"""Score audit variants with Malinois (Gosai 2024, MPRA-trained; K562, HepG2, SK-N-SH) and test
against measured allelic effects. 200 bp window centred on the variant, MPRA flanks, fwd+RC averaged.
Usage: python malinois_score.py --boda DIR --artifact DIR --genome hg38.2bit --rows CSV.GZ --audit-code DIR --out DIR"""
import argparse, importlib.util, inspect, os, sys, types
from pathlib import Path
import numpy as np, pandas as pd, torch, py2bit

ap = argparse.ArgumentParser()
for k in ['boda', 'artifact', 'genome', 'rows', 'audit_code', 'out']:
    ap.add_argument('--' + k.replace('_', '-'), dest=k, required=True)
a = ap.parse_args(); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, a.audit_code)
from audit_engine import errors  # noqa: E402
torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', '8')))

# import boda's model code without installing the package or its unrelated dependencies
lt = types.ModuleType('lightning'); pl = types.ModuleType('lightning.pytorch')
pl.LightningModule = torch.nn.Module
pl.__getattr__ = lambda name: torch.nn.Module   # Trainer, LightningDataModule, callbacks
lt.pytorch = pl
sys.modules.setdefault('lightning', lt); sys.modules.setdefault('lightning.pytorch', pl)
B = Path(a.boda) / 'boda'
for nm, sub in [('bd', ''), ('bd.common', 'common'), ('bd.model', 'model')]:
    m = types.ModuleType(nm); m.__path__ = [str(B / sub)]; sys.modules.setdefault(nm, m)
basset = importlib.import_module('bd.model.basset')
constants = importlib.import_module('bd.common.constants')
MPRA_UPSTREAM, MPRA_DOWNSTREAM = constants.MPRA_UPSTREAM, constants.MPRA_DOWNSTREAM

ck = torch.load(Path(a.artifact) / 'torch_checkpoint.pt', map_location='cpu', weights_only=False)
assert ck['model_module'] == 'BassetBranched', ck['model_module']
hp = vars(ck['model_hparams']); ok = set(inspect.signature(basset.BassetBranched.__init__).parameters)
net = basset.BassetBranched(**{k: v for k, v in hp.items() if k in ok})
net.load_state_dict(ck['model_state_dict']); net.eval()
L = hp['input_len']; lp = (L - 200) // 2; rp = (L - 200) - lp
FL, FR = MPRA_UPSTREAM[-lp:], MPRA_DOWNSTREAM[:rp]; IDX = {c: i for i, c in enumerate('ACGT')}

def onehot(seqs):
    x = np.zeros((len(seqs), 4, L), np.float32)
    for n, s in enumerate(seqs):
        for i, c in enumerate(FL + s + FR):
            j = IDX.get(c)
            if j is not None: x[n, j, i] = 1
    return torch.from_numpy(x)

def predict(seqs, bs=256):
    out = []
    with torch.inference_mode():
        for i in range(0, len(seqs), bs):
            x = onehot(seqs[i:i + bs]); out.append(((net(x) + net(x.flip(dims=[1, 2]))) / 2).numpy())
            print('scored', min(i + bs, len(seqs)), 'of', len(seqs), flush=True)
    return np.concatenate(out)

d = pd.read_csv(a.rows, low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true')].copy()
for c in d.columns:
    if c.endswith('_mean') or c.startswith('audit_') or c == 'effect':
        d[c] = pd.to_numeric(d[c], errors='coerce')
v = d.drop_duplicates('variant_id')[['variant_id', 'chrom', 'pos', 'ref', 'alt']].reset_index(drop=True)
g = py2bit.open(a.genome); refs, alts, bad = [], [], 0
for r in v.itertuples():
    p0 = int(r.pos) - 1; s0 = p0 - 150
    w = g.sequence(str(r.chrom), s0, s0 + 500).upper()
    ref, alt = str(r.ref).upper(), str(r.alt).upper(); o = p0 - s0
    if w[o:o + len(ref)] != ref: bad += 1
    wa = w[:o] + alt + w[o + len(ref):]
    refs.append(w[50:250]); c0 = 50 + (len(alt) - len(ref)) // 2; alts.append(wa[c0:c0 + 200])
g.close(); print('variants', len(v), 'ref mismatches', bad, 'input_len', L, flush=True)
skew = predict(alts) - predict(refs)
cells = ['K562', 'HepG2', 'SKNSH']
for i, c in enumerate(cells): v['malinois_' + c] = skew[:, i]
v.to_csv(OUT / 'malinois_variant_scores.csv.gz', index=False)
d = d.merge(v[['variant_id'] + ['malinois_' + c for c in cells]], on='variant_id', how='left')
d['cluster'] = np.where(d.dataset.eq('Kircher'), d.fit_group.astype(str), d.genomic_block.astype(str))

def boot_r(p, y, cl, nb=2000, seed=20261014):
    m = np.isfinite(p) & np.isfinite(y); p, y, cl = p[m], y[m], cl[m]
    u, inv = np.unique(cl, return_inverse=True); s = np.zeros((len(u), 6))
    np.add.at(s, inv, np.column_stack([np.ones_like(y), p, y, p * p, y * y, p * y]))
    W = np.random.default_rng(seed).multinomial(len(u), np.ones(len(u)) / len(u), size=nb)
    T = W @ s; n = T[:, 0]
    r = (T[:, 5] - T[:, 1] * T[:, 2] / n) / np.sqrt((T[:, 3] - T[:, 1] ** 2 / n) * (T[:, 4] - T[:, 2] ** 2 / n))
    return float(np.corrcoef(p, y)[0, 1]), *np.nanpercentile(r, [2.5, 97.5])

res = []
for ds in ['Kircher', 'Abell', 'Ulirsch']:
    t = d[d.dataset.eq(ds)].copy()
    for name in ['enformer_DNASE_mean'] + ['malinois_' + c for c in cells]:
        r, lo, hi = boot_r(t[name].to_numpy(), t.effect.to_numpy(), t.cluster.to_numpy())
        h = t.copy(); h['enformer_DNASE_mean'] = t[name]
        err = float(np.nanmean(errors(h, 'enformer', 'all', source=ds == 'Kircher')))
        res.append(dict(dataset=ds, score=name, pearson=r, ci_lo=lo, ci_hi=hi, sign_error=err, n=len(t)))
res = pd.DataFrame(res); res.to_csv(OUT / 'malinois_results.csv', index=False)
pd.set_option('display.width', 200); print(res.round(3).to_string(index=False))
