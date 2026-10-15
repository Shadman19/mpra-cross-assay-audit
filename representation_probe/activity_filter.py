"""Does the model work where the element is active? Malinois gives an MPRA-calibrated activity
estimate for the REFERENCE sequence, which is a prospective filter: it needs no measured outcome.
Stratify the external studies by that activity and test whether accuracy rises.
Usage: python activity_filter.py --boda D --artifact D --genome 2bit --rows CSV.GZ --audit-code D --out D"""
import argparse, importlib, inspect, os, sys, types
from pathlib import Path
import numpy as np, pandas as pd, torch, py2bit

ap = argparse.ArgumentParser()
for k in ['boda', 'artifact', 'genome', 'rows', 'audit_code', 'out']:
    ap.add_argument('--' + k.replace('_', '-'), dest=k, required=True)
a = ap.parse_args(); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, a.audit_code)
from audit_engine import errors  # noqa: E402
torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', '8')))
lt = types.ModuleType('lightning'); pl = types.ModuleType('lightning.pytorch')
pl.LightningModule = torch.nn.Module; pl.__getattr__ = lambda n: torch.nn.Module; lt.pytorch = pl
sys.modules.setdefault('lightning', lt); sys.modules.setdefault('lightning.pytorch', pl)
B = Path(a.boda) / 'boda'
for nm, sub in [('bd', ''), ('bd.common', 'common'), ('bd.model', 'model')]:
    m = types.ModuleType(nm); m.__path__ = [str(B / sub)]; sys.modules.setdefault(nm, m)
basset = importlib.import_module('bd.model.basset'); constants = importlib.import_module('bd.common.constants')

ck = torch.load(Path(a.artifact) / 'torch_checkpoint.pt', map_location='cpu', weights_only=False)
hp = vars(ck['model_hparams']); ok = set(inspect.signature(basset.BassetBranched.__init__).parameters)
net = basset.BassetBranched(**{k: v for k, v in hp.items() if k in ok})
net.load_state_dict(ck['model_state_dict']); net.eval()
L = hp['input_len']; lp = (L - 200) // 2; rp = (L - 200) - lp
FL, FR = constants.MPRA_UPSTREAM[-lp:], constants.MPRA_DOWNSTREAM[:rp]
IDX = {c: i for i, c in enumerate('ACGT')}
def predict(seqs, bs=256):
    out = []
    with torch.inference_mode():
        for i in range(0, len(seqs), bs):
            ch = seqs[i:i + bs]; x = np.zeros((len(ch), 4, L), np.float32)
            for n, s in enumerate(ch):
                for j, c in enumerate(FL + s + FR):
                    k = IDX.get(c)
                    if k is not None: x[n, k, j] = 1
            x = torch.from_numpy(x); out.append(((net(x) + net(x.flip(dims=[1, 2]))) / 2).numpy())
            print('ref', min(i + bs, len(seqs)), 'of', len(seqs), flush=True)
    return np.concatenate(out)

d = pd.read_csv(a.rows, low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true')].copy()
for c in ['effect', 'enformer_DNASE_mean', 'borzoi_DNASE_mean']:
    if c in d: d[c] = pd.to_numeric(d[c], errors='coerce')
v = d.drop_duplicates('variant_id')[['variant_id', 'chrom', 'pos', 'ref', 'alt']].reset_index(drop=True)
g = py2bit.open(a.genome)
refs = [g.sequence(str(r.chrom), int(r.pos) - 1 - 100, int(r.pos) - 1 + 100).upper() for r in v.itertuples()]
g.close()
act = predict(refs); cells = ['K562', 'HepG2', 'SKNSH']
for i, c in enumerate(cells): v['act_' + c] = act[:, i]
v['act_max'] = act.max(1)
v.to_csv(OUT / 'reference_activity.csv.gz', index=False)
sk = pd.read_csv(Path(a.out).parent / 'malinois' / 'malinois_variant_scores.csv.gz')
d = d.merge(v[['variant_id'] + ['act_' + c for c in cells] + ['act_max']], on='variant_id', how='left')
d = d.merge(sk[['variant_id', 'malinois_K562', 'malinois_HepG2', 'malinois_SKNSH']], on='variant_id', how='left')
d['cluster'] = np.where(d.dataset.eq('Kircher'), d.fit_group.astype(str), d.genomic_block.astype(str))

def boot_r(p, y, cl, nb=2000, seed=20261015):
    m = np.isfinite(p) & np.isfinite(y)
    p, y, cl = p[m], y[m], cl[m]
    if len(p) < 30: return np.nan, np.nan, np.nan, len(p)
    u, inv = np.unique(cl, return_inverse=True); s = np.zeros((len(u), 6))
    np.add.at(s, inv, np.column_stack([np.ones_like(y), p, y, p * p, y * y, p * y]))
    W = np.random.default_rng(seed).multinomial(len(u), np.ones(len(u)) / len(u), size=nb)
    T = W @ s; n = T[:, 0]
    r = (T[:, 5] - T[:, 1] * T[:, 2] / n) / np.sqrt((T[:, 3] - T[:, 1] ** 2 / n) * (T[:, 4] - T[:, 2] ** 2 / n))
    return float(np.corrcoef(p, y)[0, 1]), *np.nanpercentile(r, [2.5, 97.5]), len(p)

rows = []
MATCH = {'Abell': 'act_max', 'Ulirsch': 'act_K562', 'Kircher': 'act_max'}
SCORE = {'Abell': 'malinois_SKNSH', 'Ulirsch': 'malinois_K562', 'Kircher': 'malinois_K562'}
for ds in ['Kircher', 'Abell', 'Ulirsch']:
    t = d[d.dataset.eq(ds)].copy(); actcol = MATCH[ds]
    t['q'] = pd.qcut(t[actcol], 5, labels=False, duplicates='drop')
    for score in ['enformer_DNASE_mean', SCORE[ds]]:
        for q in sorted(t.q.dropna().unique()):
            s = t[t.q.eq(q)]
            r, lo, hi, n = boot_r(s[score].to_numpy(), s.effect.to_numpy(), s.cluster.to_numpy())
            h = s.copy(); h['enformer_DNASE_mean'] = s[score]
            err = float(np.nanmean(errors(h, 'enformer', 'all', source=ds == 'Kircher'))) if len(s) else np.nan
            rows.append(dict(dataset=ds, score=score, activity_quintile=int(q) + 1, pearson=r,
                             ci_lo=lo, ci_hi=hi, sign_error=err, n=n))
    for frac in [0.5, 0.25, 0.1]:
        thr = t[actcol].quantile(1 - frac); s = t[t[actcol] >= thr]
        for score in ['enformer_DNASE_mean', SCORE[ds]]:
            r, lo, hi, n = boot_r(s[score].to_numpy(), s.effect.to_numpy(), s.cluster.to_numpy())
            h = s.copy(); h['enformer_DNASE_mean'] = s[score]
            err = float(np.nanmean(errors(h, 'enformer', 'all', source=ds == 'Kircher')))
            rows.append(dict(dataset=ds, score=score, activity_quintile=f'top{int(frac * 100)}pct',
                             pearson=r, ci_lo=lo, ci_hi=hi, sign_error=err, n=n))
res = pd.DataFrame(rows); res.to_csv(OUT / 'activity_filter_results.csv', index=False)
pd.set_option('display.width', 220); print(res.round(3).to_string(index=False))
