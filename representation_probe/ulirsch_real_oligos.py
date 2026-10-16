"""Score Ulirsch with the ACTUAL assayed oligo sequences from the study design, not a reconstructed
genomic window. Each variant was assayed in several tile positions; predictions are averaged over the
tiles that the study actually used, forward and reverse complement averaged, skew = Mut minus Ref.
Usage: python ulirsch_real_oligos.py --boda D --artifact D --design TSV.GZ --rows CSV.GZ --audit-code D --out D"""
import argparse, importlib, inspect, os, sys, types
from pathlib import Path
import numpy as np, pandas as pd, torch

ap = argparse.ArgumentParser()
for k in ['boda', 'artifact', 'design', 'rows', 'audit_code', 'out']:
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
L = hp['input_len']; IDX = {c: i for i, c in enumerate('ACGT')}

def predict(seqs, bs=256):
    out = []
    with torch.inference_mode():
        for i in range(0, len(seqs), bs):
            ch = seqs[i:i + bs]; x = np.zeros((len(ch), 4, L), np.float32)
            for n, s in enumerate(ch):
                ins = s[:L]; lp = (L - len(ins)) // 2; rp = L - len(ins) - lp
                full = constants.MPRA_UPSTREAM[-lp:] + ins + constants.MPRA_DOWNSTREAM[:rp]
                for j, c in enumerate(full[:L]):
                    k = IDX.get(c)
                    if k is not None: x[n, k, j] = 1
            x = torch.from_numpy(x); out.append(((net(x) + net(x.flip(dims=[1, 2]))) / 2).numpy())
            print('oligo', min(i + bs, len(seqs)), 'of', len(seqs), flush=True)
    return np.concatenate(out)

des = pd.read_csv(a.design, sep='\t', header=None, names=['id', 'seq'])
des['base'] = des.id.str.split('_').str[0]; des['tag'] = des.id.str.split('_').str[1]
des = des[des.tag.isin(['Ref', 'Mut'])].copy()
d = pd.read_csv(a.rows, low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true') & d.dataset.eq('Ulirsch')].copy()
for c in ['effect', 'enformer_DNASE_mean']: d[c] = pd.to_numeric(d[c], errors='coerce')
# legacy_row_id "chr1:3684954:C>A"  ->  design base "1:3684954:C:A"
key = d.legacy_row_id.astype(str).str.replace('^chr', '', regex=True).str.replace('>', ':', regex=False)
d['dkey'] = key
hit = d.dkey.isin(set(des.base))
print('matched', int(hit.sum()), 'of', len(d), 'Ulirsch primary rows', flush=True)
des = des[des.base.isin(set(d.dkey))]
pred = predict(des.seq.str.upper().tolist())
cells = ['K562', 'HepG2', 'SKNSH']
for i, c in enumerate(cells): des['p_' + c] = pred[:, i]
agg = des.groupby(['base', 'tag'])[['p_' + c for c in cells]].mean().reset_index()
w = agg.pivot(index='base', columns='tag')
sk = pd.DataFrame({'dkey': w.index})
for c in cells: sk['oligo_' + c] = (w[('p_' + c, 'Mut')] - w[('p_' + c, 'Ref')]).to_numpy()
sk.to_csv(OUT / 'ulirsch_oligo_scores.csv.gz', index=False)
d = d.merge(sk, on='dkey', how='left'); d['cluster'] = d.genomic_block.astype(str)

def boot(p, y, cl, nb=2000, seed=20261016):
    m = np.isfinite(p) & np.isfinite(y); p, y, cl = p[m], y[m], cl[m]
    u, inv = np.unique(cl, return_inverse=True); s = np.zeros((len(u), 6))
    np.add.at(s, inv, np.column_stack([np.ones_like(y), p, y, p * p, y * y, p * y]))
    W = np.random.default_rng(seed).multinomial(len(u), np.ones(len(u)) / len(u), size=nb)
    T = W @ s; n = T[:, 0]
    r = (T[:, 5] - T[:, 1] * T[:, 2] / n) / np.sqrt((T[:, 3] - T[:, 1] ** 2 / n) * (T[:, 4] - T[:, 2] ** 2 / n))
    return float(np.corrcoef(p, y)[0, 1]), float(np.nanpercentile(r, 2.5)), float(np.nanpercentile(r, 97.5)), len(p), len(u)

rows = []
for name in ['enformer_DNASE_mean'] + ['oligo_' + c for c in cells]:
    g = d.dropna(subset=[name, 'effect'])
    r, lo, hi, n, b = boot(g[name].to_numpy(), g.effect.to_numpy(), g.cluster.to_numpy())
    e = float(np.nanmean(errors(g.assign(enformer_DNASE_mean=g[name]), 'enformer', 'all', source=False)))
    rows.append(dict(score=name, pearson=r, ci_lo=lo, ci_hi=hi, sign_error=e, n=n, blocks=b))
res = pd.DataFrame(rows); res.to_csv(OUT / 'ulirsch_oligo_results.csv', index=False)
pd.set_option('display.width', 200); print(res.round(3).to_string(index=False))
