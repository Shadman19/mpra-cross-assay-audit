"""Paired target-assay fine-tuning: the one recipe never tried. Pretrained Malinois backbone with
only the last convolutional block unfrozen, plus a small shared scalar head. Prediction is a true
difference, delta = f(ALT) - f(REF), so swapping alleles negates it and identical alleles give zero.
Outer grouped 5-fold CV over genomic blocks, inner grouped early stopping, three seeds averaged
without picking the best. Forward and reverse complement averaged.

Gates fixed before any run (external reviewer): pooled out-of-fold Pearson >= 0.15, at least 0.10
above the frozen-feature ridge, and above the strongest frozen Malinois comparator; sign error <= 0.45
and at least 5 points below that comparator; paired block-bootstrap intervals excluding zero.
Usage: python finetune.py --boda D --artifact D --seqs NPZ --rows CSV.GZ --audit-code D --out D"""
import argparse, importlib, inspect, os, sys, types, time
from pathlib import Path
import numpy as np, pandas as pd, torch, torch.nn as nn

ap = argparse.ArgumentParser()
for k in ['boda', 'artifact', 'seqs', 'rows', 'audit_code', 'out']:
    ap.add_argument('--' + k.replace('_', '-'), dest=k, required=True)
ap.add_argument('--dataset', default='Abell'); ap.add_argument('--epochs', type=int, default=12)
ap.add_argument('--seeds', type=int, default=3); ap.add_argument('--batch', type=int, default=128)
a = ap.parse_args(); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, a.audit_code)
from audit_engine import errors  # noqa: E402
lt = types.ModuleType('lightning'); pl = types.ModuleType('lightning.pytorch')
pl.LightningModule = nn.Module
def _stub(name):
    if name.startswith('__'): raise AttributeError(name)   # keep dunders real for torch/inspect
    return nn.Module
pl.__getattr__ = _stub; lt.pytorch = pl
sys.modules.setdefault('lightning', lt); sys.modules.setdefault('lightning.pytorch', pl)
BP = Path(a.boda) / 'boda'
for nm, sub in [('bd', ''), ('bd.common', 'common'), ('bd.model', 'model')]:
    m = types.ModuleType(nm); m.__path__ = [str(BP / sub)]; sys.modules.setdefault(nm, m)
basset = importlib.import_module('bd.model.basset')
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
ck = torch.load(Path(a.artifact) / 'torch_checkpoint.pt', map_location='cpu', weights_only=False)
HP = vars(ck['model_hparams']); SIG = set(inspect.signature(basset.BassetBranched.__init__).parameters)

class Delta(nn.Module):
    """Shared trunk scored on REF and ALT; the output is their difference."""
    def __init__(self):
        super().__init__()
        self.net = basset.BassetBranched(**{k: v for k, v in HP.items() if k in SIG})
        self.net.load_state_dict(ck['model_state_dict'])
        for p in self.net.parameters(): p.requires_grad = False
        for p in self.net.conv3.parameters(): p.requires_grad = True   # last conv block only
        with torch.no_grad():
            dim = self.net.encode(torch.zeros(1, 4, HP['input_len'])).shape[-1]
        self.head = nn.Sequential(nn.Linear(dim, 64), nn.ReLU(), nn.Dropout(0.2), nn.Linear(64, 1))
    def score(self, x):
        return self.head(self.net.encode(x)).squeeze(-1)
    def forward(self, r, m):
        return (self.score(m) - self.score(r) + self.score(m.flip(dims=[1, 2])) - self.score(r.flip(dims=[1, 2]))) / 2

z = np.load(a.seqs)
ref, alt, ids = z['ref'], z['alt'], z['row_id'].astype(str)
d = pd.read_csv(a.rows, low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true') & d.dataset.eq(a.dataset)].copy()
for c in ['effect', 'enformer_DNASE_mean']: d[c] = pd.to_numeric(d[c], errors='coerce')
pos = pd.Series(np.arange(len(ids)), index=ids); d['row_id'] = d.row_id.astype(str)
d = d[d.row_id.isin(pos.index)].dropna(subset=['effect']).reset_index(drop=True)
d['_i'] = pos.loc[d.row_id].to_numpy(); d['cluster'] = d.genomic_block.astype(str)
print(a.dataset, len(d), 'rows', d.cluster.nunique(), 'blocks', 'device', DEV, flush=True)

def onehot(idx):
    x = np.zeros((len(idx), 4, ref.shape[1]), np.float32)
    for n, row in enumerate(idx):
        for j, k in enumerate(row):
            if 0 <= k < 4: x[n, k, j] = 1
    return torch.from_numpy(x)

def predict(mdl, idx, bs=256):
    if len(idx) == 0: return np.zeros(0, np.float32)
    out = []
    with torch.inference_mode():
        for i in range(0, len(idx), bs):
            ii = d._i.to_numpy()[idx[i:i + bs]]
            out.append(mdl(onehot(ref[ii]).to(DEV), onehot(alt[ii]).to(DEV)).float().cpu().numpy())
    return np.concatenate(out)

def run_fold(tr, va, te, seed):
    torch.manual_seed(seed); np.random.seed(seed)
    mdl = Delta().to(DEV)
    opt = torch.optim.AdamW([p for p in mdl.parameters() if p.requires_grad], lr=1e-4, weight_decay=1e-4)
    y = torch.tensor(d.effect.to_numpy(), dtype=torch.float32)
    best, best_state, bad = -9e9, None, 0
    for ep in range(a.epochs):
        mdl.train(); order = np.random.permutation(tr)
        for i in range(0, len(order), a.batch):
            b = order[i:i + a.batch]; ii = d._i.to_numpy()[b]
            r = onehot(ref[ii]).to(DEV); m = onehot(alt[ii]).to(DEV)
            loss = nn.functional.mse_loss(mdl(r, m), y[b].to(DEV))
            opt.zero_grad(); loss.backward(); opt.step()
        mdl.eval(); pv = predict(mdl, va)
        r_va = float(np.corrcoef(pv, d.effect.to_numpy()[va])[0, 1]) if np.std(pv) > 0 else -9e9
        print(f'  seed {seed} ep {ep} val r {r_va:.4f}', flush=True)
        if r_va > best: best, best_state, bad = r_va, {k: v.detach().clone() for k, v in mdl.state_dict().items()}, 0
        else:
            bad += 1
            if bad >= 3: break
    if best_state: mdl.load_state_dict(best_state)
    mdl.eval(); return predict(mdl, te)

_m = Delta().to(DEV).eval()
with torch.inference_mode():
    q = onehot(ref[:4]).to(DEV); same = _m(q, q).abs().max().item()
    swap = (_m(onehot(ref[:4]).to(DEV), onehot(alt[:4]).to(DEV)) + _m(onehot(alt[:4]).to(DEV), onehot(ref[:4]).to(DEV))).abs().max().item()
print('sanity identical-allele delta', same, 'swap antisymmetry', swap, flush=True)
assert same < 1e-4 and swap < 1e-4, 'siamese difference is not antisymmetric'
del _m

ug = np.unique(d.cluster); rng = np.random.default_rng(20261016); rng.shuffle(ug)
fold = pd.Series({g: i % 5 for i, g in enumerate(ug)}).loc[d.cluster].to_numpy()
oof = np.zeros((a.seeds, len(d)))
for k in range(5):
    te = np.flatnonzero(fold == k); rest = np.flatnonzero(fold != k)
    if len(te) == 0 or len(rest) < 10: continue
    inner = np.unique(d.cluster.to_numpy()[rest]); rng.shuffle(inner)
    vg = set(inner[:max(1, len(inner) // 5)])
    va = rest[np.isin(d.cluster.to_numpy()[rest], list(vg))]
    tr = rest[~np.isin(d.cluster.to_numpy()[rest], list(vg))]
    for s in range(a.seeds):
        t0 = time.monotonic(); oof[s, te] = run_fold(tr, va, te, 20261016 + s)
        print(f'fold {k} seed {s} done in {time.monotonic()-t0:.0f}s', flush=True)
p = oof.mean(0); d['finetuned'] = p
d.to_csv(OUT / f'finetune_{a.dataset}_predictions.csv.gz', index=False)

def boot(pred, nb=4000, seed=20261016):
    y = d.effect.to_numpy(); cl = d.cluster.to_numpy()
    m = np.isfinite(pred) & np.isfinite(y); pred, y, cl = pred[m], y[m], cl[m]
    u, inv = np.unique(cl, return_inverse=True); s = np.zeros((len(u), 6))
    np.add.at(s, inv, np.column_stack([np.ones_like(y), pred, y, pred * pred, y * y, pred * y]))
    W = np.random.default_rng(seed).multinomial(len(u), np.ones(len(u)) / len(u), size=nb)
    T = W @ s; n = T[:, 0]
    return (T[:, 5] - T[:, 1] * T[:, 2] / n) / np.sqrt((T[:, 3] - T[:, 1] ** 2 / n) * (T[:, 4] - T[:, 2] ** 2 / n))

rows = []
for name in ['finetuned', 'enformer_DNASE_mean']:
    r = boot(d[name].to_numpy()); r0 = float(np.corrcoef(d[name], d.effect)[0, 1])
    e = float(np.nanmean(errors(d.assign(enformer_DNASE_mean=d[name]), 'enformer', 'all', source=False)))
    rows.append(dict(score=name, pearson=r0, ci_lo=float(np.nanpercentile(r, 2.5)),
                     ci_hi=float(np.nanpercentile(r, 97.5)), sign_error=e, n=len(d), blocks=d.cluster.nunique()))
dif = boot(d.finetuned.to_numpy()) - boot(d.enformer_DNASE_mean.to_numpy())
rows.append(dict(score='finetuned_minus_enformer', pearson=rows[0]['pearson'] - rows[1]['pearson'],
                 ci_lo=float(np.nanpercentile(dif, 2.5)), ci_hi=float(np.nanpercentile(dif, 97.5))))
res = pd.DataFrame(rows); res.to_csv(OUT / f'finetune_{a.dataset}_results.csv', index=False)
pd.set_option('display.width', 200); print(res.round(3).to_string(index=False))
print('GATE pearson>=0.15:', rows[0]['pearson'] >= 0.15, '| >=0.10 over ridge 0.041:', rows[0]['pearson'] - 0.041 >= 0.10,
      '| sign error<=0.45:', rows[0]['sign_error'] <= 0.45, flush=True)
