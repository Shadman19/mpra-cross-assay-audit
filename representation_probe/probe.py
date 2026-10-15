"""CPU analysis for the preregistered representation probe (see PREREGISTRATION.md).
Usage: python probe.py --features <dir with shardNN/*.npz> --rows <annotated_audit_variants.csv.gz>
       --audit-code <subgroup_audit/code> --out <dir>
"""
import argparse, glob, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ap = argparse.ArgumentParser()
ap.add_argument('--features', required=True); ap.add_argument('--rows', required=True)
ap.add_argument('--audit-code', required=True); ap.add_argument('--out', required=True)
a = ap.parse_args()
sys.path.insert(0, a.audit_code)
from audit_engine import errors  # noqa: E402

ALPHAS = 10.0 ** np.arange(-1, 5.5, 0.5)   # 13 values, 1e-1 .. 1e5
B_DRAWS, FOLDS, SEED = 2000, 5, 20271014
OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)

# ---- load features --------------------------------------------------------------------------
parts = {'B': [], 'C': [], 'D': []}; ids = []
files = sorted(glob.glob(str(Path(a.features) / 'shard*' / '*.npz')))
files = [f for f in files if not f.endswith('.tmp.npz')]
for f in files:
    z = np.load(f); ids.append(z['row_id'])
    for k in parts: parts[k].append(z[k].astype(np.float32))
ids = np.concatenate(ids); feats = {k: np.concatenate(v) for k, v in parts.items()}
d = pd.read_csv(a.rows, low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true')].copy()
for c in d.columns:
    if c.endswith('_mean') or c.startswith('audit_') or c == 'effect': d[c] = pd.to_numeric(d[c], errors='coerce')
pos = pd.Series(np.arange(len(ids)), index=ids.astype(str)); d['row_id'] = d.row_id.astype(str)
missing = (~d.row_id.isin(pos.index)).sum()
assert missing == 0, f'{missing} primary rows lack features; scoring is incomplete'
d = d.assign(_i=pos.loc[d.row_id].to_numpy()).reset_index(drop=True)
d['cluster'] = np.where(d.dataset.eq('Kircher'), d.fit_group.astype(str), d.genomic_block.astype(str))

# ---- ridge with grouped CV over alpha (eigendecomposition per fold) ----------------------
def fit_ridge(X, y, groups):
    mu, sd = X.mean(0), X.std(0); sd[sd == 0] = 1
    Z = (X - mu) / sd; yc = y - y.mean()
    ug = np.unique(groups); rng = np.random.default_rng(SEED); rng.shuffle(ug)
    fold_of = {g: i % FOLDS for i, g in enumerate(ug)}; fold = np.array([fold_of[g] for g in groups])
    sse = np.zeros(len(ALPHAS))
    for k in range(FOLDS):
        tr, te = fold != k, fold == k
        Zt, yt = Z[tr], yc[tr] - yc[tr].mean()
        w_, V = np.linalg.eigh(Zt.T @ Zt); proj = V.T @ (Zt.T @ yt)
        for j, al in enumerate(ALPHAS):
            beta = V @ (proj / (w_ + al))
            sse[j] += ((Z[te] @ beta + yc[tr].mean() - yc[te]) ** 2).sum()
    al = ALPHAS[int(np.argmin(sse))]
    w_, V = np.linalg.eigh(Z.T @ Z); beta = V @ ((V.T @ (Z.T @ yc)) / (w_ + al))
    return dict(mu=mu, sd=sd, beta=beta, y0=y.mean(), alpha=float(al))

def predict(m, X): return ((X - m['mu']) / m['sd']) @ m['beta'] + m['y0']

# ---- evaluation ----------------------------------------------------------------------------
def boot_pearson(g, preds):
    cl = g.cluster.to_numpy(); u, inv = np.unique(cl, return_inverse=True)
    y = g.effect.to_numpy(); W = np.random.default_rng(SEED).multinomial(len(u), np.ones(len(u)) / len(u), size=B_DRAWS)
    out = {}
    for name, p in preds.items():
        s = np.zeros((len(u), 6))
        np.add.at(s, inv, np.column_stack([np.ones_like(y), p, y, p * p, y * y, p * y]))
        T = W @ s; n = T[:, 0]
        num = T[:, 5] - T[:, 1] * T[:, 2] / n
        den = np.sqrt((T[:, 3] - T[:, 1] ** 2 / n) * (T[:, 4] - T[:, 2] ** 2 / n))
        out[name] = (pearsonr(p, y)[0], num / den)
    return out

def sign_error(g, p, source):
    h = g.copy(); h['enformer_DNASE_mean'] = p
    e = errors(h, 'enformer', 'all', source=source); return float(np.nanmean(e))

results, verdict = [], {}
for train, test in [('Kircher', 'Abell'), ('Kircher', 'Ulirsch'), ('Abell', 'Ulirsch')]:
    tr = d[d.dataset.eq(train)]; te = d[d.dataset.eq(test)].reset_index(drop=True)
    preds = {'A_DNASE_raw': te.enformer_DNASE_mean.to_numpy()}
    alphas = {}
    for F in ['B', 'C', 'D']:
        m = fit_ridge(feats[F][tr._i.to_numpy()], tr.effect.to_numpy(), tr.cluster.to_numpy())
        preds[F] = predict(m, feats[F][te._i.to_numpy()]); alphas[F] = m['alpha']
    bp = boot_pearson(te, preds)
    for name, (r, bs) in bp.items():
        lo, hi = np.nanpercentile(bs, [2.5, 97.5])
        results.append(dict(train=train, test=test, features=name, pearson=r, ci_lo=lo, ci_hi=hi,
                            sign_error=sign_error(te, preds[name], False), alpha=alphas.get(name)))
    dif = bp['C'][1] - bp['B'][1]; lo, hi = np.nanpercentile(dif, [2.5, 97.5])
    results.append(dict(train=train, test=test, features='C minus B', pearson=bp['C'][0] - bp['B'][0], ci_lo=lo, ci_hi=hi))
    if train == 'Kircher':
        verdict[test] = bool(lo > 0 and bp['C'][0] >= 0.15)

res = pd.DataFrame(results); res.to_csv(OUT / 'probe_results.csv', index=False)
success = all(verdict.values()) and len(verdict) == 2
summary = {'success': success, 'per_target': verdict, 'n_rows': int(len(d)), 'feature_files': len(files)}
(OUT / 'probe_summary.json').write_text(json.dumps(summary, indent=2))
pd.set_option('display.width', 200)
print(res.round(3).to_string(index=False))
print('PREREGISTERED SUCCESS:', success, verdict)
