"""Can these features predict MPRA allelic effects when trained on the target assay itself?
(1) within-study grouped 5-fold CV; (2) Kircher->target with an extended alpha grid; (3) few-label
learning curve on the target, with and without the Kircher-trained prediction as an extra feature.
Usage: python solve.py --features DIR --rows CSV.GZ --out DIR"""
import argparse, glob
from pathlib import Path
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
for k in ['features', 'rows', 'out']: ap.add_argument('--' + k, required=True)
a = ap.parse_args(); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261014; ALPHAS = 10.0 ** np.arange(-1, 8.01, 0.5); NB = 1000

parts = {'B': [], 'C': [], 'D': []}; ids = []
for f in sorted(glob.glob(str(Path(a.features) / 'shard*' / '*.npz'))):
    if f.endswith('.tmp.npz'): continue
    z = np.load(f); ids.append(z['row_id'].astype(str))
    for k in parts: parts[k].append(z[k].astype(np.float32))
ids = np.concatenate(ids); F = {k: np.concatenate(v) for k, v in parts.items()}
d = pd.read_csv(a.rows, low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true')].copy(); d['row_id'] = d.row_id.astype(str)
d['effect'] = pd.to_numeric(d.effect, errors='coerce'); d['enformer_DNASE_mean'] = pd.to_numeric(d.enformer_DNASE_mean, errors='coerce')
pos = pd.Series(np.arange(len(ids)), index=ids); assert d.row_id.isin(pos.index).all(), 'missing features'
d['_i'] = pos.loc[d.row_id].to_numpy()
d['cluster'] = np.where(d.dataset.eq('Kircher'), d.fit_group.astype(str), d.genomic_block.astype(str))
print({k: int((d.dataset == k).sum()) for k in ['Kircher', 'Abell', 'Ulirsch']}, flush=True)

def folds(groups, k, seed):
    ug = np.unique(groups); np.random.default_rng(seed).shuffle(ug)
    m = {g: i % k for i, g in enumerate(ug)}; return np.array([m[g] for g in groups])

def path_preds(Xtr, ytr, Xte):
    """Ridge predictions on Xte for every alpha; standardise on train; primal or dual form."""
    mu, sd = Xtr.mean(0), Xtr.std(0); sd[sd == 0] = 1
    Z = (Xtr - mu) / sd; Zt = (Xte - mu) / sd; y0 = ytr.mean(); yc = ytr - y0
    if Z.shape[0] >= Z.shape[1]:
        w, V = np.linalg.eigh(Z.T @ Z); q = V.T @ (Z.T @ yc); P = Zt @ V
        return np.stack([P @ (q / (w + al)) for al in ALPHAS], 1) + y0
    w, U = np.linalg.eigh(Z @ Z.T); q = U.T @ yc; P = (Zt @ Z.T) @ U
    return np.stack([P @ (q / (w + al)) for al in ALPHAS], 1) + y0

def fit_predict(Xtr, ytr, gtr, Xte, seed=SEED):
    k = min(5, len(np.unique(gtr))); fo = folds(gtr, k, seed); sse = np.zeros(len(ALPHAS))
    for j in range(k):
        tr, te = fo != j, fo == j
        sse += ((path_preds(Xtr[tr], ytr[tr], Xtr[te]) - ytr[te, None]) ** 2).sum(0)
    j = int(np.argmin(sse)); return path_preds(Xtr, ytr, Xte)[:, j], float(ALPHAS[j]), k

def boot_r(p, y, cl, seed=SEED):
    ok = np.isfinite(p) & np.isfinite(y); p, y, cl = p[ok], y[ok], cl[ok]
    u, inv = np.unique(cl, return_inverse=True); s = np.zeros((len(u), 6))
    np.add.at(s, inv, np.column_stack([np.ones_like(y), p, y, p * p, y * y, p * y]))
    W = np.random.default_rng(seed).multinomial(len(u), np.ones(len(u)) / len(u), size=NB); T = W @ s; n = T[:, 0]
    r = (T[:, 5] - T[:, 1] * T[:, 2] / n) / np.sqrt((T[:, 3] - T[:, 1] ** 2 / n) * (T[:, 4] - T[:, 2] ** 2 / n))
    return float(np.corrcoef(p, y)[0, 1]), *np.nanpercentile(r, [2.5, 97.5])

def X_of(g, fs): return np.hstack([F[f][g._i.to_numpy()] for f in fs])
res = []
def add(**kw): res.append(kw); print(kw, flush=True)

# (2) Kircher -> target, extended alpha grid (also gives the Kircher prior used in (3))
kir = d[d.dataset.eq('Kircher')]; prior = {}
for f in ['B', 'C', 'D']:
    Xk = X_of(kir, [f])
    for t in ['Abell', 'Ulirsch']:
        te = d[d.dataset.eq(t)]
        p, al, _ = fit_predict(Xk, kir.effect.to_numpy(), kir.cluster.to_numpy(), X_of(te, [f]))
        if f == 'B': prior[t] = p
        r, lo, hi = boot_r(p, te.effect.to_numpy(), te.cluster.to_numpy())
        add(exp='kircher_to_target', target=t, features=f, budget=None, pearson=r, ci_lo=lo, ci_hi=hi, alpha=al)
    del Xk

# (1) within-study grouped 5-fold CV
for t in ['Abell', 'Ulirsch']:
    g = d[d.dataset.eq(t)].reset_index(drop=True); y = g.effect.to_numpy(); cl = g.cluster.to_numpy()
    r, lo, hi = boot_r(g.enformer_DNASE_mean.to_numpy(), y, cl)
    add(exp='within_study', target=t, features='A_raw', budget=None, pearson=r, ci_lo=lo, ci_hi=hi, alpha=None)
    fo = folds(cl, 5, SEED)
    for fs in [['B'], ['C'], ['D'], ['B', 'C', 'D']]:
        X = X_of(g, fs); oof = np.full(len(g), np.nan); als = []
        for j in range(5):
            tr, te = fo != j, fo == j
            oof[te], al, _ = fit_predict(X[tr], y[tr], cl[tr], X[te]); als.append(al)
        r, lo, hi = boot_r(oof, y, cl)
        add(exp='within_study', target=t, features='+'.join(fs), budget=None, pearson=r, ci_lo=lo, ci_hi=hi, alpha=str(als))
        del X

# (3) few-label learning curve, feature set B, fixed 50% test blocks
for t in ['Abell', 'Ulirsch']:
    g = d[d.dataset.eq(t)].reset_index(drop=True); y = g.effect.to_numpy(); cl = g.cluster.to_numpy()
    XB = X_of(g, ['B']); pk = prior[t]; ub = np.unique(cl); rng = np.random.default_rng(SEED); rng.shuffle(ub)
    test_b = set(ub[:len(ub) // 2]); pool = ub[len(ub) // 2:]; te = np.isin(cl, list(test_b))
    add(exp='curve', target=t, features='kircher_prior_only', budget=0, pearson=float(np.corrcoef(pk[te], y[te])[0, 1]), ci_lo=None, ci_hi=None, alpha=None)
    for frac in [0.1, 0.25, 0.5, 1.0]:
        for variant in ['target_only', 'target_plus_prior']:
            rs = []
            for rep in range(5 if frac < 1 else 1):
                r2 = np.random.default_rng(SEED + rep); pick = r2.choice(pool, max(2, int(round(frac * len(pool)))), replace=False)
                tr = np.isin(cl, pick)
                X = XB if variant == 'target_only' else np.hstack([XB, pk[:, None]])
                p, al, k = fit_predict(X[tr], y[tr], cl[tr], X[te], SEED + rep)
                rs.append(np.corrcoef(p, y[te])[0, 1])
            add(exp='curve', target=t, features=variant, budget=frac, pearson=float(np.mean(rs)), ci_lo=float(np.min(rs)), ci_hi=float(np.max(rs)), alpha=f'n_blocks={len(pick)};inner_k={k}')

out = pd.DataFrame(res); out.to_csv(OUT / 'solve_results.csv', index=False)
pd.set_option('display.width', 220); print(out.round(3).to_string(index=False))
