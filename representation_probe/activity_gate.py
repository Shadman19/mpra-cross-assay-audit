"""Score the activity filter against the success gate fixed before the results were read.
Primary forecast is Enformer DNASE; Malinois heads are secondary. Primary retention is the top 25
percent by predicted reference activity (max over heads for Abell, K562 for Ulirsch). The comparator
is a magnitude rule retaining the SAME number of variants in each genomic block, ranked by predicted
effect size. Also reports per-stratum replicate noise ceilings and noise-adjusted correlation.
Usage: python activity_gate.py --rows CSV.GZ --act CSV.GZ --skew CSV.GZ --audit-code DIR --out DIR"""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
for k in ['rows', 'act', 'skew', 'audit_code', 'out']: ap.add_argument('--' + k.replace('_', '-'), dest=k, required=True)
a = ap.parse_args(); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, a.audit_code)
from audit_engine import errors  # noqa: E402

d = pd.read_csv(a.rows, low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true')].copy()
for c in ['effect', 'enformer_DNASE_mean', 'noise_variance']: d[c] = pd.to_numeric(d[c], errors='coerce')
d = d.merge(pd.read_csv(a.act)[['variant_id', 'act_K562', 'act_max']], on='variant_id', how='left')
d = d.merge(pd.read_csv(a.skew)[['variant_id', 'malinois_K562', 'malinois_SKNSH']], on='variant_id', how='left')
d['cluster'] = d.genomic_block.astype(str)
B, SEED = 4000, 20261015

def stats(g, score):
    p, y, cl = g[score].to_numpy(), g.effect.to_numpy(), g.cluster.to_numpy()
    m = np.isfinite(p) & np.isfinite(y); p, y, cl = p[m], y[m], cl[m]
    u, inv = np.unique(cl, return_inverse=True)
    e = errors(g.assign(enformer_DNASE_mean=g[score]), 'enformer', 'all', source=False)
    s = np.zeros((len(u), 6)); np.add.at(s, inv, np.column_stack([np.ones_like(y), p, y, p * p, y * y, p * y]))
    W = np.random.default_rng(SEED).multinomial(len(u), np.ones(len(u)) / len(u), size=B)
    T = W @ s; n = T[:, 0]
    r = (T[:, 5] - T[:, 1] * T[:, 2] / n) / np.sqrt((T[:, 3] - T[:, 1] ** 2 / n) * (T[:, 4] - T[:, 2] ** 2 / n))
    nv = g.noise_variance.mean(); ceil = float(np.sqrt(max(0.0, 1 - nv / np.var(y)))) if np.isfinite(nv) else np.nan
    r0 = float(np.corrcoef(p, y)[0, 1])
    return dict(n=len(p), blocks=len(u), pearson=r0, ci_lo=float(np.nanpercentile(r, 2.5)),
                ci_hi=float(np.nanpercentile(r, 97.5)), sign_error=float(np.nanmean(e)),
                effect_sd=float(np.std(y)), noise_ceiling=ceil,
                r_over_ceiling=float(r0 / ceil) if ceil and ceil > 0 else np.nan), r, (s, u, inv, W)

rows = []
CFG = {'Abell': ('act_max', 'malinois_SKNSH'), 'Ulirsch': ('act_K562', 'malinois_K562')}
for ds, (actcol, sk) in CFG.items():
    t = d[d.dataset.eq(ds)].dropna(subset=[actcol, 'effect']).copy()
    for frac in [0.25, 0.10]:
        sel = t[t[actcol] >= t[actcol].quantile(1 - frac)]
        keep = sel.groupby('cluster').size()
        mag = pd.concat([g.reindex(g.enformer_DNASE_mean.abs().sort_values(ascending=False).index)
                         .head(int(keep.get(b, 0))) for b, g in t.groupby('cluster')])
        for score in ['enformer_DNASE_mean', sk]:
            for label, g in [('all', t), (f'activity_top{int(frac*100)}', sel),
                             (f'magnitude_blockmatched_top{int(frac*100)}', mag)]:
                st, r, _ = stats(g, score)
                rows.append(dict(dataset=ds, score=score, selection=label, retention=frac, **st))
            # paired difference, activity vs block-matched magnitude, same block weights
            _, ra, _ = stats(sel, score); _, rm, _ = stats(mag, score)
            dif = ra - rm; 
            rows.append(dict(dataset=ds, score=score, selection=f'activity_minus_magnitude_top{int(frac*100)}',
                             retention=frac, n=len(sel), blocks=sel.cluster.nunique(),
                             pearson=float(np.nanmean(ra) - np.nanmean(rm)),
                             ci_lo=float(np.nanpercentile(dif, 2.5)), ci_hi=float(np.nanpercentile(dif, 97.5))))
res = pd.DataFrame(rows); res.to_csv(OUT / 'activity_gate_results.csv', index=False)
pd.set_option('display.width', 240); print(res.round(3).to_string(index=False))
