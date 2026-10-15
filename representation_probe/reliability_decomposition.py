"""Why does external performance collapse? Decompose it into measurement noise and model failure.
Each audited row carries a replicate-derived noise_variance, so per variant we can form a reliability
= 1 - noise_variance / var(effect). Stratify by reliability and report the raw correlation, the
per-stratum noise ceiling sqrt(mean reliability), and the ratio r / ceiling, which is the model's
accuracy AFTER removing what measurement noise alone can explain.

Decision rule, fixed before the output was read:
  - If r / ceiling is roughly flat across reliability strata and stays near the pooled value, the
    collapse is NOT explained by measurement noise: the model genuinely fails.
  - If r / ceiling rises steeply with reliability and approaches the Kircher value (0.581 / 0.816
    = 0.71), then measurement noise explains most of the collapse and the models are better than
    the headline numbers suggest.
  - Kircher is the positive control and must behave as the already published audit says.
Usage: python reliability_decomposition.py --rows CSV.GZ --skew CSV.GZ --out DIR"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
for k in ['rows', 'skew', 'out']: ap.add_argument('--' + k, required=True)
a = ap.parse_args(); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
d = pd.read_csv(a.rows, low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true')].copy()
for c in ['effect', 'enformer_DNASE_mean', 'borzoi_DNASE_mean', 'noise_variance']:
    if c in d: d[c] = pd.to_numeric(d[c], errors='coerce')
d = d.merge(pd.read_csv(a.skew)[['variant_id', 'malinois_K562', 'malinois_SKNSH']], on='variant_id', how='left')
d['cluster'] = np.where(d.dataset.eq('Kircher'), d.fit_group.astype(str), d.genomic_block.astype(str))
print('noise_variance present for', int(d.noise_variance.notna().sum()), 'of', len(d), 'rows', flush=True)
print(d.groupby('dataset').noise_variance.describe()[['count', 'mean', '50%', 'max']].round(3), flush=True)

def boot(p, y, cl, nv, nb=2000, seed=20261015):
    m = np.isfinite(p) & np.isfinite(y) & np.isfinite(nv)
    p, y, cl, nv = p[m], y[m], cl[m], nv[m]
    if len(p) < 50: return None
    u, inv = np.unique(cl, return_inverse=True)
    s = np.zeros((len(u), 8))
    np.add.at(s, inv, np.column_stack([np.ones_like(y), p, y, p * p, y * y, p * y, nv, np.zeros_like(y)]))
    W = np.random.default_rng(seed).multinomial(len(u), np.ones(len(u)) / len(u), size=nb)
    T = W @ s; n = T[:, 0]
    vy = T[:, 4] / n - (T[:, 2] / n) ** 2
    r = (T[:, 5] - T[:, 1] * T[:, 2] / n) / np.sqrt((T[:, 3] - T[:, 1] ** 2 / n) * (T[:, 4] - T[:, 2] ** 2 / n))
    ceil = np.sqrt(np.clip(1 - (T[:, 6] / n) / vy, 0, 1))
    ratio = r / np.where(ceil > 0, ceil, np.nan)
    r0 = float(np.corrcoef(p, y)[0, 1]); c0 = float(np.sqrt(max(0.0, 1 - np.mean(nv) / np.var(y))))
    return dict(n=len(p), blocks=len(u), pearson=r0, r_lo=float(np.nanpercentile(r, 2.5)),
                r_hi=float(np.nanpercentile(r, 97.5)), ceiling=c0,
                ratio=float(r0 / c0) if c0 > 0 else np.nan,
                ratio_lo=float(np.nanpercentile(ratio, 2.5)), ratio_hi=float(np.nanpercentile(ratio, 97.5)),
                mean_noise_var=float(np.mean(nv)), effect_sd=float(np.std(y)))

rows = []
SCORE = {'Kircher': 'malinois_K562', 'Abell': 'malinois_SKNSH', 'Ulirsch': 'malinois_K562'}
for ds in ['Kircher', 'Abell', 'Ulirsch']:
    t = d[d.dataset.eq(ds)].dropna(subset=['effect', 'noise_variance']).copy()
    if not len(t): continue
    t['rel'] = 1 - t.noise_variance / t.effect.var()
    t['q'] = pd.qcut(t.rel, 5, labels=False, duplicates='drop')
    for score in ['enformer_DNASE_mean', SCORE[ds]]:
        for label, g in [('all', t)] + [(f'reliability_q{int(q)+1}', t[t.q.eq(q)]) for q in sorted(t.q.dropna().unique())]:
            st = boot(g[score].to_numpy(), g.effect.to_numpy(), g.cluster.to_numpy(), g.noise_variance.to_numpy())
            if st: rows.append(dict(dataset=ds, score=score, stratum=label, **st))
res = pd.DataFrame(rows); res.to_csv(OUT / 'reliability_decomposition.csv', index=False)
pd.set_option('display.width', 240); print(res.round(3).to_string(index=False))
