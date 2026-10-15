"""Control for the activity filter: is the gain real, or just range restriction?
For each prospective retention level, compare the activity-selected set against (a) random subsets of
the same size and (b) random subsets matched on the spread of the measured effect, and report the
measured-effect SD of each set. Usage: python activity_control.py --rows CSV.GZ --act CSV.GZ --skew CSV.GZ --out DIR"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
for k in ['rows', 'act', 'skew', 'out']: ap.add_argument('--' + k, required=True)
a = ap.parse_args(); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
d = pd.read_csv(a.rows, low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true')].copy()
for c in ['effect', 'enformer_DNASE_mean']: d[c] = pd.to_numeric(d[c], errors='coerce')
d = d.merge(pd.read_csv(a.act)[['variant_id', 'act_K562', 'act_max']], on='variant_id', how='left')
d = d.merge(pd.read_csv(a.skew)[['variant_id', 'malinois_K562', 'malinois_SKNSH']], on='variant_id', how='left')
d['cluster'] = np.where(d.dataset.eq('Kircher'), d.fit_group.astype(str), d.genomic_block.astype(str))
rng = np.random.default_rng(20261015)

def blockboot(p, y, cl, nb=2000, seed=1):
    m = np.isfinite(p) & np.isfinite(y); p, y, cl = p[m], y[m], cl[m]
    u, inv = np.unique(cl, return_inverse=True); s = np.zeros((len(u), 6))
    np.add.at(s, inv, np.column_stack([np.ones_like(y), p, y, p * p, y * y, p * y]))
    W = np.random.default_rng(seed).multinomial(len(u), np.ones(len(u)) / len(u), size=nb)
    T = W @ s; n = T[:, 0]
    r = (T[:, 5] - T[:, 1] * T[:, 2] / n) / np.sqrt((T[:, 3] - T[:, 1] ** 2 / n) * (T[:, 4] - T[:, 2] ** 2 / n))
    return float(np.corrcoef(p, y)[0, 1]), *np.nanpercentile(r, [2.5, 97.5])

rows = []
CFG = {'Abell': ('act_max', 'malinois_SKNSH'), 'Ulirsch': ('act_K562', 'malinois_K562')}
for ds, (actcol, sk) in CFG.items():
    t = d[d.dataset.eq(ds)].dropna(subset=[actcol, 'effect']).copy()
    for score in ['enformer_DNASE_mean', sk]:
        r0, lo0, hi0 = blockboot(t[score].to_numpy(), t.effect.to_numpy(), t.cluster.to_numpy())
        rows.append(dict(dataset=ds, score=score, selection='all', retained=len(t),
                         effect_sd=float(t.effect.std()), pearson=r0, ci_lo=lo0, ci_hi=hi0))
        for frac in [0.5, 0.25, 0.1]:
            s = t[t[actcol] >= t[actcol].quantile(1 - frac)]
            r, lo, hi = blockboot(s[score].to_numpy(), s.effect.to_numpy(), s.cluster.to_numpy())
            rows.append(dict(dataset=ds, score=score, selection=f'active_top{int(frac*100)}', retained=len(s),
                             effect_sd=float(s.effect.std()), pearson=r, ci_lo=lo, ci_hi=hi))
            # (a) random subsets of the same size
            rr = [np.corrcoef(*[t.sample(len(s), random_state=int(rng.integers(1e9)))[c].to_numpy()
                                for c in (score, 'effect')])[0, 1] for _ in range(200)]
            rows.append(dict(dataset=ds, score=score, selection=f'random_sizematched_top{int(frac*100)}',
                             retained=len(s), effect_sd=float(t.effect.std()), pearson=float(np.mean(rr)),
                             ci_lo=float(np.percentile(rr, 2.5)), ci_hi=float(np.percentile(rr, 97.5))))
            # (b) random subsets matched on the spread of the measured effect (|effect| strata)
            q = pd.qcut(t.effect.abs(), 20, labels=False, duplicates='drop')
            want = pd.qcut(s.effect.abs().reindex(t.index).dropna(), 20, labels=False, duplicates='drop') if False else None
            tgt = pd.Series(pd.qcut(t.effect.abs(), 20, labels=False, duplicates='drop'), index=t.index)
            need = tgt.loc[s.index].value_counts()
            rm = []
            for _ in range(200):
                idx = []
                for k, n in need.items():
                    pool = tgt[tgt.eq(k)].index
                    idx.extend(rng.choice(pool, min(n, len(pool)), replace=False))
                u = t.loc[idx]
                rm.append(np.corrcoef(u[score], u.effect)[0, 1])
            rows.append(dict(dataset=ds, score=score, selection=f'random_sdmatched_top{int(frac*100)}',
                             retained=len(s), effect_sd=float(np.nan), pearson=float(np.mean(rm)),
                             ci_lo=float(np.percentile(rm, 2.5)), ci_hi=float(np.percentile(rm, 97.5))))
res = pd.DataFrame(rows); res.to_csv(OUT / 'activity_control_results.csv', index=False)
pd.set_option('display.width', 220); print(res.round(3).to_string(index=False))
