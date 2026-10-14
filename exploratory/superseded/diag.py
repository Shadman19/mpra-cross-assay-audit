"""Exploratory diagnostics requested in review round 2. Model outputs select; measured effects only evaluate."""
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr
rng = np.random.default_rng(7)
d = pd.read_csv('results/annotated_audit_variants.csv.gz', low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true')]
for c in ['effect', 'enformer_DNASE_mean', 'enformer_DNASE_std', 'borzoi_DNASE_mean', 'risk_enformer']:
    d[c] = pd.to_numeric(d[c], errors='coerce')
d = d.dropna(subset=['effect', 'enformer_DNASE_mean', 'borzoi_DNASE_mean'])
d['grp'] = np.where(d.dataset.eq('Kircher'), d.element.astype(str), d.genomic_block.astype(str))
P = lambda s: print('\n### ' + s)

P('A. Pooled vs within-group Pearson (Enformer DNASE); groups = element (Kircher) or genomic block (others), n>=10')
rows = []
for ds, g in d.groupby('dataset'):
    for lab, h in [('all', g), ('active |effect|>=0.5', g[g.effect.abs() >= 0.5])]:
        pooled = pearsonr(h.enformer_DNASE_mean, h.effect)[0]
        w = [(len(x), pearsonr(x.enformer_DNASE_mean, x.effect)[0]) for _, x in h.groupby('grp') if len(x) >= 10 and x.effect.std() > 0 and x.enformer_DNASE_mean.std() > 0]
        n = np.array([a for a, _ in w]); r = np.array([b for _, b in w])
        fz = np.tanh(np.sum((n - 3) * np.arctanh(np.clip(r, -.999, .999))) / np.sum(n - 3)) if len(w) else np.nan
        rows.append(dict(dataset=ds, subset=lab, n=len(h), effect_sd=h.effect.std(), pooled_r=pooled, groups=len(w), median_within_r=np.median(r) if len(w) else np.nan, weighted_within_r=fz))
A = pd.DataFrame(rows); print(A.round(3).to_string())

P('B. Sign confusion among active variants (|effect|>=0.5), Enformer DNASE')
rows = []
for ds, g in d.groupby('dataset'):
    h = g[g.effect.abs() >= 0.5]; yp = h.effect > 0; pp = h.enformer_DNASE_mean > 0
    rows.append(dict(dataset=ds, active=len(h), meas_up=int(yp.sum()), meas_down=int((~yp).sum()), pred_up_meas_up=int((pp & yp).sum()), pred_up_meas_down=int((pp & ~yp).sum()), pred_down_meas_up=int((~pp & yp).sum()), pred_down_meas_down=int((~pp & ~yp).sum()),
                     model_error=float((pp != yp).mean()), always_up_error=float((~yp).mean()), always_down_error=float(yp.mean()), balanced_acc=0.5 * ((pp & yp).sum() / max(yp.sum(), 1) + (~pp & ~yp).sum() / max((~yp).sum(), 1))))
B = pd.DataFrame(rows); print(B.round(3).to_string())

P('C. Simultaneous 95% bands (max-|z| bootstrap) for gain vs magnitude over rules x retention {10,25,50,75}')
QS = [0.1, 0.25, 0.5, 0.75]
res = []
for ds, g in d.groupby('dataset'):
    g = g.reset_index(drop=True)
    ze = g.enformer_DNASE_mean / g.enformer_DNASE_mean.std(); zb = g.borzoi_DNASE_mean / g.borzoi_DNASE_mean.std(); sd = g.enformer_DNASE_std / g.enformer_DNASE_std.std()
    rules = {'magnitude': ze.abs(), 'low dispersion': -sd, 'magnitude/dispersion': ze.abs() / (sd + 1e-3), 'cross-model agreement': -(ze - zb).abs(),
             'consensus': pd.Series(np.where(np.sign(ze) == np.sign(zb), (ze.abs() + zb.abs()) / 2, 0.0)), 'learned policy': -g.risk_enformer, 'random': pd.Series(rng.random(len(g)))}
    act = (g.effect.abs() >= 0.5).values; W = ((np.sign(g.enformer_DNASE_mean) != np.sign(g.effect)).values & act).astype(float)
    gi = g.groupby('grp').indices; blocks = np.array(list(gi))
    def keep(sc, q):
        k = np.zeros(len(g), bool); v = np.asarray(sc, float)
        for idx in gi.values():
            m = int(np.floor(q * len(idx)))
            if m: k[idx[np.argsort(-v[idx], kind='stable')[:m]]] = True
        return k & act
    K = {(r, q): keep(s, q) for r, s in rules.items() for q in QS}
    names = [(r, q) for r in rules if r != 'magnitude' for q in QS]
    def gains(sel):
        ii = np.concatenate([gi[b] for b in sel]); out = []
        for r, q in names:
            kb, kc = K[('magnitude', q)][ii], K[(r, q)][ii]
            out.append((W[ii][kb].mean() if kb.any() else np.nan) - (W[ii][kc].mean() if kc.any() else np.nan))
        return np.array(out)
    pt = gains(blocks); B_ = np.array([gains(rng.choice(blocks, len(blocks))) for _ in range(1000)])
    se = np.nanstd(B_, 0); z = np.nanmax(np.abs((B_ - pt) / se), 1); c = np.nanpercentile(z, 95)
    for (r, q), p, s in zip(names, pt, se):
        res.append(dict(dataset=ds, rule=r, retain=q, gain_pp=100 * p, sim_lo=100 * (p - c * s), sim_hi=100 * (p + c * s), crit=c))
C = pd.DataFrame(res); C['excl0'] = (C.sim_lo > 0) | (C.sim_hi < 0)
print(C[C.excl0].round(2).to_string()); print(C.groupby('dataset').agg(tests=('excl0', 'size'), exclude_zero=('excl0', 'sum'), crit=('crit', 'first')).round(2).to_string())
print(C[C.retain == 0.25].round(2).to_string())

P('D. Learned-policy risk vs observed sign error (active variants), risk quintiles within dataset')
rows = []
for ds, g in d.groupby('dataset'):
    h = g[g.effect.abs() >= 0.5].dropna(subset=['risk_enformer']).copy()
    h['err'] = (np.sign(h.enformer_DNASE_mean) != np.sign(h.effect)).astype(float); h['q'] = pd.qcut(h.risk_enformer.rank(method='first'), 5, labels=False)
    e = h.groupby('q').agg(mean_risk=('risk_enformer', 'mean'), observed_error=('err', 'mean'), n=('err', 'size'))
    rho = spearmanr(h.risk_enformer, h.err)[0]
    print(ds, 'spearman(risk,error)=%.3f' % rho); print(e.round(3).to_string())

P('E. Feature support: share of target variants inside source (Kircher) 5-95% range')
src = d[d.dataset.eq('Kircher')]
feats = {'|Enformer DNASE|': lambda x: x.enformer_DNASE_mean.abs(), 'Enformer track SD': lambda x: x.enformer_DNASE_std, '|Enformer-Borzoi|': lambda x: (x.enformer_DNASE_mean - x.borzoi_DNASE_mean).abs()}
rows = []
for ds in ['Abell', 'Ulirsch']:
    t = d[d.dataset.eq(ds)]
    for f, fn in feats.items():
        s, v = fn(src), fn(t); lo, hi = np.nanpercentile(s, [5, 95])
        rows.append(dict(target=ds, feature=f, source_median=np.nanmedian(s), target_median=np.nanmedian(v), target_inside_source_5_95=float(((v >= lo) & (v <= hi)).mean())))
print(pd.DataFrame(rows).round(4).to_string())
A.to_csv('results/diag_within_group_r.csv', index=False); B.to_csv('results/diag_sign_confusion.csv', index=False); C.to_csv('results/diag_simultaneous_bands.csv', index=False)
