"""Round-3 checks: (1) family-held-out Kircher risk vs error; (2) bootstrap CIs for within-group correlation summaries."""
import sys, numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier
sys.path.insert(0, '/home/spath004/cardiac_r21/audit2/Reliability_Audit_v2/code')
from audit import deployment_features
rng = np.random.default_rng(11)
d = pd.read_csv('results/annotated_audit_variants.csv.gz', low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true')].copy()
for c in ['effect', 'enformer_DNASE_mean', 'borzoi_DNASE_mean', 'risk_enformer']: d[c] = pd.to_numeric(d[c], errors='coerce')
d = d.dropna(subset=['effect', 'enformer_DNASE_mean', 'borzoi_DNASE_mean']).reset_index(drop=True)

print('### 1. Family-held-out Kircher risk (HGB refit per held-out family on the package feature set)')
k = d[d.dataset.eq('Kircher')].reset_index(drop=True)
X, _, _ = deployment_features(k)
act = (k.effect.abs() >= 0.5).values
y = (np.sign(k.enformer_DNASE_mean) != np.sign(k.effect)).values
fam = k.fit_group.astype(str).values
lofo = np.full(len(k), np.nan)
for f in np.unique(fam):
    tr = (fam != f) & act; te = fam == f
    if tr.sum() < 50 or len(np.unique(y[tr])) < 2: continue
    m = HistGradientBoostingClassifier(max_depth=3, max_iter=200, learning_rate=0.05, random_state=0).fit(X[tr], y[tr])
    lofo[te] = m.predict_proba(X[te])[:, 1]
h = pd.DataFrame({'risk': lofo, 'err': y.astype(float), 'act': act, 'fam': fam}).query('act').dropna()
rho = spearmanr(h.risk, h.err)[0]
h['q'] = pd.qcut(h.risk.rank(method='first'), 5, labels=False)
print('families', h.fam.nunique(), 'n', len(h), 'spearman(LOFO risk, error)=%.3f' % rho)
print(h.groupby('q').agg(mean_risk=('risk', 'mean'), observed_error=('err', 'mean'), n=('err', 'size')).round(3).to_string())
fb = h.fam.unique(); bs = []
for _ in range(1000):
    s = rng.choice(fb, len(fb)); hh = pd.concat([h[h.fam == x] for x in s]); bs.append(spearmanr(hh.risk, hh.err)[0])
print('family-bootstrap 95%% CI: %.3f to %.3f' % tuple(np.nanpercentile(bs, [2.5, 97.5])))
for ds in ['Abell', 'Ulirsch']:
    g = d[d.dataset.eq(ds) & (d.effect.abs() >= 0.5)]
    e = (np.sign(g.enformer_DNASE_mean) != np.sign(g.effect)).astype(float); r0 = spearmanr(g.risk_enformer, e)[0]
    gb = g.genomic_block.astype(str).values; ub = np.unique(gb); idx = {b: np.where(gb == b)[0] for b in ub}; bs = []
    for _ in range(1000):
        ii = np.concatenate([idx[b] for b in rng.choice(ub, len(ub))]); bs.append(spearmanr(g.risk_enformer.values[ii], e.values[ii])[0])
    print(ds, 'spearman(risk,error)=%.3f, block-bootstrap 95%% CI %.3f to %.3f' % ((r0,) + tuple(np.nanpercentile(bs, [2.5, 97.5]))))

print('\n### 2. Within-group correlation summaries with group-bootstrap 95% CIs (groups n>=10, nonzero variance)')
d['grp'] = np.where(d.dataset.eq('Kircher'), d.element.astype(str), d.genomic_block.astype(str))
for ds, g in d.groupby('dataset'):
    for lab, hh in [('all', g), ('active', g[g.effect.abs() >= 0.5])]:
        w = [pearsonr(x.enformer_DNASE_mean, x.effect)[0] for _, x in hh.groupby('grp') if len(x) >= 10 and x.effect.std() > 0 and x.enformer_DNASE_mean.std() > 0]
        w = np.array(w); bs = [np.median(rng.choice(w, len(w))) for _ in range(2000)]
        print('%-8s %-6s groups=%3d median r=%.3f  95%% CI %.3f to %.3f' % ((ds, lab, len(w), np.median(w)) + tuple(np.percentile(bs, [2.5, 97.5]))))
