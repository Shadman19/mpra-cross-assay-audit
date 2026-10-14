"""Round-4 check: direct source-target contrasts with bootstrap CIs (independent resampling per study).
(a) Spearman(risk, sign error) among active variants: Kircher LOFO risk vs Abell/Ulirsch frozen-policy risk.
(b) Median within-group Pearson r (Enformer DNASE vs effect): Kircher (elements) vs targets (blocks)."""
import sys, numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier
sys.path.insert(0, '/home/spath004/cardiac_r21/audit2/Reliability_Audit_v2/code')
import audit
from audit import deployment_features
rng = np.random.default_rng(404); B = 2000
d = pd.read_csv('results/annotated_audit_variants.csv.gz', low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true')].copy()
for c in ['effect', 'enformer_DNASE_mean', 'borzoi_DNASE_mean', 'risk_enformer']: d[c] = pd.to_numeric(d[c], errors='coerce')
d = d.dropna(subset=['effect', 'enformer_DNASE_mean', 'borzoi_DNASE_mean']).reset_index(drop=True)

k = d[d.dataset.eq('Kircher')].reset_index(drop=True)
X, _, _ = deployment_features(k)
act = (k.effect.abs() >= 0.5).values
y = (np.sign(k.enformer_DNASE_mean) != np.sign(k.effect)).values
fam = k.fit_group.astype(str).values
lofo = np.full(len(k), np.nan)
for f in np.unique(fam):
    tr = (fam != f) & act; te = fam == f
    if tr.sum() < 50 or len(np.unique(y[tr])) < 2: continue
    m = HistGradientBoostingClassifier(max_iter=120, max_leaf_nodes=7, min_samples_leaf=30, l2_regularization=5, learning_rate=.07, random_state=audit.SEED)
    m.fit(X[tr], y[tr], sample_weight=audit.balanced(fam[tr]))  # exact frozen recipe
    lofo[te] = m.predict_proba(X[te])[:, 1]
src = pd.DataFrame({'risk': lofo, 'err': y.astype(float), 'g': fam})[act].dropna()

def boot_rho(t):
    gi = t.groupby('g').indices; u = np.array(list(gi)); r, e = t.risk.values, t.err.values
    pt = spearmanr(r, e)[0]; bs = []
    for _ in range(B):
        ii = np.concatenate([gi[x] for x in rng.choice(u, len(u))]); bs.append(spearmanr(r[ii], e[ii])[0])
    return pt, np.array(bs)

rows = []
s_pt, s_bs = boot_rho(src)
print('LOFO risk-error Spearman (frozen recipe): %.3f' % spearmanr(src.risk, src.err)[0])
print('source family-bootstrap 95%% CI %.3f to %.3f' % tuple(np.nanpercentile(s_bs, [2.5, 97.5])))
src['q'] = pd.qcut(src.risk.rank(method='first'), 5, labels=False)
print(src.groupby('q').err.mean().round(3).to_string())
print('### (a) risk-error Spearman, source minus target (Kircher LOFO, family bootstrap; target block bootstrap)')
for ds in ['Abell', 'Ulirsch']:
    g = d[d.dataset.eq(ds) & (d.effect.abs() >= 0.5)].dropna(subset=['risk_enformer'])
    t = pd.DataFrame({'risk': g.risk_enformer.values, 'err': (np.sign(g.enformer_DNASE_mean) != np.sign(g.effect)).astype(float).values, 'g': g.genomic_block.astype(str).values})
    t_pt, t_bs = boot_rho(t); diff = s_bs - t_bs
    lo, hi = np.nanpercentile(diff, [2.5, 97.5])
    print('%-8s source=%.3f target=%.3f diff=%.3f 95%% CI %.3f to %.3f' % (ds, s_pt, t_pt, s_pt - t_pt, lo, hi))
    rows.append(dict(metric='risk_error_spearman', target=ds, source=s_pt, target_value=t_pt, diff=s_pt - t_pt, lo=lo, hi=hi))

print('\n### (b) median within-group Pearson r, source minus target (group bootstrap)')
d['grp'] = np.where(d.dataset.eq('Kircher'), d.element.astype(str), d.genomic_block.astype(str))
def wr(g):
    return np.array([pearsonr(x.enformer_DNASE_mean, x.effect)[0] for _, x in g.groupby('grp') if len(x) >= 10 and x.effect.std() > 0 and x.enformer_DNASE_mean.std() > 0])
for lab, flt in [('all', lambda x: x), ('active', lambda x: x[x.effect.abs() >= 0.5])]:
    ws = wr(flt(d[d.dataset.eq('Kircher')])); sb = np.array([np.median(rng.choice(ws, len(ws))) for _ in range(B)])
    for ds in ['Abell', 'Ulirsch']:
        wt = wr(flt(d[d.dataset.eq(ds)])); tb = np.array([np.median(rng.choice(wt, len(wt))) for _ in range(B)])
        diff = sb - tb; lo, hi = np.percentile(diff, [2.5, 97.5]); pt = np.median(ws) - np.median(wt)
        print('%-6s %-8s source=%.3f (%d) target=%.3f (%d) diff=%.3f 95%% CI %.3f to %.3f' % (lab, ds, np.median(ws), len(ws), np.median(wt), len(wt), pt, lo, hi))
        rows.append(dict(metric='median_within_r_' + lab, target=ds, source=np.median(ws), target_value=np.median(wt), diff=pt, lo=lo, hi=hi))
pd.DataFrame(rows).to_csv('results/diag_source_target_contrasts.csv', index=False)
