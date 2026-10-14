"""Exploratory diagnostics for the GLBIO 2027 paper (single consolidated script).

Consistency rules, matching the locked audit package (subgroup_audit/code):
  * Sign-error target: audit_engine.errors(). External studies use replicate screening
    columns to choose active variants and held-out replicate columns as sign truth;
    Kircher (no replicates) uses its measured effect for both, as in the package.
  * Kircher risks: the exported risk_enformer column, which the package produced with
    policies.source_lofo() (family-held-out fit + calibration split). This script
    recomputes source_lofo() and asserts equality before using it.
  * Selection and resampling unit: the package cluster (Kircher source family = fit_group;
    Abell/Ulirsch genomic_block). Kircher elements are carried with their family.
  * Bootstrap: B = 2000 cluster-multinomial draws everywhere, fixed seed. Intervals are
    descriptive percentile 95% intervals conditional on the frozen policies.
Model outputs select; measured effects only evaluate. No GPU, no new scoring.
Run from the repository root:  python exploratory/diagnostics.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'subgroup_audit' / 'code'))
from audit_engine import errors  # noqa: E402
from policies import source_lofo  # noqa: E402

OUT = ROOT / 'exploratory' / 'results'
OUT.mkdir(parents=True, exist_ok=True)
B = 2000
QS = [0.1, 0.25, 0.5, 0.75]
rng = np.random.default_rng(20271)

d = pd.read_csv(ROOT / 'subgroup_audit/results/annotated_audit_variants.csv.gz', low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true')].copy()
num = [c for c in d.columns if c.endswith(('_mean', '_std')) or c.startswith(('audit_', 'risk_')) or c in ('effect', 'original_enformer_risk')]
for c in num:
    d[c] = pd.to_numeric(d[c], errors='coerce')
d['cluster'] = np.where(d.dataset.eq('Kircher'), d.fit_group.astype(str), d.genomic_block.astype(str))
d['human_cell'] = ~d.cross_species_proxy.astype(str).str.lower().eq('true')
log = []
def say(s=''):
    print(s); log.append(s)

# ---- 0. Kircher risks must equal the package's own source_lofo() output -----------------
k_all = d[d.dataset.eq('Kircher')].copy()
recomputed, _ = source_lofo(k_all, 'enformer')
maxdiff = float(np.nanmax(np.abs(recomputed.to_numpy() - k_all.risk_enformer.to_numpy())))
assert maxdiff < 1e-9, f'Exported Kircher risk differs from source_lofo by {maxdiff}'
say(f'0. Kircher risk_enformer equals policies.source_lofo() (max abs diff {maxdiff:.1e})')

def prep(ds):
    g = d[d.dataset.eq(ds)].reset_index(drop=True).copy()
    g['err'] = errors(g, 'enformer', 'all', source=(ds == 'Kircher'))
    g['valid'] = np.isfinite(g.err)
    cl = sorted(g.cluster.unique()); g['ci'] = g.cluster.map({c: i for i, c in enumerate(cl)})
    W = np.random.default_rng(abs(hash(ds)) % 2**32 if False else {'Kircher': 11, 'Abell': 12, 'Ulirsch': 13}[ds]).multinomial(len(cl), np.ones(len(cl)) / len(cl), size=B).astype(float)
    return g, len(cl), W

def csum(g, v, n):  # per-cluster sums
    return np.bincount(g.ci.to_numpy(), weights=np.asarray(v, float), minlength=n)

def keep(g, score, q):  # package-style: floor(q * cluster size) lowest scores, tie-break row_id
    t = g.assign(_r=np.asarray(score, float)).sort_values(['ci', '_r', 'row_id'], kind='stable')
    k = np.zeros(len(g), bool)
    for _, h in t.groupby('ci', sort=False):
        m = int(np.floor(q * len(h)))
        if m: k[h.index[:m]] = True
    return k

def boot_spearman(x, y, ci, n, W):
    pt = spearmanr(x, y)[0]; bs = np.empty(B)
    idx = [np.where(ci == j)[0] for j in range(n)]
    for b in range(B):
        ii = np.concatenate([np.repeat(idx[j], int(W[b, j])) for j in range(n) if W[b, j] > 0])
        bs[b] = spearmanr(x[ii], y[ii])[0]
    return pt, bs

def median_within_r(g, mask, unit):
    rows = []
    for (cl, u), h in g[mask].groupby(['ci', unit]):
        if len(h) >= 10 and h.effect.std() > 0 and h.enformer_DNASE_mean.std() > 0:
            rows.append((cl, pearsonr(h.enformer_DNASE_mean, h.effect)[0]))
    return pd.DataFrame(rows, columns=['ci', 'r'])

def boot_median(rr, n, W):
    pt = float(rr.r.median()) if len(rr) else np.nan; bs = np.full(B, np.nan)
    by = [rr.r[rr.ci == j].to_numpy() for j in range(n)]
    for b in range(B):
        v = np.concatenate([np.tile(by[j], int(W[b, j])) for j in range(n) if W[b, j] > 0 and len(by[j])])
        if len(v): bs[b] = np.median(v)
    return pt, bs

def ci95(v):
    return np.nanpercentile(v, [2.5, 97.5])

S = {}
rule_rows, band_rows, curve_rows, diag_rows, sens_rows = [], [], [], [], []
for ds in ['Kircher', 'Abell', 'Ulirsch']:
    g, n, W = prep(ds)
    unit = 'element' if ds == 'Kircher' else 'genomic_block'
    v = g.valid.to_numpy(); e = np.where(v, g.err, 0.0)
    # ---- 1. learned risk vs observed error ------------------------------------------
    x = g.risk_enformer.to_numpy()[v]; y = g.err.to_numpy()[v]; ci = g.ci.to_numpy()[v]
    rho, rho_bs = boot_spearman(x, y, ci, n, W)
    q5 = pd.qcut(pd.Series(x).rank(method='first'), 5, labels=False)
    quint = pd.Series(y).groupby(q5).mean().round(3).tolist()
    # ---- 2. within-group accuracy ----------------------------------------------------
    r_all = median_within_r(g, np.ones(len(g), bool), unit)
    r_act = median_within_r(g, g.valid.to_numpy(), unit)
    m_all, m_all_bs = boot_median(r_all, n, W); m_act, m_act_bs = boot_median(r_act, n, W)
    # ---- 3. constant-sign baselines and balanced accuracy (same error machinery) ----
    def err_of(forecast):
        h = g.copy(); h['enformer_DNASE_mean'] = forecast
        return errors(h, 'enformer', 'all', source=(ds == 'Kircher'))
    up = err_of(np.ones(len(g))); down = err_of(-np.ones(len(g)))
    pos = g.enformer_DNASE_mean.to_numpy() > 0
    # class-conditional error: truth-positive share of a variant = 1 - error of 'always up'
    tp = 1 - up[v]; tn = up[v]; err_m = g.err.to_numpy()[v]
    wrong_pos = np.where(pos[v], 0, tp)          # truth up, model said down
    wrong_neg = np.where(pos[v], tn, 0)          # truth down, model said up
    bal = 1 - 0.5 * (wrong_pos.sum() / tp.sum() + wrong_neg.sum() / tn.sum())
    # ---- 4. simple fit-free rules vs magnitude ---------------------------------------
    ze = g.enformer_DNASE_mean / g.enformer_DNASE_mean.std(); zb = g.borzoi_DNASE_mean / g.borzoi_DNASE_mean.std()
    sd = g.enformer_DNASE_std / g.enformer_DNASE_std.std()
    rules = {'magnitude': -ze.abs(), 'low dispersion': sd, 'magnitude/dispersion': -ze.abs() / (sd + 1e-3),
             'cross-model agreement': (ze - zb).abs(),
             'consensus': -pd.Series(np.where(np.sign(ze) == np.sign(zb), (ze.abs() + zb.abs()) / 2, 0.0)),
             'learned policy': g.risk_enformer, 'random': pd.Series(np.random.default_rng(7).random(len(g)))}
    sums = {}
    for r, sc in rules.items():
        for q in QS:
            kk = keep(g, sc, q)
            sums[(r, q)] = (csum(g, e * (kk & v), n), csum(g, (kk & v).astype(float), n))
            curve_rows.append(dict(dataset=ds, rule=r, retain=q, sign_error=sums[(r, q)][0].sum() / sums[(r, q)][1].sum(), n_scored=int(sums[(r, q)][1].sum())))
    unf = e.sum() / v.sum()
    for r in rules:
        curve_rows.append(dict(dataset=ds, rule=r, retain=1.0, sign_error=unf, n_scored=int(v.sum())))
    names = [(r, q) for r in rules if r != 'magnitude' for q in QS]
    def gains(Wm):
        out = []
        for r, q in names:
            me, mn = Wm @ sums[('magnitude', q)][0], Wm @ sums[('magnitude', q)][1]
            ce, cn = Wm @ sums[(r, q)][0], Wm @ sums[(r, q)][1]
            out.append(me / mn - ce / cn)
        return np.array(out).T
    pt = gains(np.ones((1, n)))[0]; bs = gains(W)
    se = np.nanstd(bs, 0); crit = np.nanpercentile(np.nanmax(np.abs((bs - pt) / se), 1), 95)
    for j, (r, q) in enumerate(names):
        lo, hi = ci95(bs[:, j])
        band_rows.append(dict(dataset=ds, rule=r, retain=q, gain_pp=100 * pt[j], ci_lo=100 * lo, ci_hi=100 * hi,
                              sim_lo=100 * (pt[j] - crit * se[j]), sim_hi=100 * (pt[j] + crit * se[j]), sim_crit=crit))
    # ---- 5. feature support -----------------------------------------------------------
    S[ds] = dict(g=g, n=n, W=W, rho=rho, rho_bs=rho_bs, m_act=m_act, m_act_bs=m_act_bs, m_all=m_all, m_all_bs=m_all_bs)
    diag_rows.append(dict(dataset=ds, clusters=n, scored=int(v.sum()), model_error=err_m.mean(), always_up_error=up[v].mean(),
                          always_down_error=down[v].mean(), balanced_accuracy=bal, risk_error_spearman=rho,
                          rho_lo=ci95(rho_bs)[0], rho_hi=ci95(rho_bs)[1], risk_quintile_error=quint,
                          within_r_all=m_all, within_r_all_lo=ci95(m_all_bs)[0], within_r_all_hi=ci95(m_all_bs)[1], groups_all=len(r_all),
                          within_r_active=m_act, within_r_active_lo=ci95(m_act_bs)[0], within_r_active_hi=ci95(m_act_bs)[1], groups_active=len(r_act),
                          median_abs_pred=float(np.median(np.abs(g.enformer_DNASE_mean)))))

# ---- 6. source-minus-target contrasts (independent resampling per study) -------------
con_rows = []
for t in ['Abell', 'Ulirsch']:
    for key, lab in [('rho', 'risk_error_spearman'), ('m_act', 'median_within_r_active'), ('m_all', 'median_within_r_all')]:
        dif = S['Kircher'][key + '_bs'] - S[t][key + '_bs']; lo, hi = ci95(dif)
        con_rows.append(dict(metric=lab, target=t, source=S['Kircher'][key], target_value=S[t][key], diff=S['Kircher'][key] - S[t][key], lo=lo, hi=hi))
src_abs = np.abs(S['Kircher']['g'].enformer_DNASE_mean); lo5, hi95 = np.percentile(src_abs, [5, 95])
for t in ['Abell', 'Ulirsch']:
    a = np.abs(S[t]['g'].enformer_DNASE_mean)
    con_rows.append(dict(metric='share_outside_source_5_95_abs_pred', target=t, source=np.nan, target_value=float(((a < lo5) | (a > hi95)).mean()), diff=np.nan, lo=np.nan, hi=np.nan))

# ---- 7. Kircher sensitivity: human-cell libraries only -------------------------------
g, n, W = S['Kircher']['g'], S['Kircher']['n'], S['Kircher']['W']
for lab, msk in [('all Kircher', np.ones(len(g), bool)), ('human cells only', g.human_cell.to_numpy())]:
    h = g[msk]; v = h.valid.to_numpy()
    kk = keep(h.reset_index(drop=True), -h.enformer_DNASE_mean.abs().to_numpy(), 0.25)
    sens_rows.append(dict(subset=lab, variants=len(h), elements=h.element.nunique(), pearson=pearsonr(h.enformer_DNASE_mean, h.effect)[0],
                          median_within_element_r=median_within_r(h.reset_index(drop=True), np.ones(len(h), bool), 'element').r.median(),
                          unfiltered_error=h.err[v].mean(), magnitude_top25_error=h.err.to_numpy()[kk & v].mean(),
                          risk_error_spearman=spearmanr(h.risk_enformer[v], h.err[v])[0]))

pd.DataFrame(diag_rows).to_csv(OUT / 'diagnostics_by_dataset.csv', index=False)
pd.DataFrame(con_rows).to_csv(OUT / 'source_target_contrasts.csv', index=False)
pd.DataFrame(band_rows).to_csv(OUT / 'simple_rule_gains_bands.csv', index=False)
pd.DataFrame(curve_rows).to_csv(OUT / 'simple_rule_curves.csv', index=False)
pd.DataFrame(sens_rows).to_csv(OUT / 'kircher_human_cell_sensitivity.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
for name, df in [('diagnostics', pd.DataFrame(diag_rows)), ('contrasts', pd.DataFrame(con_rows)), ('sensitivity', pd.DataFrame(sens_rows))]:
    say(f'\n### {name}'); say(df.round(4).to_string())
bd = pd.DataFrame(band_rows); bd['sim_excl0'] = (bd.sim_lo > 0) | (bd.sim_hi < 0)
say('\n### simple rules at 25% retention'); say(bd[bd.retain == 0.25].round(2).to_string())
say('\n### simultaneous bands: comparisons excluding zero'); say(bd.groupby('dataset').agg(tests=('sim_excl0', 'size'), excl0=('sim_excl0', 'sum'), better=('sim_lo', lambda s: int((s > 0).sum())), crit=('sim_crit', 'first')).round(2).to_string())
cv = pd.DataFrame(curve_rows).pivot_table(index=['dataset', 'rule'], columns='retain', values='sign_error')
say('\n### sign error curves (%)'); say((100 * cv).round(1).to_string())
(OUT / 'diagnostics_log.txt').write_text('\n'.join(log) + '\n')
