"""Exploratory: simple label-free confidence rules vs magnitude, matched per-block retention.
Selection uses model outputs only; evaluation uses measured effect sign among active variants."""
import numpy as np, pandas as pd, sys
rng = np.random.default_rng(2027)
d = pd.read_csv('results/annotated_audit_variants.csv.gz', low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true')] if 'primary' in d else d
num = lambda c: pd.to_numeric(d[c], errors='coerce')
for c in ['effect', 'enformer_DNASE_mean', 'enformer_DNASE_std', 'borzoi_DNASE_mean', 'borzoi_DNASE_std', 'risk_enformer', 'risk_borzoi']:
    d[c] = num(c)
d = d.dropna(subset=['effect', 'enformer_DNASE_mean', 'borzoi_DNASE_mean'])
d['block'] = d['dataset'] + '|' + d['genomic_block'].astype(str)
out, curves = [], []
QS = [0.1, 0.25, 0.5, 0.75, 1.0]
for ds, g in d.groupby('dataset'):
    g = g.copy()
    ze = g.enformer_DNASE_mean / g.enformer_DNASE_mean.std(); zb = g.borzoi_DNASE_mean / g.borzoi_DNASE_mean.std()
    rules = {
        'magnitude (baseline)': ze.abs(),
        'low track dispersion': -(g.enformer_DNASE_std / g.enformer_DNASE_std.std()),
        'magnitude / dispersion': ze.abs() / (g.enformer_DNASE_std / g.enformer_DNASE_std.std() + 1e-3),
        'cross-model agreement': -(ze - zb).abs(),
        'consensus (signs agree)': np.where(np.sign(ze) == np.sign(zb), (ze.abs() + zb.abs()) / 2, 0.0),
        'learned policy (HGB)': -g.risk_enformer if g.risk_enformer.notna().all() else None,
        'random': pd.Series(rng.random(len(g)), index=g.index),
    }
    active = g.effect.abs() >= 0.5
    wrong = (np.sign(g.enformer_DNASE_mean) != np.sign(g.effect)) & active
    blocks = g.block.values
    def err(score, q):
        s = pd.Series(np.asarray(score, float), index=g.index)
        keep = np.zeros(len(g), bool)
        for b, idx in g.groupby('block').indices.items():
            k = int(np.floor(q * len(idx))) if q < 1 else len(idx)
            if k: keep[idx[np.argsort(-s.values[idx], kind='stable')[:k]]] = True
        return keep
    keeps = {r: {q: err(sc, q) for q in QS} for r, sc in rules.items() if sc is not None}
    ub = np.unique(blocks); bidx = {b: np.where(blocks == b)[0] for b in ub}
    for r in keeps:
        for q in QS:
            k = keeps[r][q] & active.values
            curves.append(dict(dataset=ds, rule=r, retain=q, n_active=int(k.sum()), sign_error=float(wrong.values[k].mean()) if k.sum() else np.nan))
    # paired block bootstrap at q=0.25 vs magnitude
    base = keeps['magnitude (baseline)'][0.25] & active.values
    W = wrong.values.astype(float)
    for r in keeps:
        if r.startswith('magnitude (b'): continue
        cand = keeps[r][0.25] & active.values
        def gain(sel):
            ii = np.concatenate([bidx[b] for b in sel])
            eb, ec = base[ii], cand[ii]
            return (W[ii][eb].mean() if eb.sum() else np.nan) - (W[ii][ec].mean() if ec.sum() else np.nan)
        pt = gain(ub)
        bs = [gain(rng.choice(ub, len(ub))) for _ in range(2000)]
        out.append(dict(dataset=ds, rule=r, gain_pp=100 * pt, lo=100 * np.nanpercentile(bs, 2.5), hi=100 * np.nanpercentile(bs, 97.5), n_blocks=len(ub)))
pd.DataFrame(curves).to_csv('results/exploratory_simple_rule_curves.csv', index=False)
o = pd.DataFrame(out); o.to_csv('results/exploratory_simple_rule_gains.csv', index=False)
pd.set_option('display.width', 200)
print(o.round(2).to_string())
c = pd.DataFrame(curves).pivot_table(index=['dataset', 'rule'], columns='retain', values='sign_error')
print((100 * c).round(1).to_string())
