"""Check the numbers quoted in the GLBIO 2027 manuscript against machine-generated result tables.
Run from the repository root after exploratory/diagnostics.py:  python exploratory/check_paper_numbers.py
Each check compares the value as rounded in the paper with the value in a results CSV."""
import json
from pathlib import Path
import pandas as pd

R = Path(__file__).resolve().parents[1]
V2 = R / 'reliability_audit_v2/results'; SG = R / 'subgroup_audit/results'; EX = R / 'exploratory/results'
fails = []
def check(label, value, paper, nd):
    ok = round(float(value), nd) == round(float(paper), nd)
    print(('OK  ' if ok else 'FAIL'), f'{label}: data {float(value):.{nd + 2}f} -> {round(float(value), nd)}, paper {paper}')
    if not ok: fails.append(label)

# Primary transfer test (Table 3): 95% descriptive and Bonferroni 97.5% decision intervals, Enformer
p = pd.read_csv(V2 / 'primary_bootstrap_comparisons.csv').set_index('dataset')
for ds, vals in {'Abell': (-2.1, -3.6, -0.7, -3.8, -0.5), 'Ulirsch': (-0.3, -6.5, 6.7, -7.5, 7.7)}.items():
    for col, pv in zip(['error_reduction', 'ci95_low', 'ci95_high', 'ci97_5_low', 'ci97_5_high'], vals):
        check(f'primary {ds} {col} (pp)', 100 * p.loc[ds, col], pv, 1)

# Point accuracy (Table 2)
pm = pd.read_csv(SG / 'point_model_estimates.csv')
for (ds, m), pv in {('Kircher', 'enformer'): 0.581, ('Kircher', 'borzoi'): 0.552, ('Abell', 'enformer'): 0.031,
                    ('Abell', 'borzoi'): 0.029, ('Ulirsch', 'enformer'): 0.095, ('Ulirsch', 'borzoi'): 0.098}.items():
    row = pm[pm.dataset.str.startswith(ds) & pm.model.eq(m) & pm.readout.eq('DNASE')]
    check(f'Pearson {ds} {m}', row.Pearson.iloc[0], pv, 3)

# Noise ceilings and corrected correlations (Table 6, Figure 2)
nc = pd.read_csv(SG / 'noise_ceiling_by_dataset.csv').set_index('dataset')
for ds, (ceil, ce, cb) in {'Abell': (0.667, 0.047, 0.043), 'Ulirsch': (0.699, 0.137, 0.140)}.items():
    check(f'noise ceiling {ds}', nc.loc[ds, 'Pearson_noise_ceiling'], ceil, 3)
    check(f'corrected Enformer {ds}', nc.loc[ds, 'enformer_noise_corrected_Pearson'], ce, 3)
    check(f'corrected Borzoi {ds}', nc.loc[ds, 'borzoi_noise_corrected_Pearson'], cb, 3)

# Kircher family-held-out gain over magnitude
st = pd.read_csv(SG / 'subgroup_tests.csv')
for m, vals in {'enformer': (1.7, -0.2, 5.1), 'borzoi': (1.2, -0.4, 4.2)}.items():
    row = st[st.dataset.eq('Kircher_LOFO') & st.model.eq(m) & st.analysis_mode.eq('global') & st.test_kind.eq('gain') & st.dimension.eq('all')].iloc[0]
    for col, pv in zip(['estimate', 'ci95_low', 'ci95_high'], vals):
        check(f'Kircher LOFO {m} {col} (pp)', 100 * row[col], pv, 1)

# Exploratory diagnostics (Table 4, Results text)
dg = pd.read_csv(EX / 'diagnostics_by_dataset.csv').set_index('dataset')
for ds, vals in {'Kircher': (0.864, 0.244, 0.317, 0.35, 0.61), 'Abell': (0.504, 0.491, -0.008, 0.03, 0.03),
                 'Ulirsch': (0.525, 0.480, 0.008, 0.03, 0.04)}.items():
    for col, pv, nd in zip(['balanced_accuracy', 'always_down_error', 'risk_error_spearman', 'within_r_all', 'within_r_active'], vals, [3, 3, 3, 2, 2]):
        check(f'{ds} {col}', dg.loc[ds, col], pv, nd)
check('Abell unfiltered sign error (%)', 100 * dg.loc['Abell', 'model_error'], 49.6, 1)
check('Kircher unfiltered sign error (%)', 100 * dg.loc['Kircher', 'model_error'], 13.0, 1)
ct = pd.read_csv(EX / 'source_target_contrasts.csv')
for (m, t), vals in {('risk_error_spearman', 'Abell'): (0.33, 0.25, 0.39), ('risk_error_spearman', 'Ulirsch'): (0.31, 0.22, 0.40),
                     ('median_within_r_active', 'Abell'): (0.58, 0.44, 0.70), ('median_within_r_active', 'Ulirsch'): (0.57, 0.38, 0.74)}.items():
    row = ct[ct.metric.eq(m) & ct.target.eq(t)].iloc[0]
    for col, pv in zip(['diff', 'lo', 'hi'], vals):
        check(f'contrast {m} {t} {col}', row[col], pv, 2)
cv = pd.read_csv(EX / 'simple_rule_curves.csv')
check('Kircher magnitude 25% error (%)', 100 * cv[(cv.dataset == 'Kircher') & (cv.rule == 'magnitude') & (cv.retain == 0.25)].sign_error.iloc[0], 6.1, 1)
print('\nFAILED: ' + ', '.join(fails) if fails else '\nAll checked manuscript numbers match the result tables.')
raise SystemExit(1 if fails else 0)
