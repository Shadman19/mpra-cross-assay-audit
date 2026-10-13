# Locked follow-up subgroup plan

This plan was saved before any new subgroup performance calculation. The
completed primary Enformer FAIL and the reported Abell/Ulirsch aggregate results
were already known. This is an internally timestamped follow-up, not an external
registry submission or a blinded confirmatory study. No subgroup can reverse the
original FAIL. The full machine-readable plan and hash are in `config/`.

## Population and unchanged forecasts

Use the completed v2 primary variant universe, exclusion rules, frozen Enformer
and Borzoi DNASE forecasts and recorded selected track IDs. Reuse the exact old
Enformer external risk ranks. Fit a new Borzoi risk policy on Kircher only using
the frozen original source split, hyperparameters and output-only features.
For source diagnostics, remove each held family before fitting, scaling or
calibration. Source diagnostics cannot establish external replication.

No AlphaGenome outputs are read, trained on, normalized or fitted. No foundation
model is executed. Original data, scores, classifiers and primary results are
read-only. The reported N-padding repair remains in the input provenance.

## Fixed strata

| Dimension | Rule fixed before performance inspection | Guard |
|---|---|---|
| Promoter vs enhancer | Uniform coordinate overlap with ENCODE3 cCRE: PLS = promoter-like; pELS/dELS = enhancer-like | These are genomic signatures, not MPRA construct labels. Other/unannotated/ambiguous remain distinct. |
| Strong vs weak measured effect | Strong active screen views: absolute effect >=1.0 log2. Weak active: >=0.5 and <1.0. | Screening and sign truth use disjoint replicates. No effect label sets selection budgets. A variant may have views in both bins. |
| Enformer vs Borzoi | Frozen DNASE point forecasts, their own magnitude rankings and source-only reliability policies | No oracle best modality. Other modality correlations descriptive only. |
| Matched vs unmatched cells | All actual selected DNASE track identities must match strict declared-cell aliases | Tissue/other-cell/cross-species proxies unmatched; generic unresolved identity unknown. No synthetic mismatched scores. Cell identity does not guarantee treatment/timepoint match. |

The main subgroup estimand preserves the original global selection, then
stratifies. Exactly the same total candidate budget is retained by policy and
magnitude within each original block. Selected counts inside a subgroup can
differ and are displayed, alongside audited sign calls.

A separately fixed sensitivity reselects within metadata strata for regulatory
signature and cell matching. It has exactly equal within-stratum budgets.
Measured strong/weak effects are never selection strata.

## Statistical family

| External test slots | Number |
|---|---:|
| Overall policy-vs-magnitude gains: 2 studies x 2 models | 4 |
| Global-selection subgroup gains: 2 studies x 2 models x 3 dimensions x 2 levels | 24 |
| Global subgroup interactions: 2 studies x 2 models x 3 dimensions | 12 |
| Within-metadata gains: 2 studies x 2 models x 2 dimensions x 2 levels | 16 |
| Within-metadata interactions: 2 studies x 2 models x 2 dimensions | 8 |
| Paired model contrasts: 2 studies x Pearson/sign error/filtering gain | 6 |
| **Fixed family** | **70** |

Resample paired 1 Mb clusters 10,000 times within each study. Use null-centered,
two-sided bootstrap tests. Apply Holm correction to all 70 slots at family
alpha 0.05. Non-estimable slots remain as p=1; do not remove them to weaken the
correction. Report pointwise 95% CIs and conservative Bonferroni family intervals
separately. The bootstrap tests are approximate conditional inference, not exact
randomization tests; dependence or invalid marginal tests can invalidate nominal
coverage/FWER. Genomic blocks are not proven independent LD units.

Require at least 20 audited calls per policy and 10 contributing clusters for
an external gain. Missing categories, unstable bootstrap distributions and
constant point forecasts are NOT ESTIMABLE. A significant interaction does not
prove a positive gain in either subgroup. A described positive gain also needs
at least 90% of the magnitude comparator's audited opportunities. Any favorable
subgroup still needs independent replication and cannot rescue the main FAIL.

## Noise ceiling

Estimate R = 1 - mean(mean-effect noise variance) / variance(mean effect), and
Pearson ceiling = sqrt(R). For Abell use variance of three paired effects /3.
For Ulirsch use RNA-ratio variance /6 plus DNA-ratio variance /2, explicitly
acknowledging unmodeled covariance. Bootstrap re-estimates R and corrected
Pearson with the same paired cluster draws. Nonpositive R is unidentified;
corrected absolute correlations >1 are unstable and never clipped. Do not apply
this correction to Spearman or claim a maximum sign accuracy.

The full Kircher ceiling is not identifiable from its processed table. Report
that limitation. LDLR/SORT1 two-library agreement is descriptive evidence for
those two elements, not a ceiling for the entire source benchmark. No p-value
is converted into an invented variant standard error.

## Interpretation and reproducibility

The existing Ulirsch interval allows appreciable gains and losses; a null
subgroup cannot prove equivalence. Any inferred benchmark-to-new-assay gap must
use comparable source and validation metrics, retain readout/assay limitations,
and avoid clinical or universal foundation-model failure claims.

`preregistered_plan.json` and `plan_lock.json` fix the design. Statistical code,
cell aliases and the public annotation subset were hash-locked before the
public-data noise calculation. New model subgroup statistics remain pending
until the Wahab CPU job. Distribution-wide hashes accompany the package.
