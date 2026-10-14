# Exploratory diagnostics

`diagnostics.py` produces every exploratory number in the manuscript (Table 4, Table 5, Figure 1 and the
source-minus-target contrasts). It follows the locked audit package exactly:

- sign-error target: `subgroup_audit/code/audit_engine.errors()` (replicate screen and held-out replicate truth on
  Abell and Ulirsch; measured effect on Kircher, which has no replicate estimates);
- Kircher risks: the exported family-held-out scores, which the script recomputes with `policies.source_lofo()`
  and asserts are identical before use;
- selection and resampling unit: the audit cluster (Kircher source family, otherwise genomic block);
  Kircher elements are resampled with their family;
- 2,000 cluster-bootstrap draws; percentile 95% intervals conditional on the frozen policies; simultaneous
  max-|z| bands over each study's 24 rule-by-retention comparisons.

`check_paper_numbers.py` compares 58 selected values quoted in the paper (hard-coded in the script) with the result tables and exits non-zero on
any mismatch. Run both with `bash reproduce.sh` from the repository root (pinned versions in `requirements.txt`;
exact Kircher-risk reproduction needs these versions).

`superseded/` keeps earlier scripts and outputs for transparency. They used a pooled-mean error definition,
element-level units and, in `diag2.py`/`diag3.py`, refits that did not follow the package's fit/calibration split.
Their numbers are not used in the manuscript.

Zero forecasts: none of the audited variants has an Enformer DNASE prediction of exactly zero in any study, so the error function (zero counts as wrong) and the balanced-accuracy convention (zero counts as down) cannot disagree.
