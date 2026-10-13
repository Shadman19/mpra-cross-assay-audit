# Cross-assay audit of Enformer and Borzoi regulatory-variant predictions

Code, frozen analysis plans and results for the manuscript
"Cross-assay performance and confidence-filtering limits of Enformer and Borzoi regulatory-variant predictions"
(Shadman Mahmood Khan Pathan, Old Dominion University).

## Contents
- `reliability_audit_v2/` - primary audit. GPU scoring of public MPRA variants with released Enformer and Borzoi weights,
  Kircher-only reliability policies, frozen transfer test on Abell 2022 and Ulirsch 2016. Includes the
  chromosome-boundary N-padding fix in `code/score_open_models.py`.
- `subgroup_audit/` - prespecified 70-test follow-up (Holm), noise ceilings, `PREREGISTRATION.md` and logged
  protocol deviations in `DEVIATIONS.md`.
- `exploratory/` - post hoc analyses reported as exploratory in the paper: simple confidence rules, risk-retention
  curves, simultaneous bootstrap bands, within-group correlations, sign confusion, family-held-out risk diagnostics.

## Data
All inputs are public (GEO GSE126550, GSE174534, GSE87711; UCSC hg38 and liftover chain). Raw files are not
redistributed; `reliability_audit_v2/code/download_reference.py` and `prepare_data.py` download them and verify the
SHA-256 values listed in `reliability_audit_v2/README.md`. Model weights are downloaded from their public releases.

## Reproduce
See `reliability_audit_v2/README.md` (GPU scoring, Slurm) and `subgroup_audit/README.md` (CPU follow-up).
Exploratory scripts run from `subgroup_audit/` with the same Python environment, e.g.
`python ../exploratory/simple_rules.py`.

## Note
An optional AlphaGenome module present in an internal version was removed; no AlphaGenome outputs are used
anywhere in this repository or the paper.

## License
MIT (see `reliability_audit_v2/LICENSE`); third-party notices in `subgroup_audit/THIRD_PARTY_NOTICES.md`.
