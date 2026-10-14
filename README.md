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
- `exploratory/` - post hoc analyses reported as exploratory in the paper (one script, `diagnostics.py`), and
  `check_paper_numbers.py`, which checks every manuscript number against the result tables.

## Data
All inputs are public (GEO GSE126550, GSE174534, GSE87711; UCSC hg38 and liftover chain). Raw files are not
redistributed; `reliability_audit_v2/code/download_reference.py` and `prepare_data.py` download them and verify the
SHA-256 values listed in `reliability_audit_v2/README.md`. Model weights are downloaded from their public releases.

## Reproduce
CPU, from the repository root: `bash reproduce.sh`. It installs the pinned environment, regenerates all exploratory
tables and checks every manuscript number (primary, follow-up and exploratory) against the result files.

GPU scoring and the original CPU audits are described in `reliability_audit_v2/README.md` and
`subgroup_audit/README.md`. Those two READMEs and `reliability_audit_v2/FILE_SHA256.json` are the original package
documents, kept unchanged for provenance. They mention items not redistributed here: the working `preprint/` and
`conference/` folders (superseded by the GLBIO 2027 manuscript) and an optional AlphaGenome module
(`code/alpha_optional.py`, `hpc/alpha_cpu.slurm`), which was removed. Absolute Wahab paths in
`subgroup_audit/results/execution.json` are provenance records, not requirements.

## Note
An optional AlphaGenome module present in an internal version was removed; no AlphaGenome outputs are used
anywhere in this repository or the paper.

## License
MIT (see `reliability_audit_v2/LICENSE`); third-party notices in `subgroup_audit/THIRD_PARTY_NOTICES.md`.
