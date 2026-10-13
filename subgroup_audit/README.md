# Subgroup Audit

**Original primary result: FAIL.** This follow-up cannot turn a selected subgroup
into success for the original Enformer protocol. All new model subgroup results
are pending on your existing Wahab scores. No GPU, model download or API call is
needed. AlphaGenome is excluded from loading, training, fitting and statistics.

## Run on Wahab / OnDemand

Unzip this package alongside the completed `Reliability_Audit_v2` directory.
Open an OnDemand terminal and run:

```bash
cd Subgroup_Audit
export CARDIAC_PYTHON=/absolute/path/to/your/v2/.venv/bin/python
export RELIABILITY_BASE_DIR=/absolute/path/to/completed/Reliability_Audit_v2
bash hpc/submit.sh
```

Use your existing v2 interpreter; its required CPU libraries should already be
installed. If needed install `requirements.txt` in that environment. The job
uses four CPUs,16GB host memory and a2-hour limit on Wahab's default CPU
partition. There is no GPU partition or GPU allocation. Set your existing
authorized account/QOS environment flags if your ODU allocation requires them.

No changes are made inside the original directory. The job reads both complete
eight-shard journals, the original external-policy export, data/protocol hashes,
selected track IDs and model metadata. It preserves the N-padding repair's code
provenance rather than editing old hashes. Results go into this new package's
`results/`. The original Enformer external ranks are reused exactly. A new
Borzoi policy and grouped source LOFO diagnostics are fitted on Kircher only.

The14 software tests include temporary synthetic end-to-end CPU analysis,
source-family isolation, label-independent selection, exact metadata budgets,
70 fixed test slots, Holm arithmetic, absent-stratum handling and noise-estimator
edge cases. Mock results are deleted and are not scientific evidence. The actual
Wahab model subgroup statistics have not been computed in the authoring sandbox.

If interrupted, rerun the same job. Partial per-study checkpoints are preserved,
but incomplete runs recompute deterministically rather than splice changed
inputs. Complete runs with the same identity exit without recomputing. Changes
to the plan, source scores or statistical assets require a separate run version.

## Read the preregistration first

`PREREGISTRATION.md` and `config/preregistered_plan.json` define the follow-up.
The plan was internally hash-locked before new subgroup calculations. The
aggregate FAIL was known. It is not external registration or a claim of blinded
confirmatory evidence. Statistical assets and aliases were locked before the
public-data noise calculation. The formal external family has **70 fixed slots**,
including absent/non-estimable contrasts. Holm correction at0.05 and conservative
family intervals are included. Source LOFO results are descriptive only.

Global selection has identical total candidate budgets, then is stratified.
Subgroup selected counts may differ; they are reported. A separate matched
budget sensitivity operates inside available regulatory/cell metadata strata.
Measured effect strength never sets selection budgets. External strong/weak
labels use replicate screening separate from held-out sign truth, with fixed
thresholds1.0 and0.5log2. A variant can have qualifying views in both bins; views
are not independent samples. An absent cell subgroup is NOT ESTIMABLE.

**Important count correction:** the reported4765/4726 and264/217 are retained
AUDITED sign calls. They should not be described as unmatched total retention.
The package separately checks both total candidate budgets and audited calls,
and refuses to change the original Enformer point estimate or call counts.

## Noise ceiling: real calculation already included

`results/public_data_noise_ceiling.csv` was actually computed from the verified
prepared public v2 variant table with10,000 genomic-cluster bootstrap draws,
after plan and asset locking. No model scores entered this calculation.

| Study | Approximate Pearson ceiling | Pointwise95% cluster interval |
|---|---:|---:|
| Abell |0.666878|[0.653190,0.679630]|
| Ulirsch |0.698980|[0.645163,0.739262]|
| Kircher |Not identifiable|Most variant-level replicate estimates unavailable|

This is the classical additive-error approximation, not a guarantee of a
specific model's achievable performance or a maximum sign accuracy. Abell uses
three paired transfections. Ulirsch includes both RNA and DNA variance but
unmodeled library covariance remains. The intervals condition on this
measurement-noise model. The CPU run recomputes these ceilings on the actual
verified Wahab population and adds paired corrected model correlations.

Kircher's two-library LDLR/SORT1 agreement is exported descriptively. It does not
provide a ceiling for all21 elements. No standard error is invented from p-values.

## Public promoter/enhancer annotation

A coordinate-matched subset is bundled for all62,057 unique v2 contrasts.
Classification is **genomic promoter-like/enhancer-like signature**, not a
claimed experimental MPRA construct type. ENCODE3 cCRE PLS overlaps are
promoter-like; pELS/dELS overlaps are enhancer-like. Other, ambiguous and
unannotated classes remain separate. Cell match uses actual selected DNASE
track identities and strict declared-cell aliases; it does not establish
matched culture/treatment or every other modality.

Verified public raw annotation:
https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb

Actual anonymously downloaded bytes:49,809,776.
SHA-256:`42a9cf761ed87e2ce7b26d60a84b9b623848a2a8c2de9c84c246eb4c6d38417c`.

Bundled subset SHA-256:
`8523e7b42b7a4706f5e12aebe726459843b9cbfbf40129f156723f6e2db120a8`.

Exact coordinate and source hashes are in `annotations/annotation_manifest.json`.
A changed universe is rejected, not silently labeled unannotated. Rebuilding the
subset requires only the original coordinate CSV, this frozen BigBed and optional
`pyBigWig==0.3.25`; statistical runs do not need pyBigWig or a new download.

## Outputs

- `subgroup_tests.csv`:70 corrected external tests plus source diagnostics;
  estimates, pointwise/family CIs, raw/Holm p-values, opportunity ratios/statuses.
- `retention_checks.csv`: global and metadata-stratum exact candidate budgets.
- `original_primary_reconciliation.csv`: original point estimates/counts unchanged.
- `point_model_estimates.csv`, `descriptive_point_accuracy_by_modality.csv`:
  frozen DNASE estimates and a non-oracle, descriptive weak-readout check.
- `noise_ceiling_by_dataset.csv`, `Kircher_two_library_agreement.csv`.
- `cell_match_provenance.csv`, `subgroup_population_coverage.csv`,
  `assay_effect_view_coverage.csv` (includes negligible/noisy views), source LOFO and
  Borzoi source-policy manifests, input/asset hashes and complete input identities.
- `annotated_audit_variants.csv.gz`, `follow_up_summary.json`, `SUMMARY.txt`.
- Subgroup forests, model/noise-ceiling plot and interactions, each PNG320dpi,
  vector PDF and SVG. Empty contrasts are not plotted as zero-effect estimates.

Null-centered cluster-bootstrap p-values are approximate conditional tests.
Holm's nominal FWER requires valid marginal tests. Shared assay effects,
interdependent genomic blocks or fixed-training uncertainty are not eliminated
by the correction. Non-significance is not equivalence. A favorable subgroup
cannot rescue the original global FAIL or establish clinical use.

## Full peer-reviewed conference target

**GLBIO2027 full-paper track. Deadline:January7,2027,23:59AoE.**
Use OUP **Modern Large** and author-year references, a structured abstract<=250
words and<=9 inclusive pages. `conference/Conference_Recommendation.md` gives
the verified official sources and the honest negative-result framing.

The compiled `conference/GLBIO2027_manuscript.pdf` rewrites the old workshop draft
as a conference-format content outline. Original primary numbers are labeled
user-reported, the noise calculation is real, and all new model subgroup panels
are visibly pending. `Conference_Format_Outline.md` gives an8-page full-paper
budget plus the9-page inclusive ceiling. No article has been submitted or accepted.

The OUP class and supporting template archive are unchanged and retain their
LPPL licenses. The author-year Modern Large options are journal-specific.
`update_outline.py` updates the abstract/result snippets after verified CPU
completion; then `bash conference/compile.sh` regenerates the PDF if a local
TeX installation is available. Otherwise upload the conference folder to
Overleaf. CPU results do not depend on TeX. Additional CTAN algorithm packages
are included solely to compile the public OUP class on a minimal TeX install.

Conference acceptance is peer-reviewed. Publication in the Bioinformatics
Advances special issue requires a separate subsequent journal review; conference
acceptance does not guarantee journal publication. Add a permanent public code
archive/DOI, contributing authors and truthful funding/disclosure statements
before submission. The9-page format must not be achieved by shrinking fonts.

No R21 document, original result or Mayo personnel/resource section is modified.
