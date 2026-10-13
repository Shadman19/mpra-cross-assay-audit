# Reliability Audit v2

**PENDING: no new Enformer/Borzoi efficacy result has been generated.** Public
assay preparation and software checks have actually run. Full released-weight
inference, its GPU reproducibility test, and external reliability verdicts run on
ODU Wahab. Prior AlphaGenome results are author-reported, explicitly separated.

This package studies **deployment-time reliability of regulatory variant-effect
models**. It makes no cardiac, 3D-genome, clinical or grant-fundability claim.

## Five-line run guide

1. Upload/unzip on Wahab in OnDemand, open a terminal, and `cd Reliability_Audit_v2`.
2. Run `bash hpc/setup_env.sh`, then `export CARDIAC_PYTHON="$PWD/.venv/bin/python"`.
3. Run `bash hpc/submit.sh` to launch the `gpu` array and its dependent CPU audit.
4. Check `squeue` and `logs/`; rerun `bash hpc/submit.sh` for exit-75/partial shards.
5. Read `results/SUMMARY.txt` and `preprint/Workshop_Preprint.pdf`; AlphaGenome fitting stays gated.

## Actual evidence delivered

All six raw resources below were downloaded anonymously, parsed and checksum
verified. Full UCSC hg38 2bit was downloaded separately and checked against its
publisher MD5. Prepared CSVs are included. The genome and model weights are
large, reproducibly downloadable resources, not embedded in the ZIP. Data QC:

| Study | Acquired contrasts | Primary unique SNVs | Primary 1 Mb blocks |
|---|---:|---:|---:|
| Kircher 2019 | 44,647 | 28,932 | 19 |
| Ulirsch 2016 | 2,433 | 2,375 | 79 |
| Abell 2022 | 26,959 | 22,037 | 525 |

Kircher contains 21 displayed elements / 29 conditions, grouped into 19 source
families; raw GP1BA maps to GP1BB promoter coordinates. HBG2 is not a separate
source element. Low-quality contrasts and single-base deletions are acquired;
only preregistered canonical quality-filtered SNVs train the primary policy.
All mapped acquisition contrasts receive actual network passes. Separate
conditions may share a genomic contrast and a prediction pass, never a label.

The real-data tests check RNA/DNA effect and variance reconstruction,
independent locus exclusions, one canonical construct per Abell variant,
source/metadata integrity, track-strand permutations, exact retention,
external-label poisoning, and the closed AlphaGenome fitting gate. A mock
end-to-end audit writes its CSVs/figures only inside a temporary directory and
is deleted. These are software validations, not evidence for improved accuracy.
A tiny untrained Enformer native forward pass was also checked in the authoring
sandbox. Released-weight GPU output contracts are checked on Wahab before any
scientific scoring. Read `results/validation_report.json` for precise status.

## Exact public downloads and observed checksums

| Resource / accession | Anonymous download URL | SHA-256 of actual bytes |
|---|---|---|
| data/raw/kircher_elements.tsv.gz / GSE126550 | https://raw.githubusercontent.com/kircherlab/MPRA_SaturationMutagenesis/05d2ffb965090d3f5dd27dfb038cec493a15ab35/data/elements.tsv.gz | `fec2eed91fe27af3aae07ebce2eca65e9bad4bb6abba5d8c27f478887dd7b134` |
| data/raw/ulirsch_counts.tsv.gz / GSE87711 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87711/suppl/GSE87711_RBC_MPRA_minP_raw.txt.gz | `94f74f5de695f5389c6c9f42e163f8abdeb3544c1309145a9fbd44fa74f462fb` |
| data/raw/ulirsch_design.tsv.gz / GSE87711 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87711/suppl/GSE87711_design.txt.gz | `141eafcb023c9555dc129b11e5cf45985930d8cce2acb5193f511795a2750114` |
| data/raw/abell_counts.tsv.gz / GSE174534 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174534/suppl/GSE174534_OligoCounts.txt.gz | `576ffc492a455528c2fa2cd48de08bf9d4516ed4f8ebf10356bf3368a5723de4` |
| data/raw/abell_hg_conversion.tsv / GSE174534 | https://raw.githubusercontent.com/nsabell/mpra-v2/d62fca86f10b61ab8dbc64f2face45327ae040d0/annotation/hgConversionTable.txt | `12916af586dad818ac3292f580334f4f2eac926272a5a1d2f2c72282de122c20` |
| data/raw/hg19ToHg38.over.chain.gz / UCSC chain | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz | `5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d` |

Raw bytes are bundled under `data/raw/`; `data/source_manifest.json` records sizes.
`prepare_data.py` reconstructs Kircher directly from its raw table and derives
validation effects from count data. Ulirsch's design file is also included.

**Genome:** https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
(835,393,456 bytes; publisher MD5 `dcc3ea27079aa6dc3f9deccd7275e0f8`).
Publisher checksum list: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt.
The exact list downloaded is under `data/model_metadata/hg38_md5.txt`.
`code/download_reference.py` verifies the MD5 after download to `runtime/`.

## Replicates, effects and conservative scope

**Ulirsch:** GEO **GSE87711**, K562 CTRL only, centered 1/2 insert. Do not use
GSE70531, which is a different microarray submission. Pool barcode counts per
variant/allele, retaining all six RNA libraries and two DNA libraries. A shared
DNA denominator makes RNA-only replicate agreement optimistic, so the noise
estimate includes RNA and DNA variance separately. Independence of all six RNA
libraries at the biological-transfection level is not fully documented in the
retrieved processed resource. The approximation and unidentifiable shared
covariance are reported, rather than presenting six RNA libraries as six
independent complete effect measurements. GATA1 perturbations are excluded.

**Abell:** GEO **GSE174534**, GM12878, three reported paired independent RNA/DNA
transfections. Use 1KGTest/1KGControl `Base` canonical SNVs; exclude haplotype
edits and reused Sabeti controls. Coordinates use the pinned author hg19/hg38
conversion. Retain a lexicographically first source identifier per QC-eligible
identical genomic variant, never the most significant or largest-effect
construct. The other acquired constructs remain available for sensitivity work.

Apply a fixed DNA-count minimum of 10 in each reference/alternate DNA sample,
verify hg38 REF, exclude source-overlapping 1 Mb blocks, and exclude Abell blocks
appearing in Ulirsch. All exclusions and replicate components are saved. These
are independent studies and disjoint audit loci, not proven model-pretraining
holdouts or proven independent LD blocks.

ALT/REF log2 RNA ratio minus ALT/REF log2 DNA ratio uses pseudocount 0.5. Library
size normalization cancels inside an allelic ratio. These reprocessed estimates
are not DESeq2 publication coefficients or clinical disease classifications.
Source Kircher coefficients remain as published; no SE is manufactured from
p-values. Source/target effects can differ in scale and biology.

## Models and licensing

Primary inference uses these published, ungated checkpoints. Metadata licenses
are **CC BY 4.0**, with attribution obligations. No AlphaGenome-style prohibition
on output-based fitting was found in the checkpoint model-card licenses. Their
code licenses and copyright notices are included unchanged. This statement is
about these pinned public releases, not a promise about every future model.

| Model | Pinned weights URL | SHA-256 from public Hub LFS metadata |
|---|---|---|
| Enformer PyTorch conversion | https://huggingface.co/EleutherAI/enformer-official-rough/resolve/affe5713ae9017460706a44108289b13c5fee16c/pytorch_model.bin | `99b09d602e195d89c7d4debe144bb2f43907ba0d74006e97098e99d9171c439c` |
| Borzoi replicate 0 | https://huggingface.co/johahi/borzoi-replicate-0/resolve/521db103372769581b8d1d616e6a3b496da490c5/model.safetensors | `d839593c3350bc3aea93b7a6d1fc3ce2b1629e0f4322ecb3c6b6743d14932f75` |

Weights are 1,005,149,571 and 743,692,788 bytes respectively. Model-card and LFS
metadata were downloaded anonymously and bundled. The weight files themselves
were NOT downloaded/run in the authoring sandbox; the hashes are verified
publisher/Hub metadata, and each actual weight download is checked on Wahab.

Enformer code: https://github.com/lucidrains/enformer-pytorch/tree/29cd2529f31cc41456da7511e720d730f0bdd19d (MIT).
Borzoi code: https://github.com/johahi/borzoi-pytorch/tree/8a05eb15dd56871806771b2870c98a0c1dd9e0c8 (Apache 2.0 / retained third-party notices).
The pinned native implementations are vendored **unmodified**, including
Enformer's precomputed corrected TensorFlow gamma tensor. Source manifests give
actual source-file checksums. Vendoring avoids incompatible upstream declared
Transformers pins; shared Transformers 4.57.6 native imports were tested.

Enformer input is **196,608 bp**, matching the corrected conversion's 1,536
attention positions. This is not an exact claim to reproduce the original
TensorFlow 393,216 bp execution. Borzoi input is **524,288 bp**. Full context is
retained; only its final pointwise output embedding is cropped to 512 bins.
A GPU contract verifies identical center predictions before scoring. Both
models use real ALT forward passes, forward/reverse-complement averaging,
float32, no dropout, fixed seed and explicit strand remapping. No linear
surrogate or fabricated outputs are substituted.

Frozen metadata chooses DNASE/ATAC/CAGE/H3K27ac tracks in ordered cell/proxy
contexts. Track files:
- Enformer: https://raw.githubusercontent.com/calico/basenji/06ce5d387e20b47184d05433b3983163c5f923cd/manuscripts/cross2020/targets_human.txt
  SHA-256 `d90233175a9fec389e5ac04cc4aac8cfc08c4aa86217e2e5c10221de78023580`.
- Borzoi: https://raw.githubusercontent.com/johahi/borzoi-pytorch/8a05eb15dd56871806771b2870c98a0c1dd9e0c8/borzoi_pytorch/precomputed/targets.txt
  SHA-256 `cb0c1d4baf514f44de458364b2ae0b85f722a80b72cbe66c56bd700faa481785`.

No exact cell match is promised where the source assay or model lacks it. Mouse
cell/human insert Kircher assays use declared human proxy tracks. Missing
modalities remain missing, with indicators; required DNASE absence stops a run.
`results/track_metadata_preflight.csv` lists every selected context.

## Locked statistical question and decision

Read `config/protocol.json` before scoring; its SHA lock is checked on every
resume/audit. `protocol_initial.json` and a pre-output source-QC amendment are
retained. This is an **internally locked preregistered plan on retrospective
public data**, not an externally registered or blinded clinical study. The
existing author-reported AlphaGenome findings are known and recorded separately.
Do not edit thresholds, features, folds, datasets or modality after new outputs.

Fit on 14 Kircher grouped families, calibrate on 5, freeze, and evaluate both
independent datasets unchanged. The point forecast stays **Enformer DNASE**;
Borzoi DNASE point accuracy is also reported, with all four modality point correlations exported as diagnostics. The primary DNASE proxy is not a direct reporter-expression output; a negative audit does not establish that every model readout fails. Neither model is selected using
external labels. Primary reliability is a source-trained sign-error risk
classifier with fixed hyperparameters and source-only isotonic calibration.
Ranking uses raw risk; monotonic probability calibration is not credited with
a ranking improvement.

The risk features include model-output spread within each fixed candidate block; this requires a candidate set and makes no claim for an isolated variant.

Compare magnitude-only, cross-track dispersion, signed Enformer/Borzoi
disagreement, and deterministic random ranking. Same per-block retained count
for every policy at 10%, 25%, 50%, 75%. A conditional test also matches source-
frozen magnitude-quintile selection budgets to magnitude-only. This controls
magnitude coarsely; residual continuous-magnitude confounding remains possible.

Fixed measured-effect threshold is **0.5 log2**, used ONLY to define audit
labels, not select deployment variants. Screen and evaluate signs using disjoint
replicates, rotate views, and average each variant once. Audited opportunities
are reported explicitly, so a method cannot claim success just by rejecting all
measurably active variants. Noise-ceiling correction is the classical additive
measurement-error approximation, not an established unbiased biological truth.
Report unidentified or unstable corrections; never clip them into [0,1].

**PASS only if BOTH studies** meet all locked conditions: >=20 retained audited
variants and >=10 contributing blocks; positive 97.5% paired block-bootstrap
lower CI for risk-vs-magnitude error reduction; positive 95% lower CI under
conditional magnitude-bin matching; >=90% of the baseline's audited
opportunities; positive full-retention point-correlation 95% lower CI. Bootstrap
10,000 times, seed 812071. Adequate negative evidence is FAIL; incomplete,
missing-model or insufficient evidence is INCONCLUSIVE. Secondary results
cannot rescue a primary FAIL. Missing Borzoi outputs do not silently switch to a
more favorable Enformer-only primary: resume scoring or declare INCONCLUSIVE.

Partial rank correlations audit prediction spread/disagreement vs corrected
block accuracy after controlling measured spread and mean replicate-derived
SE. These **audit** covariates never enter the deployable score. Conditional
associations are exploratory and can still contain residual confounding.

## Wahab / OnDemand operational details

Use ODU OnDemand -> Files to upload the ZIP, then the terminal/shell to unzip it.
Do not score on the login node. `hpc/setup_env.sh` creates an isolated Python
3.11/3.12 environment from your `CARDIAC_PYTHON` interpreter, installs Torch 2.6.0
CUDA 11.8 plus pinned requirements, downloads/checks the genome, regenerates
assay tables and runs software/metadata tests. It never overwrites your old env.
In the submitting shell explicitly point `CARDIAC_PYTHON` to this new `.venv`.

The array uses Wahab's **gpu** partition, one GPU/task, 64 GB host RAM, four CPU
threads. Eight stable hash shards per model give 16 tasks; concurrency defaults
to two GPUs. Increase only the array concurrency if your ODU allocation permits.
Native GPU preflight uses actual weights and hardware before scoring. GPU model
runs are not benchmarked here; no total six-hour guarantee is claimed for this
large all-variant workload. Each submission has a 5h59 wall limit and stops
scoring early to flush resumable checkpoints. A single long shard may require
multiple submissions. No 'larger GPU unavailable' team limitation is asserted;
Wahab is the compute resource, and its actual allocation is reported by preflight.

**Resume:** exit 75 means a valid partial journal, not success. Resubmit the same
array with `bash hpc/submit.sh`. Locks prevent concurrent writes. Complete
shards validate their identities and exit without more inference. If a job is
killed during a row write, only a final partial JSON line is discarded and the
variant is rescored. Reproducibility hashes cover data, protocol, code, weights,
Torch, CUDA/cuDNN and output dtype. Device type is recorded for every scored row; bitwise repeatability is checked on each device type, not promised across GPU architectures. Do not change source files or environments mid-run.
If a previously submitted CPU dependency is permanently unsatisfied after a
partial array, cancel that old CPU job ID using `scancel JOB_ID`; the new submit
creates a fresh dependency. The automatic CPU audit starts only after ALL array
tasks exit successfully and independently verifies every expected row.

If your account requires Slurm account/QOS flags, use your existing authorized
ODU values via `SBATCH_ACCOUNT`/`SBATCH_QOS`; no account name is invented here.

## Optional AlphaGenome, permission gate closed

The primary analysis needs no AlphaGenome key. Optional raw acquisition:
`sbatch --export=ALL hpc/alpha_cpu.slurm` (CPU, resumable, <=5h59/job).
Reads `$HOME/.config/alphagenome/key` privately; never prints the key. Uses
single-variant Atlas queries for scattered validation SNPs, with live ALL_FOLDS
for explicit library REF mismatches. It does not query whole 1 Mb windows
containing millions of unneeded substitutions. Recommended raw DNASE/ATAC/CAGE
masks are 501 bp, histone 2001 bp; these are not exactly aligned to open-model
512 bp score bins. Metadata/allele checks prevent silent cohort mixing.

`code/alpha_optional.py acquire --compare-v1-features /path/old/features.csv`
can record raw Kircher score differences against the existing v1 acquisition;
it reacquires all outputs instead of silently splicing different server runs.

ANY source-fitted AlphaGenome normalization, risk training/calibration or
partial-regression audit requires applicable output-use permission and the
explicit `ALPHAGENOME_CALIBRATION_PERMISSION=1` gate. **Do not set the variable
merely because you have API access.** If that permission exists, the command is:
`"$CARDIAC_PYTHON" code/alpha_optional.py gated-fit`. It fits a separate source-
frozen three-model secondary analysis into `results/alphagenome_optional/`.
No Slurm script automatically opens the gate, and optional outputs never replace
the primary paper or verdict. Without permission, leave this step unrun.

## Outputs and automatic manuscript

After complete open-model scoring, the dependent CPU job writes:
- `results/SUMMARY.txt`, `audit_summary.json`: first-line PASS / FAIL / INCONCLUSIVE.
- `source_point_accuracy.csv`, `point_accuracy_by_dataset.csv`.
- `variance_robust_accuracy_by_block.csv`: replicate noise reliability and status.
- `matched_retention_by_block.csv`, `precision_with_bootstrap_CIs.csv`.
- `primary_bootstrap_comparisons.csv`, `conditional_magnitude_selection.csv`.
- `partial_correlations.csv`, `probability_calibration.csv`, `external_predictions.csv`.
- Four figures in PNG (320 dpi), PDF and SVG: independent gain, risk-retention,
  adjusted spread/disagreement association, and external probability calibration.
- `preprint/Workshop_Preprint.pdf` and editable `Workshop_Preprint.md`, exactly
  four PDF pages, with verified numbers/figures populated automatically.

The included manuscript is a genuine working draft with visibly pending figure
panels, actual QC counts and known negative pilot evidence. It is not submission-
ready evidence for a new reliability method. The PDF builder verifies completed
model states and result-file hashes; it never fills from temporary test mocks.
Calibration probabilities concern the effectful source assay target and need
external checking, not a universal correctness guarantee. Observational MPRA
validation is insufficient for clinical interpretation.

`FILE_SHA256.json` records all distributed file hashes; locally generated
runtime files and later result versions are naturally different from delivery.
No API key, large weights, genome, environment, old proposal or R21 document is
included or modified. This ZIP is self-contained apart from the listed public
model/genome downloads and your authorized Wahab account/API service.
