# Can anything recover the lost signal? Three preregistered attempts

The main audit shows that Enformer and Borzoi DNASE predictions correlate about 0.58 with Kircher
saturation-mutagenesis effects but only 0.03 (Abell) and 0.10 (Ulirsch) on the two independent
studies, far below their noise ceilings of 0.67 and 0.70. Confidence filtering does not recover it.
This directory holds three attempts to recover it. Only the first was preregistered; attempts 2 and 3
were run afterwards and are exploratory. None of them recovers the signal, and attempt 3 produces one
modest positive on a single study. Results are in `results/`.

## 1. Internal representations (preregistered 14 October 2026, before any feature was computed)

`PREREGISTRATION.md` fixes the design; `score_features.py` (GPU) extracts features and `probe.py`
(CPU) fits and evaluates. Feature sets, all ALT minus REF over a 512 bp centre window with forward
and reverse-complement passes averaged: **A** the audit's DNASE score (no probe, reported as the
floor); **B** the full 5,313-track output vector; **C** the 3,072-dimension trunk embedding over the
four centre bins; **D** a 6-mer count difference in a 21 bp window (sequence-only control). Ridge,
alpha by grouped CV inside the training study, trained on Kircher and tested on the other two.

Preregistered success required C minus B above zero and C at least 0.15 on both transfers.
**Result: failed on both.** C reaches 0.010 on Abell and 0.014 on Ulirsch, below B (0.041, 0.106).
The internal representation carries less usable signal than the output head, not more.

## 2. Training on the target study itself (`solve.py`)

Within-study grouped 5-fold CV, and a learning curve giving the model 10, 25, 50 and 100 percent of
the target study's genomic blocks, with and without the Kircher-trained prediction as an extra
feature. The alpha grid is extended to 1e8 because the first run selected the grid boundary.

**Result: failed.** Within-study CV reaches 0.041 (Abell) and 0.038 (Ulirsch). The full learning
curve ends at 0.041 and 0.105. Adding the Kircher prior changes nothing. The limit is not which
study the probe is trained on.

## 3. An MPRA-trained model (exploratory, `malinois_score.py`, `paired_contrast.py`)

Malinois (Gosai et al., Nature 2024) is trained on episomal MPRA activity in K562, HepG2 and SK-N-SH.
All 53,344 audit variants were scored as ALT minus REF on a 200 bp genomic window carrying the
published MPRA flanks, forward and reverse-complement averaged (0 reference mismatches). The
checkpoint checksum is in `results/malinois_checkpoint.md5`.

**Result: a modest gain on Abell, none on Ulirsch.** All three heads are reported; none was selected
after the fact. Abell rises from 0.031 [0.018, 0.044] to 0.069 (K562), 0.071 (HepG2), 0.075
(SK-N-SH); paired block-bootstrap differences against Enformer over the 525 Abell blocks are
+0.038 [0.022, 0.054], +0.039 [0.024, 0.055] and +0.043 [0.027, 0.059]. Ulirsch stays at
0.104 [0.047, 0.152] for K562, difference +0.009 [-0.039, 0.044]; the other two heads are slightly
negative. Sign error stays near 0.49 throughout, and 0.07 remains far below the 0.67 that the Abell
replicates themselves reach.

Interpretation is deliberately limited. Abell is GM12878, which Malinois has no head for, so the
three heads cannot be read as a cell-matching contrast, and their similarity does not separate
training objective, architecture or cell identity as the cause of the gain. Ulirsch is K562, a cell
type Malinois was trained on, and shows no gain, but with only 79 blocks that comparison is weak.

## 4. Filtering on predicted ACTIVITY rather than predicted effect (exploratory)

`activity_filter.py` and `activity_control.py`. Kircher samples densely inside 21 curated active
elements; Abell and Ulirsch are sparse variants spread over 525 and 79 genomic blocks. Where the
reference sequence is not active in a reporter there is no allelic effect to measure, only noise.
Every filter tested so far, including the confidence filtering in the paper, keys on the predicted
EFFECT. This one keys on the predicted ACTIVITY of the reference 200 bp window, scored with Malinois.
It is prospective: it never touches the measured outcome.

**Result: accuracy does rise steeply with predicted activity, but the filter FAILS its decision gate.**
A success gate was fixed by an external reviewer before the output was read: primary forecast Enformer
DNASE, primary retention the top 25 percent by predicted activity, and a required improvement over a
magnitude rule retaining the SAME number of variants in each genomic block, with an interval excluding
zero. `activity_gate.py` scores it.

Retaining the most active variants, Pearson against the measured effect:

| retained | Abell, Enformer | Abell, Malinois | Ulirsch, Enformer | Ulirsch, Malinois |
|---|---|---|---|---|
| all | 0.031 | 0.075 | 0.095 | 0.104 |
| top 50% | 0.043 | 0.097 | 0.136 | 0.133 |
| top 25% | 0.065 | 0.125 | 0.213 | 0.225 |
| top 10% | 0.095 [0.037, 0.153] | 0.170 [0.126, 0.212] | 0.335 [0.121, 0.511] | 0.348 [0.136, 0.507] |

The rise is not range restriction: the SD of the measured effect is flat or lower in the retained sets
(Abell 0.729 to 0.724; Ulirsch 0.514 to 0.450), random subsets of the same size give about 0.00, and
random subsets matched on measured effect magnitude sit at the unfiltered value.

**But the block-matched magnitude comparator does just as well.** Paired differences, activity minus
magnitude, 4,000 block bootstraps: Abell Enformer -0.002 [-0.019, 0.015] at top 25 percent and
-0.014 [-0.051, 0.021] at top 10; Ulirsch Enformer +0.021 [-0.018, 0.078] and +0.026 [-0.042, 0.106].
Every interval includes zero, and on Abell the magnitude rule is marginally ahead. Under the gate this
is descriptive enrichment, not a better filter: **what helps is selecting a few variants per genomic
block at all, not selecting them by predicted activity.**

The per-stratum replicate noise ceiling also rises with selection (Ulirsch 0.699 unfiltered to 0.777
for the activity set and 0.803 for the magnitude set), so a large part of the apparent gain is
recovered signal-to-noise rather than better discrimination.

One finding does survive, and it is a nuance the paper should carry: the confidence filtering in the
paper applies global thresholds and does not beat the unfiltered baseline, whereas selecting within
genomic blocks lifts Ulirsch from 0.095 to above 0.30 by either criterion. The unit of selection, not
the selection variable, is what matters here. That needs its own preregistered test before it is
claimed.

## 5. Why does it fail? Measurement noise versus model failure (`reliability_decomposition.py`)

The previous four attempts establish that the models do not transfer; none of them explains why. The
most common explanation in the field is that the external assays are too noisy to predict. Each
audited external row carries a replicate-derived `noise_variance`, so variants can be stratified by
reliability and the correlation divided by the per-stratum noise ceiling, giving model accuracy after
removing what measurement noise alone can account for. The decision rule was written into the script
before the output was read: a flat ratio means noise does not explain the collapse, a steeply rising
ratio approaching the Kircher value means it does.

**Result on Abell: measurement noise does NOT explain the failure.** In the most reliably measured
fifth of Abell variants the noise ceiling is 0.971, so the measurements are nearly noise-free, and
Enformer's correlation there is 0.015 [-0.014, 0.042]. The noise-adjusted ratio is flat and near zero
across all five strata (0.052, 0.058, 0.045, 0.015 for quintiles 2 to 5). Malinois behaves the same
way (ratio 0.102, 0.106, 0.102, 0.083). On Abell the models are simply wrong, not merely measured
badly. This is the clearest single result in this directory and it closes off the standard excuse.

**Result on Ulirsch: suggestive, not conclusive.** The most reliable fifth reaches 0.285 [0.064,
0.485] for Enformer and 0.308 [0.169, 0.434] for Malinois against a ceiling of 0.912, so ratios of
0.31 and 0.34 against a pooled 0.14 and 0.15. But the pattern is not monotonic (quintile 3 is
negative), each stratum holds only 475 variants in about 65 blocks, and the intervals are wide.

**Caveat that weakens the comparison.** The preregistered positive control could not be run: the
released Kircher table carries no `noise_variance`, so Kircher's noise-adjusted ratio cannot be
computed from this export and the external ratios have no internal reference point. The Kircher
figure quoted as a target in the script header (0.71) is derived from the published pooled ceiling,
not from the same per-variant procedure, so it is indicative only.

## 6. Were the inputs wrong? Two sensitivity checks (`malinois_score.py --insert-len`, `ulirsch_real_oligos.py`)

An external reviewer raised a fair objection: Malinois was applied to reconstructed genomic windows,
while these studies assayed short synthetic inserts. A reconstructed window contains sequence the
assay never tested, so the comparison could have been unfair to the model.

**Insert length makes no difference.** Rescoring every variant with a 150 bp variant-centred window
instead of 200 bp gives Abell SK-N-SH 0.074 against 0.075, Abell HepG2 0.070 against 0.071, and
Ulirsch K562 0.103 against 0.104.

**The real assayed sequences make no difference either.** The Ulirsch study design file in the audit
package carries the actual 145 bp oligos, and all 2,375 primary Ulirsch rows join to it on
`legacy_row_id`. Scoring those exact sequences, averaging over the tile positions the study used,
gives K562 0.110 [0.061, 0.154] against 0.104 for the reconstructed window and 0.095 for Enformer
DNASE, with sign error 0.488. Abell has no sequences in this package, only coordinates.

So the weak external performance is not an artifact of how the input sequence was constructed. This
is a negative control rather than a result, and it removes an objection rather than adding a finding.

## Reading

Within the models and probes tested here, none of the first three attempts recovers cross-study
performance, and one MPRA-trained model gives a modest improvement on one of the two external
studies. The fourth attempt changes which variants the model is asked about rather than the model itself.
Accuracy rises steeply with predicted activity, but it fails its decision gate because an equally
simple magnitude rule applied within the same genomic blocks does the same thing. These experiments bound what these particular models and readouts achieve. They do not
establish that the targets are unlearnable from sequence: a replicate-derived ceiling says the
measurements are reproducible, not that sequence alone can reach them.

## Known limitations in this directory

- `probe.py` standardises features on the whole training study before the inner alpha folds. The
  leak is confined to the training study and affects only the alpha choice, not the held-out test,
  and C trails B by far more than the effect could explain, but the fold-local version is correct.
- In `solve.py` the learning-curve rows report `rep_min` and `rep_max` across repeats, which are the
  spread over repeats and not confidence intervals. The 100 percent budget uses a single repeat.
- With the alpha grid extended to 1e8, Kircher-to-Ulirsch C rises from 0.014 to 0.091. This does not
  change the preregistered outcome, which compared C against B under the original grid.
