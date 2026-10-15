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

**Result: accuracy rises steeply and monotonically with predicted activity, on both studies.**
Retaining the most active variants, Pearson against the measured effect:

| retained | Abell, Enformer | Abell, Malinois | Ulirsch, Enformer | Ulirsch, Malinois |
|---|---|---|---|---|
| all | 0.031 | 0.075 | 0.095 | 0.104 |
| top 50% | 0.043 | 0.097 | 0.136 | 0.133 |
| top 25% | 0.065 | 0.125 | 0.213 | 0.225 |
| top 10% | 0.095 [0.037, 0.153] | 0.170 [0.126, 0.212] | 0.335 [0.121, 0.511] | 0.348 [0.136, 0.507] |

Sign error also moves off chance for the first time on an external study: Ulirsch Enformer falls from
0.476 unfiltered to 0.360 in the top 10 percent. On Abell it stays near 0.47.

The obvious artifact is range restriction, since selecting on activity could simply concentrate large
effects. It does not. The SD of the measured effect is flat or lower in the retained sets (Abell
0.729 to 0.724; Ulirsch 0.514 to 0.450), and two controls at each retention level separate the
selection from chance: random subsets of the same size give about 0.00, and random subsets matched on
the distribution of the measured effect magnitude sit at the unfiltered value (Abell Enformer 0.031
[-0.009, 0.067] against 0.095 selected; Ulirsch Enformer 0.100 [-0.023, 0.265] against 0.335).

This is exploratory and was not preregistered. Thresholds and controls were fixed before the output
was read, and the retained counts are small at the strictest level (2,204 Abell, 238 Ulirsch rows).

## Reading

Within the models and probes tested here, none of the first three attempts recovers cross-study
performance, and one MPRA-trained model gives a modest improvement on one of the two external
studies. The fourth attempt is different in kind: it does not change the model, it changes which
variants the model is asked about, and it works on both studies. These experiments bound what these particular models and readouts achieve. They do not
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
