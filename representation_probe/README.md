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

## Reading

Within the models and probes tested here, none of the three attempts recovers cross-study
performance, and one MPRA-trained model gives a modest improvement on one of the two external
studies. These experiments bound what these particular models and readouts achieve. They do not
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
