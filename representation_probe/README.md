# Can anything recover the lost signal? Three preregistered attempts

The main audit shows that Enformer and Borzoi DNASE predictions correlate about 0.58 with Kircher
saturation-mutagenesis effects but only 0.03 (Abell) and 0.10 (Ulirsch) on the two independent
studies, far below their noise ceilings of 0.67 and 0.70. Confidence filtering does not recover it.
This directory holds three attempts to recover it. All three fail. Results are in `results/`.

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

## 3. An MPRA-trained model in the matched cell type (`malinois_score.py`)

Malinois (Gosai et al., Nature 2024) is trained directly on lentiMPRA in K562, HepG2 and SK-N-SH.
All 53,344 audit variants were scored as ALT minus REF on a 200 bp genomic window with the MPRA
flanks, forward and reverse-complement averaged (0 reference mismatches).

**Result: a small real gain on Abell, none on Ulirsch.** Abell rises from 0.031 to 0.069 (K562),
0.071 (HepG2), 0.075 (SK-N-SH); the paired block-bootstrap difference against Enformer is
+0.038 [0.022, 0.054] for K562 and similar for the others. Ulirsch stays at 0.104 [0.047, 0.152],
difference +0.009 [-0.039, 0.044]. Two observations matter more than the gain itself: all three
Malinois cell types perform equally on Abell, which uses none of them, so the gain comes from
MPRA-style training rather than cell matching; and Ulirsch is K562, the matched cell type, where
there is no gain at all. Sign error stays near 0.49 throughout, and 0.07 remains far below the 0.67
the Abell measurements themselves permit.

## Reading

Three controlled attempts, each with its decision rule written down first, do not recover the
signal. The constraint is the training data rather than the architecture or the readout: predicting
these effects appears to need reporter measurements in the relevant cell type at sufficient scale,
not a different way of reading an existing model.
