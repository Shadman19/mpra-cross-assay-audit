# Representation probe: preregistered test (frozen 14 October 2026, before any feature is computed)

Question: do Enformer internal representations predict independent MPRA allelic effects better than
Enformer's full output-track vector, when a probe is trained on one study and tested on another?

Variants: the primary audit rows (Kircher 28,932; Abell 22,037; Ulirsch 2,375), unchanged.
Target: the signed measured effect (`effect`), as in the audit.

Feature sets (all ALT minus REF, 512 bp centre window, forward and reverse-complement averaged):
- A. Enformer DNASE score already used in the audit (one number, no probe; reported as the floor).
- B. Full Enformer human output vector: log2(ALT+1) - log2(REF+1) over all 5,313 tracks.
- C. Enformer trunk embedding (3,072 dimensions) averaged over the four centre bins.
- D. Sequence-only control: 6-mer count difference in the 21 bp window centred on the variant (4,096 dims).

Probe: standardise features on the training study only; ridge regression with alpha chosen from
10^-1 to 10^5 (13 values) by grouped cross-validation inside the training study (groups: Kircher
`fit_group`; Abell or Ulirsch `genomic_block`). Same recipe for B, C and D. No other tuning.

Transfers: primary Kircher -> Abell and Kircher -> Ulirsch. Secondary Abell -> Ulirsch.
Metric: Pearson correlation between prediction and measured effect on the test study, with 2,000
cluster-bootstrap draws over test genomic blocks; paired contrasts resample the same blocks.
Sign error on active variants uses `audit_engine.errors()`, as in the audit.

Success (all required): on both primary transfers, C minus B has a 95% interval lower bound above
zero, and C reaches Pearson >= 0.15. Anything else is reported as a negative result.
Overlapping loci were already excluded between studies by the audit.
