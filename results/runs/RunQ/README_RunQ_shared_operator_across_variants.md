# Protocol D — Run Q (operator invariance across variants)

## Question
Does the Quaia selection-driven drift use the same operator geometry (axis ŝ) when we switch variant
(baseline vs 2MRS_clean)?

## Model (O2)
E = D_obs - D_kin ≈ β w(z) ĝ + t * A(x) ŝ
A(x)=a0 + a_b*(bmin-30) + a_m*(mag-20)
t=1 for Quaia only.

## Tests
1) Per pack (boehme, wagenveld): joint fit across variants with shared operator (ŝ and a),
   but separate core (ĝ,β) per variant.
   Key statistic: AICc_joint_minus_sep < 0 supports invariance.

2) Global: one shared operator across all four datasets (two packs × two variants),
   with separate core per dataset.

## Outputs
- RunQ_joint_shared_operator_per_pack_across_variants.csv
- RunQ_global_joint_shared_operator_all4.csv
- plots/: AICc comparisons, axis angles, skyplot of ŝ.
