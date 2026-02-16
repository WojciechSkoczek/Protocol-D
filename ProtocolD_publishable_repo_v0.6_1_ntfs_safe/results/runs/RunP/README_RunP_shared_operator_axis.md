# Protocol D — Run P (shared operator axis test)

## Purpose
Run O showed that the best operator amplitude model is O2_bmin_mag:
A(x) = a0 + a_b*(bmin-30) + a_m*(mag-20)

Run P tests whether the *operator axis* ŝ and coefficients are **invariant** across measurement packs
(boehme vs wagenveld), which would support the 'map operator' interpretation.

## Part 1 — single-dataset O2 fits
File: RunP_O2_single_fits.csv
- Each dataset (pack,variant) fitted independently.
- Compare ŝ directions and (a_b, a_m).

## Part 2 — joint fit with shared operator (ŝ and a) across packs
File: RunP_O2_joint_shared_operator.csv
Model for combined data:
E_i ≈ β_p w(z_i) ĝ_p  +  t_i * A(x_i) ŝ
where p indicates pack (separate core ĝ_p,β_p), but ŝ and A(x) are shared.

### Key readout
- AICc_joint_minus_sep < 0 means the shared-operator hypothesis is supported.

## Plots
- sky_operator_axes.png: ŝ directions (single fits + shared fits) vs CMB.
- angle_s_between_packs.png: angle between ŝ from independent fits.
- aicc_joint_vs_sep.png: joint vs separate AICc.
- maskscan_shared_trend_*.png: predicted vs measured maskscan trend under shared operator.
