# Protocol D — Run Lq (operator manifold; Quaia-only low-rank)

## Why this run
Run K showed that Quaia directions shift strongly with mask and magnitude selection.
Run Lq expands the *operator manifold*: we add many Quaia-only 'measurements' produced by changing selection parameters
(|b| cut, mag cut, proper-motion cuts), then fit a model that treats these as *different outputs of the same survey operator*.

## Model
Work with excess vectors E_i = D_obs_i - D_kin.

Core term (shared):  E_i ≈ β w(z_i) ĝ
Operator term (Quaia-only): + c1_i ŝ1   (low-rank k=1), with c1_i fitted per Quaia observation.

Key constraint: **the operator term is applied only to Quaia points** (non-Quaia coefficients fixed to 0).
This prevents the operator subspace from 'explaining away' the radio/IR tracers.

## What to look at
- RunLq_summary.csv: compares directions to CMB and to non-Quaia; shows χ² collapse when operator term is allowed.
- plots/sky_*.png: CMB vs (non-Quaia ĝ) vs (BAO-all ĝ) vs (RunLq core ĝ) vs (RunLq ŝ1).
- plots/c1_vs_bmin_*.png: how the fitted operator coefficient changes with |b| cut; pmcut points sit at bmin=30.
- RunLq_coeff_regression.csv + plots/c1_pred_vs_meas_*.png: linear feature model predicting c1 from (bmin, mag, pm cuts, z).
- quaia_zbins_angle_to_CMB.png + RunLq_quaia_zbins_table.csv: 'plateau-style' check across the three z-bins for G20.0 and G20.5.

## Interpretation shortcut
If this framing is right, the cosmology-like core ĝ should stay near the non-Quaia direction,
while most selection-driven drift gets absorbed into the Quaia-only coefficient c1 along ŝ1.
The small 'CMB distance to plane {ĝ,ŝ1}' indicates CMB often lies close to the 2D subspace spanned by the core and the operator axis.
