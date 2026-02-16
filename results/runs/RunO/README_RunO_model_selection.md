# Protocol D — Run O (Run M model selection / operator feature ladder)

## Goal
Run M showed that a *single* Quaia-only operator axis explains most selection-driven dipole drift.
Run O asks: **which operator features are actually needed?**

We fit the same model as Run M but with different feature subsets for the operator amplitude A(x):
E_i = D_obs - D_kin ≈ β w(z) ĝ + t_i * A(x_i) ŝ  (t_i=1 for Quaia only)

## Ladder
- O0_intercept: A = a0
- O1_bmin: A = a0 + a_b*(bmin-30)
- O2_bmin_mag: + a_mag*(mag-20)
- O3_bmin_mag_pm: + proper-motion terms (pm_sig, pm_tot, pm_none, is_pmcut)
- O4_full_plus_z: + redshift term (z-1)

## Reading shortcut
For each (pack,variant) look at aicc_ladder_*.png.
The best model is the lowest AICc.

## Main result
Across all datasets here, **O2_bmin_mag** is the AICc winner.
Adding PM terms and z improves chi2 slightly but is penalised by extra parameters.

See RunO_best_model_per_dataset.csv and RunO_aggregate_by_model.csv.
