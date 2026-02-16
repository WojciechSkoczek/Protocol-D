# Protocol D — Run N (two-axis parametric operator)

## Model
Excess vectors: E_i = D_obs - D_kin

E_i ≈ β w(z_i) ĝ  +  t_i * [A1(x_i) ŝ1 + A2(x_i) ŝ2]

t_i=1 for Quaia points only; for non-Quaia, operator term is forced to zero.

A1(x)=a·x, A2(x)=b·x, where x includes (bmin-30, mag-20, pm_sig, pm_tot, pm_none, z-1, is_pmcut).

ŝ2 is orthogonalised to ŝ1 by Gram–Schmidt during optimisation (keeps the problem identifiable).

## Outputs
- RunN_summary.csv: directions, chi2, AICc, coefficients (a and b), and prediction diagnostics for both axes.
- RunN_diagnostics.csv: per observation, measured projections on ŝ1, ŝ2 vs predicted amplitudes.
- plots/: sky_*, axis1_pred_vs_meas_*, axis2_pred_vs_meas_*, maskscan trends for both axes.

## Quick reading
If Run N improves AICc over Run M (dAICc_N_minus_M < 0), the second operator axis is justified.
If not, Run M already captures the operator manifold sufficiently with one axis.
