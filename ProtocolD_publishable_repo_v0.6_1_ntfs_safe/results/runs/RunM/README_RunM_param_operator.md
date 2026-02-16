# Protocol D — Run M (parametric operator model)

## What Run M adds vs Run K / Run Lq
Run Lq allowed a free Quaia-only coefficient c1 for each selection variant.
Run M replaces that freedom with a *parametric* operator amplitude:

E_i = D_obs - D_kin ≈ β w(z_i) ĝ  +  t_i * A(x_i) ŝ

A(x_i) = a0 + a_b*(bmin-30) + a_mag*(mag-20) + a_pmSig*(pm_sig_max)
         + a_pmTot*(pm_tot_max) + a_pmNone*(pm_is_none) + a_z*(z-1) + a_isPmcut*(is_pmcut)

Operator term applies only to Quaia points (t_i=1), keeping non-Quaia tracers 'protected'.
We fit (ĝ, β, ŝ, {a_j}) with full 3D covariance weighting; reported AICc uses pure chi2.

## Files
- RunM_summary.csv: per pack/variant, directions, ΔAICc, fitted a_j, and Quaia prediction diagnostics.
- RunM_diagnostics.csv: per observation, measured projection on ŝ vs predicted amplitude.
- plots/: sky_*, pred_vs_meas_*, maskscan_trend_*.

## Reading shortcut
1) Core ĝ should sit close to non-Quaia direction (radio/IR), while BAO-all may drift toward CMB when Quaia is added.
2) High corr(pred,meas) indicates 'map operator' explanation (selection → dipole drift) is internally consistent.
