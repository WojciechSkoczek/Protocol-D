# Protocol D — Run K (operator / Banach framing)

## Idea
Treat dipole measurements as outputs of a *survey operator* (mask, selection function, magnitude cut). Instead of forcing one 3D vector to fit everything (BAO-only), we model an additional *operator bias axis* ŝ whose amplitude depends on mask parameters.

Model in excess-vector space E_i = D_obs_i - D_kin:
- Core (cosmology-like):  E_i ≈ β w(z_i) ĝ
- Operator term (Quaia-only):  + t_i * (a0 + a1*(bmin-30) + a2*(magcut-20)) ŝ
where t_i=1 for Quaia points, else 0.

We fit (ĝ, β, ŝ, a0, a1, a2) using the same 3D covariance weighting as earlier runs, with a Huber penalty on whitened residuals.

## Outputs
- RunK_operator_summary.csv: core ĝ and sys ŝ per pack/variant + AICc comparison.
- RunK_operator_diagnostics.csv: per-observation residual projection on ŝ vs predicted amplitude.
- plots/: sky_{pack}__{variant}.png, maskscan_pred_vs_meas_*.png, bmin_trend_*.png.

## Key readouts (what to look at)
1) Core direction ĝ (Run K) should match the *non-Quaia* direction, while BAO-only ĝ may drift when Quaia is added.
2) The maskscan points should satisfy proj(residual on ŝ) ≈ predicted amplitude (very high correlation).
3) Coefficients a1<0 means increasing |b| cut reduces the operator bias; a2>0 means deeper magnitude cut increases it.
