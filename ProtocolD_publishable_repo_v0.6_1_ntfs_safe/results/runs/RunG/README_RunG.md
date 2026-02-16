# Protocol D — Run G (operator-template extension)

This bundle contains quick diagnostics for extending BAO_P1 with a second coherent mode tied to a **template** (operator proxy).

Model family tested (Galactic coordinates):
- NULL: D_obs = D_kin
- BAO_P1: D_obs = D_kin + beta * w(z)^p * g_hat
- BAO_P1+SYS: D_obs = D_kin + beta * w(z)^p * g_hat + gamma * t(template) * s_hat

Here:
- g_hat is the fitted "excess" direction
- s_hat is an additional fitted axis (systematic / operator mode)
- template choices used in this run:
  - `radio` = 1 for radio tracers, else 0
  - `quaia` = 1 for Quaia z-bins, else 0
- AICc comparisons are always vs NULL on the same pack+variant.

Files:
- `RunG_fit_summary.csv` — best-fit parameters and LOO stability metrics
- `RunG_loo_details.csv` — leave-one-out refits (direction sensitivity)
- `RunG_direction_plots/*.png` — Aitoff views of fitted axes (CMB vs g_hat and (if present) s_hat)

Notes:
- The 2MRS_clean vs baseline distinction is meaningful for Wagenveld packs; for Boehme CatWISE-only it is identical (same input rows).
- For N=3 (CatWISE-only packs), BAO_P1+SYS (k=6) is under-supported and AICc penalises it heavily. Treat those as a sanity check.
