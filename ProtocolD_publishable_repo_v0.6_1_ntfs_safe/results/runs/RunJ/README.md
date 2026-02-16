# Protocol D — Run J (Quaia systematics stress test)

Goal:
- Quantify how strongly BAO_P1 fitted direction ĝ depends on how we treat Quaia directional uncertainty and mask choices.

What varies:
- We start from the non-Quaia packs (NVSS + CatWISE + LoTSS/RACS).
- We append ONE Quaia row in several variants (or use the 3-bin Quaia packs), then refit BAO_P1.
- We repeat with 'amplitude-only' downweighting: sigma_dir -> 90 deg for Quaia (scenario suffix _sigma90).

Scenarios (key):
- no_quaia: baseline reference (3 tracers).
- quaia_std_G20p5_b20_selcorr: small errors (direction 5 deg, mag err ~0.0016) — can dominate the fit.
- quaia_b30_debias: Quaia debiased with random_10x, |b|>=30, sigma_dir from mask sensitivity.
- quaia_conservative_sys: conservative errors (sigma_dir 20 deg, mag err 0.003).
- quaia_maskscan_bmin20/25/30: uses l,b,D from maskscan with a shared conservative sigma_dir + D_err.
- zbins_default / zbins_sigma90: 3-bin Quaia (G<20.0, pmSig<2) as provided, with/without sigma90 downweighting.

Files:
- RunJ_summary_long.csv: all scenarios for each pack/variant.
- RunJ_summary_key.csv: the subset used in plots.
- RunJ_pivot_ang_to_CMB.csv: quick view (degrees).
- RunJ_pivot_dAICc.csv: quick view (ΔAICc).
- RunJ_contrasts_zbins_sigma90.csv: the key comparison for the z-bin case.
- plots/: per-pack plots for angle-to-CMB and ΔAICc.

Interpretation guide:
- If 'CMB alignment' is robust cosmology, ĝ should not change much when Quaia direction is downweighted.
- If alignment is driven by mixing/selection effects, inflating sigma_dir for Quaia should pull ĝ back toward the non-Quaia direction.
