# Protocol D — Run I (lite): Quaia inclusion ladder

This run uses the existing 3 Quaia z-bins as a ladder to test how the fitted excess direction moves as we include/exclude Quaia bins.

Cases keep the non-Quaia tracers fixed (CatWISE, LoTSS, NVSS) and add selected Quaia bins (low/mid/high).
It is a plateau proxy given we do not yet have >3 Quaia bins.

## Outputs
- RunI_quaia_inclusion_ladder_summary.csv
- plots/angle_to_cmb_*.png, dAICc_*.png, loo_mean_*.png

## How to read it
- If adding Quaia bins pulls ĝ sharply toward CMB while non-Quaia-only stays ~20–40° away, that supports 'mixing' rather than a single cosmological axis.
- Compare LOO stability across cases.

## Run I2 (added): 2-axis SYS plane diagnostics
Added RunI2_sys_plane_summary.csv and plots/cmb_plane_*.png.
This tests whether the CMB direction lies close to the plane spanned by {ĝ, ŝ} for each inclusion case.
