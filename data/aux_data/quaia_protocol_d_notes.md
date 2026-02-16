# Quaia (G<20.5) — dipole sanity check via |b| masks (python-only)

## Dipole vs |b| cut (mask sensitivity)
Computed with debias using provided random_10x catalogue (same |b| cut applied to random).  
Key point: direction is **sensitive** to the Galactic latitude cut.

For **G<20.5**:
- |b|≥20°: D_obs=0.027559  (l,b)=(343.64,46.26)
- |b|≥25°: D_obs=0.024128  (l,b)=(341.28,54.25)
- |b|≥30°: D_obs=0.022130  (l,b)=(338.27,63.74)

Direction change (20→30): ~17.74°  
Max pairwise (20,25,30): **17.74°**  → used as empirical sigma_dir

Amplitude spread std(20,25,30): **0.002746**  
Shot-noise proxy (sqrt(3)/sqrt(N)): **0.001808**  
Combined D_obs_err (quadrature): **0.003287**

## Protocol D row (adopted)
Adopted baseline: **G<20.5, |b|≥30°** (for consistency with CatWISE |b|≥30°).

Row fields:
- D_obs=0.022130 ± 0.003287
- (l,b)=(338.27,63.74)
- sigma_dir_deg_approx=17.74
- z_mean=1.556

## BAO_P1 harness runs (4 tracers)
### Wagenveld pack (NVSS+RACS+CatWISE+Quaia)
baseline: ΔAICc = -54.35
g_hat ≈ (247.39,46.49)
coherence angles: NVSS 16.4°, RACS 40.7°, CatWISE 38.4°, Quaia 57.4°

### Boehme pack (NVSS+LoTSS+CatWISE+Quaia)
baseline: ΔAICc = -46.90
g_hat ≈ (240.39,37.17)
coherence angles: NVSS 10.2°, LoTSS 20.5°, CatWISE 27.7°, Quaia 68.0°

## Practical take
- Quaia adds leverage in z (A_eff falls with z in BAO_P1), but its **direction is still a clear outlier** (angles ~56–68°).
- Given the strong |b| sensitivity, it is safer to treat Quaia as **amplitude‑only / very low‑weight direction** until we model Gaia systematics more explicitly (e.g., scanning‑law residuals).

Files produced:
- quaia_dipole_maskscan.csv
- protocol_d_row_quaia_G20p5_bscan.csv
- updated input packs with Quaia appended (2 files)
- harness outputs (4 txt)
