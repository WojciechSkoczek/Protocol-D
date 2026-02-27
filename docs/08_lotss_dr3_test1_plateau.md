# LoTSS DR3 — Test 1: Number-count dipole plateau (fit_nside-invariant)

## Purpose
This test is a simple, tracer-specific sanity anchor for Protocol D. It measures a *pure number-count dipole* from the LoTSS DR3 source catalogue and checks **stability (plateau)** under:
- flux threshold scans (bright cuts),
- Galactic latitude cuts (|b| >= 10, 20, 30 deg),
- an explicit survey footprint mask.

If a recovered dipole is physical, it should not change violently under these reasonable perturbations. If it does, the instability is diagnostic of survey / selection effects.

## Input
- LoTSS DR3 v1.0 source catalogue (FITS): `LoTSS_DR3_v1.0.srl.fits`

Store the raw FITS locally (do not commit):
- `data/raw/LoTSS_DR3_v1.0.srl.fits`

## Method (summary)
1. Convert (RA, DEC) -> Galactic (l, b).
2. Build a footprint mask from sources above a *footprint flux threshold* and the same |b| cut.
3. For each flux cut, count sources per HEALPix pixel and fit:
   `N = a0 + ax*x + ay*y + az*z`
   Dipole vector `d = (ax, ay, az) / a0`.
4. Report amplitude, direction, and step-to-step deltas vs previous flux cut.

### Important implementation detail (why `fit_nside`)
A naive dipole fit on a very high-resolution pixelisation can become **resolution-dependent** when maps are sparse (bright cuts -> many empty pixels). To prevent estimator artefacts, we fit on a fixed `fit_nside` (default 64). This makes the plateau test resolution-invariant.

### Footprint clamping
If `footprint_flux_mjy` is set higher than the lowest plateau flux cut, the code clamps it down to the lowest plateau cut. This avoids an unstable footprint built from ultra-rare sources.

## How to run (Windows)
PowerShell:
```powershell
.\scripts\run_lotss_dr3_plateau.ps1
```

## Canonical artefact (to commit)
- `results/lotss_dr3/plateau_dr3_counts_fit64.csv`

## Interpreting the plateau
Empirically, LoTSS DR3 looks cleanest for |b| >= 20 deg at bright cuts (e.g. 60–80 mJy). |b| >= 10 deg tends to be more sensitive to residual systematics near the Galactic plane / footprint edges.
