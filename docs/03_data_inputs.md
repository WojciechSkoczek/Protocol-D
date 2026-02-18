# Data inputs and provenance

## Derived input packs (public-friendly)
The analysis uses **derived packs** under `data/packs/`. Each row represents one dipole measurement and includes:

- D_obs, Ïƒ_mag, Ïƒ_dir, (l,b), z_mean
- D_kin computed from (x, Î±, Î²_assumed) (recorded as `D_kin_from_x_alpha`)
- bibliographic reference and notes

These are sufficient to reproduce all harness-level fits.

## Additional derived inputs (Quaia systematics)
Additional helper tables and notes live in:
- `data/aux_data/` (raw exports and scans)
- `results/aux/` (tables and plots copied from the project workspace)

## What is *not* redistributed here
Raw catalogues, masks, and selection-function FITS that may be subject to third-party licences are intentionally not bundled.  
If you want full end-to-end reproduction, create a private `data/raw/` folder and follow the upstream data licences.



## Windows-safe pack filenames

The derived packs are stored under short IDs in `data/packs/`.
Use `data/pack_ids.csv` to map IDs to the original full filenames.
The original full-named pack files are archived as `data/packs_full.zip`.


## Windows NTFS note (reserved names)

On Windows, `aux` is a reserved device name, which can break Git and some tools even if the folder is visible. For this reason the repo uses `data/aux_data/` (formerly `data/aux/`).

