# Data provenance and layout

This repository is structured so that a reader can (a) reproduce the run bookkeeping and (b) inspect the archived outputs.

## Repository data layout (high level)

- `configs/`
  - Default configuration files (analysis parameters, thresholds, paths).
- `data/`
  - Auxiliary inputs used by the pipeline (masks, reference tables, small helper assets).
  - Run input packs that define reproducible runs (cuts, catalog selection, seeds, etc.).
- `results/`
  - Archived outputs produced by the runs (tables/CSVs/PNGs) and run indices.
- `scripts/`
  - Convenience scripts for reproducing and validating runs.
- `src/`
  - The Python implementation.

## External datasets

Large, external catalogues are not re-hosted here unless explicitly stated in `data/`.
If a run depends on an external catalogue, the run definition should:
- reference the catalogue by name,
- document where to obtain it,
- include the exact version / release tag when possible.

## Reproducibility notes

- All run definitions should be captured in `data/run_input_packs/` (or equivalent) so a run can be re-launched.
- Outputs should be accompanied by a minimal metadata record: timestamp, git commit hash, config hash, and random seed.
