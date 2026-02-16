# Protocol D — BAO-anchored dipole harness (v0.1)

This is a lightweight reproduction repository for the *harness* used in Protocol D.

What it does:
- compares NULL vs BAO_P1 models for a set of tracer dipoles
- fits a global excess direction `g_hat` and BAO-anchored effective amplitudes `A_eff`
- uses an approximate directional covariance (sigma_mag + sigma_dir proxy)

What it does *not* do:
- it does not re-measure dipoles from raw catalogues (that requires a consistent estimator and detailed mask handling)

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

bash scripts/repro_bao_only.sh
```

Outputs are written to `REPRO_OUT/<timestamp>/`.

## Inputs

Derived input packs are stored under `data/packs/`. They include:
- baseline radio packs (z_radio=1.0)
- packs extended with CatWISE + Quaia (zbins)

## Results

Selected run outputs and summary tables are stored under `results/`.

See `docs/data_provenance.md` for recommended provenance notes to include before making the repo public.


## Operator runs (Banach framing)

```bash
bash scripts/repro_operator_runs.sh
```

See `docs/02_models_and_assumptions.md` and `docs/04_runs.md`.

## What to verify (publishable checks)

- Key claims + where to verify: `docs/07_key_claims.md`
- One-stop summary tables: `docs/06_master_summary.md`

## Release checklist

See `RELEASE.md` for a publishable release checklist (tagging, Zenodo DOI, and verification).

## Windows compatibility note

To avoid Windows path-length issues, `results/runs/` uses short run folder names (RunG…RunR) and short pack IDs.
See `results/runs/README.md` and `data/pack_ids.csv`.
Original full-named pack files are archived as `data/packs_full.zip`.



### Windows note

On Windows, `aux` is a reserved device name (Win32/NTFS), which can break Git and some tools. This repo therefore uses `data/aux_data/` for auxiliary tables and plots.
