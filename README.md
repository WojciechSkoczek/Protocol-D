# Protocol D - BAO-anchored dipole harness (v0.6.2)

Lightweight reproduction repository for the analysis harness used in Protocol D.

## What it does
- compares NULL vs BAO_P1 models for a set of tracer dipoles
- fits a global excess direction g_hat and BAO-anchored effective amplitudes A_eff
- uses an approximate directional covariance (sigma_mag + sigma_dir proxy)

## What it does not do
- does not re-measure dipoles from raw catalogues (estimator + mask handling live upstream)

## Quick start
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
bash scripts/repro_bao_only.sh
```
Outputs are written to REPRO_OUT/.

## What to verify (publishable checks)
- Key claims + where to verify: docs/07_key_claims.md
- One-stop summary tables: docs/06_master_summary.md

## Release checklist
See RELEASE.md (tagging, Zenodo DOI, verification).

## Windows note
On Windows, aux is a reserved device name (Win32/NTFS). This repo therefore uses data/aux_data/.

## AI assistance
See `docs/ai_assistance.md`.

Context for researchers: docs/context_summary.txt

