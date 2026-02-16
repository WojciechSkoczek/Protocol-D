# Protocol D — publishable repository (v0.2)

This repository packages the **Protocol D** work completed in this project:

- the *BAO-anchored dipole harness* (reproduction from derived dipole packs)
- the *operator / Banach framing* runs (K–R): treating the survey pipeline as an operator that can induce coherent “dipole-like” structure
- canonical **run outputs** (CSV summaries, diagnostics, plots) for Runs G–R

## Scope (what this repo is and is not)

🟢 **In scope:** reproducible analysis given *derived* dipole measurements (per-tracer rows), including model fits, AICc comparisons, leave-one-out checks, and operator identifiability tests.

🟡 **Not in scope (by design):** re-measuring dipoles from the full raw catalogues and masks. That is estimator- and licensing-dependent. This repo is meant to be transparent about the *analysis layer* and its assumptions, while keeping the raw-catalogue layer out of band.

## What to read first
1. `README.md` (quick start)
2. `docs/02_models_and_assumptions.md`
3. `docs/04_runs.md` (what each Run means and where its outputs are)



See also: `docs/06_master_summary.md` for one-stop summary tables.


See also: `docs/07_key_claims.md` for a publishable checklist of claims and verification paths.
