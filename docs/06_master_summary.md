# Master summary tables (Protocol D)

This folder contains **one-stop tables** that summarise the key outcomes of Protocol D.

## 1) `master_fits.csv`
A consolidated view of the *main fit-level outcomes*:
- Run G: baseline BAO harness (and BAO+SYS where relevant)
- Run O: best model family per dataset (model selection)
- Run P: joint/shared-operator improvement within each variant
- Run Q: single O2 fits per dataset + global shared-operator fit across all four datasets

Recommended columns to inspect:
- `aicc`, `dAICc`
- `g_to_CMB_deg` (core direction vs CMB)
- `s_to_CMB_deg` (operator axis vs CMB, where applicable)

## 2) Subspace + systematics:
- `master_subspace_RunH.csv`  (Run H)
- `master_quaia_ladder_RunI.csv` (Run I)
- `master_systematics_RunJ.csv` (Run J)

These tables are meant to answer:
- Is the directional scatter mainly tied to the Quaia subspace?
- How do different Quaia scenarios shift the preferred core direction?

## 3) Invariance
- `master_invariance_axis_angles.csv` (Run Q) — stability of the operator axis across baseline vs 2MRS-clean and per pack.
- `master_invariance_RunP_sys_axis_comparison.csv` (Run P, if present)

## 4) Identifiability (“clumps”)
- `master_identifiability_thresholds_RunR.csv` — estimated injection amplitude where rank-2 becomes detectable.
- `master_identifiability_sweep_RunR.csv` — full sweep results.

All tables are derived from the canonical run outputs in `results/runs/`.
