# Key claims and where to verify them

This repo is designed so that a reader can **verify the core Protocol D conclusions** without rerunning the entire pipeline from raw sky catalogues.

Below, each claim includes:
- **what it says**
- **where the evidence lives** (tables/columns)
- **a minimal verification step** (what to look for)

> Note: All claims are about the **analysis layer** (fits over *derived dipole packs*).  
> End-to-end remeasurement from raw catalogues is intentionally out of scope here (see `docs/03_data_inputs.md`).

---

## Claim 1 — Directional scatter is largely tied to the Quaia subspace (not universal across tracers)

**Statement:** Once the core BAO term is stabilised, the dominant directional scatter behaves as a tracer/subspace effect rather than a global cosmological axis drift.

**Verify in:**
- `results/summary/master_subspace_RunH.csv`

**What to check:**
- Look for systematic differences between:
  - "full set" fits vs "Quaia-only" (or "non-Quaia") subspace fits
  - direction / stability metrics reported in Run H

---

## Claim 2 — A low-rank “selection operator” (O2) captures Quaia-specific excess structure

**Statement:** Adding a Quaia-only operator term with a fixed axis and an amplitude depending on selection variables provides a compact description of the residual behaviour.

**Verify in:**
- `results/summary/master_fits.csv` (Run O provides “best model per dataset”)
- Canonical outputs: `results/runs/RunK/` and `results/runs/RunO/`

**What to check:**
- In `master_fits.csv`, filter `run == "O"` and confirm which model family is preferred per dataset.
- Cross-check the underlying per-run summaries in the Run K / Run O folders.

---

## Claim 3 — The operator axis is invariant across baseline vs 2MRS-clean (per pack)

**Statement:** The inferred operator axis \(\hat{s}\) is stable under the baseline → 2MRS-clean variant change, indicating a repeatable mapping distortion (operator), not a fragile artefact.

**Verify in:**
- `results/summary/master_invariance_axis_angles.csv` (Run Q)

**What to check:**
- The angle between \(\hat{s}_{\mathrm{baseline}}\) and \(\hat{s}_{\mathrm{2MRS}}\) is small per pack (see the angle column).
- Confirm the per-pack values in the canonical Run Q output folder.

---

## Claim 4 — Sharing one operator across packs and/or variants is preferred by information criteria

**Statement:** Joint fits that enforce a **shared operator axis and coefficients** are preferred relative to fitting each dataset independently (AICc improves).

**Verify in:**
- `results/summary/master_fits.csv`
  - Run P row(s): `model == "O2_joint_shared_operator"` with `dAICc` representing *joint minus separate*
  - Run Q row(s): `model == "O2_global_joint_shared_operator_all4"` with `dAICc` representing *global minus separate*

**What to check:**
- In `master_fits.csv`, locate these rows and verify:
  - `dAICc < 0` (joint/shared is better than separate)
  - shared operator axis columns (`s_l`, `s_b`, `s_to_CMB_deg`) are present for those rows

---

## Claim 5 — Rank-2 (“clumps”) becomes identifiable only above a non-trivial threshold

**Statement:** A genuinely independent second-axis component (rank-2) would need to exceed a threshold amplitude before the data would prefer a two-axis operator model.

**Verify in:**
- `results/summary/master_identifiability_thresholds_RunR.csv`
- `results/summary/master_identifiability_sweep_RunR.csv`

**What to check:**
- The threshold table reports the approximate injection amplitude where rank-2 becomes AICc-preferred.
- The sweep table shows the full ΔAICc curve across injected \(\lambda\).

---

## Claim 6 — Systematics scenarios shift preferred directions in predictable ways

**Statement:** Alternative Quaia constructions (systematics / conservative floors / mask scans) move the preferred core direction and operator fits in a structured way, consistent with an operator (map) effect.

**Verify in:**
- `results/summary/master_systematics_RunJ.csv`
- Canonical outputs: `results/runs/RunJ/`

**What to check:**
- Compare direction metrics and AICc across scenarios; look for consistent patterns rather than random scatter.

---

## Minimal “how to check quickly” (no code required)

Open these CSVs:
- `results/summary/master_fits.csv`
- `results/summary/master_invariance_axis_angles.csv`
- `results/summary/master_identifiability_thresholds_RunR.csv`

And verify:
1) shared-operator rows have `dAICc < 0`
2) the baseline vs 2MRS operator-axis angles are small
3) the rank-2 threshold is non-trivial (not ~0)

For the full context, read:
- `docs/02_models_and_assumptions.md`
- `docs/04_runs.md`
- `docs/06_master_summary.md`
