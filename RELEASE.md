# RELEASE.md â€” Protocol D (publishable repo)

This checklist is for producing a clean, citable, reproducible public release (e.g., GitHub + Zenodo).

---

## 0) Decide release scope

- [ ] **Scope statement is explicit** (analysis layer over derived dipole packs; no raw-catalogue remeasurement).  
      See: `docs/03_data_inputs.md`
- [ ] **Public vs private data** is clearly separated:
  - public: `data/packs/`, `data/aux_data/`, `results/bundles/`, `results/runs/`, `results/summary/`
  - private (optional, not committed): `data/raw/` (if ever used)

---

## 1) Pre-release sanity checks (must pass)

- [ ] Fresh clone on a clean machine/VM
- [ ] Create venv and install requirements:
  ```bash
  python -m venv .venv
  . .venv/bin/activate
  pip install -r requirements.txt
  ```
- [ ] Reproduction smoke tests:
  ```bash
  make repro_bao
  make repro_operator
  ```
- [ ] Confirm outputs are created under `REPRO_OUT/` and the scripts run without manual edits.
- [ ] Confirm key tables exist:
  - `results/summary/master_fits.csv`
  - `results/summary/master_invariance_axis_angles.csv`
  - `results/summary/master_identifiability_thresholds_RunR.csv`
- [ ] Confirm the â€œpublishable checklistâ€ reads well:
  - `docs/07_key_claims.md`

---

## 2) Pin and document environment

Minimum standard:
- [ ] Ensure `requirements.txt` contains the packages needed for reproduction.
- [ ] Optional but recommended: record the exact environment used for the release:
  ```bash
  pip freeze > requirements-lock.txt
  python --version > python-version.txt
  ```
  Commit these two files for the release tag.

---

## 3) Update metadata for citation

- [ ] `CITATION.cff`:
  - [ ] version (e.g., 1.0.0)
  - [ ] date-released (YYYY-MM-DD)
  - [ ] repository-code (URL)
  - [ ] preferred-citation fields if you want a paper-style citation
- [ ] `LICENSE` present (code licence).  
- [ ] If you want a separate **data licence**, add `DATA_LICENSE.md` (optional).

---

## 4) Data availability statement (template)

Use in a paper/preprint:

> The analysis code and derived dipole input packs used in Protocol D are available in the public repository release.  
> The repository contains derived measurement packs, summary tables, and canonical run outputs (Runs Gâ€“R).  
> Raw sky catalogues and survey masks are not redistributed in this repository due to thirdâ€‘party licensing and estimator-dependence; users may reproduce the analysis layer from the provided derived packs, or reconstruct the full pipeline privately using upstream data under their respective licences.

(If you add `data/raw/` instructions later, link to them here.)

---

## 5) Create the release (GitHub)

- [ ] Merge to `main` (or your release branch)
- [ ] Tag:
  ```bash
  git tag -a v1.0.0 -m "Protocol D v1.0.0 (publishable release)"
  git push origin v1.0.0
  ```
- [ ] Create a GitHub Release for the tag:
  - Title: `Protocol D v1.0.0`
  - Attach the source zip (optional, GitHub does it automatically)
  - Include â€œKey claimsâ€ link (docs/07) and â€œMaster tablesâ€ link (docs/06)

---

## 6) Archive and mint DOI (Zenodo)

- [ ] Enable Zenodo-GitHub integration for the repository
- [ ] Trigger a release in GitHub (Zenodo picks it up)
- [ ] Confirm Zenodo record shows:
  - [ ] correct title and authors
  - [ ] version matches tag
  - [ ] licence
  - [ ] DOI minted
- [ ] Update `CITATION.cff` `repository-code` and add the DOI to your README (optional but recommended).

---

## 7) Post-release verification (the â€œsomeone elseâ€ test)

- [ ] A colleague (or â€œfuture youâ€) downloads the release archive and follows only:
  - `README.md`
  - `docs/06_master_summary.md`
  - `docs/07_key_claims.md`
- [ ] They can reproduce at least one run end-to-end (BAO and/or operator) without guessing.

---

## 8) Minimal release notes (suggested)

Include in GitHub Release:

- Whatâ€™s new in this version (e.g., master tables + key claims)
- Known limitations (no raw-catalogue remeasurement)
- How to verify (point to docs/07)

