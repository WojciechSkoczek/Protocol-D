# RELEASE.md - Protocol D (publishable repo)

Checklist for producing a clean, citable, reproducible public release (GitHub + Zenodo).

## Pre-release sanity checks
- Fresh clone on a clean machine or VM.
- Create venv and install requirements.
- Smoke tests: `make repro_bao` and `make repro_operator`.
- Confirm outputs under `REPRO_OUT/`.

## Metadata
- `CITATION.cff`: version, date-released, repo URL.
- `LICENSE` is present and consistent with `CITATION.cff`.

## Release (GitHub)
- Create tag and push tags.
- Publish GitHub Release.

## Archive and DOI (Zenodo)
- Enable Zenodo-GitHub integration for this repo.
- Trigger a GitHub Release and confirm Zenodo record (version, license, DOI).
- Optionally add the DOI back to `README.md` and `CITATION.cff`.

## Windows (Git Bash) note
- If `bash scripts/*.sh` cannot find `python`, ensure your Windows Python is on PATH (Microsoft Store installs often live under `WindowsApps`). The scripts try to add `WindowsApps` automatically when running under Git Bash.
