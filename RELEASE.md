# RELEASE.md - Protocol D (publishable repo)

Checklist for producing a clean, citable, reproducible public release (GitHub + Zenodo).

## Pre-release sanity checks
- Fresh clone on a clean machine/VM
- Create venv + install requirements
- Repro smoke tests: make repro_bao ; make repro_operator
- Confirm outputs under REPRO_OUT/

## Metadata
- CITATION.cff: version + date-released + repo URL
- LICENSE present and consistent with CITATION

## Release (GitHub)
- Tag + push tags
- Publish GitHub Release

## Archive + DOI (Zenodo)
- Enable Zenodo-GitHub integration for this repo
- Trigger a GitHub Release and confirm Zenodo record (version, license, DOI)
