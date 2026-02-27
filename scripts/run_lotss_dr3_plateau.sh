#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

FITS="${ROOT}/data/raw/LoTSS_DR3_v1.0.srl.fits"
OUT="${ROOT}/results/lotss_dr3/plateau_dr3_counts.csv"

mkdir -p "$(dirname "$FITS")" "$(dirname "$OUT")"

# Python chooser (Windows Git Bash + Linux/macOS)
PY="python"
PY_ARGS=()
if command -v python >/dev/null 2>&1; then
  PY="python"
elif command -v python3 >/dev/null 2>&1; then
  PY="python3"
elif command -v py >/dev/null 2>&1; then
  PY="py"
  PY_ARGS=(-3)
else
  echo "[ERROR] Python not found in PATH (tried: python, python3, py)."
  exit 127
fi

# Optional: quick dependency check (astropy-healpix)
"${PY}" "${PY_ARGS[@]}" - << 'PY'
try:
    from astropy_healpix import healpy as _hp
    import astropy
except Exception as e:
    raise SystemExit("Missing deps. Install: pip install astropy astropy-healpix") from e
print("OK: astropy + astropy-healpix available")
PY

# Download catalogue once (big file)
if [ ! -f "$FITS" ]; then
  URL="https://lofar-surveys.org/public/DR3/catalogues/LoTSS_DR3_v1.0.srl.fits"
  echo "[INFO] Downloading DR3 source catalogue to: $FITS"
  if command -v wget >/dev/null 2>&1; then
    wget -O "$FITS" "$URL"
  elif command -v curl >/dev/null 2>&1; then
    curl -L -o "$FITS" "$URL"
  else
    echo "[ERROR] Need wget or curl to download the FITS file."
    echo "Download manually from:"
    echo "  $URL"
    exit 2
  fi
fi

export PYTHONPATH="${ROOT}/src"

# NOTE:
# - footprint_flux_mjy is set LOW (10 mJy), so the footprint mask reflects survey geometry,
#   not only the sparsest ultra-bright sources.
# - flux plateau excludes 80 mJy (can be re-added later once footprint is independent).
"${PY}" "${PY_ARGS[@]}" "${ROOT}/src/lotss_dr3_plateau.py" \
  --fits "$FITS" \
  --out "$OUT" \
  --nside 64 \
  --flux_cuts_mjy "10,20,30,40,60" \
  --bmins_deg "10,20,30" \
  --footprint_flux_mjy 10

echo "[OK] Output: $OUT"