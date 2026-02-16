#!/usr/bin/env bash
set -euo pipefail

STAMP=$(date +"%Y%m%d_%H%M%S")
OUT="REPRO_OUT/${STAMP}"
mkdir -p "${OUT}"

echo "[INFO] Reproducing BAO harness outputs from packs -> ${OUT}"

python src/index_runs.py --packs_dir data/packs --out_dir "${OUT}"

echo "[OK] Done. Inspect ${OUT}"
