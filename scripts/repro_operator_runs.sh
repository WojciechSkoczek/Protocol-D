#!/usr/bin/env bash
set -euo pipefail
STAMP=$(date +"%Y%m%d_%H%M%S")
OUT="REPRO_OUT/${STAMP}_operator"
mkdir -p "${OUT}"

python -m protocol_d.run_operator --config configs/defaults.yaml --out_dir "${OUT}"

echo "[OK] Operator reproduction done: ${OUT}"
