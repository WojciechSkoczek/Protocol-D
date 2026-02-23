#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date +"%Y%m%d_%H%M%S")"
OUT="${ROOT}/REPRO_OUT/${STAMP}"
mkdir -p "${OUT}"

# On Git Bash/MSYS, Microsoft Store Python is often under WindowsApps (not always on PATH).
if command -v cygpath >/dev/null 2>&1 && [ -n "${LOCALAPPDATA:-}" ]; then
  WINAPPS="$(cygpath -u "$LOCALAPPDATA")/Microsoft/WindowsApps"
  if [ -d "$WINAPPS" ]; then
    export PATH="$PATH:$WINAPPS"
  fi
fi

# Python chooser (Linux/macOS + Git Bash on Windows)
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
  echo "        On Windows: install Python and/or add it to PATH."
  exit 127
fi

echo "[INFO] Reproducing BAO harness outputs from packs -> ${OUT}"

PACKS_DIR="${ROOT}/data/packs"
SCRIPT="${ROOT}/src/protocol_d_run_input_pack_any_v2_bao.py"
INDEXER="${ROOT}/src/index_runs.py"

to_py_path() {
  local p="$1"
  if command -v cygpath >/dev/null 2>&1; then
    # C:/... (Windows-friendly, avoids backslash escaping issues)
    cygpath -m "$p"
  else
    echo "$p"
  fi
}

SCRIPT_PY="$(to_py_path "$SCRIPT")"
INDEXER_PY="$(to_py_path "$INDEXER")"
OUT_PY="$(to_py_path "$OUT")"

shopt -s nullglob
packs=( "${PACKS_DIR}"/*.csv )
if [ ${#packs[@]} -eq 0 ]; then
  echo "[ERROR] No packs found under: ${PACKS_DIR}"
  exit 2
fi

for pack in "${packs[@]}"; do
  pack_py="$(to_py_path "$pack")"
  for variant in baseline 2MRS_clean; do
    (
      cd "${OUT}"
      "${PY}" "${PY_ARGS[@]}" "${SCRIPT_PY}"         --input "${pack_py}"         --variant "${variant}"         --model bao_p1
    )
  done
done

"${PY}" "${PY_ARGS[@]}" "${INDEXER_PY}"   --runs_dir "${OUT_PY}"   --out "$(to_py_path "${OUT}/run_inventory.csv")"

echo "[OK] Done. Inspect ${OUT}"
