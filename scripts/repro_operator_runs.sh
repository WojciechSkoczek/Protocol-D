#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date +"%Y%m%d_%H%M%S")"
OUT="${ROOT}/REPRO_OUT/${STAMP}_operator"
mkdir -p "${OUT}"

# On Git Bash/MSYS, Microsoft Store Python is often under WindowsApps (not always on PATH).
if command -v cygpath >/dev/null 2>&1 && [ -n "${LOCALAPPDATA:-}" ]; then
  WINAPPS="$(cygpath -u "$LOCALAPPDATA")/Microsoft/WindowsApps"
  if [ -d "$WINAPPS" ]; then
    export PATH="$PATH:$WINAPPS"
  fi
fi

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

to_py_path() {
  local p="$1"
  if command -v cygpath >/dev/null 2>&1; then
    cygpath -m "$p"
  else
    echo "$p"
  fi
}

# Allow running module without installing the package (src/ layout)
export PYTHONPATH="${ROOT}/src"

"${PY}" "${PY_ARGS[@]}" -m protocol_d.run_operator   --config "$(to_py_path "${ROOT}/configs/defaults.yaml")"   --out_dir "$(to_py_path "${OUT}")"

echo "[OK] Operator reproduction done: ${OUT}"
