#!/usr/bin/env python3
"""Index Protocol D harness run outputs into a single CSV.

This is a small bookkeeping helper intended to prevent the 'mixed folders' issue:
- it scans run-result *.txt files
- extracts key fitted quantities (g_hat, ΔAICc, beta)
- computes both CV(A_eff) and the BAO-anchored invariant CV(chi*A_eff)

Usage:
  python src/index_runs.py --runs_dir runs --out results/run_inventory.csv

You can point --runs_dir either to a folder containing run files, or to the
repo root; it will recursively search for `protocol_d_run_results_*.txt`.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np


RUN_GLOB = "protocol_d_run_results_*.txt"


def _float(m: Optional[re.Match], idx: int = 1) -> Optional[float]:
    if not m:
        return None
    try:
        return float(m.group(idx))
    except Exception:
        return None


def parse_run_text(text: str) -> Dict[str, object]:
    # Variant + input pack
    variant = None
    m = re.search(r"^variant:\s*(\S+)\s*$", text, flags=re.MULTILINE)
    if m:
        variant = m.group(1)

    input_pack = None
    m = re.search(r"^input:\s*(\S+)\s*$", text, flags=re.MULTILINE)
    if m:
        input_pack = m.group(1)

    # ΔAICc
    dAICc = _float(re.search(r"ΔAICc \(BAO_P1-null\)\s*=\s*([\-\d\.]+)", text))

    # g_hat
    m = re.search(r"g_hat \(l,b\) = \(([\d\.]+)°\,\s*([\d\.]+)°\)", text)
    g_l = _float(m, 1)
    g_b = _float(m, 2)

    # beta
    beta = _float(re.search(r"\bbeta\s*=\s*([\d\.Ee\+\-]+)", text))

    # Per-tracer chi and A_eff
    chis: List[float] = []
    aeffs: List[float] = []
    for line in text.splitlines():
        if ("chi=" in line) and ("A_eff=" in line):
            mm = re.search(r"chi=\s*([\d\.]+)\s+.*A_eff=\s*([\-\d\.Ee\+\-]+)", line)
            if mm:
                try:
                    chis.append(float(mm.group(1)))
                    aeffs.append(float(mm.group(2)))
                except Exception:
                    pass

    chis_np = np.array(chis, dtype=float) if chis else np.array([], dtype=float)
    aeff_np = np.array(aeffs, dtype=float) if aeffs else np.array([], dtype=float)

    cv_Aeff = None
    if len(aeff_np) > 1 and float(aeff_np.mean()) != 0.0:
        cv_Aeff = float(aeff_np.std(ddof=0) / aeff_np.mean() * 100.0)

    cv_chiAeff = None
    mean_chiAeff = None
    if len(aeff_np) > 1 and len(chis_np) == len(aeff_np):
        prod = chis_np * aeff_np
        mean = float(prod.mean())
        mean_chiAeff = mean
        if mean != 0.0:
            cv_chiAeff = float(prod.std(ddof=0) / mean * 100.0)

    return {
        "variant": variant,
        "input_pack": input_pack,
        "dAICc": dAICc,
        "g_l": g_l,
        "g_b": g_b,
        "beta": beta,
        "n_tracers": int(len(aeff_np)),
        "cv_Aeff_pct": cv_Aeff,
        "cv_chiAeff_pct": cv_chiAeff,
        "mean_chiAeff": mean_chiAeff,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs_dir", required=True, help="Folder to scan (recursive) for run outputs")
    ap.add_argument("--out", required=True, help="Output CSV path")
    args = ap.parse_args()

    runs_root = Path(args.runs_dir)
    files = sorted(runs_root.rglob(RUN_GLOB)) if runs_root.is_dir() else []

    rows: List[Dict[str, object]] = []
    for p in files:
        try:
            txt = p.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue
        row = parse_run_text(txt)
        row["run_file"] = p.name
        row["run_path"] = str(p)
        # simple tags to prevent accidental mixing
        inp = str(row.get("input_pack") or "")
        row["pack_type"] = "quaia_zbins" if "quaia" in inp else "catwise_only"
        row["dataset"] = "wagenveld" if "wagenveld" in inp else ("lotssswap" if "boehme" in inp else "unknown")
        rows.append(row)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        out_path.write_text("", encoding="utf-8")
        print(f"No run files found under: {runs_root}")
        return 0

    # stable column order
    cols = [
        "run_file",
        "run_path",
        "pack_type",
        "dataset",
        "variant",
        "input_pack",
        "g_l",
        "g_b",
        "dAICc",
        "beta",
        "n_tracers",
        "cv_Aeff_pct",
        "cv_chiAeff_pct",
        "mean_chiAeff",
    ]

    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c) for c in cols})

    print(f"Wrote {len(rows)} rows to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
