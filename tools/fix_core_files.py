from pathlib import Path
import subprocess

repo = Path(subprocess.check_output(["git","rev-parse","--show-toplevel"]).decode().strip())

def write(rel, content):
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write(content)

legacy = repo / ".gitattributes.txt"
if legacy.exists():
    legacy.unlink()

write("CITATION.cff", """cff-version: 1.2.0
message: "If you use this repository, please cite it."
title: "Protocol D"
type: software
version: "0.6.2"
date-released: "2026-02-21"
license: MIT
repository-code: "https://github.com/WojciechSkoczek/Protocol-D"
url: "https://github.com/WojciechSkoczek/Protocol-D"
authors:
  - family-names: Skoczek
    given-names: Wojciech
keywords:
  - cosmology
  - dipole
  - BAO
  - survey systematics
  - operator model
""")

write(".gitattributes", """* text=auto eol=lf

*.ps1 text eol=crlf
*.bat text eol=crlf
*.cmd text eol=crlf

*.png binary
*.jpg binary
*.jpeg binary
*.gif binary
*.pdf binary
*.zip binary
*.fits binary
*.npy binary
*.npz binary
*.pkl binary
*.parquet binary
""")

write("Makefile", """.PHONY: help setup repro_bao repro_operator

help:
\t@echo "Targets:"
\t@echo "  setup         Create venv + install requirements"
\t@echo "  repro_bao      Reproduce BAO harness outputs"
\t@echo "  repro_operator Reproduce operator O2 fits (structure of Runs K–Q)"

setup:
\tpython -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt

repro_bao:
\tbash scripts/repro_bao_only.sh

repro_operator:
\tbash scripts/repro_operator_runs.sh
""")

write("configs/defaults.yaml", """cosmology:
  H0: 70.0
  Om0: 0.3
  r_bao_mpc: 147.0
  bao_p: 1.0
  cmb_dipole_galactic_deg:
    l: 264.021
    b: 48.253

covariance:
  n_samples: 4000
  floor: 1.0e-10
  robust_loss:
    use_huber: true
    huber_delta: 2.5

operator_model_O2:
  t_definition: "t=1 for Quaia points, else 0"
  amplitude: "a0 + a_bmin*(bmin-30) + a_mag*(magcut-20.0)"
  default_quaia_bmin_for_zbins: 30

maskscan_error_floors:
  sigma_mag: 0.003
  sigma_dir_deg: 20.0
""")

targets = ["CITATION.cff","Makefile","configs/defaults.yaml",".gitattributes"]
for f in targets:
    b = (repo / f).read_bytes()
    print(f"{f}: LF={b.count(b'\\n')} CR={b.count(b'\\r')} BOM={b.startswith(b'\\xef\\xbb\\xbf')}")
print("OK: rewritten with UTF-8 (no BOM) and LF newlines.")