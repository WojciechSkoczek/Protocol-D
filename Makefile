# Protocol D Makefile
.RECIPEPREFIX := >
.PHONY: help setup repro_bao repro_operator

help:
>echo "Targets:"
>echo "  setup           Create venv (you still need to activate it) + show install command"
>echo "  repro_bao       Reproduce BAO harness outputs"
>echo "  repro_operator  Reproduce operator O2 fits"

setup:
>python -m venv .venv
>@echo "Activate the venv and run: pip install -r requirements.txt"

repro_bao:
>bash scripts/repro_bao_only.sh

repro_operator:
>bash scripts/repro_operator_runs.sh
