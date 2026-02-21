.PHONY: help setup repro_bao repro_operator

help:
	@echo "Targets:"
	@echo "  setup         Create venv + install requirements"
	@echo "  repro_bao      Reproduce BAO harness outputs"
	@echo "  repro_operator Reproduce operator O2 fits (structure of Runs K–Q)"

setup:
	python -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt

repro_bao:
	bash scripts/repro_bao_only.sh

repro_operator:
	bash scripts/repro_operator_runs.sh
