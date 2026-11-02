.PHONY: help install toy test lint pipeline clean

help:
	@echo "install   install package in editable mode"
	@echo "toy       regenerate the synthetic dataset"
	@echo "test      run the test suite"
	@echo "lint      run ruff"
	@echo "pipeline  run the full toy pipeline (needs hmmer/clustalo/iqtree)"

install:
	pip install -e ".[dev]"

toy:
	python -m genefamily toy --outdir data/toy

test:
	pytest -q

lint:
	ruff check src tests

pipeline: toy
	python -m genefamily run --config config/toy.yaml

clean:
	rm -rf results data/toy .pytest_cache .ruff_cache
