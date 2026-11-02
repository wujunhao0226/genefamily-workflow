#!/usr/bin/env bash
# End-to-end smoke test on the synthetic dataset.
# Set SKIP_EXTERNAL=1 to run only the pure-Python steps (no HMMER/Clustal/IQ-TREE).
set -euo pipefail

cd "$(dirname "$0")/.."

python -m genefamily toy --outdir data/toy

if [[ "${SKIP_EXTERNAL:-0}" == "1" ]]; then
  python -m genefamily run --config config/toy.yaml --skip-external
else
  python -m genefamily run --config config/toy.yaml
fi

echo "--- outputs ---"
find results -type f | sort | sed 's/^/  /'
