#!/usr/bin/env bash
# Reproduce every result PNG in results/ by running the demo notebook headlessly.
set -euo pipefail

uv run jupyter nbconvert \
  --to notebook \
  --execute notebooks/01_demo.ipynb \
  --output 01_demo_executed.ipynb \
  --ExecutePreprocessor.timeout=900
