#!/usr/bin/env bash
# One-command pipeline: train on the post-fix OpenQP DB + contributions,
# distill the model, and sync the runtime copy used by converger_type=ml.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
TOOL_DIR=$(cd -- "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd -- "${TOOL_DIR}/../.." && pwd)

cd "${TOOL_DIR}"
python train/train_selector.py --db data/db_openqp_postfix.csv "data/contrib/*.csv" \
  --geom geometries --out model/scf_selector_model.py | tee model/metrics.txt

cp model/scf_selector_model.py "${REPO_ROOT}/pyoqp/oqp/library/scf_selector_model.py"
echo "model -> tools/scf-converger-ml/model/scf_selector_model.py"
echo "runtime copy -> pyoqp/oqp/library/scf_selector_model.py"
