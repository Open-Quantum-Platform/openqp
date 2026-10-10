"""Make the helpers beside the tests importable regardless of invocation path."""

import sys
from pathlib import Path

TEST_ROOT = Path(__file__).resolve().parent
DEVKIT_ROOT = TEST_ROOT.parent

# The integrated repository also has an engine-level tools/ directory. Put the
# devkit root first so imports such as tools.nac_lagrangian continue to resolve
# to devkit/tools exactly as they did in the standalone repository.
sys.path.insert(0, str(DEVKIT_ROOT))
sys.path.insert(0, str(TEST_ROOT))
