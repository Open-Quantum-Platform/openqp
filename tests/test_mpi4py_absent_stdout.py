"""Without mpi4py, ``import oqp`` must leave stdout to the caller.

``import oqp`` instantiates ``MPIManager``; without mpi4py it falls back to
serial mode and announces that.  The notice went to stdout, so every process
whose stdout is data -- ``python -c "import oqp...; print(json.dumps(...))"``
in tests/test_namd_droplet_boundary.py, for one -- failed with
JSONDecodeError on hosts without mpi4py.  The notice is a diagnostic and
belongs on stderr; the serial fallback itself is unchanged.

mpi4py is hidden by ``sys.modules["mpi4py"] = None``, which makes
``from mpi4py import MPI`` raise ImportError exactly as on a host without it.
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

from oqp_module_isolation import child_env

ROOT = Path(__file__).resolve().parents[1]
MPI_UTILS = ROOT / "pyoqp" / "oqp" / "utils" / "mpi_utils.py"
NOTICE = "Failed to import mpi4py"
SENTINEL = {"stdout": "data only"}

_HIDE_MPI4PY = 'import sys\nsys.modules["mpi4py"] = None\n'


def _run(body, hide_mpi4py=True):
    # child_env: a checkout without an installed oqp imports the source tree
    # through conftest's sys.path fallback, which a child does not inherit.
    return subprocess.run(
        [sys.executable, "-c", (_HIDE_MPI4PY if hide_mpi4py else "") + body],
        cwd=ROOT, env=child_env(), capture_output=True, text=True, check=False,
    )


def test_mpi_manager_serial_fallback_notice_goes_to_stderr():
    """Runs without the compiled backend: load mpi_utils.py by path."""
    result = _run(
        "import importlib.util, json\n"
        f"spec = importlib.util.spec_from_file_location('mpi_utils', {str(MPI_UTILS)!r})\n"
        "mod = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(mod)\n"
        "m = mod.MPIManager()\n"
        "assert not mod.mpi4py_available\n"
        "assert (m.rank, m.size, m.use_mpi, m.world_rank, m.world_size) == (0, 1, 0, 0, 1)\n"
        f"print(json.dumps({SENTINEL!r}))\n"
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(result.stdout) == SENTINEL
    assert NOTICE in result.stderr


def test_import_oqp_library_namd_keeps_stdout_parseable():
    """The droplet-smoke route: import oqp.library.namd, print JSON."""
    result = _run(
        "import json\n"
        "import oqp.library.namd\n"
        f"print(json.dumps({SENTINEL!r}))\n"
    )
    if result.returncode != 0:
        # Skip only when the import fails with mpi4py visible too, i.e. the
        # compiled backend itself is unavailable; otherwise hiding mpi4py
        # broke the import and that is a failure.
        probe = _run("import oqp.library.namd\n", hide_mpi4py=False)
        if probe.returncode != 0:
            last = (probe.stderr.strip().splitlines() or ["?"])[-1]
            pytest.skip(f"compiled OpenQP runtime is not available: {last}")
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(result.stdout) == SENTINEL
    assert NOTICE in result.stderr
