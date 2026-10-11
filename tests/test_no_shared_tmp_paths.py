"""Static guard: tests and native code must not use fixed paths under /tmp.

A fixed path such as /tmp/hess1_selftest.out or /tmp/oqp_cam_hess_rks_analytical
is shared by every user and every concurrent pytest run on a host.  A second
user gets PermissionError on a directory the first one created, and two suites
running at once overwrite each other's work directories and self-test reports.
Tests take a per-test temporary directory instead, and native self-tests write
their reports to $OQP_SELFTEST_DIR (source/io/selftest_output.F90), which the
Python drivers point at that directory.

Fixed /tmp strings that are only data -- a parameter value that is echoed back
or passed to a mock and never opened -- are listed in ALLOWED with the reason.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
THIS_FILE = Path(__file__).resolve()

# A quoted string literal that starts with /tmp/ (or is exactly "/tmp").
PY_TMP_LITERAL = re.compile(r"""['"]/tmp['"/]""")
# A Fortran string literal that starts with /tmp/, e.g. open(file='/tmp/x.out').
F90_TMP_LITERAL = re.compile(r"""['"]/tmp/""", re.IGNORECASE)

# (file relative to the repository root, exact substring of the line, reason)
ALLOWED = (
    ("tests/test_docker_distribution.py", 'Path("/tmp/owner.so")',
     "argument to a mocked subprocess.run (ldd); never opened"),
    ("tests/test_openqp_api.py", '"/tmp/minimal_hh.opdftb"',
     "DFTB parameter_path echoed back by the input builder; never opened"),
    ("tests/test_openqp_api.py", '"/tmp/gfn1.opxtb"',
     "xTB parameter_path echoed back by the input builder; never opened"),
    ("tests/test_openqp_dftb_abi1.py", '"/tmp/minimal.opdftb"',
     "parameter_path of a fake molecule driving a mocked native library; never opened"),
    ("tests/test_openqp_dftb_abi1.py", 'Path("/tmp/libopenqp_dftb_c.fake")',
     "library path handed to a mocked ctypes.CDLL; never loaded"),
    ("tests/test_openqp_dftb_logging.py", '"/tmp/minimal.opdftb"',
     "parameter_path in a config only checked by input_checker; never opened"),
    ("tests/test_openqp_dftb_schema_hooks.py", '"/tmp/minimal.opdftb"',
     "parameter_path in a config only checked by input_checker; never opened"),
    ("tests/test_openqp_xtb_schema_hooks.py", '"/tmp/gfn1.opxtb"',
     "parameter_path in a config only checked by input_checker; never opened"),
    ("tests/test_openqp_xtb_schema_hooks.py", '"/tmp/minimal.opdftb"',
     "parameter_path in a config only checked by input_checker; never opened"),
    ("tests/test_state_labels.py", 'executable="/tmp/openqp_dftb_probe"',
     "executable name only formatted into a settings summary"),
    ("tests/test_symmetry_metadata.py", "molecule.log = '/tmp/placeholder.log'",
     "log name on a stub Molecule; get_results() does not open it"),
)


def _python_files():
    return sorted(p for p in (ROOT / "tests").rglob("*.py") if p.resolve() != THIS_FILE)


def _fortran_files():
    files = []
    for top in (ROOT / "source", ROOT / "tests"):
        for pattern in ("*.F90", "*.f90"):
            files.extend(top.rglob(pattern))
    return sorted(set(files))


def _allowed(rel, line):
    return any(rel == f and sub in line for f, sub, _ in ALLOWED)


def _scan(files, pattern, comment):
    hits = []
    for path in files:
        rel = path.relative_to(ROOT).as_posix()
        for lineno, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            if line.lstrip().startswith(comment):
                continue
            if pattern.search(line) and not _allowed(rel, line):
                hits.append(f"{rel}:{lineno}: {line.strip()}")
    return hits


def test_python_tests_use_no_fixed_tmp_paths():
    hits = _scan(_python_files(), PY_TMP_LITERAL, "#")
    assert not hits, (
        "fixed /tmp paths in tests (use tmp_path / tempfile.TemporaryDirectory and "
        "OQP_SELFTEST_DIR; or add a data-only string to ALLOWED with a reason):\n"
        + "\n".join(hits))


def test_fortran_writes_no_fixed_tmp_paths():
    hits = _scan(_fortran_files(), F90_TMP_LITERAL, "!")
    assert not hits, (
        "fixed /tmp paths in Fortran (use selftest_output_path from "
        "source/io/selftest_output.F90):\n" + "\n".join(hits))


def test_every_allow_list_entry_still_matches():
    # A stale entry would silently allow a future fixed path in that file.
    stale = []
    for rel, sub, _ in ALLOWED:
        path = ROOT / rel
        if not path.is_file() or sub not in path.read_text(encoding="utf-8"):
            stale.append(f"{rel}: {sub}")
    assert not stale, "ALLOWED entries that match nothing:\n" + "\n".join(stale)
