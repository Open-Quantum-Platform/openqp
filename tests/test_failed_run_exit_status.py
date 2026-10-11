"""Regression guard: a failed PyOQP run must end with a non-zero exit status.

The TD Z-vector non-convergence branch of ``Gradient.tddft_grad`` (and the TD
Davidson branch of ``excitation``, the QM/MM Z-vector branch, the basis-set
check and the QM/MM link-atom check) stopped the run with a bare ``exit()``.  ``exit()`` raises ``SystemExit(None)``,
which the interpreter reports as status 0, so a gradient whose Z-vector never
converged looked like a successful run to every shell script, scheduler and
workflow manager that checks ``$?``.

Two guards:

* ``BareExitSourceTests`` parses every module under ``pyoqp/oqp`` and fails on
  any ``exit()``/``quit()``/``sys.exit()``/``os._exit()`` call with no argument
  or a literal ``None``/``0`` argument -- the root cause, independent of which
  branch happens to be reached.
* ``ZVectorFailureExitStatusTests`` runs ``python -m oqp.pyoqp`` on a gradient
  input whose Z-vector cannot converge (``maxit_zv=1``, ``zvconv=1e-30``) and
  asserts a non-zero status, the log marker, and the raised error.
"""
import ast
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "pyoqp" / "oqp"

NATIVE_TESTS_REQUIRED = os.getenv("OQP_REQUIRE_NATIVE_TESTS") == "1"
_RUNTIME_ERROR = ""

EXIT_NAMES = {"exit", "quit"}
EXIT_ATTRS = {("sys", "exit"), ("os", "_exit")}

GRAD_RPA_INPUT = """[input]
system=
  8    0.000000000    0.000000000   -0.041061554
  1   -0.533194329    0.533194329   -0.614469223
  1    0.533194329   -0.533194329   -0.614469223
charge=0
method=tdhf
basis=sto-3g
runtype=grad
functional=
d4=False

[guess]
type=huckel
save_mol=False

[scf]
type=rhf
multiplicity=1
conv=1.0e-8
save_molden=False

[properties]
grad=1

[tdhf]
type=rpa
multiplicity=1
nstate=2
conv=1.0e-8
maxit_zv=1
zvconv=1.0e-30
"""



def _silent_exit_calls(path):
    """(line, source) of exit calls that report status 0 on a failure path."""
    tree = ast.parse(path.read_text(encoding="utf-8"), str(path))
    hits = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if isinstance(func, ast.Name):
            is_exit = func.id in EXIT_NAMES
        elif isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name):
            is_exit = (func.value.id, func.attr) in EXIT_ATTRS
        else:
            is_exit = False
        if not is_exit or node.keywords:
            continue
        if not node.args:
            hits.append((node.lineno, ast.unparse(node)))
        elif (len(node.args) == 1 and isinstance(node.args[0], ast.Constant)
              and node.args[0].value in (None, 0, False)):
            hits.append((node.lineno, ast.unparse(node)))
    return hits


class BareExitSourceTests(unittest.TestCase):
    def test_scanner_detects_bare_exit(self):
        with tempfile.TemporaryDirectory() as tmp:
            probe = Path(tmp) / "probe.py"
            probe.write_text("import sys, os\n"
                             "exit()\nquit(None)\nsys.exit(0)\nos._exit(0)\n"
                             "sys.exit(1)\nexit('message')\nsys.exit(status)\n")
            self.assertEqual([line for line, _ in _silent_exit_calls(probe)],
                             [2, 3, 4, 5])

    def test_no_status_zero_exit_in_package(self):
        offenders = []
        for path in sorted(PKG.rglob("*.py")):
            for line, src in _silent_exit_calls(path):
                offenders.append(f"{path.relative_to(ROOT)}:{line}: {src}")
        self.assertEqual(offenders, [],
                         "exit calls that end a failed run with status 0; "
                         "raise an exception instead:\n" + "\n".join(offenders))


def _runtime_available():
    global _RUNTIME_ERROR
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception as exc:
        _RUNTIME_ERROR = f"{type(exc).__name__}: {exc}"
        return False


RUNTIME_AVAILABLE = _runtime_available()


class FailedRunExitStatusBuildGate(unittest.TestCase):
    """Turn a native-import failure into a FAILURE under CI, not a silent skip."""

    @unittest.skipUnless(
        NATIVE_TESTS_REQUIRED,
        "source-only run; set OQP_REQUIRE_NATIVE_TESTS=1 after building OpenQP",
    )
    def test_required_runtime(self):
        self.assertTrue(RUNTIME_AVAILABLE,
                        f"compiled OpenQP runtime is required: {_RUNTIME_ERROR}")


@unittest.skipUnless(RUNTIME_AVAILABLE, "compiled OpenQP runtime not available")
class ZVectorFailureExitStatusTests(unittest.TestCase):
    def _run(self, name, text):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        work = Path(tmp.name)
        inp = work / f"{name}.inp"
        inp.write_text(text)
        env = dict(os.environ)
        env.setdefault("OMP_NUM_THREADS", "1")
        proc = subprocess.run(
            [sys.executable, "-m", "oqp.pyoqp", inp.name, "--nompi", "--silent"],
            cwd=work, env=env, capture_output=True, text=True, timeout=600)
        log = work / f"{name}.log"
        log_text = log.read_text(errors="replace") if log.exists() else ""
        return proc, log_text

    def _assert_zvector_failure(self, proc, log_text):
        detail = (f"returncode={proc.returncode}\n--- stdout ---\n{proc.stdout[-3000:]}"
                  f"\n--- stderr ---\n{proc.stderr[-3000:]}")
        # The run must fail *because of* the Z-vector, not for another reason.
        self.assertIn("PyOQP: TD Z-vector is not converged", log_text, detail)
        self.assertIn("TD Z-vector did not converge", proc.stderr, detail)
        self.assertNotEqual(proc.returncode, 0, detail)

    def test_tddft_grad_zvector_failure_is_nonzero(self):
        proc, log_text = self._run("h2o_rpa_grad_zv_fail", GRAD_RPA_INPUT)
        self._assert_zvector_failure(proc, log_text)



if __name__ == "__main__":
    unittest.main()
