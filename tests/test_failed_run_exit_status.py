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
  input whose Z-vector cannot converge (``maxit_zv=1``, ``zvconv=1e-30``) and on
  an energy input whose TD Davidson cannot converge (``maxit=1``,
  ``conv=1e-30``), and asserts a non-zero status, the log marker, and the
  raised error.
* ``SameProcessRecoveryTests`` runs, in ONE interpreter, a failing gradient,
  a good gradient, a failing TD Davidson (larger basis), a good TD energy and
  the good gradient again, and requires the good results to equal fresh-process
  runs bit for bit (at one thread): a caught failure must leave native and
  Python state reusable, with the basis growing and shrinking between calls.
* ``SetBasisFailureTests`` drives ``set_basis`` with the native basis loader
  mocked to report ``basis_set_issue`` and asserts a ``RuntimeError``, not a
  ``SystemExit``.

The QM/MM Z-vector and link-atom branches need OpenMM and are covered by the
source scan only.
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



# cc-pVDZ, not STO-3G: the TD Davidson auto-restarts twice with a doubled
# subspace (maxit 1 -> 2 -> 4).  The STO-3G symmetry blocks hold only one or
# two singles, so the subspace becomes complete and the solve is exact even at
# conv=1e-30; with cc-pVDZ (95 singles) the final residual stays ~1e-7.
DAVIDSON_TDA_INPUT = """[input]
system=
  8    0.000000000    0.000000000   -0.041061554
  1   -0.533194329    0.533194329   -0.614469223
  1    0.533194329   -0.533194329   -0.614469223
charge=0
method=tdhf
basis=cc-pvdz
runtype=energy
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

[tdhf]
type=tda
multiplicity=1
nstate=2
maxit=1
conv=1.0e-30
"""


GOOD_GRAD_INPUT = (GRAD_RPA_INPUT.replace("maxit_zv=1\n", "maxit_zv=50\n")
                   .replace("zvconv=1.0e-30\n", "zvconv=1.0e-10\n"))
GOOD_TDA_INPUT = (DAVIDSON_TDA_INPUT.replace("maxit=1\n", "maxit=50\n")
                  .replace("conv=1.0e-30\n", "conv=1.0e-10\n"))

# Driver executed in a fresh interpreter: runs each (name, expect) step through
# oqp.pyoqp.Runner in the SAME process and prints one JSON line of results.
_RUNNER_DRIVER = r"""
import json, sys, traceback
import numpy as np
from oqp.pyoqp import Runner
out = []
for name, expect in json.loads(sys.argv[1]):
    rec = {"name": name, "expect": expect}
    try:
        r = Runner(project=name, input_file=name + ".inp", log=name + ".log",
                   silent=1, usempi=False)
        r.run()
        rec["energies"] = [float(e) for e in np.ravel(r.mol.energies)]
        grads = getattr(r.mol, "grads", None)
        rec["grads"] = None if grads is None else [float(g) for g in np.ravel(grads)]
        rec["error"] = None
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {exc}"
        rec["trace"] = traceback.format_exc()[-2000:]
    out.append(rec)
print("RESULT " + json.dumps(out))
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

    def _assert_failure(self, proc, log_text, log_marker, error_text):
        detail = (f"returncode={proc.returncode}\n--- stdout ---\n{proc.stdout[-3000:]}"
                  f"\n--- stderr ---\n{proc.stderr[-3000:]}")
        # The run must fail *because of* the intended branch, not another reason.
        self.assertIn(log_marker, log_text, detail)
        self.assertIn(error_text, proc.stderr, detail)
        self.assertNotEqual(proc.returncode, 0, detail)

    def test_tddft_grad_zvector_failure_is_nonzero(self):
        proc, log_text = self._run("h2o_rpa_grad_zv_fail", GRAD_RPA_INPUT)
        self._assert_failure(proc, log_text, "PyOQP: TD Z-vector is not converged",
                             "TD Z-vector did not converge")

    def test_td_davidson_failure_is_nonzero(self):
        proc, log_text = self._run("h2o_tda_davidson_fail", DAVIDSON_TDA_INPUT)
        self._assert_failure(proc, log_text, "PyOQP: TD energy is not converged",
                             "TD Davidson did not converge")


@unittest.skipUnless(RUNTIME_AVAILABLE, "compiled OpenQP runtime not available")
class SameProcessRecoveryTests(unittest.TestCase):
    INPUTS = {
        "grad_fail": GRAD_RPA_INPUT,
        "grad_good": GOOD_GRAD_INPUT,
        "dav_fail": DAVIDSON_TDA_INPUT,
        "tda_good": GOOD_TDA_INPUT,
    }

    def _drive(self, work, steps):
        import json
        env = dict(os.environ)
        # One thread on purpose: with several threads the reduction order makes
        # two FRESH processes differ by up to ~6e-9 Eh (measured at 4 threads),
        # which would hide a small state leak.  At one thread fresh runs, repeated
        # runs and runs after a caught failure are bit-identical, so exact
        # equality is the contract checked here.
        env["OMP_NUM_THREADS"] = "1"
        proc = subprocess.run(
            [sys.executable, "-c", _RUNNER_DRIVER, json.dumps(steps)],
            cwd=work, env=env, capture_output=True, text=True, timeout=900)
        lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
        detail = (f"returncode={proc.returncode}\n--- stdout ---\n{proc.stdout[-3000:]}"
                  f"\n--- stderr ---\n{proc.stderr[-3000:]}")
        self.assertEqual(proc.returncode, 0, detail)
        self.assertEqual(len(lines), 1, detail)
        return {(i, r["name"]): r for i, r in enumerate(json.loads(lines[0][7:]))}

    def _workdir(self, tag):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        work = Path(tmp.name)
        for name, text in self.INPUTS.items():
            # A unique project name per step keeps logs apart.
            (work / f"{name}.inp").write_text(text)
        return work

    def assertSame(self, got, ref, what):
        self.assertIsNotNone(got, what)
        self.assertEqual(len(got), len(ref), what)
        diff = max(abs(a - b) for a, b in zip(got, ref))
        self.assertEqual(got, ref, f"{what}: not bit-identical, max |diff| = {diff:.3e}")

    def test_failed_runs_leave_the_process_reusable(self):
        # Fresh-process references, one interpreter per good input.
        ref_grad = self._drive(self._workdir("ref_grad"), [["grad_good", "ok"]])[(0, "grad_good")]
        ref_tda = self._drive(self._workdir("ref_tda"), [["tda_good", "ok"]])[(0, "tda_good")]
        for ref in (ref_grad, ref_tda):
            self.assertIsNone(ref["error"], ref.get("trace"))

        # One interpreter: fail, succeed, fail (larger basis), succeed, repeat.
        seq = [["grad_fail", "zv"], ["grad_good", "ok"], ["dav_fail", "dav"],
               ["tda_good", "ok"], ["grad_good", "ok"]]
        res = self._drive(self._workdir("seq"), seq)
        self.assertIn("RuntimeError: TD Z-vector did not converge",
                      res[(0, "grad_fail")]["error"] or "")
        self.assertIn("RuntimeError: TD Davidson did not converge",
                      res[(2, "dav_fail")]["error"] or "")
        for key in ((1, "grad_good"), (3, "tda_good"), (4, "grad_good")):
            self.assertIsNone(res[key]["error"], res[key].get("trace"))
        for key in ((1, "grad_good"), (4, "grad_good")):
            self.assertSame(res[key]["energies"], ref_grad["energies"], f"{key} energies")
            self.assertSame(res[key]["grads"], ref_grad["grads"], f"{key} gradient")
        self.assertSame(res[(3, "tda_good")]["energies"], ref_tda["energies"],
                        "TDA energies after failures")


@unittest.skipUnless(RUNTIME_AVAILABLE, "compiled OpenQP runtime not available")
class SetBasisFailureTests(unittest.TestCase):
    def _call(self, issue):
        import importlib
        from unittest import mock
        # oqp.library re-exports the set_basis *function* under the module's
        # name, so import the module itself.
        sb = importlib.import_module("oqp.library.set_basis")

        class FakeMol:
            config = {"input": {"basis": "sto-3g"}}
            data = {"basis_set_issue": issue}

        mol = FakeMol()
        with mock.patch.object(sb, "BasisData"), \
                mock.patch.object(sb.oqp, "apply_basis"), \
                mock.patch.object(sb, "dump_log") as dump_log:
            try:
                sb.set_basis(mol)
            except SystemExit as exc:  # the defect: status-0 SystemExit
                self.fail(f"set_basis raised SystemExit({exc.code!r})")
        return dump_log

    def test_basis_set_issue_raises_runtime_error(self):
        with self.assertRaisesRegex(RuntimeError, "basis set 'sto-3g' is not set properly"):
            self._call(1)

    def test_no_basis_set_issue_returns_normally(self):
        dump_log = self._call(0)
        dump_log.assert_not_called()



if __name__ == "__main__":
    unittest.main()
