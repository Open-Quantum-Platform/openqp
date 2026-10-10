"""The block CPHF solver and its one-right-hand-side-at-a-time fallback.

cphf_solve advances all closed-shell right-hand sides together by default and
solves them one at a time with OQP_CPHF_SERIAL=1.  Per right-hand side both
run the same preconditioned CG recurrences, so the analytic RHF and RKS
Hessians (3N right-hand sides) and the analytic IR/Raman tensors (three field
right-hand sides) must agree to well below the 1e-9 solver tolerance, and both
solvers must report every right-hand side converged.

Skipped unless the compiled OpenQP runtime is importable.
"""

import os
import re
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

INPUT = """[input]
system=
   8   0.000000000   0.000000000   0.117300000
   1   0.000000000   0.757200000  -0.469200000
   1   0.000000000  -0.757200000  -0.469200000
charge=0
runtype=hess
basis=6-31g*
method=hf
{dft}
[guess]
type=huckel
[scf]
type=rhf
multiplicity=1
conv=1.0e-10
[symmetry]
enabled=false
[hess]
type=analytical
state=0
"""


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "2")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


def _run(tmp, name, dft, serial):
    from oqp.pyoqp import Runner
    previous = os.environ.get("OQP_CPHF_SERIAL")
    os.environ["OQP_CPHF_SERIAL"] = "1" if serial else "0"
    try:
        inp = Path(tmp) / f"{name}.inp"
        inp.write_text(INPUT.format(dft=dft))
        log = Path(tmp) / f"{name}.log"
        runner = Runner(project=name, input_file=str(inp), log=str(log), silent=1, usempi=False)
        runner.run()
    finally:
        if previous is None:
            os.environ.pop("OQP_CPHF_SERIAL", None)
        else:
            os.environ["OQP_CPHF_SERIAL"] = previous
    mol = runner.mol
    tensors = [np.array(mol.data[tag], dtype=float).ravel()
               for tag in ("OQP::hf_dipole_derivatives", "OQP::hf_polarizability_derivatives")]
    return np.array(mol.hessian, dtype=float), tensors, log.read_text()


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class CphfBlockMatchesSerial(unittest.TestCase):
    def test_block_and_serial_solvers_agree(self):
        summary = re.compile(r"converged\s+(\d+) of\s+(\d+) right-hand sides")
        with tempfile.TemporaryDirectory() as tmp:
            for label, dft in (("rhf", ""), ("rks", "functional=b3lyp5")):
                with self.subTest(case=label):
                    h_blk, t_blk, log_blk = _run(tmp, f"{label}_block", dft, serial=False)
                    h_ser, t_ser, log_ser = _run(tmp, f"{label}_serial", dft, serial=True)
                    for log in (log_blk, log_ser):
                        counts = summary.findall(log)
                        self.assertTrue(counts, log[-2000:])
                        for done, total in counts:
                            self.assertEqual(done, total)
                    self.assertLess(np.abs(h_blk - h_ser).max(), 1.0e-8)
                    for a, b in zip(t_blk, t_ser):
                        self.assertLess(np.abs(a - b).max(), 1.0e-7 * max(1.0, np.abs(b).max()))


if __name__ == "__main__":
    unittest.main()
