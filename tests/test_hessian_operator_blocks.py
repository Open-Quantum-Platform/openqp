"""Coordinate blocks of the blocked derivative-ERI operator.

eri_derivative_operator_mo holds the Cartesian operator for one block of
nuclear coordinates at a time; the block size follows
OQP_HESS_OPERATOR_MEM_MB.  The default budget keeps a small molecule in one
block, and OQP_HESS_OPERATOR_MEM_MB=0 forces one atom per block, i.e. one
derivative-ERI traversal per atom.  Both must give the same analytic RHF and
UHF Hessians (the UHF path builds J and K operators separately), with one and
with two OpenMP threads, and the blocked run must report its block layout.

Each case runs in a subprocess so OMP_NUM_THREADS and the budget are read by a
fresh runtime.  Skipped unless the compiled OpenQP runtime is importable.
"""

import json
import os
import subprocess
import sys
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
charge={charge}
runtype=hess
basis=6-31g*
method=hf
[guess]
type=huckel
[scf]
type={scf}
multiplicity={mult}
conv=1.0e-10
[symmetry]
enabled=false
[hess]
type=analytical
state=0
"""

CASES = {"rhf": dict(charge=0, scf="rhf", mult=1),
         "uhf": dict(charge=1, scf="uhf", mult=2)}

CHILD = """
import json, sys
import numpy as np
from oqp.pyoqp import Runner
inp, log, out = sys.argv[1:4]
runner = Runner(project='blk', input_file=inp, log=log, silent=1, usempi=False)
runner.run()
json.dump(np.asarray(runner.mol.hessian, dtype=float).tolist(), open(out, 'w'))
"""


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


def _hessian(tmp, case, threads, budget):
    tag = f"{case}_t{threads}_{'default' if budget is None else budget}"
    work = Path(tmp) / tag
    work.mkdir()
    inp = work / "blk.inp"
    inp.write_text(INPUT.format(**CASES[case]))
    env = {k: v for k, v in os.environ.items() if not k.startswith("OQP_")}
    env["OMP_NUM_THREADS"] = str(threads)
    if budget is not None:
        env["OQP_HESS_OPERATOR_MEM_MB"] = str(budget)
    log, out = work / "blk.log", work / "hess.json"
    proc = subprocess.run([sys.executable, "-c", CHILD, str(inp), str(log), str(out)],
                          cwd=work, env=env, capture_output=True, text=True, timeout=900)
    text = log.read_text() if log.exists() else ""
    if proc.returncode != 0 or not out.exists():
        raise AssertionError(f"{tag} failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}\n{text[-2000:]}")
    return np.array(json.loads(out.read_text())), text


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class OperatorCoordinateBlocks(unittest.TestCase):
    def test_atom_blocks_match_one_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            for case in CASES:
                ref, ref_log = _hessian(tmp, case, 1, None)
                self.assertNotIn("derivative operator:", ref_log)
                self.assertGreater(np.abs(ref).max(), 1.0e-2)
                for threads in (1, 2):
                    with self.subTest(case=case, threads=threads):
                        h, log = _hessian(tmp, case, threads, 0)
                        self.assertIn("derivative operator:      9 coordinates in blocks of     3", log)
                        self.assertLess(np.abs(h - ref).max(), 1.0e-8)


if __name__ == "__main__":
    unittest.main()
