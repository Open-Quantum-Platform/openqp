"""Finite-difference check of the xi^alpha Kohn-Sham matrix (tau replacement).

Drives the ``xi_fock_selftest`` bind(C) harness after a converged M06-2X SCF
with the fractional-derivative ingredient in the tau slot.  The harness
compares (E_xc[P+eps D]-E_xc[P-eps D])/2eps with sum F_ij D_ij for plain tau
and for xi^alpha; the two ratios must agree.  Skipped without a compiled
OpenQP runtime.
"""

import os
import unittest
from pathlib import Path

SELFTEST_OUT = Path("/tmp/xi_fock_selftest.out")

INPUT = """[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g
method=hf
functional=m06-2x
d4=False
[guess]
type=huckel
[scf]
multiplicity={mult}
type={scf}
conv=1.0e-8
maxit=60
[dftgrid]
rad_npts=48
ang_npts=110
pruned=
xi_mode=1
xi_alpha={alpha}
xi_p={p}
xi_scale={scale}
"""


def _runtime_available():
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return hasattr(oqp, "xi_fock_selftest")
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class XiFockSelfTest(unittest.TestCase):
    def _run(self, tag, scf, mult, alpha, p, scale):
        import oqp
        from oqp.pyoqp import Runner

        workdir = Path("/tmp/oqp_xi_fock_test")
        workdir.mkdir(exist_ok=True)
        inp = workdir / f"h2o_{tag}.inp"
        inp.write_text(INPUT.format(scf=scf, mult=mult, alpha=alpha, p=p, scale=scale))
        log = workdir / f"h2o_{tag}.log"
        if SELFTEST_OUT.exists():
            SELFTEST_OUT.unlink()
        runner = Runner(project=f"h2o_xi_{tag}", input_file=str(inp), log=str(log))
        runner.run()
        oqp.xi_fock_selftest(runner.mol)
        self.assertTrue(SELFTEST_OUT.exists(), "self-test produced no output file")
        result = SELFTEST_OUT.read_text()
        self.assertIn("XI_FOCK_SELFTEST PASS", result, msg=result)

    def test_rhf_alpha_half(self):
        self._run("rhf_a05", "rhf", 1, 0.5, -1, 0)

    def test_rhf_alpha_minus_half_scaled(self):
        self._run("rhf_am05", "rhf", 1, -0.5, -1, 1)

    def test_uhf_triplet_alpha_half(self):
        self._run("uhf_a05", "uhf", 3, 0.5, -1, 0)


if __name__ == "__main__":
    unittest.main()
