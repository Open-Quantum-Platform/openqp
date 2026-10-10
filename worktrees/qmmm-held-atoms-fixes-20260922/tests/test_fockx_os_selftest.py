"""Integration test for the open-shell two-density derivative-Fock contraction.

Drives the ``fockx_os_selftest`` bind(C) harness, the open-shell analog of
``fockx_selftest``.  It checks the native contraction
``fock_deriv_mod::fock_deriv_contract_os``, which forms

    g_x = sum_uv M_uv ( J^x_uv[Pa+Pb] - c_x K^x_uv[Pa] ) = d/dx Tr[M . G^a]

(Coulomb from the total density, exchange from the spin density -- the spin-s
Fock that scf_addons::fock_jk builds for scftype>=2), against a central finite
difference of Tr[M . G^a] at frozen densities.  The reference is non-iterative,
so a PASS validates the open-shell exchange/Coulomb prefactors that the UHF
analytic Hessian response relies on.

Skipped unless the compiled OpenQP runtime is importable.
"""

import os
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELFTEST_OUT = Path("/tmp/fockx_os_selftest.out")

INPUT_TMPL = """[input]
system=
   8   0.000000000   0.000000000   0.000000000
   1   0.000000000   0.000000000   0.970000000
charge=0
runtype=energy
basis=6-31g*
method=hf
[guess]
type=huckel
[scf]
multiplicity=2
type={scftype}
"""


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class FockxOsSelfTest(unittest.TestCase):
    def _run_case(self, scftype):
        import oqp
        from oqp.pyoqp import Runner

        workdir = Path(f"/tmp/oqp_fockx_os_{scftype}")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "oh.inp"
        inp.write_text(INPUT_TMPL.format(scftype=scftype))

        if SELFTEST_OUT.exists():
            SELFTEST_OUT.unlink()

        runner = Runner(project=f"oh_{scftype}", input_file=str(inp),
                        log=str(workdir / "oh.log"))
        runner.run()
        oqp.fockx_os_selftest(runner.mol)

        self.assertTrue(SELFTEST_OUT.exists(), "self-test produced no output file")
        result = SELFTEST_OUT.read_text()
        self.assertIn("FOCKX_OS_SELFTEST PASS", result, msg=result)

    def test_uhf_contraction_matches_finite_difference(self):
        self._run_case("uhf")

    def test_rohf_contraction_matches_finite_difference(self):
        self._run_case("rohf")


if __name__ == "__main__":
    unittest.main()
