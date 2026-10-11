"""Validate the open-shell (ROHF) CPHF solver via the static polarizability.

The ROHF CPHF solver (cphf_mod::cphf_solve_rohf) replicates the validated TRAH
ROHF orbital-Hessian action (scf_converger::calc_h_op) over the docc/socc/virt
rotation space (layout: rohf_pack_trial): per spin the orbital-energy-difference
part  Fvv x - x Foo  (full MO Fock blocks) plus the response Fock from the trial
rotation density (get_response_packed; scftype>=2 -> Coulomb from the spin-summed
density, exchange same-spin).

Correctness check (no external reference needed): for a closed-shell molecule run
as ROHF (multiplicity 1, offset = n_socc = 0) the rotation space reduces to the
virt-docc block and the ROHF orbital Hessian reduces to (twice) the RHF one, so
the ROHF static dipole polarizability must equal the validated closed-shell (RHF)
cphf_static_polarizability. A match confirms the solver, the operator and the
docc/socc/virt packing.

Skipped unless the compiled OpenQP runtime is importable.
"""

import os
import re
import tempfile
import unittest
from pathlib import Path
from unittest import mock


INPUT_TMPL = """[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g
method=hf
[guess]
type=huckel
[scf]
multiplicity=1
type={scftype}
"""


def _runtime_available():
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


def _isotropic(text):
    m = re.search(r"isotropic\s*=\s*(-?\d+\.\d+)", text)
    assert m, f"no isotropic value in:\n{text}"
    return float(m.group(1))


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class CphfRohfPolarizability(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="oqp_cphf_rohf_polar_")
        self.addCleanup(tmp.cleanup)
        self.tmpdir = Path(tmp.name)
        # The native self-test writes its report to $OQP_SELFTEST_DIR
        # (source/io/selftest_output.F90); keep it private to this test.
        env = mock.patch.dict(os.environ, {"OQP_SELFTEST_DIR": str(self.tmpdir)})
        env.start()
        self.addCleanup(env.stop)

    def _run(self, scftype, selftest_name, report_name):
        import oqp
        from oqp.pyoqp import Runner

        workdir = self.tmpdir / scftype
        workdir.mkdir()
        inp = workdir / "h2o.inp"
        inp.write_text(INPUT_TMPL.format(scftype=scftype))
        runner = Runner(project=f"h2o_{scftype}", input_file=str(inp),
                        log=str(workdir / "h2o.log"))
        runner.run()
        getattr(oqp, selftest_name)(runner.mol)
        out_path = self.tmpdir / report_name
        self.assertTrue(out_path.exists(), f"{selftest_name} produced no output")
        return out_path.read_text()

    def test_rohf_matches_rhf_on_closed_shell(self):
        rhf_text = self._run("rhf", "cphf_polarizability_selftest", "cphf_polar.out")
        rohf_text = self._run("rohf", "cphf_rohf_polarizability_selftest", "cphf_rohf_polar.out")
        self.assertAlmostEqual(
            _isotropic(rohf_text), _isotropic(rhf_text), places=4,
            msg=f"ROHF:\n{rohf_text}\nRHF:\n{rhf_text}",
        )


if __name__ == "__main__":
    unittest.main()
