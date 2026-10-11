"""Validate the open-shell (UHF) CPHF solver via the static polarizability.

The open-shell CPHF solver (cphf_mod::cphf_solve_uhf) builds the UHF
orbital-Hessian action  (M U)^s_ia = (e^s_a - e^s_i) U^s_ia + [C^s^T dF^s C^s]_ia
with  dF^s = J[dPa+dPb] - cx K[dP^s]  (scf_addons::fock_jk, scftype>=2), i.e.
Coulomb from the spin-summed trial density and exchange from the same spin.

Correctness check (no external reference needed): for a closed-shell molecule
run as UHF with multiplicity 1, the converged Pa = Pb and the UHF CPHF static
dipole polarizability must equal the closed-shell (RHF) polarizability produced
by the already-validated cphf_static_polarizability. A match confirms the spin
coupling and the normalization (alpha_pq = -2 sum_s sum_ia mu^p,s_ia U^q,s_ia).

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
class CphfUhfPolarizability(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="oqp_cphf_polar_")
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

    def test_uhf_matches_rhf_on_closed_shell(self):
        rhf_text = self._run("rhf", "cphf_polarizability_selftest", "cphf_polar.out")
        uhf_text = self._run("uhf", "cphf_uhf_polarizability_selftest", "cphf_uhf_polar.out")
        rhf_iso = _isotropic(rhf_text)
        uhf_iso = _isotropic(uhf_text)
        self.assertAlmostEqual(
            uhf_iso, rhf_iso, places=4,
            msg=f"UHF iso={uhf_iso} vs RHF iso={rhf_iso}\nRHF:\n{rhf_text}\nUHF:\n{uhf_text}",
        )


if __name__ == "__main__":
    unittest.main()
