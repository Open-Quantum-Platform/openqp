"""Validate the SCF-type dispatch of the cphf_static_polarizability entry point.

The C entry point cphf_static_polarizability (used by the vibrational
Raman-activity path in single_point._compute_vibrational_intensities) must
route to the matching CPHF response kernel for the converged reference:
RHF/RKS -> closed-shell, UHF/UKS -> per-spin UHF, ROHF/ROKS -> single-MO ROHF.
Previously it always ran the closed-shell kernel, silently producing wrong
Raman activities for open-shell Hessian runs.

Correctness check (no external reference needed): for a closed-shell molecule
run as RHF, UHF (multiplicity 1) and ROHF (multiplicity 1), all three dispatched
tensors must agree, and the open-shell routes must reproduce the validated
closed-shell result.

Skipped unless the compiled OpenQP runtime is importable.
"""

import os
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

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
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class PolarizabilityDispatch(unittest.TestCase):
    def _alpha(self, scftype):
        import oqp
        from oqp.pyoqp import Runner

        workdir = Path(f"/tmp/oqp_polar_dispatch_{scftype}")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "h2o.inp"
        inp.write_text(INPUT_TMPL.format(scftype=scftype))
        runner = Runner(project=f"h2o_disp_{scftype}", input_file=str(inp),
                        log=str(workdir / "h2o.log"))
        runner.run()
        alpha = np.zeros((3, 3), dtype=np.float64)
        oqp.cphf_static_polarizability(
            runner.mol, oqp.ffi.cast("double *", oqp.ffi.from_buffer(alpha)))
        return alpha

    def test_dispatch_matches_closed_shell(self):
        a_rhf = self._alpha("rhf")
        a_uhf = self._alpha("uhf")
        a_rohf = self._alpha("rohf")
        # sanity: a real, nonzero tensor
        self.assertGreater(np.trace(a_rhf) / 3.0, 0.1)
        np.testing.assert_allclose(a_uhf, a_rhf, atol=2.0e-4,
                                   err_msg="UHF-dispatched polarizability != RHF")
        np.testing.assert_allclose(a_rohf, a_rhf, atol=2.0e-4,
                                   err_msg="ROHF-dispatched polarizability != RHF")


if __name__ == "__main__":
    unittest.main()
