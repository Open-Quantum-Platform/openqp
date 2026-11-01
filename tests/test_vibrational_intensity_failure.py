"""A failing displaced child run must stop the Hessian, not empty IR/Raman.

Hessian._compute_vibrational_intensities forms dipole and polarizability
derivatives for the TDDFT Hessian by central differences: each displacement is
a child Runner built from a copy of mol.config. The loop used to catch every
exception, record vibrational_intensity_metadata = {'status': 'failed'} and
return, so the run finished with exit status 0 and empty IR/Raman arrays. On
2026-10-17 that hid an input-schema rejection ("Unknown option: tdhf.cv4os")
of a key injected into mol.config; only the regression comparison of
H2_RHF_RPA_ANA_HESS noticed the missing intensities.

The test injects the same kind of key into the parent's config after it was
validated, so only the child Runner sees it, and requires the Hessian run to
raise with the child's own message. Before the fix runner.run() returns
normally and the test fails.

Skipped unless the compiled OpenQP runtime is importable.
"""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np


# H2 / STO-3G RPA analytic Hessian: the ground-state analytic dipole and
# polarizability derivatives are not used for a TDDFT state, so the
# finite-difference child-run path of the intensities is the one executed.
INPUT = """[input]
system=
  H  0.0000000000  0.0000000000 -0.3700000000
  H  0.0000000000  0.0000000000  0.3700000000
charge=0
method=tdhf
basis=sto-3g
runtype=hess
functional=

[scf]
type=rhf
multiplicity=1
conv=1.0e-10

[tdhf]
type=rpa
multiplicity=1
nstate=1
conv=1.0e-10
zvconv=1.0e-10

[hess]
type=analytical
state=1
clean=True
"""


def _runtime_available():
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class VibrationalIntensityFailure(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="oqp_vib_int_fail_")
        self.addCleanup(tmp.cleanup)
        self.tmpdir = Path(tmp.name)
        cwd = os.getcwd()
        os.chdir(self.tmpdir)
        self.addCleanup(os.chdir, cwd)

    def _runner(self, name):
        from oqp.pyoqp import Runner

        inp = self.tmpdir / f"{name}.inp"
        inp.write_text(INPUT)
        return Runner(project=name, input_file=str(inp),
                      log=str(self.tmpdir / f"{name}.log"), silent=1,
                      usempi=False)

    def test_child_config_error_is_raised(self):
        runner = self._runner("h2_vib_child_error")
        # Accepted by nothing: the parent already validated its input, so
        # only the child Runner built for the displaced geometries sees it.
        runner.mol.config["tdhf"]["cv4os"] = "true"
        with self.assertRaisesRegex(ValueError, r"Unknown option: tdhf\.cv4os"):
            runner.run()

    def test_failure_after_success_leaves_no_stale_intensities(self):
        # A caller that reuses the molecule and catches the exception must
        # not read the previous run's intensities as this run's result.
        runner = self._runner("h2_vib_rerun")
        runner.run()
        self.assertEqual(runner.mol.vibrational_intensity_metadata.get("status"), "computed")
        self.assertGreater(len(runner.mol.infrared_intensities), 0)
        runner.mol.config["tdhf"]["cv4os"] = "true"
        with self.assertRaisesRegex(ValueError, r"Unknown option: tdhf\.cv4os"):
            runner.run()
        mol = runner.mol
        self.assertEqual(mol.vibrational_intensity_metadata, {})
        self.assertEqual(mol.infrared_intensities.size, 0)
        self.assertEqual(mol.raman_activities.size, 0)
        self.assertEqual(mol.infrared_mode_dipole_derivatives.shape, (0, 3))
        self.assertEqual(mol.raman_mode_polarizability_derivatives.shape, (0, 3, 3))

    def test_clean_config_computes_intensities(self):
        # Control: the same run without the injected key must reach the
        # finite-difference backend and produce intensities, so the test
        # above fails for the injected key and not for the input.
        runner = self._runner("h2_vib_clean")
        runner.run()
        meta = runner.mol.vibrational_intensity_metadata
        self.assertEqual(meta.get("status"), "computed")
        self.assertEqual(meta.get("backend"), "native_openqp_finite_difference")
        self.assertTrue(np.all(np.isfinite(runner.mol.infrared_intensities)))
        self.assertGreater(len(runner.mol.infrared_intensities), 0)

    def test_mode_shape_mismatch_is_raised(self):
        from oqp.library.single_point import Hessian

        runner = self._runner("h2_vib_shape")
        runner.mol.config["input"]["runtype"] = "energy"
        runner.run()
        hess = Hessian(runner.mol)
        with self.assertRaisesRegex(ValueError, r"expected normal modes with shape"):
            hess._compute_vibrational_intensities(np.zeros((1, 5)))


if __name__ == "__main__":
    unittest.main()
