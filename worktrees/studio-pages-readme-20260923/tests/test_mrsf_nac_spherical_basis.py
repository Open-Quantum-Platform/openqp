"""The analytic MRSF NAC must respect molecular symmetry in a SPHERICAL basis.

Why this test exists
--------------------
grd2 drives the two-particle-density digests with CARTESIAN shell extents.
A digest that indexes spherical densities with those extents reads the wrong
elements. That is invisible in a Cartesian basis, where the two index spaces
coincide, so every analytic-NAC example and test in this repository -- all of
which used 6-31G or 6-31G* -- passed while the spherical path was wrong:

  * d functions (5 vs 6): neighbouring AO elements, wrong values, broken
    molecular symmetry, and a result that moved with the thread count because
    it amplified the 1e-13 SCF reduction-order noise to O(1).
  * f functions (7 vs 10): writes past the end of the array, heap corruption,
    SIGABRT.

The detector here is symmetry, so it needs no external reference and no
finite-difference run. The H2O geometry below has a C2 axis along z, and the
two hydrogens are exchanged by it. A derivative coupling between two states of
definite symmetry must therefore reproduce itself, up to one overall sign,
under "rotate every vector by C2z and swap the two hydrogens". A misindexed
digest scrambles AO orientation and breaks that immediately.

Both AO conventions are checked: 6-31G* is Cartesian and is the control that
passed even when the spherical path was broken; cc-pVDZ is spherical.
"""
import os
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

try:
    import oqp                      # noqa: F401  (must precede numpy: ILP64)
    import numpy as np
    from oqp.pyoqp import Runner
    HAVE_OQP = True
except Exception:                   # pragma: no cover - source-only checkout
    HAVE_OQP = False


# O on the C2 axis; the two H are exchanged by C2z: (x,y,z) -> (-x,-y,z).
GEOMETRY = """   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223"""

INPUT = """[input]
system=
{geometry}
charge=0
runtype=nac
basis={basis}
functional=bhhlyp
method=tdhf

[guess]
type=huckel

[scf]
type=rohf
multiplicity=3
conv=1.0e-10
maxit=200

[tdhf]
type=mrsf
multiplicity=1
nstate={nstate}
conv=1.0e-10
zvconv=1.0e-10
maxit=100

[nac]
type=analytical
states=1 2
"""


NSTATE = 3


def _c2z_then_swap_hydrogens(d):
    """Apply C2z to every atom's vector and exchange the two hydrogens."""
    rotated = d * np.array([-1.0, -1.0, 1.0])
    return np.array([rotated[0], rotated[2], rotated[1]])


@unittest.skipUnless(HAVE_OQP, "compiled oqp package not importable")
class MRSFNACSphericalBasisTests(unittest.TestCase):

    def _derivative_coupling(self, basis, nstate=NSTATE, pair=(1, 2)):
        with tempfile.TemporaryDirectory(prefix="nac_spherical_") as tmp:
            path = os.path.join(tmp, "h2o.inp")
            with open(path, "w") as handle:
                handle.write(INPUT.format(geometry=GEOMETRY, basis=basis,
                                          nstate=nstate))
            runner = Runner(input_file=path, log=path.replace(".inp", ".log"))
            runner.run()
            # OQP::nac_dcv is reserved as (3*natom, nstate, nstate) in Fortran
            # order. Check every retained pair before selecting one for the
            # molecular-symmetry test.
            dcv = np.array(runner.mol.data["OQP::nac_dcv"], copy=True)
            self.assertEqual(dcv.size, 9 * nstate * nstate)
            dcv = dcv.reshape(-1).reshape((9, nstate, nstate), order="F")
            self.assertTrue(np.all(np.isfinite(dcv)))
            np.testing.assert_allclose(dcv, -dcv.swapaxes(1, 2), atol=1e-12)
        return dcv[:, pair[0], pair[1]].reshape(-1, 3)

    def _assert_c2v(self, basis, nstate=NSTATE, pair=(1, 2)):
        d = self._derivative_coupling(basis, nstate, pair)
        self.assertTrue(np.all(np.isfinite(d)),
                        f"{basis}: derivative coupling is not finite:\n{d}")
        scale = np.abs(d).max()
        self.assertGreater(scale, 1e-8,
                           f"{basis}: derivative coupling is identically zero")

        image = _c2z_then_swap_hydrogens(d)
        symmetric = np.abs(image - d).max()
        antisymmetric = np.abs(image + d).max()
        residual = min(symmetric, antisymmetric) / scale
        self.assertLess(
            residual, 1e-6,
            f"{basis}: the analytic derivative coupling does not transform "
            f"under C2z as either representation of the C2v point group "
            f"(relative residual {residual:.3e}).\n"
            f"coupling:\n{d}\nC2z image with the hydrogens swapped:\n{image}")

    def test_cartesian_basis_respects_the_c2_axis(self):
        # Control: this passed even while the spherical path was broken.
        self._assert_c2v("6-31g*")

    def test_spherical_basis_respects_the_c2_axis(self):
        self._assert_c2v("cc-pvdz")

    def test_two_state_small_basis_allocation(self):
        self._assert_c2v("6-31g", nstate=2, pair=(0, 1))

    def test_four_state_spherical_allocation(self):
        self._assert_c2v("cc-pvdz", nstate=4, pair=(0, 1))


    def test_repeated_calls_restore_state_and_clear_unselected_pairs(self):
        from oqp.library.nac_analytic import analytic_nac

        # Grow and shrink the response and AO dimensions in one process.
        for basis, nstate in (("6-31g", 2), ("cc-pvdz", 4), ("6-31g", 2)):
            with self.subTest(basis=basis, nstate=nstate):
                with tempfile.TemporaryDirectory(prefix="nac_lifetime_") as tmp:
                    path = Path(tmp) / "h2o.inp"
                    path.write_text(INPUT.format(geometry=GEOMETRY, basis=basis,
                                                 nstate=nstate))
                    runner = Runner(input_file=str(path), log=str(path.with_suffix(".log")))
                    runner.run()
                    mol = runner.mol
                    original = np.array(mol.data["OQP::td_bvec_mo"], copy=True)
                    cutoff = mol.data._data.control.int2e_cutoff
                    full_h, full_d = analytic_nac(mol)
                    for pair in ((1, 2), (2, 1), (1, nstate)):
                        h, d = analytic_nac(mol, pair=pair)
                        i, j = pair[0] - 1, pair[1] - 1
                        np.testing.assert_allclose(d[i, j], full_d[i, j], atol=1e-9)
                        np.testing.assert_allclose(h[i, j], full_h[i, j], atol=1e-9)
                        mask = np.ones((nstate, nstate), dtype=bool)
                        mask[i, j] = mask[j, i] = False
                        self.assertTrue(np.all(d[mask] == 0.0))
                        self.assertTrue(np.all(h[mask] == 0.0))
                        np.testing.assert_array_equal(mol.data["OQP::td_bvec_mo"], original)
                        self.assertEqual(mol.data._data.control.int2e_cutoff, cutoff)
                    # A Python failure after the native call still restores state.
                    with patch.dict(os.environ, {"NAC_ANALYTIC_DEBUG": str(path)}):
                        with patch("numpy.savez", side_effect=OSError("test export failure")):
                            with self.assertRaisesRegex(OSError, "test export failure"):
                                analytic_nac(mol, pair=(1, 2))
                    np.testing.assert_array_equal(mol.data["OQP::td_bvec_mo"], original)
                    self.assertEqual(mol.data._data.control.int2e_cutoff, cutoff)
                    _, repeated = analytic_nac(mol)
                    np.testing.assert_allclose(repeated, full_d, atol=1e-9)

    def test_overlap_diagnostics_survive_tag_reallocation(self):
        with tempfile.TemporaryDirectory(prefix="nac_overlap_") as tmp:
            path = Path(tmp) / "h2o.inp"
            path.write_text(INPUT.format(geometry=GEOMETRY, basis="cc-pvdz", nstate=3))
            runner = Runner(input_file=str(path), log=str(path.with_suffix(".log")))
            runner.run()
            mol = runner.mol
            original = np.array(mol.data["OQP::td_bvec_mo"], copy=True)
            oqp.mrsf_nac_metric_data(mol)
            oqp.mrsf_nac_overlap(mol)
            expected = np.array(mol.data["OQP::nac_overlap"], copy=True)
            self.assertTrue(np.all(np.isfinite(expected)))
            self.assertGreater(np.max(np.abs(expected)), 1e-10)
            with patch.dict(os.environ, {"NAC_DUMP_DS": "1"}):
                for _ in range(3):
                    oqp.mrsf_nac_overlap(mol)
                    np.testing.assert_allclose(mol.data["OQP::nac_overlap"], expected,
                                               atol=1e-12)
            np.testing.assert_array_equal(mol.data["OQP::td_bvec_mo"], original)

    def test_polarization_restores_amplitudes_and_closes_its_log(self):
        import re

        with tempfile.TemporaryDirectory(prefix="nac_polarization_") as tmp:
            path = Path(tmp) / "h2o.inp"
            path.write_text(INPUT.format(geometry=GEOMETRY, basis="6-31g", nstate=3))
            log = path.with_suffix(".log")
            runner = Runner(input_file=str(path), log=str(log))
            runner.run()
            mol = runner.mol
            original = np.array(mol.data["OQP::td_bvec_mo"], copy=True)
            target = mol.data._data.tddft.target_state
            original_directory = Path.cwd()
            try:
                os.chdir(tmp)
                with patch.dict(os.environ, {"OQP_NAC_SELFTEST": "1"}):
                    oqp.mrsf_nac_polarize(mol, 1, 2)
            finally:
                os.chdir(original_directory)
            np.testing.assert_array_equal(mol.data["OQP::td_bvec_mo"], original)
            self.assertEqual(mol.data._data.tddft.target_state, target)
            text = log.read_text()
            differences = re.findall(r"max \|de_nac\(I=J\) - de_prod\|\s*=\s*([\d.E+-]+)", text)
            self.assertEqual(len(differences), 4, text[-3000:])
            self.assertLess(max(abs(float(x)) for x in differences), 1e-10)
            self.assertIn("NAC polarization", text)
            self.assertNotIn("Z-Vector breakdown", text)
            self.assertTrue(np.all(np.isfinite(mol.data["OQP::nac_amp_polar"])))
            self.assertFalse((Path(tmp) / "fort.6").exists())
            # On Linux, verify that the final native return released the log.
            fd_directory = Path("/proc/self/fd")
            if fd_directory.exists():
                targets = []
                for fd in fd_directory.iterdir():
                    try:
                        targets.append(fd.resolve(strict=True))
                    except FileNotFoundError:
                        pass
                self.assertNotIn(log, targets)
