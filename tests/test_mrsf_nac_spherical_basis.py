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
nstate=3
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

    def _derivative_coupling(self, basis):
        with tempfile.TemporaryDirectory(prefix="nac_spherical_") as tmp:
            path = os.path.join(tmp, "h2o.inp")
            with open(path, "w") as handle:
                handle.write(INPUT.format(geometry=GEOMETRY, basis=basis))
            runner = Runner(input_file=path, log=path.replace(".inp", ".log"))
            runner.run()
            # OQP::nac_dcv is reserved as (3*natom, nstate, nstate) in Fortran
            # order; take the (state 0, state 1) pair and fold the coordinate
            # axis back into (natom, 3).
            dcv = np.array(runner.mol.data["OQP::nac_dcv"], copy=True)
            dcv = dcv.reshape(-1).reshape((-1, NSTATE, NSTATE), order="F")
        # Fortran state 1/2 are stored at the 0-based positions 1 and 2; fold
        # the coordinate axis back into (natom, 3).
        return dcv[:, 1, 2].reshape(-1, 3)

    def _assert_c2v(self, basis):
        d = self._derivative_coupling(basis)
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
