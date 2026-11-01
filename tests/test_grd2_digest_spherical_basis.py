"""The analytic MRSF gradient must respect molecular symmetry in a SPHERICAL basis.

Why this test exists
--------------------
grd2 drives the two-particle-density digests with CARTESIAN shell extents.
A digest that indexes spherical densities with those extents reads the wrong
elements. That is invisible in a Cartesian basis, where the two index spaces
coincide:

  * d functions (5 vs 6): neighbouring AO elements, wrong values and broken
    molecular symmetry.
  * f functions (7 vs 10): writes past the end of the array, heap corruption,
    SIGABRT.

The detector here is symmetry, so it needs no external reference and no
finite-difference run. The H2O geometry below has a C2 axis along z, and the
two hydrogens are exchanged by it. The energy gradient of any state must
reproduce itself under "rotate every vector by C2z and swap the two hydrogens".
A misindexed digest scrambles AO orientation and breaks that immediately.

Both AO conventions are checked: 6-31G* is Cartesian and is the control that
passed even when the spherical path was broken; cc-pVDZ is spherical.

This is the AGENTS.md rule-6 detector for every grd2 two-electron gradient
digest, so a change to any file that defines one must run it.
GUARDED_NATIVE_SOURCES names them; CI's pull-request test selection
(.github/scripts/select_pr_tests.py) selects this file by those names, and
tests/test_select_pr_tests.py fails if a digest file is missing here.
"""
import os
import tempfile
import unittest
from pathlib import Path

GUARDED_NATIVE_SOURCES = (
    "grd2",
    "caspt2_gradient",
    "casscf_ao_gradient",
    "casscf_gradient",
    "fock_deriv",
    "hf_gradient",
    "mp2_gradient",
    "nevpt2_gradient",
    "tdhf_gradient",
    "tdhf_hessian_z_rhs",
    "tdhf_mrsf_gradient",
    "tdhf_sf_gradient",
)

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
runtype=grad
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

[properties]
grad={state}
"""


NSTATE = 3


def _c2z_then_swap_hydrogens(g):
    """Apply C2z to every atom's vector and exchange the two hydrogens."""
    rotated = g * np.array([-1.0, -1.0, 1.0])
    return rotated[[0, 2, 1]]


@unittest.skipUnless(HAVE_OQP, "compiled oqp package not available")
class GradientDigestSphericalBasisTests(unittest.TestCase):
    def _gradient(self, basis, nstate=NSTATE, state=1):
        with tempfile.TemporaryDirectory(prefix="digest_spherical_") as tmp:
            path = os.path.join(tmp, "h2o.inp")
            Path(path).write_text(INPUT.format(geometry=GEOMETRY, basis=basis,
                                               nstate=nstate, state=state))
            runner = Runner(input_file=path, log=path.replace(".inp", ".log"))
            runner.run()
            return np.asarray(runner.mol.grads, dtype=float)[state].reshape(-1, 3).copy()

    def _assert_c2v(self, basis, nstate=NSTATE, state=1):
        g = self._gradient(basis, nstate, state)
        self.assertTrue(np.all(np.isfinite(g)),
                        f"{basis}: gradient is not finite:\n{g}")
        scale = np.abs(g).max()
        self.assertGreater(scale, 1e-6, f"{basis}: gradient is identically zero")
        image = _c2z_then_swap_hydrogens(g)
        residual = np.abs(image - g).max() / scale
        self.assertLess(
            residual, 1e-6,
            f"{basis}: the analytic gradient breaks C2z "
            f"(relative residual {residual:.3e}).\n"
            f"gradient:\n{g}\nC2z image with the hydrogens swapped:\n{image}")

    def test_cartesian_basis_respects_the_c2_axis(self):
        # Control: this passed even while the spherical path was broken.
        self._assert_c2v("6-31g*")

    def test_spherical_basis_respects_the_c2_axis(self):
        self._assert_c2v("cc-pvdz")

    def test_higher_root_in_a_spherical_basis(self):
        self._assert_c2v("cc-pvdz", nstate=4, state=2)


if __name__ == "__main__":
    unittest.main()
