"""Analytic IR and Raman tensors of the native ground-state Hessian.

The analytic Hessian stores the nuclear derivatives of the dipole moment
(OQP::hf_dipole_derivatives) and of the static polarizability
(OQP::hf_polarizability_derivatives). They replace 6N displaced SCF + CPHF
runs, so they are checked here against exactly that finite-difference path
(Hessian._native_property_tensors_at, central differences, h = 1e-3 bohr):

  * RHF in a spherical basis (cc-pVDZ, ispher on),
  * RHF with an ECP (HBr / LANL2DZ, 28 core electrons on Br),
  * UHF and ROHF with an ECP: the CH2Br radical (non-degenerate doublet) with
    LANL2DZ, so the open-shell hf_polder_uhf / hf_polder_rohf paths are
    checked with the core-removing ECP derivative integrals in their
    right-hand sides.

In addition, closed-shell water run through the UHF and ROHF kernels must
reproduce the RHF tensors, and the H2O+ UHF and ROHF tensors must rotate with
the molecule, sum_A (n x R_A).dalpha/dR_A = [Omega_n, alpha].

Everything runs in one process, one molecule after another, so state left
behind by a previous molecule (basis size, atom count, spin) would show up as
a mismatch.

Before these tensors existed the tags are absent and the test fails. The
relative tolerances sit above the CPHF solver tolerance (sqrt(1e-9)) amplified
by the finite difference, and far below a missing or wrongly signed term.

Skipped unless the compiled OpenQP runtime is importable.
"""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

WATER = """   8   0.000000000   0.000000000   0.117300000
   1   0.000000000   0.757200000  -0.469200000
   1   0.000000000  -0.757200000  -0.469200000"""
HBR = """   35  0.000000000   0.000000000   0.000000000
    1  0.100000000   0.050000000   1.414000000"""
CH2BR = """    6  0.000  0.000  0.000
   35  0.050  0.030  1.880
    1  0.930  0.000 -0.520
    1 -0.930  0.060 -0.540"""

FD_CASES = [
    # name, geometry, charge, basis, scf type, multiplicity, extra input lines
    ("rhf_ccpvdz_spherical", WATER, 0, "cc-pvdz", "rhf", 1, "ispher=true"),
    ("rhf_hbr_lanl2dz_ecp", HBR, 0, "lanl2dz", "rhf", 1, ""),
    ("uhf_ch2br_lanl2dz_ecp", CH2BR, 0, "lanl2dz", "uhf", 2, "ispher=true"),
    ("rohf_ch2br_lanl2dz_ecp", CH2BR, 0, "lanl2dz", "rohf", 2, "ispher=true"),
]

INPUT_TMPL = """[input]
system=
{geometry}
charge={charge}
runtype=hess
basis={basis}
method=hf
{extra}
[guess]
type=huckel
[scf]
type={scftype}
multiplicity={mult}
conv=1.0e-10
maxit=200
[symmetry]
enabled=false
[hess]
type=analytical
state=0
"""


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "2")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


def _tag(mol, name, shape):
    """Fortran (3,3N) / (3,3,3N) array from the tagarray, in that index order."""
    raw = np.array(mol.data[name], dtype=float).reshape(-1)
    return raw.reshape(shape[::-1]).T


def _run(tmp, name, geometry, charge, basis, scftype, mult, extra):
    from oqp.pyoqp import Runner
    inp = Path(tmp) / f"{name}.inp"
    inp.write_text(INPUT_TMPL.format(geometry=geometry, charge=charge, basis=basis,
                                     scftype=scftype, mult=mult, extra=extra))
    runner = Runner(project=name, input_file=str(inp), log=str(Path(tmp) / f"{name}.log"),
                    silent=1, usempi=False)
    runner.run()
    mol = runner.mol
    ncoord = np.asarray(mol.get_system(), dtype=float).size
    return (mol, _tag(mol, "OQP::hf_dipole_derivatives", (3, ncoord)),
            _tag(mol, "OQP::hf_polarizability_derivatives", (3, 3, ncoord)))


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class AnalyticIrRaman(unittest.TestCase):
    def test_analytic_tensors(self):
        from oqp.library.single_point import Hessian

        h = 1.0e-3
        with tempfile.TemporaryDirectory() as tmp:
            for case in FD_CASES:
                name = case[0]
                with self.subTest(case=name):
                    mol, dmu, dalpha = _run(tmp, *case)
                    coord = np.asarray(mol.get_system(), dtype=float).reshape(-1, 3)
                    ncoord = coord.size
                    hess = Hessian(mol)
                    flat = coord.reshape(-1)
                    fd_mu = np.zeros((3, ncoord))
                    fd_alpha = np.zeros((3, 3, ncoord))
                    for i in range(ncoord):
                        step = np.zeros(ncoord)
                        step[i] = h
                        mu_p, a_p = hess._native_property_tensors_at((flat + step).reshape(coord.shape))
                        mu_m, a_m = hess._native_property_tensors_at((flat - step).reshape(coord.shape))
                        fd_mu[:, i] = (mu_p - mu_m) / (2.0 * h)
                        fd_alpha[:, :, i] = (a_p - a_m) / (2.0 * h)
                    rel_mu = np.abs(dmu - fd_mu).max() / np.abs(fd_mu).max()
                    rel_alpha = np.abs(dalpha - fd_alpha).max() / np.abs(fd_alpha).max()
                    self.assertLess(rel_mu, 5.0e-4, f"{name}: dmu/dR rel. error {rel_mu:.2e}")
                    self.assertLess(rel_alpha, 1.0e-3, f"{name}: dalpha/dR rel. error {rel_alpha:.2e}")
                    self.assertLess(np.abs(dalpha - dalpha.transpose(1, 0, 2)).max(), 1.0e-10)
                    self.assertEqual(mol.vibrational_intensity_metadata.get("backend"),
                                     "native_openqp_analytic")

            # open-shell kernels on a closed-shell molecule reproduce RHF
            _, mu_r, al_r = _run(tmp, "closed_rhf", WATER, 0, "6-31g", "rhf", 1, "")
            for scftype in ("uhf", "rohf"):
                with self.subTest(case=f"closed_{scftype}"):
                    _, mu_o, al_o = _run(tmp, f"closed_{scftype}", WATER, 0, "6-31g", scftype, 1, "")
                    self.assertLess(np.abs(mu_o - mu_r).max(), 1.0e-4)
                    self.assertLess(np.abs(al_o - al_r).max() / np.abs(al_r).max(), 1.0e-4)

            # open-shell tensors rotate with the molecule
            eps = np.zeros((3, 3, 3))
            eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1.0
            eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1.0
            for scftype in ("uhf", "rohf"):
                with self.subTest(case=f"cation_{scftype}_rotation"):
                    mol, _, dalpha = _run(tmp, f"cation_{scftype}", WATER, 1, "6-31g", scftype, 2, "")
                    coord = np.asarray(mol.get_system(), dtype=float).reshape(-1, 3)
                    _, alpha0 = Hessian(mol)._native_property_tensors_at(coord)
                    for k in range(3):
                        n = np.zeros(3)
                        n[k] = 1.0
                        omega = -np.einsum("ijk,k->ij", eps, n)
                        lhs = np.einsum("abx,x->ab", dalpha, np.cross(n, coord).reshape(-1))
                        rhs = omega @ alpha0 - alpha0 @ omega
                        self.assertLess(np.abs(lhs - rhs).max() / np.abs(alpha0).max(), 1.0e-3)


if __name__ == "__main__":
    unittest.main()
