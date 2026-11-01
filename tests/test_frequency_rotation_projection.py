"""Rigid-rotation projection in the harmonic frequency analysis.

normal_mode() removes translations and rotations from the mass-weighted
Hessian before diagonalising it. The second rotational vector was built as

    - m^1/2 * (P_x * X_z + P_z * X_x)

where P are coordinates in the principal-axis frame and X the principal axes.
A rotation about the principal y axis is P_z * X_x - P_x * X_z, so the second
term had the wrong sign. The resulting vector is not a rigid motion: projecting
it out removed a genuine vibration and left a rotation in the vibrational space.
The term vanishes when every atom has P_z = 0 or P_x = 0, which is why planar
molecules lying in a principal plane and linear molecules were unaffected, while
every nonplanar molecule had wrong frequencies, zero-point energies and
thermochemistry. NH3 PBE0/6-31G* gave 1132, 1180, 1732, 2920, 3487, 3630 cm^-1
instead of 1132, 1732, 1732, 3487, 3629, 3629 cm^-1.

The tests use two independent references:

* a central-force (pair-spring) Hessian at its equilibrium geometry, which is
  exactly invariant to translation and rotation, so the vibrational
  frequencies are the nonzero eigenvalues of the full mass-weighted Hessian and
  carry the exact point-group degeneracies;
* an OpenQP NH3 PBE0/6-31G* analytic Hessian compared with PySCF 2.14
  frequencies at the same geometry (tests/data/nh3_pbe0_631gs_hessian.json).
"""

import importlib.util
import json
import sys
import types
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FREQ_TO_INV_CM = 5140.489195376594


def load_frequency_module():
    stub_names = ("oqp", "oqp.utils", "oqp.utils.constants")
    saved_modules = {name: sys.modules.get(name) for name in stub_names}
    try:
        oqp = types.ModuleType("oqp")
        oqp.__path__ = []
        sys.modules["oqp"] = oqp

        utils = types.ModuleType("oqp.utils")
        utils.__path__ = []
        sys.modules["oqp.utils"] = utils

        constants = types.ModuleType("oqp.utils.constants")
        constants.SPEED_OF_LIGHT = 2.99792458e10
        constants.ATMOS = 101.325
        constants.BOHR = 0.52917721090299996e-10
        constants.FREQ_TO_INV_CM = FREQ_TO_INV_CM
        constants.AMU_to_KG = 1.66053886E-27
        constants.J_TO_AU = 1 / (4.184 * 627.509541 * 1000.0)
        constants.GAS_CONSTANT = 8.3144621
        constants.PLANCK_CONSTANT = 6.62606957e-34
        constants.BOLTZMANN_CONSTANT = 1.3806488e-23
        constants.AVOGADRO_CONSTANT = 6.0221415e23
        sys.modules["oqp.utils.constants"] = constants

        spec = importlib.util.spec_from_file_location(
            "frequency_rotation_under_test",
            ROOT / "pyoqp/oqp/library/frequency.py",
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules["frequency_rotation_under_test"] = module
        spec.loader.exec_module(module)
        return module
    finally:
        for name, module in saved_modules.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


def random_rotation(seed):
    rng = np.random.default_rng(seed)
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def pair_spring_hessian(xyz, force_constants):
    """Hessian of sum_ab k_ab/2 (|r_a - r_b| - |r_a^0 - r_b^0|)^2 at r = r^0."""
    natom = len(xyz)
    hess = np.zeros((3 * natom, 3 * natom))
    for a in range(natom):
        for b in range(a + 1, natom):
            u = xyz[a] - xyz[b]
            u /= np.linalg.norm(u)
            block = force_constants(a, b) * np.outer(u, u)
            sa, sb = slice(3 * a, 3 * a + 3), slice(3 * b, 3 * b + 3)
            hess[sa, sa] += block
            hess[sb, sb] += block
            hess[sa, sb] -= block
            hess[sb, sa] -= block
    return hess


def nonzero_mass_weighted_frequencies(hess, mass, nexternal):
    mr = np.repeat(mass, 3)
    eig = np.linalg.eigvalsh(hess / np.sqrt(np.outer(mr, mr)))
    eig = np.sort(eig)[nexternal:]
    return np.sign(eig) * np.sqrt(np.abs(eig)) * FREQ_TO_INV_CM


def rigid_rotations(xyz, mass):
    """Mass-weighted infinitesimal rotations e_k x (r - R_com), k = x, y, z."""
    com = np.average(xyz, axis=0, weights=mass)
    rot = np.zeros((3 * len(mass), 3))
    for k in range(3):
        axis = np.zeros(3)
        axis[k] = 1.0
        rot[:, k] = (np.cross(axis, xyz - com) * np.sqrt(mass)[:, None]).ravel()
    return rot


# Nonplanar molecules, oriented away from the lab axes so that no principal
# coordinate of any atom vanishes. Coordinates in bohr, masses in amu.
_r, _h = 1.770792, -0.745096  # H ring radius and N-H height, bohr
NH3 = (
    np.array([[0.0, 0.0, 0.0]]
             + [[_r * np.cos(phi), _r * np.sin(phi), _h]
                for phi in (np.pi / 2, np.pi / 2 + 2 * np.pi / 3,
                            np.pi / 2 + 4 * np.pi / 3)]),
    np.array([14.00307, 1.007825, 1.007825, 1.007825]),
    [2, 2],  # two E pairs
)
_t = 2.0598 / np.sqrt(3.0)
CH4 = (
    np.array([[0.0, 0.0, 0.0], [_t, _t, _t], [-_t, -_t, _t],
              [-_t, _t, -_t], [_t, -_t, -_t]]),
    np.array([12.0, 1.007825, 1.007825, 1.007825, 1.007825]),
    [2, 3, 3],  # E and two T2 multiplets
)


def multiplets(freqs, tol):
    """Sizes of the runs of frequencies separated by less than tol."""
    sizes, run = [], 1
    for lower, upper in zip(freqs[:-1], freqs[1:]):
        if upper - lower < tol:
            run += 1
        else:
            sizes.append(run)
            run = 1
    sizes.append(run)
    return sorted(size for size in sizes if size > 1)


class FrequencyRotationProjectionTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.frequency = load_frequency_module()

    def _check_invariant_model(self, xyz0, mass, expected_multiplets, seed):
        xyz = xyz0 @ random_rotation(seed).T + np.array([0.3, -0.7, 1.1])
        # Distinct heavy-light and light-light springs keep every symmetry
        # species at a different frequency.
        heavy = int(np.argmax(mass))
        hess = pair_spring_hessian(
            xyz, lambda a, b: 0.45 if heavy in (a, b) else 0.06)

        freqs, modes, _ = self.frequency.normal_mode(xyz.ravel(), mass, hess)

        reference = nonzero_mass_weighted_frequencies(hess, mass, 6)
        np.testing.assert_allclose(freqs, reference, rtol=1e-9, atol=1e-6)
        self.assertEqual(multiplets(freqs, 1e-6), expected_multiplets)

        # Every vibrational mode is orthogonal, in the mass-weighted metric, to
        # the three rigid rotations built independently with cross products.
        q = modes * np.sqrt(np.repeat(mass, 3))[None, :]
        q /= np.linalg.norm(q, axis=1)[:, None]
        rot = rigid_rotations(xyz, mass)
        rot /= np.linalg.norm(rot, axis=0)[None, :]
        self.assertLess(np.abs(q @ rot).max(), 1e-10)

    def test_nh3_invariant_model_hessian(self):
        for seed in (1, 2, 3):
            with self.subTest(seed=seed):
                self._check_invariant_model(*NH3, seed)

    def test_ch4_invariant_model_hessian(self):
        for seed in (4, 5, 6):
            with self.subTest(seed=seed):
                self._check_invariant_model(*CH4, seed)

    def test_nh3_pbe0_hessian_matches_pyscf(self):
        data = json.loads(
            (ROOT / "tests/data/nh3_pbe0_631gs_hessian.json").read_text())
        coord = np.array(data["coord_bohr"], dtype=float).ravel()
        mass = np.array(data["mass_amu"], dtype=float)
        hess = np.array(data["hessian_hartree_per_bohr2"], dtype=float)
        hess = hess.reshape(coord.size, coord.size)
        reference = np.array(data["reference"]["frequencies_cm1"])

        freqs, _, _ = self.frequency.normal_mode(coord, mass, hess)

        # The two Hessians differ by DFT grid and SCF convergence; the largest
        # measured deviation is 1.5 cm^-1. The old projection was off by
        # 48 cm^-1 (1180 vs 1732) and 709 cm^-1 (2920 vs 3629).
        np.testing.assert_allclose(freqs, reference, atol=2.5)


if __name__ == "__main__":
    unittest.main()
