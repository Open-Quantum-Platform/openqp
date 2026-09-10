"""Finite-difference consistency of the analytic ESPF QM/MM gradient.

The ESPF fitting grid is atom-centred and translates rigidly with its parent
atom, so the analytic gradient of the embedding energy must include the
derivative through the moving grid points (integral, kernel and switching
terms).  Without those terms the QM gradient is neither translationally
invariant nor the derivative of the energy (residual ~1e-4 Ha/bohr, an
energy drift in NVE dynamics).

Checked here at a FIXED MM potential (the classical dphi/dR coupling force is
the driver's job and is exact by construction): displace each QM atom of a
water molecule embedded in three point charges, central-difference the
embedded SCF energy, and compare with hf_gradient + grad_esp_qmmm.  The sum of
the analytic gradient over the QM atoms must vanish (rigid translation of the
QM region and its grid leaves the fixed-potential energy unchanged).

Requires the compiled OpenQP runtime (OPENQP_ROOT); skipped otherwise.
"""
import os
import unittest

import numpy as np

try:
    import oqp
    from oqp.openqp import OPENQP
    from oqp.library.qmmm_driver import (
        unpack_lower_tri_single, unpack_lower_tri_multi, pack_lower_tri_single)
    _HAVE_OQP = True
except Exception:  # pragma: no cover - uncompiled backend / missing OPENQP_ROOT
    _HAVE_OQP = False

ANG2BOHR = 1.8897259886

QM_XYZ_ANG = np.array([[0.000, 0.000, 0.000],
                       [0.957, 0.000, 0.000],
                       [-0.239, 0.927, 0.000]])
MM_XYZ_BOHR = np.array([[2.500, 0.300, 0.100],
                        [3.100, 0.900, 0.100],
                        [2.900, -0.40, -0.30]]) * ANG2BOHR
MM_CHARGES = np.array([-0.834, 0.417, 0.417])


def _potmm(qm_xyz_bohr):
    pot = np.zeros(len(qm_xyz_bohr))
    for a, ra in enumerate(qm_xyz_bohr):
        d = MM_XYZ_BOHR - ra
        pot[a] = np.sum(MM_CHARGES / np.sqrt(np.einsum("ij,ij->i", d, d)))
    return pot


@unittest.skipUnless(_HAVE_OQP and os.environ.get("OPENQP_ROOT"),
                     "compiled OpenQP runtime not available")
class TestEspfMovingGridGradient(unittest.TestCase):
    def setUp(self):
        sys_str = "\n" + "\n".join(
            f"{l} {x:.10f} {y:.10f} {z:.10f}"
            for l, (x, y, z) in zip(["O", "H", "H"], QM_XYZ_ANG))
        cfg = {"input.system": sys_str, "input.charge": 0, "input.runtype": "grad",
               "input.basis": "6-31g", "input.method": "hf", "guess.type": "huckel",
               "scf.multiplicity": 1, "scf.type": "rhf", "scf.conv": 1e-11}
        self.op = OPENQP(cfg, True)
        self.mol, self.sp = self.op.mol, self.op.sp
        self.sp._prep_guess()
        self.nat = self.mol.data["natom"]
        self.nbf = self.mol.data.get_basis()["nbf"]
        self.r0 = np.asarray(self.mol.get_atoms2("coords")).reshape(self.nat, 3).copy()
        self.pot0 = _potmm(self.r0)

    def _energy(self, coords):
        mol, sp, nat, nbf = self.mol, self.sp, self.nat, self.nbf
        mol.data.set_scf_conv(1e-12)
        mol.data.set_scf_maxit(300)
        mol.update_system(coords.ravel())
        sp._prep_guess()
        mol.data["OQP::POTMM"] = self.pot0.copy()
        mol.data["OQP::POTQM"] = np.zeros((nat, nat))
        oqp.espf_op_corr(mol)
        espf = unpack_lower_tri_multi(mol.data["OQP::ESPF_CORR"], nbf, nat)
        hcore = unpack_lower_tri_single(mol.get_hcore(), nbf)
        hcore += np.einsum("ijk,i->jk", espf, self.pot0)
        mol.set_hcore(pack_lower_tri_single(hcore))
        sp.scf()
        return float(mol.get_scf_energy())

    def test_fixed_potential_gradient_matches_finite_difference(self):
        mol, nat = self.mol, self.nat
        self._energy(self.r0)
        oqp.form_esp_charges(mol)
        oqp.hf_gradient(mol)
        g_an = np.asarray(mol.get_grad()).reshape(nat, 3).copy()
        oqp.grad_esp_qmmm(mol)
        g_an += np.asarray(mol.data["OQP::ESPF_GRAD"]).reshape(nat, 3)

        # Rigid translation of QM atoms + their grid at fixed potential: dE = 0.
        self.assertLess(np.abs(g_an.sum(axis=0)).max(), 1.0e-9)

        h = 1.0e-3
        g_num = np.zeros((nat, 3))
        for a in range(nat):
            for c in range(3):
                rp = self.r0.copy(); rp[a, c] += h
                rm = self.r0.copy(); rm[a, c] -= h
                g_num[a, c] = (self._energy(rp) - self._energy(rm)) / (2 * h)
        self.assertLess(np.abs(g_an - g_num).max(), 2.0e-6)


if __name__ == "__main__":
    unittest.main()
