"""After an ROHF SCF the OQP::DM_A / OQP::DM_B records must hold the alpha and
beta densities on EVERY exit path.  The final write-back used the working
density as if it held the total density, which is true only after the ROHF
combination step; a DIIS loop stopping at its iteration limit (or the
SOSCF/TRAH paths) left the ALPHA density there, so DM_A became the spin
density alpha-beta and every post-SCF consumer (ESPF charges, Mulliken,
MRSF relaxed densities) worked on a wrong density.  Needs the compiled
runtime."""
import os
import unittest

import numpy as np

try:
    from oqp.openqp import OPENQP
    _HAVE = os.environ.get("OPENQP_ROOT") is not None
except Exception:  # pragma: no cover
    _HAVE = False

H2O = "O 0.0 0.0 0.1173; H 0.0 0.7572 -0.4692; H 0.0 -0.7572 -0.4692"


def _unpack(p, n):
    m = np.zeros((n, n)); k = 0
    for i in range(n):
        for j in range(i + 1):
            m[i, j] = m[j, i] = p[k]; k += 1
    return m


def _traces(mol):
    S = np.asarray(mol.data["OQP::SM"], dtype=float).ravel()
    n = int((np.sqrt(8 * len(S) + 1) - 1) / 2)
    Sm = _unpack(S, n)
    ta = float(np.sum(_unpack(np.asarray(mol.data["OQP::DM_A"], dtype=float).ravel(), n) * Sm))
    tb = float(np.sum(_unpack(np.asarray(mol.data["OQP::DM_B"], dtype=float).ravel(), n) * Sm))
    return ta, tb


@unittest.skipUnless(_HAVE, "compiled OpenQP runtime (OPENQP_ROOT) unavailable")
class TestROHFDensityRecords(unittest.TestCase):
    def _run(self, maxit, converger="diis"):
        cfg = {"input.runtype": "energy", "input.method": "hf", "input.basis": "6-31g",
               "input.charge": 0, "input.system": H2O, "scf.type": "rohf",
               "scf.multiplicity": 3, "scf.maxit": maxit, "scf.conv": 1e-8,
               "scf.converger_type": converger, "guess.type": "huckel"}
        op = OPENQP(cfg, True)
        op.sp._prep_guess()
        op.sp.scf()                       # plain SCF, no escalation: exercises the maxit exit
        return op.mol

    def test_alpha_beta_on_converged_and_iteration_limit_exits(self):
        # triplet H2O: 6 alpha, 4 beta electrons
        for maxit in (50, 2):
            mol = self._run(maxit)
            ta, tb = _traces(mol)
            self.assertAlmostEqual(ta, 6.0, places=6, msg=f"maxit={maxit}: Tr[DM_A S]={ta}")
            self.assertAlmostEqual(tb, 4.0, places=6, msg=f"maxit={maxit}: Tr[DM_B S]={tb}")

    def test_alpha_beta_after_trah(self):
        mol = self._run(50, converger="trah")
        ta, tb = _traces(mol)
        self.assertAlmostEqual(ta, 6.0, places=6)
        self.assertAlmostEqual(tb, 4.0, places=6)
