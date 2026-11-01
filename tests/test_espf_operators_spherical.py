"""ESPF atomic-charge operators in spherical (pure d/f) bases.

The ESPF operators are formed in one pass over the fitting grid for a block
of QM atoms (int1::omp_qmmm_all).  That path has its own Cartesian-to-
spherical transformation of every shell block, and the shipped QM/MM decks
use only Pople (Cartesian) bases, so a pure-d/f indexing or normalization
error would not show up in them.

Pinned here against the per-atom omp_qmmm path of gateway main 3ff4121fa,
for cc-pVDZ (pure d) and cc-pVTZ (pure d, f): weighted sums of every packed
operator element and the ESPF charges of a non-symmetric open-shell
4-atom system, so a transposed index, a swapped atom column or a sign error
changes the numbers.  The operators must not depend on the atom block size
(ESPF_ATOM_BLOCK), which bounds their memory.

Requires the compiled OpenQP runtime (OPENQP_ROOT); skipped otherwise.
"""
import os
import re
import unittest
from pathlib import Path

import numpy as np

try:
    import oqp
    from oqp.openqp import OPENQP
    _HAVE_OQP = True
except Exception:  # pragma: no cover - uncompiled backend / missing OPENQP_ROOT
    _HAVE_OQP = False

XYZ = [("O", 0.000, 0.000, 0.110), ("H", 0.951, 0.040, -0.420),
       ("H", -0.262, 0.905, -0.390), ("N", 0.300, -1.900, 0.700)]

# Per-atom omp_qmmm path, gateway main 3ff4121fa (macOS arm64, Accelerate ILP64).
REFERENCE = {
    "cc-pvdz": {
        "nbf": 38,
        "w1": [5.5037242227924965, -3.908000375268504, 1.0904652692427823, -7.239394117052844],
        "w2": [25.9617942413099, -9.86118820856804, 0.06610688025147704, -19.87190305153546],
        "charges": [0.2511358914208983, 0.360140796377751, 0.3399942526201357, -0.9512709404187849],
    },
    "cc-pvtz": {
        "nbf": 88,
        "w1": [-6.958596333291371, -6.78348344311753, -1.107102705311045, -3.023961145596322],
        "w2": [14.553120953077181, -21.86098884641638, 14.707561955806431, 14.945419965635113],
        "charges": [0.2539559785146289, 0.3581647006387373, 0.3398848618531485, -0.9520055410065148],
    },
}


def _espf(basis, block=None):
    if block is None:
        os.environ.pop("ESPF_ATOM_BLOCK", None)
    else:
        os.environ["ESPF_ATOM_BLOCK"] = str(block)
    try:
        sys_str = "\n" + "\n".join(f"{l} {x:.6f} {y:.6f} {z:.6f}" for l, x, y, z in XYZ)
        cfg = {"input.system": sys_str, "input.charge": 0, "input.runtype": "energy",
               "input.basis": basis, "input.method": "hf", "guess.type": "huckel",
               "scf.multiplicity": 2, "scf.type": "uhf"}
        op = OPENQP(cfg, True)
        mol = op.mol
        op.sp._prep_guess()
        nat = mol.data["natom"]
        oqp.espf_op_corr(mol)
        ops = np.array(mol.data["OQP::ESPF_CORR"]).reshape(nat, -1)
        oqp.form_esp_charges(mol)
        charges = np.array(mol.data["OQP::partial_charges"]).ravel()
        return int(mol.data.get_basis()["nbf"]), ops, charges
    finally:
        os.environ.pop("ESPF_ATOM_BLOCK", None)


def _weighted(ops):
    k = np.arange(ops.shape[1])
    w1 = np.cos(0.37 * k + 0.1)
    w2 = np.sin(1.31 * k) * (1.0 + (k % 7))
    return ops @ w1, ops @ w2


@unittest.skipUnless(_HAVE_OQP and os.environ.get("OPENQP_ROOT"),
                     "compiled OpenQP runtime not available")
class TestEspfOperatorsSpherical(unittest.TestCase):
    def setUp(self):
        # The ESPF grid and weights read ESPF_* controls from the environment,
        # and the QM/MM driver sets some of them in-process (ESPF_SWSCALE for a
        # covalent boundary), so an earlier test can leave them behind.  The
        # reference is the shipped default.
        self._espf_env = {k: v for k, v in os.environ.items() if k.startswith("ESPF_")}
        for k in self._espf_env:
            del os.environ[k]

    def tearDown(self):
        for k in [k for k in os.environ if k.startswith("ESPF_")]:
            del os.environ[k]
        os.environ.update(self._espf_env)

    def _check_reference(self, basis):
        ref = REFERENCE[basis]
        nbf, ops, charges = _espf(basis)
        self.assertEqual(nbf, ref["nbf"])
        w1, w2 = _weighted(ops)
        np.testing.assert_allclose(w1, ref["w1"], rtol=0, atol=1e-9)
        np.testing.assert_allclose(w2, ref["w2"], rtol=0, atol=1e-9)
        np.testing.assert_allclose(charges, ref["charges"], rtol=0, atol=1e-8)
        return ops

    def test_cc_pvdz_matches_per_atom_operators(self):
        self._check_reference("cc-pvdz")

    def test_cc_pvtz_matches_per_atom_operators(self):
        self._check_reference("cc-pvtz")

    def test_operators_do_not_depend_on_the_atom_block(self):
        _, ops, charges = _espf("cc-pvdz")
        for block in (1, 3):
            _, ops_b, charges_b = _espf("cc-pvdz", block=block)
            self.assertEqual(np.abs(ops_b - ops).max(), 0.0, f"block={block}")
            self.assertEqual(np.abs(charges_b - charges).max(), 0.0, f"block={block}")

    def test_unusable_block_settings_keep_the_default(self):
        # empty, null values, text, zero and negative: the default block
        _, ops, charges = _espf("cc-pvdz")
        for block in ("", ",", "/", "abc", "0", "-3"):
            _, ops_b, charges_b = _espf("cc-pvdz", block=block)
            self.assertEqual(np.abs(ops_b - ops).max(), 0.0, f"block={block!r}")
            self.assertEqual(np.abs(charges_b - charges).max(), 0.0, f"block={block!r}")


QMMM_F90 = Path(__file__).resolve().parents[1] / "source" / "modules" / "qmmm.F90"


class TestEspfAtomBlockParser(unittest.TestCase):
    """ESPF_ATOM_BLOCK may hold a list-directed null value (',' or '/'), which
    a successful read leaves unassigned, or text that fails to read.  The
    requested size must then be zero (keep the memory-bounded default), never
    an undefined integer; Fortran does not short-circuit `ios == 0 .and. nreq > 0`."""

    def test_requested_size_is_defined_before_it_is_tested(self):
        src = QMMM_F90.read_text()
        m = re.search(r"integer function espf_atom_block[\s\S]*?end function espf_atom_block", src)
        self.assertIsNotNone(m)
        body = m.group(0)
        self.assertRegex(body, r"nreq = 0\s*\n\s*read\(env_b, \*, iostat=ios\) nreq\s*\n"
                               r"\s*if \(ios /= 0\) nreq = 0\s*\n\s*if \(nreq > 0\) nblk = min\(nat, nreq\)")
        self.assertNotRegex(body, r"ios == 0 \.and\. nreq")


if __name__ == "__main__":
    unittest.main()
