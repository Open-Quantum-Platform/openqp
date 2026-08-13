"""Regression tests for Hückel guesses with an empty beta channel."""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]

INPUT_TEMPLATE = """[input]
system=
 H  0.0  0.0  0.0
charge=0
runtype=energy
basis=6-31g
method=hf

[guess]
type={guess_type}
save_mol=false

[scf]
type={scf_type}
multiplicity=2
maxit=50
conv=1.0e-9
stability=false
save_molden=false
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
class NativeGuessZeroBetaTest(unittest.TestCase):
    def _run_case(self, guess_type, scf_type):
        import oqp
        from oqp.library.single_point import SinglePoint
        from oqp.pyoqp import Runner

        with tempfile.TemporaryDirectory(prefix=f"oqp_huckel_zero_beta_{scf_type}_") as tmp:
            root = Path(tmp)
            input_file = root / "h.inp"
            input_file.write_text(
                INPUT_TEMPLATE.format(guess_type=guess_type, scf_type=scf_type)
            )
            runner = Runner(
                project=f"h_{scf_type}",
                input_file=str(input_file),
                log=str(root / "h.log"),
                silent=1,
                usempi=False,
            )
            runner.mol.data["OQP::log_filename"] = str(root / "h.log")
            oqp.oqp_banner(runner.mol)
            SinglePoint(runner.mol)._prep_guess()

            self.assertEqual(int(runner.mol.data["nelec_B"]), 0)

            beta_density = np.asarray(runner.mol.data["OQP::DM_B"], dtype=float)
            alpha_orbitals = np.asarray(
                runner.mol.data["OQP::VEC_MO_A"], dtype=float
            )
            beta_orbitals = np.asarray(
                runner.mol.data["OQP::VEC_MO_B"], dtype=float
            )
            alpha_energies = np.asarray(
                runner.mol.data["OQP::E_MO_A"], dtype=float
            )
            beta_energies = np.asarray(
                runner.mol.data["OQP::E_MO_B"], dtype=float
            )
            nbf = int(runner.mol.data.get_basis()["nbf"])
            alpha_orbitals = alpha_orbitals.reshape((nbf, nbf))
            beta_orbitals = beta_orbitals.reshape((nbf, nbf))

            self.assertTrue(np.all(np.isfinite(beta_orbitals)))
            self.assertEqual(np.linalg.matrix_rank(beta_orbitals), nbf)
            self.assertGreater(np.linalg.norm(beta_orbitals), 0.0)
            self.assertEqual(np.linalg.norm(beta_density), 0.0)
            np.testing.assert_allclose(
                beta_orbitals, alpha_orbitals, rtol=0.0, atol=0.0
            )
            np.testing.assert_allclose(
                beta_energies, alpha_energies, rtol=0.0, atol=0.0
            )

    def test_uhf_huckel_preserves_empty_beta_orbital_basis(self):
        self._run_case("huckel", "uhf")

    def test_rohf_huckel_preserves_empty_beta_orbital_basis(self):
        self._run_case("huckel", "rohf")

    def test_all_native_guesses_support_empty_beta_channel(self):
        for guess_type in ("hcore", "huckel", "modhuckel", "minao", "sap"):
            for scf_type in ("uhf", "rohf"):
                with self.subTest(guess_type=guess_type, scf_type=scf_type):
                    self._run_case(guess_type, scf_type)


if __name__ == "__main__":
    unittest.main()
