"""Compare shared and private response storage with spherical d and f shells.

The input route selects shared storage first. Repeated response calculations on
the same SCF reference then select private and shared storage. Check the output
orbital coefficients with independent signs before comparing states. Copy each native
array before replacing its owner, and compare amplitudes and transition dipoles
with the same state signs. Run at OMP_NUM_THREADS=1, 2 and 4.
"""
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


class SharedSphericalBasisTests(unittest.TestCase):
    def test_spherical_d_and_f_shells_match_private_response(self):
        import oqp
        import numpy as np
        from oqp.pyoqp import Runner

        workers = int(os.environ.get("OMP_NUM_THREADS", "1"))
        self.assertGreater(workers, 0)
        examples = Path(__file__).resolve().parents[1] / "examples" / "MRSF-TDDFT"

        def snapshot(mol):
            result = {key: np.array(mol.data["OQP::" + key], copy=True)
                      for key in ("td_energies", "td_bvec_mo", "td_trans_dipole",
                                  "VEC_MO_A", "VEC_MO_B")}
            for value in result.values():
                self.assertGreater(value.size, 0)
                self.assertTrue(np.isfinite(value).all())
            return result

        def compare(private, shared, method):
            np.testing.assert_allclose(shared["td_energies"], private["td_energies"],
                                       atol=1e-8, rtol=0)
            energies = private["td_energies"].ravel()
            self.assertEqual(energies.size, 3)
            self.assertGreater(np.min(np.abs(np.diff(energies))), 1e-7,
                               "this case requires distinct state energies")
            a = private["td_bvec_mo"].reshape(3, -1)
            b = shared["td_bvec_mo"].reshape(3, -1)
            orbital_signs = []
            for key in ("VEC_MO_A", "VEC_MO_B"):
                p, s = private[key], shared[key]
                product = np.einsum("ij,ij->i", p, s)
                self.assertTrue(np.all(product != 0))
                phase = np.where(product < 0, -1.0, 1.0)
                np.testing.assert_allclose(s * phase[:, None], p, atol=1e-6, rtol=0)
                orbital_signs.append(phase)
            if method == "mrsf":
                np.testing.assert_array_equal(*orbital_signs)
            # Ten electrons, triplet reference: six alpha and four beta.
            phase = (orbital_signs[0][:6, None] * orbital_signs[1][None, 4:]).ravel(order="F")
            self.assertEqual(b.shape[1], phase.size)
            b = b * phase
            overlap = np.einsum("ij,ij->i", a, b)
            self.assertTrue(np.all(np.abs(overlap) > 0))
            signs = np.where(overlap < 0, -1.0, 1.0)
            np.testing.assert_allclose(b * signs[:, None], a, atol=1e-5, rtol=0)
            d = private["td_trans_dipole"].ravel().reshape((3, 3, 3), order="F")
            e = shared["td_trans_dipole"].ravel().reshape((3, 3, 3), order="F")
            np.testing.assert_allclose(e * signs[None, :, None] * signs[None, None, :],
                                       d, atol=1e-6, rtol=0)

        for method in ("mrsf", "umrsf"):
            native_name = "tdhf_" + method + "_energy"
            native = getattr(oqp, native_name)
            for basis, angular in (("CC-PVDZ", 2), ("CC-PVTZ", 3)):
                with self.subTest(method=method, basis=basis, workers=workers):
                    with tempfile.TemporaryDirectory(prefix="shared_spherical_") as tmp:
                        original = examples / f"H2O_BHHLYP-{method.upper()}_SHARED_F3_{basis}.inp"
                        inp = Path(tmp) / "water.inp"
                        inp.write_text(original.read_text().replace(
                            "omp_threads = 2", f"omp_threads = {workers}"))
                        entered = []

                        def selected(mol):
                            self.assertTrue(mol.data._data.control.mrsf_shared_f3)
                            entered.append(True)
                            result = native(mol)
                            self.assertTrue(mol.mol_energy.Davidson_converged)
                            return result

                        runner = Runner(input_file=str(inp), log=str(inp.with_suffix(".log")),
                                        silent=1, usempi=False)
                        with patch.object(oqp, native_name, selected):
                            runner.run()
                        self.assertEqual(len(entered), 1, "shared energy route was not called")
                        mol = runner.mol
                        shell = mol.data.get_basis()
                        self.assertEqual(max(shell["angs"]), angular)
                        self.assertTrue(any(s and l == angular for s, l in
                                            zip(shell["spherical"], shell["angs"])))
                        first_shared = snapshot(mol)
                        mol.data.set_tdhf_shared_f3(False)
                        self.assertFalse(mol.data._data.control.mrsf_shared_f3)
                        native(mol)
                        self.assertTrue(mol.mol_energy.Davidson_converged)
                        private = snapshot(mol)
                        compare(private, first_shared, method)
                        mol.data.set_tdhf_shared_f3(True)
                        self.assertTrue(mol.data._data.control.mrsf_shared_f3)
                        native(mol)
                        self.assertTrue(mol.mol_energy.Davidson_converged)
                        compare(private, snapshot(mol), method)


if __name__ == "__main__":
    unittest.main()
