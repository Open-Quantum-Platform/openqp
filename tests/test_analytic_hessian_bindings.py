import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(relpath):
    return (ROOT / relpath).read_text()


class TestAnalyticHessianBindings(unittest.TestCase):
    def test_c_header_declares_planned_hessian_entry_points(self):
        header = read("include/oqp.h")

        for symbol in [
            "hf_hessian",
            "tdhf_hessian",
            "tdhf_sf_hessian",
        ]:
            self.assertIn(f"void {symbol}(struct oqp_handle_t *inf);", header)

    def test_python_dispatch_mentions_native_hessian_entry_points(self):
        source = read("pyoqp/oqp/library/single_point.py")

        for symbol in [
            "oqp.hf_hessian",
            "oqp.tdhf_hessian",
            "oqp.tdhf_sf_hessian",
        ]:
            self.assertIn(symbol, source)

    def test_hf_native_dispatch_is_explicitly_not_a_numerical_fallback(self):
        source = read("pyoqp/oqp/library/single_point.py")

        self.assertIn("native_hess_func", source)
        self.assertIn("no numerical fallback", source.lower())
        self.assertIn("Native HF/DFT analytic Hessian kernels are not available", source)
        self.assertIn("native_hess(self.mol)", source)

    def test_public_sf_dispatch_is_separate_from_private_mrsf_path(self):
        source = read("pyoqp/oqp/library/single_point.py")

        self.assertIn("def analytical_sf_hess", source)
        self.assertIn("td_type == 'sf'", source)
        self.assertIn("return self.analytical_sf_hess()", source)

    def test_hf_hessian_fortran_scaffold_exports_c_abi_without_claiming_support(self):
        source = read("source/modules/hf_hessian.F90")

        self.assertIn("module hf_hessian_mod", source)
        self.assertIn('bind(C, name="hf_hessian")', source)
        self.assertIn("subroutine hf_hessian_C", source)
        self.assertIn("subroutine hf_hessian", source)
        self.assertIn("Analytic HF/DFT Hessian kernel scaffold reached", source)
        self.assertIn("WITH_ABORT", source)

    def test_tdhf_hessian_fortran_scaffolds_export_c_abi_without_claiming_support(self):
        expected = {
            "tdhf_hessian.F90": [
                "module tdhf_hessian_mod",
                'bind(C, name="tdhf_hessian")',
                "subroutine tdhf_hessian_C",
                "subroutine tdhf_hessian",
                "Analytic TDDFT Hessian kernel scaffold reached",
            ],
            "tdhf_sf_hessian.F90": [
                "module tdhf_sf_hessian_mod",
                'bind(C, name="tdhf_sf_hessian")',
                "subroutine tdhf_sf_hessian_C",
                "subroutine tdhf_sf_hessian",
                "Analytic SF-TDDFT Hessian kernel scaffold reached",
            ],
        }

        for filename, needles in expected.items():
            with self.subTest(filename=filename):
                source = read(f"source/modules/{filename}")
                for needle in needles:
                    self.assertIn(needle, source)
                self.assertIn("WITH_ABORT", source)

    def test_molecule_has_single_hessian_storage_helper_with_asymmetry_metadata(self):
        source = read("pyoqp/oqp/molecule/molecule.py")

        self.assertIn("def set_hessian_result", source)
        self.assertIn("hessian_metadata", source)
        self.assertIn("max_asymmetry", source)
        self.assertIn("0.5 * (hessian + hessian.T)", source)
        self.assertIn("def get_hess(self):", source)
        self.assertNotIn("def get_hess(self):\n        \"\"\"\n        Get hessian results\n        \"\"\"\n\n        return []", source)

    def test_saved_hessian_json_uses_inertia_not_modes_for_inertia_field(self):
        source = read("pyoqp/oqp/molecule/molecule.py")

        self.assertIn("'inertia': self.inertia.tolist()", source)
        self.assertNotIn("'inertia': self.modes.tolist()", source)

    def test_read_hessian_json_restores_hessian_metadata(self):
        source = read("pyoqp/oqp/molecule/molecule.py")

        self.assertIn("self.hessian_metadata = data.get('hessian_metadata', {})", source)


if __name__ == "__main__":
    unittest.main()
