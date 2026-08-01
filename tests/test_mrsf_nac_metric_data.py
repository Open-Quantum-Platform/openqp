"""Static contract for resident exact-TLF MRSF NAC metric data."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class MRSFNACMetricDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = (
            ROOT / "source" / "modules" / "mrsf_nac_metric_data.F90"
        ).read_text()
        cls.header = (ROOT / "include" / "oqp.h").read_text()

    def test_resident_fortran_entry_is_exported(self):
        self.assertIn("public :: mrsf_nac_metric_data", self.source)
        self.assertIn(
            'bind(C, name="mrsf_nac_metric_data")', self.source
        )
        self.assertIn(
            "void mrsf_nac_metric_data(struct oqp_handle_t *inf);",
            self.header,
        )

    def test_uses_resident_amplitudes_and_production_unfolding(self):
        self.assertIn("OQP_td_bvec_mo", self.source)
        self.assertIn("call mrsfxvec(", self.source)
        self.assertIn("noca - nocb /= 2", self.source)

    def test_differentiates_the_exact_state_overlap_in_closed_form(self):
        self.assertIn(
            "subroutine build_raw_overlap_sensitivities", self.source
        )
        self.assertIn("raw_sij(pi,ri,oi,ni)", self.source)
        self.assertIn("raw_sab(qi-nocb,si-nocb,oi,ni)", self.source)
        self.assertIn("raw_sia(pi,si-nocb,oi,ni)", self.source)
        self.assertIn("! Block 7:", self.source)
        for forbidden in (
            "call compute_states_overlap", "generator_step", "set_plane_rotation",
            "sin(angle)", "cos(angle)", "kernel_ij - kernel_ji",
        ):
            self.assertNotIn(forbidden, self.source)

    def test_uses_direct_sparse_determinant_cofactors(self):
        self.assertIn("function identity_subdet", self.source)
        self.assertIn("subroutine accumulate_identity_cofactor", self.source)
        self.assertIn(
            "cofactor = identity_subdet(rows, cols, unmatched_row, "
            "unmatched_col)",
            self.source,
        )
        self.assertIn("missing_rows > 1", self.source)
        self.assertNotIn("determinant_direction", self.source)

    def test_column_normalization_is_differentiated_exactly(self):
        self.assertIn(
            "*raw_overlap(kstate,jstate)/(column_norm(jstate)**3)",
            self.source,
        )
        self.assertIn("+ 1.0_dp/column_norm(jstate)", self.source)
        self.assertIn("raw_sij(:,:,kstate,jstate)", self.source)

    def test_gamma_slots_are_orbital_antisymmetric_but_state_ordered(self):
        self.assertIn(
            "half_derivative = 0.5_dp*(mo_gradient(p,q)-"
            "mo_gradient(q,p))",
            self.source,
        )
        self.assertIn(
            "gamma_tlf(p + (q-1)*nbf, istate, jstate) = half_derivative",
            self.source,
        )
        self.assertIn(
            "gamma_tlf(q + (p-1)*nbf, istate, jstate) = -half_derivative",
            self.source,
        )
        self.assertIn("do jstate = 1, nstate", self.source)
        self.assertIn("do istate = 1, nstate", self.source)
        self.assertNotIn("orbital_space", self.source)

    def test_four_index_work_is_not_nested_in_an_orbital_sweep(self):
        one_pass = self.source.index(
            "call build_raw_overlap_sensitivities"
        )
        state_loop = self.source.index("do jstate = 1, nstate", one_pass)
        self.assertLess(one_pass, state_loop)
        self.assertNotIn("build_identity_minor_direction", self.source)

    def test_exports_the_consumer_shape_without_nuclear_displacements(self):
        self.assertIn('tag_gamma = "OQP::nac_gamma_tlf"', self.source)
        self.assertIn("(/ nbf*nbf, nstate, nstate /)", self.source)
        self.assertIn("No displaced geometry", self.source)
        forbidden = (
            "atoms%xyz", "get_structures_ao_overlap", "der_overlap_matrix",
            "generator_step", "set_plane_rotation",
        )
        for token in forbidden:
            self.assertNotIn(token, self.source)


if __name__ == "__main__":
    unittest.main()
