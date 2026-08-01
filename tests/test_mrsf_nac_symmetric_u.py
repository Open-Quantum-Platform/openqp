"""Static regression for the resident overlap-fixed MRSF NAC response."""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


class MRSFNACSymmetricUTests(unittest.TestCase):
    def test_diagonal_same_space_and_cross_space_terms_are_fortran_resident(self):
        source = (
            ROOT / "source" / "modules" / "mrsf_nac_interchange.F90"
        ).read_text()
        body = source.split(
            "subroutine mrsf_nac_rohf_pair_overlap(infos, metric_only)", 1
        )[1].split("end subroutine mrsf_nac_rohf_pair_overlap", 1)[0]
        self.assertIn("overlap_weight(p,p) = -0.5_dp*xmat(p,p)", body)
        self.assertIn(
            "coefficient = -0.5_dp*(xmat(p,q)+xmat(q,p))", body
        )
        self.assertIn("coefficient = -xmat(lo,hi)", body)
        self.assertIn("value = sum(overlap_weight_ao*dsfull", body)
        self.assertIn("gsk = sum(gamma_ao*dsket", body)
        self.assertNotIn("call ao_to_mo", body)

        production = (
            ROOT / "pyoqp" / "oqp" / "library" / "nac_analytic.py"
        ).read_text()
        driver = (
            ROOT / "source" / "modules" / "mrsf_nac_driver.F90"
        ).read_text()
        self.assertIn("oqp.mrsf_nac_lagrangian(mol)", production)
        self.assertNotIn("oqp.mrsf_nac_rohf_pair_overlap(mol)", production)
        self.assertIn(
            "call mrsf_nac_rohf_pair_overlap(infos)", driver
        )
        self.assertNotIn("def _symmetric_u_contraction", production)

    def test_fortran_wpair_keeps_diagonal_generator_derivative(self):
        source = (
            ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"
        ).read_text()
        body = source.split("subroutine mrsf_nac_wpair_batch_impl", 1)[1].split(
            "end subroutine mrsf_nac_wpair_batch_impl", 1
        )[0]
        self.assertNotIn("mt(k,k) = 0.0_dp", body)


if __name__ == "__main__":
    unittest.main()
