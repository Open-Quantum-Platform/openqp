"""Static contract for the state-pair ROHF NAC Z-vector API."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class MRSFNACZVectorAPITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = (
            ROOT / "source" / "modules" / "mrsf_nac_interchange.F90"
        ).read_text()
        cls.header = (ROOT / "include" / "oqp.h").read_text()
        cls.production = (
            ROOT / "pyoqp" / "oqp" / "library" / "nac_analytic.py"
        ).read_text()
        cls.driver = (
            ROOT / "source" / "modules" / "mrsf_nac_driver.F90"
        ).read_text()

    def test_state_pair_zvector_is_the_public_canonical_entry(self):
        self.assertIn("public :: mrsf_nac_rohf_zvector", self.source)
        self.assertIn("public :: mrsf_nac_rohf_zvector_batch", self.source)
        self.assertIn(
            'bind(C, name="mrsf_nac_rohf_zvector")', self.source
        )
        self.assertIn(
            "void mrsf_nac_rohf_zvector(struct oqp_handle_t *inf);",
            self.header,
        )
        self.assertIn("call mrsf_nac_rohf_zvector_batch(", self.driver)
        self.assertIn("oqp.mrsf_nac_lagrangian(mol)", self.production)
        self.assertNotIn("oqp.mrsf_nac_rohf_solve(mol)", self.production)

    def test_public_compatibility_path_still_solves_one_adjoint_rhs(self):
        one_body = self.source.split(
            "subroutine mrsf_nac_rohf_zvector(infos)", 1
        )[1].split("end subroutine mrsf_nac_rohf_zvector", 1)[0]
        batch_body = self.source.split(
            "subroutine mrsf_nac_rohf_zvector_batch(infos", 1
        )[1].split("end subroutine mrsf_nac_rohf_zvector_batch", 1)[0]
        self.assertIn("allocate(rhs(ltot,1), solution(ltot,1))", one_body)
        self.assertIn(
            "call mrsf_nac_rohf_zvector_batch(infos, rhs, solution)", one_body
        )
        self.assertIn("call cphf_solve_rohf(infos, nrhs,", batch_body)
        self.assertIn("minres_solver=.true.", batch_body)
        self.assertIn("3n forward", self.source.lower())

    def test_production_batches_antisymmetric_unordered_pairs(self):
        self.assertIn("npair = nstate*(nstate - 1)/2", self.driver)
        self.assertIn(
            "gamma_pair = pair_sign*gamma_column(:,istate)",
            self.driver,
        )
        self.assertIn(
            "rhs_batch(:,ipair) = rhs_batch(:,ipair) + rhs_in",
            self.driver,
        )
        self.assertIn("metric_only=.true.", self.driver)
        self.assertIn("integer, parameter :: z_batch_width = 3", self.driver)
        self.assertIn("do z_first = 1, npair, z_batch_width", self.driver)
        self.assertIn("call mrsf_nac_rohf_zvector_batch(", self.driver)
        self.assertIn("rhs_batch(:,z_first:z_last)", self.driver)
        self.assertIn("solution_batch(:,z_first:z_last)", self.driver)

    def test_pair_adjoint_requests_an_actual_1e_minus_10_residual_norm(self):
        body = self.source.split(
            "subroutine mrsf_nac_rohf_zvector_batch(infos", 1
        )[1].split("end subroutine mrsf_nac_rohf_zvector_batch", 1)[0]
        self.assertIn("tol=1.0e-20_dp", body)
        self.assertIn("||H z - rhs||_2 <= 1e-10", body)
        self.assertIn("do irhs = 1, nrhs", body)
        self.assertIn(".not. converged(irhs)", body)

    def test_old_solve_name_is_a_logic_free_abi_alias(self):
        self.assertIn(
            'bind(C, name="mrsf_nac_rohf_solve")', self.source
        )
        self.assertIn(
            "void mrsf_nac_rohf_solve(struct oqp_handle_t *inf);",
            self.header,
        )
        alias = self.source.split(
            "subroutine mrsf_nac_rohf_solve(infos)", 1
        )[1].split("end subroutine mrsf_nac_rohf_solve", 1)[0]
        self.assertIn("call mrsf_nac_rohf_zvector(infos)", alias)
        self.assertNotIn("cphf_solve_rohf", alias)

    def test_production_rejects_loose_scf_and_mrsf_stationarity(self):
        self.assertIn("mol.config['scf']['conv']", self.production)
        self.assertIn("mol.config['tdhf']['conv']", self.production)
        self.assertGreaterEqual(self.production.count("> 1.0e-8"), 2)

    def test_expensive_pair_kernels_are_pair_specific_and_fortran_resident(self):
        gradient = (
            ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"
        ).read_text()
        energy = (
            ROOT / "source" / "modules" / "tdhf_mrsf_energy.F90"
        ).read_text()
        self.assertIn('bind(C, name="mrsf_nac_amp_pair")', gradient)
        self.assertIn(
            "void mrsf_nac_amp_pair(struct oqp_handle_t *inf, "
            "int32_t istate, int32_t jstate);",
            self.header,
        )
        self.assertIn("call mrsf_nac_amp(infos, istate, jstate)", self.driver)
        self.assertNotIn("oqp.mrsf_nac_amp(mol)", self.production)
        self.assertIn('"OQP::nac_mt_response"', energy)
        self.assertIn("call mrsf_nac_rohf_pair_overlap(infos)", self.driver)
        self.assertNotIn("for I in range(nstate)", self.production)


if __name__ == "__main__":
    unittest.main()
