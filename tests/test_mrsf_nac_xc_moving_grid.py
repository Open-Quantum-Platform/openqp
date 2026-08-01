"""Static contract for the interstate XC moving-grid derivative."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class MRSFNACXCMovingGridTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.consumer = (
            ROOT / "source" / "dftlib" / "dft_gridint_tdxc_grad.F90"
        ).read_text()
        cls.esum = (
            ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"
        ).read_text()
        cls.adjoint = (
            ROOT / "source" / "modules" / "mrsf_nac_interchange.F90"
        ).read_text()
        cls.partfunc = (
            ROOT / "source" / "dftlib" / "dft_partfunc.F90"
        ).read_text()
        cls.oqpdata = (
            ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"
        ).read_text()

    def test_linear_probe_api_disables_ground_state_and_enables_weights(self):
        self.assertIn("include_ground_state", self.consumer)
        self.assertIn("include_weight_derivative", self.consumer)
        for production_source in (self.esum, self.adjoint):
            self.assertIn("include_ground_state=.false.", production_source)
            self.assertIn(
                "include_weight_derivative=.true.", production_source
            )

    def test_production_grid_response_has_no_environment_kill_switch(self):
        self.assertNotIn("NAC_ESUM_NO_GRID_DERIV", self.esum)
        esum_body = self.esum.split(
            "subroutine mrsf_nac_esum(infos", 1
        )[1].split("end subroutine mrsf_nac_esum", 1)[0]
        self.assertIn("include_weight_derivative=.true.", esum_body)

    def test_production_esum_has_no_coordinate_fd_debug_branch(self):
        esum_body = self.esum.split(
            "subroutine mrsf_nac_esum(infos", 1
        )[1].split("end subroutine mrsf_nac_esum", 1)[0]
        self.assertNotIn("NAC_ESUM_FDTEST", esum_body)
        self.assertNotIn("/tmp/nac_esum_", esum_body)

    def test_diagonal_gradient_adds_probe_only_grid_correction(self):
        self.assertIn("weight_derivative_only", self.consumer)
        self.assertIn("requested_weight_only", self.consumer)
        self.assertIn("if (.not. requested_weight_only) then", self.consumer)
        main = self.esum.split(
            "subroutine tdhf_mrsf_gradient(infos)", 1
        )[1].split("end subroutine tdhf_mrsf_gradient", 1)[0]
        self.assertIn("grid_correction", main)
        self.assertIn("include_ground_state=.false.", main)
        self.assertIn("include_weight_derivative=.true.", main)
        self.assertIn("weight_derivative_only=.true.", main)

    def test_combined_ground_probe_call_does_not_request_grid_response(self):
        main = self.esum.split(
            "subroutine tdhf_mrsf_gradient(infos)", 1
        )[1].split("grid_correction", 1)[0]
        self.assertNotIn("include_weight_derivative=.true.", main)

    def test_partition_and_owner_motion_are_both_present(self):
        body = self.consumer.split(
            "subroutine add_partition_weight_gradient", 1
        )[1].split("end subroutine add_partition_weight_gradient", 1)[0]
        self.assertIn("partfunc%deriv(mu)", body)
        self.assertIn("dlog(:,b,owner)", body)
        self.assertIn("merge(1,0,b == owner)", body)
        self.assertIn("sum(tmpGrad, dim=1)", self.consumer)

    def test_partition_derivatives_include_coordinate_chain_rules(self):
        self.assertIn("(1.0_fp+x*x)*frac*frac", self.partfunc)
        self.assertEqual(self.partfunc.count("df = -SCALEF"), 4)
        self.assertIn("partfunc%limit = 0.73_fp", self.partfunc)

    def test_python_partition_names_match_fortran_type_ids(self):
        self.assertIn("{'ssf': 0, 'erf': 1, 'becke': 2", self.oqpdata)

    def test_moving_grid_rejects_fxc_mode_and_preserves_probe_axis(self):
        self.assertIn("requested_weight_derivative .and. doFxc", self.consumer)
        self.assertNotIn(
            "requested_weight_derivative .and. nMtx /= 1", self.consumer
        )
        self.assertIn("dedft_mtx", self.consumer)
        self.assertIn("self%probe_value(ipt,imtx,mythread)", self.consumer)
        self.assertIn("self%nucgrad(:,b,imtx,mythread)", self.consumer)
        self.assertIn(
            "dat%do_weight_derivative .and. dat%do_ground_state",
            self.consumer,
        )

    def test_partition_work_is_preallocated_and_quadratic_in_atom_count(self):
        body = self.consumer.split(
            "subroutine add_partition_weight_gradient", 1
        )[1].split("end subroutine add_partition_weight_gradient", 1)[0]
        self.assertNotIn("allocate(", body)
        self.assertIn("part_dlog", self.consumer)
        self.assertIn("O(Ngrid*Natom**3)", self.consumer)

    def test_partition_thread_scratch_uses_cache_line_atom_stride(self):
        self.assertIn("nat_stride = ((nat+7)/8)*8", self.consumer)
        self.assertIn("part_dist(nat_stride,nthreads)", self.consumer)
        self.assertIn("part_cells(nat_stride,nthreads)", self.consumer)
        self.assertIn("part_dlog(3,nat,nat_stride,nthreads)", self.consumer)

    def test_batch_gradient_scratch_uses_pruned_ao_extent(self):
        self.assertIn(
            "tmpGrad(1:numAOs,1:3,1:nMtx)", self.consumer
        )
        self.assertIn(
            "self%tmpGrad_(1:numAOs*3*nMtx,myThread)", self.consumer
        )
        self.assertNotIn(
            "tmpGrad => self%tmpGrad_(:,:,:,myThread)", self.consumer
        )

    def test_probe_weight_term_is_not_ground_state_xc_energy(self):
        self.assertIn("probe_value", self.consumer)
        self.assertIn("dot_product(d_r, rhoab)", self.consumer)
        self.assertIn("dot_product(d_s, sigma)", self.consumer)
        self.assertIn("dot_product(d_t, tauab)", self.consumer)

    def test_zero_probe_subtraction_is_absent_from_production_paths(self):
        esum_body = self.esum.split("subroutine mrsf_nac_esum(infos", 1)[1]
        esum_body = esum_body.split("end subroutine mrsf_nac_esum", 1)[0]
        adjoint_body = self.adjoint.split(
            "subroutine mrsf_nac_xc_adjoint(infos)", 1
        )[1].split("end subroutine mrsf_nac_xc_adjoint", 1)[0]
        self.assertNotIn("gxc0", esum_body)
        self.assertNotIn("gxc0", adjoint_body)

    def test_adjoint_builds_the_xc_kernel_without_redundant_jk_focks(self):
        body = self.adjoint.split(
            "subroutine mrsf_nac_xc_adjoint(infos)", 1
        )[1].split("end subroutine mrsf_nac_xc_adjoint", 1)[0]
        self.assertIn("call utddft_fxc(", body)
        self.assertNotIn("call fock_jk(", body)
        self.assertNotIn("call get_response_packed(", body)
        self.assertNotIn("call pack_matrix(", body)

    def test_production_batches_pair_probes_in_one_moving_grid_call(self):
        body = self.adjoint.split(
            "subroutine mrsf_nac_xc_adjoint_batch(infos", 1
        )[1].split("end subroutine mrsf_nac_xc_adjoint_batch", 1)[0]
        self.assertEqual(body.count("call dft_initialize("), 1)
        self.assertEqual(body.count("call utddft_xc_gradient("), 1)
        self.assertIn("nMtx=nrhs", body)
        self.assertIn("dedft_mtx=gxc_vectors", body)
        self.assertIn("call utddft_fxc(", body)
        self.assertIn("call ao_to_mo_occ(dsa(:,:,cart,atom)", body)
        self.assertIn("nocc, nocc, nbf", body)
        self.assertNotIn("call fock_jk(", body)


if __name__ == "__main__":
    unittest.main()
