"""Regression gates for MRSF Z/W Hessian intermediates."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_hessian_intermediates.F90"
ERI_DERIVATIVE = ROOT / "source/modules/tdhf_mrsf_hessian_eri_derivative.F90"
FOCK_DERIVATIVE = ROOT / "source/modules/fock_deriv.F90"
GRD2_RYS = ROOT / "source/integrals/grd2_rys.F90"
GRD2 = ROOT / "source/integrals/grd2.F90"


def test_apb_explicit_eri_derivative_uses_symmetric_density_convention():
    text = SOURCE.read_text()
    assert "int_apb=.true." in text
    assert "2.0_dp*explicit_pair(:,:,cart,atom,2*probe-1)" in text
    assert "2.0_dp*explicit_pair(:,:,cart,atom,2*probe)" in text
    assert "ordinary open-shell Fock derivative is therefore doubled" in text


def test_all_explicit_t_and_z_spin_focks_share_one_derivative_eri_traversal():
    text = SOURCE.read_text()
    start = text.index("  subroutine build_hf_response_fock_batch")
    stop = text.index("  end subroutine build_hf_response_fock_batch", start)
    routine = text[start:stop]
    assert "explicit_pcoul(nbf,nbf,2*nprobe)" in routine
    assert routine.count("call fock_deriv_matrix_os_batch") == 1
    assert "explicit_pexch(:,:,2*probe-1)=sym_a" in routine
    assert "explicit_pexch(:,:,2*probe)=sym_b" in routine


def test_response_fock_contracts_all_nuclear_rhs_in_one_eri_pass():
    text = SOURCE.read_text()
    start = text.index("  subroutine build_hf_response_fock_batch")
    stop = text.index("  end subroutine build_hf_response_fock_batch", start)
    routine = text[start:stop]
    assert "density(nbf,nbf,2*(ncoord+1)*nprobe)" in routine
    assert "density(:,:,pair_index+2*coordinate+1)=dpa(:,:,coordinate,probe)" in routine
    assert "density(:,:,pair_index+2*coordinate+2)=dpb(:,:,coordinate,probe)" in routine
    assert routine.count("call int2_driver%run(int2_data)") == 1


def test_t_and_z_probes_share_one_exact_response_fock_traversal():
    text = SOURCE.read_text()
    start = text.index("  subroutine build_mrsf_hf_z_intermediates")
    stop = text.index("  end subroutine build_mrsf_hf_z_intermediates", start)
    routine = text[start:stop]
    assert "probe_pa(nbf,nbf,2)" in routine
    assert routine.count("call build_hf_response_fock_batch") == 1
    assert "probe_pa(:,:,1)=ta" in routine
    assert "probe_pa(:,:,2)=matmul(mo_a" in routine


def test_reconstructed_gradient_intermediates_fail_closed():
    text = SOURCE.read_text()
    assert "max(baseline_error_a,baseline_error_b)>1.0e-8_dp" in text
    assert "data%z_rhs_hxa-(stored_hxa-2.0_dp*" in text
    assert "data%ppija-stored_ppija" in text
    assert "data%ppijb-stored_ppijb" in text


def test_no_unconditional_debug_files_or_diagnostic_stderr():
    text = SOURCE.read_text()
    assert "/private/tmp" not in text
    assert "error_unit" not in text


def test_mrsf_ordered_eri_derivative_is_direct_and_never_uses_ao_probes():
    caller = ERI_DERIVATIVE.read_text()
    start = caller.index("  subroutine build_mrsf_explicit_eri_derivative")
    stop = caller.index("  end subroutine build_mrsf_explicit_eri_derivative", start)
    routine = caller[start:stop]
    assert routine.count("call fock_deriv_matrix_mrsf_scaled_batch") == 1
    assert "call fock_deriv_matrix_mrsf_scaled(" not in routine

    direct = FOCK_DERIVATIVE.read_text()
    start = direct.index("  subroutine fock_deriv_matrix_mrsf_scaled_batch")
    stop = direct.index("  end subroutine fock_deriv_matrix_mrsf_scaled_batch", start)
    routine = direct[start:stop]
    assert routine.count("call grd2_fock_deriv_mrsf_driver_batch") == 1
    assert "call fock_deriv_matrix_mrsf_scaled(" not in routine


def test_direct_ordered_rys_adjoint_keeps_all_four_coulomb_and_eight_exchange_targets():
    text = GRD2_RYS.read_text()
    start = text.index("              else\n                ! Exact ordered-density adjoint")
    stop = text.index("              end if\n            end do", start)
    ordered = text[start:stop]
    assert ordered.count("call add_target_fock") == 12
    assert ordered.count(".false.)") == 12


def test_derivative_fock_openmp_accumulators_are_heap_allocated():
    text = GRD2.read_text()
    for name in (
        "grd2_fock_deriv_os_driver",
        "grd2_fock_deriv_os_driver_batch",
        "grd2_fock_deriv_mrsf_driver_batch",
    ):
        start = text.index(f"  subroutine {name}")
        stop = text.index(f"  end subroutine {name}", start)
        routine = text[start:stop]
        assert "reduction(+:skip1,skip2,numint,fmat)" not in routine
        assert "fmat_threads" in routine
        assert "omp_get_max_threads" in routine
        assert "omp_get_thread_num()+1" in routine
        assert "source=0.0_dp" in routine
        assert "fmat=fmat+fmat_threads" in routine
        assert "deallocate(fmat_threads)" in routine
