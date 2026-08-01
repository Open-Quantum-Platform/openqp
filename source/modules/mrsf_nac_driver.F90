!> Resident production driver for the analytic MRSF-TDDFT Lagrangian NAC.
!>
!> Every scientific operation for an ordered state pair remains in Fortran:
!> the exact eigenvector response y_IJ=X_I/(Omega_J-Omega_I), closed exact-TLF
!> metric column, amplitude/Fock skeletons, ROHF/ROKS one-RHS adjoint Z-vector,
!> coordinate contractions, pair accumulation, antisymmetrization, and gap
!> scaling.  Python invokes one C entry point and only reshapes the final data.
module mrsf_nac_driver_mod

  use precision, only: dp

  implicit none

  private
  public :: mrsf_nac_lagrangian

  character(len=*), parameter :: module_name = "mrsf_nac_driver_mod"

contains

!###############################################################################

  subroutine mrsf_nac_lagrangian_C(c_handle) &
      bind(C, name="mrsf_nac_lagrangian")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use io_constants, only: iw
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    logical :: log_was_open

    inf => oqp_handle_get_info(c_handle)
    ! The energy driver normally closes IW before NAC begins.  Keep it open
    ! around every resident kernel so no failure-free path creates fort.6.
    inquire(unit=iw, opened=log_was_open)
    if (.not. log_was_open) &
      open(unit=iw, file=inf%log_filename, position='append')
    call mrsf_nac_lagrangian(inf)
    if (.not. log_was_open) close(iw)
  end subroutine mrsf_nac_lagrangian_C

!###############################################################################

  subroutine mrsf_nac_lagrangian(infos)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use types, only: information
    use oqp_tagarray_driver, only: data_has_tags, tagarray_get_data, &
      OQP_td_bvec_mo, OQP_td_energies, TA_TYPE_REAL64
    use messages, only: show_message, WITH_ABORT
    use mrsf_nac_metric_data_mod, only: mrsf_nac_metric_column
    use tdhf_mrsf_gradient_mod, only: mrsf_nac_amp, mrsf_nac_esum
    use tdhf_mrsf_energy_mod, only: mrsf_nac_response
    use mrsf_nac_interchange_mod, only: &
      mrsf_nac_pair_accumulator_init, mrsf_nac_pair_accumulate, &
      mrsf_nac_pair_finalize, mrsf_nac_rohf_pair_overlap, &
      mrsf_nac_rohf_zvector, mrsf_nac_rohf_hf_adjoint, &
      mrsf_nac_xc_adjoint

    implicit none

    interface
      subroutine mrsf_nac_wpair_impl(infos, istate, jstate)
        use types, only: information
        type(information), target, intent(inout) :: infos
        integer, intent(in) :: istate, jstate
      end subroutine mrsf_nac_wpair_impl
    end interface

    character(len=*), parameter :: subroutine_name = &
      "mrsf_nac_lagrangian"
    character(len=*), parameter :: tag_ytil = "OQP::nac_ytil"
    character(len=*), parameter :: tag_xstate = "OQP::nac_xstate"
    character(len=*), parameter :: tag_gamma = "OQP::nac_gamma_pair"
    character(len=*), parameter :: tag_solution = &
      "OQP::nac_rohf_solution"
    character(len=*), parameter :: tag_z = "OQP::nac_rohf_z"
    character(len=*), parameter :: tags_required(2) = (/ character(len=80) :: &
      OQP_td_bvec_mo, OQP_td_energies /)

    type(information), target, intent(inout) :: infos
    real(kind=dp), contiguous, pointer :: bvec_mo(:,:), energies(:)
    real(kind=dp), contiguous, pointer :: solution(:)
    real(kind=dp), pointer :: ytil_tag(:), xstate_tag(:), gamma_tag(:), z_tag(:)
    real(kind=dp), allocatable :: bvec_saved(:,:), energies_saved(:), ytil(:)
    real(kind=dp), allocatable :: gamma_column(:,:)
    real(kind=dp), allocatable :: gamma_pair(:), z_work(:)
    real(kind=dp) :: gap, gap_floor, energy_scale, cutoff_saved
    integer(c_int64_t) :: nstate64, natom64, nbf64, noca64, nocb64
    integer(c_int64_t) :: nvirb64, nij64, nbfsq64, ncoord64
    integer(c_int64_t) :: state_pair_size64, default_int_limit64
    integer :: nbf, noca, nocb, nij, nstate, natom
    integer :: istate, jstate, redundant_index

    if (infos%control%scftype /= 3) then
      call show_message( &
        'Analytic MRSF NAC requires an ROHF/ROKS reference.', WITH_ABORT)
    end if
    if (infos%tddft%umrsf) then
      call show_message('Analytic MRSF NAC does not support UMRSF.', WITH_ABORT)
    end if
    if (infos%tddft%mult /= 1) then
      call show_message( &
        'Analytic MRSF NAC currently supports singlet states only.', WITH_ABORT)
    end if
    if (infos%control%conv > 1.0e-8_dp .or. &
        infos%tddft%cnvtol > 1.0e-8_dp) then
      call show_message( &
        'Analytic MRSF NAC requires SCF and TD thresholds <= 1e-8.', WITH_ABORT)
    end if
    if (.not. infos%mol_energy%SCF_converged .or. &
        .not. infos%mol_energy%Davidson_converged) then
      call show_message( &
        'Analytic MRSF NAC requires converged SCF and MRSF states.', WITH_ABORT)
    end if

    default_int_limit64 = int(huge(nstate),c_int64_t)
    nstate64 = infos%tddft%nstate
    if (nstate64 < 2_c_int64_t .or. &
        nstate64 > default_int_limit64) then
      call show_message('Invalid actual MRSF state count for NAC.', WITH_ABORT)
    end if
    natom64 = infos%mol_prop%natom
    noca64 = infos%mol_prop%nelec_A
    nocb64 = infos%mol_prop%nelec_B
    nbf64 = int(infos%basis%nbf,c_int64_t)
    if (natom64 < 1_c_int64_t .or. natom64 > default_int_limit64) then
      call show_message('Analytic MRSF NAC requires at least one atom.', &
                        WITH_ABORT)
    end if
    if (nbf64 < 1_c_int64_t .or. noca64 < 1_c_int64_t .or. &
        nocb64 < 1_c_int64_t .or. nocb64 >= noca64 .or. &
        noca64 > nbf64 .or. noca64-nocb64 /= 2_c_int64_t) then
      call show_message( &
        'Analytic MRSF NAC requires 1 <= nocb < noca <= nbf and two SOMOs.', &
        WITH_ABORT)
    end if
    if (natom64 > default_int_limit64/3_c_int64_t) then
      call show_message( &
        'MRSF NAC coordinate count exceeds the default-integer limit.', &
        WITH_ABORT)
    end if
    ncoord64 = 3_c_int64_t*natom64
    if (nstate64 > default_int_limit64/nstate64) then
      call show_message( &
        'MRSF NAC state-pair count exceeds the default-integer limit.', &
        WITH_ABORT)
    end if
    state_pair_size64 = nstate64*nstate64
    if (ncoord64 > default_int_limit64/state_pair_size64) then
      call show_message( &
        'MRSF NAC output exceeds the TagArray default-integer limit.', &
        WITH_ABORT)
    end if
    nvirb64 = nbf64-nocb64
    nij64 = noca64*nvirb64
    nbfsq64 = nbf64*nbf64
    if (nij64 > default_int_limit64 .or. &
        nbfsq64 > default_int_limit64) then
      call show_message( &
        'MRSF NAC orbital matrix exceeds the default-integer limit.', &
        WITH_ABORT)
    end if
    if (nstate64 > default_int_limit64/nij64 .or. &
        nstate64 > default_int_limit64/nbfsq64) then
      call show_message( &
        'MRSF NAC resident state data exceed the TagArray integer limit.', &
        WITH_ABORT)
    end if

    nstate = int(nstate64)
    natom = int(natom64)
    nbf = int(nbf64)
    noca = int(noca64)
    nocb = int(nocb64)
    nij = int(nij64)
    call data_has_tags(infos%dat, tags_required, module_name, &
                       subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    call tagarray_get_data(infos%dat, OQP_td_energies, energies)
    if (size(bvec_mo,1) /= nij .or. size(bvec_mo,2) /= nstate .or. &
        size(energies) /= nstate) then
      call show_message( &
        'Resident MRSF amplitudes/energies disagree with the actual state count.', &
        WITH_ABORT)
    end if
    if (.not. all(ieee_is_finite(bvec_mo)) .or. &
        .not. all(ieee_is_finite(energies))) then
      call show_message('Non-finite MRSF amplitude or energy in NAC driver.', &
                        WITH_ABORT)
    end if

    redundant_index = (noca-nocb-1)*noca+noca
    if (redundant_index < 1 .or. redundant_index > nij) then
      call show_message('Invalid redundant MRSF response coordinate.', &
                        WITH_ABORT)
    end if
    allocate(bvec_saved(nij,nstate), energies_saved(nstate), ytil(nij), &
             gamma_column(nbf*nbf,nstate), gamma_pair(nbf*nbf))
    bvec_saved = bvec_mo
    ! TagArray reserve/remove operations below may invalidate every cached
    ! record pointer, not only the record being changed.  Keep an owned copy
    ! of the state energies for the complete ordered-pair traversal.
    energies_saved = energies

    call infos%dat%remove_records((/ character(len=80) :: &
      tag_ytil, tag_xstate, tag_gamma, tag_z /))
    call infos%dat%reserve_data(tag_ytil, TA_TYPE_REAL64, nij, (/ nij /), &
      comment='streamed MRSF ordered-pair eigenvector response')
    call infos%dat%reserve_data(tag_xstate, TA_TYPE_REAL64, nij, (/ nij /), &
      comment='streamed MRSF right-state amplitude')
    call infos%dat%reserve_data(tag_gamma, TA_TYPE_REAL64, nbf*nbf, &
      (/ nbf*nbf /), comment='streamed exact-TLF pair metric source')

    cutoff_saved = infos%control%int2e_cutoff
    infos%control%int2e_cutoff = 1.0e-20_dp
    call mrsf_nac_pair_accumulator_init(infos)

    ! A fixed target column J shares its normalized-overlap denominator.  Build
    ! that O(nstate*nbf**2) metric column once, then consume each I immediately.
    do jstate = 1, nstate
      ! Ensure the metric always sees unmodified resident eigenvectors.
      call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
      bvec_mo = bvec_saved
      call mrsf_nac_metric_column(infos, jstate, gamma_column)

      do istate = 1, nstate
        if (istate == jstate) cycle
        gap = energies_saved(jstate)-energies_saved(istate)
        energy_scale = max(1.0_dp, abs(energies_saved(istate)), &
                           abs(energies_saved(jstate)))
        gap_floor = 128.0_dp*epsilon(1.0_dp)*energy_scale
        if (.not. ieee_is_finite(gap) .or. abs(gap) <= gap_floor) then
          call show_message( &
            'MRSF NAC state-response gap is zero or numerically unresolved.', &
            WITH_ABORT)
        end if
        ytil = bvec_saved(:,istate)/gap
        ytil(redundant_index) = 0.0_dp
        if (.not. all(ieee_is_finite(ytil))) then
          call show_message('Non-finite MRSF ordered-pair response.', &
                            WITH_ABORT)
        end if

        call tagarray_get_data(infos%dat, tag_ytil, ytil_tag)
        call tagarray_get_data(infos%dat, tag_xstate, xstate_tag)
        ytil_tag = ytil
        xstate_tag = bvec_saved(:,jstate)
        call mrsf_nac_wpair_impl(infos, istate, jstate)

        ! The pair amplitude engine reads the selected left response from its
        ! normal TD slot.  Reacquire the TagArray pointer before injection and
        ! again after kernels that reserve/remove other records.
        call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
        bvec_mo = bvec_saved
        bvec_mo(:,istate) = ytil
        call mrsf_nac_amp(infos, istate, jstate)
        call mrsf_nac_esum(infos, istate, jstate)
        call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
        bvec_mo = bvec_saved

        call mrsf_nac_response(infos)
        gamma_pair = gamma_column(:,istate)
        call tagarray_get_data(infos%dat, tag_gamma, gamma_tag)
        gamma_tag = gamma_pair
        call mrsf_nac_rohf_pair_overlap(infos)
        call mrsf_nac_rohf_zvector(infos)

        call tagarray_get_data(infos%dat, tag_solution, solution)
        if (.not. allocated(z_work)) allocate(z_work(size(solution)))
        if (size(z_work) /= size(solution)) then
          call show_message('Inconsistent streamed ROHF Z-vector dimension.', &
                            WITH_ABORT)
        end if
        z_work = solution
        call infos%dat%remove_records((/ character(len=80) :: tag_z /))
        call infos%dat%reserve_data(tag_z, TA_TYPE_REAL64, size(z_work), &
          (/ size(z_work) /), comment='current streamed pair ROHF adjoint')
        call tagarray_get_data(infos%dat, tag_z, z_tag)
        z_tag = z_work

        call mrsf_nac_rohf_hf_adjoint(infos)
        call mrsf_nac_xc_adjoint(infos)
        call mrsf_nac_pair_accumulate(infos, istate, jstate)
      end do
    end do

    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    bvec_mo = bvec_saved
    infos%control%int2e_cutoff = cutoff_saved
    call mrsf_nac_pair_finalize(infos)

    deallocate(bvec_saved, energies_saved, ytil, gamma_column, gamma_pair)
    if (allocated(z_work)) deallocate(z_work)
  end subroutine mrsf_nac_lagrangian

end module mrsf_nac_driver_mod
