!> Resident production driver for the analytic MRSF-TDDFT Lagrangian NAC.
!>
!> Every scientific operation for an ordered state pair remains in Fortran:
!> the exact eigenvector response y_IJ=X_I/(Omega_J-Omega_I), closed exact-TLF
!> metric column and amplitude/Fock skeletons.  The two ordered sources for
!> each physical pair are reduced to their exact half-difference.  Unordered
!> ROHF/ROKS adjoints and their HF/XC contractions are processed in bounded
!> batches of at most three pairs; the production three-state case therefore
!> uses one batch for all three physical pairs.  Python invokes one C entry
!> point and only reshapes the final data.
module mrsf_nac_driver_mod

  use precision, only: dp

  implicit none

  private
  public :: mrsf_nac_lagrangian, mrsf_nac_lagrangian_fused_buffered

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

  subroutine mrsf_nac_lagrangian_fused_buffered(infos)
    use types, only: information
    use mrsf_nac_fusion_buffer_mod, only: mrsf_nac_fusion_get_rhs, &
      mrsf_nac_fusion_set_solution

    type(information), target, intent(inout) :: infos
    real(kind=dp), allocatable :: gradient_rhs(:), gradient_solution(:)
    integer :: nbf, nocca, noccb, ltot

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    ltot = (nocca-noccb)*noccb + (nbf-nocca)*noccb + &
           (nbf-nocca)*(nocca-noccb)
    allocate(gradient_rhs(ltot), gradient_solution(ltot))
    call mrsf_nac_fusion_get_rhs(gradient_rhs)
    call mrsf_nac_lagrangian(infos,gradient_rhs,gradient_solution)
    call mrsf_nac_fusion_set_solution(gradient_solution)
    deallocate(gradient_rhs,gradient_solution)
  end subroutine mrsf_nac_lagrangian_fused_buffered

!###############################################################################

  subroutine mrsf_nac_lagrangian(infos, gradient_rhs, gradient_solution)
    use mrsf_nac_fusion_buffer_mod, only: mrsf_nac_fusion_get_tolerance
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use types, only: information
    use io_constants, only: iw
    use oqp_tagarray_driver, only: tagarray_reserve_data, data_has_tags, tagarray_get_data, &
      OQP_td_bvec_mo, OQP_td_energies, TA_TYPE_REAL64
    use messages, only: show_message, WITH_ABORT
    use cphf_mod, only: cphf_solve_rohf
    use mrsf_nac_metric_data_mod, only: mrsf_nac_metric_column
    use tdhf_mrsf_gradient_mod, only: mrsf_nac_amp, mrsf_nac_esum
    use tdhf_mrsf_energy_mod, only: mrsf_nac_response
    use mrsf_nac_interchange_mod, only: &
      mrsf_nac_pair_accumulator_init, &
      mrsf_nac_pair_accumulate_antisym, &
      mrsf_nac_pair_finalize, mrsf_nac_rohf_pair_overlap, &
      mrsf_nac_rohf_zvector_batch, mrsf_nac_rohf_hf_adjoint_batch, &
      mrsf_nac_xc_adjoint_batch

    implicit none

    interface
      subroutine mrsf_nac_wpair_batch_impl(infos, ytil_batch, &
                                           xstate_batch, mt_batch)
        use types, only: information
        use precision, only: dp
        type(information), target, intent(inout) :: infos
        real(kind=dp), contiguous, intent(in) :: ytil_batch(:,:), &
                                                 xstate_batch(:,:)
        real(kind=dp), contiguous, intent(out) :: mt_batch(:,:,:)
      end subroutine mrsf_nac_wpair_batch_impl
    end interface

    character(len=*), parameter :: subroutine_name = &
      "mrsf_nac_lagrangian"
    integer, parameter :: z_batch_width = 3
    integer, parameter :: wpair_batch_width = 3
    integer, parameter :: hf_batch_width = 3
    integer, parameter :: xc_batch_width = 3
    character(len=*), parameter :: tag_ytil = "OQP::nac_ytil"
    character(len=*), parameter :: tag_xstate = "OQP::nac_xstate"
    character(len=*), parameter :: tag_gamma = "OQP::nac_gamma_pair"
    character(len=*), parameter :: tag_rhs = "OQP::nac_rohf_rhs"
    character(len=*), parameter :: tag_amp = "OQP::nac_amp"
    character(len=*), parameter :: tag_esum = "OQP::nac_esum"
    character(len=*), parameter :: tag_overlap = "OQP::nac_pair_overlap"
    character(len=*), parameter :: tag_mt_frozen = "OQP::nac_mt_frozen"
    character(len=*), parameter :: tag_z = "OQP::nac_rohf_z"
    character(len=*), parameter :: tag_hf = "OQP::nac_rohf_hf_adjoint"
    character(len=*), parameter :: tag_xc = "OQP::nac_rohf_xc_adjoint"
    character(len=*), parameter :: tag_predictor_dcv = &
      "OQP::nac_predictor_dcv"
    character(len=*), parameter :: tag_predictor_nacv = &
      "OQP::nac_predictor_nacv"
    character(len=*), parameter :: tags_required(2) = (/ character(len=80) :: &
      OQP_td_bvec_mo, OQP_td_energies /)

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in), optional :: gradient_rhs(:)
    real(kind=dp), intent(out), optional :: gradient_solution(:)
    real(kind=dp), contiguous, pointer :: bvec_mo(:,:), energies(:)
    real(kind=dp), contiguous, pointer :: rhs_in(:), amp(:,:,:), esum(:,:), &
      pair_overlap(:,:)
    real(kind=dp), pointer :: ytil_tag(:), xstate_tag(:), gamma_tag(:), &
      mt_frozen_tag(:), z_tag(:), hf_tag(:,:), xc_tag(:,:), &
      predictor_dcv_tag(:,:,:), predictor_nacv_tag(:,:,:)
    real(kind=dp), allocatable :: bvec_saved(:,:), energies_saved(:), ytil(:)
    real(kind=dp), allocatable :: gamma_column(:,:)
    real(kind=dp), allocatable :: gamma_pair(:), rhs_batch(:,:), &
      solution_batch(:,:), predictor_batch(:,:), nonz_batch(:,:), &
      hf_batch(:,:,:), xc_batch(:,:,:), hf_predictor(:,:,:), &
      xc_predictor(:,:,:), predictor_initial_residual(:), &
      predictor_final_residual(:)
    real(kind=dp), allocatable :: fused_rhs(:,:), fused_solution(:,:), &
      fused_residual(:), fused_tolerance(:)
    real(kind=dp), allocatable :: wpair_ytil(:,:), wpair_xstate(:,:), &
      wpair_mt(:,:,:)
    integer, allocatable :: pair_i(:), pair_j(:)
    integer, allocatable :: predictor_iterations(:)
    logical, allocatable :: predictor_available(:), predictor_accepted(:), &
      fused_converged(:)
    integer, allocatable :: fused_iterations(:)
    real(kind=dp) :: gap, gap_floor, energy_scale, cutoff_saved, pair_sign, &
      fused_scaled_residual
    real(kind=dp) :: profile_total, profile_metric, profile_wpair
    real(kind=dp) :: profile_amp, profile_esum, profile_response
    real(kind=dp) :: profile_overlap, profile_zvector, profile_hf
    real(kind=dp) :: profile_xc, profile_accumulate, profile_finalize
    integer(c_int64_t) :: nstate64, natom64, nbf64, noca64, nocb64
    integer(c_int64_t) :: nvirb64, nij64, nbfsq64, ncoord64
    integer(c_int64_t) :: state_pair_size64, default_int_limit64
    integer :: nbf, noca, nocb, nij, nstate, natom, ncoord
    integer :: nvira, offset, ltot, npair, ipair, z_first, z_last, &
      hf_first, hf_last, xc_first, xc_last
    integer :: wpair_first, wpair_last, wpair_count, wpair_index, batch_pair
    integer :: istate, jstate, redundant_index, atom, cart, coord
    integer(c_int64_t) :: profile_start, profile_stop, profile_rate
    integer :: profile_status
    character(len=16) :: profile_value, audit_value
    logical :: profile_enabled, audit_enabled
    real(kind=dp) :: gradient_tolerance

    if (present(gradient_rhs) .neqv. present(gradient_solution)) then
      call show_message('Gradient/NAC fusion requires both RHS and solution arrays.', &
                        WITH_ABORT)
    end if

    profile_value = ''
    call get_environment_variable('OQP_NAC_PROFILE', profile_value, &
                                  status=profile_status)
    profile_enabled = profile_status == 0 .and. &
      len_trim(profile_value) > 0 .and. trim(profile_value) /= '0'
    audit_value = ''
    call get_environment_variable('OQP_MRSF_NAC_ZV_AUDIT',audit_value)
    audit_enabled = len_trim(audit_value) > 0 .and. &
                    trim(audit_value) /= '0'
    profile_total = 0.0_dp
    profile_metric = 0.0_dp
    profile_wpair = 0.0_dp
    profile_amp = 0.0_dp
    profile_esum = 0.0_dp
    profile_response = 0.0_dp
    profile_overlap = 0.0_dp
    profile_zvector = 0.0_dp
    profile_hf = 0.0_dp
    profile_xc = 0.0_dp
    profile_accumulate = 0.0_dp
    profile_finalize = 0.0_dp
    if (profile_enabled) then
      call system_clock(profile_start, profile_rate)
      if (profile_rate <= 0_c_int64_t) profile_enabled = .false.
    end if

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
        nocb64 < 0_c_int64_t .or. nocb64 >= noca64 .or. &
        noca64 > nbf64 .or. noca64-nocb64 /= 2_c_int64_t) then
      call show_message( &
        'Analytic MRSF NAC requires 0 <= nocb < noca <= nbf and two SOMOs.', &
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
    ncoord = 3*natom
    nvira = nbf - noca
    offset = noca - nocb
    ltot = nocb*(offset + nvira) + offset*nvira
    npair = nstate*(nstate - 1)/2
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
             gamma_column(nbf*nbf,nstate), gamma_pair(nbf*nbf), &
             rhs_batch(ltot,npair), solution_batch(ltot,npair), &
             predictor_batch(ltot,npair), &
             nonz_batch(ncoord,npair), hf_batch(3,natom,npair), &
             xc_batch(3,natom,npair), hf_predictor(3,natom,npair), &
             xc_predictor(3,natom,npair), &
             predictor_initial_residual(npair), &
             predictor_final_residual(npair), predictor_iterations(npair), &
             predictor_available(npair), predictor_accepted(npair), &
             wpair_ytil(nij,wpair_batch_width), &
             wpair_xstate(nij,wpair_batch_width), &
             wpair_mt(nbf,nbf,wpair_batch_width), &
             pair_i(npair), pair_j(npair))
    bvec_saved = bvec_mo
    ! TagArray reserve/remove operations below may invalidate every cached
    ! record pointer, not only the record being changed.  Keep an owned copy
    ! of the state energies for the complete ordered-pair traversal.
    energies_saved = energies
    rhs_batch = 0.0_dp
    solution_batch = 0.0_dp
    predictor_batch = 0.0_dp
    nonz_batch = 0.0_dp
    hf_batch = 0.0_dp
    xc_batch = 0.0_dp
    hf_predictor = 0.0_dp
    xc_predictor = 0.0_dp
    predictor_initial_residual = 0.0_dp
    predictor_final_residual = 0.0_dp
    predictor_iterations = 0
    predictor_available = .false.
    predictor_accepted = .false.
    wpair_first = 1
    wpair_last = 0
    ipair = 0
    do jstate = 2, nstate
      do istate = 1, jstate - 1
        ipair = ipair + 1
        pair_i(ipair) = istate
        pair_j(ipair) = jstate
      end do
    end do

    call infos%dat%erase((/ character(len=80) :: &
      tag_ytil, tag_xstate, tag_gamma, tag_z /))
    call tagarray_reserve_data(infos%dat, tag_ytil, TA_TYPE_REAL64, nij, (/ nij /), &
      comment='streamed MRSF ordered-pair eigenvector response')
    call tagarray_reserve_data(infos%dat, tag_xstate, TA_TYPE_REAL64, nij, (/ nij /), &
      comment='streamed MRSF right-state amplitude')
    call tagarray_reserve_data(infos%dat, tag_gamma, TA_TYPE_REAL64, nbf*nbf, &
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
      if (profile_enabled) call system_clock(profile_stop)
      call mrsf_nac_metric_column(infos, jstate, gamma_column)
      if (profile_enabled) call profile_add(profile_metric, profile_stop)

      do istate = 1, nstate
        if (istate == jstate) cycle
        ipair = unordered_pair_index(istate, jstate)
        pair_sign = merge(0.5_dp, -0.5_dp, istate < jstate)
        gamma_pair = pair_sign*gamma_column(:,istate)

        ! Every direct pair kernel is a symmetric bilinear form in its left
        ! and right MRSF amplitudes.  Reversing the ordered pair also reverses
        ! the energy gap, so its direct source is exactly the negative of this
        ! canonical I<J source.  Evaluate that expensive source only once.
        if (istate < jstate) then
          if (ipair > wpair_last) then
            wpair_first = ipair
            wpair_last = min(npair, wpair_first + wpair_batch_width - 1)
            wpair_count = wpair_last - wpair_first + 1
            do wpair_index = 1, wpair_count
              batch_pair = wpair_first + wpair_index - 1
              gap = energies_saved(pair_j(batch_pair)) - &
                    energies_saved(pair_i(batch_pair))
              energy_scale = max(1.0_dp, &
                abs(energies_saved(pair_i(batch_pair))), &
                abs(energies_saved(pair_j(batch_pair))))
              gap_floor = 128.0_dp*epsilon(1.0_dp)*energy_scale
              if (.not. ieee_is_finite(gap) .or. &
                  abs(gap) <= gap_floor) then
                call show_message( &
                  'MRSF NAC state-response gap is zero or numerically unresolved.', &
                  WITH_ABORT)
              end if
              wpair_ytil(:,wpair_index) = &
                bvec_saved(:,pair_i(batch_pair))/gap
              wpair_ytil(redundant_index,wpair_index) = 0.0_dp
              wpair_xstate(:,wpair_index) = &
                bvec_saved(:,pair_j(batch_pair))
              if (.not. all(ieee_is_finite( &
                    wpair_ytil(:,wpair_index)))) then
                call show_message('Non-finite MRSF ordered-pair response.', &
                                  WITH_ABORT)
              end if
            end do
            if (profile_enabled) call system_clock(profile_stop)
            call mrsf_nac_wpair_batch_impl( &
              infos, wpair_ytil(:,1:wpair_count), &
              wpair_xstate(:,1:wpair_count), &
              wpair_mt(:,:,1:wpair_count))
            if (profile_enabled) call profile_add(profile_wpair, profile_stop)
          end if
          wpair_index = ipair - wpair_first + 1
          ytil = wpair_ytil(:,wpair_index)

          call tagarray_get_data(infos%dat, tag_ytil, ytil_tag)
          call tagarray_get_data(infos%dat, tag_xstate, xstate_tag)
          ytil_tag = ytil
          xstate_tag = bvec_saved(:,jstate)
          ! Publish only the current pair.  Downstream overlap assembly sees
          ! the same record at the same point as in the scalar implementation.
          call infos%dat%erase((/ character(len=80) :: &
            tag_mt_frozen /))
          call tagarray_reserve_data(infos%dat, tag_mt_frozen, TA_TYPE_REAL64, &
            nbf*nbf, (/ nbf*nbf /), &
            comment='current batched MRSF frozen pair orbital source')
          call tagarray_get_data(infos%dat, tag_mt_frozen, mt_frozen_tag)
          mt_frozen_tag = reshape(wpair_mt(:,:,wpair_index), (/ nbf*nbf /))

          ! The pair amplitude engine reads the selected left response from its
          ! normal TD slot. Reacquire the TagArray pointer before injection and
          ! again after kernels that reserve/remove other records.
          call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
          bvec_mo = bvec_saved
          bvec_mo(:,istate) = ytil
          if (profile_enabled) call system_clock(profile_stop)
          call mrsf_nac_amp(infos, istate, jstate)
          if (profile_enabled) call profile_add(profile_amp, profile_stop)
          if (profile_enabled) call system_clock(profile_stop)
          call mrsf_nac_esum(infos, istate, jstate)
          if (profile_enabled) call profile_add(profile_esum, profile_stop)
          call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
          bvec_mo = bvec_saved

          if (profile_enabled) call system_clock(profile_stop)
          call mrsf_nac_response(infos)
          if (profile_enabled) call profile_add(profile_response, profile_stop)
        end if

        call tagarray_get_data(infos%dat, tag_gamma, gamma_tag)
        gamma_tag = gamma_pair
        if (profile_enabled) call system_clock(profile_stop)
        if (istate < jstate) then
          ! D + gamma_IJ/2, where D is the canonical direct source.
          call mrsf_nac_rohf_pair_overlap(infos)
        else
          ! The reverse direct source is -D and has already been folded into
          ! the canonical contribution. Add only -gamma_JI/2 here.
          call mrsf_nac_rohf_pair_overlap(infos, metric_only=.true.)
        end if
        if (profile_enabled) call profile_add(profile_overlap, profile_stop)

        ! The scaled metric contributions are folded as each target column is
        ! streamed.  No nbf**2 unordered-pair tensor is materialized.
        call tagarray_get_data(infos%dat, tag_rhs, rhs_in)
        call tagarray_get_data(infos%dat, tag_overlap, pair_overlap)
        if (size(rhs_in) /= ltot .or. size(pair_overlap,1) /= 3 .or. &
            size(pair_overlap,2) /= natom) then
          call show_message( &
            'Ordered MRSF NAC sources have inconsistent dimensions.', &
            WITH_ABORT)
        end if
        rhs_batch(:,ipair) = rhs_batch(:,ipair) + rhs_in
        if (istate < jstate) then
          call tagarray_get_data(infos%dat, tag_amp, amp)
          call tagarray_get_data(infos%dat, tag_esum, esum)
          if (size(amp,1) /= ncoord .or. size(amp,2) /= nstate .or. &
              size(amp,3) /= nstate .or. size(esum,1) /= 3 .or. &
              size(esum,2) /= natom) then
            call show_message( &
              'Direct MRSF NAC sources have inconsistent dimensions.', &
              WITH_ABORT)
          end if
          do atom = 1, natom
            do cart = 1, 3
              coord = (atom - 1)*3 + cart
              nonz_batch(coord,ipair) = nonz_batch(coord,ipair) + &
                amp(coord,istate,jstate) + esum(cart,atom) + &
                pair_overlap(cart,atom)
            end do
          end do
        else
          do atom = 1, natom
            do cart = 1, 3
              coord = (atom - 1)*3 + cart
              nonz_batch(coord,ipair) = nonz_batch(coord,ipair) + &
                pair_overlap(cart,atom)
            end do
          end do
        end if
      end do
    end do

    ! Bound the solver-owned nbf**2 multi-vector/Fock workspace for callers
    ! requesting many states.  The production three-state case remains one
    ! nrhs=3 solve, in place of six independent ordered-pair solves.
    if (profile_enabled) call system_clock(profile_stop)
    if (present(gradient_rhs)) then
      call mrsf_nac_fusion_get_tolerance(gradient_tolerance)
      if (size(gradient_rhs) /= ltot .or. size(gradient_solution) /= ltot) then
        call show_message('Gradient/NAC fusion has the wrong ROHF rotation dimension.', &
                          WITH_ABORT)
      end if
      allocate(fused_rhs(ltot,npair+1), fused_solution(ltot,npair+1), &
               fused_residual(npair+1), fused_tolerance(npair+1), &
               fused_converged(npair+1), &
               fused_iterations(npair+1))
      ! The native Hessian H is twice the legacy MRSF-gradient Hessian A in
      ! the shared SD/DV/SV coordinates. H*x_g=2*b_g therefore returns the
      ! exact legacy gradient multiplier, while the remaining columns are the
      ! native unordered-pair adjoints. One synchronized multi-RHS traversal
      ! serves the active-state gradient and all three-state NAC pairs.
      fused_rhs(:,1) = 2.0_dp*gradient_rhs
      fused_rhs(:,2:npair+1) = rhs_batch
      fused_tolerance = 1.0e-20_dp
      ! Keep the established gradient residual threshold while the NAC pair
      ! columns retain the tighter property threshold.  Once the gradient is
      ! converged, its density/Fock column drops out of the shared traversal.
      fused_tolerance(1) = max(1.0e-20_dp, gradient_tolerance)
      call cphf_solve_rohf(infos,npair+1,fused_rhs,fused_solution, &
                           tol=1.0e-20_dp, &
                           maxit=max(int(infos%control%maxit_zv),ltot+5), &
                           converged=fused_converged, &
                           residual=fused_residual, minres_solver=.true., &
                           iterations=fused_iterations, &
                           rhs_tolerances=fused_tolerance)
      do z_first = 1, npair+1
        fused_scaled_residual = fused_residual(z_first)/max(1.0_dp, &
          sum(fused_rhs(:,z_first)*fused_rhs(:,z_first)))
        write(iw,'(A,1X,I0,2(1X,ES16.8),1X,I0)') &
          'NAC_GRADIENT_Z_FUSION',z_first,fused_residual(z_first), &
          fused_scaled_residual,fused_iterations(z_first)
        if (.not. fused_converged(z_first) .and. fused_scaled_residual > &
            (10.0_dp*sqrt(epsilon(1.0_dp)))**2) then
          call show_message('Fused gradient/NAC ROHF Z-vector failed residual gate.', &
                            WITH_ABORT)
        end if
      end do
      gradient_solution = fused_solution(:,1)
      solution_batch = fused_solution(:,2:npair+1)
      predictor_available = .false.
      predictor_accepted = .false.
      predictor_batch = 0.0_dp
      predictor_initial_residual = 0.0_dp
      predictor_final_residual = 0.0_dp
      predictor_iterations = 0
      deallocate(fused_rhs,fused_solution,fused_residual,fused_tolerance, &
                 fused_converged,fused_iterations)
    else
      do z_first = 1, npair, z_batch_width
        z_last = min(npair, z_first + z_batch_width - 1)
        call mrsf_nac_rohf_zvector_batch( &
          infos, rhs_batch(:,z_first:z_last), &
          solution_batch(:,z_first:z_last), pair_offset=z_first-1, &
          predictor=predictor_batch(:,z_first:z_last), &
          predictor_available=predictor_available(z_first:z_last), &
          predictor_accepted=predictor_accepted(z_first:z_last), &
          initial_residual_out= &
            predictor_initial_residual(z_first:z_last), &
          final_residual_out=predictor_final_residual(z_first:z_last), &
          iterations_out=predictor_iterations(z_first:z_last))
      end do
    end if
    if (profile_enabled) call profile_add(profile_zvector, profile_stop)

    if (profile_enabled) call system_clock(profile_stop)
    ! Bound the AO response-density workspace for unusually many states.  The
    ! production three-state case shares all ground-density, one-electron
    ! derivative and AO->MO work across its three physical pairs.
    do hf_first = 1, npair, hf_batch_width
      hf_last = min(npair, hf_first + hf_batch_width - 1)
      call mrsf_nac_rohf_hf_adjoint_batch( &
        infos, solution_batch(:,hf_first:hf_last), &
        hf_batch(:,:,hf_first:hf_last))
    end do
    if (profile_enabled) call profile_add(profile_hf, profile_stop)

    if (profile_enabled) call system_clock(profile_stop)
    ! Bound the grid-consumer workspace for callers requesting many states.
    ! The production three-state case still traverses the grid only once.
    do xc_first = 1, npair, xc_batch_width
      xc_last = min(npair, xc_first + xc_batch_width - 1)
      call mrsf_nac_xc_adjoint_batch( &
        infos, solution_batch(:,xc_first:xc_last), &
        xc_batch(:,:,xc_first:xc_last))
    end do
    if (profile_enabled) call profile_add(profile_xc, profile_stop)

    ! The optional audit evaluates the transported or extrapolated Z vector
    ! before MINRES correction with the same analytic adjoint contractions.
    ! It is an observational comparison only: the exact solution above remains
    ! the source of the production NAC and trajectory coupling.
    if (audit_enabled .and. any(predictor_available)) then
      do hf_first = 1, npair, hf_batch_width
        hf_last = min(npair, hf_first + hf_batch_width - 1)
        call mrsf_nac_rohf_hf_adjoint_batch( &
          infos, predictor_batch(:,hf_first:hf_last), &
          hf_predictor(:,:,hf_first:hf_last))
      end do
      do xc_first = 1, npair, xc_batch_width
        xc_last = min(npair, xc_first + xc_batch_width - 1)
        call mrsf_nac_xc_adjoint_batch( &
          infos, predictor_batch(:,xc_first:xc_last), &
          xc_predictor(:,:,xc_first:xc_last))
      end do
    end if

    call infos%dat%erase((/ character(len=80) :: tag_z, tag_hf, &
                                                        tag_xc /))
    call tagarray_reserve_data(infos%dat, tag_z, TA_TYPE_REAL64, ltot, (/ ltot /), &
      comment='current antisymmetric unordered-pair ROHF adjoint')
    call tagarray_reserve_data(infos%dat, tag_hf, TA_TYPE_REAL64, 3*natom, &
      (/ 3, natom /), comment='batched native ROHF NAC analytic HF adjoint')
    call tagarray_reserve_data(infos%dat, tag_xc, TA_TYPE_REAL64, 3*natom, &
      (/ 3, natom /), comment='batched native ROHF NAC analytic XC adjoint')
    do ipair = 1, npair
      ! Every adjoint contraction is linear in z.  Applying it once to the
      ! half-difference solution is therefore exactly the half-difference of
      ! the two ordered adjoints, up to the solver's certified residual.
      call tagarray_get_data(infos%dat, tag_z, z_tag)
      z_tag = solution_batch(:,ipair)
      call tagarray_get_data(infos%dat, tag_hf, hf_tag)
      hf_tag = hf_batch(:,:,ipair)
      call tagarray_get_data(infos%dat, tag_xc, xc_tag)
      xc_tag = xc_batch(:,:,ipair)
      if (profile_enabled) call system_clock(profile_stop)
      call mrsf_nac_pair_accumulate_antisym( &
        infos, pair_i(ipair), pair_j(ipair), nonz_batch(:,ipair))
      if (profile_enabled) call profile_add(profile_accumulate, profile_stop)
    end do

    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    bvec_mo = bvec_saved
    infos%control%int2e_cutoff = cutoff_saved
    if (profile_enabled) call system_clock(profile_stop)
    call mrsf_nac_pair_finalize(infos)
    call infos%dat%erase((/ character(len=80) :: tag_predictor_dcv, &
      tag_predictor_nacv /))
    if (audit_enabled .and. any(predictor_available)) then
      call tagarray_reserve_data(infos%dat,tag_predictor_dcv, &
        TA_TYPE_REAL64,ncoord*nstate*nstate,(/ ncoord,nstate,nstate /), &
        comment='uncorrected transported/extrapolated MRSF derivative coupling')
      call tagarray_reserve_data(infos%dat,tag_predictor_nacv, &
        TA_TYPE_REAL64,ncoord*nstate*nstate,(/ ncoord,nstate,nstate /), &
        comment='uncorrected transported/extrapolated gap-scaled coupling')
      call tagarray_get_data(infos%dat,tag_predictor_dcv,predictor_dcv_tag)
      call tagarray_get_data(infos%dat,tag_predictor_nacv,predictor_nacv_tag)
      predictor_dcv_tag = 0.0_dp
      predictor_nacv_tag = 0.0_dp
      do ipair = 1, npair
        if (.not. predictor_available(ipair)) cycle
        gap = energies_saved(pair_j(ipair))-energies_saved(pair_i(ipair))
        do coord = 1, ncoord
          predictor_dcv_tag(coord,pair_i(ipair),pair_j(ipair)) = &
            nonz_batch(coord,ipair) + &
            hf_predictor(mod(coord-1,3)+1,(coord-1)/3+1,ipair) + &
            xc_predictor(mod(coord-1,3)+1,(coord-1)/3+1,ipair)
          predictor_dcv_tag(coord,pair_j(ipair),pair_i(ipair)) = &
            -predictor_dcv_tag(coord,pair_i(ipair),pair_j(ipair))
          predictor_nacv_tag(coord,pair_i(ipair),pair_j(ipair)) = &
            gap*predictor_dcv_tag(coord,pair_i(ipair),pair_j(ipair))
          predictor_nacv_tag(coord,pair_j(ipair),pair_i(ipair)) = &
            predictor_nacv_tag(coord,pair_i(ipair),pair_j(ipair))
        end do
      end do
    end if
    if (profile_enabled) then
      call profile_add(profile_finalize, profile_stop)
      call system_clock(profile_stop)
      profile_total = real(profile_stop-profile_start,dp)/real(profile_rate,dp)
      write(iw,'(A,12(1X,A,"=",F12.6))') 'NAC_PROFILE', &
        'total', profile_total, 'metric', profile_metric, &
        'wpair', profile_wpair, 'amp', profile_amp, &
        'esum', profile_esum, 'response', profile_response, &
        'overlap', profile_overlap, 'zvector', profile_zvector, &
        'hf', profile_hf, 'xc', profile_xc, &
        'accumulate', profile_accumulate, 'finalize', profile_finalize
      flush(iw)
    end if

    deallocate(bvec_saved, energies_saved, ytil, gamma_column, gamma_pair, &
               rhs_batch, solution_batch, predictor_batch, nonz_batch, &
               hf_batch, xc_batch, hf_predictor, xc_predictor, &
               predictor_initial_residual, predictor_final_residual, &
               predictor_iterations, predictor_available, &
               predictor_accepted, &
               wpair_ytil, wpair_xstate, wpair_mt, &
               pair_i, pair_j)
  contains
    pure integer function unordered_pair_index(left_state, right_state) &
        result(index)
      integer, intent(in) :: left_state, right_state
      integer :: lo, hi

      lo = min(left_state, right_state)
      hi = max(left_state, right_state)
      index = (hi - 1)*(hi - 2)/2 + lo
    end function unordered_pair_index

    subroutine profile_add(accumulator, start_count)
      real(kind=dp), intent(inout) :: accumulator
      integer(c_int64_t), intent(in) :: start_count
      integer(c_int64_t) :: end_count

      call system_clock(end_count)
      accumulator = accumulator + &
        real(end_count-start_count,dp)/real(profile_rate,dp)
    end subroutine profile_add
  end subroutine mrsf_nac_lagrangian

end module mrsf_nac_driver_mod

! External bridge used by the legacy MRSF gradient module without introducing
! a Fortran module-dependency cycle (the NAC driver itself consumes resident
! kernels from tdhf_mrsf_gradient_mod).
subroutine mrsf_nac_lagrangian_fused_external(infos)
  use types, only: information
  use mrsf_nac_driver_mod, only: mrsf_nac_lagrangian_fused_buffered
  implicit none
  type(information), target, intent(inout) :: infos
  call mrsf_nac_lagrangian_fused_buffered(infos)
end subroutine mrsf_nac_lagrangian_fused_external
