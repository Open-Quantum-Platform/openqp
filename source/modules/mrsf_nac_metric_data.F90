!> Closed-form resident construction of the MRSF state-overlap metric source.
!>
!> The exact ndtlf=0 overlap contains seven determinant-block contractions.
!> This module differentiates those contractions, their determinant minors,
!> and the final column normalization analytically at the identity MO overlap.
!> No displaced geometry or finite-difference sweep is used.
module mrsf_nac_metric_data_mod

  use precision, only: dp

  implicit none

  private
  public :: mrsf_nac_metric_data
  public :: mrsf_nac_metric_column

  character(len=*), parameter :: module_name = "mrsf_nac_metric_data_mod"
  character(len=*), parameter :: tag_gamma = "OQP::nac_gamma_tlf"

contains

!###############################################################################

  subroutine mrsf_nac_metric_data_C(c_handle) &
      bind(C, name="mrsf_nac_metric_data")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call mrsf_nac_metric_data(inf)
  end subroutine mrsf_nac_metric_data_C

!###############################################################################

!> Diagnostic C wrapper for comparing the streamed production column against
!> the all-pair exact metric oracle.  Production calls the internal routine
!> directly and never materializes OQP::nac_gamma_tlf.
  subroutine mrsf_nac_metric_column_C(c_handle, jstate) &
      bind(C, name="mrsf_nac_metric_column")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    use oqp_tagarray_driver, only: tagarray_get_data, TA_TYPE_REAL64
    use, intrinsic :: iso_c_binding, only: c_int32_t

    type(oqp_handle_t) :: c_handle
    integer(c_int32_t), intent(in), value :: jstate
    type(information), pointer :: inf
    real(kind=dp), allocatable :: column(:,:)
    real(kind=dp), pointer :: exported(:,:)
    integer :: nbf, nstate

    inf => oqp_handle_get_info(c_handle)
    nbf = inf%basis%nbf
    nstate = inf%tddft%nstate
    allocate(column(nbf*nbf,nstate))
    call mrsf_nac_metric_column(inf, int(jstate), column)
    call inf%dat%remove_records((/ character(len=80) :: &
      'OQP::nac_gamma_column' /))
    call inf%dat%reserve_data('OQP::nac_gamma_column', TA_TYPE_REAL64, &
      nbf*nbf*nstate, (/ nbf*nbf, nstate /), &
      comment='diagnostic streamed exact-TLF metric column')
    call tagarray_get_data(inf%dat, 'OQP::nac_gamma_column', exported)
    exported = column
    deallocate(column)
  end subroutine mrsf_nac_metric_column_C

!###############################################################################

!> Build gamma^IJ_pq for every ordered off-diagonal state pair from the
!> resident MRSF amplitudes.  For the antisymmetric MO generator
!>
!>   K_pq = +1, K_qp = -1,
!>
!> the stored convention is
!>
!>   d S_IJ / d theta = sum_rs gamma^IJ_rs K_rs.
!>
!> Hence one half of the directional derivative is placed in each of the two
!> antisymmetric orbital slots.  State labels are deliberately not
!> antisymmetrized: the exact overlap is one-sided and column-normalized, so
!> gamma^IJ and -gamma^JI are not interchangeable identities.
  subroutine mrsf_nac_metric_data(infos)
    use types, only: information
    use oqp_tagarray_driver, only: data_has_tags, tagarray_get_data, &
      OQP_td_bvec_mo, TA_TYPE_REAL64
    use tdhf_mrsf_lib, only: mrsfxvec
    use messages, only: show_message, WITH_ABORT

    implicit none

    character(len=*), parameter :: subroutine_name = "mrsf_nac_metric_data"
    character(len=*), parameter :: tags_required(1) = (/ character(len=80) :: &
      OQP_td_bvec_mo /)

    type(information), target, intent(inout) :: infos
    real(kind=dp), contiguous, pointer :: bvec_mo(:,:), gamma_tlf(:,:,:)
    real(kind=dp), allocatable :: xvec(:,:), coeff(:,:,:), coeff_generic(:,:,:)
    real(kind=dp), allocatable :: sij(:,:), sab(:,:), sia(:,:)
    real(kind=dp), allocatable :: raw_overlap(:,:), column_norm(:)
    real(kind=dp), allocatable :: raw_sij(:,:,:,:), raw_sab(:,:,:,:)
    real(kind=dp), allocatable :: raw_sia(:,:,:,:)
    real(kind=dp), allocatable :: normalized_sij(:,:), normalized_sab(:,:)
    real(kind=dp), allocatable :: normalized_sia(:,:), mo_gradient(:,:)
    real(kind=dp) :: normalization_weight, half_derivative
    integer :: nbf, noca, nocb, nvirb, xvec_dim, nstate
    integer :: istate, jstate, kstate, p, q, ok

    nbf = infos%basis%nbf
    noca = infos%mol_prop%nelec_A
    nocb = infos%mol_prop%nelec_B
    nvirb = nbf - nocb
    xvec_dim = noca*nvirb
    nstate = infos%tddft%nstate

    if (noca - nocb /= 2) then
      call show_message( &
        'mrsf_nac_metric_data requires the two-SOMO MRSF reference.', &
        WITH_ABORT)
    end if
    if (nocb < 0) then
      call show_message('Invalid negative MRSF closed-shell count.', WITH_ABORT)
    end if
    if (infos%tddft%mult /= 1 .and. infos%tddft%mult /= 3) then
      call show_message( &
        'mrsf_nac_metric_data supports singlet or triplet MRSF states.', &
        WITH_ABORT)
    end if

    call data_has_tags(infos%dat, tags_required, module_name, &
                       subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)

    if (size(bvec_mo,1) < xvec_dim .or. size(bvec_mo,2) < nstate) then
      call show_message( &
        'OQP::td_bvec_mo has dimensions inconsistent with the MRSF space.', &
        WITH_ABORT)
    end if

    allocate(xvec(xvec_dim,nstate), &
             coeff(noca,nvirb,nstate), &
             coeff_generic(noca,nvirb,nstate), &
             sij(noca,noca), sab(nvirb,nvirb), sia(noca,nvirb), &
             raw_overlap(nstate,nstate), column_norm(nstate), &
             raw_sij(noca,noca,nstate,nstate), &
             raw_sab(nvirb,nvirb,nstate,nstate), &
             raw_sia(noca,nvirb,nstate,nstate), &
             normalized_sij(noca,noca), normalized_sab(nvirb,nvirb), &
             normalized_sia(noca,nvirb), mo_gradient(nbf,nbf), &
             source=0.0_dp, stat=ok)
    if (ok /= 0) then
      call show_message('Cannot allocate MRSF NAC metric workspace.', WITH_ABORT)
    end if

    ! Use the same resident amplitude unfolding as the production MRSF
    ! overlap driver, including the LR1/LR2 singlet/triplet pairing.
    do istate = 1, nstate
      call mrsfxvec(infos, bvec_mo(:,istate), xvec(:,istate))
      coeff(:,:,istate) = reshape(xvec(:,istate), (/ noca, nvirb /))
    end do

    ! The two-SOMO by two-SOMO coefficients occur only in the explicit MRSF
    ! blocks.  The generic determinant contractions omit them exactly as
    ! compute_states_overlap does.
    coeff_generic = coeff
    coeff_generic(nocb+1:noca,1:2,:) = 0.0_dp

    call build_identity_minors(noca, nocb, nbf, sij, sab, sia)
    call build_raw_overlap_sensitivities(coeff, coeff_generic, &
      sij, sab, sia, nocb, raw_overlap, raw_sij, raw_sab, raw_sia)

    do jstate = 1, nstate
      column_norm(jstate) = norm2(raw_overlap(:,jstate))
      if (column_norm(jstate) <= sqrt(tiny(1.0_dp))) then
        call show_message( &
          'Zero exact-overlap column in MRSF NAC metric construction.', &
          WITH_ABORT)
      end if
    end do

    call infos%dat%remove_records((/ character(len=80) :: tag_gamma /))
    call infos%dat%reserve_data(tag_gamma, TA_TYPE_REAL64, &
      nbf*nbf*nstate*nstate, (/ nbf*nbf, nstate, nstate /), &
      comment='closed-form exact-tlf MRSF state-overlap orbital source')
    call tagarray_get_data(infos%dat, tag_gamma, gamma_tlf)
    gamma_tlf = 0.0_dp

    ! Reverse the exact column normalization first.  For target S_IJ and raw
    ! overlap R_KJ, dS_IJ/dR_KJ is delta_IK/n_j-R_IJ*R_KJ/n_j^3.  The three
    ! resulting minor sensitivities are then contracted once with their
    ! sparse identity cofactors to obtain dS_IJ/dM.  This avoids an orbital-
    ! generator loop around the four-index contractions.
    do jstate = 1, nstate
      do istate = 1, nstate
        if (istate == jstate) cycle

        normalized_sij = 0.0_dp
        normalized_sab = 0.0_dp
        normalized_sia = 0.0_dp
        do kstate = 1, nstate
          normalization_weight = -raw_overlap(istate,jstate) &
            *raw_overlap(kstate,jstate)/(column_norm(jstate)**3)
          if (kstate == istate) then
            normalization_weight = normalization_weight &
              + 1.0_dp/column_norm(jstate)
          end if
          normalized_sij = normalized_sij + normalization_weight &
            *raw_sij(:,:,kstate,jstate)
          normalized_sab = normalized_sab + normalization_weight &
            *raw_sab(:,:,kstate,jstate)
          normalized_sia = normalized_sia + normalization_weight &
            *raw_sia(:,:,kstate,jstate)
        end do

        call accumulate_minor_cofactors(noca, nocb, nbf, normalized_sij, &
          normalized_sab, normalized_sia, mo_gradient)

        do q = 1, nbf - 1
          do p = q + 1, nbf
            half_derivative = 0.5_dp*(mo_gradient(p,q)-mo_gradient(q,p))
            gamma_tlf(p + (q-1)*nbf, istate, jstate) = half_derivative
            gamma_tlf(q + (p-1)*nbf, istate, jstate) = -half_derivative
          end do
        end do
      end do
    end do

    deallocate(xvec, coeff, coeff_generic, sij, sab, sia, raw_overlap, &
               column_norm, raw_sij, raw_sab, raw_sia, normalized_sij, &
               normalized_sab, normalized_sia, mo_gradient)
  end subroutine mrsf_nac_metric_data

!###############################################################################

!> Build one normalized exact-overlap metric column for the production NAC
!> driver.  Only R_KJ and dR_KJ (K=1..nstate) are resident at once, reducing
!> the metric workspace from O(nstate**2*nbf**2) to O(nstate*nbf**2).
!> gamma_column(:,I) contains gamma^(I,J); the diagonal column is zero.
  subroutine mrsf_nac_metric_column(infos, jstate, gamma_column)
    use types, only: information
    use oqp_tagarray_driver, only: data_has_tags, tagarray_get_data, &
      OQP_td_bvec_mo
    use tdhf_mrsf_lib, only: mrsfxvec
    use messages, only: show_message, WITH_ABORT

    implicit none

    character(len=*), parameter :: subroutine_name = &
      "mrsf_nac_metric_column"
    character(len=*), parameter :: tags_required(1) = (/ character(len=80) :: &
      OQP_td_bvec_mo /)

    type(information), target, intent(inout) :: infos
    integer, intent(in) :: jstate
    real(kind=dp), intent(out) :: gamma_column(:,:)
    real(kind=dp), contiguous, pointer :: bvec_mo(:,:)
    real(kind=dp), allocatable :: xvec(:,:), coeff(:,:,:), coeff_generic(:,:,:)
    real(kind=dp), allocatable :: sij(:,:), sab(:,:), sia(:,:)
    real(kind=dp), allocatable :: raw_overlap(:), raw_sij(:,:,:)
    real(kind=dp), allocatable :: raw_sab(:,:,:), raw_sia(:,:,:)
    real(kind=dp), allocatable :: normalized_sij(:,:), normalized_sab(:,:)
    real(kind=dp), allocatable :: normalized_sia(:,:), mo_gradient(:,:)
    real(kind=dp) :: column_norm, normalization_weight, half_derivative
    integer :: nbf, noca, nocb, nvirb, xvec_dim, nstate
    integer :: istate, kstate, p, q, ok

    nbf = infos%basis%nbf
    noca = infos%mol_prop%nelec_A
    nocb = infos%mol_prop%nelec_B
    nvirb = nbf - nocb
    xvec_dim = noca*nvirb
    nstate = infos%tddft%nstate

    if (jstate < 1 .or. jstate > nstate) then
      call show_message('Invalid state in MRSF NAC metric column.', WITH_ABORT)
    end if
    if (noca - nocb /= 2 .or. nocb < 0) then
      call show_message( &
        'MRSF NAC metric column requires a two-SOMO reference.', &
        WITH_ABORT)
    end if
    if (infos%tddft%mult /= 1 .and. infos%tddft%mult /= 3) then
      call show_message( &
        'MRSF NAC metric column supports singlet or triplet states.', &
        WITH_ABORT)
    end if
    if (size(gamma_column,1) /= nbf*nbf .or. &
        size(gamma_column,2) /= nstate) then
      call show_message('MRSF NAC metric column has inconsistent dimensions.', &
                        WITH_ABORT)
    end if

    call data_has_tags(infos%dat, tags_required, module_name, &
                       subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    if (size(bvec_mo,1) < xvec_dim .or. size(bvec_mo,2) < nstate) then
      call show_message( &
        'OQP::td_bvec_mo is inconsistent with the MRSF metric column.', &
        WITH_ABORT)
    end if

    allocate(xvec(xvec_dim,nstate), coeff(noca,nvirb,nstate), &
             coeff_generic(noca,nvirb,nstate), &
             sij(noca,noca), sab(nvirb,nvirb), sia(noca,nvirb), &
             raw_overlap(nstate), raw_sij(noca,noca,nstate), &
             raw_sab(nvirb,nvirb,nstate), &
             raw_sia(noca,nvirb,nstate), &
             normalized_sij(noca,noca), normalized_sab(nvirb,nvirb), &
             normalized_sia(noca,nvirb), mo_gradient(nbf,nbf), &
             source=0.0_dp, stat=ok)
    if (ok /= 0) then
      call show_message('Cannot allocate streamed MRSF metric workspace.', &
                        WITH_ABORT)
    end if

    do istate = 1, nstate
      call mrsfxvec(infos, bvec_mo(:,istate), xvec(:,istate))
      coeff(:,:,istate) = reshape(xvec(:,istate), (/ noca, nvirb /))
    end do
    coeff_generic = coeff
    coeff_generic(nocb+1:noca,1:2,:) = 0.0_dp
    call build_identity_minors(noca, nocb, nbf, sij, sab, sia)

    do kstate = 1, nstate
      call build_raw_overlap_pair_sensitivities(coeff, coeff_generic, &
        sij, sab, sia, nocb, kstate, jstate, raw_overlap(kstate), &
        raw_sij(:,:,kstate), raw_sab(:,:,kstate), raw_sia(:,:,kstate))
    end do
    column_norm = norm2(raw_overlap)
    if (column_norm <= sqrt(tiny(1.0_dp))) then
      call show_message('Zero exact-overlap column in streamed MRSF metric.', &
                        WITH_ABORT)
    end if

    gamma_column = 0.0_dp
    do istate = 1, nstate
      if (istate == jstate) cycle
      normalized_sij = 0.0_dp
      normalized_sab = 0.0_dp
      normalized_sia = 0.0_dp
      do kstate = 1, nstate
        normalization_weight = -raw_overlap(istate)*raw_overlap(kstate) &
          /(column_norm**3)
        if (kstate == istate) normalization_weight = &
          normalization_weight + 1.0_dp/column_norm
        normalized_sij = normalized_sij + normalization_weight &
          *raw_sij(:,:,kstate)
        normalized_sab = normalized_sab + normalization_weight &
          *raw_sab(:,:,kstate)
        normalized_sia = normalized_sia + normalization_weight &
          *raw_sia(:,:,kstate)
      end do

      call accumulate_minor_cofactors(noca, nocb, nbf, normalized_sij, &
        normalized_sab, normalized_sia, mo_gradient)
      do q = 1, nbf - 1
        do p = q + 1, nbf
          half_derivative = 0.5_dp*(mo_gradient(p,q)-mo_gradient(q,p))
          gamma_column(p + (q-1)*nbf,istate) = half_derivative
          gamma_column(q + (p-1)*nbf,istate) = -half_derivative
        end do
      end do
    end do

    deallocate(xvec, coeff, coeff_generic, sij, sab, sia, raw_overlap, &
      raw_sij, raw_sab, raw_sia, normalized_sij, normalized_sab, &
      normalized_sia, mo_gradient)
  end subroutine mrsf_nac_metric_column

!###############################################################################

!> Exact ndtlf=0 determinant minors at the identity MO overlap.
  subroutine build_identity_minors(noca, nocb, nbf, sij, sab, sia)
    implicit none

    integer, intent(in) :: noca, nocb, nbf
    real(kind=dp), intent(out) :: sij(:,:), sab(:,:), sia(:,:)

    integer :: i1, i2, j1, j2, noc, nvirb
    integer :: rows(noca-1), cols(noca-1)
    real(kind=dp) :: sign_factor

    noc = noca - 1
    nvirb = nbf - nocb

    do i2 = 1, noca
      do i1 = 1, noca
        call build_ij_maps(i1, i2, noca, rows, cols, sign_factor)
        sij(i1,i2) = sign_factor*identity_subdet(rows, cols, 0, 0)
      end do
    end do

    do j2 = 1, nvirb
      do j1 = 1, nvirb
        call build_ab_maps(j1, j2, nocb, rows, cols)
        sab(j1,j2) = identity_subdet(rows, cols, 0, 0)
      end do
    end do

    do j1 = 1, nvirb
      do i1 = 1, noca
        call build_ia_maps(i1, j1, noca, nocb, rows, cols)
        sia(i1,j1) = identity_subdet(rows, cols, 0, 0)
      end do
    end do

    if (noc /= size(rows)) error stop 'Invalid exact-minor work dimension.'
  end subroutine build_identity_minors

!###############################################################################

!> Contract analytic sensitivities to the three minor families with the
!> corresponding determinant cofactors at M=I.  Some exact minors are
!> singular at I, so Jacobi's inverse formula is invalid; direct cofactors are
!> accumulated instead.
  subroutine accumulate_minor_cofactors(noca, nocb, nbf, weight_sij, &
                                         weight_sab, weight_sia, gradient)
    implicit none

    integer, intent(in) :: noca, nocb, nbf
    real(kind=dp), intent(in) :: weight_sij(:,:), weight_sab(:,:)
    real(kind=dp), intent(in) :: weight_sia(:,:)
    real(kind=dp), intent(out) :: gradient(:,:)

    integer :: i1, i2, j1, j2, nvirb
    integer :: rows(noca-1), cols(noca-1)
    real(kind=dp) :: sign_factor

    nvirb = nbf - nocb
    gradient = 0.0_dp

    do i2 = 1, noca
      do i1 = 1, noca
        call build_ij_maps(i1, i2, noca, rows, cols, sign_factor)
        call accumulate_identity_cofactor(rows, cols, sign_factor, &
          weight_sij(i1,i2), gradient)
      end do
    end do

    do j2 = 1, nvirb
      do j1 = 1, nvirb
        call build_ab_maps(j1, j2, nocb, rows, cols)
        call accumulate_identity_cofactor(rows, cols, 1.0_dp, &
          weight_sab(j1,j2), gradient)
      end do
    end do

    do j1 = 1, nvirb
      do i1 = 1, noca
        call build_ia_maps(i1, j1, noca, nocb, rows, cols)
        call accumulate_identity_cofactor(rows, cols, 1.0_dp, &
          weight_sia(i1,j1), gradient)
      end do
    end do
  end subroutine accumulate_minor_cofactors

!###############################################################################

!> Row and column maps of ov_exact itype=1.
  pure subroutine build_ij_maps(i1, i2, noca, rows, cols, sign_factor)
    implicit none

    integer, intent(in) :: i1, i2, noca
    integer, intent(out) :: rows(:), cols(:)
    real(kind=dp), intent(out) :: sign_factor

    integer :: k, orbital, imin, imax

    k = 0
    if (i1 == i2) then
      do orbital = 1, noca
        if (orbital == i1) cycle
        k = k + 1
        rows(k) = orbital
        cols(k) = orbital
      end do
      sign_factor = 1.0_dp
    else
      imin = min(i1,i2)
      imax = max(i1,i2)
      do orbital = 1, noca
        if (orbital == imin .or. orbital == imax) cycle
        k = k + 1
        rows(k) = orbital
        cols(k) = orbital
      end do
      k = k + 1
      rows(k) = i2
      cols(k) = i1
      sign_factor = -1.0_dp
    end if
  end subroutine build_ij_maps

!###############################################################################

!> Row and column maps of ov_exact itype=2.
  pure subroutine build_ab_maps(j1, j2, nocb, rows, cols)
    implicit none

    integer, intent(in) :: j1, j2, nocb
    integer, intent(out) :: rows(:), cols(:)

    integer :: k

    do k = 1, nocb
      rows(k) = k
      cols(k) = k
    end do
    rows(nocb+1) = nocb + j1
    cols(nocb+1) = nocb + j2
  end subroutine build_ab_maps

!###############################################################################

!> Net row and column maps of the literal ov_exact itype=3 layout.
  pure subroutine build_ia_maps(i1, j1, noca, nocb, rows, cols)
    implicit none

    integer, intent(in) :: i1, j1, noca, nocb
    integer, intent(out) :: rows(:), cols(:)

    integer :: k, noc, orbital

    noc = noca - 1
    ! For a two-electron MRSF reference noc=1.  In the literal ov_exact
    ! case(3) write order, the sole determinant element is the final (4,4)
    ! assignment M(noc+1,nocb+j1); the preceding nominal (3,*) blocks address
    ! row zero and are an F77-era out-of-bounds artifact.  State the resulting
    ! 1x1 minor explicitly so the resident cofactor kernel has the same
    ! well-defined limiting semantics without invalid indices.
    if (noc == 1) then
      rows(1) = noc + 1
      cols(1) = nocb + j1
      return
    end if
    do k = 1, noc - 2
      orbital = k
      if (k > i1 - 1) orbital = k + 1
      rows(k) = orbital
      cols(k) = orbital
    end do
    rows(noc-1) = noc
    rows(noc) = noc + 1
    cols(noc-1) = i1
    cols(noc) = nocb + j1
  end subroutine build_ia_maps

!###############################################################################

!> Determinant of an identity submatrix after optional local row/column
!> deletion.  A submatrix of I is a partial permutation matrix, so its
!> determinant is exactly 0, +1, or -1 and requires no numerical factorization.
  pure real(kind=dp) function identity_subdet(rows, cols, skip_row, skip_col) &
      result(det_value)
    implicit none

    integer, intent(in) :: rows(:), cols(:), skip_row, skip_col

    integer :: permutation(size(rows))
    integer :: irow, icol, compact_row, compact_col
    integer :: found, inversions, nrow, ncol, i, j

    det_value = 0.0_dp
    nrow = size(rows)
    ncol = size(cols)
    if (skip_row > 0) nrow = nrow - 1
    if (skip_col > 0) ncol = ncol - 1
    if (nrow /= ncol) return
    if (nrow == 0) then
      det_value = 1.0_dp
      return
    end if

    compact_row = 0
    do irow = 1, size(rows)
      if (irow == skip_row) cycle
      compact_row = compact_row + 1
      compact_col = 0
      found = 0
      do icol = 1, size(cols)
        if (icol == skip_col) cycle
        compact_col = compact_col + 1
        if (rows(irow) == cols(icol)) then
          found = found + 1
          permutation(compact_row) = compact_col
        end if
      end do
      if (found /= 1) return
    end do

    inversions = 0
    do i = 1, nrow - 1
      do j = i + 1, nrow
        if (permutation(i) == permutation(j)) return
        if (permutation(i) > permutation(j)) inversions = inversions + 1
      end do
    end do

    det_value = 1.0_dp
    if (mod(inversions,2) /= 0) det_value = -1.0_dp
  end function identity_subdet

!###############################################################################

!> Add weight*d(det M(rows,cols))/dM to a full orbital gradient.  At the
!> identity, a minor with equal row/column label sets is a permutation matrix
!> and has one cofactor per row.  If the sets differ by one label, only the
!> unmatched row/column cofactor survives.  With two or more unmatched labels
!> every first cofactor is zero.  This is the sparse exact form of
!> d(det A)=sum_ab C_ab dA_ab.
  pure subroutine accumulate_identity_cofactor(rows, cols, sign_factor, &
                                                weight, gradient)
    implicit none

    integer, intent(in) :: rows(:), cols(:)
    real(kind=dp), intent(in) :: sign_factor, weight
    real(kind=dp), intent(inout) :: gradient(:,:)

    integer :: irow, icol, matched_col
    integer :: missing_rows, missing_cols, unmatched_row, unmatched_col
    real(kind=dp) :: cofactor, determinant_value, scale

    if (abs(weight) <= tiny(weight)) return

    missing_rows = 0
    unmatched_row = 0
    do irow = 1, size(rows)
      matched_col = 0
      do icol = 1, size(cols)
        if (rows(irow) == cols(icol)) matched_col = icol
      end do
      if (matched_col == 0) then
        missing_rows = missing_rows + 1
        unmatched_row = irow
      end if
    end do

    missing_cols = 0
    unmatched_col = 0
    do icol = 1, size(cols)
      matched_col = 0
      do irow = 1, size(rows)
        if (cols(icol) == rows(irow)) matched_col = irow
      end do
      if (matched_col == 0) then
        missing_cols = missing_cols + 1
        unmatched_col = icol
      end if
    end do

    if (missing_rows /= missing_cols .or. missing_rows > 1) return
    scale = sign_factor*weight

    if (missing_rows == 0) then
      determinant_value = identity_subdet(rows, cols, 0, 0)
      do irow = 1, size(rows)
        matched_col = 0
        do icol = 1, size(cols)
          if (rows(irow) == cols(icol)) matched_col = icol
        end do
        gradient(rows(irow),cols(matched_col)) = &
          gradient(rows(irow),cols(matched_col)) + scale*determinant_value
      end do
    else
      cofactor = identity_subdet(rows, cols, unmatched_row, unmatched_col)
      if (mod(unmatched_row+unmatched_col,2) /= 0) cofactor = -cofactor
      gradient(rows(unmatched_row),cols(unmatched_col)) = &
        gradient(rows(unmatched_row),cols(unmatched_col)) + scale*cofactor
    end if
  end subroutine accumulate_identity_cofactor

!###############################################################################

!> Value and analytic minor sensitivities of all seven exact-overlap
!> contractions before column normalization.  Each product contributes to
!> raw and to d(raw)/dS_ij, d(raw)/dS_ab, or d(raw)/dS_ia in the same loop.
  subroutine build_raw_overlap_sensitivities(coeff, coeff_generic, &
      sij, sab, sia, nocb, raw, raw_sij, raw_sab, raw_sia)
    implicit none

    real(kind=dp), intent(in) :: coeff(:,:,:), coeff_generic(:,:,:)
    real(kind=dp), intent(in) :: sij(:,:), sab(:,:), sia(:,:)
    integer, intent(in) :: nocb
    real(kind=dp), intent(out) :: raw(:,:)
    real(kind=dp), intent(out) :: raw_sij(:,:,:,:), raw_sab(:,:,:,:)
    real(kind=dp), intent(out) :: raw_sia(:,:,:,:)

    real(kind=dp), parameter :: rs = 1.0_dp/sqrt(2.0_dp)
    real(kind=dp) :: acc, left_coefficient, right_coefficient
    integer :: noca, nvirb, nbf, nstate
    integer :: oi, ni, i, j, a, b, pi, qi, ri, si

    noca = size(coeff,1)
    nvirb = size(coeff,2)
    nstate = size(coeff,3)
    nbf = nocb + nvirb
    raw = 0.0_dp
    raw_sij = 0.0_dp
    raw_sab = 0.0_dp
    raw_sia = 0.0_dp

    do ni = 1, nstate
      do oi = 1, nstate
        acc = 0.0_dp

        ! Block 1: SOMO--SOMO against SOMO--SOMO, S_ij*S_ab.
        do pi = nocb + 1, noca
          do qi = nocb + 1, noca
            do ri = nocb + 1, noca
              do si = nocb + 1, noca
                left_coefficient = coeff(pi,qi-nocb,oi)
                right_coefficient = coeff(ri,si-nocb,ni)
                acc = acc + left_coefficient*right_coefficient &
                  *sij(pi,ri)*sab(qi-nocb,si-nocb)
                raw_sij(pi,ri,oi,ni) = raw_sij(pi,ri,oi,ni) &
                  + left_coefficient*right_coefficient &
                  *sab(qi-nocb,si-nocb)
                raw_sab(qi-nocb,si-nocb,oi,ni) = &
                  raw_sab(qi-nocb,si-nocb,oi,ni) &
                  + left_coefficient*right_coefficient*sij(pi,ri)
              end do
            end do
          end do
        end do

        ! Blocks 2 and 3: left SOMO--SOMO against right generic.
        do pi = nocb + 1, noca
          do qi = nocb + 1, noca
            do ri = 1, noca
              do si = nocb + 1, nbf
                if (ri >= nocb+1 .and. si <= noca) cycle
                left_coefficient = coeff(pi,qi-nocb,oi)
                right_coefficient = coeff(ri,si-nocb,ni)
                acc = acc + rs*left_coefficient*right_coefficient*( &
                  sij(pi,ri)*sab(qi-nocb,si-nocb) &
                  + sia(pi,si-nocb)*sia(ri,qi-nocb))
                raw_sij(pi,ri,oi,ni) = raw_sij(pi,ri,oi,ni) &
                  + rs*left_coefficient*right_coefficient &
                  *sab(qi-nocb,si-nocb)
                raw_sab(qi-nocb,si-nocb,oi,ni) = &
                  raw_sab(qi-nocb,si-nocb,oi,ni) &
                  + rs*left_coefficient*right_coefficient*sij(pi,ri)
                raw_sia(pi,si-nocb,oi,ni) = raw_sia(pi,si-nocb,oi,ni) &
                  + rs*left_coefficient*right_coefficient &
                  *sia(ri,qi-nocb)
                raw_sia(ri,qi-nocb,oi,ni) = raw_sia(ri,qi-nocb,oi,ni) &
                  + rs*left_coefficient*right_coefficient &
                  *sia(pi,si-nocb)
              end do
            end do
          end do
        end do

        ! Blocks 4 and 5: left generic against right SOMO--SOMO.
        do pi = 1, noca
          do qi = nocb + 1, nbf
            if (pi >= nocb+1 .and. qi <= noca) cycle
            do ri = nocb + 1, noca
              do si = nocb + 1, noca
                left_coefficient = coeff(pi,qi-nocb,oi)
                right_coefficient = coeff(ri,si-nocb,ni)
                acc = acc + rs*left_coefficient*right_coefficient*( &
                  sij(pi,ri)*sab(qi-nocb,si-nocb) &
                  + sia(pi,si-nocb)*sia(ri,qi-nocb))
                raw_sij(pi,ri,oi,ni) = raw_sij(pi,ri,oi,ni) &
                  + rs*left_coefficient*right_coefficient &
                  *sab(qi-nocb,si-nocb)
                raw_sab(qi-nocb,si-nocb,oi,ni) = &
                  raw_sab(qi-nocb,si-nocb,oi,ni) &
                  + rs*left_coefficient*right_coefficient*sij(pi,ri)
                raw_sia(pi,si-nocb,oi,ni) = raw_sia(pi,si-nocb,oi,ni) &
                  + rs*left_coefficient*right_coefficient &
                  *sia(ri,qi-nocb)
                raw_sia(ri,qi-nocb,oi,ni) = raw_sia(ri,qi-nocb,oi,ni) &
                  + rs*left_coefficient*right_coefficient &
                  *sia(pi,si-nocb)
              end do
            end do
          end do
        end do

        ! Block 6: generic amplitudes coupled by S_ij*S_ab.
        do i = 1, noca
          do a = 1, nvirb
            do j = 1, noca
              do b = 1, nvirb
                left_coefficient = coeff_generic(i,a,oi)
                right_coefficient = coeff_generic(j,b,ni)
                acc = acc + left_coefficient*right_coefficient &
                  *sij(i,j)*sab(a,b)
                raw_sij(i,j,oi,ni) = raw_sij(i,j,oi,ni) &
                  + left_coefficient*right_coefficient*sab(a,b)
                raw_sab(a,b,oi,ni) = raw_sab(a,b,oi,ni) &
                  + left_coefficient*right_coefficient*sij(i,j)
              end do
            end do
          end do
        end do

        ! Block 7: generic amplitudes coupled by the two S_ia factors.
        do i = 1, noca
          do a = 1, nvirb
            do j = 1, noca
              do b = 1, nvirb
                left_coefficient = coeff_generic(i,a,oi)
                right_coefficient = coeff_generic(j,b,ni)
                acc = acc + left_coefficient*right_coefficient &
                  *sia(j,a)*sia(i,b)
                raw_sia(j,a,oi,ni) = raw_sia(j,a,oi,ni) &
                  + left_coefficient*right_coefficient*sia(i,b)
                raw_sia(i,b,oi,ni) = raw_sia(i,b,oi,ni) &
                  + left_coefficient*right_coefficient*sia(j,a)
              end do
            end do
          end do
        end do

        raw(oi,ni) = acc
      end do
    end do
  end subroutine build_raw_overlap_sensitivities

!###############################################################################

!> One (old-state,new-state) member of the seven-term exact-overlap value and
!> minor sensitivities.  This is the streamed production counterpart of
!> build_raw_overlap_sensitivities; keeping the algebra in a pair-sized output
!> avoids four-dimensional nstate**2 work arrays in the production driver.
  subroutine build_raw_overlap_pair_sensitivities(coeff, coeff_generic, &
      sij, sab, sia, nocb, oi, ni, raw, raw_sij, raw_sab, raw_sia)
    implicit none

    real(kind=dp), intent(in) :: coeff(:,:,:), coeff_generic(:,:,:)
    real(kind=dp), intent(in) :: sij(:,:), sab(:,:), sia(:,:)
    integer, intent(in) :: nocb, oi, ni
    real(kind=dp), intent(out) :: raw
    real(kind=dp), intent(out) :: raw_sij(:,:), raw_sab(:,:), raw_sia(:,:)

    real(kind=dp), parameter :: rs = 1.0_dp/sqrt(2.0_dp)
    real(kind=dp) :: left_coefficient, right_coefficient
    integer :: noca, nvirb, nbf
    integer :: i, j, a, b, pi, qi, ri, si

    noca = size(coeff,1)
    nvirb = size(coeff,2)
    nbf = nocb + nvirb
    raw = 0.0_dp
    raw_sij = 0.0_dp
    raw_sab = 0.0_dp
    raw_sia = 0.0_dp

    ! Block 1: SOMO--SOMO against SOMO--SOMO, S_ij*S_ab.
    do pi = nocb + 1, noca
      do qi = nocb + 1, noca
        do ri = nocb + 1, noca
          do si = nocb + 1, noca
            left_coefficient = coeff(pi,qi-nocb,oi)
            right_coefficient = coeff(ri,si-nocb,ni)
            raw = raw + left_coefficient*right_coefficient &
              *sij(pi,ri)*sab(qi-nocb,si-nocb)
            raw_sij(pi,ri) = raw_sij(pi,ri) &
              + left_coefficient*right_coefficient*sab(qi-nocb,si-nocb)
            raw_sab(qi-nocb,si-nocb) = raw_sab(qi-nocb,si-nocb) &
              + left_coefficient*right_coefficient*sij(pi,ri)
          end do
        end do
      end do
    end do

    ! Blocks 2 and 3: left SOMO--SOMO against right generic.
    do pi = nocb + 1, noca
      do qi = nocb + 1, noca
        do ri = 1, noca
          do si = nocb + 1, nbf
            if (ri >= nocb+1 .and. si <= noca) cycle
            left_coefficient = coeff(pi,qi-nocb,oi)
            right_coefficient = coeff(ri,si-nocb,ni)
            raw = raw + rs*left_coefficient*right_coefficient*( &
              sij(pi,ri)*sab(qi-nocb,si-nocb) &
              + sia(pi,si-nocb)*sia(ri,qi-nocb))
            raw_sij(pi,ri) = raw_sij(pi,ri) &
              + rs*left_coefficient*right_coefficient &
              *sab(qi-nocb,si-nocb)
            raw_sab(qi-nocb,si-nocb) = raw_sab(qi-nocb,si-nocb) &
              + rs*left_coefficient*right_coefficient*sij(pi,ri)
            raw_sia(pi,si-nocb) = raw_sia(pi,si-nocb) &
              + rs*left_coefficient*right_coefficient*sia(ri,qi-nocb)
            raw_sia(ri,qi-nocb) = raw_sia(ri,qi-nocb) &
              + rs*left_coefficient*right_coefficient*sia(pi,si-nocb)
          end do
        end do
      end do
    end do

    ! Blocks 4 and 5: left generic against right SOMO--SOMO.
    do pi = 1, noca
      do qi = nocb + 1, nbf
        if (pi >= nocb+1 .and. qi <= noca) cycle
        do ri = nocb + 1, noca
          do si = nocb + 1, noca
            left_coefficient = coeff(pi,qi-nocb,oi)
            right_coefficient = coeff(ri,si-nocb,ni)
            raw = raw + rs*left_coefficient*right_coefficient*( &
              sij(pi,ri)*sab(qi-nocb,si-nocb) &
              + sia(pi,si-nocb)*sia(ri,qi-nocb))
            raw_sij(pi,ri) = raw_sij(pi,ri) &
              + rs*left_coefficient*right_coefficient &
              *sab(qi-nocb,si-nocb)
            raw_sab(qi-nocb,si-nocb) = raw_sab(qi-nocb,si-nocb) &
              + rs*left_coefficient*right_coefficient*sij(pi,ri)
            raw_sia(pi,si-nocb) = raw_sia(pi,si-nocb) &
              + rs*left_coefficient*right_coefficient*sia(ri,qi-nocb)
            raw_sia(ri,qi-nocb) = raw_sia(ri,qi-nocb) &
              + rs*left_coefficient*right_coefficient*sia(pi,si-nocb)
          end do
        end do
      end do
    end do

    ! Block 6: generic amplitudes coupled by S_ij*S_ab.
    do i = 1, noca
      do a = 1, nvirb
        do j = 1, noca
          do b = 1, nvirb
            left_coefficient = coeff_generic(i,a,oi)
            right_coefficient = coeff_generic(j,b,ni)
            raw = raw + left_coefficient*right_coefficient*sij(i,j)*sab(a,b)
            raw_sij(i,j) = raw_sij(i,j) &
              + left_coefficient*right_coefficient*sab(a,b)
            raw_sab(a,b) = raw_sab(a,b) &
              + left_coefficient*right_coefficient*sij(i,j)
          end do
        end do
      end do
    end do

    ! Block 7: generic amplitudes coupled by the two S_ia factors.
    do i = 1, noca
      do a = 1, nvirb
        do j = 1, noca
          do b = 1, nvirb
            left_coefficient = coeff_generic(i,a,oi)
            right_coefficient = coeff_generic(j,b,ni)
            raw = raw + left_coefficient*right_coefficient*sia(j,a)*sia(i,b)
            raw_sia(j,a) = raw_sia(j,a) &
              + left_coefficient*right_coefficient*sia(i,b)
            raw_sia(i,b) = raw_sia(i,b) &
              + left_coefficient*right_coefficient*sia(j,a)
          end do
        end do
      end do
    end do
  end subroutine build_raw_overlap_pair_sensitivities

end module mrsf_nac_metric_data_mod
