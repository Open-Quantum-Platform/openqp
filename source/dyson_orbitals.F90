! =====================================================================
!  dyson_orbitals.F90
!  OpenQP: Dyson orbitals via EKT (Eq. 5) with TagArray I/O
!  - Provides two modules:
!      * oqp_ekt_kernel        : small helpers for EKT assembly
!      * dyson_orbitals_mod    : public driver, compute, store, print
!  - TagArray usage follows OpenQP wiki:
!        remove_records -> reserve_data -> tagarray_get_data -> write
! =====================================================================
module dyson_orbitals_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t
  use precision,      only: dp
  use io_constants,   only: iw
  use types,          only: information
  use printing,       only: print_module_info
  use eigen,          only: diag_symm_full
  use blas_wrap,      only: oqp_dgemm_i64

  ! >>> IMPORTANT: bring 'container_t' into scope so generics resolve <<<
  use tagarray                                   ! type(container_t), container%methods

  ! TagArray driver: tag names, TA_TYPE_*, and helpers
  use oqp_tagarray_driver, only: tagarray_get_data, data_has_tags, &
                                 TA_TYPE_REAL64, TA_TYPE_INT64, &
                                 OQP_VEC_DO, OQP_E_DO, OQP_DO_Strength, OQP_DO_Type, OQP_DO_Count, &
                                 OQP_VEC_DO_comment, OQP_E_DO_comment, OQP_DO_Strength_comment, OQP_DO_Type_comment, OQP_DO_Count_comment

  implicit none
  private

  public :: dyson_orbital_driver
  public :: compute_dyson_ekt
  public :: store_dyson_results_tagarray
  public :: print_dyson_results

  real(dp), parameter :: HARTREE_TO_EV = 27.21139_dp

contains

  !--------------------------------------------------------------------
  ! Driver: given relaxed density D and Lagrangian W, run EKT and export
  !--------------------------------------------------------------------
  subroutine dyson_orbital_driver(infos, D, W)
    type(information), intent(inout) :: infos
    real(dp), intent(in) :: D(:,:), W(:,:)

    integer :: nbf, nroots, j
    real(dp), allocatable :: dyson_AO(:,:), energies(:), strengths(:)
    real(dp), allocatable :: dyson_pad(:,:,:), e_pad(:,:), s_pad(:,:)

    call print_module_info('Dyson Orbitals (EKT) module', 'Dyson Orbitals (EKT) module')

    nbf = size(D,1)
    if (size(D,2) /= nbf .or. size(W,1) /= nbf .or. size(W,2) /= nbf) then
      write(iw,'(5x,"[EKT] Dimension mismatch in D/W.")')
      return
    end if

    call compute_dyson_ekt(infos, D, W, dyson_AO, energies, strengths, nroots)

    ! Pad to (nbf x nbf x 1) for downstream convention
    allocate(dyson_pad(nbf,nbf,1)); dyson_pad = 0.0_dp
    allocate(e_pad(nbf,1));         e_pad      = 0.0_dp
    allocate(s_pad(nbf,1));         s_pad      = 0.0_dp
    do j = 1, min(nbf, nroots)
      dyson_pad(:,j,1) = dyson_AO(:,j)
      e_pad(j,1)       = energies(j)
      s_pad(j,1)       = strengths(j)
    end do

    call store_dyson_results_tagarray(infos, dyson_pad, e_pad, s_pad, nbf, 1, 1)
    call print_dyson_results(infos)
  end subroutine dyson_orbital_driver


  !--------------------------------------------------------------------
  ! Core EKT compute: implements Eq.(5)
  !--------------------------------------------------------------------
  subroutine compute_dyson_ekt(infos, D_in, W_in, dyson_orbs_AO, binding_energies, pole_strengths, nroots)
    type(information), intent(in) :: infos
    real(dp), intent(in) :: D_in(:,:), W_in(:,:)

    real(dp), allocatable, intent(out) :: dyson_orbs_AO(:,:)      ! (nbf x nroots)
    real(dp), allocatable, intent(out) :: binding_energies(:)     ! (nroots)
    real(dp), allocatable, intent(out) :: pole_strengths(:)       ! (nroots)
    integer, intent(out) :: nroots

    integer :: nbf, info, i, j, k, nkeep
    real(dp) :: def_tol, tiny, normj
    real(dp), allocatable :: D(:,:), W(:,:), U(:,:), occ(:)
    real(dp), allocatable :: Ukeep(:,:), occ_keep(:)
    real(dp), allocatable :: M(:,:), Y(:,:), a_NO(:,:)

    nbf = size(D_in,1)
    tiny = 1.0e-16_dp
    def_tol = max(infos%tddft%dyson_deflation_tol, 1.0e-12_dp)

    allocate(D(nbf,nbf)); D = D_in
    allocate(W(nbf,nbf)); W = W_in
    call symmetrize_inplace(nbf, D)
    call symmetrize_inplace(nbf, W)

    ! Diagonalize D: eigenvectors in U (columns), eigenvalues in occ
    allocate(U(nbf,nbf)); U = D          ! copy D into U; diag overwrites U with eigenvectors
    allocate(occ(nbf));   occ = 0.0_dp
    info = 0
    call diag_symm_full(0, nbf, U, nbf, occ, info)

    ! Sort NOs by descending occupation and reorder U accordingly
    call sort_descending_with_vectors(nbf, occ, U)

    ! Deflation: keep NOs with occ > def_tol
    nkeep = count(occ > def_tol)
    if (nkeep <= 0) then
      write(iw,'(5x,"[EKT] No natural orbitals above deflation threshold.")')
      allocate(dyson_orbs_AO(0,0), binding_energies(0), pole_strengths(0))
      nroots = 0
      deallocate(D,W,U,occ)
      return
    end if

    allocate(Ukeep(nbf,nkeep)); Ukeep = U(:,1:nkeep)
    allocate(occ_keep(nkeep));  occ_keep = occ(1:nkeep)

    ! M = D^{-1/2} (U^T W U) D^{-1/2} in NO space
    allocate(M(nkeep,nkeep)); M = 0.0_dp
    call form_ekt_M_in_NO(nbf, nkeep, Ukeep, occ_keep, W, M)

    ! Diagonalize M: eigenvectors in Y (columns), eigenvalues in binding_energies
    allocate(Y(nkeep,nkeep)); Y = M
    allocate(binding_energies(nkeep)); binding_energies = 0.0_dp
    info = 0
    call diag_symm_full(0, nkeep, Y, nkeep, binding_energies, info)

    ! Dyson amplitudes in NO basis: a = D^{1/2} * Y
    allocate(a_NO(nkeep,nkeep)); a_NO = 0.0_dp
    do j = 1, nkeep
      do i = 1, nkeep
        a_NO(i,j) = sqrt(max(occ_keep(i), tiny)) * Y(i,j)
      end do
    end do

    ! Transform to AO: dyson_AO = Ukeep * a_NO
    allocate(dyson_orbs_AO(nbf, nkeep)); dyson_orbs_AO = 0.0_dp
    call oqp_dgemm_i64('N','N', nbf, nkeep, nkeep, 1.0_dp, Ukeep, nbf, a_NO, nkeep, 0.0_dp, dyson_orbs_AO, nbf)

    ! Pole strengths = ||a||^2 ; normalize AO vectors by sqrt(strength)
    allocate(pole_strengths(nkeep)); pole_strengths = 0.0_dp
    do j = 1, nkeep
      normj = 0.0_dp
      do i = 1, nkeep
        normj = normj + a_NO(i,j)*a_NO(i,j)
      end do
      pole_strengths(j) = normj
      if (normj > tiny) then
        normj = 1.0_dp / sqrt(normj)
        do k = 1, nbf
          dyson_orbs_AO(k,j) = dyson_orbs_AO(k,j) * normj
        end do
      end if
    end do

    nroots = nkeep
    deallocate(D,W,U,occ,Ukeep,occ_keep,M,Y,a_NO)
  end subroutine compute_dyson_ekt


  !--------------------------------------------------------------------
  ! TagArray store: follows wiki pattern
  !--------------------------------------------------------------------
  subroutine store_dyson_results_tagarray(infos, dyson_orbs, binding_energies, &
                                          pole_strengths, nbf, nstates, target_state)
    use messages, only: show_message, WITHOUT_ABORT
    type(information), intent(inout) :: infos
    real(dp), intent(in) :: dyson_orbs(nbf, nbf, nstates)
    real(dp), intent(in) :: binding_energies(nbf, nstates)
    real(dp), intent(in) :: pole_strengths(nbf, nstates)
    integer,  intent(in) :: nbf, nstates, target_state

    real(dp), allocatable :: vec_do(:,:), e_do(:), do_strength(:)
    integer(8), allocatable :: do_type(:)
    integer :: i, idx, n_sig
    real(dp) :: thr
    character(len=80) :: tags_to_remove(5)

    ! Pointers to TagArray storage
    real(8),    pointer :: p_vec_do(:,:)
    real(8),    pointer :: p_e_do(:)
    real(8),    pointer :: p_do_strength(:)
    integer(8), pointer :: p_do_type(:)
    integer(8), pointer :: p_count(:)

    thr   = infos%tddft%dyson_pole_threshold
    n_sig = count(pole_strengths(:,target_state) > thr)

    ! Remove any previous records for our tags (equal-length constructor avoids warnings)
    tags_to_remove = [ character(len=80) :: OQP_DO_Count, OQP_VEC_DO, OQP_E_DO, OQP_DO_Strength, OQP_DO_Type ]
    call infos%dat%remove_records(tags_to_remove)

    ! DO_Count: 1-element INT64
    call infos%dat%reserve_data(OQP_DO_Count, TA_TYPE_INT64, 1, (/1/), comment=OQP_DO_Count_comment)
    call tagarray_get_data(infos%dat, OQP_DO_Count, p_count)
    p_count(1) = int(n_sig, 8)

    if (n_sig <= 0) then
      write(iw,'(5x,"[Dyson] No orbitals above pole-strength threshold.")')
      return
    end if

    ! Compact arrays (nbf x n_sig) and (n_sig)
    allocate(vec_do(nbf, n_sig), e_do(n_sig), do_strength(n_sig), do_type(n_sig))
    idx = 0
    do i = 1, nbf
      if (pole_strengths(i, target_state) > thr) then
        idx = idx + 1
        vec_do(:, idx)   = dyson_orbs(:, i, target_state)
        e_do(idx)        = binding_energies(i, target_state)
        do_strength(idx) = pole_strengths(i, target_state)
        do_type(idx)     = merge(1_8, -1_8, e_do(idx) > 0.0_dp)   ! +1=IP, -1=EA
      end if
    end do

    ! Reserve TagArray buffers and map pointers
    call infos%dat%reserve_data(OQP_VEC_DO,      TA_TYPE_REAL64, nbf*n_sig, (/ nbf, n_sig /), comment=OQP_VEC_DO_comment)
    call infos%dat%reserve_data(OQP_E_DO,        TA_TYPE_REAL64, n_sig,     (/ n_sig /),      comment=OQP_E_DO_comment)
    call infos%dat%reserve_data(OQP_DO_Strength, TA_TYPE_REAL64, n_sig,     (/ n_sig /),      comment=OQP_DO_Strength_comment)
    call infos%dat%reserve_data(OQP_DO_Type,     TA_TYPE_INT64,  n_sig,     (/ n_sig /),      comment=OQP_DO_Type_comment)

    call tagarray_get_data(infos%dat, OQP_VEC_DO,      p_vec_do)
    call tagarray_get_data(infos%dat, OQP_E_DO,        p_e_do)
    call tagarray_get_data(infos%dat, OQP_DO_Strength, p_do_strength)
    call tagarray_get_data(infos%dat, OQP_DO_Type,     p_do_type)

    ! Copy data
    p_vec_do(:,:)    = vec_do
    p_e_do(:)        = e_do
    p_do_strength(:) = do_strength
    p_do_type(:)     = do_type

    deallocate(vec_do, e_do, do_strength, do_type)
  end subroutine store_dyson_results_tagarray


  !--------------------------------------------------------------------
  ! Pretty-print results from TagArray
  !--------------------------------------------------------------------
  subroutine print_dyson_results(infos)
    use messages, only: show_message, WITH_ABORT
    type(information), intent(in) :: infos

    real(8),    pointer :: vec_do(:,:)
    real(8),    pointer :: e_do(:)
    real(8),    pointer :: do_strength(:)
    integer(8), pointer :: do_type(:)
    integer(8), pointer :: p_count(:)
    integer :: n, i, n_print
    character(len=80) :: tags(5)

    ! Ensure tags exist (will print a message if any is missing)
    tags = [ character(len=80) :: OQP_DO_Count, OQP_VEC_DO, OQP_E_DO, OQP_DO_Strength, OQP_DO_Type ]
    call data_has_tags(infos%dat, tags, 'dyson_orbitals_mod', 'print_dyson_results', WITH_ABORT)

    call tagarray_get_data(infos%dat, OQP_DO_Count, p_count)
    n = int(p_count(1))

    call tagarray_get_data(infos%dat, OQP_VEC_DO,      vec_do)
    call tagarray_get_data(infos%dat, OQP_E_DO,        e_do)
    call tagarray_get_data(infos%dat, OQP_DO_Strength, do_strength)
    call tagarray_get_data(infos%dat, OQP_DO_Type,     do_type)

    n_print = min(n, 20)
    write(iw,'(/,5x,"Dyson results (top ",i0,"):",/)') n_print
    do i = 1, n_print
      if (do_type(i) == 1_8) then
        write(iw,'(5x,"#",i3,":  IP  = ",f12.6," Ha  (",f10.4," eV)   |  Strength = ",f8.5)') &
             i, e_do(i),  e_do(i)*HARTREE_TO_EV, do_strength(i)
      else
        write(iw,'(5x,"#",i3,":  EA  = ",f12.6," Ha  (",f10.4," eV)   |  Strength = ",f8.5)') &
             i, e_do(i), -e_do(i)*HARTREE_TO_EV, do_strength(i)
      end if
    end do
    write(iw,'( )')
  end subroutine print_dyson_results


  !--------------------------------------------------------------------
  ! Utility: sort eigenvalues (desc) and reorder eigenvectors accordingly
  !--------------------------------------------------------------------
  subroutine sort_descending_with_vectors(n, eval, V)
    integer, intent(in) :: n
    real(dp), intent(inout) :: eval(n)
    real(dp), intent(inout) :: V(n,n)
    integer :: i, j, kmax
    real(dp) :: vmax, tmp
    real(dp), allocatable :: tmpcol(:)

    allocate(tmpcol(n))
    do i = 1, n-1
      kmax = i
      vmax = eval(i)
      do j = i+1, n
        if (eval(j) > vmax) then
          vmax = eval(j)
          kmax = j
        end if
      end do
      if (kmax /= i) then
        tmp        = eval(i)
        eval(i)    = eval(kmax)
        eval(kmax) = tmp
        tmpcol(:)    = V(:,i)
        V(:,i)       = V(:,kmax)
        V(:,kmax)    = tmpcol(:)
      end if
    end do
    deallocate(tmpcol)
  end subroutine sort_descending_with_vectors

  subroutine form_ekt_M_in_NO(nbf, nkeep, Ukeep, occ_keep, W, M)
    ! Build M = D^{-1/2} (U^T W U) D^{-1/2} in the NO subspace
    integer, intent(in) :: nbf, nkeep
    real(dp), intent(in) :: Ukeep(nbf,nkeep), occ_keep(nkeep), W(nbf,nbf)
    real(dp), intent(out) :: M(nkeep,nkeep)
    real(dp), allocatable :: T(:,:), A(:,:)
    integer :: i, j
    real(dp) :: di, dj, tiny

    allocate(T(nbf,nkeep)); T = 0.0_dp
    allocate(A(nkeep,nkeep)); A = 0.0_dp
    M = 0.0_dp

    ! T = W * Ukeep
    call oqp_dgemm_i64('N','N', nbf, nkeep, nbf, 1.0_dp, W, nbf, Ukeep, nbf, 0.0_dp, T, nbf)
    ! A = Ukeep^T * T
    call oqp_dgemm_i64('T','N', nkeep, nkeep, nbf, 1.0_dp, Ukeep, nbf, T, nbf, 0.0_dp, A, nkeep)

    tiny = 1.0e-16_dp
    do j = 1, nkeep
      dj = sqrt(max(occ_keep(j), tiny))
      do i = 1, nkeep
        di = sqrt(max(occ_keep(i), tiny))
        M(i,j) = A(i,j) / (di*dj)
      end do
    end do

    deallocate(T, A)
  end subroutine form_ekt_M_in_NO

end module dyson_orbitals_mod
