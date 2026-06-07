!> @file dk_scalar.F90
!> @brief Douglas-Kroll scalar relativistic correction to the one-electron Hamiltonian
!>
!> @details Implements the DK1 and DK2 Douglas-Kroll-Hess transformations that replace
!>          the non-relativistic H_core = T + V with the scalar relativistic H^DK.
!>          The DK transformation is carried out in the momentum (p) representation
!>          obtained by diagonalising the kinetic energy matrix T.
!>
!>          Pipeline (called once per SCF):
!>            1. Compute pVp integrals  <mu|p(-sum_A Z_A/r_A)p|nu>
!>            2. Build p-space basis:   S^{-1/2} -> XU, SXU, p^2 eigenvalues
!>            3. Compute kinematic factors: E_p, A, R
!>            4. Transform V and pVp to p-space
!>            5. Build H^DK1 in p-space (DK1 correction)
!>            6. Add H^DK2 correction in p-space (DK2 correction)
!>            7. Back-transform to AO basis -> overwrite OQP::Hcore
!>
!> @author Vladimir Makhnev
!> @date   March 2026

module dk_scalar_mod

  implicit none

  character(len=*), parameter :: module_name = "dk_scalar_mod"

  !> Set to .true. to enable diagnostic output from DK routines
  logical :: dk_debug = .false.

  private compute_and_check_pvp
  public dk_scalar

contains

  !> @brief C-interop wrapper: unpack the OQP handle and call dk_scalar
  subroutine dk_scalar_C(c_handle) bind(C, name="dk_scalar")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call dk_scalar(inf)
  end subroutine dk_scalar_C

  !> @brief Apply scalar relativistic Douglas-Kroll correction to H_core
  !>
  !> @details Reads OQP::SM (overlap), OQP::TM (kinetic energy), and
  !>          OQP::Hcore (= T + V) from the tagarray, performs the DK1+DK2
  !>          transformation, and overwrites OQP::Hcore with H^DK.
  !>
  !> @param[inout] infos  OQP information struct (basis, atoms, tagarray, log)
  subroutine dk_scalar(infos)

    use types,               only: information
    use oqp_tagarray_driver
    use precision,           only: dp
    use io_constants,        only: iw
    use basis_tools,         only: basis_set
    use messages,            only: show_message, WITH_ABORT
    use printing,            only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "dk_scalar"

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis

    integer :: nbf    ! number of AO basis functions
    integer :: nbf2   ! triangular size nbf*(nbf+1)/2
    integer :: ok     ! allocation status
    integer :: i, j, idx

    ! --- tagarray pointers (no copy: point into tagarray storage) ---
    real(kind=dp), contiguous, pointer :: hcore(:), tmat(:), smat(:)

    character(len=*), parameter :: tags_required(3) = (/ character(len=80) :: &
      OQP_SM, OQP_TM, OQP_Hcore /)

    ! --- working arrays ---
    real(kind=dp), allocatable :: &
      pvp(:),    &  ! <mu|p V p|nu>,  packed triangular (nbf2)
      XU(:,:),   &  ! X*U:  columns are the p-space basis vectors in AO rep.
      SXU(:,:),  &  ! S*X*U: used for the back-transformation to AO basis
      psq(:),    &  ! p_i^2: eigenvalues of 2T in the orthonormal basis
      Ep(:),     &  ! relativistic kinetic energy  E_p = c*sqrt(p^2 + c^2)
      Akin(:),   &  ! kinematic factor  A_i = sqrt((E_p+c^2)/(2*E_p))
      Rkin(:),   &  ! kinematic factor  R_i = c/(E_p+c^2)
      hdk(:)        ! H^DK in AO basis, packed triangular (nbf2)
    real(kind=dp), allocatable :: hdkp(:)  ! H^DK in p-space, packed triangular

    real(kind=dp), allocatable :: &
        Vp(:),    &  ! V = Hcore-T transformed to p-space, packed triangular
        PVPp(:)      ! pVp transformed to p-space, packed triangular
    integer :: qrnk  ! effective rank after removing linear dependencies in S

    dk_debug = (infos%control%verbose > 1)

    open(unit=iw, file=infos%log_filename, position="append")

    call print_module_info('DK_SCALAR', 'Douglas-Kroll Scalar Relativistic Correction')

    ! --- check requested DK order ---
    select case (infos%control%scal_rel)
    case (0)
      write(iw, '(1x,a)') 'scal_rel = 0: scalar relativistic correction disabled, skipping.'
      close(iw)
      return
    case (1)
      write(iw, '(1x,a)') 'scal_rel = 1: applying DK1 correction.'
    case (2)
      write(iw, '(1x,a)') 'scal_rel = 2: applying DK1 + DK2 correction.'
    case default
      write(iw, '(1x,a,i0)') 'WARNING: unknown scal_rel value: ', infos%control%scal_rel
      write(iw, '(1x,a)') 'Defaulting to DK2.'
    end select

    basis => infos%basis
    basis%atoms => infos%atoms

    nbf  = basis%nbf
    nbf2 = nbf*(nbf+1)/2

    ! --- verify required tags are present ---
    call data_has_tags(infos%dat, tags_required, &
                       module_name, subroutine_name, WITH_ABORT)

    call tagarray_get_data(infos%dat, OQP_SM,    smat)
    call tagarray_get_data(infos%dat, OQP_TM,    tmat)
    call tagarray_get_data(infos%dat, OQP_Hcore, hcore)

    ! --- allocate main working arrays ---
    allocate( pvp(nbf2),       &
              XU(nbf,nbf),     &
              SXU(nbf,nbf),    &
              psq(nbf),        &
              Ep(nbf),         &
              Akin(nbf),       &
              Rkin(nbf),       &
              hdk(nbf2),       &
              stat=ok)
    if (ok /= 0) call show_message('dk_scalar: cannot allocate', WITH_ABORT)

    ! --- Step 1: compute <mu|pVp|nu> integrals ---
    call compute_and_check_pvp(basis, infos, pvp)

    if (dk_debug) then
      write(iw, '(/,a)') '  Hcore diagonal d-block (14-19):'
      do i = 14, 19
        write(iw, '(2x,i5,es16.6)') i, hcore(i*(i-1)/2 + i)
      end do
      write(iw, '(a)') '  Hcore(17,17), (18,18), (19,19) vs (14,14):'
      write(iw, '(3es16.6)') hcore(17*16/2+17), hcore(14*13/2+14)
    end if

    ! --- Step 2: build p-space basis ---
    call build_p_space(smat, tmat, nbf, XU, SXU, psq, qrnk)

    if (dk_debug) then
      write(iw, '(a,2i5)') '  nbf, qrnk = ', nbf, qrnk
    end if

    call check_p_space(smat, tmat, nbf, qrnk, XU, SXU, psq)

    allocate( Vp(qrnk*(qrnk+1)/2), PVPp(qrnk*(qrnk+1)/2), hdkp(qrnk*(qrnk+1)/2), stat=ok )
    if (ok /= 0) call show_message('dk_scalar: cannot allocate p-space arrays', WITH_ABORT)

    ! --- Step 3: kinematic factors E_p, A, R ---
    call compute_kinematic_factors(psq, qrnk, Ep, Akin, Rkin)

    if (dk_debug) then
      write(iw, '(a,es12.4)') '  PVP min diagonal: ', minval([(pvp(i*(i+1)/2), i=1,nbf)])
      write(iw, '(a,es16.6)') '  Trace pvp: ',        sum([(pvp(i*(i+1)/2), i=1,nbf)])
    end if

    ! --- Step 4: transform V and pVp to p-space ---
    call transform_to_p_space(hcore, tmat, pvp, XU, nbf, qrnk, Vp, PVPp)

    if (dk_debug) then
      write(iw, '(/,a)') '  Vp diagonal (first 5):'
      do i = 1, min(5, qrnk)
        write(iw, '(2x,i5,es16.6)') i, Vp(i*(i+1)/2)
      end do
    end if

    ! --- Steps 5-6: build H^DK1, then add H^DK2 correction if requested ---
    call build_hdk_p(Ep, Akin, Rkin, Vp, PVPp, qrnk, hdkp)
    if (infos%control%scal_rel >= 2) &
      call build_hdk2_p(Ep, Akin, Rkin, psq, Vp, PVPp, qrnk, hdkp)

    ! --- Step 7: back-transform to AO basis ---
    call back_transform_hdk(hdkp, SXU, nbf, qrnk, hdk)

    if (dk_debug) then
      write(iw, '(/,a)') '  === NR limit check: hdk vs hcore ==='
      write(iw, '(a,es12.4)') '  Max |hdk - hcore|: ', maxval(abs(hdk(1:nbf2) - hcore(1:nbf2)))

      write(iw, '(/,a)') '  hcore diagonal (first 5):'
      do i = 1, min(5, nbf)
        write(iw, '(2x,i5,es16.6)') i, hcore(i*(i+1)/2)
      end do

      write(iw, '(/,a)') '  H_DK matrix:'
      do i = 1, nbf
        do j = 1, i
          idx = i*(i-1)/2 + j
          write(iw, '(2i5, f20.10)') i, j, hdk(idx)
        end do
      end do
    end if

    ! --- overwrite OQP::Hcore with H^DK ---
    hcore(:) = hdk(:)

    deallocate(pvp, XU, SXU, psq, Ep, Akin, Rkin, hdk, Vp, PVPp, hdkp)

    write(iw,'(/1X,"...... End Of DK Scalar Correction ......"/)')
    close(iw)

  end subroutine dk_scalar

  !> @brief Compute the <mu|pVp|nu> integrals and apply AO normalisation
  !>
  !> @details Loops over shell pairs and accumulates the momentum-weighted
  !>          nuclear attraction integrals:
  !>
  !>            (pVp)_{mu nu} = <mu| p * (-sum_A Z_A^eff/r_A) * p |nu>
  !>
  !>          The result is stored in packed triangular form and normalised
  !>          with bas_norm_matrix (accounts for the sqrt(3) factor for
  !>          Cartesian d-functions).  When dk_debug is enabled, symmetry
  !>          and positivity of the diagonal are verified.
  !>
  !> @param[in]    basis   Basis set descriptor
  !> @param[inout] infos   OQP information struct (atoms, ECP charges, log)
  !> @param[out]   pvp     pVp matrix, packed lower-triangular (nbf*(nbf+1)/2)
  subroutine compute_and_check_pvp(basis, infos, pvp)

    use types,             only: information
    use basis_tools,       only: basis_set, bas_norm_matrix
    use mod_1e_primitives, only: comp_pvp_int1_prim, update_triang_matrix
    use mod_shell_tools,   only: shell_t, shpair_t
    use constants,         only: tol_int
    use messages,          only: show_message, with_abort
    use precision,         only: dp
    use io_constants,      only: iw
    use printing,          only: print_sym_labeled

    implicit none

    type(basis_set),   intent(in)    :: basis
    type(information), intent(inout) :: infos
    real(dp),          intent(out)   :: pvp(:)  ! packed triangular, size nbf2

    type(shell_t)  :: shi, shj
    type(shpair_t) :: cntp

    integer  :: ii, jj, iat, nat, nbf, nbf2, ig, ok
    real(dp) :: tol
    integer, parameter :: blocksize = 28*28  ! max Cartesian functions per shell pair

    real(dp), allocatable :: pvpmat(:)    ! accumulated pVp, packed triangular
    real(dp), allocatable :: pvpfull(:,:) ! unpacked pVp for symmetry check (debug)

    real(dp) :: pvpblk(blocksize)   ! primitive-level buffer for one shell pair

    real(dp) :: sym_err, diag_min
    integer  :: mu, nu, idx_mu_nu

    nbf  = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nat  = ubound(infos%atoms%zn, 1)
    tol  = log(10.0_dp) * tol_int

    allocate(pvpmat(nbf2), stat=ok)
    if (ok /= 0) call show_message('compute_and_check_pvp: cannot allocate', with_abort)
    pvpmat = 0.0_dp

    ! --- loop over shell pairs, accumulate pVp contributions from all nuclei ---
    call cntp%alloc(basis)

    do ii = basis%nshell, 1, -1
        call shi%fetch_by_id(basis, ii)
        do jj = 1, ii
            call shj%fetch_by_id(basis, jj)

            call cntp%shell_pair(basis, shi, shj, tol)
            if (cntp%numpairs == 0) cycle

            pvpblk = 0.0_dp

            do iat = 1, nat
                do ig = 1, cntp%numpairs
                    ! charge weight: -(Z - Z_ecp) so that V = -sum_A Z_eff/r_A
                    call comp_pvp_int1_prim( &
                        cntp, ig,            &
                        infos%atoms%xyz(:, iat), &
                        -(infos%atoms%zn(iat) - infos%basis%ecp_zn_num(iat)), &
                        pvpblk)
                end do
            end do

            call update_triang_matrix(shi, shj, pvpblk, pvpmat)
        end do
    end do

    if (dk_debug) then
      ! --- check 1: pVp must be symmetric ---
      allocate(pvpfull(nbf, nbf), stat=ok)
      if (ok /= 0) call show_message('compute_and_check_pvp: cannot allocate pvpfull', with_abort)
      pvpfull = 0.0_dp
      do mu = 1, nbf
        do nu = 1, mu
          idx_mu_nu = mu*(mu-1)/2 + nu
          pvpfull(mu, nu) = pvpmat(idx_mu_nu)
          pvpfull(nu, mu) = pvpmat(idx_mu_nu)
        end do
      end do

      sym_err = 0.0_dp
      do mu = 1, nbf
        do nu = 1, nbf
          sym_err = max(sym_err, abs(pvpfull(mu,nu) - pvpfull(nu,mu)))
        end do
      end do

      write(iw, '(/,a)') '  === PVP integral checks ==='
      write(iw, '(a,es12.4)') '  Max symmetry error (should be 0): ', sym_err

      ! --- check 2: diagonal elements must be non-negative ---
      diag_min = huge(1.0_dp)
      do mu = 1, nbf
        idx_mu_nu = mu*(mu-1)/2 + mu
        diag_min = min(diag_min, pvpmat(idx_mu_nu))
      end do
      write(iw, '(a,es12.4)') '  Min diagonal element (should be >= 0): ', diag_min

      deallocate(pvpfull)
    end if

    ! --- apply AO normalisation (accounts for sqrt(3) on d-shell cross terms) ---
    call bas_norm_matrix(pvpmat, basis%bfnrm, nbf)

    if (dk_debug) then
      write(iw, '(/,a)') '  PVP matrix:'
      call print_sym_labeled(pvpmat, nbf, basis)
    end if

    pvp(:) = pvpmat(:)
    deallocate(pvpmat)

  end subroutine compute_and_check_pvp

  !> @brief Build the p-space transformation matrices from S and T
  !>
  !> @details Constructs the momentum-space basis by canonical orthogonalisation
  !>          of the overlap followed by diagonalisation of 2T.  The result is
  !>          a set of vectors X~ such that:
  !>
  !>            X~^T * S * X~ = I_qrnk
  !>            X~^T * 2T * X~ = diag(p_i^2),  i = 1..qrnk
  !>
  !>          Procedure:
  !>            1. X  = S^{-1/2}               (canonical orthogonalisation)
  !>            2. T~ = X^T * 2T * X            (2T in orthonormal basis)
  !>            3. T~ * U = U * diag(p^2)       (diagonalise T~)
  !>            4. XU  = X * U                  (p-space basis in AO rep.)
  !>            5. SXU = S * XU                 (needed for back-transform)
  !>
  !> @param[in]  smat   Overlap matrix, packed triangular (nbf*(nbf+1)/2)
  !> @param[in]  tmat   Kinetic energy matrix, packed triangular
  !> @param[in]  nbf    Number of AO basis functions
  !> @param[out] xu     XU matrix (nbf x nbf);  columns 1:qrnk are valid
  !> @param[out] sxu    SXU = S * XU (nbf x nbf);  columns 1:qrnk are valid
  !> @param[out] psq    p_i^2 eigenvalues (nbf); entries 1:qrnk are valid
  !> @param[out] qrnk   Effective rank (nbf minus linear dependencies)
  subroutine build_p_space(smat, tmat, nbf, xu, sxu, psq, qrnk)

    use mathlib,  only: matrix_invsqrt, orthogonal_transform_sym, unpack_f90
    use eigen,    only: diag_symm_packed
    use messages, only: show_message, with_abort
    use precision, only: dp

    implicit none

    real(dp), intent(in)  :: smat(*), tmat(*)  ! packed triangular, nbf*(nbf+1)/2
    integer,  intent(in)  :: nbf
    real(dp), intent(out) :: xu(nbf, nbf)
    real(dp), intent(out) :: sxu(nbf, nbf)
    real(dp), intent(out) :: psq(nbf)
    integer,  intent(out) :: qrnk

    real(dp), allocatable :: x(:,:)      ! S^{-1/2}  (nbf x nbf)
    real(dp), allocatable :: ttilde(:)   ! X^T * 2T * X, packed triangular
    real(dp), allocatable :: u(:,:)      ! eigenvectors of ttilde
    real(dp), allocatable :: t2(:)       ! 2*T, packed triangular
    real(dp), allocatable :: sfull(:,:)  ! S in full storage (for dsymm)

    integer :: nbf2, ok, ierr

    nbf2 = nbf*(nbf+1)/2

    allocate(x(nbf, nbf),     &
             ttilde(nbf2),    &
             u(nbf, nbf),     &
             t2(nbf2),        &
             sfull(nbf, nbf), &
             stat=ok)
    if (ok /= 0) call show_message('build_p_space: cannot allocate', with_abort)

    ! --- step 1: X = S^{-1/2} ---
    call matrix_invsqrt(smat, x, nbf, qrnk)

    ! --- step 2: T~ = X^T * 2T * X ---
    t2(1:nbf2) = 2.0_dp * tmat(1:nbf2)
    call orthogonal_transform_sym(nbf, qrnk, t2, x, nbf, ttilde)

    ! --- step 3: diagonalise T~ -> p^2 eigenvalues and eigenvectors ---
    call diag_symm_packed(1, qrnk, qrnk, qrnk, ttilde, psq, u, ierr)
    if (ierr /= 0) call show_message('build_p_space: diag_symm_packed failed', with_abort)

    ! --- step 4: XU = X * U ---
    call dgemm('n', 'n', nbf, qrnk, qrnk, &
               1.0_dp, x, nbf, u, qrnk,   &
               0.0_dp, xu, nbf)

    ! --- step 5: SXU = S * XU ---
    call unpack_f90(smat, sfull, 'u')
    call dsymm('l', 'u', nbf, qrnk,          &
                1.0_dp, sfull, nbf, xu, nbf,  &
                0.0_dp, sxu, nbf)

    deallocate(x, ttilde, u, t2, sfull)

  end subroutine build_p_space

  !> @brief Verify the p-space orthonormality conditions (debug only)
  !>
  !> @details When dk_debug is .true., checks:
  !>            - XU^T * S * XU  = I_qrnk   (orthonormality)
  !>            - XU^T * 2T * XU = diag(psq) (eigenvalue condition)
  !>          and prints the first five p^2 values.  Returns immediately
  !>          if dk_debug is .false. (zero cost in production runs).
  !>
  !> @param[in] smat   Overlap matrix, packed triangular
  !> @param[in] tmat   Kinetic energy matrix, packed triangular
  !> @param[in] nbf    Number of AO basis functions
  !> @param[in] qrnk   Effective rank
  !> @param[in] xu     p-space basis vectors in AO rep. (nbf x nbf)
  !> @param[in] sxu    S * XU  (nbf x nbf)
  !> @param[in] psq    p^2 eigenvalues (qrnk)
  subroutine check_p_space(smat, tmat, nbf, qrnk, xu, sxu, psq)

    use precision,    only: dp
    use io_constants, only: iw
    use messages,     only: show_message, with_abort
    use mathlib,      only: unpack_f90

    implicit none

    real(dp), intent(in) :: smat(*), tmat(*)
    integer,  intent(in) :: nbf, qrnk
    real(dp), intent(in) :: xu(nbf, nbf), sxu(nbf, nbf), psq(nbf)

    real(dp), allocatable :: check(:,:), t2full(:,:), tmp(:,:)
    real(dp) :: err1, err2
    integer  :: i, j, ok

    if (.not. dk_debug) return

    allocate(check(nbf,nbf), t2full(nbf,nbf), tmp(nbf,nbf), stat=ok)
    if (ok /= 0) call show_message('check_p_space: cannot allocate', with_abort)

    write(iw, '(/,a)') '  === build_p_space checks ==='

    ! --- test 1: XU^T * S * XU = I ---
    call dgemm('t', 'n', qrnk, qrnk, nbf, &
               1.0_dp, sxu, nbf, xu, nbf,  &
               0.0_dp, check, nbf)

    err1 = 0.0_dp
    do i = 1, qrnk
        do j = 1, qrnk
            if (i == j) then
                err1 = max(err1, abs(check(i,j) - 1.0_dp))
            else
                err1 = max(err1, abs(check(i,j)))
            end if
        end do
    end do
    write(iw, '(a,es12.4)') '  XU^T*S*XU = I,            max error: ', err1

    ! --- test 2: XU^T * 2T * XU = diag(psq) ---
    t2full = 0.0_dp
    call unpack_f90(tmat, t2full, 'u')
    t2full = 2.0_dp * t2full

    call dgemm('n', 'n', nbf, qrnk, nbf, &
               1.0_dp, t2full, nbf, xu, nbf, &
               0.0_dp, tmp, nbf)

    call dgemm('t', 'n', qrnk, qrnk, nbf, &
               1.0_dp, xu, nbf, tmp, nbf, &
               0.0_dp, check, nbf)

    err2 = 0.0_dp
    do i = 1, qrnk
        do j = 1, qrnk
            if (i == j) then
                err2 = max(err2, abs(check(i,j) - psq(i)))
            else
                err2 = max(err2, abs(check(i,j)))
            end if
        end do
    end do
    write(iw, '(a,es12.4)') '  XU^T*2T*XU = diag(psq),  max error: ', err2

    write(iw, '(a)') '  First 5 p^2 values:'
    do i = 1, min(5, qrnk)
        write(iw, '(2x,i5,es16.6)') i, psq(i)
    end do

    deallocate(check, t2full, tmp)

  end subroutine check_p_space

  !> @brief Compute the DK kinematic factors for each p-space eigenvalue
  !>
  !> @details For each momentum-space eigenvalue p_i^2 computes:
  !>
  !>            E_p(i) = c * sqrt(p_i^2 + c^2)   (relativistic energy)
  !>            A(i)   = sqrt((E_p + c^2) / (2*E_p))
  !>            R(i)   = c / (E_p + c^2)
  !>
  !>          where c = 137.0359895 a.u. (speed of light).
  !>
  !> @param[in]  psq    p_i^2 eigenvalues (qrnk)
  !> @param[in]  qrnk   Number of p-space basis vectors
  !> @param[out] ep     Relativistic kinetic energy E_p (qrnk)
  !> @param[out] akin   Kinematic factor A (qrnk)
  !> @param[out] rkin   Kinematic factor R (qrnk)
  subroutine compute_kinematic_factors(psq, qrnk, ep, akin, rkin)

    use precision,    only: dp
    use io_constants, only: iw

    implicit none

    real(dp), intent(in)  :: psq(*)
    integer,  intent(in)  :: qrnk
    real(dp), intent(out) :: ep(*), akin(*), rkin(*)

    real(dp), parameter :: clight  = 137.0359895_dp
    real(dp), parameter :: clight2 = clight * clight

    integer :: i

    do i = 1, qrnk
        ep(i)   = clight * sqrt(psq(i) + clight2)
        akin(i) = sqrt((ep(i) + clight2) / (2.0_dp * ep(i)))
        rkin(i) = clight / (ep(i) + clight2)
    end do

    if (dk_debug) then
      write(iw, '(/,a)') '  === Kinematic factors (first 5) ==='
      write(iw, '(2x,a5,3a16)') 'i', 'p^2', 'E_p - c^2', 'A_i'
      do i = 1, min(5, qrnk)
        write(iw, '(2x,i5,3es16.6)') i, psq(i), ep(i)-clight2, akin(i)
      end do
      write(iw, '(a)') '  Last 5:'
      do i = max(1, qrnk-4), qrnk
        write(iw, '(2x,i5,3es16.6)') i, psq(i), ep(i)-clight2, akin(i)
      end do
    end if

  end subroutine compute_kinematic_factors

  !> @brief Transform V and pVp from AO basis to p-space
  !>
  !> @details Applies the congruence transformation X~^T * M * X~
  !>          to both the potential V = Hcore - T and the pVp matrix:
  !>
  !>            V^p   = XU^T * V   * XU
  !>            pVp^p = XU^T * pVp * XU
  !>
  !>          Both results are stored as packed lower-triangular arrays
  !>          of size qrnk*(qrnk+1)/2.
  !>
  !> @param[in]  hcore  H_core = T+V, packed triangular AO (nbf2)
  !> @param[in]  tmat   Kinetic energy T, packed triangular AO (nbf2)
  !> @param[in]  pvp    pVp integrals, packed triangular AO (nbf2)
  !> @param[in]  xu     p-space basis XU (nbf x nbf)
  !> @param[in]  nbf    Number of AO basis functions
  !> @param[in]  qrnk   Effective rank
  !> @param[out] vp     V in p-space, packed triangular (qrnk*(qrnk+1)/2)
  !> @param[out] pvpp   pVp in p-space, packed triangular (qrnk*(qrnk+1)/2)
  subroutine transform_to_p_space(hcore, tmat, pvp, xu, nbf, qrnk, vp, pvpp)

    use mathlib,      only: orthogonal_transform_sym
    use messages,     only: show_message, with_abort
    use precision,    only: dp
    use io_constants, only: iw

    implicit none

    real(dp), intent(in)  :: hcore(*), tmat(*), pvp(*)  ! packed triangular AO
    real(dp), intent(in)  :: xu(nbf, nbf)
    integer,  intent(in)  :: nbf, qrnk
    real(dp), intent(out) :: vp(*)    ! V in p-space, packed triangular
    real(dp), intent(out) :: pvpp(*)  ! pVp in p-space, packed triangular

    integer :: nbf2, ok
    real(dp), allocatable :: vao(:)  ! V = Hcore - T in AO basis

    nbf2 = nbf*(nbf+1)/2

    allocate(vao(nbf2), stat=ok)
    if (ok /= 0) call show_message('transform_to_p_space: cannot allocate', with_abort)

    ! --- V = Hcore - T ---
    vao(1:nbf2) = hcore(1:nbf2) - tmat(1:nbf2)

    ! --- V^p = XU^T * V * XU ---
    call orthogonal_transform_sym(nbf, qrnk, vao, xu, nbf, vp)

    ! --- (pVp)^p = XU^T * pVp * XU ---
    call orthogonal_transform_sym(nbf, qrnk, pvp, xu, nbf, pvpp)

    if (dk_debug) then
      write(iw, '(/,a)') '  === transform_to_p_space done ==='
    end if

    deallocate(vao)

  end subroutine transform_to_p_space

  !> @brief Build the DK1 Hamiltonian in p-space
  !>
  !> @details Constructs the first-order Douglas-Kroll Hamiltonian matrix
  !>          in the momentum representation (packed lower-triangular):
  !>
  !>            H^DK1_{ij} = (E_p(i) - c^2) * delta_{ij}
  !>                       + A(i) * V^p_{ij} * A(j)
  !>                       + A(i)*R(i) * (pVp)^p_{ij} * R(j)*A(j)
  !>
  !> @param[in]  ep     Relativistic energy E_p (qrnk)
  !> @param[in]  akin   Kinematic factor A (qrnk)
  !> @param[in]  rkin   Kinematic factor R (qrnk)
  !> @param[in]  vp     V in p-space, packed triangular (qrnk*(qrnk+1)/2)
  !> @param[in]  pvpp   pVp in p-space, packed triangular
  !> @param[in]  qrnk   Effective rank
  !> @param[out] hdkp   H^DK1 in p-space, packed triangular
  subroutine build_hdk_p(ep, akin, rkin, vp, pvpp, qrnk, hdkp)

    use precision,    only: dp
    use io_constants, only: iw

    implicit none

    real(dp), intent(in)  :: ep(*), akin(*), rkin(*)
    real(dp), intent(in)  :: vp(*), pvpp(*)
    integer,  intent(in)  :: qrnk
    real(dp), intent(out) :: hdkp(*)

    real(dp), parameter :: clight  = 137.0359895_dp
    real(dp), parameter :: clight2 = clight * clight

    integer  :: i, j, ij
    real(dp) :: ar_i, ar_j

    ij = 0
    do i = 1, qrnk
        ar_i = akin(i) * rkin(i)
        do j = 1, i
            ij = ij + 1
            ar_j = akin(j) * rkin(j)

            hdkp(ij) = akin(i) * vp(ij) * akin(j)   &   ! A * V^p * A
                     + ar_i * pvpp(ij) * ar_j             ! A*R * (pVp)^p * R*A

            ! add kinetic energy contribution on diagonal
            if (i == j) hdkp(ij) = hdkp(ij) + ep(i) - clight2

        end do
    end do

    if (dk_debug) then
      write(iw, '(/,a)') '  === H^DK in p-space, diagonal (first 5) ==='
      do i = 1, min(5, qrnk)
        write(iw, '(2x,i5,3es16.6)') i, ep(i)-clight2, hdkp(i*(i+1)/2)
      end do
    end if

  end subroutine build_hdk_p

  !> @brief Back-transform H^DK from p-space to AO basis
  !>
  !> @details Applies the two-sided transformation:
  !>
  !>            H^DK_AO = SXU * H^DK_p * SXU^T
  !>
  !>          where SXU = S * XU.  The result is symmetrised and packed
  !>          into a lower-triangular array.
  !>
  !> @param[in]  hdkp   H^DK in p-space, packed triangular (qrnk*(qrnk+1)/2)
  !> @param[in]  sxu    S * XU matrix (nbf x nbf);  columns 1:qrnk are used
  !> @param[in]  nbf    Number of AO basis functions
  !> @param[in]  qrnk   Effective rank
  !> @param[out] hdk    H^DK in AO basis, packed triangular (nbf*(nbf+1)/2)
  subroutine back_transform_hdk(hdkp, sxu, nbf, qrnk, hdk)

    use mathlib,      only: unpack_f90, pack_f90
    use messages,     only: show_message, with_abort
    use precision,    only: dp
    use io_constants, only: iw

    implicit none

    real(dp), intent(in)  :: hdkp(*)
    real(dp), intent(in)  :: sxu(nbf, nbf)
    integer,  intent(in)  :: nbf, qrnk
    real(dp), intent(out) :: hdk(*)

    real(dp), allocatable :: hdkp_full(:,:)  ! H^DK_p in full storage
    real(dp), allocatable :: tmp(:,:)        ! intermediate  SXU * H^DK_p
    real(dp), allocatable :: hdk_full(:,:)   ! H^DK_AO in full storage

    integer :: ok, i, j

    allocate(hdkp_full(qrnk, qrnk), &
             tmp(nbf, qrnk),        &
             hdk_full(nbf, nbf),    &
             stat=ok)
    if (ok /= 0) call show_message('back_transform_hdk: cannot allocate', with_abort)

    ! --- unpack H^DK_p and symmetrise ---
    call unpack_f90(hdkp, hdkp_full, 'u')
    do i = 1, qrnk
        do j = i+1, qrnk
            hdkp_full(j,i) = hdkp_full(i,j)
        end do
    end do

    ! --- tmp = SXU * H^DK_p ---
    call dgemm('n', 'n', nbf, qrnk, qrnk, &
               1.0_dp, sxu, nbf, hdkp_full, qrnk, &
               0.0_dp, tmp, nbf)

    ! --- H^DK_AO = tmp * SXU^T ---
    call dgemm('n', 't', nbf, nbf, qrnk, &
               1.0_dp, tmp, nbf, sxu, nbf, &
               0.0_dp, hdk_full, nbf)

    if (dk_debug) then
      write(iw, '(/,a)') '  hdk_full diagonal (first 5):'
      do i = 1, min(5, nbf)
        write(iw, '(2x,i5,3es16.6)') i, hdk_full(i,i)
      end do
    end if

    ! --- pack to lower-triangular ---
    call pack_f90(hdk_full, hdk(1:nbf*(nbf+1)/2), 'u')

    if (dk_debug) then
      write(iw, '(/,a)') '  === back_transform_hdk done ==='
      write(iw, '(a)') '  H^DK diagonal (first 5):'
      do i = 1, min(5, nbf)
        write(iw, '(2x,i5,es16.6)') i, hdk(i*(i+1)/2)
      end do
    end if

    deallocate(hdkp_full, tmp, hdk_full)

  end subroutine back_transform_hdk

  !> @brief Add the DK2 second-order correction to H^DK in p-space
  !>
  !> @details Computes and accumulates the DK2 correction following the
  !>          direct W_1^2 approach (equivalent to GAMESS DK2X):
  !>
  !>            H^DK2 += -1/2 (E * W_1^2 + W_1^2 * E)  -  W_1 * E * W_1
  !>
  !>          where W_1 is the first-order unitary generator.  Rather than
  !>          forming W_1 explicitly, W_1^2 and W_1*E*W_1 are assembled
  !>          from four matrix-matrix products each, using scaled V^p and
  !>          pVp^p matrices:
  !>
  !>            V~_{ij}   = V^p_{ij}   / (E_p(i) + E_p(j))
  !>            pVp~_{ij} = pVp^p_{ij} / (E_p(i) + E_p(j))
  !>
  !>          Block 1 - W_1^2 (four terms):
  !>            +  (AR*pVp~*RA) * (A*V~*A)
  !>            +  (A*V~*A)     * (AR*pVp~*RA)
  !>            -  (A*V~*p^2*R^2*A) * (A*V~*A)
  !>            -  (AR*pVp~*A/p^2)  * (A*pVp~*RA)
  !>
  !>          Block 2 - W_1*E*W_1 (four analogous terms with extra E factors).
  !>
  !>          The result is symmetrised and added to hdkp in-place.
  !>
  !> @param[in]    ep     Relativistic energy E_p (qrnk)
  !> @param[in]    akin   Kinematic factor A (qrnk)
  !> @param[in]    rkin   Kinematic factor R (qrnk)
  !> @param[in]    psq    p^2 eigenvalues (qrnk)
  !> @param[in]    vp     V in p-space, packed triangular
  !> @param[in]    pvpp   pVp in p-space, packed triangular
  !> @param[in]    qrnk   Effective rank
  !> @param[inout] hdkp   H^DK in p-space (DK2 correction accumulated in-place)
  subroutine build_hdk2_p(ep, akin, rkin, psq, vp, pvpp, qrnk, hdkp)

    use precision, only: dp
    use messages,  only: show_message, with_abort
    implicit none

    real(dp), intent(in)    :: ep(*), akin(*), rkin(*), psq(*)
    real(dp), intent(in)    :: vp(*), pvpp(*)
    integer,  intent(in)    :: qrnk
    real(dp), intent(inout) :: hdkp(*)

    real(dp), allocatable :: vps(:,:)    ! V~   = V^p  / (Ei+Ej),  full
    real(dp), allocatable :: pvpps(:,:)  ! pVp~ = pVp^p/ (Ei+Ej),  full
    real(dp), allocatable :: ma(:,:)     ! scratch matrix A
    real(dp), allocatable :: mb(:,:)     ! scratch matrix B
    real(dp), allocatable :: w1sq(:,:)   ! W_1^2  accumulator
    real(dp), allocatable :: hdk2(:,:)   ! DK2 correction (full, before pack)

    integer  :: i, j, ij, ok

    allocate(vps(qrnk,qrnk), pvpps(qrnk,qrnk), &
             ma(qrnk,qrnk),  mb(qrnk,qrnk),    &
             w1sq(qrnk,qrnk), hdk2(qrnk,qrnk), &
             stat=ok)
    if (ok /= 0) call show_message('build_hdk2_p: cannot allocate', with_abort)

    ! --- scale V^p and pVp^p by 1/(Ei+Ej) ---
    ij = 0
    do i = 1, qrnk
        do j = 1, i
            ij = ij + 1
            vps(i,j)   = vp(ij)   / (ep(i) + ep(j))
            vps(j,i)   = vps(i,j)
            pvpps(i,j) = pvpp(ij) / (ep(i) + ep(j))
            pvpps(j,i) = pvpps(i,j)
        end do
    end do

    ! =========================================================
    ! Block 1: build W_1^2 from four terms
    ! =========================================================
    w1sq = 0.0_dp

    ! --- term 1: (AR*pVp~*RA) * (A*V~*A) ---
    do i = 1, qrnk
        do j = 1, i
            ma(i,j) = akin(i)*rkin(i) * pvpps(i,j) * rkin(j)*akin(j)
            ma(j,i) = ma(i,j)
            mb(i,j) = akin(i) * vps(i,j) * akin(j)
            mb(j,i) = mb(i,j)
        end do
    end do
    call dgemm('n','n', qrnk,qrnk,qrnk, 1.0_dp, ma,qrnk, mb,qrnk, 0.0_dp, w1sq,qrnk)

    ! --- term 2: (A*V~*A) * (AR*pVp~*RA) ---
    call dgemm('n','n', qrnk,qrnk,qrnk, 1.0_dp, mb,qrnk, ma,qrnk, 1.0_dp, w1sq,qrnk)

    ! --- term 3: -(A*V~*p^2*R^2*A) * (A*V~*A) ---
    do i = 1, qrnk
        do j = 1, i
            ma(i,j) = akin(i) * vps(i,j) * akin(j) * psq(j)*rkin(j)**2
            ma(j,i) = akin(j) * vps(i,j) * akin(i) * psq(i)*rkin(i)**2
        end do
    end do
    call dgemm('n','n', qrnk,qrnk,qrnk, -1.0_dp, ma,qrnk, mb,qrnk, 1.0_dp, w1sq,qrnk)

    ! --- term 4: -(AR*pVp~*A/p^2) * (A*pVp~*RA) ---
    do i = 1, qrnk
        do j = 1, i
            ma(i,j) = akin(i)*rkin(i) * pvpps(i,j) * akin(j) / psq(j)
            ma(j,i) = akin(j)*rkin(j) * pvpps(i,j) * akin(i) / psq(i)
            mb(i,j) = akin(i) * pvpps(i,j) * rkin(j)*akin(j)
            mb(j,i) = akin(j) * pvpps(i,j) * rkin(i)*akin(i)
        end do
    end do
    call dgemm('n','n', qrnk,qrnk,qrnk, -1.0_dp, ma,qrnk, mb,qrnk, 1.0_dp, w1sq,qrnk)

    ! --- accumulate -1/2*(E*W1sq + W1sq*E) into hdk2 ---
    ! E is diagonal: (E*W1sq)_{ij} = ep(i)*w1sq(i,j)
    do i = 1, qrnk
        do j = 1, qrnk
            hdk2(i,j) = -0.5_dp * (ep(i)*w1sq(i,j) + w1sq(i,j)*ep(j))
        end do
    end do

    ! =========================================================
    ! Block 2: -W_1*E*W_1 from four terms (same structure, extra ep factor)
    ! =========================================================

    ! --- term 1: -(AR*pVp~*RA*E) * (A*V~*A) ---
    do i = 1, qrnk
        do j = 1, i
            ma(i,j) = akin(i)*rkin(i) * pvpps(i,j) * rkin(j)*akin(j) * ep(j)
            ma(j,i) = akin(j)*rkin(j) * pvpps(i,j) * rkin(i)*akin(i) * ep(i)
            mb(i,j) = akin(i) * vps(i,j) * akin(j)
            mb(j,i) = mb(i,j)
        end do
    end do
    call dgemm('n','n', qrnk,qrnk,qrnk, -1.0_dp, ma,qrnk, mb,qrnk, 1.0_dp, hdk2,qrnk)

    ! --- term 2: -(A*V~*A) * (AR*pVp~*RA*E)^T ---
    call dgemm('n','t', qrnk,qrnk,qrnk, -1.0_dp, mb,qrnk, ma,qrnk, 1.0_dp, hdk2,qrnk)

    ! --- term 3: +(A*V~*p^2*R^2*E*A) * (A*V~*A) ---
    do i = 1, qrnk
        do j = 1, i
            ma(i,j) = akin(i) * vps(i,j) * akin(j) * psq(j)*rkin(j)**2 * ep(j)
            ma(j,i) = akin(j) * vps(i,j) * akin(i) * psq(i)*rkin(i)**2 * ep(i)
        end do
    end do
    call dgemm('n','n', qrnk,qrnk,qrnk, 1.0_dp, ma,qrnk, mb,qrnk, 1.0_dp, hdk2,qrnk)

    ! --- term 4: +(AR*pVp~*A*E/p^2) * (A*pVp~*RA) ---
    do i = 1, qrnk
        do j = 1, i
            ma(i,j) = akin(i)*rkin(i) * pvpps(i,j) * akin(j) * ep(j) / psq(j)
            ma(j,i) = akin(j)*rkin(j) * pvpps(i,j) * akin(i) * ep(i) / psq(i)
            mb(i,j) = akin(i) * pvpps(i,j) * rkin(j)*akin(j)
            mb(j,i) = akin(j) * pvpps(i,j) * rkin(i)*akin(i)
        end do
    end do
    call dgemm('n','n', qrnk,qrnk,qrnk, 1.0_dp, ma,qrnk, mb,qrnk, 1.0_dp, hdk2,qrnk)

    ! --- symmetrise and accumulate into hdkp (packed triangular) ---
    ij = 0
    do i = 1, qrnk
        do j = 1, i
            ij = ij + 1
            hdkp(ij) = hdkp(ij) + 0.5_dp * (hdk2(i,j) + hdk2(j,i))
        end do
    end do

    deallocate(vps, pvpps, ma, mb, w1sq, hdk2)

  end subroutine build_hdk2_p


end module dk_scalar_mod