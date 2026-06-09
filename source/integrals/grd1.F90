module grd1

   use io_constants, only: iw
   use precision, only: dp
   use types, only: information
   use atomic_structure_m, only: atomic_structure

   use basis_tools, only: basis_set, &
       bas_norm_matrix, bas_denorm_matrix, build_cart_density
   use constants, only: HARMONIC_ACTIVE

   use mod_1e_primitives, only: &
       comp_coulomb_der1, comp_coulomb_helfeyder1, comp_kinetic_der1, &
       comp_overlap_der1, &
       comp_overlap_der2, comp_kinetic_der2, comp_coulomb_der2_braC, &
       comp_overlap_der1_block, comp_kinetic_der1_block, &
       comp_coulomb_der1_block, comp_coulomb_helfeyder1_block, &
       comp_ewaldlr_der1, comp_ewaldlr_helfeyder1, &
       density_ordered

   use mod_shell_tools, only: shell_t, shpair_t
   use mathlib, only: unpack_matrix
   use ecp_tool, only: add_ecpder
   implicit none

   character(len=*), parameter :: module_name = "grd1"

   real(kind=dp), parameter :: tol_default = log(10.0d0)*20

   private
   public print_gradient
   public eijden
   public grad_nn
   public hess_nn
   public grad_ee_overlap
   public grad_ee_kinetic
   public hess_ee_overlap
   public hess_ee_kinetic
   public hess_en
   public der_overlap_matrix
   public der_kinetic_matrix
   public der_nucattr_matrix
   public grad_en_hellman_feynman
   public grad_en_pulay
   public grad_1e_ecp

contains

!-------------------------------------------------------------------------------

!> @brief Compute "energy weighted density matrix"
!> @note This quantity is actually the Lagrangian matrix,
!   backtransformed into the AO basis.
  subroutine eijden(eps, nbf, infos)
    use oqp_tagarray_driver
    use mathlib, only: orthogonal_transform_sym
    use messages, only: show_message, with_abort

    implicit none

    character(len=*), parameter :: subroutine_name = "eijden"

    type(information), intent(inout) :: infos
    integer :: nbf
    integer :: i, j, ij, ne, ok
    real(kind=dp) :: eps(:)
    real(kind=dp), allocatable :: c(:,:), scr(:),tempd(:)

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      mo_energy_a(:), mo_a(:,:), &
      fock_a(:), fock_b(:), dmat_a(:), dmat_b(:)
    character(len=*), parameter :: tags_alpha(2) = (/ character(len=80) :: &
      OQP_E_MO_A, OQP_VEC_MO_A /)
    character(len=*), parameter :: tags_beta(4) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_DM_A, OQP_FOCK_B, OQP_DM_B /)

    if (infos%control%scftype>1) then
       allocate(c(nbf,nbf), scr(8*nbf), tempd(nbf*(nbf+1)/2), stat=ok)
    end if

    ne = infos%mol_prop%nelec/2

    select case (infos%control%scftype)
!   RHF case
    case (1)

      call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
      call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

       ij = 0
       do i = 1, nbf
          do j = 1, i
             ij = ij+1
             eps(ij) = -2*sum(mo_energy_a(1:ne)&
                             *mo_a(i,1:ne)&
                             *mo_a(j,1:ne))
          end do
       end do

!   U/ROHF case
    case (2:)

      call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
      call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
      call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)

!     Alpha part
      call unpack_matrix(dmat_a,c,nbf,'U')
      call orthogonal_transform_sym(nbf, nbf, fock_a, c, nbf, tempd)

!     Beta part
      call unpack_matrix(dmat_b,c,nbf,'U')
      call orthogonal_transform_sym(nbf, nbf, fock_b, c, nbf, eps)
      eps = -eps - tempd

!     Half the diagonal
      ij = 0
      do i = 1, nbf
         ij = ij+i
         eps(ij) = 0.5d0*eps(ij)
      end do

    end select

  end subroutine eijden

!-------------------------------------------------------------------------------

!> @brief Print energy gradient vector
  subroutine grad_max_rms(n,de,gmax,grms)

    implicit none
!    type(information), intent(in) :: infos
    real(kind=dp), intent(out) :: gmax, grms
    real(kind=dp) :: de(3,n)
    integer :: n, i

  ! Calculate maximum value
  gmax = maxval(abs(de))

  ! Calculate root mean square (RMS)
  grms = 0.0
  do i = 1, n
    grms = grms + de(1,i)**2 + de(2,i) ** 2 + de(3,i) ** 2
  end do
  grms = sqrt(grms / real(n * 3))

  end subroutine grad_max_rms

!-------------------------------------------------------------------------------

!> @brief Print energy gradient vector
  subroutine print_gradient(infos)

    implicit none
    type(information), intent(in) :: infos
    real(kind=dp) :: gmax, grms
    integer :: i

    write(iw, fmt="(&
              &/25X,23('=')&
              &/25X,'Gradient (Hartree/Bohr)'&
              &/25X,23('=')&
              &/8X,'ATOM     ZNUC',9X,'dE/dX',10X,'dE/dY',10X,'dE/dZ'&
              &/6X,62('-'))")

    do i = 1, infos%mol_prop%natom
       write(iw,'(7X,I4,5X,F4.1,3X,3F15.9)') &
               i,infos%atoms%zn(i), infos%atoms%grad(:,i)
    end do

!   Compute Maximum and RMS Gradient
    call grad_max_rms(infos%mol_prop%natom,infos%atoms%grad,gmax,grms)
    write(iw,fmt="(/10X,'Maximum Gradient =',F10.7,4X,&
          &'RMS Gradient =',F10.7/)") gmax, grms

  end subroutine print_gradient

!-------------------------------------------------------------------------------

!> @brief Print energy gradient vector, without `atoms`
  subroutine print_grd(zn, de)
    implicit none
    real(kind=dp) :: zn(:), de(:,:)
    real(kind=dp) :: gmax, grms
    integer :: i, j

    write(iw, fmt="(&
              &/25X,23('=')&
              &/25X,'Gradient (Hartree/Bohr)'&
              &/25X,23('=')&
              &/8X,'ATOM     ZNUC',9X,'dE/dX',10X,'dE/dY',10X,'dE/dZ'&
              &/6X,62('-'))")

    do i = 1, ubound(de,2)
       write(iw,'(7x,i4,5x,f4.1,3x,3f15.9)') i, zn(i), (de(j,i),j=1,3)
    end do

!   Compute Maximum and RMS Gradient
    call grad_max_rms(ubound(de,2),de,gmax,grms)
    write(iw,fmt="(/10X,'Maximum Gradient =',F10.7,4X,&
          &'RMS Gradient =',F10.7/)") gmax, grms

  end subroutine


!-------------------------------------------------------------------------------

!> @brief Compute overlap energy derivative contribution to gradient
!
!> @author Vladimir Mironov
!
!   REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @brief Unpack + bfnrm-fold a packed density and, under HARMONIC_ACTIVE,
!>        expand it to the Cartesian-effective density the derivative kernels
!>        contract. Returns the (possibly Cartesian-sized) full density and
!>        the matching per-shell AO offsets. With the gate off this is the
!>        former inline unpack/bas_norm and off = basis%ao_offset.
 SUBROUTINE prepare_grad_density(basis, denab, dens, off)
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: denab(:)
    real(kind=dp), allocatable, intent(out) :: dens(:,:)
    integer, allocatable, intent(out) :: off(:)
    real(kind=dp), allocatable :: dcart(:,:)
    integer, allocatable :: cart_off(:)
    integer :: nbf_cart

    allocate(dens(basis%nbf, basis%nbf), source=0.0d0)
    call unpack_matrix(denab, dens)
    call bas_norm_matrix(dens, basis%bfnrm, basis%nbf)
    if (HARMONIC_ACTIVE) then
      call build_cart_density(basis, dens, dcart, cart_off, nbf_cart)
      call move_alloc(dcart, dens)
      call move_alloc(cart_off, off)
    else
      allocate(off(basis%nshell))
      off = basis%ao_offset(1:basis%nshell)
    end if
 END SUBROUTINE

!> @param[in,out]   denab   density matrix in packed format, remains unchanged on return
 SUBROUTINE grad_ee_overlap(basis, denab, de, logtol)

    implicit none

    type(basis_set), intent(inout) :: basis
    REAL(kind=dp), INTENT(INOUT) :: denab(:)
    real(kind=dp), optional :: logtol

    REAL(kind=dp) :: de(:,:)

    INTEGER :: ii, jj

    REAL(kind=dp) :: tol

    REAL(kind=dp) :: de_atom(3)

    LOGICAL :: norm

    REAL(kind=dp), ALLOCATABLE :: de_priv(:,:), dens(:,:)
    INTEGER, ALLOCATABLE :: off(:)

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    if (present(logtol)) then
        tol = logtol
    else
        tol = tol_default
    end if

    allocate(de_priv, mold=de)
    de_priv = 0.0d0

    call prepare_grad_density(basis, denab, dens, off)

!   Initialize parallel
!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       shi, shj, cntp, &
!$omp       de_atom &
!$omp   ) &
!$omp   reduction(+:de_priv)

    CALL cntp%alloc(basis)

!$omp do schedule(dynamic)
!   I shell
    DO ii = 1, basis%nshell

        de_atom = 0.0

        CALL shi%fetch_by_id(basis, ii)

!       J shell
        DO jj = 1, basis%nshell

            CALL shj%fetch_by_id(basis, jj)

            CALL cntp%shell_pair(basis, shi, shj, tol, dup=.false.)
            IF (cntp%numpairs==0) CYCLE

            CALL comp_overlap_der1(cntp, dens(off(ii):, off(jj):), de_atom)
        END DO

        ! Update gradient
        de_priv(:,shi%atid) = de_priv(:,shi%atid) + 2*de_atom

    END DO
!$omp end do
!   End of shell loops
!$omp end parallel

    de = de + de_priv
    DEALLOCATE(de_priv)

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Basis function derivative contributions to gradient
!> @details Compute derivative integrals of type <ii'|h|jj> = <ii'|t+v|jj>
!> @note No relativistic methods available
!
!> @author Vladimir Mironov
!
!   REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   denab   density matrix in packed format, remains unchanged on return
 SUBROUTINE grad_ee_kinetic(basis, denab, de, logtol)

    REAL(kind=dp), INTENT(INOUT) :: denab(:)
    type(basis_set), intent(inout) :: basis

    REAL(kind=dp) :: de(:,:)

    INTEGER :: ii, jj

    REAL(kind=dp), optional :: logtol

    REAL(kind=dp) :: de_atom(3)
    REAL(kind=dp), ALLOCATABLE :: de_priv(:,:), dens(:,:)
    INTEGER, ALLOCATABLE :: off(:)

    REAL(kind=dp) :: tol
    LOGICAL :: out, dbg, norm

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    dbg = .false.
    out = .false.

    IF (dbg) WRITE(iw,'(/10X,38(1H-)/10X,"GRADIENT INCLUDING AO DERIVATIVE TERMS"/10X,38(1H-))')

    if (present(logtol)) then
        tol = logtol
    else
        tol = tol_default
    end if

    call prepare_grad_density(basis, denab, dens, off)

!   temporary storage for 1e gradient
    ALLOCATE(de_priv, mold=de)
    de_priv = 0.0d0

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       shi, shj, cntp, &
!$omp       de_atom &
!$omp   ) &
!$omp   reduction(+:de_priv)

    CALL cntp%alloc(basis)

!$omp do schedule(dynamic)
!   I shell
    DO ii = 1, basis%nshell


        CALL shi%fetch_by_id(basis, ii)
        de_atom = 0.0

!       J shell
        DO jj = 1, basis%nshell

            CALL shj%fetch_by_id(basis, jj)

            CALL cntp%shell_pair(basis, shi, shj, tol, dup=.false.)
            IF (cntp%numpairs==0) CYCLE

            CALL comp_kinetic_der1(cntp, dens(off(ii):, off(jj):), de_atom)

        END DO

        de_priv(:,shi%atid) = de_priv(:,shi%atid) + 2*de_atom(:)

    END DO
!$omp end do
!   End of shell loops
!$omp end parallel

    de = de + de_priv
    DEALLOCATE(de_priv)

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Overlap second-derivative contribution to the Cartesian Hessian.
!> @details Accumulates  sum_uv M_uv d2 S_uv / dR_a dR_b  into the (3N,3N)
!>   Hessian, where M is the matrix passed in packed form (the energy-weighted
!>   density W for the HF Hessian). Mirrors grad_ee_overlap: it loops ordered
!>   shell pairs, evaluates the bra-center second derivative (comp_overlap_der2)
!>   with the same factor-2 convention as the gradient, and uses translational
!>   invariance of the two-center integral (d/dB = -d/dA) for the cross block.
 SUBROUTINE hess_ee_overlap(basis, denab, hess, logtol)
    implicit none
    type(basis_set), intent(inout) :: basis
    REAL(kind=dp), INTENT(INOUT) :: denab(:)
    real(kind=dp), intent(inout) :: hess(:,:)
    real(kind=dp), optional :: logtol

    INTEGER :: ii, jj, a, b, ai, bi
    REAL(kind=dp) :: tol, de2(3,3)
    LOGICAL :: norm
    REAL(kind=dp), ALLOCATABLE :: hess_priv(:,:), dens(:,:)
    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    if (present(logtol)) then
        tol = logtol
    else
        tol = tol_default
    end if

    norm = .true.
    allocate(dens(basis%nbf,basis%nbf), source=0.0d0)
    call unpack_matrix(denab, dens)
    IF (norm) CALL bas_norm_matrix(dens, basis%bfnrm, basis%nbf)

    allocate(hess_priv, mold=hess)
    hess_priv = 0.0d0

!$omp parallel &
!$omp   private(ii, jj, a, b, ai, bi, shi, shj, cntp, de2) &
!$omp   reduction(+:hess_priv)
    CALL cntp%alloc(basis)
!$omp do schedule(dynamic)
    DO ii = 1, basis%nshell
        CALL shi%fetch_by_id(basis, ii)
        DO jj = 1, basis%nshell
            CALL shj%fetch_by_id(basis, jj)
            CALL cntp%shell_pair(basis, shi, shj, tol, dup=.false.)
            IF (cntp%numpairs==0) CYCLE
            de2 = 0.0d0
            CALL comp_overlap_der2(cntp, &
                dens(basis%ao_offset(ii):, basis%ao_offset(jj):), de2)
            ! d/dR_C of  G_A = 2*sum_bra-on-A ...  : C=A gives +2*de2, C=B gives -2*de2
            DO a = 1, 3
                ai = 3*(shi%atid-1) + a
                DO b = 1, 3
                    hess_priv(3*(shi%atid-1)+b, ai) = &
                        hess_priv(3*(shi%atid-1)+b, ai) + 2*de2(b,a)
                    hess_priv(3*(shi%atid-1)+b, 3*(shj%atid-1)+a) = &
                        hess_priv(3*(shi%atid-1)+b, 3*(shj%atid-1)+a) - 2*de2(b,a)
                END DO
            END DO
        END DO
    END DO
!$omp end do
!$omp end parallel

    hess = hess + hess_priv
    DEALLOCATE(hess_priv)
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Kinetic-energy second-derivative contribution to the Cartesian Hessian.
!> @details Accumulates  sum_uv M_uv d2 T_uv / dR_a dR_b  into the (3N,3N)
!>   Hessian. Same structure and conventions as hess_ee_overlap, using
!>   comp_kinetic_der2.
 SUBROUTINE hess_ee_kinetic(basis, denab, hess, logtol)
    implicit none
    type(basis_set), intent(inout) :: basis
    REAL(kind=dp), INTENT(INOUT) :: denab(:)
    real(kind=dp), intent(inout) :: hess(:,:)
    real(kind=dp), optional :: logtol

    INTEGER :: ii, jj, a, b, ai
    REAL(kind=dp) :: tol, de2(3,3)
    LOGICAL :: norm
    REAL(kind=dp), ALLOCATABLE :: hess_priv(:,:), dens(:,:)
    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    if (present(logtol)) then
        tol = logtol
    else
        tol = tol_default
    end if

    norm = .true.
    allocate(dens(basis%nbf,basis%nbf), source=0.0d0)
    call unpack_matrix(denab, dens)
    IF (norm) CALL bas_norm_matrix(dens, basis%bfnrm, basis%nbf)

    allocate(hess_priv, mold=hess)
    hess_priv = 0.0d0

!$omp parallel &
!$omp   private(ii, jj, a, b, ai, shi, shj, cntp, de2) &
!$omp   reduction(+:hess_priv)
    CALL cntp%alloc(basis)
!$omp do schedule(dynamic)
    DO ii = 1, basis%nshell
        CALL shi%fetch_by_id(basis, ii)
        DO jj = 1, basis%nshell
            CALL shj%fetch_by_id(basis, jj)
            CALL cntp%shell_pair(basis, shi, shj, tol, dup=.false.)
            IF (cntp%numpairs==0) CYCLE
            de2 = 0.0d0
            CALL comp_kinetic_der2(cntp, &
                dens(basis%ao_offset(ii):, basis%ao_offset(jj):), de2)
            DO a = 1, 3
                ai = 3*(shi%atid-1) + a
                DO b = 1, 3
                    hess_priv(3*(shi%atid-1)+b, ai) = &
                        hess_priv(3*(shi%atid-1)+b, ai) + 2*de2(b,a)
                    hess_priv(3*(shi%atid-1)+b, 3*(shj%atid-1)+a) = &
                        hess_priv(3*(shi%atid-1)+b, 3*(shj%atid-1)+a) - 2*de2(b,a)
                END DO
            END DO
        END DO
    END DO
!$omp end do
!$omp end parallel

    hess = hess + hess_priv
    DEALLOCATE(hess_priv)
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Nuclear-attraction second-derivative contribution to the Cartesian
!>   Hessian (electron-nucleus 1e term), on a fully bra-validated derivative path.
!> @warning WORK IN PROGRESS, not yet validated. The mixed bra-charge block
!>   (p_AC, d2/dA dC) computed from the bra derivative of the Hellmann-Feynman
!>   term currently disagrees with finite differences (see hess1_selftest,
!>   nucattr WIP line). This routine is not wired into any production path; the
!>   native hf_hessian kernel remains guarded.
!> @details For each ordered shell pair (bra atom A, ket atom B) and nucleus C,
!>   comp_coulomb_der2_braC is called twice -- as (A,B) giving the bra blocks
!>   p_AA=d2/dA2, p_AC=d2/dA dC, and swapped as (B,A) giving p_BB=d2/dB2,
!>   p_BC=d2/dB dC. Translational invariance d/dC = -(d/dA + d/dB) then yields
!>     p_AB = -(p_AA + p_AC),  p_CC = p_AA + p_AB + p_AB^T + p_BB,
!>   and all nine atom blocks are scattered into the (3N,3N) Hessian. Summing
!>   over ordered pairs reproduces the true d2E/dRdR (no extra symmetry factor).
 SUBROUTINE hess_en(basis, coord, zq, denab, hess, logtol, hess_cc)
    implicit none
    type(basis_set), intent(inout) :: basis
    real(kind=dp), contiguous, intent(in) :: coord(:,:), zq(:)
    REAL(kind=dp), INTENT(INOUT) :: denab(:)
    real(kind=dp), intent(inout) :: hess(:,:)
    real(kind=dp), optional :: logtol
    real(kind=dp), optional, intent(inout) :: hess_cc(:,:)

    INTEGER :: ii, jj, ic, a, b, n, nat, pat, qat
    REAL(kind=dp) :: tol
    REAL(kind=dp) :: p_AA(3,3), p_AC(3,3), p_BB(3,3), p_BC(3,3)
    REAL(kind=dp) :: bAB(3,3), blocks(3,3,9)
    INTEGER :: atP(9), atQ(9)
    LOGICAL :: norm
    REAL(kind=dp), ALLOCATABLE :: hess_priv(:,:), dens(:,:)
    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cab, cba

    if (present(logtol)) then
        tol = logtol
    else
        tol = tol_default
    end if

    nat = ubound(coord, 2)

    norm = .true.
    allocate(dens(basis%nbf,basis%nbf), source=0.0d0)
    call unpack_matrix(denab, dens)
    IF (norm) CALL bas_norm_matrix(dens, basis%bfnrm, basis%nbf)

    allocate(hess_priv, mold=hess)
    hess_priv = 0.0d0

!$omp parallel &
!$omp   private(ii, jj, ic, a, b, n, pat, qat, shi, shj, cab, cba, &
!$omp           p_AA, p_AC, p_BB, p_BC, bAB, blocks, atP, atQ) &
!$omp   reduction(+:hess_priv)
    CALL cab%alloc(basis)
    CALL cba%alloc(basis)
!$omp do schedule(dynamic)
    DO ii = 1, basis%nshell
        CALL shi%fetch_by_id(basis, ii)
        DO jj = 1, basis%nshell
            CALL shj%fetch_by_id(basis, jj)
            CALL cab%shell_pair(basis, shi, shj, tol, dup=.false.)
            IF (cab%numpairs==0) CYCLE
            CALL cba%shell_pair(basis, shj, shi, tol, dup=.false.)
            DO ic = 1, nat
                p_AA = 0.0d0; p_AC = 0.0d0; p_BB = 0.0d0; p_BC = 0.0d0
                CALL comp_coulomb_der2_braC(cab, coord(:,ic), -zq(ic), &
                    dens(basis%ao_offset(ii):, basis%ao_offset(jj):), p_AA, p_AC)
                CALL comp_coulomb_der2_braC(cba, coord(:,ic), -zq(ic), &
                    dens(basis%ao_offset(jj):, basis%ao_offset(ii):), p_BB, p_BC)

                bAB = -(p_AA + p_AC)

                atP(1)=shi%atid; atQ(1)=shi%atid; blocks(:,:,1)=p_AA
                atP(2)=shj%atid; atQ(2)=shj%atid; blocks(:,:,2)=p_BB
                atP(3)=shi%atid; atQ(3)=shj%atid; blocks(:,:,3)=bAB
                atP(4)=shj%atid; atQ(4)=shi%atid; blocks(:,:,4)=transpose(bAB)
                atP(5)=shi%atid; atQ(5)=ic;       blocks(:,:,5)=p_AC
                atP(6)=ic;       atQ(6)=shi%atid; blocks(:,:,6)=transpose(p_AC)
                atP(7)=shj%atid; atQ(7)=ic;       blocks(:,:,7)=p_BC
                atP(8)=ic;       atQ(8)=shj%atid; blocks(:,:,8)=transpose(p_BC)
                atP(9)=ic;       atQ(9)=ic;       blocks(:,:,9)=p_AA+bAB+transpose(bAB)+p_BB

                DO n = 1, 9
                    pat = atP(n); qat = atQ(n)
                    DO a = 1, 3
                        DO b = 1, 3
                            hess_priv(3*(pat-1)+a, 3*(qat-1)+b) = &
                                hess_priv(3*(pat-1)+a, 3*(qat-1)+b) + blocks(a,b,n)
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO
!$omp end do
!$omp end parallel

    hess = hess + hess_priv
    DEALLOCATE(hess_priv)
    if (present(hess_cc)) then
        call cab%alloc(basis)
        call cba%alloc(basis)
        do ii = 1, basis%nshell
            call shi%fetch_by_id(basis, ii)
            do jj = 1, basis%nshell
                call shj%fetch_by_id(basis, jj)
                call cab%shell_pair(basis, shi, shj, tol, dup=.false.)
                if (cab%numpairs == 0) cycle
                call cba%shell_pair(basis, shj, shi, tol, dup=.false.)
                do ic = 1, nat
                    p_AA = 0.0d0; p_AC = 0.0d0; p_BB = 0.0d0; p_BC = 0.0d0
                    call comp_coulomb_der2_braC(cab, coord(:,ic), -zq(ic), &
                        dens(basis%ao_offset(ii):, basis%ao_offset(jj):), p_AA, p_AC)
                    call comp_coulomb_der2_braC(cba, coord(:,ic), -zq(ic), &
                        dens(basis%ao_offset(jj):, basis%ao_offset(ii):), p_BB, p_BC)
                    bAB = -(p_AA + p_AC)
                    do a = 1, 3
                        do b = 1, 3
                            hess_cc(3*(ic-1)+a, 3*(ic-1)+b) = &
                                hess_cc(3*(ic-1)+a, 3*(ic-1)+b) + &
                                p_AA(a,b) + bAB(a,b) + bAB(b,a) + p_BB(a,b)
                        end do
                    end do
                end do
            end do
        end do
    end if
    DEALLOCATE(dens)
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Build the AO overlap first-derivative matrices dS_uv/dR for every
!>   nuclear coordinate (a CPHF right-hand-side building block).
!> @details Returns dS(nbf, nbf, 3, natom) where dS(:,:,c,A) = dS/dR_{A,c}. For
!>   each ordered shell pair the bra-center derivative block is scattered to the
!>   bra atom (+) and, by translational invariance of the two-center overlap
!>   (d/dB = -d/dA), to the ket atom (-). The integrals are in the same
!>   unnormalized convention as the stored overlap matrix (basis normalization
!>   is applied to the contracting density, as in grad_ee_overlap), so
!>   sum_uv (bfnrm_u bfnrm_v M_uv) dS(u,v,c,A) reproduces grad_ee_overlap(M).
 SUBROUTINE der_overlap_matrix(basis, dS, logtol)
    implicit none
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(out) :: dS(:,:,:,:)   ! (nbf, nbf, 3, natom)
    real(kind=dp), optional :: logtol

    INTEGER :: ii, jj, c, i, j, gi, gj, A_at, B_at, oi, oj
    REAL(kind=dp) :: tol
    REAL(kind=dp), ALLOCATABLE :: dblk(:,:,:)
    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    if (present(logtol)) then
        tol = logtol
    else
        tol = tol_default
    end if

    dS = 0.0d0

    CALL cntp%alloc(basis)
    DO ii = 1, basis%nshell
        CALL shi%fetch_by_id(basis, ii)
        A_at = shi%atid
        oi = basis%ao_offset(ii) - 1
        DO jj = 1, basis%nshell
            CALL shj%fetch_by_id(basis, jj)
            B_at = shj%atid
            oj = basis%ao_offset(jj) - 1
            CALL cntp%shell_pair(basis, shi, shj, tol, dup=.false.)
            IF (cntp%numpairs==0) CYCLE
            allocate(dblk(cntp%inao, cntp%jnao, 3), source=0.0d0)
            CALL comp_overlap_der1_block(cntp, dblk)
            DO c = 1, 3
                DO i = 1, cntp%inao
                    gi = oi + i
                    DO j = 1, cntp%jnao
                        gj = oj + j
                        dS(gi, gj, c, A_at) = dS(gi, gj, c, A_at) + dblk(i,j,c)
                        dS(gi, gj, c, B_at) = dS(gi, gj, c, B_at) - dblk(i,j,c)
                    END DO
                END DO
            END DO
            deallocate(dblk)
        END DO
    END DO
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Build the AO nuclear-attraction first-derivative matrices dV_uv/dR for
!>   every nuclear coordinate (a CPHF right-hand-side building block).
!> @details V_uv = sum_C (-Z_C) <u|1/|r-C||v> depends on three centers: bra atom
!>   A, ket atom B, and each charge atom C. For every ordered shell pair and
!>   nucleus C, the bra-center derivative (comp_coulomb_der1_block) is scattered
!>   to A, the charge-center derivative (comp_coulomb_helfeyder1_block) to C, and
!>   the ket-center derivative follows from translational invariance of the
!>   integral, d/dB = -(d/dA + d/dC), scattered to B. Contracting dV with the
!>   normalized density reproduces grad_en_pulay + grad_en_hellman_feynman.
 SUBROUTINE der_nucattr_matrix(basis, coord, zq, dV, logtol)
    implicit none
    type(basis_set), intent(inout) :: basis
    real(kind=dp), contiguous, intent(in) :: coord(:,:), zq(:)
    real(kind=dp), intent(out) :: dV(:,:,:,:)   ! (nbf, nbf, 3, natom)
    real(kind=dp), optional :: logtol

    INTEGER :: ii, jj, ic, c, i, j, gi, gj, A_at, B_at, oi, oj, nat
    REAL(kind=dp) :: tol, dba, dbc
    REAL(kind=dp), ALLOCATABLE :: dA(:,:,:), dC(:,:,:)
    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    if (present(logtol)) then
        tol = logtol
    else
        tol = tol_default
    end if

    nat = ubound(coord, 2)
    dV = 0.0d0

    CALL cntp%alloc(basis)
    DO ii = 1, basis%nshell
        CALL shi%fetch_by_id(basis, ii)
        A_at = shi%atid
        oi = basis%ao_offset(ii) - 1
        DO jj = 1, basis%nshell
            CALL shj%fetch_by_id(basis, jj)
            B_at = shj%atid
            oj = basis%ao_offset(jj) - 1
            CALL cntp%shell_pair(basis, shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE
            allocate(dA(cntp%inao, cntp%jnao, 3), dC(cntp%inao, cntp%jnao, 3))
            DO ic = 1, nat
                dA = 0.0d0; dC = 0.0d0
                CALL comp_coulomb_der1_block(cntp, coord(:,ic), -zq(ic), dA)
                CALL comp_coulomb_helfeyder1_block(cntp, coord(:,ic), -zq(ic), dC)
                DO c = 1, 3
                    DO i = 1, cntp%inao
                        gi = oi + i
                        DO j = 1, cntp%jnao
                            gj = oj + j
                            dba = dA(i,j,c); dbc = dC(i,j,c)
                            dV(gi, gj, c, A_at) = dV(gi, gj, c, A_at) + dba
                            dV(gi, gj, c, ic)   = dV(gi, gj, c, ic)   + dbc
                            dV(gi, gj, c, B_at) = dV(gi, gj, c, B_at) - (dba + dbc)
                        END DO
                    END DO
                END DO
            END DO
            deallocate(dA, dC)
        END DO
    END DO
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Build the AO kinetic-energy first-derivative matrices dT_uv/dR for
!>   every nuclear coordinate (a CPHF right-hand-side building block). Same
!>   structure and conventions as der_overlap_matrix, using comp_kinetic_der1_block.
 SUBROUTINE der_kinetic_matrix(basis, dT, logtol)
    implicit none
    type(basis_set), intent(inout) :: basis
    real(kind=dp), intent(out) :: dT(:,:,:,:)   ! (nbf, nbf, 3, natom)
    real(kind=dp), optional :: logtol

    INTEGER :: ii, jj, c, i, j, gi, gj, A_at, B_at, oi, oj
    REAL(kind=dp) :: tol
    REAL(kind=dp), ALLOCATABLE :: dblk(:,:,:)
    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    if (present(logtol)) then
        tol = logtol
    else
        tol = tol_default
    end if

    dT = 0.0d0

    CALL cntp%alloc(basis)
    DO ii = 1, basis%nshell
        CALL shi%fetch_by_id(basis, ii)
        A_at = shi%atid
        oi = basis%ao_offset(ii) - 1
        DO jj = 1, basis%nshell
            CALL shj%fetch_by_id(basis, jj)
            B_at = shj%atid
            oj = basis%ao_offset(jj) - 1
            CALL cntp%shell_pair(basis, shi, shj, tol, dup=.false.)
            IF (cntp%numpairs==0) CYCLE
            allocate(dblk(cntp%inao, cntp%jnao, 3), source=0.0d0)
            CALL comp_kinetic_der1_block(cntp, dblk)
            DO c = 1, 3
                DO i = 1, cntp%inao
                    gi = oi + i
                    DO j = 1, cntp%jnao
                        gj = oj + j
                        dT(gi, gj, c, A_at) = dT(gi, gj, c, A_at) + dblk(i,j,c)
                        dT(gi, gj, c, B_at) = dT(gi, gj, c, B_at) - dblk(i,j,c)
                    END DO
                END DO
            END DO
            deallocate(dblk)
        END DO
    END DO
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Basis function derivative contributions to gradient
!> @details Compute derivative integrals of type <ii'|h|jj> = <ii'|t+v|jj>
!> @note No relativistic methods available
!
!> @author Vladimir Mironov
!
!   REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   denab   density matrix in packed format, remains unchanged on return
 SUBROUTINE grad_en_pulay(basis, coord, zq, denab, de, logtol)

    REAL(kind=dp), INTENT(INOUT) :: denab(:)
    type(basis_set), intent(inout) :: basis
    real(kind=dp), contiguous, intent(in) :: coord(:,:), zq(:)

    REAL(kind=dp) :: de(:,:)


    INTEGER :: ii, jj, ic

    REAL(kind=dp), optional :: logtol

    REAL(kind=dp) :: de1(3)
    REAL(kind=dp), ALLOCATABLE :: de_priv(:,:), dens(:,:)
    INTEGER, ALLOCATABLE :: off(:)

    REAL(kind=dp) :: dernuc(3), tol
    LOGICAL :: out, dbg, norm

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    INTEGER :: nat

    dbg = .false.
    out = .false.

    IF (dbg) WRITE(iw,'(/10X,38(1H-)/10X,"GRADIENT INCLUDING AO DERIVATIVE TERMS"/10X,38(1H-))')

    nat = ubound(de, 2)

    if (present(logtol)) then
        tol = logtol
    else
        tol = tol_default
    end if

    call prepare_grad_density(basis, denab, dens, off)

!   temporary storage for 1e gradient
    ALLOCATE(de_priv, mold=de)
    de_priv = 0.0d0

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, ic, &
!$omp       shi, shj, cntp, &
!$omp       dernuc, &
!$omp       de1 &
!$omp   ) &
!$omp   reduction(+:de_priv)

    CALL cntp%alloc(basis)

!$omp do schedule(dynamic)
!   I shell
    DO ii = 1, basis%nshell


        CALL shi%fetch_by_id(basis, ii)
        de1 = 0.0

!       J shell
        DO jj = 1, basis%nshell

            CALL shj%fetch_by_id(basis, jj)

            CALL cntp%shell_pair(basis, shi, shj, tol, dup=.false.)
            IF (cntp%numpairs==0) CYCLE

!           Nuclear attraction derivative
            DO ic = 1, nat
                CALL comp_coulomb_der1(cntp, coord(:,ic), -zq(ic), dens(off(ii):, off(jj):), dernuc)
                de1 = de1 + 2*dernuc(1:3)
            END DO
!           End of primitive loops

#ifdef DEBUG
            IF (dbg) THEN
               WRITE(iw,'(1X,"TVDER: SHELLS II,JJ=",2I5)') ii,jj
               CALL print_grd(zq, de)
            END IF
#endif

        END DO

        de_priv(:,shi%atid) = de_priv(:,shi%atid) + de1(:)

    END DO
!$omp end do
!   End of shell loops
!$omp end parallel

    de = de + de_priv
    DEALLOCATE(de_priv)

    IF (out) THEN
       WRITE(iw,'(/10X,38(1H-)/10X,"GRADIENT INCLUDING AO DERIVATIVE TERMS"/10X,38(1H-))')
       CALL print_grd(zq,de)
    END IF

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Basis function derivative contributions to gradient from external charges
!> @details Compute derivative integrals of type <ii'|h|jj> = <ii'|t+v|jj>
!> @note No relativistic methods available
!
!> @author Vladimir Mironov
!
!   REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   denab   density matrix in packed format, remains unchanged on return
 SUBROUTINE omp_extder(basis, de, denab, ext_charges, de_mm, logtol, alpha)

    type(basis_set), intent(inout) :: basis
    REAL(kind=dp), INTENT(INOUT) :: denab(:)

    REAL(kind=dp) :: de(:,:)
    REAL(kind=dp) :: ext_charges(:,:)
    REAL(kind=dp), optional :: de_mm(:,:)
    REAL(kind=dp), optional :: logtol
    REAL(kind=dp), optional :: alpha


    INTEGER :: ii, jj, ic

    REAL(kind=dp) :: tol, znuc

    REAL(kind=dp) :: de1(3), de2(3)

    REAL(kind=dp) :: dernuc(3), cxyz(3)
    LOGICAL :: out, dbg, norm

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    REAL(kind=dp), ALLOCATABLE :: de_priv(:,:), dens(:,:)
    INTEGER, ALLOCATABLE :: off(:)

    dbg = .false.
    out = .false.

    IF (dbg) WRITE(iw,'(/10X,38(1H-)/10X,"GRADIENT INCLUDING AO DERIVATIVE TERMS"/10X,38(1H-))')

    if (present(logtol)) then
        tol = logtol
    else
        tol = tol_default
    end if

    norm = .true.

    !IF (dbg) write(iw,'(A,5I10)') "OMP 1E GRD (EXTERNAL CHARGES)", basis%nshell, natfmo, nps, nchmat, nqmmatm

    call prepare_grad_density(basis, denab, dens, off)

!   temporary storage for 1e gradient
    ALLOCATE(de_priv, mold=de)
    de_priv(:,:) = 0.0

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, ic, &
!$omp       shi, shj, cntp, &
!$omp       dernuc, &
!$omp       znuc, cxyz, &
!$omp       de1, de2 &
!$omp   ) &
!$omp   reduction(+:de_priv)

    CALL cntp%alloc(basis)
!   I shell
    DO ii = 1, basis%nshell

        CALL shi%fetch_by_id(basis, ii)
        de1(:) = 0.0
        de2(:) = 0.0

!       J shell
        DO jj = 1, basis%nshell

            CALL shj%fetch_by_id(basis, jj)

            CALL cntp%shell_pair(basis, shi, shj, tol, dup=.false.)
            IF (cntp%numpairs==0) CYCLE

!           External charges (QM/MM)
!$omp do schedule(dynamic,4)
            DO ic = 1, ubound(ext_charges, 2)
                cxyz =  ext_charges(1:3,ic)
                znuc = -ext_charges(4,ic)
!                IF ( doscr &
!                    .AND. (znuc**2 < scrthr**2 * sum((cxyz-shi%r)**2)) ) CYCLE

                CALL comp_coulomb_der1(cntp, cxyz, znuc, dens(off(ii):, off(jj):), dernuc)

                ! Ewald screening
                IF (present(alpha)) THEN
                    CALL comp_ewaldlr_der1(cntp, cxyz, -znuc, dens(basis%ao_offset(ii):, basis%ao_offset(jj):), alpha, dernuc)
                END IF

                de1 = de1 + 2*dernuc(1:3)

                ! Add gradient contribution to MM atoms?
                if (present(de_mm)) then
                  de_mm(:,ic) = de_mm(:,ic) - 2*dernuc(1:3)
                end if
            END DO
!$omp end do nowait

        END DO

        de_priv(1:3,shi%atid) = de_priv(1:3,shi%atid) + de1(1:3)
    END DO
!   End of shell loops

!$omp end parallel

    de = de + de_priv
    DEALLOCATE(de_priv)

!    IF (out) THEN
!      WRITE(iw,'(/10X,38(1H-)/10X,"GRADIENT INCLUDING AO DERIVATIVE TERMS"/10X,38(1H-))')
!      CALL print_grd(zq,de)
!    END IF

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Hellmann-Feynman force
!> @details Compute derivative contributions due to the Hamiltonian
!>   operator change w.r.t. shifts of nuclei. The contribution
!>   of the form <i|T'+V'|j> is evaluated by Gauss-Rys quadrature.
!>   This version handles spdfg and L shells.
!> @note No relativistic methods available
!
!> @author Vladimir Mironov
!
!   REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   denab   density matrix in packed format, remains unchanged on return
 SUBROUTINE grad_en_hellman_feynman(basis, coord, zq, denab, de, logtol)

!$  use omp_lib, ONLY: omp_get_max_threads

    type(basis_set), intent(inout) :: basis
    real(kind=dp), contiguous, intent(in) :: coord(:,:), zq(:)
    REAL(kind=dp), INTENT(INOUT) :: denab(:)
    REAL(kind=dp) :: de(:,:)

    INTEGER :: ii, jj, ic

    REAL(kind=dp), optional :: logtol

    REAL(kind=dp) :: dernuc(3), tol
    LOGICAL :: out, dbg, norm

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    REAL(kind=dp), ALLOCATABLE :: de_priv(:,:), dens(:,:)
    INTEGER, ALLOCATABLE :: off(:)
    INTEGER :: nat

    dbg = .false.
    out = .false.

    IF (dbg) WRITE(iw,'(/10X,22("-")/10X,"HELLMANN-FEYNMAN FORCE"/10X,22("-"))')

    nat = ubound(de, 2)

    if (present(logtol)) then
        tol = logtol
    else
        tol = tol_default
    end if

    call prepare_grad_density(basis, denab, dens, off)

!   temporary storage for 1e gradient
    ALLOCATE(de_priv, mold=de)
    de_priv = 0.0d0

!   Initialize parallel
!$omp parallel &
!$omp   num_threads(min(omp_get_max_threads(), nat)) &
!$omp   reduction(+:de_priv) &
!$omp   private( &
!$omp       ii, jj, ic, &
!$omp       shi, shj, cntp, &
!$omp       dernuc &
!$omp   )

    CALL cntp%alloc(basis)

!   I shell
!$omp do schedule(dynamic)
    DO ii = 1, basis%nshell

        CALL shi%fetch_by_id(basis, ii)

!       J shell
        DO jj = 1, basis%nshell

            CALL shj%fetch_by_id(basis, jj)

            CALL cntp%shell_pair(basis, shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

!           Hellmann-Feynman term
atoms:      DO ic = 1, nat

                CALL comp_coulomb_helfeyder1(cntp, coord(:,ic), -zq(ic), dens(off(ii):, off(jj):), dernuc)

                de_priv(:,ic) = de_priv(:,ic) + dernuc(:3)

            END DO atoms

#ifdef DEBUG
            IF (dbg) THEN
               WRITE(iw,'(1X,"HELFEY: SHELLS II,JJ=",2I5)') ii,jj
               CALL print_grd(zq,de_priv)
            END IF
#endif

        END DO
    END DO
!$omp end do
!   End of shell loops
!$omp end parallel

    de(:,1:nat) = de(:,1:nat) + de_priv(:3,1:nat)

    IF (out) THEN
       WRITE(iw,'(/10X,22(1H-)/10X,"HELLMANN-FEYNMAN FORCE"/10X,22(1H-))')
       CALL print_grd(zq,de)
    END IF

    DEALLOCATE(de_priv)

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Gradient of nuclear repulsion energy
  subroutine grad_nn(atoms, ecp_el)
    implicit none
    type(atomic_structure), intent(inout) :: atoms
    integer, intent(in) :: ecp_el(:)

    integer :: k, l
    real(kind=dp) :: pkl(3), rkl3, de1(3)

    do k = 2, ubound(atoms%zn, 1)
        do l = 1, k-1
            if (k==l) cycle
            pkl = atoms%xyz(:,k)-atoms%xyz(:,l)
            rkl3 = norm2(pkl)**3
            de1 = -(atoms%zn(k)-ecp_el(k))*(atoms%zn(l)-ecp_el(l))*pkl/rkl3
            atoms%grad(:,k) = atoms%grad(:,k) + de1
            atoms%grad(:,l) = atoms%grad(:,l) - de1
        end do
    end do

  end subroutine grad_nn

!-------------------------------------------------------------------------------

!> @brief Nuclear-repulsion contribution to the Cartesian Hessian.
!> @details Accumulates the second derivatives of the nuclear-repulsion energy
!>   E_nn = sum_{k>l} Zk*Zl / r_kl into the (3N, 3N) Hessian in OpenQP
!>   atom-major coordinate order (x1,y1,z1,x2,...). For each pair (k,l) the
!>   3x3 block is  Zk*Zl * (3 p_a p_b / r^5 - delta_ab / r^3), with p = r_k-r_l;
!>   it is added to the (k,k) and (l,l) diagonal blocks and subtracted from the
!>   (k,l) and (l,k) off-diagonal blocks. Effective nuclear charges use the same
!>   (zn - ecp_el) convention as grad_nn. The result is added in place so the
!>   routine composes with the electronic Hessian terms.
  subroutine hess_nn(atoms, ecp_el, hess)
    implicit none
    type(atomic_structure), intent(in) :: atoms
    integer, intent(in) :: ecp_el(:)
    real(kind=dp), intent(inout) :: hess(:,:)

    integer :: k, l, a, b, ka, lb
    real(kind=dp) :: pkl(3), r, r2, r3, zz, blk(3,3)

    do k = 2, ubound(atoms%zn, 1)
        do l = 1, k-1
            pkl = atoms%xyz(:,k) - atoms%xyz(:,l)
            r = norm2(pkl)
            r2 = r*r
            r3 = r*r2
            zz = (atoms%zn(k)-ecp_el(k))*(atoms%zn(l)-ecp_el(l))
            do a = 1, 3
                do b = 1, 3
                    blk(b,a) = zz * 3.0_dp*pkl(b)*pkl(a) / (r2*r3)
                    if (a == b) blk(b,a) = blk(b,a) - zz/r3
                end do
            end do
            do a = 1, 3
                ka = 3*(k-1) + a
                lb = 3*(l-1) + a
                do b = 1, 3
                    hess(3*(k-1)+b, ka) = hess(3*(k-1)+b, ka) + blk(b,a)
                    hess(3*(l-1)+b, lb) = hess(3*(l-1)+b, lb) + blk(b,a)
                    hess(3*(k-1)+b, 3*(l-1)+a) = hess(3*(k-1)+b, 3*(l-1)+a) - blk(b,a)
                    hess(3*(l-1)+b, 3*(k-1)+a) = hess(3*(l-1)+b, 3*(k-1)+a) - blk(b,a)
                end do
            end do
        end do
    end do

  end subroutine hess_nn

!-------------------------------------------------------------------------------

!> @brief Effective core potential gradient
  subroutine grad_1e_ecp(infos,basis, coord, denab, de, logtol)
    use types, only: information
    use parallel, only: par_env_t

    type(information), target, intent(inout) :: infos
    type(par_env_t) :: pe
    REAL(kind=dp), INTENT(INOUT) :: denab(:)
    type(basis_set), intent(inout) :: basis
    real(kind=dp), contiguous, intent(in) :: coord(:,:)
    REAL(kind=dp) :: de(:,:)

    REAL(kind=dp), optional :: logtol

    call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    if (pe%rank == 0) then
        call add_ecpder(basis, coord, denab, de)
    end if

    call pe%bcast(de, size(de))

  end subroutine grad_1e_ecp

  !-------------------------------------------------------------------------------

end module grd1
