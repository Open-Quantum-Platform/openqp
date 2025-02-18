module grd1

   use io_constants, only: iw
   use precision, only: dp
   use types, only: information
   use atomic_structure_m, only: atomic_structure

   use basis_tools, only: basis_set, &
       bas_norm_matrix, bas_denorm_matrix

   use mod_1e_primitives, only: &
       comp_coulomb_der1, comp_coulomb_helfeyder1, comp_kinetic_der1, &
       comp_overlap_der1, &
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
   public grad_ee_overlap
   public grad_ee_kinetic
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

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    if (present(logtol)) then
        tol = logtol
    else
        tol = tol_default
    end if

    allocate(de_priv, mold=de)
    de_priv = 0.0d0

    norm = .true.
    allocate(dens(basis%nbf,basis%nbf), source=0.0d0)
    call unpack_matrix(denab, dens)

    IF (norm) THEN
        CALL bas_norm_matrix(dens, basis%bfnrm, basis%nbf)
    END IF

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

            CALL comp_overlap_der1(cntp, dens(basis%ao_offset(ii):, basis%ao_offset(jj):), de_atom)
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

    norm = .true.
    allocate(dens(basis%nbf,basis%nbf), source=0.0d0)
    call unpack_matrix(denab, dens)

    IF (norm) THEN
        CALL bas_norm_matrix(dens, basis%bfnrm, basis%nbf)
    END IF

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

            CALL comp_kinetic_der1(cntp, dens(basis%ao_offset(ii):, basis%ao_offset(jj):), de_atom)

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

    norm = .true.
    allocate(dens(basis%nbf,basis%nbf), source=0.0d0)
    call unpack_matrix(denab, dens)

    IF (norm) THEN
        CALL bas_norm_matrix(dens, basis%bfnrm, basis%nbf)
    END IF

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
                CALL comp_coulomb_der1(cntp, coord(:,ic), -zq(ic), dens(basis%ao_offset(ii):, basis%ao_offset(jj):), dernuc)
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

    norm = .true.
    allocate(dens(basis%nbf,basis%nbf), source=0.0d0)
    call unpack_matrix(denab, dens)

    IF (norm) THEN
        CALL bas_norm_matrix(dens, basis%bfnrm, basis%nbf)
    END IF

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

                CALL comp_coulomb_der1(cntp, cxyz, znuc, dens(basis%ao_offset(ii):, basis%ao_offset(jj):), dernuc)

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

    norm = .true.
    allocate(dens(basis%nbf,basis%nbf), source=0.0d0)
    call unpack_matrix(denab, dens)

    IF (norm) THEN
        CALL bas_norm_matrix(dens, basis%bfnrm, basis%nbf)
    END IF

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

                CALL comp_coulomb_helfeyder1(cntp, coord(:,ic), -zq(ic), dens(basis%ao_offset(ii):, basis%ao_offset(jj):), dernuc)

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
