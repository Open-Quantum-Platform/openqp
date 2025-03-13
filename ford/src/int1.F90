!#define DEBUG 1
!> @author  Vladimir Mironov
!
!> @brief This module contains subroutines for 1-electron integrals
!>  calculation.
!
!  REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
module int1

    use, intrinsic :: iso_fortran_env, only: real64
    use basis_tools, only: basis_set, &
            bas_norm_matrix, &
            bas_denorm_matrix
    use mod_1e_primitives, only: &
        update_triang_matrix, &
        update_rectangular_matrix, &
        density_ordered, &
        comp_coulomb_int1_prim, &
        comp_ewaldlr_int1_prim, &
        comp_kin_ovl_int1_prim, &
        comp_lz_int1_prim, &
        MAX_EL_MOM, &
        comp_mult_int1_prim, &
        comp_allmult_int1_prim, &
        comp_coulpot_prim

    use mod_shell_tools, only: shell_t, shpair_t
    use messages, only: show_message, with_abort

    implicit none

!<  size of shell pair block (square of max.num. basis functions in max.ang.m.)
    integer, parameter :: blocksize = 28*28

    integer, parameter :: mult_bs(0:MAX_EL_MOM) = [1, 3, 6, 10]
    integer, parameter :: mult_all_bs(MAX_EL_MOM) = [3, 9, 19]

    interface int1_coul
       module procedure int1_coul_xyzc
       module procedure int1_coul_x_y_z_c
       module procedure int1_coul_xyz_c
    end interface

    private
    public omp_hst
    public multipole_integrals
    public electrostatic_potential
    public basis_overlap
    public overlap

contains
!> @brief Driver for conventional h, S, and T integrals
!
!> @details  Compute one electron integrals and core Hamiltonian,
!>  - S is evaluated by Gauss-Hermite quadrature,
!>  - T is an overlap with -2,0,+2 angular momentum shifts,
!>  - V is evaluated by Gauss-Rys quadrature, then \f$ h = T+V \f$
!>  Also, do \f$ L_z \f$ integrals if requested
!
!> @note Based on `HSANDT` subroutine from file `INT1.SRC`
!
!> @author Vladimir Mironov
!
!   REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   h       one-electron Hamiltonian matrix in packet format
!> @param[in,out]   s       packed matrix of overlap integrals
!> @param[in,out]   t       packed matrix of kinetic energy integrals
!> @param[in,out]   z       packed matrix of z-angular momentum (Lz) integrals
!> @param[in]       dbug    flag for debug output
 subroutine omp_hst(basis, coord, zq, h, s, t, z, debug, logtol, comm, usempi)

    use io_constants, only: iw
    use precision, only: dp
    use basis_tools, only: basis_set
    use printing, only: print_sym_labeled
    use ecp_tool, only: add_ecpint
    use parallel, only: par_env_t
    use iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32

    type(basis_set), intent(in) :: basis
    real(real64), contiguous, intent(in) :: coord(:,:), zq(:)
    real(real64), contiguous, intent(inout) :: h(:), s(:), t(:)
    real(real64), contiguous, optional, intent(inout) :: z(:)
    real(real64), optional, intent(in) :: logtol
    logical, optional, intent(in) :: debug
    integer :: ii

    real(real64) :: tol
    logical :: lzint, dbug

    integer :: nbf, nbf_tri
    type(par_env_t) :: pe

    integer(kind=int32) :: comm
    logical(c_bool), intent(in) :: usempi

    call pe%init(comm, usempi)


    lzint = present(z)
    dbug = .false.
    if (present(debug)) dbug = debug

    if (present(logtol)) then
        tol = logtol
    else
        tol = log(10.0_dp)*20
    end if

!   Exclude 1e potential in ESDIM, because density is used,
!   not point charges

    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2

!    Zero out all arrays
     s = 0.0
     t = 0.0
     h = 0.0
     if (lzint) z = 0.0

    call kin_ovl_ints(s, t, basis, tol)

    call nuc_ints(basis, coord(:,:), zq, h, tol)

!   Add effective core potential
    if(pe%rank == 0) then
        call add_ecpint(basis,coord(:,:),h)
    end if

    call pe%bcast(h, nbf_tri)

!    IF (exterior%num_chg/=0) THEN
!        SELECT CASE (pbc%method)
!        CASE (OQP_PBC_METHOD_EWALD)
!            IF (dbug) THEN
!                WRITE(*,*) 'Computing Erfc-attenuated Coulomb 1e-integrals'
!                WRITE(*,*) 'alpha=', pbc%alpha
!            END IF
!            CALL int1_coul_ext_chg_ewaldsr(h, basis, &
!                    exterior%num_chg, &
!                    exterior%chg(:,1), &
!                    exterior%chg(:,2), &
!                    exterior%chg(:,3), &
!                    exterior%chg(:,4), &
!                    tol, 1.0d-8, pbc%alpha)
!        CASE (OQP_PBC_METHOD_OFF)
!            IF (dbug) THEN
!                WRITE(*,*) 'Computing regular Coulomb 1e-integrals'
!            END IF
!            CALL int1_coul_ext_chg(h, basis, &
!                    exterior%num_chg, &
!                    exterior%chg(:,1), &
!                    exterior%chg(:,2), &
!                    exterior%chg(:,3), &
!                    exterior%chg(:,4), &
!                    tol, 1.0d-8)
!        CASE DEFAULT
!            WRITE (iw,*) 'Unknown PBC method selected'
!            CALL abrt
!        END SELECT
!    END IF

    if (lzint) call lzints(z, basis, tol)

!   Normalize 1-e integrals all at once
    call bas_norm_matrix(h, basis%bfnrm, nbf)
    call bas_norm_matrix(s, basis%bfnrm, nbf)
    call bas_norm_matrix(t, basis%bfnrm, nbf)
    if (lzint)  call bas_norm_matrix(z, basis%bfnrm, nbf)

!   Form one electron Hamiltonian
!   Hcore = Vne + Te
    h = h + t

!   Optional debug printout
    if (dbug) then
       write(iw,*) 'Overlap matrix (S)'
       call print_sym_labeled(s,nbf,basis)
       write(iw,*) 'Bare nucleus Hamiltonian integrals (H=T+V)'
       call print_sym_labeled(h,nbf,basis)
       write(iw,*) 'Kinetic energy integrals (T)'
       call print_sym_labeled(t,nbf,basis)
       if (lzint) then
          write(iw,*) 'Z-angular momentum integrals'
          call print_sym_labeled(z,nbf,basis)
       end if
    end if

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Driver for multipole integrals
!
!> @details  Compute one electron multipole integrals
!>  Integrals are evaluated by Gauss-Hermite quadrature,
!
!> @author Vladimir Mironov
!
!   REVISION HISTORY:
!> @date _Feb, 2023_ Initial release
!>
!> @param[in,out]   ints    integrals, packed format
!> @param[in]       dbug    flag for debug output
 subroutine multipole_integrals(basis, ints, r, mxmom, debug, logtol)

    use io_constants, only: iw
    use precision, only: dp
    use basis_tools, only: basis_set
    use printing, only: print_sym_labeled

    type(basis_set), intent(in) :: basis
    real(real64), contiguous, intent(inout) :: ints(:,:)
    real(real64), intent(in) :: r(:)
    integer, intent(in) :: mxmom
    real(real64), optional, intent(in) :: logtol
    logical, optional, intent(in) :: debug
    character(2) :: mxmom_str

    real(real64) :: tol
    logical :: dbug

    character(len=*), parameter :: labels(19) = [&
        'X  ', 'Y  ', 'Z  ', &
        'XX ', 'YY ', 'ZZ ', 'XY ', 'XZ ', 'YZ ', &
        'XXX', 'YYY', 'ZZZ', &
        'XXY', 'XXZ', &
        'YYX', 'YYZ', &
        'ZZX', 'ZZY', &
        'XYZ' &
        ]

    integer :: nbf
    integer :: i

    if (mxmom > 3) then
      write(mxmom_str,'(I2)') MAX_EL_MOM
      call show_message('Maximum order of multipole integrals is'//mxmom_str, with_abort)
    end if

    if (ubound(ints,2) < mult_all_bs(mxmom)) then
      write(iw,*) 'Insufficient space for multipole moment integrals: [', ubound(ints), ']'
      write(mxmom_str,'(I2)') mult_all_bs(MAX_EL_MOM)
      call show_message('Required:'//mxmom_str, with_abort)
    end if

    dbug = .false.
    if (present(debug)) dbug = debug

    tol = log(10.0_dp)*20
    if (present(logtol)) tol = logtol
    nbf = basis%nbf

!   Zero out all arrays
    ints = 0.0

    call mult_all_ints(ints, mxmom, r, basis, tol)

!   Normalize 1-e integrals all at once
    do i = 1, mult_all_bs(mxmom)
      call bas_norm_matrix(ints(:,i), basis%bfnrm, nbf)
    end do

!   Optional debug printout
    if (dbug) then
       do i = 1, mult_all_bs(mxmom)
         write(iw,*) 'Multipole moment integrals ('//trim(labels(i))//')'
         call print_sym_labeled(ints(:,i),nbf,basis)
       end do
    end if

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute electronic contribution to electrostatic potential on a grid
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2023_ Initial release
!>
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       x       array of X grid pts
!> @param[in]       y       array of Y grid pts
!> @param[in]       z       array of Z grid pts
!> @param[in]       wt      array of grid weights
!> @param[in]       d       density matrix
!> @param[in]       tol     1-e exponential prefactor tolerance
!> @param[out]      pot     electrostatic potential on a grid
 subroutine electrostatic_potential(basis, x, y, z, wt, d, pot, logtol)

    use precision, only: dp
    implicit none
    type(basis_set), intent(inout)          :: basis
    real(real64), contiguous, intent(in)    :: x(:), y(:), z(:), wt(:)
    real(real64), contiguous, intent(inout) :: d(:)
    real(real64), contiguous, intent(out)   :: pot(:)
    real(real64), optional, intent(in)      :: logtol
    real(real64) :: tol

    call bas_norm_matrix(d, basis%bfnrm, basis%nbf)

    tol = log(10.0_dp)*20
    if (present(logtol)) tol = logtol

    call int1_el_pot(basis, x, y, z, d, pot, tol)
    pot = pot*wt

    call bas_denorm_matrix(d, basis%bfnrm, basis%nbf)

 end subroutine


!-------------------------------------------------------------------------------

!> @brief Compute overlap matrix between two basis sets
!
!> @details Overlap integrals are computed using Gauss-Hermite quadrature formula
!
!> @author   Igor S. Gerasimov
!
!     REVISION HISTORY:
!> @date _Oct, 2022_ Initial release
!>
!> @param[in,out]   s       unpacked matrix of overlap integrals
!> @param[in]       basis1  basis w/ SP-shells separated
!> @param[in]       basis2  basis w/ SP-shells separated
!> @param[in]       tol     1-e exponential prefactor tolerance (should be ~tol_int*log(10.0_dp))
 SUBROUTINE basis_overlap(s, basis1, basis2, tol)

    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: s(:,:)
    TYPE(basis_set), INTENT(IN) :: basis1, basis2 ! basis without sp-shells
    REAL(REAL64),    INTENT(IN) :: tol

    INTEGER :: ii, jj

    LOGICAL, PARAMETER :: dokinetic = .false.

    REAL(REAL64), DIMENSION(BLOCKSIZE) :: sblk
    REAL(REAL64), DIMENSION(1) :: tblk ! should be zero, but...
!dir$ attributes align : 64 :: sblk
!dir$ attributes align : 64 :: tblk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    s = 0

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       sblk, tblk, &
!$omp       shi, shj, cntp &
!$omp   )

    CALL cntp%alloc2(basis1, basis2)

!   I shell
    DO ii = basis1%nshell, 1, -1

        CALL shi%fetch_by_id(basis1, ii)

!       J shell
!$omp do schedule(dynamic)
        DO jj = basis2%nshell, 1, -1

            CALL shj%fetch_by_id(basis2, jj)

            CALL cntp%shell_pair2(basis1, basis2, shi, shj, tol)
            IF (cntp%numpairs == 0) CYCLE

            sblk = 0.0

            CALL int1_kin_ovl(cntp, dokinetic, sblk, tblk)

            CALL update_rectangular_matrix(shi, shj, sblk, s)

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute overlap and integrals
!
!> @details Overlap integrals
!>  are computed using Gauss-Hermite quadrature formula
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
!>
!> @param[in,out]   s       packed matrix of overlap integrals
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       tol     1-e exponential prefactor tolerance
 SUBROUTINE overlap(s, basis, tol)

    REAL(REAL64), CONTIGUOUS,  INTENT(INOUT)  :: s(:)
    TYPE(basis_set), INTENT(IN)     :: basis ! basis without sp-shells
    REAL(REAL64),   INTENT(IN)     :: tol

    INTEGER :: &
        ii, jj

    logical, parameter :: dokinetic = .false.

    REAL(REAL64), DIMENSION(BLOCKSIZE) :: sblk
    REAL(REAL64), DIMENSION(1) :: tblk
!dir$ attributes align : 64 :: tblk, sblk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       sblk, &
!$omp       shi, shj, cntp &
!$omp   )

    CALL cntp%alloc(basis)

!   I shell
    DO ii = basis%nshell, 1, -1

        CALL shi%fetch_by_id(basis,ii)

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(basis,jj)

            CALL cntp%shell_pair(basis,shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            sblk = 0.0

            CALL int1_kin_ovl(cntp, dokinetic, sblk, tblk)

            CALL update_triang_matrix(shi, shj, sblk, s)

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute overlap and kinetic integrals
!
!> @details Overlap and electron kinetic energy integrals
!>  are computed using Gauss-Hermite quadrature formula
!>  Kinetic energy integrals are actually overlap integrals with +2, -2 angular
!>  momentum shifts
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   s       packed matrix of overlap integrals
!> @param[in,out]   t       packed matrix of kinetic energy integrals
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       tol     1-e exponential prefactor tolerance
 SUBROUTINE kin_ovl_ints(s, t, basis, tol)

    REAL(REAL64), CONTIGUOUS,  INTENT(INOUT)  :: s(:), t(:)
    TYPE(basis_set), INTENT(IN)     :: basis ! basis without sp-shells
    REAL(REAL64),   INTENT(IN)     :: tol

    INTEGER :: &
        ii, jj

    LOGICAL :: dokinetic

    REAL(REAL64), DIMENSION(BLOCKSIZE) :: tblk, sblk
!dir$ attributes align : 64 :: tblk, sblk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       sblk, tblk, &
!$omp       shi, shj, cntp, dokinetic &
!$omp   )

    dokinetic = .true.
    CALL cntp%alloc(basis)

!   I shell
    DO ii = basis%nshell, 1, -1

        CALL shi%fetch_by_id(basis,ii)

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(basis,jj)

            CALL cntp%shell_pair(basis,shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            sblk = 0.0
            tblk = 0.0

            CALL int1_kin_ovl(cntp, dokinetic, sblk, tblk)

            CALL update_triang_matrix(shi, shj, sblk, s)
            CALL update_triang_matrix(shi, shj, tblk, t)

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute multipole moment integrals
!> @author   Vladimir Mironov
!
 SUBROUTINE mult_all_ints(ints, mxmom, r, basis, tol)

    REAL(REAL64), CONTIGUOUS,  INTENT(INOUT)  :: ints(:,:)
    TYPE(basis_set), INTENT(IN)     :: basis ! basis without sp-shells
    real(real64), contiguous, intent(in) :: r(:)
    integer, intent(in) :: mxmom
    REAL(REAL64),   INTENT(IN)     :: tol

    INTEGER :: ii, jj, m

    REAL(REAL64), DIMENSION(BLOCKSIZE,19) :: blk
!dir$ attributes align : 64 :: blk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       blk, &
!$omp       shi, shj, cntp &
!$omp   )

    CALL cntp%alloc(basis)

!   I shell
    !DO ii = basis%nshell, 1, -1
    DO ii = 1, basis%nshell

        CALL shi%fetch_by_id(basis,ii)

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(basis,jj)

            CALL cntp%shell_pair(basis,shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            blk = 0.0

            CALL int1_allmul(cntp, r, mxmom, blk)

            do m = 1, mult_all_bs(mxmom)
              CALL update_triang_matrix(shi, shj, blk(:,m), ints(:,m))
            end do

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute multipole moment integrals
!> @author   Vladimir Mironov
!
 SUBROUTINE mult_ints(ints, mom, r, basis, tol)

    REAL(REAL64), CONTIGUOUS,  INTENT(INOUT)  :: ints(:,:)
    TYPE(basis_set), INTENT(IN)     :: basis ! basis without sp-shells
    real(real64), contiguous, intent(in) :: r(:)
    integer, intent(in) :: mom
    REAL(REAL64),   INTENT(IN)     :: tol

    INTEGER :: ii, jj, m

    REAL(REAL64), DIMENSION(BLOCKSIZE,10) :: blk
!dir$ attributes align : 64 :: blk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       blk, &
!$omp       shi, shj, cntp &
!$omp   )

    CALL cntp%alloc(basis)

!   I shell
    !DO ii = basis%nshell, 1, -1
    DO ii = 1, basis%nshell

        CALL shi%fetch_by_id(basis,ii)

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(basis,jj)

            CALL cntp%shell_pair(basis,shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            blk = 0.0

            CALL int1_mul(cntp, r, mom, blk)

            do m = 1, mult_bs(mom)
              CALL update_triang_matrix(shi, shj, blk(:,m), ints(:,m))
            end do

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops
 END SUBROUTINE


!-------------------------------------------------------------------------------

!> @brief Compute \f$ L_z \f$ integrals
!
!> @details \f$ L_z \f$ are actually overlap integrals with +1, -1 angular
!>  momentum shifts
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   z       packed matrix of Lz integrals
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       tol     1-e exponential prefactor tolerance
 SUBROUTINE lzints(z, basis, tol)

    TYPE(basis_set), INTENT(IN)  :: basis ! basis without sp-shells
    REAL(REAL64), CONTIGUOUS,  INTENT(OUT) :: z(:)
    REAL(REAL64),   INTENT(IN)  :: tol

    INTEGER :: &
        ii, jj

    REAL(REAL64), DIMENSION(BLOCKSIZE) :: zblk
!dir$ attributes align : 64 :: zblk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t)  :: cntp

    CALL cntp%alloc(basis)

!   I shell
    DO ii = 1, basis%nshell
        CALL shi%fetch_by_id(basis, ii)

!       J shell
        DO jj = 1, ii

            CALL shj%fetch_by_id(basis, jj)

            CALL cntp%shell_pair(basis, shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            zblk = 0.0
            CALL int1_lz(cntp, zblk)
            CALL update_triang_matrix(shi, shj, zblk, z)

        END DO
    END DO
!   End of shell loops
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute nuclear attraction integrals
!
!> @details Nuclear attaction integrals are computed using Gauss-Rys quadrature
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   h       core Hamiltonian matrix
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       tol     1-e exponential prefactor tolerance
 subroutine nuc_ints(basis, coord, zq, h, tol)

    type(basis_set), intent(in)     :: basis ! basis without sp-shells
    real(real64),   intent(in)     :: tol
    real(real64), contiguous,  intent(in)  :: coord(:,:), zq(:)
    real(real64), contiguous,  intent(inout)  :: h(:)

    integer :: nat, ii, jj

    real(real64), dimension(blocksize) :: vblk
!dir$ attributes align : 64 :: vblk

    type(shell_t) :: shi, shj
    type(shpair_t) :: cntp
    nat = ubound(zq, 1)

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       vblk, &
!$omp       shi, shj, cntp &
!$omp   )

    call cntp%alloc(basis)

!   I shell
    do ii = basis%nshell, 1, -1

        call shi%fetch_by_id(basis, ii)

!       J shell
!$omp do schedule(dynamic)
        do jj = 1, ii

            call shj%fetch_by_id(basis, jj)

            call cntp%shell_pair(basis, shi, shj, tol)
            if (cntp%numpairs==0) cycle

            vblk = 0.0d0

            call int1_coul(cntp, coord, zq, nat, 0.0d0, vblk)

            call update_triang_matrix(shi, shj, vblk, h)

        end do
!$omp end do nowait
    end do
!$omp end parallel
!   end of shell loops

 end subroutine

!-------------------------------------------------------------------------------

!> @brief General way to compute integrals of charge interaction, charge
!>  data are stored in a structure of arrays x(:),y(:),z(:),charge(:)
!
!> @details Electron-charge interaction integrals are computed using
!>  Gauss-Rys quadrature
!> @note This case is a variation of nuclear attraction case. They differ in
!>  data representation and in screening logic
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   h       packed matrix of one-electon Coulomb integrals
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       nat     number of atoms
!> @param[in]       x       array of X particle coordinates
!> @param[in]       y       array of Y particle coordinates
!> @param[in]       z       array of Z particle coordinates
!> @param[in]       chg     array of particle charges
!> @param[in]       tol     1-e exponential prefactor tolerance
!> @param[in]       chgtol  tolerance for particle charge
 SUBROUTINE int1_coul_ext_chg(h, basis, nat, x, y, z, chg, tol, chgtol)

    REAL(REAL64), CONTIGUOUS,  INTENT(INOUT)  :: h(:)
    TYPE(basis_set), INTENT(IN)     :: basis ! basis without sp-shells
    INTEGER,         INTENT(IN)     :: nat
    REAL(REAL64), CONTIGUOUS,  INTENT(IN)     :: x(:), y(:), z(:), chg(:)
    REAL(REAL64),   INTENT(IN)     :: tol, chgtol

    INTEGER :: &
        ii, jj

    REAL(REAL64), DIMENSION(BLOCKSIZE) :: vblk
!dir$ attributes align : 64 :: vblk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       vblk, &
!$omp       shi, shj, cntp &
!$omp   )
    CALL cntp%alloc(basis)

!   I shell
    DO ii = basis%nshell, 1, -1

        CALL shi%fetch_by_id(basis, ii)

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(basis, jj)

            CALL cntp%shell_pair(basis, shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            vblk = 0.0d0

            CALL int1_coul(cntp, x, y, z, chg, nat, chgtol, vblk)

            CALL update_triang_matrix(shi, shj, vblk, h)

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Ewald summation scheme, short-range part
!
!> @details Compute 1e integrals using modified Coulomb potential:
!>  \f$ \frac{Erfc(\omega^{1/2}|r-r_C|)}{|r-r_C|} \f$
!>  First, regular integrals are computed, then, long-range Ewald
!>  term is subtracted. Integrals are computed using Gauss-Rys quadrature.
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   h       packed matrix of one-electon Coulomb integrals
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       nat     number of atoms
!> @param[in]       x       array of X particle coordinates
!> @param[in]       y       array of Y particle coordinates
!> @param[in]       z       array of Z particle coordinates
!> @param[in]       chg     array of particle charges
!> @param[in]       tol     1-e exponential prefactor tolerance
!> @param[in]       chgtol  tolerance for particle charge
!> @param[in]       omega   Ewald splitting parameter
 SUBROUTINE int1_coul_ext_chg_ewaldsr(h, basis, nat, x, y, z, chg, tol, chgtol, omega)

    REAL(REAL64), CONTIGUOUS,  INTENT(INOUT)  :: h(:)
    TYPE(basis_set), INTENT(IN)     :: basis ! basis without sp-shells
    INTEGER,         INTENT(IN)     :: nat
    REAL(REAL64), CONTIGUOUS,  INTENT(IN)     :: x(:), y(:), z(:), chg(:)
    REAL(REAL64),   INTENT(IN)     :: tol, chgtol
    REAL(REAL64),   INTENT(IN)     :: omega

    INTEGER :: ii, jj

    REAL(REAL64), DIMENSION(BLOCKSIZE) :: vblk
!dir$ attributes align : 64 :: vblk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       vblk, &
!$omp       shi, shj, cntp &
!$omp   )
    CALL cntp%alloc(basis)

!   I shell
    DO ii = basis%nshell, 1, -1

        CALL shi%fetch_by_id(basis, ii)

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(basis, jj)

            CALL cntp%shell_pair(basis, shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            vblk = 0.0d0

            CALL int1_coul(cntp, x, y, z, chg, nat, chgtol, vblk)

!           Subtract long-range Ewald term
            CALL int1_ewald(cntp, x, y, z, chg, nat, chgtol, omega, vblk)

            CALL update_triang_matrix(shi, shj, vblk, h)

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute electronic contribution to electrostatic potential on a grid
!
!> @details Integrals are computed using Gauss-Rys quadrature
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2023_ Initial release
!>
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       x       array of X grid pts
!> @param[in]       y       array of Y grid pts
!> @param[in]       z       array of Z grid pts
!> @param[in]       d       density matrix
!> @param[in]       tol     1-e exponential prefactor tolerance
!> @param[out]      pot     electrostatic potential on a grid
 subroutine int1_el_pot(basis, x, y, z, d, pot, tol)

    type(basis_set), intent(in)           :: basis
    real(real64), contiguous, intent(in)  :: x(:), y(:), z(:), d(:)
    real(real64), contiguous, intent(out) :: pot(:)
    real(real64), intent(in)              :: tol

    integer :: npts, ii, jj, n

    real(real64), dimension(BLOCKSIZE) :: den
!dir$ attributes align : 64 :: den

    type(shell_t) :: shi, shj
    type(shpair_t) :: cntp

    npts = ubound(x,1)

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       den, &
!$omp       shi, shj, cntp &
!$omp   ) &
!$omp   reduction(+:pot)
    call cntp%alloc(basis)

!   i shell
    do ii = basis%nshell, 1, -1

      call shi%fetch_by_id(basis, ii)

!     j shell
!$omp do schedule(dynamic)
      do jj = 1, ii

        call shj%fetch_by_id(basis, jj)

        call cntp%shell_pair(basis, shi, shj, tol)
        if (cntp%numpairs==0) cycle

        call density_ordered(shi, shj, den, d)
        do n = 1, npts
          pot(n) = pot(n) + int1_epoten(cntp, x(n), y(n), z(n), den)
        end do

      end do
!$omp end do nowait
    end do
!$omp end parallel

 end subroutine

!--------------------------------------------------------------------------------
!       ONE-ELECTRON INTEGRALS CALCULATION (CONTRACTED SHELLS)
!--------------------------------------------------------------------------------

!> @brief Compute contracted block of kinetic energy and overlap 1e integrals
!> @param[in]       cntp        shell pair data
!> @param[in]       dokinetic   if `.FALSE.` compute only overlap integrals
!> @param[out]      sblk        block of overlap integrals
!> @param[out]      tblk        block of kinetic energy integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_kin_ovl(cntp, dokinetic, sblk, tblk)
!dir$ attributes inline :: int1_kin_ovl
    TYPE(shpair_t), INTENT(IN) :: cntp
    LOGICAL, INTENT(IN) :: dokinetic
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: sblk(:), tblk(:)

    INTEGER :: ig

!dir$ assume_aligned sblk : 64
!dir$ assume_aligned tblk : 64

    DO ig = 1, cntp%numpairs
        CALL comp_kin_ovl_int1_prim(cntp, ig, dokinetic, sblk, tblk)
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of Coulomb 1e integrals
!> @param[in]       cntp        shell pair data
!> @param[in]       xyz         coordinates of particles
!> @param[in]       c           charges of particles
!> @param[in]       nat         number of particles
!> @param[in]       chgtol      cut-off for charge
!> @param[inout]    blk         block of 1e Coulomb integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_coul_xyz_c(cntp, xyz, c, nat, chgtol, blk)
!dir$ attributes inline :: int1_coulxyz_c
    TYPE(shpair_t), INTENT(IN) :: cntp
    REAL(REAL64), CONTIGUOUS, INTENT(IN) :: xyz(:,:), c(:)
    REAL(REAL64), INTENT(IN) :: chgtol
    INTEGER, INTENT(IN) :: nat
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: blk(:)

    INTEGER :: ig, iat

!dir$ assume_aligned blk : 64

!   Interaction with point charge
    DO iat = 1, nat
        IF (abs(c(iat))<chgtol) CYCLE
        DO ig = 1, cntp%numpairs
            CALL comp_coulomb_int1_prim(cntp, ig, xyz(:,iat), -c(iat), blk)
        END DO
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of Coulomb 1e integrals
!> @param[in]       cntp        shell pair data
!> @param[in]       xyzc        coordinates and charges of particles
!> @param[in]       nat         number of particles
!> @param[in]       chgtol      cut-off for charge
!> @param[inout]    blk         block of 1e Coulomb integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 subroutine int1_coul_xyzc(cntp, xyzc, nat, chgtol, blk)
!dir$ attributes inline :: int1_coul_xyzc
    type(shpair_t), intent(in) :: cntp
    real(real64), contiguous, intent(in) :: xyzc(:)
    real(real64), intent(in) :: chgtol
    integer, intent(in) :: nat
    real(real64), contiguous, intent(inout) :: blk(:)

    integer :: ig, iat
    real(real64) :: c(3), znuc
!dir$ assume_aligned blk : 64

!   Interaction with point charge
    do iat = 1, nat
        c = xyzc((iat-1)*4+1:iat*4-1)
        znuc = xyzc(iat*4)
        if (abs(znuc)<chgtol) cycle
        do ig = 1, cntp%numpairs
            call comp_coulomb_int1_prim(cntp, ig, c(:), -znuc, blk)
        end do
    end do

 end subroutine

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of Coulomb 1e integrals
!> @param[in]       cntp        shell pair data
!> @param[in]       x           `X` coordinates of charged particles
!> @param[in]       y           `Y` coordinates of charged particles
!> @param[in]       z           `Z` coordinates of charged particles
!> @param[in]       c           charges of particles
!> @param[in]       nat         number of particles
!> @param[in]       chgtol      cut-off for charge
!> @param[inout]    blk         block of 1e Coulomb integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_coul_x_y_z_c(cntp, x, y, z, c, nat, chgtol, blk)
!dir$ attributes inline :: int1_coul_x_y_z_c
    TYPE(shpair_t), INTENT(IN) :: cntp
    REAL(REAL64), CONTIGUOUS, INTENT(IN) :: x(:), y(:), z(:), c(:)
    REAL(REAL64), INTENT(IN) :: chgtol
    INTEGER, INTENT(IN) :: nat
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: blk(:)

    INTEGER :: ig, iat
    REAL(REAL64) :: crd(3)
!dir$ assume_aligned blk : 64

!   Interaction with point charge
    DO iat = 1, nat
        IF (abs(c(iat))<chgtol) CYCLE
        crd(1) = x(iat)
        crd(2) = y(iat)
        crd(3) = z(iat)
        DO ig = 1, cntp%numpairs
            CALL comp_coulomb_int1_prim(cntp, ig, crd, -c(iat), blk)
        END DO
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Substract long-range Ewald contribution from the block of regular
!>  Coulomb 1e integrals
!> @param[in]       cntp        shell pair data
!> @param[in]       x           `X` coordinates of charged particles
!> @param[in]       y           `Y` coordinates of charged particles
!> @param[in]       z           `Z` coordinates of charged particles
!> @param[in]       c           charges of particles
!> @param[in]       nat         number of particles
!> @param[in]       chgtol      cut-off for charge
!> @param[in]       omega       Ewald splitting parameter
!> @param[inout]    blk         block of 1e Coulomb integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_ewald(cntp, x, y, z, c, nat, chgtol, omega, blk)
!dir$ attributes inline :: int1_ewald
    TYPE(shpair_t), INTENT(IN) :: cntp
    REAL(REAL64), CONTIGUOUS, INTENT(IN) :: x(:), y(:), z(:), c(:)
    REAL(REAL64), INTENT(IN) :: chgtol
    REAL(REAL64), INTENT(IN) :: omega
    INTEGER, INTENT(IN) :: nat
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: blk(:)

    INTEGER :: ig, iat
    REAL(REAL64) :: crd(3)
!dir$ assume_aligned blk : 64

!   Interaction with point charge
    DO iat = 1, nat
        IF (abs(c(iat))<chgtol) CYCLE
        crd(1) = x(iat)
        crd(2) = y(iat)
        crd(3) = z(iat)
        DO ig = 1, cntp%numpairs
!           By passing the charge as is, we effectively
!           subtract the long-range Ewald term
            CALL comp_ewaldlr_int1_prim(cntp, ig, crd, c(iat), omega, blk)
        END DO
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute sum of Coulomb integrals over pair of contracted shells
!> @param[in]       cntp        shell pair data
!> @param[in]       x           `X` coordinate of the charged particle
!> @param[in]       y           `Y` coordinate of the charged particle
!> @param[in]       z           `Z` coordinate of the charged particle
!> @param[in]       den         normalized density matrix block
!> @return          sum of Coulomb integrals over shell pair
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Oct, 2018_ Initial release
!
 function int1_epoten(cntp, x, y, z, den) result(vsum)
        !dir$ attributes inline :: int1_epoten
    type(shpair_t), intent(in) :: cntp
    real(real64), intent(in) :: x, y, z, den(:)
    real(real64) :: vsum

    integer :: ig
    real(real64) :: crd(3)

    vsum = 0.0
    crd(1) = x
    crd(2) = y
    crd(3) = z

    do ig = 1, cntp%numpairs
!       Interaction with unit charge
        call comp_coulpot_prim(cntp, ig, crd, den, vsum)
    end do

 end function

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of Z-angular momentum integrals
!> @param[in]       cntp        shell pair data
!> @param[inout]    blk         block of 1e Coulomb Lz-integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_lz(cntp, blk)
!dir$ attributes inline :: int1_lz
    TYPE(shpair_t), INTENT(IN) :: cntp
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: blk(:)

    INTEGER :: ig
!dir$ assume_aligned blk : 64

    DO ig = 1, cntp%numpairs
        CALL comp_lz_int1_prim(cntp, ig, blk)
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of multipole moment 1e integrals
!> @param[in]       cntp        shell pair data
!> @param[in]       r           point in space to compute integrals
!> @param[in]       mom         multiplole moment order (1-dipole, 2-quadrupole, 3-octopole)
!> @param[out]      blk         block of overlap integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_mul(cntp, r, mom, blk)
!dir$ attributes inline :: int1_kin_ovl
    type(shpair_t), intent(in) :: cntp
    real(real64), contiguous, intent(in) :: r(:)
    integer, intent(in) :: mom
    real(real64), contiguous, intent(inout) :: blk(:,:)

    integer :: ig

!dir$ assume_aligned blk : 64

    do ig = 1, cntp%numpairs
        call comp_mult_int1_prim(cntp, ig, r, mom, blk)
    end do

 end subroutine

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of multipole moment 1e integrals
!> @param[in]       cntp        shell pair data
!> @param[in]       r           point in space to compute integrals
!> @param[in]       mom         multiplole moment order (1-dipole, 2-quadrupole, 3-octopole)
!> @param[out]      blk         block of overlap integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_allmul(cntp, r, mxmom, blk)
!dir$ attributes inline :: int1_kin_ovl
    type(shpair_t), intent(in) :: cntp
    real(real64), contiguous, intent(in) :: r(:)
    integer, intent(in) :: mxmom
    real(real64), contiguous, intent(inout) :: blk(:,:)

    integer :: ig

!dir$ assume_aligned blk : 64

    do ig = 1, cntp%numpairs
        call comp_allmult_int1_prim(cntp, ig, r, mxmom, blk)
    end do

 end subroutine

end module
