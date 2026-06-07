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
        comp_amom_int1_prim, &
        comp_giao_overlap_deriv_prim, &
        comp_giao_h10_core_prim, &
        comp_nmr_dia_int1_prim, &
        comp_pso_int1_prim, &
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
    public omp_qmmm
    public multipole_integrals
    public angular_momentum_integrals
    public giao_overlap_derivative
    public giao_h10_core
    public nmr_dia_shielding
    public giao_a11part_corr
    public giao_a01gp_contract
    public pso_integrals
    public electrostatic_potential
    public electrostatic_potential_unweighted
    public external_charge_potential
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

!> @brief Driver for conventional ESP QM/MM integrals on a grid around QM atoms
!
!> @details  Compute one electron integrals and core Hamiltonian,
!>  - V is evaluated by Gauss-Rys quadrature, then \f$ h = T+V \f$
!>  Also, do \f$ L_z \f$ integrals for atoms and linear cases.
!>  Also, do FMO ESP integrals if needed.
!>  This subroutine is capable to do integrals in parallel using
!>  both OpenMP and MPI. It it helpful when running large FMO jobs.
!
!> @note Based on `HSANDT` subroutine from file `INT1.SRC`
!
!> @author Miquel Huix-Rotllant
!
!   REVISION HISTORY:
!> @date _Jul, 2024_ Initial release
!>
!> @param[in]       i       QM center index
!> @param[in]       ttt     ESP integral weight
!> @param[in,out]   chg_op  one-electron ESP atomic charge operator in packet format
 subroutine omp_qmmm(basis, i, coord, ttt, chg_op, nat, logtol)

    use precision, only: dp
    use basis_tools, only: basis_set

    type(basis_set), intent(in) :: basis
    real(real64), contiguous, intent(in) :: coord(:,:), ttt(:,:)
    real(real64), contiguous, intent(inout) :: chg_op(:)
    real(real64), optional, intent(in) :: logtol
    integer, intent(in) :: i, nat
    real(real64) :: tol

    integer :: l1, l2

    if (present(logtol)) then
        tol = logtol
    else
        tol = log(10.0_dp)*20
    end if

!   Exclude 1e potential in ESDIM, because density is used,
!   not point charges

    l1 = basis%nbf
    l2 = l1*(l1+1)/2

    chg_op(:)=0
    call nuc_ints(basis, coord, ttt(i,:), chg_op(:), tol)
    call bas_norm_matrix(chg_op(:), basis%bfnrm, l1)

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

!> @brief Compute the three angular-momentum 1e integral matrices about a
!>        gauge origin `o`, in packed (lower-triangular) storage.
!> @details The orbital angular momentum operator is anti-Hermitian, so in a
!>  real AO basis the matrices A_x, A_y, A_z returned here are antisymmetric
!>  (A_qp = -A_pq, zero diagonal). Only the unique lower triangle is stored;
!>  the caller is responsible for applying the antisymmetry when expanding to a
!>  full square matrix. The physical angular momentum is L = -i * A.
!
!> @param[in]       basis   basis set (without sp-shells)
!> @param[in,out]   ints    packed integrals, dimension (nbf2, 3) for x,y,z
!> @param[in]       o       gauge origin
!> @param[in]       debug   optional flag for debug printout
!> @param[in]       logtol  optional screening tolerance
 subroutine angular_momentum_integrals(basis, ints, o, debug, logtol)

    use io_constants, only: iw
    use precision, only: dp
    use basis_tools, only: basis_set
    use printing, only: print_sym_labeled

    type(basis_set), intent(in) :: basis
    real(real64), contiguous, intent(inout) :: ints(:,:)
    real(real64), intent(in) :: o(:)
    real(real64), optional, intent(in) :: logtol
    logical, optional, intent(in) :: debug

    character(len=*), parameter :: labels(3) = ['Lx', 'Ly', 'Lz']
    real(real64) :: tol
    logical :: dbug
    integer :: nbf, i

    if (ubound(ints,2) < 3) then
      call show_message('Insufficient space for angular momentum integrals', with_abort)
    end if

    dbug = .false.
    if (present(debug)) dbug = debug

    tol = log(10.0_dp)*20
    if (present(logtol)) tol = logtol
    nbf = basis%nbf

    ints = 0.0

    call amom_ints(ints, o, basis, tol)

!   Normalize 1-e integrals
    do i = 1, 3
      call bas_norm_matrix(ints(:,i), basis%bfnrm, nbf)
    end do

    if (dbug) then
       do i = 1, 3
         write(iw,*) 'Angular momentum integrals ('//trim(labels(i))//'), lower triangle'
         call print_sym_labeled(ints(:,i),nbf,basis)
       end do
    end if

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute the GIAO/London AO overlap magnetic derivative S10.
!> @details Returns the real coefficient of the imaginary first magnetic-field
!>  derivative of the overlap matrix for the three Cartesian magnetic-field
!>  components.  This is a native one-electron GIAO building block and remains
!>  disconnected from production NMR shielding until h10, two-electron derivative
!>  contractions, and GIAO CPHF/CPKS terms are implemented and benchmarked.
 subroutine giao_overlap_derivative(basis, ints, debug, logtol)

    use io_constants, only: iw
    use precision, only: dp
    use basis_tools, only: basis_set
    use printing, only: print_sym_labeled

    type(basis_set), intent(in) :: basis
    real(real64), contiguous, intent(inout) :: ints(:,:)
    real(real64), optional, intent(in) :: logtol
    logical, optional, intent(in) :: debug

    character(len=*), parameter :: labels(3) = ['Sx', 'Sy', 'Sz']
    real(real64) :: tol
    logical :: dbug
    integer :: nbf, i

    if (ubound(ints,2) < 3) then
      call show_message('Insufficient space for GIAO overlap derivative integrals', with_abort)
    end if

    dbug = .false.
    if (present(debug)) dbug = debug

    tol = log(10.0_dp)*20
    if (present(logtol)) tol = logtol
    nbf = basis%nbf

    ints = 0.0d0
    call giao_overlap_deriv_ints(ints, basis, tol)

    do i = 1, 3
      call bas_norm_matrix(ints(:,i), basis%bfnrm, nbf)
    end do

    if (dbug) then
       do i = 1, 3
         write(iw,*) 'GIAO overlap derivative integrals ('//trim(labels(i))//'), lower triangle'
         call print_sym_labeled(ints(:,i),nbf,basis)
       end do
    end if

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute the one-electron part of the RHF GIAO h10 magnetic derivative.
!> @details Returns packed lower-triangular real coefficients of the imaginary
!>  first-order GIAO core-Hamiltonian derivative for x/y/z magnetic-field
!>  components.  This routine intentionally contains only h10 one-electron
!>  kinetic+nuclear-attraction terms; it does not include the GIAO two-electron
!>  Fock derivative, CPHF/CPKS response, shielding assembly, or GIAO ungating.
 subroutine giao_h10_core(basis, coord, zq, ints, debug, logtol)

    use io_constants, only: iw
    use precision, only: dp
    use basis_tools, only: basis_set
    use printing, only: print_sym_labeled

    type(basis_set), intent(in) :: basis
    real(real64), contiguous, intent(in) :: coord(:,:), zq(:)
    real(real64), contiguous, intent(inout) :: ints(:,:)
    real(real64), optional, intent(in) :: logtol
    logical, optional, intent(in) :: debug

    character(len=*), parameter :: labels(3) = ['Hx', 'Hy', 'Hz']
    real(real64) :: tol
    logical :: dbug
    integer :: nbf, i

    if (ubound(ints,2) < 3) then
      call show_message('Insufficient space for GIAO h10 core derivative integrals', with_abort)
    end if

    dbug = .false.
    if (present(debug)) dbug = debug

    tol = log(10.0_dp)*20
    if (present(logtol)) tol = logtol
    nbf = basis%nbf

    ints = 0.0d0
    call giao_h10_core_ints(ints, basis, coord, zq, size(zq), tol)

    do i = 1, 3
      call bas_norm_matrix(ints(:,i), basis%bfnrm, nbf)
    end do

    if (dbug) then
       do i = 1, 3
         write(iw,*) 'GIAO h10 core derivative integrals ('//trim(labels(i))//'), lower triangle'
         call print_sym_labeled(ints(:,i),nbf,basis)
       end do
    end if

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Density-contracted NMR diamagnetic shielding integrals, all nuclei.
!> @details Returns g_ab(N) = sum_{mu,nu} D_{mu,nu} <mu|(r-o)_a (r-c_N)_b/|r-c_N|^3|nu>
!>  for every nucleus N. The caller assembles the diamagnetic shielding tensor as
!>    sigma^dia_{ts}(N) = (alpha^2/2) [ delta_ts (g_xx+g_yy+g_zz) - g_{s,t} ].
!> @param[in]   basis    basis set
!> @param[in]   denab    total density matrix, packed (lower triangle)
!> @param[in]   o        gauge origin
!> @param[in]   coords   nuclear coordinates (3, nat)
!> @param[in]   nat      number of nuclei
!> @param[out]  gdia     contracted integrals (3, 3, nat)
!> @param[in]   logtol   optional screening tolerance
 subroutine nmr_dia_shielding(basis, denab, o, coords, nat, gdia, logtol)
    use precision, only: dp
    use mathlib, only: unpack_matrix

    type(basis_set), intent(in) :: basis
    real(real64), contiguous, intent(in) :: denab(:)
    real(real64), intent(in) :: o(:)
    real(real64), contiguous, intent(in) :: coords(:,:)
    integer, intent(in) :: nat
    real(real64), intent(out) :: gdia(3,3,nat)
    real(real64), optional, intent(in) :: logtol

    real(real64), allocatable :: dens(:,:)
    real(real64) :: tol
    integer :: ii, jj, ic, nbf
    type(shell_t) :: shi, shj
    type(shpair_t) :: cntp

    tol = log(10.0_dp)*20
    if (present(logtol)) tol = logtol
    nbf = basis%nbf

    allocate(dens(nbf,nbf), source=0.0d0)
    call unpack_matrix(denab, dens, nbf, 'U')
    call bas_norm_matrix(dens, basis%bfnrm, nbf)

    gdia = 0.0d0

    call cntp%alloc(basis)

    do ii = 1, basis%nshell
      call shi%fetch_by_id(basis, ii)
      do jj = 1, basis%nshell
        call shj%fetch_by_id(basis, jj)
        call cntp%shell_pair(basis, shi, shj, tol)
        if (cntp%numpairs==0) cycle
        do ic = 1, nat
          call comp_nmr_dia_int1_prim(cntp, coords(:,ic), o, &
                 dens(basis%ao_offset(ii):, basis%ao_offset(jj):), gdia(:,:,ic))
        end do
      end do
    end do

    deallocate(dens)

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute PSO (paramagnetic spin-orbit) integral matrices for one nucleus.
!> @details Returns the three matrices A_a = [(r-c) x grad]_a/|r-c|^3 as FULL
!>  (nbf x nbf) antisymmetric matrices; the physical PSO operator is -i*A.
!>  The raw field+ket-derivative product acquires a small spurious symmetric
!>  component for nuclei not centered on a basis function. Since the exact PSO
!>  operator is anti-Hermitian (its real representation is antisymmetric with a
!>  zero diagonal), the full block is assembled and the symmetric part is removed
!>  via A = (M - M^T)/2, which is exact and discards only the spurious error.
!> @param[in]   basis    basis set
!> @param[in]   c        nucleus coordinates
!> @param[inout] ints    full integrals, dimension (nbf, nbf, 3), antisymmetric
!> @param[in]   logtol   optional screening tolerance
!> @brief GIAO a11part London correction, density-contracted.
!> @details The GIAO diamagnetic a11part integral satisfies (verified vs libcint)
!>   <mu|giao_a11part_{a,b}|nu> = <mu|cg_a11part(O=0)_{a,b}|nu>
!>                              + 0.5 * <mu|(r-R_N)_a/|r-R_N|^3|nu> * R_nu,b
!>  where R_nu is the KET shell center.  This routine returns the density-
!>  contracted correction tensor
!>   corr_{a,b}(N) = 0.5 * sum_{mu,nu} <mu|(r-R_N)_a/|r-R_N|^3|nu> * R_nu,b * D_{mu,nu}
!>  using the validated Hellmann-Feynman field integral (comp_coulomb_helfeyder1)
!>  on the density whose columns are pre-scaled by the ket-shell-center component.
!>  The caller adds this (trace-corrected) to the CGO diamagnetic at gauge origin
!>  0 (nmr_dia_shielding with o=0) to form the full GIAO a11part contribution.
 subroutine giao_a11part_corr(basis, denab, coords, nat, corr, logtol)
    use precision, only: dp
    use mathlib, only: unpack_matrix
    use mod_1e_primitives, only: comp_coulomb_helfeyder1

    type(basis_set), intent(in) :: basis
    real(real64), contiguous, intent(in) :: denab(:)
    real(real64), contiguous, intent(in) :: coords(:,:)
    integer, intent(in) :: nat
    real(real64), intent(out) :: corr(3,3,nat)
    real(real64), optional, intent(in) :: logtol

    real(real64), allocatable :: dens(:,:), densb(:,:), aoc(:,:)
    real(real64) :: tol, der(3)
    integer :: ii, jj, ic, nbf, b, ish, ao, k
    type(shell_t) :: shi, shj
    type(shpair_t) :: cntp

    tol = log(10.0_dp)*20
    if (present(logtol)) tol = logtol
    nbf = basis%nbf

    allocate(dens(nbf,nbf), densb(nbf,nbf), aoc(nbf,3), source=0.0d0)
    call unpack_matrix(denab, dens, nbf, 'U')
    call bas_norm_matrix(dens, basis%bfnrm, nbf)

    ! AO -> shell center map
    do ish = 1, basis%nshell
      do k = 1, basis%naos(ish)
        ao = basis%ao_offset(ish) + k - 1
        aoc(ao,1:3) = basis%shell_centers(ish,1:3)
      end do
    end do

    corr = 0.0d0
    call cntp%alloc(basis)

    do b = 1, 3
      ! scale ket (column) by its shell-center b-component
      do ao = 1, nbf
        densb(:,ao) = dens(:,ao)*aoc(ao,b)
      end do
      do ii = 1, basis%nshell
        call shi%fetch_by_id(basis, ii)
        do jj = 1, basis%nshell
          call shj%fetch_by_id(basis, jj)
          call cntp%shell_pair(basis, shi, shj, tol)
          if (cntp%numpairs==0) cycle
          do ic = 1, nat
            der = 0.0d0
            call comp_coulomb_helfeyder1(cntp, coords(:,ic), 1.0d0, &
                   densb(basis%ao_offset(ii):, basis%ao_offset(jj):), der)
            corr(1:3,b,ic) = corr(1:3,b,ic) + 0.5d0*der(1:3)
          end do
        end do
      end do
    end do

    deallocate(dens, densb, aoc)

 end subroutine

!> @brief GIAO a01gp gauge-correction, density-contracted (9 comp -> 3x3).
!> @details Returns e2_{a,col}(N) = sum_{mu,nu} <mu|a01gp_{a,col}|nu> D_{mu,nu}
!>  for each nucleus N, with a01gp the GIAO derivative of the PSO operator
!>  (comp_giao_a01gp_prim).  cvec = R_bra - R_ket per shell pair.  The caller
!>  adds this (NOT trace-corrected, per the standard diamagnetic decomposition) to the trace-corrected
!>  a11part contribution.
 subroutine giao_a01gp_contract(basis, denab, coords, nat, e2, logtol)
    use precision, only: dp
    use mathlib, only: unpack_matrix
    use mod_1e_primitives, only: comp_giao_a01gp_prim

    type(basis_set), intent(in) :: basis
    real(real64), contiguous, intent(in) :: denab(:)
    real(real64), contiguous, intent(in) :: coords(:,:)
    integer, intent(in) :: nat
    real(real64), intent(out) :: e2(3,3,nat)
    real(real64), optional, intent(in) :: logtol

    real(real64), allocatable :: dens(:,:)
    real(real64) :: tol, cvec(3), blk(blocksize,9)
    integer :: ii, jj, ic, nbf, a, col, i, j, ij, oi, oj
    type(shell_t) :: shi, shj
    type(shpair_t) :: cntp

    tol = log(10.0_dp)*20
    if (present(logtol)) tol = logtol
    nbf = basis%nbf

    allocate(dens(nbf,nbf), source=0.0d0)
    call unpack_matrix(denab, dens, nbf, 'U')
    call bas_norm_matrix(dens, basis%bfnrm, nbf)

    e2 = 0.0d0
    call cntp%alloc(basis)

    do ii = 1, basis%nshell
      call shi%fetch_by_id(basis, ii)
      oi = basis%ao_offset(ii)
      do jj = 1, basis%nshell
        call shj%fetch_by_id(basis, jj)
        oj = basis%ao_offset(jj)
        call cntp%shell_pair(basis, shi, shj, tol)
        if (cntp%numpairs==0) cycle
        cvec = basis%shell_centers(ii,1:3) - basis%shell_centers(jj,1:3)
        do ic = 1, nat
          blk = 0.0d0
          call comp_giao_a01gp_prim(cntp, coords(:,ic), cvec, blk)
          ! contract: e2(a,col) += sum_ij den(bra,ket) * blk(ij, (a-1)*3+col)
          ij = 0
          do i = 1, cntp%inao
            do j = 1, cntp%jnao
              ij = ij + 1
              do a = 1, 3
                do col = 1, 3
                  e2(a,col,ic) = e2(a,col,ic) &
                    + dens(oi+i-1, oj+j-1) * blk(ij,(a-1)*3+col)
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    deallocate(dens)

 end subroutine

 subroutine pso_integrals(basis, c, ints, logtol)
    use precision, only: dp
    type(basis_set), intent(in) :: basis
    real(real64), intent(in) :: c(:)
    real(real64), contiguous, intent(inout) :: ints(:,:,:)
    real(real64), optional, intent(in) :: logtol

    real(real64) :: tol
    integer :: ii, jj, m, nbf, p, q
    type(shell_t) :: shi, shj
    type(shpair_t) :: cntp
    real(real64), dimension(blocksize,3) :: blk
    real(real64) :: aij

    tol = log(10.0_dp)*20
    if (present(logtol)) tol = logtol
    nbf = basis%nbf

    ints = 0.0d0
    call cntp%alloc(basis)

!   Assemble the full (both-triangle) matrix M_a[bra,ket] for all shell pairs.
    do ii = 1, basis%nshell
      call shi%fetch_by_id(basis, ii)
      do jj = 1, basis%nshell
        call shj%fetch_by_id(basis, jj)
        call cntp%shell_pair(basis, shi, shj, tol)
        if (cntp%numpairs==0) cycle
        blk = 0.0d0
        call comp_pso_int1_prim(cntp, c, blk)
        do m = 1, 3
          ! blk is ordered (bra=shi outer, ket=shj inner); update_rectangular_matrix
          ! then writes ints(ket_global, bra_global) = <bra|A|ket>.
          call update_rectangular_matrix(shi, shj, blk(:,m), ints(:,:,m))
        end do
      end do
    end do

!   Normalize, then antisymmetrize A = (M - M^T)/2 (exact for the PSO operator).
    do m = 1, 3
      call bas_norm_matrix(ints(:,:,m), basis%bfnrm, nbf)
    end do
    ! update_rectangular_matrix stored ints(a,b) = <b|A|a>; antisymmetrise into
    ! the <bra|A|ket> convention used by the paramagnetic assembly.
    do m = 1, 3
      do p = 1, nbf
        do q = 1, p
          aij = 0.5d0*(ints(q,p,m) - ints(p,q,m))
          ints(p,q,m) =  aij
          ints(q,p,m) = -aij
        end do
      end do
    end do

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

!> @brief Compute unweighted electronic electrostatic potential on arbitrary points.
!
!> @details This is the safe public wrapper around the internal `int1_el_pot`
!> kernel. Unlike `electrostatic_potential`, this routine does not multiply by
!> quadrature weights. ddX expects `phi_cav` to be the unweighted electric
!> potential at cavity points, so this is the intended OpenQP entry point for
!> building ddX primal RHS data from an AO density.
!>
!> @param[inout] basis basis with SP-shells separated
!> @param[in]    x     x coordinates of evaluation points, in Bohr
!> @param[in]    y     y coordinates of evaluation points, in Bohr
!> @param[in]    z     z coordinates of evaluation points, in Bohr
!> @param[inout] d     packed AO density matrix; restored to input normalization
!> @param[out]   pot   unweighted electronic potential on points
!> @param[in]    logtol optional 1-e exponential prefactor tolerance
!>
 subroutine electrostatic_potential_unweighted(basis, x, y, z, d, pot, logtol)

    use precision, only: dp
    implicit none
    type(basis_set), intent(in)             :: basis
    real(real64), contiguous, intent(in)    :: x(:), y(:), z(:)
    real(real64), contiguous, intent(inout) :: d(:)
    real(real64), contiguous, intent(out)   :: pot(:)
    real(real64), optional, intent(in)      :: logtol
    real(real64) :: tol
    real(real64), allocatable :: invnrm(:)

    call bas_norm_matrix(d, basis%bfnrm, basis%nbf)

    tol = log(10.0_dp)*20
    if (present(logtol)) tol = logtol

    pot = 0.0_real64
    call int1_el_pot(basis, x, y, z, d, pot, tol)

    ! Restore the input normalization of d. Use a local inverse of the basis
    ! norms rather than bas_denorm_matrix, which would transiently mutate
    ! basis%bfnrm and so force an intent(inout) basis on this otherwise
    ! read-only routine (it is called from the intent(in) SCF Fock build).
    invnrm = 1.0_real64 / basis%bfnrm
    call bas_norm_matrix(d, invnrm, basis%nbf)

 end subroutine electrostatic_potential_unweighted

!-------------------------------------------------------------------------------

!> @brief Compute packed one-electron Coulomb potential from external point charges.
!
!> @details This is the normalized public wrapper around the internal
!> `int1_coul_ext_chg` kernel. It is intended for environment/solvent reaction
!> fields such as ddX apparent charges: given point charges q_k at coordinates
!> r_k, return the packed AO matrix sum_k q_k <mu|1/|r-r_k||nu>.
!>
!> @param[in]     basis  basis with SP-shells separated
!> @param[out]    v      packed normalized AO potential matrix
!> @param[in]     x      x coordinates of point charges, in Bohr
!> @param[in]     y      y coordinates of point charges, in Bohr
!> @param[in]     z      z coordinates of point charges, in Bohr
!> @param[in]     chg    point charges
!> @param[in]     logtol optional 1-e exponential prefactor tolerance
!> @param[in]     chgtol optional charge screening threshold
!>
 subroutine external_charge_potential(basis, v, x, y, z, chg, logtol, chgtol)

    use precision, only: dp
    implicit none
    type(basis_set), intent(in)              :: basis
    real(real64), contiguous, intent(out)    :: v(:)
    real(real64), contiguous, intent(in)     :: x(:), y(:), z(:), chg(:)
    real(real64), optional, intent(in)       :: logtol, chgtol
    real(real64) :: tol, qtol

    tol = log(10.0_dp)*20
    if (present(logtol)) tol = logtol

    qtol = 1.0d-12
    if (present(chgtol)) qtol = chgtol

    v = 0.0_real64
    call int1_coul_ext_chg(v, basis, size(chg), x, y, z, chg, tol, qtol)
    call bas_norm_matrix(v, basis%bfnrm, basis%nbf)

 end subroutine external_charge_potential

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

!  I shell
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
!> @author   Miquel Huix-Rotllant
!
!     REVISION HISTORY:
!> @date _Jul, 2024_ Initial release
!
 SUBROUTINE int1_coul_xyz_u(cntp, xyz, c, blk)
!dir$ attributes inline :: int1_coulxyz_u
    TYPE(shpair_t), INTENT(IN) :: cntp
    REAL(REAL64), CONTIGUOUS, INTENT(IN) :: xyz(:)
    REAL(REAL64), INTENT(IN) :: c
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: blk(:)

    INTEGER :: ig

!dir$ assume_aligned blk : 64

!   Interaction with point charge
    DO ig = 1, cntp%numpairs
        CALL comp_coulomb_int1_prim(cntp, ig, xyz, -c, blk)
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

!> @brief Compute angular momentum integrals about gauge origin `o`
!> @author   Generated for NMR shielding (CGO)
!
 SUBROUTINE amom_ints(ints, o, basis, tol)

    REAL(REAL64), CONTIGUOUS,  INTENT(INOUT)  :: ints(:,:)
    TYPE(basis_set), INTENT(IN)     :: basis
    real(real64), contiguous, intent(in) :: o(:)
    REAL(REAL64),   INTENT(IN)     :: tol

    INTEGER :: ii, jj, m

    REAL(REAL64), DIMENSION(BLOCKSIZE,3) :: blk
!dir$ attributes align : 64 :: blk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, m, &
!$omp       blk, &
!$omp       shi, shj, cntp &
!$omp   )

    CALL cntp%alloc(basis)

    DO ii = 1, basis%nshell

        CALL shi%fetch_by_id(basis,ii)

!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(basis,jj)

            CALL cntp%shell_pair(basis,shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            blk = 0.0

            CALL int1_amom(cntp, o, blk)

            do m = 1, 3
              CALL update_triang_matrix(shi, shj, blk(:,m), ints(:,m))
            end do

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of angular momentum 1e integrals
!> @param[in]       cntp        shell pair data
!> @param[in]       o           gauge origin
!> @param[inout]    blk         block of 1e angular momentum integrals (:,1:3)
!> @author   Generated for NMR shielding (CGO)
!
 SUBROUTINE int1_amom(cntp, o, blk)
!dir$ attributes inline :: int1_amom
    type(shpair_t), intent(in) :: cntp
    real(real64), contiguous, intent(in) :: o(:)
    real(real64), contiguous, intent(inout) :: blk(:,:)

    integer :: ig
!dir$ assume_aligned blk : 64

    do ig = 1, cntp%numpairs
        call comp_amom_int1_prim(cntp, ig, o, blk)
    end do

 end subroutine

!--------------------------------------------------------------------------------

!> @brief Compute GIAO/London overlap magnetic derivative integrals.
 SUBROUTINE giao_overlap_deriv_ints(ints, basis, tol)

    REAL(REAL64), CONTIGUOUS,  INTENT(INOUT)  :: ints(:,:)
    TYPE(basis_set), INTENT(IN)     :: basis
    REAL(REAL64),   INTENT(IN)     :: tol

    INTEGER :: ii, jj, m
    REAL(REAL64), DIMENSION(BLOCKSIZE,3) :: blk
!dir$ attributes align : 64 :: blk
    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, m, &
!$omp       blk, &
!$omp       shi, shj, cntp &
!$omp   )

    CALL cntp%alloc(basis)

    DO ii = 1, basis%nshell
        CALL shi%fetch_by_id(basis,ii)
!$omp do schedule(dynamic)
        DO jj = 1, ii
            CALL shj%fetch_by_id(basis,jj)
            CALL cntp%shell_pair(basis,shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE
            blk = 0.0d0
            CALL int1_giao_overlap_deriv(cntp, blk)
            do m = 1, 3
              CALL update_triang_matrix(shi, shj, blk(:,m), ints(:,m))
            end do
        END DO
!$omp end do nowait
    END DO
!$omp end parallel
 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of GIAO/London overlap derivative integrals.
 SUBROUTINE int1_giao_overlap_deriv(cntp, blk)
!dir$ attributes inline :: int1_giao_overlap_deriv
    type(shpair_t), intent(in) :: cntp
    real(real64), contiguous, intent(inout) :: blk(:,:)

    integer :: ig
!dir$ assume_aligned blk : 64

    do ig = 1, cntp%numpairs
        call comp_giao_overlap_deriv_prim(cntp, ig, blk)
    end do

 end subroutine

!--------------------------------------------------------------------------------

!> @brief Compute one-electron GIAO h10 core derivative integrals.
 SUBROUTINE giao_h10_core_ints(ints, basis, coord, zq, nat, tol)

    REAL(REAL64), CONTIGUOUS,  INTENT(INOUT)  :: ints(:,:)
    TYPE(basis_set), INTENT(IN)     :: basis
    real(real64), contiguous, intent(in) :: coord(:,:), zq(:)
    integer, intent(in) :: nat
    REAL(REAL64),   INTENT(IN)     :: tol

    INTEGER :: ii, jj, m
    REAL(REAL64), DIMENSION(BLOCKSIZE,3) :: blk
!dir$ attributes align : 64 :: blk
    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, m, &
!$omp       blk, &
!$omp       shi, shj, cntp &
!$omp   )

    CALL cntp%alloc(basis)

    DO ii = 1, basis%nshell
        CALL shi%fetch_by_id(basis,ii)
!$omp do schedule(dynamic)
        DO jj = 1, ii
            CALL shj%fetch_by_id(basis,jj)
            CALL cntp%shell_pair(basis,shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE
            blk = 0.0d0
            CALL int1_giao_h10_core(cntp, coord, zq, nat, blk)
            do m = 1, 3
              CALL update_triang_matrix(shi, shj, blk(:,m), ints(:,m))
            end do
        END DO
!$omp end do nowait
    END DO
!$omp end parallel
 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of one-electron GIAO h10 derivative integrals.
!> @details Matches the one-electron libcint convention
!>  h10_onee = -0.5*int1e_giao_irjxp - int1e_ignuc(asym) - int1e_igkin;
!>  the two-electron GIAO magnetic Fock derivative is intentionally absent here.
 SUBROUTINE int1_giao_h10_core(cntp, coord, zq, nat, blk)
!dir$ attributes inline :: int1_giao_h10_core
    type(shpair_t), intent(in) :: cntp
    real(real64), contiguous, intent(in) :: coord(:,:), zq(:)
    integer, intent(in) :: nat
    real(real64), contiguous, intent(inout) :: blk(:,:)

    integer :: ig
    real(real64), dimension(BLOCKSIZE,3) :: amom_blk
!dir$ assume_aligned blk : 64
!dir$ assume_aligned amom_blk : 64

    amom_blk = 0.0_real64
    do ig = 1, cntp%numpairs
        call comp_giao_h10_core_prim(cntp, ig, coord, zq, nat, blk)
        call comp_amom_int1_prim(cntp, ig, cntp%rj, amom_blk)
    end do
    ! libcint/the reference int1e_giao_irjxp is the negative transpose of the
    ! angular-momentum block returned by comp_amom_int1_prim for the current
    ! (bra=shi, ket=shj) shell-pair convention.  The packed lower-triangle
    ! h10 contribution is therefore -0.5*irjxp = +0.5*amom_blk.
    blk = blk + 0.5_real64*amom_blk

 end subroutine

!--------------------------------------------------------------------------------


!--------------------------------------------------------------------------------


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
