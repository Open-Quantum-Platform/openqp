module cphf_mod
!> @brief Native coupled-perturbed Hartree-Fock / Kohn-Sham (CPHF/CPKS) solver
!>   for closed-shell (RHF/RKS) references.
!>
!>   The static CPHF A-matrix is the orbital Hessian (A+B)_{ia,jb}, the same
!>   operator the TDDFT Z-vector solver applies. This module reuses that exact
!>   operator -- built from the native Rys 2e engine via int2_td_data_t plus the
!>   DFT XC kernel (tddft_fxc) -- so it has no libint dependency. It drives the
!>   existing pcg solver with:
!>     update  : U(MO,occ-vir) -> AO density (iatogen) -> response Fock (A+B)
!>               -> MO occ-vir (mntoia) + orbital-energy diagonal (e_a-e_i) U
!>     precond : diagonal 1/(e_a-e_i)
!>   to solve  A U = B  for an arbitrary occ-vir right-hand side B.
!>
!>   cphf_solve is the reusable entry point (used by the analytic Hessian for the
!>   nuclear-perturbation response). cphf_polarizability_selftest validates the
!>   solver end to end against a known property: it solves with the dipole
!>   right-hand side and forms the static dipole polarizability, written to a file
!>   for comparison with PySCF (no libint, no geometry derivatives required).

  use precision, only: dp
  use iso_c_binding, only: c_ptr, c_loc, c_f_pointer
  use types, only: information
  use basis_tools, only: basis_set
  use int2_compute, only: int2_compute_t, int2_fock_data_t
  use tdhf_lib, only: int2_td_data_t, iatogen, mntoia
  use mod_dft_molgrid, only: dft_grid_t
  use pcg_mod, only: pcg_t, PCG_OK, PCG_CONVERGED
  use io_constants, only: iw

  implicit none

  character(len=*), parameter :: module_name = "cphf_mod"

  !> Opaque data passed to the PCG callbacks (the A-matrix action).
  type :: cphf_cg_data
    type(information), pointer :: infos => null()
    type(int2_compute_t), pointer :: int2_driver => null()
    class(int2_fock_data_t), pointer :: int2_data => null()
    type(dft_grid_t), pointer :: molgrid => null()
    real(kind=dp), pointer :: wrk(:,:) => null()
    real(kind=dp), pointer :: mo(:,:) => null()
    real(kind=dp), pointer :: pa(:,:,:) => null()
    real(kind=dp), pointer :: xm(:) => null()      ! (e_a - e_i), length nocc*nvir
    real(kind=dp), pointer :: xminv(:) => null()   ! 1/(e_a - e_i)
    integer :: nbf = 0
    integer :: nocc = 0
    logical :: dft = .false.
  end type

  !> Opaque data for the open-shell (UHF) A-matrix action.  The rotation vector
  !> is the concatenation of the alpha occ-vir block (length la = nocca*nvira)
  !> and the beta occ-vir block (length lb = noccb*nvirb).
  type :: cphf_cg_data_uhf
    type(information), pointer :: infos => null()
    type(basis_set), pointer :: basis => null()
    type(dft_grid_t), pointer :: molgrid => null()
    real(kind=dp), pointer :: moa(:,:) => null()
    real(kind=dp), pointer :: mob(:,:) => null()
    real(kind=dp), pointer :: xm(:) => null()      ! (e_a - e_i) for [alpha; beta]
    real(kind=dp), pointer :: xminv(:) => null()   ! 1/(e_a - e_i)
    real(kind=dp), pointer :: wrka(:,:) => null()  ! nbf x nbf scratch (alpha)
    real(kind=dp), pointer :: wrkb(:,:) => null()  ! nbf x nbf scratch (beta)
    integer :: nbf = 0
    integer :: nocca = 0
    integer :: noccb = 0
    integer :: la = 0
    integer :: lb = 0
    real(kind=dp) :: scale_exch = 1.0_dp
    logical :: dft = .false.
  end type

  private
  public :: cphf_solve
  public :: cphf_solve_uhf
  public :: cphf_static_polarizability
  public :: cphf_static_polarizability_C
  public :: cphf_polarizability_selftest
  public :: cphf_polarizability_selftest_C
  public :: cphf_uhf_polarizability_selftest
  public :: cphf_uhf_polarizability_selftest_C

contains

!###############################################################################

!> @brief Solve A U = B for closed-shell CPHF, B and U in MO occ-vir layout
!>   (nocc*nvir, nrhs), matching the iatogen/mntoia convention.
!> @param[in]    infos   system/control information (must have a converged RHF/RKS)
!> @param[in]    nrhs    number of right-hand sides
!> @param[in]    bvec    (nocc*nvir, nrhs) right-hand sides
!> @param[out]   uvec    (nocc*nvir, nrhs) solutions
!> @param[in]    tol     CG tolerance (optional)
!> @param[in]    maxit   max CG iterations (optional)
  subroutine cphf_solve(infos, nrhs, bvec, uvec, tol, maxit)
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_E_MO_A, OQP_VEC_MO_A
    use dft, only: dft_initialize
    real(kind=dp), parameter :: default_tol = 1.0d-9
    type(information), target, intent(inout) :: infos
    integer, intent(in) :: nrhs
    real(kind=dp), intent(in) :: bvec(:,:)
    real(kind=dp), intent(out) :: uvec(:,:)
    real(kind=dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxit

    type(basis_set), pointer :: basis
    type(dft_grid_t), target :: molgrid
    type(int2_compute_t), target :: int2_driver
    type(int2_td_data_t), target :: int2_data
    type(cphf_cg_data), target :: cgdata
    type(pcg_t) :: pcg

    real(kind=dp), contiguous, pointer :: mo_a(:,:), mo_energy_a(:)
    real(kind=dp), allocatable, target :: wrk1(:,:), pa(:,:,:), xm(:), xminv(:)
    real(kind=dp), pointer :: pxm(:,:)
    integer :: nbf, nocc, nvir, lexc, i, j, irhs, iter, mxit
    integer :: clock_rate, clock_start, clock_stop, rhs_clock_start, rhs_clock_stop
    logical :: dft
    real(kind=dp) :: cnv, scale_exch
    real(kind=dp) :: cpu_start, cpu_stop, rhs_cpu_start, rhs_cpu_stop, rhs_wall

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nocc = infos%mol_prop%nocc
    nvir = nbf - nocc
    lexc = nocc*nvir
    dft = infos%control%hamilton == 20
    cnv = default_tol; if (present(tol)) cnv = tol
    mxit = 100; if (present(maxit)) mxit = maxit
    if (mxit < lexc + 5) mxit = lexc + 5

    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

    if (dft) call dft_initialize(infos, basis, molGrid)

    allocate(wrk1(nbf,nbf), pa(nbf,nbf,1), xm(lexc), xminv(lexc), source=0.0_dp)

    ! orbital-energy difference diagonal (e_a - e_i), occ-vir layout
    pxm(1:nocc,1:nvir) => xm(1:)
    do i = 1, nvir
      do j = 1, nocc
        pxm(j,i) = mo_energy_a(nocc+i) - mo_energy_a(j)
      end do
    end do
    xminv = 1.0_dp/xm

    scale_exch = 1.0_dp
    if (dft) scale_exch = infos%dft%HFscale

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    int2_data = int2_td_data_t(d2=pa, &
            int_apb=.true., int_amb=.false., &
            tamm_dancoff=.false., scale_exchange=scale_exch)

    cgdata%infos => infos
    cgdata%int2_driver => int2_driver
    cgdata%int2_data => int2_data
    cgdata%molgrid => molgrid
    cgdata%wrk => wrk1
    cgdata%mo => mo_a
    cgdata%pa => pa
    cgdata%xm => xm
    cgdata%xminv => xminv
    cgdata%nbf = nbf
    cgdata%nocc = nocc
    cgdata%dft = dft

    call system_clock(count_rate=clock_rate)
    call system_clock(clock_start)
    call cpu_time(cpu_start)
    write(iw,'(/3x,60("-"))')
    write(iw,'(6x,"CPHF/CPKS iterative solver")')
    write(iw,'(6x,"right-hand sides =",I5,3x,"nocc =",I5,3x,"nvir =",I5)') &
            nrhs, nocc, nvir
    write(iw,'(6x,"tolerance =",1P,E10.3,3x,"max iterations =",I6)') cnv, mxit
    write(iw,'(3x,60("-"))')

    do irhs = 1, nrhs
      call system_clock(rhs_clock_start)
      call cpu_time(rhs_cpu_start)
      call pcg%init(b=bvec(:,irhs), update=cphf_apbx, precond=cphf_precond, &
                    dat=cgdata, tol=sqrt(abs(cnv)))
      write(iw,'(" INITIAL CPHF ERROR RHS",I5," =",3X,' // &
               '1P,E10.3,1X,"/",1P,E10.3)') &
              irhs, pcg%error**2, cnv
      do iter = 1, mxit
        if (pcg%errcode /= PCG_OK) exit
        call pcg%step()
        write(iw,'(" CPHF ITER RHS",I5," ITER#",I4," ERROR =",3X,' // &
                 '1P,E10.3,1X,"/",1P,E10.3)') &
                irhs, iter, pcg%error**2, cnv
        call flush(iw)
      end do
      call system_clock(rhs_clock_stop)
      call cpu_time(rhs_cpu_stop)
      rhs_wall = real(rhs_clock_stop - rhs_clock_start, kind=dp) / real(clock_rate, kind=dp)
      write(iw,'(" CPHF RHS",I5," completed in",I5," iterations;",' // &
               '" CPU time =",F10.3," s; wall time =",F10.3," s")') &
              irhs, iter - 1, rhs_cpu_stop - rhs_cpu_start, rhs_wall
      call flush(iw)
      uvec(:,irhs) = pcg%x
      call pcg%clean()
    end do

    call system_clock(clock_stop)
    call cpu_time(cpu_stop)
    write(iw,'(6x,"CPHF wall time =",F10.3," s; CPU time =",F10.3," s"/)') &
            real(clock_stop - clock_start, kind=dp) / real(clock_rate, kind=dp), cpu_stop - cpu_start
    call flush(iw)

    call int2_driver%clean()
    deallocate(wrk1, pa, xm, xminv)
  end subroutine cphf_solve

!###############################################################################

!> @brief A-matrix action y = (A+B) x, mirroring tdhf_z_vector::compute_apbx.
  subroutine cphf_apbx(y, x, dat)
    use mathlib, only: symmetrize_matrix, orthogonal_transform
    use mod_dft_gridint_fxc, only: tddft_fxc
    real(kind=dp) :: x(:)
    real(kind=dp) :: y(:)
    type(c_ptr) :: dat
    type(cphf_cg_data), pointer :: p
    real(kind=dp), pointer :: apb(:,:,:)

    call c_f_pointer(dat, p)
    associate( wrk => p%wrk, nocc => p%nocc, nbf => p%nbf, mo => p%mo, &
               pa => p%pa, int2_driver => p%int2_driver, int2_data => p%int2_data, &
               infos => p%infos, molgrid => p%molgrid, dft => p%dft, xm => p%xm )

      call iatogen(x, wrk, nocc, nocc)
      call symmetrize_matrix(wrk, nbf)
      call orthogonal_transform('t', nbf, mo, wrk, pa(:,:,1))

      call int2_driver%run(int2_data, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%dft%cam_alpha, beta=infos%dft%cam_beta, mu=infos%dft%cam_mu)

      select type (int2_data)
      type is (int2_td_data_t)
        apb => int2_data%apb(:,:,:,1)
      end select
      apb = apb*0.5_dp

      if (dft) then
        call tddft_fxc(basis=infos%basis, molGrid=molGrid, isVecs=.true., wf=mo, &
                       fx=apb(:,:,1:1), dx=pa(:,:,1:1), nmtx=1, threshold=0.0d0, infos=infos)
      end if

      call mntoia(apb(:,:,1), y, mo, mo, nocc, nocc)
      y = y + xm*x
    end associate
  end subroutine cphf_apbx

!###############################################################################

  subroutine cphf_precond(y, x, dat)
    real(kind=dp) :: x(:)
    real(kind=dp) :: y(:)
    type(c_ptr) :: dat
    type(cphf_cg_data), pointer :: p
    call c_f_pointer(dat, p)
    y = p%xminv*x
  end subroutine cphf_precond

!###############################################################################

  subroutine cphf_polarizability_selftest_C(c_handle) bind(C, name="cphf_polarizability_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call cphf_polarizability_selftest(inf)
  end subroutine cphf_polarizability_selftest_C

  subroutine cphf_static_polarizability_C(c_handle, alpha) bind(C, name="cphf_static_polarizability")
    use iso_c_binding, only: c_double
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    type(oqp_handle_t) :: c_handle
    real(c_double), intent(out) :: alpha(3,3)
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call cphf_static_polarizability(inf, alpha)
  end subroutine cphf_static_polarizability_C

!> @brief Compute native closed-shell static dipole polarizability.
!>   For each Cartesian q, the perturbation is the dipole operator; the MO
!>   occ-vir RHS is B^q_{ia} = -<i|q|a> (in MO basis). Solving A U^q = B^q gives
!>   the orbital response, and alpha_pq = -4 sum_{ia} mu^p_{ia} U^q_{ia}.
  subroutine cphf_static_polarizability(infos, alpha)
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A
    use int1, only: multipole_integrals
    use mathlib, only: unpack_matrix
    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(out) :: alpha(3,3)

    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: mo_a(:,:)
    real(kind=dp), allocatable :: mints(:,:), dipfull(:,:), dip_mo(:,:)
    real(kind=dp), allocatable :: bvec(:,:), uvec(:,:), scr(:,:)
    real(kind=dp) :: origin(3)
    integer :: nbf, nbf2, nocc, nvir, lexc, q, i, a, ia

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nocc = infos%mol_prop%nocc
    nvir = nbf - nocc
    lexc = nocc*nvir

    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

    ! dipole integrals about the origin (first 3 of the multipole set: X,Y,Z)
    allocate(mints(nbf2,19), source=0.0_dp)
    origin = 0.0_dp
    call multipole_integrals(basis, mints, origin, 3)

    allocate(dipfull(nbf,nbf), dip_mo(nbf,nbf), scr(nbf,nbf))
    allocate(bvec(lexc,3), uvec(lexc,3), source=0.0_dp)

    ! Build MO-basis dipole and the occ-vir RHS B^q_{ia} = -mu^q_{ia}
    do q = 1, 3
      call unpack_matrix(mints(:,q), dipfull)
      ! MO transform: dip_mo = C^T (dipfull) C
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mo_a, nbf, dipfull, nbf, 0.0_dp, scr, nbf)
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, scr, nbf, mo_a, nbf, 0.0_dp, dip_mo, nbf)
      ia = 0
      do a = 1, nvir
        do i = 1, nocc
          ia = ia + 1
          bvec(ia,q) = -dip_mo(i, nocc+a)
        end do
      end do
    end do

    call cphf_solve(infos, 3, bvec, uvec)

    ! alpha_pq = -4 sum_ia mu^p_ia U^q_ia  (closed shell)
    alpha = 0.0_dp
    do q = 1, 3
      do i = 1, 3
        ! recompute mu^p_ia from bvec (= -mu) : mu = -bvec
        alpha(i,q) = -4.0_dp * sum( (-bvec(:,i)) * uvec(:,q) )
      end do
    end do

    deallocate(mints, dipfull, dip_mo, scr, bvec, uvec)
  end subroutine cphf_static_polarizability

!> @brief Validate the CPHF solver via the reusable static dipole polarizability.
!>   Writes the 3x3 tensor to /tmp/cphf_polar.out for comparison with PySCF.
  subroutine cphf_polarizability_selftest(infos)
    type(information), target, intent(inout) :: infos

    real(kind=dp) :: alpha(3,3)
    integer :: i, u

    call cphf_static_polarizability(infos, alpha)

    open(newunit=u, file='/tmp/cphf_polar.out', status='replace', action='write')
    write(u,'(a)') 'CPHF static dipole polarizability (a.u.):'
    do i = 1, 3
      write(u,'(3f16.8)') alpha(i,1:3)
    end do
    write(u,'(a,f16.8)') 'isotropic = ', (alpha(1,1)+alpha(2,2)+alpha(3,3))/3.0_dp
    close(u)
  end subroutine cphf_polarizability_selftest

!###############################################################################
!  Open-shell (UHF) CPHF solver
!###############################################################################

!> @brief Solve the open-shell (UHF) CPHF equations  M U = B.
!>
!>   The unknown/RHS vectors are laid out as the concatenation of the alpha
!>   occ-vir block (length la = nocca*nvira) followed by the beta occ-vir block
!>   (length lb = noccb*nvirb), each in the iatogen/mntoia (occ-major) order.
!>
!>   The UHF orbital-Hessian action on a trial rotation U is
!>       (M U)^sigma_ia = (e^sigma_a - e^sigma_i) U^sigma_ia
!>                        + [ C^sigma^T  dF^sigma  C^sigma ]_ia ,
!>       dF^sigma = J[dP^alpha + dP^beta] - c_x K[dP^sigma]  (+ f_xc for KS),
!>       dP^sigma_mn = sum_ia ( C^s_mi U^s_ia C^s_na + C^s_ma U^s_ia C^s_ni ).
!>   The Coulomb response is built from the spin-summed trial density and the
!>   exchange response from the same-spin trial density, exactly the open-shell
!>   two-electron Fock that scf_addons::fock_jk assembles for scftype>=2.
!>
!>   This is the genuine static CPHF operator (not the TDDFT A+B), so it serves
!>   the open-shell analytic Hessian nuclear-perturbation response and the
!>   open-shell static dipole polarizability on the same footing.
  subroutine cphf_solve_uhf(infos, nrhs, bvec, uvec, tol, maxit)
    use oqp_tagarray_driver, only: tagarray_get_data, &
        OQP_E_MO_A, OQP_VEC_MO_A, OQP_E_MO_B, OQP_VEC_MO_B
    use dft, only: dft_initialize
    real(kind=dp), parameter :: default_tol = 1.0d-9
    type(information), target, intent(inout) :: infos
    integer, intent(in) :: nrhs
    real(kind=dp), intent(in) :: bvec(:,:)
    real(kind=dp), intent(out) :: uvec(:,:)
    real(kind=dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxit

    type(basis_set), pointer :: basis
    type(dft_grid_t), target :: molgrid
    type(cphf_cg_data_uhf), target :: cgdata
    type(pcg_t) :: pcg

    real(kind=dp), contiguous, pointer :: moa(:,:), mob(:,:), epsa(:), epsb(:)
    real(kind=dp), allocatable, target :: wrka(:,:), wrkb(:,:), xm(:), xminv(:)
    real(kind=dp), pointer :: pxm(:,:)
    integer :: nbf, nocca, noccb, nvira, nvirb, la, lb, ltot
    integer :: i, j, irhs, iter, mxit, off
    logical :: dft
    real(kind=dp) :: cnv, scale_exch

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    nvirb = nbf - noccb
    la = nocca*nvira
    lb = noccb*nvirb
    ltot = la + lb
    dft = infos%control%hamilton == 20
    cnv = default_tol; if (present(tol)) cnv = tol
    mxit = 100; if (present(maxit)) mxit = maxit
    if (mxit < ltot + 5) mxit = ltot + 5

    call tagarray_get_data(infos%dat, OQP_E_MO_A, epsa)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, moa)
    call tagarray_get_data(infos%dat, OQP_E_MO_B, epsb)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mob)

    if (dft) call dft_initialize(infos, basis, molGrid)

    allocate(wrka(nbf,nbf), wrkb(nbf,nbf), xm(ltot), xminv(ltot), source=0.0_dp)

    ! orbital-energy difference diagonal (e_a - e_i), occ-vir (occ-major) layout
    if (la > 0) then
      pxm(1:nocca,1:nvira) => xm(1:la)
      do i = 1, nvira
        do j = 1, nocca
          pxm(j,i) = epsa(nocca+i) - epsa(j)
        end do
      end do
    end if
    if (lb > 0) then
      pxm(1:noccb,1:nvirb) => xm(la+1:ltot)
      do i = 1, nvirb
        do j = 1, noccb
          pxm(j,i) = epsb(noccb+i) - epsb(j)
        end do
      end do
    end if
    xminv = 1.0_dp/xm

    scale_exch = 1.0_dp
    if (dft) scale_exch = infos%dft%HFscale

    cgdata%infos => infos
    cgdata%basis => basis
    cgdata%molgrid => molgrid
    cgdata%moa => moa
    cgdata%mob => mob
    cgdata%xm => xm
    cgdata%xminv => xminv
    cgdata%wrka => wrka
    cgdata%wrkb => wrkb
    cgdata%nbf = nbf
    cgdata%nocca = nocca
    cgdata%noccb = noccb
    cgdata%la = la
    cgdata%lb = lb
    cgdata%scale_exch = scale_exch
    cgdata%dft = dft

    write(iw,'(/3x,60("-"))')
    write(iw,'(6x,"open-shell (UHF) CPHF iterative solver")')
    write(iw,'(6x,"right-hand sides =",I5,3x,"la =",I6,3x,"lb =",I6)') nrhs, la, lb
    write(iw,'(6x,"tolerance =",1P,E10.3,3x,"max iterations =",I6)') cnv, mxit
    write(iw,'(3x,60("-"))')

    off = 0
    do irhs = 1, nrhs
      call pcg%init(b=bvec(:,irhs), update=cphf_apbx_uhf, precond=cphf_precond_uhf, &
                    dat=cgdata, tol=sqrt(abs(cnv)))
      do iter = 1, mxit
        if (pcg%errcode /= PCG_OK) exit
        call pcg%step()
      end do
      write(iw,'(" UHF CPHF RHS",I5," completed in",I5," iterations; error =",1P,E10.3)') &
              irhs, iter - 1, pcg%error**2
      call flush(iw)
      uvec(:,irhs) = pcg%x
      call pcg%clean()
    end do

    deallocate(wrka, wrkb, xm, xminv)
  end subroutine cphf_solve_uhf

!###############################################################################

!> @brief Open-shell (UHF) A-matrix action  y = M x  (see cphf_solve_uhf).
  subroutine cphf_apbx_uhf(y, x, dat)
    use mathlib, only: symmetrize_matrix, orthogonal_transform, pack_matrix, unpack_matrix
    use mod_dft_gridint_fxc, only: utddft_fxc
    use scf_addons, only: fock_jk
    real(kind=dp) :: x(:)
    real(kind=dp) :: y(:)
    type(c_ptr) :: dat
    type(cphf_cg_data_uhf), pointer :: p

    real(kind=dp), allocatable :: pa_ao(:,:), pb_ao(:,:), dpack(:,:), fpack(:,:)
    real(kind=dp), allocatable :: ga(:,:,:), gb(:,:,:), dxa(:,:,:), dxb(:,:,:)
    integer :: nbf, nbf2, nocca, noccb, la, lb

    call c_f_pointer(dat, p)
    nbf = p%nbf; nbf2 = nbf*(nbf+1)/2
    nocca = p%nocca; noccb = p%noccb; la = p%la; lb = p%lb

    allocate(pa_ao(nbf,nbf), pb_ao(nbf,nbf), source=0.0_dp)
    allocate(ga(nbf,nbf,1), gb(nbf,nbf,1), source=0.0_dp)
    allocate(dpack(nbf2,2), fpack(nbf2,2), source=0.0_dp)

    ! Trial AO densities from the occ-vir rotation amplitudes (per spin).
    if (la > 0) then
      call iatogen(x(1:la), p%wrka, nocca, nocca)
      call symmetrize_matrix(p%wrka, nbf)
      call orthogonal_transform('t', nbf, p%moa, p%wrka, pa_ao)
    end if
    if (lb > 0) then
      call iatogen(x(la+1:la+lb), p%wrkb, noccb, noccb)
      call symmetrize_matrix(p%wrkb, nbf)
      call orthogonal_transform('t', nbf, p%mob, p%wrkb, pb_ao)
    end if

    call pack_matrix(pa_ao, dpack(:,1))
    call pack_matrix(pb_ao, dpack(:,2))

    ! Open-shell two-electron response Fock: dF^s = J[dPa+dPb] - cx K[dP^s].
    call fock_jk(p%basis, d=dpack, f=fpack, scale_exch=p%scale_exch, infos=p%infos)
    call unpack_matrix(fpack(:,1), ga(:,:,1))
    call unpack_matrix(fpack(:,2), gb(:,:,1))

    ! XC response kernel (UKS): spin-resolved f_xc on the trial spin densities.
    if (p%dft) then
      allocate(dxa(nbf,nbf,1), dxb(nbf,nbf,1))
      dxa(:,:,1) = pa_ao; dxb(:,:,1) = pb_ao
      call utddft_fxc(basis=p%infos%basis, molGrid=p%molgrid, isVecs=.true., &
                      wfa=p%moa, wfb=p%mob, fxa=ga, fxb=gb, dxa=dxa, dxb=dxb, &
                      nmtx=1, threshold=0.0d0, infos=p%infos)
      deallocate(dxa, dxb)
    end if

    ! Project back to MO occ-vir and add the orbital-energy diagonal.
    if (la > 0) call mntoia(ga(:,:,1), y(1:la), p%moa, p%moa, nocca, nocca)
    if (lb > 0) call mntoia(gb(:,:,1), y(la+1:la+lb), p%mob, p%mob, noccb, noccb)
    y = y + p%xm*x

    deallocate(pa_ao, pb_ao, ga, gb, dpack, fpack)
  end subroutine cphf_apbx_uhf

!###############################################################################

  subroutine cphf_precond_uhf(y, x, dat)
    real(kind=dp) :: x(:)
    real(kind=dp) :: y(:)
    type(c_ptr) :: dat
    type(cphf_cg_data_uhf), pointer :: p
    call c_f_pointer(dat, p)
    y = p%xminv*x
  end subroutine cphf_precond_uhf

!###############################################################################

  subroutine cphf_uhf_polarizability_selftest_C(c_handle) bind(C, name="cphf_uhf_polarizability_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call cphf_uhf_polarizability_selftest(inf)
  end subroutine cphf_uhf_polarizability_selftest_C

!> @brief Validate the open-shell CPHF solver via the static dipole
!>   polarizability.  Built per spin: B^sigma_ia = -<i|q|a>^sigma, solve
!>   M U^q = B^q, and alpha_pq = -2 sum_sigma sum_ia mu^p,sigma_ia U^q,sigma_ia.
!>   For a closed-shell system run as UHF (multiplicity 1) the tensor must equal
!>   the closed-shell (RHF) cphf_static_polarizability, which is the unambiguous
!>   correctness check for the spin coupling and normalization.  Written to
!>   /tmp/cphf_uhf_polar.out.
  subroutine cphf_uhf_polarizability_selftest(infos)
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A, OQP_VEC_MO_B
    use int1, only: multipole_integrals
    use mathlib, only: unpack_matrix
    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: moa(:,:), mob(:,:)
    real(kind=dp), allocatable :: mints(:,:), dipfull(:,:), dmo(:,:), scr(:,:)
    real(kind=dp), allocatable :: bvec(:,:), uvec(:,:), mua(:,:), mub(:,:)
    real(kind=dp) :: origin(3), alpha(3,3)
    integer :: nbf, nbf2, nocca, noccb, nvira, nvirb, la, lb, ltot
    integer :: q, i, a, ia, pq, uu

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    nvirb = nbf - noccb
    la = nocca*nvira
    lb = noccb*nvirb
    ltot = la + lb

    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, moa)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mob)

    allocate(mints(nbf2,19), source=0.0_dp)
    origin = 0.0_dp
    call multipole_integrals(basis, mints, origin, 3)

    allocate(dipfull(nbf,nbf), dmo(nbf,nbf), scr(nbf,nbf))
    allocate(bvec(ltot,3), uvec(ltot,3), source=0.0_dp)
    allocate(mua(la,3), mub(lb,3), source=0.0_dp)

    do q = 1, 3
      call unpack_matrix(mints(:,q), dipfull)
      ! alpha MO dipole and RHS
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, moa, nbf, dipfull, nbf, 0.0_dp, scr, nbf)
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, scr, nbf, moa, nbf, 0.0_dp, dmo, nbf)
      ia = 0
      do a = 1, nvira
        do i = 1, nocca
          ia = ia + 1
          mua(ia,q) = dmo(i, nocca+a)
          bvec(ia,q) = -dmo(i, nocca+a)
        end do
      end do
      ! beta MO dipole and RHS
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mob, nbf, dipfull, nbf, 0.0_dp, scr, nbf)
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, scr, nbf, mob, nbf, 0.0_dp, dmo, nbf)
      ia = 0
      do a = 1, nvirb
        do i = 1, noccb
          ia = ia + 1
          mub(ia,q) = dmo(i, noccb+a)
          bvec(la+ia,q) = -dmo(i, noccb+a)
        end do
      end do
    end do

    call cphf_solve_uhf(infos, 3, bvec, uvec)

    alpha = 0.0_dp
    do q = 1, 3
      do pq = 1, 3
        alpha(pq,q) = -2.0_dp*( sum(mua(:,pq)*uvec(1:la,q)) &
                              + sum(mub(:,pq)*uvec(la+1:ltot,q)) )
      end do
    end do

    open(newunit=uu, file='/tmp/cphf_uhf_polar.out', status='replace', action='write')
    write(uu,'(a)') 'open-shell (UHF) CPHF static dipole polarizability (a.u.):'
    do i = 1, 3
      write(uu,'(3f16.8)') alpha(i,1:3)
    end do
    write(uu,'(a,f16.8)') 'isotropic = ', (alpha(1,1)+alpha(2,2)+alpha(3,3))/3.0_dp
    close(uu)

    deallocate(mints, dipfull, dmo, scr, bvec, uvec, mua, mub)
  end subroutine cphf_uhf_polarizability_selftest

end module cphf_mod
