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

  !> Opaque data for the open-shell (ROHF) orbital-Hessian action.  ROHF uses a
  !> SINGLE MO set with a docc / socc / virt partition, so the rotation vector is
  !> laid out over the three non-redundant blocks (socc-docc, virt-docc,
  !> virt-socc) exactly as scf_converger::pack_rohf_trial.  The action mirrors
  !> the validated TRAH ROHF orbital Hessian (scf_converger::calc_h_op):
  !>   y = pack( Fvv^a xa - xa Foo^a + [C^a^T G^a C^a]_vo ,
  !>             Fvv^b xb - xb Foo^b + [C^b^T G^b C^b]_vo )
  !> with Foo/Fvv the occ-occ / vir-vir blocks of the converged spin Fock
  !> matrices in the MO basis and G^s the response Fock from get_response_packed.
  type :: cphf_cg_data_rohf
    type(information), pointer :: infos => null()
    type(basis_set), pointer :: basis => null()
    type(dft_grid_t), pointer :: molgrid => null()
    real(kind=dp), pointer :: mo(:,:) => null()
    real(kind=dp), pointer :: foo(:,:) => null()     ! alpha occ-occ Fock (MO)
    real(kind=dp), pointer :: fvv(:,:) => null()     ! alpha vir-vir Fock (MO)
    real(kind=dp), pointer :: foo_b(:,:) => null()   ! beta  occ-occ Fock (MO)
    real(kind=dp), pointer :: fvv_b(:,:) => null()   ! beta  vir-vir Fock (MO)
    real(kind=dp), pointer :: xminv(:) => null()     ! diagonal preconditioner
    integer :: nbf = 0
    integer :: nocca = 0, noccb = 0, nvira = 0, nvirb = 0, offset = 0, ltot = 0
    real(kind=dp) :: scale_exch = 1.0_dp
    logical :: dft = .false.
  end type

  private
  public :: cphf_solve
  public :: cphf_solve_uhf
  public :: cphf_solve_rohf
  public :: rohf_pack_trial, rohf_unpack_trial
  public :: cphf_static_polarizability
  public :: cphf_static_polarizability_C
  public :: cphf_polarizability_selftest
  public :: cphf_polarizability_selftest_C
  public :: cphf_uhf_polarizability_selftest
  public :: cphf_uhf_polarizability_selftest_C
  public :: cphf_rohf_polarizability_selftest
  public :: cphf_rohf_polarizability_selftest_C

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

!###############################################################################
!  Open-shell (ROHF) CPHF solver
!###############################################################################

!> @brief Pack ROHF alpha/beta vir-occ rotation matrices into a single vector.
!>   Layout (nocc_a >= nocc_b, offset = nocc_a - nocc_b = n_socc):
!>     block 1 (socc-docc): xb(1:offset, 1:noccb)
!>     block 2 (virt-docc): xa(1:nvira, 1:noccb) + xb(offset+1:, 1:noccb)
!>     block 3 (virt-socc): xa(1:nvira, noccb+1:nocca)
!>   Mirrors scf_converger::pack_rohf_trial.
  subroutine rohf_pack_trial(x, xa, xb, nbf, nocca, noccb)
    real(kind=dp), intent(out) :: x(:)
    real(kind=dp), intent(in)  :: xa(:,:), xb(:,:)
    integer, intent(in) :: nbf, nocca, noccb
    integer :: nvira, offset, k, iv, a
    nvira = nbf - nocca
    offset = nocca - noccb
    x = 0.0_dp
    k = 0
    if (offset > 0) then
      do iv = 1, offset
        do a = 1, noccb
          k = k + 1; x(k) = xb(iv, a)
        end do
      end do
    end if
    do iv = 1, nvira
      do a = 1, noccb
        k = k + 1; x(k) = xa(iv, a) + xb(offset + iv, a)
      end do
    end do
    if (offset > 0) then
      do iv = 1, nvira
        do a = 1, offset
          k = k + 1; x(k) = xa(iv, noccb + a)
        end do
      end do
    end if
  end subroutine rohf_pack_trial

!> @brief Inverse of rohf_pack_trial (scf_converger::unpack_rohf_trial).
  subroutine rohf_unpack_trial(x, xa, xb, nbf, nocca, noccb)
    real(kind=dp), intent(in)  :: x(:)
    real(kind=dp), intent(out) :: xa(:,:), xb(:,:)
    integer, intent(in) :: nbf, nocca, noccb
    integer :: nvira, offset, k, iv, a
    nvira = nbf - nocca
    offset = nocca - noccb
    xa = 0.0_dp; xb = 0.0_dp
    k = 0
    if (offset > 0) then
      do iv = 1, offset
        do a = 1, noccb
          k = k + 1; xb(iv, a) = x(k)
        end do
      end do
    end if
    do iv = 1, nvira
      do a = 1, noccb
        k = k + 1
        xa(iv, a) = x(k)
        xb(offset + iv, a) = x(k)
      end do
    end do
    if (offset > 0) then
      do iv = 1, nvira
        do a = 1, offset
          k = k + 1; xa(iv, noccb + a) = x(k)
        end do
      end do
    end if
  end subroutine rohf_unpack_trial

!> @brief Solve the open-shell (ROHF) CPHF equations  H theta = B  over the
!>   docc/socc/virt rotation space (layout: rohf_pack_trial).  The orbital
!>   Hessian action replicates the validated TRAH ROHF operator
!>   (scf_converger::calc_h_op): per spin the orbital-energy-difference part
!>   Fvv x - x Foo (full MO Fock blocks, so non-canonical orbitals are handled)
!>   plus the response Fock from the trial rotation density (get_response_packed,
!>   scftype>=2 -> Coulomb from the spin-summed density, exchange same-spin).
  subroutine cphf_solve_rohf(infos, nrhs, bvec, uvec, tol, maxit)
    use oqp_tagarray_driver, only: tagarray_get_data, &
        OQP_VEC_MO_A, OQP_FOCK_A, OQP_FOCK_B
    use mathlib, only: unpack_matrix
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
    type(cphf_cg_data_rohf), target :: cgdata
    type(pcg_t) :: pcg

    real(kind=dp), contiguous, pointer :: mo(:,:), focka(:), fockb(:)
    real(kind=dp), allocatable, target :: foo(:,:), fvv(:,:), foo_b(:,:), fvv_b(:,:)
    real(kind=dp), allocatable, target :: xminv(:)
    real(kind=dp), allocatable :: fao(:,:), w2(:,:), w3(:,:)
    integer :: nbf, nocca, noccb, nvira, nvirb, offset, ltot
    integer :: i, a, k, irhs, iter, mxit
    logical :: dft
    real(kind=dp) :: cnv, scale_exch, d

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    nvirb = nbf - noccb
    offset = nocca - noccb
    ltot = noccb*(offset + nvira) + offset*nvira
    dft = infos%control%hamilton == 20
    cnv = default_tol; if (present(tol)) cnv = tol
    mxit = 100; if (present(maxit)) mxit = maxit
    if (mxit < ltot + 5) mxit = ltot + 5

    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, focka)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fockb)

    if (dft) call dft_initialize(infos, basis, molGrid)

    ! MO-basis occ-occ / vir-vir blocks of the converged spin Fock matrices
    allocate(foo(nocca,nocca), fvv(nvira,nvira), foo_b(noccb,noccb), fvv_b(nvirb,nvirb))
    allocate(fao(nbf,nbf), w2(nbf,nbf), w3(nbf,nbf))
    call unpack_matrix(focka, fao)
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, fao, nbf, mo, nbf, 0.0_dp, w2, nbf)
    call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mo, nbf, w2, nbf, 0.0_dp, w3, nbf)
    foo = w3(1:nocca,1:nocca); fvv = w3(nocca+1:nbf,nocca+1:nbf)
    call unpack_matrix(fockb, fao)
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, fao, nbf, mo, nbf, 0.0_dp, w2, nbf)
    call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mo, nbf, w2, nbf, 0.0_dp, w3, nbf)
    foo_b = w3(1:noccb,1:noccb); fvv_b = w3(noccb+1:nbf,noccb+1:nbf)

    ! diagonal preconditioner (orbital-energy-difference gaps, response neglected)
    allocate(xminv(ltot))
    k = 0
    if (offset > 0) then
      do i = 1, offset                       ! socc-docc
        do a = 1, noccb
          k = k + 1; d = fvv_b(i,i) - foo_b(a,a); xminv(k) = 1.0_dp/sign(max(abs(d),1.0d-8), d)
        end do
      end do
    end if
    do i = 1, nvira                          ! virt-docc (alpha + beta share)
      do a = 1, noccb
        k = k + 1
        d = (fvv(i,i) - foo(a,a)) + (fvv_b(offset+i,offset+i) - foo_b(a,a))
        xminv(k) = 1.0_dp/sign(max(abs(d),1.0d-8), d)
      end do
    end do
    if (offset > 0) then
      do i = 1, nvira                        ! virt-socc
        do a = 1, offset
          k = k + 1; d = fvv(i,i) - foo(noccb+a,noccb+a); xminv(k) = 1.0_dp/sign(max(abs(d),1.0d-8), d)
        end do
      end do
    end if

    scale_exch = 1.0_dp
    if (dft) scale_exch = infos%dft%HFscale

    cgdata%infos => infos
    cgdata%basis => basis
    cgdata%molgrid => molgrid
    cgdata%mo => mo
    cgdata%foo => foo; cgdata%fvv => fvv
    cgdata%foo_b => foo_b; cgdata%fvv_b => fvv_b
    cgdata%xminv => xminv
    cgdata%nbf = nbf
    cgdata%nocca = nocca; cgdata%noccb = noccb
    cgdata%nvira = nvira; cgdata%nvirb = nvirb
    cgdata%offset = offset; cgdata%ltot = ltot
    cgdata%scale_exch = scale_exch
    cgdata%dft = dft

    write(iw,'(/3x,60("-"))')
    write(iw,'(6x,"open-shell (ROHF) CPHF iterative solver")')
    write(iw,'(6x,"right-hand sides =",I5,3x,"rotation dim =",I6)') nrhs, ltot
    write(iw,'(6x,"tolerance =",1P,E10.3,3x,"max iterations =",I6)') cnv, mxit
    write(iw,'(3x,60("-"))')

    do irhs = 1, nrhs
      call pcg%init(b=bvec(:,irhs), update=cphf_apbx_rohf, precond=cphf_precond_rohf, &
                    dat=cgdata, tol=sqrt(abs(cnv)))
      do iter = 1, mxit
        if (pcg%errcode /= PCG_OK) exit
        call pcg%step()
      end do
      write(iw,'(" ROHF CPHF RHS",I5," completed in",I5," iterations; error =",1P,E10.3)') &
              irhs, iter - 1, pcg%error**2
      call flush(iw)
      uvec(:,irhs) = pcg%x
      call pcg%clean()
    end do

    deallocate(foo, fvv, foo_b, fvv_b, xminv, fao, w2, w3)
  end subroutine cphf_solve_rohf

!###############################################################################

!> @brief ROHF orbital-Hessian action y = H x (see cphf_solve_rohf).
  subroutine cphf_apbx_rohf(y, x, dat)
    use mathlib, only: pack_matrix, unpack_matrix
    use scf_addons, only: get_response_packed
    real(kind=dp) :: x(:)
    real(kind=dp) :: y(:)
    type(c_ptr) :: dat
    type(cphf_cg_data_rohf), pointer :: p

    real(kind=dp), allocatable :: xa(:,:), xb(:,:), x2a(:,:), x2b(:,:)
    real(kind=dp), allocatable :: work2(:,:), work3(:,:), dm(:,:), v(:,:)
    real(kind=dp), allocatable :: dm_tri(:,:), pfock(:,:)
    integer :: nbf, nbf2, nocca, noccb, nvira, nvirb, i, j

    call c_f_pointer(dat, p)
    nbf = p%nbf; nbf2 = nbf*(nbf+1)/2
    nocca = p%nocca; noccb = p%noccb; nvira = p%nvira; nvirb = p%nvirb

    allocate(xa(nvira,nocca), xb(nvirb,noccb), x2a(nvira,nocca), x2b(nvirb,noccb))
    allocate(work2(nbf,nbf), work3(nbf,nbf), dm(nbf,nbf), v(nbf,nbf))
    allocate(dm_tri(nbf2,2), pfock(nbf2,2), source=0.0_dp)

    call rohf_unpack_trial(x, xa, xb, nbf, nocca, noccb)

    ! orbital-energy-difference part:  Fvv x - x Foo  (per spin)
    call dgemm('n','n', nvira, nocca, nvira, 1.0_dp, p%fvv, nvira, xa, nvira, 0.0_dp, x2a, nvira)
    call dgemm('n','n', nvira, nocca, nocca, -1.0_dp, xa, nvira, p%foo, nocca, 1.0_dp, x2a, nvira)
    call dgemm('n','n', nvirb, noccb, nvirb, 1.0_dp, p%fvv_b, nvirb, xb, nvirb, 0.0_dp, x2b, nvirb)
    call dgemm('n','n', nvirb, noccb, noccb, -1.0_dp, xb, nvirb, p%foo_b, noccb, 1.0_dp, x2b, nvirb)

    ! orbital-rotation density (alpha):  dm = Cv xa Co^T + (Cv xa Co^T)^T
    work2 = 0.0_dp
    call dgemm('n','n', nbf, nocca, nvira, 1.0_dp, p%mo(:,nocca+1:nbf), nbf, xa, nvira, 0.0_dp, work2, nbf)
    call dgemm('n','t', nbf, nbf, nocca, 1.0_dp, work2, nbf, p%mo(:,1:nocca), nbf, 0.0_dp, work3, nbf)
    do i = 1, nbf
      do j = 1, nbf
        dm(i,j) = work3(i,j) + work3(j,i)
      end do
    end do
    call pack_matrix(dm, dm_tri(:,1))
    ! beta
    work2 = 0.0_dp
    call dgemm('n','n', nbf, noccb, nvirb, 1.0_dp, p%mo(:,noccb+1:nbf), nbf, xb, nvirb, 0.0_dp, work2, nbf)
    call dgemm('n','t', nbf, nbf, noccb, 1.0_dp, work2, nbf, p%mo(:,1:noccb), nbf, 0.0_dp, work3, nbf)
    do i = 1, nbf
      do j = 1, nbf
        dm(i,j) = work3(i,j) + work3(j,i)
      end do
    end do
    call pack_matrix(dm, dm_tri(:,2))

    ! response Fock from the trial density (open-shell: J[dPa+dPb] - cx K[dP^s])
    call get_response_packed(p%basis, p%infos, p%molgrid, p%mo, dm_tri, pfock, p%mo)

    ! add the MO vir-occ block of the response Fock (alpha)
    call unpack_matrix(pfock(:,1), v)
    call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, p%mo, nbf, v, nbf, 0.0_dp, work2, nbf)
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, work2, nbf, p%mo, nbf, 0.0_dp, work3, nbf)
    x2a = x2a + work3(nocca+1:nbf, 1:nocca)
    ! beta
    call unpack_matrix(pfock(:,2), v)
    call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, p%mo, nbf, v, nbf, 0.0_dp, work2, nbf)
    call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, work2, nbf, p%mo, nbf, 0.0_dp, work3, nbf)
    x2b = x2b + work3(noccb+1:nbf, 1:noccb)

    call rohf_pack_trial(y, x2a, x2b, nbf, nocca, noccb)

    deallocate(xa, xb, x2a, x2b, work2, work3, dm, v, dm_tri, pfock)
  end subroutine cphf_apbx_rohf

!###############################################################################

  subroutine cphf_precond_rohf(y, x, dat)
    real(kind=dp) :: x(:)
    real(kind=dp) :: y(:)
    type(c_ptr) :: dat
    type(cphf_cg_data_rohf), pointer :: p
    call c_f_pointer(dat, p)
    y = p%xminv*x
  end subroutine cphf_precond_rohf

!###############################################################################

  subroutine cphf_rohf_polarizability_selftest_C(c_handle) bind(C, name="cphf_rohf_polarizability_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call cphf_rohf_polarizability_selftest(inf)
  end subroutine cphf_rohf_polarizability_selftest_C

!> @brief Validate the ROHF CPHF solver via the static dipole polarizability.
!>   For a closed-shell molecule run as ROHF (multiplicity 1, offset=0) the
!>   rotation space reduces to the virt-docc block and the ROHF orbital Hessian
!>   reduces to (twice) the RHF one; the resulting static polarizability must
!>   equal the validated closed-shell cphf_static_polarizability.  This is the
!>   unambiguous check for the solver plumbing, the operator and the packing.
!>   Written to /tmp/cphf_rohf_polar.out.
  subroutine cphf_rohf_polarizability_selftest(infos)
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A
    use int1, only: multipole_integrals
    use mathlib, only: unpack_matrix
    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: mo(:,:)
    real(kind=dp), allocatable :: mints(:,:), dipfull(:,:), dmo(:,:), scr(:,:)
    real(kind=dp), allocatable :: xa(:,:), xb(:,:), bvec(:,:), uvec(:,:)
    real(kind=dp) :: origin(3), alpha(3,3)
    integer :: nbf, nbf2, nocca, noccb, nvira, nvirb, offset, ltot
    integer :: q, pq, i, a, uu

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    nvirb = nbf - noccb
    offset = nocca - noccb
    ltot = noccb*(offset + nvira) + offset*nvira

    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)

    allocate(mints(nbf2,19), source=0.0_dp)
    origin = 0.0_dp
    call multipole_integrals(basis, mints, origin, 3)

    allocate(dipfull(nbf,nbf), dmo(nbf,nbf), scr(nbf,nbf))
    allocate(xa(nvira,nocca), xb(nvirb,noccb))
    allocate(bvec(ltot,3), uvec(ltot,3), source=0.0_dp)

    ! dipole RHS over the rotation space (single ROHF MO set; vir-occ blocks)
    do q = 1, 3
      call unpack_matrix(mints(:,q), dipfull)
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mo, nbf, dipfull, nbf, 0.0_dp, scr, nbf)
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, scr, nbf, mo, nbf, 0.0_dp, dmo, nbf)
      do i = 1, nocca
        do a = 1, nvira
          xa(a,i) = -dmo(nocca+a, i)
        end do
      end do
      do i = 1, noccb
        do a = 1, nvirb
          xb(a,i) = -dmo(noccb+a, i)
        end do
      end do
      call rohf_pack_trial(bvec(:,q), xa, xb, nbf, nocca, noccb)
    end do

    call cphf_solve_rohf(infos, 3, bvec, uvec)

    ! alpha_pq = -2 sum over rotation space of mu^p . theta^q  (mu = -bvec)
    alpha = 0.0_dp
    do q = 1, 3
      do pq = 1, 3
        alpha(pq,q) = -2.0_dp*sum( (-bvec(:,pq)) * uvec(:,q) )
      end do
    end do

    open(newunit=uu, file='/tmp/cphf_rohf_polar.out', status='replace', action='write')
    write(uu,'(a)') 'open-shell (ROHF) CPHF static dipole polarizability (a.u.):'
    do i = 1, 3
      write(uu,'(3f16.8)') alpha(i,1:3)
    end do
    write(uu,'(a,f16.8)') 'isotropic = ', (alpha(1,1)+alpha(2,2)+alpha(3,3))/3.0_dp
    close(uu)

    deallocate(mints, dipfull, dmo, scr, xa, xb, bvec, uvec)
  end subroutine cphf_rohf_polarizability_selftest

end module cphf_mod
