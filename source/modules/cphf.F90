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

  private
  public :: cphf_solve
  public :: cphf_static_polarizability
  public :: cphf_polarizability_selftest
  public :: cphf_polarizability_selftest_C

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
    integer :: nbf, nocc, nvir, lexc, i, j, irhs, iter, ok, mxit
    logical :: dft
    real(kind=dp) :: cnv, scale_exch

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

    do irhs = 1, nrhs
      call pcg%init(b=bvec(:,irhs), update=cphf_apbx, precond=cphf_precond, &
                    dat=cgdata, tol=sqrt(abs(cnv)))
      do iter = 1, mxit
        if (pcg%errcode /= PCG_OK) exit
        call pcg%step()
      end do
      uvec(:,irhs) = pcg%x
      call pcg%clean()
    end do

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

end module cphf_mod
