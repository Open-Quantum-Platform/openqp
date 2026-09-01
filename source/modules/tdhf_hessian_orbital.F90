module tdhf_hessian_orbital_mod

  use precision, only: dp

  implicit none

  private
  public :: build_tdhf_ground_orbital_response
  logical, parameter :: enable_tddft_orbital_xc = .true.
  logical, parameter :: enable_tddft_canonical_fxc_response = .true.

contains

!###############################################################################

  subroutine build_tdhf_ground_orbital_response(infos, sx_mo, umat, eps_deriv, dp_ao, &
                                                  dxc_skeleton)
    ! Closed-shell HF orbital response to all Cartesian nuclear perturbations.
    ! This supplies the S^K, U^K, epsilon^K, and P^K quantities entering the
    ! differentiated TD response operators.  DFT requires the additional
    ! moving-grid XC derivative and is rejected until that term is supplied.

    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, &
      OQP_VEC_MO_A, OQP_E_MO_A
    use mathlib, only: unpack_matrix, pack_matrix
    use grd1, only: der_overlap_matrix, der_kinetic_matrix, der_nucattr_matrix
    use ecp_tool, only: ecp_deriv_ints
    use fock_deriv_mod, only: fock_deriv_matrix
    use scf_addons, only: fock_jk
    use cphf_mod, only: cphf_solve
    use tdhf_hessian_response_mod, only: complete_rhf_orbital_response, &
      tdhf_reference_has_degenerate_subspace
    use messages, only: show_message, WITH_ABORT

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(out) :: sx_mo(:,:,:), umat(:,:,:)
    real(kind=dp), intent(out) :: eps_deriv(:,:), dp_ao(:,:,:)
    real(kind=dp), intent(out) :: dxc_skeleton(:,:,:)

    type(basis_set), pointer :: basis
    real(kind=dp), contiguous, pointer :: dmat(:), mo(:,:), eps(:)
    real(kind=dp), allocatable, target :: pfull(:,:)
    real(kind=dp), allocatable :: ds(:,:,:,:), dt(:,:,:,:), dv(:,:,:,:), decp(:,:,:,:)
    real(kind=dp), allocatable :: dg(:,:,:,:), f0ao(:,:), f0mo(:,:), gd0mo(:,:)
    real(kind=dp), allocatable :: d0(:,:), d0pack(:,:), gpack(:,:), gfull(:,:)
    real(kind=dp), allocatable :: probe1(:,:), bvec(:,:), uvec(:,:)
    real(kind=dp), allocatable :: dmo_occ(:,:)
    real(kind=dp), allocatable :: dxc_explicit(:,:,:)
    integer :: a, c, i, ia, k, mu, nbf, nbf2, natom, ncart, nocc, nvir, p, q, status
    real(kind=dp) :: scale_exch
    real(kind=dp), parameter :: xc_deriv_step = 1.0e-3_dp
    logical :: converged

    if (infos%control%scftype /= 1 .or. infos%mol_prop%mult /= 1) then
      call show_message('TD Hessian orbital response requires a closed-shell restricted reference.', &
                        WITH_ABORT)
    end if

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    natom = size(basis%atoms%xyz,2)
    ncart = 3*natom
    nocc = infos%mol_prop%nocc
    nvir = nbf-nocc
    scale_exch = 1.0_dp
    if (infos%control%hamilton >= 20) scale_exch = infos%dft%hfscale
    if (any(shape(sx_mo) /= [nbf,nbf,ncart]) .or. &
        any(shape(umat) /= [nbf,nbf,ncart]) .or. &
        any(shape(eps_deriv) /= [nbf,ncart]) .or. &
        any(shape(dp_ao) /= [nbf,nbf,ncart]) .or. &
        any(shape(dxc_skeleton) /= [nbf,nbf,ncart])) then
      call show_message('TD Hessian orbital-response output has the wrong shape.', WITH_ABORT)
    end if

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, eps)
    if (tdhf_reference_has_degenerate_subspace(eps, nocc, 1.0e-10_dp)) then
      call show_message('Analytic TD Hessian does not yet support degenerate occupied '// &
                        'or virtual canonical MO subspaces; use a numerical Hessian.', WITH_ABORT)
    end if
    allocate(pfull(nbf,nbf)); call unpack_matrix(dmat, pfull)
    dxc_skeleton = 0.0_dp
    allocate(ds(nbf,nbf,3,natom), dt(nbf,nbf,3,natom), &
             dv(nbf,nbf,3,natom), decp(nbf,nbf,3,natom), &
             dg(nbf,nbf,3,natom))
    call der_overlap_matrix(basis, ds)
    call der_kinetic_matrix(basis, dt)
    call der_nucattr_matrix(basis, basis%atoms%xyz, &
      basis%atoms%zn-basis%ecp_zn_num, dv)
    do k = 1, natom
      do c = 1, 3
        do i = 1, nbf
          do mu = 1, nbf
            ds(mu,i,c,k) = ds(mu,i,c,k)*basis%bfnrm(mu)*basis%bfnrm(i)
            dt(mu,i,c,k) = dt(mu,i,c,k)*basis%bfnrm(mu)*basis%bfnrm(i)
            dv(mu,i,c,k) = dv(mu,i,c,k)*basis%bfnrm(mu)*basis%bfnrm(i)
          end do
        end do
      end do
    end do
    call ecp_deriv_ints(basis, basis%atoms%xyz, decp)
    dv = dv + decp
    call fock_deriv_matrix(infos, basis, pfull, scale_exch, dg)

    allocate(f0ao(nbf,nbf), f0mo(nbf,nbf), gd0mo(nbf,nbf), &
             d0(nbf,nbf), d0pack(nbf2,1), gpack(nbf2,1), gfull(nbf,nbf), &
             probe1(nbf,nbf), bvec(nocc*nvir,ncart), &
             uvec(nocc*nvir,ncart), dmo_occ(nbf,nocc), &
             dxc_explicit(nbf,nbf,ncart), source=0.0_dp)
    ia = 0
    do k = 1, natom
      do c = 1, 3
        ia = ia + 1
        call ao_to_mo(ds(:,:,c,k), mo, sx_mo(:,:,ia), probe1)
        f0ao = dt(:,:,c,k) + dv(:,:,c,k) + dg(:,:,c,k)
        if (infos%control%hamilton == 20 .and. enable_tddft_orbital_xc) then
          block
            use mod_dft, only: dft_initialize, dftclean, dftexcor
            use mod_dft_molgrid, only: dft_grid_t
            type(dft_grid_t) :: grid
            real(dp),allocatable :: dcorb(:,:),mop(:,:),fp(:),fm(:),dvxc(:,:)
            real(dp) :: exc,tele,tkin
            integer :: io,jo
            allocate(dcorb(nbf,nocc),mop(nbf,nbf),fp(nbf2),fm(nbf2),dvxc(nbf,nbf))
            dcorb=0.0_dp
            do io=1,nocc
              do jo=1,nocc
                dcorb(:,io)=dcorb(:,io)-0.5_dp*mo(:,jo)*sx_mo(jo,io,ia)
              end do
            end do
            basis%atoms%xyz(c,k)=basis%atoms%xyz(c,k)+xc_deriv_step
            call basis%init_shell_centers(); call dft_initialize(infos,basis,grid)
            mop=mo; mop(:,1:nocc)=mo(:,1:nocc)+xc_deriv_step*dcorb; fp=0.0_dp
            call dftexcor(basis,grid,1,fp,fp,mop,mop,nbf,nbf2,exc,tele,tkin,infos)
            call dftclean(infos)
            basis%atoms%xyz(c,k)=basis%atoms%xyz(c,k)-2.0_dp*xc_deriv_step
            call basis%init_shell_centers(); call dft_initialize(infos,basis,grid)
            mop=mo; mop(:,1:nocc)=mo(:,1:nocc)-xc_deriv_step*dcorb; fm=0.0_dp
            call dftexcor(basis,grid,1,fm,fm,mop,mop,nbf,nbf2,exc,tele,tkin,infos)
            call dftclean(infos)
            basis%atoms%xyz(c,k)=basis%atoms%xyz(c,k)+xc_deriv_step
            call basis%init_shell_centers()
            call unpack_symmetric((fp-fm)/(2.0_dp*xc_deriv_step),dvxc,nbf)
            f0ao=f0ao+dvxc
            dxc_explicit(:,:,ia)=dvxc
            deallocate(dcorb,mop,fp,fm,dvxc)
          end block
        end if
        call ao_to_mo(f0ao, mo, f0mo, probe1)
        d0 = -2.0_dp*matmul(matmul(mo(:,1:nocc), sx_mo(1:nocc,1:nocc,ia)), &
                            transpose(mo(:,1:nocc)))
        call pack_matrix(d0, d0pack(:,1))
        gpack = 0.0_dp
        call fock_jk(basis, d=d0pack, f=gpack, scale_exch=scale_exch, infos=infos)
        call unpack_symmetric(gpack(:,1), gfull, nbf)
        call ao_to_mo(gfull, mo, gd0mo, probe1)
        do a = 1, nvir
          do i = 1, nocc
            bvec((a-1)*nocc+i,ia) = -f0mo(i,nocc+a) &
              + eps(i)*sx_mo(i,nocc+a,ia) - gd0mo(i,nocc+a)
          end do
        end do
      end do
    end do

    call cphf_solve(infos, ncart, bvec, uvec, converged=converged)
    if (.not. converged) then
      call show_message('Ground-state CPHF did not converge for the TD Hessian.', WITH_ABORT)
    end if
    do ia = 1, ncart
      call complete_rhf_orbital_response(mo, nocc, sx_mo(:,:,ia), uvec(:,ia), &
        umat(:,:,ia), dmo_occ, dp_ao(:,:,ia), status)
      if (status /= 0) call show_message('Failed to complete the TD Hessian orbital response.', &
                                         WITH_ABORT)
      c = mod(ia-1,3)+1
      k = (ia-1)/3+1
      f0ao = dt(:,:,c,k) + dv(:,:,c,k) + dg(:,:,c,k)
      d0 = -2.0_dp*matmul(matmul(mo(:,1:nocc), sx_mo(1:nocc,1:nocc,ia)), &
                          transpose(mo(:,1:nocc)))
      call pack_matrix(dp_ao(:,:,ia), d0pack(:,1))
      gpack = 0.0_dp
      call fock_jk(basis, d=d0pack, f=gpack, scale_exch=scale_exch, infos=infos)
      call unpack_symmetric(gpack(:,1), gfull, nbf)
      if (infos%control%hamilton == 20 .and. enable_tddft_orbital_xc) then
        block
          use mod_dft, only: dft_initialize, dftclean
          use mod_dft_molgrid, only: dft_grid_t
          use mod_dft_gridint_fxc, only: tddft_fxc
          type(dft_grid_t) :: grid
          real(dp),allocatable,target :: dx(:,:,:)
          real(dp),allocatable :: fx(:,:,:)
          allocate(dx(nbf,nbf,2),fx(nbf,nbf,2),source=0.0_dp)
          dx(:,:,1)=0.5_dp*dp_ao(:,:,ia)
          dx(:,:,2)=0.5_dp*d0
          call dft_initialize(infos,basis,grid)
          call tddft_fxc(basis,grid,.true.,mo,fx,dx,2,1.0e-14_dp,infos)
          call dftclean(infos)
          if (enable_tddft_canonical_fxc_response) gfull=gfull+fx(:,:,1)
          f0ao=f0ao+dxc_explicit(:,:,ia)-fx(:,:,2)
          ! The one-electron response row contracts the relaxed excited-state
          ! density with HMO = dH + dVxc_skeleton.  The displaced-grid
          ! derivative above also contains the metric-density response; remove
          ! that fxc contribution here, exactly as for the CPKS skeleton.
          dxc_skeleton(:,:,ia)=dxc_explicit(:,:,ia)-fx(:,:,2)
          deallocate(dx,fx)
        end block
      end if
      call ao_to_mo(f0ao, mo, f0mo, probe1)
      call ao_to_mo(gfull, mo, gd0mo, probe1)
      f0mo = f0mo + gd0mo
      eps_deriv(:,ia) = [(f0mo(i,i)-eps(i)*sx_mo(i,i,ia), i=1,nbf)]
      ! Canonical orbitals fix the occupied-occupied and virtual-virtual
      ! rotations as well as the occupied-virtual response.  The former must
      ! not be replaced by the symmetric -S^R/2 gauge: differentiated TDHF
      ! matrices depend on the canonical frame within both subspaces.
      do q = 1, nbf
        do p = 1, nbf
          if (p /= q .and. abs(eps(p)-eps(q)) > 1.0e-10_dp) then
            umat(p,q,ia) = (sx_mo(p,q,ia)*eps(q)-f0mo(p,q))/(eps(p)-eps(q))
          end if
        end do
      end do
    end do

    ! A rigid translation cannot change the canonical Fock matrix, orbital
    ! energies, MO connection, or explicit XC operator.  Atom-centred grid
    ! differentiation leaves a small owner-atom residual in these quantities;
    ! restore the exact acoustic sum rule on the first atom.  This is applied
    ! only after every independently evaluated nuclear derivative is present.
    do c = 1, 3
      umat(:,:,c) = 0.0_dp
      eps_deriv(:,c) = 0.0_dp
      dxc_skeleton(:,:,c) = 0.0_dp
      do k = 2, natom
        ia = 3*(k-1)+c
        umat(:,:,c) = umat(:,:,c)-umat(:,:,ia)
        eps_deriv(:,c) = eps_deriv(:,c)-eps_deriv(:,ia)
        dxc_skeleton(:,:,c) = dxc_skeleton(:,:,c)-dxc_skeleton(:,:,ia)
      end do
    end do

    deallocate(pfull, ds, dt, dv, decp, dg, f0ao, f0mo, gd0mo, d0, &
      d0pack, gpack, gfull, probe1, bvec, uvec, dmo_occ, dxc_explicit)
  contains
    subroutine ao_to_mo(ao, coeff, mo_mat, work)
      real(kind=dp), intent(in) :: ao(:,:), coeff(:,:)
      real(kind=dp), intent(out) :: mo_mat(:,:), work(:,:)
      work = matmul(ao, coeff)
      mo_mat = matmul(transpose(coeff), work)
    end subroutine ao_to_mo

    subroutine unpack_symmetric(gpk, gfu, n)
      real(kind=dp), intent(in) :: gpk(:)
      real(kind=dp), intent(out) :: gfu(:,:)
      integer, intent(in) :: n
      integer :: ii, ij, jj
      ij = 0
      do ii = 1, n
        do jj = 1, ii
          ij = ij + 1
          gfu(ii,jj) = gpk(ij)
          gfu(jj,ii) = gpk(ij)
        end do
      end do
    end subroutine unpack_symmetric
  end subroutine build_tdhf_ground_orbital_response

end module tdhf_hessian_orbital_mod
