!> @brief Nuclear gradient from a FACTORIZED AO two-particle density.
!>
!> This module owns no theory.  It is the derivative-integral contraction that
!> the SA-CASSCF gradients of `pyoqp/oqp/library/casscf_sa_gradient.py` need
!> and that `grd1`/`grd2` already implement, exposed with a two-particle
!> density general enough for the relaxed (Z-vector) case:
!>
!>     dE/dx = sum_{mu nu} D_{mu nu} h^x_{mu nu}
!>           + 1/2 sum_{mu nu la si} Gamma_{mu nu la si} (mu nu|la si)^x
!>           - sum_{mu nu} X_{mu nu} S^x_{mu nu}
!>           + dV_NN/dx.
!>
!> `D` and `X` cross the boundary as ordinary symmetric AO matrices.  The
!> caller passes `X` itself (the energy-weighted density of the term
!> `-sum X S^x`); the sign flip the overlap-derivative kernel wants is applied
!> here, so no caller has to remember it.
!>
!> Why `Gamma` is factorized rather than dense
!> -------------------------------------------
!> `grd2` asks for a two-particle density block per shell quartet and assumes
!> the full eight-fold permutational symmetry of `(mu nu|la si)`; under
!> HARMONIC_ACTIVE it also needs the block in the Cartesian-effective basis,
!> which `build_cart_density` supplies for a symmetric MATRIX and for nothing
!> else.  A dense `nbf^4` AO tensor satisfies neither cheaply.  Every
!> two-particle density this stack builds is instead a short sum of products of
!> ordinary symmetric AO matrices, and that form is both manifestly eight-fold
!> symmetric and expandable matrix by matrix:
!>
!>     Gamma_{mu nu la si} =
!>         sum_t c_t [ P^t_{mu nu} P^t_{la si}
!>                     - 1/4 (P^t_{mu la} P^t_{nu si} + P^t_{nu la} P^t_{mu si}) ]
!>       + sum_k lam_k A^k_{mu nu} A^k_{la si}.
!>
!> The first family is the Hartree-Fock ("separable") form evaluated at an
!> arbitrary symmetric matrix, and is what carries every separable and
!> BILINEAR piece: a cross term `Gamma^cross(A, B)` -- the one-index transform
!> of a separable density, and the inactive-Fock folding of a transition
!> density -- is the polarization `Gamma^sep(A+B) - Gamma^sep(A) - Gamma^sep(B)`
!> and therefore three entries in this list, not a new kernel.  The second
!> family is the pure Coulomb-type product that the eigendecomposition of an
!> all-active correction produces, and it absorbs bilinear low-rank terms by
!> the same polarization.
!>
!> Everything about the FACTORS is the caller's business.  This module folds
!> `bfnrm` into each matrix, expands each to the Cartesian-effective basis when
!> HARMONIC_ACTIVE, contracts, and adds the one-electron, nuclear-repulsion and
!> Pulay terms.  The state-specific gradient of `casscf_gradient.F90` is the
!> `nsep = 1, c_1 = 1` special case of what is contracted here; it keeps its own
!> entry point so its numbers cannot move.
!>
!> The symmetry petite list is deliberately NOT used, for the reason spelled
!> out in `casscf_gradient.F90`: it is valid only for a totally symmetric
!> contracted density, and neither an individual SA root's density nor a
!> relaxed Z-vector density carries that guarantee.
module casscf_ao_gradient_mod

  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  use precision, only: dp
  use types, only: information
  use basis_tools, only: basis_set, bas_norm_matrix, build_cart_density
  use constants, only: HARMONIC_ACTIVE, NUM_CART_BF
  use grd2, only: grd2_driver, grd2_compute_data_t

  implicit none
  private

  integer, parameter :: i8 = c_int64_t

  public :: casscf_ao_gradient

  !> `info` slots handed back to the caller (all double):
  integer, parameter :: CAS_AG_NSEP  = 0  !< separable factors contracted
  integer, parameter :: CAS_AG_NVEC  = 1  !< low-rank factors contracted
  integer, parameter :: CAS_AG_NINFO = 2

  ! Status codes.  0 is success.
  integer(i8), parameter :: CASAG_ERR_INPUT = -30_i8  !< sizes disagree with the handle
  integer(i8), parameter :: CASAG_ERR_ALLOC = -31_i8  !< out of memory

  !> Two-particle density for the derivative-ERI driver.
  !>
  !> `sepm(t,:,:)` / `sepc(t)` are the separable factors and their
  !> coefficients; `av(k,:,:)` / `lam(k)` the Coulomb-type low-rank ones.  The
  !> `_cart` copies are the bfnrm-folded, Cartesian-expanded versions the
  !> derivative kernels contract under HARMONIC_ACTIVE.  The factor index is
  !> stored FIRST in both so the innermost contraction loop walks contiguous
  !> runs, exactly as `casscf_gradient.F90` stores its own.
  type, extends(grd2_compute_data_t) :: grd2_casao_compute_data_t
    real(kind=dp), allocatable :: sepm(:,:,:), sepc(:)
    real(kind=dp), allocatable :: av(:,:,:), lam(:)
    real(kind=dp), allocatable :: sepm_cart(:,:,:), av_cart(:,:,:)
    integer, allocatable :: cart_off(:)
    integer :: nbf = 0
    integer :: nbf_cart = 0
    integer :: nsep = 0
    integer :: nvec = 0
  contains
    procedure :: init => grd2_casao_init
    procedure :: clean => grd2_casao_clean
    procedure :: get_density => grd2_casao_get_density
    procedure :: build_cart => grd2_casao_build_cart
  end type

contains

!###############################################################################

  !> C-bound entry point.  One call produces the complete nuclear gradient in
  !> `infos%atoms%grad`.
  !>
  !> @param[in]  c_handle  the OQP handle (`infos%atoms%grad` is overwritten)
  !> @param[in]  nbf       basis size, checked against the handle
  !> @param[in]  dmat      one-particle AO density, C-order [nbf,nbf], symmetric
  !> @param[in]  xmat      energy-weighted AO density X of `-sum X S^x`,
  !>                       C-order [nbf,nbf], symmetric
  !> @param[in]  nsep      separable factors
  !> @param[in]  sepc      their coefficients, [nsep]
  !> @param[in]  sepm      their matrices, C-order [nsep,nbf,nbf], each symmetric
  !> @param[in]  nvec      low-rank Coulomb-type factors
  !> @param[in]  lam       their coefficients, [nvec]
  !> @param[in]  avm       their matrices, C-order [nvec,nbf,nbf], each symmetric
  !> @param[out] info      CAS_AG_NINFO diagnostics
  !> @return     0, or a negative status
  function casscf_ao_gradient_C(c_handle, nbf, dmat, xmat, nsep, sepc, sepm, &
                                nvec, lam, avm, info) result(status) &
      bind(C, name="casscf_ao_gradient")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    type(oqp_handle_t) :: c_handle
    integer(c_int32_t), value :: nbf, nsep, nvec
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: dmat(0:*), xmat(0:*)
    real(dp), intent(in) :: sepc(0:*), sepm(0:*)
    real(dp), intent(in) :: lam(0:*), avm(0:*)
    real(dp), intent(inout) :: info(0:*)
    integer(i8) :: status
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    status = casscf_ao_gradient(inf, int(nbf), dmat, xmat, int(nsep), sepc, &
                                sepm, int(nvec), lam, avm, info)
  end function casscf_ao_gradient_C

!###############################################################################

  !> The contraction itself, reachable from an in-process Fortran caller.
  function casscf_ao_gradient(infos, nbf, dmat, xmat, nsep, sepc, sepm, &
                              nvec, lam, avm, info) result(status)
    use io_constants, only: iw
    use grd1, only: print_gradient
    use util, only: measure_time
    use printing, only: print_module_info

    type(information), target, intent(inout) :: infos
    integer, intent(in) :: nbf, nsep, nvec
    real(dp), intent(in) :: dmat(0:*), xmat(0:*)
    real(dp), intent(in) :: sepc(0:*), sepm(0:*)
    real(dp), intent(in) :: lam(0:*), avm(0:*)
    real(dp), intent(inout) :: info(0:*)
    integer(i8) :: status

    type(basis_set), pointer :: basis
    type(grd2_casao_compute_data_t) :: gcomp
    real(dp), allocatable :: dm_packed(:), wm_packed(:)
    real(dp), allocatable :: sym(:,:)
    integer :: n, nbf_tri, k, p, q, ierr
    integer(i8) :: n2

    status = 0_i8
    do k = 0, CAS_AG_NINFO - 1
      info(k) = 0.0_dp
    end do

    n = infos%basis%nbf
    ! `nsep >= 1` is required, not merely usual: under HARMONIC_ACTIVE the
    ! shell-offset table `grd2` indexes with comes out of the factor expansion,
    ! so a request with no factors at all would have no basis to expand.  Every
    ! two-particle density this stack builds has at least the separable term.
    if (nbf /= n .or. n <= 0 .or. nsep < 1 .or. nvec < 0) then
      status = CASAG_ERR_INPUT
      return
    end if
    nbf_tri = n * (n + 1) / 2
    n2 = int(n, i8) * int(n, i8)

    open (unit=IW, file=infos%log_filename, position="append")
    call print_module_info('CASSCF_AO_Gradient', &
                           'Nuclear Gradient From A Factorized AO Density')

    ! ---- packed one-particle and energy-weighted densities.  The gradient
    ! term is `-sum X S^x` while the overlap-derivative kernel ADDS
    ! `sum W S^x`, so W = -X.  The negation lives here rather than in every
    ! caller.
    allocate(dm_packed(nbf_tri), wm_packed(nbf_tri), sym(n, n), stat=ierr)
    if (ierr /= 0) then
      status = CASAG_ERR_ALLOC
      close(iw)
      return
    end if
    do p = 0, n - 1
      do q = 0, n - 1
        sym(p + 1, q + 1) = 0.5_dp * (dmat(int(p, i8)*int(n, i8) + int(q, i8)) &
                                    + dmat(int(q, i8)*int(n, i8) + int(p, i8)))
      end do
    end do
    call pack_symmetric(n, sym, dm_packed)
    do p = 0, n - 1
      do q = 0, n - 1
        sym(p + 1, q + 1) = -0.5_dp * (xmat(int(p, i8)*int(n, i8) + int(q, i8)) &
                                     + xmat(int(q, i8)*int(n, i8) + int(p, i8)))
      end do
    end do
    call pack_symmetric(n, sym, wm_packed)
    deallocate(sym)

    ! ---- factorized two-particle density.  Each factor is symmetrized on the
    ! way in: `grd2` contracts one representative per permutation orbit, so a
    ! factor that is only symmetric to rounding would break the eight-fold
    ! symmetry the driver assumes.
    gcomp%nbf = n
    gcomp%nsep = nsep
    gcomp%nvec = nvec
    allocate(gcomp%sepc(max(nsep, 1)), gcomp%sepm(max(nsep, 1), n, n), &
             gcomp%lam(max(nvec, 1)), gcomp%av(max(nvec, 1), n, n), stat=ierr)
    if (ierr /= 0) then
      status = CASAG_ERR_ALLOC
      close(iw)
      return
    end if
    do k = 1, nsep
      gcomp%sepc(k) = sepc(k - 1)
      do p = 0, n - 1
        do q = 0, n - 1
          gcomp%sepm(k, p + 1, q + 1) = 0.5_dp * &
              (sepm(int(k - 1, i8)*n2 + int(p, i8)*int(n, i8) + int(q, i8)) &
             + sepm(int(k - 1, i8)*n2 + int(q, i8)*int(n, i8) + int(p, i8)))
        end do
      end do
    end do
    do k = 1, nvec
      gcomp%lam(k) = lam(k - 1)
      do p = 0, n - 1
        do q = 0, n - 1
          gcomp%av(k, p + 1, q + 1) = 0.5_dp * &
              (avm(int(k - 1, i8)*n2 + int(p, i8)*int(n, i8) + int(q, i8)) &
             + avm(int(k - 1, i8)*n2 + int(q, i8)*int(n, i8) + int(p, i8)))
        end do
      end do
    end do

    info(CAS_AG_NSEP) = real(nsep, dp)
    info(CAS_AG_NVEC) = real(nvec, dp)

    basis => infos%basis
    basis%atoms => infos%atoms

    write(iw,"(/' ..... Beginning CASSCF AO-Density Gradient ...'/)")
    write(iw,"(' Separable two-particle factors      ',I8)") nsep
    write(iw,"(' Low-rank two-particle factors       ',I8)") nvec
    call flush(iw)

    call casao_1e_grad(infos, basis, dm_packed, wm_packed)
    write(iw,"(' ..... End Of 1-Electron Gradient ......')")
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

    call casao_2e_grad(infos, basis, gcomp, status)
    if (status < 0_i8) then
      close(iw)
      return
    end if
    write(iw, fmt="(' ...... End Of 2-Electron Gradient ......')")
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

    call print_gradient(infos)
    call measure_time(print_total=1, log_unit=iw)
    close(iw)

    call gcomp%clean()
    deallocate(dm_packed, wm_packed)
  end function casscf_ao_gradient

!###############################################################################

  !> Pack the upper triangle of a symmetric matrix in the layout the gradient
  !> kernels read (`grd1::prepare_grad_density` unpacks it with the same
  !> default convention `mathlib::pack_matrix` writes).
  subroutine pack_symmetric(n, a, ap)
    use mathlib, only: pack_matrix
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n, n)
    real(dp), intent(out) :: ap(:)
    call pack_matrix(a, ap)
  end subroutine pack_symmetric

!###############################################################################

  !> Nuclear-repulsion, overlap (Pulay), kinetic, nuclear-attraction and ECP
  !> terms.  Structurally identical to `casscf_gradient_mod::cas_1e_grad`.
  subroutine casao_1e_grad(infos, basis, dm_packed, wm_packed)
    use constants, only: tol_int
    use grd1, only: grad_nn, grad_ee_overlap, grad_ee_kinetic, &
                    grad_en_hellman_feynman, grad_en_pulay, grad_1e_ecp

    type(information), intent(inout) :: infos
    type(basis_set), intent(inout) :: basis
    real(dp), intent(inout) :: dm_packed(:), wm_packed(:)

    real(dp) :: tol

    tol = tol_int * log(10.0_dp)

    associate(grad => infos%atoms%grad, &
              xyz  => infos%atoms%xyz, &
              zn   => infos%atoms%zn - infos%basis%ecp_zn_num)

      grad = 0.0_dp

      call grad_nn(infos%atoms, infos%basis%ecp_zn_num)
      call grad_ee_overlap(basis, wm_packed, grad, logtol=tol)
      call grad_ee_kinetic(basis, dm_packed, grad, logtol=tol)
      call grad_en_hellman_feynman(basis, xyz, zn, dm_packed, grad, logtol=tol)
      call grad_en_pulay(basis, xyz, zn, dm_packed, grad, logtol=tol)
      call grad_1e_ecp(infos, basis, xyz, dm_packed, grad, logtol=tol)

    end associate
  end subroutine casao_1e_grad

!###############################################################################

  !> Two-electron term: contract the factorized two-particle density against
  !> the derivative ERIs.
  subroutine casao_2e_grad(infos, basis, gcomp, status)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    type(grd2_casao_compute_data_t), intent(inout) :: gcomp
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: de(:,:)
    integer :: ierr

    allocate(de(3, ubound(infos%atoms%zn, 1)), source=0.0_dp, stat=ierr)
    if (ierr /= 0) then
      status = CASAG_ERR_ALLOC
      return
    end if

    ! Neither family of factors is a Coulomb or an exchange contribution in the
    ! range-separated sense, so the scales are pinned and the CAM path -- which
    ! would run the driver twice with different ones -- is switched off
    ! explicitly.  CASSCF requires Hartree-Fock integrals anyway.
    gcomp%hfscale = 1.0_dp
    gcomp%coulscale = 1.0_dp
    call gcomp%init()
    call gcomp%build_cart(basis)
    call grd2_driver(infos, basis, de, gcomp, cam=.false.)
    infos%atoms%grad = infos%atoms%grad + de

    deallocate(de)
  end subroutine casao_2e_grad

!###############################################################################

  subroutine grd2_casao_init(this)
    class(grd2_casao_compute_data_t), target, intent(inout) :: this
    this%nbf_cart = 0
  end subroutine grd2_casao_init

  !> Fold `bfnrm` into every factor and, under HARMONIC_ACTIVE, expand each to
  !> the Cartesian-effective basis the derivative kernels iterate.  Every
  !> factor is an ordinary symmetric AO matrix, so each takes exactly the route
  !> a density does.
  subroutine grd2_casao_build_cart(this, basis)
    class(grd2_casao_compute_data_t), intent(inout) :: this
    type(basis_set), intent(in) :: basis

    real(dp), allocatable :: tmp(:,:), expanded(:,:)
    integer, allocatable :: off(:)
    integer :: k, ncart

    if (.not. HARMONIC_ACTIVE) return

    ! `build_cart_density`'s outputs are intent(out) allocatables, so each call
    ! reallocates them; the shell-offset table it rebuilds is a function of the
    ! basis alone and is kept from the last call.
    allocate(tmp(this%nbf, this%nbf))
    do k = 1, this%nsep
      tmp = this%sepm(k, :, :)
      call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
      call build_cart_density(basis, tmp, expanded, off, ncart)
      if (.not. allocated(this%sepm_cart)) then
        this%nbf_cart = ncart
        allocate(this%sepm_cart(this%nsep, ncart, ncart))
      end if
      this%sepm_cart(k, :, :) = expanded
    end do
    do k = 1, this%nvec
      tmp = this%av(k, :, :)
      call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
      call build_cart_density(basis, tmp, expanded, off, ncart)
      if (.not. allocated(this%av_cart)) then
        this%nbf_cart = ncart
        allocate(this%av_cart(this%nvec, ncart, ncart))
      end if
      this%av_cart(k, :, :) = expanded
    end do
    if (allocated(off)) call move_alloc(off, this%cart_off)
    deallocate(tmp)
    if (allocated(expanded)) deallocate(expanded)
  end subroutine grd2_casao_build_cart

  subroutine grd2_casao_clean(this)
    class(grd2_casao_compute_data_t), target, intent(inout) :: this
    if (allocated(this%sepm)) deallocate(this%sepm)
    if (allocated(this%sepc)) deallocate(this%sepc)
    if (allocated(this%av)) deallocate(this%av)
    if (allocated(this%lam)) deallocate(this%lam)
    if (allocated(this%sepm_cart)) deallocate(this%sepm_cart)
    if (allocated(this%av_cart)) deallocate(this%av_cart)
    if (allocated(this%cart_off)) deallocate(this%cart_off)
  end subroutine grd2_casao_clean

!###############################################################################

  !> Two-particle density for one shell quartet.
  !>
  !> Returns `4 x Gamma`, the scale the Hartree-Fock path already hands `grd2`
  !> (`coulfact = 4`, `xcfact = 1` there is `4 x [P P - 1/4 (P P + P P)]`), so
  !> both families carry the same factor and the caller's coefficients are the
  !> plain ones of the module header.
  subroutine grd2_casao_get_density(this, basis, id, dab, dabmax)

    class(grd2_casao_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: acc, bfn, cf
    logical :: usecart
    integer :: i, j, k, l, v
    integer :: loc(4), nbf(4)
    integer :: i1, j1, k1, l1
    real(kind=dp), pointer :: ab(:,:,:,:)
    real(kind=dp), pointer :: sm(:,:,:), av(:,:,:)

    usecart = HARMONIC_ACTIVE
    if (usecart) then
      loc = this%cart_off(id) - 1
      nbf = NUM_CART_BF(basis%am(id))
    else
      loc = basis%ao_offset(id) - 1
      nbf = basis%naos(id)
    end if
    sm => null()
    av => null()
    if (this%nsep > 0) then
      if (usecart) then
        sm => this%sepm_cart
      else
        sm => this%sepm
      end if
    end if
    if (this%nvec > 0) then
      if (usecart) then
        av => this%av_cart
      else
        av => this%av
      end if
    end if

    dabmax = 0
    ab(1:nbf(4),1:nbf(3),1:nbf(2),1:nbf(1)) => dab(1:product(nbf))

    do i = 1, nbf(1)
      i1 = loc(1) + i

      do j = 1, nbf(2)
        j1 = loc(2) + j

        do k = 1, nbf(3)
          k1 = loc(3) + k

          do l = 1, nbf(4)
            l1 = loc(4) + l
            acc = 0.0_dp
            ! The factor index is stored first in both stacks, so each inner
            ! loop walks contiguous runs; these are the hot loops of the build.
            do v = 1, this%nsep
              cf = this%sepc(v)
              acc = acc + cf * (4.0_dp*sm(v,i1,j1)*sm(v,k1,l1) &
                                - sm(v,i1,k1)*sm(v,j1,l1) &
                                - sm(v,i1,l1)*sm(v,j1,k1))
            end do
            do v = 1, this%nvec
              acc = acc + 4.0_dp * this%lam(v)*av(v,i1,j1)*av(v,k1,l1)
            end do
            dabmax = max(dabmax, abs(acc))
            bfn = 1.0_dp
            if (.not. usecart) bfn = product(basis%bfnrm([i1,j1,k1,l1]))
            ab(l,k,j,i) = acc*bfn
          end do
        end do
      end do
    end do
  end subroutine grd2_casao_get_density

end module casscf_ao_gradient_mod
