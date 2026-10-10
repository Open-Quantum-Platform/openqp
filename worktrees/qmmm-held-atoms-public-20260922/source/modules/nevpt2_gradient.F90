!> @file nevpt2_gradient.F90
!>
!> @brief Derivative-integral contraction for the analytic SC-NEVPT2 nuclear
!>        gradient.
!>
!> The SC-NEVPT2 Lagrangian, its orbital/CI response equations and the relaxed
!> densities are formed in `pyoqp/oqp/library/nevpt2_gradient.py`, which owns the
!> perturbation theory.  What is left is the part that needs derivative
!> integrals, and that is all this module does: given the AO-basis relaxed
!> one-particle density, the AO-basis energy-weighted density and the AO-basis
!> two-particle density of a STATIONARY Lagrangian,
!>
!>     dE/dx = sum_{mu nu} D_{mu nu} h^x_{mu nu}
!>           + 1/2 sum_{mu nu la si} Gamma_{mu nu la si} (mu nu|la si)^x
!>           - sum_{mu nu} X_{mu nu} S^x_{mu nu}
!>           + dV_NN/dx,
!>
!> it contracts them with the existing grd1/grd2 engines and writes
!> `infos%atoms%grad`.  No finite differences and no perturbation theory enter
!> here; the routine is a pure quadrature of supplied densities, and it is
!> pinned by feeding it the plain CASSCF densities and reproducing
!> `casscf_gradient` (see `tests/test_nevpt2_gradient.py`).
!>
!> Conventions (identical to `casscf_gradient.F90`)
!> ------------------------------------------------
!> * `dm_ao` and `x_ao` are full square symmetric AO matrices; they are packed
!>   here, and `x_ao` is NEGATED before it reaches `grad_ee_overlap`, which ADDS
!>   `sum W S^x` while the gradient term is `-sum X S^x`.
!> * `g2_ao` is the PHYSICAL two-particle density in the `E = 1/2 sum Gamma g`
!>   convention, carrying the full eight-fold permutational symmetry of a real
!>   AO ERI.  `grd2`'s unique-quartet convention consumes four times that, which
!>   is the factor applied in `get_density` below (the same factor
!>   `mp2_gradient` applies to its amplitude-linear density).
!> * Because `g2_ao` is eight-fold symmetric it is invariant under complete
!>   index reversal, `Gamma_{si la nu mu} = Gamma_{mu nu la si}`, so the caller's
!>   C-ordered buffer may be read directly as a Fortran-ordered rank-4 array.
!>   The caller symmetrizes; this module verifies rather than assumes it.
!>
!> The dense nbf^4 two-particle density is deliberate.  It keeps every index
!> convention independently auditable against the Lagrangian equations at the
!> validation-grade system sizes SC-NEVPT2 itself is limited to, and it can be
!> replaced by a factorized or streamed density without changing this entry
!> point's contract.
module nevpt2_gradient_mod

  use precision, only: dp
  use types, only: information
  use basis_tools, only: basis_set, bas_norm_matrix, build_cart_density
  use constants, only: HARMONIC_ACTIVE, NUM_CART_BF
  use grd2, only: grd2_driver, grd2_compute_data_t
  use iso_c_binding, only: c_int32_t, c_int64_t

  implicit none
  private

  public :: nevpt2_gradient

  !> Status codes returned to the Python driver.  Mirrors the NEVPT2_G_* block
  !> in include/oqp.h.
  integer(c_int64_t), parameter :: NEVPT2G_OK          =  0_c_int64_t
  integer(c_int64_t), parameter :: NEVPT2G_ERR_NBF     = -1_c_int64_t
  integer(c_int64_t), parameter :: NEVPT2G_ERR_ALLOC   = -2_c_int64_t
  integer(c_int64_t), parameter :: NEVPT2G_ERR_SYMM    = -3_c_int64_t

  !> Largest accepted deviation from eight-fold permutational symmetry in the
  !> supplied two-particle density, relative to its largest element.  The
  !> derivative-ERI driver assumes that symmetry, so a density lacking it would
  !> be contracted against the wrong integrals rather than merely inaccurately.
  real(dp), parameter :: SYMMETRY_TOL = 1.0e-9_dp

  type, extends(grd2_compute_data_t) :: grd2_nevpt2_compute_data_t
    real(dp), allocatable :: g2(:,:,:,:)
    real(dp), allocatable :: g2_cart(:,:,:,:)
    integer, allocatable :: cart_off(:)
    integer :: nbf = 0
    integer :: nbf_cart = 0
  contains
    procedure :: init => grd2_nevpt2_init
    procedure :: clean => grd2_nevpt2_clean
    procedure :: get_density => grd2_nevpt2_get_density
    procedure :: build_cart => grd2_nevpt2_build_cart
  end type grd2_nevpt2_compute_data_t

contains

!###############################################################################

  function nevpt2_gradient_C(c_handle, nbf, dm_ao, x_ao, g2_ao) &
      result(status) bind(C, name="nevpt2_gradient")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    type(oqp_handle_t) :: c_handle
    integer(c_int32_t), value :: nbf
    ! bind(C) hands over bare pointers: array dummies must be assumed-size.
    real(dp), intent(in) :: dm_ao(*), x_ao(*), g2_ao(*)
    integer(c_int64_t) :: status
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    status = nevpt2_gradient(inf, int(nbf), dm_ao, x_ao, g2_ao)
  end function nevpt2_gradient_C

!###############################################################################

  !> @brief Contract supplied AO relaxed densities with the derivative integrals.
  !>
  !> @param[in,out] infos   OpenQP information handle; `infos%atoms%grad` is
  !>                        overwritten with the completed gradient.
  !> @param[in]     nbf     number of spherical basis functions (must match the
  !>                        handle's basis).
  !> @param[in]     dm_ao   relaxed one-particle density, nbf x nbf, symmetric.
  !> @param[in]     x_ao    energy-weighted density (the orthonormality
  !>                        multiplier), nbf x nbf, symmetric.  Enters as
  !>                        `-sum X S^x`.
  !> @param[in]     g2_ao   relaxed two-particle density, nbf^4, eight-fold
  !>                        symmetric, in the `E = 1/2 sum Gamma g` convention.
  function nevpt2_gradient(infos, nbf, dm_ao, x_ao, g2_ao) result(status)
    use io_constants, only: iw
    use printing, only: print_module_info
    use constants, only: tol_int
    use mathlib, only: pack_matrix
    use grd1, only: print_gradient, grad_nn, grad_ee_overlap, &
                    grad_ee_kinetic, grad_en_hellman_feynman, &
                    grad_en_pulay, grad_1e_ecp

    type(information), target, intent(inout) :: infos
    integer, intent(in) :: nbf
    real(dp), intent(in) :: dm_ao(*), x_ao(*), g2_ao(*)
    integer(c_int64_t) :: status

    type(basis_set), pointer :: basis
    type(grd2_nevpt2_compute_data_t) :: gcomp
    real(dp), allocatable :: dsq(:,:), xsq(:,:), dpack(:), wpack(:), de2(:,:)
    integer :: n, p, q, r, s, ierr
    real(dp) :: tol, gmax, asym, ref

    status = NEVPT2G_OK
    basis => infos%basis
    basis%atoms => infos%atoms
    n = basis%nbf
    if (nbf /= n) then
      status = NEVPT2G_ERR_NBF
      return
    end if

    allocate(dsq(n,n), xsq(n,n), dpack(n*(n+1)/2), wpack(n*(n+1)/2), &
             gcomp%g2(n,n,n,n), source=0.0_dp, stat=ierr)
    if (ierr /= 0) then
      status = NEVPT2G_ERR_ALLOC
      return
    end if

    ! Symmetric square matrices are layout-independent, so the caller's C-order
    ! buffers are read directly; the explicit symmetrization below is what makes
    ! that safe rather than assumed.
    do q = 1, n
      do p = 1, n
        dsq(p,q) = dm_ao((p-1)*n + q)
        xsq(p,q) = x_ao((p-1)*n + q)
      end do
    end do
    dsq = 0.5_dp*(dsq + transpose(dsq))
    xsq = 0.5_dp*(xsq + transpose(xsq))

    do s = 1, n
      do r = 1, n
        do q = 1, n
          do p = 1, n
            gcomp%g2(p,q,r,s) = g2_ao(((((s-1)*n + (r-1))*n + (q-1))*n) + p)
          end do
        end do
      end do
    end do

    ! Verify the eight-fold symmetry the derivative-ERI driver relies on.
    gmax = maxval(abs(gcomp%g2))
    ref = max(gmax, 1.0_dp)
    asym = 0.0_dp
    do s = 1, n
      do r = 1, n
        do q = 1, n
          do p = 1, n
            asym = max(asym, abs(gcomp%g2(p,q,r,s) - gcomp%g2(q,p,r,s)))
            asym = max(asym, abs(gcomp%g2(p,q,r,s) - gcomp%g2(p,q,s,r)))
            asym = max(asym, abs(gcomp%g2(p,q,r,s) - gcomp%g2(r,s,p,q)))
          end do
        end do
      end do
    end do
    if (asym > SYMMETRY_TOL*ref) then
      status = NEVPT2G_ERR_SYMM
      call gcomp%clean()
      return
    end if
    gcomp%nbf = n

    call pack_matrix(dsq, dpack)
    ! `grad_ee_overlap` ADDS sum W S^x while the gradient term is -sum X S^x,
    ! so the negated multiplier is what goes in -- the same sign convention
    ! `grd1::eijden` uses for RHF and `casscf_gradient` reuses.
    xsq = -xsq
    call pack_matrix(xsq, wpack)

    open(unit=iw, file=infos%log_filename, position='append')
    call print_module_info('NEVPT2_Gradient', &
                           'Analytic SC-NEVPT2 Nuclear Gradient')
    write(iw,'(3X,A,I0)') 'dense AO 2-particle density dimension = ', n
    write(iw,'(3X,A,ES12.4)') 'two-particle density asymmetry        = ', asym

    tol = tol_int*log(10.0_dp)
    associate(grad => infos%atoms%grad, xyz => infos%atoms%xyz, &
              zn => infos%atoms%zn - infos%basis%ecp_zn_num)
      grad = 0.0_dp
      call grad_nn(infos%atoms, infos%basis%ecp_zn_num)
      call grad_ee_overlap(basis, wpack, grad, logtol=tol)
      call grad_ee_kinetic(basis, dpack, grad, logtol=tol)
      call grad_en_hellman_feynman(basis, xyz, zn, dpack, grad, logtol=tol)
      call grad_en_pulay(basis, xyz, zn, dpack, grad, logtol=tol)
      call grad_1e_ecp(infos, basis, xyz, dpack, grad, logtol=tol)
    end associate

    allocate(de2(3, ubound(infos%atoms%zn,1)), source=0.0_dp, stat=ierr)
    if (ierr /= 0) then
      status = NEVPT2G_ERR_ALLOC
      close(iw)
      call gcomp%clean()
      return
    end if
    gcomp%coulscale = 1.0_dp
    gcomp%hfscale = 1.0_dp
    call gcomp%init()
    call gcomp%build_cart(basis)
    ! The petite (symmetry) reduction is NOT opted into: it is valid only for a
    ! totally symmetric contracted density, and the relaxed two-particle density
    ! of an arbitrary `[casscf] root` carries no such guarantee.  This is the
    ! same choice `casscf_gradient` makes, and for the same reason.
    call grd2_driver(infos, basis, de2, gcomp, cam=.false.)
    infos%atoms%grad = infos%atoms%grad + de2

    call print_gradient(infos)
    close(iw)

    deallocate(de2)
    call gcomp%clean()
  end function nevpt2_gradient

!###############################################################################

  subroutine grd2_nevpt2_init(this)
    class(grd2_nevpt2_compute_data_t), target, intent(inout) :: this
    this%nbf_cart = 0
  end subroutine grd2_nevpt2_init

  subroutine grd2_nevpt2_clean(this)
    class(grd2_nevpt2_compute_data_t), target, intent(inout) :: this
    if (allocated(this%g2)) deallocate(this%g2)
    if (allocated(this%g2_cart)) deallocate(this%g2_cart)
    if (allocated(this%cart_off)) deallocate(this%cart_off)
  end subroutine grd2_nevpt2_clean

  !> Fold the basis normalization in and, under HARMONIC_ACTIVE, expand the
  !> two-particle density to the Cartesian-effective basis the derivative
  !> kernels iterate.  Structurally identical to `mp2_gradient`'s expansion:
  !> normalize all four spherical indices exactly once, then expand the AO pairs
  !> one after the other.
  subroutine grd2_nevpt2_build_cart(this, basis)
    class(grd2_nevpt2_compute_data_t), intent(inout) :: this
    type(basis_set), intent(in) :: basis
    real(dp), allocatable :: tmp(:,:), expanded(:,:), pair12(:,:,:,:)
    integer, allocatable :: off(:)
    integer :: nc, p, q, mu, nu

    if (.not. HARMONIC_ACTIVE) return

    ! The offsets and the Cartesian dimension come from expanding any symmetric
    ! AO matrix; the first ERI pair block serves.
    tmp = this%g2(:,:,1,1)
    call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
    call build_cart_density(basis, tmp, expanded, this%cart_off, this%nbf_cart)
    deallocate(expanded)

    allocate(pair12(this%nbf_cart,this%nbf_cart,this%nbf,this%nbf), &
             this%g2_cart(this%nbf_cart,this%nbf_cart,this%nbf_cart,this%nbf_cart), &
             source=0.0_dp)
    do p = 1, this%nbf
      do q = 1, this%nbf
        tmp = this%g2(:,:,p,q)
        call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
        tmp = tmp*basis%bfnrm(p)*basis%bfnrm(q)
        call build_cart_density(basis, tmp, expanded, off, nc)
        pair12(:,:,p,q) = expanded
        deallocate(expanded, off)
      end do
    end do
    do mu = 1, this%nbf_cart
      do nu = 1, this%nbf_cart
        tmp = pair12(mu,nu,:,:)
        call build_cart_density(basis, tmp, expanded, off, nc)
        this%g2_cart(mu,nu,:,:) = expanded
        deallocate(expanded, off)
      end do
    end do
    deallocate(pair12)
  end subroutine grd2_nevpt2_build_cart

  !> Two-particle density block for one shell quartet.  `grd2`'s unique-quartet
  !> convention consumes four times the physical symmetrized density.
  subroutine grd2_nevpt2_get_density(this, basis, id, dab, dabmax)
    class(grd2_nevpt2_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(dp), target, intent(out) :: dab(*)
    real(dp), intent(out) :: dabmax
    real(dp), pointer :: ab(:,:,:,:), g2(:,:,:,:)
    integer :: nb(4), loc(4), i,j,k,l,ii,jj,kk,ll
    real(dp) :: val, bfn

    if (HARMONIC_ACTIVE) then
      g2 => this%g2_cart
      loc = this%cart_off(id)-1
      nb = NUM_CART_BF(basis%am(id))
    else
      g2 => this%g2
      loc = basis%ao_offset(id)-1
      nb = basis%naos(id)
    end if
    ab(1:nb(4),1:nb(3),1:nb(2),1:nb(1)) => dab(1:product(nb))
    dabmax = 0.0_dp
    do i=1,nb(1); ii=loc(1)+i
      do j=1,nb(2); jj=loc(2)+j
        do k=1,nb(3); kk=loc(3)+k
          do l=1,nb(4); ll=loc(4)+l
            val = 4.0_dp*g2(ii,jj,kk,ll)
            bfn = 1.0_dp
            if (.not. HARMONIC_ACTIVE) bfn = product(basis%bfnrm([ii,jj,kk,ll]))
            ab(l,k,j,i) = val*bfn
            dabmax = max(dabmax, abs(val))
          end do
        end do
      end do
    end do
  end subroutine grd2_nevpt2_get_density

end module nevpt2_gradient_mod
