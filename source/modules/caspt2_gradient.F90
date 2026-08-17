!> @brief Derivative-integral contraction for the analytic CASPT2 /
!>        XMS-CASPT2 nuclear gradient.
!>
!> Division of labour
!> ------------------
!> `oqp.library.caspt2_gradient` (Python) owns the *theory*: it rebuilds the
!> CASPT2 state, solves the amplitude, CI, XMS-rotation and orbital-response
!> equations, and reduces the whole Lagrangian to three ordinary AO-basis
!> objects.  This module owns the *integrals*: it contracts those three objects
!> with the derivative one- and two-electron integrals and writes
!> `infos%atoms%grad`.  That split follows the CASPT2 energy path itself, where
!> the driver is Python and liboqp computes; derivative integrals exist only
!> here, so this is the part that has to be Fortran.
!>
!> The three objects
!> -----------------
!>     dE/dx = sum_{mu nu}      D^AO_{mu nu}    h^x_{mu nu}
!>           + 1/2 sum_{mu nu la si} Gamma^AO   (mu nu|la si)^x
!>           - sum_{mu nu}      X^AO_{mu nu}    S^x_{mu nu}
!>           + dV_NN/dx
!>
!>   * `dm_packed` -- the relaxed one-particle density `D^AO`, packed triangle.
!>   * `wm_packed` -- MINUS the energy-weighted density `X^AO`, packed triangle
!>     (the overlap-derivative kernel ADDS its argument, so the caller negates,
!>     exactly as `grd1::eijden` does for RHF and `casscf_gradient` does for CAS).
!>   * `lam`, `avec` -- the relaxed two-particle density in the separable form
!>
!>         Gamma^AO_{mu nu la si} = sum_k lam_k A^k_{mu nu} A^k_{la si},
!>
!>     with every `A^k` an ordinary symmetric AO matrix.  The caller obtains it
!>     from the eigendecomposition of the eight-fold-symmetrized `Gamma` over the
!>     composite index `(mu nu)`, so the representation is exact, not a fit, and
!>     each term is manifestly eight-fold symmetric.  Storage is
!>     `O(nvec * nbf^2)` instead of `nbf^4`, and every `A^k` takes the same
!>     `bfnrm` folding and spherical->Cartesian expansion the density takes.
!>
!> Unlike `casscf_gradient`, NO Hartree-Fock-like separable part is split off:
!> the CASPT2 two-particle density is not separable in any orbital block (the
!> first-order amplitudes correlate inactive, active and virtual space at once),
!> so the whole of it goes through the factorization.  `coulscale` and `hfscale`
!> are therefore zero and `get_density` evaluates the factorized sum alone.
!>
!> The symmetry petite list is deliberately NOT opted into, for the same reason
!> as in `casscf_gradient`: it is valid only for a totally symmetric contracted
!> density, and neither an excited CASPT2 root nor a relaxed response density
!> carries that guarantee.
module caspt2_gradient_mod

  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  use precision, only: dp
  use types, only: information
  use basis_tools, only: basis_set, bas_norm_matrix, build_cart_density
  use constants, only: HARMONIC_ACTIVE, NUM_CART_BF
  use grd2, only: grd2_driver, grd2_compute_data_t

  implicit none
  private

  integer, parameter :: i8 = c_int64_t

  public :: caspt2_gradient

  !> `info` slots handed back to the caller (all double):
  integer, parameter :: PT2_G_NVEC  = 0  !< factorization rank actually used
  integer, parameter :: PT2_G_TRACE = 1  !< trace of the AO density (electron count probe)
  integer, parameter :: PT2_G_NINFO = 2

  ! Status codes.  0 is success.
  integer(i8), parameter :: PT2G_ERR_DFT   = -21_i8  !< non-HF Hamiltonian
  integer(i8), parameter :: PT2G_ERR_ALLOC = -22_i8  !< out of memory
  integer(i8), parameter :: PT2G_ERR_SIZE  = -25_i8  !< nbf / nvec inconsistent

  !> Two-particle density for the derivative-ERI driver: the factorized
  !> `sum_k lam_k A^k (x) A^k` and nothing else.
  !>
  !> `lav` is `lam(k) * av(k,:,:)`, formed once so the inner loop does one
  !> multiply and one load less per contracted element.  It is built from the
  !> SCALED-and-EXPANDED `av_cart`, not scaled before the expansion, because
  !> floating-point multiplication does not associate.
  type, extends(grd2_compute_data_t) :: grd2_pt2_compute_data_t
    real(kind=dp), allocatable :: av(:,:,:)
    real(kind=dp), allocatable :: lav(:,:,:)
    real(kind=dp), allocatable :: lam(:)
    real(kind=dp), allocatable :: av_cart(:,:,:)
    real(kind=dp), allocatable :: lav_cart(:,:,:)
    integer, allocatable :: cart_off(:)
    integer :: nbf = 0
    integer :: nbf_cart = 0
    integer :: nvec = 0
  contains
    procedure :: init => grd2_pt2_init
    procedure :: clean => grd2_pt2_clean
    procedure :: get_density => grd2_pt2_get_density
    procedure :: build_cart => grd2_pt2_build_cart
  end type

contains

!###############################################################################

  !> C-bound entry point.  One call produces the complete nuclear gradient in
  !> `infos%atoms%grad`.
  !>
  !> @param[in]     c_handle   the OQP handle (`infos%atoms%grad` is overwritten)
  !> @param[in]     nbf        number of (spherical) basis functions
  !> @param[in]     dm_packed  packed AO one-particle density, nbf*(nbf+1)/2
  !> @param[in]     wm_packed  packed NEGATED energy-weighted density
  !> @param[in]     nvec       number of factorization vectors
  !> @param[in]     avec       `A^k`, C-order [nbf][nbf][nvec]
  !> @param[in]     lam        factorization eigenvalues, [nvec]
  !> @param[out]    info       PT2_G_NINFO diagnostics
  !> @return        0, or a negative status
  function caspt2_gradient_C(c_handle, nbf, dm_packed, wm_packed, nvec, avec, &
                             lam, info) result(status) bind(C, name="caspt2_gradient")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    type(oqp_handle_t) :: c_handle
    integer(c_int32_t), value :: nbf, nvec
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: dm_packed(0:*), wm_packed(0:*), avec(0:*), lam(0:*)
    real(dp), intent(inout) :: info(0:*)
    integer(i8) :: status
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    status = caspt2_gradient(inf, int(nbf), dm_packed, wm_packed, int(nvec), &
                             avec, lam, info)
  end function caspt2_gradient_C

!###############################################################################

  !> The run itself, reachable from an in-process Fortran caller.
  function caspt2_gradient(infos, nbf, dm_packed, wm_packed, nvec, avec, lam, &
                           info) result(status)
    use io_constants, only: iw
    use grd1, only: print_gradient
    use util, only: measure_time
    use printing, only: print_module_info

    type(information), target, intent(inout) :: infos
    integer, intent(in) :: nbf, nvec
    real(dp), intent(in) :: dm_packed(0:*), wm_packed(0:*), avec(0:*), lam(0:*)
    real(dp), intent(inout) :: info(0:*)
    integer(i8) :: status

    type(basis_set), pointer :: basis
    type(grd2_pt2_compute_data_t) :: gcomp
    real(dp), allocatable :: dm(:), wm(:)
    integer :: nbf_tri, k, p, q, ierr
    integer(i8) :: off
    real(dp) :: trace

    status = 0_i8
    do k = 0, PT2_G_NINFO - 1
      info(k) = 0.0_dp
    end do

    if (nbf <= 0 .or. nvec < 0) then
      status = PT2G_ERR_SIZE
      return
    end if
    ! The CASPT2 energy expression has no exchange-correlation term; the energy
    ! path already refuses a functional, this only pins it.
    if (infos%control%hamilton >= 20) then
      status = PT2G_ERR_DFT
      return
    end if

    nbf_tri = nbf * (nbf + 1) / 2
    allocate(dm(nbf_tri), wm(nbf_tri), stat=ierr)
    if (ierr /= 0) then
      status = PT2G_ERR_ALLOC
      return
    end if
    do k = 1, nbf_tri
      dm(k) = dm_packed(k - 1)
      wm(k) = wm_packed(k - 1)
    end do

    ! Electron-count probe: the packed triangle stores the strict upper
    ! triangle once, so the trace is the diagonal alone.
    trace = 0.0_dp
    do p = 1, nbf
      trace = trace + dm(p * (p + 1) / 2)
    end do
    info(PT2_G_TRACE) = trace
    info(PT2_G_NVEC) = real(nvec, dp)

    gcomp%nbf = nbf
    gcomp%nvec = nvec
    allocate(gcomp%lam(max(nvec, 1)), gcomp%av(max(nvec, 1), nbf, nbf), &
             gcomp%lav(max(nvec, 1), nbf, nbf), stat=ierr)
    if (ierr /= 0) then
      status = PT2G_ERR_ALLOC
      deallocate(dm, wm)
      return
    end if
    gcomp%lam = 0.0_dp
    gcomp%av = 0.0_dp
    gcomp%lav = 0.0_dp
    do k = 1, nvec
      gcomp%lam(k) = lam(k - 1)
    end do
    ! `avec` arrives as the C-order tensor [nbf][nbf][nvec], which is exactly
    ! the Fortran (nvec, nbf, nbf) layout: the vector index is contiguous, which
    ! is what the inner loop of `get_density` walks.
    do q = 1, nbf
      do p = 1, nbf
        do k = 1, nvec
          off = int(k - 1, i8) + int(nvec, i8) * (int(p - 1, i8) &
                + int(nbf, i8) * int(q - 1, i8))
          gcomp%av(k, p, q) = avec(off)
          gcomp%lav(k, p, q) = gcomp%lam(k) * avec(off)
        end do
      end do
    end do

    basis => infos%basis
    basis%atoms => infos%atoms

    open (unit=IW, file=infos%log_filename, position="append")
    call print_module_info('CASPT2_Gradient', &
                           'Analytic CASPT2 / XMS-CASPT2 Nuclear Gradient')
    write(iw,"(/' ..... Beginning CASPT2 Gradient Calculation...'/)")
    write(iw,"(' Relaxed 2-PDM factorization vectors ',I8)") nvec
    write(iw,"(' Relaxed AO density trace            ',ES16.8)") trace
    call flush(iw)

    call pt2_1e_grad(infos, basis, dm, wm)
    write(iw,"(' ..... End Of 1-Electron Gradient ......')")
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

    call pt2_2e_grad(infos, basis, gcomp, status)
    if (status < 0_i8) then
      close(iw)
      deallocate(dm, wm)
      call gcomp%clean()
      return
    end if
    write(iw, fmt="(' ...... End Of 2-Electron Gradient ......')")
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

    call print_gradient(infos)
    call measure_time(print_total=1, log_unit=iw)
    close(iw)

    deallocate(dm, wm)
    call gcomp%clean()
  end function caspt2_gradient

!###############################################################################

  !> Nuclear-repulsion, overlap (Pulay), kinetic, nuclear-attraction and ECP
  !> terms.  Structurally identical to `casscf_gradient_mod::cas_1e_grad`; only
  !> the two densities differ -- here they are the RELAXED CASPT2 ones.
  subroutine pt2_1e_grad(infos, basis, dm_packed, wm_packed)
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

      ! Energy-weighted density (already negated by the caller)
      call grad_ee_overlap(basis, wm_packed, grad, logtol=tol)

      call grad_ee_kinetic(basis, dm_packed, grad, logtol=tol)
      call grad_en_hellman_feynman(basis, xyz, zn, dm_packed, grad, logtol=tol)
      call grad_en_pulay(basis, xyz, zn, dm_packed, grad, logtol=tol)
      call grad_1e_ecp(infos, basis, xyz, dm_packed, grad, logtol=tol)

    end associate
  end subroutine pt2_1e_grad

!###############################################################################

  !> Two-electron term: contract the relaxed two-particle density against the
  !> derivative ERIs.
  subroutine pt2_2e_grad(infos, basis, gcomp, status)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    type(grd2_pt2_compute_data_t), intent(inout) :: gcomp
    integer(i8), intent(inout) :: status

    real(dp), allocatable :: de(:,:)
    integer :: ierr

    allocate(de(3, ubound(infos%atoms%zn, 1)), source=0.0_dp, stat=ierr)
    if (ierr /= 0) then
      status = PT2G_ERR_ALLOC
      return
    end if

    ! The entire two-particle density is carried by the factorization, so there
    ! is no Coulomb or exchange scale to apply.  `cam=.false.` explicitly: the
    ! range-separated path runs the driver twice with different scales and the
    ! factorized term is neither Coulomb nor exchange, so it would be added
    ! unscaled in both passes.
    gcomp%hfscale = 0.0_dp
    gcomp%coulscale = 0.0_dp
    call gcomp%init()
    call gcomp%build_cart(basis)
    call grd2_driver(infos, basis, de, gcomp, cam=.false.)
    infos%atoms%grad = infos%atoms%grad + de

    deallocate(de)
  end subroutine pt2_2e_grad

!###############################################################################

  subroutine grd2_pt2_init(this)
    class(grd2_pt2_compute_data_t), target, intent(inout) :: this
    this%nbf_cart = 0
  end subroutine grd2_pt2_init

  !> Fold the basis normalization into each factor and, under HARMONIC_ACTIVE,
  !> expand to the Cartesian-effective basis the derivative kernels iterate.
  !> Each `A^k` is an ordinary symmetric AO matrix, so it takes exactly the
  !> route the density takes in the Hartree-Fock and CASSCF paths.
  subroutine grd2_pt2_build_cart(this, basis)
    class(grd2_pt2_compute_data_t), intent(inout) :: this
    type(basis_set), intent(in) :: basis

    real(dp), allocatable :: tmp(:,:), expanded(:,:)
    integer, allocatable :: dummy_off(:)
    integer :: k, nv, nc_tmp

    if (.not. HARMONIC_ACTIVE) return

    ! `av` is always allocated with at least one (zero) vector, so the same
    ! path serves nvec = 0 and gives `get_density` a valid target either way.
    nv = max(this%nvec, 1)
    allocate(tmp(this%nbf, this%nbf))

    do k = 1, nv
      tmp = this%av(k, :, :)
      call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
      if (k == 1) then
        call build_cart_density(basis, tmp, expanded, this%cart_off, &
                                this%nbf_cart)
        allocate(this%av_cart(nv, this%nbf_cart, this%nbf_cart), &
                 this%lav_cart(nv, this%nbf_cart, this%nbf_cart))
      else
        call build_cart_density(basis, tmp, expanded, dummy_off, nc_tmp)
        if (allocated(dummy_off)) deallocate(dummy_off)
      end if
      this%av_cart(k, :, :) = expanded
      ! Scaled AFTER the norm folding and the Cartesian expansion, so the factor
      ! multiplied by `lam` is bit-for-bit the one the hot loop reads.
      this%lav_cart(k, :, :) = this%lam(k) * expanded
      deallocate(expanded)
    end do
    deallocate(tmp)
  end subroutine grd2_pt2_build_cart

  subroutine grd2_pt2_clean(this)
    class(grd2_pt2_compute_data_t), target, intent(inout) :: this
    if (allocated(this%av)) deallocate(this%av)
    if (allocated(this%lav)) deallocate(this%lav)
    if (allocated(this%lam)) deallocate(this%lam)
    if (allocated(this%av_cart)) deallocate(this%av_cart)
    if (allocated(this%lav_cart)) deallocate(this%lav_cart)
    if (allocated(this%cart_off)) deallocate(this%cart_off)
  end subroutine grd2_pt2_clean

!###############################################################################

  !> Two-particle density for one shell quartet.
  !>
  !> `grd2` contracts `dab` against the derivative ERI with the eight-fold
  !> multiplicity of the quartet, and the convention (see the Hartree-Fock and
  !> CASSCF paths) is `dab = 4 * Gamma^sym`.  With the whole density carried by
  !> the factorization that is `4 sum_k lam_k A^k_{ij} A^k_{kl}`.
  !>
  !> This is the one hot routine of the module: `grd2` calls it once per shell
  !> quartet surviving the coarse Schwarz screen, BEFORE the fine screen.  Both
  !> factors are stored with the vector index FIRST, so the inner loop walks two
  !> contiguous runs.
  subroutine grd2_pt2_get_density(this, basis, id, dab, dabmax)

    class(grd2_pt2_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: corr, df1
    logical :: usecart
    integer :: i, j, k, l, v, nv
    integer :: loc(4), nbf(4)
    integer :: i1, j1, k1, l1
    real(kind=dp), pointer :: ab(:,:,:,:)
    real(kind=dp), contiguous, pointer :: av(:,:,:)
    real(kind=dp), contiguous, pointer :: lavij(:), avkl(:)

    ! Local copy: `this` is polymorphic and has the TARGET attribute, so a
    ! component reference inside the loop nest cannot be assumed loop-invariant.
    nv = this%nvec

    usecart = HARMONIC_ACTIVE
    if (usecart) then
      loc = this%cart_off(id) - 1
      nbf = NUM_CART_BF(basis%am(id))
      av => this%av_cart
    else
      loc = basis%ao_offset(id) - 1
      nbf = basis%naos(id)
      av => this%av
    end if

    dabmax = 0
    ab(1:nbf(4),1:nbf(3),1:nbf(2),1:nbf(1)) => dab(1:product(nbf))

    do i = 1, nbf(1)
      i1 = loc(1) + i
      do j = 1, nbf(2)
        j1 = loc(2) + j
        if (usecart) then
          lavij => this%lav_cart(:,i1,j1)
        else
          lavij => this%lav(:,i1,j1)
        end if
        do k = 1, nbf(3)
          k1 = loc(3) + k
          do l = 1, nbf(4)
            l1 = loc(4) + l
            corr = 0.0_dp
            avkl => av(:,k1,l1)
            do v = 1, nv
              corr = corr + lavij(v)*avkl(v)
            end do
            df1 = 4*corr
            dabmax = max(dabmax, abs(df1))
            if (usecart) then
              ab(l,k,j,i) = df1
            else
              ab(l,k,j,i) = df1*product(basis%bfnrm([i1,j1,k1,l1]))
            end if
          end do
        end do
      end do
    end do
  end subroutine grd2_pt2_get_density

end module caspt2_gradient_mod
