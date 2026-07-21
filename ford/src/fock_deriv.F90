module fock_deriv_mod
!> @brief Two-electron derivative-Fock contraction  tr(M . F^x[P])  for the
!>   CPHF nuclear right-hand side, built natively on top of the validated 2e
!>   gradient driver (grd2_driver) with no changes to the Rys internals and no
!>   libint.
!>
!>   The 2e gradient driver contracts the derivative ERIs d(uv|ls)/dx with a
!>   four-index density product supplied by a grd2_compute_data_t extension. The
!>   standard (energy-gradient) extension forms D (x) D. Here we instead form a
!>   MIXED product M (x) P, so the same driver returns, for each nuclear
!>   coordinate x,
!>     g_x = sum_{uvls} d(uv|ls)/dx * [ 4 c M_uv P_ls - x_hf ( M_ul P_vs + M_us P_vl ) ]
!>   which is exactly  sum_uv M_uv F^x_uv[P]  for the closed-shell response Fock
!>   F^x[P] = J^x[P] - 1/2 K^x[P] (Coulomb scaled by c, exchange by x_hf=HFscale),
!>   summed over the two equivalent index orderings that the driver already
!>   exploits. M is the "probe" matrix; for a CPHF RHS element B^x_{ia} the probe
!>   is the symmetric AO matrix C_{.,i} C_{.,a}^T + C_{.,a} C_{.,i}^T.
!>
!>   This is the F^x building block of the native CPHF chain. It is validated by
!>   the trace identity tr(P . F^x[P]) = (2e part of dE/dx), i.e. against the
!>   already-validated grd2_driver energy gradient (exact, non-iterative).

  use precision, only: dp
  use grd2, only: grd2_driver, grd2_compute_data_t
  use basis_tools, only: basis_set, bas_norm_matrix, build_cart_density
  use constants, only: HARMONIC_ACTIVE, NUM_CART_BF
  use types, only: information

  implicit none

  character(len=*), parameter :: module_name = "fock_deriv_mod"

  !> grd2 compute-data extension forming the mixed two-density product M (x) P.
  type, extends(grd2_compute_data_t) :: grd2_fockprobe_data_t
    real(kind=dp), pointer :: pmat(:,:) => null()   !< density P (nbf,nbf), full
    real(kind=dp), pointer :: mmat(:,:) => null()   !< probe  M (nbf,nbf), full (symmetric)
    ! Cartesian-effective (bfnrm-folded) copies + Cartesian offsets, used under
    ! HARMONIC_ACTIVE so the spherical probe/density contract with Cartesian
    ! derivative ERIs (set by prepare_cart).
    real(kind=dp), allocatable :: pmat_cart(:,:), mmat_cart(:,:)
    integer, allocatable :: cart_off(:)
    integer :: nbf = 0
  contains
    procedure :: init => grd2_fockprobe_init
    procedure :: clean => grd2_fockprobe_clean
    procedure :: get_density => grd2_fockprobe_get_density
  end type

  !> Open-shell (UHF/ROHF) extension: the probe M is contracted against a
  !> Coulomb density (the spin-summed total) and a SEPARATE exchange density
  !> (one spin), so the contraction returns the genuine open-shell derivative-
  !> Fock trace
  !>    g_x = sum_uv M_uv ( J^x_uv[pcoul] - c_x K^x_uv[pexch] )
  !> with the full (not 1/2) open-shell exchange factor, matching the spin-s
  !> Fock that scf_addons::fock_jk assembles for scftype>=2.  Setting
  !> pcoul == pexch == P does NOT reduce to the closed-shell grd2_fockprobe_data_t
  !> (which carries the closed-shell 1/2 K factor); this object is for the
  !> open-shell response only.
  type, extends(grd2_compute_data_t) :: grd2_fockprobe_os_data_t
    real(kind=dp), pointer :: pcoul(:,:) => null()  !< Coulomb density (total = Pa+Pb)
    real(kind=dp), pointer :: pexch(:,:) => null()  !< exchange density (one spin)
    real(kind=dp), pointer :: mmat(:,:) => null()   !< probe M (nbf,nbf), full (symmetric)
    real(kind=dp), allocatable :: pcoul_cart(:,:), pexch_cart(:,:), mmat_cart(:,:)
    integer, allocatable :: cart_off(:)
    integer :: nbf = 0
  contains
    procedure :: init => grd2_fockprobe_os_init
    procedure :: clean => grd2_fockprobe_os_clean
    procedure :: get_density => grd2_fockprobe_os_get_density
  end type

  private
  public :: grd2_fockprobe_data_t
  public :: grd2_fockprobe_os_data_t
  public :: fock_deriv_contract
  public :: fock_deriv_contract_os

contains

!###############################################################################

!> @brief Compute g_x = sum_uv M_uv F^x_uv[P] for every nuclear coordinate.
!> @param[in]  infos    system info (converged SCF)
!> @param[in]  basis    basis set
!> @param[in]  pmat     density P (nbf,nbf) full, AO basis (alpha density for RHF)
!> @param[in]  mmat     probe M (nbf,nbf) full, symmetric, AO basis
!> @param[in]  hfscale  HF exchange scale (1.0 for HF; HFscale for hybrids)
!> @param[out] gx       (3, natom) contraction per nuclear coordinate
  subroutine fock_deriv_contract(infos, basis, pmat, mmat, hfscale, gx)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), target, intent(in) :: pmat(:,:), mmat(:,:)
    real(kind=dp), intent(in) :: hfscale
    real(kind=dp), intent(out) :: gx(:,:)

    type(grd2_fockprobe_data_t) :: gcomp

    integer, allocatable :: off_dummy(:)
    integer :: ncart

    gcomp%pmat => pmat
    gcomp%mmat => mmat
    gcomp%nbf = basis%nbf
    gcomp%coulscale = 1.0_dp
    gcomp%hfscale = hfscale
    gcomp%hfscale2 = hfscale

    if (HARMONIC_ACTIVE) then
      call fockprobe_cart(basis, pmat, gcomp%pmat_cart, gcomp%cart_off, ncart)
      call fockprobe_cart(basis, mmat, gcomp%mmat_cart, off_dummy, ncart)
    end if

    gx = 0.0_dp
    call grd2_driver(infos, basis, gx, gcomp)
  end subroutine fock_deriv_contract

!###############################################################################

!> @brief Cartesian-effective (bfnrm-folded) copy of a full AO matrix +
!>   Cartesian per-shell offsets, for contraction with Cartesian derivative
!>   ERIs under HARMONIC_ACTIVE.
  subroutine fockprobe_cart(basis, m, m_cart, cart_off, nbf_cart)
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: m(:,:)
    real(kind=dp), allocatable, intent(out) :: m_cart(:,:)
    integer, allocatable, intent(out) :: cart_off(:)
    integer, intent(out) :: nbf_cart
    real(kind=dp), allocatable :: tmp(:,:)
    tmp = m
    call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
    call build_cart_density(basis, tmp, m_cart, cart_off, nbf_cart)
  end subroutine fockprobe_cart

!###############################################################################

  subroutine grd2_fockprobe_init(this)
    class(grd2_fockprobe_data_t), target, intent(inout) :: this
    ! pmat/mmat are full matrices supplied by the caller; nothing to unpack.
  end subroutine grd2_fockprobe_init

!###############################################################################

  subroutine grd2_fockprobe_clean(this)
    class(grd2_fockprobe_data_t), target, intent(inout) :: this
    this%pmat => null()
    this%mmat => null()
  end subroutine grd2_fockprobe_clean

!###############################################################################

!> @brief Mixed two-density product for the shell quartet, matching the layout
!>   and normalization of grd2_rhf_compute_data_t_get_density but replacing the
!>   second density factor with the probe M.  The energy-gradient routine forms
!>     4 c D_ij D_kl - x_hf ( D_ik D_jl + D_il D_jk ).
!>   To contract d(uv|ls)/dx with M on the (i,j)=(u,v) pair and P on the
!>   (k,l)=(l,s) pair (Coulomb), and M/P spread across the exchange index
!>   pairings symmetrically, we form
!>     4 c M_ij P_kl
!>     - x_hf/2 ( M_ik P_jl + M_il P_jk + P_ik M_jl + P_il M_jk ).
!>   The 1/2 with the four symmetric exchange terms reproduces the same total as
!>   the energy routine's x_hf ( D_ik D_jl + D_il D_jk ) when M = P, so the trace
!>   identity tr(P . F^x[P]) = (2e gradient) holds exactly.
  subroutine grd2_fockprobe_get_density(this, basis, id, dab, dabmax)
    class(grd2_fockprobe_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: coulfact, xcfact, df1, dq1, bfn
    integer :: i, j, k, l, i1, j1, k1, l1
    integer :: loc(4), nbf(4)
    real(kind=dp), pointer :: ab(:,:,:,:)
    real(kind=dp), pointer :: pmat(:,:), mmat(:,:)
    logical :: usecart

    coulfact = 4*this%coulscale
    xcfact = this%hfscale

    usecart = HARMONIC_ACTIVE
    if (usecart) then
      pmat => this%pmat_cart;  mmat => this%mmat_cart
      loc = this%cart_off(id) - 1
      nbf = NUM_CART_BF(basis%am(id))
    else
      pmat => this%pmat;  mmat => this%mmat
      loc = basis%ao_offset(id) - 1
      nbf = basis%naos(id)
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
            ! Coulomb: symmetrized over the (ij)<->(kl) permutation the driver
            ! exploits: 2 c (M_ij P_kl + M_kl P_ij) (reduces to 4 c P_ij P_kl at M=P).
            df1 = 0.5_dp*coulfact*( mmat(i1,j1)*pmat(k1,l1) &
                                  + mmat(k1,l1)*pmat(i1,j1) )
            if (xcfact/=0.0_dp) then
              ! Exchange: symmetrized 4-term (reduces to x(P_ik P_jl+P_il P_jk) at M=P).
              dq1 = 0.5_dp*( mmat(i1,k1)*pmat(j1,l1) &
                          + mmat(i1,l1)*pmat(j1,k1) &
                          + pmat(i1,k1)*mmat(j1,l1) &
                          + pmat(i1,l1)*mmat(j1,k1) )
              df1 = df1 - xcfact*dq1
            end if
            dabmax = max(dabmax, abs(df1))
            bfn = 1.0_dp
            if (.not. usecart) bfn = product(basis%bfnrm([i1,j1,k1,l1]))
            ab(l,k,j,i) = df1*bfn
          end do
        end do
      end do
    end do
  end subroutine grd2_fockprobe_get_density

!###############################################################################
!  Open-shell (UHF/ROHF) two-density derivative-Fock contraction
!###############################################################################

!> @brief Compute g_x = sum_uv M_uv ( J^x_uv[pcoul] - c_x K^x_uv[pexch] ) for
!>   every nuclear coordinate, i.e. the trace of the probe M against the
!>   open-shell spin-s derivative Fock with Coulomb from the total density and
!>   exchange from the spin density.
!> @param[in]  infos    system info (converged SCF)
!> @param[in]  basis    basis set
!> @param[in]  pcoul    Coulomb density (total Pa+Pb), full (nbf,nbf), AO basis
!> @param[in]  pexch    exchange density (one spin), full (nbf,nbf), AO basis
!> @param[in]  mmat     probe M (nbf,nbf) full, symmetric, AO basis
!> @param[in]  hfscale  HF exchange scale (1.0 for HF; HFscale for hybrids)
!> @param[out] gx       (3, natom) contraction per nuclear coordinate
  subroutine fock_deriv_contract_os(infos, basis, pcoul, pexch, mmat, hfscale, gx)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), target, intent(in) :: pcoul(:,:), pexch(:,:), mmat(:,:)
    real(kind=dp), intent(in) :: hfscale
    real(kind=dp), intent(out) :: gx(:,:)

    type(grd2_fockprobe_os_data_t) :: gcomp

    integer, allocatable :: off_dummy(:)
    integer :: ncart

    gcomp%pcoul => pcoul
    gcomp%pexch => pexch
    gcomp%mmat => mmat
    gcomp%nbf = basis%nbf
    gcomp%coulscale = 1.0_dp
    gcomp%hfscale = hfscale
    gcomp%hfscale2 = hfscale

    if (HARMONIC_ACTIVE) then
      call fockprobe_cart(basis, pcoul, gcomp%pcoul_cart, gcomp%cart_off, ncart)
      call fockprobe_cart(basis, pexch, gcomp%pexch_cart, off_dummy, ncart)
      call fockprobe_cart(basis, mmat,  gcomp%mmat_cart,  off_dummy, ncart)
    end if

    gx = 0.0_dp
    call grd2_driver(infos, basis, gx, gcomp)
  end subroutine fock_deriv_contract_os

!###############################################################################

  subroutine grd2_fockprobe_os_init(this)
    class(grd2_fockprobe_os_data_t), target, intent(inout) :: this
    ! densities/probe are full matrices supplied by the caller; nothing to do.
  end subroutine grd2_fockprobe_os_init

!###############################################################################

  subroutine grd2_fockprobe_os_clean(this)
    class(grd2_fockprobe_os_data_t), target, intent(inout) :: this
    this%pcoul => null()
    this%pexch => null()
    this%mmat => null()
  end subroutine grd2_fockprobe_os_clean

!###############################################################################

!> @brief Open-shell mixed-density product for the shell quartet.
!>   The closed-shell grd2_fockprobe_get_density forms
!>     2 c (M_ij P_kl + M_kl P_ij) - x_hf/2 (M_ik P_jl + M_il P_jk
!>                                            + P_ik M_jl + P_il M_jk),
!>   which the validated builder maps to  1/2 Tr[M J^x[P]] - 1/4 c_x Tr[M K^x[P]]
!>   (the closed-shell Fock J - 1/2 K).  Here we want the FULL open-shell trace
!>     Tr[M J^x[pcoul]] - c_x Tr[M K^x[pexch]],
!>   i.e. twice the Coulomb coefficient and four times the exchange coefficient,
!>   with Coulomb taking the total density (pcoul) and exchange the spin density
!>   (pexch).  Validated by fock_deriv_os_selftest against a finite difference of
!>   Tr[M . fock_jk_spin] at frozen densities.
  subroutine grd2_fockprobe_os_get_density(this, basis, id, dab, dabmax)
    class(grd2_fockprobe_os_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: ccoef, xcoef, df1, dq1, bfn
    integer :: i, j, k, l, i1, j1, k1, l1
    integer :: loc(4), nbf(4)
    real(kind=dp), pointer :: ab(:,:,:,:)
    real(kind=dp), pointer :: mmat(:,:), pcoul(:,:), pexch(:,:)
    logical :: usecart

    ccoef = 4*this%coulscale     ! 2x the closed-shell Coulomb -> full Tr[M J^x[pcoul]]
    xcoef = 2*this%hfscale       ! 4x the closed-shell exchange -> full c_x Tr[M K^x[pexch]]

    usecart = HARMONIC_ACTIVE
    if (usecart) then
      mmat => this%mmat_cart;  pcoul => this%pcoul_cart;  pexch => this%pexch_cart
      loc = this%cart_off(id) - 1
      nbf = NUM_CART_BF(basis%am(id))
    else
      mmat => this%mmat;  pcoul => this%pcoul;  pexch => this%pexch
      loc = basis%ao_offset(id) - 1
      nbf = basis%naos(id)
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
            ! Coulomb: M against the TOTAL density (symmetrized over (ij)<->(kl)).
            df1 = ccoef*( mmat(i1,j1)*pcoul(k1,l1) &
                        + mmat(k1,l1)*pcoul(i1,j1) )
            if (xcoef/=0.0_dp) then
              ! Exchange: M against the SPIN density (4-term symmetrized).
              dq1 = mmat(i1,k1)*pexch(j1,l1) &
                  + mmat(i1,l1)*pexch(j1,k1) &
                  + pexch(i1,k1)*mmat(j1,l1) &
                  + pexch(i1,l1)*mmat(j1,k1)
              df1 = df1 - xcoef*dq1
            end if
            dabmax = max(dabmax, abs(df1))
            bfn = 1.0_dp
            if (.not. usecart) bfn = product(basis%bfnrm([i1,j1,k1,l1]))
            ab(l,k,j,i) = df1*bfn
          end do
        end do
      end do
    end do
  end subroutine grd2_fockprobe_os_get_density

end module fock_deriv_mod
