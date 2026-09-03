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
    real(kind=dp), pointer :: mmat(:,:) => null()   !< probe M (nbf,nbf), full
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

  ! General-density derivative probe matching int2_mrsf_data_t_update exactly.
  ! The ordinary closed-shell mixed probe above is intentionally symmetric in
  ! its two density factors and cannot reconstruct the ordered Fock images of
  ! the nonsymmetric MRSF transition densities.
  type, extends(grd2_compute_data_t) :: grd2_mrsf_fockprobe_data_t
    real(kind=dp), pointer :: pmat(:,:) => null()
    real(kind=dp), pointer :: mmat(:,:) => null()
    real(kind=dp), allocatable :: pmat_cart(:,:),mmat_cart(:,:)
    integer, allocatable :: cart_off(:)
    integer :: nbf=0
  contains
    procedure :: init => grd2_mrsf_fockprobe_init
    procedure :: clean => grd2_mrsf_fockprobe_clean
    procedure :: get_density => grd2_mrsf_fockprobe_get_density
  end type grd2_mrsf_fockprobe_data_t

  private
  public :: grd2_fockprobe_data_t
  public :: grd2_fockprobe_os_data_t
  public :: fock_deriv_contract
  public :: fock_deriv_contract_scaled
  public :: fock_deriv_matrix
  public :: fock_deriv_matrix_general
  public :: fock_deriv_matrix_general_scaled
  public :: fock_deriv_matrix_mrsf_scaled
  public :: fock_deriv_matrix_mrsf_scaled_batch
  public :: fock_deriv_contract_os
  public :: fock_deriv_contract_os_batch
  public :: fock_deriv_matrix_os
  public :: fock_deriv_matrix_os_batch

contains

!###############################################################################

!> @brief Compute g_x = sum_uv M_uv F^x_uv[P] for every nuclear coordinate.
!> @param[in]  infos    system info (converged SCF)
!> @param[in]  basis    basis set
!> @param[in]  pmat     density P (nbf,nbf) full, AO basis (alpha density for RHF)
!> @param[in]  mmat     probe M (nbf,nbf) full, AO basis
!> @param[in]  hfscale  HF exchange scale (1.0 for HF; HFscale for hybrids)
!> @param[out] gx       (3, natom) contraction per nuclear coordinate
  subroutine fock_deriv_contract(infos, basis, pmat, mmat, hfscale, gx)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), target, intent(in) :: pmat(:,:), mmat(:,:)
    real(kind=dp), intent(in) :: hfscale
    real(kind=dp), intent(out) :: gx(:,:)

    call fock_deriv_contract_scaled(infos,basis,pmat,mmat,1.0_dp,hfscale,gx)
  end subroutine fock_deriv_contract

!###############################################################################

!> @brief General derivative response-Fock contraction with independent
!>        Coulomb and exchange scales.
!>
!> MRSF uses seven AO response densities. Channels 1:4 carry Coulomb and
!> exchange, whereas channels 5:7 carry exchange only; both are multiplied by
!> the response exact-exchange scale before the channel-specific SPC factors
!> are applied. Keeping the two scales explicit provides the derivative-ERI
!> action without changing representation or reconstructing four-index
!> derivative integrals in the MRSF driver.
  subroutine fock_deriv_contract_scaled(infos,basis,pmat,mmat,coulscale, &
      exchangescale,gx)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), target, intent(in) :: pmat(:,:),mmat(:,:)
    real(kind=dp), intent(in) :: coulscale,exchangescale
    real(kind=dp), intent(out) :: gx(:,:)

    type(grd2_fockprobe_data_t) :: gcomp
    integer, allocatable :: off_dummy(:)
    integer :: ncart

    gcomp%pmat=>pmat
    gcomp%mmat=>mmat
    gcomp%nbf=basis%nbf
    gcomp%coulscale=coulscale
    gcomp%hfscale=exchangescale
    gcomp%hfscale2=exchangescale
    if(HARMONIC_ACTIVE) then
      call fockprobe_cart(basis,pmat,gcomp%pmat_cart,gcomp%cart_off,ncart)
      call fockprobe_cart(basis,mmat,gcomp%mmat_cart,off_dummy,ncart)
    end if
    gx=0.0_dp
    call grd2_driver(infos,basis,gx,gcomp)
    call gcomp%clean()
  end subroutine fock_deriv_contract_scaled

!###############################################################################

!> @brief Build the complete symmetric two-electron derivative Fock matrix
!>        F^x[P] for every nuclear coordinate.
!>
!> @details The validated scalar contraction above is evaluated on the
!>          orthonormal basis of symmetric AO probe matrices.  A diagonal
!>          probe has M_uu=1, while an off-diagonal probe has
!>          M_uv=M_vu=1/2.  The closed-shell contraction is one half of the
!>          SCF response-Fock trace, so twice its value is the corresponding
!>          independent matrix element.  This
!>          reference implementation favors a direct, auditable relation to
!>          fock_deriv_contract; a later blocked-quartet implementation may
!>          replace it without changing the result or the calling convention.
!>
!> @param[in]  infos    system information for the converged SCF state
!> @param[in]  basis    basis set
!> @param[in]  pmat     closed-shell alpha density in the AO basis
!> @param[in]  hfscale  exact-exchange scale
!> @param[out] fmat     (nbf,nbf,3,natom) derivative two-electron Fock matrices
  subroutine fock_deriv_matrix(infos, basis, pmat, hfscale, fmat)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), target, intent(in) :: pmat(:,:)
    real(kind=dp), intent(in) :: hfscale
    real(kind=dp), intent(out) :: fmat(:,:,:,:)

    real(kind=dp), allocatable, target :: probe(:,:)
    real(kind=dp), allocatable :: gx(:,:), phalf(:,:)
    integer :: mu, nu, natom, nbf

    nbf = basis%nbf
    natom = size(basis%atoms%xyz, 2)
    if (any(shape(pmat) /= [nbf, nbf])) &
      error stop 'fock_deriv_matrix: density shape does not match the basis'
    if (any(shape(fmat) /= [nbf, nbf, 3, natom])) &
      error stop 'fock_deriv_matrix: output shape does not match the system'

    ! Direct single-traversal assembly.  The closed-shell response Fock
    !   F^x[P] = J^x[P] - (hfscale/2) K^x[P]
    ! is the open-shell derivative Fock J^x[pcoul] - hfscale K^x[pexch] with
    ! pcoul = P and pexch = P/2, so the exact AO-target driver builds every
    ! matrix element during one derivative-integral shell traversal.  The
    ! elementwise probe construction below costs nbf(nbf+1)/2 traversals and is
    ! kept only for the range-separated case, which the direct driver does not
    ! yet cover.
    if (.not. infos%dft%cam_flag) then
      allocate(phalf(nbf,nbf))
      phalf = 0.5_dp*pmat
      call fock_deriv_matrix_os(infos, basis, pmat, phalf, hfscale, fmat)
      deallocate(phalf)
      return
    end if

    allocate(probe(nbf,nbf), gx(3,natom))
    fmat = 0.0_dp
    do nu = 1, nbf
      do mu = nu, nbf
        probe = 0.0_dp
        if (mu == nu) then
          probe(mu,nu) = 1.0_dp
        else
          probe(mu,nu) = 0.5_dp
          probe(nu,mu) = 0.5_dp
        end if
        call fock_deriv_contract(infos, basis, pmat, probe, hfscale, gx)
        fmat(mu,nu,:,:) = 2.0_dp*gx
        fmat(nu,mu,:,:) = 2.0_dp*gx
      end do
    end do
  end subroutine fock_deriv_matrix

!###############################################################################

!> @brief Build F^x[P] for a general, possibly nonsymmetric AO density.
!>
!> @details Each ordered AO matrix unit is used as a probe.  Unlike
!>          fock_deriv_matrix, this routine does not identify transposed
!>          elements.  It is required for the separate P and P^T transition-
!>          density contractions in differentiated TDHF response operators.
!>
!> @param[in]  infos    system information for the converged SCF state
!> @param[in]  basis    basis set
!> @param[in]  pmat     general AO density
!> @param[in]  hfscale  exact-exchange scale
!> @param[out] fmat     (nbf,nbf,3,natom) derivative Fock matrices
  subroutine fock_deriv_matrix_general(infos, basis, pmat, hfscale, fmat)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), target, intent(in) :: pmat(:,:)
    real(kind=dp), intent(in) :: hfscale
    real(kind=dp), intent(out) :: fmat(:,:,:,:)

    call fock_deriv_matrix_general_scaled(infos,basis,pmat,1.0_dp,hfscale,fmat)
  end subroutine fock_deriv_matrix_general

!###############################################################################

!> @brief Build the ordered derivative response-Fock matrix with independent
!>        Coulomb and exchange scales for every nuclear coordinate.
  subroutine fock_deriv_matrix_general_scaled(infos,basis,pmat,coulscale, &
      exchangescale,fmat)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), target, intent(in) :: pmat(:,:)
    real(kind=dp), intent(in) :: coulscale,exchangescale
    real(kind=dp), intent(out) :: fmat(:,:,:,:)

    real(kind=dp), allocatable, target :: probe(:,:)
    real(kind=dp), allocatable :: gx(:,:)
    integer :: mu,nu,natom,nbf

    nbf=basis%nbf
    natom=size(basis%atoms%xyz,2)
    if(any(shape(pmat)/=[nbf,nbf])) &
      error stop 'fock_deriv_matrix_general_scaled: density shape mismatch'
    if(any(shape(fmat)/=[nbf,nbf,3,natom])) &
      error stop 'fock_deriv_matrix_general_scaled: output shape mismatch'
    allocate(probe(nbf,nbf),gx(3,natom))
    do nu=1,nbf
      do mu=1,nbf
        probe=0.0_dp
        probe(mu,nu)=1.0_dp
        call fock_deriv_contract_scaled(infos,basis,pmat,probe,coulscale, &
          exchangescale,gx)
        fmat(mu,nu,:,:)=2.0_dp*gx
      end do
    end do
    deallocate(probe,gx)
  end subroutine fock_deriv_matrix_general_scaled

!###############################################################################

!> Build the ordered derivative Fock matrix generated by the exact production
!> MRSF quartet update.  This differs from the closed-shell general probe above:
!> every ordered Coulomb and exchange target is retained for a nonsymmetric
!> transition density.
  subroutine fock_deriv_matrix_mrsf_scaled(infos,basis,pmat,coulscale, &
      exchangescale,fmat)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), target, intent(in) :: pmat(:,:)
    real(kind=dp), intent(in) :: coulscale,exchangescale
    real(kind=dp), intent(out) :: fmat(:,:,:,:)

    integer, parameter :: probe_block_size=64
    real(kind=dp), allocatable, target :: pmat_batch(:,:,:),probe_batch(:,:,:)
    real(kind=dp), allocatable :: gx_batch(:,:,:),coul_batch(:),exch_batch(:)
    integer, allocatable :: mu_index(:),nu_index(:)
    integer :: mu,nu,natom,nbf,nordered,probe,first,nactive,slot

    nbf=basis%nbf
    natom=size(basis%atoms%xyz,2)
    if(any(shape(pmat)/=[nbf,nbf])) &
      error stop 'fock_deriv_matrix_mrsf_scaled: density shape mismatch'
    if(any(shape(fmat)/=[nbf,nbf,3,natom])) &
      error stop 'fock_deriv_matrix_mrsf_scaled: output shape mismatch'
    nordered=nbf*nbf
    allocate(pmat_batch(nbf,nbf,probe_block_size), &
      probe_batch(nbf,nbf,probe_block_size), &
      gx_batch(3,natom,probe_block_size), &
      coul_batch(probe_block_size),exch_batch(probe_block_size), &
      mu_index(nordered),nu_index(nordered))
    fmat=0.0_dp
    probe=0
    do nu=1,nbf
      do mu=1,nbf
        probe=probe+1
        mu_index(probe)=mu
        nu_index(probe)=nu
      end do
    end do
    do first=1,nordered,probe_block_size
      nactive=min(probe_block_size,nordered-first+1)
      probe_batch=0.0_dp
      do slot=1,nactive
        mu=mu_index(first+slot-1)
        nu=nu_index(first+slot-1)
        pmat_batch(:,:,slot)=pmat
        probe_batch(mu,nu,slot)=1.0_dp
        coul_batch(slot)=coulscale
        exch_batch(slot)=exchangescale
      end do
      call fock_deriv_contract_mrsf_scaled_batch(infos,basis, &
        pmat_batch(:,:,1:nactive),probe_batch(:,:,1:nactive), &
        coul_batch(1:nactive),exch_batch(1:nactive), &
        gx_batch(:,:,1:nactive))
      do slot=1,nactive
        mu=mu_index(first+slot-1)
        nu=nu_index(first+slot-1)
        fmat(mu,nu,:,:)=gx_batch(:,:,slot)
      end do
    end do
    deallocate(pmat_batch,probe_batch,gx_batch,coul_batch,exch_batch, &
      mu_index,nu_index)
  end subroutine fock_deriv_matrix_mrsf_scaled

!###############################################################################

!> Exact direct derivative Fock matrices for ordered MRSF density channels.
!> This is the AO-target adjoint of the production MRSF quartet contraction;
!> no dense AO probe basis is formed.
  subroutine fock_deriv_matrix_mrsf_scaled_batch(infos,basis,pmat,coulscale, &
      exchangescale,fmat)
    use grd2, only: grd2_fock_deriv_mrsf_driver_batch
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: pmat(:,:,:),coulscale(:),exchangescale(:)
    real(kind=dp), intent(out) :: fmat(:,:,:,:,:)

    real(kind=dp), allocatable :: pmat_cart(:,:,:),pmat_one(:,:), &
      fcart(:,:,:,:,:)
    integer, allocatable :: cart_off(:),off_one(:)
    integer :: natom,nbf,ncart,ncart_one,nprobe,probe

    nbf=basis%nbf
    natom=size(basis%atoms%xyz,2)
    nprobe=size(pmat,3)
    if(nprobe<=0 .or. any(shape(pmat)/=[nbf,nbf,nprobe]) .or. &
       size(coulscale)/=nprobe .or. size(exchangescale)/=nprobe .or. &
       any(shape(fmat)/=[nbf,nbf,3,natom,nprobe])) error stop &
      'fock_deriv_matrix_mrsf_scaled_batch: inconsistent dimensions'
    fmat=0.0_dp
    if(.not.HARMONIC_ACTIVE) then
      ncart=nbf
      allocate(pmat_cart(ncart,ncart,nprobe),cart_off(basis%nshell))
      pmat_cart=pmat
      cart_off=basis%ao_offset
    else
      do probe=1,nprobe
        call fockprobe_cart(basis,pmat(:,:,probe),pmat_one,off_one,ncart_one)
        if(probe==1) then
          ncart=ncart_one
          allocate(pmat_cart(ncart,ncart,nprobe),cart_off(size(off_one)))
          cart_off=off_one
        else if(ncart_one/=ncart .or. any(off_one/=cart_off)) then
          error stop 'fock_deriv_matrix_mrsf_scaled_batch: layout changed'
        end if
        pmat_cart(:,:,probe)=pmat_one
        deallocate(pmat_one,off_one)
      end do
    end if
    allocate(fcart(ncart,ncart,3,natom,nprobe),source=0.0_dp)
    call grd2_fock_deriv_mrsf_driver_batch(infos,basis,pmat_cart,cart_off, &
      coulscale,exchangescale,fcart)
    if(HARMONIC_ACTIVE) then
      do probe=1,nprobe
        call reduce_cart_fock_deriv(basis,cart_off,fcart(:,:,:,:,probe), &
          fmat(:,:,:,:,probe))
      end do
    else
      fmat=fcart
    end if
    deallocate(pmat_cart,fcart,cart_off)
  end subroutine fock_deriv_matrix_mrsf_scaled_batch

!###############################################################################

!> Exact simultaneous MRSF derivative contractions.  The response channels
!> and ordered AO probes share each derivative-integral shell traversal while
!> retaining their independent Coulomb and exchange coefficients.
  subroutine fock_deriv_contract_mrsf_scaled_batch(infos,basis,pmat,mmat, &
      coulscale,exchangescale,gx)
    use grd2, only: grd2_driver_batch
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), target, intent(in) :: pmat(:,:,:),mmat(:,:,:)
    real(kind=dp), intent(in) :: coulscale(:),exchangescale(:)
    real(kind=dp), intent(out) :: gx(:,:,:)

    type(grd2_mrsf_fockprobe_data_t), allocatable :: gcomp(:)
    integer, allocatable :: off_dummy(:)
    integer :: nbf,natom,nprobe,probe,ncart

    nbf=basis%nbf
    natom=size(basis%atoms%xyz,2)
    nprobe=size(pmat,3)
    if(any(shape(pmat)/=[nbf,nbf,nprobe]) .or. &
       any(shape(mmat)/=[nbf,nbf,nprobe]) .or. &
       size(coulscale)/=nprobe .or. size(exchangescale)/=nprobe .or. &
       any(shape(gx)/=[3,natom,nprobe])) &
      error stop 'fock_deriv_contract_mrsf_scaled_batch: inconsistent dimensions'
    allocate(gcomp(nprobe))
    do probe=1,nprobe
      gcomp(probe)%pmat=>pmat(:,:,probe)
      gcomp(probe)%mmat=>mmat(:,:,probe)
      gcomp(probe)%nbf=nbf
      gcomp(probe)%coulscale=coulscale(probe)
      gcomp(probe)%hfscale=exchangescale(probe)
      gcomp(probe)%hfscale2=exchangescale(probe)
      if(HARMONIC_ACTIVE) then
        call fockprobe_cart(basis,pmat(:,:,probe), &
          gcomp(probe)%pmat_cart,gcomp(probe)%cart_off,ncart)
        call fockprobe_cart(basis,mmat(:,:,probe), &
          gcomp(probe)%mmat_cart,off_dummy,ncart)
      end if
    end do
    gx=0.0_dp
    call grd2_driver_batch(infos,basis,gx,gcomp,preserve_scales=.true.)
    do probe=1,nprobe
      call gcomp(probe)%clean()
    end do
    deallocate(gcomp)
  end subroutine fock_deriv_contract_mrsf_scaled_batch

!###############################################################################

  subroutine fock_deriv_contract_mrsf_scaled(infos,basis,pmat,mmat, &
      coulscale,exchangescale,gx)
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), target, intent(in) :: pmat(:,:),mmat(:,:)
    real(kind=dp), intent(in) :: coulscale,exchangescale
    real(kind=dp), intent(out) :: gx(:,:)

    type(grd2_mrsf_fockprobe_data_t) :: gcomp
    integer, allocatable :: off_dummy(:)
    integer :: ncart

    gcomp%pmat=>pmat
    gcomp%mmat=>mmat
    gcomp%nbf=basis%nbf
    gcomp%coulscale=coulscale
    gcomp%hfscale=exchangescale
    gcomp%hfscale2=exchangescale
    if(HARMONIC_ACTIVE) then
      call fockprobe_cart(basis,pmat,gcomp%pmat_cart,gcomp%cart_off,ncart)
      call fockprobe_cart(basis,mmat,gcomp%mmat_cart,off_dummy,ncart)
    end if
    gx=0.0_dp
    call grd2_driver(infos,basis,gx,gcomp,petite=.false.)
    call gcomp%clean()
  end subroutine fock_deriv_contract_mrsf_scaled

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
!>   identity tr(P . F^x[P]) = (2e gradient) holds exactly.  This bilinear
!>   expression retains the ordered M and P elements and therefore also
!>   applies to nonsymmetric transition densities.
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

  subroutine grd2_mrsf_fockprobe_init(this)
    class(grd2_mrsf_fockprobe_data_t), target, intent(inout) :: this
  end subroutine grd2_mrsf_fockprobe_init

!###############################################################################

  subroutine grd2_mrsf_fockprobe_clean(this)
    class(grd2_mrsf_fockprobe_data_t), target, intent(inout) :: this
    this%pmat=>null()
    this%mmat=>null()
    if(allocated(this%pmat_cart)) deallocate(this%pmat_cart)
    if(allocated(this%mmat_cart)) deallocate(this%mmat_cart)
    if(allocated(this%cart_off)) deallocate(this%cart_off)
  end subroutine grd2_mrsf_fockprobe_clean

!###############################################################################

!> Exact scalar contraction M:F[P] for one canonical ERI quartet, using the
!> same four Coulomb and eight exchange targets as int2_mrsf_data_t_update.
!> Repeated indices are deliberately not collapsed: the production update also
!> executes every target statement when two or more quartet indices coincide.
  subroutine grd2_mrsf_fockprobe_get_density(this,basis,id,dab,dabmax)
    class(grd2_mrsf_fockprobe_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: ccoef,xcoef,df1,bfn
    integer :: i,j,k,l,i1,j1,k1,l1
    integer :: loc(4),nbf(4)
    real(kind=dp), pointer :: ab(:,:,:,:),pmat(:,:),mmat(:,:)
    logical :: usecart

    ccoef=this%coulscale
    ! grd2_driver may replace hfscale for a DFT energy gradient.  The MRSF
    ! response scale is carried independently in hfscale2.
    xcoef=this%hfscale2
    usecart=HARMONIC_ACTIVE
    if(usecart) then
      pmat=>this%pmat_cart
      mmat=>this%mmat_cart
      loc=this%cart_off(id)-1
      nbf=NUM_CART_BF(basis%am(id))
    else
      pmat=>this%pmat
      mmat=>this%mmat
      loc=basis%ao_offset(id)-1
      nbf=basis%naos(id)
    end if

    dabmax=0.0_dp
    ab(1:nbf(4),1:nbf(3),1:nbf(2),1:nbf(1))=>dab(1:product(nbf))
    do i=1,nbf(1)
      i1=loc(1)+i
      do j=1,nbf(2)
        j1=loc(2)+j
        do k=1,nbf(3)
          k1=loc(3)+k
          do l=1,nbf(4)
            l1=loc(4)+l
            df1=ccoef*((mmat(i1,j1)+mmat(j1,i1))* &
                         (pmat(k1,l1)+pmat(l1,k1))+ &
                       (mmat(k1,l1)+mmat(l1,k1))* &
                         (pmat(i1,j1)+pmat(j1,i1)))
            if(xcoef/=0.0_dp) then
              df1=df1-xcoef*( &
                mmat(i1,k1)*pmat(j1,l1)+mmat(k1,i1)*pmat(l1,j1)+ &
                mmat(i1,l1)*pmat(j1,k1)+mmat(l1,i1)*pmat(k1,j1)+ &
                mmat(j1,k1)*pmat(i1,l1)+mmat(k1,j1)*pmat(l1,i1)+ &
                mmat(j1,l1)*pmat(i1,k1)+mmat(l1,j1)*pmat(k1,i1))
            end if
            bfn=1.0_dp
            if(.not.usecart) bfn=product(basis%bfnrm([i1,j1,k1,l1]))
            ab(l,k,j,i)=df1*bfn
            dabmax=max(dabmax,abs(ab(l,k,j,i)))
          end do
        end do
      end do
    end do
  end subroutine grd2_mrsf_fockprobe_get_density

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

!> Exact simultaneous open-shell derivative-Fock contractions.  All probes
!> share one derivative-integral shell traversal and one Rys recurrence per
!> shell quartet; the final density contractions remain independent.
  subroutine fock_deriv_contract_os_batch(infos,basis,pcoul,pexch,mmat, &
      hfscale,gx)
    use grd2, only: grd2_driver_batch
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), target, intent(in) :: pcoul(:,:,:),pexch(:,:,:),mmat(:,:,:)
    real(kind=dp), intent(in) :: hfscale
    real(kind=dp), intent(out) :: gx(:,:,:)

    type(grd2_fockprobe_os_data_t), allocatable :: gcomp(:)
    integer, allocatable :: off_dummy(:)
    integer :: nbf,natom,nprobe,probe,ncart

    nbf=basis%nbf
    natom=size(basis%atoms%xyz,2)
    nprobe=size(mmat,3)
    if(any(shape(pcoul)/=[nbf,nbf,nprobe]) .or. &
       any(shape(pexch)/=[nbf,nbf,nprobe]) .or. &
       any(shape(gx)/=[3,natom,nprobe])) &
      error stop 'fock_deriv_contract_os_batch: inconsistent batch dimensions'
    allocate(gcomp(nprobe))
    do probe=1,nprobe
      gcomp(probe)%pcoul=>pcoul(:,:,probe)
      gcomp(probe)%pexch=>pexch(:,:,probe)
      gcomp(probe)%mmat=>mmat(:,:,probe)
      gcomp(probe)%nbf=nbf
      gcomp(probe)%coulscale=1.0_dp
      gcomp(probe)%hfscale=hfscale
      gcomp(probe)%hfscale2=hfscale
      if(HARMONIC_ACTIVE) then
        call fockprobe_cart(basis,pcoul(:,:,probe), &
          gcomp(probe)%pcoul_cart,gcomp(probe)%cart_off,ncart)
        call fockprobe_cart(basis,pexch(:,:,probe), &
          gcomp(probe)%pexch_cart,off_dummy,ncart)
        call fockprobe_cart(basis,mmat(:,:,probe), &
          gcomp(probe)%mmat_cart,off_dummy,ncart)
      end if
    end do
    gx=0.0_dp
    call grd2_driver_batch(infos,basis,gx,gcomp)
    do probe=1,nprobe
      call gcomp(probe)%clean()
    end do
    deallocate(gcomp)
  end subroutine fock_deriv_contract_os_batch

!###############################################################################

!> @brief Build the full symmetric explicit derivative of an open-shell spin
!>        Fock matrix J^x[P_alpha+P_beta]-c_x K^x[P_spin].
!>
!> @details fock_deriv_contract_os returns the genuine matrix trace against a
!>          symmetric probe.  A unit diagonal probe therefore yields F_uu,
!>          while a half-weighted symmetric off-diagonal probe yields F_uv for
!>          a symmetric derivative Fock.  This validation-grade implementation
!>          avoids reconstructing derivative four-index integrals and is used
!>          by the first MRSF Hessian driver for modest molecular systems.
  subroutine fock_deriv_matrix_os(infos,basis,pcoul,pexch,hfscale,fmat)
    use grd2, only: grd2_fock_deriv_os_driver
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), target, intent(in) :: pcoul(:,:),pexch(:,:)
    real(kind=dp), intent(in) :: hfscale
    real(kind=dp), intent(out) :: fmat(:,:,:,:)

    integer, parameter :: probe_block_size=64
    real(kind=dp), allocatable, target :: pcoul_batch(:,:,:), &
      pexch_batch(:,:,:),probe_batch(:,:,:)
    real(kind=dp), allocatable :: gx_batch(:,:,:),pcoul_cart(:,:), &
      pexch_cart(:,:),fcart(:,:,:,:)
    integer, allocatable :: mu_index(:),nu_index(:),cart_off(:),off_dummy(:)
    integer :: mu,nu,natom,nbf,nunique,probe,first,nactive,slot,ncart

    nbf=basis%nbf
    natom=size(basis%atoms%xyz,2)
    if(any(shape(pcoul)/=[nbf,nbf]) .or. any(shape(pexch)/=[nbf,nbf])) &
      error stop 'fock_deriv_matrix_os: density shape does not match the basis'
    if(any(shape(fmat)/=[nbf,nbf,3,natom])) &
      error stop 'fock_deriv_matrix_os: output shape does not match the system'

    if(HARMONIC_ACTIVE) then
      call fockprobe_cart(basis,pcoul,pcoul_cart,cart_off,ncart)
      call fockprobe_cart(basis,pexch,pexch_cart,off_dummy,ncart)
      allocate(fcart(ncart,ncart,3,natom),source=0.0_dp)
      call grd2_fock_deriv_os_driver(infos,basis,pcoul_cart,pexch_cart, &
        cart_off,hfscale,fcart)
      call reduce_cart_fock_deriv(basis,cart_off,fcart,fmat)
      deallocate(pcoul_cart,pexch_cart,fcart,cart_off,off_dummy)
      return
    end if

    nunique=nbf*(nbf+1)/2
    allocate(pcoul_batch(nbf,nbf,probe_block_size), &
      pexch_batch(nbf,nbf,probe_block_size), &
      probe_batch(nbf,nbf,probe_block_size), &
      gx_batch(3,natom,probe_block_size),mu_index(nunique),nu_index(nunique))
    fmat=0.0_dp
    probe=0
    do nu=1,nbf
      do mu=nu,nbf
        probe=probe+1
        mu_index(probe)=mu
        nu_index(probe)=nu
      end do
    end do
    do first=1,nunique,probe_block_size
      nactive=min(probe_block_size,nunique-first+1)
      probe_batch=0.0_dp
      do slot=1,nactive
        mu=mu_index(first+slot-1)
        nu=nu_index(first+slot-1)
        pcoul_batch(:,:,slot)=pcoul
        pexch_batch(:,:,slot)=pexch
        if(mu==nu) then
          probe_batch(mu,nu,slot)=1.0_dp
        else
          probe_batch(mu,nu,slot)=0.5_dp
          probe_batch(nu,mu,slot)=0.5_dp
        end if
      end do
      call fock_deriv_contract_os_batch(infos,basis, &
        pcoul_batch(:,:,1:nactive),pexch_batch(:,:,1:nactive), &
        probe_batch(:,:,1:nactive),hfscale,gx_batch(:,:,1:nactive))
      do slot=1,nactive
        mu=mu_index(first+slot-1)
        nu=nu_index(first+slot-1)
        fmat(mu,nu,:,:)=gx_batch(:,:,slot)
        fmat(nu,mu,:,:)=gx_batch(:,:,slot)
      end do
    end do
    deallocate(pcoul_batch,pexch_batch,probe_batch,gx_batch,mu_index,nu_index)
  end subroutine fock_deriv_matrix_os

!###############################################################################

!> Build exact explicit derivative Fock matrices for several open-shell
!> density pairs in one derivative-ERI traversal.  The single-probe fallback
!> is retained for non-harmonic builds; production harmonic builds share all
!> Rys recurrence and derivative-integral work without an integral model.
  subroutine fock_deriv_matrix_os_batch(infos,basis,pcoul,pexch,hfscale,fmat)
    use grd2, only: grd2_fock_deriv_os_driver_batch
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: pcoul(:,:,:),pexch(:,:,:)
    real(kind=dp), intent(in) :: hfscale
    real(kind=dp), intent(out) :: fmat(:,:,:,:,:)

    real(kind=dp), allocatable :: pcoul_cart(:,:,:),pexch_cart(:,:,:), &
      pcoul_one(:,:),pexch_one(:,:),fcart(:,:,:,:,:)
    integer, allocatable :: cart_off(:),off_one(:),off_exchange(:)
    integer :: natom,nbf,ncart,ncart_one,nprobe,probe

    nbf=basis%nbf
    natom=size(basis%atoms%xyz,2)
    nprobe=size(pcoul,3)
    if(nprobe<=0 .or. any(shape(pcoul)/=[nbf,nbf,nprobe]) .or. &
       any(shape(pexch)/=[nbf,nbf,nprobe]) .or. &
       any(shape(fmat)/=[nbf,nbf,3,natom,nprobe])) &
      error stop 'fock_deriv_matrix_os_batch: inconsistent dimensions'
    fmat=0.0_dp
    if(.not.HARMONIC_ACTIVE) then
      do probe=1,nprobe
        call fock_deriv_matrix_os(infos,basis,pcoul(:,:,probe), &
          pexch(:,:,probe),hfscale,fmat(:,:,:,:,probe))
      end do
      return
    end if

    do probe=1,nprobe
      call fockprobe_cart(basis,pcoul(:,:,probe),pcoul_one,off_one,ncart_one)
      call fockprobe_cart(basis,pexch(:,:,probe),pexch_one,off_exchange, &
        ncart_one)
      if(any(off_exchange/=off_one)) &
        error stop 'fock_deriv_matrix_os_batch: Cartesian layout changed'
      if(probe==1) then
        ncart=ncart_one
        allocate(pcoul_cart(ncart,ncart,nprobe), &
          pexch_cart(ncart,ncart,nprobe),cart_off(size(off_one)))
        cart_off=off_one
      else if(ncart_one/=ncart .or. any(off_one/=cart_off) .or. &
              any(off_exchange/=cart_off)) then
        error stop 'fock_deriv_matrix_os_batch: Cartesian layout changed'
      end if
      pcoul_cart(:,:,probe)=pcoul_one
      pexch_cart(:,:,probe)=pexch_one
      deallocate(pcoul_one,pexch_one,off_one,off_exchange)
    end do
    allocate(fcart(ncart,ncart,3,natom,nprobe),source=0.0_dp)
    call grd2_fock_deriv_os_driver_batch(infos,basis,pcoul_cart,pexch_cart, &
      cart_off,hfscale,fcart)
    do probe=1,nprobe
      call reduce_cart_fock_deriv(basis,cart_off,fcart(:,:,:,:,probe), &
        fmat(:,:,:,:,probe))
    end do
    deallocate(pcoul_cart,pexch_cart,fcart,cart_off)
  end subroutine fock_deriv_matrix_os_batch

!###############################################################################

!> Transform a Cartesian derivative Fock matrix back to the native AO basis.
!> This is the covariant counterpart of fockprobe_cart and is applied to every
!> nuclear coordinate after the direct derivative-ERI build.
  subroutine reduce_cart_fock_deriv(basis,cart_off,cart,sph)
    use cart2sph, only: cart2sph_mat
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: cart_off(:)
    real(kind=dp), intent(in) :: cart(:,:,:,:)
    real(kind=dp), intent(out) :: sph(:,:,:,:)

    real(kind=dp), allocatable :: blk(:)
    integer :: atom,axis,si,sj,ii,jj,nci,ncj,nsi,nsj
    integer :: coi,coj,soi,soj,idx

    sph=0.0_dp
    do atom=1,size(sph,4)
      do axis=1,3
        do sj=1,basis%nshell
          ncj=NUM_CART_BF(basis%am(sj))
          nsj=basis%naos(sj)
          coj=cart_off(sj)
          soj=basis%ao_offset(sj)
          do si=1,basis%nshell
            nci=NUM_CART_BF(basis%am(si))
            nsi=basis%naos(si)
            coi=cart_off(si)
            soi=basis%ao_offset(si)
            allocate(blk(nci*ncj))
            idx=0
            do jj=1,ncj
              do ii=1,nci
                idx=idx+1
                blk(idx)=cart(coi+ii-1,coj+jj-1,axis,atom)
              end do
            end do
            if(basis%harmonic(si)==1 .or. basis%harmonic(sj)==1) &
              call cart2sph_mat(blk,basis%am(si),basis%harmonic(si), &
                basis%am(sj),basis%harmonic(sj))
            idx=0
            do jj=1,nsj
              do ii=1,nsi
                idx=idx+1
                sph(soi+ii-1,soj+jj-1,axis,atom)=blk(idx)
              end do
            end do
            deallocate(blk)
          end do
        end do
        call bas_norm_matrix(sph(:,:,axis,atom),basis%bfnrm,basis%nbf)
      end do
    end do
  end subroutine reduce_cart_fock_deriv

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
