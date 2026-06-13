module hf_gradient_mod

  use precision, only: dp
  use grd2, only: grd2_driver, grd2_compute_data_t
  use basis_tools, only: basis_set, bas_norm_matrix, build_cart_density
  use constants, only: HARMONIC_ACTIVE, NUM_CART_BF
  use types, only: information

  implicit none

  character(len=*), parameter :: module_name = "hf_gradient_mod"

!###############################################################################

  type, extends(grd2_compute_data_t), abstract :: grd2_hf_compute_data_t
    real(kind=dp), pointer :: da(:) => null()
    real(kind=dp), pointer :: db(:) => null()
    real(kind=dp), allocatable :: d2a(:,:), d2b(:,:)
    ! Cartesian-effective (bfnrm-folded) densities + Cartesian offsets, built
    ! by build_cart and used by get_density under HARMONIC_ACTIVE.
    real(kind=dp), allocatable :: d2a_cart(:,:), d2b_cart(:,:)
    integer, allocatable :: cart_off(:)
    integer :: nbf_cart = 0
    integer :: nbf = 0
  contains
    procedure :: build_cart => grd2_hf_build_cart
  end type

!###############################################################################

  type, extends(grd2_hf_compute_data_t) :: grd2_rhf_compute_data_t
  contains
    procedure :: init => grd2_rhf_compute_data_t_init
    procedure :: clean => grd2_rhf_compute_data_t_clean
    procedure :: get_density => grd2_rhf_compute_data_t_get_density
  end type

!###############################################################################

  type, extends(grd2_hf_compute_data_t) :: grd2_uhf_compute_data_t
  contains
    procedure :: init => grd2_uhf_compute_data_t_init
    procedure :: clean => grd2_uhf_compute_data_t_clean
    procedure :: get_density => grd2_uhf_compute_data_t_get_density
  end type

!###############################################################################

  private

  public hf_gradient
  public grd2_rhf_compute_data_t
  public grd2_uhf_compute_data_t

!###############################################################################

contains

!###############################################################################

  subroutine hf_gradient_C(c_handle) bind(C, name="hf_gradient")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call hf_gradient(inf)
  end subroutine hf_gradient_C

!###############################################################################

  subroutine hf_gradient(infos)
    use io_constants, only: iw
    use grd1, only: print_gradient
    use dft, only: dft_initialize, dftclean, dftder
    use strings, only: Cstring, fstring
    use mod_dft_molgrid, only: dft_grid_t
    use util, only: measure_time
    use printing, only: print_module_info

    implicit none

    type(information), target, intent(inout) :: infos

    type(basis_Set), pointer :: basis
    logical :: dft
    type(dft_grid_t) :: molGrid

!   3. LOG: Write: Main output file
    open (unit=IW, file=infos%log_filename, position="append")
!
    call print_module_info('HF_DFT_Gradient','Computing Gradient of HF/DFT')

    dft = infos%control%hamilton == 20

!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

!   Compute 1e gradient

    write(iw,"(/' ..... Beginning Gradient Calculation...'/)")
    call flush(iw)

    call hf_1e_grad(infos, basis)

    write(iw,"(' ..... End Of 1-Eelectron Gradient ......')")

    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

!   Calculate DFT contribution
    if (dft) then
      write(iw,"(/' ..... Gradient XC terms...'/)")
      call dft_initialize(infos, basis, molGrid, verbose=.true.)
      call dftder(basis, infos, molGrid)
      call dftclean(infos)
      write(iw,"(' ..... End of gradient XC terms......')")
      call measure_time(print_total=1, log_unit=iw)
      call flush(iw)
    end if

!   Compute 2e gradient
    call hf_2e_grad(basis, infos)

    write(iw, fmt="(' ...... End Of 2-Electron Gradient ......')")

    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)


!   Print out gradient
    call print_gradient(infos)

!   Print timings
    call measure_time(print_total=1, log_unit=iw)
    close(iw)

  end subroutine hf_gradient

!-------------------------------------------------------------------------------

  subroutine hf_1e_grad(infos, basis)

    use types, only: information
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use util, only: measure_time
    use messages, only: show_message, WITH_ABORT
    use constants, only: tol_int
    use grd1, only: eijden, grad_nn, grad_ee_overlap, &
            grad_ee_kinetic, grad_en_hellman_feynman, &
            grad_en_pulay, grad_1e_ecp
    use qmmm_mod, only: grad_esp_qmmm

    implicit none

    character(len=*), parameter :: subroutine_name = "hf_1e_grad"

    type(information), intent(inout) :: infos
    type(basis_set), intent(inout) :: basis

    real(kind=dp), allocatable :: dens(:)
    real(kind=dp) :: tol
    integer :: nbf, nbf_tri, ok
    integer(4) :: status

    real(kind=dp), pointer :: dmat_a(:), dmat_b(:)

    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2

    tol = tol_int*log(10.0_dp)

!   initial memory allocation
    allocate(dens(nbf_tri), source=0.0_dp, stat=ok)
    if(ok/=0) call show_message('Cannot allocate memory', WITH_ABORT)

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a, status)
    call check_status(status, module_name, subroutine_name, OQP_DM_A)

    if (infos%control%scftype>=2) then
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b, status)
      call check_status(status, module_name, subroutine_name, OQP_DM_B)
    endif

    associate(grad  => infos%atoms%grad, &
              xyz   => infos%atoms%xyz, &
              zn    => infos%atoms%zn - infos%basis%ecp_zn_num)

!     Zero out gradient
      grad = 0.0d0

!     Nuclear repulsion force
      call grad_nn(infos%atoms, infos%basis%ecp_zn_num)

!     Obtain Lagrangian matrix (`dens`)
      call eijden(dens, nbf, infos)

!     Overlap gradient
      call grad_ee_overlap(basis, dens, grad, logtol=tol)


!     Compute total density matrix, discard Lagrangian
      dens = dmat_a
      if (infos%control%scftype>=2) dens = dens + dmat_b

!     Kinetic gradient
      call grad_ee_kinetic(basis, dens, grad, logtol=tol)

!     Hellmann-Feynman force
      call grad_en_hellman_feynman(basis, xyz, zn, dens, grad, logtol=tol)

!     Pulay force
      call grad_en_pulay(basis, xyz, zn, dens, grad, logtol=tol)

!     Effective core potential gradient
      call grad_1e_ecp(infos, basis, xyz, dens, grad, logtol=tol)
!     QM/MM force
!      if(infos%control%qmmm_flag) call grad_esp_qmmm(infos, dens, grad,logtol=tol)

    end associate

   end subroutine

!###############################################################################

!> @brief The driver for the two electron gradient
  subroutine hf_2e_grad(basis, infos)

    use basis_tools, only: basis_set
    use precision, only: dp
    use oqp_tagarray_driver
    use messages, only: show_message, WITH_ABORT
    use types, only: information

    implicit none

    character(len=*), parameter :: subroutine_name = "hf_2e_grad"

    type(information), target, intent(inout) :: infos
    type(basis_set) :: basis

    logical :: urohf
    real(kind=dp) :: hfscale

    integer :: ok
    real(kind=dp), allocatable :: de(:,:)
    class(grd2_compute_data_t), allocatable :: gcomp

    ! tagarray
    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:)
    integer(4) :: status

    urohf = infos%control%scftype >= 2
    hfscale = 1.0d0
    if(infos%control%hamilton.ge.20) hfscale = infos%dft%hfscale

    allocate(de(3,ubound(infos%atoms%zn,1)), &
            source=0.0d0, &
            stat=ok)

    if(ok/=0) call show_message('cannot allocate memory', WITH_ABORT)

    if (urohf) then

      call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a, status)
      call check_status(status, module_name, subroutine_name, OQP_DM_A)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b, status)
      call check_status(status, module_name, subroutine_name, OQP_DM_B)

      gcomp = &
        grd2_uhf_compute_data_t( da = dmat_a &
                               , db = dmat_b &
                               , hfscale = hfscale &
                               , nbf = basis%nbf &
        )
    else

      call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a, status)
      call check_status(status, module_name, subroutine_name, OQP_DM_A)

      gcomp = &
        grd2_rhf_compute_data_t( da = dmat_a &
                               , hfscale = hfscale &
                               , nbf = basis%nbf &
        )
    end if

    call gcomp%init()

    select type (gcomp)
    class is (grd2_hf_compute_data_t)
      call gcomp%build_cart(basis)
    end select

    call grd2_driver(infos, basis, de, gcomp)
    infos%atoms%grad = infos%atoms%grad + de

    call gcomp%clean()

  end subroutine

!###############################################################################

  !> @brief Build Cartesian-effective (bfnrm-folded) densities + Cartesian
  !>        offsets for the 2e gradient under HARMONIC_ACTIVE. The 2-particle
  !>        density factorizes, so get_density uses the SAME formula on the
  !>        expanded density (no separate bfnrm factor).
  subroutine grd2_hf_build_cart(this, basis)
    class(grd2_hf_compute_data_t), intent(inout) :: this
    type(basis_set), intent(in) :: basis
    real(kind=dp), allocatable :: tmp(:,:)
    integer, allocatable :: dummy_off(:)
    integer :: ncb

    if (.not. HARMONIC_ACTIVE) return

    tmp = this%d2a
    call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
    call build_cart_density(basis, tmp, this%d2a_cart, this%cart_off, this%nbf_cart)

    if (allocated(this%d2b)) then
      tmp = this%d2b
      call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
      call build_cart_density(basis, tmp, this%d2b_cart, dummy_off, ncb)
    end if
  end subroutine grd2_hf_build_cart

!###############################################################################

  subroutine grd2_rhf_compute_data_t_init(this)
    use messages, only: show_message, WITH_ABORT
    use mathlib, only: unpack_matrix
    implicit none
    class(grd2_rhf_compute_data_t), target, intent(inout) :: this
    integer :: iok

    call this%clean()

    allocate(this%d2a(this%nbf,this%nbf), stat=iok, source=0.0d0)
    if (iok/=0) call show_message('cannot allocate memory', WITH_ABORT)
    call unpack_matrix(this%da, this%d2a)

  end subroutine

!###############################################################################

  subroutine grd2_uhf_compute_data_t_init(this)
    use messages, only: show_message, WITH_ABORT
    use mathlib, only: unpack_matrix
    implicit none
    class(grd2_uhf_compute_data_t), target, intent(inout) :: this
    integer :: iok

    call this%clean()

    allocate(this%d2a(this%nbf,this%nbf), this%d2b(this%nbf,this%nbf), stat=iok, source=0.0d0)
    if (iok/=0) call show_message('cannot allocate memory', WITH_ABORT)
    call unpack_matrix(this%da, this%d2a)
    call unpack_matrix(this%db, this%d2b)
    this%d2a = this%d2a + this%d2b
    this%d2b = this%d2a - 2*this%d2b

  end subroutine

!###############################################################################

  subroutine grd2_rhf_compute_data_t_clean(this)
    implicit none
    class(grd2_rhf_compute_data_t), target, intent(inout) :: this
    if (allocated(this%d2a)) deallocate(this%d2a)
  end subroutine

!###############################################################################

  subroutine grd2_uhf_compute_data_t_clean(this)
    implicit none
    class(grd2_uhf_compute_data_t), target, intent(inout) :: this
    if (allocated(this%d2a)) deallocate(this%d2a)
    if (allocated(this%d2b)) deallocate(this%d2b)
  end subroutine

!###############################################################################

!> @brief This routine forms the product of density
!>        matrices for use in forming the two electron
!>        gradient. Valid for closed and open shell SCF.
!> @note  dabmax is computed from the unnormalized density products (the
!>        historic screening convention); the basis normalization enters only
!>        the stored block, through norm factors hoisted out of the inner
!>        loops.
  subroutine grd2_rhf_compute_data_t_get_density(this, basis, id, dab, dabmax)

    implicit none

    class(grd2_rhf_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: coulfact, xcfact, df1, dq1, bfn
    logical :: do_exchange
    integer :: i, j, k, l
    integer :: loc(4)
    integer :: nbf(4)
    real(kind=dp), pointer :: ab(:,:,:,:)
    real(kind=dp), pointer :: d2a(:,:)
    integer :: i1, j1, k1, l1
    logical :: usecart

    coulfact = 4*this%coulscale
    xcfact = this%hfscale
    do_exchange = xcfact/=0.0_dp

    ! Under HARMONIC_ACTIVE the derivative ERIs are Cartesian, so contract
    ! against the Cartesian-effective density with Cartesian offsets; the
    ! bfnrm factor is already folded into d2a_cart. Otherwise use the
    ! spherical/Cartesian density with shell offsets and apply bfnrm.
    usecart = HARMONIC_ACTIVE
    if (usecart) then
      d2a => this%d2a_cart
      loc = this%cart_off(id) - 1
      nbf = NUM_CART_BF(basis%am(id))
    else
      d2a => this%d2a
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
            df1 = coulfact*d2a(i1,j1)*d2a(k1,l1)
            if (do_exchange) then
              dq1 = d2a(i1,k1)*d2a(j1,l1) &
                  + d2a(i1,l1)*d2a(j1,k1)
              df1 = df1-xcfact*dq1
            end if
            dabmax = max(dabmax, abs(df1))
            bfn = 1.0_dp
            if (.not. usecart) bfn = product(basis%bfnrm([i1,j1,k1,l1]))
            ab(l,k,j,i) = df1*bfn
          end do
        end do
      end do
    end do
  end subroutine grd2_rhf_compute_data_t_get_density

!###############################################################################

!> @brief This routine forms the product of density
!>        matrices for use in forming the two electron
!>        gradient. Valid for closed and open shell SCF.
!> @note  dabmax is computed from the unnormalized density products (the
!>        historic screening convention); the basis normalization enters only
!>        the stored block, through norm factors hoisted out of the inner
!>        loops.
  subroutine grd2_uhf_compute_data_t_get_density(this, basis, id, dab, dabmax)

    implicit none

    class(grd2_uhf_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: coulfact, xcfact, df1, dq1, bfn
    logical :: do_exchange
    integer :: i, j, k, l
    integer :: loc(4)
    integer :: nbf(4)
    real(kind=dp), pointer :: ab(:,:,:,:)
    real(kind=dp), pointer :: d2a(:,:), d2b(:,:)
    integer :: i1, j1, k1, l1
    logical :: usecart

    coulfact = 4*this%coulscale
    xcfact = this%hfscale
    do_exchange = xcfact/=0.0_dp

    usecart = HARMONIC_ACTIVE
    if (usecart) then
      d2a => this%d2a_cart
      d2b => this%d2b_cart
      loc = this%cart_off(id) - 1
      nbf = NUM_CART_BF(basis%am(id))
    else
      d2a => this%d2a
      d2b => this%d2b
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
            df1 = coulfact*d2a(i1,j1)*d2a(k1,l1)
            if (do_exchange) then
              dq1 = d2a(i1,k1)*d2a(j1,l1) &
                  + d2a(i1,l1)*d2a(j1,k1) &
                  + d2b(i1,k1)*d2b(j1,l1) &
                  + d2b(i1,l1)*d2b(j1,k1)
              df1 = df1-xcfact*dq1
            end if
            dabmax = max(dabmax, abs(df1))
            bfn = 1.0_dp
            if (.not. usecart) bfn = product(basis%bfnrm([i1,j1,k1,l1]))
            ab(l,k,j,i) = df1*bfn
          end do
        end do
      end do
    end do
  end subroutine grd2_uhf_compute_data_t_get_density

end module hf_gradient_mod
