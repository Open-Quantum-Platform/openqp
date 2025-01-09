module tdhf_gradient_mod

  use precision, only: dp
  use grd2, only: grd2_driver, grd2_compute_data_t
  use basis_tools, only: basis_set
  use types, only: information
  use io_constants, only: iw

!###############################################################################

  implicit none

!###############################################################################

  character(len=*), parameter :: module_name = "tdhf_gradient_mod"

!###############################################################################

  public tdhf_gradient

!###############################################################################

  type, extends(grd2_compute_data_t) :: grd2_tdhf_compute_data_t
    real(kind=dp), pointer :: d2(:,:,:) => null()
    real(kind=dp), pointer :: p2(:,:,:) => null()
    real(kind=dp), pointer :: xpy2(:,:,:) => null()
    real(kind=dp), pointer :: xmy2(:,:,:) => null()
    integer :: nbf = 0
  contains
    procedure :: init => grd2_tdhf_compute_data_t_init
    procedure :: clean => grd2_tdhf_compute_data_t_clean
    procedure :: get_density => grd2_tdhf_compute_data_t_get_density
  end type

contains

!###############################################################################

  subroutine tdhf_gradient_C(c_handle) bind(C, name="tdhf_gradient")

    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use strings, only: Cstring

    implicit none

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call tdhf_gradient(inf)

  end subroutine tdhf_gradient_c

!###############################################################################

  subroutine tdhf_gradient(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver

    use types, only: information
    use strings, only: Cstring, fstring
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort

    use grd1, only: eijden, print_gradient
    use util, only: measure_time
    use tdhf_lib, only: &
      iatogen, mntoia
    use dft, only: dft_initialize, dftclean
    use mathlib, only: symmetrize_matrix, orthogonal_transform
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_tdxc_grad, only: tddft_xc_gradient
    use mathlib, only: unpack_matrix
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_gradient"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos

    integer :: nbf, nocc

    type(dft_grid_t) :: molGrid

  ! General data
    logical :: dft
    integer :: scf_type, mol_mult

    real(kind=dp), allocatable :: p(:,:,:), d(:,:,:), xpy2(:,:,:), xmy2(:,:,:), wrk1(:,:), wrk2(:,:)

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      dmat_a(:), mo_a(:,:), td_p(:,:), xpy(:), xmy(:)
    character(len=*), parameter :: tags_general(*) = [ character(len=80) :: &
      OQP_TD_XPY, OQP_TD_XMY, OQP_DM_A, OQP_VEC_MO_A, OQP_TD_P ]

    mol_mult = infos%mol_prop%mult
    if (mol_mult/=1) call show_message(&
            'RPA-TDDFT are only available for RHF reference', with_abort)

    scf_type = infos%control%scftype
    if (scf_type/=1) error stop

    dft = infos%control%hamilton == 20

  ! Files open
    open (unit=IW, file=infos%log_filename, position="append")
  !
    call print_module_info('TDHF_Grad','Computing Grdient of TDDFT')
!
    write(iw,'(/5X,"Gradient options"/&
                &5X,18("-")/&
                &5X,"Target State: ",I8/)') infos%tddft%target_state
   ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)
    call tagarray_get_data(infos%dat, OQP_td_xpy, xpy)
    call tagarray_get_data(infos%dat, OQP_td_xmy, xmy)

  ! Allocate H, S ,T and D matrices
    nbf = basis%nbf
    nocc = infos%mol_prop%nocc

!   Compute 1e gradient
    call flush(iw)

    call tdhf_1e_grad(infos, basis)
    !call print_gradient(infos)

    write(iw,"(' ..... End Of 1-Eelectron Gradient ......')")
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

    allocate(d(nbf,nbf,2), &
             p(nbf,nbf,2), &
             xpy2(nbf,nbf,2), &
             xmy2(nbf,nbf,2), &
             wrk1(nbf,nbf), &
             wrk2(nbf,nbf), &
             source=0.0d0)

    call iatogen(xpy, wrk1, nocc, nocc)
    call symmetrize_matrix(wrk1, nbf)
    wrk1 = wrk1*0.5
    call orthogonal_transform('t', nbf, mo_a, wrk1, xpy2(:,:,1), wrk2)

    call iatogen(xmy, wrk1, nocc, nocc)
    call orthogonal_transform('t', nbf, mo_a, wrk1, xmy2(:,:,1), wrk2)


    call unpack_matrix(td_p(:,1), p(:,:,1))
    call unpack_matrix(dmat_a, d(:,:,1))

!   Compute xc gradient
    if (dft) then
      call dft_initialize(infos, basis, molGrid)

      d(:,:,2) = d(:,:,1)
      p(:,:,2) = p(:,:,1)
      xpy2(:,:,2) = xpy2(:,:,1)
      call tddft_xc_gradient(basis=basis, &
             molGrid=molGrid, &
             dedft=infos%atoms%grad, &
             da=d(:,:,1), &
             pa=p(:,:,1:1), &
             xa=xpy2(:,:,1:1), &
             nmtx=1, &
             threshold=1.0d-14, &
             infos=infos)
      call dftclean(infos)
      call measure_time(print_total=1, log_unit=iw)
      call flush(iw)
    end if

!   Compute 2e gradient
    call tdhf_2e_grad(basis, infos, d(:,:,1:1), p(:,:,1:1), xpy2(:,:,1:1), xmy2(:,:,1:1))

    call print_gradient(infos)

!   Print timings
    call measure_time(print_total=1, log_unit=iw)

    close(iw)

  end subroutine tdhf_gradient

!###############################################################################

  subroutine tdhf_1e_grad(infos, basis)

    use oqp_tagarray_driver
    use types, only: information
    use basis_tools, only: basis_set
    use util, only: measure_time
    use messages, only: show_message, WITH_ABORT
    use precision, only: dp
    use constants, only: tol_int
    use grd1, only: eijden, print_gradient, &
            grad_nn, grad_ee_overlap, &
            grad_ee_kinetic, grad_en_hellman_feynman, grad_en_pulay, grad_1e_ecp

    use mathlib, only: symmetrize_matrix

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_1e_grad"

    type(information), intent(inout) :: infos
    type(basis_set), intent(inout) :: basis

    real(kind=dp), allocatable :: dens(:)
    real(kind=dp) :: tol
    integer :: nbf, nbf2, ok

    ! tagarray
    real(kind=dp), pointer :: dmat_a(:), wao(:), td_p(:,:)
    character(len=*), parameter :: tags_general(*) = (/ character(len=80) :: &
      OQP_DM_A, OQP_WAO, OQP_TD_P /)

    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2

    tol = tol_int*log(10.0_dp)

!   initial memory allocation
    allocate(dens(nbf2), source=0.0_dp, stat=ok)
    if(ok/=0) call show_message('Cannot allocate memory', WITH_ABORT)

    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_WAO, wao)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)

    associate( grad   => infos%atoms%grad &
             , xyz    => infos%atoms%xyz &
             , zn     => infos%atoms%zn &
             , p      => td_p &
             , w      => wao &
             , urohf  => infos%control%scftype>=2 &
        )

!     Zero out gradient
      grad = 0.0d0
!     Nuclear repulsion force
      call grad_nn(infos%atoms)

!     Obtain Lagrangian matrix (`dens`)
      call eijden(dens, nbf, infos)

!     Add W matrix:
      dens = dens - 2*w

!     Overlap gradient
      call grad_ee_overlap(basis, dens, grad, logtol=tol)

!     Compute total density matrix, discard Lagrangian
      dens = dmat_a + 2*p(:,1)

!     Hellmann-Feynman force
      call grad_en_hellman_feynman(basis, xyz, zn, dens, grad, logtol=tol)

!     KE gradient
      call grad_ee_kinetic(basis, dens, grad, logtol=tol)

!     Pulay force
      call grad_en_pulay(basis, xyz, zn, dens, grad, logtol=tol)
!     Effective core potential gradient
      call grad_1e_ecp(basis, xyz, dens, grad, logtol=tol)

    end associate

   end subroutine

!###############################################################################

!> @brief The driver for the two electron gradient
  subroutine tdhf_2e_grad(basis, infos, d, p, xpy, xmy)

    use basis_tools, only: basis_set
    use precision, only: dp
    use messages, only: show_message, WITH_ABORT
    use types, only: information

    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set) :: basis
    real(kind=dp), contiguous, target :: p(:,:,:), d(:,:,:), xpy(:,:,:), xmy(:,:,:)

    logical :: urohf
    real(kind=dp) :: hfscale

    integer :: ok
    real(kind=dp), allocatable :: de(:,:)
    class(grd2_compute_data_t), allocatable :: gcomp

    urohf = infos%control%scftype >= 2
    hfscale = 1.0d0
    if(infos%control%hamilton.ge.20) hfscale = infos%dft%hfscale

    allocate(de(3,ubound(infos%atoms%zn,1)), &
            source=0.0d0, &
            stat=ok)

    if(ok/=0) call show_message('cannot allocate memory', WITH_ABORT)

    gcomp =  grd2_tdhf_compute_data_t( d2 = d &
                                     , p2 = p &
                                     , xpy2 = xpy &
                                     , xmy2 = xmy &
                                     , hfscale = hfscale &
                                     , nbf = basis%nbf &
      )

    call gcomp%init()

    call grd2_driver(infos, basis, de, gcomp)
    infos%atoms%grad = infos%atoms%grad + de

    call gcomp%clean()

  end subroutine

!###############################################################################

  subroutine grd2_tdhf_compute_data_t_init(this)
    !use messages, only: show_message, WITH_ABORT
    implicit none
    class(grd2_tdhf_compute_data_t), target, intent(inout) :: this

    call this%clean()


  end subroutine

!###############################################################################

  subroutine grd2_tdhf_compute_data_t_clean(this)
    implicit none
    class(grd2_tdhf_compute_data_t), target, intent(inout) :: this
  end subroutine

!###############################################################################

!> @brief Compute density factors for \Gamma term of 2-electron contribution
!>        to TD-DFT gradients
  subroutine grd2_tdhf_compute_data_t_get_density(this, basis, id, dab, dabmax)

    implicit none

    class(grd2_tdhf_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: coul, coul2, coul3, exch, exch2, exch3, exch4
    real(kind=dp) :: df1
    real(kind=dp) :: coulfact, xcfact
    integer :: i_, j_, k_, l_
    integer :: i, j, k, l
    integer :: loc(4)
    integer :: nbf(4)
    real(kind=dp), pointer :: ab(:,:,:,:)

    coulfact = 4*this%coulscale
    xcfact = this%hfscale
    dabmax = 0
    loc = basis%ao_offset(id)-1

    nbf = basis%naos(id)

    ab(1:nbf(4),1:nbf(3),1:nbf(2),1:nbf(1)) => dab(1:product(nbf))

    associate( p => this%p2(:,:,1) &
             , d => this%d2(:,:,1) &
             , xpy => this%xpy2(:,:,1) &
             , xmy => this%xmy2(:,:,1) &
             )

    do i_ = 1, nbf(1)
      i = loc(1) + i_

      do j_ = 1, nbf(2)
        j = loc(2) + j_

        do k_ = 1, nbf(3)
          k = loc(3) + k_

          do l_ = 1, nbf(4)
            l = loc(4) + l_

            coul = &
            +   p(j,i) * d(l,k) &
            +   p(l,k) * d(j,i)

            coul2 = &
            + xpy(i,j) * xpy(k,l)

            coul3 = d(i,j) *   d(k,l)

            coul = 2*coul + 8*coul2 + coul3

            df1 = coulfact*coul

            if (xcfact/=0.0_dp) then
              exch = &
              +   p(i,k) *   d(l,j) &
              +   p(j,k) *   d(l,i) &
              +   p(l,i) *   d(k,j) &
              +   p(l,j) *   d(k,i)

              exch2 = &
              + xpy(k,i) * xpy(l,j) &
              + xpy(l,i) * xpy(k,j)

              exch3 = &
              + (xmy(i,l)-xmy(l,i)) * (xmy(j,k)-xmy(k,j)) &
              + (xmy(j,l)-xmy(l,j)) * (xmy(i,k)-xmy(k,i))

              exch4 =    d(i,k) *   d(j,l) &
                     +   d(i,l) *   d(j,k)

              exch = 2*exch + 8*exch2 + 2*exch3 + exch4

              df1 = df1 - xcfact*exch
            end if
            dabmax = max(dabmax, abs(df1))
            ab(l_,k_,j_,i_) = df1*product(basis%bfnrm([i,j,k,l]))
          end do
        end do
      end do
    end do

    end associate

  end subroutine grd2_tdhf_compute_data_t_get_density

!###############################################################################

end module tdhf_gradient_mod
