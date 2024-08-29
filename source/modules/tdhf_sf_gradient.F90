module tdhf_sf_gradient_mod

  use precision, only: dp
  use grd2, only: grd2_driver, grd2_compute_data_t
  use basis_tools, only: basis_set
  use types, only: information

  implicit none

  character(len=*), parameter :: module_name = "tdhf_sf_gradient_mod"

  public tdhf_sf_gradient

  type, extends(grd2_compute_data_t) :: grd2_sf_compute_data_t
    real(kind=dp), pointer :: d2(:,:,:) => null()
    real(kind=dp), pointer :: p2(:,:,:) => null()
    real(kind=dp), pointer :: v2(:,:) => null()
    integer :: nbf = 0
  contains
    procedure :: init => grd2_sf_compute_data_t_init
    procedure :: clean => grd2_sf_compute_data_t_clean
    procedure :: get_density => grd2_sf_compute_data_t_get_density
  end type

contains

  subroutine sf_gradient_C(c_handle) bind(C, name="tdhf_sf_gradient")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_sf_gradient(inf)
  end subroutine sf_gradient_c

  subroutine tdhf_sf_gradient(infos)
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
    use tdhf_sf_lib, only: &
      sfrorhs, &
      sfromcal, sfrogen, sfrolhs, &
      pcgb, sfropcal, sfrowcal
    use dft, only: dft_initialize, dftclean
    use mathlib, only: symmetrize_matrix
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_tdxc_grad, only: utddft_xc_gradient
    use mathlib, only: unpack_matrix
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_sf_gradient"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos

    integer :: s_size

    integer :: nbf, nbf_tri
    logical :: roref = .false.

    type(dft_grid_t) :: molGrid

  ! General data
    logical :: dft
    integer :: scf_type, mol_mult

    real(kind=dp), allocatable :: p(:,:,:), v(:,:,:), d(:,:,:)

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      dmat_a(:), dmat_b(:), td_abxc(:,:), td_p(:,:)
    character(len=*), parameter :: tags_general(*) = (/ character(len=80) :: &
      OQP_DM_A, OQP_DM_B, OQP_td_abxc, OQP_td_p /)

    mol_mult = infos%mol_prop%mult
    if (mol_mult/=3) call show_message(&
            'SF-TDDFT are available for ROHF/UHF ref.&
            &with ONLY triplet multiplicity(mult=3)', with_abort)

    scf_type = infos%control%scftype
    if (scf_type==3) roref = .true.

    dft = infos%control%hamilton == 20

  ! Files open
    open (unit=IW, file=infos%log_filename, position="append")
  !
    call print_module_info('SF_Grad','Computing Gradient of SF-TDDFT')
!
    write(iw,'(/5X,"Gradient options"/&
                &5X,18("-")/&
                &5X,"Target State: ",I8/&
                &,5X,"*Note that the ground state of SF-TDDFT is 1.*"/)')&
                & infos%tddft%target_state

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    call tagarray_get_data(infos%dat, OQP_td_abxc, td_abxc)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)

  ! Allocate H, S ,T and D matrices
    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2
    s_size = (basis%nshell**2+basis%nshell)/2

!   Compute 1e gradient
    call flush(iw)

    call sf_1e_grad(infos, basis)

    write(iw,"(' ..... End Of 1-Eelectron Gradient ......')")
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

    allocate(d(nbf,nbf,2), source=0.0d0)
    allocate(p(nbf,nbf,2), source=0.0d0)

    call unpack_matrix(td_p(:,1), p(:,:,1))
    call unpack_matrix(td_p(:,2), p(:,:,2))

    call unpack_matrix(dmat_a, d(:,:,1))
    call unpack_matrix(dmat_b, d(:,:,2))

!   Compute xc gradient
    if (dft) then
      call dft_initialize(infos, basis, molGrid, verbose=.true.)

      call utddft_xc_gradient(basis=basis, &
           molGrid=molGrid, &
           dedft=infos%atoms%grad, &
           da=d(:,:,1), &
           db=d(:,:,2), &
           pa=p(:,:,1:1), &
           pb=p(:,:,2:2), &
           nmtx=1, &
           !threshold=1.0d-15, &
           threshold=0.0d0, &
           infos=infos)

      call dftclean(infos)
      call measure_time(print_total=1, log_unit=iw)
      call flush(iw)
    end if

!   Compute 2e gradient
    allocate(v(nbf,nbf,2), source=0.0d0)
    v(:,:,1) = td_abxc
    call sf_2e_grad(basis, infos, d, p, v(:,:,1))

    call print_gradient(infos)

!   Print timings
    call measure_time(print_total=1, log_unit=iw)

    close(iw)

  end subroutine tdhf_sf_gradient

!###############################################################################

  subroutine sf_1e_grad(infos, basis)

    use oqp_tagarray_driver
    use types, only: information
    use basis_tools, only: basis_set
    use util, only: measure_time
    use messages, only: show_message, WITH_ABORT
    use precision, only: dp
    use constants, only: tol_int
    use grd1, only: eijden, print_gradient, &
            grad_nn, grad_ee_overlap, &
            grad_ee_kinetic, grad_en_hellman_feynman, grad_en_pulay

    use mathlib, only: symmetrize_matrix

    implicit none

    character(len=*), parameter :: subroutine_name = "sf_1e_grad"

    type(information), intent(inout) :: infos
    type(basis_set), intent(inout) :: basis

    real(kind=dp), allocatable :: dens(:)
    real(kind=dp) :: tol
    integer :: nbf, nbf_tri, ok

    ! tagarray
    real(kind=dp), pointer :: dmat_a(:), dmat_b(:), wao(:), td_p(:,:)
    character(len=*), parameter :: tags_general(4) = (/ character(len=80) :: &
      OQP_DM_A, OQP_DM_B, OQP_WAO, OQP_td_p /)

    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2

    tol = tol_int*log(10.0_dp)

!   initial memory allocation
    allocate(dens(nbf_tri), source=0.0_dp, stat=ok)
    if(ok/=0) call show_message('Cannot allocate memory', WITH_ABORT)

    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
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
      dens = dens + 2*w

!     Overlap gradient
      call grad_ee_overlap(basis, dens, grad, logtol=tol)

!     Compute total density matrix, discard Lagrangian
      dens = dmat_a + p(:,1)
      if (infos%control%scftype>=2) then
        dens = dens + dmat_b + p(:,2)
      end if

!     Hellmann-Feynman force
      call grad_en_hellman_feynman(basis, xyz, zn, dens, grad, logtol=tol)

!     KE gradient
      call grad_ee_kinetic(basis, dens, grad, logtol=tol)

!     Pulay force
      call grad_en_pulay(basis, xyz, zn, dens, grad, logtol=tol)

    end associate

   end subroutine

!###############################################################################

!> @brief The driver for the two electron gradient
  subroutine sf_2e_grad(basis, infos, d, p, v)

    use basis_tools, only: basis_set
    use precision, only: dp
    use messages, only: show_message, WITH_ABORT
    use types, only: information

    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set) :: basis
    real(kind=dp), contiguous, target :: p(:,:,:), d(:,:,:), v(:,:)

    logical :: urohf, dft
    real(kind=dp) :: scale_exch  !> HF scale in Reference
    real(kind=dp) :: scale_exch2 !> HF scale in Response

    integer :: ok
    real(kind=dp), allocatable :: de(:,:)
    class(grd2_compute_data_t), allocatable :: gcomp

    dft = infos%control%hamilton == 20 ! dft or hf
    urohf = infos%control%scftype >= 2

    scale_exch = 1.0_dp
    scale_exch2 = 1.0_dp
    if (dft) then
      scale_exch = infos%dft%HFscale
      scale_exch2 = infos%tddft%HFscale
    end if

    allocate(de(3,ubound(infos%atoms%zn,1)), &
            source=0.0d0, &
            stat=ok)

    if(ok/=0) call show_message('cannot allocate memory', WITH_ABORT)

    write(*, '(/7x,"Fitting parameters")')
    if (.not.infos%dft%cam_flag) then
      write(*, '(10x,"Exact HF exchange:")')
      write(*, '(5x,"Reference: |", t20, f6.3, t29, "|")') scale_exch
      write(*, '(5x,"Response:  |", t20, f6.3, t29, "|")') scale_exch2
    else
      write(*, '(10x,"CAM parametres:")')
      write(*, '(16x,"|   alpha   |    beta   |     mu    |")')
      write(*, '(5x,"Reference: |", t20, f6.3, t29, "|", t32, f6.3, t41, "|", t44, f6.3, t53, "|")') &
         infos%dft%cam_alpha, infos%dft%cam_beta, infos%dft%cam_mu
      write(*, '(5x,"Response:  |", t20, f6.3, t29, "|", t32, f6.3, t41, "|", t44, f6.3, t53, "|")') &
         infos%tddft%cam_alpha, infos%tddft%cam_beta, infos%tddft%cam_mu
    end if
    write(*, '(10x,"Spin-pair coupling parametres:")')
    write(*, '(16x,"|   CO-CO   |   OV-OV   |   CO-OV   |")')
    write(*, '(16x,"|", t20, f6.3, t29, "|", t32, f6.3, t41, "|", t44, f6.3, t53, "|")') &
       infos%tddft%HFscale, infos%tddft%HFscale, infos%tddft%HFscale

    gcomp =  grd2_sf_compute_data_t( d2 = d &
                                   , p2 = p &
                                   , v2 = v &
                                   , nbf = basis%nbf )

    call gcomp%init()

    call grd2_driver(infos, basis, de, gcomp, &
                     cam = dft.and.infos%dft%cam_flag, &
                     alpha = infos%tddft%cam_alpha, &
                     beta = infos%tddft%cam_beta, &
                     mu = infos%tddft%cam_mu)

    infos%atoms%grad = infos%atoms%grad + de

    call gcomp%clean()

  end subroutine

!###############################################################################

  subroutine grd2_sf_compute_data_t_init(this)
    implicit none
    class(grd2_sf_compute_data_t), target, intent(inout) :: this

    call this%clean()

    this%d2(:,:,1) = this%d2(:,:,1) +   this%d2(:,:,2)
    this%d2(:,:,2) = this%d2(:,:,1) - 2*this%d2(:,:,2)

    this%p2(:,:,1) = this%p2(:,:,1) +   this%p2(:,:,2)
    this%p2(:,:,2) = this%p2(:,:,1) - 2*this%p2(:,:,2)

  end subroutine

!###############################################################################

  subroutine grd2_sf_compute_data_t_clean(this)
    implicit none
    class(grd2_sf_compute_data_t), target, intent(inout) :: this
  end subroutine

!###############################################################################

!> @brief This routine forms the product of density
!>        matrices for use in forming the two electron
!>        gradient. Valid for closed and open shell SCF.
  subroutine grd2_sf_compute_data_t_get_density(this, basis, id, dab, dabmax)

    implicit none

    class(grd2_sf_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: df1, dq1, dt2
    real(kind=dp) :: coulfact, xcfact, xcfact2
    integer :: i, j, k, l
    integer :: loc(4)
    integer :: nbf(4)
    real(kind=dp), pointer :: ab(:,:,:,:)
    integer :: i1, j1, k1, l1

    coulfact = 4*this%coulscale
    xcfact = this%hfscale
    xcfact2 = this%hfscale2
    dabmax = 0
    loc = basis%ao_offset(id)-1

    nbf = basis%naos(id)

    ab(1:nbf(4),1:nbf(3),1:nbf(2),1:nbf(1)) => dab(1:product(nbf))

    do i = 1, nbf(1)
      i1 = loc(1) + i

      do j = 1, nbf(2)
        j1 = loc(2) + j

        do k = 1, nbf(3)
          k1 = loc(3) + k

          do l = 1, nbf(4)
            l1 = loc(4) + l
            df1 = (this%d2(i1,j1,1)+this%p2(i1,j1,1))*this%d2(k1,l1,1) &
                +  this%d2(i1,j1,1)                  *this%p2(k1,l1,1)
            df1 = df1 * coulfact

            if (xcfact /= 0.0_dp .or. xcfact2 /= 0.0_dp) then
              dq1 = (this%d2(i1,k1,1)+this%p2(i1,k1,1))*this%d2(j1,l1,1) &
                  +  this%d2(i1,k1,1)                  *this%p2(j1,l1,1) &
                  + (this%d2(i1,l1,1)+this%p2(i1,l1,1))*this%d2(j1,k1,1) &
                  +  this%d2(i1,l1,1)                  *this%p2(j1,k1,1) &
                  + (this%d2(i1,k1,2)+this%p2(i1,k1,2))*this%d2(j1,l1,2) &
                  +  this%d2(i1,k1,2)                  *this%p2(j1,l1,2) &
                  + (this%d2(i1,l1,2)+this%p2(i1,l1,2))*this%d2(j1,k1,2) &
                  +  this%d2(i1,l1,2)                  *this%p2(j1,k1,2)
              dt2 = this%v2(i1,k1)*this%v2(j1,l1) &
                  + this%v2(k1,i1)*this%v2(l1,j1) &
                  + this%v2(i1,l1)*this%v2(j1,k1) &
                  + this%v2(l1,i1)*this%v2(k1,j1)

              df1 = df1-xcfact*dq1-xcfact2*2.0_dp*dt2
            end if
            dabmax = max(dabmax, abs(df1))
            ab(l,k,j,i) = df1*product(basis%bfnrm([i1,j1,k1,l1]))
          end do
        end do
      end do
    end do
  end subroutine grd2_sf_compute_data_t_get_density

!###############################################################################

end module tdhf_sf_gradient_mod
