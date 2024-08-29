module tdhf_mrsf_gradient_mod

  use precision, only: dp
  use grd2, only: grd2_driver, grd2_compute_data_t
  use basis_tools, only: basis_set
  use printing, only: print_module_info

  implicit none

  character(len=*), parameter :: module_name = "tdhf_mrsf_gradient_mod"

  public tdhf_mrsf_gradient

  type, extends(grd2_compute_data_t) :: grd2_mrsf_compute_data_t
    real(kind=dp), pointer :: d2(:,:,:) => null()
    real(kind=dp), pointer :: p2(:,:,:) => null()
    real(kind=dp), pointer :: spc2(:,:,:) => null()
    integer :: nbf = 0
    integer :: mrst = 1
    real(kind=dp), dimension(3) :: spcscale = [0.0_dp, 0.0_dp, 0.0_dp]
  contains
    procedure :: init => grd2_mrsf_compute_data_t_init
    procedure :: clean => grd2_mrsf_compute_data_t_clean
    procedure :: get_density => grd2_mrsf_compute_data_t_get_density
  end type

contains

  subroutine tdhf_mrsf_gradient_C(c_handle) bind(C, name="tdhf_mrsf_gradient")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_mrsf_gradient(inf)
  end subroutine tdhf_mrsf_gradient_c

  subroutine tdhf_mrsf_gradient(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver

    use types, only: information
    use strings, only: Cstring, fstring
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort

    use grd1, only: eijden, print_gradient
    use util, only: measure_time
    use dft, only: dft_initialize, dftclean
    use mathlib, only: symmetrize_matrix
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_tdxc_grad, only: utddft_xc_gradient
    use mathlib, only: unpack_matrix
    use tdhf_sf_gradient_mod, only: sf_1e_grad, sf_2e_grad

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_mrsf_gradient"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos

    integer :: s_size

    integer :: nbf, nbf2
    integer :: mrst
    logical :: roref = .false.

    type(dft_grid_t) :: molGrid

  ! General data
    logical :: dft
    integer :: scf_type, mol_mult

    real(kind=dp), allocatable :: p(:,:,:), v(:,:,:), d(:,:,:), spc(:,:,:)

    ! tagarray
    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:), td_mrsf_density(:,:,:), td_abxc(:,:), td_p(:,:)
    character(len=*), parameter :: tags_general(*) = (/ character(len=80) :: &
      OQP_DM_A, OQP_DM_B, OQP_td_abxc, OQP_td_p /)
    character(len=*), parameter :: tags_mrsf(1) = (/ character(len=80) :: &
      OQP_td_mrsf_density /)

    mol_mult = infos%mol_prop%mult
    if (mol_mult/=3) call show_message(&
            'MRSF-TDDFT are available for ROHF/UHF ref.&
            &with ONLY triplet multiplicity(mult=3)', with_abort)

    scf_type = infos%control%scftype
    if (scf_type==3) roref = .true.

    dft = infos%control%hamilton == 20

  ! Files open
    open (unit=iw, file=infos%log_filename, position="append")
  !
    call print_module_info('MRSF_Grad','Computing Gradient of MRSF-TDDFT')
!
    write(iw,'(/5X,"Gradient options"/&
                &5X,18("-")/&
                &5X,"Target State: ",I8/&
                &,5X,"Note: MRSF-TDDFT can calculate ground state with state number 1."/)')&
                & infos%tddft%target_state

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

   ! Input parameters
    mrst = infos%tddft%mult

  ! Allocate H, S ,T and D matrices
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    s_size = (basis%nshell**2+basis%nshell)/2

!   Compute 1e gradient
    call flush(iw)

    call sf_1e_grad(infos, basis)

    write(iw,"(' ..... End Of 1-Eelectron Gradient ......')")
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

    allocate(v(nbf,nbf,2), source=0.0d0)
    allocate(d(nbf,nbf,2), source=0.0d0)
    allocate(p(nbf,nbf,2), source=0.0d0)
    allocate(spc(7,nbf,nbf), source=0.0d0)

    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    call tagarray_get_data(infos%dat, OQP_td_abxc, td_abxc)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)
    if (mrst==1 .or. mrst==3) then
      call data_has_tags(infos%dat, tags_mrsf, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_td_mrsf_density, td_mrsf_density)
    endif

    call unpack_matrix(td_p(:,1), p(:,:,1))
    call unpack_matrix(td_p(:,2), p(:,:,2))

    call unpack_matrix(dmat_a, d(:,:,1))
    call unpack_matrix(dmat_b, d(:,:,2))

    v(:,:,1) = td_abxc
    if (mrst==1 .or. mrst==3) then
      spc(1:7,:,:) = td_mrsf_density
    end if

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
    if (mrst==1 .or. mrst==3) then
      call mrsf_2e_grad(basis, infos, d, p, spc, v(:,:,1))
    else if (mrst==5) then
      call sf_2e_grad(basis, infos, d, p, v(:,:,1))
    end if

    call print_gradient(infos)

!   Print timings
    call measure_time(print_total=1, log_unit=iw)

    close(iw)

  end subroutine tdhf_mrsf_gradient

!###############################################################################

!> @brief The driver for the two electron gradient
  subroutine mrsf_2e_grad(basis, infos, d, p, spc, v)

    use basis_tools, only: basis_set
    use precision, only: dp
    use messages, only: show_message, WITH_ABORT
    use types, only: information

    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set) :: basis
    real(kind=dp), contiguous, target :: p(:,:,:), d(:,:,:), spc(:,:,:), v(:,:)

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

    write(*, '(/7x,"Fitting parameters for MRSF-TDDFT")')
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
       infos%tddft%spc_coco, infos%tddft%spc_ovov, infos%tddft%spc_coov

    gcomp = grd2_mrsf_compute_data_t( d2 = d &
                                    , p2 = p &
                                    , spc2 = spc &
                                    , nbf = basis%nbf &
                                    , hfscale = scale_exch &
                                    , hfscale2 = scale_exch2 &
                                    , spcscale = [infos%tddft%spc_coco, &
                                                  infos%tddft%spc_ovov, &
                                                  infos%tddft%spc_coov] &
                                    , mrst = infos%tddft%mult )

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

  subroutine grd2_mrsf_compute_data_t_init(this)
    implicit none
    class(grd2_mrsf_compute_data_t), target, intent(inout) :: this

    call this%clean()

    this%d2(:,:,1) = this%d2(:,:,1) +   this%d2(:,:,2)
    this%d2(:,:,2) = this%d2(:,:,1) - 2*this%d2(:,:,2)

    this%p2(:,:,1) = this%p2(:,:,1) +   this%p2(:,:,2)
    this%p2(:,:,2) = this%p2(:,:,1) - 2*this%p2(:,:,2)

  end subroutine

!###############################################################################

  subroutine grd2_mrsf_compute_data_t_clean(this)
    implicit none
    class(grd2_mrsf_compute_data_t), target, intent(inout) :: this
  end subroutine

!###############################################################################

!> @brief This routine forms the product of density
!>        matrices for use in forming the two electron
!>        gradient. Valid for closed and open shell SCF.
  subroutine grd2_mrsf_compute_data_t_get_density(this, basis, id, dab, dabmax)

    implicit none

    class(grd2_mrsf_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: xcfact, xcfact2, coulfact, df1, dq1, dt2
    real(kind=dp) :: qfspcp1, qfspcp2, qfspcp3, sgnk
    real(kind=dp) :: db1, db2, dc1, dc2, dc3, dc4, dd1, dd2, dd3, dd4
    real(kind=dp), pointer, dimension(:,:) :: &
      ball, bo2v, bo1v, bco1, bco2, co12, o21v
    integer :: i, j, k, l
    integer :: loc(4)
    integer :: nbf(4)
    real(kind=dp), pointer :: ab(:,:,:,:)
    integer :: i1, j1, k1, l1

    ! mrsf density-like
    ball => this%spc2(7,:,:)
    ! Spin-pair coupling density-like
    bo2v => this%spc2(1,:,:)
    bo1v => this%spc2(2,:,:)
    bco1 => this%spc2(3,:,:)
    bco2 => this%spc2(4,:,:)
    o21v => this%spc2(5,:,:)
    co12 => this%spc2(6,:,:)

    coulfact = 4*this%coulscale
    xcfact = this%hfscale
    xcfact2 = this%hfscale2
    qfspcp1 = this%spcscale(1)
    qfspcp2 = this%spcscale(2)
    qfspcp3 = this%spcscale(3)

    sgnk = 1.0_dp
    if (this%mrst==3) sgnk = -1.0_dp
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
              dt2 = ball(i1,k1)*ball(j1,l1) &
                  + ball(k1,i1)*ball(l1,j1) &
                  + ball(i1,l1)*ball(j1,k1) &
                  + ball(l1,i1)*ball(k1,j1)

              df1 = df1-xcfact*dq1-xcfact2*2.0_dp*dt2
            end if

            if (qfspcp1 /= 0.0_dp) then
              db1 =  co12(i1,k1)*co12(l1,j1) &
                   + co12(i1,l1)*co12(k1,j1) &
                   + co12(j1,k1)*co12(l1,i1) &
                   + co12(j1,l1)*co12(k1,i1) &
                   + co12(l1,j1)*co12(i1,k1) &
                   + co12(k1,j1)*co12(i1,l1) &
                   + co12(l1,i1)*co12(j1,k1) &
                   + co12(k1,i1)*co12(j1,l1)

              df1 = df1 + sgnk*qfspcp1*db1
            end if

            if (qfspcp2 /= 0.0_dp) then
              db2 =  o21v(i1,k1)*o21v(l1,j1) &
                   + o21v(i1,l1)*o21v(k1,j1) &
                   + o21v(j1,k1)*o21v(l1,i1) &
                   + o21v(j1,l1)*o21v(k1,i1) &
                   + o21v(l1,j1)*o21v(i1,k1) &
                   + o21v(k1,j1)*o21v(i1,l1) &
                   + o21v(l1,i1)*o21v(j1,k1) &
                   + o21v(k1,i1)*o21v(j1,l1)

              df1 = df1 - sgnk*qfspcp2*db2
            end if

            if (qfspcp3 /= 0.0_dp) then
              dc1 =  bco1(i1,k1)*bo2v(j1,l1) &
                   + bco1(i1,l1)*bo2v(j1,k1) &
                   + bco1(j1,k1)*bo2v(i1,l1) &
                   + bco1(j1,l1)*bo2v(i1,k1) &
                   + bco1(l1,j1)*bo2v(k1,i1) &
                   + bco1(k1,j1)*bo2v(l1,i1) &
                   + bco1(l1,i1)*bo2v(k1,j1) &
                   + bco1(k1,i1)*bo2v(l1,j1)

              dc2 =  bco2(i1,k1)*bo1v(j1,l1) &
                   + bco2(i1,l1)*bo1v(j1,k1) &
                   + bco2(j1,k1)*bo1v(i1,l1) &
                   + bco2(j1,l1)*bo1v(i1,k1) &
                   + bco2(l1,j1)*bo1v(k1,i1) &
                   + bco2(k1,j1)*bo1v(l1,i1) &
                   + bco2(l1,i1)*bo1v(k1,j1) &
                   + bco2(k1,i1)*bo1v(l1,j1)

              dc3 =  bo2v(i1,k1)*bco1(j1,l1) &
                   + bo2v(i1,l1)*bco1(j1,k1) &
                   + bo2v(j1,k1)*bco1(i1,l1) &
                   + bo2v(j1,l1)*bco1(i1,k1) &
                   + bo2v(l1,j1)*bco1(k1,i1) &
                   + bo2v(k1,j1)*bco1(l1,i1) &
                   + bo2v(l1,i1)*bco1(k1,j1) &
                   + bo2v(k1,i1)*bco1(l1,j1)

              dc4 =  bo1v(i1,k1)*bco2(j1,l1) &
                   + bo1v(i1,l1)*bco2(j1,k1) &
                   + bo1v(j1,k1)*bco2(i1,l1) &
                   + bo1v(j1,l1)*bco2(i1,k1) &
                   + bo1v(l1,j1)*bco2(k1,i1) &
                   + bo1v(k1,j1)*bco2(l1,i1) &
                   + bo1v(l1,i1)*bco2(k1,j1) &
                   + bo1v(k1,i1)*bco2(l1,j1)

              dd1 =  bco1(i1,j1)*bo2v(l1,k1) &
                   + bco1(i1,j1)*bo2v(k1,l1) &
                   + bco1(j1,i1)*bo2v(l1,k1) &
                   + bco1(j1,i1)*bo2v(k1,l1) &
                   + bco1(l1,k1)*bo2v(i1,j1) &
                   + bco1(k1,l1)*bo2v(i1,j1) &
                   + bco1(l1,k1)*bo2v(j1,i1) &
                   + bco1(k1,l1)*bo2v(j1,i1)

              dd2 =  bco2(i1,j1)*bo1v(l1,k1) &
                   + bco2(i1,j1)*bo1v(k1,l1) &
                   + bco2(j1,i1)*bo1v(l1,k1) &
                   + bco2(j1,i1)*bo1v(k1,l1) &
                   + bco2(l1,k1)*bo1v(i1,j1) &
                   + bco2(k1,l1)*bo1v(i1,j1) &
                   + bco2(l1,k1)*bo1v(j1,i1) &
                   + bco2(k1,l1)*bo1v(j1,i1)

              dd3 =  bo2v(i1,j1)*bco1(l1,k1) &
                   + bo2v(i1,j1)*bco1(k1,l1) &
                   + bo2v(j1,i1)*bco1(l1,k1) &
                   + bo2v(j1,i1)*bco1(k1,l1) &
                   + bo2v(l1,k1)*bco1(i1,j1) &
                   + bo2v(k1,l1)*bco1(i1,j1) &
                   + bo2v(l1,k1)*bco1(j1,i1) &
                   + bo2v(k1,l1)*bco1(j1,i1)

              dd4 =  bo1v(i1,j1)*bco2(l1,k1) &
                   + bo1v(i1,j1)*bco2(k1,l1) &
                   + bo1v(j1,i1)*bco2(l1,k1) &
                   + bo1v(j1,i1)*bco2(k1,l1) &
                   + bo1v(l1,k1)*bco2(i1,j1) &
                   + bo1v(k1,l1)*bco2(i1,j1) &
                   + bo1v(l1,k1)*bco2(j1,i1) &
                   + bo1v(k1,l1)*bco2(j1,i1)

              df1  = df1 + sgnk*qfspcp3*(-dc1-dc2-dc3-dc4 &
                                         +dd1+dd2+dd3+dd4)
            end if

            dabmax = max(dabmax, abs(df1))
            ab(l,k,j,i) = df1*product(basis%bfnrm([i1,j1,k1,l1]))
          end do
        end do
      end do
    end do
  end subroutine grd2_mrsf_compute_data_t_get_density

!###############################################################################

end module tdhf_mrsf_gradient_mod
