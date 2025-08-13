module tdhf_z_vector_mod
  use types, only: information
  use precision, only: dp
  use basis_tools, only: basis_set
  use int2_compute, only: int2_compute_t
  use io_constants, only: iw
  use tdhf_lib, only: int2_td_data_t, &
    int2_fock_data_t, int2_tdgrd_data_t
  use mod_dft_molgrid, only: dft_grid_t
  use oqp_linalg

  implicit none

  character(len=*), parameter :: module_name = "tdhf_z_vector_mod"

  private
  public tdhf_z_vector_C
  public oqp_tdhf_z_vector

  type :: tdhf_cg_data
    type(information), pointer :: infos
    type(int2_compute_t), pointer :: int2_driver
    class(int2_fock_data_t), pointer :: int2_data
    type(dft_grid_t), pointer :: molGrid

    real(kind=dp), pointer :: wrk(:,:)
    real(kind=dp), pointer :: mo(:,:)
    real(kind=dp), pointer :: pa(:,:,:)
    real(kind=dp), pointer :: xm(:)
    real(kind=dp), pointer :: xminv(:)

    integer :: nbf
    integer :: nocc
    logical :: dft

  end type

contains

!###############################################################################

  subroutine tdhf_z_vector_C(c_handle) bind(C, name="tdhf_z_vector")
    use types, only: information
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use strings, only: Cstring
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call oqp_tdhf_z_vector(inf)
  end subroutine

!###############################################################################

  subroutine oqp_tdhf_z_vector(infos)
    use, intrinsic :: iso_c_binding, only: c_int32_t
    use oqp_tagarray_driver

    use strings, only: Cstring, fstring
    use messages, only: show_message, with_abort
    use util, only: measure_time

    use tdhf_lib, only: iatogen, mntoia, esum, &
      tdhf_unrelaxed_density
    use tdhf_sf_lib, only: sfrorhs, &
      sfromcal, sfrogen, sfrolhs, &
      pcgb, sfropcal, sfrowcal
    use dft, only: dft_initialize, dftclean
    use mod_dft_gridint_fxc, only: tddft_fxc
    use mod_dft_gridint_gxc, only: tddft_gxc
    use mathlib, only: symmetrize_matrix, orthogonal_transform
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: pack_matrix, unpack_matrix
    use mathlib, only: triangular_to_full
    use tdhf_lib, only: int2_rpagrd_data_t
    use pcg_mod
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "oqp_tdhf_z_vector"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos


    type(dft_grid_t), target :: molGrid
    type(int2_compute_t), target :: int2_driver
    class(int2_fock_data_t), allocatable, target :: int2_data
    type(tdhf_cg_data) :: cgdata
    type(pcg_t) :: pcg

    real(kind=dp), allocatable :: hpp(:,:,:), hpt(:,:,:), hmm(:,:,:), gxp(:,:,:)
    real(kind=dp), allocatable, target :: wrk1(:,:), &
            rhs(:), xm(:), zvec(:), xminv(:), pa(:,:,:)

    real(kind=dp), contiguous, pointer :: ppa(:,:,:,:), &
            pxm(:,:), prhs(:,:), wmo(:,:)


    logical :: dft, tda
    integer :: scf_type, mol_mult
    integer :: i, j, iter, ok
    integer :: nocc, nvir, lexc, nbf, nbf2
    real(kind=dp) :: cnvtol, scale_exch

    ! tagarray
    integer(c_int32_t) :: stat
    real(kind=dp), contiguous, pointer :: &
      mo_a(:,:), mo_energy_a(:), wao(:), td_p(:,:), td_t(:,:), &
      ta(:), xpy(:,:), xmy(:,:), td_energies(:)

    character(len=*), parameter :: &
      tags_required(*) = [  character(len=80) :: &
        OQP_E_MO_A, OQP_VEC_MO_A, OQP_TD_T, OQP_TD_XPY, OQP_TD_XMY, &
        OQP_td_energies &
      ]

  ! Log file
    open (unit=IW, file=infos%log_filename, position="append")
  !
    call print_module_info('TDHF_Z_Vector','Solving Z-Vector for TDDFT')
    mol_mult = infos%mol_prop%mult
    scf_type = infos%control%scftype
    dft = infos%control%hamilton == 20
    tda = infos%tddft%tda

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2

    if (dft) call dft_initialize(infos, basis, molGrid)

  ! Parameter it should be inputed later
  ! convergence tolerance in the iterative TD-DFT step.
    cnvtol = infos%tddft%zvconv

    nocc = infos%mol_prop%nocc
    nvir = nbf-nocc
    lexc = nocc*nvir

    allocate(&
      rhs(lexc), &
      pa(nbf,nbf,1), &
      wrk1(nbf,nbf), &
      stat=ok, &
      source=0.0_dp)

    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    call infos%dat%erase([character(80) :: OQP_WAO, OQP_TD_P])

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_TD_T, td_t)
    call tagarray_get_data(infos%dat, OQP_TD_XPY, xpy)
    call tagarray_get_data(infos%dat, OQP_TD_XMY, xmy)
    call tagarray_get_data(infos%dat, OQP_td_energies, td_energies)
    ta => td_t(:,1)

    write(iw,'(/1x,71("-")&
             &/19x,"TD-DFT ENERGY GRADIENT CALCULATION"&
             &/1x,71("-")/)')

    write(iw,fmt='(5x,a/&
                  &5x,16("-")/&
                  &5x,a,x,i0,x,f17.10,x,"Hartree"/&
                  &5x,a,x,i0/&
                  &5x,a,x,e10.4/&
                  &5x,a,x,i0)') &
        'Z-vector options' &
      , 'Target state       is', infos%tddft%target_state, infos%mol_energy%energy+td_energies(infos%tddft%target_state) &
      , 'Multiplicity       is', infos%tddft%mult &
      , 'Convergence        is', infos%tddft%zvconv &
      , 'Maximum iterations is', infos%control%maxit_zv
    call flush(iw)

!   1. Compute right-hand side of Z-vector equation
    allocate(hpp(nbf,nbf,1), & ! H+[X+Y]
             hpt(nbf,nbf,1), & ! H+[T]
             hmm(nbf,nbf,1), & ! H-[X-Y]
             gxp(nbf,nbf,1), & ! g_xc[X+Y][X+Y]
             source=0.0d0)

    call tdhf_unrelaxed_density(xmy(:,infos%tddft%target_state), xpy(:,infos%tddft%target_state), mo_a, td_t(:,1), nocc, tda)

    ! Initialize ERI calculations
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    ! Compute H+[X+Y], H+[T], H-[X-Y], and G_xc
    call compute_r_terms(infos, basis, int2_driver, molGrid, mo_a, xpy(:,infos%tddft%target_state), &
                         xmy(:,infos%tddft%target_state), ta, hpp, hpt, hmm, gxp)

    ! Transform Gxc to MO basis
    call orthogonal_transform('n', nbf, mo_a, gxp(:,:,1))

    ! Assemble RHS in MO basis
    prhs(1:nocc,1:nvir) => rhs(1:)
    call compute_r_mo(prhs, mo_a, xpy(:,infos%tddft%target_state), xmy(:,infos%tddft%target_state), hpt, hpp, hmm, gxp(:,:,1))
    rhs = -rhs

!   2. Initialize CG for Z-vector solution
    write(iw,'(/3x,25("-")&
             &/6x,"START Z-VECTOR LOOP"&
             &/3x,25("-")/)')
    call flush(iw)

    allocate(xminv(lexc), &
             xm(lexc), &
             stat=ok, &
             source=0.0_dp)

    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    scale_exch = 1.0_dp
    if (dft) scale_exch = infos%dft%HFscale
    int2_data = int2_td_data_t(d2=pa, &
            int_apb=.true., &
            int_amb=.false., &
            tamm_dancoff=tda, &
            scale_exchange=scale_exch)

    pxm(1:nocc, 1:nvir) => xm(1:)
    do i = 1, nvir
      do j = 1, nocc
        pxm(j,i) = mo_energy_a(nocc+i) - mo_energy_a(j)
      end do
    end do
    xminv = 1.0d0/xm

    cgdata = tdhf_cg_data( &
        infos=infos, int2_driver=int2_driver, &
        int2_data=int2_data, molgrid=molgrid, &
        wrk=wrk1, mo=mo_a, pa=pa, xm=xm, xminv=xminv, &
        nbf=nbf, nocc = nocc, dft = dft &
      )

    call pcg%init(b=rhs, &
      update=compute_apbx, precond=precond, &
      dat=cgdata, tol=sqrt(abs(cnvtol)))

    write(iw,'(" INITIAL ERROR =",3X,1P,E10.3,1X,"/",1P,E10.3)') pcg%error**2, cnvtol

    ! Begin CG iterations
    do iter = 1, infos%control%maxit_zv
      if (pcg%errcode /= PCG_OK) exit

      call pcg%step()

      write(iw,'(" ITER#",I2," ERROR =",3X,1P,E10.3,1X,"/",1P,E10.3)') &
              iter, pcg%error**2, cnvtol
      call flush(iw)
    end do

    select case (pcg%errcode)
    case (PCG_CONVERGED)
      write(iw,'(/3x,24("-")&
               &/6x,"Z-Vector converged"&
               &/3x,24("-")/)')
      infos%mol_energy%Z_Vector_converged=.true.
    case default
      write(iw,'(/3x,24("-")&
               &/6x,"Z-Vector not converged"&
               &/3x,24("-")/)')
      infos%mol_energy%Z_Vector_converged=.false.
    end select
    call flush(iw)

!   Save Z-vector and clean up memory
    deallocate(xminv, rhs, xm)
    call int2_data%clean()
    deallocate(int2_data)
    allocate(zvec, source=pcg%x)
    call pcg%clean()

!   3. Now, compute relaxed energy-weighted difference density matrix W

    ! Convert Z-vector from MO to AO and assemble P = T+Z
    pa = 0
    call iatogen(zvec, pa(:,:,1), nocc, nocc)
    call symmetrize_matrix(pa(:,:,1), nbf)
    pa(:,:,1) = 0.5*pa(:,:,1)
    call orthogonal_transform('t', nbf, mo_a, pa(:,:,1))

    call unpack_matrix(ta, wrk1) ! T
    pa(:,:,1) = pa(:,:,1) + wrk1 ! T+Z

    ! Store relaxed difference density matrix P to global memory
    stat = infos%dat%create(OQP_td_p, TA_TYPE_REAL64, [nbf2, 1], description=OQP_td_p_comment)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)
    call pack_matrix(pa(:,:,1), td_p(:,1))

    ! Compute H+[P]
    ppa(1:nbf,1:nbf, 1:1, 1:1) => pa
    int2_data = int2_rpagrd_data_t(&
            xpy=null(), &
            xmy=null(), &
            t=ppa, &
            nspin = 1, &
            tamm_dancoff=tda, &
            scale_exchange=scale_exch)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
    select type (int2_data)
    type is (int2_rpagrd_data_t)
      hpt = int2_data%hpt(:,:,:,1,1)
    end select

    call int2_data%clean()
    deallocate(int2_data)

    if (dft) then
      pa = pa*2
      call tddft_fxc(basis=basis, &
             molGrid=molGrid, &
             isVecs=.true., &
             wf=MO_A, &
             fx=hpt(:,:,1:1), &
             dx=pa(:,:,1:1), &
             nmtx=1, &
             !threshold=1.0d-15, &
             threshold=0.0d0, &
             infos=infos)
    end if

    ! Transform H+[P], H+[X+Y], H-[X-Y] from AO to MO basis
    call orthogonal_transform('n', nbf, mo_a, hpt(:,:,1), wrk=wrk1) ! P
    call orthogonal_transform('n', nbf, mo_a, hpp(:,:,1), wrk=wrk1) ! X+Y
    call orthogonal_transform('n', nbf, mo_a, hmm(:,:,1), wrk=wrk1) ! X-Y

    ! Calculate W in MO basis
    wmo(1:nbf,1:nbf) => wrk1
    call compute_w_mo(wmo, xpy(:,infos%tddft%target_state), xmy(:,infos%tddft%target_state), zvec, &
            hpt(:,:,1), hpp(:,:,1), hmm(:,:,1), &
            mo_energy_a, td_energies(infos%tddft%target_state), nocc, gxp(:,:,1))

    ! Transform W from MO to AO basis
    call orthogonal_transform('t', nbf, mo_a, wmo)

    ! Store W to global memory
    stat = infos%dat%create(OQP_WAO, TA_TYPE_REAL64, (/ nbf2 /), description=OQP_WAO_comment)
    call tagarray_get_data(infos%dat, OQP_WAO, wao)
    call pack_matrix(wmo, wao)
    wao = wao*0.5_dp

    ! Cleanup
    call int2_driver%clean()
    if (dft) call dftclean(infos)
    call measure_time(print_total=1, log_unit=iw)
    close(iw)

  end subroutine oqp_tdhf_z_vector

!###############################################################################

!> @brief Compute H+[X+Y], H+[T], and H-[X-Y]
  subroutine compute_r_terms(infos, basis, int2_driver, molGrid, mo_a, xpy, xmy, ta, hpp, hpt, hmm, gxp)

    use mod_dft_gridint_fxc, only: tddft_fxc
    use mod_dft_gridint_gxc, only: tddft_gxc
    use tdhf_lib, only: int2_rpagrd_data_t
    use mathlib, only: symmetrize_matrix, orthogonal_transform
    use mathlib, only: unpack_matrix
    use tdhf_lib, only: iatogen

    type(information) :: infos
    type(basis_set) :: basis
    type(dft_grid_t) :: molGrid
    type(int2_compute_t), target :: int2_driver
    real(kind=dp), intent(in) :: mo_a(:,:), xpy(:), xmy(:), ta(:)
    real(kind=dp), target :: hpp(:,:,:), hpt(:,:,:), hmm(:,:,:), gxp(:,:,:)
    real(kind=dp), allocatable, target :: xpy2(:,:,:,:), xmy2(:,:,:,:)
    real(kind=dp), allocatable, target :: hp(:,:,:)
    real(kind=dp) :: scale_exch
    integer :: nocc, nvir, nbf
    logical :: tda, dft
    type(int2_rpagrd_data_t), target :: int2_data

    nbf = basis%nbf
    nocc = infos%mol_prop%nocc
    nvir = nbf-nocc
    tda = infos%tddft%tda
    dft = infos%control%hamilton == 20

    allocate(xpy2(nbf,nbf,2,1), &
             xmy2(nbf,nbf,1,1), &
             source=0.0_dp)

    ! Prepare unrelaxed density, X+Y, and X-Y vectors in AO basis set
    call iatogen(xpy, xpy2(:,:,1,1), nocc, nocc)
    call symmetrize_matrix(xpy2(:,:,1,1), nbf)
    xpy2(:,:,1,1) = xpy2(:,:,1,1)*0.5
    call orthogonal_transform('t', nbf, mo_a, xpy2(:,:,1,1))

    call iatogen(xmy, xmy2(:,:,1,1), nocc, nocc)
    call orthogonal_transform('t', nbf, mo_a, xmy2(:,:,1,1))

    call unpack_matrix(ta, xpy2(:,:,2,1))

    ! Initialize ERI calculations
    scale_exch = 1.0_dp
    if (dft) scale_exch = infos%dft%HFscale

    if (dft) then
      call tddft_gxc(basis=basis, &
             molGrid=molGrid, &
             isVecs=.true., &
             wf=MO_A, &
             fx=gxp(:,:,1:1), &
             dx=xpy2(:,:,1:1,1), &
             nmtx=1, &
             !threshold=1.0d-15, &
             threshold=0.0d0, &
             infos=infos)

      allocate(hp(nbf,nbf,2), source=0.0_dp)
      call tddft_fxc(basis=basis, &
             molGrid=molGrid, &
             isVecs=.true., &
             wf=MO_A, &
             fx=hp(:,:,1:2), &
             dx=xpy2(:,:,1:2,1), &
             nmtx=2, &
             !threshold=1.0d-15, &
             threshold=0.0d0, &
             infos=infos)

      hpp(:,:,1) = hpp(:,:,1) + hp(:,:,1)
      hpt(:,:,1) = hpt(:,:,1) + hp(:,:,2)
      deallocate(hp)

    end if

    ! Compute H+[X+Y], H+[T], and H-[X-Y]
    int2_data = int2_rpagrd_data_t(&
            xpy=xpy2, &
            xmy=xmy2, &
            t=xpy2(:,:,2:2,1:1), &
            nspin = 1, &
            tamm_dancoff=tda, &
            scale_exchange=scale_exch)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)

      hpp(:,:,:) = 2*hpp(:,:,:) + int2_data%hpp(:,:,:,1,1) ! H+[X+Y]
      hpt(:,:,:) = 2*hpt(:,:,:) + int2_data%hpt(:,:,:,1,1) ! H+[T]
      hmm(:,:,:) = int2_data%hmm(:,:,:,1,1) ! H-[X-Y]

    call int2_data%clean()

  end subroutine

!###############################################################################

  subroutine compute_w_mo(wmo, xpy, xmy, z, hpt, hpp, hmm, e_mo, e_state, nocc, gxp)

    use mathlib, only: triangular_to_full, orthogonal_transform
    implicit none

    real(kind=dp), intent(in) :: xpy(nocc,*), xmy(nocc,*), z(nocc,*), e_mo(*), hpt(:,:), hpp(:,:), hmm(:,:)
    real(kind=dp), optional, intent(in) :: gxp(:,:)
    real(kind=dp), intent(in) :: e_state
    real(kind=dp), intent(out) :: wmo(:,:)
    integer, intent(in) :: nocc

    real(kind=dp), allocatable, target :: wrk1(:,:), wrk_ia(:,:)

    integer :: nbf, nvir, i

    nbf = ubound(wmo,1)
    nvir = nbf-nocc

    allocate(wrk1(nbf,nbf), &
             wrk_ia(nocc,nvir), &
             source=0.0d0)

    wmo = 0

    ! W_ij (occ x occ)
    call dsyr2k('u','n',nocc,nvir, &
            e_state, xpy, nocc, &
                     xmy, nocc, &
            1.0d0,   wmo, nbf)

    do i = 1, nvir
      wrk_ia(:,i) = xpy(:,i) * e_mo(nocc+i)
    end do
    call dgemm('n','t', nocc, nocc, nvir, &
            -1.0d0, xpy,    nocc, &
                    wrk_ia, nocc, &
             1.0d0, wmo, nbf)

    do i = 1, nvir
      wrk_ia(:,i) = xmy(:,i) * e_mo(nocc+i)
    end do
    call dgemm('n', 't', nocc, nocc, nvir, &
           -1.0d0, xmy,     nocc, &
                   wrk_ia,  nocc, &
            1.0d0, wmo, nbf)

    wmo(:nocc,:nocc) = wmo(:nocc,:nocc) + hpt(:nocc,:nocc)

    if (present(gxp)) then
      wmo(:nocc,:nocc) = wmo(:nocc,:nocc) + 2*gxp(:nocc,:nocc)
    end if

    ! W_ab (vir x vir)
    call dsyr2k('u','t',nvir,nocc, &
            e_state, xpy, nocc, &
                     xmy, nocc, &
            1.0d0,   wmo(nocc+1:,nocc+1), nbf)

    do i = 1, nocc
      wrk_ia(i,:) = xpy(i,:nvir) * e_mo(i)
    end do
    call dgemm('t', 'n', nvir, nvir, nocc, &
            1.0d0, xpy,    nocc, &
                   wrk_ia, nocc, &
            1.0d0, wmo(nocc+1:,nocc+1), nbf)

    do i = 1, nocc
      wrk_ia(i,:) = xmy(i,:nvir) * e_mo(i)
    end do
    call dgemm('t', 'n', nvir, nvir, nocc, &
            1.0d0, xmy,    nocc, &
                   wrk_ia, nocc, &
            1.0d0, wmo(nocc+1:,nocc+1), nbf)

    ! W_ia (occ x vir)
    call dgemm('t', 'n', nocc, nvir, nocc,  &
               1.0_dp, hpp, nbf, &
                       xpy, nocc, &
               1.0_dp, wmo(:,nocc+1), nbf)

    call dgemm('t', 'n', nocc, nvir, nocc,  &
               1.0_dp, hmm, nbf, &
                       xmy, nocc, &
               1.0_dp, wmo(:,nocc+1), nbf)

    do i = 1, nocc
      wmo(i,nocc+1:) = wmo(i,nocc+1:) + z(i,:nvir)*e_mo(i)
    end do

    call triangular_to_full(wmo,nbf,'u')

  end subroutine

!###############################################################################

  subroutine compute_r_mo(rhs, mo, xpy, xmy, hpt, hpp, hmm, gxp)

    use mathlib, only: orthogonal_transform
    implicit none

    real(kind=dp), contiguous, target, intent(out) :: rhs(:,:)
    real(kind=dp), intent(in) :: xpy(*), xmy(*)
    real(kind=dp), contiguous, intent(in) :: mo(:,:)
    real(kind=dp), contiguous, intent(inout) :: hpt(:,:,:), hpp(:,:,:), hmm(:,:,:), gxp(:,:)

    real(kind=dp), allocatable :: wrk1(:,:)

    integer :: nbf, nocc, nvir

    nbf = ubound(hpt,1)
    nocc = ubound(rhs,1)
    nvir = ubound(rhs,2)

    allocate(wrk1(nbf,nbf))

    rhs = 0
    call orthogonal_transform('n', nbf, mo, hpp(:,:,1), wrk1)
    call dgemm('n', 't', nocc, nvir, nvir,  &
               1.0_dp, xpy,  nocc, &
                       wrk1(nocc+1,nocc+1), nbf, &
               1.0_dp, rhs, nocc)
    call dgemm('t', 'n', nocc, nvir, nocc,  &
              -1.0_dp, wrk1, nbf, &
                       xpy,  nocc, &
               1.0_dp, rhs, nocc)

    call orthogonal_transform('n', nbf, mo, hmm(:,:,1), wrk1)
    call dgemm('n', 't', nocc, nvir, nvir,  &
               1.0_dp, xmy,  nocc, &
                       wrk1(nocc+1,nocc+1), nbf, &
               1.0_dp, rhs, nocc)
    call dgemm('t', 'n', nocc, nvir, nocc,  &
              -1.0_dp, wrk1, nbf, &
                       xmy,  nocc, &
               1.0_dp, rhs, nocc)

    call orthogonal_transform('n', nbf, mo, hpt(:,:,1), wrk1)

    rhs = rhs + wrk1(1:nocc,nocc+1:)

    rhs = rhs + 2*gxp(1:nocc,nocc+1:)

  end subroutine

!###############################################################################

  subroutine compute_apbx(y, x, dat)

    use iso_c_binding, only: c_ptr, c_f_pointer
    use tdhf_lib, only: iatogen, mntoia
    use mathlib, only: symmetrize_matrix, orthogonal_transform
    use mod_dft_gridint_fxc, only: tddft_fxc
    implicit none

    real(kind=dp) :: x(:)
    real(kind=dp) :: y(:)
    type(c_ptr) :: dat
    type(tdhf_cg_data), pointer :: p
    real(kind=dp), pointer :: apb(:,:,:)

    call c_f_pointer(dat, p)

    associate ( wrk => p%wrk, nocc => p%nocc, nbf  => p%nbf &
              , mo  => p%mo, pa  => p%pa &
              , int2_driver => p%int2_driver &
              , int2_data => p%int2_data &
              , infos => p%infos &
              , molgrid => p%molgrid &
              , dft => p%dft &
              , xm => p%xm &
              )

      call iatogen(x, wrk, nocc, nocc)
      call symmetrize_matrix(wrk, nbf)
      call orthogonal_transform('t', nbf, mo, wrk, pa(:,:,1))

!     (A+B)*PK
      call int2_driver%run(int2_data, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%dft%cam_alpha, &
              beta=infos%dft%cam_beta,&
              mu=infos%dft%cam_mu)

      select type (int2_data)
      type is (int2_td_data_t)
        apb => int2_data%apb(:,:,:,1)
      end select
      apb = apb * 0.5

      if (dft) then
        call tddft_fxc(basis=infos%basis, &
               molGrid=molGrid, &
               isVecs=.true., &
               wf=mo, &
               fx=apb(:,:,1:1), &
               dx=pa(:,:,1:1), &
               nmtx=1, &
               !threshold=1.0d-15, &
               threshold=0.0d0, &
               infos=infos)
      end if

      call mntoia(apb(:,:,1), y, mo, mo, nocc, nocc)

      y = y + xm*x

    end associate

  end subroutine

!###############################################################################

  subroutine precond(y, x, dat)

    use iso_c_binding, only: c_ptr, c_f_pointer
    implicit none

    real(kind=dp) :: x(:)
    real(kind=dp) :: y(:)
    type(c_ptr) :: dat
    type(tdhf_cg_data), pointer :: p

    call c_f_pointer(dat, p)

    ! diagonal approximation:
    ! (A+B)^{-1} ~ 1/(e_vir-e_occ)
    y = p%xminv*x

  end subroutine
end module
