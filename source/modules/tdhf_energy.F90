module tdhf_energy_mod
  implicit none

  character(len=*), parameter :: module_name = "tdhf_energy_mod"

contains

  subroutine tdhf_energy_C(c_handle) bind(C, name="tdhf_energy")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_energy(inf)
  end subroutine tdhf_energy_C

  subroutine tdhf_energy(infos)
    use, intrinsic :: iso_c_binding, only: c_int32_t
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use strings, only: Cstring, fstring
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort

    use precision, only: dp
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_td_data_t
    use tdhf_lib, only: &
      inivec, iatogen, mntoia, rparedms, rpaeig, rpavnorm, &
      rpaechk, rpaprint, rparesvec, rpaexpndv, rpanewb, esum, &
      tdhf_unrelaxed_density
    use dft, only: dft_initialize, dftclean
    use mod_dft_gridint_fxc, only: tddft_fxc
    use util, only: measure_time
    use mathlib, only: symmetrize_matrix, orthogonal_transform
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: unpack_matrix
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_energy"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos

    integer :: ok

    real(kind=dp), allocatable :: scr2(:)
    real(kind=dp), allocatable :: ab1_mo(:,:), ab2_mo(:,:), abxc(:,:), scr3(:,:)
    real(kind=dp), allocatable :: eex(:)
    real(kind=dp), allocatable :: abm_l(:,:), abp_r(:,:), amb(:,:), &
                                  apb(:,:), vlo(:,:), vro(:,:)
    real(kind=dp), allocatable, target :: vl(:), vr(:)
    real(kind=dp), pointer :: vl_p(:,:), vr_p(:,:)
    real(kind=dp), allocatable :: xm(:)
    real(kind=dp), allocatable :: bvec_mo(:,:)
    real(kind=dp), allocatable, target :: bvec(:,:,:)
    real(kind=dp), allocatable :: scr1(:,:)
    real(kind=dp), allocatable :: errors(:)
    real(kind=dp), pointer :: ab2(:,:,:)
    real(kind=dp), pointer :: ab1(:,:,:)
    real(kind=dp), allocatable :: dip(:,:)

    integer :: nocc, nvir
    integer :: nbf, nbf2, lexc
    integer :: ndsr, mxvec, nmax, ist, iend, nvec, novec
    integer :: iter, istart, nv, iv, ivec
    integer :: mxiter
    logical :: do_apb,do_amb, tamm_dancoff   ! SWITCH FOR a-b, a+b, tANDAMM COFF
    integer :: imax
    logical :: converged
    integer :: ierr
    real(kind=dp) :: mxerr, cnvtol, scale_exch, norm
    integer :: maxvec, nstates

    type(int2_compute_t) :: int2_driver
    type(int2_td_data_t), target :: int2_data
    type(dft_grid_t) :: molGrid

    logical :: dft

    ! tagarray
    integer(c_int32_t) :: stat
    real(kind=dp), contiguous, pointer :: &
      mo_energy_a(:), mo_a(:,:), td_t(:,:), &
      xpy(:,:), xmy(:,:), td_energies(:)
    character(len=*), parameter :: tags_alloc(4) = (/ character(len=80) :: &
      OQP_td_t, OQP_td_xpy, OQP_td_xmy, OQP_td_energies /)
    character(len=*), parameter :: tags_alpha(2) = (/ character(len=80) :: &
      OQP_E_MO_A, OQP_VEC_MO_A /)

    dft = infos%control%hamilton == 20

  ! Files open
  ! 3. LOG: Write: Main output file
    open (unit=IW, file=infos%log_filename, position="append")
  !
    call print_module_info('THDF_Energy','Computing Energy of TDDFT')

  ! Readings

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

  ! Allocate H, S ,T and D matrices
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2

    call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

    if (dft) call dft_initialize(infos, basis, molGrid, verbose=.true.)

    nstates = infos%tddft%nstate
    maxvec = infos%tddft%maxvec
    cnvtol = infos%tddft%cnvtol

    do_apb = .true.  ! Always do A+B part
    do_amb = .true.  ! = isGGA ! Do A-B part if a hybrid functional.
    tamm_dancoff = infos%tddft%tda  ! Tamm/Dancoff approximation

    nocc = infos%mol_prop%nocc
    nvir = nbf-nocc
    lexc = nocc*nvir

    mxvec = min(maxvec*nstates,lexc)
    nvec = min(nstates,mxvec)
    ndsr = nvec
    nmax = 2*nvec
    !nvec = 2*nvec

  ! Allocate TDDFT variables
    allocate(xm(lexc), &
             abxc(nbf,nbf), &
             bvec_mo(lexc,mxvec), &
             ab1_mo(lexc,mxvec), &
             ab2_mo(lexc,mxvec), &
             bvec(nbf,nbf,nmax), &
             eex(mxvec), &
             errors(mxvec), &
             apb(mxvec,mxvec), &
             amb(mxvec,mxvec), &
             vr(mxvec*mxvec), &
             vro(lexc,mxvec), &
             vl(mxvec*mxvec), &
             vlo(lexc,mxvec), &
             abp_r(lexc,mxvec), &
             abm_l(lexc,mxvec), &
             dip(3,nstates), &
             scr1(nbf,nbf), &
             scr2(mxvec*mxvec), &
             scr3(lexc,nmax), &
             source=0.0d0, &
             stat=ok)
    if( ok/=0 ) &
      call show_message('Cannot allocate memory', with_abort)


  ! Initialize ERI calculations
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    scale_exch = 1.0_dp
    if (dft) scale_exch = infos%dft%HFscale

  ! Construct TD trial vector
    call inivec(mo_energy_a,mo_energy_a,bvec_mo,xm,nocc,nocc,nvec)

    ist = 1
    istart = 1
    iend = nvec
    iter = 0
    mxiter = infos%control%maxit_dav
    ierr = 0

    do iter = 1, mxiter
      nv = iend-ist+1

!$omp parallel private(iv,scr1,abxc,ivec)
!$omp do schedule(static)
      do ivec = ist, iend
        iv = ivec-ist+1
        call iatogen(bvec_mo(:,ivec),abxc,nocc,nocc)
        call orthogonal_transform('t', nbf, mo_a, abxc, bvec(:,:,iv), scr1)
      end do
!$omp end do
!$omp end parallel

      int2_data = int2_td_data_t(d2=bvec(:,:,:nv), &
              int_apb=do_apb, int_amb=do_amb, tamm_dancoff=tamm_dancoff, &
              tamm_dancoff_coulomb=.true., &
              scale_exchange=scale_exch)

      call int2_driver%run(int2_data, &
                cam=dft.and.infos%dft%cam_flag, &
                alpha=infos%dft%cam_alpha, &
                beta=infos%dft%cam_beta,&
                mu=infos%dft%cam_mu)
      ab1 => int2_data%apb(:,:,:,1)
      ab2 => int2_data%amb(:,:,:,1)

      if (dft) then
        do ivec = ist, iend
          iv = ivec-ist+1
          call symmetrize_matrix(bvec(:,:,iv), nbf)
        end do
        call tddft_fxc(basis=basis, &
               molGrid=molGrid, &
               isVecs=.true., &
               wf=mo_a, &
               fx=ab1(:,:,:iv), &
               dx=bvec(:,:,:iv), &
               nmtx=iv, &
               !threshold=1.0d-15, &
               threshold=0.0d0, &
               infos=infos)
      end if

!$omp parallel private(iv,ivec)
!$omp do schedule(static)
      do ivec = ist, iend
        iv = ivec-ist+1
  !     Convert AO to MO basis
        if (tamm_dancoff) then
          ab1(:,:,iv) = 0.5_dp*ab1(:,:,iv) + ab2(:,:,iv)
        end if
        call mntoia(ab1(:,:,iv), ab1_mo(:,ivec), mo_a, mo_a, nocc, nocc)

  !     Product (A+B)*X
        call esum(mo_energy_a,ab1_mo,bvec_mo,nocc,ivec)

        if (.not.tamm_dancoff) then
!         Product (A-B)*X
          call mntoia(ab2(:,:,iv),ab2_mo(:,ivec),mo_a,mo_a,nocc,nocc)
          call esum(mo_energy_a,ab2_mo,bvec_mo,nocc,ivec)
        end if
      end do
!$omp end do
!$omp end parallel

      call rparedms(bvec_mo,ab1_mo,ab2_mo,apb,amb,nvec,tamm_dancoff)

      vl_p(1:nvec, 1:nvec) => vl(1:nvec*nvec)
      vr_p(1:nvec, 1:nvec) => vr(1:nvec*nvec)
      call rpaeig(eex,vl_p,vr_p,apb,amb,scr2,tamm_dancoff)
      call rpavnorm(vr_p,vl_p,tamm_dancoff)
      call rpaexpndv(vr_p,vl_p,vro,vlo,bvec_mo,bvec_mo,ndsr,tamm_dancoff)
      call rpaexpndv(vr_p,vl_p,abp_r,abm_l,ab1_mo,ab2_mo,ndsr,tamm_dancoff)

      call rpaechk(eex,nvec,ndsr,imax,tamm_dancoff)
!     Residual vectors W
!     Get perturbed vectors Q if required
      call rparesvec(scr3,abp_r,abm_l,vlo,vro,eex,xm,ndsr,errors,cnvtol,imax,tamm_dancoff)

      mxerr = maxval(errors(imax+1:ndsr))

      call rpaprint(eex, errors, cnvtol, iter, imax, ndsr)

!     Check convergence
      converged = mxerr<=cnvtol
      if (converged) exit

!     No space left for new vectors, exit
      if (nvec==mxvec) ierr = 1
      if (ierr/=0) exit

      call rpanewb(ndsr,bvec_mo,scr3,novec,nvec,ierr,tamm_dancoff)

  !   ierr=1 nvec over mxvec: not converged case
      if (ierr/=0) exit

      ist = novec+1
      iend = nvec
    end do

    if (iter >= mxiter .and. .not. converged) ierr = -1

    select case (ierr)
    case (-1)
      write(*,'(/,2X,"TD-DFT energies NOT CONVERGED after ",I4," iterations"/)') mxiter
      infos%mol_energy%Davidson_converged=.false.
!      call show_message("Aborting. Try to increase maxit or check your system.", WITH_ABORT)
    case (0)
      write(*,'(/,2X,"TD-DFT energies converged in ",I4," iterations"/)') iter
      infos%mol_energy%Davidson_converged=.true.
    case (1)
      write(*,'(/,2X,"..something is wrong.. nvec = mxvec")')
      infos%mol_energy%Davidson_converged=.false.
    case (2)
      write(*,'(/,2x,"..something is wrong..  nvec > mxvec")')
      write(*,'(3x,"nvec/mxvec =",I4,"/",I4)') nvec, mxvec
      infos%mol_energy%Davidson_converged=.false.
    case (3)
      write(*,'(/,2x,"..something is wrong.. No vectors were added")')
      infos%mol_energy%Davidson_converged=.false.
    end select

    call get_td_transition_dipole(basis, dip, mo_a, vro, vlo, nstates, nocc)

    call td_print_results(eex, dip, nstates)

    call measure_time(print_total=1, log_unit=iw)

    ! Bi-orthogonalize X+Y and X-Y:
    ! (X+Y)_i \dot (X-Y)_j = \delta_ij
    ! only normalization is needed
    ! RHF reference case: alpha == beta, additional factor of 1/sqrt(2)
    do ist = 1, nstates
        norm = 1/sqrt(2*dot_product(vlo(:,ist), vro(:,ist)))
        vlo(:,ist) = vlo(:,ist) * norm
        vro(:,ist) = vro(:,ist) * norm
    end do

    call infos%dat%erase(tags_alloc)

    stat = infos%dat%create(OQP_td_t, &
            TA_TYPE_REAL64, &
            [nbf2, 1], &
            description=OQP_td_t_comment)

    stat = infos%dat%create(OQP_td_xpy, &
            TA_TYPE_REAL64, &
            [lexc, nstates], &
            description="(X+Y) vector for target state in TD-DFT calculations")

    stat = infos%dat%create(OQP_td_xmy, &
            TA_TYPE_REAL64, &
            [lexc, nstates], &
            description="(X-Y) vector for target state in TD-DFT calculations")

    stat = infos%dat%create(OQP_td_energies, &
            TA_TYPE_REAL64, &
            [nstates], &
            description=OQP_td_energies_comment)

    call tagarray_get_data(infos%dat, OQP_td_t, td_t)
    call tagarray_get_data(infos%dat, OQP_td_xpy, xpy)
    call tagarray_get_data(infos%dat, OQP_td_xmy, xmy)
    xpy = vro(:,:nstates)
    xmy = vlo(:,:nstates)

    call tagarray_get_data(infos%dat, OQP_td_energies, td_energies)
    td_energies = eex(:nstates)
    infos%mol_energy%excited_energy = td_energies(infos%tddft%target_state)

    call int2_driver%clean()

    if (dft) call dftclean(infos)

    close(iw)

  end subroutine tdhf_energy

  subroutine get_td_transition_dipole(basis, dip, v, vr, vl, nstates, nocc)
    use int1
    use types, only: information
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use tdhf_lib, only: iatogen
    use mathlib, only: orthogonal_transform, symmetrize_matrix, traceprod_sym_packed
    use mathlib, only: pack_matrix

    implicit none

    type(basis_set), intent(in) :: basis
    real(kind=8) :: vl(:,:), vr(:,:), v(:,:)
    real(kind=8) :: dip(:,:)
    integer :: nstates, nocc

    real(kind=8) :: com(3)
    real(kind=8), allocatable :: mints(:,:), trden(:,:), trden_ao(:,:)
    real(kind=8), allocatable, target :: tmp(:,:)
    real(kind=8), pointer :: tmp2(:)
    integer :: nbf, nbf2, ok
    integer :: i, j, s1

    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2

    allocate(mints(nbf2,3), &
             trden(nbf,nbf), &
             trden_ao(nbf,nbf), &
             tmp(nbf,nbf), &
             source=0.0d0, stat=ok)

    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    ! Compute dipole integrals at the center of mass
    com = basis%atoms%center(weight='mass')

    call multipole_integrals(basis, mints, com, 1)

    ! Compute transition dipole between the ground state and all excited
    do s1 = 1, nstates
        ! Compute transition density
        ! unpack X+Y
        call iatogen(vr(:,s1), trden, nocc, nocc)
        ! unpack X-Y
        call iatogen(vl(:,s1), tmp, nocc, nocc)
        do i = nocc+1, nbf
          do j = 1, nocc
            trden(i,j) = trden(j,i) - tmp(j,i) ! Y: vir -> occ
            trden(j,i) = trden(j,i) + tmp(j,i) ! X: occ -> vir
          end do
        end do

        ! Convert transition density from MO to AO basis
        call orthogonal_transform('t', nbf, v, trden, trden_ao, tmp)

        ! Symmetrize and pack transition density to triangular format
        tmp2(1:nbf2) => tmp
        call symmetrize_matrix(trden_ao, nbf)
        call pack_matrix(trden_ao, tmp2)

        ! Compute dipole moment:
        ! D_i = Tr(T * dipole_ints_i), i = x, y, z
        dip(1,s1) = -traceprod_sym_packed(tmp2, mints(:,1), nbf)/(2.0d0*sqrt(2.0d0))
        dip(2,s1) = -traceprod_sym_packed(tmp2, mints(:,2), nbf)/(2.0d0*sqrt(2.0d0))
        dip(3,s1) = -traceprod_sym_packed(tmp2, mints(:,3), nbf)/(2.0d0*sqrt(2.0d0))
    end do

  end subroutine

  subroutine td_print_results(energies, dip, nstates)

    use io_constants, only: iw
    use physical_constants, only: UNITS_EV
    implicit none

    real(kind=8) :: energies(:)
    real(kind=8) :: dip(:,:)
    real(kind=8) :: f
    integer :: nstates

    integer :: i
    write(iw,'(/x,78("^"))')
    write(iw,'(/8x,a/)') "Summary of the TD-DFT calculation"
    write(iw,'(A12,A16,3A12,a14)') 'Transition', 'dE(eV)', 'DX', 'DY', 'DZ', 'Osc.str.'

    do i = 1, nstates
      f = 2.0d0 / 3.0d0 * energies(i) * dot_product(dip(:,i), dip(:,i))
      write(iw,'(3X,"0 -> ",G0,t15,F15.6, 4F12.4)') i, energies(i)/UNITS_EV, dip(1:3,i), f
    end do
    write(iw,'(x,78("=")/)')

  end subroutine

end module tdhf_energy_mod
