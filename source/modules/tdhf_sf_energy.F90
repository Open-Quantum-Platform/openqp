module tdhf_sf_energy_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_sf_energy_mod"

contains

  subroutine tdhf_sf_energy_C(c_handle) bind(C, name="tdhf_sf_energy")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_sf_energy(inf)
  end subroutine tdhf_sf_energy_C

  subroutine tdhf_sf_energy(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver

    use types, only: information
    use strings, only: Cstring, fstring
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use util, only: measure_time

    use precision, only: dp
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_td_data_t
    use tdhf_lib, only: &
      inivec, iatogen, mntoia, rparedms, rpaeig, rpavnorm, &
      rpaechk, rpaprint, rpanewb
    use tdhf_sf_lib, only: &
      sfroesum, sfresvec, sfqvec, sfesum, sfdmat, trfrmb, &
      print_results, get_spin_square, get_transition_density, &
      get_transitions, get_transition_dipole
    use mathlib, only: orthogonal_transform_sym
    use mathlib, only: unpack_matrix
    use oqp_linalg
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_sf_energy"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos

    integer :: s_size, ok

    real(kind=dp), allocatable :: scr2(:)
    real(kind=dp), allocatable :: ab2_mo(:,:), scr3(:,:)
    real(kind=dp), allocatable :: eex(:), spin_square(:)
    real(kind=dp), allocatable :: amb(:,:), &
                                  apb(:,:)
    real(kind=dp), allocatable, target :: vl(:), vr(:)
    real(kind=dp), pointer :: vl_p(:,:), vr_p(:,:)
    real(kind=dp), allocatable :: xm(:)
    real(kind=dp), allocatable :: bvec_mo(:,:), for_trnsf_b_vec(:,:)
    real(kind=dp), allocatable, target :: bvec(:,:,:)
    real(kind=dp), allocatable, target :: scr1(:,:)
    real(kind=dp), pointer :: ab2(:,:,:)
    real(kind=dp), allocatable, dimension(:,:) :: fa, fb
    real(kind=dp), allocatable, dimension(:) :: rnorm
    real(kind=dp), allocatable, dimension(:,:,:,:) :: trden
    integer, allocatable, dimension(:,:) :: trans
    real(kind=dp), pointer :: scr1t(:)
    real(kind=dp), allocatable :: dip(:,:,:), abxc(:,:)

    integer :: nocca, nvira, noccb, nvirb
    integer :: nbf, nbf2, xvec_dim
    integer :: nstates, mxvec, nmax, ist, iend, nvec, novec
    integer :: iter, istart, nv, iv, ivec
    integer :: mxiter
    integer :: imax
    logical :: converged
    integer :: ierr
    real(kind=dp) :: mxerr, cnvtol, scale_exch
    integer :: maxvec, target_state
    logical :: roref = .false.

    type(int2_compute_t) :: int2_driver
    type(int2_td_data_t), target :: int2_data

    logical :: dft
    integer :: scf_type, mol_mult

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      fock_a(:), dmat_a(:), mo_a(:,:), mo_energy_a(:), &
      fock_b(:), dmat_b(:), mo_b(:,:), mo_energy_b(:), &
      smat(:), ta(:), tb(:), bvec_mo_out(:,:), td_t(:,:), &
      sf_energies(:)
    character(len=*), parameter :: tags_alloc(3) = (/ character(len=80) :: &
      OQP_td_bvec_mo, OQP_td_t, OQP_td_energies /)
    character(len=*), parameter :: tags_required(9) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_DM_A, OQP_E_MO_A, OQP_VEC_MO_A, OQP_FOCK_B, OQP_DM_B, OQP_E_MO_B, OQP_VEC_MO_B, OQP_SM /)

    mol_mult = infos%mol_prop%mult
    if (mol_mult/=3) call show_message('SF-TDDFT are available for ROHF/UHF ref. with ONLY triplet multiplicity(mult=3)',with_abort)

    scf_type = infos%control%scftype
    if (scf_type==3) roref = .true.

    dft = infos%control%hamilton == 20

  ! Files open
  ! 3. LOG: Write: Main output file
    open (unit=IW, file=infos%log_filename, position="append")
  !
    call print_module_info('SF_TDHF_Energy','Computing Energy of SF-TDDFT')
  ! Readings

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

  ! Allocate H, S ,T and D matrices
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    s_size = (basis%nshell**2+basis%nshell)/2

  ! Allocate temporary matrices for diagonalization
    allocate (FA(nbf, nbf), &
              FB(nbf, nbf), &
              stat=ok)
    if( ok/=0 ) call show_message('Cannot allocate memory',with_abort)

    nstates = infos%tddft%nstate
    target_state = infos%tddft%target_state
    maxvec = infos%tddft%maxvec
    cnvtol = infos%tddft%cnvtol

    nocca = infos%mol_prop%nelec_A
    nvira = nbf-noccA
    noccb = infos%mol_prop%nelec_B
    nvirb = nbf-noccB
    xvec_dim = nocca*nvirb

    mxvec = min(maxvec*nstates, xvec_dim)
    nstates = min(nstates, mxvec)
    nvec = nstates
    nvec = min(max(2*nstates, 5), mxvec)
    nmax = nvec

    call infos%dat%remove_records(tags_alloc)

    call infos%dat%reserve_data(OQP_td_bvec_mo, TA_TYPE_REAL64, &
        xvec_dim*nstates, (/xvec_dim, nstates/), comment=OQP_td_bvec_mo_comment)
    call infos%dat%reserve_data(OQP_td_t, TA_TYPE_REAL64, nbf2*2, (/ nbf2, 2 /), comment=OQP_td_t_comment)
    call infos%dat%reserve_data(OQP_td_energies, TA_TYPE_REAL64, nstates, comment=OQP_td_energies_comment)

    call data_has_tags(infos%dat, tags_alloc, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo_out)
    call tagarray_get_data(infos%dat, OQP_td_t, td_t)
    call tagarray_get_data(infos%dat, OQP_td_energies, sf_energies)

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_SM, smat)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)

    allocate(xm(xvec_dim), &
             trden(nbf,nbf,nstates,nstates), &
             bvec_mo(xvec_dim,mxvec), &
             abxc(nbf,nbf), &
             ab2_mo(xvec_dim,mxvec), &
             bvec(nbf,nbf,nmax), &
             eex(mxvec), &
             spin_square(nstates), &
             apb(mxvec,mxvec), &
             amb(mxvec,mxvec), &
             for_trnsf_b_vec(mxvec,mxvec), & !
             dip(3,nstates,nstates), &
             vr(mxvec*mxvec), &
             vl(mxvec*mxvec), &
             scr1(nbf,nbf), &
             scr2(mxvec*mxvec), &
             scr3(xvec_dim,nstates), &
             rnorm(nstates), &
             source=0.0_dp,stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    allocate(trans(xvec_dim,2), &
             source=0,stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    scale_exch = 1.0_dp
    if (infos%tddft%HFscale == -1.0_dp) &
          infos%tddft%HFscale = infos%dft%HFscale

    if (infos%dft%cam_flag) then
      if (infos%tddft%cam_alpha == -1.0_dp) &
            infos%tddft%cam_alpha = infos%dft%cam_alpha
      infos%tddft%HFscale = infos%tddft%cam_alpha
      if (infos%tddft%cam_beta == -1.0_dp) &
            infos%tddft%cam_beta = infos%dft%cam_beta
      if (infos%tddft%cam_mu == -1.0_dp) &
            infos%tddft%cam_mu = infos%dft%cam_mu
    end if
    if (dft) scale_exch = infos%tddft%HFscale

    if(.true.)then
      write(*,'(/,5x,"Input parameters:")')
      write(*,'(5x,"Number of states:                 ",1x,I0)') nstates
      write(*,'(5x,"Number of single excitations:     ",1x,I0)') xvec_dim
      write(*,'(5x,"Number of atomic orbitals:        ",1x,I0)') nbf
      write(*,'(5x,"Number of electrons:              ",1x,I0)') nocca+noccb
      write(*,'(5x,"Number of occupied alpha orbitals:",1x,I0)') nocca
      write(*,'(5x,"Number of occupied beta orbitals: ",1x,I0)') noccb
      write(*,'(5x,"Number of virtual alpha orbitals: ",1x,I0)') nvira
      write(*,'(5x,"Number of virtual beta orbitals:  ",1x,I0)') nvirb
      write(*,'(5x,"Maximum vectors:                  ",1x,I0)') mxvec
      write(*,'(5x,"Initial vectors:                  ",1x,I0)') nvec
      write(*, '(/7x,"Fitting parameters for SF-TDDFT")')
      if (.not.infos%dft%cam_flag) then
        write(*, '(10x,"Exact HF exchange:")')
        write(*, '(5x,"Reference: |", t20, f6.3, t29, "|")') infos%dft%HFscale
        write(*, '(5x,"Response:  |", t20, f6.3, t29, "|")') infos%tddft%HFscale
      else
        write(*, '(10x,"CAM parametres:")')
        write(*, '(16x,"|   alpha   |    beta   |     mu    |")')
        write(*, '(5x,"Reference: |", t20, f6.3, t29, "|", t32, f6.3, t41, "|", t44, f6.3, t53, "|")') &
           infos%dft%cam_alpha, infos%dft%cam_beta, infos%dft%cam_mu
        write(*, '(5x,"Response:  |", t20, f6.3, t29, "|", t32, f6.3, t41, "|", t44, f6.3, t53, "|")') &
           infos%tddft%cam_alpha, infos%tddft%cam_beta, infos%tddft%cam_mu
      end if
    end if

    write(*,'(/,5x,46("="))')
    write(*,'(5X,"Davidson algorithm for Spin-flip TDDFT")')
    write(*,'(5x,46("="))')

    ta          => td_t(:,1)
    tb          => td_t(:,2)

  ! Initialize ERI calculations
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

  ! Prepare for ROHF
    if( roref )then
      scr1t(1:nbf*nbf) => scr1(:,:)
  !   Alpha
      call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, scr1)
      call unpack_matrix(scr1t,fa)

  !   Beta
      call orthogonal_transform_sym(nbf, nbf, fock_b, mo_b, nbf, scr1)
      call unpack_matrix(scr1t,fb)
    end if

  ! Construct TD trial vector
    call inivec(mo_energy_a,mo_energy_b,bvec_mo,xm, &
                nocca,noccb,nvec)

    ist = 1
    istart = 1
    iend = nvec
    iter = 0
    mxiter = infos%control%maxit_dav
    ierr = 0

    do iter = 1, mxiter
      nv = iend-ist+1

      do ivec = ist, iend
        iv = ivec-ist+1
        call iatogen(bvec_mo(:,ivec),abxc,nocca,noccb)
        call dgemm('n','n',nbf,nbf,nbf, &
                   1.0_dp,mo_a,nbf,abxc,nbf, &
                   0.0_dp,scr1,nbf)
        call dgemm('n','t',nbf,nbf,nbf, &
                   1.0_dp,scr1,nbf,mo_b,nbf, &
                   0.0_dp,bvec(1,1,iv),nbf)
      end do

      int2_data = int2_td_data_t(d2=bvec(:,:,:nv), &
              int_apb=.false., int_amb=.false., tamm_dancoff=.true., &
              scale_exchange=scale_exch)

      call int2_driver%run(int2_data, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%tddft%cam_alpha, &
              beta=infos%tddft%cam_beta,&
              mu=infos%tddft%cam_mu)
      ab2 => int2_data%amb(:,:,:,1)

      do ivec = ist, iend
        iv = ivec-ist+1
  !     Product (A-B)*X
        call mntoia(ab2(:,:,iv),ab2_mo(:,ivec),mo_a,mo_b,nocca,noccb)
        if (roref) then
          call iatogen(bvec_mo(:,ivec),abxc,nocca,noccb)

          ! faz
          call dgemm('n','n',nocca,nbf,nocca, &
                     1.0_dp,fa,nbf,abxc,nbf, &
                     0.0_dp,scr1,nbf)
          ! zfb - faz
          call dgemm('n','n',nocca,nbf,nbf, &
                     1.0_dp,abxc,nbf,fb,nbf, &
                    -1.0_dp,scr1,nbf)

          call sfroesum(scr1,ab2_mo,nocca,noccb,ivec)
        else
          call sfesum(mo_energy_a,mo_energy_b,ab2_mo,bvec_mo,nocca,noccb,ivec)
        endif
      end do

      vl_p(1:nvec, 1:nvec) => vl(1:nvec*nvec)
      vr_p(1:nvec, 1:nvec) => vr(1:nvec*nvec)
      call rparedms(bvec_mo,ab2_mo,ab2_mo,apb,amb,nvec,tamm_dancoff=.true.)
      call rpaeig(eex,vl_p,vr_p,apb,amb,scr2,tamm_dancoff=.true.)
      call rpavnorm(vr_p,vl_p,tamm_dancoff=.true.)

      call rpaechk(eex,nvec,nstates,imax,tamm_dancoff=.true.)

      for_trnsf_b_vec = vr_p
      call sfresvec(scr3,bvec_mo,ab2_mo,vr_p,eex,nvec,rnorm,nstates)
      call sfqvec(scr3,xm,eex,nstates)

      call rpaprint(eex, rnorm, cnvtol, iter, imax, nstates, do_neg=.true.)

      mxerr = maxval(rnorm)

!     Check convergence
      converged = mxerr<=cnvtol
      if (converged) exit

!     No space left for new vectors, exit
      if (nvec==mxvec) ierr = 1
      if (ierr/=0) exit

      call rpanewb(nstates,bvec_mo,scr3,novec,nvec,ierr,tamm_dancoff=.true.)

  !   ierr=1 nvec over mxvec: not converged case
      if (ierr/=0) exit

      ist = novec+1
      iend = nvec
    end do

    if (iter >= mxiter .and. .not. converged) ierr = -1

    select case (ierr)
    case (-1)
      write(*,'(/,2X,"SF-TD-DFT energies NOT CONVERGED after ",I4," iterations"/)') mxiter
!      call show_message("Aborting. Try to increase maxit or check your system.", WITH_ABORT)
      infos%mol_energy%Davidson_converged=.false.
    case (0)
      write(*,'(/,2X,"SF-TD-DFT energies converged in ",I4," iterations"/)') iter
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
    call flush(iw)

    call trfrmb(bvec_mo,for_trnsf_b_vec,nvec,nstates)

    call get_transition_density(trden, bvec_mo, nbf, noccb, nocca, nstates)

    call get_transition_dipole(basis, dip, mo_a, trden, nstates)

    do ist = 1, nstates
      call sfdmat(bvec_mo(:,ist),abxc,mo_a,ta,tb,nocca,noccb)
      spin_square(ist) = get_spin_square(dmat_a,dmat_b,ta,tb,abxc,Smat,noccb)
    end do

    call get_transitions(trans, nocca, noccb, nbf)

    write(*,'(2x,35("="),/,2x,"Alpha -> Beta spin-flip excitations",/,2x,35("="))')

    sf_energies = eex(:nstates)
    bvec_mo_out = bvec_mo(:,:nstates)
    infos%mol_energy%excited_energy = sf_energies(infos%tddft%target_state)
    call print_results(infos, bvec_mo, eex, trans, dip, spin_square, nstates)
    call flush(iw)

    call int2_driver%clean()

    call measure_time(print_total=1, log_unit=iw)
    close(iw)

  end subroutine tdhf_sf_energy

end module tdhf_sf_energy_mod
