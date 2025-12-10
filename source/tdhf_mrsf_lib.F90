module tdhf_mrsf_lib

    use precision, only : dp
    use int2_compute, only: int2_fock_data_t, int2_storage_t
    use basis_tools, only: basis_set
    use oqp_linalg

    type, extends(int2_fock_data_t) :: int2_mrsf_data_t

      real(kind=dp), allocatable :: f3(:,:,:,:,:)
      real(kind=dp), pointer :: d3(:,:,:,:) => null()
      logical :: tamm_dancoff = .true. !< Tamm-Dancoff approximation

    contains

        procedure :: parallel_start => int2_mrsf_data_t_parallel_start
        procedure :: parallel_stop => int2_mrsf_data_t_parallel_stop
        procedure :: init_screen => int2_mrsf_data_t_init_screen
        procedure :: update => int2_mrsf_data_t_update
        procedure :: clean => int2_mrsf_data_t_clean

    end type

contains

!###############################################################################

  subroutine int2_mrsf_data_t_parallel_start(this, basis, nthreads)

    implicit none

    class(int2_mrsf_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nthreads
    integer :: nbf, nsh, nmatrix

    nbf = basis%nbf
    this%fockdim = nbf*(nbf+1) / 2
    this%nfocks = ubound(this%d3,1)
    this%nthreads = nthreads
    nsh = basis%nshell
    nmatrix = 7 ! spin pair copuling A': bo2v, bo1v, bco1, bco2, o21v, co12; A: ball

    if (this%cur_pass == 1) then
      if (allocated(this%f3)) deallocate(this%f3)
      if (allocated(this%dsh)) deallocate(this%dsh)

      allocate(this%f3(this%nfocks, nmatrix, nbf, nbf, nthreads), &
               this%dsh(nsh,nsh), &
               source=0.0d0)
    end if

    call this%init_screen(basis)

  end subroutine

!###############################################################################

  subroutine int2_mrsf_data_t_parallel_stop(this)

    implicit none

    integer :: f3last
    class(int2_mrsf_data_t), intent(inout) :: this

    if (this%cur_pass /= this%num_passes) return

    f3last = size(shape(this%f3))

    if (this%nthreads /= 1) then
      this%f3(:,:,:,:,lbound(this%f3, f3last)) = sum(this%f3, dim=f3last)
    end if

    call this%pe%allreduce(this%f3(:,:,:,:,1), &
              size(this%f3(:,:,:,:,1)))
    this%nthreads = 1

  end subroutine

!###############################################################################

  subroutine int2_mrsf_data_t_clean(this)

    implicit none

    class(int2_mrsf_data_t), intent(inout) :: this

    deallocate(this%f3)
    deallocate(this%dsh)
    nullify(this%d3)

  end subroutine

!###############################################################################

  subroutine int2_mrsf_data_t_init_screen(this, basis)

    implicit none

    class(int2_mrsf_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis

!   Form shell density
    call shell_den_screen_mrsf(this%dsh, this%d3(:,7,:,:), basis)
    this%max_den = maxval(abs(this%dsh))

  end subroutine

!###############################################################################

  subroutine shell_den_screen_mrsf(dsh,da,basis)

    use types, only: information
    use basis_tools, only: basis_set

    implicit none

    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(out) :: dsh(:,:)
    real(kind=dp), intent(in), dimension(:,:,:) :: da

    integer :: ish, jsh, maxi, maxj, mini, minj

    ! RHF
    do ish = 1, basis%nshell
      mini = basis%ao_offset(ish)
      maxi = mini + basis%naos(ish)-1
      do jsh = 1, ish
        minj = basis%ao_offset(jsh)
        maxj = minj+basis%naos(jsh)-1
        dsh(ish,jsh) = maxval(abs(da(:,minj:maxj,mini:maxi)))
        dsh(jsh,ish) = dsh(ish,jsh)
      end do
    end do

  end subroutine shell_den_screen_mrsf

!###############################################################################

  subroutine int2_mrsf_data_t_update(this, buf)

    implicit none

    class(int2_mrsf_data_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    integer :: i, j, k, l, n
    real(kind=dp) :: val, xval, cval
    integer :: mythread

    mythread = buf%thread_id

    if (.not.this%tamm_dancoff) return

    associate ( f3 => this%f3(:,:,:,:,mythread), &
                d3 => this%d3, &
                nf => this%nfocks &
      )

      do n = 1, buf%ncur
        i = buf%ids(1,n)
        j = buf%ids(2,n)
        k = buf%ids(3,n)
        l = buf%ids(4,n)
        val = buf%ints(n)

        xval = val * this%scale_exchange
        cval = val * this%scale_coulomb

        ! f3(nF,1:7,:,:) !> 1=ado2v, 2=ado1v, 3=adco1, 4=adco2, 5=ao21v, 6=aco12, 7=agdlr
        ! d3(nF,1:7,:,:) !> 1= bo2v, 2= bo1v, 3= bco1, 4= bco2, 5= o21v, 6= co12, 7= ball
        if (this%cur_pass==1) then
          f3(:nf,:4,i,j) = f3(:nf,:4,i,j) + cval*d3(:nf,:4,k,l)! (ij|lk)
          f3(:nf,:4,k,l) = f3(:nf,:4,k,l) + cval*d3(:nf,:4,i,j)! (kl|ji)
          f3(:nf,:4,i,j) = f3(:nf,:4,i,j) + cval*d3(:nf,:4,l,k)! (ij|kl)
          f3(:nf,:4,l,k) = f3(:nf,:4,l,k) + cval*d3(:nf,:4,i,j)! (lk|ji)
          f3(:nf,:4,j,i) = f3(:nf,:4,j,i) + cval*d3(:nf,:4,k,l)! (ji|lk)
          f3(:nf,:4,k,l) = f3(:nf,:4,k,l) + cval*d3(:nf,:4,j,i)! (kl|ij)
          f3(:nf,:4,j,i) = f3(:nf,:4,j,i) + cval*d3(:nf,:4,l,k)! (ji|kl)
          f3(:nf,:4,l,k) = f3(:nf,:4,l,k) + cval*d3(:nf,:4,j,i)! (lk|ij)

          f3(:nf,:7,i,k) = f3(:nf,:7,i,k) - xval*d3(:nf,:7,j,l)
          f3(:nf,:7,k,i) = f3(:nf,:7,k,i) - xval*d3(:nf,:7,l,j)
          f3(:nf,:7,i,l) = f3(:nf,:7,i,l) - xval*d3(:nf,:7,j,k)
          f3(:nf,:7,l,i) = f3(:nf,:7,l,i) - xval*d3(:nf,:7,k,j)
          f3(:nf,:7,j,k) = f3(:nf,:7,j,k) - xval*d3(:nf,:7,i,l)
          f3(:nf,:7,k,j) = f3(:nf,:7,k,j) - xval*d3(:nf,:7,l,i)
          f3(:nf,:7,j,l) = f3(:nf,:7,j,l) - xval*d3(:nf,:7,i,k)
          f3(:nf,:7,l,j) = f3(:nf,:7,l,j) - xval*d3(:nf,:7,k,i)
        else if (this%cur_pass==2) then
          f3(1:nf,7,i,k) = f3(1:nf,7,i,k) - xval*d3(1:nf,7,j,l)
          f3(1:nf,7,k,i) = f3(1:nf,7,k,i) - xval*d3(1:nf,7,l,j)
          f3(1:nf,7,i,l) = f3(1:nf,7,i,l) - xval*d3(1:nf,7,j,k)
          f3(1:nf,7,l,i) = f3(1:nf,7,l,i) - xval*d3(1:nf,7,k,j)
          f3(1:nf,7,j,k) = f3(1:nf,7,j,k) - xval*d3(1:nf,7,i,l)
          f3(1:nf,7,k,j) = f3(1:nf,7,k,j) - xval*d3(1:nf,7,l,i)
          f3(1:nf,7,j,l) = f3(1:nf,7,j,l) - xval*d3(1:nf,7,i,k)
          f3(1:nf,7,l,j) = f3(1:nf,7,l,j) - xval*d3(1:nf,7,k,i)
        end if
      end do
    end associate

    buf%ncur = 0

  end subroutine

!###############################################################################
!###############################################################################

  subroutine mrinivec(infos,ea,eb,bvec_mo,xm,nvec)

    use precision, only: dp
    use io_constants, only: iw
    use types, only: information

    implicit none

    type(information), intent(in) :: infos
    real(kind=dp), intent(in), dimension(:) :: ea, eb
    real(kind=dp), intent(out), dimension(:,:) :: bvec_mo
    real(kind=dp), intent(out), dimension(:) :: xm
    integer, intent(in) :: nvec

    logical :: debug_mode
    real(kind=dp) :: xmj
    integer :: nocca, nbf, i, ij, j, k, xvec_dim, lr1, lr2, mrst
    integer :: itmp(nvec)
    real(kind=dp) :: xtmp(nvec)

    debug_mode = infos%tddft%debug_mode
    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_A
    xvec_dim = ubound(xm, 1)
    lr1 = nocca-1
    lr2 = nocca

    ! For Singlet or Triplet
    mrst = infos%tddft%mult

    ! Set xm(xvec_dim)
    ij = 0
    do j = lr1, nbf
      do i = 1, lr2
        if (i==lr1 .and. j==lr1) then
          ij = ij+1
          xm(ij) = (eb(lr1)-ea(lr1)+eb(lr2)-ea(lr2))*0.5_dp
          cycle
        end if

        ij = ij+1
        xm(ij) = eb(j)-ea(i)

        if (i==nocca .and. j==nocca) then
          xm(ij) = huge(1.0d0)
        else if (mrst==3 .and. i==nocca .and. j==nocca-1) then
          xm(ij) = huge(1.0d0)
        else if (mrst==3 .and. i==nocca-1 .and. j==nocca) then
          xm(ij) = huge(1.0d0)
        end if
      end do
    end do

    ! Find indices of the first `nvec` smallest
    ! Values in the `xm` array
    itmp = 0 ! indices
    xtmp = huge(1.0d0) ! values
    do i = 1, xvec_dim
      do j = 1, nvec
        if (xtmp(j) > xm(i)) exit
      end do
      if (j <= nvec) then
        ! new small value found,
        ! insert it into temporary arrays
        xtmp(j+1:nvec) = xtmp(j:nvec-1)
        itmp(j+1:nvec) = itmp(j:nvec-1)
        xtmp(j) = xm(i)
        itmp(j) = i
      end if
    end do

    ! Ordering xm(xvec_dim): xm(small) <= xm(large)
    ! Get smaller diagonal values
    do j = 1, xvec_dim-1
      do i = j+1, xvec_dim
        if (xm(j)<=xm(i)) cycle
        xmj = xm(j)
        xm(j) = xm(i)
        xm(i) = xmj
      end do
    end do

    if (debug_mode) then
      write(iw,'("print xm(xvec_dim) ordering")')
      do i = 1, xvec_dim
        write(iw,'(a,i5,f20.10,i5)') 'i,xm(ij)=', i, xm(i)
      end do
    end if

    ! Get initial vectors: bvec(xvec_dim, nvec)
    bvec_mo = 0.0_dp
    do k = 1, nvec
      bvec_mo(itmp(k),k) = 1.0_dp
    end do

    ! set xm(xvec_dim) again
    ij = 0
    do j = lr1, nbf
      do i = 1, lr2
        if (i==lr1 .and. j==lr1) then
          ij = ij+1
          xm(ij) = (eb(lr1)-ea(lr1)+eb(lr2)-ea(lr2))*0.5_dp
          cycle
        endif

        ij = ij+1
        xm(ij) = eb(j)-ea(i)

        if (i==nocca .and. j==nocca) then
          xm(ij) = 9d99
        else if (mrst==3 .and. i==nocca .and. j==nocca-1) then
          xm(ij) = 9d99
        else if (mrst==3 .and. i==nocca-1 .and. j==nocca) then
          xm(ij) = 9d99
        end if
      end do
    end do

    if (debug_mode) then
      write(iw,'("print xm(xvec_dim) ordering")')
      do i = 1, xvec_dim
        write(iw,'(a,i5,f20.10,i5)') 'i,xm(ij)=', i, xm(i)
      end do
    end if

    return

  end subroutine mrinivec

!> Transform MRSF response vectors from MO to AO basis
!>
!> This subroutine performs the transformation of MRSF-TDDFT response amplitudes
!> X^(k), where k is singlet or triplet, from molecular orbital (MO) representation
!> to atomic orbital (AO) basis. The AO-basis response matrices are needed
!> for contraction with two-electron integrals and Fock matrix contributions
!> in the response equations.
!>
!> Physical context:
!> In MRSF-TDDFT, the response space is constructed from MS=+/-1 triplet references
!> to eliminate spin contamination in target singlet and triplet excited states.
!>
!> Orbital spaces in MRSF-TDDFT:
!> - C (Closed): Doubly-occupied orbitals (indices i,j,k,l)
!> - O (Open): Singly-occupied orbitals O1=HOMO-1, O2=HOMO (indices u,v,w,z)
!> - V (Virtual): Unoccupied orbitals (indices a,b,c,d)
!>
!> For each reference state k, the response amplitudes X^(k)_pq represent orbital
!> excitations between different orbital spaces. The response configurations are:
!> - Type I (OO): Open-to-Open transitions (O1<->O2)
!> - Type II (CO): Closed-to-Open transitions (C->O1, C->O2)
!> - Type III (OV): Open-to-Virtual transitions (O1->V, O2->V)
!> - Type IV (CV): Closed-to-Virtual transitions (C->V)
!>
!> The six response components in this subroutine correspond to:
!> 1. bo2v: O2(HOMO, alpha) -> V(beta) - OV block
!> 2. bo1v: O1(HOMO-1, alpha) -> V(beta) - OV block
!> 3. bco1: C(alpha) -> O1(HOMO-1, beta) - CO block
!> 4. bco2: C(alpha) -> O2(HOMO, beta) - CO block
!> 5. o21v: Mixed OV component coupling O1 and O2 with V
!> 6. co12: Mixed CO component coupling C with O1 and O2
!>
!> Transformation scheme:
!> For each component, we transform X^(k)_pq (MO basis) to P^(k)_(mu,nu) (AO basis):
!>   P^(k)_(mu,nu) = sum_pq C_(mu,p) X^(k)_pq C_(nu,q)
!> where C are MO coefficient matrices (va for alpha-spin, vb for beta-spin).
!>
!> Spin-pairing coupling between MS=+1 and MS=-1 reference states is realized
!> through specific linear combinations of MO coefficients va and vb (with proper
!> signs and 1/sqrt(2) normalization), following Slater-Condon rules. This enables
!> proper description of singlet and triplet target states from the mixed-reference formalism.
!>
!> Reference: Lee et al., J. Chem. Phys. 150, 184111 (2019), Eq. 2.11-2.18
!>
!> \author  Konstantin Komarov (constlike@gmail.com)
!>
  subroutine mrsfcbc(infos,va,vb,bvec,fmrsf)

    use messages, only: show_message, with_abort
    use types, only: information
    use io_constants, only: iw
    use precision, only: dp

    implicit none

    type(information), intent(in) :: infos
    real(kind=dp), intent(in), dimension(:,:) :: &
      va, vb, bvec
    real(kind=dp), intent(inout), target, dimension(:,:,:) :: &
      fmrsf

    real(kind=dp), allocatable, dimension(:,:) :: &
      tmp
    real(kind=dp), pointer, dimension(:,:) :: &
      bo2v, bo1v, bco1, bco2, ball, co12, o21v
    integer :: nocca, noccb, mrst, i, j, m, nbf, lr1, lr2, ok
    logical :: debug_mode
    real(kind=dp), parameter :: isqrt2 = 1.0_dp/sqrt(2.0_dp)

    ball => fmrsf(7,:,:)
    bo2v => fmrsf(1,:,:)
    bo1v => fmrsf(2,:,:)
    bco1 => fmrsf(3,:,:)
    bco2 => fmrsf(4,:,:)
    o21v => fmrsf(5,:,:)
    co12 => fmrsf(6,:,:)

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    mrst = infos%tddft%mult
    debug_mode = infos%tddft%debug_mode

    lr1 = nocca-1
    lr2 = nocca

    allocate(tmp(nbf,max(1,noccb)), source=0.0_dp, stat=ok)
    if( ok/=0 ) call show_message('Cannot allocate memory',with_abort)

    !-----------------------------------------------------------------------
    ! Component 1: bo2v - O2(HOMO, alpha) -> V(beta) excitations (OV block)
    !-----------------------------------------------------------------------
    ! Physical meaning: This component represents Open-to-Virtual (OV) excitations
    ! from the HOMO orbital of alpha-spin (O2, lr2 = nocca) to virtual orbitals
    ! of beta-spin. In MRSF theory, this corresponds to Type III response
    ! configurations.
    !
    ! MO->AO transformation: P^bo2v_(mu,nu) = C^alpha_(mu,HOMO) * X_(HOMO,a) * C^beta_(nu,a)
    ! where a runs over virtual beta-orbitals (nocca+1:nbf)
    !
    ! Step 1: Intermediate vector tmp = sum_a C^beta_(mu,a) X_(HOMO,a)
    !   tmp_mu = sum_{a in virt_beta} C^beta_(mu,a) * X_(HOMO,a)
    call dgemm('n', 't', nbf, 1, nbf-nocca, &
               1.0_dp, vb(:,nocca+1), nbf, &
                       bvec(lr2:lr2,nocca+1), nbf, &
               0.0_dp, tmp(:,1), nbf)

    ! Step 2: Outer product to form AO-basis matrix
    !   P^bo2v_(mu,nu) += C^alpha_(mu,HOMO) * tmp_nu
    call dgemm('n', 't', nbf, nbf, 1, &
               1.0_dp, va(:,lr2:lr2), nbf, &
                       tmp(:,1:1), nbf, &
               1.0_dp, bo2v, nbf)

    !-----------------------------------------------------------------------
    ! Component 2: bo1v - O1(HOMO-1, alpha) -> V(beta) excitations (OV block)
    !-----------------------------------------------------------------------
    ! Physical meaning: This component represents Open-to-Virtual (OV) excitations
    ! from HOMO-1 orbital of alpha-spin (O1, lr1 = nocca-1) to virtual orbitals
    ! of beta-spin. Together with bo2v, this forms the complete set of Type III
    ! spin-flip excitations from the two singly-occupied MOs (O1 and O2) in the
    ! MS=+/-1 triplet reference.
    !
    ! MO->AO transformation: P^bo1v_(mu,nu) = C^alpha_(mu,HOMO-1) * X_(HOMO-1,a) * C^beta_(nu,a)
    !
    ! Step 1: Intermediate vector tmp = sum_a C^beta_(mu,a) X_(HOMO-1,a)
    !   tmp_mu = sum_{a in virt_beta} C^beta_(mu,a) * X_(HOMO-1,a)
    call dgemm('n', 't', nbf, 1, nbf-nocca, &
               1.0_dp, vb(:,nocca+1), nbf, &
                       bvec(lr1:lr1,nocca+1), nbf, &
               0.0_dp, tmp(:,1), nbf)

    ! Step 2: Outer product to form AO-basis matrix
    !   P^bo1v_(mu,nu) += C^alpha_(mu,HOMO-1) * tmp_nu
    call dgemm('n', 't', nbf, nbf, 1, &
               1.0_dp, va(:,lr1:lr1), nbf, &
                       tmp(:,1:1), nbf, &
               1.0_dp, bo1v, nbf)

    !-----------------------------------------------------------------------
    ! Component 3: bco1 - C(alpha) -> O1(HOMO-1, beta) excitations (CO block)
    !-----------------------------------------------------------------------
    ! Physical meaning: This component represents Closed-to-Open (CO) excitations
    ! from doubly-occupied alpha-orbitals (C, 1:noccb) to the HOMO-1 orbital of
    ! beta-spin (O1, lr1). In MRSF theory, this corresponds to Type II response
    ! configurations. These are spin-flip de-excitations that complement the
    ! bo1v/bo2v excitations, maintaining the symmetry of the response space.
    !
    ! MO->AO transformation: P^bco1_(mu,nu) = C^alpha_(mu,i) * X_(i,HOMO-1) * C^beta_(nu,HOMO-1)
    ! where i runs over doubly-occupied orbitals (1:noccb)
    !
    ! Step 1: Intermediate vector tmp = sum_i C^alpha_(mu,i) X_(i,HOMO-1)
    !   tmp_mu = sum_{i in occ_alpha} C^alpha_(mu,i) * X_(i,HOMO-1)
    call dgemm('n', 'n', nbf, 1, noccb, &
               1.0_dp, va, nbf, &
                       bvec(1:noccb,lr1:lr1), nbf, &
               0.0_dp, tmp(:,1), nbf)

    ! Step 2: Outer product to form AO-basis matrix
    !   P^bco1_(mu,nu) += tmp_mu * C^beta_(nu,HOMO-1)
    call dgemm('n', 't', nbf, nbf, 1, &
               1.0_dp, tmp(:,1:1), nbf, &
                       vb(:,lr1:lr1), nbf, &
               1.0_dp, bco1, nbf)

    !-----------------------------------------------------------------------
    ! Component 4: bco2 - C(alpha) -> O2(HOMO, beta) excitations (CO block)
    !-----------------------------------------------------------------------
    ! Physical meaning: This component represents Closed-to-Open (CO) excitations
    ! from doubly-occupied alpha-orbitals (C, 1:noccb) to the HOMO orbital of
    ! beta-spin (O2, lr2). In MRSF theory, this corresponds to Type II response
    ! configurations. Together with bco1, this completes the set of spin-flip
    ! de-excitations from doubly-occupied orbitals to the two singly-occupied MOs.
    !
    ! MO->AO transformation: P^bco2_(mu,nu) = C^alpha_(mu,i) * X_(i,HOMO) * C^beta_(nu,HOMO)
    ! where i runs over doubly-occupied orbitals (1:noccb)
    !
    ! Step 1: Intermediate vector tmp = sum_i C^alpha_(mu,i) X_(i,HOMO)
    !   tmp_mu = sum_{i in occ_alpha} C^alpha_(mu,i) * X_(i,HOMO)
    call dgemm('n', 'n', nbf, 1, noccb, &
               1.0_dp, va, nbf, &
                       bvec(1:noccb,lr2:lr2), nbf, &
               0.0_dp, tmp(:,1), nbf)

    ! Step 2: Outer product to form AO-basis matrix
    !   P^bco2_(mu,nu) += tmp_mu * C^beta_(nu,HOMO)
    call dgemm('n', 't', nbf, nbf, 1, &
               1.0_dp, tmp(:,1:1), nbf, &
                       vb(:,lr2:lr2), nbf, &
               1.0_dp, bco2, nbf)

    !-----------------------------------------------------------------------
    ! Component 5: o21v - Mixed (O1<->O2)(alpha) x V(beta) (OV block)
    !-----------------------------------------------------------------------
    ! Physical meaning: This is a mixed Open-to-Virtual (OV) component coupling
    ! both HOMO and HOMO-1 alpha-orbitals (O1 and O2) with virtual beta-orbitals.
    ! The subtraction ensures proper antisymmetry and represents coherent
    ! superpositions of spin-flip excitations. This component is essential for
    ! the correct description of spin-adapted states in MRSF theory, arising
    ! from the coupling between different OV response configurations.
    !
    ! MO->AO transformation:
    !   P^o21v_(mu,nu) = sum_a [
    !       C^alpha_(mu,HOMO-1) * X_(HOMO,a) - C^alpha_(mu,HOMO) * X_(HOMO-1,a)
    !                          ] * C^beta_(nu,a)
    !
    ! Step 1: Intermediate vector from HOMO -> virt_beta amplitudes
    !   tmp_mu = sum_{a in virt_beta} C^beta_(mu,a) * X_(HOMO,a)
    call dgemm('n', 't', nbf, 1, nbf-nocca, &
               1.0_dp, vb(:,nocca+1), nbf, &
                       bvec(lr2:lr2,nocca+1), nbf, &
               0.0_dp, tmp(:,1), nbf)

    ! Combined outer products with subtraction
    !   P^o21v_(mu,nu) += C^alpha_(mu,HOMO-1) * tmp_nu (positive contribution)
    call dgemm('n', 't', nbf, nbf, 1, &
               1.0_dp, tmp(:,1:1), nbf, &
                       va(:,lr1:lr1), nbf, &
               1.0_dp, o21v, nbf)

    ! Step 1: Intermediate vector from HOMO-1 -> virt_beta amplitudes
    !   tmp_mu = sum_{a in virt_beta} C^beta_(mu,a) * X_(HOMO-1,a)
    call dgemm('n', 't', nbf, 1, nbf-nocca, &
               1.0_dp, vb(:,nocca+1), nbf, &
                       bvec(lr1:lr1,nocca+1), nbf, &
               0.0_dp, tmp(:,1), nbf)

    !   P^o21v_(mu,nu) -= C^alpha_(mu,HOMO) * tmp_nu (negative contribution)
    call dgemm('n', 't', nbf, nbf, 1, &
              -1.0_dp, tmp(:,1:1), nbf, &
                       va(:,lr2:lr2), nbf, &
               1.0_dp, o21v, nbf)

    !-----------------------------------------------------------------------
    ! Component 6: co12 - C(alpha) x Mixed (O1<->O2)(beta) (CO block)
    !-----------------------------------------------------------------------
    ! Physical meaning: This is a mixed Closed-to-Open (CO) component coupling
    ! doubly-occupied alpha-orbitals with both HOMO and HOMO-1 beta-orbitals
    ! (O1 and O2). The subtraction ensures proper antisymmetry and represents
    ! coherent superpositions of spin-flip de-excitations. Together with o21v,
    ! this maintains the full symmetry of the MRSF response space under orbital
    ! permutations, arising from coupling between different CO configurations.
    !
    ! MO->AO transformation:
    !   P^co12_(mu,nu) = sum_i [
    !       C^beta_(mu,HOMO) * X_(i,HOMO-1) - C^beta_(mu,HOMO-1) * X_(i,HOMO)
    !                          ] * C^alpha_(nu,i)
    !
    ! Step 1: Intermediate vector from occ_alpha -> HOMO-1_beta amplitudes
    !   tmp_mu = sum_{i in occ_alpha} C^alpha_(mu,i) * X_(i,HOMO-1)
    call dgemm('n', 'n', nbf, 1, noccb, &
               1.0_dp, va, nbf, &
                       bvec(1:noccb,lr1:lr1), nbf, &
               0.0_dp, tmp(:,1), nbf)

    ! Combined outer products with subtraction
    !   P^co12_(mu,nu) += C^beta_(mu,HOMO) * tmp_nu (positive contribution)
    call dgemm('n', 't', nbf, nbf, 1, &
               1.0_dp, vb(:,lr2:lr2), nbf, &
                       tmp(:,1:1), nbf, &
               1.0_dp, co12, nbf)

    ! Step 1: Intermediate vector from occ_alpha -> HOMO_beta amplitudes
    !   tmp_mu = sum_{i in occ_alpha} C^alpha_(mu,i) * X_(i,HOMO)
    call dgemm('n', 'n', nbf, 1, noccb, &
               1.0_dp, va, nbf, &
                       bvec(1:noccb,lr2:lr2), nbf, &
               0.0_dp, tmp(:,1), nbf)

    !   P^co12_(mu,nu) -= C^beta_(mu,HOMO-1) * tmp_nu (negative contribution)
    call dgemm('n', 't', nbf, nbf, 1, &
              -1.0_dp, vb(:,lr1:lr1), nbf, &
                       tmp(:,1:1), nbf, &
               1.0_dp, co12, nbf)

    !-----------------------------------------------------------------------
    ! Sum the four primary components into the total response matrix
    !-----------------------------------------------------------------------
    ! Physical meaning: Combine bo2v, bo1v, bco1, and bco2 (the four components
    ! without mixed character) into the total AO-basis response matrix ball.
    ! These represent the OV and CO blocks (Type II and Type III configurations).
    ! The mixed components o21v and co12 are not included here as they are
    ! handled separately in the spin-dependent sections below.
    ball = ball + bo2v + bo1v + bco1 + bco2

    !-----------------------------------------------------------------------
    ! Additional general contribution: C(alpha) x V(beta) block (CV)
    !-----------------------------------------------------------------------
    ! Physical meaning: Transform the general Closed-to-Virtual (CV) block of
    ! response amplitudes (doubly-occupied_alpha -> virtual_beta) to AO basis.
    ! This represents Type IV spin-flip excitations from the doubly-occupied
    ! core orbitals. These excitations do not involve the singly-occupied
    ! frontier orbitals (O1, O2) and represent the CV response configurations.
    !
    ! Transformation: P^ball_(mu,nu) += sum_ia C^alpha_(mu,i) * X_(i,a) * C^beta_(nu,a)
    ! where i in doubly-occupied (1:noccb), a in virtual_beta (nocca+1:nbf)
    !
    ! Step 1: Intermediate tmp_(mu,i) = sum_a C^beta_(mu,a) * X_(i,a)
    call dgemm('n', 't', nbf, noccb, nbf-nocca, &
               1.0_dp, vb(:,nocca+1), nbf, &
                       bvec(:,nocca+1), nbf, &
               0.0_dp, tmp(:,1:noccb), nbf)

    ! Step 2: Outer product P^ball_(mu,nu) += sum_i C^alpha_(mu,i) * tmp_(nu,i)
    call dgemm('n', 't', nbf, nbf, noccb, &
               1.0_dp, va, nbf, &
                       tmp(:,1:noccb), nbf, &
               1.0_dp, ball, nbf)

    !-----------------------------------------------------------------------
    ! Spin-dependent corrections for (O1<->O2) x (O1<->O2) block (OO)
    !-----------------------------------------------------------------------
    ! Physical meaning: Add corrections for the Open-to-Open (OO) Type I response
    ! configurations that depend on the target spin state (singlet mrst=1 or
    ! triplet mrst=3). These involve the special elements X_(HOMO,HOMO-1),
    ! X_(HOMO-1,HOMO), and X_(HOMO-1,HOMO-1) that couple the two singly-occupied
    ! MOs (O1 and O2). The 1/sqrt(2) factor ensures proper normalization of
    ! spin-adapted states in the MRSF formalism.
    !
    ! Spin-pairing coupling between MS=+1 and MS=-1 triplet references (following
    ! Slater-Condon rules) is implemented via specific linear combinations of
    ! MO coefficients va and vb. The different sign patterns for singlet/triplet
    ! reflect the different spin symmetries:
    ! - Singlet: antisymmetric spatial wavefunction (subtraction)
    ! - Triplet: symmetric spatial wavefunction (addition)
    if (mrst==1) then
      ! Singlet state corrections (mrst=1):
      ! Three terms that couple HOMO and HOMO-1 orbitals:
      ! 1. X_(HOMO,HOMO-1) * C^alpha_HOMO (outer) C^beta_HOMO-1
      ! 2. X_(HOMO-1,HOMO) * C^alpha_HOMO-1 (outer) C^beta_HOMO
      ! 3. X_(HOMO-1,HOMO-1) * (C^alpha_HOMO-1 (outer) C^beta_HOMO-1 - C^alpha_HOMO (outer) C^beta_HOMO) / sqrt(2)
      ! The subtraction in term 3 ensures proper singlet spin coupling.
      do m = 1, nbf
        ball(:,m) = ball(:,m) &
          +va(:,lr2)*bvec(lr2,lr1)*vb(m,lr1) &
          +va(:,lr1)*bvec(lr1,lr2)*vb(m,lr2) &
          +(va(:,lr1)*vb(m,lr1)-va(:,lr2)*vb(m,lr2)) &
             *bvec(lr1,lr1)*isqrt2
      end do
    else if (mrst==3) then
      ! Triplet state corrections (mrst=3):
      ! Single term for diagonal HOMO-1,HOMO-1 element:
      ! X_(HOMO-1,HOMO-1) * (C^alpha_HOMO-1 (outer) C^beta_HOMO-1 + C^alpha_HOMO (outer) C^beta_HOMO) / sqrt(2)
      ! The addition ensures proper triplet spin coupling.
      ! Off-diagonal OO terms are zero for triplet states.
      do m = 1, nbf
        ball(:,m) = ball(:,m) &
          +(va(:,lr1)*vb(m,lr1)+va(:,lr2)*vb(m,lr2)) &
             *bvec(lr1,lr1)*isqrt2
      end do
    end if

    if (debug_mode) then
      write(iw,*) 'Check sum = va', sum(abs(va))
      write(iw,*) 'Check sum = vb', sum(abs(vb))
      write(iw,*) 'Check sum = bvec', sum(abs(bvec))
      write(iw,*) 'Check sum = ball', sum(abs(ball))
      write(iw,*) 'Check sum = o21v', sum(abs(o21v))
      write(iw,*) 'Check sum = co12', sum(abs(co12))
      write(iw,*) 'Check sum = bo2v', sum(abs(bo2v))
      write(iw,*) 'Check sum = bo1v', sum(abs(bo1v))
      write(iw,*) 'Check sum = bco1', sum(abs(bco1))
      write(iw,*) 'Check sum = bco2', sum(abs(bco2))
    end if

    return

  end subroutine mrsfcbc

!> Transform MRSF Fock-like matrices from AO to MO basis
!>
!> This subroutine performs the transformation of MRSF-TDDFT Fock-like matrices
!> (or generalized density contributions) from atomic orbital (AO) basis back to
!> molecular orbital (MO) representation.
!>
!> Physical context:
!> After contracting the response amplitudes (in AO basis) with two-electron
!> integrals and other operators, we obtain Fock-like matrices P^(k)_(mu,nu) in
!> AO basis. These must be transformed back to MO basis to extract elements
!> corresponding to specific orbital transitions in the MRSF response space.
!>
!> Orbital spaces (same as in mrsfcbc):
!> - C (Closed): Doubly-occupied orbitals
!> - O (Open): Singly-occupied O1=HOMO-1, O2=HOMO
!> - V (Virtual): Unoccupied orbitals
!>
!> Transformation scheme:
!> For the general contribution: F^MO_pq = sum_(mu,nu) C^alpha_(mu,p) P^AO_(mu,nu) C^beta_(nu,q)
!>
!> Then, specific corrections are added for each of the six response components
!> corresponding to different blocks of the MRSF response space:
!> 1. Section 3: Corrections from ado1v (OV) and aco12 (CO mixed) -> C x O2 block
!> 2. Section 4: Corrections from ado2v (OV) and aco12 (CO mixed) -> C x O1 block
!> 3. Section 5: Corrections from adco2 (CO) and ao21v (OV mixed) -> O1 x V block
!> 4. Section 6: Corrections from adco1 (CO) and ao21v (OV mixed) -> O2 x V block
!>
!> Reference: Lee et al., J. Chem. Phys. 150, 184111 (2019), Eq. 2.11-2.18
!>
!> \author  Konstantin Komarov (constlike@gmail.com)
!>
   subroutine mrsfmntoia(infos, fmrsf, pmo, va, vb, ivec)

    use precision, only: dp
    use types, only: information
    use messages, only: show_message, with_abort
    use io_constants, only: iw

    implicit none

    type(information), intent(in) :: infos
    real(kind=dp), intent(in), target, dimension(:,:,:) :: &
      fmrsf
    real(kind=dp), intent(out), dimension(:,:) :: pmo
    real(kind=dp), intent(in), dimension(:,:) :: va, vb
    integer, intent(in) :: ivec

    real(kind=dp), allocatable :: &
      scr(:,:), tmp(:), wrk(:,:)
    real(kind=dp), pointer, dimension(:,:) :: &
      adco1, adco2, ado1v, ado2v, agdlr, aco12, ao21v
    integer :: noca, nocb, mrst, i, ij, &
      j, lr1, lr2, nbf, ok
    real(kind=dp), parameter :: zero = 0.0_dp
    real(kind=dp), parameter :: one = 1.0_dp
    real(kind=dp), parameter :: sqrt2 = 1.0_dp/sqrt(2.0_dp)
    logical :: debug_mode

    nbf = infos%basis%nbf
    mrst = infos%tddft%mult
    noca = infos%mol_prop%nelec_a
    nocb = infos%mol_prop%nelec_b
    debug_mode = infos%tddft%debug_mode

    allocate(tmp(nbf), scr(nbf,nbf), wrk(nbf,nbf), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', with_abort)

    agdlr => fmrsf(7,:,:)
    ado2v => fmrsf(1,:,:)
    ado1v => fmrsf(2,:,:)
    adco1 => fmrsf(3,:,:)
    adco2 => fmrsf(4,:,:)
    ao21v => fmrsf(5,:,:)
    aco12 => fmrsf(6,:,:)

    lr1 = noca-1
    lr2 = noca

    !-----------------------------------------------------------------------
    ! Initial AO->MO transformation of general Fock-like contribution
    !-----------------------------------------------------------------------
    ! Physical meaning: Transform the general (summed) Fock-like matrix agdlr
    ! from AO basis to MO basis. This represents the baseline contribution
    ! before adding specific corrections for each response component.
    !
    ! Transformation: F^MO_pq = C^alpha^T * P^AO * C^beta
    !   Step 1: wrk = C^alpha^T * agdlr (transform first index)
    !   Step 2: scr = wrk * C^beta (transform second index)
    !
    ! Result: scr contains the general MO-basis matrix that will be corrected
    ! by the specific response component contributions in sections 3-6.
    call dgemm('t','n',nbf,nbf,nbf, &
               one, va, nbf, &
                    agdlr,nbf, &
               zero, wrk, nbf)
    call dgemm('n','n',nbf,nbf,nbf, &
               one, wrk, nbf, &
                    vb, nbf, &
               zero, scr, nbf)

! 1
!   ----- (m,n) to (i+,n) -----
    wrk = scr

    !-----------------------------------------------------------------------
    ! Section 3: Corrections for C(alpha) -> O2(HOMO, beta) response element
    !-----------------------------------------------------------------------
    ! Physical meaning: Add corrections to F^MO_(i,HOMO) from two sources:
    ! 1. ado1v: O1(HOMO-1, alpha) -> V(beta) excitation component (OV block)
    ! 2. aco12: C(alpha) x Mixed (O1<->O2)(beta) component (CO mixed block)
    !
    ! These represent the backflow of MO to AO, accounting for coupling
    ! between different excitation types (OV and CO blocks). 
    ! The result corrects wrk(i,lr2) for doubly-occupied
    ! orbitals i in the Closed space.
    !
    ! Combined transformation:
    !   tmp = ado1v * C^beta_HOMO + aco12 * C^beta_HOMO-1
    !   F^MO_(i,HOMO) += C^alpha_i^T * tmp  (for i=1:noca-2)
    !
    ! Step 1: Contract ado1v with HOMO_beta MO coefficient
    !   tmp_mu = sum_nu P^ado1v_(mu,nu) * C^beta_(nu,HOMO)
    call dgemm('n','n',nbf,1,nbf, &
               one, ado1v, nbf, &
                    vb(:,lr2:lr2), nbf, &
              zero, tmp, nbf)
    ! Step 2: Add contribution from aco12 with HOMO-1_beta MO coefficient
    !   tmp_mu += sum_nu P^aco12_(mu,nu) * C^beta_(nu,HOMO-1)
    call dgemm('n','n',nbf,1,nbf, &
               one, aco12, nbf, &
                    vb(:,lr1:lr1), nbf, &
               one, tmp, nbf)
    ! Step 3: Project onto doubly-occupied alpha-orbitals
    !   F^MO_(i,HOMO) += sum_mu C^alpha_(mu,i) * tmp_mu  (i=1:noca-2)
    call dgemm('t','n',noca-2,1,nbf, &
               one, va, nbf, &
                    tmp, nbf, &
               one, wrk(1:noca-2,lr2:lr2), noca-2)
    !-----------------------------------------------------------------------
    ! Section 4: Corrections for C(alpha) -> O1(HOMO-1, beta) response element
    !-----------------------------------------------------------------------
    ! Physical meaning: Add corrections to F^MO_(i,HOMO-1) from two sources:
    ! 1. ado2v: O2(HOMO, alpha) -> V(beta) excitation component (OV block, positive)
    ! 2. aco12: C(alpha) x Mixed (O1<->O2)(beta) component (CO mixed, negative)
    !
    ! The subtraction in this section ensures proper antisymmetry between
    ! HOMO and HOMO-1 contributions, maintaining consistency with the mixed
    ! character of the aco12 and ao21v components. This antisymmetry arises from
    ! the spin-pairing coupling realized via linear combinations of va and vb.
    !
    ! Combined transformation:
    !   tmp = ado2v * C^beta_HOMO-1 - aco12 * C^beta_HOMO
    !   F^MO_(i,HOMO-1) += C^alpha_i^T * tmp  (for i=1:noca-2)
    !
    ! Step 1: Contract aco12 with HOMO_beta MO coefficient
    !   tmp_mu = sum_nu P^aco12_(mu,nu) * C^beta_(nu,HOMO)
    call dgemm('n','n',nbf,1,nbf, &
               one, aco12, nbf, &
                    vb(:,lr2:lr2), nbf, &
              zero, tmp, nbf)
    ! Step 2: Add ado2v contribution and subtract aco12 contribution
    !   tmp_mu = sum_nu P^ado2v_(mu,nu) * C^beta_(nu,HOMO-1) - tmp_mu
    call dgemm('n','n',nbf,1,nbf, &
               one, ado2v, nbf, &
                    vb(:,lr1:lr1), nbf, &
              -one, tmp, nbf)
    ! Step 3: Project onto doubly-occupied alpha-orbitals
    !   F^MO_(i,HOMO-1) += sum_mu C^alpha_(mu,i) * tmp_mu  (i=1:noca-2)
    call dgemm('t','n',noca-2,1,nbf, &
               one, va, nbf, &
                    tmp, nbf, &
               one, wrk(1:noca-2,lr1:lr1), noca-2)
    !-----------------------------------------------------------------------
    ! Section 5: Corrections for O1(HOMO-1, alpha) -> V(beta) response element
    !-----------------------------------------------------------------------
    ! Physical meaning: Add corrections to F^MO_(HOMO-1,a) from two sources:
    ! 1. adco2: C(alpha) -> O2(HOMO, beta) excitation component (CO block)
    ! 2. ao21v: Mixed (O1<->O2)(alpha) x V(beta) component (OV mixed block)
    !
    ! This section computes contributions to the O1(HOMO-1, alpha) -> V(beta) block
    ! by contracting the transposed AO matrices with alpha-spin MO coefficients,
    ! then projecting onto virtual beta-orbitals.
    !
    ! Combined transformation:
    !   tmp = adco2^T * C^alpha_HOMO-1 + ao21v^T * C^alpha_HOMO
    !   F^MO_(HOMO-1,a) += C^beta_a^T * tmp  (for a in virt_beta)
    !
    ! Step 1: Contract adco2^T with HOMO-1_alpha MO coefficient
    !   tmp_mu = sum_nu P^adco2_(nu,mu) * C^alpha_(nu,HOMO-1)
    call dgemm('t','n',nbf,1,nbf, &
               one, adco2, nbf, &
                    va(:,lr1:lr1), nbf, &
              zero, tmp, nbf)
    ! Step 2: Add contribution from ao21v^T with HOMO_alpha MO coefficient
    !   tmp_mu += sum_nu P^ao21v_(nu,mu) * C^alpha_(nu,HOMO)
    call dgemm('t','n',nbf,1,nbf, &
               one, ao21v, nbf, &
                    va(:,lr2:lr2), nbf, &
               one, tmp, nbf)
    ! Step 3: Project onto virtual beta-orbitals
    !   F^MO_(HOMO-1,a) += sum_mu C^beta_(mu,a) * tmp_mu  (a=noca+1:nbf)
    call dgemm('t','n',nbf-noca,1,nbf, &
               one, vb(:,noca+1), nbf, &
                    tmp, nbf, &
               one, wrk(lr1:lr1,noca+1:nbf), nbf-noca)
    !-----------------------------------------------------------------------
    ! Section 6: Corrections for O2(HOMO, alpha) -> V(beta) response element
    !-----------------------------------------------------------------------
    ! Physical meaning: Add corrections to F^MO_(HOMO,a) from two sources:
    ! 1. adco1: C(alpha) -> O1(HOMO-1, beta) excitation component (CO block, positive)
    ! 2. ao21v: Mixed (O1<->O2)(alpha) x V(beta) component (OV mixed, negative)
    !
    ! The subtraction in this section ensures proper antisymmetry between
    ! HOMO and HOMO-1 contributions, complementing Section 5 and maintaining
    ! consistency with the mixed character of ao21v. This antisymmetry arises from
    ! the spin-pairing coupling realized via linear combinations of va and vb.
    !
    ! Combined transformation:
    !   tmp = adco1^T * C^alpha_HOMO - ao21v^T * C^alpha_HOMO-1
    !   F^MO_(HOMO,a) += C^beta_a^T * tmp  (for a in virt_beta)
    !
    ! Step 1: Contract ao21v^T with HOMO-1_alpha MO coefficient
    !   tmp_mu = sum_nu P^ao21v_(nu,mu) * C^alpha_(nu,HOMO-1)
    call dgemm('t','n',nbf,1,nbf, &
               one, ao21v, nbf, &
                    va(:,lr1:lr1), nbf, &
              zero, tmp, nbf)
    ! Step 2: Add adco1 contribution and subtract ao21v contribution
    !   tmp_mu = sum_nu P^adco1_(nu,mu) * C^alpha_(nu,HOMO) - tmp_mu
    call dgemm('t','n',nbf,1,nbf, &
               one, adco1, nbf, &
                    va(:,lr2:lr2), nbf, &
              -one, tmp, nbf)
    ! Step 3: Project onto virtual beta-orbitals
    !   F^MO_(HOMO,a) += sum_mu C^beta_(mu,a) * tmp_mu  (a=noca+1:nbf)
    call dgemm('t','n',nbf-noca,1,nbf, &
               one, vb(:,noca+1), nbf, &
                    tmp, nbf, &
               one, wrk(lr2:lr2,noca+1:nbf), nbf-noca)

    !-----------------------------------------------------------------------
    ! Spin-dependent corrections for (O1,O1) diagonal element (OO block)
    !-----------------------------------------------------------------------
    ! Physical meaning: Apply spin-state-dependent transformations to the
    ! F^MO_(HOMO-1,HOMO-1) diagonal element based on whether we are computing
    ! a singlet (mrst=1) or triplet (mrst=3) excited state.
    !
    ! The 1/sqrt(2) factor and the addition/subtraction of diagonal elements ensure
    ! proper normalization and spin coupling for the MRSF response equations.
    ! This reflects the spin-pairing coupling between MS=+1 and MS=-1 triplet
    ! reference states (via linear combinations of va and vb), following Slater-Condon rules.
    if (mrst==1) then
      ! Singlet state (mrst=1):
      ! F^MO_(HOMO-1,HOMO-1) = (scr_(HOMO-1,HOMO-1) - scr_(HOMO,HOMO)) / sqrt(2)
      ! The subtraction reflects the antisymmetric spin coupling in singlets.
      ! The (HOMO,HOMO) element is zeroed as it's not part of the singlet response space.
      wrk(lr1,lr1) = (scr(lr1,lr1)-scr(lr2,lr2))*sqrt2
      wrk(lr2,lr2) = 0.0_dp
    else if (mrst==3) then
      ! Triplet state (mrst=3):
      ! F^MO_(HOMO-1,HOMO-1) = (scr_(HOMO-1,HOMO-1) + scr_(HOMO,HOMO)) / sqrt(2)
      ! The addition reflects the symmetric spin coupling in triplets.
      ! All other OO coupling elements (HOMO,HOMO-1), (HOMO-1,HOMO), and
      ! (HOMO,HOMO) are zeroed as they're outside the triplet response space.
      wrk(lr1,lr1) = (scr(lr1,lr1)+scr(lr2,lr2))*sqrt2
      wrk(lr2,lr1) = 0.0_dp
      wrk(lr1,lr2) = 0.0_dp
      wrk(lr2,lr2) = 0.0_dp
    end if

    pmo(:,ivec) = 0.0_dp

    ij = 0
    do j = nocb+1, nbf
      do i = 1, noca
        ij = ij+1
        pmo(ij,ivec) = pmo(ij,ivec) + wrk(i,j)
      end do
    end do

    if (debug_mode) then
      write(iw,*) 'Check sum = ivec', ivec
      write(iw,*) 'Check sum = agdlr', sum(abs(agdlr))
      write(iw,*) 'Check sum = ao21v', sum(abs(ao21v))
      write(iw,*) 'Check sum = aco12', sum(abs(aco12))
      write(iw,*) 'Check sum = ado2v', sum(abs(ado2v))
      write(iw,*) 'Check sum = ado1v', sum(abs(ado1v))
      write(iw,*) 'Check sum = adco1', sum(abs(adco1))
      write(iw,*) 'Check sum = adco2', sum(abs(adco2))
      write(iw,*) 'Check sum = pmo', sum(abs(pmo(:,ivec)))
    end if

    return

  end subroutine mrsfmntoia

  subroutine mrsfesum(infos, wrk, fij, fab, pmo, iv)

    use precision, only: dp
    use types, only: information
    use messages, only: show_message, with_abort
    use io_constants, only: iw

    implicit none

    type(information), intent(in) :: infos
    real(kind=dp), intent(in), dimension(:,:) :: &
      wrk, fij, fab
    real(kind=dp), intent(inout), dimension(:,:) :: &
      pmo
    integer, intent(in) :: iv

    real(kind=dp), allocatable, dimension(:,:) :: scr, tmp1, wrk1
    real(kind=dp) :: dumn, xlr
    integer :: nbf, nocca, noccb, mrst, i, ij, j, lr1, lr2, ok
    real(kind=dp), parameter :: sqrt2 = 1.0_dp/sqrt(2.0_dp)
    logical :: debug_mode

    nbf = infos%basis%nbf
    nocca = infos%mol_prop%nelec_a
    noccb = infos%mol_prop%nelec_b
    mrst = infos%tddft%mult
    debug_mode = infos%tddft%debug_mode

    allocate(scr(nbf,nbf), tmp1(nbf,nbf), wrk1(nbf,nbf), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', with_abort)

    lr1 = nocca-1
    lr2 = nocca
    scr = wrk
    scr(lr1,lr1) = 0.0_dp
    scr(lr2,lr2) = 0.0_dp

    ! Contraction 1
    call dgemm('n', 't', nocca, nbf-noccb, nbf-noccb, &
               1.0_dp, scr(1,noccb+1), nbf, &
                       fab(noccb+1:,noccb+1), nbf, &
               0.0_dp, tmp1(1,noccb+1), nbf)

    ! Contraction 2
    call dgemm('n', 'n', nocca, nbf-noccb, nocca, &
              -1.0_dp, fij(:,1), nbf, &
                       scr(:,noccb+1), nbf, &
               1.0_dp, tmp1(1,noccb+1), nbf)

    xlr = wrk(lr1,lr1)

    if (mrst==1) then

      do j = noccb+1, nbf
        do i = 1, nocca
          wrk1(i,j) = wrk1(i,j)+tmp1(i,j)
          if(i==lr1) wrk1(i,j) = wrk1(i,j)+fab(j,lr1)*xlr*sqrt2
          if(i==lr2) wrk1(i,j) = wrk1(i,j)-fab(j,lr2)*xlr*sqrt2
          if(j==lr1) wrk1(i,j) = wrk1(i,j)-fij(i,lr1)*xlr*sqrt2
          if(j==lr2) wrk1(i,j) = wrk1(i,j)+fij(i,lr2)*xlr*sqrt2
        end do
      end do

      dumn = - dot_product(fij(lr1,1:nocca),scr(1:nocca,lr1)) &
             + dot_product(fij(lr2,1:nocca),scr(1:nocca,lr2)) &
             + dot_product(fab(lr1,noccb+1:nbf),scr(lr1,noccb+1:nbf)) &
             - dot_product(fab(lr2,noccb+1:nbf),scr(lr2,noccb+1:nbf))

      wrk1(lr1,lr1) = dumn*sqrt2 &
                    + xlr*(fab(lr1,lr1)+fab(lr2,lr2) &
                          -fij(lr1,lr1)-fij(lr2,lr2))*0.5_dp

    elseif (mrst==3) then

      do j = noccb+1, nbf
        do i = 1, nocca
          wrk1(i,j) = wrk1(i,j)+tmp1(i,j)
          if(i==lr1) wrk1(i,j) = wrk1(i,j)+fab(j,lr1)*xlr*sqrt2
          if(i==lr2) wrk1(i,j) = wrk1(i,j)+fab(j,lr2)*xlr*sqrt2
          if(j==lr1) wrk1(i,j) = wrk1(i,j)-fij(i,lr1)*xlr*sqrt2
          if(j==lr2) wrk1(i,j) = wrk1(i,j)-fij(i,lr2)*xlr*sqrt2
        end do
      end do

      dumn = dot_product(-fij(lr1,1:nocca),scr(1:nocca,lr1)) &
           + dot_product(-fij(lr2,1:nocca),scr(1:nocca,lr2)) &
           + dot_product( fab(lr1,noccb+1:nbf),scr(lr1,noccb+1:nbf)) &
           + dot_product( fab(lr2,noccb+1:nbf),scr(lr2,noccb+1:nbf))

      wrk1(lr1,lr1) = dumn*sqrt2 &
                    + xlr*(fab(lr1,lr1)+fab(lr2,lr2) &
                          -fij(lr1,lr1)-fij(lr2,lr2))*0.5_dp

    end if

    if (mrst==1) then
      wrk1(lr2,lr2) = 0.0_dp
    else if (mrst==3) then
      wrk1(lr2,lr1) = 0.0_dp
      wrk1(lr1,lr2) = 0.0_dp
      wrk1(lr2,lr2) = 0.0_dp
    end if

    ij = 0
    do j = noccb+1, nbf
      do i = 1, nocca
        ij = ij+1
        pmo(ij,iv) = pmo(ij,iv)+wrk1(i,j)
      end do
    end do

    if (debug_mode) then
      write(iw,*) 'Check sum = ivec', iv
      write(iw,*) 'Check sum = pmo', sum(abs(pmo(:,iv)))
    end if

    return

  end subroutine mrsfesum

  subroutine mrsfqroesum(fbzzfa,pmo,noca,nocb,nbf,ivec)

    use precision, only: dp

    implicit none

    real(kind=dp), target, intent(in) :: fbzzfa(*)
    real(kind=dp), intent(inout) :: pmo(:,:)
    integer, intent(in) :: noca, nocb, nbf, ivec

    real(kind=dp), pointer :: wrk(:,:)
    integer :: i, ij, j

    wrk(1:nocb,1:nbf) => fbzzfa(1:nocb*nbf)

    ij = 0
    do j = noca+1, nbf
      do i = 1, nocb
        ij = ij + 1
        pmo(ij,ivec) = pmo(ij,ivec)+wrk(i,j)
      end do
    end do

  end subroutine mrsfqroesum

  subroutine get_mrsf_transitions(trans, noca, nocb, nbf)

    implicit none

    integer, intent(out), dimension(:,:) :: trans
    integer, intent(in) :: noca, nocb, nbf

    integer :: ij, i, j, lr1, lr2

    lr1 = nocb+1
    lr2 = noca
    ij = 0
    do j = lr1, nbf
       do i = 1, lr2
          ij = ij+1
          trans(ij,1) = i
          trans(ij,2) = j
       end do
    end do

  end subroutine get_mrsf_transitions

!> @details This subroutine transforms Multi-Reference Spin-Flip (MRSF) response vectors
!>          from a compressed representation to an expanded form. It handles both
!>          singlet (mrst=1) and triplet (mrst=3) cases.
!>
!> @param[in]     infos  Information structure containing system parameters
!> @param[in]     xv     Input compressed MRSF response vector
!> @param[out]    xv12   Output expanded MRSF response vector
!>
!> @date Aug 2024
!> @author Konstantin Komarov
  subroutine mrsfxvec(infos,xv,xv12)

    use precision, only: dp
    use types, only: information
    use messages, only: show_message, with_abort

    implicit none

    type(information), intent(in) :: infos

    real(kind=dp), intent(in), dimension(:) :: xv
    real(kind=dp), intent(inout), dimension(:) :: xv12

    integer :: noca, nocb, nbf, mrst
    integer :: i, ij, ijd, ijg, ijlr1, ijlr2, j, xvec_dim, ok
    real(kind=dp), parameter :: sqrt2 = 1.0_dp/sqrt(2.0_dp)
    real(kind=dp), allocatable, dimension(:) :: tmp

    nbf = infos%basis%nbf
    noca = infos%mol_prop%nelec_A
    nocb = infos%mol_prop%nelec_B
    mrst = infos%tddft%mult

    ijlr1 = (noca-1-nocb-1)*noca+noca-1
    ijg   = (noca-1-nocb-1)*noca+noca
    ijd   = (noca  -nocb-1)*noca+noca-1
    ijlr2 = (noca  -nocb-1)*noca+noca

    xvec_dim = noca*(nbf-nocb)

    allocate(tmp(xvec_dim), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', with_abort)

    if (mrst==1) then

      do i = 1, noca
        do j = nocb+1, nbf
          ij = (j-nocb-1)*noca+i
          if(ij==ijlr1) then
            tmp(ij) = xv(ijlr1)*sqrt2
            cycle
          else if(ij==ijlr2) then
            tmp(ij) = -xv(ijlr1)*sqrt2
            cycle
          end if
          tmp(ij) = xv(ij)
        end do
      end do

    else if (mrst==3) then

      do i = 1, noca
        do j = nocb+1,nbf
          ij = (j-nocb-1)*noca+i
          if(ij==ijlr1) then
            tmp(ij) = xv(ijlr1)*sqrt2
            cycle
          else if(ij==ijg) then
            tmp(ij) = 0.0_dp
            cycle
          else if(ij==ijd) then
            tmp(ij) = 0.0_dp
            cycle
          else if(ij==ijlr2) then
            tmp(ij) = xv(ijlr1)*sqrt2
            cycle
          end if
          tmp(ij) = xv(ij)
        end do
      end do

    end if

    xv12(:) = tmp(:)

    return

  end subroutine mrsfxvec

!>    @brief    Spin-pairing parts
!>              of singlet and triplet MRSF Lagrangian
!>
  subroutine mrsfsp(xhxa, xhxb, ca, cb, xv, fmrsf, noca, nocb)

    use precision, only: dp
    use messages, only: show_message, with_abort
    implicit none

    real(kind=dp), intent(out), dimension(:,:) :: xhxa, xhxb
    real(kind=dp), intent(in), dimension(:,:) :: ca, cb, xv
    real(kind=dp), intent(in), target, dimension(:,:,:) :: fmrsf
    integer, intent(in) :: noca, nocb

    integer :: nbf, i, j, lr1, lr2, ok

    real(kind=dp), allocatable :: scr(:,:), scr2(:,:)
    real(kind=dp), pointer, dimension(:,:) :: &
      adco1, adco2, ado1v, ado2v, aco12, ao21v

    ado2v => fmrsf(1,:,:)
    ado1v => fmrsf(2,:,:)
    adco1 => fmrsf(3,:,:)
    adco2 => fmrsf(4,:,:)
    ao21v => fmrsf(5,:,:)
    aco12 => fmrsf(6,:,:)

    nbf = ubound(ca, 1)
    lr1 = nocb+1
    lr2 = noca

    allocate(scr(nbf,nbf), &
             scr2(nbf,nbf), &
             source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', with_abort)
  ! Spin-pairing coupling contributions of xhxa

  ! o1v
    call dgemm('t', 'n', nbf, nbf, nbf, &
               1.0_dp, ca, nbf,  &
                       ao21v, nbf,  &
               0.0_dp, scr2, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, &
               2.0_dp, scr2, nbf, &
                       cb, nbf, &
               0.0_dp, scr, nbf)

    do j = noca+1, nbf
      xhxa(:,lr2) = xhxa(:,lr2)+scr(:,j)*xv(lr1,j)
      xhxa(:,lr1) = xhxa(:,lr1)-scr(:,j)*xv(lr2,j)
    end do

    ! co1
    call dgemm('t', 'n', nbf, nbf, nbf, &
              -1.0_dp, ca, nbf, &
                       aco12, nbf, &
               0.0_dp, scr2, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, &
               2.0_dp, scr2, nbf, &
                       cb, nbf, &
               0.0_dp, scr, nbf)

    do i = 1, nocb
      xhxa(:,i) = xhxa(:,i)+scr(:,lr2)*xv(i,lr1)
      xhxa(:,i) = xhxa(:,i)-scr(:,lr1)*xv(i,lr2)
    end do

    call dgemm('t', 'n', nbf, nbf, nbf, &
               1.0_dp, ca, nbf, &
                       adco2, nbf, &
               0.0_dp, scr2, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, &
               2.0_dp, scr2, nbf, &
                       cb, nbf, &
               0.0_dp, scr, nbf)

    do j = noca+1, nbf
      xhxa(:,lr1) = xhxa(:,lr1)+scr(:,j)*xv(lr1,j)
    end do

    call dgemm('t', 'n', nbf, nbf, nbf, &
               1.0_dp, ca, nbf, &
                       adco1, nbf, &
               0.0_dp, scr2, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, &
               2.0_dp, scr2, nbf, &
                       cb, nbf, &
               0.0_dp, scr, nbf)

    do j = noca+1, nbf
      xhxa(:,lr2) = xhxa(:,lr2)+scr(:,j)*xv(lr2,j)
    end do

    call dgemm('t', 'n', nbf, nbf, nbf, &
               1.0_dp, ca, nbf, &
                       ado2v, nbf, &
               0.0_dp, scr2, nbf)
    call dgemm('n', 'n', nbf, 1, nbf, &
               2.0_dp, scr2, nbf, &
                       cb(:,lr1), nbf, &
               0.0_dp, scr, nbf)

    do i = 1, nocb
      xhxa(:,i) = xhxa(:,i)+scr(:,1)*xv(i,lr1)
    end do

  ! co2
    call dgemm('t', 'n', nbf, nbf, nbf, &
               1.0_dp, ca, nbf, &
                       ado1v, nbf, &
               0.0_dp, scr2, nbf)
    call dgemm('n', 'n', nbf, 1, nbf, &
               2.0_dp, scr2, nbf, &
                       cb(:,lr2), nbf, &
               0.0_dp, scr, nbf)

    do i = 1, nocb
      xhxa(:,i) = xhxa(:,i)+scr(:,1)*xv(i,lr2)
    end do

   ! Spin-pairing coupling contributions of xhxb

    call dgemm('t', 'n', nbf, nbf, nbf, &
               1.0_dp, ca, nbf, &
                       ao21v, nbf,&
               0.0_dp, scr2, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, &
               2.0_dp, scr2, nbf, &
                       cb, nbf, &
               0.0_dp, scr, nbf)

    do j = noca+1, nbf
      xhxb(:,j) = xhxb(:,j)+scr(lr2,:)*xv(lr1,j)
    end do

    call dgemm('t', 'n', nbf, nbf, nbf, &
              -1.0_dp, ca, nbf, &
                       ao21v, nbf, &
               0.0_dp, scr2, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, &
               2.0_dp, scr2, nbf, &
                       cb, nbf, &
               0.0_dp, scr, nbf)

    do j = noca+1, nbf
      xhxb(:,j) = xhxb(:,j)+scr(lr1,:)*xv(lr2,j)
    end do

  ! co1
    call dgemm('t', 'n', nbf, nbf, nbf, &
              -1.0_dp, ca, nbf, &
                       aco12, nbf, &
               0.0_dp, scr2, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, &
               2.0_dp, scr2, nbf, &
                       cb, nbf, &
               0.0_dp, scr, nbf)

    do i = 1, nocb
      xhxb(:,lr2) = xhxb(:,lr2)+scr(i,:)*xv(i,lr1)
    end do

    call dgemm('t', 'n', nbf, nbf, nbf, &
                         1.0_dp, ca, nbf, &
                                 aco12, nbf, &
                         0.0_dp, scr2, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, &
                         2.0_dp, scr2, nbf, &
                                 cb, nbf, &
                         0.0_dp, scr, nbf)

    do i = 1, nocb
      xhxb(:,lr1) = xhxb(:,lr1)+scr(i,:)*xv(i,lr2)
    end do

    call dgemm('t', 'n', nbf, nbf, nbf, &
               1.0_dp, ca, nbf, &
                       ado2v, nbf, &
               0.0_dp, scr2, nbf)
    call dgemm('n' ,'n', nbf, nbf, nbf, &
               2.0_dp, scr2, nbf, &
                       cb, nbf, &
               0.0_dp, scr, nbf)

    do i = 1, nocb
      xhxb(:,lr1) = xhxb(:,lr1)+scr(i,:)*xv(i,lr1)
    end do

    call dgemm('t', 'n', nbf, nbf, nbf, &
               1.0_dp, ca, nbf, &
                       ado1v, nbf, &
               0.0_dp, scr2, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, &
               2.0_dp, scr2, nbf, &
                       cb, nbf, &
               0.0_dp, scr, nbf)

    do i = 1, nocb
      xhxb(:,lr2) = xhxb(:,lr2)+scr(i,:)*xv(i,lr2)
    end do

  ! O1V
    call dgemm('t', 'n', 1, nbf, nbf, &
               1.0_dp, ca(:, lr1), nbf, &
                       adco2, nbf,  &
               0.0_dp, scr2, 1)
    call dgemm('n', 'n', 1, nbf, nbf, &
               2.0_dp, scr2, 1, &
                       cb, nbf, &
               0.0_dp, scr, 1)

    do j = noca+1, nbf
      xhxb(:,j) = xhxb(:,j)+scr(:,1)*xv(lr1,j)
    end do

    call dgemm('t', 'n', 1, nbf, nbf, &
               1.0_dp, ca(:, lr2), nbf, &
                       adco1, nbf,  &
               0.0_dp, scr2, 1)
    call dgemm('n',  'n', 1, nbf, nbf, &
               2.0_dp, scr2, 1, &
                       cb, nbf, &
               0.0_dp, scr, 1)

    do j = noca+1, nbf
      xhxb(:,j) = xhxb(:,j)+scr(:,1)*xv(lr2,j)
    end do

    return

  end subroutine mrsfsp

  subroutine mrsfrowcal(wmo, mo_energy_a, fa, fb, xk, &
                        xhxa, xhxb, hppija, hppijb, noca, nocb)

    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:,:) :: wmo
    real(kind=dp), intent(in), dimension(:) :: mo_energy_a
    real(kind=dp), intent(in), dimension(:,:) :: fa, fb
    real(kind=dp), intent(in), dimension(:) :: xk
    real(kind=dp), intent(in), dimension(:,:) :: xhxa, xhxb
    real(kind=dp), intent(in), dimension(:,:) :: hppija, hppijb
    integer, intent(in) :: noca, nocb

    real(kind=dp), allocatable, dimension(:,:) :: wrk, scr
    integer :: i, a, k, x, y, j, b, ij, nbf, lr1, lr2


    nbf = ubound(fa, 1)
    lr1 = nocb+1
    lr2 = noca

    allocate(wrk(nbf,nbf), &
             scr(nbf,nbf), &
             source=0.0_dp)

  ! Unpack  xk
    ij = 0
    do i = lr1, lr2
      do j = 1, nocb
        ij = ij+1
        scr(j,i) = xk(ij)
      end do
    end do

    do i = noca+1, nbf
      do j = 1, nocb
        ij = ij+1
        scr(j,i) = xk(ij)
      end do
    end do

    do k = noca+1, nbf
      do i = lr1, lr2
        ij = ij+1
        scr(i,k) = xk(ij)
      end do
    end do

  ! W_ix
    do x = 1, nocb
      do k = 1, nocb
        wrk(x,1:2) = wrk(x,1:2)-fa(k,x)*scr(k,lr1:lr2)
      end do
    end do
    do x = 1, nocb
      do k = 1, nbf-noca
        wrk(x,1:2) = wrk(x,1:2)+scr(lr1:lr2,noca+k)*fa(noca+k,x)
      end do
    end do

    wmo(1:nocb,lr1:lr2) = wrk(1:nocb,1:2)*0.5_dp &
                        + xhxa(1:nocb,lr1:lr2) &
                        + xhxb(1:nocb,lr1:lr2) &
                        + hppija(1:nocb,lr1:lr2)
    wmo(1:nocb,lr1) = wmo(1:nocb,lr1) &
                    + mo_energy_a(1:nocb)*scr(1:nocb,lr1)
    wmo(1:nocb,lr2) = wmo(1:nocb,lr2) &
                    + mo_energy_a(1:nocb)*scr(1:nocb,lr2)

!   ----- W_IA -----
    wrk = 0.0_dp
    do i = 1, nocb
      do a = 1, nbf-noca
        wrk(i,a) = wrk(i,a)+fa(lr1,i)*scr(lr1,noca+a) &
                           +fa(lr2,i)*scr(lr2,noca+a)
      end do
    end do

    wmo(1:nocb,noca+1:nbf) = wrk(1:nocb,1:nbf-noca)*0.5_dp &
                          + xhxb(1:nocb,noca+1:nbf)

    do a = 1, nbf-noca
      wmo(1:nocb,noca+a) = wmo(1:nocb,noca+a) &
                         + mo_energy_a(1:nocb)*scr(1:nocb,noca+a)
    end do

!   ----- W_XA -----
    wrk = 0.0_dp
    do a = 1, nbf-noca
      do k = 1, nocb
        wrk(1,a) = wrk(1,a)+fa(k,lr1)*scr(k,noca+a)
        wrk(2,a) = wrk(2,a)+fa(k,lr2)*scr(k,noca+a)
      end do
    end do

    do a = 1, nbf-noca
      wrk(1,a) = wrk(1,a)-fb(lr1,lr1)*scr(lr1,noca+a) &
                         -fb(lr2,lr1)*scr(lr2,noca+a)
      wrk(2,a) = wrk(2,a)-fb(lr1,lr2)*scr(lr1,noca+a) &
                         -fb(lr2,lr2)*scr(lr2,noca+a)
    end do

    wmo(lr1:lr2,noca+1:nbf) = wrk(1:2,1:nbf-noca)*0.5_dp &
                           + xhxb(lr1:lr2,noca+1:nbf)
    do a = noca+1, nbf
      wmo(lr1:lr2,a) = wmo(lr1:lr2,a) &
                     + mo_energy_a(lr1:lr2)*scr(lr1:lr2,a)
    end do

  !  W_ij
    do i = 1, nocb
      do j = 1, i
        wmo(i,j) = hppija(i,j)+hppijb(i,j)+xhxa(j,i)
      end do
    end do

  ! W_xy
    do x = nocb+1, noca
      do y = nocb+1, x
        wmo(x,y) = xhxa(y,x)+xhxb(y,x)+hppija(x,y)
      end do
    end do

  ! W_ab
    do a = noca+1, nbf
      do b = noca+1, a
        wmo(a,b) = xhxb(b,a)
      end do
    end do

  ! Scale diagonal elements
    do i = 1, nbf
      wmo(i,i) = wmo(i,i)*0.5_dp
    end do

    wmo = -wmo

    return

  end subroutine mrsfrowcal

  subroutine mrsfqrorhs(rhs, xhxa, xhxb, hpta, hptb, tab, tij, fa, fb, noca,  &
                        nocb)

    use precision, only: dp
    implicit none

    real(kind=dp), intent(out), dimension(:) :: rhs
    real(kind=dp), intent(inout), dimension(:,:) :: xhxa, xhxb
    real(kind=dp), intent(in), dimension(:,:) :: hpta
    real(kind=dp), intent(in), dimension(:,:) :: hptb
    real(kind=dp), intent(in), dimension(:,:) :: tij
    real(kind=dp), intent(in), dimension(:,:) :: tab
    real(kind=dp), intent(in), dimension(:,:) :: fa, fb
    integer, intent(in) :: noca, nocb

    real(kind=dp), allocatable, dimension(:,:) :: scr
    integer :: nbf, i, j, ij, a, x, nconf

    nbf = ubound(fa, 1)

    allocate(scr(nbf,nbf), &
             source=0.0_dp)

  ! Alpha
  ! hxa+= 2*fa(p+,a+)*ta(a+,b+)
    do j = noca+1, nbf
      do i = noca+1, nbf
        scr(i,j) = tab(i-noca,j-noca)
      end do
    end do
    call dgemm('n', 'n', nbf, nbf, nbf, &
               2.0_dp, fa, nbf, &
                       scr, nbf, &
               1.0_dp, xhxa, nbf)

  ! Beta
  ! xhxb+= 2*fb(p-, i-)*tb(i-, j-)
    call dgemm('n', 'n', nbf, nocb, nocb, &
               2.0_dp, fb, nbf, &
                       tij, nocb, &
               1.0_dp, xhxb, nbf)

    rhs = 0.0_dp

  ! doc-socc
    ij = 0
    do x = nocb+1, noca
      do i = 1, nocb
        ij = ij+1
        rhs(ij) = hptb(i,x-nocb)+xhxb(x,i)
      end do
    end do

  ! doc-virt
    do a = noca+1, nbf
      do i = 1, nocb
        ij = ij+1
        rhs(ij) = hpta(i,a-noca)+hptb(i,a-nocb) &
                + xhxb(a,i)-xhxa(i,a)
      end do
    end do

  ! soc-virt
    do a = noca+1, nbf
      do x = nocb+1, noca
        ij = ij+1
        rhs(ij) = hpta(x,a-noca)-xhxa(x,a)
      end do
    end do

    nconf = ij
    rhs(1:nconf) = -rhs(1:nconf)

    return

  end subroutine mrsfqrorhs

  subroutine mrsfqropcal(pa, pb, tab, tij, z, noca, nocb)

    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:,:) :: pa, pb
    real(kind=dp), intent(in), dimension(:,:) :: tab, tij
    real(kind=dp), intent(in), dimension(:) :: z
    integer, intent(in) :: noca, nocb

    integer :: nbf, i, j, a, x, ij

    nbf = ubound(pa, 1)

  ! alpha
    pa = 0.0_dp
    do j=noca+1, nbf
      do i=noca+1, nbf
        pa(i, j) = tab(i-noca, j-noca)
      end do
    end do
  ! beta
    pb = 0.0_dp
    do j=1, nocb
      do i=1, nocb
        pb(i, j) = tij(i, j)
      end do
    end do

  ! doc-socc
    ij = 0
    do x = nocb+1, noca
      do i = 1, nocb
        ij = ij+1
        pb(i,x) = pb(i,x)+z(ij)*0.5_dp
      end do
    end do

  ! doc-virt
    do a = noca+1, nbf
      do i = 1, nocb
        ij = ij + 1
        pa(i,a) = pa(i,a)+z(ij)*0.5_dp
        pb(i,a) = pb(i,a)+z(ij)*0.5_dp
      end do
    end do

  ! socc-virt
    do a = noca+1, nbf
      do x = nocb+1, noca
        ij = ij+1
        pa(x,a) = pa(x,a)+z(ij)*0.5_dp
      end do
    end do

    return

  end subroutine mrsfqropcal

  subroutine mrsfqrowcal(w, mo_energy_a, fa, fb, z,  &
                         xhxa, xhxb, hppija, hppijb, noca, nocb)

    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:,:) :: w
    real(kind=dp), intent(in), dimension(:) :: mo_energy_a
    real(kind=dp), intent(in), dimension(:,:) :: fa, fb
    real(kind=dp), intent(in), dimension(:) :: z
    real(kind=dp), intent(in), dimension(:,:) :: xhxa, xhxb
    real(kind=dp), intent(in), dimension(:,:) :: hppija, hppijb
    integer, intent(in) :: noca, nocb

    real(kind=dp), allocatable, dimension(:,:) :: scr, wrk
    integer :: i, a, k, x, y, j, b, ij, nbf, lr1, lr2

    nbf = ubound(fa, 1)
    lr1 = nocb+1
    lr2 = noca

    allocate(wrk(nbf,nbf), &
             scr(nbf,nbf), &
             source=0.0_dp)

    ij = 0
    do x = nocb+1, noca
      do i = 1, nocb
        ij = ij+1
        scr(i, x) = z(ij)
      end do
    end do

    do a = noca+1, nbf
      do i = 1, nocb
        ij = ij+1
        scr(i, a) = z(ij)
      end do
    end do

    do a = noca+1, nbf
      do x = nocb+1, noca
        ij = ij+1
        scr(x, a) = z(ij)
      end do
    end do

  ! w_ix
    do i = 1, nocb
      do k = 1, nocb
        wrk(i,1:2) = wrk(i,1:2)-fa(k,i)*scr(k,lr1:lr2)
      end do
    end do

    do i = 1, nocb
      do y = 1, nbf-noca
        wrk(i,1:2) = wrk(i,1:2)+fa(noca+y,i)*scr(lr1:lr2,noca+y)
      end do
    end do

    w(1:nocb,lr1:lr2) = 0.5_dp*wrk(1:nocb,1:2)+hppija(1:nocb,lr1:lr2)

    do i = 1,nocb
      w(i,lr1:lr2) = w(i,lr1:lr2) &
                   + mo_energy_a(i)*scr(i,lr1:lr2)
    end do

  ! w_ia
    wrk = 0.0_dp
    do a = 1, nbf-noca
      do i = 1, nocb
        wrk(i,a) = wrk(i,a)+fa(lr1,i)*scr(lr1,noca+a)
        wrk(i,a) = wrk(i,a)+fa(lr2,i)*scr(lr2,noca+a)
      end do
    end do

    do a = 1, nbf-noca
      w(1:nocb,noca+a) = mo_energy_a(1:nocb)*scr(1:nocb,noca+a) &
                       + 0.5_dp*wrk(1:nocb,a) &
                       + xhxa(1:nocb,noca+a)
    end do

  ! w_xa
    wrk = 0.0_dp
    do a = 1, nbf-noca
      do k = 1, nocb
        wrk(1:2,a) = wrk(1:2,a)+fa(k, lr1:lr2)*scr(k, noca+a)
      end do

    end do
    do a = 1, nbf-noca
       wrk(1:2,a) = wrk(1:2,a)-fb(lr1, lr1:lr2)*scr(lr1, noca+a)
       wrk(1:2,a) = wrk(1:2,a)-fb(lr2, lr1:lr2)*scr(lr2, noca+a)
    end do

    do a = 1, nbf-noca
      w(lr1:lr2, noca+a) = mo_energy_a(lr1:lr2)*scr(lr1:lr2, noca+a) &
                         + 0.5_dp*wrk(1:2,a) &
                         + xhxa(lr1:lr2, noca+a)
    end do

  ! w_ij
    do i = 1, nocb
      do j = 1, i
        w(i, j) = hppija(i, j)+hppijb(i, j)+xhxb(j, i)
      end do
    end do

  ! w_xy
    do x = nocb+1, noca
      do y = nocb+1, x
        w(x, y) = hppija(x, y)
      end do
    end do

  ! w_ab
    do a = noca+1, nbf
      do b = noca+1, a
        w(a, b) = xhxa(b, a)
      end do
    end do

  ! scale diagonal elements
    do i = 1, nbf
      w(i, i) = 0.5_dp*w(i, i)
    end do

    w = -w

    return

  end subroutine mrsfqrowcal

  subroutine get_mrsf_transition_density(infos, trden, bvec_mo, ist, jst)

    use precision, only: dp
    use messages, only: show_message, with_abort
    use types, only: information

    implicit none

    type(information), intent(in) :: infos
    real(kind=dp), intent(out), dimension(:,:) :: trden
    real(kind=dp), intent(in), dimension(:,:) :: bvec_mo
    integer, intent(in) :: ist, jst

    real(kind=dp), allocatable :: xv12i(:,:), xv12j(:,:)
    integer :: nbf, noca, nocb, nvirb, ok
    integer :: i, ij, ijd, ijg, ijlr1, ijlr2, j, mrst
    real(kind=dp), parameter :: rsqrt = 1.0_dp/sqrt(2.0_dp)

    nbf = infos%basis%nbf
    noca = infos%mol_prop%nelec_a
    mrst = infos%tddft%mult
    nocb = noca-2
    nvirb = nbf-nocb

    allocate(xv12i(noca,nvirb), &
             xv12j(noca,nvirb), &
             source=0.0_dp, stat=ok)

    if (ok /= 0) call show_message('Cannot allocate memory', with_abort)

    ijlr1 = (noca-1-nocb-1)*noca+noca-1
    ijg   = (noca-1-nocb-1)*noca+noca
    ijd   = (noca-nocb-1)*noca+noca-1
    ijlr2 = (noca-nocb-1)*noca+noca

    if (mrst==1) then

      do i = 1, noca
        do j = nocb + 1, nbf
          ij = (j-nocb-1)*noca+i
          if (ij==ijlr1) then
            xv12j(i,j-nocb) = bvec_mo(ijlr1,jst)*rsqrt
            xv12i(i,j-nocb) = bvec_mo(ijlr1,ist)*rsqrt
            cycle
          else if (ij==ijlr2) then
            xv12j(i,j-nocb) = -bvec_mo(ijlr1,jst)*rsqrt
            xv12i(i,j-nocb) = -bvec_mo(ijlr1,ist)*rsqrt
            cycle
          end if
          xv12j(i,j-nocb) = bvec_mo(ij,jst)
          xv12i(i,j-nocb) = bvec_mo(ij,ist)
        end do
      end do

    else if (mrst==3) then

      do i = 1, noca
        do j = nocb+1, nbf
          ij = (j-nocb-1)*noca+i
          if (ij==ijlr1) then
            xv12j(i,j-nocb) = bvec_mo(ijlr1,jst)*rsqrt
            xv12i(i,j-nocb) = bvec_mo(ijlr1,ist)*rsqrt
            cycle
          else if (ij==ijlr2) then
            xv12j(i,j-nocb) = bvec_mo(ijlr1,jst)*rsqrt
            xv12i(i,j-nocb) = bvec_mo(ijlr1,ist)*rsqrt
            cycle
          else if (ij==ijg) then
            xv12j(i,j-nocb) = 0.0_dp
            xv12i(i,j-nocb) = 0.0_dp
            cycle
          else if (ij==ijd) then
            xv12j(i,j-nocb) = 0.0_dp
            xv12i(i,j-nocb) = 0.0_dp
            cycle
          end if
          xv12j(i,j-nocb) = bvec_mo(ij,jst)
          xv12i(i,j-nocb) = bvec_mo(ij,ist)
        end do
      end do

    end if

    call get_trans_den(trden, xv12i, xv12j, noca, nocb, nvirb)

    return

  end subroutine get_mrsf_transition_density

  subroutine get_trans_den(trden, xv12i, xv12j, noca, nocb, nvirb)

    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:,:) :: trden
    real(kind=dp), intent(in), dimension(:,:) :: xv12i, xv12j
    integer, intent(in) :: noca, nocb, nvirb

    logical :: ioo, joo
    integer :: i, j, k
    real(kind=dp) :: tmp, scal
    real(kind=dp), parameter :: sqrt2 = sqrt(2.0_dp)

    ! Revised trden(vir/vir)
    do i = 1, nvirb
      do j = 1, nvirb
        tmp = 0.0_dp
        do k = 1, noca
          ioo = .false.
          joo = .false.
          if ((i==1 .or. i==2) .and. (k==noca-1 .or. k==noca)) ioo = .true.
          if ((j==1 .or. j==2) .and. (k==noca-1 .or. k==noca)) joo = .true.
          scal = 1.0_dp
          if (ioo .and. .not. joo) scal = sqrt2
          if (.not. ioo .and. joo) scal = sqrt2

          tmp = tmp+scal*xv12j(k,i)*xv12i(k,j)
        end do
        trden(i+nocb,j+nocb) = trden(i+nocb,j+nocb)+tmp
      end do
    end do

    ! Revised trden(occ/occ)
    do i = 1, noca
      do j = 1, noca
        tmp = 0.0_dp
        do k = 1, nvirb
          ioo = .false.
          joo = .false.
          if ((i==noca-1 .or. i==noca) .and. (k==1 .or. k==2)) ioo = .true.
          if ((j==noca-1 .or. j==noca) .and. (k==1 .or. k==2)) joo = .true.
          scal = 1.0_dp
          if (ioo .and. .not. joo) scal = sqrt2
          if (.not. ioo .and. joo) scal = sqrt2

          tmp = tmp + scal*xv12j(j,k)*xv12i(i,k)
        end do
        trden(i,j) = trden(i,j)-tmp
      end do
    end do

    ! Trden(occ/vir) elements are zero
    return

  end subroutine

end module tdhf_mrsf_lib
