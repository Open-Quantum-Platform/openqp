module tdhf_lib

    use precision, only : dp
    use int2_compute, only: int2_fock_data_t, int2_storage_t
    use basis_tools, only: basis_set
    use oqp_linalg

    type, extends(int2_fock_data_t) :: int2_td_data_t
      real(kind=dp), allocatable :: apb(:,:,:,:)
      real(kind=dp), allocatable :: amb(:,:,:,:)
      real(kind=dp), pointer :: d2(:,:,:) => null()
      logical :: int_apb = .true. !< do A+B part
      logical :: int_amb = .false. !< do A-B part, needed for hybrid functionals
      logical :: tamm_dancoff = .false. !< Tamm-Dancoff approximation
      logical :: tamm_dancoff_coulomb = .false. !< Whetehr to includ Coulomb terms with TDA, not needed in SFDFT
    contains
        procedure :: parallel_start => int2_td_data_t_parallel_start
        procedure :: parallel_stop => int2_td_data_t_parallel_stop
        procedure :: init_screen => int2_td_data_t_init_screen
        procedure :: update => int2_td_data_t_update
        procedure :: clean => int2_td_data_t_clean
    end type

    type, extends(int2_td_data_t) :: int2_tdgrd_data_t
    contains
      procedure :: update => int2_tdgrd_data_t_update
    end type

    type, extends(int2_fock_data_t) :: int2_rpagrd_data_t
      real(kind=dp), allocatable :: hpp(:,:,:,:,:)  ! H+[X+Y]
      real(kind=dp), allocatable :: hpt(:,:,:,:,:)  ! H+[T]
      real(kind=dp), allocatable :: hmm(:,:,:,:,:)  ! H-[X-Y]
      real(kind=dp), pointer :: xpy(:,:,:,:) => null() ! X+Y
      real(kind=dp), pointer :: xmy(:,:,:,:) => null() ! X-Y
      real(kind=dp), pointer :: t(:,:,:,:) => null() ! T
      integer :: np = 0
      integer :: nm = 0
      integer :: nt = 0
      integer :: nspin = 0
      integer :: nbf = 0
      logical :: tamm_dancoff = .false. !< Tamm-Dancoff approximation
    contains
        procedure :: parallel_start => int2_rpagrd_data_t_parallel_start
        procedure :: parallel_stop  => int2_rpagrd_data_t_parallel_stop
        procedure :: init_screen    => int2_rpagrd_data_t_init_screen
        procedure :: update         => int2_rpagrd_data_t_update
        procedure :: clean          => int2_rpagrd_data_t_clean

        procedure, non_overridable :: hplus  => int2_rpagrd_data_t_update_hplus
        procedure, non_overridable :: hminus => int2_rpagrd_data_t_update_hminus
    end type
contains

!###############################################################################

  subroutine int2_td_data_t_parallel_start(this, basis, nthreads)
    implicit none
    class(int2_td_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nthreads
    integer :: nbf, nsh

    nbf = basis%nbf
    this%fockdim = nbf*(nbf+1) / 2
    this%nfocks = ubound(this%d2,size(shape(this%d2)))
    this%nthreads = nthreads
    nsh = basis%nshell

    if (this%cur_pass == 1) then
      if (allocated(this%apb)) deallocate(this%apb)
      if (allocated(this%amb)) deallocate(this%amb)
      if (allocated(this%dsh)) deallocate(this%dsh)

      allocate(this%apb(nbf, nbf, this%nfocks, nthreads), &
               this%amb(nbf, nbf, this%nfocks, nthreads), &
               this%dsh(nsh,nsh), &
               source=0.0d0)
    end if

    call this%init_screen(basis)

  end subroutine

!###############################################################################

  subroutine int2_td_data_t_parallel_stop(this)
    use mathlib, only: symmetrize_matrix
    implicit none
    integer :: flast, amblast, nbf, i
    class(int2_td_data_t), intent(inout) :: this
    if (this%cur_pass /= this%num_passes) return
    flast  = size(shape(this%apb))
    amblast = size(shape(this%amb))
    nbf = ubound(this%amb, 1)
    if (this%nthreads /= 1) then
      this%apb(:,:,:,lbound(this%apb, flast)) = sum(this%apb, dim=flast)
      this%amb(:,:,:,lbound(this%amb, amblast)) = sum(this%amb, dim=amblast)
    end if
    call this%pe%allreduce(this%apb(:,:,:,1), &
                       size(this%apb(:,:,:,1)))
    call this%pe%allreduce(this%amb(:,:,:,1), &
                         size(this%amb(:,:,:,1)))

    do i = lbound(this%apb,3), ubound(this%apb,3)
      call symmetrize_matrix(this%apb(:,:,i,1), nbf)
    end do
    this%nthreads = 1
  end subroutine

!###############################################################################

  subroutine int2_td_data_t_clean(this)
    implicit none
    class(int2_td_data_t), intent(inout) :: this
    deallocate(this%apb)
    deallocate(this%amb)
    deallocate(this%dsh)
    nullify(this%d)
    nullify(this%d2)
  end subroutine

!###############################################################################

  subroutine int2_td_data_t_init_screen(this, basis)
    implicit none
    class(int2_td_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis

!   Form shell density
    call shltd(this%dsh, this%d2, basis)
    this%max_den = maxval(abs(this%dsh))

  end subroutine

!###############################################################################

  subroutine int2_td_data_t_update(this, buf)
    implicit none
    class(int2_td_data_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    integer :: i, j, k, l, n
    real(kind=dp) :: xval1, cval2, val2c, cval4, &
                     val, val1, val4c
    integer :: ifock, mythread

    xval1 = this%scale_exchange
    cval2 = 2 * this%scale_coulomb
    cval4 = 4 * this%scale_coulomb

    mythread = buf%thread_id

    associate (&
                apb => this%apb(:,:,:,mythread), &
                amb => this%amb(:,:,:,mythread), &
                d2 => this%d2 &
                )

      do ifock = 1, this%nfocks
        do n = 1, buf%ncur
          i  = buf%ids(1,n)
          j  = buf%ids(2,n)
          k  = buf%ids(3,n)
          l  = buf%ids(4,n)
          val = buf%ints(n)

          if (this%tamm_dancoff) then
            val1 = val*xval1
            val2c = val*cval2
            !A
            amb(i,k,ifock) = amb(i,k,ifock) - val1 * d2(j,l,ifock)
            amb(k,i,ifock) = amb(k,i,ifock) - val1 * d2(l,j,ifock)
            amb(i,l,ifock) = amb(i,l,ifock) - val1 * d2(j,k,ifock)
            amb(l,i,ifock) = amb(l,i,ifock) - val1 * d2(k,j,ifock)
            amb(j,k,ifock) = amb(j,k,ifock) - val1 * d2(i,l,ifock)
            amb(k,j,ifock) = amb(k,j,ifock) - val1 * d2(l,i,ifock)
            amb(j,l,ifock) = amb(j,l,ifock) - val1 * d2(i,k,ifock)
            amb(l,j,ifock) = amb(l,j,ifock) - val1 * d2(k,i,ifock)
            if (this%tamm_dancoff_coulomb) then
              amb(i,j,ifock) = amb(i,j,ifock) + val2c * (d2(k,l,ifock)+d2(l,k,ifock))
              amb(j,i,ifock) = amb(j,i,ifock) + val2c * (d2(k,l,ifock)+d2(l,k,ifock))
              amb(k,l,ifock) = amb(k,l,ifock) + val2c * (d2(i,j,ifock)+d2(j,i,ifock))
              amb(l,k,ifock) = amb(l,k,ifock) + val2c * (d2(i,j,ifock)+d2(j,i,ifock))
            end if
          else
            val1 = val*xval1
            val4c = val*cval4

            if (this%int_apb) then
              ! A+B
              ! Coulomb
              apb(i,j,ifock) = apb(i,j,ifock) + val4c * (d2(k,l,ifock)+d2(l,k,ifock))
              apb(k,l,ifock) = apb(k,l,ifock) + val4c * (d2(i,j,ifock)+d2(j,i,ifock))

              ! Exchange
              apb(i,k,ifock) = apb(i,k,ifock) - val1 * (d2(j,l,ifock)+d2(l,j,ifock))
              apb(i,l,ifock) = apb(i,l,ifock) - val1 * (d2(j,k,ifock)+d2(k,j,ifock))
              apb(j,k,ifock) = apb(j,k,ifock) - val1 * (d2(i,l,ifock)+d2(l,i,ifock))
              apb(j,l,ifock) = apb(j,l,ifock) - val1 * (d2(i,k,ifock)+d2(k,i,ifock))
            end if

            if (this%int_amb) then
              ! A-B
              amb(i,k,ifock) = amb(i,k,ifock) + val1 * (d2(l,j,ifock)-d2(j,l,ifock))
              amb(i,l,ifock) = amb(i,l,ifock) + val1 * (d2(k,j,ifock)-d2(j,k,ifock))
              amb(j,k,ifock) = amb(j,k,ifock) + val1 * (d2(l,i,ifock)-d2(i,l,ifock))
              amb(j,l,ifock) = amb(j,l,ifock) + val1 * (d2(k,i,ifock)-d2(i,k,ifock))

              amb(k,i,ifock) = amb(k,i,ifock) - val1 * (d2(l,j,ifock)-d2(j,l,ifock))
              amb(l,i,ifock) = amb(l,i,ifock) - val1 * (d2(k,j,ifock)-d2(j,k,ifock))
              amb(k,j,ifock) = amb(k,j,ifock) - val1 * (d2(l,i,ifock)-d2(i,l,ifock))
              amb(l,j,ifock) = amb(l,j,ifock) - val1 * (d2(k,i,ifock)-d2(i,k,ifock))
            end if
          end if
        end do
      end do

    end associate

    buf%ncur = 0

  end subroutine

!###############################################################################

  subroutine int2_tdgrd_data_t_update(this, buf)
    implicit none
    class(int2_tdgrd_data_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    integer :: i, j, k, l, n
    real(kind=dp) :: xval1, xval2, &
                     val, val1, val2
    integer :: mythread

    xval1 = 1 * this%scale_exchange
    xval2 = 2 * this%scale_coulomb
    mythread = buf%thread_id

    associate (&
                apb => this%apb(:,:,:,mythread), &
                amb => this%amb(:,:,:,mythread), &
                d2 => this%d2 &
                )

      do n = 1, buf%ncur
        i  = buf%ids(1,n)
        j  = buf%ids(2,n)
        k  = buf%ids(3,n)
        l  = buf%ids(4,n)
        val = buf%ints(n)

        val1 = val*xval1
        val2 = val*xval2

        if (this%int_apb) then
          ! A+B
          ! Coulomb
          apb(i,j,1) = apb(i,j,1)+val2*(d2(k,l,1)+d2(l,k,1)+d2(k,l,2)+d2(l,k,2))
          apb(k,l,1) = apb(k,l,1)+val2*(d2(i,j,1)+d2(j,i,1)+d2(i,j,2)+d2(j,i,2))

          apb(i,j,2) = apb(i,j,2)+val2*(d2(k,l,1)+d2(l,k,1)+d2(k,l,2)+d2(l,k,2))
          apb(k,l,2) = apb(k,l,2)+val2*(d2(i,j,1)+d2(j,i,1)+d2(i,j,2)+d2(j,i,2))

!         ! Exchange
          apb(i,k,1) = apb(i,k,1)-val1*(d2(j,l,1)+d2(l,j,1))
          apb(i,l,1) = apb(i,l,1)-val1*(d2(j,k,1)+d2(k,j,1))
          apb(j,k,1) = apb(j,k,1)-val1*(d2(i,l,1)+d2(l,i,1))
          apb(j,l,1) = apb(j,l,1)-val1*(d2(i,k,1)+d2(k,i,1))
!
          apb(i,k,2) = apb(i,k,2)-val1*(d2(j,l,2)+d2(l,j,2))
          apb(i,l,2) = apb(i,l,2)-val1*(d2(j,k,2)+d2(k,j,2))
          apb(j,k,2) = apb(j,k,2)-val1*(d2(i,l,2)+d2(l,i,2))
          apb(j,l,2) = apb(j,l,2)-val1*(d2(i,k,2)+d2(k,i,2))
        end if

        if (this%int_amb) then
          ! A-B
          amb(i,k,1) = amb(i,k,1)+val1*(d2(l,j,1)-d2(j,l,1))
          amb(i,l,1) = amb(i,l,1)+val1*(d2(k,j,1)-d2(j,k,1))
          amb(j,k,1) = amb(j,k,1)+val1*(d2(l,i,1)-d2(i,l,1))
          amb(j,l,1) = amb(j,l,1)+val1*(d2(k,i,1)-d2(i,k,1))
          amb(k,i,1) = amb(k,i,1)-val1*(d2(l,j,1)-d2(j,l,1))
          amb(l,i,1) = amb(l,i,1)-val1*(d2(k,j,1)-d2(j,k,1))
          amb(k,j,1) = amb(k,j,1)-val1*(d2(l,i,1)-d2(i,l,1))
          amb(l,j,1) = amb(l,j,1)-val1*(d2(k,i,1)-d2(i,k,1))
        end if
      end do

    end associate

    buf%ncur = 0

  end subroutine

!###############################################################################
!###############################################################################

  subroutine shltd(dsh,da,basis)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set

    implicit none

    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(out) :: dsh(:,:)
    real(kind=dp), intent(in), dimension(:,:,:) :: da

    integer :: ish, jsh, maxi, maxj, mini, &
               minj

!   RHF
    do ish = 1, basis%nshell
      mini = basis%ao_offset(ish)
      maxi = mini + basis%naos(ish) - 1
      do jsh = 1, ish
        minj = basis%ao_offset(jsh)
        maxj = minj + basis%naos(jsh) - 1
        dsh(ish,jsh) = maxval(abs(da(minj:maxj,mini:maxi,:)))
        dsh(jsh,ish) = dsh(ish,jsh)
      end do
    end do
  end subroutine shltd


  subroutine inivec(eiga,eigb,bvec_mo,xm, &
                    nocca,noccb,nvec)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(in), dimension(:) :: eiga, eigb
    real(kind=dp), intent(out), dimension(:,:) :: bvec_mo
    real(kind=dp), intent(out), dimension(:) :: xm
    integer, intent(in) :: nocca, noccb
    integer, intent(in) :: nvec

    integer :: i, ij, j, k, nbf, mxvec

    integer :: itmp(nvec)
    real(kind=dp) :: xtmp(nvec)

    nbf = ubound(eiga, 1)
    mxvec = ubound(bvec_mo, 2)

  ! -- Set xm(xvec_dim)
    do j = noccb+1, nbf
      do i = 1, nocca
        ij = (j-noccb-1)*nocca+i
        xm(ij) = eigb(j)-eiga(i)
      end do
    end do

!   Find indices of the first `nvec` smallest values in the `xm` array
    itmp = 0 ! indices
    xtmp = huge(1.0d0) ! values

    do i = 1, ubound(xm,1)
      do j = 1, nvec
        if (xtmp(j) > xm(i)) exit
      end do
      if (j <= nvec) then
        ! new small value found, insert it into temporary arrays
        xtmp(j+1:nvec) = xtmp(j:nvec-1)
        itmp(j+1:nvec) = itmp(j:nvec-1)
        xtmp(j) = xm(i)
        itmp(j) = i
      end if
    end do

  ! -- Get initial vectors: bvec(xvec_dim,nvec)
    bvec_mo = 0.0_dp
    do k = 1, nvec
      bvec_mo(itmp(k),k) = 1.0_dp
    end do

  end subroutine inivec

  subroutine iatogen(pv,av,nocca,noccb)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(in), contiguous, target :: pv(:)
    real(kind=dp), intent(out) :: av(:,:)
    integer, intent(in) :: nocca, noccb

    real(kind=dp), pointer :: ppv(:,:)
    integer :: nbf

    nbf = ubound(av, 1)
    ppv(1:nocca, noccb+1:nbf) => pv(:)

    av = 0.0_dp
    av(1:nocca, noccb+1:nbf) = ppv(1:nocca, noccb+1:nbf)

  end subroutine iatogen

  subroutine mntoia(pao,pmo,va,vb,nocca,noccb)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(in), dimension(:,:) :: pao
    real(kind=dp), intent(out), dimension(*) :: pmo
    real(kind=dp), intent(in), target, dimension(:,:) :: va, vb
    integer, intent(in) :: nocca, noccb

    integer :: nbf
    real(kind=dp), allocatable :: scr(:,:)
    real(kind=dp), pointer :: vap(:,:), vbp(:,:)

    nbf = ubound(pao, 1)

    allocate(scr(nocca,nbf))

    vap => va(:,1:nocca)
    vbp => vb(:,noccb+1:)

    call dgemm('t','n',nocca,nbf,nbf, &
               1.0_dp,vap,nbf,pao,nbf, &
               0.0_dp,scr,nocca)

    call dgemm('n','n',nocca,nbf-noccb,nbf, &
               1.0_dp,scr,nocca,vbp,nbf, &
               0.0_dp,pmo,nocca)

    deallocate(scr)
  end subroutine mntoia

  subroutine rparedms(b,ap_b,am_b,xm_p,xm_m,nvec,tamm_dancoff)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(in), dimension(:,:) :: b
    real(kind=dp), intent(in), dimension(:,:) :: am_b, ap_b
    real(kind=dp), intent(out), dimension(:,:) :: xm_m, xm_p
    integer, intent(in) :: nvec
    logical :: tamm_dancoff

    integer :: xvec_dim

    xvec_dim = ubound(b, 1)

    if (tamm_dancoff) then
      call dgemm('t','n',nvec,nvec,xvec_dim, &
                 1.0_dp,b,xvec_dim,ap_b,xvec_dim, &
                 0.0_dp,xm_p,nvec)
    else
      call dgemm('t','n',nvec,nvec,xvec_dim, &
                 1.0_dp,b,xvec_dim,ap_b,xvec_dim, &
                 0.0_dp,xm_p,nvec)
      call dgemm('t','n',nvec,nvec,xvec_dim, &
                 1.0_dp,b,xvec_dim,am_b,xvec_dim, &
                 0.0_dp,xm_m,nvec)
    end if
  end subroutine rparedms

!> @brief Diagonalize small reduced RPA matrix
  subroutine rpaeig(ee,vl,vr,apb,amb,scr,tamm_dancoff)

    use precision, only: dp
    use eigen, only: diag_symm_packed

    implicit none

    real(kind=dp), intent(out), dimension(:) :: ee
    real(kind=dp), intent(out), dimension(:,:) :: vl, vr
    real(kind=dp), intent(in), dimension(:,:) :: apb
    real(kind=dp), intent(inout), dimension(:,:) :: amb
    real(kind=dp), intent(inout), dimension(:) :: scr
    logical, intent(in) :: tamm_dancoff

    integer :: ierr, j, nvec

    nvec = ubound(vl, 2)

    if (tamm_dancoff) then

  !   Diagonailze A: VR
      if (nvec==1) then
        ee(1) = vr(1,1)
        vr(1,1) = 1.0_dp
      else
        call dtrttp('u', nvec, apb, nvec, scr, ierr)
        call diag_symm_packed(1,nvec,nvec,nvec,scr,ee,vr,ierr)
      end if
      return
    end if

!   sqrt(A-B)
    if (nvec==1) then
      amb(1,1) = sqrt(amb(1,1))
    else
      call dtrttp('u', nvec, amb, nvec, scr, ierr)
      call diag_symm_packed(1,nvec,nvec,nvec,scr,ee,vr,ierr)
      ee(1:nvec) = sign(sqrt(abs(ee(1:nvec))),ee(1:nvec))
      do j = 1, nvec
        vl(:,j) = vr(:,j)*ee(j)
      end do
      call dgemm('n','t',nvec,nvec,nvec, &
                 1.0_dp,vr,nvec,vl,nvec, &
                 0.0_dp,amb,nvec)
    end if

!   Form sqrt(A-B)*(A+B)*sqrt(A-B)
!   (A+B)*sqrt(A-B) : VL
    call dgemm('n','n',nvec,nvec,nvec, &
               1.0_dp,apb,nvec,amb,nvec, &
               0.0_dp,vl,nvec)
!   sqrt(A-B)*(A+B)*sqrt(A-B) : VR
    call dgemm('n','n',nvec,nvec,nvec, &
               1.0_dp,amb,nvec,vl,nvec, &
               0.0_dp,vr,nvec)

!   Diagonailze sqrt(A-B)*(A+B)*sqrt(A-B) : VR
    if (nvec==1) then
      ee(1) = vr(1,1)
      vl(1,1) = 1.0_dp
    else
      call dtrttp('u', nvec, vr, nvec, scr, ierr)
      call diag_symm_packed(1,nvec,nvec,nvec,scr,ee,vl,ierr)
    end if

!   Current vector VL is sqrt(1/(A-B)))|X+Y>.

!   VL into right eigenvector  VR = |X+Y>
    call dgemm('n','n',nvec,nvec,nvec, &
               1.0_dp,amb,nvec,vl,nvec, &
               0.0_dp,vr,nvec)

!   Left eigenvector   VL = = 1/E (A+B)|X+Y>
    call dgemm('n','n',nvec,nvec,nvec, &
               1.0_dp,apb,nvec,vr,nvec, &
               0.0_dp,vl,nvec)

    do j = 1, nvec
      vl(1:nvec,j) = vl(1:nvec,j)/sign(sqrt(abs(ee(j))),ee(j))
    end do

  end subroutine rpaeig

!> @brief Normalize `V1` and `V2` by biorthogonality condition
  subroutine rpavnorm(vr,vl,tamm_dancoff)

    use precision, only: dp

    implicit none

    real(kind=dp), intent(out), dimension(:,:) :: vr, vl
    logical, intent(in) :: tamm_dancoff

    real(kind=dp) :: scal, vrl, vrr
    integer :: ivec, nvec

    nvec = ubound(vl, 2)

    if (tamm_dancoff) then
      do ivec = 1, nvec
        vrr = dot_product(vr(:,ivec),vr(:,ivec))
        scal = sqrt(1.0_dp/vrr)
        vr(:,ivec) = vr(:,ivec)*scal
      end do
    else
      do ivec = 1, nvec
        vrl = dot_product(vr(:,ivec),vl(1:nvec,ivec))
        scal = sqrt(1.0D+00/vrl)
        vr(:,ivec)=vr(:,ivec)*scal
        vl(:,ivec)=vl(:,ivec)*scal
      end do
    end if
  end subroutine rpavnorm

!> @brief Remove negative eigenvalues
  subroutine rpaechk(ee,nvec,ndsr,imax,tamm_dancoff)

    use precision, only: dp

    implicit none

    real(kind=dp), intent(inout), dimension(:) :: ee
    integer, intent(in) :: nvec, ndsr
    integer, intent(inout) :: imax
    logical, intent(in) :: tamm_dancoff

  ! Number of negative eigenvalues : imax
    imax = count(ee(1:nvec)<0.0_dp)

    if (.not.tamm_dancoff) ee(1:ndsr) = sqrt(abs(ee(1:ndsr)))
  end subroutine rpaechk

!> @brief Print current excitation energies and errors
  subroutine rpaprint(ee,err,cnvtol,iter,imax,ndsr, do_neg)

    use precision, only: dp
    use physical_constants, only: UNITS_EV
    use io_constants, only: iw

    implicit none

    real(kind=dp), intent(in) :: ee(:)
    real(kind=dp), intent(in) :: cnvtol
    real(kind=dp), intent(in) :: err(:)
    integer, intent(in) :: iter
    integer, intent(in) :: imax
    integer, intent(in) :: ndsr
    logical, optional :: do_neg

    integer :: istat
    logical :: neg
    integer :: first

    neg = .false.
    if (present(do_neg)) neg = do_neg

    write(*,fmt='(/,4X,"Davidson iteration #",I4)') iter

    if (imax/=0.and..not.neg)  write(*,'(4X,"Number of negative eigenvalues =",I4)') imax
    first = imax + 1
    if (neg) first = 1

    do istat = 1, ndsr
      write(*,fmt='(4X,"State ",I4, &
             &3x,"E =",F12.6," eV",&
             &4x,"err. =",F10.6)') &
        istat, ee(istat)/UNITS_EV, err(istat)
    end do
    write(*,'(10X,"Max error =",1X,1P,E10.3,1X,"/",1P,E10.3)') maxval(err(first:ndsr)), cnvtol
    call flush(iw)
  end subroutine rpaprint

!> @brief Expand reduced vectors to real size space
  subroutine rpaexpndv(vr,vl,vro,vlo,br,bl,ndsr,tamm_dancoff)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(in) :: vr(:,:), vl(:,:)
    real(kind=dp), intent(out) :: vro(:,:), vlo(:,:)
    real(kind=dp), intent(in) :: br(:,:), bl(:,:)
    integer, intent(in) :: ndsr
    logical, intent(in) :: tamm_dancoff

    integer :: xvec_dim
    integer :: nvec

    nvec = ubound(vr, 2)
    xvec_dim = ubound(vro, 1)

    call dgemm('n','n',xvec_dim,ndsr,nvec, &
               1.0_dp,br,xvec_dim,vr,nvec, &
               0.0_dp,vro,xvec_dim)

    if (tamm_dancoff) then
      vlo(:,1:ndsr) = vro(:,1:ndsr)
    else
      call dgemm('n','n',xvec_dim,ndsr,nvec, &
                 1.0_dp,bl,xvec_dim,vl,nvec, &
                 0.0_dp,vlo,xvec_dim)
    end if

  end subroutine rpaexpndv

!> @brief Construct residual vectors and check convergence
  subroutine rparesvec(q,w_l,w_r,v_l,v_r,ee,abd, &
                       ndsr,errors,tol,imax,tamm_dancoff)

    use precision, only: dp

    implicit none

    real(kind=dp), intent(inout), dimension(:,:) :: w_l, w_r, v_l, v_r
    real(kind=dp), intent(out), &
      dimension(ubound(w_l,1),ndsr,*) :: q
    real(kind=dp), intent(in), dimension(:) :: ee
    real(kind=dp), intent(in), dimension(:) :: abd
    integer, intent(in) :: ndsr
    real(kind=dp), intent(inout) :: errors(:)
    real(kind=dp), intent(in) :: tol
    integer, intent(in) :: imax
    logical :: tamm_dancoff

    real(kind=dp) :: er_l, er_r
    integer :: ivec, xvec_dim

    xvec_dim = ubound(w_l,1)

!   Residual vector W_l and W_r
    do ivec = 1, ndsr
      w_r(1:xvec_dim,ivec) = w_r(1:xvec_dim,ivec) - ee(ivec)*v_r(1:xvec_dim,ivec)
    end do
    if (.not.tamm_dancoff) then
      do ivec = 1, ndsr
        w_l(1:xvec_dim,ivec) = w_l(1:xvec_dim,ivec) - ee(ivec)*v_l(1:xvec_dim,ivec)
      end do
    end if

!   Norms of W
    errors = 0.0_dp

    if (tamm_dancoff) then
      do ivec = imax+1, ndsr
        errors(ivec) = dot_product(w_r(:,ivec),w_r(:,ivec))
        if (errors(ivec)>1.0_dp) then
          write(*,'(4X,"Large error detected")')
          write(*,'(4X,"State#=",I4)') ivec
          write(*,'(4X,"Error right =",1X,1P,E10.3)') errors(ivec)
          errors(ivec) = 0.0_dp
        end if

!       Get new vectors Q
        if (errors(ivec)>tol) then
          q(:,ivec,1) = w_r(:,ivec)/(ee(ivec)-abd(:))
        else
          q(:,ivec,1) = 0.0_dp
        end if
      end do
    else
      do ivec = imax+1, ndsr
        er_l = dot_product(w_l(:,ivec),w_l(:,ivec))
        er_r = dot_product(w_r(:,ivec),w_r(:,ivec))
        errors(ivec) = max(er_l,er_r)
        if (errors(ivec)>1.0_dp) then
          write(*,'(4X,"Large error detected")')
          write(*,'(4X,"State#=",I4)') ivec
          write(*,'(4X,"Error left/right =",1X,1P,E10.3,"/",1X,1P,E10.3)') er_l, er_r
          errors(ivec) = 0.0_dp
        end if

!       Get new vectors Q
        if (errors(ivec)>tol) then
          q(:,ivec,1) = 1.0d0/(ee(ivec)-abd(:))
          q(:,ivec,2) = q(:,ivec,1)

          q(:,ivec,1) = q(:,ivec,1)*w_l(:,ivec)
          q(:,ivec,2) = q(:,ivec,2)*w_r(:,ivec)
        else
          q(:,ivec,1:2) = 0.0_dp
        end if
      end do
    end if

  end subroutine rparesvec

!> @brief Orthonormalize q(xvec_dim,ndsr*2) and append to bvec
  subroutine rpanewb(ndsr,bvec,q,novec,nvec,ick,tamm_dancoff)
    use precision, only: dp

    implicit none

    integer, intent(in) :: ndsr
    real(kind=dp), intent(out), dimension(:,:) :: bvec
    real(kind=dp), intent(inout), dimension(:,:) :: q
    integer, intent(inout) :: novec, nvec
    integer, intent(out) :: ick
    logical, intent(in) :: tamm_dancoff

    real(kind=dp) :: bq, fnorm
    integer :: istat, k, ms, ndsrt, mxvec
    real(kind=dp), parameter :: norm_threshold = 1.0D-09

    mxvec = ubound(bvec, 2)

!   Save nvec as novec
    novec = nvec
    ick = 0

!   Modified Gram-Schmidt orthonormalization
    ndsrt = ndsr*2
    if (tamm_dancoff) ndsrt = ndsr

    do k = 1, ndsrt
!     MGS: orthonormalize next vector w.r.t. all
!     previous vectors
      do istat = 1, nvec
        bq = dot_product(bvec(:,istat),q(:,k))
        q(:,k) = q(:,k) - bq*bvec(:,istat)
      end do

      fnorm = norm2(q(:,k))

!     Possible linear dependency, skip this vector
      if (fnorm<norm_threshold) cycle

      if (nvec==mxvec) then
!       Error termination, no space left for new vectors
        ick = 2
        return
      end if

!     Append new b vector
      nvec = nvec+1
      bvec(:,nvec) = q(:,k)/fnorm
    end do

!   Error termination, no vectors added
    ms = nvec-novec
    if (ms==0) ick = 3
  end subroutine rpanewb

!> @breif Add (E_a-E_i)*Z_ai to Pmo
  subroutine esum(e,pmo,z,nocc,ivec)

    use precision, only: dp

    implicit none

    real(kind=dp), intent(in) :: e(:)
    real(kind=dp), intent(inout) :: pmo(:,:)
    real(kind=dp), intent(in) :: z(:,:)
    integer, intent(in) :: nocc, ivec

    integer :: i, ij, j, nbf

    nbf = ubound(e, 1)

    do j = nocc+1, nbf
      do i = 1, nocc
        ij = (j-nocc-1)*nocc + i
        pmo(ij,ivec) = pmo(ij,ivec) + (e(j)-e(i))*z(ij,ivec)
      end do
    end do
  end subroutine esum

!> @brief Compute unrelaxed difference density matrix for R-TDDFT
  subroutine tdhf_unrelaxed_density(xmy, xpy, mo, t, nocc, tda)
    use precision, only: dp
    use mathlib, only: orthogonal_transform
    use mathlib, only: pack_matrix

    implicit none

    real(kind=dp), intent(in) :: xpy(*), xmy(*)
    real(kind=dp), intent(in) :: mo(:,:)
    real(kind=dp), intent(out) :: t(*)
    logical, intent(in) :: tda
    integer, intent(in) :: nocc

    integer :: nbf, nvir
    real(kind=dp), allocatable :: t_mo(:,:), t_ao(:,:), scr(:,:)

    nbf= ubound(mo,1)
    allocate(t_mo(nbf,nbf), t_ao(nbf,nbf), scr(nbf,nbf), source=0.0_dp)

    nvir = nbf-nocc

    if (tda) then
      ! vir->vib block
      call dgemm('t', 'n', nvir, nvir, nocc, &
                  1.0d0, xpy, nocc, &
                         xpy, nocc, &
                  0.0d0, t_mo(nocc+1:,nocc+1), nbf)

      ! occ->occ block
      call dgemm('n', 't', nocc, nocc, nvir,  &
                 -1.0d0, xpy, nocc, &
                         xpy, nocc,  &
                  0.0d0, t_mo,  nbf)
    else
      ! vir->vib block
      call dgemm('t', 'n', nvir, nvir, nocc, &
                 0.5d0, xpy, nocc, &
                        xpy, nocc,  &
                 0.0d0, t_mo(nocc+1:, nocc+1), nbf)
      call dgemm('t', 'n', nvir, nvir, nocc,  &
                 0.5d0, xmy, nocc, &
                        xmy, nocc,  &
                 1.0d0, t_mo(nocc+1:, nocc+1), nbf)

      ! occ->occ block
      call dgemm('n', 't', nocc, nocc, nvir,  &
                 -0.5d0, xpy, nocc, &
                         xpy, nocc,  &
                  0.0d0, t_mo, nbf)
      call dgemm('n', 't', nocc, nocc, nvir,  &
                 -0.5d0, xmy, nocc, &
                         xmy, nocc,  &
                  1.0d0, t_mo, nbf)
    endif

    ! MO->AO
    call orthogonal_transform('t', nbf, mo, t_mo, t_ao, scr)

    call pack_matrix(t_ao, t(1:(nbf+1)*nbf/2) )

  end subroutine

!###############################################################################

  subroutine int2_rpagrd_data_t_parallel_start(this, basis, nthreads)
    implicit none
    class(int2_rpagrd_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nthreads
    integer :: nbf, nsh

    this%nbf = basis%nbf
    nbf = this%nbf
    this%np = 0
    this%nm = 0
    this%nt = 0
    if (associated(this%xpy)) this%np = ubound(this%xpy, size(shape(this%xpy)))
    if (associated(this%xmy)) this%nm = ubound(this%xmy, size(shape(this%xmy)))
    if (associated(this%t))   this%nt = ubound(this%t, size(shape(this%t)))
    this%nthreads = nthreads
    nsh = basis%nshell

    if (this%cur_pass == 1) then
      if (allocated(this%hpp)) deallocate(this%hpp)
      if (allocated(this%hpt)) deallocate(this%hpt)
      if (allocated(this%hmm)) deallocate(this%hmm)
      if (allocated(this%dsh)) deallocate(this%dsh)

      if (this%np > 0) &
        allocate(this%hpp(nbf, nbf, this%nspin, this%np, nthreads), &
                 source=0.0d0)

      if (this%nt > 0) &
        allocate(this%hpt(nbf, nbf, this%nspin, this%nt, nthreads),  &
                 source=0.0d0)

      if (this%nm > 0) &
        allocate(this%hmm(nbf, nbf, this%nspin, this%nm, nthreads),  &
                 source=0.0d0)

      allocate(this%dsh(nsh,nsh), source=0.0d0)
    end if

    call this%init_screen(basis)

  end subroutine

!###############################################################################

  subroutine int2_rpagrd_data_t_parallel_stop(this)
    implicit none
    integer :: plast, tlast, mlast, nbf
    class(int2_rpagrd_data_t), intent(inout) :: this

    if (this%cur_pass /= this%num_passes) return
    plast = size(shape(this%hpp))
    tlast = size(shape(this%hpt))
    mlast = size(shape(this%hmm))
    nbf = this%nbf

    if (this%nthreads /= 1) then
      if (this%np > 0) this%hpp(:,:,:,:,lbound(this%hpp, plast )) = sum(this%hpp, dim=plast)
      if (this%nt > 0) this%hpt(:,:,:,:,lbound(this%hpt, tlast )) = sum(this%hpt, dim=tlast)
      if (this%nm > 0) this%hmm(:,:,:,:,lbound(this%hmm, mlast )) = sum(this%hmm, dim=mlast)
      this%nthreads = 1
    end if
    if (this%np > 0) call this%pe%allreduce(this%hpp(:,:,:,:,1),&
            size(this%hpp(:,:,:,:,1)))
    if (this%nt > 0) call this%pe%allreduce(this%hpt(:,:,:,:,1),&
            size(this%hpt(:,:,:,:,1)))
    if (this%nm > 0) call this%pe%allreduce(this%hmm(:,:,:,:,1),&
            size(this%hmm(:,:,:,:,1)))


    if (this%cur_pass == this%num_passes) then
      if (this%np > 0) call symmetrize_matrices(this%hpp, nbf, this%np*this%nspin)
      if (this%nt > 0) call symmetrize_matrices(this%hpt, nbf, this%nt*this%nspin)
    end if

  end subroutine

!###############################################################################

  subroutine int2_rpagrd_data_t_clean(this)
    implicit none
    class(int2_rpagrd_data_t), intent(inout) :: this
    if (allocated(this%hpp)) deallocate(this%hpp)
    if (allocated(this%hpt)) deallocate(this%hpt)
    if (allocated(this%hmm)) deallocate(this%hmm)
    if (allocated(this%dsh)) deallocate(this%dsh)
    nullify(this%xpy)
    nullify(this%xmy)
    nullify(this%t)
  end subroutine

!###############################################################################

  subroutine int2_rpagrd_data_t_init_screen(this, basis)
    implicit none
    class(int2_rpagrd_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis

!   Form shell density
    this%dsh = 0
    if (this%np > 0) call shlrpagrd(this%dsh, this%xpy, basis)
    if (this%np > 0) call shlrpagrd(this%dsh, this%xmy, basis)
    if (this%nt > 0) call shlrpagrd(this%dsh, this%t, basis)
    this%max_den = maxval(abs(this%dsh))

  end subroutine

!###############################################################################

  subroutine int2_rpagrd_data_t_update(this, buf)
    implicit none
    class(int2_rpagrd_data_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    integer :: mythread
    integer :: i

    mythread = buf%thread_id

    ! H+[X+Y]
    do i = 1, this%np
      call this%hplus(buf, this%hpp(:,:,:,i, mythread), this%xpy(:,:,:,i))
    end do

    ! H+[T]
    do i = 1, this%nt
      call this%hplus(buf, this%hpt(:,:,:,i,mythread), this%t(:,:,:,i))
    end do

    ! H-[X-Y]
    do i = 1, this%nm
      call this%hminus(buf, this%hmm(:,:,:,i,mythread), this%xmy(:,:,:,i))
    end do

    buf%ncur = 0

  end subroutine

!###############################################################################

!> @brief Compute H^+[V] over a set of integrals
  subroutine int2_rpagrd_data_t_update_hplus(this, buf, hp, v)
    implicit none
    class(int2_rpagrd_data_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    real(kind=dp), intent(inout) :: hp(:,:,:)
    real(kind=dp), intent(in) :: v(:,:,:)
    integer :: i, j, k, l, n
    real(kind=dp) :: xfact, cfact, &
                     val, xval, cval
    integer :: mythread

    mythread = buf%thread_id

    if (this%nspin == 1) then

      xfact = 2 * this%scale_exchange
      cfact = 8 * this%scale_coulomb

      do n = 1, buf%ncur
        i  = buf%ids(1,n)
        j  = buf%ids(2,n)
        k  = buf%ids(3,n)
        l  = buf%ids(4,n)
        val = buf%ints(n)

        xval = val*xfact
        cval = val*cfact
        ! Coulomb
        hp(i,j,1) = hp(i,j,1)+cval*v(l,k,1)
        hp(k,l,1) = hp(k,l,1)+cval*v(j,i,1)

!       Exhange
        hp(i,k,1) = hp(i,k,1)-xval*v(l,j,1)
        hp(i,l,1) = hp(i,l,1)-xval*v(k,j,1)
        hp(j,k,1) = hp(j,k,1)-xval*v(l,i,1)
        hp(j,l,1) = hp(j,l,1)-xval*v(k,i,1)
      end do

    else if (this%nspin == 2) then

      xfact = 1 * this%scale_exchange
      cfact = 2 * this%scale_coulomb

      do n = 1, buf%ncur
        i  = buf%ids(1,n)
        j  = buf%ids(2,n)
        k  = buf%ids(3,n)
        l  = buf%ids(4,n)
        val = buf%ints(n)

        xval = val*xfact
        cval = val*cfact
        ! Coulomb
        hp(i,j,1) = hp(i,j,1)+cval*(v(k,l,1)+v(l,k,1)+v(k,l,2)+v(l,k,2))
        hp(k,l,1) = hp(k,l,1)+cval*(v(i,j,1)+v(j,i,1)+v(i,j,2)+v(j,i,2))

        hp(i,j,2) = hp(i,j,2)+cval*(v(k,l,1)+v(l,k,1)+v(k,l,2)+v(l,k,2))
        hp(k,l,2) = hp(k,l,2)+cval*(v(i,j,1)+v(j,i,1)+v(i,j,2)+v(j,i,2))

!       Exhange
        hp(i,k,1) = hp(i,k,1)-xval*(v(j,l,1)+v(l,j,1))
        hp(i,l,1) = hp(i,l,1)-xval*(v(j,k,1)+v(k,j,1))
        hp(j,k,1) = hp(j,k,1)-xval*(v(i,l,1)+v(l,i,1))
        hp(j,l,1) = hp(j,l,1)-xval*(v(i,k,1)+v(k,i,1))
!
        hp(i,k,2) = hp(i,k,2)-xval*(v(j,l,2)+v(l,j,2))
        hp(i,l,2) = hp(i,l,2)-xval*(v(j,k,2)+v(k,j,2))
        hp(j,k,2) = hp(j,k,2)-xval*(v(i,l,2)+v(l,i,2))
        hp(j,l,2) = hp(j,l,2)-xval*(v(i,k,2)+v(k,i,2))
      end do

    end if

  end subroutine

!###############################################################################

!> @brief Compute H^-[V] over a set of integrals
  subroutine int2_rpagrd_data_t_update_hminus(this, buf, hm, v)
    implicit none
    class(int2_rpagrd_data_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    real(kind=dp), intent(inout) :: hm(:,:,:)
    real(kind=dp), intent(in) :: v(:,:,:)
    integer :: i, j, k, l, n
    real(kind=dp) :: xfact, val, xval
    integer :: mythread

    xfact = this%scale_exchange
    mythread = buf%thread_id

    do n = 1, buf%ncur
      i  = buf%ids(1,n)
      j  = buf%ids(2,n)
      k  = buf%ids(3,n)
      l  = buf%ids(4,n)
      val = buf%ints(n)

      xval = val*xfact

!     Exhange
      hm(i,k,1) = hm(i,k,1)+xval*(v(l,j,1)-v(j,l,1))
      hm(i,l,1) = hm(i,l,1)+xval*(v(k,j,1)-v(j,k,1))
      hm(j,k,1) = hm(j,k,1)+xval*(v(l,i,1)-v(i,l,1))
      hm(j,l,1) = hm(j,l,1)+xval*(v(k,i,1)-v(i,k,1))

!     Exhange
      hm(k,i,1) = hm(k,i,1)-xval*(v(l,j,1)-v(j,l,1))
      hm(l,i,1) = hm(l,i,1)-xval*(v(k,j,1)-v(j,k,1))
      hm(k,j,1) = hm(k,j,1)-xval*(v(l,i,1)-v(i,l,1))
      hm(l,j,1) = hm(l,j,1)-xval*(v(k,i,1)-v(i,k,1))
    end do
  end subroutine

!###############################################################################

  subroutine shlrpagrd(dsh,d,basis)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set

    implicit none

    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(out) :: dsh(:,:)
    real(kind=dp), intent(in), dimension(:,:,:,:) :: d

    integer :: ish, jsh, maxi, maxj, mini, minj
    real(kind=dp) :: mxv

!   RHF
    do ish = 1, basis%nshell
      mini = basis%ao_offset(ish)
      maxi = mini + basis%naos(ish)-1
      do jsh = 1, ish
        minj = basis%ao_offset(jsh)
        maxj = minj + basis%naos(jsh) - 1
        mxv = maxval(abs(d(minj:maxj,mini:maxi,:,:)))
        dsh(ish,jsh) = max(dsh(ish,jsh), mxv)
        dsh(jsh,ish) = dsh(ish,jsh)
      end do
    end do
  end subroutine shlrpagrd

!###############################################################################

  subroutine symmetrize_matrices(a, lda, nmtx)
    use mathlib, only: symmetrize_matrix
    implicit none
    real(kind=dp), intent(inout) :: a(lda,lda,*)
    integer, intent(in) :: lda, nmtx
    integer :: i

    do i = 1, nmtx
        call symmetrize_matrix(a(:,:,i), lda)
    end do
  end subroutine

!###############################################################################

end module tdhf_lib
