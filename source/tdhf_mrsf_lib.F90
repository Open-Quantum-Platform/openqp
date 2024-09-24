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
      tmp, tmp1, tmp2
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

    allocate(tmp(nbf,nbf), &
             tmp1(nbf,4), &
             tmp2(nbf,4), &
             source=0.0_dp, stat=ok)
    if( ok/=0 ) call show_message('Cannot allocate memory',with_abort)

    do j = nocca+1, nbf
      tmp1(:,1) = tmp1(:,1)+vb(:,j)*bvec(lr2,j)
      tmp1(:,2) = tmp1(:,2)+vb(:,j)*bvec(lr1,j)
      tmp1(:,3) = tmp1(:,3)+vb(:,j)*bvec(lr2,j)
      tmp1(:,4) = tmp1(:,4)+vb(:,j)*bvec(lr1,j)
    end do

    do i = 1, noccb
      tmp2(:,1) = tmp2(:,1)+va(:,i)*bvec(i,lr1)
      tmp2(:,2) = tmp2(:,2)+va(:,i)*bvec(i,lr2)
      tmp2(:,3) = tmp2(:,3)+va(:,i)*bvec(i,lr1)
      tmp2(:,4) = tmp2(:,4)+va(:,i)*bvec(i,lr2)
    end do

    do m = 1, nbf
      bo2v(:,m) = bo2v(:,m)+va(:,lr2)*tmp1(m,1)
      bo1v(:,m) = bo1v(:,m)+va(:,lr1)*tmp1(m,2)
      bco1(:,m) = bco1(:,m)+vb(m,lr1)*tmp2(:,1)
      bco2(:,m) = bco2(:,m)+vb(m,lr2)*tmp2(:,2)
    end do

    do m = 1, nbf
      o21v(m,:) = o21v(m,:)+va(:,lr1)*tmp1(m,3) &
                           -va(:,lr2)*tmp1(m,4)
      co12(m,:) = co12(m,:)+vb(m,lr2)*tmp2(:,3) &
                           -vb(m,lr1)*tmp2(:,4)
    end do

    ball = ball + bo2v + bo1v + bco1 + bco2

    call dgemm('n', 't', nbf, noccb, nbf-nocca, &
               1.0_dp, vb(:,nocca+1), nbf, &
                       bvec(:,nocca+1), nbf, &
               0.0_dp, tmp, nbf)

    call dgemm('n', 't', nbf, nbf, noccb, &
               1.0_dp, va, nbf, &
                       tmp, nbf, &
               1.0_dp, ball, nbf)

    if (mrst==1) then
      do m = 1, nbf
        ball(:,m) = ball(:,m) &
          +va(:,lr2)*bvec(lr2,lr1)*vb(m,lr1) &
          +va(:,lr1)*bvec(lr1,lr2)*vb(m,lr2) &
          +(va(:,lr1)*vb(m,lr1)-va(:,lr2)*vb(m,lr2)) &
             *bvec(lr1,lr1)*isqrt2
      end do
    else if (mrst==3) then
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
! subroutine umrsfcbc(infos, va, vb, bvec, mrsf_density)

!   use messages, only: show_message, with_abort
!   use types, only: information
!   use io_constants, only: iw
!   use precision, only: dp

!   implicit none

!   type(information), intent(in) :: infos
!   real(kind=dp), intent(in), dimension(:,:) :: &
!     va, vb, bvec
!   real(kind=dp), intent(inout), target, dimension(:,:,:) :: &
!     mrsf_density

!   real(kind=dp), allocatable, dimension(:,:) :: &
!     tmp, tmp1, tmp2
!   real(kind=dp), pointer, dimension(:,:) :: &
!     bo2va, bo2vb, bo1va, bo1vb, bco1a, bco1b, &
!     bco2a, bco2b, ball, co12, o21v
!   integer :: nocca, noccb, mrst, i, j, m, nbf, lr1, lr2, ok
!   logical :: debug_mode
!   real(kind=dp), parameter :: isqrt2 = 1.0_dp/sqrt(2.0_dp)

!   ball => mrsf_density(11,:,:)
!   bo2va => mrsf_density(1,:,:)
!   bo2vb => mrsf_density(2,:,:)
!   bo1va => mrsf_density(3,:,:)
!   bo1vb => mrsf_density(4,:,:)
!   bco1a => mrsf_density(5,:,:)
!   bco1b => mrsf_density(6,:,:)
!   bco2a => mrsf_density(7,:,:)
!   bco2b => mrsf_density(8,:,:)
!   o21v => mrsf_density(9,:,:)
!   co12 => mrsf_density(10,:,:)

!   nbf = infos%basis%nbf
!   nocca = infos%mol_prop%nelec_A
!   noccb = infos%mol_prop%nelec_B
!   mrst = infos%tddft%mult
!   debug_mode = infos%tddft%debug_mode

!   lr1 = nocca-1
!   lr2 = nocca
!   allocate(tmp(nbf,nbf), &
!            tmp1(nbf,4), &
!            tmp2(nbf,4), &
!            source=0.0_dp, stat=ok)
!   if( ok/=0 ) call show_message('Cannot allocate memory',with_abort)

!   do j = nocca+1, nbf
!     tmp1(:,1) = tmp1(:,1)+va(:,j)*bvec(nocca,j)
!     tmp1(:,2) = tmp1(:,2)+vb(:,j)*bvec(nocca,j)
!     tmp1(:,3) = tmp1(:,3)+va(:,j)*bvec(nocca-1,j)
!     tmp1(:,4) = tmp1(:,4)+vb(:,j)*bvec(nocca-1,j)
!   end do

!   do i = 1, nocca-2
!     tmp2(:,1) = tmp2(:,1)+va(:,i)*bvec(i,nocca-1)
!     tmp2(:,2) = tmp2(:,2)+vb(:,i)*bvec(i,nocca-1)
!     tmp2(:,3) = tmp2(:,3)+va(:,i)*bvec(i,nocca)
!     tmp2(:,4) = tmp2(:,4)+vb(:,i)*bvec(i,nocca)
!   end do

!   do m = 1, nbf
!     ball(:,m) = ball(:,m)+tmp1(m,2)*va(:,nocca) &
!                          +tmp1(m,4)*va(:,nocca-1) &
!                          +tmp2(:,1)*vb(m,nocca-1) &
!                          +tmp2(:,3)*vb(m,nocca)

!     bo2va(:,m) = bo2va(:,m)+tmp1(m,1)*va(:,nocca)
!     bo2vb(:,m) = bo2vb(:,m)+tmp1(m,2)*vb(:,nocca)
!     bo1va(:,m) = bo1va(:,m)+tmp1(m,3)*va(:,nocca-1)
!     bo1vb(:,m) = bo1vb(:,m)+tmp1(m,4)*vb(:,nocca-1)
!     o21v(:,m) = o21v(:,m)+tmp1(m,2)*va(:,nocca-1) &
!                          -tmp1(m,4)*va(:,nocca)

!     bco1a(:,m) = bco1a(:,m)+tmp2(:,1)*va(m,nocca-1)
!     bco1b(:,m) = bco1b(:,m)+tmp2(:,2)*vb(m,nocca-1)
!     bco2a(:,m) = bco2a(:,m)+tmp2(:,3)*va(m,nocca)
!     bco2b(:,m) = bco2b(:,m)+tmp2(:,4)*vb(m,nocca)
!     co12(:,m) = co12(:,m)+tmp2(:,1)*vb(m,nocca) &
!                          -tmp2(:,3)*vb(m,nocca-1)
!   end do

!   call dgemm('n','t', nbf, nocca-2, nbf-nocca, &
!              1.0_dp, vb(:,nocca+1), nbf, &
!                      bvec(:,nocca+1), nbf, &
!              0.0_dp, tmp, nbf)

!   call dgemm('n','t', nbf, nbf, nocca-2, &
!              1.0_dp, va, nbf, &
!                      tmp, nbf, &
!              1.0_dp, ball, nbf)

!   if (mrst==1) then
!     do m = 1, nbf
!       ball(:,m) = ball(:,m) &
!       +va(:,noca)*bvec(noca,noca-1)*vb(m,noca-1) &
!       +va(:,noca-1)*bvec(noca-1,noca)*vb(m,noca) &
!       +(va(:,noca-1)*vb(m,noca-1)-va(:,noca)*vb(m,noca)) &
!       *bvec(noca-1,noca-1)*isqrt2
!     end do
!   elseif (mrsft) then
!     do m = 1, nbf
!       ball(:,m) = ball(:,m)
!       +(va(:,noca-1)*vb(m,noca-1)+va(:,noca)*vb(m,noca))
!       *bvec(noca-1,noca-1)*sqrt2
!     end do
!   endif

! return
! end subroutine umrsfcbc

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

! 3
    call dgemv('n',nbf,nbf,one,ado1v,nbf,vb(:,lr2),1,zero,tmp,1)
    call dgemv('n',nbf,nbf,one,aco12,nbf,vb(:,lr1),1,one,tmp,1)
    do i = 1, noca-2
      wrk(i,lr2) = wrk(i,lr2) + dot_product(va(:,i),tmp(:))
    end do
! 4
    call dgemv('n',nbf,nbf,one,aco12,nbf,vb(:,lr2),1,zero,tmp,1)
    call dgemv('n',nbf,nbf,one,ado2v,nbf,vb(:,lr1),1,-one,tmp,1)
    do i = 1, noca-2
      wrk(i,lr1) = wrk(i,lr1) + dot_product(va(:,i),tmp(:))
    end do
! 5
    call dgemv('t',nbf,nbf,one,adco2,nbf,va(:,lr1),1,zero,tmp,1)
    call dgemv('t',nbf,nbf,one,ao21v,nbf,va(:,lr2),1, one,tmp,1)
    do j = noca+1, nbf
      wrk(lr1,j) = wrk(lr1,j) + dot_product(vb(:,j),tmp(:))
    end do
! 6
    call dgemv('t',nbf,nbf,one,ao21v,nbf,va(:,lr1),1,zero,tmp,1)
    call dgemv('t',nbf,nbf,one,adco1,nbf,va(:,lr2),1,-one,tmp,1)
    do j = noca+1, nbf
       wrk(lr2,j) = wrk(lr2,j) + dot_product(vb(:,j),tmp(:))
    end do

    if (mrst==1) then
      wrk(lr1,lr1) = (scr(lr1,lr1)-scr(lr2,lr2))*sqrt2
      wrk(lr2,lr2) = 0.0_dp
    else if (mrst==3) then
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
