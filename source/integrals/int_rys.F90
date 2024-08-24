module int2e_rys
    use precision, only: dp
    use basis_tools, only: basis_set

    integer, private :: iii
    integer, parameter :: INT2_RYS_MXANG = 7
    integer, parameter :: INT2_RYS_MXCART = INT2_RYS_MXANG*(INT2_RYS_MXANG+1)/2
    integer, parameter :: ijkln(INT2_RYS_MXANG) = [ (iii*(iii+1)/2, iii = 1, INT2_RYS_MXANG) ]
    integer, parameter :: MAXCONTR = 120
   !< power of X in Cartesian Gaussian basis functions
    integer, parameter :: &
   ijklx(INT2_RYS_MXCART,INT2_RYS_MXANG) = reshape([ &
      [0,                                                               (0, iii = ijkln(1)+1, INT2_RYS_MXCART)], &
      [1, 0, 0,                                                         (0, iii = ijkln(2)+1, INT2_RYS_MXCART)], &
      [2, 0, 0, 1, 1, 0,                                                (0, iii = ijkln(3)+1, INT2_RYS_MXCART)], &
      [3, 0, 0, 2, 2, 1, 0, 1, 0, 1,                                    (0, iii = ijkln(4)+1, INT2_RYS_MXCART)], &
      [4, 0, 0, 3, 3, 1, 0, 1, 0, 2, 2, 0, 2, 1, 1,                     (0, iii = ijkln(5)+1, INT2_RYS_MXCART)], &
      [5, 0, 0, 4, 4, 1, 0, 1, 0, 3, 3, 2, 0, 2, 0, 3, 1, 1, 2, 2, 1,   (0, iii = ijkln(6)+1, INT2_RYS_MXCART)], &
      [6, 0, 0, 5, 5, 1, 0, 1, 0, 4, 4, 2, 0, 2, 0, 4, 1, 1, 3, 3, 0, 3, 3, 2, 1, 2, 1, 2]  &
      ], shape(ijklx))

   !< power of Y in Cartesian Gaussian basis functions
    integer, parameter :: &
   ijkly(INT2_RYS_MXCART,INT2_RYS_MXANG) = reshape([ &
      [0,                                                              (0, iii = ijkln(1)+1, INT2_RYS_MXCART)], &
      [0, 1, 0,                                                        (0, iii = ijkln(2)+1, INT2_RYS_MXCART)], &
      [0, 2, 0, 1, 0, 1,                                               (0, iii = ijkln(3)+1, INT2_RYS_MXCART)], &
      [0, 3, 0, 1, 0, 2, 2, 0, 1, 1,                                   (0, iii = ijkln(4)+1, INT2_RYS_MXCART)], &
      [0, 4, 0, 1, 0, 3, 3, 0, 1, 2, 0, 2, 1, 2, 1,                    (0, iii = ijkln(5)+1, INT2_RYS_MXCART)], &
      [0, 5, 0, 1, 0, 4, 4, 0, 1, 2, 0, 3, 3, 0, 2, 1, 3, 1, 2, 1, 2,  (0, iii = ijkln(6)+1, INT2_RYS_MXCART)], &
      [0, 6, 0, 1, 0, 5, 5, 0, 1, 2, 0, 4, 4, 0, 2, 1, 4, 1, 3, 0, 3, 2, 1, 3, 3, 1, 2, 2] &
      ], shape(ijkly))

   !< power of Z in Cartesian Gaussian basis functions
    integer, parameter :: &
   ijklz(INT2_RYS_MXCART,INT2_RYS_MXANG) = reshape([ &
      [0,                                                              (0, iii = ijkln(1)+1, INT2_RYS_MXCART)], &
      [0, 0, 1,                                                        (0, iii = ijkln(2)+1, INT2_RYS_MXCART)], &
      [0, 0, 2, 0, 1, 1,                                               (0, iii = ijkln(3)+1, INT2_RYS_MXCART)], &
      [0, 0, 3, 0, 1, 0, 1, 2, 2, 1,                                   (0, iii = ijkln(4)+1, INT2_RYS_MXCART)], &
      [0, 0, 4, 0, 1, 0, 1, 3, 3, 0, 2, 2, 1, 1, 2,                    (0, iii = ijkln(5)+1, INT2_RYS_MXCART)], &
      [0, 0, 5, 0, 1, 0, 1, 4, 4, 0, 2, 0, 2, 3, 3, 1, 1, 3, 1, 2, 2,  (0, iii = ijkln(6)+1, INT2_RYS_MXCART)], &
      [0, 0, 6, 0, 1, 0, 1, 5, 5, 0, 2, 0, 2, 4, 4, 1, 1, 4, 0, 3, 3, 1, 2, 1, 2, 3, 3, 2] &
      ], shape(ijklz))

    type int2_rys_data_t
      integer :: id(4)
      integer :: at(4)
      integer :: am(4)
      integer :: nbf(4)
      integer :: flips(4)
      integer :: nroots
      logical :: iandj, kandl, same
      real(kind=dp), allocatable :: gijkl(:)
      real(kind=dp), allocatable :: gnkl (:)
      real(kind=dp), allocatable :: gnm  (:)
      real(kind=dp), allocatable :: dij  (:,:)
      real(kind=dp), allocatable :: dkl  (:,:)
      real(kind=dp), allocatable :: b00  (:)
      real(kind=dp), allocatable :: b01  (:)
      real(kind=dp), allocatable :: b10  (:)
      real(kind=dp), allocatable :: c00  (:)
      real(kind=dp), allocatable :: d00  (:)
      real(kind=dp), allocatable :: abv  (:,:)
      real(kind=dp), allocatable :: PQ   (:,:)
      real(kind=dp), allocatable :: PB   (:,:)
      real(kind=dp), allocatable :: QD   (:,:)
      real(kind=dp), allocatable :: rw   (:,:)
      integer :: ijklxyz(4,INT2_RYS_MXCART,4)
      real(kind=dp) :: quartet_cutoff
      contains
        procedure :: init => gdat_init
        procedure :: clean => gdat_clean
        procedure :: set_ids => gdat_set_ids
    end type

    private
    public :: int2_rys_data_t
    public :: int2_rys_compute
    public :: rys_to_ghondo, rys_print_eri
    public :: INT2_RYS_MXANG

contains

  subroutine gdat_init(gdat, maxang, &
                        cutoffs, stat)

    use int2_pairs, only: int2_pair_storage, int2_cutoffs_t
    implicit none

    class(int2_rys_data_t), intent(inout) :: gdat
    type(int2_cutoffs_t) :: cutoffs
    integer, intent(in) :: maxang
    integer, intent(out) :: stat
    integer :: mxbra, mxcart, mxrys

    gdat%quartet_cutoff = cutoffs%pair_cutoff_squared

    mxrys = (4*maxang + 2)/2
    mxcart = maxang+1
    mxbra = 2*mxcart-1

    allocate(&
      gdat%gijkl(mxcart**4      *MAXCONTR*3), &
      gdat%gnkl (mxcart**2*mxbra*MAXCONTR*3), &
      gdat%gnm  (mxbra**2       *MAXCONTR*3), &
      gdat%dij  (3,mxcart**2    *MAXCONTR  ), &
      gdat%dkl  (3,mxbra        *MAXCONTR  ), &

      gdat%b00  (mxrys*MAXCONTR  ), &
      gdat%b01  (mxrys*MAXCONTR  ), &
      gdat%b10  (mxrys*MAXCONTR  ), &
      gdat%c00  (mxrys*MAXCONTR*3), &
      gdat%d00  (mxrys*MAXCONTR*3), &
      gdat%abv  (  6, MAXCONTR  ), &
      gdat%PQ   (  3, MAXCONTR  ), &
      gdat%PB   (  3, MAXCONTR  ), &
      gdat%QD   (  3, MAXCONTR  ), &
      gdat%rw   (mxrys*2, MAXCONTR  ), &

      stat=stat)
  end subroutine gdat_init

  subroutine gdat_clean(gdat)

    implicit none

    class(int2_rys_data_t), intent(inout) :: gdat
    if (allocated(gdat%gijkl)) deallocate(gdat%gijkl)
    if (allocated(gdat%gnkl )) deallocate(gdat%gnkl )
    if (allocated(gdat%gnm  )) deallocate(gdat%gnm  )
    if (allocated(gdat%dij  )) deallocate(gdat%dij  )
    if (allocated(gdat%dkl  )) deallocate(gdat%dkl  )
    if (allocated(gdat%b00  )) deallocate(gdat%b00  )
    if (allocated(gdat%b01  )) deallocate(gdat%b01  )
    if (allocated(gdat%b10  )) deallocate(gdat%b10  )
    if (allocated(gdat%c00  )) deallocate(gdat%c00  )
    if (allocated(gdat%d00  )) deallocate(gdat%d00  )
    if (allocated(gdat%abv  )) deallocate(gdat%abv  )
    if (allocated(gdat%PQ   )) deallocate(gdat%PQ   )
    if (allocated(gdat%PB   )) deallocate(gdat%PB   )
    if (allocated(gdat%QD   )) deallocate(gdat%QD   )
    if (allocated(gdat%rw   )) deallocate(gdat%rw   )
  end subroutine gdat_clean

  subroutine gdat_set_ids(gdat, basis, id)

    implicit none

    class(int2_rys_data_t), intent(inout) :: gdat
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)

    integer :: am(4)

    ! Permute shells so L_i < L_j, L_k < L_l, L_i < L_j
    gdat%flips = [1,2,3,4]
    am = basis%ktype(id)
    if (am(1) > am(2)) then
      gdat%flips(1:2) = [2,1]
      am(1:2) = am([2,1])
    end if
    if (am(3) > am(4)) then
      gdat%flips(3:4) = [4,3]
      am(3:4) = am([4,3])
    end if
    if (am(1)+am(2) > am(3)+am(4)) then
      gdat%flips = gdat%flips([3,4,1,2])
      am = am([3,4,1,2])
    end if

    gdat%id = id(gdat%flips)
    gdat%am = am

    gdat%at = basis%katom(gdat%id)

  end subroutine gdat_set_ids

  subroutine int2_rys_compute(ints, gdat, ppairs, zero_shq, mu2)

    use int2_pairs, only: int2_pair_storage
    implicit none

    type(int2_rys_data_t) :: gdat
    type(int2_pair_storage), intent(in) :: ppairs
    real(kind=dp), intent(inout) :: ints(*)
    real(kind=dp), intent(in), optional :: mu2
    logical :: zero_shq

    integer :: ijg, klg, maxgg, mmax, ng
    integer :: nmax
    real(kind=dp) :: aa, ab, aandb1, bb, da, db, test
    real(kind=dp) :: pfac, rho
    real(kind=dp) :: p(3), q(3)
    logical :: first
    integer :: id1, id2, ppid_p, ppid_q, npp_p, npp_q
    real(kind=dp) :: mu2_1

    ! Range-separation parameter for Erfc-attenuated integrals
    mu2_1 = 0
    if (present(mu2)) mu2_1 = 1.0d0/mu2

!   Prepare shell block
    call set_shells(gdat)


    id1 = maxval(gdat%id(1:2))
    id2 = minval(gdat%id(1:2))
    npp_p = ppairs%ppid(1,id1*(id1-1)/2+id2)
    ppid_p = ppairs%ppid(2,id1*(id1-1)/2+id2)

    id1 = maxval(gdat%id(3:4))
    id2 = minval(gdat%id(3:4))
    npp_q = ppairs%ppid(1,id1*(id1-1)/2+id2)
    ppid_q = ppairs%ppid(2,id1*(id1-1)/2+id2)

    zero_shq = npp_p*npp_q == 0
    if (zero_shq) return

    nmax = gdat%am(1)+gdat%am(2)-1
    mmax = gdat%am(3)+gdat%am(4)-1

    maxgg = MAXCONTR/gdat%nroots

!   Pair of k,l primitives

    first = .true.
    zero_shq = .true.
    ng = 0

    do klg = 1, npp_q
      db = ppairs%k(ppid_q-1+klg)*ppairs%ginv(ppid_q-1+klg)
      bb = ppairs%g(ppid_q-1+klg)
      q = ppairs%P(:,ppid_q-1+klg)

!     Pair of i,j primitives
      do ijg = 1, npp_p
        da = ppairs%k(ppid_p-1+ijg)*ppairs%ginv(ppid_p-1+ijg)
        aa = ppairs%g(ppid_p-1+ijg)
        p = ppairs%P(:,ppid_p-1+ijg)

        ! 2nd term is used for Erfc-attenuated integrals
        ! It is zero for regular integrals
        ab = (aa+bb) + aa*bb*mu2_1

        pfac = da*db
        test = pfac*pfac

        if (test<gdat%quartet_cutoff*ab) cycle

        ng = ng+1

        aandb1 = 1.0_dp/ab
        rho = aa*bb*aandb1

        gdat%abv(1,ng) = ppairs%ginv(ppid_p-1+ijg)
        gdat%abv(2,ng) = ppairs%ginv(ppid_q-1+klg)
        gdat%abv(3,ng) = rho
        gdat%abv(4,ng) = pfac*sqrt(aandb1)
        gdat%abv(5,ng) = aandb1
        gdat%abv(6,ng) = rho*sum((p-q)**2)

        gdat%PQ(:,ng) = p-q

        if (nmax>1) gdat%pb(:,ng) = ppairs%PB(:,ppid_p-1+ijg)
        if (mmax>1) gdat%qd(:,ng) = ppairs%PB(:,ppid_q-1+klg)

        gdat%dij(:,ng) = ppairs%PA(:,ppid_p-1+ijg)-ppairs%PB(:,ppid_p-1+ijg)
        gdat%dkl(:,ng) = ppairs%PA(:,ppid_q-1+klg)-ppairs%PB(:,ppid_q-1+klg)

        if (ng==maxgg) then
          if (ng/=0) then
            if (first) call clear_ints(gdat, ints)
            first = .false.
            call compute(gdat, ng, nmax, mmax, ints)
            ng = 0
            zero_shq = .false.
          end if
        end if

      end do
    end do
    if (ng/=0) then
      if (first) call clear_ints(gdat, ints)
      first = .false.
      call compute(gdat, ng, nmax, mmax, ints)
      ng = 0
      zero_shq = .false.
    end if


  end subroutine int2_rys_compute

  subroutine compute(gdat, ng, nmax, mmax, ints)

    type(int2_rys_data_t), intent(inout) :: gdat

    integer :: nmax, mmax, ng
    real(kind=dp), intent(inout) :: ints(*)

!   Compute roots and weights for quadrature
    call compute_rys_rw(gdat, gdat%rw, ng)

!   Compute coefficients for recursion formulae
    call compute_coefficients(gdat%b00, gdat%b01, gdat%b10, &
                gdat%c00, gdat%d00, gdat%gnm, &
                gdat%abv, gdat%pq, gdat%pb, gdat%qd, gdat%rw, nmax, mmax, ng, gdat%nroots)

!   Compute x, y, z integrals (2 centers, 2-d )
    call compute_xyz_p0q0(gdat%gnm,ng*gdat%nroots,nmax,mmax, &
                gdat%b00, gdat%b01, gdat%b10, gdat%c00, gdat%d00)

!   Compute x, y, z integrals (4 centers, 2-d)
    call compute_xyz_ijkl(gdat%gijkl, gdat%gnkl, gdat%gnm, &
                ng, gdat%nroots, nmax, mmax, &
                gdat%am(1), gdat%am(2), gdat%am(3), gdat%am(4), &
                gdat%dij,gdat%dkl)

!   compute integrals
    call compute_ints(gdat, ng*gdat%nroots, gdat%ijklxyz, gdat%gijkl, ints)
  end subroutine

  subroutine set_shells(gdat)

    implicit none

    type(int2_rys_data_t) :: gdat
    integer :: ish, jsh, ksh, lsh

    ish = gdat%id(1)
    jsh = gdat%id(2)
    ksh = gdat%id(3)
    lsh = gdat%id(4)

    gdat%iandj = ish == jsh
    gdat%kandl = ksh == lsh
    gdat%same = ish == ksh.and.jsh == lsh

    gdat%nbf = ijkln(gdat%am)

!   Set number of quadrature points
    gdat%nroots = (sum(gdat%am)-2 )/2

!   Prepare indices for pairs of (i,j) functions
    call prepare_xyz_ids(gdat)

  end subroutine set_shells

  subroutine prepare_xyz_ids(gdat)
    implicit none
    class(int2_rys_data_t), intent(inout) :: gdat
    integer :: i, nj, nk, nl, njkl, nkl

    nj = gdat%am(2)
    nk = gdat%am(3)
    nl = gdat%am(4)
    njkl = nl*nk*nj
    do i = 1, gdat%nbf(1)
      gdat%ijklxyz(1,i,1) = ijklx(i,gdat%am(1))*njkl
      gdat%ijklxyz(2,i,1) = ijkly(i,gdat%am(1))*njkl
      gdat%ijklxyz(3,i,1) = ijklz(i,gdat%am(1))*njkl
    end do

    nkl = nl*nk
    do i = 1, gdat%nbf(2)
      gdat%ijklxyz(1,i,2) = ijklx(i,gdat%am(2))*nkl
      gdat%ijklxyz(2,i,2) = ijkly(i,gdat%am(2))*nkl
      gdat%ijklxyz(3,i,2) = ijklz(i,gdat%am(2))*nkl
    end do

!   Prepare indices for pairs of (k,l) functions
    do i = 1, gdat%nbf(3)
      gdat%ijklxyz(1,i,3) = ijklx(i,gdat%am(3))*nl
      gdat%ijklxyz(2,i,3) = ijkly(i,gdat%am(3))*nl
      gdat%ijklxyz(3,i,3) = ijklz(i,gdat%am(3))*nl
    end do

    do i = 1, gdat%nbf(4)
      gdat%ijklxyz(1,i,4) = ijklx(i,gdat%am(4))+1
      gdat%ijklxyz(2,i,4) = ijkly(i,gdat%am(4))+1
      gdat%ijklxyz(3,i,4) = ijklz(i,gdat%am(4))+1
    end do
  end subroutine

  subroutine compute_rys_rw(gdat, rwv, numg)
    use rys, only: rys_root_t
    implicit none


    class(int2_rys_data_t), intent(inout) :: gdat
    real(kind=dp), intent(out) :: rwv(2,numg,*)
    integer, intent(in) :: numg

    type(rys_root_t) :: root
    integer :: ng

    root%nroots = gdat%nroots
    do ng = 1, numg
      root%x = gdat%abv(6,ng)
      call root%evaluate
      rwv(1,ng,1:gdat%nroots) = root%u(1:gdat%nroots)
      rwv(2,ng,1:gdat%nroots) = root%w(1:gdat%nroots)
    end do

  end subroutine compute_rys_rw

  subroutine compute_coefficients(b00,b01,b10,c00,d00,gnm, &
                  abv,pq,pb,qd,rwv,nmax,mmax,numg,nroots)

    implicit none

    real(kind=dp) :: b00(numg,*),b01(numg,*),b10(numg,*) !numg*nroots
    real(kind=dp) :: c00(numg,nroots,*) !numg*nroots*3
    real(kind=dp) :: d00(numg,nroots,*) !numg*nroots*3
    real(kind=dp) :: gnm(numg,nroots,*) !numg*nroots*3

    real(kind=dp) :: abv(6,*), pq(3,*), pb(3,*), qd(3,*)
    real(kind=dp) :: rwv(2,numg,*) !rwv(2,numg,nroots)
    integer :: mmax, nmax
    integer :: numg, nroots

    integer :: nr, ng
    real(kind=dp) :: a1, b1, ab1
    real(kind=dp) :: pfac, rho, t2, t2ar, t2br, uu, ww

    do nr = 1, nroots
      do ng = 1, numg
        a1  = abv(1,ng)
        b1  = abv(2,ng)
        rho = abv(3,ng)
        pfac = abv(4,ng)
        ab1 = abv(5,ng)
        uu  = rwv(1,ng,nr)
        ww  = rwv(2,ng,nr)

!       G(0,0)
        gnm(ng,nr,1) = ww*pfac
        gnm(ng,nr,2) = 1.0_dp
        gnm(ng,nr,3) = 1.0_dp

        t2    = uu/(uu+1)
        t2ar  = t2*rho*a1
        t2br  = t2*rho*b1
        b00(ng,nr) = 0.5_dp*ab1*t2
        b01(ng,nr) = 0.5_dp*b1*(1.0_dp-t2br)
        b10(ng,nr) = 0.5_dp*a1*(1.0_dp-t2ar)

        if (mmax>1) then
          d00(ng,nr,1) = qd(1,ng) + t2br*pq(1,ng)
          d00(ng,nr,2) = qd(2,ng) + t2br*pq(2,ng)
          d00(ng,nr,3) = qd(3,ng) + t2br*pq(3,ng)
        end if

        if (nmax>1) then
          c00(ng,nr,1) = pb(1,ng) - t2ar*pq(1,ng)
          c00(ng,nr,2) = pb(2,ng) - t2ar*pq(2,ng)
          c00(ng,nr,3) = pb(3,ng) - t2ar*pq(3,ng)
        end if

      end do
    end do

  end subroutine compute_coefficients

  subroutine compute_xyz_p0q0(gnm,ng,nmax,mmax,b00,b01,b10,c00,d00)

    implicit none

    real(kind=dp) :: gnm(ng,3,nmax,*)
    real(kind=dp) :: c00(ng,*),d00(ng,*)
    real(kind=dp) :: b00(*),b01(*),b10(*)
    integer :: ng, nmax, mmax

    integer :: m, n, xyz

    if ( max(nmax,mmax) == 1 ) return

    if (nmax>1) then
!     g(1,0) = c00 * g(0,0)
      gnm(:ng,1:3,2,1) = c00(:ng,1:3)*gnm(:ng,1:3,1,1)
    end if

    if (mmax>1) then
!     g(0,1) = d00 * g(0,0)
      gnm(:ng,1:3,1,2) = d00(:ng,1:3)*gnm(:ng,1:3,1,1)

      if (nmax>1) then
!       g(1,1) = b00 * g(0,0) + d00 * g(1,0)
        do xyz = 1, 3
          gnm(:ng,xyz,2,2) = b00(:ng)    *gnm(:ng,xyz,1,1) &
                           + d00(:ng,xyz)*gnm(:ng,xyz,2,1)
        end do
      end if
    end if

    if (nmax>2) then
!     g(n+1,0) = n * b10 * g(n-1,0) + c00 * g(n,0)
      do n = 2, nmax-1
        do xyz = 1, 3
          gnm(:ng,xyz,n+1,1) = (n-1)*b10(:ng)    *gnm(:ng,xyz,n-1,1) &
                             +       c00(:ng,xyz)*gnm(:ng,xyz,n  ,1)
        end do
      end do

      if (mmax>1) then

!       g(n,1) = n * b00 * g(n-1,0) + d00 * g(n,0)
        do n = 2, nmax-1
          do xyz = 1, 3
            gnm(:ng,xyz,n+1,2) = n*b00(:ng)    *gnm(:ng,xyz,n  ,1) &
                               +   d00(:ng,xyz)*gnm(:ng,xyz,n+1,1)
          end do
        end do
      end if

    end if

    if (mmax<3) return

!   g(0,m+1) = m * b01 * g(0,m-1) + d00 * g(o,m)
    do m = 2, mmax-1
      do xyz = 1, 3
        gnm(:ng,xyz,1,m+1) = (m-1)*b01(:ng)    *gnm(:ng,xyz,1,m-1) &
                           +       d00(:ng,xyz)*gnm(:ng,xyz,1,m  )
      end do
    end do

    if (nmax<2) return

!   g(1,m) = m * b00 * g(0,m-1) + c00 * g(0,m)
    do m = 2, mmax-1
      do xyz = 1, 3
        gnm(:ng,xyz,2,m+1) = m*b00(:ng)    *gnm(:ng,xyz,1,m  ) &
                       +       c00(:ng,xyz)*gnm(:ng,xyz,1,m+1)
      end do
    end do

    if (nmax<3) return
!   g(n+1,m) = n * b10 * g(n-1,m  )
!            +     c00 * g(n  ,m  )
!            + m * b00 * g(n  ,m-1)

    do m = 2, mmax-1
      do n = 2, nmax-1
        do xyz = 1, 3
          gnm(:,xyz,n+1,m+1) = (n-1)*b10(:ng)    *gnm(:ng,xyz,n-1,m+1) &
                         +           c00(:ng,xyz)*gnm(:ng,xyz,n  ,m+1) &
                         +         m*b00(:ng)    *gnm(:ng,xyz,n  ,m  )
        end do
      end do
    end do

  end subroutine compute_xyz_p0q0

  subroutine compute_xyz_ijkl(ijkl, gnkl, gnm, ng, nr,  &
                  nmax, mmax, nimax, njmax, nkmax, nlmax, dij, dkl)

    implicit none

    real(kind=dp) :: ijkl(ng,nr,3,nlmax,nkmax,njmax,*)
    real(kind=dp) :: gnkl(ng,nr,3,nlmax,nkmax,*)
    real(kind=dp) ::   gnm(ng,nr,3,nmax,*) ! gnm(ng,nmax,mmax)
    real(kind=dp) ::   dij(3,*) ! dij(ng)
    real(kind=dp) ::   dkl(3,*) ! dkl(ng)
    integer :: ng,nr,nmax,mmax,nimax,njmax,nkmax,nlmax

    integer :: ni, nk, nl, ig, m1, n1, xyz

!   g(n,k,l)
    do nk=1, nkmax
      do nl=1,nlmax
        gnkl(:,:,:,nl,nk,:nmax) = gnm(:,:,:,:,nl)
      end do
      if(nk == nkmax) exit
      m1 = mmax-nk
      do xyz = 1, 3
        do ig = 1, ng
          gnm(ig,:,xyz,:,1:m1) = dkl(xyz,ig)*gnm(ig,:,xyz,:,1:m1) &
                               +             gnm(ig,:,xyz,:,2:m1+1)
        end do
      end do
    end do

!   g(i,j,k,l)
    do ni = 1, nimax
      ijkl(:,:,:,:,:,1:njmax,ni) = gnkl(:,:,:,:,:,1:njmax)
      if (ni == nimax) exit
      n1 = nmax-ni
      do xyz = 1, 3
        do ig = 1, ng
          gnkl(ig,:,xyz,:,:,1:n1) = dij(xyz,ig)*gnkl(ig,:,xyz,:,:,1:n1) &
                                  +             gnkl(ig,:,xyz,:,:,2:n1+1)
        end do
      end do
    end do

  end subroutine compute_xyz_ijkl

  subroutine clear_ints(gdat,ints)

    implicit none

    type(int2_rys_data_t) :: gdat
    real(kind=dp) :: ints(*)
    ints(1:product(gdat%nbf)) = 0

  end subroutine clear_ints

  subroutine compute_ints(gdat,ngnr,ijklxyz,g0,ints)

    implicit none

    type(int2_rys_data_t) :: gdat
    integer :: ngnr
    integer :: ijklxyz(:,:,:)
    real(kind=dp) ::  g0(ngnr,3,*)
    real(kind=dp), target :: ints(*)

    integer :: i, j, k, l
    integer :: nx, ny, nz
    integer :: am(4)
    real(kind=dp), pointer :: p(:,:,:,:)

    p(1:gdat%nbf(4),1:gdat%nbf(3),1:gdat%nbf(2),1:gdat%nbf(1)) => ints(1:product(gdat%nbf))
    am = gdat%am

    do i = 1, gdat%nbf(1)
      do j = 1, gdat%nbf(2)
        do k = 1, gdat%nbf(3)
          do l = 1, gdat%nbf(4)
            nx = ijklxyz(1,i,1)+ijklxyz(1,j,2)+ijklxyz(1,k,3)+ijklxyz(1,l,4)
            ny = ijklxyz(2,i,1)+ijklxyz(2,j,2)+ijklxyz(2,k,3)+ijklxyz(2,l,4)
            nz = ijklxyz(3,i,1)+ijklxyz(3,j,2)+ijklxyz(3,k,3)+ijklxyz(3,l,4)

            associate ( x  => g0(:,1,nx) &
                      , y  => g0(:,2,ny) &
                      , z  => g0(:,3,nz) &
                      )
                p(l,k,j,i) = p(l,k,j,i) + sum(x*y*z)
            end associate

          end do
        end do
      end do
    end do

  end subroutine compute_ints

  subroutine rys_to_ghondo(gdat, ints, ghondo)
    use constants, only: shells_pnrm2
    implicit none
    type(int2_rys_data_t), intent(in) :: gdat
    real(kind=dp), contiguous, intent(inout) :: ints(:,:,:,:), ghondo(:,:,:,:)

    integer :: am(4), n0(4), am0(4), ids(4), inv(4), na, nb, nc, nd
    real(kind=dp), pointer:: pnorma(:), pnormb(:), pnormc(:), pnormd(:)
    real(kind=dp) :: scaleab, scalecd

    inv(gdat%flips) = [1,2,3,4]
    am = gdat%am
    am0 = am(inv)
    n0 = am0*(am0+1)/2

    pnorma => shells_pnrm2(:,am0(1)-1)
    pnormb => shells_pnrm2(:,am0(2)-1)
    pnormc => shells_pnrm2(:,am0(3)-1)
    pnormd => shells_pnrm2(:,am0(4)-1)

    do na = 1, n0(1)
    do nb = 1, n0(2)
    scaleab =  pnorma(na) * pnormb(nb)
    do nc = 1, n0(3)
    do nd = 1, n0(4)
      ids = [na, nb, nc, nd]
      ids = ids(gdat%flips)
      scalecd =  pnormc(nc) * pnormd(nd)
      ghondo(nd,nc,nb,na) = &
        ints(ids(4),ids(3),ids(2),ids(1)) * scaleab * scalecd
    end do
    end do
    end do
    end do

  end subroutine

  subroutine rys_print_eri(gdat, ints)
    use constants, only: shells_pnrm2
    implicit none
    type(int2_rys_data_t), intent(in) :: gdat
    real(kind=dp), intent(in) :: ints(:,:,:,:)
    integer :: n(4), n0(4), na, nb, nc, nd
    integer :: am(4), am0(4)
    integer :: ids(4), inv(4)
    real(kind=dp), pointer:: pnorma(:), pnormb(:), pnormc(:), pnormd(:)

    am = gdat%am

    n = am*(am+1)/2
    inv(gdat%flips) = [1,2,3,4]

    am0 = am(inv)
    n0 = am0*(am0+1)/2

    pnorma => shells_pnrm2(:,am0(1)-1)
    pnormb => shells_pnrm2(:,am0(2)-1)
    pnormc => shells_pnrm2(:,am0(3)-1)
    pnormd => shells_pnrm2(:,am0(4)-1)

    do na = 1, n0(1)
    do nb = 1, n0(2)
    do nc = 1, n0(3)
    do nd = 1, n0(4)
      ids = [na, nb, nc, nd]
      ids = ids(gdat%flips)
      write (*, "(a6, 2i3, a, 2i3, a, es30.15)") &
        "elem (", na, nb, " |", nc, nd,  ") = ", &
        ints(ids(4),ids(3),ids(2),ids(1)) * &
        pnorma(na)*pnormb(nb)*pnormc(nc)*pnormd(nd)
    end do
    end do
    end do
    end do

  end subroutine
end module int2e_rys
