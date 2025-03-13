module grd2_rys
    use precision, only: dp
    use basis_tools, only: basis_set
    use constants, only: bas_mxang, bas_mxcart, num_cart_bf, cart_x, cart_y, cart_z

    integer, private :: iii
    integer, parameter :: MAXCONTR = 120

    logical, parameter :: skips(4,16) = reshape([&
            .true. , .true. , .true. , .true. , &  !  1      ( 1, 1 | 1, 1 )
            .true. , .true. , .true. , .false., &  !  2      ( 1, 1 | 1, 2 )
            .true. , .true. , .false., .true. , &  !  3      ( 1, 1 | 2, 1 )
            .true. , .true. , .false., .false., &  !  4      ( 1, 1 | 2, 2 )
            .true. , .true. , .false., .false., &  !  5      ( 1, 1 | 2, 3 )
            .true. , .false., .true. , .true. , &  !  6      ( 1, 2 | 1, 1 )
            .true. , .false., .true. , .false., &  !  7      ( 1, 2 | 1, 2 )
            .true. , .false., .true. , .false., &  !  8      ( 1, 2 | 1, 3 )
            .true. , .false., .false., .true. , &  !  9      ( 1, 2 | 2, 1 )
            .true. , .false., .false., .true. , &  !  10     ( 1, 2 | 3, 1 )
            .false., .true. , .true. , .true. , &  !  11     ( 1, 2 | 2, 2 )
            .false., .true. , .true. , .false., &  !  12     ( 1, 2 | 2, 3 )
            .false., .true. , .false., .true. , &  !  13     ( 1, 2 | 3, 2 )
            .false., .false., .true. , .true. , &  !  14     ( 1, 2 | 3, 3 )
            .false., .false., .false., .true. , &  !  15     ( 1, 2 | 3, 4 )
            .false., .false., .false., .false.  &  !  16
    ], shape(skips))

    type grd2_int_data_t
      logical :: skip(4)
      integer :: id(4)
      integer :: at(4)
      integer :: am(4)
      integer :: nbf(4)
      integer :: der(4)
      integer :: nder
      integer :: nroots
      integer :: invtyp
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
      real(kind=dp), allocatable :: f00  (:)
      real(kind=dp), allocatable :: abv  (:,:)
      real(kind=dp), allocatable :: PQ   (:,:)
      real(kind=dp), allocatable :: PB   (:,:)
      real(kind=dp), allocatable :: QD   (:,:)
      real(kind=dp), allocatable :: rw   (:,:)
      real(kind=dp), allocatable :: ai   (:)
      real(kind=dp), allocatable :: aj   (:)
      real(kind=dp), allocatable :: ak   (:)
      real(kind=dp), allocatable :: al   (:)
      real(kind=dp), allocatable :: fi   (:)
      real(kind=dp), allocatable :: fj   (:)
      real(kind=dp), allocatable :: fk   (:)
      real(kind=dp), allocatable :: fl   (:)
      integer :: ijklxyz(4,BAS_MXCART,4)
      real(kind=dp) :: fd(3,4)
      real(kind=dp) :: dtol
      real(kind=dp) :: dabcut
      contains
        procedure :: init => gdat_init
        procedure :: clean => gdat_clean
        procedure :: set_ids => gdat_set_ids
    end type

    logical :: dbg = .false.

    private
    public :: grd2_int_data_t
    public :: grd2_rys_compute

contains

  subroutine gdat_init(gdat, maxang, nder, &
                        dtol, dabcut, &
                        stat)

    implicit none

    class(grd2_int_data_t), intent(inout) :: gdat
    integer, intent(in) :: maxang, nder
    real(kind=dp), intent(in) :: dtol, dabcut
    integer, intent(out) :: stat
    integer :: mxbra, mxcart, mxrys

    gdat%nder = nder
    gdat%dtol = dtol
    gdat%dabcut = dabcut**2

    mxrys = (4*maxang + 2 + gdat%nder)/2
    mxcart = maxang+1 + nder
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
      gdat%f00  (mxrys*MAXCONTR*3), &
      gdat%abv  (  6, MAXCONTR  ), &
      gdat%PQ   (  3, MAXCONTR  ), &
      gdat%PB   (  3, MAXCONTR  ), &
      gdat%QD   (  3, MAXCONTR  ), &
      gdat%rw   (mxrys*2, MAXCONTR  ), &
      gdat%ai  (MAXCONTR), &
      gdat%aj  (MAXCONTR), &
      gdat%ak  (MAXCONTR), &
      gdat%al  (MAXCONTR), &

      gdat%fi   ( mxcart**4 *MAXCONTR*3), &
      gdat%fj   ( mxcart**4 *MAXCONTR*3), &
      gdat%fk   ( mxcart**4 *MAXCONTR*3), &
      gdat%fl   ( mxcart**4 *MAXCONTR*3), &
      stat=stat)
  end subroutine gdat_init

  subroutine gdat_clean(gdat)

    implicit none

    class(grd2_int_data_t), intent(inout) :: gdat
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
    if (allocated(gdat%f00  )) deallocate(gdat%f00  )
    if (allocated(gdat%abv  )) deallocate(gdat%abv  )
    if (allocated(gdat%PQ   )) deallocate(gdat%PQ   )
    if (allocated(gdat%PB   )) deallocate(gdat%PB   )
    if (allocated(gdat%QD   )) deallocate(gdat%QD   )
    if (allocated(gdat%rw   )) deallocate(gdat%rw   )
    if (allocated(gdat%ai   )) deallocate(gdat%ai   )
    if (allocated(gdat%aj   )) deallocate(gdat%aj   )
    if (allocated(gdat%ak   )) deallocate(gdat%ak   )
    if (allocated(gdat%al   )) deallocate(gdat%al   )
    if (allocated(gdat%fi   )) deallocate(gdat%fi   )
    if (allocated(gdat%fj   )) deallocate(gdat%fj   )
    if (allocated(gdat%fk   )) deallocate(gdat%fk   )
    if (allocated(gdat%fl   )) deallocate(gdat%fl   )
  end subroutine gdat_clean

  subroutine gdat_set_ids(gdat, basis, i, j, k, l)

    implicit none

    class(grd2_int_data_t), intent(inout) :: gdat
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: i, j, k, l

    integer :: id(4), am(4), flips(4), tmp(2), same

    ! Permute shells so L_i < L_j, L_k < L_l, L_i < L_j
    flips = [1,2,3,4]
    id = [i,j,k,l]
    am = basis%am([i,j,k,l])
    if (am(1) > am(2)) flips(1:2) = [2,1]
    if (am(3) > am(4)) flips(3:4) = [4,3]
    if (max(am(1),am(2)) > max(am(3),am(4))) then
      tmp = flips(1:2)
      flips(1:2) = flips(3:4)
      flips(3:4) = tmp
    end if

    gdat%id = id(flips)
    gdat%am = am(flips)

    gdat%at = basis%origin(gdat%id)

    ! Compute translation invariance class of the integral:
    same = 0
    if (gdat%at(1)/=gdat%at(2)) same = same + 32
    if (gdat%at(1)/=gdat%at(3)) same = same + 16
    if (gdat%at(1)/=gdat%at(4)) same = same + 8
    if (gdat%at(2)/=gdat%at(3)) same = same + 4
    if (gdat%at(2)/=gdat%at(4)) same = same + 2
    if (gdat%at(3)/=gdat%at(4)) same = same + 1

    select case (same)
      case(0 ); gdat%invtyp = 1  !  ( 1, 1 | 1, 1 )
      case(11); gdat%invtyp = 2  !  ( 1, 1 | 1, 2 )
      case(21); gdat%invtyp = 3  !  ( 1, 1 | 2, 1 )
      case(30); gdat%invtyp = 4  !  ( 1, 1 | 2, 2 )
      case(31); gdat%invtyp = 5  !  ( 1, 1 | 2, 3 )
      case(38); gdat%invtyp = 6  !  ( 1, 2 | 1, 1 )
      case(45); gdat%invtyp = 7  !  ( 1, 2 | 1, 2 )
      case(47); gdat%invtyp = 8  !  ( 1, 2 | 1, 3 )
      case(51); gdat%invtyp = 9  !  ( 1, 2 | 2, 1 )
      case(55); gdat%invtyp = 10 !  ( 1, 2 | 3, 1 )
      case(56); gdat%invtyp = 11 !  ( 1, 2 | 2, 2 )
      case(59); gdat%invtyp = 12 !  ( 1, 2 | 2, 3 )
      case(61); gdat%invtyp = 13 !  ( 1, 2 | 3, 2 )
      case(62); gdat%invtyp = 14 !  ( 1, 2 | 3, 3 )
      !case(63); gdat%invtyp = 15 !  ( 1, 2 | 3, 4 )
      case default; gdat%invtyp = 15
    end select

!   For debugging purposes calculate all terms
    if (dbg) gdat%invtyp = 16

    gdat%skip = skips(:,gdat%invtyp)

  end subroutine gdat_set_ids

  subroutine grd2_rys_compute(gdat, ppairs, dab, dabmax, mu2)

    use int2_pairs, only: int2_pair_storage, int2_cutoffs_t
    implicit none

    type(grd2_int_data_t) :: gdat
    type(int2_pair_storage), intent(in) :: ppairs
    real(kind=dp), intent(in) :: dabmax
    real(kind=dp), intent(in) :: dab(*)
    real(kind=dp), intent(in), optional :: mu2

    integer :: ijg, klg, maxgg, mmax, ng
    integer :: nimax, njmax, nkmax, nlmax, nmax
    real(kind=dp) :: aa, ab, aandb1, bb, da, db, test
    real(kind=dp) :: pfac, rho
    real(kind=dp) :: p(3), q(3)
    logical :: last
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
    if (npp_p*npp_q == 0) return

    nimax = gdat%am(1) + gdat%der(1) + 1
    njmax = gdat%am(2) + gdat%der(2) + 1
    nkmax = gdat%am(3) + gdat%der(3) + 1
    nlmax = gdat%am(4) + gdat%der(4) + 1

    nmax = gdat%am(1)+gdat%am(2)+1 + min(gdat%der(1)+gdat%der(2),gdat%nder)
    mmax = gdat%am(3)+gdat%am(4)+1 + min(gdat%der(3)+gdat%der(4),gdat%nder)

    maxgg = MAXCONTR/gdat%nroots

!   Pair of k,l primitives

    gdat%fd = 0
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

        if (test<gdat%dtol*ab) cycle
        if (test*dabmax*dabmax<gdat%dabcut*ab) cycle

        aandb1 = 1.0_dp/ab
        rho = aa*bb*aandb1

        ng = ng+1

        gdat%abv(1,ng) = ppairs%ginv(ppid_p-1+ijg)
        gdat%abv(2,ng) = ppairs%ginv(ppid_q-1+klg)
        gdat%abv(3,ng) = rho
        gdat%abv(4,ng) = pfac*sqrt(aandb1)
        gdat%abv(5,ng) = aandb1
        gdat%abv(6,ng) = rho*sum((p-q)**2)

        gdat%ai(ng) = 2*ppairs%alpha_a(ppid_p-1+ijg)
        gdat%aj(ng) = 2*ppairs%alpha_b(ppid_p-1+ijg)
        gdat%ak(ng) = 2*ppairs%alpha_a(ppid_q-1+klg)
        gdat%al(ng) = 2*ppairs%alpha_b(ppid_q-1+klg)

        gdat%PQ(:,ng) = p-q

        if (nmax>1) gdat%pb(:,ng) = ppairs%PB(:,ppid_p-1+ijg)
        if (mmax>1) gdat%qd(:,ng) = ppairs%PB(:,ppid_q-1+klg)

        gdat%dij(:,ng) = ppairs%PA(:,ppid_p-1+ijg)-ppairs%PB(:,ppid_p-1+ijg)
        gdat%dkl(:,ng) = ppairs%PA(:,ppid_q-1+klg)-ppairs%PB(:,ppid_q-1+klg)

        last = klg==npp_q .and. ijg==npp_p

        if (ng==maxgg .or. last) then
          if (ng==0) return

          call compute_grd_ints(gdat, dab, &
                  ng, nmax, mmax, nimax, njmax, nkmax, nlmax)

          ng = 0
        end if

      end do
    end do

!   Process derivative integrals
    call apply_translation_invariance(gdat)

  end subroutine grd2_rys_compute

  subroutine compute_grd_ints(gdat, dab, ng, nmax, mmax, nimax, njmax, nkmax, nlmax)

    type(grd2_int_data_t), intent(inout) :: gdat
    real(kind=dp), intent(in) :: dab(*)

    integer :: mmax, ng
    integer :: nimax, njmax, nkmax, nlmax, nmax

!   Compute roots and weights for quadrature
    call compute_rys_rw(gdat, gdat%rw, ng)

!   Compute coefficients for recursion formulae
    call compute_coefficients(gdat%b00, gdat%b01, gdat%b10, &
                gdat%c00, gdat%d00, gdat%f00, &
                gdat%abv, gdat%pq, gdat%pb, gdat%qd, gdat%rw, nmax, mmax, ng, gdat%nroots)

!   Compute x, y, z integrals (2 centers, 2-d )
    call compute_xyz_p0q0(gdat%gnm,ng*gdat%nroots,nmax,mmax, &
                gdat%b00, gdat%b01, gdat%b10, gdat%c00, gdat%d00, gdat%f00)

!   Compute x, y, z integrals (4 centers, 2-d)
    call compute_xyz_ijkl(gdat%gijkl, gdat%gnkl, gdat%gnm, &
                ng, gdat%nroots, nmax, mmax, nimax, njmax, nkmax,nlmax,   &
                gdat%dij,gdat%dkl)

!   Compute x, y, z integrals for derivatives
    call compute_der_xyz_ijkl(gdat, gdat%gijkl, &
                ng, gdat%nroots*3, nimax, njmax, nkmax, nlmax, &
                gdat%ai, gdat%aj, gdat%ak, gdat%al, gdat%fi, gdat%fj, gdat%fk, gdat%fl)

!   compute derivative integrals
    call compute_der_ijkl(gdat, ng*gdat%nroots, gdat%ijklxyz, gdat%gijkl, &
                gdat%fi, gdat%fj, gdat%fk, gdat%fl, dab, gdat%fd)
  end subroutine

  subroutine set_shells(gdat)

    implicit none

    type(grd2_int_data_t) :: gdat
    integer :: ish, jsh, ksh, lsh

    ish = gdat%id(1)
    jsh = gdat%id(2)
    ksh = gdat%id(3)
    lsh = gdat%id(4)

    gdat%iandj = ish == jsh
    gdat%kandl = ksh == lsh
    gdat%same = ish == ksh.and.jsh == lsh

    gdat%nbf = num_cart_bf(gdat%am)

    gdat%der(:) = gdat%nder
    if (gdat%skip(1)) gdat%der(1) = 0
    if (gdat%skip(2)) gdat%der(2) = 0
    if (gdat%skip(3)) gdat%der(3) = 0
    if (gdat%skip(4)) gdat%der(4) = 0

!   Set number of quadrature points
    gdat%nroots = (sum(gdat%am)+2 + gdat%nder )/2

!   Prepare indices for pairs of (i,j) functions
    call prepare_xyz_ids(gdat)

  end subroutine set_shells

  subroutine prepare_xyz_ids(gdat)
    implicit none
    class(grd2_int_data_t), intent(inout) :: gdat
    integer :: i, nj, nk, nl, njkl, nkl

    nj = gdat%am(2) + gdat%der(2) + 1
    nk = gdat%am(3) + gdat%der(3) + 1
    nl = gdat%am(4) + gdat%der(4) + 1
    njkl = nl*nk*nj
    do i = 1, gdat%nbf(1)
      gdat%ijklxyz(1,i,1) = cart_x(i,gdat%am(1))*njkl
      gdat%ijklxyz(2,i,1) = cart_y(i,gdat%am(1))*njkl
      gdat%ijklxyz(3,i,1) = cart_z(i,gdat%am(1))*njkl
    end do

    nkl = nl*nk
    do i = 1, gdat%nbf(2)
      gdat%ijklxyz(1,i,2) = cart_x(i,gdat%am(2))*nkl
      gdat%ijklxyz(2,i,2) = cart_y(i,gdat%am(2))*nkl
      gdat%ijklxyz(3,i,2) = cart_z(i,gdat%am(2))*nkl
    end do

!   Prepare indices for pairs of (k,l) functions
    do i = 1, gdat%nbf(3)
      gdat%ijklxyz(1,i,3) = cart_x(i,gdat%am(3))*nl
      gdat%ijklxyz(2,i,3) = cart_y(i,gdat%am(3))*nl
      gdat%ijklxyz(3,i,3) = cart_z(i,gdat%am(3))*nl
    end do

    do i = 1, gdat%nbf(4)
      gdat%ijklxyz(1,i,4) = cart_x(i,gdat%am(4))+1
      gdat%ijklxyz(2,i,4) = cart_y(i,gdat%am(4))+1
      gdat%ijklxyz(3,i,4) = cart_z(i,gdat%am(4))+1
    end do
  end subroutine

  subroutine compute_rys_rw(gdat, rwv, numg)
    use rys, only: rys_root_t
    implicit none


    class(grd2_int_data_t), intent(inout) :: gdat
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

  subroutine compute_coefficients(b00,b01,b10,c00,d00,f00,abv,pq,pb,qd,rwv,nmax,mmax,numg,nroots)

    implicit none

    real(kind=dp) :: b00(numg,*),b01(numg,*),b10(numg,*) !numg*nroots
    real(kind=dp) :: c00(numg,nroots,*) !numg*nroots*3
    real(kind=dp) :: d00(numg,nroots,*) !numg*nroots*3
    real(kind=dp) :: f00(numg,nroots,*) !numg*nroots*3

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

        f00(ng,nr,1) = ww*pfac
        f00(ng,nr,2) = 1.0_dp
        f00(ng,nr,3) = 1.0_dp

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

  subroutine compute_xyz_p0q0(gnm,ng,nmax,mmax,b00,b01,b10,c00,d00,f00)

    implicit none

    real(kind=dp) :: gnm(ng,3,nmax,*)
    real(kind=dp) :: c00(ng,*),d00(ng,*),f00(ng,*)
    real(kind=dp) :: b00(*),b01(*),b10(*)
    integer :: ng, nmax, mmax

    integer :: m, n, xyz

!   G(0,0)
    gnm(:ng,:,1,1) = f00(:,:3)
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

  subroutine compute_xyz_ijkl(ijkl,gnkl,gnm,ng,nr,nmax,mmax,nimax,njmax,nkmax,nlmax,dij,dkl)

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
      if(nk == nkmax) cycle
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
      if (ni == nimax) cycle
      n1 = nmax-ni
      do xyz = 1, 3
        do ig = 1, ng
          gnkl(ig,:,xyz,:,:,1:n1) = dij(xyz,ig)*gnkl(ig,:,xyz,:,:,1:n1) &
                                  +             gnkl(ig,:,xyz,:,:,2:n1+1)
        end do
      end do
    end do

  end subroutine compute_xyz_ijkl

  subroutine compute_der_xyz_ijkl(gdat,g, &
       ng,nr3,nimax,njmax,nkmax,nlmax,aai,aaj,aak,aal,fi,fj,fk,fl)

      implicit none

      type(grd2_int_data_t) :: gdat
      integer :: ng, nr3, nimax, njmax, nkmax, nlmax
      real(kind=dp) :: g(ng,nr3,nlmax,nkmax,njmax,*)
      real(kind=dp) ::  aai(*) !   aai(ng)
      real(kind=dp) ::  aaj(*) !   aaj(ng)
      real(kind=dp) ::  aak(*) !   aak(ng)
      real(kind=dp) ::  aal(*) !   aal(ng)
      real(kind=dp) :: fi(ng,nr3,nlmax,nkmax,njmax,*) ! fi(ng,nlmax,nkmax,njmax,nimax)
      real(kind=dp) :: fj(ng,nr3,nlmax,nkmax,njmax,*) ! fj(ng,nlmax,nkmax,njmax,nimax)
      real(kind=dp) :: fk(ng,nr3,nlmax,nkmax,njmax,*) ! fk(ng,nlmax,nkmax,njmax,nimax)
      real(kind=dp) :: fl(ng,nr3,nlmax,nkmax,njmax,*) ! fl(ng,nlmax,nkmax,njmax,nimax)

      integer :: i, j, k, l, n
      integer :: ni, nj, nk, nl

      ni = gdat%am(1) + 1
      nj = gdat%am(2) + 1
      nk = gdat%am(3) + 1
      nl = gdat%am(4) + 1

!     First derivatives only
      if (.not.gdat%skip(1)) then

!       FI only
        do n = 1, ng
          fi(n,:,:,:,:,1) = g(n,:,:,:,:,2)*aai(n)
        end do

        if (ni/=1) then
          do i = 2, ni
            do n = 1, ng
              fi(n,:,:,:,:,i)= g(n,:,:,:,:,i+1)*aai(n) &
                             - g(n,:,:,:,:,i-1)*(i-1)
             end do
          end do
        end if

      end if

      if (.not.gdat%skip(2)) then

!       FJ only
        do i = 1, nimax
          do n = 1, ng
            fj(n,:,:,:,1,i) = g(n,:,:,:,2,i)*aaj(n)
          end do
        end do

        if (nj/=1) then

          do i = 1, nimax
            do j = 2, nj
              do n = 1, ng
                fj(n,:,:,:,j,i)= g(n,:,:,:,j+1,i)*aaj(n) &
                               - g(n,:,:,:,j-1,i)*(j-1)
              end do
            end do
          end do
        end if

      end if

      if (.not.gdat%skip(3)) then

!       FK only
        do i = 1, nimax
          do n = 1, ng
            fk(n,:,:,1,:,i) = g(n,:,:,2,:,i)*aak(n)
          end do
        end do

        if (nk/=1) then

          do i = 1, nimax
            do k = 2, nk
              do n = 1, ng
                fk(n,:,:,k,:,i) = g(n,:,:,k+1,:,i)*aak(n) &
                                - g(n,:,:,k-1,:,i)*(k-1)
              end do
            end do
          end do
        end if

      end if

      if (.not.gdat%skip(4)) then

!       FL and SLL
        do i = 1, nimax
          do n = 1, ng
            fl(n,:,1,:,:,i) = g(n,:,2,:,:,i)*aal(n)
          end do
        end do

        if (nl/=1) then
        do i = 1, nimax
          do l = 2, nl
            do n = 1, ng
              fl(n,:,l,:,:,i) = g(n,:,l+1,:,:,i)*aal(n) &
                              - g(n,:,l-1,:,:,i)*(l-1)
            end do
          end do
        end do
        end if

      end if

  end subroutine compute_der_xyz_ijkl

  subroutine compute_der_ijkl(gdat,ngnr,&
                  ijklxyz,g0,fi,fj,fk,fl,den,fd)

    implicit none

    type(grd2_int_data_t) :: gdat
    integer :: ngnr
    integer :: ijklxyz(:,:,:)
    real(kind=dp), target :: den(*)
    real(kind=dp) ::  g0(ngnr,3,*)
    real(kind=dp) ::  fi(ngnr,3,*), fj(ngnr,3,*)
    real(kind=dp) ::  fk(ngnr,3,*), fl(ngnr,3,*)
    real(kind=dp) :: fd(3,4)

    integer :: i, j, k, l
    integer :: nx, ny, nz
    real(kind=dp), pointer :: pd(:,:,:,:)

    pd(1:gdat%nbf(4), 1:gdat%nbf(3), 1:gdat%nbf(2), 1:gdat%nbf(1)) => den(1:product(gdat%nbf))

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
                      , df => pd(l,k,j,i) &
                      )

              if (.not.gdat%skip(1)) then
                fd(1,1) = fd(1,1) + df * sum(fi(:,1,nx)*y*z)
                fd(2,1) = fd(2,1) + df * sum(fi(:,2,ny)*x*z)
                fd(3,1) = fd(3,1) + df * sum(fi(:,3,nz)*x*y)
              end if
              if (.not.gdat%skip(2)) then
                fd(1,2) = fd(1,2) + df * sum(fj(:,1,nx)*y*z)
                fd(2,2) = fd(2,2) + df * sum(fj(:,2,ny)*x*z)
                fd(3,2) = fd(3,2) + df * sum(fj(:,3,nz)*x*y)
              end if
              if (.not.gdat%skip(3)) then
                fd(1,3) = fd(1,3) + df * sum(fk(:,1,nx)*y*z)
                fd(2,3) = fd(2,3) + df * sum(fk(:,2,ny)*x*z)
                fd(3,3) = fd(3,3) + df * sum(fk(:,3,nz)*x*y)
              end if
              if (.not.gdat%skip(4)) then
                fd(1,4) = fd(1,4) + df * sum(fl(:,1,nx)*y*z)
                fd(2,4) = fd(2,4) + df * sum(fl(:,2,ny)*x*z)
                fd(3,4) = fd(3,4) + df * sum(fl(:,3,nz)*x*y)
              end if
            end associate

          end do
        end do
      end do
    end do

  end subroutine compute_der_ijkl

  subroutine apply_translation_invariance(gdat)

    implicit none

    type(grd2_int_data_t) :: gdat

    if (gdat%nder == 0) return

!   Translational invariance for gradient elements
    associate(fd => gdat%fd)

    if (gdat%iandj) fd = 0.5*fd
    if (gdat%kandl) fd = 0.5*fd
    if (gdat%same) fd = 0.5*fd

    select case (gdat%invtyp)
    case (2)    ; fd(:,1) = - fd(:,4)
    case (3)    ; fd(:,1) = - fd(:,3)
    case (4,5)  ; fd(:,1) = -(fd(:,3)+fd(:,4))
    case (6)    ; fd(:,1) = - fd(:,2)
    case (7,8)  ; fd(:,1) = -(fd(:,2)+fd(:,4))
    case (9,10) ; fd(:,1) = -(fd(:,2)+fd(:,3))
    case (11)   ; fd(:,2) = - fd(:,1)
    case (12)   ; fd(:,2) = -(fd(:,1)+fd(:,4))
    case (13)   ; fd(:,2) = -(fd(:,1)+fd(:,3))
    case (14)   ; fd(:,3) = -(fd(:,1)+fd(:,2))
    case (15)   ; fd(:,4) = -(fd(:,1)+fd(:,2)+fd(:,3))
    end select
end associate

  end subroutine apply_translation_invariance
end module grd2_rys
