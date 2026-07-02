module grd2_rys
    use precision, only: dp
    use basis_tools, only: basis_set
    use constants, only: bas_mxang, bas_mxcart, num_cart_bf, cart_x, cart_y, cart_z

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
      real(kind=dp) :: fd2(3,4,3,4)
      ! Second-derivative 1D integral arrays (allocated only when nder>=2).
      ! f2tmp is scratch; f2_cc' holds d^2/d(center c)d(center c') of the
      ! 4-center 1D integrals, stored in the same layout as fi/fj/fk/fl.
      real(kind=dp), allocatable :: f2tmp(:)
      real(kind=dp), allocatable :: f2_11(:), f2_12(:), f2_13(:), f2_14(:)
      real(kind=dp), allocatable :: f2_22(:), f2_23(:), f2_24(:)
      real(kind=dp), allocatable :: f2_33(:), f2_34(:)
      real(kind=dp), allocatable :: f2_44(:)
      real(kind=dp) :: dtol
      real(kind=dp) :: dabcut
      contains
        procedure :: init => gdat_init
        procedure :: clean => gdat_clean
        procedure :: set_ids => gdat_set_ids
    end type

type soc2e_int_data_t
      integer :: id(4)
      integer :: at(4)
      integer :: am(4)
      integer :: nbf(4)
      integer :: nder = 1   ! derivatives only on IJ (electron 1)
      integer :: nroots
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
      integer :: ijklxyz(4,BAS_MXCART,4)
      integer :: ao_offset(4)
      real(kind=dp) :: dtol
      real(kind=dp) :: dabcut
    contains
      procedure :: init => soc2e_gdat_init
      procedure :: clean => soc2e_gdat_clean
      procedure :: set_ids => soc2e_gdat_set_ids
    end type

    logical :: dbg = .false.

    private
    public :: grd2_int_data_t
    public :: grd2_rys_compute
    public :: soc2e_int_data_t
    public :: soc2e_rys_compute
    public :: soc2e_driver
    public :: grd2_rys_hess_compute

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

!   Second-derivative work arrays (Hessian path only)
    if (nder >= 2) then
      allocate(&
        gdat%f2tmp( mxcart**4 *MAXCONTR*3), &
        gdat%f2_11( mxcart**4 *MAXCONTR*3), &
        gdat%f2_12( mxcart**4 *MAXCONTR*3), &
        gdat%f2_13( mxcart**4 *MAXCONTR*3), &
        gdat%f2_14( mxcart**4 *MAXCONTR*3), &
        gdat%f2_22( mxcart**4 *MAXCONTR*3), &
        gdat%f2_23( mxcart**4 *MAXCONTR*3), &
        gdat%f2_24( mxcart**4 *MAXCONTR*3), &
        gdat%f2_33( mxcart**4 *MAXCONTR*3), &
        gdat%f2_34( mxcart**4 *MAXCONTR*3), &
        gdat%f2_44( mxcart**4 *MAXCONTR*3), &
        stat=stat)
    end if
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
    if (allocated(gdat%f2tmp)) deallocate(gdat%f2tmp)
    if (allocated(gdat%f2_11)) deallocate(gdat%f2_11)
    if (allocated(gdat%f2_12)) deallocate(gdat%f2_12)
    if (allocated(gdat%f2_13)) deallocate(gdat%f2_13)
    if (allocated(gdat%f2_14)) deallocate(gdat%f2_14)
    if (allocated(gdat%f2_22)) deallocate(gdat%f2_22)
    if (allocated(gdat%f2_23)) deallocate(gdat%f2_23)
    if (allocated(gdat%f2_24)) deallocate(gdat%f2_24)
    if (allocated(gdat%f2_33)) deallocate(gdat%f2_33)
    if (allocated(gdat%f2_34)) deallocate(gdat%f2_34)
    if (allocated(gdat%f2_44)) deallocate(gdat%f2_44)
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

  !> @brief Compute Rys quadrature roots and weights for the 2e SOC integrals
  !> @details
  !>  Wrapper around the Rys root/weight evaluator for the SOC case.
  !>  The number of roots is increased by one relative to the standard 2e case
  !>  to accommodate the higher-order numerator arising from the angular momentum
  !>  operator: ⟨μν|r^{-3} L|λσ⟩ ~ ⟨∂μ|r^{-1}|∂σ⟩.
  !>
  !> @param[inout] gdat   SOC integral data structure (nroots is read from gdat)
  !> @param[out]   rwv    Roots and weights array (2*nroots x ng)
  !> @param[in]    numg   Number of primitive pairs in the current batch
  subroutine soc2e_compute_rys_rw(gdat, rwv, numg)
    use rys, only: rys_root_t
    implicit none

    type(soc2e_int_data_t), intent(inout) :: gdat
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

  end subroutine soc2e_compute_rys_rw

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

    integer :: ni, nk, nl, ig, m1, n1, xyz, n, m, ir
    real(kind=dp) :: dt(ng,3)

!   Transposed pair-distance coefficients: contiguous along the batch index,
!   so the recursions below run on stride-1 vectors of length ng.
    do ig = 1, ng
      dt(ig,1) = dkl(1,ig)
      dt(ig,2) = dkl(2,ig)
      dt(ig,3) = dkl(3,ig)
    end do

!   g(n,k,l)
    do nk=1, nkmax
      do nl=1,nlmax
        gnkl(:,:,:,nl,nk,:nmax) = gnm(:,:,:,:,nl)
      end do
      if(nk == nkmax) cycle
      m1 = mmax-nk
!     m ascending: iteration m reads column m+1 before it is overwritten,
!     reproducing the whole-array statement semantics.
      do m = 1, m1
        do n = 1, nmax
          do xyz = 1, 3
            do ir = 1, nr
              gnm(:,ir,xyz,n,m) = dt(:,xyz)*gnm(:,ir,xyz,n,m) &
                                +           gnm(:,ir,xyz,n,m+1)
            end do
          end do
        end do
      end do
    end do

    do ig = 1, ng
      dt(ig,1) = dij(1,ig)
      dt(ig,2) = dij(2,ig)
      dt(ig,3) = dij(3,ig)
    end do

!   g(i,j,k,l)
    do ni = 1, nimax
      ijkl(:,:,:,:,:,1:njmax,ni) = gnkl(:,:,:,:,:,1:njmax)
      if (ni == nimax) cycle
      n1 = nmax-ni
      do n = 1, n1
        do nk = 1, nkmax
          do nl = 1, nlmax
            do xyz = 1, 3
              do ir = 1, nr
                gnkl(:,ir,xyz,nl,nk,n) = dt(:,xyz)*gnkl(:,ir,xyz,nl,nk,n) &
                                       +           gnkl(:,ir,xyz,nl,nk,n+1)
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine compute_xyz_ijkl

  subroutine compute_der_xyz_ijkl(gdat,g, &
       ng,nr3,nimax,njmax,nkmax,nlmax,aai,aaj,aak,aal,fi,fj,fk,fl)

      implicit none

      type(grd2_int_data_t) :: gdat
      integer :: ng, nr3, nimax, njmax, nkmax, nlmax
      real(kind=dp) :: g(*)
      real(kind=dp) ::  aai(*) !   aai(ng)
      real(kind=dp) ::  aaj(*) !   aaj(ng)
      real(kind=dp) ::  aak(*) !   aak(ng)
      real(kind=dp) ::  aal(*) !   aal(ng)
      real(kind=dp) :: fi(*), fj(*), fk(*), fl(*)

      integer :: ni, nj, nk, nl

      ni = gdat%am(1) + 1
      nj = gdat%am(2) + 1
      nk = gdat%am(3) + 1
      nl = gdat%am(4) + 1

!     Apply the first-derivative operator
!         d/d(center) = 2*alpha * raise - power * lower
!     along each center's 1D-integral index. The shared kernel views the
!     (ng,nr3,nlmax,nkmax,njmax,nimax) arrays as (ng, mid, axis, post) so all
!     updates run on contiguous stride-1 vectors of length ng.

      if (.not.gdat%skip(1)) &  ! FI: axis = i (outermost)
        call deriv_axis_1d(g, fi, ng, nr3*nlmax*nkmax*njmax, nimax, 1, &
                           ni, aai)

      if (.not.gdat%skip(2)) &  ! FJ: axis = j
        call deriv_axis_1d(g, fj, ng, nr3*nlmax*nkmax, njmax, nimax, &
                           nj, aaj)

      if (.not.gdat%skip(3)) &  ! FK: axis = k
        call deriv_axis_1d(g, fk, ng, nr3*nlmax, nkmax, njmax*nimax, &
                           nk, aak)

      if (.not.gdat%skip(4)) &  ! FL: axis = l (innermost after nr3)
        call deriv_axis_1d(g, fl, ng, nr3, nlmax, nkmax*njmax*nimax, &
                           nl, aal)

  end subroutine compute_der_xyz_ijkl

!> @brief Apply the single-center first-derivative operator along one axis of
!>        a 4-center 1D-integral array, vectorized over the contiguous
!>        primitive-batch index.
!> @details g and f are interpreted as (ng, mid, nax, npost), where `mid` and
!>        `npost` are the products of the extents below/above the derivative
!>        axis. f(:,m,a,p) = g(:,m,a+1,p)*aa - (a-1)*g(:,m,a-1,p) for
!>        a = 1..nphys (the a=1 lowering term vanishes).
  subroutine deriv_axis_1d(g, f, ng, mid, nax, npost, nphys, aa)
      implicit none
      integer, intent(in) :: ng, mid, nax, npost, nphys
      real(kind=dp), intent(in)  :: g(ng, mid, nax, npost)
      real(kind=dp), intent(out) :: f(ng, mid, nax, npost)
      real(kind=dp), intent(in)  :: aa(ng)

      integer :: m, a, p

      do p = 1, npost
        do m = 1, mid
          f(:,m,1,p) = g(:,m,2,p)*aa
        end do
        do a = 2, nphys
          do m = 1, mid
            f(:,m,a,p) = g(:,m,a+1,p)*aa - g(:,m,a-1,p)*(a-1)
          end do
        end do
      end do

  end subroutine deriv_axis_1d

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
    real(kind=dp) :: yz(ngnr), xz(ngnr), xy(ngnr)

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

!             Hoist the pair products: they are shared by all four centers.
              yz = y*z
              xz = x*z
              xy = x*y

              if (.not.gdat%skip(1)) then
                fd(1,1) = fd(1,1) + df * sum(fi(:,1,nx)*yz)
                fd(2,1) = fd(2,1) + df * sum(fi(:,2,ny)*xz)
                fd(3,1) = fd(3,1) + df * sum(fi(:,3,nz)*xy)
              end if
              if (.not.gdat%skip(2)) then
                fd(1,2) = fd(1,2) + df * sum(fj(:,1,nx)*yz)
                fd(2,2) = fd(2,2) + df * sum(fj(:,2,ny)*xz)
                fd(3,2) = fd(3,2) + df * sum(fj(:,3,nz)*xy)
              end if
              if (.not.gdat%skip(3)) then
                fd(1,3) = fd(1,3) + df * sum(fk(:,1,nx)*yz)
                fd(2,3) = fd(2,3) + df * sum(fk(:,2,ny)*xz)
                fd(3,3) = fd(3,3) + df * sum(fk(:,3,nz)*xy)
              end if
              if (.not.gdat%skip(4)) then
                fd(1,4) = fd(1,4) + df * sum(fl(:,1,nx)*yz)
                fd(2,4) = fd(2,4) + df * sum(fl(:,2,ny)*xz)
                fd(3,4) = fd(3,4) + df * sum(fl(:,3,nz)*xy)
              end if
            end associate

          end do
        end do
      end do
    end do

  end subroutine compute_der_ijkl

!###############################################################################
!   Second-derivative (Hessian) skeleton: analytic 2e ERI second derivatives
!###############################################################################

  subroutine grd2_rys_hess_compute(gdat, ppairs, dab, dabmax, mu2)
!   Mirror of grd2_rys_compute, but accumulates the per-quartet second
!   derivative block gdat%fd2(3,4,3,4) instead of the gradient gdat%fd(3,4).
!   All four centers are differentiated explicitly (no translation-invariance
!   recovery); the only post-processing is the permutational multiplicity
!   halving identical to the gradient path.
    use int2_pairs, only: int2_pair_storage
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

    mu2_1 = 0
    if (present(mu2)) mu2_1 = 1.0d0/mu2

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

    gdat%fd2 = 0
    ng = 0

    do klg = 1, npp_q
      db = ppairs%k(ppid_q-1+klg)*ppairs%ginv(ppid_q-1+klg)
      bb = ppairs%g(ppid_q-1+klg)
      q = ppairs%P(:,ppid_q-1+klg)

      do ijg = 1, npp_p
        da = ppairs%k(ppid_p-1+ijg)*ppairs%ginv(ppid_p-1+ijg)
        aa = ppairs%g(ppid_p-1+ijg)
        p = ppairs%P(:,ppid_p-1+ijg)

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
          call compute_grd2_ints(gdat, dab, &
                  ng, nmax, mmax, nimax, njmax, nkmax, nlmax)
          ng = 0
        end if

      end do
    end do

!   Permutational multiplicity (same correction as the gradient path)
    if (gdat%iandj) gdat%fd2 = 0.5_dp*gdat%fd2
    if (gdat%kandl) gdat%fd2 = 0.5_dp*gdat%fd2
    if (gdat%same)  gdat%fd2 = 0.5_dp*gdat%fd2

  end subroutine grd2_rys_hess_compute

  subroutine compute_grd2_ints(gdat, dab, ng, nmax, mmax, nimax, njmax, nkmax, nlmax)
    type(grd2_int_data_t), intent(inout) :: gdat
    real(kind=dp), intent(in) :: dab(*)
    integer :: mmax, ng
    integer :: nimax, njmax, nkmax, nlmax, nmax

    call compute_rys_rw(gdat, gdat%rw, ng)

    call compute_coefficients(gdat%b00, gdat%b01, gdat%b10, &
                gdat%c00, gdat%d00, gdat%f00, &
                gdat%abv, gdat%pq, gdat%pb, gdat%qd, gdat%rw, nmax, mmax, ng, gdat%nroots)

    call compute_xyz_p0q0(gdat%gnm,ng*gdat%nroots,nmax,mmax, &
                gdat%b00, gdat%b01, gdat%b10, gdat%c00, gdat%d00, gdat%f00)

    call compute_xyz_ijkl(gdat%gijkl, gdat%gnkl, gdat%gnm, &
                ng, gdat%nroots, nmax, mmax, nimax, njmax, nkmax,nlmax,   &
                gdat%dij,gdat%dkl)

!   First derivatives (all four centers)
    call compute_der_xyz_ijkl(gdat, gdat%gijkl, &
                ng, gdat%nroots*3, nimax, njmax, nkmax, nlmax, &
                gdat%ai, gdat%aj, gdat%ak, gdat%al, gdat%fi, gdat%fj, gdat%fk, gdat%fl)

!   Second-derivative 1D arrays (10 unique center pairs)
    call compute_der2_xyz_ijkl(gdat, &
                ng, gdat%nroots*3, nimax, njmax, nkmax, nlmax)

!   Contract with density into the per-quartet second-derivative block
    call compute_der2_ijkl(gdat, ng*gdat%nroots, gdat%ijklxyz, gdat%gijkl, &
                gdat%fi, gdat%fj, gdat%fk, gdat%fl, &
                gdat%f2_11, gdat%f2_12, gdat%f2_13, gdat%f2_14, &
                gdat%f2_22, gdat%f2_23, gdat%f2_24, &
                gdat%f2_33, gdat%f2_34, gdat%f2_44, &
                dab, gdat%fd2)
  end subroutine compute_grd2_ints

  subroutine der_center(src, dst, ng, nr3, nlmax, nkmax, njmax, nimax, &
                        active, aa, imax, jmax, kmax, lmax)
!   Apply the single-center first-derivative operator
!       d/d(center) = 2*alpha * raise  -  power * lower
!   to a 4-center 1D integral array (g layout), filling dst over the
!   requested physical index ranges. Composing this operator twice (or on
!   two different centers) yields the second-derivative arrays.
    implicit none
    integer, intent(in) :: ng, nr3, nlmax, nkmax, njmax, nimax
    integer, intent(in) :: active, imax, jmax, kmax, lmax
    real(kind=dp), intent(in)  :: src(ng,nr3,nlmax,nkmax,njmax,nimax)
    real(kind=dp), intent(out) :: dst(ng,nr3,nlmax,nkmax,njmax,nimax)
    real(kind=dp), intent(in)  :: aa(*)
    integer :: i, j, k, l, n

    select case (active)
    case (1)  ! d/d center 1 (i index)
      do i = 1, imax
       do j = 1, jmax
        do k = 1, kmax
         do l = 1, lmax
          if (i==1) then
            do n = 1, ng
              dst(n,:,l,k,j,1) = src(n,:,l,k,j,2)*aa(n)
            end do
          else
            do n = 1, ng
              dst(n,:,l,k,j,i) = src(n,:,l,k,j,i+1)*aa(n) - src(n,:,l,k,j,i-1)*(i-1)
            end do
          end if
         end do
        end do
       end do
      end do
    case (2)  ! d/d center 2 (j index)
      do i = 1, imax
       do j = 1, jmax
        do k = 1, kmax
         do l = 1, lmax
          if (j==1) then
            do n = 1, ng
              dst(n,:,l,k,1,i) = src(n,:,l,k,2,i)*aa(n)
            end do
          else
            do n = 1, ng
              dst(n,:,l,k,j,i) = src(n,:,l,k,j+1,i)*aa(n) - src(n,:,l,k,j-1,i)*(j-1)
            end do
          end if
         end do
        end do
       end do
      end do
    case (3)  ! d/d center 3 (k index)
      do i = 1, imax
       do j = 1, jmax
        do k = 1, kmax
         do l = 1, lmax
          if (k==1) then
            do n = 1, ng
              dst(n,:,l,1,j,i) = src(n,:,l,2,j,i)*aa(n)
            end do
          else
            do n = 1, ng
              dst(n,:,l,k,j,i) = src(n,:,l,k+1,j,i)*aa(n) - src(n,:,l,k-1,j,i)*(k-1)
            end do
          end if
         end do
        end do
       end do
      end do
    case (4)  ! d/d center 4 (l index)
      do i = 1, imax
       do j = 1, jmax
        do k = 1, kmax
         do l = 1, lmax
          if (l==1) then
            do n = 1, ng
              dst(n,:,1,k,j,i) = src(n,:,2,k,j,i)*aa(n)
            end do
          else
            do n = 1, ng
              dst(n,:,l,k,j,i) = src(n,:,l+1,k,j,i)*aa(n) - src(n,:,l-1,k,j,i)*(l-1)
            end do
          end if
         end do
        end do
       end do
      end do
    end select
  end subroutine der_center

  subroutine compute_der2_xyz_ijkl(gdat, ng, nr3, nimax, njmax, nkmax, nlmax)
!   Build the 10 unique center-pair second-derivative 1D integral arrays by
!   composing the single-center first-derivative operator (der_center).
    implicit none
    type(grd2_int_data_t) :: gdat
    integer, intent(in) :: ng, nr3, nimax, njmax, nkmax, nlmax
    integer :: ni, nj, nk, nl

    ni = gdat%am(1) + 1
    nj = gdat%am(2) + 1
    nk = gdat%am(3) + 1
    nl = gdat%am(4) + 1

!   Same-center second derivatives: apply the operator twice (via scratch).
!   f2_11 = d/di ( d/di g )
    call der_center(gdat%gijkl, gdat%f2tmp, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    1, gdat%ai, ni+1, nj, nk, nl)
    call der_center(gdat%f2tmp, gdat%f2_11, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    1, gdat%ai, ni, nj, nk, nl)
!   f2_22 = d/dj ( d/dj g )
    call der_center(gdat%gijkl, gdat%f2tmp, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    2, gdat%aj, ni, nj+1, nk, nl)
    call der_center(gdat%f2tmp, gdat%f2_22, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    2, gdat%aj, ni, nj, nk, nl)
!   f2_33 = d/dk ( d/dk g )
    call der_center(gdat%gijkl, gdat%f2tmp, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    3, gdat%ak, ni, nj, nk+1, nl)
    call der_center(gdat%f2tmp, gdat%f2_33, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    3, gdat%ak, ni, nj, nk, nl)
!   f2_44 = d/dl ( d/dl g )
    call der_center(gdat%gijkl, gdat%f2tmp, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    4, gdat%al, ni, nj, nk, nl+1)
    call der_center(gdat%f2tmp, gdat%f2_44, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    4, gdat%al, ni, nj, nk, nl)

!   Mixed-center second derivatives: apply the second center's operator to the
!   already-formed first-derivative array (fi/fj/fk/fl).
    call der_center(gdat%fi, gdat%f2_12, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    2, gdat%aj, ni, nj, nk, nl)
    call der_center(gdat%fi, gdat%f2_13, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    3, gdat%ak, ni, nj, nk, nl)
    call der_center(gdat%fi, gdat%f2_14, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    4, gdat%al, ni, nj, nk, nl)
    call der_center(gdat%fj, gdat%f2_23, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    3, gdat%ak, ni, nj, nk, nl)
    call der_center(gdat%fj, gdat%f2_24, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    4, gdat%al, ni, nj, nk, nl)
    call der_center(gdat%fk, gdat%f2_34, ng, nr3, nlmax, nkmax, njmax, nimax, &
                    4, gdat%al, ni, nj, nk, nl)
  end subroutine compute_der2_xyz_ijkl

  subroutine compute_der2_ijkl(gdat, ngnr, ijklxyz, g0, fi, fj, fk, fl, &
                  f2_11, f2_12, f2_13, f2_14, f2_22, f2_23, f2_24, &
                  f2_33, f2_34, f2_44, den, fd2)
!   Contract the 1D first/second derivative arrays with the 2-body density to
!   form the per-quartet second-derivative block fd2(a1,c1,a2,c2).
    implicit none
    type(grd2_int_data_t) :: gdat
    integer :: ngnr
    integer :: ijklxyz(:,:,:)
    real(kind=dp), target :: den(*)
    real(kind=dp) :: g0(ngnr,3,*)
    real(kind=dp) :: fi(ngnr,3,*), fj(ngnr,3,*), fk(ngnr,3,*), fl(ngnr,3,*)
    real(kind=dp) :: f2_11(ngnr,3,*), f2_12(ngnr,3,*), f2_13(ngnr,3,*)
    real(kind=dp) :: f2_14(ngnr,3,*), f2_22(ngnr,3,*), f2_23(ngnr,3,*)
    real(kind=dp) :: f2_24(ngnr,3,*), f2_33(ngnr,3,*), f2_34(ngnr,3,*)
    real(kind=dp) :: f2_44(ngnr,3,*)
    real(kind=dp) :: fd2(3,4,3,4)

    integer :: i, j, k, l
    integer :: noff(3)
    integer :: c1, c2, a1, a2, a3, o1, o2
    real(kind=dp) :: df, val
    real(kind=dp), pointer :: pd(:,:,:,:)

    pd(1:gdat%nbf(4), 1:gdat%nbf(3), 1:gdat%nbf(2), 1:gdat%nbf(1)) => den(1:product(gdat%nbf))

    do i = 1, gdat%nbf(1)
      do j = 1, gdat%nbf(2)
        do k = 1, gdat%nbf(3)
          do l = 1, gdat%nbf(4)
            noff(1) = ijklxyz(1,i,1)+ijklxyz(1,j,2)+ijklxyz(1,k,3)+ijklxyz(1,l,4)
            noff(2) = ijklxyz(2,i,1)+ijklxyz(2,j,2)+ijklxyz(2,k,3)+ijklxyz(2,l,4)
            noff(3) = ijklxyz(3,i,1)+ijklxyz(3,j,2)+ijklxyz(3,k,3)+ijklxyz(3,l,4)
            df = pd(l,k,j,i)

            do c1 = 1, 4
             do a1 = 1, 3
              do c2 = 1, 4
               do a2 = 1, 3
                if (a1 == a2) then
                  ! both derivatives act on the same 1D direction factor
                  o1 = mod(a1,  3) + 1
                  o2 = mod(a1+1,3) + 1
                  val = df * sum( f2pick(c1,c2,a1,noff(a1)) &
                                * g0(:,o1,noff(o1)) * g0(:,o2,noff(o2)) )
                else
                  ! different direction factors: product of first derivatives
                  a3 = 6 - a1 - a2
                  val = df * sum( f1pick(c1,a1,noff(a1)) &
                                * f1pick(c2,a2,noff(a2)) * g0(:,a3,noff(a3)) )
                end if
                fd2(a1,c1,a2,c2) = fd2(a1,c1,a2,c2) + val
               end do
              end do
             end do
            end do

          end do
        end do
      end do
    end do

  contains

    function f1pick(c, d, o) result(v)
      integer, intent(in) :: c, d, o
      real(kind=dp) :: v(ngnr)
      select case (c)
      case (1); v = fi(:,d,o)
      case (2); v = fj(:,d,o)
      case (3); v = fk(:,d,o)
      case (4); v = fl(:,d,o)
      end select
    end function f1pick

    function f2pick(ca, cb, d, o) result(v)
      integer, intent(in) :: ca, cb, d, o
      real(kind=dp) :: v(ngnr)
      integer :: lo, hi
      lo = min(ca,cb); hi = max(ca,cb)
      select case (lo*10+hi)
      case (11); v = f2_11(:,d,o)
      case (12); v = f2_12(:,d,o)
      case (13); v = f2_13(:,d,o)
      case (14); v = f2_14(:,d,o)
      case (22); v = f2_22(:,d,o)
      case (23); v = f2_23(:,d,o)
      case (24); v = f2_24(:,d,o)
      case (33); v = f2_33(:,d,o)
      case (34); v = f2_34(:,d,o)
      case (44); v = f2_44(:,d,o)
      end select
    end function f2pick

  end subroutine compute_der2_ijkl

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

!> @brief Allocate and initialise the 2e SOC integral data structure
!> @details
!>  Extends gdat_init for the 2e mean-field SOC case. Allocates the standard
!>  Rys quadrature work arrays (b00, b01, b10, c00, d00, f00, rw, abv, PQ, PB, QD,
!>  gnm, gnkl, gijkl, dij, dkl) plus fi and fj which hold the derivative-shifted
!>  integrals needed by the SOC recurrence (fi = d/dA phi_i, fj = d/dB phi_j).
!>  Sizes are determined by maxang and gdat%nder.
!>
!> @param[inout] gdat    SOC integral data structure (soc2e_int_data_t)
!> @param[in]    maxang  Maximum angular momentum in the basis
!> @param[in]    dtol    Distance screening threshold
!> @param[in]    dabcut  |AB|^2 cut-off for shell-pair prescreening
!> @param[out]   stat    Allocate status (0 = success)
subroutine soc2e_gdat_init(gdat, maxang, dtol, dabcut, stat)

    implicit none

    class(soc2e_int_data_t), intent(inout) :: gdat
    integer, intent(in) :: maxang
    real(kind=dp), intent(in) :: dtol, dabcut
    integer, intent(out) :: stat
    integer :: mxbra, mxcart, mxrys

    gdat%dtol   = dtol
    gdat%dabcut = dabcut**2

    mxrys  = (4*maxang + 2 + gdat%nder)/2
    mxcart = maxang + 1 + gdat%nder
    mxbra  = 2*mxcart - 1

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
      gdat%rw   (mxrys*2, MAXCONTR), &
      gdat%ai   (MAXCONTR), &
      gdat%aj   (MAXCONTR), &
      gdat%ak   (MAXCONTR), &
      gdat%al   (MAXCONTR), &
      gdat%fi   (mxcart**4 *MAXCONTR*3), &
      gdat%fj   (mxcart**4 *MAXCONTR*3), &
      stat=stat)

  end subroutine soc2e_gdat_init

!> @brief Deallocate all arrays in the 2e SOC integral data structure
!> @param[inout] gdat  SOC integral data structure to be cleaned
subroutine soc2e_gdat_clean(gdat)

    implicit none

    class(soc2e_int_data_t), intent(inout) :: gdat

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

end subroutine soc2e_gdat_clean

!> @brief Store shell quartet (i,j,k,l) metadata in the SOC integral data structure
!> @details
!>  Records angular momenta, AO offsets, and contraction degrees for the four
!>  shells of the current quartet. Called once per quartet before soc2e_rys_compute.
!>
!> @param[inout] gdat   SOC integral data structure
!> @param[in]    basis  Basis set descriptor
!> @param[in]    i,j,k,l  Shell indices of the current quartet (1-based)
subroutine soc2e_gdat_set_ids(gdat, basis, i, j, k, l)

    implicit none

    class(soc2e_int_data_t), intent(inout) :: gdat
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: i, j, k, l

    integer :: id(4), am(4), flips(4), tmp(2)

    flips = [1,2,3,4]
    id = [i,j,k,l]
    am = basis%am([i,j,k,l])
    ! Sort within IJ pair and within KL pair only
    ! Never swap IJ with KL: derivatives are on IJ side only
    if (am(1) > am(2)) flips(1:2) = [2,1]
    if (am(3) > am(4)) flips(3:4) = [4,3]

    gdat%id = id(flips)
    gdat%am = am(flips)
    gdat%at = basis%origin(gdat%id)
    gdat%ao_offset = basis%ao_offset(gdat%id)

  end subroutine soc2e_gdat_set_ids

  !> @brief Set shell parameters for a given quartet in the SOC integral data structure
  !> @details
  !>  Mirrors gdat_set_ids but with a SOC-specific shell ordering constraint:
  !>  shells i,j (bra) may be swapped to put lower-am shell first, and similarly
  !>  for k,l (ket), but the bra-ket pair is never exchanged because the derivative
  !>  recurrence (d/dA phi_i) acts only on the bra side.
  !>
  !> @param[inout] gdat  SOC integral data structure
  subroutine soc2e_set_shells(gdat)

    implicit none

    type(soc2e_int_data_t) :: gdat

    gdat%nbf    = num_cart_bf(gdat%am)
    gdat%nroots = (sum(gdat%am) + 2 + gdat%nder) / 2

    call soc2e_prepare_xyz_ids(gdat)

  end subroutine soc2e_set_shells

  !> @brief Prepare Cartesian index tables for the 2e SOC integral recurrence
  !> @details
  !>  Builds the ijklxyz index array that maps each Cartesian component (nx,ny,nz)
  !>  of each shell function to a flat position in the integral array. Also sets
  !>  nroots (number of Rys roots) from the total angular momentum of the quartet.
  !>  Analogous to prepare_xyz_ids but extended to include the derivative index
  !>  dimension needed for d/dA phi_i and d/dB phi_j.
  !>
  !> @param[inout] gdat  SOC integral data structure
  subroutine soc2e_prepare_xyz_ids(gdat)

    implicit none

    type(soc2e_int_data_t) :: gdat
    integer :: i, nj, nk, nl, njkl, nkl
    integer :: nd

    nd = gdat%nder

    nj   = gdat%am(2) + nd + 1
    nk   = gdat%am(3) + 1        ! no derivative on KL
    nl   = gdat%am(4) + 1        ! no derivative on KL
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

  end subroutine soc2e_prepare_xyz_ids

!> @brief Compute 2e SOC AO integrals for a shell quartet (i,j,k,l) via Rys quadrature
!> @details
!>  Outer loop over primitive pairs on the bra (ij) and ket (kl) sides.
!>  For each pair batch, calls:
!>    soc2e_compute_rys_rw   -- Rys roots and weights
!>    compute_coefficients   -- recursion coefficients (b00, b01, b10, c00, d00, f00)
!>    compute_xyz_p0q0       -- base integrals [p0|q0] in gnm
!>    compute_soc2e_xyz      -- apply derivative recurrence to get Cartesian ints
!>    compute_soc2e_ao       -- contract with density and accumulate into wao
!>  Prescreening is applied via the Schwarz bound gmax * dabmax.
!>
!> @param[inout] gdat   SOC integral data structure (set for the current quartet)
!> @param[in]    ppairs Shell-pair list with pre-screened primitive pairs
!> @param[in]    gmax   Maximum Schwarz estimate over all ket pairs (for screening)
!> @param[in]    den    ROHF density matrix (nbf x nbf)
!> @param[inout] wao    2e SOC AO contribution matrix (3 x nbf x nbf); accumulated
subroutine soc2e_rys_compute(gdat, ppairs, gmax, den, wao)

    use int2_pairs, only: int2_pair_storage
    implicit none

    type(soc2e_int_data_t), intent(inout) :: gdat
    type(int2_pair_storage), intent(in) :: ppairs
    real(kind=dp), intent(in) :: gmax
    real(kind=dp), intent(in) :: den(:,:)
    real(kind=dp), intent(inout) :: wao(:,:,:)

    integer :: ijg, klg, maxgg, mmax, ng
    integer :: nimax, njmax, nkmax, nlmax, nmax
    real(kind=dp) :: aa, ab, aandb1, bb, da, db, test
    real(kind=dp) :: pfac, rho
    real(kind=dp) :: p(3), q(3)
    logical :: last
    integer :: id1, id2, ppid_p, ppid_q, npp_p, npp_q
    integer :: mk, ml, ao_k, ao_l
    real(kind=dp) :: dabmax

    call soc2e_set_shells(gdat)

    ! Density screening bound.  It must cover ALL density blocks contracted
    ! in compute_soc2e_ao: the Coulomb term uses D(K,L), while the exchange
    ! terms use the cross blocks D(I,L), D(J,L), D(I,K), D(J,K).  Screening
    ! on the ket block alone silently drops the exchange contributions
    ! whenever D(K,L) ~ 0 although the cross blocks are finite (e.g. the
    ! s-p blocks of a spherical atom are exactly zero), which under-screens
    ! the SOC (C atom 3P spacing 22.3 instead of 16.4 cm-1).
    dabmax = 0.0_dp
    do mk = 1, gdat%nbf(3)
      ao_k = gdat%ao_offset(3) - 1 + mk
      do ml = 1, gdat%nbf(4)
        dabmax = max(dabmax, abs(den(gdat%ao_offset(4)-1+ml, ao_k)))
      end do
      do ml = 1, gdat%nbf(1)
        dabmax = max(dabmax, abs(den(gdat%ao_offset(1)-1+ml, ao_k)))
      end do
      do ml = 1, gdat%nbf(2)
        dabmax = max(dabmax, abs(den(gdat%ao_offset(2)-1+ml, ao_k)))
      end do
    end do
    do mk = 1, gdat%nbf(4)
      ao_l = gdat%ao_offset(4) - 1 + mk
      do ml = 1, gdat%nbf(1)
        dabmax = max(dabmax, abs(den(gdat%ao_offset(1)-1+ml, ao_l)))
      end do
      do ml = 1, gdat%nbf(2)
        dabmax = max(dabmax, abs(den(gdat%ao_offset(2)-1+ml, ao_l)))
      end do
    end do

    if (dabmax * gmax < 5.0d-11) return

    id1 = maxval(gdat%id(1:2))
    id2 = minval(gdat%id(1:2))
    npp_p = ppairs%ppid(1,id1*(id1-1)/2+id2)
    ppid_p = ppairs%ppid(2,id1*(id1-1)/2+id2)

    id1 = maxval(gdat%id(3:4))
    id2 = minval(gdat%id(3:4))
    npp_q = ppairs%ppid(1,id1*(id1-1)/2+id2)
    ppid_q = ppairs%ppid(2,id1*(id1-1)/2+id2)
    if (npp_p*npp_q == 0) return

    nimax = gdat%am(1) + gdat%nder + 1
    njmax = gdat%am(2) + gdat%nder + 1
    nkmax = gdat%am(3) + 1                ! no derivative on KL
    nlmax = gdat%am(4) + 1                ! no derivative on KL

    nmax = gdat%am(1) + gdat%am(2) + 1 + gdat%nder
    mmax = gdat%am(3) + gdat%am(4) + 1   ! no extra for KL

    maxgg = MAXCONTR / gdat%nroots

    ng = 0

    do klg = 1, npp_q
      db = ppairs%k(ppid_q-1+klg)*ppairs%ginv(ppid_q-1+klg)
      bb = ppairs%g(ppid_q-1+klg)
      q  = ppairs%P(:,ppid_q-1+klg)

      do ijg = 1, npp_p
        da = ppairs%k(ppid_p-1+ijg)*ppairs%ginv(ppid_p-1+ijg)
        aa = ppairs%g(ppid_p-1+ijg)
        p  = ppairs%P(:,ppid_p-1+ijg)

        ab   = aa + bb
        pfac = da*db
        test = pfac*pfac

        if (test < gdat%dtol*ab) cycle
        if (test*dabmax*dabmax < gdat%dabcut*ab) cycle

        aandb1 = 1.0_dp/ab
        rho    = aa*bb*aandb1

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

        gdat%PQ(:,ng) = p - q

        if (nmax>1) gdat%PB(:,ng) = ppairs%PB(:,ppid_p-1+ijg)
        if (mmax>1) gdat%QD(:,ng) = ppairs%PB(:,ppid_q-1+klg)

        gdat%dij(:,ng) = ppairs%PA(:,ppid_p-1+ijg) - ppairs%PB(:,ppid_p-1+ijg)
        gdat%dkl(:,ng) = ppairs%PA(:,ppid_q-1+klg) - ppairs%PB(:,ppid_q-1+klg)

        last = klg==npp_q .and. ijg==npp_p

        if (ng==maxgg .or. last) then
          if (ng==0) return

          call compute_soc2e_ints(gdat, den, wao, &
                  ng, nmax, mmax, nimax, njmax, nkmax, nlmax)

          ng = 0
        end if

      end do
    end do

end subroutine soc2e_rys_compute

!> @brief Evaluate 2e SOC integrals for one primitive batch after Rys setup
!> @details
!>  Called from soc2e_rys_compute after Rys roots/weights and recursion
!>  coefficients are ready in gdat. Evaluates the base integrals [p0|q0]
!>  via compute_xyz_p0q0, then calls compute_soc2e_xyz (derivative recurrence)
!>  and compute_soc2e_ao (density contraction and accumulation into wao).
!>
!> @param[inout] gdat              SOC integral data structure
!> @param[in]    den               ROHF density matrix (nbf x nbf)
!> @param[inout] wao               2e SOC AO matrix (3 x nbf x nbf)
!> @param[in]    ng                Number of primitive pairs in this batch
!> @param[in]    nmax, mmax        Maximum angular momentum orders for the recursion
!> @param[in]    nimax..nlmax      Per-shell Cartesian dimension bounds
subroutine compute_soc2e_ints(gdat, den, wao, &
                  ng, nmax, mmax, nimax, njmax, nkmax, nlmax)

    type(soc2e_int_data_t), intent(inout) :: gdat
    real(kind=dp), intent(in) :: den(:,:)
    real(kind=dp), intent(inout) :: wao(:,:,:)

    integer, intent(in) :: ng, nmax, mmax
    integer, intent(in) :: nimax, njmax, nkmax, nlmax

!   Compute roots and weights for quadrature
    call soc2e_compute_rys_rw(gdat, gdat%rw, ng)

!   Compute coefficients for recursion formulae
    call compute_coefficients(gdat%b00, gdat%b01, gdat%b10, &
                gdat%c00, gdat%d00, gdat%f00, &
                gdat%abv, gdat%pq, gdat%pb, gdat%qd, gdat%rw, nmax, mmax, ng, gdat%nroots)

!   Compute x, y, z integrals (2 centers)
    call compute_xyz_p0q0(gdat%gnm, ng*gdat%nroots, nmax, mmax, &
                gdat%b00, gdat%b01, gdat%b10, gdat%c00, gdat%d00, gdat%f00)

!   Compute x, y, z integrals (4 centers)
    call compute_xyz_ijkl(gdat%gijkl, gdat%gnkl, gdat%gnm, &
                ng, gdat%nroots, nmax, mmax, nimax, njmax, nkmax, nlmax, &
                gdat%dij, gdat%dkl)

!   Compute SOC derivative integrals (IJ side only)
    call compute_soc2e_xyz(gdat, gdat%gijkl, &
                ng, gdat%nroots*3, nimax, njmax, nkmax, nlmax, &
                gdat%ai, gdat%aj, gdat%fi, gdat%fj)

    call compute_soc2e_ao(gdat, ng, gdat%nroots*3, gdat%ijklxyz, &
                gdat%gijkl, gdat%fi, gdat%fj, den, wao)
end subroutine compute_soc2e_ints

!> @brief Apply derivative recurrence to build Cartesian 2e SOC integrals
!> @details
!>  Uses the relation d/dA phi_i(A) = i*phi_{i-1}(A) - 2*alpha*phi_{i+1}(A)
!>  (identical to GAMESS XYZ2E, routines XINTI/XINTJ) to build the fi and fj
!>  arrays from the base integrals g0. These are subsequently contracted in
!>  compute_soc2e_ao to form the mean-field SOC AO matrix element.
!>
!> @param[inout] gdat             SOC integral data structure
!> @param[inout] g                Base Cartesian integral array (in/out)
!> @param[in]    ng, nr3          Batch size and Cartesian xyz dimension
!> @param[in]    nimax..nlmax     Per-shell Cartesian dimension bounds
!> @param[in]    aai, aaj         Exponents for shells i and j
!> @param[out]   fi               d/dA integrals (derivative on shell i side)
!> @param[out]   fj               d/dB integrals (derivative on shell j side)
subroutine compute_soc2e_xyz(gdat, g, &
       ng, nr3, nimax, njmax, nkmax, nlmax, &
       aai, aaj, fi, fj)

    implicit none

    type(soc2e_int_data_t) :: gdat
    integer, intent(in) :: ng, nr3, nimax, njmax, nkmax, nlmax
    real(kind=dp) :: g(ng,nr3,nlmax,nkmax,njmax,*)
    real(kind=dp) :: aai(*), aaj(*)
    real(kind=dp) :: fi(ng,nr3,nlmax,nkmax,njmax,*)
    real(kind=dp) :: fj(ng,nr3,nlmax,nkmax,njmax,*)

    integer :: i, j, n
    integer :: ni, nj

    ni = gdat%am(1) + 1
    nj = gdat%am(2) + 1

!   FI: derivative over I index (electron 1, shell i)
    do n = 1, ng
      fi(n,:,:,:,:,1) = -g(n,:,:,:,:,2)*aai(n)
    end do
    if (ni/=1) then
      do i = 2, ni
        do n = 1, ng
          fi(n,:,:,:,:,i) = g(n,:,:,:,:,i-1)*(i-1) &
                          - g(n,:,:,:,:,i+1)*aai(n)
        end do
      end do
    end if

!   FJ: derivative over J index (electron 1, shell j)
    do i = 1, nimax
      do n = 1, ng
        fj(n,:,:,:,1,i) = -g(n,:,:,:,2,i)*aaj(n)
      end do
    end do
    if (nj/=1) then
      do i = 1, nimax
        do j = 2, nj
          do n = 1, ng
            fj(n,:,:,:,j,i) = g(n,:,:,:,j-1,i)*(j-1) &
                             - g(n,:,:,:,j+1,i)*aaj(n)
          end do
        end do
      end do
    end if

  end subroutine compute_soc2e_xyz

!  subroutine compute_soc2e_ao(gdat, ngnr, ijklxyz, g0, fi, fj, den, wao)
!> @brief Contract 2e SOC integrals with the density matrix and accumulate into wao
!> @details
!>  Implements the mean-field contraction of the two-electron SOC integrals:
!>    W_x(mu,nu) += sum_{lambda,sigma} P(lambda,sigma) * [d_y phi_mu  | r^{-1} | d_z phi_sigma] * D(nu,lambda)
!>                                                      - [d_z phi_mu  | r^{-1} | d_y phi_sigma] * D(nu,lambda)
!>  (and cyclic permutations for Wy, Wz). The derivative integrals fi (on shell i)
!>  and fj (on shell j) are provided by compute_soc2e_xyz. The routine exploits
!>  8-fold permutation symmetry (IJ <-> JI, KL <-> LK, IJ <-> KL) to reduce cost.
!>
!> @param[inout] gdat    SOC integral data structure (shell quartet metadata)
!> @param[in]    ng      Number of primitive pairs in this batch
!> @param[in]    nr3     xyz dimension of the integral arrays (3 for x,y,z)
!> @param[in]    ijklxyz Cartesian index table from soc2e_prepare_xyz_ids
!> @param[in]    g0      Base (unshifted) Cartesian integrals
!> @param[in]    fi      Derivative integrals d/dA on shell i
!> @param[in]    fj      Derivative integrals d/dB on shell j
!> @param[in]    den     ROHF density matrix (nbf x nbf)
!> @param[inout] wao     2e SOC AO matrix (3 x nbf x nbf); Lx,Ly,Lz accumulated
subroutine compute_soc2e_ao(gdat, ng, nr3, ijklxyz, g0, fi, fj, den, wao)
    implicit none
    integer, intent(in) :: ng, nr3
    real(kind=dp), intent(in) :: g0(ng,nr3,*)
    real(kind=dp), intent(in) :: fi(ng,nr3,*)
    real(kind=dp), intent(in) :: fj(ng,nr3,*)

    type(soc2e_int_data_t), intent(inout) :: gdat
    integer, intent(in) :: ijklxyz(:,:,:)
    real(kind=dp), intent(in) :: den(:,:)
    real(kind=dp), intent(inout) :: wao(:,:,:)
    integer :: r, rx, ry, rz
    integer :: i, j, k, l
    integer :: nx, ny, nz
    integer :: ao_I, ao_J, ao_K, ao_L
    integer :: pi, pj, pk, pl
    integer :: jmax
    real(kind=dp) :: sol(3), val2, val3, val4
    logical :: same_KL, same_IJ

    ao_I = gdat%ao_offset(1) - 1
    ao_J = gdat%ao_offset(2) - 1
    ao_K = gdat%ao_offset(3) - 1
    ao_L = gdat%ao_offset(4) - 1
!   Check if KL shells are the same — determines if off-diagonal KL factor applies
    same_KL = (gdat%id(3) == gdat%id(4))
    same_IJ = (gdat%id(1) == gdat%id(2))

    do i = 1, gdat%nbf(1)
      pi = ao_I + i
        if (same_IJ) then
            jmax = i - 1
        else
            jmax = gdat%nbf(2)
        end if
      do j = 1, jmax!gdat%nbf(2)
        pj = ao_J + j
        do k = 1, gdat%nbf(3)
          pk = ao_K + k
          do l = 1, gdat%nbf(4)
            pl = ao_L + l

            nx = ijklxyz(1,i,1)+ijklxyz(1,j,2)+ijklxyz(1,k,3)+ijklxyz(1,l,4) !+ 1
            ny = ijklxyz(2,i,1)+ijklxyz(2,j,2)+ijklxyz(2,k,3)+ijklxyz(2,l,4) !+ 1
            nz = ijklxyz(3,i,1)+ijklxyz(3,j,2)+ijklxyz(3,k,3)+ijklxyz(3,l,4) !+ 1
            sol = 0.0_dp
            do r = 1, nr3/3
              rx = r
              ry = nr3/3 + r
              rz = 2*(nr3/3) + r
              sol(1) = sol(1) + sum((fi(:,ry,ny)*fj(:,rz,nz) - fj(:,ry,ny)*fi(:,rz,nz)) * g0(:,rx,nx))
              sol(2) = sol(2) + sum((fi(:,rz,nz)*fj(:,rx,nx) - fj(:,rz,nz)*fi(:,rx,nx)) * g0(:,ry,ny))
              sol(3) = sol(3) + sum((fi(:,rx,nx)*fj(:,ry,ny) - fj(:,rx,nx)*fi(:,ry,ny)) * g0(:,rz,nz))
            end do
            sol = -sol

!           ---- Coulomb 21 ----
!           W(I,J) += 2*(1+delta_KL)*D(K,L)*VAL
!           W(J,I) -= 2*(1+delta_KL)*D(K,L)*VAL
            val2 = den(pk, pl) * 2.0_dp
            if (.not. same_KL) val2 = val2 * 2.0_dp
            wao(:, pi, pj) = wao(:, pi, pj) + val2 * sol
            wao(:, pj, pi) = wao(:, pj, pi) - val2 * sol

!           ---- Exchange 12 ----
!           W(K,J) += -3*D(I,L)*VAL
!           W(K,I) -= -3*D(J,L)*VAL
            val3 = -3.0_dp
            wao(:, pk, pj) = wao(:, pk, pj) + val3 * den(pi, pl) * sol
            wao(:, pk, pi) = wao(:, pk, pi) - val3 * den(pj, pl) * sol
            if (.not. same_KL) then
!             W(L,J) += -3*D(I,K)*VAL
!             W(L,I) -= -3*D(J,K)*VAL
              wao(:, pl, pj) = wao(:, pl, pj) + val3 * den(pi, pk) * sol
              wao(:, pl, pi) = wao(:, pl, pi) - val3 * den(pj, pk) * sol
            end if

!           ---- Exchange 21 ----
!           W(I,L) += -3*D(K,J)*VAL
!           W(J,L) -= -3*D(K,I)*VAL
            val4 = -3.0_dp
            wao(:, pi, pl) = wao(:, pi, pl) + val4 * den(pk, pj) * sol
            wao(:, pj, pl) = wao(:, pj, pl) - val4 * den(pk, pi) * sol
            if (.not. same_KL) then
!             W(I,K) += -3*D(L,J)*VAL
!             W(J,K) -= -3*D(L,I)*VAL
              wao(:, pi, pk) = wao(:, pi, pk) + val4 * den(pl, pj) * sol
              wao(:, pj, pk) = wao(:, pj, pk) - val4 * den(pl, pi) * sol
            end if

          end do
        end do
      end do
    end do
end subroutine compute_soc2e_ao

  !> @brief Top-level driver for the 2e mean-field SOC correction
  !> @details
  !>  Computes the mean-field two-electron spin-orbit coupling matrix in the AO basis:
  !>    W_x(mu,nu) = sum_{lambda,sigma} P(lambda,sigma) *
  !>                 [<d_y phi_mu | r^{-1} | d_z phi_sigma> - <d_z phi_mu | r^{-1} | d_y phi_sigma>]
  !>  (and cyclic permutations for Wy, Wz). The result is stored in wao(3, nbf, nbf).
  !>
  !>  Procedure:
  !>    1. Build shell-pair list with Schwarz prescreening (cutoff = 1e-10)
  !>    2. Compute Schwarz estimates for screening
  !>    3. Loop over shell quartet (ij, kl); for each quartet call soc2e_rys_compute
  !>
  !>  The physical identity used is:
  !>    <mu | Z * L_x / r^3 | nu>  =  <d_y mu | 1/r | d_z nu> - <d_z mu | 1/r | d_y nu>
  !>
  !> @param[inout] infos  OQP information struct (basis, atoms, MPI info)
  !> @param[in]    basis  Basis set descriptor
  !> @param[in]    den    ROHF density matrix (nbf x nbf)
  !> @param[inout] wao    Output: 2e SOC AO matrix (3 x nbf x nbf), zero-initialised by caller
  subroutine soc2e_driver(infos, basis, den, wao)

    use types, only: information
    use int2_pairs, only: int2_pair_storage, int2_cutoffs_t
    use int2_compute, only: ints_exchange
    use parallel, only: par_env_t
    use constants, only: tol_int

    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: den(:,:)
    real(kind=dp), intent(inout) :: wao(:,:,:)

    type(soc2e_int_data_t) :: gdat
    type(int2_pair_storage) :: ppairs
    type(int2_cutoffs_t) :: cutoffs
    type(par_env_t) :: pe

    real(kind=dp), allocatable :: schwarz_ints(:,:)
    real(kind=dp), allocatable :: den_phys(:,:)


    real(kind=dp) :: cutoff, dabcut, gmax
    real(kind=dp) :: dtol, rtol, zbig
    integer :: i, j, k, l, ij, kl, iok, mpi_ij

    integer iao, jao

    call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    cutoff = 1.0d-10

    zbig = maxval(basis%ex)
    dabcut = 1.0d-11
    if (zbig>1.0d+06) dabcut = dabcut/10
    if (zbig>1.0d+07) dabcut = dabcut/10

    dtol = 10.0d0**(-tol_int)
    rtol = log(10.0_dp)*tol_int

    call cutoffs%set(&
            cutoff_integral_value=dabcut, &
            cutoff_exp=rtol, &
            cutoff_prefactor_pq=dtol, &
            cutoff_prefactor_p=dtol)
    call ppairs%alloc(basis, cutoffs)
    call ppairs%compute(basis, cutoffs)

    allocate(schwarz_ints(basis%nshell, basis%nshell))
    call ints_exchange(basis, schwarz_ints)

    dtol = dtol*dtol

    allocate(den_phys(basis%nbf, basis%nbf))

    do jao = 1, basis%nbf
      do iao = 1, basis%nbf
        den_phys(iao, jao) = den(iao, jao) * basis%bfnrm(iao) * basis%bfnrm(jao)
      end do
    end do

!$omp parallel &
!$omp   firstprivate(gdat) &
!$omp   private(i, j, k, l, ij, kl, gmax, iok, mpi_ij) &
!$omp   reduction(+:wao)

    call gdat%init(basis%mxam, dtol, dabcut, iok)

!$omp barrier
    if (infos%mpiinfo%usempi) mpi_ij = 0

    do i = 1, basis%nshell
      do j = 1, i
        ij = i*(i-1)/2+j
        if (ppairs%ppid(1,ij)==0) cycle
        if (infos%mpiinfo%usempi) then
          mpi_ij = mpi_ij+1
          if (mod(mpi_ij, pe%size) /= pe%rank) cycle
        end if

!$omp do schedule(dynamic)
        do k = 1, basis%nshell
          do l = 1, k
            kl = k*(k-1)/2+l
            if (ppairs%ppid(1,kl)==0) cycle

            gmax = schwarz_ints(i,j)*schwarz_ints(k,l)
            if (gmax < cutoff) cycle
            call gdat%set_ids(basis, i, j, k, l)
            call soc2e_rys_compute(gdat, ppairs, gmax, den_phys, wao)
          end do
        end do
!$omp end do

      end do
    end do
    call gdat%clean()
!$omp end parallel

    call pe%allreduce(wao, size(wao))
    do jao = 1, basis%nbf
      do iao = 1, basis%nbf
        wao(:, iao, jao) = wao(:, iao, jao) * basis%bfnrm(iao) * basis%bfnrm(jao)
      end do
    end do
    call ppairs%clean()
    deallocate(schwarz_ints)
    deallocate(den_phys)
  end subroutine soc2e_driver
end module grd2_rys
