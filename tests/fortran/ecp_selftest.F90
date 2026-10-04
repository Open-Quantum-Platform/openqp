!> Self-test of the native ECP integrals (source/ecp.F90) against
!> references that do not use its expansion machinery.
!>
!> System: one O p primitive at the origin; one s and one p primitive on Br at
!> (-0.4117, 1.7831, 0) Angstrom carrying the LANL2DZ Br ECP (18 terms).
!> Coefficients multiply the unnormalised monomials, as basis%cc does.
!>
!>   err(1)  O py,py: |native - 0.16376847398766514|, an independent
!>           mpmath (type 1) + 2-D Gauss-Legendre/adaptive radial (type 2)
!>           evaluation
!>   err(2)  Br s,s on the ECP centre: closed form
!>             4 pi c^2 sum_{t in L, l=0} d_t Gamma((n_t+1)/2) / (2 q_t^((n_t+1)/2))
!>   err(3)  Br px,px on the ECP centre: closed form
!>             (4 pi/3) c^2 sum_{t in L, l=1} d_t Gamma((n_t+3)/2) / (2 q_t^((n_t+3)/2))
!>   err(4)  max |first derivative - five-point FD of the value|, h = 1e-3 bohr
!>   err(5)  max |second derivative - five-point FD of the first|, h = 1e-3 bohr
!>   err(6)  max |sum over atoms of the first derivatives| (translational invariance)
!>
!> System 2: an f (zeta 1.148) and a g (zeta 2.376) primitive on F at
!> (0.1, 0.2, 1.91) Angstrom, with the def2 iodine ECP (29 terms) at the origin.
!> libecpint 1.0.7 returned -4.0e-5 for the element of err(7).
!>
!>   err(7)  F f_zzz, g_zzzz: |native - (-2.699075746770718e-4)|, independent
!>           mpmath (type 1) + numerical projection/adaptive radial (type 2)
!>   err(8), err(9)  as err(4), err(5) for system 2
!>   err(10) max over both systems of |ecp_hess_contract - explicit contraction
!>           of the stored second-derivative matrices| / max |Hessian|, for a
!>           fixed symmetric density
module ecp_selftest_mod
  use ecp_tool, only: ecp_raw_ints, ecp_hess_start, ecp_hess_contract
  use basis_tools, only: basis_set
  use precision, only: dp
  use, intrinsic :: iso_c_binding, only: c_double
  implicit none

contains

  subroutine oqp_ecp_selftest(err) bind(C, name='oqp_ecp_selftest')
    real(c_double), intent(out) :: err(10)

    real(dp), parameter :: ANG = 1.0_dp/0.529177210903_dp
    real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
    integer, parameter :: NT = 18
    integer, parameter :: ECP_L(NT) = [3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
    integer, parameter :: ECP_N(NT) = [1, 2, 2, 2, 0, 1, 2, 2, 0, 1, 2, 2, 2, 0, 1, 2, 2, 2]
    real(dp), parameter :: ECP_Z(NT) = [213.6143969_dp, 41.058538_dp, 8.708653_dp, 2.6074661_dp, &
        54.1980682_dp, 32.9053558_dp, 13.674489_dp, 3.0341152_dp, &
        54.256334_dp, 26.0095593_dp, 28.2012995_dp, 9.4341061_dp, 2.5321764_dp, &
        87.6328721_dp, 61.7373377_dp, 32.4385104_dp, 8.7537199_dp, 1.6633189_dp]
    real(dp), parameter :: ECP_D(NT) = [-28.0_dp, -134.9268852_dp, -41.9271913_dp, -5.933642_dp, &
        3.0_dp, 27.3430642_dp, 118.8028847_dp, 43.4354876_dp, &
        5.0_dp, 25.0504252_dp, 92.6157463_dp, 95.8249016_dp, 26.2684983_dp, &
        3.0_dp, 22.5533557_dp, 178.1241988_dp, 76.9924162_dp, 9.481827_dp]

    type(basis_set) :: b
    real(dp) :: coord(3, 2), ref, q
    real(dp), allocatable :: v(:), d1(:)
    integer, parameter :: NR = 7
    integer(8), parameter :: NN = NR*NR
    integer :: t, k

    coord(:, 1) = 0.0_dp
    coord(:, 2) = [-0.4117_dp, 1.7831_dp, 0.0_dp]*ANG
    b%nshell = 3
    allocate(b%am(3), b%origin(3), b%ncontr(3), b%g_offset(3), b%ex(3), b%cc(3))
    b%am = [1, 0, 1]
    b%origin = [1, 2, 2]
    b%ncontr = [1, 1, 1]
    b%g_offset = [1, 2, 3]
    b%ex = [0.2137_dp, 0.1905_dp, 0.1377_dp]
    b%cc = [0.20710745912085646_dp, 0.2055092277258375_dp, 0.11956594844924859_dp]
    b%ecp_params%is_ecp = .true.
    allocate(b%ecp_params%n_expo(1), b%ecp_params%ecp_coord(3))
    b%ecp_params%n_expo = NT
    b%ecp_params%ecp_coord = coord(:, 2)
    allocate(b%ecp_params%ecp_am(NT), source=ECP_L)
    allocate(b%ecp_params%ecp_r_ex(NT), source=ECP_N)
    allocate(b%ecp_params%ecp_ex(NT), source=ECP_Z)
    allocate(b%ecp_params%ecp_cc(NT), source=ECP_D)

    call ecp_raw_ints(b, coord, 0, v)
    ! raw order: O px py pz | Br s | Br px py pz
    err(1) = abs(v(1*NR + 2) - 0.16376847398766514_dp)

    ref = 0.0_dp
    do t = 1, NT
      if (ECP_L(t) /= 3 .and. ECP_L(t) /= 0) cycle
      q = 2.0_dp*b%ex(2) + ECP_Z(t)
      ref = ref + ECP_D(t)*gamma(0.5_dp*(ECP_N(t) + 1))/(2.0_dp*q**(0.5_dp*(ECP_N(t) + 1)))
    end do
    err(2) = abs(v(3*NR + 4) - 4.0_dp*PI*b%cc(2)**2*ref)

    ref = 0.0_dp
    do t = 1, NT
      if (ECP_L(t) /= 3 .and. ECP_L(t) /= 1) cycle
      q = 2.0_dp*b%ex(3) + ECP_Z(t)
      ref = ref + ECP_D(t)*gamma(0.5_dp*(ECP_N(t) + 3))/(2.0_dp*q**(0.5_dp*(ECP_N(t) + 3)))
    end do
    err(3) = abs(v(4*NR + 5) - 4.0_dp*PI/3.0_dp*b%cc(3)**2*ref)

    ! finite differences: the ECP moves with atom 2
    call fd_check(b, coord, 2, err(4), err(5))
    call ecp_raw_ints(b, coord, 1, d1)
    err(6) = 0.0_dp
    do k = 1, 3
      err(6) = max(err(6), maxval(abs(d1((k - 1)*NN + 1:k*NN) + d1((k + 2)*NN + 1:(k + 3)*NN))))
    end do

    call contract_check(b, coord, err(10))
    call system2(err(7), err(8), err(9), q)
    err(10) = max(err(10), q)
  end subroutine oqp_ecp_selftest

  subroutine system2(eref, e1, e2, ec)
    real(c_double), intent(out) :: eref, e1, e2
    real(dp), intent(out) :: ec
    real(dp), parameter :: ANG = 1.0_dp/0.529177210903_dp
    integer, parameter :: NT = 29
    integer, parameter :: ECP_L(NT) = [3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, &
                                       2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    real(dp), parameter :: ECP_Z(NT) = [19.458609_dp, 19.34926_dp, 4.823767_dp, 4.884315_dp, &
        40.015835_dp, 17.429747_dp, 9.005484_dp, 19.458609_dp, 19.34926_dp, 4.823767_dp, 4.884315_dp, &
        15.355466_dp, 14.971833_dp, 8.960164_dp, 8.259096_dp, 19.458609_dp, 19.34926_dp, 4.823767_dp, &
        4.884315_dp, 15.068908_dp, 14.555322_dp, 6.718647_dp, 6.456393_dp, 1.191779_dp, 1.291157_dp, &
        19.458609_dp, 19.34926_dp, 4.823767_dp, 4.884315_dp]
    real(dp), parameter :: ECP_D(NT) = [-21.84204_dp, -28.468191_dp, -0.243713_dp, -0.320804_dp, &
        49.994293_dp, 281.025317_dp, 61.573326_dp, 21.84204_dp, 28.468191_dp, 0.243713_dp, 0.320804_dp, &
        67.442841_dp, 134.881137_dp, 14.675051_dp, 29.375666_dp, 21.84204_dp, 28.468191_dp, 0.243713_dp, &
        0.320804_dp, 35.439529_dp, 53.176057_dp, 9.067195_dp, 13.206937_dp, 0.089335_dp, 0.05238_dp, &
        21.84204_dp, 28.468191_dp, 0.243713_dp, 0.320804_dp]
    type(basis_set) :: b
    real(dp) :: coord(3, 2)
    real(dp), allocatable :: v(:)
    integer :: ie

    coord(:, 1) = 0.0_dp
    coord(:, 2) = [0.1_dp, 0.2_dp, 1.91_dp]*ANG
    b%nshell = 2
    allocate(b%am(2), b%origin(2), b%ncontr(2), b%g_offset(2), b%ex(2), b%cc(2))
    b%am = [3, 4]
    b%origin = [2, 2]
    b%ncontr = [1, 1]
    b%g_offset = [1, 2]
    b%ex = [1.148_dp, 2.376_dp]
    b%cc = [7.778024868236892_dp, 123.19916794689166_dp]
    b%ecp_params%is_ecp = .true.
    allocate(b%ecp_params%n_expo(1), b%ecp_params%ecp_coord(3))
    b%ecp_params%n_expo = NT
    b%ecp_params%ecp_coord = coord(:, 1)
    allocate(b%ecp_params%ecp_am(NT), source=ECP_L)
    allocate(b%ecp_params%ecp_r_ex(NT), source=[(2, ie = 1, NT)])
    allocate(b%ecp_params%ecp_ex(NT), source=ECP_Z)
    allocate(b%ecp_params%ecp_cc(NT), source=ECP_D)

    call ecp_raw_ints(b, coord, 0, v)
    ! raw order: f (10 components, zzz last) | g (15 components, zzzz last); nraw = 25
    eref = abs(v((10 - 1)*25 + 25) - (-2.699075746770718e-4_dp))
    call fd_check(b, coord, 1, e1, e2)
    call contract_check(b, coord, ec)
  end subroutine system2

  !> The contracted second derivatives against the stored matrices of
  !> ecp_raw_ints (deriv 2) contracted with the same symmetric density
  subroutine contract_check(b, coord, e)
    type(basis_set), intent(in) :: b
    real(dp), intent(in) :: coord(:,:)
    real(dp), intent(out) :: e
    real(dp), allocatable :: d2(:), dens(:,:), ref(:,:), h(:,:)
    integer :: natm, nraw, i, j, ia, ib, k, kb, gi, gj, mat
    integer(8) :: nn
    real(dp) :: val

    natm = size(coord, 2)
    call ecp_raw_ints(b, coord, 2, d2)
    nraw = sum([(merge(3, merge(1, (b%am(i) + 1)*(b%am(i) + 2)/2, b%am(i) == 0), &
                 b%am(i) == 1), i = 1, b%nshell)])
    nn = int(nraw, 8)*nraw
    allocate(dens(nraw, nraw), ref(3*natm, 3*natm), h(3*natm, 3*natm))
    do j = 1, nraw
      do i = 1, nraw
        dens(i, j) = 1.0_dp/real(i + j - 1, dp) + 0.1_dp*cos(real(i*j, dp))
      end do
    end do
    ref = 0.0_dp
    do ia = 1, natm
      do k = 1, 3
        gi = 3*(ia - 1) + k
        do ib = 1, natm
          do kb = 1, 3
            gj = 3*(ib - 1) + kb
            if (gi > gj) cycle
            if (ia == ib) then
              mat = ecp_hess_start(ia - 1, ia - 1, natm) + 3 &
                  + (k - 1)*3 - ((k - 1)*(k - 2))/2 + (kb - k)
            else
              mat = ecp_hess_start(ia - 1, ib - 1, natm) + (k - 1)*3 + (kb - 1)
            end if
            val = sum(reshape(d2(int(mat, 8)*nn + 1:int(mat + 1, 8)*nn), [nraw, nraw])*dens)
            ref(gi, gj) = val
            ref(gj, gi) = val
          end do
        end do
      end do
    end do
    h = 0.0_dp
    call ecp_hess_contract(b, coord, dens, h)
    e = maxval(abs(h - ref))/max(1.0_dp, maxval(abs(ref)))
  end subroutine contract_check

  !> Five-point finite differences (h = 1e-3 bohr) of the value and of the
  !> first derivative against the first and second derivatives.  The ECP
  !> centre moves with atom iecp.
  subroutine fd_check(b, coord, iecp, e1, e2)
    type(basis_set), intent(inout) :: b
    real(dp), intent(inout) :: coord(:,:)
    integer, intent(in) :: iecp
    real(c_double), intent(out) :: e1, e2
    real(dp), parameter :: H = 1.0e-3_dp
    real(dp), allocatable :: v(:), d1(:), d2(:), fp(:), f2(:), c0(:,:)
    real(dp) :: fdw(4)
    integer :: natm, ia, k, ib, kb, gi, gj, mat, step, isgn, nraw, nmat1
    integer(8) :: nn

    natm = size(coord, 2)
    call ecp_raw_ints(b, coord, 0, v)
    nn = size(v, kind=8)
    nraw = nint(sqrt(real(nn, dp)))
    nmat1 = 3*natm
    call ecp_raw_ints(b, coord, 1, d1)
    call ecp_raw_ints(b, coord, 2, d2)
    c0 = coord
    fdw = [1.0_dp, -8.0_dp, 8.0_dp, -1.0_dp]/(12.0_dp*H)
    e1 = 0.0_dp
    e2 = 0.0_dp
    allocate(fp(nn), f2(nn*nmat1))
    do ia = 1, natm
      do k = 1, 3
        fp = 0.0_dp
        f2 = 0.0_dp
        do step = 1, 4
          isgn = merge(step - 3, step - 2, step <= 2)    ! -2, -1, +1, +2
          coord = c0
          coord(k, ia) = coord(k, ia) + isgn*H
          b%ecp_params%ecp_coord = coord(:, iecp)
          call ecp_raw_ints(b, coord, 0, v)
          fp = fp + fdw(step)*v
          call ecp_raw_ints(b, coord, 1, v)
          f2 = f2 + fdw(step)*v
        end do
        coord = c0
        b%ecp_params%ecp_coord = coord(:, iecp)
        gi = 3*(ia - 1) + k
        e1 = max(e1, maxval(abs(fp - d1((gi - 1)*nn + 1:gi*nn))))
        do ib = 1, natm
          do kb = 1, 3
            gj = 3*(ib - 1) + kb
            if (gi > gj) cycle
            if (ia == ib) then
              mat = ecp_hess_start(ia - 1, ia - 1, natm) + 3 &
                  + (k - 1)*3 - ((k - 1)*(k - 2))/2 + (kb - k)
            else
              mat = ecp_hess_start(ia - 1, ib - 1, natm) + (k - 1)*3 + (kb - 1)
            end if
            e2 = max(e2, maxval(abs(f2((gj - 1)*nn + 1:gj*nn) &
                                    - d2(int(mat, 8)*nn + 1:int(mat + 1, 8)*nn))))
          end do
        end do
      end do
    end do
  end subroutine fd_check

end module ecp_selftest_mod
