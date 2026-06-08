!> @brief   Cartesian -> pure spherical-harmonic (c2s) transforms for shells.
!>
!> @details Provides the per-shell transform matrices B(l) that map the
!>          unit-normalized Cartesian components of a shell (OpenQP's
!>          canonical bf_names order; see constants::CART_X/Y/Z) onto the
!>          2l+1 real solid harmonics in CCA/libint order (m = -l..+l).
!>
!>          Convention is identical to the validated Python reference in
!>          pyoqp/oqp/library/symmetry.py (_solid_harmonic_coefficients):
!>          each column of C2S_x holds the Cartesian coefficients of one
!>          spherical component, orthonormal against the intra-shell metric
!>          S of unit-normalized Cartesian Gaussians (B S B^T = I). The
!>          matrices below were generated from that reference and are
!>          re-verified at runtime by c2s_selftest().
!>
!>          The transform is applied to integrals that are already in the
!>          unit-normalized Cartesian basis (e.g. 2e blocks AFTER the
!>          rotation/Rys/libint normalization in int2::shellquartet, where
!>          all backends agree). s and p shells are passed through unchanged
!>          (Cartesian == spherical up to the trivial 1:1 / 3:3 mapping).
module cart2sph

  use precision, only: dp
  use constants, only: NUM_CART_BF, NUM_SPH_BF, BAS_MXANG

  implicit none
  private

  public :: c2s_ncomp
  public :: cart2sph_eri
  public :: cart2sph_mat
  public :: c2s_selftest

  ! l=2 (D): Cart(6) -> Sph(5); column i = Cartesian coeffs of spherical i (m=-l..+l)
  real(dp), parameter :: C2S_D(6,5) = reshape([ &
        0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  1.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  &
        0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  1.000000000000000e+00_dp,  &
       -4.999999999999999e-01_dp, -4.999999999999999e-01_dp,  9.999999999999999e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  &
        0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  1.000000000000000e+00_dp,  0.000000000000000e+00_dp,  &
        8.660254037844386e-01_dp, -8.660254037844386e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp &
     ], shape=[6,5])

  ! l=3 (F): Cart(10) -> Sph(7); column i = Cartesian coeffs of spherical i (m=-l..+l)
  real(dp), parameter :: C2S_F(10,7) = reshape([ &
        0.000000000000000e+00_dp, -7.905694150420950e-01_dp,  0.000000000000000e+00_dp,  1.060660171779821e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  &
        0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  1.000000000000000e+00_dp,  &
        0.000000000000000e+00_dp, -6.123724356957946e-01_dp,  0.000000000000000e+00_dp, -2.738612787525830e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  1.095445115010332e+00_dp,  0.000000000000000e+00_dp,  &
        0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  1.000000000000000e+00_dp,  0.000000000000000e+00_dp, -6.708203932499369e-01_dp,  0.000000000000000e+00_dp, -6.708203932499369e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  &
       -6.123724356957946e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp, -2.738612787525830e-01_dp,  0.000000000000000e+00_dp,  1.095445115010332e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  &
        0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  8.660254037844385e-01_dp,  0.000000000000000e+00_dp, -8.660254037844385e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  &
        7.905694150420950e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp, -1.060660171779821e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp &
     ], shape=[10,7])

  ! l=4 (G): Cart(15) -> Sph(9); column i = Cartesian coeffs of spherical i (m=-l..+l)
  real(dp), parameter :: C2S_G(15,9) = reshape([ &
        0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  1.118033988749895e+00_dp,  0.000000000000000e+00_dp, -1.118033988749895e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  &
        0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp, -7.905694150420950e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  1.060660171779821e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  &
        0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp, -4.225771273642583e-01_dp,  0.000000000000000e+00_dp, -4.225771273642583e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  1.133893419027681e+00_dp,  &
        0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp, -8.964214570007953e-01_dp,  0.000000000000000e+00_dp,  1.195228609334394e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp, -4.008918628686365e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  &
        3.749999999999999e-01_dp,  3.749999999999999e-01_dp,  9.999999999999998e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  2.195775164134199e-01_dp, -8.783100656536798e-01_dp, -8.783100656536798e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  &
        0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp, -8.964214570007953e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  1.195228609334394e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp, -4.008918628686365e-01_dp,  0.000000000000000e+00_dp,  &
       -5.590169943749475e-01_dp,  5.590169943749475e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  9.819805060619657e-01_dp, -9.819805060619657e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  &
        0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  7.905694150420950e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp, -1.060660171779821e+00_dp,  0.000000000000000e+00_dp,  &
        7.395099728874520e-01_dp,  7.395099728874520e-01_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp, -1.299038105676658e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp,  0.000000000000000e+00_dp &
     ], shape=[15,9])

contains

  !> @brief AO component count for a shell, honoring its harmonic flag.
  elemental integer function c2s_ncomp(l, pure) result(n)
    integer, intent(in) :: l, pure
    if (pure == 1 .and. l >= 2) then
      n = NUM_SPH_BF(l)
    else
      n = NUM_CART_BF(l)
    end if
  end function c2s_ncomp

  !> @brief Return the c2s matrix B(l) (ncart x nsph) for a pure shell.
  !> @details Only l = 2,3,4 carry a non-trivial transform. Caller guarantees
  !>          l >= 2 (s/p never reach here because c2s_ncomp keeps them
  !>          Cartesian). The result is a copy sized (NUM_CART_BF(l), 2l+1).
  subroutine c2s_get(l, b)
    integer, intent(in) :: l
    real(dp), allocatable, intent(out) :: b(:,:)
    select case (l)
    case (2)
      b = C2S_D
    case (3)
      b = C2S_F
    case (4)
      b = C2S_G
    case default
      ! Should never happen; leave unallocated so callers fail loudly.
      allocate(b(0,0))
    end select
  end subroutine c2s_get

  !> @brief Contract one index of a 3-way-folded block: out(il,is,ir) =
  !>        sum_ic B(ic,is) * a(il,ic,ir), with the index laid out as
  !>        (left, n_cart, right) in column-major order.
  subroutine contract_index(a, left, ncart, right, b, nsph, out)
    integer, intent(in) :: left, ncart, right, nsph
    real(dp), intent(in) :: a(left, ncart, right)
    real(dp), intent(in) :: b(ncart, nsph)
    real(dp), intent(out) :: out(left, nsph, right)
    integer :: il, is, ic, ir
    real(dp) :: bval
    out = 0.0_dp
    do ir = 1, right
      do is = 1, nsph
        do ic = 1, ncart
          bval = b(ic, is)
          if (bval == 0.0_dp) cycle
          do il = 1, left
            out(il, is, ir) = out(il, is, ir) + bval * a(il, ic, ir)
          end do
        end do
      end do
    end do
  end subroutine contract_index

  !> @brief Transform a 2e shell-quartet block from unit-normalized Cartesian
  !>        to pure spherical for any index whose shell is flagged harmonic.
  !>
  !> @param[inout] ints  flat ERI buffer; on entry holds the Cartesian block
  !>                     with storage dims [nbf(1)..nbf(4)] (column-major,
  !>                     index 1 fastest); on exit the spherical block with
  !>                     dims [nbf_out(1)..nbf_out(4)].
  !> @param[in]    am    angular momentum of the four shells, in storage order
  !> @param[in]    pure  per-shell harmonic flag (1=spherical), storage order
  !> @param[in]    nbf   Cartesian component counts, storage order
  !> @param[out]   nbf_out spherical component counts, storage order
  !>
  !> Storage order means dimension k of `ints` corresponds to am(k)/pure(k).
  !> Callers in int2 pass these already in the flipped (stored) order.
  subroutine cart2sph_eri(ints, am, pure, nbf, nbf_out)
    real(dp), intent(inout) :: ints(:)
    integer, intent(in) :: am(4), pure(4), nbf(4)
    integer, intent(out) :: nbf_out(4)

    integer :: dims(4)        ! running per-index sizes (Cartesian -> spherical)
    integer :: p, left, right, k, nc, ns
    real(dp), allocatable :: b(:,:), src(:), dst(:)

    nbf_out = nbf
    do p = 1, 4
      if (pure(p) /= 1 .or. am(p) < 2) cycle
      nbf_out(p) = NUM_SPH_BF(am(p))
    end do

    if (all(nbf_out == nbf)) return   ! nothing pure -> leave Cartesian block

    dims = nbf
    allocate(src(product(nbf)))
    src(1:product(nbf)) = ints(1:product(nbf))

    do p = 1, 4
      if (pure(p) /= 1 .or. am(p) < 2) cycle
      nc = dims(p)
      ns = NUM_SPH_BF(am(p))
      left  = product(dims(1:p-1))
      right = product(dims(p+1:4))
      call c2s_get(am(p), b)
      allocate(dst(left * ns * right))
      call contract_index(src, left, nc, right, b, ns, dst)
      call move_alloc(dst, src)
      dims(p) = ns
      deallocate(b)
    end do

    k = product(dims)
    ints(1:k) = src(1:k)
    deallocate(src)
  end subroutine cart2sph_eri

  !> @brief Transform a 1e shell-pair block (unit-normalized Cartesian) to
  !>        pure spherical for any harmonic-flagged shell.
  !>
  !> @details The block is laid out with the "fast" shell varying quickest,
  !>          i.e. blk(nn), nn over (slow outer, fast inner) -- the order
  !>          consumed by update_triang_matrix/update_rectangular_matrix
  !>          (fast = shj, slow = shi). On exit blk(1:n_out) holds the
  !>          spherical block in the same fast/slow layout and n_fast_out/
  !>          n_slow_out give its extents.
  subroutine cart2sph_mat(blk, l_fast, pure_fast, l_slow, pure_slow, n_fast_out, n_slow_out)
    use constants, only: shells_pnrm2
    real(dp), intent(inout) :: blk(:)
    integer, intent(in) :: l_fast, pure_fast, l_slow, pure_slow
    integer, intent(out), optional :: n_fast_out, n_slow_out

    integer :: ncf, ncs, nsf, nss, k, ic, ir
    real(dp), allocatable :: b(:,:), src(:), dst(:)

    ncf = NUM_CART_BF(l_fast)
    ncs = NUM_CART_BF(l_slow)
    nsf = c2s_ncomp(l_fast, pure_fast)
    nss = c2s_ncomp(l_slow, pure_slow)
    if (present(n_fast_out)) n_fast_out = nsf
    if (present(n_slow_out)) n_slow_out = nss

    if (nsf == ncf .and. nss == ncs) return   ! nothing pure -> Cartesian block

    allocate(src(ncf*ncs))
    src(1:ncf*ncs) = blk(1:ncf*ncs)

    ! 1e blocks arrive in the pure-power Cartesian normalization (the
    ! shells_pnrm2 per-component factors are applied later by bas_norm_matrix
    ! for Cartesian shells). B is defined for unit-normalized Cartesians, so
    ! for each index we transform, first fold in shells_pnrm2 along that index;
    ! the resulting spherical components are unit-normalized (set_bfnorms then
    ! uses bfnrm = 1 for them). Non-pure indices are left pure-power untouched.

    ! Contract the fast index (storage layout (ncf, ncs), fast contiguous).
    if (pure_fast == 1 .and. l_fast >= 2) then
      do ir = 1, ncs
        do ic = 1, ncf
          src((ir-1)*ncf + ic) = src((ir-1)*ncf + ic) * shells_pnrm2(ic, l_fast)
        end do
      end do
      call c2s_get(l_fast, b)
      allocate(dst(nsf*ncs))
      call contract_index(src, 1, ncf, ncs, b, nsf, dst)
      call move_alloc(dst, src)
      deallocate(b)
    end if
    ! Contract the slow index ((nsf, ncs, 1); the slow index is the outer one).
    if (pure_slow == 1 .and. l_slow >= 2) then
      do ir = 1, ncs
        do ic = 1, nsf
          src((ir-1)*nsf + ic) = src((ir-1)*nsf + ic) * shells_pnrm2(ir, l_slow)
        end do
      end do
      call c2s_get(l_slow, b)
      allocate(dst(nsf*nss))
      call contract_index(src, nsf, ncs, 1, b, nss, dst)
      call move_alloc(dst, src)
      deallocate(b)
    end if

    k = nsf*nss
    blk(1:k) = src(1:k)
    deallocate(src)
  end subroutine cart2sph_mat

  !> @brief Self-test: rebuild the intra-shell metric S of unit-normalized
  !>        Cartesian Gaussians from the canonical exponents and verify
  !>        B(l) S B(l)^T = I for l = 2,3,4. Returns the worst deviation.
  subroutine c2s_selftest(max_err)
    use constants, only: CART_X, CART_Y, CART_Z
    real(dp), intent(out) :: max_err
    integer :: l, nc, ns, i, j, is, js
    real(dp), allocatable :: b(:,:), s(:,:), g(:,:)
    real(dp) :: err

    max_err = 0.0_dp
    do l = 2, 4
      nc = NUM_CART_BF(l)
      ns = NUM_SPH_BF(l)
      call c2s_get(l, b)
      allocate(s(nc, nc))
      do i = 1, nc
        do j = 1, nc
          s(i,j) = cart_overlap(CART_X(i,l), CART_Y(i,l), CART_Z(i,l), &
                                CART_X(j,l), CART_Y(j,l), CART_Z(j,l))
        end do
      end do
      ! g = B^T S B  (ns x ns), should be identity
      allocate(g(ns, ns))
      g = matmul(matmul(transpose(b), s), b)
      do is = 1, ns
        do js = 1, ns
          err = abs(g(is,js) - merge(1.0_dp, 0.0_dp, is == js))
          if (err > max_err) max_err = err
        end do
      end do
      deallocate(b, s, g)
    end do
  end subroutine c2s_selftest

  !> @brief Overlap of two unit-normalized Cartesian Gaussians of the same
  !>        shell (same exponent), i.e. the intra-shell metric element.
  !>        For a normalized x^a y^b z^c: <i|j> = prod_k (a_k+b_k-1)!! /
  !>        sqrt((2a_k-1)!! (2b_k-1)!!) over k in {x,y,z}; zero if any sum odd.
  pure real(dp) function cart_overlap(ax, ay, az, bx, by, bz) result(s)
    integer, intent(in) :: ax, ay, az, bx, by, bz
    if (mod(ax+bx,2) /= 0 .or. mod(ay+by,2) /= 0 .or. mod(az+bz,2) /= 0) then
      s = 0.0_dp
      return
    end if
    s = ratio(ax, bx) * ratio(ay, by) * ratio(az, bz)
  end function cart_overlap

  !> @brief (a+b-1)!! / sqrt((2a-1)!! (2b-1)!!) for one Cartesian axis.
  pure real(dp) function ratio(a, b) result(r)
    integer, intent(in) :: a, b
    r = real(idfact(a+b-1), dp) / sqrt(real(idfact(2*a-1), dp) * real(idfact(2*b-1), dp))
  end function ratio

  !> @brief Integer double factorial n!! with (-1)!! = 0!! = 1.
  pure integer function idfact(n) result(r)
    integer, intent(in) :: n
    integer :: k
    r = 1
    k = n
    do while (k > 1)
      r = r * k
      k = k - 2
    end do
  end function idfact

end module cart2sph
