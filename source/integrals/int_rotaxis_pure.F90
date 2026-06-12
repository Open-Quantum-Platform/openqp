!> @brief Direct pure-spherical output for the Ishimura rotated-axis ERI engine.
!>
!> @details The Cartesian path (genr22) rotates the rotated-frame block to
!>   the lab frame index-by-index (r30s1d) and, for harmonic-flagged shells,
!>   projects 6d Cartesian components to 5d pure components in a second pass
!>   (genr22_reduce_pure). Both steps are linear maps acting on one shell
!>   index at a time, so they are fused here: each index is transformed once
!>   by T = C * R^T, where R is the per-shell rotated->lab rotation in the
!>   engine's component convention and C the Cartesian->pure projection.
!>   Pure 5d blocks are written directly; the 6d lab-frame Cartesian block
!>   is never materialized, and s-shell (identity) indices are skipped.
!>
!>   Component conventions (must match r30s1d_NN and the projection tables):
!>   d order xx,yy,zz,xy,xz,yz; rotated-frame cross components carry no
!>   sqrt(3) normalization while lab cross components do - the sqrt(3) is
!>   folded into the rotation, exactly as in r30s1d_07.
submodule (int2e_rotaxis) int2e_rotaxis_pure

  ! Sparse pure-d projection C(out, cart) for unit-normalized Cartesians,
  ! cart order xx,yy,zz,xy,xz,yz, pure order m = -2..+2. The values must
  ! stay identical to int2_pure_generated::load_l2 (asserted by
  ! tests/test_ispher_rotaxis_direct_pure.py).
  integer, parameter :: CD_NTERM(6) = [2, 2, 1, 1, 1, 1]
  integer, parameter :: CD_OUT(2,6) = reshape([3, 5,  3, 5,  3, 0,  1, 0,  4, 0,  2, 0], [2,6])
  real(dp), parameter :: CD_COEF(2,6) = reshape([ &
      -4.999999999999999e-01_dp,  8.660254037844386e-01_dp, &
      -4.999999999999999e-01_dp, -8.660254037844386e-01_dp, &
       9.999999999999999e-01_dp,  0.0_dp, &
       1.000000000000000e+00_dp,  0.0_dp, &
       1.000000000000000e+00_dp,  0.0_dp, &
       1.000000000000000e+00_dp,  0.0_dp], [2,6])

contains

  module subroutine genr22_pure(basis, ppairs, grotspd, shell_ids, flips, cutoffs, nbf, emu2)
    use constants, only: NUM_CART_BF
    implicit none

    type(basis_set), intent(in) :: basis
    type(int2_pair_storage), intent(in) :: ppairs
    real(kind=dp), intent(inout) :: grotspd(*)
    integer, intent(in) :: shell_ids(4)
    integer, intent(out) :: flips(4)
    type(int2_cutoffs_t), intent(in) :: cutoffs
    integer, intent(out) :: nbf(4)
    real(kind=dp), optional :: emu2

    real(kind=dp) :: prot(3,3)
    real(kind=dp) :: t(6,6,4)
    integer :: jtype
    integer :: ids(4), am(4), nin(4)
    integer :: s

    call genr22_core(basis, ppairs, grotspd, shell_ids, flips, cutoffs, prot, jtype, emu2)

    ids = shell_ids(flips)
    am = basis%am(ids)

    do s = 1, 4
      nin(s) = NUM_CART_BF(am(s))
      call build_pure_rotation(am(s), basis%harmonic(ids(s)), prot, t(:,:,s), nbf(s))
    end do

    call apply_index_transforms(grotspd, am, nin, nbf, t)

  end subroutine genr22_pure

  !> Build the fused rotated->lab(->pure) transform for one shell index.
  !> tmat(o,nu): coefficient of rotated-frame component nu in output
  !> component o. For l=1 this is the plain rotation (r30s1d_02); for l=2
  !> it is the d rotation of r30s1d_07 composed with the constant pure
  !> projection (or the rotation alone for Cartesian-flagged shells).
  !> l=0 indices are identity and skipped by the caller.
  subroutine build_pure_rotation(l, pure, prot, tmat, nout)
    implicit none

    integer, intent(in) :: l, pure
    real(kind=dp), intent(in) :: prot(3,3)
    real(kind=dp), intent(out) :: tmat(6,6)
    integer, intent(out) :: nout

    real(kind=dp) :: q(6,6)
    integer :: c, k, nu

    select case (l)
    case (0)
      nout = 1

    case (1)
      ! f_lab(i) = sum_nu f_rot(nu) * prot(nu,i), as in r30s1d_02
      nout = 3
      do nu = 1, 3
        tmat(1:3,nu) = prot(nu,1:3)
      end do

    case (2)
      ! q(nu,i): rotated d component nu -> lab d component i (r30s1d_07)
      q(1,1:3) = prot(1,1:3)*prot(1,1:3)
      q(2,1:3) = prot(2,1:3)*prot(2,1:3)
      q(3,1:3) = prot(3,1:3)*prot(3,1:3)
      q(4,1:3) = prot(1,1:3)*prot(2,1:3)*2.0_dp
      q(5,1:3) = prot(1,1:3)*prot(3,1:3)*2.0_dp
      q(6,1:3) = prot(2,1:3)*prot(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (prot(1,1)*prot(1,2:3))
      q(2,4:5) = sqrt3 * (prot(2,1)*prot(2,2:3))
      q(3,4:5) = sqrt3 * (prot(3,1)*prot(3,2:3))
      q(4,4:5) = sqrt3 * (prot(1,1)*prot(2,2:3)+prot(2,1)*prot(1,2:3))
      q(5,4:5) = sqrt3 * (prot(1,1)*prot(3,2:3)+prot(3,1)*prot(1,2:3))
      q(6,4:5) = sqrt3 * (prot(2,1)*prot(3,2:3)+prot(3,1)*prot(2,2:3))

      q(1,6) = sqrt3 * (prot(1,2)*prot(1,3))
      q(2,6) = sqrt3 * (prot(2,2)*prot(2,3))
      q(3,6) = sqrt3 * (prot(3,2)*prot(3,3))
      q(4,6) = sqrt3 * (prot(1,2)*prot(2,3)+prot(2,2)*prot(1,3))
      q(5,6) = sqrt3 * (prot(1,2)*prot(3,3)+prot(3,2)*prot(1,3))
      q(6,6) = sqrt3 * (prot(2,2)*prot(3,3)+prot(3,2)*prot(2,3))

      if (pure == 1) then
        ! tmat = C * Q^T from the constant sparse d projection
        nout = 5
        tmat(1:5,1:6) = 0.0_dp
        do c = 1, 6
          do k = 1, CD_NTERM(c)
            tmat(CD_OUT(k,c),1:6) = tmat(CD_OUT(k,c),1:6) + CD_COEF(k,c) * q(1:6,c)
          end do
        end do
      else
        nout = 6
        do c = 1, 6
          tmat(c,1:6) = q(1:6,c)
        end do
      end if

    case default
      error stop 'genr22_pure: rotated-axis engine supports l <= 2 only'
    end select

  end subroutine build_pure_rotation

  !> Apply the per-index transforms to the quartet block. Storage order is
  !> B(n4,n3,n2,n1): storage dimension k holds canonical shell slot 5-k.
  !> Sequential one-index contractions; l=0 (identity) indices are skipped,
  !> the first active stage reads from f and the last writes back to f, so
  !> no boundary copies are made unless only one index is active.
  subroutine apply_index_transforms(f, am, nin, nout, t)
    implicit none

    real(kind=dp), intent(inout) :: f(*)
    integer, intent(in) :: am(4), nin(4), nout(4)
    real(kind=dp), intent(in) :: t(6,6,4)

    real(kind=dp) :: work(1296,2)
    integer :: dims(4)
    integer :: k, pos, j, ntot
    integer :: nleft, nright
    integer :: src_id, dst_id, last_k

    dims = nin([4,3,2,1])

    last_k = 0
    do k = 1, 4
      if (am(5-k) > 0) last_k = k
    end do
    if (last_k == 0) return

    src_id = 0                     ! 0 = f, 1/2 = work columns
    do k = 1, 4
      pos = 5 - k
      if (am(pos) == 0) cycle      ! identity index

      nleft = 1
      do j = 1, k-1
        nleft = nleft*dims(j)
      end do
      nright = 1
      do j = k+1, 4
        nright = nright*dims(j)
      end do

      if (k == last_k .and. src_id /= 0) then
        dst_id = 0
      else
        dst_id = merge(2, 1, src_id == 1)
      end if

      if (src_id == 0) then
        call transform_one_dim(f, work(:,dst_id), nleft, nin(pos), nout(pos), nright, t(:,:,pos))
      else if (dst_id == 0) then
        call transform_one_dim(work(:,src_id), f, nleft, nin(pos), nout(pos), nright, t(:,:,pos))
      else
        call transform_one_dim(work(:,src_id), work(:,dst_id), nleft, nin(pos), nout(pos), nright, t(:,:,pos))
      end if

      dims(k) = nout(pos)
      src_id = dst_id
    end do

    if (src_id /= 0) then          ! single active index: copy back
      ntot = product(dims)
      f(1:ntot) = work(1:ntot,src_id)
    end if

  end subroutine apply_index_transforms

  !> dst(:,o,b) = sum_i tmat(o,i) * src(:,i,b)
  subroutine transform_one_dim(src, dst, nleft, ni, no, nright, tmat)
    implicit none

    integer, intent(in) :: nleft, ni, no, nright
    real(kind=dp), intent(in) :: src(nleft, ni, nright)
    real(kind=dp), intent(out) :: dst(nleft, no, nright)
    real(kind=dp), intent(in) :: tmat(6,6)

    integer :: b, o, i

    do b = 1, nright
      do o = 1, no
        dst(:,o,b) = tmat(o,1) * src(:,1,b)
        do i = 2, ni
          dst(:,o,b) = dst(:,o,b) + tmat(o,i) * src(:,i,b)
        end do
      end do
    end do

  end subroutine transform_one_dim

end submodule int2e_rotaxis_pure
