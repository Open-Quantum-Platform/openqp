!> @brief Direct pure-spherical output for the Ishimura rotated-axis ERI engine.
!>
!> @details The Cartesian path (genr22) rotates the rotated-frame block to
!>   the lab frame index-by-index (r30s1d) and, for harmonic-flagged shells,
!>   projects 6d Cartesian components to 5d pure components in a second pass
!>   (genr22_reduce_pure). Both steps are linear maps acting on one shell
!>   index at a time, so they are fused here: each index is transformed once
!>   by T = C * R^T, where R is the per-shell rotated->lab rotation in the
!>   engine's component convention and C the Cartesian->pure projection from
!>   int2_pure_generated. Pure 5d blocks are written directly; the 6d
!>   lab-frame Cartesian block is never materialized.
!>
!>   Component conventions (must match r30s1d_NN and the projection tables):
!>   d order xx,yy,zz,xy,xz,yz; rotated-frame cross components carry no
!>   sqrt(3) normalization while lab cross components do - the sqrt(3) is
!>   folded into the rotation, exactly as in r30s1d_07.
submodule (int2e_rotaxis) int2e_rotaxis_pure

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
    integer :: ids(4), am(4), pure(4), nin(4)
    integer :: s

    call genr22_core(basis, ppairs, grotspd, shell_ids, flips, cutoffs, prot, jtype, emu2)

    ids = shell_ids(flips)
    am = basis%am(ids)
    pure = basis%harmonic(ids)

    do s = 1, 4
      nin(s) = NUM_CART_BF(am(s))
      call build_pure_rotation(am(s), pure(s), prot, t(:,:,s), nbf(s))
    end do

    call apply_index_transforms(grotspd, nin, nbf, t)

  end subroutine genr22_pure

  !> Build the fused rotated->lab(->pure) transform for one shell index.
  !> tmat(o,nu): coefficient of rotated-frame component nu in output
  !> component o. For l=0/1 this is the plain rotation; for l=2 it is the
  !> d rotation of r30s1d_07 composed with the sparse pure projection
  !> (identity for Cartesian-flagged shells).
  subroutine build_pure_rotation(l, pure, prot, tmat, nout)
    use int2_pure_generated, only: int2_shell_projection_t, int2_init_shell_projection
    implicit none

    integer, intent(in) :: l, pure
    real(kind=dp), intent(in) :: prot(3,3)
    real(kind=dp), intent(out) :: tmat(6,6)
    integer, intent(out) :: nout

    type(int2_shell_projection_t) :: proj
    real(kind=dp) :: q(6,6)
    integer :: c, k, o, nu

    tmat = 0.0_dp

    select case (l)
    case (0)
      nout = 1
      tmat(1,1) = 1.0_dp

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

      ! Compose with the pure projection: tmat = C * Q^T. The identity
      ! projection of Cartesian-flagged shells reduces this to Q^T.
      call int2_init_shell_projection(l, pure, proj)
      nout = proj%nout
      do c = 1, proj%ncart
        do k = 1, proj%nterm(c)
          o = proj%out_idx(k,c)
          tmat(o,1:6) = tmat(o,1:6) + proj%coeff(k,c) * q(1:6,c)
        end do
      end do

    case default
      error stop 'genr22_pure: rotated-axis engine supports l <= 2 only'
    end select

  end subroutine build_pure_rotation

  !> Apply the four per-index transforms to the quartet block in place.
  !> Storage order is B(n4,n3,n2,n1): storage dimension k holds canonical
  !> shell slot 5-k. Sequential one-index contractions with ping-pong work
  !> buffers; intermediate dimensions shrink as pure indices are produced.
  subroutine apply_index_transforms(f, nin, nout, t)
    implicit none

    real(kind=dp), intent(inout) :: f(*)
    integer, intent(in) :: nin(4), nout(4)
    real(kind=dp), intent(in) :: t(6,6,4)

    real(kind=dp) :: work(1296,2)
    integer :: dims(4)
    integer :: k, pos, cur, nxt, ntot
    integer :: nleft, nright, j

    dims = nin([4,3,2,1])
    ntot = product(dims)
    work(1:ntot,1) = f(1:ntot)
    cur = 1

    do k = 1, 4
      pos = 5 - k
      nleft = 1
      do j = 1, k-1
        nleft = nleft*dims(j)
      end do
      nright = 1
      do j = k+1, 4
        nright = nright*dims(j)
      end do
      nxt = 3 - cur
      call transform_one_dim(work(:,cur), work(:,nxt), nleft, nin(pos), nout(pos), nright, t(:,:,pos))
      dims(k) = nout(pos)
      cur = nxt
    end do

    ntot = product(dims)
    f(1:ntot) = work(1:ntot,cur)

  end subroutine apply_index_transforms

  !> dst(a,o,b) = sum_i tmat(o,i) * src(a,i,b)
  subroutine transform_one_dim(src, dst, nleft, ni, no, nright, tmat)
    implicit none

    real(kind=dp), intent(in) :: src(*)
    real(kind=dp), intent(out) :: dst(*)
    integer, intent(in) :: nleft, ni, no, nright
    real(kind=dp), intent(in) :: tmat(6,6)

    integer :: a, b, i, o
    integer :: src0, dst0
    real(kind=dp) :: acc

    do b = 1, nright
      do o = 1, no
        dst0 = (b-1)*nleft*no + (o-1)*nleft
        do a = 1, nleft
          acc = 0.0_dp
          src0 = (b-1)*nleft*ni + a
          do i = 1, ni
            acc = acc + tmat(o,i)*src(src0 + (i-1)*nleft)
          end do
          dst(dst0 + a) = acc
        end do
      end do
    end do

  end subroutine transform_one_dim

end submodule int2e_rotaxis_pure
