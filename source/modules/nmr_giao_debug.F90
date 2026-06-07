module nmr_giao_debug_mod

  implicit none
  private
  public giao_h10_twoe_matrix

contains


!> @brief Native RHF GIAO two-electron magnetic-field-derivative Fock matrices.
!> @details Returns the real coefficients of the imaginary first-order (wrt the
!>  uniform field B) two-electron Fock contributions for x/y/z, contracted with
!>  the ground-state density dm: the Coulomb-like image vj, the exchange-like
!>  image vk, and the RHF combination h10 = vj - 0.5*vk.  Each block is
!>  antisymmetric in its AO indices.  Validated against an independent
!>  two-electron GIAO-derivative reference.  This is the production-reusable core
!>  of nmr_giao_h10_twoe_debug; it computes integrals only and emits nothing.
  subroutine giao_h10_twoe_matrix(basis, infos, dm, vj, vk, h10)

    use basis_tools, only: basis_set
    use constants, only: cart_x, cart_y, cart_z, num_cart_bf
    use int2_compute, only: int2_compute_t
    use int2e_rys, only: int2_rys_compute_ordered_am, int2_rys_data_t
    use messages, only: show_message, with_abort
    use precision, only: dp
    use types, only: information

    implicit none

    type(basis_set), intent(in) :: basis
    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: dm(:,:)
    real(kind=dp), intent(out) :: vj(:,:,:), vk(:,:,:), h10(:,:,:)

    type(int2_compute_t), target :: int2_driver
    type(int2_rys_data_t) :: gdat
    real(kind=dp), allocatable :: vk_pre(:,:,:)
    real(kind=dp), allocatable, target :: eri0(:), erir(:)
    real(kind=dp), pointer :: p0(:,:,:,:), pr(:,:,:,:)
    real(kind=dp) :: g(3), d(3), rij(3), norm4, ket_fac, base_val
    integer :: i, j, m, nbf, nshell, maxang, maxcart
    integer :: si, sj, sk, sl, ni, nj, nk, nl, mu, nu, kap, lam
    integer :: ii, jj, kk, ll, axis, mapr, ids(4), am0(4), amr(4), ok
    logical :: zero_shq

    nbf = basis%nbf
    nshell = basis%nshell
    maxang = maxval(basis%am)
    if (maxang + 1 > 6) call show_message('GIAO two-electron: angular momentum too high', with_abort)
    maxcart = num_cart_bf(maxang + 1)

    vj = 0.0_dp
    vk = 0.0_dp
    h10 = 0.0_dp
    allocate(vk_pre(3,nbf,nbf), source=0.0_dp)
    allocate(eri0(maxcart**4), erir(maxcart**4), source=0.0_dp)

    call int2_driver%init(basis, infos)
    ! The magnetic derivative is ordered/antisymmetric in the bra pair.  Rebuild
    ! pair data without angular-momentum reordering so the Rys recurrence sees
    ! PA/PB vectors in the same shell order as the explicit ids below.
    call int2_driver%ppairs%compute(basis, int2_driver%cutoffs, noswap=.true.)
    call gdat%init(maxang + 1, int2_driver%cutoffs, ok)
    if (ok /= 0) call show_message('GIAO two-electron: cannot allocate Rys workspace', with_abort)

    do si = 1, nshell
      ni = basis%naos(si)
      do sj = 1, si
        nj = basis%naos(sj)
        rij = basis%shell_centers(si,1:3) - basis%shell_centers(sj,1:3)
        if (maxval(abs(rij)) < 1.0d-14) cycle
        do sk = 1, nshell
          nk = basis%naos(sk)
          do sl = 1, sk
            nl = basis%naos(sl)
            ket_fac = merge(1.0_dp, 2.0_dp, sk == sl)
            ids = [si, sj, sk, sl]
            am0 = basis%am(ids)
            eri0 = 0.0_dp
            call int2_rys_compute_ordered_am(eri0, gdat, int2_driver%ppairs, ids, am0, zero_shq)
            if (zero_shq) cycle
            p0(1:nl,1:nk,1:nj,1:ni) => eri0(1:nl*nk*nj*ni)
            ! The raised-bra quartet depends only on (ids, am0 + 1 on the bra);
            ! it is identical for all three axes and all (ii,jj,kk,ll), so
            ! compute it ONCE per shell quartet.  (It was previously recomputed
            ! 3*ni*nj*nk*nl times inside the loops below, which made d/f-basis
            ! GIAO NMR intractable.)  Only the Cartesian component index mapr
            ! depends on (ii, axis).
            amr = am0
            amr(1) = amr(1) + 1
            erir = 0.0_dp
            call int2_rys_compute_ordered_am(erir, gdat, int2_driver%ppairs, ids, amr, zero_shq)
            pr(1:nl,1:nk,1:nj,1:num_cart_bf(amr(1))) => erir(1:nl*nk*nj*num_cart_bf(amr(1)))
            do ii = 1, ni
              mu = basis%ao_offset(si) + ii - 1
              do jj = 1, nj
                nu = basis%ao_offset(sj) + jj - 1
                do kk = 1, nk
                  kap = basis%ao_offset(sk) + kk - 1
                  do ll = 1, nl
                    lam = basis%ao_offset(sl) + ll - 1
                    norm4 = basis%bfnrm(mu)*basis%bfnrm(nu)*basis%bfnrm(kap)*basis%bfnrm(lam)
                    base_val = p0(ll,kk,jj,ii)
                    do axis = 1, 3
                      mapr = raised_cart_index(ii, am0(1), axis)
                      d(axis) = (pr(ll,kk,jj,mapr) + basis%shell_centers(si,axis)*base_val) * norm4
                    end do
                    g(1) = -0.5_dp*(rij(2)*d(3) - rij(3)*d(2))
                    g(2) = -0.5_dp*(rij(3)*d(1) - rij(1)*d(3))
                    g(3) = -0.5_dp*(rij(1)*d(2) - rij(2)*d(1))
                    do m = 1, 3
                      vj(m,mu,nu) = vj(m,mu,nu) - ket_fac*g(m)*dm(lam,kap)
                      vk_pre(m,mu,lam) = vk_pre(m,mu,lam) + g(m)*dm(nu,kap)
                      if (sk /= sl) vk_pre(m,mu,kap) = vk_pre(m,mu,kap) + g(m)*dm(nu,lam)
                      if (si /= sj) then
                        vj(m,nu,mu) = vj(m,nu,mu) + ket_fac*g(m)*dm(lam,kap)
                        vk_pre(m,nu,lam) = vk_pre(m,nu,lam) - g(m)*dm(mu,kap)
                        if (sk /= sl) vk_pre(m,nu,kap) = vk_pre(m,nu,kap) - g(m)*dm(mu,lam)
                      end if
                    end do
                  end do
                end do
              end do
            end do
            nullify(p0)
            nullify(pr)
          end do
        end do
      end do
    end do

    do m = 1, 3
      do i = 1, nbf
        do j = 1, nbf
          vk(m,i,j) = -(vk_pre(m,i,j) - vk_pre(m,j,i))
          h10(m,i,j) = vj(m,i,j) - 0.5_dp*vk(m,i,j)
        end do
      end do
    end do

    call gdat%clean()
    call int2_driver%clean()
    deallocate(vk_pre, eri0, erir)

  end subroutine giao_h10_twoe_matrix


  integer function raised_cart_index(idx, am, axis) result(match)
    use constants, only: cart_x, cart_y, cart_z, num_cart_bf
    integer, intent(in) :: idx, am, axis
    integer :: i, tx, ty, tz
    tx = cart_x(idx, am)
    ty = cart_y(idx, am)
    tz = cart_z(idx, am)
    if (axis == 1) tx = tx + 1
    if (axis == 2) ty = ty + 1
    if (axis == 3) tz = tz + 1
    match = 0
    do i = 1, num_cart_bf(am + 1)
      if (cart_x(i, am + 1) == tx .and. cart_y(i, am + 1) == ty .and. cart_z(i, am + 1) == tz) then
        match = i
        return
      end if
    end do
    error stop 'raised_cart_index failed'
  end function raised_cart_index

end module nmr_giao_debug_mod
