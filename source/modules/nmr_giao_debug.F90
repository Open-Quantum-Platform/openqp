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

    use basis_tools, only: basis_set, bas_norm_matrix, build_cart_density
    use constants, only: cart_x, cart_y, cart_z, HARMONIC_ACTIVE, num_cart_bf
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
    real(kind=dp), allocatable :: dm_norm(:,:), dm_cart(:,:), dm_work(:,:)
    real(kind=dp), allocatable :: vj_work(:,:,:), vk_work(:,:,:), h10_work(:,:,:), vk_pre(:,:,:)
    real(kind=dp), allocatable, target :: eri0(:), erir(:)
    real(kind=dp), pointer :: p0(:,:,:,:), pr(:,:,:,:)
    real(kind=dp) :: g(3), d(3), rij(3), norm4, ket_fac, base_val
    integer, allocatable :: loc(:), cart_off(:)
    integer :: i, j, m, nbf, nshell, maxang, maxcart, nbf_work
    integer :: si, sj, sk, sl, ni, nj, nk, nl, mu, nu, kap, lam
    integer :: ii, jj, kk, ll, axis, mapr, ids(4), am0(4), amr(4), ok
    logical :: zero_shq, usecart

    nbf = basis%nbf
    nshell = basis%nshell
    maxang = maxval(basis%am)
    if (maxang + 1 > 6) call show_message('GIAO two-electron: angular momentum too high', with_abort)
    maxcart = num_cart_bf(maxang + 1)

    usecart = HARMONIC_ACTIVE
    if (usecart) then
      dm_norm = dm
      call bas_norm_matrix(dm_norm, basis%bfnrm, basis%nbf)
      call build_cart_density(basis, dm_norm, dm_cart, cart_off, nbf_work)
      dm_work = dm_cart
      loc = cart_off
    else
      nbf_work = nbf
      dm_work = dm
      allocate(loc(nshell))
      loc = basis%ao_offset
    end if

    vj = 0.0_dp
    vk = 0.0_dp
    h10 = 0.0_dp
    allocate(vj_work(3,nbf_work,nbf_work), source=0.0_dp)
    allocate(vk_work(3,nbf_work,nbf_work), h10_work(3,nbf_work,nbf_work), source=0.0_dp)
    allocate(vk_pre(3,nbf_work,nbf_work), source=0.0_dp)
    allocate(eri0(maxcart**4), erir(maxcart**4), source=0.0_dp)

    call int2_driver%init(basis, infos)
    int2_driver%rys_only = .true.   ! NMR is Rys-only by construction
    ! The magnetic derivative is ordered/antisymmetric in the bra pair.  Rebuild
    ! pair data without angular-momentum reordering so the Rys recurrence sees
    ! PA/PB vectors in the same shell order as the explicit ids below.
    call int2_driver%ppairs%compute(basis, int2_driver%cutoffs, noswap=.true.)
    call gdat%init(maxang + 1, int2_driver%cutoffs, ok)
    if (ok /= 0) call show_message('GIAO two-electron: cannot allocate Rys workspace', with_abort)

    do si = 1, nshell
      if (usecart) then
        ni = NUM_CART_BF(basis%am(si))
      else
        ni = basis%naos(si)
      end if
      do sj = 1, si
        if (usecart) then
          nj = NUM_CART_BF(basis%am(sj))
        else
          nj = basis%naos(sj)
        end if
        rij = basis%shell_centers(si,1:3) - basis%shell_centers(sj,1:3)
        if (maxval(abs(rij)) < 1.0d-14) cycle
        do sk = 1, nshell
          if (usecart) then
            nk = NUM_CART_BF(basis%am(sk))
          else
            nk = basis%naos(sk)
          end if
          do sl = 1, sk
            if (usecart) then
              nl = NUM_CART_BF(basis%am(sl))
            else
              nl = basis%naos(sl)
            end if
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
              mu = loc(si) + ii - 1
              do jj = 1, nj
                nu = loc(sj) + jj - 1
                do kk = 1, nk
                  kap = loc(sk) + kk - 1
                  do ll = 1, nl
                    lam = loc(sl) + ll - 1
                    norm4 = 1.0_dp
                    if (.not. usecart) norm4 = basis%bfnrm(mu)*basis%bfnrm(nu)*basis%bfnrm(kap)*basis%bfnrm(lam)
                    base_val = p0(ll,kk,jj,ii)
                    do axis = 1, 3
                      mapr = raised_cart_index(ii, am0(1), axis)
                      d(axis) = (pr(ll,kk,jj,mapr) + basis%shell_centers(si,axis)*base_val) * norm4
                    end do
                    g(1) = -0.5_dp*(rij(2)*d(3) - rij(3)*d(2))
                    g(2) = -0.5_dp*(rij(3)*d(1) - rij(1)*d(3))
                    g(3) = -0.5_dp*(rij(1)*d(2) - rij(2)*d(1))
                    do m = 1, 3
                      vj_work(m,mu,nu) = vj_work(m,mu,nu) - ket_fac*g(m)*dm_work(lam,kap)
                      vk_pre(m,mu,lam) = vk_pre(m,mu,lam) + g(m)*dm_work(nu,kap)
                      if (sk /= sl) vk_pre(m,mu,kap) = vk_pre(m,mu,kap) + g(m)*dm_work(nu,lam)
                      if (si /= sj) then
                        vj_work(m,nu,mu) = vj_work(m,nu,mu) + ket_fac*g(m)*dm_work(lam,kap)
                        vk_pre(m,nu,lam) = vk_pre(m,nu,lam) - g(m)*dm_work(mu,kap)
                        if (sk /= sl) vk_pre(m,nu,kap) = vk_pre(m,nu,kap) - g(m)*dm_work(mu,lam)
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
      do i = 1, nbf_work
        do j = 1, nbf_work
          vk_work(m,i,j) = -(vk_pre(m,i,j) - vk_pre(m,j,i))
          h10_work(m,i,j) = vj_work(m,i,j) - 0.5_dp*vk_work(m,i,j)
        end do
      end do
    end do

    if (usecart) then
      call reduce_cart_giao_matrix(basis, cart_off, vj_work, vj)
      call reduce_cart_giao_matrix(basis, cart_off, vk_work, vk)
      call reduce_cart_giao_matrix(basis, cart_off, h10_work, h10)
    else
      vj = vj_work
      vk = vk_work
      h10 = h10_work
    end if

    call gdat%clean()
    call int2_driver%clean()
    deallocate(vj_work, vk_work, h10_work, vk_pre, eri0, erir)

  end subroutine giao_h10_twoe_matrix


  subroutine reduce_cart_giao_matrix(basis, cart_off, cart, sph)

    use basis_tools, only: basis_set, bas_norm_matrix
    use cart2sph, only: cart2sph_mat
    use constants, only: num_cart_bf
    use precision, only: dp

    implicit none

    type(basis_set), intent(in) :: basis
    integer, intent(in) :: cart_off(:)
    real(kind=dp), intent(in) :: cart(:,:,:)
    real(kind=dp), intent(out) :: sph(:,:,:)

    real(kind=dp), allocatable :: blk(:)
    integer :: m, si, sj, ii, jj, nci, ncj, nsi, nsj, coi, coj, soi, soj, idx

    sph = 0.0_dp
    do m = 1, 3
      do sj = 1, basis%nshell
        ncj = NUM_CART_BF(basis%am(sj))
        nsj = basis%naos(sj)
        coj = cart_off(sj)
        soj = basis%ao_offset(sj)
        do si = 1, basis%nshell
          nci = NUM_CART_BF(basis%am(si))
          nsi = basis%naos(si)
          coi = cart_off(si)
          soi = basis%ao_offset(si)
          allocate(blk(nci*ncj))
          idx = 0
          do jj = 1, ncj
            do ii = 1, nci
              idx = idx + 1
              blk(idx) = cart(m, coi + ii - 1, coj + jj - 1)
            end do
          end do
          if (basis%harmonic(si) == 1 .or. basis%harmonic(sj) == 1) &
              call cart2sph_mat(blk, basis%am(si), basis%harmonic(si), basis%am(sj), basis%harmonic(sj))
          idx = 0
          do jj = 1, nsj
            do ii = 1, nsi
              idx = idx + 1
              sph(m, soi + ii - 1, soj + jj - 1) = blk(idx)
            end do
          end do
          deallocate(blk)
        end do
      end do
      call bas_norm_matrix(sph(m,:,:), basis%bfnrm, basis%nbf)
    end do

  end subroutine reduce_cart_giao_matrix


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
