module nmr_giao_debug_mod

  implicit none
  private
  public nmr_giao_h10_debug
  public nmr_giao_h10_twoe_debug

contains

  subroutine nmr_giao_h10_debug_C(c_handle) bind(C, name="nmr_giao_h10_debug")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call nmr_giao_h10_debug(inf)
  end subroutine nmr_giao_h10_debug_C

  subroutine nmr_giao_h10_twoe_debug_C(c_handle) bind(C, name="nmr_giao_h10_twoe_debug")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call nmr_giao_h10_twoe_debug(inf)
  end subroutine nmr_giao_h10_twoe_debug_C

!> @brief Test-only/debug emitter for native GIAO one-electron h10 matrices.
!> @details This is deliberately separate from production NMR dispatch.  It
!>  computes the packed lower-triangle AO matrices for the one-electron GIAO
!>  h10 block and writes machine-parseable records to the calculation log so the
!>  kernel can be compared to an external PySCF/libcint oracle without ungating
!>  `nmr_gauge=giao` or routing production shieldings through partial GIAO code.
  subroutine nmr_giao_h10_debug(infos)

    use basis_tools, only: basis_set
    use constants, only: tol_int
    use int1, only: giao_h10_core, giao_h10_core_terms_full
    use io_constants, only: iw

    use precision, only: dp
    use types, only: information

    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis
    character(len=24), parameter :: term_name(8) = [ &
      'raw_irjxp               ', '-0.5_irjxp              ', &
      'raw_igkin               ', '-igkin                  ', &
      'raw_ignuc               ', 'ignuc_asym              ', &
      '-ignuc_asym             ', 'final_h10               ']
    real(kind=dp), allocatable :: h10(:,:), terms(:,:,:,:)
    real(kind=dp) :: tol
    integer :: i, j, m, ij, nbf, nbf2

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    allocate(h10(nbf2,3), terms(nbf,nbf,3,8))

    tol = log(10.0d0)*tol_int
    call giao_h10_core(basis, infos%atoms%xyz, &
         infos%atoms%zn - infos%basis%ecp_zn_num, h10, debug=.false., logtol=tol)
    call giao_h10_core_terms_full(basis, infos%atoms%xyz, &
         infos%atoms%zn - infos%basis%ecp_zn_num, size(infos%atoms%zn), terms, tol)

    open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/,A)') 'GIAO_H10_DEBUG_BEGIN packed_lower one-electron h10 real-coefficient'
    write(iw,'(A,1X,I0)') 'GIAO_H10_DEBUG_NBF', nbf
    write(iw,'(A,1X,I0)') 'GIAO_H10_DEBUG_NATOM', size(infos%atoms%zn)
    do i = 1, size(infos%atoms%zn)
      write(iw,'(A,1X,I0,1X,F6.1,3(1X,ES24.16))') 'GIAO_H10_DEBUG_ATOM_BOHR', &
        i, infos%atoms%zn(i), infos%atoms%xyz(1:3,i)
    end do
    ij = 0
    do i = 1, nbf
      do j = 1, i
        ij = ij + 1
        do m = 1, 3
          write(iw,'(A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') 'GIAO_H10_DEBUG_PACKED', &
            m, i, j, h10(ij,m)
        end do
      end do
    end do
    write(iw,'(A)') 'GIAO_H10_DEBUG_TERMS_BEGIN full_ao current-native-term-construction'
    do ij = 1, 8
      write(iw,'(A,1X,I0,1X,A)') 'GIAO_H10_DEBUG_TERM', ij, trim(term_name(ij))
      do m = 1, 3
        do i = 1, nbf
          do j = 1, nbf
            write(iw,'(A,1X,I0,1X,A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') &
              'GIAO_H10_DEBUG_FULL', ij, trim(term_name(ij)), m, i, j, terms(i,j,m,ij)
          end do
        end do
      end do
    end do
    write(iw,'(A)') 'GIAO_H10_DEBUG_TERMS_END'
    write(iw,'(A)') 'GIAO_H10_DEBUG_END'
    close(iw)

    deallocate(h10, terms)

  end subroutine nmr_giao_h10_debug

!> @brief Test-only/debug emitter for native RHF GIAO two-electron h10 terms.
!> @details RED/GREEN scaffold for the two-electron magnetic Fock derivative.
!>  The output format is intentionally machine-parseable and remains gated away
!>  from production shielding until it matches the PySCF/libcint get_jk oracle.
  subroutine nmr_giao_h10_twoe_debug(infos)

    use basis_tools, only: basis_set
    use constants, only: cart_x, cart_y, cart_z, num_cart_bf
    use int2_compute, only: int2_compute_t
    use int2e_rys, only: int2_rys_compute_ordered_am, int2_rys_data_t
    use io_constants, only: iw
    use messages, only: show_message, with_abort
    use oqp_tagarray_driver
    use precision, only: dp
    use types, only: information

    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis
    type(int2_compute_t), target :: int2_driver
    type(int2_rys_data_t) :: gdat
    real(kind=dp), allocatable :: vj(:,:,:), vk(:,:,:), vk_pre(:,:,:), h10(:,:,:)
    real(kind=dp), allocatable :: dm(:,:), dmat_packed(:)
    real(kind=dp), allocatable, target :: eri0(:), erir(:)
    real(kind=dp), pointer :: p0(:,:,:,:), pr(:,:,:,:)
    real(kind=dp), contiguous, pointer :: dmat_a(:)
    real(kind=dp) :: g(3), d(3), rij(3), norm4, ket_fac
    real(kind=dp) :: base_val, raised_val
    integer :: i, j, k, l, m, nbf, nbf2, nshell, maxang, maxcart
    integer :: si, sj, sk, sl, ni, nj, nk, nl, mu, nu, kap, lam
    integer :: ii, jj, kk, ll, axis, mapr, ids(4), am0(4), amr(4), ok
    integer(4) :: status
    logical :: zero_shq

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nshell = basis%nshell
    maxang = maxval(basis%am)
    if (maxang + 1 > 6) call show_message('GIAO two-electron debug: angular momentum too high', with_abort)
    maxcart = num_cart_bf(maxang + 1)

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a, status)
    if (status /= 0) call show_message('GIAO two-electron debug: OQP::DM_A missing', with_abort)
    allocate(dm(nbf,nbf), dmat_packed(nbf2), source=0.0_dp)
    dmat_packed = dmat_a
    call unpack_symmetric_density(dmat_packed, dm, nbf)

    allocate(vj(3,nbf,nbf), vk(3,nbf,nbf), vk_pre(3,nbf,nbf), h10(3,nbf,nbf), source=0.0_dp)
    allocate(eri0(maxcart**4), erir(maxcart**4), source=0.0_dp)

    call int2_driver%init(basis, infos)
    ! The magnetic derivative is ordered/antisymmetric in the bra pair.  Rebuild
    ! pair data without angular-momentum reordering so the Rys recurrence sees
    ! PA/PB vectors in the same shell order as the explicit ids below.
    call int2_driver%ppairs%compute(basis, int2_driver%cutoffs, noswap=.true.)
    call gdat%init(maxang + 1, int2_driver%cutoffs, ok)
    if (ok /= 0) call show_message('GIAO two-electron debug: cannot allocate Rys workspace', with_abort)

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
                      amr = am0
                      amr(1) = amr(1) + 1
                      erir = 0.0_dp
                      call int2_rys_compute_ordered_am(erir, gdat, int2_driver%ppairs, ids, amr, zero_shq)
                      pr(1:nl,1:nk,1:nj,1:num_cart_bf(amr(1))) => erir(1:nl*nk*nj*num_cart_bf(amr(1)))
                      mapr = raised_cart_index(ii, am0(1), axis)
                      d(axis) = (pr(ll,kk,jj,mapr) + basis%shell_centers(si,axis)*base_val) * norm4
                      nullify(pr)
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

    open(unit=iw, file=infos%log_filename, position="append")
    write(iw,'(/,A)') 'GIAO_H10_TWOE_DEBUG_BEGIN full_ao rhf-two-electron-h10'
    write(iw,'(A,1X,I0)') 'GIAO_H10_TWOE_DEBUG_NBF', nbf
    do m = 1, 3
      do i = 1, nbf
        do j = 1, nbf
          write(iw,'(A,1X,A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') &
            'GIAO_H10_TWOE_DEBUG_FULL', 'vj', m, i, j, vj(m,i,j)
          write(iw,'(A,1X,A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') &
            'GIAO_H10_TWOE_DEBUG_FULL', 'vk', m, i, j, vk(m,i,j)
          write(iw,'(A,1X,A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') &
            'GIAO_H10_TWOE_DEBUG_FULL', 'twoe_h10', m, i, j, h10(m,i,j)
        end do
      end do
    end do
    write(iw,'(A)') 'GIAO_H10_TWOE_DEBUG_END'
    close(iw)

    call gdat%clean()
    call int2_driver%clean()
    deallocate(vj, vk, vk_pre, h10, dm, dmat_packed, eri0, erir)

  end subroutine nmr_giao_h10_twoe_debug

  subroutine unpack_symmetric_density(packed, full, n)
    use precision, only: dp
    real(kind=dp), intent(in) :: packed(:)
    real(kind=dp), intent(out) :: full(:,:)
    integer, intent(in) :: n
    integer :: i, j, ij
    full = 0.0_dp
    do i = 1, n
      do j = 1, i
        ij = j + i*(i-1)/2
        full(i,j) = packed(ij)
        full(j,i) = packed(ij)
      end do
    end do
  end subroutine unpack_symmetric_density

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
