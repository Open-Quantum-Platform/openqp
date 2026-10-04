!> @brief ECP (effective core potential) one-electron integrals for OpenQP.
!> @detail Contracts the raw Cartesian ECP integrals and nuclear derivatives of
!>         @ref ecp_native with densities, or returns them in OpenQP AO order
!>         (canonical Cartesian order, spherical transformation).
!> @author Mohsen Mazaherifar
!> @date January 2025
module ecp_tool
    use iso_c_binding, only: c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only: real64
    use basis_tools, only: basis_set
    use precision, only: dp
    use constants, only: HARMONIC_ACTIVE, NUM_CART_BF

    implicit none

    private
    public add_ecpint
    public add_ecpder
    public add_ecphess
    public ecp_deriv_ints

contains
    !> @brief Add ECP one-electron contribution to the AO-core Hamiltonian (packed).
    !> @detail Computes scalar ECP integrals with ecp_native (deriv order 0),
    !>         remaps them into OpenQP AO ordering via @ref transform_ecp_matrix,
    !>         and accumulates into upper-triangular packed Hcore.
    !> @param[in]  basis   Basis set (contains ECP params and AO metadata).
    !> @param[in]  coord   Nuclear coordinates (3×natm).
    !> @param[inout] hcore Upper-triangular packed AO core Hamiltonian (size nbf*(nbf+1)/2).
    !> @note No-op if basis%ecp_params%is_ecp == .false.
    !> @author Mohsen Mazaherifar
    !> @date January 2025
    subroutine add_ecpint(basis, coord, hcore)
        real(real64), contiguous, intent(in) :: coord(:,:)
        type(basis_set), intent(in) :: basis
        real(real64), contiguous, intent(inout) :: hcore(:)
        real(c_double), allocatable :: raw_res(:)
        real(c_double), allocatable :: ecp_mat(:)
        integer :: i, j, c
        integer(c_int) :: driv_order

        if (.not.(basis%ecp_params%is_ecp)) then
            return
        end if
        driv_order = 0

        call ecp_raw(basis, coord, driv_order, raw_res)


        call transform_ecp_matrix(basis, raw_res, ecp_mat)

        c = 0
        do i = 1, basis%nbf
            do j = 1, i
                c = c + 1
                hcore(c) = ecp_mat((i - 1) * basis%nbf + j) + hcore(c)
            end do
        end do

        deallocate(raw_res)
        deallocate(ecp_mat)


    end subroutine add_ecpint

    !> @brief Add ECP force contribution (first derivatives) to nuclear gradients.
    !> @detail Computes dV_ECP/dR_A in AO full-square form for each atom using
    !>         ecp_native (deriv order 1), transforms to OpenQP AO ordering, and
    !>         contracts with the symmetric density `denab` (packed) to accumulate
    !>         into atomic gradient components `de(:,A)`.
    !> @param[in]    basis  Basis set (with ECP params).
    !> @param[in]    coord  Nuclear coordinates (3×natm).
    !> @param[inout] denab  Packed AO density (size nbf*(nbf+1)/2).
    !> @param[inout] de     Nuclear gradients (3×natm), incremented by ECP part.
    !> @note No-op if basis%ecp_params%is_ecp == .false.
    !> @author Mohsen Mazaherifar
    !> @date January 2025
    subroutine add_ecpder(basis, coord, denab, de)

        real(real64), contiguous, intent(in) :: coord(:,:)
        type(basis_set), intent(in) :: basis
        REAL(kind=dp), INTENT(INOUT) :: denab(:)
        REAL(kind=dp), intent(INOUT) :: de(:,:)

        real(c_double), allocatable :: raw_res(:)
        real(c_double), allocatable :: raw_block(:), ecp_mat(:)
        real(real64), allocatable :: deloc(:,:)
        integer :: i, j, c, n, natm, prim, cc, nbf_raw
        ! 64-bit: slice offsets reach 3*natm*nbf^2 and overflow default integers
        integer(c_int64_t) :: full_size
        integer(c_int) :: driv_order

        if (.not.(basis%ecp_params%is_ecp)) then
            return
        end if

        driv_order = 1

        nbf_raw = ecp_cart_nbf(basis)
        full_size = int(nbf_raw, c_int64_t) * nbf_raw
        allocate(raw_block(full_size))

        natm = size(coord, dim=2)

        allocate(deloc(3, natm))
        deloc = 0

        call ecp_raw(basis, coord, driv_order, raw_res)


        do n = 1, natm
            do cc = 1, 3
                raw_block = raw_res(full_size * (3 * (n - 1) + cc - 1) + 1 : &
                                       full_size * (3 * (n - 1) + cc))
                call transform_ecp_matrix(basis, raw_block, ecp_mat)

                do j = 1, basis%nbf
                    do i = 1, j
                        c = j * (j - 1) / 2 + i

                        if (i == j) then
                            prim = 1
                        else
                            prim = 2
                        end if

                        deloc(cc, n) = deloc(cc, n) + prim * ecp_mat((i - 1) * basis%nbf + j) * denab(c)
                    end do
                end do
            end do

        end do

        de(:, 1:natm) = de(:, 1:natm) + deloc(:, 1:natm)

        deallocate(raw_res)
        if (allocated(ecp_mat)) deallocate(ecp_mat)
        deallocate(raw_block)


    end subroutine add_ecpder

    !> @brief Return ECP one-electron first-derivative integrals (uncontracted).
    !> @detail Computes dV_ECP_{mu,nu}/dR_{I,c} for every atom I and Cartesian
    !>         direction c using ecp_native (deriv order 1), transforms each block
    !>         to OpenQP AO ordering, and stores the full-square AO matrices into
    !>         `dVecp(mu,nu,c,I)`.  These are the response counterpart of
    !>         @ref add_ecpder (which contracts the same integrals with a density);
    !>         the analytic Hessian adds them into the core-Hamiltonian derivative
    !>         dHcore/dR so the ECP enters the CPHF right-hand side and the
    !>         orbital-relaxation response, exactly as nuclear attraction does.
    !>         Like @ref add_ecpint, the integrals are returned in the OpenQP
    !>         normalized (density/Hcore) convention, so callers must NOT apply an
    !>         additional bfnrm scaling.
    !> @param[in]  basis  Basis set (with ECP params).
    !> @param[in]  coord  Nuclear coordinates (3 x natm).
    !> @param[out] dVecp  ECP derivative integrals (nbf x nbf x 3 x natm).
    !> @note Returns zeros if basis%ecp_params%is_ecp == .false.
    subroutine ecp_deriv_ints(basis, coord, dVecp)

        real(real64), contiguous, intent(in) :: coord(:,:)
        type(basis_set), intent(in) :: basis
        real(kind=dp), intent(out) :: dVecp(:,:,:,:)

        real(c_double), allocatable :: raw_res(:)
        real(c_double), allocatable :: raw_block(:), ecp_mat(:)
        integer :: nbf, nbf_raw, natm, n, cc, i, j
        ! 64-bit: slice offsets reach 3*natm*nbf^2 and overflow default integers
        integer(c_int64_t) :: full_size
        integer(c_int) :: driv_order

        dVecp = 0.0_dp
        if (.not.(basis%ecp_params%is_ecp)) then
            return
        end if

        driv_order = 1
        nbf = basis%nbf
        nbf_raw = ecp_cart_nbf(basis)
        full_size = int(nbf_raw, c_int64_t) * nbf_raw
        natm = size(coord, dim=2)
        allocate(raw_block(full_size))

        call ecp_raw(basis, coord, driv_order, raw_res)

        do n = 1, natm
            do cc = 1, 3
                raw_block = raw_res(full_size*(3*(n - 1) + cc - 1) + 1 : &
                                       full_size*(3*(n - 1) + cc))
                call transform_ecp_matrix(basis, raw_block, ecp_mat)
                do j = 1, nbf
                    do i = 1, nbf
                        dVecp(i, j, cc, n) = ecp_mat((i - 1)*nbf + j)
                    end do
                end do
            end do
        end do

        deallocate(raw_res)
        if (allocated(ecp_mat)) deallocate(ecp_mat)

        deallocate(raw_block)

    end subroutine ecp_deriv_ints

    !> @brief Add ECP second-derivative contribution to the nuclear Hessian.
    !> @detail Computes d^2 V_ECP/dR_I dR_J in AO full-square form for every atom
    !>         pair using ecp_native (deriv order 2), transforms each block to
    !>         OpenQP AO ordering, and contracts with the symmetric density
    !>         `denab` (packed) to accumulate the fixed-density ECP skeleton into
    !>         the Cartesian Hessian `hess` (3*natm x 3*natm, atom-major layout
    !>         hess(3*(I-1)+a, 3*(J-1)+b)).
    !>
    !>         ecp_native returns the packed upper triangle of atom-coordinate
    !>         pairs: matrix index H_START(I,J,natm) (0-based) starts each (I<=J)
    !>         atom block.  Diagonal blocks (I==J) store 6 matrices in the order
    !>         {xx,xy,xz,yy,yz,zz}; off-diagonal blocks (I<J) store 9 matrices in
    !>         row-major {xx,xy,xz,yx,yy,yz,zx,zy,zz} (first index = coordinate of
    !>         atom I, second = coordinate of atom J).  Each block is the AO matrix
    !>         packed M(k,l) = (k-1)*nbf + l.  We scatter symmetrically so the
    !>         returned Hessian is exactly symmetric.
    !> @param[in]    basis  Basis set (with ECP params).
    !> @param[in]    coord  Nuclear coordinates (3 x natm).
    !> @param[in]    denab  Packed AO density (size nbf*(nbf+1)/2), upper triangle.
    !> @param[inout] hess   Cartesian Hessian (3*natm x 3*natm), incremented by ECP.
    !> @note No-op if basis%ecp_params%is_ecp == .false.
    subroutine add_ecphess(basis, coord, denab, hess)

        real(real64), contiguous, intent(in) :: coord(:,:)
        type(basis_set), intent(in) :: basis
        real(kind=dp), intent(in) :: denab(:)
        real(kind=dp), intent(inout) :: hess(:,:)

        real(c_double), allocatable :: raw_res(:)
        real(c_double), allocatable :: raw_block(:), ecp_mat(:)
        integer :: nbf, nbf_raw, natm
        integer :: iat, jat, ia0, ja0, hstart, base, ncomp, n
        integer :: a, b, i, j, c, prim
        ! 64-bit: block offsets reach 3N(3N+1)/2 * nbf^2 and overflow default integers
        integer(c_int64_t) :: mat_sz
        integer(c_int) :: driv_order
        integer :: amap(9), bmap(9)
        real(real64) :: val

        if (.not.(basis%ecp_params%is_ecp)) then
            return
        end if

        driv_order = 2
        nbf = basis%nbf
        nbf_raw = ecp_cart_nbf(basis)
        mat_sz = int(nbf_raw, c_int64_t) * nbf_raw
        natm = size(coord, dim=2)
        allocate(raw_block(mat_sz))

        call ecp_raw(basis, coord, driv_order, raw_res)

        do iat = 1, natm
            do jat = iat, natm
                ia0 = iat - 1
                ja0 = jat - 1
                ! 0-based starting matrix index of the (iat,jat) atom block
                hstart = 9*ja0 + 3*(3*natm - 1)*ia0 - (9*ia0*(ia0 + 1))/2 - 3
                if (iat == jat) then
                    base = hstart + 3
                    ncomp = 6
                    amap(1:6) = [1, 1, 1, 2, 2, 3]
                    bmap(1:6) = [1, 2, 3, 2, 3, 3]
                else
                    base = hstart
                    ncomp = 9
                    amap(1:9) = [1, 1, 1, 2, 2, 2, 3, 3, 3]
                    bmap(1:9) = [1, 2, 3, 1, 2, 3, 1, 2, 3]
                end if

                do n = 1, ncomp
                    raw_block = raw_res((base + n - 1)*mat_sz + 1 : (base + n - 1)*mat_sz + mat_sz)
                    call transform_ecp_matrix(basis, raw_block, ecp_mat)

                    val = 0.0_dp
                    do j = 1, nbf
                        do i = 1, j
                            c = j*(j - 1)/2 + i
                            if (i == j) then
                                prim = 1
                            else
                                prim = 2
                            end if
                            val = val + prim * ecp_mat((i - 1)*nbf + j) * denab(c)
                        end do
                    end do

                    a = amap(n)
                    b = bmap(n)
                    hess(3*(iat - 1) + a, 3*(jat - 1) + b) = &
                        hess(3*(iat - 1) + a, 3*(jat - 1) + b) + val
                    ! symmetric partner (skip if it is the same matrix element)
                    if (.not. (iat == jat .and. a == b)) then
                        hess(3*(jat - 1) + b, 3*(iat - 1) + a) = &
                            hess(3*(jat - 1) + b, 3*(iat - 1) + a) + val
                    end if
                end do
            end do
        end do

        deallocate(raw_res)
        if (allocated(ecp_mat)) deallocate(ecp_mat)

        deallocate(raw_block)

    end subroutine add_ecphess

    !> @brief Raw Cartesian ECP integrals or nuclear derivatives (see ecp_native).
    subroutine ecp_raw(basis, coord, deriv_order, raw)
        use ecp_native, only: ecp_native_raw
        type(basis_set), intent(in) :: basis
        real(real64), contiguous, intent(in) :: coord(:,:)
        integer(c_int), intent(in) :: deriv_order
        real(c_double), allocatable, intent(out) :: raw(:)

        call ecp_native_raw(basis, coord, int(deriv_order), raw)

    end subroutine ecp_raw
  !> @brief Build AO index remapping from the raw Cartesian order (x power descending,
  !>        then y) to OpenQP AO order.
  !> @detail Fills `label_map(i_old)=i_new` using shell origins and angular-momentum
  !>         layout so that full-square AO matrices can be permuted consistently.
  !> @param[in]    basis     Basis set (AO layout and shell metadata).
  !> @param[inout] label_map Integer array of length nbf receiving the permutation.
  !> @see transform_ecp_matrix
  !> @author Mohsen Mazaherifar
  !> @date January 2025
  integer function ecp_cart_nbf(basis) result(nbf_cart)

    type(basis_set), intent(in) :: basis
    integer :: ish

    nbf_cart = 0
    do ish = 1, basis%nshell
      nbf_cart = nbf_cart + NUM_CART_BF(basis%am(ish))
    end do

  end function ecp_cart_nbf

  subroutine ecp_cart_offsets(basis, cart_off, nbf_cart)

    type(basis_set), intent(in) :: basis
    integer, allocatable, intent(out) :: cart_off(:)
    integer, intent(out) :: nbf_cart
    integer :: ish

    allocate(cart_off(basis%nshell))
    nbf_cart = 0
    do ish = 1, basis%nshell
      cart_off(ish) = nbf_cart + 1
      nbf_cart = nbf_cart + NUM_CART_BF(basis%am(ish))
    end do

  end subroutine ecp_cart_offsets

  subroutine raw_cart_map(basis, cart_off, label_map)

    use basis_tools, only: basis_set
    use constants, only: map_canonical

    type(basis_set), intent(in) :: basis
    integer, dimension(:), intent(in) :: cart_off
    integer, dimension(:), intent(inout):: label_map
    integer :: ish, i, old

    label_map = 0
    do ish = 1, basis%nshell
      do i = 1, NUM_CART_BF(basis%am(ish))
        old = cart_off(ish) + i - 1
        label_map(old + map_canonical(i, basis%am(ish))) = old
      end do
    end do

  end  subroutine raw_cart_map

  !> @brief Permute a full AO square matrix into OpenQP AO ordering.
  !> @detail Applies the mapping from @ref raw_cart_map to reorder rows/cols
  !>         of `matrix` in-place (via a temporary copy). Expects size nbf×nbf.
  !> @param[in]    basis   Basis set (provides AO label map).
  !> @param[inout] matrix  Full AO square matrix flattened (size nbf*nbf).
  !> @throws Stops if `size(matrix) != nbf*nbf`.
  !> @see raw_cart_map
  !> @author Mohsen Mazaherifar
  !> @date January 2025
  subroutine transform_ecp_matrix(basis, raw_matrix, matrix)

    use basis_tools, only: basis_set
    use cart2sph, only: cart2sph_mat
    type(basis_set), intent(in) :: basis
    real(c_double), dimension(:), intent(in) :: raw_matrix
    real(c_double), dimension(:), allocatable, intent(out) :: matrix
    real(c_double), dimension(:), allocatable :: cart_matrix
    real(c_double), dimension(:), allocatable :: blk
    integer, dimension(:), allocatable :: label_map
    integer, allocatable :: cart_off(:)
    integer :: i, j, row, col, nbf_raw, nbf_sph
    integer :: ish, jsh, nci, ncj, nsi, nsj, coi, coj, soi, soj
    integer :: si, sj, max_blk, pure_i, pure_j

    call ecp_cart_offsets(basis, cart_off, nbf_raw)
    nbf_sph = basis%nbf
    allocate(label_map(nbf_raw))

    if (size(raw_matrix) /= nbf_raw * nbf_raw) then
      print *, "Error: original_matrix size does not match labels."
      stop
    end if

    call raw_cart_map(basis, cart_off, label_map)
    allocate(cart_matrix(nbf_raw * nbf_raw))
    allocate(matrix(nbf_sph * nbf_sph))

    cart_matrix = 0.0_dp
    matrix = 0.0_dp
    do i = 1, nbf_raw
      do j = 1, nbf_raw
        row = label_map(i)
        col = label_map(j)
        cart_matrix((row - 1) * nbf_raw + col) = raw_matrix((i - 1) * nbf_raw + j)
      end do
    end do

    max_blk = 0
    do ish = 1, basis%nshell
      do jsh = 1, basis%nshell
        max_blk = max(max_blk, NUM_CART_BF(basis%am(ish)) * NUM_CART_BF(basis%am(jsh)))
      end do
    end do
    allocate(blk(max_blk))

    do ish = 1, basis%nshell
      nci = NUM_CART_BF(basis%am(ish))
      nsi = basis%naos(ish)
      coi = cart_off(ish)
      soi = basis%ao_offset(ish)
      if (HARMONIC_ACTIVE) then
        pure_i = basis%harmonic(ish)
      else
        pure_i = 0
      end if

      do jsh = 1, basis%nshell
        ncj = NUM_CART_BF(basis%am(jsh))
        nsj = basis%naos(jsh)
        coj = cart_off(jsh)
        soj = basis%ao_offset(jsh)
        if (HARMONIC_ACTIVE) then
          pure_j = basis%harmonic(jsh)
        else
          pure_j = 0
        end if

        do si = 1, nci
          do sj = 1, ncj
            blk((si - 1) * ncj + sj) = cart_matrix((coi + si - 2) * nbf_raw + coj + sj - 1)
          end do
        end do

        ! Raw ECP blocks are in the same pure-power Cartesian convention as
        ! the native 1e primitives (bas_norm_matrix folds shells_pnrm2 for
        ! Cartesian shells later, but bfnrm = 1 for pure shells), so the
        ! transform must fold shells_pnrm2 along each pure index itself.
        call cart2sph_mat(blk, basis%am(jsh), pure_j, basis%am(ish), pure_i)

        do si = 1, nsi
          do sj = 1, nsj
            matrix((soi + si - 2) * nbf_sph + soj + sj - 1) = blk((si - 1) * nsj + sj)
          end do
        end do
      end do
    end do

  end subroutine transform_ecp_matrix


end module ecp_tool
