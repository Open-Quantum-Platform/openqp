!> @brief Effective-core-potential (ECP) one-electron integrals and their nuclear
!>        derivatives, and their contraction or transformation to OpenQP AO order.
!> @detail Semilocal ECP
!>
!>           U(r) = U_L(r_C) + sum_{l<L} sum_m |lm> U_l(r_C) <lm|,
!>           U_l(r) = sum_t d_t r^(n_t-2) exp(-zeta_t r^2),
!>
!>         between Cartesian Gaussian primitives.  The angular integrations are
!>         analytic: each primitive is expanded about the ECP centre C with
!>
!>           exp(k.r) = 4 pi sum_lam i_lam(k r) sum_mu Y_lam,mu(k^) Y_lam,mu(r^),
!>
!>         and the products of monomials and real spherical harmonics are
!>         integrated over the unit sphere in closed form.  The remaining radial
!>         integrals are evaluated with a fixed Gauss-Legendre rule on the window
!>         of the Gaussian that multiplies them.  After the exponential growth of
!>         i_lam is folded into that Gaussian, the rest of the integrand
!>         (r^N e^(-x) i_lam(x)) is smooth on the window, so no adaptive
!>         convergence test is needed and the value, first and second
!>         derivatives are all evaluated on the same abscissae.
!>
!>         Type 1 (local channel L): the AO pair is combined into one Gaussian at
!>         P = (a A + b B)/(a+b) before the expansion.  Type 2 (l < L): each AO is
!>         projected onto Y_lm separately.
!>
!>         Nuclear derivatives follow from
!>           d/dA_k [(x-A)^t e^(-a|r-A|^2)] = -t_k (x-A)^(t-e_k) e^.. + 2a (x-A)^(t+e_k) e^..,
!>         and translational invariance for the ECP centre.  Contractions are
!>         carried with the primitive weights 1, a, b, a^2, ab, b^2, so the
!>         exponent-independent angular work is done once per shell pair.
!>
!>         Method: the semilocal-projector expansion of McMurchie & Davidson,
!>         J. Comput. Phys. 44, 289 (1981), in the form of Flores-Moreno et al.,
!>         J. Comput. Chem. 27, 1009 (2006) and Shaw & Hill, J. Chem. Phys. 147,
!>         074108 (2017), with a different radial quadrature (above).
!>
!>         Raw result layout (ecp_raw_ints): Cartesian components ordered
!>         xx, xy, xz, yy, yz, zz (x power descending, then y), shells in basis
!>         order, basis%cc as the multipliers of the unnormalised monomials, one
!>         nraw x nraw matrix per derivative component:
!>           order 0: one matrix;
!>           order 1: 3*natm matrices, atom-major (x, y, z per atom);
!>           order 2: the upper triangle of atom-coordinate pairs, atom blocks
!>                    (I<=J); diagonal blocks {xx,xy,xz,yy,yz,zz}, off-diagonal
!>                    blocks row-major in (coordinate of I, coordinate of J).
!>
!>         add_ecpint, add_ecpder, ecp_deriv_ints and add_ecphess contract these
!>         raw matrices with densities or return them in OpenQP AO order
!>         (canonical Cartesian order, spherical transformation).
!> @author Mohsen Mazaherifar (AO-order drivers, January 2025)
module ecp_tool

  use iso_c_binding, only: c_double, c_int, c_int64_t
  use, intrinsic :: iso_fortran_env, only: real64
  use precision, only: dp
  use basis_tools, only: basis_set
  use constants, only: HARMONIC_ACTIVE, NUM_CART_BF
  use messages, only: show_message, WITH_ABORT

  implicit none

  private
  public add_ecpint
  public add_ecpder
  public add_ecphess
  public ecp_deriv_ints
  public ecp_raw_ints
  public ecp_hess_start
  public ecp_hess_contract

  !> Gauss-Legendre points per radial window
  integer, parameter :: NGL = 64
  !> The window spans |r - r0| <= sqrt(WIN2/q) for a Gaussian exp(-q (r-r0)^2)
  real(dp), parameter :: WIN2 = 60.0_dp
  !> Primitive/term contributions whose Gaussian amplitude is below exp(LOGSKIP)
  !> are skipped
  real(dp), parameter :: LOGSKIP = -90.0_dp
  !> Distance (sum of |dx|) below which a shell centre coincides with the ECP centre
  real(dp), parameter :: ONCENTRE = 1.0e-12_dp
  real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
  real(dp), parameter :: FOURPI = 4.0_dp*PI

  !> Per-atom angular projection table for one ECP centre (type 2)
  type gtab_t
    logical :: used = .false.
    logical :: oncentre = .false.
    integer :: namax = -1, lammax = -1
    real(dp), allocatable :: g(:,:,:,:)   ! (tuple, lm, N, lam)
  end type gtab_t

  !> Tables shared by all shell pairs of one call
  type tables_t
    integer :: lb = 0, le = 0, ly = 0, pm = 0
    real(dp), allocatable :: ycoef(:,:)       ! (tuple, lm) real orthonormal Y_lm
    real(dp), allocatable :: mom(:,:,:)       ! unit-sphere monomial integrals
    real(dp), allocatable :: om1(:,:)         ! (tuple, lam mu)
    real(dp), allocatable :: om2(:,:,:)       ! (tuple, lm, lam mu)
    real(dp) :: glx(NGL), glw(NGL)
    integer, allocatable :: tx(:), ty(:), tz(:), tn(:)   ! components and degree of each tuple
  end type tables_t

contains
    !> @brief Add ECP one-electron contribution to the AO-core Hamiltonian (packed).
    !> @detail Computes scalar ECP integrals with ecp_raw_ints (deriv order 0),
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

        call ecp_raw_ints(basis, coord, int(driv_order), raw_res)


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
    !>         ecp_raw_ints (deriv order 1), transforms to OpenQP AO ordering, and
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

        call ecp_raw_ints(basis, coord, int(driv_order), raw_res)


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
    !>         direction c using ecp_raw_ints (deriv order 1), transforms each block
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

        call ecp_raw_ints(basis, coord, int(driv_order), raw_res)

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
    !> @detail Accumulates the fixed-density ECP skeleton
    !>
    !>           hess(3(I-1)+a, 3(J-1)+b) += sum_{mu nu} D_{mu nu} d^2 V_{mu nu} / dR_Ia dR_Jb
    !>
    !>         without forming the second-derivative matrices.  The packed AO
    !>         density is taken once to the raw Cartesian order of ecp_raw_ints
    !>         with the adjoint of @ref transform_ecp_matrix (per shell pair
    !>         D_cart = B'_i D B'_j^T with B' = B * shells_pnrm2 for pure shells,
    !>         then the inverse canonical-order permutation), and contracted per
    !>         shell pair by ecp_hess_contract.  Memory is O(nbf^2) instead of
    !>         O(nbf^2 * natm^2).
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

        real(dp), allocatable :: draw(:,:)

        if (.not.(basis%ecp_params%is_ecp)) then
            return
        end if

        call ecp_density_to_raw(basis, denab, draw)
        call ecp_hess_contract(basis, coord, draw, hess)

    end subroutine add_ecphess

    !> @brief Packed OpenQP AO density -> symmetric density in the raw Cartesian
    !>        order of ecp_raw_ints, such that
    !>        sum_kl draw(k,l) V_raw(k,l) = sum_{mu nu} D(mu,nu) V(mu,nu)
    !>        for V = transform_ecp_matrix(V_raw).
    subroutine ecp_density_to_raw(basis, denab, draw)
        use cart2sph, only: c2s_expand_block

        type(basis_set), intent(in) :: basis
        real(kind=dp), intent(in) :: denab(:)
        real(dp), allocatable, intent(out) :: draw(:,:)

        real(dp), allocatable :: dfull(:,:), dcart(:,:)
        integer, allocatable :: cart_off(:), label_map(:)
        integer :: nbf, nbf_raw, i, j, ish, jsh, nci, ncj, nsi, nsj
        integer :: coi, coj, soi, soj, pure_i, pure_j

        nbf = basis%nbf
        allocate(dfull(nbf, nbf))
        do j = 1, nbf
            do i = 1, j
                dfull(i, j) = denab(j*(j - 1)/2 + i)
                dfull(j, i) = dfull(i, j)
            end do
        end do

        call ecp_cart_offsets(basis, cart_off, nbf_raw)
        allocate(dcart(nbf_raw, nbf_raw))
        do ish = 1, basis%nshell
            nci = NUM_CART_BF(basis%am(ish))
            nsi = basis%naos(ish)
            coi = cart_off(ish)
            soi = basis%ao_offset(ish)
            pure_i = 0
            if (HARMONIC_ACTIVE) pure_i = basis%harmonic(ish)
            do jsh = 1, basis%nshell
                ncj = NUM_CART_BF(basis%am(jsh))
                nsj = basis%naos(jsh)
                coj = cart_off(jsh)
                soj = basis%ao_offset(jsh)
                pure_j = 0
                if (HARMONIC_ACTIVE) pure_j = basis%harmonic(jsh)
                call c2s_expand_block(dfull(soi:soi + nsi - 1, soj:soj + nsj - 1), &
                                      dcart(coi:coi + nci - 1, coj:coj + ncj - 1), &
                                      basis%am(ish), pure_i, basis%am(jsh), pure_j)
            end do
        end do
        deallocate(dfull)

        allocate(label_map(nbf_raw))
        call raw_cart_map(basis, cart_off, label_map)
        allocate(draw(nbf_raw, nbf_raw))
        do j = 1, nbf_raw
            do i = 1, nbf_raw
                draw(i, j) = dcart(label_map(i), label_map(j))
            end do
        end do

    end subroutine ecp_density_to_raw

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




!###############################################################################
! Indexing helpers
!###############################################################################

  !> Number of Cartesian tuples with total degree <= n
  pure integer function ntup_upto(n)
    integer, intent(in) :: n
    if (n < 0) then
      ntup_upto = 0
    else
      ntup_upto = (n + 1)*(n + 2)*(n + 3)/6
    end if
  end function ntup_upto

  !> 1-based global index of the tuple (i,j,k); within a degree the order is
  !> x power descending, then y power descending
  pure integer function tid(i, j, k)
    integer, intent(in) :: i, j, k
    integer :: n
    n = i + j + k
    tid = n*(n + 1)*(n + 2)/6 + (n - i)*(n - i + 1)/2 + (n - i - j) + 1
  end function tid

  !> The tuple with global index t
  pure subroutine tuple_of(t, i, j, k)
    integer, intent(in) :: t
    integer, intent(out) :: i, j, k
    integer :: n, loc
    n = 0
    do while (ntup_upto(n) < t)
      n = n + 1
    end do
    loc = t - ntup_upto(n - 1) - 1
    i = n
    do while (loc > n - i)
      loc = loc - (n - i + 1)
      i = i - 1
    end do
    j = n - i - loc
    k = n - i - j
  end subroutine tuple_of

  !> Index of (l,m) in an l-major list of real spherical harmonics, 1-based
  pure integer function lmi(l, m)
    integer, intent(in) :: l, m
    lmi = l*l + l + m + 1
  end function lmi

  pure real(dp) function ipow(x, n)
    real(dp), intent(in) :: x
    integer, intent(in) :: n
    if (n == 0) then
      ipow = 1.0_dp
    else
      ipow = x**n
    end if
  end function ipow

  pure real(dp) function binom(n, k)
    integer, intent(in) :: n, k
    integer :: i
    binom = 1.0_dp
    if (k < 0 .or. k > n) then
      binom = 0.0_dp
      return
    end if
    do i = 1, k
      binom = binom*real(n - k + i, dp)/real(i, dp)
    end do
  end function binom

  pure real(dp) function fact(n)
    integer, intent(in) :: n
    integer :: i
    fact = 1.0_dp
    do i = 2, n
      fact = fact*real(i, dp)
    end do
  end function fact

  pure real(dp) function dfact(n)
    integer, intent(in) :: n
    integer :: i
    dfact = 1.0_dp
    i = n
    do while (i > 1)
      dfact = dfact*real(i, dp)
      i = i - 2
    end do
  end function dfact

  !> 0-based index of the first matrix of the (I,J) atom block, I <= J (0-based atoms)
  pure integer function ecp_hess_start(i0, j0, natm)
    integer, intent(in) :: i0, j0, natm
    ecp_hess_start = 9*j0 + 3*(3*natm - 1)*i0 - (9*i0*(i0 + 1))/2 - 3
  end function ecp_hess_start

!###############################################################################
! Tables
!###############################################################################

  !> Real orthonormal spherical harmonics as homogeneous polynomials in the
  !> direction cosines (Helgaker, Joergensen & Olsen, eq. 6.4.47, scaled by
  !> sqrt((2l+1)/4pi))
  subroutine build_ylm(tab, lmax)
    type(tables_t), intent(inout) :: tab
    integer, intent(in) :: lmax
    integer :: l, m, am, t, u, v2, v2min, px, py, pz, sgn
    real(dp) :: nrm, c

    allocate(tab%ycoef(ntup_upto(lmax), (lmax + 1)**2), source=0.0_dp)
    do l = 0, lmax
      do m = -l, l
        am = abs(m)
        nrm = 1.0_dp/(2.0_dp**am*fact(l)) &
            * sqrt(2.0_dp*fact(l + am)*fact(l - am)/merge(2.0_dp, 1.0_dp, m == 0))
        v2min = merge(0, 1, m >= 0)
        do t = 0, (l - am)/2
          do u = 0, t
            do v2 = v2min, am, 2
              sgn = t + (v2 - v2min)/2
              c = merge(-1.0_dp, 1.0_dp, mod(sgn, 2) == 1) * 0.25_dp**t &
                * binom(l, t)*binom(l - t, am + t)*binom(t, u)*binom(am, v2)
              px = 2*t + am - 2*u - v2
              py = 2*u + v2
              pz = l - 2*t - am
              tab%ycoef(tid(px, py, pz), lmi(l, m)) = &
                  tab%ycoef(tid(px, py, pz), lmi(l, m)) + nrm*c
            end do
          end do
        end do
        tab%ycoef(:, lmi(l, m)) = tab%ycoef(:, lmi(l, m))*sqrt(real(2*l + 1, dp)/FOURPI)
      end do
    end do
  end subroutine build_ylm

  !> Integrals over the unit sphere of x^p y^q z^r
  subroutine build_moments(tab, pmax)
    type(tables_t), intent(inout) :: tab
    integer, intent(in) :: pmax
    integer :: p, q, r
    allocate(tab%mom(0:pmax, 0:pmax, 0:pmax), source=0.0_dp)
    do r = 0, pmax, 2
      do q = 0, pmax, 2
        do p = 0, pmax, 2
          tab%mom(p, q, r) = FOURPI*dfact(p - 1)*dfact(q - 1)*dfact(r - 1)/dfact(p + q + r + 1)
        end do
      end do
    end do
  end subroutine build_moments

  !> Gauss-Legendre abscissae and weights on [-1,1]
  subroutine build_gl(tab)
    type(tables_t), intent(inout) :: tab
    integer :: i, j, it
    real(dp) :: z, z1, p1, p2, p3, pp
    do i = 1, NGL/2 + mod(NGL, 2)
      z = cos(PI*(real(i, dp) - 0.25_dp)/(real(NGL, dp) + 0.5_dp))
      do it = 1, 100
        p1 = 1.0_dp
        p2 = 0.0_dp
        do j = 1, NGL
          p3 = p2
          p2 = p1
          p1 = (real(2*j - 1, dp)*z*p2 - real(j - 1, dp)*p3)/real(j, dp)
        end do
        pp = real(NGL, dp)*(z*p1 - p2)/(z*z - 1.0_dp)
        z1 = z
        z = z1 - p1/pp
        if (abs(z - z1) < 1.0e-15_dp) exit
      end do
      tab%glx(i) = -z
      tab%glx(NGL + 1 - i) = z
      tab%glw(i) = 2.0_dp/((1.0_dp - z*z)*pp*pp)
      tab%glw(NGL + 1 - i) = tab%glw(i)
    end do
  end subroutine build_gl

  !> Angular tables.
  !>   om1(s, lam mu)    = int x^s Y_lam,mu dOmega,                deg s <= 2 lb
  !>   om2(s, lm, lam mu) = int x^s Y_lm Y_lam,mu dOmega, l < le,   deg s <= lb
  subroutine build_tables(tab, lb, le)
    type(tables_t), intent(inout) :: tab
    integer, intent(in) :: lb, le
    integer :: lammax2, ly, nq, s, t1, t2, i, j, k, a, b, c, l, m, lam, mu, lm, lmu, n1
    real(dp), allocatable :: q(:,:)

    tab%lb = lb
    tab%le = le
    lammax2 = max(le - 1, 0) + lb
    ly = max(2*lb, lammax2, le)
    tab%ly = ly
    tab%pm = max(4*lb, 2*lb + 2*le) + 2
    call build_ylm(tab, ly)
    call build_moments(tab, tab%pm)
    allocate(tab%tx(ntup_upto(tab%pm)), tab%ty(ntup_upto(tab%pm)), tab%tz(ntup_upto(tab%pm)), &
             tab%tn(ntup_upto(tab%pm)))
    do s = 1, ntup_upto(tab%pm)
      call tuple_of(s, tab%tx(s), tab%ty(s), tab%tz(s))
      tab%tn(s) = tab%tx(s) + tab%ty(s) + tab%tz(s)
    end do
    call build_gl(tab)

    ! om1
    allocate(tab%om1(ntup_upto(2*lb), (2*lb + 1)**2), source=0.0_dp)
    do s = 1, ntup_upto(2*lb)
      call tuple_of(s, i, j, k)
      do lam = 0, i + j + k
        do mu = -lam, lam
          lmu = lmi(lam, mu)
          do t2 = ntup_upto(lam - 1) + 1, ntup_upto(lam)
            if (tab%ycoef(t2, lmu) == 0.0_dp) cycle
            call tuple_of(t2, a, b, c)
            tab%om1(s, lmu) = tab%om1(s, lmu) + tab%ycoef(t2, lmu)*tab%mom(i + a, j + b, k + c)
          end do
        end do
      end do
    end do

    if (le < 1) return
    ! q(s, lm) = int x^s Y_lm, deg s <= lb + lammax2
    nq = ntup_upto(lb + lammax2)
    allocate(q(nq, le*le), source=0.0_dp)
    do s = 1, nq
      call tuple_of(s, i, j, k)
      do l = 0, le - 1
        do m = -l, l
          lm = lmi(l, m)
          do t1 = ntup_upto(l - 1) + 1, ntup_upto(l)
            if (tab%ycoef(t1, lm) == 0.0_dp) cycle
            call tuple_of(t1, a, b, c)
            q(s, lm) = q(s, lm) + tab%ycoef(t1, lm)*tab%mom(i + a, j + b, k + c)
          end do
        end do
      end do
    end do
    allocate(tab%om2(ntup_upto(lb), le*le, (lammax2 + 1)**2), source=0.0_dp)
    do s = 1, ntup_upto(lb)
      call tuple_of(s, i, j, k)
      n1 = i + j + k
      do lm = 1, le*le
        do lam = 0, lammax2
          do mu = -lam, lam
            lmu = lmi(lam, mu)
            do t2 = ntup_upto(lam - 1) + 1, ntup_upto(lam)
              if (tab%ycoef(t2, lmu) == 0.0_dp) cycle
              call tuple_of(t2, a, b, c)
              tab%om2(s, lm, lmu) = tab%om2(s, lm, lmu) &
                  + tab%ycoef(t2, lmu)*q(tid(i + a, j + b, k + c), lm)
            end do
          end do
        end do
      end do
    end do
  end subroutine build_tables

  !> Values of Y_lm(d) for l <= lmax at the unit vector d
  subroutine ylm_at(tab, lmax, d, y)
    type(tables_t), intent(in) :: tab
    integer, intent(in) :: lmax
    real(dp), intent(in) :: d(3)
    real(dp), intent(out) :: y(:)
    integer :: t, i, j, k, l, lm
    real(dp) :: mono
    y(1:(lmax + 1)**2) = 0.0_dp
    do l = 0, lmax
      do t = ntup_upto(l - 1) + 1, ntup_upto(l)
        i = tab%tx(t)
        j = tab%ty(t)
        k = tab%tz(t)
        mono = ipow(d(1), i)*ipow(d(2), j)*ipow(d(3), k)
        do lm = l*l + 1, (l + 1)**2
          y(lm) = y(lm) + tab%ycoef(t, lm)*mono
        end do
      end do
    end do
  end subroutine ylm_at

!###############################################################################
! Scaled modified spherical Bessel functions  K_l(x) = exp(-x) i_l(x)
!###############################################################################

  pure subroutine kbessel(lmax, x, kb)
    integer, intent(in) :: lmax
    real(dp), intent(in) :: x
    real(dp), intent(out) :: kb(0:lmax)
    integer :: l, l0, j
    real(dp) :: e2, k0, ex, term, s
    ! fixed size: an automatic array here would be heap-allocated on every
    ! call (once per quadrature point); l0 + 1 <= lmax + 17 + max(16, 2 lmax)
    integer, parameter :: KB_LMAX = 40
    real(dp) :: f(0:3*KB_LMAX + 34)

    if (lmax > KB_LMAX) error stop "ecp: kbessel order above KB_LMAX"
    if (x <= 0.0_dp) then
      kb = 0.0_dp
      kb(0) = 1.0_dp
      return
    end if
    if (x < 2.0_dp) then
      ! power series for the two highest orders, then the downward recurrence
      ! i_(l-1) = i_(l+1) + (2l+1)/x i_l
      ex = exp(-x)
      do l = lmax, max(lmax - 1, 0), -1
        term = ex
        do j = 1, l
          term = term*x/real(2*j + 1, dp)
        end do
        s = term
        e2 = term
        do j = 1, 100
          e2 = e2*(0.5_dp*x*x)/(real(j, dp)*real(2*l + 2*j + 1, dp))
          s = s + e2
          if (e2 <= 1.0e-17_dp*s) exit
        end do
        kb(l) = s
      end do
      do l = lmax - 1, 1, -1
        kb(l - 1) = kb(l + 1) + real(2*l + 1, dp)/x*kb(l)
      end do
    else if (x < max(16.0_dp, 2.0_dp*real(lmax, dp))) then
      ! Miller: the same downward recurrence from an arbitrary start above lmax,
      ! normalised to the exact K_0
      k0 = (1.0_dp - exp(-2.0_dp*x))/(2.0_dp*x)
      l0 = lmax + 16 + int(x)
      f(l0 + 1) = 0.0_dp
      f(l0) = 1.0e-280_dp
      do l = l0, 1, -1
        f(l - 1) = f(l + 1) + real(2*l + 1, dp)/x*f(l)
        if (abs(f(l - 1)) > 1.0e250_dp) f(l - 1:l0) = f(l - 1:l0)*1.0e-250_dp
      end do
      kb(0:lmax) = f(0:lmax)*(k0/f(0))
    else
      e2 = exp(-2.0_dp*x)
      kb(0) = (1.0_dp - e2)/(2.0_dp*x)
      if (lmax >= 1) kb(1) = (0.5_dp*(1.0_dp + e2) - kb(0))/x
      do l = 1, lmax - 1
        kb(l + 1) = kb(l - 1) - real(2*l + 1, dp)/x*kb(l)
      end do
    end if
  end subroutine kbessel

!###############################################################################
! Main driver
!###############################################################################

  !> @brief Raw Cartesian ECP integrals and nuclear derivatives (layout in the module header).
  !> @param[in]  basis  Basis set with ECP parameters.
  !> @param[in]  coord  Nuclear coordinates (3 x natm), Bohr.
  !> @param[in]  deriv  0, 1 or 2.
  !> @param[out] res    Concatenated nraw x nraw matrices (see module header).
  subroutine ecp_raw_ints(basis, coord, deriv, res)
    type(basis_set), intent(in) :: basis
    real(dp), intent(in) :: coord(:,:)
    integer, intent(in) :: deriv
    real(dp), allocatable, intent(out) :: res(:)

    integer, allocatable :: off(:)
    integer :: natm, nraw, nmat

    natm = size(coord, 2)
    call raw_offsets(basis, off, nraw)
    select case (deriv)
    case (0)
      nmat = 1
    case (1)
      nmat = 3*natm
    case default
      nmat = (3*natm*(3*natm + 1))/2
    end select
    allocate(res(int(nraw, 8)*int(nraw, 8)*nmat), source=0.0_dp)
    if (.not. basis%ecp_params%is_ecp) return

    call ecp_centres(basis, coord, deriv, off, nraw, res)
  end subroutine ecp_raw_ints

  !> @brief Second nuclear derivatives of the ECP energy, contracted per shell pair.
  !> @detail hess(3(I-1)+a, 3(J-1)+b) += sum_kl dens(k,l) d^2 V_kl / dR_Ia dR_Jb,
  !>         with dens in the raw Cartesian order of @ref ecp_raw_ints.  The
  !>         derivative matrices are never stored, so the memory is that of dens
  !>         and one 3N x 3N accumulator per thread.
  !> @param[in]    dens  Symmetric density, nraw x nraw, raw Cartesian order.
  !> @param[inout] hess  Cartesian Hessian (3*natm x 3*natm), incremented.
  subroutine ecp_hess_contract(basis, coord, dens, hess)
    type(basis_set), intent(in) :: basis
    real(dp), intent(in) :: coord(:,:)
    real(dp), intent(in) :: dens(:,:)
    real(dp), intent(inout) :: hess(:,:)

    integer, allocatable :: off(:)
    integer :: nraw
    real(dp) :: dummy(1)

    if (.not. basis%ecp_params%is_ecp) return
    call raw_offsets(basis, off, nraw)
    call ecp_centres(basis, coord, 2, off, nraw, dummy, dens, hess)
  end subroutine ecp_hess_contract

  !> Offsets of the shells in the raw Cartesian order and its dimension
  subroutine raw_offsets(basis, off, nraw)
    type(basis_set), intent(in) :: basis
    integer, allocatable, intent(out) :: off(:)
    integer, intent(out) :: nraw
    integer :: ish

    allocate(off(basis%nshell))
    nraw = 0
    do ish = 1, basis%nshell
      off(ish) = nraw
      nraw = nraw + NUM_CART_BF(basis%am(ish))
    end do
  end subroutine raw_offsets

  !> Loop over ECP centres and shell pairs.  Without dens/hess the derivative
  !> matrices are added to res; with them (deriv = 2 only) they are contracted
  !> into hess instead and res is not referenced.
  subroutine ecp_centres(basis, coord, deriv, off, nraw, res, dens, hess)
    type(basis_set), intent(in) :: basis
    real(dp), intent(in) :: coord(:,:)
    integer, intent(in) :: deriv, nraw
    integer, intent(in) :: off(:)
    real(dp), intent(inout) :: res(:)
    real(dp), intent(in), optional :: dens(:,:)
    real(dp), intent(inout), optional :: hess(:,:)

    type(tables_t) :: tab
    type(gtab_t), allocatable :: gt(:)
    integer, allocatable :: toff(:)
    integer :: natm, nsh, lbmax, le, lecp, ic, nc, it, iat, iatc, s1, s2
    real(dp) :: c(3)
    real(dp), allocatable :: hloc(:,:)
    logical :: contract

    natm = size(coord, 2)
    nsh = basis%nshell
    contract = present(hess)
    if (contract) then
      allocate(hloc(3*natm, 3*natm), source=0.0_dp)
    else
      allocate(hloc(1, 1), source=0.0_dp)
    end if

    nc = size(basis%ecp_params%n_expo)
    allocate(toff(nc + 1))
    toff(1) = 0
    do ic = 1, nc
      toff(ic + 1) = toff(ic) + basis%ecp_params%n_expo(ic)
    end do
    le = 0
    do it = 1, toff(nc + 1)
      le = max(le, basis%ecp_params%ecp_am(it))
    end do
    lbmax = maxval(basis%am(1:nsh)) + deriv
    call build_tables(tab, lbmax, le)

    allocate(gt(natm))
    do ic = 1, nc
      c = basis%ecp_params%ecp_coord(3*ic - 2:3*ic)
      iatc = 0
      do iat = 1, natm
        if (sum(abs(coord(:, iat) - c)) < 1.0e-4_dp) then
          iatc = iat
          exit
        end if
      end do
      ! derivatives must be attributed to the atom that carries the ECP
      if (iatc == 0 .and. deriv > 0) call show_message( &
          'ECP centre does not coincide with any atom; cannot assign its nuclear derivatives', &
          WITH_ABORT)
      lecp = 0
      do it = toff(ic) + 1, toff(ic + 1)
        lecp = max(lecp, basis%ecp_params%ecp_am(it))
      end do

      call build_gtabs(tab, basis, coord, c, lecp, deriv, gt)

      if (contract) then
        !$omp parallel do schedule(dynamic) private(s1, s2) reduction(+:hloc)
        do s1 = 1, nsh
          do s2 = 1, s1
            call shell_pair(tab, basis, coord, c, iatc, toff(ic) + 1, toff(ic + 1), lecp, &
                            deriv, gt, s1, s2, off, nraw, natm, res, dens, hloc)
          end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do schedule(dynamic) private(s1, s2)
        do s1 = 1, nsh
          do s2 = 1, s1
            call shell_pair(tab, basis, coord, c, iatc, toff(ic) + 1, toff(ic + 1), lecp, &
                            deriv, gt, s1, s2, off, nraw, natm, res)
          end do
        end do
        !$omp end parallel do
      end if

      do iat = 1, natm
        if (allocated(gt(iat)%g)) deallocate(gt(iat)%g)
        gt(iat)%used = .false.
      end do
    end do

    if (contract) hess(1:3*natm, 1:3*natm) = hess(1:3*natm, 1:3*natm) + hloc
  end subroutine ecp_centres

  !> Type-2 projection tables of every atom that carries shells, for one ECP
  !>   g(t, lm, N, lam) = sum_{s <= t, deg s = N} c_t(s) 4pi sum_mu Y_lam,mu(A^) om2(s, lm, lam mu)
  !> where c_t(s) expands (r-A)^t in monomials r^s about the ECP centre.
  subroutine build_gtabs(tab, basis, coord, c, lecp, deriv, gt)
    type(tables_t), intent(in) :: tab
    type(basis_set), intent(in) :: basis
    real(dp), intent(in) :: coord(:,:), c(3)
    integer, intent(in) :: lecp, deriv
    type(gtab_t), intent(inout) :: gt(:)
    integer :: iat, ish, na, lammax, t, i, j, k, u, v, w, s, lm, lam, mu, nn
    real(dp) :: a(3), an, ahat(3), cu, cv, cw
    real(dp), allocatable :: y(:), h(:,:,:)

    if (lecp < 1) return
    do iat = 1, size(gt)
      gt(iat)%namax = -1
    end do
    do ish = 1, basis%nshell
      iat = basis%origin(ish)
      gt(iat)%namax = max(gt(iat)%namax, basis%am(ish) + deriv)
    end do
    allocate(y((tab%ly + 1)**2))
    do iat = 1, size(gt)
      na = gt(iat)%namax
      if (na < 0) cycle
      a = coord(:, iat) - c
      an = sqrt(sum(a*a))
      gt(iat)%oncentre = sum(abs(a)) < ONCENTRE
      if (gt(iat)%oncentre) then
        lammax = 0
        ahat = [0.0_dp, 0.0_dp, 1.0_dp]
      else
        lammax = lecp - 1 + na
        ahat = a/an
      end if
      gt(iat)%lammax = lammax
      call ylm_at(tab, lammax, ahat, y)
      ! h(s, lm, lam)
      allocate(h(ntup_upto(na), lecp*lecp, 0:lammax), source=0.0_dp)
      do s = 1, ntup_upto(na)
        do lm = 1, lecp*lecp
          do lam = 0, lammax
            do mu = -lam, lam
              h(s, lm, lam) = h(s, lm, lam) + FOURPI*y(lmi(lam, mu))*tab%om2(s, lm, lmi(lam, mu))
            end do
          end do
        end do
      end do
      allocate(gt(iat)%g(ntup_upto(na), lecp*lecp, 0:na, 0:lammax), source=0.0_dp)
      do t = 1, ntup_upto(na)
        call tuple_of(t, i, j, k)
        do u = 0, i
          cu = binom(i, u)*ipow(-a(1), i - u)
          if (cu == 0.0_dp) cycle
          do v = 0, j
            cv = cu*binom(j, v)*ipow(-a(2), j - v)
            if (cv == 0.0_dp) cycle
            do w = 0, k
              cw = cv*binom(k, w)*ipow(-a(3), k - w)
              if (cw == 0.0_dp) cycle
              s = tid(u, v, w)
              nn = u + v + w
              gt(iat)%g(t, :, nn, :) = gt(iat)%g(t, :, nn, :) + cw*h(s, :, :)
            end do
          end do
        end do
      end do
      deallocate(h)
      gt(iat)%used = .true.
    end do
  end subroutine build_gtabs

  !> Expansion of (r-A)^t about the ECP centre:  sum_s cs(s) r^s
  subroutine expand_about(t, a, ns, sidx, scoef)
    integer, intent(in) :: t
    real(dp), intent(in) :: a(3)
    integer, intent(out) :: ns
    integer, intent(out) :: sidx(:)
    real(dp), intent(out) :: scoef(:)
    integer :: i, j, k, u, v, w
    real(dp) :: cu, cv, cw
    call tuple_of(t, i, j, k)
    ns = 0
    do u = 0, i
      cu = binom(i, u)*ipow(-a(1), i - u)
      if (cu == 0.0_dp) cycle
      do v = 0, j
        cv = cu*binom(j, v)*ipow(-a(2), j - v)
        if (cv == 0.0_dp) cycle
        do w = 0, k
          cw = cv*binom(k, w)*ipow(-a(3), k - w)
          if (cw == 0.0_dp) cycle
          ns = ns + 1
          sidx(ns) = tid(u, v, w)
          scoef(ns) = cw
        end do
      end do
    end do
  end subroutine expand_about

  !> Type-1 angular factors  t(s, lam) = 4 pi sum_mu Y_lam,mu(d) om1(s, lam mu)
  !> for every tuple s with deg s <= nmax (zero unless lam <= deg s, same parity)
  subroutine angular_type1(tab, nmax, d, y, t)
    type(tables_t), intent(in) :: tab
    integer, intent(in) :: nmax
    real(dp), intent(in) :: d(3)
    real(dp), intent(inout) :: y(:)
    real(dp), intent(out) :: t(:, 0:)
    integer :: s, n, lam, mu
    real(dp) :: e
    call ylm_at(tab, nmax, d, y)
    t = 0.0_dp
    do s = 1, ntup_upto(nmax)
      n = tab%tn(s)
      do lam = mod(n, 2), n, 2
        e = 0.0_dp
        do mu = -lam, lam
          e = e + y(lmi(lam, mu))*tab%om1(s, lmi(lam, mu))
        end do
        t(s, lam) = FOURPI*e
      end do
    end do
  end subroutine angular_type1

  !> One shell pair with one ECP centre: sizes the work arrays for
  !> shell_pair_core, which evaluates it
  subroutine shell_pair(tab, basis, coord, c, iatc, t1, t2, lecp, deriv, gt, s1, s2, off, nraw, natm, res, &
                        dens, hess)
    type(tables_t), intent(in) :: tab
    type(basis_set), intent(in) :: basis
    real(dp), intent(in) :: coord(:,:), c(3)
    integer, intent(in) :: iatc, t1, t2, lecp, deriv, s1, s2, nraw, natm
    type(gtab_t), intent(in) :: gt(:)
    integer, intent(in) :: off(:)
    real(dp), intent(inout) :: res(:)
    real(dp), intent(in), optional :: dens(:,:)
    real(dp), intent(inout), optional :: hess(:,:)
    integer :: na, nb, lamA, lamB

    na = basis%am(s1) + deriv
    nb = basis%am(s2) + deriv
    lamA = 0
    lamB = 0
    if (lecp >= 1) then
      if (sum(abs(coord(:, basis%origin(s1)) - c)) >= ONCENTRE) lamA = lecp - 1 + na
      if (sum(abs(coord(:, basis%origin(s2)) - c)) >= ONCENTRE) lamB = lecp - 1 + nb
    end if
    call shell_pair_core(tab, basis, coord, c, iatc, t1, t2, lecp, deriv, gt, s1, s2, off, nraw, natm, res, &
                         na, nb, na + nb, ntup_upto(na), ntup_upto(nb), merge(1, merge(3, 6, deriv == 1), deriv == 0), &
                         ntup_upto(max(na, nb)), lamA, lamB, dens, hess)
  end subroutine shell_pair

  !> One shell pair with one ECP centre: weighted primitive sums, then the value
  !> and derivative blocks, scattered into res
  subroutine shell_pair_core(tab, basis, coord, c, iatc, t1, t2, lecp, deriv, gt, s1, s2, off, nraw, natm, res, &
                             na, nb, nab, ntA, ntB, nw, nsmax, lamA, lamB, dens, hess)
    type(tables_t), intent(in) :: tab
    type(basis_set), intent(in) :: basis
    real(dp), intent(in) :: coord(:,:), c(3)
    integer, intent(in) :: iatc, t1, t2, lecp, deriv, s1, s2, nraw, natm
    type(gtab_t), intent(in) :: gt(:)
    integer, intent(in) :: off(:)
    real(dp), intent(inout) :: res(:)
    integer, intent(in) :: na, nb, nab, ntA, ntB, nw, nsmax, lamA, lamB
    real(dp), intent(in), optional :: dens(:,:)
    real(dp), intent(inout), optional :: hess(:,:)

    integer :: la, lb, iata, iatb, ipa, ipb, it, l, n
    integer :: ta, tb, ia, ib, lm, mm, nA_, nB_, lam1, lam2
    real(dp) :: a(3), b(3), an, bn, ea, eb, ca, cb, wts(6), pv(3), pn, p, pref, phat(3)
    real(dp) :: logamp, q, r0, dw, lo, hi, rr, g, zt, dc, prod
    logical :: onA, onB, havelocal
    integer :: nterm
    logical :: fixdir
    integer :: ig
    ! work arrays sized by the shell pair (automatic: no heap traffic per pair)
    real(dp) :: iw(nw, ntA, ntB)
    integer :: sidxa(nsmax, ntA), sidxb(nsmax, ntB), nsa(ntA), nsb(ntB)
    real(dp) :: scoa(nsmax, ntA), scob(nsmax, ntB)
    real(dp) :: sw(nw, ntup_upto(nab)), s(ntup_upto(nab)), r1(0:nab, 0:nab), y((nab + 1)**2)
    real(dp) :: tang(ntup_upto(nab), 0:nab), kb1(0:nab)
    real(dp) :: r2(0:nab, 0:lamA, 0:lamB), r2w(nw, 0:nab, 0:lamA, 0:lamB), x(ntA, 0:nb, 0:lamB)
    real(dp) :: kba(0:lamA), kbb(0:lamB)

    la = basis%am(s1)
    lb = basis%am(s2)
    iata = basis%origin(s1)
    iatb = basis%origin(s2)
    a = coord(:, iata) - c
    b = coord(:, iatb) - c
    an = sqrt(sum(a*a))
    bn = sqrt(sum(b*b))
    onA = sum(abs(a)) < ONCENTRE
    onB = sum(abs(b)) < ONCENTRE
    iw = 0.0_dp

    ! expansions of (r-A)^t and (r-B)^t about the ECP centre
    do ta = 1, ntA
      call expand_about(ta, a, nsa(ta), sidxa(:, ta), scoa(:, ta))
    end do
    do tb = 1, ntB
      call expand_about(tb, b, nsb(tb), sidxb(:, tb), scob(:, tb))
    end do

    !------------------------------------------------------------------ type 1
    havelocal = .false.
    do it = t1, t2
      if (basis%ecp_params%ecp_am(it) == lecp) havelocal = .true.
    end do
    if (havelocal) then
      sw = 0.0_dp
      fixdir = iata == iatb .or. onA .or. onB
      if (fixdir) then
        if (.not. onA) then
          phat = a/an
        else if (.not. onB) then
          phat = b/bn
        else
          phat = [0.0_dp, 0.0_dp, 1.0_dp]
        end if
        call angular_type1(tab, nab, phat, y, tang)
      end if
      do ipa = 0, basis%ncontr(s1) - 1
        ea = basis%ex(basis%g_offset(s1) + ipa)
        ca = basis%cc(basis%g_offset(s1) + ipa)
        do ipb = 0, basis%ncontr(s2) - 1
          eb = basis%ex(basis%g_offset(s2) + ipb)
          cb = basis%cc(basis%g_offset(s2) + ipb)
          p = ea + eb
          pv = (ea*a + eb*b)/p
          pn = sqrt(sum(pv*pv))
          pref = -ea*an*an - eb*bn*bn + p*pn*pn
          r1 = 0.0_dp
          nterm = 0
          do it = t1, t2
            if (basis%ecp_params%ecp_am(it) /= lecp) cycle
            zt = basis%ecp_params%ecp_ex(it)
            dc = basis%ecp_params%ecp_cc(it)
            n = basis%ecp_params%ecp_r_ex(it)
            q = p + zt
            logamp = pref - p*zt*pn*pn/q
            if (logamp < LOGSKIP) cycle
            nterm = nterm + 1
            r0 = p*pn/q
            dw = sqrt(WIN2/q)
            lo = max(0.0_dp, r0 - dw)
            hi = r0 + dw
            do ig = 1, NGL
              rr = 0.5_dp*(hi - lo)*tab%glx(ig) + 0.5_dp*(hi + lo)
              g = dc*0.5_dp*(hi - lo)*tab%glw(ig) &
                * exp(-zt*rr*rr - p*(rr - pn)**2 + pref)*ipow(rr, n)
              call kbessel(nab, 2.0_dp*p*pn*rr, kb1)
              do l = 0, nab
                do lam1 = mod(l, 2), l, 2
                  r1(l, lam1) = r1(l, lam1) + g*kb1(lam1)
                end do
                g = g*rr
              end do
            end do
          end do
          if (nterm == 0) cycle
          ! angular factors at P^; the direction is fixed for the whole shell
          ! pair when both shells share an atom or one sits on the ECP centre
          if (.not. fixdir) then
            if (pn > 1.0e-12_dp) then
              phat = pv/pn
            else
              phat = [0.0_dp, 0.0_dp, 1.0_dp]
            end if
            call angular_type1(tab, nab, phat, y, tang)
          end if
          do ta = 1, ntup_upto(nab)
            n = tab%tn(ta)
            prod = 0.0_dp
            do lam1 = mod(n, 2), n, 2
              prod = prod + tang(ta, lam1)*r1(n, lam1)
            end do
            s(ta) = prod
          end do
          call weights(nw, ea, eb, wts)
          do ia = 1, nw
            sw(ia, :) = sw(ia, :) + wts(ia)*ca*cb*s(:)
          end do
        end do
      end do
      do tb = 1, ntB
        do ta = 1, ntA
          do ia = 1, nsa(ta)
            do ib = 1, nsb(tb)
              iw(:, ta, tb) = iw(:, ta, tb) + scoa(ia, ta)*scob(ib, tb) &
                  * sw(:, tid(tab%tx(sidxa(ia, ta)) + tab%tx(sidxb(ib, tb)), &
                              tab%ty(sidxa(ia, ta)) + tab%ty(sidxb(ib, tb)), &
                              tab%tz(sidxa(ia, ta)) + tab%tz(sidxb(ib, tb))))
            end do
          end do
        end do
      end do
    end if

    !------------------------------------------------------------------ type 2
    if (lecp >= 1 .and. gt(iata)%used .and. gt(iatb)%used) then
      do l = 0, lecp - 1
        r2w = 0.0_dp
        nterm = 0
        do ipa = 0, basis%ncontr(s1) - 1
          ea = basis%ex(basis%g_offset(s1) + ipa)
          ca = basis%cc(basis%g_offset(s1) + ipa)
          do ipb = 0, basis%ncontr(s2) - 1
            eb = basis%ex(basis%g_offset(s2) + ipb)
            cb = basis%cc(basis%g_offset(s2) + ipb)
            r2 = 0.0_dp
            do it = t1, t2
              if (basis%ecp_params%ecp_am(it) /= l) cycle
              zt = basis%ecp_params%ecp_ex(it)
              dc = basis%ecp_params%ecp_cc(it)
              n = basis%ecp_params%ecp_r_ex(it)
              q = ea + eb + zt
              logamp = -(ea*an*an + eb*bn*bn) + (ea*an + eb*bn)**2/q
              if (logamp < LOGSKIP) cycle
              nterm = nterm + 1
              r0 = (ea*an + eb*bn)/q
              dw = sqrt(WIN2/q)
              lo = max(0.0_dp, r0 - dw)
              hi = r0 + dw
              do ig = 1, NGL
                rr = 0.5_dp*(hi - lo)*tab%glx(ig) + 0.5_dp*(hi + lo)
                g = dc*0.5_dp*(hi - lo)*tab%glw(ig) &
                  * exp(-zt*rr*rr - ea*(rr - an)**2 - eb*(rr - bn)**2)*ipow(rr, n)
                call kbessel(lamA, 2.0_dp*ea*an*rr, kba)
                call kbessel(lamB, 2.0_dp*eb*bn*rr, kbb)
                do mm = 0, nab
                  do lam2 = 0, lamB
                    do lam1 = 0, lamA
                      r2(mm, lam1, lam2) = r2(mm, lam1, lam2) + g*kba(lam1)*kbb(lam2)
                    end do
                  end do
                  g = g*rr
                end do
              end do
            end do
            call weights(nw, ea, eb, wts)
            do ia = 1, nw
              r2w(ia, :, :, :) = r2w(ia, :, :, :) + wts(ia)*ca*cb*r2
            end do
          end do
        end do
        if (nterm == 0) cycle
        do ia = 1, nw
          do lm = l*l + 1, (l + 1)**2
            ! x(ta, Nb, lam2) = sum_{Na, lam1} gA(ta, lm, Na, lam1) r2w(Na+Nb, lam1, lam2)
            x = 0.0_dp
            do lam2 = 0, lamB
              do nB_ = 0, nb
                do lam1 = 0, lamA
                  do nA_ = 0, na
                    if (r2w(ia, nA_ + nB_, lam1, lam2) == 0.0_dp) cycle
                    x(:, nB_, lam2) = x(:, nB_, lam2) &
                        + gt(iata)%g(1:ntA, lm, nA_, lam1)*r2w(ia, nA_ + nB_, lam1, lam2)
                  end do
                end do
              end do
            end do
            do tb = 1, ntB
              do lam2 = 0, lamB
                do nB_ = 0, nb
                  if (gt(iatb)%g(tb, lm, nB_, lam2) == 0.0_dp) cycle
                  iw(ia, :, tb) = iw(ia, :, tb) + x(:, nB_, lam2)*gt(iatb)%g(tb, lm, nB_, lam2)
                end do
              end do
            end do
          end do
        end do
      end do
    end if

    call scatter(iw, nw, la, lb, deriv, iata, iatb, iatc, s1, s2, off, nraw, natm, res, dens, hess)

  end subroutine shell_pair_core

  !> Primitive weights: 1, a, b, a^2, ab, b^2
  pure subroutine weights(nw, ea, eb, w)
    integer, intent(in) :: nw
    real(dp), intent(in) :: ea, eb
    real(dp), intent(out) :: w(6)
    w = [1.0_dp, ea, eb, ea*ea, ea*eb, eb*eb]
    if (nw < 6) w(nw + 1:) = 0.0_dp
  end subroutine weights

  !> Value and derivative blocks from the weighted tuple integrals, added to res;
  !> with dens and hess (deriv = 2), the second-derivative elements are instead
  !> contracted with dens and added to hess
  subroutine scatter(iw, nw, la, lb, deriv, iata, iatb, iatc, s1, s2, off, nraw, natm, res, dens, hess)
    integer, intent(in) :: nw, la, lb, deriv, iata, iatb, iatc, s1, s2, nraw, natm
    real(dp), intent(in) :: iw(:,:,:)
    integer, intent(in) :: off(:)
    real(dp), intent(inout) :: res(:)
    real(dp), intent(in), optional :: dens(:,:)
    real(dp), intent(inout), optional :: hess(:,:)

    integer, parameter :: W1 = 1, WA = 2, WB = 3, WAA = 4, WAB = 5, WBB = 6
    integer :: ca_, cb_, ta, tb, ka(3), kb(3), k, k2, row, col, nca, ncb
    integer :: x, y, ix, iy, kk, gi, gj, iat, jat, mat, n0
    integer(8) :: nn
    real(dp) :: val, dA(3), dB(3), AA(3,3), AB(3,3), BB(3,3), m9(9,9), wd
    integer :: atomof(3)
    logical :: contract

    contract = present(hess)

    nca = NUM_CART_BF(la)
    ncb = NUM_CART_BF(lb)
    nn = int(nraw, 8)*int(nraw, 8)
    atomof = [iata, iatb, iatc]
    do cb_ = 1, ncb
      tb = ntup_upto(lb - 1) + cb_
      call tuple_of(tb, kb(1), kb(2), kb(3))
      do ca_ = 1, nca
        ta = ntup_upto(la - 1) + ca_
        call tuple_of(ta, ka(1), ka(2), ka(3))
        row = off(s1) + ca_
        col = off(s2) + cb_
        if (s1 == s2 .and. ca_ < cb_) cycle

        if (deriv == 0) then
          val = iw(W1, ta, tb)
          call put(0, row, col, val)
          cycle
        end if

        do k = 1, 3
          dA(k) = 2.0_dp*iw(WA, shifted(ka, k, 1), tb)
          if (ka(k) > 0) dA(k) = dA(k) - real(ka(k), dp)*iw(W1, shifted(ka, k, -1), tb)
          dB(k) = 2.0_dp*iw(WB, ta, shifted(kb, k, 1))
          if (kb(k) > 0) dB(k) = dB(k) - real(kb(k), dp)*iw(W1, ta, shifted(kb, k, -1))
        end do

        if (deriv == 1) then
          do k = 1, 3
            call put(3*(iata - 1) + k - 1, row, col, dA(k))
            call put(3*(iatb - 1) + k - 1, row, col, dB(k))
            call put(3*(iatc - 1) + k - 1, row, col, -(dA(k) + dB(k)))
          end do
          cycle
        end if

        do k = 1, 3
          do k2 = 1, 3
            AA(k, k2) = dd_same(ka, k, k2, .true.)
            BB(k, k2) = dd_same(kb, k, k2, .false.)
            AB(k, k2) = dd_mixed(k, k2)
          end do
        end do
        ! 9x9 over (centre A,B,C) x (x,y,z); C = -(A+B)
        m9(1:3, 1:3) = AA
        m9(1:3, 4:6) = AB
        m9(4:6, 1:3) = transpose(AB)
        m9(4:6, 4:6) = BB
        m9(7:9, 1:6) = -(m9(1:3, 1:6) + m9(4:6, 1:6))
        m9(1:9, 7:9) = -(m9(1:9, 1:3) + m9(1:9, 4:6))
        do x = 1, 3
          do y = 1, 3
            iat = atomof(x)
            jat = atomof(y)
            do ix = 1, 3
              do iy = 1, 3
                gi = 3*(iat - 1) + ix
                gj = 3*(jat - 1) + iy
                if (gi > gj) cycle
                if (contract) then
                  ! put() adds v at (row,col) and (col,row); dens is symmetric
                  if (row == col) then
                    wd = dens(row, col)
                  else
                    wd = 2.0_dp*dens(row, col)
                  end if
                  val = wd*m9(3*(x - 1) + ix, 3*(y - 1) + iy)
                  hess(gi, gj) = hess(gi, gj) + val
                  if (gi /= gj) hess(gj, gi) = hess(gj, gi) + val
                  cycle
                end if
                if (iat == jat) then
                  n0 = ecp_hess_start(iat - 1, iat - 1, natm) + 3
                  ! {xx,xy,xz,yy,yz,zz}
                  kk = (ix - 1)*3 - ((ix - 1)*(ix - 2))/2 + (iy - ix)
                  mat = n0 + kk
                else
                  n0 = ecp_hess_start(iat - 1, jat - 1, natm)
                  mat = n0 + (ix - 1)*3 + (iy - 1)
                end if
                call put(mat, row, col, m9(3*(x - 1) + ix, 3*(y - 1) + iy))
              end do
            end do
          end do
        end do
      end do
    end do

  contains

    subroutine put(mat_, r, cc, v)
      integer, intent(in) :: mat_, r, cc
      real(dp), intent(in) :: v
      integer(8) :: base
      base = int(mat_, 8)*nn
      res(base + int(r - 1, 8)*nraw + cc) = res(base + int(r - 1, 8)*nraw + cc) + v
      if (r /= cc) res(base + int(cc - 1, 8)*nraw + r) = res(base + int(cc - 1, 8)*nraw + r) + v
    end subroutine put

    integer function shifted(t, k_, d)
      integer, intent(in) :: t(3), k_, d
      integer :: u(3)
      u = t
      u(k_) = u(k_) + d
      shifted = tid(u(1), u(2), u(3))
    end function shifted

    !> d/dX_k2 d/dX_k of one side (A if isA, else B)
    real(dp) function dd_same(t, k_, k2_, isA)
      integer, intent(in) :: t(3), k_, k2_
      logical, intent(in) :: isA
      integer :: u(3), v(3), w1_, wx, wxx
      real(dp) :: f1
      integer :: s1_, s2_
      if (isA) then
        w1_ = W1; wx = WA; wxx = WAA
      else
        w1_ = W1; wx = WB; wxx = WBB
      end if
      dd_same = 0.0_dp
      ! D_k t = -t_k [t - e_k] + 2x [t + e_k]; apply D_k2 to each term
      do s1_ = -1, 1, 2
        u = t
        u(k_) = u(k_) + s1_
        if (s1_ == -1) then
          if (t(k_) == 0) cycle
          f1 = -real(t(k_), dp)
        else
          f1 = 2.0_dp
        end if
        do s2_ = -1, 1, 2
          v = u
          v(k2_) = v(k2_) + s2_
          if (s2_ == -1) then
            if (u(k2_) == 0) cycle
            if (s1_ == -1) then
              dd_same = dd_same + f1*(-real(u(k2_), dp))*side(v, w1_, isA)
            else
              dd_same = dd_same + f1*(-real(u(k2_), dp))*side(v, wx, isA)
            end if
          else
            if (s1_ == -1) then
              dd_same = dd_same + f1*2.0_dp*side(v, wx, isA)
            else
              dd_same = dd_same + f1*2.0_dp*side(v, wxx, isA)
            end if
          end if
        end do
      end do
    end function dd_same

    real(dp) function side(v, w_, isA)
      integer, intent(in) :: v(3), w_
      logical, intent(in) :: isA
      if (isA) then
        side = iw(w_, tid(v(1), v(2), v(3)), tb)
      else
        side = iw(w_, ta, tid(v(1), v(2), v(3)))
      end if
    end function side

    !> d/dA_k d/dB_k2
    real(dp) function dd_mixed(k_, k2_)
      integer, intent(in) :: k_, k2_
      integer :: u(3), v(3), sa, sb, wsel
      real(dp) :: fa, fb
      dd_mixed = 0.0_dp
      do sa = -1, 1, 2
        u = ka
        u(k_) = u(k_) + sa
        if (sa == -1) then
          if (ka(k_) == 0) cycle
          fa = -real(ka(k_), dp)
        else
          fa = 2.0_dp
        end if
        do sb = -1, 1, 2
          v = kb
          v(k2_) = v(k2_) + sb
          if (sb == -1) then
            if (kb(k2_) == 0) cycle
            fb = -real(kb(k2_), dp)
          else
            fb = 2.0_dp
          end if
          if (sa == -1 .and. sb == -1) then
            wsel = W1
          else if (sa == 1 .and. sb == -1) then
            wsel = WA
          else if (sa == -1 .and. sb == 1) then
            wsel = WB
          else
            wsel = WAB
          end if
          dd_mixed = dd_mixed + fa*fb*iw(wsel, tid(u(1), u(2), u(3)), tid(v(1), v(2), v(3)))
        end do
      end do
    end function dd_mixed

  end subroutine scatter

end module ecp_tool
