!> @brief ECP (effective core potential) interface built on libecpint.
!> @detail Provides ECP one-electron integrals and first derivatives, handling
!>         AO-label remapping and shell-origin geometry. Wraps libecpint’s C API
!>         and exposes simple Fortran-callable routines for OpenQP.
!> @author Mohsen Mazaherifar
!> @date January 2025
module ecp_tool
    use iso_c_binding, only: c_double, c_ptr, c_int,&
            c_f_pointer, C_LOC, c_null_ptr
    use, intrinsic :: iso_fortran_env, only: real64
    use libecpint_wrapper
    use libecp_result, only : ecp_result
    use basis_tools, only: basis_set
    use precision, only: dp

    implicit none

    private
    public add_ecpint
    public add_ecpder

contains
    !> @brief Add ECP one-electron contribution to the AO-core Hamiltonian (packed).
    !> @detail Computes scalar ECP integrals with libecpint (deriv order 0),
    !>         remaps them into OpenQP AO ordering via @ref transform_matrix,
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
        type(c_ptr) :: integrator
        type(ecp_result) :: result_ptr
        real(c_double), pointer :: libecp_res(:)
        integer :: i, j, c
        integer(c_int) :: driv_order

        if (.not.(basis%ecp_params%is_ecp)) then
            return
        end if
        driv_order = 0

        call set_integrator(integrator, basis, coord, driv_order)

        result_ptr = compute_integrals(integrator)
        call c_f_pointer(result_ptr%data, libecp_res, [result_ptr%size])


        call transform_matrix(basis, libecp_res)

        c = 0
        do i = 1, basis%nbf
            do j = 1, i
                c = c + 1
                hcore(c) = libecp_res((i - 1) * basis%nbf + j) + hcore(c)
            end do
        end do

        result_ptr%data = c_null_ptr
        result_ptr%size = 0
        nullify(libecp_res)

        call free_integrator(integrator)
        call free_result(result_ptr)

    end subroutine add_ecpint

    !> @brief Add ECP force contribution (first derivatives) to nuclear gradients.
    !> @detail Computes dV_ECP/dR_A in AO full-square form for each atom using
    !>         libecpint (deriv order 1), transforms to OpenQP AO ordering, and
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

        type(ecp_result) :: result_ptr
        type(c_ptr) :: integrator
        real(c_double), pointer :: libecp_res(:)
        real(real64), allocatable :: res_x(:), res_y(:), res_z(:)
        real(real64), allocatable :: deloc(:,:)
        integer :: i, j, c, n, natm, prim
        integer :: tri_size, full_size
        integer(c_int) :: driv_order

        if (.not.(basis%ecp_params%is_ecp)) then
            return
        end if

        driv_order = 1

        tri_size = basis%nbf * (basis%nbf + 1) / 2
        full_size = basis%nbf * basis%nbf

        allocate(res_x(full_size))
        allocate(res_y(full_size))
        allocate(res_z(full_size))

        natm = size(coord, dim=2)

        allocate(deloc(3, natm))
        deloc = 0

        call set_integrator(integrator, basis, coord, driv_order)

        result_ptr = compute_first_derivs(integrator)

        call c_f_pointer(result_ptr%data, libecp_res, [result_ptr%size])


        do n = 1, natm
            res_x = libecp_res(full_size * (3 * n - 3)+1:full_size * (3 * n - 2))
            res_y = libecp_res(full_size * (3 * n - 2)+1:full_size * (3 * n - 1))
            res_z = libecp_res(full_size * (3 * n - 1)+1:full_size * (3 * n))
            call transform_matrix(basis, res_x)
            call transform_matrix(basis, res_y)
            call transform_matrix(basis, res_z)

            do j = 1, basis%nbf
                do i = 1, j
                    c = j * (j - 1) / 2 + i

                    if (i == j) then
                        prim = 1
                    else
                        prim = 2
                    end if

                    deloc(1, n) = deloc(1, n) + prim * res_x((i - 1) * basis%nbf + j) * denab(c)
                    deloc(2, n) = deloc(2, n) + prim * res_y((i - 1) * basis%nbf + j) * denab(c)
                    deloc(3, n) = deloc(3, n) + prim * res_z((i - 1) * basis%nbf + j) * denab(c)
                end do
            end do

        end do

        de(:, 1:natm) = de(:, 1:natm) + deloc(:, 1:natm)

        result_ptr%data = c_null_ptr
        result_ptr%size = 0
        nullify(libecp_res)

        call free_integrator(integrator)
        call free_result(result_ptr)

    end subroutine add_ecpder

    !> @brief Construct and initialize a libecpint integrator instance.
    !> @detail Marshals Gaussian basis (centers, exponents, contractions, AMs) and
    !>         ECP basis (centers, exponents, coefficients, AMs, powers) from
    !>         OpenQP’s `basis_set` into libecpint arrays, assigns the ECP data,
    !>         and finalizes the integrator for the requested derivative order.
    !> @param[out] integrator    Opaque libecpint handle (C pointer).
    !> @param[in]  basis         Basis + ECP data.
    !> @param[in]  coord         Nuclear coordinates (3×natm).
    !> @param[in]  deriv_order   0 = value, 1 = first derivatives.
    !> @pre `basis%ecp_params` fields are allocated when is_ecp is true.
    !> @author Mohsen Mazaherifar
    !> @date January 2025
    subroutine set_integrator(integrator, basis, coord, deriv_order)

        real(c_double), intent(in), contiguous :: coord(:,:)
        type(basis_set), intent(in) :: basis
        integer(c_int), intent(in) :: deriv_order

        type(c_ptr) :: integrator
        real(c_double), allocatable :: g_coords(:), g_exps(:), g_coefs(:)
        integer(c_int), allocatable :: g_ams(:), g_lengths(:)
        real(c_double), allocatable :: u_coords(:), u_exps(:), u_coefs(:)
        integer(c_int), allocatable :: u_ams(:), u_ns(:), u_lengths(:)
        integer(c_int) :: num_ecps, num_gaussians, n_coord, f_expo_len
        integer :: tri_size, full_size, natm


        tri_size = basis%nbf * (basis%nbf + 1) / 2
        full_size = basis%nbf * basis%nbf

        f_expo_len = sum(basis%ecp_params%n_expo)
        natm = size(coord, dim=2)

        num_gaussians = basis%nshell
        n_coord = num_gaussians * 3

        allocate(g_coords(n_coord), g_exps(basis%nprim), g_coefs(basis%nprim))
        allocate(g_ams(basis%nshell), g_lengths(basis%nshell))

        allocate(u_coords(size(basis%ecp_params%ecp_coord)), u_exps(f_expo_len))
        allocate(u_coefs(f_expo_len), u_ams(f_expo_len))
        allocate(u_ns(f_expo_len), u_lengths(size(basis%ecp_params%n_expo)))

        call libecp_g_coords(basis, coord, g_coords)
        g_exps = real(basis%ex, kind=c_double)
        g_coefs = real(basis%cc, kind=c_double)
        g_ams = int(basis%am, kind=c_int)
        g_lengths = int(basis%ncontr, kind=c_int)

        num_ecps = int(size(basis%ecp_params%n_expo), kind=c_int)
        u_coords = real(basis%ecp_params%ecp_coord, kind=c_double)
        u_exps = real(basis%ecp_params%ecp_ex, kind=c_double)
        u_coefs = real(basis%ecp_params%ecp_cc, kind=c_double)
        u_ams = int(basis%ecp_params%ecp_am, kind=c_int)
        u_ns = int(basis%ecp_params%ecp_r_ex, kind=c_int)
        u_lengths = int(basis%ecp_params%n_expo, kind=c_int)


        integrator = init_integrator(num_gaussians, g_coords, g_exps, g_coefs, &
                                     g_ams, g_lengths)

        call set_ecp_basis(integrator, num_ecps, u_coords, u_exps, u_coefs, &
                           u_ams, u_ns, u_lengths)

        call init_integrator_instance(integrator, deriv_order)

    end subroutine set_integrator
  !> @brief Build AO index remapping from libecpint canonical order to OpenQP AO order.
  !> @detail Fills `label_map(i_old)=i_new` using shell origins and angular-momentum
  !>         layout so that full-square AO matrices can be permuted consistently.
  !> @param[in]    basis     Basis set (AO layout and shell metadata).
  !> @param[inout] label_map Integer array of length nbf receiving the permutation.
  !> @see transform_matrix
  !> @author Mohsen Mazaherifar
  !> @date January 2025
  subroutine libecpint_map(basis, label_map)

    use basis_tools, only: basis_set
    use constants, only: map_canonical

    type(basis_set), intent(in) :: basis
    integer, dimension(:), intent(inout):: label_map
    integer :: i,j

    do i = 1, basis%nbf
      j = basis%bf_to_shell(i)
      label_map(i + map_canonical(i-basis%ao_offset(j)+1, basis%am(j))) = i
    end do

  end  subroutine libecpint_map
  !> @brief Pack Gaussian-center coordinates per shell for libecpint.
  !> @detail Writes (x,y,z) per shell index using `basis%origin(shell)` to select
  !>         the parent atom for the shell center as expected by libecpint.
  !> @param[in]  basis    Basis set.
  !> @param[in]  coord    Nuclear coordinates (3×natm).
  !> @param[out] g_coords Flat array of size 3*nshell: [x1,y1,z1, x2,y2,z2, ...].
  !> @note Coordinates are cast to C double precision for the C API.
  !> @author Mohsen Mazaherifar
  !> @date January 2025
  subroutine libecp_g_coords(basis, coord, g_coords)

      type(basis_set), intent(in) :: basis
      real(real64), intent(in) :: coord(:,:)
      real(c_double), intent(out) :: g_coords(:)

      integer :: shell

      do shell = 1, basis%nshell
          g_coords(3*shell-2:3*shell) = real(coord(1:3, basis%origin(shell)),c_double)
      end do
  end subroutine libecp_g_coords

  !> @brief Permute a full AO square matrix into OpenQP AO ordering.
  !> @detail Applies the mapping from @ref libecpint_map to reorder rows/cols
  !>         of `matrix` in-place (via a temporary copy). Expects size nbf×nbf.
  !> @param[in]    basis   Basis set (provides AO label map).
  !> @param[inout] matrix  Full AO square matrix flattened (size nbf*nbf).
  !> @throws Stops if `size(matrix) != nbf*nbf`.
  !> @see libecpint_map
  !> @author Mohsen Mazaherifar
  !> @date January 2025
  subroutine transform_matrix(basis, matrix)

    use basis_tools, only: basis_set
    type(basis_set), intent(in) :: basis
    real(c_double), dimension(:), intent(inout) :: matrix
    real(c_double), dimension(:), allocatable :: tmp_matrix
    integer, dimension(:), allocatable :: label_map
    integer :: i, j, row, col, n

    n = basis%nbf
    allocate(label_map(n))

    if (size(matrix) /= n * n) then
      print *, "Error: original_matrix size does not match labels."
      stop
    end if

    call libecpint_map(basis, label_map)
    allocate(tmp_matrix(n * n))

    tmp_matrix = matrix

    do i = 1, n
      do j = 1, n
        row = label_map(i)
        col = label_map(j)
        matrix((row - 1) * n + col) = tmp_matrix((i - 1) * n + j)
      end do
    end do
  end subroutine transform_matrix


end module ecp_tool
