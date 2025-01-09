module ecp_tool
    use iso_c_binding, only: c_double, c_ptr, c_int,&
            c_f_pointer, C_LOC, c_null_ptr
    use, intrinsic :: iso_fortran_env, only: real64
    use libecpint_wrapper
    use libecp_result, only : ecp_result
    use basis_api, only : ecp_head
    use basis_tools, only: basis_set
    use precision, only: dp

    implicit none

contains

subroutine add_ecpint(basis, coord, hcore)
    real(real64), contiguous, intent(in) :: coord(:,:)
    type(basis_set), intent(in) :: basis
    real(real64), contiguous, intent(inout) :: hcore(:)
    type(c_ptr) :: integrator
    type(ecp_result) :: result_ptr
    real(c_double), pointer :: result1(:)
    real(c_double), dimension(:), allocatable :: result2
    integer :: i, j, c
    integer(c_int) :: driv_order

    ! Early exit if no ECP data is present
    if (ecp_head%element_id == 0) then
        return
    end if
    driv_order = 0

    ! Use the set_integrator function to initialize the integrator
    call set_integrator(integrator, basis, coord, driv_order)

    result_ptr = compute_integrals(integrator)
    call c_f_pointer(result_ptr%data, result1, [result_ptr%size])

    call transform_matrix(result1, basis, result2)

    c = 0
    do i = 1, basis%nbf
        do j = 1, i
            c = c + 1
            hcore(c) = result2((i - 1) * basis%nbf + j) + hcore(c)
        end do
    end do

    result_ptr%data = c_null_ptr
    result_ptr%size = 0
    nullify(result1)

    call free_integrator(integrator)
    call free_result(result_ptr)

end subroutine add_ecpint

subroutine add_ecpder(basis, coord, denab, de)

    real(real64), contiguous, intent(in) :: coord(:,:)
    type(basis_set), intent(in) :: basis
    REAL(kind=dp), INTENT(INOUT) :: denab(:)
    REAL(kind=dp), intent(INOUT) :: de(:,:)

    type(ecp_result) :: result_ptr
    type(c_ptr) :: integrator
    real(c_double), pointer :: result1(:)
    real(real64), allocatable :: deloc(:,:)
    integer :: i, j, c, n, natm
    integer :: tri_size, full_size
    integer(c_int) :: driv_order

    ! Early exit if no ECP data is present
    if (ecp_head%element_id == 0) then
        return
    end if
    driv_order = 1

    tri_size = basis%nbf * (basis%nbf + 1) / 2
    full_size = basis%nbf * basis%nbf
    natm = size(coord, dim=2)

    allocate(deloc(3, natm))
    deloc = 0

    ! Use the set_integrator function to initialize the integrator
    call set_integrator(integrator, basis, coord, driv_order)

    result_ptr = compute_first_derivs(integrator)

    call c_f_pointer(result_ptr%data, result1, [result_ptr%size])

    ! Compute derivatives
    do i = 1, basis%nbf
        do j = 1, basis%nbf
            do n = 1, natm
                if (i <= j) then
                    c = j * (j - 1) / 2 + i
                else
                    c = i * (i - 1) / 2 + j
                end if

                deloc(1, n) = deloc(1, n) + &
                        result1(full_size * (3 * n - 3) + (i - 1) * basis%nbf + j) * denab(c)
                deloc(2, n) = deloc(2, n) + &
                        result1(full_size * (3 * n - 2) + (i - 1) * basis%nbf + j) * denab(c)
                deloc(3, n) = deloc(3, n) + &
                        result1(full_size * (3 * n - 1) + (i - 1) * basis%nbf + j) * denab(c)

            end do
        end do
    end do

    do n = 1, natm
        de(1, n) = de(1, n) + deloc(1, n)
        de(2, n) = de(2, n) + deloc(2, n)
        de(3, n) = de(3, n) + deloc(3, n)
    end do

    result_ptr%data = c_null_ptr
    result_ptr%size = 0
    nullify(result1)

    call free_integrator(integrator)
    call free_result(result_ptr)

end subroutine add_ecpder

    subroutine set_integrator(integrator, basis, coord, deriv_order)
    
        ! Input arguments
        real(c_double), intent(in), contiguous :: coord(:,:)
        type(basis_set), intent(in) :: basis
        integer(c_int), intent(in) :: deriv_order
    
        ! Local variables
        type(c_ptr) :: integrator
        real(c_double), allocatable :: g_coords(:), g_exps(:), g_coefs(:)
        integer(c_int), allocatable :: g_ams(:), g_lengths(:)
        real(c_double), allocatable :: u_coords(:), u_exps(:), u_coefs(:)
        integer(c_int), allocatable :: u_ams(:), u_ns(:), u_lengths(:)
        integer(c_int) :: num_ecps, num_gaussians, n_coord, f_expo_len
        integer :: tri_size, full_size, natm
    
    
        ! Compute sizes and allocate memory
        tri_size = basis%nbf * (basis%nbf + 1) / 2
        full_size = basis%nbf * basis%nbf
    
        f_expo_len = sum(ecp_head%n_exponents)
        natm = size(coord, dim=2)
    
        num_gaussians = basis%nshell
        n_coord = num_gaussians * 3
    
        allocate(g_coords(n_coord), g_exps(basis%nprim), g_coefs(basis%nprim))
        allocate(g_ams(basis%nshell), g_lengths(basis%nshell))
    
        allocate(u_coords(size(ecp_head%ecp_coord)), u_exps(f_expo_len))
        allocate(u_coefs(f_expo_len), u_ams(f_expo_len))
        allocate(u_ns(f_expo_len), u_lengths(ecp_head%element_id))
    
        ! Initialize Gaussian basis data
        call libecp_g_coords(basis, coord, g_coords)
        g_exps = real(basis%ex, kind=c_double)
        g_coefs = real(basis%cc, kind=c_double)
        g_ams = int(basis%am, kind=c_int)
        g_lengths = int(basis%ncontr, kind=c_int)
    
        ! Initialize ECP data
        num_ecps = int(size(ecp_head%ecp_coord) / 3, kind=c_int)
        u_coords = real(ecp_head%ecp_coord, kind=c_double)
        u_exps = real(ecp_head%exponents, kind=c_double)
        u_coefs = real(ecp_head%coefficient, kind=c_double)
        u_ams = int(ecp_head%ecp_am, kind=c_int)
        u_ns = int(ecp_head%ecp_r_expo, kind=c_int)
        u_lengths = ecp_head%n_exponents
    
        ! Initialize integrator
        integrator = init_integrator(num_gaussians, g_coords, g_exps, g_coefs, &
                                     g_ams, g_lengths)
    
        call set_ecp_basis(integrator, num_ecps, u_coords, u_exps, u_coefs, &
                           u_ams, u_ns, u_lengths)
    
        call init_integrator_instance(integrator, deriv_order)
    
    end subroutine set_integrator

  subroutine libecpint_map(basis, label_map)

    use basis_tools, only: basis_set
    type(basis_set), intent(in) :: basis

    integer, dimension(:), intent(inout):: label_map
    character(len=8) :: extact_label

    integer :: i
    do i = 1, basis%nbf
      extact_label = basis%bf_label(i)
      select case (trim(extact_label(6:8)))
        case ('YY')
          label_map(i+2) = i
        case ('ZZ')
          label_map(i+3) = i
        case ('XY')
          label_map(i-2) = i
        case ('XZ')
          label_map(i-2) = i
        case ('YZ')
          label_map(i-1) = i
        case default
          label_map(i) = i
      end select
    end do

  end  subroutine libecpint_map

  subroutine libecp_g_coords(basis, coord, g_coords)

      type(basis_set), intent(in) :: basis
      real(real64), intent(in) :: coord(:,:)
      real(c_double), intent(out) :: g_coords(:)

      integer :: shell

      do shell = 1, basis%nshell
          g_coords(3*shell-2:3*shell) = real(coord(1:3, basis%origin(shell)),c_double)
      end do
  end subroutine libecp_g_coords

  subroutine transform_matrix(original_matrix, basis, transformed_matrix)

    use basis_tools, only: basis_set
    type(basis_set), intent(in) :: basis
    real(c_double), dimension(:), intent(in) :: original_matrix
    real(c_double), dimension(:), allocatable, intent(out) :: transformed_matrix
    integer, dimension(:), allocatable :: label_map
    integer :: i, j, row, col, n

    n = basis%nbf
    allocate(label_map(n))

    if (size(original_matrix) /= n * n) then
      print *, "Error: original_matrix size does not match labels."
      stop
    end if

    call libecpint_map(basis, label_map)
    allocate(transformed_matrix(n * n))

    transformed_matrix = 0.0_c_double

    do i = 1, n
      do j = 1, n
        row = label_map(i)
        col = label_map(j)
        transformed_matrix((row - 1) * n + col) = original_matrix((i - 1) * n + j)
      end do
    end do
  end subroutine transform_matrix


end module ecp_tool
