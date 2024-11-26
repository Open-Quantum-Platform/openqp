module ecp_tool
    use iso_c_binding, only: c_double, c_ptr, c_int, c_f_pointer, C_LOC
    use, intrinsic :: iso_fortran_env, only: real64
    use libecpint_wrapper
    use libecp_result, only : ecp_result
    use basis_api, only : ecp_head
    use basis_tools, only: basis_set

    implicit none

contains
    subroutine add_ecpint(basis,coord ,hcore)
!        use types, only: information
!        type(information), intent(in) :: info
        real(real64), contiguous,  intent(in)  :: coord(:,:)
        type(basis_set), intent(in) :: basis
        type(c_ptr) :: integrator
        real(real64), contiguous, intent(inout) :: hcore(:)
        type(ecp_result) :: result_ptr
        real(c_double), pointer :: result1(:)
        real(c_double), dimension(:), allocatable :: result2
        integer :: i, j, c
        real(real64) :: tmp

        real(c_double), allocatable :: g_coords(:)
        real(c_double), allocatable :: g_exps(:)
        real(c_double), allocatable :: g_coefs(:)
        integer(c_int), allocatable :: g_ams(:), g_lengths(:)
        real(c_double), allocatable :: u_coords(:)
        real(c_double), allocatable :: u_exps(:), u_coefs(:)
        integer(c_int), allocatable  :: u_ams(:), u_ns(:)
        integer(c_int),  allocatable :: u_lengths(:)
        integer(c_int) :: num_ecps, num_gaussians, n_coord

        if (ecp_head%element_id .EQ. 0) then
            return
        end if

        num_gaussians = basis%nshell
        n_coord = num_gaussians*3
        allocate(g_coords(n_coord))
        allocate(g_exps(basis%nprim))
        allocate(g_coefs(basis%nprim))
        allocate(g_ams(basis%nshell))
        allocate(g_lengths(basis%nshell))

        allocate(u_coords(size(ecp_head%ecp_coord)))
        allocate(u_exps(ecp_head%n_exponents))
        allocate(u_coefs(ecp_head%n_exponents))
        allocate(u_ams(ecp_head%n_exponents))
        allocate(u_ns(ecp_head%n_exponents))
        allocate(u_lengths(ecp_head%element_id))

        call libecp_g_coords(basis, coord ,g_coords)

        g_exps = real(basis%ex, kind=c_double)
        g_coefs = real(basis%cc, kind=c_double)
        g_ams = int(basis%am, kind=c_int)
        g_lengths = int(basis%ncontr, kind=c_int)


        num_ecps = int(size(ecp_head%ecp_coord)/3, kind =c_int)
        u_coords =  real(ecp_head%ecp_coord, kind=c_double)
        u_exps =  real(ecp_head%exponents, kind=c_double)
        u_coefs =  real(ecp_head%coefficient, kind=c_double)
        u_ams =  int(ecp_head%ecp_am, kind=c_int)
        u_ns =  int(ecp_head%ecp_r_expo, kind=c_int)
        u_lengths(1) = ecp_head%n_exponents


        integrator = init_integrator(num_gaussians, g_coords, g_exps, g_coefs,  &
                                g_ams, g_lengths)

        call set_ecp_basis(integrator, num_ecps,u_coords, u_exps, u_coefs,    &
                 u_ams, u_ns, u_lengths)
        call init_integrator_instance(integrator)

        result_ptr = compute_integrals(integrator)

        call c_f_pointer(result_ptr%data, result1, [result_ptr%size])
        do i = 1, basis%nbf
            print *, basis%bf_label(i)
        end do


        call transform_matrix(result1, basis, result2)
        c=0
        do i = 1 ,basis%nbf!size(g_coefs)
            do j = 1, i!size(g_coefs)
                c =c+1
                tmp = hcore(c)
                hcore(c) = result2((i-1)*basis%nbf +j) + hcore(c)
                print *, c, result2((i-1)*basis%nbf +j),tmp, hcore(c)
            end do
        end do

        call free_integrator(integrator)
        call free_result(result_ptr)
    end subroutine add_ecpint

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
      
      use physical_constants, only: BOHR_TO_ANGSTROM
      type(basis_set), intent(in) :: basis
      real(real64), intent(in) :: coord(:,:)
      real(c_double), intent(out) :: g_coords(:)

      integer :: shell 

      do shell = 1, basis%nshell
          g_coords(3*shell-2:3*shell) = real(coord(1:3, basis%origin(shell)),c_double)!*BOHR_TO_ANGSTROM
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
