module basis_api
    use iso_c_binding, only: c_f_pointer, c_ptr, c_double
    use iso_fortran_env, only: real64
    implicit none

    type :: electron_shell
        integer :: id
        integer :: element_id
        character(len=20) :: function_type
        character(len=20) :: region
        integer :: angular_momentum
        integer :: n_exponents
        real, pointer :: exponents(:)   ! Pointer to an array of exponents
        real, pointer :: coefficient(:) ! Pointer to an array of coefficients
        type(electron_shell), pointer :: next => null()  ! Fortran pointer for linked list
    end type electron_shell

    type(electron_shell), pointer :: head => null()  ! Head of the linked list

contains

    subroutine append_shell_C(c_handle) bind(C, name="append_shell")
        use c_interop, only: oqp_handle_t, oqp_handle_get_info
        use types, only: information
        type(oqp_handle_t) :: c_handle
        type(information), pointer :: inf
        inf => oqp_handle_get_info(c_handle)
        call oqp_append_shell(inf)
    end subroutine append_shell_C

    ! Append a new electron shell to the linked list
    subroutine oqp_append_shell(info)
        use types, only: information
        type(information), intent(in) :: info
        type(electron_shell), pointer :: new_node, temp
        real(c_double), pointer :: expo_ptr(:), coef_ptr(:)
        ! Map c_ptr fields from info to Fortran pointers
        call c_f_pointer(info%elshell%expo, expo_ptr, [info%elshell%num_expo])
        call c_f_pointer(info%elshell%coef, coef_ptr, [info%elshell%num_expo])

        allocate(new_node)
        new_node%id = info%elshell%id
        new_node%element_id = info%elshell%element_id
        new_node%angular_momentum = info%elshell%ang_mom
        allocate(new_node%exponents(info%elshell%num_expo))
        allocate(new_node%coefficient(info%elshell%num_expo))
        new_node%n_exponents = info%elshell%num_expo
        new_node%exponents = expo_ptr
        new_node%coefficient = coef_ptr
        new_node%next => null()
        if (.not. associated(head)) then
            head => new_node
        else
            temp => head
            do while (associated(temp%next))
                temp => temp%next
            end do
            temp%next => new_node
        end if
    end subroutine oqp_append_shell

    ! Print all shells in the linked list
    subroutine print_all_shells() bind(C, name="print_all_shells")
        type(electron_shell), pointer :: temp

        temp => head
        print *, "Printing all shells:"
        do while (associated(temp))

            print *, "Shell ID: ", temp%id
            print *, "Element ID: ", temp%element_id
            print *, "Angular Momentum: ", temp%angular_momentum
            print *, "Number of Exponents: ", temp%n_exponents
            print *, "Exponents: ", temp%exponents
            print *, "Coefficients: ", temp%coefficient
            print *, "----------------------"
            temp => temp%next

        end do
        print *, "----------------------"
    end subroutine print_all_shells

    subroutine map_shell2basis_set(basis)
        use basis_tools, only: basis_set
        class(basis_set) ,intent(inout):: basis
        type(electron_shell), pointer :: temp
        type(electron_shell), pointer :: temp1
        integer :: nbf, nshell, nprim, mxcontr, mxam, ii
        integer :: n1,n2
        real, dimension(:), allocatable :: ex

        temp => head
        mxam = 0
        mxcontr = 0
        nbf = 0
        nshell = 0
        nprim = 0  ! Initialize nprim
        ii = 0

        do while (associated(temp))

            mxcontr = max(mxcontr, temp%n_exponents)
            mxam = max(mxam, temp%angular_momentum)

            nshell = temp%id
            nprim = nprim + temp%n_exponents

            select case (temp%angular_momentum)
                case (0)
                    nbf = nbf + 1
                case (1)
                    nbf = nbf + 3
                case (2)
                    nbf = nbf + 6
                case (3)
                    nbf = nbf + 9
            end select
            temp => temp%next  ! Move to the next shell
        end do
        temp1 => head
        basis%mxam = mxam
        basis%mxcontr = mxcontr
        basis%nbf = nbf
        basis%nshell = nshell
        basis%nprim = nprim
        ! Allocate arrays based on the received sizes (on all processes)
        if (.not. allocated(basis%ex)) allocate(basis%ex(nprim))
        if (.not. allocated(basis%cc)) allocate(basis%cc(nprim))
        if (.not. allocated(basis%bfnrm)) allocate(basis%bfnrm(nbf))
        if (.not. allocated(basis%g_offset)) allocate(basis%g_offset(nshell))
        if (.not. allocated(basis%origin)) allocate(basis%origin(nshell))
        if (.not. allocated(basis%am)) allocate(basis%am(nshell))
        if (.not. allocated(basis%ncontr)) allocate(basis%ncontr(nshell))
        if (.not. allocated(basis%ao_offset)) allocate(basis%ao_offset(nshell))
        if (.not. allocated(basis%naos)) allocate(basis%naos(nshell))
        if (.not. allocated(basis%at_mx_dist2)) allocate(basis%at_mx_dist2(nbf))
        if (.not. allocated(basis%prim_mx_dist2)) allocate(basis%prim_mx_dist2(nprim))
        if (.not. allocated(basis%shell_mx_dist2)) allocate(basis%shell_mx_dist2(nshell))
        if (.not. allocated(basis%shell_centers)) allocate(basis%shell_centers(nshell, 3))

        if (.not. allocated(ex)) allocate(ex(nprim))



        do while (associated(temp1))
            ii = ii + 1
            n2 = temp1%n_exponents
            basis%ncontr(ii) = n2

            if (ii == 1) then
               basis%g_offset(1) = 1
               basis%ao_offset(ii) = 1
            else
               basis%g_offset(ii) = basis%g_offset(ii-1) + n1
               basis%ao_offset(ii) = basis%naos(ii-1) + basis%ao_offset(ii-1)
            end if

            basis%ex(basis%g_offset(ii):(basis%g_offset(ii) + n2 - 1)) = real(temp1%exponents, kind=real64)
            basis%cc(basis%g_offset(ii):(basis%g_offset(ii) + n2 - 1)) = real(temp1%coefficient,  kind=real64)



            basis%origin(ii) = temp1%element_id
            basis%am(ii) = temp1%angular_momentum
            select case (temp1%angular_momentum)
                case (0)
                    basis%naos(ii) = 1
                case (1)
                    basis%naos(ii) = 3
                case (2)
                    basis%naos(ii) = 6
                case (3)
                    basis%naos(ii) = 9
            end select

            n1 = temp1%n_exponents

            temp1 => temp1%next

        end do
        call basis%set_bfnorms()
        call basis%normalize_primitives()


    end subroutine map_shell2basis_set

    ! Delete all shells in the linked list
    subroutine delete_all_shells()
        type(electron_shell), pointer :: to_delete

        do while (associated(head))
            to_delete => head
            head => head%next
            deallocate(to_delete%exponents)
            deallocate(to_delete%coefficient)
            deallocate(to_delete)
        end do
        print *, "All shells deleted."
    end subroutine delete_all_shells

end module basis_api

