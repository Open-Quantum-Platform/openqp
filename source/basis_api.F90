
module basis_api
    use iso_c_binding, only: c_int, c_double, c_char, c_ptr, c_f_pointer, c_null_ptr, c_associated, c_loc
    implicit none

    type, bind(C) :: electron_shell
        integer(c_int) :: id
        integer(c_int) :: element_id
        character(c_char), dimension(20) :: function_type
        character(c_char), dimension(20) :: region
        integer(c_int) :: angular_momentum
        integer(c_int) :: n_exponents
        type(c_ptr) :: exponents
        type(c_ptr) :: coefficient
        type(c_ptr) :: next
    end type electron_shell

    type(c_ptr) :: head = c_null_ptr 

contains
    subroutine append_shell_C(c_handle) bind(C, name="append_shell")
        use c_interop, only: oqp_handle_t, oqp_handle_get_info
        use types, only: information
        type(oqp_handle_t) :: c_handle
        type(information), pointer :: inf
        inf => oqp_handle_get_info(c_handle)
        call oqp_append_shell(inf)
    end subroutine append_shell_C

    subroutine oqp_append_shell(info)
        use types, only: information
        type(information), intent(in) :: info
!        type(electron_shell), intent(in) :: new_shell
        type(electron_shell), pointer :: new_node
        type(c_ptr) :: current, new_node_ptr
        type(electron_shell), pointer :: temp

        allocate(new_node)
        new_node%id = info%elshell%id
        new_node%element_id = info%elshell%element_id
!        new_node%function_type = info%elshell%function_type
!        new_node%region = info%elshell%region
        new_node%angular_momentum = info%elshell%ang_mom
        new_node%n_exponents = info%elshell%num_expo
!        call c_f_pointer(infos%elshell%expo, , [infos%elshell%num_expo])
        new_node%exponents = info%elshell%expo
        new_node%coefficient = info%elshell%coef
        new_node%next = c_null_ptr

        new_node_ptr = c_loc(new_node)

        if (.not. c_associated(head)) then
            head = new_node_ptr
        else
            current = head
            do
                call c_f_pointer(current, temp)
                if (.not. c_associated(temp%next)) exit
                current = temp%next
            end do
            call c_f_pointer(current, temp)
            temp%next = new_node_ptr
        end if
        call print_all_shells()
    end subroutine oqp_append_shell

    subroutine print_all_shells() bind(C, name="print_all_shells")
        type(c_ptr) :: current
        type(electron_shell), pointer :: temp
        real(c_double), pointer :: exponents_ptr(:), coefficients_ptr(:)

        current = head
        print *, "Printing all shells:"
        do while (c_associated(current))
            call c_f_pointer(current, temp)
            print *, "Shell ID: ", temp%id
            print *, "element ID: ", temp%element_id
            print *, "Angular Momentum: ", temp%angular_momentum
            print *, "Number of Exponents: ", temp%n_exponents

            call c_f_pointer(temp%exponents, exponents_ptr, [temp%n_exponents])
            call c_f_pointer(temp%coefficient, coefficients_ptr, [temp%n_exponents])

            print *, "Exponents: "
            print *, exponents_ptr

            print *, "Coefficients: "
            print *, coefficients_ptr
            print *, "----------------------"
            current = temp%next
        end do
        print *, "----------------------"
        print *, "----------------------"
        print *, "----------------------"
    end subroutine print_all_shells

    subroutine delete_all_shells() bind(C, name="delete_all_shells")
        type(c_ptr) :: temp
        type(electron_shell), pointer :: to_delete

        do while (c_associated(head))
            temp = head
            call c_f_pointer(temp, to_delete)
            head = to_delete%next
            deallocate(to_delete)
        end do
        print *, "All shells deleted."
    end subroutine delete_all_shells

end module basis_api

