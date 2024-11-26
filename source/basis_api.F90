module basis_api
    use iso_c_binding, only: c_f_pointer, c_ptr, c_double
    use iso_fortran_env, only: real64
    use libecpint_wrapper
    implicit none

    type, abstract :: base_shell
        integer :: id
        integer :: element_id
        integer :: n_exponents
        real, pointer :: exponents(:)
        real, pointer :: coefficient(:)
    end type base_shell

    type, extends(base_shell) :: electron_shell
        integer :: angular_momentum
        type(electron_shell), pointer :: next => null()
    end type electron_shell

    type, extends(base_shell) :: ecpdata
        integer :: n_angular_m
        integer, pointer :: ecp_zn(:)
        integer, pointer :: ecp_r_expo(:)
        integer, pointer :: ecp_am(:)
        real, pointer :: ecp_coord(:)
!        type(ecpdata), pointer :: next => null()
    end type ecpdata

    type(electron_shell), pointer :: head => null()
    type(ecpdata) :: ecp_head
!    type(ecpdata), pointer :: ecp_head => null()

contains

    subroutine append_shell_C(c_handle) bind(C, name="append_shell")
        use c_interop, only: oqp_handle_t, oqp_handle_get_info
        use types, only: information
        type(oqp_handle_t) :: c_handle
        type(information), pointer :: inf
        inf => oqp_handle_get_info(c_handle)
        call oqp_append_shell(inf)
    end subroutine append_shell_C

    subroutine append_ecp_C(c_handle) bind(C, name="append_ecp")
        use c_interop, only: oqp_handle_t, oqp_handle_get_info
        use types, only: information
        type(oqp_handle_t) :: c_handle
        type(information), pointer :: inf
        inf => oqp_handle_get_info(c_handle)
        call oqp_append_ecp(inf)
    end subroutine append_ecp_C

    ! Append a new ECP  to the linked list
    subroutine oqp_append_ecp(info)
        use types, only: information
        type(information), intent(in) :: info
 !       type(ecpdata), pointer :: ecp_node, ecp_temp
        real(c_double), pointer :: expo_ptr(:), coef_ptr(:), rexpo_ptr(:),&
                am_ptr(:), coord_ptr(:)
        integer(c_int) , pointer :: ecp_zn_ptr(:)
        integer :: natm
!        call c_f_pointer(info%elshell%ecp_zn, ecp_zn_ptr, [natm])
        print *, "info%mol_prop%natom", info%mol_prop%natom
        natm = info%mol_prop%natom

!        natm = get_number_atoms()
        call c_f_pointer(info%elshell%ecp_zn, ecp_zn_ptr, [natm])
        print *, "ecp_zn_ptr", ecp_zn_ptr        
        allocate(ecp_head%ecp_zn(natm))
        ecp_head%ecp_zn = ecp_zn_ptr
        print *,"segmentation", ecp_head%ecp_zn

        if (info%elshell%element_id .EQ. 0) then
            ecp_head%element_id = info%elshell%element_id
            print *,"HIHIHIHI"
            return
        end if
        call c_f_pointer(info%elshell%expo, expo_ptr, [info%elshell%num_expo])
        call c_f_pointer(info%elshell%coef, coef_ptr, [info%elshell%num_expo])
        call c_f_pointer(info%elshell%ecp_rex, rexpo_ptr, [info%elshell%num_expo])
        call c_f_pointer(info%elshell%ecp_am, am_ptr, [info%elshell%ecp_nam])
        call c_f_pointer(info%elshell%ecp_coord, coord_ptr, [3*info%elshell%element_id])
!        call c_f_pointer(info%elshell%ecp_zn, ecp_zn_ptr, [natm])

        ecp_head%element_id = info%elshell%element_id
        ecp_head%n_exponents = info%elshell%num_expo
        ecp_head%n_angular_m = info%elshell%ecp_nam
        allocate(ecp_head%exponents(info%elshell%num_expo))
        allocate(ecp_head%coefficient(info%elshell%num_expo))
        allocate(ecp_head%ecp_r_expo(info%elshell%num_expo))
        allocate(ecp_head%ecp_am(info%elshell%ecp_nam))
        allocate(ecp_head%ecp_coord(3 * info%elshell%element_id))
!        allocate(ecp_head%ecp_zn(natm))
        print *, "head%element_id", ecp_zn_ptr
        ecp_head%exponents = expo_ptr
        ecp_head%coefficient = coef_ptr
        ecp_head%ecp_r_expo = rexpo_ptr
        ecp_head%ecp_am = am_ptr
        ecp_head%ecp_coord = coord_ptr
        !ecp_head%ecp_zn = ecp_zn_ptr
        print *,"aaaa", ecp_head%ecp_zn
!        print *, "ecp_node%n_exponents",ecp_head%n_exponents
!        print *, "ecp_node%ecp_am", ecp_head%ecp_am
!        print *, "ecp_node%ecp_r_expo", ecp_head%ecp_r_expo
!        allocate(ecp_node)
!        ecp_node%id = info%elshell%id
!        ecp_node%element_id = info%elshell%element_id
!        ecp_node%n_exponents = info%elshell%num_expo
!        ecp_node%n_angular_m = info%elshell%ecp_nam
!        allocate(ecp_node%exponents(info%elshell%num_expo))
!        allocate(ecp_node%coefficient(info%elshell%num_expo))
!        allocate(ecp_node%ecp_r_expo(info%elshell%num_expo))
!        allocate(ecp_node%ecp_am(info%elshell%ecp_nam))
!        allocate(ecp_node%ecp_coord(3 * info%elshell%ecp_nam))
!        ecp_node%exponents = expo_ptr
!        ecp_node%coefficient = coef_ptr
!        ecp_node%ecp_r_expo = rexpo_ptr
!        ecp_node%ecp_am = am_ptr
!        ecp_node%ecp_coord = coord_ptr
!        print *, "ecp_node%n_exponents",ecp_node%n_exponents
!        print *, "ecp_node%ecp_am", ecp_node%ecp_am
!        print *, "ecp_node%ecp_r_expo", ecp_node%ecp_r_expo

!        ecp_node%next => null()
!        if (.not. associated(ecp_head)) then
!            ecp_head => ecp_node
!        else
!            ecp_temp => ecp_head
!            do while (associated(ecp_temp%next))
!                ecp_temp => ecp_temp%next
!            end do
!            ecp_temp%next => ecp_node
!        end if
    end subroutine oqp_append_ecp

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

    function get_number_atoms() result(natm)
        type(electron_shell), pointer :: temp
        integer:: natm
        temp => head
        natm = 0
        do while (associated(temp%next))
            temp => temp%next
        end do
        if (associated(temp)) then
            natm = temp%element_id
        end if

        nullify(temp)
     end function get_number_atoms


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

        print *, "map_shell"

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
        if (.not. allocated(basis%ecp_zn_num)) allocate(basis%ecp_zn_num(maxval(basis%origin)))

        basis%ecp_zn_num = ecp_head%ecp_zn 
        print *, "ecp_head", basis%ecp_zn_num
        print *, "ecp_head", maxval(basis%origin)

        call basis%set_bfnorms()
        call basis%normalize_primitives()
!        call basis%init_shell_centers()


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

