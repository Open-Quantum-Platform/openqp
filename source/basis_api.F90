module basis_api
    use iso_c_binding, only: c_f_pointer, c_ptr, c_double
    use iso_fortran_env, only: real64
    use physical_constants, only: UNITS_ANGSTROM
    use libecpint_wrapper
    implicit none

!###############################################################################

    type, abstract :: base_shell
        integer :: id
        integer :: element_id
        integer, pointer :: n_exponents(:)
        real, pointer :: exponents(:)
        real, pointer :: coefficient(:)
    contains
        procedure(base_shell_clear), deferred, pass :: clear
    end type base_shell

    abstract interface
        subroutine base_shell_clear(this)
            import base_shell
            class(base_shell), intent(inout) :: this
        end subroutine base_shell_clear
    end interface

!###############################################################################

    type, extends(base_shell) :: electron_shell
        integer :: angular_momentum
        type(electron_shell), pointer :: next => null()
   contains
      procedure :: clear => electron_shell_clear
    end type electron_shell

!###############################################################################

    type, extends(base_shell) :: ecpdata
        integer :: n_angular_m
        integer, pointer :: ecp_zn(:)
        integer, pointer :: ecp_r_expo(:)
        integer, pointer :: ecp_am(:)
        real, pointer :: ecp_coord(:)
   contains
      procedure :: clear => ecpdata_clear
    end type ecpdata

    type(electron_shell), pointer :: head => null()
    type(ecpdata) :: ecp_head

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

    subroutine oqp_append_ecp(info)
        use types, only: information
        type(information), intent(in) :: info
        real(c_double), pointer :: expo_ptr(:), coef_ptr(:), rexpo_ptr(:),&
                am_ptr(:), coord_ptr(:)
        integer(c_int) , pointer :: n_expo_ptr(:), ecp_zn_ptr(:)
        integer :: natm, f_expo_len
        natm = info%mol_prop%natom

        call c_f_pointer(info%elshell%ecp_zn, ecp_zn_ptr, [natm])
        allocate(ecp_head%ecp_zn(natm))
        ecp_head%ecp_zn = ecp_zn_ptr

        if (info%elshell%element_id .EQ. 0) then
            ecp_head%element_id = info%elshell%element_id
            return
        end if
        call c_f_pointer(info%elshell%num_expo, n_expo_ptr, [info%elshell%element_id])

        allocate(ecp_head%n_exponents(info%elshell%element_id))

        ecp_head%n_exponents = n_expo_ptr

        f_expo_len = sum(ecp_head%n_exponents)

        call c_f_pointer(info%elshell%expo, expo_ptr, [f_expo_len])
        call c_f_pointer(info%elshell%coef, coef_ptr, [f_expo_len])
        call c_f_pointer(info%elshell%ecp_rex, rexpo_ptr, [f_expo_len])
        call c_f_pointer(info%elshell%ecp_am, am_ptr, [info%elshell%ecp_nam])
        call c_f_pointer(info%elshell%ecp_coord, coord_ptr, [3*info%elshell%element_id])

        ecp_head%element_id = info%elshell%element_id
        ecp_head%n_angular_m = info%elshell%ecp_nam


        allocate(ecp_head%exponents(f_expo_len))
        allocate(ecp_head%coefficient(f_expo_len))
        allocate(ecp_head%ecp_r_expo(f_expo_len))
        allocate(ecp_head%ecp_am(info%elshell%ecp_nam))
        allocate(ecp_head%ecp_coord(3 * info%elshell%element_id))

        ecp_head%exponents = expo_ptr
        ecp_head%coefficient = coef_ptr
        ecp_head%ecp_r_expo = rexpo_ptr
        ecp_head%ecp_am = am_ptr
        ecp_head%ecp_coord = coord_ptr !* UNITS_ANGSTROM
    end subroutine oqp_append_ecp

    subroutine oqp_append_shell(info)
        use types, only: information
        type(information), intent(in) :: info
        type(electron_shell), pointer :: new_node, temp
        real(c_double), pointer :: expo_ptr(:), coef_ptr(:)
        integer(c_int), pointer ::n_expo_ptr(:)
        integer :: n_expo

        call c_f_pointer(info%elshell%num_expo, n_expo_ptr, [1])
        n_expo = n_expo_ptr(1)

        call c_f_pointer(info%elshell%expo, expo_ptr, [n_expo])
        call c_f_pointer(info%elshell%coef, coef_ptr, [n_expo])

        allocate(new_node)
        new_node%id = info%elshell%id
        new_node%element_id = info%elshell%element_id
        new_node%angular_momentum = info%elshell%ang_mom
        allocate(new_node%exponents(n_expo))
        allocate(new_node%coefficient(n_expo))
        allocate(new_node%n_exponents(1))
        new_node%n_exponents = n_expo_ptr
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

    subroutine map_shell2basis_set(infos, basis)
        use basis_tools, only: basis_set
        use types, only: information

        type(information), intent(inout) :: infos
        class(basis_set) ,intent(inout):: basis
        type(electron_shell), pointer :: temp
        type(electron_shell), pointer :: temp1
        integer :: nbf, nshell, nprim, mxcontr, mxam, ii
        integer :: n1,n2
        integer :: f_expo_len
        real, dimension(:), allocatable :: ex

        infos%control%basis_set_issue = .false.

        temp => head
        mxam = 0
        mxcontr = 0
        nbf = 0
        nshell = 0
        nprim = 0  ! Initialize nprim
        ii = 0


        do while (associated(temp))

            mxcontr = max(mxcontr, temp%n_exponents(1))
            mxam = max(mxam, temp%angular_momentum)

            nshell = temp%id
            nprim = nprim + temp%n_exponents(1)

            select case (temp%angular_momentum)
                case (0)
                    nbf = nbf + 1
                case (1)
                    nbf = nbf + 3
                case (2)
                    nbf = nbf + 6
                case (3)
                    nbf = nbf + 10
            end select
            temp => temp%next  ! Move to the next shell
        end do
        temp1 => head
        basis%mxam = mxam
        basis%mxcontr = mxcontr
        basis%nbf = nbf
        basis%nshell = nshell
        basis%nprim = nprim

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
            n2 = temp1%n_exponents(1)
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
                    basis%naos(ii) = 10
            end select

            n1 = temp1%n_exponents(1)

            temp1 => temp1%next

        end do
        if (.not. allocated(basis%ecp_zn_num)) allocate(basis%ecp_zn_num(maxval(basis%origin)))

        basis%ecp_zn_num = ecp_head%ecp_zn

        if (sum(basis%ecp_zn_num) > 0) then
            infos%mol_prop%nelec = infos%mol_prop%nelec - sum(basis%ecp_zn_num)
            infos%mol_prop%nelec_A = infos%mol_prop%nelec_A - sum(basis%ecp_zn_num)/2
            infos%mol_prop%nelec_B = infos%mol_prop%nelec_B - sum(basis%ecp_zn_num)/2
            infos%mol_prop%nocc = max(infos%mol_prop%nelec_A,infos%mol_prop%nelec_B)
        end if

        call basis%set_bfnorms()
        call basis%normalize_primitives()

        call head%clear()
        nullify(head)
        nullify(temp)
        nullify(temp1)


        if (ecp_head%element_id == 0) then
            basis%ecp_params%is_ecp = .false.
            return
        end if

        basis%ecp_params%is_ecp = .true.
        f_expo_len = sum(ecp_head%n_exponents)

        allocate(basis%ecp_params%ecp_ex(f_expo_len))
        allocate(basis%ecp_params%ecp_cc(f_expo_len))
        allocate(basis%ecp_params%ecp_coord(size(ecp_head%ecp_coord)))
        allocate(basis%ecp_params%ecp_r_ex(f_expo_len))
        allocate(basis%ecp_params%ecp_am(size(ecp_head%ecp_am)))
        allocate(basis%ecp_params%n_expo(size(ecp_head%n_exponents)))

        basis%ecp_params%ecp_ex = ecp_head%exponents
        basis%ecp_params%ecp_cc = ecp_head%coefficient
        basis%ecp_params%ecp_r_ex = ecp_head%ecp_r_expo
        basis%ecp_params%ecp_coord = ecp_head%ecp_coord
        basis%ecp_params%ecp_am = ecp_head%ecp_am
        basis%ecp_params%n_expo = ecp_head%n_exponents

        call ecp_head%clear()

    end subroutine map_shell2basis_set

    subroutine print_basis(infos)
        use types, only: information
        use elements, only: ELEMENTS_SHORT_NAME


        type(information), intent(inout) :: infos
        integer :: iw, i, j, atom, elem, end_i
        character(len=1) :: orbit

        open (newunit=iw, file=infos%log_filename, position="append")

        write(iw, '(/,5X,"====================== Basis Set Details ======================")')
        atom = 0

        do j = 1, infos%basis%nshell
            if (atom .NE. infos%basis%origin(j)) then
                elem = nint(infos%atoms%zn(infos%basis%origin(j)))
                write(iw, '(5X, A2)') ELEMENTS_SHORT_NAME(elem)
            end if
            select case (infos%basis%am(j))
                case (0)
                  orbit = 'S'
                case (1)
                  orbit = 'P'
                case (2)
                  orbit = 'D'
                case (3)
                  orbit = 'F'
                case default
                  orbit = '?'
            end select

            write(iw, '(10X, A1)') orbit

            end_i = infos%basis%g_offset(j) + infos%basis%ncontr(j) - 1

            do i = infos%basis%g_offset(j), end_i
                write(iw, '(15X, ES12.5, 15X, ES12.5)') infos%basis%ex(i),&
                        infos%basis%cc(i)
            end do
            atom = infos%basis%origin(j)
        end do
        if (infos%basis%ecp_params%is_ecp) then
            do i=1, infos%mol_prop%natom
                if (infos%basis%ecp_zn_num(i) .EQ. 0) cycle
                elem = nint(infos%atoms%zn(i))
                write(iw, '(5X, A2, A4)') ELEMENTS_SHORT_NAME(elem), '-ECP'
                write(iw, '(5X, A, I5)') 'Core Electrons Removed:', infos%basis%ecp_zn_num(i)
                call ecp_printing(infos%basis, iw, i)
            end do
        end if
        write(iw, '(/,5X,"==================== End of Basis Set Data ====================")')
        close(iw)

    end subroutine print_basis

    subroutine ecp_printing(basis,iw,j)

        use basis_tools, only: basis_set

        class(basis_set) ,intent(in):: basis
        integer, intent(in) :: iw, j
        integer :: i, start_i, end_i

        if (j > 1) then
            start_i = sum(basis%ecp_params%n_expo(1:j-1)) + 1
            end_i   = sum(basis%ecp_params%n_expo(1:j))
        else
            start_i = 1
            end_i   = basis%ecp_params%n_expo(j)
        end if
        do i = start_i, end_i
            write(iw, '(5X, I5, 5X, ES12.5, 5X, I5, 5X, ES12.5)') basis%ecp_params%ecp_am(i),&
                    basis%ecp_params%ecp_ex(i), basis%ecp_params%ecp_r_ex(i),&
                    basis%ecp_params%ecp_cc(i)
        end do

    end subroutine ecp_printing

    subroutine electron_shell_clear(this)
        class(electron_shell), intent(inout) :: this

        if (associated(this%n_exponents))  deallocate(this%n_exponents)
        if (associated(this%exponents))    deallocate(this%exponents)
        if (associated(this%coefficient))  deallocate(this%coefficient)

        this%angular_momentum = 0
        this%id         = 0
        this%element_id = 0

        if (associated(this%next)) then
            call this%next%clear()
            nullify(this%next)
        end if

    end subroutine electron_shell_clear

    subroutine ecpdata_clear(this)
        class(ecpdata), intent(inout) :: this

        if (associated(this%n_exponents))  deallocate(this%n_exponents)
        if (associated(this%exponents))    deallocate(this%exponents)
        if (associated(this%coefficient))  deallocate(this%coefficient)

        if (associated(this%ecp_zn))     deallocate(this%ecp_zn)
        if (associated(this%ecp_r_expo)) deallocate(this%ecp_r_expo)
        if (associated(this%ecp_am))     deallocate(this%ecp_am)
        if (associated(this%ecp_coord))  deallocate(this%ecp_coord)

        this%n_angular_m = 0
        this%id          = 0
        this%element_id  = 0
    end subroutine ecpdata_clear

    subroutine delete_all_shells()
        type(electron_shell), pointer :: to_delete

        do while (associated(head))
            to_delete => head
            head => head%next
            deallocate(to_delete%exponents)
            deallocate(to_delete%coefficient)
            deallocate(to_delete)
        end do
    end subroutine delete_all_shells

end module basis_api

