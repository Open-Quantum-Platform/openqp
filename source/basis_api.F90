!> @brief Basis ingestion API bridging external handles to OpenQP basis_set.
!> @detail Collects electron shells and ECP data from C/handles, builds the
!>         internal basis_set (cartesian AO layout), normalizes primitives,
!>         and prints a compact basis/ECP summary to the log.
!> @author Mohsen Mazaherifar
!> @date November 2025
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

    private
    public append_shell
    public append_ecp
    public map_shell2basis_set
    public print_basis 

contains
   !> @brief Append the current electron shell from an external handle.
   !> @detail Pulls one shell from `oqp_handle_get_info(c_handle)%elshell` and
   !>         pushes it to the internal linked list (`head`).
   !> @param[in] c_handle  Foreign handle carrying an `information` pointer.
   !> @note C binding: name="append_shell".
   !> @author Mohsen
   !> @date November 2025
    subroutine append_shell(c_handle) bind(C, name="append_shell")
        use c_interop, only: oqp_handle_t, oqp_handle_get_info
        use types, only: information
        type(oqp_handle_t) :: c_handle
        type(information), pointer :: inf
        inf => oqp_handle_get_info(c_handle)
        call oqp_append_shell(inf)
    end subroutine append_shell
    !> @brief Append ECP metadata from an external handle.
    !> @detail Copies global ECP arrays (per-element exponents/coeffs, AM, radii,
    !>         coords, and removed core Z) into `ecp_head` buffers.
    !> @param[in] c_handle  Foreign handle carrying an `information` pointer.
    !> @note C binding: name="append_ecp". No-op if element_id==0 (no ECP). 
    !> @author Mohsen
    !> @date November 2025
    subroutine append_ecp(c_handle) bind(C, name="append_ecp")
        use c_interop, only: oqp_handle_t, oqp_handle_get_info
        use types, only: information
        type(oqp_handle_t) :: c_handle
        type(information), pointer :: inf
        inf => oqp_handle_get_info(c_handle)
        call oqp_append_ecp(inf)
    end subroutine append_ecp
    !> @brief Internal: stage ECP data from `information%elshell` into `ecp_head`.
    !> @param[in] info  Read-only `information` snapshot with ECP fields populated.
    !> @author Mohsen
    !> @date November 2025
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
    !> @brief Internal: stage one electron shell from `information%elshell` to list.
    !> @param[in] info  Read-only `information` snapshot with shell fields populated.
    !> @note Maintains insertion order; computes n_exponents per shell. 
    !> @date November 2025
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
    !> @brief Build `basis_set` from staged shells/ECP and finalize normalization.
    !> @detail Computes sizes (nbf, nshell, nprim), allocates arrays, copies
    !>         exponents/coeffs, AO offsets, origins, AM; normalizes primitives,
    !>         sets ECP params (if present), and clears staging lists.
    !> @param[inout] infos  Provides target `basis_set` (primary or alt by flag).
    !> @note Sets `basis%ecp_params%is_ecp` and `basis%ecp_zn_num` when ECP present.
    !> @throws Sets `infos%control%basis_set_issue=.false.` on entry (reserved).
    !> @date November 2025
    subroutine map_shell2basis_set(infos)
        use basis_tools, only: basis_set
        use types, only: information
        use constants, only: NUM_CART_BF

        type(information), target, intent(inout) :: infos
        class(basis_set), pointer :: basis
        type(electron_shell), pointer :: temp
        type(electron_shell), pointer :: temp1
        integer :: nbf, nshell, nprim, mxcontr, mxam, ii
        integer :: n1,n2
        integer :: f_expo_len
        real, dimension(:), allocatable :: ex

        infos%control%basis_set_issue = .false.
        if (infos%control%active_basis == 0) then
            basis => infos%basis
        else
            basis => infos%alt_basis
        end if
        if (allocated(basis%ex)) call basis%destroy()

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

            nbf = nbf + NUM_CART_BF(temp%angular_momentum)

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
            basis%naos(ii) = NUM_CART_BF(temp1%angular_momentum)
            n1 = temp1%n_exponents(1)

            temp1 => temp1%next

        end do
        if (.not. allocated(basis%ecp_zn_num)) allocate(basis%ecp_zn_num(maxval(basis%origin)))

        basis%ecp_zn_num = ecp_head%ecp_zn


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

        if (.not. allocated(basis%ecp_params%ecp_ex)) allocate(basis%ecp_params%ecp_ex(f_expo_len))
        if (.not. allocated(basis%ecp_params%ecp_cc)) allocate(basis%ecp_params%ecp_cc(f_expo_len))
        if (.not. allocated(basis%ecp_params%ecp_coord)) allocate(basis%ecp_params%ecp_coord(size(ecp_head%ecp_coord)))
        if (.not. allocated(basis%ecp_params%ecp_r_ex)) allocate(basis%ecp_params%ecp_r_ex(f_expo_len))
        if (.not. allocated(basis%ecp_params%ecp_am)) allocate(basis%ecp_params%ecp_am(size(ecp_head%ecp_am)))
        if (.not. allocated(basis%ecp_params%n_expo)) allocate(basis%ecp_params%n_expo(size(ecp_head%n_exponents)))


        basis%ecp_params%ecp_ex = ecp_head%exponents
        basis%ecp_params%ecp_cc = ecp_head%coefficient
        basis%ecp_params%ecp_r_ex = ecp_head%ecp_r_expo
        basis%ecp_params%ecp_coord = ecp_head%ecp_coord
        basis%ecp_params%ecp_am = ecp_head%ecp_am
        basis%ecp_params%n_expo = ecp_head%n_exponents

        call ecp_head%clear()

    end subroutine map_shell2basis_set
    !> @brief Pretty-print the active basis (and ECP, if any) to the log file.
    !> @detail Lists element, shell AM, primitives with normalized coefficients;
    !>         for ECP atoms, prints removed core Z and per-term parameters.
    !> @param[inout] infos  Supplies basis, atom labels, and output filename.
    !> @date November 2025
    subroutine print_basis(infos)
        use types, only: information
        use elements, only: ELEMENTS_SHORT_NAME
        use basis_tools, only: basis_set

        type(information), target, intent(inout) :: infos
        type(basis_set), pointer :: basis

        integer :: iw, i, j, atom, elem, end_i
        character(len=1) :: orbit
        character(len=1), dimension(5) :: orbital_types = ['S', 'P', 'D', 'F', 'G']
        if (infos%control%active_basis == 0) then
           basis => infos%basis
        else
           basis => infos%alt_basis
        end if
        open (newunit=iw, file=infos%log_filename, position="append")

        write(iw, '(/,5X,"====================== Basis Set Details ======================")')
        atom = 0

        write(iw, '(/,17X, A,13X,A)') 'Exponent', 'Normalized Coefficient'
        do j = 1, basis%nshell
            if (atom .NE. basis%origin(j)) then
                elem = nint(infos%atoms%zn(infos%basis%origin(j)))
                write(iw, '(5X, A2)') ELEMENTS_SHORT_NAME(elem)
            end if
            orbit = orbital_types(basis%am(j)+1)

            write(iw, '(10X, A1)') orbit

            end_i = basis%g_offset(j) + basis%ncontr(j) - 1

            do i = basis%g_offset(j), end_i
                write(iw, '(15X, ES12.5, 15X, ES12.5)') basis%ex(i),&
                        basis%cc(i)
            end do
            atom = basis%origin(j)
        end do
        if (basis%ecp_params%is_ecp) then
            do i=1, infos%mol_prop%natom
                if (basis%ecp_zn_num(i) .EQ. 0) cycle
                elem = nint(infos%atoms%zn(i))
                write(iw, '(5X, A2, A4)') ELEMENTS_SHORT_NAME(elem), '-ECP'
                write(iw, '(5X, A, I5)') 'Core Electrons Removed:', basis%ecp_zn_num(i)
                call ecp_printing(basis, iw, i)
            end do
        end if
        write(iw, '(/,5X,"==================== End of Basis Set Data ====================")')
        close(iw)

    end subroutine print_basis
    !> @brief Print ECP terms for atom j to the open unit `iw`.
    !> @param[in] basis  Basis with `ecp_params` populated.
    !> @param[in] iw     Fortran unit already opened for append.
    !> @param[in] j      Atom index (1-based).
    !> @date November 2025
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
                    basis%ecp_params%ecp_cc(i), basis%ecp_params%ecp_r_ex(i),&
                    basis%ecp_params%ecp_ex(i)
        end do

    end subroutine ecp_printing
    !> @brief Clear an electron_shell node and recursively clear its `next` chain.
    !> @date November 2025
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
    !> @brief Clear staged ECP buffers in `ecpdata`.
    !> @date November 2025
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

end module basis_api

