module libecp_result
    use iso_c_binding
    implicit none

    type, bind(c) :: ecp_result
        type(c_ptr) :: data
        integer(c_int) :: size
    end type ecp_result
end module libecp_result

module libecpint_wrapper
    use iso_c_binding, only : c_int, c_double, c_ptr
!    use ecpresult, only : ecp_result
    implicit none

    interface
        function init_integrator(num_gaussians, g_coords, g_exps, g_coefs, &
                                 g_ams, g_lengths) bind(c, name="init_integrator")
            use iso_c_binding, only : c_int, c_double, c_ptr
            type(c_ptr) :: init_integrator
            integer(c_int), value :: num_gaussians
            real(c_double), dimension(*), intent(in) :: g_coords
            real(c_double), dimension(*), intent(in) :: g_exps
            real(c_double), dimension(*), intent(in) :: g_coefs
!            real(c_double), dimension(*) :: g_coords, g_exps, g_coefs
            integer(c_int), dimension(*) :: g_ams, g_lengths
        end function init_integrator

        subroutine set_ecp_basis(integrator, num_ecps, u_coords, u_exps, &
                                 u_coefs, u_ams, u_ns, u_lengths) &
                                 bind(c, name="set_ecp_basis")
            use iso_c_binding, only : c_int, c_double, c_ptr
            type(c_ptr), value :: integrator
            integer(c_int), value :: num_ecps
            real(c_double), dimension(*) :: u_coords, u_exps, u_coefs
            integer(c_int), dimension(*) :: u_ams, u_ns, u_lengths
        end subroutine set_ecp_basis

        subroutine init_integrator_instance(integrator) &
                   bind(c, name="init_integrator_instance")
            use iso_c_binding, only : c_ptr
            type(c_ptr), value :: integrator
        end subroutine init_integrator_instance

        function compute_integrals(integrator) bind(c, name="compute_integrals")
            use iso_c_binding, only : c_ptr
            use libecp_result
            type(ecp_result) :: compute_integrals
            type(c_ptr), value :: integrator
        end function compute_integrals

        subroutine free_integrator(integrator) bind(c, name="free_integrator")
            use iso_c_binding, only : c_ptr
            type(c_ptr), value :: integrator
        end subroutine free_integrator

        subroutine free_result(result) bind(c, name="free_result")
            use libecp_result
            type(ecp_result), value :: result
        end subroutine free_result
    end interface
end module libecpint_wrapper
