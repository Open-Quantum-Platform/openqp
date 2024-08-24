! MODULE MX_LIMITS
!>    @author  Vladimir Mironov
!
!>    @brief   Contains parameters scattered throughout all
!>             of the code that define sizes of static
!>             memory arrays
!
!     REVISION HISTORY:
!>    @date _Jan, 2017_ Initial release
!
module openqp_config

  implicit none

  integer, parameter :: &
    basis_max_angular_momentum = 7, &
    basis_max_contraction = 30         !< Max. degree of contraction

end module openqp_config
