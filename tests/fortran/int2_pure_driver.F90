program int2_pure_driver
  implicit none
  interface
    subroutine OQP_ERI_SELFTEST() bind(C)
    end subroutine OQP_ERI_SELFTEST
  end interface
  call OQP_ERI_SELFTEST()
end program int2_pure_driver
