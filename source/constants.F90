module constants

   use precision, only: dp

   implicit none

   real(kind=dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
   integer, parameter :: tol_int = 20

!  atomic number
   integer, private :: i

!  angular momentum labels
   character(len=7) :: angular_label = 'SPDFGHI'

!  number of cartesian bf for each shell kind
   integer, parameter :: num_cart_bf(7) = [(i*(i+1)/2, i = 1, 7)]

!  ao symbol
   character(len=4), parameter :: bfnam(36)=(/'  S ','  X ','  Y ','  Z ', &
                                 ' XX ',' YY ',' ZZ ',' XY ',' XZ ',' YZ ', &
                                 ' XXX',' YYY',' ZZZ',' XXY',' XXZ', &
                                 ' YYX',' YYZ',' ZZX',' ZZY',' XYZ', &
                                 'XXXX','YYYY','ZZZZ','XXXY','XXXZ', &
                                 'YYYX','YYYZ','ZZZX','ZZZY','XXYY', &
                                 'XXZZ','YYZZ','XXYZ','YYXZ','ZZXY', '????'/)

! normalization constants
  real(kind=dp), target, save :: shells_pnrm2(28,0:6) = reshape([ &
    [1.0_dp, &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0],  & ! s-shell

    [1.0_dp, 1.0_dp, 1.0_dp, &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], & ! p-shell

    [1.0_dp, 1.0_dp, 1.0_dp, sqrt(3.0_dp), sqrt(3.0_dp), sqrt(3.0_dp), &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], & ! d-shell

    [1.0_dp, 1.0_dp, 1.0_dp, sqrt(5.0_dp), sqrt(5.0_dp), sqrt(5.0_dp), sqrt(5.0_dp), sqrt(5.0_dp), &
     sqrt(5.0_dp), sqrt(5.0_dp)*sqrt(3.0_dp), &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
     0.0d0, 0.0d0, 0.0d0], & ! f-shell

    [1.0_dp, 1.0_dp, 1.0_dp, sqrt(7.0_dp), sqrt(7.0_dp), sqrt(7.0_dp), sqrt(7.0_dp), &
     sqrt(7.0_dp), sqrt(7.0_dp), sqrt(7.0_dp)*sqrt(5.0_dp)/sqrt(3.0_dp), &
     sqrt(7.0_dp)*sqrt(5.0_dp)/sqrt(3.0_dp), sqrt(7.0_dp)*sqrt(5.0_dp)/sqrt(3.0_dp), &
     sqrt(7.0_dp)*sqrt(5.0_dp), sqrt(7.0_dp)*sqrt(5.0_dp), sqrt(7.0_dp)*sqrt(5.0_dp), &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], & ! g-shell

    [1.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, sqrt(7.0_dp)*sqrt(3.0_dp), &
     sqrt(7.0_dp)*sqrt(3.0_dp), sqrt(7.0_dp)*sqrt(3.0_dp), sqrt(7.0_dp)*sqrt(3.0_dp), &
     sqrt(7.0_dp)*sqrt(3.0_dp), sqrt(7.0_dp)*sqrt(3.0_dp), sqrt(7.0_dp)*3.0_dp, &
     sqrt(7.0_dp)*3.0_dp, sqrt(7.0_dp)*3.0_dp, sqrt(3.0_dp)*sqrt(5.0_dp)*sqrt(7.0_dp), &
     sqrt(3.0_dp)*sqrt(5.0_dp)*sqrt(7.0_dp), sqrt(3.0_dp)*sqrt(5.0_dp)*sqrt(7.0_dp), &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], & ! h-shell

    [1.0_dp, 1.0_dp, 1.0_dp, sqrt(11.0_dp), sqrt(11.0_dp), sqrt(11.0_dp), sqrt(11.0_dp), &
     sqrt(11.0_dp), sqrt(11.0_dp), sqrt(3.0_dp)*sqrt(11.0_dp), &
     sqrt(3.0_dp)*sqrt(11.0_dp), sqrt(3.0_dp)*sqrt(11.0_dp), sqrt(3.0_dp)*sqrt(11.0_dp), &
     sqrt(3.0_dp)*sqrt(11.0_dp), sqrt(3.0_dp)*sqrt(11.0_dp), &
     3.0_dp*sqrt(11.0_dp), 3.0_dp*sqrt(11.0_dp), 3.0_dp*sqrt(11.0_dp), &
     sqrt(3.0_dp)*sqrt(7.0_dp)*sqrt(11.0_dp)/sqrt(5.0_dp), &
     sqrt(3.0_dp)*sqrt(7.0_dp)*sqrt(11.0_dp)/sqrt(5.0_dp), &
     sqrt(3.0_dp)*sqrt(7.0_dp)*sqrt(11.0_dp)/sqrt(5.0_dp), sqrt(3.0_dp)*sqrt(7.0_dp)*sqrt(11.0_dp), &
     sqrt(3.0_dp)*sqrt(7.0_dp)*sqrt(11.0_dp), sqrt(3.0_dp)*sqrt(7.0_dp)*sqrt(11.0_dp), &
     sqrt(3.0_dp)*sqrt(7.0_dp)*sqrt(11.0_dp), sqrt(3.0_dp)*sqrt(7.0_dp)*sqrt(11.0_dp), &
     sqrt(3.0_dp)*sqrt(7.0_dp)*sqrt(11.0_dp), sqrt(5.0_dp)*sqrt(7.0_dp)*sqrt(11.0_dp)] &! i-shell
    ], shape = [28,7])

  contains

end module constants
