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

   integer, parameter :: BAS_MXANG = 6
   integer, parameter :: BAS_MXCONTR = 30
   integer, parameter :: BAS_MXCART = (BAS_MXANG+1)*(BAS_MXANG+2)/2
   integer, parameter :: NUM_CART_BF(0:BAS_MXANG) = [((i+1)*(i+2)/2, i = 0, BAS_MXANG)]
   !< powers of X,Y,Z in Cartesian Gaussian basis functions
   integer, parameter :: &
    CART_X(BAS_MXCART,0:BAS_MXANG) = reshape([ &
       [0,                                                               (0, i = NUM_CART_BF(0)+1, BAS_MXCART)], &
       [1, 0, 0,                                                         (0, i = NUM_CART_BF(1)+1, BAS_MXCART)], &
       [2, 0, 0, 1, 1, 0,                                                (0, i = NUM_CART_BF(2)+1, BAS_MXCART)], &
       [3, 0, 0, 2, 2, 1, 0, 1, 0, 1,                                    (0, i = NUM_CART_BF(3)+1, BAS_MXCART)], &
       [4, 0, 0, 3, 3, 1, 0, 1, 0, 2, 2, 0, 2, 1, 1,                     (0, i = NUM_CART_BF(4)+1, BAS_MXCART)], &
       [5, 0, 0, 4, 4, 1, 0, 1, 0, 3, 3, 2, 0, 2, 0, 3, 1, 1, 2, 2, 1,   (0, i = NUM_CART_BF(5)+1, BAS_MXCART)], &
       [6, 0, 0, 5, 5, 1, 0, 1, 0, 4, 4, 2, 0, 2, 0, 4, 1, 1, 3, 3, 0, 3, 3, 2, 1, 2, 1, 2]  &
       ], shape(CART_X))
   integer, parameter :: &
    CART_Y(BAS_MXCART,0:BAS_MXANG) = reshape([ &
       [0,                                                              (0, i = NUM_CART_BF(0)+1, BAS_MXCART)], &
       [0, 1, 0,                                                        (0, i = NUM_CART_BF(1)+1, BAS_MXCART)], &
       [0, 2, 0, 1, 0, 1,                                               (0, i = NUM_CART_BF(2)+1, BAS_MXCART)], &
       [0, 3, 0, 1, 0, 2, 2, 0, 1, 1,                                   (0, i = NUM_CART_BF(3)+1, BAS_MXCART)], &
       [0, 4, 0, 1, 0, 3, 3, 0, 1, 2, 0, 2, 1, 2, 1,                    (0, i = NUM_CART_BF(4)+1, BAS_MXCART)], &
       [0, 5, 0, 1, 0, 4, 4, 0, 1, 2, 0, 3, 3, 0, 2, 1, 3, 1, 2, 1, 2,  (0, i = NUM_CART_BF(5)+1, BAS_MXCART)], &
       [0, 6, 0, 1, 0, 5, 5, 0, 1, 2, 0, 4, 4, 0, 2, 1, 4, 1, 3, 0, 3, 2, 1, 3, 3, 1, 2, 2] &
       ], shape(CART_Y))
   integer, parameter :: &
    CART_Z(BAS_MXCART,0:BAS_MXANG) = reshape([ &
       [0,                                                              (0, i = NUM_CART_BF(0)+1, BAS_MXCART)], &
       [0, 0, 1,                                                        (0, i = NUM_CART_BF(1)+1, BAS_MXCART)], &
       [0, 0, 2, 0, 1, 1,                                               (0, i = NUM_CART_BF(2)+1, BAS_MXCART)], &
       [0, 0, 3, 0, 1, 0, 1, 2, 2, 1,                                   (0, i = NUM_CART_BF(3)+1, BAS_MXCART)], &
       [0, 0, 4, 0, 1, 0, 1, 3, 3, 0, 2, 2, 1, 1, 2,                    (0, i = NUM_CART_BF(4)+1, BAS_MXCART)], &
       [0, 0, 5, 0, 1, 0, 1, 4, 4, 0, 2, 0, 2, 3, 3, 1, 1, 3, 1, 2, 2,  (0, i = NUM_CART_BF(5)+1, BAS_MXCART)], &
       [0, 0, 6, 0, 1, 0, 1, 5, 5, 0, 2, 0, 2, 4, 4, 1, 1, 4, 0, 3, 3, 1, 2, 1, 2, 3, 3, 2] &
       ], shape(CART_Z))

!  ao symbol
   integer, private :: iii
   character(len=4), parameter :: bf_names(15,0:5) = reshape([&
                ['  S ', &
                ('    ', iii = num_cart_bf(0)+1, 15)], &
                ['  X ','  Y ','  Z ', &
                ('    ', iii = num_cart_bf(1)+1, 15)], &
                [' XX ',' YY ',' ZZ ',' XY ',' XZ ',' YZ ', &
                ('    ', iii = num_cart_bf(2)+1, 15)], &
                [' XXX',' YYY',' ZZZ',' XXY',' XXZ', &
                 ' YYX',' YYZ',' ZZX',' ZZY',' XYZ', &
                 ('    ', iii = num_cart_bf(3)+1, 15)], &
                ['XXXX','YYYY','ZZZZ','XXXY','XXXZ', &
                 'YYYX','YYYZ','ZZZX','ZZZY','XXYY', &
                 'XXZZ','YYZZ','XXYZ','YYXZ','ZZXY', &
                 ('    ', iii = num_cart_bf(4)+1, 15)], &
                ['????', &
                ('    ', iii = 2, 15)] &
             ], shape(bf_names))

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
