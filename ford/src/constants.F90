module constants

   use precision, only: dp

   implicit none

   real(kind=dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
   integer, parameter :: tol_int = 20

!  atomic number
   integer, private :: i

!  angular momentum labels
   character(len=7) :: angular_label = 'SPDFGHI'

!  Boltzmann Constant
   real(kind=dp), parameter :: kB_HaK = 3.166811563e-6_dp

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

!  canonical order is achieved by ordering angular momentum components (x,y,z)
!  in descending order. For a given total angular momentum L, the components
!  are generated as follows:
!  Example for L = 2:
!     do x = L, 0, -1          ! x descends from L to 0
!         do y = L-x, 0, -1    ! y descends from remaining momentum (L-x) to 0
!             z = L - x - y     ! z takes the remaining momentum
!             ! This generates ordered triplets (x,y,z) where x >= y >= z
!             ! and x + y + z = L
!         end do
!     end do 
!  Last Modified: 2025-02-03 06:00:22 UTC

   integer, parameter :: &
      map_canonical(BAS_MXCART,0:BAS_MXANG) = reshape([ &
      [0,                                                               (0, i = NUM_CART_BF(0)+1, BAS_MXCART)], & ! l = 0 (S)
      [0, 0, 0,                                                         (0, i = NUM_CART_BF(1)+1, BAS_MXCART)], & ! l = 1 (P)
      [0, 2, 3, -2, -2, -1,                                             (0, i = NUM_CART_BF(2)+1, BAS_MXCART)], & ! l = 2 (D)
      [0, 5, 7, -2, -2, -2, 1, -2, 0, -5,                               (0, i = NUM_CART_BF(3)+1, BAS_MXCART)], & ! l = 3 (F)
      [0, 9, 12, -2, -2, 1, 5, 2, 5, -6, -5, 1, -8, -6, -6,             (0, i = NUM_CART_BF(4)+1, BAS_MXCART)], & ! l = 4 (G)
      [0,14,18,-2,-2,5,10,7,11,-6,-5,-5,5,-4,4,-11,-5,-4,-11,-11,-8,    (0, i = NUM_CART_BF(5)+1, BAS_MXCART)], & ! l = 4 (H)
      [0,20,25,-2,-2,10,16,13,18,-6,-5,-1,11,1,11,-11,0,2,-12,-10,4,-14,-14,-12,-7,-12,-8,-15] &
      ],  shape(map_canonical))


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
