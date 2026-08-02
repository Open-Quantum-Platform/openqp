MODULE AFQMC_Random_Normal
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: AFQMC_init_random_seed, AFQMC_rand_normal, AFQMC_rand_normal_array

! 数学的背景
! 一様乱数 u1, u2 ∈(0,1) から
! z = \sqrt{-2 \ln u_1}\cos(2\pi u_2)
! を作ると、z ∼ N(0,1)  
CONTAINS

  !--------------------------------------------
  ! 乱数シード初期化（時刻ベース）
  !--------------------------------------------
  SUBROUTINE AFQMC_init_random_seed(ISeed)
    INTEGER, INTENT(IN) :: ISeed
    INTEGER :: n, i
    INTEGER, ALLOCATABLE :: seed(:)
    ! シードの格納に必要なサイズを取得する
    CALL Random_Seed(SIZE = n)
    ALLOCATE(seed(n))
    !call system_clock(count = clock)
    !seed = clock + 37 * [(i - 1, i = 1, n)]
    !seed = clock + 37 * [(i - 1, i = 1, n)]
    !CALL Random_Seed(get = seed)
    ! シード値を設定
    !seed = seed + ISeed
    ! Fill the complete implementation-defined seed vector.  Assigning the
    ! same integer to every element produces unnecessarily correlated streams
    ! on some compilers, especially for neighboring MPI ranks.
    DO i=1,n
       seed(i) = MODULO(ABS(ISeed) + 104729*i + 7919*i*i, HUGE(1)-1)
       IF (seed(i) == 0) seed(i) = i
    ENDDO
    !write(6,'("seeds=",10I20)')seed(1:n)
    CALL Random_Seed(put = seed)
    DEALLOCATE(seed)
  END SUBROUTINE 

  !--------------------------------------------
  ! 正規分布乱数 N(mu, sigma^2)
  !--------------------------------------------
  DOUBLE PRECISION FUNCTION AFQMC_rand_normal(mu, sigma) RESULT(randn)
  !
  !  Box--Muller 
  !    x = 
  !
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: mu, sigma
  DOUBLE PRECISION             :: u1, u2
  DOUBLE PRECISION, PARAMETER  :: Pi=3.141592653589793238D0
  DOUBLE PRECISION, PARAMETER  :: Zero=0D+00
  DOUBLE PRECISION, PARAMETER  :: Small=TINY(1D+00)
  DOUBLE PRECISION, PARAMETER  :: Two=2D+00
  CALL random_number(u1)
  CALL random_number(u2)

  ! log(0) 回避
  !IF (u1 <= 0.0D0) u1 = TINY(1.D0)
  u1    = MAX(u1,Small)
  randn = mu + sigma * SQRT(-Two * LOG(u1)) &
                     * COS(Two * Pi * u2)

  !randn = mu + sigma * SQRT(-2.D0 * LOG(u1)) &
  !               * COS(2.0D0 * ACOS(-1.D0) * u2)
  END FUNCTION

  !--------------------------------------------
  ! OpenQP port: block version of AFQMC_rand_normal.
  !
  ! The auxiliary fields are drawn one scalar at a time before the threaded
  ! walker loop, which is the single largest serial cost per step. This draws
  ! the whole block with one RANDOM_NUMBER call and vectorises the Box-Muller
  ! transform.
  !
  ! It consumes the uniform stream in exactly the same order as N successive
  ! AFQMC_rand_normal calls (u1 then u2, per variate) and applies the identical
  ! expression, so the numbers produced are bit-identical to the scalar path.
  !--------------------------------------------
  SUBROUTINE AFQMC_rand_normal_array(mu, sigma, n, z)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)  :: mu, sigma
    INTEGER,          INTENT(IN)  :: n
    DOUBLE PRECISION, INTENT(OUT) :: z(1:n)
    DOUBLE PRECISION, PARAMETER   :: Pi=3.141592653589793238D0
    DOUBLE PRECISION, PARAMETER   :: Small=TINY(1D+00)
    DOUBLE PRECISION, PARAMETER   :: Two=2D+00
    DOUBLE PRECISION, ALLOCATABLE :: u(:)
    INTEGER :: i
    IF (n <= 0) RETURN
    ALLOCATE(u(1:2*n))
    CALL random_number(u)
    DO i = 1, n
       z(i) = mu + sigma * SQRT(-Two * LOG(MAX(u(2*i-1),Small))) &
                         * COS(Two * Pi * u(2*i))
    ENDDO
    DEALLOCATE(u)
  END SUBROUTINE

END MODULE


!program test_normal_random
!  use normal_random
!  implicit none
!
!  integer :: i
!  real(8) :: x
!
!  call init_random_seed()
!
!  do i = 1, 1000000
!     x = randn(0.0d0, 1.0d0)   ! 平均0、標準偏差1
!     write(*,'(F10.6)') x
!  end do
!end program test_normal_random
