SUBROUTINE Merge_D(N,IMid,IA,IB,A,B)
    IMPLICIT NONE
    INTEGER, INTENT(IN)             :: N
    INTEGER, INTENT(IN)             :: IMid
    INTEGER, INTENT(INOUT)          :: IA(1:N)
    INTEGER, INTENT(OUT)            :: IB(1:N)
    DOUBLE PRECISION, INTENT(INOUT) :: A(1:N)
    DOUBLE PRECISION, INTENT(OUT)   :: B(1:N)
    INTEGER                         :: I, J, K
    I = 1
    J = IMid+1 
    K = 1
    DO
       IF (I > IMid .OR. J > N) EXIT
       IF (A(I) <= A(J)) THEN
          B(K)  = A(I)
          IB(K) = IA(I)
          I = I + 1
       ELSE
          B(K)  = A(J)
          IB(K) = IA(J)
          J = J + 1
       ENDIF
       K = K + 1
    END DO
    IF (I > IMid) THEN
       B(K:)  = A(J:N)
       IB(K:) = IA(J:N)
    ELSE
       B(K:)  = A(I:IMid)
       IB(K:) = IA(I:IMid)
    ENDIF
    A(1:N)  = B(1:N)
    IA(1:N) = IB(1:N)
END SUBROUTINE

RECURSIVE SUBROUTINE Merge_Sort_D(N,IA,IB,A,B) 
    IMPLICIT NONE
    INTEGER, INTENT(IN)             :: N
    INTEGER, INTENT(INOUT)          :: IA(1:N)
    INTEGER, INTENT(OUT)            :: IB(1:N)
    DOUBLE PRECISION, INTENT(INOUT) :: A(1:N)
    DOUBLE PRECISION, INTENT(OUT)   :: B(1:N)
    INTEGER                         :: IMid
    INTEGER                         :: NL
    INTEGER                         :: NR
    IF (N <= 1) RETURN
    IMid  = N/2 
    NL    = IMid 
    NR    = N - NL
    CALL Merge_Sort_D(NL,IA(1),IB(1),A(1),B(1))
    CALL Merge_Sort_D(NR,IA(IMid+1),IB(IMid+1),A(IMid+1),B(IMid+1))
    CALL Merge_D(N,IMid,IA,IB,A,B)
END SUBROUTINE

!PROGRAM Main
!IMPLICIT NONE
!  INTEGER :: i
!  INTEGER, PARAMETER :: n=10
!  DOUBLE PRECISION   :: A(1:n) 
!  DOUBLE PRECISION   :: B(1:n) 
!  INTEGER  :: IA(1:n)
!  INTEGER  :: IB(1:n)
!  A(1:n) = [8D+0,4D+0,7D+0,2D+0,1D+0,3D+0,5D+0,6D+0,9D+0,10D+0]
!  B = 0.D0
!  IA(1:n) = [(i,i=1,n)]
!  WRITE(6,*)A
!  CALL Merge_Sort_D(n, IA, IB, A, B)
!  WRITE(6,*)A
!  WRITE(6,*)IA
!END PROGRAM
