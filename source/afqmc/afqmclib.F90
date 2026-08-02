SUBROUTINE AFQMC_TFT(M,N,ShapeA,A,LdA,T,LdT,B,LdB)
   !
   !
   !   B(1:M,1:M) = TRANSPOSE(T(1:N,1:M)).A(1:N,1:N).T(1:N,1:M)
   !   LdB > M
   !   LdT > N
   !   LdA > N
   !
   IMPLICIT NONE
   DOUBLE PRECISION, PARAMETER   :: Zero=0D+00, One=1D+00
   INTEGER,          INTENT(IN)  :: M
   INTEGER,          INTENT(IN)  :: N
   CHARACTER(LEN=*), INTENT(IN)  :: ShapeA
   INTEGER,          INTENT(IN)  :: LdA
   DOUBLE PRECISION, INTENT(IN)  :: A(1:*)
   INTEGER,          INTENT(IN)  :: LdT
   DOUBLE PRECISION, INTENT(IN)  :: T(1:LdT,1:M)
   INTEGER,          INTENT(IN)  :: LdB
   DOUBLE PRECISION, INTENT(OUT) :: B(1:*)
   INTEGER :: LUse, Last, Need
   INTEGER :: LWrk
   IF (LdB < M) ERROR STOP
   IF (LdT < N) ERROR STOP
   IF (LdA < N) ERROR STOP

   SELECT CASE (ShapeA(1:1))
   CASE ("T","t")
      CALL AFQMC_TFT_TRI(M,N,A,LdA,T,LdT,B,LdB)
   CASE ("S","s")
      CALL AFQMC_TFT_SQR(M,N,A,LdA,T,LdT,B,LdB)
   END SELECT
END SUBROUTINE

SUBROUTINE AFQMC_TFT_SQR(M,N,A,LdA,T,LdT,B,LdB)
   !
   !
   !   B(1:M,1:M) = TRANSPOSE(T(1:N,1:M)).A(1:N,1:N).T(1:N,1:M)
   !   LdB > M
   !   LdT > N
   !   LdA > N
   !
   IMPLICIT NONE
   DOUBLE PRECISION, PARAMETER   :: Zero=0D+00, One=1D+00
   INTEGER,          INTENT(IN)  :: M
   INTEGER,          INTENT(IN)  :: N
   INTEGER,          INTENT(IN)  :: LdA
   DOUBLE PRECISION, INTENT(IN)  :: A(1:LdA,1:N)
   INTEGER,          INTENT(IN)  :: LdT
   DOUBLE PRECISION, INTENT(IN)  :: T(1:LdT,1:M)
   INTEGER,          INTENT(IN)  :: LdB
   DOUBLE PRECISION, INTENT(OUT) :: B(1:LdB,1:M)
   ! OpenQP port: the Gellan global scratch (Core_Memory::X) is replaced by a
   ! local allocatable.  X is a fixed 10^6-word static array shared by every
   ! thread, so it both caps the basis size and races once the walker loop is
   ! threaded.
   DOUBLE PRECISION, ALLOCATABLE :: Wrk(:,:)
   ALLOCATE(Wrk(1:M,1:N))
   CALL DGeMM('T','N',M,N,N,One,T,LdT,A,LdA,Zero,Wrk,M)
   CALL DGeMM('N','N',M,M,N,One,Wrk,M,T,LdT,Zero,B,  M)
   DEALLOCATE(Wrk)
END SUBROUTINE
 
SUBROUTINE AFQMC_TFT_TRI(M,N,A,LdA,T,LdT,B,LdB)
   !
   !
   !   B(1:M,1:M) = TRANSPOSE(T(1:N,1:M)).A(1:N,1:N).T(1:N,1:M)
   !   LdB > M
   !   LdT > N
   !   LdA > N
   !
   IMPLICIT NONE
   DOUBLE PRECISION, PARAMETER   :: Zero=0D+00, One=1D+00
   INTEGER,          INTENT(IN)  :: M
   INTEGER,          INTENT(IN)  :: N
   INTEGER,          INTENT(IN)  :: LdA
   DOUBLE PRECISION, INTENT(IN)  :: A(1:*)
   INTEGER,          INTENT(IN)  :: LdT
   DOUBLE PRECISION, INTENT(IN)  :: T(1:LdT,1:M)
   INTEGER,          INTENT(IN)  :: LdB
   DOUBLE PRECISION, INTENT(OUT) :: B(1:*)
   ! OpenQP port: local allocatables in place of the Gellan global scratch.
   DOUBLE PRECISION, ALLOCATABLE :: AFull(:,:), Wrk(:,:), BFull(:,:)
   ALLOCATE(AFull(1:N,1:N), Wrk(1:M,1:N), BFull(1:M,1:M))
   CALL Expnd(A,AFull,N,0)
   CALL DGeMM('T','N',M,N,N,One,T,LdT,AFull,N,Zero,Wrk,M)
   CALL DGeMM('N','N',M,M,N,One,Wrk,M,T,LdT,Zero,BFull,M)
   ! BFull holds the transformed M-by-M matrix.  Using N here can overrun both
   ! the workspace and the packed output when the AO and retained-MO
   ! dimensions differ.
   CALL Cmprss(BFull,B,M)
   DEALLOCATE(AFull,Wrk,BFull)
END SUBROUTINE

SUBROUTINE Cmprss(A,B,N)
   IMPLICIT NONE
   INTEGER,          INTENT(IN)  :: N
   DOUBLE PRECISION, INTENT(IN)  :: A(1:N,1:N)
   DOUBLE PRECISION, INTENT(OUT) :: B(1:*)
   INTEGER  :: I, J, IJ
   IJ = 0
   DO J = 1,N
      DO I = 1,J
         IJ = IJ + 1
         B(IJ) = (A(I,J) + A(J,I))/2.D0
      ENDDO
   ENDDO
END SUBROUTINE

SUBROUTINE AFQMC_Greens_Functions_D(NVar,NSpin,NOccA,NOccB,U,V,Det,G,G_Half)
   !
   !  G_pq = <Psi_T|ap+.aq|Psi_W>/<Psi_T|Psi_W> = [CW.(CT+.CW)-1.CT+]_qp
   !
   !           |Psi_T> = U|Psi>  (U=CT)
   !           |Psi_W> = V|Psi>  (V=CW)
   !     であるとすれば
   !        <Psi_T|ap+.aq|Psi_W>/<Psi_T|Psi_W> 
   !        = <Psi|U+.ap+.aq.V|Psi> / <Psi|U+.V|Psi> 
   !        = <Psi|[U+.a_p+.U].U+.V.[V+.a_q.V]|Psi> / <Psi|U+.V|Psi> 
   !        = <Psi|U*_pr.a_r+.[U+.V].a_s.V_qs|Psi> / <Psi|U+.V|Psi>
   !        = U*_pr V_qs <Psi|a_r+.[U+.V].a_s|Psi> / <Psi|U+.V|Psi>
   !        = U*_pr V_qs [U+.V]^{-1}_sr / 
   !        = [V.(U+.V)^{-1}.U+]_qp
   !        = [U.(V+.U)^{-1}.V+]*_pq
   !
   !  p and q are occupied orbitals
   ! 
   !  ここで CTは実数、CWは複素数を想定しないといけない
   !
   USE IO_Files,     ONLY: IW
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   INTEGER,          INTENT(IN)  :: NVar
   INTEGER,          INTENT(IN)  :: NSpin
   INTEGER,          INTENT(IN)  :: NOccA
   INTEGER,          INTENT(IN)  :: NOccB
#if 1
   DOUBLE PRECISION, INTENT(IN)  :: U(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(IN)  :: V(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(OUT) :: Det
   DOUBLE PRECISION, INTENT(OUT) :: G(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(OUT) :: G_Half(1:NVar,1:NVar,1:NSpin)
#else
   DOUBLE COMPLEX,   INTENT(IN)  :: U(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX,   INTENT(IN)  :: V(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX,   INTENT(OUT) :: Det
   DOUBLE COMPLEX,   INTENT(OUT) :: G(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX,   INTENT(OUT) :: G_Half(1:NVar,1:NVar,1:NSpin)
#endif
   DOUBLE PRECISION              :: SA(1:NOccA,1:NOccA)
   DOUBLE PRECISION              :: SB(1:NOccA,1:NOccB)
   DOUBLE PRECISION              :: InvSA(1:NOccA,1:NOccA)
   DOUBLE PRECISION              :: InvSB(1:NOccB,1:NOccB)
   DOUBLE PRECISION              :: Wrk(1:NVar**2)
   DOUBLE PRECISION              :: DetA, DetB
   DOUBLE PRECISION, PARAMETER   :: Zero= 0D+00
   DOUBLE PRECISION, PARAMETER   :: One = 1D+00
   INTEGER                       :: i

   !
   !--- Initialize ---
   !
   G_Half = Zero 
   G      = Zero
   !
   !Note we assume real wavefunction, so we treat Hetmite conjugate and transpose in the same manner
   !
   !---  S= [Vt.U*] ---
   !
   !     #
   !     #   SA = CWA+.CAT  <=> V+.U
   !     #
   !     ovlp = xp.dot(walker_batch.phia[iw].T, trial.psi0a.conj())
#if 1
   CALL DGeMM('T','N',NOccA,NOccA,NVar,One,V(1,1,1),NVar,U(1,1,1),NVar,Zero,SA,NOccA)
#else
   CALL ZGeMM('T','N',NOccA,NOccA,NVar,One,V(1,1,1),NVar,U(1,1,1),NVar,Zero,SA,NOccA)
#endif
   !
   !---  S-1= [Vt.U*]-1 ---
   !
   !     #
   !     #   SA-1 = [CWA+.CTA]^{-1}   <=> [V+.U]^{-1}
   !     #
   !     ovlp_inv = xp.linalg.inv(ovlp)
   WRITE(IW,*)"SA Matrix"
   CALL PrSq(SA,NOccA,NOccA,NOccA)
   CALL GetInvLR_D(NOccA,SA,InvSA,DetA)
   WRITE(IW,*)"SA-1 Matrix"
   CALL PrSq(InvSA,NOccA,NOccA,NOccA)
   WRITE(IW,*)"V Matrix"
   CALL PrSq(V,NVar,NVar,NVar)
   !
   !--- G_Half =  [(Vt.U*)^{-1}].[Vt] ---
   !
   !     #
   !     #  Green_Half_alpha = SA-1.CWA+  <=? [V+.U]^{-1}.V+
   !     #
   !     walker_batch.Ghalfa[iw] = xp.dot(ovlp_inv, walker_batch.phia[iw].T)
   CALL DGeMM('N','T',NOccA,NVar,NVar,One,InvSA,NOccA,V(1,1,1),NVar,Zero,G_Half(1,1,1),NVar)
   WRITE(IW,*)"G_Half  (Occ x NVar)  (Vir x NVar = 0)"
   CALL PrSq(G_Half,NVar,NVar,NVar)
   !
   !--- G =  [U*.(Vt.U*)^{-1}.Vt]+ ---
   !     #
   !     #  [Green_alpha]_pq = [CTA.SA-1.CWA+]_pq = [CWA.[CTA+.CWA]^{-1}/CTA]_qp
   !     #
   !     walker_batch.Ga[iw] = xp.dot(trial.psi0a.conj(), walker_batch.Ghalfa[iw])
   CALL DGeMM('N','N',NOccA,NOccA,NVar,One,U(1,1,1),NVar,G_Half(1,1,1),NVar,Zero,G(1,1,1),NVar)
   !
   !
   !
   CALL DGeMM('T','N',NOccB,NOccB,NVar,One,V(1,1,2),NVar,U(1,1,2),NVar,Zero,SB,NOccB)
   CALL GetInvLR_D(NOccB,SB,InvSB,DetB)
   CALL DGeMM('N','T',NOccB,NVar,NVar,One,InvSB,NOccB,V(1,1,2),NVar,Zero,G_Half(1,1,2),NVar)
   CALL DGeMM('N','N',NOccB,NOccB,NVar,One,U(1,1,2),NVar,G_Half(1,1,2),NVar,Zero,G(1,1,2),NVar)
   
   WRITE(IW,'("Determinants of the overlap ")')
   WRITE(IW,'("Det A =",F20.12)')DetA
   WRITE(IW,'("Det B =",F20.12)')DetB
   WRITE(IW,'("Half rotated Greens function ")')
   DO i=1,2
      CALL PrSq(G(1,1,i),NVar,NVar,NVar)
   ENDDO
   Det = DetA*DetB

END SUBROUTINE 


SUBROUTINE AFQMC_Read_MRSFCIS_Header(FileName,NDet,NOcc,NSpin)
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN) :: FileName
   INTEGER, INTENT(OUT)         :: NDet
   INTEGER, INTENT(IN)          :: NSpin
   INTEGER, INTENT(IN)          :: NOcc(1:NSpin)
   INTEGER                      :: IU, IOS, NA, NB

   IF (NSpin /= 2) ERROR STOP "MRSF-CIS trials require explicit alpha and beta blocks"
   OPEN(NEWUNIT=IU,FILE=TRIM(FileName),STATUS='OLD',ACTION='READ',IOSTAT=IOS)
   IF (IOS /= 0) ERROR STOP "Unable to open the MRSF-CIS trial file"
   READ(IU,*,IOSTAT=IOS) NDet,NA,NB
   CLOSE(IU)
   IF (IOS /= 0 .OR. NDet < 1) ERROR STOP "Invalid MRSF-CIS trial header"
   IF (NA /= NOcc(1) .OR. NB /= NOcc(2)) THEN
      ERROR STOP "MRSF-CIS electron counts do not match MULT/NCH"
   ENDIF
END SUBROUTINE


SUBROUTINE AFQMC_Read_MRSFCIS_Trial(FileName,NVar,NSpin,NOcc,NDet,Coeff,Phi,DetOcc,IDom)
   USE IO_Files,     ONLY: IW
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN)  :: FileName
   INTEGER,          INTENT(IN)  :: NVar, NSpin, NDet
   INTEGER,          INTENT(IN)  :: NOcc(1:NSpin)
   DOUBLE COMPLEX,   INTENT(OUT) :: Coeff(1:NDet)
   DOUBLE COMPLEX,   INTENT(OUT) :: Phi(1:NVar,1:NVar,1:NSpin,1:NDet)
   INTEGER,          INTENT(OUT) :: DetOcc(1:NVar,1:NSpin,1:NDet)
   INTEGER,          INTENT(OUT) :: IDom
   INTEGER :: IU, IOS, NDFile, NA, NB, D, S, I, J, Orb
   DOUBLE PRECISION :: CR, CI, Norm2, MaxCoeff

   IF (NSpin /= 2) ERROR STOP "MRSF-CIS trials require NSpin=2"
   OPEN(NEWUNIT=IU,FILE=TRIM(FileName),STATUS='OLD',ACTION='READ',IOSTAT=IOS)
   IF (IOS /= 0) ERROR STOP "Unable to open the MRSF-CIS trial file"
   READ(IU,*,IOSTAT=IOS) NDFile,NA,NB
   IF (IOS /= 0 .OR. NDFile /= NDet) ERROR STOP "Inconsistent MRSF-CIS trial header"

   Phi = (0.D0,0.D0)
   DetOcc = 0
   DO D=1,NDet
      READ(IU,*,IOSTAT=IOS) CR,CI, &
     &   (DetOcc(I,1,D),I=1,NOcc(1)),(DetOcc(I,2,D),I=1,NOcc(2))
      IF (IOS /= 0) ERROR STOP "Invalid determinant record in MRSF-CIS trial file"
      Coeff(D) = CMPLX(CR,CI,KIND=KIND(Coeff))
      DO S=1,NSpin
         DO I=1,NOcc(S)
            Orb = DetOcc(I,S,D)
            IF (Orb < 1 .OR. Orb > NVar) ERROR STOP "MRSF-CIS orbital index is out of range"
            DO J=1,I-1
               IF (DetOcc(J,S,D) == Orb) ERROR STOP "Duplicate orbital in trial determinant"
            ENDDO
            ! Column ordering is part of the determinant phase convention.
            Phi(Orb,I,S,D) = (1.D0,0.D0)
         ENDDO
      ENDDO
   ENDDO
   CLOSE(IU)

   Norm2 = SUM(ABS(Coeff)**2)
   IF (Norm2 <= TINY(1.D0)) ERROR STOP "All MRSF-CIS coefficients are zero"
   Coeff = Coeff/DSQRT(Norm2)
   IDom = 1
   MaxCoeff = ABS(Coeff(1))
   DO D=2,NDet
      IF (ABS(Coeff(D)) > MaxCoeff) THEN
         IDom = D
         MaxCoeff = ABS(Coeff(D))
      ENDIF
   ENDDO
   IF (main_rank) THEN
      WRITE(IW,'(1X,"Spin-flip CI trial: ",I0," determinants; dominant determinant ",I0)')NDet,IDom
   ENDIF
END SUBROUTINE


SUBROUTINE AFQMC_Report_Trial_Spin_Z(NVar,NSpin,NOcc,NDet,Coeff,DetOcc,TrialKind)
   ! Evaluate <S**2> from ||S_+ Psi||**2 + M_S(M_S+1).  This expression is
   ! exact for a fixed-M_S CI expansion in a common orthonormal spatial-orbital
   ! basis.  Occupation lists are canonicalized here so arbitrary determinant
   ! column order (and its associated phase) remains valid.
   USE IO_Files,     ONLY: IW
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   INTEGER,          INTENT(IN) :: NVar,NSpin,NDet,TrialKind
   INTEGER,          INTENT(IN) :: NOcc(1:NSpin)
   INTEGER,          INTENT(IN) :: DetOcc(1:NVar,1:NSpin,1:NDet)
   DOUBLE COMPLEX,   INTENT(IN) :: Coeff(1:NDet)
   INTEGER :: MaxTerm,NTerm,D,I,J,K,P,Exponent,NLess,Parity,Pos
   INTEGER, ALLOCATABLE :: OutA(:,:),OutB(:,:)
   INTEGER :: A(1:NVar),B(1:NVar),ATmp(1:NVar),BTmp(1:NVar)
   DOUBLE COMPLEX, ALLOCATABLE :: Amp(:)
   DOUBLE COMPLEX :: CanonCoeff,Term
   DOUBLE PRECISION :: MS,S2,RaiseNorm,SEff
   LOGICAL :: Found,Occupied
   CHARACTER(LEN=8) :: Label

   IF (NSpin /= 2) ERROR STOP "Spin diagnostic requires alpha and beta blocks"
   MaxTerm=MAX(1,NDet*NOcc(2))
   ALLOCATE(OutA(1:NVar,1:MaxTerm),OutB(1:NVar,1:MaxTerm),Amp(1:MaxTerm))
   OutA=0
   OutB=0
   Amp=(0.D0,0.D0)
   NTerm=0

   DO D=1,NDet
      A=0
      B=0
      A(1:NOcc(1))=DetOcc(1:NOcc(1),1,D)
      B(1:NOcc(2))=DetOcc(1:NOcc(2),2,D)
      Parity=1
      DO I=1,NOcc(1)-1
         DO J=I+1,NOcc(1)
            IF (A(I) > A(J)) THEN
               P=A(I); A(I)=A(J); A(J)=P
               Parity=-Parity
            ENDIF
         ENDDO
      ENDDO
      DO I=1,NOcc(2)-1
         DO J=I+1,NOcc(2)
            IF (B(I) > B(J)) THEN
               P=B(I); B(I)=B(J); B(J)=P
               Parity=-Parity
            ENDIF
         ENDDO
      ENDDO
      CanonCoeff=DBLE(Parity)*Coeff(D)

      DO J=1,NOcc(2)
         P=B(J)
         Occupied=.FALSE.
         DO I=1,NOcc(1)
            IF (A(I) == P) Occupied=.TRUE.
         ENDDO
         IF (Occupied) CYCLE

         NLess=COUNT(A(1:NOcc(1)) < P)
         Exponent=NOcc(1)+(J-1)+NLess
         Term=CanonCoeff
         IF (MOD(Exponent,2) /= 0) Term=-Term

         ATmp=0
         BTmp=0
         Pos=1
         DO I=1,NOcc(1)
            IF (A(I) < P) THEN
               ATmp(Pos)=A(I); Pos=Pos+1
            ENDIF
         ENDDO
         ATmp(Pos)=P; Pos=Pos+1
         DO I=1,NOcc(1)
            IF (A(I) > P) THEN
               ATmp(Pos)=A(I); Pos=Pos+1
            ENDIF
         ENDDO
         Pos=1
         DO I=1,NOcc(2)
            IF (I /= J) THEN
               BTmp(Pos)=B(I); Pos=Pos+1
            ENDIF
         ENDDO

         Found=.FALSE.
         DO K=1,NTerm
            IF (ALL(OutA(1:NOcc(1)+1,K) == ATmp(1:NOcc(1)+1)) .AND. &
           &    ALL(OutB(1:MAX(0,NOcc(2)-1),K) == BTmp(1:MAX(0,NOcc(2)-1)))) THEN
               Amp(K)=Amp(K)+Term
               Found=.TRUE.
               EXIT
            ENDIF
         ENDDO
         IF (.NOT.Found) THEN
            NTerm=NTerm+1
            IF (NTerm > MaxTerm) ERROR STOP "Internal spin-diagnostic workspace exceeded"
            OutA(:,NTerm)=ATmp
            OutB(:,NTerm)=BTmp
            Amp(NTerm)=Term
         ENDIF
      ENDDO
   ENDDO

   RaiseNorm=0.D0
   IF (NTerm > 0) RaiseNorm=SUM(ABS(Amp(1:NTerm))**2)
   MS=0.5D0*DBLE(NOcc(1)-NOcc(2))
   S2=MS*(MS+1.D0)+RaiseNorm
   SEff=0.5D0*(SQRT(MAX(0.D0,1.D0+4.D0*S2))-1.D0)
   ! OpenQP port: the branch was inverted -- kind 1 is SF-CIS, kind 2 is
   ! MRSF-CIS (as the trial-file selection elsewhere already assumes), so an
   ! MRSF trial was reported as SF-CIS and vice versa.
   IF (TrialKind == 1) THEN
      Label="SF-CIS  "
   ELSE
      Label="MRSF-CIS"
   ENDIF
   IF (main_rank) THEN
      WRITE(IW,'(1X,A," trial <S^2> = ",F16.10,"; effective S = ",F12.8)')&
     &     TRIM(Label),S2,SEff
   ENDIF
   DEALLOCATE(OutA,OutB,Amp)
END SUBROUTINE


SUBROUTINE AFQMC_Read_OO_Orbitals_Z(FileName,NVar,NSpin,C)
   USE IO_Files,     ONLY: IW
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN)  :: FileName
   INTEGER,          INTENT(IN)  :: NVar, NSpin
   DOUBLE COMPLEX,   INTENT(OUT) :: C(1:NVar,1:NVar,1:NSpin)
   INTEGER :: IU, IOS, NVFile, NSFile, S, P, Q, I, J
   DOUBLE PRECISION :: CR, CI, Err, MaxErr
   DOUBLE COMPLEX :: Dot

   OPEN(NEWUNIT=IU,FILE=TRIM(FileName),STATUS='OLD',ACTION='READ',IOSTAT=IOS)
   IF (IOS /= 0) ERROR STOP "Unable to open the OO orbital file"
   READ(IU,*,IOSTAT=IOS) NVFile,NSFile
   IF (IOS /= 0 .OR. NVFile /= NVar) ERROR STOP "Invalid OO orbital header"
   IF (NSFile /= 1 .AND. NSFile /= NSpin) THEN
      ERROR STOP "OO orbital file must contain one spatial block or NSpin blocks"
   ENDIF
   C = (0.D0,0.D0)
   DO S=1,NSFile
      DO Q=1,NVar
         DO P=1,NVar
            READ(IU,*,IOSTAT=IOS) CR,CI
            IF (IOS /= 0) ERROR STOP "Invalid matrix element in OO orbital file"
            C(P,Q,S) = CMPLX(CR,CI,KIND=KIND(C))
         ENDDO
      ENDDO
   ENDDO
   CLOSE(IU)
   IF (NSFile == 1) THEN
      DO S=2,NSpin
         C(:,:,S) = C(:,:,1)
      ENDDO
   ENDIF

   MaxErr = 0.D0
   DO S=1,NSpin
      DO I=1,NVar
         DO J=1,NVar
            Dot = SUM(CONJG(C(:,I,S))*C(:,J,S))
            IF (I == J) Dot = Dot-(1.D0,0.D0)
            Err = ABS(Dot)
            MaxErr = MAX(MaxErr,Err)
         ENDDO
      ENDDO
   ENDDO
   IF (MaxErr > 1.D-7) ERROR STOP "OO orbitals are not unitary in the canonical MO basis"
   IF (main_rank) WRITE(IW,'(1X,"OO trial orbitals loaded; unitarity error = ",ES12.4)')MaxErr
END SUBROUTINE


SUBROUTINE AFQMC_Apply_Orbital_Transform_Z(NVar,NSpin,NDet,C,Phi)
   IMPLICIT NONE
   INTEGER,        INTENT(IN)    :: NVar, NSpin, NDet
   DOUBLE COMPLEX, INTENT(IN)    :: C(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX, INTENT(INOUT) :: Phi(1:NVar,1:NVar,1:NSpin,1:NDet)
   DOUBLE COMPLEX :: Tmp(1:NVar,1:NVar)
   DOUBLE COMPLEX, PARAMETER :: One=(1.D0,0.D0), Zero=(0.D0,0.D0)
   INTEGER :: D, S

   DO D=1,NDet
      DO S=1,NSpin
         CALL ZGEMM('N','N',NVar,NVar,NVar,One,C(1,1,S),NVar,&
        &           Phi(1,1,S,D),NVar,Zero,Tmp,NVar)
         Phi(:,:,S,D) = Tmp
      ENDDO
   ENDDO
END SUBROUTINE


SUBROUTINE AFQMC_Read_Lower_States_Z(FileName,NDet,NLow,LowerCoeff)
   USE IO_Files,     ONLY: IW
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN)  :: FileName
   INTEGER,          INTENT(IN)  :: NDet, NLow
   DOUBLE COMPLEX,   INTENT(OUT) :: LowerCoeff(1:NDet,1:NLow)
   DOUBLE PRECISION :: Record(1:2*NDet), Norm2
   DOUBLE COMPLEX :: Projection
   INTEGER :: IU, IOS, NDFile, NLFile, K, J, D

   OPEN(NEWUNIT=IU,FILE=TRIM(FileName),STATUS='OLD',ACTION='READ',IOSTAT=IOS)
   IF (IOS /= 0) ERROR STOP "Unable to open lower-state coefficient file"
   READ(IU,*,IOSTAT=IOS) NLFile,NDFile
   IF (IOS /= 0 .OR. NLFile /= NLow .OR. NDFile /= NDet) THEN
      ERROR STOP "Lower-state header does not match NLOW/NDet"
   ENDIF
   DO K=1,NLow
      READ(IU,*,IOSTAT=IOS) (Record(J),J=1,2*NDet)
      IF (IOS /= 0) ERROR STOP "Invalid lower-state coefficient record"
      DO D=1,NDet
         LowerCoeff(D,K) = CMPLX(Record(2*D-1),Record(2*D),KIND=KIND(LowerCoeff))
      ENDDO
   ENDDO
   CLOSE(IU)

   ! Modified Gram-Schmidt gives an orthonormal basis for the supplied lower
   ! state span without assuming the upstream roots are numerically exact.
   DO K=1,NLow
      DO J=1,K-1
         Projection = DOT_PRODUCT(LowerCoeff(:,J),LowerCoeff(:,K))
         LowerCoeff(:,K) = LowerCoeff(:,K)-Projection*LowerCoeff(:,J)
      ENDDO
      Norm2 = SUM(ABS(LowerCoeff(:,K))**2)
      IF (Norm2 <= 1.D-14) ERROR STOP "Linearly dependent lower-state trials"
      LowerCoeff(:,K) = LowerCoeff(:,K)/DSQRT(Norm2)
   ENDDO
   IF (main_rank) WRITE(IW,'(1X,I0," lower-state trial vectors loaded")')NLow
END SUBROUTINE


SUBROUTINE AFQMC_Project_Lower_States_Z(NDet,NLow,LowerCoeff,Coeff)
   USE IO_Files,     ONLY: IW
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   INTEGER,        INTENT(IN)    :: NDet, NLow
   DOUBLE COMPLEX, INTENT(IN)    :: LowerCoeff(1:NDet,1:NLow)
   DOUBLE COMPLEX, INTENT(INOUT) :: Coeff(1:NDet)
   DOUBLE COMPLEX :: Projection
   DOUBLE PRECISION :: Norm2, Removed
   INTEGER :: K

   Removed = 0.D0
   DO K=1,NLow
      Projection = DOT_PRODUCT(LowerCoeff(:,K),Coeff)
      Removed = Removed+ABS(Projection)**2
      Coeff = Coeff-Projection*LowerCoeff(:,K)
   ENDDO
   Norm2 = SUM(ABS(Coeff)**2)
   IF (Norm2 <= 1.D-12) ERROR STOP "Target trial lies in the supplied lower-state span"
   Coeff = Coeff/DSQRT(Norm2)
   IF (main_rank) WRITE(IW,'(1X,"Lower-state component removed from target = ",ES12.4)')Removed
END SUBROUTINE


SUBROUTINE AFQMC_Dominant_Coefficient_Z(NDet,Coeff,IDom)
   IMPLICIT NONE
   INTEGER,        INTENT(IN)  :: NDet
   DOUBLE COMPLEX, INTENT(IN)  :: Coeff(1:NDet)
   INTEGER,        INTENT(OUT) :: IDom
   INTEGER :: D

   IDom = 1
   DO D=2,NDet
      IF (ABS(Coeff(D)) > ABS(Coeff(IDom))) IDom = D
   ENDDO
END SUBROUTINE


SUBROUTINE AFQMC_MultiDet_Overlap_Z(NVar,NSpin,NOcc,NDet,Coeff,TrialPhi,WalkerPhi,Ovlp)
   IMPLICIT NONE
   INTEGER,        INTENT(IN)  :: NVar, NSpin, NDet
   INTEGER,        INTENT(IN)  :: NOcc(1:NSpin)
   DOUBLE COMPLEX, INTENT(IN)  :: Coeff(1:NDet)
   DOUBLE COMPLEX, INTENT(IN)  :: TrialPhi(1:NVar,1:NVar,1:NSpin,1:NDet)
   DOUBLE COMPLEX, INTENT(IN)  :: WalkerPhi(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX, INTENT(OUT) :: Ovlp
   DOUBLE COMPLEX :: OvlpD
   INTEGER :: D

   Ovlp = (0.D0,0.D0)
   DO D=1,NDet
      CALL AFQMC_Overlap_Z(NVar,NSpin,NOcc,TrialPhi(1,1,1,D),WalkerPhi,OvlpD)
      Ovlp = Ovlp + CONJG(Coeff(D))*OvlpD
   ENDDO
END SUBROUTINE


SUBROUTINE AFQMC_Lower_State_Diagnostic_Z(NVar,NSpin,NOcc,NWlk,Walker,&
&                  NDet,NLow,LowerCoeff,TrialPhi,IStep,LowMax,LowCap,Enforce)
   USE AFQMC_Module, ONLY: AFQMC_Walker_Type
   USE IO_Files,     ONLY: IW
! OpenQP port: the upstream #ifdef AFQMC_STANDALONE selected between Gellan's
! real MPI layer and a serial stub set. liboqp always uses the serial stubs in
! afqmc_host.F90 (it parallelises with OpenMP), so the guard is resolved here.
   USE MPI_Parallel, ONLY: my_rank,NProc,MPI_Comm_World,MPI_SUM,main_rank,Gel_AllReduce
   IMPLICIT NONE
   INTEGER,                  INTENT(IN) :: NVar,NSpin,NWlk,NDet,NLow,IStep
   INTEGER,                  INTENT(IN) :: NOcc(1:NSpin)
   TYPE(AFQMC_Walker_Type),  INTENT(IN) :: Walker(1:NWlk)
   DOUBLE COMPLEX,           INTENT(IN) :: LowerCoeff(1:NDet,1:NLow)
   DOUBLE COMPLEX,           INTENT(IN) :: TrialPhi(1:NVar,1:NVar,1:NSpin,1:NDet)
   DOUBLE PRECISION,         INTENT(IN) :: LowMax
   DOUBLE PRECISION,         INTENT(IN) :: LowCap
   LOGICAL,                  INTENT(IN) :: Enforce
   DOUBLE PRECISION :: LeakageRMS(1:NLow), DummyRMS(1:NLow)
   DOUBLE COMPLEX :: LeakageAmp(1:NLow), DummyAmp(1:NLow), Ratio,RobustRatio
   DOUBLE PRECISION :: TotalWeight(1), DummyWeight(1)
   DOUBLE PRECISION :: W
   DOUBLE COMPLEX :: LowerOvlp
   INTEGER :: IWlk, K

   LeakageAmp = (0.D0,0.D0)
   LeakageRMS = 0.D0
   TotalWeight(1) = 0.D0
   DO IWlk=my_rank+1,NWlk,NProc
      W = ABS(Walker(IWlk)%Weight)
      IF (W <= 0.D0) CYCLE
      IF (ABS(Walker(IWlk)%Ovlp) <= 1.D-14) THEN
         ERROR STOP "Finite-weight walker has zero target-trial overlap"
      ENDIF
      TotalWeight(1) = TotalWeight(1)+W
      DO K=1,NLow
         CALL AFQMC_MultiDet_Overlap_Z(NVar,NSpin,NOcc,NDet,&
        &        LowerCoeff(1,K),TrialPhi,Walker(IWlk)%UW,LowerOvlp)
         Ratio = LowerOvlp/Walker(IWlk)%Ovlp
         ! Importance sampling represents the projected state with
         ! w_i |phi_i>/<Psi_T|phi_i>.  The coherent weighted mean therefore
         ! estimates its lower-state amplitude; the RMS is also useful as a
         ! variance/conditioning diagnostic but must not trigger collapse by
         ! itself because legitimate cancellation occurs between walkers.
         RobustRatio=Ratio
         IF (LowCap > 0.D0 .AND. ABS(RobustRatio) > LowCap) THEN
            RobustRatio=LowCap*RobustRatio/ABS(RobustRatio)
         ENDIF
         LeakageAmp(K) = LeakageAmp(K)+W*RobustRatio
         LeakageRMS(K) = LeakageRMS(K)+W*ABS(Ratio)**2
      ENDDO
   ENDDO
   CALL Gel_AllReduce(LeakageAmp,DummyAmp,NLow,'Z',MPI_SUM,MPI_Comm_World,1)
   CALL Gel_AllReduce(LeakageRMS,DummyRMS,NLow,'D',MPI_SUM,MPI_Comm_World,1)
   CALL Gel_AllReduce(TotalWeight,DummyWeight,1,'D',MPI_SUM,MPI_Comm_World,1)
   IF (TotalWeight(1) <= 0.D0) ERROR STOP "No finite-weight walkers in leakage diagnostic"
   LeakageAmp = LeakageAmp/TotalWeight(1)
   LeakageRMS = SQRT(LeakageRMS/TotalWeight(1))
   IF (main_rank) THEN
      WRITE(IW,'(1X,"Lower-state coherent leakage at step ",I0,":",*(1X,ES12.4))')&
     &      IStep,ABS(LeakageAmp)
      WRITE(IW,'(1X,"Lower-state walker RMS ratios at step ",I0,":",*(1X,ES12.4))')&
     &      IStep,LeakageRMS
   ENDIF
   IF (Enforce .AND. LowMax > 0.D0 .AND. MAXVAL(ABS(LeakageAmp)) > LowMax) THEN
      ERROR STOP "Lower-state leakage exceeded LOWMAX; excited-state projection stopped"
   ENDIF
END SUBROUTINE


SUBROUTINE AFQMC_TrialEnsemble_Green_Z(NVar,NSpin,NOcc,NDet,Coeff,TrialPhi,G)
   ! Stable diagonal-ensemble density used only for the AFQMC mean-field
   ! subtraction.  The coherent MRSF-CIS trial is retained in overlap, force
   ! bias, and energy calculations below.
   IMPLICIT NONE
   INTEGER,        INTENT(IN)  :: NVar, NSpin, NDet
   INTEGER,        INTENT(IN)  :: NOcc(1:NSpin)
   DOUBLE COMPLEX, INTENT(IN)  :: Coeff(1:NDet)
   DOUBLE COMPLEX, INTENT(IN)  :: TrialPhi(1:NVar,1:NVar,1:NSpin,1:NDet)
   DOUBLE COMPLEX, INTENT(OUT) :: G(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX :: GD(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX :: GHalf(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX :: OvlpD
   INTEGER :: D

   G = (0.D0,0.D0)
   DO D=1,NDet
      CALL AFQMC_Greens_Functions_Z(NVar,NSpin,NOcc,TrialPhi(1,1,1,D),&
     &                              TrialPhi(1,1,1,D),OvlpD,GD,GHalf)
      G = G + ABS(Coeff(D))**2*GD
   ENDDO
END SUBROUTINE


SUBROUTINE AFQMC_MultiDet_Greens_Z(NVar,NSpin,NOcc,NDet,Coeff,&
&                                   TrialPhi,WalkerPhi,Ovlp,G)
   ! Mixed Green's function for a coherent multideterminant trial:
   ! G = sum_D c_D^* <D|W> G_DW / sum_D c_D^* <D|W>.
   IMPLICIT NONE
   INTEGER,        INTENT(IN)  :: NVar, NSpin, NDet
   INTEGER,        INTENT(IN)  :: NOcc(1:NSpin)
   DOUBLE COMPLEX, INTENT(IN)  :: Coeff(1:NDet)
   DOUBLE COMPLEX, INTENT(IN)  :: TrialPhi(1:NVar,1:NVar,1:NSpin,1:NDet)
   DOUBLE COMPLEX, INTENT(IN)  :: WalkerPhi(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX, INTENT(OUT) :: Ovlp
   DOUBLE COMPLEX, INTENT(OUT) :: G(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX :: GNumerator(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX :: OvlpD
   INTEGER :: D

   Ovlp = (0.D0,0.D0)
   G = (0.D0,0.D0)
   DO D=1,NDet
      CALL AFQMC_Overlap_Green_Numerator_Z(NVar,NSpin,NOcc,&
     &             TrialPhi(1,1,1,D),WalkerPhi,OvlpD,GNumerator)
      Ovlp = Ovlp + CONJG(Coeff(D))*OvlpD
      G = G + CONJG(Coeff(D))*GNumerator
   ENDDO
   IF (ABS(Ovlp) <= 100.D0*TINY(1.D0)) ERROR STOP "Zero total MRSF-CIS trial overlap"
   G = G/Ovlp
END SUBROUTINE


SUBROUTINE AFQMC_MultiDet_Local_Energy_Z(NVar,NSpin,NOcc,NChol,NDet,&
&                                         Coeff,TrialPhi,WalkerPhi,ENu,Hc,Chol,E)
   IMPLICIT NONE
   INTEGER,          INTENT(IN)  :: NVar, NSpin, NChol, NDet
   INTEGER,          INTENT(IN)  :: NOcc(1:NSpin)
   DOUBLE COMPLEX,   INTENT(IN)  :: Coeff(1:NDet)
   DOUBLE COMPLEX,   INTENT(IN)  :: TrialPhi(1:NVar,1:NVar,1:NSpin,1:NDet)
   DOUBLE COMPLEX,   INTENT(IN)  :: WalkerPhi(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(IN)  :: ENu
   DOUBLE PRECISION, INTENT(IN)  :: Hc(1:NVar*(NVar+1)/2)
   DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
   DOUBLE COMPLEX,   INTENT(OUT) :: E
   DOUBLE COMPLEX :: Ovlp, OvlpD, WeightD, Numerator, ED
   DOUBLE COMPLEX :: SingularNumerator
   DOUBLE COMPLEX :: G(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX :: GHalf(1:NVar,1:NVar,1:NSpin)
   INTEGER :: D

   Ovlp = (0.D0,0.D0)
   Numerator = (0.D0,0.D0)
   DO D=1,NDet
      CALL AFQMC_Overlap_Z(NVar,NSpin,NOcc,TrialPhi(1,1,1,D),WalkerPhi,OvlpD)
      WeightD = CONJG(Coeff(D))*OvlpD
      Ovlp = Ovlp + WeightD
      IF (ABS(OvlpD) > 100.D0*TINY(1.D0)) THEN
         CALL AFQMC_Greens_Functions_Z(NVar,NSpin,NOcc,TrialPhi(1,1,1,D),&
        &                              WalkerPhi,OvlpD,G,GHalf)
         CALL AFQMC_Energy_Estimator_Z(NVar,NSpin,NChol,G,GHalf,ENu,Hc,Chol,ED,.TRUE.)
         Numerator = Numerator + WeightD*ED
      ELSEIF (ABS(Coeff(D)) > 100.D0*TINY(1.D0)) THEN
         ! OvlpD*ED has a finite Slater--Condon limit even when the overlap is
         ! exactly zero.  Evaluate that Hamiltonian numerator directly rather
         ! than silently dropping a singly/doubly connected determinant.
         CALL AFQMC_Singular_Local_Energy_Numerator_Z(NVar,NSpin,NOcc,NChol,&
        &       TrialPhi(1,1,1,D),WalkerPhi,ENu,Hc,Chol,SingularNumerator)
         Numerator = Numerator + CONJG(Coeff(D))*SingularNumerator
      ENDIF
   ENDDO
   IF (ABS(Ovlp) <= 100.D0*TINY(1.D0)) ERROR STOP "Zero total MRSF-CIS trial overlap"
   E = Numerator/Ovlp
END SUBROUTINE


SUBROUTINE AFQMC_Singular_Local_Energy_Numerator_Z(NVar,NSpin,NOcc,NChol,&
&                                                   U,V,ENu,Hc,Chol,HNumerator)
   ! Exact <U|H|V> for the singular-overlap fallback.  For each Cholesky
   ! one-body operator L, derivatives of
   !   det[U^+ exp(t L) V]
   ! give <U|rho_L|V> and <U|rho_L^2|V> without forming S^{-1}.  Column-
   ! replacement determinant derivatives remain valid at any rank of S.
   IMPLICIT NONE
   INTEGER,          INTENT(IN)  :: NVar,NSpin,NChol
   INTEGER,          INTENT(IN)  :: NOcc(1:NSpin)
   DOUBLE COMPLEX,   INTENT(IN)  :: U(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX,   INTENT(IN)  :: V(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(IN)  :: ENu
   DOUBLE PRECISION, INTENT(IN)  :: Hc(1:NVar*(NVar+1)/2)
   DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
   DOUBLE COMPLEX,   INTENT(OUT) :: HNumerator
   DOUBLE COMPLEX :: GNumerator(1:NVar,1:NVar,1:NSpin),Ovlp
   DOUBLE COMPLEX :: S(1:NVar,1:NVar),S1(1:NVar,1:NVar),S2(1:NVar,1:NVar)
   DOUBLE COMPLEX :: LV(1:NVar,1:NVar),L2V(1:NVar,1:NVar)
   DOUBLE COMPLEX :: F0(1:NSpin),F1(1:NSpin),F2(1:NSpin),Q1(1:NSpin)
   DOUBLE COMPLEX :: F2Total,QTotal,ProductOther,EOne,ETwo
   INTEGER :: I,J,K,P,Q,M,ISpin,JSpin,NO,IJ

   CALL AFQMC_Overlap_Green_Numerator_Z(NVar,NSpin,NOcc,U,V,Ovlp,GNumerator)
   EOne = CMPLX(ENu,0.D0,KIND=KIND(EOne))*Ovlp
   IJ = 0
   DO I=1,NVar
      DO J=1,I-1
         IJ = IJ+1
         DO ISpin=1,NSpin
            EOne = EOne+Hc(IJ)*(GNumerator(J,I,ISpin)+GNumerator(I,J,ISpin))
         ENDDO
      ENDDO
      IJ = IJ+1
      DO ISpin=1,NSpin
         EOne = EOne+Hc(IJ)*GNumerator(I,I,ISpin)
      ENDDO
   ENDDO

   ETwo = (0.D0,0.D0)
   DO M=1,NChol
      F0 = (1.D0,0.D0)
      F1 = (0.D0,0.D0)
      F2 = (0.D0,0.D0)
      Q1 = (0.D0,0.D0)
      DO ISpin=1,NSpin
         NO = NOcc(ISpin)
         IF (NO == 0) CYCLE
         LV = (0.D0,0.D0)
         L2V = (0.D0,0.D0)
         DO J=1,NO
            DO P=1,NVar
               DO Q=1,NVar
                  LV(P,J) = LV(P,J)+Chol(P,Q,M)*V(Q,J,ISpin)
               ENDDO
            ENDDO
            DO P=1,NVar
               DO Q=1,NVar
                  L2V(P,J) = L2V(P,J)+Chol(P,Q,M)*LV(Q,J)
               ENDDO
            ENDDO
         ENDDO
         S = (0.D0,0.D0)
         S1 = (0.D0,0.D0)
         S2 = (0.D0,0.D0)
         DO I=1,NO
            DO J=1,NO
               DO P=1,NVar
                  S(I,J) = S(I,J)+CONJG(U(P,I,ISpin))*V(P,J,ISpin)
                  S1(I,J) = S1(I,J)+CONJG(U(P,I,ISpin))*LV(P,J)
                  S2(I,J) = S2(I,J)+CONJG(U(P,I,ISpin))*L2V(P,J)
               ENDDO
            ENDDO
         ENDDO
         CALL AFQMC_Determinant_Derivatives_Z(NO,S,NVar,S1,NVar,S2,NVar,&
        &                                      F0(ISpin),F1(ISpin),F2(ISpin),Q1(ISpin))
      ENDDO

      F2Total = (0.D0,0.D0)
      QTotal = (0.D0,0.D0)
      DO ISpin=1,NSpin
         ProductOther = (1.D0,0.D0)
         DO JSpin=1,NSpin
            IF (JSpin /= ISpin) ProductOther = ProductOther*F0(JSpin)
         ENDDO
         F2Total = F2Total+F2(ISpin)*ProductOther
         QTotal = QTotal+Q1(ISpin)*ProductOther
      ENDDO
      DO ISpin=1,NSpin
         DO JSpin=ISpin+1,NSpin
            ProductOther = (1.D0,0.D0)
            DO K=1,NSpin
               IF (K /= ISpin .AND. K /= JSpin) ProductOther=ProductOther*F0(K)
            ENDDO
            F2Total = F2Total+2.D0*F1(ISpin)*F1(JSpin)*ProductOther
         ENDDO
      ENDDO
      ETwo = ETwo+0.5D0*(F2Total-QTotal)
   ENDDO
   HNumerator = EOne+ETwo
END SUBROUTINE


SUBROUTINE AFQMC_Determinant_Derivatives_Z(N,S,LdS,S1,LdS1,S2,LdS2,&
&                                           F0,F1,F2,Q1)
   ! For A(t)=S+t*S1+t^2*S2/2, return det(A), its first and second
   ! derivatives at zero, and the part of the second derivative containing
   ! S2.  Multilinearity of the determinant makes this exact without S^{-1}.
   IMPLICIT NONE
   INTEGER,        INTENT(IN)  :: N,LdS,LdS1,LdS2
   DOUBLE COMPLEX, INTENT(IN)  :: S(1:LdS,1:MAX(1,N))
   DOUBLE COMPLEX, INTENT(IN)  :: S1(1:LdS1,1:MAX(1,N))
   DOUBLE COMPLEX, INTENT(IN)  :: S2(1:LdS2,1:MAX(1,N))
   DOUBLE COMPLEX, INTENT(OUT) :: F0,F1,F2,Q1
   DOUBLE COMPLEX :: Work(1:MAX(1,N),1:MAX(1,N)),DetTerm
   INTEGER :: I,J

   IF (N == 0) THEN
      F0=(1.D0,0.D0); F1=(0.D0,0.D0); F2=(0.D0,0.D0); Q1=(0.D0,0.D0)
      RETURN
   ENDIF
   CALL AFQMC_Determinant_Z(N,S,LdS,F0)
   F1=(0.D0,0.D0)
   Q1=(0.D0,0.D0)
   DO I=1,N
      Work(1:N,1:N)=S(1:N,1:N)
      Work(1:N,I)=S1(1:N,I)
      CALL AFQMC_Determinant_Z(N,Work,N,DetTerm)
      F1=F1+DetTerm
      Work(1:N,1:N)=S(1:N,1:N)
      Work(1:N,I)=S2(1:N,I)
      CALL AFQMC_Determinant_Z(N,Work,N,DetTerm)
      Q1=Q1+DetTerm
   ENDDO
   F2=Q1
   DO I=1,N
      DO J=I+1,N
         Work(1:N,1:N)=S(1:N,1:N)
         Work(1:N,I)=S1(1:N,I)
         Work(1:N,J)=S1(1:N,J)
         CALL AFQMC_Determinant_Z(N,Work,N,DetTerm)
         F2=F2+2.D0*DetTerm
      ENDDO
   ENDDO
END SUBROUTINE


SUBROUTINE AFQMC_Determinant_Z(N,A,LdA,Det)
   IMPLICIT NONE
   INTEGER,        INTENT(IN)  :: N, LdA
   DOUBLE COMPLEX, INTENT(IN)  :: A(1:LdA,1:MAX(1,N))
   DOUBLE COMPLEX, INTENT(OUT) :: Det
   DOUBLE COMPLEX :: Work(1:MAX(1,N),1:MAX(1,N))
   INTEGER :: IPiv(1:MAX(1,N)), Info, I

   Det = (1.D0,0.D0)
   IF (N == 0) RETURN
   Work(1:N,1:N) = A(1:N,1:N)
   CALL ZGETRF(N,N,Work,N,IPiv,Info)
   IF (Info < 0) ERROR STOP "Invalid argument in AFQMC_Determinant_Z"
   IF (Info > 0) THEN
      Det = (0.D0,0.D0)
      RETURN
   ENDIF
   DO I=1,N
      Det = Det*Work(I,I)
      IF (IPiv(I) /= I) Det = -Det
   ENDDO
END SUBROUTINE


SUBROUTINE AFQMC_Adjugate_Z(N,A,LdA,Adj,LdAdj)
   IMPLICIT NONE
   INTEGER,        INTENT(IN)  :: N, LdA, LdAdj
   DOUBLE COMPLEX, INTENT(IN)  :: A(1:LdA,1:MAX(1,N))
   DOUBLE COMPLEX, INTENT(OUT) :: Adj(1:LdAdj,1:MAX(1,N))
   DOUBLE COMPLEX :: Minor(1:MAX(1,N-1),1:MAX(1,N-1)), DetMinor
   INTEGER :: I, J, R, C, RR, CC

   Adj = (0.D0,0.D0)
   IF (N == 0) RETURN
   IF (N == 1) THEN
      Adj(1,1) = (1.D0,0.D0)
      RETURN
   ENDIF
   DO I=1,N
      DO J=1,N
         RR = 0
         DO R=1,N
            IF (R == J) CYCLE
            RR = RR + 1
            CC = 0
            DO C=1,N
               IF (C == I) CYCLE
               CC = CC + 1
               Minor(RR,CC) = A(R,C)
            ENDDO
         ENDDO
         CALL AFQMC_Determinant_Z(N-1,Minor,MAX(1,N-1),DetMinor)
         IF (MOD(I+J,2) == 0) THEN
            Adj(I,J) = DetMinor
         ELSE
            Adj(I,J) = -DetMinor
         ENDIF
      ENDDO
   ENDDO
END SUBROUTINE


SUBROUTINE AFQMC_Overlap_Green_Numerator_Z(NVar,NSpin,NOcc,U,V,Ovlp,GNumerator)
   ! Return <U|V> and the unnormalized one-body transition density
   ! <U|a_p^+ a_q|V>.  The adjugate fallback remains well-defined when an
   ! individual determinant overlap is exactly singular, which occurs at the
   ! start of an MRSF calculation with orthogonal canonical determinants.
   IMPLICIT NONE
   INTEGER,        INTENT(IN)  :: NVar, NSpin
   INTEGER,        INTENT(IN)  :: NOcc(1:NSpin)
   DOUBLE COMPLEX, INTENT(IN)  :: U(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX, INTENT(IN)  :: V(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX, INTENT(OUT) :: Ovlp
   DOUBLE COMPLEX, INTENT(OUT) :: GNumerator(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX :: S(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX :: Adj(1:NVar,1:NVar), UConj(1:NVar,1:NVar)
   DOUBLE COMPLEX :: Tmp(1:NVar,1:NVar), G(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX :: GHalf(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX :: DetSpin(1:NSpin), DetOther, DetCheck
   DOUBLE COMPLEX, PARAMETER :: One=(1.D0,0.D0), Zero=(0.D0,0.D0)
   INTEGER :: ISpin, JSpin, NO

   S = Zero
   DetSpin = One
   DO ISpin=1,NSpin
      NO = NOcc(ISpin)
      IF (NO == 0) CYCLE
      UConj(1:NVar,1:NO) = CONJG(U(1:NVar,1:NO,ISpin))
      CALL ZGEMM('T','N',NO,NO,NVar,One,V(1,1,ISpin),NVar,&
     &           UConj,NVar,Zero,S(1,1,ISpin),NVar)
      CALL AFQMC_Determinant_Z(NO,S(1,1,ISpin),NVar,DetSpin(ISpin))
   ENDDO
   Ovlp = PRODUCT(DetSpin)

   IF (ABS(Ovlp) > 100.D0*TINY(1.D0)) THEN
      CALL AFQMC_Greens_Functions_Z(NVar,NSpin,NOcc,U,V,DetCheck,G,GHalf)
      GNumerator = Ovlp*G
      RETURN
   ENDIF

   GNumerator = Zero
   DO ISpin=1,NSpin
      NO = NOcc(ISpin)
      IF (NO == 0) CYCLE
      DetOther = One
      DO JSpin=1,NSpin
         IF (JSpin /= ISpin) DetOther = DetOther*DetSpin(JSpin)
      ENDDO
      IF (ABS(DetOther) <= 100.D0*TINY(1.D0)) CYCLE
      CALL AFQMC_Adjugate_Z(NO,S(1,1,ISpin),NVar,Adj,NVar)
      UConj(1:NVar,1:NO) = CONJG(U(1:NVar,1:NO,ISpin))
      CALL ZGEMM('N','N',NVar,NO,NO,One,UConj,NVar,Adj,NVar,Zero,Tmp,NVar)
      CALL ZGEMM('N','T',NVar,NVar,NO,DetOther,Tmp,NVar,V(1,1,ISpin),&
     &           NVar,Zero,GNumerator(1,1,ISpin),NVar)
   ENDDO
END SUBROUTINE


SUBROUTINE AFQMC_Greens_Functions_Z(NVar,NSpin,NOcc,U,V,Det,G,G_Half)
   !
   !  G_pq = <Psi_T|ap+.aq|Psi_W>/<Psi_T|Psi_W> = [CW.(CT+.CW)-1.CT+]_qp
   !
   !           |Psi_T> = U|Psi>  (U=CT)
   !           |Psi_W> = V|Psi>  (V=CW)
   !     であるとすれば
   !        <Psi_T|ap+.aq|Psi_W>/<Psi_T|Psi_W> 
   !        = <Psi|U+.ap+.aq.V|Psi> / <Psi|U+.V|Psi> 
   !        = <Psi|[U+.a_p+.U].U+.V.[V+.a_q.V]|Psi> / <Psi|U+.V|Psi> 
   !        = <Psi|U*_pr.a_r+.[U+.V].a_s.V_qs|Psi> / <Psi|U+.V|Psi>
   !        = U*_pr V_qs <Psi|a_r+.[U+.V].a_s|Psi> / <Psi|U+.V|Psi>
   !        = U*_pr V_qs [U+.V]^{-1}_sr / 
   !        = [V.(U+.V)^{-1}.U+]_qp
   !        = [U.(V+.U)^{-1}.V+]*_pq
   !
   !  p and q are occupied orbitals
   ! 
   !  ここで CTは実数、CWは複素数を想定しないといけない
   !
   USE IO_Files,     ONLY: IW
   USE AFQMC_Module, ONLY: PrintT
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   INTEGER,          INTENT(IN)  :: NVar
   INTEGER,          INTENT(IN)  :: NSpin
   INTEGER,          INTENT(IN)  :: NOcc(1:NSpin)
   DOUBLE COMPLEX,   INTENT(IN)  :: U(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX,   INTENT(IN)  :: V(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX,   INTENT(OUT) :: Det
   DOUBLE COMPLEX,   INTENT(OUT) :: G(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX,   INTENT(OUT) :: G_Half(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX,   TARGET      :: Wrk1(1:NVar**2)
   DOUBLE COMPLEX,   TARGET      :: Wrk2(1:NVar**2)
   DOUBLE COMPLEX,   TARGET      :: Wrk3(1:NVar**2)
   DOUBLE COMPLEX                :: DetX(1:NSpin)
   DOUBLE COMPLEX, PARAMETER     :: Zero = CMPLX(0D+00, 0D+00)
   DOUBLE COMPLEX, PARAMETER     :: One = CMPLX(1D+00, 0D+00)
   INTEGER                       :: I
   INTEGER                       :: NO
   DOUBLE COMPLEX, POINTER       :: S(:,:)
   DOUBLE COMPLEX, POINTER       :: SInv(:,:)
   DOUBLE COMPLEX, POINTER       :: U_Cnj(:,:)
   !
   !
   !
   IF (main_rank .AND. PrintT%Walker > 1) THEN 
      WRITE(IW,'(3X,"AFQMC_Greens_Functions_Z")')
   ENDIF
   !
   !--- Initialize ---
   !
   G_Half = Zero
   G      = Zero 
   Det    = One
   !
   !Note we assume real wavefunction, so we treat Hetmite conjugate and transpose in the same manner
   !
   !---  S= [Vt.U*] ---
   !
   !     #
   !     #   SA = CWA+.CAT  <=> V+.U
   !     #
   !     ovlp = xp.dot(walker_batch.phia[iw].T, trial.psi0a.conj())
   DO I=1,NSpin
      NO = NOcc(I)
      IF (NO == 0) CYCLE
      !U_Cnj(1:NVar,1:NVar)  =  CONJG(U(1:NVar,1:NVar,I))
      U_Cnj(1:NVar,1:NO) => Wrk1(1:NVar*NO)
      S(1:NO,1:NO)       => Wrk2(1:NO**2) 
      SInv(1:NO,1:NO)    => Wrk3(1:NO**2) 
      U_Cnj(1:NVar,1:NO)  =  CONJG(U(1:NVar,1:NO,I))
      CALL ZGeMM('T','N',NO,NO,NVar,One,V(1,1,I),NVar,U_Cnj,NVar,Zero,S,NO)
      !
      !---  S-1= [Vt.U*]-1 ---
      !
      !     #
      !     #   SA-1 = [CWA+.CTA]^{-1}   <=> [V+.U]^{-1}
      !     #
      !     ovlp_inv = xp.linalg.inv(ovlp)
      CALL GetInvLR_Z(NO,S,SInv,DetX(I))
#if 0
      WRITE(IW,*)"V Matrix"
      CALL PrSqZ(V,NVar,NVar,NVar)
      WRITE(IW,*)"SA Matrix"
      CALL PrSqZ(S,NO,NO,NO)
      WRITE(IW,*)"SA-1 Matrix"
      CALL PrSqZ(SInv,NO,NO,NO)
#endif
       !
       !--- G_Half =  [(Vt.U*)^{-1}].[Vt] ---
       !
       !     #
       !     #  Green_Half_alpha = SA-1.CWA+  <=? [V+.U]^{-1}.V+
       !     #
       !     walker_batch.Ghalfa[iw] = xp.dot(ovlp_inv, walker_batch.phia[iw].T)
       !CALL ZGeMM('N','T',NOccA,NVar,NVar,One,InvSA,NOccA,V(1,1,1),NVar,Zero,G_Half(1,1,1),NVar)
       ! 
       !     [X_oo, 0_ov]      [ [V_oo]T,  [V_vo]T ]      [ X_oo.[V_oo]T,  X_oo.[V_vo]T ]
       !     [0_vo, 0_vv]   x  [ [V_ov]T,  [V_vv]T ]   =  [ 0_vo,          0_vv         ]
       !
       CALL ZGeMM('N','T',NO,NVar,NO,One,SInv,NO,V(1,1,I),NVar,Zero,G_Half(1,1,I),NVar)
       !  
       !   [GH_oo, GH_ov]   [ X_oo.[V_oo]T,  X_oo.[V_vo]T ]
       !   [ 0_vo, 0_vv ] = [ 0_vo,          0_vv         ]
       ! 
#if 0
       WRITE(IW,*)"G_Half  (Occ x NVar)  (Vir x NVar = 0)"
       CALL PrSqZ(G_Half(1,1,I),NVar,NVar,NVar)
#endif
       !
       !--- G =  [U*.(Vt.U*)^{-1}.Vt]+ ---
       !     #
       !     #  [Green_alpha]_pq = [CTA.SA-1.CWA+]_pq = [CWA.[CTA+.CWA]^{-1}/CTA]_qp
       !     #
       !     walker_batch.Ga[iw] = xp.dot(trial.psi0a.conj(), walker_batch.Ghalfa[iw])
       ! 
       !     [ U*_oo,  U*_ov ]    [GH_oo, GH_ov]      [ U*_oo.GH_oo, U*_oo.GH_ov]
       !     [ U*_vo,  U*_vv ]  x [ 0_vo,  0_vv]   =  [ U*_vo.GH_oo, U*_vo.GH_ov]
       !
       ! OpenQP port: the contraction index is NO, not NVar. U_Cnj is a pointer
       ! remap of only the first NVar*NO words of Wrk1, so K=NVar read NVar*(NVar-NO)
       ! uninitialized words past its end. The extra columns multiply rows
       ! NO+1..NVar of G_Half, which are zero, so the result was unchanged
       ! whenever the garbage happened to be finite -- but a NaN there poisons G
       ! through NaN*0. K=NO is both in bounds and exactly equivalent (and cheaper).
       CALL ZGeMM('N','N',NVar,NVar,NO,One,U_Cnj,NVar,G_Half(1,1,I),NVar,Zero,G(1,1,I),NVar)
       Det = Det*DetX(I)
   ENDDO
   !
   !
   !
   !CALL ZGeMM('T','N',NOccB,NOccB,NVar,One,V(1,1,2),NVar,U_Cnj(1,1,2),NVar,Zero,SB,NOccB)
   !CALL GetInvLR_Z(NOccB,SB,InvSB,DetB)
   !!CALL ZGeMM('N','T',NOccB,NVar,NVar,One,InvSB,NOccB,V(1,1,2),NVar,Zero,G_Half(1,1,2),NVar)
   !CALL ZGeMM('N','T',NOccB,NVar,NOccB,One,InvSB,NOccB,V(1,1,2),NVar,Zero,G_Half(1,1,2),NVar)
   !CALL ZGeMM('N','N',NOccB,NOccB,NVar,One,U_Cnj(1,1,2),NVar,G_Half(1,1,2),NVar,Zero,G(1,1,2),NVar)
   
   IF (main_rank .AND. PrintT%Walker > 1) THEN 
      WRITE(IW,'("Determinants of the overlap ")')
      WRITE(IW,'("Det (=DetA*DetB):",6F20.12)')Det,DetX(1:NSpin)
   ENDIF

   IF (main_rank .AND. PrintT%Walker > 2) THEN 
      WRITE(IW,'("Half rotated Greens function ")')
      DO I=1,NSpin
         CALL PrSqZ(G_Half(1,1,I),NVar,NVar,NVar)
      ENDDO
      WRITE(IW,'("Full rotated Greens function ")')
      DO I=1,NSpin
         CALL PrSqZ(G(1,1,I),NVar,NVar,NVar)
      ENDDO
   ENDIF
END SUBROUTINE 


SUBROUTINE AFQMC_Init_Wfn_D(NWlk,NVar,NSpin,NOccA,NOccB,U)
   USE IO_Files,     ONLY: IW
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   INTEGER,          INTENT(IN)  :: NWlk
   INTEGER,          INTENT(IN)  :: NVar
   INTEGER,          INTENT(IN)  :: NSpin
   INTEGER,          INTENT(IN)  :: NOccA
   INTEGER,          INTENT(IN)  :: NOccB
   DOUBLE PRECISION, INTENT(OUT) :: U(1:NVar,1:NVar,1:NSpin,1:NWlk)
   DOUBLE PRECISION, PARAMETER   :: Zero= 0D+00
   DOUBLE PRECISION, PARAMETER   :: One = 1D+00
   INTEGER                       :: ISpin
   INTEGER                       :: IMO
   INTEGER                       :: IWlk
 
   DO IWlk=1,NWlk
      DO ISpin=1,NSpin
         U(1:NVar,1:NVar,ISpin,IWlk) = Zero
         IF (ISpin==1) THEN
            DO IMO=1,NOccA
               U(IMO,IMO,ISpin,IWlk) = One
            ENDDO
         ELSEIF (ISpin==2) THEN
            DO IMO=1,NOccB
               U(IMO,IMO,ISpin,IWlk) = One
            ENDDO
         ENDIF
      ENDDO
   ENDDO
END SUBROUTINE

SUBROUTINE AFQMC_Init_Wfn_Z(NWlk,NVar,NSpin,NOcc,U)
   USE IO_Files,     ONLY: IW
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   INTEGER,          INTENT(IN)  :: NWlk
   INTEGER,          INTENT(IN)  :: NVar
   INTEGER,          INTENT(IN)  :: NSpin
   INTEGER,          INTENT(IN)  :: NOcc(1:NSpin)
   DOUBLE COMPLEX,   INTENT(OUT) :: U(1:NVar,1:NVar,1:NSpin,1:NWlk)
   DOUBLE PRECISION, PARAMETER   :: Zero= 0D+00
   DOUBLE PRECISION, PARAMETER   :: One = 1D+00
   INTEGER                       :: ISpin
   INTEGER                       :: IMO
   INTEGER                       :: IWlk
   INTEGER                       :: NO

   U(1:NVar,1:NVar,1:NSpin,1:NWlk) = CMPLX(Zero, Zero)
   DO IWlk=1,NWlk
      DO ISpin=1,NSpin
         NO = NOcc(ISpin)
         DO IMO=1,NO
            U(IMO,IMO,ISpin,IWlk)  = CMPLX(One, Zero)
         ENDDO
      ENDDO
   ENDDO
END SUBROUTINE


SUBROUTINE AFQMC_Orthogonalize_Walkers_QR_D(NVar,NSpin,U_old,U_new,Det)
!
!  see ipie/walkers/uhf_walkers.py
!
!  self.ovlp = trial.calc_greens_function(self)
!            = Det ([Green_alpha]_pq) 
!            = Det ([CTA.SA-1.CWA+]_pq) 
!            = Det ([CWA.[CTA+.CWA]^{-1}/CTA]_qp)
!
   USE IO_Files,     ONLY: IW
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NVar
   INTEGER, INTENT(IN) :: NSpin
   DOUBLE PRECISION, INTENT(IN)   :: U_old(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(OUT)  :: U_new(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(OUT)  :: Det
   DOUBLE PRECISION     :: Q(1:NVar,1:NVar)
   DOUBLE PRECISION     :: R(1:NVar,1:NVar)
   DOUBLE PRECISION     :: V(1:NVar,1:NVar)
   INTEGER              :: i
   INTEGER              :: ISpin

   Det = 1.D0
   DO ISpin=1,NSpin
      CALL QR_Decompostion_D(NVar,U_old(1,1,ISpin),Q,R)
      Det =  Det*R(i,i)
      U_new(1:NVar,1:NVar,ISpin) = Q(1:NVar,1:NVar)
   ENDDO

END SUBROUTINE

SUBROUTINE QR_Decompostion_D(n,H,Q,R)
!     The matrix Q is  represented  as  a  product  of  elementary
!     reflectors
!
!        Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!     Each H(i) has the form
!        H(i) = I - tau * v * v'
!
!     where tau is a real scalar, and v is a real vector with
!     v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is  stored  on  exit  in
!     A(i+1:m,i), and tau in TAU(i).
IMPLICIT NONE
DOUBLE PRECISION, PARAMETER :: Zero = 0D+00
DOUBLE PRECISION, PARAMETER :: One  = 1D+00
INTEGER,          INTENT(IN)  :: n
DOUBLE PRECISION, INTENT(IN)  :: H(1:n,1:n)
DOUBLE PRECISION, INTENT(OUT) :: Q(1:n,1:n)
DOUBLE PRECISION, INTENT(OUT) :: R(1:n,1:n)
INTEGER            :: i,j
INTEGER            :: m
INTEGER            :: lda
INTEGER            :: lwork
INTEGER            :: Info
DOUBLE PRECISION   :: Det
DOUBLE PRECISION   :: DetQ
DOUBLE PRECISION   :: A(1:n,1:n)
DOUBLE PRECISION   :: W(1:n,1:n)
DOUBLE PRECISION   :: T(1:n,1:n)
DOUBLE PRECISION   :: O(1:n,1:n)
DOUBLE PRECISION   :: V(1:n)
DOUBLE PRECISION   :: Tau(1:n)
DOUBLE PRECISION   :: Work(1:n*64)

m     =n
lda   =m
lwork =n*64

IF (LdA < MAX(1, m)) THEN
   ERROR STOP "Error in QR_Decompostion: LdA < MAX(1, m)"
ENDIF

A = H
CALL DGeQRF(m, n, a, lda, tau, work, lwork, info)
!CALL PrSq(A,n,n,n)

DO i=1,n
   DO j=1,i-1
      R(i,j)=Zero
   ENDDO
   DO j=i,n
      R(i,j)=A(i,j)
   ENDDO
ENDDO
write(6,*)"R-mat"
CALL PrSq(R,n,n,n)

O = Zero
DO i=1,n
   O(i,i) = One
ENDDO
write(6,*)"O-mat"
CALL PrSq(O,n,n,n)

T = O
DO i=1,n
   W        = O
   V(1:i-1) = Zero
   V(i)     = One
   V(i+1:n) = A(i+1:n,i)
   CALL DGEMM('N','T',n,n,1,-tau(i),V,n,V,n,One,W,n) ! I - tau*V.V+
   CALL DGEMM('N','N',n,n,n,One,T,n,W,n,Zero,Q,n)     ! Q = Ti.W
   T = Q                                              ! Ti+1 = Q
ENDDO
write(6,*)"T-mat"
CALL PrSq(T,n,n,n)
write(6,*)"Q-mat"
CALL PrSq(Q,n,n,n)

! OpenQP port: the upstream call named the non-existent generic GetInvLR. The
! arguments here are real, so it resolves to GetInvLR_D. Upstream only linked
! because this debug block sits in a routine the standalone build dead-strips;
! liboqp keeps every symbol, so the dangling reference had to be fixed.
CALL GetInvLR_D(n,Q,T,DetQ)
write(6,*)"Det Q=",DetQ
write(6,*)"Q-1"
CALL PrSq(T,n,n,n)
   CALL DGEMM('N','N',n,n,n,One,Q,n,T,n,Zero,O,n)     ! Q = Ti.W

write(6,*)"Q-1.Q"
CALL PrSq(O,n,n,n)
!DO i=1,n
!   IF (R(i,i) < Zero) O(i,i) = -One
!ENDDO
!write(6,*)"O-mat"
!CALL PrSq(O,n,n,n)
!CALL DGEMM('N','N',n,n,n,One,T,n,O,n,Zero,Q,n)        ! Q = Ti.W
END SUBROUTINE


SUBROUTINE AFQMC_Orthogonalize_Walkers_QR_Z(NVar,NSpin,NOcc,U_old,U_new,Det,Ovlp)
!
!  see ipie/walkers/uhf_walkers.py
!
!  self.ovlp = trial.calc_greens_function(self)
!            = Det ([Green_alpha]_pq) 
!            = Det ([CTA.SA-1.CWA+]_pq) 
!            = Det ([CWA.[CTA+.CWA]^{-1}/CTA]_qp)
!
   USE IO_Files,     ONLY: IW
   USE AFQMC_Module, ONLY: PrintT
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NVar
   INTEGER, INTENT(IN) :: NSpin
   INTEGER, INTENT(IN) :: NOcc(1:NSpin)
   DOUBLE COMPLEX,   PARAMETER    :: Zero = CMPLX(0D+00, 0D+00)
   DOUBLE COMPLEX,   PARAMETER    :: One  = CMPLX(1D+00, 0D+00)
   DOUBLE COMPLEX,   INTENT(IN)     :: U_old(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX,   INTENT(OUT)    :: U_new(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(OUT)    :: Det
   DOUBLE COMPLEX,   INTENT(INOUT)  :: Ovlp
   DOUBLE COMPLEX  :: Q(1:NVar,1:NVar)
   DOUBLE COMPLEX, POINTER  :: R(:,:)
   DOUBLE COMPLEX, TARGET   :: Wrk(1:NVar**2)
   DOUBLE COMPLEX  :: Phase
   INTEGER         :: i
   INTEGER         :: NO
   INTEGER         :: ISpin
   DOUBLE PRECISION  :: AbsR

   Det = 1.D0
   U_new = Zero
#if 0
   WRITE(IW,'(6X,"AFQMC_Orthogonalize_Walkers_QR_Z")')
   WRITE(IW,'(6X,"Input Walker.Det =",2F16.8)')Det
#endif
   DO ISpin=1,NSpin
      NO = NOcc(ISpin)
      IF (NO == 0) CYCLE
      R(1:NO,1:NO) => Wrk(1:NO**2)
      CALL QR_Decompostion_Z(NVar,NO,U_old(1,1,ISpin),Q,R)
      !
      !  Q: NVar x NVar
      !  R: NOcc x NOcc
      !
#if 0
      WRITE(IW,'("Q")')
      CALL PrSqZ(Q,NVar,NVar,NVar)
      WRITE(IW,'("R")')
      CALL PrSqZ(R,NO,NO,NO)
      CALL ZGeMM('N','N',NVar,NO,NO,One,Q,NVar,R,NO,Zero,V,NVar)
      WRITE(IW,'("Q.R = U_old")')
      CALL PrSqZ(V,NVar,NVar,NVar)
#endif
      DO i=1,NO
         AbsR = ABS(R(i,i))
         IF (AbsR <= TINY(1.D0)) THEN
            ERROR STOP "Singular walker encountered during QR stabilization"
         ENDIF
         ! Rotate each Q column so the corresponding R diagonal is real and
         ! positive.  The removed determinant is then a positive real scale,
         ! and no complex phase is lost when the stored overlap is rescaled.
         Phase = R(i,i)/CMPLX(AbsR,0.D0,KIND=KIND(R))
         U_new(1:NVar,i,ISpin) = Q(1:NVar,i)*Phase
         Det = Det*AbsR
      ENDDO
      !U_new(1:NVar,1:NVar,ISpin) = Q(1:NVar,1:NVar)
      !U_new(1:NVar,1:NO,ISpin) = Q(1:NVar,1:NO)
#if 0
      CALL ZGeMM('C','N',NVar,NVar,NVar,One,Q,NVar,Q,NVar,Zero,V,NVar)
      WRITE(IW,'("Check orthogonality Q+.Q = I")')
      CALL PrSqZ(V,NVar,NVar,NVar)
#endif
      !  Note:
      !       A = Q.R
      !    and always |det Q| = 1.
      !    Therefore, |det A| = |det Q| |det R| = |det R|
      !    R matrix is upper trianglar matrix so,
      !     det R = Prod_i  R_ii
#if 0
      WRITE(IW,'(6X,"ISpin =",I4)')ISpin
      WRITE(IW,'(6X,"Current      Det =",2F16.8)')DetX
#endif
      !Det = Det*DetX
   ENDDO
   !WRITE(IW,'(6X,"Input  Walker.Ovlp =",3F16.8)')Ovlp,ABS(Ovlp)
   Ovlp = Ovlp/CMPLX(Det,0D+00,KIND=KIND(Ovlp))
   !WRITE(IW,'(6X,"Output Walker.Det  =", F16.8)')Det
   !WRITE(IW,'(6X,"Output Walker.Ovlp =",4F16.8)')Ovlp,ABS(Ovlp),Theta

#if 0
   WRITE(IW,'(4X,"U_new")')
   DO ISpin=1,NSpin
      WRITE(IW,'(4X,"ISpin=",I4)')ISpin
      !NO = NOcc(ISpin)
      !U_new(1:NVar,NO+1:NVar,ISpin) = Zero
      CALL PrSqZ(U_new(1,1,ISpin),NVar,NVar,NVar)
   ENDDO
#endif

END SUBROUTINE

SUBROUTINE QR_Decompostion_Z(m,n,H,Q,R)
!
!    法線ベクトルnの（原点を通る）超平面に対して対称に点を写す鏡映変換をHausholder変換という。
!  この変換行列をHで表し、||n||=1とする。
!   H.x = x − 2(n, x)n = (I − 2n.nT).x
! これは直交行列であり、かつ対称行列である。
!     R = H_n. H_n-1. ... H_2. H_1. A 
!  とすれば
!     A = (H_n. H_n-1. ... H_2. H_1.)^{-1}.R
!       = H_1T.H_2T. ....  H_nT.R 
!       = H_1.H_2. ....  H_n.R
!       = Q.R 
!
!     The matrix Q is  represented  as  a  product  of  elementary
!     reflectors
!
!        Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!     Each H(i) has the form
!        H(i) = I - tau * v * v'
!
!     where tau is a real scalar, and v is a real vector with
!     v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is  stored  on  exit  in
!     A(i+1:m,i), and tau in TAU(i).
!
!  see https://qiita.com/Sharkkii/items/cebba1bc538fdaaa9fb5
   USE IO_Files, ONLY: IW
   IMPLICIT NONE
   DOUBLE COMPLEX,   PARAMETER :: Zero = CMPLX(0D+00, 0D+00)
   DOUBLE COMPLEX,   PARAMETER :: One  = CMPLX(1D+00, 0D+00)
   INTEGER,          INTENT(IN)  :: m
   INTEGER,          INTENT(IN)  :: n
   DOUBLE COMPLEX,   INTENT(IN)  :: H(1:m,1:m)
   DOUBLE COMPLEX,   INTENT(OUT) :: Q(1:m,1:m)
   DOUBLE COMPLEX,   INTENT(OUT) :: R(1:n,1:n)
   INTEGER            :: i,j
   INTEGER            :: lda
   INTEGER            :: lwork
   INTEGER            :: Info
   DOUBLE COMPLEX     :: Det
   DOUBLE COMPLEX     :: DetQ
   DOUBLE COMPLEX     :: A(1:m,1:m)
   DOUBLE COMPLEX     :: Tau(1:m)
   DOUBLE COMPLEX     :: Work(1:m*64)
   
   !m     =n
   
   
#if 0
   WRITE(IW,'("H-Mat")')
   CALL PrSqZ(H,m,m,m)
#endif
   A = H
   !
   !
   !
   ! ZGEQRF computes a QR factorization of a complex M-by-N matrix A:
   !
   !    A = Q * ( R ),
   !            ( 0 )
   !
   ! where:
   !    Q is a M-by-M orthogonal matrix;
   !    R is an upper-triangular N-by-N matrix;
   !    0 is a (M-N)-by-N zero matrix, if M > N.
   !
   !
   !        A is COMPLEX*16 array, dimension (LDA,N)
   !        On entry, the M-by-N matrix A.
   !        On exit, the elements on and above the diagonal of the array
   !        contain the min(M,N)-by-N upper trapezoidal matrix R (R is
   !        upper triangular if m >= n); the elements below the diagonal,
   !        with the array TAU, represent the unitary matrix Q as a
   !        product of min(m,n) elementary reflectors (see Further
   !        Details).
   !
   LdA = m
   lwork =m*64
   IF (LdA < MAX(1, n)) THEN
      ERROR STOP "Error in QR_Decompostion_Z: LdA < MAX(1, n)"
   ENDIF
   CALL ZGeQRF(m, n, A, LdA, tau, work, lwork, info)
   IF (Info /= 0) THEN
      ERROR STOP "ZGEQRF failed in QR_Decompostion_Z"
   ENDIF
#if 0
   WRITE(IW,'("A-Mat")')
   CALL PrSqZ(A,n,n,n)
#endif
   !   
   !   
   !   
   DO i=1,n
      DO j=1,i-1
         R(i,j)=Zero
      ENDDO
      DO j=i,n
         R(i,j)=A(i,j)
      ENDDO
   ENDDO
#if 0
   WRITE(IW,'("R-Mat")')
   CALL PrSqZ(R,n,n,n)
#endif
   !   
   !   
   !   
   ! Generate the first n columns of Q directly from LAPACK.  The previous
   ! hand-written loop consumed tau(n+1:m), although ZGEQRF only defines n
   ! reflectors for an m-by-n walker matrix.
   CALL ZUngQR(m,n,n,A,LdA,tau,work,lwork,info)
   IF (Info /= 0) THEN
      ERROR STOP "ZUNGQR failed in QR_Decompostion_Z"
   ENDIF
   Q = Zero
   Q(1:m,1:n) = A(1:m,1:n)
   RETURN

#if 0
   WRITE(IW,'("T-Mat")')
   CALL PrSqZ(T,n,n,n)
   WRITE(IW,'("Q-Mat")')
   CALL PrSqZ(Q,n,n,n)
#endif

#if 0
   WRITE(IW,'("Q+.Q = I")')
   CALL ZGEMM('C','N',n,n,n,1.D0,Q,n,Q,n,0.D0,T,n)    ! Q = Ti.W
   CALL PrSqZ(T,n,n,n)
#endif
   !
   !
   !
#if 0
   CALL GetInvLR_Z(n,Q,T,DetQ)
   WRITE(IW,'("Det Q = ",2F20.12)')DetQ
   WRITE(IW,'("Q-1")')
   CALL PrSqZ(T,n,n,n)
   CALL ZGEMM('N','N',n,n,n,One,Q,n,T,n,Zero,O,n)     
   WRITE(IW,'("Q-1.Q")')
   CALL PrSqZ(O,n,n,n)
   DO i=1,n
      !IF (R(i,i) < Zero) O(i,i) = -One
      IF (DBLE(R(i,i)) < 0D+0) O(i,i) = -One
   ENDDO
   WRITE(IW,'("O-Mat")')
   CALL PrSqZ(O,n,n,n)
   CALL ZGEMM('N','N',n,n,n,One,T,n,O,n,Zero,Q,n)        ! Q = Ti.W
#endif

END SUBROUTINE


SUBROUTINE PrEig(E,V,M,N,NDim)
!======================================================================!
!----------------------------------------------------------------------!
!                                                                      !
!     Function: Print eigensystem                                      !
!               V is N rows by M columns, with leading dimension NDim  !
!                                                                      !
!     Author(s): Motoyuki Uejima (2023)                                !
!                                                                      !
!     (c) Copyright MOLFEX, Inc. All rights reserved.                  !
!----------------------------------------------------------------------!
!======================================================================!
   IMPLICIT NONE
   INTEGER, PARAMETER :: IW = 6
   INTEGER,          INTENT(IN) :: M
   INTEGER,          INTENT(IN) :: N
   INTEGER,          INTENT(IN) :: NDim
   DOUBLE PRECISION, INTENT(IN) :: E(1:M)
   DOUBLE PRECISION, INTENT(IN) :: V(1:NDim,1:M)
   INTEGER, PARAMETER           :: NMax = 10
   INTEGER                      :: IMin, IMax, I, J
   IMax = 0
   DO
      IMin = IMax+1
      IMax = IMax+NMax
      IF (IMax > M) IMax = M
      WRITE (IW,'(1X)')
      WRITE (IW,'(16X,10(4X,I6,4X))')(I,I = IMin,IMax)
      !WRITE (Iw,'(2X,"Eigenvalue",2X,10F14.8)')(E(I),I = IMin,IMax)
      WRITE (IW,'(2X,"Eigenvalue",2X)',ADVANCE='NO')
      DO I=IMin,IMax
         IF (ABS(E(I)) < 1D+04) THEN
            WRITE (IW,'(F14.8)',ADVANCE='NO')E(I)
         ELSE
            WRITE (IW,'(E14.6)',ADVANCE='NO')E(I)
         ENDIF
      ENDDO
      WRITE (IW,'(A)')
      WRITE (IW,'(A)')
      DO J = 1,N
         WRITE (Iw,'(I8,6X,10F14.8)')J,(V(J,I),I = IMin,IMax)
      ENDDO
      IF (IMax >= M) EXIT
   ENDDO
END SUBROUTINE

SUBROUTINE GetInvLR_Z(n,A0,InvA,DetA)
!======================================================================!
!----------------------------------------------------------------------!
!                                                                      !
!     Function:  Calculate the inverse of a matrix and its determinant !
!                 using the LR decomposion                             !
!                                                                      !
!     Author(s): Motoyuki Uejima (2023)                                !
!                                                                      !
!     (c) Copyright MOLFEX, Inc. All rights reserved.                  !
!----------------------------------------------------------------------!
!======================================================================!
   IMPLICIT NONE
   INTEGER, PARAMETER :: IW = 6
   INTEGER,          INTENT(IN)  :: n
   DOUBLE COMPLEX,   INTENT(IN)  :: A0(1:n,1:n)
   DOUBLE COMPLEX,   INTENT(OUT) :: InvA(1:n,1:n)
   DOUBLE COMPLEX,   INTENT(OUT) :: DetA
   !--- Local ---
   DOUBLE COMPLEX                :: A(1:n,1:n)
   DOUBLE COMPLEX,   ALLOCATABLE :: Wrk(:)
   INTEGER                       :: m
   INTEGER                       :: i, Info, LdA, LWrk, LdtV, LdU
#if defined(IFORT) 
   !-- MKL---
   INTEGER                       :: IPiv(n)
#else
   INTEGER                       :: IPiv(n)
   !INTEGER                       :: IPiv(n)
#endif
   DOUBLE COMPLEX                :: tV(1:n,1:n)
   DOUBLE COMPLEX                :: S(1:n)
   DOUBLE COMPLEX                :: U(1:n,1:n)
   DOUBLE COMPLEX, PARAMETER     :: One  = CMPLX(1D+0,0D+0)
   DOUBLE COMPLEX, PARAMETER     :: Zero = CMPLX(0D+0,0D+0)

!#if defined(IFORT) 
!   write(0,*)"IFort"
!#else
!   write(0,*)"Not IFort"
!#endif
   m    = n
   LdA  = m
   LdU  = m
   LdtV = n
   LWrk = 64*n
   A    = A0
  
   ALLOCATE(Wrk(1:n*LWrk))
 
   IPiv  = 0
   InvA  = A0
   LdA   = n
   Info  = 0
   CALL ZGeTrF(n,n,InvA,LdA,IPiv,Info)
#if 0
   WRITE(IW,'("Pivot")')
   DO i=1, n
      WRITE(IW,'(2I6)')i,IPiv(i)
   ENDDO
#endif
   IF (Info /= 0) THEN
      WRITE(IW,'(2X,"ZGeTrF: The factor U is singular")')
      ERROR STOP
   ENDIF
 
   DetA = One
   DO i=1,n
      DetA = DetA*InvA(i,i)
   ENDDO
 
   DO i=1,n
      IF (IPiv(i) /= i) DetA = -DetA
   ENDDO
#if 0
   WRITE(IW,'("DetA=",2F20.12)')DetA
#endif

#if DEBUG
   WRITE(IW,'("Pivot")')
   DO i=1, n
      WRITE(IW,'(2I6)')i,IPiv(i)
   ENDDO
#endif

   Wrk  = Zero
   Info = 0
   LdA  = n
   CALL ZGeTrI(n,InvA,LdA,IPiv,Wrk,LWrk,Info)
 
   IF (Info /= 0) THEN
      WRITE(IW,'("Failed in generating inversion matrix")')
      ERROR STOP
   ENDIF

#if DEBUG
   WRITE(IW,'("A")')
   CALL PrSqZ(A0,n,n,n)
   WRITE(IW,'("A-1")')
   CALL PrSqZ(InvA,n,n,n)
   WRITE(IW,'("Check Inversion A.A-1")')
   CALL PrSqZ(MatMul(A0,InvA),n,n,n)
#endif

   DEALLOCATE(Wrk)
END SUBROUTINE


SUBROUTINE GetInvLR_D(n,A0,InvA,DetA)
!======================================================================!
!----------------------------------------------------------------------!
!                                                                      !
!     Function:  Calculate the inverse of a matrix and its determinant !
!                 using the LR decomposion                             !
!                                                                      !
!     Author(s): Motoyuki Uejima (2023)                                !
!                                                                      !
!     (c) Copyright MOLFEX, Inc. All rights reserved.                  !
!----------------------------------------------------------------------!
!======================================================================!
   IMPLICIT NONE
   INTEGER, PARAMETER :: IW = 6
   INTEGER,          INTENT(IN)  :: n
   DOUBLE PRECISION, INTENT(IN)  :: A0(1:n,1:n)
   DOUBLE PRECISION, INTENT(OUT) :: InvA(1:n,1:n)
   DOUBLE PRECISION, INTENT(OUT) :: DetA
   !--- Local ---
   DOUBLE PRECISION              :: A(1:n,1:n)
   DOUBLE PRECISION, ALLOCATABLE :: Wrk(:)
   INTEGER                       :: m
   INTEGER                       :: i, Info, LdA, LWrk, LdtV, LdU
#if defined(IFORT) 
   !-- MKL---
   INTEGER                       :: IPiv(n)
#else
   INTEGER                       :: IPiv(n)
   !INTEGER                       :: IPiv(n)
#endif
   DOUBLE PRECISION              :: tV(1:n,1:n)
   DOUBLE PRECISION              :: S(1:n)
   DOUBLE PRECISION              :: U(1:n,1:n)

#if defined(IFORT) 
   write(0,*)"IFort"
#else
   write(0,*)"Not IFort"
#endif
   m    = n
   LdA  = m
   LdU  = m
   LdtV = n
   LWrk = 64*n
   A    = A0
  
   ALLOCATE(Wrk(1:n*LWrk))
 
   IPiv  = 0
   InvA  = A0
   LdA   = n
   Info  = 0
   CALL DGeTrF(n,n,InvA,LdA,IPiv,Info)
#if 0
   WRITE(IW,'("Pivot")')
   DO i=1, n
      WRITE(IW,'(2I6)')i,IPiv(i)
   ENDDO
#endif
!
 
   IF (Info /= 0) THEN
      WRITE(IW,'(2X,"DGeTrF: The factor U is singular")')
      ERROR STOP
   ENDIF
 
   DetA = 1.D0
   DO i=1,n
      DetA   = DetA*InvA(i,i)
   ENDDO
 
   DO i=1,n
      IF (IPiv(i) /= i) DetA = -DetA
   ENDDO
   !WRITE(6,*)"DetA=",DetA

#if DEBUG
   WRITE(IW,'("Pivot")')
   DO i=1, n
      WRITE(IW,'(2I6)')i,IPiv(i)
   ENDDO
#endif

   Wrk  = 0.D0
   Info = 0
   LdA  = n
   CALL DGeTrI(n,InvA,LdA,IPiv,Wrk,LWrk,Info)
 
   IF (Info /= 0) THEN
      WRITE(IW,'("Failed in generating inversion matrix")')
      ERROR STOP
   ENDIF

#if DEBUG
   WRITE(IW,'("InvA")')
   CALL PrSq(InvA,n,n,n)
   WRITE(IW,'("Check Inversion A0.InvA")')
   CALL PrSq(MatMul(A0,InvA),n,n,n)
#endif

   DEALLOCATE(Wrk)
END SUBROUTINE


SUBROUTINE ExpMat(n,A,ExpA)
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: n
  DOUBLE PRECISION, INTENT(IN)  :: A(1:n,1:n)
  DOUBLE PRECISION, INTENT(OUT) :: ExpA(1:n,1:n)
  ! OpenQP port: local allocatables in place of the Gellan global scratch.
  DOUBLE PRECISION, ALLOCATABLE :: E(:)
  DOUBLE PRECISION, ALLOCATABLE :: V(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: W(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: ExpE(:,:)
  INTEGER :: i
  ALLOCATE(E(1:n),V(1:n,1:n),W(1:n,1:n),ExpE(1:n,1:n))
  CALL diagMRRR(n,A,E,V)

  ExpE=0.D0
  DO i=1,n
     ExpE(i,i) = DEXP(E(i))
  ENDDO
  !CALL PrSq(ExpE,n,n,n)
  CALL DGeMM('N','N',n,n,n,1.D0,V,n,ExpE,n,0.D0,W,n)
  CALL DGeMM('N','T',n,n,n,1.D0,W,n,V,n,0.D0,ExpA,n)
  !     A.V = V.E
  !  V+.A.V = E
  !  V+.A^n.V = (V+.A.V)^n = E^n
  !  A^n = V.E^n.V+
  DEALLOCATE(E,V,W,ExpE)
END SUBROUTINE



SUBROUTINE PrSqZ(A,M,N,LdA)
!======================================================================!
!----------------------------------------------------------------------!
!                                                                      !
!     Function: Print square matrix A                                  !
!               A is N rows by M columns, with leading dimension NDim  !
!                                                                      !
!     Author(s): Motoyuki Uejima (2023)                                !
!                                                                      !
!----------------------------------------------------------------------!
!======================================================================!
   USE IO_Files,     ONLY: IW
   IMPLICIT DOUBLE PRECISION(A-H,O-Z)
   INTEGER,          INTENT(IN) :: M
   INTEGER,          INTENT(IN) :: N
   INTEGER,          INTENT(IN) :: LdA
   DOUBLE COMPLEX,   INTENT(IN) :: A(1:LdA,1:M)
   INTEGER,          PARAMETER  :: NMax = 5

   IMax = 0
   DO
      IMin = IMax+1
      IMax = IMax+NMax
      IF (IMax > M) IMax = M
      WRITE (IW,'(1X)')
      WRITE (IW,'(10X,5(7X,I6,15X))')(I,I = IMin,IMax)
      DO J = 1,N
         WRITE (IW,'(I5,1X)',ADVANCE='NO')J
         DO I=IMin,IMax
            IF (ABS(A(J,I)) < 1D+03) THEN
               WRITE (IW,'(F14.7)',ADVANCE='NO')REAL(A(J,I))
               WRITE (IW,'(F13.7,"i")',ADVANCE='NO')AIMAG(A(J,I))
            ELSE
               WRITE (IW,'(ES14.6)',ADVANCE='NO')REAL(A(J,I))
               WRITE (IW,'(ES13.6,"i")',ADVANCE='NO')AIMAG(A(J,I))
            ENDIF
         ENDDO
         WRITE (IW,'(A)')
      ENDDO
      IF (IMax >= M) EXIT
   ENDDO
END SUBROUTINE

SUBROUTINE ExpMat_Z(n,nterm,A,ExpA)
!======================================================================!
!----------------------------------------------------------------------!
!                                                                      !
!     Function: ExpA = Exp(A)                                          !
!                                                                      !
!     Author(s): Motoyuki Uejima (2017)                                !
!                                                                      !
!     (c) Copyright Ten-no Research Group. All rights reserved.        !
!----------------------------------------------------------------------!
!======================================================================!
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: n
  INTEGER,          INTENT(IN)  :: nterm
  DOUBLE COMPLEX,   INTENT(IN)  :: A(1:n,1:n)
  DOUBLE COMPLEX,   INTENT(OUT) :: ExpA(1:n,1:n)
  DOUBLE COMPLEX                :: I(1:n,1:n)
  DOUBLE COMPLEX                :: Wrk(1:n,1:n)
  !DOUBLE PRECISION, PARAMETER   :: Zero=0D+00, One=1D+00
  DOUBLE COMPLEX, PARAMETER     :: Zero = CMPLX(0D+00, 0D+00)
  DOUBLE COMPLEX, PARAMETER     :: One = CMPLX(1D+00, 0D+00)
  INTEGER                       :: j

  I = Zero
  DO j=1,n
    I(j,j) = One
  ENDDO
  ExpA = I

  Wrk = I
  DO j=nterm,1,-1
     CALL ZGeMM('N','N',n,n,n,One/CMPLX(DBLE(j),0D+00),A,n,ExpA,n,One,Wrk,n)
     !
     !  Wrk = I + A/n
     !  Wrk = I + (I + A/n)A/n-1 = I + A/n-1 +  A^2/n(n-1)
     !  Wrk = I + (I + A/n-1 +  A^2/n(n-1))A/n-2 = I + A/n-2 +  A^2/(n-1)(n-2) + A^3/(n-1)(n-2)(n-3)
     !
     ExpA = Wrk
     Wrk  = I
  ENDDO

END SUBROUTINE

SUBROUTINE AFQMC_Overlap_Z(NVar,NSpin,NOcc,U,V,Det)
   !
   USE IO_Files,     ONLY: IW
   USE AFQMC_Module, ONLY: PrintT
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   INTEGER,        INTENT(IN)  :: NVar
   INTEGER,        INTENT(IN)  :: NSpin
   INTEGER,        INTENT(IN)  :: NOcc(1:NSpin)
   DOUBLE COMPLEX, INTENT(IN)  :: U(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX, INTENT(IN)  :: V(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX, INTENT(OUT) :: Det
   DOUBLE COMPLEX, POINTER     :: S(:,:)
   DOUBLE COMPLEX              :: DetX(1:NSpin)
   DOUBLE COMPLEX              :: U_Cnj(1:NVar,1:NVar)
   DOUBLE COMPLEX, TARGET      :: Wrk1(1:NVar**2)
   DOUBLE COMPLEX, PARAMETER   :: Zero = CMPLX(0D+00, 0D+00)
   DOUBLE COMPLEX, PARAMETER   :: One  = CMPLX(1D+00, 0D+00)
   INTEGER                     :: IPiv(1:NVar)
   INTEGER                     :: I, J, NO, Info
   !
   !  U: Trial
   !  V: Walker
   !
   Det = One
   DO I=1,NSpin
      NO = NOcc(I)
      IF (NO == 0) CYCLE
      !U_Cnj(1:NVar,1:NVar)  =  CONJG(U(1:NVar,1:NVar,I))
      U_Cnj(1:NVar,1:NO)  =  CONJG(U(1:NVar,1:NO,I))
      S(1:NO,1:NO)    => Wrk1(1:NO**2) 
      CALL ZGeMM('T','N',NO,NO,NVar,One,V(1,1,I),NVar,U_Cnj,NVar,Zero,S,NO)
      ! Only the determinant is needed here.  Avoid the old inverse call so a
      ! trial component with exactly zero overlap can contribute zero to a
      ! multi-determinant MRSF-CIS overlap without aborting the run.
      CALL ZGeTrF(NO,NO,S,NO,IPiv,Info)
      IF (Info < 0) THEN
         ERROR STOP "Invalid argument passed to ZGETRF in AFQMC_Overlap_Z"
      ELSEIF (Info > 0) THEN
         DetX(I) = Zero
      ELSE
         DetX(I) = One
         DO J=1,NO
            DetX(I) = DetX(I)*S(J,J)
            IF (IPiv(J) /= J) DetX(I) = -DetX(I)
         ENDDO
      ENDIF
      Det = Det*DetX(I)
   ENDDO
   ! 
   ! ndown = walkers.ndown
   ! ovlp_a = xp.einsum("wmi,mj->wij", walkers.phia, trial.psi0a.conj(), optimize=True)
   ! 
   !       S(w,i,j) = V(w,m,i).U*(m,j)
   !                = [VT.U*]w,ij
   ! 
   ! sign_a, log_ovlp_a = xp.linalg.slogdet(ovlp_a)
   ! 
   ! if ndown > 0 and not walkers.rhf:
   !     ovlp_b = xp.einsum("wmi,mj->wij", walkers.phib, trial.psi0b.conj(), optimize=True)
   !     sign_b, log_ovlp_b = xp.linalg.slogdet(ovlp_b)
   !     # log_shift=0
   !     #print(f"log_shift={walkers.log_shift}")
   !     ot = sign_a * sign_b * xp.exp(log_ovlp_a + log_ovlp_b - walkers.log_shift)
   ! elif ndown > 0 and walkers.rhf:
   !     ot = sign_a * sign_a * xp.exp(log_ovlp_a + log_ovlp_a - walkers.log_shift)
   ! elif ndown == 0:
   !     ot = sign_a * xp.exp(log_ovlp_a - walkers.log_shift)
   ! 
   ! synchronize()
   ! 
   ! return ot
   ! 
END SUBROUTINE 
