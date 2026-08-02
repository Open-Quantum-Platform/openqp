!> @file afqmc_hamiltonian.F90
!>
!> @brief AFQMC Fock build and energy estimators (vendored from
!>        Open-Quantum-Platform/AFQMC).
!>
!> The upstream file also carries Gellan integral-file loaders and the
!> Gellan-side Hamiltonian constructor.  Those are dropped here: OpenQP builds
!> the MO core Hamiltonian and the Cholesky vectors in memory from the
!> tagarray, in afqmc_oqp_hamiltonian.F90.


SUBROUTINE AFQMC_Construct_Fock(NVar,NSpin,NChol,NOcc,Hc,Chol,Fc)
   USE IO_Files,     ONLY: IW
   USE MPI_Parallel, ONLY: main_rank
   USE AFQMC_Module, ONLY: PrintT
   IMPLICIT NONE
   INTEGER,          INTENT(IN)  :: NVar
   INTEGER,          INTENT(IN)  :: NSpin
   INTEGER,          INTENT(IN)  :: NChol
   INTEGER,          INTENT(IN)  :: NOcc(1:NSpin)
   DOUBLE PRECISION, INTENT(IN)  :: Hc(1:NVar*(NVar+1)/2)
   DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
   DOUBLE PRECISION, INTENT(OUT) :: Fc(1:NVar,1:NVar,1:NSpin)
   INTEGER                       :: i,j,k,m
   IF (main_rank .AND. PrintT%Integral > 1) THEN
      WRITE(IW,'("AFQMC_Form_Fock")')
   ENDIF
   IF (NSpin==1) THEN  
      CALL Expnd(Hc,Fc(1,1,1),NVar,0)
      DO m=1,NChol
         DO j=1,NVar
            DO i=1,NVar
               DO k=1,NOcc(1)
                 ! 2(ij|kk) - (ik|jk)
                  Fc(i,j,1) = Fc(i,j,1) + 2.D0*Chol(i,j,m)*Chol(k,k,m)
                  Fc(i,j,1) = Fc(i,j,1) - 1.D0*Chol(i,k,m)*Chol(j,k,m)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      IF (main_rank .AND. PrintT%Integral > 1) THEN
         WRITE(IW,'("Fock Matrix")')
         CALL PrSq(Fc(1,1,1),NVar,NVar,NVar)
      ENDIF
   ELSE
   ENDIF
END SUBROUTINE 
!============
!
!   Evaluator
!
!============
SUBROUTINE AFQMC_Mean_Field_Energy(NVar,NChol,NEl,ENu,Hc,Chol)
   USE IO_Files,     ONLY: IW
   USE MPI_Parallel, ONLY: main_rank
   USE AFQMC_Module, ONLY: PrintT
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NVar
   INTEGER, INTENT(IN) :: NChol
   INTEGER, INTENT(IN) :: NEl
   DOUBLE PRECISION, INTENT(IN)  :: ENu
   DOUBLE PRECISION, INTENT(IN)  :: Hc(1:NVar*(NVar+1)/2)
   DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
   DOUBLE PRECISION, EXTERNAL    :: DDOT
   INTEGER              :: i,j,m,ii
   DOUBLE PRECISION     :: E0, E1, E2
   DOUBLE PRECISION     :: EJ, EK
   INTEGER              :: NOcc
   !
   !E2 = 2*(ij|ij) - (ij|ji)
   !
   NOcc = NEl/2
   E0 = ENu
   E1 = 0.D0
   ii = 0
   DO i=1,NOcc
      ii = ii + i
      E1 = E1 + 2.D0*Hc(ii)
   ENDDO

   E2 = 0.D0
   EJ = 0.D0
   EK = 0.D0
   DO m=1,NChol
     DO i=1,NOcc
       DO j=1,NOcc
            EJ = EJ + Chol(i,i,m)*Chol(j,j,m)
            EK = EK + Chol(i,j,m)*Chol(j,i,m)
         ENDDO
      ENDDO
   ENDDO
   E2 = 2.D0*EJ - EK
   IF (main_rank .AND. PrintT%Energy > 1) THEN
      WRITE(IW,'("EJ   =",F20.12)')EJ
      WRITE(IW,'("EK   =",F20.12)')EK
      WRITE(IW,'("E0   =",F20.12)')E0
      WRITE(IW,'("E1   =",F20.12)')E1
      WRITE(IW,'("E2   =",F20.12)')E2
      WRITE(IW,'("ETot =",F20.12)')E0+E1+E2
   ENDIF
END SUBROUTINE 

SUBROUTINE AFQMC_Energy_Estimate(NVar,NSpin,NOcc,NChol,NWlk,Walker,&
&                                NDet,Trial_Coeff,Trial_Phi,ENu,Hc,Chol,E_Shift)
   USE AFQMC_Module, ONLY: AFQMC_Walker_Type
! OpenQP port: the upstream #ifdef AFQMC_STANDALONE selected between Gellan's
! real MPI layer and a serial stub set. liboqp always uses the serial stubs in
! afqmc_host.F90 (it parallelises with OpenMP), so the guard is resolved here.
   USE MPI_Parallel, ONLY: my_rank,NProc,MPI_Comm_World,MPI_SUM,Gel_AllReduce
   IMPLICIT NONE
   INTEGER,          INTENT(IN)  :: NVar
   INTEGER,          INTENT(IN)  :: NSpin
   INTEGER,          INTENT(IN)  :: NOcc(1:NSpin)
   INTEGER,          INTENT(IN)  :: NChol
   INTEGER,          INTENT(IN)  :: NWlk
   INTEGER,          INTENT(IN)  :: NDet
   DOUBLE PRECISION, INTENT(IN)  :: ENu
   DOUBLE COMPLEX,   INTENT(IN)  :: Trial_Coeff(1:NDet)
   DOUBLE COMPLEX,   INTENT(IN)  :: Trial_Phi(1:NVar,1:NVar,1:NSpin,1:NDet)
   TYPE(AFQMC_Walker_Type), INTENT(INOUT)  :: Walker(1:NWlk)
   DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
   DOUBLE PRECISION, INTENT(IN)  :: Hc(1:NVar*(NVar+1)/2)
   DOUBLE COMPLEX,   INTENT(OUT) :: E_Shift
   DOUBLE COMPLEX                :: E_Loc
   DOUBLE PRECISION              :: Total_Weight
   DOUBLE PRECISION              :: Dummy
   DOUBLE COMPLEX                :: ZDummy
   INTEGER                       :: IWlk
   
   ! OpenQP port: the per-walker local energies are the second-largest cost in a
   ! run (a full multideterminant contraction each) and were evaluated serially
   ! while the propagation itself was threaded. They are computed in parallel
   ! here and accumulated afterwards in walker order, so the sum is bit-identical
   ! whatever OMP_NUM_THREADS is -- a reduction clause would not be, and the
   ! benchmarks pin thread-independence to 1e-12.
   E_Shift = (0D+00,0D+00)
   Total_Weight = 0.D0
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(IWlk,E_Loc)
   DO IWlk = my_rank + 1, NWlk, NProc
      CALL AFQMC_MultiDet_Local_Energy_Z(NVar,NSpin,NOcc,NChol,NDet,&
     &             Trial_Coeff,Trial_Phi,Walker(IWlk)%UW,ENu,Hc,Chol,E_Loc)
      Walker(IWlk)%E_Loc = E_Loc
   ENDDO
!$OMP END PARALLEL DO
   DO IWlk = my_rank + 1, NWlk, NProc
      Total_Weight =  Total_Weight + Walker(IWlk)%Weight
      E_Shift = E_Shift + Walker(IWlk)%Weight*Walker(IWlk)%E_Loc
   ENDDO
   CALL Gel_AllReduce(E_Shift,ZDummy,1,'Z',MPI_SUM,MPI_Comm_World,1)
   CALL Gel_AllReduce(Total_Weight,Dummy,1,'D',MPI_SUM,MPI_Comm_World,1)
   E_Shift = E_Shift/Total_Weight
END SUBROUTINE

SUBROUTINE AFQMC_Energy_Estimator_Z(NVar,NSpin,NChol,G,G_Half,ENu,Hc,Chol,E,HalfRot)
   !
   !  G_pq = <Psi_T|ap+.aq|Psi_W>/<Psi_T|Psi_W> = [CW.(CT+.CW)-1.CT+]_qp
   !
   !           |Psi_T> = U|Psi>  (U=CT)
   !           |Psi_W> = V|Psi>  (V=CW)
   !
   !  Theta_pq = [CW.(CT+.CW)-1]_pq = [V.(U+.V)-1]_pq
   !
   !     であるとすれば
   !
   !        <Psi_T|a_p+.a_q|Psi_W>/<Psi_T|Psi_W> h_pq
   !        = <Psi|U+.a_p+.a_q.V|Psi> / <Psi|U+.V|Psi> h_pq
   !        = <Psi|[U+.a_p+.U].U+.V.[V+.a_q.V]|Psi> / <Psi|U+.V|Psi> h_pq
   !        = U*_pi V_qj <Psi|a_i+.[U+.V].a_j|Psi> / <Psi|U+.V|Psi> h_pq
   !        = U*_pi V_qj [(U+.V)-1]_ji h_pq 
   !        = h_iq. Theta_qi 
   !
   !     #
   !     #  [Green_half]_pq = [U.S-1.V+]_pq = [U.[U+.V]^{-1}.V+]_qp
   !     #   Theta_pq       = [V.(U+.V)-1]_pq = GH_qp
   !
   !        <Psi_T|a_p+.a_q+.a_s.a_r|Psi_W>/<Psi_T|Psi_W> (pr|qs)
   !        = <Psi|U+.a_p+.a_q+.a_s.a_r.V|Psi> / <Psi|U+.V|Psi> (pr|qs)
   !        = <Psi|[U+.a_p+.U].[U+.a_q+.U].U+.V.[V+.a_s.V].[V+.a_r.V]|Psi> / <Psi|U+.V|Psi> (pr|qs)
   !        = U*_pi U*_qj V_sl V_rk <Psi|a_i+.a_j+.[U+.V].a_l.a_k|Psi> / <Psi|U+.V|Psi> (pr|qs)
   !        = U*_pi U*_qj V_sl V_rk ( [(U+.V)-1]_ki  [(U+.V)-1]_lj - [(U+.V)-1]_kj [(U+.V)-1]_li) (pr|qs) 
   !        = (ir|js) ( [V.(U+.V)-1]_ri [V.(U+.V)-1]_sj - [V.(U+.V)-1]_rj [V.(U+.V)-1]_si )
   !        = (ir|js) ( θ_ri θ_sj - θ_rj θ_si )
   !        = (ir|js) ( GH_ir GH_js - GH_jr GH_is )
   !
   USE IO_Files,     ONLY: IW
   USE AFQMC_Module, ONLY: PrintT
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NVar
   INTEGER, INTENT(IN) :: NSpin
   INTEGER, INTENT(IN) :: NChol
   DOUBLE COMPLEX,   INTENT(IN)  :: G(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX,   INTENT(IN)  :: G_Half(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(IN)  :: ENu
   DOUBLE PRECISION, INTENT(IN)  :: Hc(1:NVar*(NVar+1)/2)
   DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
   DOUBLE COMPLEX,   INTENT(OUT) :: E 
   LOGICAL,          INTENT(IN)  :: HalfRot
   INTEGER                       :: i,j,k,l,m,ii,ij,ispin
   DOUBLE COMPLEX                :: E0, E1, E2, EJ, EK
   DOUBLE COMPLEX                :: X
   DOUBLE COMPLEX                :: Y(1:NVar,1:NVar)
   DOUBLE COMPLEX                :: Z(1:NVar,1:NVar)
   DOUBLE PRECISION, PARAMETER   :: Zero= 0D+00
   DOUBLE PRECISION, PARAMETER   :: One = 1D+00
   DOUBLE PRECISION, PARAMETER   :: Two = 2D+00
   DOUBLE COMPLEX,   PARAMETER   :: ZZero = (0D+00,0D+00)
   DOUBLE COMPLEX,   PARAMETER   :: ZHalf = (5D-01,0D+00)
   DOUBLE COMPLEX,   PARAMETER   :: ZTwo  = (2D+00,0D+00)

   E0 = ENu
   E1 = ZZero
   E2 = ZZero
   EJ = ZZero
   EK = ZZero
   IF (HalfRot) THEN
      DO ispin=1,NSpin
         ij = 0
         DO i=1,NVar
            DO j=1,i-1
               ij = ij +1
               E1 = E1 + Hc(ij)*G(j,i,ISpin)
               E1 = E1 + Hc(ij)*G(i,j,ISpin)
            ENDDO
            ij = ij +1
            ii = ij
            E1 = E1 + Hc(ii)*G(i,i,ISpin)
         ENDDO
      ENDDO

         ! (ij|kl) GH_ij GH_kl
      DO m=1,NChol
         X= ZZero
         DO ispin=1,NSpin
            DO i=1,NVar
               DO j=1,NVar
                  X = X + Chol(i,j,m)*G(i,j,ispin)
               ENDDO
            ENDDO
         ENDDO
         EJ = EJ + ZHalf*X**2
      ENDDO
      !
      ! (ij|kl) GH_kj GH_il  
      ! j-spin = k-spin <= GH_jk
      ! i-spin = l-spin <= GH_il
      !
      DO ispin=1,NSpin
         DO m=1,NChol
            Y = ZZero
            DO i=1,NVar
               DO j=1,NVar
                  DO k=1,NVar
                     Y(i,k) = Y(i,k) + Chol(i,j,m)*G(k,j,ispin)
                  ENDDO
               ENDDO
            ENDDO
            Z = ZZero
            DO i=1,NVar
               DO k=1,NVar
                  DO l=1,NVar
                     Z(k,i) = Z(k,i) + Chol(k,l,m)*G(i,l,ispin)
                  ENDDO
               ENDDO
            ENDDO
            DO i=1,NVar
               DO k=1,NVar
                  EK = EK + Y(i,k)*Z(k,i)*ZHalf
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF
   E2 = EJ - EK
   E  = E0+E1+E2

   IF (main_rank .AND. PrintT%Energy > 2) THEN 
      WRITE(IW,'(3X,"AFQMC_Energy_Estimator_Z")')
      WRITE(IW,'(3X,"EJ    =",2F20.12)')EJ
      WRITE(IW,'(3X,"EK    =",2F20.12)')EK
      WRITE(IW,'(3X,"E0    =",2F20.12)')E0
      WRITE(IW,'(3X,"E1    =",2F20.12)')E1
      WRITE(IW,'(3X,"E0+E1 =",2F20.12)')E0+E1
      WRITE(IW,'(3X,"E2    =",2F20.12)')E2
      WRITE(IW,'(3X,"E_Loc =",2F20.12)')E
   ENDIF
END SUBROUTINE 

SUBROUTINE AFQMC_Energy_Estimator_D(NVar,NSpin,NChol,GH,GF,ENu,Hc,Chol,E,HalfRot)
   !
   !  G_pq = <Psi_T|ap+.aq|Psi_W>/<Psi_T|Psi_W> = [CW.(CT+.CW)-1.CT+]_qp
   !
   !           |Psi_T> = U|Psi>  (U=CT)
   !           |Psi_W> = V|Psi>  (V=CW)
   !
   !  Theta_pq = [CW.(CT+.CW)-1]_pq = [V.(U+.V)-1]_pq
   !
   !     であるとすれば
   !
   !        <Psi_T|a_p+.a_q|Psi_W>/<Psi_T|Psi_W> h_pq
   !        = <Psi|U+.a_p+.a_q.V|Psi> / <Psi|U+.V|Psi> h_pq
   !        = <Psi|[U+.a_p+.U].U+.V.[V+.a_q.V]|Psi> / <Psi|U+.V|Psi> h_pq
   !        = U*_pi V_qj <Psi|a_i+.[U+.V].a_j|Psi> / <Psi|U+.V|Psi> h_pq
   !        = U*_pi V_qj [(U+.V)-1]_ji h_pq 
   !        = h_iq. Theta_qi 
   !
   !     #
   !     #  [Green_half]_pq = [U.S-1.V+]_pq = [U.[U+.V]^{-1}.V+]_qp
   !     #   Theta_pq       = [V.(U+.V)-1]_pq = GH_qp
   !
   !        <Psi_T|a_p+.a_q+.a_s.a_r|Psi_W>/<Psi_T|Psi_W> (pr|qs)
   !        = <Psi|U+.a_p+.a_q+.a_s.a_r.V|Psi> / <Psi|U+.V|Psi> (pr|qs)
   !        = <Psi|[U+.a_p+.U].[U+.a_q+.U].U+.V.[V+.a_s.V].[V+.a_r.V]|Psi> / <Psi|U+.V|Psi> (pr|qs)
   !        = U*_pi U*_qj V_sl V_rk <Psi|a_i+.a_j+.[U+.V].a_l.a_k|Psi> / <Psi|U+.V|Psi> (pr|qs)
   !        = U*_pi U*_qj V_sl V_rk ( [(U+.V)-1]_ki  [(U+.V)-1]_lj - [(U+.V)-1]_kj [(U+.V)-1]_li) (pr|qs) 
   !        = (ir|js) ( [V.(U+.V)-1]_ri [V.(U+.V)-1]_sj - [V.(U+.V)-1]_rj [V.(U+.V)-1]_si )
   !        = (ir|js) ( θ_ri θ_sj - θ_rj θ_si )
   !        = (ir|js) ( GH_ir GH_js - GH_jr GH_is )
   !
   USE IO_Files,     ONLY: IW
   USE AFQMC_Module, ONLY: PrintT
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NVar
   INTEGER, INTENT(IN) :: NSpin
   INTEGER, INTENT(IN) :: NChol
   DOUBLE PRECISION, INTENT(IN)  :: GF(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(IN)  :: GH(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(IN)  :: ENu
   DOUBLE PRECISION, INTENT(IN)  :: Hc(1:NVar*(NVar+1)/2)
   DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
   DOUBLE PRECISION, INTENT(OUT) :: E 
   LOGICAL,          INTENT(IN)  :: HalfRot
   INTEGER                       :: i,j,k,l,m,ii,ij,ispin
   DOUBLE PRECISION               :: E0, E1, E2, EJ, EK
   DOUBLE PRECISION              :: X
   DOUBLE PRECISION              :: Y(1:NVar,1:NVar)
   DOUBLE PRECISION              :: Z(1:NVar,1:NVar)
   DOUBLE PRECISION, PARAMETER   :: Zero = 0D+00
   DOUBLE PRECISION, PARAMETER   :: One  = 1D+00
   DOUBLE PRECISION, PARAMETER   :: Two  = 2D+00
   DOUBLE PRECISION, PARAMETER   :: Half = 5D-01
   !
   !  G_pq = <Psi_T|ap+.aq|Psi_W>/<Psi_T|Psi_W> = [CW.(CT+.CW)-1.CT+]_qp
   !
   !           |Psi_T> = U|Psi>  (U=CT)
   !           |Psi_W> = V|Psi>  (V=CW)
   !
   !  Theta_pq = [CW.(CT+.CW)-1]_pq = [V.(U+.V)-1]_pq
   !
   !     であるとすれば
   !
   !        <Psi_T|a_p+.a_q|Psi_W>/<Psi_T|Psi_W> h_pq
   !        = <Psi|U+.a_p+.a_q.V|Psi> / <Psi|U+.V|Psi> h_pq
   !        = <Psi|[U+.a_p+.U].U+.V.[V+.a_q.V]|Psi> / <Psi|U+.V|Psi> h_pq
   !        = U*_pi V_qj <Psi|a_i+.[U+.V].a_j|Psi> / <Psi|U+.V|Psi> h_pq
   !        = U*_pi V_qj [(U+.V)-1]_ji h_pq 
   !        = h_iq. Theta_qi 
   !
   !     #
   !     #  [Green_half]_pq = [U.S-1.V+]_pq = [U.[U+.V]^{-1}.V+]_qp
   !     #   Theta_pq       = [V.(U+.V)-1]_pq = GH_qp
   !
   !        <Psi_T|a_p+.a_q+.a_s.a_r|Psi_W>/<Psi_T|Psi_W> (pr|qs)
   !        = <Psi|U+.a_p+.a_q+.a_s.a_r.V|Psi> / <Psi|U+.V|Psi> (pr|qs)
   !        = <Psi|[U+.a_p+.U].[U+.a_q+.U].U+.V.[V+.a_s.V].[V+.a_r.V]|Psi> / <Psi|U+.V|Psi> (pr|qs)
   !        = U*_pi U*_qj V_sl V_rk <Psi|a_i+.a_j+.[U+.V].a_l.a_k|Psi> / <Psi|U+.V|Psi> (pr|qs)
   !        = U*_pi U*_qj V_sl V_rk ( [(U+.V)-1]_ki  [(U+.V)-1]_lj - [(U+.V)-1]_kj [(U+.V)-1]_li) (pr|qs) 
   !        = (ir|js) ( [V.(U+.V)-1]_ri [V.(U+.V)-1]_sj - [V.(U+.V)-1]_rj [V.(U+.V)-1]_si )
   !        = (ir|js) ( θ_ri θ_sj - θ_rj θ_si )
   !        = (ir|js) ( GH_ir GH_js - GH_jr GH_is )
   !
   E0 = ENu
   E1 = Zero
   E2 = Zero
   EJ = Zero
   EK = Zero
   IF (HalfRot) THEN
      DO ispin=1,NSpin
         ij = 0
         DO i=1,NVar
            DO j=1,i-1
               ij = ij +1
               E1 = E1 + Hc(ij)*GH(j,i,ISpin)
               E1 = E1 + Hc(ij)*GH(i,j,ISpin)
            ENDDO
            ij = ij +1
            ii = ij
            E1 = E1 + Hc(ii)*GH(i,i,ISpin)
         ENDDO
      ENDDO
      !
      !  (ij|kl).GH_ij.GH_kl
      !
      DO m=1,NChol
         X= Zero
         DO ispin=1,NSpin
            DO i=1,NVar
               DO j=1,NVar
                  X = X + Chol(i,j,m)*GH(i,j,ispin)
               ENDDO
            ENDDO
         ENDDO
         EJ = EJ + Half*X**2
      ENDDO
      !
      ! (ij|kl) GH_kj GH_il  
      ! j-spin = k-spin <= GH_jk
      ! i-spin = l-spin <= GH_il
      !
      DO ispin=1,NSpin
         DO m=1,NChol
            Y = Zero
            DO i=1,NVar
               DO j=1,NVar
                  DO k=1,NVar
                     Y(i,k) = Y(i,k) + Chol(i,j,m)*GH(k,j,ispin)
                  ENDDO
               ENDDO
            ENDDO
            Z = Zero
            DO i=1,NVar
               DO k=1,NVar
                  DO l=1,NVar
                     Z(k,i) = Z(k,i) + Chol(k,l,m)*GH(i,l,ispin)
                  ENDDO
               ENDDO
            ENDDO
            DO i=1,NVar
               DO k=1,NVar
                  EK = EK + Y(i,k)*Z(k,i)*Half
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF
   E2 = EJ - EK
   E  = E0+E1+E2

   IF (main_rank .AND. PrintT%Energy > 2) THEN 
      WRITE(IW,'(3X,"AFQMC_Energy_Estimator_D")')
      WRITE(IW,'(3X,"EJ    =",1F20.12)')EJ
      WRITE(IW,'(3X,"EK    =",1F20.12)')EK
      WRITE(IW,'(3X,"E0    =",1F20.12)')E0
      WRITE(IW,'(3X,"E1    =",1F20.12)')E1
      WRITE(IW,'(3X,"E0+E1 =",1F20.12)')E0+E1
      WRITE(IW,'(3X,"E2    =",1F20.12)')E2
      WRITE(IW,'(3X,"E_Loc =",1F20.12)')E
   ENDIF
END SUBROUTINE 
