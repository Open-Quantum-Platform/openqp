SUBROUTINE Survial_Pair_Branch(NVar,NSpin,NWlk,Walker,Weight_Min,Weight_Max)
   USE AFQMC_Random_Normal
! OpenQP port: the upstream #ifdef AFQMC_STANDALONE selected between Gellan's
! real MPI layer and a serial stub set. liboqp always uses the serial stubs in
! afqmc_host.F90 (it parallelises with OpenMP), so the guard is resolved here.
   USE MPI_Parallel, ONLY: MPI_Comm_World,MPI_SUM,my_rank,NProc,main_rank,&
                          Gel_AllReduce,Gel_Barrier,Gel_BCast,Gel_Sendrecv
   USE AFQMC_Module,       ONLY: AFQMC_Walker_Type, PrintT
   USE IO_Files,           ONLY: IW
   IMPLICIT NONE
   INTEGER,                 INTENT(IN)    :: NVar
   INTEGER,                 INTENT(IN)    :: NSpin
   INTEGER,                 INTENT(IN)    :: NWlk
   TYPE(AFQMC_Walker_Type), INTENT(INOUT) :: Walker(1:NWlk)
   DOUBLE PRECISION,        INTENT(IN)    :: Weight_Min
   DOUBLE PRECISION,        INTENT(IN)    :: Weight_Max
   DOUBLE PRECISION,        ALLOCATABLE   :: W(:)
   DOUBLE PRECISION,        ALLOCATABLE   :: Wrk(:)
   INTEGER,                 ALLOCATABLE   :: ID(:)
   INTEGER,                 ALLOCATABLE   :: IWrk(:)
   INTEGER,                 ALLOCATABLE   :: List(:,:)
   DOUBLE PRECISION,        ALLOCATABLE   :: WPair(:)
   DOUBLE PRECISION                       :: Ratio
   INTEGER                                :: I, J, K
   INTEGER                                :: NList
   INTEGER                                :: IWlk
   INTEGER                                :: ISource, IDest
   INTEGER                                :: ISender, IRecver
   DOUBLE PRECISION                       :: R
   DOUBLE PRECISION                       :: Dummy
   TYPE(AFQMC_Walker_Type), ALLOCATABLE   :: SendBuf
   TYPE(AFQMC_Walker_Type), ALLOCATABLE   :: RecvBuf
   !
   !
   ! 
   ALLOCATE(W(1:NWlk))
   ALLOCATE(ID(1:NWlk))
   ALLOCATE(Wrk(1:NWlk))
   ALLOCATE(IWrk(1:NWlk))
   ALLOCATE(List(1:2,1:NWlk))
   ALLOCATE(WPair(1:NWlk))
   ID(1:NWlk) = [(I,I=1,NWlk)]
   W = 0.D0
   DO IWlk=my_rank+1,NWlk, NProc
      W(IWlk) = ABS(Walker(IWlk)%Weight)
   ENDDO
   CALL Gel_AllReduce(W,Dummy,NWlk,'D',MPI_SUM,MPI_Comm_World,1)
   CALL Merge_Sort_D(NWlk, ID, IWrk, W, Wrk)
   I = 0
   J = NWlk+1
   IF (PrintT%PopCntrl > 1 .AND. main_rank) THEN
      WRITE(IW,'(4X,"------------------")')
      WRITE(IW,'(4X,"Population Control")')
      WRITE(IW,'(4X,"------------------")')
      WRITE(IW,'(4X,"Weight_Min=",F20.12)')Weight_Min
      WRITE(IW,'(4X,"Weight_Max=",F20.12)')Weight_Max
   ENDIF
   NList = 0
   IF (main_rank) THEN
      K = 0
      DO 
         I = I + 1
         J = J - 1
         IF (PrintT%PopCntrl > 2 .AND. main_rank) THEN
            WRITE(IW,'("I=",I4," J=",I4," WI=",F20.12," WJ=",F20.12," L=",L4)') &
         &    I,J,W(I),W(J), (W(I) >= Weight_Min .AND. W(J) <= Weight_Max)
         ENDIF
         IF (I >= J) EXIT
         IF (W(I) >= Weight_Min .AND. W(J) <= Weight_Max) CYCLE
         CALL Random_Number(r)
         Ratio  = W(J)/(W(I) + W(J))
         IF (PrintT%PopCntrl > 1 .AND. main_rank) THEN
            WRITE(IW,'(4I4,2ES15.6,2F8.3,L)')I,J,ID(I),ID(J),W(I),W(J),r,Ratio,(r < Ratio)
         ENDIF
         K = K + 1
         WPair(K) = 0.5D0*(W(I) + W(J))
         IF (r < Ratio) THEN
            List(1,K) =  J   ! Sender
            List(2,K) =  I   ! Receiver
         ELSE
            List(1,K) =  I   ! Sender
            List(2,K) =  J   ! Receiver
         ENDIF
      ENDDO
      NList = K
      IF (PrintT%PopCntrl > 2 .AND. main_rank) THEN
         DO K=1,NList
            WRITE(IW,'(2X,"K=",I4,2X,"List=",2I4)')K,List(1:2,K)
         ENDDO
      ENDIF
   ENDIF
   ! The pairing list is constructed by the global rank zero, so it must be
   ! distributed through the global communicator (not only within each node).
   CALL Gel_Barrier(MPI_Comm_World)
   CALL Gel_BCast(NList,1,'I',0,MPI_Comm_World)
   IF (NList > 0) THEN
      CALL Gel_BCast(List,2*NList,'I',0,MPI_Comm_World)
      CALL Gel_BCast(WPair,NList,'D',0,MPI_Comm_World)
   ENDIF
   DO K=1,NList
      !
      !--- Clone large weight walker (to smal) ---
      !
      I = List(1,K)   ! Sender
      J = List(2,K)   ! Receiver
      ALLOCATE(SendBuf,SOURCE=Walker(ID(I)))
      ALLOCATE(RecvBuf,SOURCE=Walker(ID(J)))
      ISender   = MOD(ID(I)-1,NProc)  !Sender
      IRecver   = MOD(ID(J)-1,NProc)  !Receiver
      !WRITE(IW,'(2X,"A ID(I)",I4,2X,"ID(J)",I4)')ID(I),ID(J)
      !WRITE(IW,'(2X,"Sender",I4," Recver=",I4," MyRank=",I4)') &
      ! &    ISender,IRecver,my_rank
      !Walker(ID(I)) = Walker(ID(J))
      IF (my_rank == ISender .OR. my_rank == IRecver) THEN
         IF (my_rank == IRecver) THEN !もし自分が受信ランクなら送信(IDest)は関係ない
            IF (PrintT%PopCntrl > 2) THEN
               WRITE(IW,'(2X,"K=",I4," Recver=",I4," JWalk=",I4," WJ=",F20.12," Wave=",F20.12)') &
        &         K,IRecver,ID(J),Walker(ID(J))%Weight,WPair(K)
            ENDIF
            IDest   = ISender !Any rank is okay
            ISource = ISender
         ELSE !もし自分が送信ランクなら受信(ISource)は関係ない
            IF (PrintT%PopCntrl > 2) THEN
               WRITE(IW,'(2X,"K=",I4," Sender",I4," IWalk=",I4," WI=",F20.12," Wave=",F20.12)') &
        &         K,ISender,ID(I),Walker(ID(I))%Weight,WPair(K)
            ENDIF
            IDest   = IRecver
            ISource = IRecver !Any rank is okay
         ENDIF
         CALL Gel_Sendrecv(SendBuf%UW,NVar**2*NSpin,"Z",IDest,0, &
             &             RecvBuf%UW,NVar**2*NSpin,ISource,0,MPI_Comm_World)
         CALL Gel_Sendrecv(SendBuf%E_Hyb,1,"Z",IDest,0, &
             &             RecvBuf%E_Hyb,1,ISource,0,MPI_Comm_World)
         CALL Gel_Sendrecv(SendBuf%E_Loc,1,"Z",IDest,0, &
             &             RecvBuf%E_Loc,1,ISource,0,MPI_Comm_World)
         CALL Gel_Sendrecv(SendBuf%Ovlp,1,"Z",IDest,0, &
             &             RecvBuf%Ovlp,1,ISource,0,MPI_Comm_World)
         CALL Gel_Sendrecv(SendBuf%Det,1,"D",IDest,0,   &
             &             RecvBuf%Det,1,ISource,0,MPI_Comm_World)
      ENDIF
      IF (my_rank == ISender .OR. my_rank == IRecver) RecvBuf%Weight = WPair(K)
      ! もし自分が受信ランクならWalkerをUpdate
      IF (my_rank == IRecver ) Walker(ID(J))        = RecvBuf
      IF (my_rank == ISender ) Walker(ID(I))%Weight = WPair(K)
      DEALLOCATE(SendBuf)
      DEALLOCATE(RecvBuf)
   ENDDO
   DEALLOCATE(W)
   DEALLOCATE(ID)
   DEALLOCATE(Wrk)
   DEALLOCATE(IWrk)
   DEALLOCATE(List)
   DEALLOCATE(WPair)
END SUBROUTINE
