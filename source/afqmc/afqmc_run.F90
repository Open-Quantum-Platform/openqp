SUBROUTINE AFQMC_Init(AFQMC_Cntrl)
   USE AFQMC_Module, ONLY: AFQMC_Control_Type
   IMPLICIT NONE
   TYPE(AFQMC_Control_Type), INTENT(INOUT) :: AFQMC_Cntrl
   AFQMC_Cntrl%NStep          = 1000
   AFQMC_Cntrl%NStep_Stblz    = 5
   AFQMC_Cntrl%NStep_PopCntrl = 5
   AFQMC_Cntrl%NStep_Estimate = 25
   AFQMC_Cntrl%NStep_Accum    = 10000
   AFQMC_Cntrl%NWlk           = 640
   AFQMC_Cntrl%dT             = 5D-03
   AFQMC_Cntrl%Weight_Min     = 0.1D0
   AFQMC_Cntrl%Weight_Max     = 4.0D0
   AFQMC_Cntrl%E_Bound        = SQRT(2.D0/AFQMC_Cntrl%dT)
   AFQMC_Cntrl%x_bound        = 1.D0
   AFQMC_Cntrl%Chol_Tol       = 1D-06
   AFQMC_Cntrl%ISeed          = 1
END SUBROUTINE

! AFQMC_Set (the Gellan SetInp/ChkDlt keyword reader) is not vendored:
! OpenQP fills AFQMC_Control_Type from the [afqmc] input section in
! afqmc_oqp.F90.

SUBROUTINE AFQMC_Run(NVar,NSpin,NOcc,NChol,ENu,Hc,Chol,AFQMC_Cntrl,&
&                    MRSF_CIS,MRSF_File,OO_Orb,OO_File,NLow,Low_File,Low_Max,&
&                    Low_Cap,Low_Start,Results,Ext_Coeff,Ext_Occ)
! OpenQP port, two optional arguments:
!   Results(4)         returns the final [E_Shift, E_Ave, E_Shift_Accum,
!                      E_Ave_Accum] -- the mixed and hybrid estimators and their
!                      post-equilibration running averages -- so the library
!                      entry point can report an energy instead of only printing
!                      one.
!   Ext_Coeff/Ext_Occ  supply the trial determinant expansion in memory, which
!                      is how an OpenQP MRSF-TDDFT root is used as the trial.
!                      When present they take precedence over MRSF_File.
  USE IO_Files, ONLY: IW
  USE AFQMC_Module, ONLY: AFQMC_Walker_Type, AFQMC_Control_Type, PrintT

! OpenQP port: the upstream #ifdef AFQMC_STANDALONE selected between Gellan's
! real MPI layer and a serial stub set. liboqp always uses the serial stubs in
! afqmc_host.F90 (it parallelises with OpenMP), so the guard is resolved here.
  USE MPI_Parallel, ONLY: my_rank,NProc,MPI_Comm_World,my_node_comm,&
                         MPI_SUM,main_rank,Gel_AllReduce,Gel_Barrier
  USE AFQMC_Random_Normal
!
!        Walker(IWlk)%UW(1:NVar,1:NVar)  Complex
!        Walker(IWlk)%E_Hyb              Complex
!        Walker(IWlk)%E_loc              Complex
!        Walker(IWlk)%Weight             Real

  IMPLICIT NONE
  INTEGER,                  INTENT(IN) :: NVar
  INTEGER,                  INTENT(IN) :: NSpin
  INTEGER,                  INTENT(IN) :: NChol
  INTEGER,                  INTENT(IN) :: NOcc(1:NSpin)
  DOUBLE PRECISION,         INTENT(IN) :: ENu
  DOUBLE PRECISION,         INTENT(IN) :: Hc(1:NVar*(NVar+1)/2)
  DOUBLE PRECISION,         INTENT(IN) :: Chol(1:NVar,1:NVar,1:NChol)
  TYPE(AFQMC_Control_Type), INTENT(IN) :: AFQMC_Cntrl
  INTEGER,                  INTENT(IN) :: MRSF_CIS
  CHARACTER(LEN=8),         INTENT(IN) :: MRSF_File
  INTEGER,                  INTENT(IN) :: OO_Orb, NLow
  CHARACTER(LEN=8),         INTENT(IN) :: OO_File, Low_File
   DOUBLE PRECISION,         INTENT(IN) :: Low_Max
   DOUBLE PRECISION,         INTENT(IN) :: Low_Cap
  INTEGER,                  INTENT(IN) :: Low_Start
  DOUBLE PRECISION, OPTIONAL, INTENT(OUT) :: Results(1:4)
  DOUBLE COMPLEX,   OPTIONAL, INTENT(IN)  :: Ext_Coeff(:)
  INTEGER,          OPTIONAL, INTENT(IN)  :: Ext_Occ(:,:,:)

  DOUBLE PRECISION :: H1(1:NVar,1:NVar,1:NSpin)
  DOUBLE PRECISION :: ExpH1(1:NVar,1:NVar,1:NSpin)
  DOUBLE COMPLEX   :: GF(1:NVar,1:NVar,1:NSpin)
  DOUBLE COMPLEX   :: GH(1:NVar,1:NVar,1:NSpin)
  DOUBLE COMPLEX   :: Ovlp
  DOUBLE COMPLEX   :: Ovlp_FB
  DOUBLE PRECISION :: dtheta
  DOUBLE PRECISION :: dt

  INTEGER :: IStep 
  INTEGER :: IWlk
  INTEGER :: IField
  INTEGER :: NWlk
  INTEGER :: NField
  INTEGER :: NDet
  INTEGER :: IDom
  INTEGER :: NStep 
  INTEGER :: NStep_Stblz
  INTEGER :: NStep_PopCntrl
  INTEGER :: NStep_Estimate
  INTEGER :: NStep_Accum

  DOUBLE COMPLEX   :: E_Hyb
  DOUBLE COMPLEX   :: E_Ave

  DOUBLE PRECISION, ALLOCATABLE :: x(:,:)
  DOUBLE COMPLEX                :: X_Shift(1:NChol)
  DOUBLE COMPLEX                :: X_Bar(1:NChol)
  DOUBLE COMPLEX,   ALLOCATABLE :: MF_Shift(:)
  DOUBLE COMPLEX,   ALLOCATABLE :: UW_Tmp(:,:,:)
  DOUBLE COMPLEX,   ALLOCATABLE :: Trial_Coeff(:)
  DOUBLE COMPLEX,   ALLOCATABLE :: Trial_Phi(:,:,:,:)
  INTEGER,          ALLOCATABLE :: Trial_Occ(:,:,:)
  DOUBLE COMPLEX,   ALLOCATABLE :: OO_Coeff(:,:,:)
  DOUBLE COMPLEX,   ALLOCATABLE :: Lower_Coeff(:,:)
  

  TYPE(AFQMC_Walker_Type)              :: Trial
  TYPE(AFQMC_Walker_Type), ALLOCATABLE :: Walker(:)
  DOUBLE PRECISION          :: x_bound
  DOUBLE PRECISION          :: E_bound
  DOUBLE PRECISION          :: W_bound
  DOUBLE PRECISION          :: Weight_Min
  DOUBLE PRECISION          :: Weight_Max
  DOUBLE PRECISION          :: Weight_Scale
  DOUBLE PRECISION          :: Total_Weight
  DOUBLE PRECISION          :: Target_Weight
  DOUBLE COMPLEX, PARAMETER :: ZHalf=(5D-01,0D+00)
  DOUBLE COMPLEX, PARAMETER :: ZZero=(0D+00,0D+00)
  DOUBLE PRECISION, PARAMETER :: DZero=0D+00
  DOUBLE COMPLEX            :: cmf1,cmf2,cmf,cfb,v2
  !DOUBLE COMPLEX, EXTERNAL  :: ZDOTU
  DOUBLE COMPLEX            :: omega
  DOUBLE COMPLEX            :: Important_Function
  DOUBLE COMPLEX            :: E_Shift
  DOUBLE PRECISION          :: E_Max, E_Min
  DOUBLE PRECISION          :: E_Real, E_Imag
  DOUBLE PRECISION          :: Weight0
  DOUBLE PRECISION          :: Dummy
  DOUBLE COMPLEX            :: ZDummy
  DOUBLE COMPLEX            :: E_Shift_Accum
  DOUBLE COMPLEX            :: E_Ave_Accum
  INTEGER                   :: IStep_Accum
  INTEGER                   :: ISeed
  INTEGER                   :: I, ISpin

  NField         = NChol
  dT             = AFQMC_Cntrl%dT
  NWlk           = AFQMC_Cntrl%NWlk
  ISeed          = AFQMC_Cntrl%ISeed
  !
  !--- Interval parameters ---
  !
  NStep          = AFQMC_Cntrl%NStep
  NStep_Stblz    = AFQMC_Cntrl%NStep_Stblz 
  NStep_PopCntrl = AFQMC_Cntrl%NStep_PopCntrl
  NStep_Estimate = AFQMC_Cntrl%NStep_Estimate
  NStep_Accum    = AFQMC_Cntrl%NStep_Accum
  !
  !--- Clipping parameters ---
  !
  Weight_Min = AFQMC_Cntrl%Weight_Min
  Weight_Max = AFQMC_Cntrl%Weight_Max
  E_Bound    = AFQMC_Cntrl%E_Bound 
  x_bound    = AFQMC_Cntrl%X_Bound
  !
  !
  !
  Total_Weight  = DBLE(NWlk)
  Target_Weight = DBLE(NWlk)
  !
  !
  PrintT%QMCLog = 1
  PrintT%PopCntrl    = 0
  !PrintT%QMCLog = 2
#if 0
  PrintT%Walker      = 5
  PrintT%RandNum     = 5
  PrintT%Energy      = 5
  PrintT%Hamiltonian = 5
  PrintT%Misc        = 5
  PrintT%PopCntrl    = 5
  PrintT%QMCLog      = 5
#endif

  IF (PRESENT(Ext_Coeff)) THEN
     ! OpenQP port: the determinant expansion was handed in directly (e.g. built
     ! from an MRSF-TDDFT root in memory) instead of read from a trial file.
     NDet = SIZE(Ext_Coeff)
  ELSE IF (MRSF_CIS /= 0) THEN
     CALL AFQMC_Read_MRSFCIS_Header(MRSF_File,NDet,NOcc,NSpin)
  ELSE
     NDet = 1
  ENDIF
  ALLOCATE(Trial_Coeff(1:NDet))
  ALLOCATE(Trial_Phi(1:NVar,1:NVar,1:NSpin,1:NDet))
  ALLOCATE(Trial_Occ(1:NVar,1:NSpin,1:NDet))
  IF (PRESENT(Ext_Coeff)) THEN
     Trial_Coeff(1:NDet) = Ext_Coeff(1:NDet)
     Trial_Occ(1:NVar,1:NSpin,1:NDet) = Ext_Occ(1:NVar,1:NSpin,1:NDet)
     CALL AFQMC_Trial_From_Occupations(NVar,NSpin,NOcc,NDet,Trial_Coeff,&
    &                                  Trial_Occ,Trial_Phi,IDom)
  ELSE IF (MRSF_CIS /= 0) THEN
     CALL AFQMC_Read_MRSFCIS_Trial(MRSF_File,NVar,NSpin,NOcc,NDet,&
    &                           Trial_Coeff,Trial_Phi,Trial_Occ,IDom)
  ELSE
     Trial_Coeff(1) = (1.D0,0.D0)
     CALL AFQMC_Init_Wfn_Z(1,NVar,NSpin,NOcc,Trial_Phi)
     Trial_Occ = 0
     DO ISpin=1,NSpin
        DO I=1,NOcc(ISpin)
           Trial_Occ(I,ISpin,1)=I
        ENDDO
     ENDDO
     IDom = 1
  ENDIF

  IF (OO_Orb /= 0) THEN
     ALLOCATE(OO_Coeff(1:NVar,1:NVar,1:NSpin))
     CALL AFQMC_Read_OO_Orbitals_Z(OO_File,NVar,NSpin,OO_Coeff)
     CALL AFQMC_Apply_Orbital_Transform_Z(NVar,NSpin,NDet,OO_Coeff,Trial_Phi)
  ENDIF

  IF (NLow > 0) THEN
     IF (MRSF_CIS == 0) ERROR STOP "NLOW requires an MRSF-CIS determinant basis"
     ALLOCATE(Lower_Coeff(1:NDet,1:NLow))
     CALL AFQMC_Read_Lower_States_Z(Low_File,NDet,NLow,Lower_Coeff)
     CALL AFQMC_Project_Lower_States_Z(NDet,NLow,Lower_Coeff,Trial_Coeff)
     CALL AFQMC_Dominant_Coefficient_Z(NDet,Trial_Coeff,IDom)
     IF (main_rank) THEN
        WRITE(IW,'(1X,"Lower-state collapse threshold/cap/start = ",2ES12.4,1X,I0)')&
       &     Low_Max,Low_Cap,Low_Start
     ENDIF
  ENDIF

  IF (MRSF_CIS /= 0) THEN
     CALL AFQMC_Report_Trial_Spin_Z(NVar,NSpin,NOcc,NDet,Trial_Coeff,&
    &                               Trial_Occ,MRSF_CIS)
  ENDIF

  ALLOCATE(Walker(1:NWlk))
  ALLOCATE(Trial%UW(1:NVar,1:NVar,1:NSpin))
  DO IWlk=1,NWlk
     ALLOCATE(Walker(IWlk)%UW(1:NVar,1:NVar,1:NSpin))
     Walker(IWlk)%E_Hyb  = ZZero
     Walker(IWlk)%E_Loc  = ZZero
     Walker(IWlk)%Weight = DZero
     Walker(IWlk)%Det    = DZero
     Walker(IWlk)%Ovlp   = ZZero
  ENDDO
  ALLOCATE(UW_Tmp(1:NVar,1:NVar,1:NSpin))
  ALLOCATE(x(1:NField,1:NWlk))
  ALLOCATE(mf_shift(1:NField))
  !
  !------------------------------------
  !   1. Set up the tiral wavefunction 
  !------------------------------------
  !
  !--- 1-1. Initiate the tiral wavefunction ---
  !
  Trial%UW = Trial_Phi(:,:,:,IDom)
  !
  !--- 1-2. Obtain density matrix for trial to generate mean-field (MF) substruction  ---
  !
  CALL AFQMC_Greens_Functions_Z(NVar,NSpin,NOcc,Trial%UW,Trial%UW,Ovlp,GF,GH)
  !
  !--- 1-3. Evaluate the dominant-determinant energy used to start E_Shift ---
  !
  CALL AFQMC_Energy_Estimator_Z(NVar,NSpin,NChol, &
     &                          GF,GH,ENu,Hc,Chol,E_Shift,.TRUE.)      
  !
  !--- 1-4. Construct the mixed-reference density and mean-field shift ---
  !
  IF (MRSF_CIS /= 0) THEN
     CALL AFQMC_TrialEnsemble_Green_Z(NVar,NSpin,NOcc,NDet,&
    &                                 Trial_Coeff,Trial_Phi,GF)
  ENDIF
  CALL AFQMC_Mean_Field_Shift_Z(NVar,NSpin,NChol,Chol,GF,MF_Shift)
  IF (PrintT%Misc > 1) THEN
     WRITE(IW,'("Initial energy=",2F20.12)')E_Shift
  ENDIF
  !
  !--------------------------------------
  !   2. Set up the walker wavefunctions
  !--------------------------------------
  !
  DO IWlk = my_rank+1, NWlk, NProc
     Walker(IWlk)%UW = Trial%UW
     Walker(IWlk)%Weight = 1D+00
     Walker(IWlk)%E_Hyb  = E_Shift
     CALL AFQMC_MultiDet_Overlap_Z(NVar,NSpin,NOcc,NDet,Trial_Coeff,&
    &                              Trial_Phi,Walker(IWlk)%UW,Walker(IWlk)%Ovlp)
  ENDDO
  !
  !-------------------------
  !   3. Set up Hamiltonians
  !-------------------------
  !
  !
  !--- 3-1. Calculate one-body H1 ---
  !
  !      v0_{pq} = H1_{pq} ap+.aq 
  !      H1_{pq} = h_{pq} - 1/2 L_{pr,m} L_{qr,m} 
  !              = h_{pq} - 1/2 (pr|qr) = h_{pq} - 1/2 <pr|rq>
  !
  !      h_{pq}: Core Hamiltonioan
  !
  CALL AFQMC_Construct_One_Body_Hamiltonian(NVar,NSpin,NChol,Hc,Chol,H1)
  !
  !--- 3-2. Construct one body propagator ---
  !
  CALL Construct_One_Body_Propagator(NVar,NSpin,NChol,H1,Chol,MF_Shift,dt,ExpH1)
  !
  !
  !------------------------------
  !   4. Initial AFQMC parameters
  !-------------------------------
  CALL AFQMC_Init_random_seed(ISeed + my_rank)
  !
  !--- Run AFQMC ---
  !
  E_Shift_Accum = ZZero
  E_Ave_Accum = ZZero
  IStep_Accum = 0
  ! OpenQP port: label the columns of the per-step table below.
  IF (main_rank) THEN
     WRITE(IW,'(/,1X,A)')"AFQMC propagation"
     WRITE(IW,'(1X,A6,3A20,A1,2A20)')"step","E(mixed)","E(hybrid)",&
    &     "total weight","|","E(mixed) avg","E(hybrid) avg"
  ENDIF
  DO IStep = 1, NStep
     IF (IStep > NStep_Accum) THEN
        IStep_Accum = IStep_Accum + 1
     ENDIF
     IF (PrintT%QMCLog > 1) THEN
        WRITE(IW,'(/,3X,10("="))')
        WRITE(IW,'(3X,"IStep=",I4)')IStep
        WRITE(IW,'(3X,10("="),/)')
     ENDIF
     IF (MOD(IStep,NStep_Stblz) == 0) THEN
        ! OpenQP port: the re-orthogonalisation is independent per walker, so it
        ! is threaded like the propagation itself. The scratch copy moved into
        ! AFQMC_Stabilize_Walker_Z as an automatic array, because the shared
        ! UW_Tmp used here would race. Each walker's QR is untouched, so the
        ! result does not depend on the thread count.
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(IWlk)
        DO IWlk = my_rank+1, NWlk, NProc
           CALL AFQMC_Stabilize_Walker_Z(NVar,NSpin,NOcc,Walker(IWlk)%UW,&
          &                              Walker(IWlk)%Det,Walker(IWlk)%Ovlp)
        ENDDO
!$OMP END PARALLEL DO
        CALL Gel_Barrier(my_node_comm)
     ENDIF

     ! Generate random fields outside the threaded region.  This avoids races
     ! in the RNG state and makes a run reproducible when OMP_NUM_THREADS is
     ! changed (for a fixed MPI layout and seed).
     DO IWlk = my_rank+1, NWlk, NProc
        ! OpenQP port: one block draw instead of NField scalar calls. The
        ! uniform stream is consumed in the same order, so the fields are
        ! bit-identical to the scalar path.
        CALL AFQMC_Rand_Normal_Array(0.0D0, 1.0D0, NField, x(1,IWlk))
        IF (PrintT%RandNum > 1 .AND. main_rank) THEN
           WRITE(IW,'(/,4X,32("-"))')
           WRITE(IW,'(4X,"IField",2X,"IWlk",2X,"Aux. Field Coord.")')
           WRITE(IW,'(4X,32("-"))')
           DO IField=1,NField
              WRITE(IW,'(4X,2I6,F20.12)') IField,IWlk,x(IField,IWlk)
           ENDDO
           WRITE(IW,'(4X,32("-"),/)')
        ENDIF
     ENDDO

!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) &
!$OMP& IF(PrintT%QMCLog <= 1 .AND. PrintT%Walker <= 1 .AND. PrintT%RandNum <= 1) &
!$OMP& PRIVATE(IWlk,IField,GF,GH,Ovlp,Ovlp_FB,X_Shift,X_Bar,E_Hyb,dTheta) &
!$OMP& PRIVATE(cmf1,cmf2,cmf,cfb,v2,omega,Important_Function) &
!$OMP& PRIVATE(E_Min,E_Max,E_Real,E_Imag,Weight0)
     DO IWlk = my_rank+1, NWlk, NProc
        !
        !--- 1. Propagate one independently sampled walker ---
        !
        IF (PrintT%QMCLog > 1 .AND. main_rank) THEN
           WRITE(IW,'(/,5X,"--------------------------")')
           WRITE(IW,'(5X," 2. Phaseless propagation ")')
           WRITE(IW,'(5X,"--------------------------")')
           WRITE(IW,'(5X,"IWlk=",I4,/)')IWlk
        ENDIF
        !
        !--- 2-1. Calculate Green's functions ---
        !
        IF (MRSF_CIS /= 0) THEN
           CALL AFQMC_MultiDet_Greens_Z(NVar,NSpin,NOcc,NDet,Trial_Coeff,&
          &                             Trial_Phi,Walker(IWlk)%UW,Ovlp_FB,GF)
        ELSE
           CALL AFQMC_Greens_Functions_Z(NVar,NSpin,NOcc,&
          &                              Trial%UW,Walker(IWlk)%UW,Ovlp_FB,GF,GH)
        ENDIF
        !
        !--- 2-2. Calculate force bias ---
        !
        ! The MRSF force bias and phaseless overlap both use the complete
        ! determinant expansion.  The mean-field subtraction uses a stable
        ! coefficient-weighted diagonal ensemble of the imported trial.
        CALL AFQMC_Force_Bias_Z(NVar,NSpin,NChol,dT,GF, &
        &                       Chol,MF_Shift,X_Bar)
        DO IField=1,NField
           IF (ABS(x_bar(IField)) > x_bound) THEN
              IF (ABS(x_bar(IField)) > 1D-13) THEN
                 x_bar(IField) = x_bound*x_bar(IField)/ABS(x_bar(IField))
              ENDIF
           ENDIF
        ENDDO
        x_shift(1:NField) = x(1:NField,IWlk) - x_bar(1:NField)       !x_shift: Complex
        IF (PrintT%QMCLog > 2 .AND. main_rank) THEN
           WRITE(IW,'(5X,"--------------------")')
           WRITE(IW,'(5X," IField  X Shift    ")')
           WRITE(IW,'(5X,"--------------------")')
           DO IField=1,NField
               WRITE(IW,'(5X,I6,2F20.12)')IField,x_shift(IField)
           ENDDO
        ENDIF
        !
        !--- 2-3. Update wavefunctions ---
        ! 
        CALL AFQMC_Propagate_One_Body_Z(NVar,NSpin,ExpH1,Walker(IWlk)%UW)
        CALL AFQMC_Propagate_Two_Body_Z(NVar,NSpin,NChol,dT,x_shift,Chol,Walker(IWlk)%UW)
        CALL AFQMC_Propagate_One_Body_Z(NVar,NSpin,ExpH1,Walker(IWlk)%UW)
        !
        !   ===================================================
        !   2-4--2-13: ipie.propagation.phaseless_baseupdate_weight
        !   ===================================================
        !
        !--- 2-4. Calculate constant factors arising from force bias and mean field shift ---
        !
        !    cmf = -self.sqrt_dt * xp.einsum("wx,x->w", xshifted, self.mf_shift)
        !    cmf[w] = -√Δτ \sum_m (x_wm - xbar_wm) . mf_m
        !
        ! ZDOTU did not work ...
        !
        !cmf1 = - CMPLX(SQRT(dt),0D+00) * ZDOTU(NField,x_shift,1,mf_shift,1)         ! cm1: Complex
        !v2   = ZDOTU(NField,mf_shift,1,mf_shift,1)                     ! v2 : Negative Double
        !
        ! Constant factor arising from shifting the propability distribution.
        !    cfb = xp.einsum("wx,wx->w", xi, xbar) - 0.5 * xp.einsum("wx,wx->w", xbar, xbar)
        !
        !    cfb[w] = \sum_m [ x_{wm} xbar_{wm} - 1/2 xbar_{wm}^2 ]
        !
        !cfb = ZDOTU(NField,x2,1,x_bar,1) - Half*ZDOTU(NField,x_bar,1,x_bar,1)  ! Complex
        cmf1 = ZZero
        cfb  = ZZero 
        v2   = ZZero 
        DO IField=1,NField
           cmf1 = cmf1 - DSQRT(dt) * x_shift(IField)*mf_shift(IField)
           v2   = v2 + mf_shift(IField)*mf_shift(IField)
           cfb  = cfb + x(IField,IWlk)*x_bar(IField) - ZHalf*x_bar(IField)*x_bar(IField)
        ENDDO
        cmf2 = (- CMPLX(ENu,0D+0,KIND=KIND(cmf2)) - ZHalf*v2) * &
       &       CMPLX(dT,0D+0,KIND=KIND(cmf2))
        cmf  = cmf1 + cmf2                                           
        !
        !--- 2-5. Calculate overlap ratio ---
        !
        CALL AFQMC_MultiDet_Overlap_Z(NVar,NSpin,NOcc,NDet,Trial_Coeff,&
       &                              Trial_Phi,Walker(IWlk)%UW,Ovlp)
        omega = Ovlp/Walker(IWlk)%Ovlp
        IF (PrintT%QMCLog > 1 .AND. main_rank) THEN
           WRITE(IW,'(5X,"--------------------------")')
           WRITE(IW,'(5X," 2-5. Overlap ratio ")')
           WRITE(IW,'(5X,"--------------------------")')
           WRITE(IW,'(5X,"Current  Ovlp. =",3F20.12)')Ovlp,ABS(Ovlp)
           WRITE(IW,'(5X,"Previous Ovlp. =",3F20.12)')Walker(IWlk)%Ovlp,ABS(Ovlp)
           WRITE(IW,'(5X,"omega          =",2F20.12)')omega
        ENDIF 
        !
        !--- 2-6. Calculate hybrid energy ---
        !
        E_Hyb = - (LOG(omega) + cfb + cmf)/dT 
        IF (PrintT%QMCLog > 1 .AND. main_rank) THEN
           WRITE(IW,'(5X,"------------------------------")')
           WRITE(IW,'(5X," 2-6. Calculate hybrid energy ")')
           WRITE(IW,'(5X,"------------------------------")')
           WRITE(IW,'(5X,"E_Hyb      =",2F20.12)')E_Hyb
           WRITE(IW,'(5X,"ENu        =", F20.12)')ENu
           WRITE(IW,'(5X,"log(omega) =",2F20.12)')LOG(omega)
           WRITE(IW,'(5X,"v2         =",2F20.12)')v2
           WRITE(IW,'(5X,"cmf1       =",2F20.12)')cmf1
           WRITE(IW,'(5X,"cmf2       =",2F20.12)')cmf2
           WRITE(IW,'(5X,"cmf        =",2F20.12)')cmf
           WRITE(IW,'(5X,"cfb        =",2F20.12)')cfb
        ENDIF

        !
        !--- 2-7. Apply E_hyb ---
        !
        IF (ABS(E_Shift) > 1D-10) THEN
           E_Min  = DBLE(E_Shift) - E_Bound
           E_Max  = DBLE(E_Shift) + E_Bound
           E_Real = MIN(MAX(DBLE(E_Hyb), E_Min),E_Max)
           E_Imag = AIMAG(E_Hyb)
           E_Hyb  = CMPLX(E_Real,E_Imag,KIND=KIND(E_Hyb))
        ENDIF
        !   
        !--- 2-8. Importance Sampling ---
        ! 
        Important_Function = EXP(-dT*(0.5D0*(E_Hyb + Walker(IWlk)%E_Hyb) - E_Shift))
        !
        ! importance_function = xp.exp(
        !    -self.dt * (0.5 * (hybrid_energy + walkers.hybrid_energy) - eshift)
        !)
        !Bug
        !Important_Function = EXP(-0.5D0*dT*(E_Hyb + Walker(IWlk)%E_Hyb - E_Shift))
        IF (PrintT%QMCLog > 1 .AND. main_rank) THEN
           WRITE(IW,'(5X,"--------------------------")')
           WRITE(IW,'(5X," 2-8. Importance Sampling ")')
           WRITE(IW,'(5X,"--------------------------")')
           WRITE(IW,'(5X,"E_Hyb       =",2F20.12)')E_Hyb
           WRITE(IW,'(5X,"Walk.E_Hyb  =",2F20.12)')Walker(IWlk)%E_Hyb
           WRITE(IW,'(5X,"E_Shift     =",2F20.12)')E_Shift
           WRITE(IW,'(5X,"Imp. Func.  =",2F20.12)')Important_Function
        ENDIF
        !
        !--- 2-9. Phase change ---
        ! 
        dTheta = AIMAG(-dT*E_hyb - cfb)
        !dtheta = (-self.dt * hybrid_energy - cfb).imag
        !cosine_fac = xp.cos(dtheta)

        IF (PrintT%QMCLog > 1 .AND. main_rank) THEN
           WRITE(IW,'(5X,"--------------------------")')
           WRITE(IW,'(5X," 2-9. Phase change ")')
           WRITE(IW,'(5X,"--------------------------")')
           WRITE(IW,'(5X,"dTheta       =", F20.12)')dTheta
           WRITE(IW,'(5X,"Cos(dTheta)  =", F20.12)')DCOS(dTheta)
           WRITE(IW,'(5X,"|Imp. Func.| =", F20.12)')ABS(Important_Function)
        ENDIF
        !
        !--- 2-10. Update weight ---
        !
        Weight0 = Walker(IWlk)%Weight
        Walker(IWlk)%Weight = ABS(Important_Function)*MAX(0D+00,DCOS(dTheta))*Weight0
        !
        !walkers.weight = walkers.weight * magn * cosine_fac
        !
        IF (PrintT%QMCLog > 1) THEN
           WRITE(IW,'(5X,"--------------------------")')
           WRITE(IW,'(5X," 2-10. Update weight ")')
           WRITE(IW,'(5X,"--------------------------")')
           WRITE(IW,'(5X, "Weight (old) =", F20.12)') Weight0
           WRITE(IW,'(5X, "Weight (new) =", F20.12)') Walker(IWlk)%Weight
        ENDIF
        !
        !--- 2-11. Update hybrid energy ---
        !
        Walker(IWlk)%E_Hyb = E_Hyb
        !
        !--- 2-12. Update overlap ---
        !
        Walker(IWlk)%Ovlp   = Ovlp
     ENDDO
!$OMP END PARALLEL DO
     !
     !--- 3. Control walker's weight ---
     !
     IF (IStep > 1) THEN
        Total_Weight = 0.D0
        DO IWlk=my_rank+1,NWlk, NProc
           Total_Weight = Total_Weight + ABS(Walker(IWlk)%Weight)
        ENDDO
        CALL Gel_AllReduce(Total_Weight,Dummy,1,'D',MPI_SUM,MPI_Comm_World,1)
        !IF (MOD(IStep,NStep_Stblz) == 0) THEN
        W_Bound = Total_Weight * 0.1D0
        DO IWlk=my_rank+1,NWlk, NProc
           Weight0 = Walker(IWlk)%Weight
           Weight0 = MAX(-W_Bound, Weight0)
           Weight0 = MIN( W_Bound, Weight0)
           Walker(IWlk)%Weight = Weight0
        ENDDO
        !ENDIF
     ENDIF
     !
     !--- 4. Population controll ---
     !
     IF (MOD(IStep,NStep_PopCntrl) == 0) THEN
        Total_Weight = 0.D0
        DO IWlk=my_rank+1, NWlk, NProc
           Total_Weight = Total_Weight + ABS(Walker(IWlk)%Weight)
        ENDDO
        CALL Gel_AllReduce(Total_Weight,Dummy,1,'D',MPI_SUM,MPI_Comm_World,1)
        IF (Total_Weight < 1D-08) THEN
           ERROR STOP "Walkers' weight is significantly small"
        ENDIF
        !
        !--- 4-1. Renormalize the walker's weight ---
        !
        Weight_Scale = Total_Weight/Target_Weight
        DO IWlk=my_rank+1, NWlk, NProc
           Walker(IWlk)%Weight = Walker(IWlk)%Weight/Weight_Scale 
        ENDDO
        !
        !--- 4-2. Death/Cloning step ---
        !
        CALL Survial_Pair_Branch(NVar,NSpin,NWlk,Walker,Weight_Min,Weight_Max)
     ENDIF
     !
     !--- 5. Update accumlators ---
     !
     E_Ave        = ZZero
     Total_Weight = DZero
     DO IWlk=my_rank+1,NWlk,NProc
        E_Ave        = E_Ave + Walker(IWlk)%Weight*Walker(IWlk)%E_Hyb 
        Total_Weight = Total_Weight + Walker(IWlk)%Weight
     ENDDO
     CALL Gel_AllReduce(E_Ave,ZDummy,1,'Z',MPI_SUM,MPI_Comm_World,1)
     CALL Gel_AllReduce(Total_Weight,Dummy,1,'D',MPI_SUM,MPI_Comm_World,1)
     E_Ave   = E_Ave/Total_Weight
     !
     !--- 6. Estimators ---
     !
     IF (MOD(IStep,NStep_Estimate) == 0) THEN
        CALL AFQMC_Energy_Estimate(NVar,NSpin,NOcc,NChol,NWlk,Walker,&
       &                           NDet,Trial_Coeff,Trial_Phi,ENu,Hc,Chol,E_Shift)
        IF (NLow > 0) THEN
           CALL AFQMC_Lower_State_Diagnostic_Z(NVar,NSpin,NOcc,NWlk,Walker,&
          &       NDet,NLow,Lower_Coeff,Trial_Phi,IStep,Low_Max,&
          &       Low_Cap,IStep >= Low_Start)
        ENDIF
     ENDIF

     IF (IStep_Accum > 0) THEN
        E_Ave_Accum   = (E_Ave_Accum*DBLE(IStep_Accum-1) + E_Ave)/DBLE(IStep_Accum)
        E_Shift_Accum = (E_Shift_Accum*DBLE(IStep_Accum-1) + E_Shift)/DBLE(IStep_Accum)
     ENDIF
     !
     !--- 7. Uodate energy estimate ---
     !
     IF (main_rank) THEN
        IF (IStep > NStep_Accum) THEN
           WRITE(IW,'(1X,I6,3F20.12,"|",2F20.12)')IStep,DBLE(E_Shift),DBLE(E_Ave),&
          &     Total_Weight,DBLE(E_Shift_Accum),DBLE(E_Ave_Accum)
        ELSE
           WRITE(IW,'(1X,I6,3F20.12)')IStep,DBLE(E_Shift),DBLE(E_Ave),Total_Weight
        ENDIF
     ENDIF
  ENDDO

  IF (PRESENT(Results)) THEN
     Results(1) = DBLE(E_Shift)
     Results(2) = DBLE(E_Ave)
     Results(3) = DBLE(E_Shift_Accum)
     Results(4) = DBLE(E_Ave_Accum)
  ENDIF

  DEALLOCATE(x)
  DEALLOCATE(mf_shift)
  DEALLOCATE(UW_Tmp)
  DEALLOCATE(Trial_Coeff)
  DEALLOCATE(Trial_Phi)
  IF (ALLOCATED(OO_Coeff)) DEALLOCATE(OO_Coeff)
  IF (ALLOCATED(Lower_Coeff)) DEALLOCATE(Lower_Coeff)
  DO IWlk=1,NWlk
     DEALLOCATE(Walker(IWlk)%UW)
  ENDDO
  DEALLOCATE(Walker)
  DEALLOCATE(Trial%UW)
 
END SUBROUTINE
