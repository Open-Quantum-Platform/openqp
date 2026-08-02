SUBROUTINE AFQMC_Construct_One_Body_Hamiltonian(NVar,NSpin,NChol,Hc,Chol,H1)
!
!   ipie.hamiltonians.utils.get_generic_integrals                     [5]
!      ipie.hamiltonians.generic.construct_h1e_mod                    [6]
!
!   Subtract one-body bit following reordering of 2-body operators.
!   Eqn (17) of [Motta17].
!
!    H = H0 + H1 + H2
!      = H0 + \sum_pq hpq.ap+.aq + 1/2 \sum_pqrs <pq|rs>ap+.aq+.as.ar
!    To use the Habburd Stranovich transformation, two-body part is reordered as one-body product 
!    H = H0 + \sum_pqr (hpq -1/2 <pr|rq>) ap+.aq - 1/2 \sum_pqrs <pq|rs> ap+.ar.aq+.as
!      = H0 + \hat{v}_0 - 1/2 \sum_m \sum_pr [Lpr,m ap+.ar] \sum_qs [Lqs,m aq+.aq]
!   
!    where 
!         \hat{v}_0 = [h_pq - 1/2 \sum_r (pr|qr)] ap+.aq
!                   = [h_pq - 1/2 \sum_{mr} Lpr,m Lqr,m] ap+.aq
!
  USE IO_Files, ONLY: IW
  USE AFQMC_Module, ONLY: PrintT
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: NVar
  INTEGER,          INTENT(IN)  :: NSpin
  INTEGER,          INTENT(IN)  :: NChol
  DOUBLE PRECISION, INTENT(IN)  :: Hc(1:NVar*(NVar+1)/2)
  DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
  DOUBLE PRECISION, INTENT(OUT) :: H1(1:NVar,1:NVar,1:NSpin)
  DOUBLE PRECISION              :: Shift(1:NVar,1:NVar)  
  !                                 = 1/2 <pr|rq> = 1/2 sum_{rm} L_{pr,m} L_{qr,m}
  INTEGER :: ISpin
  CALL DGeMM('N','T',NVar,NVar,NVar*NChol,0.5D0,Chol,NVar,Chol,NVar,0.D0,Shift,NVar)
  IF (PrintT%Hamiltonian > 2) THEN
     WRITE(IW,'("AFQMC_Construct_One_Body_Hamiltonian")')
     WRITE(IW,'(" H1_{pq} = h_{pq}  -  0.5 L_{pr,m} L_{qr,m}",&
               &"=  h_{pq} - 1/2 (pr|qr) = h_{pq} - 1/2 <pr|rq>")')
     CALL PrSq(Shift,NVar,NVar,NVar)
  ENDIF
  DO ISpin=1,NSpin
     CALL Expnd(Hc,H1(1,1,ISpin),NVar,0)
     H1(1:NVar,1:NVar,ISpin) = H1(1:NVar,1:NVar,ISpin) - Shift(1:NVar,1:NVar)
  ENDDO
  IF (PrintT%Hamiltonian > 2) THEN
     WRITE(IW,'(" H1(α): ipie.hamiltonians.generic.construct_h1e_mod.py")')
     CALL PrSq(H1(1,1,1),NVar,NVar,NVar)
  ENDIF
END SUBROUTINE


! self.mf_shift = construct_mean_field_shift(hamiltonian, trial)
! self.expH1 = construct_one_body_propagator(hamiltonian, self.mf_shift, self.dt)

SUBROUTINE Construct_One_Body_Propagator(NVar,NSpin,NChol,H1,Chol,MF_Shift,dt,ExpH1)
!
!    Construct mean-field shifted one-body propagator.
!        
!    .. math::
!        
!        H1 \rightarrow H1 - v0
!        v0_{ik} = \sum_n v_{(ik),n} \bar{v}_n
!           
!    Parameters
!    ----------
!    hamiltonian : hamiltonian class.
!        Generic hamiltonian object.
!    mf_shift : xp.ndarray
!        Average value of Choleskies with respect to the trial wavefunction.
!    dt : float
!        Timestep.
!    """
!
!    H = H0 + H1 + H2
!      = H0 + \sum_pq hpq.ap+.aq + 1/2 \sum_pqrs <pq|rs>ap+.aq+.as.ar
!    Using fermionic anticommutation relations, 
!        aq+.as = Iqs - as.aq+
!    the two‑body term can be rewritten in density–density form:
!    H = H0 + \sum_pqr (hpq -1/2 <pr|rq>) ap+.aq - 1/2 \sum_pqrs <pq|rs> ap+.ar.aq+.as
!      = H0 + \hat{v}_0                          - 1/2 \sum_m \sum_pr [Lpr,m ap+.ar] \sum_qs [Lqs,m aq+.as]
!    where 
!     [\hat{v}_0]_pq = [h_pq - 1/2 \sum_r (pr|qr)] ap+.aq
!                    = [h_pq - 1/2 \sum_{mr} Lpr,m Lqr,m] ap+.aq
!    is defined as the effective one‑body operator.
!    Now, we consider the mean-field substraction.
!    H = H0 + \sum_pqr (hpq -1/2 <pr|rq>) ap+.aq - 1/2 \sum_pqrs <pq|rs> ap+.ar.aq+.as
!      = H0 + \hat{v}_0 - 1/2 \sum_m \sum_pr [Lpr,m ap+.ar] \sum_qs [Lqs,m aq+.aq]
!      = H0 + \hat{v}_0 - 1/2 \sum_m \hat{v}_m^2
!   Mean-field shift value \bar{v}_m is definded as
!       \bar{v}_m = <Psi_T| \hat{v}_m |Psi_T> / <Psi_T|Psi_T>
!   The shift may reduce the fluctuation of the auxualy fields
!     H = H0 + \hat{v}_0  -1/2 \sum_m (\hat{v}_m - \bar{v}_m)^2
!           - \sum_m \hat{v}_m \bar{v}_m + 1/2 \sum_m (\bar{v}_m)^2
!   Accordingly,        
!     H0' = H0 + 1/2 \sum_m (\bar{v}_m)^2
!
!     H1' = \hat{v}_0 - \sum_m \hat{v}_m \bar{v}_m
!     (H1')_pq = [\hat{v}_0]_pq - \sum_m  \bar{v}_m [\hat{v}_m]_{pq}
!
!     H2' = 1/2 \sum_m  (i \hat{v}_m - i \bar{v}_m)^2
!
!   See also  
!        H1 \rightarrow H1 - v0
!        v0_{ik} = \sum_n v_{(ik),n} \bar{v}_n
!           
!   shift = 1j * numpy.einsum("mx,x->m", hamiltonian.chol, mf_shift).reshape(nb, nb)
!
!
! ipie.propagation.phaseless_base.construct_one_body_propagator 
!
!  walkers.phia = propagate_one_body(walkers.phia, self.expH1[0])
!
!
!    self.expH1 = construct_one_body_propagator(hamiltonian, self.mf_shift, self.dt)
!
!
!   expH1 = numpy.array(
!        [scipy.linalg.expm(-0.5 * dt * H1[0]), scipy.linalg.expm(-0.5 * dt * H1[1])]
!          Alpha block                           Beta block
!        )
!   shift = 1j * numpy.einsum("mx,x->m", hamiltonian.chol, mf_shift).reshape(nb, nb)
   USE IO_Files, ONLY: IW
   USE AFQMC_Module, ONLY: PrintT
   IMPLICIT NONE
   INTEGER,          INTENT(IN)  :: NVar
   INTEGER,          INTENT(IN)  :: NSpin
   INTEGER,          INTENT(IN)  :: NChol
   DOUBLE PRECISION, INTENT(IN)  :: H1(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
   DOUBLE COMPLEX,   INTENT(IN)  :: MF_Shift(1:NChol) ! Imaginary quantity
   DOUBLE PRECISION, INTENT(IN)  :: dt
   DOUBLE PRECISION, INTENT(OUT) :: ExpH1(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION              :: Shift1(1:NVar,1:NVar)
   DOUBLE PRECISION              :: Shift2(1:NVar,1:NVar)
   DOUBLE PRECISION              :: H1m(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION              :: A(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION              :: ExpA(1:NVar,1:NVar,1:NSpin)
   INTEGER :: i
   !
   !   H1' = H1 - v_p <v_p>_T +  1/2 <v_p>_T^2
   !     <v_p>_T = [MF_shfit]_x
   !      v_p    = [Chol]_mx
   !
   !   Scalar term  "1/2 <v_p>_T^2" is added afterhand.
   !     
   !
   ! ---
   !
   Shift1(1:NVar,1:NVar) = 0.D0
   DO i=1,NChol
      Shift1(1:NVar,1:NVar) = Shift1(1:NVar,1:NVar) + Chol(1:NVar,1:NVar,i)*AIMAG(MF_Shift(i))
   ENDDO
   IF (PrintT%Hamiltonian > 2) THEN
      WRITE(IW,'(3X,"Construct_One_Body_Propagator")')
      WRITE(IW,'(4X,"Shift1 = L_{ij,m} v_{m}")')
      CALL PrSq(Shift1,NVar,NVar,NVar)
   ENDIF
   !
   !
   ! 
   DO I=1,NSpin
      H1m(:,:,I) = H1(:,:,I) + Shift1(:,:)
   ENDDO
   IF (PrintT%Hamiltonian > 2) THEN
      WRITE(IW,'(4X,"H1m_{pq}   = H1_{pq} - L_{pq,m} hat{v}_{m}")')
      WRITE(IW,'(4X,"hat{v}_{m} = h_{pq} ap+.aq -  1/2 sum_r L_{pr,m} L_{qr,m}ap+.aq")')
      CALL PrSq(H1m,NVar,NVar,NVar)
   ENDIF
   !
   !
   ! 
   A = -0.5D0 * dt * H1m
   DO I=1,NSpin
      CALL AFQMC_ExpMat(NVar,-1,A(1,1,I),ExpH1(1,1,I))
   ENDDO
   IF (PrintT%Hamiltonian > 2) THEN
     WRITE(IW,'( "ExpH1m = [exp(-Δτ/2 H1[α]), exp(-Δτ/2 H1[β] ]")')
     CALL PrSq(ExpH1,NVar,NVar,NVar)
   ENDIF
END SUBROUTINE

SUBROUTINE AFQMC_Mean_Field_Shift_D(NVar,NSpin,NChol,Chol,G,MF_Shift)
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: NVar
  INTEGER,          INTENT(IN)  :: NSpin
  INTEGER,          INTENT(IN)  :: NChol
  DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
  DOUBLE PRECISION, INTENT(IN)  :: G(1:NVar,1:NVar,1:NSpin)
  DOUBLE PRECISION, INTENT(OUT) :: MF_Shift(1:NChol)
  DOUBLE PRECISION, EXTERNAL    :: DDOT
  INTEGER :: i
  WRITE(6,'("AFQMC_Mean_Field_Shift")')
  WRITE(6,'(" ipie.propagation.phaseless_base.construct_mean_field_shift (B)")')
  WRITE(6,'("  called by ipie.propagation.phaseless_base.build.")')
  IF (NSpin==1) THEN
     DO i=1,NChol
        MF_Shift(i) = 2.D0*DDOT(NVar*NVar,Chol(1,1,i),1,G(1,1,1),1)
        WRITE(6,'("MF_Shift(",I0,")=",F20.12,"i")')i,MF_Shift(i)
     ENDDO
  ELSE
     DO i=1,NChol
        MF_Shift(i) = DDOT(NVar*NVar,Chol(1,1,i),1,G(1,1,1),1)  &
      &             + DDOT(NVar*NVar,Chol(1,1,i),1,G(1,1,2),1)
        WRITE(6,'("MF_Shift(",I0,")=",F20.12,"i")')i,MF_Shift(i)
     ENDDO
  ENDIF
END SUBROUTINE


SUBROUTINE AFQMC_Mean_Field_Shift_Z(NVar,NSpin,NChol,Chol,G,MF_Shift)
  USE IO_Files,     ONLY: IW
  USE AFQMC_Module, ONLY: PrintT
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: NVar
  INTEGER,          INTENT(IN)  :: NSpin
  INTEGER,          INTENT(IN)  :: NChol
  DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
  DOUBLE COMPLEX,   INTENT(IN)  :: G(1:NVar,1:NVar,1:NSpin)
  DOUBLE COMPLEX,   INTENT(OUT) :: MF_Shift(1:NChol)
  DOUBLE PRECISION              :: GR(1:NVar,1:NVar,1:NSpin)
  DOUBLE PRECISION              :: GI(1:NVar,1:NVar,1:NSpin)
  DOUBLE PRECISION, EXTERNAL    :: DDOT
  DOUBLE PRECISION              :: DR, DI
  DOUBLE PRECISION, PARAMETER   :: Zero = 0D+00
  DOUBLE PRECISION, PARAMETER   :: One  = 1D+00
  DOUBLE PRECISION, PARAMETER   :: Two  = 2D+00
  INTEGER :: i, ISp
  !
  ! MF_Shift(m) 
  ! =   <Psi_T| hat{v}_m|Psi_T>
  ! =   sum_pq <Psi_T|i L^m_pq cp+.cq |Psi_T>/<Psi_T|Psi_T>
  ! =   i sum_pq G_pq L^m_pq
  ! =   sum_pq [ - IMG(G_pq) L^m_pq ] + i [ REAL(G_pq) L^m_pq ] 
  ! 
  IF (PrintT%Misc > 0) THEN 
     WRITE(IW,'(/,4X,"AFQMC_Mean_Field_Shift")')
     WRITE(IW,'(4X," ipie.propagation.phaseless_base.construct_mean_field_shift (B)")')
     WRITE(IW,'(4X,"    called by ipie.propagation.phaseless_base.build.",/)')
  ENDIF
  IF (PrintT%Misc > 1) THEN 
     WRITE(IW,'(4X,46("-"))')
     WRITE(IW,'(4X," IField          Mean-field shift")')
     WRITE(IW,'(4X,46("-"))')
  ENDIF
  GR = DBLE(G)
  GI = AIMAG(G)
  IF (NSpin==1) THEN
     DO i=1,NChol
        DR =  - Two*DDOT(NVar*NVar,Chol(1,1,i),1,GI(1,1,1),1)
        DI =  + Two*DDOT(NVar*NVar,Chol(1,1,i),1,GR(1,1,1),1)
        MF_Shift(i) = CMPLX(DR,DI,KIND=KIND(MF_Shift))
        IF (PrintT%Misc > 1) THEN 
           WRITE(IW,'(4X,I6,2F20.12)')i,MF_Shift(i)
        ENDIF
     ENDDO
  ELSE
     DO i=1,NChol
        DR = Zero
        DI = Zero
        DO ISp=1,NSpin
           DR = DR - DDOT(NVar*NVar,Chol(1,1,i),1,GI(1,1,ISp),1)
           DI = DI + DDOT(NVar*NVar,Chol(1,1,i),1,GR(1,1,ISp),1)
        ENDDO
        MF_Shift(i) = CMPLX(DR,DI,KIND=KIND(MF_Shift))
        IF (PrintT%Misc > 1) THEN 
           WRITE(IW,'(I6,2F20.12)')i,MF_Shift(i)
        ENDIF
     ENDDO
  ENDIF
  IF (PrintT%Misc > 1) THEN 
     WRITE(IW,'(4X,46("-"),/)')
  ENDIF
END SUBROUTINE

SUBROUTINE AFQMC_ExpMat(n,nterm,A,ExpA)
!
!  Infinite polynomial order
!    ExpMat   => afqmclib.F90
!
!  Finite polynomial order
!    ExpMatQC => mcscflib.F90
!
!--
  IMPLICIT NONE
  INTEGER, INTENT(IN)           :: n
  INTEGER, INTENT(IN)           :: nterm
  DOUBLE PRECISION, INTENT(IN)  :: A(1:n,1:n)
  DOUBLE PRECISION, INTENT(OUT) :: ExpA(1:n,1:n)
  IF (nterm < 0) THEN
     CALL ExpMat(n,A,ExpA)
  ELSE
     CALL ExpMatQC(n,nterm,0.D0,A,ExpA)
  ENDIF
  !CALL PrSq(ExpA,n,n,n)
END SUBROUTINE

SUBROUTINE AFQMC_Propagate_One_Body_Z(NVar,NSpin,ExpH1,UW)
!
!    Propagate by the kinetic term by direct matrix multiplication.
!    Only one spin component. Assuming phi is a batch.
!    For use with the continuus algorithm and free propagation.
!    Parameters
!    ----------
!    walker : :class:`pie.walker.Walker`
!        Walker object to be updated. on output we have acted on
!        :math:`|\phi_i\rangle` by :math:`B_{T/2}` and updated the weight
!        appropriately.  updates inplace.
!    state : :class:`pie.state.State`
!        Simulation state.
!
!   phi = xp.einsum("ik,wkj->wij", bt2, phi, optimize=True)
!      phi_{w,i,j} =  \sum_k [bt2]_{i,k}  [phi]_{w,k,j} 
!
!
!
!   hat{U}|Ψ> 
!       |ψ_p'> = hat{U}|ψ_p> = \sum_q |ψ_q><ψ_q|hat{U}|ψ_p> 
!              = \sum_q |ψ_q> U_qp
!
!   exp(-Δτ/2 H1')|Ψ'> 
!   =
!   exp(-Δτ/2 H1') \hat{U}|Ψ> 

   USE IO_Files,     ONLY: IW
   USE AFQMC_Module, ONLY: PrintT
   IMPLICIT NONE
   INTEGER,          INTENT(IN)    :: NVar
   INTEGER,          INTENT(IN)    :: NSpin
   DOUBLE PRECISION, INTENT(IN)    :: ExpH1(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX,   INTENT(INOUT) :: UW(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION                :: UR(1:NVar,1:NVar)
   DOUBLE PRECISION                :: UI(1:NVar,1:NVar)
   DOUBLE PRECISION                :: VR(1:NVar,1:NVar)
   DOUBLE PRECISION                :: VI(1:NVar,1:NVar)
   DOUBLE PRECISION, PARAMETER     :: Zero = 0D+00 
   DOUBLE PRECISION, PARAMETER     :: One = 1D+00 
   INTEGER                         :: ISpin
   !
   !
   !
   IF (PrintT%Walker > 1) THEN 
      WRITE(IW,'(3X,"AFQMC_Propagate_One_Body_Z")')
   ENDIF
   !
   !
   !
   DO ISpin=1,NSpin
      UR(1:NVar,1:NVar) = DBLE(UW(1:NVar,1:NVar,ISpin))
      UI(1:NVar,1:NVar) = AIMAG(UW(1:NVar,1:NVar,ISpin))
      ! conventional
#if 1
      CALL DGeMM('N','N',NVar,NVar,NVar,One,ExpH1(1,1,ISpin),NVar,UR,NVar,Zero,VR,NVar)
      CALL DGeMM('N','N',NVar,NVar,NVar,One,ExpH1(1,1,ISpin),NVar,UI,NVar,Zero,VI,NVar)
#else
      ! test
      CALL DGeMM('N','N',NVar,NVar,NVar,One,UR,NVar,ExpH1(1,1,ISpin),NVar,Zero,VR,NVar)
      CALL DGeMM('N','N',NVar,NVar,NVar,One,UI,NVar,ExpH1(1,1,ISpin),NVar,Zero,VI,NVar)
#endif
      UW(1:NVar,1:NVar,ISpin)=CMPLX(VR,VI,KIND=KIND(UW))
   ENDDO
   !
   !
   !
   IF (PrintT%Walker > 2) THEN 
      WRITE(IW,'(3X,A)')"ExpH1  = exp[-Δτ/2 H1']"
      CALL PrSq(ExpH1,NVar,NVar,NVar)
      WRITE(IW,'(3X,A)')"U(New) = exp[-Δτ/2 H1']U(old)"
      CALL PrSqZ(UW,NVar,NVar,NVar)
   ENDIF
END SUBROUTINE


SUBROUTINE AFQMC_Propagate_Two_Body_Z(NVar,NSpin,NChol,dT,x,Chol,UW)
   USE IO_Files,     ONLY: IW
   USE AFQMC_Module, ONLY: PrintT
   IMPLICIT NONE
   INTEGER,          INTENT(IN)    :: NVar
   INTEGER,          INTENT(IN)    :: NSpin
   INTEGER,          INTENT(IN)    :: NChol
   DOUBLE PRECISION, INTENT(IN)    :: dT
   DOUBLE COMPLEX,   INTENT(IN)    :: x(1:NChol)
   !DOUBLE COMPLEX,   INTENT(IN)    :: GH(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(IN)    :: Chol(1:NVar,1:NVar,1:NChol)
   DOUBLE COMPLEX,   INTENT(INOUT) :: UW(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX                  :: VHS(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX                  :: H2(1:NVar,1:NVar)
   DOUBLE COMPLEX                  :: NewUW(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX                  :: ExpH2(1:NVar,1:NVar)
   DOUBLE COMPLEX, PARAMETER       :: Zero = CMPLX(0D+00, 0D+00)
   DOUBLE COMPLEX, PARAMETER       :: One  = CMPLX(1D+00, 0D+00)
   DOUBLE PRECISION                :: XRI(1:NChol,1:2)
   DOUBLE PRECISION                :: H2Buf(1:NVar*NVar,1:2)
   DOUBLE PRECISION                :: SqDt
   DOUBLE PRECISION, PARAMETER     :: DZero = 0D+00, DOne = 1D+00
   INTEGER                         :: ISpin
   INTEGER                         :: I
   INTEGER                         :: NField
   !
   !
   !
   IF (PrintT%Walker > 1) THEN 
      WRITE(IW,'(3X,"AFQMC_Propagate_Two_Body_Z")')
   ENDIF
   !
   !
   !
   NField = NChol
   !
   !--- 1. Calculate Two-electron potential <Psi_T|V_m|Psi_T> ---
   ! 
   ! OpenQP port: the VHS build was a loop of NField calls, each doing two
   ! NVar^2 array expressions over one Cholesky vector -- NField passes over the
   ! tensor per walker per step, with temporaries created every call. The whole
   ! sum is a single matrix-vector product over Chol viewed as (NVar^2, NChol),
   ! done here as one DGEMM with the (-sqrt(dt)*Im x, +sqrt(dt)*Re x) pair as
   ! two right-hand sides: one pass over Chol, in BLAS.
   ! (AFQMC_Construct_Two_Body_Hamiltonian is kept for the upstream API; it is
   ! the same expression one vector at a time.)
   SqDt = DSQRT(dT)
   DO I = 1, NField
      XRI(I,1) = -SqDt*AIMAG(X(I))
      XRI(I,2) =  SqDt*DBLE(X(I))
   ENDDO
   CALL DGeMM('N','N',NVar*NVar,2,NField,DOne,Chol,NVar*NVar,XRI,NChol,&
  &           DZero,H2Buf,NVar*NVar)
   H2 = RESHAPE(CMPLX(H2Buf(:,1),H2Buf(:,2),KIND=KIND(H2)),(/NVar,NVar/))
   IF (PrintT%Walker > 2) THEN 
      WRITE(IW,'(3X,A)')" Δτ.H2' (=VHS)"
      CALL PrSqZ(H2,NVar,NVar,NVar)
   ENDIF
   !
   !--- 2. Form propagator ---
   !
   CALL ExpMat_Z(NVar,6,H2,ExpH2)
   IF (PrintT%Walker > 2) THEN 
      WRITE(IW,'(3X,A)')"Exp[Δτ H2']"
      CALL PrSqZ(ExpH2,NVar,NVar,NVar)
   ENDIF
   !
   !--- 3. Orbital rotation ---
   !
   DO ISpin=1,NSpin
      ! conventional
#if 1
      CALL ZGeMM('N','N',NVar,NVar,NVar,One,ExpH2,NVar,UW(1,1,ISpin),NVar,Zero,NewUW(1,1,ISpin),NVar)
#else
      ! test
      CALL ZGeMM('N','N',NVar,NVar,NVar,One,UW(1,1,ISpin),NVar,ExpH2,NVar,Zero,NewUW(1,1,ISpin),NVar)
#endif
   ENDDO
   !
   !--- 4. Update walkers ---
   !
   UW = NewUW
   IF (PrintT%Walker > 2) THEN 
      DO ISpin=1,NSpin
         WRITE(IW,'(3X,"ISpin=",I4)')ISpin
         WRITE(IW,'(3X,A)')"U(New) = exp[Δτ H2']"
         CALL PrSqZ(UW(1,1,ISpin),NVar,NVar,NVar)
      ENDDO
   ENDIF
END SUBROUTINE




SUBROUTINE AFQMC_Propagate_One_Body(NVar,NSpin,UW, ExpH1t2)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NVar
  INTEGER, INTENT(IN) :: NSpin
  DOUBLE PRECISION, INTENT(INOUT) :: UW(1:NVar,1:NVar,1:NSpin)
  DOUBLE PRECISION, INTENT(IN)    :: ExpH1t2(1:NVar,1:NVar,1:NSpin)
  DOUBLE PRECISION    :: Wrk(1:NVar,1:NVar)
!
!    Propagate by the kinetic term by direct matrix multiplication.
!    Only one spin component. Assuming phi is a batch.
!    For use with the continuus algorithm and free propagation.
!    Parameters
!    ----------
!    walker : :class:`pie.walker.Walker`
!        Walker object to be updated. on output we have acted on
!        :math:`|\phi_i\rangle` by :math:`B_{T/2}` and updated the weight
!        appropriately.  updates inplace.
!    state : :class:`pie.state.State`
!        Simulation state.
!
!   phi = xp.einsum("ik,wkj->wij", bt2, phi, optimize=True)
  INTEGER :: ISpin

    Wrk = 0.D0
    DO ISpin=1,NSpin
      CALL DGeMM('N','N',NVar,NVar,NVar,1.D0,ExpH1t2(1,1,ISpin),NVar,UW(1,1,ISpin),NVar,0.D0,Wrk,NVar)
      UW(1:NVar,1:NVar,ISpin)=Wrk(1:NVar,1:NVar)
    ENDDO
    WRITE(6,*)"U(New) = exp[Δτ/2 H1']"
    CALL PrSq(UW,NVar,NVar,NVar)
#if 0
       DO i=1,
          DO k=1,
             DO j=1,
                Wrk(i,j,m) = Wrk(i,j,m) + Bt2(i,k)*Phi(k,j,m)
             ENDDO
          ENDDO
       ENDDO
#endif
END SUBROUTINE


!SUBROUTINE AFQMC_Construct_Two_Body_Hamiltonian(NVar,NSpin,dT,X,GH,Chol,H2)
SUBROUTINE AFQMC_Construct_Two_Body_Hamiltonian(NVar,NSpin,dT,X,Chol,H2)
!
!   ipie.propagation.phaseless_base.PhaselessBase.propagate_walkers_two_body     [4]
!     ipie.propagation.phaseless_base.PhaselessBase.apply_VHS
!      <= ipie.propagation.phaseless_generic.PhaselessGeneric.apply_VHS          [5]
!          ipie.propagation.phaseless_generic.PhaselessGeneric.construct_VHS     [6]
!              [V_HS]_pq = i √Δτ [ Sum_m Re(x_m - xbar_m).Lpq,m  + j*Im(x_m - xbar_m).Lpq,m ]
!              v_pq,m = 
!          ipie.propagation.operations.apply_exponential [6]
!               phi(new) = Exp(V_HS).phi(old)
!
!
!
!   [V_m]_pq = <Psi_T|hat{V}_m|Psi_w>/<Psi_T|Psi_w>  
!            = i G^w_pq L^m_pq
!   を用いてるわけではない.むしろ
!   <V_m>_pq = <Psi_T|hat{V}_m|Psi_T>
!            = i L^m_pq
!   を用いている (2026.5.25)
!
!   したがってGHは不要
!
!  [x] [V_HS]_pq = √dt <Psi_T|1/2 (x_m - xbar_m) hat{V}_m|Psi_T> = i L^m_pq 
!  [o] [V_HS]_pq = √dt <Psi_T| (x_m - xbar_m) hat{V}_m|Psi_T> = i L^m_pq 
!  と思えるが、実際は[x]がipieに一致。なぜ1/2の因子が必要？？ (2026.5.30)
!
   IMPLICIT NONE
   INTEGER,          INTENT(IN)    :: NVar
   INTEGER,          INTENT(IN)    :: NSpin
   DOUBLE PRECISION, INTENT(IN)    :: dT
   DOUBLE COMPLEX,   INTENT(IN)    :: x
  ! DOUBLE COMPLEX,   INTENT(IN)    :: GH(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(IN)    :: Chol(1:NVar,1:NVar)
   DOUBLE COMPLEX,   INTENT(INOUT) :: H2(1:NVar,1:NVar)
   DOUBLE PRECISION                :: H2R(1:NVar,1:NVar)
   DOUBLE PRECISION                :: H2I(1:NVar,1:NVar)
   INTEGER                         :: ISpin
   !
   ! V_HS = i SQRT(dt)*x_shift(IField)*Chol(1:NVar,1:NVar,IField))       
   !
   H2R(1:NVar,1:NVar) = -DSQRT(dT)*AIMAG(x)*Chol(1:NVar,1:NVar)
   H2I(1:NVar,1:NVar) =  DSQRT(dT)*DBLE(x)*Chol(1:NVar,1:NVar)
   !
   !---  Add alpha and beta parts ---
   !
   ! Why -> Rotation alpha and beta respectively
   H2(1:NVar,1:NVar) = H2(1:NVar,1:NVar) + &
  & CMPLX(H2R(1:NVar,1:NVar),H2I(1:NVar,1:NVar),KIND=KIND(H2))
   !H2(1:NVar,1:NVar) = H2(1:NVar,1:NVar) + 2.0*0.5D0*CMPLX(H2R(1:NVar,1:NVar), H2I(1:NVar,1:NVar))
   !H2(1:NVar,1:NVar) = H2(1:NVar,1:NVar) + 2.D0*CMPLX(H2R(1:NVar,1:NVar), H2I(1:NVar,1:NVar))
END SUBROUTINE
