SUBROUTINE AFQMC_Force_Bias_D(NVar,NSpin,NChol,dT,GH,Chol,MF_Shift,Xbar)
!
!   xbar =  - √Δτ [ <ψ_T|hat{v}_m|ψ_W>/<ψ_T|ψ_W> - <ψ_T|hat{v}_m|ψ_T>/<ψ_T|ψ_T> ]
!       
!
!
   USE IO_Files, ONLY: IW
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NSpin
   INTEGER, INTENT(IN) :: NVar
   INTEGER, INTENT(IN) :: NChol
   DOUBLE PRECISION, INTENT(IN)  :: dT
   DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
   DOUBLE PRECISION, INTENT(IN)  :: GH(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(IN)  :: MF_Shift(1:NChol)
   DOUBLE PRECISION, INTENT(OUT) :: XBar(1:NChol)
   DOUBLE PRECISION              :: Vm(1:NChol)
   DOUBLE PRECISION              :: VmS
   DOUBLE PRECISION, EXTERNAL    :: DDOT
   INTEGER :: IChol, ISpin, i, j
   DO IChol=1,NChol 
      VmS = 0.D0
      DO ISpin=1,NSpin 
         VmS = VmS + DDOT(NVar*NVar,GH(1,1,ISpin),1,Chol(1,1,IChol),1)
         !DO i=1,NVar
         !   DO j=1,NVar
         !      XBar(IChol) = XBar(IChol) + GH(i,j,ISpin)*Chol(j,i,IChol)
         !   ENDDO
         !ENDDO
      ENDDO
      Vm(IChol) = VmS
   ENDDO
   !
   !  xbar = -self.sqrt_dt * (1j * self.vbias - self.mf_shift)
   !
   XBar(1:NChol) = (Vm(1:NChol)-MF_Shift(1:NChol))*DSQRT(dT)
   WRITE(IW,'(3X,"Force bias")')
   DO IChol=1,NChol
      WRITE(IW,'(I4,F20.12)')IChol,XBar(IChol)
   ENDDO
END SUBROUTINE

SUBROUTINE AFQMC_Force_Bias_Z(NVar,NSpin,NChol,dT,GH,Chol,MF_Shift,Xbar)
! ipie.propagation.phaseless_base.PhaselessBase.propagate_walkers                [3]
!   #
!   # (2.b) Apply two-body
!   #
!   ipie.propagation.phaseless_base.PhaselessBase.propagate_walkers_two_body     [4]
!     ipie.trial_wavefunction.single_det.SingleDet.calc_force_bias_A             [5]
!       ipie.propagation.force_bias.construct_force_bias_batch_single_det (Real) [6]
!         各Walkerに対してHalf Rotate Greens FunctionとCholVecを縮約
!            v_wm = GH_ri.Lri,m
!
!  Mean-field substraction is included 
!
!  xbar(m,w) = - sqrt(dt) [<Psi_T|v_m|Psi_w>/<Psi_T|Psi_w> -<Psi_T|v_m|Psi_T>/<Psi_T|Psi_T> ]
!            = - sqrt(dt) [ i G^w_pq L^m_pq - MF_Shift(m)] 
!
   USE IO_Files,     ONLY: IW
   USE AFQMC_Module, ONLY: PrintT
   IMPLICIT NONE
   INTEGER,          INTENT(IN)  :: NSpin
   INTEGER,          INTENT(IN)  :: NVar
   INTEGER,          INTENT(IN)  :: NChol
   DOUBLE PRECISION, INTENT(IN)  :: dT
   DOUBLE PRECISION, INTENT(IN)  :: Chol(1:NVar,1:NVar,1:NChol)
   DOUBLE COMPLEX,   INTENT(IN)  :: GH(1:NVar,1:NVar,1:NSpin)
   DOUBLE COMPLEX,   INTENT(IN)  :: MF_Shift(1:NChol)
   DOUBLE COMPLEX,   INTENT(OUT) :: XBar(1:NChol)
   DOUBLE COMPLEX                :: Vm(1:NChol)
   DOUBLE PRECISION              :: GSum(1:NVar*NVar,1:2)
   DOUBLE PRECISION              :: Proj(1:NChol,1:2)
   DOUBLE PRECISION, PARAMETER   :: Zero=0D+00, One=1D+00
   INTEGER :: IChol, ISpin
   !
   ! OpenQP port: this used to be 2*NSpin*NChol separate DDOT calls, i.e. one
   ! BLAS-1 pass over the whole Cholesky tensor per (vector, spin, re/im). Chol
   ! is NVar^2*NChol doubles and is re-streamed by every walker at every step,
   ! so this kernel is memory-bandwidth bound and the traffic was the cost.
   !
   ! Summing the spins first and stacking (Re, Im) as two right-hand sides turns
   ! the whole contraction into ONE DGEMM over Chol viewed as (NVar^2, NChol):
   ! a single pass instead of 2*NSpin, and BLAS-2/3 instead of BLAS-1.
   !
   ! The spin sum now happens before the contraction rather than after, so the
   ! rounding differs from the previous order at the 1e-16 level. It stays
   ! deterministic and independent of the thread count.
   !
   GSum(:,1) = RESHAPE(DBLE(GH(:,:,1)), (/NVar*NVar/))
   GSum(:,2) = RESHAPE(AIMAG(GH(:,:,1)), (/NVar*NVar/))
   DO ISpin = 2, NSpin
      GSum(:,1) = GSum(:,1) + RESHAPE(DBLE(GH(:,:,ISpin)), (/NVar*NVar/))
      GSum(:,2) = GSum(:,2) + RESHAPE(AIMAG(GH(:,:,ISpin)), (/NVar*NVar/))
   ENDDO
   CALL DGeMM('T','N',NChol,2,NVar*NVar,One,Chol,NVar*NVar,GSum,NVar*NVar,&
  &           Zero,Proj,NChol)
   DO IChol=1,NChol
      Vm(IChol) = CMPLX(-Proj(IChol,2),Proj(IChol,1),KIND=KIND(Vm))
   ENDDO
   !bar{x}_W = - √Δτ <ψ_T|hat{v}-<v>|ψ_W>/<ψ_T|ψ_W>
   XBar(1:NChol) = -(Vm(1:NChol)-MF_Shift(1:NChol))*DSQRT(dT)

   IF (PrintT%Walker > 1) THEN 
      WRITE(IW,'(3X,"AFQMC_Force_Bias_Z")')
      WRITE(IW,'(6X,"  Force bias: bar{x}_W = - √Δτ <ψ_T|hat{v}-<v>|ψ_W>/<ψ_T|ψ_W>")')
      WRITE(IW,'(A)')
      WRITE(IW,'(3X,"--------------------")')
      WRITE(IW,'(3X," IField  Force Bias ")')
      WRITE(IW,'(3X,"--------------------")')
      DO IChol=1,NChol
         WRITE(IW,'(I6,2F20.12)')IChol,XBar(IChol)
      ENDDO
      WRITE(IW,'(3X,"--------------------")')
      WRITE(IW,'(A)')
   ENDIF
END SUBROUTINE
