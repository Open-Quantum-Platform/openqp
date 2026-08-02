!> @file afqmc_host.F90
!>
!> @brief Host-environment shims for the vendored AFQMC kernels.
!>
!> The AFQMC sources under source/afqmc/ are vendored from
!> Open-Quantum-Platform/AFQMC, where they are written against the Gellan
!> host program.  This file supplies the handful of host modules and utility
!> subroutines they `USE`/`CALL`, implemented against OpenQP instead:
!>
!>   IO_Files      -> OpenQP's log unit
!>   MPI_Parallel  -> serial no-ops (liboqp parallelises with OpenMP)
!>   AFQMC_Module  -> the walker/control/print derived types
!>   Expnd/PrSq/DiagMRRR/ExpMatQC -> Gellan mathlib utilities
!>
!> Keeping the shims here rather than editing every vendored file lets the
!> AFQMC sources stay close to their upstream form, so they can be re-synced.

module IO_Files
  implicit none
  integer :: IW = 6
end module IO_Files

!> Serial stand-ins for Gellan's MPI layer.  liboqp is an OpenMP library; the
!> AFQMC walker loop is threaded, not distributed, so every collective here is
!> a copy and every rank query answers "rank 0 of 1".
module MPI_Parallel
  implicit none
  integer :: my_rank = 0, NProc = 1, MPI_Comm_World = 0, my_node_comm = 0, MPI_SUM = 0
  logical :: main_rank = .true., writ_rank = .true.
  interface Gel_AllReduce
    module procedure AllReduce_Real_Scalar, AllReduce_Real_Vector
    module procedure AllReduce_Real_Vector_Scalar
    module procedure AllReduce_Complex_Scalar, AllReduce_Complex_Vector
  end interface
  interface Gel_BCast
    module procedure BCast_Integer_Scalar, BCast_Integer_Vector
    module procedure BCast_Integer_Matrix, BCast_Real_Vector
  end interface
  interface Gel_Sendrecv
    module procedure Sendrecv_Complex_Vector, Sendrecv_Complex_Cube
    module procedure Sendrecv_Complex_Scalar
    module procedure Sendrecv_Real_Scalar
  end interface
contains
  subroutine AllReduce_Real_Scalar(A, B, N, TypeCode, Op, Comm, Flag)
    double precision, intent(in) :: A
    double precision, intent(out) :: B
    integer, intent(in) :: N, Op, Comm, Flag
    character(len=*), intent(in) :: TypeCode
    B = A
  end subroutine
  subroutine AllReduce_Real_Vector(A, B, N, TypeCode, Op, Comm, Flag)
    double precision, intent(in) :: A(*)
    double precision, intent(out) :: B(*)
    integer, intent(in) :: N, Op, Comm, Flag
    character(len=*), intent(in) :: TypeCode
    B(1:N) = A(1:N)
  end subroutine
  subroutine AllReduce_Real_Vector_Scalar(A, B, N, TypeCode, Op, Comm, Flag)
    double precision, intent(inout) :: A(:)
    double precision, intent(out) :: B
    integer, intent(in) :: N, Op, Comm, Flag
    character(len=*), intent(in) :: TypeCode
    B = sum(A(1:N))
  end subroutine
  subroutine AllReduce_Complex_Scalar(A, B, N, TypeCode, Op, Comm, Flag)
    double complex, intent(in) :: A
    double complex, intent(out) :: B
    integer, intent(in) :: N, Op, Comm, Flag
    character(len=*), intent(in) :: TypeCode
    B = A
  end subroutine
  subroutine AllReduce_Complex_Vector(A, B, N, TypeCode, Op, Comm, Flag)
    double complex, intent(in) :: A(*)
    double complex, intent(out) :: B(*)
    integer, intent(in) :: N, Op, Comm, Flag
    character(len=*), intent(in) :: TypeCode
    B(1:N) = A(1:N)
  end subroutine
  subroutine Gel_Barrier(Comm)
    integer, intent(in) :: Comm
  end subroutine
  subroutine BCast_Integer_Scalar(A, N, TypeCode, Root, Comm)
    integer, intent(inout) :: A
    integer, intent(in) :: N, Root, Comm
    character(len=*), intent(in) :: TypeCode
  end subroutine
  subroutine BCast_Integer_Vector(A, N, TypeCode, Root, Comm)
    integer, intent(inout) :: A(*)
    integer, intent(in) :: N, Root, Comm
    character(len=*), intent(in) :: TypeCode
  end subroutine
  subroutine BCast_Integer_Matrix(A, N, TypeCode, Root, Comm)
    integer, intent(inout) :: A(:,:)
    integer, intent(in) :: N, Root, Comm
    character(len=*), intent(in) :: TypeCode
  end subroutine
  subroutine BCast_Real_Vector(A, N, TypeCode, Root, Comm)
    double precision, intent(inout) :: A(*)
    integer, intent(in) :: N, Root, Comm
    character(len=*), intent(in) :: TypeCode
  end subroutine
  subroutine Sendrecv_Complex_Vector(Send, N, TypeCode, IDest, SendTag, Recv, NRecv, &
                                     ISource, RecvTag, Comm)
    double complex, intent(in) :: Send(*)
    double complex, intent(out) :: Recv(*)
    integer, intent(in) :: N, NRecv, IDest, SendTag, ISource, RecvTag, Comm
    character(len=*), intent(in) :: TypeCode
    Recv(1:min(N,NRecv)) = Send(1:min(N,NRecv))
  end subroutine
  subroutine Sendrecv_Complex_Cube(Send, N, TypeCode, IDest, SendTag, Recv, NRecv, &
                                   ISource, RecvTag, Comm)
    double complex, intent(in) :: Send(:,:,:)
    double complex, intent(out) :: Recv(:,:,:)
    integer, intent(in) :: N, NRecv, IDest, SendTag, ISource, RecvTag, Comm
    character(len=*), intent(in) :: TypeCode
    Recv = Send
  end subroutine
  subroutine Sendrecv_Complex_Scalar(Send, N, TypeCode, IDest, SendTag, Recv, NRecv, &
                                     ISource, RecvTag, Comm)
    double complex, intent(in) :: Send
    double complex, intent(out) :: Recv
    integer, intent(in) :: N, NRecv, IDest, SendTag, ISource, RecvTag, Comm
    character(len=*), intent(in) :: TypeCode
    Recv = Send
  end subroutine
  subroutine Sendrecv_Real_Scalar(Send, N, TypeCode, IDest, SendTag, Recv, NRecv, &
                                  ISource, RecvTag, Comm)
    double precision, intent(in) :: Send
    double precision, intent(out) :: Recv
    integer, intent(in) :: N, NRecv, IDest, SendTag, ISource, RecvTag, Comm
    character(len=*), intent(in) :: TypeCode
    Recv = Send
  end subroutine
end module MPI_Parallel

!> Walker/control/print derived types shared by every AFQMC kernel.
module AFQMC_Module
  implicit none
  type AFQMC_Walker_Type
    double complex, allocatable :: UW(:,:,:)
    double complex :: E_Hyb = (0.d0,0.d0), E_Loc = (0.d0,0.d0), Ovlp = (0.d0,0.d0)
    double precision :: Weight = 0.d0, Det = 0.d0
  end type
  type AFQMC_Control_Type
    integer :: NStep, NStep_Stblz, NStep_PopCntrl, NStep_Estimate
    integer :: NStep_Accum, NWlk, ISeed
    double precision :: dT, Weight_Min, Weight_Max, E_Bound, x_bound, Chol_Tol
  end type
  type Print_Type
    integer :: QMCLog = 0, PopCntrl = 0, Walker = 0, RandNum = 0, Energy = 0
    integer :: Hamiltonian = 0, Misc = 0, Integral = 0
  end type
  type(Print_Type) :: PrintT
end module AFQMC_Module

!> Expand a column-wise packed upper triangle into a full symmetric matrix.
subroutine Expnd(Packed, Full, N, Mode)
  implicit none
  integer, intent(in) :: N, Mode
  double precision, intent(in) :: Packed(*)
  double precision, intent(out) :: Full(N,N)
  integer :: I, J, IJ
  IJ = 0
  do J = 1, N
    do I = 1, J
      IJ = IJ + 1
      Full(I,J) = Packed(IJ)
      Full(J,I) = Packed(IJ)
    end do
  end do
end subroutine Expnd

!> Debug matrix print. Every call site is behind an AFQMC PrintT guard that is
!> off in production, so this stays cheap and unconditional-format.
subroutine PrSq(A, NRow, NCol, LdA)
  use IO_Files, only: IW
  implicit none
  integer, intent(in) :: NRow, NCol, LdA
  double precision, intent(in) :: A(LdA,*)
  integer :: I, J
  do I = 1, NRow
    write(IW,'(1x,10(1x,f12.6))') (A(I,J), J = 1, NCol)
  end do
end subroutine PrSq

!> Symmetric eigensolve, Gellan's name for it. A is not overwritten.
subroutine DiagMRRR(N, A, E, V)
  use IO_Files, only: IW
  implicit none
  integer, intent(in) :: N
  double precision, intent(in) :: A(N,N)
  double precision, intent(out) :: E(N), V(N,N)
  double precision, allocatable :: Work(:)
  integer :: LWork, Info
  LWork = max(1, 34*N)
  allocate(Work(LWork))
  V = A
  call DSYEV('V', 'U', N, V, N, E, Work, LWork, Info)
  deallocate(Work)
  if (Info /= 0) then
    write(IW,'(1x,"AFQMC DiagMRRR: DSYEV info=",i0," n=",i0," lwork=",i0)') Info, N, LWork
    write(IW,'(1x,"matrix finite: ",l1,"  max|a|=",es12.4)') &
      all(A == A), maxval(abs(A))
    error stop 'AFQMC: DSYEV failed in DiagMRRR'
  end if
end subroutine DiagMRRR

!> Truncated Taylor series for exp(A), Gellan's mcscflib entry point.
subroutine ExpMatQC(N, NTerm, Threshold, A, ExpA)
  implicit none
  integer, intent(in) :: N, NTerm
  double precision, intent(in) :: Threshold, A(N,N)
  double precision, intent(out) :: ExpA(N,N)
  double precision, allocatable :: Term(:,:), Next(:,:)
  integer :: I, K
  allocate(Term(N,N), Next(N,N))
  ExpA = 0.d0
  Term = 0.d0
  do I = 1, N
    ExpA(I,I) = 1.d0
    Term(I,I) = 1.d0
  end do
  do K = 1, NTerm
    call DGEMM('N', 'N', N, N, N, 1.d0/dble(K), Term, N, A, N, 0.d0, Next, N)
    ExpA = ExpA + Next
    Term = Next
    if (Threshold > 0.d0 .and. maxval(abs(Term)) <= Threshold) exit
  end do
  deallocate(Term, Next)
end subroutine ExpMatQC
