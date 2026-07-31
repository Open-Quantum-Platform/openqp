#define UNUSED_DUMMY(x) if (.false.) then ; if (size(shape(x))<0) continue ; end if
!> @file cc_ao2mo.F90
!>
!> @brief AO -> MO two-electron integral transformation for the coupled-cluster
!>        module, plus extraction of the occupied/virtual blocks CCSD(T) needs.
!>
!> OpenQP's two-electron engine is consumer-driven: `int2_compute_t` walks the
!> shell quartets and hands each thread a buffer of unique integrals, and the
!> consumer decides what to do with them.  The Fock consumers contract the
!> buffer on the fly; here we instead scatter it into a dense AO array, because
!> a conventional CCSD needs random access to every (pq|rs).
!>
!> Two details of the buffer contract matter:
!>
!>   * `storeints` emits ONE canonical representative per 8-fold orbit and
!>     pre-scales it by 1/2 for each index coincidence (ii=jj, kk=ll, and the
!>     bra/ket pair exchange).  That convention exists so the Fock builders can
!>     use a flat six-term expansion.  A verbatim tensor needs the undressed
!>     value, so the scaling is undone here before the orbit is expanded.
!>   * Orbits are disjoint, so expanding each representative to its 8 permuted
!>     positions writes every AO address exactly once.  The writes are plain
!>     assignments rather than accumulations, which makes the scatter safe from
!>     the OpenMP threads of the integral driver without any locking.
!>
!> The transformation itself is the textbook N^5 sequence of four quarter
!> transformations, each a single DGEMM over a flattened index triple.
module cc_ao2mo

  use precision, only: dp
  use int2_compute, only: int2_compute_data_t, int2_storage_t
  use basis_tools, only: basis_set

  implicit none

  private
  public :: cc_eri_collect_t
  public :: cc_build_mo_blocks

  !> Integral consumer that materialises the full AO integral tensor.
  type, extends(int2_compute_data_t) :: cc_eri_collect_t
    real(kind=dp), pointer :: ao(:,:,:,:) => null()
    integer :: nbf = 0
  contains
    procedure :: parallel_start => cc_eri_parallel_start
    procedure :: parallel_stop  => cc_eri_parallel_stop
    procedure :: update         => cc_eri_update
    procedure :: clean          => cc_eri_clean
    procedure :: screen_ij      => cc_eri_screen_ij
    procedure :: screen_ijkl    => cc_eri_screen_ijkl
  end type cc_eri_collect_t

contains

!###############################################################################

  subroutine cc_eri_parallel_start(this, basis, nthreads)
    class(cc_eri_collect_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nthreads
    this%nbf = basis%nbf
    if (associated(this%ao)) this%ao = 0.0_dp
    UNUSED_DUMMY(nthreads)
  end subroutine cc_eri_parallel_start

!###############################################################################

  subroutine cc_eri_parallel_stop(this)
    class(cc_eri_collect_t), intent(inout) :: this
    ! Shell pairs are distributed across ranks by the driver, so each rank holds
    ! a disjoint slice of the tensor; one reduction completes it everywhere.
    call this%pe%barrier()
    if (this%pe%size > 1) then
      call this%pe%allreduce(this%ao, size(this%ao))
    end if
    call this%pe%barrier()
  end subroutine cc_eri_parallel_stop

!###############################################################################

  subroutine cc_eri_clean(this)
    class(cc_eri_collect_t), intent(inout) :: this
    nullify(this%ao)
  end subroutine cc_eri_clean

!###############################################################################

!> Schwarz bound on the bra pair.  The base class does no screening at all; a
!> density-free |(ij|ij)|^1/2 bound is the correct one for a bare integral list.
  function cc_eri_screen_ij(this, xints, i, j) result(res)
    class(cc_eri_collect_t), intent(in) :: this
    real(kind=dp), contiguous, intent(in) :: xints(:,:)
    real(kind=dp) :: res
    integer, intent(in) :: i, j
    res = xints(i,j)
    UNUSED_DUMMY(this)
  end function cc_eri_screen_ij

  function cc_eri_screen_ijkl(this, xints, i, j, k, l) result(res)
    class(cc_eri_collect_t), intent(in) :: this
    real(kind=dp), contiguous, intent(in) :: xints(:,:)
    real(kind=dp) :: res
    integer, intent(in) :: i, j, k, l
    res = xints(i,j)*xints(k,l)
    UNUSED_DUMMY(this)
  end function cc_eri_screen_ijkl

!###############################################################################

!> Scatter one buffer into the AO tensor, undoing the storeints coincidence
!> scaling and expanding the 8-fold permutational orbit.
  subroutine cc_eri_update(this, buf)
    class(cc_eri_collect_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf

    integer :: n, i, j, k, l
    real(kind=dp) :: val

    do n = 1, buf%ncur
      i = int(buf%ids(1,n))
      j = int(buf%ids(2,n))
      k = int(buf%ids(3,n))
      l = int(buf%ids(4,n))
      val = buf%ints(n)

      ! Undo the 1/2 factors applied per index coincidence in storeints.
      if (i == j) val = val*2.0_dp
      if (k == l) val = val*2.0_dp
      if (i == k .and. j == l) val = val*2.0_dp

      this%ao(i,j,k,l) = val
      this%ao(j,i,k,l) = val
      this%ao(i,j,l,k) = val
      this%ao(j,i,l,k) = val
      this%ao(k,l,i,j) = val
      this%ao(l,k,i,j) = val
      this%ao(k,l,j,i) = val
      this%ao(l,k,j,i) = val
    end do

    buf%ncur = 0

  end subroutine cc_eri_update

!###############################################################################

!> @brief Transform the AO integrals to the MO basis and slice out the blocks
!>        consumed by the coupled-cluster equations.
!>
!> @param[in]  nbf   number of basis functions
!> @param[in]  nmo   number of MOs kept (frozen orbitals already removed)
!> @param[in]  no    correlated occupied orbitals (first @p no columns of cmo)
!> @param[in]  cmo   MO coefficients of the retained orbitals, (nbf, nmo)
!> @param[in]  ao    full AO integral tensor (nbf^4), chemist notation
!>
!> Blocks are returned in chemist notation with the index orders assumed by
!> cc_lib: oooo(i,j,k,l)=(ij|kl), ooov(i,j,k,a)=(ij|ka), oovv(i,j,a,b)=(ij|ab),
!> ovov(i,a,j,b)=(ia|jb), ovvv(i,a,b,c)=(ia|bc), vvvv(a,b,c,d)=(ab|cd).
  subroutine cc_build_mo_blocks(nbf, nmo, no, cmo, ao, &
                                oooo, ooov, oovv, ovov, ovvv, vvvv)

    integer, intent(in) :: nbf, nmo, no
    real(dp), intent(in) :: cmo(nbf, nmo)
    real(dp), intent(in) :: ao(nbf, nbf, nbf, nbf)
    real(dp), intent(out) :: oooo(no,no,no,no)
    real(dp), intent(out) :: ooov(no,no,no,nmo-no)
    real(dp), intent(out) :: oovv(no,no,nmo-no,nmo-no)
    real(dp), intent(out) :: ovov(no,nmo-no,no,nmo-no)
    real(dp), intent(out) :: ovvv(no,nmo-no,nmo-no,nmo-no)
    real(dp), intent(out) :: vvvv(nmo-no,nmo-no,nmo-no,nmo-no)

    real(dp), allocatable :: x(:,:,:,:), y(:,:,:,:)
    integer :: i, j, k, l, a, b, c, d, nv

    nv = nmo - no

    ! Four quarter transformations.  Each contracts the leading AO index and
    ! leaves the new MO index trailing, so after four passes the chemist index
    ! order is restored -- the same rotation trick used for the AO tensor build.
    allocate(x(nbf,nbf,nbf,nmo))
    call dgemm('t','n', nbf*nbf*nbf, nmo, nbf, 1.0_dp, ao, nbf, cmo, nbf, &
               0.0_dp, x, nbf*nbf*nbf)

    allocate(y(nbf,nbf,nmo,nmo))
    call dgemm('t','n', nbf*nbf*nmo, nmo, nbf, 1.0_dp, x, nbf, cmo, nbf, &
               0.0_dp, y, nbf*nbf*nmo)
    deallocate(x)

    allocate(x(nbf,nmo,nmo,nmo))
    call dgemm('t','n', nbf*nmo*nmo, nmo, nbf, 1.0_dp, y, nbf, cmo, nbf, &
               0.0_dp, x, nbf*nmo*nmo)
    deallocate(y)

    allocate(y(nmo,nmo,nmo,nmo))
    call dgemm('t','n', nmo*nmo*nmo, nmo, nbf, 1.0_dp, x, nbf, cmo, nbf, &
               0.0_dp, y, nmo*nmo*nmo)
    deallocate(x)

    ! --- slice the blocks -------------------------------------------------
    !$omp parallel do collapse(3) private(i,j,k,l) schedule(static)
    do l = 1, no
      do k = 1, no
        do j = 1, no
          do i = 1, no
            oooo(i,j,k,l) = y(i,j,k,l)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do collapse(3) private(i,j,k,a) schedule(static)
    do a = 1, nv
      do k = 1, no
        do j = 1, no
          do i = 1, no
            ooov(i,j,k,a) = y(i,j,k,no+a)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do collapse(3) private(i,j,a,b) schedule(static)
    do b = 1, nv
      do a = 1, nv
        do j = 1, no
          do i = 1, no
            oovv(i,j,a,b) = y(i,j,no+a,no+b)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do collapse(3) private(i,j,a,b) schedule(static)
    do b = 1, nv
      do j = 1, no
        do a = 1, nv
          do i = 1, no
            ovov(i,a,j,b) = y(i,no+a,j,no+b)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do collapse(3) private(i,a,b,c) schedule(static)
    do c = 1, nv
      do b = 1, nv
        do a = 1, nv
          do i = 1, no
            ovvv(i,a,b,c) = y(i,no+a,no+b,no+c)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do collapse(3) private(a,b,c,d) schedule(static)
    do d = 1, nv
      do c = 1, nv
        do b = 1, nv
          do a = 1, nv
            vvvv(a,b,c,d) = y(no+a,no+b,no+c,no+d)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    deallocate(y)

  end subroutine cc_build_mo_blocks

!###############################################################################

end module cc_ao2mo
