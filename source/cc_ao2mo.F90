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
!> transformations driven through the symmetric pair index, so neither the AO
!> integrals nor the intermediate is ever stored as a dense nbf^4 tensor.
module cc_ao2mo

  use precision, only: dp
  use int2_compute, only: int2_compute_data_t, int2_storage_t
  use basis_tools, only: basis_set

  implicit none

  private
  public :: cc_eri_collect_t
  public :: cc_build_mo_blocks
  public :: cc_packed_length
  public :: cc_build_full_mo

  !> Integral consumer that materialises the AO integrals in packed form.
  !>
  !> Only the canonical eighth of the tensor is kept: g(idx(P,Q)) with
  !> P = pair(mu,nu) >= Q = pair(lambda,sigma) is (mu nu | lambda sigma).
  !> The dense nbf^4 tensor the transformation used to build is by far the
  !> largest array in a CCSD(T) run -- 12 GB at nbf = 200 against 1.5 GB here --
  !> and the driver never needs an element that this form cannot supply.
  type, extends(int2_compute_data_t) :: cc_eri_collect_t
    real(kind=dp), pointer :: g(:) => null()
    integer :: nbf = 0
    integer :: npair = 0
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
    this%npair = basis%nbf*(basis%nbf+1)/2
    if (associated(this%g)) this%g = 0.0_dp
    UNUSED_DUMMY(nthreads)
  end subroutine cc_eri_parallel_start

!###############################################################################

  subroutine cc_eri_parallel_stop(this)
    class(cc_eri_collect_t), intent(inout) :: this
    ! Shell pairs are distributed across ranks by the driver, so each rank holds
    ! a disjoint slice of the tensor; one reduction completes it everywhere.
    call this%pe%barrier()
    if (this%pe%size > 1) then
      ! Reduce in chunks: the packed length passes the 2^31 that the
      ! allreduce count argument holds at around nbf = 400, and a wrapped
      ! count would silently reduce the wrong extent rather than fail.
      block
        integer(8) :: ntot, off
        integer(8), parameter :: CHUNK = 134217728_8   ! 2^27 elements, 1 GB
        integer :: nthis
        ntot = size(this%g, kind=8)
        off = 0
        do while (off < ntot)
          nthis = int(min(CHUNK, ntot - off))
          call this%pe%allreduce(this%g(off+1:off+nthis), nthis)
          off = off + nthis
        end do
      end block
    end if
    call this%pe%barrier()
  end subroutine cc_eri_parallel_stop

!###############################################################################

  subroutine cc_eri_clean(this)
    class(cc_eri_collect_t), intent(inout) :: this
    nullify(this%g)
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

!> Transform the packed AO integrals into a FULL spatial MO tensor,
!>   eri(p,q,r,s) = (pq|rs),
!> with the bra pair (pq) over @p cmo_bra and the ket pair (rs) over
!> @p cmo_ket.  Passing the same coefficients twice gives a same-spin tensor;
!> passing alpha then beta gives the mixed-spin one.
!>
!> Same two half transformations as cc_build_mo_blocks -- only the second half
!> writes a full tensor instead of scattering into the six coupled-cluster
!> blocks, because the spin-orbital assembly needs every (pq|rs).
subroutine cc_build_full_mo(nbf, nmo, cmo_bra, cmo_ket, g, eri)

  integer, intent(in) :: nbf, nmo
  real(dp), intent(in) :: cmo_bra(nbf,nmo), cmo_ket(nbf,nmo)
  real(dp), intent(in) :: g(*)
  real(dp), intent(out) :: eri(nmo,nmo,nmo,nmo)

  real(dp), allocatable :: half(:,:)
  integer, allocatable :: prow(:), pcol(:), qrow(:), qcol(:)
  integer :: npair, nmopair, pq, p, q, ip

  npair   = nbf*(nbf+1)/2
  nmopair = nmo*(nmo+1)/2

  allocate(prow(npair), pcol(npair), qrow(nmopair), qcol(nmopair))
  do p = 1, nbf
    do q = 1, p
      ip = p*(p-1)/2 + q
      prow(ip) = p; pcol(ip) = q
    end do
  end do
  do p = 1, nmo
    do q = 1, p
      ip = p*(p-1)/2 + q
      qrow(ip) = p; qcol(ip) = q
    end do
  end do

  allocate(half(nmopair, npair))

  !$omp parallel default(shared) private(pq, p, q, ip)
  block
    real(dp), allocatable :: d(:,:), scr(:,:), m(:,:)
    allocate(d(nbf,nbf), scr(nbf,nmo), m(nmo,nmo))
    !$omp do schedule(static)
    do pq = 1, npair
      do ip = 1, npair
        p = prow(ip); q = pcol(ip)
        d(p,q) = g(packed_index(ip, pq))
        d(q,p) = d(p,q)
      end do
      call dgemm('n','n', nbf, nmo, nbf, 1.0_dp, d, nbf, cmo_bra, nbf, 0.0_dp, scr, nbf)
      call dgemm('t','n', nmo, nmo, nbf, 1.0_dp, cmo_bra, nbf, scr, nbf, 0.0_dp, m, nmo)
      do ip = 1, nmopair
        half(ip, pq) = m(qrow(ip), qcol(ip))
      end do
    end do
    !$omp end do
    deallocate(d, scr, m)
  end block
  !$omp end parallel

  !$omp parallel default(shared) private(pq, p, q, ip)
  block
    real(dp), allocatable :: e(:,:), scr(:,:), n(:,:)
    allocate(e(nbf,nbf), scr(nbf,nmo), n(nmo,nmo))
    !$omp do schedule(static)
    do pq = 1, nmopair
      do ip = 1, npair
        p = prow(ip); q = pcol(ip)
        e(p,q) = half(pq, ip)
        e(q,p) = e(p,q)
      end do
      call dgemm('n','n', nbf, nmo, nbf, 1.0_dp, e, nbf, cmo_ket, nbf, 0.0_dp, scr, nbf)
      call dgemm('t','n', nmo, nmo, nbf, 1.0_dp, cmo_ket, nbf, scr, nbf, 0.0_dp, n, nmo)
      p = qrow(pq); q = qcol(pq)
      eri(p,q,:,:) = n
      if (p /= q) eri(q,p,:,:) = n
    end do
    !$omp end do
    deallocate(e, scr, n)
  end block
  !$omp end parallel

  deallocate(half, prow, pcol, qrow, qcol)

end subroutine cc_build_full_mo

!###############################################################################

!> Canonical index of the AO pair (mu,nu).
  pure integer function pair_index(mu, nu) result(p)
    integer, intent(in) :: mu, nu
    if (mu >= nu) then
      p = mu*(mu-1)/2 + nu
    else
      p = nu*(nu-1)/2 + mu
    end if
  end function pair_index

!> Offset of the pair-pair element (P,Q) in the packed triangular store.
!> int64 throughout -- the packed length passes 2^31 near nbf = 400.
  pure integer(8) function packed_index(p, q) result(idx)
    integer, intent(in) :: p, q
    integer(8) :: hi, lo
    hi = int(max(p,q), 8)
    lo = int(min(p,q), 8)
    idx = hi*(hi-1)/2 + lo
  end function packed_index

!> Number of elements in the packed AO store for a given basis size.
  pure integer(8) function cc_packed_length(nbf) result(n)
    integer, intent(in) :: nbf
    integer(8) :: np
    np = int(nbf, 8)*int(nbf+1, 8)/2
    n = np*(np+1)/2
  end function cc_packed_length

!###############################################################################

!> Store one buffer, undoing the storeints coincidence scaling.  The driver
!> hands over each canonical quartet once and storeints has already sorted the
!> indices into the canonical order, so a single packed slot receives it and
!> the 8-fold orbit is implicit in the storage.
  subroutine cc_eri_update(this, buf)
    class(cc_eri_collect_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf

    integer :: n, i, j, k, l, p, q
    real(kind=dp) :: val

    do n = 1, buf%ncur
      i = int(buf%ids(1,n))
      j = int(buf%ids(2,n))
      k = int(buf%ids(3,n))
      l = int(buf%ids(4,n))
      val = buf%ints(n)

      p = pair_index(i, j)
      q = pair_index(k, l)

      ! Undo the 1/2 factors applied per index coincidence in storeints.
      if (i == j) val = val*2.0_dp
      if (k == l) val = val*2.0_dp
      if (p == q) val = val*2.0_dp

      this%g(packed_index(p, q)) = val
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
!> @param[in]  g     packed AO integrals from cc_eri_collect_t, chemist notation
!>
!> Blocks are returned in chemist notation with the index orders assumed by
!> cc_lib: oooo(i,j,k,l)=(ij|kl), ooov(i,j,k,a)=(ij|ka), oovv(i,j,a,b)=(ij|ab),
!> ovov(i,a,j,b)=(ia|jb), ovvv(i,a,b,c)=(ia|bc), vvvv(a,b,c,d)=(ab|cd).
  subroutine cc_build_mo_blocks(nbf, nmo, no, cmo, g, &
                                oooo, ooov, oovv, ovov, ovvv, vvvv)

    integer, intent(in) :: nbf, nmo, no
    real(dp), intent(in) :: cmo(nbf, nmo)
    real(dp), intent(in) :: g(*)
    real(dp), intent(out) :: oooo(no,no,no,no)
    real(dp), intent(out) :: ooov(no,no,no,nmo-no)
    real(dp), intent(out) :: oovv(no,no,nmo-no,nmo-no)
    real(dp), intent(out) :: ovov(no,nmo-no,no,nmo-no)
    real(dp), intent(out) :: ovvv(no,nmo-no,nmo-no,nmo-no)
    real(dp), intent(out) :: vvvv(nmo-no,nmo-no,nmo-no,nmo-no)

    real(dp), allocatable :: half(:,:)
    integer, allocatable :: prow(:), pcol(:), qrow(:), qcol(:)
    integer :: npair, nmopair, nv, pq, p, q, ip

    nv = nmo - no
    npair   = nbf*(nbf+1)/2
    nmopair = nmo*(nmo+1)/2

    ! Pair -> (row,col) tables, so the unpack loops below index the packed
    ! stores by a single running counter instead of recomputing a triangular
    ! root per element.
    allocate(prow(npair), pcol(npair), qrow(nmopair), qcol(nmopair))
    do p = 1, nbf
      do q = 1, p
        ip = p*(p-1)/2 + q
        prow(ip) = p; pcol(ip) = q
      end do
    end do
    do p = 1, nmo
      do q = 1, p
        ip = p*(p-1)/2 + q
        qrow(ip) = p; qcol(ip) = q
      end do
    end do

    ! Two half transformations instead of four quarter transformations.  A
    ! quarter transformation has to hold a full nbf^3 x nmo intermediate and
    ! its successor at the same time; going through the symmetric pair index
    ! keeps both the input and the intermediate at an eighth and a quarter of
    ! the dense tensor, which is what lets the whole run fit in memory.
    allocate(half(nmopair, npair))

    ! --- first half: (mu nu | Q) -> (p q | Q) for every ket pair Q ---------
    !$omp parallel default(shared) private(pq, p, q, ip)
    block
      real(dp), allocatable :: d(:,:), scr(:,:), m(:,:)
      allocate(d(nbf,nbf), scr(nbf,nmo), m(nmo,nmo))
      !$omp do schedule(static)
      do pq = 1, npair
        do ip = 1, npair
          p = prow(ip); q = pcol(ip)
          d(p,q) = g(packed_index(ip, pq))
          d(q,p) = d(p,q)
        end do
        call dgemm('n','n', nbf, nmo, nbf, 1.0_dp, d, nbf, cmo, nbf, 0.0_dp, scr, nbf)
        call dgemm('t','n', nmo, nmo, nbf, 1.0_dp, cmo, nbf, scr, nbf, 0.0_dp, m, nmo)
        do ip = 1, nmopair
          half(ip, pq) = m(qrow(ip), qcol(ip))
        end do
      end do
      !$omp end do
      deallocate(d, scr, m)
    end block
    !$omp end parallel

    ! --- second half: (p q | lambda sigma) -> (p q | r s) -----------------
    !$omp parallel default(shared) private(pq, p, q, ip)
    block
      real(dp), allocatable :: e(:,:), scr(:,:), n(:,:)
      allocate(e(nbf,nbf), scr(nbf,nmo), n(nmo,nmo))
      !$omp do schedule(static)
      do pq = 1, nmopair
        do ip = 1, npair
          p = prow(ip); q = pcol(ip)
          e(p,q) = half(pq, ip)
          e(q,p) = e(p,q)
        end do
        call dgemm('n','n', nbf, nmo, nbf, 1.0_dp, e, nbf, cmo, nbf, 0.0_dp, scr, nbf)
        call dgemm('t','n', nmo, nmo, nbf, 1.0_dp, cmo, nbf, scr, nbf, 0.0_dp, n, nmo)

        p = qrow(pq); q = qcol(pq)
        call scatter_pair(p, q, n)
        if (p /= q) call scatter_pair(q, p, n)
      end do
      !$omp end do
      deallocate(e, scr, n)
    end block
    !$omp end parallel

    deallocate(half, prow, pcol, qrow, qcol)

  contains

    !> File the bra pair (p,q) of a completed (pq|rs) panel into the block that
    !> stores it.  A virtual-occupied bra needs no case of its own: every block
    !> keeping a mixed pair stores it occupied-first, so the transposed call
    !> covers it.
    subroutine scatter_pair(p, q, n)
      integer, intent(in) :: p, q
      real(dp), intent(in) :: n(nmo,nmo)
      integer :: i, j, a, b

      if (p <= no .and. q <= no) then
        i = p; j = q
        oooo(i,j,:,:) = n(1:no,      1:no)
        ooov(i,j,:,:) = n(1:no,      no+1:nmo)
        oovv(i,j,:,:) = n(no+1:nmo,  no+1:nmo)
      else if (p <= no .and. q > no) then
        i = p; a = q - no
        ovov(i,a,:,:) = n(1:no,      no+1:nmo)
        ovvv(i,a,:,:) = n(no+1:nmo,  no+1:nmo)
      else if (p > no .and. q > no) then
        a = p - no; b = q - no
        vvvv(a,b,:,:) = n(no+1:nmo,  no+1:nmo)
      end if
    end subroutine scatter_pair

  end subroutine cc_build_mo_blocks

!###############################################################################

end module cc_ao2mo
