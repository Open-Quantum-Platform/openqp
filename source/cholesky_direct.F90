!> @file cholesky_direct.F90
!>
!> @brief Integral-direct pivoted Cholesky: factorise without ever storing
!>        the two-electron integrals.
!>
!> `cholesky_eri` reads a packed AO array costing nbf^4/8 to hold, so it
!> removes the v^4 ladder block but leaves that store as the remaining N^4
!> object.  This version never builds it.
!>
!> Pivoted Cholesky needs only the diagonal (mu nu | mu nu) and one column
!> (mu nu | pivot) per accepted vector.  Both can be produced on demand, so
!> the working set is the vectors themselves, O(npair * nchol).
!>
!> Extracting a column needs no change to the integral engine.
!> `int2_compute_t` asks its consumer whether a shell quartet is worth
!> computing, through `screen_ij` and `screen_ijkl`, and skips whatever the
!> consumer values below the cutoff.  A consumer that values only the quartets
!> touching the pivot's shell pair therefore receives exactly that column.
!>
!> Pivoting is by shell pair: one pass yields every column of that pair at
!> once, so a shell pair of n functions costs one pass rather than n.
!>
!> @note The AO extents are built from `num_ao()` rather than read from
!>   `basis%ao_offset` / `basis%naos`.  Those are Cartesian per-shell offsets
!>   (see `build_cart_density`), while the integral buffers index the AO space
!>   the calculation actually runs in.  `num_ao()` is the single place that
!>   resolves 5d against 6d, from the shell's pure flag and the global gate
!>   `[input] ispher` sets, so both conventions are followed automatically and
!>   a basis with d or higher is addressed correctly either way.
module cholesky_direct

  use precision, only: dp
  use int2_compute, only: int2_compute_data_t, int2_storage_t
  use basis_tools, only: basis_set

  implicit none

  private
  public :: cholesky_direct_decompose

  !> Per-shell extents in the AO space the integrals are delivered in.
  type :: ao_map_t
    integer, allocatable :: first(:)   !< first AO of each shell
    integer, allocatable :: count(:)   !< AOs in each shell
    integer, allocatable :: shell(:)   !< shell owning each AO
  end type ao_map_t

  !> Consumer keeping only (mu nu | mu nu), the pivoting diagonal.
  type, extends(int2_compute_data_t) :: diag_collect_t
    real(kind=dp), pointer :: diag(:) => null()
  contains
    procedure :: parallel_start => diag_start
    procedure :: parallel_stop  => diag_stop
    procedure :: update         => diag_update
    procedure :: clean          => diag_clean
    procedure :: screen_ij      => diag_screen_ij
    procedure :: screen_ijkl    => diag_screen_ijkl
  end type diag_collect_t

  !> Consumer keeping only the columns belonging to one shell pair.
  type, extends(int2_compute_data_t) :: column_collect_t
    real(kind=dp), pointer :: col(:,:) => null()
    integer :: want_i = 0, want_j = 0
    integer :: ao_i = 0, n_i = 0, ao_j = 0, n_j = 0
  contains
    procedure :: parallel_start => col_start
    procedure :: parallel_stop  => col_stop
    procedure :: update         => col_update
    procedure :: clean          => col_clean
    procedure :: screen_ij      => col_screen_ij
    procedure :: screen_ijkl    => col_screen_ijkl
  end type column_collect_t

contains

!###############################################################################
! Diagonal consumer
!###############################################################################

  subroutine diag_start(this, basis, nthreads)
    class(diag_collect_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nthreads
    if (.false.) then
      if (size(shape(basis%nbf))<0) continue
      if (size(shape(nthreads))<0) continue
    end if
  end subroutine diag_start

  subroutine diag_stop(this)
    class(diag_collect_t), intent(inout) :: this
    call this%pe%barrier()
    if (this%pe%size > 1) call this%pe%allreduce(this%diag, size(this%diag))
    call this%pe%barrier()
  end subroutine diag_stop

  subroutine diag_clean(this)
    class(diag_collect_t), intent(inout) :: this
    nullify(this%diag)
  end subroutine diag_clean

  function diag_screen_ij(this, xints, i, j) result(res)
    class(diag_collect_t), intent(in) :: this
    real(kind=dp), contiguous, intent(in) :: xints(:,:)
    real(kind=dp) :: res
    integer, intent(in) :: i, j
    res = xints(i,j)
    if (.false.) then; if (size(shape(this%pe%size))<0) continue; end if
  end function diag_screen_ij

  !> Only diagonal shell quartets can carry diagonal function pairs.
  function diag_screen_ijkl(this, xints, i, j, k, l) result(res)
    class(diag_collect_t), intent(in) :: this
    real(kind=dp), contiguous, intent(in) :: xints(:,:)
    real(kind=dp) :: res
    integer, intent(in) :: i, j, k, l
    if (i == k .and. j == l) then
      res = xints(i,j)*xints(k,l)
    else
      res = 0.0_dp
    end if
    if (.false.) then; if (size(shape(this%pe%size))<0) continue; end if
  end function diag_screen_ijkl

  subroutine diag_update(this, buf)
    class(diag_collect_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    integer :: n, i, j, k, l, p, q
    real(kind=dp) :: val
    do n = 1, buf%ncur
      i = int(buf%ids(1,n)); j = int(buf%ids(2,n))
      k = int(buf%ids(3,n)); l = int(buf%ids(4,n))
      p = pair_index(i, j); q = pair_index(k, l)
      if (p /= q) cycle
      val = buf%ints(n)
      if (i == j) val = val*2.0_dp
      if (k == l) val = val*2.0_dp
      val = val*2.0_dp                     ! the p == q coincidence
      this%diag(p) = val
    end do
    buf%ncur = 0
  end subroutine diag_update

!###############################################################################
! Column consumer
!###############################################################################

  subroutine col_start(this, basis, nthreads)
    class(column_collect_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nthreads
    if (associated(this%col)) this%col = 0.0_dp
    if (.false.) then
      if (size(shape(basis%nbf))<0) continue
      if (size(shape(nthreads))<0) continue
    end if
  end subroutine col_start

  subroutine col_stop(this)
    class(column_collect_t), intent(inout) :: this
    call this%pe%barrier()
    if (this%pe%size > 1) call this%pe%allreduce(this%col(:,1), size(this%col))
    call this%pe%barrier()
  end subroutine col_stop

  subroutine col_clean(this)
    class(column_collect_t), intent(inout) :: this
    nullify(this%col)
  end subroutine col_clean

  !> The wanted pair can be the ket of any bra, so no bra may be rejected here
  !> on the pivot's account; only the ordinary Schwarz bound applies.
  function col_screen_ij(this, xints, i, j) result(res)
    class(column_collect_t), intent(in) :: this
    real(kind=dp), contiguous, intent(in) :: xints(:,:)
    real(kind=dp) :: res
    integer, intent(in) :: i, j
    res = xints(i,j)
    if (.false.) then; if (size(shape(this%want_i))<0) continue; end if
  end function col_screen_ij

  !> Where the column is selected: a quartet survives only if one of its pairs
  !> is the wanted one.  Everything else is valued at zero and skipped before
  !> any integral is computed.
  function col_screen_ijkl(this, xints, i, j, k, l) result(res)
    class(column_collect_t), intent(in) :: this
    real(kind=dp), contiguous, intent(in) :: xints(:,:)
    real(kind=dp) :: res
    integer, intent(in) :: i, j, k, l
    if ((i == this%want_i .and. j == this%want_j) .or. &
        (k == this%want_i .and. l == this%want_j)) then
      res = xints(i,j)*xints(k,l)
    else
      res = 0.0_dp
    end if
  end function col_screen_ijkl

  subroutine col_update(this, buf)
    class(column_collect_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    integer :: n, i, j, k, l, p, q, m
    real(kind=dp) :: val
    do n = 1, buf%ncur
      i = int(buf%ids(1,n)); j = int(buf%ids(2,n))
      k = int(buf%ids(3,n)); l = int(buf%ids(4,n))
      val = buf%ints(n)
      p = pair_index(i, j); q = pair_index(k, l)
      if (i == j) val = val*2.0_dp
      if (k == l) val = val*2.0_dp
      if (p == q) val = val*2.0_dp
      m = slot_of(this, i, j)
      if (m > 0) this%col(q, m) = val
      m = slot_of(this, k, l)
      if (m > 0) this%col(p, m) = val
    end do
    buf%ncur = 0
  end subroutine col_update

  !> Column index within the block for AO pair (mu,nu), or 0 if outside it.
  pure integer function slot_of(this, mu, nu) result(m)
    class(column_collect_t), intent(in) :: this
    integer, intent(in) :: mu, nu
    m = 0
    if (in_range(mu, this%ao_i, this%n_i) .and. &
        in_range(nu, this%ao_j, this%n_j)) then
      m = (mu - this%ao_i + 1) + (nu - this%ao_j)*this%n_i
    else if (in_range(nu, this%ao_i, this%n_i) .and. &
             in_range(mu, this%ao_j, this%n_j)) then
      m = (nu - this%ao_i + 1) + (mu - this%ao_j)*this%n_i
    end if
  end function slot_of

  pure logical function in_range(x, first, n) result(y)
    integer, intent(in) :: x, first, n
    y = x >= first .and. x < first + n
  end function in_range

  pure integer function pair_index(mu, nu) result(p)
    integer, intent(in) :: mu, nu
    if (mu >= nu) then
      p = mu*(mu-1)/2 + nu
    else
      p = nu*(nu-1)/2 + mu
    end if
  end function pair_index

!###############################################################################
! Driver
!###############################################################################

!> Per-shell AO extents in the space the integral buffers use.
!>
!> Aborts if they do not sum to nbf: that would mean this mapping and the
!> engine disagree about what an AO index means, and every column would be
!> scattered to the wrong place -- silently.
  subroutine build_ao_map(basis, map)
    use constants, only: num_ao
    use messages, only: show_message, with_abort
    type(basis_set), intent(in) :: basis
    type(ao_map_t), intent(out) :: map
    integer :: s, n, k

    allocate(map%first(basis%nshell), map%count(basis%nshell))
    allocate(map%shell(basis%nbf))

    k = 0
    do s = 1, basis%nshell
      ! num_ao() is the one place that knows the 5d/6d rule: it returns 2l+1
      ! only when the shell is flagged pure AND the global gate that [input]
      ! ispher sets is on.  Deciding it here instead would be a second copy of
      ! that rule, free to drift from the engine's.
      n = num_ao(basis%am(s), basis%harmonic(s))
      map%first(s) = k + 1
      map%count(s) = n
      if (k + n > basis%nbf) exit
      map%shell(k+1:k+n) = s
      k = k + n
    end do

    if (k /= basis%nbf) then
      call show_message('cholesky_direct: shell AO extents do not sum to nbf; &
                        &the AO mapping and the integral engine disagree', &
                        with_abort)
    end if
  end subroutine build_ao_map

  pure subroutine pair_to_ao(p, mu, nu)
    integer, intent(in) :: p
    integer, intent(out) :: mu, nu
    ! Largest mu with mu(mu-1)/2 < p.  Start from the real root and walk; the
    ! walk is at most one step, but it is what makes this exact.
    mu = int((sqrt(8.0_dp*real(p,dp) + 1.0_dp) - 1.0_dp)*0.5_dp) + 1
    do while (mu > 1 .and. mu*(mu-1)/2 >= p)
      mu = mu - 1
    end do
    do while ((mu+1)*mu/2 < p)
      mu = mu + 1
    end do
    nu = p - mu*(mu-1)/2
  end subroutine pair_to_ao

!> @brief Integral-direct pivoted Cholesky.
!>
!> @param[inout] basis    AO basis
!> @param[inout] infos    calculation state, for the integral driver
!> @param[in]    tol      stop once the largest remaining diagonal is below it
!> @param[in]    maxchol  capacity of @p lvec
!> @param[out]   lvec     the vectors, lvec(P,J)
!> @param[out]   nchol    how many were produced
!> @param[out]   err      largest remaining diagonal on exit
!> @param[out]   truncated .true. if capacity ran out before the tolerance
!> @param[out]   npass    integral passes used: one for the diagonal, then one
!>                        per pivot block
  subroutine cholesky_direct_decompose(basis, infos, tol, maxchol, lvec, &
                                       nchol, err, truncated, npass)

    use types, only: information
    use int2_compute, only: int2_compute_t

    type(basis_set), intent(inout) :: basis
    type(information), target, intent(inout) :: infos
    real(dp), intent(in) :: tol
    integer, intent(in) :: maxchol
    real(dp), intent(out) :: lvec(:,:)
    integer, intent(out) :: nchol
    real(dp), intent(out) :: err
    logical, intent(out) :: truncated
    integer, intent(out) :: npass

    type(int2_compute_t) :: drv
    type(diag_collect_t), target :: dcons
    type(column_collect_t), target :: ccons
    type(ao_map_t) :: map
    real(dp), allocatable, target :: colblk(:,:), diag(:)
    integer :: npair, p, j, pivot, mu, nu, si, sj, sw, nblk, slot
    real(dp) :: vmax, scale_

    npair = basis%nbf*(basis%nbf+1)/2
    nchol = 0
    npass = 0
    truncated = .false.

    call build_ao_map(basis, map)
    allocate(diag(npair), source=0.0_dp)

    call drv%init(basis, infos)
    call drv%set_screening()

    dcons%diag => diag
    call drv%run(dcons)
    call dcons%clean()
    npass = 1

    err = maxval(diag)

    outer: do
      pivot = maxloc(diag, dim=1)
      vmax = diag(pivot)
      if (vmax < tol) exit outer
      if (nchol >= maxchol) then
        truncated = .true.
        exit outer
      end if

      ! The shell pair owning this pivot.  Because the block is chosen FROM
      ! the pivot, the pivot is always inside it -- which is what makes the
      ! inner loop guaranteed to accept at least one vector and the outer loop
      ! guaranteed to make progress.
      call pair_to_ao(pivot, mu, nu)
      si = map%shell(mu)
      sj = map%shell(nu)
      if (si < sj) then
        sw = si; si = sj; sj = sw
      end if

      ccons%want_i = si
      ccons%want_j = sj
      ccons%ao_i = map%first(si); ccons%n_i = map%count(si)
      ccons%ao_j = map%first(sj); ccons%n_j = map%count(sj)
      nblk = ccons%n_i*ccons%n_j

      if (allocated(colblk)) deallocate(colblk)
      allocate(colblk(npair, nblk), source=0.0_dp)
      ccons%col => colblk

      call drv%run(ccons)
      call ccons%clean()
      npass = npass + 1

      ! Drain this block: its columns are already in hand, so taking every
      ! vector it still supports is free relative to fetching another.
      inner: do
        pivot = maxloc(diag, dim=1)
        vmax = diag(pivot)
        if (vmax < tol) exit outer
        if (nchol >= maxchol) then
          truncated = .true.
          exit outer
        end if

        call pair_to_ao(pivot, mu, nu)
        slot = slot_of(ccons, mu, nu)
        if (slot == 0) exit inner     ! the best pivot moved to another block

        nchol = nchol + 1
        do p = 1, npair
          lvec(p, nchol) = colblk(p, slot)
        end do
        do j = 1, nchol-1
          call daxpy(npair, -lvec(pivot,j), lvec(1,j), 1, lvec(1,nchol), 1)
        end do
        scale_ = 1.0_dp/sqrt(vmax)
        do p = 1, npair
          lvec(p,nchol) = lvec(p,nchol)*scale_
          diag(p) = max(diag(p) - lvec(p,nchol)*lvec(p,nchol), 0.0_dp)
        end do
        diag(pivot) = 0.0_dp
      end do inner
    end do outer

    err = maxval(diag)
    call drv%clean()
    if (allocated(colblk)) deallocate(colblk)
    deallocate(diag)

  end subroutine cholesky_direct_decompose

end module cholesky_direct
