module fci_integrals_mod

  use precision, only: dp
  use int2_compute, only: int2_compute_t, int2_compute_data_t, int2_storage_t
  use basis_tools, only: basis_set

  implicit none

  character(len=*), parameter :: module_name = "fci_integrals_mod"

  private

  public fci_ao_integrals

  type, extends(int2_compute_data_t) :: fci_eri_data_t
    integer :: nbf = 0
    real(kind=dp), pointer :: eri(:) => null()
  contains
    procedure :: parallel_start => fci_eri_parallel_start
    procedure :: parallel_stop => fci_eri_parallel_stop
    procedure :: update => fci_eri_update
    procedure :: clean => fci_eri_clean
  end type fci_eri_data_t

contains

  subroutine fci_ao_integrals_C(c_handle) bind(C, name="fci_ao_integrals")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call fci_ao_integrals(inf)
  end subroutine fci_ao_integrals_C

  subroutine fci_ao_integrals(infos)
    use types, only: information
    use messages, only: show_message, WITH_ABORT
    use oqp_tagarray_driver, only: &
      OQP_AO_ERI, OQP_AO_ERI_comment, TA_OK

    implicit none

    character(len=*), parameter :: subroutine_name = "fci_ao_integrals"

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis
    type(int2_compute_t) :: int2_driver
    type(fci_eri_data_t) :: eri_data
    integer :: nbf
    integer :: neri
    integer :: ta_status
    integer(8) :: neri64
    real(kind=dp) :: eri_mb
    character(len=256) :: msg

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf

    if (nbf <= 0) call show_message("FCI AO ERI requested before basis setup", WITH_ABORT)

    ! Dense AO ERI storage scales as nbf**4. Guard against 32-bit overflow of the
    ! element count and against an unreasonable allocation: dense FCI is meant for
    ! small systems only, so a few-GiB ceiling is the right backstop here.
    neri64 = int(nbf, 8)**4
    eri_mb = real(neri64, dp) * 8.0_dp / (1024.0_dp**2)
    if (neri64 > int(huge(neri), 8) .or. eri_mb > 4096.0_dp) then
      write(msg, '(A,I0,A,F0.1,A)') &
        "FCI dense AO ERI request is too large (nbf=", nbf, ", ", eri_mb, &
        " MiB); dense FCI is intended for small active spaces only."
      call show_message(trim(msg), WITH_ABORT)
    end if
    neri = int(neri64)

    ! tagarray 1.0.0 API: alloc binds the pointer and override=.true. replaces any
    ! existing record (subsumes the old remove_records + reserve_data + get trio).
    ta_status = infos%dat%alloc(OQP_AO_ERI, [neri], eri_data%eri, &
                                description=OQP_AO_ERI_comment, override=.true.)
    if (ta_status /= TA_OK) &
      call show_message("FCI: failed to allocate the AO ERI tagarray record", WITH_ABORT)
    eri_data%eri = 0.0_dp
    eri_data%nbf = nbf

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    call int2_driver%run(eri_data)
    call int2_driver%clean()
    call eri_data%clean()

  end subroutine fci_ao_integrals

  subroutine fci_eri_parallel_start(this, basis, nthreads)
    implicit none
    class(fci_eri_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nthreads

    this%nbf = basis%nbf

    if (.false.) then
      if (nthreads < 0) continue
    end if
  end subroutine fci_eri_parallel_start

  subroutine fci_eri_parallel_stop(this)
    implicit none
    class(fci_eri_data_t), intent(inout) :: this

    call this%pe%barrier()
    call this%pe%allreduce(this%eri, size(this%eri))
    call this%pe%barrier()
  end subroutine fci_eri_parallel_stop

  subroutine fci_eri_clean(this)
    implicit none
    class(fci_eri_data_t), intent(inout) :: this

    nullify(this%eri)
    this%nbf = 0
  end subroutine fci_eri_clean

  subroutine fci_eri_update(this, buf)
    implicit none
    class(fci_eri_data_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf

    integer :: n
    integer :: ii, jj, kk, ll
    real(kind=dp) :: val

!$omp critical(fci_eri_update_region)
    do n = 1, buf%ncur
      ii = int(buf%ids(1, n))
      jj = int(buf%ids(2, n))
      kk = int(buf%ids(3, n))
      ll = int(buf%ids(4, n))
      val = buf%ints(n)

      ! int2 stores permutation-unique integrals with symmetry weights for
      ! Fock builds. Undo those weights before expanding to the full tensor.
      if (ii == jj) val = 2.0_dp * val
      if (kk == ll) val = 2.0_dp * val
      if (ii == kk .and. jj == ll) val = 2.0_dp * val

      call fci_eri_set_permutations(this, ii, jj, kk, ll, val)
    end do
!$omp end critical(fci_eri_update_region)

    buf%ncur = 0
  end subroutine fci_eri_update

  subroutine fci_eri_set_permutations(this, i, j, k, l, val)
    implicit none
    class(fci_eri_data_t), intent(inout) :: this
    integer, intent(in) :: i, j, k, l
    real(kind=dp), intent(in) :: val

    call fci_eri_set(this, i, j, k, l, val)
    call fci_eri_set(this, j, i, k, l, val)
    call fci_eri_set(this, i, j, l, k, val)
    call fci_eri_set(this, j, i, l, k, val)
    call fci_eri_set(this, k, l, i, j, val)
    call fci_eri_set(this, l, k, i, j, val)
    call fci_eri_set(this, k, l, j, i, val)
    call fci_eri_set(this, l, k, j, i, val)
  end subroutine fci_eri_set_permutations

  subroutine fci_eri_set(this, i, j, k, l, val)
    implicit none
    class(fci_eri_data_t), intent(inout) :: this
    integer, intent(in) :: i, j, k, l
    real(kind=dp), intent(in) :: val
    integer :: idx

    idx = i &
        + (j - 1) * this%nbf &
        + (k - 1) * this%nbf * this%nbf &
        + (l - 1) * this%nbf * this%nbf * this%nbf
    this%eri(idx) = val
  end subroutine fci_eri_set

end module fci_integrals_mod
