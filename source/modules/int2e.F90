module int2e_mod

  use precision, only: dp
  use int2_compute, only: int2_compute_data_t

  implicit none

  character(len=*), parameter :: module_name = "int2e_mod"

  private
  public int2e

!> @brief Consumer that scatters computed shell-quartet ERIs into a dense
!>        (nbf,nbf,nbf,nbf) AO tensor, applying the full 8-fold permutational
!>        symmetry. Mirrors the int2_rhf_data_t consumer but accumulates the
!>        raw integrals instead of contracting them into a Fock matrix.
  type, extends(int2_compute_data_t) :: int2_dump_data_t
    integer :: nbf = 0
    real(kind=dp), pointer :: eri(:,:,:,:) => null()
  contains
    procedure :: parallel_start => dump_parallel_start
    procedure :: parallel_stop => dump_parallel_stop
    procedure :: update => dump_update
    procedure :: clean => dump_clean
  end type

contains

  subroutine int2e_C(c_handle) bind(C, name="int2e")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call int2e(inf)
  end subroutine int2e_C

!> @brief Compute all two-electron repulsion integrals (mu nu|la si) in the AO
!>        basis (chemist notation) and store them in the OQP::ERI_AO tag as a
!>        full nbf**4 array, so the Python layer can build a FCIDUMP / qubit
!>        Hamiltonian. This is the conventional (in-core) path: memory grows as
!>        nbf**4, so it is intended for small active systems, not production SCF.
  subroutine int2e(infos)

    use types, only: information
    use oqp_tagarray_driver
    use precision, only: dp
    use io_constants, only: iw
    use basis_tools, only: basis_set
    use printing, only: print_module_info
    use messages, only: show_message, WITH_ABORT
    use int2_compute, only: int2_compute_t

    implicit none

    character(len=*), parameter :: subroutine_name = "int2e"

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis

    type(int2_compute_t) :: int2_driver
    type(int2_dump_data_t) :: dump
    real(kind=dp), contiguous, pointer :: eri_flat(:)
    integer :: nbf
    integer(8) :: nbf4
    real(kind=dp) :: mem_mb

    open(unit=iw, file=infos%log_filename, position="append")

    basis => infos%basis
    basis%atoms => infos%atoms

    call print_module_info('int2e', &
      'Computing Two-Electron Repulsion Integrals (AO ERIs)')

    nbf = basis%nbf
    nbf4 = int(nbf,8)**4
    mem_mb = real(nbf4,dp) * 8.0d0 / (1024.0d0*1024.0d0)

    write(iw,'(/1x,"AO basis functions (nbf): ",i0)') nbf
    write(iw,'(1x,"In-core ERI tensor size : ",i0," elements (",f0.1," MB)")') &
      nbf4, mem_mb

!   Guard against an accidental, ruinous allocation. nbf**4 doubles is the
!   conventional in-core cost; refuse clearly above ~16 GB rather than thrash.
    if (mem_mb > 16384.0d0) then
      call show_message( &
        "int2e: in-core AO ERI tensor exceeds 16 GB; this routine targets "// &
        "small active systems for FCIDUMP export, not full production basis "// &
        "sets.", WITH_ABORT)
    end if

!   Allocate and zero the destination tag. Screened (negligible) integrals are
!   left at zero, matching standard quantum-chemistry practice.
    call infos%dat%remove_records([character(len=80) :: OQP_ERI_AO])
    call infos%dat%reserve_data(OQP_ERI_AO, TA_TYPE_REAL64, nbf4, &
                                comment=OQP_ERI_AO_comment)
    call tagarray_get_data(infos%dat, OQP_ERI_AO, eri_flat)
    eri_flat = 0.0d0

!   Remap the flat storage to a rank-4 view for convenient scatter. The full
!   8-fold symmetry of the integrals makes the C/Fortran index-order difference
!   irrelevant: the Python side reshapes the same bytes to (pq|rs) directly.
    dump%nbf = nbf
    dump%eri(1:nbf,1:nbf,1:nbf,1:nbf) => eri_flat

!   Drive the conventional two-electron engine with the dump consumer. No CAM
!   attenuation: we want the bare 1/r12 Coulomb integrals.
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    call int2_driver%run(dump)
    call int2_driver%clean()

    dump%eri => null()

    write(iw,"(/1x,'...... End Of Two-Electron Integrals ......'/)")
    close(iw)

  end subroutine int2e

!###############################################################################

  subroutine dump_parallel_start(this, basis, nthreads)
    use basis_tools, only: basis_set
    implicit none
    class(int2_dump_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nthreads
!   Nothing to set up: distinct shell quartets map to disjoint canonical AO
!   index sets, so threads never write the same destination element.
    if (.false.) then
      this%nbf = this%nbf
      if (basis%nbf < 0 .or. nthreads < 0) continue
    end if
  end subroutine dump_parallel_start

  subroutine dump_parallel_stop(this)
    implicit none
    class(int2_dump_data_t), intent(inout) :: this

!   int2_compute_t distributes shell-quartet work across MPI ranks.  Each rank
!   scatters only its local quartets into the dense tensor, so combine the full
!   tensor before returning it through OQP::ERI_AO.  For non-MPI runs this is a
!   no-op through par_env_t%allreduce.
    call this%pe%barrier()
    if (associated(this%eri)) then
      call this%pe%allreduce(this%eri, size(this%eri))
    end if
    call this%pe%barrier()
  end subroutine dump_parallel_stop

  subroutine dump_clean(this)
    implicit none
    class(int2_dump_data_t), intent(inout) :: this
    if (.false.) this%nbf = this%nbf
  end subroutine dump_clean

!> @brief Scatter one buffer of unique integrals into the dense tensor using
!>        the 8-fold permutational symmetry of (ij|kl).
  subroutine dump_update(this, buf)
    use int2_compute, only: int2_storage_t
    implicit none
    class(int2_dump_data_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    integer :: n, i, j, k, l
    real(kind=dp) :: v

    do n = 1, buf%ncur
      i = buf%ids(1,n)
      j = buf%ids(2,n)
      k = buf%ids(3,n)
      l = buf%ids(4,n)
      v = buf%ints(n)

!     storeints (int2_compute_data_t_storeints) pre-scales the buffered value
!     by 0.5 for each "diagonal" coincidence so the Fock build can apply
!     uniform Coulomb/exchange factors. Undo that scaling here to recover the
!     true integral (mu nu|la si) before scattering it into the dense tensor.
      if (i == j) v = v*2.0d0
      if (k == l) v = v*2.0d0
      if (i == k .and. j == l) v = v*2.0d0

      this%eri(i,j,k,l) = v
      this%eri(j,i,k,l) = v
      this%eri(i,j,l,k) = v
      this%eri(j,i,l,k) = v
      this%eri(k,l,i,j) = v
      this%eri(l,k,i,j) = v
      this%eri(k,l,j,i) = v
      this%eri(l,k,j,i) = v
    end do

    buf%ncur = 0
  end subroutine dump_update

end module int2e_mod
