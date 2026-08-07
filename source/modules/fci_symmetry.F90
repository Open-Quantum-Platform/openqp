!> @file fci_symmetry.F90
!>
!> @brief Point-group classification of determinants and CI roots.
!>
!> A determinant is a product of occupied spin orbitals, so in an abelian
!> group its irrep is the direct product of the irreps of those orbitals.
!> With every irrep encoded as the bitmask of operations whose character is
!> -1 (which pyoqp stages as OQP::sym_irrep_xor), that product is a bitwise
!> XOR -- one `ieor` per occupied orbital, over a determinant that is already
!> a bit pattern. Classifying the whole determinant basis is therefore
!> O(ndet * nact) integer work, negligible beside the sigma builds.
!>
!> Why this is worth having is not speed. The CI solver returns the lowest
!> roots of the active-space Hamiltonian, and those roots need not be the
!> state the user asked for: a CASSCF that starts in one irrep can converge
!> to another, and the run reports a converged energy for the wrong state
!> with no diagnostic. The spin channel already handles the analogous problem
!> -- solve extra roots, classify by <S^2>, keep the ones with the requested
!> multiplicity (`select_by_multiplicity`) -- and this is the same filter for
!> the spatial quantum number, deliberately shaped the same way.
!>
!> Nothing here changes a result on its own; a caller that does not ask for a
!> target irrep never reaches it.
module fci_symmetry

  use precision, only: dp
  use iso_fortran_env, only: int64

  implicit none

  private
  public :: fci_det_irrep_code
  public :: fci_det_irrep_codes
  public :: fci_classify_root_irreps
  public :: fci_select_by_irrep
  public :: fci_code_to_index

  integer, parameter :: i8 = int64

contains

!###############################################################################

!> @brief XOR code of one determinant.
!>
!> @param[in] nact      active orbitals
!> @param[in] orb_code  XOR code of each active orbital, 1-based, [nact]
!> @param[in] det       determinant bit pattern: alpha string in bits
!>                      0..nact-1, beta string in bits nact..2*nact-1, which
!>                      is the layout fci_sigma_strings builds and reads
!>
!> Both spin strings contribute: a beta electron carries the same spatial
!> irrep as an alpha electron in the same orbital.
  pure integer function fci_det_irrep_code(nact, orb_code, det) result(code)
    integer, intent(in) :: nact
    integer, intent(in) :: orb_code(:)
    integer(i8), intent(in) :: det

    integer :: p

    code = 0
    do p = 1, nact
      ! bit p-1 = alpha occupation, bit nact+p-1 = beta occupation
      if (btest(det, p-1)) code = ieor(code, orb_code(p))
      if (btest(det, nact+p-1)) code = ieor(code, orb_code(p))
    end do
  end function fci_det_irrep_code

!###############################################################################

!> @brief XOR code of every determinant in the basis.
  subroutine fci_det_irrep_codes(nact, orb_code, ndet, dets, codes)
    integer, intent(in) :: nact
    integer, intent(in) :: orb_code(:)
    integer(i8), intent(in) :: ndet
    integer(i8), intent(in) :: dets(*)
    integer, intent(out) :: codes(*)

    integer(i8) :: k

    !$omp parallel do default(shared) private(k) schedule(static)
    do k = 1_i8, ndet
      codes(k) = fci_det_irrep_code(nact, orb_code, dets(k))
    end do
    !$omp end parallel do
  end subroutine fci_det_irrep_codes

!###############################################################################

!> @brief Map an XOR code back to a 1-based irrep index. 0 if unknown.
  pure integer function fci_code_to_index(xor_code, nirrep, code) result(idx)
    integer, intent(in) :: xor_code(:)   !< xor_code(i) for irrep i, 1-based
    integer, intent(in) :: nirrep
    integer, intent(in) :: code

    integer :: i

    idx = 0
    do i = 1, nirrep
      if (xor_code(i) == code) then
        idx = i
        return
      end if
    end do
  end function fci_code_to_index

!###############################################################################

!> @brief Dominant irrep of each CI root, with the weight carried there.
!>
!> @param[in]  nact     active orbitals
!> @param[in]  ndet     determinants
!> @param[in]  dets     determinant bit patterns, [ndet]
!> @param[in]  orb_code XOR code per active orbital, [nact]
!> @param[in]  xor_code XOR code per irrep, 1-based, [nirrep]
!> @param[in]  nirrep   number of irreps
!> @param[in]  nvec     roots to classify
!> @param[in]  civecs   CI vectors, column-major [ndet, nvec]
!> @param[out] irr_win  dominant 1-based irrep index per root, [nvec]
!> @param[out] purity   weight in that irrep per root, [nvec]
!>
!> An exact eigenvector of a Hamiltonian that commutes with the group is pure,
!> so `purity` should come out at 1 to solver accuracy. It is returned rather
!> than asserted because it is the honest diagnostic: a value below 1 means
!> the orbitals are not symmetry-adapted (a 'mixed' label, a broken-symmetry
!> reference, or a geometry that only nearly has the group), and a caller that
!> filters on irrep should decline rather than pick a dominant component out
!> of a genuinely mixed vector.
  subroutine fci_classify_root_irreps(nact, ndet, dets, orb_code, &
                                      xor_code, nirrep, nvec, civecs, &
                                      irr_win, purity)
    integer, intent(in) :: nact, nirrep, nvec
    integer(i8), intent(in) :: ndet
    integer(i8), intent(in) :: dets(*)
    integer, intent(in) :: orb_code(:)
    integer, intent(in) :: xor_code(:)
    real(dp), intent(in) :: civecs(ndet, nvec)
    integer, intent(out) :: irr_win(*)
    real(dp), intent(out) :: purity(*)

    integer, allocatable :: codes(:), idx(:)
    real(dp), allocatable :: weight(:)
    integer(i8) :: k
    integer :: v, i, best
    real(dp) :: total, c

    allocate(codes(ndet), idx(ndet), weight(nirrep))

    call fci_det_irrep_codes(nact, orb_code, ndet, dets, codes)
    do k = 1_i8, ndet
      idx(k) = fci_code_to_index(xor_code, nirrep, codes(k))
    end do

    do v = 1, nvec
      weight = 0.0_dp
      total = 0.0_dp
      do k = 1_i8, ndet
        c = civecs(k, v)*civecs(k, v)
        total = total + c
        i = idx(k)
        ! A determinant whose code matches no irrep cannot be attributed;
        ! it still counts in the total, so it lowers the reported purity
        ! instead of quietly vanishing from it.
        if (i >= 1) weight(i) = weight(i) + c
      end do
      best = maxloc(weight, 1)
      irr_win(v) = best
      if (total > 0.0_dp) then
        purity(v) = weight(best)/total
      else
        purity(v) = 0.0_dp
      end if
    end do

    deallocate(codes, idx, weight)

  end subroutine fci_classify_root_irreps

!###############################################################################

!> @brief Keep the first @p nroot roots whose dominant irrep is @p target.
!>
!> Deliberately the same shape as select_by_multiplicity: a prefix of solved
!> roots comes in, the ones that match come out in energy order, and the
!> caller grows the prefix and retries when too few match.
!>
!> @param[in]  irr_win  dominant irrep index per root, [limit]
!> @param[in]  purity   weight in that irrep per root, [limit]
!> @param[in]  limit    roots classified so far
!> @param[in]  target   1-based irrep index to keep
!> @param[in]  min_pure reject a root whose dominant irrep holds less than
!>                      this fraction -- such a root is not a symmetry
!>                      eigenstate and calling it one would be a guess
!> @param[in]  nroot    how many are wanted
!> @param[out] keep     indices of the kept roots, [nroot]
!> @param[out] nsel     how many were found
  subroutine fci_select_by_irrep(irr_win, purity, limit, target, min_pure, &
                                 nroot, keep, nsel)
    integer, intent(in) :: irr_win(*)
    real(dp), intent(in) :: purity(*)
    integer, intent(in) :: limit, target, nroot
    real(dp), intent(in) :: min_pure
    integer, allocatable, intent(inout) :: keep(:)
    integer, intent(out) :: nsel

    integer :: k

    if (allocated(keep)) deallocate(keep)
    allocate(keep(nroot))
    nsel = 0
    do k = 1, limit
      if (irr_win(k) /= target) cycle
      if (purity(k) < min_pure) cycle
      nsel = nsel + 1
      keep(nsel) = k
      if (nsel == nroot) exit
    end do
  end subroutine fci_select_by_irrep

end module fci_symmetry
