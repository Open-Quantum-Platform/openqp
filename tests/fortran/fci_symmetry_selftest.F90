module fci_symmetry_selftest_mod
!> @brief Check determinant/root irrep classification against hand-built cases.
!>
!>   The properties worth pinning are the ones a wrong implementation would
!>   still look plausible under:
!>
!>     * the closed-shell reference determinant is totally symmetric, whatever
!>       the orbital irreps are -- every orbital enters twice, so the XOR
!>       cancels. A version that folded only the alpha string would pass a
!>       symmetric test case and fail here.
!>     * a single alpha excitation carries the product of the two orbitals it
!>       connects, with the C2v products spelled out.
!>     * a CI vector confined to one irrep classifies as that irrep with
!>       purity 1; a deliberate 50/50 mixture classifies with purity 0.5, so
!>       a caller filtering on irrep can tell the two apart.
!>     * selection returns matching roots in energy order and skips impure
!>       ones.
!>
!>   Results are written to /tmp/fci_symmetry_selftest.out.

  implicit none
  character(len=*), parameter :: module_name = "fci_symmetry_selftest_mod"

contains

  subroutine fci_symmetry_selftest_C() bind(C, name="fci_symmetry_selftest")
    call fci_symmetry_selftest()
  end subroutine fci_symmetry_selftest_C

  subroutine fci_symmetry_selftest()

    use precision, only: dp
    use iso_fortran_env, only: int64
    use fci_symmetry, only: fci_det_irrep_code, fci_classify_root_irreps, &
                            fci_select_by_irrep, fci_code_to_index

    implicit none

    integer, parameter :: i8 = int64
    ! C2v: a1 a2 b1 b2 -> XOR codes over (C2, sigma_v, sigma_v')
    integer, parameter :: NIRR = 4
    integer :: xor_code(NIRR)
    integer, parameter :: A1 = 1, A2 = 2, B1 = 3, B2 = 4

    integer, parameter :: NACT = 4
    integer :: orb_irrep(NACT), orb_code(NACT)
    integer(i8) :: det_ref, det_ex, dets(4)
    integer :: code, idx, p, iu, nsel
    integer :: irr_win(3)
    real(dp) :: purity(3), civecs(4, 3)
    integer, allocatable :: keep(:)
    logical :: ok_ref, ok_ex, ok_pure, ok_mix, ok_sel

    ! Characters over 3 non-identity ops; code = bitmask of chi = -1.
    xor_code(A1) = 0                       ! (+ + +)
    xor_code(A2) = 1 + 2                   ! (+ - -)
    xor_code(B1) = 1 + 4                   ! (- + -)  (bit0, bit2)
    xor_code(B2) = 2 + 4                   ! (- - +)  (bit1, bit2)

    ! Active orbitals: a1 b2 a1 b1
    orb_irrep = [A1, B2, A1, B1]
    do p = 1, NACT
      orb_code(p) = xor_code(orb_irrep(p))
    end do

    ! ---- closed-shell reference: orbitals 1,2 doubly occupied -------------
    ! alpha bits 0,1 ; beta bits NACT+0, NACT+1
    det_ref = 0_i8
    det_ref = ibset(det_ref, 0)
    det_ref = ibset(det_ref, 1)
    det_ref = ibset(det_ref, NACT + 0)
    det_ref = ibset(det_ref, NACT + 1)
    code = fci_det_irrep_code(NACT, orb_code, det_ref)
    ok_ref = (fci_code_to_index(xor_code, NIRR, code) == A1)

    ! ---- single alpha excitation 2 -> 4, i.e. b2 -> b1 --------------------
    ! b2 x b1 = a2 in C2v; the beta pair still cancels.
    det_ex = det_ref
    det_ex = ibclr(det_ex, 1)
    det_ex = ibset(det_ex, 3)
    code = fci_det_irrep_code(NACT, orb_code, det_ex)
    ok_ex = (fci_code_to_index(xor_code, NIRR, code) == A2)

    ! ---- root classification ---------------------------------------------
    ! Four determinants: two a1, two a2, by construction above.
    dets(1) = det_ref
    dets(2) = det_ref
    dets(3) = det_ex
    dets(4) = det_ex

    civecs = 0.0_dp
    civecs(1, 1) = 1.0_dp                       ! pure a1
    civecs(3, 2) = 1.0_dp                       ! pure a2
    civecs(1, 3) = sqrt(0.5_dp)                 ! 50/50 mixture
    civecs(3, 3) = sqrt(0.5_dp)

    call fci_classify_root_irreps(NACT, 4_i8, dets, orb_code, xor_code, &
                                  NIRR, 3, civecs, irr_win, purity)

    ok_pure = (irr_win(1) == A1 .and. abs(purity(1) - 1.0_dp) < 1.0e-12_dp) &
        .and. (irr_win(2) == A2 .and. abs(purity(2) - 1.0_dp) < 1.0e-12_dp)
    ok_mix = abs(purity(3) - 0.5_dp) < 1.0e-12_dp

    ! ---- selection --------------------------------------------------------
    ! Asking for a2 with a purity floor must take root 2 and reject root 3.
    call fci_select_by_irrep(irr_win, purity, 3, A2, 0.9_dp, 1, keep, nsel)
    ok_sel = (nsel == 1)
    if (ok_sel) ok_sel = (keep(1) == 2)

    open(newunit=iu, file='/tmp/fci_symmetry_selftest.out', status='replace')
    write(iu,'(A,L1)') 'closed_shell_is_totally_symmetric = ', ok_ref
    write(iu,'(A,L1)') 'b2_to_b1_excitation_is_a2         = ', ok_ex
    write(iu,'(A,L1)') 'pure_roots_classify_with_purity_1 = ', ok_pure
    write(iu,'(A,L1)') 'mixed_root_reports_purity_half    = ', ok_mix
    write(iu,'(A,L1)') 'selection_skips_impure_root       = ', ok_sel
    write(iu,'(A,F6.3)') 'purity_mixed                    = ', purity(3)
    write(iu,'(A,L1)') 'ALL_PASS                          = ', &
        (ok_ref .and. ok_ex .and. ok_pure .and. ok_mix .and. ok_sel)
    close(iu)

    if (allocated(keep)) deallocate(keep)

  end subroutine fci_symmetry_selftest

end module fci_symmetry_selftest_mod
