"""S^2 resolution inside a degenerate CI cluster (issue #339).

A CI root's multiplicity is read from round(sqrt(1 + 4<S^2>)), which identifies
the spin only when the eigenvector IS an S^2 eigenstate.  Inside a degenerate
energy subspace the eigensolver may return any orthogonal mixture, and a
mixture of different-spin states is not an S^2 eigenstate at all.

These tests use the smallest system where that happens: two electrons in two
orbitals at Ms = 0.  Determinant keys put alpha in bits 0..norb-1 and beta in
bits norb..2*norb-1, which is the layout fci_sigma_strings builds and reads.
"""

import numpy as np
import pytest

import oqp.library.fci as fci_module
from oqp.library.fci import (
    ActiveSpacePlan,
    FCISettings,
    compute_s2,
    degenerate_clusters,
    s2_matrix,
    solve_active_ci,
    spin_label_diagnosis,
    spin_purify_degenerate_clusters,
)

NORB = 2
NELEC = (1, 1)
#  |0a 0b>   |0a 1b>   |1a 0b>   |1a 1b>
DETS = [5, 9, 6, 10]
OPEN_A, OPEN_B = 1, 2          # positions of |0a 1b> and |1a 0b> in DETS


def _unit(position):
    vector = np.zeros(len(DETS))
    vector[position] = 1.0
    return vector


def _open_shell_block():
    return np.column_stack([_unit(OPEN_A), _unit(OPEN_B)])


def test_bare_open_shell_determinant_reports_an_impossible_multiplicity():
    """The failure the issue describes, before anything is done about it."""
    for position in (OPEN_A, OPEN_B):
        s2 = compute_s2(_unit(position), DETS, NORB, NELEC)
        assert s2 == pytest.approx(1.0)
        multiplicity = int(round((1.0 + 4.0 * s2) ** 0.5))
        assert multiplicity == 2          # a "doublet" from two electrons
        problems = spin_label_diagnosis([s2], [multiplicity], NELEC)
        assert problems, "the impossible label must be flagged"
        assert "parity" in problems[0][3]


def test_s2_matrix_exposes_the_off_diagonal_the_scalar_form_cannot():
    """<a|S^2|b> is what resolves the subspace; <a|S^2|a> alone cannot."""
    matrix = s2_matrix(_open_shell_block(), DETS, NORB, NELEC)
    assert matrix == pytest.approx(np.array([[1.0, -1.0], [-1.0, 1.0]]))


def test_degenerate_cluster_is_resolved_into_spin_pure_roots():
    energies = np.array([-1.0, -1.0])
    coeffs, s2, multiplicity, changed = spin_purify_degenerate_clusters(
        energies, _open_shell_block(), DETS, NORB, NELEC)

    assert changed
    assert sorted(np.round(s2, 10)) == [0.0, 2.0]
    assert sorted(multiplicity.tolist()) == [1, 3]
    assert not spin_label_diagnosis(s2, multiplicity, NELEC)

    # The textbook combinations, up to the sign convention below.
    root = 1.0 / np.sqrt(2.0)
    for column in range(2):
        weights = np.abs(coeffs[:, column])
        assert weights[OPEN_A] == pytest.approx(root)
        assert weights[OPEN_B] == pytest.approx(root)
        assert weights[0] == pytest.approx(0.0, abs=1e-12)
        assert weights[3] == pytest.approx(0.0, abs=1e-12)

    # Ascending S^2 inside the cluster, so the order does not depend on what
    # the eigensolver happened to produce.
    assert s2[0] < s2[1]


def test_sign_convention_makes_the_rotation_reproducible():
    """The largest-magnitude coefficient of every rotated root is positive."""
    coeffs, _s2, _mult, changed = spin_purify_degenerate_clusters(
        np.array([-1.0, -1.0]), _open_shell_block(), DETS, NORB, NELEC)
    assert changed
    for column in range(coeffs.shape[1]):
        col = coeffs[:, column]
        assert col[int(np.argmax(np.abs(col)))] > 0.0


def test_non_degenerate_roots_are_left_alone():
    block = _open_shell_block()
    coeffs, _s2, _mult, changed = spin_purify_degenerate_clusters(
        np.array([-1.0, -0.5]), block, DETS, NORB, NELEC)
    assert not changed
    assert coeffs == pytest.approx(block)


def test_already_pure_degenerate_cluster_is_not_rotated():
    """A degenerate pair the solver already returned spin-pure stays put.

    This is the common case -- every committed WF example that reaches a
    degenerate cluster is in it -- and rotating there would churn references
    for nothing.
    """
    root = 1.0 / np.sqrt(2.0)
    singlet = root * (_unit(OPEN_A) + _unit(OPEN_B))
    triplet = root * (_unit(OPEN_A) - _unit(OPEN_B))
    block = np.column_stack([singlet, triplet])
    coeffs, s2, _mult, changed = spin_purify_degenerate_clusters(
        np.array([-1.0, -1.0]), block, DETS, NORB, NELEC)
    assert not changed
    assert coeffs == pytest.approx(block)
    assert sorted(np.round(s2, 10)) == [0.0, 2.0]


def test_native_engine_results_are_purified_too():
    """solve_active_ci purifies the NATIVE result, not just the Python fallback.

    The production path returns from _lib_fci_solve; if purification ran only
    after the solve_fci fallback, the compiled engine would still publish
    spin-mixed roots (the first Codex P1 on PR #367).  The engine is faked
    here to return the issue's mixed pair, in the determinant order
    _determinants produces -- the same order the real engine is pinned to by
    tests/test_fci_solve.py.
    """
    plan = ActiveSpacePlan(norb=NORB, active=(0, 1), core=(), nelec=NELEC,
                           metadata={})
    settings = FCISettings(nroot=2, active_electrons=NELEC,
                           active_orbitals=NORB)
    mixed = _open_shell_block()

    def fake_native(h1e, eri, plan_, spec, *, nthreads, want_s2,
                    use_target_spin):
        return (np.array([-2.0, -2.0]), np.array(mixed, copy=True),
                np.array([1.0, 1.0]) if want_s2 else None)

    original = fci_module._lib_fci_solve
    fci_module._lib_fci_solve = fake_native
    try:
        energies, coeffs, s2 = solve_active_ci(
            -np.eye(NORB), np.zeros((NORB,) * 4), plan, 0.0, settings,
            nroot=2, want_s2=True)
    finally:
        fci_module._lib_fci_solve = original

    assert energies == pytest.approx([-2.0, -2.0])
    assert sorted(np.round(s2, 10)) == [0.0, 2.0]
    assert not spin_label_diagnosis(
        s2, np.rint(np.sqrt(1.0 + 4.0 * np.asarray(s2))).astype(int), NELEC)
    root = 1.0 / np.sqrt(2.0)
    for column in range(2):
        weights = np.abs(coeffs[:, column])
        assert weights[OPEN_A] == pytest.approx(root)
        assert weights[OPEN_B] == pytest.approx(root)


def test_ambiguous_root_beyond_the_selection_does_not_refuse():
    """Refuse only for labels that can change WHICH roots come back.

    Selection stops once the request is filled, so an untrustworthy label on a
    higher inspected root cannot affect the answer; refusing on it turned an
    unambiguous request into an error (the second Codex P2 on PR #367).
    """
    four = (2, 2)                             # singlet/triplet/quintet allowed
    energies = np.array([-1.0, -0.5, -0.2])
    coeffs = np.eye(3)
    s2 = np.array([0.0, 1.0, 2.0])            # root 1 carries the bad label
    mult = np.array([1, 2, 3])                # "doublet" is impossible here

    # Root 1 is untrustworthy but sits past the singlet the request is filled
    # by, so it cannot change the answer.
    kept_e, _c, kept_s2, kept_m, idx = fci_module._filter_roots_by_target_spin(
        energies, coeffs, s2, mult,
        target_spin="singlet", requested_nroot=1,
        ci_label="FCI", ci_section="[fci]", nelec=four)
    assert kept_e == pytest.approx([-1.0])
    assert kept_s2 == pytest.approx([0.0])
    assert kept_m.tolist() == [1] and idx.tolist() == [0]

    # The same bad label now lies inside the span the selection walks.
    with pytest.raises(fci_module.SpinLabelAmbiguityError):
        fci_module._filter_roots_by_target_spin(
            energies, coeffs, s2, mult,
            target_spin="triplet", requested_nroot=1,
            ci_label="FCI", ci_section="[fci]", nelec=four)

    # And when the request cannot be filled, every inspected root could have
    # supplied a missing match, so the ambiguity still decides it.
    with pytest.raises(fci_module.SpinLabelAmbiguityError):
        fci_module._filter_roots_by_target_spin(
            energies, coeffs, s2, mult,
            target_spin="singlet", requested_nroot=2,
            ci_label="FCI", ci_section="[fci]", nelec=four)


def test_target_the_electron_count_forbids_is_rejected_as_such():
    """An impossible target is answered as impossible, not as bad labels.

    Two electrons have no quintet; refusing that request because some root's
    <S^2> is untrustworthy names the wrong problem, and no amount of extra
    roots would ever satisfy it.
    """
    with pytest.raises(ValueError, match="found no matching roots") as raised:
        fci_module._filter_roots_by_target_spin(
            np.array([-1.0, -1.0]), _open_shell_block(),
            np.array([1.0, 1.0]), np.array([2, 2]),
            target_spin="quintet", requested_nroot=1,
            ci_label="FCI", ci_section="[fci]", nelec=NELEC)
    assert not isinstance(raised.value, fci_module.SpinLabelAmbiguityError)


def test_cluster_grouping_uses_the_tolerance():
    assert degenerate_clusters([-1.0, -1.0 + 1e-12, -0.5]) == [[0, 1], [2]]
    assert degenerate_clusters([-1.0, -0.9, -0.5]) == [[0], [1], [2]]
    assert degenerate_clusters([]) == []
