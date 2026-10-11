"""MS-CASPT2 state-to-state coefficients from independent Slater overlaps."""
from itertools import combinations

import numpy as np
import pytest

from oqp.library import caspt2_dyall as cp


def _slater_overlap(dst, src, nelec):
    """<D_dst|D_src> in a common orthonormal orbital frame, alpha-major."""
    alpha = list(combinations(range(src.shape[1]), nelec[0]))
    beta = list(combinations(range(src.shape[1]), nelec[1]))
    dets = [(a, b) for a in alpha for b in beta]
    overlap = np.empty((len(dets), len(dets)))
    for row, (a_dst, b_dst) in enumerate(dets):
        for col, (a_src, b_src) in enumerate(dets):
            # Derive each spin overlap directly from occupied orbital columns;
            # do not reuse the production minor or rotation helper.
            oa = dst[:, a_dst].T @ src[:, a_src]
            ob = dst[:, b_dst].T @ src[:, b_src]
            overlap[row, col] = np.linalg.det(oa) * np.linalg.det(ob)
    return overlap


@pytest.mark.parametrize("nelec", [(1, 0), (2, 1), (2, 2)])
@pytest.mark.parametrize("backend", ["python", "native"])
def test_rotate_coefficients_follow_orbital_columns(monkeypatch, nelec, backend):
    if backend == "python":
        monkeypatch.setattr(cp, "_pt2_lib", lambda: None)
    else:
        probe = cp._pt2_lib()
        assert probe is not None and hasattr(probe[0], "pt2_minor_transform"), \
            "This test requires the installed native minor-transform engine"
    src = np.eye(4)
    src[:2, :2] = [[5/13, -12/13], [12/13, 5/13]]
    dst = np.eye(4)
    dst[:2, :2] = [[3/5, -4/5], [4/5, 3/5]]
    dst[2, 2] = -1
    rotation = src.T @ dst
    assert np.max(np.abs(rotation - rotation.T)) > 0.1
    expected_transform = _slater_overlap(dst, src, nelec)
    vec = (-1.0) ** np.arange(len(expected_transform)) * \
        np.arange(1, len(expected_transform)+1, dtype=float)
    vec /= np.linalg.norm(vec)
    saved = vec.copy(), rotation.copy()
    result = cp._rotate_det_vector(vec, rotation, 4, nelec)
    np.testing.assert_allclose(result, expected_transform @ vec, rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        cp._rotate_det_vector(result, rotation.T, 4, nelec), vec, rtol=0, atol=1e-12)
    np.testing.assert_array_equal(vec, saved[0])
    np.testing.assert_array_equal(rotation, saved[1])


@pytest.mark.parametrize("frozen", [0, 1])
def test_multiset_cross_rotations_match_independent_slater_overlaps(
        tmp_path, monkeypatch, frozen):
    import oqp
    from oqp.pyoqp import Runner

    assert oqp.BACKEND_AVAILABLE, "The installed native OpenQP backend is required"
    probe = cp._pt2_lib()
    assert probe is not None and hasattr(probe[0], "pt2_minor_transform")
    original_rotate = cp._rotate_det_vector
    original_driver = cp._multistate_multiset
    state_orbitals = []
    cross_calls = []
    direction_sensitive = []
    pairs = [(i, j) for i in range(3) for j in range(3) if i != j]

    def rotate(vec, rotation, norb, nelec):
        saved = vec.copy(), rotation.copy()
        result = original_rotate(vec, rotation, norb, nelec)
        np.testing.assert_array_equal(vec, saved[0])
        np.testing.assert_array_equal(rotation, saved[1])
        if len(state_orbitals) < 3:
            # Intra-state calls transform C_ref -> C_ref U_I, so the passed
            # rotation supplies each state's orbitals in the common MO frame.
            state_orbitals.append(rotation.copy())
        else:
            index = len(cross_calls)
            i, j = pairs[index // 2]
            source, destination = (j, i) if index % 2 == 0 else (i, j)
            src, dst = state_orbitals[source], state_orbitals[destination]
            np.testing.assert_allclose(src.T @ src, np.eye(norb), rtol=0, atol=1e-12)
            np.testing.assert_allclose(dst.T @ dst, np.eye(norb), rtol=0, atol=1e-12)
            expected = _slater_overlap(dst, src, nelec) @ vec
            inverse = _slater_overlap(src, dst, nelec) @ vec
            direction_sensitive.append(float(np.max(np.abs(expected - inverse))))
            # Signed coefficients, not an energy or a saved producer output,
            # distinguish J->I from I->J even if spectra happen to agree.
            np.testing.assert_allclose(result, expected, rtol=0, atol=1e-12)
            cross_calls.append((source, destination))
        return result

    def driver(*args, **kwargs):
        assert list(args[10]) == [0, 1, 2]
        state_orbitals.clear()
        cross_calls.clear()
        direction_sensitive.clear()
        result = original_driver(*args, **kwargs)
        assert len(state_orbitals) == 3 and len(cross_calls) == 12
        assert max(direction_sensitive) > 1e-6, "Fixture must distinguish inverse rotations"
        return result

    monkeypatch.setattr(cp, "_rotate_det_vector", rotate)
    monkeypatch.setattr(cp, "_multistate_multiset", driver)
    runner = Runner(project=f"ms_rotation_{frozen}", input_file=None,
                    log=str(tmp_path / "ms.log"), silent=1, usempi=False,
                    input_dict={
        "input": {"method": "ms-caspt2", "runtype": "energy", "basis": "sto-3g",
                  "system": "\nH 0 0 0\nH 0 0 0.74\nH 0.13 0.07 1.55\nH -0.08 0.02 2.27"},
        "guess": {"type": "hcore"},
        "scf": {"type": "rhf", "multiplicity": "1", "conv": "1e-10",
                "maxit": "80", "forced_attempt": "3", "save_molden": "false"},
        "properties": {"scf_prop": ""},
        "cas": {"active_electrons": "2", "active_orbitals": "2",
                "frozen_core": "1", "max_det": "2000"},
        "ci": {"nroot": "3", "target_spin": "singlet", "solver": "dense",
               "eig_tol": "1e-10", "integral_backend": "native"},
        "pt2": {"reference": "casci", "target_roots": "0,1,2", "nroot": "3",
                "ipea_shift": "0", "frozen": str(frozen)},
        "tests": {"exception": "true"},
    })
    runner.run(test_mod=True)
    assert len(cross_calls) == 12
