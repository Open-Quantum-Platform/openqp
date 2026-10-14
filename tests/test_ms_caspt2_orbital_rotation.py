"""MS-CASPT2 state-to-state coefficients from independent Slater overlaps."""
from itertools import combinations, product
from inspect import signature

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
        inputs = signature(original_driver).bind(*args, **kwargs).arguments
        assert list(inputs["roots"]) == [0, 1, 2]
        state_orbitals.clear()
        cross_calls.clear()
        direction_sensitive.clear()
        result = original_driver(*args, **kwargs)
        assert len(state_orbitals) == 3 and len(cross_calls) == 12
        assert max(direction_sensitive) > 1e-6, "Fixture must distinguish inverse rotations"
        return result

    monkeypatch.setattr(cp, "_rotate_det_vector", rotate)
    monkeypatch.setattr(cp, "_multistate_multiset", driver)
    runner = _asymmetric_h4_runner(tmp_path, frozen)
    runner.run(test_mod=True)
    assert len(cross_calls) == 12


def _asymmetric_h4_runner(tmp_path, frozen):
    from oqp.pyoqp import Runner

    return Runner(project=f"ms_rotation_{frozen}", input_file=None,
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


def _fermion_action(det, operations):
    """Apply creation/annihilation in right-to-left order, including parity."""
    phase = 1
    for orbital, create in operations:
        occupied = bool(det & (1 << orbital))
        if occupied == create:
            return None
        phase *= (-1) ** bin(det & ((1 << orbital) - 1)).count("1")
        det ^= 1 << orbital
    return det, phase


def _common_frame_hamiltonian(hcore_ao, eri_ao, coeff, enuc, nelec):
    """Independent <D|H|D'> in C_ref, with no production integral/FCI helper.

    H = enuc + h[p,q] a_p^+ a_q
        + (1/2) (p r|q s) a_p^+ a_q^+ a_s a_r.
    Alpha spin orbitals precede beta; spatial ERIs use chemists' notation.
    """
    norb = coeff.shape[1]
    h = np.einsum("up,uv,vq->pq", coeff, hcore_ao, coeff, optimize=True)
    eri = np.einsum("up,vq,wr,xs,uvwx->pqrs", coeff, coeff, coeff, coeff,
                    eri_ao, optimize=True)
    alpha = combinations(range(norb), nelec[0])
    beta = list(combinations(range(norb, 2*norb), nelec[1]))
    dets = [sum(1 << p for p in a+b) for a in alpha for b in beta]
    rows = {d: i for i, d in enumerate(dets)}
    ham = enuc * np.eye(len(dets))
    for col, det in enumerate(dets):
        for p, q in product(range(2*norb), repeat=2):
            if p // norb != q // norb:
                continue
            action = _fermion_action(det, [(q, False), (p, True)])
            if action is not None:
                target, phase = action
                ham[rows[target], col] += phase * h[p % norb, q % norb]
        for p, q, r, t in product(range(2*norb), repeat=4):
            if p // norb != r // norb or q // norb != t // norb:
                continue
            action = _fermion_action(det, [(r, False), (t, False),
                                           (q, True), (p, True)])
            if action is not None:
                target, phase = action
                ham[rows[target], col] += 0.5 * phase * eri[
                    p % norb, r % norb, q % norb, t % norb]
    np.testing.assert_allclose(ham, ham.T, rtol=0, atol=1e-12)
    return ham


@pytest.mark.parametrize("frozen", [0, 1])
def test_multiset_final_energies_and_states_match_common_frame_reference(
        tmp_path, monkeypatch, frozen):
    import oqp

    assert oqp.BACKEND_AVAILABLE, "The installed native OpenQP backend is required"
    probe = cp._pt2_lib()
    assert probe is not None and hasattr(probe[0], "pt2_minor_transform")
    original_driver = cp._multistate_multiset
    original_rotate = cp._rotate_det_vector
    original_lift = cp._lift_to_full
    expected = {}

    def driver(*args, **kwargs):
        inputs = signature(original_driver).bind(*args, **kwargs).arguments
        roots = list(inputs["roots"])
        assert roots == [0, 1, 2]
        norb = inputs["nbf"]
        nelec = tuple(inputs["ncore"] + n for n in inputs["active_nelec"])
        state_orbitals, refs_common, psi_local = [], [], []
        saved_inputs = {key: inputs[key].copy() for key in
                        ("hcore_ao", "eri_ao", "coeff_ref")}

        def rotate(vec, rotation, *rest):
            if len(state_orbitals) < len(roots):
                # These are upstream inputs BEFORE inter-state contractions.
                state_orbitals.append(rotation.copy())
                refs_common.append(vec.copy())
            return original_rotate(vec, rotation, *rest)

        def lift(*lift_args, **lift_kwargs):
            vec = original_lift(*lift_args, **lift_kwargs)
            psi_local.append(vec.copy())
            return vec

        with monkeypatch.context() as capture:
            capture.setattr(cp, "_rotate_det_vector", rotate)
            capture.setattr(cp, "_lift_to_full", lift)
            result = original_driver(*args, **kwargs)
        assert len(state_orbitals) == len(psi_local) == len(roots)
        for key, saved in saved_inputs.items():
            np.testing.assert_array_equal(inputs[key], saved)
        identity = np.eye(norb)
        to_common = [_slater_overlap(identity, u, nelec) for u in state_orbitals]
        for u in state_orbitals:
            np.testing.assert_allclose(u.T @ u, identity, rtol=0, atol=1e-12)
        refs = np.column_stack(refs_common)
        np.testing.assert_allclose(refs.T @ refs, np.eye(len(roots)), rtol=0, atol=1e-12)
        psi = np.column_stack([t @ v for t, v in zip(to_common, psi_local)])
        for state, local in enumerate(psi_local):
            np.testing.assert_allclose(np.linalg.norm(psi[:, state]), np.linalg.norm(local),
                                       rtol=0, atol=1e-12)
        np.testing.assert_allclose(refs.T @ psi, 0, rtol=0, atol=1e-12)
        # Use actual prepared AO/MO inputs, not result['heff'], energies, e2,
        # mixing, production minor matrices or per-state Hamiltonian outputs.
        ham = _common_frame_hamiltonian(inputs["hcore_ao"], inputs["eri_ao"],
                                        inputs["coeff_ref"], inputs["enuc"], nelec)
        signed = psi.T @ ham @ refs
        ref_energies = np.diag(refs.T @ ham @ refs).copy()
        corrections = np.diag(signed).copy()
        heff = 0.5 * (signed + signed.T) + np.diag(ref_energies)
        energies, mixing = np.linalg.eigh(heff)
        assert np.min(np.diff(energies)) > 1e-4, "Fixture needs distinct state projectors"

        # Independently evaluate the old inverse-direction coupling on EXACTLY
        # the same upstream tuple. Sensitivity is required at the final-energy
        # boundary, not merely at a rotation coefficient or call count.
        wrong = heff.copy()
        for i, j in combinations(range(len(roots)), 2):
            local_i = _slater_overlap(state_orbitals[i], identity, nelec) @ refs[:, i]
            local_j = _slater_overlap(state_orbitals[j], identity, nelec) @ refs[:, j]
            wrong_j = to_common[i] @ (
                _slater_overlap(state_orbitals[j], state_orbitals[i], nelec) @ local_j)
            wrong_i = to_common[j] @ (
                _slater_overlap(state_orbitals[i], state_orbitals[j], nelec) @ local_i)
            wrong[i, j] = wrong[j, i] = 0.5 * (
                psi[:, i] @ ham @ wrong_j + psi[:, j] @ ham @ wrong_i)
        assert np.max(np.abs(heff - wrong)) > 1e-8, "Fixture needs signed Heff sensitivity"
        assert np.max(np.abs(energies - np.linalg.eigvalsh(wrong))) > 1e-8, \
            "Fixture must fail for the old final energies"
        np.testing.assert_allclose(result["heff"], heff, rtol=0, atol=1e-11)
        np.testing.assert_allclose(result["ms_energies"], energies, rtol=0, atol=1e-10)
        np.testing.assert_allclose(result["ref_energies"], ref_energies, rtol=0, atol=1e-11)
        np.testing.assert_allclose(result["e2"], corrections, rtol=0, atol=1e-11)
        np.testing.assert_allclose(result["ss_energies"], np.diag(heff), rtol=0, atol=1e-11)
        actual_mixing = result["mixing"]
        np.testing.assert_allclose(actual_mixing.T @ actual_mixing, np.eye(3),
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(heff @ actual_mixing, actual_mixing * energies,
                                   rtol=0, atol=1e-11)
        # Each projector is invariant to the arbitrary eigenvector sign, while
        # signed Heff above retains the phases of the actual reference roots.
        for state in range(len(roots)):
            actual = result["mixing"][:, state]
            reference = mixing[:, state]
            np.testing.assert_allclose(np.outer(actual, actual),
                                       np.outer(reference, reference), rtol=0, atol=1e-9)
        expected.update(heff=heff, energies=energies, mixing=mixing, refs=refs,
                        ref_energies=ref_energies, corrections=corrections)
        return result

    monkeypatch.setattr(cp, "_multistate_multiset", driver)
    runner = _asymmetric_h4_runner(tmp_path, frozen)
    runner.run(test_mod=True)
    assert expected, "The actual molecular multiset driver must execute"
    np.testing.assert_allclose(runner.mol.energies, expected["energies"], rtol=0, atol=1e-10)
    np.testing.assert_allclose(runner.mol.mol_energy.energy, expected["energies"][0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(runner.mol.data["OQP::CASPT2_ENERGIES"],
                               expected["energies"], rtol=0, atol=1e-10)
    np.testing.assert_allclose(runner.mol.data["OQP::CASPT2_EFFECTIVE_HAMILTONIAN"],
                               expected["heff"], rtol=0, atol=1e-11)
    for key, reference in [
            ("REFERENCE_ENERGIES", expected["ref_energies"]),
            ("STATE_SPECIFIC_CORRECTIONS", expected["corrections"]),
            ("SS_ENERGIES", np.diag(expected["heff"]))]:
        np.testing.assert_allclose(runner.mol.data["OQP::CASPT2_" + key], reference,
                                   rtol=0, atol=1e-11)
    np.testing.assert_array_equal(runner.mol.data["OQP::CASPT2_TARGET_ROOTS"], [0, 1, 2])
    mixing = np.asarray(runner.mol.data["OQP::CASPT2_MIXING_VECTORS"])
    for state in range(3):
        actual = expected["refs"] @ mixing[:, state]
        reference = expected["refs"] @ expected["mixing"][:, state]
        np.testing.assert_allclose(np.outer(actual, actual), np.outer(reference, reference),
                                   rtol=0, atol=1e-9)
