"""Signed bra-side Hamiltonian actions and the public MS-CASPT2 route."""
import weakref

import numpy as np
import pytest

from oqp.library import caspt2_dyall as cp


@pytest.mark.parametrize("ndet", [1, 9, 36])
def test_action_pins_original_signed_contraction(ndet):
    rng = np.random.default_rng(723 + ndet)
    # Non-symmetric H catches a lost transpose; unequal bra/ket catches signs.
    h = rng.normal(size=(ndet, ndet))
    bra = rng.normal(size=ndet)
    ket = rng.normal(size=ndet)
    saved = h.copy(), bra.copy(), ket.copy()
    action = cp._multiset_coupling_action(h, bra)
    np.testing.assert_allclose(action @ ket, bra @ (h @ ket), rtol=0, atol=1e-12)
    np.testing.assert_array_equal(h, saved[0])
    np.testing.assert_array_equal(bra, saved[1])
    np.testing.assert_array_equal(ket, saved[2])
    assert not np.shares_memory(action, h)
    assert not np.shares_memory(action, bra)
    href = weakref.ref(h)
    del h, saved
    assert href() is None
    assert np.isfinite(action).all()


@pytest.mark.parametrize("frozen", [0, 1])
def test_multiset_reaches_action_and_pins_old_cross_terms(tmp_path, monkeypatch, frozen):
    import oqp
    from oqp.pyoqp import Runner

    assert oqp.BACKEND_AVAILABLE, "This regression requires the installed native backend"
    actions, rotations, references = [], [], []
    original_action = cp._multiset_coupling_action
    original_rotate = cp._rotate_det_vector
    original_driver = cp._multistate_multiset

    def action(h, bra):
        # A previous state's dense matrix must already be gone while this
        # driver's per-state records are still alive, not just after return.
        assert all(ref() is None for _, _, ref in actions)
        actions.append((h.copy(), bra.copy(), weakref.ref(h)))
        return original_action(h, bra)

    def rotate(vec, *args):
        if len(actions) == 3:
            # Cross-state rotations start after the last action.  Retaining
            # even the final local H until driver return must fail this test.
            assert all(ref() is None for _, _, ref in actions)
        result = original_rotate(vec, *args)
        rotations.append(result.copy())
        return result

    def driver(*args, **kwargs):
        rotations.clear()
        result = original_driver(*args, **kwargs)
        nstate = len(actions)
        assert nstate == 3
        # First nstate rotations form each reference; the following rotations
        # are the two bra/ket directions for each ordered off-diagonal pair.
        cross_rotations = iter(rotations[nstate:])
        pinned = result["heff"].copy()
        for i in range(nstate):
            for j in range(nstate):
                if i == j:
                    continue
                vji, vij = next(cross_rotations), next(cross_rotations)
                hi, bi, _ = actions[i]
                hj, bj, _ = actions[j]
                pinned[i, j] = 0.5 * (bi @ (hi @ vji) + bj @ (hj @ vij))
        pinned = 0.5 * (pinned + pinned.T)
        np.testing.assert_allclose(result["heff"], pinned, rtol=0, atol=1e-12)
        np.testing.assert_allclose(result["ms_energies"], np.linalg.eigvalsh(pinned),
                                   rtol=0, atol=1e-12)
        # Snapshots are owned by the test; the production matrix itself must
        # no longer be retained after each state's action has been formed.
        assert all(ref() is None for _, _, ref in actions)
        references.extend(result["ms_energies"])
        return result

    monkeypatch.setattr(cp, "_multiset_coupling_action", action)
    monkeypatch.setattr(cp, "_rotate_det_vector", rotate)
    monkeypatch.setattr(cp, "_multistate_multiset", driver)
    runner = Runner(project=f"ms_action_{frozen}", input_file=None,
                    log=str(tmp_path / "ms.log"), silent=1, usempi=False,
                    input_dict={
        "input": {"method": "ms-caspt2", "runtype": "energy", "basis": "sto-3g",
                  "system": "\nH 0 0 0\nH 0 0 0.74\nH 0 0 1.48\nH 0 0 2.22"},
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
    assert len(references) == 3
    np.testing.assert_allclose(runner.mol.energies, references, rtol=0, atol=1e-12)
