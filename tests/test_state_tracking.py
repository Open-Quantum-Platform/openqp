"""Regressions for global MO/root assignment and transported phase metadata."""

from __future__ import annotations

import numpy as np
import pytest


oqp = pytest.importorskip("oqp")

from oqp.library.state_tracking import (
    diagonal_phase_tracking,
    maximum_overlap_assignment,
)
from oqp.library import single_point
from oqp.molecule.molecule import Molecule
from oqp.pyoqp import Runner


def test_global_assignment_beats_row_greedy_counterexample():
    overlap = np.array([[0.90, 0.80], [0.85, 0.10]])

    order, signs, matched, margins = maximum_overlap_assignment(overlap)

    # Row-greedy chooses [0, 1] with total 1.00.  The globally optimal
    # one-to-one state map is [1, 0] with total 1.65.
    assert order.tolist() == [1, 0]
    assert signs.tolist() == [1.0, 1.0]
    assert np.sum(matched) == pytest.approx(1.65)
    assert margins.tolist() == pytest.approx([-0.10, 0.75])


def test_assignment_tracks_multi_root_exchange_and_phase():
    overlap = np.array([
        [0.02, -0.97, 0.01],
        [0.03, 0.01, 0.96],
        [-0.95, 0.02, 0.04],
    ])

    order, signs, matched, margins = maximum_overlap_assignment(overlap)

    assert order.tolist() == [1, 2, 0]
    assert signs.tolist() == [-1.0, 1.0, -1.0]
    assert matched.tolist() == pytest.approx([0.97, 0.96, 0.95])
    assert np.all(margins > 0.90)


def test_zero_selected_overlap_has_deterministic_positive_phase():
    order, signs, matched, _ = maximum_overlap_assignment(np.zeros((2, 2)))

    assert order.tolist() == [0, 1]
    assert signs.tolist() == [1.0, 1.0]
    assert matched.tolist() == [0.0, 0.0]


def test_diagonal_phase_preserves_adiabatic_order_at_root_exchange():
    overlap = np.array([
        [0.10, -0.90],
        [0.95, 0.20],
    ])

    signs, matched, margins = diagonal_phase_tracking(overlap)

    # A global assignment would exchange the roots here.  Phase-only NAMD
    # must not borrow that permutation unless it transports every state label.
    assert signs.tolist() == [1.0, 1.0]
    assert matched.tolist() == pytest.approx([0.10, 0.20])
    assert margins.tolist() == pytest.approx([-0.80, -0.75])


def test_phase_mode_mo_tracking_preserves_orbital_order_at_exchange():
    tracker = single_point.BasisOverlap.__new__(single_point.BasisOverlap)
    tracker.align_type = "phase"
    overlap = np.array([
        [0.10, -0.90],
        [0.95, 0.20],
    ])

    order, signs, matched, margins = tracker.find_mo_order(overlap, diagnostics=True)

    assert order.tolist() == [0, 1]
    assert signs.tolist() == [1.0, 1.0]
    assert matched.tolist() == pytest.approx([0.10, 0.20])
    assert margins.tolist() == pytest.approx([-0.80, -0.75])


def test_namd_align_x_preserves_energy_root_order_and_records_phase(monkeypatch):
    class DummyMol:
        def __init__(self):
            self.data = {
                "OQP::td_bvec_mo_old": np.eye(3),
                "OQP::td_bvec_mo": np.array([
                    [0.0, -1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [1.0, 0.0, 0.0],
                ]),
                "OQP::state_tracking_lineage": np.array([10, 11, 12]),
                "OQP::state_tracking_phase_initial": np.array([1.0, -1.0, 1.0]),
            }

    monkeypatch.setattr(single_point, "dump_log", lambda *_args, **_kwargs: None)
    tracker = single_point.NACME.__new__(single_point.NACME)
    tracker.mol = DummyMol()

    tracker.align_x()

    data = tracker.mol.data
    assert tracker.mol._state_tracking_fresh is True
    assert data["OQP::state_tracking_order"].tolist() == [0, 1, 2]
    assert data["OQP::state_tracking_lineage"].tolist() == [10, 11, 12]
    assert data["OQP::state_tracking_phase_step"].tolist() == [1.0, 1.0, 1.0]
    assert data["OQP::state_tracking_phase_initial"].tolist() == [1.0, 1.0, 1.0]
    assert data["OQP::state_tracking_previous_phase_initial"].tolist() == [1.0, -1.0, 1.0]
    assert np.allclose(data["OQP::td_bvec_mo"], np.array([
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
    ]))


def test_numerical_worker_reorders_exchanged_roots_to_central_order(monkeypatch):
    class DummyMol:
        def __init__(self):
            self.data = {
                "OQP::td_bvec_mo_old": np.eye(3),
                "OQP::td_bvec_mo": np.array([
                    [0.0, -1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [1.0, 0.0, 0.0],
                ]),
                "OQP::state_tracking_lineage": np.array([10, 11, 12]),
                "OQP::state_tracking_phase_initial": np.ones(3),
            }

    monkeypatch.setattr(single_point, "dump_log", lambda *_args, **_kwargs: None)
    tracker = single_point.NACME.__new__(single_point.NACME)
    tracker.mol = DummyMol()

    tracker.align_x(reorder=True)

    data = tracker.mol.data
    assert data["OQP::state_tracking_raw_order"].tolist() == [1, 2, 0]
    assert data["OQP::state_tracking_order"].tolist() == [0, 1, 2]
    assert data["OQP::state_tracking_lineage"].tolist() == [10, 11, 12]
    assert data["OQP::state_tracking_output_reordered"].tolist() == [1]
    assert np.allclose(data["OQP::td_bvec_mo"], np.eye(3))


def test_align_x_tracks_nonsquare_state_major_bridge_buffer(monkeypatch):
    class DummyMol:
        def __init__(self):
            documented_previous = np.array([
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 0.0],
            ])
            documented_current = documented_previous[:, [2, 1, 0]]
            # Emulate the state-major contiguous tag-array bridge view.
            previous = documented_previous.T.ravel().reshape(4, 3)
            current = documented_current.T.ravel().reshape(4, 3)
            self.data = {
                "OQP::td_bvec_mo_old": previous,
                "OQP::td_bvec_mo": current,
                "OQP::state_tracking_lineage": np.array([10, 11, 12]),
                "OQP::state_tracking_phase_initial": np.ones(3),
            }

    monkeypatch.setattr(single_point, "dump_log", lambda *_args, **_kwargs: None)
    tracker = single_point.NACME.__new__(single_point.NACME)
    tracker.mol = DummyMol()

    tracker.align_x(reorder=True)

    data = tracker.mol.data
    assert data["OQP::state_tracking_raw_order"].tolist() == [2, 1, 0]
    assert data["OQP::state_tracking_order"].tolist() == [0, 1, 2]
    assert np.allclose(data["OQP::td_bvec_mo"], data["OQP::td_bvec_mo_old"])


def test_initial_gauge_sign_is_not_double_counted_across_json_restart(monkeypatch):
    class DummyMol:
        def __init__(self):
            self.data = {
                # The previous JSON stores already transported response vectors.
                "OQP::td_bvec_mo_old": np.eye(3),
                "OQP::td_bvec_mo": np.diag([1.0, -1.0, 1.0]),
                "OQP::state_tracking_lineage_old": np.arange(3),
                "OQP::state_tracking_phase_initial_old": np.array([-1.0, 1.0, -1.0]),
            }

    monkeypatch.setattr(single_point, "dump_log", lambda *_args, **_kwargs: None)
    tracker = single_point.NACME.__new__(single_point.NACME)
    tracker.mol = DummyMol()

    tracker.align_x()

    data = tracker.mol.data
    # The current raw-to-initial correction is measured directly against the
    # corrected previous vectors.  Multiplying by the recorded previous sign
    # would double count its gauge operation.
    assert data["OQP::state_tracking_phase_initial"].tolist() == [1.0, -1.0, 1.0]
    assert data["OQP::state_tracking_previous_phase_initial"].tolist() == [-1.0, 1.0, -1.0]
    assert np.allclose(data["OQP::td_bvec_mo"], np.eye(3))


def test_public_tracking_contract_is_available_to_memory_drivers():
    class DummyMol:
        _state_tracking_fresh = True
        data = {
            "OQP::state_tracking_order": np.array([1, 0]),
            "OQP::state_tracking_raw_order": np.array([1, 0]),
            "OQP::state_tracking_output_reordered": np.array([0]),
            "OQP::state_tracking_lineage": np.array([7, 4]),
            "OQP::state_tracking_phase_step": np.array([-1.0, 1.0]),
            "OQP::state_tracking_phase_initial": np.array([-1.0, -1.0]),
            "OQP::state_tracking_previous_phase_initial": np.array([1.0, -1.0]),
            "OQP::state_tracking_overlap": np.array([0.98, 0.97]),
            "OQP::state_tracking_margin": np.array([0.90, 0.88]),
        }

    tracking = Molecule.get_state_tracking(DummyMol())

    assert tracking == {
        "schema_version": 1,
        "index_base": 0,
        "order_semantics": "current_root_to_previous_root",
        "order": [1, 0],
        "raw_order": [1, 0],
        "output_reordered": False,
        "lineage": [7, 4],
        "phase_step": [-1.0, 1.0],
        "phase_initial": [-1.0, -1.0],
        "previous_phase_initial": [1.0, -1.0],
        "matched_overlap": [0.98, 0.97],
        "margin": [0.90, 0.88],
    }


def test_tracking_loaded_from_guess_is_not_republished_as_current_result():
    tracking_tags = {
        "OQP::state_tracking_order": [1, 0],
        "OQP::state_tracking_raw_order": [1, 0],
        "OQP::state_tracking_output_reordered": [0],
        "OQP::state_tracking_lineage": [7, 4],
        "OQP::state_tracking_phase_step": [-1.0, 1.0],
        "OQP::state_tracking_phase_initial": [-1.0, -1.0],
        "OQP::state_tracking_previous_phase_initial": [1.0, -1.0],
        "OQP::state_tracking_overlap": [0.98, 0.97],
        "OQP::state_tracking_margin": [0.90, 0.88],
    }

    class DummyMol:
        tag = list(tracking_tags)
        config_tag = {}
        config = {}
        data = {}
        _state_tracking_fresh = True

    molecule = DummyMol()
    Molecule.put_data(molecule, tracking_tags)

    assert molecule._state_tracking_fresh is False
    assert Molecule.get_state_tracking(molecule) is None


def test_runner_results_forwards_public_tracking_contract():
    tracking = {"schema_version": 1, "phase_initial": [-1.0, 1.0]}

    class DummyMol:
        energies = [0.0]
        grads = []
        dcm = []
        nac = []
        soc = []

        @staticmethod
        def get_atoms():
            return ["H"]

        @staticmethod
        def get_system():
            return [[0.0, 0.0, 0.0]]

        @staticmethod
        def get_data():
            return {"OQP::VEC_MO_A": [1.0]}

        @staticmethod
        def get_state_tracking():
            return tracking

    runner = Runner.__new__(Runner)
    runner.mol = DummyMol()

    assert runner.results()["state_tracking"] is tracking
