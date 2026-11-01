"""NAMD trial-step history and energy audits."""
import copy

import numpy as np
import types

import pytest

from oqp.library import namd as mod


def test_energy_retry_restores_audit_history():
    driver = mod.NAMD.__new__(mod.NAMD)
    accepted = np.array([[0.0, 1.0], [-1.0, 0.0]])
    driver._nacme_reference_tdc = accepted.copy()
    driver._nacme_reference_source = 1
    saved = driver._energy_retry_state()

    driver._nacme_reference_tdc = 99.0*np.ones((2, 2))
    driver._nacme_reference_source = 127
    driver._restore_energy_retry_state(saved)

    np.testing.assert_array_equal(driver._nacme_reference_tdc, accepted)
    assert driver._nacme_reference_source == 1


def _nve_driver():
    driver = mod.NAMD.__new__(mod.NAMD)
    driver.ensemble = 'nve'
    driver.nve_gate = 'off'
    driver.nve_gate_abs_tol = 0.004
    driver.nve_gate_step_tol = 0.001
    driver.nve_gate_transition_tol = 1.0e-6
    driver.nve_gate_consecutive = 3
    driver._nve_reference_energy = None
    driver._nve_previous_energy = None
    driver._nve_gate_failures = 0
    driver.dt_adaptive = False
    driver._time_origin_fs = 0.0
    driver.dt_fs = 0.5
    driver._disc_energy_absorbed = 0.0
    driver._step_numerical_correction = 0.0
    return driver


def test_nve_gate_audits_the_energy_a_numerical_correction_absorbed():
    driver = _nve_driver()
    driver._update_nve_gate(0, -1.0, 0.0)
    # The dynamics gained 0.01 Ha; the numerical correction then rescaled the
    # velocities so that the corrected total equals the previous total.
    driver._step_numerical_correction = 0.01
    driver._disc_energy_absorbed = 0.01
    result = driver._update_nve_gate(1, -1.0, 0.0)
    assert result['step_change'] == pytest.approx(0.01)
    assert result['drift'] == pytest.approx(0.01)
    assert result['step_failure'] and result['drift_failure']

    # Without a correction the same corrected totals pass.
    clean = _nve_driver()
    clean._update_nve_gate(0, -1.0, 0.0)
    result = clean._update_nve_gate(1, -1.0, 0.0)
    assert not result['step_failure'] and not result['drift_failure']


def test_restart_histories_keep_nvt_energy_baseline():
    driver = mod.NAMD.__new__(mod.NAMD)
    driver.nstate = 2
    histories = {
        'ba_energy_left': None, 'ba_energy_center': None,
        'ba_tdc_left': None, 'ba_dt_left': None,
        # NVT keeps no NVE history, but energy recovery still needs a baseline.
        'nve_reference_energy': None, 'nve_previous_energy': None,
        'etot_prev': -76.25, 'disc_energy_absorbed': 0.003,
    }
    validated = driver._validate_restart_histories(
        copy.deepcopy(histories), context='test')
    assert validated['etot_prev'] == pytest.approx(-76.25)
    assert validated['disc_energy_absorbed'] == pytest.approx(0.003)

    payload = {}
    for name, value in validated.items():
        mod.NAMD._checkpoint_optional(payload, name, value)
    assert int(payload['has_etot_prev'][0]) == 1
    with pytest.raises(RuntimeError, match='invalid etot_prev'):
        driver._validate_restart_histories(
            dict(histories, etot_prev=float('nan')), context='test')


def _trajectory_driver(qm_atoms=None):
    driver = mod.NAMD.__new__(mod.NAMD)
    if qm_atoms is not None:
        driver.qm_atoms = np.asarray(qm_atoms, dtype=int)
    driver.mol = None
    return driver


def test_gas_phase_hop_direction_is_recorded_unchanged():
    driver = _trajectory_driver()
    driver._last_hop_direction = np.arange(9.0).reshape(3, 3)
    np.testing.assert_allclose(
        driver._trajectory_hop_direction(np.zeros((3, 3))),
        driver._last_hop_direction)


@pytest.mark.parametrize('qm_atoms', ([0, 1, 2, 3], [8, 9, 16, 17]))
def test_qmmm_hop_direction_maps_onto_the_full_atom_record(qm_atoms):
    """A QM/MM trajectory records every atom; d_IJ covers the QM atoms only.

    Writing the QM-sized vector into the full-length field raised
    "could not broadcast input array from shape (4,3) into shape (1,19,3)"
    and aborted every QM/MM NAMD run at its first trajectory frame.
    """
    driver = _trajectory_driver(qm_atoms)
    driver._last_hop_direction = np.arange(12.0).reshape(4, 3)
    coords = np.zeros((19, 3))
    mapped = driver._trajectory_hop_direction(coords)
    assert mapped.shape == coords.shape
    np.testing.assert_allclose(mapped[np.asarray(qm_atoms)],
                               driver._last_hop_direction)
    others = [i for i in range(19) if i not in qm_atoms]
    assert not mapped[others].any()


def test_unmappable_hop_direction_is_zeroed_not_fatal(monkeypatch):
    logged = []
    monkeypatch.setattr(mod, 'dump_log',
                        lambda *a, **k: logged.append(k.get('title', '')))
    driver = _trajectory_driver([0, 1])
    driver._last_hop_direction = np.arange(12.0).reshape(4, 3)
    mapped = driver._trajectory_hop_direction(np.zeros((19, 3)))
    assert mapped.shape == (19, 3) and not mapped.any()
    assert any('hop-direction shape' in title for title in logged)
