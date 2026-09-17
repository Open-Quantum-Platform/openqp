"""NAMD analytic-NAC route selection, trial-step history and energy audits."""
import copy

import numpy as np
import pytest

from oqp.library import namd as mod


def _mrsf_singlet_config(**sections):
    config = {
        'input': {'method': 'tdhf', 'qmmm_flag': False},
        'scf': {'type': 'rohf', 'multiplicity': 3, 'conv': 1.0e-8},
        'tdhf': {'type': 'mrsf', 'multiplicity': 1, 'conv': 1.0e-8},
        'md': {'soc': False},
    }
    for section, values in sections.items():
        config[section] = dict(config.get(section, {}), **values)
    return config


def test_analytic_nac_route_is_available_only_where_it_is_defined():
    assert mod.analytic_nac_route_issue(_mrsf_singlet_config()) is None
    unsupported = {
        'legacy default convergence': _mrsf_singlet_config(scf={'conv': 1.0e-6}),
        'loose response convergence': _mrsf_singlet_config(tdhf={'conv': 1.0e-6}),
        'triplet MRSF states': _mrsf_singlet_config(tdhf={'multiplicity': 3}),
        'closed-shell reference': _mrsf_singlet_config(scf={'type': 'rhf', 'multiplicity': 1}),
        'linear-response TDDFT': _mrsf_singlet_config(tdhf={'type': 'rpa'}),
        'spin-orbit NAMD': _mrsf_singlet_config(md={'soc': True}),
        'QM/MM NAMD': _mrsf_singlet_config(input={'qmmm_flag': True}),
        'tight binding': _mrsf_singlet_config(input={'method': 'dftb'}),
    }
    for label, config in unsupported.items():
        assert mod.analytic_nac_route_issue(config) is not None, label
    # SOC and QM/MM keep their own NotImplementedError guards; the model
    # check reports only the electronic-structure limitation.
    assert mod.analytic_nac_model_issue(
        _mrsf_singlet_config(md={'soc': True})) is None
    assert mod.analytic_nac_model_issue(
        _mrsf_singlet_config(input={'qmmm_flag': True})) is None


def test_energy_retry_restores_analytic_audit_history():
    driver = mod.NAMD.__new__(mod.NAMD)
    accepted = np.array([[0.0, 1.0], [-1.0, 0.0]])
    driver._analytic_tdc_previous = accepted.copy()
    driver._analytic_tdc_centered = accepted.copy()
    driver._nacme_reference_tdc = accepted.copy()
    driver._nacme_reference_source = 2
    saved = driver._energy_retry_state()

    rejected = 99.0*np.ones((2, 2))
    driver._analytic_tdc_previous = rejected
    driver._analytic_tdc_centered = rejected
    driver._nacme_reference_tdc = rejected
    driver._nacme_reference_source = 127
    driver._restore_energy_retry_state(saved)

    np.testing.assert_array_equal(driver._analytic_tdc_previous, accepted)
    np.testing.assert_array_equal(driver._analytic_tdc_centered, accepted)
    np.testing.assert_array_equal(driver._nacme_reference_tdc, accepted)
    assert driver._nacme_reference_source == 2


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


def test_restart_histories_keep_nvt_energy_baseline_and_analytic_history():
    driver = mod.NAMD.__new__(mod.NAMD)
    driver.nstate = 2
    histories = {
        'ba_energy_left': None, 'ba_energy_center': None,
        'ba_tdc_left': None, 'ba_dt_left': None,
        # NVT keeps no NVE history, but energy recovery still needs a baseline.
        'nve_reference_energy': None, 'nve_previous_energy': None,
        'analytic_tdc_previous': np.array([[0.0, 0.2], [-0.2, 0.0]]),
        'etot_prev': -76.25, 'disc_energy_absorbed': 0.003,
    }
    validated = driver._validate_restart_histories(
        copy.deepcopy(histories), context='test')
    assert validated['etot_prev'] == pytest.approx(-76.25)
    assert validated['disc_energy_absorbed'] == pytest.approx(0.003)
    np.testing.assert_array_equal(
        validated['analytic_tdc_previous'], histories['analytic_tdc_previous'])

    payload = {}
    for name, value in validated.items():
        mod.NAMD._checkpoint_optional(payload, name, value)
    assert int(payload['has_etot_prev'][0]) == 1
    assert int(payload['has_analytic_tdc_previous'][0]) == 1
    with pytest.raises(RuntimeError, match='invalid etot_prev'):
        driver._validate_restart_histories(
            dict(histories, etot_prev=float('nan')), context='test')
