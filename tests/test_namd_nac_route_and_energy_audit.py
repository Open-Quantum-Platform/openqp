"""NAMD analytic-NAC route selection, trial-step history and energy audits."""
import copy

import numpy as np
import types

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
    assert mod.analytic_nac_route_issue(
        _mrsf_singlet_config(input={'qmmm_flag': True})) is None
    unsupported = {
        'legacy default convergence': _mrsf_singlet_config(scf={'conv': 1.0e-6}),
        'loose response convergence': _mrsf_singlet_config(tdhf={'conv': 1.0e-6}),
        'triplet MRSF states': _mrsf_singlet_config(tdhf={'multiplicity': 3}),
        'closed-shell reference': _mrsf_singlet_config(scf={'type': 'rhf', 'multiplicity': 1}),
        'linear-response TDDFT': _mrsf_singlet_config(tdhf={'type': 'rpa'}),
        'spin-orbit NAMD': _mrsf_singlet_config(md={'soc': True}),
        'tight binding': _mrsf_singlet_config(input={'method': 'dftb'}),
    }
    for label, config in unsupported.items():
        assert mod.analytic_nac_route_issue(config) is not None, label
    # SOC keeps its own NotImplementedError guard; the model check reports only
    # the electronic-structure limitation.  QM/MM is a route, not a different
    # MRSF electronic model.
    assert mod.analytic_nac_model_issue(
        _mrsf_singlet_config(md={'soc': True})) is None
    assert mod.analytic_nac_model_issue(
        _mrsf_singlet_config(input={'qmmm_flag': True})) is None


def test_qmmm_analytic_nac_rejects_link_atoms_until_direction_is_projected():
    driver = mod.NAMD_QMMM.__new__(mod.NAMD_QMMM)
    driver.mol = types.SimpleNamespace(
        config={'md': {'rescale': 'hop_analytic_nac'}})      # named by the input
    driver.tdc_provider = 'npi'
    driver.rescale_provider = 'hop_analytic_nac'
    driver.nacme_check = 'off'
    driver.link_atoms = [object()]
    with pytest.raises(NotImplementedError, match='without link atoms'):
        driver._validate_qmmm_analytic_nac_scope()

    driver.link_atoms = []
    driver._validate_qmmm_analytic_nac_scope()


def _link_atom_driver(requested, provider, tdc='npi', nacme_check='off'):
    driver = mod.NAMD_QMMM.__new__(mod.NAMD_QMMM)
    driver.mol = types.SimpleNamespace(config={'md': {'rescale': requested}})
    driver.tdc_provider = tdc
    driver.rescale_provider = provider
    driver.nacme_check = nacme_check
    driver._rescale_auto_issue = None
    driver.link_atoms = [object()]
    return driver


def test_rescale_auto_falls_back_to_isotropic_for_a_link_atom_region():
    """rescale=auto is a default, not a request for the analytic route.

    It is resolved in NAMD.__init__ by a gate that only looks at the
    electronic model, before the QM/MM subclass knows about link atoms.  With
    scf/tdhf conv <= 1e-8 it therefore became hop_analytic_nac, which the
    link-atom guard then rejected -- so tightening the convergence turned a
    working link-atom trajectory into a startup error.
    """
    driver = _link_atom_driver('auto', 'hop_analytic_nac')
    driver._validate_qmmm_analytic_nac_scope()
    assert driver.rescale_provider == 'isotropic'
    assert driver._rescale_auto_issue == 'QM region has link atoms'

    # the default spelling (key absent) behaves the same
    driver = _link_atom_driver('auto', 'hop_analytic_nac')
    driver.mol.config = {}
    driver._validate_qmmm_analytic_nac_scope()
    assert driver.rescale_provider == 'isotropic'

    # without link atoms auto keeps the analytic route it resolved to
    driver = _link_atom_driver('auto', 'hop_analytic_nac')
    driver.link_atoms = []
    driver._validate_qmmm_analytic_nac_scope()
    assert driver.rescale_provider == 'hop_analytic_nac'


@pytest.mark.parametrize("requested, provider, tdc, nacme_check", (
    ('hop_analytic_nac', 'hop_analytic_nac', 'npi', 'off'),   # named rescale
    ('analytic_nac', 'analytic_nac', 'npi', 'off'),
    ('auto', 'hop_analytic_nac', 'analytic', 'off'),          # named analytic TDC
    ('auto', 'hop_analytic_nac', 'npi', 'analytic'),          # named analytic audit
))
def test_named_analytic_modes_are_still_rejected_with_link_atoms(
        requested, provider, tdc, nacme_check):
    driver = _link_atom_driver(requested, provider, tdc, nacme_check)
    with pytest.raises(NotImplementedError, match='without link atoms'):
        driver._validate_qmmm_analytic_nac_scope()


def test_qmmm_analytic_nac_contracts_the_current_qm_velocity():
    driver = mod.NAMD_QMMM.__new__(mod.NAMD_QMMM)
    current = np.array([[1.0, 2.0, 3.0]])
    seen = {}
    driver.vel = np.zeros_like(current)
    driver._qm_velocities = lambda: current.copy()

    def state_overlap(step):
        seen['step'] = step
        seen['velocity'] = driver.vel.copy()
        return 'overlap'

    driver._state_overlap = state_overlap
    assert driver._prepare_qmmm_hop_couplings(7) == 'overlap'
    assert seen['step'] == 7
    np.testing.assert_array_equal(seen['velocity'], current)


def test_minimal_mrsf_namd_input_resolves_auto_to_hop_analytic_nac():
    from oqp.utils import oqp_input

    spec = oqp_input.parse_canonical_oqp(
        'mrsf/bhhlyp/6-31g* geom="thymine.xyz" '
        'namd(S2,scheme=TDC_NAC,nstep=1000,dt=0.5,'
        'velocity="velocity.au")')
    config = oqp_input.lower_to_legacy(spec)
    values = oqp_input._effective_config(config, oqp_input._load_schema_defaults())
    effective = {}
    for (section, key), value in values.items():
        effective.setdefault(section, {})[key] = value
    assert effective['md']['rescale'] == 'hop_analytic_nac'
    assert mod.analytic_nac_route_issue(effective) is None
    effective['scf']['conv'] = 1e-6
    assert mod.analytic_nac_route_issue(effective) is not None


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
