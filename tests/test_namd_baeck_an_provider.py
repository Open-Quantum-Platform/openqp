"""Causal signed Baeck-An TDC provider tests."""

import numpy as np


def test_baeck_an_sign_comes_only_from_phase_tracked_overlap():
    from oqp.library.namd import NAMD

    magnitude = np.array([
        [0.0, -0.4, 0.2],
        [0.4, 0.0, -0.3],
        [-0.2, 0.3, 0.0],
    ])
    overlap = np.array([
        [0.0, 0.7, 0.0],
        [-0.7, 0.0, -0.1],
        [0.0, 0.1, 0.0],
    ])
    signed = NAMD._signed_baeck_an_tdc(magnitude, overlap)

    np.testing.assert_allclose(signed, np.array([
        [0.0, 0.4, 0.0],
        [-0.4, 0.0, -0.3],
        [0.0, 0.3, 0.0],
    ]))
    np.testing.assert_array_equal(signed + signed.T, np.zeros((3, 3)))


def test_baeck_an_provider_uses_npi_only_for_phase_transport():
    from oqp.library.namd import NAMD

    angle = 0.2
    overlap = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)],
    ])
    driver = NAMD.__new__(NAMD)
    driver.tdc_scheme = 3
    driver.dt = 2.0

    tdc = driver._compute_tdc(overlap)
    np.testing.assert_allclose(
        tdc, np.array([[0.0, -angle/2.0], [angle/2.0, 0.0]]),
        atol=1.0e-13, rtol=0.0)


def test_three_energy_points_activate_one_step_lagged_baeck_an():
    from types import SimpleNamespace
    from oqp.library.namd import NAMD

    driver = NAMD.__new__(NAMD)
    driver.nacme_check = 'off'
    driver.tdc_provider = 'baeck_an'
    driver.nstate = 2
    driver.dt = 1.0
    driver.ba_gap_max = 1.0
    driver._ba_energy_left = None
    driver._ba_energy_center = None
    driver._ba_tdc_left = None
    driver._ba_dt_left = None
    driver._ba_last = None
    driver._nacme_gate_failures = 0
    driver._nacme_gate_last = None
    driver._nacme_candidate_tdc = None
    driver._nacme_reference_tdc = None
    driver._nacme_reference_mask = None
    driver._nacme_reference_source = 0
    driver._pending_nacme_gate_error = None
    overlap_tdc = np.array([[0.0, 0.1], [-0.1, 0.0]])
    driver._compute_tdc = lambda _overlap: overlap_tdc.copy()
    driver.mol = SimpleNamespace(data={
        'OQP::td_energies_old': np.array([0.0, 0.03]),
        'OQP::td_energies': np.array([0.0, 0.02]),
    })

    driver._update_baeck_an_check(1, np.eye(2))
    assert driver._last_baeck_an_tdc is None
    assert driver._last_tdc_source == 1

    driver.mol.data['OQP::td_energies_old'] = np.array([0.0, 0.02])
    driver.mol.data['OQP::td_energies'] = np.array([0.0, 0.03])
    driver._update_baeck_an_check(2, np.eye(2))
    np.testing.assert_allclose(
        driver._last_baeck_an_tdc,
        np.array([[0.0, 0.5], [-0.5, 0.0]]))
    assert driver._last_tdc_source == 3
    assert driver._ba_last is None


def test_hop_uses_lagged_baeck_an_matrix_and_isotropic_rescaling(monkeypatch):
    import oqp
    from oqp.library.namd import NAMD, _P_RESCALE, _P_TDC

    class Mol:
        def __init__(self):
            self.data = {
                'OQP::td_states_overlap': np.eye(2),
                'OQP::td_energies': np.array([-0.5, -0.4]),
            }

    driver = NAMD.__new__(NAMD)
    driver.mol = Mol()
    driver.nstate = 2
    driver.natom = 1
    driver.coef = np.array([1.0 + 0.0j, 0.0 + 0.0j])
    driver.vel = np.array([[0.1, 0.2, 0.3]])
    driver.dt_fs = 0.5
    driver.substep = 4
    driver.thrshe = 1.0
    driver.active = 1
    driver.decoherence = 0
    driver.edc_c = 0.1
    driver.tdc_scheme = 3
    driver.tdc_provider = 'baeck_an'
    driver.rescale_provider = 'isotropic'
    driver.trivial = 0
    driver.trivial_thresh = 0.5
    driver._last_baeck_an_tdc = np.array([[0.0, 0.7], [-0.7, 0.0]])
    driver._last_overlap_tdc = np.array([[0.0, 0.1], [-0.1, 0.0]])
    driver._hop_random = lambda: 0.25

    captured = {}

    def fake_hop(mol):
        captured['tdc'] = np.array(mol.data['OQP::namd_tdc'], copy=True)
        captured['params'] = np.array(mol.data['OQP::namd_params'], copy=True)
        results = np.zeros(2*2 + 8)
        results[2*2 + 5] = 0.0
        mol.data['OQP::namd_results'] = results

    monkeypatch.setattr(oqp, 'mrsf_namd_hop', fake_hop)
    new_active, hopped = driver._hop()

    np.testing.assert_allclose(
        captured['tdc'], driver._last_baeck_an_tdc.reshape(-1))
    assert captured['params'][_P_TDC] == 3.0
    assert captured['params'][_P_RESCALE] == 0.0
    assert driver._last_tdc_source == 3
    assert new_active == 1
    assert hopped is False


def test_baeck_an_warmup_uses_overlap_npi_and_records_source(monkeypatch):
    import oqp
    from oqp.library.namd import NAMD

    class Mol:
        def __init__(self):
            self.data = {
                'OQP::td_states_overlap': np.eye(2),
                'OQP::td_energies': np.array([-0.5, -0.4]),
            }

    driver = NAMD.__new__(NAMD)
    driver.mol = Mol()
    driver.nstate = 2
    driver.natom = 1
    driver.coef = np.array([1.0 + 0.0j, 0.0 + 0.0j])
    driver.vel = np.array([[0.1, 0.2, 0.3]])
    driver.dt_fs = 0.5
    driver.substep = 4
    driver.thrshe = 1.0
    driver.active = 1
    driver.decoherence = 0
    driver.edc_c = 0.1
    driver.tdc_scheme = 3
    driver.tdc_provider = 'baeck_an'
    driver.rescale_provider = 'isotropic'
    driver.trivial = 0
    driver.trivial_thresh = 0.5
    driver._last_baeck_an_tdc = None
    driver._last_overlap_tdc = np.array([[0.0, 0.1], [-0.1, 0.0]])
    driver._hop_random = lambda: 0.25

    def fake_hop(mol):
        results = np.zeros(2*2 + 8)
        results[2*2 + 5] = 0.0
        mol.data['OQP::namd_results'] = results

    monkeypatch.setattr(oqp, 'mrsf_namd_hop', fake_hop)
    driver._hop()
    assert driver._last_tdc_source == 1
