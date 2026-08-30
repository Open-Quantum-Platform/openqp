"""Analytic-NAC TDC and directional FSSH rescaling regression gates."""

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_directional_rescaling_contract_is_native_and_explicit():
    source = (ROOT / "source" / "modules" / "namd.F90").read_text()
    tags = (ROOT / "source" / "tagarray_driver.F90").read_text()
    header = (ROOT / "include" / "oqp.h").read_text()
    driver = (ROOT / "pyoqp" / "oqp" / "library" / "namd.py").read_text()

    assert "subroutine namd_rescale_velocities_directional(" in source
    assert 'bind(C, name="oqp_namd_rescale_directional")' in source
    assert "int oqp_namd_rescale_directional(" in header
    assert "OQP_namd_dcv" in tags
    assert "direction_vectors(:,:,i,ncrst)" in source
    assert "if (tdc_scheme == 0 .and. all(abs(tdc) < 1.0e-30_dp))" in source
    assert "np.einsum('ijac,ac->ij'" in driver
    assert "self.tdc_provider == 'analytic'" in driver
    assert "self.rescale_provider == 'analytic_nac'" in driver
    assert "if (mode == 2) then" in source
    assert "hop_analytic_nac" in driver


def test_built_directional_kernel_conserves_energy_and_is_gauge_invariant():
    script = r"""
import json
import numpy as np
import oqp

def call(velocity, mass, direction, delta_e):
    velocity = np.array(velocity, dtype=np.float64, order='C', copy=True)
    mass = np.ascontiguousarray(mass, dtype=np.float64)
    direction = np.ascontiguousarray(direction, dtype=np.float64)
    initial = velocity.copy()
    gamma = np.zeros(1, dtype=np.float64)
    disc = np.zeros(1, dtype=np.float64)
    status = oqp.oqp_namd_rescale_directional(
        mass.size,
        oqp.ffi.cast('double *', velocity.ctypes.data),
        oqp.ffi.cast('double *', mass.ctypes.data),
        oqp.ffi.cast('double *', direction.ctypes.data),
        delta_e,
        oqp.ffi.cast('double *', gamma.ctypes.data),
        oqp.ffi.cast('double *', disc.ctypes.data),
    )
    e0 = 0.5*np.sum(mass[:, None]*initial**2)
    e1 = 0.5*np.sum(mass[:, None]*velocity**2)
    return {
        'status': int(status), 'velocity': velocity.tolist(),
        'gamma': float(gamma[0]), 'disc': float(disc[0]),
        'energy_residual': float(e1 + delta_e - e0),
        'changed': float(np.max(np.abs(velocity-initial))),
    }

mass = np.array([1836.0, 29156.0])
velocity = np.array([[0.003, -0.001, 0.0005], [-0.001, 0.002, 0.001]])
direction = np.array([[0.8, -0.3, 0.2], [-0.4, 0.6, -0.5]])
base = call(velocity, mass, direction, 5.0e-4)
sign = call(velocity, mass, -direction, 5.0e-4)
scale = call(velocity, mass, 7.5*direction, 5.0e-4)
frustrated = call(np.zeros_like(velocity), mass, direction, 5.0e-4)
downhill = call(np.zeros_like(velocity), mass, direction, -5.0e-4)
zero = call(velocity, mass, np.zeros_like(direction), 0.0)
nonfinite = call(velocity, mass, np.full_like(direction, np.nan), 0.0)
print('DIRECTIONAL=' + json.dumps({
    'base': base, 'sign': sign, 'scale': scale,
    'frustrated': frustrated, 'downhill': downhill,
    'zero': zero, 'nonfinite': nonfinite,
}))
"""
    env = os.environ.copy()
    pythonpath = str(ROOT / "pyoqp")
    if env.get("PYTHONPATH"):
        pythonpath += os.pathsep + env["PYTHONPATH"]
    env["PYTHONPATH"] = pythonpath
    result = subprocess.run(
        [sys.executable, "-c", script], cwd=ROOT, env=env,
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        if ("No module named" in result.stderr or "cannot load" in result.stderr
                or "oqp_namd_rescale_directional" in result.stderr):
            pytest.skip("matching compiled OpenQP runtime is not available")
        pytest.fail(result.stdout + result.stderr)
    marker = next(
        (line for line in result.stdout.splitlines()
         if line.startswith("DIRECTIONAL=")), None)
    assert marker is not None, result.stdout + result.stderr
    values = json.loads(marker.removeprefix("DIRECTIONAL="))

    base = values['base']
    assert base['status'] == 0
    assert abs(base['energy_residual']) < 2.0e-14
    np.testing.assert_allclose(
        values['sign']['velocity'], base['velocity'], rtol=0.0, atol=2.0e-15)
    np.testing.assert_allclose(
        values['scale']['velocity'], base['velocity'], rtol=0.0, atol=2.0e-15)
    assert values['frustrated']['status'] == 1
    assert values['frustrated']['changed'] == 0.0
    assert values['downhill']['status'] == 0
    assert abs(values['downhill']['energy_residual']) < 2.0e-14
    for key in ('zero', 'nonfinite'):
        assert values[key]['status'] == 1
        assert values[key]['changed'] == 0.0


def test_analytic_tdc_contraction_is_signed_and_antisymmetric():
    rng = np.random.default_rng(81)
    raw = rng.normal(size=(3, 3, 4, 3))
    dcv = 0.5*(raw - raw.swapaxes(0, 1))
    velocity = rng.normal(size=(4, 3))
    tdc = np.einsum('ijac,ac->ij', dcv, velocity)
    np.testing.assert_allclose(np.diag(tdc), 0.0, atol=1.0e-15)
    np.testing.assert_allclose(tdc + tdc.T, 0.0, atol=1.0e-15)


def test_analytic_provider_keeps_endpoint_and_centered_tdc_separate(monkeypatch):
    import oqp.library.nac_analytic as nac_module
    from oqp.library.namd import NAMD

    dcv = np.zeros((2, 2, 1, 3))
    dcv[0, 1, 0] = [2.0, -1.0, 0.5]
    dcv[1, 0] = -dcv[0, 1]
    monkeypatch.setattr(
        nac_module, 'analytic_nac',
        lambda _mol: (np.zeros_like(dcv), dcv.copy()))

    driver = NAMD.__new__(NAMD)
    driver.mol = object()
    driver.nstate = 2
    driver.natom = 1
    driver.vel = np.array([[0.2, 0.3, -0.4]])
    driver.nacme_check = 'analytic'
    driver._last_overlap_tdc = np.array([[0.0, 0.25], [-0.25, 0.0]])
    driver._analytic_tdc_previous = np.array([[0.0, -0.10], [0.10, 0.0]])
    driver._nacme_reference_tdc = None
    driver._nacme_reference_mask = None
    driver._nacme_reference_source = 0
    captured = {}
    driver._run_nacme_gate = lambda candidate, reference, **kwargs: (
        captured.update(candidate=candidate.copy(), reference=reference.copy(),
                        kwargs=kwargs) or {'verdict': 'pass'})
    monkeypatch.setattr(
        'oqp.library.namd.dump_log', lambda *_args, **_kwargs: None)

    driver._update_analytic_nac(4, compare_overlap=True)
    endpoint = np.array([[0.0, -0.1], [0.1, 0.0]])
    centered = np.array([[0.0, -0.1], [0.1, 0.0]])
    np.testing.assert_allclose(driver._last_analytic_tdc, endpoint)
    np.testing.assert_allclose(driver._analytic_tdc_centered, centered)
    np.testing.assert_allclose(captured['candidate'], driver._last_overlap_tdc)
    np.testing.assert_allclose(captured['reference'], centered)
    assert captured['kwargs']['signed'] is True
    assert captured['kwargs']['source'] == 'analytic'


def test_hop_passes_analytic_tdc_and_full_direction_tensor(monkeypatch):
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
    driver.tdc_scheme = 2
    driver.tdc_provider = 'analytic'
    driver.rescale_provider = 'analytic_nac'
    driver.trivial = 0
    driver.trivial_thresh = 0.5
    driver._last_analytic_tdc = np.array([[0.0, 0.7], [-0.7, 0.0]])
    driver._last_analytic_dcv = np.arange(12.0).reshape(2, 2, 1, 3)
    driver._hop_random = lambda: 0.25

    captured = {}
    def fake_hop(mol):
        captured['tdc'] = np.array(mol.data['OQP::namd_tdc'], copy=True)
        captured['dcv'] = np.array(mol.data['OQP::namd_dcv'], copy=True)
        captured['params'] = np.array(mol.data['OQP::namd_params'], copy=True)
        results = np.zeros(2*2 + 8)
        results[2*2 + 5] = 1.0
        results[2*2 + 6] = 0.2
        results[2*2 + 7] = 0.3
        mol.data['OQP::namd_results'] = results
    monkeypatch.setattr(oqp, 'mrsf_namd_hop', fake_hop)

    new_active, hopped = driver._hop()
    np.testing.assert_allclose(
        captured['tdc'], driver._last_analytic_tdc.reshape(-1))
    np.testing.assert_allclose(
        captured['dcv'], driver._last_analytic_dcv.reshape(-1))
    assert captured['params'][_P_TDC] == 2.0
    assert captured['params'][_P_RESCALE] == 1.0
    assert new_active == 1
    assert hopped is False


def _hop_triggered_driver():
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
    driver.mass = np.array([1836.0])
    driver.coef = np.array([1.0 + 0.0j, 0.0 + 0.0j])
    driver.vel = np.array([[0.1, 0.2, 0.3]])
    driver.dt_fs = 0.5
    driver.dt = 0.5 * 41.3413745758
    driver.substep = 4
    driver.thrshe = 1.0
    driver.active = 1
    driver.decoherence = 0
    driver.edc_c = 0.1
    driver.tdc_scheme = 0
    driver.tdc_provider = 'fd'
    driver.rescale_provider = 'hop_analytic_nac'
    driver.trivial = 0
    driver.trivial_thresh = 0.5
    driver._last_analytic_tdc = None
    driver._last_analytic_dcv = None
    driver._analytic_tdc_previous = None
    driver._analytic_tdc_centered = None
    driver._nacme_reference_tdc = None
    driver._nacme_reference_mask = None
    driver._nacme_reference_source = 0
    driver._last_rescale_source = 2
    driver._last_rescale_gamma = np.nan
    driver._last_rescale_discriminant = np.nan
    driver._last_hop_direction = np.zeros((1, 3))
    return driver


def _fake_native_hop_result(mol, *, target):
    params = np.array(mol.data['OQP::namd_params'], copy=True)
    n = int(round(params[12]))
    params[10] = 0.0
    params[11] = float(target)
    mol.data['OQP::namd_params'] = params
    results = np.zeros(n*n + 8)
    results[n*n + 1] = float(target)
    results[n*n + 5] = 2.0
    results[n*n + 6] = 0.0
    results[n*n + 7] = -1.0
    mol.data['OQP::namd_results'] = results


def test_hop_triggered_analytic_evaluates_exact_nac_once_for_candidate(monkeypatch):
    import oqp
    from oqp.library.namd import _P_RESCALE

    driver = _hop_triggered_driver()
    calls = {'random': 0, 'exact': 0, 'rescale': 0}

    def random_once():
        calls['random'] += 1
        return 0.25

    def fake_native(mol):
        assert np.all(np.asarray(mol.data['OQP::namd_dcv']) == 0.0)
        assert mol.data['OQP::namd_params'][_P_RESCALE] == 2.0
        _fake_native_hop_result(mol, target=2)

    def fake_exact(_istep, *, compare_overlap=False):
        assert compare_overlap is False
        calls['exact'] += 1
        dcv = np.zeros((2, 2, 1, 3))
        dcv[0, 1, 0] = [1.0, 0.0, 0.0]
        dcv[1, 0, 0] = -dcv[0, 1, 0]
        driver._last_analytic_dcv = dcv

    def fake_rescale(_natom, velocity_ptr, _mass_ptr, direction_ptr,
                     _delta_e, gamma_ptr, discriminant_ptr):
        calls['rescale'] += 1
        velocity = np.frombuffer(oqp.ffi.buffer(velocity_ptr, 3*8), dtype=np.float64)
        direction = np.frombuffer(oqp.ffi.buffer(direction_ptr, 3*8), dtype=np.float64)
        gamma = np.frombuffer(oqp.ffi.buffer(gamma_ptr, 8), dtype=np.float64)
        discriminant = np.frombuffer(
            oqp.ffi.buffer(discriminant_ptr, 8), dtype=np.float64)
        np.testing.assert_allclose(direction, [1.0, 0.0, 0.0])
        velocity[0] += 0.01
        gamma[0] = 0.125
        discriminant[0] = 0.5
        return 0

    driver._hop_random = random_once
    driver._update_analytic_nac = fake_exact
    monkeypatch.setattr(oqp, 'mrsf_namd_hop', fake_native)
    monkeypatch.setattr(oqp, 'oqp_namd_rescale_directional', fake_rescale)

    new_active, hopped = driver._hop(allow_hop=True, istep=7)
    assert (new_active, hopped) == (2, True)
    assert calls == {'random': 1, 'exact': 1, 'rescale': 1}
    assert driver.vel[0, 0] == pytest.approx(0.11)
    assert driver._last_rescale_source == 2
    assert driver._last_rescale_gamma == pytest.approx(0.125)
    assert driver._last_rescale_discriminant == pytest.approx(0.5)
    np.testing.assert_allclose(driver._last_hop_direction, [[1.0, 0.0, 0.0]])


def test_hop_triggered_analytic_skips_exact_nac_without_candidate(monkeypatch):
    import oqp

    driver = _hop_triggered_driver()
    driver._hop_random = lambda: 0.75
    driver._update_analytic_nac = lambda *_args, **_kwargs: pytest.fail(
        'exact analytic NAC was evaluated without an FSSH candidate')
    monkeypatch.setattr(
        oqp, 'mrsf_namd_hop',
        lambda mol: _fake_native_hop_result(mol, target=1))

    new_active, hopped = driver._hop(allow_hop=True, istep=8)
    assert (new_active, hopped) == (1, False)


def test_hop_triggered_record_clear_prevents_stale_exact_vector():
    driver = _hop_triggered_driver()
    driver._last_analytic_dcv = np.ones((2, 2, 1, 3))
    driver._last_analytic_tdc = np.ones((2, 2))
    driver._analytic_tdc_previous = np.ones((2, 2))
    driver._analytic_tdc_centered = np.ones((2, 2))
    driver._nacme_reference_tdc = np.ones((2, 2))
    driver._nacme_reference_mask = np.ones((2, 2), dtype=np.int32)
    driver._nacme_reference_source = 2
    driver._clear_hop_triggered_analytic_record()
    assert driver._last_analytic_dcv is None
    assert driver._last_analytic_tdc is None
    assert driver._analytic_tdc_previous is None
    assert driver._analytic_tdc_centered is None
    assert driver._nacme_reference_tdc is None
    assert driver._nacme_reference_mask is None
    assert driver._nacme_reference_source == 0


def test_built_hop_triggered_kernel_returns_uncommitted_candidate(tmp_path):
    """Exercise native rescale mode 2, not the Python test double above."""
    script = r"""
import json
import numpy as np
import oqp
from oqp.pyoqp import Runner

runner = Runner(
    project='native_ht_contract',
    input_file='examples/ODP/H2CO_BHHLYP-MRSF-NAMD-ODP.inp',
    log=LOG_PATH, silent=1, usempi=False)
mol = runner.mol
n = 2
natom = int(mol.data['natom'])
coef = np.array([2**-0.5, 0.0, 2**-0.5, 0.0])
velocity = np.linspace(1.0e-4, 1.2e-3, 3*natom)
params = np.zeros(16)
params[:15] = [0.5, 1, 1.0, 1.0e-8, 1, 0, 0.1, 2,
               0, 0.5, 0, 1, 2, 1, 2]
mol.data['OQP::namd_coef'] = coef.copy()
mol.data['OQP::namd_velocity'] = velocity.copy()
mol.data['OQP::namd_params'] = params
mol.data['OQP::namd_tdc'] = np.array([0.0, -0.05, 0.05, 0.0])
mol.data['OQP::namd_eabs'] = np.array([0.0, 0.0])
mol.data['OQP::namd_stas'] = np.eye(n).reshape(-1)
mol.data['OQP::namd_dcv'] = np.zeros(n*n*natom*3)
oqp.mrsf_namd_hop(mol)
params_out = np.asarray(mol.data['OQP::namd_params'])
velocity_out = np.asarray(mol.data['OQP::namd_velocity'])
results = np.asarray(mol.data['OQP::namd_results'])
print('HT_NATIVE=' + json.dumps({
    'active': int(round(params_out[4])),
    'hopped': int(round(params_out[10])),
    'target': int(round(params_out[11])),
    'result_target': int(round(results[n*n + 1])),
    'rescale_source': int(round(results[n*n + 5])),
    'velocity_unchanged': bool(np.array_equal(velocity_out, velocity)),
    'candidate_probability': float(results[2]),
}))
""".replace('LOG_PATH', repr(str(tmp_path / 'native_ht_contract.log')))
    env = os.environ.copy()
    pythonpath = str(ROOT / 'pyoqp')
    if env.get('PYTHONPATH'):
        pythonpath += os.pathsep + env['PYTHONPATH']
    env['PYTHONPATH'] = pythonpath
    result = subprocess.run(
        [sys.executable, '-c', script], cwd=ROOT, env=env,
        capture_output=True, text=True, check=False)
    if result.returncode != 0:
        if ('No module named' in result.stderr or 'cannot load' in result.stderr
                or 'mode 2' in result.stderr):
            pytest.skip('matching compiled HT-NAC runtime is not available')
        pytest.fail(result.stdout + result.stderr)
    marker = next(
        line for line in result.stdout.splitlines()
        if line.startswith('HT_NATIVE='))
    values = json.loads(marker.removeprefix('HT_NATIVE='))
    assert values == {
        'active': 1,
        'hopped': 0,
        'target': 2,
        'result_target': 2,
        'rescale_source': 2,
        'velocity_unchanged': True,
        'candidate_probability': pytest.approx(0.5097010658),
    }
