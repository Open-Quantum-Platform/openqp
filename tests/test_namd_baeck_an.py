"""Resident TD-Baeck-An NACME diagnostic tests."""

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_baeck_an_kernel_is_fortran_resident_and_c_interoperable():
    source = (ROOT / "source" / "modules" / "namd.F90").read_text()
    header = (ROOT / "include" / "oqp.h").read_text()
    driver = (ROOT / "pyoqp" / "oqp" / "library" / "namd.py").read_text()

    assert "function namd_baeck_an_tdc(" in source
    assert 'bind(C, name="oqp_namd_baeck_an_tdc")' in source
    assert "int oqp_namd_baeck_an_tdc(" in header
    assert "oqp.oqp_namd_baeck_an_tdc(" in driver
    assert "magnitude diagnostic" in driver
    assert "function namd_nacme_gate(" in source
    assert 'bind(C, name="oqp_namd_nacme_gate")' in source
    assert "int oqp_namd_nacme_gate(" in header
    assert "def _run_nacme_gate(" in driver
    assert "signed=True" in driver


def test_built_baeck_an_kernel_matches_quadratic_gap_oracle():
    script = r"""
import json
import numpy as np
import oqp

def run(old, center, current, dt_left=1.0, dt_right=1.0, gap_max=1.0):
    old = np.ascontiguousarray(old, dtype=np.float64)
    center = np.ascontiguousarray(center, dtype=np.float64)
    current = np.ascontiguousarray(current, dtype=np.float64)
    out = np.zeros((old.size, old.size), dtype=np.float64)
    status = oqp.oqp_namd_baeck_an_tdc(
        old.size, dt_left, dt_right, gap_max,
        oqp.ffi.cast('double *', old.ctypes.data),
        oqp.ffi.cast('double *', center.ctypes.data),
        oqp.ffi.cast('double *', current.ctypes.data),
        oqp.ffi.cast('double *', out.ctypes.data),
    )
    return status, out.tolist()

def gate(candidate, reference, signed=False):
    candidate = np.ascontiguousarray(candidate, dtype=np.float64)
    reference = np.ascontiguousarray(reference, dtype=np.float64)
    n = candidate.shape[0]
    mask = np.ones((n, n), dtype=np.int32)
    np.fill_diagonal(mask, 0)
    metrics = np.zeros(7, dtype=np.float64)
    counts = np.zeros(3, dtype=np.int64)
    status = oqp.oqp_namd_nacme_gate(
        n,
        oqp.ffi.cast('double *', candidate.ctypes.data),
        oqp.ffi.cast('double *', reference.ctypes.data),
        oqp.ffi.cast('int *', mask.ctypes.data),
        int(signed), 1.0e-12, 1.0e-8, 0.01,
        oqp.ffi.cast('double *', metrics.ctypes.data),
        oqp.ffi.cast('int64_t *', counts.ctypes.data),
    )
    return status, metrics.tolist(), counts.tolist()

result = {
    'uniform': run([0.0, 0.03], [0.0, 0.02], [0.0, 0.03]),
    'nonuniform': run([0.0, 0.06], [0.0, 0.02], [0.0, 0.03], 2.0, 1.0),
    'negative_radicand': run([0.0, 0.01], [0.0, 0.02], [0.0, 0.01]),
    'gap_cutoff': run([0.0, 0.03], [0.0, 0.02], [0.0, 0.03], 1.0, 1.0, 0.01),
    'gate_pass': gate([[0.0, -0.5], [0.5, 0.0]],
                      [[0.0, -0.5], [0.5, 0.0]]),
    'gate_magnitude': gate([[0.0, 0.5], [-0.5, 0.0]],
                           [[0.0, -0.5], [0.5, 0.0]]),
    'gate_signed': gate([[0.0, 0.5], [-0.5, 0.0]],
                        [[0.0, -0.5], [0.5, 0.0]], signed=True),
    'gate_invariant': gate([[0.0, -0.5], [0.4, 0.0]],
                           [[0.0, -0.5], [0.5, 0.0]]),
}
print('BAECK_AN=' + json.dumps(result))
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
        if "No module named" in result.stderr or "cannot load" in result.stderr:
            pytest.skip("compiled OpenQP runtime is not available")
        pytest.fail(result.stdout + result.stderr)
    marker = next(
        (line for line in result.stdout.splitlines() if line.startswith("BAECK_AN=")),
        None,
    )
    assert marker is not None, result.stdout + result.stderr
    values = json.loads(marker.removeprefix("BAECK_AN="))

    for key in ("uniform", "nonuniform"):
        status, matrix = values[key]
        assert status == 0
        np.testing.assert_allclose(matrix, [[0.0, -0.5], [0.5, 0.0]], atol=1e-14)
    for key in ("negative_radicand", "gap_cutoff"):
        status, matrix = values[key]
        assert status == 0
        assert matrix == [[0.0, 0.0], [0.0, 0.0]]

    for key in ('gate_pass', 'gate_magnitude'):
        status, _metrics, counts = values[key]
        assert status == 0
        assert counts == [1, 0, 0]
    status, _metrics, counts = values['gate_signed']
    assert status == 0
    assert counts == [1, 0, 1]
    status, metrics, counts = values['gate_invariant']
    assert status == 0
    assert counts[1] == 1
    assert metrics[1] == pytest.approx(0.1)


def test_dense_trajectory_is_appendable_and_restart_manifest_is_runnable(tmp_path):
    script = r"""
import json
import os
from types import SimpleNamespace
import numpy as np
import oqp.library.namd as namd_module
from oqp.library.namd import (
    NAMD, read_namd_trajectory, _validate_distinct_output_paths,
    _validate_nacme_gate_activation,
)
namd_module.dump_log = lambda *_args, **_kwargs: None

root = os.environ['OQP_NAMD_TEST_ROOT']
try:
    _validate_distinct_output_paths(
        trajectory_file=os.path.join(root, 'collision'),
        nacme_audit_file=os.path.join(root, '.', 'collision'),
    )
except ValueError as error:
    output_collision_rejected = 'must be distinct' in str(error)
else:
    output_collision_rejected = False
inactive_gate_rejected = []
for policy in ('warn', 'error'):
    try:
        _validate_nacme_gate_activation('off', policy)
    except ValueError as error:
        inactive_gate_rejected.append('nacme_check is off' in str(error))
    else:
        inactive_gate_rejected.append(False)
_validate_nacme_gate_activation('off', 'off')
_validate_nacme_gate_activation('baeck_an', 'error')

class Mol:
    log = os.path.join(root, 'job.log')
    oqp_input_source = os.path.join(root, 'source', 'request.oqp')
    oqp_canonical_input = (
        'mrsf(nstate=2)/bhhlyp/6-31g*\n'
        'namd(S1,nstep=20,dt=0.5)\ngeom="h2o.xyz"\n'
    )
    config = {
        'input': {'method': 'tdhf', 'functional': 'bhhlyp', 'basis': '6-31g*',
                  'charge': 0},
        'scf': {'type': 'rohf', 'multiplicity': 3},
        'tdhf': {'type': 'mrsf', 'nstate': 2, 'tlf': 2, 'multiplicity': 3},
    }
    data = {'OQP::td_energies': np.array([0.0, 0.1])}
    @staticmethod
    def get_state_tracking():
        return {'order': [0, 1], 'phase_step': [1.0, -1.0],
                'matched_overlap': [0.99, 0.98], 'margin': [0.2, 0.3]}
    def put_data(self, data):
        self.loaded = data

d = NAMD.__new__(NAMD)
d.mol = Mol(); d.nstate = 2; d.dt_fs = 0.5; d.dt_adaptive = False
d._t_fs = 0.0; d.seed = 1; d.rng_stream = 2; d.trajectory_interval = 1
d.trajectory_file = os.path.join(root, 'job.namd.trj')
d.restart_file = os.path.join(root, 'job.namd.restart.npz')
d.nacme_audit_file = os.path.join(root, 'job.namd.nacme.tsv')
d.restart_manifest_file = os.path.join(root, 'restart.oqp')
d._restart_system_identity = {'kind': 'test', 'sha256': 'system-a'}
d.active = 1; d._last_hop_random = 0.25; d.coef = np.array([1+0j, 0+0j])
d.init_temp = 300.0; d.velocity_source = '/provided/velocity.dat'
d.mass = np.ones(3); d.vel = np.ones((3, 3))*0.01
d._last_state_overlap = np.eye(2)
d._last_overlap_tdc = np.array([[0.0, 0.1], [-0.1, 0.0]])
d._nacme_reference_tdc = np.array([[0.0, 0.11], [-0.11, 0.0]])
d._nacme_reference_mask = np.array([[0, 1], [1, 0]], dtype=np.int32)
d._nacme_reference_source = 1
d._nacme_gate_last = {
    'center_step': 0, 'verdict': 'pass', 'compared_pairs': 1,
    'invariant_failures': 0, 'reference_failures': 0,
    'candidate_diagonal_max': 0.0, 'candidate_antisymmetry_max': 0.0,
    'reference_diagonal_max': 0.0, 'reference_antisymmetry_max': 0.0,
    'pair_rms_error': 0.01, 'pair_max_error': 0.01,
    'max_tolerance_ratio': 0.2,
}
d._pending_nacme_gate_error = None
d.nve_gate = 'warn'; d.nve_gate_abs_tol = 0.005
d.nve_gate_step_tol = 0.001; d.nve_gate_transition_tol = 1.0e-6
d.nve_gate_consecutive = 3; d._nve_reference_energy = None
d._nve_previous_energy = None; d._nve_gate_failures = 0; d._nve_gate_last = None
d._update_nve_gate(0, -1.0, 0.1)
d._update_nve_gate(1, -1.0, 0.1)
d._write_md_trajectory(1, np.zeros((3, 3)), -1.0, 0.1, False)
d.restart_interval = 1; d.prev_xyz = np.arange(9.0)
d.prev_data = {
    'OQP::VEC_MO_A': np.arange(6.0).reshape(2, 3),
    'OQP::td_bvec_mo': np.arange(4.0).reshape(2, 2),
    'OQP::state_tracking_phase_initial': np.array([1.0, -1.0]),
}
d._ba_energy_left = np.array([0.0, 0.2]); d._ba_energy_center = np.array([0.0, 0.1])
d._ba_tdc_left = np.array([[0.0, 0.1], [-0.1, 0.0]]); d._ba_dt_left = 2.0
d._nacme_gate_failures = 2; d._rng_step = 1
d._write_nacme_audit_row({
    'center_step': 0, 'source': 'baeck_an', 'verdict': 'pass',
    'signed_comparison': False, 'compared_pairs': 1,
})
d._save_restart(1, np.zeros((3, 3)), d.vel, np.ones((3, 3))*0.001)
d._update_nve_gate(2, -0.9, 0.1, 0.002)
d._write_md_trajectory(2, np.ones((3, 3)), -0.9, 0.1, True)
d._write_nacme_audit_row({
    'center_step': 1, 'source': 'baeck_an', 'verdict': 'fail',
    'signed_comparison': True, 'compared_pairs': 1,
})

d2 = NAMD.__new__(NAMD); d2.mol = Mol(); d2.nstate = 2; d2.dt_fs = 0.5
d2.seed = 1; d2.rng_stream = 2; d2.restart_requested = True
d2._restart_system_identity = {'kind': 'test', 'sha256': 'system-a'}
d2.restart_file = d.restart_file; d2.trajectory_file = d.trajectory_file
d2.nacme_audit_file = d.nacme_audit_file
loaded = d2._load_restart()
d_missing = NAMD.__new__(NAMD); d_missing.mol = Mol()
d_missing.nstate = 2; d_missing.dt_fs = 0.5
d_missing.seed = 1; d_missing.rng_stream = 2
d_missing.restart_requested = True
d_missing._restart_system_identity = {'kind': 'test', 'sha256': 'system-a'}
d_missing.restart_file = d.restart_file
d_missing.trajectory_file = os.path.join(root, 'missing.namd.trj')
d_missing.nacme_audit_file = d.nacme_audit_file
try:
    d_missing._load_restart()
except FileNotFoundError as error:
    missing_trajectory_rejected = 'trajectory not found' in str(error)
else:
    missing_trajectory_rejected = False
d3 = NAMD.__new__(NAMD); d3.mol = Mol(); d3.nstate = 2; d3.dt_fs = 0.5
d3.seed = 1; d3.rng_stream = 2; d3.restart_requested = True
d3._restart_system_identity = {'kind': 'test', 'sha256': 'system-b'}
d3.restart_file = d.restart_file; d3.trajectory_file = d.trajectory_file
try:
    d3._load_restart()
except ValueError as error:
    system_mismatch = 'does not match the current run' in str(error)
else:
    system_mismatch = False
d4 = NAMD.__new__(NAMD); d4.mol = Mol(); d4.nstate = 2; d4.dt_fs = 0.5
d4.mol.config = json.loads(json.dumps(Mol.config))
d4.mol.config['input']['charge'] = 1
d4.seed = 1; d4.rng_stream = 2; d4.restart_requested = True
d4._restart_system_identity = {'kind': 'test', 'sha256': 'system-a'}
d4.restart_file = d.restart_file; d4.trajectory_file = d.trajectory_file
try:
    d4._load_restart()
except ValueError as error:
    electronic_mismatch = 'does not match the current run' in str(error)
else:
    electronic_mismatch = False
d_d4 = NAMD.__new__(NAMD); d_d4.mol = Mol()
d_d4.mol.config = json.loads(json.dumps(Mol.config))
d_d4.mol.config['input']['d4'] = True
d_d4.nstate = 2; d_d4.dt_fs = 0.5
d_d4.seed = 1; d_d4.rng_stream = 2; d_d4.restart_requested = True
d_d4._restart_system_identity = {'kind': 'test', 'sha256': 'system-a'}
d_d4.restart_file = d.restart_file; d_d4.trajectory_file = d.trajectory_file
try:
    d_d4._load_restart()
except ValueError as error:
    d4_mismatch = 'does not match the current run' in str(error)
else:
    d4_mismatch = False
d_pcm = NAMD.__new__(NAMD); d_pcm.mol = Mol()
d_pcm.mol.config = json.loads(json.dumps(Mol.config))
d_pcm.mol.config['pcm'] = {
    'enabled': True, 'backend': 'ddx', 'mode': 'reference_scf',
    'model': 'ddpcm', 'solvent': 'water', 'epsilon': 40.0, 'radii': 'uff',
}
d_pcm.nstate = 2; d_pcm.dt_fs = 0.5
d_pcm.seed = 1; d_pcm.rng_stream = 2; d_pcm.restart_requested = True
d_pcm._restart_system_identity = {'kind': 'test', 'sha256': 'system-a'}
d_pcm.restart_file = d.restart_file; d_pcm.trajectory_file = d.trajectory_file
try:
    d_pcm._load_restart()
except ValueError as error:
    pcm_mismatch = 'does not match the current run' in str(error)
else:
    pcm_mismatch = False
d_gate = NAMD.__new__(NAMD); d_gate.mol = Mol()
d_gate.mol.config = json.loads(json.dumps(Mol.config))
d_gate.mol.config['md'] = {'nve_gate': 'error'}
d_gate.nstate = 2; d_gate.dt_fs = 0.5
d_gate.seed = 1; d_gate.rng_stream = 2; d_gate.restart_requested = True
d_gate._restart_system_identity = {'kind': 'test', 'sha256': 'system-a'}
d_gate.restart_file = d.restart_file; d_gate.trajectory_file = d.trajectory_file
try:
    d_gate._load_restart()
except ValueError as error:
    gate_mismatch = 'does not match the current run' in str(error)
else:
    gate_mismatch = False
header, records = read_namd_trajectory(d.trajectory_file)
manifest = open(d.restart_manifest_file, encoding='utf-8').read()
manifest_geometry = os.path.realpath(
    os.path.join(root, 'source', 'h2o.xyz')) in manifest
with open(d.nacme_audit_file, encoding='utf-8') as stream:
    audit_lines = [line.rstrip('\n').split('\t') for line in stream]
audit_header, audit_row = audit_lines
audit = dict(zip(audit_header, audit_row))
d.nacme_check = 'baeck_an'; d.dt = 1.0
d._ba_energy_center = np.array([0.0, 0.1])
d._nacme_gate_last = {'verdict': 'fail'}
d._nacme_reference_tdc = np.ones((2, 2))
d._nacme_reference_mask = np.ones((2, 2), dtype=np.int32)
d._nacme_reference_source = 1; d._nacme_gate_failures = 2
d.mol.data['OQP::td_energies_old'] = np.array([0.0, 0.2])
d.mol.data['OQP::td_energies'] = np.array([0.0, 0.3])
d._compute_tdc = lambda overlap: np.asarray(overlap, dtype=float)
d._update_baeck_an_check(3, np.eye(2))
reseed_cleared = (
    d._nacme_gate_last is None
    and d._nacme_reference_tdc is None
    and d._nacme_reference_mask is None
    and d._nacme_reference_source == 0
    and d._nacme_gate_failures == 0
)
d.nacme_gate = 'warn'; d.nacme_gate_invariant_tol = 1.0e-12
d.nacme_gate_abs_tol = 1.0e-8; d.nacme_gate_rel_tol = 0.01
d.nacme_gate_consecutive = 3
invariant_without_pairs = d._run_nacme_gate(
    np.array([[1.0e-3, 0.0], [0.0, 0.0]]),
    np.zeros((2, 2)), reference_mask=np.zeros((2, 2), dtype=np.int32),
    source='analytic', center_step=4, signed=True,
)
d.nve_gate = 'error'; d.nve_gate_consecutive = 1
d._nve_reference_energy = -0.9; d._nve_previous_energy = -0.9
d._nve_gate_failures = 0
d._update_nve_gate(3, -0.8, 0.1)
deferred_error = d._pending_nve_gate_error is not None
d.trajectory_interval = 10
d.trajectory_file = os.path.join(root, 'failure.namd.trj')
d._write_md_trajectory(3, np.zeros((3, 3)), -0.8, 0.1, False)
_, failure_records = read_namd_trajectory(d.trajectory_file)
forced_failure_steps = failure_records['step'].tolist()
try:
    d._enforce_nve_gate()
except RuntimeError:
    enforced_error = True
else:
    enforced_error = False
d.nacme_gate = 'error'; d.nacme_gate_consecutive = 1
d._pending_nacme_gate_error = None
d._run_nacme_gate(
    np.array([[1.0e-3, 0.0], [0.0, 0.0]]),
    np.zeros((2, 2)), reference_mask=np.zeros((2, 2), dtype=np.int32),
    source='analytic', center_step=5, signed=True,
)
pending_nacme_error = d._pending_nacme_gate_error
d.trajectory_file = os.path.join(root, 'failure-nacme.namd.trj')
d._write_md_trajectory(5, np.ones((3, 3)), -0.8, 0.1, False)
_, nacme_failure_records = read_namd_trajectory(d.trajectory_file)
try:
    d._enforce_nacme_gate()
except RuntimeError as error:
    enforced_nacme_error = error is pending_nacme_error
else:
    enforced_nacme_error = False
d.nacme_gate = 'warn'; d._pending_nacme_gate_error = None
nonfinite_gate = d._run_nacme_gate(
    np.array([[np.nan, 0.0], [0.0, 0.0]]),
    np.zeros((2, 2)), reference_mask=np.zeros((2, 2), dtype=np.int32),
    source='analytic', center_step=6, signed=True,
)
pending_nonfinite_error = d._pending_nacme_gate_error
d.trajectory_file = os.path.join(root, 'failure-nonfinite.namd.trj')
d._write_md_trajectory(6, np.ones((3, 3)), -0.8, 0.1, False)
_, nonfinite_failure_records = read_namd_trajectory(d.trajectory_file)
try:
    d._enforce_nacme_gate()
except RuntimeError as error:
    enforced_nonfinite_error = (
        error is pending_nonfinite_error and 'status=-2' in str(error))
else:
    enforced_nonfinite_error = False
print('DENSE=' + json.dumps({
    'shape': records.shape, 'steps': records['step'].tolist(),
    'phase': records['tracking_phase'].tolist(), 'nstate': header['nstate'],
    'natom': header['natom'], 'restart': 'restart=true' in manifest,
    'manifest_geometry': manifest_geometry,
    'checkpoint': 'restart_file="job.namd.restart.npz"' in manifest,
    'loaded_step': loaded['step'],
    'system_mismatch': system_mismatch,
    'electronic_mismatch': electronic_mismatch,
    'd4_mismatch': d4_mismatch,
    'pcm_mismatch': pcm_mismatch,
    'gate_mismatch': gate_mismatch,
    'missing_trajectory_rejected': missing_trajectory_rejected,
    'output_collision_rejected': output_collision_rejected,
    'inactive_gate_rejected': inactive_gate_rejected,
    'reseed_cleared': reseed_cleared,
    'audit_steps': [int(row[0]) for row in audit_lines[1:]],
    'audit_signed_comparison': audit['signed_comparison'],
    'invariant_without_pairs': invariant_without_pairs['verdict'],
    'phase_history': d2.mol.loaded['OQP::state_tracking_phase_initial'].tolist(),
        'gate_failures': d2._nacme_gate_failures,
        'nve_failures': d2._nve_gate_failures,
        'nve_verdict': records['nve_verdict'].tolist(),
        'deferred_error': deferred_error,
        'enforced_error': enforced_error,
        'forced_failure_steps': forced_failure_steps,
        'forced_nacme_failure_steps': nacme_failure_records['step'].tolist(),
        'forced_nacme_verdict': nacme_failure_records['gate_verdict'].tolist(),
        'enforced_nacme_error': enforced_nacme_error,
        'nonfinite_native_status': nonfinite_gate['native_status'],
        'forced_nonfinite_steps': nonfinite_failure_records['step'].tolist(),
        'forced_nonfinite_verdict': nonfinite_failure_records[
            'gate_verdict'].tolist(),
        'forced_nonfinite_candidate_recorded': bool(np.isnan(
            nonfinite_failure_records['overlap_tdc_au'][0, 0, 0])),
        'enforced_nonfinite_error': enforced_nonfinite_error,
        'initial_temperature_kelvin': header['initial_temperature_kelvin'],
        'initial_temperature_dof': header[
            'initial_temperature_degrees_of_freedom'],
        'requested_initial_temperature_kelvin': header[
            'requested_initial_temperature_kelvin'],
        'initial_velocity_source': header['initial_velocity_source'],
}))
"""
    env = os.environ.copy()
    env["OQP_NAMD_TEST_ROOT"] = str(tmp_path)
    pythonpath = str(ROOT / "pyoqp")
    if env.get("PYTHONPATH"):
        pythonpath += os.pathsep + env["PYTHONPATH"]
    env["PYTHONPATH"] = pythonpath
    result = subprocess.run(
        [sys.executable, "-c", script], cwd=ROOT, env=env,
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        if "No module named" in result.stderr or "cannot load" in result.stderr:
            pytest.skip("compiled OpenQP runtime is not available")
        pytest.fail(result.stdout + result.stderr)
    marker = next(
        (line for line in result.stdout.splitlines() if line.startswith("DENSE=")),
        None,
    )
    assert marker is not None, result.stdout + result.stderr
    values = json.loads(marker.removeprefix("DENSE="))
    assert values == {
        'shape': [1], 'steps': [1],
        'phase': [[1.0, -1.0]],
        'nstate': 2, 'natom': 3, 'restart': True,
        'manifest_geometry': True, 'checkpoint': True,
        'loaded_step': 1, 'phase_history': [1.0, -1.0], 'gate_failures': 2,
        'system_mismatch': True, 'audit_signed_comparison': 'False',
        'electronic_mismatch': True, 'reseed_cleared': True,
        'd4_mismatch': True,
        'pcm_mismatch': True,
        'gate_mismatch': True,
        'missing_trajectory_rejected': True,
        'output_collision_rejected': True,
        'inactive_gate_rejected': [True, True], 'audit_steps': [0],
        'invariant_without_pairs': 'fail',
        'nve_failures': 0, 'nve_verdict': [1],
        'deferred_error': True, 'enforced_error': True,
        'forced_failure_steps': [3],
        'forced_nacme_failure_steps': [5],
        'forced_nacme_verdict': [2],
        'enforced_nacme_error': True,
        'nonfinite_native_status': -2,
        'forced_nonfinite_steps': [6],
        'forced_nonfinite_verdict': [2],
        'forced_nonfinite_candidate_recorded': True,
        'enforced_nonfinite_error': True,
        'initial_temperature_kelvin': pytest.approx(
            9.0e-4/(6*3.166811563e-6)),
        'initial_temperature_dof': 6,
        'requested_initial_temperature_kelvin': 300.0,
        'initial_velocity_source': 'file',
    }
