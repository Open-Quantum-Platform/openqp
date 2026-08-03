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
    NAMD, read_namd_trajectory, _validate_gate_tolerances,
    _validate_namd_sidecar_paths, _validate_soc_control_compatibility,
    _validate_thermostat_parameters,
)
namd_module.dump_log = lambda *_args, **_kwargs: None

root = os.environ['OQP_NAMD_TEST_ROOT']
class Mol:
    log = os.path.join(root, 'job.log')
    oqp_canonical_input = (
        'mrsf(nstate=2)/bhhlyp/6-31g*\n'
        'namd(S1,nstep=20,dt=0.5)\ngeom="h2o.xyz"\n'
    )
    config = {
        'input': {'method': 'tdhf', 'functional': 'bhhlyp', 'basis': '6-31g*'},
        'tdhf': {'type': 'mrsf', 'nstate': 2, 'tlf': 2},
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
d.active = 1; d._last_hop_random = 0.25; d.coef = np.array([1+0j, 0+0j])
d.vel = np.ones((3, 3))*0.01; d._last_state_overlap = np.eye(2)
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
d.nve_gate = 'warn'; d.nve_gate_abs_tol = 0.005
d.nve_gate_step_tol = 0.001; d.nve_gate_transition_tol = 1.0e-6
d.nve_gate_consecutive = 3; d._nve_reference_energy = None
d._nve_previous_energy = None; d._nve_gate_failures = 0; d._nve_gate_last = None
d._droplet_energy = 0.02; d._droplet_max_penetration = 0.3
d._droplet_active_count = 2; d._solute_com_energy = 0.004
d._solute_com_displacement = 0.15; d._conservative_restraint_energy = 0.024
d._thermostat_exchange = 0.001; d._thermostat_exchange_cumulative = 0.003
d._update_nve_gate(0, -1.0, 0.1)
d._update_nve_gate(1, -1.0, 0.1)
d._write_md_trajectory(1, np.zeros((3, 3)), -1.0, 0.1, False)
d._update_nve_gate(2, -0.9, 0.1, 0.002)
d._write_md_trajectory(2, np.ones((3, 3)), -0.9, 0.1, True)
d.restart_interval = 1; d.prev_xyz = np.arange(9.0)
d.prev_data = {
    'OQP::VEC_MO_A': np.arange(6.0).reshape(2, 3),
    'OQP::td_bvec_mo': np.arange(4.0).reshape(2, 2),
    'OQP::state_tracking_phase_initial': np.array([1.0, -1.0]),
}
d._ba_energy_left = np.array([0.0, 0.2]); d._ba_energy_center = np.array([0.0, 0.1])
d._ba_tdc_left = np.array([[0.0, 0.1], [-0.1, 0.0]]); d._ba_dt_left = 2.0
d._nacme_gate_failures = 2; d._rng_step = 1
audit = {
    'center_step': 0, 'source': 'analytic', 'verdict': 'pass',
    'signed_comparison': True, 'compared_pairs': 1,
    'invariant_failures': 0, 'reference_failures': 0,
}
d._write_nacme_audit_row(dict(audit, evaluation_step=1))
d._write_nacme_audit_row(dict(audit, evaluation_step=2, center_step=1,
                              signed_comparison=False))
d._save_restart(1, np.zeros((3, 3)), d.vel, np.ones((3, 3))*0.001)
last_good_step = 1
d.coef[0] = np.nan + 0j
try:
    d._save_restart(2, np.zeros((3, 3)), d.vel, np.ones((3, 3))*0.001)
except RuntimeError:
    invalid_electronic_checkpoint_rejected = True
else:
    invalid_electronic_checkpoint_rejected = False
with np.load(d.restart_file, allow_pickle=False) as last_good:
    last_good_step = int(last_good['step'][0])
d.coef = np.array([1+0j, 0+0j])

d2 = NAMD.__new__(NAMD); d2.mol = Mol(); d2.nstate = 2; d2.dt_fs = 0.5
d2.seed = 1; d2.rng_stream = 2; d2.restart_requested = True
d2.restart_file = d.restart_file; d2.trajectory_file = d.trajectory_file
d2.nacme_audit_file = d.nacme_audit_file
loaded = d2._load_restart()
header, records = read_namd_trajectory(d.trajectory_file)
manifest = open(d.restart_manifest_file, encoding='utf-8').read()
audit_lines = open(d.nacme_audit_file, encoding='utf-8').read().splitlines()
audit_columns = audit_lines[0].split('\t')
audit_values = audit_lines[1].split('\t')
audit_row = dict(zip(audit_columns, audit_values))

pdb_path = os.path.join(root, 'identity.pdb')
ff_path = os.path.join(root, 'identity.xml')
with open(pdb_path, 'w', encoding='utf-8') as stream:
    stream.write('MODEL 1\nENDMDL\n')
with open(ff_path, 'w', encoding='utf-8') as stream:
    stream.write('<ForceField/>\n')
class QMol:
    input_file = os.path.join(root, 'identity.oqp')
    config = {
        'input': {'method': 'tdhf', 'functional': 'bhhlyp', 'basis': '6-31g*',
                  'qmmm_flag': True},
        'tdhf': {'type': 'mrsf', 'nstate': 2, 'tlf': 2},
        'md': {},
        'qmmm': {'pdb_file': pdb_path, 'forcefield_files': ff_path,
                 'qm_atoms': '0-3', 'cutoff': 'NoCutoff',
                 'embedding': 'electrostatic', 'rigidwater': False,
                 'frontier_scheme': 'none'},
    }
q = NAMD.__new__(NAMD); q.mol = QMol(); q.dt_fs = 0.5
q.seed = 1; q.rng_stream = 2
signature_before = q._restart_signature()
with open(pdb_path, 'a', encoding='utf-8') as stream:
    stream.write('REMARK topology changed\n')
q._qmmm_restart_identity_cache = None
signature_after = q._restart_signature()

class SystemMol:
    config = {
        'input': {'method': 'tdhf', 'functional': 'bhhlyp', 'basis': '6-31g*',
                  'charge': 0, 'system': 'H 0 0 0\nH 0 0 0.7'},
        'scf': {'multiplicity': 1},
        'tdhf': {'type': 'mrsf', 'nstate': 2, 'tlf': 2},
        'md': {},
    }
    @staticmethod
    def get_atoms():
        return np.array([1, 1])
    @staticmethod
    def get_mass():
        return np.array([1.00784, 1.00784])
system_driver = NAMD.__new__(NAMD); system_driver.mol = SystemMol()
system_driver.dt_fs = 0.5; system_driver.seed = 1; system_driver.rng_stream = 2
system_signature_before = system_driver._restart_signature()
SystemMol.config['input']['charge'] = 1
system_driver._molecular_restart_identity_cache = None
system_signature_after = system_driver._restart_signature()

rejected_tolerances = []
for name, values in (
        ('nan', (1.0e-10, np.nan, 1.0)),
        ('inf', (np.inf, 1.0e-4, 1.0))):
    try:
        _validate_gate_tolerances('test gate', *values)
    except ValueError:
        rejected_tolerances.append(name)
try:
    _validate_thermostat_parameters(300.0, 0.0, True)
except ValueError:
    zero_friction_rejected = True
else:
    zero_friction_rejected = False
soc_rejections = []
for name, nve_gate, dense in (
        ('nve', 'warn', (False,)),
        ('dense', 'off', (True,))):
    try:
        _validate_soc_control_compatibility(True, nve_gate, dense)
    except NotImplementedError:
        soc_rejections.append(name)
try:
    _validate_namd_sidecar_paths(
        trajectory_file=os.path.join(root, 'collision.dat'),
        nacme_audit_file=os.path.join(root, '.', 'collision.dat'),
        restart_file=os.path.join(root, 'checkpoint.npz'),
    )
except ValueError:
    sidecar_collision_rejected = True
else:
    sidecar_collision_rejected = False

asset_dir = os.path.join(root, 'assets')
output_dir = os.path.join(root, 'outputs')
os.makedirs(asset_dir); os.makedirs(output_dir)
geom_path = os.path.join(asset_dir, 'geom.xyz')
topology_path = os.path.join(asset_dir, 'system.pdb')
local_ff_path = os.path.join(asset_dir, 'local.xml')
for path, contents in (
        (geom_path, '1\n\nH 0 0 0\n'),
        (topology_path, 'MODEL 1\nENDMDL\n'),
        (local_ff_path, '<ForceField/>\n')):
    with open(path, 'w', encoding='utf-8') as stream:
        stream.write(contents)
source_input = os.path.join(asset_dir, 'source.oqp')
canonical_manifest_input = (
    'mrsf(nstate=2)/bhhlyp/6-31g*\n'
    'namd(T0,nstep=2,dt=0.5)\n'
    'qmmm(pdb_file="system.pdb",forcefield_files="local.xml amber14-all.xml",'
    'qm_atoms="0-3",cutoff=NoCutoff)\n'
    'geom="geom.xyz"\n'
)
with open(source_input, 'w', encoding='utf-8') as stream:
    stream.write(canonical_manifest_input)
class ManifestMol:
    log = os.path.join(output_dir, 'job.log')
    oqp_input_source = source_input
    oqp_canonical_input = canonical_manifest_input
manifest_driver = NAMD.__new__(NAMD); manifest_driver.mol = ManifestMol()
manifest_driver.restart_file = os.path.join(output_dir, 'job.restart.npz')
manifest_driver.trajectory_file = os.path.join(output_dir, 'job.trj')
manifest_driver.nacme_audit_file = os.path.join(output_dir, 'job.tsv')
manifest_driver.restart_manifest_file = os.path.join(output_dir, 'restart.oqp')
manifest_driver._write_restart_manifest()
rebased_manifest = open(
    manifest_driver.restart_manifest_file, encoding='utf-8').read()
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
print('DENSE=' + json.dumps({
    'shape': records.shape, 'steps': records['step'].tolist(),
    'phase': records['tracking_phase'].tolist(), 'nstate': header['nstate'],
    'natom': header['natom'], 'restart': 'restart=true' in manifest,
    'controls': sorted(header['independent_controls']),
    'droplet_energy': records['droplet_energy_hartree'].tolist(),
    'droplet_penetration': records['droplet_max_penetration_bohr'].tolist(),
    'droplet_count': records['droplet_active_count'].tolist(),
    'thermostat_exchange': records['thermostat_exchange_hartree'].tolist(),
    'loaded_droplet': d2._droplet_energy,
    'loaded_thermostat_cumulative': d2._thermostat_exchange_cumulative,
    'checkpoint': 'restart_file="job.namd.restart.npz"' in manifest,
    'audit_rows': len(audit_lines) - 1,
    'audit_signed': audit_row['signed'],
    'qmmm_signature_changed': signature_before != signature_after,
    'system_signature_changed': system_signature_before != system_signature_after,
    'rejected_tolerances': rejected_tolerances,
    'zero_friction_rejected': zero_friction_rejected,
    'soc_rejections': soc_rejections,
    'invalid_electronic_checkpoint_rejected': invalid_electronic_checkpoint_rejected,
    'last_good_step': last_good_step,
    'sidecar_collision_rejected': sidecar_collision_rejected,
    'manifest_geom_rebased': geom_path in rebased_manifest,
    'manifest_topology_rebased': topology_path in rebased_manifest,
    'manifest_forcefield_rebased': local_ff_path in rebased_manifest,
    'manifest_builtin_forcefield': 'amber14-all.xml' in rebased_manifest,
    'loaded_step': loaded['step'],
    'phase_history': d2.mol.loaded['OQP::state_tracking_phase_initial'].tolist(),
        'gate_failures': d2._nacme_gate_failures,
        'nve_failures': d2._nve_gate_failures,
        'nve_verdict': records['nve_verdict'].tolist(),
        'deferred_error': deferred_error,
        'enforced_error': enforced_error,
        'forced_failure_steps': forced_failure_steps,
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
        'nstate': 2, 'natom': 3, 'restart': True, 'checkpoint': True,
        'controls': ['droplet', 'solute_com', 'thermostat'],
        'droplet_energy': [0.02], 'droplet_penetration': [0.3],
        'droplet_count': [2], 'thermostat_exchange': [0.001],
        'loaded_droplet': 0.02, 'loaded_thermostat_cumulative': 0.003,
        'audit_rows': 1, 'audit_signed': 'True',
        'qmmm_signature_changed': True,
        'system_signature_changed': True,
        'rejected_tolerances': ['nan', 'inf'],
        'zero_friction_rejected': True,
        'soc_rejections': ['nve', 'dense'],
        'invalid_electronic_checkpoint_rejected': True,
        'last_good_step': 1,
        'sidecar_collision_rejected': True,
        'manifest_geom_rebased': True, 'manifest_topology_rebased': True,
        'manifest_forcefield_rebased': True,
        'manifest_builtin_forcefield': True,
        'loaded_step': 1, 'phase_history': [1.0, -1.0], 'gate_failures': 2,
        'nve_failures': 1, 'nve_verdict': [1],
        'deferred_error': True, 'enforced_error': True,
        'forced_failure_steps': [3],
    }
