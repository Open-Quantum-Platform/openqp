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
    example_inp = (
        ROOT / "examples" / "QMMM" /
        "H2CO-water_BHHLYP-MRSF-NAMD-QMMM.inp"
    ).read_text()
    example_oqp = (
        ROOT / "examples" / "QMMM" /
        "H2CO-water_BHHLYP-MRSF-NAMD-QMMM.oqp"
    ).read_text()
    assert "nstep=2" in example_inp
    assert "nstep=2" in example_oqp
    assert "nacme_check=baeck_an" in example_inp
    assert "nacme_check=baeck_an" in example_oqp


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


def test_soc_namd_rejects_unimplemented_nve_and_dense_trajectory_controls(tmp_path):
    script = r"""
import json
import os
from oqp.library.namd import NAMD

root = os.environ['OQP_NAMD_TEST_ROOT']
base_md = {
    'nstep': 2, 'dt': 0.5, 'active': 1, 'substep': 10,
    'decoherence': 'edc', 'edc_c': 0.1, 'thrshe': 0.01, 'tdc': 'fd',
    'trivial': False, 'trivial_thresh': 0.1, 'init_temp': 0.0,
    'velocity': 'zero', 'seed': 1, 'soc': True,
}

class Mol:
    log = os.path.join(root, 'soc.log')
    def __init__(self, **overrides):
        md = dict(base_md)
        md.update(overrides)
        self.config = {
            'input': {'method': 'tdhf', 'functional': 'bhhlyp',
                      'basis': '6-31g*'},
            'tdhf': {'type': 'mrsf', 'nstate': 2, 'tlf': 2},
            'properties': {}, 'nac': {}, 'md': md,
        }

def rejected(**overrides):
    try:
        NAMD(Mol(**overrides))
    except NotImplementedError:
        return True
    return False

print('SOC_GATES=' + json.dumps({
    'nve': rejected(nve_gate='warn'),
    'trajectory_file': rejected(trajectory_file='soc.trj'),
    'trajectory_interval': rejected(trajectory_interval=2),
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
        (line for line in result.stdout.splitlines()
         if line.startswith("SOC_GATES=")),
        None,
    )
    assert marker is not None, result.stdout + result.stderr
    assert json.loads(marker.removeprefix("SOC_GATES=")) == {
        'nve': True, 'trajectory_file': True, 'trajectory_interval': True,
    }


def test_dense_trajectory_is_appendable_and_restart_manifest_is_runnable(tmp_path):
    script = r"""
import json
import os
from types import SimpleNamespace
import numpy as np
import oqp.library.namd as namd_module
from oqp.library.namd import NAMD, read_namd_trajectory
namd_module.dump_log = lambda *_args, **_kwargs: None

root = os.environ['OQP_NAMD_TEST_ROOT']
input_root = os.path.join(root, 'input')
os.makedirs(input_root)
for filename in ('h2o.xyz', 'vel.dat', 'water.pdb', 'local.xml'):
    with open(os.path.join(input_root, filename), 'w', encoding='utf-8') as stream:
        stream.write('test\n')
class Mol:
    log = os.path.join(root, 'job.log')
    oqp_input_source = os.path.join(input_root, 'request.oqp')
    oqp_canonical_input = (
        'mrsf(nstate=2)/bhhlyp/6-31g*\n'
        'namd(S1,nstep=20,dt=0.5,velocity="vel.dat")\n'
        'qmmm(pdb_file="water.pdb",forcefield_files="local.xml tip3p.xml",'
        'qm_atoms="0-2")\ngeom="h2o.xyz"\n'
    )
    config = {
        'input': {'method': 'tdhf', 'functional': 'bhhlyp', 'basis': '6-31g*'},
        'tdhf': {'type': 'mrsf', 'nstate': 2, 'tlf': 2},
    }
    data = {'OQP::td_energies': np.array([0.0, 0.1])}
    atoms = np.array([8, 1, 1])
    masses = np.array([15.999, 1.008, 1.008])
    def get_atoms(self):
        return self.atoms.copy()
    def get_mass(self):
        return self.masses.copy()
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
d.restart_manifest_file = d._restart_manifest_path()
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
d._save_restart(1, np.zeros((3, 3)), d.vel, np.ones((3, 3))*0.001)

# Rebased local force-field paths must retain the same restart identity, while
# OpenMM built-in names remain stable symbolic identities.
relative_ff_identity = d._qmmm_forcefield_identity('local.xml tip3p.xml')
absolute_ff_identity = d._qmmm_forcefield_identity(
    os.path.join(input_root, 'local.xml') + ' tip3p.xml')
forcefield_identity_stable = relative_ff_identity == absolute_ff_identity

# Electronic arrays must retain their exact vector shape, and tracking
# histories must remain one entry per state.
try:
    d._validate_restart_state(
        (np.zeros((3, 3)), d.vel, np.ones((3, 3))*0.001),
        np.array([[1+0j], [0+0j]]), d.active, d.prev_xyz, d.prev_data,
        context='shape test')
except RuntimeError:
    coefficient_shape_rejected = True
else:
    coefficient_shape_rejected = False
bad_tracking = dict(d.prev_data)
bad_tracking['OQP::state_tracking_phase_initial'] = np.array([1.0])
try:
    d._validate_restart_state(
        (np.zeros((3, 3)), d.vel, np.ones((3, 3))*0.001),
        d.coef, d.active, d.prev_xyz, bad_tracking, context='tracking test')
except RuntimeError:
    tracking_shape_rejected = True
else:
    tracking_shape_rejected = False

# A non-finite electronic state must not replace the last-good checkpoint.
d.coef = np.array([np.nan+0j, 0+0j])
try:
    d._save_restart(2, np.zeros((3, 3)), d.vel, np.ones((3, 3))*0.001)
except RuntimeError:
    invalid_save_rejected = True
else:
    invalid_save_rejected = False
d.coef = np.array([1+0j, 0+0j])
with np.load(d.restart_file, allow_pickle=False) as saved:
    checkpoint_after_rejection = int(saved['step'][0])

# The audit has one committed row, one uncommitted row, and a torn final row.
audit_base = {
    'source': 'TD-Baeck-An', 'verdict': 'pass', 'signed_comparison': False,
    'compared_pairs': 1, 'invariant_failures': 0, 'reference_failures': 0,
    'consecutive_reference_failures': 0,
}
d._write_nacme_audit_row(dict(audit_base, center_step=1))
d._write_nacme_audit_row(dict(audit_base, center_step=2,
                              signed_comparison=True))
with open(d.nacme_audit_file, 'a', encoding='utf-8') as stream:
    stream.write('3\tpartial')

# Simulate interruption in the middle of the final packed trajectory record.
with open(d.trajectory_file, 'ab') as stream:
    stream.write(b'partial-record')

d2 = NAMD.__new__(NAMD); d2.mol = Mol(); d2.nstate = 2; d2.dt_fs = 0.5
d2.seed = 1; d2.rng_stream = 2; d2.restart_requested = True
d2.restart_file = d.restart_file; d2.trajectory_file = d.trajectory_file
d2.nacme_audit_file = d.nacme_audit_file
loaded = d2._load_restart()
header, records = read_namd_trajectory(d.trajectory_file)
manifest = open(d.restart_manifest_file, encoding='utf-8').read()
audit_lines = open(d.nacme_audit_file, encoding='utf-8').read().splitlines()

# Externally corrupted electronic payloads are rejected during load too.
with np.load(d.restart_file, allow_pickle=False) as saved:
    corrupt = {key: np.array(saved[key], copy=True) for key in saved.files}
corrupt['coef_real'][0] = np.nan
bad_restart = os.path.join(root, 'bad.namd.restart.npz')
np.savez_compressed(bad_restart, **corrupt)
d_bad = NAMD.__new__(NAMD); d_bad.mol = Mol(); d_bad.nstate = 2
d_bad.dt_fs = 0.5; d_bad.seed = 1; d_bad.rng_stream = 2
d_bad.restart_requested = True; d_bad.restart_file = bad_restart
d_bad.trajectory_file = d.trajectory_file
d_bad.nacme_audit_file = d.nacme_audit_file
try:
    d_bad._load_restart()
except RuntimeError:
    corrupt_load_rejected = True
else:
    corrupt_load_rejected = False

# Same method/count but different atom identities must not accept the file.
wrong_mol = Mol(); wrong_mol.atoms = np.array([7, 1, 2])
d_wrong = NAMD.__new__(NAMD); d_wrong.mol = wrong_mol; d_wrong.nstate = 2
d_wrong.dt_fs = 0.5; d_wrong.seed = 1; d_wrong.rng_stream = 2
d_wrong.restart_requested = True; d_wrong.restart_file = d.restart_file
d_wrong.trajectory_file = d.trajectory_file
d_wrong.nacme_audit_file = d.nacme_audit_file
try:
    d_wrong._load_restart()
except ValueError:
    molecule_mismatch_rejected = True
else:
    molecule_mismatch_rejected = False

# QM selections participate in the signature even with identical QM atoms.
d_qm1 = NAMD.__new__(NAMD); d_qm1.mol = Mol(); d_qm1.nstate = 2
d_qm1.dt_fs = 0.5; d_qm1.seed = 1; d_qm1.rng_stream = 2
d_qm1.qm_atoms = np.array([0, 1])
d_qm2 = NAMD.__new__(NAMD); d_qm2.mol = Mol(); d_qm2.nstate = 2
d_qm2.dt_fs = 0.5; d_qm2.seed = 1; d_qm2.rng_stream = 2
d_qm2.qm_atoms = np.array([1, 2])
qm_selection_bound = d_qm1._restart_signature() != d_qm2._restart_signature()

# Defaults derive a unique manifest from the log stem.
d_other = NAMD.__new__(NAMD)
d_other.mol = SimpleNamespace(log=os.path.join(root, 'other.log'))
unique_manifests = d._restart_manifest_path() != d_other._restart_manifest_path()
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
    'checkpoint': 'restart_file="job.namd.restart.npz"' in manifest,
    'manifest_name': os.path.basename(d.restart_manifest_file),
    'geom_rebased': os.path.join(input_root, 'h2o.xyz') in manifest,
    'velocity_rebased': os.path.join(input_root, 'vel.dat') in manifest,
    'pdb_rebased': os.path.join(input_root, 'water.pdb') in manifest,
    'forcefield_rebased': os.path.join(input_root, 'local.xml') in manifest,
    'builtin_forcefield_preserved': 'tip3p.xml' in manifest,
    'loaded_step': loaded['step'],
    'phase_history': d2.mol.loaded['OQP::state_tracking_phase_initial'].tolist(),
        'gate_failures': d2._nacme_gate_failures,
        'nve_failures': d2._nve_gate_failures,
        'nve_verdict': records['nve_verdict'].tolist(),
        'deferred_error': deferred_error,
        'enforced_error': enforced_error,
        'forced_failure_steps': forced_failure_steps,
        'invalid_save_rejected': invalid_save_rejected,
        'checkpoint_after_rejection': checkpoint_after_rejection,
        'corrupt_load_rejected': corrupt_load_rejected,
        'molecule_mismatch_rejected': molecule_mismatch_rejected,
        'qm_selection_bound': qm_selection_bound,
        'unique_manifests': unique_manifests,
        'forcefield_identity_stable': forcefield_identity_stable,
        'coefficient_shape_rejected': coefficient_shape_rejected,
        'tracking_shape_rejected': tracking_shape_rejected,
        'audit_rows': audit_lines,
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
        'manifest_name': 'job.namd.restart.oqp',
        'geom_rebased': True, 'velocity_rebased': True,
        'pdb_rebased': True, 'forcefield_rebased': True,
        'builtin_forcefield_preserved': True,
        'loaded_step': 1, 'phase_history': [1.0, -1.0], 'gate_failures': 2,
        'nve_failures': 1, 'nve_verdict': [1],
        'deferred_error': True, 'enforced_error': True,
        'forced_failure_steps': [3],
        'invalid_save_rejected': True, 'checkpoint_after_rejection': 1,
        'corrupt_load_rejected': True, 'molecule_mismatch_rejected': True,
        'qm_selection_bound': True, 'unique_manifests': True,
        'forcefield_identity_stable': True,
        'coefficient_shape_rejected': True, 'tracking_shape_rejected': True,
        'audit_rows': [
            'center_step\tsource\tverdict\tsigned\tcompared_pairs\t'
            'invariant_failures\treference_failures\t'
            'consecutive_reference_failures\tcandidate_diagonal_max\t'
            'candidate_antisymmetry_max\treference_diagonal_max\t'
            'reference_antisymmetry_max\tpair_rms_error\tpair_max_error\t'
            'max_tolerance_ratio',
            '1\tTD-Baeck-An\tpass\tFalse\t1\t0\t0\t0\t\t\t\t\t\t\t',
        ],
    }
