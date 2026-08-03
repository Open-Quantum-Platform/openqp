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
    tracer = (ROOT / "tools" / "diagnostics" / "trace_namd_hop.py").read_text()

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
    assert "if self.restart_requested else self._init_velocities()" in driver
    assert "istep, self.r_all, epot, ekin, hopped" in driver
    assert 'for state in range(nstate)' in tracer
    assert 'self.coef[2]' not in tracer
    assert "does not support SOC NAMD logging paths" in tracer
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
    for keyword in (
            "trajectory_interval=2", "restart_interval=2",
            "trajectory_file=", "nacme_audit_file=", "restart_file="):
        assert keyword in example_inp
        assert keyword in example_oqp
    for keyword in (
            "nacme_gate_invariant_tol=", "nacme_gate_abs_tol=",
            "nacme_gate_rel_tol=", "nacme_gate_consecutive=",
            "nve_gate_abs_tol=", "nve_gate_step_tol=",
            "nve_gate_transition_tol=", "nve_gate_consecutive="):
        assert keyword in example_inp
        assert keyword in example_oqp


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

def invalid(**overrides):
    try:
        NAMD(Mol(**overrides))
    except ValueError:
        return True
    return False

print('SOC_GATES=' + json.dumps({
    'nve': rejected(nve_gate='warn'),
    'trajectory_file': rejected(trajectory_file='soc.trj'),
    'trajectory_interval': rejected(trajectory_interval=2),
    'ba_gap_nan': invalid(ba_gap_max=float('nan')),
    'nacme_nan': invalid(nacme_gate_abs_tol=float('nan')),
    'nacme_inf': invalid(nacme_gate_rel_tol=float('inf')),
    'nve_nan': invalid(nve_gate_abs_tol=float('nan')),
    'nve_inf': invalid(nve_gate_step_tol=float('inf')),
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
        'ba_gap_nan': True, 'nacme_nan': True, 'nacme_inf': True,
        'nve_nan': True, 'nve_inf': True,
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
os.makedirs(os.path.join(input_root, 'parameters'))
with open(os.path.join(input_root, 'parameters', 'H-H.skf'),
          'w', encoding='utf-8') as stream:
    stream.write('parameter\n')
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
        'input': {'method': 'tdhf', 'functional': 'bhhlyp', 'basis': '6-31g*',
                  'charge': 0, 'multiplicity': 3},
        'scf': {'type': 'rohf', 'multiplicity': 3},
        'tdhf': {'type': 'mrsf', 'nstate': 2, 'tlf': 2, 'multiplicity': 3},
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
                'raw_order': [1, 0], 'lineage': [7, 4],
                'phase_initial': [-1.0, 1.0],
                'previous_phase_initial': [1.0, 1.0],
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
    'OQP::state_tracking_lineage': np.array([7, 4]),
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
builtin_forcefield_fingerprinted = (
    relative_ff_identity[1].get('builtin') == 'tip3p.xml'
    and len(relative_ff_identity[1].get('sha256', '')) == 64)

# Tight-binding model settings and parameter contents define the Hamiltonian.
tb1 = Mol(); tb1.config = {
    section: dict(settings) for section, settings in Mol.config.items()}
tb1.config['input']['method'] = 'dftb'
tb1.config['dftb'] = {'backend': 'native', 'model': 'mio',
                      'parameter_path': 'local.xml',
                      'library_path': 'local.xml', 'scc_mixing': 0.35}
tb2 = Mol(); tb2.config = {
    section: dict(settings) for section, settings in tb1.config.items()}
tb2.config['dftb']['model'] = '3ob'
d_tb1 = NAMD.__new__(NAMD); d_tb1.mol = tb1; d_tb1.nstate = 2
d_tb1.dt_fs = 0.5; d_tb1.seed = 1; d_tb1.rng_stream = 2
d_tb2 = NAMD.__new__(NAMD); d_tb2.mol = tb2; d_tb2.nstate = 2
d_tb2.dt_fs = 0.5; d_tb2.seed = 1; d_tb2.rng_stream = 2
tight_binding_bound = d_tb1._restart_signature() != d_tb2._restart_signature()

# Empty TB paths must fingerprint the environment/wheel artifacts that the
# runtime resolver actually selects, not merely the empty input strings.
os.environ['OPENQP_DFTB_PARAMETER_PATH'] = os.path.join(input_root, 'local.xml')
os.environ['OPENQP_DFTB_LIBRARY'] = os.path.join(input_root, 'local.xml')
tb_default = Mol(); tb_default.config = {
    section: dict(settings) for section, settings in Mol.config.items()}
tb_default.config['input']['method'] = 'dftb'
tb_default.config['dftb'] = {'backend': 'native', 'model': 'mio'}
d_tb_default1 = NAMD.__new__(NAMD); d_tb_default1.mol = tb_default
d_tb_default1.nstate = 2; d_tb_default1.dt_fs = 0.5
d_tb_default1.seed = 1; d_tb_default1.rng_stream = 2
default_tb_signature1 = d_tb_default1._restart_signature()
with open(os.path.join(input_root, 'local.xml'), 'a', encoding='utf-8') as stream:
    stream.write('changed\n')
d_tb_default2 = NAMD.__new__(NAMD); d_tb_default2.mol = tb_default
d_tb_default2.nstate = 2; d_tb_default2.dt_fs = 0.5
d_tb_default2.seed = 1; d_tb_default2.rng_stream = 2
default_tb_artifact_bound = (
    default_tb_signature1 != d_tb_default2._restart_signature())
os.environ.pop('OPENQP_DFTB_PARAMETER_PATH')
os.environ.pop('OPENQP_DFTB_LIBRARY')

# The pre-calculation empty DFTB model and the runtime-materialized default
# must describe one restart identity.
tb_model_default = Mol(); tb_model_default.config = {
    section: dict(settings) for section, settings in Mol.config.items()}
tb_model_default.config['input']['method'] = 'dftb'
tb_model_default.config['dftb'] = {
    'backend': 'native', 'type': 'mrsf', 'model': '',
    'parameter_path': 'local.xml', 'library_path': 'local.xml'}
d_tb_model_default = NAMD.__new__(NAMD); d_tb_model_default.mol = tb_model_default
d_tb_model_default.nstate = 2; d_tb_model_default.dt_fs = 0.5
d_tb_model_default.seed = 1; d_tb_model_default.rng_stream = 2
default_model_signature = d_tb_model_default._restart_signature()
default_model_not_mutated = tb_model_default.config['dftb']['model'] == ''
tb_model_resolved = Mol(); tb_model_resolved.config = {
    section: dict(settings) for section, settings in tb_model_default.config.items()}
tb_model_resolved.config['dftb']['model'] = 'dtcam'
d_tb_model_resolved = NAMD.__new__(NAMD); d_tb_model_resolved.mol = tb_model_resolved
d_tb_model_resolved.nstate = 2; d_tb_model_resolved.dt_fs = 0.5
d_tb_model_resolved.seed = 1; d_tb_model_resolved.rng_stream = 2
default_tb_model_stable = (
    default_model_not_mutated
    and default_model_signature == d_tb_model_resolved._restart_signature())

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

# For checkpoint k, only audits centered before k are committed: the
# center_step=k row requires energies from step k+1 and must be regenerated.
audit_base = {
    'source': 'TD-Baeck-An', 'verdict': 'pass', 'signed_comparison': False,
    'compared_pairs': 1, 'invariant_failures': 0, 'reference_failures': 0,
    'consecutive_reference_failures': 0,
}
d._write_nacme_audit_row(dict(audit_base, center_step=0))
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

# Independent real/imaginary vectors must not exploit NumPy broadcasting to
# synthesize an apparently valid nstate-vector during checkpoint loading.
with np.load(d.restart_file, allow_pickle=False) as saved:
    broadcast_corrupt = {
        key: np.array(saved[key], copy=True) for key in saved.files}
broadcast_corrupt['coef_real'] = np.array([0.0])
broadcast_corrupt['coef_imag'] = np.array([1.0, 0.0])
broadcast_restart = os.path.join(root, 'broadcast.namd.restart.npz')
np.savez_compressed(broadcast_restart, **broadcast_corrupt)
d_broadcast = NAMD.__new__(NAMD); d_broadcast.mol = Mol()
d_broadcast.nstate = 2; d_broadcast.dt_fs = 0.5
d_broadcast.seed = 1; d_broadcast.rng_stream = 2
d_broadcast.restart_requested = True
d_broadcast.restart_file = broadcast_restart
d_broadcast.trajectory_file = d.trajectory_file
d_broadcast.nacme_audit_file = d.nacme_audit_file
try:
    d_broadcast._load_restart()
except RuntimeError:
    broadcast_load_rejected = True
else:
    broadcast_load_rejected = False

# Malformed optional native histories must be rejected before restoration.
with np.load(d.restart_file, allow_pickle=False) as saved:
    history_corrupt = {
        key: np.array(saved[key], copy=True) for key in saved.files}
history_corrupt['ba_energy_left'] = np.array([0.0])
history_restart = os.path.join(root, 'history.namd.restart.npz')
np.savez_compressed(history_restart, **history_corrupt)
d_history = NAMD.__new__(NAMD); d_history.mol = Mol(); d_history.nstate = 2
d_history.dt_fs = 0.5; d_history.seed = 1; d_history.rng_stream = 2
d_history.restart_requested = True; d_history.restart_file = history_restart
d_history.trajectory_file = d.trajectory_file
d_history.nacme_audit_file = d.nacme_audit_file
try:
    d_history._load_restart()
except RuntimeError:
    history_load_rejected = True
else:
    history_load_rejected = False

# Gate streaks are control-flow state: malformed, negative, or impossible
# values must not silently disable or prematurely trigger a restarted gate.
streak_rejections = []
for index, replacement in enumerate((
        np.array([-1], dtype=np.int64),
        np.array([[1]], dtype=np.int64),
        np.array([99], dtype=np.int64))):
    with np.load(d.restart_file, allow_pickle=False) as saved:
        streak_corrupt = {
            key: np.array(saved[key], copy=True) for key in saved.files}
    streak_corrupt['gate_failures' if index != 1 else 'nve_failures'] = replacement
    streak_restart = os.path.join(root, f'streak-{index}.namd.restart.npz')
    np.savez_compressed(streak_restart, **streak_corrupt)
    d_streak = NAMD.__new__(NAMD); d_streak.mol = Mol(); d_streak.nstate = 2
    d_streak.dt_fs = 0.5; d_streak.seed = 1; d_streak.rng_stream = 2
    d_streak.restart_requested = True; d_streak.restart_file = streak_restart
    d_streak.trajectory_file = d.trajectory_file
    d_streak.nacme_audit_file = d.nacme_audit_file
    try:
        d_streak._load_restart()
    except RuntimeError:
        streak_rejections.append(True)
    else:
        streak_rejections.append(False)
gate_streak_metadata_rejected = all(streak_rejections)

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

# Charge and reference/excited-state spin settings define the electronic
# Hamiltonian and therefore participate in checkpoint compatibility.
d_charge = NAMD.__new__(NAMD); d_charge.mol = Mol(); d_charge.nstate = 2
d_charge.dt_fs = 0.5; d_charge.seed = 1; d_charge.rng_stream = 2
charge_mol = Mol()
charge_mol.config = {
    section: dict(settings) for section, settings in Mol.config.items()}
charge_mol.config['input']['charge'] = 1
d_wrong_charge = NAMD.__new__(NAMD); d_wrong_charge.mol = charge_mol
d_wrong_charge.nstate = 2; d_wrong_charge.dt_fs = 0.5
d_wrong_charge.seed = 1; d_wrong_charge.rng_stream = 2
charge_bound = d_charge._restart_signature() != d_wrong_charge._restart_signature()
spin_mol = Mol()
spin_mol.config = {
    section: dict(settings) for section, settings in Mol.config.items()}
spin_mol.config['tdhf']['multiplicity'] = 1
d_wrong_spin = NAMD.__new__(NAMD); d_wrong_spin.mol = spin_mol
d_wrong_spin.nstate = 2; d_wrong_spin.dt_fs = 0.5
d_wrong_spin.seed = 1; d_wrong_spin.rng_stream = 2
spin_bound = d_charge._restart_signature() != d_wrong_spin._restart_signature()
basis_mol = Mol()
basis_mol.config = {
    section: dict(settings) for section, settings in Mol.config.items()}
basis_mol.config['input']['library'] = 'H 6-31g; O aug-cc-pvdz'
basis_mol.config['input']['ispher'] = 0
d_wrong_basis = NAMD.__new__(NAMD); d_wrong_basis.mol = basis_mol
d_wrong_basis.nstate = 2; d_wrong_basis.dt_fs = 0.5
d_wrong_basis.seed = 1; d_wrong_basis.rng_stream = 2
basis_definition_bound = (
    d_charge._restart_signature() != d_wrong_basis._restart_signature())
pcm_mol = Mol(); pcm_mol.config = {
    section: dict(settings) for section, settings in Mol.config.items()}
pcm_mol.config['pcm'] = {'enabled': True, 'epsilon': 40.0, 'model': 'ddpcm'}
d_pcm = NAMD.__new__(NAMD); d_pcm.mol = pcm_mol; d_pcm.nstate = 2
d_pcm.dt_fs = 0.5; d_pcm.seed = 1; d_pcm.rng_stream = 2
pcm_bound = d_charge._restart_signature() != d_pcm._restart_signature()
operator_mol = Mol(); operator_mol.config = {
    section: dict(settings) for section, settings in Mol.config.items()}
operator_mol.config['tdhf']['cam_alpha'] = 0.25
d_operator = NAMD.__new__(NAMD); d_operator.mol = operator_mol
d_operator.nstate = 2; d_operator.dt_fs = 0.5
d_operator.seed = 1; d_operator.rng_stream = 2
tdhf_operator_bound = (
    d_charge._restart_signature() != d_operator._restart_signature())
grid_mol = Mol(); grid_mol.config = {
    section: dict(settings) for section, settings in Mol.config.items()}
grid_mol.config['dftgrid'] = {'rad_npts': 120, 'ang_npts': 434,
                              'hfscale': 0.55}
d_grid = NAMD.__new__(NAMD); d_grid.mol = grid_mol; d_grid.nstate = 2
d_grid.dt_fs = 0.5; d_grid.seed = 1; d_grid.rng_stream = 2
dftgrid_bound = d_charge._restart_signature() != d_grid._restart_signature()
policy_mol = Mol(); policy_mol.config = {
    section: dict(settings) for section, settings in Mol.config.items()}
policy_mol.config['md'] = {'nacme_gate': 'error',
                           'nacme_gate_consecutive': 5}
d_policy = NAMD.__new__(NAMD); d_policy.mol = policy_mol; d_policy.nstate = 2
d_policy.dt_fs = 0.5; d_policy.seed = 1; d_policy.rng_stream = 2
gate_policy_bound = d_charge._restart_signature() != d_policy._restart_signature()
scf_mol = Mol(); scf_mol.config = {
    section: dict(settings) for section, settings in Mol.config.items()}
scf_mol.config['scf']['scal_rel'] = 'x2c'
d_scf = NAMD.__new__(NAMD); d_scf.mol = scf_mol; d_scf.nstate = 2
d_scf.dt_fs = 0.5; d_scf.seed = 1; d_scf.rng_stream = 2
scf_settings_bound = d_charge._restart_signature() != d_scf._restart_signature()

# Rank-zero reconciliation failures must be broadcast before any rank raises;
# no unmatched explicit barrier is permitted on either simulated rank.
class FakeMPI:
    status = None
    barrier_calls = 0
    def __init__(self, rank):
        self.rank = rank
        self.use_mpi = 1
    def bcast(self, value, root=0):
        if self.rank == root:
            FakeMPI.status = value
        return FakeMPI.status
    def barrier(self):
        FakeMPI.barrier_calls += 1

broken_trajectory = os.path.join(root, 'broken.namd.trj')
with open(broken_trajectory, 'wb') as stream:
    stream.write(b'not-a-valid-header')
mpi_errors = []
for rank in (0, 1):
    d_mpi = NAMD.__new__(NAMD)
    d_mpi.mol = SimpleNamespace(mpi_manager=FakeMPI(rank))
    d_mpi.trajectory_file = broken_trajectory
    try:
        d_mpi._reconcile_trajectory_with_restart(1)
    except ValueError as error:
        mpi_errors.append(str(error))
collective_error_propagated = (
    len(mpi_errors) == 2 and mpi_errors[0] == mpi_errors[1]
    and FakeMPI.barrier_calls == 0)

# Checkpoint writes use the same collective failure propagation path.
FakeMPI.status = None
save_errors = []
for rank in (0, 1):
    d_save_mpi = NAMD.__new__(NAMD)
    d_save_mpi.mol = SimpleNamespace(mpi_manager=FakeMPI(rank))
    d_save_mpi.restart_interval = 1
    d_save_mpi._save_restart_on_io_rank = (
        lambda *_args: (_ for _ in ()).throw(ValueError('checkpoint failure')))
    try:
        d_save_mpi._save_restart(1, None, None, None)
    except ValueError as error:
        save_errors.append(str(error))
collective_save_error = len(save_errors) == 2 and save_errors[0] == save_errors[1]

# Fresh starts invalidate old runnable checkpoints, while path aliases are
# rejected before any sidecar is opened.
d_fresh = NAMD.__new__(NAMD); d_fresh.mol = SimpleNamespace(
    log=os.path.join(root, 'fresh.log'))
d_fresh.restart_requested = False
d_fresh.trajectory_file = os.path.join(root, 'fresh.trj')
d_fresh.nacme_audit_file = os.path.join(root, 'fresh.tsv')
d_fresh.restart_file = os.path.join(root, 'fresh.npz')
d_fresh.restart_manifest_file = os.path.join(root, 'fresh.oqp')
for path in (d_fresh.trajectory_file, d_fresh.nacme_audit_file,
             d_fresh.restart_file, d_fresh.restart_manifest_file):
    with open(path, 'w', encoding='utf-8') as stream:
        stream.write('stale')
d_fresh._prepare_md_outputs()
fresh_outputs_invalidated = (
    os.path.getsize(d_fresh.trajectory_file) == 0
    and os.path.getsize(d_fresh.nacme_audit_file) == 0
    and not os.path.exists(d_fresh.restart_file)
    and not os.path.exists(d_fresh.restart_manifest_file))
d_collision = NAMD.__new__(NAMD)
d_collision.mol = SimpleNamespace(log=os.path.join(root, 'collision.log'))
d_collision.trajectory_file = d_fresh.trajectory_file
d_collision.nacme_audit_file = d_fresh.trajectory_file
d_collision.restart_file = d_fresh.restart_file
d_collision.restart_manifest_file = d_fresh.restart_manifest_file
try:
    d_collision._validate_sidecar_paths()
except ValueError:
    sidecar_collision_rejected = True
else:
    sidecar_collision_rejected = False
d_log_collision = NAMD.__new__(NAMD)
d_log_collision.mol = SimpleNamespace(log=d_fresh.trajectory_file)
d_log_collision.trajectory_file = d_fresh.trajectory_file
d_log_collision.nacme_audit_file = d_fresh.nacme_audit_file
d_log_collision.restart_file = d_fresh.restart_file
d_log_collision.restart_manifest_file = d_fresh.restart_manifest_file
try:
    d_log_collision._validate_sidecar_paths()
except ValueError:
    log_collision_rejected = True
else:
    log_collision_rejected = False

def input_collision(candidate, *, velocity='zero', qmmm=None, config=None):
    probe = NAMD.__new__(NAMD)
    probe.mol = SimpleNamespace(
        log=os.path.join(root, 'input-collision.log'),
        oqp_input_source=os.path.join(input_root, 'request.oqp'),
        config=(config or {
            'input': {'method': 'tdhf'}, 'qmmm': qmmm or {}}))
    probe.velocity_source = velocity
    probe.trajectory_file = candidate
    probe.nacme_audit_file = os.path.join(root, 'input-collision.tsv')
    probe.restart_file = os.path.join(root, 'input-collision.npz')
    probe.restart_manifest_file = os.path.join(root, 'input-collision.oqp')
    try:
        probe._validate_sidecar_paths()
    except ValueError:
        return True
    return False

input_collisions_rejected = all((
    input_collision(os.path.join(input_root, 'request.oqp')),
    input_collision(os.path.join(input_root, 'vel.dat'), velocity='vel.dat'),
    input_collision(
        os.path.join(input_root, 'water.pdb'),
        qmmm={'pdb_file': 'water.pdb'}),
    input_collision(
        os.path.join(input_root, 'local.xml'),
        qmmm={'forcefield_files': 'local.xml tip3p.xml'}),
    input_collision(
        os.path.join(input_root, 'local.xml'), config={
            'input': {'method': 'dftb'},
            'dftb': {'backend': 'native', 'type': 'mrsf', 'model': 'mio',
                     'parameter_path': 'local.xml',
                     'library_path': 'local.xml'}}),
    input_collision(
        os.path.join(input_root, 'parameters', 'H-H.skf'), config={
            'input': {'method': 'dftb'},
            'dftb': {'backend': 'native', 'type': 'mrsf', 'model': 'mio',
                     'parameter_path': 'parameters',
                     'library_path': 'local.xml'}}),
))

# A discontinuity invalidates the previous gate result before the current
# dense record is written.
d_ba = NAMD.__new__(NAMD); d_ba.nacme_check = 'baeck_an'; d_ba.nstate = 2
d_ba.dt = 1.0; d_ba._ba_energy_left = np.array([0.0, 0.1])
d_ba._ba_energy_center = np.array([0.0, 0.2]); d_ba._ba_tdc_left = np.eye(2)
d_ba._ba_dt_left = 1.0; d_ba._ba_last = {'stale': True}
d_ba._nacme_gate_failures = 2
d_ba._nacme_gate_last = {'stale': True}; d_ba._nacme_reference_tdc = np.eye(2)
d_ba._nacme_reference_mask = np.ones((2, 2)); d_ba._nacme_reference_source = 1
d_ba._compute_tdc = lambda _overlap: np.zeros((2, 2))
d_ba.mol = SimpleNamespace(data={
    'OQP::td_energies_old': np.array([0.0, 0.3]),
    'OQP::td_energies': np.array([0.0, 0.4]),
})
d_ba._update_baeck_an_check(3, np.eye(2))
stale_gate_cleared = (
    d_ba._ba_last is None and d_ba._nacme_gate_last is None
    and d_ba._nacme_reference_tdc is None
    and d_ba._nacme_reference_mask is None
    and d_ba._nacme_reference_source == 0
    and d_ba._nacme_gate_failures == 0)
bad_live = NAMD.__new__(NAMD); bad_live.nacme_check = 'baeck_an'
bad_live.nstate = 2; bad_live.mol = SimpleNamespace(data={
    'OQP::td_energies_old': np.array([0.0]),
    'OQP::td_energies': np.array([0.0, 0.2]),
})
try:
    bad_live._update_baeck_an_check(2, np.eye(2))
except RuntimeError:
    short_live_energy_rejected = True
else:
    short_live_energy_rejected = False
bad_hop = NAMD.__new__(NAMD); bad_hop.nstate = 2
bad_hop.mol = SimpleNamespace(data={'OQP::td_energies': np.array([0.0])})
try:
    bad_hop._validated_td_energies('OQP::td_energies')
except RuntimeError:
    hop_energy_rejected = True
else:
    hop_energy_rejected = False

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
    'phase': records['tracking_phase'].tolist(),
    'phase_initial': records['tracking_phase_initial'].tolist(),
    'gate_streak': records['gate_streak'].tolist(),
    'lineage': records['tracking_lineage'].tolist(),
    'raw_order': records['tracking_raw_order'].tolist(),
    'nstate': header['nstate'],
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
        'broadcast_load_rejected': broadcast_load_rejected,
        'history_load_rejected': history_load_rejected,
        'gate_streak_metadata_rejected': gate_streak_metadata_rejected,
        'molecule_mismatch_rejected': molecule_mismatch_rejected,
        'qm_selection_bound': qm_selection_bound,
        'charge_bound': charge_bound, 'spin_bound': spin_bound,
        'basis_definition_bound': basis_definition_bound,
        'pcm_bound': pcm_bound, 'tdhf_operator_bound': tdhf_operator_bound,
        'dftgrid_bound': dftgrid_bound, 'gate_policy_bound': gate_policy_bound,
        'collective_error_propagated': collective_error_propagated,
        'collective_save_error': collective_save_error,
        'fresh_outputs_invalidated': fresh_outputs_invalidated,
        'sidecar_collision_rejected': sidecar_collision_rejected,
        'log_collision_rejected': log_collision_rejected,
        'input_collisions_rejected': input_collisions_rejected,
        'stale_gate_cleared': stale_gate_cleared,
        'short_live_energy_rejected': short_live_energy_rejected,
        'hop_energy_rejected': hop_energy_rejected,
        'tight_binding_bound': tight_binding_bound,
        'default_tb_artifact_bound': default_tb_artifact_bound,
        'default_tb_model_stable': default_tb_model_stable,
        'scf_settings_bound': scf_settings_bound,
        'unique_manifests': unique_manifests,
        'forcefield_identity_stable': forcefield_identity_stable,
        'builtin_forcefield_fingerprinted': builtin_forcefield_fingerprinted,
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
        'phase_initial': [[-1.0, 1.0]], 'lineage': [[7, 4]],
        'gate_streak': [0],
        'raw_order': [[1, 0]],
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
        'broadcast_load_rejected': True,
        'history_load_rejected': True,
        'gate_streak_metadata_rejected': True,
        'qm_selection_bound': True, 'charge_bound': True, 'spin_bound': True,
        'basis_definition_bound': True,
        'pcm_bound': True, 'tdhf_operator_bound': True,
        'dftgrid_bound': True, 'gate_policy_bound': True,
        'collective_error_propagated': True, 'unique_manifests': True,
        'collective_save_error': True, 'fresh_outputs_invalidated': True,
        'sidecar_collision_rejected': True, 'stale_gate_cleared': True,
        'log_collision_rejected': True, 'short_live_energy_rejected': True,
        'input_collisions_rejected': True,
        'hop_energy_rejected': True,
        'tight_binding_bound': True,
        'default_tb_artifact_bound': True, 'scf_settings_bound': True,
        'default_tb_model_stable': True,
        'forcefield_identity_stable': True,
        'builtin_forcefield_fingerprinted': True,
        'coefficient_shape_rejected': True, 'tracking_shape_rejected': True,
        'audit_rows': [
            'center_step\tsource\tverdict\tsigned\tcompared_pairs\t'
            'invariant_failures\treference_failures\t'
            'consecutive_reference_failures\tcandidate_diagonal_max\t'
            'candidate_antisymmetry_max\treference_diagonal_max\t'
            'reference_antisymmetry_max\tpair_rms_error\tpair_max_error\t'
            'max_tolerance_ratio',
            '0\tTD-Baeck-An\tpass\tFalse\t1\t0\t0\t0\t\t\t\t\t\t\t',
        ],
    }
