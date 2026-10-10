"""Smaller-dt continuation preserves a checkpoint without editing its source."""
import json
import struct
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from oqp.library import namd as mod


def driver(tmp_path, dt=0.1, name='source'):
    d = mod.NAMD.__new__(mod.NAMD)
    d.dt_fs = dt
    d.dt_adaptive = False
    d._time_origin_fs = 0.0
    d._continuation_provenance = None
    d.nstate = 2
    d.nstep = 9
    d.active = 2
    d.first_hop_step = 1
    d.seed, d.rng_stream, d._rng_step = 12, 3, 7
    d.restart_requested = False
    d.continuation_checkpoint = d.continuation_trajectory = ''
    d.trajectory_file = str(tmp_path / (name + '.trj'))
    d.restart_file = str(tmp_path / (name + '.npz'))
    d.restart_manifest_file = str(tmp_path / (name + '.restart.oqp'))
    d.zpredict_audit_file = str(tmp_path / (name + '.tsv'))
    d._restart_signature = lambda: json.dumps({'dt_fs': d.dt_fs, 'nstate': 2,
        'tdhf_settings': {'zvconv': 1e-7}, 'seed': d.seed, 'guess_settings': {}})
    d.mol = SimpleNamespace(log=str(tmp_path / (name + '.log')), put_data=lambda data: None)
    d._molecular_identity = lambda: {'atomic_numbers': [1]}
    d._independent_settings_record = lambda: {}
    d._validate_sidecar_paths = lambda: None
    d._write_restart_manifest = lambda: None
    d.coef = np.array([0.6+0.2j, np.sqrt(0.6)*1j])
    d.prev_xyz = np.array([1., 2., 3.])
    d.prev_data = {k: np.eye(2) for k in ('OQP::VEC_MO_A', 'OQP::VEC_MO_B',
        'OQP::E_MO_A', 'OQP::E_MO_B', 'OQP::DM_A', 'OQP::DM_B',
        'OQP::FOCK_A', 'OQP::FOCK_B', 'OQP::SM', 'OQP::td_bvec_mo')}
    d.prev_data.update({'OQP::td_energies': np.array([0., .1]),
        'OQP::state_tracking_phase_initial': np.array([-1., 1.]),
        'OQP::state_tracking_lineage': np.array([1, 0])})
    d._ba_energy_left = d._ba_energy_center = d._ba_tdc_left = d._ba_dt_left = None
    d._nve_reference_energy = d._nve_previous_energy = -0.8
    d._nacme_gate_failures = d._nve_gate_failures = 0
    d._trajectory_prefix_hasher = None
    d._trajectory_prefix_stat = None
    return d


def source_checkpoint(tmp_path, configure=None):
    d = driver(tmp_path)
    if configure is not None:
        configure(d)
    dtype = mod._namd_trajectory_dtype(2, 1, 0)
    record = np.zeros(1, dtype=dtype)
    record['step'], record['time_fs'], record['active'] = 7, .7, 2
    record['coordinates_bohr'] = d.prev_xyz.reshape(1, 3)
    record['velocities_au'] = np.array([[.01, .02, .03]])
    record['coef_real'], record['coef_imag'] = d.coef.real, d.coef.imag
    header = {'schema_version': mod.NAMD_TRAJECTORY_SCHEMA_VERSION,
              'nstate': 2, 'natom': 1, 'ncv': 0, 'record_bytes': dtype.itemsize,
              'signature': d._restart_signature()}
    encoded = json.dumps(header).encode()
    Path(d.trajectory_file).write_bytes(mod.NAMD_TRAJECTORY_MAGIC +
        struct.pack('<Q', len(encoded)) + encoded + record.tobytes())
    d._save_restart_on_io_rank(7, d.prev_xyz.reshape(1, 3),
        np.array([[.01, .02, .03]]), np.array([[.001, .002, .003]]))
    # An uncommitted failed point must remain unchanged in the source.
    record['step'], record['time_fs'] = 8, .8
    with open(d.trajectory_file, 'ab') as stream:
        stream.write(record.tobytes())
    return d


def child(tmp_path, source):
    d = driver(tmp_path, .05, 'child')
    d.continuation_checkpoint = source.restart_file
    d.continuation_trajectory = source.trajectory_file
    return d


def test_continuation_preserves_full_state_source_rng_and_time(tmp_path, monkeypatch):
    monkeypatch.delattr(mod.hashlib, "file_digest", raising=False)
    monkeypatch.setattr(mod, 'dump_log', lambda *args, **kwargs: None)
    source = source_checkpoint(tmp_path)
    original = {p: Path(p).read_bytes() for p in (source.restart_file, source.trajectory_file)}
    d = child(tmp_path, source)
    loaded = []
    d.mol.put_data = loaded.append
    state = d._load_restart()
    assert state['step'] == 7 and d._rng_step == 7 and d.active == 2
    assert np.array_equal(d.coef, source.coef)
    assert np.array_equal(state['acceleration'], [[.001, .002, .003]])
    for k, value in source.prev_data.items():
        assert np.array_equal(loaded[0][k], value)
    assert d._etot_prev == -.8 and d._nve_reference_energy == -.8
    assert d._physical_time_fs(7) == pytest.approx(.7)
    assert d._physical_time_fs(8) == pytest.approx(.75)
    assert d._physical_time_fs(9) == pytest.approx(.8)
    d._prepare_hop_step(8)
    assert d._rng_step == 8
    assert all(Path(p).read_bytes() == value for p, value in original.items())
    header, records = mod.read_namd_trajectory(d.trajectory_file)
    assert records['step'].tolist() == [7]
    assert records['time_fs'].tolist() == [.7]
    assert header['continuation']['source_dt_fs'] == .1
    del records
    # A normal restart of the child restores its nonzero physical time origin.
    resumed = driver(tmp_path, .05, 'child')
    resumed.restart_requested = True
    resumed._load_restart()
    assert resumed._physical_time_fs(9) == pytest.approx(.8)
    assert resumed._etot_prev == -.8
    assert np.array_equal(resumed.coef, d.coef)


def test_restart_restores_nvt_energy_recovery_baseline(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, 'dump_log', lambda *args, **kwargs: None)

    def nvt(d):
        # NVT keeps no NVE history; energy recovery uses its own baseline.
        d._nve_reference_energy = d._nve_previous_energy = None
        d._etot_prev = -.75
        d._disc_energy_absorbed = .002

    source_checkpoint(tmp_path, configure=nvt)
    resumed = driver(tmp_path)
    resumed.restart_requested = True
    resumed._nve_reference_energy = resumed._nve_previous_energy = None
    resumed._etot_prev = None
    resumed._disc_energy_absorbed = 0.0
    resumed._load_restart()
    assert resumed._etot_prev == pytest.approx(-.75)
    assert resumed._disc_energy_absorbed == pytest.approx(.002)


@pytest.mark.parametrize('change', ['larger_dt', 'equal_dt', 'nan_dt', 'seed', 'solver'])
def test_continuation_rejects_unrequested_identity_changes(tmp_path, change):
    d = driver(tmp_path, .05)
    saved = json.loads(d._restart_signature())
    saved['dt_fs'] = .1
    if change == 'larger_dt': d.dt_fs = .2
    if change == 'equal_dt': d.dt_fs = .1
    if change == 'nan_dt': d.dt_fs = float('nan')
    if change == 'seed': d.seed += 1
    if change == 'solver': saved['tdhf_settings']['zvconv'] = 1e-16
    assert not d._restart_signature_matches(json.dumps(saved), allow_smaller_dt=True)
    assert not d._restart_signature_matches(json.dumps(saved)) or change == 'equal_dt'


@pytest.mark.parametrize('alias', ['direct', 'symlink', 'hardlink', 'existing'])
def test_source_output_isolation(tmp_path, alias):
    source = source_checkpoint(tmp_path)
    d = child(tmp_path, source)
    if alias == 'direct': d.trajectory_file = source.trajectory_file
    elif alias == 'symlink': Path(d.trajectory_file).symlink_to(source.trajectory_file)
    elif alias == 'hardlink': Path(d.trajectory_file).hardlink_to(source.trajectory_file)
    else: Path(d.trajectory_file).write_text('prior output')
    before = Path(source.trajectory_file).read_bytes()
    with pytest.raises(ValueError): d._load_continuation_on_io_rank()
    assert Path(source.trajectory_file).read_bytes() == before


def test_tampered_prefix_is_rejected_before_output_creation(tmp_path):
    source = source_checkpoint(tmp_path)
    d = child(tmp_path, source)
    path = Path(source.trajectory_file)
    raw = bytearray(path.read_bytes());raw[100] ^= 1;path.write_bytes(raw)
    with pytest.raises((ValueError, KeyError)):
        d._load_continuation_on_io_rank()
    assert not Path(d.trajectory_file).exists()


def test_physical_time_metadata_must_match_saved_step(tmp_path):
    source = source_checkpoint(tmp_path)
    with np.load(source.restart_file, allow_pickle=False) as f:
        data = {k: f[k] for k in f.files}
    data['time_fs'] = np.array([99.])
    np.savez_compressed(source.restart_file, **data)
    with pytest.raises(ValueError, match='physical time'):
        child(tmp_path, source)._load_continuation_on_io_rank()


@pytest.mark.parametrize('alias', ['direct', 'symlink', 'hardlink'])
def test_runner_rejects_log_alias_before_logging(tmp_path, alias):
    import ast
    import os
    path = Path(__file__).parents[1] / 'pyoqp/oqp/pyoqp.py'
    tree = ast.parse(path.read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                    and n.name == '_protect_continuation_log')
    namespace = {'os': os}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(path), 'exec'), namespace)
    source = tmp_path / 'source.npz';source.write_bytes(b'checkpoint')
    output = tmp_path / 'child.log'
    if alias == 'direct': output = source
    elif alias == 'symlink': output.symlink_to(source)
    else: output.hardlink_to(source)
    with pytest.raises(ValueError, match='calculation log'):
        namespace['_protect_continuation_log']({'md': {'continuation_checkpoint': str(source)}}, str(output))
    assert source.read_bytes() == b'checkpoint'
    assert path.read_text().index('_protect_continuation_log(self.mol.config') < path.read_text().index("dump_log(self.mol, title='', section='start'")


def test_child_restart_manifest_removes_continuation_options(tmp_path):
    from oqp.utils.oqp_input import parse_canonical_oqp
    d = driver(tmp_path, .05, 'child')
    example = Path(__file__).parents[1] / 'examples/namd_local_continuation/source.continuation.oqp'
    d.mol.oqp_canonical_input = example.read_text()
    d.mol.oqp_input_source = str(example)
    d._rebase_restart_spec_paths = lambda spec, source_dir: spec
    mod.NAMD._write_restart_manifest(d)
    spec = parse_canonical_oqp(Path(d.restart_manifest_file).read_text())
    assert spec.driver.kwargs['restart'] is True
    assert 'continuation_checkpoint' not in spec.driver.kwargs
    assert 'continuation_trajectory' not in spec.driver.kwargs
    assert spec.driver.kwargs['nstep'] == 3


@pytest.mark.parametrize('restored_dt', [.075, .1])
def test_return_toward_original_dt_preserves_state_and_time(tmp_path, monkeypatch, restored_dt):
    monkeypatch.setattr(mod, 'dump_log', lambda *args, **kwargs: None)
    source = source_checkpoint(tmp_path)
    fine = child(tmp_path, source)
    fine._load_restart()
    returning = driver(tmp_path, restored_dt, 'return')
    returning.continuation_checkpoint = fine.restart_file
    returning.continuation_trajectory = fine.trajectory_file
    before = Path(fine.restart_file).read_bytes()
    returning._load_restart()
    assert returning._physical_time_fs(7) == pytest.approx(.7)
    assert returning._physical_time_fs(8) == pytest.approx(.7 + restored_dt)
    assert returning._rng_step == 7
    assert np.array_equal(returning.coef, fine.coef)
    assert Path(fine.restart_file).read_bytes() == before
    # Ordinary restart must still reject a changed dt.
    returning.dt_fs = .06
    with pytest.raises(ValueError, match='mismatch'):
        returning._load_restart_on_io_rank()


def test_return_cannot_exceed_original_dt(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, 'dump_log', lambda *args, **kwargs: None)
    source = source_checkpoint(tmp_path)
    fine = child(tmp_path, source)
    fine._load_restart()
    returning = driver(tmp_path, .2, 'return')
    returning.continuation_checkpoint = fine.restart_file
    returning.continuation_trajectory = fine.trajectory_file
    with pytest.raises(ValueError, match='mismatch'):
        returning._load_restart()
    assert not Path(returning.trajectory_file).exists()


@pytest.mark.parametrize("content", [b"", b"abc", b"x" * (1024 * 1024 + 13)])
def test_checkpoint_hash_matches_sha256_without_python311_api(content, monkeypatch):
    import io
    import hashlib
    monkeypatch.delattr(hashlib, "file_digest", raising=False)
    assert mod._sha256_stream(io.BytesIO(content)) == hashlib.sha256(content).hexdigest()
