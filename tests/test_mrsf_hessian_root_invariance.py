"""Isolated-root and raw-invariance gates for numerical MRSF Hessians."""

from __future__ import annotations

import importlib.util
import json
import sys
import types
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load_source_module(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _single_point_without_native_runtime():
    """Load single_point.py with the same backend-free boundary as CI."""

    names = [
        'oqp', 'oqp.molecule', 'oqp.utils', 'oqp.utils.mpi_utils',
        'oqp.library', 'oqp.library.frequency', 'oqp.library.state_tracking',
        'oqp.library.nac_utils', 'oqp.utils.tb_backends',
        'oqp.utils.file_utils', 'oqp.utils.state_labels', 'oqp.utils.qmmm',
    ]
    snapshot = {name: sys.modules.get(name) for name in names}
    oqp = types.ModuleType('oqp')
    oqp.lib = types.SimpleNamespace()
    oqp.ffi = types.SimpleNamespace()
    sys.modules['oqp'] = oqp
    molecule = types.ModuleType('oqp.molecule')
    molecule.Molecule = object
    sys.modules['oqp.molecule'] = molecule
    utils = types.ModuleType('oqp.utils')
    sys.modules['oqp.utils'] = utils
    mpi = types.ModuleType('oqp.utils.mpi_utils')
    mpi.MPIManager = type('MPIManager', (), {'use_mpi': False, 'rank': 0})
    mpi.MPIPool = object
    sys.modules['oqp.utils.mpi_utils'] = mpi
    library = types.ModuleType('oqp.library')
    sys.modules['oqp.library'] = library
    frequency = types.ModuleType('oqp.library.frequency')
    frequency.normal_mode = lambda *_args, **_kwargs: ([], [], [])
    frequency.thermal_analysis = lambda *_args, **_kwargs: {}
    sys.modules['oqp.library.frequency'] = frequency
    native_tracking = types.ModuleType('oqp.library.state_tracking')
    native_tracking.diagonal_phase_tracking = lambda *_args, **_kwargs: None
    native_tracking.maximum_overlap_assignment = lambda *_args, **_kwargs: None
    sys.modules['oqp.library.state_tracking'] = native_tracking
    nac_utils = types.ModuleType('oqp.library.nac_utils')
    for name in (
            'canonical_state_overlap', 'hst_derivative_coupling',
            'interstate_coupling', 'load_numerical_nac_cache',
            'write_numerical_nac_cache_marker'):
        setattr(nac_utils, name, lambda *_args, **_kwargs: None)
    sys.modules['oqp.library.nac_utils'] = nac_utils
    tb = types.ModuleType('oqp.utils.tb_backends')
    tb.is_tb_method = lambda *_args, **_kwargs: False
    tb.make_tb_adapter = lambda *_args, **_kwargs: None
    tb.tb_config = lambda *_args, **_kwargs: {}
    sys.modules['oqp.utils.tb_backends'] = tb
    files = types.ModuleType('oqp.utils.file_utils')
    for name in ('dump_log', 'dump_data', 'write_config', 'write_xyz'):
        setattr(files, name, lambda *_args, **_kwargs: None)
    sys.modules['oqp.utils.file_utils'] = files
    labels = types.ModuleType('oqp.utils.state_labels')
    labels.is_mrsf = lambda *_args, **_kwargs: True
    labels.public_state_label = lambda *_args, **_kwargs: 'S1'
    sys.modules['oqp.utils.state_labels'] = labels
    sys.modules['oqp.utils.qmmm'] = types.ModuleType('oqp.utils.qmmm')
    module = _load_source_module(
        'single_point_mrsf_hessian_invariance',
        'pyoqp/oqp/library/single_point.py')
    return module, snapshot


def _restore_modules(snapshot):
    sys.modules.pop('single_point_mrsf_hessian_invariance', None)
    for name, module in snapshot.items():
        if module is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = module


def _molecule_module_without_native_runtime():
    # Other collected tests may already have cached a differently stubbed
    # Molecule module. Force this source checkout to load, then restore every
    # pre-existing oqp module so the fixture is order-independent.
    snapshot = {
        name: module for name, module in sys.modules.items()
        if name == 'oqp' or name.startswith('oqp.')
    }
    sys.modules.pop('oqp.molecule.molecule', None)
    helper = _load_source_module(
        'mrsf_hessian_symmetry_metadata_loader',
        'tests/test_symmetry_metadata.py')
    # The helper stubs the oqp package, so pin the matching source-tree JSON
    # adapter explicitly instead of inheriting an installed OpenQP copy.
    _load_source_module(
        'oqp.utils.json_utils',
        'pyoqp/oqp/utils/json_utils.py')
    module = helper.load_molecule_module()
    for name in list(sys.modules):
        if (name == 'oqp' or name.startswith('oqp.')) and name not in snapshot:
            sys.modules.pop(name, None)
    sys.modules.update(snapshot)
    return module


def _assignment(overlap):
    """Small permutation-only oracle for the tracker unit fixtures."""

    matrix = np.asarray(overlap, dtype=float)
    order = np.argmax(np.abs(matrix), axis=1).astype(np.int32)
    assert len(set(order.tolist())) == matrix.shape[0]
    chosen = matrix[np.arange(matrix.shape[0]), order]
    signs = np.where(chosen < 0.0, -1.0, 1.0)
    matched = np.abs(chosen)
    margins = np.empty(matrix.shape[0])
    for row, column in enumerate(order):
        other = np.delete(np.abs(matrix[row]), column)
        margins[row] = matched[row] - (np.max(other) if other.size else 0.0)
    return order, signs, matched, margins


class DummyMRSFMolecule:
    def __init__(self):
        # Current raw order is old [1, 2, 0].  The public S1 root is therefore
        # raw displaced root 3 and must be transported back to slot 1.
        central_rows = np.zeros((3, 6))
        central_rows[0, 0] = 1.0
        central_rows[1, 1] = 1.0
        central_rows[2, 2] = 1.0
        current_rows = central_rows[[1, 2, 0]]
        # Native tag shape is (packed, state), but the contiguous buffer is
        # state-major.  Keep this deliberately non-square so a layout mistake
        # cannot be hidden by an identity matrix.
        self.central_bridge = central_rows.reshape((6, 3))
        self.data = {
            'OQP::td_bvec_mo_old': self.central_bridge.copy(),
            'OQP::td_bvec_mo': current_rows.reshape((6, 3)),
            'OQP::td_energies_old': np.array([0.10, 0.20, 0.30]),
            'OQP::td_energies': np.array([0.20, 0.30, 0.10]),
            'nelec_A': np.array([2]),
            'nelec_B': np.array([0]),
        }
        self.config = {
            'tdhf': {'type': 'mrsf', 'multiplicity': 1},
            'nac': {'dt': 1.0},
        }
        self.energies = np.array([-10.0, -9.80, -9.70, -9.90])
        self._previous_symmetry_metadata = {
            'status': 'active',
            'detected_subgroup': 'c2v',
            'state_labels': {
                'status': 'ok',
                'labels': ['a1', 'b1', 'b2'],
                'terms': ['1A1', '1B1', '1B2'],
            },
        }
        self.symmetry_metadata = {
            'status': 'active',
            'detected_subgroup': 'c2v',
            'state_labels': {
                'status': 'ok',
                'labels': ['b1', 'b2', 'a1'],
                'terms': ['1B1', '1B2', '1A1'],
            },
        }
        self._state_tracking_fresh = False


@pytest.fixture
def single_point():
    module, snapshot = _single_point_without_native_runtime()
    try:
        yield module
    finally:
        _restore_modules(snapshot)


@pytest.fixture
def tracker(single_point, monkeypatch):
    monkeypatch.setattr(single_point, 'maximum_overlap_assignment', _assignment)
    monkeypatch.setattr(single_point, 'dump_log', lambda *_args, **_kwargs: None)
    obj = single_point.NACME.__new__(single_point.NACME)
    obj.mol = DummyMRSFMolecule()
    return obj


def test_isolated_root_exchange_transports_amplitude_energy_and_irrep(tracker):
    report = tracker.track_isolated_mrsf_hessian_root(1)

    assert report['selected_raw_state'] == 3
    assert report['transported_state'] == 1
    assert report['matched_overlap'] == pytest.approx(1.0)
    assert report['symmetry_check'] == 'same_subgroup_same_irrep'
    assert report['spin_check'] == 'fixed_spin_adapted_mrsf_manifold'
    assert report['representation'] == 'spin_adapted_mrsf_co_ov_cv_oo'
    assert report['response_overlap_metric'] == \
        'spin_adapted_mrsf_physical_packed'
    assert report['expanded_response_dimension'] == 6
    assert report['physical_response_dimension'] == 5
    assert report['degenerate_manifold_solver'] is False
    np.testing.assert_allclose(
        tracker.mol.data['OQP::td_bvec_mo'], tracker.mol.central_bridge)
    np.testing.assert_allclose(
        tracker.mol.data['OQP::td_energies'], [0.10, 0.20, 0.30])
    np.testing.assert_allclose(
        tracker.mol.energies, [-10.0, -9.90, -9.80, -9.70])
    assert tracker.mol.symmetry_metadata['state_labels']['labels'] == [
        'a1', 'b1', 'b2']


def test_subgroup_lowering_keeps_pure_reference_overlap_gate(tracker):
    tracker.mol.symmetry_metadata['detected_subgroup'] = 'c1'
    tracker.mol.symmetry_metadata['state_labels']['labels'] = ['a', 'a', 'a']

    report = tracker.track_isolated_mrsf_hessian_root(1)

    assert report['reference_irrep'] == 'a1'
    assert report['selected_raw_irrep'] == 'a'
    assert report['symmetry_check'] == 'pure_reference_subgroup_lowered'


def test_same_subgroup_irrep_change_is_rejected(tracker):
    tracker.mol.symmetry_metadata['state_labels']['labels'][2] = 'b2'

    with pytest.raises(RuntimeError, match='changed spatial irrep'):
        tracker.track_isolated_mrsf_hessian_root(1)


def test_missing_requested_symmetry_metadata_is_rejected(tracker):
    tracker.mol._previous_symmetry_metadata = None

    with pytest.raises(RuntimeError, match='central symmetry metadata'):
        tracker.track_isolated_mrsf_hessian_root(1)


def test_weak_overlap_is_rejected_before_gradient(tracker):
    angle = np.arccos(0.98)
    current_rows = tracker.mol.data['OQP::td_bvec_mo'].reshape((3, 6)).copy()
    current_rows[2, :] = 0.0
    current_rows[2, 0] = np.cos(angle)
    current_rows[2, 4] = np.sin(angle)
    tracker.mol.data['OQP::td_bvec_mo'] = current_rows.reshape((6, 3))

    with pytest.raises(RuntimeError, match='root overlap'):
        tracker.track_isolated_mrsf_hessian_root(1)


def test_h2_spin_adapted_metric_normalizes_the_folded_oo_response(tracker):
    current_rows = tracker.mol.data['OQP::td_bvec_mo'].reshape((3, 6)).copy()
    current_rows[2, :] *= 7.0
    tracker.mol.data['OQP::td_bvec_mo'] = current_rows.reshape((6, 3))

    report = tracker.track_isolated_mrsf_hessian_root(1)

    assert report['matched_overlap'] == pytest.approx(1.0)


def test_h2_nearby_spin_adapted_response_has_unit_limit(tracker):
    angle = 1.0e-3
    current_rows = tracker.mol.data['OQP::td_bvec_mo'].reshape((3, 6)).copy()
    current_rows[2, :] = 0.0
    current_rows[2, 0] = np.cos(angle)
    current_rows[2, 4] = np.sin(angle)
    tracker.mol.data['OQP::td_bvec_mo'] = current_rows.reshape((6, 3))

    report = tracker.track_isolated_mrsf_hessian_root(1)

    assert report['matched_overlap'] == pytest.approx(np.cos(angle))
    assert report['matched_overlap'] > 0.999999


def test_h2_spin_adapted_metric_rejects_redundant_r_contamination(tracker):
    current_rows = tracker.mol.data['OQP::td_bvec_mo'].reshape((3, 6)).copy()
    current_rows[2, 3] = 1.0e-4
    tracker.mol.data['OQP::td_bvec_mo'] = current_rows.reshape((6, 3))

    with pytest.raises(RuntimeError, match='nonphysical displaced OO'):
        tracker.track_isolated_mrsf_hessian_root(1)


def test_h2_triplet_metric_retains_only_physical_oo_l(single_point):
    previous = np.zeros((1, 6))
    current = np.zeros((1, 6))
    previous[0, 0] = 1.0
    current[0, 0] = -3.0

    overlap, physical_dimension = (
        single_point._mrsf_spin_adapted_response_overlap(
            current, previous, nalpha=2, nbeta=0, multiplicity=3)
    )

    assert physical_dimension == 3
    np.testing.assert_allclose(overlap, [[-1.0]])

    current[0, 1] = 1.0e-4
    with pytest.raises(RuntimeError, match='nonphysical displaced OO'):
        single_point._mrsf_spin_adapted_response_overlap(
            current, previous, nalpha=2, nbeta=0, multiplicity=3)


def test_near_degenerate_root_requires_projector_treatment(tracker):
    tracker.mol.data['OQP::td_energies_old'][1] = 0.100001

    with pytest.raises(RuntimeError, match='near-degenerate'):
        tracker.track_isolated_mrsf_hessian_root(1)


def test_highest_solved_root_is_not_claimed_isolated(tracker):
    with pytest.raises(RuntimeError, match='not bracketed'):
        tracker.track_isolated_mrsf_hessian_root(3)


def test_unvalidated_spin_extension_is_rejected(tracker):
    tracker.mol.config['tdhf']['multiplicity'] = 5

    with pytest.raises(RuntimeError, match='singlet and triplet'):
        tracker.track_isolated_mrsf_hessian_root(1)


def test_tracking_sidecar_is_bound_to_state_and_displacement(single_point, tmp_path):
    path = tmp_path / 'disp.tracking.json'
    report = {
        'status': 'accepted',
        'reference_state': 1,
        'displacement_tag': 'c0p',
    }
    single_point._write_mrsf_hessian_tracking(path, report)

    assert single_point._read_mrsf_hessian_tracking(
        path, state=1, tag='c0p') is not None
    assert single_point._read_mrsf_hessian_tracking(
        path, state=2, tag='c0p') is None
    assert single_point._read_mrsf_hessian_tracking(
        path, state=1, tag='c0m') is None


def test_numerical_hessian_worker_tracks_before_gradient_and_binds_restart():
    runfunc = (ROOT / 'pyoqp/oqp/library/runfunc.py').read_text(encoding='utf-8')
    single = (ROOT / 'pyoqp/oqp/library/single_point.py').read_text(encoding='utf-8')

    tracking_call = 'track_isolated_mrsf_hessian_root('
    assert runfunc.index(tracking_call) < runfunc.index(
        'Gradient(mol).gradient()', runfunc.index(tracking_call))
    assert "config['guess']['file2'] = guess_file" in single
    assert "config['nac']['align'] = 'reorder'" in single
    assert 'missing or invalid MRSF Hessian root-tracking report' in single
    assert "'OQP_NUM_HESS_WORKER': '1'" in single


def test_basis_overlap_keeps_displaced_symmetry_and_captures_central_metadata(
        single_point, tmp_path):
    central_meta = {
        'status': 'active', 'detected_subgroup': 'c2v',
        'state_labels': {'status': 'ok', 'labels': ['a1', 'b1']},
    }
    displaced_meta = {
        'status': 'active', 'detected_subgroup': 'c1',
        'state_labels': {'status': 'ok', 'labels': ['a', 'a']},
    }
    previous_file = tmp_path / 'central.json'
    previous_file.write_text(
        json.dumps({'symmetry_metadata': central_meta}),
        encoding='utf-8')

    class Data(dict):
        mol2 = np.array([])

    class Mol:
        def __init__(self):
            self.system = np.array([0.0, 0.0, 0.01])
            self.symmetry_metadata = dict(displaced_meta)
            self.config = {
                'guess': {
                    'file': str(tmp_path / 'displaced.json'),
                    'file2': str(previous_file),
                    'continue_geom': False,
                },
                'input': {'method': 'tdhf'},
            }
            self.data = Data({'natom': 1})
            self._arrays = {
                'OQP::VEC_MO_A': np.eye(2), 'OQP::VEC_MO_B': np.eye(2),
                'OQP::E_MO_A': np.array([-1.0, 0.2]),
                'OQP::E_MO_B': np.array([-1.0, 0.2]),
                'OQP::td_bvec_mo': np.eye(2),
                'OQP::td_energies': np.array([0.1, 0.2]),
            }

        def get_system(self):
            return self.system.copy()

        def get_data(self):
            return {key: np.array(value).copy() for key, value in self._arrays.items()}

        def load_data(self):
            self.system = np.zeros(3)
            self.symmetry_metadata = dict(central_meta)

        def update_system(self, value):
            self.system = np.asarray(value, dtype=float).copy()

        def put_data(self, _data):
            return None

    tracker = single_point.BasisOverlap.__new__(single_point.BasisOverlap)
    tracker.mol = Mol()
    tracker.back_door = False
    tracker.load_previous_data()

    assert tracker.mol.symmetry_metadata['detected_subgroup'] == 'c1'
    assert tracker.mol._previous_symmetry_metadata[
        'detected_subgroup'] == 'c2v'


def _diagnostic_molecule():
    Molecule = _molecule_module_without_native_runtime().Molecule
    mol = Molecule.__new__(Molecule)
    mol.data = {'natom': 2}
    mol.get_system = lambda: np.array([-0.5, 0.0, 0.0, 0.5, 0.0, 0.0])
    mol.get_mass = lambda: np.array([1.0, 1.0])
    return mol


def test_raw_hessian_rigid_motion_residuals_precede_symmetrization():
    mol = _diagnostic_molecule()
    Molecule = type(mol)
    hessian = np.zeros((6, 6))
    hessian[0, 0] = hessian[3, 3] = 2.0
    hessian[0, 3] = hessian[3, 0] = -2.0

    final = Molecule.set_hessian_result(
        mol, hessian, reference_gradient=np.zeros(6))
    raw = mol.hessian_metadata['pre_symmetrization_invariance']

    np.testing.assert_allclose(final, hessian)
    assert raw['stage'] == 'raw_cartesian_before_symmetrization_or_rigid_motion_projection'
    assert raw['max_abs_asymmetry'] == 0.0
    assert raw['translation_max_abs_residual'] == 0.0
    assert raw['rotation_max_abs_residual'] == 0.0
    assert raw['rotation_gradient_covariance_correction'] is True
    assert raw['rigid_motion_projection_applied'] is False


def test_raw_asymmetry_is_recorded_before_the_final_average():
    mol = _diagnostic_molecule()
    Molecule = type(mol)
    raw = np.zeros((6, 6))
    raw[0, 1] = 0.2

    final = Molecule.set_hessian_result(mol, raw, asymmetry_tol=np.inf)

    assert mol.hessian_metadata['max_asymmetry'] == pytest.approx(0.2)
    assert mol.hessian_metadata['pre_symmetrization_invariance'][
        'max_abs_asymmetry'] == pytest.approx(0.2)
    assert final[0, 1] == pytest.approx(0.1)
    assert final[1, 0] == pytest.approx(0.1)


def test_native_mrsf_processed_hessian_metadata_is_not_labeled_raw():
    mol = _diagnostic_molecule()
    Molecule = type(mol)
    hessian = np.zeros((6, 6))
    hessian[0, 0] = hessian[3, 3] = 2.0
    hessian[0, 3] = hessian[3, 0] = -2.0

    final = Molecule.set_hessian_result(
        mol,
        hessian,
        producer_stage=(
            'native_mrsf_after_response_row_symmetrization_and_'
            'translation_projection'),
        upstream_symmetrization_applied=True,
        translation_projection_applied=True,
        rotation_projection_applied=False,
    )

    np.testing.assert_allclose(final, hessian)
    metadata = mol.hessian_metadata
    assert metadata['upstream_symmetrization_applied'] is True
    assert metadata['python_symmetrization_applied'] is False
    assert metadata['translation_projection_applied'] is True
    assert metadata['rotation_projection_applied'] is False
    assert 'pre_symmetrization_invariance' not in metadata
    stored = metadata['stored_matrix_invariance']
    assert stored['stage'] == (
        'native_mrsf_after_response_row_symmetrization_and_'
        'translation_projection')
    assert stored['upstream_symmetrization_applied'] is True
    assert stored['translation_projection_applied'] is True
    assert stored['rotation_projection_applied'] is False
