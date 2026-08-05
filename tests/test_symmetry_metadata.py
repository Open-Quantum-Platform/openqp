import importlib.util
import json
import sys
import types
import tempfile
from pathlib import Path
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]



def load_molecule_module():
    if 'oqp.molecule.molecule' in sys.modules:
        return sys.modules['oqp.molecule.molecule']

    # Minimal runtime stubs so we can import the molecule module without a compiled library.
    package_root = ROOT / 'pyoqp'

    oqp_mod = types.ModuleType('oqp')
    oqp_mod.__path__ = [str(package_root)]

    utils_mod = types.ModuleType('oqp.utils')
    utils_mod.__path__ = [str(package_root / 'oqp' / 'utils')]

    mol_pkg = types.ModuleType('oqp.molecule')
    mol_pkg.__path__ = [str(package_root / 'oqp' / 'molecule')]

    mpi_mod = types.ModuleType('oqp.utils.mpi_utils')

    input_parser_mod = types.ModuleType('oqp.utils.input_parser')

    molden_mod = types.ModuleType('oqp.molden')
    molden_writer_mod = types.ModuleType('oqp.molden.moldenwriter')

    class DummyFFI:
        def buffer(self, *_args, **_kwargs):
            return b""

        def sizeof(self, *_args, **_kwargs):
            return 8

    class DummyLib:
        def oqp_init(self):
            return types.SimpleNamespace()

        def oqp_clean(self, *_args, **_kwargs):
            return None

    class OQPConfigParser:
        def __init__(self, *_, **__):
            self._data = {}

        def read(self, *_args, **_kwargs):
            return None

        def load_dict(self, data):
            self._data = data

        def print_config(self):
            return None

        def validate(self):
            return self._data

    class MPIManager:
        size = 1
        use_mpi = False

    def mpi_get_attr(func):
        return func

    def mpi_dump(func):
        return func

    class MoldenWriter:
        def __init__(self, *_args, **_kwargs):
            pass

    class OQPData:
        def __init__(self, *_, **__):
            self._store = {}
            self._data = types.SimpleNamespace(grad=np.array([]))

        def __getitem__(self, key):
            return self._store[key]

        def __setitem__(self, key, value):
            self._store[key] = np.array(value)

        def apply_config(self, config):
            self.config = config

        @property
        def natom(self):
            return self._store.get('natom', 0)

    class OQP_DATA:
        pass

    setattr(input_parser_mod, 'OQPConfigParser', OQPConfigParser)
    setattr(mpi_mod, 'MPIManager', MPIManager)
    setattr(mpi_mod, 'mpi_get_attr', mpi_get_attr)
    setattr(mpi_mod, 'mpi_dump', mpi_dump)
    setattr(molden_writer_mod, 'MoldenWriter', MoldenWriter)

    setattr(oqp_mod, 'ffi', DummyFFI())
    setattr(oqp_mod, 'lib', DummyLib())
    setattr(oqp_mod, 'utils', utils_mod)
    setattr(utils_mod, 'mpi_utils', mpi_mod)
    setattr(oqp_mod, 'molden', molden_mod)
    setattr(molden_mod, 'moldenwriter', molden_writer_mod)

    oqpdata_stub = types.ModuleType('oqp.molecule.oqpdata')
    setattr(oqpdata_stub, 'OQPData', OQPData)
    setattr(oqpdata_stub, 'OQP_DATA', OQP_DATA)
    setattr(
        oqpdata_stub,
        'OQP_CONFIG_SCHEMA',
        {
            'symmetry': {
                'enabled': {'type': str, 'default': 'false'},
                'point_group': {'type': str, 'default': 'auto'},
                'subgroup': {'type': str, 'default': 'auto'},
                'label_mo': {'type': bool, 'default': 'True'},
                'label_states': {'type': bool, 'default': 'True'},
                'label_modes': {'type': bool, 'default': 'True'},
                'use_integral_symmetry': {'type': bool, 'default': 'False'},
                'use_response_symmetry': {'type': bool, 'default': 'False'},
                'tolerance': {'type': float, 'default': '1.0e-5'},
                'strict': {'type': bool, 'default': 'False'},
            }
        },
    )

    sys.modules.update(
        {
            'oqp': oqp_mod,
            'oqp.utils': utils_mod,
            'oqp.utils.mpi_utils': mpi_mod,
            'oqp.utils.input_parser': input_parser_mod,
            'oqp.molden': molden_mod,
            'oqp.molden.moldenwriter': molden_writer_mod,
            'oqp.molecule': mol_pkg,
            'oqp.molecule.oqpdata': oqpdata_stub,
        }
    )

    spec = importlib.util.spec_from_file_location(
        'oqp.molecule.molecule',
        ROOT / 'pyoqp/oqp/molecule/molecule.py',
        submodule_search_locations=[str(ROOT / 'pyoqp/oqp/molecule')],
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    module.__package__ = 'oqp.molecule'
    sys.modules['oqp.molecule.molecule'] = module
    spec.loader.exec_module(module)
    return module


class TestSymmetryMetadata(unittest.TestCase):
    def test_symmetry_metadata_defaults_to_c1(self):
        molecule_module = load_molecule_module()
        molecule = molecule_module.Molecule.__new__(molecule_module.Molecule)
        molecule.config = {}
        molecule.symmetry_metadata = {}

        metadata = molecule.initialize_symmetry_metadata()

        self.assertEqual(metadata['point_group'], 'c1')
        self.assertEqual(metadata['subgroup'], 'c1')
        self.assertEqual(metadata['status'], 'disabled')
        self.assertFalse(metadata['use_integral_symmetry'])
        self.assertFalse(metadata['use_response_symmetry'])
        self.assertEqual(metadata['label_mo'], True)
        self.assertEqual(metadata['label_states'], True)
        self.assertEqual(metadata['label_modes'], True)

    def test_symmetry_metadata_can_request_auto_c2v(self):
        molecule_module = load_molecule_module()
        molecule = molecule_module.Molecule.__new__(molecule_module.Molecule)
        molecule.config = {
            'symmetry': {
                'enabled': 'auto',
                'point_group': 'c2v',
                'subgroup': 'c1',
                'label_mo': False,
                'label_states': False,
                'label_modes': False,
                'use_integral_symmetry': 'False',
                'use_response_symmetry': 'False',
                'tolerance': 2e-4,
                'strict': True,
            }
        }
        molecule.symmetry_metadata = {}

        metadata = molecule.initialize_symmetry_metadata()

        self.assertEqual(metadata['status'], 'auto')
        self.assertEqual(metadata['requested_point_group'], 'c2v')
        self.assertEqual(metadata['detected_point_group'], 'c2v')
        self.assertEqual(metadata['requested_subgroup'], 'c1')
        self.assertEqual(metadata['detected_subgroup'], 'c1')
        self.assertFalse(metadata['use_integral_symmetry'])
        self.assertFalse(metadata['use_response_symmetry'])
        self.assertFalse(metadata['label_mo'])
        self.assertFalse(metadata['label_states'])
        self.assertFalse(metadata['label_modes'])
        self.assertTrue(metadata['strict'])
        self.assertAlmostEqual(metadata['tolerance'], 2e-4)

    def _molecule_for_put_data(self, coord):
        molecule_module = load_molecule_module()
        molecule = molecule_module.Molecule.__new__(molecule_module.Molecule)
        molecule.tag = []
        molecule.config = {}
        molecule.config_tag = {}
        molecule._state_tracking_fresh = True
        molecule.get_system = lambda: coord
        molecule.symmetry_metadata = {
            'status': 'disabled',
            'point_group': 'c1',
            'subgroup': 'c1',
            'requested_point_group': 'auto',
            'requested_subgroup': 'auto',
            'label_mo': True,
            'label_states': True,
            'label_modes': True,
            'use_integral_symmetry': False,
            'use_response_symmetry': False,
            'strict': False,
            'tolerance': 1e-5,
            'raw': {'enabled': 'false'},
        }
        return molecule

    def test_symmetry_metadata_round_trips_at_the_same_geometry(self):
        coord = [0.0, 0.0, 0.0, 0.0, 0.0, 1.4]
        molecule = self._molecule_for_put_data(coord)

        molecule.put_data({
            'coord': coord,
            'symmetry_metadata': {'point_group': 'c2v',
                                  'detection': {'point_group': 'c2v'}},
        })

        self.assertEqual(molecule.symmetry_metadata['point_group'], 'c2v')
        self.assertEqual(
            molecule.symmetry_metadata['detection']['point_group'], 'c2v')

    def test_geometry_derived_metadata_does_not_follow_a_guess_to_a_new_geometry(self):
        """A guess file is read routinely at a *different* geometry.

        Every optimiser step and every finite-difference displacement reads the
        reference job's guess, and a displaced geometry generally has lower
        symmetry than the reference -- so carrying the producer's point group,
        operations and standard-frame transform over would assert symmetry the
        molecule does not have.
        """
        molecule = self._molecule_for_put_data([0.0, 0.0, 0.0, 0.0, 0.0, 1.4])

        molecule.put_data({
            'coord': [0.0, 0.0, 0.0, 0.0, 0.0, 1.405],
            'symmetry_metadata': {'point_group': 'c2v',
                                  'detection': {'point_group': 'c2v'},
                                  'integral_symmetry': {'status': 'reoriented'}},
        })

        # Falls back to this job's own detection rather than vanishing --
        # consumers index these keys directly.
        self.assertEqual(molecule.symmetry_metadata['point_group'], 'c1')
        self.assertNotIn('detection', molecule.symmetry_metadata)
        self.assertNotIn('integral_symmetry', molecule.symmetry_metadata)

    def test_staging_state_is_never_restored_from_a_guess_file(self):
        """These have a Fortran-side half that put_data does not restore.

        Inheriting the Python half alone claims a reduction that is not staged
        in this process: symmetrize_gradient gates on status == 'active' and
        would project with the producer's operations, and _petite_is_staged
        would let the stability path switch the Fortran flag on with no maps.
        A job that never asked for integral symmetry could pick both up.
        """
        coord = [0.0, 0.0, 0.0, 0.0, 0.0, 1.4]
        molecule = self._molecule_for_put_data(coord)

        molecule.put_data({
            'coord': coord,  # same geometry -- the permissive case
            'symmetry_metadata': {
                'integral_symmetry': {'status': 'active'},
                'reduction_maps': {'n_operations': 8},
                'reduction_maps_full': {'n_operations': 24},
                'detection': {'point_group': 'c2v'},
            },
        })

        for key in ('integral_symmetry', 'reduction_maps',
                    'reduction_maps_full'):
            self.assertNotIn(key, molecule.symmetry_metadata, key)
        # Descriptive geometry data still comes across when it matches.
        self.assertEqual(
            molecule.symmetry_metadata['detection']['point_group'], 'c2v')

    def test_guess_file_never_overrides_this_job_symmetry_switches(self):
        """The bug this pins: a numerical-Hessian child reoriented itself even
        though its own input said use_integral_symmetry=false, because the
        parent's guess JSON carried True and won. Water/6-31G HF frequencies
        went from 1828.86/3906.50/4001.54 to -2675.60/697.64/2728.85 cm^-1.
        """
        coord = [0.0, 0.0, 0.0, 0.0, 0.0, 1.4]
        molecule = self._molecule_for_put_data(coord)

        molecule.put_data({
            'coord': coord,
            'symmetry_metadata': {'use_integral_symmetry': True,
                                  'use_response_symmetry': True,
                                  'status': 'enabled'},
        })

        self.assertFalse(molecule.symmetry_metadata['use_integral_symmetry'])
        self.assertFalse(molecule.symmetry_metadata['use_response_symmetry'])
        self.assertEqual(molecule.symmetry_metadata['status'], 'disabled')

    def test_get_results_and_hess_json_include_symmetry_metadata(self):
        molecule_module = load_molecule_module()
        molecule = molecule_module.Molecule.__new__(molecule_module.Molecule)
        molecule.symmetry_metadata = {'status': 'disabled', 'point_group': 'c1', 'subgroup': 'c1',
                                    'requested_point_group': 'auto', 'requested_subgroup': 'auto',
                                    'label_mo': True, 'label_states': True, 'label_modes': True,
                                    'use_integral_symmetry': False, 'use_response_symmetry': False,
                                    'strict': False, 'tolerance': 1e-5}
        molecule.mol_energy = types.SimpleNamespace(energy=-1.23)
        molecule.log = '/tmp/placeholder.log'
        molecule.idx = 1
        molecule.config = {}
        molecule.mrsf_ekt_results_by_kind = {}

        class _StubData:
            def __getitem__(self, key):
                return np.array([])

        molecule.data = _StubData()
        molecule.energies = np.array([-1.23])
        molecule.hessian = np.eye(9)
        molecule.hessian_metadata = {}
        molecule.freqs = np.array([1.0])
        molecule.modes = np.ones((1, 9))
        molecule.inertia = np.ones(3)
        molecule.infrared_intensities = np.zeros(0)
        molecule.raman_activities = np.zeros(0)
        molecule.vibrational_intensity_metadata = {}
        molecule.infrared_mode_dipole_derivatives = np.zeros((0, 3))
        molecule.raman_mode_polarizability_derivatives = np.zeros((0, 3, 3))

        molecule.get_atoms = lambda: np.array([1, 1, 8], dtype=int)
        molecule.get_system = lambda: np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        molecule.get_mass = lambda: np.array([1.0, 1.0, 16.0])
        molecule.get_grad = lambda: []
        molecule.get_nac = lambda: []
        molecule.get_soc = lambda: []
        molecule.get_hess = lambda: []

        results = molecule.get_results()

        self.assertIn('symmetry_metadata', results)
        self.assertEqual(results['symmetry_metadata']['point_group'], 'c1')

        molecule.get_data = lambda: {}
        molecule.get_grad = lambda: np.array([])
        molecule.get_nac = lambda: []
        molecule.get_soc = lambda: []
        molecule.get_hess = lambda: []
        with tempfile.TemporaryDirectory() as tmp:
            molecule.log = str(Path(tmp) / 'run.log')
            molecule.save_freqs(0)
            with open(str(Path(tmp) / 'run.hess.json'), 'r', encoding='utf-8') as f:
                data = json.load(f)

            molecule.hessian = np.zeros((9, 9))
            molecule.freqs = np.zeros(1)
            # Thermochemistry may read a Hessian sidecar before any SinglePoint
            # calculation has populated Molecule.energies.
            molecule.energies = None
            molecule.read_freqs()
            self.assertTrue(np.array_equal(molecule.hessian, np.eye(9)))
            self.assertTrue(np.array_equal(molecule.freqs, np.array([1.0])))

            molecule.get_mass = lambda: np.array([2.0, 1.0, 16.0])
            with self.assertRaisesRegex(ValueError, 'isotopic masses'):
                molecule.read_freqs()
            molecule.get_mass = lambda: np.array([1.0, 1.0, 16.0])

            molecule.config = {'input': {'basis': 'different'}}
            with self.assertRaisesRegex(ValueError, 'electronic-model configuration/state'):
                molecule.read_freqs()

            molecule.config = {}
            corrupt = dict(data)
            corrupt['freqs'] = [float('nan')]
            with open(str(Path(tmp) / 'run.hess.json'), 'w', encoding='utf-8') as f:
                json.dump(corrupt, f)
            with self.assertRaisesRegex(ValueError, 'cached frequencies'):
                molecule.read_freqs()

            unsigned = dict(data)
            unsigned.pop('hessian_cache_version')
            unsigned.pop('hessian_request')
            with open(str(Path(tmp) / 'run.hess.json'), 'w', encoding='utf-8') as f:
                json.dump(unsigned, f)
            with self.assertRaisesRegex(ValueError, 'current versioned model-configuration/state'):
                molecule.read_freqs()

        self.assertIn('symmetry_metadata', data)
        self.assertEqual(data['symmetry_metadata']['subgroup'], 'c1')
        self.assertIn('hessian_request', data)
        self.assertEqual(data['hessian_cache_version'], 2)


if __name__ == '__main__':
    unittest.main()
