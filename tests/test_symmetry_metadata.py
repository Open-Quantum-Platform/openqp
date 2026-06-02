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

    def test_symmetry_metadata_round_trips_through_data(self):
        molecule_module = load_molecule_module()
        molecule = molecule_module.Molecule.__new__(molecule_module.Molecule)
        molecule.tag = []
        molecule.config = {}
        molecule.config_tag = {}
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

        payload = {'symmetry_metadata': {'point_group': 'c2v', 'raw': {'point_group': 'c2v'}}}
        molecule.put_data(payload)

        self.assertEqual(molecule.symmetry_metadata['point_group'], 'c2v')

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

        class _StubData:
            def __getitem__(self, key):
                return np.array([])

        molecule.data = _StubData()
        molecule.energies = np.array([-1.23])
        molecule.hessian = np.array([1.0])
        molecule.freqs = np.array([1.0])
        molecule.modes = np.array([1.0])

        molecule.get_atoms = lambda: np.array([1, 1, 8], dtype=int)
        molecule.get_system = lambda: np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        molecule.get_mass = lambda: np.array([1.0])
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

        self.assertIn('symmetry_metadata', data)
        self.assertEqual(data['symmetry_metadata']['subgroup'], 'c1')


if __name__ == '__main__':
    unittest.main()
