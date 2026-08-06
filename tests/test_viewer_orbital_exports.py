from io import StringIO

import numpy as np

from oqp.molden.moldenwriter import MoldenWriter, write_frequency
from oqp.molecule.molecule import Molecule


class FakeData:
    def __init__(self):
        self.basis = {
            'centers': np.array([0, 0]),
            'angs': np.array([0, 0]),
            'ncontr': np.array([1, 1]),
            'alpha': np.array([1.0, 0.5]),
            'coef': np.array([1.0, 0.8]),
            'nsh': 2,
            'nbf': 2,
        }
        self.values = {
            'OQP::VEC_MO_A': np.array([1.0, 2.0, 3.0, 4.0]),
            'OQP::E_MO_A': np.array([-0.5, 0.2]),
            'nocc': 1,
        }

    def get_basis(self):
        return self.basis

    def __getitem__(self, key):
        if key in self.values:
            return self.values[key]
        raise AttributeError(key)


def fake_molecule():
    molecule = Molecule.__new__(Molecule)
    molecule.data = FakeData()
    molecule.config = {
        'scf': {'type': 'rhf'},
        'tdhf': {'target': 2},
        'ekt': {'ip': True, 'ea': False},
    }
    molecule.mrsf_ekt_results_by_kind = {
        'ip': {
            'eigenvalues_hartree': [-0.4],
            'ebe_ev': [10.8845544983952],
            'pole_strengths': [0.75],
            'dyson_orbitals_mo': [[0.5, 0.25]],
        }
    }
    molecule.freqs = np.array([1234.5])
    molecule.modes = np.array([[0.1, 0.2, 0.3]])
    molecule.infrared_intensities = np.array([12.5])
    molecule.raman_activities = np.array([34.75])
    return molecule


def test_portable_json_exports_scf_and_state_specific_dyson_ao_coefficients():
    molecule = fake_molecule()
    portable = molecule._viewer_orbital_data()
    assert portable['basis_set']['basis_function_order'] == 'molden'
    assert portable['molecular_orbitals']['alpha']['coefficients'] == [[1.0, 2.0], [3.0, 4.0]]

    states = molecule._viewer_dyson_data()['dyson_orbitals']['states']
    assert len(states) == 1
    assert states[0]['label'] == 'Dyson IP state 1'
    assert states[0]['parent_state'] == 2
    assert states[0]['pole_strength'] == 0.75
    assert np.allclose(states[0]['coefficients'], [1.25, 2.0])


def test_molden_writer_combines_mo_dyson_labels_and_frequency_sections():
    molecule = fake_molecule()
    output = StringIO()
    writer = MoldenWriter(output)
    basis = molecule.data.get_basis()
    writer.write_atoms(1, np.array([1]), np.array([[0.0, 0.0, 0.0]]))
    writer.write_basis(1, basis)
    writer.write_mo(
        basis,
        np.array([[1.25, 2.0]]),
        [-0.4],
        [0.75],
        spin='Alpha',
        symmetries=['Dyson-IP-state-1'],
        already_reordered=True,
    )
    molecule.get_atoms = lambda: np.array([1])
    molecule.get_system = lambda: np.array([0.0, 0.0, 0.0])
    writer.write_frequency(molecule, molecule.freqs, molecule.modes)
    text = output.getvalue()
    assert text.count('[Molden Format]') == 1
    assert '\n1 0\n' in text
    assert '\ns 1 1.00\n' in text
    assert 'Sym= Dyson-IP-state-1' in text
    assert '[GTO]' in text and '[MO]' in text
    assert '[FREQ]' in text and '[FR-NORM-COORD]' in text
    assert '[INT]\n     12.50000000\n' in text
    assert '[RAMAN]\n     34.75000000\n' in text
    assert '12.50000000      34.75000000' not in text


def test_basis_writer_does_not_mutate_openqp_primitive_coefficients():
    basis = fake_molecule().data.get_basis()
    original = basis['coef'].copy()
    writer = MoldenWriter(StringIO())
    writer.write_basis(1, basis)
    assert np.array_equal(basis['coef'], original)


def test_standard_mo_records_always_include_a_symmetry_field():
    basis = fake_molecule().data.get_basis()
    output = StringIO()
    writer = MoldenWriter(output)
    writer.write_mo(basis, np.eye(2), [-0.5, 0.2], [2.0, 0.0], spin='Alpha')
    assert output.getvalue().count('Sym= A') == 2


def test_frequency_only_writer_keeps_one_canonical_molden_header():
    molecule = fake_molecule()
    molecule.get_atoms = lambda: np.array([1])
    molecule.get_system = lambda: np.array([0.0, 0.0, 0.0])
    text = write_frequency(molecule, molecule.freqs, molecule.modes)
    assert text.startswith('[Molden Format]\n')
    assert text.count('[Molden Format]') == 1


def test_portable_json_omits_unsupported_h_shells_without_aborting_save():
    molecule = fake_molecule()
    nbf = 21
    molecule.data.basis = {
        'centers': np.array([0]),
        'angs': np.array([5]),
        'ncontr': np.array([1]),
        'alpha': np.array([1.0]),
        'coef': np.array([1.0]),
        'nsh': 1,
        'nbf': nbf,
    }
    molecule.data.values['OQP::VEC_MO_A'] = np.eye(nbf).reshape(-1)
    molecule.data.values['OQP::E_MO_A'] = np.arange(nbf, dtype=float)
    molecule.data.values['nocc'] = 1

    assert molecule.has_molden_orbitals() is False
    assert molecule._viewer_basis_data() == {}
    assert molecule._viewer_orbital_data() == {}
    assert molecule._viewer_dyson_data() == {}


def test_portable_json_omits_mixed_cartesian_spherical_basis_without_bad_indices():
    molecule = fake_molecule()
    nbf = 11  # one Cartesian D shell (6) plus one pure D shell (5)
    molecule.data.basis = {
        'centers': np.array([0, 0]),
        'angs': np.array([2, 2]),
        'ncontr': np.array([1, 1]),
        'alpha': np.array([1.0, 0.5]),
        'coef': np.array([1.0, 0.8]),
        'nsh': 2,
        'nbf': nbf,
    }
    molecule.data.values['OQP::VEC_MO_A'] = np.eye(nbf).reshape(-1)
    molecule.data.values['OQP::E_MO_A'] = np.arange(nbf, dtype=float)
    molecule.data.values['nocc'] = 1

    assert MoldenWriter._harmonic_representation(molecule.data.basis) is None
    assert molecule.has_molden_orbitals() is False
    assert molecule._viewer_orbital_data() == {}


def test_legacy_direct_ekt_uses_tdhf_type_for_dyson_kind():
    molecule = fake_molecule()
    records = molecule.mrsf_ekt_results_by_kind['ip']
    molecule.mrsf_ekt_results_by_kind = {}
    molecule._read_mrsf_ekt_records = lambda: records
    molecule.config['tdhf']['type'] = 'mrsf_ekt_ea'
    # Preserve the normal [ekt] defaults that previously misclassified EA as IP.
    molecule.config['ekt'] = {'ip': True, 'ea': False}

    states = molecule._viewer_dyson_data()['dyson_orbitals']['states']
    assert states[0]['kind'] == 'EA'
    assert states[0]['label'] == 'Dyson EA state 1'


def test_non_ekt_restart_does_not_export_stale_raw_dyson_records():
    molecule = fake_molecule()
    molecule.mrsf_ekt_results_by_kind = {}
    molecule.config['input'] = {'runtype': 'energy'}
    molecule.config['tdhf']['type'] = 'mrsf'

    assert molecule._viewer_dyson_data() == {}


def test_write_molden_requires_explicit_post_kernel_dyson_append(tmp_path):
    molecule = fake_molecule()
    molecule.usempi = False
    molecule.elem = np.array([1])
    molecule.xyz = np.zeros((1, 3))
    molecule.data.values['natom'] = 1

    scf_file = tmp_path / 'scf.molden'
    molecule.write_molden(scf_file)
    assert 'Dyson-IP-state-1' not in scf_file.read_text()

    dyson_file = tmp_path / 'dyson.molden'
    molecule.write_molden(dyson_file, include_dyson=True)
    assert 'Sym= Dyson-IP-state-1' in dyson_file.read_text()
