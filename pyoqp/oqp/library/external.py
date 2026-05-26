"""OQP external quantum chemical program"""
import os
import sys
import json
from pathlib import Path

import numpy as np
import oqp
from oqp.utils.constants import ANGSTROM_TO_BOHR
from oqp.periodic_table import SYMBOL_MAP, ELEMENTS_NAME

try:
    from pyscf import gto, dft, lib

except ModuleNotFoundError:
    gto = None
    dft = None
    lib = None

pyscf_functional = {
    "pbe0": "pbe0",
    "b3lyp": "b3lyp",
    "bhhlyp": "bhandhlyp",
    "cam-b3lyp": "camb3lyp",
    "camb3lyp": "camb3lyp",
}


def _json_array(value):
    """Convert NumPy/scalar values into JSON-serializable Python values."""

    return np.asarray(value).tolist()


def _spin_pair(value):
    """Return alpha/beta arrays from PySCF RHF/ROHF/UHF-style data."""

    if isinstance(value, (tuple, list)) and len(value) == 2:
        return np.asarray(value[0]), np.asarray(value[1])
    array = np.asarray(value)
    if array.ndim >= 3 and array.shape[0] == 2:
        return array[0], array[1]
    return array, array


def _pack_lower_triangle(matrix):
    """Pack a square matrix in OpenQP's lower-triangular restart layout."""

    mat = np.asarray(matrix)
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise ValueError(f"Expected a square density matrix, got shape={mat.shape}")
    idx = np.tril_indices(mat.shape[0])
    return mat[idx].tolist()


def _density_pair_from_mf(mf):
    """Return alpha/beta AO density matrices from a PySCF mean-field object."""

    mo_coeff = mf.mo_coeff
    mo_occ = getattr(mf, 'mo_occ', None)
    if isinstance(mo_coeff, (tuple, list)) and len(mo_coeff) == 2:
        return _spin_pair(mf.make_rdm1())
    if isinstance(mo_occ, (tuple, list)) and len(mo_occ) == 2:
        return _spin_pair(mf.make_rdm1())

    coeff = np.asarray(mo_coeff)
    occ = np.asarray(mo_occ) if mo_occ is not None else None
    if coeff.ndim == 2 and occ is not None and occ.ndim == 1:
        occ_a = np.where(occ > 0.0, 1.0, 0.0)
        occ_b = np.where(occ > 1.0, 1.0, 0.0)
        dm_a = np.dot(coeff * occ_a, coeff.T)
        dm_b = np.dot(coeff * occ_b, coeff.T)
        return dm_a, dm_b

    return _spin_pair(mf.make_rdm1())


def _pyscf_total_energy(mf):
    """Read the total energy from a PySCF mean-field object if available."""

    try:
        return float(mf.energy_tot())
    except Exception:
        pass
    try:
        return float(mf.e_tot)
    except Exception:
        pass
    try:
        return float(mf.mol.energy_tot())
    except Exception:
        return 0.0


def export_pyscf_guess_to_openqp_json(mf, mol, filename='pyscf_wfn.json'):
    """Export a PySCF mean-field object as an OpenQP JSON restart guess.

    This native exporter replaces the previous MOKIT ``py2openqp`` bridge for
    PySCF-backed guess modes. It writes the OpenQP fields consumed by
    ``Molecule.load_data``: MO coefficients, orbital energies, packed density
    matrices, atoms, coordinates, and basis/library metadata.
    """

    filename = Path(filename)
    coeff_a, coeff_b = _spin_pair(mf.mo_coeff)
    energy_a, energy_b = _spin_pair(mf.mo_energy)
    dm_a, dm_b = _density_pair_from_mf(mf)

    data = {
        'OQP::DM_A': _pack_lower_triangle(dm_a),
        'OQP::DM_B': _pack_lower_triangle(dm_b),
        'OQP::E_MO_A': _json_array(energy_a),
        'OQP::E_MO_B': _json_array(energy_b),
        'OQP::VEC_MO_A': _json_array(coeff_a),
        'OQP::VEC_MO_B': _json_array(coeff_b),
        'atoms': _json_array(mol.get_atoms()),
        'coord': _json_array(mol.get_system()),
        'energy': _pyscf_total_energy(mf),
        'td_energies': [0.0],
        'grad': [],
        'nac': [],
        'soc': [],
        'hess': [],
        'json': {
            'scf_type': mol.config['scf']['type'],
            'basis': mol.config['input']['basis'],
            'library': mol.config['input']['library'],
        },
    }

    with open(filename, 'w') as outdata:
        json.dump(data, outdata, indent=2)
    return filename


def load_pyscf_guess(mf, mol, filename='pyscf_wfn.json'):
    """Export a PySCF guess, load it into OpenQP, and project the basis."""

    target_basis = mol.config['input']['basis']
    target_library = mol.config['input']['library']
    guess_file = export_pyscf_guess_to_openqp_json(mf, mol, filename)
    mol.config['guess']['file'] = str(guess_file)
    mol.load_data()

    mol.data.set_scf_active_basis(1)
    mol.config['input']['basis'] = mol.config['json']['basis']
    mol.config['input']['library'] = mol.config['json']['library']
    oqp.library.set_basis(mol)

    mol.data.set_scf_active_basis(0)
    mol.config['input']['basis'] = target_basis
    mol.config['input']['library'] = target_library
    oqp.library.project_basis(mol)


def guess_from_pyscf(mol):
    """Set up pySCF calculation"""

    if gto is None:
        print(f'\n   PyOQP cannot find PySCF\n')
        ## pyscf run mpi with a different lib, which might cause conflict with mpi4py, disabling mpi now for future dev
        sys.exit(1)

    if mol.usempi:
        print(f'\n   PyOQP cannot run pySCF with MPI\n')
        ## pyscf run mpi with a different lib, which might cause conflict with mpi4py, disabling mpi now for future dev
        sys.exit(1)

    atoms = []
    coord = mol.get_system().reshape((-1, 3)) * ANGSTROM_TO_BOHR
    for n, at in enumerate(mol.get_atoms()):
        atoms.append([ELEMENTS_NAME[SYMBOL_MAP[int(at)]], coord[n]])

    mole = gto.Mole()
    mole.atom = atoms
    mole.unit = 'Bohr'
    mole.basis = mol.config["input"]["basis"]
    mole.charge = mol.config["input"]["charge"]
    mole.spin = mol.config["scf"]["multiplicity"] - 1
    mole.output = '%s.pyscf' % mol.project_name
    mole.verbose = 4
    mole.build(cart=True)
    functional = mol.config["input"]["functional"]

    if mol.config["scf"]["type"] == 'rohf':
        mf = dft.ROKS(mole)
    elif mol.config["scf"]["type"] == 'uhf':
        mf = dft.UKS(mole)
    elif mol.config["scf"]["type"] == 'rhf':
        mf = dft.RKS(mole)
    else:
        print(f'\n   PyOQP only support pySCF with rhf, uhf, and rohf\n')
        sys.exit(1)

    try:
        xcfun = pyscf_functional[functional]
    except KeyError:
        xcfun = 'pbe0'

    mf.xc = xcfun
    mf.max_cycle = mol.config["scf"]["maxit"]
    mf.level_shift = mol.config["scf"]["vshift"]

    try:
        os.environ["PYSCF_MAX_MEMORY"]
    except KeyError:
        mf.max_memory = 1000  # MB
    try:
        os.environ["OMP_NUM_THREADS"]
    except KeyError:
        lib.num_threads(1)

    mf.kernel()

    if mf.converged is False:
        mf = mf.newton()
        mf.kernel()

    load_pyscf_guess(mf, mol)

    return

def guess_from_pyscf_initial_density(mol, guess_type):
    """Build a PySCF atomic-density/potential initial guess for OpenQP.

    The resulting PySCF orbitals are exported through OpenQP's native JSON
    restart writer. ``sad`` uses PySCF's atomic-density guess. ``sap`` requests
    PySCF's superposition-of-atomic-potentials guess when available.
    """

    if gto is None:
        print(f'\n   PyOQP cannot find PySCF\n')
        sys.exit(1)

    if mol.usempi:
        print(f'\n   PyOQP cannot run pySCF with MPI\n')
        sys.exit(1)

    guess_type = guess_type.lower()
    if guess_type not in {'sad', 'sap'}:
        raise ValueError(f'Unsupported PySCF initial-density guess={guess_type}')

    atoms = []
    coord = mol.get_system().reshape((-1, 3)) * ANGSTROM_TO_BOHR
    for n, at in enumerate(mol.get_atoms()):
        atoms.append([ELEMENTS_NAME[SYMBOL_MAP[int(at)]], coord[n]])

    mole = gto.Mole()
    mole.atom = atoms
    mole.unit = 'Bohr'
    mole.basis = mol.config["input"]["basis"]
    mole.charge = mol.config["input"]["charge"]
    mole.spin = mol.config["scf"]["multiplicity"] - 1
    mole.output = '%s.%s.pyscf' % (mol.project_name, guess_type)
    mole.verbose = 4
    mole.build(cart=True)

    if mol.config["scf"]["type"] == 'rohf':
        mf = dft.ROKS(mole)
    elif mol.config["scf"]["type"] == 'uhf':
        mf = dft.UKS(mole)
    elif mol.config["scf"]["type"] == 'rhf':
        mf = dft.RKS(mole)
    else:
        print(f'\n   PyOQP only support pySCF with rhf, uhf, and rohf\n')
        sys.exit(1)

    functional = mol.config["input"]["functional"]
    try:
        xcfun = pyscf_functional[functional]
    except KeyError:
        xcfun = 'pbe0'

    mf.xc = xcfun
    mf.max_cycle = 0
    mf.level_shift = mol.config["scf"]["vshift"]

    try:
        os.environ["PYSCF_MAX_MEMORY"]
    except KeyError:
        mf.max_memory = 1000  # MB
    try:
        os.environ["OMP_NUM_THREADS"]
    except KeyError:
        lib.num_threads(1)

    pyscf_guess_key = 'atom' if guess_type == 'sad' else 'sap'
    try:
        dm0 = mf.get_init_guess(mole, key=pyscf_guess_key)
    except Exception as exc:
        if guess_type == 'sap':
            print(f'\n   PyOQP/PySCF SAP guess failed ({exc}); falling back to SAD\n')
            dm0 = mf.get_init_guess(mole, key='atom')
        else:
            raise

    # max_cycle=0 performs the initial Fock build/diagonalization from dm0,
    # giving the native exporter a consistent MO coefficient set to pass to OpenQP.
    mf.kernel(dm0=dm0)

    load_pyscf_guess(mf, mol)

    return

