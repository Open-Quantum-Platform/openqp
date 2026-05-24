"""OQP external quantum chemical program"""
import re
import os
import sys
import oqp
from oqp.utils.constants import ANGSTROM_TO_BOHR
from oqp.periodic_table import SYMBOL_MAP, ELEMENTS_NAME

try:
    from pyscf import gto, dft, lib

except ModuleNotFoundError:
    gto = None
    dft = None
    lib = None

try:
    from mokit.lib import py2openqp

except ModuleNotFoundError:
    py2openqp = None

pyscf_functional = {
    "pbe0": "pbe0",
    "b3lyp": "b3lyp",
    "bhhlyp": "bhandhlyp",
    "cam-b3lyp": "camb3lyp",
    "camb3lyp": "camb3lyp",
}


def guess_from_pyscf(mol):
    """Set up pySCF calculation"""

    if gto is None:
        print(f'\n   PyOQP cannot find PySCF\n')
        ## pyscf run mpi with a different lib, which might cause conflict with mpi4py, disabling mpi now for future dev
        sys.exit(1)

    if py2openqp is None:
        print(f'\n   PyOQP cannot find Mokit\n')
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

    use_mrsf_export = (
        mol.config["input"]["method"] == 'tdhf'
        and mol.config["tdhf"]["type"] in {'mrsf', 'umrsf'}
    )
    py2openqp(mf, 'mokit.inp', mrsf=use_mrsf_export)

    target_basis = mol.config['input']['basis']
    target_library = mol.config['input']['library']
    mol.config['guess']['file'] = 'mokit_wfn.json'
    mol.load_data()

    mol.data.set_scf_active_basis(1)
    mol.config['input']['basis'] = mol.config['json']['basis']
    mol.config['input']['library'] = mol.config['json']['library']
    oqp.library.set_basis(mol)

    mol.data.set_scf_active_basis(0)
    mol.config['input']['basis'] = target_basis
    mol.config['input']['library'] = target_library
    oqp.library.project_basis(mol)

    return

def guess_from_pyscf_initial_density(mol, guess_type):
    """Build a PySCF atomic-density/potential initial guess for OpenQP.

    The resulting PySCF orbitals are exported through the same Mokit bridge used
    by ``guess.type=pyscf``. ``sad`` uses PySCF's atomic-density guess. ``sap``
    requests PySCF's superposition-of-atomic-potentials guess when available.
    """

    if gto is None:
        print(f'\n   PyOQP cannot find PySCF\n')
        sys.exit(1)

    if py2openqp is None:
        print(f'\n   PyOQP cannot find Mokit\n')
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
    # giving Mokit a consistent MO coefficient set to export to OpenQP.
    mf.kernel(dm0=dm0)

    use_mrsf_export = (
        mol.config["input"]["method"] == 'tdhf'
        and mol.config["tdhf"]["type"] in {'mrsf', 'umrsf'}
    )
    py2openqp(mf, 'mokit.inp', mrsf=use_mrsf_export)

    target_basis = mol.config['input']['basis']
    target_library = mol.config['input']['library']
    mol.config['guess']['file'] = 'mokit_wfn.json'
    mol.load_data()

    mol.data.set_scf_active_basis(1)
    mol.config['input']['basis'] = mol.config['json']['basis']
    mol.config['input']['library'] = mol.config['json']['library']
    oqp.library.set_basis(mol)

    mol.data.set_scf_active_basis(0)
    mol.config['input']['basis'] = target_basis
    mol.config['input']['library'] = target_library
    oqp.library.project_basis(mol)

    return

