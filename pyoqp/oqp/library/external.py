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
    mole.basis = mol.config["input"]["basis"]
    mole.charge = mol.config["input"]["charge"]
    mole.spin = mol.config["scf"]["multiplicity"] - 1
    mole.output = '%s.pyscf' % mol.project_name
    mole.verbose = 4
    mole.build(cart=False)
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

    py2openqp(mf, 'mokit.inp', mrsf=True)

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
