import oqp
import numpy as np
from oqp.utils.qmmm import openmm_potential,openmm_energy

def ints_1e(mol):
    """Compute a set of one-electron integrals"""
    if mol.config['input']['qmmm_flag']:
       current_xyz = mol.get_system().reshape((-1, 3))
       mol.data["OQP::mm_potential"]=np.array(openmm_potential(current_xyz)).tolist()
       mol.data["OQP::mm_energy"]=openmm_energy()
    oqp.int1e(mol)
