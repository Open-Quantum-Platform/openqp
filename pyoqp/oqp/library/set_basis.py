"""Set up basis set for the molecule"""

import oqp
from oqp.utils.file_utils import try_basis, dump_log

def set_basis(mol):
    """Set up basis set for the molecule"""
    basis_file = try_basis(mol.config["input"]["basis"])
    mol.data["OQP::basis_filename"] = basis_file

    oqp.apply_basis(mol)

    if mol.control.basis_set_issue:
        dump_log(mol, title='PyOQP: basis set is not set properly', section='end')
        exit()
