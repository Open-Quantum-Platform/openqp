"""Set up basis set for the molecule"""

import oqp
from oqp.utils.file_utils import try_basis, dump_log
import numpy as np
import cffi

def get_basis_data(elements, basis_name='STO-3G', shell_num=0):
    """Get basis set data for the given elements using the Basis Set Exchange library."""
    import basis_set_exchange as bse

    shells_data = []
    basis = bse.get_basis(basis_name, elements=elements)

    for element_id in basis['elements']:
        for shell in basis['elements'][element_id]['electron_shells']:
            exponents = list(map(float, shell['exponents']))

            for ang_mom, coefficients in zip(shell['angular_momentum'], shell['coefficients']):
                shell_num = shell_num + 1
                shell_dict = {
                    "id" : shell_num,
                    "element_id": element_id,
                    "function_type": shell['function_type'],
                    "region": shell['region'],
                    "angular_momentum": ang_mom,
                    "exponents": exponents,
                    "coefficients": list(map(float, coefficients))
                }
                shells_data.append(shell_dict)

    return shells_data, shell_num

def set_basis_data(mol):

    ffi = cffi.FFI()
    shells_data,shell_num = get_basis_data(["O"], mol.config["input"]["basis"])
    temp,shell_num = get_basis_data(["H"], mol.config["input"]["basis"],shell_num)
    shells_data.extend(temp)
    temp,shell_num = get_basis_data(["H"], mol.config["input"]["basis"],shell_num)
    shells_data.extend(temp)
    print(shells_data)

    for shell in shells_data:
        mol.data["id"] = int(shell["id"])
        mol.data["element_id"] = int(shell["element_id"])
        mol.data["num_expo"] = len(shell["exponents"])
        mol.data["ang_mom"] = shell["angular_momentum"]

        expo_array = np.array(shell["exponents"], dtype=np.float64)
        mol.data["expo"] = ffi.cast("double*", ffi.from_buffer(expo_array))

        coef_array = np.array(shell["coefficients"], dtype=np.float64)
        mol.data["coef"] = ffi.cast("double*", ffi.from_buffer(coef_array))

        oqp.append_shell(mol)


def set_basis(mol):
    """Set up basis set for the molecule"""
    basis_file = try_basis(mol.config["input"]["basis"])
    mol.data["OQP::basis_filename"] = basis_file

    set_basis_data(mol)
    oqp.apply_basis(mol)

    if mol.data["basis_set_issue"]:
        dump_log(mol, title='PyOQP: basis set is not set properly', section='end')
        exit()
