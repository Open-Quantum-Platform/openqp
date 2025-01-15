"""Set up basis set for the molecule"""

import oqp
from oqp.utils.file_utils import try_basis, dump_log
import numpy as np
from oqp import ffi
import basis_set_exchange as bse
import json

class BasisData:
    def __init__(self, mol):
        self.mol = mol
        self.shell_num = 0
        self.num_atoms = mol.data["natom"]
        self.atoms = mol.data["qn"]
        self.atom_xyz = mol.data["xyz"]
        print(type(self.atom_xyz))
        print("xyz",list(self.atom_xyz[0][0:3]))
        self.basis_names = []
        self.shells_data = []
        self.ecp = {}
        self.ecp["ang"] = []
        self.ecp["r_expo"] = []
        self.ecp["g_expo"] = []
        self.ecp["coef"] = []
        self.ecp["coord"] = []
        self.ecp["element_id"] = []
        self.ecp["ecp_electron"] = []
        self.ecp["num_expo"] = []


    def get_basis_data(self, elements, basis_name='STO-3G',el_index=0):

        shells_data = []
        basis = self.read_basis_fmt(basis_name, elements)

        for shell in basis['elements'][str(elements)]['electron_shells']:
            ang_ii = 0

            for coefficients in shell['coefficients']:
                shell['angular_momentum']
                if len(shell['angular_momentum']) > 1:
                    ang_mom = shell['angular_momentum'][ang_ii]
                else:
                    ang_mom = shell['angular_momentum'][0]
                self.shell_num = self.shell_num + 1
                f_coefficients =  list(map(float, coefficients))
                indices = [i for i, value in enumerate(f_coefficients) if value != 0]
                shell_dict = {
                "id" : self.shell_num,
                "element_id":  el_index+1,
                "angular_momentum": ang_mom,
                "exponents": list(map(float, [shell['exponents'][i] for i in indices])),
                "coefficients": list(map(float, [coefficients[i] for i in indices]))
            }

                shells_data.append(shell_dict)
                ang_ii+=1

        if 'ecp_potentials' in basis['elements'][str(elements)]:
            ecp_list = basis['elements'][str(elements)]['ecp_potentials']
            self.ecp["ecp_electron"].extend([basis['elements'][str(elements)]['ecp_electrons']])
            self.ecp["element_id"] +=1
            ecp_num_expo = 0
            for term in ecp_list:
                ecp_num_expo += len(term['r_exponents'])
                self.ecp["ang"].extend((term['angular_momentum']*len(term['gaussian_exponents'])))
                self.ecp["r_expo"].extend(term['r_exponents'])
                self.ecp["g_expo"].extend(term['gaussian_exponents'])
                self.ecp["coef"].extend(term['coefficients'][0])

            self.ecp["num_expo"].extend([ecp_num_expo])
            self.ecp["coord"].extend(list(self.atom_xyz[0][0:3]))
        else:
            self.ecp["ecp_electron"].extend([0])

        return shells_data



    def get_basislist(self):

        self.basis_names = self.mol.config["input"]["basis"].split(',')


        if len(self.basis_names) == 1:
            self.basis_names = [self.basis_names[0]] * self.num_atoms
        return self.basis_names

    def create_shell_data(self):

        num_atoms = self.num_atoms
        basis_list = self.get_basislist()

        self.ecp["element_id"] = 0
        for el_index in range(0, num_atoms):
            element_shells = self.get_basis_data(int(self.atoms[el_index]), basis_list[el_index], el_index)
            self.shells_data.extend(element_shells)
        return self.shells_data

    def set_ecp_data(self):


        self.mol.data["element_id"] = int(self.ecp["element_id"])
        self.mol.data["ecp_nam"] = len(self.ecp["ang"])

        n_expo_array = np.array(self.ecp["num_expo"], dtype=np.int32)
        self.mol.data["num_expo"] = ffi.cast("int*", ffi.from_buffer(n_expo_array))

        expo_array = np.array(self.ecp["g_expo"], dtype=np.float64)
        self.mol.data["expo"] = ffi.cast("double*", ffi.from_buffer(expo_array))
        coef_array = np.array(self.ecp["coef"], dtype=np.float64)
        self.mol.data["coef"] = ffi.cast("double*", ffi.from_buffer(coef_array))

        r_expo_array = np.array(self.ecp["r_expo"], dtype=np.float64)
        self.mol.data["ecp_rex"] =  ffi.cast("int*", ffi.from_buffer(r_expo_array))
        coord_array = np.array(self.ecp["coord"], dtype=np.float64)

        self.mol.data["ecp_am"] = ffi.cast("int*", ffi.from_buffer(np.array(self.ecp["ang"], dtype=np.float64)))
        self.mol.data["ecp_zn"] = ffi.cast("int*", ffi.from_buffer(np.array(self.ecp["ecp_electron"], dtype=np.int32)))
        self.mol.data["ecp_coord"] = ffi.cast("double*", ffi.from_buffer(coord_array))
        oqp.append_ecp(self.mol)


    def set_basis_data(self):

        shells_data = self.create_shell_data()


        for shell in shells_data:
            self.mol.data["id"] = int(shell["id"])
            self.mol.data["element_id"] = int(shell["element_id"])
            self.mol.data["ang_mom"] = shell["angular_momentum"]

            n_expo_array = np.array([len(shell["exponents"])], dtype=np.int32)
            self.mol.data["num_expo"] = ffi.cast("int*", ffi.from_buffer(n_expo_array))

            expo_array = np.array(shell["exponents"], dtype=np.float64)
            self.mol.data["expo"] = ffi.cast("double*", ffi.from_buffer(expo_array))

            coef_array = np.array(shell["coefficients"], dtype=np.float64)
            self.mol.data["coef"] = ffi.cast("double*", ffi.from_buffer(coef_array))

            oqp.append_shell(self.mol)

        self.set_ecp_data()

    def read_basis_fmt(self, url: str, elements: str) -> dict:
        """
        Reads and returns the 'elements' dictionary from a basis set specification.

        :param url: A string indicating either a local file path ('file:/path/...').
        :param basis_name: Optional name of the basis set, used by bse.get_basis if not a local file.
        :param elements:
        :return: Dictionary of basis data of element.
        """
#        url = "file:631g.json"
        if url.startswith("file:"):
            file_path = url[len("file:"):]

            if file_path.endswith('.json'):
                with open(file_path, "r") as f:
                    basis_data = json.load(f)
            else:
                try:
                    basis_data = bse.read_formatted_basis_file(file_path)
                except Exception as e:
                    raise ValueError(f"Could not read basis file '{file_path}': {e}")

            if 'elements' in basis_data:
                return basis_data
            else:
                raise KeyError(f"No 'elements' key found in basis data for '{file_path}'")
        else:
            basis_data = bse.get_basis(url, elements=elements)
            return basis_data


def set_basis(mol):
    """Set up basis set for the molecule"""
    basis_file = mol.config["input"]["basis"]
    mol.data["OQP::basis_filename"] = basis_file

    basis_data= BasisData(mol)

    basis_data.set_basis_data()
    oqp.apply_basis(mol)

    if mol.data["basis_set_issue"]:
        dump_log(mol, title='PyOQP: basis set is not set properly', section='end')
        exit()
