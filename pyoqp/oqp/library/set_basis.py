"""Set up basis set for the molecule"""

import oqp
from oqp.utils.file_utils import try_basis, dump_log
import numpy as np
from oqp import ffi
import basis_set_exchange as bse

class BasisData:
    def __init__(self, mol):
        self.mol = mol
        self.shell_num = 0
        self.num_atoms = 0
        self.atoms = []
        self.basis_names = []
        self.shells_data = []
        self.ecp = {}
        self.ecp["ang"] = []
        self.ecp["r_expo"] = []
        self.ecp["g_expo"] = []
        self.ecp["coef"] = []
        self.ecp["coord"] = []
        self.ecp["element_id"] = []


    def get_basis_data(self, elements, basis_name='STO-3G',el_index=0):

        shells_data = []
        basis = bse.get_basis(basis_name, elements=elements)
        print(basis['elements'])

        for element_id in basis['elements']:
            for shell in basis['elements'][element_id]['electron_shells']:
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
                    "function_type": shell['function_type'],
                    "region": shell['region'],
                    "angular_momentum": ang_mom,
                    "exponents": list(map(float, [shell['exponents'][i] for i in indices])),
                    "coefficients": list(map(float, [coefficients[i] for i in indices]))
                }
#                    print("shellshell",shell_dict)
                    shells_data.append(shell_dict)
                    ang_ii+=1

        self.ecp["element_id"] = 0
        if 'ecp_potentials' in basis['elements'][element_id]:
            ecp_list = basis['elements'][element_id]['ecp_potentials']

            self.ecp["element_id"] +=1
            for term in ecp_list:
                self.ecp["ang"].extend((term['angular_momentum']*len(term['gaussian_exponents'])))
                self.ecp["r_expo"].extend(term['r_exponents'])
                self.ecp["g_expo"].extend(term['gaussian_exponents'])
                self.ecp["coef"].extend(term['coefficients'][0])
            self.ecp["coord"].extend(self.atoms[el_index,1:])
        else:
            self.ecp["ang"] = [0]
            self.ecp["r_expo"] = [0]
            self.ecp["g_expo"] = [0]
            self.ecp["coef"] = [0]
            self.ecp["coord"] = [0]
#        print("ECP:::::",self.ecp)

        return shells_data



    def get_basislist(self):

        self.basis_names = self.mol.config["input"]["basis"].split(',')

        if len(self.basis_names) == 1:
            self.basis_names = [self.basis_names[0]] * self.num_atoms
        return self.basis_names

    def get_num_atom(self):

        system = self.mol.config["input"]["system"].strip().split('\n')
        self.atoms = np.array([list(map(float, line.split())) for line in system])
        self.num_atoms = self.atoms.shape[0]
        return self.num_atoms

    def create_shell_data(self):

        num_atoms = self.get_num_atom()
        basis_list = self.get_basislist()
        for el_index in range(0, num_atoms):
            temp_shell = self.get_basis_data(int(self.atoms[el_index,0]), basis_list[el_index], el_index)
            self.shells_data.extend(temp_shell)
        return self.shells_data

    def set_ecp_data(self):


        self.mol.data["element_id"] = int(self.ecp["element_id"])
        self.mol.data["num_expo"] = len(self.ecp["g_expo"])
        self.mol.data["ecp_nam"] = len(self.ecp["ang"])

        expo_array = np.array(self.ecp["g_expo"], dtype=np.float64)
        self.mol.data["expo"] = ffi.cast("double*", ffi.from_buffer(expo_array))
        coef_array = np.array(self.ecp["coef"], dtype=np.float64)
        self.mol.data["coef"] = ffi.cast("double*", ffi.from_buffer(coef_array))

        r_expo_array = np.array(self.ecp["r_expo"], dtype=np.float64)
        self.mol.data["ecp_rex"] =  ffi.cast("int*", ffi.from_buffer(r_expo_array))
        coord_array = np.array(self.ecp["coord"], dtype=np.float64)
        
        self.mol.data["ecp_am"] = ffi.cast("int*", ffi.from_buffer(np.array(self.ecp["ang"], dtype=np.float64)))
        self.mol.data["ecp_coord"] = ffi.cast("double*", ffi.from_buffer(coord_array))
        oqp.append_ecp(self.mol)


    def set_basis_data(self):

        shells_data = self.create_shell_data()

        for shell in shells_data:
            self.mol.data["id"] = int(shell["id"])
            self.mol.data["element_id"] = int(shell["element_id"])
            self.mol.data["num_expo"] = len(shell["exponents"])
            self.mol.data["ang_mom"] = shell["angular_momentum"]

            expo_array = np.array(shell["exponents"], dtype=np.float64)
            self.mol.data["expo"] = ffi.cast("double*", ffi.from_buffer(expo_array))

            coef_array = np.array(shell["coefficients"], dtype=np.float64)
            self.mol.data["coef"] = ffi.cast("double*", ffi.from_buffer(coef_array))

            oqp.append_shell(self.mol)
        self.set_ecp_data()




def set_basis(mol):
    """Set up basis set for the molecule"""
#    basis_file = try_basis(mol.config["input"]["basis"])
#    mol.data["OQP::basis_filename"] = basis_file
#    print(mol.config["input"])

    basis_data= BasisData(mol)

    basis_data.set_basis_data()
    oqp.apply_basis(mol)

    if mol.data["basis_set_issue"]:
        dump_log(mol, title='PyOQP: basis set is not set properly', section='end')
        exit()
