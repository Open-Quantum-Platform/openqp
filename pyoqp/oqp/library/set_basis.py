"""Set up basis set for the molecule"""

import oqp
import os
import sys
from oqp.utils.file_utils import try_basis, dump_log
import numpy as np
from oqp import ffi
import basis_set_exchange as bse
from oqp.utils.mpi_utils import MPIManager
import json

class BasisData:
    def __init__(self, mol):
        self.mol = mol
        self.shell_num = 0
        self.num_atoms = mol.data["natom"]
        self.atoms = mol.data["qn"]
        self.atom_xyz = mol.data["xyz"]
        self.basis_names = []
        self.shells_data = []
        self.ecp = {
            "ang": [],
            "r_expo": [],
            "g_expo": [],
            "coef": [],
            "coord": [],
            "element_id": [],
            "ecp_electron": [],
            "num_expo": []
        }

    def get_basis_data(self, elements, basis_name='STO-3G',el_index=0):
        """
        Retrieves and organizes electron shell data (and ECP data if present) for a specified element.

        :param elements: Atomic number (int) of the element.
        :param basis_name: Name of the basis set (default: 'STO-3G').
        :param el_index: Index of the element in the molecule.
        :return: A list of dictionaries containing shell data.
        """

        shells_data = []
        basis = self.read_basis_fmt(basis_name, elements)
        element_key = str(elements)
        if 'electron_shells' not in basis['elements'].get(element_key, {}):
            print(
                f"Error: The basis set '{basis_name}' does not support the element '{element_key}'.\n"
                "Please choose a valid basis set. For more information, visit:\n"
                "www.basissetexchange.org\n"
                "Aborting."
            )
            sys.exit(1)

        for shell in basis['elements'][element_key]['electron_shells']:
            ang_ii = 0

            for coefficients in shell['coefficients']:
                shell['angular_momentum']
                if len(shell['angular_momentum']) > 1:
                    ang_mom = shell['angular_momentum'][ang_ii]
                else:
                    ang_mom = shell['angular_momentum'][0]

                self.shell_num += 1

                float_coeffs =  list(map(float, coefficients))
                nonzero_indices = [i for i, value in enumerate(float_coeffs) if value != 0]
                shell_dict = {
                "id" : self.shell_num,
                "element_id":  el_index+1,
                "angular_momentum": ang_mom,
                "exponents": [float(x) for i, x in enumerate(shell['exponents']) if i in nonzero_indices],
                "coefficients": list(map(float, [coefficients[i] for i in nonzero_indices]))
            }

                shells_data.append(shell_dict)
                ang_ii += 1

        if 'ecp_potentials' in basis['elements'][element_key]:
            self.ecp["ecp_electron"].append(basis['elements'][element_key]['ecp_electrons'])
            self.ecp["element_id"] +=1

            ecp_list = basis['elements'][element_key]['ecp_potentials']
            ecp_num_expo = 0
            for term in ecp_list:
                ecp_num_expo += len(term['r_exponents'])
                self.ecp["ang"].extend((term['angular_momentum']*len(term['gaussian_exponents'])))
                self.ecp["r_expo"].extend(term['r_exponents'])
                self.ecp["g_expo"].extend(term['gaussian_exponents'])
                self.ecp["coef"].extend(term['coefficients'][0])

            self.ecp["num_expo"].append(ecp_num_expo)
            self.ecp["coord"].extend(list(self.atom_xyz[el_index][0:3]))
        else:
            self.ecp["ecp_electron"].extend([0])

        return shells_data



    def get_basislist(self):
        """
        Retrieves a list of basis names for structure.

        :return: A list of basis names (strings).
        """
        basis_tags = []
        system = self.mol.config["input"]["system"]
        system = system.split("\n")
        if system[0]:
            if not os.path.exists(system[0]):
                raise FileNotFoundError("XYZ file %s is not found!" % system[0])

            with open(system[0], 'r') as xyzfile:
                system = xyzfile.read().splitlines()

            num_atoms = int(system[0])
            system = system[2: 2 + num_atoms]
        else:
            system = system[1:]

        for i, line in enumerate(system):
            parts = line.split()
            if len(parts) < 5:
                basis_tags.clear()
                break
            basis_tags.append(parts[4])

        if basis_tags:
            basis_list_str = self.mol.config["basis_set"]["library"].strip()
            if basis_list_str:
                basis_dict = {}
                for line in basis_list_str.splitlines():
                    library_parts = line.split()
                    if len(library_parts) >= 2:
                        key = library_parts[0]
                        basis_info = " ".join(library_parts[1:])
                        basis_dict[key] = basis_info

                self.basis_names = [basis_dict.get(tag, "UNKNOWN") for tag in basis_tags]
                return self.basis_names


        self.basis_names = self.mol.config["input"]["basis"].split(',')
        if len(self.basis_names) == 1:
            self.basis_names = [self.basis_names[0]] * self.num_atoms

        return self.basis_names

    def create_shell_data(self):
        """
        Creates the shells data for each atom in the molecule based on provided basis sets.

        :return: A list of all shells data across all atoms.
        """

        basis_list = self.get_basislist()

        self.ecp["element_id"] = 0

        for el_index in range(self.num_atoms):
            element = int(self.atoms[el_index])
            basis_name = basis_list[el_index]
            element_shells = self.get_basis_data(element, basis_name,el_index)
            self.shells_data.extend(element_shells)
        return self.shells_data

    def set_ecp_data(self):
        """
        Sets the ECP-related data in the molecule's data structure and appends it via oqp.
        """

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
        """
        Sets the basis data for the molecule's shells, then sets the ECP data.
        """

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

    if(MPIManager().rank == 0):
        basis_data= BasisData(mol)
        basis_data.set_basis_data()

    oqp.apply_basis(mol)

    if mol.data["basis_set_issue"]:
        dump_log(mol, title='PyOQP: basis set is not set properly', section='end')
        exit()
