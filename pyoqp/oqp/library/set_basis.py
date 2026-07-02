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
from oqp.molecule.oqpdata import compute_alpha_beta_electrons, ispher_mode


class BasisData:
    def __init__(self, mol):
        self.mol = mol
        self.shell_num = 0
        self.mpi_manager = MPIManager()
        self.num_atoms = mol.data["natom"]
        self.atoms = mol.data["qn"]
        self.atom_xyz = mol.data["xyz"]
        self.basis_names = []
        self.shells_data = []
        self._ffi_buffer_refs = []
        self.mol._ffi_buffer_refs = self._ffi_buffer_refs
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
        self.use_ecp = False

    def _ispher_mode(self):
        return ispher_mode(self.mol.config['input'].get('ispher', 'auto'))

    def _buffer_ptr(self, values, dtype, ctype):
        array = np.ascontiguousarray(np.asarray(values, dtype=dtype))
        if not hasattr(self, "_ffi_buffer_refs"):
            self._ffi_buffer_refs = []
        if getattr(self.mol, "_ffi_buffer_refs", None) is not self._ffi_buffer_refs:
            self.mol._ffi_buffer_refs = self._ffi_buffer_refs
        self._ffi_buffer_refs.append(array)
        return ffi.cast(ctype, ffi.from_buffer(array))

    def get_basis_data(self, elements, basis_name='STO-3G', el_index=0):
        """
        Retrieves and organizes electron shell data (and ECP data if present) for a specified element.

        :param elements: Atomic number (int) of the element.
        :param basis_name: Name of the basis set (default: 'STO-3G').
        :param el_index: Index of the element in the molecule.
        :return: A list of dictionaries containing shell data.
        """

        shells_data = []
        if self.mol.usempi:
            basis = self.read_basis_fmt(basis_name, elements) if self.mpi_manager.rank == 0 else None
            basis = self.mpi_manager.bcast(basis, root=0)
        else:
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

            # BSE tags the harmonic type per shell: 'gto_spherical' |
            # 'gto_cartesian' | 'gto' (s/p, where Cartesian == spherical).
            # ispher=auto selects the AO convention the basis was published
            # with, e.g. 6-31G* -> Cartesian (6d), cc-pVDZ/def2 -> spherical
            # (5d); ispher=true forces pure spherical for every shell (GAMESS
            # ISPHER=1); ispher=false deactivates the harmonic gate globally.
            shell_ft = shell.get('function_type', 'gto')
            shell_is_spherical = shell_ft.endswith('spherical')
            force_spherical = self._ispher_mode() == 'true'

            for coefficients in shell['coefficients']:
                shell['angular_momentum']
                if len(shell['angular_momentum']) > 1:
                    ang_mom = shell['angular_momentum'][ang_ii]
                else:
                    ang_mom = shell['angular_momentum'][0]
                shell_harmonic = 1 if (shell_is_spherical or force_spherical) else 0

                self.shell_num += 1

                float_coeffs = list(map(float, coefficients))
                nonzero_indices = [i for i, value in enumerate(float_coeffs) if value != 0]
                shell_dict = {
                    "id": self.shell_num,
                    "element_id": el_index + 1,
                    "angular_momentum": ang_mom,
                    "harmonic": shell_harmonic,
                    "exponents": [float(x) for i, x in enumerate(shell['exponents']) if i in nonzero_indices],
                    "coefficients": list(map(float, [coefficients[i] for i in nonzero_indices]))
                }

                shells_data.append(shell_dict)
                ang_ii += 1

        if 'ecp_potentials' in basis['elements'][element_key]:
            self.use_ecp = True
            self.ecp["ecp_electron"].append(basis['elements'][element_key]['ecp_electrons'])
            self.ecp["element_id"] += 1

            ecp_list = basis['elements'][element_key]['ecp_potentials']
            ecp_num_expo = 0
            for term in ecp_list:
                ecp_num_expo += len(term['r_exponents'])
                self.ecp["ang"].extend((term['angular_momentum'] * len(term['gaussian_exponents'])))
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
        if self.mol.config["input"]["basis"] == 'library':
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
                    raise FileNotFoundError(f"Please correctly add a tag for each atom. ({parts})")
                    break
                basis_tags.append(parts[4])

            basis_list_str = self.mol.config["input"]["library"].strip()
            if not basis_list_str:
                raise FileNotFoundError("Please ensure the necessary library for tags is correctly added.")
            basis_dict = {}
            for line in basis_list_str.splitlines():
                library_parts = line.split()
                if len(library_parts) >= 2:
                    key = library_parts[0]
                    basis_info = " ".join(library_parts[1:])
                    basis_dict[key] = basis_info

            self.basis_names = [basis_dict.get(tag, "UNKNOWN") for tag in basis_tags]
            return self.basis_names

        self.basis_names = self.mol.config["input"]["basis"].split(';')
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
            element_shells = self.get_basis_data(element, basis_name, el_index)
            self.shells_data.extend(element_shells)
        return self.shells_data

    def set_ecp_data(self):
        """
        Sets the ECP-related data in the molecule's data structure and appends it via oqp.
        """

        self.mol.data["element_id"] = int(self.ecp["element_id"])
        self.mol.data["ecp_nam"] = len(self.ecp["ang"])

        self.mol.data["num_expo"] = self._buffer_ptr(self.ecp["num_expo"], np.int32, "int*")
        self.mol.data["expo"] = self._buffer_ptr(self.ecp["g_expo"], np.float64, "double*")
        self.mol.data["coef"] = self._buffer_ptr(self.ecp["coef"], np.float64, "double*")
        self.mol.data["ecp_rex"] = self._buffer_ptr(self.ecp["r_expo"], np.float64, "int*")
        self.mol.data["ecp_am"] = self._buffer_ptr(self.ecp["ang"], np.float64, "int*")
        self.mol.data["ecp_zn"] = self._buffer_ptr(self.ecp["ecp_electron"], np.int32, "int*")
        self.mol.data["ecp_coord"] = self._buffer_ptr(self.ecp["coord"], np.float64, "double*")

        oqp.append_ecp(self.mol)

        if self.use_ecp:

            molecule = self.mol.data._data
            natom = molecule.mol_prop.natom
            charge = molecule.mol_prop.charge
            nelec_base = sum(int(molecule.qn[i]) for i in range(natom)) - charge
            ecp_all_n = 0
            ecp_all_n = sum(int(v) for v in self.ecp["ecp_electron"])

            nelec_qm = nelec_base - ecp_all_n
            if nelec_qm < 0:
                raise ValueError(f"ECP removes more electrons than available: nelec_qm={nelec_qm}")

            mult = molecule.mol_prop.mult
            na, nb = compute_alpha_beta_electrons(nelec_qm, mult)

            molecule.mol_prop.nelec     = nelec_qm
            molecule.mol_prop.nelec_A   = na
            molecule.mol_prop.nelec_B   = nb
            molecule.mol_prop.nocc      = max(na, nb)

    def set_basis_data(self):
        """
        Sets the basis data for the molecule's shells, then sets the ECP data.
        """

        shells_data = self.create_shell_data()

        for shell in shells_data:
            self.mol.data["id"] = int(shell["id"])
            self.mol.data["element_id"] = int(shell["element_id"])
            self.mol.data["ang_mom"] = shell["angular_momentum"]
            self.mol.data["harmonic"] = int(shell.get("harmonic", 0))

            self.mol.data["num_expo"] = self._buffer_ptr([len(shell["exponents"])], np.int32, "int*")
            self.mol.data["expo"] = self._buffer_ptr(shell["exponents"], np.float64, "double*")
            self.mol.data["coef"] = self._buffer_ptr(shell["coefficients"], np.float64, "double*")

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
            file_name = url[len("file:"):]
            directory = os.path.dirname(self.mol.input_file)
            file_path = os.path.join(directory, file_name)

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

    basis_data = BasisData(mol)
    basis_data.set_basis_data()

    oqp.apply_basis(mol)

    if mol.data["basis_set_issue"]:
        dump_log(mol, title='PyOQP: basis set is not set properly', section='end')
        exit()
