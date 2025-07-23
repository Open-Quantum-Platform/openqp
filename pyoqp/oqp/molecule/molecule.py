"""OQP molecule class"""
import os
import copy
import json
import platform
import numpy as np
import oqp
from oqp.utils.input_parser import OQPConfigParser
from oqp.molden.moldenwriter import MoldenWriter
from .oqpdata import OQPData, OQP_CONFIG_SCHEMA
from oqp.utils.mpi_utils import MPIManager
from oqp.utils.mpi_utils import mpi_get_attr, mpi_dump
from oqp import ffi
class Molecule:
    """
    OQP molecule representation in python
    """

    def __init__(self, project_name, input_file, log,
                 xyz=None, elem=None, mass=None, charge=0, mult=1, silent=0, idx=1):
        self.mpi_manager = MPIManager()
        self.usempi = True
        self.silent = silent
        self.idx = idx

        self.xyz = xyz
        self.elem = elem
        self.mass = mass
        self.charge = charge
        self.mult = mult

        self.control = None
        self.mol_energy = None

        self.data = None
        self.data_allocate()

        self.config = {}
        self.project_name = project_name
        self.input_file = input_file
        self.log = log
        self.log_path = os.path.dirname(log)
        self.energies = None
        self.grads = None
        self.dcm = []  # Nstate, Nstate
        self.nac = []  # Npairs, 3, Natom,
        self.soc = []  # Npairs, 1,
        self.freqs = np.zeros(0)  # 3Natom-6
        self.hessian = np.zeros(0)  # 3Natom, 3Natom
        self.modes = np.zeros(0)  # 3Natom-6, 3Natom
        self.inertia = np.zeros(0)  # 3

        self.tag = [
            'OQP::DM_A', 'OQP::DM_B',
            'OQP::FOCK_A', 'OQP::FOCK_B',
            'OQP::E_MO_A', 'OQP::E_MO_B',
            'OQP::VEC_MO_A', 'OQP::VEC_MO_B',
            'OQP::Hcore', 'OQP::SM', 'OQP::TM', 'OQP::WAO',
            'OQP::td_abxc', 'OQP::td_bvec_mo', 'OQP::td_mrsf_density', 'OQP::td_energies',
            'OQP::td_states_overlap',
            'OQP::dc_matrix', 'OQP::nac_matrix',
        ]
        self.skip_tag = {"rhf": ['OQP::DM_B', 'OQP::FOCK_B', 'OQP::E_MO_B', 'OQP::VEC_MO_B'],
                         "rohf": [],
                         "uhf": []
                         }
        self.config_tag = {
            'json': ['scf_type', 'basis', 'library']
        }
        self.start_time = None
        self.back_door = None

    def get_atoms(self):
        """
        Get read-only atoms
        """
        natom = self.data["natom"]
        atoms = np.frombuffer(
            oqp.ffi.buffer(self.elem,
                           natom * oqp.ffi.sizeof("double"))
        ).astype(int)

        return copy.deepcopy(atoms)

    def get_mass(self):
        """
        Get read-only molar mass
        """
        natom = self.data["natom"]
        atoms = np.frombuffer(
            oqp.ffi.buffer(self.mass,
                           natom * oqp.ffi.sizeof("double"))
        ).astype(float)

        return copy.deepcopy(atoms)

    def get_system(self):
        """
        Get read-only coordinates
        """
        natom = self.data['natom']
        coord = np.frombuffer(
            oqp.ffi.buffer(self.xyz, 3 * natom * oqp.ffi.sizeof("double")),
            dtype=np.double)

        return copy.deepcopy(coord)

    def get_grad(self):
        """
        Get gradient in Hartree/Bohr
        """
        natom = self.data['natom']
        grad = np.frombuffer(
            oqp.ffi.buffer(self.data._data.grad, 3 * natom * oqp.ffi.sizeof("double"))
        )

        return copy.deepcopy(grad)

    def get_nac(self):
        """
        Get non-adiabatic couping in Hartree/Bohr
        """

        return []

    def get_soc(self):
        """
        Get spin-orbit coupling in cm-1
        """

        return []

    def get_hess(self):
        """
        Get hessian results
        """

        return []

    def get_data(self):
        """
        Extract data from mol to dict
        """
        scf_type = self.config['scf']['type']
        data = {}
        for key in self.tag:
            if key in self.skip_tag[scf_type]:
                continue 
            try:
                data[key] = np.array(self.data[key]).tolist()

            except AttributeError:
                continue

        return data

    def get_data_from_back_door(self):
        """
        Extract mol data for nacme calculation
        """
        if isinstance(self.back_door, tuple):
            return self.back_door
        else:
            # previous data is not available, return current data to bypass nacme calculation
            return self.get_system(), self.get_data()

    def get_results(self):
        """
        Collect computed results to dict
        """
        data = {
            'atoms': self.get_atoms().tolist(),
            'coord': self.get_system().tolist(),
            'energy': self.mol_energy.energy,
        }

        # save td energies if available
        try:
            data['td_energies'] = np.array(self.data['OQP::td_energies']
                                           ).tolist()
        except AttributeError:
            data['td_energies'] = np.array([0]).tolist()

        # save gradients if available
        data['grad'] = np.array(self.get_grad()).tolist()
        data['nac'] = np.array(self.get_nac()).tolist()
        data['soc'] = np.array(self.get_soc()).tolist()
        data['hess'] = np.array(self.get_hess()).tolist()

        return data

    @mpi_get_attr
    def get_coord(self, coordinates):
        return coordinates

    def update_system(self, coordinates):
        """
        Modify coordinates in memory
        """
        coordinates = self.get_coord(coordinates)
        coordinates = coordinates.reshape((-1, 3))
        natom = self.data['natom']
        coord = np.frombuffer(
            oqp.ffi.buffer(self.xyz, 3 * natom * oqp.ffi.sizeof("double")),
            dtype=np.double).reshape((natom, 3))
        for at in range(len(coordinates)):
            for c in range(3):
                coord[at, c] = np.float64(coordinates[at, c])

    def update_mol(self, ref_mol):
        """
        Pass data from ref_mol to current mol
        """
        for key in ref_mol.tag:
            try:
                self.data[key] = copy.deepcopy(ref_mol.data[key])

            except AttributeError:
                continue

    def check(self, info):
        """
        Check internal data
        """
        if self.data._data.qn[0] != self.elem[0]:
            raise ValueError(info, self.data._data.qn[0], self.elem[0],
                             'var changed!')
        else:
            print(info, 'var checked!')

    def data_allocate(self):
        """Allocate new oqp data object"""
        if not self.data:
            self.data = OQPData(silent=self.silent)

    def data_deallocate(self):
        """Deallocate oqp data object"""
        self.data = None

    @mpi_get_attr
    def get_config(self, input_source):
        parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA, allow_no_value=True)

        # Determine the type of the input source and process accordingly
        if isinstance(input_source, str):  # Assuming input is a filename
            parser.read(input_source)
        elif isinstance(input_source, dict):  # Assuming input is a dictionary
            parser.load_dict(input_source)
        else:
            raise ValueError("Input must be a filename (str) or a configuration dictionary (dict)")

        # Print configuration if not in silent mode
        if not self.silent:
            parser.print_config()

        # Validate the configuration and apply it
        config = parser.validate()

        return config

    def load_config(self, input_source):
        """
        Load calculation parameters from a file or a dictionary based on the input type.

        input_source: filename (str) or config dictionary (dict)
        """
        self.mpi_manager.set_mpi_comm(self.data)
        self.config = self.get_config(input_source)
        self.data.apply_config(self.config)
        self.data['usempi'] = int(self.usempi)
        self.xyz = self.data._data.xyz
        self.elem = self.data._data.qn
        self.mass = self.data._data.mass
        self.mol_energy = self.data._data.mol_energy

        return self

    @mpi_dump
    def write_molden(self, filename):
        """Write calculation results in Molden format"""

        with open(filename, mode='w', encoding='ascii') as fout:
            basis = self.data.get_basis()
            nat = self.data['natom']
            nbf = basis['nbf']
            mdw = MoldenWriter(fout)
            mdw.write_atoms(nat, self.elem, self.xyz, angstrom=False)
            mdw.write_basis(nat, basis)

            if self.config['scf']['type'] == 'rhf':
                # alpha only
                orbitals = self.data['OQP::VEC_MO_A'].reshape([nbf, nbf])
                eorbitals = self.data['OQP::E_MO_A']
                nocc = self.data['nocc']
                occupancies = (2.0 if i < nocc else 0.0 for i in range(nbf))
                mdw.write_mo(basis, orbitals, eorbitals,
                             occupancies, spin='Alpha')
            else:
                # alpha
                orbitals = self.data['OQP::VEC_MO_A'].reshape([nbf, nbf])
                eorbitals = self.data['OQP::E_MO_A']
                nocc = self.data['nelec_A']
                occupancies = (1.0 if i < nocc else 0.0 for i in range(nbf))
                mdw.write_mo(basis, orbitals, eorbitals,
                             occupancies, spin='Alpha')
                # beta
                orbitals = self.data['OQP::VEC_MO_B'].reshape([nbf, nbf])
                eorbitals = self.data['OQP::E_MO_B']
                nocc = self.data['nelec_B']
                occupancies = (1.0 if i < nocc else 0.0 for i in range(nbf))
                mdw.write_mo(basis, orbitals, eorbitals,
                             occupancies, spin='Beta', header=False)

    def set_log(self):
        """
        Set up log file
        """
        if not self.log:
            if platform.uname()[0] == "Windows":
                log_file = b'NUL'
            elif platform.uname()[0] == "Linux":
                log_file = b'/dev/null'
            elif platform.uname()[0] == "Darwin":
                log_file = b'/dev/null'
            else:
                log_file = b'/dev/null'
        else:
            log_file = bytes(str(self.log), encoding='ascii')

        _log_c = oqp.ffi.new("char[]", log_file)
        log_c = oqp.ffi.new("struct Cstring *", [len(log_file), _log_c])

        return log_c

    def set_config_json(self):
        data = {}
        data['json'] = {
            'scf_type': self.config['scf']['type'],
            'basis': self.config['input']['basis'],
            'library': self.config['input']['library']
        }
        return data

    @mpi_dump
    def save_data(self):
        """
        Save mol data and computed results to json
        """
        if self.idx != 1:
            jsonfile = self.log.replace('.log', f'_{self.idx}.json')
        else:
            jsonfile = self.log.replace('.log', '.json')
        data = self.get_data()
        data.update(self.get_results())
        data.update(self.set_config_json())

        with open(jsonfile, 'w') as outdata:
            json.dump(data, outdata, indent=2)

    @mpi_dump
    def save_freqs(self, state):
        jsonfile = self.log.replace('.log', '.hess.json')
        data = {
            'atoms': self.get_atoms().tolist(),
            'coord': self.get_system().tolist(),
            'mass': self.get_mass().tolist(),
            'energy': self.energies[state],
            'hessian': self.hessian.tolist(),
            'freqs': self.freqs.tolist(),
            'modes': self.modes.tolist(),
            'inertia': self.modes.tolist(),
        }

        with open(jsonfile, 'w') as outdata:
            json.dump(data, outdata, indent=2)

    def load_data(self):
        # load data from json to mol
        guess_geom = self.config['guess']['continue_geom']
        guess_file = self.config['guess']['file']

        if not os.path.exists(guess_file):
            exit(f'mol object {guess_file} does not exist')

        with open(guess_file, 'r') as indata:
            data = json.load(indata)

        in_atoms = self.get_atoms()
        ld_atoms = np.array(data['atoms'])

        if len(in_atoms) != len(ld_atoms):
            exit('loading data from json, the number of atoms does not match!')

        if np.amax(np.abs(in_atoms - ld_atoms)) > 0:
            exit('loading data from json, the types of atoms does not match!')

        self.put_data(data)
        self.update_config_json()

        if guess_geom:
            self.update_system(np.array(data['coord']))

    def update_config_json(self):
        # Update the configuration from JSON
        config = self.config
        if config['guess']['type'] != 'json':
            return
        if (config['input']['basis'] == config['json']['basis'] and
                config['scf']['init_library'] == config['json']['library']):
            return
        self.config['json']['do_init'] = 'yes'
        self.config['scf']['init_scf'] = self.config['json']['scf_type']
        self.config['scf']['init_basis'] = self.config['json']['basis']
        self.config['scf']['init_library'] = self.config['json']['library']

    def put_data(self, data):
        # convert list to data
        for key in self.tag:
            try:
                self.data[key] = np.array(data[key])

            except KeyError:
                continue
        for key in self.config_tag.keys():
            for item in self.config_tag[key]:
                try:
                    self.config[key][item] = data[key][item]
                except KeyError:
                    print(f"Warning: Key {key} not found in data")
                except Exception as e:
                    print(f"Error: {e}")

    def read_freqs(self):
        jsonfile = self.log.replace('.log', '.hess.json')

        if not os.path.exists(jsonfile):
            exit(f'hess file {jsonfile} does not exist')

        with open(jsonfile, 'r') as indata:
            data = json.load(indata)

        energy = data['energy']
        hessian = data['hessian']
        freqs = data['freqs']
        modes = data['modes']
        inertia = data['inertia']

        return energy, hessian, freqs, modes, inertia

    def check_ref(self):
        # compare test data with ref data
        runtype = self.config['input']['runtype']
        ref_file = self.input_file.replace('.inp', '.json')
        runtime_data = self.get_data()
        runtime_data.update(self.get_results())
        skip_keys = [
            'OQP::VEC_MO_A', 'OQP::VEC_MO_B',
            'OQP::td_abxc', 'OQP::td_bvec_mo', 'OQP::td_mrsf_density',
            'OQP::td_states_overlap', 'OQP::state_sign', 'OQP::td_states_phase',
            'OQP::dc_matrix', 'OQP::nac_matrix', 'OQP::DM_A', 'OQP::DM_B', 'OQP::DM_B', 'E_MO_A', 'OQP::Hcore',
            'OQP::SM', 'OQP::TM', 'OQP::FOCK_A', 'OQP::FOCK_B', 'OQP::E_MO_A', 'OQP::E_MO_B', 'OQP::WAO', 'json'
        ]

        if runtype in ['energy']:
            skip_keys.append('grad')
            skip_keys.append('hess')

        if runtype in ['grad', 'optimize', 'meci', 'mep']:
            skip_keys.append('hess')

        if runtype in ['hess', 'nacme', 'nac']:
            skip_keys.append('grad')

        message = ''
        total_diff = 0

        if os.path.exists(ref_file):
            message += f'   PyOQP reference data {ref_file}\n'

            with open(ref_file, 'r') as indata:
                ref_data = json.load(indata)

            for key, value in ref_data.items():
                if key in skip_keys:
                    continue
                flag, diff = compare_data(runtime_data[key], value)
                total_diff += diff
                message += f'   PyOQP checking {key:<20} ... {flag} ({diff:.8f})\n'
        else:
            message += '   PyOQP reference data is not found (skip and save data)\n'
            self.save_data()

        return message, total_diff


def compare_data(data_1, data_2):
    """
    Compute the numerical differences between two arrays
    """
    diff = np.sum(np.abs(np.array(data_1) - np.array(data_2)))
    if np.round(diff, 4) > 0:
        return 'failed', diff

    return 'passed', diff


def get_coord(xyz, nat):
    """Get coordinate"""
    return np.frombuffer(oqp.ffi.buffer(xyz, 3 * nat * oqp.ffi.sizeof("double")),
                         dtype=np.double).reshape((nat, 3))


def string_config(config):
    # convert dict value to strings
    str_config = {}
    for section in config.keys():
        str_config[section] = {}
        for option, value in config[section].items():
            if isinstance(value, list) or isinstance(value, tuple):
                value = list2string(value)
            else:
                value = str(value)

            str_config[section][option] = value

    return str_config


def list2string(in_list):
    # convert list to str
    # [1, 2, 3, 4] -> '1, 2, 3, 4'
    # [[1 2], [3, 4]] -> '1 2, 3 4'
    # do not support three-layer list

    if len(in_list) == 0:
        return ''

    str_list = []
    for item in in_list:
        if isinstance(item, list):
            item = ' '.join([str(x) for x in item])
            if isinstance(item[0], list):
                raise ValueError('do not support three-layer list %s' % in_list)
        else:
            item = str(item)

        str_list.append(item)

    str_list = ','.join(str_list)

    return str_list
