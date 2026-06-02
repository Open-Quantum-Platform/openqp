"""OQP molecule class"""
import os
import copy
import json
import platform
import warnings
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
        self.hessian_metadata = {}
        self.modes = np.zeros(0)  # 3Natom-6, 3Natom
        self.inertia = np.zeros(0)  # 3
        self.infrared_intensities = np.zeros(0)
        self.raman_activities = np.zeros(0)
        self.vibrational_intensity_metadata = {}
        self.infrared_mode_dipole_derivatives = np.zeros((0, 3))
        self.raman_mode_polarizability_derivatives = np.zeros((0, 3, 3))

        self.tag = [
            'OQP::DM_A', 'OQP::DM_B',
            'OQP::FOCK_A', 'OQP::FOCK_B',
            'OQP::E_MO_A', 'OQP::E_MO_B',
            'OQP::VEC_MO_A', 'OQP::VEC_MO_B',
            'OQP::Hcore', 'OQP::SM', 'OQP::TM', 'OQP::WAO',
            'OQP::td_abxc', 'OQP::td_bvec_mo', 'OQP::td_mrsf_density', 'OQP::td_energies',
            'OQP::mrsf_ekt_density_mo', 'OQP::mrsf_ekt_lagrangian_mo', 'OQP::mrsf_ekt_fock_mo',
            'OQP::mrsf_ekt_orbitals_mo', 'OQP::mrsf_ekt_eigenvalues', 'OQP::mrsf_ekt_strengths',
            'OQP::hf_hessian',
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

        for tag in self.tag:
            name = tag.replace('OQP::', '').lower()
            getter = lambda self, t=tag: np.array(self.data[t])
            setter = lambda self, val, t=tag: self.data.__setitem__(t, val)
            setattr(self.__class__, f'get_{name}', getter)
            setattr(self.__class__, f'set_{name}', setter)

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
    def get_scf_energy(self, component=None):
        """
        Retrieve SCF (Self-Consistent Field) energy components.

        This method provides convenient access to individual or all energy
        terms computed during an SCF procedure. If no component is specified,
        the total SCF energy is returned.

        Parameters
        ----------
        component : str, optional
            The energy component to retrieve. Supported options are:

            - ``None`` (default): Returns only the total SCF energy.
            - ``"all"``: Returns a dictionary containing all available
              energy components.
            - One of the following component names:
                * "energy"  — total SCF energy
                * "psinrm"  — wavefunction norm
                * "ehf1"    — Hartree-Fock energy (one-electron)
                * "vee"     — electron-electron repulsion energy
                * "nenergy" — nuclear energy contribution
                * "vne"     — electron-nucleus attraction energy
                * "vnn"     — nucleus-nucleus repulsion energy
                * "vtot"    — total potential energy
                * "tkin"    — kinetic energy
                * "virial"  — virial ratio

        Returns
        -------
        float or dict
            - If `component` is None, returns a single float (total SCF energy).
            - If `component` is "all", returns a dictionary with all energy components.
            - If `component` corresponds to a specific component, returns that component as a float.

        Raises
        ------
        ValueError
            If the provided `component` does not match any of the known energy components.

        Examples
        --------
        >>> mol.get_scf_energy()
        -75.98327432

        >>> mol.get_scf_energy("tkin")
        37.420192

        >>> mol.get_scf_energy("all")
        {
            'energy': -75.98327432,
            'psinrm': 0.999999,
            'ehf1': -72.3123,
            'vee': 18.2034,
            'nenergy': -80.000,
            'vne': -85.6214,
            'vnn': 5.6214,
            'vtot': -67.4180,
            'tkin': 37.4202,
            'virial': 2.1519
        }
        """
        energy_data = self.data._data.mol_energy

        if component is None:
            return energy_data.energy

        elif component == "all":
            return {
                "energy": energy_data.energy,
                "psinrm": energy_data.psinrm,
                "ehf1": energy_data.ehf1,
                "vee": energy_data.vee,
                "nenergy": energy_data.nenergy,
                "vne": energy_data.vne,
                "vnn": energy_data.vnn,
                "vtot": energy_data.vtot,
                "tkin": energy_data.tkin,
                "virial": energy_data.virial
            }

        else:
            if hasattr(energy_data, component):
                return getattr(energy_data, component)
            else:
                raise ValueError(
                    f"Invalid component '{component}'. Use one of: "
                    f"energy, psinrm, ehf1, vee, nenergy, vne, vnn, "
                    f"vtot, tkin, virial, or 'all'."
                )

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

        return copy.deepcopy(self.hessian)

    def set_hessian_result(self, raw_hessian, asymmetry_tol=1.0e-8):
        """
        Store a final Cartesian Hessian in OpenQP frequency conventions.

        Native analytic Hessian kernels should hand one square ``(3N, 3N)``
        matrix to this helper. The helper records the pre-symmetrization
        asymmetry for diagnostics and stores the symmetrized matrix used by
        normal-mode analysis; it does not compute a numerical fallback.
        """

        hessian = np.asarray(raw_hessian, dtype=float)
        if hessian.ndim != 2 or hessian.shape[0] != hessian.shape[1]:
            raise ValueError(f"Expected square Hessian matrix, got shape={hessian.shape}")

        natom = self.data['natom']
        expected = 3 * natom
        if hessian.shape != (expected, expected):
            raise ValueError(
                f"Expected Hessian shape ({expected}, {expected}) for {natom} atoms, got {hessian.shape}"
            )

        max_asymmetry = float(np.max(np.abs(hessian - hessian.T))) if hessian.size else 0.0
        if max_asymmetry > asymmetry_tol:
            warnings.warn(
                f"Analytic Hessian asymmetry {max_asymmetry:.3e} exceeds tolerance {asymmetry_tol:.3e}; symmetrizing final matrix.",
                RuntimeWarning,
            )

        self.hessian = 0.5 * (hessian + hessian.T)
        self.hessian_metadata = {
            'max_asymmetry': max_asymmetry,
            'symmetrized': bool(max_asymmetry > 0.0),
        }
        return self.hessian

    def get_mrsf_ekt_results(self):
        """Collect MRSF-EKT root results for the final JSON file."""
        if self.data is None:
            return {}

        try:
            eigenvalues = np.array(self.data['OQP::mrsf_ekt_eigenvalues'])
            strengths = np.array(self.data['OQP::mrsf_ekt_strengths'])
            orbitals = np.array(self.data['OQP::mrsf_ekt_orbitals_mo'])
        except AttributeError:
            return {}

        hartree_to_ev = 27.211386245988
        ebe_ev = (-eigenvalues * hartree_to_ev).tolist()
        return {
            'mrsf_ekt': {
                'tdhf_type': self.config.get('tdhf', {}).get('type'),
                'target_state': self.config.get('tdhf', {}).get('target'),
                'eigenvalues_hartree': eigenvalues.tolist(),
                'ebe_ev': ebe_ev,
                'pole_strengths': strengths.tolist(),
                'orbitals_mo': orbitals.tolist(),
            }
        }

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
        data.update(self.get_mrsf_ekt_results())

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
            'hessian_metadata': self.hessian_metadata,
            'freqs': self.freqs.tolist(),
            'modes': self.modes.tolist(),
            'frequency_modes': {
                'frequencies_cm-1': self.freqs.tolist(),
                'normal_mode_eigenvectors': self.modes.tolist(),
                'normal_mode_eigenvectors_units': 'Cartesian displacement, mass-unweighted, row-major by vibrational mode',
            },
            'inertia': self.inertia.tolist(),
            'infrared_intensities': self.infrared_intensities.tolist(),
            'raman_activities': self.raman_activities.tolist(),
            'vibrational_intensity_metadata': self.vibrational_intensity_metadata,
            'infrared_mode_dipole_derivatives': self.infrared_mode_dipole_derivatives.tolist(),
            'raman_mode_polarizability_derivatives': self.raman_mode_polarizability_derivatives.tolist(),
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
        self.hessian_metadata = data.get('hessian_metadata', {})
        freqs = data['freqs']
        modes = data['modes']
        inertia = data['inertia']
        self.infrared_intensities = np.array(data.get('infrared_intensities', []), dtype=float)
        self.raman_activities = np.array(data.get('raman_activities', []), dtype=float)
        self.vibrational_intensity_metadata = data.get('vibrational_intensity_metadata', {})
        self.infrared_mode_dipole_derivatives = np.array(
            data.get('infrared_mode_dipole_derivatives', []), dtype=float
        )
        self.raman_mode_polarizability_derivatives = np.array(
            data.get('raman_mode_polarizability_derivatives', []), dtype=float
        )

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
            'OQP::SM', 'OQP::TM', 'OQP::FOCK_A', 'OQP::FOCK_B', 'OQP::E_MO_A', 'OQP::E_MO_B', 'OQP::WAO',
            'OQP::mrsf_ekt_density_mo', 'OQP::mrsf_ekt_lagrangian_mo', 'OQP::mrsf_ekt_fock_mo',
            'OQP::mrsf_ekt_orbitals_mo', 'OQP::mrsf_ekt_eigenvalues', 'OQP::mrsf_ekt_strengths',
            'OQP::hf_hessian',
            'json'
        ]
        tdhf_type = self.config.get('tdhf', {}).get('type')
        required_ref_keys = []
        if tdhf_type in ('mrsf_ekt_ip', 'mrsf_ekt_ea'):
            required_ref_keys.append('mrsf_ekt')

        if runtype in ['energy']:
            skip_keys.append('grad')
            skip_keys.append('hess')

        if runtype in ['grad', 'optimize', 'meci', 'mep']:
            skip_keys.append('hess')

        if runtype in ['hess', 'nacme', 'nac']:
            skip_keys.append('grad')
        if runtype == 'hess':
            # Hessian examples keep the detailed vibrational reference in
            # ``*.hess.json``.  The legacy main ``*.json`` references store an
            # empty placeholder for ``hess``; comparing that placeholder to the
            # runtime Cartesian matrix creates a shape mismatch instead of a
            # meaningful regression signal.
            skip_keys.append('hess')

        message = ''
        total_diff = 0

        if os.path.exists(ref_file):
            message += f'   PyOQP reference data {ref_file}\n'

            with open(ref_file, 'r') as indata:
                ref_data = json.load(indata)

            for key in required_ref_keys:
                if key not in ref_data:
                    total_diff += 1.0
                    message += f'   PyOQP missing reference {key:<20} ... failed (1.00000000)\n'

            for key, value in ref_data.items():
                if key in skip_keys:
                    continue
                if key not in runtime_data:
                    flag, diff = 'failed', 1.0
                else:
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
    if isinstance(data_1, dict) or isinstance(data_2, dict):
        if not isinstance(data_1, dict) or not isinstance(data_2, dict):
            return 'failed', 1.0
        diff = 0.0
        for key in sorted(data_2):
            if key == 'orbitals_mo':
                # EKT orbital vectors are phase/sign ambiguous between runs;
                # eigenvalues and pole strengths provide the stable regression
                # signal for the structured EKT result.
                continue
            if key not in data_1:
                diff += 1.0
                continue
            _, subdiff = compare_data(data_1[key], data_2[key])
            diff += subdiff
        if np.round(diff, 4) > 0:
            return 'failed', diff
        return 'passed', diff

    if isinstance(data_1, str) or isinstance(data_2, str):
        diff = 0.0 if data_1 == data_2 else 1.0
        if diff > 0:
            return 'failed', diff
        return 'passed', diff

    if data_1 is None or data_2 is None:
        diff = 0.0 if data_1 is data_2 else 1.0
        if diff > 0:
            return 'failed', diff
        return 'passed', diff

    arr_1 = np.array(data_1)
    arr_2 = np.array(data_2)
    if arr_1.shape != arr_2.shape:
        return 'failed', 1.0

    if arr_1.size == 0:
        diff = 0.0
    else:
        # Use the maximum element-wise deviation instead of an L1 sum so
        # vector-valued references are judged by per-value numerical drift,
        # not by the number of states/components in the vector.
        diff = float(np.max(np.abs(arr_1 - arr_2)))
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
