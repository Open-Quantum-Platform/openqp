"""OQP single point class"""
import os
import oqp
import copy
import time
import shutil
import platform
import subprocess
import multiprocessing
import numpy as np
from numpy import linalg as la
from oqp.molecule import Molecule
from oqp.utils.mpi_utils import MPIManager, MPIPool

try:
    from dftd4.interface import DampingParam, DispersionModel

    dftd_installed = 'dftd4'
except ModuleNotFoundError:
    from oqp.utils.matrix import DampingParam, DispersionModel

    dftd_installed = 'not installed'

from oqp.library.frequency import normal_mode, thermal_analysis
from oqp.utils.file_utils import dump_log, dump_data, write_config, write_xyz


class Calculator:
    """
    OQP calculator base class

    """

    def __init__(self, mol):
        self.mol = mol
        self.save_mol = mol.config['guess']['save_mol']
        self.export = mol.config['properties']['export']
        self.export_title = mol.config['properties']['title']
        self.exception = mol.config['tests']['exception']
        self.mpi_manager = MPIManager()


class LastStep(Calculator):
    """
    OQP last step single point calculation class

    """

    def __init__(self, mol, param=None):
        super().__init__(mol)

        self.functional = mol.config['input']['functional']
        if len(self.functional) == 0:
            self.functional = 'hf'

        self.do_d4 = mol.config['input']['d4']
        self.res = None
        self.set_param(param)
        self.natom = 0

        dump_log(
            mol,
            title='PyOQP: Dispersion Correction',
            section='dftd',
            info={'type': dftd_installed, 'd4': self.do_d4}
        )

    def get_dispersion(self, mol, grad_list):
        if grad_list:
            do_grad = True
        else:
            do_grad = False

        atoms = mol.get_atoms()
        coordinates = mol.get_system().reshape((-1, 3))
        natom = len(atoms)
        if self.do_d4:
            model = DispersionModel(atoms, coordinates)
            res = model.get_dispersion(DampingParam(method=self.functional),
                                       grad=do_grad)

            energy = res['energy']
            if do_grad:
                grad = res['gradient']
            else:
                grad = np.zeros((natom, 3))
        else:
            energy = 0.0
            grad = np.zeros((natom, 3))

        return energy, grad

    def set_param(self, param):
        # TODO pass user-defined parameters to dftd4
        pass

    def compute(self, mol, grad_list=None):
        # do dftd4
        energy, grad = self.get_dispersion(mol, grad_list)
        energies = self.final_energy(energy)

        if not grad_list:
            return energies

        grads = self.final_grad(grad, grad_list)
        return energies, grads

    def final_energy(self, d4_energy):
        # add dftd4 energy
        el = self.mol.energies
        self.mol.energies = [x + d4_energy for x in el]

        dump_log(
            self.mol,
            title='PyOQP: Final Energy',
            section='energy',
            info={'el': el, 'd4': d4_energy}
        )

        # export data
        if self.export:
            dump_data(self.mol, (self.mol.energies, self.export_title), title='ENERGY', fpath=self.mol.log_path)

        # save mol
        if self.save_mol:
            self.mol.save_data()

        return self.mol.energies

    def final_grad(self, d4_grad, grad_list):
        # add dftd4 grad
        el = self.mol.grads
        self.mol.grads = [x + d4_grad for x in el]

        dump_log(
            self.mol,
            title='PyOQP: Final Gradient',
            section='grad',
            info={'el': el, 'd4': d4_grad, 'grad_list': grad_list}
        )

        # export data
        if self.export:
            dump_data(self.mol, (self.mol.grads, self.export_title, grad_list), title='GRADIENT',
                      fpath=self.mol.log_path)

        # save mol
        if self.save_mol:
            self.mol.save_data()

        return self.mol.grads


class SinglePoint(Calculator):
    """
    OQP single point calculation class

    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.method = mol.config['input']['method']
        self.functional = mol.config['input']['functional']
        self.basis = mol.config['input']['basis']
        self.library = mol.config['input']['library']
        self.scf_type = mol.config['scf']['type']
        self.scf_maxit = mol.config['scf']['maxit']
        self.forced_attempt = mol.config['scf']['forced_attempt']
        self.alternative_scf = mol.config["scf"]["alternative_scf"]
        self.scf_mult = mol.config['scf']['multiplicity']
        self.init_scf = mol.config['scf']['init_scf']
        self.init_it = mol.config['scf']['init_it']
        self.init_basis = mol.config['scf']['init_basis']
        self.init_library = mol.config['scf']['init_library']
        self.init_conv = mol.config['scf']['init_conv']
        self.conv = mol.config['scf']['conv']
        self.save_molden = mol.config['scf']['save_molden']
        self.td = mol.config['tdhf']['type']
        self.nstate = mol.config['tdhf']['nstate']
        self.energy_func = {
            'hf': oqp.hf_energy,
            'rpa': oqp.tdhf_energy,
            'tda': oqp.tdhf_energy,
            'sf': oqp.tdhf_sf_energy,
            'mrsf': oqp.tdhf_mrsf_energy,
        }

        # initialize state sign
        self.mol.data["OQP::state_sign"] = np.ones(self.nstate)

    def _prep_guess(self):
        oqp.library.set_basis(self.mol)
        oqp.library.ints_1e(self.mol)
        oqp.library.guess(self.mol)

    def _project_basis(self):
        oqp.library.project_basis(self.mol)

    def _init_convergence(self):
        init_calc = self.energy_func['hf']
        target_basis = self.basis
        target_library = self.library
        if self.init_basis == 'none':
            init_basis = target_basis
            init_library = target_library
        else:
            init_basis = self.init_basis
            init_library = self.init_library

        init_converger = self.mol.config['scf']['init_converger']
        target_converger = self.mol.config['scf']['soscf_type']
        self.mol.data.set_scf_soscf_type(init_converger)
        self.mol.data.set_scf_conv(self.init_conv)

        if init_basis:
            self.mol.config['input']['basis'] = init_basis
            self.mol.config['input']['library'] = init_library

        self.mol.data.set_scf_maxit(self.init_it)

        if self.init_scf == 'rhf':
            self.mol.config['input']['functional'] = ''
            self.mol.data.set_dft_functional('')
            self.mol.data.set_scf_type('rhf')
            self.mol.data.set_mol_multiplicity(1)

        elif self.init_scf == 'uhf':
            self.mol.config['input']['functional'] = ''
            self.mol.data.set_dft_functional('')
            self.mol.data.set_scf_type('uhf')
            self.mol.data.set_mol_multiplicity(3)

        elif self.init_scf == 'rohf':
            self.mol.config['input']['functional'] = ''
            self.mol.data.set_dft_functional('')
            self.mol.data.set_scf_type('rohf')
            self.mol.data.set_mol_multiplicity(3)

        elif self.init_scf == 'rks':
            self.mol.data.set_scf_type('rhf')
            self.mol.data.set_mol_multiplicity(1)

        elif self.init_scf == 'uks':
            self.mol.data.set_scf_type('uhf')
            self.mol.data.set_mol_multiplicity(3)

        elif self.init_scf == 'roks':
            self.mol.data.set_scf_type('rohf')
            self.mol.data.set_mol_multiplicity(3)
        else:
            raise ValueError(f'Unknown initial scf method {self.init_scf}')

        dump_log(self.mol, title='PyOQP: Initial SCF steps', section='scf')
        self._prep_guess()
        if self.mol.config['guess']['type'] != 'json':
            init_calc(self.mol)
        # save initially converge orbitals
        if self.save_molden:
            guess_file = self.pack_molden_name('init', self.init_scf, self.mol.config['input']['functional'])
            self.mol.write_molden(guess_file)

        if init_basis:
            dump_log(self.mol, title='OQP: Applying Basis Sets for Overlap calculation', section='scf')
            self.mol.data.set_scf_active_basis(1)
            oqp.library.set_basis(self.mol)
            self.mol.data.set_scf_active_basis(0)
            self.mol.config['input']['basis'] = target_basis
            self.mol.config['input']['library'] = target_library
            oqp.library.set_basis(self.mol)

        # set parameters back to normal scf
        self.mol.config['input']['basis'] = target_basis
        self.mol.config['input']['functional'] = self.functional
        self.mol.data.set_scf_soscf_type(target_converger)
        self.mol.data.set_dft_functional(self.functional)
        self.mol.data.set_scf_type(self.scf_type)
        self.mol.data.set_scf_maxit(self.scf_maxit)
        self.mol.data.set_scf_conv(self.conv)
        self.mol.data.set_mol_multiplicity(self.scf_mult)
        self._project_basis()
        #        oqp.library.update_guess(self.mol)

        dump_log(self.mol, title='PyOQP: Initial SCF steps done, switching back to normal SCF', section='scf')

    def pack_molden_name(self, cal_type, scf_type, functional):
        """Add information to molden file name"""
        if len(functional) == 0:
            functional = 'hf'
        basis = self.basis.replace('*', 's')
        if self.mol.idx != 1:
            guess_file = self.mol.log.replace('.log', '_%s_%s_%s_%s_%s.molden' % (
                self.mol.idx, cal_type, scf_type, functional, basis))
        else:
            guess_file = self.mol.log.replace('.log', '_%s_%s_%s_%s.molden' % (
                cal_type, scf_type, functional, basis))

        return guess_file

    @staticmethod
    def molden_unique_name(base_filename):
        """Check if the base file exists"""
        print('kk: check base name:', base_filename)
        if not os.path.exists(base_filename):
            # The file does not exist
            final_filename = base_filename
        else:
            # The file exists -- generate a new filename with '_i'
            i = 1
            while True:
                new_filename = base_filename.replace('.molden', f'_{i}.molden')
                # Check if this new filename exists
                if not os.path.exists(new_filename):
                    # If it doesn't exist, use this filename
                    final_filename = new_filename
                    break
                i += 1
        return final_filename

    # IXCORE for XAS (X-ray absorption spectroscopy)
    # Fock matrix is in AO here, so we need to shift it in Fortran after transform it into MO
    def ixcore_shift(self):
        from oqp import ffi
        ixcore = self.mol.config["tdhf"]["ixcore"]
        if ixcore == "-1":  # if default
            return
        ixcore_array = np.array(ixcore.split(','), dtype=np.int32)
        # shift MO energies 
        noccB = self.mol.data['nelec_B']
        tmp = self.mol.data["OQP::E_MO_A"]
        for i in range(noccB + 1):  # up to HOMO-1
            if i not in ixcore_array:
                tmp[i - 1] = -100000  # shift the MO energy down

    def energy(self, do_init_scf=True):
        # check method
        if self.method not in ['hf', 'tdhf']:
            raise ValueError(f'Unknown method type {self.method}')

        # compute reference
        ref_energy = self.reference(do_init_scf=do_init_scf)

        # ixcore
        self.ixcore_shift()

        # compute excitations
        if self.method == 'tdhf':
            energies = self.excitation(ref_energy)
        else:
            energies = ref_energy

        return energies

    def swapmo(self):
        # swap MO energy and AO coefficient depending on user's request
        swapmo = self.mol.config["guess"]["swapmo"]
        if swapmo:  # if not default (empty)
            swapmo_array = [int(x.strip()) for x in swapmo.split(',')]

            # Initial MO energy and coefficient
            og_val = self.mol.data["OQP::E_MO_A"]
            og_vec = self.mol.data["OQP::VEC_MO_A"]
            # It only takes pairs. If it is not pair, it will be ignored.
            for i, j in zip(swapmo_array[::2], swapmo_array[1::2]):
                og_val[[i - 1, j - 1]] = og_val[[j - 1, i - 1]]
                og_vec[[i - 1, j - 1]] = og_vec[[j - 1, i - 1]]

    def reference(self, do_init_scf=True):
        dump_log(self.mol, title='PyOQP: Entering Electronic Energy Calculation', section='input')

        if self.init_scf != 'no' and do_init_scf:
            # do initial scf iteration to help convergence
            self._init_convergence()
        else:
            self._prep_guess()

        self.swapmo()

        scf_flag = False
        itr = 0
        energy = 0
        for itr in range(self.forced_attempt):
            self.scf()
            energy = [self.mol.mol_energy.energy]

            # check convergence
            scf_flag = self.mol.mol_energy.SCF_converged

            if scf_flag:
                dump_log(self.mol, title='PyOQP: SCF energy is converged after %s attempts' % (itr + 1), section='')
                break
            else:
                dump_log(self.mol, title='PyOQP: SCF energy is not converged after %s attempts' % (itr + 1), section='')

        if not scf_flag:
            dump_log(self.mol, title='PyOQP: SCF energy is not converged', section='end')

            if self.alternative_scf:
                dump_log(self.mol, title='PyOQP: Enable the SOSCF flag in SCF to improve convergence.', section='input')
                self.mol.data.set_scf_soscf_type(1)
                self.scf()
                energy = [self.mol.mol_energy.energy]
                scf_flag = self.mol.mol_energy.SCF_converged
                if scf_flag:
                    dump_log(self.mol, title='PyOQP: SCF energy is converged after %s attempts' % (itr + 1), section='')

            if self.exception is True:
                raise SCFnotConverged()
            else:
                exit()

        if self.save_molden:
            guess_file = self.pack_molden_name('scf', self.scf_type, self.functional)
            self.mol.write_molden(guess_file)

        self.mol.energies = energy

        return energy

    def excitation(self, ref_energy):
        self.tddft()
        energies = ref_energy + [ex + ref_energy[0] for ex in self.mol.data['OQP::td_energies']]

        # check convergence
        td_flag = self.mol.mol_energy.Davidson_converged

        if self.method != 'hf' and not td_flag:
            dump_log(self.mol, title='PyOQP: TD energy is not converged', section='end')

            if self.exception is True:
                raise TDnotConverged()
            else:
                exit()

        self.mol.energies = energies

        return energies

    def scf(self):
        # do SCF
        dump_log(self.mol, title='PyOQP: Normal SCF steps', section='scf')
        self.energy_func['hf'](self.mol)

    def tddft(self):
        # check td type
        if self.td not in ['rpa', 'tda', 'sf', 'mrsf']:
            raise ValueError(f'Unknown tdhf type {self.td}')

        # do TDDFT
        dump_log(self.mol, title='PyOQP: TDDFT steps', section='tdhf')
        self.energy_func[self.td](self.mol)


class Gradient(Calculator):
    """
    OQP gradient calculation class

    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.method = mol.config["input"]["method"]
        self.td = mol.config["tdhf"]["type"]
        self.grads = mol.config["properties"]["grad"]
        self.natom = mol.data["natom"]
        self.nstate = mol.config['tdhf']['nstate']

        self.zvec_func = {
            'rpa': oqp.tdhf_z_vector,
            'tda': oqp.tdhf_z_vector,
            'sf': oqp.tdhf_sf_z_vector,
            'mrsf': oqp.tdhf_mrsf_z_vector,
        }

        self.grad_func = {
            'hf': oqp.hf_gradient,
            'rpa': oqp.tdhf_gradient,
            'tda': oqp.tdhf_gradient,
            'sf': oqp.tdhf_sf_gradient,
            'mrsf': oqp.tdhf_mrsf_gradient,
        }

    def gradient(self):
        # check method
        if self.method not in ['hf', 'tdhf']:
            raise ValueError(f'Unknown method type {self.method}')

        dump_log(self.mol, title='PyOQP: Entering Gradient Calculation')

        # compute gradients
        grads = []
        if self.method == 'hf':
            grads = self.scf_grad()
        elif self.method == 'tdhf':
            grads = self.tddft_grad()

        self.mol.grads = grads

        return grads

    def scf_grad(self):
        dump_log(self.mol, title='PyOQP: Gradient of Root 0')
        self.grad_func['hf'](self.mol)
        grad = self.mol.get_grad()
        grads = np.array([grad.copy()]).reshape((1, self.natom, 3))

        return grads

    def tddft_grad(self):
        if self.td not in ['rpa', 'tda', 'sf', 'mrsf']:
            raise ValueError(f'Unknown tdhf type {self.td}')

        if self.nstate < max(self.grads):
            raise ValueError(f'Gradient requested state {max(self.grads)} > the highest computed state {self.nstate}')

        if self.nstate == max(self.grads):
            print(f'\nGradient requested state {max(self.grads)} == the highest computed state {self.nstate}')
            print(f'It is recommended to compute {self.nstate + 1} states to avoid missing degenerate states\n')

        grads = np.zeros((self.nstate + 1, self.natom, 3))
        for i in self.grads:
            dump_log(self.mol, title='PyOQP: Gradient of Root %s' % i)
            self.mol.data.set_tdhf_target(i)
            self.zvec_func[self.td](self.mol)

            # check convergence
            z_flag = self.mol.mol_energy.Z_Vector_converged

            if not z_flag:
                dump_log(self.mol, title='PyOQP: TD Z-vector is not converged', section='end')

                if self.exception is True:
                    raise ZVnotConverged()
                else:
                    exit()

            self.grad_func[self.td](self.mol)
            grad = self.mol.get_grad().reshape((self.natom, 3))
            grads[i] = grad.copy()

        return grads


"""
Complete Dyson class
"""
class Dyson(Calculator):
    """
    OQP Dyson orbital calculation class
    Computes Dyson orbitals using Extended Koopmans' Theorem (EKT)
    Requires gradient calculation first to obtain relaxed density and Lagrangian
    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.method = mol.config["input"]["method"]
        self.td = mol.config["tdhf"]["type"]
        self.dyson_config = mol.config.get('dyson', {})
        
        # Dyson calculation parameters
        self.compute_ip = self.dyson_config.get('ip', True)
        self.compute_ea = self.dyson_config.get('ea', False)
        self.target_state = self.dyson_config.get('target_state', 0)
        self.threshold = self.dyson_config.get('threshold', 1.0e-4)
        self.deflation_tol = self.dyson_config.get('deflation_tol', 1.0e-8)
        self.max_orbitals = self.dyson_config.get('max_orbitals', 0)
        self.print_level = self.dyson_config.get('print_level', 1)
        self.max_print = self.dyson_config.get('max_print', 10)
        self.export_dyson = self.dyson_config.get('export', False)
        self.export_format = self.dyson_config.get('export_format', 'text')
        self.export_name = self.dyson_config.get('export_name', 'dyson_results')
        self.save_cube = self.dyson_config.get('save_cube', False)
        self.spectrum = self.dyson_config.get('spectrum', False)
        self.orbital_analysis = self.dyson_config.get('orbital_analysis', False)
        
        # If target_state is 0, use TDHF target
        if self.target_state == 0:
            self.target_state = mol.config.get('tdhf', {}).get('target', 1)
        
        # Function mappings for different TD methods
        self.dyson_func = {
            'rpa': oqp.dyson_orbital_calculation,
            'tda': oqp.dyson_orbital_calculation,
            'sf': oqp.dyson_mrsf,
            'mrsf': oqp.dyson_mrsf,
        }

    def dyson(self):
        """
        Main entry point for Dyson orbital calculation.
        Workflow:
        1. Validate input
        2. Set parameters for Fortran
        3. Call appropriate Dyson calculation
        4. Retrieve and process results
        5. Export/save if requested
        """
        # Check method
        if self.method != 'tdhf':
            raise ValueError(f'Dyson calculation requires method=tdhf, but found {self.method}')
        
        # Check TD type
        if self.td not in ['rpa', 'tda', 'sf', 'mrsf']:
            raise ValueError(f'Unknown tdhf type {self.td} for Dyson calculation')
        

        # Check that Z-vector was computed
        if not self.mol.mol_energy.Z_Vector_converged:
           dump_log(self.mol, title='PyOQP: Error - Z-vector not computed', section='error')
           raise ValueError('Z-vector must be computed before Dyson calculation')

        dump_log(self.mol, title='PyOQP: Entering Dyson Orbital Calculation')
        dump_log(self.mol, title='', section='dyson_setup', info={
            'method': self.td,
            'target_state': self.target_state,
            'compute_ip': self.compute_ip,
            'compute_ea': self.compute_ea,
            'threshold': self.threshold,
            'deflation_tol': self.deflation_tol,
            'max_orbitals': self.max_orbitals,
            'print_level': self.print_level
        })
        
        # Set Dyson parameters in the data structure for Fortran to access
        self.set_dyson_parameters()
        
        # Call appropriate Dyson calculation based on TD method
        dump_log(self.mol, title=f'PyOQP: Calling Fortran Dyson calculation ({self.td})')
        self.dyson_func[self.td](self.mol)
        
        # Process and retrieve results from tagarray including orbitals
        results = self.retrieve_results_with_orbitals()
        
        # Store results in mol for access
        self.mol.dyson_results = results
        
        # Export if requested
        if self.export_dyson and results and results['n_orbitals'] > 0:
            self.export_results(results)
        
        # Save cube files if requested
        if self.save_cube and results and results['n_orbitals'] > 0:
            self.save_cube_files(results)
        
        # Generate spectrum if requested
        if self.spectrum and results and results['n_orbitals'] > 0:
            self.generate_spectrum(results)
        
        # Save mol data if requested
        if self.save_mol:
            self.mol.save_data()
        
        return results
    
    def set_dyson_parameters(self):
        """
        Set Dyson parameters in the data structure for Fortran access.
        Uses the tddft structure that's already bound in C.
        """
        # These parameters map to tddft_parameters type in types.F90
        self.mol.data._data.tddft.dyson_compute_ip = self.compute_ip
        self.mol.data._data.tddft.dyson_compute_ea = self.compute_ea
        self.mol.data._data.tddft.dyson_target_state = self.target_state
        self.mol.data._data.tddft.dyson_pole_threshold = self.threshold
        self.mol.data._data.tddft.dyson_deflation_tol = self.deflation_tol
        self.mol.data._data.tddft.dyson_max_orbitals = self.max_orbitals
        self.mol.data._data.tddft.dyson_print_level = self.print_level
        self.mol.data._data.tddft.dyson_max_print = self.max_print
        self.mol.data._data.tddft.dyson_export = self.export_dyson
        self.mol.data._data.tddft.dyson_orbital_analysis = self.orbital_analysis
        
        # Set export format as integer
        format_map = {'text': 1, 'csv': 2, 'json': 3, 'numpy': 4}
        self.mol.data._data.tddft.dyson_export_format = format_map.get(self.export_format, 1)
        
        # Enable Dyson calculation globally
        self.mol.data._data.control.dyson_enabled = True
    
    def retrieve_results_with_orbitals(self):
        """
        Retrieve complete Dyson calculation results from tagarray including orbitals.
        Results are stored in tagarray by Fortran code.
        """
        try:
            # Get the number of significant orbitals stored
            try:
                n_orbitals = self.mol.data["OQP::DO_Count"]
            except (KeyError, AttributeError):
                n_orbitals = 0
            
            if n_orbitals == 0:
                dump_log(self.mol, title='PyOQP: No significant Dyson orbitals found', section='warning')
                return self.empty_results()
            
            # Access Dyson orbital coefficients (like MO coefficients)
            vec_do = self.mol.data["OQP::VEC_DO"]  # Shape: (nbf, n_orbitals)
            
            # Access Dyson orbital energies (binding energies)
            e_do = self.mol.data["OQP::E_DO"]  # Shape: (n_orbitals,)
            
            # Access pole strengths
            do_strength = self.mol.data["OQP::DO_Strength"]  # Shape: (n_orbitals,)
            
            # Access orbital types (1=IP, -1=EA)
            do_type = self.mol.data["OQP::DO_Type"]  # Shape: (n_orbitals,)
            
            # Convert to numpy arrays if needed
            vec_do = np.array(vec_do) if not isinstance(vec_do, np.ndarray) else vec_do
            e_do = np.array(e_do) if not isinstance(e_do, np.ndarray) else e_do
            do_strength = np.array(do_strength) if not isinstance(do_strength, np.ndarray) else do_strength
            do_type = np.array(do_type) if not isinstance(do_type, np.ndarray) else do_type
            
            # Reshape vec_do if necessary (Fortran column-major to Python)
            if vec_do.ndim == 1 and n_orbitals > 0:
                nbf = self.mol.data["nbf"]
                vec_do = vec_do.reshape((nbf, n_orbitals), order='F')
            
            # Process results
            results = self.process_results_with_orbitals(vec_do, e_do, do_strength, do_type)
            
            # Log summary
            self.print_summary(results)
            
            return results
            
        except (KeyError, AttributeError) as e:
            dump_log(self.mol, title=f'PyOQP: Error retrieving Dyson results: {e}', section='error')
            return self.empty_results()
    
    def process_results_with_orbitals(self, vec_do, e_do, do_strength, do_type):
        """
        Process Dyson orbital data into structured results.
        Separates IPs and EAs and converts to appropriate units.
        """
        hartree_to_ev = 27.211396
        
        # Ensure arrays are numpy arrays
        do_type = np.array(do_type, dtype=int)
        
        # Separate IPs and EAs based on do_type
        ip_mask = do_type == 1
        ea_mask = do_type == -1
        
        results = {
            # Full orbital data
            'dyson_orbitals': vec_do,  # All Dyson orbital coefficients
            'binding_energies': e_do,   # In Hartree
            'pole_strengths': do_strength,
            'orbital_types': do_type,
            'n_orbitals': len(e_do),
            'nbf': vec_do.shape[0] if vec_do.ndim > 1 else 0,
            'target_state': self.target_state,
            
            # IP-specific data
            'ips': None,
            'ip_poles': None,
            'ip_orbitals': None,
            'n_ip': 0,
            'homo': None,
            
            # EA-specific data
            'eas': None,
            'ea_poles': None,
            'ea_orbitals': None,
            'n_ea': 0,
            'lumo': None,
            
            # Gap
            'gap': None
        }
        
        # Process IPs
        if self.compute_ip and np.any(ip_mask):
            results['ips'] = e_do[ip_mask] * hartree_to_ev  # Convert to eV
            results['ip_poles'] = do_strength[ip_mask]
            results['ip_orbitals'] = vec_do[:, ip_mask] if vec_do.ndim > 1 else None
            results['n_ip'] = int(np.sum(ip_mask))
            if len(results['ips']) > 0:
                results['homo'] = np.min(results['ips'])
        
        # Process EAs
        if self.compute_ea and np.any(ea_mask):
            results['eas'] = -e_do[ea_mask] * hartree_to_ev  # Convert to eV (flip sign)
            results['ea_poles'] = do_strength[ea_mask]
            results['ea_orbitals'] = vec_do[:, ea_mask] if vec_do.ndim > 1 else None
            results['n_ea'] = int(np.sum(ea_mask))
            if len(results['eas']) > 0:
                results['lumo'] = np.min(results['eas'])
        
        # Calculate HOMO-LUMO gap
        if results['homo'] is not None and results['lumo'] is not None:
            results['gap'] = results['homo'] + results['lumo']
        
        return results
    
    def empty_results(self):
        """Return empty results structure"""
        return {
            'dyson_orbitals': np.array([]),
            'binding_energies': np.array([]),
            'pole_strengths': np.array([]),
            'orbital_types': np.array([]),
            'n_orbitals': 0,
            'nbf': 0,
            'target_state': self.target_state,
            'ips': None,
            'ip_poles': None,
            'ip_orbitals': None,
            'n_ip': 0,
            'homo': None,
            'eas': None,
            'ea_poles': None,
            'ea_orbitals': None,
            'n_ea': 0,
            'lumo': None,
            'gap': None
        }
    
    def print_summary(self, results):
        """
        Print summary of Dyson calculation results including orbital info.
        """
        if self.print_level == 0:
            return
        
        dump_log(self.mol, title='PyOQP: Dyson Orbital Results', section='dyson_results')
        
        # Basic summary
        dump_log(self.mol, title='', section='dyson_summary', info={
            'n_orbitals': results['n_orbitals'],
            'nbf': results['nbf'],
            'n_ip': results['n_ip'],
            'n_ea': results['n_ea'],
            'target_state': results['target_state']
        })
        
        # IP details
        if results['n_ip'] > 0 and self.print_level >= 1:
            ip_info = {
                'n_ip': results['n_ip'],
                'homo': f"{results['homo']:.4f} eV" if results['homo'] is not None else None,
                'range': [f"{np.min(results['ips']):.4f}", f"{np.max(results['ips']):.4f}"] if results['n_ip'] > 0 else None
            }
            
            # Add detailed list if verbose
            if self.print_level >= 2:
                n_show = min(self.max_print, results['n_ip'])
                ip_details = []
                for i in range(n_show):
                    ip_details.append(f"  {i+1:3d}: {results['ips'][i]:8.4f} eV (pole: {results['ip_poles'][i]:.6f})")
                ip_info['details'] = ip_details
            
            dump_log(self.mol, title='Ionization Potentials', section='dyson_ip', info=ip_info)
        
        # EA details
        if results['n_ea'] > 0 and self.print_level >= 1:
            ea_info = {
                'n_ea': results['n_ea'],
                'lumo': f"{results['lumo']:.4f} eV" if results['lumo'] is not None else None,
                'range': [f"{np.min(results['eas']):.4f}", f"{np.max(results['eas']):.4f}"] if results['n_ea'] > 0 else None
            }
            
            # Add detailed list if verbose
            if self.print_level >= 2:
                n_show = min(self.max_print, results['n_ea'])
                ea_details = []
                for i in range(n_show):
                    ea_details.append(f"  {i+1:3d}: {results['eas'][i]:8.4f} eV (pole: {results['ea_poles'][i]:.6f})")
                ea_info['details'] = ea_details
            
            dump_log(self.mol, title='Electron Affinities', section='dyson_ea', info=ea_info)
        
        # Gap
        if results['gap'] is not None:
            dump_log(self.mol, title='', section='dyson_gap', info={'gap': f"{results['gap']:.4f} eV"})
    
    def save_cube_files(self, results):
        """
        Save Dyson orbitals as cube files for visualization.
        """
        if results['n_orbitals'] == 0:
            return
        
        # Get molecular geometry
        atoms = self.mol.get_atoms()
        coords = self.mol.get_system().reshape((-1, 3))
        
        # Check if we can write cube files
        try:
            from oqp.utils.cube_writer import write_cube
            cube_available = True
        except ImportError:
            dump_log(self.mol, title='PyOQP: Cube writer not available, skipping cube file export', section='warning')
            cube_available = False
            return
        
        # Save IP orbitals
        if results['ip_orbitals'] is not None and results['n_ip'] > 0:
            n_save = min(results['n_ip'], 5)  # Save up to 5 IP orbitals
            for i in range(n_save):
                filename = f"{self.mol.log_path}/dyson_ip_{i+1}.cube"
                comment = f"Dyson IP orbital {i+1}, BE={results['ips'][i]:.4f} eV, pole={results['ip_poles'][i]:.6f}"
                write_cube(filename, atoms, coords, results['ip_orbitals'][:, i], comment=comment)
                dump_log(self.mol, title=f'Saved IP Dyson orbital {i+1} to {filename}', section='export')
        
        # Save EA orbitals
        if results['ea_orbitals'] is not None and results['n_ea'] > 0:
            n_save = min(results['n_ea'], 5)  # Save up to 5 EA orbitals
            for i in range(n_save):
                filename = f"{self.mol.log_path}/dyson_ea_{i+1}.cube"
                comment = f"Dyson EA orbital {i+1}, EA={results['eas'][i]:.4f} eV, pole={results['ea_poles'][i]:.6f}"
                write_cube(filename, atoms, coords, results['ea_orbitals'][:, i], comment=comment)
                dump_log(self.mol, title=f'Saved EA Dyson orbital {i+1} to {filename}', section='export')
    
    def generate_spectrum(self, results):
        """
        Generate a spectrum plot of IPs and EAs.
        """
        try:
            import matplotlib.pyplot as plt
            
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
            
            # IP spectrum
            if results['n_ip'] > 0:
                ax1.stem(results['ips'], results['ip_poles'], linefmt='b-', markerfmt='bo', basefmt=' ')
                ax1.set_xlabel('Ionization Potential (eV)')
                ax1.set_ylabel('Pole Strength')
                ax1.set_title(f'IP Spectrum (HOMO = {results["homo"]:.4f} eV)')
                ax1.grid(True, alpha=0.3)
            
            # EA spectrum
            if results['n_ea'] > 0:
                ax2.stem(results['eas'], results['ea_poles'], linefmt='r-', markerfmt='ro', basefmt=' ')
                ax2.set_xlabel('Electron Affinity (eV)')
                ax2.set_ylabel('Pole Strength')
                ax2.set_title(f'EA Spectrum (LUMO = {results["lumo"]:.4f} eV)')
                ax2.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            spectrum_file = self.dyson_config.get('spectrum_file', 'dyson_spectrum.png')
            spectrum_path = f"{self.mol.log_path}/{spectrum_file}"
            plt.savefig(spectrum_path, dpi=300, bbox_inches='tight')
            
            if self.dyson_config.get('show_spectrum', False):
                plt.show()
            else:
                plt.close()
            
            dump_log(self.mol, title=f'Saved spectrum to {spectrum_path}', section='export')
            
        except ImportError:
            dump_log(self.mol, title='PyOQP: Matplotlib not available, skipping spectrum generation', section='warning')
    
    def export_results(self, results):
        """
        Export Dyson results to file in various formats.
        """
        if self.export_format == 'text':
            self.export_text(results)
        elif self.export_format == 'numpy':
            self.export_numpy(results)
        elif self.export_format == 'json':
            self.export_json(results)
        elif self.export_format == 'csv':
            self.export_csv(results)
        
        # Also export using dump_data for consistency with other OQP exports
        if self.export:
            dump_data(self.mol, (results, self.export_title), title='DYSON', fpath=self.mol.log_path)
    
    def export_text(self, results):
        """Export results as text file"""
        filename = f"{self.mol.log_path}/{self.export_name}.txt"
        with open(filename, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write(" " * 25 + "DYSON ORBITAL RESULTS\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Method: {self.td.upper()}\n")
            f.write(f"Target State: {results['target_state']}\n")
            f.write(f"Pole Strength Threshold: {self.threshold}\n")
            f.write(f"Total Significant Orbitals: {results['n_orbitals']}\n")
            f.write(f"Basis Functions: {results['nbf']}\n\n")
            
            if results['n_ip'] > 0:
                f.write("-" * 70 + "\n")
                f.write("IONIZATION POTENTIALS\n")
                f.write("-" * 70 + "\n")
                f.write(f"Number of IPs: {results['n_ip']}\n")
                f.write(f"HOMO: {results['homo']:.6f} eV\n")
                f.write(f"Range: {np.min(results['ips']):.6f} - {np.max(results['ips']):.6f} eV\n\n")
                f.write("  No.      IP (eV)    Pole Strength\n")
                f.write("  ---    ----------   -------------\n")
                for i in range(results['n_ip']):
                    f.write(f"  {i+1:3d}    {results['ips'][i]:10.6f}   {results['ip_poles'][i]:12.8f}\n")
                f.write("\n")
            
            if results['n_ea'] > 0:
                f.write("-" * 70 + "\n")
                f.write("ELECTRON AFFINITIES\n")
                f.write("-" * 70 + "\n")
                f.write(f"Number of EAs: {results['n_ea']}\n")
                f.write(f"LUMO: {results['lumo']:.6f} eV\n")
                f.write(f"Range: {np.min(results['eas']):.6f} - {np.max(results['eas']):.6f} eV\n\n")
                f.write("  No.      EA (eV)    Pole Strength\n")
                f.write("  ---    ----------   -------------\n")
                for i in range(results['n_ea']):
                    f.write(f"  {i+1:3d}    {results['eas'][i]:10.6f}   {results['ea_poles'][i]:12.8f}\n")
                f.write("\n")
            
            if results['gap'] is not None:
                f.write("-" * 70 + "\n")
                f.write(f"HOMO-LUMO Gap: {results['gap']:.6f} eV\n")
                f.write("=" * 70 + "\n")
        
        dump_log(self.mol, title=f'PyOQP: Dyson results exported to {filename}', section='export')
    
    def export_numpy(self, results):
        """Export results as NumPy npz file"""
        filename = f"{self.mol.log_path}/{self.export_name}.npz"
        np.savez(filename, **results)
        dump_log(self.mol, title=f'PyOQP: Dyson results exported to {filename}', section='export')
    
    def export_json(self, results):
        """Export results as JSON file"""
        import json
        filename = f"{self.mol.log_path}/{self.export_name}.json"
        
        # Convert numpy arrays to lists for JSON serialization
        json_results = {}
        for key, value in results.items():
            if isinstance(value, np.ndarray):
                json_results[key] = value.tolist()
            elif value is None:
                json_results[key] = None
            else:
                json_results[key] = value
        
        with open(filename, 'w') as f:
            json.dump(json_results, f, indent=2)
        
        dump_log(self.mol, title=f'PyOQP: Dyson results exported to {filename}', section='export')
    
    def export_csv(self, results):
        """Export results as CSV file"""
        filename = f"{self.mol.log_path}/{self.export_name}.csv"
        
        with open(filename, 'w') as f:
            # Header
            f.write("Type,Index,Energy_eV,Pole_Strength\n")
            
            # IPs
            if results['n_ip'] > 0:
                for i in range(results['n_ip']):
                    f.write(f"IP,{i+1},{results['ips'][i]:.6f},{results['ip_poles'][i]:.8f}\n")
            
            # EAs
            if results['n_ea'] > 0:
                for i in range(results['n_ea']):
                    f.write(f"EA,{i+1},{results['eas'][i]:.6f},{results['ea_poles'][i]:.8f}\n")
        
        dump_log(self.mol, title=f'PyOQP: Dyson results exported to {filename}', section='export')

class Hessian(Calculator):
    """
    OQP frequence calculation class

    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.hess_type = mol.config['hess']['type']
        self.state = mol.config['hess']['state']
        self.read = mol.config['hess']['read']
        self.restart = mol.config['hess']['restart']
        self.temperature = mol.config['hess']['temperature']
        self.clean = mol.config['hess']['clean']

        if self.hess_type == 'analytical':
            self.hess_func = self.analytical_hess
        else:
            self.hess_func = self.numerical_hess

        method = mol.config['input']['method']
        scf_mult = mol.config['scf']['multiplicity']
        td_mult = mol.config['tdhf']['multiplicity']
        if method == 'hf':
            self.hess_mult = scf_mult
        else:
            self.hess_mult = td_mult

    def hessian(self):
        dump_log(self.mol, title='PyOQP: Entering Hessian Calculation')

        if self.read:
            # read .hess file
            dump_log(self.mol, title='', section='read_hess')
            energy, hessian, freqs, modes, inertia = self.mol.read_freqs()

        else:
            # compute hessian
            energy = self.mol.energies[self.state]
            hessian, flags = self.hess_func()
            if 'failed' in flags:
                dump_log(self.mol, title='PyOQP: numerical hessian calculations failed')
                return None
            else:
                freqs, modes, inertia = normal_mode(self.mol.get_system(), self.mol.get_mass(), hessian)
                self.mol.freqs = freqs
                self.mol.hessian = hessian
                self.mol.modes = modes
                self.mol.inertia = inertia

                self.mol.save_freqs(self.state)
                dump_data(self.mol, (self.mol, freqs, modes), title='FREQ', fpath=self.mol.log_path)

                # save mol
                if self.save_mol:
                    self.mol.save_data()

        dump_log(self.mol, title='PyOQP: Frequencies', section='freq', info=freqs)

        for t in self.temperature:
            thermal_data = thermal_analysis(
                energy=energy,
                atoms=self.mol.get_atoms(),
                mass=self.mol.get_mass(),
                freqs=freqs,
                inertia=inertia,
                temperature=t,
                mult=self.hess_mult,
            )
            dump_log(self.mol, title='PyOQP: Thermochemistry at %-10.2f K' % t, section='thermo', info=thermal_data)

    def analytical_hess(self):
        exit('analytical hessian is not available yet, choose numerical')

    def numerical_hess(self):
        dir_hess = f'{self.mol.log_path}/{self.mol.project_name}_num_hess'
        nproc = self.mol.config['hess']['nproc']
        dx = self.mol.config['hess']['dx']
        origin_coord = self.mol.get_system()

        # prepare scratch folder
        os.makedirs(dir_hess, exist_ok=True)

        # shift origin 3N coord with 6N displacement
        ncoord = len(origin_coord)
        shift = np.diag(np.ones(ncoord) * dx).reshape(ncoord, ncoord)
        shifted_coord = np.concatenate((origin_coord + shift, origin_coord - shift), axis=0)
        ndim = len(shifted_coord)

        # prepare grad calculations
        self.mol.save_data()
        self.mpi_manager.barrier()
        atoms = self.mol.get_atoms()
        guess_file = self.mol.log.replace('.log', '.json')
        variables_wrapper = [
            {
                'idx': idx,
                'atoms': atoms,
                'coord': coord,
                'dir_hess': dir_hess,
                'project_name': self.mol.project_name,
                'config': copy.deepcopy(self.mol.config),
                'guess_file': guess_file,
                'state': self.state,
                'restart': self.restart,
            }
            for idx, coord in enumerate(shifted_coord)
        ]

        ## adjust multiprocessing if necessary
        if self.mpi_manager.use_mpi:
            ncpu = np.amin([ndim, self.mpi_manager.comm.size])
            pool = MPIPool(processes=ncpu)
        else:
            ncpu = np.amin([ndim, nproc])
            pool = multiprocessing.Pool(processes=ncpu)

        dump_log(self.mol,
                 title='',
                 section='num_hess',
                 info=[self.state, ndim, dx, self.restart, len(variables_wrapper), ncpu, os.environ['OMP_NUM_THREADS']]
                 )

        ## start multiprocessing
        grads = [[] for _ in range(ndim)]
        flags = []
        n = 0
        for val in pool.imap_unordered(grad_wrapper, variables_wrapper):
            if self.mpi_manager.rank == 0:
                n += 1
                idx, grad, flag, timing = val
                grads[idx] = grad
                flags.append(flag)
                dump_log(self.mol, title=None, section='hess_worker', info=[n, idx, flag, timing])
                dump_data(self.mol, (n, idx, np.sum(grad ** 2) ** 2, timing), title='NUM_HESS', fpath=self.mol.log_path)

        pool.close()

        grads = self.mpi_manager.bcast(grads)
        # compute hessian
        forward = np.array(grads[0:ncoord])
        backward = np.array(grads[ncoord:])
        hessian = (forward - backward) / (2 * dx)

        # symmetrize hessian
        hessian = (hessian + hessian.T) / 2

        # delete scratch folder
        if 'failed' not in flags and self.clean and self.mpi_manager.rank == 0:
            shutil.rmtree(dir_hess)

        return hessian, flags


def grad_wrapper(key_dict):
    start_time = time.time()
    rank = MPIManager().rank
    threads = os.environ['OMP_NUM_THREADS']
    host = platform.node()
    # unpack variables
    idx = key_dict['idx']
    atoms = key_dict['atoms']
    coord = key_dict['coord']
    dir_hess = key_dict['dir_hess']
    project_name = key_dict['project_name']
    config = key_dict['config']
    guess_file = key_dict['guess_file']
    state = key_dict['state']
    restart = key_dict['restart']

    # prepare log files
    inp = f'{dir_hess}/{project_name}.{idx}.tmp.inp'
    xyz = f'{dir_hess}/{project_name}.{idx}.tmp.xyz'
    dat = f'{dir_hess}/{project_name}.{idx}.grad_{state}'
    log = f'{dir_hess}/{project_name}.{idx}.tmp.log'

    # attempt to read computed data
    if restart and os.path.exists(dat):
        status = 'loaded'
    else:
        status = 'computed'

        # modify config
        config['input']['runtype'] = 'grad'
        config['input']['system'] = xyz
        config['guess']['type'] = 'json'
        config['guess']['file'] = guess_file
        config['guess']['continue_geom'] = 'false'
        config['properties']['grad'] = config['hess']['state']
        config['properties']['export'] = 'True'
        config['properties']['title'] = f'{project_name}.{idx}'
        config['hess']['temperature'] = ','.join([str(x) for x in config['hess']['temperature']])
        config['tests']['exception'] = 'false'

        # save config
        input_xyz = write_xyz(atoms, coord, [idx])
        input_file, input_dict = write_config(config)

        with open(xyz, 'w') as out:
            out.write(input_xyz)

        with open(inp, 'w') as out:
            out.write(input_file)

        if not MPIManager().use_mpi:
            # run grad calculation externally
            subprocess.run(['openqp', inp, '--silent'])
        else:
            # run grad calculation internally
            start_time = time.time()
            mol = Molecule(project_name, inp, log, silent=1)
            mol.usempi = False
            mol.load_config(input_dict)
            mol.load_data()
            mol.start_time = start_time
            dump_log(mol, title='', section='start')
            mol.data["OQP::log_filename"] = log
            oqp.oqp_banner(mol)
            SinglePoint(mol).energy()
            Gradient(mol).gradient()
            LastStep(mol).compute(mol, grad_list=mol.config['properties']['grad'])
            dump_log(mol, title='', section='end')
    try:
        grad = np.loadtxt(dat).reshape(-1)

    except FileNotFoundError:
        grad = np.zeros_like(coord)
        status = 'failed'

    end_time = time.time()

    return idx, grad, status, (start_time, end_time, rank, threads, host)


class BasisOverlap(Calculator):
    """
    OQP basis overlap calculation class
    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.natom = mol.data["natom"]
        self.nstate = mol.config['tdhf']['nstate']
        self.back_door = self.mol.config['properties']['back_door']
        self.align_type = self.mol.config['nac']['align']
        self.overlap_func = oqp.get_structures_ao_overlap

    def overlap(self):
        # load previous data
        self.load_previous_data()

        # compute basis overlap
        self.overlap_func(self.mol)

        # align mo before tdhf
        if self.align_type != 'no':
            self.align_mo()

    def load_previous_data(self):
        dump_log(self.mol, title='PyOQP: Loading Previous Data')
        if self.back_door:
            # get data externally
            previous_xyz, previous_data = self.mol.get_data_from_back_door()
        else:
            # record data from the current step
            current_file = self.mol.config["guess"]["file"]
            current_xyz = copy.deepcopy(self.mol.get_system())
            current_data = copy.deepcopy(self.mol.get_data())

            # check data from the previous step
            previous_coord = self.mol.data.mol2
            previous_file = self.mol.config['guess']['file2']

            if previous_file:
                self.mol.config['guess']['file'] = previous_file
                self.mol.config['guess']['continue_geom'] = True
                self.mol.load_data()
                previous_xyz = copy.deepcopy(self.mol.get_system())
                previous_data = copy.deepcopy(self.mol.get_data())
            else:
                if len(previous_coord) > 1:
                    # compute data for previous step
                    self.mol.idx = 2
                    self.mol.update_system(previous_coord)
                    oqp.library.ints_1e(self.mol)
                    oqp.library.guess(self.mol)
                    SinglePoint(self.mol).energy()
                    LastStep(self.mol).compute(self.mol)
                    previous_xyz = previous_coord
                    previous_data = copy.deepcopy(self.mol.get_data())
                else:
                    previous_xyz = None
                    previous_data = None
                    exit(f'\nmolecule loads previous step data cannot find [guess] file2 or [input] system2')

            # restore current data
            self.mol.idx = 1
            self.mol.config['guess']['file'] = current_file
            self.mol.update_system(current_xyz)
            self.mol.put_data(current_data)

        # copy previous data to old tags
        natom = self.mol.data["natom"]
        self.mol.data["OQP::xyz_old"] = previous_xyz.reshape((3, natom))
        self.mol.data["OQP::VEC_MO_A_old"] = previous_data["OQP::VEC_MO_A"]
        self.mol.data["OQP::VEC_MO_B_old"] = previous_data["OQP::VEC_MO_B"]
        self.mol.data["OQP::E_MO_A_old"] = previous_data["OQP::E_MO_A"]
        self.mol.data["OQP::E_MO_B_old"] = previous_data["OQP::E_MO_B"]

        # check if td data are available
        try:
            self.mol.data["OQP::td_bvec_mo_old"] = previous_data["OQP::td_bvec_mo"]
            self.mol.data["OQP::td_energies_old"] = previous_data["OQP::td_energies"]
        except KeyError:
            pass

    def align_mo(self):
        dump_log(self.mol, title='PyOQP: Entering Overlap Calculation', section='basis_overlap')
        current_mo = copy.deepcopy(self.mol.data["OQP::VEC_MO_A"])
        current_energy = copy.deepcopy(self.mol.data["OQP::E_MO_A"])
        nocc = self.mol.data['nelec_A']

        # get MO overlap data
        mo_overlap_matrix = self.mol.data["OQP::overlap_mo_non_orthogonal"]

        # current MO in row, previous MO in column
        occ_order, occ_sign = self.find_vec_order(mo_overlap_matrix[:nocc - 2, : nocc - 2])
        somo_order, somo_sign = self.find_vec_order(mo_overlap_matrix[nocc - 2: nocc, nocc - 2: nocc])
        vir_order, vir_sign = self.find_vec_order(mo_overlap_matrix[nocc:, nocc:])

        mo_order = np.concatenate((occ_order, somo_order + nocc - 2, vir_order + nocc))
        mo_sign = np.concatenate((occ_sign, somo_sign, vir_sign)).reshape((-1, 1))

        # apply sign correction
        current_mo *= mo_sign

        # reorder mo if requested
        if self.align_type == 'reorder':
            current_mo = current_mo[np.argsort(mo_order)]
            current_energy = current_energy[np.argsort(mo_order)]

        # update mo
        self.mol.data["OQP::VEC_MO_A"] = current_mo
        self.mol.data["OQP::VEC_MO_B"] = current_mo
        self.mol.data["OQP::E_MO_A"] = current_energy
        self.mol.data["OQP::E_MO_B"] = current_energy

        dump_log(self.mol, title='PyOQP: Aligning MOs')

        # compute new basis overlap
        self.overlap_func(self.mol)

    @staticmethod
    def find_vec_order(overlap_matrix):
        overlap_matrix = copy.deepcopy(overlap_matrix)
        vec_order = []
        vec_sign = []
        for s_i in overlap_matrix:
            # find the matched vec and sign
            n_i = np.argmax(np.abs(s_i))
            p_i = np.sign(s_i)[n_i]

            # record the vec index and sign
            vec_order.append(n_i)
            vec_sign.append(p_i)

            # remove the matched vec
            overlap_matrix[:, n_i] = 0

        return np.array(vec_order), np.array(vec_sign)


class NACME(BasisOverlap):
    """
    Class to calculate Non-Adiabatic Coupling (NAC) matrix elements
    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.dt = mol.config['nac']['dt']
        self.nac_dim = self.nstate

    def align_x(self):
        # use current td data if previous td data is unavailable
        try:
            self.mol.data["OQP::td_bvec_mo_old"]
        except AttributeError:
            self.mol.data["OQP::td_bvec_mo_old"] = self.mol.data["OQP::td_bvec_mo"]
            self.mol.data["OQP::td_energies_old"] = self.mol.data["OQP::td_energies"]

        previous_x = copy.deepcopy(self.mol.data['OQP::td_bvec_mo_old'])
        current_x = copy.deepcopy(self.mol.data['OQP::td_bvec_mo'])
        x_shape = current_x.shape

        # reshape Fortran data into Python style
        current_x = current_x.reshape((x_shape[1], x_shape[0]))
        previous_x = previous_x.reshape((x_shape[1], x_shape[0]))

        # current x in row, previous x in column
        x_overlap_matrix = np.matmul(current_x, previous_x.T)
        x_order, x_sign = self.find_vec_order(x_overlap_matrix)

        # apply sign correction
        current_x *= x_sign.reshape((-1, 1))

        # update x in Fortran data shape
        self.mol.data['OQP::td_bvec_mo'] = current_x.reshape((x_shape[0], x_shape[1]))

        dump_log(self.mol, title='PyOQP: Aligning X amplitudes')

    def nacme(self):
        """
        Calculates the non-adiabatic coupling (NAC) matrix elements
        between the two geometries.

        Currently, adapted only for MRSF-TDDFT approach
        """
        dump_log(self.mol, title='PyOQP: Entering State Overlap Calculation')

        # align X amplitudes
        self.align_x()

        # compute state overlap
        oqp.get_states_overlap(self.mol)
        state_overlap = self.mol.data["OQP::td_states_overlap"]

        # compute time-derivative nac
        dc_matrix = (state_overlap - state_overlap.T) / self.dt
        e_i = np.array(self.mol.energies[1:]).reshape((-1, 1))
        e_j = np.array(self.mol.energies[1:]).reshape((1, -1))
        gap = e_j - e_i
        nac_matrix = dc_matrix * gap
        nac_matrix = nac_matrix[0:self.nac_dim, 0:self.nac_dim]
        self.mol.data["OQP::dc_matrix"] = dc_matrix
        self.mol.data["OQP::nac_matrix"] = nac_matrix

        dump_log(self.mol, title='PyOQP: Non-Adiabatic Coupling Matrix Calculation', section='nacme')
        dump_log(self.mol, title='PyOQP: phase corrected state overlap (s_ij)', section='nacm', info=state_overlap)
        dump_log(self.mol, title='PyOQP: phase corrected derivative coupling (d_ij)', section='nacm', info=dc_matrix)
        dump_log(self.mol, title='PyOQP: state energy gap (e_ji)', section='nacm', info=gap)
        dump_log(self.mol, title='PyOQP: phase corrected non-adiabatic coupling (h_ij)', section='nacm',
                 info=nac_matrix)

        # save corrected data
        self.mol.save_data()
        self.mol.dcm = dc_matrix

        # export data
        if self.export:
            dump_data(self.mol, (nac_matrix, self.export_title), title='NACME', fpath=self.mol.log_path)
            dump_data(self.mol, (dc_matrix, self.export_title), title='DCME', fpath=self.mol.log_path)

        return dc_matrix, nac_matrix


class NAC(Calculator):
    """
    Class to calculate Non-Adiabatic Coupling (NAC) vector
    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.nac_type = mol.config["nac"]["type"]
        self.natom = mol.data["natom"]
        self.nstate = mol.config['tdhf']['nstate']
        self.restart = mol.config['nac']['restart']
        self.clean = mol.config['nac']['clean']
        self.nac_states = mol.config['nac']['states']
        self.bp = mol.config['nac']['bp']
        self.nac_func = oqp.get_states_overlap

        if self.nac_type == 'analytical':
            self.nac_func = self.analytical_nac
        else:
            self.nac_func = self.numerical_nac

    def nac(self):
        dump_log(self.mol, title='PyOQP: Entering NAC Vector Calculation')

        # compute nacv
        nacv, dcv, flags = self.nac_func()
        if 'failed' in flags:
            dump_log(self.mol, title='PyOQP: numerical nac calculations failed')
            return None
        else:
            self.mol.nac = nacv

            # save mol
            if self.save_mol:
                self.mol.save_data()

            # export data
            if self.export:
                dump_data(self.mol, (self.mol, nacv, self.export_title), title='NACV', fpath=self.mol.log_path)

        dump_log(self.mol, title='PyOQP: NAC Vector (h_ij)', section='nacv', info=nacv)
        dump_log(self.mol, title='PyOQP: DC Vector (d_ij)', section='dcv', info=dcv)

        # compute gradients for branching plane
        if self.bp:
            dump_log(self.mol, title='PyOQP: Entering Branching Plane Calculation')
            grad = Gradient(self.mol)
            grad.grads = np.unique(self.nac_states)
            grads = grad.gradient()
            for ij in self.nac_states:
                i, j = np.sort(ij)
                g1 = grads[i]
                g2 = grads[j]
                h = nacv[i - 1, j - 1]
                x, y, sx, sy, pitch, tilt, peak, bifu = self.compute_bp(h, g1, g2)
                dump_log(self.mol,
                         title='PyOQP: Final Gradient',
                         section='grad',
                         info={'el': grads, 'd4': np.zeros_like(g1), 'grad_list': [i, j]}
                         )
                dump_log(self.mol,
                         title='PyOQP: Branching Plane Info %s - %s' % (i, j),
                         section='bp',
                         info=(g1, g2, h, x, y, sx, sy, pitch, tilt, peak, bifu)
                         )
                dump_data(self.mol, (self.mol, g1, g2, h, x, y, i, j), title='BP', fpath=self.mol.log_path)

    @staticmethod
    def compute_bp(h, g1, g2):
        s = (g2 + g1) / 2
        g = (g2 - g1) / 2

        tan2b = 2 * np.sum(g * h) / (np.sum(g ** 2) - np.sum(h ** 2))
        b = np.arctan(tan2b) / 2
        sinb = np.sin(b)
        cosb = np.cos(b)
        gt = g * cosb + h * sinb
        ht = h * cosb - g * sinb

        x = gt / la.norm(gt)
        y = ht / la.norm(ht)

        ngt = np.sum(gt * gt)
        nht = np.sum(ht * ht)

        pitch = (0.5 * (ngt + nht)) ** 0.5
        tilt = (ngt - nht) / (ngt + nht)

        sx = np.sum(s * x)
        sy = np.sum(s * y)

        th = ((sx / pitch) ** 2 + (sy / pitch) ** 2) ** 0.5
        ts = np.arctan(sy / sx)

        peak = th ** 2 / (1 - tilt ** 2) * (1 - tilt * np.cos(2 * ts))
        bifu = (th ** 2 / (4 * tilt ** 2)) ** (1 / 3) * \
               (
                       ((1 + tilt) * np.cos(ts) ** 2) ** (1 / 3) +
                       ((1 - tilt) * np.sin(ts) ** 2) ** (1 / 3)
               )
        return x, y, sx, sy, pitch, tilt, peak, bifu

    def analytical_nac(self):
        exit('analytical nac vector is not available yet, choose numerical')

    def numerical_nac(self):
        dir_nacv = f'{self.mol.log_path}/{self.mol.project_name}_num_nacv'
        nproc = self.mol.config['nac']['nproc']
        dx = self.mol.config['nac']['dx']
        origin_coord = self.mol.get_system()

        # prepare scratch folder
        os.makedirs(dir_nacv, exist_ok=True)

        # shift origin 3N coord with 6N displacement
        ncoord = len(origin_coord)
        shift = np.diag(np.ones(ncoord) * dx).reshape(ncoord, ncoord)
        shifted_coord = np.concatenate((origin_coord + shift, origin_coord - shift), axis=0)
        ndim = len(shifted_coord)

        # prepare grad calculations
        self.mol.save_data()
        self.mpi_manager.barrier()
        atoms = self.mol.get_atoms()
        guess_file = self.mol.log.replace('.log', '.json')
        variables_wrapper = [
            {
                'idx': idx,
                'dx': dx,
                'atoms': atoms,
                'coord': coord,
                'nstate': self.nstate,
                'dir_nacv': dir_nacv,
                'project_name': self.mol.project_name,
                'config': copy.deepcopy(self.mol.config),
                'guess_file': guess_file,
                'restart': self.restart,
            }
            for idx, coord in enumerate(shifted_coord)
        ]

        ## adjust multiprocessing if necessary
        if self.mpi_manager.use_mpi:
            ncpu = np.amin([ndim, self.mpi_manager.comm.size])
            pool = MPIPool(processes=ncpu)
        else:
            ncpu = np.amin([ndim, nproc])
            pool = multiprocessing.Pool(processes=ncpu)

        dump_log(self.mol,
                 title='',
                 section='num_nacv',
                 info=[ndim, dx, self.restart, len(variables_wrapper), ncpu, os.environ['OMP_NUM_THREADS']]
                 )

        ## start multiprocessing
        dcm = [[] for _ in range(ndim)]
        flags = []
        n = 0
        for val in pool.imap_unordered(nacme_wrapper, variables_wrapper):
            if self.mpi_manager.rank == 0:
                n += 1
                idx, dcme, flag, timing = val
                dcm[idx] = dcme.reshape(-1)  # 1D nstate x nstate
                flags.append(flag)
                dump_log(self.mol, title=None, section='nacv_worker', info=[n, idx, flag, timing])
                dump_data(self.mol, (n, idx, np.sum(dcme ** 2) ** 2, timing), title='NUM_NACV', fpath=self.mol.log_path)

        pool.close()

        dcm = self.mpi_manager.bcast(dcm)
        # compute nacv (natom x 3, nstate x nstate)
        forward = np.array(dcm[0:ncoord])
        backward = np.array(dcm[ncoord:])
        dcm = (forward - backward) / 2
        e_i = np.array(self.mol.energies[1:]).reshape((1, -1))
        e_j = np.array(self.mol.energies[1:]).reshape((-1, 1))
        gap = e_j - e_i
        np.fill_diagonal(gap, 1)
        gap = gap.reshape((1, -1))
        nacm = dcm * gap

        # reshape matrix -> (nstate x nstate, natom x 3) -> (nstate, nstate, natom, 3)
        dcv = dcm.T.reshape((self.nstate, self.nstate, self.natom, 3))
        nacv = nacm.T.reshape((self.nstate, self.nstate, self.natom, 3))

        # delete scratch folder
        if 'failed' not in flags and self.clean and self.mpi_manager.rank == 0:
            shutil.rmtree(dir_nacv)

        return nacv, dcv, flags


def nacme_wrapper(key_dict):
    start_time = time.time()
    rank = MPIManager().rank
    threads = os.environ['OMP_NUM_THREADS']
    host = platform.node()
    # unpack variables
    idx = key_dict['idx']
    dx = key_dict['dx']
    atoms = key_dict['atoms']
    coord = key_dict['coord']
    nstate = key_dict['nstate']
    dir_nacv = key_dict['dir_nacv']
    project_name = key_dict['project_name']
    config = key_dict['config']
    guess_file = key_dict['guess_file']
    restart = key_dict['restart']

    # prepare log files
    inp = f'{dir_nacv}/{project_name}.{idx}.tmp.inp'
    xyz = f'{dir_nacv}/{project_name}.{idx}.tmp.xyz'
    dat = f'{dir_nacv}/{project_name}.{idx}.dcme'
    log = f'{dir_nacv}/{project_name}.{idx}.tmp.log'

    # attempt to read computed data
    if restart and os.path.exists(dat):
        status = 'loaded'
    else:
        status = 'computed'

        # modify config
        config['input']['runtype'] = 'nacme'
        config['input']['system'] = xyz
        config['guess']['type'] = 'json'
        config['guess']['file'] = guess_file
        config['guess']['file2'] = guess_file
        config['guess']['continue_geom'] = 'false'
        config['properties']['export'] = 'true'
        config['properties']['title'] = f'{project_name}.{idx}'
        config['nac']['dt'] = str(dx)
        config['tests']['exception'] = 'false'

        # save config
        input_xyz = write_xyz(atoms, coord, [idx])
        input_file, input_dict = write_config(config)

        with open(xyz, 'w') as out:
            out.write(input_xyz)

        with open(inp, 'w') as out:
            out.write(input_file)

        if not MPIManager().use_mpi:
            # run nac calculation externally
            subprocess.run(['openqp', inp, '--silent'])
        else:
            # run nac calculation internally
            start_time = time.time()
            mol = Molecule(project_name, inp, log, silent=1)
            mol.usempi = False
            mol.load_config(input_dict)
            mol.load_data()
            mol.start_time = start_time
            dump_log(mol, title='', section='start')
            mol.data["OQP::log_filename"] = log
            oqp.oqp_banner(mol)
            sp = SinglePoint(mol)
            ref_energy = sp.reference()
            BasisOverlap(mol).overlap()
            sp.excitation(ref_energy)
            LastStep(mol).compute(mol)
            NACME(mol).nacme()
            dump_log(mol, title='', section='end')

    try:
        dcme = np.loadtxt(dat).reshape(-1)

    except FileNotFoundError:
        dcme = np.zeros_like(nstate * nstate)
        status = 'failed'

    end_time = time.time()

    return idx, dcme, status, (start_time, end_time, rank, threads, host)


class SCFnotConverged(Exception):
    pass


class TDnotConverged(Exception):
    pass


class ZVnotConverged(Exception):
    pass
