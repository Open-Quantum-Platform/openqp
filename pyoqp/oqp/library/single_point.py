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
            dump_data(self.mol, (self.mol.grads, self.export_title, grad_list), title='GRADIENT', fpath=self.mol.log_path)

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
        self.init_conv =  mol.config['scf']['init_conv']
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
        else :
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
        ixcore= self.mol.config["tdhf"]["ixcore"]
        
        if ixcore == "-1":  # if default
            # Pass pointer 
            ixcore_array = np.array([-1], dtype=np.int32)
            self.mol.data['ixcore'] = ffi.cast("int32_t *", ffi.from_buffer(ixcore_array))
            self.mol.data['ixcore_len'] = ixcore_array.size
            return

        # Pass pointer to Fortran via C
        ixcore_array = np.array(ixcore.split(','), dtype=np.int32)
        self.mol.data['ixcore'] = ffi.cast("int32_t *", ffi.from_buffer(ixcore_array))
        self.mol.data['ixcore_len'] = ixcore_array.size

        # shift MO energies 
        noccB = self.mol.data['nelec_B']
        tmp = self.mol.data["OQP::E_MO_A"]
        for i in range(noccB+1):  # up to HOMO-1
            if i not in ixcore_array:
                tmp[i-1] = -100000  # shift the MO energy down

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
        swapmo= self.mol.config["guess"]["swapmo"]
        if swapmo:   # if not default (empty)
            swapmo_array = [int(x.strip()) for x in swapmo.split(',')]

            # Initial MO energy and coefficient
            og_val = self.mol.data["OQP::E_MO_A"]
            og_vec = self.mol.data["OQP::VEC_MO_A"]
            # It only takes pairs. If it is not pair, it will be ignored.
            for i, j in zip(swapmo_array[::2], swapmo_array[1::2]):
                og_val[[i-1, j-1]] = og_val[[j-1, i-1]]
                og_vec[[i-1, j-1]] = og_vec[[j-1, i-1]]

    def reference(self, do_init_scf=True):
        dump_log(self.mol, title='PyOQP: Entering Electronic Energy Calculation', section='input')

        if self.init_scf != 'no' and do_init_scf:
            # do initial scf iteration to help convergence
            self._init_convergence()
        else:
            self._prep_guess()

        self.swapmo()

        scf_flag = False
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
        dump_log(self.mol, title='PyOQP: phase corrected non-adiabatic coupling (h_ij)', section='nacm', info=nac_matrix)

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
        ts = np.arctan(sy/sx)

        peak = th ** 2 / (1 - tilt ** 2) * (1 - tilt * np.cos(2 * ts))
        bifu = (th ** 2 / (4 * tilt ** 2)) ** (1/3) *\
               (
                       ((1 + tilt) * np.cos(ts) ** 2) ** (1/3) +
                       ((1 - tilt) * np.sin(ts) ** 2) ** (1/3)
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
