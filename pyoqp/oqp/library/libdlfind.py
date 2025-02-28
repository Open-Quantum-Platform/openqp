"""OQP mini-app with DL-FIND for geometry optimization"""

import os
import time
import signal
import oqp
import functools
import numpy as np
from oqp.library.libscipy import Optimizer
from oqp.utils.file_utils import dump_data, dump_log
from libdlfind import dl_find
from libdlfind.callback import (
    dlf_get_gradient_wrapper,
    dlf_get_multistate_gradients_wrapper,
    dlf_put_coords_wrapper,
    make_dlf_get_params,
)


@dlf_get_gradient_wrapper
def get_grad(coordinates, *args, calculator=None, stopper=None):
    try:
        energy, grad = calculator(coordinates)
        return energy, grad

    except StopIteration:
        stopper()


@dlf_get_multistate_gradients_wrapper
def get_msg(coordinates, *args, calculator=None, stopper=None):
    try:
        energies, grads = calculator(coordinates)
        return energies[0], energies[1], grads[0], grads[1], None

    except StopIteration:
        stopper()


@dlf_put_coords_wrapper
def store_results(*args):
    # we implemented our-own result printing methods
    pass

class DLFinder(Optimizer):
    """
    OQP DL-Find optimization class
    user should use the subclass to define the dlf_get_params and dlf_get_gradient function

    """

    def __init__(self, mol):
        super().__init__(mol)
        self.printl = mol.config['dlfind']['printl']
        self.icoord = mol.config['dlfind']['icoord']
        self.iopt = mol.config['dlfind']['iopt']
        self.ims = mol.config['dlfind']['ims']

        self.dlf_get_params = make_dlf_get_params(
            coords=self.pre_coord.reshape((self.natom, 3)),
            printl=self.printl,
            icoord=self.icoord,
            iopt=self.iopt,
            imultistate=self.ims,
            maxcycle=self.maxit,
            tolerance=1e-20,
            tolerance_e=1e-20,
        )

    def optimize(self):
        dl_find(
            nvarin=int(self.natom * 3),
            dlf_get_params=self.dlf_get_params,
            dlf_get_gradient=functools.partial(get_grad, calculator=self.opt_func, stopper=self.early_stop),
            dlf_put_coords=store_results,
        )

    def early_stop(self):
        # we implemented our-own convergence checker, just early stop the dlfind optimization
        # I have to kill the pid since the fortran subroutine ignores the exception from python callback
        self.mol.end_time = time.time()
        dump_log(self.mol, title='', section='end')
        pid = os.getpid()
        print('OQP early stops DL-FIND')
        os.kill(pid, signal.SIGKILL)

class DLFindMin(DLFinder):
    """
    OQP DL-Find state specific optimization class

    """
    def __init__(self, mol):
        super().__init__(mol)

        dump_log(self.mol, title='PyOQP: DL-FIND Local Minimum Optimization', section='dlf')

    def one_step(self, coordinates):
        # flatten coord
        coordinates = coordinates.reshape(-1)

        # add iteration
        self.itr += 1

        dump_log(self.mol, title='PyOQP: Geometry Optimization Step %s' % self.itr)

        if self.itr == 1:
            do_init_scf = True
        else:
            do_init_scf = self.init_scf

        # update coordinates
        self.mol.update_system(coordinates)

        # compute energy
        energies = self.sp.energy(do_init_scf=do_init_scf)

        # compute gradient
        self.grad.grads = [self.istate]
        grads = self.grad.gradient()
        self.mol.energies = energies
        self.mol.grads = grads

        # compute dftd4
        energies, grads = self.ls.compute(self.mol, grad_list=self.grad.grads)
        self.mol.energies = energies
        self.mol.grads = grads

        # flatten data
        energy = energies[self.istate]
        grad = grads[self.istate].reshape(-1)

        # evaluate metrics
        de = energy - self.pre_energy
        rmsd_step = np.mean((coordinates - self.pre_coord) ** 2) ** 0.5
        max_step = np.amax(np.abs(coordinates - self.pre_coord))
        rmsd_grad = np.mean(grad ** 2) ** 0.5
        max_grad = np.amax(np.abs(grad))
        self.metrics['itr'] = self.itr
        self.metrics['de'] = de
        self.metrics['rmsd_step'] = rmsd_step
        self.metrics['max_step'] = max_step
        self.metrics['rmsd_grad'] = rmsd_grad
        self.metrics['max_grad'] = max_grad

        # store energy and coordinates
        self.pre_energy = energy
        self.pre_coord = coordinates.copy()
        dump_data(self.mol,
                  (self.itr, self.atoms, coordinates, energy, de, rmsd_step, max_step, rmsd_grad, max_grad),
                  title='OPTIMIZATION', fpath=self.mol.log_path)

        return energy, grad.reshape((self.natom, 3))

class DLFindTS(DLFindMin):
    """
    OQP DL-Find transition state optimization class

    """
    def __init__(self, mol):
        super().__init__(mol)

        dump_log(self.mol, title='PyOQP: DL-FIND Transition State Optimization', section='dlf')

class DLFindMECI(DLFinder):
    """
    OQP DL-Find MECI optimization class

    """
    def __init__(self, mol):
        super().__init__(mol)

        dump_log(self.mol, title='PyOQP: DL-FIND Minimum Energy Conical Intersection Optimization', section='dlf')

    def one_step(self, coordinates):
        # flatten coord
        coordinates = coordinates.reshape(-1)

        # add iteration
        self.itr += 1

        dump_log(self.mol, title='PyOQP: Geometry Optimization Step %s' % self.itr)

        if self.itr == 1:
            do_init_scf = True
        else:
            do_init_scf = self.init_scf

        # update coordinates
        self.mol.update_system(coordinates)

        # compute energy
        energies = self.sp.energy(do_init_scf=do_init_scf)

        # compute gradient
        self.grad.grads = [self.istate]
        grads = self.grad.gradient()
        self.mol.energies = energies
        self.mol.grads = grads

        # compute dftd4
        energies, grads = self.ls.compute(self.mol, grad_list=self.grad.grads)
        self.mol.energies = energies
        self.mol.grads = grads

        # flatten data
        energy_1 = energies[self.istate]
        energy_2 = energies[self.jstate]
        grad_1 = grads[self.istate].reshape(-1)
        grad_2 = grads[self.jstate].reshape(-1)
        sum_e = (energy_1 + energy_2) / 2
        gap_e = energy_2 - energy_1
        gap_g = grad_1 - grad_2
        dg = 2 * gap_e * gap_g
        # evaluate metrics
        de = sum_e - self.pre_energy
        rmsd_step = np.mean((coordinates - self.pre_coord) ** 2) ** 0.5
        max_step = np.amax(np.abs(coordinates - self.pre_coord))
        rmsd_grad = np.mean(dg) ** 0.5
        max_grad = np.amax(np.abs(dg))
        self.metrics['itr'] = self.itr
        self.metrics['de'] = de
        self.metrics['gap'] = gap_e
        self.metrics['rmsd_step'] = rmsd_step
        self.metrics['max_step'] = max_step
        self.metrics['rmsd_grad'] = rmsd_grad
        self.metrics['max_grad'] = max_grad

        # store energy and coordinates
        self.pre_energy = sum_e
        self.pre_coord = coordinates.copy()
        dump_data(self.mol,
                  (self.itr, self.atoms, coordinates, sum_e, de, gap_e, rmsd_step, max_step, rmsd_grad, max_grad, 0, 0),
                  title='MECI', fpath=self.mol.log_path)

        return (energy_1, energy_2), (grad_1.reshape((self.natom, 3)), grad_2.reshape((self.natom, 3)))

    def optimize(self):
        dl_find(
            nvarin=int(self.natom * 3),
            dlf_get_params=self.dlf_get_params,
            dlf_get_multistate_gradients=functools.partial(get_msg, calculator=self.one_step, stopper=self.early_stop),
            dlf_put_coords=store_results,
        )

    def check_convergence(self):
        # write convergence to log
        dump_log(
            self.mol, title='Geometry Optimization Convergence %s' % self.itr,
            section='dlf-ci',
            info=self.metrics
        )

        if np.abs(self.metrics['de']) <= self.energy_shift and \
                self.metrics['gap'] <= self.energy_gap and \
                self.metrics['rmsd_step'] <= self.rmsd_step and \
                self.metrics['max_step'] <= self.rmsd_step and \
                self.metrics['rmsd_grad'] <= self.rmsd_grad and \
                self.metrics['max_grad'] <= self.max_grad:
            dump_log(self.mol, title='PyOQP: Geometry Optimization Has Converged')
            raise StopIteration

        else:
            if self.itr == self.maxit:
                dump_log(self.mol,
                         title='PyOQP: Geometry Optimization Has Not Converged. Reached The Maximum Iteration')
                raise StopIteration
