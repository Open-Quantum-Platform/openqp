"""OQP single point class"""
import copy
import oqp
import numpy as np
import scipy as sc
from oqp.library.single_point import SinglePoint, Gradient, LastStep
from oqp.utils.file_utils import dump_log, dump_data
import oqp.utils.qmmm as qmmm

class Optimizer:
    """
    OQP optimization class

    user should use the subclass to define the one_step optimization function
    see TemplateOpt class for example

    """

    def __init__(self, mol):
        self.mol = mol
        self.optimizer = mol.config['optimize']['optimizer']
        self.step_size = mol.config['optimize']['step_size']
        self.step_tol = mol.config['optimize']['step_tol']
        self.maxit = mol.config['optimize']['maxit']
        self.mep_maxit = mol.config['optimize']['mep_maxit']
        self.rmsd_grad = mol.config['optimize']['rmsd_grad']
        self.rmsd_step = mol.config['optimize']['rmsd_step']
        self.max_grad = mol.config['optimize']['max_grad']
        self.max_step = mol.config['optimize']['max_step']
        self.istate = mol.config['optimize']['istate']
        self.jstate = mol.config['optimize']['jstate']
        self.kstate = mol.config['optimize']['kstate']
        self.imult = mol.config['optimize']['imult']
        self.jmult = mol.config['optimize']['jmult']
        self.energy_shift = mol.config['optimize']['energy_shift']
        self.energy_gap = mol.config['optimize']['energy_gap']
        self.init_scf = mol.config['optimize']['init_scf']
        self.nstate = mol.config['tdhf']['nstate']
        self.natom = mol.data['natom']
        self.atoms = mol.get_atoms().reshape((self.natom, 1))
        self.sp = SinglePoint(mol)
        self.grad = Gradient(mol)
        self.ls = LastStep(mol)
        self.itr = 0
        self.pre_energy = 0
        self.pre_coord = mol.get_system()
        if mol.config['optimize']['lib'] == 'scipy':
            self.optimizer = mol.config['optimize']['optimizer']
        else:
            self.optimizer = mol.config['optimize']['lib']

        # check optimizer
        if self.optimizer not in ['cg', 'bfgs', 'l-bfgs-b', 'newton-cg', 'dlfind']:
            raise ValueError(f'Unknown optimizer {self.optimizer}')

        self.metrics = {
            'nstate': self.nstate,
            'istate': self.istate,
            'jstate': self.jstate,
            'itr': self.itr,
            'de': 0,
            'gap': 0,
            'rmsd_step': 0,
            'max_step': 0,
            'rmsd_grad': 0,
            'max_grad': 0,
            'energy_shift': self.energy_shift,
            'energy_gap': self.energy_gap,
            'target_rmsd_step': self.rmsd_step,
            'target_max_step': self.max_step,
            'target_rmsd_grad': self.rmsd_grad,
            'target_max_grad': self.max_grad,
        }

        dump_log(self.mol, title='PyOQP: Entering Geometry Optimization (%s)' % self.optimizer)

    def optimize(self):
        try:
            sc.optimize.minimize(
                fun=self.opt_func, x0=self.pre_coord, method=self.optimizer, jac=True, tol=1.0e-20,
                options={'maxiter': 9999},
            )
        except StopIteration:
            pass

    def one_step(self, coordinates):
        # user defined energy and gradient calculation function in the subclass
        raise NameError('Optimizer one_step functon is not defined')

    def opt_func(self, coordinates):
        # wrapper function to call one_step then check convergence
        energy, grad = self.one_step(coordinates)
        self.check_convergence()

        return energy, grad

    def check_convergence(self):
        # write convergence to log
        dump_log(self.mol, title='Geometry Optimization Convergence %s' % self.itr, section='opt', info=self.metrics)

        if np.abs(self.metrics['de']) <= self.energy_shift and \
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


class TemplateOpt(Optimizer):
    """
    OQP abstract class for optimization

    """

    def __init__(self, mol):
        super().__init__(mol)

    def one_step(self, coordinates):
        # add iteration
        self.itr += 1

        dump_log(self.mol, title='PyOQP: Geometry Optimization Step %s' % self.itr)

        # update coordinates
        coordinates = coordinates.reshape((self.natom, 3))
        self.mol.update_system(coordinates)

        # do energy and gradient calculation

        return None, None


class ConstrainOpt(Optimizer):
    """
    OQP Constrained optimization class

    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.ref_coord = mol.get_system()
        self.metrics['radius'] = 0
        self.metrics['distance'] = 0
        self.metrics['step_size'] = self.step_size
        self.metrics['step_tol'] = self.step_tol
        self.mass = np.repeat(mol.get_mass(), 3)
        self.tmass = np.sum(self.mass)
        self.init_energy = 0
        self.last_energy = 0
        self.init_xyz = None
        self.last_xyz = None
        self.last_radius = 0
        self.opt_status = 0
        self.message = 'not converged'

    def one_step(self, coordinates):
        # add iteration
        self.itr += 1

        dump_log(self.mol, title='PyOQP: Constrained Geometry Optimization Step %s' % self.itr)

        # compute mass-weighted radius and Cartesian distance
        if self.itr == 1:
            radius = 0
            distance = 0
            f = 1
            df = np.zeros_like(coordinates)

            do_init_scf = True
        else:
            dcoord = (coordinates - self.ref_coord) * (self.mass / self.tmass) ** 0.5
            radius = np.sum(dcoord ** 2) ** 0.5
            distance = np.sum((coordinates - self.ref_coord) ** 2) ** 0.5
            fs = 1
            if radius > 0:
                f = fs * (radius - self.step_size) ** 2 / self.step_size ** 2
                df = fs * 2 * (radius - self.step_size) / radius * dcoord / self.step_size ** 2
            else:
                f = 1
                df = np.zeros_like(coordinates)

            do_init_scf = self.init_scf

        # update coordinates
        self.mol.update_system(coordinates)

        # compute 1e integral
        oqp.library.ints_1e(self.mol)

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
        en = energy  # copy electronic energy

        # apply RELU to constrain
        f = ReLU(f, radius - self.step_size)
        df = ReLU(df, radius - self.step_size)

        # project out grad from df and nf
        df -= np.sum(df * grad) / np.sum(grad ** 2) * grad
        nf = df * self.step_size ** 2

        """
        print(
            'R %8.4f T %8.4f E %14.6f F %14.6f De %14.6f G %14.6f D %14.6f N %14.6f ' % (
                radius, self.step_size, energy, f, f + energy - self.pre_energy,
                np.mean(grad ** 2) ** 0.5, np.mean(df ** 2) ** 0.5, np.mean(nf ** 2) ** 0.5
                )
        )
        """

        # apply constrain
        energy += f
        grad += df

        # evaluate metrics
        de = energy - self.pre_energy
        rmsd_step = np.mean((coordinates - self.pre_coord) ** 2) ** 0.5
        max_step = np.amax(np.abs(coordinates - self.pre_coord))
        rmsd_grad = np.mean(nf ** 2) ** 0.5
        max_grad = np.amax(np.abs(nf))

        self.metrics['itr'] = self.itr
        self.metrics['de'] = de
        self.metrics['rmsd_step'] = rmsd_step
        self.metrics['max_step'] = max_step
        self.metrics['rmsd_grad'] = rmsd_grad
        self.metrics['max_grad'] = max_grad
        self.metrics['radius'] = radius
        self.metrics['distance'] = distance

        # store energy and coordinates
        self.pre_energy = energy
        self.pre_coord = coordinates.copy()
        dump_data(
            self.mol,
            (self.itr, self.atoms, coordinates, energy, de, rmsd_step, max_step, rmsd_grad, max_grad, radius, distance),
            title='CONS_SPHERE',
            fpath=self.mol.log_path,
        )

        # store info for mep
        if self.itr == 1:
            self.init_energy = en
            self.init_xyz = coordinates

        self.last_energy = en
        self.last_xyz = coordinates
        self.last_radius = radius

        return energy, grad

    def check_convergence(self):
        # write convergence to log
        dump_log(
            self.mol, title='Geometry Optimization Convergence %s' % self.itr,
            section='cons_sphere', info=self.metrics
        )

        if np.abs(self.metrics['de']) <= self.energy_shift and \
                self.metrics['rmsd_step'] <= self.rmsd_step and \
                self.metrics['max_step'] <= self.rmsd_step and \
                self.metrics['rmsd_grad'] <= self.rmsd_grad and \
                self.metrics['max_grad'] <= self.max_grad and \
                self.metrics['radius'] > self.step_tol:
            dump_log(self.mol, title='PyOQP: Geometry Optimization Has Converged')
            self.opt_status = 1
            self.message = 'converged'
            raise StopIteration

        else:
            if self.itr == self.maxit:
                dump_log(self.mol,
                         title='PyOQP: Geometry Optimization Has Not Converged. Reached The Maximum Iteration')
                raise StopIteration


class MECIOpt(Optimizer):
    """
    OQP MECI optimization class

    upb or penalty method
    """

    def __init__(self, mol):
        super().__init__(mol)

        self.meci_search = mol.config['optimize']['meci_search']
        self.sigma = mol.config['optimize']['pen_sigma']
        self.alpha = mol.config['optimize']['pen_alpha']
        self.incre = mol.config['optimize']['pen_incre']
        self.weights = mol.config['optimize']['gap_weight']
        self.metrics['meci_search'] = self.meci_search
        self.metrics['sigma'] = self.sigma
        self.metrics['alpha'] = self.alpha
        self.metrics['incre'] = self.incre
        self.x = np.zeros(0)  # record dgv
        self.y = np.zeros(0)  # record cgv

        # check method
        self.method = mol.config['input']['method']
        if self.method == 'hf':
            raise ValueError(f'MECI optimization require method=tdhf, but found method={self.method}')

        # choose meci search method
        func_dict = {
            'penalty': self.penalty,
            'ubp': self.ubp,
            'hybrid': self.hybrid,
        }
        self.work_func = func_dict[self.meci_search]

    def one_step(self, coordinates):
        # check state order:
        if self.jstate <= self.istate:
            raise ValueError(f'MECI state i {self.istate} is equal to or higher than state j {self.jstate}')

        # add iteration
        self.itr += 1

        dump_log(self.mol, title='PyOQP: Geometry Optimization Step %s' % self.itr)

        if self.itr == 1:
            do_init_scf = True
        else:
            do_init_scf = self.init_scf

        # update coordinates
        self.mol.update_system(coordinates)

        # compute 1e integral
        oqp.library.ints_1e(self.mol)

        # compute energy
        energies = self.sp.energy(do_init_scf=do_init_scf)

        # compute gradient
        self.grad.grads = [self.istate, self.jstate]
        grads = self.grad.gradient()
        self.mol.energies = energies
        self.mol.grads = grads

        # compute dftd4
        energies, grads = self.ls.compute(self.mol, grad_list=self.grad.grads)
        self.mol.energies = energies
        self.mol.grads = grads

        f, df = self.work_func(coordinates, energies, grads)

        return f, df

    def hybrid(self, coordinates, energies, grads):
        """
        first use penalty method
        switch to ubp when gap is < self.energy_gap
        """

        if self.itr > 1 and self.metrics['gap'] < self.energy_gap:
            self.meci_search = 'ubp'
            f, df = self.ubp(coordinates, energies, grads)
        else:
            self.meci_search = 'penalty'
            f, df = self.penalty(coordinates, energies, grads)

        return f, df

    def ubp(self, coordinates, energies, grads):
        """
        update branching plane method
        for state i and j,
        j = i + 1

        in the first step
        x_0 = [G(j) - G(i)] / np.sum([G(j) - G(i)] ** 2) ** 2
        y_0 = G(i) - np.sum(G(i) * x_0) * x_0
        y_0 = y_0 / np.sum(y_0 ** 2), so x_0 and y_0 are orthonormal

        in the k-th step
        y_k = [(x_k-1 * x_k) * y_k-1 - (y_k-1 * x_k) * x_k-1] / [(y_k-1 * x_k) ** 2 + (x_k-1 * x_k) ** 2] ** 0.5
        note the sign of x_k-1 and y_k-1 are opposite to the paper

        F = [E(j) + E(i)] * 0.5 + [E(j) - E(i)] ** 2 / sigma
        dF = PE * 0.5 + 2 * [E(j) - E(i)] / sigma * [G(j) - G(i)]
        PE = G(j) + G(i) - np.sum([G(j) + G(i)] * x_k) * x_k - np.sum([G(j) + G(i)] * y_k) * y_k

        in the original paper
        sigma = np.sum([G(j) - G(i)] ** 2) ** 0.5

        in this implementation
        we set sigma = np.sum([G(j) - G(i)] ** 2)

        reference: J. Chem. Theory Comput.2010,6,1538–1545
        """
        # flatten data
        energy_i = energies[self.istate]
        energy_j = energies[self.jstate]

        grad_i = grads[self.istate].reshape(-1)
        grad_j = grads[self.jstate].reshape(-1)

        # compute energy and gradient
        sum_e = energy_j + energy_i
        gap_e = energy_j - energy_i
        sum_g = grad_j + grad_i
        gap_g = grad_j - grad_i

        # compute x, dgv
        alpha = np.sum(gap_g ** 2) ** 0.5
        x = gap_g / alpha

        # compute y, cgv
        if len(self.y) > 0:
            t1 = np.sum(self.x * x) * self.y - np.sum(self.y * x) * self.x
            t2 = np.sum(self.y * x) ** 2 + np.sum(self.x * x) ** 2
            y = t1 / t2 ** 0.5
        else:
            y = sum_g * 0.5
            y -= np.sum(y * x) * x
            y /= np.sum(y ** 2) ** 0.5

        # record x and y
        self.x = x
        self.y = y

        # check orthonormal
        self.metrics['norm'] = np.sum(y * y)
        self.metrics['orth'] = np.sum(x * y)

        # compute function and gradient
        f = sum_e * 0.5 + self.weights * gap_e ** 2 / alpha ** 2
        df_1 = 0.5 * (sum_g - np.sum(sum_g * x) * x - np.sum(sum_g * y) * y)
        df_2 = 2 * gap_e * x
        df = df_1 + self.weights * df_2 / alpha

        # evaluate metrics
        de = f - self.pre_energy
        rmsd_step = np.mean((coordinates - self.pre_coord) ** 2) ** 0.5
        max_step = np.amax(np.abs(coordinates - self.pre_coord))
        rmsd_grad = np.mean(df_2 ** 2) ** 0.5
        max_grad = np.amax(np.abs(df_2))
        rmsd_df_1 = np.mean(df_1 ** 2) ** 0.5
        max_df_1 = np.amax(np.abs(df_1))
        self.metrics['itr'] = self.itr
        self.metrics['de'] = de
        self.metrics['gap'] = gap_e
        self.metrics['rmsd_step'] = rmsd_step
        self.metrics['max_step'] = max_step
        self.metrics['rmsd_grad'] = rmsd_grad
        self.metrics['max_grad'] = max_grad

        # store energy and coordinates
        self.pre_energy = f
        self.pre_coord = coordinates.copy()
        dump_data(
            self.mol,
            (self.itr, self.atoms, coordinates, f, de, gap_e, rmsd_step, max_step, rmsd_grad, max_grad, rmsd_df_1, max_df_1),
            title='MECI',
            fpath=self.mol.log_path,
        )

        return f, df

    def penalty(self, coordinates, energies, grads):
        """
        penalty method
        for state i and j,
        j = i + 1
        F = [E(j) + E(i)] * 0.5 + sigma * P
        P = [E(j) - E(i)] ** 2 / [E(j) - E(i) + alpha]
        dF = [G(j) + G(i)] * 0.5 + sigma * dP
        dP = {[E(j) - E(i)] ** 2 + 2 * alpha * [E(j) - E(i)]} / [E(j) - E(i) + alpha] ** 2 * [G(j) - G(i)]

        in the original paper
        sigma is a growing factor, 3.5
        alpha is a constant factor, 0.02

        in this implementation
        we set a constant sigma to 1
        we replace alpha with np.mean([G(j) - G(i)] ** 2) ** 0.5

        reference: J. Phys. Chem. B 2008,112,405-413
        """
        # flatten data
        energy_i = energies[self.istate]
        energy_j = energies[self.jstate]

        grad_i = grads[self.istate].reshape(-1)
        grad_j = grads[self.jstate].reshape(-1)

        # compute function and gradient
        sum_e = energy_j + energy_i
        gap_e = energy_j - energy_i
        sum_g = grad_j + grad_i
        gap_g = grad_j - grad_i

        # gradient of the energy gap function
        dg = 2 * gap_e * gap_g

        if self.alpha == 0:
            alpha = np.mean(dg ** 2) ** 0.5
        else:
            alpha = self.alpha

        self.sigma *= self.incre ** self.itr

        f = sum_e * 0.5 + self.weights * self.sigma * gap_e ** 2 / (gap_e + alpha)
        df_1 = sum_g * 0.5
        df_2 = (gap_e ** 2 + 2 * alpha * gap_e) / (gap_e + alpha) ** 2 * gap_g
        df_1 = df_1 - np.sum(df_1 * df_2) / np.sum(df_2 ** 2) * df_2
        df = df_1 + self.weights * self.sigma * df_2

        """
        print(
            'E %14.6f  P %14.6f G %14.6f S %14.6f G %14.6f D %14.6f A %14.6f' % (
                sum_e * 0.5,
                self.sigma * gap_e ** 2 / (gap_e + self.alpha),
                gap_e,
                self.sigma,
                np.mean(df_1 ** 2) ** 0.5,
                np.mean(df_2 ** 2) ** 0.5,
                self.alpha,
            )
        )
        """

        # evaluate metrics
        de = f - self.pre_energy
        rmsd_step = np.mean((coordinates - self.pre_coord) ** 2) ** 0.5
        max_step = np.amax(np.abs(coordinates - self.pre_coord))
        rmsd_grad = np.mean(np.array(dg ** 2)) ** 0.5
        max_grad = np.amax(np.abs(dg))
        rmsd_df = np.mean(df ** 2) ** 0.5
        max_df = np.amax(np.abs(df))

        self.metrics['itr'] = self.itr
        self.metrics['sigma'] = self.sigma
        self.metrics['de'] = de
        self.metrics['gap'] = gap_e
        self.metrics['rmsd_step'] = rmsd_step
        self.metrics['max_step'] = max_step
        self.metrics['rmsd_grad'] = rmsd_grad
        self.metrics['max_grad'] = max_grad

        # store energy and coordinates
        self.pre_energy = f
        self.pre_coord = coordinates.copy()
        dump_data(
            self.mol,
            (self.itr, self.atoms, coordinates, f, de, gap_e, rmsd_step, max_step, rmsd_grad, max_grad, rmsd_df, max_df),
            title='MECI',
            fpath=self.mol.log_path,
        )

        return f, df

    def check_convergence(self):
        # write convergence to log
        dump_log(
            self.mol, title='Geometry Optimization Convergence %s' % self.itr,
            section=self.meci_search,
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


class MECPOpt(Optimizer):
    """
    OQP MECP optimization class

    quadratic method
    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mecp_search = 'quad'
        self.weights = mol.config['optimize']['gap_weight']
        self.metrics['mecp_search'] = self.mecp_search

        # check method
        self.method = mol.config['input']['method']
        if self.method == 'hf':
            raise ValueError(f'MECP optimization require method=tdhf, but found method={self.method}')

        # choose meci search method
        func_dict = {
            'quad': self.quad,
        }
        self.work_func = func_dict[self.mecp_search]

    def one_step(self, coordinates):
        # check state multiplicity:
        if self.imult == self.jmult:
            raise ValueError(f'MECP state i multiplicity {self.istate} is equal to state j multiplicity {self.jstate}')
        # add iteration
        self.itr += 1

        dump_log(self.mol, title='PyOQP: Geometry Optimization Step %s' % self.itr)

        if self.itr == 1:
            do_init_scf = True
        else:
            do_init_scf = self.init_scf

        # update coordinates
        self.mol.update_system(coordinates)

        # compute 1e integral
        oqp.library.ints_1e(self.mol)

        # compute reference
        ref_energy = self.sp.reference(do_init_scf=do_init_scf)

        # set multiplicity for state i
        dump_log(self.mol, title='PyOQP: PES 1 Mult = %s Root = %s' % (self.imult, self.istate), section='input')
        self.mol.data.set_tdhf_multiplicity(self.imult)
        energies_1 = self.sp.excitation(ref_energy)

        # compute gradient for state i
        self.grad.grads = [self.istate]
        grads_1 = self.grad.gradient()

        # set multiplicity for state j
        dump_log(self.mol, title='PyOQP: PES 2 Mult = %s Root = %s' % (self.jmult, self.jstate), section='input')
        self.mol.data.set_tdhf_multiplicity(self.jmult)
        energies_2 = self.sp.excitation(ref_energy)

        # compute gradient for state j
        self.grad.grads = [self.jstate]
        grads_2 = self.grad.gradient()

        # compute dftd4
        self.mol.energies = np.concatenate((energies_1, energies_2[1:]))
        self.mol.grads = np.concatenate((grads_1, grads_2[1:]))
        energies, grads = self.ls.compute(self.mol, grad_list=[self.istate, self.jstate + self.nstate])
        self.mol.energies = energies
        self.mol.grads = grads

        f, df = self.work_func(coordinates, energies, grads)

        return f, df

    def quad(self, coordinates, energies, grads):
        """
        quadratic optimization

        f = E1 + (E2 - E1) ** 2
        dg = dE2 - dE1
        df1 = dE1 - dg/norm(dg) * np.sum(dE1 * dg/norm(dg))
        df2 = 2 * (E2 - E1) * dg
        df = df1 + df2
        """

        # flatten data
        energy_i = energies[self.istate]
        energy_j = energies[self.jstate + self.nstate]

        grad_i = grads[self.istate].reshape(-1)
        grad_j = grads[self.jstate + self.nstate].reshape(-1)

        # compute energy and gradient
        gap_e = energy_j - energy_i
        gap_g = grad_j - grad_i

        f = energy_j + self.weights * gap_e ** 2
        df1 = grad_j - gap_g * np.sum(grad_j * gap_g) / np.sum(gap_g ** 2)
        df2 = 2 * gap_e * gap_g
        df = df1 + self.weights * df2

        # evaluate metrics
        de = f - self.pre_energy
        rmsd_step = np.mean((coordinates - self.pre_coord) ** 2) ** 0.5
        max_step = np.amax(np.abs(coordinates - self.pre_coord))
        rmsd_grad = np.mean(np.array(df ** 2)) ** 0.5
        max_grad = np.amax(np.abs(df))
        rmsd_df2 = np.mean(df2 ** 2) ** 0.5
        max_df2 = np.amax(np.abs(df2))
        self.metrics['itr'] = self.itr
        self.metrics['de'] = de
        self.metrics['gap'] = gap_e
        self.metrics['rmsd_step'] = rmsd_step
        self.metrics['max_step'] = max_step
        self.metrics['rmsd_grad'] = rmsd_grad
        self.metrics['max_grad'] = max_grad

        # store energy and coordinates
        self.pre_energy = f
        self.pre_coord = coordinates.copy()
        dump_data(
            self.mol,
            (self.itr, self.atoms, coordinates, f, de, gap_e, rmsd_step, max_step, rmsd_grad, max_grad, rmsd_df2, max_df2),
            title='MECP',
            fpath=self.mol.log_path,
        )

        return f, df

    def check_convergence(self):
        # write convergence to log
        dump_log(
            self.mol, title='Geometry Optimization Convergence %s' % self.itr,
            section=self.mecp_search,
            info=self.metrics
        )

        if np.abs(self.metrics['de']) <= self.energy_shift and \
                np.abs(self.metrics['gap']) <= self.energy_gap and \
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


class StateSpecificOpt(Optimizer):
    """
    OQP state specific optimization class

    """

    def __init__(self, mol):
        super().__init__(mol)

        # reset istate for hf methods
        self.method = mol.config['input']['method']
        if self.method == 'hf':
            self.istate = 0

    def one_step(self, coordinates):
        # add iteration
        self.itr += 1

        dump_log(self.mol, title='PyOQP: Geometry Optimization Step %s' % self.itr)

        if self.itr == 1:
            do_init_scf = True
        else:
            do_init_scf = self.init_scf

        # update coordinates
        self.mol.update_system(coordinates)

        # compute 1e integral
        oqp.library.ints_1e(self.mol)

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
        dump_data(
            self.mol, (self.itr, self.atoms, coordinates, energy, de, rmsd_step, max_step, rmsd_grad, max_grad),
            title='OPTIMIZATION', fpath=self.mol.log_path)

        return energy, grad


class MEP:
    """
    OQP minimum energy path calculation class

    """

    def __init__(self, mol):
        self.mol = mol
        self.istate = mol.config['optimize']['istate']
        self.mep_maxit = mol.config['optimize']['mep_maxit']
        self.atoms = mol.get_atoms()
        self.mo = None
        self.mep_itr = 0
        self.mep_energies = []
        self.optimizer = ConstrainOpt(self.mol)

    def optimize(self):
        # record the first geometry
        dump_log(self.mol, title='PyOQP: Entering Minimum Energy Path Calculations')

        for i in range(self.mep_maxit):
            self.mep_itr += 1
            dump_log(self.mol, title='PyOQP: MEP Step %s' % self.mep_itr)

            # constrained geometry optimization
            self.optimizer.optimize()

            # check mep results
            status = self.check_mep()

            if status > 0:
                break

            # update optimizer
            self.optimizer.itr = 0
            self.optimizer.pre_coord = self.optimizer.last_xyz
            self.optimizer.ref_coord = self.optimizer.last_xyz

    def check_mep(self):
        self.mep_energies.append(self.optimizer.last_energy)

        if self.mep_itr == 1:
            dump_data(self.mol,
                      (0, self.atoms, self.optimizer.init_xyz, self.optimizer.init_energy, self.optimizer.init_energy),
                      title='MEP',
                      fpath=self.mol.log_path)
            de = self.optimizer.init_energy - self.optimizer.last_energy
        else:
            de = self.mep_energies[-2] - self.mep_energies[-1]

        # write output
        results = {
            'itr': self.optimizer.itr,
            'istate': self.istate,
            'status': self.optimizer.message,
            'radius': self.optimizer.last_radius,
            'energy': self.optimizer.last_energy,
            'de': de,
        }

        dump_log(self.mol, title='MEP Convergence %s' % self.mep_itr, section='mep', info=results)
        dump_data(
            self.mol,
            (self.mep_itr, self.atoms, self.optimizer.last_xyz, self.optimizer.last_energy, de),
            title='MEP',
            fpath=self.mol.log_path,
        )

        # not converged
        if self.optimizer.opt_status == 0:
            dump_log(self.mol, title='PyOQP: MEP Stopped Due to Constrained Geometry Optimization Has Not Converged')
            return 1

        # energy increase
        if de < 0:
            dump_log(self.mol, title='PyOQP: MEP Stopped Due to Energy Are Not Decreased')

            return 2

        # normal termination
        if 0 < de < 1e-6:
            dump_log(self.mol, title='PyOQP: MEP Has Converged to Minimum')

            return 3

        # reach max iteration
        if self.mep_itr == self.mep_maxit:
            dump_log(self.mol, title='PyOQP: MEP Stopped At The Max Step')

            return 4

        return 0


def ReLU(x, y):
    if y > 0:
        return x
    else:
        return 0

class QMMMOpt(Optimizer):
    """
    OQP QM/MM optimization class

    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.metrics['radius'] = 0
        self.metrics['distance'] = 0
        self.metrics['step_size'] = self.step_size
        self.metrics['step_tol'] = self.step_tol
        self.ref_coord,self.mass,self.natom,self.atoms = qmmm.openmm_info()
        self.pre_coord = self.ref_coord
        self.tmass = np.sum(self.mass)
        self.init_energy = 0
        self.last_energy = 0
        self.init_xyz = None
        self.last_xyz = None
        self.last_radius = 0
        self.opt_status = 0
        self.message = 'not converged'

        if not self.mol.config['input']['qmmm_flag']:
           exit(f"QM/MM optimizations require an active qmmm_flag")

    def one_step(self, coordinates):
        # add iteration
        self.itr += 1

        dump_log(self.mol, title='PyOQP: QM/MM Geometry Optimization Step %s' % self.itr)


        # update coordinates
        qmmm.openmm_update_system(coordinates)
        current_xyz=coordinates.reshape((-1,3))
#       for i in qmmm.qm_atoms:
#           coordinates_qm=np.append(coordinates_qm,current_xyz[i][0])
#           coordinates_qm=np.append(coordinates_qm,current_xyz[i][1])
#           coordinates_qm=np.append(coordinates_qm,current_xyz[i][2])

        num_atoms, x, y, z, q, mass = qmmm.openmm_system()
        coordinates_qm=np.array(x + y + z).reshape((3, num_atoms)).T.reshape(-1)

        self.mol.update_system(coordinates_qm)

        # compute 1e integral
        if self.itr == 1:
            do_init_scf = True
        else:
            do_init_scf = self.init_scf
            oqp.library.ints_1e(self.mol)

        # compute energy
        energies = self.sp.energy(do_init_scf=do_init_scf)

        # compute QM/MM gradient
        current_xyz = self.mol.get_system().reshape((-1, 3))
        gradient_qm,gradient_mm=qmmm.openmm_gradient(current_xyz,self.mol.data["OQP::partial_charges"])
        self.mol.data["OQP::mm_gradient"]=np.transpose(gradient_qm).tolist()

        self.grad.grads = [self.istate]
        grads = self.grad.gradient()
        self.mol.energies = energies
        self.mol.grads = grads

        qmmm.gradient_qmmm=qmmm.form_gradient_qmmm(grads,gradient_mm)

        # flatten data
        energy = energies[self.istate]
        grad = qmmm.gradient_qmmm.reshape(-1)

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
        dump_data((self.itr, self.atoms, coordinates, energy, de, rmsd_step, max_step, rmsd_grad, max_grad),
                  title='QM/MM OPTIMIZATION', fpath=self.mol.log_path)

        return energy, grad

    def check_convergence(self):
        # write convergence to log
        dump_log(
            self.mol, title='Geometry Optimization Convergence %s' % self.itr,
            section='QM/MM', info=self.metrics
        )

        if np.abs(self.metrics['de']) <= self.energy_shift and \
                self.metrics['rmsd_step'] <= self.rmsd_step and \
                self.metrics['max_step'] <= self.rmsd_step and \
                self.metrics['rmsd_grad'] <= self.rmsd_grad and \
                self.metrics['max_grad'] <= self.max_grad and \
                self.metrics['radius'] > self.step_tol:
            dump_log(self.mol, title='PyOQP: Geometry Optimization Has Converged')
            self.opt_status = 1
            self.message = 'converged'
            raise StopIteration

        else:
            if self.itr == self.maxit:
                dump_log(self.mol,
                         title='PyOQP: Geometry Optimization Has Not Converged. Reached The Maximum Iteration')
                raise StopIteration



