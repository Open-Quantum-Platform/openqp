"""OQP single point class"""
import copy
import oqp
import numpy as np
import scipy as sc
from oqp.library.baeka import BaekAState
from oqp.library.single_point import SinglePoint, Gradient, LastStep
from oqp.utils.file_utils import dump_log, dump_data
from oqp.utils.state_labels import is_mrsf, public_state_label
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
        if self.optimizer not in ['cg', 'bfgs', 'l-bfgs-b', 'newton-cg', 'geometric', 'oqp']:
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
                self.metrics['max_step'] <= self.max_step and \
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
                self.metrics['max_step'] <= self.max_step and \
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
        # ``auto`` is orchestrated by the native OQP backend.  Other backends
        # retain their established two-state penalty objective rather than
        # failing while constructing the shared MECI objective.
        if self.meci_search == 'auto':
            self.meci_search = 'penalty'
        self.sigma = mol.config['optimize']['pen_sigma']
        self.alpha = mol.config['optimize']['pen_alpha']
        self.incre = mol.config['optimize']['pen_incre']
        self.weights = mol.config['optimize']['gap_weight']
        self.mu = mol.config['optimize']['gap_sigma']
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
            'auglag': self.auglag,
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

        # evaluate metrics.  MECI convergence must use the full
        # crossing-objective gradient, not only the energy-gap term;
        # individual state gradients are finite at a CI.
        de = f - self.pre_energy
        rmsd_step = np.mean((coordinates - self.pre_coord) ** 2) ** 0.5
        max_step = np.amax(np.abs(coordinates - self.pre_coord))
        rmsd_grad = np.mean(df ** 2) ** 0.5
        max_grad = np.amax(np.abs(df))
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

    def auglag(self, coordinates, energies, grads):
        """
        augmented Lagrangian method
        for state i and j,
        j = i + 1

        The branching plane (x, y) is updated exactly as in ubp, so no
        derivative coupling is required.  On top of that projection a
        least-squares multiplier

            lam = -[G(j) + G(i)] * 0.5 . [G(j) - G(i)] / |G(j) - G(i)| ** 2

        removes the mean-gradient component along the gap direction, and the
        remaining gap term is scaled to match:

            mu = 2 * gap_sigma / |G(j) - G(i)|
            F = [E(j) + E(i)] * 0.5 + lam * [E(j) - E(i)]
                + 0.5 * mu * [E(j) - E(i)] ** 2
            dF = PE + mu * [E(j) - E(i)] * [G(j) - G(i)]

        where PE is the branching-plane-projected mean gradient.  The two
        contributions of dF are orthogonal, so both must vanish separately and
        the stationary point is a genuine intersection.  Unlike ubp, the
        reported objective is an augmented Lagrangian value rather than a
        ratio that diverges as the gradient difference becomes small.
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
        mean_g = sum_g * 0.5
        mu = 2.0 * self.mu / alpha
        lagrange = -np.sum(mean_g * gap_g) / alpha ** 2
        df_1 = mean_g - np.sum(mean_g * x) * x - np.sum(mean_g * y) * y
        df_2 = mu * gap_e * gap_g
        df = df_1 + df_2
        f = sum_e * 0.5 + lagrange * gap_e + 0.5 * mu * gap_e ** 2

        # evaluate metrics
        de = f - self.pre_energy
        rmsd_step = np.mean((coordinates - self.pre_coord) ** 2) ** 0.5
        max_step = np.amax(np.abs(coordinates - self.pre_coord))
        rmsd_grad = np.mean(df ** 2) ** 0.5
        max_grad = np.amax(np.abs(df))
        rmsd_df_2 = np.mean(df_2 ** 2) ** 0.5
        max_df_2 = np.amax(np.abs(df_2))
        self.metrics['itr'] = self.itr
        self.metrics['de'] = de
        self.metrics['gap'] = gap_e
        self.metrics['rmsd_step'] = rmsd_step
        self.metrics['max_step'] = max_step
        self.metrics['rmsd_grad'] = rmsd_grad
        self.metrics['max_grad'] = max_grad
        self.metrics['lagrange'] = lagrange

        # store energy and coordinates
        self.pre_energy = f
        self.pre_coord = coordinates.copy()
        dump_data(
            self.mol,
            (self.itr, self.atoms, coordinates, f, de, gap_e, rmsd_step, max_step, rmsd_grad, max_grad, rmsd_df_2, max_df_2),
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

        if self.alpha == 0:
            # Use the coupling-gradient scale described for the penalty
            # method.  Do not scale alpha by the energy gap; otherwise alpha
            # collapses to zero near the crossing and the penalty objective
            # becomes ill-conditioned.
            alpha = np.mean(gap_g ** 2) ** 0.5
        else:
            alpha = self.alpha

        self.sigma *= self.incre

        f = sum_e * 0.5 + self.weights * self.sigma * gap_e ** 2 / (gap_e + alpha)
        df_1 = sum_g * 0.5
        df_2 = (gap_e ** 2 + 2 * alpha * gap_e) / (gap_e + alpha) ** 2 * gap_g
        # True penalty-function gradient.  The average-gradient component
        # must not be projected out when a scalar optimizer such as scipy or
        # geomeTRIC is minimizing the penalty objective.
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

        # evaluate metrics.  MECI convergence uses the full penalty-objective
        # gradient.  Keep the penalty contribution separately for diagnostics.
        de = f - self.pre_energy
        rmsd_step = np.mean((coordinates - self.pre_coord) ** 2) ** 0.5
        max_step = np.amax(np.abs(coordinates - self.pre_coord))
        rmsd_grad = np.mean(df ** 2) ** 0.5
        max_grad = np.amax(np.abs(df))
        rmsd_df = np.mean(df_2 ** 2) ** 0.5
        max_df = np.amax(np.abs(df_2))

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
                self.metrics['max_step'] <= self.max_step and \
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

    The two states differ in spin multiplicity, so the derivative coupling
    vanishes by spin symmetry and the branching space is the single
    gradient-difference direction.  The crossing algorithms that converge the
    energy gap exactly are therefore available here without any coupling
    calculation, which is why the fixed-weight quadratic penalty is no longer
    the default.

    auglag   augmented Lagrangian (default).  A least-squares multiplier
             absorbs the mean-gradient component along the branching
             direction, leaving a gap term and a projected mean gradient that
             are orthogonal, so a vanishing total gradient forces both to
             vanish separately and the stationary point is the true crossing.
    ubp      the same objective; at pen_sigma = 1 the multiplier form is the
             Bearpark gradient projection.
    penalty  Levine-Martinez smooth penalty, matching the MECI path.
    baeka    adaptive additive-sigma penalty, driven by the shared BaekA
             kernel.
    quad     legacy fixed-weight quadratic penalty.  Its stationary condition
             balances the mean gradient against the gap term, and because the
             two are not orthogonal they cancel at a residual gap of order
             1/gap_weight.  The energy_gap criterion is then unreachable in
             practice.  Retained only to reproduce earlier runs.
    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mecp_search = mol.config['optimize']['mecp_search']
        self.weights = mol.config['optimize']['gap_weight']
        self.sigma = mol.config['optimize']['pen_sigma']
        self.incre = mol.config['optimize']['pen_incre']
        alpha = mol.config['optimize']['pen_alpha']
        # zero is the historical sentinel for "use the published value"
        self.alpha = 0.02 if alpha == 0 else alpha
        # augmented-Lagrangian multiplier and penalty parameter
        self.lagrange = 0.0
        self.mu = mol.config['optimize']['gap_sigma']
        self.metrics['mecp_search'] = self.mecp_search

        # check method
        self.method = mol.config['input']['method']
        if self.method == 'hf':
            raise ValueError(f'MECP optimization require method=tdhf, but found method={self.method}')

        # choose mecp search method
        func_dict = {
            'auglag': self.auglag,
            # the multiplier form reduces to the Bearpark gradient projection
            # at pen_sigma = 1, so ubp names the same objective here
            'ubp': self.auglag,
            'penalty': self.penalty,
            'baeka': self.baeka,
            'quad': self.quad,
        }
        if self.mecp_search not in func_dict:
            raise ValueError(
                f'Unknown MECP search algorithm {self.mecp_search}, '
                f'use one of {", ".join(sorted(func_dict))}'
            )
        self.work_func = func_dict[self.mecp_search]

        if self.mecp_search == 'baeka':
            self.baeka_state = BaekAState(
                sigma=self.sigma,
                alpha=self.alpha,
                delta_beta=mol.config['optimize']['pen_delta'],
                beta_schedule=mol.config['optimize']['pen_jump'],
                gap_tol=self.energy_gap,
                tol_f=self.energy_shift,
                tol_g=self.rmsd_grad,
                gap_weight=self.weights,
            )

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

        # compute reference
        ref_energy = self.sp.reference(do_init_scf=do_init_scf)

        td_config = self.mol.config['tdhf']
        saved_mult = td_config.get('multiplicity', self.imult)

        def select_multiplicity(mult):
            # Synchronize the Python configuration used by log/state-label
            # helpers with the native target actually being evaluated.
            td_config['multiplicity'] = mult
            self.mol.data.set_tdhf_multiplicity(mult)

        try:
            # set multiplicity for state i
            pes1 = (public_state_label(self.mol.config, self.istate, self.imult)
                    if is_mrsf(self.mol.config)
                    else 'multiplicity %s, root %s' % (self.imult, self.istate))
            dump_log(self.mol, title='PyOQP: PES 1 = %s' % pes1, section='input')
            select_multiplicity(self.imult)
            energies_1 = self.sp.excitation(ref_energy)

            # compute gradient for state i
            self.grad.grads = [self.istate]
            grads_1 = self.grad.gradient()

            # set multiplicity for state j
            pes2 = (public_state_label(self.mol.config, self.jstate, self.jmult)
                    if is_mrsf(self.mol.config)
                    else 'multiplicity %s, root %s' % (self.jmult, self.jstate))
            dump_log(self.mol, title='PyOQP: PES 2 = %s' % pes2, section='input')
            select_multiplicity(self.jmult)
            energies_2 = self.sp.excitation(ref_energy)

            # compute gradient for state j
            self.grad.grads = [self.jstate]
            grads_2 = self.grad.gradient()
        finally:
            select_multiplicity(saved_mult)

        # compute dftd4
        self.mol.energies = np.concatenate((energies_1, energies_2[1:]))
        self.mol.grads = np.concatenate((grads_1, grads_2[1:]))
        energies, grads = self.ls.compute(self.mol, grad_list=[self.istate, self.jstate + self.nstate])
        self.mol.energies = energies
        self.mol.grads = grads

        f, df = self.work_func(coordinates, energies, grads)

        return f, df

    def state_pair(self, energies, grads):
        """Return the two crossing states as (E_i, E_j, G_i, G_j).

        State j lives in the second multiplicity block, which one_step appends
        after the nstate energies of the first.
        """
        energy_i = energies[self.istate]
        energy_j = energies[self.jstate + self.nstate]

        grad_i = grads[self.istate].reshape(-1)
        grad_j = grads[self.jstate + self.nstate].reshape(-1)

        return energy_i, energy_j, grad_i, grad_j

    def ordered_pair(self, energies, grads):
        """Return the crossing states sorted by energy, lower state first.

        The smooth penalties are defined for a non-negative gap, while the two
        MECP states come from independent multiplicity solves and may arrive in
        either order.
        """
        energy_i, energy_j, grad_i, grad_j = self.state_pair(energies, grads)

        if energy_j >= energy_i:
            return energy_i, energy_j, grad_i, grad_j

        return energy_j, energy_i, grad_j, grad_i

    def record(self, coordinates, f, df, gap_e, gap_term):
        """Store metrics, append the status row, and return (f, df)."""
        de = f - self.pre_energy
        rmsd_step = np.mean((coordinates - self.pre_coord) ** 2) ** 0.5
        max_step = np.amax(np.abs(coordinates - self.pre_coord))
        rmsd_grad = np.mean(df ** 2) ** 0.5
        max_grad = np.amax(np.abs(df))
        rmsd_gap = np.mean(gap_term ** 2) ** 0.5
        max_gap = np.amax(np.abs(gap_term))
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
            (self.itr, self.atoms, coordinates, f, de, gap_e, rmsd_step, max_step, rmsd_grad, max_grad, rmsd_gap, max_gap),
            title='MECP',
            fpath=self.mol.log_path,
        )

        return f, df

    def auglag(self, coordinates, energies, grads):
        """
        augmented Lagrangian optimization

        x   = (G2 - G1) / |G2 - G1|
        lam = -(0.5 * (G1 + G2) . (G2 - G1)) / |G2 - G1| ** 2
        mu  = 2 * sigma / |G2 - G1|
        F   = 0.5 * (E1 + E2) + lam * (E2 - E1) + 0.5 * mu * (E2 - E1) ** 2
        dF  = [I - x x^T] 0.5 * (G1 + G2) + mu * (E2 - E1) * (G2 - G1)

        A plain penalty stalls because 0.5 * (G1 + G2) has a component along
        the branching direction that only an infinite weight can overcome.  The
        least-squares multiplier lam removes exactly that component, so the two
        remaining contributions are orthogonal and dF = 0 forces E2 - E1 = 0
        and a vanishing projected mean gradient separately, for any positive
        mu.  Scaling mu by 1/|G2 - G1| keeps the gap term dimensionally
        matched to the projected gradient; at pen_sigma = 1 it reduces to the
        Bearpark gradient projection 2 (E2 - E1) x.

        reference: Chem. Phys. Lett. 1994, 223, 269-274 (gradient projection);
        Nocedal and Wright, Numerical Optimization, ch. 17 (multiplier form)
        """

        energy_i, energy_j, grad_i, grad_j = self.state_pair(energies, grads)

        gap_e = energy_j - energy_i
        gap_g = grad_j - grad_i
        mean_g = (grad_i + grad_j) * 0.5

        norm = np.sum(gap_g ** 2) ** 0.5
        if norm == 0:
            raise ValueError(
                'MECP branching direction is undefined: the two states have '
                'identical gradients'
            )
        x = gap_g / norm
        mu = 2.0 * self.mu / norm
        lagrange = -np.sum(mean_g * gap_g) / norm ** 2

        gap_term = mu * gap_e * gap_g
        df = (mean_g - np.sum(mean_g * x) * x) + gap_term
        f = (energy_i + energy_j) * 0.5 + lagrange * gap_e + 0.5 * mu * gap_e ** 2

        self.lagrange = lagrange
        self.metrics['lagrange'] = lagrange

        return self.record(coordinates, f, df, gap_e, gap_term)

    def penalty(self, coordinates, energies, grads):
        """
        smooth penalty optimization

        F = 0.5 * (E1 + E2) + w * sigma * D ** 2 / (D + alpha)
        dF = 0.5 * (G1 + G2)
             + w * sigma * (D ** 2 + 2 * alpha * D) / (D + alpha) ** 2 * dG
        where D is the non-negative gap and dG the matching gradient difference.

        reference: J. Phys. Chem. B 2008, 112, 405-413
        """

        energy_lo, energy_up, grad_lo, grad_up = self.ordered_pair(energies, grads)
        energy_i, energy_j, grad_i, grad_j = self.state_pair(energies, grads)

        gap_e = energy_j - energy_i
        span = energy_up - energy_lo
        gap_g = grad_up - grad_lo
        mean_g = (grad_i + grad_j) * 0.5

        scale = self.weights * self.sigma
        slope = (span ** 2 + 2 * self.alpha * span) / (span + self.alpha) ** 2
        f = (energy_i + energy_j) * 0.5 + scale * span ** 2 / (span + self.alpha)
        gap_term = scale * slope * gap_g
        df = mean_g + gap_term

        self.sigma *= self.incre
        self.metrics['sigma'] = self.sigma

        return self.record(coordinates, f, df, gap_e, gap_term)

    def baeka(self, coordinates, energies, grads):
        """
        adaptive additive-sigma penalty optimization

        Reuses the shared BaekA kernel, which needs only ascending energies and
        their gradients and makes no reference to spin, so the spin-crossing
        pair is a valid two-state input.

        reference: Y. S. Baek, S. Lee, M. Filatov, and C. H. Choi,
        J. Phys. Chem. A 125, 1994-2006 (2021)
        """

        energy_lo, energy_up, grad_lo, grad_up = self.ordered_pair(energies, grads)
        energy_i, energy_j, _, _ = self.state_pair(energies, grads)

        gap_e = energy_j - energy_i
        evaluation = self.baeka_state.evaluate(
            [energy_lo, energy_up], np.array([grad_lo, grad_up])
        )
        data = evaluation.data

        f = data.objective
        df = data.gradient
        gap_term = data.effective_sigma * data.penalty_gradient
        self.metrics['sigma'] = evaluation.next_sigma

        return self.record(coordinates, f, df, gap_e, gap_term)

    def quad(self, coordinates, energies, grads):
        """
        quadratic penalty optimization

        F = 0.5 * (E1 + E2) + w * (E2 - E1) ** 2
        dF = 0.5 * (G1 + G2) + 2 * w * (E2 - E1) * (G2 - G1)

        Legacy objective.  The stationary condition cancels the mean gradient
        against the gap term, leaving a residual gap of order 1/w, so this
        cannot satisfy a tight energy_gap criterion.  Kept to reproduce runs
        made before the converging objectives were added.
        """

        energy_i, energy_j, grad_i, grad_j = self.state_pair(energies, grads)

        # compute energy and gradient
        gap_e = energy_j - energy_i
        gap_g = grad_j - grad_i

        f = (energy_i + energy_j) * 0.5 + self.weights * gap_e ** 2
        df1 = (grad_i + grad_j) * 0.5
        df2 = 2 * gap_e * gap_g
        df = df1 + self.weights * df2

        return self.record(coordinates, f, df, gap_e, df2)

    def check_convergence(self):
        # write convergence to log
        dump_log(
            self.mol, title='Geometry Optimization Convergence %s' % self.itr,
            section='mecp',
            info=self.metrics
        )

        if np.abs(self.metrics['de']) <= self.energy_shift and \
                np.abs(self.metrics['gap']) <= self.energy_gap and \
                self.metrics['rmsd_step'] <= self.rmsd_step and \
                self.metrics['max_step'] <= self.max_step and \
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

        target_converger = self.mol.config['scf']['converger_type']
        try:
            self.grad.grads = [self.istate]
            if self.method == 'dftb':
                # OpenQP-DFTB's analytic-gradient C call also returns the
                # complete, matching state-energy spectrum.  Use that single
                # solve for optimization so an independent energy-only ROKS
                # cycle cannot select a different open-shell branch.
                grads = self.grad.gradient()
                energies = np.asarray(self.mol.energies, dtype=float)
            else:
                # If SCF falls back to alternative_scf (e.g. TRAH), keep that
                # recovered reference through the matching excited-state
                # gradient/Z-vector evaluation for this geometry only.
                energies = self.sp.energy(
                    do_init_scf=do_init_scf,
                    restore_scf_converger=False,
                )
                grads = self.grad.gradient()
            self.mol.energies = energies
            self.mol.grads = grads

            # compute dftd4
            energies, grads = self.ls.compute(self.mol, grad_list=self.grad.grads)
            self.mol.energies = energies
            self.mol.grads = grads
        finally:
            # The fallback converger is only for the current geometry-cycle
            # recovery. The next geometry cycle must start from the original
            # user-requested SCF optimizer.
            self.mol.data.set_scf_converger_type(target_converger)

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

#        if not self.mol.config['input']['qmmm_flag']:
#           exit(f"QM/MM optimizations require an active qmmm_flag")

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
