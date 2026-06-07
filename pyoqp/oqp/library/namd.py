"""
Nonadiabatic molecular dynamics (NAMD) — Tully fewest-switches surface hopping
(FSSH) on MRSF-TDDFT adiabatic states.

Design
------
The numerically critical surface-hopping physics lives in Fortran
(``source/modules/namd.F90`` -> C entry ``mrsf_namd_hop``); this Python driver
only *sequences* the existing electronic-structure drivers and integrates the
nuclei with velocity Verlet:

  per step:  position update (Verlet)
             -> SinglePoint.reference()  (SCF)
             -> BasisOverlap.overlap()   (MO overlap vs previous step, phase aligned)
             -> SinglePoint.excitation() (MRSF energies + response vectors)
             -> Gradient (active state)
             -> velocity 2nd half-kick
             -> NACME.nacme()            (state overlap S = <i(t-dt)|j(t)>)
             -> oqp.mrsf_namd_hop()      (TDC, RK4 amplitude propagation, EDC,
                                          trivial-crossing follow, FSSH hop +
                                          isotropic velocity rescaling)
             -> on hop: recompute gradient on the new active surface
             -> output / restart

Units: coordinates in bohr, masses in electron masses, energies in Hartree,
velocities in bohr / atomic-time.  Time step is given in fs and converted.

This is the gas-phase (all-QM) path; QM/MM and PBC are layered on later.
"""

import os
import copy
import numpy as np

import oqp
from oqp.library.single_point import SinglePoint, Gradient, LastStep, BasisOverlap, NACME
from oqp.utils.file_utils import dump_log

# 1 fs in atomic units of time
FS_TO_AU = 41.341374575751
# Boltzmann constant in Hartree / Kelvin
KB_HARTREE = 3.166811563e-6
# 1 atomic mass unit (Dalton) in electron masses
AMU_TO_AU = 1822.888486209

# OQP::namd_params packed-scalar indices (0-based here; 1-based in the contract)
_P_DT_FS = 0
_P_NSUB = 1
_P_THRSHE = 2
_P_RAND = 3
_P_ACTIVE = 4
_P_DECO = 5
_P_EDC_C = 6
_P_TDC = 7
_P_TRIV = 8
_P_TRIV_THR = 9
_P_HOPPED = 10
_P_TARGET = 11
_NPARAMS = 16


class NAMD:
    """Driver for FSSH nonadiabatic molecular dynamics."""

    def __init__(self, mol):
        self.mol = mol
        cfg = mol.config
        md = cfg['md']

        self.nstep = int(md['nstep'])
        self.dt_fs = float(md['dt'])
        self.dt = self.dt_fs * FS_TO_AU
        self.active = int(md['active'])            # 1-based excited-state index
        self.nstate = int(cfg['tdhf']['nstate'])
        self.substep = int(md['substep'])
        self.decoherence = 1 if str(md['decoherence']).lower() in ('edc', 'on', 'true', '1') else 0
        self.edc_c = float(md['edc_c'])
        self.thrshe = float(md['thrshe'])
        self.tdc_scheme = 1 if str(md['tdc']).lower() == 'npi' else 0
        self.trivial = 1 if str(md['trivial']).lower() in ('true', '1', 'on', 'yes') else 0
        self.trivial_thresh = float(md['trivial_thresh'])
        self.init_temp = float(md['init_temp'])
        self.seed = int(md['seed'])
        self.velocity_source = str(md['velocity'])

        self.natom = mol.data['natom']
        # get_mass() returns atomic masses in amu; the integrator works in
        # atomic units, so convert to electron masses.
        self.mass = mol.get_mass() * AMU_TO_AU     # (natom,) electron masses
        self.rng = np.random.default_rng(self.seed)

        # electronic amplitudes (complex), one per excited state
        self.coef = np.zeros(self.nstate, dtype=complex)
        self.coef[self.active - 1] = 1.0 + 0.0j

        # velocities (natom, 3) in atomic units
        self.vel = self._init_velocities()

        # previous-step payload for the overlap (back_door carry)
        self.prev_xyz = None
        self.prev_data = None

        # force the per-step electronic pipeline to use the in-memory previous
        # step rather than recomputing it
        cfg['properties']['back_door'] = True
        # NACME needs a dt; reuse the MD step (atomic units) for the TDC scale
        cfg['nac']['dt'] = self.dt

    # ------------------------------------------------------------------ #
    # initialisation
    # ------------------------------------------------------------------ #
    def _init_velocities(self):
        src = self.velocity_source.lower()
        if src in ('zero', 'none', '0'):
            return np.zeros((self.natom, 3))
        if src in ('maxwell', 'boltzmann', 'random'):
            sigma = np.sqrt(KB_HARTREE * self.init_temp / self.mass)  # (natom,)
            v = self.rng.normal(0.0, 1.0, size=(self.natom, 3)) * sigma[:, None]
            return self._remove_com_motion(v)
        # otherwise treat as a file path: "vx vy vz" per atom (atomic units)
        if os.path.isfile(self.velocity_source):
            v = np.loadtxt(self.velocity_source).reshape((self.natom, 3))
            return self._remove_com_motion(v)
        raise ValueError(f"[md] velocity='{self.velocity_source}' is not zero/maxwell or a readable file")

    def _remove_com_motion(self, v):
        p = (self.mass[:, None] * v).sum(axis=0)        # total momentum
        v = v - p / self.mass.sum()
        return v

    # ------------------------------------------------------------------ #
    # electronic structure for one geometry
    # ------------------------------------------------------------------ #
    def _electronic(self, with_overlap):
        """Run SCF + (optional overlap) + MRSF excitation at the current geometry."""
        mol = self.mol
        sp = SinglePoint(mol)
        ref_energy = sp.reference()
        if with_overlap:
            mol.back_door = (self.prev_xyz, self.prev_data)
            BasisOverlap(mol).overlap()
        sp.excitation(ref_energy)
        LastStep(mol).compute(mol)

    def _active_gradient(self):
        """Compute and return the gradient (natom,3) on the current active state."""
        mol = self.mol
        mol.config['properties']['grad'] = [self.active]
        Gradient(mol).gradient()
        return np.array(mol.grads[self.active]).reshape((self.natom, 3))

    def _state_overlap(self):
        """Compute the phase-corrected state overlap S(i,j)=<i(t-dt)|j(t)>."""
        NACME(self.mol).nacme()
        return np.array(self.mol.data["OQP::td_states_overlap"])

    # ------------------------------------------------------------------ #
    # Fortran FSSH hop
    # ------------------------------------------------------------------ #
    def _hop(self):
        """Call the Fortran FSSH kernel; updates amplitudes, velocities, active state."""
        mol = self.mol
        n = self.nstate

        # amplitudes: flat 1-D, interleaved [re1, im1, re2, im2, ...]
        coef_io = np.zeros(2 * n)
        coef_io[0::2] = self.coef.real
        coef_io[1::2] = self.coef.imag
        mol.data["OQP::namd_coef"] = coef_io
        # velocities: flat 1-D, [vx1, vy1, vz1, vx2, ...] (atom-major)
        mol.data["OQP::namd_velocity"] = self.vel.reshape(-1).copy()

        params = np.zeros(_NPARAMS)
        params[_P_DT_FS] = self.dt_fs
        params[_P_NSUB] = float(self.substep)
        params[_P_THRSHE] = self.thrshe
        params[_P_RAND] = float(self.rng.random())
        params[_P_ACTIVE] = float(self.active)
        params[_P_DECO] = float(self.decoherence)
        params[_P_EDC_C] = self.edc_c
        params[_P_TDC] = float(self.tdc_scheme)
        params[_P_TRIV] = float(self.trivial)
        params[_P_TRIV_THR] = self.trivial_thresh
        mol.data["OQP::namd_params"] = params

        oqp.mrsf_namd_hop(mol)

        # read back
        coef_io = np.array(mol.data["OQP::namd_coef"]).reshape(-1)
        self.coef = coef_io[0::2] + 1j * coef_io[1::2]
        self.vel = np.array(mol.data["OQP::namd_velocity"]).reshape((self.natom, 3)).copy()
        params = np.array(mol.data["OQP::namd_params"])
        new_active = int(round(params[_P_ACTIVE]))
        hopped = int(round(params[_P_HOPPED])) == 1
        return new_active, hopped

    # ------------------------------------------------------------------ #
    # main loop
    # ------------------------------------------------------------------ #
    def run(self):
        mol = self.mol
        dump_log(mol, title='PyOQP: Tully FSSH Nonadiabatic Molecular Dynamics')

        r = mol.get_system().reshape((self.natom, 3))   # bohr

        # initial electronic structure + force on the active state
        self._electronic(with_overlap=False)
        accel = -self._active_gradient() / self.mass[:, None]
        self._record_previous(r)
        self._log_step(0, r)

        for istep in range(1, self.nstep + 1):
            # velocity-Verlet position update
            r = r + self.vel * self.dt + 0.5 * accel * self.dt ** 2
            mol.update_system(r.reshape(-1))

            # electronic structure at the new geometry (with overlap vs previous)
            self._electronic(with_overlap=True)
            accel_new = -self._active_gradient() / self.mass[:, None]

            # velocity-Verlet velocity update
            self.vel = self.vel + 0.5 * (accel + accel_new) * self.dt

            # state overlap (couplings) and FSSH hop
            self._state_overlap()
            new_active, hopped = self._hop()

            if hopped:
                self.active = new_active
                # force for the next step is on the new active surface
                accel_new = -self._active_gradient() / self.mass[:, None]

            accel = accel_new
            self._record_previous(r)
            self._log_step(istep, r, hopped=hopped)

        dump_log(mol, title='PyOQP: NAMD trajectory complete')

    # ------------------------------------------------------------------ #
    # helpers
    # ------------------------------------------------------------------ #
    def _record_previous(self, r):
        self.prev_xyz = copy.deepcopy(r.reshape(-1))
        self.prev_data = copy.deepcopy(self.mol.get_data())

    def _log_step(self, istep, r, hopped=False):
        mol = self.mol
        e = np.array(mol.energies)
        ekin = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
        epot = float(e[self.active])
        pops = np.abs(self.coef) ** 2
        dump_log(
            mol,
            title=(f'NAMD step {istep:6d}  t={istep*self.dt_fs:9.3f} fs  '
                   f'active={self.active}  E_tot={ekin+epot:.8f}  '
                   f'E_pot={epot:.8f}  E_kin={ekin:.8f}  '
                   f'hop={hopped}  pop={np.array2string(pops, precision=4)}'),
        )
