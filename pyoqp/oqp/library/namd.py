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
# unit conversions for the QM/MM (OpenMM <-> atomic units) coupling
BOHR_TO_NM = 0.052917721090
NM_TO_BOHR = 1.0 / BOHR_TO_NM
# 1 Hartree/bohr in kJ/mol/nm  (2625.499639 kJ/mol per Ha / 0.0529177 nm per bohr)
HABOHR_TO_KJMOLNM = 2625.499639 / BOHR_TO_NM
KJMOL_TO_HARTREE = 1.0 / 2625.499639

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
_P_NSTATE = 12          # number of states for the hop (0 -> tddft.nstate)
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
    # time-derivative couplings
    # ------------------------------------------------------------------ #
    def _compute_tdc(self, s):
        """Time-derivative coupling matrix from the phase-corrected state
        overlap s(i,j)=<i(t-dt)|j(t)>.

        'fd'  : Hammes-Schiffer/Tully finite difference  (s - s^T)/(2 dt)
        'npi' : norm-preserving interpolation (Meek & Levine, JPCL 5, 2351
                (2014)) in its rigorous matrix form -- the real antisymmetric
                logarithm of the Loewdin-orthonormalised step overlap,
                T = logm(s (s^T s)^{-1/2}) / dt, which reduces to the exact
                two-state identity T*dt = arcsin(s_10) and to the finite
                difference in the weak-coupling limit.
        """
        if self.tdc_scheme == 1:
            from scipy.linalg import logm, sqrtm
            m = s.T @ s
            u = s @ np.linalg.inv(np.real(sqrtm(m)))     # nearest orthogonal (Loewdin)
            t = np.real(logm(u))
            t = 0.5 * (t - t.T)                          # enforce antisymmetry
            return t / self.dt
        return (s - s.T) / (2.0 * self.dt)

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
        params[_P_NSTATE] = float(n)
        mol.data["OQP::namd_params"] = params

        # state overlap + time-derivative couplings (FD or NPI), passed to the
        # Fortran hop as flat row-major (n x n) matrices; absolute state
        # energies via namd_eabs. (Same-spin MRSF path.)
        s = np.array(mol.data["OQP::td_states_overlap"]).reshape((n, n))
        tdc = self._compute_tdc(s)
        mol.data["OQP::namd_tdc"] = tdc.reshape(-1).copy()
        mol.data["OQP::namd_stas"] = s.reshape(-1).copy()
        mol.data["OQP::namd_eabs"] = np.array(mol.data["OQP::td_energies"]).reshape(-1)[:n].copy()

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


def _parse_int_list(spec):
    """Parse '0-2,5,7-8' into [0,1,2,5,7,8]."""
    out = []
    for tok in str(spec).replace(',', ' ').split():
        if '-' in tok:
            a, b = tok.split('-')
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(tok))
    return out


class NAMD_QMMM(NAMD):
    """FSSH NAMD with electrostatic ESPF QM/MM embedding (non-periodic).

    The QM region is the OpenQP Molecule; the MM region + QM/MM coupling are
    handled by OpenMM via the OpenQpQMMM driver.  Per step:
      * sync positions (QM Molecule + OpenMM context),
      * MM electrostatic potential at QM atoms (POTMM),
      * embedded SCF + MRSF excitation (all states),
      * active-state embedded gradient = Gradient + grad_esp_qmmm_excited,
      * ESPF QM charges -> MM forces (forces_mm),
      * full-system velocity Verlet (QM+MM, atomic units),
      * QM-only FSSH hop (rescale only QM velocities, as in GAMESS RESCALV).
    """

    def __init__(self, mol):
        super().__init__(mol)
        import openmm as mm
        import openmm.app as app
        import openmm.unit as u
        from oqp.library.qmmm_driver import OpenQpQMMM
        self._mm = mm
        self._app = app
        self._u = u

        from oqp.library.qmmm_md import _resolve_cutoff

        q = mol.config['qmmm']
        pdb_file = q['pdb_file']
        ff_files = [s for s in str(q['forcefield_files']).replace(',', ' ').split() if s]
        self.qm_atoms = np.array(_parse_int_list(q['qm_atoms']), dtype=int)
        self.cutoff = _resolve_cutoff(str(q['cutoff']).strip())   # NoCutoff | PME | Ewald | ...
        self.periodic = self.cutoff is not app.NoCutoff
        embedding = str(q['embedding']).strip()

        self.pdb = app.PDBFile(pdb_file)
        self.forcefield = app.ForceField(*ff_files)
        self.driver = OpenQpQMMM(
            positions=self.pdb.positions,
            topology=self.pdb.topology,
            forcefield=self.forcefield,
            qm_atoms=self.qm_atoms,
            mol=mol,
            Cutoff=self.cutoff,
            Embedding=embedding,
        )
        self.mm = self.driver.mm_systems

        # full-system state (atomic units)
        self.natom_all = self.pdb.topology.getNumAtoms()
        pos_nm = np.array(self.pdb.positions.value_in_unit(u.nanometer))
        self.r_all = pos_nm * NM_TO_BOHR                       # (natom_all, 3) bohr
        sys0 = self.mm["sys0"]
        self.m_all = np.array([
            sys0.getParticleMass(i).value_in_unit(u.dalton) for i in range(self.natom_all)
        ]) * AMU_TO_AU                                          # electron masses

        # full-system Maxwell-Boltzmann velocities (a.u.), COM removed
        sig = np.sqrt(KB_HARTREE * self.init_temp / self.m_all)
        self.v_all = self.rng.normal(0.0, 1.0, size=(self.natom_all, 3)) * sig[:, None]
        p = (self.m_all[:, None] * self.v_all).sum(axis=0)
        self.v_all -= p / self.m_all.sum()

        # sync the QM Molecule geometry from the pdb QM atoms
        self._sync_positions()
        # QM-region masses for the hop (already set by super from mol.get_mass())
        self.qm_mass = self.mass.copy()

    # ------------------------------------------------------------------ #
    def _sync_positions(self):
        """Push the current full-system positions into OpenMM + the QM Molecule."""
        u = self._u
        pos_q = (self.r_all * BOHR_TO_NM) * u.nanometer
        self.driver.positions = pos_q
        self.mm["sim0"].context.setPositions(pos_q)
        # periodic contexts used by the Ewald QM-QM correction (electrostatic_potential)
        if self.periodic:
            for key in ("simew", "simor"):
                sim = self.mm.get(key)
                if sim is not None:
                    sim.context.setPositions(pos_q)
        # QM Molecule coords (bohr) from the pdb-indexed positions
        self.mol.update_system(self.r_all[self.qm_atoms].reshape(-1))

    def _electronic_qmmm(self, with_overlap):
        """Embedded SCF + MRSF excitation; returns (potmm, potqm)."""
        from oqp.library.qmmm_driver import (
            unpack_lower_tri_single, unpack_lower_tri_multi, pack_lower_tri_single)
        mol = self.mol
        potmm, potqm = self.driver.electrostatic_potential()

        sp = SinglePoint(mol)
        sp._prep_guess()
        nat = mol.data["natom"]
        nbf = mol.data.get_basis()["nbf"]
        mol.data["OQP::POTMM"] = potmm
        mol.data["OQP::POTQM"] = potqm if potqm is not None else np.zeros((nat, nat))
        oqp.espf_op_corr(mol)
        espf = unpack_lower_tri_multi(mol.data["OQP::ESPF_CORR"], nbf, nat)
        hcore = unpack_lower_tri_single(mol.get_hcore(), nbf)
        hcore += np.einsum("ijk,i->jk", espf, potmm)
        mol.set_hcore(pack_lower_tri_single(hcore))
        sp.scf()
        ref = [mol.get_scf_energy()]
        if with_overlap:
            mol.back_door = (self.prev_xyz, self.prev_data)
            BasisOverlap(mol).overlap()
        sp.excitation(ref)
        LastStep(mol).compute(mol)
        return potmm, potqm

    def _qm_gradient(self):
        """Embedded active-state gradient (Hartree/bohr) incl. ESPF force."""
        mol = self.mol
        mol.config['properties']['grad'] = [self.active]
        Gradient(mol).gradient()
        g = np.array(mol.grads[self.active]).reshape(-1, 3)
        oqp.grad_esp_qmmm_excited(mol)
        g = g + np.array(mol.data["OQP::ESPF_GRAD"]).reshape(-1, 3)
        return g

    def _potqm_force(self, pchg, sign=1.0):
        """QM-QM Ewald self-interaction correction force (a.u.), QM atoms only.

        Mirrors electrostatic_potential's energy construction but for forces:
        with the QM atoms carrying their ESPF charges and all MM charges zeroed,
        the (Ewald - direct) force difference is the periodic QM-QM correction
        force.  Returns a (natom_all, 3) array (nonzero only on QM atoms).
        """
        u = self._u
        f = np.zeros((self.natom_all, 3))
        simew = self.mm.get("simew")
        simor = self.mm.get("simor")
        if simew is None or simor is None:
            return f
        nbew = next(x for x in simew.system.getForces() if isinstance(x, self._mm.NonbondedForce))
        nbor = next(x for x in simor.system.getForces() if isinstance(x, self._mm.NonbondedForce))
        # save + set charges: QM -> pchg, everything else -> 0
        saved = []
        for i in range(self.natom_all):
            c_e, s_e, e_e = nbew.getParticleParameters(i)
            c_o, s_o, e_o = nbor.getParticleParameters(i)
            saved.append((c_e, s_e, e_e, c_o, s_o, e_o))
            q = 0.0
            if i in self.qm_atoms:
                q = float(pchg[np.where(self.qm_atoms == i)[0][0]])
            nbew.setParticleParameters(i, q * u.elementary_charge, 0.0, 0.0)
            nbor.setParticleParameters(i, q * u.elementary_charge, 0.0, 0.0)
        nbew.updateParametersInContext(simew.context)
        nbor.updateParametersInContext(simor.context)
        kj = u.kilojoule_per_mole / u.nanometer
        f_ew = np.array(simew.context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(kj))
        f_or = np.array(simor.context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(kj))
        f_corr = sign * (f_ew - f_or) / HABOHR_TO_KJMOLNM      # a.u.
        for k, i in enumerate(self.qm_atoms):
            f[i] = f_corr[i]
        # restore
        for i in range(self.natom_all):
            c_e, s_e, e_e, c_o, s_o, e_o = saved[i]
            nbew.setParticleParameters(i, c_e, s_e, e_e)
            nbor.setParticleParameters(i, c_o, s_o, e_o)
        nbew.updateParametersInContext(simew.context)
        nbor.updateParametersInContext(simor.context)
        return f

    def _total_force(self, potmm):
        """Assemble full-system force (a.u.) and total potential energy (Ha)."""
        mol = self.mol
        u = self._u
        # active-state embedded QM gradient (Ha/bohr). The z-vector step inside
        # the gradient already forms the excited-state ESPF charges, so
        # OQP::partial_charges holds the active state's QM charges afterwards.
        gqm = self._qm_gradient()
        pchg = np.array(mol.data["OQP::partial_charges"])
        # MM forces with embedded QM charges (OpenMM units)
        emm_q, gmm_q = self.driver.forces_mm(pchg)
        gmm = np.array(gmm_q.value_in_unit(u.kilojoule_per_mole / u.nanometer)) / HABOHR_TO_KJMOLNM
        emm = emm_q.value_in_unit(u.kilojoule_per_mole) * KJMOL_TO_HARTREE

        # total force = MM forces; on QM atoms subtract the QM gradient
        f_all = gmm.copy()
        for k, i in enumerate(self.qm_atoms):
            f_all[i] = f_all[i] - gqm[k]
        # periodic QM-QM Ewald self-interaction correction force (QM atoms only;
        # physically correct but small for large boxes -- NOT the dominant
        # source of the remaining periodic force-energy drift, which is the PME
        # embedding consistency term, still under development)
        if self.periodic:
            f_all = f_all + self._potqm_force(pchg, sign=1.0)
        # remove net (COM) force
        f_all -= f_all.mean(axis=0)

        # total potential energy (matches compute_force bookkeeping)
        eqm = float(mol.energies[self.active])
        znuc = np.array(mol.get_atoms2("charge"))
        eqm -= np.dot(pchg - znuc, potmm)
        epot = eqm + emm
        return f_all, epot

    # ------------------------------------------------------------------ #
    def run(self):
        mol = self.mol
        dump_log(mol, title='PyOQP: QM/MM Tully FSSH Nonadiabatic Molecular Dynamics')

        # initial electronic structure + force
        self._sync_positions()
        potmm0, _ = self._electronic_qmmm(with_overlap=False)
        f_all, epot = self._total_force(potmm0)
        accel = f_all / self.m_all[:, None]
        self.prev_xyz = copy.deepcopy(self.r_all[self.qm_atoms].reshape(-1))
        self.prev_data = copy.deepcopy(mol.get_data())
        self._log_qmmm(0, epot)

        for istep in range(1, self.nstep + 1):
            # velocity-Verlet position update (all atoms)
            self.r_all = self.r_all + self.v_all * self.dt + 0.5 * accel * self.dt ** 2
            self._sync_positions()

            # embedded electronic structure at the new geometry
            potmm, _ = self._electronic_qmmm(with_overlap=True)
            f_all, epot = self._total_force(potmm)
            accel_new = f_all / self.m_all[:, None]

            # velocity-Verlet velocity update (all atoms)
            self.v_all = self.v_all + 0.5 * (accel + accel_new) * self.dt

            # couplings + QM-only FSSH hop
            self._state_overlap()
            self.vel = self.v_all[self.qm_atoms].copy()       # hop sees QM velocities
            new_active, hopped = self._hop()
            self.v_all[self.qm_atoms] = self.vel              # write back rescaled QM velocities
            if hopped:
                self.active = new_active
                f_all, epot = self._total_force(potmm)
                accel_new = f_all / self.m_all[:, None]

            accel = accel_new
            self.prev_xyz = copy.deepcopy(self.r_all[self.qm_atoms].reshape(-1))
            self.prev_data = copy.deepcopy(mol.get_data())
            self._log_qmmm(istep, epot, hopped=hopped)

        dump_log(mol, title='PyOQP: QM/MM NAMD trajectory complete')

    def _log_qmmm(self, istep, epot, hopped=False):
        ekin = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
        pops = np.abs(self.coef) ** 2
        dump_log(
            self.mol,
            title=(f'QMMM-NAMD step {istep:6d}  t={istep*self.dt_fs:9.3f} fs  '
                   f'active={self.active}  E_tot={ekin+epot:.8f}  '
                   f'E_pot={epot:.8f}  E_kin={ekin:.8f}  '
                   f'hop={hopped}  pop={np.array2string(pops, precision=4)}'),
        )


HA_TO_WAVENUM = 219474.6313708


class NAMD_SOC(NAMD):
    """SOC-NAMD (intersystem crossing) on the SHARC spin-adiabatic representation.

    soc_mrsf builds and diagonalises H = diag(E_MCH) + H_SOC, giving the
    spin-adiabatic energies (OQP::soc_eval, cm^-1 relative to the lowest
    excitation) and eigenvectors U (OQP::soc_evec_*).  Surface hopping is done
    on these spin-adiabatic states.

    Nuclei propagate on the active spin-adiabatic surface using the SHARC
    weighted-MCH diagonal gradient (sum_i |U_i,a|^2 grad E_i^MCH, triplet Ms
    sublevels sharing one gradient), which is correct through an S/T crossing
    where SOC mixes states strongly.  The spin-adiabatic FSSH hopping layer uses
    the local-diabatization propagator with substep energy integration and
    U-phase tracking (block-Procrustes within degenerate Ms groups) on the
    block-diagonal MCH state overlap.
    """

    def __init__(self, mol):
        # set up nuclear/velocity state via the base class, then override the
        # electronic state space to the spin-adiabatic dimension.
        super().__init__(mol)
        self.nstate_mrsf = int(mol.config['tdhf']['nstate'])      # per multiplicity (ns = nt)
        self.ns = self.nstate_mrsf
        self.nt = self.nstate_mrsf
        self.nstate_soc = self.ns + 3 * self.nt
        # active spin-adiabatic state (1-based); [md] active is reused
        self.active = int(mol.config['md']['active'])
        # electronic amplitudes over the spin-adiabatic states
        self.coef = np.zeros(self.nstate_soc, dtype=complex)
        self.coef[self.active - 1] = 1.0 + 0.0j
        self.e_ref = 0.0
        self.e0 = 0.0
        # weight threshold for the SHARC weighted-MCH diagonal gradient
        try:
            self.grad_wthr = float(mol.config['md'].get('grad_wthr', 0.05))
        except Exception:
            self.grad_wthr = 0.05

    # ------------------------------------------------------------------ #
    def _electronic_soc(self, with_overlap=False):
        """SCF + singlet MRSF + triplet MRSF + soc_mrsf; returns (eval_ha, U).

        Stores the current singlet/triplet response vectors (for the MCH state
        overlap) and, when with_overlap, the MO overlap vs the previous geometry.
        """
        mol = self.mol
        sp = SinglePoint(mol)
        ref = sp.reference()
        self.e_ref = float(ref[0])

        if with_overlap:
            mol.back_door = (self.prev_xyz, self.prev_data)
            BasisOverlap(mol).overlap()                     # sets OQP::overlap_mo

        mol.data.set_tdhf_multiplicity(1)
        sing = sp.excitation(ref)
        self.sbvec = np.array(mol.data['OQP::td_bvec_mo']).copy()
        mol.data['OQP::td_singlet_energies'] = mol.data['OQP::td_energies'].copy()
        mol.data['OQP::td_bvec_mo_s'] = mol.data['OQP::td_bvec_mo'].copy()

        mol.data.set_tdhf_multiplicity(3)
        trip = sp.excitation(ref)
        self.tbvec = np.array(mol.data['OQP::td_bvec_mo']).copy()
        mol.data['OQP::td_triplet_energies'] = mol.data['OQP::td_energies'].copy()
        mol.data['OQP::td_bvec_mo_t'] = mol.data['OQP::td_bvec_mo'].copy()

        oqp.soc_mrsf(mol)

        eval_wn = np.array(mol.data['OQP::soc_eval']).reshape(-1)           # cm^-1 rel e0
        u = (np.array(mol.data['OQP::soc_evec_re'])
             + 1j * np.array(mol.data['OQP::soc_evec_im'])).reshape(self.nstate_soc, self.nstate_soc)
        self.e0 = float(min(np.array(sing[1:]).min() - self.e_ref,
                            np.array(trip[1:]).min() - self.e_ref)) if len(sing) > 1 else 0.0
        eval_ha = eval_wn / HA_TO_WAVENUM                                   # Hartree rel e0
        return eval_ha, u

    # ------------------------------------------------------------------ #
    # spin-adiabatic couplings (SHARC scheme)
    # ------------------------------------------------------------------ #
    def _mch_overlap(self):
        """Block-diagonal MCH state overlap S(t-dt,t) over the spin-adiabatic
        basis (ns singlets + 3*nt triplet Ms sublevels). Singlet and triplet
        spatial overlaps come from get_states_overlap; triplet Ms sublevels
        share the spatial overlap and are spin-orthogonal across Ms; singlet-
        triplet blocks vanish (different spin)."""
        mol = self.mol
        ns, nt, n = self.ns, self.nt, self.nstate_soc

        mol.data.set_tdhf_multiplicity(1)
        mol.data['OQP::td_bvec_mo'] = self.sbvec.copy()
        mol.data['OQP::td_bvec_mo_old'] = self.prev_sbvec.copy()
        oqp.get_states_overlap(mol)
        s_s = np.array(mol.data['OQP::td_states_overlap']).reshape((ns, ns))

        mol.data.set_tdhf_multiplicity(3)
        mol.data['OQP::td_bvec_mo'] = self.tbvec.copy()
        mol.data['OQP::td_bvec_mo_old'] = self.prev_tbvec.copy()
        oqp.get_states_overlap(mol)
        s_t = np.array(mol.data['OQP::td_states_overlap']).reshape((nt, nt))

        s = np.zeros((n, n))
        s[:ns, :ns] = s_s
        for m in range(3):
            for i in range(nt):
                for j in range(nt):
                    s[ns + i * 3 + m, ns + j * 3 + m] = s_t[i, j]
        return s

    @staticmethod
    def _phase_track(u, u_prev, s_mch, eval_cur, tol=5.0e-5):
        """Align the freshly diagonalised U to the previous step on the
        spin-adiabatic overlap T = U_prev^dag S_MCH U, using orthogonal
        Procrustes WITHIN each (near-)degenerate energy group only.

        Diagonalisation returns eigenvectors with arbitrary phase AND arbitrary
        rotation within degenerate subspaces (e.g. the three triplet Ms
        sublevels).  Restricting the alignment to degenerate blocks (which are
        adjacent since soc_eval is energy-sorted) removes that artifact while
        preserving the energy<->state correspondence (a global rotation would
        mix non-degenerate states and desynchronise eval).  Singleton groups
        reduce to a phase fix."""
        t = u_prev.conj().T @ s_mch @ u
        n = u.shape[1]
        w = np.eye(n, dtype=complex)
        i = 0
        while i < n:
            j = i + 1
            while j < n and abs(eval_cur[j] - eval_cur[i]) < tol:
                j += 1
            g = list(range(i, j))
            sub = t[np.ix_(g, g)]
            a_mat, _, bh = np.linalg.svd(sub)
            w[np.ix_(g, g)] = bh.conj().T @ a_mat.conj().T
            i = j
        u_aligned = u @ w
        t_aligned = u_prev.conj().T @ s_mch @ u_aligned
        return u_aligned, t_aligned

    def _propagate_and_hop(self, eval_prev, eval_cur, t):
        """Local-diabatization (SHARC) propagation of the spin-adiabatic
        amplitudes + fewest-switches hop + isotropic velocity rescaling.

        The orthonormalised spin-adiabatic overlap T (T[I,J]=<I(t)|J(t+dt)>) is
        used directly as the basis-change propagator, which is unitary and
        therefore robust to the arbitrary within-subspace rotation of degenerate
        states (e.g. the triplet Ms sublevels):

            P = diag(e^{-i E(t+dt) dt/2}) . T_u^dag . diag(e^{-i E(t) dt/2})
            c(t+dt) = P c(t)

        Hop probabilities (SHARC) attribute the active-state population loss to
        the states it flowed into through P.
        """
        from scipy.linalg import sqrtm, logm, expm
        n = self.nstate_soc
        a = self.active - 1
        dt = self.dt
        nsub = max(1, self.substep)

        tu = t @ np.linalg.inv(sqrtm(t.conj().T @ t))       # nearest unitary
        # substep local diabatization: split the basis rotation into nsub equal
        # fractional rotations (tu^{1/nsub}) and integrate the energy phase with
        # linearly interpolated diagonal energies.  The net full-step propagator
        # p is accumulated and used for the SHARC flux hop probabilities.
        # Reduces exactly to the single-step LD propagator when nsub = 1.
        kgen = logm(tu)                                     # skew-Hermitian generator
        rsub_dag = expm(-kgen / nsub)                       # (tu^{1/nsub})^dagger
        dtau = dt / nsub
        p = np.eye(n, dtype=complex)
        for s in range(nsub):
            ea = eval_prev + (eval_cur - eval_prev) * (s / nsub)
            eb = eval_prev + (eval_cur - eval_prev) * ((s + 1) / nsub)
            d1 = np.exp(-1j * ea * dtau / 2.0)
            d2 = np.exp(-1j * eb * dtau / 2.0)
            psub = (d2[:, None]) * rsub_dag * (d1[None, :])
            p = psub @ p                                   # propagator c(t+dt)=P c(t)

        c_old = self.coef.copy()
        c_new = p @ c_old
        nrm = np.linalg.norm(c_new)
        if nrm > 0:
            c_new = c_new / nrm

        # SHARC hop probabilities: distribute the active-state population loss
        rho_a = abs(c_old[a]) ** 2
        dp = rho_a - abs(c_new[a]) ** 2
        cmhp = np.zeros(n)
        if dp > 0.0 and rho_a > 1e-30:
            flux = np.array([max(0.0, np.real(np.conj(c_new[j]) * p[j, a] * c_old[a]))
                             for j in range(n)])
            flux[a] = 0.0
            fsum = flux.sum()
            if fsum > 1e-30:
                cmhp = (dp / rho_a) * flux / fsum

        # energy-based decoherence correction (Granucci-Persico)
        if self.decoherence == 1:
            ekin = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
            if ekin > 0:
                p_others = 0.0
                for k in range(n):
                    if k == a:
                        continue
                    gap = abs(eval_cur[k] - eval_cur[a])
                    if gap < 1e-12:
                        p_others += abs(c_new[k]) ** 2
                        continue
                    tau = (1.0 / gap) * (1.0 + self.edc_c / ekin)
                    c_new[k] *= np.exp(-dt / tau)
                    p_others += abs(c_new[k]) ** 2
                pa = abs(c_new[a]) ** 2
                if pa > 1e-30:
                    c_new[a] *= np.sqrt(max(0.0, 1.0 - p_others) / pa)
        self.coef = c_new

        # fewest-switches hop decision
        rand = float(self.rng.random())
        hopped = False
        lower = 0.0
        for j in range(n):
            if j == a:
                continue
            upper = lower + cmhp[j]
            if lower < rand < upper:
                de = eval_cur[a] - eval_cur[j]             # E_old - E_new
                ekin = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
                if de < 0.0 and ekin < abs(de):
                    break                                  # frustrated hop
                if abs(de) > self.thrshe:
                    break
                scale = np.sqrt(max(0.0, 1.0 + de / ekin)) if ekin > 0 else 1.0
                self.vel = scale * self.vel                # isotropic rescale
                self.active = j + 1
                hopped = True
                break
            lower = upper
        return hopped

    def _dominant_component(self, u, active):
        """Largest |U|^2 MCH component of the active spin-adiabatic state.
        Returns (multiplicity, state_index_1based, weight)."""
        col = np.abs(u[:, active - 1]) ** 2
        k = int(np.argmax(col))
        if k < self.ns:
            return 1, k + 1, col[k]                       # singlet S(k)
        kk = k - self.ns
        return 3, kk // 3 + 1, col[k]                     # triplet T(kk//3), Ms = kk%3

    def _soc_gradient(self, u, active, eval_ha):
        """Weighted-MCH (SHARC-diagonal) gradient of the active spin-adiabatic
        state:

            dE_diag,a/dR  ~  sum_i |U_i,a|^2  dE_i^MCH/dR

        neglecting the off-diagonal NAC terms and the (slowly varying) SOC
        derivative -- the standard SHARC diagonal-gradient approximation.  Only
        MCH components with weight above grad_wthr contribute, and the three
        triplet Ms sublevels of a spatial triplet share a single gradient (their
        weights are summed).  This is exact in the weak-mixing limit (one MCH
        component dominates) and, unlike the dominant-component approximation,
        gives the correct averaged force through an S/T crossing where SOC mixes
        states ~50/50.

        Returns (grad[natom,3] Hartree/bohr, E_diag absolute Hartree,
        dom_mult, dom_state, dom_weight) where the dominant labels are for
        logging only."""
        mol = self.mol
        col = np.abs(u[:, active - 1]) ** 2

        # collapse to unique MCH spatial states, summing triplet Ms weights
        wmap = {}
        for k in range(self.nstate_soc):
            if col[k] < self.grad_wthr:
                continue
            if k < self.ns:
                key = (1, k + 1)                          # singlet S(k+1)
            else:
                key = (3, (k - self.ns) // 3 + 1)         # triplet T(.), any Ms
            wmap[key] = wmap.get(key, 0.0) + col[k]

        if not wmap:                                      # safety: fall back to dominant
            kdom = int(np.argmax(col))
            key = (1, kdom + 1) if kdom < self.ns else (3, (kdom - self.ns) // 3 + 1)
            wmap[key] = float(col[kdom])

        wtot = sum(wmap.values())
        g = np.zeros((self.natom, 3))
        for (mult, state), w in wmap.items():
            mol.data.set_tdhf_multiplicity(mult)
            SinglePoint(mol).excitation([self.e_ref])     # set td vectors for this multiplicity
            mol.config['properties']['grad'] = [state]
            Gradient(mol).gradient()
            g += (w / wtot) * np.array(mol.grads[state]).reshape((self.natom, 3))

        dom_mult, dom_state, dom_w = self._dominant_component(u, active)
        e_diag = self.e_ref + self.e0 + float(eval_ha[active - 1])   # absolute (Hartree)
        return g, e_diag, dom_mult, dom_state, dom_w

    # ------------------------------------------------------------------ #
    def run(self):
        mol = self.mol
        dump_log(mol, title='PyOQP: SOC-NAMD (ISC, SHARC spin-adiabatic FSSH)')

        r = mol.get_system().reshape((self.natom, 3))

        # initial electronic structure + active-surface force
        eval_ha, u = self._electronic_soc(with_overlap=False)
        grad, e_pure, mult, state, w = self._soc_gradient(u, self.active, eval_ha)
        accel = -grad / self.mass[:, None]
        self._ulog = u
        self._store_prev(r, u, eval_ha)
        self._log_soc(0, e_pure, mult, state, w, False)

        for istep in range(1, self.nstep + 1):
            # velocity-Verlet position update
            r = r + self.vel * self.dt + 0.5 * accel * self.dt ** 2
            mol.update_system(r.reshape(-1))

            # electronic structure (+ MO overlap vs previous geometry)
            eval_ha, u = self._electronic_soc(with_overlap=True)

            # spin-adiabatic overlap: MCH overlap -> Procrustes-align U -> T
            s_mch = self._mch_overlap()
            u, t = self._phase_track(u, self.prev_u, s_mch, eval_ha)

            # active-surface force (weighted-MCH diagonal gradient) + vel update
            grad, e_pure, mult, state, w = self._soc_gradient(u, self.active, eval_ha)
            accel_new = -grad / self.mass[:, None]
            self.vel = self.vel + 0.5 * (accel + accel_new) * self.dt

            # local-diabatization propagation + fewest-switches hop
            hopped = self._propagate_and_hop(self.prev_eval, eval_ha, t)
            if hopped:
                grad, e_pure, mult, state, w = self._soc_gradient(u, self.active, eval_ha)
                accel_new = -grad / self.mass[:, None]

            accel = accel_new
            self._ulog = u
            self._store_prev(r, u, eval_ha)
            self._log_soc(istep, e_pure, mult, state, w, hopped)

        dump_log(mol, title='PyOQP: SOC-NAMD trajectory complete')

    def _store_prev(self, r, u, eval_ha):
        self.prev_xyz = copy.deepcopy(r.reshape(-1))
        self.prev_data = copy.deepcopy(self.mol.get_data())
        self.prev_u = u.copy()
        self.prev_eval = eval_ha.copy()
        self.prev_sbvec = self.sbvec.copy()
        self.prev_tbvec = self.tbvec.copy()

    def _log_soc(self, istep, e_pure, mult, state, w, hopped):
        ekin = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
        # manifold-summed populations via the MCH projection (U c): the spin
        # character is in the MCH basis, where the first ns components are
        # singlets and the rest triplet Ms sublevels. The adiabatic states are
        # energy-sorted mixtures, so summing adiabatic amplitudes by index is
        # not the spin character.
        mch = self._ulog @ self.coef
        pmch = np.abs(mch) ** 2
        pop_s = float(pmch[:self.ns].sum())
        pop_t = float(pmch[self.ns:].sum())
        dump_log(
            self.mol,
            title=(f'SOC-NAMD step {istep:6d}  t={istep*self.dt_fs:9.3f} fs  '
                   f'active={self.active}  E_tot={ekin+e_pure:.8f}  '
                   f'E_pure={e_pure:.8f}  E_kin={ekin:.8f}  hop={hopped}  '
                   f'dom=({"S" if mult==1 else "T"}{state},w={w:.3f})  '
                   f'pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )
