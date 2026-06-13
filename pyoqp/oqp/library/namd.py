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
        # adaptive (variable) timestep: shrink dt when the fastest atom would
        # move more than dx_max in one step (resolves stiff/hot modes without a
        # globally small dt). dt_max is the configured dt; clamped to dt_min.
        self.dt_max = self.dt
        _da = md.get('dt_adaptive', False)
        self.dt_adaptive = (_da is True) or (str(_da).lower() in ('true', '1', 'on', 'yes'))
        self.dt_min = float(md.get('dt_min', 0.05)) * FS_TO_AU
        self.dx_max = float(md.get('dx_max', 0.02))
        self._t_fs = 0.0                            # cumulative time (fs) for variable-dt logging
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

        # electronic amplitudes (complex), one per excited state. For SOC-NAMD
        # the active index runs over the larger spin-adiabatic manifold and the
        # subclass overwrites coef; guard the base indexing against that.
        self.coef = np.zeros(self.nstate, dtype=complex)
        if 1 <= self.active <= self.nstate:
            self.coef[self.active - 1] = 1.0 + 0.0j
        else:
            self.coef[0] = 1.0 + 0.0j

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

    def _adaptive_dt(self, vel, accel):
        """Return the timestep for this step (atomic units). If dt_adaptive,
        shrink dt so the largest predicted atomic displacement
        |v*dt + 1/2 a*dt^2| stays below dx_max; clamp to [dt_min, dt_max]."""
        if not getattr(self, 'dt_adaptive', False):
            return self.dt_max
        disp = np.abs(vel * self.dt_max + 0.5 * accel * self.dt_max ** 2)
        dmax = float(disp.max()) if disp.size else 0.0
        dt = self.dt_max * (self.dx_max / dmax) if dmax > self.dx_max else self.dt_max
        return float(max(self.dt_min, min(self.dt_max, dt)))

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
            active_old = self.active
            new_active, hopped = self._hop()

            active_changed = new_active != active_old
            if active_changed:
                self.active = new_active
                # force for the next step is on the new active surface. This
                # also covers trivial-crossing following, where the Fortran
                # kernel can update ACTIVE without marking HOPPED.
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
            title=(f'NAMD step {istep:6d}  t={(self._t_fs if self.dt_adaptive else istep*self.dt_fs):9.3f} fs  '
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
        # rigid-water (SHAKE/RATTLE) constraints for the MM region
        self._build_constraints()

    # ------------------------------------------------------------------ #
    def _build_constraints(self):
        """Collect the MM rigid-water bond/angle constraints (O-H, O-H, H-H per
        TIP3P water) from an OpenMM rigidWater system, as (i, j, d_bohr).  QM
        atoms are never constrained (they move under the QM forces).  Enables a
        normal MD timestep (~0.5-1 fs) despite the stiff O-H stretch."""
        u = self._u
        qm = set(int(i) for i in self.qm_atoms)
        ref = self.forcefield.createSystem(
            self.pdb.topology, nonbondedMethod=self.cutoff,
            constraints=None, rigidWater=True)
        ci, cj, cd = [], [], []
        for k in range(ref.getNumConstraints()):
            p1, p2, dist = ref.getConstraintParameters(k)
            if p1 in qm or p2 in qm:
                continue
            ci.append(p1); cj.append(p2)
            cd.append(dist.value_in_unit(u.nanometer) * NM_TO_BOHR)
        self._ci = np.array(ci, dtype=int)
        self._cj = np.array(cj, dtype=int)
        self._cd2 = np.array(cd) ** 2
        self._inv_m = 1.0 / self.m_all
        self._has_constraints = len(ci) > 0

    def _shake(self, r_old, r, v, dt, tol=1.0e-9, maxit=500):
        """Constrain bond lengths after the position update (SHAKE), and apply
        the implied velocity correction.  r is modified in place; v gets
        += (r_constrained - r_unconstrained)/dt."""
        if not self._has_constraints:
            return
        ci, cj, d2, inv = self._ci, self._cj, self._cd2, self._inv_m
        r_unc = r.copy()
        rij0 = r_old[ci] - r_old[cj]
        for _ in range(maxit):
            rij = r[ci] - r[cj]
            diff = np.einsum('ij,ij->i', rij, rij) - d2
            if np.max(np.abs(diff)) < tol:
                break
            denom = 2.0 * np.einsum('ij,ij->i', rij, rij0) * (inv[ci] + inv[cj])
            g = diff / denom
            dr = g[:, None] * rij0
            np.add.at(r, ci, -inv[ci][:, None] * dr)
            np.add.at(r, cj,  inv[cj][:, None] * dr)
        v += (r - r_unc) / dt

    def _thermalize_initial(self):
        """Rescale the full-system velocities to init_temp using the CONSTRAINED
        degrees of freedom (3N - n_constraints - 3 COM).  Called after the
        initial RATTLE has projected out the rigid-water internal velocities,
        which otherwise leaves the system below the target temperature.  Uniform
        scaling preserves both the RATTLE projection and zero COM momentum."""
        ncon = len(self._ci) if self._has_constraints else 0
        ndof = 3 * self.natom_all - ncon - 3
        if ndof <= 0:
            return
        ke = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
        if ke <= 0:
            return
        t_cur = 2.0 * ke / (ndof * KB_HARTREE)
        if t_cur > 0:
            self.v_all *= np.sqrt(self.init_temp / t_cur)

    def _rattle(self, r, v, tol=1.0e-9, maxit=500):
        """Project velocities onto the constraint manifold (RATTLE): make the
        relative velocity along each constrained bond zero.  v modified in place."""
        if not self._has_constraints:
            return
        ci, cj, inv = self._ci, self._cj, self._inv_m
        rij = r[ci] - r[cj]
        rr = np.einsum('ij,ij->i', rij, rij)
        for _ in range(maxit):
            vij = v[ci] - v[cj]
            rv = np.einsum('ij,ij->i', rij, vij)
            if np.max(np.abs(rv)) < tol:
                break
            k = rv / (rr * (inv[ci] + inv[cj]))
            dv = k[:, None] * rij
            np.add.at(v, ci, -inv[ci][:, None] * dv)
            np.add.at(v, cj,  inv[cj][:, None] * dv)

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
        # Zero POTQM: POTMM (PME) already captures the periodic MM embedding and
        # has the QM self-image removed; the residual QM-QM periodic image term
        # is negligible for solvation-size boxes, and the OpenMM correction was
        # buggy (over-corrected E by ~5 Ha, force-inconsistent -- verified by
        # finite difference). See pme_fd_diag.py.
        mol.data["OQP::POTQM"] = np.zeros((nat, nat))
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
        import os
        mol = self.mol
        mol.config['properties']['grad'] = [self.active]
        Gradient(mol).gradient()
        g = np.array(mol.grads[self.active]).reshape(-1, 3)
        # ESPF_ROHF=1: use the ROHF reference density for ESPF charges and
        # gradient, matching GAMESS which always fits ESPF from the SCF (ROHF)
        # density regardless of the target excited state.  Combined with
        # ESPF_GAMESS=1 this reproduces GAMESS QM/MM energy conservation for
        # direct validation.  Default (ESPF_ROHF unset): physically correct
        # S1 relaxed density via grad_esp_qmmm_excited.
        if os.environ.get('ESPF_ROHF', '').strip() in ('1', 'on'):
            oqp.form_esp_charges(mol)   # partial_charges from ROHF density
            oqp.grad_esp_qmmm(mol)      # ESPF gradient from ROHF density
        else:
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
        # No POTQM force: the QM-QM periodic image self-interaction is neglected
        # (POTQM zeroed in the embedded SCF; see _electronic_qmmm). Adding the
        # _potqm_force here without the matching energy term would reintroduce a
        # force-energy inconsistency.
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
        self._rattle(self.r_all, self.v_all)          # project initial MM velocities onto constraints
        self._thermalize_initial()                    # rescale to init_temp on the constrained DOF
        self.prev_xyz = copy.deepcopy(self.r_all[self.qm_atoms].reshape(-1))
        self.prev_data = copy.deepcopy(mol.get_data())
        self._log_qmmm(0, epot)

        for istep in range(1, self.nstep + 1):
            # velocity-Verlet position update (all atoms) + SHAKE (rigid MM water)
            # (fixed dt: the same-spin path uses the Fortran hop kernel with dt_fs)
            r_old = self.r_all.copy()
            self.r_all = self.r_all + self.v_all * self.dt + 0.5 * accel * self.dt ** 2
            self._shake(r_old, self.r_all, self.v_all, self.dt)
            self._sync_positions()

            # embedded electronic structure at the new geometry
            potmm, _ = self._electronic_qmmm(with_overlap=True)
            f_all, epot = self._total_force(potmm)
            accel_new = f_all / self.m_all[:, None]

            # velocity-Verlet velocity update (all atoms) + RATTLE (rigid MM water)
            self.v_all = self.v_all + 0.5 * (accel + accel_new) * self.dt
            self._rattle(self.r_all, self.v_all)

            # couplings + QM-only FSSH hop
            self._state_overlap()
            self.vel = self.v_all[self.qm_atoms].copy()       # hop sees QM velocities
            active_old = self.active
            new_active, hopped = self._hop()
            self.v_all[self.qm_atoms] = self.vel              # write back rescaled QM velocities
            active_changed = new_active != active_old
            if active_changed:
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
            title=(f'QMMM-NAMD step {istep:6d}  t={(self._t_fs if self.dt_adaptive else istep*self.dt_fs):9.3f} fs  '
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
            self.grad_wthr = float(mol.config['md'].get('grad_wthr', 0.001))
        except Exception:
            self.grad_wthr = 0.05
        # optional: choose the initial active state by MCH spin character
        # (e.g. 'S1') instead of a fixed spin-adiabatic index -- robust when the
        # spin-adiabatic ordering is ambiguous at the start (S/T near-degeneracy)
        self.init_state = str(mol.config['md'].get('init_state', '') or '').strip()
        _ev = mol.config['md'].get('econs', False)
        self.econs = (_ev is True) or (str(_ev).lower() in ('true', '1', 'on', 'yes'))
        _du = mol.config['md'].get('soc_du_dt_corr', False)
        self.soc_du_dt_corr = (_du is True) or (str(_du).lower() in ('true', '1', 'on', 'yes'))
        _tdcg = mol.config['md'].get('soc_tdc_grad_corr', False)
        self.soc_tdc_grad_corr = (_tdcg is True) or (str(_tdcg).lower() in ('true', '1', 'on', 'yes'))

    # ------------------------------------------------------------------ #
    def _resolve_initial_active(self, u):
        """If [md] init_state names an MCH state (S0/S1/.../T1/T2/...), set the
        active spin-adiabatic state to the adiabat with the largest character of
        that MCH state at t=0 (summing the three Ms sublevels for a triplet).
        Otherwise keep the configured integer active index."""
        label = self.init_state
        if not label:
            return
        mult = 1 if label[0].upper() == 'S' else 3
        n = int(label[1:])
        if mult == 1:
            mch_idx = [n]                                  # singlet root: S0=0, S1=1, ...
        else:
            base = self.ns + (n - 1) * 3
            mch_idx = [base, base + 1, base + 2]           # triplet Ms sublevels
        char = (np.abs(u[mch_idx, :]) ** 2).sum(axis=0)    # character per adiabat
        a = int(np.argmax(char))
        self.active = a + 1
        self.coef = np.zeros(self.nstate_soc, dtype=complex)
        self.coef[a] = 1.0 + 0.0j
        dump_log(self.mol, title=(
            f'SOC-NAMD: initial active set to adiabat {self.active} by {label} '
            f'character ({char[a]*100:.1f}% {label})'))

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
        self.sing_energies = np.array(sing, dtype=float)
        self.sbvec = np.array(mol.data['OQP::td_bvec_mo']).copy()
        mol.data['OQP::td_singlet_energies'] = mol.data['OQP::td_energies'].copy()
        mol.data['OQP::td_bvec_mo_s'] = mol.data['OQP::td_bvec_mo'].copy()

        mol.data.set_tdhf_multiplicity(3)
        trip = sp.excitation(ref)
        self.trip_energies = np.array(trip, dtype=float)
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

    def _mch_target(self, k):
        """Map an MCH (diabatic) basis index k to its (multiplicity, MRSF grad
        target).  MRSF roots are 1-based with the LOWEST root being the ground
        state: root 1 = S0, root 2 = S1, ...  (S0 is the lowest eigenvalue of
        the MRSF orbital-Hessian response, so it has a normal MRSF gradient.)
        Hence the singlet block index k maps to target k+1 (k=0->S0, k=1->S1).
        The triplet block is also 1-based per multiplicity (T1=target 1, ...);
        the three Ms sublevels of a spatial triplet share one target."""
        if k < self.ns:
            return 1, k + 1                               # singlet root: S0=1, S1=2, ...
        return 3, (k - self.ns) // 3 + 1                  # triplet root: T1=1, T2=2, ...

    @staticmethod
    def _mch_label(mult, target):
        """Human-readable MCH state name for a (mult, MRSF target): singlet
        target t -> S(t-1) (so target 1 = S0), triplet target t -> T(t)."""
        return f'S{target - 1}' if mult == 1 else f'T{target}'

    def _mch_energies_abs(self):
        """Absolute MCH energies expanded over singlet + triplet Ms sublevels."""
        e = []
        for target in range(1, self.ns + 1):
            e.append(float(self.sing_energies[target]))
        for target in range(1, self.nt + 1):
            e.extend([float(self.trip_energies[target])] * 3)
        return np.array(e)

    def _mch_hamiltonian_from_u(self, u, eval_ha):
        """MCH-basis Hamiltonian, relative to the common e0 shift, in Hartree."""
        h = u @ np.diag(eval_ha) @ u.conj().T
        return 0.5 * (h + h.conj().T)

    def _build_wmap(self, col):
        """{(mult,target): weight} map of MCH components contributing to the
        active spin-adiabatic state's gradient, keeping components above
        grad_wthr; triplet Ms sublevels share a target (weights summed).
        Falls back to the dominant component if none clear the threshold."""
        wmap = {}
        for k in range(self.nstate_soc):
            if col[k] < self.grad_wthr:
                continue
            key = self._mch_target(k)
            wmap[key] = wmap.get(key, 0.0) + col[k]
        if not wmap:
            kdom = int(np.argmax(col))
            wmap[self._mch_target(kdom)] = float(col[kdom])
        return wmap

    def _dominant_component(self, u, active):
        """Largest |U|^2 MCH component of the active spin-adiabatic state.
        Returns (multiplicity, MRSF grad target, weight)."""
        col = np.abs(u[:, active - 1]) ** 2
        k = int(np.argmax(col))
        mult, target = self._mch_target(k)
        return mult, target, col[k]

    def _du_dt_gradient_correction(self, u, active, eval_ha, vel):
        """Option 2: projected finite-difference dU/dt correction to dE/dR.

        The phase-tracked active-column derivative gives dU/dt along the actual
        nuclear displacement.  The minimum-norm spatial projection is
        dU/dR ~= dU/dt * v / |v|^2, yielding a gradient correction in
        Hartree/bohr with the same shape as the QM gradient.
        """
        if not getattr(self, 'soc_du_dt_corr', False):
            return np.zeros((self.natom, 3))
        if not hasattr(self, 'prev_u') or self.prev_u is None or self.dt <= 0:
            return np.zeros((self.natom, 3))

        v = np.array(vel, dtype=float).reshape((self.natom, 3))
        v2 = float(np.sum(v * v))
        if v2 < 1.0e-30:
            return np.zeros((self.natom, 3))

        a = active - 1
        du_dt = (u - self.prev_u) / self.dt
        coeff = 0.0
        for i in range(self.nstate_soc):
            for j in range(self.nstate_soc):
                if i == j:
                    continue
                coeff += 2.0 * np.real(
                    u[i, a].conj() * u[j, a] * (eval_ha[j] - eval_ha[i]) * du_dt[i, a]
                )
        g_corr = coeff * v / v2
        self._du_dt_corr_norm = float(np.linalg.norm(g_corr))
        return g_corr

    def _tdc_gradient_correction(self, u, active, s_mch, vel):
        """Approximate the off-diagonal MCH derivative-Hamiltonian force term.

        The exact SHARC diagonal gradient contains MCH NAC vectors through

            G_ij = (E_j - E_i) d_ij,  i != j.

        We already have overlap-derived time-derivative couplings for the TDSE,
        tau_ij = d_ij dot v.  This correction uses the minimum-norm projection
        d_ij ~= tau_ij * v / |v|^2, giving an approximate vector correction
        without additional QM calls.  The SOC derivative term is still omitted.
        """
        if not getattr(self, 'soc_tdc_grad_corr', False):
            return np.zeros((self.natom, 3))
        if s_mch is None:
            return np.zeros((self.natom, 3))

        v = np.array(vel, dtype=float).reshape((self.natom, 3))
        v2 = float(np.sum(v * v))
        if v2 < 1.0e-30:
            return np.zeros((self.natom, 3))

        a = active - 1
        tau = self._compute_tdc(np.array(s_mch, dtype=float).reshape((self.nstate_soc, self.nstate_soc)))
        e_mch = self._mch_energies_abs()
        coeff = 0.0
        for i in range(self.nstate_soc):
            for j in range(self.nstate_soc):
                if i == j:
                    continue
                coeff += np.real(
                    u[i, a].conj() * u[j, a] * (e_mch[j] - e_mch[i]) * tau[i, j]
                )
        g_corr = coeff * v / v2
        self._tdc_grad_corr_norm = float(np.linalg.norm(g_corr))
        return g_corr

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

        # collapse to unique MCH spatial states (summing triplet Ms weights)
        wmap = self._build_wmap(col)
        wtot = sum(wmap.values())
        g = np.zeros((self.natom, 3))
        for (mult, state), w in wmap.items():
            mol.data.set_tdhf_multiplicity(mult)
            SinglePoint(mol).excitation([self.e_ref])     # set td vectors for this multiplicity
            mol.config['properties']['grad'] = [state]
            Gradient(mol).gradient()
            g += (w / wtot) * np.array(mol.grads[state]).reshape((self.natom, 3))
        g += self._du_dt_gradient_correction(u, active, eval_ha, self.vel)
        g += self._tdc_gradient_correction(
            u, active, getattr(self, '_last_s_mch', None), self.vel)

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
        self._resolve_initial_active(u)
        grad, e_pure, mult, state, w = self._soc_gradient(u, self.active, eval_ha)
        accel = -grad / self.mass[:, None]
        self._ulog = u
        self._store_prev(r, u, eval_ha)
        self._log_soc(0, e_pure, mult, state, w, False)
        self._e_ref_tot = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2) + e_pure

        for istep in range(1, self.nstep + 1):
            # adaptive timestep + velocity-Verlet position update
            self.dt = self._adaptive_dt(self.vel, accel)
            self._t_fs += self.dt / FS_TO_AU
            r = r + self.vel * self.dt + 0.5 * accel * self.dt ** 2
            mol.update_system(r.reshape(-1))

            # electronic structure (+ MO overlap vs previous geometry)
            eval_ha, u = self._electronic_soc(with_overlap=True)

            # spin-adiabatic overlap: MCH overlap -> Procrustes-align U -> T
            s_mch = self._mch_overlap()
            self._last_s_mch = s_mch
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
            if self.econs:                                 # temporary E_tot-conservation rescale
                ke = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
                ket = self._e_ref_tot - e_pure
                if ke > 0 and ket > 0:
                    self.vel *= np.sqrt(ket / ke)
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
            title=(f'SOC-NAMD step {istep:6d}  t={(self._t_fs if self.dt_adaptive else istep*self.dt_fs):9.3f} fs  '
                   f'active={self.active}  E_tot={ekin+e_pure:.8f}  '
                   f'E_pure={e_pure:.8f}  E_kin={ekin:.8f}  hop={hopped}  '
                   f'dom=({self._mch_label(mult, state)},w={w:.3f})  '
                   f'pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )


class NAMD_SOC_MCH(NAMD_SOC):
    """SOC-NAMD in the MCH (spin-pure) basis.

    The active state is a single MCH basis function (singlet root or one
    triplet Ms component), so the nuclear force is the exact MCH root gradient.
    Electronic amplitudes are propagated by the SOC Hamiltonian in that MCH
    basis instead of the spin-adiabatic local-diabatization propagator.
    """

    def __init__(self, mol):
        super().__init__(mol)
        self.coef = np.zeros(self.nstate_soc, dtype=complex)
        self.coef[self.active - 1] = 1.0 + 0.0j

    def _resolve_initial_mch_active(self):
        label = self.init_state
        if not label:
            return
        mult = 1 if label[0].upper() == 'S' else 3
        target = int(label[1:])
        if mult == 1:
            active = target + 1                         # S0 -> MCH basis 1
        else:
            active = self.ns + (target - 1) * 3 + 1      # choose Ms=-1 member
        if not (1 <= active <= self.nstate_soc):
            raise ValueError(f"[md] init_state='{label}' is outside the SOC MCH basis")
        self.active = active
        self.coef[:] = 0.0
        self.coef[self.active - 1] = 1.0 + 0.0j
        dump_log(self.mol, title=f'SOC-MCH-NAMD: initial active set to MCH state {self._mch_active_label(self.active)}')

    def _mch_active_label(self, active):
        k = active - 1
        mult, state = self._mch_target(k)
        if mult == 1:
            return self._mch_label(mult, state)
        ms = (k - self.ns) % 3 - 1
        return f'{self._mch_label(mult, state)}(ms={ms:+d})'

    def _mch_exact_gradient(self, active):
        mol = self.mol
        mult, state = self._mch_target(active - 1)
        mol.data.set_tdhf_multiplicity(mult)
        SinglePoint(mol).excitation([self.e_ref])
        mol.config['properties']['grad'] = [state]
        Gradient(mol).gradient()
        g = np.array(mol.grads[state]).reshape((self.natom, 3))
        e = self._mch_energies_abs()[active - 1]
        return g, e, mult, state

    def _mch_propagate_and_hop(self, h_mch, e_mch):
        from scipy.linalg import expm
        n = self.nstate_soc
        a = self.active - 1
        dt = self.dt

        c_old = self.coef.copy()
        p = expm(-1j * h_mch * dt)
        c_new = p @ c_old
        nrm = np.linalg.norm(c_new)
        if nrm > 0:
            c_new /= nrm

        rho_a = abs(c_old[a]) ** 2
        cmhp = np.zeros(n)
        if rho_a > 1.0e-30:
            for j in range(n):
                if j == a:
                    continue
                # TDSE in a.u.: c_dot = -i H c. Positive loss of active
                # population through channel a->j becomes a hop probability.
                loss = 2.0 * np.real(1j * c_old[a].conj() * h_mch[a, j] * c_old[j])
                cmhp[j] = max(0.0, dt * loss / rho_a)

        if self.decoherence == 1:
            ekin = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
            if ekin > 0:
                p_others = 0.0
                for k in range(n):
                    if k == a:
                        continue
                    gap = abs(e_mch[k] - e_mch[a])
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

        rand = float(self.rng.random())
        hopped = False
        lower = 0.0
        for j in range(n):
            if j == a:
                continue
            upper = lower + cmhp[j]
            if lower < rand < upper:
                de = e_mch[a] - e_mch[j]
                ekin = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
                if de < 0.0 and ekin < abs(de):
                    break
                if abs(de) > self.thrshe:
                    break
                scale = np.sqrt(max(0.0, 1.0 + de / ekin)) if ekin > 0 else 1.0
                self.vel = scale * self.vel
                self.active = j + 1
                hopped = True
                break
            lower = upper
        return hopped

    def run(self):
        mol = self.mol
        dump_log(mol, title='PyOQP: SOC-NAMD (ISC, MCH-basis FSSH)')

        r = mol.get_system().reshape((self.natom, 3))
        eval_ha, u = self._electronic_soc(with_overlap=False)
        self._resolve_initial_mch_active()
        h_mch = self._mch_hamiltonian_from_u(u, eval_ha)
        e_mch = self._mch_energies_abs()
        grad, e_pure, mult, state = self._mch_exact_gradient(self.active)
        accel = -grad / self.mass[:, None]
        self._store_prev(r, u, eval_ha)
        self._log_mch(0, e_pure, mult, state, False)
        self._e_ref_tot = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2) + e_pure

        for istep in range(1, self.nstep + 1):
            self.dt = self._adaptive_dt(self.vel, accel)
            self._t_fs += self.dt / FS_TO_AU
            r = r + self.vel * self.dt + 0.5 * accel * self.dt ** 2
            mol.update_system(r.reshape(-1))

            eval_ha, u = self._electronic_soc(with_overlap=False)
            h_mch = self._mch_hamiltonian_from_u(u, eval_ha)
            e_mch = self._mch_energies_abs()
            grad, e_pure, mult, state = self._mch_exact_gradient(self.active)
            accel_new = -grad / self.mass[:, None]
            self.vel = self.vel + 0.5 * (accel + accel_new) * self.dt

            hopped = self._mch_propagate_and_hop(h_mch, e_mch)
            if hopped:
                grad, e_pure, mult, state = self._mch_exact_gradient(self.active)
                accel_new = -grad / self.mass[:, None]

            accel = accel_new
            if self.econs:
                ke = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
                ket = self._e_ref_tot - e_pure
                if ke > 0 and ket > 0:
                    self.vel *= np.sqrt(ket / ke)
            self._store_prev(r, u, eval_ha)
            self._log_mch(istep, e_pure, mult, state, hopped)

        dump_log(mol, title='PyOQP: SOC-MCH-NAMD trajectory complete')

    def _log_mch(self, istep, e_pure, mult, state, hopped):
        ekin = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
        pmch = np.abs(self.coef) ** 2
        pop_s = float(pmch[:self.ns].sum())
        pop_t = float(pmch[self.ns:].sum())
        dump_log(
            self.mol,
            title=(f'SOC-MCH-NAMD step {istep:6d}  t={(self._t_fs if self.dt_adaptive else istep*self.dt_fs):9.3f} fs  '
                   f'active={self.active}:{self._mch_active_label(self.active)}  '
                   f'E_tot={ekin+e_pure:.8f}  E_pure={e_pure:.8f}  '
                   f'E_kin={ekin:.8f}  hop={hopped}  '
                   f'grad={self._mch_label(mult, state)}  pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )


class NAMD_SOC_QMMM(NAMD_QMMM):
    """SOC-NAMD (intersystem crossing) with electrostatic ESPF QM/MM embedding.

    Union of the SHARC spin-adiabatic SOC-NAMD machinery (NAMD_SOC) and the
    ESPF/OpenMM embedding (NAMD_QMMM).  Per step:
      * sync positions (QM Molecule + OpenMM context),
      * embedded SCF + singlet MRSF + triplet MRSF + soc_mrsf -> (E_diag, U),
      * spin-adiabatic MCH overlap -> U-phase tracking -> overlap T,
      * active-surface force = weighted-MCH diagonal gradient, each MCH
        component carrying its own ESPF embedding gradient,
      * ESPF QM charges (of the dominant MCH component) -> MM forces,
      * full-system velocity Verlet (QM+MM, atomic units),
      * local-diabatization propagation + spin-adiabatic fewest-switches hop,
        rescaling QM velocities only (as in GAMESS RESCALV).

    The SOC electronic/hopping kernels are borrowed from NAMD_SOC via explicit
    NAMD_SOC.<method>(self, ...) calls so this class can inherit the QM/MM
    embedding plumbing from NAMD_QMMM.
    """

    # borrow the small SOC helpers so they resolve via self inside the borrowed
    # NAMD_SOC methods (this class inherits NAMD_QMMM, not NAMD_SOC)
    _mch_target = NAMD_SOC._mch_target
    _mch_label = staticmethod(NAMD_SOC._mch_label)
    _build_wmap = NAMD_SOC._build_wmap
    _dominant_component = NAMD_SOC._dominant_component
    _mch_energies_abs = NAMD_SOC._mch_energies_abs
    _mch_hamiltonian_from_u = NAMD_SOC._mch_hamiltonian_from_u

    def __init__(self, mol):
        super().__init__(mol)                                  # NAMD_QMMM: OpenMM + QM masses + v_all
        # spin-adiabatic electronic state space (ns singlets + 3 nt triplet Ms)
        self.nstate_mrsf = int(mol.config['tdhf']['nstate'])
        self.ns = self.nstate_mrsf
        self.nt = self.nstate_mrsf
        self.nstate_soc = self.ns + 3 * self.nt
        self.active = int(mol.config['md']['active'])
        self.coef = np.zeros(self.nstate_soc, dtype=complex)
        self.coef[self.active - 1] = 1.0 + 0.0j
        self.e_ref = 0.0
        self.e0 = 0.0
        try:
            self.grad_wthr = float(mol.config['md'].get('grad_wthr', 0.001))
        except Exception:
            self.grad_wthr = 0.05
        self.init_state = str(mol.config['md'].get('init_state', '') or '').strip()
        _ev = mol.config['md'].get('econs', False)
        self.econs = (_ev is True) or (str(_ev).lower() in ('true', '1', 'on', 'yes'))
        _du = mol.config['md'].get('soc_du_dt_corr', False)
        self.soc_du_dt_corr = (_du is True) or (str(_du).lower() in ('true', '1', 'on', 'yes'))
        _tdcg = mol.config['md'].get('soc_tdc_grad_corr', False)
        self.soc_tdc_grad_corr = (_tdcg is True) or (str(_tdcg).lower() in ('true', '1', 'on', 'yes'))

    # ------------------------------------------------------------------ #
    def _electronic_soc_qmmm(self, with_overlap):
        """Embedded SCF + singlet MRSF + triplet MRSF + soc_mrsf.
        Returns (eval_ha rel e0, U, potmm, potqm)."""
        from oqp.library.qmmm_driver import (
            unpack_lower_tri_single, unpack_lower_tri_multi, pack_lower_tri_single)
        mol = self.mol
        potmm, potqm = self.driver.electrostatic_potential()

        sp = SinglePoint(mol)
        sp._prep_guess()
        nat = mol.data["natom"]
        nbf = mol.data.get_basis()["nbf"]
        mol.data["OQP::POTMM"] = potmm
        # Zero POTQM: POTMM (PME) already captures the periodic MM embedding and
        # has the QM self-image removed; the residual QM-QM periodic image term
        # is negligible for solvation-size boxes, and the OpenMM correction was
        # buggy (over-corrected E by ~5 Ha, force-inconsistent -- verified by
        # finite difference). See pme_fd_diag.py.
        mol.data["OQP::POTQM"] = np.zeros((nat, nat))
        oqp.espf_op_corr(mol)
        espf = unpack_lower_tri_multi(mol.data["OQP::ESPF_CORR"], nbf, nat)
        hcore = unpack_lower_tri_single(mol.get_hcore(), nbf)
        hcore += np.einsum("ijk,i->jk", espf, potmm)
        mol.set_hcore(pack_lower_tri_single(hcore))
        sp.scf()
        ref = [mol.get_scf_energy()]
        self.e_ref = float(ref[0])

        if with_overlap:
            mol.back_door = (self.prev_xyz, self.prev_data)
            BasisOverlap(mol).overlap()

        mol.data.set_tdhf_multiplicity(1)
        sing = sp.excitation(ref)
        self.sing_energies = np.array(sing, dtype=float)
        self.sbvec = np.array(mol.data['OQP::td_bvec_mo']).copy()
        mol.data['OQP::td_singlet_energies'] = mol.data['OQP::td_energies'].copy()
        mol.data['OQP::td_bvec_mo_s'] = mol.data['OQP::td_bvec_mo'].copy()

        mol.data.set_tdhf_multiplicity(3)
        trip = sp.excitation(ref)
        self.trip_energies = np.array(trip, dtype=float)
        self.tbvec = np.array(mol.data['OQP::td_bvec_mo']).copy()
        mol.data['OQP::td_triplet_energies'] = mol.data['OQP::td_energies'].copy()
        mol.data['OQP::td_bvec_mo_t'] = mol.data['OQP::td_bvec_mo'].copy()

        oqp.soc_mrsf(mol)

        eval_wn = np.array(mol.data['OQP::soc_eval']).reshape(-1)            # cm^-1 rel e0
        u = (np.array(mol.data['OQP::soc_evec_re'])
             + 1j * np.array(mol.data['OQP::soc_evec_im'])).reshape(self.nstate_soc, self.nstate_soc)
        self.e0 = float(min(np.array(sing[1:]).min() - self.e_ref,
                            np.array(trip[1:]).min() - self.e_ref)) if len(sing) > 1 else 0.0
        eval_ha = eval_wn / HA_TO_WAVENUM                                    # Hartree rel e0
        return eval_ha, u, potmm, potqm

    # ------------------------------------------------------------------ #
    def _soc_gradient_qmmm(self, u, active, eval_ha):
        """Weighted-MCH diagonal gradient with ESPF embedding force per MCH
        component, plus the dominant component's ESPF QM charges for the MM
        forces.  Returns (grad_qm[natom,3], E_diag, dom_mult, dom_state,
        dom_w, pchg_dominant)."""
        mol = self.mol
        col = np.abs(u[:, active - 1]) ** 2

        wmap = NAMD_SOC._build_wmap(self, col)
        dom_mult, dom_state, dom_w = NAMD_SOC._dominant_component(self, u, active)
        dom_key = (dom_mult, dom_state)
        wtot = sum(wmap.values())
        g = np.zeros((self.natom, 3))
        pchg_dom = None
        for (mult, state), w in wmap.items():
            mol.data.set_tdhf_multiplicity(mult)
            SinglePoint(mol).excitation([self.e_ref])
            mol.config['properties']['grad'] = [state]
            Gradient(mol).gradient()
            gi = np.array(mol.grads[state]).reshape((self.natom, 3))
            # ESPF_ROHF=1: use ROHF reference density for ESPF gradient across
            # all SOC MCH components, matching GAMESS behaviour and ensuring
            # the ESPF energy is constant across state hops.
            if os.environ.get('ESPF_ROHF', '').strip() in ('1', 'on'):
                oqp.form_esp_charges(mol)
                oqp.grad_esp_qmmm(mol)
            else:
                oqp.grad_esp_qmmm_excited(mol)
            gi = gi + np.array(mol.data["OQP::ESPF_GRAD"]).reshape((self.natom, 3))
            g += (w / wtot) * gi
            if (mult, state) == dom_key:
                pchg_dom = np.array(mol.data["OQP::partial_charges"]).copy()
        g += NAMD_SOC._du_dt_gradient_correction(
            self, u, active, eval_ha, self.v_all[self.qm_atoms])
        g += NAMD_SOC._tdc_gradient_correction(
            self, u, active, getattr(self, '_last_s_mch', None), self.v_all[self.qm_atoms])

        if pchg_dom is None:                                  # dominant below threshold: take last
            pchg_dom = np.array(mol.data["OQP::partial_charges"]).copy()

        e_diag = self.e_ref + self.e0 + float(eval_ha[active - 1])
        return g, e_diag, dom_mult, dom_state, dom_w, pchg_dom

    # ------------------------------------------------------------------ #
    def _total_force_soc(self, potmm, g_qm, e_diag, pchg):
        """Assemble full-system force (a.u.) and total potential energy (Ha)."""
        mol = self.mol
        u = self._u
        emm_q, gmm_q = self.driver.forces_mm(pchg)
        gmm = np.array(gmm_q.value_in_unit(u.kilojoule_per_mole / u.nanometer)) / HABOHR_TO_KJMOLNM
        emm = emm_q.value_in_unit(u.kilojoule_per_mole) * KJMOL_TO_HARTREE

        f_all = gmm.copy()
        for k, i in enumerate(self.qm_atoms):
            f_all[i] = f_all[i] - g_qm[k]
        # No POTQM force: the QM-QM periodic image self-interaction is neglected
        # (POTQM zeroed in the embedded SCF; see _electronic_qmmm). Adding the
        # _potqm_force here without the matching energy term would reintroduce a
        # force-energy inconsistency.
        f_all -= f_all.mean(axis=0)

        eqm = float(e_diag)
        znuc = np.array(mol.get_atoms2("charge"))
        eqm -= np.dot(pchg - znuc, potmm)
        epot = eqm + emm
        return f_all, epot

    # ------------------------------------------------------------------ #
    def run(self):
        mol = self.mol
        dump_log(mol, title='PyOQP: SOC-NAMD QM/MM (ISC, SHARC spin-adiabatic FSSH + ESPF embedding)')

        self._sync_positions()
        eval_ha, u, potmm, _ = self._electronic_soc_qmmm(with_overlap=False)
        NAMD_SOC._resolve_initial_active(self, u)
        g_qm, e_diag, mult, state, w, pchg = self._soc_gradient_qmmm(u, self.active, eval_ha)
        f_all, epot = self._total_force_soc(potmm, g_qm, e_diag, pchg)
        accel = f_all / self.m_all[:, None]
        self._rattle(self.r_all, self.v_all)          # project initial MM velocities onto constraints
        self._thermalize_initial()                    # rescale to init_temp on the constrained DOF
        self._ulog = u
        r_qm = self.r_all[self.qm_atoms].reshape((self.natom, 3))
        NAMD_SOC._store_prev(self, r_qm, u, eval_ha)
        self._log_soc_qmmm(0, epot, mult, state, w, False)
        self._e_ref_tot = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2) + epot

        for istep in range(1, self.nstep + 1):
            # adaptive timestep + velocity-Verlet position update + SHAKE
            self.dt = self._adaptive_dt(self.v_all, accel)
            self._t_fs += self.dt / FS_TO_AU
            r_old = self.r_all.copy()
            self.r_all = self.r_all + self.v_all * self.dt + 0.5 * accel * self.dt ** 2
            self._shake(r_old, self.r_all, self.v_all, self.dt)
            self._sync_positions()

            # embedded spin-adiabatic electronic structure (+ MO overlap)
            eval_ha, u, potmm, _ = self._electronic_soc_qmmm(with_overlap=True)
            s_mch = NAMD_SOC._mch_overlap(self)
            self._last_s_mch = s_mch
            u, t = NAMD_SOC._phase_track(u, self.prev_u, s_mch, eval_ha)

            # active-surface force (weighted-MCH diagonal gradient + ESPF) + vel update
            g_qm, e_diag, mult, state, w, pchg = self._soc_gradient_qmmm(u, self.active, eval_ha)
            f_all, epot = self._total_force_soc(potmm, g_qm, e_diag, pchg)
            accel_new = f_all / self.m_all[:, None]
            self.v_all = self.v_all + 0.5 * (accel + accel_new) * self.dt
            self._rattle(self.r_all, self.v_all)

            # local-diabatization propagation + spin-adiabatic hop (QM velocities only)
            active_old = self.active                           # save for ESPF correction below
            epot_old = epot                                    # total E_pot before hop
            self.vel = self.v_all[self.qm_atoms].copy()
            hopped = NAMD_SOC._propagate_and_hop(self, self.prev_eval, eval_ha, t)
            self.v_all[self.qm_atoms] = self.vel
            if hopped:
                g_qm, e_diag, mult, state, w, pchg = self._soc_gradient_qmmm(u, self.active, eval_ha)
                f_all, epot = self._total_force_soc(potmm, g_qm, e_diag, pchg)
                accel_new = f_all / self.m_all[:, None]
                # Correct velocity rescaling for ESPF energy change at hop.
                # _propagate_and_hop accounts for ΔE_QM only (eval_ha gap). When the
                # ESPF density switches at an ISC hop the ESPF electrostatic energy
                # also changes by ΔE_ESPF = (epot_new - epot_old) - ΔE_QM. Apply an
                # additional isotropic rescaling to all atoms so total energy is
                # conserved. For ESPF_ROHF=1, charges are state-independent → 0.
                de_espf = ((epot_old - epot) +
                           (eval_ha[self.active - 1] - eval_ha[active_old - 1]))
                if abs(de_espf) > 1e-10:
                    ekin_all = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
                    if ekin_all > 0:
                        self.v_all *= np.sqrt(max(0.0, 1.0 + de_espf / ekin_all))

            accel = accel_new
            if self.econs:                                 # temporary E_tot-conservation rescale
                ke = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
                ket = self._e_ref_tot - epot
                if ke > 0 and ket > 0:
                    self.v_all *= np.sqrt(ket / ke)
            self._ulog = u
            NAMD_SOC._store_prev(self, self.r_all[self.qm_atoms].reshape((self.natom, 3)), u, eval_ha)
            self._log_soc_qmmm(istep, epot, mult, state, w, hopped)

        dump_log(mol, title='PyOQP: SOC-NAMD QM/MM trajectory complete')

    def _log_soc_qmmm(self, istep, epot, mult, state, w, hopped):
        ekin = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
        mch = self._ulog @ self.coef
        pmch = np.abs(mch) ** 2
        pop_s = float(pmch[:self.ns].sum())
        pop_t = float(pmch[self.ns:].sum())
        dump_log(
            self.mol,
            title=(f'SOC-QMMM-NAMD step {istep:6d}  t={(self._t_fs if self.dt_adaptive else istep*self.dt_fs):9.3f} fs  '
                   f'active={self.active}  E_tot={ekin+epot:.8f}  '
                   f'E_pot={epot:.8f}  E_kin={ekin:.8f}  hop={hopped}  '
                   f'dom=({NAMD_SOC._mch_label(mult, state)},w={w:.3f})  '
                   f'pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )


class NAMD_SOC_MCH_QMMM(NAMD_SOC_QMMM):
    """QM/MM SOC-NAMD in the MCH basis with exact active-root QM gradient."""

    _resolve_initial_mch_active = NAMD_SOC_MCH._resolve_initial_mch_active
    _mch_active_label = NAMD_SOC_MCH._mch_active_label
    _mch_propagate_and_hop = NAMD_SOC_MCH._mch_propagate_and_hop
    _log_mch = NAMD_SOC_MCH._log_mch

    def __init__(self, mol):
        super().__init__(mol)
        self.coef = np.zeros(self.nstate_soc, dtype=complex)
        self.coef[self.active - 1] = 1.0 + 0.0j

    def _mch_exact_gradient_qmmm(self, active):
        mol = self.mol
        mult, state = self._mch_target(active - 1)
        mol.data.set_tdhf_multiplicity(mult)
        SinglePoint(mol).excitation([self.e_ref])
        mol.config['properties']['grad'] = [state]
        Gradient(mol).gradient()
        g = np.array(mol.grads[state]).reshape((self.natom, 3))
        if os.environ.get('ESPF_ROHF', '').strip() in ('1', 'on'):
            oqp.form_esp_charges(mol)
            oqp.grad_esp_qmmm(mol)
        else:
            oqp.grad_esp_qmmm_excited(mol)
        g = g + np.array(mol.data["OQP::ESPF_GRAD"]).reshape((self.natom, 3))
        pchg = np.array(mol.data["OQP::partial_charges"]).copy()
        e = self._mch_energies_abs()[active - 1]
        return g, e, mult, state, pchg

    def run(self):
        mol = self.mol
        dump_log(mol, title='PyOQP: SOC-NAMD QM/MM (ISC, MCH-basis FSSH + ESPF embedding)')

        self._sync_positions()
        eval_ha, u, potmm, _ = self._electronic_soc_qmmm(with_overlap=False)
        self._resolve_initial_mch_active()
        h_mch = self._mch_hamiltonian_from_u(u, eval_ha)
        e_mch = self._mch_energies_abs()
        g_qm, e_pure, mult, state, pchg = self._mch_exact_gradient_qmmm(self.active)
        f_all, epot = self._total_force_soc(potmm, g_qm, e_pure, pchg)
        accel = f_all / self.m_all[:, None]
        self._rattle(self.r_all, self.v_all)
        self._thermalize_initial()
        r_qm = self.r_all[self.qm_atoms].reshape((self.natom, 3))
        NAMD_SOC._store_prev(self, r_qm, u, eval_ha)
        self._log_mch_qmmm(0, epot, mult, state, False)
        self._e_ref_tot = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2) + epot

        for istep in range(1, self.nstep + 1):
            self.dt = self._adaptive_dt(self.v_all, accel)
            self._t_fs += self.dt / FS_TO_AU
            r_old = self.r_all.copy()
            self.r_all = self.r_all + self.v_all * self.dt + 0.5 * accel * self.dt ** 2
            self._shake(r_old, self.r_all, self.v_all, self.dt)
            self._sync_positions()

            eval_ha, u, potmm, _ = self._electronic_soc_qmmm(with_overlap=False)
            h_mch = self._mch_hamiltonian_from_u(u, eval_ha)
            e_mch = self._mch_energies_abs()
            g_qm, e_pure, mult, state, pchg = self._mch_exact_gradient_qmmm(self.active)
            f_all, epot = self._total_force_soc(potmm, g_qm, e_pure, pchg)
            accel_new = f_all / self.m_all[:, None]
            self.v_all = self.v_all + 0.5 * (accel + accel_new) * self.dt
            self._rattle(self.r_all, self.v_all)

            active_old = self.active
            epot_old = epot
            self.vel = self.v_all[self.qm_atoms].copy()
            hopped = self._mch_propagate_and_hop(h_mch, e_mch)
            self.v_all[self.qm_atoms] = self.vel
            if hopped:
                g_qm, e_pure, mult, state, pchg = self._mch_exact_gradient_qmmm(self.active)
                f_all, epot = self._total_force_soc(potmm, g_qm, e_pure, pchg)
                accel_new = f_all / self.m_all[:, None]
                de_espf = ((epot_old - epot) +
                           (e_mch[self.active - 1] - e_mch[active_old - 1]))
                if abs(de_espf) > 1e-10:
                    ekin_all = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
                    if ekin_all > 0:
                        self.v_all *= np.sqrt(max(0.0, 1.0 + de_espf / ekin_all))

            accel = accel_new
            if self.econs:
                ke = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
                ket = self._e_ref_tot - epot
                if ke > 0 and ket > 0:
                    self.v_all *= np.sqrt(ket / ke)
            NAMD_SOC._store_prev(self, self.r_all[self.qm_atoms].reshape((self.natom, 3)), u, eval_ha)
            self._log_mch_qmmm(istep, epot, mult, state, hopped)

        dump_log(mol, title='PyOQP: SOC-MCH-QMMM-NAMD trajectory complete')

    def _log_mch_qmmm(self, istep, epot, mult, state, hopped):
        ekin = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
        pmch = np.abs(self.coef) ** 2
        pop_s = float(pmch[:self.ns].sum())
        pop_t = float(pmch[self.ns:].sum())
        dump_log(
            self.mol,
            title=(f'SOC-MCH-QMMM-NAMD step {istep:6d}  t={(self._t_fs if self.dt_adaptive else istep*self.dt_fs):9.3f} fs  '
                   f'active={self.active}:{self._mch_active_label(self.active)}  '
                   f'E_tot={ekin+epot:.8f}  E_pot={epot:.8f}  '
                   f'E_kin={ekin:.8f}  hop={hopped}  '
                   f'grad={NAMD_SOC._mch_label(mult, state)}  pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )
