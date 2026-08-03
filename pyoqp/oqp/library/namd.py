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
import hashlib
import json
import re
import struct
import tempfile
import numpy as np

import oqp
from oqp.library.single_point import SinglePoint, Gradient, LastStep, BasisOverlap, NACME
from oqp.library.nac_utils import canonical_state_overlap
from oqp.library.odp import odp_from_config
from oqp.utils.tb_backends import is_tb_method, make_tb_adapter, tb_section_name
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
INT64_MIN = -(1 << 63)
INT64_MAX = (1 << 63) - 1
NAMD_RESTART_SCHEMA_VERSION = 3
NAMD_TRAJECTORY_SCHEMA_VERSION = 3
NAMD_TRAJECTORY_MAGIC = b'OQPNTRJ1'


def _restart_identity_digest(array_parts=(), text_parts=()):
    """Return a stable digest for static molecular and topology identity."""
    digest = hashlib.sha256()
    for label, value, dtype in array_parts:
        name = str(label).encode('utf-8')
        array = np.ascontiguousarray(np.asarray(value, dtype=dtype))
        digest.update(struct.pack('<Q', len(name)))
        digest.update(name)
        shape = json.dumps(array.shape, separators=(',', ':')).encode('ascii')
        digest.update(struct.pack('<Q', len(shape)))
        digest.update(shape)
        digest.update(array.dtype.str.encode('ascii'))
        digest.update(array.tobytes(order='C'))
    for label, value in text_parts:
        name = str(label).encode('utf-8')
        payload = str(value).encode('utf-8')
        digest.update(struct.pack('<Q', len(name)))
        digest.update(name)
        digest.update(struct.pack('<Q', len(payload)))
        digest.update(payload)
    return digest.hexdigest()


def _validate_gate_tolerances(label, values):
    """Reject NaN/Inf as well as negative validation tolerances."""
    if not all(np.isfinite(value) and value >= 0.0 for value in values):
        raise ValueError(
            f"[md] {label} gate tolerances must be finite and non-negative")


def _validate_distinct_output_paths(**paths):
    """Reject NAMD outputs that would overwrite or corrupt one another."""
    aliases = {}
    for label, path in paths.items():
        normalized = os.path.normcase(os.path.realpath(os.fspath(path)))
        aliases.setdefault(normalized, []).append(label)
    collisions = [labels for labels in aliases.values() if len(labels) > 1]
    if collisions:
        rendered = "; ".join(" = ".join(labels) for labels in collisions)
        raise ValueError(
            "[md] NAMD output paths must be distinct; collision: " + rendered)


def _namd_trajectory_dtype(nstate, natom, ncv=0):
    """Fixed-width, appendable binary record used by ``*.namd.trj``."""
    matrix = (nstate, nstate)
    vectors = (natom, 3)
    return np.dtype([
        ('step', '<i8'), ('time_fs', '<f8'), ('active', '<i4'),
        ('hopped', 'i1'), ('rng', '<f8'),
        ('e_unbiased_pot_hartree', '<f8'), ('e_pot_hartree', '<f8'),
        ('e_kin_hartree', '<f8'),
        ('e_tot_hartree', '<f8'), ('state_energies', '<f8', (nstate,)),
        ('populations', '<f8', (nstate,)), ('coef_real', '<f8', (nstate,)),
        ('coef_imag', '<f8', (nstate,)), ('coordinates_bohr', '<f8', vectors),
        ('velocities_au', '<f8', vectors), ('state_overlap', '<f8', matrix),
        ('overlap_tdc_au', '<f8', matrix), ('reference_tdc_au', '<f8', matrix),
        ('reference_mask', 'u1', matrix), ('reference_source', 'i1'),
        ('gate_center_step', '<i8'), ('gate_verdict', 'i1'),
        ('gate_counts', '<i8', (3,)), ('gate_metrics', '<f8', (7,)),
        ('nve_verdict', 'i1'), ('nve_streak', '<i8'),
        ('nve_metrics', '<f8', (4,)),
        ('odp_window', '<i8'), ('odp_xi', '<f8'),
        ('odp_cv_raw', '<f8', (ncv,)), ('odp_cv_scaled', '<f8', (ncv,)),
        ('odp_cv_perpendicular', '<f8', (ncv,)),
        ('odp_perpendicular_norm', '<f8'),
        ('odp_bias_parallel_hartree', '<f8'),
        ('odp_bias_perpendicular_hartree', '<f8'),
        ('odp_bias_hartree', '<f8'),
        ('tracking_valid', 'i1'), ('tracking_order', '<i8', (nstate,)),
        ('tracking_phase', '<f8', (nstate,)),
        ('tracking_overlap', '<f8', (nstate,)),
        ('tracking_margin', '<f8', (nstate,)),
    ], align=False)


def read_namd_trajectory(path, mmap_mode='r'):
    """Return ``(header, records)`` for a dense OpenQP NAMD trajectory.

    ``records`` is a NumPy memmap by default, so even a long ensemble can be
    analysed without loading all coordinates and coupling matrices into RAM.
    """
    with open(path, 'rb') as stream:
        if stream.read(8) != NAMD_TRAJECTORY_MAGIC:
            raise ValueError(f'not an OpenQP dense NAMD trajectory: {path}')
        header_size = struct.unpack('<Q', stream.read(8))[0]
        header = json.loads(stream.read(header_size).decode('utf-8'))
        offset = 16 + header_size
    if int(header.get('schema_version', -1)) != NAMD_TRAJECTORY_SCHEMA_VERSION:
        raise ValueError('unsupported OpenQP NAMD trajectory schema')
    dtype = _namd_trajectory_dtype(
        int(header['nstate']), int(header['natom']), int(header.get('ncv', 0)))
    payload_size = os.path.getsize(path) - offset
    if payload_size < 0 or payload_size % dtype.itemsize:
        raise ValueError('truncated OpenQP NAMD trajectory record')
    count = payload_size // dtype.itemsize
    records = np.memmap(path, dtype=dtype, mode=mmap_mode,
                        offset=offset, shape=(count,))
    return header, records


def read_odp_wham_series(path):
    """Extract the complete per-window ODP series from one packed NAMD TRJ."""
    header, records = read_namd_trajectory(path)
    provenance = header.get('odp', {})
    if not provenance.get('enabled', False):
        raise ValueError('trajectory does not contain an enabled ODP umbrella')
    required = {
        'window', 'cv', 'cv_atom_indexing', 'cv_native_units', 'scale',
        'reference_r', 'reference_p', 'scaled_path_length', 'center',
        'k_parallel_hartree', 'k_perpendicular_hartree', 'projection',
        'perpendicular_restraint',
    }
    missing = sorted(required.difference(provenance))
    if missing:
        raise ValueError(
            f'trajectory ODP provenance is missing {", ".join(missing)}')
    ncv = int(header.get('ncv', 0))
    scale = np.asarray(provenance['scale'], dtype=np.float64).reshape(-1)
    reference_r = np.asarray(
        provenance['reference_r'], dtype=np.float64).reshape(-1)
    reference_p = np.asarray(
        provenance['reference_p'], dtype=np.float64).reshape(-1)
    center = float(provenance['center'])
    k_parallel = float(provenance['k_parallel_hartree'])
    k_perpendicular = float(provenance['k_perpendicular_hartree'])
    if not (ncv > 0 and len(provenance['cv']) == ncv
            and len(provenance['cv_native_units']) == ncv
            and scale.size == reference_r.size == reference_p.size == ncv
            and np.all(np.isfinite(scale)) and np.all(scale > 0.0)
            and np.all(np.isfinite(reference_r))
            and np.all(np.isfinite(reference_p))
            and np.isfinite(center) and np.isfinite(k_parallel)
            and np.isfinite(k_perpendicular) and k_parallel > 0.0
            and k_perpendicular >= 0.0
            and provenance['cv_atom_indexing'] == '1-based'
            and provenance['projection']
            == 'signed_scaled_dot_over_path_norm_squared'
            and bool(provenance['perpendicular_restraint'])
            == (k_perpendicular > 0.0)):
        raise ValueError('trajectory ODP provenance has an invalid CV metric')
    expected_window = int(provenance['window'])
    recorded_windows = np.asarray(records['odp_window'], dtype=np.int64)
    if np.any(recorded_windows != expected_window):
        raise ValueError(
            'trajectory ODP window records disagree with header provenance')
    raw = np.asarray(records['odp_cv_raw'], dtype=np.float64)
    scaled = np.asarray(records['odp_cv_scaled'], dtype=np.float64)
    xi = np.asarray(records['odp_xi'], dtype=np.float64)
    perpendicular = np.asarray(
        records['odp_cv_perpendicular'], dtype=np.float64)
    perpendicular_norm = np.asarray(
        records['odp_perpendicular_norm'], dtype=np.float64)
    direction = scale*(reference_p - reference_r)
    path_length_squared = float(np.dot(direction, direction))
    if not np.isfinite(path_length_squared) or path_length_squared <= 1.0e-24:
        raise ValueError('trajectory ODP provenance has a degenerate R/P path')
    if not np.isclose(
            float(provenance['scaled_path_length']),
            np.sqrt(path_length_squared), rtol=2.0e-12, atol=2.0e-12):
        raise ValueError(
            'trajectory ODP scaled path length disagrees with its metric')
    displacement = scaled - scale*reference_r
    expected_xi = displacement @ direction/path_length_squared
    expected_perpendicular = displacement - expected_xi[:, None]*direction
    expected_perpendicular_norm = np.linalg.norm(
        expected_perpendicular, axis=1)
    expected_bias = (
        0.5*k_parallel*(xi - center)**2
        + 0.5*k_perpendicular*perpendicular_norm**2
    )
    expected_bias_parallel = 0.5*k_parallel*(xi - center)**2
    expected_bias_perpendicular = 0.5*k_perpendicular*perpendicular_norm**2
    checks = (
        (scaled, raw*scale, 'scaled CV'),
        (xi, expected_xi, 'signed projection'),
        (perpendicular, expected_perpendicular, 'perpendicular CV'),
        (perpendicular_norm, expected_perpendicular_norm,
         'perpendicular norm'),
        (np.asarray(records['odp_bias_parallel_hartree'], dtype=np.float64),
         expected_bias_parallel, 'parallel bias energy'),
        (np.asarray(
            records['odp_bias_perpendicular_hartree'], dtype=np.float64),
         expected_bias_perpendicular, 'perpendicular bias energy'),
        (np.asarray(records['odp_bias_hartree'], dtype=np.float64),
         expected_bias, 'bias energy'),
    )
    for actual, expected, label in checks:
        if (not np.all(np.isfinite(actual))
                or not np.allclose(actual, expected, rtol=2.0e-12,
                                   atol=2.0e-12)):
            raise ValueError(
                f'trajectory ODP {label} records disagree with provenance')
    result = {
        'provenance': provenance,
        'ensemble': header.get('ensemble'),
        'step': np.array(records['step'], copy=True),
        'time_fs': np.array(records['time_fs'], copy=True),
        'window': np.array(records['odp_window'], copy=True),
        'xi': np.array(xi, copy=True),
        'cv_raw': np.array(raw, copy=True),
        'cv_scaled': np.array(scaled, copy=True),
        'perpendicular_norm': np.array(perpendicular_norm, copy=True),
        'bias_hartree': np.array(records['odp_bias_hartree'], copy=True),
        'unbiased_potential_hartree': np.array(
            records['e_unbiased_pot_hartree'], copy=True),
        'total_conservative_hartree': np.array(
            records['e_tot_hartree'], copy=True),
    }
    del records
    return result

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


def _parse_soc_init_state(label, ns, nt, *, public_labels=False):
    """Return ``(multiplicity, zero_based_root, public_label)`` for SOC-NAMD.

    Canonical ``.oqp`` inputs use the public MRSF labels ``S0/T0`` for the
    lowest roots.  Historical ``[md] init_state`` values, however, numbered
    triplets as ``T1, T2, ...``.  Keep that legacy spelling working while also
    accepting ``T0`` as an unambiguous alias for the first legacy triplet.
    Bounds are checked here so an invalid label never reaches a NumPy index.
    """
    text = str(label or "").strip().upper()
    match = re.fullmatch(r"([ST])(\d+)", text)
    if match is None:
        raise ValueError(
            "[md] init_state must be an MRSF state label such as S0, S1, T0, or T1"
        )

    manifold, requested = match.group(1), int(match.group(2))
    if manifold == "S":
        mult = 1
        root = requested
        count = int(ns)
    else:
        mult = 3
        count = int(nt)
        # .oqp state names are uniformly zero based.  Legacy INI/API decks
        # historically used T1 for the first triplet; T0 is accepted there as
        # a transition-friendly alias rather than producing a negative index.
        root = requested if public_labels else max(0, requested - 1)

    if root < 0 or root >= count:
        last = max(0, count - 1)
        if manifold == "T" and not public_labels:
            valid = "T1" if count == 1 else "T1-T%d" % count
            valid += " (T0 is also accepted for the first root)"
        else:
            valid = "%s0" % manifold if count == 1 else "%s0-%s%d" % (
                manifold, manifold, last)
        raise ValueError(
            "[md] init_state='%s' is outside the SOC MCH basis; available %s states: %s"
            % (text, manifold, valid)
        )

    return mult, root, "%s%d" % (manifold, root)


def _select_response_manifold(mol, multiplicity):
    """Synchronize native, config, DFTB, and log-facing target spin state."""
    mult = int(multiplicity)
    mol.config['tdhf']['multiplicity'] = mult
    mol.data.set_tdhf_multiplicity(mult)
    # Both tight-binding backends (dftb and xtb) mirror the target multiplicity
    # into their own [dftb]/[xtb] section so the native library selects the
    # matching MRSF response manifold.
    if is_tb_method(mol.config['input']['method']):
        mol.config[tb_section_name(mol.config)]['target_multiplicity'] = mult


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
        self.rng_stream = int(md.get('rng_stream', 0))
        self.first_hop_step = int(md.get('first_hop_step', 2))
        self.nacme_check = str(md.get('nacme_check', 'off')).strip().lower().replace('-', '_')
        if self.nacme_check == 'tdba':
            self.nacme_check = 'baeck_an'
        if self.nacme_check not in ('off', 'baeck_an'):
            raise ValueError("[md] nacme_check must be off or baeck_an")
        self.ba_gap_max = float(md.get('ba_gap_max', 0.0734986443513))
        if not np.isfinite(self.ba_gap_max) or self.ba_gap_max <= 0.0:
            raise ValueError("[md] ba_gap_max must be positive and finite")
        self.nacme_gate = str(md.get('nacme_gate', 'warn')).strip().lower()
        if self.nacme_gate not in ('off', 'warn', 'error'):
            raise ValueError("[md] nacme_gate must be off, warn, or error")
        self.nacme_gate_invariant_tol = float(
            md.get('nacme_gate_invariant_tol', 1.0e-10))
        self.nacme_gate_abs_tol = float(md.get('nacme_gate_abs_tol', 1.0e-4))
        self.nacme_gate_rel_tol = float(md.get('nacme_gate_rel_tol', 1.0))
        self.nacme_gate_consecutive = int(md.get('nacme_gate_consecutive', 3))
        _validate_gate_tolerances('NACME', (
            self.nacme_gate_invariant_tol, self.nacme_gate_abs_tol,
            self.nacme_gate_rel_tol,
        ))
        if self.nacme_gate_consecutive < 1:
            raise ValueError("[md] nacme_gate_consecutive must be at least 1")
        self.nve_gate = str(md.get('nve_gate', 'off')).strip().lower()
        if self.nve_gate not in ('off', 'warn', 'error'):
            raise ValueError("[md] nve_gate must be off, warn, or error")
        self.nve_gate_abs_tol = float(md.get('nve_gate_abs_tol', 5.0e-3))
        self.nve_gate_step_tol = float(md.get('nve_gate_step_tol', 1.0e-3))
        self.nve_gate_transition_tol = float(
            md.get('nve_gate_transition_tol', 1.0e-6))
        self.nve_gate_consecutive = int(md.get('nve_gate_consecutive', 3))
        _validate_gate_tolerances('NVE', (
            self.nve_gate_abs_tol, self.nve_gate_step_tol,
            self.nve_gate_transition_tol,
        ))
        if self.nve_gate_consecutive < 1:
            raise ValueError("[md] nve_gate_consecutive must be at least 1")
        self.trajectory_interval = int(md.get('trajectory_interval', 1))
        self.restart_interval = int(md.get('restart_interval', 1))
        if self.trajectory_interval < 1 or self.restart_interval < 1:
            raise ValueError("[md] trajectory/restart intervals must be at least 1")
        self.restart_requested = self._as_bool(md.get('restart', False))
        self.trajectory_file = self._md_output_path(
            md.get('trajectory_file', ''), '.namd.trj')
        self.nacme_audit_file = self._md_output_path(
            md.get('nacme_audit_file', ''), '.namd.nacme.tsv')
        self.restart_file = self._md_output_path(
            md.get('restart_file', ''), '.namd.restart.npz')
        self.restart_manifest_file = os.path.join(
            os.path.dirname(os.path.abspath(self.mol.log)), 'restart.oqp')
        _validate_distinct_output_paths(
            log_file=self.mol.log,
            trajectory_file=self.trajectory_file,
            nacme_audit_file=self.nacme_audit_file,
            restart_file=self.restart_file,
            restart_manifest_file=self.restart_manifest_file,
        )
        self.odp = odp_from_config(cfg)
        self._odp_last = None
        self._unbiased_potential_energy = np.nan
        _soc = md.get('soc', False)
        soc_requested = (_soc is True) or (str(_soc).lower() in ('true', '1', 'on', 'yes'))
        if soc_requested and self.nacme_check != 'off':
            raise NotImplementedError(
                "[md] nacme_check currently supports same-spin NAMD only"
            )
        if soc_requested and self.restart_requested:
            raise NotImplementedError(
                "[md] restart currently supports same-spin NAMD only"
            )
        if soc_requested and self.nve_gate != 'off':
            raise NotImplementedError(
                "[md] nve_gate currently supports same-spin NAMD only"
            )
        if soc_requested and self.odp is not None:
            raise NotImplementedError(
                "[odp] currently supports same-spin NVE NAMD only"
            )
        if not INT64_MIN <= self.seed <= INT64_MAX:
            raise ValueError("[md] seed must fit in a signed 64-bit integer")
        if not 0 <= self.rng_stream <= INT64_MAX:
            raise ValueError(
                "[md] rng_stream must be a non-negative signed 64-bit integer"
            )
        if self.first_hop_step < 1:
            raise ValueError("[md] first_hop_step must be at least 1")
        self.velocity_source = str(md['velocity'])

        self.natom = mol.data['natom']
        if self.odp is not None and not self._as_bool(
                cfg.get('input', {}).get('qmmm_flag', False)):
            self.odp.validate_atom_count(self.natom)
        # get_mass() returns atomic masses in amu; the integrator works in
        # atomic units, so convert to electron masses.
        self.mass = mol.get_mass() * AMU_TO_AU     # (natom,) electron masses
        self._restart_system_identity = self._qm_restart_system_identity()
        self._rng_step = 0
        self._last_hop_random = np.nan
        self._ba_energy_left = None
        self._ba_energy_center = None
        self._ba_tdc_left = None
        self._ba_dt_left = None
        self._ba_last = None
        self._nacme_gate_failures = 0
        self._nacme_gate_last = None
        self._nacme_reference_tdc = None
        self._nacme_reference_mask = None
        self._nacme_reference_source = 0
        self._last_state_overlap = None
        self._last_overlap_tdc = None
        self._nve_reference_energy = None
        self._nve_previous_energy = None
        self._nve_gate_failures = 0
        self._nve_gate_last = None
        self._pending_nve_gate_error = None

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

    @staticmethod
    def _as_bool(value):
        return (value is True) or (str(value).lower() in ('true', '1', 'on', 'yes'))

    def _qm_restart_system_identity(self):
        """Bind restarts to the ordered atoms, masses, and starting geometry."""
        digest = _restart_identity_digest(array_parts=(
            ('atomic_numbers', self.mol.get_atoms(), '<i8'),
            ('masses_electron', self.mass, '<f8'),
            ('initial_coordinates_bohr', self.mol.get_system(), '<f8'),
        ))
        return {'kind': 'qm', 'natom': int(self.natom), 'sha256': digest}

    def _md_output_path(self, configured, suffix):
        """Resolve NAMD sidecars beside the main log, not the process CWD."""
        value = str(configured or '').strip()
        log_dir = os.path.dirname(os.path.abspath(self.mol.log))
        if value:
            return value if os.path.isabs(value) else os.path.join(log_dir, value)
        stem = os.path.splitext(os.path.basename(self.mol.log))[0]
        return os.path.join(log_dir, stem + suffix)

    def _is_io_rank(self):
        manager = getattr(self.mol, 'mpi_manager', None)
        return manager is None or int(getattr(manager, 'rank', 0)) == 0

    def _io_barrier(self):
        manager = getattr(self.mol, 'mpi_manager', None)
        if manager is not None and hasattr(manager, 'barrier'):
            manager.barrier()

    def _prepare_md_outputs(self):
        """Start fresh sidecars or preserve them when explicitly restarting."""
        if self._is_io_rank() and not self.restart_requested:
            for path in (self.trajectory_file, self.nacme_audit_file):
                with open(path, 'w', encoding='utf-8'):
                    pass
        self._io_barrier()
        dump_log(
            self.mol,
            title=(f'NAMD files: trajectory={self.trajectory_file} '
                   f'nacme_audit={self.nacme_audit_file} '
                   f'restart={self.restart_file} '
                   f'manifest={self.restart_manifest_file}'),
        )
        if self.odp is not None:
            dump_log(
                self.mol,
                title=('ODP umbrella: '
                       f'window={self.odp.window} center={self.odp.center:g} '
                       f'k_parallel={self.odp.k_parallel:g} Ha '
                       f'k_perpendicular={self.odp.k_perpendicular:g} Ha '
                       f'CVs={"; ".join(self.odp.cv_labels)}'),
            )

    def _odp_provenance(self):
        if getattr(self, 'odp', None) is None:
            return {'enabled': False}
        return self.odp.provenance()

    def _evaluate_odp(self, coordinates):
        if self.odp is None:
            self._odp_last = None
            return None
        self._odp_last = self.odp.evaluate(coordinates)
        return self._odp_last

    def _prepare_hop_step(self, istep):
        """Bind the physical MD step to the stateless hop RNG.

        Returning ``False`` suppresses both electronic propagation and the
        stochastic FSSH decision.  The default first_hop_step=2 reproduces the
        KNU-GAMESS/TLF2 initialisation convention without consuming a random
        value at the skipped first interval.
        """
        self._rng_step = int(istep)
        self._last_hop_random = np.nan
        return self._rng_step >= self.first_hop_step

    def _hop_random(self):
        """Return the Fortran counter-RNG value for this trajectory and step."""
        if np.isfinite(self._last_hop_random):
            return self._last_hop_random
        override = getattr(self, '_hop_random_override', None)
        if override is not None:
            value = float(override.random())
        else:
            value = float(oqp.oqp_namd_counter_random(
                self.seed, self.rng_stream, self._rng_step))
        if not 0.0 <= value < 1.0:
            raise RuntimeError(
                "NAMD hop random value must lie in [0,1); got "
                f"{value!r} at step {self._rng_step}"
            )
        self._last_hop_random = value
        return value

    def _hop_rng_log(self):
        value = self._last_hop_random
        rendered = "skipped" if not np.isfinite(value) else f"{value:.17g}"
        return f"rng_step={self._rng_step} random={rendered}"

    # ------------------------------------------------------------------ #
    # initialisation
    # ------------------------------------------------------------------ #
    def _counter_normals(self, shape):
        """Return resident-Fortran standard normals for this trajectory."""
        size = int(np.prod(shape, dtype=np.int64))
        normals = np.empty(size, dtype=np.float64)
        oqp.oqp_namd_counter_normal_fill(
            self.seed,
            self.rng_stream,
            size,
            oqp.ffi.cast("double *", normals.ctypes.data),
        )
        return normals.reshape(shape)

    def _init_velocities(self):
        src = self.velocity_source.lower()
        if src in ('zero', 'none', '0'):
            return np.zeros((self.natom, 3))
        if src in ('maxwell', 'boltzmann', 'random'):
            sigma = np.sqrt(KB_HARTREE * self.init_temp / self.mass)  # (natom,)
            v = self._counter_normals((self.natom, 3)) * sigma[:, None]
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
        gradient = np.array(mol.grads[self.active]).reshape((self.natom, 3))
        odp = self._evaluate_odp(mol.get_system().reshape((self.natom, 3)))
        if odp is not None:
            gradient = gradient - odp['force']
        return gradient

    def _state_overlap(self, istep=None):
        """Compute the phase-corrected state overlap S(i,j)=<i(t-dt)|j(t)>."""
        NACME(self.mol).nacme()
        state_overlap = canonical_state_overlap(
            self.mol.data["OQP::td_states_overlap"]
        )
        self._last_state_overlap = np.array(state_overlap, copy=True)
        self._last_overlap_tdc = np.array(self._compute_tdc(state_overlap), copy=True)
        self._update_baeck_an_check(istep, state_overlap)
        return state_overlap

    def _update_baeck_an_check(self, istep, state_overlap):
        """Compare overlap TDC magnitudes with a centred TD-Baeck-An estimate.

        TD-BA is phase-free and therefore cannot validate the signed gauge.
        It is retained only as an independent energy-curvature diagnostic.
        """
        if self.nacme_check != 'baeck_an':
            return

        n = self.nstate
        data = self.mol.data
        energies_old = np.ascontiguousarray(
            np.asarray(data["OQP::td_energies_old"], dtype=np.float64).reshape(-1)[:n]
        )
        energies_current = np.ascontiguousarray(
            np.asarray(data["OQP::td_energies"], dtype=np.float64).reshape(-1)[:n]
        )
        tdc_current = np.ascontiguousarray(
            self._compute_tdc(state_overlap), dtype=np.float64
        )
        dt_right = float(self.dt)

        if self._ba_energy_center is None:
            self._ba_energy_left = energies_old.copy()
            self._ba_energy_center = energies_current.copy()
            self._ba_tdc_left = tdc_current.copy()
            self._ba_dt_left = dt_right
            return

        if not np.allclose(
                energies_old, self._ba_energy_center, rtol=0.0, atol=1.0e-12):
            dump_log(
                self.mol,
                title='NACME check: Baeck-An history discontinuity; reseeding',
            )
            self._ba_energy_left = energies_old.copy()
            self._ba_energy_center = energies_current.copy()
            self._ba_tdc_left = tdc_current.copy()
            self._ba_dt_left = dt_right
            self._reset_nacme_gate_state()
            return

        ba_tdc = np.zeros((n, n), dtype=np.float64)
        status = oqp.oqp_namd_baeck_an_tdc(
            n,
            self._ba_dt_left,
            dt_right,
            self.ba_gap_max,
            oqp.ffi.cast("double *", self._ba_energy_left.ctypes.data),
            oqp.ffi.cast("double *", self._ba_energy_center.ctypes.data),
            oqp.ffi.cast("double *", energies_current.ctypes.data),
            oqp.ffi.cast("double *", ba_tdc.ctypes.data),
        )
        if status != 0:
            raise RuntimeError(f"native Baeck-An NACME check failed (status={status})")

        dt_sum = self._ba_dt_left + dt_right
        overlap_center = (
            dt_right*self._ba_tdc_left + self._ba_dt_left*tdc_current
        )/dt_sum
        center_step = None if istep is None else int(istep) - 1
        ba_mask = np.asarray(np.abs(ba_tdc) > 0.0, dtype=np.int32)
        gate = self._run_nacme_gate(
            overlap_center,
            ba_tdc,
            reference_mask=ba_mask,
            source='TD-Baeck-An',
            center_step=center_step,
            signed=False,
        )
        self._ba_last = {
            'center_step': center_step,
            'baeck_an_tdc': ba_tdc.copy(),
            'overlap_tdc_centered': overlap_center.copy(),
            'magnitude_rms_error': gate['pair_rms_error'],
            'magnitude_max_error': gate['pair_max_error'],
            'gate': gate,
        }
        dump_log(
            self.mol,
            title='NACME check: TD-Baeck-An TDC (magnitude diagnostic)',
            section='nacm',
            info=ba_tdc,
        )
        dump_log(
            self.mol,
            title='NACME check: centered overlap TDC',
            section='nacm',
            info=overlap_center,
        )

        self._ba_energy_left = self._ba_energy_center.copy()
        self._ba_energy_center = energies_current.copy()
        self._ba_tdc_left = tdc_current.copy()
        self._ba_dt_left = dt_right

    def _reset_nacme_gate_state(self):
        """Discard gate evidence that cannot cross a history discontinuity."""
        self._nacme_gate_last = None
        self._nacme_reference_tdc = None
        self._nacme_reference_mask = None
        self._nacme_reference_source = 0
        self._nacme_gate_failures = 0

    def _run_nacme_gate(self, candidate_tdc, reference_tdc, *,
                        reference_mask=None, source='reference',
                        center_step=None, signed=False):
        """Run the common resident-Fortran NACME validation gate.

        Future analytic NAC support should contract the phase-aligned analytic
        coupling vector with the nuclear velocity at the same time point, then
        call this method with ``signed=True``.  TD-Baeck-An calls it with
        ``signed=False`` because an energy-only estimate has no wavefunction
        gauge.  Thus the invariant and policy machinery is shared without
        treating the approximate TD-BA sign as physical.
        """
        n = self.nstate
        candidate = np.ascontiguousarray(
            np.asarray(candidate_tdc, dtype=np.float64).reshape((n, n)))
        reference = np.ascontiguousarray(
            np.asarray(reference_tdc, dtype=np.float64).reshape((n, n)))
        if reference_mask is None:
            mask = np.ones((n, n), dtype=np.int32)
            np.fill_diagonal(mask, 0)
        else:
            mask = np.ascontiguousarray(
                np.asarray(reference_mask, dtype=np.int32).reshape((n, n)))

        metrics = np.zeros(7, dtype=np.float64)
        counts = np.zeros(3, dtype=np.int64)
        status = oqp.oqp_namd_nacme_gate(
            n,
            oqp.ffi.cast("double *", candidate.ctypes.data),
            oqp.ffi.cast("double *", reference.ctypes.data),
            oqp.ffi.cast("int *", mask.ctypes.data),
            int(bool(signed)),
            self.nacme_gate_invariant_tol,
            self.nacme_gate_abs_tol,
            self.nacme_gate_rel_tol,
            oqp.ffi.cast("double *", metrics.ctypes.data),
            oqp.ffi.cast("int64_t *", counts.ctypes.data),
        )
        if status != 0:
            raise RuntimeError(
                f"native NACME validation gate failed for {source} (status={status})"
            )

        compared_pairs = int(counts[0])
        invariant_failures = int(counts[1])
        reference_failures = int(counts[2])
        if reference_failures:
            self._nacme_gate_failures += 1
        else:
            self._nacme_gate_failures = 0
        if invariant_failures or reference_failures:
            verdict = 'fail'
        elif compared_pairs == 0:
            verdict = 'not_evaluable'
        else:
            verdict = 'pass'

        result = {
            'source': source,
            'center_step': center_step,
            'signed_comparison': bool(signed),
            'verdict': verdict,
            'compared_pairs': compared_pairs,
            'invariant_failures': invariant_failures,
            'reference_failures': reference_failures,
            'consecutive_reference_failures': self._nacme_gate_failures,
            'candidate_diagonal_max': float(metrics[0]),
            'candidate_antisymmetry_max': float(metrics[1]),
            'reference_diagonal_max': float(metrics[2]),
            'reference_antisymmetry_max': float(metrics[3]),
            'pair_rms_error': float(metrics[4]),
            'pair_max_error': float(metrics[5]),
            'max_tolerance_ratio': float(metrics[6]),
        }
        self._nacme_gate_last = result
        self._nacme_reference_tdc = reference.copy()
        self._nacme_reference_mask = mask.copy()
        self._nacme_reference_source = {
            'TD-Baeck-An': 1,
            'analytic': 2,
        }.get(source, 127)
        table = (
            "   center  source          verdict          pairs  inv  ref  "
            "diag_max      anti_max      rms_error     max_error     ratio   streak\n"
            "   ------  --------------  ---------------  -----  ---  ---  "
            "------------  ------------  ------------  ------------  ------  ------\n"
            f"   {str(center_step):>6}  {source[:14]:<14}  {verdict:<15}  "
            f"{compared_pairs:5d}  {invariant_failures:3d}  {reference_failures:3d}  "
            f"{metrics[0]:12.4e}  {metrics[1]:12.4e}  {metrics[4]:12.4e}  "
            f"{metrics[5]:12.4e}  {metrics[6]:6.2f}  "
            f"{self._nacme_gate_failures:6d}"
        )
        dump_log(
            self.mol,
            title='NACME validation gate',
            section='text',
            info={'text': table},
        )
        self._write_nacme_audit_row(result)

        if self.nacme_gate == 'error':
            if invariant_failures:
                raise RuntimeError(
                    f"NACME invariant gate failed for {source} at step {center_step}"
                )
            if self._nacme_gate_failures >= self.nacme_gate_consecutive:
                raise RuntimeError(
                    f"NACME reference gate failed {self._nacme_gate_failures} "
                    f"consecutive times for {source} at step {center_step}"
                )
        return result

    def _write_nacme_audit_row(self, result):
        """Append one machine-readable gate row for ensemble/post-MD audits."""
        if not self._is_io_rank():
            self._io_barrier()
            return
        path = self.nacme_audit_file
        needs_header = not os.path.exists(path) or os.path.getsize(path) == 0
        columns = (
            'center_step', 'source', 'verdict', 'signed_comparison',
            'compared_pairs',
            'invariant_failures', 'reference_failures',
            'consecutive_reference_failures', 'candidate_diagonal_max',
            'candidate_antisymmetry_max', 'reference_diagonal_max',
            'reference_antisymmetry_max', 'pair_rms_error', 'pair_max_error',
            'max_tolerance_ratio',
        )
        values = [result.get(name, '') for name in columns]
        with open(path, 'a', encoding='utf-8') as stream:
            if needs_header:
                stream.write('\t'.join(columns) + '\n')
            stream.write('\t'.join(str(value) for value in values) + '\n')
            stream.flush()
            os.fsync(stream.fileno())
        self._io_barrier()

    def _update_nve_gate(self, istep, epot, ekin, transition_energy_jump=np.nan):
        """Audit microcanonical energy conservation for same-spin FSSH."""
        self._pending_nve_gate_error = None
        total = float(epot + ekin)
        if not np.isfinite(total):
            drift = np.nan
            step_change = np.nan
        elif self._nve_reference_energy is None:
            self._nve_reference_energy = total
            self._nve_previous_energy = total
            drift = 0.0
            step_change = 0.0
        else:
            drift = total - self._nve_reference_energy
            step_change = total - self._nve_previous_energy
        transition_jump = float(transition_energy_jump)
        time_fs = self._t_fs if self.dt_adaptive else istep*self.dt_fs
        drift_rate = drift/time_fs if time_fs > 0.0 else 0.0
        transition_failure = (
            np.isfinite(transition_jump)
            and abs(transition_jump) > self.nve_gate_transition_tol
        )
        finite_failure = not np.isfinite(total)
        drift_failure = finite_failure or abs(drift) > self.nve_gate_abs_tol
        step_failure = finite_failure or (
            istep > 0 and abs(step_change) > self.nve_gate_step_tol)
        failed = finite_failure or transition_failure or drift_failure or step_failure
        if failed:
            self._nve_gate_failures += 1
        else:
            self._nve_gate_failures = 0
        verdict = 'off' if self.nve_gate == 'off' else ('fail' if failed else 'pass')
        result = {
            'step': int(istep), 'verdict': verdict,
            'total_energy': total, 'drift': drift, 'step_change': step_change,
            'transition_energy_jump': transition_jump, 'drift_rate': drift_rate,
            'drift_failure': bool(drift_failure),
            'step_failure': bool(step_failure),
            'transition_failure': bool(transition_failure),
            'consecutive_failures': self._nve_gate_failures,
        }
        self._nve_gate_last = result
        if np.isfinite(total):
            self._nve_previous_energy = total

        if self.nve_gate != 'off':
            table = (
                "   step  verdict  E_total(Ha)       drift(Ha)       step_dE(Ha)     "
                "transition_dE    drift(Ha/fs)   streak\n"
                "   ----  -------  ----------------  --------------  --------------  "
                "---------------  -------------  ------\n"
                f"   {istep:4d}  {verdict:<7}  {total:16.9f}  {drift:14.6e}  "
                f"{step_change:14.6e}  {transition_jump:15.6e}  "
                f"{drift_rate:13.5e}  {self._nve_gate_failures:6d}"
            )
            dump_log(
                self.mol, title='NVE energy validation gate', section='text',
                info={'text': table},
            )
        if self.nve_gate == 'error':
            if finite_failure:
                self._pending_nve_gate_error = RuntimeError(
                    f'non-finite NAMD total energy at step {istep}'
                )
            elif transition_failure:
                self._pending_nve_gate_error = RuntimeError(
                    f'NVE transition-energy gate failed at step {istep}: '
                    f'{transition_jump:.8e} Ha'
                )
            elif self._nve_gate_failures >= self.nve_gate_consecutive:
                self._pending_nve_gate_error = RuntimeError(
                    f'NVE energy gate failed {self._nve_gate_failures} '
                    f'consecutive times at step {istep}'
                )
        return result

    def _enforce_nve_gate(self):
        """Stop only after the failing point has been written to the dense TRJ."""
        error = self._pending_nve_gate_error
        self._pending_nve_gate_error = None
        if error is not None:
            raise error

    def _write_md_trajectory(self, istep, coordinates, epot, ekin, hopped):
        """Append one lossless, fixed-width record to the dense binary TRJ."""
        gate_failure = getattr(self, '_pending_nve_gate_error', None) is not None
        if istep % self.trajectory_interval != 0 and not gate_failure:
            return
        if not self._is_io_rank():
            self._io_barrier()
            return
        coords = np.asarray(coordinates, dtype=np.float64).reshape((-1, 3))
        if hasattr(self, 'r_all') and len(coords) == len(self.r_all):
            velocities = np.asarray(self.v_all, dtype=np.float64).reshape(coords.shape)
        else:
            velocities = np.asarray(self.vel, dtype=np.float64).reshape(coords.shape)
        ncv = self.odp.ncv if getattr(self, 'odp', None) is not None else 0
        trajectory_nstate = int(np.asarray(self.coef).size)
        dtype = _namd_trajectory_dtype(trajectory_nstate, len(coords), ncv)
        if not os.path.exists(self.trajectory_file) or os.path.getsize(self.trajectory_file) == 0:
            header = {
                'schema_version': NAMD_TRAJECTORY_SCHEMA_VERSION,
                'nstate': trajectory_nstate,
                'natom': len(coords),
                'ncv': ncv,
                'record_bytes': dtype.itemsize,
                'signature': self._restart_signature(),
                'electronic_representation': getattr(
                    self, '_trajectory_representation', 'same_spin_adiabatic'),
                'ensemble': 'NVE',
                'initial_temperature_kelvin': getattr(self, 'init_temp', None),
                'odp': self._odp_provenance(),
                'wham': {
                    'reaction_coordinate_field': 'odp_xi',
                    'window_field': 'odp_window',
                    'bias_field': 'odp_bias_hartree',
                    'unbiased_potential_field': 'e_unbiased_pot_hartree',
                    'energy_unit': 'hartree',
                    'temperature_note': (
                        'NVE trajectory: initial_temperature_kelvin is not an '
                        'NVT thermostat temperature'),
                },
                'units': {
                    'time': 'fs', 'coordinates': 'bohr',
                    'velocities': 'bohr/atomic_time', 'energies': 'hartree',
                    'tdc': 'atomic_time^-1',
                },
                'reference_source': {'0': 'none', '1': 'TD-Baeck-An',
                                     '2': 'analytic', '127': 'other'},
                'gate_metrics': [
                    'candidate_diagonal_max', 'candidate_antisymmetry_max',
                    'reference_diagonal_max', 'reference_antisymmetry_max',
                    'pair_rms_error', 'pair_max_error', 'max_tolerance_ratio',
                ],
                'gate_counts': [
                    'compared_pairs', 'invariant_failures', 'reference_failures',
                ],
                'gate_verdict': {'-1': 'none', '0': 'not_evaluable',
                                 '1': 'pass', '2': 'fail'},
                'nve_metrics': [
                    'total_energy_drift', 'step_energy_change',
                    'transition_energy_jump', 'drift_rate_hartree_per_fs',
                ],
                'nve_verdict': {'-1': 'off', '1': 'pass', '2': 'fail'},
            }
            encoded = json.dumps(header, sort_keys=True).encode('utf-8')
            with open(self.trajectory_file, 'wb') as stream:
                stream.write(NAMD_TRAJECTORY_MAGIC)
                stream.write(struct.pack('<Q', len(encoded)))
                stream.write(encoded)
                stream.flush()
                os.fsync(stream.fileno())
        else:
            header, existing = read_namd_trajectory(self.trajectory_file)
            del existing
            if (int(header['nstate']) != trajectory_nstate
                    or int(header['natom']) != len(coords)
                    or int(header.get('ncv', 0)) != ncv
                    or int(header['record_bytes']) != dtype.itemsize
                    or header.get('signature') != self._restart_signature()):
                raise ValueError('NAMD trajectory schema/model mismatch on append')

        record = np.zeros(1, dtype=dtype)
        for field in (
                'rng', 'e_unbiased_pot_hartree', 'e_pot_hartree',
                'e_kin_hartree', 'e_tot_hartree',
                'state_energies', 'populations', 'coef_real', 'coef_imag',
                'coordinates_bohr', 'velocities_au', 'state_overlap',
                'overlap_tdc_au', 'reference_tdc_au', 'gate_metrics',
                'nve_metrics',
                'odp_xi', 'odp_cv_raw', 'odp_cv_scaled',
                'odp_cv_perpendicular', 'odp_perpendicular_norm',
                'odp_bias_parallel_hartree',
                'odp_bias_perpendicular_hartree', 'odp_bias_hartree',
                'tracking_phase', 'tracking_overlap', 'tracking_margin'):
            record[field] = np.nan
        record['tracking_order'] = -1
        record['gate_center_step'] = -1
        record['gate_verdict'] = -1
        record['nve_verdict'] = -1
        record['odp_window'] = -1
        gate = self._nacme_gate_last or {}
        time_fs = self._t_fs if self.dt_adaptive else istep*self.dt_fs
        record['step'] = istep
        record['time_fs'] = time_fs
        record['active'] = self.active
        record['hopped'] = int(bool(hopped))
        record['rng'] = self._last_hop_random
        record['e_unbiased_pot_hartree'] = getattr(
            self, '_unbiased_potential_energy', epot)
        record['e_pot_hartree'] = epot
        record['e_kin_hartree'] = ekin
        record['e_tot_hartree'] = epot + ekin
        record['populations'] = np.abs(self.coef)**2
        record['coef_real'] = self.coef.real
        record['coef_imag'] = self.coef.imag
        record['coordinates_bohr'] = coords
        record['velocities_au'] = velocities
        if getattr(self, 'odp', None) is not None:
            if getattr(self, '_odp_last', None) is None:
                raise RuntimeError('ODP trajectory record has no native evaluation')
            record['odp_window'] = self.odp.window
            record['odp_xi'] = self._odp_last['xi']
            record['odp_cv_raw'] = self._odp_last['cv_raw']
            record['odp_cv_scaled'] = self._odp_last['cv_scaled']
            record['odp_cv_perpendicular'] = self._odp_last['cv_perpendicular']
            record['odp_perpendicular_norm'] = self._odp_last['perpendicular_norm']
            record['odp_bias_parallel_hartree'] = self._odp_last['energy_parallel']
            record['odp_bias_perpendicular_hartree'] = (
                self._odp_last['energy_perpendicular'])
            record['odp_bias_hartree'] = self._odp_last['energy']
        state_energies = getattr(self, '_trajectory_state_energies', None)
        if state_energies is None:
            try:
                state_energies = self.mol.data['OQP::td_energies']
            except (KeyError, TypeError):
                state_energies = None
        if state_energies is not None:
            values = np.asarray(state_energies, dtype=float).reshape(-1)
            if values.size >= trajectory_nstate:
                record['state_energies'] = values[:trajectory_nstate]
        if (self._last_state_overlap is not None
                and np.shape(self._last_state_overlap)
                == (trajectory_nstate, trajectory_nstate)):
            record['state_overlap'] = self._last_state_overlap
        if (self._last_overlap_tdc is not None
                and np.shape(self._last_overlap_tdc)
                == (trajectory_nstate, trajectory_nstate)):
            record['overlap_tdc_au'] = self._last_overlap_tdc
        if (self._nacme_reference_tdc is not None
                and np.shape(self._nacme_reference_tdc)
                == (trajectory_nstate, trajectory_nstate)):
            record['reference_tdc_au'] = self._nacme_reference_tdc
            record['reference_mask'] = self._nacme_reference_mask
            record['reference_source'] = self._nacme_reference_source
        if gate:
            record['gate_center_step'] = (-1 if gate.get('center_step') is None
                                          else int(gate['center_step']))
            record['gate_verdict'] = {
                'not_evaluable': 0, 'pass': 1, 'fail': 2,
            }.get(gate.get('verdict'), -1)
            record['gate_counts'] = (
                gate.get('compared_pairs', 0), gate.get('invariant_failures', 0),
                gate.get('reference_failures', 0),
            )
            record['gate_metrics'] = (
                gate.get('candidate_diagonal_max', np.nan),
                gate.get('candidate_antisymmetry_max', np.nan),
                gate.get('reference_diagonal_max', np.nan),
                gate.get('reference_antisymmetry_max', np.nan),
                gate.get('pair_rms_error', np.nan),
                gate.get('pair_max_error', np.nan),
                gate.get('max_tolerance_ratio', np.nan),
            )
        nve = self._nve_gate_last or {}
        if nve and nve.get('verdict') != 'off':
            record['nve_verdict'] = {'pass': 1, 'fail': 2}.get(
                nve.get('verdict'), -1)
        if nve:
            record['nve_streak'] = nve.get('consecutive_failures', 0)
            record['nve_metrics'] = (
                nve.get('drift', np.nan), nve.get('step_change', np.nan),
                nve.get('transition_energy_jump', np.nan),
                nve.get('drift_rate', np.nan),
            )
        tracking = self.mol.get_state_tracking()
        if (tracking is not None
                and all(np.asarray(tracking[name]).size == trajectory_nstate
                        for name in ('order', 'phase_step', 'matched_overlap',
                                     'margin'))):
            record['tracking_valid'] = 1
            record['tracking_order'] = np.asarray(tracking['order'], dtype=np.int64)
            record['tracking_phase'] = np.asarray(tracking['phase_step'], dtype=float)
            record['tracking_overlap'] = np.asarray(tracking['matched_overlap'], dtype=float)
            record['tracking_margin'] = np.asarray(tracking['margin'], dtype=float)

        with open(self.trajectory_file, 'ab') as stream:
            stream.write(record.tobytes(order='C'))
            stream.flush()
            os.fsync(stream.fileno())
        self._io_barrier()

    def _restart_signature(self):
        cfg = self.mol.config
        md = cfg.get('md', {})
        identity = {
            'method': cfg['input'].get('method', ''),
            'charge': cfg['input'].get('charge', ''),
            'functional': cfg['input'].get('functional', ''),
            'basis': cfg['input'].get('basis', ''),
            'scf_type': cfg.get('scf', {}).get('type', ''),
            'scf_multiplicity': cfg.get('scf', {}).get('multiplicity', ''),
            'tdhf_type': cfg['tdhf'].get('type', ''),
            'tdhf_multiplicity': cfg['tdhf'].get('multiplicity', ''),
            'nstate': cfg['tdhf'].get('nstate', ''),
            'tlf': cfg['tdhf'].get('tlf', ''),
            'dt_fs': self.dt_fs, 'seed': self.seed,
            'rng_stream': self.rng_stream,
            'substep': md.get('substep', ''),
            'decoherence': md.get('decoherence', ''),
            'edc_c': md.get('edc_c', ''), 'thrshe': md.get('thrshe', ''),
            'tdc': md.get('tdc', ''), 'trivial': md.get('trivial', ''),
            'trivial_thresh': md.get('trivial_thresh', ''),
            'first_hop_step': md.get('first_hop_step', ''),
            'soc': md.get('soc', ''), 'soc_basis': md.get('soc_basis', ''),
            'trajectory_representation': getattr(
                self, '_trajectory_representation', 'same_spin_adiabatic'),
            'odp': self._odp_provenance(),
            'system': getattr(
                self, '_restart_system_identity', {'kind': 'unavailable'}),
        }
        return json.dumps(identity, sort_keys=True, separators=(',', ':'))

    @staticmethod
    def _checkpoint_optional(payload, name, value):
        if value is None:
            payload[f'has_{name}'] = np.array([0], dtype=np.int8)
            payload[name] = np.empty(0, dtype=np.float64)
        else:
            payload[f'has_{name}'] = np.array([1], dtype=np.int8)
            payload[name] = np.asarray(value)

    def _save_restart(self, istep, coordinates, velocities, acceleration):
        """Atomically save all state needed for phase-continuous continuation."""
        if istep % self.restart_interval != 0:
            return
        if not self._is_io_rank():
            self._io_barrier()
            return
        if self.prev_data is None or self.prev_xyz is None:
            self._io_barrier()
            return
        nuclear_state = tuple(np.asarray(value, dtype=float) for value in (
            coordinates, velocities, acceleration))
        if (nuclear_state[0].shape != nuclear_state[1].shape
                or nuclear_state[0].shape != nuclear_state[2].shape
                or nuclear_state[0].size == 0
                or not all(np.all(np.isfinite(value)) for value in nuclear_state)):
            raise RuntimeError(
                f'refusing to overwrite the last-good NAMD restart with an '
                f'invalid nuclear state at step {istep}')
        prev_keys = sorted(self.prev_data)
        payload = {
            'schema_version': np.array([NAMD_RESTART_SCHEMA_VERSION], dtype=np.int64),
            'signature': np.array([self._restart_signature()]),
            'step': np.array([istep], dtype=np.int64),
            'time_fs': np.array([self._t_fs if self.dt_adaptive else istep*self.dt_fs]),
            'active': np.array([self.active], dtype=np.int64),
            'rng_step': np.array([self._rng_step], dtype=np.int64),
            'gate_failures': np.array([self._nacme_gate_failures], dtype=np.int64),
            'nve_failures': np.array([self._nve_gate_failures], dtype=np.int64),
            'coordinates': np.asarray(coordinates, dtype=np.float64),
            'velocities': np.asarray(velocities, dtype=np.float64),
            'acceleration': np.asarray(acceleration, dtype=np.float64),
            'coef_real': np.asarray(self.coef.real, dtype=np.float64),
            'coef_imag': np.asarray(self.coef.imag, dtype=np.float64),
            'prev_xyz': np.asarray(self.prev_xyz, dtype=np.float64),
            'prev_keys': np.asarray(prev_keys, dtype=np.str_),
            'nacme_audit_bytes': np.array([
                os.path.getsize(self.nacme_audit_file)
                if os.path.isfile(self.nacme_audit_file) else 0
            ], dtype=np.int64),
            'odp_provenance': np.array([
                json.dumps(self._odp_provenance(), sort_keys=True,
                           separators=(',', ':'))
            ]),
        }
        for index, key in enumerate(prev_keys):
            value = np.asarray(self.prev_data[key])
            if value.dtype == object:
                raise TypeError(f'NAMD restart cannot serialize ragged tag {key!r}')
            payload[f'prev_{index}'] = value
        for name, value in (
                ('ba_energy_left', self._ba_energy_left),
                ('ba_energy_center', self._ba_energy_center),
                ('ba_tdc_left', self._ba_tdc_left),
                ('ba_dt_left', self._ba_dt_left),
                ('nve_reference_energy', self._nve_reference_energy),
                ('nve_previous_energy', self._nve_previous_energy)):
            self._checkpoint_optional(payload, name, value)

        directory = os.path.dirname(self.restart_file) or '.'
        descriptor, temporary = tempfile.mkstemp(
            prefix='.namd-restart-', suffix='.tmp', dir=directory)
        try:
            with os.fdopen(descriptor, 'wb') as stream:
                np.savez_compressed(stream, **payload)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, self.restart_file)
            self._write_restart_manifest()
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        self._io_barrier()

    def _write_restart_manifest(self):
        """Write a directly runnable ``restart.oqp`` beside the checkpoint."""
        canonical = str(getattr(self.mol, 'oqp_canonical_input', '') or '').strip()
        if not canonical:
            source = getattr(self.mol, 'oqp_input_source', None)
            if source and str(source).lower().endswith('.oqp') and os.path.isfile(source):
                with open(source, 'r', encoding='utf-8') as stream:
                    canonical = stream.read().strip()
        if not canonical:
            dump_log(
                self.mol,
                title=('NAMD checkpoint saved, but restart.oqp was not generated: '
                       'the run did not originate from canonical .oqp input'),
            )
            return
        from oqp.utils.oqp_input import (
            CallSpec, CalculationSpec, parse_canonical_oqp, render_canonical_oqp,
        )
        spec = parse_canonical_oqp(canonical)
        if spec.driver.name != 'namd':
            raise ValueError('cannot create restart.oqp from a non-NAMD request')
        directory = os.path.dirname(self.restart_manifest_file) or '.'
        kwargs = dict(spec.driver.kwargs)
        kwargs.update({
            'restart': True,
            'restart_file': os.path.relpath(self.restart_file, directory),
            'trajectory_file': os.path.relpath(self.trajectory_file, directory),
            'nacme_audit_file': os.path.relpath(self.nacme_audit_file, directory),
        })
        driver = CallSpec(
            spec.driver.name, spec.driver.args, kwargs, spec.driver.explicit)
        restart_spec = CalculationSpec(
            spec.model, spec.functional, spec.basis, spec.model_options,
            spec.options, driver, spec.modifiers, spec.source_text)
        rendered = render_canonical_oqp(restart_spec)
        descriptor, temporary = tempfile.mkstemp(
            prefix='.restart-oqp-', suffix='.tmp', dir=directory)
        try:
            with os.fdopen(descriptor, 'w', encoding='utf-8') as stream:
                stream.write(rendered)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, self.restart_manifest_file)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)

    def _load_restart(self):
        """Load and validate a same-spin NAMD restart without pickle data."""
        if not self.restart_requested:
            return None
        if not os.path.isfile(self.restart_file):
            raise FileNotFoundError(f'NAMD restart file not found: {self.restart_file}')
        with np.load(self.restart_file, allow_pickle=False) as saved:
            version = int(saved['schema_version'][0])
            if version != NAMD_RESTART_SCHEMA_VERSION:
                raise ValueError(f'unsupported NAMD restart schema {version}')
            signature = str(saved['signature'][0])
            if signature != self._restart_signature():
                raise ValueError(
                    'NAMD restart configuration/system identity does not '
                    'match the current run')
            saved_odp = str(saved['odp_provenance'][0])
            current_odp = json.dumps(
                self._odp_provenance(), sort_keys=True, separators=(',', ':'))
            if saved_odp != current_odp:
                raise ValueError('NAMD restart ODP definition/metric mismatch')
            keys = [str(key) for key in saved['prev_keys']]
            self.prev_data = {
                key: np.array(saved[f'prev_{index}'], copy=True)
                for index, key in enumerate(keys)
            }
            self.prev_xyz = np.array(saved['prev_xyz'], copy=True)
            self.mol.put_data(self.prev_data)
            self.active = int(saved['active'][0])
            self.coef = (np.array(saved['coef_real'], copy=True)
                         + 1j*np.array(saved['coef_imag'], copy=True))
            self._rng_step = int(saved['rng_step'][0])
            self._last_hop_random = np.nan
            self._nacme_gate_failures = int(saved['gate_failures'][0])
            self._nve_gate_failures = int(saved['nve_failures'][0])
            self._t_fs = float(saved['time_fs'][0])
            for name in ('ba_energy_left', 'ba_energy_center', 'ba_tdc_left',
                         'ba_dt_left', 'nve_reference_energy',
                         'nve_previous_energy'):
                value = (np.array(saved[name], copy=True)
                         if int(saved[f'has_{name}'][0]) else None)
                if name in ('ba_dt_left', 'nve_reference_energy',
                            'nve_previous_energy') and value is not None:
                    value = float(value.reshape(-1)[0])
                setattr(self, f'_{name}', value)
            result = {
                'step': int(saved['step'][0]),
                'coordinates': np.array(saved['coordinates'], copy=True),
                'velocities': np.array(saved['velocities'], copy=True),
                'acceleration': np.array(saved['acceleration'], copy=True),
                'nacme_audit_bytes': int(saved['nacme_audit_bytes'][0]),
            }
        self._reconcile_trajectory_with_restart(result['step'])
        self._reconcile_nacme_audit_with_restart(
            result['step'], result['nacme_audit_bytes'])
        dump_log(
            self.mol,
            title=(f'NAMD restart loaded: step={result["step"]} '
                   f'file={self.restart_file} phase_history=restored '
                   f'rng=({self.seed},{self.rng_stream},step)'),
        )
        return result

    def _reconcile_trajectory_with_restart(self, checkpoint_step):
        """Discard only uncommitted TRJ tail records newer than checkpoint."""
        if not self._is_io_rank():
            self._io_barrier()
            return
        if not os.path.isfile(self.trajectory_file) or os.path.getsize(self.trajectory_file) == 0:
            self._io_barrier()
            return
        with open(self.trajectory_file, 'rb') as stream:
            if stream.read(8) != NAMD_TRAJECTORY_MAGIC:
                raise ValueError('restart trajectory is not an OpenQP dense TRJ')
            header_size = struct.unpack('<Q', stream.read(8))[0]
        offset = 16 + header_size
        header, records = read_namd_trajectory(self.trajectory_file)
        if header.get('signature') != self._restart_signature():
            raise ValueError('restart trajectory and checkpoint model mismatch')
        steps = np.array(records['step'], copy=True)
        keep = int(np.count_nonzero(steps <= int(checkpoint_step)))
        if keep < len(records):
            del records
            dtype = _namd_trajectory_dtype(
                int(header['nstate']), int(header['natom']),
                int(header.get('ncv', 0)))
            with open(self.trajectory_file, 'r+b') as stream:
                stream.truncate(offset + keep*dtype.itemsize)
                stream.flush()
                os.fsync(stream.fileno())
            dump_log(
                self.mol,
                title=(f'NAMD restart removed {len(steps) - keep} uncommitted '
                       f'trajectory record(s) after step {checkpoint_step}'),
            )
        else:
            del records
        self._io_barrier()

    def _reconcile_nacme_audit_with_restart(self, checkpoint_step,
                                            checkpoint_bytes):
        """Restore the exact machine-readable audit prefix at a checkpoint."""
        if not self._is_io_rank():
            self._io_barrier()
            return
        checkpoint_bytes = int(checkpoint_bytes)
        if checkpoint_bytes < 0:
            raise ValueError('restart checkpoint has an invalid NACME audit size')
        current_bytes = (
            os.path.getsize(self.nacme_audit_file)
            if os.path.isfile(self.nacme_audit_file) else 0
        )
        if current_bytes < checkpoint_bytes:
            raise ValueError(
                'NAMD restart NACME audit is missing committed checkpoint data')
        if current_bytes > checkpoint_bytes:
            with open(self.nacme_audit_file, 'r+b') as stream:
                stream.truncate(checkpoint_bytes)
                stream.flush()
                os.fsync(stream.fileno())
            dump_log(
                self.mol,
                title=(f'NAMD restart removed {current_bytes - checkpoint_bytes} '
                       f'uncommitted NACME audit byte(s) after checkpoint step '
                       f'{checkpoint_step}'),
            )
        self._io_barrier()

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
        params[_P_RAND] = self._hop_random()
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
        s = canonical_state_overlap(
            np.asarray(mol.data["OQP::td_states_overlap"]).reshape((n, n))
        )
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
        self._prepare_md_outputs()
        restart = self._load_restart()
        if restart is None:
            r = mol.get_system().reshape((self.natom, 3))   # bohr
            # initial electronic structure + force on the active state
            self._electronic(with_overlap=False)
            accel = -self._active_gradient() / self.mass[:, None]
            self._record_previous(r)
            self._log_step(0, r)
            self._save_restart(0, r, self.vel, accel)
            start_step = 0
        else:
            r = restart['coordinates'].reshape((self.natom, 3))
            self.vel = restart['velocities'].reshape((self.natom, 3))
            accel = restart['acceleration'].reshape((self.natom, 3))
            mol.update_system(r.reshape(-1))
            start_step = restart['step']

        for istep in range(start_step + 1, self.nstep + 1):
            # velocity-Verlet position update
            r = r + self.vel * self.dt + 0.5 * accel * self.dt ** 2
            mol.update_system(r.reshape(-1))

            # electronic structure at the new geometry (with overlap vs previous)
            self._electronic(with_overlap=True)
            accel_new = -self._active_gradient() / self.mass[:, None]

            # velocity-Verlet velocity update
            self.vel = self.vel + 0.5 * (accel + accel_new) * self.dt

            # state overlap (couplings) and FSSH hop
            self._state_overlap(istep)
            active_old = self.active
            odp = self._evaluate_odp(r)
            bias_energy = 0.0 if odp is None else odp['energy']
            energy_before_transition = (
                0.5*np.sum(self.mass[:, None]*self.vel**2)
                + float(np.asarray(mol.energies)[active_old])
                + bias_energy
            )
            if self._prepare_hop_step(istep):
                new_active, hopped = self._hop()
            else:
                new_active, hopped = self.active, False

            active_changed = new_active != active_old
            if active_changed:
                self.active = new_active
                # force for the next step is on the new active surface. This
                # also covers trivial-crossing following, where the Fortran
                # kernel can update ACTIVE without marking HOPPED.
                accel_new = -self._active_gradient() / self.mass[:, None]
                energy_after_transition = (
                    0.5*np.sum(self.mass[:, None]*self.vel**2)
                    + float(np.asarray(mol.energies)[self.active])
                    + bias_energy
                )
                transition_energy_jump = (
                    energy_after_transition - energy_before_transition)
            else:
                transition_energy_jump = np.nan

            accel = accel_new
            self._record_previous(r)
            self._log_step(
                istep, r, hopped=hopped,
                transition_energy_jump=transition_energy_jump)
            self._save_restart(istep, r, self.vel, accel)

        dump_log(mol, title='PyOQP: NAMD trajectory complete')

    # ------------------------------------------------------------------ #
    # helpers
    # ------------------------------------------------------------------ #
    def _record_previous(self, r):
        self.prev_xyz = copy.deepcopy(r.reshape(-1))
        self.prev_data = copy.deepcopy(self.mol.get_data())

    def _log_step(self, istep, r, hopped=False, transition_energy_jump=np.nan):
        mol = self.mol
        e = np.array(mol.energies)
        ekin = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
        self._unbiased_potential_energy = float(e[self.active])
        odp = self._evaluate_odp(r)
        epot = self._unbiased_potential_energy + (
            0.0 if odp is None else odp['energy'])
        pops = np.abs(self.coef) ** 2
        self._update_nve_gate(istep, epot, ekin, transition_energy_jump)
        dump_log(
            mol,
            title=(f'NAMD step {istep:6d}  t={(self._t_fs if self.dt_adaptive else istep*self.dt_fs):9.3f} fs  '
                   f'active={self.active}  E_tot={ekin+epot:.8f}  '
                   f'E_pot={epot:.8f}  E_kin={ekin:.8f}  '
                   f'U_ODP={(0.0 if odp is None else odp["energy"]):.8f}  '
                   f'hop={hopped}  {self._hop_rng_log()}  '
                   f'pop={np.array2string(pops, precision=4)}'),
        )
        self._write_md_trajectory(istep, r, epot, ekin, hopped)
        self._enforce_nve_gate()


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
      * QM-only FSSH hop with rescaling of QM velocities only.
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

        # Resolve QM/MM auxiliary files (PDB, local force-field XMLs) relative to
        # the input file's directory when a bare/relative name does not resolve
        # against the current working directory. This lets the examples run from
        # any CWD (e.g. `openqp --run_tests`, which executes each example by its
        # full path). OpenMM built-in force fields (amber14-all.xml, ...) are left
        # untouched: the join is only used when it actually points at a file.
        inp_dir = os.path.dirname(os.path.abspath(mol.input_file)) \
            if getattr(mol, 'input_file', None) else ''

        def _resolve_aux(name):
            if name and inp_dir and not os.path.isabs(name) \
                    and not os.path.exists(name):
                cand = os.path.join(inp_dir, name)
                if os.path.exists(cand):
                    return cand
            return name

        q = mol.config['qmmm']
        pdb_file = _resolve_aux(q['pdb_file'])
        ff_files = [_resolve_aux(s) for s in
                    str(q['forcefield_files']).replace(',', ' ').split() if s]
        self.qm_atoms = np.array(_parse_int_list(q['qm_atoms']), dtype=int)
        self.cutoff = _resolve_cutoff(str(q['cutoff']).strip())   # NoCutoff | PME | Ewald | ...
        self.periodic = self.cutoff is not app.NoCutoff
        embedding = str(q['embedding']).strip()
        frontier_scheme = str(q.get('frontier_scheme', 'none')).strip()

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
            frontier_scheme=frontier_scheme,
        )
        self.mm = self.driver.mm_systems

        # full-system state (atomic units)
        self.natom_all = self.pdb.topology.getNumAtoms()
        if self.odp is not None:
            self.odp.validate_atom_count(self.natom_all)
        pos_nm = np.array(self.pdb.positions.value_in_unit(u.nanometer))
        self.r_all = pos_nm * NM_TO_BOHR                       # (natom_all, 3) bohr
        sys0 = self.mm["sys0"]
        self.m_all = np.array([
            sys0.getParticleMass(i).value_in_unit(u.dalton) for i in range(self.natom_all)
        ]) * AMU_TO_AU                                          # electron masses
        self._restart_system_identity = self._qmmm_restart_system_identity(
            sys0, q)

        # full-system Maxwell-Boltzmann velocities (a.u.), COM removed
        sig = np.sqrt(KB_HARTREE * self.init_temp / self.m_all)
        self.v_all = self._counter_normals((self.natom_all, 3)) * sig[:, None]
        p = (self.m_all[:, None] * self.v_all).sum(axis=0)
        self.v_all -= p / self.m_all.sum()

        # sync the QM Molecule geometry from the pdb QM atoms
        self._sync_positions()
        # QM-region masses for the hop (already set by super from mol.get_mass())
        self.qm_mass = self.mass.copy()
        # rigid-water (SHAKE/RATTLE) constraints for the MM region
        self._build_constraints()

    # ------------------------------------------------------------------ #
    def _qmmm_restart_system_identity(self, system, qmmm_config):
        """Bind QM/MM restarts to atoms, topology, selection, and force field."""
        atoms = list(self.pdb.topology.atoms())
        atomic_numbers = [
            0 if atom.element is None else atom.element.atomic_number
            for atom in atoms
        ]
        atom_metadata = [
            (atom.name, atom.residue.name, atom.residue.index,
             atom.residue.chain.id)
            for atom in atoms
        ]
        bonds = np.asarray(sorted(
            (min(atom1.index, atom2.index), max(atom1.index, atom2.index))
            for atom1, atom2 in self.pdb.topology.bonds()
        ), dtype='<i8').reshape((-1, 2))
        system_xml = self._mm.XmlSerializer.serialize(system)
        digest = _restart_identity_digest(
            array_parts=(
                ('atomic_numbers', atomic_numbers, '<i8'),
                ('masses_electron', self.m_all, '<f8'),
                ('initial_coordinates_bohr', self.r_all, '<f8'),
                ('qm_atoms', self.qm_atoms, '<i8'),
                ('bonds', bonds, '<i8'),
            ),
            text_parts=(
                ('atom_metadata', json.dumps(
                    atom_metadata, separators=(',', ':'))),
                ('qmmm_config', json.dumps(
                    qmmm_config, sort_keys=True, separators=(',', ':'),
                    default=str)),
                ('openmm_system', system_xml),
            ),
        )
        return {
            'kind': 'qmmm', 'natom': int(self.natom_all),
            'nqm': int(len(self.qm_atoms)), 'sha256': digest,
        }

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

        if is_tb_method(str(mol.config['input']['method'])):
            # DFTB electrostatic embedding: the openqp-dftb library folds the
            # per-atom MM potential (Hartree/e) directly into the SCC
            # Hamiltonian, so there is no ESPF operator / hcore mutation here.
            # POTQM stays zeroed (same resolution as the native path below).
            # Only the full-ESPF scheme is supported: the legacy 'split'
            # scheme would double-count the coupling already inside the
            # embedded DFTB state energy.
            if not getattr(self.driver, 'espf_full', False):
                raise NotImplementedError(
                    "NAMD QM/MM with method=dftb/xtb requires [qmmm] embedding="
                    "electrostatic/espf (full-ESPF scheme).")
            mol.dftb_external_potential = np.asarray(potmm, dtype=float)
            sp = SinglePoint(mol)
            ref = sp.reference()
            if with_overlap:
                mol.back_door = (self.prev_xyz, self.prev_data)
                BasisOverlap(mol).overlap()
            sp.excitation(ref)
            LastStep(mol).compute(mol)
            return potmm, potqm

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
        if is_tb_method(str(mol.config['input']['method'])):
            # The DFTB analytic gradient is d(E_embedded)/dR at FIXED potential
            # values: the charge-response (Pulay-type) coupling term is already
            # inside it, so the native grad_esp_qmmm(_excited)/OQP::ESPF_GRAD
            # addition is skipped. The adapter has also just published
            # OQP::partial_charges (relaxed net charges of the active state)
            # for the classical coupling forces in _total_force.
            return g
        # ESPF_ROHF=1: use the ROHF reference density for ESPF charges and
        # gradient regardless of the target excited state.  Combined with
        # ESPF_HARD_GRID=1 this provides stable QM/MM energy conservation for
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

        if getattr(self.driver, "espf_full", False):
            force, epot = self._total_force_espf(potmm, gqm, pchg)
            return self._apply_odp_to_force_energy(force, epot)

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
        return self._apply_odp_to_force_energy(f_all, epot)

    def _apply_odp_to_force_energy(self, force, epot):
        """Add the conservative native ODP term to a full QM/MM force."""
        self._unbiased_potential_energy = float(epot)
        odp = self._evaluate_odp(self.r_all)
        if odp is None:
            return force, float(epot)
        return force + odp['force'], float(epot) + odp['energy']

    def _total_force_espf(self, potmm, gqm, pchg):
        """Full-ESPF force/energy assembly (mirrors OpenQpQMMM._assemble_force_espf,
        in atomic units). QM-MM electrostatics go entirely through ESPF: OpenMM
        is pure MM-MM, the embedded SCF carries the electronic coupling, the
        nuclear-MM term sum_A Z_A phi_A is added explicitly, and the analytic
        coupling force is applied to QM and MM atoms."""
        mol = self.mol
        u = self._u
        nqm = len(self.qm_atoms)
        if len(self.driver.link_atoms) and gqm.shape[0] == nqm:
            raise NotImplementedError(
                "NAMD full-ESPF across a covalent boundary needs the QM link "
                "atoms in the NAMD mol (only the driver has them). Use "
                "QMMM_MD for covalent-boundary QM/MM, or add link atoms to the "
                "NAMD QM geometry.")

        emm_q, gmm_q = self.driver.forces_mm(pchg)     # pure MM-MM
        gmm = np.array(gmm_q.value_in_unit(u.kilojoule_per_mole / u.nanometer)) / HABOHR_TO_KJMOLNM
        emm = emm_q.value_in_unit(u.kilojoule_per_mole) * KJMOL_TO_HARTREE

        fq, fm, mm_idx = self.driver._coupling_forces(pchg)   # a.u. (Ha/bohr)

        f_all = gmm.copy()
        for a, link in enumerate(self.driver.link_atoms):
            gl, fl = gqm[nqm + a], fq[nqm + a]
            gqm[link.host_row] = gqm[link.host_row] + (1.0 - link.g) * gl
            fq[link.host_row] = fq[link.host_row] + (1.0 - link.g) * fl
            f_all[link.mm_index] = f_all[link.mm_index] - link.g * gl + link.g * fl
        for k, i in enumerate(self.qm_atoms):
            f_all[i] = f_all[i] - gqm[k] + fq[k]
        for j, m in enumerate(mm_idx):
            f_all[m] = f_all[m] + fm[j]
        f_all -= f_all.mean(axis=0)

        if is_tb_method(str(mol.config['input']['method'])):
            # The DFTB embedded state energy is already COMPLETE: the library
            # folds the MM potential into the SCC Hamiltonian and the returned
            # energy includes the full net-charge coupling
            #   E_ext = sum_A q_A phi_A   (q_A = NET atomic charge, cores +
            # electrons; dE/dphi_A = +q_A, verified by finite differences).
            # The native path below instead carries only the electronic
            # coupling in the embedded SCF and must add the nuclear term
            # sum_A Z_A phi_A -- do NOT add it for dftb.
            eqm = float(mol.energies[self.active])
        else:
            eqm = float(mol.energies[self.active]) + float(
                np.dot(np.array(mol.get_atoms2("charge")), potmm))
        return f_all, eqm + emm

    # ------------------------------------------------------------------ #
    def run(self):
        mol = self.mol
        dump_log(mol, title='PyOQP: QM/MM Tully FSSH Nonadiabatic Molecular Dynamics')
        self._prepare_md_outputs()
        restart = self._load_restart()
        if restart is None:
            # initial electronic structure + force
            self._sync_positions()
            potmm0, _ = self._electronic_qmmm(with_overlap=False)
            f_all, epot = self._total_force(potmm0)
            accel = f_all / self.m_all[:, None]
            self._rattle(self.r_all, self.v_all)      # constrained velocities
            self._thermalize_initial()
            self.prev_xyz = copy.deepcopy(self.r_all[self.qm_atoms].reshape(-1))
            self.prev_data = copy.deepcopy(mol.get_data())
            self._log_qmmm(0, epot)
            self._save_restart(0, self.r_all, self.v_all, accel)
            start_step = 0
        else:
            self.r_all = restart['coordinates'].reshape((self.natom_all, 3))
            self.v_all = restart['velocities'].reshape((self.natom_all, 3))
            accel = restart['acceleration'].reshape((self.natom_all, 3))
            self.vel = self.v_all[self.qm_atoms].copy()
            self._sync_positions()
            start_step = restart['step']

        for istep in range(start_step + 1, self.nstep + 1):
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
            self._state_overlap(istep)
            self.vel = self.v_all[self.qm_atoms].copy()       # hop sees QM velocities
            active_old = self.active
            energy_before_transition = (
                0.5*np.sum(self.m_all[:, None]*self.v_all**2) + epot)
            if self._prepare_hop_step(istep):
                new_active, hopped = self._hop()
            else:
                new_active, hopped = self.active, False
            self.v_all[self.qm_atoms] = self.vel              # write back rescaled QM velocities
            active_changed = new_active != active_old
            if active_changed:
                self.active = new_active
                f_all, epot = self._total_force(potmm)
                accel_new = f_all / self.m_all[:, None]
                energy_after_transition = (
                    0.5*np.sum(self.m_all[:, None]*self.v_all**2) + epot)
                transition_energy_jump = (
                    energy_after_transition - energy_before_transition)
            else:
                transition_energy_jump = np.nan

            accel = accel_new
            self.prev_xyz = copy.deepcopy(self.r_all[self.qm_atoms].reshape(-1))
            self.prev_data = copy.deepcopy(mol.get_data())
            self._log_qmmm(
                istep, epot, hopped=hopped,
                transition_energy_jump=transition_energy_jump)
            self._save_restart(istep, self.r_all, self.v_all, accel)

        dump_log(mol, title='PyOQP: QM/MM NAMD trajectory complete')

    def _log_qmmm(self, istep, epot, hopped=False,
                  transition_energy_jump=np.nan):
        ekin = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
        pops = np.abs(self.coef) ** 2
        self._update_nve_gate(istep, epot, ekin, transition_energy_jump)
        dump_log(
            self.mol,
            title=(f'QMMM-NAMD step {istep:6d}  t={(self._t_fs if self.dt_adaptive else istep*self.dt_fs):9.3f} fs  '
                   f'active={self.active}  E_tot={ekin+epot:.8f}  '
                   f'E_pot={epot:.8f}  E_kin={ekin:.8f}  '
                   f'U_ODP={(0.0 if self._odp_last is None else self._odp_last["energy"]):.8f}  '
                   f'hop={hopped}  {self._hop_rng_log()}  '
                   f'pop={np.array2string(pops, precision=4)}'),
        )
        self._write_md_trajectory(
            istep, self.r_all, epot, ekin, hopped)
        self._enforce_nve_gate()


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
        self._trajectory_representation = 'soc_adiabatic'
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
        # (e.g. 'S1' or 'T0') instead of a fixed spin-adiabatic index -- robust when the
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
        """If [md] init_state names an MCH state (S0/S1/.../T0/T1/...), set the
        active spin-adiabatic state to the adiabat with the largest character of
        that MCH state at t=0 (summing the three Ms sublevels for a triplet).
        Otherwise keep the configured integer active index."""
        label = self.init_state
        if not label:
            return
        mult, n, public_label = _parse_soc_init_state(
            label,
            self.ns,
            self.nt,
            public_labels=bool(getattr(self.mol, 'oqp_public_state_labels', False)),
        )
        if mult == 1:
            mch_idx = [n]                                  # singlet root: S0=0, S1=1, ...
        else:
            base = self.ns + n * 3
            mch_idx = [base, base + 1, base + 2]           # triplet Ms sublevels
        char = (np.abs(u[mch_idx, :]) ** 2).sum(axis=0)    # character per adiabat
        a = int(np.argmax(char))
        self.active = a + 1
        self.coef = np.zeros(self.nstate_soc, dtype=complex)
        self.coef[a] = 1.0 + 0.0j
        requested = str(label).strip().upper()
        state_note = (public_label if requested == public_label else
                      f'{public_label} (legacy {requested})')
        dump_log(self.mol, title=(
            f'SOC-NAMD: initial active set to adiabat {self.active} by {state_note} '
            f'character ({char[a]*100:.1f}% {public_label})'))

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

        is_dftb = is_tb_method(mol.config['input']['method'])

        _select_response_manifold(mol, 1)
        sing = sp.excitation(ref)
        self.sing_energies = np.array(sing, dtype=float)
        self.sbvec = np.array(mol.data['OQP::td_bvec_mo']).copy()
        mol.data['OQP::td_singlet_energies'] = mol.data['OQP::td_energies'].copy()
        mol.data['OQP::td_bvec_mo_s'] = mol.data['OQP::td_bvec_mo'].copy()

        _select_response_manifold(mol, 3)
        trip = sp.excitation(ref)
        self.trip_energies = np.array(trip, dtype=float)
        self.tbvec = np.array(mol.data['OQP::td_bvec_mo']).copy()
        mol.data['OQP::td_triplet_energies'] = mol.data['OQP::td_energies'].copy()
        mol.data['OQP::td_bvec_mo_t'] = mol.data['OQP::td_bvec_mo'].copy()

        if is_dftb:
            _dftb_soc_tags(mol)
        else:
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

        _select_response_manifold(mol, 1)
        mol.data['OQP::td_bvec_mo'] = self.sbvec.copy()
        mol.data['OQP::td_bvec_mo_old'] = self.prev_sbvec.copy()
        if is_tb_method(mol.config['input']['method']):
            _dftb_spatial_overlap(mol, 1)
        else:
            oqp.get_states_overlap(mol)
        s_s = canonical_state_overlap(
            np.asarray(mol.data['OQP::td_states_overlap']).reshape((ns, ns))
        )

        _select_response_manifold(mol, 3)
        mol.data['OQP::td_bvec_mo'] = self.tbvec.copy()
        mol.data['OQP::td_bvec_mo_old'] = self.prev_tbvec.copy()
        if is_tb_method(mol.config['input']['method']):
            _dftb_spatial_overlap(mol, 3)
        else:
            oqp.get_states_overlap(mol)
        s_t = canonical_state_overlap(
            np.asarray(mol.data['OQP::td_states_overlap']).reshape((nt, nt))
        )

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
        rand = self._hop_random()
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
        The triplet block is also 1-based internally (T0=target 1, ...);
        the three Ms sublevels of a spatial triplet share one target."""
        if k < self.ns:
            return 1, k + 1                               # singlet root: S0=1, S1=2, ...
        return 3, (k - self.ns) // 3 + 1                  # triplet root: T0=1, T1=2, ...

    @staticmethod
    def _mch_label(mult, target):
        """Human-readable zero-based MCH state name for an internal target."""
        return f'S{target - 1}' if mult == 1 else f'T{target - 1}'

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
            _select_response_manifold(mol, mult)
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
        self._prepare_md_outputs()

        r = mol.get_system().reshape((self.natom, 3))

        # initial electronic structure + active-surface force
        eval_ha, u = self._electronic_soc(with_overlap=False)
        self._resolve_initial_active(u)
        grad, e_pure, mult, state, w = self._soc_gradient(u, self.active, eval_ha)
        accel = -grad / self.mass[:, None]
        self._ulog = u
        self._store_prev(r, u, eval_ha)
        self._log_soc(
            0, e_pure, mult, state, w, False,
            self.e_ref + self.e0 + eval_ha)
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
            if self._prepare_hop_step(istep):
                hopped = self._propagate_and_hop(self.prev_eval, eval_ha, t)
            else:
                hopped = False
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
            self._log_soc(
                istep, e_pure, mult, state, w, hopped,
                self.e_ref + self.e0 + eval_ha)

        dump_log(mol, title='PyOQP: SOC-NAMD trajectory complete')

    def _store_prev(self, r, u, eval_ha):
        self.prev_xyz = copy.deepcopy(r.reshape(-1))
        self.prev_data = copy.deepcopy(self.mol.get_data())
        self.prev_u = u.copy()
        self.prev_eval = eval_ha.copy()
        self.prev_sbvec = self.sbvec.copy()
        self.prev_tbvec = self.tbvec.copy()

    def _log_soc(self, istep, e_pure, mult, state, w, hopped,
                 state_energies):
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
                   f'{self._hop_rng_log()}  '
                   f'dom=({self._mch_label(mult, state)},w={w:.3f})  '
                   f'pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )
        self._unbiased_potential_energy = float(e_pure)
        self._trajectory_state_energies = np.asarray(
            state_energies, dtype=float).reshape(-1)
        self._write_md_trajectory(
            istep, self.mol.get_system(), e_pure, ekin, hopped)


class NAMD_SOC_MCH(NAMD_SOC):
    """SOC-NAMD in the MCH (spin-pure) basis.

    The active state is a single MCH basis function (singlet root or one
    triplet Ms component), so the nuclear force is the exact MCH root gradient.
    Electronic amplitudes are propagated by the SOC Hamiltonian in that MCH
    basis instead of the spin-adiabatic local-diabatization propagator.
    """

    def __init__(self, mol):
        super().__init__(mol)
        self._trajectory_representation = 'soc_mch'
        self.coef = np.zeros(self.nstate_soc, dtype=complex)
        self.coef[self.active - 1] = 1.0 + 0.0j

    def _resolve_initial_mch_active(self):
        label = self.init_state
        if not label:
            return
        mult, target, public_label = _parse_soc_init_state(
            label,
            self.ns,
            self.nt,
            public_labels=bool(getattr(self.mol, 'oqp_public_state_labels', False)),
        )
        if mult == 1:
            active = target + 1                         # S0 -> MCH basis 1
        else:
            active = self.ns + target * 3 + 1            # T0 -> first triplet Ms member
        self.active = active
        self.coef[:] = 0.0
        self.coef[self.active - 1] = 1.0 + 0.0j
        requested = str(label).strip().upper()
        state_note = (public_label if requested == public_label else
                      f'{public_label} (legacy {requested})')
        dump_log(self.mol, title=(
            f'SOC-MCH-NAMD: initial active set to MCH state '
            f'{self._mch_active_label(self.active)} from {state_note}'))

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
        _select_response_manifold(mol, mult)
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

        rand = self._hop_random()
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
        self._prepare_md_outputs()

        r = mol.get_system().reshape((self.natom, 3))
        eval_ha, u = self._electronic_soc(with_overlap=False)
        self._resolve_initial_mch_active()
        h_mch = self._mch_hamiltonian_from_u(u, eval_ha)
        e_mch = self._mch_energies_abs()
        grad, e_pure, mult, state = self._mch_exact_gradient(self.active)
        accel = -grad / self.mass[:, None]
        self._store_prev(r, u, eval_ha)
        self._log_mch(0, e_pure, mult, state, False, e_mch)
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

            if self._prepare_hop_step(istep):
                hopped = self._mch_propagate_and_hop(h_mch, e_mch)
            else:
                hopped = False
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
            self._log_mch(istep, e_pure, mult, state, hopped, e_mch)

        dump_log(mol, title='PyOQP: SOC-MCH-NAMD trajectory complete')

    def _log_mch(self, istep, e_pure, mult, state, hopped, state_energies):
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
                   f'{self._hop_rng_log()}  '
                   f'grad={self._mch_label(mult, state)}  pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )
        self._unbiased_potential_energy = float(e_pure)
        self._trajectory_state_energies = np.asarray(
            state_energies, dtype=float).reshape(-1)
        self._write_md_trajectory(
            istep, self.mol.get_system(), e_pure, ekin, hopped)


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
        rescaling QM velocities only.

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
        self._trajectory_representation = 'soc_adiabatic'
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

        is_dftb = is_tb_method(mol.config['input']['method'])

        _select_response_manifold(mol, 1)
        sing = sp.excitation(ref)
        self.sing_energies = np.array(sing, dtype=float)
        self.sbvec = np.array(mol.data['OQP::td_bvec_mo']).copy()
        mol.data['OQP::td_singlet_energies'] = mol.data['OQP::td_energies'].copy()
        mol.data['OQP::td_bvec_mo_s'] = mol.data['OQP::td_bvec_mo'].copy()

        _select_response_manifold(mol, 3)
        trip = sp.excitation(ref)
        self.trip_energies = np.array(trip, dtype=float)
        self.tbvec = np.array(mol.data['OQP::td_bvec_mo']).copy()
        mol.data['OQP::td_triplet_energies'] = mol.data['OQP::td_energies'].copy()
        mol.data['OQP::td_bvec_mo_t'] = mol.data['OQP::td_bvec_mo'].copy()

        if is_dftb:
            _dftb_soc_tags(mol)
        else:
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
            _select_response_manifold(mol, mult)
            SinglePoint(mol).excitation([self.e_ref])
            mol.config['properties']['grad'] = [state]
            Gradient(mol).gradient()
            gi = np.array(mol.grads[state]).reshape((self.natom, 3))
            # ESPF_ROHF=1: use ROHF reference density for ESPF gradient across
            # all SOC MCH components, ensuring the ESPF energy is constant
            # across state hops.
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

        if getattr(self.driver, "espf_full", False):
            # Full-ESPF: pure MM-MM + nuclear-MM energy + analytic coupling force
            # (see NAMD_QMMM._total_force_espf / OpenQpQMMM._assemble_force_espf).
            nqm = len(self.qm_atoms)
            if len(self.driver.link_atoms) and g_qm.shape[0] == nqm:
                raise NotImplementedError(
                    "SOC-NAMD full-ESPF across a covalent boundary needs the QM "
                    "link atoms in the NAMD mol; use QMMM_MD for boundary QM/MM.")
            emm_q, gmm_q = self.driver.forces_mm(pchg)
            gmm = np.array(gmm_q.value_in_unit(u.kilojoule_per_mole / u.nanometer)) / HABOHR_TO_KJMOLNM
            emm = emm_q.value_in_unit(u.kilojoule_per_mole) * KJMOL_TO_HARTREE
            fq, fm, mm_idx = self.driver._coupling_forces(pchg)
            f_all = gmm.copy()
            for a, link in enumerate(self.driver.link_atoms):
                gl, fl = g_qm[nqm + a], fq[nqm + a]
                g_qm[link.host_row] = g_qm[link.host_row] + (1.0 - link.g) * gl
                fq[link.host_row] = fq[link.host_row] + (1.0 - link.g) * fl
                f_all[link.mm_index] = f_all[link.mm_index] - link.g * gl + link.g * fl
            for k, i in enumerate(self.qm_atoms):
                f_all[i] = f_all[i] - g_qm[k] + fq[k]
            for j, m in enumerate(mm_idx):
                f_all[m] = f_all[m] + fm[j]
            f_all -= f_all.mean(axis=0)
            eqm = float(e_diag) + float(
                np.dot(np.array(mol.get_atoms2("charge")), potmm))
            return f_all, eqm + emm

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
        self._prepare_md_outputs()

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
        self._log_soc_qmmm(
            0, epot, mult, state, w, False,
            self.e_ref + self.e0 + eval_ha)
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
            if self._prepare_hop_step(istep):
                hopped = NAMD_SOC._propagate_and_hop(self, self.prev_eval, eval_ha, t)
            else:
                hopped = False
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
            self._log_soc_qmmm(
                istep, epot, mult, state, w, hopped,
                self.e_ref + self.e0 + eval_ha)

        dump_log(mol, title='PyOQP: SOC-NAMD QM/MM trajectory complete')

    def _log_soc_qmmm(self, istep, epot, mult, state, w, hopped,
                      state_energies):
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
                   f'{self._hop_rng_log()}  '
                   f'dom=({NAMD_SOC._mch_label(mult, state)},w={w:.3f})  '
                   f'pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )
        self._unbiased_potential_energy = float(epot)
        self._trajectory_state_energies = np.asarray(
            state_energies, dtype=float).reshape(-1)
        self._write_md_trajectory(
            istep, self.r_all, epot, ekin, hopped)


class NAMD_SOC_MCH_QMMM(NAMD_SOC_QMMM):
    """QM/MM SOC-NAMD in the MCH basis with exact active-root QM gradient."""

    _resolve_initial_mch_active = NAMD_SOC_MCH._resolve_initial_mch_active
    _mch_active_label = NAMD_SOC_MCH._mch_active_label
    _mch_propagate_and_hop = NAMD_SOC_MCH._mch_propagate_and_hop
    _log_mch = NAMD_SOC_MCH._log_mch

    def __init__(self, mol):
        super().__init__(mol)
        self._trajectory_representation = 'soc_mch'
        self.coef = np.zeros(self.nstate_soc, dtype=complex)
        self.coef[self.active - 1] = 1.0 + 0.0j

    def _mch_exact_gradient_qmmm(self, active):
        mol = self.mol
        mult, state = self._mch_target(active - 1)
        _select_response_manifold(mol, mult)
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
        self._prepare_md_outputs()

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
        self._log_mch_qmmm(0, epot, mult, state, False, e_mch)
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
            if self._prepare_hop_step(istep):
                hopped = self._mch_propagate_and_hop(h_mch, e_mch)
            else:
                hopped = False
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
            self._log_mch_qmmm(
                istep, epot, mult, state, hopped, e_mch)

        dump_log(mol, title='PyOQP: SOC-MCH-QMMM-NAMD trajectory complete')

    def _log_mch_qmmm(self, istep, epot, mult, state, hopped,
                      state_energies):
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
                   f'{self._hop_rng_log()}  '
                   f'grad={NAMD_SOC._mch_label(mult, state)}  pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )
        self._unbiased_potential_energy = float(epot)
        self._trajectory_state_energies = np.asarray(
            state_energies, dtype=float).reshape(-1)
        self._write_md_trajectory(
            istep, self.r_all, epot, ekin, hopped)


def _dftb_soc_tags(mol):
    """Build OQP::soc_* tags for method=dftb/xtb (one-center SOC + numpy eigh)."""
    from oqp.library.openqp_dftb import HA_TO_WAVENUMBER, FINE_STRUCTURE
    adapter = make_tb_adapter(mol)
    data = mol.data
    dims = np.asarray(data['OQP::dftb_wf_dims']).ravel()
    nbf, noca, nocb = (int(round(v)) for v in dims[:3])
    x_s = np.asarray(data['OQP::td_bvec_mo_s'])
    x_t = np.asarray(data['OQP::td_bvec_mo_t'])
    hsoc_re, hsoc_im = adapter.soc_matrix(
        np.asarray(data['OQP::VEC_MO_A']).ravel(), x_s.ravel(), x_t.ravel(),
        noca=noca, nocb=nocb)
    e_s = np.asarray(data['OQP::td_singlet_energies']).ravel()
    e_t = np.asarray(data['OQP::td_triplet_energies']).ravel()
    e0 = min(e_s[0], e_t[0])
    diag = np.concatenate([e_s - e0, np.repeat(e_t - e0, 3)]) * HA_TO_WAVENUMBER
    dfac = 0.5 * FINE_STRUCTURE ** 2 * HA_TO_WAVENUMBER
    h_total = np.diag(diag).astype(complex) + (hsoc_re + 1j * hsoc_im) * dfac
    eigenvalues, eigenvectors = np.linalg.eigh(h_total)
    fortran_tag = adapter._fortran_tag
    data['OQP::soc_eval'] = np.ascontiguousarray(eigenvalues.real)
    data['OQP::soc_evec_re'] = fortran_tag(np.ascontiguousarray(eigenvectors.real))
    data['OQP::soc_evec_im'] = fortran_tag(np.ascontiguousarray(eigenvectors.imag))
    data['OQP::soc_hsoc_re'] = fortran_tag(np.ascontiguousarray(hsoc_re))
    data['OQP::soc_hsoc_im'] = fortran_tag(np.ascontiguousarray(hsoc_im))


def _dftb_spatial_overlap(mol, multiplicity):
    """TB (dftb/xtb) spatial state overlap for the current td_bvec_mo(_old) tags."""
    adapter = make_tb_adapter(mol)
    data = mol.data
    dims = np.asarray(data['OQP::dftb_wf_dims']).ravel()
    nbf, noca, nocb = (int(round(v)) for v in dims[:3])
    tlf = int(mol.config.get('tdhf', {}).get('tlf', 2))
    _, s_st = adapter.states_overlap(
        np.asarray(data['OQP::xyz_old']).ravel(),
        np.asarray(mol.get_system(), dtype=float).ravel(),
        np.asarray(data['OQP::VEC_MO_A_old']).ravel(),
        np.asarray(data['OQP::VEC_MO_A']).ravel(),
        np.asarray(data['OQP::td_bvec_mo_old']).ravel(),
        np.asarray(data['OQP::td_bvec_mo']).ravel(),
        noca=noca, nocb=nocb, multiplicity=multiplicity, tlf_order=tlf)
    data['OQP::td_states_overlap'] = s_st
