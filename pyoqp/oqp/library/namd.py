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
             -> TDC provider             (overlap FD/NPI or resident analytic
                                          sigma_ij = v dot d_ij)
             -> oqp.mrsf_namd_hop()      (RK4 amplitude propagation, EDC,
                                          trivial-crossing follow, FSSH hop +
                                          isotropic or analytic-NAC directional
                                          velocity rescaling)
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
from datetime import date
from importlib import resources
import numpy as np

from oqp.library.qmmm_active import (
    freeze_constrained_partners, held_atoms, resolve_active_set)

import oqp
from oqp.library.ints_1e import ints_1e
from oqp.library.single_point import SinglePoint, Gradient, LastStep, BasisOverlap, NACME, SCFnotConverged
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
ANGSTROM_TO_BOHR = 0.1 * NM_TO_BOHR
# 1 Hartree/bohr in kJ/mol/nm  (2625.499639 kJ/mol per Ha / 0.0529177 nm per bohr)
HABOHR_TO_KJMOLNM = 2625.499639 / BOHR_TO_NM
KJMOL_TO_HARTREE = 1.0 / 2625.499639
KCALMOL_TO_HARTREE = 1.0 / 627.5094740631
KCALMOLANG2_TO_HARTREEBOHR2 = (
    KCALMOL_TO_HARTREE / ANGSTROM_TO_BOHR**2
)
INT64_MIN = -(1 << 63)
INT64_MAX = (1 << 63) - 1
NAMD_RESTART_SCHEMA_VERSION = 8
NAMD_TRAJECTORY_SCHEMA_VERSION = 7
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


def _normalize_identity_value(value):
    """Convert parsed configuration values into deterministic JSON values."""
    if isinstance(value, dict):
        return {
            str(key): _normalize_identity_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_normalize_identity_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return _normalize_identity_value(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def _electronic_config_identity(config):
    """Capture every normalized option that can change NAMD surfaces."""
    input_keys = (
        'charge', 'basis', 'library', 'functional', 'method', 'ispher', 'd4',
        'soc_2e',
    )
    identity = {
        'input': {
            key: config.get('input', {}).get(key)
            for key in input_keys
        },
    }
    for section in (
            'scf', 'tdhf', 'dftgrid', 'pcm', 'dftb', 'xtb', 'symmetry'):
        identity[section] = config.get(section, {})
    identity['nac'] = {'align': config.get('nac', {}).get('align')}
    return _normalize_identity_value(identity)


def _validate_gate_tolerances(label, *values):
    """Reject NaN/Inf as well as negative validation tolerances."""
    if len(values) == 1 and not np.isscalar(values[0]):
        values = tuple(values[0])
    tolerances = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(tolerances)) or np.any(tolerances < 0.0):
        raise ValueError(
            f"[md] {label} gate tolerances must be finite and non-negative")


def _validate_thermostat_parameters(temperature, friction, enabled):
    """Validate Langevin parameters, including a nonzero coupling rate."""
    values = np.asarray((temperature, friction), dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError(
            '[md] thermostat_temperature and thermostat_friction must be '
            'finite and non-negative')
    if enabled and friction <= 0.0:
        raise ValueError(
            '[md] thermostat_friction must be positive for Langevin NVT')


def _validate_nacme_gate_activation(check, gate):
    """Reject a policy that cannot observe any NACME reference diagnostic."""
    if check == 'off' and gate != 'off':
        raise ValueError(
            "[md] nacme_gate must be off when nacme_check is off; enable "
            "nacme_check=baeck_an before selecting warn or error")


def _resolve_namd_seed(value, current_date=None):
    """Resolve the zero default to a reproducible YYYYMMDD run seed."""
    seed = int(value)
    if seed != 0:
        return seed
    run_date = date.today() if current_date is None else current_date
    return int(run_date.strftime('%Y%m%d'))


def _validate_distinct_output_paths(*, protected_paths=(), **paths):
    """Reject NAMD outputs that overwrite one another or an input deck."""
    def aliases(first, second):
        first = os.fspath(first)
        second = os.fspath(second)
        if (os.path.normcase(os.path.realpath(first))
                == os.path.normcase(os.path.realpath(second))):
            return True
        try:
            return os.path.samefile(first, second)
        except (FileNotFoundError, OSError):
            return False

    outputs = list(paths.items())
    collisions = []
    for index, (left_label, left_path) in enumerate(outputs):
        for right_label, right_path in outputs[index + 1:]:
            if aliases(left_path, right_path):
                collisions.append((left_label, right_label))
    if collisions:
        rendered = "; ".join(" = ".join(labels) for labels in collisions)
        raise ValueError(
            "[md] NAMD output paths must be distinct; collision: " + rendered)

    protected = [path for path in protected_paths if path]
    for label, output_path in outputs:
        if any(aliases(output_path, input_path) for input_path in protected):
            raise ValueError(
                f"[md] NAMD output {label} must not alias the input deck")


_ANALYTIC_NAC_CONV_MAX = 1.0e-8

_TRAJECTORY_CONTROL_KEYS = (
    'mo_reuse', 'scf_fail', 'scf_guess_retry', 'ref_follow',
    'ref_switch_rescale', 'somo_tol', 'frustrated', 'disc_rescale',
    'disc_tol', 'disc_substeps')


def _normalized_md_control(key, value):
    """Canonical text of an [md] control, typed by its schema entry."""
    from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
    kind = OQP_CONFIG_SCHEMA['md'][key]['type']
    text = str(value).strip().lower()
    try:
        if kind is bool:
            return 'true' if text in ('true', '1', 'on', 'yes') else 'false'
        if kind is int:
            return repr(int(float(text)))
        if kind is float:
            return repr(float(text))
    except ValueError:
        pass
    return text


def _non_default_md_controls(md, keys):
    """Return {key: canonical value} for controls that differ from defaults."""
    from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
    controls = {}
    for key in keys:
        if key not in md:
            continue
        value = _normalized_md_control(key, md[key])
        default = _normalized_md_control(key, OQP_CONFIG_SCHEMA['md'][key]['default'])
        if value != default:
            controls[key] = value
    return controls


def _config_flag(value):
    return (value is True) or (str(value).strip().lower() in ('true', '1', 'on', 'yes'))


def analytic_nac_model_issue(config):
    """Return why the resident analytic MRSF NAC cannot serve this model.

    ``None`` means the electronic model is the two-SOMO ROHF/ROKS MRSF singlet
    response that ``oqp.library.nac_analytic.analytic_nac`` implements, with
    SCF and response thresholds tight enough for its accuracy guard.  The
    spin-orbit and QM/MM boundaries are enforced separately by their drivers.
    """
    inp = config.get('input', {})
    scf = config.get('scf', {})
    tdhf = config.get('tdhf', {})
    if is_tb_method(str(inp.get('method', ''))):
        return 'tight-binding models have no analytic NAC'
    if str(tdhf.get('type', '')).strip().lower() != 'mrsf':
        return 'tdhf.type must be mrsf'
    try:
        tdhf_mult = int(tdhf.get('multiplicity', 1))
        scf_mult = int(scf.get('multiplicity', 1))
    except (TypeError, ValueError):
        return 'invalid multiplicity'
    if tdhf_mult != 1:
        return 'analytic NAC implements singlet MRSF states only'
    if str(scf.get('type', '')).strip().lower() != 'rohf' or scf_mult != 3:
        return 'analytic NAC requires a two-SOMO ROHF/ROKS triplet reference'
    for section, values in (('scf', scf), ('tdhf', tdhf)):
        try:
            conv = float(values.get('conv', 1.0e-6))
        except (TypeError, ValueError):
            return '%s.conv is not a number' % section
        if not (np.isfinite(conv) and 0.0 < conv <= _ANALYTIC_NAC_CONV_MAX):
            return '%s.conv=%g exceeds %g' % (section, conv, _ANALYTIC_NAC_CONV_MAX)
    return None


def analytic_nac_route_issue(config):
    """Return why a NAMD route cannot use analytic NAC, including SOC/QM/MM."""
    if bool(config.get('md', {}).get('soc', False)):
        return 'spin-orbit NAMD'
    if _config_flag(config.get('input', {}).get('qmmm_flag', False)):
        return 'QM/MM NAMD'
    return analytic_nac_model_issue(config)


def _sha256_stream(stream):
    """Hash a binary stream in bounded chunks, including on Python 3.10."""
    digest = hashlib.sha256()
    for chunk in iter(lambda: stream.read(1024 * 1024), b''):
        digest.update(chunk)
    return digest.hexdigest()


def _restart_manifest_path(log_path):
    """Return a job-specific runnable restart manifest beside the main log."""
    absolute = os.path.abspath(os.fspath(log_path))
    stem = os.path.splitext(os.path.basename(absolute))[0]
    return os.path.join(os.path.dirname(absolute), stem + '.namd.restart.oqp')


def _validate_odp_boundary_conditions(odp, periodic):
    """Reject periodic ODP until native CVs implement minimum images."""
    if odp is not None and bool(periodic):
        raise NotImplementedError(
            "[odp] periodic QM/MM is unsupported because ODP CVs do not yet "
            "apply the minimum-image convention; use qmmm.cutoff=NoCutoff")


def _namd_trajectory_dtype(nstate, natom, ncv=0):
    """Fixed-width, appendable binary record used by ``*.namd.trj``."""
    matrix = (nstate, nstate)
    vectors = (natom, 3)
    return np.dtype([
        ('step', '<i8'), ('time_fs', '<f8'), ('active', '<i4'),
        ('hopped', 'i1'), ('rng', '<f8'),
        ('e_unbiased_pot_hartree', '<f8'), ('e_pot_hartree', '<f8'),
        ('e_kin_hartree', '<f8'),
        ('e_tot_hartree', '<f8'),
        ('droplet_energy_hartree', '<f8'),
        ('droplet_max_penetration_bohr', '<f8'),
        ('droplet_active_count', '<i8'),
        ('solute_com_energy_hartree', '<f8'),
        ('solute_com_displacement_bohr', '<f8'),
        ('conservative_restraint_energy_hartree', '<f8'),
        ('thermostat_exchange_hartree', '<f8'),
        ('thermostat_exchange_cumulative_hartree', '<f8'),
        ('thermostat_adjusted_energy_hartree', '<f8'),
        ('droplet_force_hartree_per_bohr', '<f8', vectors),
        ('solute_com_force_hartree_per_bohr', '<f8', vectors),
        ('state_energies', '<f8', (nstate,)),
        ('populations', '<f8', (nstate,)), ('coef_real', '<f8', (nstate,)),
        ('coef_imag', '<f8', (nstate,)), ('coordinates_bohr', '<f8', vectors),
        ('velocities_au', '<f8', vectors), ('state_overlap', '<f8', matrix),
        ('state_overlap_imag', '<f8', matrix),
        ('overlap_tdc_au', '<f8', matrix),
        ('overlap_tdc_imag_au', '<f8', matrix),
        ('gate_candidate_tdc_au', '<f8', matrix),
        ('reference_tdc_au', '<f8', matrix),
        ('reference_mask', 'u1', matrix), ('reference_source', 'i1'),
        ('tdc_source', 'i1'), ('rescale_source', 'i1'),
        ('rescale_gamma', '<f8'), ('rescale_discriminant', '<f8'),
        ('hop_direction', '<f8', vectors),
        ('gate_center_step', '<i8'), ('gate_verdict', 'i1'),
        ('gate_counts', '<i8', (3,)), ('gate_streak', '<i8'),
        ('gate_metrics', '<f8', (7,)),
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
        ('tracking_raw_order', '<i8', (nstate,)),
        ('tracking_lineage', '<i8', (nstate,)),
        ('tracking_phase', '<f8', (nstate,)),
        ('tracking_phase_initial', '<f8', (nstate,)),
        ('tracking_previous_phase_initial', '<f8', (nstate,)),
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
    if count == 0:
        return header, np.empty(0, dtype=dtype)
    records = np.memmap(path, dtype=dtype, mode=mmap_mode,
                        offset=offset, shape=(count,))
    return header, records


def read_odp_wham_series(path):
    """Extract the complete per-window ODP series from one packed NAMD TRJ."""
    header, records = read_namd_trajectory(path)
    snapshot_bytes = int(records.offset + records.nbytes)
    provenance = header.get('odp', {})
    if not provenance.get('enabled', False):
        raise ValueError('trajectory does not contain an enabled ODP umbrella')
    try:
        restart_identity = json.loads(header['signature'])
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError(
            'trajectory does not contain a valid NAMD system identity') from exc
    identity_keys = (
        'method', 'charge', 'functional', 'basis', 'scf_type',
        'scf_multiplicity', 'tdhf_type', 'tdhf_multiplicity', 'nstate', 'tlf',
        'd4', 'pcm', 'electronic_config', 'trajectory_representation',
    )
    system_identity = {
        key: restart_identity.get(key) for key in identity_keys
    }
    independent_controls = restart_identity.get('independent_controls', {})
    conservative_restraints = {}
    for name in ('droplet', 'solute_com'):
        control = independent_controls.get(name, {})
        conservative_restraints[name] = (
            control if bool(control.get('enabled', False))
            else {'enabled': False}
        )
    system_identity['conservative_restraints'] = conservative_restraints
    system = header.get('wham_system_identity')
    system_identity['system'] = system
    if (not isinstance(system, dict) or not system.get('sha256')
            or system_identity.get('method') is None
            or system_identity.get('nstate') is None):
        raise ValueError('trajectory has an incomplete NAMD system identity')
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
        'system_identity': system_identity,
        'ensemble': header.get('ensemble'),
        'thermostat_temperature_kelvin': (
            header.get('independent_controls', {})
            .get('thermostat', {}).get('temperature_kelvin')
        ),
        'snapshot_bytes': snapshot_bytes,
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
_P_ALLOW_HOP = 13       # +1 permit state changes; -1 propagate coefficients only
_P_RESCALE = 14         # 0 isotropic; 1 resident analytic NAC; 2 hop-triggered NAC
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
        # Freeze the user-resolved electronic configuration before SOC/TB
        # helpers temporarily switch response manifolds during propagation.
        self._electronic_config_identity = _electronic_config_identity(cfg)

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
        self.tdc_provider = str(md['tdc']).strip().lower().replace('-', '_')
        if self.tdc_provider in ('ba', 'tdba'):
            self.tdc_provider = 'baeck_an'
        if self.tdc_provider not in ('fd', 'npi', 'analytic', 'baeck_an'):
            raise ValueError(
                "[md] tdc must be fd, npi, analytic, or baeck_an")
        self.tdc_scheme = {
            'fd': 0, 'npi': 1, 'analytic': 2, 'baeck_an': 3,
        }[self.tdc_provider]
        self.rescale_provider = str(md.get('rescale', 'auto')).strip().lower().replace('-', '_')
        if self.rescale_provider in ('analytic', 'nac'):
            self.rescale_provider = 'analytic_nac'
        if self.rescale_provider in ('hop_analytic', 'hop_nac', 'ht_nac'):
            self.rescale_provider = 'hop_analytic_nac'
        if self.rescale_provider not in (
                'auto', 'isotropic', 'analytic_nac', 'hop_analytic_nac'):
            raise ValueError(
                "[md] rescale must be auto, isotropic, analytic_nac, or "
                "hop_analytic_nac")
        # rescale=auto selects hop-triggered analytic NAC only where the
        # resident analytic NAC is defined.  Resolving this once at startup
        # keeps every other route on isotropic rescaling instead of failing at
        # the first stochastic hop candidate, possibly deep into a trajectory.
        self._rescale_auto_issue = None
        if self.rescale_provider == 'auto':
            self._rescale_auto_issue = analytic_nac_route_issue(cfg)
            self.rescale_provider = (
                'hop_analytic_nac' if self._rescale_auto_issue is None
                else 'isotropic')
        # Carry the converged orbitals of the previous geometry into the SCF
        # of the next geometry (guess type 'previous' after the first step).
        # Reuse the previous orbitals by default. Repeating the configured
        # guess at every step can select a different ROHF solution or
        # orbital ordering and collapse the electronic-state overlap.
        self.mo_reuse = str(md.get('mo_reuse', 'true')).strip().lower() in (
            'true', '1', 'on', 'yes')
        self._overlap_collapse_steps = 0
        # scf_fail=restart: no in-step converger escalation.  When the primary
        # converger fails from the previous-step orbitals, re-solve the
        # reference from a fresh guess with SOSCF (the KNU-GAMESS restart
        # procedure) and treat the step as a restart boundary: the electronic
        # coefficients are frozen and no hop is attempted for that step, so a
        # bra from the old SCF branch is never combined with a ket from the new
        # one.  scf_fail=escalate keeps the SinglePoint SOSCF/TRAH ladder.
        self.scf_fail = str(md.get('scf_fail', 'escalate')).strip().lower()
        if self.scf_fail not in ('escalate', 'restart'):
            raise ValueError("[md] scf_fail must be escalate or restart")
        self.scf_guess_retry = str(md.get('scf_guess_retry', 'true')).strip().lower() in (
            'true', '1', 'on', 'yes')
        self._restart_boundary = False
        self._scf_restart_steps = 0
        # Reference (SOMO) continuity controls.  ref_follow selects the SCF
        # continuation converger for steps after the first: SOSCF or DIIS
        # with a 0.2 Hartree level shift both keep the previous-step SOMO
        # configuration where plain C-DIIS can jump to a different ROHF
        # triplet configuration.  The SOMO block of the aligned MO overlap
        # detects a configuration change (reference switch event).  With
        # ref_switch_rescale the velocities are rescaled isotropically at such
        # a step so that the total energy is conserved across the jump of the
        # active-state MRSF energy, and the jump is recorded.
        self.ref_follow = str(md.get('ref_follow', 'soscf')).strip().lower().replace('-', '_')
        if self.ref_follow not in ('off', 'soscf', 'diis_vshift'):
            raise ValueError("[md] ref_follow must be off, soscf, or diis_vshift")
        self.ref_switch_rescale = str(md.get('ref_switch_rescale', 'true')).strip().lower() in (
            'true', '1', 'on', 'yes')
        self.somo_tol = float(md.get('somo_tol', 0.5))
        # Frustrated-hop treatment for derivative-coupling (directional)
        # rescaling: 'none' leaves the velocity unchanged (Tully 1990);
        # 'reflect' reverses the momentum component along d_IJ
        # (Hammes-Schiffer & Tully 1994).
        self.frustrated = str(md.get('frustrated', 'reflect')).strip().lower()
        if self.frustrated not in ('none', 'reflect'):
            raise ValueError("[md] frustrated must be none or reflect")
        self._frustrated_reflect_count = 0
        # Numerical energy correction is an optional last resort after
        # finer nuclear integration, separate from a physical surface hop.
        self.disc_rescale = str(md.get('disc_rescale', 'true')).strip().lower() in (
            'true', '1', 'on', 'yes')
        self.disc_tol = float(md.get('disc_tol', 0.002))
        if not np.isfinite(self.disc_tol) or self.disc_tol <= 0.0:
            raise ValueError("[md] disc_tol must be positive")
        self._disc_event_count = 0
        self._disc_energy_absorbed = 0.0
        self._step_numerical_correction = 0.0
        # Energy-guarded nuclear substepping: when the pre-hop total-energy
        # jump of a step exceeds disc_tol, the step is repeated from the
        # previous phase point and electronic state with disc_substeps
        # velocity-Verlet substeps (electronic structure and active-state
        # force at every substep, orbitals carried along).  The state overlap,
        # couplings and hop decision are then evaluated once between the
        # start and the end of the full step as usual.  Only the residual jump
        # remaining after refinement may be treated by disc_rescale.
        self.disc_substeps = int(md.get('disc_substeps', 10))
        if self.disc_substeps < 0:
            raise ValueError("[md] disc_substeps must be >= 0")
        if self.disc_rescale or self.ref_switch_rescale:
            self.disc_substeps = max(2, self.disc_substeps)
        self._disc_substep_events = 0
        self._window_leak_step = False
        self._somo_switch_step = False
        self._somo_switch_count = 0
        self._window_leak_count = 0
        self._scf_fallback_steps = 0
        self._etot_prev = None
        self._ref_switch_jump = np.nan
        self.trivial = 1 if str(md['trivial']).lower() in ('true', '1', 'on', 'yes') else 0
        self.trivial_thresh = float(md['trivial_thresh'])
        self.init_temp = float(md['init_temp'])
        self.seed = _resolve_namd_seed(md.get('seed', 0))
        self.rng_stream = int(md.get('rng_stream', 1))
        self.first_hop_step = int(md.get('first_hop_step', 1))
        self.nacme_check = str(md.get(
            'nacme_check', 'off')).strip().lower().replace('-', '_')
        if self.nacme_check == 'tdba':
            self.nacme_check = 'baeck_an'
        if self.nacme_check not in ('off', 'baeck_an', 'analytic'):
            raise ValueError("[md] nacme_check must be off, baeck_an, or analytic")
        self.ba_gap_max = float(md.get('ba_gap_max', 0.0734986443513))
        if not np.isfinite(self.ba_gap_max) or self.ba_gap_max <= 0.0:
            raise ValueError("[md] ba_gap_max must be positive and finite")
        self.nacme_gate = str(md.get('nacme_gate', 'off')).strip().lower()
        if self.nacme_gate not in ('off', 'warn', 'error'):
            raise ValueError("[md] nacme_gate must be off, warn, or error")
        _validate_nacme_gate_activation(self.nacme_check, self.nacme_gate)
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
        self.nve_gate = str(md.get('nve_gate', 'warn')).strip().lower()
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
        self.ensemble = str(md.get('ensemble', 'nve')).strip().lower()
        self.thermostat = str(md.get('thermostat', 'off')).strip().lower()
        if self.ensemble not in ('nve', 'nvt'):
            raise ValueError("[md] ensemble must be nve or nvt")
        if self.thermostat not in ('off', 'langevin'):
            raise ValueError("[md] thermostat must be off or langevin")
        if self.ensemble == 'nve' and self.thermostat != 'off':
            raise ValueError(
                "[md] thermostat is independent of NVE; use ensemble=nvt or thermostat=off"
            )
        if self.ensemble == 'nvt' and self.thermostat == 'off':
            raise ValueError("[md] ensemble=nvt requires thermostat=langevin")
        if self.ensemble == 'nvt' and self.nve_gate == 'warn':
            # ``warn`` is the NVE safety default.  NVT records thermostat work
            # instead, so an inherited default becomes contextually inactive.
            self.nve_gate = 'off'
        if self.ensemble == 'nvt' and self.nve_gate != 'off':
            raise ValueError(
                "[md] nve_gate does not apply to NVT; thermostat exchange is recorded separately"
            )
        self.thermostat_temperature = float(
            md.get('thermostat_temperature', self.init_temp))
        self.thermostat_friction = float(md.get('thermostat_friction', 1.0))
        _validate_thermostat_parameters(
            self.thermostat_temperature, self.thermostat_friction,
            self.thermostat == 'langevin')
        self.trajectory_interval_input = int(md.get('trajectory_interval', 1))
        self.restart_interval_input = int(md.get('restart_interval', 10))
        self.trajectory_interval = self._output_interval_steps(
            self.trajectory_interval_input, self.dt_fs)
        self.restart_interval = self._output_interval_steps(
            self.restart_interval_input, self.dt_fs)
        self.restart_requested = self._as_bool(md.get('restart', False))
        self.continuation_checkpoint = str(md.get('continuation_checkpoint', '')).strip()
        self.continuation_trajectory = str(md.get('continuation_trajectory', '')).strip()
        if bool(self.continuation_checkpoint) != bool(self.continuation_trajectory):
            raise ValueError('continuation_checkpoint and continuation_trajectory are both required')
        if self.continuation_checkpoint:
            if self.restart_requested or type(self) is not NAMD:
                raise ValueError('local continuation requires a new same-spin NAMD run')
            if self.tdc_provider != 'analytic':
                raise ValueError('local continuation currently requires analytic TDC')
            self.continuation_checkpoint = os.path.abspath(os.path.expanduser(self.continuation_checkpoint))
            self.continuation_trajectory = os.path.abspath(os.path.expanduser(self.continuation_trajectory))
        self._time_origin_fs = 0.0
        self._continuation_provenance = None
        self.trajectory_file = self._md_output_path(
            md.get('trajectory_file', ''), '.namd.trj')
        # Observational sidecar for the transported/extrapolated NAC adjoint.
        # Its path is deliberately derived rather than added to the public
        # input schema while the approximation remains experimental.
        self.zpredict_audit_file = self._md_output_path(
            '', '.namd.zpredict.tsv')
        self.restart_file = self._md_output_path(
            md.get('restart_file', ''), '.namd.restart.npz')
        self.restart_manifest_file = self._restart_manifest_path()
        self._restart_manifest_written = False
        self.velocity_source = str(md['velocity'])
        self._validate_sidecar_paths()
        # Capture external guess inputs before the first electronic step. A
        # save_mol target can be rewritten at every geometry, but that output
        # mutation must not change trajectory or checkpoint identity.
        self._restart_guess_identity = self._guess_settings_identity()
        self.odp = odp_from_config(cfg)
        self._odp_last = None
        self._unbiased_potential_energy = np.nan
        _soc = md.get('soc', False)
        # Keep this guard identical to compute_namd's ``if soc:`` dispatch.
        # Programmatic callers are not necessarily constrained by the input
        # schema and may supply another truthy value.
        soc_requested = bool(_soc)
        if (soc_requested and self.nacme_check == 'baeck_an'
                and self.nacme_gate == 'off'):
            # Baeck-An is a real same-spin magnitude diagnostic.  SOC records
            # the full complex spin-adiabatic overlap and anti-Hermitian TDC.
            self.nacme_check = 'off'
        if soc_requested and (
                self.nacme_check != 'off'
                or self.tdc_provider in ('analytic', 'baeck_an')
                or self.rescale_provider in (
                    'analytic_nac', 'hop_analytic_nac')):
            raise NotImplementedError(
                "analytic NAC TDC/rescaling/check currently supports same-spin NAMD only"
            )
        if (self._needs_analytic_nac()
                or self.rescale_provider == 'hop_analytic_nac'):
            model_issue = analytic_nac_model_issue(cfg)
            if model_issue is not None:
                raise ValueError(
                    "[md] analytic NAC TDC/rescaling/check was requested, but "
                    "this electronic model cannot provide it: %s. Use "
                    "rescale=auto or isotropic, or tighten the model." % model_issue)
        if soc_requested and self.odp is not None:
            raise NotImplementedError(
                "[odp] currently supports same-spin NVE NAMD only"
            )
        if self.odp is not None and self.ensemble != 'nve':
            raise NotImplementedError(
                "[odp] currently supports same-spin NVE NAMD only"
            )
        if self.dt_adaptive and not soc_requested:
            raise NotImplementedError(
                "[md] dt_adaptive currently supports SOC-NAMD only"
            )
        if not INT64_MIN <= self.seed <= INT64_MAX:
            raise ValueError("[md] seed must fit in a signed 64-bit integer")
        if not 0 <= self.rng_stream <= INT64_MAX:
            raise ValueError(
                "[md] rng_stream must be a non-negative signed 64-bit integer"
            )
        if self.first_hop_step < 1:
            raise ValueError("[md] first_hop_step must be at least 1")
        self.natom = mol.data['natom']
        if self.odp is not None and not self._as_bool(
                cfg.get('input', {}).get('qmmm_flag', False)):
            self.odp.validate_atom_count(self.natom)
        # get_mass() returns atomic masses in amu; the integrator works in
        # atomic units, so convert to electron masses.
        self.mass = mol.get_mass() * AMU_TO_AU     # (natom,) electron masses
        self._restart_system_identity = self._qm_restart_system_identity()
        self._wham_system_identity = self._qm_wham_system_identity()
        self._init_independent_controls(cfg)
        if soc_requested and (self.droplet_enabled or self.solute_com_enabled
                              or self.thermostat != 'off'):
            raise NotImplementedError(
                "droplet/solute_com/NVT controls currently support same-spin NAMD only"
            )
        if (self.droplet_enabled or self.solute_com_enabled) and not self._as_bool(
                cfg.get('input', {}).get('qmmm_flag', False)):
            self._setup_gas_restraint_targets()
        self._rng_step = 0
        self._last_hop_random = np.nan
        self._last_hop_probabilities = None
        self._ba_energy_left = None
        self._ba_energy_center = None
        self._ba_tdc_left = None
        self._ba_dt_left = None
        self._ba_last = None
        self._last_baeck_an_tdc = None
        # Baeck-An needs three energy points.  Its first interval and every
        # reseeded interval therefore use the phase-tracked overlap TDC.
        self._last_tdc_source = (
            1 if self.tdc_provider == 'baeck_an' else self.tdc_scheme)
        self._nacme_gate_failures = 0
        self._nacme_gate_last = None
        self._pending_nacme_gate_error = None
        self._nacme_candidate_tdc = None
        self._nacme_reference_tdc = None
        self._nacme_reference_mask = None
        self._nacme_reference_source = 0
        self._last_state_overlap = None
        self._last_overlap_tdc = None
        self._last_analytic_dcv = None
        self._last_analytic_tdc = None
        self._last_analytic_pair = None
        self._last_analytic_step = None
        self._analytic_tdc_previous = None
        self._analytic_tdc_centered = None
        self._last_rescale_source = {
            'isotropic': 0, 'analytic_nac': 1,
            'hop_analytic_nac': 2,
        }[self.rescale_provider]
        self._last_rescale_gamma = np.nan
        self._last_rescale_discriminant = np.nan
        self._last_hop_direction = np.zeros((self.natom, 3), dtype=float)
        self._nve_reference_energy = None
        self._nve_previous_energy = None
        self._nve_gate_failures = 0
        self._nve_gate_last = None
        self._pending_nve_gate_error = None
        self._thermostat_exchange = 0.0
        self._thermostat_exchange_cumulative = 0.0
        self._pending_nacme_gate_error = None
        self._trajectory_prefix_hasher = None
        self._trajectory_prefix_bytes = 0
        self._trajectory_prefix_last_step = None
        self._trajectory_prefix_stat = None

        # electronic amplitudes (complex), one per excited state. For SOC-NAMD
        # the active index runs over the larger spin-adiabatic manifold and the
        # subclass overwrites coef; guard the base indexing against that.
        self.coef = np.zeros(self.nstate, dtype=complex)
        if 1 <= self.active <= self.nstate:
            self.coef[self.active - 1] = 1.0 + 0.0j
        else:
            self.coef[0] = 1.0 + 0.0j

        # velocities (natom, 3) in atomic units
        # A checkpoint is the authoritative velocity state.  Do not retain a
        # transport-time dependency on the original velocity input file.
        self.vel = (np.zeros((self.natom, 3), dtype=float)
                    if self.restart_requested or self.continuation_checkpoint else self._init_velocities())

        # previous-step payload for the overlap (back_door carry)
        self.prev_xyz = None
        self.prev_data = None

        # force the per-step electronic pipeline to use the in-memory previous
        # step rather than recomputing it
        cfg['properties']['back_door'] = True
        # NACME needs a dt; reuse the MD step (atomic units) for the TDC scale
        cfg['nac']['dt'] = self.dt

    @staticmethod
    def _output_interval_steps(configured, dt_fs):
        """Resolve zero to an approximately 10 fs fixed-step output cadence."""
        configured = int(configured)
        dt_fs = float(dt_fs)
        if configured < 0:
            raise ValueError(
                "[md] trajectory/restart intervals must be zero or positive")
        if not np.isfinite(dt_fs) or dt_fs <= 0.0:
            raise ValueError("[md] dt must be finite and positive")
        return configured or max(1, int(round(10.0 / dt_fs)))

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

    def _qm_wham_system_identity(self):
        """Bind WHAM windows to atom identity while allowing new coordinates."""
        digest = _restart_identity_digest(array_parts=(
            ('atomic_numbers', self.mol.get_atoms(), '<i8'),
            ('masses_electron', self.mass, '<f8'),
        ))
        return {'kind': 'qm', 'natom': int(self.natom), 'sha256': digest}

    @staticmethod
    def _vector3(value, label):
        if isinstance(value, str):
            values = [float(item) for item in value.replace(',', ' ').split()]
        else:
            values = [float(item) for item in value]
        if len(values) != 3 or not np.all(np.isfinite(values)):
            raise ValueError(f"{label} must contain three finite Cartesian values")
        return np.asarray(values, dtype=np.float64)

    def _init_independent_controls(self, cfg):
        """Parse droplet, solute-COM, and thermostat records once.

        User-facing droplet lengths are angstrom and force constants are
        kcal mol^-1 angstrom^-2.  Everything below this method is atomic units.
        Neither control is inferred from ODP, QM/MM, or the NVT thermostat.
        """
        droplet = cfg.get('droplet', {})
        self.droplet_enabled = self._as_bool(droplet.get('enabled', False))
        self.droplet_center_angstrom = self._vector3(
            droplet.get('center', (0.0, 0.0, 0.0)), '[droplet] center')
        self.droplet_center = (
            self.droplet_center_angstrom * ANGSTROM_TO_BOHR)
        self.droplet_radius_angstrom = float(droplet.get('radius', 20.0))
        self.droplet_buffer_angstrom = float(droplet.get('buffer', 1.0))
        self.droplet_force_constant_input = float(
            droplet.get('force_constant', 10.0))
        self.droplet_max_penetration_angstrom = float(
            droplet.get('max_penetration', 10.0))
        self.droplet_radius = self.droplet_radius_angstrom * ANGSTROM_TO_BOHR
        self.droplet_buffer = self.droplet_buffer_angstrom * ANGSTROM_TO_BOHR
        self.droplet_force_constant = (
            self.droplet_force_constant_input *
            KCALMOLANG2_TO_HARTREEBOHR2)
        self.droplet_max_penetration = (
            self.droplet_max_penetration_angstrom * ANGSTROM_TO_BOHR)
        self.droplet_target = str(
            droplet.get('target', 'water_com')).strip().lower().replace('-', '_')
        if self.droplet_target in ('com', 'molecule_com', 'water_molecule_com'):
            self.droplet_target = 'water_com'
        if self.droplet_target in ('o', 'water_oxygen'):
            self.droplet_target = 'oxygen'
        if self.droplet_target not in ('water_com', 'oxygen', 'atoms'):
            raise ValueError(
                "[droplet] target must be water_com, oxygen, or atoms")
        self.droplet_atoms_spec = str(droplet.get('atoms', '') or '').strip()
        names = droplet.get(
            'water_resnames', ('hoh', 'wat', 'sol', 'tip3', 'tip3p'))
        if isinstance(names, str):
            names = names.replace(',', ' ').split()
        self.droplet_water_resnames = tuple(
            str(name).strip().lower() for name in names if str(name).strip())
        if self.droplet_enabled:
            values = (self.droplet_radius_angstrom,
                      self.droplet_buffer_angstrom,
                      self.droplet_force_constant_input,
                      self.droplet_max_penetration_angstrom)
            if not np.all(np.isfinite(values)):
                raise ValueError("[droplet] numeric settings must be finite")
            if (self.droplet_radius_angstrom <= 0.0
                    or self.droplet_buffer_angstrom < 0.0
                    or self.droplet_force_constant_input <= 0.0
                    or self.droplet_max_penetration_angstrom < 0.0):
                raise ValueError(
                    "[droplet] radius/force_constant must be positive and "
                    "buffer/max_penetration non-negative")
        self._droplet_group_index = None
        self._droplet_group_count = 0
        self._droplet_energy = 0.0
        self._droplet_max_penetration = 0.0
        self._droplet_active_count = 0
        self._droplet_force = None
        self._droplet_force_max = 0.0

        solute = cfg.get('solute_com', {})
        self.solute_com_enabled = self._as_bool(solute.get('enabled', False))
        self.solute_com_center_angstrom = self._vector3(
            solute.get('center', (0.0, 0.0, 0.0)), '[solute_com] center')
        self.solute_com_center = (
            self.solute_com_center_angstrom * ANGSTROM_TO_BOHR)
        self.solute_com_force_constant_input = float(
            solute.get('force_constant', 5.0))
        self.solute_com_force_constant = (
            self.solute_com_force_constant_input *
            KCALMOLANG2_TO_HARTREEBOHR2)
        self.solute_com_atoms_spec = str(solute.get('atoms', '') or '').strip()
        if (self.solute_com_enabled
                and (not np.isfinite(self.solute_com_force_constant_input)
                     or self.solute_com_force_constant_input <= 0.0)):
            raise ValueError("[solute_com] force_constant must be positive and finite")
        self._solute_com_selected = None
        self._solute_com_energy = 0.0
        self._solute_com_displacement = 0.0
        self._solute_com_force = None
        self._conservative_restraint_energy = 0.0
        self._conservative_restraint_force = None

    @staticmethod
    def _single_atom_groups(natom, atom_indices, label):
        indices = np.asarray(sorted(set(int(i) for i in atom_indices)), dtype=int)
        if indices.size == 0:
            raise ValueError(f"{label} selected no atoms")
        if indices[0] < 0 or indices[-1] >= int(natom):
            raise ValueError(
                f"{label} atom indices must be zero-based and within 0..{int(natom)-1}")
        groups = np.zeros(int(natom), dtype=np.int64)
        groups[indices] = np.arange(1, len(indices) + 1, dtype=np.int64)
        return groups

    def _setup_gas_restraint_targets(self):
        if self.droplet_enabled:
            if self.droplet_target == 'water_com':
                raise ValueError(
                    "[droplet] target=water_com requires qmmm(...) topology; "
                    "use target=oxygen or target=atoms for all-QM NAMD")
            if self.droplet_target == 'oxygen':
                indices = np.flatnonzero(
                    np.asarray(self.mol.get_atoms(), dtype=int) == 8)
            else:
                indices = _parse_int_list(self.droplet_atoms_spec)
            self._droplet_group_index = self._single_atom_groups(
                self.natom, indices, '[droplet]')
            self._droplet_group_count = int(self._droplet_group_index.max())
        if self.solute_com_enabled:
            if self.solute_com_atoms_spec:
                indices = _parse_int_list(self.solute_com_atoms_spec)
            else:
                indices = range(self.natom)
            selected = self._single_atom_groups(
                self.natom, indices, '[solute_com]')
            self._solute_com_selected = (selected > 0).astype(np.int64)

    def _setup_qmmm_restraint_targets(self):
        if self.droplet_enabled:
            if self.periodic:
                raise ValueError(
                    "[droplet] finite spherical containment requires qmmm(cutoff=NoCutoff)"
                )
            if self.droplet_target == 'atoms':
                groups = self._single_atom_groups(
                    self.natom_all, _parse_int_list(self.droplet_atoms_spec),
                    '[droplet]')
            else:
                groups = np.zeros(self.natom_all, dtype=np.int64)
                group = 0
                for residue in self.pdb.topology.residues():
                    if str(residue.name).strip().lower() not in self.droplet_water_resnames:
                        continue
                    atoms = list(residue.atoms())
                    if self.droplet_target == 'oxygen':
                        atoms = [atom for atom in atoms if (
                            getattr(getattr(atom, 'element', None), 'symbol', '') == 'O'
                            or str(atom.name).strip().upper().startswith('O'))]
                        atoms = atoms[:1]
                    if not atoms:
                        continue
                    group += 1
                    for atom in atoms:
                        groups[int(atom.index)] = group
                if group == 0:
                    raise ValueError(
                        "[droplet] found no target waters; check target and water_resnames")
            self._droplet_group_index = groups
            self._droplet_group_count = int(groups.max())
        if self.solute_com_enabled:
            indices = (_parse_int_list(self.solute_com_atoms_spec)
                       if self.solute_com_atoms_spec else self.qm_atoms)
            selected = self._single_atom_groups(
                self.natom_all, indices, '[solute_com]')
            self._solute_com_selected = (selected > 0).astype(np.int64)

    def _independent_settings_record(self):
        # Keep low-level trajectory/checkpoint helpers usable in focused tests
        # that construct a driver with ``__new__`` instead of running __init__.
        if not hasattr(self, 'droplet_enabled'):
            return {
                'droplet': {'enabled': False},
                'solute_com': {'enabled': False},
                'thermostat': {
                    'ensemble': getattr(self, 'ensemble', 'nve'),
                    'type': getattr(self, 'thermostat', 'off'),
                    'temperature_kelvin': getattr(
                        self, 'thermostat_temperature', 0.0),
                    'friction_ps_inverse': getattr(
                        self, 'thermostat_friction', 0.0),
                },
            }
        return {
            'droplet': {
                'enabled': bool(self.droplet_enabled),
                'center_angstrom': self.droplet_center_angstrom.tolist(),
                'radius_angstrom': self.droplet_radius_angstrom,
                'buffer_angstrom': self.droplet_buffer_angstrom,
                'force_constant_kcal_mol_angstrom2':
                    self.droplet_force_constant_input,
                'target': self.droplet_target,
                'atoms_zero_based': self.droplet_atoms_spec,
                'water_resnames': list(self.droplet_water_resnames),
                'max_penetration_angstrom':
                    self.droplet_max_penetration_angstrom,
                'group_count': int(self._droplet_group_count),
            },
            'solute_com': {
                'enabled': bool(self.solute_com_enabled),
                'center_angstrom': self.solute_com_center_angstrom.tolist(),
                'force_constant_kcal_mol_angstrom2':
                    self.solute_com_force_constant_input,
                'atoms_zero_based': self.solute_com_atoms_spec,
            },
            'thermostat': {
                'ensemble': self.ensemble,
                'type': self.thermostat,
                'temperature_kelvin': self.thermostat_temperature,
                'friction_ps_inverse': self.thermostat_friction,
            },
        }

    def _evaluate_conservative_restraints(self, coordinates, masses):
        coords = np.ascontiguousarray(coordinates, dtype=np.float64).reshape((-1, 3))
        mass = np.ascontiguousarray(masses, dtype=np.float64).reshape(-1)
        if len(coords) != len(mass):
            raise ValueError("restraint coordinate/mass sizes do not match")
        total_force = np.zeros_like(coords)
        self._droplet_energy = 0.0
        self._droplet_max_penetration = 0.0
        self._droplet_active_count = 0
        self._droplet_force = np.zeros_like(coords)
        self._droplet_force_max = 0.0
        if self.droplet_enabled:
            if (self._droplet_group_index is None
                    or len(self._droplet_group_index) != len(coords)):
                raise RuntimeError("droplet target groups were not initialized")
            force = np.zeros_like(coords)
            energy = np.zeros(1, dtype=np.float64)
            penetration = np.zeros(1, dtype=np.float64)
            active = np.zeros(1, dtype=np.int64)
            status = int(oqp.oqp_namd_droplet_boundary(
                len(coords), self._droplet_group_count,
                oqp.ffi.cast("double *", coords.ctypes.data),
                oqp.ffi.cast("double *", mass.ctypes.data),
                oqp.ffi.cast("int64_t *", self._droplet_group_index.ctypes.data),
                oqp.ffi.cast("double *", self.droplet_center.ctypes.data),
                self.droplet_radius, self.droplet_buffer,
                self.droplet_force_constant, self.droplet_max_penetration,
                oqp.ffi.cast("double *", energy.ctypes.data),
                oqp.ffi.cast("double *", force.ctypes.data),
                oqp.ffi.cast("double *", penetration.ctypes.data),
                oqp.ffi.cast("int64_t *", active.ctypes.data),
            ))
            self._droplet_max_penetration = float(penetration[0])
            self._droplet_active_count = int(active[0])
            if status == 1:
                raise RuntimeError(
                    "droplet boundary failsafe: maximum penetration "
                    f"{penetration[0]/ANGSTROM_TO_BOHR:.6f} angstrom exceeds "
                    f"{self.droplet_max_penetration_angstrom:.6f} angstrom")
            if status != 0:
                raise RuntimeError(
                    f"native droplet boundary rejected coordinates (status={status})")
            self._droplet_energy = float(energy[0])
            self._droplet_force = force.copy()
            self._droplet_force_max = float(
                np.max(np.linalg.norm(force, axis=1))) if len(force) else 0.0
            total_force += force

        self._solute_com_energy = 0.0
        self._solute_com_displacement = 0.0
        self._solute_com_force = np.zeros_like(coords)
        if self.solute_com_enabled:
            if (self._solute_com_selected is None
                    or len(self._solute_com_selected) != len(coords)):
                raise RuntimeError("solute COM target atoms were not initialized")
            force = np.zeros_like(coords)
            energy = np.zeros(1, dtype=np.float64)
            displacement = np.zeros(1, dtype=np.float64)
            status = int(oqp.oqp_namd_com_restraint(
                len(coords), oqp.ffi.cast("double *", coords.ctypes.data),
                oqp.ffi.cast("double *", mass.ctypes.data),
                oqp.ffi.cast("int64_t *", self._solute_com_selected.ctypes.data),
                oqp.ffi.cast("double *", self.solute_com_center.ctypes.data),
                self.solute_com_force_constant,
                oqp.ffi.cast("double *", energy.ctypes.data),
                oqp.ffi.cast("double *", force.ctypes.data),
                oqp.ffi.cast("double *", displacement.ctypes.data),
            ))
            if status != 0:
                raise RuntimeError(
                    f"native solute COM restraint rejected coordinates (status={status})")
            self._solute_com_energy = float(energy[0])
            self._solute_com_displacement = float(displacement[0])
            self._solute_com_force = force.copy()
            total_force += force
        self._conservative_restraint_energy = (
            self._droplet_energy + self._solute_com_energy)
        self._conservative_restraint_force = total_force
        return total_force, self._conservative_restraint_energy

    def _add_last_conservative_restraints(self, force, potential_energy):
        if self._conservative_restraint_force is None:
            raise RuntimeError("conservative restraints have not been evaluated")
        return (np.asarray(force) + self._conservative_restraint_force,
                float(potential_energy) + self._conservative_restraint_energy)

    def _langevin_update(self, velocities, masses, istep):
        values = np.ascontiguousarray(velocities, dtype=np.float64)
        mass = np.ascontiguousarray(masses, dtype=np.float64).reshape(-1)
        heat = np.zeros(1, dtype=np.float64)
        friction_au = self.thermostat_friction/(1000.0*FS_TO_AU)
        status = int(oqp.oqp_namd_langevin_thermostat(
            len(mass), self.dt, self.thermostat_temperature, friction_au,
            self.seed, self.rng_stream, int(istep),
            oqp.ffi.cast("double *", mass.ctypes.data),
            oqp.ffi.cast("double *", values.ctypes.data),
            oqp.ffi.cast("double *", heat.ctypes.data),
        ))
        if status != 0:
            raise RuntimeError(
                f"native Langevin thermostat rejected state (status={status})")
        return values, float(heat[0])

    def _apply_thermostat(self, istep):
        self._thermostat_exchange = 0.0
        if self.thermostat == 'off':
            return
        self.vel, self._thermostat_exchange = self._langevin_update(
            self.vel, self.mass, istep)
        self._thermostat_exchange_cumulative += self._thermostat_exchange

    def _md_output_path(self, configured, suffix):
        """Resolve NAMD sidecars beside the main log, not the process CWD."""
        value = str(configured or '').strip()
        log_dir = os.path.dirname(os.path.abspath(self.mol.log))
        if value:
            return value if os.path.isabs(value) else os.path.join(log_dir, value)
        stem = os.path.splitext(os.path.basename(self.mol.log))[0]
        return os.path.join(log_dir, stem + suffix)

    def _restart_manifest_path(self):
        """Return a per-job manifest path that cannot collide in an ensemble."""
        log_dir = os.path.dirname(os.path.abspath(self.mol.log))
        stem = os.path.splitext(os.path.basename(self.mol.log))[0]
        return os.path.join(log_dir, stem + '.namd.restart.oqp')

    def _resolved_velocity_file(self):
        """Resolve a file velocity source exactly as the runtime consumes it."""
        velocity = str(getattr(self, 'velocity_source', '') or '').strip()
        if velocity.lower() in (
                'zero', 'none', '0', 'maxwell', 'boltzmann', 'random'):
            return None
        # Relative velocity paths have historically been interpreted from the
        # process working directory by _init_velocities, not from the input
        # file directory. Keep validation and loading on one resolver.
        return os.path.abspath(os.path.expanduser(velocity))

    def _resolve_qmmm_aux_file(self, name):
        """Resolve a QM/MM auxiliary path exactly as NAMD_QMMM does."""
        value = str(name or '')
        input_file = getattr(self.mol, 'input_file', None)
        input_dir = (os.path.dirname(os.path.abspath(input_file))
                     if input_file else '')
        if (value and input_dir and not os.path.isabs(value)
                and not os.path.exists(value)):
            candidate = os.path.join(input_dir, value)
            if os.path.exists(candidate):
                return candidate
        return value

    def _resolved_basis_definition_file(self, value):
        """Resolve a ``file:`` basis exactly as BasisData.read_basis_fmt."""
        if not isinstance(value, str) or not value.startswith('file:'):
            return None
        filename = value[len('file:'):]
        input_file = getattr(self.mol, 'input_file', None)
        directory = os.path.dirname(input_file) if input_file else ''
        return os.path.join(directory, filename)

    def _validate_sidecar_paths(self):
        """Reject aliases between NAMD sidecars and simulation inputs."""
        outputs = {
            'log_file': self.mol.log,
            'trajectory_file': self.trajectory_file,
            'zpredict_audit_file': self.zpredict_audit_file,
            'restart_file': self.restart_file,
            'restart_manifest_file': self.restart_manifest_file,
        }
        inputs = {}
        for name in ('continuation_checkpoint', 'continuation_trajectory'):
            if getattr(self, name, ''):
                inputs[name] = getattr(self, name)
        original_source = getattr(self.mol, 'oqp_input_source', None)
        resolved_input = getattr(self.mol, 'input_file', None)
        source = original_source or resolved_input
        source_dir = (os.path.dirname(os.path.abspath(source))
                      if source else os.getcwd())
        if original_source:
            inputs['input_source'] = original_source
        if resolved_input:
            inputs['input_file'] = resolved_input

        input_config = getattr(self.mol, 'config', {}).get('input', {})
        for key in ('system', 'system2'):
            geometry = input_config.get(key, '')
            if not isinstance(geometry, str) or not geometry.strip():
                continue
            candidate = geometry.strip()
            if ('\n' in candidate or '\r' in candidate):
                continue
            expanded = os.path.expanduser(candidate)
            path = (expanded if os.path.isabs(expanded)
                    else os.path.join(source_dir, expanded))
            if (os.path.isfile(path)
                    or os.path.splitext(candidate)[1].lower() in ('.xyz', '.pdb')):
                inputs[f'input_{key}'] = path

        basis_inputs = {
            'input_basis': input_config.get('basis', ''),
            'scf_init_basis': getattr(self.mol, 'config', {}).get(
                'scf', {}).get('init_basis', 'none'),
        }
        for name, value in basis_inputs.items():
            path = self._resolved_basis_definition_file(value)
            if path is not None:
                inputs[name] = path

        guess = getattr(self.mol, 'config', {}).get('guess', {})
        for key in ('file', 'file2'):
            value = str(guess.get(key, '') or '').strip()
            if value:
                inputs[f'guess_{key}'] = os.path.abspath(
                    os.path.expanduser(value))

        velocity_file = self._resolved_velocity_file()
        if velocity_file is not None:
            inputs['velocity_file'] = velocity_file

        qmmm = getattr(self.mol, 'config', {}).get('qmmm', {})
        pdb_file = str(qmmm.get('pdb_file', '') or '').strip()
        if pdb_file:
            inputs['qmmm_pdb_file'] = self._resolve_qmmm_aux_file(pdb_file)
        qm_atoms_xyz = str(qmmm.get('qm_atoms_xyz', '') or '').strip()
        if qm_atoms_xyz:
            expanded = os.path.expanduser(qm_atoms_xyz)
            inputs['qmmm_qm_atoms_xyz'] = (
                expanded if os.path.isabs(expanded)
                else os.path.join(source_dir, expanded))
        forcefields = str(qmmm.get('forcefield_files', '')
                          or qmmm.get('forcefield', '') or '')
        for index, item in enumerate(forcefields.replace(',', ' ').split()):
            candidate = self._resolve_qmmm_aux_file(os.path.expanduser(item))
            is_builtin = False
            if not os.path.isfile(candidate):
                try:
                    resource = resources.files('openmm.app').joinpath(
                        'data', *item.split('/'))
                    is_builtin = resource.is_file()
                except (ImportError, ModuleNotFoundError):
                    pass
            # Existing OpenMM resources are package-owned. Any other token is
            # a local input path (including a currently missing input, which a
            # destructive sidecar open must never manufacture accidentally).
            if not is_builtin:
                inputs[f'qmmm_forcefield_file_{index}'] = candidate

        for name, path in self._resolved_tight_binding_artifacts().items():
            inputs[f'tight_binding_{name}'] = path

        resolved = {}
        for name, path in outputs.items():
            canonical = os.path.normcase(os.path.realpath(os.path.abspath(path)))
            if canonical in resolved:
                raise ValueError(
                    f'[md] {name} and {resolved[canonical]} resolve to the '
                    f'same NAMD sidecar path: {path}')
            resolved[canonical] = name
        for name, path in inputs.items():
            canonical = os.path.normcase(os.path.realpath(os.path.abspath(path)))
            if canonical in resolved:
                raise ValueError(
                    f'[md] {resolved[canonical]} aliases simulation input '
                    f'{name}: {path}')
            if os.path.isdir(canonical):
                for output_path, output_name in resolved.items():
                    try:
                        inside = os.path.commonpath(
                            (canonical, output_path)) == canonical
                    except ValueError:
                        inside = False
                    if inside:
                        raise ValueError(
                            f'[md] {output_name} lies inside simulation input '
                            f'directory {name}: {path}')

    def _is_io_rank(self):
        manager = getattr(self.mol, 'mpi_manager', None)
        return manager is None or int(getattr(manager, 'rank', 0)) == 0

    def _run_io_collective(self, operation):
        """Run rank-zero I/O and propagate any failure to every MPI rank."""
        manager = getattr(self.mol, 'mpi_manager', None)
        if manager is None or not bool(getattr(manager, 'use_mpi', False)):
            return operation() if self._is_io_rank() else None

        status = None
        result = None
        if self._is_io_rank():
            try:
                result = operation()
                status = (True, '', '')
            except Exception as error:  # broadcast before raising on rank zero
                status = (False, type(error).__name__, str(error))
        status = manager.bcast(status, root=0)
        if not status[0]:
            exception_type = {
                'FileNotFoundError': FileNotFoundError,
                'OSError': OSError,
                'RuntimeError': RuntimeError,
                'TypeError': TypeError,
                'ValueError': ValueError,
            }.get(status[1], RuntimeError)
            raise exception_type(status[2])
        return result

    def _run_io_collective_result(self, operation):
        """Run rank-zero I/O and broadcast its validated result or failure."""
        manager = getattr(self.mol, 'mpi_manager', None)
        if manager is None or not bool(getattr(manager, 'use_mpi', False)):
            return operation() if self._is_io_rank() else None

        message = None
        if self._is_io_rank():
            try:
                message = (True, operation(), '', '')
            except Exception as error:
                message = (False, None, type(error).__name__, str(error))
        message = manager.bcast(message, root=0)
        if not message[0]:
            exception_type = {
                'FileNotFoundError': FileNotFoundError,
                'OSError': OSError,
                'RuntimeError': RuntimeError,
                'TypeError': TypeError,
                'ValueError': ValueError,
            }.get(message[2], RuntimeError)
            raise exception_type(message[3])
        return message[1]

    def _prepare_md_outputs(self):
        """Start fresh sidecars or preserve them when explicitly restarting."""
        self._run_io_collective(self._prepare_md_outputs_on_io_rank)
        dump_log(
            self.mol,
            title=(f'NAMD files: trajectory={self.trajectory_file} '
                   f'trajectory_interval={self.trajectory_interval}step '
                   f'zpredict_audit={self.zpredict_audit_file} '
                   f'restart={self.restart_file} '
                   f'restart_interval={self.restart_interval}step '
                   f'manifest={self.restart_manifest_file}'),
        )
        if getattr(self, 'odp', None) is not None:
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

    def _prepare_md_outputs_on_io_rank(self):
        if self.restart_requested or getattr(self, 'continuation_checkpoint', ''):
            return
        # Invalidate the runnable stale manifest first.  A failed fresh start
        # must never leave a launchable checkpoint from an older trajectory.
        for path in (self.restart_manifest_file, self.restart_file):
            if os.path.lexists(path):
                os.unlink(path)
        with open(self.trajectory_file, 'w', encoding='utf-8'):
            pass
        with open(self.zpredict_audit_file, 'w', encoding='utf-8'):
            pass
        self._trajectory_prefix_hasher = None
        self._trajectory_prefix_bytes = 0
        self._trajectory_prefix_last_step = None
        self._trajectory_prefix_stat = None

    def _prepare_hop_step(self, istep):
        """Bind the physical MD step and report whether transitions are allowed.

        Electronic coefficients and hop probabilities are propagated at every
        overlap-defined interval, including step 1.  Returning ``False`` only
        suppresses active-state changes, velocity rescaling, and RNG use.
        """
        self._rng_step = int(istep)
        self._last_hop_random = np.nan
        self._last_hop_probabilities = None
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
        if self.restart_requested:
            return np.zeros((self.natom, 3))
        src = self.velocity_source.lower()
        if src in ('zero', 'none', '0'):
            return np.zeros((self.natom, 3))
        if src in ('maxwell', 'boltzmann', 'random'):
            sigma = np.sqrt(KB_HARTREE * self.init_temp / self.mass)  # (natom,)
            v = self._counter_normals((self.natom, 3)) * sigma[:, None]
            return self._remove_com_motion(v)
        # otherwise treat as a file path: "vx vy vz" per atom (atomic units)
        velocity_file = self._resolved_velocity_file()
        if velocity_file is not None and os.path.isfile(velocity_file):
            v = np.loadtxt(velocity_file).reshape((self.natom, 3))
            return self._remove_com_motion(v)
        raise ValueError(f"[md] velocity='{self.velocity_source}' is not zero/maxwell or a readable file")

    def _initial_temperature_metadata(self):
        """Describe the actual initial kinetic temperature stored in the TRJ."""
        requested = getattr(self, 'init_temp', None)
        if hasattr(self, 'm_all') and hasattr(self, 'v_all'):
            masses = np.asarray(self.m_all, dtype=np.float64).reshape(-1)
            velocities = np.asarray(self.v_all, dtype=np.float64).reshape((-1, 3))
            constraints = len(self._ci) if getattr(self, '_has_constraints', False) else 0
            mask = getattr(self, '_move_mask', None)
            moving = len(masses) if mask is None else int(round(float(np.sum(mask))))
            dof = 3*moving - constraints - 3
            source = 'restart' if getattr(self, 'restart_requested', False) else 'maxwell'
        elif hasattr(self, 'mass') and hasattr(self, 'vel'):
            masses = np.asarray(self.mass, dtype=np.float64).reshape(-1)
            velocities = np.asarray(self.vel, dtype=np.float64).reshape((-1, 3))
            dof = 3*len(masses) - 3
            configured = str(getattr(self, 'velocity_source', '')).strip().lower()
            if getattr(self, 'restart_requested', False):
                source = 'restart'
            elif configured in ('zero', 'none', '0'):
                source = 'zero'
            elif configured in ('maxwell', 'boltzmann', 'random'):
                source = 'maxwell'
            else:
                source = 'file'
        else:
            masses = velocities = None
            dof = None
            source = 'unknown'

        measured = None
        if (dof is not None and dof > 0 and masses is not None
                and velocities.shape == (len(masses), 3)):
            kinetic = 0.5*np.sum(masses[:, None]*velocities**2)
            if np.isfinite(kinetic):
                measured = float(2.0*kinetic/(dof*KB_HARTREE))
        return {
            'measured_kelvin': measured,
            'requested_kelvin': requested,
            'dof': dof,
            'velocity_source': source,
        }

    def _trajectory_ensemble_metadata(self):
        """Return honest ensemble provenance for the packed trajectory."""
        if getattr(self, 'ensemble', 'nve') == 'nvt':
            return {
                'ensemble': 'NVT',
                'integrator': 'velocity_verlet_langevin',
                'per_step_velocity_rescaling': False,
                'velocity_rescaling_mode': 'none',
                'thermostat': True,
                'thermostat_type': getattr(self, 'thermostat', 'langevin'),
                'target_temperature_kelvin': getattr(
                    self, 'thermostat_temperature', None),
                'friction_ps_inverse': getattr(
                    self, 'thermostat_friction', None),
                'energy_exchange_field':
                    'thermostat_exchange_cumulative_hartree',
            }
        # Only SOC drivers implement econs and therefore define this attribute;
        # an irrelevant econs spelling on a same-spin deck must remain NVE.
        econs = bool(getattr(self, 'econs', False))
        if econs:
            return {
                'ensemble': 'ENERGY_CONSTRAINED_VELOCITY_RESCALING',
                'integrator': 'velocity_verlet',
                'per_step_velocity_rescaling': True,
                'velocity_rescaling_mode': 'restore_initial_total_energy',
                'thermostat': False,
            }
        return {
            'ensemble': 'NVE',
            'integrator': 'velocity_verlet',
            'per_step_velocity_rescaling': False,
            'velocity_rescaling_mode': 'none',
            'thermostat': False,
        }

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
    def _electronic(self, with_overlap, continuation=False):
        """Run SCF + (optional overlap) + MRSF excitation at the current geometry.

        ``continuation`` requests the previous-orbital guess and the
        reference-following converger without the state overlap (used for
        the intermediate points of an energy-guarded substep).
        """
        mol = self.mol
        self._restart_boundary = False
        cont = with_overlap or continuation
        if self.mo_reuse and cont:
            # Resident orbitals exist once the first geometry has converged;
            # reuse them instead of restarting from the configured guess.
            mol.config['guess']['type'] = 'previous'
        scf_saved = None
        if self.ref_follow != 'off' and cont:
            # Temporarily select the SOMO-preserving continuation converger.
            # The user configuration is restored after the SCF so that the
            # trajectory/restart signature (which echoes the scf section)
            # stays identical to the one written at step 0.
            scf_cfg = mol.config['scf']
            scf_saved = {k: scf_cfg.get(k) for k in ('converger_type', 'escalation', 'vshift')}
            # The SOMO-preserving converger is the primary; if it stalls the
            # SinglePoint ladder escalates to TRAH (warm-started from the same
            # resident orbitals).  A TRAH solution that changes the SOMO
            # configuration is caught by the SOMO check below and handled as
            # a reference switch event rather than aborting the trajectory.
            if self.ref_follow == 'soscf':
                scf_cfg['converger_type'] = 'soscf'
                scf_cfg['escalation'] = 'soscf,trah'
                mol.data.set_scf_converger_type('soscf')
            else:
                scf_cfg['converger_type'] = 'diis'
                scf_cfg['escalation'] = 'soscf,trah'
                if float(scf_cfg.get('vshift', 0.0) or 0.0) <= 0.0:
                    scf_cfg['vshift'] = 0.2
                setter = getattr(mol.data, 'set_scf_vshift', None)
                if setter is not None:
                    setter(float(scf_cfg['vshift']))
                mol.data.set_scf_converger_type('diis')
        try:
            if self.scf_fail == 'restart' and with_overlap:
                sp, ref_energy = self._reference_with_restart()
            else:
                try:
                    sp = SinglePoint(mol)
                    ref_energy = sp.reference()
                except (SCFnotConverged, RuntimeError) as exc:
                    if isinstance(exc, RuntimeError) and 'SCF did not converge' not in str(exc):
                        raise
                    if not (self.mo_reuse and cont and
                            getattr(self, 'scf_guess_retry', True)):
                        raise
                    self._scf_fallback_steps += 1
                    dump_log(mol, title='PyOQP: NAMD SCF continuation from the '
                             'previous-step orbitals failed (%s); re-solving the '
                             'reference from a fresh Huckel guess (attempt %d)'
                             % (exc, self._scf_fallback_steps), section='input')
                    sp, ref_energy = self._reference_fresh_guess()

        finally:
            if scf_saved is not None:
                scf_cfg = mol.config['scf']
                for k, v in scf_saved.items():
                    if v is None:
                        scf_cfg.pop(k, None)
                    else:
                        scf_cfg[k] = v
        self._require_converged_reference()
        if with_overlap:
            mol.back_door = (self.prev_xyz, self.prev_data)
            BasisOverlap(mol).overlap()
        sp.excitation(ref_energy)
        LastStep(mol).compute(mol)

    def _energy_retry_state(self):
        """Preserve electronic histories and diagnostics before a trial step."""
        prefixes = ('_ba_', '_last_', '_pending_', '_nacme_gate_',
                    '_analytic_tdc_', '_nacme_reference_')
        names = {'_somo_switch_step', '_window_leak_step', '_somo_switch_count',
                 '_window_leak_count', '_overlap_collapse_steps', '_restart_boundary'}
        return {name: copy.deepcopy(value) for name, value in self.__dict__.items()
                if name.startswith(prefixes) or name in names}

    def _restore_energy_retry_state(self, saved):
        for name in self._energy_retry_state():
            if name not in saved:
                delattr(self, name)
        for name, value in saved.items():
            setattr(self, name, copy.deepcopy(value))

    @staticmethod
    def _energy_refinement_counts(max_substeps):
        """Increasing subdivisions, with the configured maximum included once."""
        maximum = int(max_substeps)
        counts = []
        count = 2
        while count < maximum:
            counts.append(count)
            count *= 2
        if maximum >= 2:
            counts.append(maximum)
        return counts

    def _reference_fresh_guess(self):
        """Retry once at the same geometry and require actual SCF convergence."""
        mol = self.mol
        scf_cfg = mol.config['scf']
        guess_cfg = mol.config['guess']
        saved_scf = {k: scf_cfg.get(k) for k in ('converger_type', 'escalation')}
        saved_guess = guess_cfg.get('type')
        try:
            guess_cfg['type'] = 'huckel'
            scf_cfg['converger_type'] = 'diis'
            scf_cfg['escalation'] = 'soscf,trah'
            mol.data.set_scf_converger_type('diis')
            sp = SinglePoint(mol)
            energy = sp.reference()
            self._require_converged_reference()
            dump_log(mol, title=('NAMD SCF recovery succeeded: fresh Huckel guess; '
                                 'reference energy %.12f Hartree; requested SCF '
                                 'criterion satisfied; orbital/state continuity '
                                 'is evaluated separately' % float(energy[0] if isinstance(energy, (list, tuple)) else energy)), section='input')
            return sp, energy
        finally:
            if saved_guess is None:
                guess_cfg.pop('type', None)
            else:
                guess_cfg['type'] = saved_guess
            for key, value in saved_scf.items():
                if value is None:
                    scf_cfg.pop(key, None)
                else:
                    scf_cfg[key] = value
            mol.data.set_scf_converger_type(saved_scf['converger_type'] or 'diis')

    def _reference_with_restart(self):
        """Primary converger only; on failure perform a GAMESS-style restart.

        The SinglePoint escalation ladder is disabled by naming the primary
        converger as the whole chain.  If the primary converger does not
        converge from the resident (previous-step) orbitals, the reference is
        re-solved from a fresh Huckel guess with SOSCF, exactly as the archived
        KNU-GAMESS restart inputs do (``diis=.f. soscf=.t.``).  The step is
        then marked as a restart boundary.
        """
        mol = self.mol
        scf_cfg = mol.config['scf']
        guess_cfg = mol.config['guess']
        saved = {
            'escalation': scf_cfg.get('escalation', ''),
            'converger_type': scf_cfg.get('converger_type', 'diis'),
            'guess_type': guess_cfg.get('type', 'huckel'),
        }
        primary = str(saved['converger_type'] or 'diis')
        scf_cfg['escalation'] = primary      # chain minus primary == empty
        try:
            sp = SinglePoint(mol)
            return sp, sp.reference()
        except SCFnotConverged:
            pass
        except RuntimeError as err:
            if 'SCF did not converge' not in str(err):
                raise
        finally:
            scf_cfg['escalation'] = saved['escalation']
        self._scf_restart_steps += 1
        dump_log(
            mol,
            title=('NAMD: %s did not converge from the previous-step orbitals; '
                   'restarting the reference from a Huckel guess with SOSCF '
                   '(restart boundary %d)' % (primary, self._scf_restart_steps)),
            section='input')
        guess_cfg['type'] = 'huckel'
        scf_cfg['converger_type'] = 'soscf'
        scf_cfg['escalation'] = 'soscf'       # again no further escalation
        try:
            sp = SinglePoint(mol)
            ref_energy = sp.reference()
        finally:
            guess_cfg['type'] = saved['guess_type']
            scf_cfg['converger_type'] = saved['converger_type']
            scf_cfg['escalation'] = saved['escalation']
            mol.data.set_scf_converger_type(saved['converger_type'])
        self._restart_boundary = True
        return sp, ref_energy

    def _require_converged_reference(self):
        """Never propagate an electronic reference rejected by the SCF solver."""
        # Tight-binding adapters validate their own SCC solver and do not
        # populate the ab initio SCF status field.
        if is_tb_method(self.mol.config['input']['method']):
            return
        if not self.mol.mol_energy.SCF_converged:
            raise RuntimeError(
                'NAMD cannot continue: SCF did not converge; no force or '
                'electronic propagation is permitted for this reference.')

    def _active_gradient(self):
        """Compute and return the gradient (natom,3) on the current active state."""
        self._require_converged_reference()
        mol = self.mol
        mol.config['properties']['grad'] = [self.active]
        Gradient(mol).gradient()
        gradient = np.array(mol.grads[self.active]).reshape((self.natom, 3))
        odp = self._evaluate_odp(mol.get_system().reshape((self.natom, 3)))
        if odp is not None:
            gradient = gradient - odp['force']
        return gradient

    def _state_overlap(self, istep=None, *, update_analytic=True):
        """Compute the phase-corrected state overlap S(i,j)=<i(t-dt)|j(t)>."""
        NACME(self.mol).nacme()
        state_overlap = canonical_state_overlap(
            self.mol.data["OQP::td_states_overlap"]
        )
        self._last_state_overlap = np.array(state_overlap, copy=True)
        # Diagnose a collapsed retained-manifold overlap (every column of the
        # old->new state overlap nearly zero).  A collapsed matrix cannot
        # describe continuous states; NPI in particular then returns a large
        # spurious coupling.  Report it so the trajectory can be audited.
        column_norm = np.linalg.norm(np.asarray(state_overlap, dtype=float), axis=0)
        if np.all(np.isfinite(column_norm)) and column_norm.max() < 0.5:
            self._overlap_collapse_steps += 1
            dump_log(
                self.mol,
                title=('NAMD WARNING: collapsed state overlap at step %s '
                       '(max column norm %.3f, %d collapsed steps so far)'
                       % (istep, column_norm.max(), self._overlap_collapse_steps)),
                section='nacm', info=state_overlap)
        self._last_overlap_tdc = np.array(self._compute_tdc(state_overlap), copy=True)
        # SOMO identity check: the aligned MO overlap of the two singly
        # occupied orbitals with their previous-step counterparts.  A value
        # below somo_tol means the ROHF reference changed its open-shell
        # configuration (reference switch event).
        self._somo_switch_step = False
        self._window_leak_step = False
        try:
            mo_ov = np.abs(np.asarray(
                self.mol.data['OQP::mo_tracking_overlap'], dtype=float).ravel())
            nocc = int(self.mol.data['nelec_A'])
            somo = mo_ov[nocc - 2:nocc]
        except Exception:
            somo = None
        if somo is not None and somo.size == 2 and np.all(np.isfinite(somo)) \
                and somo.min() < self.somo_tol:
            self._somo_switch_step = True
            self._somo_switch_count += 1
            dump_log(
                self.mol,
                title=('NAMD WARNING: SOMO identity change at step %s '
                       '(SOMO overlaps with the previous step %.3f %.3f < %.2f; '
                       'reference switch event %d)'
                       % (istep, somo[0], somo[1], self.somo_tol,
                          self._somo_switch_count)),
                section='input')
        # Retained-window leakage of the active state
        try:
            col = np.asarray(state_overlap, dtype=float)[:, self.active - 1]
            leak = float(np.linalg.norm(col))
        except Exception:
            leak = np.nan
        if np.isfinite(leak) and leak < 0.7:
            self._window_leak_count += 1
            self._window_leak_step = True
            dump_log(
                self.mol,
                title=('NAMD WARNING: active-state overlap column norm %.3f < 0.7 '
                       'at step %s; %.0f%% of the state lies outside the retained '
                       'window (consider a larger nstate); event %d'
                       % (leak, istep, 100.0*(1.0 - leak**2), self._window_leak_count)),
                section='input')
        self._update_baeck_an_check(istep, state_overlap)
        if update_analytic and self._needs_analytic_nac():
            self._update_analytic_nac(istep, compare_overlap=True)
        return state_overlap

    def _validated_td_energies(self, tag):
        """Return an exact finite nstate vector safe for native pointer use."""
        energies = np.asarray(self.mol.data[tag], dtype=np.float64).reshape(-1)
        if (energies.shape != (self.nstate,)
                or not np.all(np.isfinite(energies))):
            raise RuntimeError(
                f'{tag} must be an exact finite nstate TD-energy vector')
        return np.ascontiguousarray(energies)

    def _scale_analytic_velocity_contractions(self, factor, istep):
        """Keep cached d.v couplings consistent with a uniform velocity scale.

        The analytic TDC endpoint of this step was contracted with the velocity
        before a numerical energy correction.  Both the propagated TDC and the
        endpoint history are linear in v, so they scale by the same factor.
        """
        if (getattr(self, '_last_analytic_tdc', None) is None
                or getattr(self, '_last_analytic_step', None) != int(istep)):
            return
        self._last_analytic_tdc = np.asarray(self._last_analytic_tdc) * factor
        if (getattr(self, '_last_analytic_pair', None) is None
                and getattr(self, '_analytic_tdc_previous', None) is not None):
            self._analytic_tdc_previous = (
                np.asarray(self._analytic_tdc_previous) * factor)

    def _log_rescale_resolution(self):
        """Record which velocity-rescaling direction rescale=auto selected."""
        requested = str(self.mol.config.get('md', {}).get('rescale', 'auto')).strip().lower()
        if requested != 'auto':
            return
        reason = getattr(self, '_rescale_auto_issue', None)
        dump_log(self.mol, title=(
            'NAMD velocity rescaling: rescale=auto -> %s%s' % (
                self.rescale_provider,
                '' if reason is None else ' (analytic NAC unavailable: %s)' % reason)),
            section='input')

    def _needs_analytic_nac(self):
        """Return whether this trajectory consumes the resident analytic NAC."""
        return (
            self.tdc_provider == 'analytic'
            or self.rescale_provider == 'analytic_nac'
            or self.nacme_check == 'analytic'
        )

    def _gradient_nac_fusion_enabled(self):
        return (
            self._needs_analytic_nac()
            and os.environ.get('OQP_MRSF_NAC_ZV_FUSE_GRADIENT', '')
            .strip().lower() in ('1', 'y', 'yes', 't', 'true', 'on')
        )

    def _update_analytic_nac(self, istep=None, *, compare_overlap=False,
                             pair=None):
        """Evaluate phase-aligned analytic d and contract it with velocity.

        The endpoint value is used when ``tdc=analytic``.  A trapezoidal value
        is retained separately for comparison with the overlap integrated over
        the preceding nuclear interval.

        ``pair=(I, J)`` (1-based) evaluates only that physical pair; the
        reference mask then marks every other pair as not evaluated and no
        trapezoidal history is kept, because the selected vector serves one
        hop candidate only.
        """
        from oqp.library.nac_analytic import analytic_nac, _resident_pair_cartesian

        if pair is not None and compare_overlap:
            raise ValueError('an overlap NAC check requires every analytic pair')
        if getattr(self.mol, '_nac_fused_gradient_ready', False):
            _nacv = _resident_pair_cartesian(
                self.mol.data['OQP::nac_nacv'], self.nstate, self.natom)
            dcv = _resident_pair_cartesian(
                self.mol.data['OQP::nac_dcv'], self.nstate, self.natom)
            self.mol._nac_fused_gradient_ready = False
            pair = None   # the fused solve already evaluated every pair
        elif pair is None:
            _nacv, dcv = analytic_nac(self.mol)
        else:
            _nacv, dcv = analytic_nac(self.mol, pair=pair)
        dcv = np.asarray(dcv, dtype=np.float64).reshape(
            (self.nstate, self.nstate, self.natom, 3))
        try:
            predictor_dcv = _resident_pair_cartesian(
                self.mol.data['OQP::nac_predictor_dcv'],
                self.nstate, self.natom)
            predictor_nacv = _resident_pair_cartesian(
                self.mol.data['OQP::nac_predictor_nacv'],
                self.nstate, self.natom)
        except (AttributeError, KeyError, RuntimeError, TypeError, ValueError):
            predictor_dcv = None
            predictor_nacv = None
        if predictor_dcv is not None:
            self._write_zpredict_audit_row(
                istep, dcv, predictor_dcv,
                np.asarray(_nacv, dtype=np.float64).reshape(dcv.shape),
                predictor_nacv, pair=pair)
        endpoint = np.einsum('ijac,ac->ij', dcv, self.vel, optimize=True)
        self._last_analytic_dcv = np.array(dcv, copy=True)
        self._last_analytic_tdc = np.array(endpoint, copy=True)
        self._last_analytic_pair = None if pair is None else (
            int(pair[0]), int(pair[1]))
        self._last_analytic_step = None if istep is None else int(istep)

        if pair is None:
            previous = self._analytic_tdc_previous
            centered = None if previous is None else 0.5*(previous + endpoint)
            self._analytic_tdc_centered = (
                None if centered is None else np.array(centered, copy=True))
            self._analytic_tdc_previous = np.array(endpoint, copy=True)
        else:
            # A single-pair vector is not a step-to-step series; keep no
            # trapezoidal history so a later full evaluation cannot center
            # against an incomplete matrix.
            centered = None
            self._analytic_tdc_centered = None
            self._analytic_tdc_previous = None

        # Preserve the analytic quantity in the dense trajectory even when it
        # is the production provider rather than a validation reference.
        reference = endpoint if centered is None else centered
        if pair is None:
            mask = np.ones((self.nstate, self.nstate), dtype=np.int32)
            np.fill_diagonal(mask, 0)
        else:
            # Only the selected pair was evaluated; every other zero entry is
            # "not evaluated", never a computed zero.
            mask = np.zeros((self.nstate, self.nstate), dtype=np.int32)
            i, j = pair[0] - 1, pair[1] - 1
            mask[i, j] = mask[j, i] = 1
        self._nacme_reference_tdc = np.array(reference, copy=True)
        self._nacme_reference_mask = mask
        self._nacme_reference_source = 2

        if (compare_overlap and self.nacme_check == 'analytic'
                and centered is not None):
            gate = self._run_nacme_gate(
                self._last_overlap_tdc,
                centered,
                reference_mask=mask,
                source='analytic',
                center_step=None if istep is None else int(istep),
                signed=True,
            )
            dump_log(
                self.mol,
                title='NACME check: centered analytic d_ij dot velocity',
                section='nacm',
                info=centered,
            )
            return gate
        return None

    def _write_zpredict_audit_row(self, istep, exact_dcv, predictor_dcv,
                                  exact_nacv, predictor_nacv, *, pair=None):
        """Append gauge-aligned full-vector and velocity-contraction errors."""
        if not self._is_io_rank():
            self._io_barrier()
            return
        n = self.nstate
        pairs = np.triu_indices(n, 1)
        if pair is not None:
            pairs = (np.array([min(pair) - 1]), np.array([max(pair) - 1]))
        exact_d = np.asarray(exact_dcv, dtype=float)[pairs].reshape(len(pairs[0]), -1)
        pred_d = np.asarray(predictor_dcv, dtype=float)[pairs].reshape(len(pairs[0]), -1)
        exact_h = np.asarray(exact_nacv, dtype=float)[pairs].reshape(len(pairs[0]), -1)
        pred_h = np.asarray(predictor_nacv, dtype=float)[pairs].reshape(len(pairs[0]), -1)
        delta_d = pred_d - exact_d
        delta_h = pred_h - exact_h
        dot = np.sum(pred_d*exact_d, axis=1)
        denom = np.linalg.norm(pred_d, axis=1)*np.linalg.norm(exact_d, axis=1)
        cosine = np.divide(dot, denom, out=np.ones_like(dot), where=denom > 0.0)
        exact_vd = np.einsum('pac,ac->p', exact_d.reshape(-1, self.natom, 3),
                             self.vel, optimize=True)
        pred_vd = np.einsum('pac,ac->p', pred_d.reshape(-1, self.natom, 3),
                            self.vel, optimize=True)
        delta_vd = pred_vd - exact_vd
        mode = os.environ.get('OQP_MRSF_NAC_ZV_PREDICTOR', 'off').strip()
        result = {
            'step': '' if istep is None else int(istep),
            'mode': mode,
            'production_is_predictor': int(mode.lower().endswith('_approx')),
            'eta': os.environ.get('OQP_MRSF_NAC_ZV_ETA', '1.0'),
            'exact_every': os.environ.get('OQP_MRSF_NAC_ZV_EXACT_EVERY', '0'),
            'pair_count': len(pairs[0]),
            'd_rms': float(np.sqrt(np.mean(delta_d*delta_d))),
            'd_max': float(np.max(np.abs(delta_d))),
            'd_relative_l2': float(np.linalg.norm(delta_d) /
                                   max(np.linalg.norm(exact_d), 1.0e-300)),
            'd_cosine_mean': float(np.mean(cosine)),
            'd_cosine_min': float(np.min(cosine)),
            'h_rms': float(np.sqrt(np.mean(delta_h*delta_h))),
            'h_max': float(np.max(np.abs(delta_h))),
            'h_relative_l2': float(np.linalg.norm(delta_h) /
                                   max(np.linalg.norm(exact_h), 1.0e-300)),
            'vd_rms': float(np.sqrt(np.mean(delta_vd*delta_vd))),
            'vd_max': float(np.max(np.abs(delta_vd))),
            'vd_relative_l2': float(np.linalg.norm(delta_vd) /
                                    max(np.linalg.norm(exact_vd), 1.0e-300)),
        }
        tracking = self.mol.get_state_tracking()
        result['tracking_overlap_min'] = (
            '' if tracking is None else
            float(np.min(np.abs(np.asarray(tracking['matched_overlap'], dtype=float))))
        )
        result['tracking_margin_min'] = (
            '' if tracking is None else
            float(np.min(np.asarray(tracking['margin'], dtype=float)))
        )
        columns = tuple(result)
        path = self.zpredict_audit_file
        needs_header = not os.path.exists(path) or os.path.getsize(path) == 0
        with open(path, 'a', encoding='utf-8') as stream:
            if needs_header:
                stream.write('\t'.join(columns) + '\n')
            stream.write('\t'.join(str(result[name]) for name in columns) + '\n')
            stream.flush()
            os.fsync(stream.fileno())
        self._io_barrier()

    def _update_baeck_an_check(self, istep, state_overlap):
        """Compare overlap TDC magnitudes with a centred TD-Baeck-An estimate.

        TD-BA is phase-free and therefore cannot validate the signed gauge.
        It is retained only as an independent energy-curvature diagnostic.
        """
        if (self.nacme_check != 'baeck_an'
                and getattr(self, 'tdc_provider', '') != 'baeck_an'):
            return

        n = self.nstate
        data = self.mol.data
        energies_old = self._validated_td_energies("OQP::td_energies_old")
        energies_current = self._validated_td_energies("OQP::td_energies")
        tdc_current = np.ascontiguousarray(
            self._compute_tdc(state_overlap), dtype=np.float64
        )
        dt_right = float(self.dt)

        if self._ba_energy_center is None:
            self._reset_nacme_gate_evaluation()
            self._last_baeck_an_tdc = None
            if getattr(self, 'tdc_provider', '') == 'baeck_an':
                self._last_tdc_source = 1
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
            self._reset_nacme_gate_evaluation()
            self._last_baeck_an_tdc = None
            if getattr(self, 'tdc_provider', '') == 'baeck_an':
                self._last_tdc_source = 1
            self._ba_energy_left = energies_old.copy()
            self._ba_energy_center = energies_current.copy()
            self._ba_tdc_left = tdc_current.copy()
            self._ba_dt_left = dt_right
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
        signed_ba_tdc = self._signed_baeck_an_tdc(
            ba_tdc, overlap_center)
        self._last_baeck_an_tdc = signed_ba_tdc
        if getattr(self, 'tdc_provider', '') == 'baeck_an':
            self._last_tdc_source = 3
        center_step = None if istep is None else int(istep) - 1
        if self.nacme_check == 'baeck_an':
            ba_mask = np.asarray(np.abs(ba_tdc) > 0.0, dtype=np.int32)
            gate = self._run_nacme_gate(
                overlap_center,
                ba_tdc,
                reference_mask=ba_mask,
                source='TD-Baeck-An',
                center_step=center_step,
                evaluation_step=istep,
                signed=False,
            )
            self._ba_last = {
                'center_step': center_step,
                'baeck_an_tdc': ba_tdc.copy(),
                'signed_baeck_an_tdc': signed_ba_tdc.copy(),
                'overlap_tdc_centered': overlap_center.copy(),
                'signed_pair_count': int(np.count_nonzero(
                    np.triu(signed_ba_tdc, k=1))),
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
        else:
            # Production Baeck-An dynamics needs the coupling but not the
            # optional matrix dump and comparison at every nuclear step.
            self._ba_last = None

        self._ba_energy_left = self._ba_energy_center.copy()
        self._ba_energy_center = energies_current.copy()
        self._ba_tdc_left = tdc_current.copy()
        self._ba_dt_left = dt_right

    @staticmethod
    def _signed_baeck_an_tdc(baeck_an_tdc, overlap_tdc):
        """Apply only a transported wavefunction-gauge sign to TD-BA.

        Baeck-An supplies a magnitude.  The phase-tracked overlap coupling
        supplies the sign; a pair with an exactly indeterminate sign remains
        zero.  Constructing one triangle and reflecting it makes
        antisymmetry exact rather than a floating-point postcondition.
        """
        magnitude_matrix = np.asarray(baeck_an_tdc, dtype=np.float64)
        phase_matrix = np.asarray(overlap_tdc, dtype=np.float64)
        if (magnitude_matrix.ndim != 2
                or magnitude_matrix.shape[0] != magnitude_matrix.shape[1]
                or phase_matrix.shape != magnitude_matrix.shape
                or not np.all(np.isfinite(magnitude_matrix))
                or not np.all(np.isfinite(phase_matrix))):
            raise ValueError(
                'Baeck-An magnitude and overlap sign matrices must be finite '
                'square matrices of the same shape')
        signed = np.zeros_like(magnitude_matrix)
        for i in range(magnitude_matrix.shape[0]):
            for j in range(i + 1, magnitude_matrix.shape[1]):
                magnitude = abs(float(magnitude_matrix[i, j]))
                phase_reference = float(phase_matrix[i, j])
                if magnitude > 0.0 and phase_reference != 0.0:
                    value = np.copysign(magnitude, phase_reference)
                    signed[i, j] = value
                    signed[j, i] = -value
        return signed

    def _reset_nacme_gate_evaluation(self):
        """Clear streak and record state for a non-evaluable NACME interval."""
        self._nacme_gate_failures = 0
        self._nacme_gate_last = None
        self._nacme_candidate_tdc = None
        self._nacme_reference_tdc = None
        self._nacme_reference_mask = None
        self._nacme_reference_source = 0
        self._pending_nacme_gate_error = None
        self._ba_last = None

    def _run_nacme_gate(self, candidate_tdc, reference_tdc, *,
                        reference_mask=None, source='reference',
                        center_step=None, evaluation_step=None, signed=False):
        """Run the common resident-Fortran NACME validation gate.

        The analytic NAC path contracts the phase-aligned coupling vector with
        the nuclear velocity and calls this method with ``signed=True``.
        TD-Baeck-An calls it with ``signed=False`` because an energy-only
        estimate has no wavefunction gauge. Thus the invariant and policy
        machinery is shared without treating the approximate TD-BA sign as
        physical.
        """
        self._pending_nacme_gate_error = None
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
        compared_pairs = int(counts[0])
        invariant_failures = int(counts[1])
        reference_failures = int(counts[2])
        native_error = None
        if status != 0:
            # The native kernel uses -2/-3 for non-finite candidate/reference
            # matrices. Preserve that fatal status, but populate the dense
            # diagnostic state before enforcing it after trajectory output.
            native_error = RuntimeError(
                f"native NACME validation gate failed for {source} "
                f"(status={status})"
            )
            metrics[:] = np.nan
            invariant_failures = max(1, invariant_failures)
            verdict = 'fail'
        elif invariant_failures or reference_failures:
            verdict = 'fail'
        elif compared_pairs == 0:
            verdict = 'not_evaluable'
        else:
            verdict = 'pass'
        if status == 0 and reference_failures:
            self._nacme_gate_failures += 1
        elif status == 0:
            self._nacme_gate_failures = 0

        result = {
            'source': source,
            'center_step': center_step,
            'evaluation_step': evaluation_step,
            'signed_comparison': bool(signed),
            'native_status': int(status),
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
        if self.nacme_gate == 'off':
            result['verdict'] = 'off'
            self._nacme_gate_failures = 0
            self._nacme_gate_last = None
            self._nacme_candidate_tdc = None
            self._nacme_reference_tdc = None
            self._nacme_reference_mask = None
            self._nacme_reference_source = 0
            self._pending_nacme_gate_error = None
            return result
        self._nacme_gate_last = result
        self._nacme_candidate_tdc = candidate.copy()
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
        error = native_error
        if error is None and self.nacme_gate == 'error':
            if invariant_failures:
                error = RuntimeError(
                    f"NACME invariant gate failed for {source} at step {center_step}"
                )
            elif self._nacme_gate_failures >= self.nacme_gate_consecutive:
                error = RuntimeError(
                    f"NACME reference gate failed {self._nacme_gate_failures} "
                    f"consecutive times for {source} at step {center_step}"
                )
        if error is not None and self._pending_nacme_gate_error is None:
            self._pending_nacme_gate_error = error
        return result

    def _enforce_nacme_gate(self):
        """Stop only after the failing NACME point reaches the dense TRJ."""
        error = getattr(self, '_pending_nacme_gate_error', None)
        self._pending_nacme_gate_error = None
        if error is not None:
            raise error

    def _update_nve_gate(self, istep, epot, ekin, transition_energy_jump=np.nan):
        """Audit microcanonical energy conservation for FSSH/ISC dynamics."""
        self._pending_nve_gate_error = None
        total = float(epot + ekin)
        if getattr(self, 'ensemble', 'nve') != 'nve':
            self._nve_gate_last = {
                'step': int(istep), 'verdict': 'off',
                'total_energy': total, 'drift': np.nan,
                'step_change': np.nan,
                'transition_energy_jump': float(transition_energy_jump),
                'drift_rate': np.nan, 'drift_failure': False,
                'step_failure': False, 'transition_failure': False,
                'consecutive_failures': 0,
            }
            return self._nve_gate_last
        if not np.isfinite(total):
            drift = np.nan
            step_change = np.nan
        elif self._nve_reference_energy is None:
            self._nve_reference_energy = total
            self._nve_previous_energy = total
            drift = 0.0
            step_change = 0.0
        else:
            # A numerical velocity correction forces the total energy back to
            # its previous value.  Audit the energy the dynamics actually
            # produced: add back this step's absorbed change and the
            # cumulative absorbed energy, so the gate cannot be satisfied by
            # the correction it is meant to police.
            step_correction = float(getattr(self, '_step_numerical_correction', 0.0))
            absorbed = float(getattr(self, '_disc_energy_absorbed', 0.0))
            drift = total + absorbed - self._nve_reference_energy
            step_change = total + step_correction - self._nve_previous_energy
        transition_jump = float(transition_energy_jump)
        time_fs = self._physical_time_fs(istep)
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

    def _enforce_nacme_gate(self):
        """Stop only after the failing NACME matrices are in the dense TRJ."""
        error = self._pending_nacme_gate_error
        self._pending_nacme_gate_error = None
        if error is not None:
            raise error

    def _write_md_trajectory(self, istep, coordinates, epot, ekin, hopped):
        """Append one lossless, fixed-width record to the dense binary TRJ."""
        last_nve = getattr(self, '_nve_gate_last', None)
        last_nacme = getattr(self, '_nacme_gate_last', None)
        gate_failure = (
            getattr(self, '_pending_nve_gate_error', None) is not None
            or getattr(self, '_pending_nacme_gate_error', None) is not None
            or (isinstance(last_nve, dict)
                and last_nve.get('verdict') == 'fail')
            or (isinstance(last_nacme, dict)
                and last_nacme.get('verdict') == 'fail'))
        if (istep % self.trajectory_interval != 0
                and istep != self.nstep and not gate_failure):
            return
        return self._run_io_collective(
            lambda: self._write_md_trajectory_on_io_rank(
                istep, coordinates, epot, ekin, hopped))

    def _write_md_trajectory_on_io_rank(self, istep, coordinates, epot, ekin,
                                        hopped):
        """Append one packed trajectory record on rank zero."""
        coords = np.asarray(coordinates, dtype=np.float64).reshape((-1, 3))
        if hasattr(self, 'r_all') and len(coords) == len(self.r_all):
            velocities = np.asarray(self.v_all, dtype=np.float64).reshape(coords.shape)
        else:
            velocities = np.asarray(self.vel, dtype=np.float64).reshape(coords.shape)
        ncv = self.odp.ncv if getattr(self, 'odp', None) is not None else 0
        trajectory_nstate = int(np.asarray(self.coef).size)
        dtype = _namd_trajectory_dtype(trajectory_nstate, len(coords), ncv)
        new_file = (not os.path.exists(self.trajectory_file)
                    or os.path.getsize(self.trajectory_file) == 0)
        if new_file:
            temperature = self._initial_temperature_metadata()
            ensemble = self._trajectory_ensemble_metadata()
            header = {
                'schema_version': NAMD_TRAJECTORY_SCHEMA_VERSION,
                'nstate': trajectory_nstate,
                'natom': len(coords),
                'representation': getattr(
                    self, '_trajectory_representation',
                    'same_spin_adiabatic'),
                'ncv': ncv,
                'record_bytes': dtype.itemsize,
                'signature': self._restart_signature(),
                'time_origin_fs': getattr(self, '_time_origin_fs', 0.0),
                'continuation': getattr(self, '_continuation_provenance', None),
                'wham_system_identity': getattr(
                    self, '_wham_system_identity', {'kind': 'unavailable'}),
                'electronic_representation': getattr(
                    self, '_trajectory_representation', 'same_spin_adiabatic'),
                'ensemble': ensemble['ensemble'],
                'ensemble_provenance': ensemble,
                'initial_temperature_kelvin': temperature['measured_kelvin'],
                'initial_temperature_degrees_of_freedom': temperature['dof'],
                'requested_initial_temperature_kelvin': temperature[
                    'requested_kelvin'],
                'initial_velocity_source': temperature['velocity_source'],
                'odp': self._odp_provenance(),
                'wham': {
                    'reaction_coordinate_field': 'odp_xi',
                    'window_field': 'odp_window',
                    'bias_field': 'odp_bias_hartree',
                    'unbiased_potential_field': 'e_unbiased_pot_hartree',
                    'energy_unit': 'hartree',
                    'temperature_note': (
                        'initial_temperature_kelvin is measured from the '
                        'initial kinetic energy and stated degrees of freedom; '
                        'requested_initial_temperature_kelvin is only a '
                        'velocity-generation target, not an NVT thermostat '
                        'temperature'),
                },
                'units': {
                    'time': 'fs', 'coordinates': 'bohr',
                    'velocities': 'bohr/atomic_time', 'energies': 'hartree',
                    'tdc': 'atomic_time^-1', 'penetration': 'bohr',
                    'restraint_force': 'hartree/bohr',
                    'thermostat_exchange': 'hartree (positive into system)',
                },
                'independent_controls': self._independent_settings_record(),
                'reference_source': {'0': 'none', '1': 'TD-Baeck-An',
                                     '2': 'analytic', '127': 'other'},
                'tdc_source': {'0': 'overlap_fd', '1': 'overlap_npi',
                               '2': 'analytic_endpoint',
                               '3': 'lagged_baeck_an_overlap_sign'},
                'rescale_source': {
                    '0': 'isotropic', '1': 'analytic_nac',
                    '2': 'hop_triggered_analytic_nac',
                },
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
            header_record = (NAMD_TRAJECTORY_MAGIC
                             + struct.pack('<Q', len(encoded)) + encoded)
            with open(self.trajectory_file, 'wb') as stream:
                stream.write(header_record)
                stream.flush()
                os.fsync(stream.fileno())
            self._trajectory_prefix_hasher = hashlib.sha256(header_record)
            self._trajectory_prefix_bytes = len(header_record)
            self._trajectory_prefix_last_step = None
            self._trajectory_prefix_stat = self._trajectory_stat_identity()
        else:
            cached_hasher = getattr(self, '_trajectory_prefix_hasher', None)
            if cached_hasher is None:
                scanned = self._scan_trajectory_prefix(INT64_MAX)
                if scanned['partial_bytes'] or scanned['removed_records']:
                    raise ValueError(
                        'NAMD trajectory has an incomplete append history')
                self._remember_trajectory_prefix(scanned)
            else:
                self._require_unchanged_trajectory_prefix()
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
                'droplet_energy_hartree',
                'droplet_max_penetration_bohr',
                'solute_com_energy_hartree',
                'solute_com_displacement_bohr',
                'conservative_restraint_energy_hartree',
                'thermostat_exchange_hartree',
                'thermostat_exchange_cumulative_hartree',
                'thermostat_adjusted_energy_hartree',
                'state_energies', 'populations', 'coef_real', 'coef_imag',
                'coordinates_bohr', 'velocities_au', 'state_overlap',
                'state_overlap_imag', 'overlap_tdc_au',
                'overlap_tdc_imag_au', 'gate_candidate_tdc_au',
                'reference_tdc_au', 'gate_metrics',
                'nve_metrics',
                'odp_xi', 'odp_cv_raw', 'odp_cv_scaled',
                'odp_cv_perpendicular', 'odp_perpendicular_norm',
                'odp_bias_parallel_hartree',
                'odp_bias_perpendicular_hartree', 'odp_bias_hartree',
                'tracking_phase', 'tracking_phase_initial',
                'tracking_previous_phase_initial', 'tracking_overlap',
                'tracking_margin', 'rescale_gamma',
                'rescale_discriminant'):
            record[field] = np.nan
        record['tracking_order'] = -1
        record['tracking_raw_order'] = -1
        record['tracking_lineage'] = -1
        record['gate_center_step'] = -1
        record['gate_verdict'] = -1
        record['gate_streak'] = -1
        record['nve_verdict'] = -1
        record['odp_window'] = -1
        gate = self._nacme_gate_last or {}
        time_fs = self._physical_time_fs(istep)
        record['step'] = istep
        record['time_fs'] = time_fs
        record['active'] = self.active
        record['hopped'] = int(bool(hopped))
        record['tdc_source'] = getattr(
            self, '_last_tdc_source', getattr(self, 'tdc_scheme', 0))
        default_rescale = {
            'isotropic': 0, 'analytic_nac': 1,
            'hop_analytic_nac': 2,
        }.get(getattr(self, 'rescale_provider', 'isotropic'), 0)
        record['rescale_source'] = getattr(
            self, '_last_rescale_source', default_rescale)
        record['rescale_gamma'] = getattr(self, '_last_rescale_gamma', np.nan)
        record['rescale_discriminant'] = getattr(
            self, '_last_rescale_discriminant', np.nan)
        record['hop_direction'] = getattr(
            self, '_last_hop_direction', np.zeros_like(coords))
        record['rng'] = self._last_hop_random
        record['e_unbiased_pot_hartree'] = getattr(
            self, '_unbiased_potential_energy', epot)
        record['e_pot_hartree'] = epot
        record['e_kin_hartree'] = ekin
        record['e_tot_hartree'] = epot + ekin
        record['droplet_energy_hartree'] = getattr(
            self, '_droplet_energy', 0.0)
        record['droplet_max_penetration_bohr'] = getattr(
            self, '_droplet_max_penetration', 0.0)
        record['droplet_active_count'] = getattr(
            self, '_droplet_active_count', 0)
        record['solute_com_energy_hartree'] = getattr(
            self, '_solute_com_energy', 0.0)
        record['solute_com_displacement_bohr'] = getattr(
            self, '_solute_com_displacement', 0.0)
        record['conservative_restraint_energy_hartree'] = getattr(
            self, '_conservative_restraint_energy', 0.0)
        record['thermostat_exchange_hartree'] = getattr(
            self, '_thermostat_exchange', 0.0)
        record['thermostat_exchange_cumulative_hartree'] = getattr(
            self, '_thermostat_exchange_cumulative', 0.0)
        record['thermostat_adjusted_energy_hartree'] = (
            epot + ekin - getattr(
                self, '_thermostat_exchange_cumulative', 0.0))
        droplet_force = getattr(self, '_droplet_force', None)
        if droplet_force is not None and np.asarray(droplet_force).shape == coords.shape:
            record['droplet_force_hartree_per_bohr'] = droplet_force
        solute_force = getattr(self, '_solute_com_force', None)
        if solute_force is not None and np.asarray(solute_force).shape == coords.shape:
            record['solute_com_force_hartree_per_bohr'] = solute_force
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
        if callable(state_energies):
            try:
                state_energies = state_energies()
            except (KeyError, TypeError, ValueError, AttributeError):
                state_energies = None
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
            overlap = np.asarray(self._last_state_overlap)
            record['state_overlap'] = overlap.real
            record['state_overlap_imag'] = overlap.imag
        if (self._last_overlap_tdc is not None
                and np.shape(self._last_overlap_tdc)
                == (trajectory_nstate, trajectory_nstate)):
            overlap_tdc = np.asarray(self._last_overlap_tdc)
            record['overlap_tdc_au'] = overlap_tdc.real
            record['overlap_tdc_imag_au'] = overlap_tdc.imag
        candidate_tdc = getattr(self, '_nacme_candidate_tdc', None)
        if (candidate_tdc is not None
                and np.shape(candidate_tdc)
                == (trajectory_nstate, trajectory_nstate)):
            record['gate_candidate_tdc_au'] = candidate_tdc
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
            record['gate_streak'] = gate.get(
                'consecutive_reference_failures', 0)
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
            record['tracking_raw_order'] = np.asarray(
                tracking.get('raw_order', tracking['order']), dtype=np.int64)
            record['tracking_lineage'] = np.asarray(
                tracking.get('lineage', tracking['order']), dtype=np.int64)
            record['tracking_phase'] = np.asarray(tracking['phase_step'], dtype=float)
            record['tracking_phase_initial'] = np.asarray(
                tracking.get('phase_initial', tracking['phase_step']), dtype=float)
            record['tracking_previous_phase_initial'] = np.asarray(
                tracking.get('previous_phase_initial',
                             tracking.get('phase_initial', tracking['phase_step'])),
                dtype=float)
            record['tracking_overlap'] = np.asarray(tracking['matched_overlap'], dtype=float)
            record['tracking_margin'] = np.asarray(tracking['margin'], dtype=float)

        record_bytes = record.tobytes(order='C')
        with open(self.trajectory_file, 'ab') as stream:
            stream.write(record_bytes)
            stream.flush()
            os.fsync(stream.fileno())
        self._trajectory_prefix_hasher.update(record_bytes)
        self._trajectory_prefix_bytes += len(record_bytes)
        self._trajectory_prefix_last_step = int(istep)
        self._trajectory_prefix_stat = self._trajectory_stat_identity()

    def _qmmm_forcefield_identity(self, value):
        """Canonicalize local force fields across relocated restart manifests."""
        if not isinstance(value, str):
            return value
        identity = []
        for item in value.replace(',', ' ').split():
            candidate = self._resolve_qmmm_aux_file(os.path.expanduser(item))
            if os.path.isfile(candidate):
                digest = hashlib.sha256()
                with open(candidate, 'rb') as stream:
                    for block in iter(lambda: stream.read(1024 * 1024), b''):
                        digest.update(block)
                identity.append({
                    'path': os.path.realpath(candidate),
                    'sha256': digest.hexdigest(),
                })
            else:
                try:
                    resource = resources.files('openmm.app').joinpath(
                        'data', *item.split('/'))
                    if not resource.is_file():
                        raise FileNotFoundError(item)
                    digest = hashlib.sha256(resource.read_bytes()).hexdigest()
                except (ImportError, FileNotFoundError, ModuleNotFoundError):
                    raise RuntimeError(
                        f'cannot fingerprint OpenMM force-field resource {item!r}')
                identity.append({'builtin': item, 'sha256': digest})
        return identity

    def _basis_definition_identity(self, value):
        """Fingerprint the file-backed basis definition used by BasisData."""
        path = self._resolved_basis_definition_file(value)
        if path is None:
            return value
        if not os.path.isfile(path):
            return {'path': os.path.realpath(path)}
        digest = hashlib.sha256()
        with open(path, 'rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(block)
        return {
            'path': os.path.realpath(path),
            'sha256': digest.hexdigest(),
        }

    @staticmethod
    def _external_file_identity(value):
        """Fingerprint a runtime file whose contents influence a trajectory."""
        if not isinstance(value, str) or not value.strip():
            return value
        path = os.path.realpath(os.path.abspath(
            os.path.expanduser(value.strip())))
        if not os.path.isfile(path):
            return {'path': path}
        digest = hashlib.sha256()
        with open(path, 'rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(block)
        return {'path': path, 'sha256': digest.hexdigest()}

    def _guess_settings_identity(self):
        """Return the immutable, pre-propagation guess configuration."""
        cached = getattr(self, '_restart_guess_identity', None)
        if cached is not None:
            return cached
        settings = dict(getattr(self.mol, 'config', {}).get('guess', {}))
        for key in ('file', 'file2'):
            if key in settings:
                settings[key] = self._external_file_identity(settings[key])
        self._restart_guess_identity = settings
        return settings

    def _restart_signature_matches(self, saved_signature, *, allow_smaller_dt=False,
                                   continuation_dt_limit=None):
        """Validate model identity and explicitly permitted continuation dt changes."""
        try:
            saved = json.loads(saved_signature)
            current = json.loads(self._restart_signature())
        except (TypeError, ValueError) as error:
            raise RuntimeError(
                'NAMD restart checkpoint has invalid signature metadata'
            ) from error
        if not isinstance(saved, dict) or not isinstance(current, dict):
            raise RuntimeError(
                'NAMD restart checkpoint has invalid signature metadata'
            )
        if allow_smaller_dt:
            old_dt, new_dt = saved.get('dt_fs'), current.get('dt_fs')
            if (not isinstance(old_dt, (int, float)) or isinstance(old_dt, bool)
                    or not isinstance(new_dt, (int, float)) or isinstance(new_dt, bool)
                    or not np.isfinite(old_dt) or not np.isfinite(new_dt)
                    or old_dt <= 0.0 or not 0.0 < new_dt
                    or (new_dt >= old_dt and (continuation_dt_limit is None
                        or new_dt > continuation_dt_limit or new_dt == old_dt))):
                return False
            # Explicit new-output continuation can reduce dt or return up to
            # the original dt recorded in the continuation history.
            # Every Hamiltonian, response, RNG, and acceptance setting stays bound.
            current['dt_fs'] = old_dt
        if saved == current:
            return True

        saved_guess = saved.get('guess_settings')
        current_guess = current.get('guess_settings')
        if not isinstance(saved_guess, dict) or not isinstance(current_guess, dict):
            return False
        if not self._as_bool(saved_guess.get('save_mol', False)):
            return False

        # save_mol rewrites guess.file, but never file2.  The configured path
        # must remain identical; only its live digest/existence may differ.
        saved_file = saved_guess.get('file')
        current_file = current_guess.get('file')
        if not (isinstance(saved_file, dict)
                and isinstance(current_file, dict)
                and saved_file.get('path') == current_file.get('path')):
            return False
        current_guess['file'] = saved_file
        if current != saved:
            return False

        # Subsequent trajectory reconciliation must reproduce the exact saved
        # signature rather than re-hashing the mutable result file.
        self._restart_guess_identity = copy.deepcopy(saved_guess)
        return True

    @staticmethod
    def _espf_environment_identity():
        """Return effective ESPF environment controls that alter QM/MM forces."""
        def enabled(name):
            return os.environ.get(name, '').strip() in ('1', 'on')

        def enabled_by_default(name):
            return os.environ.get(name, '').strip() not in ('0', 'off')

        def real_value(name, default):
            value = os.environ.get(name, '').strip()
            if not value:
                return float(default)
            # Fortran list-directed input accepts D exponents; normalize them
            # before producing the equivalent numeric restart identity.
            return float(value.replace('d', 'e').replace('D', 'E'))

        return {
            'rohf': enabled('ESPF_ROHF'),
            'legacy_gradient': enabled('ESPF_LEGACY'),
            'hard_grid': enabled('ESPF_HARD_GRID'),
            'keep_all': enabled('ESPF_KEEPALL'),
            'smooth': enabled_by_default('ESPF_SMOOTH'),
            'weight_derivative': enabled_by_default('ESPF_WDERIV'),
            'weight_scale': real_value('ESPF_WSCALE', 1.0),
            'switch_delta': real_value('ESPF_SWDELTA', 0.7),
            # 1.8 is the whole-molecule default; use 1.5 at a covalent QM/MM
            # boundary (openqp-devkit docs/espf_qmmm_switching.md, issue #260).
            'switch_scale': real_value('ESPF_SWSCALE', 1.8),
        }

    def _effective_tight_binding_settings(self):
        """Return TB settings after non-mutating default-model resolution."""
        method = str(self.mol.config['input'].get('method', '')).lower()
        if method not in ('dftb', 'xtb'):
            return method, {}
        config = self.mol.config
        if method == 'dftb':
            config = copy.deepcopy(config)
            from oqp.utils.input_checker import apply_dftb_model_default
            apply_dftb_model_default(config)
        return method, dict(config.get(method, {}))

    def _resolved_tight_binding_artifacts(self, settings=None):
        """Resolve the parameter and executable artifacts used at runtime."""
        config = getattr(self.mol, 'config', {})
        method = str(config.get('input', {}).get('method', '')).lower()
        if method not in ('dftb', 'xtb'):
            return {}
        if settings is None:
            _method, settings = self._effective_tight_binding_settings()
        # Resolve defaults through the same adapter methods used by the
        # calculation. This binds environment, installed-wheel, staged-lib,
        # and PATH fallbacks even when the input leaves these fields blank.
        if method == 'xtb':
            from oqp.library.openqp_xtb import OpenQPXTBAdapter
            adapter_class = OpenQPXTBAdapter
        else:
            from oqp.library.openqp_dftb import OpenQPDFTBAdapter
            adapter_class = OpenQPDFTBAdapter
        adapter = adapter_class.__new__(adapter_class)
        adapter.mol = self.mol
        adapter.dftb = settings
        backend = str(settings.get('backend', 'native')).strip().lower()
        artifacts = {'parameter_path': adapter._parameter_path()}
        if backend in ('native', 'auto'):
            artifacts['library_path'] = str(adapter._native_library_path())
        elif backend == 'probe':
            artifacts['executable'] = adapter._probe_executable()
        return artifacts

    def _tight_binding_identity(self):
        """Bind a restart to the active DFTB/xTB Hamiltonian definition."""
        cached = getattr(self, '_restart_tb_identity', None)
        if cached is not None:
            return cached
        method, settings = self._effective_tight_binding_settings()
        if method not in ('dftb', 'xtb'):
            return None
        source = (getattr(self.mol, 'oqp_input_source', None)
                  or getattr(self.mol, 'input_file', None))
        source_dir = (os.path.dirname(os.path.abspath(source))
                      if source else os.getcwd())
        resolved_defaults = self._resolved_tight_binding_artifacts(settings)

        for key in ('parameter_path', 'library_path', 'executable'):
            value = settings.get(key, '')
            resolved_value = resolved_defaults.get(key)
            if resolved_value is None and (
                    not isinstance(value, str) or not value.strip()):
                continue
            raw_value = str(resolved_value if resolved_value is not None
                            else value).strip()
            expanded = os.path.expanduser(raw_value)
            path = (expanded if os.path.isabs(expanded)
                    else os.path.join(source_dir, expanded))
            if not os.path.exists(path):
                settings[key] = {'configured': value, 'resolved': raw_value}
                continue
            digest = hashlib.sha256()
            real_path = os.path.realpath(path)
            files = ([real_path] if os.path.isfile(real_path) else [
                os.path.join(root, filename)
                for root, _dirs, names in os.walk(real_path)
                for filename in sorted(names)
            ])
            for filename in sorted(files):
                relative = os.path.relpath(filename, real_path)
                digest.update(relative.encode('utf-8') + b'\0')
                with open(filename, 'rb') as stream:
                    for block in iter(lambda: stream.read(1024 * 1024), b''):
                        digest.update(block)
            settings[key] = {'path': real_path, 'sha256': digest.hexdigest()}
        self._restart_tb_identity = {'section': method, 'settings': settings}
        return self._restart_tb_identity

    def _molecular_identity(self):
        """Return stable nuclear and QM/MM topology identity for sidecars."""
        cached = getattr(self, '_restart_molecular_identity', None)
        if cached is not None:
            return cached
        atoms = np.asarray(self.mol.get_atoms(), dtype=np.int64).reshape(-1)
        masses = np.asarray(self.mol.get_mass(), dtype=float).reshape(-1)
        if (atoms.size == 0 or masses.shape != atoms.shape
                or not np.all(np.isfinite(masses)) or np.any(masses <= 0.0)):
            raise RuntimeError('cannot identify NAMD molecule: invalid atoms/masses')
        identity = {
            'atomic_numbers': atoms.tolist(),
            'masses_amu': masses.tolist(),
        }
        if hasattr(self, 'qm_atoms'):
            identity['qm_atoms'] = np.asarray(
                self.qm_atoms, dtype=np.int64).reshape(-1).tolist()
        if hasattr(self, 'pdb'):
            topology = self.pdb.topology
            topology_atoms = []
            for atom in topology.atoms():
                element = getattr(atom, 'element', None)
                residue = getattr(atom, 'residue', None)
                topology_atoms.append({
                    'index': int(atom.index),
                    'atomic_number': int(getattr(element, 'atomic_number', 0) or 0),
                    'name': str(getattr(atom, 'name', '')),
                    'residue_index': int(getattr(residue, 'index', -1)),
                    'residue_name': str(getattr(residue, 'name', '')),
                })
            bonds = []
            for bond in topology.bonds():
                if hasattr(bond, 'atom1'):
                    atom1, atom2 = bond.atom1, bond.atom2
                else:
                    atom1, atom2 = bond[0], bond[1]
                bonds.append(sorted((int(atom1.index), int(atom2.index))))
            identity['qmmm_topology'] = {
                'atoms': topology_atoms,
                'bonds': sorted(bonds),
                'masses_amu': (
                    np.asarray(self.m_all, dtype=float).reshape(-1) / AMU_TO_AU
                ).tolist() if hasattr(self, 'm_all') else [],
            }
            box_getter = getattr(topology, 'getPeriodicBoxVectors', None)
            box = box_getter() if box_getter is not None else None
            if box is not None and hasattr(box, 'value_in_unit'):
                try:
                    from openmm import unit as openmm_unit
                except ImportError:
                    from simtk import unit as openmm_unit
                box = box.value_in_unit(openmm_unit.nanometer)
            identity['qmmm_topology']['periodic_box_vectors'] = (
                None if box is None else [
                    [float(component) for component in vector]
                    for vector in box
                ])
            qmmm = self.mol.config.get('qmmm', {})
            identity['qmmm_model'] = {
                key: str(qmmm.get(key, ''))
                for key in ('cutoff', 'embedding', 'frontier_scheme')
            }
            # These environment controls change the QM/MM charge, grid, or
            # gradient model in Python/Fortran.  Store effective values so a
            # checkpointed acceleration cannot cross a force-model boundary.
            identity['qmmm_model']['espf'] = (
                self._espf_environment_identity())
            forcefields = (qmmm.get('forcefield_files', '')
                           or qmmm.get('forcefield', ''))
            identity['qmmm_model']['forcefields'] = (
                self._qmmm_forcefield_identity(forcefields))
        self._restart_molecular_identity = identity
        return identity

    def _restart_signature(self):
        cfg = self.mol.config
        md = cfg.get('md', {})
        electronic_config = getattr(
            self, '_electronic_config_identity', None)
        if electronic_config is None:
            electronic_config = _electronic_config_identity(cfg)
        else:
            electronic_config = copy.deepcopy(electronic_config)
        # The detailed electronic snapshot must use the same canonical
        # identities as the dedicated restart fields below.  Keeping raw
        # ``file:relative`` spellings or the pre-resolution empty DFTB model
        # here would make two physically identical runs fail restart matching
        # even though their basis contents and effective Hamiltonian agree.
        electronic_input = electronic_config.setdefault('input', {})
        electronic_input['basis'] = self._basis_definition_identity(
            electronic_input.get('basis', ''))
        electronic_scf = electronic_config.setdefault('scf', {})
        if 'init_basis' in electronic_scf:
            electronic_scf['init_basis'] = self._basis_definition_identity(
                electronic_scf['init_basis'])
        method = str(electronic_input.get('method', '')).lower()
        if method in ('dftb', 'xtb'):
            if method == 'dftb':
                from oqp.utils.input_checker import apply_dftb_model_default
                apply_dftb_model_default(electronic_config)
            electronic_tb = dict(electronic_config.get(method, {}))
            # External artifacts are content-addressed by
            # ``_tight_binding_identity``; retaining path spellings here
            # would duplicate that identity non-canonically.
            for key in ('parameter_path', 'library_path', 'executable'):
                electronic_tb.pop(key, None)
            electronic_config[method] = electronic_tb
        electronic_config = _normalize_identity_value(electronic_config)
        scf_settings = dict(cfg.get('scf', {}))
        if 'init_basis' in scf_settings:
            scf_settings['init_basis'] = self._basis_definition_identity(
                scf_settings['init_basis'])
        guess_settings = self._guess_settings_identity()
        identity = {
            'molecule': self._molecular_identity(),
            'method': cfg['input'].get('method', ''),
            'functional': cfg['input'].get('functional', ''),
            'basis': self._basis_definition_identity(
                cfg['input'].get('basis', '')),
            'basis_library': cfg['input'].get('library', ''),
            'basis_ispher': cfg['input'].get('ispher', ''),
            'd4': cfg['input'].get('d4', False),
            'charge': cfg['input'].get('charge', ''),
            'input_multiplicity': cfg['input'].get('multiplicity', ''),
            'scf_type': cfg.get('scf', {}).get('type', ''),
            'scf_multiplicity': cfg.get('scf', {}).get('multiplicity', ''),
            'scf_settings': scf_settings,
            'scf_init_basis': self._basis_definition_identity(
                cfg.get('scf', {}).get('init_basis', 'none')),
            'guess_settings': guess_settings,
            'tdhf_type': cfg['tdhf'].get('type', ''),
            'tdhf_multiplicity': cfg['tdhf'].get('multiplicity', ''),
            'tdhf_settings': dict(cfg['tdhf']),
            'dftgrid_settings': dict(cfg.get('dftgrid', {})),
            'pcm_settings': dict(cfg.get('pcm', {})),
            # Bind the stable, effective symmetry controls but not the
            # geometry-dependent detected group, which may evolve during MD.
            'symmetry_input': dict(cfg.get('symmetry', {})),
            'symmetry_settings': {
                key: getattr(self.mol, 'symmetry_metadata', {}).get(
                    key, cfg.get('symmetry', {}).get(key, ''))
                for key in (
                    'status', 'enabled', 'requested_point_group',
                    'requested_subgroup', 'label_mo', 'label_states',
                    'label_modes', 'use_integral_symmetry',
                    'use_response_symmetry', 'strict', 'tolerance')
            },
            'tight_binding': self._tight_binding_identity(),
            'nstate': cfg['tdhf'].get('nstate', ''),
            'tlf': cfg['tdhf'].get('tlf', ''),
            'electronic_config': electronic_config,
            'dt_fs': self.dt_fs, 'seed': self.seed,
            'rng_stream': self.rng_stream,
            'substep': md.get('substep', ''),
            'decoherence': md.get('decoherence', ''),
            'edc_c': md.get('edc_c', ''), 'thrshe': md.get('thrshe', ''),
            # The resolved provider keeps rescale=auto checkpoints compatible
            # with runs that named the provider it resolves to.
            'tdc': md.get('tdc', ''),
            'rescale': getattr(self, 'rescale_provider', None) or md.get('rescale', ''),
            'trivial': md.get('trivial', ''),
            'trivial_thresh': md.get('trivial_thresh', ''),
            'first_hop_step': md.get('first_hop_step', ''),
            'soc_settings': {
                key: md.get(key, '') for key in (
                    'soc', 'soc_basis', 'soc_du_dt_corr',
                    'soc_tdc_grad_corr', 'grad_wthr', 'init_state',
                    'dt_adaptive', 'dt_min', 'dx_max', 'econs')
            },
            'soc': md.get('soc', ''), 'soc_basis': md.get('soc_basis', ''),
            'trajectory_representation': getattr(
                self, '_trajectory_representation', 'same_spin_adiabatic'),
            'odp': self._odp_provenance(),
            'system': getattr(
                self, '_restart_system_identity', {'kind': 'unavailable'}),
            'independent_controls': self._independent_settings_record(),
            'nac_align': cfg.get('nac', {}).get('align', ''),
            'gate_policy': {
                key: md.get(key, '') for key in (
                    'nacme_check', 'ba_gap_max', 'nacme_gate',
                    'nacme_gate_invariant_tol', 'nacme_gate_abs_tol',
                    'nacme_gate_rel_tol', 'nacme_gate_consecutive',
                    'nve_gate', 'nve_gate_abs_tol', 'nve_gate_step_tol',
                    'nve_gate_transition_tol', 'nve_gate_consecutive')
            },
        }
        # Controls that change which reference is followed, how frustrated
        # hops and energy discontinuities alter velocities, or how failed SCF
        # steps are recovered all change the trajectory, so a restart must
        # not silently continue under different settings.
        # Only non-default values are recorded, so checkpoints written before
        # these controls were bound (all at their defaults) remain loadable.
        controls = _non_default_md_controls(md, _TRAJECTORY_CONTROL_KEYS)
        if controls:
            identity['trajectory_controls'] = controls
        return json.dumps(identity, sort_keys=True, separators=(',', ':'))

    def _tracking_state_count(self):
        """Number of response roots represented by Molecule tracking tags."""
        return self.nstate

    def _trajectory_state_energies(self):
        """Return energies in the trajectory's electronic propagation basis."""
        return np.asarray(
            self.mol.data['OQP::td_energies'], dtype=float
        ).reshape(-1)[:self.nstate]

    def _restart_extra_payload(self):
        """Subclass hook for representation-specific checkpoint arrays."""
        return {}

    def _load_restart_extra(self, saved, prev_data=None):
        """Subclass hook for validating representation-specific arrays."""
        del saved, prev_data
        return {}

    def _restore_restart_extra(self, extra):
        """Subclass hook for restoring representation-specific state."""
        if extra:
            raise RuntimeError('unexpected NAMD restart representation state')

    def _validate_restart_state(self, nuclear_state, coef, active, prev_xyz,
                                prev_data, *, context):
        """Reject invalid nuclear/electronic state before checkpoint use."""
        coordinates, velocities, acceleration = (
            np.asarray(value, dtype=float) for value in nuclear_state)
        if (coordinates.ndim != 2 or coordinates.shape[1:] != (3,)
                or coordinates.shape != velocities.shape
                or coordinates.shape != acceleration.shape
                or coordinates.size == 0
                or not all(np.all(np.isfinite(value)) for value in (
                    coordinates, velocities, acceleration))):
            raise RuntimeError(
                f'{context} contains an invalid NAMD nuclear state')
        expected_natom = int(getattr(
            self, 'natom_all', len(self._molecular_identity()['atomic_numbers'])))
        if coordinates.shape[0] != expected_natom:
            raise RuntimeError(
                f'{context} nuclear atom count {coordinates.shape[0]} does not '
                f'match the current system ({expected_natom})')

        coef = np.asarray(coef, dtype=np.complex128)
        if (coef.shape != (self.nstate,)
                or not np.all(np.isfinite(coef.real))
                or not np.all(np.isfinite(coef.imag))):
            raise RuntimeError(f'{context} contains invalid electronic coefficients')
        norm = float(np.vdot(coef, coef).real)
        if not np.isfinite(norm) or abs(norm - 1.0) > 1.0e-6:
            raise RuntimeError(
                f'{context} electronic coefficient norm is {norm!r}, not 1')
        if not 1 <= int(active) <= self.nstate:
            raise RuntimeError(
                f'{context} active state {active} is outside 1..{self.nstate}')

        prev_xyz = np.asarray(prev_xyz, dtype=float).reshape(-1)
        expected_qm_size = 3 * len(self._molecular_identity()['atomic_numbers'])
        if (prev_xyz.size != expected_qm_size
                or not np.all(np.isfinite(prev_xyz))):
            raise RuntimeError(f'{context} contains invalid previous QM coordinates')
        for key, raw_value in prev_data.items():
            value = np.asarray(raw_value)
            if value.dtype == object:
                raise TypeError(f'NAMD restart cannot serialize ragged tag {key!r}')
            if value.dtype.kind in 'biufc' and not np.all(np.isfinite(value)):
                raise RuntimeError(
                    f'{context} contains non-finite previous-state tag {key!r}')
            if key.startswith('OQP::state_tracking_'):
                tracking_nstate = self._tracking_state_count()
                scalar_tags = {'OQP::state_tracking_output_reordered'}
                expected_shape = ((1,) if key in scalar_tags
                                  else (tracking_nstate,))
                if value.shape != expected_shape:
                    raise RuntimeError(
                        f'{context} contains invalid tracking tag {key!r} '
                        f'shape {value.shape}; expected {expected_shape}')
                if key in {
                        'OQP::state_tracking_order',
                        'OQP::state_tracking_raw_order'}:
                    order = np.asarray(value, dtype=np.int64)
                    if not np.array_equal(
                            np.sort(order), np.arange(tracking_nstate)):
                        raise RuntimeError(
                            f'{context} contains invalid tracking permutation '
                            f'{key!r}')
                if key == 'OQP::state_tracking_lineage':
                    if (value.dtype.kind not in 'iu'
                            or np.unique(value).size != tracking_nstate):
                        raise RuntimeError(
                            f'{context} contains invalid tracking lineage IDs')
                if key in {
                        'OQP::state_tracking_phase_step',
                        'OQP::state_tracking_phase_initial',
                        'OQP::state_tracking_previous_phase_initial'}:
                    if not np.allclose(np.abs(value), 1.0, atol=1.0e-12,
                                       rtol=0.0):
                        raise RuntimeError(
                            f'{context} contains invalid tracking phase {key!r}')
        return coordinates, velocities, acceleration, coef, prev_xyz

    @staticmethod
    def _checkpoint_optional(payload, name, value):
        if value is None:
            payload[f'has_{name}'] = np.array([0], dtype=np.int8)
            payload[name] = np.empty(0, dtype=np.float64)
        else:
            payload[f'has_{name}'] = np.array([1], dtype=np.int8)
            payload[name] = np.asarray(value)

    def _validate_restart_histories(self, histories, *, context):
        """Validate optional Baeck-An and NVE state without reshaping it."""
        shapes = {
            'ba_energy_left': (self.nstate,),
            'ba_energy_center': (self.nstate,),
            'ba_tdc_left': (self.nstate, self.nstate),
            'ba_dt_left': (),
            'nve_reference_energy': (),
            'nve_previous_energy': (),
            'analytic_tdc_previous': (self.nstate, self.nstate),
            'etot_prev': (),
            'disc_energy_absorbed': (),
        }
        validated = {}
        for name, expected in shapes.items():
            value = histories.get(name)
            if value is None:
                validated[name] = None
                continue
            array = np.asarray(value)
            if array.shape != expected or not np.all(np.isfinite(array)):
                raise RuntimeError(f'{context} contains invalid {name}')
            if name == 'ba_dt_left' and float(array) <= 0.0:
                raise RuntimeError(f'{context} contains invalid ba_dt_left')
            validated[name] = float(array) if expected == () else array.copy()
        ba_names = (
            'ba_energy_left', 'ba_energy_center', 'ba_tdc_left', 'ba_dt_left')
        ba_present = [validated[name] is not None for name in ba_names]
        if any(ba_present) and not all(ba_present):
            raise RuntimeError(
                f'{context} contains incomplete Baeck-An history')
        nve_names = ('nve_reference_energy', 'nve_previous_energy')
        nve_present = [validated[name] is not None for name in nve_names]
        if any(nve_present) and not all(nve_present):
            raise RuntimeError(f'{context} contains incomplete NVE history')
        return validated

    def _save_restart(self, istep, coordinates, velocities, acceleration):
        """Atomically save all state needed for phase-continuous continuation."""
        if istep % self.restart_interval != 0 and istep != self.nstep:
            return
        return self._run_io_collective(lambda: self._save_restart_on_io_rank(
            istep, coordinates, velocities, acceleration))

    def _save_restart_on_io_rank(self, istep, coordinates, velocities,
                                 acceleration):
        """Validate and atomically write a checkpoint on rank zero."""
        if self.prev_data is None or self.prev_xyz is None:
            return
        nuclear_state = (coordinates, velocities, acceleration)
        coordinates, velocities, acceleration, coef, prev_xyz = (
            self._validate_restart_state(
                nuclear_state, self.coef, self.active, self.prev_xyz,
                self.prev_data,
                context=(f'refusing to overwrite the last-good NAMD restart '
                         f'at step {istep}: state'),
            )
        )
        prev_keys = sorted(self.prev_data)
        histories = self._validate_restart_histories({
            'ba_energy_left': self._ba_energy_left,
            'ba_energy_center': self._ba_energy_center,
            'ba_tdc_left': self._ba_tdc_left,
            'ba_dt_left': self._ba_dt_left,
            'analytic_tdc_previous': getattr(
                self, '_analytic_tdc_previous', None),
            'nve_reference_energy': self._nve_reference_energy,
            'nve_previous_energy': self._nve_previous_energy,
            # The energy-recovery baseline is also kept for NVT, where the NVE
            # history above stays empty.
            'etot_prev': getattr(self, '_etot_prev', None),
            'disc_energy_absorbed': getattr(self, '_disc_energy_absorbed', 0.0),
        }, context=(f'refusing to overwrite the last-good NAMD restart at '
                    f'step {istep}: history'))
        trajectory_prefix = self._trajectory_checkpoint_identity(istep)
        payload = {
            'schema_version': np.array([NAMD_RESTART_SCHEMA_VERSION], dtype=np.int64),
            'signature': np.array([self._restart_signature()]),
            'time_origin_fs': np.array([getattr(self, '_time_origin_fs', 0.0)]),
            'continuation_provenance_json': np.array([json.dumps(
                getattr(self, '_continuation_provenance', None), sort_keys=True)]),
            'step': np.array([istep], dtype=np.int64),
            'time_fs': np.array([self._physical_time_fs(istep)]),
            'active': np.array([self.active], dtype=np.int64),
            'rng_step': np.array([self._rng_step], dtype=np.int64),
            'gate_failures': np.array([self._nacme_gate_failures], dtype=np.int64),
            'nve_failures': np.array([self._nve_gate_failures], dtype=np.int64),
            'independent_controls_json': np.array([
                json.dumps(self._independent_settings_record(), sort_keys=True,
                           separators=(',', ':'))
            ]),
            'droplet_energy': np.array([
                getattr(self, '_droplet_energy', 0.0)], dtype=np.float64),
            'droplet_max_penetration': np.array([
                getattr(self, '_droplet_max_penetration', 0.0)],
                dtype=np.float64),
            'droplet_active_count': np.array([
                getattr(self, '_droplet_active_count', 0)], dtype=np.int64),
            'solute_com_energy': np.array([
                getattr(self, '_solute_com_energy', 0.0)], dtype=np.float64),
            'solute_com_displacement': np.array([
                getattr(self, '_solute_com_displacement', 0.0)],
                dtype=np.float64),
            'thermostat_exchange': np.array([
                getattr(self, '_thermostat_exchange', 0.0)], dtype=np.float64),
            'thermostat_exchange_cumulative': np.array([
                getattr(self, '_thermostat_exchange_cumulative', 0.0)],
                dtype=np.float64),
            'coordinates': np.asarray(coordinates, dtype=np.float64),
            'velocities': np.asarray(velocities, dtype=np.float64),
            'acceleration': np.asarray(acceleration, dtype=np.float64),
            'coef_real': np.asarray(coef.real, dtype=np.float64),
            'coef_imag': np.asarray(coef.imag, dtype=np.float64),
            'prev_xyz': np.asarray(prev_xyz, dtype=np.float64),
            'prev_keys': np.asarray(prev_keys, dtype=np.str_),
            'odp_provenance': np.array([
                json.dumps(self._odp_provenance(), sort_keys=True,
                           separators=(',', ':'))
            ]),
            'trajectory_prefix_bytes': np.array(
                [trajectory_prefix['bytes']], dtype=np.int64),
            'trajectory_prefix_sha256': np.array(
                [trajectory_prefix['sha256']]),
        }
        for index, key in enumerate(prev_keys):
            value = np.asarray(self.prev_data[key])
            payload[f'prev_{index}'] = value
        for name, value in histories.items():
            self._checkpoint_optional(payload, name, value)
        extra_payload = self._restart_extra_payload()
        duplicate = set(payload).intersection(extra_payload)
        if duplicate:
            raise RuntimeError(
                'duplicate representation checkpoint fields: '
                + ', '.join(sorted(duplicate)))
        payload.update(extra_payload)

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
    def _rebase_restart_spec_paths(self, spec, source_dir):
        """Make input-owned paths stable when the manifest moves to log_dir."""
        from oqp.utils.oqp_input import CallSpec, CalculationSpec

        def absolute_path(value):
            if not isinstance(value, str) or not value.strip():
                return value
            expanded = os.path.expanduser(value.strip())
            if os.path.isabs(expanded):
                return os.path.normpath(expanded)
            return os.path.normpath(os.path.abspath(os.path.join(
                source_dir, expanded)))

        def runtime_path(value):
            """Resolve paths exactly as the live guess readers do."""
            if not isinstance(value, str) or not value.strip():
                return value
            return os.path.realpath(os.path.abspath(
                os.path.expanduser(value.strip())))

        def geometry_path(value):
            if not isinstance(value, str) or '\n' in value or '\r' in value:
                return value
            candidate = value.strip()
            if (os.path.splitext(candidate)[1].lower() in ('.xyz', '.pdb')
                    or os.path.isfile(os.path.join(source_dir, candidate))):
                return absolute_path(value)
            return value

        def search_path_list(value):
            if not isinstance(value, str):
                return value
            entries = [item for item in value.replace(',', ' ').split() if item]
            resolved = []
            for item in entries:
                candidate = self._resolve_qmmm_aux_file(
                    os.path.expanduser(item))
                # Preserve OpenMM built-in force-field names; only local files
                # are made absolute.  The shared resolver preserves the
                # runtime's existing-CWD-file precedence over input_dir.
                resolved.append(os.path.abspath(candidate)
                                if os.path.isfile(candidate) else item)
            return ' '.join(resolved)

        def basis_path(value):
            if not isinstance(value, str) or not value.startswith('file:'):
                return value
            return 'file:' + absolute_path(value[len('file:'):])

        options = dict(spec.options)
        for key in ('geom', 'geom2'):
            if key in options:
                options[key] = geometry_path(options[key])

        model_options = dict(spec.model_options)
        for key in ('parameter_path', 'library_path'):
            if key in model_options:
                model_options[key] = absolute_path(model_options[key])

        driver_kwargs = dict(spec.driver.kwargs)
        velocity = driver_kwargs.get('velocity')
        if (isinstance(velocity, str)
                and velocity.strip().lower() not in ('maxwell', 'zero')):
            driver_kwargs['velocity'] = absolute_path(velocity)
        driver = CallSpec(
            spec.driver.name, spec.driver.args, driver_kwargs,
            spec.driver.explicit)

        single_path_keys = {
            'input': {'system', 'system2'},
            'neb': {'product'},
            'guess': {'file', 'file2'},
            'dftb': {'parameter_path', 'library_path'},
            'geometric': {'constraints_file'},
            'oqp': {'neb_output'},
            'qmmm': {
                'pdb_file', 'qm_atoms_xyz', 'trajectory_file', 'log_file',
                'energy_file',
            },
        }
        modifiers = []
        for call in spec.modifiers:
            kwargs = dict(call.kwargs)
            for key in single_path_keys.get(call.name, set()):
                if key in kwargs:
                    kwargs[key] = (runtime_path(kwargs[key])
                                   if call.name == 'guess'
                                   else absolute_path(kwargs[key]))
            if call.name == 'qmmm':
                for key in ('forcefield', 'forcefield_files'):
                    if key in kwargs:
                        kwargs[key] = search_path_list(kwargs[key])
            if call.name == 'scf' and 'init_basis' in kwargs:
                kwargs['init_basis'] = basis_path(kwargs['init_basis'])
            modifiers.append(CallSpec(
                call.name, call.args, kwargs, call.explicit))

        return CalculationSpec(
            spec.model, spec.functional, basis_path(spec.basis),
            model_options, options,
            driver, tuple(modifiers), spec.source_text)

    def _write_restart_manifest(self):
        """Write a directly runnable per-job restart manifest."""
        if getattr(self, '_restart_manifest_written', False):
            return
        source = (getattr(self.mol, 'oqp_input_source', None)
                  or getattr(self.mol, 'input_file', None))
        canonical = str(getattr(self.mol, 'oqp_canonical_input', '') or '').strip()
        if not canonical:
            if source and str(source).lower().endswith('.oqp') and os.path.isfile(source):
                with open(source, 'r', encoding='utf-8') as stream:
                    canonical = stream.read().strip()
        if not canonical:
            if not getattr(self, '_manifest_notice_logged', False):
                dump_log(
                    self.mol,
                    title=(
                        'NAMD checkpoint saved, but restart.oqp was not generated: '
                        'the run did not originate from canonical .oqp input'),
                )
                self._manifest_notice_logged = True
            return
        from oqp.utils.oqp_input import (
            CallSpec, CalculationSpec, parse_canonical_oqp,
            rebase_calculation_paths, render_canonical_oqp,
        )
        spec = parse_canonical_oqp(canonical)
        if spec.driver.name != 'namd':
            raise ValueError('cannot create a restart manifest from a non-NAMD request')
        if source:
            spec = rebase_calculation_paths(
                spec, source_dir=os.path.dirname(os.path.abspath(source)))
        directory = os.path.dirname(self.restart_manifest_file) or '.'
        source_dir = (os.path.dirname(os.path.abspath(source))
                      if source else os.getcwd())
        spec = self._rebase_restart_spec_paths(spec, source_dir)
        kwargs = dict(spec.driver.kwargs)
        kwargs.pop('continuation_checkpoint', None)
        kwargs.pop('continuation_trajectory', None)
        kwargs.update({
            'restart': True,
            # Freeze a date-derived default so restarting on a later day keeps
            # exactly the original stochastic stream and signature.
            'seed': self.seed,
            'restart_file': os.path.relpath(self.restart_file, directory),
            'trajectory_file': os.path.relpath(self.trajectory_file, directory),
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
            self._restart_manifest_written = True
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)

    def _physical_time_fs(self, istep):
        if getattr(self, 'dt_adaptive', False):
            return self._t_fs
        return getattr(self, '_time_origin_fs', 0.0) + istep*self.dt_fs

    def _load_continuation_on_io_rank(self):
        """Validate an immutable source and start isolated changed-dt outputs."""
        self._validate_sidecar_paths()
        outputs = (self.trajectory_file, self.restart_file,
                   self.restart_manifest_file, self.zpredict_audit_file)
        sources = (self.continuation_checkpoint, self.continuation_trajectory)
        for output in outputs:
            if os.path.lexists(output):
                raise ValueError('local continuation requires new output paths: ' + output)
        for output in outputs + (self.mol.log,):
            for source in sources:
                if (os.path.realpath(output) == os.path.realpath(source)
                        or (os.path.exists(output) and os.path.exists(source)
                            and os.path.samefile(output, source))):
                    raise ValueError('continuation source aliases an output')
        with open(self.continuation_checkpoint, 'rb') as stream:
            checkpoint_hash = _sha256_stream(stream)
        payload = self._load_restart_on_io_rank(
            self.continuation_checkpoint, allow_smaller_dt=True)
        required_tags = {'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A',
                         'OQP::E_MO_B', 'OQP::DM_A', 'OQP::DM_B',
                         'OQP::FOCK_A', 'OQP::FOCK_B', 'OQP::SM',
                         'OQP::td_bvec_mo', 'OQP::td_energies',
                         'OQP::state_tracking_phase_initial',
                         'OQP::state_tracking_lineage'}
        if required_tags.difference(payload['prev_data']):
            raise ValueError('continuation checkpoint lacks reference/phase history')
        if self.nstep <= payload['step']:
            raise ValueError('continuation nstep must exceed the saved absolute step index')
        if (payload['optional']['nve_previous_energy'] is None
                and payload['optional'].get('etot_prev') is None):
            raise ValueError('local continuation requires the saved total-energy history')
        scanned = self._scan_trajectory_prefix(
            payload['step'], path=self.continuation_trajectory,
            expected_signature=payload['signature'])
        if (scanned['last_step'] != payload['step']
                or {k: scanned[k] for k in ('bytes', 'sha256')} != payload['trajectory_prefix']):
            raise ValueError('continuation trajectory does not match checkpoint prefix')
        header, records = read_namd_trajectory(self.continuation_trajectory)
        anchor = np.array(records[scanned['records']-1:scanned['records']], copy=True)
        del records
        for field, value in (
                ('coordinates_bohr', payload['coordinates']),
                ('velocities_au', payload['velocities']),
                ('coef_real', payload['coef'].real), ('coef_imag', payload['coef'].imag)):
            if not np.allclose(anchor[field][0].reshape(-1), np.asarray(value).reshape(-1),
                               rtol=0.0, atol=1e-14):
                raise ValueError('continuation checkpoint and trajectory state disagree')
        if (int(anchor['active'][0]) != payload['active']
                or not np.isclose(float(anchor['time_fs'][0]), payload['time_fs'],
                                  rtol=0.0, atol=1e-10)):
            raise ValueError('continuation checkpoint and trajectory time/state disagree')
        with open(self.continuation_checkpoint, 'rb') as stream:
            if _sha256_stream(stream) != checkpoint_hash:
                raise ValueError('continuation checkpoint changed during validation')
        provenance = {
            'checkpoint': self.continuation_checkpoint,
            'checkpoint_sha256': checkpoint_hash,
            'trajectory': self.continuation_trajectory,
            'committed_prefix': payload['trajectory_prefix'],
            'source_step': payload['step'], 'source_time_fs': payload['time_fs'],
            'source_rng_step': payload['rng_step'],
            'source_dt_fs': json.loads(payload['signature'])['dt_fs'],
            'new_dt_fs': self.dt_fs,
            'parent': payload.get('continuation_provenance'),
        }
        payload['time_origin_fs'] = payload['time_fs'] - payload['step']*self.dt_fs
        payload['continuation_provenance'] = provenance
        header.update(signature=self._restart_signature(), continuation=provenance,
                      time_origin_fs=payload['time_origin_fs'])
        encoded = json.dumps(header, sort_keys=True).encode('utf-8')
        data = NAMD_TRAJECTORY_MAGIC + struct.pack('<Q', len(encoded)) + encoded + anchor.tobytes()
        with open(self.trajectory_file, 'xb') as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        with open(self.zpredict_audit_file, 'xb'):
            pass
        self._remember_trajectory_prefix(self._scan_trajectory_prefix(payload['step']))
        return payload

    def _load_restart(self):
        """Collectively load and restore a representation-aware checkpoint."""
        continuing = bool(getattr(self, 'continuation_checkpoint', ''))
        if not self.restart_requested and not continuing:
            return None
        payload = self._run_io_collective_result(
            self._load_continuation_on_io_rank if continuing else self._load_restart_on_io_rank)
        self.prev_data = payload['prev_data']
        self.prev_xyz = payload['prev_xyz']
        self.mol.put_data(self.prev_data)
        self.active = payload['active']
        self.coef = payload['coef']
        self._rng_step = payload['rng_step']
        self._last_hop_random = np.nan
        self._nacme_gate_failures = payload['gate_failures']
        self._nve_gate_failures = payload['nve_failures']
        self._t_fs = payload['time_fs']
        self._time_origin_fs = payload.get('time_origin_fs', 0.0)
        self._continuation_provenance = payload.get('continuation_provenance')
        for name, value in payload['optional'].items():
            setattr(self, f'_{name}', value)
        self._restore_restart_extra(payload['extra'])
        for name, value in payload['independent_state'].items():
            setattr(self, f'_{name}', value)
        self._conservative_restraint_energy = (
            self._droplet_energy + self._solute_com_energy)
        if getattr(self, '_etot_prev', None) is None:
            self._etot_prev = self._nve_previous_energy
        if getattr(self, '_disc_energy_absorbed', None) is None:
            self._disc_energy_absorbed = 0.0
        if continuing:
            self._run_io_collective(lambda: self._save_restart_on_io_rank(
                payload['step'], payload['coordinates'], payload['velocities'],
                payload['acceleration']))
        else:
            self._reconcile_trajectory_with_restart(
                payload['step'], payload['trajectory_prefix'])
        dump_log(
            self.mol,
            title=(f'NAMD restart loaded: step={payload["step"]} '
                   f'file={self.restart_file} phase_history=restored '
                   f'rng=({self.seed},{self.rng_stream},step)'),
        )
        return {
            key: payload[key] for key in (
                'step', 'coordinates', 'velocities', 'acceleration')
        }

    def _load_restart_on_io_rank(self, checkpoint_file=None, *, allow_smaller_dt=False):
        """Read and validate a checkpoint on rank zero without pickle data."""
        checkpoint_file = checkpoint_file or self.restart_file
        if not os.path.isfile(checkpoint_file):
            raise FileNotFoundError(f'NAMD restart file not found: {checkpoint_file}')
        with np.load(checkpoint_file, allow_pickle=False) as saved:
            version = self._restart_integer(saved, 'schema_version')
            if version != NAMD_RESTART_SCHEMA_VERSION:
                raise ValueError(f'unsupported NAMD restart schema {version}')
            signature_array = np.asarray(saved['signature'])
            if signature_array.shape != (1,) or signature_array.dtype.kind not in 'SU':
                raise RuntimeError(
                    'NAMD restart checkpoint has invalid signature metadata')
            signature = str(signature_array[0])
            odp_array = np.asarray(saved['odp_provenance'])
            current_odp = json.dumps(
                self._odp_provenance(), sort_keys=True, separators=(',', ':'))
            if (odp_array.shape != (1,) or odp_array.dtype.kind not in 'SU'
                    or str(odp_array[0]) != current_odp):
                raise ValueError('NAMD restart ODP definition/metric mismatch')
            step = self._restart_integer(saved, 'step')
            active = self._restart_integer(saved, 'active', minimum=1)
            rng_step = self._restart_integer(saved, 'rng_step')
            gate_failures = self._restart_integer(saved, 'gate_failures')
            nve_failures = self._restart_integer(saved, 'nve_failures')
            time_fs = self._restart_float(saved, 'time_fs', minimum=0.0)
            time_origin_fs = (self._restart_float(saved, 'time_origin_fs')
                              if 'time_origin_fs' in saved else 0.0)
            continuation_provenance = None
            if 'continuation_provenance_json' in saved:
                provenance_array = np.asarray(saved['continuation_provenance_json'])
                if provenance_array.shape != (1,) or provenance_array.dtype.kind not in 'SU':
                    raise ValueError('invalid continuation provenance')
                continuation_provenance = json.loads(str(provenance_array[0]))
                if continuation_provenance is not None and not isinstance(continuation_provenance, dict):
                    raise ValueError('invalid continuation provenance')
            dt_limit = None
            if allow_smaller_dt and continuation_provenance is not None:
                # Walk the recorded chain rather than treating an arbitrary
                # increase from a fixed-dt checkpoint as a return.
                ancestor = continuation_provenance
                child_dt = json.loads(signature).get('dt_fs')
                while ancestor is not None:
                    if not isinstance(ancestor, dict):
                        raise ValueError('invalid continuation dt history')
                    old_dt = ancestor.get('source_dt_fs')
                    new_dt = ancestor.get('new_dt_fs')
                    if (not isinstance(old_dt, (int, float)) or isinstance(old_dt, bool)
                            or not np.isfinite(old_dt) or old_dt <= 0.0
                            or new_dt != child_dt):
                        raise ValueError('invalid continuation dt history')
                    dt_limit = old_dt
                    child_dt = old_dt
                    ancestor = ancestor.get('parent')
            if not self._restart_signature_matches(
                    signature, allow_smaller_dt=allow_smaller_dt,
                    continuation_dt_limit=dt_limit):
                raise ValueError('NAMD restart electronic model/RNG/time-step mismatch')
            saved_dt = json.loads(signature).get('dt_fs')
            if (not self.dt_adaptive and (not isinstance(saved_dt, (int, float))
                    or not np.isclose(time_fs, time_origin_fs + step*saved_dt,
                                      rtol=0.0, atol=1e-10))):
                raise ValueError('NAMD restart physical time is inconsistent with its step and dt')
            trajectory_prefix_bytes = self._restart_integer(
                saved, 'trajectory_prefix_bytes')
            trajectory_digest_array = np.asarray(
                saved['trajectory_prefix_sha256'])
            if (trajectory_digest_array.shape != (1,)
                    or trajectory_digest_array.dtype.kind not in 'SU'):
                raise RuntimeError(
                    'NAMD restart checkpoint has invalid dense trajectory '
                    'prefix metadata')
            trajectory_prefix_sha256 = str(
                trajectory_digest_array[0]).lower()
            if not re.fullmatch(r'[0-9a-f]{64}', trajectory_prefix_sha256):
                raise RuntimeError(
                    'NAMD restart checkpoint has invalid dense trajectory '
                    'prefix metadata')
            if rng_step > step:
                raise RuntimeError(
                    'NAMD restart checkpoint has invalid rng_step metadata')
            if gate_failures > step + 1 or nve_failures > step + 1:
                raise RuntimeError(
                    'NAMD restart checkpoint has implausible gate failure streak')
            keys = [str(key) for key in saved['prev_keys']]
            prev_data = {
                key: np.array(saved[f'prev_{index}'], copy=True)
                for index, key in enumerate(keys)
            }
            prev_xyz = np.array(saved['prev_xyz'], copy=True)
            coef_real = np.array(saved['coef_real'], copy=True)
            coef_imag = np.array(saved['coef_imag'], copy=True)
            if (coef_real.shape != (self.nstate,)
                    or coef_imag.shape != (self.nstate,)
                    or not np.all(np.isfinite(coef_real))
                    or not np.all(np.isfinite(coef_imag))):
                raise RuntimeError(
                    'NAMD restart checkpoint contains invalid serialized '
                    'electronic coefficient vectors')
            coef = coef_real + 1j*coef_imag
            nuclear_state = tuple(np.array(saved[name], copy=True) for name in (
                'coordinates', 'velocities', 'acceleration'))
            coordinates, velocities, acceleration, coef, prev_xyz = (
                self._validate_restart_state(
                    nuclear_state, coef, active, prev_xyz, prev_data,
                    context='NAMD restart checkpoint',
                )
            )
            optional = {}
            for name in ('ba_energy_left', 'ba_energy_center', 'ba_tdc_left',
                         'ba_dt_left', 'nve_reference_energy',
                         'nve_previous_energy', 'analytic_tdc_previous',
                         'etot_prev', 'disc_energy_absorbed'):
                if f'has_{name}' not in saved and name in (
                        'analytic_tdc_previous', 'etot_prev',
                        'disc_energy_absorbed'):
                    # Checkpoints written before these histories existed.
                    optional[name] = None
                    continue
                present = np.asarray(saved[f'has_{name}'])
                if (present.shape != (1,) or int(present[0]) not in (0, 1)):
                    raise RuntimeError(
                        f'NAMD restart checkpoint has invalid {name} marker')
                value = (np.array(saved[name], copy=True)
                         if int(present[0]) else None)
                optional[name] = value
            optional = self._validate_restart_histories(
                optional, context='NAMD restart checkpoint')
            extra = self._load_restart_extra(saved, prev_data)
            settings = np.asarray(saved['independent_controls_json'])
            if (settings.shape != (1,) or settings.dtype.kind not in 'SU'
                    or str(settings[0]) != json.dumps(
                        self._independent_settings_record(), sort_keys=True,
                        separators=(',', ':'))):
                raise ValueError(
                    'NAMD restart droplet/restraint/thermostat mismatch')
            independent_state = {
                'droplet_energy': self._restart_float(
                    saved, 'droplet_energy', minimum=0.0),
                'droplet_max_penetration': self._restart_float(
                    saved, 'droplet_max_penetration', minimum=0.0),
                'droplet_active_count': self._restart_integer(
                    saved, 'droplet_active_count'),
                'solute_com_energy': self._restart_float(
                    saved, 'solute_com_energy', minimum=0.0),
                'solute_com_displacement': self._restart_float(
                    saved, 'solute_com_displacement', minimum=0.0),
                'thermostat_exchange': self._restart_float(
                    saved, 'thermostat_exchange'),
                'thermostat_exchange_cumulative': self._restart_float(
                    saved, 'thermostat_exchange_cumulative'),
            }
            result = {
                'step': step,
                'coordinates': coordinates,
                'velocities': velocities,
                'acceleration': acceleration,
                'prev_data': prev_data,
                'prev_xyz': prev_xyz,
                'active': active,
                'coef': coef,
                'rng_step': rng_step,
                'gate_failures': gate_failures,
                'nve_failures': nve_failures,
                'time_fs': time_fs,
                'time_origin_fs': time_origin_fs,
                'continuation_provenance': continuation_provenance,
                'signature': signature,
                'trajectory_prefix': {
                    'bytes': trajectory_prefix_bytes,
                    'sha256': trajectory_prefix_sha256,
                },
                'optional': optional,
                'extra': extra,
                'independent_state': independent_state,
            }
        return result

    @staticmethod
    def _restart_integer(saved, name, *, minimum=0):
        """Read one exact integer checkpoint field without coercion."""
        value = np.asarray(saved[name])
        if (value.shape != (1,) or value.dtype.kind not in 'iu'
                or int(value[0]) < minimum):
            raise RuntimeError(
                f'NAMD restart checkpoint has invalid {name} metadata')
        return int(value[0])

    @staticmethod
    def _restart_float(saved, name, *, minimum=None):
        """Read one finite scalar checkpoint field without broadcasting."""
        value = np.asarray(saved[name])
        if value.shape != (1,) or value.dtype.kind not in 'fiu':
            raise RuntimeError(
                f'NAMD restart checkpoint has invalid {name} metadata')
        scalar = float(value[0])
        if not np.isfinite(scalar) or (minimum is not None and scalar < minimum):
            raise RuntimeError(
                f'NAMD restart checkpoint has invalid {name} metadata')
        return scalar

    def _trajectory_stat_identity(self):
        """Return cheap metadata that detects sidecar replacement/mutation."""
        try:
            status = os.stat(self.trajectory_file)
        except OSError as error:
            raise RuntimeError(
                'NAMD dense trajectory disappeared before checkpoint') from error
        return (
            int(status.st_dev), int(status.st_ino), int(status.st_size),
            int(status.st_mtime_ns),
        )

    def _require_unchanged_trajectory_prefix(self):
        """Reject external sidecar changes before using the cached digest."""
        expected = getattr(self, '_trajectory_prefix_stat', None)
        if expected is None or self._trajectory_stat_identity() != expected:
            raise RuntimeError(
                'NAMD dense trajectory changed outside the active writer')

    def _remember_trajectory_prefix(self, scanned):
        """Install a validated incremental digest after start or restart."""
        self._trajectory_prefix_hasher = scanned['hasher']
        self._trajectory_prefix_bytes = scanned['bytes']
        self._trajectory_prefix_last_step = scanned['last_step']
        self._trajectory_prefix_stat = self._trajectory_stat_identity()

    def _scan_trajectory_prefix(self, checkpoint_step, *, path=None, expected_signature=None):
        """Scan one committed prefix without retaining its trajectory bytes."""
        path = path or self.trajectory_file
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            return {
                'hasher': hashlib.sha256(), 'bytes': 0, 'sha256':
                hashlib.sha256(b'').hexdigest(), 'records': 0,
                'last_step': None, 'removed_records': 0,
                'partial_bytes': 0,
            }
        with open(path, 'rb') as stream:
            if stream.read(8) != NAMD_TRAJECTORY_MAGIC:
                raise ValueError('restart trajectory is not an OpenQP dense TRJ')
            size_bytes = stream.read(8)
            if len(size_bytes) != 8:
                raise ValueError('restart trajectory has a truncated header')
            header_size = struct.unpack('<Q', size_bytes)[0]
            header_bytes = stream.read(header_size)
            if len(header_bytes) != header_size:
                raise ValueError('restart trajectory has a truncated header')
            try:
                header = json.loads(header_bytes.decode('utf-8'))
            except (UnicodeDecodeError, json.JSONDecodeError) as error:
                raise ValueError(
                    'restart trajectory has an invalid header') from error
        offset = 16 + header_size
        if int(header.get('schema_version', -1)) != NAMD_TRAJECTORY_SCHEMA_VERSION:
            raise ValueError('unsupported OpenQP NAMD trajectory schema')
        if header.get('signature') != (expected_signature or self._restart_signature()):
            raise ValueError('restart trajectory and checkpoint model mismatch')
        try:
            dtype = _namd_trajectory_dtype(
                int(header['nstate']), int(header['natom']),
                int(header.get('ncv', 0)))
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(
                'restart trajectory has invalid dimensions') from error
        if int(header.get('record_bytes', -1)) != dtype.itemsize:
            raise ValueError('restart trajectory record layout mismatch')
        payload_size = os.path.getsize(path) - offset
        if payload_size < 0:
            raise ValueError('restart trajectory has an invalid header size')
        partial_bytes = payload_size % dtype.itemsize
        count = payload_size // dtype.itemsize
        if count:
            records = np.memmap(
                path, dtype=dtype, mode='r', offset=offset, shape=(count,))
            steps = np.array(records['step'], copy=True)
            del records
            if (np.any(steps < 0)
                    or (len(steps) > 1 and np.any(np.diff(steps) <= 0))):
                raise ValueError(
                    'restart trajectory contains non-monotonic step records')
            keep = int(np.searchsorted(
                steps, int(checkpoint_step), side='right'))
            last_step = int(steps[keep - 1]) if keep else None
        else:
            keep = 0
            last_step = None
        prefix_bytes = offset + keep*dtype.itemsize
        digest = hashlib.sha256()
        remaining = prefix_bytes
        with open(path, 'rb') as stream:
            while remaining:
                block = stream.read(min(1024 * 1024, remaining))
                if not block:
                    raise ValueError(
                        'restart trajectory changed while being validated')
                digest.update(block)
                remaining -= len(block)
        return {
            'hasher': digest, 'bytes': prefix_bytes,
            'sha256': digest.hexdigest(), 'records': keep,
            'last_step': last_step, 'removed_records': count - keep,
            'partial_bytes': partial_bytes,
        }

    def _trajectory_checkpoint_identity(self, checkpoint_step):
        """Return a checkpoint prefix digest in O(1) during normal MD."""
        cached = getattr(self, '_trajectory_prefix_hasher', None)
        if cached is not None:
            self._require_unchanged_trajectory_prefix()
            last_step = getattr(self, '_trajectory_prefix_last_step', None)
            if last_step is not None and last_step <= int(checkpoint_step):
                return {
                    'bytes': int(self._trajectory_prefix_bytes),
                    'sha256': cached.hexdigest(),
                }

        scanned = self._scan_trajectory_prefix(checkpoint_step)
        if scanned['records'] == 0:
            raise RuntimeError(
                'refusing to checkpoint without the committed dense '
                'trajectory record')
        if not scanned['removed_records'] and not scanned['partial_bytes']:
            self._remember_trajectory_prefix(scanned)
        return {'bytes': scanned['bytes'], 'sha256': scanned['sha256']}

    def _reconcile_trajectory_with_restart(self, checkpoint_step,
                                           expected_prefix):
        """Verify the committed TRJ prefix, then discard later bytes."""
        return self._run_io_collective(
            lambda: self._reconcile_trajectory_on_io_rank(
                checkpoint_step, expected_prefix))

    def _reconcile_trajectory_on_io_rank(self, checkpoint_step,
                                         expected_prefix):
        """Perform packed-trajectory reconciliation on rank zero."""
        scanned = self._scan_trajectory_prefix(checkpoint_step)
        if scanned['records'] == 0:
            raise ValueError(
                'restart checkpoint requires its committed dense trajectory')
        observed = {
            'bytes': scanned['bytes'],
            'sha256': scanned['sha256'],
        }
        if observed != expected_prefix:
            raise ValueError(
                'restart dense trajectory does not match the checkpoint '
                'committed prefix')
        if scanned['removed_records'] or scanned['partial_bytes']:
            with open(self.trajectory_file, 'r+b') as stream:
                stream.truncate(scanned['bytes'])
                stream.flush()
                os.fsync(stream.fileno())
            dump_log(
                self.mol,
                title=(f'NAMD restart removed '
                       f'{scanned["removed_records"]} uncommitted '
                       f'trajectory record(s) and '
                       f'{scanned["partial_bytes"]} incomplete '
                       f'byte(s) after step {checkpoint_step}'),
            )
        self._remember_trajectory_prefix(scanned)

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
        if self.tdc_scheme in (1, 3):
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
    def _hop_triggered_analytic_rescale(self, active, target, istep):
        """Evaluate one exact analytic NAC and rescale a deferred hop.

        The stochastic FSSH decision and coefficient propagation have already
        occurred in the native kernel.  This routine performs no second random
        draw and no second electronic propagation.
        """
        if self.rescale_provider == 'hop_analytic_nac':
            if self._hop_candidate_nac_is_resident(istep):
                # Every pair was already evaluated at this geometry for the
                # TDC/overlap-check consumer; reuse that vector.
                dump_log(self.mol, title=(
                    'NAMD: hop candidate %d -> %d at step %d uses the '
                    'all-pair analytic NAC already evaluated at this step'
                    % (active, target, istep)), section='input')
            else:
                pair = None if self._hop_nac_all_pairs() else (active, target)
                self._update_analytic_nac(istep, compare_overlap=False,
                                          pair=pair)
                dump_log(self.mol, title=(
                    'NAMD: hop candidate %d -> %d at step %d; analytic NAC '
                    'evaluated for %s' % (active, target, istep,
                                          'every state pair (all-pair mode)'
                                          if pair is None else
                                          'the selected pair %d-%d only'
                                          % pair)), section='input')
        direction = np.ascontiguousarray(
            np.asarray(self._last_analytic_dcv, dtype=np.float64)[
                active - 1, target - 1])
        velocity = np.ascontiguousarray(self.vel, dtype=np.float64)
        mass = np.ascontiguousarray(self.mass, dtype=np.float64)
        energies = self._validated_td_energies("OQP::td_energies")
        delta_e = float(energies[target - 1] - energies[active - 1])
        gamma = np.zeros(1, dtype=np.float64)
        discriminant = np.zeros(1, dtype=np.float64)
        status = oqp.oqp_namd_rescale_directional(
            self.natom,
            oqp.ffi.cast("double *", velocity.ctypes.data),
            oqp.ffi.cast("double *", mass.ctypes.data),
            oqp.ffi.cast("double *", direction.ctypes.data),
            delta_e,
            oqp.ffi.cast("double *", gamma.ctypes.data),
            oqp.ffi.cast("double *", discriminant.ctypes.data),
        )
        self._last_rescale_source = 1 if self.rescale_provider == 'analytic_nac' else 2
        self._last_rescale_gamma = float(gamma[0])
        self._last_rescale_discriminant = float(discriminant[0])
        if int(status) != 0:
            if self.frustrated == 'reflect':
                # Frustrated hop: reverse the momentum component along d_IJ
                # (p_a -> p_a - 2 (b/a) d_a with a = sum d_a^2/m_a,
                # b = sum v_a.d_a); the kinetic energy is unchanged.
                avec = float(np.sum(direction**2 / mass[:, None]))
                bvec = float(np.sum(self.vel * direction))
                if avec > 0.0 and np.isfinite(bvec):
                    self.vel = self.vel - 2.0 * (bvec / avec) * direction / mass[:, None]
                    self._frustrated_reflect_count += 1
                    dump_log(self.mol, title=(
                        'NAMD: frustrated hop %d -> %d at step %d (discriminant %.3e); '
                        'velocity component along d_IJ reversed (reflection %d)'
                        % (active, target, istep, float(discriminant[0]),
                           self._frustrated_reflect_count)), section='input')
            return active, False
        self.vel = velocity
        self._last_hop_direction = np.array(direction, copy=True)
        return target, True

    @staticmethod
    def _hop_nac_all_pairs():
        """Return whether hop-candidate NAC evaluation is forced to all pairs.

        ``OQP_NAMD_HOP_NAC_PAIRS=all`` reproduces the pre-2026-09-08
        behaviour (every pair at each candidate hop), which is how the
        published uracil TDC+NAC ensembles were generated.  The default
        ``selected`` evaluates only the active-candidate pair.
        """
        mode = os.environ.get('OQP_NAMD_HOP_NAC_PAIRS', 'selected')
        mode = str(mode).strip().lower()
        if mode not in ('selected', 'all'):
            raise ValueError(
                "OQP_NAMD_HOP_NAC_PAIRS must be 'selected' or 'all'")
        return mode == 'all'

    def _hop_candidate_nac_is_resident(self, istep):
        """Return whether a complete analytic NAC from this step is resident.

        True only when another consumer (``tdc=analytic``, ``rescale=
        analytic_nac`` or ``nacme_check=analytic``) evaluated every pair at
        this same step, so the candidate pair needs no second evaluation.
        """
        return (
            self._needs_analytic_nac()
            and self._last_analytic_dcv is not None
            and self._last_analytic_pair is None
            and istep is not None
            and self._last_analytic_step == int(istep)
        )

    def _clear_hop_triggered_analytic_record(self):
        """Clear exact-NAC fields before a step without a known hop candidate."""
        if self.rescale_provider != 'hop_analytic_nac':
            return
        if self._needs_analytic_nac():
            # A continuous consumer (tdc=analytic, rescale history or
            # nacme_check=analytic) re-evaluates every pair at each step and
            # owns the centered/audit history; clearing it here would leave
            # the requested analytic NACME check permanently unevaluated.
            return
        self._last_analytic_dcv = None
        self._last_analytic_pair = None
        self._last_analytic_step = None
        self._last_analytic_tdc = None
        self._analytic_tdc_previous = None
        self._analytic_tdc_centered = None
        self._nacme_reference_tdc = None
        self._nacme_reference_mask = None
        self._nacme_reference_source = 0

    def _hop(self, allow_hop=True, istep=None):
        """Propagate amplitudes in Fortran and optionally permit a state change."""
        mol = self.mol
        n = self.nstate
        active_before = self.active

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
        params[_P_RAND] = self._hop_random() if allow_hop else 0.0
        params[_P_ACTIVE] = float(self.active)
        params[_P_DECO] = float(self.decoherence)
        params[_P_EDC_C] = self.edc_c
        params[_P_TDC] = float(self.tdc_scheme)
        params[_P_TRIV] = float(self.trivial)
        params[_P_TRIV_THR] = self.trivial_thresh
        params[_P_NSTATE] = float(n)
        params[_P_ALLOW_HOP] = 1.0 if allow_hop else -1.0
        deferred_directional = (
            self.rescale_provider == 'hop_analytic_nac'
            or (self.rescale_provider == 'analytic_nac'
                and self.frustrated == 'reflect'))
        params[_P_RESCALE] = float(2 if deferred_directional else {
            'isotropic': 0, 'analytic_nac': 1,
            'hop_analytic_nac': 2,
        }[self.rescale_provider])
        mol.data["OQP::namd_params"] = params

        # state overlap + time-derivative couplings (FD or NPI), passed to the
        # Fortran hop as flat row-major (n x n) matrices; absolute state
        # energies via namd_eabs. (Same-spin MRSF path.)
        s = canonical_state_overlap(
            np.asarray(mol.data["OQP::td_states_overlap"]).reshape((n, n))
        )
        if self.tdc_provider == 'analytic':
            if self._last_analytic_tdc is None:
                raise RuntimeError(
                    "analytic TDC requested before an analytic NAC was evaluated")
            tdc = np.asarray(self._last_analytic_tdc, dtype=np.float64)
            self._last_tdc_source = 2
        elif self.tdc_provider == 'baeck_an':
            if self._last_baeck_an_tdc is None:
                if self._last_overlap_tdc is None:
                    raise RuntimeError(
                        "Baeck-An TDC requested before overlap warm-up")
                tdc = np.asarray(self._last_overlap_tdc, dtype=np.float64)
                self._last_tdc_source = 1
            else:
                tdc = np.asarray(self._last_baeck_an_tdc, dtype=np.float64)
                self._last_tdc_source = 3
        else:
            tdc = self._compute_tdc(s)
            self._last_tdc_source = self.tdc_scheme
        mol.data["OQP::namd_tdc"] = tdc.reshape(-1).copy()
        mol.data["OQP::namd_stas"] = s.reshape(-1).copy()
        mol.data["OQP::namd_eabs"] = self._validated_td_energies(
            "OQP::td_energies").copy()
        if self.rescale_provider == 'analytic_nac':
            if self._last_analytic_dcv is None:
                raise RuntimeError(
                    "analytic NAC rescaling requested before d_ij was evaluated")
            dcv = np.asarray(self._last_analytic_dcv, dtype=np.float64)
        else:
            dcv = np.zeros((n, n, self.natom, 3), dtype=np.float64)
        mol.data["OQP::namd_dcv"] = dcv.reshape(-1).copy()

        oqp.mrsf_namd_hop(mol)

        # read back
        coef_io = np.array(mol.data["OQP::namd_coef"]).reshape(-1)
        self.coef = coef_io[0::2] + 1j * coef_io[1::2]
        self.vel = np.array(mol.data["OQP::namd_velocity"]).reshape((self.natom, 3)).copy()
        params = np.array(mol.data["OQP::namd_params"])
        results = np.asarray(mol.data["OQP::namd_results"], dtype=float).reshape(-1)
        new_active = int(round(params[_P_ACTIVE]))
        hopped = int(round(params[_P_HOPPED])) == 1
        self._last_rescale_source = int(round(results[n*n + 5]))
        self._last_rescale_gamma = float(results[n*n + 6])
        self._last_rescale_discriminant = float(results[n*n + 7])
        self._last_hop_direction = np.zeros((self.natom, 3), dtype=float)
        target = int(round(results[n*n + 1]))
        blocked = int(round(results[n*n + 2])) == 1
        # In deferred mode the native kernel commits only a trivial-crossing
        # relabel, so params[_P_ACTIVE] is the post-relabel source state.  A
        # stochastic candidate is target /= that state; comparing with the
        # pre-relabel state would turn a relabel into a spurious hop.
        source = new_active
        if (allow_hop and deferred_directional
                and target != source and not blocked):
            new_active, hopped = self._hop_triggered_analytic_rescale(
                source, target, istep)
            blocked = not hopped
            params[_P_ACTIVE] = float(new_active)
            params[_P_HOPPED] = 1.0 if hopped else 0.0
            params[_P_TARGET] = float(target)
            results[n*n] = 1.0 if hopped else 0.0
            results[n*n + 2] = 1.0 if blocked else 0.0
            results[n*n + 5] = 2.0
            results[n*n + 6] = self._last_rescale_gamma
            results[n*n + 7] = self._last_rescale_discriminant
            mol.data["OQP::namd_params"] = params
            mol.data["OQP::namd_results"] = results
            mol.data["OQP::namd_velocity"] = self.vel.reshape(-1).copy()
        # Native mode 1 rescales along d(source, target) with the post-relabel
        # source; that state is not recoverable after a relabel plus a hop, so
        # only record the direction when no trivial relabel occurred.
        relabeled = int(round(results[n*n + 4])) == 1
        if (hopped and self.rescale_provider == 'analytic_nac' and not relabeled
                and 1 <= active_before <= n and 1 <= new_active <= n):
            self._last_hop_direction = np.array(
                dcv[active_before - 1, new_active - 1], copy=True)
        return new_active, hopped

    # ------------------------------------------------------------------ #
    # one nuclear (sub)step: electronic structure, active force, Verlet kick
    # ------------------------------------------------------------------ #
    def _advance_electronic_and_kick(self, istep, r, vel, accel, dt,
                                     with_overlap, continuation):
        """Electronic structure at the installed geometry ``r``, active-state
        force and the velocity-Verlet velocity update over ``dt``.  Returns
        (vel, accel_new, fused) where ``fused`` tells the caller that the
        fused gradient/NAC path was used (its velocity contraction is done by
        the caller after the kick)."""
        mol = self.mol
        if getattr(self, '_nacme_reference_source', 0) == 127:
            self._nacme_reference_source = 0
        self._electronic(with_overlap=with_overlap, continuation=continuation)
        # HT-NAC is evaluated only after the native FSSH kernel selects a
        # stochastic candidate.  Clear the preceding candidate's exact
        # vector so a no-candidate step cannot publish stale NAC data.
        self._clear_hop_triggered_analytic_record()
        fused = self._gradient_nac_fusion_enabled() and with_overlap
        if fused:
            # Fix root identity and phase before constructing the fused
            # gradient/NAC right-hand sides. The analytic velocity
            # contraction is deferred until after the Verlet half kick.
            self._state_overlap(istep, update_analytic=False)
        restraint_force, _ = self._evaluate_conservative_restraints(
            r, self.mass)
        accel_new = (-self._active_gradient() + restraint_force) / self.mass[:, None]
        vel = vel + 0.5 * (accel + accel_new) * dt
        return vel, accel_new, fused

    # ------------------------------------------------------------------ #
    # main loop
    # ------------------------------------------------------------------ #
    def run(self):
        mol = self.mol
        dump_log(mol, title='PyOQP: Tully FSSH Nonadiabatic Molecular Dynamics')
        self._log_rescale_resolution()
        self._prepare_md_outputs()
        restart = self._load_restart()
        if restart is None:
            r = mol.get_system().reshape((self.natom, 3))   # bohr
            # initial electronic structure + force on the active state
            self._electronic(with_overlap=False)
            restraint_force, _ = self._evaluate_conservative_restraints(
                r, self.mass)
            accel = (-self._active_gradient() + restraint_force) / self.mass[:, None]
            if self._needs_analytic_nac():
                self._update_analytic_nac(0, compare_overlap=False)
            self._etot_prev = (
                0.5*np.sum(self.mass[:, None]*self.vel**2)
                + float(np.asarray(mol.energies)[self.active])
                + self._conservative_restraint_energy)
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
            # phase point at the start of the step (for an energy-guarded retry)
            r_start = np.array(r, copy=True)
            vel_start = np.array(self.vel, copy=True)
            accel_start = np.array(accel, copy=True)
            retry_state = self._energy_retry_state()

            # velocity-Verlet position update
            r = r + self.vel * self.dt + 0.5 * accel * self.dt ** 2
            mol.update_system(r.reshape(-1))

            # electronic structure at the new geometry (with overlap vs
            # previous), active-state force and velocity update
            self.vel, accel_new, fused_gradient_nac = (
                self._advance_electronic_and_kick(
                    istep, r, self.vel, accel, self.dt, True, False))

            # state overlap (couplings) and FSSH hop
            if fused_gradient_nac:
                self._update_analytic_nac(istep, compare_overlap=True)
            else:
                self._state_overlap(istep)

            # Energy-guarded substepping: repeat the step from the stored
            # phase point and electronic state with finer nuclear substeps
            # when the total energy jumped by more than disc_tol.
            refinement_attempted = False
            recovery_tol = self.disc_tol
            if self.nve_gate != 'off':
                recovery_tol = min(recovery_tol, self.nve_gate_step_tol)
            if self.disc_substeps > 1 and self._etot_prev is not None:
                odp0 = self._evaluate_odp(r)
                bias0 = 0.0 if odp0 is None else odp0['energy']
                epot0 = (float(np.asarray(mol.energies)[self.active])
                         + bias0 + self._conservative_restraint_energy)
                jump0 = (epot0 + 0.5*np.sum(self.mass[:, None]*self.vel**2)
                         - self._etot_prev)
                if abs(jump0) > recovery_tol and self.prev_data is not None:
                    for nsub in self._energy_refinement_counts(self.disc_substeps):
                        refinement_attempted = True
                        dts = self.dt / nsub
                        self._restore_energy_retry_state(retry_state)
                        mol.put_data(copy.deepcopy(self.prev_data))
                        r = np.array(r_start, copy=True)
                        self.vel = np.array(vel_start, copy=True)
                        accel = np.array(accel_start, copy=True)
                        for k in range(1, nsub + 1):
                            r = r + self.vel * dts + 0.5 * accel * dts ** 2
                            mol.update_system(r.reshape(-1))
                            last = (k == nsub)
                            self.vel, accel_new, fused_gradient_nac = (
                                self._advance_electronic_and_kick(
                                    istep, r, self.vel, accel, dts, last, True))
                            if not last:
                                accel = accel_new
                        if fused_gradient_nac:
                            self._update_analytic_nac(istep, compare_overlap=True)
                        else:
                            self._state_overlap(istep)
                        odp1 = self._evaluate_odp(r)
                        bias1 = 0.0 if odp1 is None else odp1['energy']
                        epot1 = (float(np.asarray(mol.energies)[self.active])
                                 + bias1 + self._conservative_restraint_energy)
                        jump1 = (epot1 + 0.5*np.sum(self.mass[:, None]*self.vel**2)
                                 - self._etot_prev)
                        dump_log(mol, title=(
                            'NAMD energy recovery: step %d; %d nuclear substeps, '
                            'dt %.8f fs; initial change %+.8e Hartree; '
                            'remaining change %+.8e Hartree; criterion %.8e Hartree; '
                            '%s; SCF and active-state force recomputed at each substep'
                            % (istep, nsub, dts/FS_TO_AU, jump0, jump1, recovery_tol,
                               'accepted without numerical rescaling' if abs(jump1) <= recovery_tol
                               else 'criterion not satisfied')), section='input')
                        if abs(jump1) <= recovery_tol:
                            break
                    self._disc_substep_events += 1
                    dump_log(mol, title=('NAMD: total-energy jump %+.4f Hartree at step %d '
                                         'exceeded disc_tol; step repeated with %d '
                                         'substeps of %.3f fs -> residual %+.4f Hartree '
                                         '(substep events %d)'
                                         % (jump0, istep, nsub, dts/FS_TO_AU,
                                            jump1, self._disc_substep_events)),
                             section='input')
            self._last_rescale_source = {
                'isotropic': 0, 'analytic_nac': 1,
                'hop_analytic_nac': 2,
            }[self.rescale_provider]
            self._last_rescale_gamma = np.nan
            self._last_rescale_discriminant = np.nan
            self._last_hop_direction = np.zeros((self.natom, 3), dtype=float)
            active_old = self.active
            odp = self._evaluate_odp(r)
            bias_energy = 0.0 if odp is None else odp['energy']
            # Numerical correction is permitted only after finer nuclear steps
            # fail. It is separate from momentum adjustment at a surface hop.
            self._ref_switch_jump = np.nan
            self._step_numerical_correction = 0.0
            if self._etot_prev is not None:
                epot_now = (float(np.asarray(mol.energies)[self.active])
                            + bias_energy + self._conservative_restraint_energy)
                ke_now = 0.5*np.sum(self.mass[:, None]*self.vel**2)
                jump = (epot_now + ke_now) - self._etot_prev
                somo_case = self._somo_switch_step and self.ref_switch_rescale
                disc_case = self.disc_rescale and abs(jump) > recovery_tol
                if (somo_case or disc_case) and abs(jump) > recovery_tol and not refinement_attempted:
                    dump_log(mol, title=('NAMD numerical energy rescaling withheld at step %d: '
                             'no finer-step retry was completed; energy change %+.8e Hartree'
                             % (istep, jump)), section='input')
                if ((somo_case or disc_case) and refinement_attempted
                        and abs(jump) > recovery_tol):
                    self._require_converged_reference()
                    dump_log(mol, title=('NAMD last-resort numerical energy rescaling: '
                             'step %d; configured finer-step attempts exhausted; '
                             'SCF converged; remaining total-energy change %+.8e Hartree; '
                             'this correction is distinct from a physical surface hop'
                             % (istep, jump)), section='input')
                    ke_target = self._etot_prev - epot_now
                    self._ref_switch_jump = jump
                    if self._somo_switch_step:
                        kind = 'reference switch'
                    elif self._window_leak_step:
                        kind = 'window-leak discontinuity'
                    else:
                        kind = 'energy discontinuity'
                    if (np.isfinite(ke_target) and np.isfinite(ke_now)
                            and ke_target > 0.0 and ke_now > 0.0):
                        if not self._somo_switch_step:
                            self._disc_event_count += 1
                        factor = np.sqrt(ke_target/ke_now)
                        self.vel = self.vel*factor
                        self._scale_analytic_velocity_contractions(factor, istep)
                        dump_log(mol, title=('NAMD numerical correction details: kinetic energy '
                                 '%.12f -> %.12f Hartree; velocity factor %.12f; '
                                 'energy correction %+.8e Hartree'
                                 % (ke_now, ke_target, factor, -jump)), section='input')
                        dump_log(mol, title=(
                            'NAMD numerical correction: relative kinetic-energy change '
                            '%+.6f%%; next interval retries the configured dt %.8f fs; '
                            'energy criteria and finer-step recovery remain active'
                            % (100.0*(ke_target-ke_now)/ke_now, self.dt/FS_TO_AU)),
                            section='input')
                        self._disc_energy_absorbed += jump
                        self._step_numerical_correction = float(jump)
                        dump_log(mol, title=('NAMD: %s at step %d; '
                                             'active-state energy jump %+.4f Hartree '
                                             'absorbed by isotropic velocity rescaling '
                                             '(events: switch %d, other %d; absorbed %+.4f Ha)'
                                             % (kind, istep, jump,
                                                self._somo_switch_count,
                                                self._disc_event_count,
                                                self._disc_energy_absorbed)),
                                 section='input')
                    else:
                        dump_log(mol, title=('NAMD: %s at step %d; '
                                             'energy jump %+.4f Hartree exceeds the '
                                             'kinetic energy, velocities unchanged'
                                             % (kind, istep, jump)),
                                 section='input')
            energy_before_transition = (
                0.5*np.sum(self.mass[:, None]*self.vel**2)
                + float(np.asarray(mol.energies)[active_old])
                + bias_energy
                + self._conservative_restraint_energy
            )
            hop_ready = self._prepare_hop_step(istep)
            if getattr(self, '_restart_boundary', False):
                # Restart boundary (GAMESS FIRST-step behaviour): the
                # previous-step states belong to another SCF branch, so the
                # overlap-derived coupling is meaningless.  Freeze the
                # electronic coefficients and attempt no hop on this step.
                new_active, hopped = self.active, False
                self._nacme_reference_source = 127
                dump_log(mol, title=('NAMD: restart boundary at step %d; '
                                     'coefficient propagation and hopping '
                                     'skipped for this step' % istep),
                         section='input')
            elif getattr(self, '_pending_nacme_gate_error', None) is not None:
                new_active, hopped = self.active, False
            else:
                new_active, hopped = self._hop(
                    allow_hop=hop_ready, istep=istep)

            active_changed = new_active != active_old
            if active_changed:
                self.active = new_active
                # force for the next step is on the new active surface. This
                # also covers trivial-crossing following, where the Fortran
                # kernel can update ACTIVE without marking HOPPED.
                accel_new = (
                    -self._active_gradient() +
                    self._conservative_restraint_force
                ) / self.mass[:, None]
                energy_after_transition = (
                    0.5*np.sum(self.mass[:, None]*self.vel**2)
                    + float(np.asarray(mol.energies)[self.active])
                    + bias_energy
                    + self._conservative_restraint_energy
                )
                transition_energy_jump = (
                    energy_after_transition - energy_before_transition)
            else:
                transition_energy_jump = np.nan
            # The numerical correction was applied before energy_before_transition.
            # Do not test that correction against the physical-hop tolerance.
            # Its uncorrected energy change and velocity factor are logged above;
            # the ordinary step and cumulative NVE criteria still apply below.

            self._apply_thermostat(istep)
            accel = accel_new
            self._etot_prev = (
                0.5*np.sum(self.mass[:, None]*self.vel**2)
                + float(np.asarray(mol.energies)[self.active])
                + bias_energy + self._conservative_restraint_energy)
            self._record_previous(r)
            self._log_step(
                istep, r, hopped=hopped,
                transition_energy_jump=transition_energy_jump)
            self._save_restart(istep, r, self.vel, accel)

        # State-overlap evaluation normally refreshes the generic molecule
        # JSON after the force calculation.  In the fused gradient/NAC path it
        # must run before the gradient so that root phases are fixed for both
        # right-hand sides.  Save once at trajectory completion to ensure the
        # generic JSON carries the final current-geometry gradient as well as
        # the restart/trajectory records.
        if (mol.config['guess']['save_mol']
                or self._gradient_nac_fusion_enabled()):
            mol.save_data()
        dump_log(mol, title=('PyOQP: NAMD trajectory complete '
                             '(reference switches %d, window leaks %d, other energy '
                             'discontinuities %d, absorbed %+.4f Ha, frustrated-hop '
                             'reflections %d, SCF fallbacks %d)'
                             % (self._somo_switch_count, self._window_leak_count,
                                self._disc_event_count, self._disc_energy_absorbed,
                                self._frustrated_reflect_count,
                                self._scf_fallback_steps)))

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
        electronic_epot = self._unbiased_potential_energy
        epot = (electronic_epot
                + (0.0 if odp is None else odp['energy'])
                + getattr(self, '_conservative_restraint_energy', 0.0))
        pops = np.abs(self.coef) ** 2
        self._update_nve_gate(istep, epot, ekin, transition_energy_jump)
        dump_log(
            mol,
            title=(f'NAMD step {istep:6d}  t={(self._physical_time_fs(istep)):9.3f} fs  '
                   f'active={self.active}  E_tot={ekin+epot:.8f}  '
                   f'E_pot={epot:.8f}  E_elec={electronic_epot:.8f}  '
                   f'U_ODP={(0.0 if odp is None else odp["energy"]):.8f}  '
                   f'E_drop={getattr(self, "_droplet_energy", 0.0):.8f}  '
                   f'drop_max={getattr(self, "_droplet_max_penetration", 0.0):.5f}bohr  '
                   f'drop_n={getattr(self, "_droplet_active_count", 0)}  '
                   f'drop_fmax={getattr(self, "_droplet_force_max", 0.0):.3e}Ha/bohr  '
                   f'E_com={getattr(self, "_solute_com_energy", 0.0):.8f}  '
                   f'E_kin={ekin:.8f}  '
                   f'dQ_therm={getattr(self, "_thermostat_exchange", 0.0):+.3e}  '
                   f'Q_therm={getattr(self, "_thermostat_exchange_cumulative", 0.0):+.3e}  '
                   f'hop={hopped}  {self._hop_rng_log()}  '
                   f'pop={np.array2string(pops, precision=4)}'),
        )
        self._write_md_trajectory(istep, r, epot, ekin, hopped)
        self._enforce_nacme_gate()
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



def _image_field_stagnant(iteration, delta, energies, driver):
    """True when a reference-density image field whose ESPF charges only
    fluctuate at the SCF noise floor can be accepted.

    ``iteration`` is the zero-based image iteration, ``delta`` its max |dq|
    (e) and ``energies`` the reference SCF energy of every iteration so far.
    The strict ``IMAGE_TOL`` test runs first; this is the fallback for an ROHF
    reference whose soft orbital rotations amplify the SCF gradient noise into
    charge noise above ``IMAGE_TOL`` (see the constants on the QM/MM driver).
    """
    if iteration + 1 < int(driver.IMAGE_STAGNANT_MINITER) or len(energies) < 3:
        return False
    last = energies[-3:]
    return delta < driver.IMAGE_TOL_STAGNANT and max(last) - min(last) < driver.IMAGE_ETOL

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
        if (self._needs_analytic_nac()
                or self.rescale_provider == 'hop_analytic_nac'):
            raise NotImplementedError(
                "analytic NAC TDC/rescaling/check is not yet available for QM/MM NAMD")
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
        q = mol.config['qmmm']
        pdb_file = self._resolve_qmmm_aux_file(q['pdb_file'])
        ff_files = [self._resolve_qmmm_aux_file(s) for s in
                    str(q['forcefield_files']).replace(',', ' ').split() if s]
        self._qmmm_pdb_file = pdb_file
        self._qmmm_forcefield_files = ff_files
        self._qmmm_restart_identity_cache = None
        # Topology (ascending) order: the embedded driver sorts its own copy and
        # numbers link-atom host rows in that order, and the QM molecule built
        # from the PDB is in topology order, so every per-row array here (QM
        # geometry, gradients, charges, hop velocities) must use the same order.
        self.qm_atoms = np.array(sorted(_parse_int_list(q['qm_atoms'])), dtype=int)
        self.cutoff = _resolve_cutoff(str(q['cutoff']).strip())   # NoCutoff | PME | Ewald | ...
        from oqp.library.qmmm_driver import is_periodic_method
        self.periodic = is_periodic_method(self.cutoff)   # PME / Ewald / CutoffPeriodic
        _validate_odp_boundary_conditions(self.odp, self.periodic)
        embedding = str(q['embedding']).strip()
        frontier_scheme = str(q.get('frontier_scheme', 'none')).strip()
        _et = q.get('ewald_tol', None)
        ewald_tol = None if _et in (None, '', 'none', 'None') else float(_et)
        lj_switch = str(q.get('lj_switch', 'false')).strip().lower() in ('1', 'true', 'yes', 'on')
        h_lj = str(q.get('h_lj', 'false')).strip().lower() in ('1', 'true', 'yes', 'on')
        _w = q.get('mm_charge_width', None)
        mm_charge_width = None if _w in (None, '', 'none', 'None', 0, 0.0, '0') else float(_w)

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
            ewald_tol=ewald_tol,
            lj_switch=lj_switch,
            h_lj=h_lj,
            mm_charge_width=mm_charge_width,
        )
        self.mm = self.driver.mm_systems

        # Covalent QM/MM boundary: the QM Molecule must carry the real QM atoms
        # followed by one hydrogen link atom per cut bond, in the order the
        # driver detects them (sorted by (QM host, MM host) topology index).
        # ``[input] system = file.pdb <1-based QM indices>`` builds exactly that.
        self.link_atoms = list(self.driver.link_atoms)
        self.nqm = int(len(self.qm_atoms))
        self._validate_qm_molecule_layout()

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
        self._resolve_active_atoms(q)
        # rigid-water (SHAKE/RATTLE) constraints for the MM region.  Built
        # before the velocities and before the identities: constraint closure
        # can move an atom either into or out of the active set, and both need
        # the final set.  An atom that closure ACTIVATES cannot have its random
        # velocity restored afterwards, so the Maxwell draw has to come second.
        self._build_constraints()
        self._restart_system_identity = self._qmmm_restart_system_identity(
            sys0, q)
        self._wham_system_identity = self._qmmm_wham_system_identity(sys0, q)

        # full-system Maxwell-Boltzmann velocities (a.u.), COM removed
        if self.restart_requested:
            self.v_all = np.zeros((self.natom_all, 3))
        else:
            sig = np.sqrt(KB_HARTREE * self.init_temp / self.m_all)
            self.v_all = self._counter_normals((self.natom_all, 3)) * sig[:, None]
            p = (self.m_all[:, None] * self.v_all).sum(axis=0)
            self.v_all -= p / self.m_all.sum()
            # A held atom is not thermalised and carries no momentum.  With
            # nothing held the two lines above are the whole draw, unchanged.
            self._hold_velocities()
            self._remove_moving_com()

        # sync the QM Molecule geometry from the pdb QM atoms
        self._sync_positions()
        # QM-region masses for the hop (already set by super from mol.get_mass())
        self.qm_mass = self.mass.copy()
        self._setup_qmmm_restraint_targets()

    # ------------------------------------------------------------------ #
    # [qmmm] keys added with the periodic/embedding controls.  They enter the
    # restart and WHAM identities only when set to a non-default value, so
    # checkpoints written before these keys existed keep validating.
    _QMMM_OPTIONAL_HAMILTONIAN_KEYS = ('ewald_tol', 'lj_switch', 'h_lj', 'mm_charge_width')
    # The active/frozen selection decides which atoms a restart propagates, so a
    # checkpoint is bound to it -- but only once one is actually requested, so
    # checkpoints written before these keys existed keep validating.
    _QMMM_SELECTION_KEYS = ('active_atoms', 'frozen_atoms', 'active_radius', 'active_from_pdb')

    @classmethod
    def _qmmm_identity_config(cls, qmmm_config):
        """Copy of ``qmmm_config`` with the optional Hamiltonian keys dropped
        when at their default and normalised (bool / float) otherwise."""
        cfg = dict(qmmm_config)
        for key in cls._QMMM_OPTIONAL_HAMILTONIAN_KEYS:
            if key not in cfg:
                continue
            value = cfg.pop(key)
            if key in ('lj_switch', 'h_lj'):
                on = value if isinstance(value, bool) else (
                    str(value).strip().lower() in ('1', 'true', 'yes', 'on', 't'))
                if on:
                    cfg[key] = True
            else:
                text = '' if value is None else str(value).strip()
                if text.lower() in ('', 'none'):
                    continue
                try:
                    number = float(text)
                except ValueError:
                    cfg[key] = text            # let the run's own validation complain
                    continue
                if number == 0.0:              # '0', '0.0', '0.00', '0e0': point charges
                    continue
                cfg[key] = number
        for key in cls._QMMM_SELECTION_KEYS:
            if key not in cfg:
                continue
            value = cfg.pop(key)
            text = '' if value is None else str(value).strip()
            if key in ('active_atoms', 'frozen_atoms'):
                # '0' is atom zero, a real selection -- only a blank is 'unset'
                if text.lower() in ('', 'none'):
                    continue
            elif key == 'active_from_pdb':
                if text.lower() in ('', 'none', 'false', 'off', 'no', '0'):
                    continue
            else:                              # active_radius: 0 A selects nothing
                if text.lower() in ('', 'none'):
                    continue
                try:
                    if float(text) == 0.0:
                        continue
                except ValueError:
                    pass
            cfg[key] = text
        return cfg

    def _qmmm_restart_system_identity(self, system, qmmm_config):
        """Bind QM/MM restarts to atoms, topology, selection, and force field."""
        qmmm_config = self._qmmm_identity_config(qmmm_config)
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
        array_parts = [
            ('atomic_numbers', atomic_numbers, '<i8'),
            ('masses_electron', self.m_all, '<f8'),
            ('initial_coordinates_bohr', self.r_all, '<f8'),
            ('qm_atoms', self.qm_atoms, '<i8'),
            ('bonds', bonds, '<i8'),
        ]
        if not self._all_atoms_move:
            # Bind the checkpoint to the atoms this run actually propagates, so
            # a restart cannot silently resume with a different set.  The
            # resolved set is used rather than the config text because
            # active_from_pdb keeps its selection in the PDB's B-factor column,
            # which the config does not capture.  Added only when a selection
            # exists, so checkpoints written before these keys existed keep
            # validating.
            array_parts.append(('active_atoms_resolved', self.active_atoms, '<i8'))
        digest = _restart_identity_digest(
            array_parts=tuple(array_parts),
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

    def _qmmm_wham_system_identity(self, system, qmmm_config):
        """Hash QM/MM topology and Hamiltonian without initial coordinates."""
        qmmm_config = self._qmmm_identity_config(qmmm_config)
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
        hamiltonian_options = {
            key: qmmm_config.get(key)
            for key in ('embedding', 'frontier_scheme', 'cutoff',
                        'nonbondedmethod', 'ewald_tol', 'lj_switch', 'h_lj',
                        'mm_charge_width')
            if key in qmmm_config
        }
        digest = _restart_identity_digest(
            array_parts=(
                ('atomic_numbers', atomic_numbers, '<i8'),
                ('masses_electron', self.m_all, '<f8'),
                ('qm_atoms', self.qm_atoms, '<i8'),
                ('bonds', bonds, '<i8'),
            ),
            text_parts=(
                ('atom_metadata', json.dumps(
                    atom_metadata, separators=(',', ':'))),
                ('qmmm_hamiltonian', json.dumps(
                    hamiltonian_options, sort_keys=True, separators=(',', ':'),
                    default=str)),
                ('openmm_system', self._mm.XmlSerializer.serialize(system)),
            ),
        )
        return {
            'kind': 'qmmm', 'natom': int(self.natom_all),
            'nqm': int(len(self.qm_atoms)), 'sha256': digest,
        }

    # ------------------------------------------------------------------ #
    def _resolve_active_atoms(self, qmmm_config):
        """``[qmmm] active_atoms`` / ``frozen_atoms`` / ``active_radius`` /
        ``active_from_pdb``: which atoms this trajectory propagates.  With no
        selection every atom moves, which is what every deck written before
        these keys existed expects.

        Freezing changes what moves, never the physics: a held atom keeps its
        charge, its embedding field and its force contribution.  It is simply
        not integrated, so no work is done on it and energy conservation still
        holds for the atoms that do move.  The QM region is always propagated.
        """
        u = self._u
        box = self.driver._box_lengths_bohr()
        self.active_atoms, self.frozen_atoms = resolve_active_set(
            qmmm_config,
            self.pdb.topology,
            np.array(self.pdb.positions.value_in_unit(u.angstrom)),
            self.qm_atoms,
            box_ang=None if box is None else np.asarray(box) * BOHR_TO_NM * 10.0,
            default_all=True,          # dynamics propagates everything unless asked
            pdb_path=getattr(self, '_qmmm_pdb_file', None),
        )
        self._set_move_mask(self.active_atoms)
        if not self._all_atoms_move:
            print(f"[QM/MM NAMD] active atoms: {len(self.active_atoms)} of {self.natom_all} "
                  f"propagated, {self.natom_all - len(self.active_atoms)} held fixed "
                  f"(their charges and forces still act)")

    def _hold_velocities(self):
        """Zero the velocity of every atom this run does not propagate.

        A held atom that kept a velocity would never move (the position update
        is masked) but its fictitious kinetic energy would still enter the
        reported temperature, the NVE audit and the checkpoint.
        """
        if self._all_atoms_move:
            return
        self.v_all *= self._move_mask

    def _remove_moving_com(self):
        """Put the centre of mass of the propagated subsystem at rest.

        Initialisation only.  While the trajectory runs, the moving atoms can
        legitimately pick up net momentum from the held environment, and that
        translation is a real degree of freedom: removing it at every step
        would silently constrain the dynamics and charge the removed kinetic
        energy to the thermostat's reported exchange.
        """
        if self._all_atoms_move:
            return
        m_move = self.m_all * self._move_mask[:, 0]
        if m_move.sum() > 0.0:
            p = (m_move[:, None] * self.v_all).sum(axis=0)
            self.v_all -= (p / m_move.sum()) * self._move_mask

    def _set_move_mask(self, active):
        """Column mask (natom_all, 1): 0 on a held atom, 1 everywhere else.

        The mask zeroes only atoms that are genuinely held.  A virtual site is
        never "active" -- OpenMM places it from its parents -- so building the
        mask from membership in ``active`` would mark every virtual site as
        held, and a TIP4P topology that selected nothing would stop looking
        like an unselected run.  Taking the complement instead keeps the mask
        all ones exactly when nothing is held.
        """
        held = held_atoms(self.pdb.topology, active)
        mask = np.ones((self.natom_all, 1))
        if held:
            mask[np.asarray(sorted(held), dtype=int), 0] = 0.0
        self._move_mask = mask
        self._held_atoms = held
        self._all_atoms_move = not held

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
        if not self._all_atoms_move:
            # SHAKE/RATTLE cannot satisfy a constraint between a propagated atom
            # and a held one: correcting the bond would drag the held atom.  So
            # such a pair is resolved first (held when either atom was frozen,
            # propagated otherwise), and the constraints left inside the held set
            # are dropped -- those atoms never move, so there is nothing to
            # constrain, and keeping them would let a starting geometry that
            # slightly violates the bond length push a fixed atom.
            active, frozen = freeze_constrained_partners(
                list(zip(ci, cj)), self.active_atoms, self.frozen_atoms)
            self.active_atoms = np.array(sorted(active), dtype=int)
            self.frozen_atoms = frozen
            self._set_move_mask(self.active_atoms)
            # Every atom outside the active set is held, not just the ones
            # frozen_atoms named: an active_atoms / active_radius /
            # active_from_pdb selection holds everything it did not select.
            # Testing `frozen` here would leave the constraints inside those
            # unselected waters in place, and SHAKE would then move them.
            held = held_atoms(self.pdb.topology, self.active_atoms)
            keep = [k for k in range(len(ci)) if ci[k] not in held and cj[k] not in held]
            ci = [ci[k] for k in keep]
            cj = [cj[k] for k in keep]
            cd = [cd[k] for k in keep]
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
        # only the propagated atoms carry kinetic degrees of freedom
        nmove = int(round(float(np.sum(self._move_mask))))
        ndof = 3 * nmove - ncon - 3
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

    def _apply_thermostat(self, istep):
        """Thermostat the full QM/MM system and report constrained dK."""
        self._thermostat_exchange = 0.0
        if self.thermostat == 'off':
            return
        kinetic_before = 0.5*np.sum(self.m_all[:, None]*self.v_all**2)
        self.v_all, _ = self._langevin_update(
            self.v_all, self.m_all, istep)
        # The thermostat draws a new velocity for every row; a held atom must
        # not be given one, or its fictitious kinetic energy would enter the
        # reported temperature and the thermostat's energy exchange.  Only the
        # held rows are zeroed: the moving subsystem's net translation is a
        # physical degree of freedom and must survive the thermostat.
        self._hold_velocities()
        self._rattle(self.r_all, self.v_all)
        kinetic_after = 0.5*np.sum(self.m_all[:, None]*self.v_all**2)
        self._thermostat_exchange = float(kinetic_after - kinetic_before)
        self._thermostat_exchange_cumulative += self._thermostat_exchange
        self.vel = self._qm_velocities()

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
        # QM Molecule coords (bohr): real QM atoms + hydrogen link atoms
        self.mol.update_system(self._qm_positions_bohr().reshape(-1))

    # ------------------------------------------------------------------ #
    # covalent-boundary (link-atom) helpers
    #
    # The QM Molecule has natom = nqm + nlink centres.  Link atoms have no
    # dynamical degrees of freedom: their position is the fixed linear
    # combination r_L = r_QM + g (r_MM - r_QM) of two real atoms, and every
    # force on them is redistributed onto those hosts by the chain rule (see
    # _total_force_espf).  The FSSH hop rescales the REAL QM-atom velocities
    # only, so the velocity vector handed to the kernel carries zeros in the
    # link rows (kinetic energy = that of the real QM atoms); coupling-derivative
    # corrections that need the motion of every QM centre get the kinematic
    # link velocity (1-g) v_QM + g v_MM instead.
    # ------------------------------------------------------------------ #
    def _validate_qm_molecule_layout(self):
        nlink = len(self.link_atoms)
        if self.natom != self.nqm + nlink:
            if nlink:
                raise ValueError(
                    f"NAMD QM/MM across a covalent boundary: the QM/MM partition "
                    f"cuts {nlink} bond(s) but the QM molecule has {self.natom} "
                    f"atoms instead of {self.nqm} QM atoms + {nlink} link "
                    f"hydrogen(s). Build the QM molecule from the PDB "
                    f"('[input] system = file.pdb <1-based QM indices>') so the "
                    f"link atoms are appended automatically.")
            raise ValueError(
                f"NAMD QM/MM: the QM molecule has {self.natom} atoms but "
                f"[qmmm] qm_atoms selects {self.nqm}.")
        z_mol = np.asarray(self.mol.get_atoms2("charge"), dtype=float).reshape(-1)
        z_top = {a.index: (0 if a.element is None else a.element.atomic_number)
                 for a in self.pdb.topology.atoms()}
        z_expected = [z_top[int(i)] for i in self.qm_atoms] + [1] * nlink
        if any(abs(z_mol[k] - z_expected[k]) > 0.5 for k in range(self.natom)):
            raise ValueError(
                "NAMD QM/MM: the QM molecule's atoms do not match [qmmm] "
                "qm_atoms (in topology order) followed by the hydrogen link "
                f"atoms: molecule Z={z_mol.astype(int).tolist()}, expected "
                f"{z_expected}. Note '[input] system = file.pdb ...' indices are "
                "1-based while [qmmm] qm_atoms are 0-based.")

    def _qm_positions_bohr(self):
        """(natom, 3) QM-centre coordinates (bohr): real QM atoms in topology
        order, then the hydrogen link atoms on their cut bonds."""
        box = self.driver._box_lengths_bohr()       # None for a cluster
        if box is None:
            r = self.r_all[self.qm_atoms]
            qm_xyz = None
        else:
            # bonded QM fragments made whole under PBC (minimum-image bonds)
            qm_xyz = self.driver.unwrap_qm(lambda i: self.r_all[i], box)
            r = np.asarray([qm_xyz[int(i)] for i in self.qm_atoms], dtype=float)
        if not self.link_atoms:
            return r
        links = []
        for l in self.link_atoms:
            bond = self.driver._min_image(self.r_all[l.mm_index] - self.r_all[l.qm_index], box)
            host = self.r_all[l.qm_index] if qm_xyz is None else qm_xyz[int(l.qm_index)]
            links.append(host + l.g * bond)
        return np.vstack([r, np.asarray(links, dtype=float)])

    def _qm_velocities(self, kinematic=False):
        """(natom, 3) QM-centre velocities.  Link rows are zero (no dynamical
        DOF; the hop rescales real QM atoms only) unless ``kinematic`` is set,
        in which case they carry (1-g) v_QM + g v_MM."""
        v = np.zeros((self.natom, 3))
        v[:self.nqm] = self.v_all[self.qm_atoms]
        if kinematic:
            for a, l in enumerate(self.link_atoms):
                v[self.nqm + a] = ((1.0 - l.g) * self.v_all[l.qm_index]
                                   + l.g * self.v_all[l.mm_index])
        return v

    def _store_qm_velocities(self, vel):
        """Write the (possibly rescaled) real QM-atom velocities back into the
        full-system velocity array; link rows are discarded."""
        self.v_all[self.qm_atoms] = np.asarray(vel)[:self.nqm]

    def _embedding_field(self):
        """(potmm, potqm) the embedded SCF sees: the MM electrostatic potential
        at every QM centre, or a zero field for [qmmm] embedding=mechanical
        (gas-phase QM Hamiltonian; the QM-MM electrostatics is then left to
        OpenMM with the QM ESP charges, as in OpenQpQMMM.compute_force)."""
        if getattr(self.driver, "Embedding", "") == "mechanical":
            return self.driver._zero_embedding()
        return self.driver.electrostatic_potential()

    @staticmethod
    def _embedded_scf(sp, warm=False):
        """Run the embedded SCF (the ESPF term is already in hcore) through the
        same robustness ladder as a gas-phase reference: primary converger,
        then SOSCF/TRAH escalation warm-started from the current orbitals.
        A reference that still does not converge stops the run: propagating
        on an unconverged SCF gives an inconsistent energy/force pair, and a
        DIIS loop that stops at its iteration limit even leaves the density
        records in an intermediate state.

        ``warm``: the orbitals already held are a converged solution of a
        nearby Hamiltonian (the previous MD step or an image iteration).  With
        a DIIS primary the solve then starts with SOSCF: where the ROHF triplet
        is not in the order of the effective-Fock energies, the first refill of
        the orbitals swaps an occupied and an open orbital (+0.2 Hartree), DIIS
        stalls and hands over to SOSCF anyway, and starting with SOSCF skips the
        stalled DIIS stage.  An explicitly selected primary (soscf, trah, auto,
        ...) runs as requested."""
        saved_converger = sp.converger_type
        if warm and str(saved_converger).lower() == 'diis':
            sp.converger_type = 'soscf'
        try:
            converged = sp._run_scf()
        finally:
            sp.converger_type = saved_converger
        if not converged:
            raise RuntimeError(
                "NAMD QM/MM: the embedded SCF did not converge (primary "
                "converger and the SOSCF/TRAH escalation).  Raise [scf] maxit, "
                "loosen [scf] conv, or check the QM/MM contacts.")

    def _load_restart(self):
        """Restore the checkpoint and, when it carries the converged orbitals
        of the saved step (restored into the molecule with the rest of the
        previous-step data), let the next step warm-start from them exactly as
        an uninterrupted trajectory would.  _start_scf_orbitals still checks
        their size and overlap metric before using them."""
        restart = super()._load_restart()
        if restart is not None:
            data = self.prev_data if isinstance(self.prev_data, dict) else {}
            self._scf_orbitals_ready = "OQP::VEC_MO_A" in data
        return restart

    def _start_scf_orbitals(self, sp):
        """Basis, one-electron integrals and starting orbitals for this MD
        step's embedded SCF.  Once a step has converged, the next step starts
        from its orbitals, re-orthonormalised in the overlap metric of the new
        geometry (symmetric orthonormalisation keeps each orbital's character)
        with the density rebuilt from them.  A fresh guess every step starts
        the SCF ~1.3 Hartree above the solution, costs ~30 iterations, and
        lets the open-shell reference settle on a different solution from
        one step to the next.  Returns True for a warm start."""
        mol = self.mol
        if not getattr(self, "_scf_orbitals_ready", False):
            sp._prep_guess()
            return False
        c_prev = np.array(mol.data["OQP::VEC_MO_A"], dtype=float, copy=True)
        e_prev = np.array(mol.data["OQP::E_MO_A"], dtype=float, copy=True)
        oqp.library.set_basis(mol)
        ints_1e(mol)
        nbf = mol.data.get_basis()["nbf"]
        target = np.asarray(mol.data["OQP::VEC_MO_A"], dtype=float)
        if c_prev.size != nbf * nbf or target.size != nbf * nbf:
            oqp.library.guess(mol)
            return False
        c = c_prev.reshape((nbf, nbf)).T                      # C[ao, mo]
        packed_s = np.asarray(mol.data["OQP::SM"], dtype=float).ravel()
        if packed_s.size != nbf * (nbf + 1) // 2:
            oqp.library.guess(mol)
            return False
        s_ao = np.zeros((nbf, nbf))
        s_ao[np.tril_indices(nbf)] = packed_s             # row-major lower triangle
        s_ao = s_ao + s_ao.T - np.diag(np.diag(s_ao))
        w, v = np.linalg.eigh(c.T @ s_ao @ c)
        if not np.all(np.isfinite(w)) or w.min() <= 1.0e-8:
            oqp.library.guess(mol)
            return False
        c = c @ (v * (1.0 / np.sqrt(w))) @ v.T
        packed = np.ascontiguousarray(c.T.reshape(target.shape))
        mol.data["OQP::VEC_MO_A"][...] = packed
        mol.data["OQP::VEC_MO_B"][...] = packed
        mol.data["OQP::E_MO_A"][...] = e_prev.reshape(np.shape(mol.data["OQP::E_MO_A"]))
        mol.data["OQP::E_MO_B"][...] = e_prev.reshape(np.shape(mol.data["OQP::E_MO_B"]))
        oqp.guess_json(mol)
        dump_log(mol, title="PyOQP: NAMD QM/MM SCF warm start from the previous step's orbitals",
                 section='')
        return True

    def _fold_link_charges(self, pchg):
        """(nqm,) MM-facing QM charges: each link atom's ESPF charge is added
        to its QM host so the total QM charge is conserved when the QM region
        is represented by point charges on the real QM atoms only."""
        pchg = np.asarray(pchg, dtype=float)
        q = pchg[:self.nqm].copy()
        for a, l in enumerate(self.link_atoms):
            q[l.host_row] += pchg[self.nqm + a]
        return q

    def _project_link_rows(self, g):
        """Chain-rule a (natom, 3) QM-centre gradient/force onto the real QM
        atoms: returns (g_qm (nqm,3), mm_host contributions {mm_index: (3,)})."""
        g = np.asarray(g, dtype=float)
        g_qm = g[:self.nqm].copy()
        mm_part = {}
        for a, l in enumerate(self.link_atoms):
            gl = g[self.nqm + a]
            g_qm[l.host_row] += (1.0 - l.g) * gl
            mm_part[l.mm_index] = mm_part.get(l.mm_index, 0.0) + l.g * gl
        return g_qm, mm_part

    def _electronic_qmmm(self, with_overlap):
        """Embedded SCF + MRSF excitation; returns (potmm, potqm)."""
        from oqp.library.qmmm_driver import (
            unpack_lower_tri_single, unpack_lower_tri_multi, pack_lower_tri_single)
        mol = self.mol
        potmm, potqm = self._embedding_field()

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
            if self.driver._ewald() is not None:
                # The periodic QM-image self-consistency below needs the ESPF
                # charge operator of the native path; the tight-binding backend
                # exposes no equivalent, so a periodic DFTB/xTB NAMD would
                # silently drop the image term and its energy correction.
                raise NotImplementedError(
                    "NAMD QM/MM with method=dftb/xtb is implemented for non-"
                    "periodic clusters only ([qmmm] cutoff=NoCutoff); the "
                    "periodic (PME/Ewald) QM-image self-consistency is not "
                    "available for the tight-binding backend.")
            # (embedding=mechanical is rejected above and by the input
            # checker for the tight-binding NAMD path, so potmm is the field.)
            mol.dftb_external_potential = np.asarray(potmm, dtype=float)
            self._e_img, self._f_img = 0.0, None
            sp = SinglePoint(mol)
            ref = sp.reference()
            if with_overlap:
                mol.back_door = (self.prev_xyz, self.prev_data)
                BasisOverlap(mol).overlap()
            sp.excitation(ref)
            LastStep(mol).compute(mol)
            return potmm, potqm

        sp = SinglePoint(mol)
        warm_start = self._start_scf_orbitals(sp)
        nat = mol.data["natom"]
        nbf = mol.data.get_basis()["nbf"]
        # Periodic full-ESPF: the QM charges also interact with their own images
        # (Bonfrate et al. JCTC 2024, eq 8).  The image field psi_img q is added to
        # the MM potential and made self-consistent with the ground-state ESPF
        # charges; the energy carries the double-counting correction
        # -1/2 q psi_img q (see _total_force_espf).  With ground-state charges
        # the force is exact for the ground state and an approximation for
        # excited states (their charges differ slightly from the ground state).
        ewald = self.driver._ewald() if getattr(self.driver, "espf_full", False) else None
        potmm_mm = np.asarray(potmm, dtype=float).copy()
        if ewald is not None:
            psi_img, dpsi_img = ewald.qm_image_matrix(self.driver._qm_center_positions_bohr())
            q_prev = (self._q_img if getattr(self, "_q_img", None) is not None
                      and len(self._q_img) == nat else np.zeros(nat))
        else:
            psi_img = dpsi_img = None
            q_prev = None
        converged = psi_img is None
        delta, it = float("inf"), -1
        e_hist = []
        for it in range(int(self.driver.IMAGE_MAXITER)):
            potmm = potmm_mm if psi_img is None else potmm_mm + psi_img @ q_prev
            mol.data["OQP::POTMM"] = potmm
            mol.data["OQP::POTQM"] = np.zeros((nat, nat))
            oqp.espf_op_corr(mol)
            espf = unpack_lower_tri_multi(mol.data["OQP::ESPF_CORR"], nbf, nat)
            hcore = unpack_lower_tri_single(mol.get_hcore(), nbf)
            hcore += np.einsum("ijk,i->jk", espf, potmm)
            mol.set_hcore(pack_lower_tri_single(hcore))
            self._embedded_scf(sp, warm=(warm_start or it > 0))
            self._scf_orbitals_ready = True
            if psi_img is None:
                break
            oqp.form_esp_charges(mol)
            q_new = np.array(mol.data["OQP::partial_charges"], dtype=float)
            delta = float(np.abs(q_new - q_prev).max())
            e_hist.append(float(mol.get_scf_energy()))
            if delta < self.driver.IMAGE_TOL:
                q_prev = q_new
                converged = True
                break
            if _image_field_stagnant(it, delta, e_hist, self.driver):
                dump_log(mol, title=(f"PyOQP: QM-image field (reference density) accepted on "
                                     f"stagnation after {it + 1} iterations: max |dq| = {delta:.2e} e, "
                                     f"SCF energy stable to {max(e_hist[-3:]) - min(e_hist[-3:]):.1e} "
                                     f"Hartree over three iterations"), section='')
                # keep q_prev: it is the field this SCF (and the excitation that
                # follows) was computed in, so the active-state loop measures its
                # first residual against the field the electrons actually saw
                converged = True
                break
            q_prev = 0.5 * (q_new + q_prev) if it > 6 else q_new
            # Warm start: keep the converged orbitals as the guess for the next
            # image iteration and only rebuild the bare one-electron integrals
            # (the ESPF term is re-added above).  A fresh Hueckel guess every
            # iteration lets a reference with two nearby SCF solutions flip
            # between them as the image field changes, and the loop then
            # oscillates instead of converging.
            ints_1e(mol)
        if not converged:
            raise RuntimeError(
                f"Periodic ESPF QM/MM NAMD: the QM-image charge self-consistency "
                f"did not converge in {it + 1} iterations "
                f"(max |dq| = {delta:.2e} e > {self.driver.IMAGE_TOL:.0e}, and the "
                f"charges did not stagnate below {self.driver.IMAGE_TOL_STAGNANT:.0e} e "
                f"with a stable SCF energy); the "
                "energy/force would be inconsistent.  Tighten [scf] conv or "
                "check the QM/MM contacts.")
        self._grad_cache = None
        self._img_ctx = None
        ref = [mol.get_scf_energy()]
        if with_overlap:
            mol.back_door = (self.prev_xyz, self.prev_data)
            BasisOverlap(mol).overlap()
        sp.excitation(ref)
        LastStep(mol).compute(mol)
        if psi_img is None:
            self._e_img, self._f_img = 0.0, None
            return potmm, potqm
        # The reference-density loop above only seeds the field; the state
        # that is propagated gets its own self-consistent field below (and
        # again after a surface hop, see _total_force).
        self._img_ctx = dict(sp=sp, potmm_mm=potmm_mm, psi_img=psi_img, dpsi_img=dpsi_img,
                             with_overlap=with_overlap, geom=self._geometry_key())
        return self._refine_image_field(q_prev), potqm

    # -- restart: the periodic image-field seed --------------------------- #
    def _restart_extra_payload(self):
        """Checkpoint the converged active-state ESPF charges of the periodic
        image field, so a resumed trajectory warm-starts the image iteration
        exactly like an uninterrupted one (a zero seed can settle on a
        different SCF/image branch when two nearby solutions exist)."""
        q = getattr(self, "_q_img", None)
        if q is None:
            return {}
        return {"qmmm_q_img": np.asarray(q, dtype=np.float64).reshape(-1)}

    def _load_restart_extra(self, saved, prev_data=None):
        del prev_data
        if "qmmm_q_img" not in saved:
            return {}
        q = np.asarray(saved["qmmm_q_img"], dtype=float).reshape(-1)
        if q.shape != (self.natom,) or not np.all(np.isfinite(q)):
            raise RuntimeError(
                "NAMD restart checkpoint has an invalid periodic QM-image charge seed "
                f"(shape {q.shape}, expected ({self.natom},))")
        return {"q_img": q.copy()}

    def _restore_restart_extra(self, extra):
        if extra:
            self._q_img = np.array(extra["q_img"], dtype=float, copy=True)

    def _absorb_hop_field_shift(self, jump):
        """Periodic full-ESPF, after an accepted hop: the potential energy of
        the new state moved by ``jump`` (Hartree) when its image field was
        refined, after the hop kernel had already rescaled the velocities.
        Rescale the real QM-atom velocities uniformly (the hop kernel's own
        degrees of freedom) so the kinetic energy changes by -jump; if they
        do not carry enough kinetic energy, rescale all atoms (uniform scaling
        keeps the rigid-water constraints satisfied).  Returns the residual
        total-energy jump."""
        mol = self.mol
        jump = float(jump)
        if not np.isfinite(jump) or abs(jump) < 1e-14:
            return jump
        for label, idx in (("QM", np.asarray(self.qm_atoms, dtype=int)),
                           ("all", np.arange(self.natom_all))):
            ke = 0.5 * float(np.sum(self.m_all[idx, None] * self.v_all[idx] ** 2))
            if ke > jump and ke > 0.0:
                self.v_all[idx] *= np.sqrt(1.0 - jump / ke)
                dump_log(mol, title=(f"PyOQP: QM-image field refinement after the hop shifted "
                                     f"E({self.active}) by {jump:+.3e} Hartree; absorbed into "
                                     f"the {label}-atom kinetic energy"), section='')
                return 0.0
        raise RuntimeError(
            f"NAMD_QMMM: the QM-image field refinement after the hop raised the energy of "
            f"state {self.active} by {jump:.3e} Hartree, more than the available kinetic "
            f"energy; the hop is not allowed under the refined Hamiltonian.")

    def _geometry_key(self):
        """Full-system coordinates the current electronic state belongs to."""
        return np.array(self.r_all, dtype=float, copy=True)

    def _refine_image_field(self, q_prev):
        """Make the periodic QM-image field self-consistent with the relaxed
        ESPF charges of the ACTIVE state (the force is the derivative of the
        energy only if the field the SCF saw belongs to the propagated state):
        iterate SCF -> excitation -> active-state gradient (which publishes
        the relaxed charges) until those charges reproduce the field they were
        computed in.  Seeds from ``q_prev`` (the reference-density charges at a
        new geometry, or the previous state's field after a surface hop),
        stores the image energy/force terms and caches the converged gradient
        for _qm_gradient; returns the MM potential incl. the image field."""
        from oqp.library.qmmm_driver import (
            unpack_lower_tri_single, unpack_lower_tri_multi, pack_lower_tri_single)
        mol = self.mol
        ctx = self._img_ctx
        sp, potmm_mm, psi_img, dpsi_img = ctx["sp"], ctx["potmm_mm"], ctx["psi_img"], ctx["dpsi_img"]
        nat = mol.data["natom"]
        nbf = mol.data.get_basis()["nbf"]
        from oqp.library.qmmm_driver import anderson_step
        q_prev = np.array(q_prev, dtype=float)
        potmm = potmm_mm + psi_img @ q_prev
        self._grad_cache = None
        converged = False
        delta = float("inf")
        q_hist, f_hist, e_hist = [], [], []
        for k in range(int(self.driver.IMAGE_MAXITER_ACTIVE)):
            g = self._qm_gradient()
            q_act = np.array(mol.data["OQP::partial_charges"], dtype=float)
            delta = float(np.abs(q_act - q_prev).max())
            e_act = float(mol.energies[self.active])
            dump_log(mol, title=(f"PyOQP: QM-image field, active-state iteration {k + 1}: "
                                 f"max |dq| = {delta:.2e} e, E({self.active}) = "
                                 f"{e_act:.10f} Hartree"), section='')
            if delta < self.driver.IMAGE_TOL_ACTIVE:
                converged = True
                break
            # The relaxed charges carry the Z-vector residual (about 1e-4 e at
            # the default zvconv for an 18-atom indole), so |dq| can settle at
            # a noise floor above the tolerance while the state energy no
            # longer moves: accept when the residual is within ten times the
            # tolerance and the energy has been stationary over three
            # iterations.  A tighter [tdhf] zvconv lowers the floor.
            if (delta < 10.0 * self.driver.IMAGE_TOL_ACTIVE and len(e_hist) >= 2
                    and max(abs(e_act - e) for e in e_hist[-2:]) < self.driver.IMAGE_ETOL_ACTIVE):
                dump_log(mol, title=(f"PyOQP: QM-image field accepted on energy stagnation "
                                     f"(|dE| < {self.driver.IMAGE_ETOL_ACTIVE:.0e} Hartree over three "
                                     f"iterations, max |dq| = {delta:.2e} e)"), section='')
                converged = True
                break
            if len(f_hist) >= 2 and delta > np.abs(f_hist[-1]).max() > np.abs(f_hist[-2]).max():
                # residual grew twice in a row: the history is dominated by
                # noise, restart the extrapolation from the current point
                q_hist, f_hist = [], []
            q_hist.append(q_prev.copy()); f_hist.append(q_act - q_prev); e_hist.append(e_act)
            q_prev = q_act if k == 0 else anderson_step(q_hist, f_hist)
            ints_1e(mol)                                   # bare hcore, orbitals kept
            potmm = potmm_mm + psi_img @ q_prev
            mol.data["OQP::POTMM"] = potmm
            mol.data["OQP::POTQM"] = np.zeros((nat, nat))
            oqp.espf_op_corr(mol)
            espf = unpack_lower_tri_multi(mol.data["OQP::ESPF_CORR"], nbf, nat)
            hcore = unpack_lower_tri_single(mol.get_hcore(), nbf)
            hcore += np.einsum("ijk,i->jk", espf, potmm)
            mol.set_hcore(pack_lower_tri_single(hcore))
            self._embedded_scf(sp, warm=True)
            ref = [mol.get_scf_energy()]
            if ctx["with_overlap"]:
                mol.back_door = (self.prev_xyz, self.prev_data)
                BasisOverlap(mol).overlap()
            sp.excitation(ref)
            LastStep(mol).compute(mol)
        if not converged:
            raise RuntimeError(
                f"Periodic ESPF QM/MM NAMD: the QM-image field did not become "
                f"self-consistent with the active-state charges in "
                f"{self.driver.IMAGE_MAXITER_ACTIVE} iterations (max |dq| = {delta:.2e} e).")
        self._q_img = q_prev.copy()
        self._e_img = 0.5 * float(q_prev @ psi_img @ q_prev)
        self._f_img = -np.einsum("a,b,abc->ac", q_prev, q_prev, dpsi_img)
        self._grad_cache = (int(self.active), np.array(g, dtype=float), q_act.copy(), ctx["geom"])
        return potmm

    def _qm_gradient(self):
        """Embedded active-state gradient (Hartree/bohr) incl. ESPF force."""
        import os
        mol = self.mol
        cache = getattr(self, "_grad_cache", None)
        if cache is not None and cache[0] == int(self.active) \
                and np.array_equal(cache[3], self.r_all):
            # gradient already evaluated by the image self-consistency loop
            # for this state at this geometry; restore its relaxed charges
            self._grad_cache = None
            mol.data["OQP::partial_charges"] = cache[2].copy()
            return cache[1].copy()
        self._grad_cache = None
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
        ctx = getattr(self, "_img_ctx", None)
        if ctx is not None:
            # Periodic full-ESPF: the image field stored by _electronic_qmmm
            # belongs to the state it was refined for.  After a surface hop
            # the new active state must get its own self-consistent field
            # (seeded from the previous state's charges) before its force is
            # integrated; a geometry change needs a new electronic step.
            if not np.array_equal(ctx["geom"], self.r_all):
                raise RuntimeError("NAMD_QMMM._total_force: the geometry changed since the "
                                   "last electronic step; call _electronic_qmmm first.")
            cache = getattr(self, "_grad_cache", None)
            if cache is None or cache[0] != int(self.active):
                dump_log(mol, title=(f"PyOQP: QM-image field re-iterated for the new active "
                                     f"state {self.active} (surface hop)"), section='')
                potmm = self._refine_image_field(self._q_img)
        # active-state embedded QM gradient (Ha/bohr). The z-vector step inside
        # the gradient already forms the excited-state ESPF charges, so
        # OQP::partial_charges holds the active state's QM charges afterwards.
        gqm = self._qm_gradient()
        pchg = np.array(mol.data["OQP::partial_charges"])

        if getattr(self.driver, "espf_full", False):
            force, epot = self._total_force_espf(potmm, gqm, pchg)
            return self._apply_odp_to_force_energy(force, epot)

        # MM forces with embedded QM charges (OpenMM units).  Link atoms are
        # not MM particles: fold each link charge onto its QM host so the
        # charge the QM region presents to the MM electrostatics is conserved
        # (same rule as OpenQpQMMM.compute_force).
        emm_q, gmm_q = self.driver.forces_mm(self._fold_link_charges(pchg))
        gmm = np.array(gmm_q.value_in_unit(u.kilojoule_per_mole / u.nanometer)) / HABOHR_TO_KJMOLNM
        emm = emm_q.value_in_unit(u.kilojoule_per_mole) * KJMOL_TO_HARTREE

        # total force = MM forces; on QM atoms subtract the QM gradient
        # (link-atom rows chain-ruled onto their QM and MM hosts)
        f_all = gmm.copy()
        gqm, g_mm_host = self._project_link_rows(gqm)
        for k, i in enumerate(self.qm_atoms):
            f_all[i] = f_all[i] - gqm[k]
        for m, gl in g_mm_host.items():
            f_all[m] = f_all[m] - gl
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
        if getattr(self, "_f_img", None) is not None:
            fq = fq + self._f_img                               # periodic QM-image force

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
        eqm -= float(getattr(self, "_e_img", 0.0))              # image double counting
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
            restraint_force, restraint_energy = self._evaluate_conservative_restraints(
                self.r_all, self.m_all)
            f_all = f_all + restraint_force
            epot = epot + restraint_energy
            accel = f_all / self.m_all[:, None]
            self._rattle(self.r_all, self.v_all)      # constrained velocities
            self._thermalize_initial()
            self.prev_xyz = copy.deepcopy(self._qm_positions_bohr().reshape(-1))
            self.prev_data = copy.deepcopy(mol.get_data())
            self._log_qmmm(0, epot)
            self._save_restart(0, self.r_all, self.v_all, accel)
            start_step = 0
        else:
            self.r_all = restart['coordinates'].reshape((self.natom_all, 3))
            self.v_all = restart['velocities'].reshape((self.natom_all, 3))
            accel = restart['acceleration'].reshape((self.natom_all, 3))
            self.vel = self._qm_velocities()
            self._sync_positions()
            start_step = restart['step']

        for istep in range(start_step + 1, self.nstep + 1):
            # velocity-Verlet position update (all atoms) + SHAKE (rigid MM water)
            # (fixed dt: the same-spin path uses the Fortran hop kernel with dt_fs)
            r_old = self.r_all.copy()
            # Each term carries the mask separately so the sum keeps the
            # association it had before: multiplying by an exact 1.0 is exact,
            # so a run with nothing held is bit-identical to the old expression.
            self.r_all = (self.r_all + self.v_all * self.dt * self._move_mask
                          + 0.5 * accel * self.dt ** 2 * self._move_mask)
            self._shake(r_old, self.r_all, self.v_all, self.dt)
            self._sync_positions()

            # embedded electronic structure at the new geometry
            potmm, _ = self._electronic_qmmm(with_overlap=True)
            f_all, epot = self._total_force(potmm)
            restraint_force, restraint_energy = self._evaluate_conservative_restraints(
                self.r_all, self.m_all)
            f_all = f_all + restraint_force
            epot = epot + restraint_energy
            accel_new = f_all / self.m_all[:, None]

            # velocity-Verlet velocity update (all atoms) + RATTLE (rigid MM water)
            self.v_all = self.v_all + 0.5 * (accel + accel_new) * self.dt * self._move_mask
            self._rattle(self.r_all, self.v_all)

            # couplings + QM-only FSSH hop
            self._state_overlap(istep)
            self.vel = self._qm_velocities()       # hop sees QM velocities
            active_old = self.active
            energy_before_transition = (
                0.5*np.sum(self.m_all[:, None]*self.v_all**2) + epot)
            hop_ready = self._prepare_hop_step(istep)
            if getattr(self, '_pending_nacme_gate_error', None) is not None:
                new_active, hopped = self.active, False
            else:
                new_active, hopped = self._hop(allow_hop=hop_ready)
            self._store_qm_velocities(self.vel)              # write back rescaled QM velocities
            active_changed = new_active != active_old
            if active_changed:
                self.active = new_active
                f_all, epot = self._total_force(potmm)
                f_all, epot = self._add_last_conservative_restraints(
                    f_all, epot)
                accel_new = f_all / self.m_all[:, None]
                energy_after_transition = (
                    0.5*np.sum(self.m_all[:, None]*self.v_all**2) + epot)
                transition_energy_jump = (
                    energy_after_transition - energy_before_transition)
                if getattr(self, "_img_ctx", None) is not None:
                    # The hop kernel rescaled the velocities with the target
                    # energy evaluated in the previous state's image field;
                    # _total_force has since refined the field for the new
                    # state, shifting its energy.  Absorb that shift into the
                    # QM kinetic energy so the total energy is continuous
                    # across the hop.
                    transition_energy_jump = self._absorb_hop_field_shift(
                        transition_energy_jump)
            else:
                transition_energy_jump = np.nan

            self._apply_thermostat(istep)
            accel = accel_new
            self.prev_xyz = copy.deepcopy(self._qm_positions_bohr().reshape(-1))
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
            title=(f'QMMM-NAMD step {istep:6d}  t={(self._physical_time_fs(istep)):9.3f} fs  '
                   f'active={self.active}  E_tot={ekin+epot:.8f}  '
                   f'E_pot={epot:.8f}  '
                   f'U_ODP={(0.0 if self._odp_last is None else self._odp_last["energy"]):.8f}  '
                   f'E_drop={getattr(self, "_droplet_energy", 0.0):.8f}  '
                   f'drop_max={getattr(self, "_droplet_max_penetration", 0.0):.5f}bohr  '
                   f'drop_n={getattr(self, "_droplet_active_count", 0)}  '
                   f'drop_fmax={getattr(self, "_droplet_force_max", 0.0):.3e}Ha/bohr  '
                   f'E_com={getattr(self, "_solute_com_energy", 0.0):.8f}  '
                   f'E_kin={ekin:.8f}  '
                   f'dQ_therm={getattr(self, "_thermostat_exchange", 0.0):+.3e}  '
                   f'Q_therm={getattr(self, "_thermostat_exchange_cumulative", 0.0):+.3e}  '
                   f'hop={hopped}  {self._hop_rng_log()}  '
                   f'pop={np.array2string(pops, precision=4)}'),
        )
        self._write_md_trajectory(
            istep, self.r_all, epot, ekin, hopped)
        self._enforce_nacme_gate()
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
        # Generic trajectory/restart machinery follows the propagated SOC
        # basis, rather than one spatial MRSF spin manifold.
        self.nstate = self.nstate_soc
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

    def _tracking_state_count(self):
        """Molecule tracking tags contain spatial roots, not SOC sublevels."""
        return self.nstate_mrsf

    def _trajectory_state_energies(self):
        fallback = getattr(self, 'prev_eval', None)
        energies = np.asarray(
            getattr(self, '_trajectory_energies', fallback), dtype=float
        ).reshape(-1)
        if energies.shape != (self.nstate_soc,):
            raise ValueError('SOC trajectory energy vector has the wrong shape')
        return energies

    def _restart_extra_payload(self):
        """Serialize the SOC gauge and both MRSF response histories."""
        required = {
            'soc_prev_u': np.asarray(self.prev_u, dtype=np.complex128),
            'soc_prev_eval': np.asarray(self.prev_eval, dtype=np.float64),
            'soc_prev_sbvec': np.asarray(self.prev_sbvec),
            'soc_prev_tbvec': np.asarray(self.prev_tbvec),
        }
        n = self.nstate_soc
        if required['soc_prev_u'].shape != (n, n):
            raise RuntimeError('invalid SOC eigenvector history at checkpoint')
        if required['soc_prev_eval'].shape != (n,):
            raise RuntimeError('invalid SOC energy history at checkpoint')
        if (required['soc_prev_sbvec'].size == 0
                or required['soc_prev_tbvec'].size == 0):
            raise RuntimeError('empty SOC response-vector history at checkpoint')
        current_vectors = {
            'soc_prev_sbvec': np.asarray(self.sbvec),
            'soc_prev_tbvec': np.asarray(self.tbvec),
        }
        for name, current in current_vectors.items():
            previous = required[name]
            if previous.shape != current.shape or previous.dtype != current.dtype:
                raise RuntimeError(
                    f'{name} does not match the current TD response-vector '
                    'shape and dtype at checkpoint')
        for name, value in required.items():
            if value.dtype.kind in 'fc' and not np.all(np.isfinite(value)):
                raise RuntimeError(f'non-finite {name} at checkpoint')
        return {
            'soc_prev_u_real': required['soc_prev_u'].real,
            'soc_prev_u_imag': required['soc_prev_u'].imag,
            'soc_prev_eval': required['soc_prev_eval'],
            'soc_prev_sbvec': required['soc_prev_sbvec'],
            'soc_prev_tbvec': required['soc_prev_tbvec'],
        }

    def _load_restart_extra(self, saved, prev_data=None):
        """Validate SOC gauge/history arrays before exposing a checkpoint."""
        names = (
            'soc_prev_u_real', 'soc_prev_u_imag', 'soc_prev_eval',
            'soc_prev_sbvec', 'soc_prev_tbvec')
        missing = [name for name in names if name not in saved]
        if missing:
            raise RuntimeError(
                'SOC restart checkpoint lacks ' + ', '.join(missing))
        n = self.nstate_soc
        u_real = np.asarray(saved['soc_prev_u_real'], dtype=float)
        u_imag = np.asarray(saved['soc_prev_u_imag'], dtype=float)
        prev_eval = np.asarray(saved['soc_prev_eval'], dtype=float)
        sbvec = np.array(saved['soc_prev_sbvec'], copy=True)
        tbvec = np.array(saved['soc_prev_tbvec'], copy=True)
        if u_real.shape != (n, n) or u_imag.shape != (n, n):
            raise RuntimeError('SOC restart checkpoint has invalid eigenvectors')
        if prev_eval.shape != (n,):
            raise RuntimeError('SOC restart checkpoint has invalid energies')
        response_tags = {
            'singlet': (sbvec, 'OQP::td_bvec_mo_s'),
            'triplet': (tbvec, 'OQP::td_bvec_mo_t'),
        }
        if not isinstance(prev_data, dict):
            raise RuntimeError(
                'SOC restart checkpoint lacks TD response-vector metadata')
        for label, (saved_vector, tag) in response_tags.items():
            if tag not in prev_data:
                raise RuntimeError(
                    f'SOC restart checkpoint lacks {label} TD dimensions')
            expected = np.asarray(prev_data[tag])
            if (saved_vector.size == 0
                    or saved_vector.shape != expected.shape
                    or saved_vector.dtype != expected.dtype
                    or saved_vector.dtype.kind not in 'fc'):
                raise RuntimeError(
                    f'SOC restart checkpoint has incompatible {label} '
                    'response-vector shape or dtype')
        for value in (u_real, u_imag, prev_eval, sbvec, tbvec):
            if value.dtype.kind in 'fc' and not np.all(np.isfinite(value)):
                raise RuntimeError('SOC restart checkpoint has non-finite history')
        prev_u = u_real + 1j*u_imag
        gram = prev_u.conj().T @ prev_u
        if not np.allclose(gram, np.eye(n), atol=1.0e-7, rtol=0.0):
            raise RuntimeError('SOC restart checkpoint eigenvectors are not unitary')
        return {
            'prev_u': prev_u, 'prev_eval': prev_eval.copy(),
            'prev_sbvec': sbvec, 'prev_tbvec': tbvec,
        }

    def _restore_restart_extra(self, extra):
        for name in ('prev_u', 'prev_eval', 'prev_sbvec', 'prev_tbvec'):
            setattr(self, name, np.array(extra[name], copy=True))

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

    def _soc_unitary_overlap(self, t):
        """Return the nearest-unitary SOC overlap and its anti-Hermitian log.

        The complex generator is also retained in the packed trajectory.  It
        must not be forced through the real antisymmetric same-spin NACME gate:
        SOC adiabatic derivative couplings are generally complex and
        anti-Hermitian.
        """
        from scipy.linalg import sqrtm, logm
        overlap = np.asarray(t, dtype=np.complex128)
        metric_root = np.asarray(
            sqrtm(overlap.conj().T @ overlap), dtype=np.complex128)
        tu = overlap @ np.linalg.inv(metric_root)
        kgen = np.asarray(logm(tu), dtype=np.complex128)
        kgen = 0.5 * (kgen - kgen.conj().T)
        self._last_state_overlap = overlap.copy()
        self._last_overlap_tdc = kgen / self.dt
        return tu, kgen

    def _propagate_and_hop(self, eval_prev, eval_cur, t, allow_hop=True):
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
        from scipy.linalg import expm
        n = self.nstate_soc
        a = self.active - 1
        dt = self.dt
        nsub = max(1, self.substep)

        tu, kgen = self._soc_unitary_overlap(t)
        # substep local diabatization: split the basis rotation into nsub equal
        # fractional rotations (tu^{1/nsub}) and integrate the energy phase with
        # linearly interpolated diagonal energies.  The net full-step propagator
        # p is accumulated and used for the SHARC flux hop probabilities.
        # Reduces exactly to the single-step LD propagator when nsub = 1.
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
        self._last_hop_probabilities = np.zeros((n, n), dtype=float)
        self._last_hop_probabilities[a, :] = cmhp

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
        if not allow_hop:
            return False
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
        restart = self._load_restart()
        if restart is None:
            r = mol.get_system().reshape((self.natom, 3))
            eval_ha, u = self._electronic_soc(with_overlap=False)
            self._resolve_initial_active(u)
            grad, e_pure, mult, state, w = self._soc_gradient(
                u, self.active, eval_ha)
            accel = -grad / self.mass[:, None]
            self._ulog = u
            self._trajectory_energies = (
                self.e_ref + self.e0 + np.asarray(eval_ha, dtype=float))
            self._store_prev(r, u, eval_ha)
            self._log_soc(0, e_pure, mult, state, w, False)
            self._save_restart(0, r, self.vel, accel)
            start_step = 0
        else:
            r = restart['coordinates'].reshape((self.natom, 3))
            self.vel = restart['velocities'].reshape((self.natom, 3))
            accel = restart['acceleration'].reshape((self.natom, 3))
            mol.update_system(r.reshape(-1))
            start_step = restart['step']
        self._e_ref_tot = (
            self._nve_reference_energy
            if self._nve_reference_energy is not None
            else 0.5*np.sum(self.mass[:, None]*self.vel**2)
                 + (e_pure if restart is None else 0.0))

        for istep in range(start_step + 1, self.nstep + 1):
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
            self._last_state_overlap = np.array(t, copy=True)

            # active-surface force (weighted-MCH diagonal gradient) + vel update
            grad, e_pure, mult, state, w = self._soc_gradient(u, self.active, eval_ha)
            accel_new = -grad / self.mass[:, None]
            self.vel = self.vel + 0.5 * (accel + accel_new) * self.dt

            # local-diabatization propagation + fewest-switches hop
            active_old = self.active
            energy_before_transition = (
                0.5*np.sum(self.mass[:, None]*self.vel**2) + e_pure)
            allow_hop = self._prepare_hop_step(istep)
            hopped = self._propagate_and_hop(
                self.prev_eval, eval_ha, t, allow_hop=allow_hop)
            if hopped:
                grad, e_pure, mult, state, w = self._soc_gradient(u, self.active, eval_ha)
                accel_new = -grad / self.mass[:, None]
                transition_energy_jump = (
                    0.5*np.sum(self.mass[:, None]*self.vel**2) + e_pure
                    - energy_before_transition)
            else:
                transition_energy_jump = np.nan

            accel = accel_new
            if self.econs:                                 # temporary E_tot-conservation rescale
                ke = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
                ket = self._e_ref_tot - e_pure
                if ke > 0 and ket > 0:
                    self.vel *= np.sqrt(ket / ke)
            self._ulog = u
            self._trajectory_energies = (
                self.e_ref + self.e0 + np.asarray(eval_ha, dtype=float))
            self._store_prev(r, u, eval_ha)
            self._log_soc(
                istep, e_pure, mult, state, w, hopped,
                transition_energy_jump=transition_energy_jump)
            self._save_restart(istep, r, self.vel, accel)

        dump_log(mol, title='PyOQP: SOC-NAMD trajectory complete')

    def _store_prev(self, r, u, eval_ha):
        self.prev_xyz = copy.deepcopy(r.reshape(-1))
        self.prev_data = copy.deepcopy(self.mol.get_data())
        self.prev_u = u.copy()
        self.prev_eval = eval_ha.copy()
        self.prev_sbvec = self.sbvec.copy()
        self.prev_tbvec = self.tbvec.copy()

    def _log_soc(self, istep, e_pure, mult, state, w, hopped,
                 transition_energy_jump=np.nan):
        ekin = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
        self._update_nve_gate(
            istep, e_pure, ekin, transition_energy_jump)
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
            title=(f'SOC-NAMD step {istep:6d}  t={(self._physical_time_fs(istep)):9.3f} fs  '
                   f'active={self.active}  E_tot={ekin+e_pure:.8f}  '
                   f'E_pure={e_pure:.8f}  E_kin={ekin:.8f}  hop={hopped}  '
                   f'{self._hop_rng_log()}  '
                   f'dom=({self._mch_label(mult, state)},w={w:.3f})  '
                   f'pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )
        self._unbiased_potential_energy = float(e_pure)
        self._trajectory_state_energies = np.asarray(
            self._trajectory_energies, dtype=float).reshape(-1)
        self._write_md_trajectory(istep, self.mol.get_system().reshape(
            (self.natom, 3)), e_pure, ekin, hopped)
        self._enforce_nve_gate()


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

    def _mch_propagate_and_hop(self, h_mch, e_mch, allow_hop=True):
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
        self._last_hop_probabilities = np.zeros((n, n), dtype=float)
        self._last_hop_probabilities[a, :] = cmhp

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

        if not allow_hop:
            return False
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
        restart = self._load_restart()
        if restart is None:
            r = mol.get_system().reshape((self.natom, 3))
            eval_ha, u = self._electronic_soc(with_overlap=False)
            self._resolve_initial_mch_active()
            h_mch = self._mch_hamiltonian_from_u(u, eval_ha)
            e_mch = self._mch_energies_abs()
            grad, e_pure, mult, state = self._mch_exact_gradient(self.active)
            accel = -grad / self.mass[:, None]
            self._trajectory_energies = e_mch.copy()
            self._store_prev(r, u, eval_ha)
            self._log_mch(0, e_pure, mult, state, False)
            self._save_restart(0, r, self.vel, accel)
            start_step = 0
        else:
            r = restart['coordinates'].reshape((self.natom, 3))
            self.vel = restart['velocities'].reshape((self.natom, 3))
            accel = restart['acceleration'].reshape((self.natom, 3))
            mol.update_system(r.reshape(-1))
            start_step = restart['step']
        self._e_ref_tot = (
            self._nve_reference_energy
            if self._nve_reference_energy is not None
            else 0.5*np.sum(self.mass[:, None]*self.vel**2)
                 + (e_pure if restart is None else 0.0))

        for istep in range(start_step + 1, self.nstep + 1):
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

            energy_before_transition = (
                0.5*np.sum(self.mass[:, None]*self.vel**2) + e_pure)
            allow_hop = self._prepare_hop_step(istep)
            hopped = self._mch_propagate_and_hop(
                h_mch, e_mch, allow_hop=allow_hop)
            if hopped:
                grad, e_pure, mult, state = self._mch_exact_gradient(self.active)
                accel_new = -grad / self.mass[:, None]
                transition_energy_jump = (
                    0.5*np.sum(self.mass[:, None]*self.vel**2) + e_pure
                    - energy_before_transition)
            else:
                transition_energy_jump = np.nan

            accel = accel_new
            if self.econs:
                ke = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
                ket = self._e_ref_tot - e_pure
                if ke > 0 and ket > 0:
                    self.vel *= np.sqrt(ket / ke)
            self._trajectory_energies = e_mch.copy()
            self._store_prev(r, u, eval_ha)
            self._log_mch(
                istep, e_pure, mult, state, hopped,
                transition_energy_jump=transition_energy_jump)
            self._save_restart(istep, r, self.vel, accel)

        dump_log(mol, title='PyOQP: SOC-MCH-NAMD trajectory complete')

    def _log_mch(self, istep, e_pure, mult, state, hopped,
                 transition_energy_jump=np.nan):
        ekin = 0.5 * np.sum(self.mass[:, None] * self.vel ** 2)
        self._update_nve_gate(
            istep, e_pure, ekin, transition_energy_jump)
        pmch = np.abs(self.coef) ** 2
        pop_s = float(pmch[:self.ns].sum())
        pop_t = float(pmch[self.ns:].sum())
        dump_log(
            self.mol,
            title=(f'SOC-MCH-NAMD step {istep:6d}  t={(self._physical_time_fs(istep)):9.3f} fs  '
                   f'active={self.active}:{self._mch_active_label(self.active)}  '
                   f'E_tot={ekin+e_pure:.8f}  E_pure={e_pure:.8f}  '
                   f'E_kin={ekin:.8f}  hop={hopped}  '
                   f'{self._hop_rng_log()}  '
                   f'grad={self._mch_label(mult, state)}  pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )
        self._unbiased_potential_energy = float(e_pure)
        self._trajectory_state_energies = np.asarray(
            self._trajectory_energies, dtype=float).reshape(-1)
        self._write_md_trajectory(istep, self.mol.get_system().reshape(
            (self.natom, 3)), e_pure, ekin, hopped)
        self._enforce_nve_gate()


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
    _soc_unitary_overlap = NAMD_SOC._soc_unitary_overlap
    _tracking_state_count = NAMD_SOC._tracking_state_count
    _trajectory_state_energies = NAMD_SOC._trajectory_state_energies
    _restart_extra_payload = NAMD_SOC._restart_extra_payload
    _load_restart_extra = NAMD_SOC._load_restart_extra
    _restore_restart_extra = NAMD_SOC._restore_restart_extra

    def __init__(self, mol):
        super().__init__(mol)                                  # NAMD_QMMM: OpenMM + QM masses + v_all
        # spin-adiabatic electronic state space (ns singlets + 3 nt triplet Ms)
        self.nstate_mrsf = int(mol.config['tdhf']['nstate'])
        self.ns = self.nstate_mrsf
        self.nt = self.nstate_mrsf
        self.nstate_soc = self.ns + 3 * self.nt
        self.nstate = self.nstate_soc
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
        ewald = self.driver._ewald() if getattr(self.driver, "espf_full", False) else None
        if ewald is not None:
            # The periodic QM-image term must be self-consistent with the
            # charges of the propagated state; the spin-adiabatic SOC state is
            # a weighted mixture of MCH states whose relaxed charges are not
            # available per iteration, so periodic SOC-NAMD is not offered.
            # (The input checker reports the same restriction at parse time.)
            raise NotImplementedError(
                "SOC-NAMD QM/MM is implemented for non-periodic clusters only "
                "([qmmm] cutoff=NoCutoff or CutoffNonPeriodic); the periodic "
                "QM-image term is not available for the spin-mixed active state.")
        potmm, potqm = self._embedding_field()

        sp = SinglePoint(mol)
        warm_start = self._start_scf_orbitals(sp)
        nat = mol.data["natom"]
        nbf = mol.data.get_basis()["nbf"]
        potmm_mm = np.asarray(potmm, dtype=float).copy()
        if ewald is not None:
            psi_img, dpsi_img = ewald.qm_image_matrix(self.driver._qm_center_positions_bohr())
            q_prev = (self._q_img if getattr(self, "_q_img", None) is not None
                      and len(self._q_img) == nat else np.zeros(nat))
        else:
            psi_img = dpsi_img = None
            q_prev = None
        converged = psi_img is None
        delta, it = float("inf"), -1
        e_hist = []
        for it in range(int(self.driver.IMAGE_MAXITER)):
            potmm = potmm_mm if psi_img is None else potmm_mm + psi_img @ q_prev
            mol.data["OQP::POTMM"] = potmm
            mol.data["OQP::POTQM"] = np.zeros((nat, nat))
            oqp.espf_op_corr(mol)
            espf = unpack_lower_tri_multi(mol.data["OQP::ESPF_CORR"], nbf, nat)
            hcore = unpack_lower_tri_single(mol.get_hcore(), nbf)
            hcore += np.einsum("ijk,i->jk", espf, potmm)
            mol.set_hcore(pack_lower_tri_single(hcore))
            self._embedded_scf(sp, warm=(warm_start or it > 0))
            self._scf_orbitals_ready = True
            if psi_img is None:
                break
            oqp.form_esp_charges(mol)
            q_new = np.array(mol.data["OQP::partial_charges"], dtype=float)
            delta = float(np.abs(q_new - q_prev).max())
            e_hist.append(float(mol.get_scf_energy()))
            if delta < self.driver.IMAGE_TOL:
                q_prev = q_new
                converged = True
                break
            if _image_field_stagnant(it, delta, e_hist, self.driver):
                dump_log(mol, title=(f"PyOQP: QM-image field (reference density) accepted on "
                                     f"stagnation after {it + 1} iterations: max |dq| = {delta:.2e} e, "
                                     f"SCF energy stable to {max(e_hist[-3:]) - min(e_hist[-3:]):.1e} "
                                     f"Hartree over three iterations"), section='')
                # keep q_prev: it is the field this SCF (and the excitation that
                # follows) was computed in, so the active-state loop measures its
                # first residual against the field the electrons actually saw
                converged = True
                break
            q_prev = 0.5 * (q_new + q_prev) if it > 6 else q_new
            # Warm start: keep the converged orbitals as the guess for the next
            # image iteration and only rebuild the bare one-electron integrals
            # (the ESPF term is re-added above).  A fresh Hueckel guess every
            # iteration lets a reference with two nearby SCF solutions flip
            # between them as the image field changes, and the loop then
            # oscillates instead of converging.
            ints_1e(mol)
        if not converged:
            raise RuntimeError(
                f"Periodic ESPF QM/MM NAMD: the QM-image charge self-consistency "
                f"did not converge in {it + 1} iterations "
                f"(max |dq| = {delta:.2e} e > {self.driver.IMAGE_TOL:.0e}, and the "
                f"charges did not stagnate below {self.driver.IMAGE_TOL_STAGNANT:.0e} e "
                f"with a stable SCF energy); the "
                "energy/force would be inconsistent.  Tighten [scf] conv or "
                "check the QM/MM contacts.")
        if psi_img is not None:
            self._q_img = q_prev.copy()
            self._e_img = 0.5 * float(q_prev @ psi_img @ q_prev)
            self._f_img = -np.einsum("a,b,abc->ac", q_prev, q_prev, dpsi_img)
        else:
            self._e_img, self._f_img = 0.0, None
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
            self, u, active, eval_ha, self._qm_velocities(kinematic=True))
        g += NAMD_SOC._tdc_gradient_correction(
            self, u, active, getattr(self, '_last_s_mch', None), self._qm_velocities(kinematic=True))

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
            if getattr(self, "_f_img", None) is not None:
                fq = fq + self._f_img                           # periodic QM-image force
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
            eqm -= float(getattr(self, "_e_img", 0.0))          # image double counting
            return f_all, eqm + emm

        # Split scheme: link-atom charges folded onto their QM hosts for the MM
        # electrostatics, link-row gradients chain-ruled onto both hosts (same
        # bookkeeping as NAMD_QMMM._total_force).
        emm_q, gmm_q = self.driver.forces_mm(self._fold_link_charges(pchg))
        gmm = np.array(gmm_q.value_in_unit(u.kilojoule_per_mole / u.nanometer)) / HABOHR_TO_KJMOLNM
        emm = emm_q.value_in_unit(u.kilojoule_per_mole) * KJMOL_TO_HARTREE

        f_all = gmm.copy()
        g_real, g_mm_host = self._project_link_rows(g_qm)
        for k, i in enumerate(self.qm_atoms):
            f_all[i] = f_all[i] - g_real[k]
        for m, gl in g_mm_host.items():
            f_all[m] = f_all[m] - gl
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
        restart = self._load_restart()
        if restart is None:
            self._sync_positions()
            eval_ha, u, potmm, _ = self._electronic_soc_qmmm(
                with_overlap=False)
            NAMD_SOC._resolve_initial_active(self, u)
            g_qm, e_diag, mult, state, w, pchg = self._soc_gradient_qmmm(
                u, self.active, eval_ha)
            f_all, epot = self._total_force_soc(
                potmm, g_qm, e_diag, pchg)
            accel = f_all / self.m_all[:, None]
            self._rattle(self.r_all, self.v_all)
            self._thermalize_initial()
            self._ulog = u
            self._trajectory_energies = (
                self.e_ref + self.e0 + np.asarray(eval_ha, dtype=float))
            r_qm = self._qm_positions_bohr()
            NAMD_SOC._store_prev(self, r_qm, u, eval_ha)
            self._log_soc_qmmm(0, epot, mult, state, w, False)
            self._save_restart(0, self.r_all, self.v_all, accel)
            start_step = 0
        else:
            self.r_all = restart['coordinates'].reshape((self.natom_all, 3))
            self.v_all = restart['velocities'].reshape((self.natom_all, 3))
            accel = restart['acceleration'].reshape((self.natom_all, 3))
            self.vel = self._qm_velocities()
            self._sync_positions()
            start_step = restart['step']
        self._e_ref_tot = (
            self._nve_reference_energy
            if self._nve_reference_energy is not None
            else 0.5*np.sum(self.m_all[:, None]*self.v_all**2)
                 + (epot if restart is None else 0.0))

        for istep in range(start_step + 1, self.nstep + 1):
            # adaptive timestep + velocity-Verlet position update + SHAKE
            # A held atom cannot move, so neither its velocity nor the force on
            # it may shrink the step the moving atoms get.  With nothing held
            # the mask is 1 everywhere and this is the previous call.
            self.dt = self._adaptive_dt(self.v_all * self._move_mask,
                                        accel * self._move_mask)
            self._t_fs += self.dt / FS_TO_AU
            r_old = self.r_all.copy()
            # Each term carries the mask separately so the sum keeps the
            # association it had before: multiplying by an exact 1.0 is exact,
            # so a run with nothing held is bit-identical to the old expression.
            self.r_all = (self.r_all + self.v_all * self.dt * self._move_mask
                          + 0.5 * accel * self.dt ** 2 * self._move_mask)
            self._shake(r_old, self.r_all, self.v_all, self.dt)
            self._sync_positions()

            # embedded spin-adiabatic electronic structure (+ MO overlap)
            eval_ha, u, potmm, _ = self._electronic_soc_qmmm(with_overlap=True)
            s_mch = NAMD_SOC._mch_overlap(self)
            self._last_s_mch = s_mch
            u, t = NAMD_SOC._phase_track(u, self.prev_u, s_mch, eval_ha)
            self._last_state_overlap = np.array(t, copy=True)

            # active-surface force (weighted-MCH diagonal gradient + ESPF) + vel update
            g_qm, e_diag, mult, state, w, pchg = self._soc_gradient_qmmm(u, self.active, eval_ha)
            f_all, epot = self._total_force_soc(potmm, g_qm, e_diag, pchg)
            accel_new = f_all / self.m_all[:, None]
            self.v_all = self.v_all + 0.5 * (accel + accel_new) * self.dt * self._move_mask
            self._rattle(self.r_all, self.v_all)

            # local-diabatization propagation + spin-adiabatic hop (QM velocities only)
            active_old = self.active                           # save for ESPF correction below
            epot_old = epot                                    # total E_pot before hop
            energy_before_transition = (
                0.5*np.sum(self.m_all[:, None]*self.v_all**2) + epot)
            self.vel = self._qm_velocities()
            allow_hop = self._prepare_hop_step(istep)
            hopped = NAMD_SOC._propagate_and_hop(
                self, self.prev_eval, eval_ha, t, allow_hop=allow_hop)
            self._store_qm_velocities(self.vel)
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
                transition_energy_jump = (
                    0.5*np.sum(self.m_all[:, None]*self.v_all**2) + epot
                    - energy_before_transition)
            else:
                transition_energy_jump = np.nan

            accel = accel_new
            if self.econs:                                 # temporary E_tot-conservation rescale
                ke = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
                ket = self._e_ref_tot - epot
                if ke > 0 and ket > 0:
                    self.v_all *= np.sqrt(ket / ke)
            self._ulog = u
            self._trajectory_energies = (
                self.e_ref + self.e0 + np.asarray(eval_ha, dtype=float))
            NAMD_SOC._store_prev(self, self._qm_positions_bohr(), u, eval_ha)
            self._log_soc_qmmm(
                istep, epot, mult, state, w, hopped,
                transition_energy_jump=transition_energy_jump)
            self._save_restart(istep, self.r_all, self.v_all, accel)

        dump_log(mol, title='PyOQP: SOC-NAMD QM/MM trajectory complete')

    def _log_soc_qmmm(self, istep, epot, mult, state, w, hopped,
                      transition_energy_jump=np.nan):
        ekin = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
        self._update_nve_gate(
            istep, epot, ekin, transition_energy_jump)
        mch = self._ulog @ self.coef
        pmch = np.abs(mch) ** 2
        pop_s = float(pmch[:self.ns].sum())
        pop_t = float(pmch[self.ns:].sum())
        dump_log(
            self.mol,
            title=(f'SOC-QMMM-NAMD step {istep:6d}  t={(self._physical_time_fs(istep)):9.3f} fs  '
                   f'active={self.active}  E_tot={ekin+epot:.8f}  '
                   f'E_pot={epot:.8f}  E_kin={ekin:.8f}  hop={hopped}  '
                   f'{self._hop_rng_log()}  '
                   f'dom=({NAMD_SOC._mch_label(mult, state)},w={w:.3f})  '
                   f'pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )
        self._unbiased_potential_energy = float(epot)
        self._trajectory_state_energies = np.asarray(
            self._trajectory_energies, dtype=float).reshape(-1)
        self._write_md_trajectory(
            istep, self.r_all, epot, ekin, hopped)
        self._enforce_nve_gate()


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
        restart = self._load_restart()
        if restart is None:
            self._sync_positions()
            eval_ha, u, potmm, _ = self._electronic_soc_qmmm(
                with_overlap=False)
            self._resolve_initial_mch_active()
            h_mch = self._mch_hamiltonian_from_u(u, eval_ha)
            e_mch = self._mch_energies_abs()
            g_qm, e_pure, mult, state, pchg = (
                self._mch_exact_gradient_qmmm(self.active))
            f_all, epot = self._total_force_soc(
                potmm, g_qm, e_pure, pchg)
            accel = f_all / self.m_all[:, None]
            self._rattle(self.r_all, self.v_all)
            self._thermalize_initial()
            self._trajectory_energies = e_mch.copy()
            r_qm = self._qm_positions_bohr()
            NAMD_SOC._store_prev(self, r_qm, u, eval_ha)
            self._log_mch_qmmm(0, epot, mult, state, False)
            self._save_restart(0, self.r_all, self.v_all, accel)
            start_step = 0
        else:
            self.r_all = restart['coordinates'].reshape((self.natom_all, 3))
            self.v_all = restart['velocities'].reshape((self.natom_all, 3))
            accel = restart['acceleration'].reshape((self.natom_all, 3))
            self.vel = self._qm_velocities()
            self._sync_positions()
            start_step = restart['step']
        self._e_ref_tot = (
            self._nve_reference_energy
            if self._nve_reference_energy is not None
            else 0.5*np.sum(self.m_all[:, None]*self.v_all**2)
                 + (epot if restart is None else 0.0))

        for istep in range(start_step + 1, self.nstep + 1):
            # A held atom cannot move, so neither its velocity nor the force on
            # it may shrink the step the moving atoms get.  With nothing held
            # the mask is 1 everywhere and this is the previous call.
            self.dt = self._adaptive_dt(self.v_all * self._move_mask,
                                        accel * self._move_mask)
            self._t_fs += self.dt / FS_TO_AU
            r_old = self.r_all.copy()
            # Each term carries the mask separately so the sum keeps the
            # association it had before: multiplying by an exact 1.0 is exact,
            # so a run with nothing held is bit-identical to the old expression.
            self.r_all = (self.r_all + self.v_all * self.dt * self._move_mask
                          + 0.5 * accel * self.dt ** 2 * self._move_mask)
            self._shake(r_old, self.r_all, self.v_all, self.dt)
            self._sync_positions()

            eval_ha, u, potmm, _ = self._electronic_soc_qmmm(with_overlap=False)
            h_mch = self._mch_hamiltonian_from_u(u, eval_ha)
            e_mch = self._mch_energies_abs()
            g_qm, e_pure, mult, state, pchg = self._mch_exact_gradient_qmmm(self.active)
            f_all, epot = self._total_force_soc(potmm, g_qm, e_pure, pchg)
            accel_new = f_all / self.m_all[:, None]
            self.v_all = self.v_all + 0.5 * (accel + accel_new) * self.dt * self._move_mask
            self._rattle(self.r_all, self.v_all)

            active_old = self.active
            epot_old = epot
            energy_before_transition = (
                0.5*np.sum(self.m_all[:, None]*self.v_all**2) + epot)
            self.vel = self._qm_velocities()
            allow_hop = self._prepare_hop_step(istep)
            hopped = self._mch_propagate_and_hop(
                h_mch, e_mch, allow_hop=allow_hop)
            self._store_qm_velocities(self.vel)
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
                transition_energy_jump = (
                    0.5*np.sum(self.m_all[:, None]*self.v_all**2) + epot
                    - energy_before_transition)
            else:
                transition_energy_jump = np.nan

            accel = accel_new
            if self.econs:
                ke = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
                ket = self._e_ref_tot - epot
                if ke > 0 and ket > 0:
                    self.v_all *= np.sqrt(ket / ke)
            self._trajectory_energies = e_mch.copy()
            NAMD_SOC._store_prev(self, self._qm_positions_bohr(), u, eval_ha)
            self._log_mch_qmmm(
                istep, epot, mult, state, hopped,
                transition_energy_jump=transition_energy_jump)
            self._save_restart(istep, self.r_all, self.v_all, accel)

        dump_log(mol, title='PyOQP: SOC-MCH-QMMM-NAMD trajectory complete')

    def _log_mch_qmmm(self, istep, epot, mult, state, hopped,
                      transition_energy_jump=np.nan):
        ekin = 0.5 * np.sum(self.m_all[:, None] * self.v_all ** 2)
        self._update_nve_gate(
            istep, epot, ekin, transition_energy_jump)
        pmch = np.abs(self.coef) ** 2
        pop_s = float(pmch[:self.ns].sum())
        pop_t = float(pmch[self.ns:].sum())
        dump_log(
            self.mol,
            title=(f'SOC-MCH-QMMM-NAMD step {istep:6d}  t={(self._physical_time_fs(istep)):9.3f} fs  '
                   f'active={self.active}:{self._mch_active_label(self.active)}  '
                   f'E_tot={ekin+epot:.8f}  E_pot={epot:.8f}  '
                   f'E_kin={ekin:.8f}  hop={hopped}  '
                   f'{self._hop_rng_log()}  '
                   f'grad={NAMD_SOC._mch_label(mult, state)}  pop[S]={pop_s:.4f} pop[T]={pop_t:.4f}'),
        )
        self._unbiased_potential_energy = float(epot)
        self._trajectory_state_energies = np.asarray(
            self._trajectory_energies, dtype=float).reshape(-1)
        self._write_md_trajectory(
            istep, self.r_all, epot, ekin, hopped)
        self._enforce_nve_gate()


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
