"""OQP molecule class"""
import os
import copy
import json
import platform
import warnings
import numpy as np
import oqp
from oqp.utils.input_parser import OQPConfigParser
from oqp.molden.moldenwriter import MoldenWriter
from .oqpdata import OQPData, OQP_CONFIG_SCHEMA
from oqp.utils.mpi_utils import MPIManager
from oqp.utils.mpi_utils import mpi_get_attr, mpi_dump
from oqp import ffi
from oqp.utils import regression as regkeys
from oqp.utils.json_utils import json_array
from oqp.utils.state_labels import is_mrsf, public_state_label

# Environment variable that opts JSON dumps into "lean" mode: internal
# ``OQP::`` arrays (density/Fock/MO matrices, etc.) and any other non-regression
# data are dropped before the bundle is written, keeping only the keys declared
# in the regression registry (physics + identity/metadata). Used when
# (re)generating committed example *test references*. The full bundle is still
# written by default so the ``guess=json`` restart workflow keeps working.
LEAN_JSON_ENV = 'OQP_LEAN_JSON'
HESSIAN_CACHE_VERSION = 2


def _env_wants_lean_json():
    """True when ``OQP_LEAN_JSON`` is set to a truthy value."""
    return os.environ.get(LEAN_JSON_ENV, '').strip().lower() not in (
        '', '0', 'false', 'no', 'off')


def _load_json(path):
    """Load a JSON file, returning {} if it is missing or unreadable."""
    try:
        with open(path, 'r') as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return {}


def _abs_nested(value):
    """Element-wise magnitude of a (possibly nested) numeric array.

    Used for sign/phase-ambiguous regression values (e.g. non-adiabatic
    couplings) so a global sign flip between builds does not register as a diff.
    """
    return np.abs(np.array(value, dtype=float)).tolist()


class Molecule:
    """
    OQP molecule representation in python
    """

    def __init__(self, project_name, input_file, log,
                 xyz=None, elem=None, mass=None, charge=0, mult=1, silent=0, idx=1):
        self.mpi_manager = MPIManager()
        self.usempi = True
        self.silent = silent
        self.idx = idx

        self.xyz = xyz
        self.elem = elem
        self.mass = mass
        self.charge = charge
        self.mult = mult

        self.control = None
        self.mol_energy = None

        self.data = None
        self.data_allocate()

        self.config = {}
        self.project_name = project_name
        self.input_file = input_file
        self.log = log
        self.log_path = os.path.dirname(log)
        if self.log_path == '':
            self.log_path = os.path.dirname(os.path.abspath(__file__))
        self.energies = None
        self.grads = None
        self.dcm = []  # Nstate, Nstate
        self.nac = []  # Npairs, 3, Natom,
        self.soc = []  # Npairs, 1,
        self.freqs = np.zeros(0)  # 3Natom-6
        self.hessian = np.zeros(0)  # 3Natom, 3Natom
        self.hessian_metadata = {}
        self.modes = np.zeros(0)  # 3Natom-6, 3Natom
        self.inertia = np.zeros(0)  # 3
        self.infrared_intensities = np.zeros(0)
        self.raman_activities = np.zeros(0)
        self.vibrational_intensity_metadata = {}
        self.infrared_mode_dipole_derivatives = np.zeros((0, 3))
        self.raman_mode_polarizability_derivatives = np.zeros((0, 3, 3))
        self.symmetry_metadata = {}
        self.mrsf_ekt_results_by_kind = {}
        # State-tracking tags loaded from a guess are transport history, not a
        # result of the new calculation.  NACME.align_x sets this only after it
        # has refreshed the mapping and gauge in the current run.
        self._state_tracking_fresh = False

        self.tag = [
            'OQP::DM_A', 'OQP::DM_B',
            'OQP::FOCK_A', 'OQP::FOCK_B',
            'OQP::E_MO_A', 'OQP::E_MO_B',
            'OQP::VEC_MO_A', 'OQP::VEC_MO_B',
            'OQP::Hcore', 'OQP::SM', 'OQP::TM', 'OQP::WAO',
            'OQP::td_abxc', 'OQP::td_bvec_mo', 'OQP::td_mrsf_density', 'OQP::td_energies',
            'OQP::mrsf_ekt_density_mo', 'OQP::mrsf_ekt_lagrangian_mo', 'OQP::mrsf_ekt_fock_mo',
            'OQP::mrsf_ekt_orbitals_mo', 'OQP::mrsf_ekt_eigenvalues', 'OQP::mrsf_ekt_strengths',
            'OQP::hf_hessian',
            'OQP::td_states_overlap',
            'OQP::dc_matrix', 'OQP::nac_matrix',
            'OQP::hamiltonian_qmmm', 'OQP::mm_potential', 'OQP::charge_operator', 'OQP::partial_charges',
            'OQP::namd_coef', 'OQP::namd_velocity', 'OQP::namd_params', 'OQP::namd_results',
            'OQP::namd_tdc', 'OQP::namd_eabs', 'OQP::namd_stas',
            'OQP::td_singlet_energies', 'OQP::td_triplet_energies',
            'OQP::td_bvec_mo_s', 'OQP::td_bvec_mo_t',
            'OQP::mo_tracking_order', 'OQP::mo_tracking_phase',
            'OQP::mo_tracking_overlap', 'OQP::mo_tracking_margin',
            'OQP::state_tracking_order', 'OQP::state_tracking_lineage',
            'OQP::state_tracking_raw_order', 'OQP::state_tracking_output_reordered',
            'OQP::state_tracking_phase_step', 'OQP::state_tracking_phase_initial',
            'OQP::state_tracking_previous_phase_initial',
            'OQP::state_tracking_overlap', 'OQP::state_tracking_margin',
            'OQP::soc_eval',
            'OQP::soc_evec_re', 'OQP::soc_evec_im',
            'OQP::soc_hsoc_re', 'OQP::soc_hsoc_im',
        ]
        self.skip_tag = {"rhf": ['OQP::DM_B', 'OQP::FOCK_B', 'OQP::E_MO_B', 'OQP::VEC_MO_B'],
                         "rohf": [],
                         "uhf": []
                         }
        self.config_tag = {
            'json': ['scf_type', 'basis', 'library']
        }
        self.start_time = None
        self.back_door = None

        for tag in self.tag:
            name = tag.replace('OQP::', '').lower()
            getter = lambda self, t=tag: np.array(self.data[t])
            setter = lambda self, val, t=tag: self.data.__setitem__(t, val)
            setattr(self.__class__, f'get_{name}', getter)
            setattr(self.__class__, f'set_{name}', setter)

    def get_atoms(self):
        """
        Get read-only atoms
        """
        natom = self.data["natom"]
        atoms = np.frombuffer(
            oqp.ffi.buffer(self.elem,
                           natom * oqp.ffi.sizeof("double"))
        ).astype(int)

        return copy.deepcopy(atoms)

    @staticmethod
    def _parse_bool_like(value):
        if isinstance(value, bool):
            return value
        if isinstance(value, (int, np.integer)):
            return bool(value)
        if isinstance(value, str):
            lower = value.strip().lower()
            if lower in ['.true.', 'true', 'on', '1', 'yes', 'y', 't', 'full']:
                return True
            if lower in ['.false.', 'false', 'off', '0', 'no', 'n', 'f', '']:
                return False
        return False

    @staticmethod
    def _parse_enabled_mode(value):
        if isinstance(value, str) and value.strip().lower() == 'auto':
            return 'auto'
        if isinstance(value, bool):
            return bool(value)
        if isinstance(value, (int, np.integer)):
            return bool(value)
        if isinstance(value, str):
            lower = value.strip().lower()
            if lower in ['.true.', 'true', '1', 'on', 'yes']:
                return True
            if lower in ['.false.', 'false', '0', 'off', 'no', 'f']:
                return False
        return False

    def initialize_symmetry_metadata(self):
        symmetry = self.config.get('symmetry', {}) if isinstance(self.config, dict) else {}
        requested_point_group = symmetry.get('point_group', 'auto')
        requested_subgroup = symmetry.get('subgroup', 'auto')

        requested_point_group = requested_point_group if isinstance(requested_point_group, str) and requested_point_group else 'auto'
        requested_subgroup = requested_subgroup if isinstance(requested_subgroup, str) and requested_subgroup else 'auto'

        enabled = self._parse_enabled_mode(symmetry.get('enabled', 'true'))

        if enabled == 'auto':
            status = 'auto'
            point_group = requested_point_group.lower() if requested_point_group != 'auto' else 'c1'
            subgroup = requested_subgroup.lower() if requested_subgroup != 'auto' else 'c1'
        elif enabled:
            status = 'enabled'
            point_group = requested_point_group.lower() if requested_point_group != 'auto' else 'c1'
            subgroup = requested_subgroup.lower() if requested_subgroup != 'auto' else 'c1'
        else:
            status = 'disabled'
            point_group = 'c1'
            subgroup = 'c1'

        self.symmetry_metadata = {
            'status': status,
            'enabled': enabled,
            'requested_point_group': requested_point_group,
            'requested_subgroup': requested_subgroup,
            'point_group': point_group,
            'subgroup': subgroup,
            'detected_point_group': point_group,
            'detected_subgroup': subgroup,
            'label_mo': self._parse_bool_like(symmetry.get('label_mo', True)),
            'label_states': self._parse_bool_like(symmetry.get('label_states', True)),
            'label_modes': self._parse_bool_like(symmetry.get('label_modes', True)),
            'use_integral_symmetry': self._parse_bool_like(symmetry.get('use_integral_symmetry', 'False')),
            'use_response_symmetry': self._parse_bool_like(symmetry.get('use_response_symmetry', 'False')),
            'strict': self._parse_bool_like(symmetry.get('strict', False)),
            'tolerance': float(symmetry.get('tolerance', 1.0e-5)),
            'raw': {
                'enabled': symmetry.get('enabled', 'true'),
                'point_group': requested_point_group,
                'subgroup': requested_subgroup,
                'label_mo': symmetry.get('label_mo', True),
                'label_states': symmetry.get('label_states', True),
                'label_modes': symmetry.get('label_modes', True),
                'use_integral_symmetry': symmetry.get('use_integral_symmetry', 'False'),
                'use_response_symmetry': symmetry.get('use_response_symmetry', 'False'),
                'strict': symmetry.get('strict', False),
                'tolerance': symmetry.get('tolerance', 1.0e-5),
            },
        }

        if status != 'disabled':
            self._detect_symmetry_metadata()

        return self.symmetry_metadata

    def _detect_symmetry_metadata(self):
        """Geometry-based point-group detection (metadata only, non-fatal)."""
        try:
            from oqp.library.symmetry_detect import attach_detection_metadata
            atoms = np.asarray(self.get_atoms(), dtype=float).ravel()
            coords = np.asarray(self.get_system(), dtype=float).reshape(-1, 3)
            if atoms.size == 0 or coords.shape[0] != atoms.size:
                return
            attach_detection_metadata(self.symmetry_metadata, atoms, coords)
        except Exception as exc:
            # Detection must never break the run in the metadata-only phase.
            self.symmetry_metadata['detection_error'] = str(exc)
            return

        if self.symmetry_metadata.get('strict') and \
                not self.symmetry_metadata.get('requested_matches_detected', True):
            raise ValueError(
                "symmetry.strict: requested point group "
                f"'{self.symmetry_metadata['requested_point_group']}' does not match "
                f"detected '{self.symmetry_metadata['detected_point_group']}'"
            )

    def _symmetry_labeling_inputs(self):
        """Shared shells + square overlap for symmetry labeling.

        Returns (shells, smat, nbf, None) on success or
        (None, None, 0, reason) when labeling must be skipped.
        """
        from oqp.library.symmetry import _cartesian_shell_size, _spherical_shell_size

        basis = self.data.get_basis()
        if not basis:
            return None, None, 0, 'skipped_no_basis'
        nbf = int(basis['nbf'])
        pairs = [(int(at), int(l)) for at, l in zip(basis['centers'], basis['angs'])]
        if any(l > 4 for _, l in pairs):
            return None, None, 0, 'skipped_unsupported_shells_beyond_g'
        # Prefer the library's own per-shell flag over any dimension test.
        # The merged integral-path fix exports it (oqp_get_basis_spherical ->
        # basis['spherical']), and it is authoritative in a way a dimension
        # test cannot be: Cartesian and spherical sizes agree for l <= 1, so
        # the AO total alone cannot tell "all shells pure" -- which applies the
        # spherical component order to p shells and yields WRONG irrep labels
        # rather than 'mixed' ones -- from OpenQP's real convention. The
        # dimension branches below stay as a fallback for a runtime that
        # predates the export.
        spherical = basis.get('spherical')
        if spherical is not None and len(spherical) == len(pairs):
            shells = [(at, l, bool(p)) for (at, l), p in zip(pairs, spherical)]
            if any(p for p in spherical):
                self.symmetry_metadata['spherical_order_assumed'] = \
                    'cca_m_ascending'
        elif sum(_cartesian_shell_size(l) for _, l in pairs) == nbf:
            shells = [(at, l, False) for at, l in pairs]
        elif sum(_spherical_shell_size(l) for _, l in pairs) == nbf:
            # Pure spherical-harmonic basis (ISPHER=1). The component order
            # is assumed to be CCA/libint (m = -l..+l); record the
            # assumption so runs are auditable until the runtime spherical
            # path is validated end-to-end.
            # OpenQP only makes l >= 2 spherical: source/integrals/cart2sph.F90
            # c2s_ncomp transforms a shell when (pure == 1 .and. l >= 2), so s
            # and p AOs stay Cartesian (x, y, z) even in a spherical basis.
            # Tagging them pure applied the spherical component permutation
            # [1,2,0] to p shells, which produced WRONG irrep labels rather
            # than 'mixed' ones -- e.g. under C2z a p_z-dominant MO picked up
            # chi = -1 where it must be +1.  Those labels feed
            # stage_response_symmetry -> OQP::sym_pair_irrep -> the Davidson
            # seeding, so the error was silent and landed on exactly the
            # spherical families (cc-pVXZ, def2) most production runs use.
            # The dimension test above is unaffected: cart and spherical sizes
            # agree for l <= 1.
            shells = [(at, l, l >= 2) for at, l in pairs]
            self.symmetry_metadata['spherical_order_assumed'] = 'cca_m_ascending'
        else:
            return None, None, 0, 'skipped_unrecognized_basis_dimension'

        # Unpack triangular overlap to a square symmetric matrix.
        packed = np.asarray(self.data['OQP::SM'], dtype=float).ravel()
        if packed.size != nbf * (nbf + 1) // 2:
            return None, None, 0, 'skipped_overlap_shape_mismatch'
        smat = np.zeros((nbf, nbf))
        rows, cols = np.tril_indices(nbf)
        smat[rows, cols] = packed
        smat[cols, rows] = packed
        return shells, smat, nbf, None

    def _overlap_square(self, nbf):
        """Square symmetric overlap from OQP::SM, or None if unavailable."""
        try:
            packed = np.asarray(self.data['OQP::SM'], dtype=float).ravel()
            if packed.size != nbf * (nbf + 1) // 2:
                return None
            smat = np.zeros((nbf, nbf))
            rows, cols = np.tril_indices(nbf)
            smat[rows, cols] = packed
            smat[cols, rows] = packed
            return smat
        except Exception:
            return None

    def _mo_coefficients(self, tag, nbf):
        """AO x MO coefficient matrix from a Fortran column-stored tag."""
        # Fortran stores MOs column-wise; the C-order view has MOs as rows,
        # so transpose to (n_ao, n_mo).
        return np.asarray(self.data[tag], dtype=float).reshape(nbf, nbf).T

    def _labelling_overlap_deviation(self, shells, smat, detection):
        """max |T^T S T - S| over the symmetry operations, or None if it
        cannot be formed.

        Built from ``matrix_input_frame`` -- the SAME operator
        ``assign_mo_irreps`` labels with -- because that is the only thing this
        guard is entitled to validate. The overlap it is tested against is the
        run's real AO overlap, which lives in the INPUT frame: nothing
        reorients the geometry unless ``use_integral_symmetry`` is set, and
        that is off by default.

        The first version used ``build_reduction_maps``, which hard-codes
        ``op['matrix']`` -- the sign-diagonal operation in the symmetry
        STANDARD orientation. Against an input-frame overlap that is O(1)
        wrong for any molecule not already in the standard frame, so the guard
        fired, labelling was refused, ``OQP::sym_pair_irrep`` was never staged
        and the Davidson coverage fix silently did nothing. Measured over the
        shipped decks: 211 of the 241 with a non-C1 abelian subgroup, including
        every H2O deck. It looked correct only because the two geometries it
        was validated on happen to be in the standard frame already.

        The twin guard in ``stage_integral_symmetry_maps`` may keep using the
        standard-frame maps: that path returns early unless the geometry has
        actually been reoriented, so there the two frames coincide.
        """
        try:
            import numpy as np
            from oqp.library.symmetry import _ao_operator_matrix, _normalize_shells

            if smat is None:
                return None
            smat = np.asarray(smat, dtype=float)
            norm_shells = _normalize_shells(shells)
            worst = 0.0
            for op in detection['operations']:
                transform = _ao_operator_matrix(norm_shells, op,
                                                matrix_key='matrix_input_frame')
                if transform.shape[0] != smat.shape[0]:
                    return None
                # Functions transform with T as columns, so invariance of the
                # metric reads T^T S T = S -- not T S T^T, which is the
                # signed-permutation convention the integral path uses.
                worst = max(worst, float(np.max(np.abs(
                    transform.T @ smat @ transform - smat))))
            return worst
        except Exception:
            # A guard that cannot be evaluated must not block labelling.
            return None

    def label_molecular_orbitals(self):
        """Assign abelian irrep labels to converged MOs (metadata only, non-fatal).

        Stores the result under ``symmetry_metadata['mo_labels']``; never
        changes SCF/integral/response behavior.
        """
        meta = self.symmetry_metadata
        if not meta or meta.get('status', 'disabled') == 'disabled':
            return None
        if not meta.get('label_mo', True):
            return None
        detection = meta.get('detection')
        if not detection:
            return None

        try:
            from oqp.library.symmetry import assign_mo_irreps

            shells, smat, nbf, skip_reason = self._symmetry_labeling_inputs()
            if shells is None:
                meta['mo_labels'] = {'status': skip_reason}
                return None

            # Defense in depth, mirroring the integral path: the operator
            # these labels are computed from must leave the real overlap matrix
            # invariant, T^T S T = S.  Any shell-convention or frame mistake
            # shows up here.
            #
            # T^T S T, not T S T^T: T is metric-orthogonal but not orthogonal
            # once a spherical shell with l >= 2 mixes components under a
            # rotation. The two forms coincide for a signed permutation, which
            # is why a test restricted to axis-aligned frames -- or to
            # Cartesian d shells -- cannot tell them apart. Verified where they
            # DO differ: water in a generic three-angle rotated frame labels
            # identically to the standard frame under 6-31G* (Cartesian d,
            # control), cc-pVDZ (spherical d) and cc-pVTZ (spherical d and f),
            # with no orbital coming back 'mixed'.  This is not hypothetical -- tagging s and p
            # shells as spherical (fixed in the previous commit) produced
            # confidently WRONG p-shell signs, and this check is what would
            # have caught it.  Labels are metadata, so a failure records the
            # reason and declines to label rather than aborting the run; the
            # consumers (stage_response_symmetry -> OQP::sym_pair_irrep ->
            # the Davidson seeding) then stay inert, which is the safe state.
            deviation = self._labelling_overlap_deviation(shells, smat,
                                                          detection)
            if deviation is not None and deviation > 1.0e-6:
                meta['mo_labels'] = {
                    'status': 'skipped_overlap_invariance',
                    'deviation': deviation,
                }
                return None

            tolerance = float(meta.get('tolerance', 1.0e-5))
            result = {'status': 'ok'}
            spins = [('alpha', 'OQP::VEC_MO_A')]
            if self.config.get('scf', {}).get('type', 'rhf') != 'rhf':
                spins.append(('beta', 'OQP::VEC_MO_B'))
            for spin, tag in spins:
                coefficients = self._mo_coefficients(tag, nbf)
                result[spin] = assign_mo_irreps(
                    coefficients, smat, shells,
                    detection['operations'], detection['character_table'],
                    tolerance=max(tolerance, 1.0e-4),
                    matrix_key='matrix_input_frame',
                )
            state = self._label_scf_state(result)
            if state is not None:
                result['scf_state'] = state
            meta['mo_labels'] = result
            # Gate A of the reductions plan: shell/AO symmetry maps
            # (metadata only; consumed by future petite-list code).
            try:
                from oqp.library.symmetry import build_reduction_maps
                meta['reduction_maps'] = build_reduction_maps(
                    shells, detection['operations'])
            except Exception as exc:
                meta['reduction_maps'] = {'status': 'error', 'error': str(exc)}
            try:
                self._dump_mo_labels_log(result)
            except Exception:
                pass
            return result
        except Exception as exc:
            # Labeling must never break the run in the metadata-only phase.
            meta['mo_labels'] = {'status': 'error', 'error': str(exc)}
            return None

    def _label_scf_state(self, mo_label_result):
        """Total-symmetry label of the SCF determinant (metadata only).

        For an abelian group the state irrep is the direct product of the
        occupied MO irreps (closed pairs cancel, so only singly occupied
        orbitals contribute for ROHF/UHF). Returns None when occupations
        are unavailable.
        """
        try:
            from oqp.library.symmetry import product_irrep

            table = self.symmetry_metadata['detection']['character_table']
            try:
                na = int(np.asarray(self.data['nelec_A']).ravel()[0])
                nb = int(np.asarray(self.data['nelec_B']).ravel()[0])
            except Exception:
                nocc = int(np.asarray(self.data['nocc']).ravel()[0])
                na = nb = nocc
            if na <= 0 or nb < 0:
                return None

            alpha = mo_label_result['alpha']['labels']
            beta = mo_label_result.get('beta', mo_label_result['alpha'])['labels']
            if na > len(alpha) or nb > len(beta):
                return None

            irrep = product_irrep(list(alpha[:na]) + list(beta[:nb]), table)
            state = {'irrep': irrep, 'nelec_alpha': na, 'nelec_beta': nb}
            multiplicity = getattr(self, 'mult', None)
            if multiplicity:
                state['multiplicity'] = int(multiplicity)
                state['term'] = f"{int(multiplicity)}{irrep.upper()}"
            else:
                state['term'] = irrep.upper()
            return state
        except Exception:
            return None

    def reorient_for_integral_symmetry(self):
        """OpenQP reorientation to the standard frame (geometry only).

        Call before the guess/basis stage; ``stage_integral_symmetry_maps``
        completes the activation once the basis is available. No-op unless
        ``[symmetry] use_integral_symmetry`` is enabled.
        """
        meta = self.symmetry_metadata
        if not meta or not meta.get('use_integral_symmetry'):
            return False
        detection = meta.get('detection')
        if not detection:
            return False

        # [symmetry] move_to_standard_frame = false: do the reduction where the
        # molecule already is.
        #
        # The move exists to make the AO-level operator a SIGNED PERMUTATION.
        # In the standard frame each symmetry operation sends an AO to plus or
        # minus one other AO, which is a cheap scatter; in a tilted frame the
        # same operation mixes components within a shell and the operator is a
        # dense block. That is a performance choice, not a correctness
        # requirement -- and it is the reason build_reduction_maps refuses a
        # dense matrix outright.
        #
        # What the reduction actually needs is the SHELL PERMUTATION, and that
        # is frame-independent: which shell maps to which is a property of the
        # atom permutation. petite_quartet_weight reads nothing else.
        #
        # So the no-move path keeps the same shell map and stages dense
        # input-frame operator blocks for the projection instead of signed
        # permutations. Everything the move was introduced to avoid goes away
        # with it: no back-transform of outputs, no runtype allow-list, no C1
        # molecule translated several bohr for a reduction it never gets.
        if not self._parse_bool_like(
                self.config.get('symmetry', {}).get(
                    'move_to_standard_frame', True)):
            meta['integral_symmetry'] = {'status': 'input_frame'}
            return True

        # Geometry-displacing drivers (optimizers, numerical Hessians, MEP,
        # NEB, ...) must not have the frame rotated under them; the petite
        # reduction is restricted to single-point runtypes for now.
        runtype = str(self.config.get('input', {}).get('runtype', 'energy')).lower()
        if runtype not in ('energy', 'grad', 'prop', 'properties'):
            meta['integral_symmetry'] = {'status': f'skipped_runtype_{runtype}'}
            return False

        try:
            from oqp.library.symmetry_detect import attach_detection_metadata

            # Reorient.
            #
            # Detection of a rotated geometry may legitimately pick a
            # different but equivalent standard frame (degenerate axis
            # choices, e.g. the three C2 axes of d2h), so a single
            # rotate-then-redetect is NOT guaranteed to converge. Iterate
            # until detection of the current geometry returns the identity
            # frame -- only then are the stored operations valid as-is.
            atoms = np.asarray(self.get_atoms(), dtype=float).ravel()
            coords = np.asarray(self.get_system(), dtype=float).reshape(-1, 3)
            input_coords = coords.copy()
            total_rotation = np.eye(3)
            total_origin = np.zeros(3)
            converged = False
            for _ in range(4):
                attach_detection_metadata(meta, atoms, coords)
                detection = meta['detection']
                origin = np.asarray(detection['origin'], dtype=float)
                rotation = np.asarray(detection['orientation'], dtype=float)
                if (np.max(np.abs(rotation - np.eye(3))) < 1.0e-12
                        and np.max(np.abs(origin)) < 1.0e-10):
                    converged = True
                    break
                total_origin = total_origin + total_rotation.T @ origin
                total_rotation = rotation @ total_rotation
                coords = (coords - origin) @ rotation.T
                self.update_system(coords.ravel())
            if not converged:
                meta['integral_symmetry'] = {'status': 'skipped_orientation_not_converged'}
                return False

            # OpenQP contract: ALL outputs (geometry, gradients, MOs)
            # are consistently in the standard orientation. The transform
            # below maps user input axes to it:
            #   r_std = (r_input - origin) @ rotation^T
            # so input-frame vectors are recovered via v_input = v_std @ R.
            meta['integral_symmetry'] = {
                'status': 'reoriented',
                'input_to_standard': {
                    'rotation': total_rotation.tolist(),
                    'origin': total_origin.tolist(),
                },
            }
            return True
        except Exception as exc:
            meta['integral_symmetry'] = {'status': 'error', 'error': str(exc)}
            return False

    def stage_integral_symmetry_maps(self):
        """Stage petite-list maps for the Fortran SCF (requires the basis).

        Fail-safe: any inconsistency leaves the run on the C1 path with the
        reason recorded in the metadata -- and, unlike before, printed. The
        fallback produces an identical energy, so a silent one is invisible:
        the reduction was disabled on every spherical basis for as long as it
        existed and no output said so.
        """
        meta = self.symmetry_metadata
        if not meta or not meta.get('use_integral_symmetry'):
            return False
        staged = self._stage_integral_symmetry_maps(meta)
        if not staged and not meta.get('integral_symmetry', {}).get('status', '').startswith('skipped'):
            meta.setdefault('integral_symmetry', {})
            meta['integral_symmetry'].setdefault('status', 'skipped_unknown')
        try:
            self._dump_symmetry_log()
        except Exception:
            pass
        return staged

    def _stage_integral_symmetry_maps(self, meta):
        """Body of :meth:`stage_integral_symmetry_maps`; records a status and
        returns True only when the reduction is actually active."""
        detection = meta.get('detection')
        if not detection:
            meta['integral_symmetry'] = {'status': 'skipped_no_detection'}
            return False
        # reorient_for_integral_symmetry records its own reason (a runtype the
        # allow-list rejects, a frame that would not converge, ...); keep it.
        frame_status = meta.get('integral_symmetry', {}).get('status')
        if frame_status not in ('reoriented', 'input_frame'):
            return False
        # No move: the AO-level operator is dense, so the signed-permutation
        # maps cannot be built. The shell map is the same either way.
        input_frame = frame_status == 'input_frame'

        try:
            from oqp.library.symmetry import build_reduction_maps

            basis = self.data.get_basis()
            if not basis:
                meta['integral_symmetry'] = {'status': 'skipped_no_basis'}
                return False
            # Per-shell sphericity comes from the library (OQP::basis
            # 'spherical'), which reports the EFFECTIVE flag including the
            # l >= 2 rule. Passing (atom, l) pairs without it makes
            # _normalize_shells default to pure=False, so n_ao is the
            # Cartesian count and every spherical basis -- cc-pVXZ, def2 and
            # the whole 6-311G family -- failed the check below and silently
            # ran C1. Do not reintroduce a dimension-matching guess here: the
            # Cartesian and spherical sizes agree for l <= 1, so the AO total
            # cannot distinguish "all shells pure" (wrong: it applies the
            # spherical component order to p shells) from OpenQP's real
            # convention.
            spherical = basis.get('spherical')
            if spherical is None:
                meta['integral_symmetry'] = {'status': 'skipped_no_shell_purity'}
                return False
            shells = [(int(at), int(l), bool(p)) for at, l, p
                      in zip(basis['centers'], basis['angs'], spherical)]
            if input_frame:
                # Same shell permutation, dense input-frame blocks in place of
                # the signed AO maps. symmetrize_skeleton_fock already prefers
                # the blocked projector whenever OQP::sym_op_blocks is staged,
                # so nothing downstream changes.
                from oqp.library.symmetry import build_full_group_blocks
                dense = build_full_group_blocks(
                    shells, detection['operations'],
                    matrix_key='matrix_input_frame')
                if int(dense['n_ao']) != int(basis['nbf']):
                    meta['integral_symmetry'] = {
                        'status': 'skipped_basis_mismatch',
                        'n_ao': int(dense['n_ao']),
                        'nbf': int(basis['nbf']),
                    }
                    return False
                self.data['OQP::sym_shell_map'] = \
                    (np.asarray(dense['shell_permutation'], dtype=np.int64) + 1).ravel()
                self.data['OQP::sym_op_blocks'] = \
                    np.asarray(dense['blocks'], dtype=np.float64)
                # Dense operator, abelian group: screen as the standard-frame
                # abelian path does, not as the non-abelian tier.
                self.data['OQP::sym_nonabelian'] = np.array([0], dtype=np.int64)
                self.data['OQP::sym_petite_enable'] = np.array([1], dtype=np.int64)
                meta['reduction_maps'] = {
                    'n_operations': int(dense['n_operations']),
                    'n_ao': int(dense['n_ao']),
                }
                meta['integral_symmetry'] = {
                    'status': 'active',
                    'group': meta.get('subgroup'),
                    'n_operations': int(dense['n_operations']),
                    'full_group': False,
                    'reoriented': False,
                    'frame': 'input',
                }
                try:
                    self._dump_symmetry_log()
                except Exception:
                    pass
                return True

            maps = build_reduction_maps(shells, detection['operations'])
            if maps['n_ao'] != int(basis['nbf']):
                meta['integral_symmetry'] = {
                    'status': 'skipped_basis_mismatch',
                    'n_ao': int(maps['n_ao']),
                    'nbf': int(basis['nbf']),
                }
                return False

            # Defense-in-depth: the maps must leave the real overlap matrix
            # invariant (T S T^T = S); any frame/staging inconsistency shows
            # up here and falls back to C1.
            smat = self._overlap_square(int(basis['nbf']))
            if smat is not None:
                identity = np.arange(maps['n_ao'])
                for iop in range(maps['n_operations']):
                    transform = np.zeros((maps['n_ao'], maps['n_ao']))
                    transform[np.array(maps['ao_target'][iop]), identity] = \
                        np.array(maps['ao_sign'][iop], dtype=float)
                    deviation = float(np.max(np.abs(
                        transform @ smat @ transform.T - smat)))
                    if deviation > 1.0e-6:
                        meta['integral_symmetry'] = {
                            'status': 'skipped_overlap_invariance',
                            'operation': maps['operation_names'][iop],
                            'deviation': deviation,
                        }
                        return False

            # 1-based flat maps for the Fortran consumers (op-major,
            # shell/AO index fastest -- see load_petite_list /
            # symmetrize_skeleton_fock).
            self.data['OQP::sym_shell_map'] = \
                (np.asarray(maps['shell_permutation'], dtype=np.int64) + 1).ravel()
            self.data['OQP::sym_ao_target'] = \
                (np.asarray(maps['ao_target'], dtype=np.int64) + 1).ravel()
            self.data['OQP::sym_ao_sign'] = \
                np.asarray(maps['ao_sign'], dtype=np.float64).ravel()
            # Per-atom orbit weights for the XC grid reduction: orbit size
            # for the unique (lowest-index) atom of each orbit, zero for
            # its images.
            permutations = np.array(
                [op['permutation'] for op in detection['operations']], dtype=int)
            representative = permutations.min(axis=0)
            natom = permutations.shape[1]
            atom_weight = np.zeros(natom)
            for atom in range(natom):
                if representative[atom] == atom:
                    atom_weight[atom] = float(np.count_nonzero(representative == atom))
            self.data['OQP::sym_atom_weight'] = atom_weight
            self._sym_atom_weight_ops = detection['operations']

            # Non-abelian upgrade: if the FULL point group is larger than
            # the abelian subgroup, stage its shell map and dense per-shell
            # operation blocks instead -- the petite filter then keeps
            # 1/|G_full| of the quartets and the skeleton is symmetrized
            # with the block transforms (e.g. benzene: 24 ops vs 8).
            # Tier selection: 'true' = abelian subgroup (machine-exact,
            # validated to ~1e-12); 'full' = full point group (up to |G|
            # quartet reduction, e.g. 24 for D6h benzene, accurate to
            # ~1e-7 -- a residual kernel-threshold asymmetry between
            # non-abelian orbit members is still under investigation).
            want_full = str(self.config.get('symmetry', {})
                            .get('use_integral_symmetry', '')).strip().lower() == 'full'

            # The full (non-abelian) tier is incompatible with DFT, and not by
            # an oversight that can be patched.
            #
            # It reduces and projects the two halves of the Fock with DIFFERENT
            # groups: J/K with the full group, Vxc with the abelian subgroup.
            # The latter is deliberate -- Lebedev angular grids are invariant
            # under the axis-aligned octahedral operations but NOT under C3/C6.
            # So the JK half is forced exactly symmetric under the full group
            # while the XC half can only ever be abelian-symmetric, and the SCF
            # converges to a density that satisfies neither.
            #
            # Measured on benzene 6-31G* at an exact D6h geometry, |G| = 24,
            # against the C1 reference:
            #
            #     pure HF                            dE = 5.9e-09
            #     hfscale = 1.0 with a BLYP Vxc      dE = 1.6e-03
            #     BHHLYP                             dE = 3.1e-04
            #     BLYP                               dE = 1.9e-03
            #
            # The second row is the discriminating one: exchange is treated
            # bit-identically to the case that works, and merely putting a Vxc
            # into the Fock costs six orders. So the variable is the presence
            # of an abelian-projected Vxc beside a full-group-projected JK, not
            # the exchange fraction.
            #
            # All three repairs are closed. Projecting both halves with the
            # full group is impossible while the grid is not C3/C6 invariant --
            # and note this is the GRID, not the grid reduction, so withdrawing
            # the XC atom weights does not help (measured: dE unchanged to
            # every digit). Projecting both halves with the abelian group is
            # wrong because the JK skeleton was reduced with 24 operations
            # (measured: dE = 7.2e+05). Reducing JK with 8 operations IS the
            # abelian tier, which is already exact (dE = 0 for both a hybrid
            # and a pure GGA) and already delivers 3.4-3.7x.
            #
            # So: fall back to the abelian tier and say so.
            #
            # This is not free, and the first version of this comment wrongly
            # said it was. Measured on benzene HF/cc-pVTZ, whole-job wall
            # clock: C1 41.6 s, abelian 22.2 s (1.87x), full 11.6 s (3.59x).
            # The 3.6x quoted for this work is the FULL tier; abelian is 1.9x.
            # So DFT gives up about a factor of 1.9 here. The decision stands
            # anyway -- the full tier produces a WRONG DFT energy and a wrong
            # number is not a tradeoff against speed -- but a maintainer
            # weighing a symmetry-adapted DFT grid should see the real size of
            # the prize. Joint finding with @karmachoi.
            if want_full and str(self.config.get('input', {})
                                 .get('functional', '')).strip():
                want_full = False
                meta['full_group_declined'] = 'dft_grid_not_invariant'

            full_group = False
            full_ops = None
            try:
                if not want_full:
                    raise StopIteration  # stay on the exact abelian tier
                from oqp.library.symmetry import build_full_group_blocks
                from oqp.library.symmetry import _ao_operator_matrix, _normalize_shells
                from oqp.library.symmetry_detect import enumerate_full_group

                atoms_arr = np.asarray(self.get_atoms(), dtype=float).ravel()
                coords_arr = np.asarray(self.get_system(), dtype=float).reshape(-1, 3)
                tolerance = float(meta.get('tolerance', 1.0e-5))
                full_ops = enumerate_full_group(atoms_arr, coords_arr,
                                                tolerance=tolerance)
                if len(full_ops) > maps['n_operations']:
                    full = build_full_group_blocks(shells, full_ops)
                    # Same defense-in-depth as the abelian path: every
                    # operation must leave the real overlap invariant.
                    ok = True
                    if smat is not None:
                        norm_shells = _normalize_shells(shells)
                        for op in full_ops:
                            transform = _ao_operator_matrix(norm_shells, op)
                            # Operator-transform side: functions transform
                            # with T as columns, so S-invariance is
                            # T^T S T = S (T is metric-orthogonal, not
                            # orthogonal, once d shells mix under rotations).
                            deviation = float(np.max(np.abs(
                                transform.T @ smat @ transform - smat)))
                            # Tight gate: ~1e-7 overlap residuals (geometry
                            # symmetric only to input precision) become
                            # 1e-5-level Fock errors that SCF amplifies, so
                            # fall back to the exact abelian path instead.
                            if deviation > 1.0e-8:
                                ok = False
                                break
                    if ok:
                        self.data['OQP::sym_shell_map'] = \
                            (np.asarray(full['shell_permutation'],
                                        dtype=np.int64) + 1).ravel()
                        self.data['OQP::sym_op_blocks'] = \
                            np.asarray(full['blocks'], dtype=np.float64)
                        # NOTE: XC atom weights intentionally stay with the
                        # abelian (sign-operation) group: Lebedev angular
                        # grids are invariant under the axis-aligned
                        # octahedral operations but NOT under C3/C6
                        # rotations, so full-group grid reduction would be
                        # inexact.
                        self.data['OQP::sym_nonabelian'] = \
                            np.array([1], dtype=np.int64)
                        full_group = True
                        meta['reduction_maps_full'] = {
                            'n_operations': full['n_operations'],
                            'operations': full_ops,
                        }
            except Exception:
                full_group = False

            self.data['OQP::sym_petite_enable'] = np.array([1], dtype=np.int64)

            meta['reduction_maps'] = maps
            input_to_standard = meta.get('integral_symmetry', {}).get('input_to_standard')
            meta['integral_symmetry'] = {
                'status': 'active',
                'group': meta.get('subgroup'),
                'n_operations': (meta['reduction_maps_full']['n_operations']
                                 if full_group else maps['n_operations']),
                'full_group': full_group,
                'reoriented': True,
                'input_to_standard': input_to_standard,
            }
            return True
        except Exception as exc:
            # Fail safe to the C1 path.
            meta['integral_symmetry'] = {'status': 'error', 'error': str(exc)}
            try:
                self.data['OQP::sym_petite_enable'] = np.array([0], dtype=np.int64)
            except Exception:
                pass
            return False

    def _clear_response_symmetry_tags(self):
        """Drop any previously staged response-symmetry tables.

        See stage_response_symmetry: the tag store survives across steps of a
        multi-geometry job, so a stale table is worse than none.
        """
        for tag in ('OQP::sym_pair_irrep', 'OQP::sym_response_project_enable'):
            try:
                del self.data[tag]
            except Exception:
                pass

    def stage_response_symmetry(self):
        """Stage per-pair irrep indices for response-space blocking.

        Builds OQP::sym_pair_irrep (1-based irrep index per excitation
        pair, occupied index fastest) from the converged MO labels. Only
        acts when ``use_response_symmetry`` is enabled; bails to the
        unblocked solver on any 'mixed' orbital or inconsistency.
        """
        meta = self.symmetry_metadata
        # Invalidate FIRST, before any bail can return. The table is consumed
        # by the Fortran Davidson guess, and the tag store outlives a single
        # step: an optimiser, IRC or NAMD run that staged a table at step N and
        # bails at step N+1 -- lower symmetry, a 'mixed' orbital, labels
        # refused -- would otherwise leave the previous step's table in place
        # and seed the guess from irreps that belong to a different geometry.
        # Nothing downstream can detect that: the table is well formed, just
        # stale. Deleting it makes a bail mean "no table", which mrinivec
        # already handles by falling back to the historical seeds.
        self._clear_response_symmetry_tags()

        if not meta:
            return False
        detection = meta.get('detection')
        if not detection:
            return False

        td_type = str(self.config.get('tdhf', {}).get('type', '')).lower()
        if td_type not in ('tda', 'rpa', 'sf', 'mrsf'):
            if td_type:
                meta['response_symmetry'] = {
                    'status': f'skipped_unsupported_td_type_{td_type}'}
            return False

        try:
            from oqp.library.symmetry import product_irrep

            mo_labels = meta.get('mo_labels')
            if not mo_labels or mo_labels.get('status') != 'ok':
                mo_labels = self.label_molecular_orbitals()
            if not mo_labels or mo_labels.get('status') != 'ok':
                meta['response_symmetry'] = {'status': 'skipped_no_mo_labels'}
                return False

            na = int(np.asarray(self.data['nelec_A']).ravel()[0])
            nb = int(np.asarray(self.data['nelec_B']).ravel()[0])
            nbf = len(mo_labels['alpha']['labels'])

            if td_type in ('sf', 'mrsf'):
                occ_labels = mo_labels['alpha']['labels'][:na]
                vir_labels = mo_labels.get('beta', mo_labels['alpha'])['labels'][nb:]
            else:
                occ_labels = mo_labels['alpha']['labels'][:na]
                vir_labels = mo_labels['alpha']['labels'][na:]

            if 'mixed' in occ_labels or 'mixed' in vir_labels:
                meta['response_symmetry'] = {'status': 'skipped_mixed_orbitals'}
                return False

            table = detection['character_table']
            irreps = list(table.keys())
            # Fortran xvec layout: occupied index fastest within each virtual.
            pair_irrep = np.zeros(len(occ_labels)*len(vir_labels), dtype=np.int64)
            idx = 0
            for vir in vir_labels:
                for occ in occ_labels:
                    label = product_irrep([occ, vir], table)
                    if label == 'mixed':
                        meta['response_symmetry'] = {'status': 'skipped_mixed_pair'}
                        return False
                    pair_irrep[idx] = irreps.index(label) + 1
                    idx += 1

            if td_type == 'mrsf' and na - nb >= 2:
                # The three open-open slots do not carry the plain
                # Gamma_i (x) Gamma_a label of their (i, a) indices, because
                # MRSF stores that sector spin-adapted and folded (SI Fig. S2
                # and Sec. 8 of JCP 149, 104101).  With O1 = na-2, O2 = na-1
                # (0-based) and the slot convention "annihilate alpha in i,
                # create beta in a" acting on |O1a O2a>:
                #   (O1,O1) holds the folded (L -/+ R)/sqrt2, whose two
                #     determinants are both the open-shell {O1,O2} pair, so it
                #     transforms as Gamma_O1 (x) Gamma_O2;
                #   (O2,O1) is G = |O1 O1> and (O1,O2) is D = |O2 O2>, both
                #     doubly occupied and therefore totally symmetric.
                # Labelling them by the raw index product puts all three in the
                # wrong blocks, which defeats the coverage the Davidson guess
                # and the response projection both rely on.
                o1, o2 = na - 2, na - 1
                v_o1, v_o2 = o1 - nb, o2 - nb
                n_occ = len(occ_labels)
                if 0 <= v_o1 and v_o2 < len(vir_labels):
                    oo_label = product_irrep(
                        [occ_labels[o1], vir_labels[v_o2]], table)
                    if oo_label != 'mixed':
                        pair_irrep[o1 + v_o1*n_occ] = irreps.index(oo_label) + 1
                    # irreps[0] is the totally symmetric representation
                    pair_irrep[o2 + v_o1*n_occ] = 1
                    pair_irrep[o1 + v_o2*n_occ] = 1

            self.data['OQP::sym_pair_irrep'] = pair_irrep
            # The table itself is staged whenever it can be built: the Davidson
            # initial trial vectors need it to guarantee that every irreducible
            # representation is represented, without which an unseeded irrep
            # contributes no roots and every state index shifts.  Confining
            # residuals to a dominant irrep is a separate, experimental
            # behaviour and stays behind use_response_symmetry.
            self.data['OQP::sym_response_project_enable'] = np.array(
                [1 if meta.get('use_response_symmetry') else 0], dtype=np.int64)
            meta['response_symmetry'] = {
                'status': 'active',
                'td_type': td_type,
                'n_pairs': int(pair_irrep.size),
                'irreps': irreps,
            }
            return True
        except Exception as exc:
            meta['response_symmetry'] = {'status': 'error', 'error': str(exc)}
            return False

    def label_excited_states(self):
        """Assign abelian irrep labels to TD excited states (metadata only).

        Supports tda/rpa (closed-shell, occ/vir from VEC_MO_A) and sf/mrsf
        (occ from alpha, vir from beta MOs; total symmetry includes the
        direct product of the reference SOMO irreps). Stores results under
        ``symmetry_metadata['state_labels']``; never fatal.
        """
        meta = self.symmetry_metadata
        if not meta or meta.get('status', 'disabled') == 'disabled':
            return None
        if not meta.get('label_states', True):
            return None
        detection = meta.get('detection')
        if not detection:
            return None

        td_type = str(self.config.get('tdhf', {}).get('type', '')).lower()
        if td_type not in ('tda', 'rpa', 'sf', 'mrsf'):
            if td_type:
                meta['state_labels'] = {'status': f'skipped_unsupported_td_type_{td_type}'}
            return None

        try:
            from oqp.library.symmetry import assign_state_irreps

            shells, smat, nbf, skip_reason = self._symmetry_labeling_inputs()
            if shells is None:
                meta['state_labels'] = {'status': skip_reason}
                return None

            # MO labels provide the SOMO reference product for sf/mrsf.
            mo_labels = meta.get('mo_labels')
            if not mo_labels or mo_labels.get('status') != 'ok':
                mo_labels = self.label_molecular_orbitals()
            if not mo_labels or mo_labels.get('status') != 'ok':
                meta['state_labels'] = {'status': 'skipped_no_mo_labels'}
                return None

            na = int(np.asarray(self.data['nelec_A']).ravel()[0])
            nb = int(np.asarray(self.data['nelec_B']).ravel()[0])

            c_alpha = self._mo_coefficients('OQP::VEC_MO_A', nbf)
            if td_type in ('sf', 'mrsf'):
                # Spin-flip: occupied alpha -> virtual beta.
                c_beta = self._mo_coefficients('OQP::VEC_MO_B', nbf)
                occ, vir = c_alpha[:, :na], c_beta[:, nb:]
                reference_labels = mo_labels['alpha']['labels'][nb:na]
            else:
                occ, vir = c_alpha[:, :na], c_alpha[:, na:]
                reference_labels = []
            n_occ, n_vir = occ.shape[1], vir.shape[1]

            # Fortran bvec(xvec_dim, nstates), occupied index fastest; the
            # C-order buffer is state-major.
            bvec = np.asarray(self.data['OQP::td_bvec_mo'], dtype=float).ravel()
            xvec_dim = n_occ * n_vir
            if xvec_dim == 0 or bvec.size % xvec_dim != 0:
                meta['state_labels'] = {'status': 'skipped_amplitude_shape_mismatch'}
                return None
            nstates = bvec.size // xvec_dim
            amplitudes = bvec.reshape(nstates, n_vir, n_occ).transpose(0, 2, 1)

            tolerance = float(meta.get('tolerance', 1.0e-5))
            result = assign_state_irreps(
                amplitudes, occ, vir, smat, shells,
                detection['operations'], detection['character_table'],
                reference_labels=reference_labels,
                tolerance=max(tolerance, 1.0e-3),
                matrix_key='matrix_input_frame',
            )
            result = dict(result)
            result['status'] = 'ok'
            result['td_type'] = td_type
            multiplicity = self.config.get('tdhf', {}).get('mult')
            if multiplicity:
                result['terms'] = [
                    f"{int(multiplicity)}{lbl.upper()}" for lbl in result['labels']
                ]
            meta['state_labels'] = result
            try:
                self._dump_state_labels_log(result)
            except Exception:
                pass
            return result
        except Exception as exc:
            # Labeling must never break the run in the metadata-only phase.
            meta['state_labels'] = {'status': 'error', 'error': str(exc)}
            return None

    @mpi_dump
    def _dump_state_labels_log(self, result):
        """Append excited-state irrep labels to the main log (best effort)."""
        try:
            lines = [
                '',
                '   ==============================================',
                '   PyOQP: excited-state symmetry labels (metadata only)',
                '   ==============================================',
            ]
            try:
                energies = np.asarray(self.data['OQP::td_energies'], dtype=float).ravel()
            except Exception:
                energies = None
            terms = result.get('terms') or [lbl.upper() for lbl in result['labels']]
            for istate, term in enumerate(terms):
                state = (public_state_label(self.config, istate + 1)
                         if is_mrsf(self.config) else 'state %3d' % (istate + 1))
                if energies is not None and istate < energies.size:
                    lines.append(f'   {state:10s}  {energies[istate]:14.8f}  {term}')
                else:
                    lines.append(f'   {state:10s}  {"":14s}  {term}')
            lines.append('')
            with open(self.log, 'a', encoding='utf-8') as fout:
                fout.write('\n'.join(lines))
        except Exception:
            pass

    def symmetrize_gradient(self, grads):
        """Project gradients onto the totally symmetric component.

        g'_a = (1/|G|) sum_op M_op^T g_{perm_op(a)}, exact for the skeleton
        two-electron gradient and a noise-cleaner for the rest.

        The operator MUST be expressed in the frame the gradient lives in.
        ``op['matrix']`` is the sign-diagonal operation in the standard
        orientation; ``op['matrix_input_frame']`` is the same operation
        conjugated back to the input frame. Both are orthogonal 3x3 Cartesian
        matrices, so the formula is unchanged -- only the frame differs, and
        picking the wrong one is silent.

        This used ``op['matrix']`` unconditionally, which was correct while the
        reduction always reoriented the molecule and became wrong the moment
        ``move_to_standard_frame=false`` existed. On the shipped
        ``h2o_rhf_6-31g_hf`` deck -- C2v water lying on a diagonal -- the
        standard-frame operator zeroes the x component of every gradient it
        touches: |g_H| came out 0.00989131 against a reference 0.01314967, a
        25% shortfall. Not a rotation: the norm itself is wrong. Energies were
        unaffected, which is why an energy-only verification missed it, and an
        already-aligned test molecule cannot detect it at all because the two
        frames coincide there.
        """
        meta = self.symmetry_metadata
        if not meta or meta.get('integral_symmetry', {}).get('status') != 'active':
            return grads
        detection = meta.get('detection')
        if not detection:
            return grads

        try:
            operations = detection['operations']
            full = meta.get('reduction_maps_full')
            if meta.get('integral_symmetry', {}).get('full_group') and full:
                operations = full['operations']
            # Which frame is the staged reduction operating in? The staging
            # records 'input' only on the no-move path; everything else moved
            # the molecule first. reduction_maps_full is built exclusively on
            # the standard-frame branch, so the two never combine.
            matrix_key = 'matrix' \
                if meta.get('integral_symmetry', {}).get('frame') != 'input' \
                else 'matrix_input_frame'
            if any(matrix_key not in op for op in operations):
                # Refuse rather than project with an operator from the wrong
                # frame. Returning the skeleton gradient unprojected keeps a
                # small symmetry-breaking residual; projecting with the wrong
                # operator deletes whole Cartesian components.
                return grads
            arr = np.asarray(grads, dtype=float)
            shape = arr.shape
            natom = len(operations[0]['permutation'])
            flat = arr.reshape(-1, natom, 3)
            result = np.zeros_like(flat)
            for op in operations:
                matrix = np.asarray(op[matrix_key], dtype=float)
                permutation = list(op['permutation'])
                # g_{perm(a)} = M g_a  =>  contribution (M^T g)[perm[a]]
                result += np.einsum('kj,sak->saj', matrix, flat[:, permutation, :])
            result /= len(operations)
            meta.setdefault('integral_symmetry', {})['gradient_symmetrized'] = True
            return result.reshape(shape)
        except Exception:
            return grads

    def label_normal_modes(self):
        """Assign abelian irrep labels to normal modes (metadata only, non-fatal).

        Stores the result under ``symmetry_metadata['mode_labels']``.
        """
        meta = self.symmetry_metadata
        if not meta or meta.get('status', 'disabled') == 'disabled':
            return None
        if not meta.get('label_modes', True):
            return None
        detection = meta.get('detection')
        if not detection:
            return None

        try:
            from oqp.library.symmetry import assign_mode_irreps

            modes = np.asarray(self.modes, dtype=float)
            if modes.ndim != 2 or modes.size == 0:
                return None
            natom = len(detection['operations'][0]['permutation'])
            if modes.shape[1] != 3 * natom:
                meta['mode_labels'] = {'status': 'skipped_mode_shape_mismatch'}
                return None

            tolerance = float(meta.get('tolerance', 1.0e-5))
            result = assign_mode_irreps(
                modes,
                detection['operations'],
                detection['character_table'],
                tolerance=max(tolerance, 1.0e-2),
                matrix_key='matrix_input_frame',
            )
            result = dict(result)
            result['status'] = 'ok'
            meta['mode_labels'] = result
            return result
        except Exception as exc:
            # Labeling must never break the run in the metadata-only phase.
            meta['mode_labels'] = {'status': 'error', 'error': str(exc)}
            return None

    @mpi_dump
    def _dump_symmetry_log(self):
        """Append a symmetry summary block to the main log (best effort)."""
        try:
            meta = self.symmetry_metadata
            active = meta.get('integral_symmetry', {})
            full = meta.get('reduction_maps_full')
            lines = [
                '',
                '   ==============================================',
                '   PyOQP: molecular symmetry',
                '   ==============================================',
                f"   detected point group : {meta.get('detected_point_group', '?')}",
                f"   abelian subgroup     : {meta.get('detected_subgroup', '?')}",
            ]
            if full:
                lines.append(f"   full group order     : {full.get('n_operations')}")
            if active:
                status = active.get('status')
                lines.append(f"   integral reduction   : {status}"
                             + (f" (|G| = {active.get('n_operations')}"
                                + (', full group' if active.get('full_group')
                                   else ', abelian subgroup') + ')'
                                if status == 'active' else ''))
                if status != 'active':
                    # The C1 fallback gives a numerically identical answer, so
                    # without this line a user who explicitly asked for the
                    # reduction has no way to discover it never ran.
                    lines.append('   *** the petite-list reduction is NOT active:'
                                 ' this run used the full (C1) integral list ***')
                    if status == 'skipped_basis_mismatch':
                        lines.append(f"   AO count from the symmetry maps ({active.get('n_ao')})"
                                     f" disagrees with the basis ({active.get('nbf')})")
                    elif active.get('error'):
                        lines.append(f"   reason: {active.get('error')}")
                if active.get('reoriented'):
                    lines.append('   geometry reoriented to the symmetry standard orientation')
            response = meta.get('response_symmetry')
            if response:
                lines.append(f"   response blocking    : {response.get('status')}")
            lines.append('')
            with open(self.log, 'a', encoding='utf-8') as fout:
                fout.write('\n'.join(lines))
        except Exception:
            pass

    @mpi_dump
    def _dump_mo_labels_log(self, result):
        """Append MO irrep labels to the main log (best effort, non-fatal)."""
        try:
            meta = self.symmetry_metadata
            lines = [
                '',
                '   ==============================================',
                '   PyOQP: MO symmetry labels (metadata only)',
                '   ==============================================',
                f"   point group:      {meta.get('point_group', '?')}",
                f"   abelian subgroup: {meta.get('subgroup', '?')}",
            ]
            scf_state = result.get('scf_state')
            if scf_state:
                lines.append(f"   SCF state:        {scf_state['term']}")
            for spin in ('alpha', 'beta'):
                if spin not in result:
                    continue
                labels = result[spin]['labels']
                try:
                    energies = np.asarray(
                        self.data[f"OQP::E_MO_{'A' if spin == 'alpha' else 'B'}"],
                        dtype=float,
                    ).ravel()
                except Exception:
                    energies = None
                lines.append(f'   {spin} MOs:')
                for imo, label in enumerate(labels):
                    if energies is not None and imo < energies.size:
                        lines.append(f'   {imo + 1:5d}  {energies[imo]:16.8f}  {label}')
                    else:
                        lines.append(f'   {imo + 1:5d}  {"":16s}  {label}')
            lines.append('')
            with open(self.log, 'a', encoding='utf-8') as fout:
                fout.write('\n'.join(lines))
        except Exception:
            pass

    def get_mass(self):
        """
        Get read-only molar mass
        """
        natom = self.data["natom"]
        atoms = np.frombuffer(
            oqp.ffi.buffer(self.mass,
                           natom * oqp.ffi.sizeof("double"))
        ).astype(float)

        return copy.deepcopy(atoms)

    def get_system(self):
        """
        Get read-only coordinates
        """
        natom = self.data['natom']
        coord = np.frombuffer(
            oqp.ffi.buffer(self.xyz, 3 * natom * oqp.ffi.sizeof("double")),
            dtype=np.double)

        return copy.deepcopy(coord)
    def get_scf_energy(self, component=None):
        """
        Retrieve SCF (Self-Consistent Field) energy components.

        This method provides convenient access to individual or all energy
        terms computed during an SCF procedure. If no component is specified,
        the total SCF energy is returned.

        Parameters
        ----------
        component : str, optional
            The energy component to retrieve. Supported options are:

            - ``None`` (default): Returns only the total SCF energy.
            - ``"all"``: Returns a dictionary containing all available
              energy components.
            - One of the following component names:
                * "energy"  — total SCF energy
                * "psinrm"  — wavefunction norm
                * "ehf1"    — Hartree-Fock energy (one-electron)
                * "vee"     — electron-electron repulsion energy
                * "nenergy" — nuclear energy contribution
                * "vne"     — electron-nucleus attraction energy
                * "vnn"     — nucleus-nucleus repulsion energy
                * "vtot"    — total potential energy
                * "tkin"    — kinetic energy
                * "virial"  — virial ratio

        Returns
        -------
        float or dict
            - If `component` is None, returns a single float (total SCF energy).
            - If `component` is "all", returns a dictionary with all energy components.
            - If `component` corresponds to a specific component, returns that component as a float.

        Raises
        ------
        ValueError
            If the provided `component` does not match any of the known energy components.

        Examples
        --------
        >>> mol.get_scf_energy()
        -75.98327432

        >>> mol.get_scf_energy("tkin")
        37.420192

        >>> mol.get_scf_energy("all")
        {
            'energy': -75.98327432,
            'psinrm': 0.999999,
            'ehf1': -72.3123,
            'vee': 18.2034,
            'nenergy': -80.000,
            'vne': -85.6214,
            'vnn': 5.6214,
            'vtot': -67.4180,
            'tkin': 37.4202,
            'virial': 2.1519
        }
        """
        energy_data = self.data._data.mol_energy

        if component is None:
            return energy_data.energy

        elif component == "all":
            return {
                "energy": energy_data.energy,
                "psinrm": energy_data.psinrm,
                "ehf1": energy_data.ehf1,
                "vee": energy_data.vee,
                "nenergy": energy_data.nenergy,
                "vne": energy_data.vne,
                "vnn": energy_data.vnn,
                "vtot": energy_data.vtot,
                "tkin": energy_data.tkin,
                "virial": energy_data.virial
            }

        else:
            if hasattr(energy_data, component):
                return getattr(energy_data, component)
            else:
                raise ValueError(
                    f"Invalid component '{component}'. Use one of: "
                    f"energy, psinrm, ehf1, vee, nenergy, vne, vnn, "
                    f"vtot, tkin, virial, or 'all'."
                )
    def get_atoms2(self, prop=None):
        """
        Get atomic data.
        Parameters
        ----------
        prop : str or None
            If None, return full atomic data dict.
            If a string, return that specific atomic property.
            Available keys: 'natom', 'coords', 'charge', 'mass'
        Returns
        -------
        dict or np.ndarray
        """
        data = self.data.atomic_data
        if prop is None:
            return data
        if prop not in data:
            raise KeyError(f"Unknown atomic property '{prop}'. "
                           f"Available keys: {list(data.keys())}")
        return data[prop]

    def get_grad(self):
        """
        Get gradient in Hartree/Bohr
        """
        natom = self.data['natom']
        grad = np.frombuffer(
            oqp.ffi.buffer(self.data._data.grad, 3 * natom * oqp.ffi.sizeof("double"))
        )

        return copy.deepcopy(grad)

    def set_grad(self, grad):
        """Write a gradient (natom, 3) or flat 3*natom back into the library buffer.

        ``get_grad`` reads the raw buffer the Fortran layer wrote. Anything the
        Python layer does to a gradient afterwards -- notably
        ``symmetrize_gradient`` -- is invisible to every consumer that goes
        through ``get_grad``: the regression comparison
        (``get_data()['grad']``), and the QM/MM driver. Without this the
        projected gradient reaches the printed output while the raw skeleton
        gradient reaches everything else, so the log looks right and the stored
        result is wrong.
        """
        natom = int(self.data['natom'])
        arr = np.ascontiguousarray(np.asarray(grad, dtype=float).reshape(-1))
        if arr.size != 3 * natom:
            raise ValueError(
                'set_grad expects %d values, got %d' % (3 * natom, arr.size))
        buf = oqp.ffi.buffer(
            self.data._data.grad, 3 * natom * oqp.ffi.sizeof("double"))
        buf[:] = arr.tobytes()

    def get_nac(self):
        """
        Get the non-adiabatic (phase-corrected derivative) coupling matrix d_ij.

        Populated by a NACME run (``self.dcm``); empty for every other runtype.
        The elements are sign/phase ambiguous between builds, so the regression
        comparison uses magnitudes (see the ``nac`` registry entry,
        ``phase_invariant=True``).
        """
        dcm = np.asarray(self.dcm, dtype=float)
        if dcm.size == 0:
            return []
        return dcm.tolist()

    def get_soc(self):
        """
        Get spin-orbit coupling eigenvalues in cm-1
        """
        try:
            return np.array(self.data['OQP::soc_eval']).tolist()
        except AttributeError:
            return []

    def get_hess(self):
        """
        Get hessian results
        """

        return copy.deepcopy(self.hessian)

    def set_hessian_result(self, raw_hessian, asymmetry_tol=1.0e-8):
        """
        Store a final Cartesian Hessian in OpenQP frequency conventions.

        Native analytic Hessian kernels should hand one square ``(3N, 3N)``
        matrix to this helper. The helper records the pre-symmetrization
        asymmetry for diagnostics and stores the symmetrized matrix used by
        normal-mode analysis; it does not compute a numerical fallback.
        """

        hessian = np.asarray(raw_hessian, dtype=float)
        if hessian.ndim != 2 or hessian.shape[0] != hessian.shape[1]:
            raise ValueError(f"Expected square Hessian matrix, got shape={hessian.shape}")

        natom = self.data['natom']
        expected = 3 * natom
        if hessian.shape != (expected, expected):
            raise ValueError(
                f"Expected Hessian shape ({expected}, {expected}) for {natom} atoms, got {hessian.shape}"
            )

        max_asymmetry = float(np.max(np.abs(hessian - hessian.T))) if hessian.size else 0.0
        if max_asymmetry > asymmetry_tol:
            warnings.warn(
                f"Analytic Hessian asymmetry {max_asymmetry:.3e} exceeds tolerance {asymmetry_tol:.3e}; symmetrizing final matrix.",
                RuntimeWarning,
            )

        self.hessian = 0.5 * (hessian + hessian.T)
        self.hessian_metadata = {
            'max_asymmetry': max_asymmetry,
            'symmetrized': bool(max_asymmetry > 0.0),
        }
        return self.hessian

    def _read_mrsf_ekt_records(self):
        """Read the MRSF-EKT root records (eigenvalues, Dyson orbitals, pole
        strengths) from the tagarray, or None when no EKT data is present."""
        if self.data is None:
            return None
        try:
            eigenvalues = np.array(self.data['OQP::mrsf_ekt_eigenvalues'])
            strengths = np.array(self.data['OQP::mrsf_ekt_strengths'])
            orbitals = np.array(self.data['OQP::mrsf_ekt_orbitals_mo'])
        except AttributeError:
            return None

        # The Fortran tagarray reserves nbf slots, but the EKT solver may
        # skip null/near-zero eigenvalue channels and only fill a compact
        # prefix of physically stored Dyson roots. Keep the structured JSON
        # aligned with the printed/stored root ladder instead of exposing
        # trailing all-zero scratch slots.
        ekt_slot_tol = 1.0e-12
        if eigenvalues.ndim == 1 and strengths.ndim == 1 and orbitals.ndim == 2 \
                and strengths.shape[0] == eigenvalues.shape[0]:
            active = (np.abs(eigenvalues) > ekt_slot_tol) \
                | (np.abs(strengths) > ekt_slot_tol)
            if orbitals.shape[0] == eigenvalues.shape[0]:
                active = active | np.any(np.abs(orbitals) > ekt_slot_tol, axis=1)
                eigenvalues = eigenvalues[active]
                strengths = strengths[active]
                orbitals = orbitals[active, :]
            elif orbitals.shape[1] == eigenvalues.shape[0]:
                active = active | np.any(np.abs(orbitals) > ekt_slot_tol, axis=0)
                eigenvalues = eigenvalues[active]
                strengths = strengths[active]
                orbitals = orbitals[:, active]

        hartree_to_ev = 27.211386245988
        return {
            "eigenvalues_hartree": eigenvalues.tolist(),
            "ebe_ev": (-eigenvalues * hartree_to_ev).tolist(),
            "pole_strengths": strengths.tolist(),
            "dyson_orbitals_mo": orbitals.tolist(),
        }

    def snapshot_mrsf_ekt_results(self, kind):
        """Snapshot MRSF-EKT root results ('ip' or 'ea') right after the call.

        The Fortran EKT driver reuses the same OQP::mrsf_ekt_* records for IP
        and EA, so when both are requested in one runtype=ekt job the second
        call overwrites the first.  Snapshotting after each call keeps the
        Dyson orbitals and pole strengths of both kinds for the final JSON.
        """
        records = self._read_mrsf_ekt_records()
        if records is not None:
            self.mrsf_ekt_results_by_kind[kind] = records

    def get_mrsf_ekt_results(self):
        """Collect MRSF-EKT root results for the final JSON file."""
        records = self._read_mrsf_ekt_records()
        if records is None:
            return {}

        ekt_type = self.config.get('tdhf', {}).get('type')
        if self.config.get('input', {}).get('runtype') == 'ekt':
            # The records reflect the most recent EKT call (EA runs after IP).
            if self.config.get('ekt', {}).get('ea'):
                ekt_type = 'mrsf_ekt_ea'
            elif self.config.get('ekt', {}).get('ip'):
                ekt_type = 'mrsf_ekt_ip'

        result = {
            'tdhf_type': ekt_type,
            'target_state': self.config.get('tdhf', {}).get('target'),
            **records,
            # legacy key kept for backward compatibility
            'orbitals_mo': records['dyson_orbitals_mo'],
        }
        # Per-kind snapshots preserve both IP and EA Dyson orbitals and pole
        # strengths when a single runtype=ekt job requests both.
        result.update(self.mrsf_ekt_results_by_kind)
        return {'mrsf_ekt': result}

    def get_data(self):
        """
        Extract data from mol to dict
        """
        scf_type = self.config['scf']['type']
        data = {}
        for key in self.tag:
            if key in self.skip_tag[scf_type]:
                continue
            try:
                data[key] = np.array(self.data[key]).tolist()

            except AttributeError:
                continue

        return data

    def get_data_from_back_door(self):
        """
        Extract mol data for nacme calculation
        """
        if isinstance(self.back_door, tuple):
            return self.back_door
        else:
            # previous data is not available, return current data to bypass nacme calculation
            return self.get_system(), self.get_data()

    def explicit_scf_props(self):
        """Lowercased scf_prop values requested for this calculation.

        ``scf_prop`` defaults to empty, so its config value is exactly the set of
        properties the user asked for -- works identically for file-based and
        scripting-API (input_dict) runs. Regression coverage of a property is
        opt-in: only requested properties are surfaced to the JSON, required by
        the gate, and compared.
        """
        try:
            requested = self.config.get('properties', {}).get('scf_prop', []) or []
        except (AttributeError, TypeError):
            return []
        return [str(p).strip().lower() for p in requested if str(p).strip()]

    def get_results(self):
        """
        Collect computed results to dict
        """
        data = {
            'atoms': self.get_atoms().tolist(),
            'coord': self.get_system().tolist(),
            'energy': self.mol_energy.energy,
            'symmetry_metadata': self.symmetry_metadata,
        }

        # save td energies if available
        try:
            data['td_energies'] = np.array(self.data['OQP::td_energies']
                                           ).tolist()
        except AttributeError:
            data['td_energies'] = np.array([0]).tolist()

        # Surface the spin-resolved excitation ladders as public regression
        # keys when present (a SOC run computes both). The public td_energies
        # mirrors only one ladder, so without these the SOC singlet excitation
        # energies would no longer be compared once the internal OQP:: arrays
        # are trimmed from the references.
        for key, tag in (('td_singlet_energies', 'OQP::td_singlet_energies'),
                         ('td_triplet_energies', 'OQP::td_triplet_energies')):
            try:
                data[key] = np.array(self.data[tag]).tolist()
            except (AttributeError, KeyError, TypeError):
                pass

        # Transition dipoles between excited states, (3, nstates, nstates) in
        # a.u.  Without this the oscillator strengths are printed but never
        # regression-checked, which is how the UMRSF path shipped identically
        # zero transition dipoles: the energies matched their reference while
        # every dipole was zero and nothing compared them.
        try:
            raw = np.asarray(self.data['OQP::td_trans_dipole'],
                             dtype=float).ravel(order='C')
            nstates = int(round((raw.size / 3) ** 0.5))
            if 3 * nstates * nstates != raw.size:
                raise ValueError('unexpected td_trans_dipole tag shape')
            data['td_trans_dipole'] = raw.reshape(
                (3, nstates, nstates), order='F').tolist()
        except (AttributeError, KeyError, TypeError, ValueError):
            pass

        # save NMR isotropic shielding if available (CGO or GIAO).
        # Flat atom-major array -> (natom, 5) in ppm; columns =
        # [dia, para_uncoupled, para_coupled, total_uncoupled, total_coupled].
        try:
            sh = np.array(self.data['OQP::nmr_shielding']).reshape(-1, 5)
            data['nmr_shielding'] = sh.tolist()
        except (AttributeError, KeyError, TypeError, ValueError):
            pass

        # SCF properties are surfaced to the JSON only when EXPLICITLY requested
        # in the input (scf_prop defaults to 'el_mom,mulliken', which every SCF
        # run computes for the log -- we must not bloat/require those on every
        # reference). The explicit request is the regression opt-in.
        props = self.explicit_scf_props()
        charge_tags = (('mulliken', 'mulliken_charges', 'OQP::mulliken_charges'),
                       ('lowdin', 'lowdin_charges', 'OQP::lowdin_charges'),
                       ('resp', 'resp_charges', 'OQP::resp_charges'))
        for prop, key, tag in charge_tags:
            if prop not in props:
                continue
            try:
                data[key] = np.array(self.data[tag]).tolist()
            except (AttributeError, KeyError, TypeError):
                pass
        if 'el_mom' in props:
            try:
                dip = np.zeros(3, dtype=np.float64)
                oqp.electric_dipole_au(
                    self, oqp.ffi.cast("double *", oqp.ffi.from_buffer(dip)))
                data['dipole'] = dip.tolist()
            except Exception:
                pass

        # save gradients if available
        data['grad'] = np.array(self.get_grad()).tolist()
        data['nac'] = np.array(self.get_nac()).tolist()
        data['soc'] = np.array(self.get_soc()).tolist()
        data['hess'] = np.array(self.get_hess()).tolist()
        data.update(self.get_mrsf_ekt_results())

        # Keep the programmatic API and the on-disk JSON on the same public
        # state-tracking contract.  External dynamics drivers (notably
        # PyRAI2MD) must be able to tell that root/phase transport has already
        # been applied without depending on the chosen I/O mode.
        state_tracking = self.get_state_tracking()
        if state_tracking is not None:
            data['state_tracking'] = state_tracking

        # OpenQP-DFTB backend results.  The DFTB adapter is pure Python, so
        # liboqp's mol_energy struct (the 'energy' scalar above) is never
        # populated on this path; surface the adapter results explicitly:
        # total state energies, the excited-state summary (VEE in eV,
        # transition dipoles, oscillator strengths), the reference energy
        # decomposition, and Mulliken charges + point-charge dipole.
        if str(self.config.get('input', {}).get('method', '')
               ).strip().lower() == 'dftb':
            energies = getattr(self, 'energies', None)
            if energies is not None and np.asarray(energies).size:
                data['energies'] = np.asarray(
                    energies, dtype=float).ravel().tolist()
                data['energy'] = data['energies'][0]
            for key in ('dftb_excited_states', 'dftb_energy_components',
                        'dftb_resolved_options', 'dftb_mo',
                        'dftb_mulliken', 'dftb_relaxed_mulliken'):
                value = getattr(self, key, None)
                if value:
                    data[key] = value

        return data

    def get_state_tracking(self):
        """Return the JSON-safe public root/phase transport record, if any."""
        if not getattr(self, '_state_tracking_fresh', False):
            return None
        try:
            return {
                'schema_version': 1,
                'index_base': 0,
                'order_semantics': 'current_root_to_previous_root',
                'order': np.asarray(
                    self.data['OQP::state_tracking_order'], dtype=int
                ).reshape(-1).tolist(),
                'raw_order': np.asarray(
                    self.data['OQP::state_tracking_raw_order'], dtype=int
                ).reshape(-1).tolist(),
                'output_reordered': bool(np.asarray(
                    self.data['OQP::state_tracking_output_reordered'], dtype=int
                ).reshape(-1)[0]),
                'lineage': np.asarray(
                    self.data['OQP::state_tracking_lineage'], dtype=int
                ).reshape(-1).tolist(),
                'phase_step': np.asarray(
                    self.data['OQP::state_tracking_phase_step'], dtype=float
                ).reshape(-1).tolist(),
                'phase_initial': np.asarray(
                    self.data['OQP::state_tracking_phase_initial'], dtype=float
                ).reshape(-1).tolist(),
                'previous_phase_initial': np.asarray(
                    self.data['OQP::state_tracking_previous_phase_initial'], dtype=float
                ).reshape(-1).tolist(),
                'matched_overlap': np.asarray(
                    self.data['OQP::state_tracking_overlap'], dtype=float
                ).reshape(-1).tolist(),
                'margin': np.asarray(
                    self.data['OQP::state_tracking_margin'], dtype=float
                ).reshape(-1).tolist(),
            }
        except (AttributeError, KeyError, TypeError, ValueError, IndexError):
            return None

    @mpi_get_attr
    def get_coord(self, coordinates):
        return coordinates

    def update_system(self, coordinates):
        """
        Modify coordinates in memory
        """
        coordinates = self.get_coord(coordinates)
        coordinates = coordinates.reshape((-1, 3))
        natom = self.data['natom']
        coord = np.frombuffer(
            oqp.ffi.buffer(self.xyz, 3 * natom * oqp.ffi.sizeof("double")),
            dtype=np.double).reshape((natom, 3))
        for at in range(len(coordinates)):
            for c in range(3):
                coord[at, c] = np.float64(coordinates[at, c])

    def update_mol(self, ref_mol):
        """
        Pass data from ref_mol to current mol
        """
        for key in ref_mol.tag:
            try:
                self.data[key] = copy.deepcopy(ref_mol.data[key])

            except AttributeError:
                continue

    def check(self, info):
        """
        Check internal data
        """
        if self.data._data.qn[0] != self.elem[0]:
            raise ValueError(info, self.data._data.qn[0], self.elem[0],
                             'var changed!')
        else:
            print(info, 'var checked!')

    def data_allocate(self):
        """Allocate new oqp data object"""
        if not self.data:
            self.data = OQPData(silent=self.silent)

    def data_deallocate(self):
        """Deallocate oqp data object"""
        self.data = None

    @mpi_get_attr
    def get_config(self, input_source):
        parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA, allow_no_value=True)

        # Determine the type of the input source and process accordingly
        if isinstance(input_source, str):  # Assuming input is a filename
            parser.read(input_source)
        elif isinstance(input_source, dict):  # Assuming input is a dictionary
            parser.load_dict(input_source)
        else:
            raise ValueError("Input must be a filename (str) or a configuration dictionary (dict)")

        # Print configuration if not in silent mode
        if not self.silent:
            parser.print_config()

        # Validate the configuration and apply it
        config = parser.validate()

        return config

    def load_config(self, input_source):
        """
        Load calculation parameters from a file or a dictionary based on the input type.

        input_source: filename (str) or config dictionary (dict)
        """
        self.mpi_manager.set_mpi_comm(self.data)
        self.config = self.get_config(input_source)
        self._resolve_perf(input_source)
        self.data.apply_config(self.config)
        self.data['usempi'] = int(self.usempi)
        self.xyz = self.data._data.xyz
        self.elem = self.data._data.qn
        self.mass = self.data._data.mass
        self.mol_energy = self.data._data.mol_energy
        self.initialize_symmetry_metadata()

        return self

    def _resolve_perf(self, input_source):
        """Apply the `perf` preset to self.config before it is pushed to the control
        struct: fill the performance input keys left at the sentinel 'auto' (an explicit
        value always wins). Env-var free. Stores the resolution report for logging."""
        from oqp.utils import perf_levels
        self.perf_level = self.config.get("input", {}).get("perf", perf_levels.UNSET)
        self.perf_report, self.perf_warns = perf_levels.apply(
            self.config, self.perf_level,
            scf_conv=self.config.get("scf", {}).get("conv"),
            zv_conv=self.config.get("tdhf", {}).get("zvconv"))

    def has_molden_orbitals(self):
        """Return whether the current molecule has an AO basis and SCF MOs."""
        try:
            basis = self.data.get_basis()
            nbf = int(basis['nbf'])
            return (
                MoldenWriter.supports_portable_ordering(basis)
                and np.asarray(self.data['OQP::VEC_MO_A']).size == nbf * nbf
            )
        except (AttributeError, KeyError, TypeError, ValueError):
            return False

    def _viewer_basis_data(self):
        """Portable AO basis metadata used by OpenqpView JSON readers."""
        basis = self.data.get_basis()
        if not MoldenWriter.supports_portable_ordering(basis):
            return {}
        shells = []
        primitive_offset = 0
        for center, angular_momentum, nprimitive in zip(
                basis['centers'], basis['angs'], basis['ncontr']):
            angular_momentum = int(angular_momentum)
            nprimitive = int(nprimitive)
            exponents = np.asarray(
                basis['alpha'][primitive_offset:primitive_offset + nprimitive], dtype=float)
            coefficients = np.asarray(
                basis['coef'][primitive_offset:primitive_offset + nprimitive], dtype=float)
            shells.append({
                'atom_index': int(center),
                'angular_momentum': angular_momentum,
                'shell': MoldenWriter.SHELL_TYPES[angular_momentum],
                'exponents': exponents.tolist(),
                'coefficients': [
                    MoldenWriter.molden_primitive_coefficient(
                        angular_momentum, exponent, coefficient)
                    for exponent, coefficient in zip(exponents, coefficients)
                ],
            })
            primitive_offset += nprimitive

        return {
            'format': 'OpenQP portable basis v1',
            'coefficient_convention': 'molden',
            'basis_function_order': 'molden',
            'spherical_harmonics': MoldenWriter._is_spherical(basis),
            'nbf': int(basis['nbf']),
            'shells': shells,
        }

    def _viewer_orbital_data(self):
        """Portable, renderer-ready SCF orbitals in Molden AO ordering."""
        if not self.has_molden_orbitals():
            return {}

        basis = self.data.get_basis()
        if not MoldenWriter.supports_portable_ordering(basis):
            return {}
        nbf = int(basis['nbf'])
        reorder = MoldenWriter.orbital_reorder(basis)
        scf_type = self.config['scf']['type']

        def spin_data(suffix, noccupied, occupancy):
            orbitals = np.asarray(
                self.data[f'OQP::VEC_MO_{suffix}'], dtype=float).reshape((nbf, nbf))
            energies = np.asarray(self.data[f'OQP::E_MO_{suffix}'], dtype=float).reshape(-1)
            return {
                'energies_hartree': energies.tolist(),
                'occupancies': [occupancy if index < int(noccupied) else 0.0
                                for index in range(nbf)],
                'coefficients': orbitals[:, reorder].tolist(),
            }

        orbitals = {
            'format': 'OpenQP portable orbitals v1',
            'basis_function_order': 'molden',
        }
        if scf_type == 'rhf':
            orbitals['alpha'] = spin_data('A', self.data['nocc'], 2.0)
        else:
            orbitals['alpha'] = spin_data('A', self.data['nelec_A'], 1.0)
            orbitals['beta'] = spin_data('B', self.data['nelec_B'], 1.0)

        return {
            'basis_set': self._viewer_basis_data(),
            'molecular_orbitals': orbitals,
        }

    @staticmethod
    def _dyson_mo_state_rows(records, nbf):
        """Normalize saved EKT coefficients to ``(state, SCF MO)``."""
        orbitals = np.asarray(records.get('dyson_orbitals_mo', []), dtype=float)
        if orbitals.ndim != 2:
            return np.zeros((0, nbf), dtype=float)
        if orbitals.shape[1] == nbf:
            return orbitals
        if orbitals.shape[0] == nbf:
            return orbitals.T
        return np.zeros((0, nbf), dtype=float)

    def _viewer_dyson_data(self):
        """State-specific MRSF-EKT Dyson orbitals in renderer-ready AO form."""
        if not self.has_molden_orbitals():
            return {}

        records_by_kind = dict(self.mrsf_ekt_results_by_kind)
        if not records_by_kind:
            tdhf_type = str(self.config.get('tdhf', {}).get('type', '')).strip().lower()
            runtype = str(self.config.get('input', {}).get('runtype', '')).strip().lower()
            if runtype != 'ekt' and tdhf_type not in ('mrsf_ekt_ip', 'mrsf_ekt_ea'):
                return {}
            records = self._read_mrsf_ekt_records()
            if records is None:
                return {}
            if tdhf_type == 'mrsf_ekt_ea':
                kind = 'ea'
            elif tdhf_type == 'mrsf_ekt_ip':
                kind = 'ip'
            else:
                ekt = self.config.get('ekt', {})
                kind = 'ea' if ekt.get('ea') else 'ip'
            records_by_kind[kind] = records

        basis = self.data.get_basis()
        if not MoldenWriter.supports_portable_ordering(basis):
            return {}
        nbf = int(basis['nbf'])
        reorder = MoldenWriter.orbital_reorder(basis)
        scf_orbitals = np.asarray(
            self.data['OQP::VEC_MO_A'], dtype=float).reshape((nbf, nbf))
        target_state = self.config.get('tdhf', {}).get('target')
        states = []

        for kind in ('ip', 'ea'):
            records = records_by_kind.get(kind)
            if not records:
                continue
            dyson_mo = self._dyson_mo_state_rows(records, nbf)
            eigenvalues = np.asarray(records.get('eigenvalues_hartree', []), dtype=float)
            ebe = np.asarray(records.get('ebe_ev', []), dtype=float)
            strengths = np.asarray(records.get('pole_strengths', []), dtype=float)
            nstate = min(len(dyson_mo), len(eigenvalues), len(strengths))
            dyson_ao = np.matmul(dyson_mo[:nstate], scf_orbitals)[:, reorder]
            for index in range(nstate):
                state_number = index + 1
                state = {
                    'kind': kind.upper(),
                    'state_index': state_number,
                    'label': f'Dyson {kind.upper()} state {state_number}',
                    'parent_state': target_state,
                    'eigenvalue_hartree': float(eigenvalues[index]),
                    'pole_strength': float(strengths[index]),
                    'coefficients': dyson_ao[index].tolist(),
                }
                if index < len(ebe):
                    state['electron_binding_energy_ev'] = float(ebe[index])
                states.append(state)

        if not states:
            return {}
        return {
            'dyson_orbitals': {
                'format': 'OpenQP MRSF-EKT Dyson orbitals v1',
                'basis_function_order': 'molden',
                'states': states,
            }
        }

    def _viewer_frequency_data(self):
        """Portable normal modes shared by normal and Hessian JSON files."""
        frequencies = np.asarray(self.freqs, dtype=float).reshape(-1)
        if frequencies.size == 0:
            return {}
        modes = np.asarray(self.modes, dtype=float).reshape((len(frequencies), -1))
        return {
            'frequency_modes': {
                'format': 'OpenQP normal modes v1',
                'frequencies_cm-1': frequencies.tolist(),
                'normal_mode_eigenvectors': modes.tolist(),
                'normal_mode_eigenvectors_units': (
                    'Cartesian displacement, mass-unweighted, row-major by vibrational mode'
                ),
            }
        }

    @mpi_dump
    def write_molden(self, filename, freqs=None, modes=None, include_dyson=False):
        """Write SCF MOs, frequencies, and explicitly requested Dyson states.

        A basis this writer cannot represent is skipped with a warning rather
        than aborting the calculation -- the orbitals are a by-product, and
        losing them should not lose the energy too. The commonest case is a
        shell beyond g (cc-pV5Z and larger), which used to index past the end
        of SHELL_TYPES and raise IndexError after the SCF had converged.

        The gate is ``supports_portable_ordering``, which is what
        ``has_molden_orbitals`` already uses: besides the angular-momentum
        ceiling it also rejects mixed cartesian/spherical bases and shells with
        no tabulated reordering, either of which would silently permute the MO
        coefficients. Checking here rather than at the call sites covers the
        four unguarded ``write_molden`` calls in single_point.py too.

        No status is returned: this runs under @mpi_dump, which does not call
        the wrapped function at all on non-zero ranks, so any return value
        would be meaningless there.
        """

        basis = self.data.get_basis()
        if not MoldenWriter.supports_portable_ordering(basis):
            bad_l = MoldenWriter.unsupported_angular_momentum(basis)
            if bad_l is not None:
                reason = ('the Molden format defines shells only up to %s '
                          '(l=%d) and this basis contains l=%d'
                          % (MoldenWriter.SHELL_TYPES[-1],
                             MoldenWriter.MAX_ANG, bad_l))
            else:
                reason = ('this basis has no portable Molden ordering (a mixed '
                          'cartesian/spherical basis, or a shell with no '
                          'tabulated reordering)')
            warnings.warn(
                'Skipping Molden output for %s: %s. The calculation itself is '
                'unaffected; set [scf] save_molden=false to silence this.'
                % (filename, reason)
            )
            # A stale file from an earlier run would otherwise survive and look
            # like valid output for this one.  Failing to remove it must not
            # abort the run either: the whole point of skipping is that a
            # by-product should not cost the user the energy they computed.
            if os.path.exists(filename):
                try:
                    os.remove(filename)
                except OSError as err:
                    warnings.warn(
                        'Could not remove the stale Molden file %s (%s). It is '
                        'left in place and does NOT describe this calculation.'
                        % (filename, err)
                    )
            return

        with open(filename, mode='w', encoding='ascii') as fout:
            nat = self.data['natom']
            nbf = basis['nbf']
            mdw = MoldenWriter(fout)
            mdw.write_atoms(nat, self.elem, self.xyz, angstrom=False)
            mdw.write_basis(nat, basis)
            mdw.write_spherical_markers(basis)

            if self.config['scf']['type'] == 'rhf':
                # alpha only
                orbitals = self.data['OQP::VEC_MO_A'].reshape([nbf, nbf])
                eorbitals = self.data['OQP::E_MO_A']
                nocc = self.data['nocc']
                occupancies = (2.0 if i < nocc else 0.0 for i in range(nbf))
                mdw.write_mo(basis, orbitals, eorbitals,
                             occupancies, spin='Alpha')
            else:
                # alpha
                orbitals = self.data['OQP::VEC_MO_A'].reshape([nbf, nbf])
                eorbitals = self.data['OQP::E_MO_A']
                nocc = self.data['nelec_A']
                occupancies = (1.0 if i < nocc else 0.0 for i in range(nbf))
                mdw.write_mo(basis, orbitals, eorbitals,
                             occupancies, spin='Alpha')
                # beta
                orbitals = self.data['OQP::VEC_MO_B'].reshape([nbf, nbf])
                eorbitals = self.data['OQP::E_MO_B']
                nocc = self.data['nelec_B']
                occupancies = (1.0 if i < nocc else 0.0 for i in range(nbf))
                mdw.write_mo(basis, orbitals, eorbitals,
                             occupancies, spin='Beta', header=False)

            if include_dyson:
                dyson = self._viewer_dyson_data().get('dyson_orbitals', {})
                states = dyson.get('states', [])
                if states:
                    mdw.write_mo(
                        basis,
                        np.asarray([state['coefficients'] for state in states]),
                        [state['eigenvalue_hartree'] for state in states],
                        [state['pole_strength'] for state in states],
                        spin='Alpha',
                        header=False,
                        symmetries=[state['label'].replace(' ', '-') for state in states],
                        already_reordered=True,
                    )

            if freqs is not None and modes is not None:
                mdw.write_frequency(self, freqs, modes)

    def set_log(self):
        """
        Set up log file
        """
        if not self.log:
            if platform.uname()[0] == "Windows":
                log_file = b'NUL'
            elif platform.uname()[0] == "Linux":
                log_file = b'/dev/null'
            elif platform.uname()[0] == "Darwin":
                log_file = b'/dev/null'
            else:
                log_file = b'/dev/null'
        else:
            log_file = bytes(str(self.log), encoding='ascii')

        _log_c = oqp.ffi.new("char[]", log_file)
        log_c = oqp.ffi.new("struct Cstring *", [len(log_file), _log_c])

        return log_c

    def set_config_json(self):
        data = {}
        data['json'] = {
            'scf_type': self.config['scf']['type'],
            'basis': self.config['input']['basis'],
            'library': self.config['input']['library'],
            'ispher': self.config['input'].get('ispher', 'auto'),
        }
        # Report the AO dimension and whether pure spherical harmonics are in
        # use (the dimension is reduced vs Cartesian when d/f/g are spherical).
        try:
            basis = self.data.get_basis()
            angs = basis['angs']
            nbf = int(basis['nbf'])
            ncart = int(sum(int((a + 1) * (a + 2) // 2) for a in angs))
            data['json'].update({
                'nbf': nbf,
                'nbf_cartesian': ncart,
                'spherical_harmonics': bool(nbf != ncart),
            })
        except Exception:
            pass
        return data

    @mpi_dump
    def save_data(self, lean=None):
        """
        Save mol data and computed results to json

        Args:
            lean: When True, drop internal ``OQP::`` arrays (density, Fock,
                MO coefficients, overlap/kinetic, SOC/TD/MRSF scratch, etc.)
                before writing. These are never compared by ``check_ref`` and
                are not consumed by any example, so test references can omit
                them to stay ~99% smaller. When None (default) the behaviour
                follows the ``OQP_LEAN_JSON`` environment variable; otherwise
                the full bundle is written so the ``guess=json`` restart
                workflow (``load_data``/``put_data``) keeps its DM/Fock/MO data.
        """
        if self.idx != 1:
            jsonfile = self.log.replace('.log', f'_{self.idx}.json')
        else:
            jsonfile = self.log.replace('.log', '.json')
        data = self.get_data()
        # Keep get_data() in the native tag-array layout for restart, back-door,
        # and NAMD consumers.  Only the on-disk JSON presentation needs the TD
        # response-vector axes repaired.
        data = {key: json_array(key, value) for key, value in data.items()}
        data.update(self.get_results())
        data.update(self.set_config_json())
        data.update(self._viewer_orbital_data())
        data.update(self._viewer_dyson_data())
        data.update(self._viewer_frequency_data())

        if lean is None:
            lean = _env_wants_lean_json()
        if lean:
            data = {k: v for k, v in data.items() if regkeys.lean_keep(k)}

        with open(jsonfile, 'w') as outdata:
            json.dump(data, outdata, indent=2)

    @staticmethod
    def _hessian_cache_value(value):
        """Return a deterministic JSON-safe representation for cache signing."""
        if isinstance(value, dict):
            return {
                str(key): Molecule._hessian_cache_value(value[key])
                for key in sorted(value, key=str)
            }
        if isinstance(value, (list, tuple)):
            return [Molecule._hessian_cache_value(item) for item in value]
        if isinstance(value, set):
            return sorted(Molecule._hessian_cache_value(item) for item in value)
        if isinstance(value, np.ndarray):
            return Molecule._hessian_cache_value(value.tolist())
        if isinstance(value, np.generic):
            return value.item()
        if value is None or isinstance(value, (bool, int, float, str)):
            return value
        return str(value)

    def _hessian_request_signature(self, state):
        """Versioned electronic-model configuration identity for Hessian reuse."""
        config = self.config
        # Include complete model-defining configuration sections so less-visible controls
        # (custom basis libraries/ispher, grids/CAM, PCM radii/epsilon/mode,
        # DFTB parameter sets/shifts, response settings, and QM/MM embedding)
        # cannot silently reuse a Hessian produced by another Hamiltonian.
        model_sections = (
            'input', 'mp2', 'guess', 'pcm', 'dftb', 'symmetry', 'scf',
            'dftgrid', 'tdhf', 'qmmm',
        )
        model_config = {}
        for section in model_sections:
            section_config = dict(config.get(section, {}))
            if section == 'input':
                # Geometry is validated numerically below; runtype and thread
                # count select a workflow/runtime but do not define H(R).
                for key in ('system', 'system2', 'runtype', 'omp_threads'):
                    section_config.pop(key, None)
            model_config[section] = self._hessian_cache_value(section_config)
        return {
            'version': HESSIAN_CACHE_VERSION,
            'state': int(state),
            'model_config': model_config,
        }

    @mpi_dump
    def save_freqs(self, state):
        jsonfile = self.log.replace('.log', '.hess.json')
        data = {
            'atoms': self.get_atoms().tolist(),
            'coord': self.get_system().tolist(),
            'mass': self.get_mass().tolist(),
            'energy': self.energies[state],
            'hessian_cache_version': HESSIAN_CACHE_VERSION,
            'state': int(state),
            'hessian_request': self._hessian_request_signature(state),
            'hessian': self.hessian.tolist(),
            'hessian_metadata': self.hessian_metadata,
            'freqs': self.freqs.tolist(),
            'modes': self.modes.tolist(),
            'frequency_modes': {
                'format': 'OpenQP normal modes v1',
                'frequencies_cm-1': self.freqs.tolist(),
                'normal_mode_eigenvectors': self.modes.tolist(),
                'normal_mode_eigenvectors_units': 'Cartesian displacement, mass-unweighted, row-major by vibrational mode',
            },
            'inertia': self.inertia.tolist(),
            'infrared_intensities': self.infrared_intensities.tolist(),
            'raman_activities': self.raman_activities.tolist(),
            'vibrational_intensity_metadata': self.vibrational_intensity_metadata,
            'infrared_mode_dipole_derivatives': self.infrared_mode_dipole_derivatives.tolist(),
            'raman_mode_polarizability_derivatives': self.raman_mode_polarizability_derivatives.tolist(),
            'symmetry_metadata': self.symmetry_metadata,
        }
        data.update(self._viewer_orbital_data())
        data.update(self._viewer_dyson_data())

        with open(jsonfile, 'w') as outdata:
            json.dump(data, outdata, indent=2)

    def load_data(self):
        # load data from json to mol
        guess_geom = self.config['guess']['continue_geom']
        guess_file = self.config['guess']['file']

        if not os.path.exists(guess_file):
            exit(f'mol object {guess_file} does not exist')

        with open(guess_file, 'r') as indata:
            data = json.load(indata)

        in_atoms = self.get_atoms()
        ld_atoms = np.array(data['atoms'])

        if len(in_atoms) != len(ld_atoms):
            exit('loading data from json, the number of atoms does not match!')

        if np.amax(np.abs(in_atoms - ld_atoms)) > 0:
            exit('loading data from json, the types of atoms does not match!')

        self.put_data(data)
        self.update_config_json()

        if guess_geom:
            self.update_system(np.array(data['coord']))

    def update_config_json(self):
        # Update the configuration from JSON
        config = self.config
        if config['guess']['type'] != 'json':
            return
        if (config['input']['basis'] == config['json']['basis'] and
                config['scf']['init_library'] == config['json']['library']):
            return
        self.config['json']['do_init'] = 'yes'
        self.config['scf']['init_scf'] = self.config['json']['scf_type']
        self.config['scf']['init_basis'] = self.config['json']['basis']
        self.config['scf']['init_library'] = self.config['json']['library']

    def put_data(self, data):
        # convert list to data
        # Keep loaded tracking arrays available as history for a subsequent
        # overlap calculation, but never publish them as this run's result.
        self._state_tracking_fresh = False
        for key in self.tag:
            try:
                self.data[key] = np.array(data[key])

            except KeyError:
                continue
        for key in self.config_tag.keys():
            for item in self.config_tag[key]:
                try:
                    self.config[key][item] = data[key][item]
                except KeyError:
                    print(f"Warning: Key {key} not found in data")
                except Exception as e:
                    print(f"Error: {e}")
        if isinstance(data, dict) and 'symmetry_metadata' in data:
            self.symmetry_metadata = data['symmetry_metadata']


    def read_freqs(self):
        jsonfile = self.log.replace('.log', '.hess.json')

        if not os.path.exists(jsonfile):
            exit(f'hess file {jsonfile} does not exist')

        with open(jsonfile, 'r') as indata:
            data = json.load(indata)

        cached_request = data.get('hessian_request')
        if (data.get('hessian_cache_version') != HESSIAN_CACHE_VERSION
                or not isinstance(cached_request, dict)
                or cached_request.get('version') != HESSIAN_CACHE_VERSION):
            raise ValueError(
                'cached Hessian lacks a current versioned model-configuration/state identity; '
                'recompute it with hess.read=false'
            )

        cached_atoms = np.asarray(data.get('atoms', []), dtype=int).reshape(-1)
        current_atoms = np.asarray(self.get_atoms(), dtype=int).reshape(-1)
        if cached_atoms.shape != current_atoms.shape or not np.array_equal(
                cached_atoms, current_atoms):
            raise ValueError(
                'cached Hessian atoms do not match the current molecule'
            )
        cached_coord = np.asarray(data.get('coord', []), dtype=float).reshape(-1)
        current_coord = np.asarray(self.get_system(), dtype=float).reshape(-1)
        if cached_coord.shape != current_coord.shape or not np.allclose(
                cached_coord, current_coord, rtol=0.0, atol=1.0e-8):
            raise ValueError(
                'cached Hessian geometry does not match the current geometry'
            )
        cached_mass = np.asarray(data.get('mass', []), dtype=float).reshape(-1)
        current_mass = np.asarray(self.get_mass(), dtype=float).reshape(-1)
        if (cached_mass.shape != current_mass.shape
                or not np.all(np.isfinite(cached_mass))
                or not np.allclose(cached_mass, current_mass,
                                   rtol=0.0, atol=1.0e-12)):
            raise ValueError(
                'cached Hessian masses do not match the current isotopic masses'
            )
        requested_state = int(self.config.get('hess', {}).get('state', 0))
        cached_state = data.get('state')
        if cached_state is not None and int(cached_state) != requested_state:
            raise ValueError(
                'cached Hessian state %s does not match requested state %s'
                % (cached_state, requested_state)
            )
        current_request = self._hessian_request_signature(requested_state)
        if cached_request != current_request:
            cached_model = cached_request.get('model_config', {})
            current_model = current_request['model_config']
            mismatched = [
                section for section in current_model
                if cached_model.get(section) != current_model[section]
            ]
            if cached_request.get('state') != requested_state:
                mismatched.append('state')
            raise ValueError(
                'cached Hessian electronic-model configuration/state does not match the '
                'current request (%s)' % ', '.join(sorted(set(mismatched)))
            )

        try:
            energy = float(data['energy'])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError('cached Hessian energy must be a finite scalar') from exc
        if not np.isfinite(energy):
            raise ValueError('cached Hessian energy must be a finite scalar')
        hessian = np.asarray(data['hessian'], dtype=float)
        expected = (current_coord.size, current_coord.size)
        if hessian.shape != expected or not np.all(np.isfinite(hessian)):
            raise ValueError(
                'cached Hessian must be a finite matrix with shape %s' % (expected,)
            )
        if not np.allclose(hessian, hessian.T, rtol=1.0e-8, atol=1.0e-10):
            raise ValueError('cached Hessian matrix must be symmetric')
        self.hessian_metadata = data.get('hessian_metadata', {})
        freqs = np.asarray(data['freqs'], dtype=float)
        modes = np.asarray(data['modes'], dtype=float)
        inertia = np.asarray(data['inertia'], dtype=float)
        if freqs.ndim != 1 or not np.all(np.isfinite(freqs)):
            raise ValueError('cached frequencies must be a finite one-dimensional array')
        if freqs.size == 0 and modes.size == 0:
            modes = np.empty((0, current_coord.size), dtype=float)
        if (modes.shape != (freqs.size, current_coord.size)
                or not np.all(np.isfinite(modes))):
            raise ValueError(
                'cached normal modes must be a finite array with shape (%d, %d)'
                % (freqs.size, current_coord.size)
            )
        if inertia.shape != (3,) or not np.all(np.isfinite(inertia)):
            raise ValueError('cached principal moments of inertia must be three finite values')

        infrared = np.asarray(data.get('infrared_intensities', []), dtype=float)
        raman = np.asarray(data.get('raman_activities', []), dtype=float)
        for name, values in (
                ('infrared_intensities', infrared),
                ('raman_activities', raman)):
            if (values.ndim != 1 or values.size not in (0, freqs.size)
                    or not np.all(np.isfinite(values))):
                raise ValueError(
                    'cached %s must be empty or a finite value for every mode' % name
                )

        infrared_derivatives = np.asarray(
            data.get('infrared_mode_dipole_derivatives', []), dtype=float
        )
        if infrared_derivatives.size == 0:
            infrared_derivatives = np.empty((0, 3), dtype=float)
        if (infrared_derivatives.shape not in ((0, 3), (freqs.size, 3))
                or not np.all(np.isfinite(infrared_derivatives))):
            raise ValueError(
                'cached infrared mode derivatives must have shape (nmode, 3)'
            )

        raman_derivatives = np.asarray(
            data.get('raman_mode_polarizability_derivatives', []), dtype=float
        )
        if raman_derivatives.size == 0:
            raman_derivatives = np.empty((0, 3, 3), dtype=float)
        if (raman_derivatives.shape not in ((0, 3, 3), (freqs.size, 3, 3))
                or not np.all(np.isfinite(raman_derivatives))):
            raise ValueError(
                'cached Raman mode derivatives must have shape (nmode, 3, 3)'
            )
        self.hessian = hessian
        self.freqs = freqs
        self.modes = modes
        self.inertia = inertia
        if (getattr(self, 'energies', None) is not None
                and len(self.energies) > requested_state):
            self.energies[requested_state] = energy
        self.infrared_intensities = infrared
        self.raman_activities = raman
        self.vibrational_intensity_metadata = data.get('vibrational_intensity_metadata', {})
        self.infrared_mode_dipole_derivatives = infrared_derivatives
        self.raman_mode_polarizability_derivatives = raman_derivatives

        return energy, hessian, freqs, modes, inertia

    def regression_context(self):
        """(runtype, excited, props) used to resolve which registry keys apply.

        ``excited`` is True when an excited-state (TDDFT/MRSF/...) method is
        active, which is the gate for comparing ``td_energies`` (a ground-state
        run stores the placeholder [0]). ``props`` is the list of requested
        ``scf_prop`` values (e.g. 'nmr'), the gate for property keys.
        """
        runtype = self.config['input']['runtype']
        tdhf_type = self.config.get('tdhf', {}).get('type')
        excited = bool(tdhf_type) and str(tdhf_type).lower() not in ('none',)
        # Only EXPLICITLY requested properties are regression targets (the
        # default el_mom,mulliken is computed for the log but not tested).
        props = self.explicit_scf_props()
        return runtype, excited, props

    def check_ref(self):
        """Compare runtime results against the reference, driven by the
        regression registry (an allowlist: exactly the declared quantities for
        this runtype/method/properties are compared)."""
        input_stem = os.path.splitext(self.input_file)[0]
        ref_file = input_stem + '.json'
        runtime_data = self.get_data()
        runtime_data.update(self.get_results())
        runtype, excited, props = self.regression_context()

        message = ''
        total_diff = 0

        if not os.path.exists(ref_file):
            message += '   PyOQP reference data is not found (skip and save data)\n'
            # A freshly generated reference keeps only the registry keys.
            self.save_data(lean=True)
            return message, total_diff

        message += f'   PyOQP reference data {ref_file}\n'
        with open(ref_file, 'r') as indata:
            ref_data = json.load(indata)

        compare = regkeys.keys_to_compare(runtype, excited, props)

        # Vibrational quantities (hessian/freqs/IR/Raman) live in the
        # <input>.hess.json sidecar; the matching runtime sidecar was written
        # next to the run log by save_freqs(). Load both lazily.
        ref_side = runtime_side = {}
        if any(e.source == 'sidecar' for e in compare):
            ref_side = _load_json(input_stem + '.hess.json')
            runtime_side = _load_json(self.log.replace('.log', '.hess.json'))

        _MISSING = object()
        for e in compare:
            if e.source == 'sidecar':
                refv = ref_side.get(e.field(), _MISSING)
                rtv = runtime_side.get(e.field(), _MISSING)
            else:
                refv = ref_data.get(e.key, _MISSING)
                rtv = runtime_data.get(e.key, _MISSING)

            # Reference missing the key: a failure only if the key is required.
            if refv is _MISSING or refv is None or refv == []:
                if e.required:
                    total_diff += 1.0
                    message += f'   PyOQP missing reference {e.key:<20} ... failed (1.00000000)\n'
                continue
            if rtv is _MISSING or rtv is None:
                total_diff += 1.0
                message += f'   PyOQP checking {e.key:<20} ... failed (runtime missing)\n'
                continue

            a, b = rtv, refv
            if e.phase_invariant:   # sign/phase ambiguous -> compare magnitudes
                a, b = _abs_nested(a), _abs_nested(b)
            flag, diff = compare_data(a, b, skip_sub=e.skip_sub, rtol=e.rtol)
            total_diff += diff
            message += f'   PyOQP checking {e.key:<20} ... {flag} ({diff:.8f})\n'

        return message, total_diff


def compare_data(data_1, data_2, skip_sub=(), rtol=0.0):
    """
    Compute the numerical differences between two arrays

    ``skip_sub`` names dict sub-keys to ignore during comparison (e.g. the
    phase/sign-ambiguous EKT orbital vectors); the registry supplies these per
    quantity so they are declared in one place.

    ``rtol`` (0 = exact absolute compare) sets a relative tolerance: a value is
    counted as matching when ``|a-b| <= rtol*|b|``, and only the excess beyond
    that tolerance contributes to the reported diff. Used for large-magnitude
    quantities (e.g. SOC, ~1e5 cm^-1) where the absolute round(diff,4) gate
    would otherwise demand far more significant figures than ULP-level noise.
    """
    if isinstance(data_1, dict) or isinstance(data_2, dict):
        if not isinstance(data_1, dict) or not isinstance(data_2, dict):
            return 'failed', 1.0
        diff = 0.0
        for key in sorted(data_2):
            if key in skip_sub:
                continue
            if key not in data_1:
                diff += 1.0
                continue
            _, subdiff = compare_data(data_1[key], data_2[key],
                                      skip_sub=skip_sub, rtol=rtol)
            diff += subdiff
        if np.round(diff, 4) > 0:
            return 'failed', diff
        return 'passed', diff

    if isinstance(data_1, str) or isinstance(data_2, str):
        diff = 0.0 if data_1 == data_2 else 1.0
        if diff > 0:
            return 'failed', diff
        return 'passed', diff

    if data_1 is None or data_2 is None:
        diff = 0.0 if data_1 is data_2 else 1.0
        if diff > 0:
            return 'failed', diff
        return 'passed', diff

    arr_1 = np.array(data_1)
    arr_2 = np.array(data_2)
    if arr_1.shape != arr_2.shape:
        # Some references intentionally store a compact prefix of a longer
        # runtime vector (e.g. EKT roots). Compare the reference-sized prefix
        # when the remaining dimensions agree; otherwise report a clean
        # failure instead of raising a broadcasting ValueError.
        if arr_1.ndim == arr_2.ndim and arr_1.ndim > 0 \
                and arr_1.shape[0] >= arr_2.shape[0] \
                and arr_1.shape[1:] == arr_2.shape[1:]:
            arr_1 = arr_1[:arr_2.shape[0]]
        else:
            return 'failed', 1.0

    if arr_1.size == 0:
        diff = 0.0
    else:
        # Use the maximum element-wise deviation instead of an L1 sum so
        # vector-valued references are judged by per-value numerical drift,
        # not by the number of states/components in the vector.
        absdiff = np.abs(arr_1 - arr_2)
        if rtol > 0:
            # Count only the deviation in excess of the per-element relative
            # tolerance, so noise within rtol*|ref| passes while a genuine
            # regression still shows up in the reported diff.
            absdiff = np.maximum(0.0, absdiff - rtol * np.abs(arr_2))
        diff = float(np.max(absdiff))
    if np.round(diff, 4) > 0:
        return 'failed', diff

    return 'passed', diff


def get_coord(xyz, nat):
    """Get coordinate"""
    return np.frombuffer(oqp.ffi.buffer(xyz, 3 * nat * oqp.ffi.sizeof("double")),
                         dtype=np.double).reshape((nat, 3))


def string_config(config):
    # convert dict value to strings
    str_config = {}
    for section in config.keys():
        str_config[section] = {}
        for option, value in config[section].items():
            if isinstance(value, list) or isinstance(value, tuple):
                value = list2string(value)
            else:
                value = str(value)

            str_config[section][option] = value

    return str_config


def list2string(in_list):
    # convert list to str
    # [1, 2, 3, 4] -> '1, 2, 3, 4'
    # [[1 2], [3, 4]] -> '1 2, 3 4'
    # do not support three-layer list

    if len(in_list) == 0:
        return ''

    str_list = []
    for item in in_list:
        if isinstance(item, list):
            item = ' '.join([str(x) for x in item])
            if isinstance(item[0], list):
                raise ValueError('do not support three-layer list %s' % in_list)
        else:
            item = str(item)

        str_list.append(item)

    str_list = ','.join(str_list)

    return str_list
