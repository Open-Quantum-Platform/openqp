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

        enabled = self._parse_enabled_mode(symmetry.get('enabled', 'false'))

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
                'enabled': symmetry.get('enabled', 'false'),
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
        if sum(_cartesian_shell_size(l) for _, l in pairs) == nbf:
            shells = [(at, l, False) for at, l in pairs]
        elif sum(_spherical_shell_size(l) for _, l in pairs) == nbf:
            # Pure spherical-harmonic basis (ISPHER=1). The component order
            # is assumed to be CCA/libint (m = -l..+l); record the
            # assumption so runs are auditable until the runtime spherical
            # path is validated end-to-end.
            shells = [(at, l, True) for at, l in pairs]
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
        """GAMESS-style reorientation to the standard frame (geometry only).

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

        # Geometry-displacing drivers (optimizers, numerical Hessians, MEP,
        # NEB, ...) must not have the frame rotated under them; the petite
        # reduction is restricted to single-point runtypes for now.
        runtype = str(self.config.get('input', {}).get('runtype', 'energy')).lower()
        if runtype not in ('energy', 'grad', 'prop', 'properties'):
            meta['integral_symmetry'] = {'status': f'skipped_runtype_{runtype}'}
            return False

        try:
            from oqp.library.symmetry_detect import attach_detection_metadata

            # Reorient (design decision recorded in
            # docs/plans/2026-06-07-symmetry-reductions-design.md).
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

            # GAMESS-style contract: ALL outputs (geometry, gradients, MOs)
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
        reason recorded in the metadata.
        """
        meta = self.symmetry_metadata
        if not meta or not meta.get('use_integral_symmetry'):
            return False
        detection = meta.get('detection')
        if not detection:
            return False
        if meta.get('integral_symmetry', {}).get('status') != 'reoriented':
            return False

        try:
            from oqp.library.symmetry import build_reduction_maps

            basis = self.data.get_basis()
            if not basis:
                meta['integral_symmetry'] = {'status': 'skipped_no_basis'}
                return False
            shells = [(int(at), int(l)) for at, l in zip(basis['centers'], basis['angs'])]
            maps = build_reduction_maps(shells, detection['operations'])
            if maps['n_ao'] != int(basis['nbf']):
                meta['integral_symmetry'] = {'status': 'skipped_basis_mismatch'}
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
            try:
                self._dump_symmetry_log()
            except Exception:
                pass
            return True
        except Exception as exc:
            # Fail safe to the C1 path.
            meta['integral_symmetry'] = {'status': 'error', 'error': str(exc)}
            try:
                self.data['OQP::sym_petite_enable'] = np.array([0], dtype=np.int64)
            except Exception:
                pass
            return False

    def stage_response_symmetry(self):
        """Stage per-pair irrep indices for response-space blocking.

        Builds OQP::sym_pair_irrep (1-based irrep index per excitation
        pair, occupied index fastest) from the converged MO labels. Only
        acts when ``use_response_symmetry`` is enabled; bails to the
        unblocked solver on any 'mixed' orbital or inconsistency.
        """
        meta = self.symmetry_metadata
        if not meta or not meta.get('use_response_symmetry'):
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

            self.data['OQP::sym_pair_irrep'] = pair_irrep
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
                if energies is not None and istate < energies.size:
                    lines.append(f'   state {istate + 1:3d}  {energies[istate]:14.8f}  {term}')
                else:
                    lines.append(f'   state {istate + 1:3d}  {"":14s}  {term}')
            lines.append('')
            with open(self.log, 'a', encoding='utf-8') as fout:
                fout.write('\n'.join(lines))
        except Exception:
            pass

    def symmetrize_gradient(self, grads):
        """Project gradients onto the totally symmetric component.

        Valid in the standard orientation when the petite reduction is
        active: g'_a = (1/|G|) sum_op M_op^T g_{perm_op(a)}. Exact for the
        skeleton two-electron gradient and a noise-cleaner for the rest.
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
            arr = np.asarray(grads, dtype=float)
            shape = arr.shape
            natom = len(operations[0]['permutation'])
            flat = arr.reshape(-1, natom, 3)
            result = np.zeros_like(flat)
            for op in operations:
                matrix = np.asarray(op['matrix'], dtype=float)
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
                lines.append(f"   integral reduction   : {active.get('status')}"
                             + (f" (|G| = {active.get('n_operations')}"
                                + (', full group' if active.get('full_group')
                                   else ', abelian subgroup') + ')'
                                if active.get('status') == 'active' else ''))
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

    def get_grad(self):
        """
        Get gradient in Hartree/Bohr
        """
        natom = self.data['natom']
        grad = np.frombuffer(
            oqp.ffi.buffer(self.data._data.grad, 3 * natom * oqp.ffi.sizeof("double"))
        )

        return copy.deepcopy(grad)

    def get_nac(self):
        """
        Get non-adiabatic couping in Hartree/Bohr
        """

        return []

    def get_soc(self):
        """
        Get spin-orbit coupling in cm-1
        """

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

        hartree_to_ev = 27.211386245988
        return {
            'eigenvalues_hartree': eigenvalues.tolist(),
            'ebe_ev': (-eigenvalues * hartree_to_ev).tolist(),
            'pole_strengths': strengths.tolist(),
            'dyson_orbitals_mo': orbitals.tolist(),
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

        # save NMR isotropic shielding if available (CGO or GIAO).
        # Flat atom-major array -> (natom, 5) in ppm; columns =
        # [dia, para_uncoupled, para_coupled, total_uncoupled, total_coupled].
        try:
            sh = np.array(self.data['OQP::nmr_shielding']).reshape(-1, 5)
            data['nmr_shielding'] = sh.tolist()
        except (AttributeError, KeyError, TypeError, ValueError):
            pass

        # save gradients if available
        data['grad'] = np.array(self.get_grad()).tolist()
        data['nac'] = np.array(self.get_nac()).tolist()
        data['soc'] = np.array(self.get_soc()).tolist()
        data['hess'] = np.array(self.get_hess()).tolist()
        data.update(self.get_mrsf_ekt_results())

        return data

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
        self.data.apply_config(self.config)
        self.data['usempi'] = int(self.usempi)
        self.xyz = self.data._data.xyz
        self.elem = self.data._data.qn
        self.mass = self.data._data.mass
        self.mol_energy = self.data._data.mol_energy
        self.initialize_symmetry_metadata()

        return self

    @mpi_dump
    def write_molden(self, filename):
        """Write calculation results in Molden format"""

        with open(filename, mode='w', encoding='ascii') as fout:
            basis = self.data.get_basis()
            nat = self.data['natom']
            nbf = basis['nbf']
            mdw = MoldenWriter(fout)
            mdw.write_atoms(nat, self.elem, self.xyz, angstrom=False)
            mdw.write_basis(nat, basis)

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
            'library': self.config['input']['library']
        }
        return data

    @mpi_dump
    def save_data(self):
        """
        Save mol data and computed results to json
        """
        if self.idx != 1:
            jsonfile = self.log.replace('.log', f'_{self.idx}.json')
        else:
            jsonfile = self.log.replace('.log', '.json')
        data = self.get_data()
        data.update(self.get_results())
        data.update(self.set_config_json())

        with open(jsonfile, 'w') as outdata:
            json.dump(data, outdata, indent=2)

    @mpi_dump
    def save_freqs(self, state):
        jsonfile = self.log.replace('.log', '.hess.json')
        data = {
            'atoms': self.get_atoms().tolist(),
            'coord': self.get_system().tolist(),
            'mass': self.get_mass().tolist(),
            'energy': self.energies[state],
            'hessian': self.hessian.tolist(),
            'hessian_metadata': self.hessian_metadata,
            'freqs': self.freqs.tolist(),
            'modes': self.modes.tolist(),
            'frequency_modes': {
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

        energy = data['energy']
        hessian = data['hessian']
        self.hessian_metadata = data.get('hessian_metadata', {})
        freqs = data['freqs']
        modes = data['modes']
        inertia = data['inertia']
        self.infrared_intensities = np.array(data.get('infrared_intensities', []), dtype=float)
        self.raman_activities = np.array(data.get('raman_activities', []), dtype=float)
        self.vibrational_intensity_metadata = data.get('vibrational_intensity_metadata', {})
        self.infrared_mode_dipole_derivatives = np.array(
            data.get('infrared_mode_dipole_derivatives', []), dtype=float
        )
        self.raman_mode_polarizability_derivatives = np.array(
            data.get('raman_mode_polarizability_derivatives', []), dtype=float
        )

        return energy, hessian, freqs, modes, inertia

    def check_ref(self):
        # compare test data with ref data
        runtype = self.config['input']['runtype']
        ref_file = self.input_file.replace('.inp', '.json')
        runtime_data = self.get_data()
        runtime_data.update(self.get_results())
        skip_keys = [
            'OQP::VEC_MO_A', 'OQP::VEC_MO_B',
            'OQP::td_abxc', 'OQP::td_bvec_mo', 'OQP::td_mrsf_density',
            'OQP::td_states_overlap', 'OQP::state_sign', 'OQP::td_states_phase',
            'OQP::dc_matrix', 'OQP::nac_matrix', 'OQP::DM_A', 'OQP::DM_B', 'OQP::DM_B', 'E_MO_A', 'OQP::Hcore',
            'OQP::SM', 'OQP::TM', 'OQP::FOCK_A', 'OQP::FOCK_B', 'OQP::E_MO_A', 'OQP::E_MO_B', 'OQP::WAO',
            'OQP::mrsf_ekt_density_mo', 'OQP::mrsf_ekt_lagrangian_mo', 'OQP::mrsf_ekt_fock_mo',
            'OQP::mrsf_ekt_orbitals_mo', 'OQP::mrsf_ekt_eigenvalues', 'OQP::mrsf_ekt_strengths',
            'OQP::hf_hessian',
            'json', 'symmetry_metadata'
        ]
        tdhf_type = self.config.get('tdhf', {}).get('type')
        required_ref_keys = []
        if tdhf_type in ('mrsf_ekt_ip', 'mrsf_ekt_ea') or runtype == 'ekt':
            required_ref_keys.append('mrsf_ekt')

        if runtype in ['energy', 'ekt']:
            skip_keys.append('grad')
            skip_keys.append('hess')

        if runtype in ['grad', 'optimize', 'meci', 'mep']:
            skip_keys.append('hess')

        if runtype in ['hess', 'nacme', 'nac']:
            skip_keys.append('grad')

        message = ''
        total_diff = 0

        if os.path.exists(ref_file):
            message += f'   PyOQP reference data {ref_file}\n'

            with open(ref_file, 'r') as indata:
                ref_data = json.load(indata)

            # Hessian runs write detailed Hessian references to a sidecar
            # <input>.hess.json. Keep the primary restart/reference JSON
            # compact, but compare the runtime hess matrix against the
            # sidecar when the primary hess field is intentionally empty.
            if runtype == 'hess' and ref_data.get('hess') == []:
                hess_ref_file = self.input_file.replace('.inp', '.hess.json')
                if os.path.exists(hess_ref_file):
                    with open(hess_ref_file, 'r') as hess_indata:
                        hess_ref_data = json.load(hess_indata)
                    if 'hessian' in hess_ref_data:
                        ref_data['hess'] = hess_ref_data['hessian']

            for key in required_ref_keys:
                if key not in ref_data:
                    total_diff += 1.0
                    message += f'   PyOQP missing reference {key:<20} ... failed (1.00000000)\n'

            for key, value in ref_data.items():
                if key in skip_keys:
                    continue
                if key not in runtime_data:
                    flag, diff = 'failed', 1.0
                else:
                    flag, diff = compare_data(runtime_data[key], value)
                total_diff += diff
                message += f'   PyOQP checking {key:<20} ... {flag} ({diff:.8f})\n'
        else:
            message += '   PyOQP reference data is not found (skip and save data)\n'
            self.save_data()

        return message, total_diff


def compare_data(data_1, data_2):
    """
    Compute the numerical differences between two arrays
    """
    if isinstance(data_1, dict) or isinstance(data_2, dict):
        if not isinstance(data_1, dict) or not isinstance(data_2, dict):
            return 'failed', 1.0
        diff = 0.0
        for key in sorted(data_2):
            if key == 'orbitals_mo':
                # EKT orbital vectors are phase/sign ambiguous between runs;
                # eigenvalues and pole strengths provide the stable regression
                # signal for the structured EKT result.
                continue
            if key not in data_1:
                diff += 1.0
                continue
            _, subdiff = compare_data(data_1[key], data_2[key])
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
        diff = float(np.max(np.abs(arr_1 - arr_2)))
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
