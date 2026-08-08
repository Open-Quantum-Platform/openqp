"""OQP single point class"""
import os
import sys
import oqp
import copy
import time
import shutil
import platform
import subprocess
import multiprocessing
import numpy as np
from numpy import linalg as la
from oqp.molecule import Molecule
from oqp.utils.mpi_utils import MPIManager, MPIPool

try:
    from dftd4.interface import DampingParam, DispersionModel

    dftd_installed = 'dftd4'
except ModuleNotFoundError:
    from oqp.utils.matrix import DampingParam, DispersionModel

    dftd_installed = 'not installed'

from oqp.library.frequency import normal_mode, thermal_analysis
from oqp.utils.file_utils import dump_log, dump_data, write_config, write_xyz

try:
    from oqp.utils.mrsf_reference import (
        build_mrsf_reference_metadata,
        collect_block_diagonal_response,
        collect_state_interaction_response,
        ensemble_occupation_vectors,
        freeze_mrsf_reference_config,
        reference_mo_permutation,
        requires_coupled_response,
    )
except Exception:
    build_mrsf_reference_metadata = None
    collect_block_diagonal_response = None
    collect_state_interaction_response = None
    ensemble_occupation_vectors = None
    freeze_mrsf_reference_config = None
    reference_mo_permutation = None
    requires_coupled_response = None


class Calculator:
    """
    OQP calculator base class

    """

    def __init__(self, mol):
        self.mol = mol
        self.save_mol = mol.config['guess']['save_mol']
        self.export = mol.config['properties']['export']
        self.export_title = mol.config['properties']['title']
        self.exception = mol.config['tests']['exception']
        self.mpi_manager = MPIManager()


class LastStep(Calculator):
    """
    OQP last step single point calculation class

    """

    def __init__(self, mol, param=None):
        super().__init__(mol)

        self.functional = mol.config['input']['functional']
        if len(self.functional) == 0:
            self.functional = 'hf'

        self.do_d4 = mol.config['input']['d4']
        self.res = None
        self.set_param(param)
        self.natom = 0

        dump_log(
            mol,
            title='PyOQP: Dispersion Correction',
            section='dftd',
            info={'type': dftd_installed, 'd4': self.do_d4}
        )

    def get_dispersion(self, mol, grad_list):
        if grad_list:
            do_grad = True
        else:
            do_grad = False

        atoms = mol.get_atoms()
        coordinates = mol.get_system().reshape((-1, 3))
        natom = len(atoms)
        if self.do_d4:
            model = DispersionModel(atoms, coordinates)
            func_for_d4 = 'bhlyp' if self.functional.lower() in ('bhhlyp') else self.functional
            res = model.get_dispersion(DampingParam(method=func_for_d4),
                                       grad=do_grad)

            energy = res['energy']
            if do_grad:
                grad = res['gradient']
            else:
                grad = np.zeros((natom, 3))
        else:
            energy = 0.0
            grad = np.zeros((natom, 3))

        return energy, grad

    def set_param(self, param):
        # TODO pass user-defined parameters to dftd4
        pass

    def compute(self, mol, grad_list=None):
        # do dftd4
        energy, grad = self.get_dispersion(mol, grad_list)
        energies = self.final_energy(energy)

        if not grad_list:
            return energies

        grads = self.final_grad(grad, grad_list)
        return energies, grads

    def final_energy(self, d4_energy):
        # add dftd4 energy
        el = self.mol.energies
        self.mol.energies = [x + d4_energy for x in el]

        dump_log(
            self.mol,
            title='PyOQP: Final Energy',
            section='energy',
            info={'el': el, 'd4': d4_energy}
        )

        # export data
        if self.export:
            dump_data(self.mol, (self.mol.energies, self.export_title), title='ENERGY', fpath=self.mol.log_path)

        # save mol
        if self.save_mol:
            self.mol.save_data()

        return self.mol.energies

    def final_grad(self, d4_grad, grad_list):
        # add dftd4 grad
        el = self.mol.grads
        self.mol.grads = [x + d4_grad for x in el]

        dump_log(
            self.mol,
            title='PyOQP: Final Gradient',
            section='grad',
            info={'el': el, 'd4': d4_grad, 'grad_list': grad_list}
        )

        # export data
        if self.export:
            dump_data(self.mol, (self.mol.grads, self.export_title, grad_list), title='GRADIENT',
                      fpath=self.mol.log_path)

        # save mol
        if self.save_mol:
            self.mol.save_data()

        return self.mol.grads


class SinglePoint(Calculator):
    """
    OQP single point calculation class

    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.method = mol.config['input']['method']
        self.runtype = mol.config['input']['runtype']
        self.functional = mol.config['input']['functional']
        self.basis = mol.config['input']['basis']
        self.library = mol.config['input']['library']
        self.scf_type = mol.config['scf']['type']
        self.scf_maxit = mol.config['scf']['maxit']
        self.forced_attempt = mol.config['scf']['forced_attempt']
        self.alternative_scf = mol.config["scf"]["alternative_scf"]
        self.converger_type = mol.config["scf"]["converger_type"]
        self.stability = mol.config["scf"]["stability"]
        self.scf_mult = mol.config['scf']['multiplicity']
        self.init_scf = mol.config['scf']['init_scf']
        self.init_it = mol.config['scf']['init_it']
        self.init_basis = mol.config['scf']['init_basis']
        self.init_library = mol.config['scf']['init_library']
        self.init_conv = mol.config['scf']['init_conv']
        self.conv = mol.config['scf']['conv']
        self.save_molden = mol.config['scf']['save_molden']
        self.td = mol.config['tdhf']['type']
        self.nstate = mol.config['tdhf']['nstate']
        self.energy_func = {
            'hf': oqp.hf_energy,
            'rpa': oqp.tdhf_energy,
            'tda': oqp.tdhf_energy,
            'sf': oqp.tdhf_sf_energy,
            'mrsf': oqp.tdhf_mrsf_energy,
            'umrsf': oqp.tdhf_umrsf_energy,
            'mrsf_ekt_ip': oqp.tdhf_mrsf_ekt_ip,
            'mrsf_ekt_ea': oqp.tdhf_mrsf_ekt_ea,
            # QMRSF dressings: present only when liboqp is built with the new modules.
            # Degrade gracefully (None) against a pre-QMRSF liboqp; only fails if dispatched.
            'qmrsf_icpt2': getattr(oqp, 'tdhf_qmrsf_icpt2', None),
            'qmrsf_dk': getattr(oqp, 'tdhf_qmrsf_dk', None),
        }

        # initialize state sign
        self.mol.data["OQP::state_sign"] = np.ones(self.nstate)

    def _mrsf_reference_mode(self):
        section = self.mol.config.get('mrsf_ref', {}) if isinstance(self.mol.config, dict) else {}
        mode = str(section.get('mode', 'off')).strip().lower()
        if mode == 'state_average':
            return 'ensemble'
        return mode

    def _update_mrsf_reference_metadata(self):
        if build_mrsf_reference_metadata is None:
            return {}
        mode = self._mrsf_reference_mode()
        if mode == 'off':
            self.mol.mrsf_reference_metadata = {}
            return {}
        previous_scf = {}
        if isinstance(getattr(self.mol, 'mrsf_reference_metadata', None), dict):
            previous_scf = self.mol.mrsf_reference_metadata.get('scf', {})
        config = self.mol.config
        frozen_response_config = False
        if previous_scf and freeze_mrsf_reference_config is not None:
            config = freeze_mrsf_reference_config(
                self.mol.config,
                {'scf': previous_scf},
            )
            frozen_response_config = config is not self.mol.config
        metadata = build_mrsf_reference_metadata(config, self.mol.data)
        if previous_scf:
            metadata['scf'] = previous_scf
            if frozen_response_config:
                metadata['response_reference_source'] = 'scf_applied_open_pairs'
                metadata['warnings'] = list(metadata.get('warnings', [])) + [
                    'MRSF ensemble response references are pinned to the open pairs used by ensemble SCF'
                ]
        self.mol.mrsf_reference_metadata = metadata
        dump_log(self.mol, title='PyOQP: MRSF Ensemble Reference', section='mrsf_ref', info=metadata)
        return metadata

    def _prepare_mrsf_reference_scf(self):
        mode = self._mrsf_reference_mode()
        if mode != 'ensemble':
            return {}
        if build_mrsf_reference_metadata is None or ensemble_occupation_vectors is None:
            raise RuntimeError("mrsf_ref.mode=ensemble is unavailable in this build")

        metadata = build_mrsf_reference_metadata(self.mol.config, self.mol.data)
        alpha_occ, beta_occ = ensemble_occupation_vectors(metadata)
        alpha_occ = np.asarray(alpha_occ, dtype=np.float64)
        beta_occ = np.asarray(beta_occ, dtype=np.float64)

        self.mol.data["OQP::mrsf_ref_occ_a"] = alpha_occ
        self.mol.data["OQP::mrsf_ref_occ_b"] = beta_occ

        # A single-reference ensemble (one open pair) has integer occupations and
        # is just an ordinary ROHF; only a genuinely fractional (>= 2 reference)
        # ensemble must bypass the SCF stability safeguard.
        is_fractional = bool(
            np.any((alpha_occ > 1.0e-6) & (alpha_occ < 1.0 - 1.0e-6))
            or np.any((beta_occ > 1.0e-6) & (beta_occ < 1.0 - 1.0e-6))
        )

        metadata['scf'] = {
            'ensemble_occupations_applied': True,
            'fractional_occupations': is_fractional,
            'occupation_tags': ["OQP::mrsf_ref_occ_a", "OQP::mrsf_ref_occ_b"],
            'occupation_model': 'fixed_mixed_reference_ensemble',
            'applied_open_pairs': metadata.get('open_pairs', []),
            'applied_weights': metadata.get('weights', []),
            'applied_weight_model': metadata.get('weight_model', {}),
            'applied_pair_selection': metadata.get('pair_selection', {}),
            'alpha_sum': float(np.sum(alpha_occ)),
            'beta_sum': float(np.sum(beta_occ)),
            'nmo': int(alpha_occ.size),
        }
        self.mol.mrsf_reference_metadata = metadata
        dump_log(self.mol, title='PyOQP: MRSF Ensemble Reference', section='mrsf_ref', info=metadata)
        return metadata

    def _guard_mrsf_reference_mode(self):
        if self.td != 'mrsf' or requires_coupled_response is None:
            return
        if requires_coupled_response(self.mol.config):
            metadata = self._update_mrsf_reference_metadata()
            raise NotImplementedError(
                "mrsf_ref.mode=ensemble has mixed-reference ensemble SCF support, "
                "but this run path still needs the coupled ensemble-reference MRSF response matrix. "
                f"reference metadata: {metadata}"
            )

    def _run_mrsf_reference_response(self):
        if self.td != 'mrsf' or requires_coupled_response is None:
            return False
        if not requires_coupled_response(self.mol.config):
            return False
        if self.runtype != 'energy':
            raise NotImplementedError(
                "mrsf_ref.mode=ensemble currently supports MRSF energies through "
                "an uncoupled block-diagonal response prototype; gradients, NAC, and "
                "transition-property vector handoff still need the coupled response implementation."
            )
        if reference_mo_permutation is None or collect_block_diagonal_response is None:
            raise RuntimeError("mrsf_ref.mode=ensemble response helpers are unavailable in this build")

        metadata = self._update_mrsf_reference_metadata()
        ensemble = metadata.get('ensemble', {}) if isinstance(metadata, dict) else {}
        references = ensemble.get('references', []) if isinstance(ensemble, dict) else []
        coupling_model = metadata.get('response_coupling_model', 'overlap_offdiagonal')
        if not ensemble.get('available', False) or not references:
            raise RuntimeError(
                "mrsf_ref.mode=ensemble did not produce reference blocks for MRSF response: "
                f"{ensemble.get('reason', ensemble.get('status', 'unknown'))}"
            )

        requested_nstate = int(self.nstate)
        nmo = int(ensemble.get('nmo', 0))
        # Solve each block at the user's nstate.  Do NOT inflate the root count:
        # the Davidson initial guess homes on the physical S0 (the lowest singlet,
        # which for a triplet-ground-state system such as O2 lies *above* the
        # high-spin reference), exactly as ordinary MRSF does.  Requesting extra
        # roots makes the solver converge onto spin-flip ghost states that sit
        # below the reference (unphysical, below the true ground state), which is
        # not S0.  The variational floor below still drops wrong-open-shell
        # references' ghost states.
        block_nstate = requested_nstate
        scf_snapshot = self._snapshot_scf_state()
        mol_energy_snapshot = self._snapshot_mol_energy_state()
        blocks = []
        block_vectors = {}

        try:
            for reference in references:
                permutation = reference_mo_permutation(reference, nmo)
                self._restore_scf_state(scf_snapshot)
                self.mol.data.set_tdhf_nstate(block_nstate)
                self._apply_mrsf_reference_mo_permutation(permutation)
                trial_vector_info = self._apply_mrsf_reference_trial_vector_policy(
                    reference,
                    metadata,
                    permutation,
                )

                block_info = {
                    'reference_id': int(reference.get('id', len(blocks) + 1)),
                    'weight': float(reference.get('weight', 0.0)),
                    'open_pair': [int(item) for item in reference.get('open_pair', [])],
                    'closed_orbitals': [int(item) for item in reference.get('closed_orbitals', [])],
                    'mo_permutation': [int(item) for item in permutation],
                    'trial_vectors': trial_vector_info,
                }
                dump_log(
                    self.mol,
                    title='PyOQP: MRSF Ensemble Reference Block',
                    section='mrsf_ref',
                    info=block_info,
                )

                self.energy_func['mrsf'](self.mol)
                energies = np.asarray(self.mol.data['OQP::td_energies'], dtype=np.float64).ravel().copy()
                block_info['energies'] = [float(item) for item in energies]
                block_info['converged'] = bool(self.mol.mol_energy.Davidson_converged)
                blocks.append(block_info)

                try:
                    block_vectors[block_info['reference_id']] = np.asarray(
                        self.mol.data['OQP::td_bvec_mo'],
                        dtype=np.float64,
                    ).copy()
                except Exception:
                    pass
        finally:
            self._restore_scf_state(scf_snapshot)
            for attr, value in mol_energy_snapshot.items():
                try:
                    setattr(self.mol.mol_energy, attr, value)
                except Exception:
                    pass
            self.mol.data.set_tdhf_nstate(requested_nstate)

        # Anchor the ground state on the physically dominant reference (the
        # frontier ROHF open pair) so alternative references forced into the wrong
        # open-shell orbitals cannot contribute spuriously low states.
        frontier = metadata.get('frontier', {}) if isinstance(metadata, dict) else {}
        anchor_open_pair = frontier.get('current_open_pair') if isinstance(frontier, dict) else None
        if not anchor_open_pair and references:
            anchor_open_pair = references[0].get('open_pair')
        overlap_threshold = float(metadata.get('overlap_threshold', 0.85))

        block_diagonal = collect_block_diagonal_response(
            blocks,
            requested_nstate,
            anchor_open_pair=anchor_open_pair,
        )
        if block_diagonal.get('status') != 'ready':
            raise RuntimeError("mrsf_ref.mode=ensemble produced no finite MRSF response states")

        # The off-diagonal state interaction is computed for diagnostics only.
        # An energy-only overlap off-diagonal coupling is not variationally safe
        # (the non-orthogonal metric drives the lowest root below min(E_i) for
        # negative spin-flip excitations), so it does NOT define the reported
        # energies; a faithful coupling needs the native sigma-action kernel.
        state_interaction = {}
        if collect_state_interaction_response is not None and coupling_model != 'block_diagonal':
            state_interaction = collect_state_interaction_response(
                blocks,
                block_vectors,
                references,
                requested_nstate,
                nmo,
                coupling=coupling_model,
                overlap_threshold=overlap_threshold,
                anchor_open_pair=anchor_open_pair,
            )

        # Floored block-diagonal selection is the reported ground/excited-state
        # energy result.
        final_response = block_diagonal
        selected_energies = np.asarray(final_response['energies'], dtype=np.float64)
        self.mol.data['OQP::td_energies'] = selected_energies
        self._store_mrsf_reference_selected_vectors(block_diagonal, block_vectors)

        target_idx = max(0, min(int(self.mol.config['tdhf']['target']) - 1, selected_energies.size - 1))
        self.mol.mol_energy.excited_energy = float(selected_energies[target_idx])
        self.mol.mol_energy.Davidson_converged = all(bool(block.get('converged')) for block in blocks)

        metadata['response'] = {
            'status': 'implemented_energy_only',
            'model': final_response.get('model', 'block_diagonal_floored'),
            'coupling': 'none_floored_block_diagonal',
            'coupled': False,
            'full_response_kernel': False,
            'energy_only': True,
            'anchor_open_pair': final_response.get('anchor_open_pair', []),
            'variational_floor_hartree': final_response.get('variational_floor_hartree'),
            'dropped_below_floor_count': final_response.get('dropped_below_floor_count', 0),
            'block_state_floor': int(block_nstate),
            'has_offdiagonal_coupling': False,
            'offdiagonal_count': state_interaction.get('offdiagonal_count', 0),
            'max_abs_offdiagonal_hamiltonian': state_interaction.get('max_abs_offdiagonal_hamiltonian', 0.0),
            'redundancy': state_interaction.get('redundancy', {}),
            'kept_count': final_response.get('candidate_count', 0),
            'block_count': len(blocks),
            'requested_states': requested_nstate,
            'selected_states': final_response.get('selected_states', []),
            'candidate_count': final_response.get('candidate_count', 0),
            'raw_candidate_count': block_diagonal.get('raw_candidate_count', 0),
            'skipped_nonconverged_blocks': block_diagonal.get('skipped_nonconverged_blocks', []),
            'block_diagonal': block_diagonal,
            'state_interaction_diagnostic': state_interaction,
            'blocks': blocks,
            'vector_warning': (
                'energy-only floored block-diagonal prototype; '
                'do not use for gradients, NAC, or transition properties'
            ),
        }
        metadata['warnings'] = list(metadata.get('warnings', [])) + [
            'MRSF ensemble energy uses a floored block-diagonal selection (anchored on the '
            'frontier reference); the off-diagonal state interaction is a diagnostic only because '
            'an energy-only coupling is not variationally safe -- the native sigma-action kernel is not included'
        ]
        self.mol.mrsf_reference_metadata = metadata
        dump_log(self.mol, title='PyOQP: MRSF Ensemble Response', section='mrsf_ref', info=metadata)
        return True

    def _apply_mrsf_reference_mo_permutation(self, permutation):
        zero_based = [int(item) - 1 for item in permutation]
        for tag in ('OQP::VEC_MO_A', 'OQP::VEC_MO_B'):
            try:
                vec = np.asarray(self.mol.data[tag])
                # Tagarray stores Fortran MO columns; the NumPy view exposes the
                # orbital index on axis 0, matching the existing swapmo path.
                self.mol.data[tag] = np.asarray(vec[zero_based, :], dtype=vec.dtype).copy()
            except Exception:
                pass
        for tag in ('OQP::E_MO_A', 'OQP::E_MO_B'):
            try:
                values = np.asarray(self.mol.data[tag])
                self.mol.data[tag] = np.asarray(values[zero_based], dtype=values.dtype).copy()
            except Exception:
                pass

    def _apply_mrsf_reference_trial_vector_policy(self, reference, metadata, permutation):
        model = metadata.get('trial_vector_model', {}) if isinstance(metadata, dict) else {}
        mode = str(model.get('mode', 'adaptive')).strip().lower()
        if mode != 'adaptive':
            return {'mode': mode, 'shifted_orbitals': [], 'shifted_positions': []}

        ensemble = metadata.get('ensemble', {}) if isinstance(metadata, dict) else {}
        active = {int(item) for item in ensemble.get('active_open_orbitals', [])}
        closed = {int(item) for item in reference.get('closed_orbitals', [])}
        open_pair = {int(item) for item in reference.get('open_pair', [])}
        active_virtuals = sorted(active - closed - open_pair)
        if not active_virtuals:
            return {'mode': mode, 'shifted_orbitals': [], 'shifted_positions': []}

        position_by_label = {int(label): index + 1 for index, label in enumerate(permutation)}
        shift = float(model.get('active_virtual_shift_hartree', 1.0e6))
        shifted_positions = []
        for tag in ('OQP::E_MO_A', 'OQP::E_MO_B'):
            try:
                values = np.asarray(self.mol.data[tag]).copy()
            except Exception:
                continue
            for orbital in active_virtuals:
                position = position_by_label.get(orbital)
                if position is None:
                    continue
                values[position - 1] += shift
                if tag == 'OQP::E_MO_A':
                    shifted_positions.append(int(position))
            try:
                self.mol.data[tag] = values
            except Exception:
                pass

        return {
            'mode': mode,
            'active_virtual_shift_hartree': shift,
            'shifted_orbitals': [int(item) for item in active_virtuals],
            'shifted_positions': shifted_positions,
        }

    def _store_mrsf_reference_selected_vectors(self, combined, block_vectors):
        selected = combined.get('selected_states', []) if isinstance(combined, dict) else []
        vectors = []
        for item in selected:
            ref_id = int(item.get('reference_id', -1))
            state_index = int(item.get('state_index', 0)) - 1
            bvec = block_vectors.get(ref_id)
            if bvec is None or state_index < 0:
                return
            try:
                vectors.append(np.asarray(bvec)[:, state_index])
            except Exception:
                return
        if not vectors:
            return
        try:
            self.mol.data['OQP::td_bvec_mo'] = np.asarray(vectors, dtype=np.float64).T.copy()
        except Exception:
            pass

    def _prep_guess(self):
        oqp.library.set_basis(self.mol)
        oqp.library.ints_1e(self.mol)
        oqp.library.guess(self.mol)

    def _project_basis(self):
        oqp.library.project_basis(self.mol)

    def _init_convergence(self):
        init_calc = self.energy_func['hf']
        target_basis = self.basis
        target_library = self.library
        if self.init_basis == 'none':
            init_basis = target_basis
            init_library = target_library
        else:
            init_basis = self.init_basis
            init_library = self.init_library

        init_converger = self.mol.config['scf']['init_converger']
        target_converger = self.mol.config['scf']['converger_type']
        self.mol.data.set_scf_converger_type(init_converger)
        self.mol.data.set_scf_conv(self.init_conv)

        if init_basis:
            self.mol.config['input']['basis'] = init_basis
            self.mol.config['input']['library'] = init_library

        self.mol.data.set_scf_maxit(self.init_it)

        if self.init_scf == 'rhf':
            self.mol.config['input']['functional'] = ''
            self.mol.data.set_dft_functional('')
            self.mol.data.set_scf_type('rhf')
            self.mol.data.set_mol_multiplicity(1)

        elif self.init_scf == 'uhf':
            self.mol.config['input']['functional'] = ''
            self.mol.data.set_dft_functional('')
            self.mol.data.set_scf_type('uhf')
            self.mol.data.set_mol_multiplicity(3)

        elif self.init_scf == 'rohf':
            self.mol.config['input']['functional'] = ''
            self.mol.data.set_dft_functional('')
            self.mol.data.set_scf_type('rohf')
            self.mol.data.set_mol_multiplicity(3)

        elif self.init_scf == 'rks':
            self.mol.data.set_scf_type('rhf')
            self.mol.data.set_mol_multiplicity(1)

        elif self.init_scf == 'uks':
            self.mol.data.set_scf_type('uhf')
            self.mol.data.set_mol_multiplicity(3)

        elif self.init_scf == 'roks':
            self.mol.data.set_scf_type('rohf')
            self.mol.data.set_mol_multiplicity(3)
        else:
            raise ValueError(f'Unknown initial scf method {self.init_scf}')

        dump_log(self.mol, title='PyOQP: Initial SCF steps', section='scf')
        self._prep_guess()
        if self.mol.config['guess']['type'] != 'json':
            init_calc(self.mol)
        # save initially converge orbitals
        if self.save_molden:
            guess_file = self.pack_molden_name('init', self.init_scf, self.mol.config['input']['functional'])
            self.mol.write_molden(guess_file)

        if init_basis:
            dump_log(self.mol, title='OQP: Applying Basis Sets for Overlap calculation', section='scf')
            self.mol.data.set_scf_active_basis(1)
            oqp.library.set_basis(self.mol)
            self.mol.data.set_scf_active_basis(0)
            self.mol.config['input']['basis'] = target_basis
            self.mol.config['input']['library'] = target_library
            oqp.library.set_basis(self.mol)

        # set parameters back to normal scf
        self.mol.config['input']['basis'] = target_basis
        self.mol.config['input']['functional'] = self.functional
        self.mol.data.set_scf_converger_type(target_converger)
        self.mol.data.set_dft_functional(self.functional)
        self.mol.data.set_scf_type(self.scf_type)
        self.mol.data.set_scf_maxit(self.scf_maxit)
        self.mol.data.set_scf_conv(self.conv)
        self.mol.data.set_mol_multiplicity(self.scf_mult)
        self._project_basis()
        #        oqp.library.update_guess(self.mol)

        dump_log(self.mol, title='PyOQP: Initial SCF steps done, switching back to normal SCF', section='scf')

    def pack_molden_name(self, cal_type, scf_type, functional):
        """Add information to molden file name"""
        if len(functional) == 0:
            functional = 'hf'
        basis = self.basis.replace('*', 's')
        if self.mol.idx != 1:
            guess_file = self.mol.log.replace('.log', '_%s_%s_%s_%s_%s.molden' % (
                self.mol.idx, cal_type, scf_type, functional, basis))
        else:
            guess_file = self.mol.log.replace('.log', '_%s_%s_%s_%s.molden' % (
                cal_type, scf_type, functional, basis))

        return guess_file

    @staticmethod
    def molden_unique_name(base_filename):
        """Check if the base file exists"""
        print('kk: check base name:', base_filename)
        if not os.path.exists(base_filename):
            # The file does not exist
            final_filename = base_filename
        else:
            # The file exists -- generate a new filename with '_i'
            i = 1
            while True:
                new_filename = base_filename.replace('.molden', f'_{i}.molden')
                # Check if this new filename exists
                if not os.path.exists(new_filename):
                    # If it doesn't exist, use this filename
                    final_filename = new_filename
                    break
                i += 1
        return final_filename

    # IXCORE for XAS (X-ray absorption spectroscopy)
    # Fock matrix is in AO here, so we need to shift it in Fortran after transform it into MO
    def ixcore_shift(self):
        from oqp import ffi
        ixcore = self.mol.config["tdhf"]["ixcore"]
        if ixcore == "-1":  # if default
            return
        ixcore_array = np.array(ixcore.split(','), dtype=np.int32)
        # Shift occupied MO energies before building the TD trial vectors, leaving
        # the requested core orbital(s) available for ixcore excitations.
        noccB = self.mol.data['nelec_B']
        tmp = self.mol.data["OQP::E_MO_A"]
        for i in range(1, noccB + 1):  # 1-based occupied MO indices
            if i not in ixcore_array:
                tmp[i - 1] = -100000  # shift the MO energy down

    def energy(self, do_init_scf=True, restore_scf_converger=True):
        # check method
        if self.method not in ['hf', 'tdhf']:
            raise ValueError(f'Unknown method type {self.method}')

        target_converger = self.mol.config['scf']['converger_type']
        try:
            # compute reference
            ref_energy = self.reference(do_init_scf=do_init_scf)

            # ixcore
            self.ixcore_shift()

            # compute excitations
            if self.method == 'tdhf':
                if self.td in ('qmrsf_icpt2', 'qmrsf_dk'):
                    # Standalone QMRSF pathways: built on the quintet ROHF
                    # reference + their own CAS(4,4) backbone, NOT the triplet
                    # MRSF Davidson, so they bypass excitation()/td_energies.
                    energies = self.qmrsf_standalone(ref_energy)
                else:
                    energies = self.excitation(ref_energy)
            else:
                energies = ref_energy
        finally:
            if restore_scf_converger:
                self.mol.data.set_scf_converger_type(target_converger)

        return energies

    def swapmo(self):
        # swap MO energy and AO coefficient depending on user's request
        swapmo = self.mol.config["guess"]["swapmo"]
        if swapmo:  # if not default (empty)
            swapmo_array = [int(x.strip()) for x in swapmo.split(',')]

            # Initial MO energy and coefficient
            og_val = self.mol.data["OQP::E_MO_A"]
            og_vec = self.mol.data["OQP::VEC_MO_A"]
            # It only takes pairs. If it is not pair, it will be ignored.
            for i, j in zip(swapmo_array[::2], swapmo_array[1::2]):
                og_val[[i - 1, j - 1]] = og_val[[j - 1, i - 1]]
                og_vec[[i - 1, j - 1]] = og_vec[[j - 1, i - 1]]

    def reference(self, do_init_scf=True):
        dump_log(self.mol, title='PyOQP: Entering Electronic Energy Calculation', section='input')

        # Experimental petite-list reduction (no-op unless
        # [symmetry] use_integral_symmetry is enabled): reorient before the
        # guess/basis stage; stage the maps once the basis exists.
        symmetry_on = bool(getattr(self.mol, 'symmetry_metadata', None) and
                           self.mol.symmetry_metadata.get('use_integral_symmetry'))
        if symmetry_on:
            self.mol.reorient_for_integral_symmetry()

        if self.init_scf != 'no' and do_init_scf:
            # do initial scf iteration to help convergence
            self._init_convergence()
        else:
            self._prep_guess()

        self.swapmo()
        self._prepare_mrsf_reference_scf()

        if symmetry_on:
            self.mol.stage_integral_symmetry_maps()

        scf_start = time.perf_counter()
        scf_flag = self._run_scf()
        self.mol.scf_elapsed_s = time.perf_counter() - scf_start

        if not scf_flag:
            dump_log(self.mol, title='PyOQP: SCF energy is not converged', section='end')
            if self.exception is True:
                raise SCFnotConverged()
            else:
                raise RuntimeError("SCF did not converge — stopping current run.")

        if self.save_molden:
            guess_file = self.pack_molden_name('scf', self.scf_type, self.functional)
            self.mol.write_molden(guess_file)

        energy = [self.mol.mol_energy.energy]
        self.mol.energies = energy

        # Metadata-only MO irrep labels (no-op unless symmetry is enabled).
        if getattr(self.mol, 'symmetry_metadata', None):
            self.mol.label_molecular_orbitals()

        self._update_mrsf_reference_metadata()

        return energy

    def _scf_features(self):
        """Cheap per-system features for the SCF manager."""
        cfg = self.mol.config
        func = (cfg.get('input', {}).get('functional', '') or '').strip()
        scft = str(cfg.get('scf', {}).get('type', 'rhf')).lower()
        try:
            mult = int(cfg.get('scf', {}).get('multiplicity', 1))
        except (TypeError, ValueError):
            mult = 1
        try:
            z = self.mol.get_atoms()
            tm = bool(((z >= 21) & (z <= 30)).any() or ((z >= 39) & (z <= 48)).any()
                      or ((z >= 57) & (z <= 80)).any())
            natom = int(len(z))
        except Exception:
            tm, natom = False, 0
        try:
            charge = int(cfg.get('input', {}).get('charge', 0))
        except (TypeError, ValueError):
            charge = 0
        return dict(is_dft=len(func) > 0, open_shell=(scft in ('uhf', 'rohf') or mult > 1),
                    is_rohf=(scft == 'rohf'), is_uhf=(scft == 'uhf'),
                    transition_metal=tm, natom=natom, mult=mult, charge=charge)

    def _select_converger(self, mode, stability):
        """SCF manager (converger_type=auto|ml): choose the primary converger from
        cheap system features, calibrated on the SCF-convergence database.

        Data summary (OpenQP, 58 cells): C-DIIS is the best primary (wins 52/58, mean
        ~21 Fock builds vs A-DIIS 39, E-DIIS 53); the hard class is open-shell
        transition-metal systems (esp. DFT), where C-DIIS can stall/find a wrong basin
        -> keep C-DIIS but enable the TRAH stability safeguard. The reactive escalation
        ladder (DIIS->SOSCF->TRAH) remains the safety net for all classes.

        mode='ml' uses a trained, distilled model if shipped in pyoqp; otherwise it
        transparently falls back to these rules.
        """
        f = self._scf_features()
        used = 'rules'
        primary = 'diis'   # C-DIIS: best default

        if mode == 'ml':
            try:
                from oqp.library.scf_selector_model import predict as _ml_predict  # distilled, optional
                primary = _ml_predict(f)
                used = 'ML-model'
            except Exception:
                used = 'ML(no model -> rules)'

        # ROHF/UHF stability is the priority: ROKS is the MRSF-TDDFT reference, so an
        # open-shell SCF that lands on a non-lowest (unstable) solution corrupts the
        # excited-state calculation. Enable the stability safeguard for open-shell
        # references (especially DFT/ROKS and transition-metal), where it matters most.
        hard = f['open_shell'] and (f['transition_metal'] or f['is_dft'])
        if used != 'ML-model' and hard:
            stability = True

        dump_log(self.mol, section='',
                 title='PyOQP SCF manager [%s]: %s%s%s%s -> primary=%s%s' % (
                     used,
                     'open-shell ' if f['open_shell'] else 'closed-shell ',
                     'TM ' if f['transition_metal'] else '',
                     'DFT' if f['is_dft'] else 'HF',
                     ' natom=%d' % f['natom'],
                     primary, ' +stability' if (stability and hard) else ''))
        return primary, stability

    def _run_scf(self):
        """Unified robust SCF driver.

        Replaces the old ``forced_attempt`` / ``alternative_scf`` retry loop
        with a single coherent robustness ladder:

          1. **Primary converger** (``scf.converger_type``, default DIIS) —
             fast, gets most cases. ``auto``/``ml`` let the SCF manager pick it.
          2. **Escalation ladder** — if the primary converger does not converge,
             walk a chain of progressively more robust (and costlier) methods,
             each warm-started from the previous orbitals. The default chain is
             ``SOSCF`` (cheap second-order, fixes most DIIS stalls) → then
             ``scf.alternative_scf`` (default TRAH, a globally convergent
             trust-region method). Set ``scf.escalation`` to a comma-separated
             list (e.g. ``soscf,trah``) to override the chain explicitly;
             ``scf.alternative_scf`` (back-compat) sets only the final method.
             The primary converger is dropped from the chain and duplicates are
             removed while preserving order.
          3. **Stability safeguard** (``scf.stability``, default on) — seed a
             stability-following TRAH pass from the converged orbitals.  At a
             genuine minimum this is a ~0-iteration no-op; when the converged
             point is an unstable saddle it relaxes to the lowest solution.
             This catches the case where DIIS *converges* to a non-aufbau /
             non-lowest open-shell (UHF/ROHF) solution and would otherwise be
             returned silently.

        Returns
        -------
        bool
            True if a converged SCF solution was obtained.
        """
        data = self.mol.data
        scf_config = self.mol.config.get('scf', {})
        primary = getattr(self, 'converger_type', scf_config.get('converger_type', 'diis'))
        fallback = getattr(self, 'alternative_scf', scf_config.get('alternative_scf', 'trah'))
        stability = getattr(self, 'stability', scf_config.get('stability', False))
        trah_stab_default = scf_config.get('trh_stab', False)

        # --- SCF manager: converger_type=auto|ml picks the primary from system features ---
        if str(primary).lower() in ('auto', 'ml'):
            primary, stability = self._select_converger(str(primary).lower(), stability)

        # --- Stage 1: primary converger ---
        data.set_scf_converger_type(primary)
        self.scf()
        converged = self.mol.mol_energy.SCF_converged
        if converged:
            dump_log(self.mol, title='PyOQP: SCF converged with %s' % primary, section='')

        # --- Stage 2: escalate the converger on non-convergence ---
        # Escalation ladder of progressively more robust (and costlier) methods.
        # Default: SOSCF (cheap second-order, fixes most DIIS stalls) BEFORE TRAH
        # (globally convergent trust region). 'scf.escalation' overrides the
        # comma-separated chain; 'scf.alternative_scf' (back-compat) sets the
        # final method. Each stage warm-starts from the previous orbitals.
        escalation = scf_config.get('escalation', None)
        if escalation:
            chain = [c.strip().lower() for c in str(escalation).split(',') if c.strip()]
        else:
            chain = ['soscf']
            if fallback:
                chain.append(fallback)
        # Drop the primary and de-duplicate while preserving order.
        _seen = set()
        chain = [c for c in chain if c != primary and not (c in _seen or _seen.add(c))]
        for conv in chain:
            if converged:
                break
            dump_log(self.mol,
                     title='PyOQP: SCF not converged; escalating to %s' % conv,
                     section='input')
            data.set_scf_converger_type(conv)
            data.set_sd_scf(False)
            self.scf()
            converged = self.mol.mol_energy.SCF_converged

        # --- Stage 3: stability safeguard ---
        # Applied to ground-state targets (method='hf') AND to excited-state
        # reference SCFs (method='tdhf', e.g. MRSF): a DIIS-converged but
        # *unstable* open-shell solution is just as wrong a reference for MRSF as
        # it is a wrong ground state, and using it makes the MRSF reference (and
        # hence the excited states) disagree with the standalone ROHF along a
        # PES.  The safeguard only KEEPS the relaxed orbitals when TRAH finds a
        # genuinely lower solution (e_post < e_pre); an energy-invariant
        # re-canonicalization (no lowering) is reverted below by restoring the
        # snapshot, so a sensitive excited-state gradient at a genuine minimum is
        # not perturbed.
        #
        # EXCEPTION: a mixed-reference *ensemble* SCF uses fixed fractional
        # occupations that are deliberately higher in energy than any single
        # determinant.  The stability safeguard would "relax" that symmetric
        # ensemble mean field straight back to the lowest integer ROHF, defeating
        # the ensemble.  Skip it whenever the ensemble occupation tags are staged.
        ens_meta = getattr(self.mol, 'mrsf_reference_metadata', None)
        ensemble_scf = bool(
            isinstance(ens_meta, dict)
            and ens_meta.get('scf', {}).get('ensemble_occupations_applied')
            and ens_meta.get('scf', {}).get('fractional_occupations')
        )
        # No orbital-rotation space (no virtual orbitals) -> the SCF determinant is
        # trivially stable, and the TRAH stability Hessian over an empty occupied-virtual
        # block would issue an empty cblas_dgemm (LDA=0) that hard-aborts under Accelerate.
        # This occurs for fully-occupied high-spin references (e.g. a minimal-basis quintet
        # where nbf == n_occupied). Skip the safeguard in that degenerate case.
        # Two degenerate cases produce an empty occupied-virtual rotation block:
        #   (a) nbf == n_occupied      -> no virtual orbitals at all (minimal-basis quintet);
        #   (b) min(nelec_A, nelec_B) == 0 -> one spin channel has no occupied orbitals
        #       (a maximally high-spin reference, e.g. an M_s=+2 quintet with nelec_B=0),
        #       so that spin's occ-vir TRAH Hessian block is empty -> dgemm LDB=0 abort.
        # Either way there is nothing for the stability safeguard to relax; skip it.
        no_rot = False
        try:
            nbf = int(self.mol.data.get_basis()['nbf'])
            nelec_a = int(self.mol.data['nelec_A'])
            nelec_b = int(self.mol.data['nelec_B'])
            nocc_max = max(nelec_a, nelec_b)
            no_rot = (nbf <= nocc_max) or (min(nelec_a, nelec_b) == 0)
        except Exception:
            no_rot = False
        if no_rot and converged and stability:
            dump_log(self.mol, section='',
                     title='PyOQP: SCF stability safeguard skipped (empty occ-virtual '
                           'rotation block: no virtuals or a fully spin-polarized reference)')
        if (converged and stability and primary != 'trah' and not no_rot
                and self.method in ('hf', 'tdhf') and not ensemble_scf):
            e_pre = self.mol.mol_energy.energy
            mol_energy_snapshot = self._snapshot_mol_energy_state()
            # Snapshot the converged orbitals so the safeguard is a true no-op
            # at a stable minimum (TRAH may re-canonicalize/rotate orbitals
            # energy-invariantly, which would otherwise perturb sensitive
            # downstream quantities such as range-separated excited gradients).
            snapshot = self._snapshot_scf_state()
            dump_log(self.mol, title='PyOQP: Verifying SCF stability (TRAH)', section='input')

            # Stability-following explores symmetry-breaking rotations whose
            # densities are not totally symmetric, so the petite-list
            # reduction must be off during (and after, if a broken-symmetry
            # solution is kept) this stage.
            petite_staged = self._petite_is_staged()
            if petite_staged:
                self._set_petite_enabled(False)

            data.set_scf_converger_type('trah')
            data.set_trah_stability(True)
            data.set_sd_scf(False)
            self.scf()
            trah_ok = self.mol.mol_energy.SCF_converged
            e_post = self.mol.mol_energy.energy

            if trah_ok and e_post < e_pre - 1.0e-7:
                # The converged point was unstable: keep the lower solution.
                # The kept density may be symmetry-broken: petite stays off.
                if petite_staged:
                    self.mol.symmetry_metadata['integral_symmetry']['status'] = \
                        'disabled_symmetry_broken_scf'
                dump_log(self.mol,
                         title='PyOQP: SCF point was unstable; relaxed to a lower '
                               'solution (dE = %.3e Hartree)' % (e_post - e_pre),
                         section='')
            else:
                # Stable (no lower solution found) or the verification did not
                # converge: restore the original converged orbitals unchanged.
                self._restore_scf_state(snapshot)
                for attr, value in mol_energy_snapshot.items():
                    try:
                        setattr(self.mol.mol_energy, attr, value)
                    except Exception:
                        pass
                # Symmetric solution kept: the petite reduction is valid again.
                if petite_staged:
                    self._set_petite_enabled(True)
                if not trah_ok:
                    # Re-run the primary converger (warm-started) so mol_energy
                    # is consistent with the restored orbitals.
                    data.set_scf_converger_type(primary)
                    self.scf()
                    converged = self.mol.mol_energy.SCF_converged

            # restore the user-configured stability flag for later SCF calls
            data.set_trah_stability(trah_stab_default)

        # restore the primary converger for any subsequent reference() calls
        data.set_scf_converger_type(primary)
        return converged

    # Wavefunction tags that define an SCF solution (alpha + beta channels).
    _scf_state_tags = (
        'OQP::VEC_MO_A', 'OQP::E_MO_A', 'OQP::DM_A', 'OQP::FOCK_A',
        'OQP::VEC_MO_B', 'OQP::E_MO_B', 'OQP::DM_B', 'OQP::FOCK_B',
    )

    def _petite_is_staged(self):
        """True when the petite-list reduction maps are staged and active."""
        meta = getattr(self.mol, 'symmetry_metadata', None)
        if not meta:
            return False
        return meta.get('integral_symmetry', {}).get('status') == 'active'

    def _set_petite_enabled(self, enabled):
        import numpy as np
        try:
            self.mol.data['OQP::sym_petite_enable'] = \
                np.array([1 if enabled else 0], dtype=np.int64)
        except Exception:
            pass

    def _snapshot_scf_state(self):
        """Deep-copy the tags that define the current converged SCF solution."""
        snap = {}
        for tag in self._scf_state_tags:
            try:
                snap[tag] = np.array(self.mol.data[tag]).copy()
            except Exception:
                pass
        return snap

    def _restore_scf_state(self, snap):
        """Write a previously captured SCF solution back into the molecule."""
        for tag, value in snap.items():
            try:
                self.mol.data[tag] = value
            except Exception:
                pass

    _mol_energy_state_attrs = (
        'energy',
        'SCF_converged',
        'Davidson_converged',
    )

    def _snapshot_mol_energy_state(self):
        """Capture scalar energy/convergence metadata that must match SCF tags."""
        mol_energy = self.mol.mol_energy
        snap = {}
        try:
            snap.update(getattr(mol_energy, '__dict__', {}))
        except Exception:
            pass
        for attr in self._mol_energy_state_attrs:
            if attr not in snap and hasattr(mol_energy, attr):
                try:
                    snap[attr] = getattr(mol_energy, attr)
                except Exception:
                    pass
        return snap

    def qmrsf_standalone(self, ref_energy):
        """Run a standalone QMRSF pathway on the converged quintet ROHF reference.

        The Fortran routine (tdhf_qmrsf_icpt2 / tdhf_qmrsf_dk) builds the active
        integrals via the int2-reuse AO->MO transform, solves the CAS(4,4)
        backbone, applies the dynamic-correlation dressing, and writes its
        results to the log + a validation dump (qmrsf_*_live.dat). SCF has
        already run in reference(); here we just dispatch and report.
        """
        dump_log(self.mol, title='PyOQP: QMRSF pathway (%s)' % self.td, section='tdhf')
        response_start = time.perf_counter()
        fn = self.energy_func.get(self.td)
        if fn is None:
            raise RuntimeError(
                'liboqp lacks %s — rebuild liboqp with the QMRSF modules '
                '(source/modules/qmrsf_*.F90, tdhf_qmrsf_*.F90).' % self.td)
        if self.td == 'qmrsf_dk':
            gauge = self.mol.canonicalize_qmrsf_active_orbitals()
            dump_log(self.mol, title='PyOQP: QMRSF active-orbital gauge',
                     section='', info=gauge)
            # The native routine returns without writing its results file when
            # a manifold fails to diagonalize or the reference is rejected.  A
            # dump left by an earlier calculation in this directory would then
            # be read against the new reference energy, so remove it first and
            # let its absence report the failure.
            stale_dump = os.path.join(os.getcwd(), 'qmrsf_dk_full_live.dat')
            if os.path.isfile(stale_dump):
                os.remove(stale_dump)
        fn(self.mol)
        response_elapsed = time.perf_counter() - response_start

        qmrsf_energies = None
        # Post-process the icPT2 validation dump into a clean JSON + log table.
        # The Fortran routine writes 'qmrsf_icpt2_full_live.dat' into the run's
        # cwd; failure here must never abort the run (the scalar reference energy
        # is already in hand), so the whole block is guarded.
        if self.td == 'qmrsf_icpt2':
            try:
                from oqp.library.qmrsf_results import (
                    parse_qmrsf_icpt2_dump,
                    build_qmrsf_results,
                    write_qmrsf_json,
                    format_qmrsf_log_table,
                )
                dump_path = os.path.join(os.getcwd(), 'qmrsf_icpt2_full_live.dat')
                if os.path.isfile(dump_path):
                    dump = parse_qmrsf_icpt2_dump(dump_path)
                    ref_scalar = ref_energy[0] if isinstance(ref_energy, (list, tuple)) else ref_energy
                    results = build_qmrsf_results(dump, ref_scalar)
                    # JSON next to the log: replace the log extension with
                    # '.qmrsf.json' (or append it when the log has no extension).
                    log_path = self.mol.log
                    base, ext = os.path.splitext(log_path)
                    json_path = (base if ext else log_path) + '.qmrsf.json'
                    write_qmrsf_json(results, json_path)
                    self.mol.qmrsf_results = results
                    # Reported spectrum = Dyall totals (intruder-robust production
                    # denominator); EN is kept in the JSON/results as a diagnostic.
                    qmrsf_energies = [s['E_icPT2_Dyall'] for s in results['states']]
                    dump_log(self.mol, title=format_qmrsf_log_table(results), section='')
                    dump_log(self.mol,
                             title='PyOQP: QMRSF-icPT2 results written to %s' % json_path,
                             section='')
                else:
                    dump_log(self.mol,
                             title='PyOQP: QMRSF-icPT2 dump not found (%s); '
                                   'skipping results output' % dump_path,
                             section='')
            except Exception as err:
                dump_log(self.mol,
                         title='PyOQP: QMRSF-icPT2 results post-processing failed '
                               '(%s); continuing' % err,
                         section='')

        # Post-process the DK validation dump into a clean JSON + log table.
        # The Fortran routine writes 'qmrsf_dk_full_live.dat' into the run's cwd;
        # failure here must never abort the run, so the whole block is guarded.
        if self.td == 'qmrsf_dk':
            try:
                from oqp.library.qmrsf_results import (
                    parse_qmrsf_dk_dump,
                    build_qmrsf_dk_results,
                    write_qmrsf_json,
                    format_qmrsf_dk_log_table,
                )
                dump_path = os.path.join(os.getcwd(), 'qmrsf_dk_full_live.dat')
                # A dump left by an earlier calculation in the same directory
                # would be combined with this run's reference energy, so the
                # file is removed before the native call above and its absence
                # here means the response did not complete.
                if os.path.isfile(dump_path):
                    dump = parse_qmrsf_dk_dump(dump_path)
                    ref_scalar = ref_energy[0] if isinstance(ref_energy, (list, tuple)) else ref_energy
                    results = build_qmrsf_dk_results(dump, ref_scalar)
                    results['timing'] = {
                        'reference_scf_s': float(getattr(self.mol, 'scf_elapsed_s', 0.0)),
                        'qmrsf_response_s': float(response_elapsed),
                    }
                    log_path = self.mol.log
                    base, ext = os.path.splitext(log_path)
                    json_path = (base if ext else log_path) + '.qmrsf_dk.json'
                    write_qmrsf_json(results, json_path)
                    self.mol.qmrsf_results = results
                    # On a KS/ROKS reference report the DFT-dressed DK spectrum
                    # (the genuine Pathway-II value); otherwise the bare DK.
                    if results.get('is_dft_dressed') and \
                            all('E_DK_DFT' in s for s in results['states']):
                        qmrsf_energies = [s['E_DK_DFT'] for s in results['states']]
                    else:
                        qmrsf_energies = [s['E_DK'] for s in results['states']]
                    dump_log(self.mol, title=format_qmrsf_dk_log_table(results), section='')
                    dump_log(self.mol,
                             title='PyOQP: QMRSF-DK results written to %s' % json_path,
                             section='')
                else:
                    dump_log(self.mol,
                             title='PyOQP: QMRSF-DK dump not found (%s); '
                                   'skipping results output' % dump_path,
                             section='')
            except Exception as err:
                dump_log(self.mol,
                         title='PyOQP: QMRSF-DK results post-processing failed '
                               '(%s); continuing' % err,
                         section='')

        # results() returns the QMRSF spectrum when available -- the icPT2 Dyall totals
        # (intruder-robust default) or the DK dressed totals -- so the standard OpenQP
        # energy interface exposes the dressed states; otherwise it falls back to the
        # reference scalar. Full per-state detail is on self.mol.qmrsf_results + the JSON.
        self.mol.energies = qmrsf_energies if qmrsf_energies is not None else ref_energy
        return self.mol.energies

    def excitation(self, ref_energy):
        # Response-space symmetry blocking (no-op unless
        # [symmetry] use_response_symmetry is enabled).
        if getattr(self.mol, 'symmetry_metadata', None) and \
                self.mol.symmetry_metadata.get('use_response_symmetry'):
            self.mol.stage_response_symmetry()

        self.tddft()
        energies = ref_energy + [ex + ref_energy[0] for ex in self.mol.data['OQP::td_energies']]

        # check convergence
        td_flag = self.mol.mol_energy.Davidson_converged

        if self.method != 'hf' and not td_flag:
            dump_log(self.mol, title='PyOQP: TD energy is not converged', section='end')

            if self.exception is True:
                raise TDnotConverged()
            else:
                exit()

        self.mol.energies = energies

        # Metadata-only state irrep labels (no-op unless symmetry enabled).
        if getattr(self.mol, 'symmetry_metadata', None):
            self.mol.label_excited_states()

        return energies

    def scf(self):
        # do SCF
        dump_log(self.mol, title='PyOQP: Normal SCF steps', section='scf')
        self.energy_func['hf'](self.mol)

    def tddft(self):
        if self.td == 'mrsf' and requires_coupled_response is not None and \
                requires_coupled_response(self.mol.config) and self.runtype != 'energy':
            raise NotImplementedError(
                "mrsf_ref.mode=ensemble is currently scoped to energy-only MRSF. "
                "Do not use it with EKT, gradients, NAC, RT-MRSF, or transition-property workflows yet."
            )

        if self.runtype == 'ekt':
            if self.td != 'mrsf':
                raise ValueError('EKT runtype only supports MRSF-TDDFT: set [tdhf] type=mrsf')
            ekt_ip = self.mol.config['ekt']['ip']
            ekt_ea = self.mol.config['ekt']['ea']
            if not ekt_ip and not ekt_ea:
                raise ValueError('EKT runtype requires [ekt] ip=True and/or ea=True')
            dump_log(self.mol, title='PyOQP: MRSF-EKT steps', section='tdhf')
            if ekt_ip:
                self.energy_func['mrsf_ekt_ip'](self.mol)
                self.mol.snapshot_mrsf_ekt_results('ip')
            if ekt_ea:
                self.energy_func['mrsf_ekt_ea'](self.mol)
                self.mol.snapshot_mrsf_ekt_results('ea')
            return

        # check td type
        if self.td not in ['rpa', 'tda', 'sf', 'mrsf', 'umrsf', 'mrsf_ekt_ip', 'mrsf_ekt_ea']:
            raise ValueError(f'Unknown tdhf type {self.td}')

        if self._run_mrsf_reference_response():
            return

        self._guard_mrsf_reference_mode()

        # do TDDFT
        dump_log(self.mol, title='PyOQP: TDDFT steps', section='tdhf')
        self.energy_func[self.td](self.mol)

        # QMRSF dynamic-correlation pathway (additive; runs only on an MRSF backbone)
        qmrsf_pathway = str(self.mol.config['tdhf'].get('qmrsf_pathway', 'none')).strip().lower()
        if qmrsf_pathway != 'none':
            if self.td not in ('mrsf', 'umrsf'):
                raise ValueError('qmrsf_pathway requires an MRSF backbone ([tdhf] type=mrsf)')
            dump_log(self.mol, title='PyOQP: QMRSF pathway (%s)' % qmrsf_pathway, section='tdhf')
            self.energy_func['qmrsf_%s' % qmrsf_pathway](self.mol)


class Gradient(Calculator):
    """
    OQP gradient calculation class

    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.method = mol.config["input"]["method"]
        self.td = mol.config["tdhf"]["type"]
        self.grads = mol.config["properties"]["grad"]
        self.natom = mol.data["natom"]
        self.nstate = mol.config['tdhf']['nstate']
        self.td_prop = mol.config['properties']['td_prop']

        self.zvec_func = {
            'rpa': oqp.tdhf_z_vector,
            'tda': oqp.tdhf_z_vector,
            'sf': oqp.tdhf_sf_z_vector,
            'mrsf': oqp.tdhf_mrsf_z_vector,
        }

        self.grad_func = {
            'hf': oqp.hf_gradient,
            'rpa': oqp.tdhf_gradient,
            'tda': oqp.tdhf_gradient,
            'sf': oqp.tdhf_sf_gradient,
            'mrsf': oqp.tdhf_mrsf_gradient,
        }

    def gradient(self):
        # check method
        if self.method not in ['hf', 'tdhf']:
            raise ValueError(f'Unknown method type {self.method}')

        dump_log(self.mol, title='PyOQP: Entering Gradient Calculation')

        # compute gradients
        grads = []
        if self.method == 'hf':
            grads = self.scf_grad()
        elif self.method == 'tdhf':
            grads = self.tddft_grad()

        # Petite-list runs produce a skeleton two-electron gradient; project
        # onto the totally symmetric component (exact for 1-dim irreps; all
        # abelian irreps are 1-dim). No-op unless the reduction is active.
        grads = self.mol.symmetrize_gradient(grads)

        self.mol.grads = grads

        return grads

    def scf_grad(self):
        dump_log(self.mol, title='PyOQP: Gradient of Root 0')
        self.grad_func['hf'](self.mol)
        grad = self.mol.get_grad()
        grads = np.array([grad.copy()]).reshape((1, self.natom, 3))

        return grads

    def tddft_grad(self):
        if self.td == 'umrsf':
            raise NotImplementedError('UMRSF-TDDFT gradients are not implemented; run UMRSF-TDDFT with runtype=energy only.')
        if self.td not in ['rpa', 'tda', 'sf', 'mrsf']:
            raise ValueError(f'Unknown tdhf type {self.td}')

        if self.nstate < max(self.grads):
            raise ValueError(f'Gradient requested state {max(self.grads)} > the highest computed state {self.nstate}')

        grads = np.zeros((self.nstate + 1, self.natom, 3))
        for i in self.grads:
            dump_log(self.mol, title='PyOQP: Gradient of Root %s' % i)
            self.mol.data.set_tdhf_target(i)
            self.zvec_func[self.td](self.mol)

            # check convergence
            z_flag = self.mol.mol_energy.Z_Vector_converged

            if not z_flag:
                dump_log(self.mol, title='PyOQP: TD Z-vector is not converged', section='end')

                if self.exception is True:
                    raise ZVnotConverged()
                else:
                    exit()
            if self.td_prop == True:
                oqp.electric_moments_excited(self.mol)
                oqp.mulliken_excited(self.mol)

            self.grad_func[self.td](self.mol)
            grad = self.mol.get_grad().reshape((self.natom, 3))
            grads[i] = grad.copy()

        return grads


class Hessian(Calculator):
    """
    OQP frequence calculation class

    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.hess_type = mol.config['hess']['type']
        self.state = mol.config['hess']['state']
        self.read = mol.config['hess']['read']
        self.restart = mol.config['hess']['restart']
        self.temperature = mol.config['hess']['temperature']
        self.clean = mol.config['hess']['clean']

        if self.hess_type == 'analytical':
            self.hess_func = self.analytical_hess
        else:
            self.hess_func = self.numerical_hess

        # Native Hessian ABI placeholders (oqp.hf_hessian, oqp.tdhf_hessian,
        # oqp.tdhf_sf_hessian). These entries are intentionally not used as a
        # numerical fallback while kernels/storage are still being implemented.
        self.native_hess_func = {
            'hf': getattr(oqp, 'hf_hessian', None),
            'rpa': getattr(oqp, 'tdhf_hessian', None),
            'tda': getattr(oqp, 'tdhf_hessian', None),
            'sf': getattr(oqp, 'tdhf_sf_hessian', None),
        }

        method = mol.config['input']['method']
        scf_mult = mol.config['scf']['multiplicity']
        td_mult = mol.config['tdhf']['multiplicity']
        if method == 'hf':
            self.hess_mult = scf_mult
        else:
            self.hess_mult = td_mult

    def _collect_native_fort6_logs(self, mol=None, append_to_log=True):
        """Append and remove Fortran unit-6 scratch logs left by native kernels."""

        mol = mol or self.mol
        native_cphf_logs = []
        for log_dir in (getattr(mol, 'log_path', os.getcwd()), os.getcwd()):
            native_cphf_log = os.path.abspath(os.path.join(log_dir, 'fort.6'))
            if native_cphf_log not in native_cphf_logs:
                native_cphf_logs.append(native_cphf_log)
        for native_cphf_log in native_cphf_logs:
            if not os.path.exists(native_cphf_log):
                continue
            if append_to_log and hasattr(mol, 'log'):
                with open(native_cphf_log, 'r', encoding='utf-8', errors='replace') as source:
                    native_text = source.read()
                if native_text.strip():
                    with open(mol.log, 'a', encoding='utf-8') as target:
                        target.write('\n\n')
                        target.write('PyOQP: Native Fortran HF/DFT analytic Hessian log\n')
                        target.write(native_text)
                        target.write('\n')
            os.remove(native_cphf_log)

    def hessian(self):
        dump_log(self.mol, title='PyOQP: Entering Hessian Calculation')

        if self.read:
            # read .hess file
            dump_log(self.mol, title='', section='read_hess')
            energy, hessian, freqs, modes, inertia = self.mol.read_freqs()

        else:
            # compute hessian
            energy = self.mol.energies[self.state]
            hessian, flags = self.hess_func()
            if 'failed' in flags:
                dump_log(self.mol, title='PyOQP: numerical hessian calculations failed')
                return None
            else:
                freqs, modes, inertia = normal_mode(self.mol.get_system(), self.mol.get_mass(), hessian)
                self.mol.freqs = freqs
                self.mol.hessian = hessian
                self.mol.modes = modes
                self.mol.inertia = inertia
                self._compute_vibrational_intensities(modes)

                # Metadata-only mode irrep labels (no-op unless symmetry enabled).
                if getattr(self.mol, 'symmetry_metadata', None):
                    self.mol.label_normal_modes()

                self.mol.save_freqs(self.state)
                dump_data(self.mol, (self.mol, freqs, modes), title='FREQ', fpath=self.mol.log_path)

                # save mol
                if self.save_mol:
                    self.mol.save_data()

        dump_log(self.mol, title='PyOQP: Frequencies', section='freq', info=freqs)
        dump_log(
            self.mol,
            title='PyOQP: Frequency Normal Mode Eigenvectors',
            section='freq_modes',
            info=(self.mol.get_atoms(), freqs, modes),
        )

        for t in self.temperature:
            thermal_data = thermal_analysis(
                energy=energy,
                atoms=self.mol.get_atoms(),
                mass=self.mol.get_mass(),
                freqs=freqs,
                inertia=inertia,
                temperature=t,
                mult=self.hess_mult,
            )
            dump_log(self.mol, title='PyOQP: Thermochemistry at %-10.2f K' % t, section='thermo', info=thermal_data)

    def _native_property_tensors_at(self, coord_bohr):
        """Return native OpenQP dipole (a.u.) and static polarizability at displaced geometry."""

        from oqp.pyoqp import Runner

        def _config_value(value):
            if isinstance(value, bool):
                return str(value).lower()
            if isinstance(value, (list, tuple)):
                if len(value) == 1 and not isinstance(value[0], (list, tuple)):
                    return str(value[0])
                return ','.join(
                    ' '.join(str(item) for item in entry) if isinstance(entry, (list, tuple)) else str(entry)
                    for entry in value
                )
            return str(value)

        raw_config = copy.deepcopy(self.mol.config)
        config = {
            section: {key: _config_value(value) for key, value in values.items()}
            for section, values in raw_config.items()
        }
        config['input']['runtype'] = 'energy'
        if config.get('guess', {}).get('type') == 'json':
            config['guess']['type'] = 'huckel'
        config['guess']['save_mol'] = 'false'
        project = f"{self.mol.project_name}_vibprop"
        runner = Runner(
            project=project,
            input_dict=config,
            log=os.devnull,
            silent=1,
            usempi=False,
        )
        runner.mol.update_system(np.asarray(coord_bohr, dtype=float).reshape((-1, 3)))
        runner.run()

        dipole = np.zeros(3, dtype=np.float64)
        alpha = np.zeros((3, 3), dtype=np.float64)
        oqp.electric_dipole_au(runner.mol, oqp.ffi.cast("double *", oqp.ffi.from_buffer(dipole)))
        oqp.cphf_static_polarizability(runner.mol, oqp.ffi.cast("double *", oqp.ffi.from_buffer(alpha)))
        self._collect_native_fort6_logs(runner.mol, append_to_log=False)
        return dipole, alpha

    def _compute_vibrational_intensities(self, modes):
        """Compute IR/Raman intensities using native OpenQP property kernels."""

        modes = np.ascontiguousarray(np.asarray(modes, dtype=np.float64))
        coord0 = np.asarray(self.mol.get_system(), dtype=float).reshape((-1, 3))
        ncoord = coord0.size
        if modes.ndim != 2 or modes.shape[1] != ncoord:
            self.mol.vibrational_intensity_metadata = {
                'status': 'failed',
                'reason': f'Expected modes with shape (nmode, {ncoord}), got {modes.shape}',
            }
            return

        displacement = 1.0e-3
        dipole_derivs = np.zeros((3, ncoord), dtype=np.float64)
        polar_derivs = np.zeros((3, 3, ncoord), dtype=np.float64)
        flat0 = coord0.reshape(-1)
        for idx in range(ncoord):
            disp = np.zeros(ncoord, dtype=float)
            disp[idx] = displacement
            try:
                dip_plus, polar_plus = self._native_property_tensors_at((flat0 + disp).reshape(coord0.shape))
                dip_minus, polar_minus = self._native_property_tensors_at((flat0 - disp).reshape(coord0.shape))
            except Exception as exc:
                self.mol.vibrational_intensity_metadata = {
                    'status': 'failed',
                    'backend': 'native_openqp_finite_difference',
                    'reason': str(exc),
                }
                return
            dipole_derivs[:, idx] = (dip_plus - dip_minus) / (2.0 * displacement)
            polar_derivs[:, :, idx] = (polar_plus - polar_minus) / (2.0 * displacement)

        nmode = modes.shape[0]
        ir = np.zeros(nmode, dtype=np.float64)
        mode_dipoles = np.zeros((nmode, 3), dtype=np.float64)
        raman = np.zeros(nmode, dtype=np.float64)
        mode_polars = np.zeros((nmode, 3, 3), dtype=np.float64)
        oqp.vibrational_intensities_native(
            self.mol,
            np.int64(nmode),
            np.int64(ncoord),
            oqp.ffi.cast("double *", oqp.ffi.from_buffer(modes)),
            oqp.ffi.cast("double *", oqp.ffi.from_buffer(dipole_derivs)),
            oqp.ffi.cast("double *", oqp.ffi.from_buffer(polar_derivs)),
            oqp.ffi.cast("double *", oqp.ffi.from_buffer(ir)),
            oqp.ffi.cast("double *", oqp.ffi.from_buffer(mode_dipoles)),
            oqp.ffi.cast("double *", oqp.ffi.from_buffer(raman)),
            oqp.ffi.cast("double *", oqp.ffi.from_buffer(mode_polars)),
        )

        self.mol.infrared_intensities = ir
        self.mol.raman_activities = raman
        self.mol.infrared_mode_dipole_derivatives = mode_dipoles
        self.mol.raman_mode_polarizability_derivatives = mode_polars
        self.mol.vibrational_intensity_metadata = {
            'status': 'computed',
            'backend': 'native_openqp_finite_difference',
            'property_kernels': 'electric_dipole_au,cphf_static_polarizability,vibrational_intensities_native',
            'displacement_bohr': float(displacement),
            'ir_units': 'km/mol',
            'raman_units': 'a.u.',
        }

    def analytical_hess(self):
        method = self.mol.config['input']['method']
        td_type = self.mol.config['tdhf']['type']

        if method == 'hf':
            return self.analytical_ground_state_hess()
        if method == 'tdhf' and td_type in {'tda', 'rpa'}:
            return self.analytical_tddft_hess()
        if method == 'tdhf' and td_type == 'sf':
            return self.analytical_sf_hess()
        if method == 'tdhf' and td_type in {'mrsf', 'umrsf'}:
            return self.analytical_mrsf_hess()
        raise NotImplementedError(
            f"Analytic Hessian is not implemented for method={method}, tdhf.type={td_type}"
        )

    def _spherical_ao_active(self):
        """Return True when the current basis is dimension-reduced by ispher."""
        from oqp.molecule.oqpdata import ispher_mode
        if ispher_mode(self.mol.config.get('input', {}).get('ispher', 'auto')) == 'false':
            return False
        try:
            basis = self.mol.data.get_basis()
            nbf = int(basis['nbf'])
            ncart = int(sum(int((ang + 1) * (ang + 2) // 2) for ang in basis['angs']))
            return nbf != ncart
        except Exception:
            return False

    def analytical_ground_state_hess(self):
        """Run the native OpenQP HF/DFT analytic Hessian kernel and return its stored matrix."""

        native_hess_func = self.native_hess_func['hf']
        if native_hess_func is None:
            raise NotImplementedError('Native OpenQP analytic Hessian entry point oqp.hf_hessian is not available.')
        native_hess_func(self.mol)
        self._collect_native_fort6_logs(self.mol)

        try:
            raw_hessian = self.mol.data['OQP::hf_hessian']
        except (AttributeError, KeyError) as exc:
            raise RuntimeError('Native oqp.hf_hessian did not store OQP::hf_hessian.') from exc

        hessian = self.mol.set_hessian_result(raw_hessian)

        # The native electronic Hessian excludes the empirical dftd4 dispersion
        # term.  The numerical Hessian includes it implicitly (each displaced
        # gradient is dispersion-corrected), so add d2 E_disp / dR2 here to keep
        # the analytic path consistent with the numerical one.
        disp_hessian = self._dispersion_hessian()
        d4_added = np.ndim(disp_hessian) != 0
        if d4_added:
            hessian = hessian + disp_hessian
            self.mol.hessian = hessian

        metadata = dict(getattr(self.mol, 'hessian_metadata', {}) or {})
        metadata.update({
            'backend': 'native_openqp',
            'native_openqp_kernel': True,
            'native_openqp_cphf_solver_exercised': True,
            'native_openqp_final_assembly': True,
            'native_openqp_d4_dispersion': d4_added,
            'no_external_hessian_backend': True,
            'no_numerical_fallback': True,
            'shape': list(hessian.shape),
        })
        setattr(self.mol, 'hessian_metadata', metadata)
        return hessian, ['computed', 'native_openqp']

    def analytical_tddft_hess(self):
        td_type = self.mol.config['tdhf']['type']
        raise NotImplementedError(
            f'TDDFT analytic Hessian is not implemented yet for tdhf.type={td_type}.'
        )

    def analytical_sf_hess(self):
        raise NotImplementedError(
            'SF-TDDFT analytic Hessian is not implemented yet; no numerical fallback will be used.'
        )

    def analytical_mrsf_hess(self):
        td_type = self.mol.config['tdhf']['type']
        label = 'MRSF-TDDFT' if td_type == 'mrsf' else td_type.upper()
        raise NotImplementedError(
            f'{label} analytic Hessian is not implemented yet; no numerical fallback will be used.'
        )

    def _dispersion_hessian(self):
        """D4 dispersion contribution to the analytic Hessian, or 0.0 if disabled.

        dftd4 exposes the dispersion energy and gradient (not a Hessian), so we
        central-difference its analytic gradient with the same step the numerical
        Hessian uses.  dftd4 gradients are cheap (~ms), so the 6N evaluations add
        negligible cost.  Coordinates are in Bohr and the result is in
        Hartree/Bohr**2, matching the native electronic Hessian, with the same
        atom-major (x, y, z) ordering used throughout.
        """
        if not self.mol.config.get('input', {}).get('d4', False):
            return 0.0

        try:
            from dftd4.interface import DampingParam, DispersionModel
        except Exception as exc:  # pragma: no cover - exercised only without dftd4
            raise RuntimeError(
                'hess.type=analytical with input.d4=true requires the dftd4 package; '
                'install dftd4 or use hess.type=numerical.'
            ) from exc

        functional = self.mol.config['input']['functional'].lower() or 'hf'
        func_for_d4 = 'bhlyp' if functional in ('bhhlyp',) else functional
        param = DampingParam(method=func_for_d4)

        atoms = self.mol.get_atoms()
        dx = self.mol.config['hess']['dx']
        flat = np.asarray(self.mol.get_system(), dtype=float).reshape(-1)
        ncoord = flat.size

        def disp_grad(coord_flat):
            model = DispersionModel(atoms, coord_flat.reshape((-1, 3)))
            res = model.get_dispersion(param, grad=True)
            return np.asarray(res['gradient'], dtype=float).reshape(-1)

        hess = np.zeros((ncoord, ncoord))
        for i in range(ncoord):
            cp = flat.copy(); cp[i] += dx
            cm = flat.copy(); cm[i] -= dx
            hess[i, :] = (disp_grad(cp) - disp_grad(cm)) / (2.0 * dx)
        # symmetrize (the FD asymmetry is O(dx**2))
        return 0.5 * (hess + hess.T)

    def numerical_hess(self):
        dir_hess = f'{self.mol.log_path}/{self.mol.project_name}_num_hess'
        nproc = self.mol.config['hess']['nproc']
        dx = self.mol.config['hess']['dx']
        origin_coord = self.mol.get_system()

        # prepare scratch folder
        os.makedirs(dir_hess, exist_ok=True)

        # shift origin 3N coord with 6N displacement
        ncoord = len(origin_coord)
        shift = np.diag(np.ones(ncoord) * dx).reshape(ncoord, ncoord)
        shifted_coord = np.concatenate((origin_coord + shift, origin_coord - shift), axis=0)
        ndim = len(shifted_coord)

        # prepare grad calculations
        self.mol.save_data()
        self.mpi_manager.barrier()
        atoms = self.mol.get_atoms()
        guess_file = self.mol.log.replace('.log', '.json')
        variables_wrapper = [
            {
                'idx': idx,
                'atoms': atoms,
                'coord': coord,
                'dir_hess': dir_hess,
                'project_name': self.mol.project_name,
                'config': copy.deepcopy(self.mol.config),
                'guess_file': guess_file,
                'state': self.state,
                'restart': self.restart,
            }
            for idx, coord in enumerate(shifted_coord)
        ]

        ## adjust multiprocessing if necessary
        if self.mpi_manager.use_mpi:
            ncpu = np.amin([ndim, self.mpi_manager.comm.size])
            pool = MPIPool(processes=ncpu)
        else:
            ncpu = np.amin([ndim, nproc])
            pool = multiprocessing.Pool(processes=ncpu)

        dump_log(self.mol,
                 title='',
                 section='num_hess',
                 info=[self.state, ndim, dx, self.restart, len(variables_wrapper), ncpu, os.environ['OMP_NUM_THREADS']]
                 )

        ## start multiprocessing
        grads = [[] for _ in range(ndim)]
        flags = []
        n = 0
        for val in pool.imap_unordered(grad_wrapper, variables_wrapper):
            if self.mpi_manager.rank == 0:
                n += 1
                idx, grad, flag, timing = val
                grads[idx] = grad
                flags.append(flag)
                dump_log(self.mol, title=None, section='hess_worker', info=[n, idx, flag, timing])
                dump_data(self.mol, (n, idx, np.sum(grad ** 2) ** 2, timing), title='NUM_HESS', fpath=self.mol.log_path)

        pool.close()

        grads = self.mpi_manager.bcast(grads)
        # compute hessian
        forward = np.array(grads[0:ncoord])
        backward = np.array(grads[ncoord:])
        hessian = (forward - backward) / (2 * dx)

        # symmetrize hessian
        hessian = (hessian + hessian.T) / 2

        # delete scratch folder
        if 'failed' not in flags and self.clean and self.mpi_manager.rank == 0:
            shutil.rmtree(dir_hess)

        return hessian, flags


def _run_oqp_external(inp):
    # Run a calculation in a fresh process. Prefer the installed `openqp`
    # console script; fall back to invoking the same entry point
    # (oqp.pyoqp:main) via the current interpreter so this works when OpenQP
    # is run from source without a pip install.
    openqp_exe = shutil.which('openqp')
    if openqp_exe:
        cmd = [openqp_exe, inp, '--silent']
    else:
        cmd = [sys.executable, '-m', 'oqp.pyoqp', inp, '--silent']
    subprocess.run(cmd)


def grad_wrapper(key_dict):
    start_time = time.time()
    rank = MPIManager().rank
    threads = os.environ['OMP_NUM_THREADS']
    host = platform.node()
    # unpack variables
    idx = key_dict['idx']
    atoms = key_dict['atoms']
    coord = key_dict['coord']
    dir_hess = key_dict['dir_hess']
    project_name = key_dict['project_name']
    config = key_dict['config']
    guess_file = key_dict['guess_file']
    state = key_dict['state']
    restart = key_dict['restart']

    # prepare log files
    inp = f'{dir_hess}/{project_name}.{idx}.tmp.inp'
    xyz = f'{dir_hess}/{project_name}.{idx}.tmp.xyz'
    dat = f'{dir_hess}/{project_name}.{idx}.grad_{state}'
    log = f'{dir_hess}/{project_name}.{idx}.tmp.log'

    # attempt to read computed data
    if restart and os.path.exists(dat):
        status = 'loaded'
    else:
        status = 'computed'

        # modify config
        config['input']['runtype'] = 'grad'
        config['input']['system'] = xyz
        config['guess']['type'] = 'json'
        config['guess']['file'] = guess_file
        config['guess']['continue_geom'] = 'false'
        config['properties']['grad'] = config['hess']['state']
        config['properties']['export'] = 'True'
        config['properties']['title'] = f'{project_name}.{idx}'
        config['hess']['temperature'] = ','.join([str(x) for x in config['hess']['temperature']])
        config['tests']['exception'] = 'false'

        # save config
        input_xyz = write_xyz(atoms, coord, [idx])
        input_file, input_dict = write_config(config)

        with open(xyz, 'w') as out:
            out.write(input_xyz)

        with open(inp, 'w') as out:
            out.write(input_file)

        if not MPIManager().use_mpi:
            # run grad calculation externally
            _run_oqp_external(inp)
        else:
            # run grad calculation internally
            start_time = time.time()
            mol = Molecule(project_name, inp, log, silent=1)
            mol.usempi = False
            mol.load_config(input_dict)
            mol.load_data()
            mol.start_time = start_time
            dump_log(mol, title='', section='start')
            mol.data["OQP::log_filename"] = log
            oqp.oqp_banner(mol)
            SinglePoint(mol).energy()
            Gradient(mol).gradient()
            LastStep(mol).compute(mol, grad_list=mol.config['properties']['grad'])
            dump_log(mol, title='', section='end')
    try:
        grad = np.loadtxt(dat).reshape(-1)

    except FileNotFoundError:
        grad = np.zeros_like(coord)
        status = 'failed'

    end_time = time.time()

    return idx, grad, status, (start_time, end_time, rank, threads, host)


class BasisOverlap(Calculator):
    """
    OQP basis overlap calculation class
    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.natom = mol.data["natom"]
        self.nstate = mol.config['tdhf']['nstate']
        self.back_door = self.mol.config['properties']['back_door']
        self.align_type = self.mol.config['nac']['align']
        self.overlap_func = oqp.get_structures_ao_overlap

    def overlap(self):
        # load previous data
        self.load_previous_data()

        # compute basis overlap
        self.overlap_func(self.mol)

        # align mo before tdhf
        if self.align_type != 'no':
            self.align_mo()

    def load_previous_data(self):
        dump_log(self.mol, title='PyOQP: Loading Previous Data')
        if self.back_door:
            # get data externally
            previous_xyz, previous_data = self.mol.get_data_from_back_door()
        else:
            # record data from the current step
            current_file = self.mol.config["guess"]["file"]
            current_xyz = copy.deepcopy(self.mol.get_system())
            current_data = copy.deepcopy(self.mol.get_data())

            # check data from the previous step
            previous_coord = self.mol.data.mol2
            previous_file = self.mol.config['guess']['file2']

            if previous_file:
                self.mol.config['guess']['file'] = previous_file
                self.mol.config['guess']['continue_geom'] = True
                self.mol.load_data()
                previous_xyz = copy.deepcopy(self.mol.get_system())
                previous_data = copy.deepcopy(self.mol.get_data())
            else:
                if len(previous_coord) > 1:
                    # compute data for previous step
                    self.mol.idx = 2
                    self.mol.update_system(previous_coord)
                    oqp.library.ints_1e(self.mol)
                    oqp.library.guess(self.mol)
                    SinglePoint(self.mol).energy()
                    LastStep(self.mol).compute(self.mol)
                    previous_xyz = previous_coord
                    previous_data = copy.deepcopy(self.mol.get_data())
                else:
                    previous_xyz = None
                    previous_data = None
                    exit(f'\nmolecule loads previous step data cannot find [guess] file2 or [input] system2')

            # restore current data
            self.mol.idx = 1
            self.mol.config['guess']['file'] = current_file
            self.mol.update_system(current_xyz)
            self.mol.put_data(current_data)

        # copy previous data to old tags
        natom = self.mol.data["natom"]
        self.mol.data["OQP::xyz_old"] = previous_xyz.reshape((3, natom))
        self.mol.data["OQP::VEC_MO_A_old"] = previous_data["OQP::VEC_MO_A"]
        self.mol.data["OQP::VEC_MO_B_old"] = previous_data["OQP::VEC_MO_B"]
        self.mol.data["OQP::E_MO_A_old"] = previous_data["OQP::E_MO_A"]
        self.mol.data["OQP::E_MO_B_old"] = previous_data["OQP::E_MO_B"]

        # check if td data are available
        try:
            self.mol.data["OQP::td_bvec_mo_old"] = previous_data["OQP::td_bvec_mo"]
            self.mol.data["OQP::td_energies_old"] = previous_data["OQP::td_energies"]
        except KeyError:
            pass

    def align_mo(self):
        dump_log(self.mol, title='PyOQP: Entering Overlap Calculation', section='basis_overlap')
        current_mo = copy.deepcopy(self.mol.data["OQP::VEC_MO_A"])
        current_energy = copy.deepcopy(self.mol.data["OQP::E_MO_A"])
        nocc = self.mol.data['nelec_A']

        # get MO overlap data
        mo_overlap_matrix = self.mol.data["OQP::overlap_mo_non_orthogonal"]

        # current MO in row, previous MO in column
        occ_order, occ_sign = self.find_vec_order(mo_overlap_matrix[:nocc - 2, : nocc - 2])
        somo_order, somo_sign = self.find_vec_order(mo_overlap_matrix[nocc - 2: nocc, nocc - 2: nocc])
        vir_order, vir_sign = self.find_vec_order(mo_overlap_matrix[nocc:, nocc:])

        mo_order = np.concatenate((occ_order, somo_order + nocc - 2, vir_order + nocc))
        mo_sign = np.concatenate((occ_sign, somo_sign, vir_sign)).reshape((-1, 1))

        # apply sign correction
        current_mo *= mo_sign

        # reorder mo if requested
        if self.align_type == 'reorder':
            current_mo = current_mo[np.argsort(mo_order)]
            current_energy = current_energy[np.argsort(mo_order)]

        # update mo
        self.mol.data["OQP::VEC_MO_A"] = current_mo
        self.mol.data["OQP::VEC_MO_B"] = current_mo
        self.mol.data["OQP::E_MO_A"] = current_energy
        self.mol.data["OQP::E_MO_B"] = current_energy

        dump_log(self.mol, title='PyOQP: Aligning MOs')

        # compute new basis overlap
        self.overlap_func(self.mol)

    @staticmethod
    def find_vec_order(overlap_matrix):
        overlap_matrix = copy.deepcopy(overlap_matrix)
        vec_order = []
        vec_sign = []
        for s_i in overlap_matrix:
            # find the matched vec and sign
            n_i = np.argmax(np.abs(s_i))
            p_i = np.sign(s_i)[n_i]

            # record the vec index and sign
            vec_order.append(n_i)
            vec_sign.append(p_i)

            # remove the matched vec
            overlap_matrix[:, n_i] = 0

        return np.array(vec_order), np.array(vec_sign)


class NACME(BasisOverlap):
    """
    Class to calculate Non-Adiabatic Coupling (NAC) matrix elements
    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.dt = mol.config['nac']['dt']
        self.nac_dim = self.nstate

    def align_x(self):
        # use current td data if previous td data is unavailable
        try:
            self.mol.data["OQP::td_bvec_mo_old"]
        except AttributeError:
            self.mol.data["OQP::td_bvec_mo_old"] = self.mol.data["OQP::td_bvec_mo"]
            self.mol.data["OQP::td_energies_old"] = self.mol.data["OQP::td_energies"]

        previous_x = copy.deepcopy(self.mol.data['OQP::td_bvec_mo_old'])
        current_x = copy.deepcopy(self.mol.data['OQP::td_bvec_mo'])
        x_shape = current_x.shape

        # reshape Fortran data into Python style
        current_x = current_x.reshape((x_shape[1], x_shape[0]))
        previous_x = previous_x.reshape((x_shape[1], x_shape[0]))

        # current x in row, previous x in column
        x_overlap_matrix = np.matmul(current_x, previous_x.T)
        x_order, x_sign = self.find_vec_order(x_overlap_matrix)

        # apply sign correction
        current_x *= x_sign.reshape((-1, 1))

        # update x in Fortran data shape
        self.mol.data['OQP::td_bvec_mo'] = current_x.reshape((x_shape[0], x_shape[1]))

        dump_log(self.mol, title='PyOQP: Aligning X amplitudes')

    def nacme(self):
        """
        Calculates the non-adiabatic coupling (NAC) matrix elements
        between the two geometries.

        Currently, adapted only for MRSF-TDDFT approach
        """
        dump_log(self.mol, title='PyOQP: Entering State Overlap Calculation')

        # align X amplitudes
        self.align_x()

        # compute state overlap
        oqp.get_states_overlap(self.mol)
        state_overlap = self.mol.data["OQP::td_states_overlap"]

        # compute time-derivative nac
        dc_matrix = (state_overlap - state_overlap.T) / self.dt
        e_i = np.array(self.mol.energies[1:]).reshape((-1, 1))
        e_j = np.array(self.mol.energies[1:]).reshape((1, -1))
        gap = e_j - e_i
        nac_matrix = dc_matrix * gap
        nac_matrix = nac_matrix[0:self.nac_dim, 0:self.nac_dim]
        self.mol.data["OQP::dc_matrix"] = dc_matrix
        self.mol.data["OQP::nac_matrix"] = nac_matrix

        dump_log(self.mol, title='PyOQP: Non-Adiabatic Coupling Matrix Calculation', section='nacme')
        dump_log(self.mol, title='PyOQP: phase corrected state overlap (s_ij)', section='nacm', info=state_overlap)
        dump_log(self.mol, title='PyOQP: phase corrected derivative coupling (d_ij)', section='nacm', info=dc_matrix)
        dump_log(self.mol, title='PyOQP: state energy gap (e_ji)', section='nacm', info=gap)
        dump_log(self.mol, title='PyOQP: phase corrected non-adiabatic coupling (h_ij)', section='nacm',
                 info=nac_matrix)

        # save corrected data
        self.mol.save_data()
        self.mol.dcm = dc_matrix

        # export data
        if self.export:
            dump_data(self.mol, (nac_matrix, self.export_title), title='NACME', fpath=self.mol.log_path)
            dump_data(self.mol, (dc_matrix, self.export_title), title='DCME', fpath=self.mol.log_path)

        return dc_matrix, nac_matrix


class NAC(Calculator):
    """
    Class to calculate Non-Adiabatic Coupling (NAC) vector
    """

    def __init__(self, mol):
        super().__init__(mol)
        self.mol = mol
        self.nac_type = mol.config["nac"]["type"]
        self.natom = mol.data["natom"]
        self.nstate = mol.config['tdhf']['nstate']
        self.restart = mol.config['nac']['restart']
        self.clean = mol.config['nac']['clean']
        self.nac_states = mol.config['nac']['states']
        self.bp = mol.config['nac']['bp']
        self.nac_func = oqp.get_states_overlap

        if self.nac_type == 'analytical':
            self.nac_func = self.analytical_nac
        else:
            self.nac_func = self.numerical_nac

    def nac(self):
        dump_log(self.mol, title='PyOQP: Entering NAC Vector Calculation')

        # compute nacv
        nacv, dcv, flags = self.nac_func()
        if 'failed' in flags:
            dump_log(self.mol, title='PyOQP: numerical nac calculations failed')
            return None
        else:
            self.mol.nac = nacv

            # save mol
            if self.save_mol:
                self.mol.save_data()

            # export data
            if self.export:
                dump_data(self.mol, (self.mol, nacv, self.export_title), title='NACV', fpath=self.mol.log_path)

        dump_log(self.mol, title='PyOQP: NAC Vector (h_ij)', section='nacv', info=nacv)
        dump_log(self.mol, title='PyOQP: DC Vector (d_ij)', section='dcv', info=dcv)

        # compute gradients for branching plane
        if self.bp:
            dump_log(self.mol, title='PyOQP: Entering Branching Plane Calculation')
            grad = Gradient(self.mol)
            grad.grads = np.unique(self.nac_states)
            grads = grad.gradient()
            for ij in self.nac_states:
                i, j = np.sort(ij)
                g1 = grads[i]
                g2 = grads[j]
                h = nacv[i - 1, j - 1]
                x, y, sx, sy, pitch, tilt, peak, bifu = self.compute_bp(h, g1, g2)
                dump_log(self.mol,
                         title='PyOQP: Final Gradient',
                         section='grad',
                         info={'el': grads, 'd4': np.zeros_like(g1), 'grad_list': [i, j]}
                         )
                dump_log(self.mol,
                         title='PyOQP: Branching Plane Info %s - %s' % (i, j),
                         section='bp',
                         info=(g1, g2, h, x, y, sx, sy, pitch, tilt, peak, bifu)
                         )
                dump_data(self.mol, (self.mol, g1, g2, h, x, y, i, j), title='BP', fpath=self.mol.log_path)

    @staticmethod
    def compute_bp(h, g1, g2):
        s = (g2 + g1) / 2
        g = (g2 - g1) / 2

        tan2b = 2 * np.sum(g * h) / (np.sum(g ** 2) - np.sum(h ** 2))
        b = np.arctan(tan2b) / 2
        sinb = np.sin(b)
        cosb = np.cos(b)
        gt = g * cosb + h * sinb
        ht = h * cosb - g * sinb

        x = gt / la.norm(gt)
        y = ht / la.norm(ht)

        ngt = np.sum(gt * gt)
        nht = np.sum(ht * ht)

        pitch = (0.5 * (ngt + nht)) ** 0.5
        tilt = (ngt - nht) / (ngt + nht)

        sx = np.sum(s * x)
        sy = np.sum(s * y)

        th = ((sx / pitch) ** 2 + (sy / pitch) ** 2) ** 0.5
        ts = np.arctan(sy / sx)

        peak = th ** 2 / (1 - tilt ** 2) * (1 - tilt * np.cos(2 * ts))
        bifu = (th ** 2 / (4 * tilt ** 2)) ** (1 / 3) * \
               (
                       ((1 + tilt) * np.cos(ts) ** 2) ** (1 / 3) +
                       ((1 - tilt) * np.sin(ts) ** 2) ** (1 / 3)
               )
        return x, y, sx, sy, pitch, tilt, peak, bifu

    def analytical_nac(self):
        exit('analytical nac vector is not available yet, choose numerical')

    def numerical_nac(self):
        dir_nacv = f'{self.mol.log_path}/{self.mol.project_name}_num_nacv'
        nproc = self.mol.config['nac']['nproc']
        dx = self.mol.config['nac']['dx']
        origin_coord = self.mol.get_system()

        # prepare scratch folder
        os.makedirs(dir_nacv, exist_ok=True)

        # shift origin 3N coord with 6N displacement
        ncoord = len(origin_coord)
        shift = np.diag(np.ones(ncoord) * dx).reshape(ncoord, ncoord)
        shifted_coord = np.concatenate((origin_coord + shift, origin_coord - shift), axis=0)
        ndim = len(shifted_coord)

        # prepare grad calculations
        self.mol.save_data()
        self.mpi_manager.barrier()
        atoms = self.mol.get_atoms()
        guess_file = self.mol.log.replace('.log', '.json')
        variables_wrapper = [
            {
                'idx': idx,
                'dx': dx,
                'atoms': atoms,
                'coord': coord,
                'nstate': self.nstate,
                'dir_nacv': dir_nacv,
                'project_name': self.mol.project_name,
                'config': copy.deepcopy(self.mol.config),
                'guess_file': guess_file,
                'restart': self.restart,
            }
            for idx, coord in enumerate(shifted_coord)
        ]

        ## adjust multiprocessing if necessary
        if self.mpi_manager.use_mpi:
            ncpu = np.amin([ndim, self.mpi_manager.comm.size])
            pool = MPIPool(processes=ncpu)
        else:
            ncpu = np.amin([ndim, nproc])
            pool = multiprocessing.Pool(processes=ncpu)

        dump_log(self.mol,
                 title='',
                 section='num_nacv',
                 info=[ndim, dx, self.restart, len(variables_wrapper), ncpu, os.environ['OMP_NUM_THREADS']]
                 )

        ## start multiprocessing
        dcm = [[] for _ in range(ndim)]
        flags = []
        n = 0
        for val in pool.imap_unordered(nacme_wrapper, variables_wrapper):
            if self.mpi_manager.rank == 0:
                n += 1
                idx, dcme, flag, timing = val
                dcm[idx] = dcme.reshape(-1)  # 1D nstate x nstate
                flags.append(flag)
                dump_log(self.mol, title=None, section='nacv_worker', info=[n, idx, flag, timing])
                dump_data(self.mol, (n, idx, np.sum(dcme ** 2) ** 2, timing), title='NUM_NACV', fpath=self.mol.log_path)

        pool.close()

        dcm = self.mpi_manager.bcast(dcm)
        # compute nacv (natom x 3, nstate x nstate)
        forward = np.array(dcm[0:ncoord])
        backward = np.array(dcm[ncoord:])
        dcm = (forward - backward) / 2
        e_i = np.array(self.mol.energies[1:]).reshape((1, -1))
        e_j = np.array(self.mol.energies[1:]).reshape((-1, 1))
        gap = e_j - e_i
        np.fill_diagonal(gap, 1)
        gap = gap.reshape((1, -1))
        nacm = dcm * gap

        # reshape matrix -> (nstate x nstate, natom x 3) -> (nstate, nstate, natom, 3)
        dcv = dcm.T.reshape((self.nstate, self.nstate, self.natom, 3))
        nacv = nacm.T.reshape((self.nstate, self.nstate, self.natom, 3))

        # delete scratch folder
        if 'failed' not in flags and self.clean and self.mpi_manager.rank == 0:
            shutil.rmtree(dir_nacv)

        return nacv, dcv, flags


def nacme_wrapper(key_dict):
    start_time = time.time()
    rank = MPIManager().rank
    threads = os.environ['OMP_NUM_THREADS']
    host = platform.node()
    # unpack variables
    idx = key_dict['idx']
    dx = key_dict['dx']
    atoms = key_dict['atoms']
    coord = key_dict['coord']
    nstate = key_dict['nstate']
    dir_nacv = key_dict['dir_nacv']
    project_name = key_dict['project_name']
    config = key_dict['config']
    guess_file = key_dict['guess_file']
    restart = key_dict['restart']

    # prepare log files
    inp = f'{dir_nacv}/{project_name}.{idx}.tmp.inp'
    xyz = f'{dir_nacv}/{project_name}.{idx}.tmp.xyz'
    dat = f'{dir_nacv}/{project_name}.{idx}.dcme'
    log = f'{dir_nacv}/{project_name}.{idx}.tmp.log'

    # attempt to read computed data
    if restart and os.path.exists(dat):
        status = 'loaded'
    else:
        status = 'computed'

        # modify config
        config['input']['runtype'] = 'nacme'
        config['input']['system'] = xyz
        config['guess']['type'] = 'json'
        config['guess']['file'] = guess_file
        config['guess']['file2'] = guess_file
        config['guess']['continue_geom'] = 'false'
        config['properties']['export'] = 'true'
        config['properties']['title'] = f'{project_name}.{idx}'
        config['nac']['dt'] = str(dx)
        config['tests']['exception'] = 'false'

        # save config
        input_xyz = write_xyz(atoms, coord, [idx])
        input_file, input_dict = write_config(config)

        with open(xyz, 'w') as out:
            out.write(input_xyz)

        with open(inp, 'w') as out:
            out.write(input_file)

        if not MPIManager().use_mpi:
            # run nac calculation externally
            _run_oqp_external(inp)
        else:
            # run nac calculation internally
            start_time = time.time()
            mol = Molecule(project_name, inp, log, silent=1)
            mol.usempi = False
            mol.load_config(input_dict)
            mol.load_data()
            mol.start_time = start_time
            dump_log(mol, title='', section='start')
            mol.data["OQP::log_filename"] = log
            oqp.oqp_banner(mol)
            sp = SinglePoint(mol)
            ref_energy = sp.reference()
            BasisOverlap(mol).overlap()
            sp.excitation(ref_energy)
            LastStep(mol).compute(mol)
            NACME(mol).nacme()
            dump_log(mol, title='', section='end')

    try:
        dcme = np.loadtxt(dat).reshape(-1)

    except FileNotFoundError:
        dcme = np.zeros_like(nstate * nstate)
        status = 'failed'

    end_time = time.time()

    return idx, dcme, status, (start_time, end_time, rank, threads, host)


class SCFnotConverged(Exception):
    pass


class TDnotConverged(Exception):
    pass


class ZVnotConverged(Exception):
    pass
