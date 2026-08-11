"""CASCI energy driver (clean module).

Extracted from the historical 19k-line ``casci.py``: the :class:`CASCI` class plus
the four module helpers it uses, with dependencies trimmed to :mod:`oqp.library.fci`,
:mod:`oqp.library.cas_orbitals` and :mod:`oqp.library.rdm`.  The CASSCF
orbital-optimization scaffolding and the Epstein-Nesbet CASPT2 preflight hook were
retired -- CASSCF/SA-CASSCF now live in :mod:`oqp.library.casscf`, and
CASPT2/NEVPT2/MS-/XMS-CASPT2 in :mod:`oqp.library.caspt2_dyall`.  ``method=casci``
dispatches here and is validated against PySCF (H4/STO-3G CAS(2,2) -2.1147334908).
"""

from __future__ import annotations

import os
import numpy as np

import oqp
from oqp.library.cas_orbitals import load_cas_mo_coeff
from oqp.library.fci import (
    FCI,
    active_space_plan,
    apply_active_space,
    _integer_vector,
    _davidson_subspace_limit,
    _real_array,
    _save_npz_artifact,
    _transform_integrals,
    _unpack_lower_triangle,
    resolve_ci_solve,
    settings_from_casci_config,
)


def _validate_rdm_root_stacks(rdm1_roots: np.ndarray, rdm2_roots: np.ndarray) -> tuple[int, int]:
    if rdm1_roots.ndim != 3 or rdm1_roots.shape[1] != rdm1_roots.shape[2]:
        raise ValueError("RDM1 root stack must have shape (nroot, nactive, nactive)")
    nroot = int(rdm1_roots.shape[0])
    nactive = int(rdm1_roots.shape[1])
    if rdm2_roots.shape != (nroot, nactive, nactive, nactive, nactive):
        raise ValueError("RDM2 root stack shape must match the RDM1 active space")
    return nroot, nactive


def _casscf_mol_data_has(data, key: str) -> bool:
    try:
        data[key]
    except (KeyError, AttributeError, TypeError, ValueError):
        return False
    return True

def _casscf_indices_from_metadata(metadata: dict, key: str) -> tuple[int, ...]:
    raw = str(metadata.get(key, "")).strip()
    if not raw:
        return ()
    values = tuple(int(token.strip()) - 1 for token in raw.split(",") if token.strip())
    if any(index < 0 for index in values):
        raise ValueError(f"{key} metadata must contain positive 1-based orbital indices")
    return values

def _casscf_real_contiguous_copy(array: np.ndarray) -> np.ndarray:
    """Return an owned contiguous real float64 copy of a numeric tensor."""
    real_array = np.real_if_close(np.asarray(array), tol=1000)
    if np.iscomplexobj(real_array):
        raise ValueError("native dependent-state trial tensors must be real")
    return np.ascontiguousarray(np.asarray(real_array, dtype=np.float64).copy())


class CASCI(FCI):
    """CASCI energy calculation using OpenQP-native RHF orbitals and CI solver."""

    method_label = "casci"
    data_prefix = "CASCI"
    log_title = "PyOQP: Complete Active Space Configuration Interaction"
    active_section = "[cas]"
    ci_section = "[ci]"

    def __init__(self, mol):
        self.mol = mol
        self.settings = settings_from_casci_config(mol.config)
        self._native_dependent_state_trial = None

    def _check_combined_ci_memory(self, nbf, plan):
        """Guard the native CI peak together with tensors retained by CASCI."""
        spec = resolve_ci_solve(
            plan.nact,
            plan.nelec,
            ecore=0.0,
            nroot=self.settings.nroot,
            max_det=self.settings.max_det,
            max_memory=self.settings.max_memory,
            eig_tol=self.settings.eig_tol,
            integral_cutoff=self.settings.integral_cutoff,
            solver=self.settings.solver,
            davidson_maxiter=self.settings.davidson_maxiter,
            davidson_subspace=self.settings.davidson_subspace,
            target_spin="any",
            active_section=self.active_section,
            ci_section=self.ci_section,
        )
        nact = int(plan.nact)
        ndet = int(spec.ndet)
        nroot = int(spec.nroot)
        # The AO ERI record, transformed full-MO ERI, dependent-state copy of
        # that MO ERI, and its active-ERI copy all stay live across the native
        # call.  The driver then gathers one more active ERI and expands it to
        # the (2*nact)^4 spin tensor before entering either eigensolver.
        resident = 3 * 8 * int(nbf) ** 4 + 8 * nact ** 4
        integral_work = 8 * nact ** 4 + 8 * (2 * nact) ** 4
        bookkeeping = 8 * ((2 * nact) ** 2 + nact ** 2)
        bookkeeping += 24 * ndet + 8 * ndet * nroot
        if spec.solver == "dense":
            solver_work = 3 * 8 * ndet ** 2
        else:
            effective = (max(spec.davidson_subspace, nroot)
                         if spec.davidson_subspace else 0)
            max_subspace = _davidson_subspace_limit(ndet, nroot, effective)
            solver_work = 8 * ndet * (2 * max_subspace + 4 * nroot)
        peak = resident + integral_work + bookkeeping + solver_work
        budget = int(self.settings.max_memory) * 1024 ** 2
        if peak > budget:
            raise ValueError(
                "CASCI combined integral and CI working set needs ~%.2f GiB "
                "while AO/MO integrals remain resident, exceeding the %s "
                "max_memory budget of %d MiB. Reduce the basis or active "
                "space, raise %s max_memory, or use a smaller CI root window."
                % (peak / 1024 ** 3, self.active_section,
                   int(self.settings.max_memory), self.active_section)
            )

    def energy(self, ref_energy=None):
        energies = super().energy(ref_energy=ref_energy)
        self._store_native_dependent_state_trial_tensors()
        self._store_state_average_rdm_tensors()
        if self.settings.save_rdm:
            self._save_native_rdm_artifact()
        return energies

    def _native_mo_integrals(self):
        nbf = int(self.mol.data.get_basis()["nbf"])
        self._check_ao_eri_budget(nbf)
        oqp.fci_ao_integrals(self.mol)

        hcore = _unpack_lower_triangle(np.asarray(self.mol.data["OQP::Hcore"], dtype=float), nbf)
        default_coeff = np.asarray(self.mol.data["OQP::VEC_MO_A"], dtype=float).reshape((nbf, nbf)).T
        _ovl = _unpack_lower_triangle(
            np.asarray(self.mol.data["OQP::SM"], dtype=float), nbf)
        coeff, source_label = load_cas_mo_coeff(
            self.mol.config, nbf, default_coeff, overlap=_ovl,
            input_dir=os.path.dirname(os.path.abspath(self.mol.input_file or '.')))
        # Commit non-RHF orbitals to the handle.  Molecule.save_data() serializes
        # OQP::VEC_MO_A, so with orbital_source=json and guess.save_mol=true the
        # saved file carried the OLD RHF coefficients -- and feeding that file
        # back as the next orbital_file silently ran a different CASCI.  CASSCF
        # and the PT2 reference already do this; CASCI was the one path left.
        if source_label != "rhf":
            _tgt = np.asarray(self.mol.data["OQP::VEC_MO_A"], dtype=float)
            self.mol.data["OQP::VEC_MO_A"][...] = np.ascontiguousarray(
                np.asarray(coeff, dtype=float).T.reshape(_tgt.shape))
        eri_ao = np.asarray(self.mol.data["OQP::AO_ERI"], dtype=float).reshape(
            (nbf, nbf, nbf, nbf),
            order="F",
        )

        h1e, eri = _transform_integrals(hcore, eri_ao, coeff)
        nelec = (int(self.mol.data["nelec_A"]), int(self.mol.data["nelec_B"]))
        ecore = float(self.mol.mol_energy.nenergy)
        plan = active_space_plan(h1e.shape[0], nelec, self.settings)
        self._check_combined_ci_memory(nbf, plan)
        metadata = dict(plan.metadata)
        metadata["orbital_source"] = source_label
        # The CI solve itself takes the FULL MO integrals and does its own
        # gather/fold natively; the active tensors are still built here, once,
        # for the dependent-state trial payload below.
        h_active, eri_active, active_nelec, ecore_active = apply_active_space(
            h1e, eri, plan, ecore)
        active_indices = plan.active
        if len(active_indices) != h_active.shape[0]:
            raise ValueError("active orbital metadata is inconsistent with active integral shape")
        bare_active_h1e = h1e[np.ix_(active_indices, active_indices)]
        self._native_dependent_state_trial = {
            "h1e": _casscf_real_contiguous_copy(h1e),
            "eri": _casscf_real_contiguous_copy(eri),
            "active_h1e": _casscf_real_contiguous_copy(bare_active_h1e),
            "active_eri": _casscf_real_contiguous_copy(eri_active),
            "active_nelec": tuple(int(value) for value in active_nelec),
            "ecore": float(ecore_active),
            "metadata": dict(metadata),
            "mo_coeff": _casscf_real_contiguous_copy(coeff),
        }
        return h1e, eri, plan, ecore, metadata

    def _store_native_dependent_state_trial_tensors(self):
        state = self._native_dependent_state_trial
        if state is None:
            return

        from oqp.library.rdm import (
            determinant_basis,
            make_rdm1_spatial,
            make_rdm12_spatial_strings,
            make_rdm2_spatial,
            natural_orbital_occupations,
        )

        active_h1e = state["active_h1e"]
        active_nelec = state["active_nelec"]
        ci_vectors = _real_array(
            self.mol.data[f"OQP::{self.data_prefix}_CI_VECTORS"],
            f"{self.data_prefix} CI vectors",
        )
        if ci_vectors.ndim == 1:
            ci_vectors = ci_vectors[:, None]
        elif ci_vectors.ndim != 2 or ci_vectors.shape[1] < 1:
            raise ValueError("CI vectors must have shape (ndet,) or (ndet, nroot)")

        determinants = determinant_basis(active_h1e.shape[0], active_nelec)
        # Bulk per-root record storage: the string-factorized engine (round-off
        # equivalent, orders of magnitude faster at wide active spaces) with
        # the bit-pinned builders as the fallback for non-product lists.
        # These per-root records are built unconditionally, outside [cas]
        # max_memory: nroot copies of an nact^4 RDM2, and np.stack duplicates
        # them while pair_roots still holds the originals, so the transient peak
        # is twice the stack.  CAS(2,30) at nroot=100 clears the ~105 MiB native
        # CI guard and then asks for ~618 MiB of RDM2 and over 1.2 GiB at the
        # peak.  Weigh the peak before building any of it.
        _nact_r = int(active_h1e.shape[0])
        _nroot_r = int(ci_vectors.shape[1])
        _rdm_bytes = 8 * _nroot_r * (_nact_r ** 4 + _nact_r ** 2)
        _rdm_peak = 2 * _rdm_bytes          # pair_roots + the stacked copy
        _budget_r = max(1, int(self.settings.max_memory)) * 1024 ** 2
        # Raise only when the artifacts were actually asked for.  These records
        # are part of the dependent-state trial payload and are built for every
        # CASCI, so an unconditional error turned my own guard into a
        # regression: an energy-only CAS(2,30) at 100 roots fits its CI working
        # set under 256 MiB and was rejected on a stacking peak it never
        # requested.  Deferring the tensors entirely would be better still, but
        # that drops documented OQP::*_TRIAL_RDM* tags from the payload, which
        # is the author's call rather than a review fix.
        _rdm_requested = (self.settings.save_rdm
                          or getattr(self.settings, "state_average_enabled", False))
        if _rdm_requested and _rdm_peak > _budget_r:
            raise ValueError(
                "CASCI per-root RDM records for %d root(s) at nact=%d need "
                "~%.2f GiB (~%.2f GiB at the stacking peak), exceeding the "
                "%s max_memory budget of %d MiB.  Reduce [ci] nroot or the "
                "active space, or raise %s max_memory."
                % (_nroot_r, _nact_r, _rdm_bytes / 1024 ** 3,
                   _rdm_peak / 1024 ** 3, self.active_section,
                   int(self.settings.max_memory), self.active_section))
        pair_roots = [
            make_rdm12_spatial_strings(
                ci_vectors[:, root], determinants, active_h1e.shape[0])
            for root in range(ci_vectors.shape[1])
        ]
        rdm1_roots = np.stack(
            [
                pair_roots[root][0] if pair_roots[root] is not None else
                make_rdm1_spatial(ci_vectors[:, root], determinants, active_h1e.shape[0])
                for root in range(ci_vectors.shape[1])
            ],
            axis=0,
        )
        rdm2_roots = np.stack(
            [
                pair_roots[root][1] if pair_roots[root] is not None else
                make_rdm2_spatial(ci_vectors[:, root], determinants, active_h1e.shape[0])
                for root in range(ci_vectors.shape[1])
            ],
            axis=0,
        )
        natural_occupation_roots = np.stack(
            [
                natural_orbital_occupations(rdm1_roots[root])
                for root in range(ci_vectors.shape[1])
            ],
            axis=0,
        )
        prefix = self.method_label.upper()
        trial_tensors = {
            f"OQP::{prefix}_TRIAL_H1E": state["h1e"],
            f"OQP::{prefix}_TRIAL_ERI": state["eri"],
            f"OQP::{prefix}_TRIAL_ACTIVE_H1E": active_h1e,
            f"OQP::{prefix}_TRIAL_ACTIVE_ERI": state["active_eri"],
            f"OQP::{prefix}_TRIAL_ACTIVE_ORBITALS": np.asarray(
                _casscf_indices_from_metadata(state["metadata"], "active_orbital_indices"),
                dtype=np.int64,
            ),
            f"OQP::{prefix}_TRIAL_CORE_ORBITALS": np.asarray(
                _casscf_indices_from_metadata(state["metadata"], "core_orbital_indices"),
                dtype=np.int64,
            ),
            f"OQP::{prefix}_TRIAL_ACTIVE_NELEC": np.asarray(
                state["active_nelec"],
                dtype=np.int64,
            ),
            f"OQP::{prefix}_TRIAL_NUCLEAR_REPULSION": np.asarray(
                [float(self.mol.mol_energy.nenergy)],
                dtype=np.float64,
            ),
            f"OQP::{prefix}_TRIAL_RDM1": rdm1_roots[0],
            f"OQP::{prefix}_TRIAL_RDM2": rdm2_roots[0],
            f"OQP::{prefix}_TRIAL_RDM1_ROOTS": rdm1_roots,
            f"OQP::{prefix}_TRIAL_RDM2_ROOTS": rdm2_roots,
            f"OQP::{prefix}_TRIAL_NATURAL_OCCUPATIONS": natural_occupation_roots[0],
            f"OQP::{prefix}_TRIAL_NATURAL_OCCUPATIONS_ROOTS": natural_occupation_roots,
        }
        integer_metadata = {
            f"OQP::{prefix}_TRIAL_ACTIVE_ORBITALS",
            f"OQP::{prefix}_TRIAL_CORE_ORBITALS",
            f"OQP::{prefix}_TRIAL_ACTIVE_NELEC",
        }
        for key, value in trial_tensors.items():
            # tagarray rejects zero-byte records (Record.hpp asserts
            # byte_count() != 0), so an empty tensor would abort the process
            # rather than round-trip.  These tensors are published for
            # downstream consumers only, and an absent key is the natural
            # encoding for "none of these orbitals exist" -- e.g. a CAS run
            # with frozen_core=0 has no core orbitals at all.
            if np.asarray(value).size == 0:
                continue
            if key in integer_metadata:
                self.mol.data[key] = np.ascontiguousarray(value, dtype=np.int64)
            else:
                self.mol.data[key] = _casscf_real_contiguous_copy(value)

    def _store_state_average_rdm_tensors(self):
        if not self.settings.state_average_enabled:
            return
        prefix = self.method_label.upper()
        roots_key = f"OQP::{prefix}_STATE_AVERAGE_ROOTS"
        weights_key = f"OQP::{prefix}_STATE_AVERAGE_WEIGHTS"
        if not _casscf_mol_data_has(self.mol.data, roots_key) or not _casscf_mol_data_has(self.mol.data, weights_key):
            return

        roots_source = self.mol.data[roots_key]
        roots = _integer_vector(roots_source, f"{prefix} state-average roots")
        weights = _real_array(self.mol.data[weights_key], f"{prefix} state-average weights")
        rdm1_roots = _real_array(
            self.mol.data[f"OQP::{prefix}_TRIAL_RDM1_ROOTS"],
            f"{prefix} RDM1 root stack",
        )
        rdm2_roots = _real_array(
            self.mol.data[f"OQP::{prefix}_TRIAL_RDM2_ROOTS"],
            f"{prefix} RDM2 root stack",
        )
        if np.asarray(roots_source).ndim != 1 or weights.ndim != 1:
            raise ValueError("State-average roots and weights must be one-dimensional arrays of equal length")
        weights = weights.reshape(-1)
        if weights.shape != roots.shape:
            raise ValueError("State-average roots and weights must be one-dimensional arrays of equal length")
        if roots.size == 0:
            raise ValueError("State-average root list cannot be empty")
        if len(set(int(root) for root in roots)) != roots.size:
            raise ValueError("State-average roots must contain unique roots")
        if np.any(weights < 0.0):
            raise ValueError("State-average weights cannot be negative")
        if float(np.sum(weights)) <= 0.0:
            raise ValueError("State-average weights must contain a positive total weight")
        nroot, _ = _validate_rdm_root_stacks(rdm1_roots, rdm2_roots)
        if np.any(roots < 0) or np.any(roots >= nroot):
            raise ValueError("State-average roots exceed the available RDM root stack")

        state_average_rdm1 = _casscf_real_contiguous_copy(
            np.tensordot(weights, rdm1_roots[roots], axes=(0, 0))
        )
        state_average_rdm2 = _casscf_real_contiguous_copy(
            np.tensordot(weights, rdm2_roots[roots], axes=(0, 0))
        )
        self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_RDM1"] = state_average_rdm1
        self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_RDM2"] = state_average_rdm2

        from oqp.library.rdm import natural_orbital_occupations

        self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_NATURAL_OCCUPATIONS"] = (
            _casscf_real_contiguous_copy(natural_orbital_occupations(state_average_rdm1))
        )

    def _save_native_rdm_artifact(self):
        prefix = self.method_label.upper()
        rdm1_roots = _real_array(
            self.mol.data[f"OQP::{prefix}_TRIAL_RDM1_ROOTS"],
            f"{prefix} RDM1 root stack",
        )
        rdm2_roots = _real_array(
            self.mol.data[f"OQP::{prefix}_TRIAL_RDM2_ROOTS"],
            f"{prefix} RDM2 root stack",
        )
        nroot, nactive = _validate_rdm_root_stacks(rdm1_roots, rdm2_roots)
        natural_occupation_roots = _real_array(
            self.mol.data[f"OQP::{prefix}_TRIAL_NATURAL_OCCUPATIONS_ROOTS"],
            f"{prefix} natural occupation root stack",
        )
        if natural_occupation_roots.shape != (nroot, nactive):
            raise ValueError("Natural occupation root stack shape must match the RDM1 root stack")
        natural_occupations = _real_array(
            self.mol.data[f"OQP::{prefix}_TRIAL_NATURAL_OCCUPATIONS"],
            f"{prefix} natural occupations",
        )
        if natural_occupations.shape != (nactive,):
            raise ValueError("Natural occupation vector shape must match the active space")
        artifact_arrays = {
            "rdm1": np.ascontiguousarray(
                _real_array(self.mol.data[f"OQP::{prefix}_TRIAL_RDM1"], f"{prefix} RDM1"),
                dtype=np.float64,
            ),
            "rdm2": np.ascontiguousarray(
                _real_array(self.mol.data[f"OQP::{prefix}_TRIAL_RDM2"], f"{prefix} RDM2"),
                dtype=np.float64,
            ),
            "rdm1_roots": np.ascontiguousarray(rdm1_roots, dtype=np.float64),
            "rdm2_roots": np.ascontiguousarray(rdm2_roots, dtype=np.float64),
            "natural_occupations": np.ascontiguousarray(
                natural_occupations,
                dtype=np.float64,
            ),
            "natural_occupation_roots": np.ascontiguousarray(
                natural_occupation_roots,
                dtype=np.float64,
            ),
            "energies": np.ascontiguousarray(
                _real_array(
                    self.mol.data[f"OQP::{self.data_prefix}_ENERGIES"],
                    f"{self.data_prefix} energies",
                ),
                dtype=np.float64,
            ),
            "root_indices": np.ascontiguousarray(
                _integer_vector(
                    self.mol.data[f"OQP::{self.data_prefix}_ROOT_INDICES"],
                    f"{self.data_prefix} root indices",
                ),
                dtype=np.int64,
            ),
        }
        if _casscf_mol_data_has(self.mol.data, f"OQP::{prefix}_STATE_AVERAGE_RDM1"):
            state_average_natural_occupations = _real_array(
                self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_NATURAL_OCCUPATIONS"],
                f"{prefix} state-average natural occupations",
            )
            if state_average_natural_occupations.shape != (nactive,):
                raise ValueError(
                    "State-average natural occupation vector shape must match the active space"
                )
            artifact_arrays.update(
                {
                    "state_average_energy": np.ascontiguousarray(
                        _real_array(
                            self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_ENERGY"],
                            f"{prefix} state-average energy",
                        ),
                        dtype=np.float64,
                    ),
                    "state_average_weights": np.ascontiguousarray(
                        _real_array(
                            self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_WEIGHTS"],
                            f"{prefix} state-average weights",
                        ),
                        dtype=np.float64,
                    ),
                    "state_average_roots": np.ascontiguousarray(
                        _integer_vector(
                            self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_ROOTS"],
                            f"{prefix} state-average roots",
                        ),
                        dtype=np.int64,
                    ),
                    "state_average_root_indices": np.ascontiguousarray(
                        _integer_vector(
                            self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_ROOT_INDICES"],
                            f"{prefix} state-average root indices",
                        ),
                        dtype=np.int64,
                    ),
                    "state_average_rdm1": np.ascontiguousarray(
                        _real_array(
                            self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_RDM1"],
                            f"{prefix} state-average RDM1",
                        ),
                        dtype=np.float64,
                    ),
                    "state_average_rdm2": np.ascontiguousarray(
                        _real_array(
                            self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_RDM2"],
                            f"{prefix} state-average RDM2",
                        ),
                        dtype=np.float64,
                    ),
                    "state_average_natural_occupations": np.ascontiguousarray(
                        state_average_natural_occupations,
                        dtype=np.float64,
                    ),
                }
            )
        _save_npz_artifact(self.mol, f"{self.method_label}_rdms", **artifact_arrays)
