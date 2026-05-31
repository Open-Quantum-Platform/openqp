"""MRSF Hessian oracle helpers.

These helpers are backend-free diagnostics for Gate 3B. They do not enable an
MRSF Hessian by themselves; they provide the root-character checks that the
finite-difference oracle must pass before a displaced gradient/energy can be
used in a Hessian assembly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np


class MRSFRootTrackingError(RuntimeError):
    """Raised when a displaced MRSF calculation cannot be assigned safely."""


@dataclass(frozen=True)
class MRSFRootTrackingDiagnostic:
    requested_state: int
    matched_state: int
    overlap: float
    overlap_gap: float
    energy_gap: float
    root_flip: bool
    ok: bool


@dataclass(frozen=True)
class MRSFHessianStabilityDiagnostic:
    dx: float
    reference_dx: float
    shape: list[int]
    max_abs_delta: float
    rms_delta: float
    max_asymmetry: float
    ok: bool


def physical_root_label(state: int) -> str:
    """Return the OpenQP/MRSF physical root label for a state index."""

    if state <= 0:
        return "high-spin reference"
    return f"physical S{state - 1}"


def _as_2d_vectors(vectors: np.ndarray | list[list[float]], *, name: str) -> np.ndarray:
    array = np.asarray(vectors, dtype=float)
    if array.ndim != 2:
        raise ValueError(f"{name} must be a 2D array of state vectors")
    norms = np.linalg.norm(array, axis=1)
    if np.any(norms <= 0.0):
        raise ValueError(f"{name} contains a zero-norm state vector")
    return array / norms[:, None]


def _record_value(record: Mapping[str, Any], keys: tuple[str, ...], *, missing: str) -> Any:
    for key in keys:
        try:
            if hasattr(record, "get"):
                value = record.get(key)
                if value is not None:
                    return value
            else:
                return record[key]
        except (KeyError, TypeError):
            pass
    raise KeyError(missing)


def extract_mrsf_root_vectors(record: Mapping[str, Any], *, nstate: int | None = None) -> np.ndarray:
    """Extract flattened TD/MRSF state vectors from an OpenQP data record."""

    array = np.asarray(
        _record_value(
            record,
            ("OQP::td_bvec_mo", "td_bvec_mo", "OQP::td_abxc", "td_abxc"),
            missing="OpenQP MRSF root-vector data not found; expected OQP::td_bvec_mo or OQP::td_abxc",
        ),
        dtype=float,
    )

    if array.ndim < 2:
        raise ValueError(f"MRSF root-vector data must have at least 2 dimensions, got {array.shape}")
    vectors = array.reshape((array.shape[0], -1))
    if nstate is not None:
        vectors = vectors[:nstate]
    return _as_2d_vectors(vectors, name="mrsf_root_vectors")


def extract_mrsf_root_energies(record: Mapping[str, Any], *, nstate: int | None = None) -> np.ndarray:
    """Extract TD/MRSF root energies from an OpenQP data record."""

    energies = np.asarray(
        _record_value(
            record,
            ("OQP::td_energies", "td_energies"),
            missing="OpenQP MRSF root energies not found; expected OQP::td_energies",
        ),
        dtype=float,
    ).reshape(-1)
    if nstate is not None:
        energies = energies[:nstate]
    return energies


def assess_root_tracking(
    *,
    requested_state: int,
    reference_vectors: np.ndarray | list[list[float]],
    displaced_vectors: np.ndarray | list[list[float]],
    reference_energies: list[float] | np.ndarray,
    displaced_energies: list[float] | np.ndarray,
    min_overlap: float = 0.85,
    min_overlap_gap: float = 0.05,
    min_energy_gap: float = 1.0e-4,
) -> MRSFRootTrackingDiagnostic:
    """Assign a displaced MRSF root to the requested reference root.

    The requested state index is preserved. The displaced root is assigned by
    absolute vector overlap against the requested reference vector. A mismatch,
    weak/ambiguous overlap, or near-degenerate assigned root fails loudly.
    """

    if requested_state < 1:
        raise ValueError("MRSF Hessian root tracking requires a positive physical state index")

    ref_vec = _as_2d_vectors(reference_vectors, name="reference_vectors")
    disp_vec = _as_2d_vectors(displaced_vectors, name="displaced_vectors")
    ref_e = np.asarray(reference_energies, dtype=float)
    disp_e = np.asarray(displaced_energies, dtype=float)

    nstate = min(ref_vec.shape[0], disp_vec.shape[0], ref_e.size, disp_e.size)
    if requested_state >= nstate:
        raise ValueError(
            f"requested root {requested_state} is outside available roots 0..{nstate - 1}"
        )

    overlaps = np.abs(disp_vec[:nstate] @ ref_vec[requested_state])
    candidate_states = np.arange(1, nstate, dtype=int)
    if candidate_states.size == 0:
        raise ValueError("MRSF Hessian root tracking requires at least one physical MRSF root")
    candidate_overlaps = overlaps[candidate_states]
    order = candidate_states[np.argsort(candidate_overlaps)[::-1]]
    matched_state = int(order[0])
    best_overlap = float(overlaps[matched_state])
    second_overlap = float(overlaps[order[1]]) if order.size > 1 else 0.0
    overlap_gap = best_overlap - second_overlap if order.size > 1 else float("inf")

    assigned_energy_gaps = np.abs(disp_e[:nstate] - disp_e[matched_state])
    nonzero_gaps = assigned_energy_gaps[assigned_energy_gaps > 0.0]
    energy_gap = float(nonzero_gaps.min()) if nonzero_gaps.size else float("inf")

    diagnostic = MRSFRootTrackingDiagnostic(
        requested_state=requested_state,
        matched_state=matched_state,
        overlap=best_overlap,
        overlap_gap=overlap_gap,
        energy_gap=energy_gap,
        root_flip=(matched_state != requested_state),
        ok=(matched_state == requested_state and best_overlap >= min_overlap and overlap_gap >= min_overlap_gap and energy_gap >= min_energy_gap),
    )

    if diagnostic.root_flip:
        raise MRSFRootTrackingError(
            f"MRSF root flip detected: requested root {requested_state} "
            f"({physical_root_label(requested_state)}) matched root {matched_state} "
            f"with overlap {best_overlap:.6f}."
        )
    if best_overlap < min_overlap or overlap_gap < min_overlap_gap:
        raise MRSFRootTrackingError(
            f"MRSF ambiguous vector overlap for requested root {requested_state}: "
            f"best overlap {best_overlap:.6f}, overlap gap {overlap_gap:.6f}."
        )
    if energy_gap < min_energy_gap:
        raise MRSFRootTrackingError(
            f"MRSF ambiguous near-degenerate energy gap for requested root {requested_state}: "
            f"minimum displaced gap {energy_gap:.6e} Hartree."
        )

    return diagnostic


def _as_gradient_block(values, *, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 2:
        raise ValueError(f"{name} must be a 2D gradient block")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values")
    return array


def _as_nuclear_geometry(coords, charges) -> tuple[np.ndarray, np.ndarray]:
    geometry = np.asarray(coords, dtype=float)
    nuclear_charges = np.asarray(charges, dtype=float).reshape(-1)
    if geometry.ndim == 1 and geometry.size % 3 == 0:
        geometry = geometry.reshape((-1, 3))
    if geometry.ndim != 2 or geometry.shape[1] != 3:
        raise ValueError(f"coords must be an (nat, 3) Cartesian array, got {geometry.shape}")
    if nuclear_charges.size != geometry.shape[0]:
        raise ValueError(
            f"charges length {nuclear_charges.size} does not match nat {geometry.shape[0]}"
        )
    if not np.all(np.isfinite(geometry)) or not np.all(np.isfinite(nuclear_charges)):
        raise ValueError("nuclear geometry contains non-finite values")
    return geometry, nuclear_charges


def nuclear_repulsion_energy(coords, charges) -> float:
    """State-independent nuclear repulsion energy for analytical Hessian checks."""

    geometry, nuclear_charges = _as_nuclear_geometry(coords, charges)
    energy = 0.0
    for atom_b in range(1, geometry.shape[0]):
        for atom_a in range(atom_b):
            displacement = geometry[atom_a] - geometry[atom_b]
            distance = float(np.linalg.norm(displacement))
            if distance <= 1.0e-12:
                raise ValueError("nuclear repulsion is singular for coincident atoms")
            energy += nuclear_charges[atom_a] * nuclear_charges[atom_b] / distance
    return float(energy)


def mrsf_nuclear_repulsion_hessian(coords, charges) -> np.ndarray:
    """Analytical nuclear-repulsion Hessian term for MRSF/TDDFT Hessian assembly.

    This term is state independent and is the first safe analytical increment:
    it does not include electronic one-electron, two-electron, XC, or response
    terms and must not be reported as a complete MRSF Hessian.
    """

    geometry, nuclear_charges = _as_nuclear_geometry(coords, charges)
    natom = geometry.shape[0]
    hessian = np.zeros((3 * natom, 3 * natom), dtype=float)
    identity = np.eye(3)

    for atom_b in range(1, natom):
        for atom_a in range(atom_b):
            displacement = geometry[atom_a] - geometry[atom_b]
            distance = float(np.linalg.norm(displacement))
            if distance <= 1.0e-12:
                raise ValueError("nuclear repulsion Hessian is singular for coincident atoms")
            charge_product = nuclear_charges[atom_a] * nuclear_charges[atom_b]
            block = charge_product * (
                3.0 * np.outer(displacement, displacement) / distance**5
                - identity / distance**3
            )
            a_slice = slice(3 * atom_a, 3 * atom_a + 3)
            b_slice = slice(3 * atom_b, 3 * atom_b + 3)
            hessian[a_slice, a_slice] += block
            hessian[b_slice, b_slice] += block
            hessian[a_slice, b_slice] -= block
            hessian[b_slice, a_slice] -= block

    return 0.5 * (hessian + hessian.T)


def assemble_root_tracked_central_fd_hessian(
    *,
    requested_state: int,
    reference_vectors: np.ndarray | list[list[float]],
    reference_energies: list[float] | np.ndarray,
    gradients_plus,
    gradients_minus,
    displaced_vectors,
    displaced_energies,
    dx: float,
    min_overlap: float = 0.85,
    min_overlap_gap: float = 0.05,
    min_energy_gap: float = 1.0e-4,
) -> tuple[np.ndarray, list[MRSFRootTrackingDiagnostic]]:
    """Assemble a central-FD Hessian only after all displaced roots pass tracking.

    ``gradients_plus`` and ``gradients_minus`` are ordered as one flattened
    gradient vector per Cartesian coordinate. ``displaced_vectors`` and
    ``displaced_energies`` must contain the matching plus records followed by
    the matching minus records, mirroring OpenQP's numerical-Hessian worker
    ordering.
    """

    if dx <= 0.0:
        raise ValueError("dx must be positive for central finite differences")

    forward = _as_gradient_block(gradients_plus, name="gradients_plus")
    backward = _as_gradient_block(gradients_minus, name="gradients_minus")
    if forward.shape != backward.shape:
        raise ValueError(
            f"plus/minus gradient blocks must have the same shape, got {forward.shape} and {backward.shape}"
        )

    ncoord = forward.shape[0]
    if len(displaced_vectors) != 2 * ncoord or len(displaced_energies) != 2 * ncoord:
        raise ValueError("displaced root diagnostics must contain plus and minus records for every coordinate")

    diagnostics: list[MRSFRootTrackingDiagnostic] = []
    for vectors, energies in zip(displaced_vectors, displaced_energies):
        diagnostics.append(
            assess_root_tracking(
                requested_state=requested_state,
                reference_vectors=reference_vectors,
                displaced_vectors=vectors,
                reference_energies=reference_energies,
                displaced_energies=energies,
                min_overlap=min_overlap,
                min_overlap_gap=min_overlap_gap,
                min_energy_gap=min_energy_gap,
            )
        )

    hessian = (forward - backward) / (2.0 * dx)
    return (hessian + hessian.T) / 2.0, diagnostics


def compare_fd_step_hessians(
    *,
    dx: float,
    hessian,
    reference_dx: float,
    reference_hessian,
    max_step_delta: float = 5.0e-4,
    max_asymmetry: float = 5.0e-5,
) -> MRSFHessianStabilityDiagnostic:
    """Compare one FD-step Hessian to a reference-step Hessian."""

    matrix = np.asarray(hessian, dtype=float)
    reference = np.asarray(reference_hessian, dtype=float)
    if matrix.shape != reference.shape:
        raise ValueError(f"Hessian step comparison shape mismatch: {matrix.shape} vs {reference.shape}")
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Hessian must be a square 2D matrix, got {matrix.shape}")
    if not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(reference)):
        raise ValueError("Hessian step comparison contains non-finite values")

    delta = matrix - reference
    max_delta = float(np.max(np.abs(delta))) if delta.size else 0.0
    rms_delta = float(np.sqrt(np.mean(delta * delta))) if delta.size else 0.0
    asym = matrix - matrix.T
    asymmetry = float(np.max(np.abs(asym))) if asym.size else 0.0
    return MRSFHessianStabilityDiagnostic(
        dx=float(dx),
        reference_dx=float(reference_dx),
        shape=list(matrix.shape),
        max_abs_delta=max_delta,
        rms_delta=rms_delta,
        max_asymmetry=asymmetry,
        ok=(max_delta <= max_step_delta and asymmetry <= max_asymmetry),
    )
