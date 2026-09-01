"""State-tracked MRSF excited-state IR and Raman property derivatives.

This module does not introduce a determinant representation of MRSF.  It
consumes OpenQP's spin-adapted two-SOMO response vectors and state-interaction
one-particle densities.  The four orbital classes ``CO``, ``OV``, ``CV``, and
``OO`` are retained when response-vector weights are recorded for every
displacement.

The permanent-dipole path is a central nuclear finite difference of the full
target-state dipole.  At each geometry the nuclear moment is combined directly
with the electronic contraction of the complete MRSF state 1-RDM.  It is
therefore a target-state quantity, not the ROHF dipole relabelled as MRSF.

OpenQP currently has no uniform-electric-field Hamiltonian input for an MRSF
finite-field polarizability.  The only implemented Raman backend is consequently
an explicitly truncated, convergence-gated sum over MRSF states (SOS).  It is
useful for method development, but is not an analytic or complete MRSF
polarizability.  A requested finite-field backend fails closed.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import hashlib
import json
import os
import tempfile
from typing import Callable, Literal, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .vibronic import (
    ExcitedStatePropertyDerivative,
    InfraredResult,
    RamanResult,
    excited_state_ir_intensities,
    excited_state_raman_activities,
    finite_difference_excited_state_property,
)


FloatArray = NDArray[np.float64]
PolarizabilityBackend = Literal["truncated_sos", "finite_field"]

SCHEMA_VERSION = 1
SPIN_ADAPTED_LAYOUT = "two_somo_spin_adapted_CO_OV_CV_OO"


def _finite_real(name: str, values: ArrayLike) -> FloatArray:
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _sha256_json(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), default=str
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True)
class MRSFPropertyFDRequest:
    """Validated normal-coordinate finite-difference request.

    ``normal_modes`` has shape ``(nmode, 3*natom)`` and follows OpenQP's
    mass-normalized convention, so ``R(Q +/- h) = R(0) +/- h * mode`` when
    ``h`` is in ``sqrt(amu)*bohr``.
    """

    electronic_method: str
    electronic_state: str
    state_index: int
    normal_modes: FloatArray
    displacement: FloatArray
    coordinate_phase_convention: str
    minimum_state_overlap: float = 0.99
    minimum_tracking_margin: float = 0.05
    fd_relative_tolerance: float = 0.05
    fd_absolute_tolerance: float = 1.0e-6
    polarizability_backend: PolarizabilityBackend = "truncated_sos"
    sos_tail_states: int = 2
    sos_tail_relative_tolerance: float = 0.05
    sos_minimum_gap_hartree: float = 1.0e-5

    @classmethod
    def create(
        cls,
        *,
        electronic_method: str,
        electronic_state: str,
        state_index: int,
        normal_modes: ArrayLike,
        displacement: ArrayLike | float,
        coordinate_phase_convention: str,
        minimum_state_overlap: float = 0.99,
        minimum_tracking_margin: float = 0.05,
        fd_relative_tolerance: float = 0.05,
        fd_absolute_tolerance: float = 1.0e-6,
        polarizability_backend: PolarizabilityBackend = "truncated_sos",
        sos_tail_states: int = 2,
        sos_tail_relative_tolerance: float = 0.05,
        sos_minimum_gap_hartree: float = 1.0e-5,
    ) -> "MRSFPropertyFDRequest":
        modes = _finite_real("normal_modes", normal_modes)
        if modes.ndim != 2 or modes.shape[0] == 0 or modes.shape[1] == 0:
            raise ValueError("normal_modes must have shape (nmode, 3*natom)")
        if modes.shape[1] % 3:
            raise ValueError("normal_modes Cartesian dimension must be divisible by 3")
        step = _finite_real("displacement", displacement)
        if step.ndim == 0:
            step = np.full(modes.shape[0], float(step))
        if step.shape != (modes.shape[0],) or np.any(step <= 0.0):
            raise ValueError("displacement must be positive, scalar or one per mode")
        if electronic_method not in ("MRSF-TDHF", "MRSF-TDDFT"):
            raise ValueError("electronic_method must be MRSF-TDHF or MRSF-TDDFT")
        if not electronic_state.strip() or not coordinate_phase_convention.strip():
            raise ValueError(
                "electronic_state and coordinate_phase_convention must be nonempty"
            )
        if int(state_index) < 1:
            raise ValueError(
                "state_index is the one-based MRSF response root and must be >= 1"
            )
        if not 0.99 <= minimum_state_overlap <= 1.0:
            raise ValueError(
                "minimum_state_overlap must lie in [0.99, 1]; weaker isolated-root "
                "tracking requires a projector treatment"
            )
        if minimum_tracking_margin < 0.05:
            raise ValueError(
                "minimum_tracking_margin must be at least 0.05; weaker isolated-root "
                "tracking requires a projector treatment"
            )
        if fd_relative_tolerance <= 0.0 or fd_absolute_tolerance <= 0.0:
            raise ValueError("finite-difference convergence tolerances must be positive")
        if polarizability_backend not in ("truncated_sos", "finite_field"):
            raise ValueError("unsupported polarizability_backend")
        if int(sos_tail_states) < 1:
            raise ValueError("sos_tail_states must be positive")
        if sos_tail_relative_tolerance <= 0.0:
            raise ValueError("sos_tail_relative_tolerance must be positive")
        if sos_minimum_gap_hartree <= 0.0:
            raise ValueError("sos_minimum_gap_hartree must be positive")
        return cls(
            electronic_method=electronic_method,
            electronic_state=electronic_state,
            state_index=int(state_index),
            normal_modes=modes.copy(),
            displacement=step.copy(),
            coordinate_phase_convention=coordinate_phase_convention,
            minimum_state_overlap=float(minimum_state_overlap),
            minimum_tracking_margin=float(minimum_tracking_margin),
            fd_relative_tolerance=float(fd_relative_tolerance),
            fd_absolute_tolerance=float(fd_absolute_tolerance),
            polarizability_backend=polarizability_backend,
            sos_tail_states=int(sos_tail_states),
            sos_tail_relative_tolerance=float(sos_tail_relative_tolerance),
            sos_minimum_gap_hartree=float(sos_minimum_gap_hartree),
        )

    @classmethod
    def from_dict(cls, record: Mapping[str, object]) -> "MRSFPropertyFDRequest":
        version = int(record.get("schema_version", 0))
        if version != SCHEMA_VERSION:
            raise ValueError(
                f"schema_version must be {SCHEMA_VERSION}, got {version}"
            )
        if record.get("electronic_method") not in ("MRSF-TDHF", "MRSF-TDDFT"):
            raise ValueError("electronic_method must be MRSF-TDHF or MRSF-TDDFT")
        if record.get("electronic_state_role") != "target_excited_state":
            raise ValueError(
                "electronic_state_role must be 'target_excited_state'"
            )
        if record.get("response_representation") != SPIN_ADAPTED_LAYOUT:
            raise ValueError(
                "response_representation must declare the spin-adapted "
                "two-SOMO CO/OV/CV/OO layout"
            )
        if record.get("coordinate_basis") != "excited_normal":
            raise ValueError("coordinate_basis must be 'excited_normal'")
        if record.get("coordinate_unit") != "sqrt(amu)*bohr":
            raise ValueError("coordinate_unit must be 'sqrt(amu)*bohr'")
        step_scales = _finite_real(
            "finite_difference_step_scales",
            record.get("finite_difference_step_scales", [1.0, 0.5, 0.25]),
        )
        if step_scales.shape != (3,) or not np.allclose(
            step_scales, [1.0, 0.5, 0.25], rtol=0.0, atol=0.0
        ):
            raise ValueError(
                "finite_difference_step_scales must be exactly [1.0, 0.5, 0.25]"
            )
        return cls.create(
            electronic_method=str(record["electronic_method"]),
            electronic_state=str(record["electronic_state"]),
            state_index=int(record["state_index"]),
            normal_modes=record["normal_modes"],
            displacement=record["displacement"],
            coordinate_phase_convention=str(
                record["coordinate_phase_convention"]
            ),
            minimum_state_overlap=float(
                record.get("minimum_state_overlap", 0.99)
            ),
            minimum_tracking_margin=float(
                record.get("minimum_tracking_margin", 0.05)
            ),
            fd_relative_tolerance=float(
                record.get("fd_relative_tolerance", 0.05)
            ),
            fd_absolute_tolerance=float(
                record.get("fd_absolute_tolerance", 1.0e-6)
            ),
            polarizability_backend=str(
                record.get("polarizability_backend", "truncated_sos")
            ),  # type: ignore[arg-type]
            sos_tail_states=int(record.get("sos_tail_states", 2)),
            sos_tail_relative_tolerance=float(
                record.get("sos_tail_relative_tolerance", 0.05)
            ),
            sos_minimum_gap_hartree=float(
                record.get("sos_minimum_gap_hartree", 1.0e-5)
            ),
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": SCHEMA_VERSION,
            "electronic_method": self.electronic_method,
            "electronic_state": self.electronic_state,
            "electronic_state_role": "target_excited_state",
            "state_index": self.state_index,
            "response_representation": SPIN_ADAPTED_LAYOUT,
            "coordinate_basis": "excited_normal",
            "coordinate_unit": "sqrt(amu)*bohr",
            "coordinate_phase_convention": self.coordinate_phase_convention,
            "normal_modes": self.normal_modes.tolist(),
            "displacement": self.displacement.tolist(),
            "minimum_state_overlap": self.minimum_state_overlap,
            "minimum_tracking_margin": self.minimum_tracking_margin,
            "finite_difference_step_scales": [1.0, 0.5, 0.25],
            "fd_relative_tolerance": self.fd_relative_tolerance,
            "fd_absolute_tolerance": self.fd_absolute_tolerance,
            "polarizability_backend": self.polarizability_backend,
            "sos_tail_states": self.sos_tail_states,
            "sos_tail_relative_tolerance": self.sos_tail_relative_tolerance,
            "sos_minimum_gap_hartree": self.sos_minimum_gap_hartree,
        }


@dataclass(frozen=True)
class MRSFTrackedPropertySnapshot:
    """One displaced MRSF target-state evaluation in the central-state gauge."""

    state_dipole_au: FloatArray
    polarizability_au: FloatArray | None
    matched_overlap: float
    tracking_margin: float
    phase_to_central: float
    raw_root_index: int
    state_energy_hartree: float
    target_multiplicity: int
    expected_s2: float
    state_irrep: str | None
    response_block_weights: Mapping[str, float]
    polarizability_diagnostics: Mapping[str, object]
    provenance: Mapping[str, object]

    @classmethod
    def create(
        cls,
        *,
        state_dipole_au: ArrayLike,
        polarizability_au: ArrayLike | None,
        matched_overlap: float,
        tracking_margin: float,
        phase_to_central: float,
        raw_root_index: int,
        state_energy_hartree: float,
        target_multiplicity: int,
        expected_s2: float,
        state_irrep: str | None,
        response_block_weights: Mapping[str, float],
        polarizability_diagnostics: Mapping[str, object] | None = None,
        provenance: Mapping[str, object] | None = None,
    ) -> "MRSFTrackedPropertySnapshot":
        dipole = _finite_real("state_dipole_au", state_dipole_au)
        if dipole.shape != (3,):
            raise ValueError("state_dipole_au must have shape (3,)")
        polar = None
        if polarizability_au is not None:
            polar = _finite_real("polarizability_au", polarizability_au)
            if polar.shape != (3, 3):
                raise ValueError("polarizability_au must have shape (3, 3)")
            residual = float(np.max(np.abs(polar - polar.T), initial=0.0))
            if residual > 1.0e-10 * max(
                1.0, float(np.max(np.abs(polar), initial=0.0))
            ):
                raise ValueError("polarizability_au must be symmetric")
        if not np.isfinite([matched_overlap, tracking_margin, phase_to_central]).all():
            raise ValueError("tracking diagnostics must be finite")
        if not 0.0 <= matched_overlap <= 1.0 + 1.0e-8:
            raise ValueError("matched_overlap must lie in [0, 1]")
        if tracking_margin < 0.0:
            raise ValueError("tracking_margin must be non-negative")
        if abs(abs(phase_to_central) - 1.0) > 1.0e-10:
            raise ValueError("phase_to_central must be a real unit phase")
        if not np.isfinite([state_energy_hartree, expected_s2]).all():
            raise ValueError("state energy and expected S^2 must be finite")
        if int(target_multiplicity) not in (1, 3):
            raise ValueError("target_multiplicity must be singlet or triplet")
        required_s2 = 0.0 if int(target_multiplicity) == 1 else 2.0
        if abs(expected_s2 - required_s2) > 1.0e-12:
            raise ValueError("expected_s2 is inconsistent with target_multiplicity")
        if state_irrep is not None and not str(state_irrep).strip():
            raise ValueError("state_irrep must be nonempty when supplied")
        keys = set(response_block_weights)
        if keys != {"CO", "OV", "CV", "OO"}:
            raise ValueError("response_block_weights must contain CO, OV, CV, and OO")
        weights = {name: float(response_block_weights[name]) for name in sorted(keys)}
        if not np.isfinite(list(weights.values())).all() or min(weights.values()) < 0.0:
            raise ValueError("response block weights must be finite and non-negative")
        weight_sum = sum(weights.values())
        if abs(weight_sum - 1.0) > 1.0e-8:
            raise ValueError("response block weights must be normalized to one")
        if int(raw_root_index) < 0:
            raise ValueError("raw_root_index must be non-negative")
        return cls(
            state_dipole_au=dipole.copy(),
            polarizability_au=None if polar is None else polar.copy(),
            matched_overlap=float(matched_overlap),
            tracking_margin=float(tracking_margin),
            phase_to_central=float(phase_to_central),
            raw_root_index=int(raw_root_index),
            state_energy_hartree=float(state_energy_hartree),
            target_multiplicity=int(target_multiplicity),
            expected_s2=float(expected_s2),
            state_irrep=None if state_irrep is None else str(state_irrep),
            response_block_weights=weights,
            polarizability_diagnostics=dict(polarizability_diagnostics or {}),
            provenance=dict(provenance or {}),
        )


@dataclass(frozen=True)
class MRSFSpectroscopyFDResult:
    dipole_derivative: ExcitedStatePropertyDerivative
    polarizability_derivative: ExcitedStatePropertyDerivative | None
    infrared: InfraredResult
    raman: RamanResult | None
    displacement_records: tuple[Mapping[str, object], ...]
    provenance: Mapping[str, object]

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": SCHEMA_VERSION,
            "method": (
                "state-tracked MRSF normal-coordinate central differences "
                "at h, h/2, and h/4"
            ),
            "dipole_derivative": self.dipole_derivative.to_dict(),
            "polarizability_derivative": None
            if self.polarizability_derivative is None
            else self.polarizability_derivative.to_dict(),
            "infrared": self.infrared.to_dict(),
            "raman": None if self.raman is None else self.raman.to_dict(),
            "displacements": [dict(item) for item in self.displacement_records],
            "provenance": dict(self.provenance),
        }


def response_block_weights(states: object, root_index: int) -> dict[str, float]:
    """Return normalized spin-adapted CO/OV/CV/OO amplitude weights."""

    nalpha = int(getattr(states, "na"))
    nbeta = int(getattr(states, "nb"))
    if nalpha - nbeta != 2:
        raise ValueError("ordinary MRSF spectroscopy requires exactly two SOMOs")
    amplitude = _finite_real(
        "spin-adapted MRSF amplitude",
        getattr(states, "amplitude_matrix")(int(root_index)),
    )
    if amplitude.ndim != 2 or amplitude.shape[0] != nalpha:
        raise ValueError("unexpected MRSF amplitude layout")
    if amplitude.shape[1] < 2:
        raise ValueError("MRSF amplitude layout is missing the two open orbitals")
    closed = slice(0, nbeta)
    open_rows = slice(nbeta, nalpha)
    open_columns = slice(0, 2)
    virtual = slice(2, amplitude.shape[1])
    raw = {
        "CO": float(np.sum(np.abs(amplitude[closed, open_columns]) ** 2)),
        "OV": float(np.sum(np.abs(amplitude[open_rows, virtual]) ** 2)),
        "CV": float(np.sum(np.abs(amplitude[closed, virtual]) ** 2)),
        "OO": float(np.sum(np.abs(amplitude[open_rows, open_columns]) ** 2)),
    }
    total = sum(raw.values())
    if total <= 1.0e-14:
        raise ValueError("MRSF response vector has zero CO/OV/CV/OO norm")
    return {key: value / total for key, value in raw.items()}


def full_mrsf_state_dipole(states: object, root_index: int, mol: object) -> FloatArray:
    """Full target-state dipole from nuclei and the complete MRSF state 1-RDM.

    ``MRSFExcitedStates.R`` contains AO position integrals about the molecular
    center of mass.  The nuclear term is evaluated about that same origin,
    including electrons removed by an effective core potential.  No generic
    ROHF permanent-dipole kernel is called.
    """

    state_density = _finite_real(
        "MRSF state density",
        getattr(states, "state_density_ao")(int(root_index)),
    )
    dipole_integrals = getattr(states, "R")
    if len(dipole_integrals) != 3:
        raise ValueError("MRSF dipole-integral record must have three components")
    overlap = _finite_real("AO overlap", getattr(states, "S"))
    if state_density.shape != overlap.shape or state_density.ndim != 2:
        raise ValueError("MRSF state density and AO overlap shapes differ")
    electron_count = float(np.sum(state_density * overlap.T))
    expected_electrons = float(getattr(states, "n_elec"))
    if abs(electron_count - expected_electrons) > 1.0e-6:
        raise ValueError(
            f"complete MRSF state density has AO-metric trace {electron_count:.8f}, "
            f"expected {expected_electrons:.8f}"
        )

    geometry = _finite_real("OpenQP geometry", getattr(mol, "get_system")())
    atoms = _finite_real("OpenQP nuclear charges", getattr(mol, "get_atoms")())
    masses = _finite_real("OpenQP atomic masses", getattr(mol, "get_mass")())
    if geometry.size != 3 * atoms.size or masses.shape != atoms.shape:
        raise ValueError("OpenQP geometry, atoms, and masses are inconsistent")
    geometry = geometry.reshape((-1, 3))
    if np.any(masses <= 0.0):
        raise ValueError("atomic masses must be positive")
    charge = int(getattr(mol, "config", {}).get("input", {}).get("charge", 0))
    total_ecp_electrons = float(np.sum(atoms) - charge - expected_electrons)
    if (
        total_ecp_electrons < -1.0e-8
        or abs(total_ecp_electrons - round(total_ecp_electrons)) > 1.0e-8
    ):
        raise ValueError("nuclei, charge, and MRSF electron count are inconsistent")
    total_ecp_electrons = float(round(total_ecp_electrons))
    try:
        exposed_ecp = _finite_real(
            "OpenQP ECP core electrons", getattr(mol, "data")["ecp_zn"]
        ).reshape(-1)
    except (AttributeError, KeyError, TypeError):
        exposed_ecp = np.zeros(0)
    if exposed_ecp.shape == atoms.shape:
        ecp_electrons = exposed_ecp
    elif total_ecp_electrons == 0.0:
        # An all-electron calculation has an unambiguous effective charge even
        # if the optional zero-valued ECP tag was not retained.
        ecp_electrons = np.zeros_like(atoms)
    else:
        raise RuntimeError(
            "OpenQP did not expose one ECP core-electron count per atom; the "
            "target-state nuclear dipole cannot be formed safely"
        )
    if abs(float(np.sum(ecp_electrons)) - total_ecp_electrons) > 1.0e-8:
        raise ValueError("ECP core-electron counts do not match the electron count")
    effective_charges = atoms - ecp_electrons
    if np.any(effective_charges < 0.0):
        raise ValueError("ECP core-electron count exceeds a nuclear charge")
    center_of_mass = np.sum(geometry * masses[:, None], axis=0) / np.sum(masses)
    nuclear = np.sum(
        effective_charges[:, None] * (geometry - center_of_mass), axis=0
    )
    electronic = np.array(
        [-np.sum(state_density * np.asarray(component, dtype=float))
         for component in dipole_integrals],
        dtype=float,
    )
    result = nuclear + electronic
    if not np.all(np.isfinite(result)):
        raise ValueError("target-state MRSF dipole is non-finite")
    return result


def truncated_sos_polarizability(
    states: object,
    root_index: int,
    *,
    tail_states: int,
    tail_relative_tolerance: float,
    minimum_gap_hartree: float,
) -> tuple[FloatArray, dict[str, object]]:
    """Static target-state polarizability from a convergence-gated MRSF SOS.

    The finite state ladder is not a complete response calculation.  The
    result is accepted only when removing the highest ``tail_states`` changes
    the tensor by no more than ``tail_relative_tolerance`` in Frobenius norm.
    """

    energies = _finite_real("MRSF excitation energies", getattr(states, "energies"))
    if energies.ndim != 1:
        raise ValueError("MRSF excitation energies must be one-dimensional")
    root = int(root_index)
    if not 0 <= root < energies.size:
        raise IndexError("target MRSF root is outside the calculated state ladder")
    if energies.size < tail_states + 3:
        raise ValueError(
            "truncated SOS Raman requires at least tail_states + 3 MRSF roots"
        )
    if root >= energies.size - tail_states:
        raise ValueError(
            "target root lies in the SOS tail; calculate more MRSF roots"
        )

    terms: list[tuple[int, FloatArray]] = []
    minimum_observed_gap = np.inf
    for other in range(energies.size):
        if other == root:
            continue
        gap = float(energies[other] - energies[root])
        minimum_observed_gap = min(minimum_observed_gap, abs(gap))
        if abs(gap) < minimum_gap_hartree:
            raise ValueError(
                f"MRSF SOS polarizability is singular: roots {root + 1} and "
                f"{other + 1} differ by {abs(gap):.3e} Hartree"
            )
        transition = _finite_real(
            f"MRSF transition dipole {root}->{other}",
            getattr(states, "transition_dipole")(root, other),
        )
        if transition.shape != (3,):
            raise ValueError("MRSF transition dipoles must have shape (3,)")
        terms.append((other, 2.0 * np.outer(transition, transition) / gap))

    polarizability = sum((tensor for _, tensor in terms), np.zeros((3, 3)))
    tail_start = energies.size - tail_states
    tail = sum(
        (tensor for other, tensor in terms if other >= tail_start),
        np.zeros((3, 3)),
    )
    relative_tail = float(
        np.linalg.norm(tail) / max(np.linalg.norm(polarizability), 1.0e-14)
    )
    if relative_tail > tail_relative_tolerance:
        raise ValueError(
            f"truncated MRSF SOS polarizability is not converged: tail relative "
            f"norm {relative_tail:.3e} exceeds {tail_relative_tolerance:.3e}; "
            "calculate more MRSF roots or use a future finite-field backend"
        )
    polarizability = 0.5 * (polarizability + polarizability.T)
    return polarizability, {
        "backend": "truncated_mrsf_state_sos",
        "analytic": False,
        "finite_field": False,
        "number_of_states": int(energies.size),
        "tail_states": int(tail_states),
        "tail_relative_norm": relative_tail,
        "tail_relative_tolerance": float(tail_relative_tolerance),
        "minimum_gap_hartree": float(minimum_observed_gap),
        "scope": (
            "finite-state spin-adapted MRSF SOS; not a complete static "
            "polarizability and not an analytic MRSF response"
        ),
    }


def assemble_mrsf_spectroscopy_derivatives(
    request: MRSFPropertyFDRequest,
    evaluator: Callable[[int, int, FloatArray], MRSFTrackedPropertySnapshot],
    equilibrium_geometry_bohr: ArrayLike,
    *,
    electronic_method: str,
    source_provenance: Mapping[str, object] | None = None,
) -> MRSFSpectroscopyFDResult:
    """Run plus/minus normal-coordinate evaluations and assemble IR/Raman."""

    if electronic_method not in ("MRSF-TDHF", "MRSF-TDDFT"):
        raise ValueError("electronic_method must be MRSF-TDHF or MRSF-TDDFT")
    if electronic_method != request.electronic_method:
        raise ValueError("request and evaluated electronic methods differ")
    geometry = _finite_real("equilibrium_geometry_bohr", equilibrium_geometry_bohr)
    if geometry.size != request.normal_modes.shape[1]:
        raise ValueError("geometry and normal_modes Cartesian dimensions differ")
    geometry = geometry.reshape(-1)
    nmode = request.normal_modes.shape[0]
    records: list[Mapping[str, object]] = []
    have_polar = request.polarizability_backend == "truncated_sos"
    dipole_derivatives: list[ExcitedStatePropertyDerivative] = []
    polarizability_derivatives: list[
        ExcitedStatePropertyDerivative | None
    ] = []
    step_scales = (1.0, 0.5, 0.25)

    if request.polarizability_backend == "finite_field":
        raise NotImplementedError(
            "OpenQP does not yet expose a uniform-electric-field Hamiltonian "
            "for MRSF. finite_field Raman is unavailable; no ROHF/CPHF "
            "polarizability will be substituted."
        )

    phase_label = "OpenQP response-overlap aligned to the central MRSF root"

    for step_scale in step_scales:
        step = request.displacement * step_scale
        plus_dipole = np.zeros((nmode, 3))
        minus_dipole = np.zeros((nmode, 3))
        plus_polar = np.zeros((nmode, 3, 3))
        minus_polar = np.zeros((nmode, 3, 3))
        overlap = np.zeros((nmode, 2))
        polar_complete = have_polar
        for mode in range(nmode):
            mode_records = []
            for side_index, sign in enumerate((+1, -1)):
                displaced = geometry + sign * step[mode] * request.normal_modes[mode]
                snapshot = evaluator(mode, sign, displaced.reshape((-1, 3)))
                if snapshot.matched_overlap < request.minimum_state_overlap:
                    raise ValueError(
                        f"mode {mode + 1} scale {step_scale:g} {sign:+d}: state "
                        f"overlap {snapshot.matched_overlap:.6f} is below "
                        f"{request.minimum_state_overlap:.6f}"
                    )
                if snapshot.tracking_margin < request.minimum_tracking_margin:
                    raise ValueError(
                        f"mode {mode + 1} scale {step_scale:g} {sign:+d}: "
                        f"root-tracking margin {snapshot.tracking_margin:.6f} is "
                        f"below {request.minimum_tracking_margin:.6f}"
                    )
                target = plus_dipole if sign > 0 else minus_dipole
                target[mode] = snapshot.state_dipole_au
                if have_polar and snapshot.polarizability_au is not None:
                    polar_target = plus_polar if sign > 0 else minus_polar
                    polar_target[mode] = snapshot.polarizability_au
                elif have_polar:
                    # Raman remains explicitly unavailable; the independently
                    # valid target-state IR derivative is retained.
                    polar_complete = False
                overlap[mode, side_index] = snapshot.matched_overlap
                item = {
                    "mode": mode + 1,
                    "sign": sign,
                    "step_scale": step_scale,
                    "displacement_sqrt_amu_bohr": float(step[mode]),
                    "matched_overlap": snapshot.matched_overlap,
                    "tracking_margin": snapshot.tracking_margin,
                    "phase_to_central": snapshot.phase_to_central,
                    "raw_root_index": snapshot.raw_root_index,
                    "state_energy_hartree": snapshot.state_energy_hartree,
                    "state_irrep": snapshot.state_irrep,
                    "target_multiplicity": snapshot.target_multiplicity,
                    "expected_s2": snapshot.expected_s2,
                    "expected_s2_source": (
                        "fixed spin-adapted MRSF singlet/triplet sector; not an "
                        "independent numerical S^2 measurement"
                    ),
                    "response_representation": SPIN_ADAPTED_LAYOUT,
                    "response_block_weights": dict(snapshot.response_block_weights),
                    "state_dipole_au": snapshot.state_dipole_au.tolist(),
                    "state_dipole_unit": "e*bohr",
                    "polarizability_au": (
                        None
                        if snapshot.polarizability_au is None
                        else snapshot.polarizability_au.tolist()
                    ),
                    "polarizability_unit": "bohr^3",
                    "polarizability": dict(snapshot.polarizability_diagnostics),
                    "provenance": dict(snapshot.provenance),
                }
                mode_records.append(item)
                records.append(item)
            if {entry["sign"] for entry in mode_records} != {-1, +1}:
                raise RuntimeError(
                    "each normal coordinate and step requires plus and minus evaluations"
                )

        dipole_derivatives.append(
            finite_difference_excited_state_property(
                plus_dipole,
                minus_dipole,
                step,
                property_kind="state_dipole",
                coordinate_basis="excited_normal",
                coordinate_unit="sqrt(amu)*bohr",
                property_unit="e*bohr",
                electronic_state=request.electronic_state,
                electronic_phase_convention=phase_label,
                coordinate_phase_convention=request.coordinate_phase_convention,
                state_tracking_overlaps=overlap,
                minimum_overlap=request.minimum_state_overlap,
            )
        )
        polarizability_derivatives.append(
            finite_difference_excited_state_property(
                plus_polar,
                minus_polar,
                step,
                property_kind="polarizability",
                coordinate_basis="excited_normal",
                coordinate_unit="sqrt(amu)*bohr",
                property_unit="bohr^3",
                electronic_state=request.electronic_state,
                electronic_phase_convention=phase_label,
                coordinate_phase_convention=request.coordinate_phase_convention,
                state_tracking_overlaps=overlap,
                minimum_overlap=request.minimum_state_overlap,
            )
            if polar_complete
            else None
        )

    def convergence_record(
        name: str, derivatives: Sequence[ExcitedStatePropertyDerivative]
    ) -> dict[str, object]:
        if len(derivatives) != 3:
            raise RuntimeError("finite-difference validation requires h, h/2, and h/4")
        coarse, medium, fine = [item.values for item in derivatives]
        coarse_change = float(np.linalg.norm(medium - coarse))
        fine_change = float(np.linalg.norm(fine - medium))
        fine_norm = float(np.linalg.norm(fine))
        relative_change = fine_change / max(
            fine_norm, request.fd_absolute_tolerance
        )
        within_tolerance = (
            fine_change <= request.fd_absolute_tolerance
            or relative_change <= request.fd_relative_tolerance
        )
        decreasing = (
            coarse_change <= request.fd_absolute_tolerance
            or fine_change <= 1.25 * coarse_change
        )
        if not within_tolerance or not decreasing:
            raise ValueError(
                f"{name} finite difference is not converged over h, h/2, h/4: "
                f"fine change {fine_change:.3e}, relative change "
                f"{relative_change:.3e}, coarse change {coarse_change:.3e}"
            )
        return {
            "step_scales": list(step_scales),
            "coarse_to_medium_norm": coarse_change,
            "medium_to_fine_norm": fine_change,
            "medium_to_fine_relative_norm": relative_change,
            "relative_tolerance": request.fd_relative_tolerance,
            "absolute_tolerance": request.fd_absolute_tolerance,
            "accepted_step_scale": step_scales[-1],
        }

    dipole_convergence = convergence_record("MRSF state dipole", dipole_derivatives)
    dipole_derivative = dipole_derivatives[-1]
    infrared = excited_state_ir_intensities(dipole_derivative)

    polarizability_derivative = None
    polarizability_convergence = None
    polarizability_fd_failure = None
    raman = None
    if all(item is not None for item in polarizability_derivatives):
        complete_polar = [
            item for item in polarizability_derivatives if item is not None
        ]
        try:
            polarizability_convergence = convergence_record(
                "MRSF truncated-SOS polarizability", complete_polar
            )
        except ValueError as exc:
            # The permanent-dipole derivative has an independent electronic
            # source and convergence test, so an unconverged SOS Raman tensor
            # must not erase a valid target-state IR result.
            polarizability_fd_failure = str(exc)
        else:
            polarizability_derivative = complete_polar[-1]
            raman = excited_state_raman_activities(polarizability_derivative)

    provenance = {
        "schema_version": SCHEMA_VERSION,
        "electronic_method": electronic_method,
        "electronic_state": request.electronic_state,
        "state_index_one_based": request.state_index,
        "electronic_state_role": "target_excited_state",
        "response_representation": SPIN_ADAPTED_LAYOUT,
        "nuclear_derivative": (
            "two-point central finite differences at h, h/2, and h/4"
        ),
        "dipole_fd_convergence": dipole_convergence,
        "polarizability_fd_convergence": polarizability_convergence,
        "polarizability_fd_failure": polarizability_fd_failure,
        "dipole_source": (
            "nuclear moment plus complete spin-adapted MRSF state-density "
            "contraction at the molecular center of mass"
        ),
        "dipole_origin": "molecular center of mass",
        "polarizability_source": (
            "convergence-gated finite-state spin-adapted MRSF SOS; not analytic"
        ),
        "analytic_gap": (
            "analytic target-state MRSF dipole and polarizability derivatives "
            "and a uniform-field MRSF Hamiltonian are not implemented"
        ),
        "state_phase_use": (
            "response phases are tracked and recorded; permanent dipoles and "
            "static SOS polarizabilities are invariant to a unit state phase"
        ),
        "request_sha256": _sha256_json(request.to_dict()),
        **dict(source_provenance or {}),
    }
    return MRSFSpectroscopyFDResult(
        dipole_derivative=dipole_derivative,
        polarizability_derivative=polarizability_derivative,
        infrared=infrared,
        raman=raman,
        displacement_records=tuple(records),
        provenance=provenance,
    )


class OpenQPMRSFNormalModeEvaluator:
    """Evaluate displaced MRSF states with OpenQP's response-overlap tracker."""

    def __init__(self, central_mol: object, request: MRSFPropertyFDRequest):
        self.central_mol = central_mol
        self.request = request
        self.target_root = request.state_index - 1
        self.central_geometry = _finite_real(
            "central OpenQP geometry", getattr(central_mol, "get_system")()
        ).reshape((-1, 3))
        if self.central_geometry.size != request.normal_modes.shape[1]:
            raise ValueError("central OpenQP geometry and normal_modes disagree")
        config = getattr(central_mol, "config")
        self._validate_config(config, central_mol)
        nstate = int(config["tdhf"]["nstate"])
        if self.target_root >= nstate:
            raise ValueError(
                f"requested MRSF root {request.state_index} exceeds nstate={nstate}"
            )
        self.central_data = deepcopy(getattr(central_mol, "get_data")())
        self.central_symmetry_metadata = deepcopy(
            getattr(central_mol, "symmetry_metadata", None)
        )
        self.config = deepcopy(config)
        self._scratch = tempfile.TemporaryDirectory(prefix="openqp-mrsf-vibprop-")
        self._counter = 0

    @staticmethod
    def _validate_config(config: Mapping[str, Mapping[str, object]], mol: object) -> None:
        if str(config["input"]["method"]).lower() != "tdhf":
            raise ValueError("MRSF spectroscopy requires [input] method=tdhf")
        if str(config["tdhf"]["type"]).lower() != "mrsf":
            raise ValueError("MRSF spectroscopy requires [tdhf] type=mrsf")
        if str(config["scf"]["type"]).lower() != "rohf":
            raise ValueError("ordinary two-SOMO MRSF requires an ROHF reference")
        if int(config["scf"]["multiplicity"]) != 3:
            raise ValueError("ordinary two-SOMO MRSF requires a triplet ROHF reference")
        nalpha = int(np.asarray(getattr(mol, "data")["nelec_A"]).ravel()[0])
        nbeta = int(np.asarray(getattr(mol, "data")["nelec_B"]).ravel()[0])
        if nalpha - nbeta != 2:
            raise ValueError("ordinary MRSF spectroscopy requires exactly two SOMOs")
        if int(config["tdhf"]["multiplicity"]) not in (1, 3):
            raise ValueError("only isolated singlet and triplet MRSF targets are supported")

    def close(self) -> None:
        self._scratch.cleanup()

    def __enter__(self) -> "OpenQPMRSFNormalModeEvaluator":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    @property
    def electronic_method(self) -> str:
        functional = str(self.config["input"].get("functional", "")).strip()
        return "MRSF-TDDFT" if functional else "MRSF-TDHF"

    def _child_config(self) -> dict[str, dict[str, object]]:
        config = deepcopy(self.config)
        config["input"]["runtype"] = "energy"
        config["guess"]["type"] = "huckel"
        config["guess"]["save_mol"] = False
        config["properties"]["back_door"] = True
        config["properties"]["scf_prop"] = []
        config["nac"]["align"] = "reorder"
        config["symmetry"]["use_integral_symmetry"] = False
        config["symmetry"]["move_to_standard_frame"] = False

        # Runner passes a dictionary through ConfigParser, whose values must
        # have the same textual representation as an input file.  The central
        # molecule stores already-converted Python values, so serialize them
        # before constructing the displaced child calculation.
        def input_value(value: object) -> str:
            if isinstance(value, bool):
                return str(value).lower()
            if isinstance(value, (list, tuple, np.ndarray)):
                entries = value.tolist() if isinstance(value, np.ndarray) else value
                if len(entries) == 1 and not isinstance(
                    entries[0], (list, tuple, np.ndarray)
                ):
                    return str(entries[0])
                return ",".join(
                    " ".join(str(item) for item in entry)
                    if isinstance(entry, (list, tuple, np.ndarray))
                    else str(entry)
                    for entry in entries
                )
            return str(value)

        return {
            section: {key: input_value(value) for key, value in values.items()}
            for section, values in config.items()
        }

    def __call__(
        self, mode: int, sign: int, geometry_bohr: FloatArray
    ) -> MRSFTrackedPropertySnapshot:
        import oqp
        from oqp.analysis.transition_density import MRSFExcitedStates
        from oqp.library.single_point import BasisOverlap, NACME, SinglePoint
        from oqp.pyoqp import Runner

        self._counter += 1
        project = f"mrsf_vibprop_m{mode + 1}_{'p' if sign > 0 else 'm'}"
        log = os.path.join(self._scratch.name, f"{project}_{self._counter}.log")
        runner = Runner(
            project=project,
            input_dict=self._child_config(),
            log=log,
            silent=1,
            usempi=False,
        )
        mol = runner.mol
        mol.update_system(np.asarray(geometry_bohr, dtype=float).reshape((-1, 3)))
        if (getattr(mol, "symmetry_metadata", None) or {}).get(
            "status", "disabled"
        ) != "disabled":
            mol._detect_symmetry_metadata()
        mol.back_door = (self.central_geometry.copy(), deepcopy(self.central_data))

        driver = SinglePoint(mol)
        reference_energy = driver.reference()
        BasisOverlap(mol).overlap()
        driver.excitation(reference_energy)
        raw_response_vectors = np.array(mol.data["OQP::td_bvec_mo"], copy=True)
        raw_excitation_energies = np.array(
            mol.data["OQP::td_energies"], copy=True
        )
        raw_total_energies = np.array(mol.energies, copy=True)
        # BasisOverlap cannot infer geometry-derived symmetry metadata from an
        # in-memory back door, so provide the central record explicitly before
        # invoking the production isolated-root gate.
        mol._previous_symmetry_metadata = deepcopy(
            self.central_symmetry_metadata
        )
        tracking_report = NACME(mol).track_isolated_mrsf_hessian_root(
            self.request.state_index,
            overlap_tol=self.request.minimum_state_overlap,
            margin_tol=self.request.minimum_tracking_margin,
        )

        tracking = mol.get_state_tracking()
        if not tracking or not tracking.get("output_reordered"):
            raise RuntimeError("OpenQP did not produce a reordered state-tracking record")
        raw_root = int(tracking_report["selected_raw_state"]) - 1
        overlap = float(tracking_report["matched_overlap"])
        margin = float(tracking_report["assignment_margin"])
        phase = float(np.asarray(tracking["phase_step"])[self.target_root])

        # align_x reorders the response vectors for downstream NAC consumers,
        # while the state-interaction densities and energy ladder retain raw
        # Davidson order.  Restore the saved raw vector array before building
        # a self-consistent MRSFExcitedStates view.
        mol.data["OQP::td_bvec_mo"] = raw_response_vectors
        mol.data["OQP::td_energies"] = raw_excitation_energies
        mol.energies = raw_total_energies

        states = MRSFExcitedStates(mol)
        dipole = full_mrsf_state_dipole(states, raw_root, mol)
        blocks = response_block_weights(states, raw_root)
        state_irrep = tracking_report.get("selected_raw_irrep")
        target_multiplicity = int(self.config["tdhf"]["multiplicity"])
        expected_s2 = 0.0 if target_multiplicity == 1 else 2.0
        total_energies = _finite_real("MRSF total-state energies", mol.energies)
        if raw_root + 1 >= total_energies.size:
            raise RuntimeError("tracked MRSF root has no total-state energy")
        if self.request.polarizability_backend == "finite_field":
            raise NotImplementedError(
                "OpenQP has no uniform electric-field MRSF Hamiltonian; "
                "finite-field polarizability is unavailable"
            )
        try:
            polarizability, polar_diagnostics = truncated_sos_polarizability(
                states,
                raw_root,
                tail_states=self.request.sos_tail_states,
                tail_relative_tolerance=self.request.sos_tail_relative_tolerance,
                minimum_gap_hartree=self.request.sos_minimum_gap_hartree,
            )
        except ValueError as exc:
            polarizability = None
            polar_diagnostics = {
                "backend": "truncated_mrsf_state_sos",
                "status": "unavailable",
                "reason": str(exc),
                "analytic": False,
                "finite_field": False,
                "replacement_used": False,
            }
        return MRSFTrackedPropertySnapshot.create(
            state_dipole_au=dipole,
            polarizability_au=polarizability,
            matched_overlap=overlap,
            tracking_margin=margin,
            phase_to_central=phase,
            raw_root_index=raw_root,
            state_energy_hartree=float(total_energies[raw_root + 1]),
            target_multiplicity=target_multiplicity,
            expected_s2=expected_s2,
            state_irrep=state_irrep,
            response_block_weights=blocks,
            polarizability_diagnostics=polar_diagnostics,
            provenance={
                "engine": "OpenQP",
                "engine_version": str(getattr(oqp, "__version__", "unknown")),
                "child_project": project,
                "electronic_method": self.electronic_method,
                "state_tracking": "OpenQP cross-geometry MO alignment and response overlap",
                "root_mapping": "raw current root to central response root",
                "isolated_root_invariants": dict(tracking_report),
                "state_dipole": "nuclear plus complete MRSF state-density contraction",
                "response_representation": SPIN_ADAPTED_LAYOUT,
            },
        )


def run_openqp_mrsf_spectroscopy_fd(
    mol: object,
    modes: ArrayLike,
    *,
    state_index: int,
    electronic_state: str,
    coordinate_phase_convention: str,
    displacement: float = 1.0e-3,
    minimum_state_overlap: float = 0.99,
    minimum_tracking_margin: float = 0.05,
    fd_relative_tolerance: float = 0.05,
    fd_absolute_tolerance: float = 1.0e-6,
    polarizability_backend: PolarizabilityBackend = "truncated_sos",
    sos_tail_states: int = 2,
    sos_tail_relative_tolerance: float = 0.05,
    sos_minimum_gap_hartree: float = 1.0e-5,
) -> MRSFSpectroscopyFDResult:
    """Run the complete OpenQP normal-coordinate MRSF FD property path."""

    config = getattr(mol, "config")
    electronic_method = (
        "MRSF-TDDFT"
        if str(config["input"].get("functional", "")).strip()
        else "MRSF-TDHF"
    )
    request = MRSFPropertyFDRequest.create(
        electronic_method=electronic_method,
        electronic_state=electronic_state,
        state_index=state_index,
        normal_modes=modes,
        displacement=displacement,
        coordinate_phase_convention=coordinate_phase_convention,
        minimum_state_overlap=minimum_state_overlap,
        minimum_tracking_margin=minimum_tracking_margin,
        fd_relative_tolerance=fd_relative_tolerance,
        fd_absolute_tolerance=fd_absolute_tolerance,
        polarizability_backend=polarizability_backend,
        sos_tail_states=sos_tail_states,
        sos_tail_relative_tolerance=sos_tail_relative_tolerance,
        sos_minimum_gap_hartree=sos_minimum_gap_hartree,
    )
    with OpenQPMRSFNormalModeEvaluator(mol, request) as evaluator:
        return assemble_mrsf_spectroscopy_derivatives(
            request,
            evaluator,
            evaluator.central_geometry,
            electronic_method=evaluator.electronic_method,
            source_provenance={
                "openqp_input": str(getattr(mol, "input_file", "")),
                "openqp_project": str(getattr(mol, "project_name", "")),
            },
        )


__all__ = (
    "MRSFPropertyFDRequest",
    "MRSFSpectroscopyFDResult",
    "MRSFTrackedPropertySnapshot",
    "OpenQPMRSFNormalModeEvaluator",
    "SCHEMA_VERSION",
    "SPIN_ADAPTED_LAYOUT",
    "assemble_mrsf_spectroscopy_derivatives",
    "full_mrsf_state_dipole",
    "response_block_weights",
    "run_openqp_mrsf_spectroscopy_fd",
    "truncated_sos_polarizability",
)
