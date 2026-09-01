"""Harmonic vibronic overlaps and property-derivative observables.

The Franck--Condon implementation evaluates the multidimensional harmonic-
oscillator generating function for

``Q_excited = J @ Q_ground + K``.

It therefore retains frequency changes, Duschinsky rotation, and displacement
without a product-mode approximation.  The public functions require positive
frequencies and a square orthogonal ``J`` because a noncanonical coordinate
transformation does not define normalized harmonic wavefunctions.

Herzberg--Teller, excited-state infrared, and nonresonant Raman quantities are
evaluated only from explicitly supplied electronic-property derivatives.  This
module never substitutes a reference-state ROHF/SCF property for an excited-
state MRSF property.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import argparse
import json
from collections.abc import Mapping
from math import exp, lgamma, log, pi, sqrt
from pathlib import Path
from typing import Literal, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]

CM1_PER_HARTREE = 219474.6313705
CM1_TO_HARTREE = 1.0 / CM1_PER_HARTREE
AMU_TO_ELECTRON_MASS = 1822.888486209
BOLTZMANN_CM1_PER_K = 0.6950348004861274
# Integrated harmonic IR intensity for a dipole derivative expressed as
# e*bohr / (sqrt(amu)*bohr) = e / sqrt(amu).  The often quoted 42.255 factor
# instead applies to derivatives in debye / (angstrom*sqrt(amu)).
IR_INTENSITY_CONVERSION_KM_MOL = 974.88011

CoordinateUnit = Literal["sqrt(amu)*bohr", "sqrt(me)*bohr"]
OriginKind = Literal["zero_zero", "adiabatic_minima"]
BroadeningKind = Literal["gaussian", "lorentzian"]
PropertyKind = Literal[
    "transition_dipole", "state_dipole", "polarizability"
]


def _finite_real(name: str, values: ArrayLike) -> FloatArray:
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _finite_numeric(name: str, values: ArrayLike) -> ComplexArray:
    array = np.asarray(values, dtype=complex)
    if not np.all(np.isfinite(array.real)) or not np.all(np.isfinite(array.imag)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _validate_coordinate_unit(unit: str) -> CoordinateUnit:
    if unit not in ("sqrt(amu)*bohr", "sqrt(me)*bohr"):
        raise ValueError(
            "coordinate_unit must be 'sqrt(amu)*bohr' or 'sqrt(me)*bohr'"
        )
    return unit  # type: ignore[return-value]


def _coordinate_scale_to_electron_mass(unit: CoordinateUnit) -> float:
    return sqrt(AMU_TO_ELECTRON_MASS) if unit == "sqrt(amu)*bohr" else 1.0


def _trapezoid(values: FloatArray, grid: FloatArray) -> float:
    integration = getattr(np, "trapezoid", None)
    if integration is None:  # NumPy < 2.0
        integration = np.trapz
    return float(integration(values, grid))


def _validate_excited_state_role(role: str) -> None:
    if role != "target_excited_state":
        raise ValueError(
            "electronic_state_role must be 'target_excited_state'; reference, "
            "ROHF, SCF, and unspecified properties cannot represent an "
            "excited-state observable"
        )


@dataclass(frozen=True)
class HarmonicVibronicModel:
    """Two harmonic surfaces in a common mass-weighted coordinate system.

    Frequencies are positive wavenumbers in cm^-1.  ``K`` uses the explicitly
    stated mass-weighted coordinate unit.  ``coordinate_phase_convention``
    identifies the normal-mode phases and must match all supplied normal-
    coordinate property derivatives.
    """

    ground_frequencies_cm1: FloatArray
    excited_frequencies_cm1: FloatArray
    J: FloatArray
    K: FloatArray
    coordinate_unit: CoordinateUnit
    coordinate_phase_convention: str
    orthogonality_tolerance: float
    orthogonality_residual: float

    @classmethod
    def create(
        cls,
        ground_frequencies_cm1: ArrayLike,
        excited_frequencies_cm1: ArrayLike,
        J: ArrayLike,
        K: ArrayLike,
        *,
        coordinate_unit: CoordinateUnit,
        coordinate_phase_convention: str,
        orthogonality_tolerance: float = 1.0e-8,
    ) -> "HarmonicVibronicModel":
        ground = _finite_real("ground_frequencies_cm1", ground_frequencies_cm1)
        excited = _finite_real(
            "excited_frequencies_cm1", excited_frequencies_cm1
        )
        if ground.ndim != 1 or ground.size == 0:
            raise ValueError("ground_frequencies_cm1 must be a nonempty vector")
        if excited.shape != ground.shape:
            raise ValueError(
                "ground and excited frequencies must have the same shape"
            )
        if np.any(ground <= 0.0) or np.any(excited <= 0.0):
            raise ValueError(
                "harmonic vibronic spectra require strictly positive frequencies"
            )
        rotation = _finite_real("J", J)
        displacement = _finite_real("K", K)
        nmode = ground.size
        if rotation.shape != (nmode, nmode):
            raise ValueError(f"J must have shape ({nmode}, {nmode})")
        if displacement.shape != (nmode,):
            raise ValueError(f"K must have shape ({nmode},)")
        if orthogonality_tolerance < 0.0:
            raise ValueError("orthogonality_tolerance must be non-negative")
        left = rotation @ rotation.T - np.eye(nmode)
        right = rotation.T @ rotation - np.eye(nmode)
        residual = float(
            max(
                np.max(np.abs(left), initial=0.0),
                np.max(np.abs(right), initial=0.0),
            )
        )
        if residual > orthogonality_tolerance:
            raise ValueError(
                f"Duschinsky orthogonality residual {residual:.3e} exceeds "
                f"{orthogonality_tolerance:.3e}; the harmonic overlap model "
                "requires a canonical vibrational-coordinate transformation"
            )
        if not coordinate_phase_convention.strip():
            raise ValueError(
                "coordinate_phase_convention must identify the normal-mode phases"
            )
        unit = _validate_coordinate_unit(coordinate_unit)
        return cls(
            ground_frequencies_cm1=ground.copy(),
            excited_frequencies_cm1=excited.copy(),
            J=rotation.copy(),
            K=displacement.copy(),
            coordinate_unit=unit,
            coordinate_phase_convention=coordinate_phase_convention,
            orthogonality_tolerance=orthogonality_tolerance,
            orthogonality_residual=residual,
        )

    @classmethod
    def from_duschinsky(
        cls,
        result: object,
        ground_frequencies_cm1: ArrayLike,
        excited_frequencies_cm1: ArrayLike,
        *,
        coordinate_unit: CoordinateUnit = "sqrt(amu)*bohr",
        coordinate_phase_convention: str = (
            "OpenQP Duschinsky aligned normal modes"
        ),
        orthogonality_tolerance: float = 1.0e-8,
    ) -> "HarmonicVibronicModel":
        """Construct the spectral model from a ``DuschinskyResult``-like object."""

        if not hasattr(result, "J") or not hasattr(result, "K"):
            raise TypeError("result must provide J and K arrays")
        return cls.create(
            ground_frequencies_cm1,
            excited_frequencies_cm1,
            getattr(result, "J"),
            getattr(result, "K"),
            coordinate_unit=coordinate_unit,
            coordinate_phase_convention=coordinate_phase_convention,
            orthogonality_tolerance=orthogonality_tolerance,
        )

    @property
    def nmode(self) -> int:
        return int(self.ground_frequencies_cm1.size)


@dataclass(frozen=True)
class ElectronicTransitionMoment:
    """Phase-labelled equilibrium electronic transition dipole."""

    value: ComplexArray
    electronic_state: str
    electronic_phase_convention: str
    dipole_unit: str = "e*bohr"
    electronic_state_role: str = "target_excited_state"

    @classmethod
    def create(
        cls,
        value: ArrayLike,
        *,
        electronic_state: str,
        electronic_phase_convention: str,
        dipole_unit: str = "e*bohr",
        electronic_state_role: str = "target_excited_state",
    ) -> "ElectronicTransitionMoment":
        _validate_excited_state_role(electronic_state_role)
        dipole = _finite_numeric("transition_dipole", value)
        if dipole.shape != (3,):
            raise ValueError("transition_dipole must have shape (3,)")
        if not electronic_state.strip() or not electronic_phase_convention.strip():
            raise ValueError(
                "electronic_state and electronic_phase_convention must be nonempty"
            )
        if dipole_unit != "e*bohr":
            raise ValueError("only transition dipoles in e*bohr are supported")
        return cls(
            value=dipole.copy(),
            electronic_state=electronic_state,
            electronic_phase_convention=electronic_phase_convention,
            dipole_unit=dipole_unit,
            electronic_state_role=electronic_state_role,
        )


@dataclass(frozen=True)
class ExcitedStatePropertyDerivative:
    """Explicit property derivative for a tracked target excited state."""

    values: ComplexArray
    property_kind: PropertyKind
    coordinate_basis: str
    coordinate_unit: str
    property_unit: str
    electronic_state: str
    electronic_state_role: str
    electronic_phase_convention: str
    coordinate_phase_convention: str
    provenance: str
    minimum_state_overlap: float | None = None

    @classmethod
    def create(
        cls,
        values: ArrayLike,
        *,
        property_kind: PropertyKind,
        coordinate_basis: str,
        coordinate_unit: str,
        property_unit: str,
        electronic_state: str,
        electronic_phase_convention: str,
        coordinate_phase_convention: str,
        provenance: str,
        electronic_state_role: str,
        minimum_state_overlap: float | None = None,
    ) -> "ExcitedStatePropertyDerivative":
        _validate_excited_state_role(electronic_state_role)
        if property_kind not in (
            "transition_dipole",
            "state_dipole",
            "polarizability",
        ):
            raise ValueError("unsupported property_kind")
        array = _finite_numeric("property derivative", values)
        if property_kind in ("transition_dipole", "state_dipole"):
            if array.ndim != 2 or array.shape[0] != 3:
                raise ValueError(
                    f"{property_kind} derivative must have shape (3, ncoordinate)"
                )
        elif array.ndim != 3 or array.shape[:2] != (3, 3):
            raise ValueError(
                "polarizability derivative must have shape (3, 3, ncoordinate)"
            )
        required = {
            "coordinate_basis": coordinate_basis,
            "coordinate_unit": coordinate_unit,
            "property_unit": property_unit,
            "electronic_state": electronic_state,
            "electronic_phase_convention": electronic_phase_convention,
            "coordinate_phase_convention": coordinate_phase_convention,
            "provenance": provenance,
        }
        empty = [name for name, value in required.items() if not value.strip()]
        if empty:
            raise ValueError(f"nonempty metadata required for: {', '.join(empty)}")
        if minimum_state_overlap is not None and not (
            0.0 <= minimum_state_overlap <= 1.0
        ):
            raise ValueError("minimum_state_overlap must lie in [0, 1]")
        return cls(
            values=array.copy(),
            property_kind=property_kind,
            coordinate_basis=coordinate_basis,
            coordinate_unit=coordinate_unit,
            property_unit=property_unit,
            electronic_state=electronic_state,
            electronic_state_role=electronic_state_role,
            electronic_phase_convention=electronic_phase_convention,
            coordinate_phase_convention=coordinate_phase_convention,
            provenance=provenance,
            minimum_state_overlap=minimum_state_overlap,
        )

    def to_dict(self) -> dict[str, object]:
        """Return a unit- and provenance-complete JSON-compatible record."""

        return {
            "values": _complex_json(self.values),
            "property_kind": self.property_kind,
            "coordinate_basis": self.coordinate_basis,
            "coordinate_unit": self.coordinate_unit,
            "property_unit": self.property_unit,
            "electronic_state": self.electronic_state,
            "electronic_state_role": self.electronic_state_role,
            "electronic_phase_convention": self.electronic_phase_convention,
            "coordinate_phase_convention": self.coordinate_phase_convention,
            "provenance": self.provenance,
            "minimum_state_overlap": self.minimum_state_overlap,
        }


def mrsf_analytic_property_derivative_from_provider(
    provider: object,
    *,
    property_kind: PropertyKind,
    electronic_state: str,
) -> ExcitedStatePropertyDerivative:
    """Read an analytic first derivative through the MRSF-specific hook.

    The provider must implement ``get_mrsf_analytic_property_derivative`` and
    return a mapping containing every metadata field consumed below.  Generic
    reference-state dipole or polarizability accessors are intentionally never
    inspected as a fallback.
    """

    getter = getattr(provider, "get_mrsf_analytic_property_derivative", None)
    if not callable(getter):
        raise TypeError(
            "provider must implement get_mrsf_analytic_property_derivative(); "
            "generic ROHF/SCF property accessors are not accepted"
        )
    record = getter(
        property_kind=property_kind,
        electronic_state=electronic_state,
    )
    if not isinstance(record, Mapping):
        raise TypeError("MRSF analytic derivative provider must return a mapping")
    required = {
        "values",
        "property_kind",
        "coordinate_basis",
        "coordinate_unit",
        "property_unit",
        "electronic_state",
        "electronic_state_role",
        "electronic_phase_convention",
        "coordinate_phase_convention",
        "electronic_method",
        "derivative_order",
        "provenance",
    }
    missing = sorted(required - set(record))
    if missing:
        raise ValueError(
            "MRSF analytic derivative record is missing: " + ", ".join(missing)
        )
    if record["property_kind"] != property_kind:
        raise ValueError("MRSF derivative provider returned a different property")
    if record["electronic_state"] != electronic_state:
        raise ValueError("MRSF derivative provider returned a different state")
    if record["electronic_method"] not in ("MRSF-TDDFT", "MRSF-TDHF"):
        raise ValueError("electronic_method must be MRSF-TDDFT or MRSF-TDHF")
    if record["derivative_order"] != 1:
        raise ValueError("MRSF property derivative_order must be exactly 1")
    provenance = str(record["provenance"])
    if "analytic" not in provenance.lower():
        raise ValueError("MRSF analytic derivative provenance must state 'analytic'")
    return ExcitedStatePropertyDerivative.create(
        record["values"],
        property_kind=property_kind,
        coordinate_basis=str(record["coordinate_basis"]),
        coordinate_unit=str(record["coordinate_unit"]),
        property_unit=str(record["property_unit"]),
        electronic_state=electronic_state,
        electronic_state_role=str(record["electronic_state_role"]),
        electronic_phase_convention=str(record["electronic_phase_convention"]),
        coordinate_phase_convention=str(record["coordinate_phase_convention"]),
        provenance=provenance,
    )


def finite_difference_excited_state_property(
    plus_values: ArrayLike,
    minus_values: ArrayLike,
    displacements: ArrayLike,
    *,
    property_kind: PropertyKind,
    coordinate_basis: str,
    coordinate_unit: str,
    property_unit: str,
    electronic_state: str,
    electronic_phase_convention: str,
    coordinate_phase_convention: str,
    state_tracking_overlaps: ArrayLike,
    minimum_overlap: float = 0.9,
    transition_phase_factors: ArrayLike | None = None,
) -> ExcitedStatePropertyDerivative:
    """Central difference of explicitly evaluated, state-tracked properties.

    Input arrays have shape ``(ncoordinate, *property_shape)``.  Transition-
    dipole differences additionally require complex unit-modulus phase factors
    with shape ``(ncoordinate, 2)`` for the plus and minus geometries.  The
    factors multiply the raw transition moments and must align the complete
    bra--ket phase with the central geometry.  The returned derivative never
    invokes an electronic-structure calculation.  The overlap threshold is an
    isolated-state diagnostic, not a replacement for energy, symmetry, spin,
    and response-vector identity checks by the electronic-structure driver.
    """

    plus = _finite_numeric("plus_values", plus_values)
    minus = _finite_numeric("minus_values", minus_values)
    if plus.shape != minus.shape or plus.ndim < 2:
        raise ValueError("plus_values and minus_values must have the same shape")
    step = _finite_real("displacements", displacements)
    if step.shape != (plus.shape[0],) or np.any(step <= 0.0):
        raise ValueError(
            "displacements must be a positive vector with one value per coordinate"
        )
    overlaps = _finite_real("state_tracking_overlaps", state_tracking_overlaps)
    if overlaps.shape != (plus.shape[0], 2):
        raise ValueError(
            "state_tracking_overlaps must have shape (ncoordinate, 2)"
        )
    if not 0.0 <= minimum_overlap <= 1.0:
        raise ValueError("minimum_overlap must lie in [0, 1]")
    if np.any(overlaps < 0.0) or np.any(overlaps > 1.0):
        raise ValueError("state_tracking_overlaps must lie in [0, 1]")
    observed_minimum = float(np.min(overlaps))
    if observed_minimum < minimum_overlap:
        raise ValueError(
            f"state-tracking overlap {observed_minimum:.6f} is below "
            f"the required {minimum_overlap:.6f}"
        )

    if property_kind == "transition_dipole":
        if transition_phase_factors is None:
            raise ValueError(
                "transition_phase_factors are required for a transition-dipole "
                "finite difference"
            )
        phases = _finite_numeric(
            "transition_phase_factors", transition_phase_factors
        )
        if phases.shape != (plus.shape[0], 2):
            raise ValueError(
                "transition_phase_factors must have shape (ncoordinate, 2)"
            )
        if np.max(np.abs(np.abs(phases) - 1.0), initial=0.0) > 1.0e-10:
            raise ValueError("transition_phase_factors must have unit modulus")
        phase_shape = (plus.shape[0],) + (1,) * (plus.ndim - 1)
        plus = plus * phases[:, 0].reshape(phase_shape)
        minus = minus * phases[:, 1].reshape(phase_shape)
    elif transition_phase_factors is not None:
        raise ValueError(
            "transition_phase_factors apply only to transition-dipole derivatives"
        )

    derivative_coordinate_first = (plus - minus) / (
        2.0 * step.reshape((step.size,) + (1,) * (plus.ndim - 1))
    )
    derivative = np.moveaxis(derivative_coordinate_first, 0, -1)
    return ExcitedStatePropertyDerivative.create(
        derivative,
        property_kind=property_kind,
        coordinate_basis=coordinate_basis,
        coordinate_unit=coordinate_unit,
        property_unit=property_unit,
        electronic_state=electronic_state,
        electronic_phase_convention=electronic_phase_convention,
        coordinate_phase_convention=coordinate_phase_convention,
        provenance="central finite difference of state-tracked evaluations",
        electronic_state_role="target_excited_state",
        minimum_state_overlap=observed_minimum,
    )


class HarmonicOverlapEngine:
    """Analytic multidimensional harmonic-oscillator overlap evaluator."""

    def __init__(self, model: HarmonicVibronicModel):
        if not isinstance(model, HarmonicVibronicModel):
            raise TypeError("model must be a HarmonicVibronicModel")
        validated = HarmonicVibronicModel.create(
            model.ground_frequencies_cm1,
            model.excited_frequencies_cm1,
            model.J,
            model.K,
            coordinate_unit=model.coordinate_unit,
            coordinate_phase_convention=model.coordinate_phase_convention,
            orthogonality_tolerance=model.orthogonality_tolerance,
        )
        self.model = validated
        scale = _coordinate_scale_to_electron_mass(validated.coordinate_unit)
        self._K = validated.K * scale
        ground = validated.ground_frequencies_cm1 * CM1_TO_HARTREE
        excited = validated.excited_frequencies_cm1 * CM1_TO_HARTREE
        root_ground = np.sqrt(ground)
        root_excited = np.sqrt(excited)
        omega_ground = np.diag(ground)
        omega_excited = np.diag(excited)
        A = omega_ground + validated.J.T @ omega_excited @ validated.J
        b0 = -validated.J.T @ omega_excited @ self._K
        B = np.concatenate(
            (
                2.0 * np.diag(root_ground),
                2.0 * validated.J.T @ np.diag(root_excited),
            ),
            axis=1,
        )
        linear_external = np.concatenate(
            (np.zeros(validated.nmode), 2.0 * root_excited * self._K)
        )
        solved_b0 = np.linalg.solve(A, b0)
        solved_B = np.linalg.solve(A, B)
        self._d = linear_external + B.T @ solved_b0
        self._R = B.T @ solved_B - 2.0 * np.eye(2 * validated.nmode)
        determinant_sign, log_determinant = np.linalg.slogdet(A)
        if determinant_sign <= 0.0:
            raise ValueError("combined harmonic Gaussian is not positive definite")
        log_s00 = (
            0.5 * validated.nmode * log(2.0)
            + 0.25 * np.sum(np.log(ground * excited))
            - 0.5 * log_determinant
            - 0.5 * self._K @ omega_excited @ self._K
            + 0.5 * b0 @ solved_b0
        )
        self.zero_zero_overlap = float(exp(log_s00))

    @lru_cache(maxsize=None)
    def _derivative_polynomial(self, alpha: tuple[int, ...]) -> float:
        if not any(alpha):
            return 1.0
        index = next(i for i, value in enumerate(alpha) if value)
        beta = list(alpha)
        beta[index] -= 1
        beta_tuple = tuple(beta)
        value = self._d[index] * self._derivative_polynomial(beta_tuple)
        for other, count in enumerate(beta):
            if count:
                reduced = beta.copy()
                reduced[other] -= 1
                value += (
                    self._R[index, other]
                    * count
                    * self._derivative_polynomial(tuple(reduced))
                )
        return float(value)

    def overlap(self, initial: Sequence[int], final: Sequence[int]) -> float:
        """Return ``<final,excited|initial,ground>``."""

        initial_state = _validate_state("initial", initial, self.model.nmode)
        final_state = _validate_state("final", final, self.model.nmode)
        alpha = initial_state + final_state
        log_norm = 0.5 * (
            sum(alpha) * log(2.0) + sum(lgamma(value + 1.0) for value in alpha)
        )
        return self.zero_zero_overlap * self._derivative_polynomial(alpha) * exp(
            -log_norm
        )

    def ground_coordinate_overlap(
        self,
        initial: Sequence[int],
        final: Sequence[int],
        mode: int,
    ) -> float:
        """Return ``<final,excited|Q_ground_mode|initial,ground>``."""

        initial_state = list(_validate_state("initial", initial, self.model.nmode))
        if not 0 <= mode < self.model.nmode:
            raise IndexError("mode index is outside the vibrational space")
        omega = self.model.ground_frequencies_cm1[mode] * CM1_TO_HARTREE
        result = 0.0
        if initial_state[mode] > 0:
            lowered = initial_state.copy()
            lowered[mode] -= 1
            result += sqrt(initial_state[mode] / (2.0 * omega)) * self.overlap(
                lowered, final
            )
        raised = initial_state.copy()
        raised[mode] += 1
        result += sqrt((initial_state[mode] + 1.0) / (2.0 * omega)) * self.overlap(
            raised, final
        )
        return result


def _validate_state(name: str, state: Sequence[int], nmode: int) -> tuple[int, ...]:
    if len(state) != nmode:
        raise ValueError(f"{name} state must contain {nmode} occupations")
    values = tuple(int(value) for value in state)
    if any(value < 0 for value in values) or any(
        float(source) != value for source, value in zip(state, values)
    ):
        raise ValueError(f"{name} state occupations must be non-negative integers")
    return values


def enumerate_vibrational_states(
    nmode: int, max_total_quanta: int, *, max_states: int = 250000
) -> tuple[tuple[int, ...], ...]:
    """Enumerate product states with a bounded total quantum number."""

    if nmode <= 0 or max_total_quanta < 0 or max_states <= 0:
        raise ValueError("state-enumeration bounds must be non-negative")
    count = 1
    for numerator in range(nmode + 1, nmode + max_total_quanta + 1):
        count = count * numerator // (numerator - nmode)
    if count > max_states:
        raise ValueError(
            f"requested vibrational basis contains {count} states, exceeding "
            f"max_states={max_states}"
        )
    states: list[tuple[int, ...]] = []

    def append_compositions(
        remaining: int, dimensions: int, prefix: tuple[int, ...]
    ) -> None:
        if dimensions == 1:
            states.append(prefix + (remaining,))
            return
        for occupation in range(remaining + 1):
            append_compositions(
                remaining - occupation, dimensions - 1, prefix + (occupation,)
            )

    for total_quanta in range(max_total_quanta + 1):
        append_compositions(total_quanta, nmode, ())
    return tuple(states)


def _thermal_populations(
    frequencies_cm1: FloatArray,
    states: Sequence[tuple[int, ...]],
    temperature_kelvin: float,
) -> tuple[FloatArray, float]:
    if temperature_kelvin < 0.0:
        raise ValueError("temperature_kelvin must be non-negative")
    if temperature_kelvin == 0.0:
        populations = np.array(
            [1.0 if not any(state) else 0.0 for state in states], dtype=float
        )
        return populations, float(np.sum(populations))
    beta = 1.0 / (BOLTZMANN_CM1_PER_K * temperature_kelvin)
    log_partition_inverse = float(np.sum(np.log1p(-np.exp(-beta * frequencies_cm1))))
    populations = np.array(
        [
            exp(log_partition_inverse - beta * np.dot(state, frequencies_cm1))
            for state in states
        ],
        dtype=float,
    )
    return populations, float(np.sum(populations))


def _validate_transition_inputs(
    model: HarmonicVibronicModel,
    transition: ElectronicTransitionMoment | None,
    derivative: ExcitedStatePropertyDerivative | None,
) -> tuple[ComplexArray | None, ComplexArray | None]:
    if derivative is not None and transition is None:
        raise ValueError(
            "an explicit equilibrium transition dipole, including an explicit "
            "zero vector for a forbidden transition, is required with a "
            "Herzberg--Teller derivative"
        )
    if transition is None:
        return None, None
    _validate_excited_state_role(transition.electronic_state_role)
    transition_value = _finite_numeric(
        "transition_dipole", transition.value
    )
    if transition_value.shape != (3,) or transition.dipole_unit != "e*bohr":
        raise ValueError("transition dipole must have shape (3,) and units e*bohr")
    if derivative is None:
        return transition_value, None
    _validate_excited_state_role(derivative.electronic_state_role)
    if derivative.property_kind != "transition_dipole":
        raise ValueError("Herzberg--Teller input must be a transition-dipole derivative")
    if derivative.coordinate_basis != "ground_normal":
        raise ValueError(
            "Herzberg--Teller derivatives must use the ground_normal coordinate basis"
        )
    if derivative.electronic_state != transition.electronic_state:
        raise ValueError("transition dipole and derivative refer to different states")
    if (
        derivative.electronic_phase_convention
        != transition.electronic_phase_convention
    ):
        raise ValueError(
            "transition dipole and derivative have inconsistent electronic phases"
        )
    if derivative.coordinate_phase_convention != model.coordinate_phase_convention:
        raise ValueError(
            "transition-dipole derivative and harmonic model have inconsistent "
            "normal-coordinate phases"
        )
    derivative_values = _finite_numeric(
        "transition_dipole_derivative", derivative.values
    )
    if derivative_values.shape != (3, model.nmode):
        raise ValueError(
            f"transition-dipole derivative must have shape (3, {model.nmode})"
        )
    if (
        derivative.property_unit != "e*bohr"
        or derivative.coordinate_unit != model.coordinate_unit
    ):
        raise ValueError(
            "transition-dipole derivative must use e*bohr per the model coordinate unit"
        )
    scale = _coordinate_scale_to_electron_mass(model.coordinate_unit)
    return transition_value, derivative_values / scale


def vibronic_transition_moment(
    engine: HarmonicOverlapEngine,
    initial: Sequence[int],
    final: Sequence[int],
    transition: ElectronicTransitionMoment,
    derivative: ExcitedStatePropertyDerivative | None = None,
) -> ComplexArray:
    """Return the Condon plus first-order Herzberg--Teller transition moment."""

    dipole, gradient = _validate_transition_inputs(
        engine.model, transition, derivative
    )
    assert dipole is not None
    value = dipole * engine.overlap(initial, final)
    if gradient is not None:
        for mode in range(engine.model.nmode):
            value = value + gradient[:, mode] * engine.ground_coordinate_overlap(
                initial, final, mode
            )
    return np.asarray(value, dtype=complex)


@dataclass(frozen=True)
class VibronicLine:
    initial_state: tuple[int, ...]
    final_state: tuple[int, ...]
    position_cm1: float
    overlap: float
    franck_condon_factor: float
    thermal_population: float
    transition_moment: ComplexArray | None
    raw_strength: float
    normalized_strength: float


@dataclass(frozen=True)
class VibronicSpectrum:
    lines: tuple[VibronicLine, ...]
    wavenumbers_cm1: FloatArray
    intensity: FloatArray
    temperature_kelvin: float
    broadening_kind: BroadeningKind
    fwhm_cm1: float
    normalization: str
    retained_thermal_population: float
    franck_condon_completeness: float


def harmonic_vibronic_spectrum(
    model: HarmonicVibronicModel,
    *,
    electronic_origin_cm1: float,
    origin_kind: OriginKind,
    max_final_quanta: int,
    temperature_kelvin: float = 0.0,
    max_initial_quanta: int | None = None,
    transition: ElectronicTransitionMoment | None = None,
    transition_dipole_derivative: ExcitedStatePropertyDerivative | None = None,
    broadening_kind: BroadeningKind = "gaussian",
    fwhm_cm1: float = 100.0,
    normalization: Literal["none", "sum"] = "sum",
    grid_cm1: ArrayLike | None = None,
    minimum_thermal_population: float = 0.999,
    minimum_franck_condon_completeness: float = 0.999,
    max_states: int = 250000,
) -> VibronicSpectrum:
    """Calculate harmonic FC/FC--HT stick and broadened absorption spectra.

    Strengths are relative vibronic strengths, not absorption cross sections.
    With no electronic transition dipole, strengths are thermally weighted
    Franck--Condon factors.  With electronic data, they are squared transition
    moments in ``(e*bohr)^2`` before optional normalization.
    """

    engine = HarmonicOverlapEngine(model)
    model = engine.model
    if not np.isfinite(electronic_origin_cm1):
        raise ValueError("electronic_origin_cm1 must be finite")
    if origin_kind not in ("zero_zero", "adiabatic_minima"):
        raise ValueError("origin_kind must be 'zero_zero' or 'adiabatic_minima'")
    if broadening_kind not in ("gaussian", "lorentzian"):
        raise ValueError("unsupported broadening_kind")
    if fwhm_cm1 <= 0.0 or not np.isfinite(fwhm_cm1):
        raise ValueError("fwhm_cm1 must be positive and finite")
    if normalization not in ("none", "sum"):
        raise ValueError("normalization must be 'none' or 'sum'")
    if max_initial_quanta is None:
        if temperature_kelvin > 0.0:
            raise ValueError(
                "max_initial_quanta is required at nonzero temperature so that "
                "thermal truncation is explicit"
            )
        max_initial_quanta = 0
    initial_states = enumerate_vibrational_states(
        model.nmode, max_initial_quanta, max_states=max_states
    )
    final_states = enumerate_vibrational_states(
        model.nmode, max_final_quanta, max_states=max_states
    )
    populations, retained_population = _thermal_populations(
        model.ground_frequencies_cm1, initial_states, temperature_kelvin
    )
    if not 0.0 <= minimum_thermal_population <= 1.0:
        raise ValueError("minimum_thermal_population must lie in [0, 1]")
    if retained_population + 1.0e-15 < minimum_thermal_population:
        raise ValueError(
            f"retained thermal population {retained_population:.8f} is below "
            f"the required {minimum_thermal_population:.8f}; increase "
            "max_initial_quanta"
        )
    dipole, derivative_values = _validate_transition_inputs(
        model, transition, transition_dipole_derivative
    )
    provisional: list[
        tuple[
            tuple[int, ...],
            tuple[int, ...],
            float,
            float,
            float,
            ComplexArray | None,
            float,
        ]
    ] = []
    fc_weight = 0.0
    for initial, population in zip(initial_states, populations):
        if population == 0.0:
            continue
        for final in final_states:
            overlap = engine.overlap(initial, final)
            fc_factor = overlap * overlap
            fc_weight += population * fc_factor
            if origin_kind == "zero_zero":
                position = electronic_origin_cm1 + float(
                    np.dot(final, model.excited_frequencies_cm1)
                    - np.dot(initial, model.ground_frequencies_cm1)
                )
            else:
                position = electronic_origin_cm1 + float(
                    np.dot(
                        np.asarray(final) + 0.5,
                        model.excited_frequencies_cm1,
                    )
                    - np.dot(
                        np.asarray(initial) + 0.5,
                        model.ground_frequencies_cm1,
                    )
                )
            if dipole is None:
                moment = None
                raw_strength = population * fc_factor
            else:
                assert transition is not None
                if derivative_values is None:
                    moment = dipole * overlap
                else:
                    moment = vibronic_transition_moment(
                        engine,
                        initial,
                        final,
                        transition,
                        transition_dipole_derivative,
                    )
                raw_strength = population * float(np.vdot(moment, moment).real)
            provisional.append(
                (initial, final, position, overlap, population, moment, raw_strength)
            )
    completeness = fc_weight / retained_population if retained_population else 0.0
    if not 0.0 <= minimum_franck_condon_completeness <= 1.0:
        raise ValueError(
            "minimum_franck_condon_completeness must lie in [0, 1]"
        )
    if completeness + 1.0e-15 < minimum_franck_condon_completeness:
        raise ValueError(
            f"Franck-Condon completeness {completeness:.8f} is below the "
            f"required {minimum_franck_condon_completeness:.8f}; increase "
            "max_final_quanta"
        )
    total_strength = sum(item[-1] for item in provisional)
    if normalization == "sum" and total_strength <= 0.0:
        raise ValueError("spectrum cannot be normalized because total strength is zero")
    strength_scale = 1.0 / total_strength if normalization == "sum" else 1.0
    lines = tuple(
        VibronicLine(
            initial_state=initial,
            final_state=final,
            position_cm1=position,
            overlap=overlap,
            franck_condon_factor=overlap * overlap,
            thermal_population=population,
            transition_moment=moment,
            raw_strength=raw_strength,
            normalized_strength=raw_strength * strength_scale,
        )
        for initial, final, position, overlap, population, moment, raw_strength in provisional
    )
    if grid_cm1 is None:
        positions = np.array([line.position_cm1 for line in lines])
        extent = 12.0 * fwhm_cm1 if broadening_kind == "lorentzian" else 5.0 * fwhm_cm1
        grid = np.linspace(
            float(np.min(positions) - extent),
            float(np.max(positions) + extent),
            max(1001, int((np.ptp(positions) + 2.0 * extent) / (fwhm_cm1 / 25.0)) + 1),
        )
    else:
        grid = _finite_real("grid_cm1", grid_cm1)
        if grid.ndim != 1 or grid.size < 2 or np.any(np.diff(grid) <= 0.0):
            raise ValueError("grid_cm1 must be a strictly increasing vector")
    intensity = np.zeros_like(grid)
    for line in lines:
        delta = grid - line.position_cm1
        if broadening_kind == "gaussian":
            profile = sqrt(4.0 * log(2.0) / pi) / fwhm_cm1 * np.exp(
                -4.0 * log(2.0) * (delta / fwhm_cm1) ** 2
            )
        else:
            half_width = 0.5 * fwhm_cm1
            profile = half_width / (pi * (delta * delta + half_width * half_width))
        intensity += line.normalized_strength * profile
    if normalization == "sum":
        area = _trapezoid(intensity, grid)
        if area <= 0.0:
            raise ValueError("broadened spectrum has nonpositive numerical area")
        intensity /= area
    return VibronicSpectrum(
        lines=lines,
        wavenumbers_cm1=grid,
        intensity=intensity,
        temperature_kelvin=temperature_kelvin,
        broadening_kind=broadening_kind,
        fwhm_cm1=fwhm_cm1,
        normalization=normalization,
        retained_thermal_population=retained_population,
        franck_condon_completeness=completeness,
    )


@dataclass(frozen=True)
class InfraredResult:
    intensities_km_mol: FloatArray
    mode_dipole_derivatives: ComplexArray

    def to_dict(self) -> dict[str, object]:
        return {
            "model": "target-excited-state double-harmonic infrared",
            "intensities": self.intensities_km_mol.tolist(),
            "intensity_unit": "km/mol",
            "mode_dipole_derivative_unit": "e/sqrt(amu)",
            "mode_dipole_derivatives": _complex_json(
                self.mode_dipole_derivatives
            ),
        }


def excited_state_ir_intensities(
    derivative: ExcitedStatePropertyDerivative,
) -> InfraredResult:
    """Return double-harmonic IR intensities from a target-state derivative."""

    _validate_excited_state_role(derivative.electronic_state_role)
    if derivative.property_kind != "state_dipole":
        raise ValueError("IR intensities require a state_dipole derivative")
    if derivative.coordinate_basis != "excited_normal":
        raise ValueError("IR derivatives must use excited_normal coordinates")
    if (
        derivative.coordinate_unit != "sqrt(amu)*bohr"
        or derivative.property_unit != "e*bohr"
    ):
        raise ValueError(
            "IR derivatives must use e*bohr per sqrt(amu)*bohr"
        )
    derivative_values = _finite_numeric("state_dipole_derivative", derivative.values)
    if derivative_values.ndim != 2 or derivative_values.shape[0] != 3:
        raise ValueError("state-dipole derivative must have shape (3, nmode)")
    if np.max(np.abs(derivative_values.imag), initial=0.0) > 1.0e-12:
        raise ValueError("nonresonant excited-state dipole derivatives must be real")
    values = derivative_values.T
    intensities = IR_INTENSITY_CONVERSION_KM_MOL * np.einsum(
        "ij,ij->i", values.conj(), values
    ).real
    return InfraredResult(intensities.astype(float), values.copy())


@dataclass(frozen=True)
class RamanResult:
    activities_au: FloatArray
    depolarization_ratio_polarized: FloatArray
    depolarization_ratio_unpolarized: FloatArray
    isotropic_invariant_squared: FloatArray
    anisotropic_invariant_squared: FloatArray
    mode_polarizability_derivatives: ComplexArray

    def to_dict(self) -> dict[str, object]:
        return {
            "model": "target-excited-state nonresonant double-harmonic Raman",
            "activities": self.activities_au.tolist(),
            "activity_unit": "bohr^4/amu",
            "depolarization_ratio_polarized": (
                self.depolarization_ratio_polarized.tolist()
            ),
            "depolarization_ratio_unpolarized": (
                self.depolarization_ratio_unpolarized.tolist()
            ),
            "isotropic_invariant_squared": (
                self.isotropic_invariant_squared.tolist()
            ),
            "anisotropic_invariant_squared": (
                self.anisotropic_invariant_squared.tolist()
            ),
            "mode_polarizability_derivatives": _complex_json(
                self.mode_polarizability_derivatives
            ),
        }


def _raman_invariants(tensors: ComplexArray) -> tuple[FloatArray, FloatArray]:
    symmetric = 0.5 * (tensors + np.swapaxes(tensors, 1, 2))
    asymmetry = float(
        np.max(np.abs(tensors - np.swapaxes(tensors, 1, 2)), initial=0.0)
    )
    scale = max(1.0, float(np.max(np.abs(tensors), initial=0.0)))
    if asymmetry / scale > 1.0e-8:
        raise ValueError(
            "polarizability derivatives are not symmetric within 1e-8 relative residual"
        )
    axx, ayy, azz = symmetric[:, 0, 0], symmetric[:, 1, 1], symmetric[:, 2, 2]
    axy, ayz, azx = symmetric[:, 0, 1], symmetric[:, 1, 2], symmetric[:, 2, 0]
    mean = (axx + ayy + azz) / 3.0
    anisotropic = 0.5 * (
        np.abs(axx - ayy) ** 2
        + np.abs(ayy - azz) ** 2
        + np.abs(azz - axx) ** 2
    ) + 3.0 * (np.abs(axy) ** 2 + np.abs(ayz) ** 2 + np.abs(azx) ** 2)
    return np.abs(mean) ** 2, anisotropic.real


def excited_state_raman_activities(
    derivative: ExcitedStatePropertyDerivative,
) -> RamanResult:
    """Return nonresonant Raman activities and depolarization ratios."""

    _validate_excited_state_role(derivative.electronic_state_role)
    if derivative.property_kind != "polarizability":
        raise ValueError("Raman activities require a polarizability derivative")
    if derivative.coordinate_basis != "excited_normal":
        raise ValueError("Raman derivatives must use excited_normal coordinates")
    if (
        derivative.coordinate_unit != "sqrt(amu)*bohr"
        or derivative.property_unit != "bohr^3"
    ):
        raise ValueError(
            "Raman derivatives must use bohr^3 per sqrt(amu)*bohr"
        )
    derivative_values = _finite_numeric(
        "polarizability_derivative", derivative.values
    )
    if derivative_values.ndim != 3 or derivative_values.shape[:2] != (3, 3):
        raise ValueError(
            "polarizability derivative must have shape (3, 3, nmode)"
        )
    if np.max(np.abs(derivative_values.imag), initial=0.0) > 1.0e-12:
        raise ValueError(
            "nonresonant excited-state polarizability derivatives must be real"
        )
    tensors = np.moveaxis(derivative_values, -1, 0)
    isotropic, anisotropic = _raman_invariants(tensors)
    activities = 45.0 * isotropic + 7.0 * anisotropic
    polarized_denominator = 45.0 * isotropic + 4.0 * anisotropic
    unpolarized_denominator = activities
    polarized = np.divide(
        3.0 * anisotropic,
        polarized_denominator,
        out=np.zeros_like(activities),
        where=polarized_denominator > 0.0,
    )
    unpolarized = np.divide(
        6.0 * anisotropic,
        unpolarized_denominator,
        out=np.zeros_like(activities),
        where=unpolarized_denominator > 0.0,
    )
    return RamanResult(
        activities_au=activities.astype(float),
        depolarization_ratio_polarized=polarized.astype(float),
        depolarization_ratio_unpolarized=unpolarized.astype(float),
        isotropic_invariant_squared=isotropic.astype(float),
        anisotropic_invariant_squared=anisotropic.astype(float),
        mode_polarizability_derivatives=tensors.copy(),
    )


@dataclass(frozen=True)
class ResonanceRamanResult:
    final_states: tuple[tuple[int, ...], ...]
    tensors: ComplexArray
    tensor_norm_squared: FloatArray
    incident_wavenumber_cm1: float
    damping_cm1: float
    intermediate_transition_completeness: float
    approximation: str

    def to_dict(self) -> dict[str, object]:
        return {
            "model": self.approximation,
            "final_states": [list(state) for state in self.final_states],
            "tensors": _complex_json(self.tensors),
            "tensor_norm_squared": self.tensor_norm_squared.tolist(),
            "incident_wavenumber_cm1": self.incident_wavenumber_cm1,
            "damping_cm1": self.damping_cm1,
            "intermediate_transition_completeness": (
                self.intermediate_transition_completeness
            ),
            "tensor_unit": "(e*bohr)^2/cm^-1",
            "absolute_cross_section": None,
        }


def resonance_raman_fc_ht(
    model: HarmonicVibronicModel,
    *,
    electronic_origin_cm1: float,
    origin_kind: OriginKind,
    incident_wavenumber_cm1: float,
    damping_cm1: float,
    transition: ElectronicTransitionMoment,
    transition_dipole_derivative: ExcitedStatePropertyDerivative | None,
    max_intermediate_quanta: int,
    final_states: Sequence[Sequence[int]] | None = None,
    minimum_intermediate_completeness: float = 0.999,
    max_states: int = 250000,
) -> ResonanceRamanResult:
    """Resonant-term KHD sum with harmonic FC/first-order HT moments.

    This proof-of-principle omits the antiresonant KHD term, non-Condon terms
    beyond first order, Duschinsky dependence of electronic properties, and
    lifetime variation.  It returns complex scattering tensors, not absolute
    Raman cross sections.
    """

    engine = HarmonicOverlapEngine(model)
    model = engine.model
    if origin_kind not in ("zero_zero", "adiabatic_minima"):
        raise ValueError("unsupported origin_kind")
    if not np.all(
        np.isfinite(
            [electronic_origin_cm1, incident_wavenumber_cm1, damping_cm1]
        )
    ):
        raise ValueError("electronic origin, incident wavenumber, and damping must be finite")
    if incident_wavenumber_cm1 <= 0.0 or damping_cm1 <= 0.0:
        raise ValueError("incident wavenumber and damping must be positive")
    dipole, gradient = _validate_transition_inputs(
        model, transition, transition_dipole_derivative
    )
    assert dipole is not None
    if not 0.0 <= minimum_intermediate_completeness <= 1.0:
        raise ValueError("minimum_intermediate_completeness must lie in [0, 1]")
    initial = (0,) * model.nmode
    if final_states is None:
        requested_final = []
        for mode in range(model.nmode):
            state = [0] * model.nmode
            state[mode] = 1
            requested_final.append(tuple(state))
        finals = tuple(requested_final)
    else:
        finals = tuple(
            _validate_state("final", state, model.nmode) for state in final_states
        )
    intermediate = enumerate_vibrational_states(
        model.nmode, max_intermediate_quanta, max_states=max_states
    )
    tensors = np.zeros((len(finals), 3, 3), dtype=complex)
    retained_transition_strength = 0.0
    for state in intermediate:
        absorption = vibronic_transition_moment(
            engine, initial, state, transition, transition_dipole_derivative
        )
        retained_transition_strength += float(np.vdot(absorption, absorption).real)
        if origin_kind == "zero_zero":
            vibronic_energy = electronic_origin_cm1 + float(
                np.dot(state, model.excited_frequencies_cm1)
            )
        else:
            vibronic_energy = electronic_origin_cm1 + float(
                np.dot(
                    np.asarray(state) + 0.5,
                    model.excited_frequencies_cm1,
                )
                - 0.5 * np.sum(model.ground_frequencies_cm1)
            )
        denominator = (
            vibronic_energy - incident_wavenumber_cm1 - 1j * damping_cm1
        )
        for index, final in enumerate(finals):
            emission_conjugate = np.conj(
                vibronic_transition_moment(
                    engine, final, state, transition, transition_dipole_derivative
                )
            )
            tensors[index] += np.outer(emission_conjugate, absorption) / denominator
    complete_transition_strength = float(np.vdot(dipole, dipole).real)
    if gradient is not None:
        ground_omega = model.ground_frequencies_cm1 * CM1_TO_HARTREE
        complete_transition_strength += float(
            np.sum(np.abs(gradient) ** 2 / (2.0 * ground_omega[None, :]))
        )
    if complete_transition_strength == 0.0:
        intermediate_completeness = 1.0
    else:
        intermediate_completeness = (
            retained_transition_strength / complete_transition_strength
        )
    if intermediate_completeness + 1.0e-15 < minimum_intermediate_completeness:
        raise ValueError(
            f"intermediate transition-strength completeness "
            f"{intermediate_completeness:.8f} is below the required "
            f"{minimum_intermediate_completeness:.8f}; increase "
            "max_intermediate_quanta"
        )
    norms = np.einsum("ijk,ijk->i", tensors.conj(), tensors).real
    return ResonanceRamanResult(
        final_states=finals,
        tensors=tensors,
        tensor_norm_squared=norms.astype(float),
        incident_wavenumber_cm1=float(incident_wavenumber_cm1),
        damping_cm1=float(damping_cm1),
        intermediate_transition_completeness=float(intermediate_completeness),
        approximation=(
            "resonant-term Kramers-Heisenberg-Dirac harmonic FC/first-order HT sum"
        ),
    )


def _complex_json(values: ComplexArray) -> object:
    if np.max(np.abs(values.imag), initial=0.0) < 1.0e-14:
        return values.real.tolist()
    return {"real": values.real.tolist(), "imag": values.imag.tolist()}


def main(argv: Sequence[str] | None = None) -> int:
    """Calculate a small FC spectrum from an explicit JSON record."""

    parser = argparse.ArgumentParser(
        description="Harmonic Franck-Condon spectrum from explicit J, K, and frequencies"
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args(argv)
    data = json.loads(arguments.input.read_text())
    model = HarmonicVibronicModel.create(
        data["ground_frequencies_cm1"],
        data["excited_frequencies_cm1"],
        data["J"],
        data["K"],
        coordinate_unit=data["coordinate_unit"],
        coordinate_phase_convention=data["coordinate_phase_convention"],
        orthogonality_tolerance=data.get("orthogonality_tolerance", 1.0e-8),
    )
    spectrum = harmonic_vibronic_spectrum(
        model,
        electronic_origin_cm1=data["electronic_origin_cm1"],
        origin_kind=data["origin_kind"],
        max_final_quanta=data["max_final_quanta"],
        temperature_kelvin=data.get("temperature_kelvin", 0.0),
        max_initial_quanta=data.get("max_initial_quanta"),
        broadening_kind=data.get("broadening_kind", "gaussian"),
        fwhm_cm1=data.get("fwhm_cm1", 100.0),
        normalization=data.get("normalization", "sum"),
        minimum_thermal_population=data.get("minimum_thermal_population", 0.999),
        minimum_franck_condon_completeness=data.get(
            "minimum_franck_condon_completeness", 0.999
        ),
    )
    output = {
        "model": "multidimensional harmonic Franck-Condon",
        "retained_thermal_population": spectrum.retained_thermal_population,
        "franck_condon_completeness": spectrum.franck_condon_completeness,
        "normalization": spectrum.normalization,
        "lines": [
            {
                "initial_state": list(line.initial_state),
                "final_state": list(line.final_state),
                "position_cm1": line.position_cm1,
                "overlap": line.overlap,
                "franck_condon_factor": line.franck_condon_factor,
                "thermal_population": line.thermal_population,
                "transition_moment": None
                if line.transition_moment is None
                else _complex_json(line.transition_moment),
                "raw_strength": line.raw_strength,
                "normalized_strength": line.normalized_strength,
            }
            for line in spectrum.lines
        ],
        "spectrum": {
            "wavenumbers_cm1": spectrum.wavenumbers_cm1.tolist(),
            "intensity": spectrum.intensity.tolist(),
            "broadening_kind": spectrum.broadening_kind,
            "fwhm_cm1": spectrum.fwhm_cm1,
        },
    }
    arguments.output.write_text(json.dumps(output, indent=2) + "\n")
    return 0


__all__ = (
    "ElectronicTransitionMoment",
    "ExcitedStatePropertyDerivative",
    "HarmonicOverlapEngine",
    "HarmonicVibronicModel",
    "InfraredResult",
    "RamanResult",
    "ResonanceRamanResult",
    "VibronicLine",
    "VibronicSpectrum",
    "enumerate_vibrational_states",
    "excited_state_ir_intensities",
    "excited_state_raman_activities",
    "finite_difference_excited_state_property",
    "harmonic_vibronic_spectrum",
    "mrsf_analytic_property_derivative_from_provider",
    "resonance_raman_fc_ht",
    "vibronic_transition_moment",
)


if __name__ == "__main__":
    raise SystemExit(main())
