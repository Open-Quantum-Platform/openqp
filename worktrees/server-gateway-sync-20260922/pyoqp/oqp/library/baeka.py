"""Generalized Baek adaptive-penalty kernel for conical intersections.

The published three-state objective removes the redundant outer-state gap and
penalizes only the two consecutive gaps.  This module extends that construction
to ``N >= 2`` states by using the ``N - 1`` adjacent gaps, together with an SVD
projector for their independent gradient subspace.  It contains no OpenQP
backend calls and is therefore directly unit-testable.

The multiplier updates are additive, as in Baek's equations 3.12--3.13::

    sigma <- sigma + delta_beta        (ordinary macro iteration)
    sigma <- sigma + beta              (locally stationary, gap still open)

``alpha`` and all gaps are energies (Hartree); ``sigma``, ``delta_beta`` and
``beta`` are dimensionless.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Tuple, Union

import numpy as np


ArrayLike = Union[np.ndarray, Iterable[float]]


def parse_jump_schedule(value) -> Tuple[float, ...]:
    """Return a validated additive ``beta`` schedule.

    The sectioned input reader may deliver the schedule as a comma-separated
    string, while the concise parser and Python API can provide a sequence.
    A scalar is accepted as a one-jump schedule; exhaustion is deliberately an
    error rather than silently inventing a continuation rule.
    """

    if isinstance(value, str):
        fields = [field.strip() for field in value.split(",")]
        if any(not field for field in fields):
            raise ValueError("pen_jump must be a comma-separated list of numbers")
        values = tuple(float(field) for field in fields)
    elif np.isscalar(value):
        values = (float(value),)
    else:
        values = tuple(float(item) for item in value)

    if not values:
        raise ValueError("pen_jump must contain at least one additive beta")
    if not np.all(np.isfinite(values)) or any(item <= 0.0 for item in values):
        raise ValueError("pen_jump values must be finite and positive")
    return values


def penalty_subspace(
    gap_gradients: ArrayLike,
    *,
    rtol: float = 1.0e-10,
    atol: float = 1.0e-12,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return an orthonormal basis and singular values for the gap subspace.

    ``gap_gradients`` has shape ``(N - 1, 3 * natom)``.  Raw adjacent
    gap-gradient vectors are used instead of multiplying by the penalty slope:
    the column space is identical away from exact degeneracy, while the raw
    vectors remain well-defined as the gap and the penalty slope approach zero.
    """

    gaps = np.asarray(gap_gradients, dtype=float)
    if gaps.ndim != 2 or gaps.shape[0] < 1 or gaps.shape[1] < 1:
        raise ValueError("gap_gradients must have shape (N-1, ncoord)")
    if not np.all(np.isfinite(gaps)):
        raise ValueError("gap_gradients must contain only finite values")
    if not np.isfinite(rtol) or not np.isfinite(atol) or rtol < 0.0 or atol < 0.0:
        raise ValueError("projector tolerances must be finite and non-negative")

    left, singular, _ = np.linalg.svd(gaps.T, full_matrices=False)
    scale = float(singular[0]) if singular.size else 0.0
    keep = singular > max(float(atol), float(rtol) * scale)
    return left[:, keep], singular


def penalty_projector(
    gap_gradients: ArrayLike,
    *,
    rtol: float = 1.0e-10,
    atol: float = 1.0e-12,
) -> Tuple[np.ndarray, int, np.ndarray]:
    """Return ``P = U U.T``, its numerical rank, and singular values."""

    basis, singular = penalty_subspace(
        gap_gradients, rtol=rtol, atol=atol
    )
    projector = basis @ basis.T
    return projector, int(basis.shape[1]), singular


@dataclass(frozen=True)
class BaekAObjective:
    """One generalized Baek objective/gradient evaluation."""

    energies: np.ndarray
    gradients: np.ndarray
    gaps: np.ndarray
    gap_gradients: np.ndarray
    average_energy: float
    average_gradient: np.ndarray
    penalty: float
    penalty_gradient: np.ndarray
    objective: float
    gradient: np.ndarray
    span: float
    sigma: float
    effective_sigma: float
    alpha: float
    projector_basis: np.ndarray
    projector_singular_values: np.ndarray
    projector_rank: int
    parallel_gradient: np.ndarray
    perpendicular_gradient: np.ndarray
    parallel_norm: float
    parallel_scaled_norm: float
    perpendicular_norm: float


def baeka_objective(
    energies: ArrayLike,
    gradients: ArrayLike,
    *,
    sigma: float,
    alpha: float,
    gap_weight: float = 1.0,
    projector_rtol: float = 1.0e-10,
    projector_atol: float = 1.0e-12,
) -> BaekAObjective:
    r"""Evaluate the generalized adjacent-gap Baek objective.

    .. math::

       F_N = N^{-1}\sum_i E_i + w\sigma\sum_i
             \frac{(E_{i+1}-E_i)^2}{E_{i+1}-E_i+\alpha}.

    States must be supplied in ascending adiabatic/root order.  Only adjacent
    gaps are included; adding all pairs would reintroduce the redundancy that
    the three-state Baek objective was designed to remove.
    """

    energies = np.asarray(energies, dtype=float).reshape(-1)
    gradients = np.asarray(gradients, dtype=float)
    if energies.size < 2:
        raise ValueError("BaekA requires at least two states")
    if gradients.ndim < 2 or gradients.shape[0] != energies.size:
        raise ValueError("gradients must have one row per selected state")
    gradients = gradients.reshape(energies.size, -1)
    if gradients.shape[1] < 1:
        raise ValueError("gradients must contain at least one coordinate")
    if not np.all(np.isfinite(energies)) or not np.all(np.isfinite(gradients)):
        raise ValueError("BaekA energies and gradients must be finite")

    sigma = float(sigma)
    alpha = float(alpha)
    gap_weight = float(gap_weight)
    if not np.isfinite(sigma) or sigma <= 0.0:
        raise ValueError("BaekA sigma must be finite and positive")
    if not np.isfinite(alpha) or alpha <= 0.0:
        raise ValueError("BaekA alpha must be a finite positive energy")
    if not np.isfinite(gap_weight) or gap_weight != 1.0:
        raise ValueError(
            "BaekA gap_weight is fixed at 1.0; use sigma for the published "
            "penalty multiplier"
        )

    gaps = np.diff(energies)
    energy_scale = max(1.0, float(np.max(np.abs(energies))))
    ordering_tolerance = 64.0 * np.finfo(float).eps * energy_scale
    if np.any(gaps < -ordering_tolerance):
        raise ValueError("BaekA states must be in ascending energy/root order")
    # Remove harmless negative roundoff at a numerical degeneracy without
    # accepting a genuinely misordered state list.
    gaps = np.maximum(gaps, 0.0)
    gap_gradients = np.diff(gradients, axis=0)

    denominator = gaps + alpha
    penalty_terms = gaps * gaps / denominator
    slopes = (gaps * gaps + 2.0 * alpha * gaps) / (denominator * denominator)

    average_energy = float(np.mean(energies))
    average_gradient = np.mean(gradients, axis=0)
    penalty = float(np.sum(penalty_terms))
    penalty_gradient = np.sum(slopes[:, None] * gap_gradients, axis=0)
    effective_sigma = sigma * gap_weight
    objective = average_energy + effective_sigma * penalty
    gradient = average_gradient + effective_sigma * penalty_gradient

    basis, singular = penalty_subspace(
        gap_gradients, rtol=projector_rtol, atol=projector_atol
    )
    if basis.shape[1]:
        parallel = basis @ (basis.T @ gradient)
    else:
        parallel = np.zeros_like(gradient)
    perpendicular = gradient - parallel
    parallel_norm = float(np.linalg.norm(parallel))

    return BaekAObjective(
        energies=energies.copy(),
        gradients=gradients.copy(),
        gaps=gaps,
        gap_gradients=gap_gradients,
        average_energy=average_energy,
        average_gradient=average_gradient,
        penalty=penalty,
        penalty_gradient=penalty_gradient,
        objective=float(objective),
        gradient=gradient,
        span=float(np.sum(gaps)),
        sigma=sigma,
        effective_sigma=effective_sigma,
        alpha=alpha,
        projector_basis=basis,
        projector_singular_values=singular,
        projector_rank=int(basis.shape[1]),
        parallel_gradient=parallel,
        perpendicular_gradient=perpendicular,
        parallel_norm=parallel_norm,
        parallel_scaled_norm=parallel_norm / sigma,
        perpendicular_norm=float(np.linalg.norm(perpendicular)),
    )


@dataclass(frozen=True)
class BaekAEvaluation:
    """Objective data plus one transition of the additive-sigma state machine."""

    tested_data: BaekAObjective
    data: BaekAObjective
    tested_sigma: float
    next_sigma: float
    delta_f: float
    local_stationary: bool
    gap_converged: bool
    converged: bool
    action: str
    jump: Optional[float]
    jump_index: int


class BaekAJumpScheduleExhausted(RuntimeError):
    """Raised when another large additive jump is required but unavailable."""


class BaekAState:
    """State machine for the generalized Baek adaptive multiplier."""

    def __init__(
        self,
        *,
        sigma: float,
        alpha: float,
        delta_beta: float,
        beta_schedule,
        gap_tol: float,
        tol_f: float,
        tol_g: float,
        gap_weight: float = 1.0,
        projector_rtol: float = 1.0e-10,
        projector_atol: float = 1.0e-12,
    ):
        self.sigma = float(sigma)
        self.alpha = float(alpha)
        self.delta_beta = float(delta_beta)
        self.beta_schedule = parse_jump_schedule(beta_schedule)
        self.gap_tol = float(gap_tol)
        self.tol_f = float(tol_f)
        self.tol_g = float(tol_g)
        self.gap_weight = float(gap_weight)
        self.projector_rtol = float(projector_rtol)
        self.projector_atol = float(projector_atol)
        self.jump_index = 0
        self._previous = None

        if not np.isfinite(self.sigma) or self.sigma <= 0.0:
            raise ValueError("BaekA sigma must be finite and positive")
        if not np.isfinite(self.alpha) or self.alpha <= 0.0:
            raise ValueError("BaekA alpha must be a finite positive energy")
        if not np.isfinite(self.delta_beta) or self.delta_beta <= 0.0:
            raise ValueError("BaekA pen_delta must be finite and positive")
        if not np.isfinite(self.gap_tol) or self.gap_tol <= 0.0:
            raise ValueError("BaekA gap tolerance must be finite and positive")
        if not np.isfinite(self.tol_f) or self.tol_f <= 0.0:
            raise ValueError("BaekA objective tolerance must be finite and positive")
        if not np.isfinite(self.tol_g) or self.tol_g <= 0.0:
            raise ValueError("BaekA gradient tolerance must be finite and positive")
        if not np.isfinite(self.gap_weight) or self.gap_weight != 1.0:
            raise ValueError(
                "BaekA gap_weight is fixed at 1.0; use sigma for the published "
                "penalty multiplier"
            )

    def _objective(self, energies, gradients, sigma) -> BaekAObjective:
        return baeka_objective(
            energies,
            gradients,
            sigma=sigma,
            alpha=self.alpha,
            gap_weight=self.gap_weight,
            projector_rtol=self.projector_rtol,
            projector_atol=self.projector_atol,
        )

    def evaluate(self, energies: ArrayLike, gradients: ArrayLike) -> BaekAEvaluation:
        """Evaluate the objective and advance the sigma state machine once."""

        tested_sigma = self.sigma
        tested_data = self._objective(energies, gradients, tested_sigma)
        data = tested_data

        if self._previous is None:
            delta_f = float("inf")
            local_stationary = False
        else:
            previous_at_current_sigma = (
                self._previous["average_energy"]
                + tested_sigma * self.gap_weight * self._previous["penalty"]
            )
            delta_f = abs(tested_data.objective - previous_at_current_sigma)
            local_stationary = bool(
                delta_f <= self.tol_f
                and tested_data.parallel_scaled_norm <= self.tol_g
                and tested_data.perpendicular_norm <= self.tol_g
            )

        gap_converged = bool(tested_data.span <= self.gap_tol)
        jump = None

        if local_stationary and gap_converged:
            action = "converged"
            converged = True
        elif local_stationary:
            if self.jump_index >= len(self.beta_schedule):
                raise BaekAJumpScheduleExhausted(
                    "BaekA requires another large beta jump, but pen_jump "
                    "has been exhausted"
                )
            jump = self.beta_schedule[self.jump_index]
            self.jump_index += 1
            self.sigma = tested_sigma + jump
            if not np.isfinite(self.sigma):
                raise FloatingPointError("BaekA sigma overflowed after beta jump")
            # The electronic energies/gradients are unchanged at this geometry;
            # re-form the objective immediately so the optimizer's next step is
            # taken on the strengthened objective, not the stale pre-jump one.
            data = self._objective(energies, gradients, self.sigma)
            action = "jump"
            converged = False
        else:
            self.sigma = tested_sigma + self.delta_beta
            if not np.isfinite(self.sigma):
                raise FloatingPointError("BaekA sigma overflowed after pen_delta")
            action = "step"
            converged = False

        # Store sigma-independent components.  At the next geometry they can be
        # recombined at that iteration's sigma, which makes the Eq. 3.8 objective
        # change a same-objective comparison even though sigma evolves.
        self._previous = {
            "average_energy": data.average_energy,
            "penalty": data.penalty,
        }

        return BaekAEvaluation(
            tested_data=tested_data,
            data=data,
            tested_sigma=tested_sigma,
            next_sigma=self.sigma,
            delta_f=float(delta_f),
            local_stationary=local_stationary,
            gap_converged=gap_converged,
            converged=converged,
            action=action,
            jump=jump,
            jump_index=self.jump_index,
        )
