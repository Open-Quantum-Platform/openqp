"""Structural translation gate for analytic MRSF NAC debug artifacts.

The production analytic NAC is a *raw electronic* derivative coupling: no
electron-translation-factor (ETF) correction is applied.  Its atom sum is
therefore generally nonzero and must not be tested against zero.  For the
ordered-pair debug decomposition written by ``NAC_ANALYTIC_DEBUG``, rigid
translation instead requires

* the antisymmetrized ``t1``, ``z_hf``, ``z_xc``, and ``vmask`` atom sums to
  vanish; and
* the atom sum of the antisymmetrized full ``dp`` to equal that of ``gsk``
  (the ``gamma:Sk`` electronic-translation contribution).

Run::

    python tools/nac_lagrangian/translation_gate.py analytic_debug.npz

The tolerances apply to Cartesian components of the atom sums.  The reported
raw atom sum is diagnostic only and is deliberately never compared with zero.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Mapping

import numpy as np


NONOVERLAP_COMPONENTS = ("t1", "z_hf", "z_xc", "vmask")


class TranslationGateError(ValueError):
    """Raised when an analytic NAC debug artifact violates the sum rule."""


@dataclass(frozen=True)
class PairTranslationMetric:
    istate: int
    jstate: int
    raw_atom_sum: tuple[float, float, float]
    gamma_sk_atom_sum: tuple[float, float, float]
    term_atom_sums: dict[str, tuple[float, float, float]]
    max_nonoverlap_atom_sum: float
    identity_error: float


@dataclass(frozen=True)
class TranslationGateResult:
    pair_metrics: tuple[PairTranslationMetric, ...]
    max_nonoverlap_atom_sum: float
    max_identity_error: float
    max_raw_atom_sum: float


def _finite_array(value: np.ndarray, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise TranslationGateError(f"{name} contains NaN or infinity")
    return array


def _infer_nstate(payload: Mapping[str, np.ndarray]) -> int:
    if "Xf" in payload:
        amplitudes = _finite_array(payload["Xf"], "Xf")
        if amplitudes.ndim != 2:
            raise TranslationGateError(
                f"Xf must have shape (namplitude,nstate); got {amplitudes.shape}"
            )
        nstate = int(amplitudes.shape[1])
    elif "energies" in payload:
        energies = _finite_array(payload["energies"], "energies").reshape(-1)
        if energies.size < 3:
            raise TranslationGateError(
                "energies must contain a reference energy and at least two states"
            )
        nstate = int(energies.size - 1)
    else:
        raise TranslationGateError(
            "debug artifact needs Xf (preferred) or ground-plus-state energies"
        )

    if nstate < 2:
        raise TranslationGateError(f"at least two states are required; got {nstate}")
    # NAC_ANALYTIC_DEBUG currently writes keys as e.g. dp_01, without a
    # separator between state indices.  Indices >= 10 would collide.
    if nstate > 10:
        raise TranslationGateError(
            "NAC_ANALYTIC_DEBUG pair keys are ambiguous for more than 10 states"
        )
    return nstate


def _ordered_component(
    payload: Mapping[str, np.ndarray],
    component: str,
    istate: int,
    jstate: int,
) -> np.ndarray:
    key = f"{component}_{istate}{jstate}"
    if key not in payload:
        raise TranslationGateError(f"debug artifact is missing {key}")
    array = _finite_array(payload[key], key).reshape(-1)
    if array.size == 0 or array.size % 3 != 0:
        raise TranslationGateError(
            f"{key} must contain 3*natom values; got {array.size}"
        )
    return array.reshape(-1, 3)


def _antisymmetrized_component(
    payload: Mapping[str, np.ndarray],
    component: str,
    istate: int,
    jstate: int,
) -> np.ndarray:
    forward = _ordered_component(payload, component, istate, jstate)
    reverse = _ordered_component(payload, component, jstate, istate)
    if forward.shape != reverse.shape:
        raise TranslationGateError(
            f"{component} ordered-pair shapes differ for "
            f"({istate + 1},{jstate + 1}): {forward.shape} vs {reverse.shape}"
        )
    return 0.5 * (forward - reverse)


def analyze_payload(
    payload: Mapping[str, np.ndarray],
    *,
    term_atol: float = 1.0e-9,
    identity_atol: float = 1.0e-9,
) -> TranslationGateResult:
    """Check the raw-translation structure of one analytic debug payload."""
    if term_atol < 0.0 or identity_atol < 0.0:
        raise TranslationGateError("translation tolerances must be nonnegative")

    nstate = _infer_nstate(payload)
    metrics: list[PairTranslationMetric] = []
    max_nonoverlap = 0.0
    max_identity = 0.0
    max_raw = 0.0

    for istate in range(nstate):
        for jstate in range(istate + 1, nstate):
            dp = _antisymmetrized_component(
                payload, "dp", istate, jstate
            )
            gsk = _antisymmetrized_component(
                payload, "gsk", istate, jstate
            )
            if dp.shape != gsk.shape:
                raise TranslationGateError(
                    f"dp and gsk shapes differ for ({istate + 1},{jstate + 1}): "
                    f"{dp.shape} vs {gsk.shape}"
                )

            raw_atom_sum = np.sum(dp, axis=0)
            gamma_atom_sum = np.sum(gsk, axis=0)
            identity_error = float(
                np.max(np.abs(raw_atom_sum - gamma_atom_sum))
            )
            if identity_error > identity_atol:
                raise TranslationGateError(
                    "raw dp atom sum differs from gamma:Sk for pair "
                    f"({istate + 1},{jstate + 1}): {identity_error:.8e} "
                    f"> {identity_atol:.8e}"
                )

            term_sums: dict[str, tuple[float, float, float]] = {}
            pair_nonoverlap = 0.0
            for component in NONOVERLAP_COMPONENTS:
                values = _antisymmetrized_component(
                    payload, component, istate, jstate
                )
                if values.shape != dp.shape:
                    raise TranslationGateError(
                        f"{component} shape {values.shape} differs from dp "
                        f"{dp.shape} for pair ({istate + 1},{jstate + 1})"
                    )
                atom_sum = np.sum(values, axis=0)
                component_error = float(np.max(np.abs(atom_sum)))
                if component_error > term_atol:
                    raise TranslationGateError(
                        f"{component} atom sum for pair "
                        f"({istate + 1},{jstate + 1}) is "
                        f"{component_error:.8e} > {term_atol:.8e}"
                    )
                pair_nonoverlap = max(pair_nonoverlap, component_error)
                term_sums[component] = tuple(float(value) for value in atom_sum)

            raw_max = float(np.max(np.abs(raw_atom_sum)))
            max_nonoverlap = max(max_nonoverlap, pair_nonoverlap)
            max_identity = max(max_identity, identity_error)
            max_raw = max(max_raw, raw_max)
            metrics.append(
                PairTranslationMetric(
                    istate=istate + 1,
                    jstate=jstate + 1,
                    raw_atom_sum=tuple(float(value) for value in raw_atom_sum),
                    gamma_sk_atom_sum=tuple(
                        float(value) for value in gamma_atom_sum
                    ),
                    term_atom_sums=term_sums,
                    max_nonoverlap_atom_sum=pair_nonoverlap,
                    identity_error=identity_error,
                )
            )

    return TranslationGateResult(
        pair_metrics=tuple(metrics),
        max_nonoverlap_atom_sum=max_nonoverlap,
        max_identity_error=max_identity,
        max_raw_atom_sum=max_raw,
    )


def analyze_file(
    path: str | Path,
    *,
    term_atol: float = 1.0e-9,
    identity_atol: float = 1.0e-9,
) -> TranslationGateResult:
    """Load an ``NAC_ANALYTIC_DEBUG`` NPZ and apply the structural gate."""
    with np.load(Path(path), allow_pickle=False) as frozen:
        payload = {key: np.array(frozen[key], copy=True) for key in frozen.files}
    return analyze_payload(
        payload,
        term_atol=term_atol,
        identity_atol=identity_atol,
    )


def _format_vector(values: tuple[float, float, float]) -> str:
    return "[" + " ".join(f"{value:+.8e}" for value in values) + "]"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Gate raw analytic NAC translation structure (no ETF)."
    )
    parser.add_argument("debug_npz", type=Path)
    parser.add_argument("--term-atol", type=float, default=1.0e-9)
    parser.add_argument("--identity-atol", type=float, default=1.0e-9)
    args = parser.parse_args(argv)

    try:
        result = analyze_file(
            args.debug_npz,
            term_atol=args.term_atol,
            identity_atol=args.identity_atol,
        )
    except (OSError, TranslationGateError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1

    print("===== raw analytic NAC translation gate (no ETF) =====")
    for metric in result.pair_metrics:
        print(
            f"({metric.istate},{metric.jstate}) "
            f"raw atom sum={_format_vector(metric.raw_atom_sum)} "
            f"gamma:Sk={_format_vector(metric.gamma_sk_atom_sum)} "
            f"identity={metric.identity_error:.8e} "
            f"nonoverlap={metric.max_nonoverlap_atom_sum:.8e}"
        )
    print(
        "raw atom-sum magnitude (diagnostic, not zero-gated): "
        f"{result.max_raw_atom_sum:.8e}"
    )
    print(
        "PASS: max non-overlap atom sum "
        f"{result.max_nonoverlap_atom_sum:.8e}; "
        f"max dp/gamma:Sk identity error {result.max_identity_error:.8e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
