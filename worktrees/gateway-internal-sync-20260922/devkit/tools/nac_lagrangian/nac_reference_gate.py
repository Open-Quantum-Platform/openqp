"""Gauge-resolved comparison gate for frozen NAC ``npz`` artifacts.

The Davidson eigenvector phase of each state is process-random.  A valid
comparison may therefore differ by one sign per state, but the signs of the
unordered state pairs cannot be chosen independently.  This gate resolves the
state gauge, rejects inconsistent pair signs, and then applies an absolute
component tolerance to every pair.

Examples
--------
Check finite-difference convergence::

    python nac_reference_gate.py coarse_dnum.npz fine_dnum.npz \
        --component-atol 1.1e-5 --energy-atol 1e-8 --require-flags

Check an analytic result against the converged numerical reference::

    python nac_reference_gate.py fine_dnum.npz analytic.npz \
        --component-atol 3e-4

The same command applies to H2O, ethylene, and Acrolein artifacts.  Worker
flags are checked whenever present; ``--require-flags`` additionally rejects
an artifact that omits them.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys
from typing import Mapping

import numpy as np


class NACGateError(ValueError):
    """Raised when a frozen NAC artifact violates the gate contract."""


@dataclass(frozen=True)
class PairMetric:
    istate: int
    jstate: int
    phase: int
    reference_norm: float
    candidate_norm: float
    max_component_error: float
    l2_error: float


@dataclass(frozen=True)
class GateResult:
    label: str
    state_signs: tuple[int, ...]
    pair_metrics: tuple[PairMetric, ...]
    max_component_error: float
    max_energy_error: float | None
    reference_antisymmetry: float
    candidate_antisymmetry: float

    def to_jsonable(self) -> dict[str, object]:
        result = asdict(self)
        result["state_signs"] = list(self.state_signs)
        return result


def _finite_array(value: np.ndarray, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise NACGateError(f"{name} contains NaN or infinity")
    return array


def validate_dcv(dcv: np.ndarray, name: str, antisym_atol: float) -> float:
    """Validate shape, finiteness, zero diagonal, and pair antisymmetry."""
    array = _finite_array(dcv, name)
    if array.ndim < 3 or array.shape[0] != array.shape[1]:
        raise NACGateError(
            f"{name} must have shape (nstate,nstate,...); got {array.shape}"
        )
    nstate = array.shape[0]
    diagonal = max(
        float(np.max(np.abs(array[state, state])))
        for state in range(nstate)
    )
    if diagonal > antisym_atol:
        raise NACGateError(
            f"{name} diagonal max {diagonal:.8e} exceeds {antisym_atol:.8e}"
        )
    antisymmetry = float(
        np.max(np.abs(array + np.swapaxes(array, 0, 1)))
    )
    if antisymmetry > antisym_atol:
        raise NACGateError(
            f"{name} antisymmetry error {antisymmetry:.8e} exceeds "
            f"{antisym_atol:.8e}"
        )
    return antisymmetry


def _pair_phases(
    reference: np.ndarray,
    candidate: np.ndarray,
    zero_pair_atol: float,
) -> tuple[dict[tuple[int, int], int], set[tuple[int, int]]]:
    phases: dict[tuple[int, int], int] = {}
    active: set[tuple[int, int]] = set()
    for istate in range(reference.shape[0]):
        for jstate in range(istate + 1, reference.shape[0]):
            ref = reference[istate, jstate].reshape(-1)
            cand = candidate[istate, jstate].reshape(-1)
            ref_norm = float(np.linalg.norm(ref))
            cand_norm = float(np.linalg.norm(cand))
            pair = (istate, jstate)
            if max(ref_norm, cand_norm) <= zero_pair_atol:
                phases[pair] = 1
                continue
            active.add(pair)
            minus = float(np.linalg.norm(cand - ref))
            plus = float(np.linalg.norm(cand + ref))
            phases[pair] = 1 if minus <= plus else -1
    return phases, active


def resolve_state_gauge(
    phases: Mapping[tuple[int, int], int],
    active_pairs: set[tuple[int, int]],
    nstate: int,
) -> tuple[int, ...]:
    """Resolve pair phases into state signs and reject frustrated cycles."""
    state_signs: list[int | None] = [None] * nstate
    adjacency: list[list[tuple[int, int]]] = [[] for _ in range(nstate)]
    for istate, jstate in active_pairs:
        phase = int(phases[(istate, jstate)])
        adjacency[istate].append((jstate, phase))
        adjacency[jstate].append((istate, phase))

    for root in range(nstate):
        if state_signs[root] is not None:
            continue
        state_signs[root] = 1
        stack = [root]
        while stack:
            istate = stack.pop()
            assert state_signs[istate] is not None
            for jstate, phase in adjacency[istate]:
                expected = int(state_signs[istate]) * phase
                if state_signs[jstate] is None:
                    state_signs[jstate] = expected
                    stack.append(jstate)
                elif state_signs[jstate] != expected:
                    raise NACGateError(
                        "pair phases are not a state gauge: inconsistent "
                        f"cycle at pair ({istate + 1},{jstate + 1})"
                    )

    signs = tuple(int(sign) for sign in state_signs)
    for istate, jstate in active_pairs:
        expected = signs[istate] * signs[jstate]
        if phases[(istate, jstate)] != expected:
            raise NACGateError(
                "pair phases are not a state gauge: pair "
                f"({istate + 1},{jstate + 1}) has "
                f"{phases[(istate, jstate)]:+d}, expected {expected:+d}"
            )
    return signs


def _excited_energies(values: np.ndarray, nstate: int, name: str) -> np.ndarray:
    energies = _finite_array(values, name).reshape(-1)
    if energies.size == nstate:
        return energies
    if energies.size == nstate + 1:
        return energies[1:]
    raise NACGateError(
        f"{name} has {energies.size} values; expected {nstate} excited "
        f"states or ground+{nstate}"
    )


def _validate_optional_payload(
    payload: Mapping[str, np.ndarray],
    dcv: np.ndarray,
    name: str,
    identity_atol: float,
    require_flags: bool,
) -> np.ndarray | None:
    if require_flags and "flags" not in payload:
        raise NACGateError(f"{name} has no flags array")
    if "flags" in payload:
        flags = {
            flag.decode() if isinstance(flag, bytes) else str(flag)
            for flag in np.asarray(payload["flags"]).reshape(-1)
        }
        if not flags:
            raise NACGateError(f"{name} has an empty flags array")
        successful = {
            "computed",
            "loaded",
            "analytic-v2-reference",
            "analytic-v3-zvector",
        }
        unexpected = flags - successful
        if unexpected:
            raise NACGateError(
                f"{name} has unsuccessful flags: {sorted(unexpected)}"
            )

    nstate = dcv.shape[0]
    energies = None
    if "energies" in payload:
        energies = _excited_energies(payload["energies"], nstate, f"{name}.energies")

    if "nacv" in payload:
        nacv = _finite_array(payload["nacv"], f"{name}.nacv")
        if nacv.shape != dcv.shape:
            raise NACGateError(
                f"{name}.nacv shape {nacv.shape} differs from dcv {dcv.shape}"
            )
        symmetry = float(np.max(np.abs(nacv - np.swapaxes(nacv, 0, 1))))
        if symmetry > identity_atol:
            raise NACGateError(
                f"{name}.nacv symmetry error {symmetry:.8e} exceeds "
                f"{identity_atol:.8e}"
            )
        diagonal = max(
            float(np.max(np.abs(nacv[state, state])))
            for state in range(nstate)
        )
        if diagonal > identity_atol:
            raise NACGateError(
                f"{name}.nacv diagonal max {diagonal:.8e} exceeds "
                f"{identity_atol:.8e}"
            )
        if energies is not None:
            for istate in range(nstate):
                for jstate in range(nstate):
                    if istate == jstate:
                        continue
                    expected = (energies[jstate] - energies[istate]) * dcv[
                        istate, jstate
                    ]
                    error = float(np.max(np.abs(nacv[istate, jstate] - expected)))
                    if error > identity_atol:
                        raise NACGateError(
                            f"{name} gap identity error for pair "
                            f"({istate + 1},{jstate + 1}) is {error:.8e}, "
                            f"above {identity_atol:.8e}"
                        )
    return energies


def compare_payloads(
    reference_payload: Mapping[str, np.ndarray],
    candidate_payload: Mapping[str, np.ndarray],
    *,
    component_atol: float,
    reference_key: str = "dcv",
    candidate_key: str = "dcv",
    antisym_atol: float = 1.0e-10,
    energy_atol: float | None = None,
    identity_atol: float = 1.0e-8,
    zero_pair_atol: float = 1.0e-14,
    require_flags: bool = False,
    label: str = "NAC",
) -> GateResult:
    """Compare two NAC payloads after resolving their relative state gauge."""
    if component_atol < 0.0:
        raise NACGateError("component_atol must be non-negative")
    try:
        reference = _finite_array(reference_payload[reference_key], "reference.dcv")
    except KeyError as exc:
        raise NACGateError(f"reference has no {reference_key!r} array") from exc
    try:
        candidate = _finite_array(candidate_payload[candidate_key], "candidate.dcv")
    except KeyError as exc:
        raise NACGateError(f"candidate has no {candidate_key!r} array") from exc
    if reference.shape != candidate.shape:
        raise NACGateError(
            f"NAC shapes differ: reference {reference.shape}, candidate {candidate.shape}"
        )

    reference_asym = validate_dcv(reference, "reference.dcv", antisym_atol)
    candidate_asym = validate_dcv(candidate, "candidate.dcv", antisym_atol)
    reference_energies = _validate_optional_payload(
        reference_payload,
        reference,
        "reference",
        identity_atol,
        require_flags,
    )
    candidate_energies = _validate_optional_payload(
        candidate_payload,
        candidate,
        "candidate",
        identity_atol,
        require_flags,
    )

    phases, active = _pair_phases(reference, candidate, zero_pair_atol)
    signs = resolve_state_gauge(phases, active, reference.shape[0])
    metrics: list[PairMetric] = []
    largest = 0.0
    for istate in range(reference.shape[0]):
        for jstate in range(istate + 1, reference.shape[0]):
            phase = phases[(istate, jstate)]
            ref = reference[istate, jstate].reshape(-1)
            cand = candidate[istate, jstate].reshape(-1)
            error = phase * cand - ref
            max_error = float(np.max(np.abs(error)))
            largest = max(largest, max_error)
            metrics.append(
                PairMetric(
                    istate=istate + 1,
                    jstate=jstate + 1,
                    phase=phase,
                    reference_norm=float(np.linalg.norm(ref)),
                    candidate_norm=float(np.linalg.norm(cand)),
                    max_component_error=max_error,
                    l2_error=float(np.linalg.norm(error)),
                )
            )

    energy_error = None
    if energy_atol is not None and (
        reference_energies is None or candidate_energies is None
    ):
        raise NACGateError(
            "energy_atol requires energies arrays in both artifacts"
        )
    if reference_energies is not None and candidate_energies is not None:
        energy_error = float(np.max(np.abs(reference_energies - candidate_energies)))
        if energy_atol is not None and energy_error > energy_atol:
            raise NACGateError(
                f"{label} energy error {energy_error:.8e} exceeds "
                f"{energy_atol:.8e}"
            )
    if largest > component_atol:
        raise NACGateError(
            f"{label} max component error {largest:.8e} exceeds "
            f"{component_atol:.8e}"
        )

    return GateResult(
        label=label,
        state_signs=signs,
        pair_metrics=tuple(metrics),
        max_component_error=largest,
        max_energy_error=energy_error,
        reference_antisymmetry=reference_asym,
        candidate_antisymmetry=candidate_asym,
    )


def _load_npz(path: str | Path) -> dict[str, np.ndarray]:
    source = Path(path)
    if not source.is_file():
        raise NACGateError(f"artifact does not exist: {source}")
    with np.load(source, allow_pickle=False) as frozen:
        return {key: np.array(frozen[key], copy=True) for key in frozen.files}


def compare_files(
    reference_path: str | Path,
    candidate_path: str | Path,
    **kwargs: object,
) -> GateResult:
    return compare_payloads(
        _load_npz(reference_path),
        _load_npz(candidate_path),
        **kwargs,
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", help="reference/finer NAC npz")
    parser.add_argument("candidate", help="candidate/coarser NAC npz")
    parser.add_argument("--component-atol", type=float, required=True)
    parser.add_argument("--energy-atol", type=float)
    parser.add_argument("--antisym-atol", type=float, default=1.0e-10)
    parser.add_argument("--identity-atol", type=float, default=1.0e-8)
    parser.add_argument("--zero-pair-atol", type=float, default=1.0e-14)
    parser.add_argument("--reference-key", default="dcv")
    parser.add_argument("--candidate-key", default="dcv")
    parser.add_argument("--require-flags", action="store_true")
    parser.add_argument("--label", default="NAC")
    parser.add_argument("--json", action="store_true", dest="as_json")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        result = compare_files(
            args.reference,
            args.candidate,
            component_atol=args.component_atol,
            reference_key=args.reference_key,
            candidate_key=args.candidate_key,
            antisym_atol=args.antisym_atol,
            energy_atol=args.energy_atol,
            identity_atol=args.identity_atol,
            zero_pair_atol=args.zero_pair_atol,
            require_flags=args.require_flags,
            label=args.label,
        )
    except NACGateError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1

    if args.as_json:
        print(json.dumps(result.to_jsonable(), indent=2, sort_keys=True))
        return 0
    print(f"===== {result.label} gauge-resolved NAC gate =====")
    print("state gauge:", " ".join(f"{sign:+d}" for sign in result.state_signs))
    for metric in result.pair_metrics:
        print(
            f"({metric.istate},{metric.jstate}) phase={metric.phase:+d} "
            f"|ref|={metric.reference_norm:.10e} "
            f"|candidate|={metric.candidate_norm:.10e} "
            f"maxdiff={metric.max_component_error:.8e} "
            f"|diff|={metric.l2_error:.8e}"
        )
    if result.max_energy_error is not None:
        print(f"energy maxdiff={result.max_energy_error:.8e}")
    print(f"PASS: max component error={result.max_component_error:.8e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
