#!/usr/bin/env python3
"""Record the native NAMD hop kernel without changing its numerical path.

This diagnostic runner is deliberately a thin Python observer.  The TLF,
electronic propagation, and FSSH probability construction remain in the
production Fortran implementation.  ``--skip-first-hop`` emulates the
historical KNU-GAMESS first-interval convention for controlled comparisons.
"""

from __future__ import annotations

import argparse
import csv
import inspect
from pathlib import Path

import numpy as np


class SequenceRNG:
    """Minimal Generator-compatible source for a prescribed hop RNG stream."""

    def __init__(self, values: np.ndarray):
        self.values = np.asarray(values, dtype=float).reshape(-1)
        self.index = 0

    def random(self, size=None):
        if size is not None:
            raise RuntimeError("Prescribed NAMD hop RNG only supports scalar draws")
        if self.index >= self.values.size:
            raise RuntimeError(
                f"Prescribed NAMD hop RNG exhausted after {self.index} draws"
            )
        value = float(self.values[self.index])
        self.index += 1
        return value


def build_trace_row(*, istep, dt_fs, active, hopped, last_hop_random,
                    coef, params, tdc, cmhp):
    """Build one stable trace row for any NAMD state count >= 2."""
    coef = np.asarray(coef)
    tdc = np.asarray(tdc, dtype=float)
    cmhp = np.asarray(cmhp, dtype=float)
    params = np.asarray(params, dtype=float)

    def population(index):
        return abs(coef[index]) ** 2 if index < coef.size else np.nan

    def matrix_value(matrix, row, column):
        if row < matrix.shape[0] and column < matrix.shape[1]:
            return matrix[row, column]
        return np.nan

    return {
        "step": istep,
        "t_fs": istep * dt_fs,
        "kernel_called": int(np.isfinite(last_hop_random)),
        "active": active,
        "hopped": int(bool(hopped)),
        "random": params[3] if params.size > 3 else np.nan,
        "pop_1": population(0),
        "pop_2": population(1),
        "pop_3": population(2),
        "tdc_12_au": matrix_value(tdc, 0, 1),
        "tdc_13_au": matrix_value(tdc, 0, 2),
        "tdc_23_au": matrix_value(tdc, 1, 2),
        "p_31": matrix_value(cmhp, 2, 0),
        "p_32": matrix_value(cmhp, 2, 1),
    }


def install_trace(
    output: Path,
    skip_first_hop: bool,
    random_values: np.ndarray | None = None,
    namd_module=None,
) -> SequenceRNG | None:
    if namd_module is None:
        from oqp.library import namd as namd_module

    original_prepare = namd_module.NAMD._prepare_hop_step
    sequence_rng = None if random_values is None else SequenceRNG(random_values)

    def traced_prepare(self, istep):
        enabled = original_prepare(self, istep)
        if skip_first_hop and int(istep) == 1:
            enabled = False
        if sequence_rng is not None:
            self._hop_random_override = sequence_rng
        return enabled

    def record(self, istep, hopped):
        if istep == 0:
            return
        nstate = np.asarray(self.coef).size
        try:
            params = np.asarray(self.mol.data["OQP::namd_params"], dtype=float)
        except (AttributeError, KeyError):
            params = np.empty(0)
        try:
            results = np.asarray(self.mol.data["OQP::namd_results"], dtype=float)
        except (AttributeError, KeyError):
            results = np.empty(0)
        try:
            overlap = np.asarray(
                self.mol.data["OQP::td_states_overlap"], dtype=float
            ).reshape(nstate, nstate).T
            tdc = self._compute_tdc(overlap)
        except (AttributeError, KeyError, ValueError):
            # SOC variants assemble their propagator from manifold-specific
            # overlaps, so a single same-spin overlap matrix may not exist.
            tdc = np.full((nstate, nstate), np.nan)
        cmhp = np.full((nstate, nstate), np.nan)
        if results.size >= nstate * nstate:
            # The native record is a Fortran column-major flattening.
            cmhp = results[: nstate * nstate].reshape(nstate, nstate, order="F")
        row = build_trace_row(
            istep=istep, dt_fs=self.dt_fs, active=self.active,
            hopped=hopped, last_hop_random=self._last_hop_random,
            coef=self.coef, params=params, tdc=tdc, cmhp=cmhp,
        )
        write_header = not output.exists()
        with output.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=row.keys())
            if write_header:
                writer.writeheader()
            writer.writerow(row)

    def traced_logger(original):
        signature = inspect.signature(original)

        def traced(self, *args, **kwargs):
            bound = signature.bind(self, *args, **kwargs)
            bound.apply_defaults()
            result = original(self, *args, **kwargs)
            record(
                self,
                int(bound.arguments["istep"]),
                bool(bound.arguments.get("hopped", False)),
            )
            return result

        return traced

    namd_module.NAMD._prepare_hop_step = traced_prepare
    logger_specs = (
        ("NAMD", "_log_step"),
        ("NAMD_QMMM", "_log_qmmm"),
        ("NAMD_SOC", "_log_soc"),
        ("NAMD_SOC_MCH", "_log_mch"),
        ("NAMD_SOC_QMMM", "_log_soc_qmmm"),
        ("NAMD_SOC_MCH_QMMM", "_log_mch_qmmm"),
    )
    for class_name, method_name in logger_specs:
        cls = getattr(namd_module, class_name, None)
        if cls is None or method_name not in cls.__dict__:
            continue
        setattr(cls, method_name, traced_logger(cls.__dict__[method_name]))
    return sequence_rng


def main() -> None:
    from oqp.pyoqp import Runner

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="thymine.inp")
    parser.add_argument("--trace", default="hop_trace.csv")
    parser.add_argument("--skip-first-hop", action="store_true")
    parser.add_argument(
        "--random-file",
        type=Path,
        help="one prescribed uniform random value per effective hop interval",
    )
    args = parser.parse_args()
    trace = Path(args.trace)
    trace.unlink(missing_ok=True)
    random_values = None
    if args.random_file is not None:
        random_values = np.atleast_1d(np.loadtxt(args.random_file, dtype=float))
        if random_values.size == 0 or np.any((random_values < 0.0) | (random_values >= 1.0)):
            raise ValueError("Prescribed random values must lie in [0, 1)")
    sequence_rng = install_trace(trace, args.skip_first_hop, random_values)
    Runner(project="thymine", input_file=args.input, log="thymine.log", silent=1, usempi=False).run()
    if sequence_rng is not None and sequence_rng.index != sequence_rng.values.size:
        raise RuntimeError(
            "Prescribed NAMD hop RNG consumed "
            f"{sequence_rng.index} of {sequence_rng.values.size} values"
        )


if __name__ == "__main__":
    main()
