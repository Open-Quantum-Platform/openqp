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
from pathlib import Path

import numpy as np

from oqp.library import namd as namd_module
from oqp.pyoqp import Runner


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


def install_trace(
    output: Path,
    skip_first_hop: bool,
    random_values: np.ndarray | None = None,
) -> SequenceRNG | None:
    original_prepare = namd_module.NAMD._prepare_hop_step
    original_log = namd_module.NAMD._log_step
    sequence_rng = None if random_values is None else SequenceRNG(random_values)

    def traced_prepare(self, istep):
        enabled = original_prepare(self, istep)
        if skip_first_hop and int(istep) == 1:
            enabled = False
        if sequence_rng is not None:
            self._hop_random_override = sequence_rng
        return enabled

    def traced_log(
        self,
        istep,
        r,
        hopped=False,
        transition_energy_jump=np.nan,
    ):
        original_failure = None
        try:
            original_log(
                self,
                istep,
                r,
                hopped=hopped,
                transition_energy_jump=transition_energy_jump,
            )
        except Exception as error:
            original_failure = (error, error.__traceback__)
        if istep == 0:
            if original_failure is not None:
                error, traceback = original_failure
                raise error.with_traceback(traceback)
            return
        try:
            nstate = self.nstate
            try:
                params = np.asarray(self.mol.data["OQP::namd_params"], dtype=float)
            except (AttributeError, KeyError):
                params = np.empty(0)
            try:
                results = np.asarray(self.mol.data["OQP::namd_results"], dtype=float)
            except (AttributeError, KeyError):
                results = np.empty(0)
            overlap = np.asarray(
                self.mol.data["OQP::td_states_overlap"], dtype=float
            ).reshape(nstate, nstate).T
            tdc = self._compute_tdc(overlap)
            cmhp = np.full((nstate, nstate), np.nan)
            if results.size >= nstate * nstate:
                # The native record is a Fortran column-major flattening.
                cmhp = results[: nstate * nstate].reshape(
                    nstate, nstate, order="F"
                )
            row = {
                "step": istep,
                "t_fs": istep * self.dt_fs,
                "kernel_called": int(np.isfinite(self._last_hop_random)),
                "active": self.active,
                "hopped": int(bool(hopped)),
                "random": params[3] if params.size > 3 else np.nan,
            }
            for index, coefficient in enumerate(self.coef):
                row[f"pop_{index + 1}"] = abs(coefficient) ** 2
            for left in range(nstate):
                for right in range(left + 1, nstate):
                    row[f"tdc_{left + 1}{right + 1}_au"] = tdc[left, right]
            for source in range(nstate):
                for target in range(nstate):
                    if source != target:
                        row[f"p_{source + 1}{target + 1}"] = cmhp[source, target]
            write_header = not output.exists()
            with output.open("a", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=row.keys())
                if write_header:
                    writer.writeheader()
                writer.writerow(row)
        except Exception:
            if original_failure is not None:
                error, traceback = original_failure
                raise error.with_traceback(traceback) from None
            raise
        if original_failure is not None:
            error, traceback = original_failure
            raise error.with_traceback(traceback)

    namd_module.NAMD._prepare_hop_step = traced_prepare
    namd_module.NAMD._log_step = traced_log
    return sequence_rng


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--trace", default="hop_trace.csv")
    parser.add_argument("--skip-first-hop", action="store_true")
    parser.add_argument(
        "--random-file",
        type=Path,
        help="one prescribed uniform random value per effective hop interval",
    )
    args = parser.parse_args()
    random_values = None
    if args.random_file is not None:
        random_values = np.atleast_1d(np.loadtxt(args.random_file, dtype=float))
        if random_values.size == 0 or np.any((random_values < 0.0) | (random_values >= 1.0)):
            raise ValueError("Prescribed random values must lie in [0, 1)")
    project = args.input.stem
    runner = Runner(
        project=project,
        input_file=str(args.input),
        log=f"{project}.log",
        silent=1,
        usempi=False,
    )
    config = runner.mol.config
    def enabled(value):
        return (value is True) or str(value).lower() in {
            "true", "1", "on", "yes",
        }
    if str(config.get("input", {}).get("runtype", "")).strip().lower() != "namd":
        parser.error("hop tracing requires a NAMD input (input.runtype=namd)")
    if enabled(config.get("input", {}).get("qmmm_flag", False)):
        parser.error("hop tracing currently supports gas-phase NAMD only; QM/MM is rejected")
    if enabled(config.get("md", {}).get("soc", False)):
        parser.error("hop tracing currently supports same-spin NAMD only; SOC is rejected")
    trace = Path(args.trace)
    trace.unlink(missing_ok=True)
    sequence_rng = install_trace(trace, args.skip_first_hop, random_values)
    runner.run()
    if sequence_rng is not None and sequence_rng.index != sequence_rng.values.size:
        raise RuntimeError(
            "Prescribed NAMD hop RNG consumed "
            f"{sequence_rng.index} of {sequence_rng.values.size} values"
        )


if __name__ == "__main__":
    main()
