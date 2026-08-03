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
import os
from pathlib import Path

import numpy as np


def _namd_sidecar_paths(md: dict, log_path: Path) -> dict[str, Path]:
    """Resolve the production NAMD outputs without constructing the driver."""
    log_path = log_path.expanduser().resolve()
    log_dir = log_path.parent
    stem = log_path.stem

    def output(key: str, suffix: str) -> Path:
        configured = str(md.get(key, '') or '').strip()
        if configured:
            path = Path(configured).expanduser()
            return path.resolve() if path.is_absolute() else (log_dir / path).resolve()
        return (log_dir / f'{stem}{suffix}').resolve()

    return {
        'trajectory-file': output('trajectory_file', '.namd.trj'),
        'restart-file': output('restart_file', '.namd.restart.npz'),
        'restart-manifest-file': (
            log_dir / f'{stem}.namd.restart.oqp').resolve(),
    }


def _paths_alias(left: Path, right: Path) -> bool:
    """Return whether two path spellings designate the same filesystem object."""
    if os.path.normcase(os.path.realpath(left)) == os.path.normcase(
            os.path.realpath(right)):
        return True
    try:
        return os.path.samefile(left, right)
    except (FileNotFoundError, OSError):
        return False


def _reject_trace_aliases(trace: Path, protected_paths: dict[str, Path]) -> None:
    aliases = [
        name for name, path in protected_paths.items()
        if _paths_alias(trace, path)
    ]
    if aliases:
        raise ValueError(
            f"--trace must not alias the {', '.join(aliases)} path")


def validate_trace_output(
    trace: Path,
    *,
    input_path: Path,
    random_file: Path | None,
    log_path: Path,
) -> None:
    """Protect diagnostic inputs before the trace is removed or opened."""
    protected = {'input deck': input_path, 'job log': log_path}
    if random_file is not None:
        protected['random replay file'] = random_file
    _reject_trace_aliases(trace, protected)


def _trace_hop_matrices(driver, kernel_called: bool):
    """Return current overlap/TDC matrices even when the hop was skipped."""
    nstate = getattr(driver, "nstate", None)
    if nstate is None:
        nstate = np.asarray(driver.coef).size
    data = getattr(driver.mol, "data", {})

    def native_tag(name):
        if not kernel_called:
            return None
        try:
            return data[name]
        except (AttributeError, KeyError, TypeError):
            return None

    overlap = native_tag("OQP::namd_stas")
    tdc = native_tag("OQP::namd_tdc")
    if overlap is None:
        overlap = getattr(driver, "_last_state_overlap", None)
    if tdc is None:
        tdc = getattr(driver, "_last_overlap_tdc", None)
    if overlap is None:
        overlap = np.full((nstate, nstate), np.nan)
    if tdc is None:
        tdc = np.full((nstate, nstate), np.nan)
    return (
        np.asarray(overlap, dtype=float).reshape(nstate, nstate),
        np.asarray(tdc, dtype=float).reshape(nstate, nstate),
    )


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
                    coef, params, tdc, cmhp, time_fs=None):
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
        "t_fs": istep * dt_fs if time_fs is None else time_fs,
        "kernel_called": int(np.isfinite(last_hop_random)),
        "active": active,
        "hopped": int(bool(hopped)),
        "random": (last_hop_random if np.isfinite(last_hop_random)
                   else params[3] if params.size > 3 else np.nan),
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

    def record_trace(self, istep, hopped):
        if istep == 0:
            return
        nstate = np.asarray(self.coef).size
        coef = np.asarray([
            self.coef[state] for state in range(nstate)
        ])
        kernel_called = bool(np.isfinite(self._last_hop_random))
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
            overlap, tdc = _trace_hop_matrices(self, kernel_called)
        cmhp = np.full((nstate, nstate), np.nan)
        observed_probabilities = getattr(
            self, "_last_hop_probabilities", None)
        if (observed_probabilities is not None
                and np.asarray(observed_probabilities).size == nstate*nstate):
            cmhp = np.asarray(
                observed_probabilities, dtype=float).reshape(nstate, nstate)
        elif results.size >= nstate * nstate:
            # The native record is a Fortran column-major flattening.
            cmhp = results[: nstate * nstate].reshape(nstate, nstate, order="F")
        row = build_trace_row(
            istep=istep, dt_fs=self.dt_fs, active=self.active,
            hopped=hopped, last_hop_random=self._last_hop_random,
            coef=coef, params=params, tdc=tdc, cmhp=cmhp,
            time_fs=(getattr(self, '_t_fs', istep*self.dt_fs)
                     if getattr(self, 'dt_adaptive', False) else None),
        )
        row.update({
            f"pop_{state + 1}": abs(coef[state]) ** 2
            for state in range(nstate)
        })
        row.update({
            f"overlap_{left + 1}{right + 1}":
                overlap[left, right]
            for left in range(nstate)
            for right in range(nstate)
        })
        row.update({
            f"tdc_{left + 1}{right + 1}_au": tdc[left, right]
            for left in range(nstate)
            for right in range(left + 1, nstate)
        })
        row.update({
            f"p_{source + 1}{target + 1}": cmhp[source, target]
            for source in range(nstate)
            for target in range(nstate)
            if source != target
        })
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
            step = bound.arguments.get("istep")
            if step is None and args:
                step = args[0]
            original_failure = None
            result = None
            try:
                try:
                    result = original(self, *args, **kwargs)
                except Exception as error:
                    original_failure = (error, error.__traceback__)
            finally:
                try:
                    record_trace(
                        self,
                        int(step),
                        bool(bound.arguments.get(
                            "hopped", kwargs.get("hopped", False))),
                    )
                except Exception:
                    if original_failure is not None:
                        error, traceback = original_failure
                        raise error.with_traceback(traceback) from None
                    raise
            if original_failure is not None:
                error, traceback = original_failure
                raise error.with_traceback(traceback)
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
    trace = Path(args.trace).expanduser().resolve()
    input_path = args.input.expanduser().resolve()
    random_path = (None if args.random_file is None
                   else args.random_file.expanduser().resolve())
    project = input_path.stem
    log_path = input_path.with_suffix('.log')
    validate_trace_output(
        trace, input_path=input_path, random_file=random_path,
        log_path=log_path)
    from oqp.pyoqp import Runner

    runner = Runner(
        project=project, input_file=str(input_path), log=str(log_path),
        silent=1, usempi=False)
    protected_paths = _namd_sidecar_paths(
        runner.mol.config.get('md', {}), log_path)
    _reject_trace_aliases(trace, protected_paths)
    runtype = str(
        runner.mol.config.get('input', {}).get('runtype', '')
    ).strip().lower()
    if runtype != 'namd':
        raise ValueError('hop tracing requires a NAMD input')
    soc = runner.mol.config.get("md", {}).get("soc", False)
    if soc:
        raise ValueError(
            "trace_namd_hop does not support SOC NAMD logging paths"
        )
    trace.unlink(missing_ok=True)
    random_values = None
    if random_path is not None:
        random_values = np.atleast_1d(np.loadtxt(random_path, dtype=float))
        if random_values.size == 0 or np.any((random_values < 0.0) | (random_values >= 1.0)):
            raise ValueError("Prescribed random values must lie in [0, 1)")
    sequence_rng = install_trace(trace, args.skip_first_hop, random_values)
    runner.run()
    if sequence_rng is not None and sequence_rng.index != sequence_rng.values.size:
        raise RuntimeError(
            "Prescribed NAMD hop RNG consumed "
            f"{sequence_rng.index} of {sequence_rng.values.size} values"
        )


if __name__ == "__main__":
    main()
