#!/usr/bin/env python3
"""Gate the resident exact-tlf MRSF metric against two formula oracles.

This diagnostic is deliberately outside the production NAC Python path.  It
runs one H2O energy calculation, asks resident Fortran to export
``OQP::nac_gamma_tlf``, and compares that record against:

1. ``nac_formula_kernel.gamma_closed``: the cofactor-sensitivity oracle;
2. a fourth-order generator derivative of the literal exact state-overlap
   replica for every independent antisymmetric orbital generator.

The gate compares every ordered state pair I!=J independently.  It includes
the same-space dd, ss, and vv generators as well as ds, dv, and sv.  It never
replaces (J,I) from (I,J), averages the two state directions, or projects the
resident matrix onto an antisymmetric form.  Orbital antisymmetry is instead a
separate raw-record invariant that must pass on its own.

For K[p,q]=+1 and K[q,p]=-1 (p>q), the resident slot convention is

    d S_IJ / d theta_pq = 2 gamma_IJ[p,q].

Default H2O invocation::

    python3 tools/nac_lagrangian/gamma_gate.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import nac_formula_kernel as FORMULA


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE / "H2O_energy_tlf0_tight_analytic.inp"


def ordered_state_pairs(nstate):
    """Return all I!=J pairs without inferring either state direction."""
    return [(istate, jstate) for jstate in range(nstate)
            for istate in range(nstate) if istate != jstate]


def decode_resident_gamma(raw, nstate, nbf):
    """Decode TagArray (nbf**2,nstate,nstate) in its Fortran layout."""
    flat = np.asarray(raw).reshape(-1)
    expected = nbf*nbf*nstate*nstate
    if flat.size != expected:
        raise RuntimeError(
            f"OQP::nac_gamma_tlf has {flat.size} values; expected {expected}"
        )
    gamma = np.empty((nstate, nstate, nbf, nbf))
    block = nbf*nbf
    for istate in range(nstate):
        for jstate in range(nstate):
            start = (istate + jstate*nstate)*block
            gamma[istate, jstate] = flat[start:start + block].reshape(
                nbf, nbf, order="F"
            )
    if not np.all(np.isfinite(gamma)):
        raise RuntimeError("OQP::nac_gamma_tlf contains a non-finite value")
    return gamma


def generator_blocks(nocb, noca, nbf):
    """Return all six independent ROHF orbital-generator blocks."""
    if not 0 <= nocb <= noca <= nbf:
        raise ValueError("orbital spaces must satisfy 0 <= nocb <= noca <= nbf")
    spaces = {
        "d": range(0, nocb),
        "s": range(nocb, noca),
        "v": range(noca, nbf),
    }
    blocks = {}
    for low_name, high_name in (
        ("d", "d"), ("d", "s"), ("d", "v"),
        ("s", "s"), ("s", "v"), ("v", "v"),
    ):
        indices = []
        for q in spaces[low_name]:
            for p in spaces[high_name]:
                if p > q:
                    indices.append((p, q))
        blocks[low_name + high_name] = tuple(indices)
    return blocks


def _max_values(values):
    values = list(values)
    return float(np.max(values, initial=0.0))


def analyze_gamma(resident, cofactor, generator_derivative, nocb, noca):
    """Compare raw tensors without orbital or state-space projection."""
    resident = np.asarray(resident)
    cofactor = np.asarray(cofactor)
    generator_derivative = np.asarray(generator_derivative)
    if resident.ndim != 4 or resident.shape[0] != resident.shape[1]:
        raise ValueError("resident gamma must have shape (nstate,nstate,nbf,nbf)")
    if resident.shape != cofactor.shape or resident.shape != generator_derivative.shape:
        raise ValueError("resident, cofactor, and generator tensors must have one shape")
    if resident.shape[2] != resident.shape[3]:
        raise ValueError("orbital gamma dimensions must be square")

    nstate, _, nbf, _ = resident.shape
    pairs = ordered_state_pairs(nstate)
    blocks = generator_blocks(nocb, noca, nbf)
    npair, nblock = len(pairs), len(blocks)
    cofactor_pair = np.zeros(npair)
    resident_generator_pair = np.zeros(npair)
    cofactor_generator_pair = np.zeros(npair)
    orbital_antisym_pair = np.zeros(npair)
    cofactor_block = np.zeros((npair, nblock))
    resident_generator_block = np.zeros((npair, nblock))
    cofactor_generator_block = np.zeros((npair, nblock))

    for pair_index, (istate, jstate) in enumerate(pairs):
        raw = resident[istate, jstate]
        oracle = cofactor[istate, jstate]
        derivative = generator_derivative[istate, jstate]
        cofactor_pair[pair_index] = _max_values(np.abs(raw - oracle).flat)
        orbital_antisym_pair[pair_index] = _max_values(
            np.abs(raw + raw.T).flat
        )

        resident_generator_errors = []
        cofactor_generator_errors = []
        for block_index, indices in enumerate(blocks.values()):
            block_cofactor = []
            block_resident_generator = []
            block_cofactor_generator = []
            for p, q in indices:
                block_cofactor.extend((
                    abs(raw[p, q] - oracle[p, q]),
                    abs(raw[q, p] - oracle[q, p]),
                ))
                block_resident_generator.append(
                    abs(2.0*raw[p, q] - derivative[p, q])
                )
                block_cofactor_generator.append(
                    abs(2.0*oracle[p, q] - derivative[p, q])
                )
            cofactor_block[pair_index, block_index] = _max_values(block_cofactor)
            resident_generator_block[pair_index, block_index] = _max_values(
                block_resident_generator
            )
            cofactor_generator_block[pair_index, block_index] = _max_values(
                block_cofactor_generator
            )
            resident_generator_errors.extend(block_resident_generator)
            cofactor_generator_errors.extend(block_cofactor_generator)

        resident_generator_pair[pair_index] = _max_values(
            resident_generator_errors
        )
        cofactor_generator_pair[pair_index] = _max_values(
            cofactor_generator_errors
        )

    return {
        "ordered_pairs": np.asarray(pairs, dtype=int),
        "block_names": np.asarray(tuple(blocks)),
        "cofactor_pair_max_abs": cofactor_pair,
        "resident_generator_pair_max_abs": resident_generator_pair,
        "cofactor_generator_pair_max_abs": cofactor_generator_pair,
        "orbital_antisym_pair_max_abs": orbital_antisym_pair,
        "cofactor_block_max_abs": cofactor_block,
        "resident_generator_block_max_abs": resident_generator_block,
        "cofactor_generator_block_max_abs": cofactor_generator_block,
    }


def run_gate(input_file=DEFAULT_INPUT, output=None, cofactor_atol=5.0e-11,
             generator_atol=1.0e-9, orbital_antisym_atol=5.0e-13,
             generator_step=1.0e-4, log=None):
    import oqp
    from oqp.pyoqp import Runner

    input_path = Path(input_file).resolve()
    log_path = Path(log).resolve() if log else input_path.with_name(
        input_path.stem + "_gamma_gate.log"
    )
    runner = Runner(input_file=str(input_path), log=str(log_path))
    runner.run()
    mol = runner.mol
    context = FORMULA.build_context(mol)
    nstate = context["nstate"]
    nbf = context["nbf"]
    noca = context["noca"]
    nocb = context["nocb"]
    if nstate < 2:
        raise RuntimeError("the resident gamma gate needs at least two states")
    if noca - nocb != 2:
        raise RuntimeError("the resident gamma gate requires a two-SOMO MRSF reference")

    oqp.mrsf_nac_metric_data(mol)
    resident = decode_resident_gamma(
        mol.data["OQP::nac_gamma_tlf"], nstate, nbf
    )
    streamed = np.zeros_like(resident)
    for jstate in range(nstate):
        oqp.mrsf_nac_metric_column(mol, jstate + 1)
        flat_column = np.asarray(
            mol.data["OQP::nac_gamma_column"]
        ).reshape(-1)
        expected = nbf*nbf*nstate
        if flat_column.size != expected:
            raise RuntimeError(
                "OQP::nac_gamma_column has an inconsistent size"
            )
        for istate in range(nstate):
            start = istate*nbf*nbf
            streamed[istate, jstate] = flat_column[
                start:start + nbf*nbf
            ].reshape(nbf, nbf, order="F")
    streamed_max = float(np.max(np.abs(streamed - resident)))
    cofactor = FORMULA.gamma_closed(context)

    def progress(done, total, p, q):
        if done == 1 or done == total or done % 25 == 0:
            print(
                f"exact-overlap generator sweep {done}/{total} "
                f"(p={p + 1}, q={q + 1})",
                flush=True,
            )

    generator_derivative = FORMULA.generator_derivative_sweep(
        context, step=generator_step, progress=progress
    )
    report = analyze_gamma(
        resident, cofactor, generator_derivative, nocb=nocb, noca=noca
    )

    print("\n===== resident exact-tlf gamma gate =====")
    print("state directions are compared independently; no state antisymmetry is imposed")
    print("orbital blocks: dd ds dv ss sv vv (same-space dd/ss/vv included)")
    print(f"streamed-column/all-pair maxabs={streamed_max:.3e}")
    for pair_index, pair in enumerate(report["ordered_pairs"]):
        istate, jstate = pair + 1
        print(
            f"({istate},{jstate}): "
            f"resident/cofactor={report['cofactor_pair_max_abs'][pair_index]:.3e} "
            f"resident/dS={report['resident_generator_pair_max_abs'][pair_index]:.3e} "
            f"cofactor/dS={report['cofactor_generator_pair_max_abs'][pair_index]:.3e} "
            f"orbital-antisym={report['orbital_antisym_pair_max_abs'][pair_index]:.3e}"
        )
        for block_index, block_name in enumerate(report["block_names"]):
            print(
                f"  {block_name}: resident/cofactor="
                f"{report['cofactor_block_max_abs'][pair_index, block_index]:.3e} "
                f"resident/dS="
                f"{report['resident_generator_block_max_abs'][pair_index, block_index]:.3e} "
                f"cofactor/dS="
                f"{report['cofactor_generator_block_max_abs'][pair_index, block_index]:.3e}"
            )

    output_path = Path(output).resolve() if output else input_path.with_name(
        input_path.stem + "_gamma_gate.npz"
    )
    np.savez(
        output_path,
        resident_gamma=resident,
        streamed_gamma=streamed,
        streamed_column_max_abs=streamed_max,
        cofactor_gamma=cofactor,
        exact_generator_derivative=generator_derivative,
        cofactor_atol=float(cofactor_atol),
        generator_atol=float(generator_atol),
        orbital_antisym_atol=float(orbital_antisym_atol),
        generator_step=float(generator_step),
        input_path=str(input_path),
        **report,
    )
    print(f"saved {output_path}")

    cofactor_max = float(np.max(report["cofactor_pair_max_abs"]))
    resident_generator_max = float(np.max(
        report["resident_generator_pair_max_abs"]
    ))
    cofactor_generator_max = float(np.max(
        report["cofactor_generator_pair_max_abs"]
    ))
    orbital_antisym_max = float(np.max(
        report["orbital_antisym_pair_max_abs"]
    ))
    failures = []
    if streamed_max > cofactor_atol:
        failures.append(
            f"streamed/all-pair={streamed_max:.3e}>{cofactor_atol:.3e}"
        )
    if cofactor_max > cofactor_atol:
        failures.append(f"resident/cofactor={cofactor_max:.3e}>{cofactor_atol:.3e}")
    if resident_generator_max > generator_atol:
        failures.append(
            f"resident/dS={resident_generator_max:.3e}>{generator_atol:.3e}"
        )
    if cofactor_generator_max > generator_atol:
        failures.append(
            f"cofactor/dS={cofactor_generator_max:.3e}>{generator_atol:.3e}"
        )
    if orbital_antisym_max > orbital_antisym_atol:
        failures.append(
            f"orbital-antisym={orbital_antisym_max:.3e}>"
            f"{orbital_antisym_atol:.3e}"
        )
    if failures:
        raise SystemExit("FAIL: " + "; ".join(failures))
    print(
        "PASS: every ordered state pair and all six orbital blocks satisfy "
        "the resident/cofactor/generator contracts"
    )
    return output_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", default=str(DEFAULT_INPUT))
    parser.add_argument("--output")
    parser.add_argument("--cofactor-atol", type=float, default=5.0e-11)
    parser.add_argument("--generator-atol", type=float, default=1.0e-9)
    parser.add_argument("--orbital-antisym-atol", type=float, default=5.0e-13)
    parser.add_argument("--generator-step", type=float, default=1.0e-4)
    parser.add_argument("--log")
    args = parser.parse_args()
    run_gate(
        args.input,
        output=args.output,
        cofactor_atol=args.cofactor_atol,
        generator_atol=args.generator_atol,
        orbital_antisym_atol=args.orbital_antisym_atol,
        generator_step=args.generator_step,
        log=args.log,
    )


if __name__ == "__main__":
    main()
