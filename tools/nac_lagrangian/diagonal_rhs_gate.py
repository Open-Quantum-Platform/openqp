"""Diagonal Lee-gradient source gate for the pairwise MRSF NAC Lagrangian.

This is a diagnostic harness, not a production code path.  For every selected
MRSF state ``I`` it compares two representations of the same diagonal orbital
Lagrangian source:

``rhs_lee``
    The standard MRSF-gradient ``sfrorhs`` result exported through
    ``OQP::nac_zvec_rhs``.  In the Lee et al. normalization this is the
    conventional multiplier right-hand side ``-2 Rbar``.

``ell_pair``
    The current pairwise NAC source evaluated at the algebraic diagonal
    continuation ``y = X_I``, ``X_J = X_I``, and ``gamma = 0``.  It is built
    from the resident Fortran ``mrsf_nac_wpair`` and ``mrsf_nac_response``
    kernels and projected with the native ROHF dual map.  With OpenQP's native
    generator ``K(high,low)=+kappa`` this is ``+2 Rbar``.

These two records must therefore be equal and opposite, not equal:

``ell_pair + rhs_lee = 0``.

The sign is also the distinction between the conventional Lee multiplier and
the computational adjoint used by the production NAC path.  OpenQP's forward
response convention is ``H U^R = B^R`` with ``B^R=-r_R``.  Production solves
``H z=ell_pair`` and adds ``z^T B^R=ell_pair^T U^R``; hence its stored ``z`` is
the negative of the conventional Lee multiplier satisfying
``H Z_Lee=-ell_pair``.  The public operation is still one state-pair Z-vector,
not a ``3N`` forward CPHF calculation.

The pair kernels intentionally reject ``I == J`` because a physical diagonal
derivative coupling is zero.  The harness therefore copies ``X_I`` into two
*distinct* TD slots before calling the bilinear engines.  The slot labels carry
no physics here; only their equal amplitude values do.  Original amplitudes,
target state, selected TagArrays, and diagnostic environment variables are
restored before returning.

Run, for example::

    python tools/nac_lagrangian/diagonal_rhs_gate.py input.inp \
        --states 1,2,3 --output diagonal_rhs_gate.npz --atol 1e-6

Python performs only state-slot orchestration and comparison.  All MRSF,
two-electron, SPC, Fock-response, and gradient-source kernels remain resident
in Fortran.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np


def pack_symmetric(matrix):
    """Pack a symmetric matrix in OpenQP's lower-triangle column order."""
    nbf = matrix.shape[0]
    packed = np.empty(nbf * (nbf + 1) // 2)
    index = 0
    for column in range(nbf):
        for row in range(column + 1):
            packed[index] = matrix[row, column]
            index += 1
    return packed


def unpack_symmetric(packed, nbf):
    """Expand an OpenQP packed symmetric matrix."""
    matrix = np.zeros((nbf, nbf))
    index = 0
    for column in range(nbf):
        for row in range(column + 1):
            matrix[row, column] = packed[index]
            matrix[column, row] = packed[index]
            index += 1
    return matrix


def pack_rohf_dual(matrix, nocb, noca):
    """Project an MO orbital gradient into native ROHF cotangent space."""
    nbf = matrix.shape[0]
    packed = []
    for socc in range(nocb, noca):
        for docc in range(nocb):
            packed.append(matrix[socc, docc] - matrix[docc, socc])
    for virt in range(noca, nbf):
        for docc in range(nocb):
            packed.append(matrix[virt, docc] - matrix[docc, virt])
    for virt in range(noca, nbf):
        for socc in range(nocb, noca):
            packed.append(matrix[virt, socc] - matrix[socc, virt])
    return np.asarray(packed)


def diagonal_source_closure(ell_pair, rhs_lee):
    """Return the sign-aware diagonal identity ``ell_pair + rhs_lee``."""
    ell_pair = np.asarray(ell_pair)
    rhs_lee = np.asarray(rhs_lee)
    if ell_pair.shape != rhs_lee.shape:
        raise ValueError("pair and Lee sources must have the same shape")
    return ell_pair + rhs_lee


def duplicate_state_slots(raw, nstate, nij, state, auxiliary):
    """Return ``raw`` with X_state copied into two distinct TD slots."""
    if state == auxiliary:
        raise ValueError("the auxiliary TD slot must differ from the state")
    if not (0 <= state < nstate and 0 <= auxiliary < nstate):
        raise IndexError("TD state or auxiliary slot is out of range")
    flat = np.asarray(raw).reshape(-1).copy()
    if flat.size != nstate * nij:
        raise ValueError("TD amplitude record has an unexpected size")
    amplitude = flat[state * nij:(state + 1) * nij].copy()
    flat[auxiliary * nij:(auxiliary + 1) * nij] = amplitude
    flat[state * nij:(state + 1) * nij] = amplitude
    return flat.reshape(np.asarray(raw).shape), amplitude


def parse_states(text, nstate):
    """Parse a comma-separated one-based state list."""
    if text is None:
        return list(range(nstate))
    states = []
    for field in text.split(","):
        value = int(field.strip())
        if not 1 <= value <= nstate:
            raise ValueError(f"state {value} is outside 1..{nstate}")
        state = value - 1
        if state not in states:
            states.append(state)
    if not states:
        raise ValueError("at least one state is required")
    return states


def _snapshot_optional_records(mol, names):
    snapshot = {}
    for name in names:
        try:
            snapshot[name] = np.array(mol.data[name], copy=True)
        except Exception:
            pass
    return snapshot


def _restore_records(mol, snapshot):
    for name, value in snapshot.items():
        mol.data[name] = value.copy()


def run_gate(input_file, states_text=None, output=None, atol=1.0e-6, log=None):
    import oqp
    from oqp.pyoqp import Runner

    input_path = Path(input_file).resolve()
    log_path = Path(log).resolve() if log else input_path.with_name(
        input_path.stem + "_diagonal_rhs.log"
    )
    runner = Runner(input_file=str(input_path), log=str(log_path))
    runner.run()
    mol = runner.mol

    nstate = int(mol.config["tdhf"]["nstate"])
    if nstate < 2:
        raise RuntimeError("the diagonal source gate needs at least two TD slots")
    states = parse_states(states_text, nstate)

    raw_original = np.array(mol.data["OQP::td_bvec_mo"], copy=True)
    nbf = np.array(mol.data["OQP::VEC_MO_A"], copy=True).shape[0]
    noca = int(np.asarray(mol.data["nelec_A"]).reshape(-1)[0])
    nocb = int(np.asarray(mol.data["nelec_B"]).reshape(-1)[0])
    nij = raw_original.size // nstate
    expected_nij = noca * (nbf - nocb)
    if nij != expected_nij:
        raise RuntimeError(
            f"unexpected MRSF amplitude dimension {nij}; expected {expected_nij}"
        )

    c_mo = np.array(mol.data["OQP::VEC_MO_A"], copy=True).T
    occ_a = np.zeros(nbf)
    occ_b = np.zeros(nbf)
    occ_a[:noca] = 1.0
    occ_b[:nocb] = 1.0
    original_target = int(mol.config["tdhf"].get("target") or 1)
    mutable_records = _snapshot_optional_records(
        mol,
        (
            "OQP::td_p",
            "OQP::WAO",
            "OQP::td_abxc",
            "OQP::td_mrsf_density",
        ),
    )
    env_original = {
        name: os.environ.get(name)
        for name in ("NAC_DUMP_RHS", "NAC_DUMP_PIJ")
    }
    os.environ["NAC_DUMP_RHS"] = "1"
    os.environ["NAC_DUMP_PIJ"] = "1"

    rhs_lee_all = []
    ell_pair_all = []
    auxiliary_all = []
    block_errors = []

    nds = (noca - nocb) * nocb
    ndv = (nbf - noca) * nocb
    blocks = (
        ("doc-socc", slice(0, nds)),
        ("doc-virt", slice(nds, nds + ndv)),
        ("socc-virt", slice(nds + ndv, None)),
    )

    try:
        # Ensure the old interstate CPHF override cannot affect the standard
        # Lee-gradient source exported below.
        oqp.set_mrsf_nac_cphf(mol, 0, 0)
        for state in states:
            auxiliary = (state + 1) % nstate
            duplicated, amplitude = duplicate_state_slots(
                raw_original, nstate, nij, state, auxiliary
            )

            # Standard diagonal MRSF gradient source: -2 Rbar.
            mol.data["OQP::td_bvec_mo"] = raw_original.copy()
            mol.data.set_tdhf_target(state + 1)
            oqp.tdhf_mrsf_z_vector(mol)
            rhs_lee = np.array(
                mol.data["OQP::nac_zvec_rhs"], copy=True
            ).reshape(-1)

            # Algebraic diagonal continuation of the pair source.  Distinct
            # labels bypass the physical I==J NAC guard; equal amplitudes make
            # every resident bilinear exactly its diagonal quadratic form.
            mol.data["OQP::td_bvec_mo"] = duplicated
            mol.data["OQP::nac_ytil"] = amplitude.copy()
            mol.data["OQP::nac_xstate"] = amplitude.copy()
            oqp.mrsf_nac_wpair(mol, state + 1, auxiliary + 1)
            mt_frozen = np.array(
                mol.data["OQP::nac_mt_frozen"], copy=True
            ).reshape(nbf, nbf).T

            oqp.mrsf_nac_esum(mol, state + 1, auxiliary + 1)
            pij_a = np.array(
                mol.data["OQP::dbg_pij_a"], copy=True
            ).reshape(nbf, nbf).T
            pij_b = np.array(
                mol.data["OQP::dbg_pij_b"], copy=True
            ).reshape(nbf, nbf).T

            mol.data["OQP::nac_dm1_a"] = pack_symmetric(pij_a)
            mol.data["OQP::nac_dm1_b"] = pack_symmetric(pij_b)
            oqp.mrsf_nac_response(mol)
            response_a = unpack_symmetric(
                np.array(mol.data["OQP::nac_v1_a"], copy=True).reshape(-1),
                nbf,
            )
            response_b = unpack_symmetric(
                np.array(mol.data["OQP::nac_v1_b"], copy=True).reshape(-1),
                nbf,
            )
            response_mo_a = c_mo.T @ response_a @ c_mo
            response_mo_b = c_mo.T @ response_b @ c_mo
            mt_response = 2.0 * (
                response_mo_a * occ_a[None, :]
                + response_mo_b * occ_b[None, :]
            )
            ell_pair = pack_rohf_dual(
                mt_frozen + mt_response, nocb=nocb, noca=noca
            )

            if ell_pair.shape != rhs_lee.shape:
                raise RuntimeError(
                    "pair and Lee sources have different ROHF dimensions: "
                    f"{ell_pair.size} versus {rhs_lee.size}"
                )
            if not np.all(np.isfinite(ell_pair)) or not np.all(np.isfinite(rhs_lee)):
                raise RuntimeError("non-finite value in a diagonal source")

            closure = diagonal_source_closure(ell_pair, rhs_lee)
            scale = float(np.dot(ell_pair, rhs_lee)) / (
                float(np.dot(rhs_lee, rhs_lee)) + 1.0e-300
            )
            state_blocks = []
            print(f"state {state + 1} (auxiliary slot {auxiliary + 1}):")
            for label, block in blocks:
                block_delta = closure[block]
                block_max = (
                    float(np.max(np.abs(block_delta)))
                    if block_delta.size
                    else 0.0
                )
                state_blocks.append(block_max)
                print(f"  {label:9s} max|ell_pair+rhs_lee| = {block_max:.12e}")
            max_abs = float(np.max(np.abs(closure)))
            rms = float(np.sqrt(np.mean(closure * closure)))
            rel = float(np.linalg.norm(closure)) / (
                float(np.linalg.norm(rhs_lee)) + 1.0e-300
            )
            print(
                f"  all       max={max_abs:.12e} rms={rms:.12e} "
                f"relative-L2={rel:.12e} fitted-scale={scale:.12e}"
            )

            rhs_lee_all.append(rhs_lee)
            ell_pair_all.append(ell_pair)
            auxiliary_all.append(auxiliary + 1)
            block_errors.append(state_blocks)
    finally:
        mol.data["OQP::td_bvec_mo"] = raw_original.copy()
        mol.data.set_tdhf_target(original_target)
        oqp.set_mrsf_nac_cphf(mol, 0, 0)
        _restore_records(mol, mutable_records)
        for name, value in env_original.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value

    rhs_lee_array = np.stack(rhs_lee_all)
    ell_pair_array = np.stack(ell_pair_all)
    closure = diagonal_source_closure(ell_pair_array, rhs_lee_array)
    state_max = np.max(np.abs(closure), axis=1)
    overall_max = float(np.max(state_max))

    if output is None:
        output_path = input_path.with_name(input_path.stem + "_diagonal_rhs.npz")
    else:
        output_path = Path(output).resolve()
    np.savez(
        output_path,
        states=np.asarray(states, dtype=int) + 1,
        auxiliary_slots=np.asarray(auxiliary_all, dtype=int),
        ell_pair=ell_pair_array,
        # Compatibility alias for the first diagnostic artifact.
        rhs_pair=ell_pair_array,
        rhs_lee=rhs_lee_array,
        closure=closure,
        raw_difference=ell_pair_array - rhs_lee_array,
        state_max_abs=state_max,
        block_max_abs=np.asarray(block_errors),
        block_labels=np.asarray([label for label, _ in blocks]),
        input_path=str(input_path),
        scf_conv=float(mol.config["scf"]["conv"]),
        tdhf_conv=float(mol.config["tdhf"]["conv"]),
        atol=float(atol),
    )
    print(f"saved {output_path}")
    if overall_max > atol:
        raise SystemExit(
            f"FAIL: diagonal ell_pair+rhs_lee max error {overall_max:.12e} "
            f"exceeds {atol:.12e}"
        )
    print(
        f"PASS: diagonal ell_pair+rhs_lee max error {overall_max:.12e} "
        f"<= {atol:.12e}"
    )
    return output_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="converged MRSF energy input")
    parser.add_argument(
        "--states",
        help="comma-separated one-based TD state list (default: all)",
    )
    parser.add_argument("--output", help="output NPZ artifact")
    parser.add_argument("--atol", type=float, default=1.0e-6)
    parser.add_argument("--log", help="OpenQP log path")
    args = parser.parse_args()
    run_gate(
        args.input,
        states_text=args.states,
        output=args.output,
        atol=args.atol,
        log=args.log,
    )


if __name__ == "__main__":
    main()
