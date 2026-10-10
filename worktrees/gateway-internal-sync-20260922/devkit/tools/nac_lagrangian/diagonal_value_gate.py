"""Full-value Lee-gradient continuation gate for the MRSF NAC Lagrangian.

For each requested state I, copy X_I into a distinct auxiliary TD slot and
evaluate the pair kernels with y=X_I and gamma=0.  This duplicated-slot
continuation must reproduce the Lee excitation-energy gradient

    d Omega_I/dR = d E_I/dR - d E_0/dR.

The same run also compares the two multiplier conventions.  The standard
gradient path exports xk=2*Zbar, whereas the native pair adjoint stores
zeta=-Zbar, so

    zeta + xk/2 = 0.

Both records are one-dimensional native-ROHF vectors with identical block
ordering: singly-occupied/doubly-occupied (SD), virtual/doubly-occupied (DV),
then virtual/singly-occupied (SV).  No reshape, transpose, or block permutation
is allowed at this comparison boundary.

This is a diagnostic gate only.  All scientific kernels and coordinate loops
remain in resident Fortran; Python selects states, snapshots records, and
compares the returned arrays.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np


def rohf_rotation_slices(nbf, noca, nocb):
    """Return the SD/DV/SV slices in OpenQP's native ROHF vector order."""
    if not 0 <= nocb <= noca <= nbf:
        raise ValueError("ROHF occupations must satisfy 0 <= nocb <= noca <= nbf")
    nsocc = noca - nocb
    nvirt = nbf - noca
    nds = nsocc*nocb
    ndv = nvirt*nocb
    nsv = nvirt*nsocc
    return {
        "socc-docc": slice(0, nds),
        "virt-docc": slice(nds, nds + ndv),
        "virt-socc": slice(nds + ndv, nds + ndv + nsv),
    }


def require_rohf_vector(record, name, nbf, noca, nocb):
    """Validate a flat SD/DV/SV multiplier without hiding shape mistakes."""
    vector = np.asarray(record)
    slices = rohf_rotation_slices(nbf, noca, nocb)
    expected = slices["virt-socc"].stop
    if vector.ndim != 1:
        raise RuntimeError(
            f"{name} must be a one-dimensional ROHF vector; got shape "
            f"{vector.shape}"
        )
    if vector.size != expected:
        raise RuntimeError(
            f"{name} has length {vector.size}; expected {expected} "
            "for the SD/DV/SV packing"
        )
    if not np.all(np.isfinite(vector)):
        raise RuntimeError(f"{name} contains a non-finite multiplier")
    return np.array(vector, copy=True)


def zvector_convention_closure(zeta, xk, nbf, noca, nocb):
    """Return the matched OpenQP identity ``zeta + xk/2``."""
    native = require_rohf_vector(
        zeta, "OQP::nac_rohf_solution", nbf, noca, nocb
    )
    legacy = require_rohf_vector(
        xk, "OQP::nac_zvec_solution", nbf, noca, nocb
    )
    return native + 0.5*legacy


def duplicate_state_slots(raw, nstate, nij, state, auxiliary):
    """Copy one amplitude into two distinct slots and return it separately."""
    if state == auxiliary:
        raise ValueError("the auxiliary TD slot must differ from the state")
    flat = np.asarray(raw).reshape(-1).copy()
    if flat.size != nstate * nij:
        raise ValueError("TD amplitude record has an unexpected size")
    amplitude = flat[state * nij:(state + 1) * nij].copy()
    flat[auxiliary * nij:(auxiliary + 1) * nij] = amplitude
    flat[state * nij:(state + 1) * nij] = amplitude
    return flat.reshape(np.asarray(raw).shape), amplitude


def parse_states(text, nstate):
    if text is None:
        return list(range(nstate))
    states = []
    for field in text.split(","):
        state = int(field.strip()) - 1
        if not 0 <= state < nstate:
            raise ValueError(f"state {state + 1} is outside 1..{nstate}")
        if state not in states:
            states.append(state)
    if not states:
        raise ValueError("at least one state is required")
    return states


def run_gate(input_file, states_text=None, output=None, value_atol=5.0e-6,
             z_atol=5.0e-6, log=None):
    import oqp
    from oqp.pyoqp import Runner

    input_path = Path(input_file).resolve()
    log_path = Path(log).resolve() if log else input_path.with_name(
        input_path.stem + "_diagonal_value.log"
    )
    runner = Runner(input_file=str(input_path), log=str(log_path))
    runner.run()
    mol = runner.mol

    nstate = int(mol.config["tdhf"]["nstate"])
    if nstate < 2:
        raise RuntimeError("the diagonal value gate needs at least two TD slots")
    if int(mol.config["tdhf"]["z_solver"]) not in (2, 3):
        raise RuntimeError("use [tdhf] z_solver=2 or 3 for this matched gate")
    if float(mol.config["tdhf"]["zvconv"]) > 1.0e-10:
        raise RuntimeError("use [tdhf] zvconv<=1e-10 for this matched gate")

    states = parse_states(states_text, nstate)
    natom = int(mol.data["natom"])
    ncoord = 3*natom
    nbf = np.asarray(mol.data["OQP::VEC_MO_A"]).shape[0]
    noca = int(np.asarray(mol.data["nelec_A"]).reshape(-1)[0])
    nocb = int(np.asarray(mol.data["nelec_B"]).reshape(-1)[0])
    nij = noca*(nbf - nocb)
    raw_original = np.array(mol.data["OQP::td_bvec_mo"], copy=True)
    original_target = int(mol.config["tdhf"].get("target") or 1)
    old_dump = os.environ.get("NAC_DUMP_RHS")
    os.environ["NAC_DUMP_RHS"] = "1"

    # Ground-state contribution needed to isolate Lee Eq. (3.21), Omega^R.
    oqp.hf_gradient(mol)
    ground = np.array(mol.get_grad(), copy=True).reshape(-1)

    excitation_values = []
    pair_values = []
    matched_pair_values = []
    pair_components = []
    matched_pair_components = []
    legacy_solutions = []
    native_solutions = []
    z_closures = []
    z_block_errors = []
    z_scales = []
    auxiliaries = []

    try:
        for state in states:
            auxiliary = (state + 1) % nstate

            # Tightly converged standard Lee diagonal gradient and multiplier.
            mol.data["OQP::td_bvec_mo"] = raw_original.copy()
            mol.data.set_tdhf_target(state + 1)
            oqp.tdhf_mrsf_z_vector(mol)
            if not bool(mol.mol_energy.Z_Vector_converged):
                raise RuntimeError(
                    f"legacy Lee Z-vector did not converge for state {state + 1}"
                )
            xk = np.array(
                mol.data["OQP::nac_zvec_solution"], copy=True
            )
            xk = require_rohf_vector(
                xk, "OQP::nac_zvec_solution", nbf, noca, nocb
            )
            oqp.tdhf_mrsf_gradient(mol)
            total = np.array(mol.get_grad(), copy=True).reshape(-1)
            excitation = total - ground

            # Algebraic I=J continuation through a distinct auxiliary slot.
            duplicated, amplitude = duplicate_state_slots(
                raw_original, nstate, nij, state, auxiliary
            )
            mol.data["OQP::td_bvec_mo"] = duplicated
            mol.data["OQP::nac_ytil"] = amplitude.copy()
            mol.data["OQP::nac_xstate"] = amplitude.copy()
            oqp.mrsf_nac_wpair(mol, state + 1, auxiliary + 1)
            oqp.mrsf_nac_amp_pair(mol, state + 1, auxiliary + 1)
            oqp.mrsf_nac_esum(mol, state + 1, auxiliary + 1)
            mol.data["OQP::td_bvec_mo"] = raw_original.copy()
            oqp.mrsf_nac_response(mol)
            mol.data["OQP::nac_gamma_pair"] = np.zeros(nbf*nbf)
            oqp.mrsf_nac_rohf_pair_overlap(mol)
            oqp.mrsf_nac_rohf_zvector(mol)
            zeta = np.array(
                mol.data["OQP::nac_rohf_solution"], copy=True
            )
            zeta = require_rohf_vector(
                zeta, "OQP::nac_rohf_solution", nbf, noca, nocb
            )
            mol.data["OQP::nac_rohf_z"] = zeta.copy()
            oqp.mrsf_nac_rohf_hf_adjoint(mol)
            oqp.mrsf_nac_xc_adjoint(mol)

            amp = np.array(mol.data["OQP::nac_amp"], copy=True).reshape(-1)
            pair_start = (state + auxiliary*nstate)*ncoord
            amp_component = amp[pair_start:pair_start+ncoord].copy()
            esum_component = np.array(
                mol.data["OQP::nac_esum"], copy=True
            ).reshape(-1)
            overlap_component = np.array(
                mol.data["OQP::nac_pair_overlap"], copy=True
            ).reshape(-1)
            z_hf_component = np.array(
                mol.data["OQP::nac_rohf_hf_adjoint"], copy=True
            ).reshape(-1)
            z_xc_component = np.array(
                mol.data["OQP::nac_rohf_xc_adjoint"], copy=True
            ).reshape(-1)
            components = np.stack((
                amp_component,
                esum_component,
                z_hf_component,
                z_xc_component,
                overlap_component,
            ))
            pair = np.sum(components, axis=0)

            # Isolate value errors caused only by the iterative native
            # multiplier.  Re-contract the resident adjoints with the matched
            # Lee solution zeta=-xk/2 while leaving every skeleton/overlap term
            # unchanged.  This is a diagnostic comparison, not a production
            # fallback: the native solver must still meet its own residual and
            # value contracts.
            mol.data["OQP::nac_rohf_z"] = -0.5*xk
            oqp.mrsf_nac_rohf_hf_adjoint(mol)
            oqp.mrsf_nac_xc_adjoint(mol)
            matched_z_hf = np.array(
                mol.data["OQP::nac_rohf_hf_adjoint"], copy=True
            ).reshape(-1)
            matched_z_xc = np.array(
                mol.data["OQP::nac_rohf_xc_adjoint"], copy=True
            ).reshape(-1)
            matched_components = np.stack((
                amp_component,
                esum_component,
                matched_z_hf,
                matched_z_xc,
                overlap_component,
            ))
            matched_pair = np.sum(matched_components, axis=0)

            value_error = float(np.max(np.abs(pair - excitation)))
            matched_value_error = float(
                np.max(np.abs(matched_pair - excitation))
            )
            z_closure = zvector_convention_closure(
                zeta, xk, nbf, noca, nocb
            )
            z_error = float(np.max(np.abs(z_closure)))
            z_slices = rohf_rotation_slices(nbf, noca, nocb)
            z_block_error = np.asarray([
                np.max(np.abs(z_closure[block_slice]), initial=0.0)
                for block_slice in z_slices.values()
            ])
            scale = float(np.dot(pair, excitation)) / (
                float(np.dot(excitation, excitation)) + 1.0e-300
            )
            z_scale = float(np.dot(zeta, xk)) / (
                float(np.dot(xk, xk)) + 1.0e-300
            )
            print(
                f"state {state + 1}: max|pair-dOmega|={value_error:.12e} "
                f"matched-z={matched_value_error:.12e} "
                f"max|zeta+xk/2|={z_error:.12e} scale={scale:.12e} "
                f"zeta/xk={z_scale:.12e}",
                flush=True,
            )

            excitation_values.append(excitation)
            pair_values.append(pair)
            matched_pair_values.append(matched_pair)
            pair_components.append(components)
            matched_pair_components.append(matched_components)
            legacy_solutions.append(xk)
            native_solutions.append(zeta)
            z_closures.append(z_closure)
            z_block_errors.append(z_block_error)
            z_scales.append(z_scale)
            auxiliaries.append(auxiliary + 1)
    finally:
        mol.data["OQP::td_bvec_mo"] = raw_original.copy()
        mol.data.set_tdhf_target(original_target)
        if old_dump is None:
            os.environ.pop("NAC_DUMP_RHS", None)
        else:
            os.environ["NAC_DUMP_RHS"] = old_dump

    excitation_array = np.stack(excitation_values)
    pair_array = np.stack(pair_values)
    matched_pair_array = np.stack(matched_pair_values)
    pair_component_array = np.stack(pair_components)
    matched_pair_component_array = np.stack(matched_pair_components)
    legacy_array = np.stack(legacy_solutions)
    native_array = np.stack(native_solutions)
    z_closure_array = np.stack(z_closures)
    z_block_error_array = np.stack(z_block_errors)
    value_state_max = np.max(np.abs(pair_array - excitation_array), axis=1)
    z_state_max = np.max(np.abs(z_closure_array), axis=1)

    output_path = Path(output).resolve() if output else input_path.with_name(
        input_path.stem + "_diagonal_value.npz"
    )
    np.savez(
        output_path,
        states=np.asarray(states, dtype=int) + 1,
        auxiliary_slots=np.asarray(auxiliaries, dtype=int),
        ground_gradient=ground,
        excitation_gradient=excitation_array,
        pair_continuation=pair_array,
        matched_z_pair_continuation=matched_pair_array,
        pair_component_names=np.asarray((
            "amp", "esum", "z_hf", "z_xc", "overlap"
        )),
        pair_components=pair_component_array,
        matched_z_pair_components=matched_pair_component_array,
        legacy_xk=legacy_array,
        native_zeta=native_array,
        z_closure=z_closure_array,
        z_block_names=np.asarray(tuple(rohf_rotation_slices(
            nbf, noca, nocb
        ).keys())),
        z_block_max_abs=z_block_error_array,
        zeta_over_xk=np.asarray(z_scales),
        value_state_max_abs=value_state_max,
        z_state_max_abs=z_state_max,
        value_atol=float(value_atol),
        z_atol=float(z_atol),
        input_path=str(input_path),
    )
    print(f"saved {output_path}")
    value_max = float(np.max(value_state_max))
    z_max = float(np.max(z_state_max))
    if value_max > value_atol or z_max > z_atol:
        raise SystemExit(
            f"FAIL: value={value_max:.12e}/{value_atol:.12e}, "
            f"z={z_max:.12e}/{z_atol:.12e}"
        )
    print(
        f"PASS: value={value_max:.12e}<={value_atol:.12e}, "
        f"z={z_max:.12e}<={z_atol:.12e}"
    )
    return output_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("--states")
    parser.add_argument("--output")
    parser.add_argument("--value-atol", type=float, default=5.0e-6)
    parser.add_argument("--z-atol", type=float, default=5.0e-6)
    parser.add_argument("--log")
    args = parser.parse_args()
    run_gate(
        args.input,
        states_text=args.states,
        output=args.output,
        value_atol=args.value_atol,
        z_atol=args.z_atol,
        log=args.log,
    )


if __name__ == "__main__":
    main()
