"""v21: re-judge production analytic NAC after the 7.54 XC response fix.

Run:
    python v21_production_gate.py <energy.inp> <numerical-dcv.npz> \
        [out.npz] [max-error]

Davidson state phases are process-random.  The gate resolves one sign per
state, rejects inconsistent independent pair flips, and reports the aligned
error for every unordered state pair.  It also requires matching tight-run
energies and records the exact reference SHA256 in the output artifact.
"""

import os
import hashlib
import sys

import numpy as np


def main(inp, reference_npz, output_npz=None, tolerance=1.0e-3):
    from oqp.library.nac_analytic import analytic_nac
    from oqp.pyoqp import Runner
    from nac_reference_gate import (
        NACGateError,
        compare_payloads,
    )

    runner = Runner(input_file=inp, log=inp.replace('.inp', '_v21.log'))
    runner.run()
    nacv, dcv = analytic_nac(runner.mol)

    with np.load(reference_npz, allow_pickle=False) as frozen:
        reference_payload = {
            key: np.array(frozen[key], copy=True) for key in frozen.files
        }
    reference = reference_payload['dcv']
    nstate = dcv.shape[0]
    energies = np.array(runner.mol.energies, copy=True)
    flags = np.array(['analytic-v3-zvector'] * (3 * runner.mol.data['natom']))
    candidate_payload = {
        'nacv': nacv,
        'dcv': dcv,
        'energies': energies,
        'flags': flags,
    }
    try:
        result = compare_payloads(
            reference_payload,
            candidate_payload,
            component_atol=float('inf'),
            energy_atol=1.0e-8,
            require_flags=True,
            label='production analytic NAC',
        )
    except NACGateError as exc:
        raise SystemExit(f'FAIL: invalid production/reference artifact: {exc}')

    aligned = np.array(dcv, copy=True)
    for istate in range(nstate):
        for jstate in range(nstate):
            aligned[istate, jstate] *= (
                result.state_signs[istate] * result.state_signs[jstate]
            )

    print('===== v21 production analytic NAC gate =====')
    print('state gauge:', ' '.join(f'{sign:+d}' for sign in result.state_signs))
    for metric in result.pair_metrics:
        print(
            f'({metric.istate},{metric.jstate}): '
            f'|num|={metric.reference_norm:.10f} '
            f'|analytic|={metric.candidate_norm:.10f} '
            f'phase={metric.phase:+d} '
            f'maxdiff={metric.max_component_error:.8e} '
            f'|diff|={metric.l2_error:.8e}'
        )
    print(f'energy maxdiff={result.max_energy_error:.8e}')

    translation = np.sum(aligned, axis=2)
    print(
        'raw-electronic NAC atom-sum diagnostic '
        '(not a translation-invariance pass/fail test): '
        f'maxabs={np.max(np.abs(translation)):.8e}'
    )

    if output_npz is None:
        output_npz = inp.replace('.inp', '_v21.npz')
    with open(reference_npz, 'rb') as reference_file:
        reference_sha256 = hashlib.sha256(reference_file.read()).hexdigest()
    np.savez(
        output_npz,
        nacv=nacv,
        dcv=dcv,
        energies=energies,
        flags=flags,
        dcv_phase_aligned=aligned,
        dcv_reference=reference,
        input_path=os.path.abspath(inp),
        reference_path=os.path.abspath(reference_npz),
        reference_sha256=reference_sha256,
        scf_conv=float(runner.mol.config['scf']['conv']),
        tdhf_conv=float(runner.mol.config['tdhf']['conv']),
    )
    print(f'saved {os.path.abspath(output_npz)}')
    if result.max_component_error > tolerance:
        raise SystemExit(
            f'FAIL: max pair/component error {result.max_component_error:.8e} '
            f'exceeds {tolerance:.8e}'
        )
    print(
        f'PASS: max pair/component error {result.max_component_error:.8e} '
        f'<= {tolerance:.8e}'
    )


if __name__ == '__main__':
    main(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3] if len(sys.argv) > 3 else None,
        float(sys.argv[4]) if len(sys.argv) > 4 else 1.0e-3,
    )
