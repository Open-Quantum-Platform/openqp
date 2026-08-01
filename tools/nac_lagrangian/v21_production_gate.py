"""v21: re-judge production analytic NAC after the 7.54 XC response fix.

Run:
    python v21_production_gate.py <energy.inp> <numerical-dcv.npz> \
        [out.npz] [max-error]

Davidson state phases are process-random.  The gate therefore reports both
the raw cosine and the sign-resolved error for every unordered state pair.
"""

import os
import sys

import numpy as np


def main(inp, reference_npz, output_npz=None, tolerance=1.0e-3):
    from oqp.library.nac_analytic import analytic_nac
    from oqp.pyoqp import Runner

    runner = Runner(input_file=inp, log=inp.replace('.inp', '_v21.log'))
    runner.run()
    nacv, dcv = analytic_nac(runner.mol)

    frozen = np.load(reference_npz)
    reference = frozen['dcv' if 'dcv' in frozen.files else frozen.files[0]]
    nstate = dcv.shape[0]
    aligned = np.array(dcv, copy=True)
    largest_error = 0.0

    print('===== v21 production analytic NAC gate =====')
    for istate in range(nstate):
        for jstate in range(istate + 1, nstate):
            pred = dcv[istate, jstate].reshape(-1)
            ref = reference[istate, jstate].reshape(-1)
            denom = np.linalg.norm(pred) * np.linalg.norm(ref) + 1.0e-300
            cosine = float(np.dot(pred, ref) / denom)
            phase = 1.0 if np.linalg.norm(pred - ref) <= np.linalg.norm(pred + ref) else -1.0
            aligned[istate, jstate] *= phase
            aligned[jstate, istate] *= phase
            error = phase * pred - ref
            largest_error = max(largest_error, float(np.max(np.abs(error))))
            print(
                f'({istate + 1},{jstate + 1}): '
                f'|num|={np.linalg.norm(ref):.10f} '
                f'|analytic|={np.linalg.norm(pred):.10f} '
                f'cos={cosine:+.10f} phase={phase:+.0f} '
                f'maxdiff={np.max(np.abs(error)):.8e} '
                f'|diff|={np.linalg.norm(error):.8e}'
            )

    translation = np.sum(aligned, axis=2)
    print(f'translation maxabs={np.max(np.abs(translation)):.8e}')

    if output_npz is None:
        output_npz = inp.replace('.inp', '_v21.npz')
    np.savez(
        output_npz,
        nacv=nacv,
        dcv=dcv,
        dcv_phase_aligned=aligned,
        dcv_reference=reference,
    )
    print(f'saved {os.path.abspath(output_npz)}')
    if largest_error > tolerance:
        raise SystemExit(
            f'FAIL: max pair/component error {largest_error:.8e} '
            f'exceeds {tolerance:.8e}'
        )
    print(
        f'PASS: max pair/component error {largest_error:.8e} '
        f'<= {tolerance:.8e}'
    )


if __name__ == '__main__':
    main(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3] if len(sys.argv) > 3 else None,
        float(sys.argv[4]) if len(sys.argv) > 4 else 1.0e-3,
    )
