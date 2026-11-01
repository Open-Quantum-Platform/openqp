"""Zero-scaled J/K digestion in the replicated and low-memory atomic Fock paths.

int2_{rhf,urohf}_data_t_update and *_consume_shell skip the exchange updates
when scale_exchange == 0 (pure LDA/GGA/meta-GGA) and the Coulomb updates when
scale_coulomb == 0 (the short-range exchange pass of a range-separated
functional).  Both the per-thread replicated accumulator and the shared atomic
accumulator (OQP_FOCK_ATOMIC=1, the automatic low-memory mode) carry the
guard, so each functional class below is run in fresh processes with the
atomic accumulator off and on, serially and with four OpenMP threads.

* The run must enter the requested accumulator: the atomic run prints the
  "shared+atomic" banner and the replicated run does not.
* Every accumulator/thread combination must reproduce the serial replicated
  energy and final alpha/beta Fock matrices to 1e-9.
* The serial replicated energy is pinned to main before the guard, so a
  dropped or doubled J or K term fails.
* Every run must converge with its first-stage converger.  A Fock comparison
  is only a measure of accumulator rounding when the SCF solution is
  well conditioned: the former UKS case, OH (2Pi), has a twofold degenerate
  pi hole whose rotation about the bond axis is a near-null orbital mode,
  lifted only by the DFT grid.  DIIS stalled along it, SOSCF finished the
  run, and the final position along that mode -- and with it the a'/a''
  block of the AO Fock matrix, up to 5e-8 serially and 2e-6 with four
  threads -- was set by rounding, while the energy agreed to 1e-13.  The
  UKS case is now NH2 (2B1, nondegenerate SOMO).

Every run is repeated with OQP_INT2_LAYOUT=legacy (quartets digested by
*_update) and OQP_INT2_LAYOUT=shell (every quartet digested by
*_consume_shell).  That the shell override selects shell blocks for these
consumers is asserted by the production layout resolution in
tests/fortran/int2_shell_block_selftest.F90 (test_int2_shell_blocks.py),
which also pins both consume_shell guards on non-symmetric data by
linearity, F(cJ,cK) = F(cJ,0) + F(0,cK).
"""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

CHILD = r'''
import json, sys
import numpy as np
from oqp.pyoqp import Runner
system, scf, mult, functional = sys.argv[1:5]
config = {
    'input': {'system': system.replace(';', '\n'), 'method': 'hf',
              'basis': '6-31g*', 'runtype': 'energy', 'functional': functional},
    'guess': {'type': 'huckel', 'save_mol': 'False'},
    'scf': {'type': scf, 'multiplicity': mult, 'conv': '1.0e-9',
            'maxit': '100', 'save_molden': 'False'},
    'tests': {'exception': 'True'},
}
runner = Runner(project='fock', input_file=None, log='fock.log',
                input_dict=config, silent=1, usempi=False)
runner.run()
data = runner.mol.data
out = {'energy': float(runner.mol.energies[0]),
       'fock_a': np.asarray(data['OQP::FOCK_A'], dtype=float).ravel().tolist()}
if scf != 'rhf':
    out['fock_b'] = np.asarray(data['OQP::FOCK_B'], dtype=float).ravel().tolist()
print('RESULT ' + json.dumps(out))
'''

WATER = 'O 0 0 -0.0410616; H -0.5331943 0.5331943 -0.6144692; H 0.5331943 -0.5331943 -0.6144692'
NH2 = 'N 0 0 0; H 0 0.79 0.62; H 0 -0.83 0.58'

# Serial replicated energies from main without the zero-scale guard: the RKS
# cases from 84f324bf8, managed build on zeus (macOS x86_64, GCC 15,
# Accelerate ILP64); NH2 from 6808c4bb6, managed build on Ultra (macOS arm64,
# GCC 15, Accelerate ILP64), which reproduces the three RKS pins to 3e-14.
CASES = {
    'rks-pbe': (WATER, 'rhf', '1', 'pbe', -76.32090487015014),
    'uks-pbe': (NH2, 'uhf', '2', 'pbe', -55.794957862486406),
    'rks-cam-b3lyp': (WATER, 'rhf', '1', 'cam-b3lyp', -76.37952876256894),
    # Hybrid control: both J and K are digested (no guard taken).
    'rks-b3lyp': (WATER, 'rhf', '1', 'b3lypv1r', -76.40808725453417),
}


def _run(tmp_path, case, layout, atomic, threads):
    system, scf, mult, functional, _ = CASES[case]
    work = tmp_path / f'{case}-{layout}-a{atomic}-t{threads}'
    work.mkdir()
    env = {k: v for k, v in os.environ.items() if not k.startswith('OQP_')}
    env.update(OMP_NUM_THREADS=str(threads), OQP_FOCK_ATOMIC=str(atomic),
               OQP_INT2_LAYOUT=layout,
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1')
    proc = subprocess.run([sys.executable, '-c', CHILD, system, scf, mult, functional],
                          cwd=work, env=env, capture_output=True, text=True,
                          timeout=900)
    log = work / 'fock.log'
    text = proc.stdout + proc.stderr
    if log.exists():
        text += log.read_text(errors='replace')
    assert proc.returncode == 0, text[-4000:]
    banner = 'Fock accumulator: shared+atomic' in text
    assert banner == bool(atomic), f'{case} {layout} atomic={atomic} threads={threads}: wrong accumulator'
    assert 'Converger stalled' not in text and 'escalating' not in text, \
        f'{case} {layout} atomic={atomic} threads={threads}: SCF needed the escalation ladder'
    line = [x for x in proc.stdout.splitlines() if x.startswith('RESULT ')][-1]
    return json.loads(line[len('RESULT '):])


@pytest.mark.parametrize('layout', ['legacy', 'shell'])
@pytest.mark.parametrize('case', sorted(CASES))
def test_zero_scale_guard_matches_across_accumulators_and_threads(tmp_path, case, layout):
    reference = _run(tmp_path, case, layout, 0, 1)
    assert np.isfinite(reference['energy'])
    pinned = CASES[case][4]
    if pinned is not None:
        assert abs(reference['energy'] - pinned) <= 1e-8, (reference['energy'], pinned)
    # The atomic accumulator rounds differently from the replicated one (it
    # stores each product before the add, the replicated path may contract it
    # into an FMA) and four threads change the summation order.  For these
    # well-conditioned SCF solutions both perturb the final Fock matrices at
    # the 1e-14 level, so they are compared to 1e-9 in every combination.
    for atomic, threads in [(1, 1), (0, 4), (1, 4)]:
        result = _run(tmp_path, case, layout, atomic, threads)
        assert abs(result['energy'] - reference['energy']) <= 1e-9, (layout, atomic, threads)
        for key in ('fock_a', 'fock_b'):
            if key in reference:
                diff = np.max(np.abs(np.asarray(result[key]) - np.asarray(reference[key])))
                assert diff <= 1e-9, (key, layout, atomic, threads, diff)
