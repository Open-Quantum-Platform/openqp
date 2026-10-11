"""The default TRAH micro-solver (augmented-Hessian Davidson) must return a
Newton step at every gradient size.

It used to stop at a fixed absolute residual of 1e-6.  Once |g| fell below that
value the gradient-seeded start vector already passed the test, so every later
step was the diagonal-preconditioned gradient step.  That step contracts only
linearly; where diag(H) is a poor model of H (the DFT orbital Hessian) the SCF
stalled near |g| ~ 1e-8 and never reached conv=1e-11.  The model problem below
has such a Hessian and must converge in a few steps.
"""
import os
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_augmented_hessian_step_stays_newton_at_small_gradient(tmp_path):
    compiler = os.environ.get('FC') or shutil.which('gfortran-15') or shutil.which('gfortran')
    if not compiler:
        pytest.skip('Fortran compiler required for the native TRAH acceptance test')
    exe = tmp_path / 'trah_ah_newton'
    subprocess.run([
        compiler, '-fcheck=all', '-fdefault-integer-8', '-ffree-line-length-none',
        str(ROOT / 'tests/fixtures/trah_core_support.f90'),
        str(ROOT / 'source/trah_core.F90'),
        str(ROOT / 'tests/fortran/trah_ah_newton_selftest.F90'),
        '-o', str(exe),
    ], cwd=tmp_path, check=True, capture_output=True, text=True)
    result = subprocess.run([str(exe)], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'PASS: augmented-Hessian Newton' in result.stdout


CH2O_DECK = """[input]
system=
   C 0.000000 0.000000 -0.528500
   O 0.000000 0.000000  0.677500
   H 0.000000 0.935100 -1.112000
   H 0.000000 -0.935100 -1.112000
charge=0
runtype=energy
basis=6-31g*
functional=blyp
method=hf

[guess]
type=huckel

[scf]
type={scftype}
multiplicity={mult}
maxit=100
conv=1e-11
converger_type=trah
trh_impl=native
escalation=trah
"""


def _runtime_available():
    try:
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


@pytest.mark.skipif(not _runtime_available(), reason='compiled OpenQP runtime unavailable')
@pytest.mark.parametrize('scftype,mult', [('rhf', 1), ('rohf', 3), ('uhf', 3)])
def test_trah_converges_dft_to_tight_gradient(tmp_path, scftype, mult):
    """CH2O/BLYP/6-31G*: with the absolute Davidson stop TRAH reached
    |g| = 2.3e-8 (RKS), 1.2e-8 (ROKS triplet) and 3.6e-8 (UKS triplet) and
    gave up after 22-32 refinement steps; Newton steps reach 1e-11 in a couple
    of steps."""
    import re
    import sys
    deck = tmp_path / 'ch2o_trah.inp'
    deck.write_text(CH2O_DECK.format(scftype=scftype, mult=mult))
    env = dict(os.environ, OMP_NUM_THREADS='2')
    proc = subprocess.run([sys.executable, '-m', 'oqp.pyoqp', deck.name], cwd=tmp_path, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=300)
    log = (tmp_path / 'ch2o_trah.log').read_text(errors='ignore')
    diagnostic = proc.stdout.decode(errors='ignore')[-2000:] + '\n' + log[-4000:]
    assert proc.returncode == 0, diagnostic
    assert 'SCF convergence achieved' in log, diagnostic
    residuals = re.findall(r'final fresh-Fock residual\s*=\s*([0-9.E+-]+)', log)
    assert residuals and float(residuals[-1]) < 1e-11, diagnostic
    steps = re.findall(r'CONVERGED \(FP precision, (\d+) refinement steps\)', log)
    assert not steps or int(steps[-1]) <= 6, diagnostic
