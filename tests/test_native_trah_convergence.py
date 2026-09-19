"""Execute the native TRAH loop with controlled convergent/stagnant providers."""
import os
from pathlib import Path
import shutil
import subprocess
import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_native_trah_respects_requested_tolerance_and_iteration_limit(tmp_path):
    compiler = os.environ.get('FC') or shutil.which('gfortran-15') or shutil.which('gfortran')
    if not compiler:
        pytest.skip('Fortran compiler required for the native TRAH acceptance test')
    exe = tmp_path / 'trah_convergence'
    subprocess.run([
        compiler, '-fdefault-integer-8', '-ffree-line-length-none',
        str(ROOT / 'tests/fixtures/trah_core_support.f90'),
        str(ROOT / 'source/trah_core.F90'),
        str(ROOT / 'tests/fortran/trah_convergence_selftest.F90'),
        '-o', str(exe),
    ], cwd=tmp_path, check=True, capture_output=True, text=True)
    result = subprocess.run([str(exe)], check=True, capture_output=True,
                            text=True, timeout=10)
    assert 'PASS: quadratic, precision stagnation, trust collapse, maximum iterations' in result.stdout
