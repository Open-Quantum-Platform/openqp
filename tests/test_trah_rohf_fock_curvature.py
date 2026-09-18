"""Check the production ROHF Fock correction against orbital finite differences."""
import os
from pathlib import Path
import re
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_rohf_common_rotation_fock_curvature(tmp_path):
    compiler = os.environ.get('FC') or shutil.which('gfortran-15') or shutil.which('gfortran')
    if not compiler:
        pytest.skip('Fortran compiler required')
    source = (ROOT / 'source/scf_converger.F90').read_text()
    helpers = []
    for name in ('rohf_missing_fock_terms', 'skew_sym_k'):
        match = re.search(rf'  subroutine {name}\(.*?end subroutine(?: {name})?\s*\n',
                          source, re.S | re.I)
        assert match, f'missing production routine {name}'
        helpers.append(match.group())
    # Compile the actual production routines with only their data container and
    # packed-matrix adapter supplied by the fixture. No SCF runtime is needed.
    fixture = (ROOT / 'tests/fortran/trah_rohf_fock_selftest.F90').read_text()
    program = fixture.replace('! PRODUCTION_ROUTINES', '\n'.join(helpers))
    src = tmp_path / 'check.f90'
    src.write_text(program)
    exe = tmp_path / 'check'
    subprocess.run([compiler, '-fdefault-integer-8', '-ffree-line-length-none',
                    '-fcheck=all', str(src), '-o', str(exe)], cwd=tmp_path,
                   check=True, capture_output=True, text=True)
    result = subprocess.run([str(exe)], check=True, capture_output=True,
                            text=True, timeout=10)
    assert 'PASS: ROHF Fock curvature finite differences' in result.stdout
