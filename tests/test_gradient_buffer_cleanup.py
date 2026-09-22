"""Execute production cleanup on partially and repeatedly allocated gradients."""
from pathlib import Path
import re
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
CASES = [
    ('hf_gradient', 'grd2_rhf_compute_data_t', 'grd2_rhf_compute_data_t_clean'),
    ('hf_gradient', 'grd2_uhf_compute_data_t', 'grd2_uhf_compute_data_t_clean'),
    ('tdhf_gradient', 'grd2_tdhf_compute_data_t', 'grd2_tdhf_compute_data_t_clean'),
    ('tdhf_sf_gradient', 'grd2_sf_compute_data_t', 'grd2_sf_compute_data_t_clean'),
    ('tdhf_mrsf_gradient', 'grd2_mrsf_compute_data_t', 'grd2_mrsf_compute_data_t_clean'),
    ('tdhf_mrsf_gradient', 'grd2_mrsf_nac_compute_data_t', 'grd2_mrsf_nac_compute_data_t_clean'),
    ('fock_deriv', 'grd2_fockprobe_data_t', 'grd2_fockprobe_clean'),
    ('fock_deriv', 'grd2_fockprobe_os_data_t', 'grd2_fockprobe_os_clean'),
]


@pytest.mark.parametrize('module,typename,routine', CASES)
def test_cleanup_releases_owned_arrays(tmp_path, module, typename, routine):
    compiler = shutil.which('gfortran-15') or shutil.which('gfortran')
    if compiler is None:
        pytest.skip('GNU Fortran required for allocation checks')
    source = (ROOT / 'source/modules' / (module + '.F90')).read_text()
    # Preserve the real member declarations (including borrowed pointers).
    def members(name):
        match = re.search(r'^\s*type[^\n]*::\s*' + name + r'\s*\n(.*?)^\s*contains',
                          source, re.M | re.S | re.I)
        body = match.group(1)
        # A match can start with a newline; find its actual type header.
        header = next(line for line in match.group().splitlines() if '::' in line)
        parent = re.search(r'extends\((\w+)\)', header, re.I)
        if parent and parent[1] == 'grd2_hf_compute_data_t':
            body = members(parent[1]) + body
        return body
    body = members(typename)
    cleanup = re.search(r'  subroutine ' + routine + r'\(.*?  end subroutine[^\n]*',
                        source, re.S).group()
    arrays = []
    for line in body.splitlines():
        line = line.split('!')[0]
        if 'allocatable' in line.lower():
            arrays.extend((m[1], m[2].count(':')) for m in re.finditer(r'(\w+)\(([: ,]+)\)', line))
    assert arrays, 'test must exercise owned allocations'
    # Use the exact complete type body to retain member kinds and initializers.
    allocate = '\n'.join('allocate(obj%' + name + '(' + ','.join('n' for _ in range(rank)) + '))'
                         for name, rank in arrays)
    checks = '\n'.join("if(allocated(obj%" + name + ")) error stop 'retained " + name + "'"
                       for name, _ in arrays)
    program = f'''
module fixture
 implicit none
 integer,parameter::dp=kind(1d0)
 type {typename}
 {body}
 contains
 procedure::clean=>{routine}
 end type
contains
{cleanup}
end module
program check
 use fixture
 implicit none
 type({typename})::obj
 integer::n
 call obj%clean()
 do n=1,3
 {allocate}
 call obj%clean()
 {checks}
 call obj%clean()
 end do
end program
'''
    src = tmp_path / 'check.f90'
    src.write_text(program)
    exe = tmp_path / 'check'
    subprocess.run([compiler, '-fcheck=all', '-ffree-line-length-none', str(src), '-o', str(exe)],
                   cwd=tmp_path, check=True, capture_output=True, text=True)
    subprocess.run([str(exe)], check=True, capture_output=True, text=True, timeout=10)
