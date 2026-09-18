"""Execute the production rescaler at the zero-projection root tie."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_downhill_equal_roots_are_phase_independent(tmp_path):
    compiler = os.environ.get('FC') or shutil.which('gfortran-15') or shutil.which('gfortran')
    if not compiler:
        pytest.skip('Fortran compiler required')
    source = (ROOT / 'source/modules/namd.F90').read_text()
    routine = re.search(r'  subroutine namd_rescale_velocities_directional\(.*?end subroutine namd_rescale_velocities_directional', source, re.S).group()
    driver = '''
end module
program check
use fixture
implicit none
real(dp)::v(3,1),d(3,1),mass(1),g,disc,expected(3,1),scale(4),initial(3,1),e0
logical::ok
integer::i,j
mass=2;scale=[1d0,-1d0,7d0,-7d0]
do j=1,2
  initial=0
  if(j==2) initial(2,1)=0.3d0
  e0=0.5d0*sum(initial**2)*mass(1)
  do i=1,4
    d(:,1)=[scale(i),0d0,0d0];v=initial
    call namd_rescale_velocities_directional(v,mass,d,-0.5d0,ok,g,disc)
    if(.not.ok) error stop 'downhill hop rejected'
    if(abs(0.5d0*sum(v**2)*mass(1)-0.5d0-e0)>1d-12) error stop 'energy changed'
    if(i==1) expected=v
    if(maxval(abs(v-expected))>1d-12) error stop 'NAC phase or scale changed velocity'
  end do
end do
print *, 'PASS: equal-root phase invariance and energy conservation'
end program
'''
    src=tmp_path/'tie.f90'
    src.write_text('module fixture\nimplicit none\ninteger,parameter::dp=kind(1d0)\ncontains\n'+routine+driver)
    exe=tmp_path/'tie'
    subprocess.run([compiler,'-fdefault-integer-8','-ffree-line-length-none','-fcheck=all',str(src),'-o',str(exe)],cwd=tmp_path,check=True,capture_output=True,text=True)
    result=subprocess.run([str(exe)],check=True,capture_output=True,text=True,timeout=10)
    assert 'PASS: equal-root phase invariance and energy conservation' in result.stdout
