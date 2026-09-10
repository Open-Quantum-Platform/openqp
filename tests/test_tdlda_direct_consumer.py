"""Point-chain and production-dispatch checks for direct LDA TDDFT Hessians."""
from pathlib import Path
import shutil, subprocess, tempfile

ROOT=Path(__file__).resolve().parents[1]

PROGRAM=r"""
program test_lda
 use precision,only:fp
 use mod_dft_gridint_tdlda_point,only:lda_weighted_fixed_hessian,lda_weighted_response_row
 implicit none
 real(fp),parameter::s=1.2_fp,hstep=1e-4_fp
 real(fp)::dr(3,1),dp(3,1),dx(3,1),d2(3,3,1,1),dw(3,1),d2w(3,1,3,1)
 real(fp)::hh(3,3,1,1),row(3,1),af,ar,nf,nr,r,p,x,rd,rx
 r=.7_fp;p=.3_fp;x=.2_fp;rd=.11_fp;rx=-.04_fp
 dr=0;dp=0;dx=0;d2=0;dw=0;d2w=0
 dr(1,1)=.13_fp;dp(1,1)=-.06_fp;dx(1,1)=.05_fp;dw(1,1)=.03_fp
 hh=0;call lda_weighted_fixed_hessian(r**4,4*r**3,12*r**2,24*r,p,x,dr,dp,dx,d2,d2,d2,s,.8_fp,dw,d2w,hh)
 row=0;call lda_weighted_response_row(4*r**3,12*r**2,24*r,p,x,rd,rx,dr,dp,dx, &
   reshape([-.02_fp,0._fp,0._fp],[3,1]),reshape([.015_fp,0._fp,0._fp],[3,1]),s,.8_fp,dw,row)
 af=hh(1,1,1,1);ar=row(1,1)
 nf=(energy(hstep,0._fp)-2*energy(0._fp,0._fp)+energy(-hstep,0._fp))/hstep**2
 nr=(energy(hstep,hstep)-energy(hstep,-hstep)-energy(-hstep,hstep)+energy(-hstep,-hstep))/(4*hstep**2)
 if(abs(af-nf)>2e-6_fp.or.abs(ar-nr)>2e-6_fp)error stop 'LDA direct chain failed'
contains
 function energy(t,l)result(e)
  real(fp),intent(in)::t,l;real(fp)::e,rr,pp,xx,rrd,rrx,w
  rr=.7_fp+.13_fp*t+l*(.11_fp-.02_fp*t);pp=.3_fp-.06_fp*t
  xx=.2_fp+.05_fp*t+l*(-.04_fp+.015_fp*t);w=.8_fp+.03_fp*t
  e=s*w*(rr**4*pp+16*rr**3*xx**2)
 end function
end program
"""

def test_direct_lda_dispatch_and_density_factors():
    driver=(ROOT/'source/dftlib/dft_gridint_tdlda_driver.F90').read_text()
    assembly=(ROOT/'source/modules/tdhf_hessian_xc.F90').read_text()
    assert '2.0_fp*p(:,i)' in driver
    assert 'dat%dd(:,i,k)=dd(:,i,k)' in driver
    assert 'dat%du(:,i,k)=du(:,i,k)' in driver
    assert 'call tddft_lda_fixed_hessian' in assembly
    assert 'call tddft_lda_response_rows' in assembly
    assert 'eps=1e-4_fp' in driver

def test_direct_lda_point_chain(tmp_path):
    compiler=shutil.which('gfortran-15') or shutil.which('gfortran')
    if not compiler:return
    src=tmp_path/'driver.F90';src.write_text(PROGRAM);exe=tmp_path/'test'
    subprocess.run([compiler,'-ffree-line-length-none',str(ROOT/'source/precision.F90'),
      str(ROOT/'source/dftlib/dft_gridint_tdxc_deriv.F90'),
      str(ROOT/'source/dftlib/dft_gridint_tdlda_point.F90'),str(src),'-o',str(exe)],check=True,cwd=tmp_path)
    subprocess.run([str(exe)],check=True,cwd=tmp_path)
