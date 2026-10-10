"""Focused checks for the direct restricted-GGA TDHXR1G point algebra."""

from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]

DRIVER = r"""
program test_tdgga_response
  use precision, only: fp
  use mod_dft_gridint_tdxc_deriv, only: gga_tdxc_kernel_t, &
    gga_tdxc_variables_t, gga_tdxc_fixed_value
  use mod_dft_gridint_tdgga_response, only: &
    gga_build_response_direction, &
    gga_build_response_direction_derivative, gga_weighted_response_row
  implicit none
  type(gga_tdxc_kernel_t) :: kernel
  type(gga_tdxc_variables_t) :: u
  real(fp), parameter :: h=1.0e-4_fp, qscale=1.3_fp
  real(fp) :: r,p,v,gr(3),gp(3),gv(3),rd,ru,gd(3),gu(3)
  real(fp) :: dr,dp,dv,dgr(3),dgp(3),dgv(3),drd,dru,dgd(3),dgu(3)
  real(fp) :: base_du(7,3,1),duk(7),dduk(7,3,1),dw(3,1),row(3,1)
  real(fp) :: analytic,numeric

  call fields(0.0_fp,r,p,v,gr,gp,gv,rd,gd,ru,gu)
  call field_derivatives(dr,dp,dv,dgr,dgp,dgv,drd,dgd,dru,dgu)
  call kernel_at(r,dot_product(gr,gr),kernel)
  u%r=r;u%s=dot_product(gr,gr);u%p=p;u%pg=dot_product(gr,gp)
  u%v=v;u%gv=dot_product(gr,gv);u%w=dot_product(gv,gv)
  base_du=0.0_fp
  base_du(:,1,1)=[dr,2.0_fp*dot_product(gr,dgr),dp, &
    dot_product(dgr,gp)+dot_product(gr,dgp),dv, &
    dot_product(dgr,gv)+dot_product(gr,dgv),2.0_fp*dot_product(gv,dgv)]
  call gga_build_response_direction(gr,gp,gv,rd,gd,ru,gu,duk)
  call gga_build_response_direction_derivative(gr,gp,gv,gd,gu,dgr,dgp,dgv, &
    drd,dgd,dru,dgu,dduk(:,1,1))
  dduk(:,2:3,1)=0.0_fp
  dw=0.0_fp;dw(1,1)=0.04_fp;row=0.0_fp
  call gga_weighted_response_row(kernel,u,base_du,duk,dduk,qscale, &
    0.8_fp,dw,row)
  analytic=row(1,1)
  numeric=(energy(h,h)-energy(h,-h)-energy(-h,h)+energy(-h,-h))/(4*h*h)
  print '(2es24.15)',analytic,numeric
  if (abs(analytic-numeric)>2.0e-6_fp) &
    error stop 'direct GGA response-row chain rule failed'

contains
  subroutine fields(t,r,p,v,gr,gp,gv,rd,gd,ru,gu)
    real(fp),intent(in)::t
    real(fp),intent(out)::r,p,v,gr(3),gp(3),gv(3),rd,gd(3),ru,gu(3)
    r=.7_fp+.2_fp*t;p=.4_fp-.12_fp*t;v=.3_fp+.08_fp*t
    gr=[.5_fp+.04_fp*t,-.2_fp-.03_fp*t,.1_fp+.02_fp*t]
    gp=[-.1_fp-.02_fp*t,.3_fp+.01_fp*t,.2_fp+.05_fp*t]
    gv=[.2_fp+.03_fp*t,.15_fp+.02_fp*t,-.25_fp-.01_fp*t]
    rd=.11_fp-.025_fp*t;ru=-.07_fp+.018_fp*t
    gd=[.06_fp+.012_fp*t,-.035_fp+.009_fp*t,.021_fp-.006_fp*t]
    gu=[-.022_fp+.007_fp*t,.031_fp-.011_fp*t,.014_fp+.005_fp*t]
  end subroutine

  subroutine field_derivatives(dr,dp,dv,dgr,dgp,dgv,drd,dgd,dru,dgu)
    real(fp),intent(out)::dr,dp,dv,dgr(3),dgp(3),dgv(3),drd,dgd(3),dru,dgu(3)
    dr=.2_fp;dp=-.12_fp;dv=.08_fp;dgr=[.04_fp,-.03_fp,.02_fp]
    dgp=[-.02_fp,.01_fp,.05_fp];dgv=[.03_fp,.02_fp,-.01_fp]
    drd=-.025_fp;dgd=[.012_fp,.009_fp,-.006_fp]
    dru=.018_fp;dgu=[.007_fp,-.011_fp,.005_fp]
  end subroutine

  function energy(t,l) result(e)
    real(fp),intent(in)::t,l
    real(fp)::e,r,p,v,rd,ru,gr(3),gp(3),gv(3),gd(3),gu(3),w
    type(gga_tdxc_kernel_t)::k
    type(gga_tdxc_variables_t)::x
    call fields(t,r,p,v,gr,gp,gv,rd,gd,ru,gu)
    r=r+l*rd;gr=gr+l*gd;v=v+l*ru;gv=gv+l*gu
    call kernel_at(r,dot_product(gr,gr),k)
    x%r=r;x%s=dot_product(gr,gr);x%p=p;x%pg=dot_product(gr,gp)
    x%v=v;x%gv=dot_product(gr,gv);x%w=dot_product(gv,gv)
    w=.8_fp+.04_fp*t
    e=qscale*w*gga_tdxc_fixed_value(k,x)
  end function

  subroutine kernel_at(r,s,k)
    real(fp),intent(in)::r,s
    type(gga_tdxc_kernel_t),intent(out)::k
    k%fr=4*r**3+6*r*r*s+6*r*s*s+4*s**3
    k%fs=2*r**3+6*r*r*s+12*r*s*s+20*s**3
    k%frr=12*r*r+12*r*s+6*s*s;k%frs=6*r*r+12*r*s+12*s*s
    k%fss=6*r*r+24*r*s+60*s*s
    k%frrr=24*r+12*s;k%frrs=12*r+12*s
    k%frss=12*r+24*s;k%fsss=24*r+120*s
    k%frrrr=24;k%frrrs=12;k%frrss=12;k%frsss=24;k%fssss=120
  end subroutine
end program
"""


class TestTDGGAResponseConsumer(unittest.TestCase):
    def test_direct_gga_dispatch_and_normalization(self):
        driver = (ROOT / "source" / "dftlib" /
                  "dft_gridint_tdgga_response_driver.F90").read_text()
        assembly = (ROOT / "source" / "modules" /
                    "tdhf_hessian_xc.F90").read_text()
        self.assertIn("2.0_fp*density_p(:,i)", driver)
        self.assertIn("dat%density_d(:,i,k) = density_d(:,i,k)", driver)
        self.assertIn("dat%density_u(:,i,k) = density_u(:,i,k)", driver)
        self.assertIn("call tddft_gga_response_rows", assembly)
        self.assertIn("dground,dxpy,direct_rows", assembly)

    def test_point_chain_rule_against_mixed_finite_difference(self):
        compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
        if compiler is None:
            self.skipTest("GNU Fortran compiler is not available")
        with tempfile.TemporaryDirectory(prefix="oqp-tdgga-response-") as td:
            tmp = Path(td)
            source = tmp / "driver.F90"
            source.write_text(DRIVER)
            executable = tmp / "test_tdgga_response"
            subprocess.run([
                compiler,
                "-std=f2018",
                str(ROOT / "source" / "precision.F90"),
                str(ROOT / "source" / "dftlib" /
                    "dft_gridint_tdxc_deriv.F90"),
                str(ROOT / "source" / "dftlib" /
                    "dft_gridint_tdgga_response.F90"),
                str(source), "-o", str(executable),
            ], cwd=tmp, check=True)
            subprocess.run([str(executable)], check=True, cwd=tmp)


if __name__ == "__main__":
    unittest.main()
