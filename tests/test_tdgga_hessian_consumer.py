"""Finite-difference checks of the direct restricted-GGA Hessian consumer."""

from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]

DRIVER = r"""
program test_tdgga_hessian_consumer
  use precision, only: fp
  use mod_dft_gridint_tdxc_deriv, only: gga_tdxc_kernel_t, &
    gga_tdxc_variables_t, gga_tdxc_fixed_value
  use mod_dft_gridint_tdgga_hessian, only: gga_add_owner_motion, &
    gga_build_variable_derivatives, gga_weighted_point_hessian
  implicit none
  integer, parameter :: nat=1
  integer :: i
  real(fp), parameter :: h=2.0e-5_fp
  type(gga_tdxc_kernel_t) :: kernel
  type(gga_tdxc_variables_t) :: u
  real(fp) :: dr(3,nat),dp(3,nat),dv(3,nat)
  real(fp) :: dgr(3,3,nat),dgp(3,3,nat),dgv(3,3,nat)
  real(fp) :: d2r(3,3,nat,nat),d2p(3,3,nat,nat),d2v(3,3,nat,nat)
  real(fp) :: d2gr(3,3,3,nat,nat),d2gp(3,3,3,nat,nat)
  real(fp) :: d2gv(3,3,3,nat,nat),du(7,3,nat),d2u(7,3,3,nat,nat)
  real(fp) :: hh(3,3,nat,nat),fd, e0, ep, em, epp, emm
  real(fp) :: gr(3),gp(3),gv(3),r,p,v,w,dw(3,nat),d2w(3,nat,3,nat)
  real(fp) :: f1(3,2),fg1(3,3,2),f2(3,3,2,2),fg2(3,3,3,2,2)
  real(fp) :: o1(3,2),og1(3,3,2),o2(3,3,2,2),og2(3,3,3,2,2)

  call data_at_zero(r,p,v,gr,gp,gv,dr,dp,dv,dgr,dgp,dgv, &
    d2r,d2p,d2v,d2gr,d2gp,d2gv)
  call kernel_at(r,dot_product(gr,gr),kernel)
  call gga_build_variable_derivatives(r,p,v,gr,gp,gv,dr,dp,dv, &
    dgr,dgp,dgv,d2r,d2p,d2v,d2gr,d2gp,d2gv,u,du,d2u)
  w=0.8_fp; dw=0; dw(1,1)=0.04_fp; d2w=0; d2w(1,1,1,1)=0.04_fp
  hh=0
  call gga_weighted_point_hessian(kernel,u,du,d2u,1.3_fp,w,dw,d2w,hh)
  e0=energy(0.0_fp); ep=energy(h); em=energy(-h)
  epp=energy(2*h); emm=energy(-2*h)
  fd=(-epp+16*ep-30*e0+16*em-emm)/(12*h*h)
  print '(2es24.15)',hh(1,1,1,1),fd
  if (abs(hh(1,1,1,1)-fd)>2.0e-6_fp) error stop 'chain-rule Hessian failed'

  f1=reshape([(0.01_fp*i,i=1,size(f1))],shape(f1))
  fg1=reshape([(0.003_fp*i,i=1,size(fg1))],shape(fg1))
  f2=reshape([(0.0007_fp*i,i=1,size(f2))],shape(f2))
  fg2=reshape([(0.0001_fp*i,i=1,size(fg2))],shape(fg2))
  call gga_add_owner_motion(1,f1,fg1,f2,fg2,o1,og1,o2,og2)
  print '(4es16.7)',maxval(abs(sum(o1,dim=2))), &
    maxval(abs(sum(og1,dim=3))),abs(sum(o2)),abs(sum(og2))
  if (maxval(abs(sum(o1,dim=2)))>1e-13_fp .or. &
      maxval(abs(sum(og1,dim=3)))>1e-13_fp .or. &
      abs(sum(o2))>1e-13_fp .or. abs(sum(og2))>1e-13_fp) &
    error stop 'owner-motion translational sum rule failed'

contains
  subroutine data_at_zero(r,p,v,gr,gp,gv,dr,dp,dv,dgr,dgp,dgv, &
      d2r,d2p,d2v,d2gr,d2gp,d2gv)
    real(fp),intent(out)::r,p,v,gr(3),gp(3),gv(3)
    real(fp),intent(out)::dr(3,nat),dp(3,nat),dv(3,nat)
    real(fp),intent(out)::dgr(3,3,nat),dgp(3,3,nat),dgv(3,3,nat)
    real(fp),intent(out)::d2r(3,3,nat,nat),d2p(3,3,nat,nat),d2v(3,3,nat,nat)
    real(fp),intent(out)::d2gr(3,3,3,nat,nat),d2gp(3,3,3,nat,nat)
    real(fp),intent(out)::d2gv(3,3,3,nat,nat)
    r=.7_fp;p=.4_fp;v=.3_fp
    gr=[.5_fp,-.2_fp,.1_fp];gp=[-.1_fp,.3_fp,.2_fp];gv=[.2_fp,.15_fp,-.25_fp]
    dr=0;dp=0;dv=0;dgr=0;dgp=0;dgv=0;d2r=0;d2p=0;d2v=0
    d2gr=0;d2gp=0;d2gv=0
    dr(1,1)=.2_fp; dp(1,1)=-.12_fp; dv(1,1)=.08_fp
    d2r(1,1,1,1)=.06_fp; d2p(1,1,1,1)=.10_fp; d2v(1,1,1,1)=-.04_fp
    dgr(:,1,1)=[.04_fp,-.03_fp,.02_fp]
    dgp(:,1,1)=[-.02_fp,.01_fp,.05_fp]
    dgv(:,1,1)=[.03_fp,.02_fp,-.01_fp]
    d2gr(:,1,1,1,1)=[.02_fp,.01_fp,-.03_fp]
    d2gp(:,1,1,1,1)=[.01_fp,-.04_fp,.02_fp]
    d2gv(:,1,1,1,1)=[-.02_fp,.03_fp,.01_fp]
  end subroutine

  function energy(t) result(e)
    real(fp),intent(in)::t
    real(fp)::e,rr,pp,vv,gg(3),gpp(3),gvv(3),ww
    type(gga_tdxc_kernel_t)::kk
    type(gga_tdxc_variables_t)::uu
    rr=.7_fp+.2_fp*t+.03_fp*t*t
    pp=.4_fp-.12_fp*t+.05_fp*t*t
    vv=.3_fp+.08_fp*t-.02_fp*t*t
    gg=[.5_fp+.04_fp*t+.01_fp*t*t,-.2_fp-.03_fp*t+.005_fp*t*t, &
        .1_fp+.02_fp*t-.015_fp*t*t]
    gpp=[-.1_fp-.02_fp*t+.005_fp*t*t,.3_fp+.01_fp*t-.02_fp*t*t, &
         .2_fp+.05_fp*t+.01_fp*t*t]
    gvv=[.2_fp+.03_fp*t-.01_fp*t*t,.15_fp+.02_fp*t+.015_fp*t*t, &
         -.25_fp-.01_fp*t+.005_fp*t*t]
    call kernel_at(rr,dot_product(gg,gg),kk)
    uu%r=rr;uu%s=dot_product(gg,gg);uu%p=pp
    uu%pg=dot_product(gg,gpp);uu%v=vv;uu%gv=dot_product(gg,gvv)
    uu%w=dot_product(gvv,gvv)
    ww=.8_fp+.04_fp*t+.02_fp*t*t
    e=1.3_fp*ww*gga_tdxc_fixed_value(kk,uu)
  end function

  subroutine kernel_at(r,s,k)
    real(fp),intent(in)::r,s
    type(gga_tdxc_kernel_t),intent(out)::k
    k%fr=4*r**3+6*r*r*s+6*r*s*s+4*s**3
    k%fs=2*r**3+6*r*r*s+12*r*s*s+20*s**3
    k%frr=12*r*r+12*r*s+6*s*s
    k%frs=6*r*r+12*r*s+12*s*s
    k%fss=6*r*r+24*r*s+60*s*s
    k%frrr=24*r+12*s;k%frrs=12*r+12*s
    k%frss=12*r+24*s;k%fsss=24*r+120*s
    k%frrrr=24;k%frrrs=12;k%frrss=12;k%frsss=24;k%fssss=120
  end subroutine
end program
"""


class TestTDGGAHessianConsumer(unittest.TestCase):
    def test_production_gga_dispatch_and_restricted_spin_factors(self):
        driver = (ROOT / "source" / "dftlib" /
                  "dft_gridint_tdgga_hessian_driver.F90").read_text()
        assembly = (ROOT / "source" / "modules" /
                    "tdhf_hessian_xc.F90").read_text()
        self.assertIn("2.0_fp*density_p(:,i)", driver)
        self.assertIn("call tddft_gga_fixed_hessian", assembly)
        self.assertIn("if (infos%functional%needGrd) then", assembly)
        self.assertIn("call tddft_lda_fixed_hessian", assembly)
        self.assertNotIn("build_scalar_fixed_hessian", assembly)

    def test_gamess_gga_cutoff_and_fourth_kernel_step_are_preserved(self):
        driver = (ROOT / "source" / "dftlib" /
                  "dft_gridint_tdgga_hessian_driver.F90").read_text().lower()
        point = (ROOT / "source" / "dftlib" /
                 "dft_gridint_tdgga_hessian.F90").read_text().lower()
        self.assertIn("eps_r=1.0e-3_fp, eps_s=1.0e-3_fp", driver)
        self.assertIn("rho_cutoff=1.0e-8_fp", driver)
        self.assertIn("eps_r,eps_s,rho_cutoff,kernels(i)", driver)
        self.assertIn("if (r < rho_cutoff) return", point)

    def test_full_chain_rule_weight_and_owner_motion(self):
        compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
        if compiler is None:
            self.skipTest("GNU Fortran compiler is not available")
        with tempfile.TemporaryDirectory(prefix="oqp-tdgga-hess-") as td:
            tmp = Path(td)
            driver = tmp / "driver.F90"
            driver.write_text(DRIVER)
            exe = tmp / "test_tdgga_hessian_consumer"
            subprocess.run([
                compiler,
                str(ROOT / "source" / "precision.F90"),
                str(ROOT / "source" / "dftlib" / "dft_partfunc.F90"),
                str(ROOT / "source" / "dftlib" / "dft_gridint_tdxc_deriv.F90"),
                str(ROOT / "source" / "dftlib" / "dft_gga_nuclear_point.F90"),
                str(ROOT / "source" / "dftlib" / "dft_partition_hessian.F90"),
                str(ROOT / "source" / "dftlib" / "dft_gridint_tdgga_hessian.F90"),
                str(driver), "-o", str(exe),
            ], cwd=tmp, check=True)
            subprocess.run([str(exe)], check=True, cwd=tmp)


if __name__ == "__main__":
    unittest.main()
