"""Standalone checks of the pure TDDFT GGA Hessian point algebra."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "source" / "dftlib" / "dft_gridint_tdxc_deriv.F90"


DRIVER = r"""
program test_tdxc_gga_algebra
  use precision, only: fp
  use mod_dft_gridint_tdxc_deriv
  implicit none
  type(gga_tdxc_kernel_t) :: k, rp, rm, sp, sm
  type(gga_tdxc_variables_t) :: u
  real(fp) :: d1(5), d2(5,5), d3(5,5,5), fa(7), fab(7,7)
  real(fp) :: gr(3), gp(3), dgr(3), dgp(3), a, c, b(3), da, dc, db(3)
  real(fp) :: a3, c2, c3, b3(3), da3, dc2, dc3, db3(3)

  d1=1; d2=1; d3=1
  call gga_total_kernel_from_spin(d1,d2,d3,k)
  write(*,'(*(ES25.16,1X))') k%fr,k%fs,k%frr,k%frs,k%fss, &
    k%frrr,k%frrs,k%frss,k%fsss

  rp%frrr=5; rm%frrr=1
  sp%frrr=7; sm%frrr=1
  sp%frrs=8; sm%frrs=2
  sp%frss=9; sm%frss=3
  sp%fsss=10; sm%fsss=4
  call gga_fourth_from_third(2.0_fp,3.0_fp,rp,rm,sp,sm, &
    0.1_fp,0.1_fp,1.e-8_fp,k)
  write(*,'(*(ES25.16,1X))') k%frrrr,k%frrrs,k%frrss,k%frsss,k%fssss

  k%fr=.11_fp; k%fs=.12_fp; k%frr=.13_fp; k%frs=.14_fp; k%fss=.15_fp
  k%frrr=.16_fp; k%frrs=.17_fp; k%frss=.18_fp; k%fsss=.19_fp
  k%frrrr=.20_fp; k%frrrs=.21_fp; k%frrss=.22_fp
  k%frsss=.23_fp; k%fssss=.24_fp
  u=gga_tdxc_variables_t(.31_fp,.32_fp,.33_fp,.34_fp,.35_fp,.36_fp,.37_fp)
  call gga_tdxc_fixed_derivatives(k,u,fa,fab)
  write(*,'(*(ES25.16,1X))') gga_tdxc_fixed_value(k,u),fa,fab

  gr=[.2_fp,-.3_fp,.4_fp]; gp=[-.5_fp,.6_fp,.7_fp]
  dgr=[.08_fp,-.09_fp,.10_fp]; dgp=[-.11_fp,.12_fp,.13_fp]
  call gga_x2_coefficients(k,gr,.45_fp,gp,a,c,b)
  call gga_x2_coefficient_derivative(k,gr,.45_fp,gp,.07_fp,dgr, &
    -.06_fp,dgp,da,dc,db)
  write(*,'(*(ES25.16,1X))') a,c,b,da,dc,db
  call gga_x3_coefficients(k,gr,.45_fp,gp,a3,c2,c3,b3)
  call gga_x3_coefficient_derivative(k,gr,.45_fp,gp,.07_fp,dgr, &
    -.06_fp,dgp,da3,dc2,dc3,db3)
  write(*,'(*(ES25.16,1X))') a3,c2,c3,b3,da3,dc2,dc3,db3
end program test_tdxc_gga_algebra
"""


def fixed_reference():
    fr, fs, frr, frs, fss = .11, .12, .13, .14, .15
    frrr, frrs, frss, fsss = .16, .17, .18, .19
    f4a, f4b, f4c, f4d, f4e = .20, .21, .22, .23, .24
    _, _, p, pg, v, gv, w = [.31, .32, .33, .34, .35, .36, .37]
    q = fr*p + 2*fs*pg + 4*(frr*v*v + 4*frs*v*gv + 4*fss*gv*gv + 2*fs*w)
    fa = [
        frr*p+2*frs*pg+4*frrr*v*v+16*frrs*v*gv+16*frss*gv*gv+8*frs*w,
        frs*p+2*fss*pg+4*frrs*v*v+16*frss*v*gv+16*fsss*gv*gv+8*fss*w,
        fr, 2*fs, 8*frr*v+16*frs*gv, 16*frs*v+32*fss*gv, 8*fs,
    ]
    h = [[0.0]*7 for _ in range(7)]
    def sym(i, j, x): h[i][j] = h[j][i] = x
    h[4][4]=8*frr; sym(4,5,16*frs); h[5][5]=32*fss
    sym(2,0,frr); sym(2,1,frs); sym(3,0,2*frs); sym(3,1,2*fss)
    sym(6,0,8*frs); sym(6,1,8*fss)
    sym(4,0,8*frrr*v+16*frrs*gv); sym(4,1,8*frrs*v+16*frss*gv)
    sym(5,0,16*frrs*v+32*frss*gv); sym(5,1,16*frss*v+32*fsss*gv)
    h[0][0]=frrr*p+2*frrs*pg+4*f4a*v*v+16*f4b*v*gv+16*f4c*gv*gv+8*frrs*w
    sym(0,1,frrs*p+2*frss*pg+4*f4b*v*v+16*f4c*v*gv+16*f4d*gv*gv+8*frss*w)
    h[1][1]=frss*p+2*fsss*pg+4*f4c*v*v+16*f4d*v*gv+16*f4e*gv*gv+8*fsss*w
    # Fortran prints arrays in column-major order.
    return [q] + fa + [h[i][j] for j in range(7) for i in range(7)]


def coefficient_reference():
    fs, frr, frs, fss = .12, .13, .14, .15
    frrr, frrs, frss, fsss = .16, .17, .18, .19
    f4a, f4b, f4c, f4d, f4e = .20, .21, .22, .23, .24
    gr, gp = [.2, -.3, .4], [-.5, .6, .7]
    dgr, dgp = [.08, -.09, .10], [-.11, .12, .13]
    p, dr, dp = .45, .07, -.06
    dot = lambda x, y: sum(a*b for a, b in zip(x, y))
    gamma = dot(gr, gp)
    dgamma = dot(dgr, gp) + dot(gr, dgp)
    gp2, dgp2 = dot(gp, gp), 2*dot(gp, dgp)
    ds = 2*dot(gr, dgr)
    dfrr, dfrs = frrr*dr+frrs*ds, frrs*dr+frss*ds
    dfss, dfs = frss*dr+fsss*ds, frs*dr+fss*ds
    a = frr*p+2*frs*gamma
    c = 2*frs*p+4*fss*gamma
    b = [c*x+2*fs*y for x, y in zip(gr, gp)]
    da = dfrr*p+frr*dp+2*(dfrs*gamma+frs*dgamma)
    dc = 2*(dfrs*p+frs*dp)+4*(dfss*gamma+fss*dgamma)
    db = [dc*x+c*dx+2*dfs*y+2*fs*dy
          for x, dx, y, dy in zip(gr, dgr, gp, dgp)]
    a3 = frrr*p*p+4*frrs*p*gamma+4*frss*gamma*gamma+2*frs*gp2
    c2 = c
    c3 = 2*frrs*p*p+8*frss*p*gamma+8*fsss*gamma*gamma+4*fss*gp2
    b3 = [c3*x+2*c2*y for x, y in zip(gr, gp)]
    dfrrr, dfrrs = f4a*dr+f4b*ds, f4b*dr+f4c*ds
    dfrss, dfsss = f4c*dr+f4d*ds, f4d*dr+f4e*ds
    da3 = (dfrrr*p*p+2*frrr*p*dp
           +4*(dfrrs*p*gamma+frrs*(dp*gamma+p*dgamma))
           +4*(dfrss*gamma*gamma+2*frss*gamma*dgamma)
           +2*(dfrs*gp2+frs*dgp2))
    dc2 = dc
    dc3 = (2*dfrrs*p*p+4*frrs*p*dp
           +8*(dfrss*p*gamma+frss*(dp*gamma+p*dgamma))
           +8*(dfsss*gamma*gamma+2*fsss*gamma*dgamma)
           +4*(dfss*gp2+fss*dgp2))
    db3 = [dc3*x+c3*dx+2*dc2*y+2*c2*dy
           for x, dx, y, dy in zip(gr, dgr, gp, dgp)]
    return [a, c] + b + [da, dc] + db, \
           [a3, c2, c3] + b3 + [da3, dc2, dc3] + db3


class TdxcGgaAlgebraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
        if compiler is None:
            raise unittest.SkipTest("GNU Fortran compiler is not available")
        cls.tmp = tempfile.TemporaryDirectory(prefix="oqp-tdxc-gga-")
        tmp = Path(cls.tmp.name)
        driver = tmp / "driver.F90"
        driver.write_text(DRIVER)
        cls.exe = tmp / "test_tdxc_gga_algebra"
        subprocess.run([
            compiler, str(ROOT / "source" / "precision.F90"), str(MODULE),
            str(driver), "-o", str(cls.exe)
        ], cwd=tmp, check=True)
        cls.lines = [[float(x) for x in line.split()] for line in
                     subprocess.check_output([str(cls.exe)], text=True).splitlines()]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def assert_close(self, got, expected, tol=2e-13):
        self.assertEqual(len(got), len(expected))
        self.assertLessEqual(max(abs(a-b) for a, b in zip(got, expected)), tol)

    def test_tdhxkg_closed_shell_spin_factors(self):
        self.assert_close(self.lines[0], [1, .75, 1, .75, .5625, 1, .75, .5625, .421875])

    def test_tdhxfd4_central_difference(self):
        self.assert_close(self.lines[1], [10, 10, 10, 10, 10])

    def test_tdhxfab_point_value_and_derivatives(self):
        self.assert_close(self.lines[2], fixed_reference())

    def test_tdhxgpg_x2_x3_coefficients_and_derivatives(self):
        x2, x3 = coefficient_reference()
        self.assert_close(self.lines[3], x2)
        self.assert_close(self.lines[4], x3)


if __name__ == "__main__":
    unittest.main()
