"""The native TRAH objective must not depend on previous trial densities."""
import os
from pathlib import Path
import re
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_revisiting_orbitals_restores_the_same_fock_and_energy(tmp_path):
    compiler = os.environ.get("FC") or shutil.which("gfortran-15") or shutil.which("gfortran")
    if not compiler:
        pytest.skip("Fortran compiler required")
    source = (ROOT / "source/trah_converger.F90").read_text()
    rebuild = re.search(r"  subroutine rebuild_fock\(.*?  end subroutine rebuild_fock", source, re.S).group()
    # Exercise the production provider with a small Fock builder whose
    # incremental path accumulates a screened-integral error. A full build
    # returns the same objective whenever the same orbitals are revisited.
    harness = """
module fixture
 implicit none
 integer,parameter::dp=kind(1d0)
 type basis_set
  integer::nbf=1
 end type
 type control_t
  integer::scftype=3
 end type
 type information
  type(control_t)::control
  type(basis_set),pointer::basis
 end type
 type dft_grid_t
 end type
 type trah_converger
  real(dp)::dens(1,2)=0,fock_ao(1,2)=0,f_old(1,2)=0,d_old(1,2)=0
 end type
 type scf_energy_t
  real(dp)::value=0
 end type
contains
 subroutine get_ab_initio_density(da,ma,db,mb,infos,basis)
  real(dp)::da(:),ma(:,:),db(:),mb(:,:)
  type(information)::infos
  type(basis_set)::basis
  da=ma(:,1)**2;db=mb(:,1)**2
 end subroutine
 subroutine calc_fock(basis,infos,grid,fock,e,ma,dens,mb,nschwz,f_old,d_old)
  type(basis_set)::basis
  type(information)::infos
  type(dft_grid_t)::grid
  type(scf_energy_t)::e
  real(dp)::fock(:,:),ma(:,:),dens(:,:),mb(:,:)
  real(dp),optional::f_old(:,:),d_old(:,:)
  integer::nschwz
  fock=2*dens
  if(present(f_old).and.present(d_old))then
   fock=f_old+2*(dens-d_old)
   if(any(d_old/=0))fock=fock+1d-8
   f_old=fock;d_old=dens
  endif
  e%value=sum(dens*fock);nschwz=0
 end subroutine
""" + rebuild + """
end module
program check_history
 use fixture
 type(information)::infos
 type(basis_set),target::basis
 type(dft_grid_t)::grid
 type(trah_converger)::conv
 type(scf_energy_t)::e
 real(dp)::ma(1,1),mb(1,1),first_energy,first_fock(1,2)
 integer::n
 infos%basis=>basis
 ma=1;mb=1
 call rebuild_fock(infos,grid,conv,e,ma,mb,n,.true.)
 first_energy=e%value;first_fock=conv%fock_ao
 ma=0.8_dp;mb=0.8_dp
 call rebuild_fock(infos,grid,conv,e,ma,mb,n,.true.)
 ma=1;mb=1
 call rebuild_fock(infos,grid,conv,e,ma,mb,n,.true.)
 if(abs(e%value-first_energy)>1d-14)error stop 'history-dependent energy'
 if(maxval(abs(conv%fock_ao-first_fock))>1d-14)error stop 'history-dependent Fock'
 print *, 'PASS: orbital objective is independent of trial history'
end program
"""
    program = tmp_path / "history.f90"
    program.write_text(harness)
    exe = tmp_path / "history"
    subprocess.run([compiler, "-ffree-line-length-none", str(program), "-o", str(exe)],
                   cwd=tmp_path, check=True, capture_output=True, text=True)
    subprocess.run([str(exe)], check=True, capture_output=True, text=True, timeout=10)


def test_trah_refreshes_periodically_and_keeps_full_builds_in_tight_tail(tmp_path):
    compiler = os.environ.get("FC") or shutil.which("gfortran-15") or shutil.which("gfortran")
    if not compiler:
        pytest.skip("Fortran compiler required")
    source = (ROOT / "source/trah_converger.F90").read_text()
    callback = re.search(r"  subroutine scf_grad_hdiag\(.*?  end subroutine scf_grad_hdiag", source, re.S).group()
    # The phase thresholds are taken from the converger, not restated here.
    thresholds = "\n".join(re.findall(r"^\s*real\(dp\), parameter :: arh_gnorm_\w+ *=.*$", source, re.M))
    assert thresholds.count("arh_gnorm_") == 2, thresholds
    harness = """
module fixture
 implicit none
 integer,parameter::dp=kind(1d0)
 integer,parameter::IW=6
 """ + thresholds + """
 type control_t
  integer::verbose=1
 end type
 type infos_t
  type(control_t)::control
 end type
 ! ARH hooks of trah_converger: record/rebuild only, no Fock work.
 type conv_t
  real(dp)::mo_a(1,1),mo_b(1,1)
  integer::arh_nb=0,pushes=0,arh_mem=20
 contains
  procedure::arh_push,arh_set_reference
  procedure::arh_linearity_residual
 end type
 type scf_trah_provider_t
  type(infos_t)::infos
  integer::molgrid=0,energy=0,model_builds=0,fock_builds=0
  logical::full_fock=.false.,arh_tail=.false.,model_hessian=.false.,model_only=.false.
  type(conv_t)::conv
 end type
 real(dp)::residual=1d-2
 integer::calls=0
 logical::last_full=.false.
contains
 subroutine arh_push(self)
  class(conv_t)::self
  self%pushes=self%pushes+1
 end subroutine
 subroutine arh_set_reference(self)
  class(conv_t)::self
 end subroutine
 real(dp) function arh_linearity_residual(self,infos)
  class(conv_t)::self
  type(infos_t)::infos
  arh_linearity_residual=0
 end function
 subroutine build_fock_grad(infos,grid,conv,energy,ma,mb,g,h,e,full)
  type(infos_t)::infos
  integer::grid,energy
  type(conv_t)::conv
  real(dp)::ma(:,:),mb(:,:),g(:),h(:),e
  logical::full
  calls=calls+1;last_full=full
  g=residual;h=1;e=0
 end subroutine
""" + callback + """
end module
program check_policy
 use fixture
 type(scf_trah_provider_t)::p
 real(dp)::g(1),h(1),e
 integer::i,ierr
 do i=1,9
  call scf_grad_hdiag(p,g,h,e,ierr)
  if(last_full)error stop 'premature full Fock'
 enddo
 call scf_grad_hdiag(p,g,h,e,ierr)
 if(.not.last_full)error stop 'missing periodic refresh'
 call scf_grad_hdiag(p,g,h,e,ierr)
 if(last_full)error stop 'incremental descent disabled'
 residual=1d-5
 call scf_grad_hdiag(p,g,h,e,ierr)
 if(.not.p%full_fock.or..not.last_full.or.calls/=13)error stop 'missing tight refresh'
 ! A subsequent larger gradient must not mix incremental and full energies.
 residual=1d-2
 call scf_grad_hdiag(p,g,h,e,ierr)
 if(.not.last_full)error stop 'tight mode was lost'
 ! Every gradient evaluation is an ARH reference point and a Fock build each.
 if(p%conv%pushes/=13)error stop 'ARH reference not recorded per gradient'
 if(p%fock_builds/=calls)error stop 'Fock builds miscounted'
 ! Hybrid phases: the tight tail (|g| < 1e-4) is exact for good, even when
 ! the gradient grows again; above it the model is used from the start.
 if(p%model_hessian.or..not.p%arh_tail)error stop 'tail must stay exact'
 p%arh_tail=.false.;residual=1d-3
 call scf_grad_hdiag(p,g,h,e,ierr)
 if(.not.p%model_hessian)error stop 'model not used before the tail'
 residual=1d-1
 call scf_grad_hdiag(p,g,h,e,ierr)
 if(.not.p%model_hessian)error stop 'model not used far from convergence'
 ! converger_type=arh: the model is the method, so the tail keeps it too.
 p%model_only=.true.;residual=1d-6
 call scf_grad_hdiag(p,g,h,e,ierr)
 if(p%arh_tail.or..not.p%model_hessian)error stop 'ARH must keep the model in the tail'
end program
"""
    program = tmp_path / "policy.f90"
    program.write_text(harness)
    exe = tmp_path / "policy"
    subprocess.run([compiler, "-ffree-line-length-none", str(program), "-o", str(exe)],
                   cwd=tmp_path, check=True, capture_output=True, text=True)
    subprocess.run([str(exe)], check=True, capture_output=True, text=True, timeout=10)
