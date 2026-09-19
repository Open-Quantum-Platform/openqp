! Behavioral acceptance tests against the actual method-independent TRAH loop.
module trah_test_provider
 use precision, only: dp
 use trah_core_mod
 implicit none
 type,extends(trah_provider_t)::model
  real(dp)::x=1.0_dp, residual=5e-8_dp, curvature=1.0_dp
  logical::stuck=.false.
  integer::evaluations=0
 contains
  procedure::grad_hdiag=>gh
  procedure::hess_vec=>model_hv
  procedure::trial_energy=>te
  procedure::apply_step=>ap
 end type
 contains
 subroutine gh(this,g,hdiag,e,ierr)
 class(model),intent(inout)::this
 real(dp),intent(out)::g(:),hdiag(:),e
 integer,intent(out)::ierr
 this%evaluations=this%evaluations+1
 if(this%stuck)then
  g=this%residual;e=0
 else
  g=this%x;e=0.5_dp*this%x**2
 endif
 hdiag=this%curvature;ierr=0
 end subroutine
 subroutine model_hv(this,v,hv,ierr)
 class(model),intent(inout)::this
 real(dp),intent(in)::v(:)
 real(dp),intent(out)::hv(:)
 integer,intent(out)::ierr
 hv=this%curvature*v;ierr=0
 end subroutine
 subroutine te(this,p,e,ierr)
 class(model),intent(inout)::this
 real(dp),intent(in)::p(:)
 real(dp),intent(out)::e
 integer,intent(out)::ierr
 e=0
 if(.not.this%stuck)e=0.5_dp*(this%x+p(1))**2
 ierr=0
 end subroutine
 subroutine ap(this,p,ierr)
 class(model),intent(inout)::this
 real(dp),intent(in)::p(:)
 integer,intent(out)::ierr
 this%x=this%x+p(1);ierr=0
 end subroutine
end module
program check_trah_convergence
 use trah_test_provider
 implicit none
 type(model)::p
 type(trah_params_t)::par
 type(trah_result_t)::res
 p%nparam=1
 par%deterministic=.true.;par%sub_solver=2;par%conv_tol=1e-8_dp
 par%nmac=20;par%verbose=.false.
 call trah_run(p,par,res)
 if(.not.res%converged.or.res%ierr/=0.or.res%error>=par%conv_tol)error stop 'quadratic convergence'
 ! A constant small residual used to be forced below the requested tolerance.
 p%stuck=.true.;p%residual=5e-8_dp;p%evaluations=0
 call trah_run(p,par,res)
 if(res%converged.or.res%ierr==0.or.res%error<par%conv_tol)error stop 'precision false success'
 if(p%evaluations>par%nmac+1)error stop 'unbounded precision refinement'
 ! Trust collapse must not substitute the loose 1e-4 acceptance threshold.
 p%residual=5e-7_dp;p%curvature=1e-6_dp;p%evaluations=0
 call trah_run(p,par,res)
 if(res%converged.or.res%ierr==0)error stop 'trust collapse false success'
 if(abs(res%error-p%residual)>1e-15_dp)error stop 'residual changed'
 ! Exhausting the iteration budget is failure even when iter == nmac.
 p%stuck=.false.;p%x=1;p%curvature=1
 par%nmac=1;par%r0=0.01_dp
 call trah_run(p,par,res)
 if(res%converged.or.res%ierr==0.or.res%iter/=1)error stop 'maxit false success'
 if(abs(res%error-abs(p%x))>1e-15_dp)error stop 'stale final residual'
 print *, 'PASS: quadratic, precision stagnation, trust collapse, maximum iterations'
end program
