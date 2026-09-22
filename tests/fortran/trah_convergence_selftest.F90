! Behavioral acceptance tests against the actual method-independent TRAH loop.
module trah_test_provider
 use precision, only: dp
 use trah_core_mod
 implicit none
 type,extends(trah_provider_t)::model
  real(dp)::x=1.0_dp, residual=5e-8_dp, curvature=1.0_dp
  logical::stuck=.false.
  logical::check_cache=.false.,trial_cache=.false.
  integer::reject_trials=0
  integer::evaluations=0,accepted_steps=0,hess_evaluations=0
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
 this%trial_cache=.false.
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
 this%hess_evaluations=this%hess_evaluations+1
 if(this%check_cache.and.this%trial_cache)then
  ierr=77
  return
 endif
 hv=this%curvature*v;ierr=0
 end subroutine
 subroutine te(this,p,e,ierr)
 class(model),intent(inout)::this
 real(dp),intent(in)::p(:)
 real(dp),intent(out)::e
 integer,intent(out)::ierr
 e=0
 if(.not.this%stuck)e=0.5_dp*(this%x+p(1))**2
 this%trial_cache=.true.
 if(this%reject_trials>0)then
  this%reject_trials=this%reject_trials-1
  e=1.0_dp+0.5_dp*this%x**2
 endif
 ierr=0
 end subroutine
 subroutine ap(this,p,ierr)
 class(model),intent(inout)::this
 real(dp),intent(in)::p(:)
 integer,intent(out)::ierr
 this%x=this%x+p(1);ierr=0
 this%accepted_steps=this%accepted_steps+1
 end subroutine
end module
program check_trah_convergence
 use trah_test_provider
 implicit none
 type(model)::p
 type(trah_params_t)::par
 type(trah_result_t)::res
 integer::nh,hi(1),used,ierr
 real(dp)::he(3),hde(2),hg(3),hs(3),g1(1),h1(1),step1(1),pred1,expected
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
 ! SCF trial energies overwrite Fock caches. Reject one full step and all
 ! five line-search trials, then require a refreshed model before H.v.
 p%x=1;p%check_cache=.true.;p%reject_trials=6;p%trial_cache=.false.
 p%refresh_on_rejection=.true.
 par%nmac=20;par%r0=0.4_dp
 call trah_run(p,par,res)
 if(.not.res%converged.or.res%ierr/=0.or.res%error>=par%conv_tol)error stop 'rejected trial cache'
 if(p%reject_trials/=0)error stop 'rejection case not exercised'
 ! Providers with an independent accepted-point Hessian (such as CASSCF)
 ! must not rebuild it when only the trial point changes.
 p%x=1;p%check_cache=.false.;p%refresh_on_rejection=.false.;p%reject_trials=6
 p%evaluations=0;p%accepted_steps=0
 call trah_run(p,par,res)
 if(.not.res%converged.or.res%ierr/=0)error stop 'independent Hessian convergence'
 if(p%reject_trials/=0)error stop 'independent Hessian rejection not exercised'
 if(p%evaluations/=p%accepted_steps+1)error stop 'unnecessary accepted-point rebuild'
 ! Optional history may be absent or have unequal capacities.
 par%want_history=.true.;p%x=1
 nh=-1
 call trah_run(p,par,res,nhist=nh)
 if(.not.res%converged.or.nh/=0)error stop 'absent history arrays'
 p%x=1;he=-999;hde=-999;hg=-999;hs=-999
 call trah_run(p,par,res,hi,he,hde,hg,hs,nh)
 if(.not.res%converged.or.nh/=1)error stop 'history minimum capacity'
 if(any(he(2:)/=-999).or.hde(2)/=-999)error stop 'history overwritten'
 ! The augmented-Hessian solve already forms A*u. Its tail determines H*p,
 ! including the eigenvector sign and trust-radius scaling, without another
 ! Hessian-vector product.
 p%hess_evaluations=0;p%curvature=2
 par%deterministic=.false.;par%sub_solver=1;par%nrtv=1;par%nmic=6
 g1=0.2_dp;h1=2.0_dp
 call trah_micro_step(p,par,g1,h1,0.05_dp,1,step1,pred1,used,ierr)
 if(ierr/=0)error stop 'Davidson micro-solve failed'
 if(p%hess_evaluations/=used)error stop 'redundant Davidson Hessian product'
 expected=-(g1(1)*step1(1)+0.5_dp*p%curvature*step1(1)*step1(1))
 if(abs(pred1-expected)>1e-13_dp)error stop 'Davidson predicted reduction'
 if(abs(step1(1))>0.05_dp+1e-14_dp)error stop 'Davidson trust radius'
 print *, 'PASS: quadratic, precision stagnation, trust collapse, maximum iterations'
end program
