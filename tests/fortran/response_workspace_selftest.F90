! Check retained response storage across flag, thread-count and shape changes.
subroutine response_workspace_selftest(failed) bind(C)
  use iso_c_binding, only: c_int
  use precision, only: dp
  use basis_tools, only: basis_set
  use tdhf_lib, only: int2_td_data_t
  implicit none
  integer(c_int), intent(out) :: failed
  type(basis_set) :: basis
  type(int2_td_data_t) :: response
  real(dp), target :: td(16,16,2)
  integer :: pass, ns
  logical :: ap, am, tda
  failed=0
  basis%nbf=16; basis%nshell=4
  allocate(basis%ao_offset(4), basis%naos(4))
  basis%ao_offset=[1,5,9,13]; basis%naos=4
  ! Exercise allocation reuse, shape changes, TDA's amb output, and repeated
  ! cleanup. The inactive image remains a zero single-copy borrowed view.
  td=0.0_dp; response%d2=>td
  do pass=1,6
    ns=merge(3,1,mod(pass,2)==0)
    tda=pass==3; ap=pass/=4; am=pass==4.or.pass==5
    response%tamm_dancoff=tda; response%int_apb=ap; response%int_amb=am
    response%cur_pass=1
    call response%parallel_start(basis,ns)
    if(size(response%apb,4)/=merge(ns,1,ap)) failed=failed+1
    if(size(response%amb,4)/=merge(ns,1,am.or.tda)) failed=failed+1
    if(any(response%apb/=0).or.any(response%amb/=0)) failed=failed+1
    response%apb=1.0_dp; response%amb=2.0_dp
  end do
  call response%clean(); call response%clean()
end subroutine response_workspace_selftest
