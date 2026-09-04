program test_tdhf_mrsf_hessian_rows

  use precision,only: dp
  use tdhf_mrsf_hessian_row_primitives_mod,only: &
    MRSF_ROWS_SUCCESS,MRSF_ROWS_XC_INCOMPLETE, &
    MRSF_ROWS_INCONSISTENT_BALL,assemble_mrsf_one_e_primitive_rows, &
    combine_mrsf_response_row_blocks,check_mrsf_seven_density_alias

  implicit none

  real(kind=dp) :: ds(2,2,2),dh(2,2,2),drelaxed(2,2,2,2)
  real(kind=dp) :: dw(2,2,2),one_e(2,2),expected_one(2,2)
  real(kind=dp) :: two_response(2,2),two_reference(2,2),xc(2,2)
  real(kind=dp) :: rows(2,2),rows_one(2,2),rows_two(2,2),rows_xc(2,2)
  real(kind=dp) :: seven_derivative(7,2,2,2),td_abxc_derivative(2,2,2)
  integer :: status

  ds=0.0_dp
  dh=0.0_dp
  drelaxed=0.0_dp
  dw=0.0_dp

  ds(:,:,1)=reshape([1.0_dp,0.2_dp,0.2_dp,2.0_dp],[2,2])
  ds(:,:,2)=reshape([0.1_dp,0.3_dp,0.7_dp,-0.4_dp],[2,2])
  dh(:,:,1)=reshape([3.0_dp,0.1_dp,0.1_dp,-1.0_dp],[2,2])
  dh(:,:,2)=reshape([0.2_dp,0.4_dp,0.6_dp,0.3_dp],[2,2])
  dw(:,:,1)=reshape([0.2_dp,0.0_dp,0.0_dp,-0.3_dp],[2,2])
  dw(:,:,2)=reshape([0.0_dp,0.2_dp,0.4_dp,0.0_dp],[2,2])
  drelaxed(:,:,1,1)=reshape([1.0_dp,0.0_dp,0.0_dp,2.0_dp],[2,2])
  drelaxed(:,:,2,1)=reshape([-0.5_dp,0.0_dp,0.0_dp,0.25_dp],[2,2])
  drelaxed(:,:,1,2)=reshape([0.0_dp,0.2_dp,0.6_dp,0.0_dp],[2,2])
  drelaxed(:,:,2,2)=reshape([0.0_dp,0.3_dp,-0.1_dp,0.0_dp],[2,2])
  expected_one=reshape([-1.55_dp,1.055_dp,0.34_dp,1.10_dp],[2,2])

  call assemble_mrsf_one_e_primitive_rows(ds,dh,drelaxed,dw,one_e,status)
  if(status/=MRSF_ROWS_SUCCESS) &
    error stop 'one-electron primitive rejected valid spin densities'
  if(maxval(abs(one_e-expected_one))>1.0e-14_dp) &
    error stop 'one-electron contraction has an incorrect sign or spin sum'

  two_response=reshape([1.0_dp,3.0_dp,2.0_dp,4.0_dp],[2,2])
  two_reference=reshape([0.1_dp,0.3_dp,0.2_dp,0.4_dp],[2,2])
  xc=reshape([-0.5_dp,0.7_dp,0.6_dp,-0.8_dp],[2,2])
  call combine_mrsf_response_row_blocks(one_e,two_response,two_reference, &
    xc,.true.,.true.,rows,rows_one,rows_two,rows_xc,status)
  if(status/=MRSF_ROWS_SUCCESS) &
    error stop 'complete DFT row blocks were rejected'
  if(maxval(abs(rows_one-one_e))>1.0e-14_dp) &
    error stop 'one-electron block changed during assembly'
  if(maxval(abs(rows_two-two_response-two_reference))>1.0e-14_dp) &
    error stop 'two-electron mixed block was not added exactly once'
  if(maxval(abs(rows_xc-xc))>1.0e-14_dp) &
    error stop 'complete XC block was not retained'
  if(maxval(abs(rows-rows_one-rows_two-rows_xc))>1.0e-14_dp) &
    error stop 'raw response rows do not equal their physical blocks'
  if(abs(rows(1,2)-rows(2,1))<0.1_dp) &
    error stop 'raw Cartesian response rows were silently symmetrized'

  call combine_mrsf_response_row_blocks(one_e,two_response,two_reference, &
    xc,.true.,.false.,rows,rows_one,rows_two,rows_xc,status)
  if(status/=MRSF_ROWS_XC_INCOMPLETE) &
    error stop 'an incomplete DFT XC row did not fail closed'
  if(maxval(abs(rows))>0.0_dp) &
    error stop 'failed assembly returned partial rows'

  seven_derivative=0.0_dp
  td_abxc_derivative=0.0_dp
  seven_derivative(7,1,2,1)=0.25_dp
  td_abxc_derivative(1,2,1)=0.25_dp
  call check_mrsf_seven_density_alias(seven_derivative, &
    td_abxc_derivative,1.0e-14_dp,status)
  if(status/=MRSF_ROWS_SUCCESS) &
    error stop 'consistent channel-7 and td_abxc derivatives were rejected'
  td_abxc_derivative(1,2,1)=0.24_dp
  call check_mrsf_seven_density_alias(seven_derivative, &
    td_abxc_derivative,1.0e-14_dp,status)
  if(status/=MRSF_ROWS_INCONSISTENT_BALL) &
    error stop 'inconsistent channel-7 and td_abxc derivatives were accepted'

end program test_tdhf_mrsf_hessian_rows
