! Compare cache reuse with forced LibXC initialization on the same SCF densities.
! No information/functional_t copy is made: their LibXC handles have one owner.
subroutine functional_owner_cache_compare(a, b, error, failed) bind(C)
  use iso_c_binding, only: c_double, c_int, c_null_char
  use precision, only: dp, i8b
  use c_interop, only: oqp_handle_t, oqp_handle_get_info
  use types, only: information, dft_parameters, tddft_parameters
  use basis_tools, only: basis_set
  use mod_dft, only: dft_initialize, dftclean, dft_forget_setups
  use mod_dft_molgrid, only: dft_grid_t
  use mod_dft_gridint_energy, only: dmatd_density_blk
  use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A
  use mathlib, only: unpack_matrix
  implicit none
  type(oqp_handle_t), intent(in) :: a, b
  real(c_double), intent(out) :: error
  integer(c_int), intent(out) :: failed
  type(information), pointer :: ia, ib, infos
  type(dft_parameters) :: saved_a, saved_b
  type(tddft_parameters) :: saved_ta, saved_tb
  type(dft_parameters) :: cached_dft
  type(tddft_parameters) :: cached_tddft
  type(basis_set) :: basis
  type(dft_grid_t) :: grid
  real(dp), pointer :: packed(:)
  real(dp), allocatable :: da(:,:), db(:,:), fa(:), fb(:), ra(:), rb(:)
  real(dp) :: energy, electrons, kinetic, reference
  integer(i8b) :: revision, previous_revision
  integer :: step, n, nt, i
  integer(c_int) :: status
  character(len=20) :: name

  ia => oqp_handle_get_info(a)
  ib => oqp_handle_get_info(b)
  saved_a=ia%dft; saved_b=ib%dft
  saved_ta=ia%tddft; saved_tb=ib%tddft
  error=0; failed=0
  call dftclean(ia); call dftclean(ib)
  do step=1,27
    infos=>ia
    if (step==2.or.step==12) infos=>ib
    name='BHHLYP'
    if (step==1) name='SLATER'
    if (step>=6.and.step<=10) name='CAM-B3LYP'
    if (step==17.or.step==18) name='STG1X'
    if (step>=19) name='DTCAMVAEE'
    infos%dft%xc_functional_name=c_null_char
    do i=1,len_trim(name)
      infos%dft%xc_functional_name(i)=name(i:i)
    end do
    infos%dft%grid_pruned=.false.
    infos%dft%grid_rad_size=24
    infos%dft%grid_ang_size=50
    if (step==12) then
      infos%dft%grid_rad_size=48
      infos%dft%grid_ang_size=110
    end if
    ! A/B/A has identical request metadata but A still owns SLATER at step 3.
    if (step<=3) then
      infos%dft%hfscale=0.5_dp
      infos%dft%cam_flag=.false.
      infos%dft%cam_alpha=0; infos%dft%cam_beta=0; infos%dft%cam_mu=0
    end if
    if (step==4) infos%dft%hfscale=0.27_dp
    ! A small but real setting change must invalidate, without a tolerance.
    if (step==5) infos%dft%hfscale=infos%dft%hfscale+1.0e-13_dp
    if (step==6) then
      infos%dft%cam_flag=.true.
      infos%dft%cam_alpha=-1; infos%dft%cam_beta=-1; infos%dft%cam_mu=-1
    end if
    if (step==7) infos%dft%cam_alpha=0.22_dp
    if (step==8) infos%dft%cam_beta=0.38_dp
    if (step==9) infos%dft%cam_mu=0.41_dp
    if (step==10) infos%dft%cam_mu=infos%dft%cam_mu+1.0e-13_dp
    if (step==11) then
      infos%dft%cam_flag=.false.
      call infos%functional%destroy()
    end if
    if (step==13) call dft_forget_setups()
    if (step==14) then
      call dftclean(infos); call dftclean(infos)
    end if
    if (step==15) infos%dft%MP2OS_Scale=0.31_dp
    if (step==16) infos%dft%dh_flag=.true.
    if (step==18) then
      infos%tddft%spc_coco=0.21_dp
      infos%tddft%spc_coov=0.23_dp
      infos%tddft%spc_ovov=0.27_dp
    end if
    if (step==20) infos%tddft%cam_alpha=0.31_dp
    if (step==21) infos%tddft%cam_beta=0.32_dp
    if (step==22) infos%tddft%cam_mu=infos%tddft%cam_mu+0.01_dp
    if (step==23) infos%dft%MP2SS_Scale=0.29_dp
    if (step==24) infos%tddft%spc_coco=0.21_dp
    if (step==25) infos%tddft%spc_coov=0.23_dp
    if (step==26) infos%tddft%spc_ovov=0.27_dp
    ! Change only this predicate; the functional and every CAM value are kept.
    if (step==27) infos%dft%cam_flag=.not.infos%dft%cam_flag
    basis=infos%basis; basis%atoms=>infos%atoms
    previous_revision=infos%functional%cache_revision
    call dft_initialize(infos,basis,grid,verbose=.false.)
    revision=infos%functional%cache_revision
    if (previous_revision==revision) failed=failed+1
    call dft_initialize(infos,basis,grid,verbose=.false.)
    if (revision/=infos%functional%cache_revision) failed=failed+1
    n=basis%nbf; nt=n*(n+1)/2
    allocate(da(n,n),db(n,n),fa(nt),fb(nt),ra(nt),rb(nt))
    call tagarray_get_data(infos%dat,OQP_DM_A,packed,status)
    if (status/=0) then
      failed=failed+1
      deallocate(da,db,fa,fb,ra,rb)
      exit
    end if
    call unpack_matrix(packed,da)
    ! Unequal spin densities retain the signs/off-diagonals of the SCF density.
    db=0.73_dp*da
    fa=0; fb=0
    call dmatd_density_blk(basis,grid,da,db,fa,fb,energy,electrons,kinetic, &
                          basis%mxam+2,n,1.0e-15_dp,.true.,infos)
    cached_dft=infos%dft; cached_tddft=infos%tddft
    call dftclean(infos)
    call dft_initialize(infos,basis,grid,verbose=.false.)
    if (cached_dft%MP2SS_Scale/=infos%dft%MP2SS_Scale.or. &
        cached_dft%MP2OS_Scale/=infos%dft%MP2OS_Scale) failed=failed+1
    if (cached_tddft%cam_alpha/=infos%tddft%cam_alpha.or. &
        cached_tddft%cam_beta/=infos%tddft%cam_beta.or. &
        cached_tddft%cam_mu/=infos%tddft%cam_mu.or. &
        cached_tddft%spc_coco/=infos%tddft%spc_coco.or. &
        cached_tddft%spc_coov/=infos%tddft%spc_coov.or. &
        cached_tddft%spc_ovov/=infos%tddft%spc_ovov) failed=failed+1
    ra=0; rb=0
    call dmatd_density_blk(basis,grid,da,db,ra,rb,reference,electrons,kinetic, &
                          basis%mxam+2,n,1.0e-15_dp,.true.,infos)
    error=max(error,abs(energy-reference),maxval(abs(fa-ra)),maxval(abs(fb-rb)))
    if (.not.all(abs(fa-ra)<=1.0e-11_dp)) failed=failed+1
    if (.not.all(abs(fb-rb)<=1.0e-11_dp)) failed=failed+1
    if (.not.(abs(energy-reference)<=1.0e-11_dp)) failed=failed+1
    deallocate(da,db,fa,fb,ra,rb)
  end do
  call dftclean(ia); call dftclean(ib)
  ia%dft=saved_a; ib%dft=saved_b
  ia%tddft=saved_ta; ib%tddft=saved_tb
end subroutine functional_owner_cache_compare
