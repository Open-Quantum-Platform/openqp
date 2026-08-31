module tdhf_hessian_z_rhs_mod

  use precision, only: dp

  implicit none
  private
  public :: build_tdhf_z_rhs_derivative
  public :: differentiated_channel
  public :: explicit_channel_derivative_matrix
  logical, parameter :: enable_tddft_explicit_gxc = .true.

contains

  subroutine explicit_channel_derivative_matrix(infos, coeff, base, channel, result)
    use types, only: information
    use basis_tools, only: basis_set
    use grd2, only: grd2_driver
    use tdhf_gradient_mod, only: grd2_tdhf_compute_data_t
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A
    use mathlib, only: unpack_matrix
    use mod_dft, only: dft_initialize, dftclean
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_tdxc_grad, only: tddft_xc_gradient
    type(information),target,intent(inout)::infos
    real(dp),intent(in)::coeff(:,:),base(:,:)
    integer,intent(in)::channel
    real(dp),intent(out)::result(:,:,:)
    type(basis_set),pointer::basis
    type(grd2_tdhf_compute_data_t)::gc
    type(dft_grid_t)::grid
    real(dp),contiguous,pointer::dpk(:)
    real(dp),allocatable,target::d(:,:,:),p(:,:,:),xp(:,:,:),xm(:,:,:),dxc(:,:)
    real(dp),allocatable::probe(:,:),buse(:,:),quse(:,:),gp(:,:),gm(:,:),gxcp(:,:),gxcm(:,:)
    integer::i,j,nbf,ncart
    real(dp)::scale_exch
    basis=>infos%basis; basis%atoms=>infos%atoms
    nbf=size(coeff,1); ncart=3*size(basis%atoms%xyz,2)
    scale_exch=1.0_dp
    if(infos%control%hamilton>=20) scale_exch=infos%dft%hfscale
    allocate(d(nbf,nbf,1),p(nbf,nbf,1),xp(nbf,nbf,1),xm(nbf,nbf,1),dxc(nbf,nbf), &
      probe(nbf,nbf),buse(nbf,nbf),quse(nbf,nbf),gp(3,ncart/3),gm(3,ncart/3), &
      gxcp(3,ncart/3),gxcm(3,ncart/3),source=0.0_dp)
    if(infos%control%hamilton==20 .and. channel>0 .and. enable_tddft_explicit_gxc) then
      call tagarray_get_data(infos%dat,OQP_DM_A,dpk); call unpack_matrix(dpk,dxc)
      call dft_initialize(infos,basis,grid)
    end if
    result=0.0_dp
    if(channel>0) then
      buse=0.5_dp*(base+transpose(base))
    else
      buse=base
    end if
    do j=1,nbf; do i=1,nbf
      probe=spread(coeff(:,i),2,nbf)*spread(coeff(:,j),1,nbf)
      if(channel>0) then
        quse=0.5_dp*(probe+transpose(probe))
      else
        quse=probe
      end if
      xp=0.0_dp; xm=0.0_dp
      if(channel>0) then; xp(:,:,1)=buse+quse; else; xm(:,:,1)=buse+quse; end if
      gp=0.0_dp; gc=grd2_tdhf_compute_data_t(d2=d,p2=p,xpy2=xp,xmy2=xm,hfscale=scale_exch,nbf=nbf)
      call gc%init(); call gc%build_cart(basis); call grd2_driver(infos,basis,gp,gc); call gc%clean()
      if(infos%control%hamilton==20 .and. channel>0 .and. enable_tddft_explicit_gxc) then
        gxcp=0.0_dp
        call tddft_xc_gradient(basis,grid,gxcp,dxc,p,xp,1,1.0e-14_dp,infos)
      end if
      if(channel>0) then; xp(:,:,1)=buse-quse; else; xm(:,:,1)=buse-quse; end if
      gm=0.0_dp; gc=grd2_tdhf_compute_data_t(d2=d,p2=p,xpy2=xp,xmy2=xm,hfscale=scale_exch,nbf=nbf)
      call gc%init(); call gc%build_cart(basis); call grd2_driver(infos,basis,gm,gc); call gc%clean()
      if(infos%control%hamilton==20 .and. channel>0 .and. enable_tddft_explicit_gxc) then
        gxcm=0.0_dp
        call tddft_xc_gradient(basis,grid,gxcm,dxc,p,xp,1,1.0e-14_dp,infos)
        gp=gp+gxcp; gm=gm+gxcm
      end if
      result(i,j,:)=reshape(0.25_dp*(gp-gm),[ncart])
    end do; end do
    if(infos%control%hamilton==20 .and. channel>0 .and. enable_tddft_explicit_gxc) &
      call dftclean(infos)
    deallocate(d,p,xp,xm,dxc,probe,buse,quse,gp,gm,gxcp,gxcm)
  end subroutine explicit_channel_derivative_matrix

  subroutine build_tdhf_z_rhs_derivative(infos, umat, u0, v0, du, dv, drhs)
    ! Nuclear derivative of the closed-shell TDHF Z-vector right-hand side.
    ! In the OpenQP occupied-virtual ordering the undifferentiated expression is
    !
    ! R = -( U H+_vv - H+_oo U + V H-_vv - H-_oo V + H+[T]_ov ),
    !
    ! with T_vv=(U^T U+V^T V)/2 and
    ! T_oo=-(U U^T+V V^T)/2.  Every product below is differentiated, including
    ! the integral skeleton, the response density, and both MO transformations.

    use types, only: information
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A
    use messages, only: show_message, WITH_ABORT
    use tdhf_hessian_gxc_derivative_mod, only: build_gxc_and_derivative

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: umat(:,:,:), u0(:), v0(:), du(:,:), dv(:,:)
    real(kind=dp), intent(out) :: drhs(:,:)

    real(kind=dp), contiguous, pointer :: mo(:,:)
    real(kind=dp), allocatable :: um(:,:), vm(:,:), dum(:,:), dvm(:,:)
    real(kind=dp), allocatable :: tm(:,:), dtm(:,:,:)
    real(kind=dp), allocatable :: hp(:,:), hm(:,:), ht(:,:)
    real(kind=dp), allocatable :: dhp(:,:,:), dhm(:,:,:), dht(:,:,:)
    real(kind=dp), allocatable :: gxp(:,:), dgxp(:,:,:)
    real(kind=dp), allocatable :: block(:,:)
    integer :: k, nbf, ncoord, nocc, nvir, nexc

    nbf = infos%basis%nbf
    nocc = infos%mol_prop%nocc
    nvir = nbf-nocc
    nexc = nocc*nvir
    ncoord = size(umat,3)
    if (any(shape(umat) /= [nbf,nbf,ncoord]) .or. size(u0) /= nexc .or. &
        size(v0) /= nexc .or. any(shape(du) /= [nexc,ncoord]) .or. &
        any(shape(dv) /= [nexc,ncoord]) .or. &
        any(shape(drhs) /= [nexc,ncoord])) then
      call show_message('TDHF Z-RHS derivative arrays have incompatible shapes.', WITH_ABORT)
    end if

    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    allocate(um(nocc,nvir), vm(nocc,nvir), dum(nocc,nvir), dvm(nocc,nvir))
    allocate(tm(nbf,nbf), dtm(nbf,nbf,ncoord), source=0.0_dp)
    um = reshape(u0, [nocc,nvir])
    vm = reshape(v0, [nocc,nvir])
    call accumulate_transition_square(tm, um, um, 0.5_dp)
    call accumulate_transition_square(tm, vm, vm, 0.5_dp)
    do k = 1, ncoord
      dum = reshape(du(:,k), [nocc,nvir])
      dvm = reshape(dv(:,k), [nocc,nvir])
      call accumulate_transition_square(dtm(:,:,k), dum, um, 1.0_dp)
      call accumulate_transition_square(dtm(:,:,k), dvm, vm, 1.0_dp)
    end do

    allocate(hp(nbf,nbf), hm(nbf,nbf), ht(nbf,nbf))
    allocate(dhp(nbf,nbf,ncoord), dhm(nbf,nbf,ncoord), &
             dht(nbf,nbf,ncoord), block(nocc,nvir))
    call differentiated_channel(infos, mo, umat, ov_matrix(um,nbf,nocc), &
                                ov_derivatives(du,nbf,nocc), +1, hp, dhp)
    call differentiated_channel(infos, mo, umat, ov_matrix(vm,nbf,nocc), &
                                ov_derivatives(dv,nbf,nocc), -1, hm, dhm)
    call differentiated_channel(infos, mo, umat, tm, dtm, +1, ht, dht)
    allocate(gxp(nbf,nbf), dgxp(nbf,nbf,ncoord), source=0.0_dp)
    if (infos%control%hamilton == 20) &
      call build_gxc_and_derivative(infos, mo, umat, um, du, gxp, dgxp)

    do k = 1, ncoord
      dum = reshape(du(:,k), [nocc,nvir])
      dvm = reshape(dv(:,k), [nocc,nvir])
      block = matmul(dum, transpose(hp(nocc+1:,nocc+1:))) &
            + matmul(um, transpose(dhp(nocc+1:,nocc+1:,k))) &
            - matmul(transpose(dhp(1:nocc,1:nocc,k)), um) &
            - matmul(transpose(hp(1:nocc,1:nocc)), dum) &
            + matmul(dvm, transpose(hm(nocc+1:,nocc+1:))) &
            + matmul(vm, transpose(dhm(nocc+1:,nocc+1:,k))) &
            - matmul(transpose(dhm(1:nocc,1:nocc,k)), vm) &
            - matmul(transpose(hm(1:nocc,1:nocc)), dvm) &
            + dht(1:nocc,nocc+1:,k)
      drhs(:,k) = -reshape(block, [nexc])
      ! The static Z-vector RHS contains -2 Gxc[X+Y,X+Y]_ov after the
      ! enclosing sign convention is applied.  Differentiate that term too.
      drhs(:,k) = drhs(:,k) &
        - 2.0_dp*reshape(dgxp(1:nocc,nocc+1:,k), [nexc])
    end do

    deallocate(um, vm, dum, dvm, tm, dtm, hp, hm, ht, dhp, dhm, dht, &
               gxp, dgxp, block)
  end subroutine build_tdhf_z_rhs_derivative

!###############################################################################

  subroutine differentiated_channel(infos, coeff, umat, m0, dm, sign_channel, gmo, dgmo)
    use types, only: information
    use basis_tools, only: basis_set

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: coeff(:,:), umat(:,:,:), m0(:,:), dm(:,:,:)
    integer, intent(in) :: sign_channel
    real(kind=dp), intent(out) :: gmo(:,:), dgmo(:,:,:)

    type(basis_set), pointer :: basis
    real(kind=dp), allocatable, target :: p0(:,:), pt(:,:)
    real(kind=dp), allocatable :: dpao(:,:,:), dcoeff(:,:), work(:,:)
    real(kind=dp), allocatable :: deri(:,:,:)
    real(kind=dp), allocatable :: dens(:,:,:), amb(:,:,:), apb(:,:,:)
    real(kind=dp), allocatable :: gao(:,:), dgao(:,:)
    integer :: atom, cart, k, nbf, ncoord

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = size(coeff,1)
    ncoord = size(umat,3)
    allocate(p0(nbf,nbf), pt(nbf,nbf), dpao(nbf,nbf,ncoord), &
             dcoeff(nbf,nbf), work(nbf,nbf))
    p0 = matmul(matmul(coeff,m0),transpose(coeff))
    if (sign_channel > 0) p0 = 0.5_dp*(p0+transpose(p0))
    pt = transpose(p0)
    do k = 1, ncoord
      dcoeff = matmul(coeff,umat(:,:,k))
      dpao(:,:,k) = matmul(matmul(dcoeff,m0),transpose(coeff)) &
                    + matmul(matmul(coeff,dm(:,:,k)),transpose(coeff)) &
                    + matmul(matmul(coeff,m0),transpose(dcoeff))
      if (sign_channel > 0) dpao(:,:,k) = 0.5_dp*(dpao(:,:,k)+transpose(dpao(:,:,k)))
    end do

    allocate(deri(nbf,nbf,ncoord))
    call explicit_channel_derivative_matrix(infos,coeff,p0,sign_channel,deri)
    allocate(dens(nbf,nbf,ncoord+1), amb(nbf,nbf,ncoord+1), &
             apb(nbf,nbf,ncoord+1))
    dens(:,:,1) = p0
    dens(:,:,2:) = dpao
    if (sign_channel < 0) then
      call apply_gradient_channel(infos,dens,sign_channel,amb)
      apb=0.0_dp
    else
      call apply_gradient_channel(infos,dens,sign_channel,apb)
      amb=0.0_dp
    end if
    allocate(gao(nbf,nbf),dgao(nbf,nbf))
    if (sign_channel < 0) then
      gao = amb(:,:,1)
    else
      gao = apb(:,:,1)
    end if
    call ao_to_mo(gao,coeff,gmo,work)
    do k = 1, ncoord
      if (sign_channel < 0) then
        dgao = amb(:,:,k+1)
      else
        dgao = apb(:,:,k+1)
      end if
      call ao_to_mo(dgao,coeff,dgmo(:,:,k),work)
      dgmo(:,:,k) = dgmo(:,:,k) + deri(:,:,k) + matmul(transpose(umat(:,:,k)),gmo) &
                                     + matmul(gmo,umat(:,:,k))
    end do
    deallocate(p0,pt,dpao,dcoeff,work,deri,dens,amb,apb,gao,dgao)
  end subroutine differentiated_channel

!###############################################################################

  subroutine apply_gradient_channel(infos,density,sign_channel,result)
    ! Full AO H+/H- matrices in exactly the representation used by the
    ! production TDHF gradient.  The Davidson-oriented int2_td_data_t agrees
    ! in the occupied-virtual action but not in all OO/VV blocks, which are
    ! required by the differentiated Z-vector and W density.
    use types, only: information
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_td_data_t
    use tdhf_response_operator_mod, only: apply_tdhf_ao_operators
    type(information),target,intent(inout)::infos
    real(dp),intent(in),target::density(:,:,:)
    integer,intent(in)::sign_channel
    real(dp),intent(out)::result(:,:,:)
    type(int2_compute_t),target::driver
    type(int2_td_data_t),target::data
    real(dp),allocatable::amb(:,:,:),apb(:,:,:)
    if(infos%control%hamilton==20) then
      allocate(amb(size(result,1),size(result,2),size(result,3)), &
               apb(size(result,1),size(result,2),size(result,3)))
      call apply_tdhf_ao_operators(infos,density,amb,apb)
      if(sign_channel>0) then
        result=apb
      else
        result=amb
      end if
      deallocate(amb,apb)
      return
    end if
    call driver%init(infos%basis,infos); call driver%set_screening()
    if(sign_channel>0) then
      data=int2_td_data_t(d2=density,int_apb=.true.,int_amb=.false., &
                         tamm_dancoff=.false.,scale_exchange=1.0_dp)
      call driver%run(data)
      result=data%apb(:,:,:,1)
    else
      data=int2_td_data_t(d2=density,int_apb=.false.,int_amb=.true., &
                         tamm_dancoff=.false.,scale_exchange=1.0_dp)
      call driver%run(data)
      result=data%amb(:,:,:,1)
    end if
    call data%clean(); call driver%clean()
  end subroutine apply_gradient_channel

!###############################################################################

  pure subroutine accumulate_transition_square(t, a, b, factor)
    real(kind=dp), intent(inout) :: t(:,:)
    real(kind=dp), intent(in) :: a(:,:), b(:,:), factor
    integer :: no
    no = size(a,1)
    t(1:no,1:no) = t(1:no,1:no) &
      - 0.5_dp*factor*(matmul(a,transpose(b))+matmul(b,transpose(a)))
    t(no+1:,no+1:) = t(no+1:,no+1:) &
      + 0.5_dp*factor*(matmul(transpose(a),b)+matmul(transpose(b),a))
  end subroutine accumulate_transition_square

  pure function ov_matrix(z, nbf, nocc) result(m)
    real(kind=dp), intent(in) :: z(:,:)
    integer, intent(in) :: nbf, nocc
    real(kind=dp) :: m(nbf,nbf)
    m = 0.0_dp
    m(1:nocc,nocc+1:) = z
  end function ov_matrix

  pure function ov_derivatives(dz, nbf, nocc) result(dm)
    real(kind=dp), intent(in) :: dz(:,:)
    integer, intent(in) :: nbf, nocc
    real(kind=dp) :: dm(nbf,nbf,size(dz,2))
    integer :: k
    dm = 0.0_dp
    do k = 1, size(dz,2)
      dm(1:nocc,nocc+1:,k) = reshape(dz(:,k),[nocc,nbf-nocc])
    end do
  end function ov_derivatives

  subroutine ao_to_mo(ao, coeff, mo_mat, work)
    real(kind=dp), intent(in) :: ao(:,:), coeff(:,:)
    real(kind=dp), intent(out) :: mo_mat(:,:), work(:,:)
    work = matmul(ao,coeff)
    mo_mat = matmul(transpose(coeff),work)
  end subroutine ao_to_mo

end module tdhf_hessian_z_rhs_mod
