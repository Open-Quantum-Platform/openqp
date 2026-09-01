module tdhf_hessian_z_rhs_mod

  use precision, only: dp
  use grd2_rys, only: grd2_operator_consumer_t

  implicit none
  private
  public :: build_tdhf_z_rhs_derivative
  public :: differentiated_channel
  public :: explicit_channel_derivative_matrix
  public :: accumulate_tdhf_channel_quartet
  logical, parameter :: enable_tddft_explicit_gxc = .true.

  type, extends(grd2_operator_consumer_t) :: tdhf_channel_operator_consumer_t
    real(dp), pointer :: base(:,:) => null()
    real(dp), pointer :: operator(:,:,:) => null()
    integer, allocatable :: cart_off(:)
    integer :: channel = 0
    real(dp) :: coulscale = 1.0_dp
    real(dp) :: hfscale = 1.0_dp
  contains
    procedure :: accumulate => accumulate_tdhf_channel_operator
  end type tdhf_channel_operator_consumer_t

contains

  subroutine explicit_channel_derivative_matrix(infos, coeff, base, channel, result)
    use types, only: information
    use basis_tools, only: basis_set, bas_norm_matrix, build_cart_density
    use grd2, only: grd2_operator_driver
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
    type(tdhf_channel_operator_consumer_t)::consumer
    type(dft_grid_t)::grid
    real(dp),contiguous,pointer::dpk(:)
    real(dp),allocatable,target::p(:,:,:),xp(:,:,:),dxc(:,:)
    real(dp),allocatable,target::bwork(:,:),base_cart(:,:),operator_cart(:,:,:)
    real(dp),allocatable::probe(:,:),buse(:,:),quse(:,:),operator_ao(:,:,:), &
      work(:,:),gp(:,:),gm(:,:)
    integer::i,j,k,nbf,ncart,nwork
    real(dp)::scale_exch
    basis=>infos%basis; basis%atoms=>infos%atoms
    nbf=size(coeff,1); ncart=3*size(basis%atoms%xyz,2)
    scale_exch=1.0_dp
    if(infos%control%hamilton>=20) scale_exch=infos%dft%hfscale
    allocate(buse(nbf,nbf),source=base)
    if(channel>0) then
      buse=0.5_dp*(base+transpose(base))
    end if

    ! Match the exact density convention consumed by grd2 under pure
    ! spherical harmonics: bfnrm folding followed by blockwise expansion to
    ! Cartesian effective densities.
    bwork=buse
    call bas_norm_matrix(bwork,basis%bfnrm,nbf)
    call build_cart_density(basis,bwork,base_cart,consumer%cart_off,nwork)
    consumer%base=>base_cart
    allocate(operator_cart(nwork,nwork,ncart),source=0.0_dp)
    consumer%operator=>operator_cart
    consumer%channel=channel
    consumer%coulscale=1.0_dp
    consumer%hfscale=scale_exch
    call grd2_operator_driver(infos,basis,consumer)

    do k=1,ncart
      if(channel>0) then
        operator_cart(:,:,k)=0.5_dp*(operator_cart(:,:,k)+transpose(operator_cart(:,:,k)))
      else
        operator_cart(:,:,k)=0.5_dp*(operator_cart(:,:,k)-transpose(operator_cart(:,:,k)))
      end if
    end do
    allocate(operator_ao(nbf,nbf,ncart))
    call reduce_cartesian_operator(basis,consumer%cart_off,operator_cart,operator_ao)
    allocate(work(nbf,nbf))
    do k=1,ncart
      work=matmul(operator_ao(:,:,k),coeff)
      result(:,:,k)=matmul(transpose(coeff),work)
    end do

    ! The ERI term above is now one blocked quartet traversal.  The XC grid
    ! contribution uses its existing finite-difference oracle independently;
    ! it does not invoke grd2 and therefore cannot reintroduce derivative-ERI
    ! traversals.
    if(infos%control%hamilton==20 .and. channel>0 .and. enable_tddft_explicit_gxc) then
      allocate(p(nbf,nbf,1),xp(nbf,nbf,1),dxc(nbf,nbf),probe(nbf,nbf), &
        quse(nbf,nbf),gp(3,ncart/3),gm(3,ncart/3),source=0.0_dp)
      call tagarray_get_data(infos%dat,OQP_DM_A,dpk); call unpack_matrix(dpk,dxc)
      call dft_initialize(infos,basis,grid)
      do j=1,nbf; do i=1,nbf
        probe=spread(coeff(:,i),2,nbf)*spread(coeff(:,j),1,nbf)
        quse=0.5_dp*(probe+transpose(probe))
        xp(:,:,1)=buse+quse; gp=0.0_dp
        call tddft_xc_gradient(basis,grid,gp,dxc,p,xp,1,1.0e-14_dp,infos)
        xp(:,:,1)=buse-quse; gm=0.0_dp
        call tddft_xc_gradient(basis,grid,gm,dxc,p,xp,1,1.0e-14_dp,infos)
        result(i,j,:)=result(i,j,:)+reshape(0.25_dp*(gp-gm),[ncart])
      end do; end do
      call dftclean(infos)
      deallocate(p,xp,dxc,probe,quse,gp,gm)
    end if
    nullify(consumer%base,consumer%operator)
    deallocate(buse,bwork,base_cart,operator_cart,operator_ao,work)
  end subroutine explicit_channel_derivative_matrix

  subroutine accumulate_tdhf_channel_operator(this,basis,shell_ids,atom_ids, &
                                               local_ids,derivative)
    use basis_tools, only: basis_set
    class(tdhf_channel_operator_consumer_t),intent(inout)::this
    type(basis_set),intent(in)::basis
    integer,intent(in)::shell_ids(4),atom_ids(4),local_ids(4)
    real(dp),intent(in)::derivative(3,4)
    integer::gi,gj,gk,gl,center,axis,coord
    real(dp)::value

    gi=this%cart_off(shell_ids(1))+local_ids(1)-1
    gj=this%cart_off(shell_ids(2))+local_ids(2)-1
    gk=this%cart_off(shell_ids(3))+local_ids(3)-1
    gl=this%cart_off(shell_ids(4))+local_ids(4)-1
    ! A shell quartet can place more than one derivative center on the same
    ! atom.  Accumulate each center explicitly so those contributions add
    ! rather than being lost through a repeated vector subscript.
    do center=1,4
      do axis=1,3
        value=derivative(axis,center)
        if(value==0.0_dp) cycle
        coord=3*(atom_ids(center)-1)+axis
        call accumulate_tdhf_channel_quartet(this%base,this%operator(:,:,coord), &
          [gi,gj,gk,gl],this%channel,this%coulscale,this%hfscale,value)
      end do
    end do
  end subroutine accumulate_tdhf_channel_operator

  pure subroutine accumulate_tdhf_channel_quartet(base,operator,ids,channel, &
                                                   coulscale,hfscale,value)
    real(dp),intent(in)::base(:,:),coulscale,hfscale,value
    real(dp),intent(inout)::operator(:,:)
    integer,intent(in)::ids(4),channel
    integer::gi,gj,gk,gl
    real(dp)::ail,ajk,ajl,aik
    gi=ids(1);gj=ids(2);gk=ids(3);gl=ids(4)
    if(channel>0) then
      operator(gk,gl)=operator(gk,gl)+16.0_dp*coulscale*base(gi,gj)*value
      operator(gi,gj)=operator(gi,gj)+16.0_dp*coulscale*base(gk,gl)*value
      operator(gl,gj)=operator(gl,gj)-4.0_dp*hfscale*base(gk,gi)*value
      operator(gk,gi)=operator(gk,gi)-4.0_dp*hfscale*base(gl,gj)*value
      operator(gk,gj)=operator(gk,gj)-4.0_dp*hfscale*base(gl,gi)*value
      operator(gl,gi)=operator(gl,gi)-4.0_dp*hfscale*base(gk,gj)*value
    else
      ail=base(gi,gl)-base(gl,gi)
      ajk=base(gj,gk)-base(gk,gj)
      ajl=base(gj,gl)-base(gl,gj)
      aik=base(gi,gk)-base(gk,gi)
      call add_antisymmetric(operator,gj,gk,-hfscale*ail*value)
      call add_antisymmetric(operator,gi,gl,-hfscale*ajk*value)
      call add_antisymmetric(operator,gi,gk,-hfscale*ajl*value)
      call add_antisymmetric(operator,gj,gl,-hfscale*aik*value)
    end if
  end subroutine accumulate_tdhf_channel_quartet

  pure subroutine add_antisymmetric(matrix,row,col,value)
    real(dp),intent(inout)::matrix(:,:)
    integer,intent(in)::row,col
    real(dp),intent(in)::value
    matrix(row,col)=matrix(row,col)+value
    matrix(col,row)=matrix(col,row)-value
  end subroutine add_antisymmetric

  subroutine reduce_cartesian_operator(basis,cart_off,cart,sph)
    use basis_tools, only: basis_set,bas_norm_matrix
    use cart2sph, only: cart2sph_mat
    use constants, only: NUM_CART_BF
    type(basis_set),intent(in)::basis
    integer,intent(in)::cart_off(:)
    real(dp),intent(in)::cart(:,:,:)
    real(dp),intent(out)::sph(:,:,:)
    real(dp),allocatable::block(:)
    integer::coord,si,sj,ii,jj,nci,ncj,nsi,nsj,coi,coj,soi,soj,idx

    sph=0.0_dp
    do coord=1,size(cart,3)
      do sj=1,basis%nshell
        ncj=NUM_CART_BF(basis%am(sj));nsj=basis%naos(sj)
        coj=cart_off(sj);soj=basis%ao_offset(sj)
        do si=1,basis%nshell
          nci=NUM_CART_BF(basis%am(si));nsi=basis%naos(si)
          coi=cart_off(si);soi=basis%ao_offset(si)
          allocate(block(nci*ncj));idx=0
          do jj=1,ncj;do ii=1,nci
            idx=idx+1
            block(idx)=cart(coi+ii-1,coj+jj-1,coord)
          end do;end do
          if(basis%harmonic(si)==1 .or. basis%harmonic(sj)==1) &
            call cart2sph_mat(block,basis%am(si),basis%harmonic(si), &
              basis%am(sj),basis%harmonic(sj))
          idx=0
          do jj=1,nsj;do ii=1,nsi
            idx=idx+1
            sph(soi+ii-1,soj+jj-1,coord)=block(idx)
          end do;end do
          deallocate(block)
        end do
      end do
      call bas_norm_matrix(sph(:,:,coord),basis%bfnrm,basis%nbf)
    end do
  end subroutine reduce_cartesian_operator

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
    real(kind=dp), allocatable :: gao(:,:), dgao(:,:), kxc_ground(:,:,:)
    integer :: atom, cart, k, nbf, ncoord, nocc

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = size(coeff,1)
    ncoord = size(umat,3)
    nocc = infos%mol_prop%nocc
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
    allocate(kxc_ground(nbf,nbf,ncoord),source=0.0_dp)
    if (infos%control%hamilton==20 .and. sign_channel>0) then
      block
        use mod_dft, only: dft_initialize,dftclean
        use mod_dft_molgrid, only: dft_grid_t
        use mod_dft_gridint_gxc, only: tddft_gxc
        type(dft_grid_t) :: grid
        real(dp),allocatable,target :: dx(:,:,:)
        real(dp),allocatable :: fx(:,:,:),dground(:,:)
        integer :: kk
        allocate(dx(nbf,nbf,3*ncoord),fx(nbf,nbf,3*ncoord),dground(nbf,nbf),source=0.0_dp)
        do kk=1,ncoord
          dcoeff=matmul(coeff,umat(:,:,kk))
          dground=2.0_dp*(matmul(dcoeff(:,1:nocc),transpose(coeff(:,1:nocc)))+ &
            matmul(coeff(:,1:nocc),transpose(dcoeff(:,1:nocc))))
          ! The restricted grid interface consumes a one-spin ground density.
          dx(:,:,3*kk-2)=p0+0.5_dp*dground
          dx(:,:,3*kk-1)=p0
          dx(:,:,3*kk)=0.5_dp*dground
        end do
        call dft_initialize(infos,basis,grid)
        call tddft_gxc(basis,grid,.true.,coeff,fx,dx,3*ncoord,1.0e-14_dp,infos)
        call dftclean(infos)
        do kk=1,ncoord
          kxc_ground(:,:,kk)=fx(:,:,3*kk-2)-fx(:,:,3*kk-1)-fx(:,:,3*kk)
        end do
        deallocate(dx,fx,dground)
      end block
    end if
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
        dgao = apb(:,:,k+1)+kxc_ground(:,:,k)
      end if
      call ao_to_mo(dgao,coeff,dgmo(:,:,k),work)
      dgmo(:,:,k) = dgmo(:,:,k) + deri(:,:,k) + matmul(transpose(umat(:,:,k)),gmo) &
                                     + matmul(gmo,umat(:,:,k))
    end do
    deallocate(p0,pt,dpao,dcoeff,work,deri,dens,amb,apb,gao,dgao,kxc_ground)
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
