module tdhf_hessian_gxc_derivative_mod

  use precision, only: dp
  implicit none
  private
  public :: build_gxc_and_derivative

contains

  ! Gxc[X+Y,X+Y] and its total nuclear derivative in the moving MO frame.
  ! The central geometry evaluation includes grid/basis motion, while c*U and
  ! d(X+Y) carry the orbital and amplitude response required by GAMESS TDHXG1.
  subroutine build_gxc_and_derivative(infos,c,umat,xpy,dxpy,gxp,dgxp)
    use types, only: information
    use basis_tools, only: basis_set
    use mod_dft, only: dft_initialize,dftclean
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_gxc, only: tddft_gxc
    use mathlib, only: orthogonal_transform
    type(information),target,intent(inout)::infos
    real(dp),intent(in)::c(:,:),umat(:,:,:),xpy(:,:),dxpy(:,:)
    real(dp),intent(out)::gxp(:,:),dgxp(:,:,:)
    type(basis_set),pointer::basis
    type(dft_grid_t)::grid
    real(dp),allocatable,target::xao(:,:,:)
    real(dp),allocatable::cp(:,:),cm(:,:),dcp(:,:),mp(:,:),mm(:,:),gp(:,:,:),gm(:,:,:)
    real(dp),parameter::step=1.0e-3_dp
    integer::k,cc,aa,nbf,no,ncart

    basis=>infos%basis; basis%atoms=>infos%atoms
    nbf=size(c,1); no=size(xpy,1); ncart=size(umat,3)
    allocate(xao(nbf,nbf,1),cp(nbf,nbf),cm(nbf,nbf),dcp(nbf,nbf), &
      mp(nbf,nbf),mm(nbf,nbf),gp(nbf,nbf,1),gm(nbf,nbf,1),source=0.0_dp)

    call symmetric_transition(xpy,mp)
    xao(:,:,1)=matmul(matmul(c,mp),transpose(c))
    gxp=0.0_dp
    call dft_initialize(infos,basis,grid)
    gp=0.0_dp
    call tddft_gxc(basis,grid,.true.,c,gp,xao,1,0.0_dp,infos)
    call dftclean(infos)
    gxp=gp(:,:,1)
    call orthogonal_transform('n',nbf,c,gxp)

    do k=1,ncart
      dcp=matmul(c,umat(:,:,k)); cp=c+step*dcp; cm=c-step*dcp
      call symmetric_transition(xpy+step*reshape(dxpy(:,k),shape(xpy)),mp)
      call symmetric_transition(xpy-step*reshape(dxpy(:,k),shape(xpy)),mm)
      xao(:,:,1)=matmul(matmul(cp,mp),transpose(cp))
      cc=mod(k-1,3)+1; aa=(k-1)/3+1
      basis%atoms%xyz(cc,aa)=basis%atoms%xyz(cc,aa)+step
      call basis%init_shell_centers()
      call dft_initialize(infos,basis,grid)
      gp=0.0_dp
      call tddft_gxc(basis,grid,.true.,cp,gp,xao,1,0.0_dp,infos)
      call dftclean(infos)
      call orthogonal_transform('n',nbf,cp,gp(:,:,1))

      xao(:,:,1)=matmul(matmul(cm,mm),transpose(cm))
      basis%atoms%xyz(cc,aa)=basis%atoms%xyz(cc,aa)-2.0_dp*step
      call basis%init_shell_centers()
      call dft_initialize(infos,basis,grid)
      gm=0.0_dp
      call tddft_gxc(basis,grid,.true.,cm,gm,xao,1,0.0_dp,infos)
      call dftclean(infos)
      call orthogonal_transform('n',nbf,cm,gm(:,:,1))
      basis%atoms%xyz(cc,aa)=basis%atoms%xyz(cc,aa)+step
      call basis%init_shell_centers()
      dgxp(:,:,k)=(gp(:,:,1)-gm(:,:,1))/(2.0_dp*step)
    end do
    deallocate(xao,cp,cm,dcp,mp,mm,gp,gm)
  contains
    subroutine symmetric_transition(a,m)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(out)::m(:,:)
      m=0.0_dp
      m(1:no,no+1:)=0.5_dp*a
      m(no+1:,1:no)=0.5_dp*transpose(a)
    end subroutine symmetric_transition
  end subroutine build_gxc_and_derivative

end module tdhf_hessian_gxc_derivative_mod
