! Fixed-grid nuclear derivatives of density data used by GGA quadratures.
module mod_dft_gga_nuclear_point

  use precision, only: fp

  implicit none
  private

  integer, parameter :: hmap(3,3) = reshape([1,4,6, 4,2,5, 6,5,3], [3,3])
  integer, parameter :: tmap(3,3,3) = reshape([ &
    1,4,5, 4,6,10, 5,10,8, &
    4,6,10, 6,2,7, 10,7,9, &
    5,10,8, 10,7,9, 8,9,3], [3,3,3])

  public :: gga_density_nuclear_point,gga_density_nuclear_point_batch
  public :: gga_density_nuclear_point_first_batch

contains

!> Compute fixed-grid first and second nuclear derivatives of rho and grad(rho).
!>
!> AO derivatives are electronic-coordinate derivatives ordered as
!> G2=(xx,yy,zz,xy,yz,xz) and
!> G3=(xxx,yyy,zzz,xxy,xxz,yyx,yyz,zzx,zzy,xyz).
!> Moving the center of AO mu gives d(phi_mu)/dR_Aa=-delta(mu,A)d_a(phi_mu).
!> Thus this routine supplies the DRA/GDA/DGGA data underlying the GAMESS
!> DDDENCNST and TDHXG1G/TDHXGPG construction, without grid-weight or
!> grid-center translation terms.
  subroutine gga_density_nuclear_point(density, ao_atom, aov, aog1, aog2, &
                                       aog3, drho, dgrho, d2rho, d2grho)
    real(fp), intent(in) :: density(:,:)
    integer, intent(in) :: ao_atom(:)
    real(fp), intent(in) :: aov(:), aog1(:,:), aog2(:,:), aog3(:,:)
    real(fp), intent(out) :: drho(:,:), dgrho(:,:,:)
    real(fp), intent(out) :: d2rho(:,:,:,:), d2grho(:,:,:,:,:)

    integer :: mu, nu, a, b, c, atom_a, atom_b, nao, nat
    real(fp) :: p, qmu, qnu, rmu, rnu, qrmu, qrnu
    real(fp) :: gc_mu, gc_nu, qgc_mu, qgc_nu, rgc_mu, rgc_nu
    real(fp) :: qrgc_mu, qrgc_nu

    nao = size(aov)
    nat = size(drho,2)
    call check_shapes(nao, nat)

    drho = 0.0_fp
    dgrho = 0.0_fp
    d2rho = 0.0_fp
    d2grho = 0.0_fp

    do nu = 1, nao
      do mu = 1, nao
        p = density(mu,nu)
        if (p == 0.0_fp) cycle
        do atom_a = 1, nat
          do a = 1, 3
            qmu = center_d1(mu,atom_a,a)
            qnu = center_d1(nu,atom_a,a)
            drho(a,atom_a) = drho(a,atom_a) &
              + p*(qmu*aov(nu) + aov(mu)*qnu)
            do c = 1, 3
              gc_mu = aog1(mu,c)
              gc_nu = aog1(nu,c)
              qgc_mu = center_gd1(mu,atom_a,a,c)
              qgc_nu = center_gd1(nu,atom_a,a,c)
              dgrho(c,a,atom_a) = dgrho(c,a,atom_a) + p*( &
                qgc_mu*aov(nu) + gc_mu*qnu + qmu*gc_nu + aov(mu)*qgc_nu)
            end do

            do atom_b = 1, nat
              do b = 1, 3
                rmu = center_d1(mu,atom_b,b)
                rnu = center_d1(nu,atom_b,b)
                qrmu = center_d2(mu,atom_a,a,atom_b,b)
                qrnu = center_d2(nu,atom_a,a,atom_b,b)
                d2rho(a,b,atom_a,atom_b) = d2rho(a,b,atom_a,atom_b) + p*( &
                  qrmu*aov(nu) + qmu*rnu + rmu*qnu + aov(mu)*qrnu)
                do c = 1, 3
                  gc_mu = aog1(mu,c)
                  gc_nu = aog1(nu,c)
                  qgc_mu = center_gd1(mu,atom_a,a,c)
                  qgc_nu = center_gd1(nu,atom_a,a,c)
                  rgc_mu = center_gd1(mu,atom_b,b,c)
                  rgc_nu = center_gd1(nu,atom_b,b,c)
                  qrgc_mu = center_gd2(mu,atom_a,a,atom_b,b,c)
                  qrgc_nu = center_gd2(nu,atom_a,a,atom_b,b,c)
                  d2grho(c,a,b,atom_a,atom_b) = &
                    d2grho(c,a,b,atom_a,atom_b) + p*( &
                    qrgc_mu*aov(nu) + qgc_mu*rnu + rgc_mu*qnu + gc_mu*qrnu &
                    + qrmu*gc_nu + qmu*rgc_nu + rmu*qgc_nu + aov(mu)*qrgc_nu)
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  contains

    pure real(fp) function center_d1(i, atom, ixyz)
      integer, intent(in) :: i, atom, ixyz
      center_d1 = 0.0_fp
      if (ao_atom(i) == atom) center_d1 = -aog1(i,ixyz)
    end function center_d1

    pure real(fp) function center_gd1(i, atom, ixyz, cxyz)
      integer, intent(in) :: i, atom, ixyz, cxyz
      center_gd1 = 0.0_fp
      if (ao_atom(i) == atom) center_gd1 = -aog2(i,hmap(ixyz,cxyz))
    end function center_gd1

    pure real(fp) function center_d2(i, atom1, xyz1, atom2, xyz2)
      integer, intent(in) :: i, atom1, xyz1, atom2, xyz2
      center_d2 = 0.0_fp
      if (ao_atom(i) == atom1 .and. atom1 == atom2) &
        center_d2 = aog2(i,hmap(xyz1,xyz2))
    end function center_d2

    pure real(fp) function center_gd2(i, atom1, xyz1, atom2, xyz2, cxyz)
      integer, intent(in) :: i, atom1, xyz1, atom2, xyz2, cxyz
      center_gd2 = 0.0_fp
      if (ao_atom(i) == atom1 .and. atom1 == atom2) &
        center_gd2 = aog3(i,tmap(xyz1,xyz2,cxyz))
    end function center_gd2

    subroutine check_shapes(n, n_atom)
      integer, intent(in) :: n, n_atom
      if (size(density,1) /= n .or. size(density,2) /= n) &
        error stop 'gga_density_nuclear_point: density shape mismatch'
      if (size(ao_atom) /= n .or. size(aog1,1) /= n .or. &
          size(aog2,1) /= n .or. size(aog3,1) /= n) &
        error stop 'gga_density_nuclear_point: AO dimension mismatch'
      if (size(aog1,2) /= 3 .or. size(aog2,2) /= 6 .or. size(aog3,2) /= 10) &
        error stop 'gga_density_nuclear_point: derivative component mismatch'
      if (size(drho,1) /= 3 .or. size(dgrho,1) /= 3 .or. &
          size(dgrho,2) /= 3 .or. size(dgrho,3) /= n_atom) &
        error stop 'gga_density_nuclear_point: first derivative output mismatch'
      if (any(shape(d2rho) /= [3,3,n_atom,n_atom]) .or. &
          any(shape(d2grho) /= [3,3,3,n_atom,n_atom])) &
        error stop 'gga_density_nuclear_point: second derivative output mismatch'
    end subroutine check_shapes

  end subroutine gga_density_nuclear_point

!> Evaluate the same exact nuclear derivatives for several density matrices.
!> The leading matrix index is contiguous so all matrices reuse each AO-pair
!> derivative before the next pair is visited.
  subroutine gga_density_nuclear_point_batch(density,ao_atom,aov,aog1,aog2, &
      aog3,drho,dgrho,d2rho,d2grho)
    real(fp), intent(in) :: density(:,:,:)
    integer, intent(in) :: ao_atom(:)
    real(fp), intent(in) :: aov(:),aog1(:,:),aog2(:,:),aog3(:,:)
    real(fp), intent(out) :: drho(:,:,:),dgrho(:,:,:,:)
    real(fp), intent(out) :: d2rho(:,:,:,:,:),d2grho(:,:,:,:,:,:)

    integer :: mu,nu,a,b,c,atom_a,atom_b,nao,nat,nmat
    real(fp) :: pair_density(size(density,1))
    real(fp) :: contracted_value(size(density,1),size(aov))
    real(fp) :: contracted_gradient(size(density,1),3,size(aov))

    nao=size(aov)
    nmat=size(density,1)
    nat=size(drho,2)
    if(any(shape(density)/=[nmat,nao,nao]) .or. size(ao_atom)/=nao .or. &
       any(shape(aog1)/=[nao,3]) .or. any(shape(aog2)/=[nao,6]) .or. &
       any(shape(aog3)/=[nao,10]) .or. &
       any(shape(drho)/=[3,nat,nmat]) .or. &
       any(shape(dgrho)/=[3,3,nat,nmat]) .or. &
       any(shape(d2rho)/=[3,3,nat,nat,nmat]) .or. &
       any(shape(d2grho)/=[3,3,3,nat,nat,nmat])) &
      error stop 'gga_density_nuclear_point_batch: shape mismatch'
    drho=0.0_fp
    dgrho=0.0_fp
    d2rho=0.0_fp
    d2grho=0.0_fp
    contracted_value=0.0_fp
    contracted_gradient=0.0_fp

    ! Contract the symmetrized density with AO values and gradients once.
    ! The two unconditional target updates preserve the two product-rule
    ! contributions for diagonal AO pairs.
    do nu=1,nao
      do mu=1,nu
        pair_density=density(:,mu,nu)
        if(mu/=nu) pair_density=pair_density+density(:,nu,mu)
        if(all(pair_density==0.0_fp)) cycle
        contracted_value(:,mu)=contracted_value(:,mu)+pair_density*aov(nu)
        contracted_value(:,nu)=contracted_value(:,nu)+pair_density*aov(mu)
        do c=1,3
          contracted_gradient(:,c,mu)=contracted_gradient(:,c,mu)+ &
            pair_density*aog1(nu,c)
          contracted_gradient(:,c,nu)=contracted_gradient(:,c,nu)+ &
            pair_density*aog1(mu,c)
        end do
      end do
    end do

    ! Terms where one or both nuclear derivatives act on the same AO.
    do mu=1,nao
      atom_a=ao_atom(mu)
      do a=1,3
        drho(a,atom_a,:)=drho(a,atom_a,:)- &
          aog1(mu,a)*contracted_value(:,mu)
        do c=1,3
          dgrho(c,a,atom_a,:)=dgrho(c,a,atom_a,:)- &
            aog2(mu,hmap(a,c))*contracted_value(:,mu)- &
            aog1(mu,a)*contracted_gradient(:,c,mu)
        end do
        do b=1,3
          d2rho(a,b,atom_a,atom_a,:)=d2rho(a,b,atom_a,atom_a,:)+ &
            aog2(mu,hmap(a,b))*contracted_value(:,mu)
          do c=1,3
            d2grho(c,a,b,atom_a,atom_a,:)= &
              d2grho(c,a,b,atom_a,atom_a,:)+ &
              aog3(mu,tmap(a,b,c))*contracted_value(:,mu)+ &
              aog2(mu,hmap(a,b))*contracted_gradient(:,c,mu)
          end do
        end do
      end do
    end do

    ! Cross terms where the two derivatives act on the two AO factors.
    do nu=1,nao
      atom_b=ao_atom(nu)
      do mu=1,nu
        atom_a=ao_atom(mu)
        pair_density=density(:,mu,nu)
        if(mu/=nu) pair_density=pair_density+density(:,nu,mu)
        if(all(pair_density==0.0_fp)) cycle
        do a=1,3
          do b=1,3
            d2rho(a,b,atom_a,atom_b,:)=d2rho(a,b,atom_a,atom_b,:)+ &
              pair_density*aog1(mu,a)*aog1(nu,b)
            d2rho(a,b,atom_b,atom_a,:)=d2rho(a,b,atom_b,atom_a,:)+ &
              pair_density*aog1(nu,a)*aog1(mu,b)
            do c=1,3
              d2grho(c,a,b,atom_a,atom_b,:)= &
                d2grho(c,a,b,atom_a,atom_b,:)+pair_density*( &
                aog2(mu,hmap(a,c))*aog1(nu,b)+ &
                aog1(mu,a)*aog2(nu,hmap(b,c)))
              d2grho(c,a,b,atom_b,atom_a,:)= &
                d2grho(c,a,b,atom_b,atom_a,:)+pair_density*( &
                aog2(nu,hmap(a,c))*aog1(mu,b)+ &
                aog1(nu,a)*aog2(mu,hmap(b,c)))
            end do
          end do
        end do
      end do
    end do

  contains

    pure real(fp) function center_d1_batch(i,atom,ixyz)
      integer, intent(in) :: i,atom,ixyz
      center_d1_batch=0.0_fp
      if(ao_atom(i)==atom) center_d1_batch=-aog1(i,ixyz)
    end function center_d1_batch

    pure real(fp) function center_gd1_batch(i,atom,ixyz,cxyz)
      integer, intent(in) :: i,atom,ixyz,cxyz
      center_gd1_batch=0.0_fp
      if(ao_atom(i)==atom) center_gd1_batch=-aog2(i,hmap(ixyz,cxyz))
    end function center_gd1_batch

    pure real(fp) function center_d2_batch(i,atom1,xyz1,atom2,xyz2)
      integer, intent(in) :: i,atom1,xyz1,atom2,xyz2
      center_d2_batch=0.0_fp
      if(ao_atom(i)==atom1 .and. atom1==atom2) &
        center_d2_batch=aog2(i,hmap(xyz1,xyz2))
    end function center_d2_batch

    pure real(fp) function center_gd2_batch(i,atom1,xyz1,atom2,xyz2,cxyz)
      integer, intent(in) :: i,atom1,xyz1,atom2,xyz2,cxyz
      center_gd2_batch=0.0_fp
      if(ao_atom(i)==atom1 .and. atom1==atom2) &
        center_gd2_batch=aog3(i,tmap(xyz1,xyz2,cxyz))
    end function center_gd2_batch

  end subroutine gga_density_nuclear_point_batch

!> Compute fixed-grid first nuclear derivatives for a batch of spin pairs.
!>
!> The leading density index labels independent right-hand sides.  Only the
!> one or two AO centers of each AO pair can contribute, which is an exact
!> property of atom-centred Gaussian basis functions.
  pure subroutine gga_density_nuclear_point_first_batch(density,ao_atom, &
      aov,aog1,aog2,drho,dgrho)
    real(fp), intent(in) :: density(:,:,:,:),aov(:),aog1(:,:),aog2(:,:)
    integer, intent(in) :: ao_atom(:)
    real(fp), intent(out) :: drho(:,:,:,:),dgrho(:,:,:,:,:)
    real(fp) :: pair_density(size(density,1),2)
    real(fp) :: contracted_value(size(density,1),2,size(aov))
    real(fp) :: contracted_gradient(size(density,1),3,2,size(aov))
    integer :: atom,cart,mu,nu,space,spin,nao,nat,nvec

    nao=size(aov)
    nvec=size(density,1)
    nat=size(drho,2)
    if(any(shape(density)/=[nvec,2,nao,nao]) .or. &
       size(ao_atom)/=nao .or. any(shape(aog1)/=[nao,3]) .or. &
       any(shape(aog2)/=[nao,6]) .or. &
       any(shape(drho)/=[3,nat,2,nvec]) .or. &
       any(shape(dgrho)/=[3,3,nat,2,nvec])) &
      error stop 'gga_density_nuclear_point_first_batch: shape mismatch'
    drho=0.0_fp
    dgrho=0.0_fp
    contracted_value=0.0_fp
    contracted_gradient=0.0_fp

    ! First contract every symmetrized density with the AO value and gradient.
    ! Both target AOs are updated even on the diagonal; this reproduces the
    ! two product-rule terms for phi_mu*phi_mu exactly.
    do nu=1,nao
      do mu=1,nu
        pair_density=density(:,:,mu,nu)
        if(mu/=nu) pair_density=pair_density+density(:,:,nu,mu)
        do spin=1,2
          contracted_value(:,spin,mu)=contracted_value(:,spin,mu)+ &
            pair_density(:,spin)*aov(nu)
          contracted_value(:,spin,nu)=contracted_value(:,spin,nu)+ &
            pair_density(:,spin)*aov(mu)
          do space=1,3
            contracted_gradient(:,space,spin,mu)= &
              contracted_gradient(:,space,spin,mu)+ &
              pair_density(:,spin)*aog1(nu,space)
            contracted_gradient(:,space,spin,nu)= &
              contracted_gradient(:,space,spin,nu)+ &
              pair_density(:,spin)*aog1(mu,space)
          end do
        end do
      end do
    end do
    do mu=1,nao
      atom=ao_atom(mu)
      do cart=1,3
        do spin=1,2
          drho(cart,atom,spin,:)=drho(cart,atom,spin,:)- &
            aog1(mu,cart)*contracted_value(:,spin,mu)
          do space=1,3
            dgrho(space,cart,atom,spin,:)= &
              dgrho(space,cart,atom,spin,:)- &
              aog2(mu,hmap(cart,space))*contracted_value(:,spin,mu)- &
              aog1(mu,cart)*contracted_gradient(:,space,spin,mu)
          end do
        end do
      end do
    end do
  end subroutine gga_density_nuclear_point_first_batch

end module mod_dft_gga_nuclear_point
