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

  public :: gga_density_nuclear_point
  public :: gga_density_nuclear_point_reference

contains

!> Same results as gga_density_nuclear_point_reference, evaluated without the
!> atom-pair loop.  A center derivative of AO mu is nonzero only for its own
!> atom, so every term reduces to per-AO vectors and atom-block sums:
!>   u   = P phi + P^T phi,   w^c = P g^c + P^T g^c,
!>   S(X,Y)_AB = sum_{mu in A, nu in B} X_mu P_mu,nu Y_nu   (X,Y in {g, h}),
!>   drho(a,A)        = -sum_{mu in A} g^a u
!>   dgrho(c,a,A)     = -sum_{mu in A} (h^ac u + g^a w^c)
!>   d2rho(a,b,A,B)   = d_AB sum_{mu in A} h^ab u + S(g^a,g^b)_AB + S(g^b,g^a)_BA
!>   d2grho(c,a,b,A,B)= d_AB sum_{mu in A} (t^abc u + h^ab w^c)
!>                      + S(h^ac,g^b)_AB + S(h^bc,g^a)_BA
!>                      + S(g^a,h^bc)_AB + S(g^b,h^ac)_BA
!> (g = first, h = second, t = third electronic AO derivatives).  Cost per
!> point: O(nao^2) contractions plus O(nao*nat) block sums, instead of
!> O(nao^2 nat^2).
  subroutine gga_density_nuclear_point(density, ao_atom, aov, aog1, aog2, &
                                       aog3, drho, dgrho, d2rho, d2grho)
    real(fp), intent(in) :: density(:,:)
    integer, intent(in) :: ao_atom(:)
    real(fp), intent(in) :: aov(:), aog1(:,:), aog2(:,:), aog3(:,:)
    real(fp), intent(out) :: drho(:,:), dgrho(:,:,:)
    real(fp), intent(out) :: d2rho(:,:,:,:), d2grho(:,:,:,:,:)

    integer :: mu, nu, a, b, c, k, kx, ky, ia, ib, nao, nat
    real(fp) :: p
    real(fp), allocatable :: y(:,:), u(:), w(:,:), r(:,:,:), sblk(:,:,:,:)

    nao = size(aov)
    nat = size(drho,2)
    if (size(density,1) /= nao .or. size(density,2) /= nao .or. size(ao_atom) /= nao) &
      error stop 'gga_density_nuclear_point: AO dimension mismatch'
    if (any(shape(d2rho) /= [3,3,nat,nat]) .or. &
        any(shape(d2grho) /= [3,3,3,nat,nat])) &
      error stop 'gga_density_nuclear_point: second derivative output mismatch'

    ! y(:,1:3) = g^x,y,z ; y(:,3+k) = h column k (k = hmap index)
    allocate(y(nao,9), u(nao), w(nao,3), r(nao,9,nat), sblk(9,9,nat,nat))
    y(:,1:3) = aog1(:,1:3)
    y(:,4:9) = aog2(:,1:6)

    u = 0.0_fp; w = 0.0_fp; r = 0.0_fp
    do nu = 1, nao
      ib = ao_atom(nu)
      do mu = 1, nao
        p = density(mu,nu)
        if (p == 0.0_fp) cycle
        u(mu) = u(mu) + p*aov(nu)
        u(nu) = u(nu) + p*aov(mu)
        do c = 1, 3
          w(mu,c) = w(mu,c) + p*aog1(nu,c)
          w(nu,c) = w(nu,c) + p*aog1(mu,c)
        end do
        do k = 1, 9
          r(mu,k,ib) = r(mu,k,ib) + p*y(nu,k)
        end do
      end do
    end do
    sblk = 0.0_fp
    do mu = 1, nao
      ia = ao_atom(mu)
      do ib = 1, nat
        do ky = 1, 9
          if (r(mu,ky,ib) == 0.0_fp) cycle
          do kx = 1, 9
            sblk(kx,ky,ia,ib) = sblk(kx,ky,ia,ib) + y(mu,kx)*r(mu,ky,ib)
          end do
        end do
      end do
    end do

    drho = 0.0_fp; dgrho = 0.0_fp; d2rho = 0.0_fp; d2grho = 0.0_fp
    do mu = 1, nao
      ia = ao_atom(mu)
      do a = 1, 3
        drho(a,ia) = drho(a,ia) - aog1(mu,a)*u(mu)
        do c = 1, 3
          dgrho(c,a,ia) = dgrho(c,a,ia) - aog2(mu,hmap(a,c))*u(mu) - aog1(mu,a)*w(mu,c)
        end do
        do b = 1, 3
          d2rho(a,b,ia,ia) = d2rho(a,b,ia,ia) + aog2(mu,hmap(a,b))*u(mu)
          do c = 1, 3
            d2grho(c,a,b,ia,ia) = d2grho(c,a,b,ia,ia) &
              + aog3(mu,tmap(a,b,c))*u(mu) + aog2(mu,hmap(a,b))*w(mu,c)
          end do
        end do
      end do
    end do
    do ib = 1, nat
      do ia = 1, nat
        do b = 1, 3
          do a = 1, 3
            d2rho(a,b,ia,ib) = d2rho(a,b,ia,ib) + sblk(a,b,ia,ib) + sblk(b,a,ib,ia)
            do c = 1, 3
              d2grho(c,a,b,ia,ib) = d2grho(c,a,b,ia,ib) &
                + sblk(3+hmap(a,c),b,ia,ib) + sblk(3+hmap(b,c),a,ib,ia) &
                + sblk(a,3+hmap(b,c),ia,ib) + sblk(b,3+hmap(a,c),ib,ia)
            end do
          end do
        end do
      end do
    end do
    deallocate(y, u, w, r, sblk)
  end subroutine gga_density_nuclear_point

!> Compute fixed-grid first and second nuclear derivatives of rho and grad(rho).
!>
!> AO derivatives are electronic-coordinate derivatives ordered as
!> G2=(xx,yy,zz,xy,yz,xz) and
!> G3=(xxx,yyy,zzz,xxy,xxz,yyx,yyz,zzx,zzy,xyz).
!> Moving the center of AO mu gives d(phi_mu)/dR_Aa=-delta(mu,A)d_a(phi_mu).
!> Thus this routine supplies the DRA/GDA/DGGA data underlying the GAMESS
!> DDDENCNST and TDHXG1G/TDHXGPG construction, without grid-weight or
!> grid-center translation terms.
  subroutine gga_density_nuclear_point_reference(density, ao_atom, aov, aog1, aog2, &
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

  end subroutine gga_density_nuclear_point_reference

end module mod_dft_gga_nuclear_point
