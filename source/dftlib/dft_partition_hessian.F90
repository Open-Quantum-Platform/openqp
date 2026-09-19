module mod_dft_partition_hessian

  use precision, only: fp
  use mod_dft_partfunc, only: partition_function, PTYPE_SSF, PTYPE_ERF, &
    PTYPE_BECKE4, PTYPE_BECKE3, PTYPE_SMSTP2, PTYPE_SMSTP3, &
    PTYPE_SMSTP4, PTYPE_SMSTP5
  implicit none
  private

  public :: partition_weight_nuclear_derivatives
  public :: partition_function_second_derivative

contains

!> Return all normalized atom-centred partition weights and their first and
!> second nuclear derivatives for a grid point rigidly attached to owner.
!>
!> The output ordering is (Cartesian component, derivative atom, cell atom)
!> for dweight and (component, atom, component, atom, cell atom) for d2weight.
!> surface_shift(j,i) follows dft_fuzzycell/dft_gridint_grad for the pair i>j.
  subroutine partition_weight_nuclear_derivatives(atom_xyz, point, owner, &
      dummy_atom, part_fun_type, has_surface_shift, surface_shift, weight, &
      dweight, d2weight, status)

    real(fp), intent(in) :: atom_xyz(:,:), point(3), surface_shift(:,:)
    integer, intent(in) :: owner, part_fun_type
    logical, intent(in) :: dummy_atom(:), has_surface_shift
    real(fp), intent(out) :: weight(:), dweight(:,:,:), d2weight(:,:,:,:,:)
    integer, intent(out) :: status

    type(partition_function) :: partfunc
    real(fp), allocatable :: cell(:), gcell(:,:), hcell(:,:,:)
    real(fp), allocatable :: rho(:), grho(:,:), hrho(:,:,:)
    real(fp), allocatable :: gpair(:), hpair(:,:), gmu(:), hmu(:,:)
    real(fp), allocatable :: gf(:), hf(:,:), gsum(:), hsum(:,:)
    real(fp) :: rvec(3), unit(3), rijvec(3), rhat(3)
    real(fp) :: rij, mu0, mu, aij, fac, f, fp1, fp2, sumc
    integer :: nat, nq, i, j, b, c, alpha, beta, ib, ic

    status = 0
    nat = size(atom_xyz,2)
    nq = 3*nat
    if (size(atom_xyz,1) /= 3 .or. owner < 1 .or. owner > nat .or. &
        size(dummy_atom) /= nat .or. any(shape(surface_shift) /= [nat,nat]) .or. &
        size(weight) /= nat .or. any(shape(dweight) /= [3,nat,nat]) .or. &
        any(shape(d2weight) /= [3,nat,3,nat,nat])) then
      status = 1
      return
    end if

    allocate(cell(nat), gcell(nq,nat), hcell(nq,nq,nat), rho(nat), &
      grho(nq,nat), hrho(nq,nq,nat), gpair(nq), hpair(nq,nq), &
      gmu(nq), hmu(nq,nq), gf(nq), hf(nq,nq), gsum(nq), hsum(nq,nq))
    cell = 1.0_fp
    where (dummy_atom) cell = 0.0_fp
    gcell = 0.0_fp
    hcell = 0.0_fp
    grho = 0.0_fp
    hrho = 0.0_fp

    ! rho_i=|r_grid-R_i|, with r_grid translated rigidly with R_owner.
    do i = 1, nat
      rvec = point-atom_xyz(:,i)
      rho(i) = norm2(rvec)
      if (rho(i) <= sqrt(tiny(1.0_fp))) then
        status = 2
        return
      end if
      unit = rvec/rho(i)
      do b = 1, nat
        fac = real(merge(1,0,b == owner)-merge(1,0,b == i),fp)
        do alpha = 1, 3
          ib = 3*(b-1)+alpha
          grho(ib,i) = fac*unit(alpha)
          do c = 1, nat
            do beta = 1, 3
              ic = 3*(c-1)+beta
              hrho(ib,ic,i) = fac* &
                real(merge(1,0,c == owner)-merge(1,0,c == i),fp)* &
                (real(merge(1,0,alpha == beta),fp)-unit(alpha)*unit(beta))/rho(i)
            end do
          end do
        end do
      end do
    end do

    call partfunc%set(part_fun_type)
    do i = 2, nat
      if (dummy_atom(i)) cycle
      do j = 1, i-1
        if (dummy_atom(j)) cycle
        rijvec = atom_xyz(:,i)-atom_xyz(:,j)
        rij = norm2(rijvec)
        if (rij <= sqrt(tiny(1.0_fp))) then
          status = 3
          return
        end if
        rhat = rijvec/rij
        gpair = 0.0_fp
        hpair = 0.0_fp
        do b = 1, nat
          fac = real(merge(1,0,b == i)-merge(1,0,b == j),fp)
          do alpha = 1, 3
            ib = 3*(b-1)+alpha
            gpair(ib) = fac*rhat(alpha)
            do c = 1, nat
              do beta = 1, 3
                ic = 3*(c-1)+beta
                hpair(ib,ic) = fac* &
                  real(merge(1,0,c == i)-merge(1,0,c == j),fp)* &
                  (real(merge(1,0,alpha == beta),fp)-rhat(alpha)*rhat(beta))/rij
              end do
            end do
          end do
        end do

        mu0 = (rho(i)-rho(j))/rij
        do ib = 1, nq
          gmu(ib) = (grho(ib,i)-grho(ib,j)-mu0*gpair(ib))/rij
          do ic = 1, nq
            hmu(ib,ic) = (hrho(ib,ic,i)-hrho(ib,ic,j))/rij &
              -(grho(ib,i)-grho(ib,j))*gpair(ic)/(rij*rij) &
              -gpair(ib)*(grho(ic,i)-grho(ic,j))/(rij*rij) &
              -mu0*hpair(ib,ic)/rij &
              +2.0_fp*mu0*gpair(ib)*gpair(ic)/(rij*rij)
          end do
        end do

        mu = mu0
        if (has_surface_shift) then
          aij = surface_shift(j,i)
          mu = mu0+aij*(1.0_fp-mu0*mu0)
          fac = 1.0_fp-2.0_fp*aij*mu0
          hmu = fac*hmu-2.0_fp*aij*spread(gmu,2,nq)*spread(gmu,1,nq)
          gmu = fac*gmu
        end if

        f = partfunc%eval(mu)
        fp1 = partfunc%deriv(mu)
        fp2 = partition_function_second_derivative(part_fun_type,mu,status)
        if (status /= 0) return
        gf = fp1*gmu
        hf = fp1*hmu+fp2*spread(gmu,2,nq)*spread(gmu,1,nq)
        call multiply_second_order(cell(i),gcell(:,i),hcell(:,:,i),f,gf,hf)
        call multiply_second_order(cell(j),gcell(:,j),hcell(:,:,j), &
                                   1.0_fp-f,-gf,-hf)
      end do
    end do

    sumc = sum(cell)
    if (sumc <= tiny(1.0_fp)) then
      status = 4
      return
    end if
    gsum = sum(gcell,dim=2)
    hsum = sum(hcell,dim=3)
    do i = 1, nat
      weight(i) = cell(i)/sumc
      do b = 1, nat
        do alpha = 1, 3
          ib = 3*(b-1)+alpha
          dweight(alpha,b,i) = gcell(ib,i)/sumc-cell(i)*gsum(ib)/(sumc*sumc)
          do c = 1, nat
            do beta = 1, 3
              ic = 3*(c-1)+beta
              d2weight(alpha,b,beta,c,i) = hcell(ib,ic,i)/sumc &
                -(gcell(ib,i)*gsum(ic)+gsum(ib)*gcell(ic,i))/(sumc*sumc) &
                -cell(i)*hsum(ib,ic)/(sumc*sumc) &
                +2.0_fp*cell(i)*gsum(ib)*gsum(ic)/(sumc*sumc*sumc)
            end do
          end do
        end do
      end do
    end do
  end subroutine partition_weight_nuclear_derivatives

!> Multiply two scalar second-order jets without division, including at a
!> saturated partition factor of exactly zero.
  pure subroutine multiply_second_order(value, gradient, hessian, factor, &
                                        factor_gradient, factor_hessian)
    real(fp), intent(inout) :: value, gradient(:), hessian(:,:)
    real(fp), intent(in) :: factor, factor_gradient(:), factor_hessian(:,:)
    real(fp) :: old_value
    real(fp) :: old_gradient(size(gradient)), old_hessian(size(hessian,1),size(hessian,2))
    old_value = value
    old_gradient = gradient
    old_hessian = hessian
    hessian = factor*old_hessian+old_value*factor_hessian &
      +spread(old_gradient,2,size(gradient))*spread(factor_gradient,1,size(gradient)) &
      +spread(factor_gradient,2,size(gradient))*spread(old_gradient,1,size(gradient))
    gradient = factor*old_gradient+old_value*factor_gradient
    value = old_value*factor
  end subroutine multiply_second_order

!> Analytic second derivative of every partition function in dft_partfunc.
  function partition_function_second_derivative(ptype,x,status) result(d2f)
    integer, intent(in) :: ptype
    real(fp), intent(in) :: x
    integer, intent(out) :: status
    real(fp) :: d2f, y, dy, d2y, gp, gpp, limit, scale, q, qp, qpp
    integer :: iter, order

    status = 0
    select case (ptype)
    case (PTYPE_BECKE4,PTYPE_BECKE3)
      y = x
      dy = 1.0_fp
      d2y = 0.0_fp
      iter = merge(4,3,ptype == PTYPE_BECKE4)
      do order = 1, iter
        gp = 1.5_fp*(1.0_fp-y*y)
        gpp = -3.0_fp*y
        d2y = gpp*dy*dy+gp*d2y
        dy = gp*dy
        y = 0.5_fp*y*(3.0_fp-y*y)
      end do
      d2f = -0.5_fp*d2y
    case (PTYPE_SSF)
      limit = 0.64_fp
      if (abs(x) > limit) then
        d2f = 0.0_fp
      else
        y = x/limit
        d2f = 105.0_fp*y*(1.0_fp-y*y)**2/(16.0_fp*limit*limit)
      end if
    case (PTYPE_ERF)
      limit = 0.725_fp
      scale = 1.0_fp/0.3_fp
      if (abs(x) > limit) then
        d2f = 0.0_fp
      else
        q = x/(1.0_fp-x*x)
        qp = (1.0_fp+x*x)/(1.0_fp-x*x)**2
        qpp = 2.0_fp*x*(3.0_fp+x*x)/(1.0_fp-x*x)**3
        d2f = -scale/sqrt(acos(-1.0_fp))*exp(-(scale*q)**2)* &
          (qpp-2.0_fp*scale*scale*q*qp*qp)
      end if
    case (PTYPE_SMSTP2,PTYPE_SMSTP3,PTYPE_SMSTP4,PTYPE_SMSTP5)
      select case (ptype)
      case (PTYPE_SMSTP2)
        limit=0.55_fp; order=2; scale=30.0_fp
      case (PTYPE_SMSTP3)
        limit=0.62_fp; order=3; scale=140.0_fp
      case (PTYPE_SMSTP4)
        limit=0.69_fp; order=4; scale=630.0_fp
      case default
        limit=0.73_fp; order=5; scale=2772.0_fp
      end select
      if (abs(x) > limit) then
        d2f = 0.0_fp
      else
        y = 0.5_fp-0.5_fp*x/limit
        d2f = (0.5_fp/limit)**2*scale*real(order,fp)* &
          y**(order-1)*(1.0_fp-y)**(order-1)*(1.0_fp-2.0_fp*y)
      end if
    case default
      status = 5
      d2f = 0.0_fp
    end select
  end function partition_function_second_derivative

end module mod_dft_partition_hessian
