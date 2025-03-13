module int2e_rotaxis
  use precision, only: dp, qp
  use basis_tools, only: basis_set
  use boys_lut, only: fgrid, xgrid, rxinc, rfinc, rmr, tmax
  use int2_pairs, only: int2_pair_storage, int2_cutoffs_t
  use constants, only: pi
  implicit none

  private
  public genr22

  real(dp), parameter :: acy_threshold = 1.0e-10_dp
  real(dp), parameter :: sqrt3 = sqrt(3.0_dp)
  real(dp), parameter :: pi4 = pi/4

  integer, parameter :: in6(6) = [ 1, 4, 5, 2, 6, 3 ]
  integer, parameter :: max_contraction = 30 !< Max. degree of contraction

  type :: rotaxis_data_t
    logical :: lrint = .false.
    real(kind=dp) :: emu2 = 1.0d99

    integer :: ngangb
    real(kind=dp) :: cutoff = 0

    real(kind=dp) :: acy, acy2, aqx, aqx2, aqxy, y03, y04
    real(kind=dp) :: rab, x34, x43, aqz, qps, sq

    real(kind=dp) :: tx12(max_contraction**2), ty02(max_contraction**2), &
          sp(max_contraction**2)


    real(kind=dp) :: fq(0:8), fq0(5), fq1(2,13), fq2(3,16), fq3(4,16), fq4(5,16), &
                fq5(6,11), fq6(7,7), fq7(8,3), fq8( 9)

    real(kind=dp) :: r00(5,5), r01(3,40), r02(6,56), r03(10,52), r04(15,42), &
       r05(21,24), r06(28,12), r07(36,4), r08(45)

  end type

contains

! >
! >    @brief   rotated axis integration involving s,p,l,d shells
! >
! >    @details rotated axis integration involving s,p,l,d shells,
! >             by a mix of rotated axis and mcmurchie/davidson.
! >                       k.ishimura, s.nagase
! >                theoret.chem.acc. 120, 185-189(2008)
! >
! >    @author  kazuya ishimura, at the institute for molecular science,
! >             sponsored by naregi nano science project, in 2004.
! >             extensively revised by jose sierra in 2013.
! >
  subroutine genr22(basis, ppairs, grotspd, shell_ids, flips, cutoffs, emu2)

    implicit none

    type(basis_set), intent(in) :: basis
    type(int2_pair_storage), intent(in) :: ppairs
    real(kind=dp), intent(inout) :: grotspd(*)
    integer, intent(in) :: shell_ids(4)
    integer, intent(out) :: flips(4)
    type(int2_cutoffs_t), intent(in) :: cutoffs
    real(kind=dp), optional :: emu2

    real(kind=dp) :: a(3), b(3), c(3), d(3), p(3,3), t(3)

    integer :: flips_inv(4)
    integer :: angm(4)
    integer :: shl_new(4)
    integer :: itype, jtype
    integer :: i
    integer :: ji

    real(kind=dp) :: acx, acz, cq, cqx, cqz, qx, qz
    real(kind=dp) :: cosg, sing
    real(kind=dp) :: rcd
    real(kind=dp) :: tmp
    real(kind=dp) :: x03, x04
    integer :: nga, ngb, ngc, ngd

    integer :: ppid_p, ppid_q, npp_p, npp_q
    integer :: id1, id2

    integer, parameter :: jtypes(*) = [ &
       1, 2, 7, 2, 3, &
       8, 7, 8, 10, 2, &
       4, 9, 4, 5, 11, &
       9, 11, 14, 7, 9, &
       12, 9, 13, 15, 12, &
       15, 17, 2, 4, 9, &
       4, 5, 11, 9, 11, &
       14, 3, 5, 13, 5, &
       6, 16, 13, 16, 18, &
       8, 11, 15, 11, 16, &
       19, 15, 19, 20, 7, &
       9, 12, 9, 13, 15, &
       12, 15, 17, 8, 11, &
       15, 11, 16, 19, 15, &
       19, 20, 10, 14, 17, &
       14, 18, 20, 17, 20, &
       21]

    type(rotaxis_data_t) :: rdat

    angm = basis%am(shell_ids)

    itype = 1 + angm(4) + 3*angm(3) + 9*angm(2) + 27*angm(1)
    jtype = jtypes(itype)

    flips = [1,2,3,4]
    if (angm(1) > angm(2)) then
      flips(1:2) = [2,1]
      angm(1:2) = angm([2,1])
    end if
    if (angm(3) > angm(4)) then
      flips(3:4) = [4,3]
      angm(3:4) = angm([4,3])
    end if
    if (sum(angm(1:2)) > sum(angm(3:4)) .or. angm(2) > angm(4)) then
      flips = flips([3,4,1,2])
      angm = angm([3,4,1,2])
    end if

    flips_inv(flips) = [1, 2, 3, 4]
    shl_new = shell_ids(flips_inv)
! If (la/=angm(1).or.lb/=angm(2)) print *, la, angm(1), lb, angm(2)
! Empty integral summation storage

    call intclean(rdat,jtype)

    if (present(emu2)) then
      rdat%lrint = .true.
      rdat%emu2 = emu2
    end if

    id1 = maxval(shl_new([1,2]))
    id2 = minval(shl_new([1,2]))
    npp_p = ppairs%ppid(1,id1*(id1-1)/2+id2)
    ppid_p = ppairs%ppid(2,id1*(id1-1)/2+id2)

    id1 = maxval(shl_new([3,4]))
    id2 = minval(shl_new([3,4]))
    npp_q = ppairs%ppid(1,id1*(id1-1)/2+id2)
    ppid_q = ppairs%ppid(2,id1*(id1-1)/2+id2)

! If (npp_p == 0 .or. npp_q == 0) return

! Obtain information about shells: inew, knew, jnew, lnew
! Number of gaussians go into nga,... in common shllfo
! Shell angular quantum numbers la,... go into common shllfo
! Gaussian exponents go into arrays exa,exb,exc,exd in common shllfo
! Gaussian coefficients go into arrays csa,cpa,... in common shllfo

! Loop over gaussians in each shell
! First shell inew
! Coordinates of atoms associated with shells inew jnew knew and lnew

    a = ppairs%p(:,ppid_p)-ppairs%pa(:, ppid_p)
    b = ppairs%p(:,ppid_p)-ppairs%pb(:, ppid_p)

    c = ppairs%p(:,ppid_q)-ppairs%pa(:, ppid_q)
    d = ppairs%p(:,ppid_q)-ppairs%pb(:, ppid_q)

! Find direction cosines of penultimate axes from coordinates of ab
! P(1,1),p(1,2),... are direction cosines of axes at p.  z-axis along ab
! T(1),t(2),t(3)... are direction cosines of axes at q.  z-axis along cd

! Find direction cosines of ab and cd. these are local z-axes.
! If indeterminate take along space z-axis

    p(:,3) = [0, 0, 1]
    rdat%rab = ppairs%rab(ppid_p)
    if (rdat%rab>0) then
      p(:,3) = (b-a)*ppairs%uab(ppid_p)
    end if

    t = [0, 0, 1]
    rcd = ppairs%rab(ppid_q)
    if (rcd>0) then
      t = (d-c)*ppairs%uab(ppid_q)
    end if

! Find local y-axis as common perpendicular to ab and cd
! If indeterminate take perpendicular to ab and space z-axis
! If still indeterminate take perpendicular to ab and space x-axis

    cosg = dot_product(t,p(:,3))

! Modified rotation testing.
! This fix cures the small angle problem.

    p(1,2) = t(3)*p(2,3) - t(2)*p(3,3)
    p(2,2) = t(1)*p(3,3) - t(3)*p(1,3)
    p(3,2) = t(2)*p(1,3) - t(1)*p(2,3)
    if (abs(cosg)>0.9) then
       sing = norm2(p(:,2))
    else
       sing = sqrt(1-cosg*cosg)
    end if

    if (sing<1d-12) then
       if (abs(p(1,3))<sqrt(0.5d0)) then
          tmp = 1/sqrt(1-p(1,3)*p(1,3))
          p(:,2) = [0d0, p(3,3), -p(2,3)] * tmp
       else
          tmp = 1/sqrt(1-p(3,3)*p(3,3))
          p(:,2) = [p(2,3), -p(1,3), 0d0] * tmp
       end if
    else
       p(:,2)= p(:,2)/sing
    end if

    p(1,1) = p(2,2)*p(3,3)-p(3,2)*p(2,3)
    p(2,1) = p(3,2)*p(1,3)-p(1,2)*p(3,3)
    p(3,1) = p(1,2)*p(2,3)-p(2,2)*p(1,3)

! Find coordinates of c relative to local axes at a

    t = c - a
    acx = dot_product(t,p(:,1))
    rdat%acy = dot_product(t,p(:,2))
    acz = dot_product(t,p(:,3))

! Set acy= 0  if close

    if (abs(rdat%acy)<=acy_threshold) then
       rdat%acy = 0.0_dp
       rdat%acy2 = 0.0_dp
    else
       rdat%acy2 = rdat%acy*rdat%acy
    end if

! Direction cosines of cd local axes with respect to ab local axes
! ( cosg,   0,-sing )
! (    0,   1,    0 )
! ( sing,   0, cosg )


! Preliminary p loop

! Fill geompq with information about p in preliminary p-loop

    rdat%cutoff = cutoffs%quartet_cutoff_squared
    ji = 0
    do i = 0, npp_p-1
      if (abs(ppairs%k(ppid_p+i)*ppairs%ginv(ppid_p+i)) < cutoffs%quartet_cutoff) cycle
      ji = ji+1
      rdat%tx12(ji) = ppairs%g(ppid_p+i)
      rdat%ty02(ji) = ppairs%alpha_b(ppid_p+i)*ppairs%ginv(ppid_p+i)*rdat%rab
      rdat%sp(ji) = ppairs%k(ppid_p+i)*ppairs%ginv(ppid_p+i)
    end do
    rdat%ngangb = ji

! Begin q loop

    do i = 0, npp_q-1
      if (abs(ppairs%k(ppid_q+i)*ppairs%ginv(ppid_q+i)) < cutoffs%quartet_cutoff) cycle
      x03 = ppairs%alpha_a(ppid_q+i)
      x04 = ppairs%alpha_b(ppid_q+i)

      rdat%x34 = ppairs%g(ppid_q+i)
      rdat%x43 = ppairs%ginv(ppid_q+i)
      rdat%y03 = x03*rdat%x43
      rdat%y04 = x04*rdat%x43
      rdat%sq = ppairs%k(ppid_q+i)*ppairs%ginv(ppid_q+i)

! Cqx = component of cq along penultimate x-axis
! Cqz = component of cq along penultimate z-axis
      cq = rcd*rdat%y04
      cqx = cq*sing
      cqz = cq*cosg

! Find coordinates of q relative to axes at a
! Qpr is perpendicular from q to ab

      rdat%aqx = acx+cqx
      rdat%aqx2 = rdat%aqx*rdat%aqx
      rdat%aqxy = rdat%aqx*rdat%acy
      rdat%aqz = acz+cqz
      rdat%qps = rdat%aqx2+rdat%acy2

! Use special fast routine for inner loops for 0000 ... 1111

      call spdgen(jtype,rdat,1)
    end do

    qx = rcd*sing
    qz = rcd*cosg

    call mcdv_all(grotspd, rdat,qx, qz, jtype)

    p = transpose(p)

    call r30s1d(jtype, grotspd, p)

  end subroutine genr22

  subroutine intclean(rdat,jtype)
    implicit none
    type(rotaxis_data_t) :: rdat
    integer, intent(in) :: jtype
    select case (jtype)
    case(1)
       rdat%r00(1,1) = 0.0_dp
    case(2)
       call intk_02(rdat,0)
    case(3)
       call intk_03(rdat,0)
    case(4)
       call intk_04(rdat,0)
    case(5)
       call intk_05(rdat,0)
    case(6)
       call intk_06(rdat,0)
    case(7)
       call intk_07(rdat,0)
    case(8)
       call intk_08(rdat,0)
    case(9)
       call intk_09(rdat,0)
    case(10)
       call intk_10(rdat,0)
    case(11)
       call intk_11(rdat,0)
    case(12)
       call intk_12(rdat,0)
    case(13)
       call intk_13(rdat,0)
    case(14)
       call intk_14(rdat,0)
    case(15)
       call intk_15(rdat,0)
    case(16)
       call intk_16(rdat,0)
    case(17)
       call intk_17(rdat,0)
    case(18)
       call intk_18(rdat,0)
    case(19)
       call intk_19(rdat,0)
    case(20)
       call intk_20(rdat,0)
    case(21)
       call intk_21(rdat,0)
    end select
  end subroutine

  subroutine mcdv_all(grotspd, rdat, qx, qz, jtype)
    implicit none
    type(rotaxis_data_t) :: rdat
    integer, intent(in) :: jtype
    real(kind=dp), intent(in) :: qx, qz
    real(kind=dp), intent(inout) :: grotspd(*)

    select case (jtype)
    case (1)
       grotspd(1) = rdat%r00 (1,1)
       return
    case (2)
       call mcdv_02(grotspd, rdat, qx, qz)
    case (3)
       call mcdv_03(grotspd, rdat, qx, qz)
    case (4)
       call mcdv_04(grotspd, rdat, qx, qz)
    case (5)
       call mcdv_05(grotspd, rdat, qx, qz)
    case (6)
       call mcdv_06(grotspd, rdat, qx, qz)
    case( 7)
       call mcdv_07(grotspd, rdat, qx, qz)
    case( 8)
       call mcdv_08(grotspd, rdat, qx, qz)
    case( 9)
       call mcdv_09(grotspd, rdat, qx, qz)
    case(10)
       call mcdv_10(grotspd, rdat, qx, qz)
    case(11)
       call mcdv_11(grotspd, rdat, qx, qz)
    case(12)
       call mcdv_12(grotspd, rdat, qx, qz)
    case(13)
       call mcdv_13(grotspd, rdat, qx, qz)
    case(14)
       call mcdv_14(grotspd, rdat, qx, qz)
    case(15)
       call mcdv_15(grotspd, rdat, qx, qz)
    case(16)
       call mcdv_16(grotspd, rdat, qx, qz)
    case(17)
       call mcdv_17(grotspd, rdat, qx, qz)
    case(18)
       call mcdv_18(grotspd, rdat, qx, qz)
    case(19)
       call mcdv_19(grotspd, rdat, qx, qz)
    case(20)
       call mcdv_20(grotspd, rdat, qx, qz)
    case(21)
       call mcdv_21(grotspd, rdat, qx, qz)
    end select

  end subroutine

! >    @brief   s,p,d rotated axis type selection
! >    @details s,p,d rotated axis type selection
 subroutine spdgen(jtype,rdat,ikl)

    implicit none
    type(rotaxis_data_t) :: rdat
    integer :: jtype, ikl

    select case(jtype)
    case (1)
       call intj_01(rdat)
       rdat%r00 (1,1) = rdat%r00 (1,1) + rdat%fq0 (1)
    case (2)
       call intj_02(rdat)
       call intk_02 (rdat,ikl)
    case(3)
       call intj_03(rdat)
       call intk_03 (rdat,ikl)
    case (4)
       call intj_04(rdat)
       call intk_04 (rdat,ikl)
    case (5)
       call intj_05(rdat)
       call intk_05 (rdat,ikl)
    case (6)
       call intj_06(rdat)
       call intk_06 (rdat,ikl)
    case(7)
       call intj_07(rdat)
       call intk_07(rdat,ikl)
    case(8)
       call intj_08(rdat)
       call intk_08(rdat,ikl)
    case(9)
       call intj_09(rdat)
       call intk_09(rdat,ikl)
    case(10)
       call intj_10(rdat)
       call intk_10(rdat,ikl)
    case(11)
       call intj_11(rdat)
       call intk_11(rdat,ikl)
    case(12)
       call intj_12(rdat)
       call intk_12(rdat,ikl)
    case(13)
       call intj_13(rdat)
       call intk_13(rdat,ikl)
    case(14)
       call intj_14(rdat)
       call intk_14(rdat,ikl)
    case(15)
       call intj_15(rdat)
       call intk_15(rdat,ikl)
    case(16)
       call intj_16(rdat)
       call intk_16(rdat,ikl)
    case(17)
       call intj_17(rdat)
       call intk_17(rdat,ikl)
    case(18)
       call intj_18(rdat)
       call intk_18(rdat,ikl)
    case(19)
       call intj_19(rdat)
       call intk_19(rdat,ikl)
    case(20)
       call intj_20(rdat)
       call intk_20(rdat,ikl)
    case(21)
       call intj_21(rdat)
       call intk_21(rdat,ikl)
    end select

  end subroutine

! >
! >    @brief   rotate up to 1296 s,p,d integrals to space fixed axes
! >
! >    @details rotate up to 1296 s,p,d integrals to space fixed axes
! >             incoming and outgoing integrals in f, while p(1,1),...
! >             are direction cosines of space fixed axes wrt axes at p
! >
      subroutine r30s1d(jtype,f,p)

      implicit none

      integer :: jtype
      real(kind=dp) :: f(*), p(3,3)

      select case (jtype)
      case (2)
        call r30s1d_02(f, p)
      case (3)
        call r30s1d_03(f, p)
      case (4)
        call r30s1d_04(f, p)
      case (5)
        call r30s1d_05(f, p)
      case (6)
        call r30s1d_06(f, p)
      case (7)
        call r30s1d_07(f, p)
      case (8)
        call r30s1d_08(f, p)
      case (9)
        call r30s1d_09(f, p)
      case (10)
        call r30s1d_10(f, p)
      case (11)
        call r30s1d_11(f, p)
      case (12)
        call r30s1d_12(f, p)
      case (13)
        call r30s1d_13(f, p)
      case (14)
        call r30s1d_14(f, p)
      case (15)
        call r30s1d_15(f, p)
      case (16)
        call r30s1d_16(f, p)
      case (17)
        call r30s1d_17(f, p)
      case (18)
        call r30s1d_18(f, p)
      case (19)
        call r30s1d_19(f, p)
      case (20)
        call r30s1d_20(f, p)
      case (21)
        call r30s1d_21(f, p)
      end select

      end subroutine r30s1d

! >
! >    @brief   ssss case
! >
! >    @details integration of a ssss case
! >
      subroutine intj_01(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype= 1 integrals
      integer :: i, n, ip
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin

      rdat%fq0 (1) = 0.0_dp

      do i = 1, rdat%ngangb
         fqz = rdat%sp (i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr = rdat%ty02(i) - rdat%aqz
         pqs = pqr * pqr
         rho = rdat%tx12(i) * rdat%x34 * x41
         if (rdat%lrint) then
            efr = rdat%emu2 / (rdat%emu2 + rho)
            rho = rho * efr
            fqz = fqz * sqrt (efr)
         endif
         xva = (pqs + rdat%qps) * rho
         rho = rho + rho
         n = 0
         if (xva<=tmax) then

! Fm(t) evaluation

            tv = xva * rfinc (n)
            ip = nint (tv)
            fx = fgrid (4, ip, n) * tv
            fx = (fx + fgrid (3, ip, n) ) * tv
            fx = (fx + fgrid (2, ip, n) ) * tv
            fx = (fx + fgrid (1, ip, n) ) * tv
            fx = fx + fgrid (0, ip, n)

            rdat%fq (n) = fx
! T2= xva+xva
! Do m=n,1,-1
! Rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
! End do

            fqf = fqz * sqrt (x41)
            rdat%fq (0) = rdat%fq (0) * fqf
! Do 210 m=0,n
! Rdat%fq(m)= rdat%fq(m)*fqf
! 210       fqf= fqf*rho
         else
            xin = 1.0_dp / xva
            rdat%fq (0) = fqz * sqrt (pi4 * xin * x41)
! Rox= rho*xin
! Fqf= 0.5_dp*rox
! Do 220 m=1,n
! Rdat%fq(m)= rdat%fq(m-1)*fqf
! 220       fqf= fqf+rox
         endif

         rdat%fq0 (1) = rdat%fq0 (1) + rdat%fq (0)

      end do

      end subroutine intj_01

! >
! >    @brief   psss case
! >
! >    @details integration of a psss case
! >
      subroutine intj_02(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype= 2 integrals
      integer :: i, n, ip
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox

      rdat%fq0 (1) = 0.0_dp
      rdat%fq1 (1,1) = 0.0_dp
      rdat%fq1 (2,1) = 0.0_dp

! Write(iw,*) 'rdat%ngangb',rdat%ngangb
      do i = 1, rdat%ngangb
         fqz = rdat%sp (i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr = rdat%ty02(i) - rdat%aqz
         pqs = pqr * pqr
         rho = rdat%tx12(i) * rdat%x34 * x41
         if (rdat%lrint) then
            efr = rdat%emu2 / (rdat%emu2 + rho)
            rho = rho * efr
            fqz = fqz * sqrt (efr)
         endif
         xva = (pqs + rdat%qps) * rho
         rho = rho + rho
         n = 1
         if (xva<=tmax) then

! Fm(t) evaluation

            tv = xva * rfinc (0)
            ip = nint (tv)
            fx = fgrid (4, ip, 0) * tv
            fx = (fx + fgrid (3, ip, 0) ) * tv
            fx = (fx + fgrid (2, ip, 0) ) * tv
            fx = (fx + fgrid (1, ip, 0) ) * tv
            fx = fx + fgrid (0, ip, 0)

            rdat%fq (0) = fx

            tv = xva * rfinc (n)
            ip = nint (tv)
            fx = fgrid (4, ip, n) * tv
            fx = (fx + fgrid (3, ip, n) ) * tv
            fx = (fx + fgrid (2, ip, n) ) * tv
            fx = (fx + fgrid (1, ip, n) ) * tv
            fx = fx + fgrid (0, ip, n)

            rdat%fq (n) = fx
! T2= xva+xva
! Do m=n,1,-1
! Rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
! End do

            fqf = fqz * sqrt (x41)
            rdat%fq (0) = rdat%fq (0) * fqf
            fqf = fqf * rho
            rdat%fq (1) = rdat%fq (1) * fqf
! Do 210 m=0,n
! Rdat%fq(m)= rdat%fq(m)*fqf
! 210       fqf= fqf*rho
         else
            xin = 1.0_dp / xva
            rdat%fq (0) = fqz * sqrt (pi4 * xin * x41)
            rox = rho * xin
            fqf = 0.5_dp * rox
            rdat%fq (1) = rdat%fq (0) * fqf
! Do 220 m=1,n
! Rdat%fq(m)= rdat%fq(m-1)*fqf
! 220       fqf= fqf+rox
         endif

         rdat%fq0 (1) = rdat%fq0 (1) + rdat%fq (0)
         rdat%fq1 (1,1) = rdat%fq1 (1,1) + rdat%fq (1)
         rdat%fq1 (2,1) = rdat%fq1 (2,1) + rdat%fq (1) * pqr
! Write(*,*) 'intj_02',i,rdat%fq0(1),rdat%fq(0)

      end do

      end subroutine intj_02

! >
! >    @brief   ppss case
! >
! >    @details integration of a ppss case
! >
      subroutine intj_03(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype= 3 integrals
      integer :: i, n, ip
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2

      rdat%fq0 (1) = 0.0_dp
      rdat%fq1 (:,1) = 0.0_dp
      rdat%fq2 (:,1) = 0.0_dp

      do i = 1, rdat%ngangb
         fqz = rdat%sp (i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr = rdat%ty02(i) - rdat%aqz
         pqs = pqr * pqr
         rho = rdat%tx12(i) * rdat%x34 * x41
         if (rdat%lrint) then
            efr = rdat%emu2 / (rdat%emu2 + rho)
            rho = rho * efr
            fqz = fqz * sqrt (efr)
         endif
         xva = (pqs + rdat%qps) * rho
         rho = rho + rho
         n = 2
         if (xva<=tmax) then

! Fm(t) evaluation

            tv = xva * rfinc (n)
            ip = nint (tv)
            fx = fgrid (4, ip, n) * tv
            fx = (fx + fgrid (3, ip, n) ) * tv
            fx = (fx + fgrid (2, ip, n) ) * tv
            fx = (fx + fgrid (1, ip, n) ) * tv
            fx = fx + fgrid (0, ip, n)
            tv = xva * rxinc
            ip = nint (tv)
            et = xgrid (4, ip) * tv
            et = (et + xgrid (3, ip) ) * tv
            et = (et + xgrid (2, ip) ) * tv
            et = (et + xgrid (1, ip) ) * tv
            et = et + xgrid (0, ip)

            rdat%fq (n) = fx
            t2 = xva + xva
            rdat%fq (2 - 1) = (t2 * rdat%fq (2) + et) * rmr (2)
            rdat%fq (1 - 1) = (t2 * rdat%fq (1) + et) * rmr (1)
! Do m=n,1,-1
! Rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
! End do

            fqf = fqz * sqrt (x41)
            rdat%fq (0) = rdat%fq (0) * fqf
            fqf = fqf * rho
            rdat%fq (1) = rdat%fq (1) * fqf
            fqf = fqf * rho
            rdat%fq (2) = rdat%fq (2) * fqf
! Do 210 m=0,n
! Rdat%fq(m)= rdat%fq(m)*fqf
! 210       fqf= fqf*rho
         else
            xin = 1.0_dp / xva
            rdat%fq (0) = fqz * sqrt (pi4 * xin * x41)
            rox = rho * xin
            fqf = 0.5_dp * rox
            rdat%fq (1) = rdat%fq (0) * fqf
            fqf = fqf + rox
            rdat%fq (2) = rdat%fq (1) * fqf
! Do 220 m=1,n
! Rdat%fq(m)= rdat%fq(m-1)*fqf
! 220       fqf= fqf+rox
         endif

         rdat%fq0 (1) = rdat%fq0 (1) + rdat%fq (0)
         rdat%fq1 (1,1) = rdat%fq1 (1,1) + rdat%fq (1)
         rdat%fq1 (2,1) = rdat%fq1 (2,1) + rdat%fq (1) * pqr
         rdat%fq2 (1,1) = rdat%fq2 (1,1) + rdat%fq (2)
         rdat%fq2 (2,1) = rdat%fq2 (2,1) + rdat%fq (2) * pqr
         rdat%fq2 (3,1) = rdat%fq2 (3,1) + rdat%fq (2) * pqs

      end do

      end subroutine intj_03

! >
! >    @brief   psps case
! >
! >    @details integration of a psps case
! >
      subroutine intj_04(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype= 4 integrals

      integer :: i, n, ip
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: xmd1, y01, tmp1, tmp2

      rdat%fq0 (1:2) = 0.0_dp
      rdat%fq1 (:,1:3) = 0.0_dp
      rdat%fq2 (1:3,1) = 0.0_dp

      do i = 1, rdat%ngangb
         fqz = rdat%sp (i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr = rdat%ty02(i) - rdat%aqz
         pqs = pqr * pqr
         rho = rdat%tx12(i) * rdat%x34 * x41
         if (rdat%lrint) then
            efr = rdat%emu2 / (rdat%emu2 + rho)
            rho = rho * efr
            fqz = fqz * sqrt (efr)
         endif
         xva = (pqs + rdat%qps) * rho
         rho = rho + rho
         n = 2
         if (xva<=tmax) then

! Fm(t) evaluation

            tv = xva * rfinc (n)
            ip = nint (tv)
            fx = fgrid (4, ip, n) * tv
            fx = (fx + fgrid (3, ip, n) ) * tv
            fx = (fx + fgrid (2, ip, n) ) * tv
            fx = (fx + fgrid (1, ip, n) ) * tv
            fx = fx + fgrid (0, ip, n)
            tv = xva * rxinc
            ip = nint (tv)
            et = xgrid (4, ip) * tv
            et = (et + xgrid (3, ip) ) * tv
            et = (et + xgrid (2, ip) ) * tv
            et = (et + xgrid (1, ip) ) * tv
            et = et + xgrid (0, ip)

            rdat%fq (n) = fx
            t2 = xva + xva
            rdat%fq (2 - 1) = (t2 * rdat%fq (2) + et) * rmr (2)
            rdat%fq (1 - 1) = (t2 * rdat%fq (1) + et) * rmr (1)
! Do m=n,1,-1
! Rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
! End do

            fqf = fqz * sqrt (x41)
            rdat%fq (0) = rdat%fq (0) * fqf
            fqf = fqf * rho
            rdat%fq (1) = rdat%fq (1) * fqf
            fqf = fqf * rho
            rdat%fq (2) = rdat%fq (2) * fqf
! Do 210 m=0,n
! Rdat%fq(m)= rdat%fq(m)*fqf
! 210       fqf= fqf*rho
         else
            xin = 1.0_dp / xva
            rdat%fq (0) = fqz * sqrt (pi4 * xin * x41)
            rox = rho * xin
            fqf = 0.5_dp * rox
            rdat%fq (1) = rdat%fq (0) * fqf
            fqf = fqf + rox
            rdat%fq (2) = rdat%fq (1) * fqf
! Do 220 m=1,n
! Rdat%fq(m)= rdat%fq(m-1)*fqf
! 220       fqf= fqf+rox
         endif

         xmd1 = 0.5_dp/rdat%tx12(i)
         y01 = rdat%ty02(i)-rdat%rab
         tmp1 = xmd1 * pqr
         tmp2 = xmd1 * pqs

         rdat%fq0 (1) = rdat%fq0 (1) + rdat%fq (0)
         rdat%fq0 (2) = rdat%fq0 (2) + rdat%fq (0) * y01
         rdat%fq1 (1,1) = rdat%fq1 (1,1) + rdat%fq (1)
         rdat%fq1 (2,1) = rdat%fq1 (2,1) + rdat%fq (1) * pqr
         rdat%fq1 (1,2) = rdat%fq1 (1,2) + rdat%fq (1) * xmd1
         rdat%fq1 (2,2) = rdat%fq1 (2,2) + rdat%fq (1) * tmp1
         rdat%fq1 (1,3) = rdat%fq1 (1,3) + rdat%fq (1) * y01
         rdat%fq1 (2,3) = rdat%fq1 (2,3) + rdat%fq (1) * y01 * pqr
         rdat%fq2 (1,1) = rdat%fq2 (1,1) + rdat%fq (2) * xmd1
         rdat%fq2 (2,1) = rdat%fq2 (2,1) + rdat%fq (2) * tmp1
         rdat%fq2 (3,1) = rdat%fq2 (3,1) + rdat%fq (2) * tmp2

      end do

      end subroutine intj_04

! >
! >    @brief   ppps case
! >
! >    @details integration of a ppps case
! >
      subroutine intj_05(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype= 5 integrals
      integer :: i, n, ip
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: xmd1, y01
      real(kind=dp) :: work (3, 3)

      rdat%fq0 (1:2) = 0.0_dp
      rdat%fq0 (2) = 0.0_dp
      rdat%fq1 (:, 1:3) = 0.0_dp
      rdat%fq2 (:, 1:3) = 0.0_dp
      rdat%fq3 (1:4,1) = 0.0_dp

      do i = 1, rdat%ngangb
         fqz = rdat%sp (i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr = rdat%ty02(i) - rdat%aqz
         pqs = pqr * pqr
         rho = rdat%tx12(i) * rdat%x34 * x41
         if (rdat%lrint) then
            efr = rdat%emu2 / (rdat%emu2 + rho)
            rho = rho * efr
            fqz = fqz * sqrt (efr)
         endif
         xva = (pqs + rdat%qps) * rho
         rho = rho + rho
         n = 3
         if (xva<=tmax) then

! Fm(t) evaluation...downward recursion for jtype >= 5

            tv = xva * rfinc (n)
            ip = nint (tv)
            fx = fgrid (4, ip, n) * tv
            fx = (fx + fgrid (3, ip, n) ) * tv
            fx = (fx + fgrid (2, ip, n) ) * tv
            fx = (fx + fgrid (1, ip, n) ) * tv
            fx = fx + fgrid (0, ip, n)
            tv = xva * rxinc
            ip = nint (tv)
            et = xgrid (4, ip) * tv
            et = (et + xgrid (3, ip) ) * tv
            et = (et + xgrid (2, ip) ) * tv
            et = (et + xgrid (1, ip) ) * tv
            et = et + xgrid (0, ip)

            rdat%fq (n) = fx
            t2 = xva + xva
            rdat%fq (3 - 1) = (t2 * rdat%fq (3) + et) * rmr (3)
            rdat%fq (2 - 1) = (t2 * rdat%fq (2) + et) * rmr (2)
            rdat%fq (1 - 1) = (t2 * rdat%fq (1) + et) * rmr (1)
! Do m=n,1,-1
! Rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
! End do

            fqf = fqz * sqrt (x41)
            rdat%fq (0) = rdat%fq (0) * fqf
            fqf = fqf * rho
            rdat%fq (1) = rdat%fq (1) * fqf
            fqf = fqf * rho
            rdat%fq (2) = rdat%fq (2) * fqf
            fqf = fqf * rho
            rdat%fq (3) = rdat%fq (3) * fqf
! Do 210 m=0,n
! Rdat%fq(m)= rdat%fq(m)*fqf
! 210       fqf= fqf*rho
         else
            xin = 1.0_dp / xva
            rdat%fq (0) = fqz * sqrt (pi4 * xin * x41)
            rox = rho * xin
            fqf = 0.5_dp * rox
            rdat%fq (1) = rdat%fq (0) * fqf
            fqf = fqf + rox
            rdat%fq (2) = rdat%fq (1) * fqf
            fqf = fqf + rox
            rdat%fq (3) = rdat%fq (2) * fqf
! Do 220 m=1,n
! Rdat%fq(m)= rdat%fq(m-1)*fqf
! 220       fqf= fqf+rox
         endif

         xmd1 = 0.5_dp/rdat%tx12(i)
         y01 = rdat%ty02(i)-rdat%rab
         work (1, 1) = 1.0d0
         work (1, 2) = y01
         work (1, 3) = xmd1
         work (2, 1:3) = work (1, 1:3) * pqr
         work (3, 1:3) = work (1, 1:3) * pqs

         rdat%fq0 (1:2)      = rdat%fq0 (1:2) + rdat%fq (0) * work (1, 1:2)
         rdat%fq1 (1:2, 1:3) = rdat%fq1 (1:2, 1:3) + rdat%fq (1) * work (1:2, 1:3)
         rdat%fq2 (1:3, 1:3) = rdat%fq2 (1:3, 1:3) + rdat%fq (2) * work (1:3, 1:3)
         rdat%fq3 (1:3,1)    = rdat%fq3 (1:3,1) + rdat%fq (3) * work (1:3, 3)
         rdat%fq3 (4,1)      = rdat%fq3 (4,1) + rdat%fq (3) * work (3, 3) * pqr

      end do

      end subroutine intj_05

! >
! >    @brief   pppp case
! >
! >    @details integration of a pppp case
! >
      subroutine intj_06(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype= 6 integrals
      integer :: i, n, ip, m
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: xmd1, y01
      real(kind=dp) :: pqt
      real(kind=dp) :: work (4, 9)

      rdat%fq0 = 0.0_dp
      rdat%fq1 = 0.0_dp
      rdat%fq2 = 0.0_dp
      rdat%fq3 = 0.0_dp
      rdat%fq4 = 0.0_dp

      do i = 1, rdat%ngangb
         fqz = rdat%sp (i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr = rdat%ty02(i) - rdat%aqz
         pqs = pqr * pqr
         rho = rdat%tx12(i) * rdat%x34 * x41
         if (rdat%lrint) then
            efr = rdat%emu2 / (rdat%emu2 + rho)
            rho = rho * efr
            fqz = fqz * sqrt (efr)
         endif
         xva = (pqs + rdat%qps) * rho
         rho = rho + rho
         n = 4
         if (xva<=tmax) then

! Fm(t) evaluation...downward recursion for jtype >= 5

            tv = xva * rfinc (n)
            ip = nint (tv)
            fx = fgrid (4, ip, n) * tv
            fx = (fx + fgrid (3, ip, n) ) * tv
            fx = (fx + fgrid (2, ip, n) ) * tv
            fx = (fx + fgrid (1, ip, n) ) * tv
            fx = fx + fgrid (0, ip, n)
            tv = xva * rxinc
            ip = nint (tv)
            et = xgrid (4, ip) * tv
            et = (et + xgrid (3, ip) ) * tv
            et = (et + xgrid (2, ip) ) * tv
            et = (et + xgrid (1, ip) ) * tv
            et = et + xgrid (0, ip)

            rdat%fq (n) = fx
            t2 = xva + xva
            do m = n, 1, - 1
            rdat%fq (m - 1) = (t2 * rdat%fq (m) + et) * rmr (m)
            enddo

            fqf = fqz * sqrt (x41)
            do m = 0, n
               rdat%fq (m) = rdat%fq (m) * fqf
               fqf = fqf * rho
            end do
         else
            xin = 1.0_dp / xva
            rdat%fq (0) = fqz * sqrt (pi4 * xin * x41)
            rox = rho * xin
            fqf = 0.5_dp * rox
            do m = 1, n
               rdat%fq (m) = rdat%fq (m - 1) * fqf
               fqf = fqf + rox
            end do
         endif

         pqt = pqr * pqs
         xmd1 = 0.5_dp/rdat%tx12(i)
         y01 = rdat%ty02(i)-rdat%rab
         work (1, 1) = 1.0d0
         work (1, 2) = y01
         work (1, 3) = rdat%ty02(i)
         work (1, 4) = y01 * rdat%ty02(i)
         work (1, 5) = xmd1
         work (1, 6) = xmd1
         work (1, 7) = xmd1
         work (1, 8) = xmd1 * y01
         work (1, 9) = xmd1 * xmd1
         work (2, 1:9) = work (1, 1:9) * pqr
         work (3, 1:9) = work (1, 1:9) * pqs
         work (4, 5:9) = work (1, 5:9) * pqt

         rdat%fq0 (1:5) = rdat%fq0 (1:5) + rdat%fq (0) * work (1, 1:5)
         rdat%fq1 (1:2, 1:8) = rdat%fq1 (1:2, 1:8) + rdat%fq (1) * work (1:2, 1:8)
         rdat%fq1 (1, 9) = rdat%fq1 (1, 9) + rdat%fq (1) * work (1, 9)

         rdat%fq2 (1:3, 1:9) = rdat%fq2 (1:3, 1:9) + rdat%fq (2) * work (1:3, 1:9)
         rdat%fq3 (1:4, 1:5) = rdat%fq3 (1:4, 1:5) + rdat%fq (3) * work (1:4, 5:9)

         rdat%fq4 (1:4,1) = rdat%fq4 (1:4,1) + rdat%fq (4) * work (1:4, 9)
         rdat%fq4 (5,1) = rdat%fq4 (5,1) + rdat%fq (4) * work (4, 9) * pqr

      end do

      end subroutine intj_06

! >
! >    @brief   dsss case
! >
! >    @details integration of the dsss case
! >
      subroutine intj_07(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype= 7 integrals
      integer :: i, n, ip
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2

      rdat%fq0(1)= 0.0_dp
      rdat%fq1(1:2,1)= 0.0_dp
      rdat%fq2(1:3,1)= 0.0_dp

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=2
         if(xva<=tmax) then

! Fm(t) evaluation

            tv= xva*rfinc(n)
            ip= nint(tv)
            fx=    fgrid(4,ip,n) *tv
            fx=(fx+fgrid(3,ip,n))*tv
            fx=(fx+fgrid(2,ip,n))*tv
            fx=(fx+fgrid(1,ip,n))*tv
            fx= fx+fgrid(0,ip,n)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(n)= fx
            t2= xva+xva
               rdat%fq(2-1)=(t2*rdat%fq(2)+et)*rmr(2)
               rdat%fq(1-1)=(t2*rdat%fq(1)+et)*rmr(1)
! Do m=n,1,-1
! Rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
! End do

            fqf= fqz*sqrt(x41)
               rdat%fq(0)= rdat%fq(0)*fqf
            fqf= fqf*rho
               rdat%fq(1)= rdat%fq(1)*fqf
            fqf= fqf*rho
               rdat%fq(2)= rdat%fq(2)*fqf
! Do 210 m=0,n
! Rdat%fq(m)= rdat%fq(m)*fqf
! 210       fqf= fqf*rho
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
               rdat%fq(1)= rdat%fq(0)*fqf
            fqf= fqf+rox
               rdat%fq(2)= rdat%fq(1)*fqf
! Do 220 m=1,n
! Rdat%fq(m)= rdat%fq(m-1)*fqf
! 220       fqf= fqf+rox
         endif

         rdat%fq0(1)= rdat%fq0(1)+rdat%fq(0)
         rdat%fq1(1,1)= rdat%fq1(1,1)+rdat%fq(1)
         rdat%fq1(2,1)= rdat%fq1(2,1)+rdat%fq(1)*pqr
         rdat%fq2(1,1)= rdat%fq2(1,1)+rdat%fq(2)
         rdat%fq2(2,1)= rdat%fq2(2,1)+rdat%fq(2)*pqr
         rdat%fq2(3,1)= rdat%fq2(3,1)+rdat%fq(2)*pqs
      end do

      end subroutine intj_07

! >
! >    @brief   dpss case
! >
! >    @details integration of the dpss case
! >
      subroutine intj_08(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype= 8 integrals
      integer :: i, n, ip
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2

      rdat%fq0(1) = 0.0_dp
      rdat%fq1(1:2,1) = 0.0_dp
      rdat%fq2(1:3,1) = 0.0_dp
      rdat%fq3(1:4,1) = 0.0_dp

      do i = 1, rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=3
         if(xva<=tmax) then

! Fm(t) evaluation...downward recursion for m=3

            tv= xva*rfinc(n)
            ip= nint(tv)
            fx=    fgrid(4,ip,n) *tv
            fx=(fx+fgrid(3,ip,n))*tv
            fx=(fx+fgrid(2,ip,n))*tv
            fx=(fx+fgrid(1,ip,n))*tv
            fx= fx+fgrid(0,ip,n)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(n)= fx
            t2= xva+xva
               rdat%fq(3-1)=(t2*rdat%fq(3)+et)*rmr(3)
               rdat%fq(2-1)=(t2*rdat%fq(2)+et)*rmr(2)
               rdat%fq(1-1)=(t2*rdat%fq(1)+et)*rmr(1)
! Do m=n,1,-1
! Rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
! End do

            fqf= fqz*sqrt(x41)
               rdat%fq(0)= rdat%fq(0)*fqf
            fqf= fqf*rho
               rdat%fq(1)= rdat%fq(1)*fqf
            fqf= fqf*rho
               rdat%fq(2)= rdat%fq(2)*fqf
            fqf= fqf*rho
               rdat%fq(3)= rdat%fq(3)*fqf
! Do 210 m=0,n
! Rdat%fq(m)= rdat%fq(m)*fqf
! 210       fqf= fqf*rho
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
               rdat%fq(1)= rdat%fq(0)*fqf
            fqf= fqf+rox
               rdat%fq(2)= rdat%fq(1)*fqf
            fqf= fqf+rox
               rdat%fq(3)= rdat%fq(2)*fqf
! Do 220 m=1,n
! Rdat%fq(m)= rdat%fq(m-1)*fqf
! 220       fqf= fqf+rox
         endif

         rdat%fq0(1)= rdat%fq0(1)+rdat%fq(0)
         rdat%fq1(1,1)= rdat%fq1(1,1)+rdat%fq(1)
         rdat%fq1(2,1)= rdat%fq1(2,1)+rdat%fq(1)*pqr
         rdat%fq2(1,1)= rdat%fq2(1,1)+rdat%fq(2)
         rdat%fq2(2,1)= rdat%fq2(2,1)+rdat%fq(2)*pqr
         rdat%fq2(3,1)= rdat%fq2(3,1)+rdat%fq(2)*pqs
         rdat%fq3(1,1)= rdat%fq3(1,1)+rdat%fq(3)
         rdat%fq3(2,1)= rdat%fq3(2,1)+rdat%fq(3)*pqr
         rdat%fq3(3,1)= rdat%fq3(3,1)+rdat%fq(3)*pqs
         rdat%fq3(4,1)= rdat%fq3(4,1)+rdat%fq(3)*pqs*pqr
      end do

      end subroutine intj_08

! >
! >    @brief   dsps case
! >
! >    @details integration of the dsps case
! >
      subroutine intj_09(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype= 9 integrals
      integer :: i, n, ip
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: work

      dimension  work(3,3)

      rdat%fq0(1:2) = 0.0_dp
      rdat%fq1(1:2,1:3) = 0.0_dp
      rdat%fq2(1:3,1:3) = 0.0_dp
      rdat%fq3(1:4,1) = 0.0_dp

      do i = 1, rdat%ngangb
         fqz = rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=3
         if(xva<=tmax) then

! Fm(t) evaluation...downward recursion for m=3

            tv= xva*rfinc(n)
            ip= nint(tv)
            fx=    fgrid(4,ip,n) *tv
            fx=(fx+fgrid(3,ip,n))*tv
            fx=(fx+fgrid(2,ip,n))*tv
            fx=(fx+fgrid(1,ip,n))*tv
            fx= fx+fgrid(0,ip,n)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(n)= fx
            t2= xva+xva
               rdat%fq(3-1)=(t2*rdat%fq(3)+et)*rmr(3)
               rdat%fq(2-1)=(t2*rdat%fq(2)+et)*rmr(2)
               rdat%fq(1-1)=(t2*rdat%fq(1)+et)*rmr(1)
! Do m=n,1,-1
! Rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
! End do

            fqf= fqz*sqrt(x41)
               rdat%fq(0)= rdat%fq(0)*fqf
            fqf= fqf*rho
               rdat%fq(1)= rdat%fq(1)*fqf
            fqf= fqf*rho
               rdat%fq(2)= rdat%fq(2)*fqf
            fqf= fqf*rho
               rdat%fq(3)= rdat%fq(3)*fqf
! Do 210 m=0,n
! Rdat%fq(m)= rdat%fq(m)*fqf
! 210       fqf= fqf*rho
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
               rdat%fq(1)= rdat%fq(0)*fqf
            fqf= fqf+rox
               rdat%fq(2)= rdat%fq(1)*fqf
            fqf= fqf+rox
               rdat%fq(3)= rdat%fq(2)*fqf
! Do 220 m=1,n
! Rdat%fq(m)= rdat%fq(m-1)*fqf
! 220       fqf= fqf+rox
         endif

         work(1,1)= 1.0d0
         work(1,2)= rdat%ty02(i)-rdat%rab
         work(1,3)= 0.5_dp/rdat%tx12(i)
         work(2,1:3)= work(1,1:3)*pqr
         work(3,1:3)= work(1,1:3)*pqs

         rdat%fq0(1:2) = rdat%fq0(1:2)+rdat%fq(0)*work(1,1:2)
         rdat%fq1(1:2,1:3) = rdat%fq1(1:2,1:3)+rdat%fq(1)*work(1:2,1:3)
         rdat%fq2(1:3,1:3) = rdat%fq2(1:3,1:3)+rdat%fq(2)*work(1:3,1:3)

         rdat%fq3(1:3,1) = rdat%fq3(1:3,1)+rdat%fq(3)*work(1:3,3)
         rdat%fq3(4,1) = rdat%fq3(4,1)+rdat%fq(3)*work(3,3)*pqr
      end do

      end subroutine intj_09

! >
! >    @brief   ddss case
! >
! >    @details integration of the ddss case
! >
      subroutine intj_10(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype=10 integrals
      integer :: i, n, ip, m
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: pqq, pqt

      rdat%fq0(1)= 0.0_dp
      rdat%fq1(1:2,1)= 0.0_dp
      rdat%fq2(1:3,1)= 0.0_dp
      rdat%fq3(1:4,1)= 0.0_dp
      rdat%fq4(1:5,1)= 0.0_dp

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=4
         if(xva<=tmax) then

! Fm(t) evaluation...downward recursion for m=4

            tv= xva*rfinc(n)
            ip= nint(tv)
            fx=    fgrid(4,ip,n) *tv
            fx=(fx+fgrid(3,ip,n))*tv
            fx=(fx+fgrid(2,ip,n))*tv
            fx=(fx+fgrid(1,ip,n))*tv
            fx= fx+fgrid(0,ip,n)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(n)= fx
            t2= xva+xva
            do m=n,1,-1
               rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
            end do

            fqf= fqz*sqrt(x41)
            do m=0,n
               rdat%fq(m)= rdat%fq(m)*fqf
               fqf= fqf*rho
            end do
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
            do m=1,n
               rdat%fq(m)= rdat%fq(m-1)*fqf
               fqf= fqf+rox
            end do
         endif

         pqt= pqs*pqr
         pqq= pqs*pqs

         rdat%fq0(1)= rdat%fq0(1)+rdat%fq(0)
         rdat%fq1(1,1)= rdat%fq1(1,1)+rdat%fq(1)
         rdat%fq1(2,1)= rdat%fq1(2,1)+rdat%fq(1)*pqr
         rdat%fq2(1,1)= rdat%fq2(1,1)+rdat%fq(2)
         rdat%fq2(2,1)= rdat%fq2(2,1)+rdat%fq(2)*pqr
         rdat%fq2(3,1)= rdat%fq2(3,1)+rdat%fq(2)*pqs
         rdat%fq3(1,1)= rdat%fq3(1,1)+rdat%fq(3)
         rdat%fq3(2,1)= rdat%fq3(2,1)+rdat%fq(3)*pqr
         rdat%fq3(3,1)= rdat%fq3(3,1)+rdat%fq(3)*pqs
         rdat%fq3(4,1)= rdat%fq3(4,1)+rdat%fq(3)*pqt
         rdat%fq4(1,1)= rdat%fq4(1,1)+rdat%fq(4)
         rdat%fq4(2,1)= rdat%fq4(2,1)+rdat%fq(4)*pqr
         rdat%fq4(3,1)= rdat%fq4(3,1)+rdat%fq(4)*pqs
         rdat%fq4(4,1)= rdat%fq4(4,1)+rdat%fq(4)*pqt
         rdat%fq4(5,1)= rdat%fq4(5,1)+rdat%fq(4)*pqq
      end do

      end subroutine intj_10

! >
! >    @brief   dpps case
! >
! >    @details integration of the dpps case
! >
      subroutine intj_11(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype=11 integrals
      integer :: i, n, ip, j, m
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: pqt, work

      dimension  work(4,3)

      rdat%fq0(1:2)= 0.0_dp
      rdat%fq1(:,1:3)= 0.0_dp
      rdat%fq2(:,1:3)= 0.0_dp
      rdat%fq3(:,1:3)= 0.0_dp
      rdat%fq4(:,1)= 0.0_dp

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=4
         if(xva<=tmax) then

! Fm(t) evaluation...downward recursion for m=4

            tv= xva*rfinc(n)
            ip= nint(tv)
            fx=    fgrid(4,ip,n) *tv
            fx=(fx+fgrid(3,ip,n))*tv
            fx=(fx+fgrid(2,ip,n))*tv
            fx=(fx+fgrid(1,ip,n))*tv
            fx= fx+fgrid(0,ip,n)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(n)= fx
            t2= xva+xva
            do m=n,1,-1
               rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
            end do

            fqf= fqz*sqrt(x41)
            do m=0,n
               rdat%fq(m)= rdat%fq(m)*fqf
               fqf= fqf*rho
            end do
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
            do m=1,n
               rdat%fq(m)= rdat%fq(m-1)*fqf
               fqf= fqf+rox
            end do
         endif

         pqt = pqr*pqs
         work(1,1)= 1.0d0
         work(1,2)= rdat%ty02(i)-rdat%rab
         work(1,3)= 0.5_dp/rdat%tx12(i)
         do j= 1, 3
            work(2,j)= work(1,j)*pqr
            work(3,j)= work(1,j)*pqs
            work(4,j)= work(1,j)*pqt
         enddo

         rdat%fq0(1)= rdat%fq0(1)+rdat%fq(0)*work(1,1)
         rdat%fq0(2)= rdat%fq0(2)+rdat%fq(0)*work(1,2)
         rdat%fq1(:,1:3)= rdat%fq1(:,1:3)+rdat%fq(1)*work(1:2,1:3)
         rdat%fq2(:,1:3)= rdat%fq2(:,1:3)+rdat%fq(2)*work(1:3,1:3)
         rdat%fq3(:,1:3)= rdat%fq3(:,1:3)+rdat%fq(3)*work(:,1:3)
         rdat%fq4(1:4,1)= rdat%fq4(1:4,1)+rdat%fq(4)*work(1:4,3)
         rdat%fq4(5,1)= rdat%fq4(5,1)+rdat%fq(4)*work(4,3)*pqr
      end do

      end subroutine intj_11

! >
! >    @brief   dsds case
! >
! >    @details integration of the dsds case
! >
      subroutine intj_12(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype=12 integrals
      integer :: i, n, ip, j, m
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: pqt, xmd1, y01, work

      dimension  work(4,4)

         rdat%fq0(1:2)= 0.0_dp
         rdat%fq1(:,1:4)= 0.0_dp
         rdat%fq2(:,1:4)= 0.0_dp
         rdat%fq3(:,1:2)= 0.0_dp
         rdat%fq4(:,1)= 0.0_dp

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=4
         if(xva<=tmax) then

! Fm(t) evaluation...downward recursion for m=4

            tv= xva*rfinc(n)
            ip= nint(tv)
            fx=    fgrid(4,ip,n) *tv
            fx=(fx+fgrid(3,ip,n))*tv
            fx=(fx+fgrid(2,ip,n))*tv
            fx=(fx+fgrid(1,ip,n))*tv
            fx= fx+fgrid(0,ip,n)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(n)= fx
            t2= xva+xva
            do m=n,1,-1
               rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
            end do

            fqf= fqz*sqrt(x41)
            do m=0,n
               rdat%fq(m)= rdat%fq(m)*fqf
               fqf= fqf*rho
            end do
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
            do m=1,n
               rdat%fq(m)= rdat%fq(m-1)*fqf
               fqf= fqf+rox
            end do
         endif

         pqt = pqr*pqs
         xmd1= 0.5_dp/rdat%tx12(i)
         y01 = rdat%ty02(i)-rdat%rab
         work(1,1)= xmd1
         work(1,2)= y01 *y01
         work(1,3)= xmd1*y01
         work(1,4)= xmd1*xmd1
         do j= 1, 4
            work(2,j)= work(1,j)*pqr
            work(3,j)= work(1,j)*pqs
         enddo
         work(4,3)= work(1,3)*pqt
         work(4,4)= work(1,4)*pqt

         rdat%fq0(1:2)= rdat%fq0(1:2)+rdat%fq(0)*work(1,1:2)
         rdat%fq1(:,1:3)= rdat%fq1(:,1:3)+rdat%fq(1)*work(1:2,1:3)
         rdat%fq1(1,4)= rdat%fq1(1,4)+rdat%fq(1)*work(1,4)
         rdat%fq2(:,1:4)= rdat%fq2(:,1:4)+rdat%fq(2)*work(1:3,1:4)
         rdat%fq3(:,1:2)= rdat%fq3(:,1:2)+rdat%fq(3)*work(:,3:4)
         rdat%fq4(1:4,1)= rdat%fq4(1:4,1)+rdat%fq(4)*work(1:4,4)
         rdat%fq4(5,1)= rdat%fq4(5,1)+rdat%fq(4)*work(4,4)*pqr
      end do

      end subroutine intj_12

! >
! >    @brief   dspp case
! >
! >    @details integration of the dspp case
! >
      subroutine intj_13(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype=13 integrals
      integer :: i, n, ip, m
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: pqt, xmd1, y01, work

      dimension  work(4,9)

      rdat%fq0(:)= 0.0_dp
      rdat%fq1(:,1:9)= 0.0_dp
      rdat%fq2(:,1:9)= 0.0_dp
      rdat%fq3(:,1:5)= 0.0_dp
      rdat%fq4(1:5,1)= 0.0_dp

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=4
         if(xva<=tmax) then

! Fm(t) evaluation...downward recursion for m=4

            tv= xva*rfinc(n)
            ip= nint(tv)
            fx=    fgrid(4,ip,n) *tv
            fx=(fx+fgrid(3,ip,n))*tv
            fx=(fx+fgrid(2,ip,n))*tv
            fx=(fx+fgrid(1,ip,n))*tv
            fx= fx+fgrid(0,ip,n)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(n)= fx
            t2= xva+xva
            do m=n,1,-1
               rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
            end do

            fqf= fqz*sqrt(x41)
            do m=0,n
               rdat%fq(m)= rdat%fq(m)*fqf
               fqf= fqf*rho
            end do
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
            do m=1,n
               rdat%fq(m)= rdat%fq(m-1)*fqf
               fqf= fqf+rox
            end do
         endif

        pqt = pqr*pqs
        xmd1= 0.5_dp/rdat%tx12(i)
        y01 = rdat%ty02(i)-rdat%rab
        work(1,1)= 1.0d0
        work(1,2)= y01
        work(1,3)= rdat%ty02(i)
        work(1,4)= y01 *rdat%ty02(i)
        work(1,5)= xmd1
        work(1,6)= xmd1
        work(1,7)= xmd1
        work(1,8)= xmd1*y01
        work(1,9)= xmd1*xmd1
        work(2,1:9)= work(1,1:9)*pqr
        work(3,1:9)= work(1,1:9)*pqs
        work(4,5:9)= work(1,5:9)*pqt
        rdat%fq0(1:5)= rdat%fq0(1:5)+rdat%fq(0)*work(1,1:5)
        rdat%fq1(:,1:8)= rdat%fq1(:,1:8)+rdat%fq(1)*work(1:2,1:8)
        rdat%fq1(1,9)= rdat%fq1(1,9)+rdat%fq(1)*work(1,9)
        rdat%fq2(:,1:9)= rdat%fq2(:,1:9)+rdat%fq(2)*work(1:3,1:9)
        rdat%fq3(:,1:5)= rdat%fq3(:,1:5)+rdat%fq(3)*work(1:4,5:9)
        rdat%fq4(1:4,1)= rdat%fq4(1:4,1)+rdat%fq(4)*work(1:4,9)
        rdat%fq4(5,1)= rdat%fq4(5,1)+rdat%fq(4)*work(4,9)*pqr
      end do
      end subroutine intj_13

! >
! >    @brief   ddps case
! >
! >    @details integration of the ddps case
! >
      subroutine intj_14(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype=14 integrals
      integer :: i, n, ip, j, m
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: pqt, pqq, work

      dimension  work(5,3)

      rdat%fq0(1:2)= 0.0_dp
      rdat%fq1(:,1:3)= 0.0_dp
      rdat%fq2(:,1:3)= 0.0_dp
      rdat%fq3(:,1:3)= 0.0_dp
      rdat%fq4(:,1:3)= 0.0_dp
      rdat%fq5(:,1)= 0.0_dp

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=5
         if(xva<=tmax) then

! Fm(t) m=5 interpolation, generating wasted m=8,7,6 data
! Fgrid(,,x) for  x=0,1,2,3,4,5, 6, 7 holds necessary data to
! Interpolate for m=0,1,2,3,4,8,12,16.
! Here m=5, so we must generate m=8,7,6 values we don't use.
! Downward recursion is used for greater numerical stability.
! Note that we also use an interpolation for exp(-t) here.

            tv= xva*rfinc(n)
            ip= nint(tv)
            fx=    fgrid(4,ip,n) *tv
            fx=(fx+fgrid(3,ip,n))*tv
            fx=(fx+fgrid(2,ip,n))*tv
            fx=(fx+fgrid(1,ip,n))*tv
            fx= fx+fgrid(0,ip,n)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(8)= fx
            t2= xva+xva
            do m=8,1,-1
               rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
            end do

! This is other parts of the integral, not fm(t)

            fqf= fqz*sqrt(x41)
            do m=0,n
               rdat%fq(m)= rdat%fq(m)*fqf
               fqf= fqf*rho
            end do
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
            do m=1,n
               rdat%fq(m)= rdat%fq(m-1)*fqf
               fqf= fqf+rox
            end do
         endif

         pqt = pqr*pqs
         pqq = pqs*pqs
         work(1,1)= 1.0d0
         work(1,2)= rdat%ty02(i)-rdat%rab
         work(1,3)= 0.5_dp/rdat%tx12(i)
         do j= 1, 3
            work(2,j)= work(1,j)*pqr
            work(3,j)= work(1,j)*pqs
            work(4,j)= work(1,j)*pqt
            work(5,j)= work(1,j)*pqq
         enddo

         rdat%fq0(1:2)= rdat%fq0(1:2)+rdat%fq(0)*work(1,1:2)
         rdat%fq1(:,1:3)= rdat%fq1(:,1:3)+rdat%fq(1)*work(1:2,1:3)
         rdat%fq2(:,1:3)= rdat%fq2(:,1:3)+rdat%fq(2)*work(1:3,1:3)
         rdat%fq3(:,1:3)= rdat%fq3(:,1:3)+rdat%fq(3)*work(1:4,1:3)
         rdat%fq4(:,1:3)= rdat%fq4(:,1:3)+rdat%fq(4)*work(1:5,1:3)
         rdat%fq5(1:5,1)= rdat%fq5(1:5,1)+rdat%fq(5)*work(1:5,3)
         rdat%fq5(6,1)= rdat%fq5(6,1)+rdat%fq(5)*work(5,3)*pqr
      end do

      end subroutine intj_14

! >
! >    @brief   dpds case
! >
! >    @details integration of the dpds case
! >
      subroutine intj_15(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype=15 integrals
      integer :: i, n, ip, j, m
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: pqt, xmd1, y01, work

      dimension  work(5,4)

         rdat%fq0(1:2)= 0.0_dp
         rdat%fq1(:,1:4)= 0.0_dp
         rdat%fq2(:,1:4)= 0.0_dp
         rdat%fq3(:,1:4)= 0.0_dp
         rdat%fq4(:,1:2)= 0.0_dp
         rdat%fq5(:,1)= 0.0_dp

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=5
         if(xva<=tmax) then

! Fm(t) m=5 interpolation, generating wasted m=8,7,6 data

            tv= xva*rfinc(n)
            ip= nint(tv)
            fx=    fgrid(4,ip,n) *tv
            fx=(fx+fgrid(3,ip,n))*tv
            fx=(fx+fgrid(2,ip,n))*tv
            fx=(fx+fgrid(1,ip,n))*tv
            fx= fx+fgrid(0,ip,n)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(8)= fx
            t2= xva+xva
            do m=8,1,-1
               rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
            end do

            fqf= fqz*sqrt(x41)
            do m=0,n
               rdat%fq(m)= rdat%fq(m)*fqf
               fqf= fqf*rho
            end do
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
            do m=1,n
               rdat%fq(m)= rdat%fq(m-1)*fqf
               fqf= fqf+rox
            end do
         endif

         pqt = pqr*pqs
         xmd1= 0.5_dp/rdat%tx12(i)
         y01 = rdat%ty02(i)-rdat%rab
         work(1,1)= xmd1
         work(1,2)= y01 *y01
         work(1,3)= xmd1*y01
         work(1,4)= xmd1*xmd1
         do j= 1, 4
            work(2,j)= work(1,j)*pqr
            work(3,j)= work(1,j)*pqs
            work(4,j)= work(1,j)*pqt
         enddo
         work(5,3)= work(4,3)*pqr
         work(5,4)= work(4,4)*pqr

         rdat%fq0(1)= rdat%fq0(1)+rdat%fq(0)*work(1,1)
         rdat%fq0(2)= rdat%fq0(2)+rdat%fq(0)*work(1,2)
         do j= 1, 3
            rdat%fq1(1,j)= rdat%fq1(1,j)+rdat%fq(1)*work(1,j)
            rdat%fq1(2,j)= rdat%fq1(2,j)+rdat%fq(1)*work(2,j)
         enddo
            rdat%fq1(1,4)= rdat%fq1(1,4)+rdat%fq(1)*work(1,4)
         do j= 1, 4
            rdat%fq2(1,j)= rdat%fq2(1,j)+rdat%fq(2)*work(1,j)
            rdat%fq2(2,j)= rdat%fq2(2,j)+rdat%fq(2)*work(2,j)
            rdat%fq2(3,j)= rdat%fq2(3,j)+rdat%fq(2)*work(3,j)

            rdat%fq3(1,j)= rdat%fq3(1,j)+rdat%fq(3)*work(1,j)
            rdat%fq3(2,j)= rdat%fq3(2,j)+rdat%fq(3)*work(2,j)
            rdat%fq3(3,j)= rdat%fq3(3,j)+rdat%fq(3)*work(3,j)
            rdat%fq3(4,j)= rdat%fq3(4,j)+rdat%fq(3)*work(4,j)
         enddo
         do j= 1, 2
            rdat%fq4(1,j)= rdat%fq4(1,j)+rdat%fq(4)*work(1,j+ 2)
            rdat%fq4(2,j)= rdat%fq4(2,j)+rdat%fq(4)*work(2,j+ 2)
            rdat%fq4(3,j)= rdat%fq4(3,j)+rdat%fq(4)*work(3,j+ 2)
            rdat%fq4(4,j)= rdat%fq4(4,j)+rdat%fq(4)*work(4,j+ 2)
            rdat%fq4(5,j)= rdat%fq4(5,j)+rdat%fq(4)*work(5,j+ 2)
         enddo
         rdat%fq5(1:5,1)= rdat%fq5(1:5,1)+rdat%fq(5)*work(1:5,4)
         rdat%fq5(6,1)= rdat%fq5(6,1)+rdat%fq(5)*work(5,4)*pqr
      end do

      end subroutine intj_15

! >
! >    @brief   dppp case
! >
! >    @details integration of the dppp case
! >
      subroutine intj_16(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype=16 integrals
      integer :: i, n, ip, j, m, k
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: pqt, xmd1, y01, pqq, work

      dimension  work(5,9)

      rdat%fq0(1:5)= 0.0_dp
      rdat%fq1(:,1:9)= 0.0_dp
      rdat%fq2(:,1:9)= 0.0_dp
      rdat%fq3(:,1:9)= 0.0_dp
      rdat%fq4(:,1:5)= 0.0_dp
      rdat%fq5(:,1)= 0.0_dp

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=5
         if(xva<=tmax) then

! Fm(t) m=5 interpolation, generating wasted m=8,7,6 data

            tv= xva*rfinc(n)
            ip= nint(tv)
            fx=    fgrid(4,ip,n) *tv
            fx=(fx+fgrid(3,ip,n))*tv
            fx=(fx+fgrid(2,ip,n))*tv
            fx=(fx+fgrid(1,ip,n))*tv
            fx= fx+fgrid(0,ip,n)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(8)= fx
            t2= xva+xva
            do m=8,1,-1
               rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
            end do

            fqf= fqz*sqrt(x41)
            do m=0,n
               rdat%fq(m)= rdat%fq(m)*fqf
               fqf= fqf*rho
            end do
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
            do m=1,n
               rdat%fq(m)= rdat%fq(m-1)*fqf
               fqf= fqf+rox
            end do
         endif

         pqt = pqr*pqs
         pqq = pqs*pqs
         xmd1= 0.5_dp/rdat%tx12(i)
         y01 = rdat%ty02(i)-rdat%rab
         work(1,1)= 1.0d0
         work(1,2)= y01
         work(1,3)= rdat%ty02(i)
         work(1,4)= y01 *rdat%ty02(i)
         work(1,5)= xmd1
         work(1,6)= xmd1
         work(1,7)= xmd1
         work(1,8)= xmd1*y01
         work(1,9)= xmd1*xmd1
         do j= 1, 9
            work(2,j)= work(1,j)*pqr
            work(3,j)= work(1,j)*pqs
            work(4,j)= work(1,j)*pqt
         enddo
         do j= 5, 9
            work(5,j)= work(1,j)*pqq
         enddo

         do j= 1, 5
            rdat%fq0(j)= rdat%fq0(j)+rdat%fq(0)*work(1,j)
         enddo
         do j= 1, 8
            rdat%fq1(1,j)= rdat%fq1(1,j)+rdat%fq(1)*work(1,j)
            rdat%fq1(2,j)= rdat%fq1(2,j)+rdat%fq(1)*work(2,j)
         enddo
            rdat%fq1(1,9)= rdat%fq1(1,9)+rdat%fq(1)*work(1,9)
         do j= 1, 9
            rdat%fq2(1,j)= rdat%fq2(1,j)+rdat%fq(2)*work(1,j)
            rdat%fq2(2,j)= rdat%fq2(2,j)+rdat%fq(2)*work(2,j)
            rdat%fq2(3,j)= rdat%fq2(3,j)+rdat%fq(2)*work(3,j)

            rdat%fq3(1,j)= rdat%fq3(1,j)+rdat%fq(3)*work(1,j)
            rdat%fq3(2,j)= rdat%fq3(2,j)+rdat%fq(3)*work(2,j)
            rdat%fq3(3,j)= rdat%fq3(3,j)+rdat%fq(3)*work(3,j)
            rdat%fq3(4,j)= rdat%fq3(4,j)+rdat%fq(3)*work(4,j)
         enddo
         do j= 1, 5
            rdat%fq4(1,j)= rdat%fq4(1,j)+rdat%fq(4)*work(1,j+ 4)
            rdat%fq4(2,j)= rdat%fq4(2,j)+rdat%fq(4)*work(2,j+ 4)
            rdat%fq4(3,j)= rdat%fq4(3,j)+rdat%fq(4)*work(3,j+ 4)
            rdat%fq4(4,j)= rdat%fq4(4,j)+rdat%fq(4)*work(4,j+ 4)
            rdat%fq4(5,j)= rdat%fq4(5,j)+rdat%fq(4)*work(5,j+ 4)
         enddo
         do k= 1, 5
            rdat%fq5(k,1)= rdat%fq5(k,1)+rdat%fq(5)*work(k,9)
         enddo
            rdat%fq5(6,1)= rdat%fq5(6,1)+rdat%fq(5)*work(5,9)*pqr
      end do

      end subroutine intj_16

! >
! >    @brief   ddds case
! >
! >    @details integration of the ddds case
! >
      subroutine intj_17(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype=17 integrals
      integer :: i, n, ip, j, m, k
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: pqt, xmd1, y01, pqq, work

      dimension  work(6,4)

         rdat%fq0(1)= 0.0_dp
         rdat%fq0(2)= 0.0_dp
         do j= 1, 4
            rdat%fq1(1,j)= 0.0_dp
            rdat%fq1(2,j)= 0.0_dp

            rdat%fq2(1,j)= 0.0_dp
            rdat%fq2(2,j)= 0.0_dp
            rdat%fq2(3,j)= 0.0_dp

            rdat%fq3(1,j)= 0.0_dp
            rdat%fq3(2,j)= 0.0_dp
            rdat%fq3(3,j)= 0.0_dp
            rdat%fq3(4,j)= 0.0_dp

            rdat%fq4(1,j)= 0.0_dp
            rdat%fq4(2,j)= 0.0_dp
            rdat%fq4(3,j)= 0.0_dp
            rdat%fq4(4,j)= 0.0_dp
            rdat%fq4(5,j)= 0.0_dp
         enddo
         do j= 1, 2
            do k= 1, 6
               rdat%fq5(k,j)= 0.0_dp
            enddo
         enddo
         rdat%fq6(:,1)= 0.0_dp

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=6
         if(xva<=tmax) then

! Fm(t) m=6 interpolation, generating wasted m=8,7 data

            m=5
            tv= xva*rfinc(m)
            ip= nint(tv)
            fx=    fgrid(4,ip,m) *tv
            fx=(fx+fgrid(3,ip,m))*tv
            fx=(fx+fgrid(2,ip,m))*tv
            fx=(fx+fgrid(1,ip,m))*tv
            fx= fx+fgrid(0,ip,m)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(8)= fx
            t2= xva+xva
            do m=8,1,-1
               rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
            end do

            fqf= fqz*sqrt(x41)
            do m=0,n
               rdat%fq(m)= rdat%fq(m)*fqf
               fqf= fqf*rho
            end do
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
            do m=1,n
               rdat%fq(m)= rdat%fq(m-1)*fqf
               fqf= fqf+rox
            end do
         endif

         pqt = pqr*pqs
         pqq = pqs*pqs
         xmd1= 0.5_dp/rdat%tx12(i)
         y01 = rdat%ty02(i)-rdat%rab
         work(1,1)= xmd1
         work(1,2)= y01 *y01
         work(1,3)= xmd1*y01
         work(1,4)= xmd1*xmd1
         do j= 1, 4
            work(2,j)= work(1,j)*pqr
            work(3,j)= work(1,j)*pqs
            work(4,j)= work(1,j)*pqt
            work(5,j)= work(1,j)*pqq
         enddo
         work(6,3)= work(5,3)*pqr
         work(6,4)= work(5,4)*pqr

         rdat%fq0(1)= rdat%fq0(1)+rdat%fq(0)*work(1,1)
         rdat%fq0(2)= rdat%fq0(2)+rdat%fq(0)*work(1,2)
         do j= 1, 3
            rdat%fq1(1,j)= rdat%fq1(1,j)+rdat%fq(1)*work(1,j)
            rdat%fq1(2,j)= rdat%fq1(2,j)+rdat%fq(1)*work(2,j)
         enddo
            rdat%fq1(1,4)= rdat%fq1(1,4)+rdat%fq(1)*work(1,4)
         do j= 1, 4
            rdat%fq2(1,j)= rdat%fq2(1,j)+rdat%fq(2)*work(1,j)
            rdat%fq2(2,j)= rdat%fq2(2,j)+rdat%fq(2)*work(2,j)
            rdat%fq2(3,j)= rdat%fq2(3,j)+rdat%fq(2)*work(3,j)

            rdat%fq3(1,j)= rdat%fq3(1,j)+rdat%fq(3)*work(1,j)
            rdat%fq3(2,j)= rdat%fq3(2,j)+rdat%fq(3)*work(2,j)
            rdat%fq3(3,j)= rdat%fq3(3,j)+rdat%fq(3)*work(3,j)
            rdat%fq3(4,j)= rdat%fq3(4,j)+rdat%fq(3)*work(4,j)

            rdat%fq4(1,j)= rdat%fq4(1,j)+rdat%fq(4)*work(1,j)
            rdat%fq4(2,j)= rdat%fq4(2,j)+rdat%fq(4)*work(2,j)
            rdat%fq4(3,j)= rdat%fq4(3,j)+rdat%fq(4)*work(3,j)
            rdat%fq4(4,j)= rdat%fq4(4,j)+rdat%fq(4)*work(4,j)
            rdat%fq4(5,j)= rdat%fq4(5,j)+rdat%fq(4)*work(5,j)
         enddo
         do j= 1, 2
            do k= 1, 6
               rdat%fq5(k,j)= rdat%fq5(k,j)+rdat%fq(5)*work(k,j+ 2)
            enddo
         enddo
         do k= 1, 6
            rdat%fq6(k,1)= rdat%fq6(k,1)+rdat%fq(6)*work(k,4)
         enddo
            rdat%fq6(7,1)= rdat%fq6(7,1)+rdat%fq(6)*work(6,4)*pqr
      end do

      end subroutine intj_17

! >
! >    @brief   ddpp case
! >
! >    @details integration of the ddpp case
! >
      subroutine intj_18(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype=18 integrals
      integer :: i, n, ip, j, m, k
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: pqt, xmd1, y01, pqq, pq5, work

      dimension  work(6,9)

         do j= 1, 5
            rdat%fq0(j)= 0.0_dp
         enddo
         do j= 1, 9
            rdat%fq1(1,j)= 0.0_dp
            rdat%fq1(2,j)= 0.0_dp

            rdat%fq2(1,j)= 0.0_dp
            rdat%fq2(2,j)= 0.0_dp
            rdat%fq2(3,j)= 0.0_dp

            rdat%fq3(1,j)= 0.0_dp
            rdat%fq3(2,j)= 0.0_dp
            rdat%fq3(3,j)= 0.0_dp
            rdat%fq3(4,j)= 0.0_dp

            rdat%fq4(1,j)= 0.0_dp
            rdat%fq4(2,j)= 0.0_dp
            rdat%fq4(3,j)= 0.0_dp
            rdat%fq4(4,j)= 0.0_dp
            rdat%fq4(5,j)= 0.0_dp
         enddo
         do j= 1, 5
            do k= 1, 6
               rdat%fq5(k,j)= 0.0_dp
            enddo
         enddo
         rdat%fq6(:,1)= 0.0_dp

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=6
         if(xva<=tmax) then

! Fm(t) m=6 interpolation, generating wasted m=8,7 data

            m=5
            tv= xva*rfinc(m)
            ip= nint(tv)
            fx=    fgrid(4,ip,m) *tv
            fx=(fx+fgrid(3,ip,m))*tv
            fx=(fx+fgrid(2,ip,m))*tv
            fx=(fx+fgrid(1,ip,m))*tv
            fx= fx+fgrid(0,ip,m)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(8)= fx
            t2= xva+xva
            do m=8,1,-1
               rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
            end do

            fqf= fqz*sqrt(x41)
            do m=0,n
               rdat%fq(m)= rdat%fq(m)*fqf
               fqf= fqf*rho
            end do
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
            do m=1,n
               rdat%fq(m)= rdat%fq(m-1)*fqf
               fqf= fqf+rox
            end do
         endif

         pqt = pqr*pqs
         pqq = pqs*pqs
         pq5 = pqt*pqs
         xmd1= 0.5_dp/rdat%tx12(i)
         y01 = rdat%ty02(i)-rdat%rab
         work(1,1)= 1.0d0
         work(1,2)= y01
         work(1,3)= rdat%ty02(i)
         work(1,4)= y01 *rdat%ty02(i)
         work(1,5)= xmd1
         work(1,6)= xmd1
         work(1,7)= xmd1
         work(1,8)= xmd1*y01
         work(1,9)= xmd1*xmd1
         do j= 1, 9
            work(2,j)= work(1,j)*pqr
            work(3,j)= work(1,j)*pqs
            work(4,j)= work(1,j)*pqt
            work(5,j)= work(1,j)*pqq
         enddo
         do j= 5, 9
            work(6,j)= work(1,j)*pq5
         enddo

         do j= 1, 5
            rdat%fq0(j)= rdat%fq0(j)+rdat%fq(0)*work(1,j)
         enddo
         do j= 1, 8
            rdat%fq1(1,j)= rdat%fq1(1,j)+rdat%fq(1)*work(1,j)
            rdat%fq1(2,j)= rdat%fq1(2,j)+rdat%fq(1)*work(2,j)
         enddo
            rdat%fq1(1,9)= rdat%fq1(1,9)+rdat%fq(1)*work(1,9)
         do j= 1, 9
            rdat%fq2(1,j)= rdat%fq2(1,j)+rdat%fq(2)*work(1,j)
            rdat%fq2(2,j)= rdat%fq2(2,j)+rdat%fq(2)*work(2,j)
            rdat%fq2(3,j)= rdat%fq2(3,j)+rdat%fq(2)*work(3,j)

            rdat%fq3(1,j)= rdat%fq3(1,j)+rdat%fq(3)*work(1,j)
            rdat%fq3(2,j)= rdat%fq3(2,j)+rdat%fq(3)*work(2,j)
            rdat%fq3(3,j)= rdat%fq3(3,j)+rdat%fq(3)*work(3,j)
            rdat%fq3(4,j)= rdat%fq3(4,j)+rdat%fq(3)*work(4,j)

            rdat%fq4(1,j)= rdat%fq4(1,j)+rdat%fq(4)*work(1,j)
            rdat%fq4(2,j)= rdat%fq4(2,j)+rdat%fq(4)*work(2,j)
            rdat%fq4(3,j)= rdat%fq4(3,j)+rdat%fq(4)*work(3,j)
            rdat%fq4(4,j)= rdat%fq4(4,j)+rdat%fq(4)*work(4,j)
            rdat%fq4(5,j)= rdat%fq4(5,j)+rdat%fq(4)*work(5,j)
         enddo
         do j= 1, 5
            do k= 1, 6
               rdat%fq5(k,j)= rdat%fq5(k,j)+rdat%fq(5)*work(k,j+ 4)
            enddo
         enddo
         do k= 1, 6
            rdat%fq6(k,1)= rdat%fq6(k,1)+rdat%fq(6)*work(k,9)
         enddo
            rdat%fq6(7,1)= rdat%fq6(7,1)+rdat%fq(6)*work(6,9)*pqr
      end do

      end subroutine intj_18

! >
! >    @brief   dpdp case
! >
! >    @details integration of the dpdp case
! >
      subroutine intj_19(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype=19 integrals
      integer :: i, n, ip, j, m, k
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: pqt, xmd1, y01, pqq, xmd2, y11, pq5, work

      dimension  work(6,11)

         do j= 1, 5
            rdat%fq0(j)= 0.0_dp
         enddo
         do j= 1,10
            rdat%fq1(1,j)= 0.0_dp
            rdat%fq1(2,j)= 0.0_dp
         enddo
         do j= 1,11
            rdat%fq2(1,j)= 0.0_dp
            rdat%fq2(2,j)= 0.0_dp
            rdat%fq2(3,j)= 0.0_dp

            rdat%fq3(1,j)= 0.0_dp
            rdat%fq3(2,j)= 0.0_dp
            rdat%fq3(3,j)= 0.0_dp
            rdat%fq3(4,j)= 0.0_dp
         enddo
         do j= 1, 7
            rdat%fq4(1,j)= 0.0_dp
            rdat%fq4(2,j)= 0.0_dp
            rdat%fq4(3,j)= 0.0_dp
            rdat%fq4(4,j)= 0.0_dp
            rdat%fq4(5,j)= 0.0_dp
         enddo
         do j= 1, 4
            do k= 1, 6
               rdat%fq5(k,j)= 0.0_dp
            enddo
         enddo
         rdat%fq6(:,1)= 0.0_dp

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=6
         if(xva<=tmax) then

! Fm(t) m=6 interpolation, generating wasted m=8,7 data

            m=5
            tv= xva*rfinc(m)
            ip= nint(tv)
            fx=    fgrid(4,ip,m) *tv
            fx=(fx+fgrid(3,ip,m))*tv
            fx=(fx+fgrid(2,ip,m))*tv
            fx=(fx+fgrid(1,ip,m))*tv
            fx= fx+fgrid(0,ip,m)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(8)= fx
            t2= xva+xva
            do m=8,1,-1
               rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
            end do

            fqf= fqz*sqrt(x41)
            do m=0,n
               rdat%fq(m)= rdat%fq(m)*fqf
               fqf= fqf*rho
            end do
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
            do m=1,n
               rdat%fq(m)= rdat%fq(m-1)*fqf
               fqf= fqf+rox
            end do
         endif

         pqt = pqr*pqs
         pqq = pqs*pqs
         pq5 = pqt*pqs
         xmd1= 0.5_dp/rdat%tx12(i)
         xmd2= xmd1*xmd1
         y01 = rdat%ty02(i)-rdat%rab
         y11 = y01 *y01
         work(1, 1)= xmd1
         work(1, 2)= y11
         work(1, 3)= y11 *rdat%ty02(i)
         work(1, 4)= xmd1*rdat%ty02(i)
         work(1, 5)= xmd1*y01
         work(1, 6)= xmd1*y01
         work(1, 7)= xmd1*y11
         work(1, 8)= xmd2
         work(1, 9)= xmd2
         work(1,10)= xmd2*y01
         work(1,11)= xmd2*xmd1
         do j= 1,11
            work(2,j)= work(1,j)*pqr
            work(3,j)= work(1,j)*pqs
            work(4,j)= work(1,j)*pqt
         enddo
         do j= 5,11
            work(5,j)= work(1,j)*pqq
         enddo
         do j= 8,11
            work(6,j)= work(1,j)*pq5
         enddo

         do j= 1, 5
            rdat%fq0(j)= rdat%fq0(j)+rdat%fq(0)*work(1,j)
         enddo
         do j= 1, 8
            rdat%fq1(1,j)= rdat%fq1(1,j)+rdat%fq(1)*work(1,j)
            rdat%fq1(2,j)= rdat%fq1(2,j)+rdat%fq(1)*work(2,j)
         enddo
            rdat%fq1(1, 9)= rdat%fq1(1, 9)+rdat%fq(1)*work(1, 9)
            rdat%fq1(1,10)= rdat%fq1(1,10)+rdat%fq(1)*work(1,10)
         do j= 1,11
            rdat%fq2(1,j)= rdat%fq2(1,j)+rdat%fq(2)*work(1,j)
            rdat%fq2(2,j)= rdat%fq2(2,j)+rdat%fq(2)*work(2,j)
            rdat%fq2(3,j)= rdat%fq2(3,j)+rdat%fq(2)*work(3,j)

            rdat%fq3(1,j)= rdat%fq3(1,j)+rdat%fq(3)*work(1,j)
            rdat%fq3(2,j)= rdat%fq3(2,j)+rdat%fq(3)*work(2,j)
            rdat%fq3(3,j)= rdat%fq3(3,j)+rdat%fq(3)*work(3,j)
            rdat%fq3(4,j)= rdat%fq3(4,j)+rdat%fq(3)*work(4,j)
         enddo
         do j= 1, 7
            rdat%fq4(1,j)= rdat%fq4(1,j)+rdat%fq(4)*work(1,j+4)
            rdat%fq4(2,j)= rdat%fq4(2,j)+rdat%fq(4)*work(2,j+4)
            rdat%fq4(3,j)= rdat%fq4(3,j)+rdat%fq(4)*work(3,j+4)
            rdat%fq4(4,j)= rdat%fq4(4,j)+rdat%fq(4)*work(4,j+4)
            rdat%fq4(5,j)= rdat%fq4(5,j)+rdat%fq(4)*work(5,j+4)
         enddo
         do j= 1, 4
            do k= 1, 6
               rdat%fq5(k,j)= rdat%fq5(k,j)+rdat%fq(5)*work(k,j+ 7)
            enddo
         enddo
         do k= 1, 6
            rdat%fq6(k,1)= rdat%fq6(k,1)+rdat%fq(6)*work(k,11)
         enddo
            rdat%fq6(7,1)= rdat%fq6(7,1)+rdat%fq(6)*work(6,11)*pqr
      end do

      end subroutine intj_19

! >
! >    @brief   dddp case
! >
! >    @details integration of the dddp case
! >
      subroutine intj_20(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype=20 integrals
      integer :: i, n, ip, j, m, k
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: pqt, xmd1, y01, pqq, xmd2, y11, pq5, pq6, work

      dimension  work(7,11)

         do j= 1, 5
            rdat%fq0(j)= 0.0_dp
         enddo
         do j= 1,10
            rdat%fq1(1,j)= 0.0_dp
            rdat%fq1(2,j)= 0.0_dp
         enddo
         do j= 1,11
            rdat%fq2(1,j)= 0.0_dp
            rdat%fq2(2,j)= 0.0_dp
            rdat%fq2(3,j)= 0.0_dp

            rdat%fq3(1,j)= 0.0_dp
            rdat%fq3(2,j)= 0.0_dp
            rdat%fq3(3,j)= 0.0_dp
            rdat%fq3(4,j)= 0.0_dp

            rdat%fq4(1,j)= 0.0_dp
            rdat%fq4(2,j)= 0.0_dp
            rdat%fq4(3,j)= 0.0_dp
            rdat%fq4(4,j)= 0.0_dp
            rdat%fq4(5,j)= 0.0_dp
         enddo
         do j= 1, 7
            do k= 1, 6
               rdat%fq5(k,j)= 0.0_dp
            enddo
         enddo
         do j= 1, 4
            do k= 1, 7
               rdat%fq6(k,j)= 0.0_dp
            enddo
         enddo
         rdat%fq7(:,1)= 0.0_dp

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=7
         if(xva<=tmax) then

! Fm(t) m=7 interpolation, generating wasted m=8 data

            m=5
            tv= xva*rfinc(m)
            ip= nint(tv)
            fx=    fgrid(4,ip,m) *tv
            fx=(fx+fgrid(3,ip,m))*tv
            fx=(fx+fgrid(2,ip,m))*tv
            fx=(fx+fgrid(1,ip,m))*tv
            fx= fx+fgrid(0,ip,m)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(8)= fx
            t2= xva+xva
            do m=8,1,-1
               rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
            end do

            fqf= fqz*sqrt(x41)
            do m=0,n
               rdat%fq(m)= rdat%fq(m)*fqf
               fqf= fqf*rho
            end do
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
            do m=1,n
               rdat%fq(m)= rdat%fq(m-1)*fqf
               fqf= fqf+rox
            end do
         endif

         pqt = pqr*pqs
         pqq = pqs*pqs
         pq5 = pqt*pqs
         pq6 = pqq*pqs
         xmd1= 0.5_dp/rdat%tx12(i)
         xmd2= xmd1*xmd1
         y01 = rdat%ty02(i)-rdat%rab
         y11 = y01 *y01
         work(1, 1)= xmd1
         work(1, 2)= y11
         work(1, 3)= y11 *rdat%ty02(i)
         work(1, 4)= xmd1*rdat%ty02(i)
         work(1, 5)= xmd1*y01
         work(1, 6)= xmd1*y01
         work(1, 7)= xmd1*y11
         work(1, 8)= xmd2
         work(1, 9)= xmd2
         work(1,10)= xmd2*y01
         work(1,11)= xmd2*xmd1
         do j= 1,11
            work(2,j)= work(1,j)*pqr
            work(3,j)= work(1,j)*pqs
            work(4,j)= work(1,j)*pqt
            work(5,j)= work(1,j)*pqq
         enddo
         do j= 5,11
            work(6,j)= work(1,j)*pq5
         enddo
         do j= 8,11
            work(7,j)= work(1,j)*pq6
         enddo

         do j= 1, 5
            rdat%fq0(j)= rdat%fq0(j)+rdat%fq(0)*work(1,j)
         enddo
         do j= 1, 8
            rdat%fq1(1,j)= rdat%fq1(1,j)+rdat%fq(1)*work(1,j)
            rdat%fq1(2,j)= rdat%fq1(2,j)+rdat%fq(1)*work(2,j)
         enddo
            rdat%fq1(1, 9)= rdat%fq1(1, 9)+rdat%fq(1)*work(1, 9)
            rdat%fq1(1,10)= rdat%fq1(1,10)+rdat%fq(1)*work(1,10)
         do j= 1,11
            rdat%fq2(1,j)= rdat%fq2(1,j)+rdat%fq(2)*work(1,j)
            rdat%fq2(2,j)= rdat%fq2(2,j)+rdat%fq(2)*work(2,j)
            rdat%fq2(3,j)= rdat%fq2(3,j)+rdat%fq(2)*work(3,j)

            rdat%fq3(1,j)= rdat%fq3(1,j)+rdat%fq(3)*work(1,j)
            rdat%fq3(2,j)= rdat%fq3(2,j)+rdat%fq(3)*work(2,j)
            rdat%fq3(3,j)= rdat%fq3(3,j)+rdat%fq(3)*work(3,j)
            rdat%fq3(4,j)= rdat%fq3(4,j)+rdat%fq(3)*work(4,j)

            rdat%fq4(1,j)= rdat%fq4(1,j)+rdat%fq(4)*work(1,j)
            rdat%fq4(2,j)= rdat%fq4(2,j)+rdat%fq(4)*work(2,j)
            rdat%fq4(3,j)= rdat%fq4(3,j)+rdat%fq(4)*work(3,j)
            rdat%fq4(4,j)= rdat%fq4(4,j)+rdat%fq(4)*work(4,j)
            rdat%fq4(5,j)= rdat%fq4(5,j)+rdat%fq(4)*work(5,j)
         enddo
         do j= 1, 7
            do k= 1, 6
               rdat%fq5(k,j)= rdat%fq5(k,j)+rdat%fq(5)*work(k,j+ 4)
            enddo
         enddo
         do j= 1, 4
            do k= 1, 7
               rdat%fq6(k,j)= rdat%fq6(k,j)+rdat%fq(6)*work(k,j+ 7)
            enddo
         enddo
         do k= 1, 7
            rdat%fq7(k,1)= rdat%fq7(k,1)+rdat%fq(7)*work(k,11)
         enddo
            rdat%fq7(8,1)= rdat%fq7(8,1)+rdat%fq(7)*work(7,11)*pqr
      end do

      end subroutine intj_20

! >
! >    @brief   dddd case
! >
! >    @details integration of the dddd case
! >
      subroutine intj_21(rdat)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype=21 integrals
      integer :: i, n, ip, j, m, k
      real(kind=dp) :: fqz, x41, pqr, pqs, rho, xva
      real(kind=dp) :: tv, fx, efr, fqf, xin, rox, et, t2
      real(kind=dp) :: xmd1, y01, xmd2, y11, pqwk, xmd3, y12, y22, work

      dimension pqwk(2:8),work(8,16)

         do j= 1, 5
            rdat%fq0(j)= 0.0_dp
         enddo
         do j= 1,13
            rdat%fq1(1,j)= 0.0_dp
            rdat%fq1(2,j)= 0.0_dp
         enddo
         do j= 1,16
            rdat%fq2(1,j)= 0.0_dp
            rdat%fq2(2,j)= 0.0_dp
            rdat%fq2(3,j)= 0.0_dp

            rdat%fq3(1,j)= 0.0_dp
            rdat%fq3(2,j)= 0.0_dp
            rdat%fq3(3,j)= 0.0_dp
            rdat%fq3(4,j)= 0.0_dp

            rdat%fq4(1,j)= 0.0_dp
            rdat%fq4(2,j)= 0.0_dp
            rdat%fq4(3,j)= 0.0_dp
            rdat%fq4(4,j)= 0.0_dp
            rdat%fq4(5,j)= 0.0_dp
         enddo
         do j= 1,11
            do k= 1, 6
               rdat%fq5(k,j)= 0.0_dp
            enddo
         enddo
         do j= 1, 7
            do k= 1, 7
               rdat%fq6(k,j)= 0.0_dp
            enddo
         enddo
         do j= 1, 3
            do k= 1, 8
               rdat%fq7(k,j)= 0.0_dp
            enddo
         enddo
         do k= 1, 9
            rdat%fq8(k)= 0.0_dp
         enddo

      do i=1,rdat%ngangb
         fqz= rdat%sp(i) * rdat%sq
         x41 = rdat%tx12(i) + rdat%x34
         if (fqz*fqz<rdat%cutoff*x41) cycle
         x41 = 1 / x41
         pqr= rdat%ty02(i)-rdat%aqz
         pqs= pqr*pqr
         rho= rdat%tx12(i)*rdat%x34*x41
         if(rdat%lrint) then
            efr= rdat%emu2/(rdat%emu2+rho)
            rho= rho*efr
            fqz= fqz*sqrt(efr)
         endif
         xva=(pqs+rdat%qps)*rho
         rho= rho+rho
         n=8
         if(xva<=tmax) then

! Fm(t) m=8 interpolation

            m=5
            tv= xva*rfinc(m)
            ip= nint(tv)
            fx=    fgrid(4,ip,m) *tv
            fx=(fx+fgrid(3,ip,m))*tv
            fx=(fx+fgrid(2,ip,m))*tv
            fx=(fx+fgrid(1,ip,m))*tv
            fx= fx+fgrid(0,ip,m)
            tv= xva*rxinc
            ip= nint(tv)
            et=    xgrid(4,ip) *tv
            et=(et+xgrid(3,ip))*tv
            et=(et+xgrid(2,ip))*tv
            et=(et+xgrid(1,ip))*tv
            et= et+xgrid(0,ip)

            rdat%fq(8)= fx
            t2= xva+xva
            do m=8,1,-1
               rdat%fq(m-1)=(t2*rdat%fq(m)+et)*rmr(m)
            end do

            fqf= fqz*sqrt(x41)
            do m=0,n
               rdat%fq(m)= rdat%fq(m)*fqf
               fqf= fqf*rho
            end do
         else
            xin= 1.0_dp/xva
            rdat%fq(0)= fqz*sqrt(pi4*xin*x41)
            rox= rho*xin
            fqf= 0.5_dp*rox
            do m=1,n
               rdat%fq(m)= rdat%fq(m-1)*fqf
               fqf= fqf+rox
            end do
         endif

         pqwk(2)= pqr
         pqwk(3)= pqs
         pqwk(4)= pqs*pqr
         pqwk(5)= pqs*pqs
         pqwk(6)= pqwk(4)*pqs
         pqwk(7)= pqwk(5)*pqs
         pqwk(8)= pqwk(6)*pqs
         xmd1= 0.5_dp/rdat%tx12(i)
         xmd2= xmd1*xmd1
         xmd3= xmd2*xmd1
         y01 = rdat%ty02(i)-rdat%rab
         y11 = y01 *y01
         y12 = y01 *rdat%ty02(i)
         y22 = rdat%ty02(i) *rdat%ty02(i)
         work(1, 1)= xmd2
         work(1, 2)= xmd1*y11
         work(1, 3)= xmd1*y12
         work(1, 4)= xmd1*y22
         work(1, 5)= y11 *y22
         work(1, 6)= xmd2*y01
         work(1, 7)= xmd2*rdat%ty02(i)
         work(1, 8)= xmd1*y11*rdat%ty02(i)
         work(1, 9)= xmd1*y12*rdat%ty02(i)
         work(1,10)= xmd3
         work(1,11)= xmd2*y11
         work(1,12)= xmd2*y12
         work(1,13)= xmd2*y22
         work(1,14)= xmd3*y01
         work(1,15)= xmd3*rdat%ty02(i)
         work(1,16)= xmd3*xmd1
         do j= 1, 5
            do k= 2, 5
               work(k,j)= work(1,j)*pqwk(k)
            enddo
         enddo
         do j= 6, 9
            do k= 2, 6
               work(k,j)= work(1,j)*pqwk(k)
            enddo
         enddo
         do j=10,13
            do k= 2, 7
               work(k,j)= work(1,j)*pqwk(k)
            enddo
         enddo
         do j=14,16
            do k= 2, 8
               work(k,j)= work(1,j)*pqwk(k)
            enddo
         enddo

         do j= 1, 5
            rdat%fq0(j)= rdat%fq0(j)+rdat%fq(0)*work(1,j)
         enddo
         do j= 1,13
            rdat%fq1(1,j)= rdat%fq1(1,j)+rdat%fq(1)*work(1,j)
            rdat%fq1(2,j)= rdat%fq1(2,j)+rdat%fq(1)*work(2,j)
         enddo
         do j= 1,16
            rdat%fq2(1,j)= rdat%fq2(1,j)+rdat%fq(2)*work(1,j)
            rdat%fq2(2,j)= rdat%fq2(2,j)+rdat%fq(2)*work(2,j)
            rdat%fq2(3,j)= rdat%fq2(3,j)+rdat%fq(2)*work(3,j)

            rdat%fq3(1,j)= rdat%fq3(1,j)+rdat%fq(3)*work(1,j)
            rdat%fq3(2,j)= rdat%fq3(2,j)+rdat%fq(3)*work(2,j)
            rdat%fq3(3,j)= rdat%fq3(3,j)+rdat%fq(3)*work(3,j)
            rdat%fq3(4,j)= rdat%fq3(4,j)+rdat%fq(3)*work(4,j)

            rdat%fq4(1,j)= rdat%fq4(1,j)+rdat%fq(4)*work(1,j)
            rdat%fq4(2,j)= rdat%fq4(2,j)+rdat%fq(4)*work(2,j)
            rdat%fq4(3,j)= rdat%fq4(3,j)+rdat%fq(4)*work(3,j)
            rdat%fq4(4,j)= rdat%fq4(4,j)+rdat%fq(4)*work(4,j)
            rdat%fq4(5,j)= rdat%fq4(5,j)+rdat%fq(4)*work(5,j)
         enddo
         do j= 1,11
            do k= 1, 6
               rdat%fq5(k,j)= rdat%fq5(k,j)+rdat%fq(5)*work(k,j+ 5)
            enddo
         enddo
         do j= 1, 7
            do k= 1, 7
               rdat%fq6(k,j)= rdat%fq6(k,j)+rdat%fq(6)*work(k,j+ 9)
            enddo
         enddo
         do j= 1, 3
            do k= 1, 8
               rdat%fq7(k,j)= rdat%fq7(k,j)+rdat%fq(7)*work(k,j+13)
            enddo
         enddo
         do k= 1, 8
            rdat%fq8(k)= rdat%fq8(k)+rdat%fq(8)*work(k,16)
         enddo
            rdat%fq8(9)= rdat%fq8(9)+rdat%fq(8)*work(8,16)*pqr
      end do

      end subroutine intj_21

! >
! >    @brief   psss case
! >
! >    @details integration of a psss case
! >
      subroutine intk_02 (rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd2, xmdt

! Generate jtype= 2 integrals


      if (ikl == 0) then
         rdat%r00 (1:2,1) = 0.0_dp
         rdat%r01 (1, 1) = 0.0_dp
         rdat%r01 (2, 1) = 0.0_dp
         rdat%r01 (3, 1) = 0.0_dp

         return
      endif

      xmd2 = rdat%x43 * 0.5d+00
      xmdt = xmd2

      rdat%r00 (1,1) = rdat%r00 (1,1) + rdat%fq0 (1)
      rdat%r00 (2,1) = rdat%r00 (2,1) - rdat%fq0 (1) * rdat%y03
! Write(iw,*) 'intk_02 rdat%r00',rdat%r00(1),rdat%r00(2)
! Write(iw,*) 'rdat%fq0 rdat%sq',rdat%fq0(1),rdat%sq
! Write(iw,*) 'rdat%r00',rdat%sq,rdat%y03

      rdat%r01 (1, 1) = rdat%r01 (1, 1) - rdat%fq1 (1, 1) * rdat%aqx * xmdt
      rdat%r01 (2, 1) = rdat%r01 (2, 1) - rdat%fq1 (1, 1) * rdat%acy * xmdt
      rdat%r01 (3, 1) = rdat%r01 (3, 1) + rdat%fq1 (2, 1) * xmdt
      end subroutine intk_02

! >
! >    @brief   ppss case
! >
! >    @details integration of a ppss case
! >
      subroutine intk_03 (rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat

! Generate jtype= 3 integrals

      integer :: ikl,i
      real(dp) :: fq11,fq12,fq13, work(4), xmd2,xmd3,xmdt

      if (ikl == 0) then
         rdat%r00 (1:5,1) = 0.0_dp
         rdat%r01 (1:3,1:4) = 0.0_dp
         rdat%r02 (1:6,1) = 0.0_dp
         return
      endif

      xmd2 = rdat%x43 * 0.5d+00
      xmd3 = xmd2
      xmdt = xmd3 * xmd2
      work (1) = xmd2
      work (2) = xmd2
      work (3) = - xmd3 * rdat%y03
      work (4) = xmd3 * rdat%y04

      rdat%r00 (1,1) = rdat%r00 (1,1) + rdat%fq0 (1)
      rdat%r00 (2,1) = rdat%r00 (2,1) - rdat%fq0 (1) * rdat%y03
      rdat%r00 (3,1) = rdat%r00 (3,1) + rdat%fq0 (1) * rdat%y04
      rdat%r00 (4,1) = rdat%r00 (4,1) + rdat%fq0 (1) * xmd3
      rdat%r00 (5,1) = rdat%r00 (5,1) - rdat%fq0 (1) * rdat%y03 * rdat%y04

      fq11 = - rdat%fq1 (1, 1) * rdat%aqx
      fq12 = - rdat%fq1 (1, 1) * rdat%acy
      fq13 = rdat%fq1 (2, 1)
      do i = 1, 4
      rdat%r01 (1, i) = rdat%r01 (1, i) + fq11 * work (i)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fq12 * work (i)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fq13 * work (i)
      enddo

      rdat%r02 (1,1) = rdat%r02 (1,1) + (rdat%fq2 (1,1) * rdat%aqx2 - rdat%fq1 (1, 1) ) * xmdt
      rdat%r02 (2,1) = rdat%r02 (2,1) + (rdat%fq2 (1,1) * rdat%acy2 - rdat%fq1 (1, 1) ) * xmdt
      rdat%r02 (3,1) = rdat%r02 (3,1) + (rdat%fq2 (3,1) - rdat%fq1 (1, 1) ) * xmdt
      rdat%r02 (4,1) = rdat%r02 (4,1) + rdat%fq2 (1,1) * rdat%aqxy * xmdt
      rdat%r02 (5,1) = rdat%r02 (5,1) - rdat%fq2 (2,1) * rdat%aqx * xmdt
      rdat%r02 (6,1) = rdat%r02 (6,1) - rdat%fq2 (2,1) * rdat%acy * xmdt

      end subroutine intk_03

! >
! >    @brief   psps case
! >
! >    @details integration of a psps case
! >
      subroutine intk_04 (rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd2
      real(kind=dp) :: work
      integer :: i


! Generate jtype= 4 integrals

      dimension work (4)

      if (ikl == 0) then
         rdat%r00 (1:4,1) = 0.0_dp
         rdat%r01 (1:3,1:4) = 0.0_dp
         rdat%r02 (1:6,1) = 0.0_dp
         return
      endif

      xmd2 = rdat%x43 * 0.5d+00
      work (1) = xmd2
      work (2) = 1.0d0
      work (3) = work (1)
      work (4) = - rdat%y03

      rdat%r00 (1,1) = rdat%r00 (1,1) + rdat%fq0 (1) * work (2)
      rdat%r00 (2,1) = rdat%r00 (2,1) + rdat%fq0 (1) * work (4)
      rdat%r00 (3,1) = rdat%r00 (3,1) + rdat%fq0 (2) * work (2)
      rdat%r00 (4,1) = rdat%r00 (4,1) + rdat%fq0 (2) * work (4)

      do i = 1, 3
      rdat%r01 (1, i) = rdat%r01 (1, i) - rdat%fq1 (1, i) * rdat%aqx * work (i)
      rdat%r01 (2, i) = rdat%r01 (2, i) - rdat%fq1 (1, i) * rdat%acy * work (i)
      rdat%r01 (3, i) = rdat%r01 (3, i) + rdat%fq1 (2, i) * work (i)
      enddo
      rdat%r01 (1, 4) = rdat%r01 (1, 4) - rdat%fq1 (1, 2) * rdat%aqx * work (4)
      rdat%r01 (2, 4) = rdat%r01 (2, 4) - rdat%fq1 (1, 2) * rdat%acy * work (4)
      rdat%r01 (3, 4) = rdat%r01 (3, 4) + rdat%fq1 (2, 2) * work (4)

      rdat%r02 (1,1) = rdat%r02 (1,1) + (rdat%fq2 (1,1) * rdat%aqx2 - rdat%fq1 (1, 2) ) * work (1)
      rdat%r02 (2,1) = rdat%r02 (2,1) + (rdat%fq2 (1,1) * rdat%acy2 - rdat%fq1 (1, 2) ) * work (1)
      rdat%r02 (3,1) = rdat%r02 (3,1) + (rdat%fq2 (3,1) - rdat%fq1 (1, 2) ) * work (1)
      rdat%r02 (4,1) = rdat%r02 (4,1) + rdat%fq2 (1,1) * rdat%aqxy * work (1)
      rdat%r02 (5,1) = rdat%r02 (5,1) - rdat%fq2 (2,1) * rdat%aqx * work (1)
      rdat%r02 (6,1) = rdat%r02 (6,1) - rdat%fq2 (2,1) * rdat%acy * work (1)

      end subroutine intk_04

! >
! >    @brief   ppps case
! >
! >    @details integration of a ppps case
! >
      subroutine intk_05 (rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd2, xmd3, xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: work, fwk
      integer :: i, j

! Generate jtype= 5 integrals

      dimension work (9), fwk (6, 3)

      if (ikl == 0) then
         do i = 1, 2
         rdat%r00 (1, i) = 0.0_dp
         rdat%r00 (2, i) = 0.0_dp
         rdat%r00 (3, i) = 0.0_dp
         rdat%r00 (4, i) = 0.0_dp
         rdat%r00 (5, i) = 0.0_dp
         enddo
         do i = 1, 13
         rdat%r01 (1, i) = 0.0_dp
         rdat%r01 (2, i) = 0.0_dp
         rdat%r01 (3, i) = 0.0_dp
         enddo
         do i = 1, 6
         rdat%r02 (1, i) = 0.0_dp
         rdat%r02 (2, i) = 0.0_dp
         rdat%r02 (3, i) = 0.0_dp
         rdat%r02 (4, i) = 0.0_dp
         rdat%r02 (5, i) = 0.0_dp
         rdat%r02 (6, i) = 0.0_dp
         enddo
         do j = 1, 10
         rdat%r03 (j, 1) = 0.0_dp
         enddo

         return
      endif

      xmd2 = rdat%x43 * 0.5d+00
      xmd3 = xmd2
      xmdt = xmd3 * xmd2

      xmdty = - xmdt * rdat%acy
      xmdtx = - xmdt * rdat%aqx
      xmdtxy = xmdt * rdat%aqxy

      work (1) = 1.0d0
      work (2) = - rdat%y03
      work (3) = rdat%y04
      work (4) = - rdat%y03 * rdat%y04
      work (5) = xmd3
      work (6) = xmd2
      work (7) = xmd2
      work (8) = - xmd3 * rdat%y03
      work (9) = xmd3 * rdat%y04

      do i = 1, 2
        rdat%r00 (1:5, i) = rdat%r00 (1:5, i) + rdat%fq0 (i) * work (1:5)
      end do

      do i = 1, 3
      fwk (1, i) = - rdat%fq1 (1, i) * rdat%aqx
      fwk (2, i) = - rdat%fq1 (1, i) * rdat%acy
      fwk (3, i) = rdat%fq1 (2, i)
      enddo
      do i = 1, 4
      rdat%r01 (1, i) = rdat%r01 (1, i) + fwk (1, 1) * work (i + 5)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fwk (2, 1) * work (i + 5)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fwk (3, 1) * work (i + 5)
      enddo
      do i = 5, 8
      rdat%r01 (1, i) = rdat%r01 (1, i) + fwk (1, 2) * work (i + 1)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fwk (2, 2) * work (i + 1)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fwk (3, 2) * work (i + 1)
      enddo
      do i = 9, 13
      rdat%r01 (1, i) = rdat%r01 (1, i) + fwk (1, 3) * work (i - 8)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fwk (2, 3) * work (i - 8)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fwk (3, 3) * work (i - 8)
      enddo

      do i = 1, 2
      rdat%r02 (1, i) = rdat%r02 (1, i) + (rdat%fq2 (1, i) * rdat%aqx2 - rdat%fq1 (1, i) )       &
      * xmdt
      rdat%r02 (2, i) = rdat%r02 (2, i) + (rdat%fq2 (1, i) * rdat%acy2 - rdat%fq1 (1, i) )       &
      * xmdt
      rdat%r02 (3, i) = rdat%r02 (3, i) + (rdat%fq2 (3, i) - rdat%fq1 (1, i) ) * xmdt
      rdat%r02 (4, i) = rdat%r02 (4, i) + rdat%fq2 (1, i) * xmdtxy
      rdat%r02 (5, i) = rdat%r02 (5, i) + rdat%fq2 (2, i) * xmdtx
      rdat%r02 (6, i) = rdat%r02 (6, i) + rdat%fq2 (2, i) * xmdty
      enddo
      fwk (1, 3) = rdat%fq2 (1, 3) * rdat%aqx2 - rdat%fq1 (1, 3)
      fwk (2, 3) = rdat%fq2 (1, 3) * rdat%acy2 - rdat%fq1 (1, 3)
      fwk (3, 3) = rdat%fq2 (3, 3) - rdat%fq1 (1, 3)
      fwk (4, 3) = rdat%fq2 (1, 3) * rdat%aqxy
      fwk (5, 3) = - rdat%fq2 (2, 3) * rdat%aqx
      fwk (6, 3) = - rdat%fq2 (2, 3) * rdat%acy
      do i = 3, 6
      do j = 1, 6
      rdat%r02 (j, i) = rdat%r02 (j, i) + fwk (j, 3) * work (i + 3)
      enddo
      enddo

      rdat%r03 (1, 1) = rdat%r03 (1, 1) + (rdat%fq3 (1, 1) * rdat%aqx2 - rdat%fq2 (1, 3) *  3 )  &
      * xmdtx
      rdat%r03 (2, 1) = rdat%r03 (2, 1) + (rdat%fq3 (1, 1) * rdat%aqx2 - rdat%fq2 (1, 3) )       &
      * xmdty
      rdat%r03 (3, 1) = rdat%r03 (3, 1) + (rdat%fq3 (2, 1) * rdat%aqx2 - rdat%fq2 (2, 3) )       &
      * xmdt
      rdat%r03 (4, 1) = rdat%r03 (4, 1) + (rdat%fq3 (1, 1) * rdat%acy2 - rdat%fq2 (1, 3) )       &
      * xmdtx
      rdat%r03 (5, 1) = rdat%r03 (5, 1) + rdat%fq3 (2, 1) * xmdtxy
      rdat%r03 (6, 1) = rdat%r03 (6, 1) + (rdat%fq3 (3, 1) - rdat%fq2 (1, 3) ) * xmdtx
      rdat%r03 (7, 1) = rdat%r03 (7, 1) + (rdat%fq3 (1, 1) * rdat%acy2 - rdat%fq2 (1, 3) *  3 )  &
      * xmdty
      rdat%r03 (8, 1) = rdat%r03 (8, 1) + (rdat%fq3 (2, 1) * rdat%acy2 - rdat%fq2 (2, 3) )       &
      * xmdt
      rdat%r03 (9, 1) = rdat%r03 (9, 1) + (rdat%fq3 (3, 1) - rdat%fq2 (1, 3) ) * xmdty
      rdat%r03 (10, 1) = rdat%r03 (10, 1) + (rdat%fq3 (4, 1) - rdat%fq2 (2, 3) *  3 )       &
      * xmdt

      end subroutine intk_05

! >
! >    @brief   pppp case
! >
! >    @details integration of a pppp case
! >
      subroutine intk_06 (rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd2, xmd3, xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: tmp0, tmp1, tmp2, tmp3, aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: work, fwk, fqw
      integer :: i, j

! Generate jtype= 6 integrals

      dimension work (9), fwk (10), fqw (6, 9)

      if (ikl == 0) then
         do i = 1, 5
         rdat%r00 (1, i) = 0.0_dp
         rdat%r00 (2, i) = 0.0_dp
         rdat%r00 (3, i) = 0.0_dp
         rdat%r00 (4, i) = 0.0_dp
         rdat%r00 (5, i) = 0.0_dp
         enddo
         do i = 1, 40
         rdat%r01 (1, i) = 0.0_dp
         rdat%r01 (2, i) = 0.0_dp
         rdat%r01 (3, i) = 0.0_dp
         enddo
         do i = 1, 26
         rdat%r02 (1, i) = 0.0_dp
         rdat%r02 (2, i) = 0.0_dp
         rdat%r02 (3, i) = 0.0_dp
         rdat%r02 (4, i) = 0.0_dp
         rdat%r02 (5, i) = 0.0_dp
         rdat%r02 (6, i) = 0.0_dp
         enddo
         do i = 1, 8
         do j = 1, 10
         rdat%r03 (j, i) = 0.0_dp
         enddo
         enddo
         do j = 1, 15
         rdat%r04 (j,1) = 0.0_dp
         enddo

         return
      endif

      xmd2 = rdat%x43 * 0.5d+00
      xmd3 = xmd2
      xmdt = xmd3 * xmd2

      xmdty = - xmdt * rdat%acy
      xmdtx = - xmdt * rdat%aqx
      xmdtxy = xmdt * rdat%aqxy

      work (1) = 1.0d0
      work (2) = - rdat%y03
      work (3) = rdat%y04
      work (4) = - rdat%y03 * rdat%y04
      work (5) = xmd3
      work (6) = xmd2
      work (7) = xmd2
      work (8) = - xmd3 * rdat%y03
      work (9) = xmd3 * rdat%y04

      do j = 1, 5
         do i = 1, 5
            rdat%r00 (i, j) = rdat%r00 (i, j) + rdat%fq0 (j) * work (i)
         end do
      end do

      do i = 1, 8
      fqw (1, i) = - rdat%fq1 (1, i) * rdat%aqx
      fqw (2, i) = - rdat%fq1 (1, i) * rdat%acy
      fqw (3, i) = rdat%fq1 (2, i)
      enddo
      tmp1 = rdat%fq1 (1, 5) * rdat%rab + rdat%fq1 (1, 8)
      tmp2 = rdat%fq1 (2, 5) * rdat%rab + rdat%fq1 (2, 8)
      fqw (1, 9) = - tmp1 * rdat%aqx
      fqw (2, 9) = - tmp1 * rdat%acy
      fqw (3, 9) = tmp2
      do i = 1, 4
      rdat%r01 (1, i) = rdat%r01 (1, i) + fqw (1, 1) * work (i + 5)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fqw (2, 1) * work (i + 5)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fqw (3, 1) * work (i + 5)
      enddo
      do i = 5, 8
      rdat%r01 (1, i) = rdat%r01 (1, i) + fqw (1, 2) * work (i + 1)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fqw (2, 2) * work (i + 1)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fqw (3, 2) * work (i + 1)
      enddo
      do i = 9, 12
      rdat%r01 (1, i) = rdat%r01 (1, i) + fqw (1, 3) * work (i - 3)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fqw (2, 3) * work (i - 3)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fqw (3, 3) * work (i - 3)
      enddo
      do i = 13, 16
      rdat%r01 (1, i) = rdat%r01 (1, i) + fqw (1, 4) * work (i - 7)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fqw (2, 4) * work (i - 7)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fqw (3, 4) * work (i - 7)
      enddo
      do i = 17, 20
      rdat%r01 (1, i) = rdat%r01 (1, i) + fqw (1, 5) * work (i - 11)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fqw (2, 5) * work (i - 11)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fqw (3, 5) * work (i - 11)
      enddo
      do i = 21, 25
      rdat%r01 (1, i) = rdat%r01 (1, i) + fqw (1, 6) * work (i - 20)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fqw (2, 6) * work (i - 20)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fqw (3, 6) * work (i - 20)
      enddo
      do i = 26, 30
      rdat%r01 (1, i) = rdat%r01 (1, i) + fqw (1, 7) * work (i - 25)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fqw (2, 7) * work (i - 25)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fqw (3, 7) * work (i - 25)
      enddo
      do i = 31, 35
      rdat%r01 (1, i) = rdat%r01 (1, i) + fqw (1, 8) * work (i - 30)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fqw (2, 8) * work (i - 30)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fqw (3, 8) * work (i - 30)
      enddo
      do i = 36, 40
      rdat%r01 (1, i) = rdat%r01 (1, i) + fqw (1, 9) * work (i - 35)
      rdat%r01 (2, i) = rdat%r01 (2, i) + fqw (2, 9) * work (i - 35)
      rdat%r01 (3, i) = rdat%r01 (3, i) + fqw (3, 9) * work (i - 35)
      enddo

      do i = 1, 5
      rdat%r02 (1, i) = rdat%r02 (1, i) + (rdat%fq2 (1, i) * rdat%aqx2 - rdat%fq1 (1, i) )       &
      * xmdt
      rdat%r02 (2, i) = rdat%r02 (2, i) + (rdat%fq2 (1, i) * rdat%acy2 - rdat%fq1 (1, i) )       &
      * xmdt
      rdat%r02 (3, i) = rdat%r02 (3, i) + (rdat%fq2 (3, i) - rdat%fq1 (1, i) ) * xmdt
      rdat%r02 (4, i) = rdat%r02 (4, i) + rdat%fq2 (1, i) * xmdtxy
      rdat%r02 (5, i) = rdat%r02 (5, i) + rdat%fq2 (2, i) * xmdtx
      rdat%r02 (6, i) = rdat%r02 (6, i) + rdat%fq2 (2, i) * xmdty
      enddo
      tmp0 = rdat%fq1 (1, 5) * rdat%rab + rdat%fq1 (1, 8)
      tmp1 = rdat%fq2 (1, 5) * rdat%rab + rdat%fq2 (1, 8)
      tmp2 = rdat%fq2 (2, 5) * rdat%rab + rdat%fq2 (2, 8)
      tmp3 = rdat%fq2 (3, 5) * rdat%rab + rdat%fq2 (3, 8)
      fqw (1, 5) = tmp1 * rdat%aqx2 - tmp0
      fqw (2, 5) = tmp1 * rdat%acy2 - tmp0
      fqw (3, 5) = tmp3 - tmp0
      fqw (4, 5) = tmp1 * rdat%aqxy
      fqw (5, 5) = - tmp2 * rdat%aqx
      fqw (6, 5) = - tmp2 * rdat%acy
      do i = 6, 9
      fqw (1, i) = rdat%fq2 (1, i) * rdat%aqx2 - rdat%fq1 (1, i)
      fqw (2, i) = rdat%fq2 (1, i) * rdat%acy2 - rdat%fq1 (1, i)
      fqw (3, i) = rdat%fq2 (3, i) - rdat%fq1 (1, i)
      fqw (4, i) = rdat%fq2 (1, i) * rdat%aqxy
      fqw (5, i) = - rdat%fq2 (2, i) * rdat%aqx
      fqw (6, i) = - rdat%fq2 (2, i) * rdat%acy
      enddo
      do i = 6, 9
         do j = 1, 6
            rdat%r02 (j, i) = rdat%r02 (j, i) + fqw (j, 6) * work (i)
         end do
      end do
      do i = 10, 13
         do j = 1, 6
            rdat%r02 (j, i) = rdat%r02 (j, i) + fqw (j, 7) * work (i - 4)
         end do
      end do
      do i = 14, 17
         do j = 1, 6
            rdat%r02 (j, i) = rdat%r02 (j, i) + fqw (j, 8) * work (i - 8)
         end do
      end do
      do i = 18, 21
         do j = 1, 6
            rdat%r02 (j, i) = rdat%r02 (j, i) + fqw (j, 5) * work (i - 12)
         end do
      end do
      do i = 22, 26
         do j = 1, 6
            rdat%r02 (j, i) = rdat%r02 (j, i) + fqw (j, 9) * work (i - 21)
         end do
      end do

      rdat%fq2 (1, 5) = rdat%fq2 (1, 5) * rdat%rab + rdat%fq2 (1, 8)
      rdat%fq2 (2, 5) = rdat%fq2 (2, 5) * rdat%rab + rdat%fq2 (2, 8)
      rdat%fq3 (1, 1) = rdat%fq3 (1, 1) * rdat%rab + rdat%fq3 (1, 4)
      rdat%fq3 (2, 1) = rdat%fq3 (2, 1) * rdat%rab + rdat%fq3 (2, 4)
      rdat%fq3 (3, 1) = rdat%fq3 (3, 1) * rdat%rab + rdat%fq3 (3, 4)
      rdat%fq3 (4, 1) = rdat%fq3 (4, 1) * rdat%rab + rdat%fq3 (4, 4)
      do i = 1, 4
      rdat%r03 (1, i) = rdat%r03 (1, i) + (rdat%fq3 (1, i) * rdat%aqx2 - rdat%fq2 (1, i + 4)     &
      *  3 ) * xmdtx
      rdat%r03 (2, i) = rdat%r03 (2, i) + (rdat%fq3 (1, i) * rdat%aqx2 - rdat%fq2 (1, i + 4) )   &
      * xmdty
      rdat%r03 (3, i) = rdat%r03 (3, i) + (rdat%fq3 (2, i) * rdat%aqx2 - rdat%fq2 (2, i + 4) )   &
      * xmdt
      rdat%r03 (4, i) = rdat%r03 (4, i) + (rdat%fq3 (1, i) * rdat%acy2 - rdat%fq2 (1, i + 4) )   &
      * xmdtx
      rdat%r03 (5, i) = rdat%r03 (5, i) + rdat%fq3 (2, i) * xmdtxy
      rdat%r03 (6, i) = rdat%r03 (6, i) + (rdat%fq3 (3, i) - rdat%fq2 (1, i + 4) ) * xmdtx
      rdat%r03 (7, i) = rdat%r03 (7, i) + (rdat%fq3 (1, i) * rdat%acy2 - rdat%fq2 (1, i + 4)     &
      *  3 ) * xmdty
      rdat%r03 (8, i) = rdat%r03 (8, i) + (rdat%fq3 (2, i) * rdat%acy2 - rdat%fq2 (2, i + 4) )   &
      * xmdt
      rdat%r03 (9, i) = rdat%r03 (9, i) + (rdat%fq3 (3, i) - rdat%fq2 (1, i + 4) ) * xmdty
      rdat%r03 (10, i) = rdat%r03 (10, i) + (rdat%fq3 (4, i) - rdat%fq2 (2, i + 4) *  3 )   &
      * xmdt
      enddo
      fwk (1) = - (rdat%fq3 (1, 5) * rdat%aqx2 - rdat%fq2 (1, 9) *  3 ) * rdat%aqx
      fwk (2) = - (rdat%fq3 (1, 5) * rdat%aqx2 - rdat%fq2 (1, 9) ) * rdat%acy
      fwk (3) = rdat%fq3 (2, 5) * rdat%aqx2 - rdat%fq2 (2, 9)
      fwk (4) = - (rdat%fq3 (1, 5) * rdat%acy2 - rdat%fq2 (1, 9) ) * rdat%aqx
      fwk (5) = rdat%fq3 (2, 5) * rdat%aqxy
      fwk (6) = - (rdat%fq3 (3, 5) - rdat%fq2 (1, 9) ) * rdat%aqx
      fwk (7) = - (rdat%fq3 (1, 5) * rdat%acy2 - rdat%fq2 (1, 9) *  3 ) * rdat%acy
      fwk (8) = rdat%fq3 (2, 5) * rdat%acy2 - rdat%fq2 (2, 9)
      fwk (9) = - (rdat%fq3 (3, 5) - rdat%fq2 (1, 9) ) * rdat%acy
      fwk (10) = rdat%fq3 (4, 5) - rdat%fq2 (2, 9) *  3
      do i = 5, 8
         do j = 1, 10
            rdat%r03 (j, i) = rdat%r03 (j, i) + fwk (j) * work (i + 1)
         end do
      end do

      aqx4 = rdat%aqx2 * rdat%aqx2
      acy4 = rdat%acy2 * rdat%acy2
      x2y2 = rdat%aqx2 * rdat%acy2
      q2c2 = rdat%aqx2 + rdat%acy2
      rdat%r04 (1,1) = rdat%r04 (1,1) + (rdat%fq4 (1,1) * aqx4 - rdat%fq3 (1, 5) *  6  * rdat%aqx2 +   &
      rdat%fq2 (1, 9) *  3 ) * xmdt
      rdat%r04 (2,1) = rdat%r04 (2,1) + (rdat%fq4 (1,1) * rdat%aqx2 - rdat%fq3 (1, 5) *  3 ) * xmdtxy
      rdat%r04 (3,1) = rdat%r04 (3,1) + (rdat%fq4 (2,1) * rdat%aqx2 - rdat%fq3 (2, 5) *  3 ) * xmdtx
      rdat%r04 (4,1) = rdat%r04 (4,1) + (rdat%fq4 (1,1) * x2y2 - rdat%fq3 (1, 5) * q2c2 + rdat%fq2 (1, &
      9) ) * xmdt
      rdat%r04 (5,1) = rdat%r04 (5,1) + (rdat%fq4 (2,1) * rdat%aqx2 - rdat%fq3 (2, 5) ) * xmdty
      rdat%r04 (6,1) = rdat%r04 (6,1) + (rdat%fq4 (3,1) * rdat%aqx2 - rdat%fq3 (1, 5) * rdat%aqx2 - rdat%fq3 (3, &
      5) + rdat%fq2 (1, 9) ) * xmdt
      rdat%r04 (7,1) = rdat%r04 (7,1) + (rdat%fq4 (1,1) * rdat%acy2 - rdat%fq3 (1, 5) *  3 ) * xmdtxy
      rdat%r04 (8,1) = rdat%r04 (8,1) + (rdat%fq4 (2,1) * rdat%acy2 - rdat%fq3 (2, 5) ) * xmdtx
      rdat%r04 (9,1) = rdat%r04 (9,1) + (rdat%fq4 (3,1) - rdat%fq3 (1, 5) ) * xmdtxy
      rdat%r04 (10,1) = rdat%r04 (10,1) + (rdat%fq4 (4,1) - rdat%fq3 (2, 5) *  3 ) * xmdtx
      rdat%r04 (11,1) = rdat%r04 (11,1) + (rdat%fq4 (1,1) * acy4 - rdat%fq3 (1, 5) *  6  * rdat%acy2 + &
      rdat%fq2 (1, 9) *  3 ) * xmdt
      rdat%r04 (12,1) = rdat%r04 (12,1) + (rdat%fq4 (2,1) * rdat%acy2 - rdat%fq3 (2, 5) *  3 ) * xmdty
      rdat%r04 (13,1) = rdat%r04 (13,1) + (rdat%fq4 (3,1) * rdat%acy2 - rdat%fq3 (1, 5) * rdat%acy2 - rdat%fq3 ( &
      3, 5) + rdat%fq2 (1, 9) ) * xmdt
      rdat%r04 (14,1) = rdat%r04 (14,1) + (rdat%fq4 (4,1) - rdat%fq3 (2, 5) *  3 ) * xmdty
      rdat%r04 (15,1) = rdat%r04 (15,1) + (rdat%fq4 (5,1) - rdat%fq3 (3, 5) *  6  + rdat%fq2 (1, 9)    &
      *  3 ) * xmdt

      end subroutine intk_06

! >
! >    @brief   dsss case
! >
! >    @details integration of the dsss case
! >
      subroutine intk_07(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd1, xmd2, xmdt

! Generate jtype= 7 integrals

      if(ikl == 0) then
         rdat%r00(1:2,1)= 0.0_dp
         rdat%r01(:,1)= 0.0_dp
         rdat%r02(:,1)= 0.0_dp
         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd1= xmd2
      xmdt= xmd2*xmd1
      xmd2=-xmd1*rdat%y03

      rdat%r00(1,1)= rdat%r00(1,1)+rdat%fq0(1)*xmd1
      rdat%r00(2,1)= rdat%r00(2,1)+rdat%fq0(1)*rdat%y03 *rdat%y03

      rdat%r01(1,1)= rdat%r01(1,1)-rdat%fq1(1,1)*rdat%aqx *xmd2
      rdat%r01(2,1)= rdat%r01(2,1)-rdat%fq1(1,1)*rdat%acy *xmd2
      rdat%r01(3,1)= rdat%r01(3,1)+rdat%fq1(2,1)     *xmd2

      rdat%r02(1,1)= rdat%r02(1,1)+(rdat%fq2(1,1)*rdat%aqx2-rdat%fq1(1,1))*xmdt
      rdat%r02(2,1)= rdat%r02(2,1)+(rdat%fq2(1,1)*rdat%acy2-rdat%fq1(1,1))*xmdt
      rdat%r02(3,1)= rdat%r02(3,1)+(rdat%fq2(3,1)     -rdat%fq1(1,1))*xmdt
      rdat%r02(4,1)= rdat%r02(4,1)+ rdat%fq2(1,1)*rdat%aqxy         *xmdt
      rdat%r02(5,1)= rdat%r02(5,1)- rdat%fq2(2,1)*rdat%aqx          *xmdt
      rdat%r02(6,1)= rdat%r02(6,1)- rdat%fq2(2,1)*rdat%acy          *xmdt

      end subroutine intk_07

! >
! >    @brief   dpss case
! >
! >    @details integration of the dpss case
! >
      subroutine intk_08(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd1, xmd2, xmd3, xmd4, xmd5, xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: y33, y334, y34
      real(kind=dp) :: work, fwk
      integer :: i, j

! Generate jtype= 8 integrals

      dimension  work(4),fwk(6)

      if(ikl == 0) then
         rdat%r00(1,1)= 0.0_dp
         rdat%r00(2,1)= 0.0_dp
         rdat%r00(3,1)= 0.0_dp
         rdat%r00(4,1)= 0.0_dp
         rdat%r00(5,1)= 0.0_dp
         do i= 1, 4
            rdat%r01(1,i)= 0.0_dp
            rdat%r01(2,i)= 0.0_dp
            rdat%r01(3,i)= 0.0_dp
         enddo
         do i= 1, 3
            rdat%r02(1,i)= 0.0_dp
            rdat%r02(2,i)= 0.0_dp
            rdat%r02(3,i)= 0.0_dp
            rdat%r02(4,i)= 0.0_dp
            rdat%r02(5,i)= 0.0_dp
            rdat%r02(6,i)= 0.0_dp
         enddo
         do j= 1,10
            rdat%r03(j,1)= 0.0_dp
         enddo

         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd1= xmd2
      xmd3= xmd2
      xmd4= xmd2*xmd2
      xmd5= xmd4
      xmdt= xmd4*xmd3

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      y33 = rdat%y03 *rdat%y03
      y34 =-rdat%y03 *rdat%y04
      y334= y33 *rdat%y04
      work(1)=-xmd1*rdat%y03
      work(2)= xmd5
      work(3)= xmd3*y33
      work(4)= xmd3*y34

      rdat%r00(1,1)= rdat%r00(1,1)+rdat%fq0(1)*xmd1
      rdat%r00(2,1)= rdat%r00(2,1)+rdat%fq0(1)*y33
      rdat%r00(3,1)= rdat%r00(3,1)-rdat%fq0(1)*xmd3*rdat%y03
      rdat%r00(4,1)= rdat%r00(4,1)+rdat%fq0(1)*xmd3*rdat%y04
      rdat%r00(5,1)= rdat%r00(5,1)+rdat%fq0(1)*y334

         fwk(1)=-rdat%fq1(1,1)*rdat%aqx
         fwk(2)=-rdat%fq1(1,1)*rdat%acy
         fwk(3)= rdat%fq1(2,1)
      do i= 1, 4
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1)*work(i)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2)*work(i)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3)*work(i)
      enddo

      work(1)= xmd4
      work(2)=-xmd5*rdat%y03
      work(3)= xmd5*rdat%y04
         fwk(1)= rdat%fq2(1,1)*rdat%aqx2-rdat%fq1(1,1)
         fwk(2)= rdat%fq2(1,1)*rdat%acy2-rdat%fq1(1,1)
         fwk(3)= rdat%fq2(3,1)     -rdat%fq1(1,1)
         fwk(4)= rdat%fq2(1,1)*rdat%aqxy
         fwk(5)=-rdat%fq2(2,1)*rdat%aqx
         fwk(6)=-rdat%fq2(2,1)*rdat%acy
      do i=1,3
         do j=1,6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j)*work(i)
         end do
      end do

      rdat%r03( 1,1)= rdat%r03( 1,1)+(rdat%fq3(1,1)*rdat%aqx2-rdat%fq2(1,1)* 3 )*xmdtx
      rdat%r03( 2,1)= rdat%r03( 2,1)+(rdat%fq3(1,1)*rdat%aqx2-rdat%fq2(1,1)    )*xmdty
      rdat%r03( 3,1)= rdat%r03( 3,1)+(rdat%fq3(2,1)*rdat%aqx2-rdat%fq2(2,1)    )*xmdt
      rdat%r03( 4,1)= rdat%r03( 4,1)+(rdat%fq3(1,1)*rdat%acy2-rdat%fq2(1,1)    )*xmdtx
      rdat%r03( 5,1)= rdat%r03( 5,1)+ rdat%fq3(2,1)                    *xmdtxy
      rdat%r03( 6,1)= rdat%r03( 6,1)+(rdat%fq3(3,1)     -rdat%fq2(1,1)    )*xmdtx
      rdat%r03( 7,1)= rdat%r03( 7,1)+(rdat%fq3(1,1)*rdat%acy2-rdat%fq2(1,1)* 3 )*xmdty
      rdat%r03( 8,1)= rdat%r03( 8,1)+(rdat%fq3(2,1)*rdat%acy2-rdat%fq2(2,1)    )*xmdt
      rdat%r03( 9,1)= rdat%r03( 9,1)+(rdat%fq3(3,1)     -rdat%fq2(1,1)    )*xmdty
      rdat%r03(10,1)= rdat%r03(10,1)+(rdat%fq3(4,1)     -rdat%fq2(2,1)* 3 )*xmdt

      end subroutine intk_08

! >
! >    @brief   dsps case
! >
! >    @details integration of the dsps case
! >
      subroutine intk_09(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd2, xmd3, xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: work
      integer :: i

! Generate jtype= 9 integrals

      dimension  work(3)

      if(ikl == 0) then
         rdat%r00(1:4, 1)= 0.0_dp
         rdat%r01(:,1:4)= 0.0_dp
         rdat%r02(:,1:3)= 0.0_dp
         rdat%r03(:,1)= 0.0_dp
         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd3= xmd2
      xmdt= xmd3*xmd2

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      work(1)= xmd3
      work(2)= rdat%y03 *rdat%y03
      work(3)=-xmd3*rdat%y03

      rdat%r00(1,1)= rdat%r00(1,1)+rdat%fq0(1)*work(1)
      rdat%r00(2,1)= rdat%r00(2,1)+rdat%fq0(1)*work(2)
      rdat%r00(3,1)= rdat%r00(3,1)+rdat%fq0(2)*work(1)
      rdat%r00(4,1)= rdat%r00(4,1)+rdat%fq0(2)*work(2)

      do i= 1, 2
         rdat%r01(1,i)= rdat%r01(1,i)-rdat%fq1(1,i)*rdat%aqx *work(3)
         rdat%r01(2,i)= rdat%r01(2,i)-rdat%fq1(1,i)*rdat%acy *work(3)
         rdat%r01(3,i)= rdat%r01(3,i)+rdat%fq1(2,i)     *work(3)
      enddo
      do i= 3, 4
         rdat%r01(1,i)= rdat%r01(1,i)-rdat%fq1(1,3)*rdat%aqx *work(i- 2)
         rdat%r01(2,i)= rdat%r01(2,i)-rdat%fq1(1,3)*rdat%acy *work(i- 2)
         rdat%r01(3,i)= rdat%r01(3,i)+rdat%fq1(2,3)     *work(i- 2)
      enddo

      work(1)= xmdt
      work(2)= xmdt
      do i= 1, 3
         rdat%r02(1,i)= rdat%r02(1,i)+(rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i))*work(i)
         rdat%r02(2,i)= rdat%r02(2,i)+(rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i))*work(i)
         rdat%r02(3,i)= rdat%r02(3,i)+(rdat%fq2(3,i)     -rdat%fq1(1,i))*work(i)
         rdat%r02(4,i)= rdat%r02(4,i)+ rdat%fq2(1,i)*rdat%aqxy           *work(i)
         rdat%r02(5,i)= rdat%r02(5,i)- rdat%fq2(2,i)*rdat%aqx            *work(i)
         rdat%r02(6,i)= rdat%r02(6,i)- rdat%fq2(2,i)*rdat%acy            *work(i)
      enddo

      rdat%r03( 1,1)= rdat%r03( 1,1)+(rdat%fq3(1,1)*rdat%aqx2-rdat%fq2(1,3)* 3 )*xmdtx
      rdat%r03( 2,1)= rdat%r03( 2,1)+(rdat%fq3(1,1)*rdat%aqx2-rdat%fq2(1,3)    )*xmdty
      rdat%r03( 3,1)= rdat%r03( 3,1)+(rdat%fq3(2,1)*rdat%aqx2-rdat%fq2(2,3)    )*xmdt
      rdat%r03( 4,1)= rdat%r03( 4,1)+(rdat%fq3(1,1)*rdat%acy2-rdat%fq2(1,3)    )*xmdtx
      rdat%r03( 5,1)= rdat%r03( 5,1)+ rdat%fq3(2,1)                    *xmdtxy
      rdat%r03( 6,1)= rdat%r03( 6,1)+(rdat%fq3(3,1)     -rdat%fq2(1,3)    )*xmdtx
      rdat%r03( 7,1)= rdat%r03( 7,1)+(rdat%fq3(1,1)*rdat%acy2-rdat%fq2(1,3)* 3 )*xmdty
      rdat%r03( 8,1)= rdat%r03( 8,1)+(rdat%fq3(2,1)*rdat%acy2-rdat%fq2(2,3)    )*xmdt
      rdat%r03( 9,1)= rdat%r03( 9,1)+(rdat%fq3(3,1)     -rdat%fq2(1,3)    )*xmdty
      rdat%r03(10,1)= rdat%r03(10,1)+(rdat%fq3(4,1)     -rdat%fq2(2,3)* 3 )*xmdt

      end subroutine intk_09

! >
! >    @brief   ddss case
! >
! >    @details integration of the ddss case
! >
      subroutine intk_10(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd2, xmd3, xmd4, xmd6
      real(kind=dp) :: xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: y33, y334, y34, y344, y44
      real(kind=dp) :: work, fwk
      integer :: i, j

! Generate jtype=10 integrals

      dimension  work(4),fwk(10)

      if(ikl == 0) then
         rdat%r00(:,1)= 0.0_dp
         rdat%r01(:,1:4)= 0.0_dp
         rdat%r02(:,1:4)= 0.0_dp
         rdat%r03(:,1:2)= 0.0_dp
         rdat%r04(:,1)= 0.0_dp

         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd3= xmd2
      xmd4= xmd3*xmd2
      xmd6= xmd4*xmd2
      xmdt= xmd6*xmd2
      xmd2= xmd3

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      y33 = rdat%y03 *rdat%y03
      y34 =-rdat%y03 *rdat%y04
      y44 = rdat%y04 *rdat%y04
      y334= y33 *rdat%y04
      y344= y34 *rdat%y04
      work(1)=-xmd4*rdat%y03
      work(2)= xmd4*rdat%y04
      work(3)= xmd2*y334
      work(4)= xmd2*y344

      rdat%r00(1,1)= rdat%r00(1,1)+rdat%fq0(1)*xmd4
      rdat%r00(2,1)= rdat%r00(2,1)+rdat%fq0(1)*xmd2*y33
      rdat%r00(3,1)= rdat%r00(3,1)+rdat%fq0(1)*xmd2*y34
      rdat%r00(4,1)= rdat%r00(4,1)+rdat%fq0(1)*xmd2*y44
      rdat%r00(5,1)= rdat%r00(5,1)+rdat%fq0(1)*y33 *y44

         fwk(1)=-rdat%fq1(1,1)*rdat%aqx
         fwk(2)=-rdat%fq1(1,1)*rdat%acy
         fwk(3)= rdat%fq1(2,1)
      do i= 1, 4
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1)*work(i)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2)*work(i)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3)*work(i)
      enddo

      work(1)= xmd6
      work(2)= xmd4*y33
      work(3)= xmd4*y34
      work(4)= xmd4*y44
         fwk(1)= rdat%fq2(1,1)*rdat%aqx2-rdat%fq1(1,1)
         fwk(2)= rdat%fq2(1,1)*rdat%acy2-rdat%fq1(1,1)
         fwk(3)= rdat%fq2(3,1)     -rdat%fq1(1,1)
         fwk(4)= rdat%fq2(1,1)*rdat%aqxy
         fwk(5)=-rdat%fq2(2,1)*rdat%aqx
         fwk(6)=-rdat%fq2(2,1)*rdat%acy
      do i= 1, 4
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j)*work(i)
         end do
      end do

      work(1)=-xmd6*rdat%y03
      work(2)= xmd6*rdat%y04
         fwk( 1)=-(rdat%fq3(1,1)*rdat%aqx2-rdat%fq2(1,1)* 3 )*rdat%aqx
         fwk( 2)=-(rdat%fq3(1,1)*rdat%aqx2-rdat%fq2(1,1)    )*rdat%acy
         fwk( 3)=  rdat%fq3(2,1)*rdat%aqx2-rdat%fq2(2,1)
         fwk( 4)=-(rdat%fq3(1,1)*rdat%acy2-rdat%fq2(1,1)    )*rdat%aqx
         fwk( 5)=  rdat%fq3(2,1)                  *rdat%aqxy
         fwk( 6)=-(rdat%fq3(3,1)     -rdat%fq2(1,1)    )*rdat%aqx
         fwk( 7)=-(rdat%fq3(1,1)*rdat%acy2-rdat%fq2(1,1)* 3 )*rdat%acy
         fwk( 8)=  rdat%fq3(2,1)*rdat%acy2-rdat%fq2(2,1)
         fwk( 9)=-(rdat%fq3(3,1)     -rdat%fq2(1,1)    )*rdat%acy
         fwk(10)=  rdat%fq3(4,1)     -rdat%fq2(2,1)* 3
      do i= 1, 2
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j)*work(i)
         end do
      end do

      aqx4= rdat%aqx2*rdat%aqx2
      acy4= rdat%acy2*rdat%acy2
      x2y2= rdat%aqx2*rdat%acy2
      q2c2= rdat%aqx2+rdat%acy2
      rdat%r04( 1,1)= rdat%r04( 1,1)+(rdat%fq4(1,1)*aqx4-rdat%fq3(1,1)* 6 *rdat%aqx2                   &
     &                                    +rdat%fq2(1,1)* 3            )*xmdt
      rdat%r04( 2,1)= rdat%r04( 2,1)+(rdat%fq4(1,1)*rdat%aqx2-rdat%fq3(1,1)* 3            )*xmdtxy
      rdat%r04( 3,1)= rdat%r04( 3,1)+(rdat%fq4(2,1)*rdat%aqx2-rdat%fq3(2,1)* 3            )*xmdtx
      rdat%r04( 4,1)= rdat%r04( 4,1)+(rdat%fq4(1,1)*x2y2-rdat%fq3(1,1)*q2c2+rdat%fq2(1,1))*xmdt
      rdat%r04( 5,1)= rdat%r04( 5,1)+(rdat%fq4(2,1)*rdat%aqx2-rdat%fq3(2,1)               )*xmdty
      rdat%r04( 6,1)= rdat%r04( 6,1)+(rdat%fq4(3,1)*rdat%aqx2-rdat%fq3(1,1)*rdat%aqx2-rdat%fq3(3,1)               &
     &                                    +rdat%fq2(1,1)               )*xmdt
      rdat%r04( 7,1)= rdat%r04( 7,1)+(rdat%fq4(1,1)*rdat%acy2-rdat%fq3(1,1)* 3            )*xmdtxy
      rdat%r04( 8,1)= rdat%r04( 8,1)+(rdat%fq4(2,1)*rdat%acy2-rdat%fq3(2,1)               )*xmdtx
      rdat%r04( 9,1)= rdat%r04( 9,1)+(rdat%fq4(3,1)     -rdat%fq3(1,1)               )*xmdtxy
      rdat%r04(10,1)= rdat%r04(10,1)+(rdat%fq4(4,1)     -rdat%fq3(2,1)* 3            )*xmdtx
      rdat%r04(11,1)= rdat%r04(11,1)+(rdat%fq4(1,1)*acy4-rdat%fq3(1,1)* 6 *rdat%acy2                   &
     &                                    +rdat%fq2(1,1)* 3            )*xmdt
      rdat%r04(12,1)= rdat%r04(12,1)+(rdat%fq4(2,1)*rdat%acy2-rdat%fq3(2,1)* 3            )*xmdty
      rdat%r04(13,1)= rdat%r04(13,1)+(rdat%fq4(3,1)*rdat%acy2-rdat%fq3(1,1)*rdat%acy2-rdat%fq3(3,1)               &
     &                                    +rdat%fq2(1,1)               )*xmdt
      rdat%r04(14,1)= rdat%r04(14,1)+(rdat%fq4(4,1)     -rdat%fq3(2,1)* 3            )*xmdty
      rdat%r04(15,1)= rdat%r04(15,1)+(rdat%fq4(5,1)  -rdat%fq3(3,1)* 6 +rdat%fq2(1,1)* 3 )*xmdt

      end subroutine intk_10

! >
! >    @brief   dpps case
! >
! >    @details integration of the dpps case
! >
      subroutine intk_11(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd1, xmd2, xmd3, xmd4, xmd5
      real(kind=dp) :: xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: y33, y334, y34
      real(kind=dp) :: work, fwk
      integer :: i, j

! Generate jtype=11 integrals

      dimension  work(12),fwk(10,3)

      if(ikl == 0) then
         rdat%r00(:,1:2)= 0.0_dp
         rdat%r01(:,1:13)= 0.0_dp
         rdat%r02(:,1:10)= 0.0_dp
         rdat%r03(:,1:5)= 0.0_dp
         rdat%r04(:,1)= 0.0_dp
         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd1= xmd2
      xmd3= xmd2
      xmd4= xmd2*xmd2
      xmd5= xmd4
      xmdt= xmd4*xmd3

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      y33 = rdat%y03 *rdat%y03
      y34 =-rdat%y03 *rdat%y04
      y334= y33 *rdat%y04

      work( 1)= xmd1
      work( 2)= y33
      work( 3)=-xmd3*rdat%y03
      work( 4)= xmd3*rdat%y04
      work( 5)= y334
      work( 6)=-xmd1*rdat%y03
      work( 7)= xmd5
      work( 8)= xmd3*y33
      work( 9)= xmd2*y34
      work(10)= xmd4
      work(11)=-xmd5*rdat%y03
      work(12)= xmd5*rdat%y04

      do i= 1, 2
         do j= 1, 5
            rdat%r00(j,i)= rdat%r00(j,i)+rdat%fq0(i)*work(j)
         end do
      end do

      do i= 1, 3
         fwk(1,i)=-rdat%fq1(1,i)*rdat%aqx
         fwk(2,i)=-rdat%fq1(1,i)*rdat%acy
         fwk(3,i)= rdat%fq1(2,i)
      enddo
      do i= 1, 4
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,1)*work(i+ 5)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,1)*work(i+ 5)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,1)*work(i+ 5)
      enddo
      do i= 5, 8
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,2)*work(i+ 1)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,2)*work(i+ 1)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,2)*work(i+ 1)
      enddo
      do i= 9,13
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,3)*work(i- 8)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,3)*work(i- 8)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,3)*work(i- 8)
      enddo

      do i= 1, 3
         fwk(1,i)= rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i)
         fwk(2,i)= rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i)
         fwk(3,i)= rdat%fq2(3,i)     -rdat%fq1(1,i)
         fwk(4,i)= rdat%fq2(1,i)*rdat%aqxy
         fwk(5,i)=-rdat%fq2(2,i)*rdat%aqx
         fwk(6,i)=-rdat%fq2(2,i)*rdat%acy
      enddo
      do i= 1, 3
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,1)*work(i+ 9)
         enddo
      enddo
      do i= 4, 6
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,2)*work(i+ 6)
         enddo
      enddo
      do i= 7,10
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,3)*work(i- 1)
         enddo
      enddo

      do i= 1, 2
         rdat%r03( 1,i)= rdat%r03( 1,i)+(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)* 3 )*xmdtx
         rdat%r03( 2,i)= rdat%r03( 2,i)+(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)    )*xmdty
         rdat%r03( 3,i)= rdat%r03( 3,i)+(rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i)    )*xmdt
         rdat%r03( 4,i)= rdat%r03( 4,i)+(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)    )*xmdtx
         rdat%r03( 5,i)= rdat%r03( 5,i)+ rdat%fq3(2,i)                    *xmdtxy
         rdat%r03( 6,i)= rdat%r03( 6,i)+(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*xmdtx
         rdat%r03( 7,i)= rdat%r03( 7,i)+(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)* 3 )*xmdty
         rdat%r03( 8,i)= rdat%r03( 8,i)+(rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i)    )*xmdt
         rdat%r03( 9,i)= rdat%r03( 9,i)+(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*xmdty
         rdat%r03(10,i)= rdat%r03(10,i)+(rdat%fq3(4,i)     -rdat%fq2(2,i)* 3 )*xmdt
      enddo
         fwk( 1,3)=-(rdat%fq3(1,3)*rdat%aqx2-rdat%fq2(1,3)* 3 )*rdat%aqx
         fwk( 2,3)=-(rdat%fq3(1,3)*rdat%aqx2-rdat%fq2(1,3)    )*rdat%acy
         fwk( 3,3)=  rdat%fq3(2,3)*rdat%aqx2-rdat%fq2(2,3)
         fwk( 4,3)=-(rdat%fq3(1,3)*rdat%acy2-rdat%fq2(1,3)    )*rdat%aqx
         fwk( 5,3)=  rdat%fq3(2,3)                    *rdat%aqxy
         fwk( 6,3)=-(rdat%fq3(3,3)     -rdat%fq2(1,3)    )*rdat%aqx
         fwk( 7,3)=-(rdat%fq3(1,3)*rdat%acy2-rdat%fq2(1,3)* 3 )*rdat%acy
         fwk( 8,3)=  rdat%fq3(2,3)*rdat%acy2-rdat%fq2(2,3)
         fwk( 9,3)=-(rdat%fq3(3,3)     -rdat%fq2(1,3)    )*rdat%acy
         fwk(10,3)=  rdat%fq3(4,3)     -rdat%fq2(2,3)* 3
      do i= 3, 5
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,3)*work(i+ 7)
         end do
      end do

      aqx4= rdat%aqx2*rdat%aqx2
      acy4= rdat%acy2*rdat%acy2
      x2y2= rdat%aqx2*rdat%acy2
      q2c2= rdat%aqx2+rdat%acy2
      rdat%r04( 1,1)= rdat%r04( 1,1)+(rdat%fq4(1,1)*aqx4-rdat%fq3(1,3)* 6 *rdat%aqx2                 &
     &                                    +rdat%fq2(1,3)* 3            )*xmdt
      rdat%r04( 2,1)= rdat%r04( 2,1)+(rdat%fq4(1,1)*rdat%aqx2-rdat%fq3(1,3)* 3            )*xmdtxy
      rdat%r04( 3,1)= rdat%r04( 3,1)+(rdat%fq4(2,1)*rdat%aqx2-rdat%fq3(2,3)* 3            )*xmdtx
      rdat%r04( 4,1)= rdat%r04( 4,1)+(rdat%fq4(1,1)*x2y2-rdat%fq3(1,3)*q2c2+rdat%fq2(1,3))*xmdt
      rdat%r04( 5,1)= rdat%r04( 5,1)+(rdat%fq4(2,1)*rdat%aqx2-rdat%fq3(2,3)               )*xmdty
      rdat%r04( 6,1)= rdat%r04( 6,1)+(rdat%fq4(3,1)*rdat%aqx2-rdat%fq3(1,3)*rdat%aqx2-rdat%fq3(3,3)           &
     &                                    +rdat%fq2(1,3)               )*xmdt
      rdat%r04( 7,1)= rdat%r04( 7,1)+(rdat%fq4(1,1)*rdat%acy2-rdat%fq3(1,3)* 3            )*xmdtxy
      rdat%r04( 8,1)= rdat%r04( 8,1)+(rdat%fq4(2,1)*rdat%acy2-rdat%fq3(2,3)               )*xmdtx
      rdat%r04( 9,1)= rdat%r04( 9,1)+(rdat%fq4(3,1)     -rdat%fq3(1,3)               )*xmdtxy
      rdat%r04(10,1)= rdat%r04(10,1)+(rdat%fq4(4,1)     -rdat%fq3(2,3)* 3            )*xmdtx
      rdat%r04(11,1)= rdat%r04(11,1)+(rdat%fq4(1,1)*acy4-rdat%fq3(1,3)* 6 *rdat%acy2                 &
     &                                    +rdat%fq2(1,3)* 3            )*xmdt
      rdat%r04(12,1)= rdat%r04(12,1)+(rdat%fq4(2,1)*rdat%acy2-rdat%fq3(2,3)* 3            )*xmdty
      rdat%r04(13,1)= rdat%r04(13,1)+(rdat%fq4(3,1)*rdat%acy2-rdat%fq3(1,3)*rdat%acy2-rdat%fq3(3,3)           &
     &                                    +rdat%fq2(1,3)               )*xmdt
      rdat%r04(14,1)= rdat%r04(14,1)+(rdat%fq4(4,1)     -rdat%fq3(2,3)* 3            )*xmdty
      rdat%r04(15,1)= rdat%r04(15,1)+(rdat%fq4(5,1)  -rdat%fq3(3,3)* 6 +rdat%fq2(1,3)* 3 )*xmdt

      end subroutine intk_11

! >
! >    @brief   dsds case
! >
! >    @details integration of the dsds case
! >
      subroutine intk_12(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd1, xmd2, xmd3
      real(kind=dp) :: xmd3y, xmd3x, xmd3xy
      real(kind=dp) :: xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: work, fwk
      integer :: i, j

! Generate jtype=12 integrals

      dimension  work(5),fwk(10)

      if(ikl == 0) then
         rdat%r00(1:4,1)= 0.0_dp
         rdat%r01(:,1:4)= 0.0_dp
         rdat%r02(:,1:5)= 0.0_dp
         rdat%r03(:,1:2)= 0.0_dp
         rdat%r04(:,1)= 0.0_dp
         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd1= xmd2
      xmd3=-xmd1*rdat%y03
      xmdt= xmd1*xmd2

      xmd3y=-xmd3*rdat%acy
      xmd3x=-xmd3*rdat%aqx
      xmd3xy=xmd3*rdat%aqxy
      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      work(1)= xmd1
      work(2)= rdat%y03 *rdat%y03
      work(3)= xmdt
      work(4)= xmdt
      work(5)= xmd3

      rdat%r00(1,1)= rdat%r00(1,1)+rdat%fq0(1)*work(1)
      rdat%r00(2,1)= rdat%r00(2,1)+rdat%fq0(1)*work(2)
      rdat%r00(3,1)= rdat%r00(3,1)+rdat%fq0(2)*work(1)
      rdat%r00(4,1)= rdat%r00(4,1)+rdat%fq0(2)*work(2)

      do i= 1, 2
         rdat%r01(1,i)= rdat%r01(1,i)+rdat%fq1(1,i)*xmd3x
         rdat%r01(2,i)= rdat%r01(2,i)+rdat%fq1(1,i)*xmd3y
         rdat%r01(3,i)= rdat%r01(3,i)+rdat%fq1(2,i)*xmd3
      enddo
         fwk(1)=-rdat%fq1(1,3)*rdat%aqx
         fwk(2)=-rdat%fq1(1,3)*rdat%acy
         fwk(3)= rdat%fq1(2,3)
      do i= 3, 4
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1)*work(i- 2)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2)*work(i- 2)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3)*work(i- 2)
      enddo

      do i= 1, 3
         rdat%r02(1,i)= rdat%r02(1,i)+(rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i))*work(i+ 2)
         rdat%r02(2,i)= rdat%r02(2,i)+(rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i))*work(i+ 2)
         rdat%r02(3,i)= rdat%r02(3,i)+(rdat%fq2(3,i)     -rdat%fq1(1,i))*work(i+ 2)
         rdat%r02(4,i)= rdat%r02(4,i)+ rdat%fq2(1,i)*rdat%aqxy           *work(i+ 2)
         rdat%r02(5,i)= rdat%r02(5,i)- rdat%fq2(2,i)*rdat%aqx            *work(i+ 2)
         rdat%r02(6,i)= rdat%r02(6,i)- rdat%fq2(2,i)*rdat%acy            *work(i+ 2)
      enddo
         fwk(1)= rdat%fq2(1,4)*rdat%aqx2-rdat%fq1(1,4)
         fwk(2)= rdat%fq2(1,4)*rdat%acy2-rdat%fq1(1,4)
         fwk(3)= rdat%fq2(3,4)     -rdat%fq1(1,4)
         fwk(4)= rdat%fq2(1,4)*rdat%aqxy
         fwk(5)=-rdat%fq2(2,4)*rdat%aqx
         fwk(6)=-rdat%fq2(2,4)*rdat%acy
      do i= 4, 5
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j)*work(i-3)
         end do
      end do

      rdat%r03( 1,1)= rdat%r03( 1,1)+(rdat%fq3(1,1)*rdat%aqx2-rdat%fq2(1,3)* 3 )*xmdtx
      rdat%r03( 2,1)= rdat%r03( 2,1)+(rdat%fq3(1,1)*rdat%aqx2-rdat%fq2(1,3)    )*xmdty
      rdat%r03( 3,1)= rdat%r03( 3,1)+(rdat%fq3(2,1)*rdat%aqx2-rdat%fq2(2,3)    )*xmdt
      rdat%r03( 4,1)= rdat%r03( 4,1)+(rdat%fq3(1,1)*rdat%acy2-rdat%fq2(1,3)    )*xmdtx
      rdat%r03( 5,1)= rdat%r03( 5,1)+ rdat%fq3(2,1)                    *xmdtxy
      rdat%r03( 6,1)= rdat%r03( 6,1)+(rdat%fq3(3,1)     -rdat%fq2(1,3)    )*xmdtx
      rdat%r03( 7,1)= rdat%r03( 7,1)+(rdat%fq3(1,1)*rdat%acy2-rdat%fq2(1,3)* 3 )*xmdty
      rdat%r03( 8,1)= rdat%r03( 8,1)+(rdat%fq3(2,1)*rdat%acy2-rdat%fq2(2,3)    )*xmdt
      rdat%r03( 9,1)= rdat%r03( 9,1)+(rdat%fq3(3,1)     -rdat%fq2(1,3)    )*xmdty
      rdat%r03(10,1)= rdat%r03(10,1)+(rdat%fq3(4,1)     -rdat%fq2(2,3)* 3 )*xmdt

      rdat%r03( 1,2)= rdat%r03( 1,2)+(rdat%fq3(1,2)*rdat%aqx2-rdat%fq2(1,4)* 3 )*xmd3x
      rdat%r03( 2,2)= rdat%r03( 2,2)+(rdat%fq3(1,2)*rdat%aqx2-rdat%fq2(1,4)    )*xmd3y
      rdat%r03( 3,2)= rdat%r03( 3,2)+(rdat%fq3(2,2)*rdat%aqx2-rdat%fq2(2,4)    )*xmd3
      rdat%r03( 4,2)= rdat%r03( 4,2)+(rdat%fq3(1,2)*rdat%acy2-rdat%fq2(1,4)    )*xmd3x
      rdat%r03( 5,2)= rdat%r03( 5,2)+ rdat%fq3(2,2)                    *xmd3xy
      rdat%r03( 6,2)= rdat%r03( 6,2)+(rdat%fq3(3,2)     -rdat%fq2(1,4)    )*xmd3x
      rdat%r03( 7,2)= rdat%r03( 7,2)+(rdat%fq3(1,2)*rdat%acy2-rdat%fq2(1,4)* 3 )*xmd3y
      rdat%r03( 8,2)= rdat%r03( 8,2)+(rdat%fq3(2,2)*rdat%acy2-rdat%fq2(2,4)    )*xmd3
      rdat%r03( 9,2)= rdat%r03( 9,2)+(rdat%fq3(3,2)     -rdat%fq2(1,4)    )*xmd3y
      rdat%r03(10,2)= rdat%r03(10,2)+(rdat%fq3(4,2)     -rdat%fq2(2,4)* 3 )*xmd3

      aqx4= rdat%aqx2*rdat%aqx2
      acy4= rdat%acy2*rdat%acy2
      x2y2= rdat%aqx2*rdat%acy2
      q2c2= rdat%aqx2+rdat%acy2
      rdat%r04( 1,1)= rdat%r04( 1,1)+(rdat%fq4(1,1)*aqx4-rdat%fq3(1,2)* 6 *rdat%aqx2                 &
     &                                    +rdat%fq2(1,4)* 3            )*xmdt
      rdat%r04( 2,1)= rdat%r04( 2,1)+(rdat%fq4(1,1)*rdat%aqx2-rdat%fq3(1,2)* 3            )*xmdtxy
      rdat%r04( 3,1)= rdat%r04( 3,1)+(rdat%fq4(2,1)*rdat%aqx2-rdat%fq3(2,2)* 3            )*xmdtx
      rdat%r04( 4,1)= rdat%r04( 4,1)+(rdat%fq4(1,1)*x2y2-rdat%fq3(1,2)*q2c2+rdat%fq2(1,4))*xmdt
      rdat%r04( 5,1)= rdat%r04( 5,1)+(rdat%fq4(2,1)*rdat%aqx2-rdat%fq3(2,2)               )*xmdty
      rdat%r04( 6,1)= rdat%r04( 6,1)+(rdat%fq4(3,1)*rdat%aqx2-rdat%fq3(1,2)*rdat%aqx2-rdat%fq3(3,2)           &
     &                                    +rdat%fq2(1,4)               )*xmdt
      rdat%r04( 7,1)= rdat%r04( 7,1)+(rdat%fq4(1,1)*rdat%acy2-rdat%fq3(1,2)* 3            )*xmdtxy
      rdat%r04( 8,1)= rdat%r04( 8,1)+(rdat%fq4(2,1)*rdat%acy2-rdat%fq3(2,2)               )*xmdtx
      rdat%r04( 9,1)= rdat%r04( 9,1)+(rdat%fq4(3,1)     -rdat%fq3(1,2)               )*xmdtxy
      rdat%r04(10,1)= rdat%r04(10,1)+(rdat%fq4(4,1)     -rdat%fq3(2,2)* 3            )*xmdtx
      rdat%r04(11,1)= rdat%r04(11,1)+(rdat%fq4(1,1)*acy4-rdat%fq3(1,2)* 6 *rdat%acy2                 &
     &                                    +rdat%fq2(1,4)* 3            )*xmdt
      rdat%r04(12,1)= rdat%r04(12,1)+(rdat%fq4(2,1)*rdat%acy2-rdat%fq3(2,2)* 3            )*xmdty
      rdat%r04(13,1)= rdat%r04(13,1)+(rdat%fq4(3,1)*rdat%acy2-rdat%fq3(1,2)*rdat%acy2-rdat%fq3(3,2)           &
     &                                    +rdat%fq2(1,4)               )*xmdt
      rdat%r04(14,1)= rdat%r04(14,1)+(rdat%fq4(4,1)     -rdat%fq3(2,2)* 3            )*xmdty
      rdat%r04(15,1)= rdat%r04(15,1)+(rdat%fq4(5,1)  -rdat%fq3(3,2)* 6 +rdat%fq2(1,4)* 3 )*xmdt

      end subroutine intk_12

! >
! >    @brief   dspp case
! >
! >    @details integration of the dspp case
! >
      subroutine intk_13(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd1, xmd2, xmd3
      real(kind=dp) :: xmd3y, xmd3x, xmd3xy
      real(kind=dp) :: xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: fqd11, fqd12
      real(kind=dp) :: work, fwk
      integer :: i, j

! Generate jtype=13 integrals

      dimension  work(2),fwk(6,4)

      if(ikl == 0) then
         rdat%r00(:,1:2)= 0.0_dp
         rdat%r01(:,1:13)= 0.0_dp
         rdat%r02(:,1:11)= 0.0_dp
         rdat%r03(:,1:5)= 0.0_dp
         rdat%r04(:,1)= 0.0_dp
         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd1= xmd2
      xmd3=-xmd1*rdat%y03
      xmdt= xmd1*xmd2

      xmd3y=-xmd3*rdat%acy
      xmd3x=-xmd3*rdat%aqx
      xmd3xy=xmd3*rdat%aqxy
      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      work(1)= xmd1
      work(2)= rdat%y03 *rdat%y03
      do i= 1, 2
         do j= 1, 5
            rdat%r00(j,i)= rdat%r00(j,i)+rdat%fq0(j)*work(i)
         end do
      end do

      do i= 1, 5
         rdat%r01(1,i)= rdat%r01(1,i)+rdat%fq1(1,i)*xmd3x
         rdat%r01(2,i)= rdat%r01(2,i)+rdat%fq1(1,i)*xmd3y
         rdat%r01(3,i)= rdat%r01(3,i)+rdat%fq1(2,i)*xmd3
      enddo
      do i= 1, 3
         fwk(1,i)=-rdat%fq1(1,i+5)*rdat%aqx
         fwk(2,i)=-rdat%fq1(1,i+5)*rdat%acy
         fwk(3,i)= rdat%fq1(2,i+5)
      enddo
      fqd11 = rdat%fq1(1,5)*rdat%rab +rdat%fq1(1,8)
      fqd12 = rdat%fq1(2,5)*rdat%rab +rdat%fq1(2,8)
         fwk(1,4)=-fqd11*rdat%aqx
         fwk(2,4)=-fqd11*rdat%acy
         fwk(3,4)= fqd12
      do i= 6, 7
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,1)*work(i- 5)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,1)*work(i- 5)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,1)*work(i- 5)
      enddo
      do i= 8, 9
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,2)*work(i- 7)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,2)*work(i- 7)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,2)*work(i- 7)
      enddo
      do i=10,11
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,3)*work(i- 9)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,3)*work(i- 9)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,3)*work(i- 9)
      enddo
      do i=12,13
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,4)*work(i-11)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,4)*work(i-11)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,4)*work(i-11)
      enddo

      do i= 1, 5
         rdat%r02(1,i)= rdat%r02(1,i)+(rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i))*xmdt
         rdat%r02(2,i)= rdat%r02(2,i)+(rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i))*xmdt
         rdat%r02(3,i)= rdat%r02(3,i)+(rdat%fq2(3,i)     -rdat%fq1(1,i))*xmdt
         rdat%r02(4,i)= rdat%r02(4,i)+ rdat%fq2(1,i)                *xmdtxy
         rdat%r02(5,i)= rdat%r02(5,i)+ rdat%fq2(2,i)                *xmdtx
         rdat%r02(6,i)= rdat%r02(6,i)+ rdat%fq2(2,i)                *xmdty
      enddo
      do i= 6, 8
         rdat%r02(1,i)= rdat%r02(1,i)+(rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i))*xmd3
         rdat%r02(2,i)= rdat%r02(2,i)+(rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i))*xmd3
         rdat%r02(3,i)= rdat%r02(3,i)+(rdat%fq2(3,i)     -rdat%fq1(1,i))*xmd3
         rdat%r02(4,i)= rdat%r02(4,i)+ rdat%fq2(1,i)                *xmd3xy
         rdat%r02(5,i)= rdat%r02(5,i)+ rdat%fq2(2,i)                *xmd3x
         rdat%r02(6,i)= rdat%r02(6,i)+ rdat%fq2(2,i)                *xmd3y
      enddo
      rdat%fq2(1,5)= rdat%fq2(1,5)*rdat%rab +rdat%fq2(1,8)
      rdat%fq2(2,5)= rdat%fq2(2,5)*rdat%rab +rdat%fq2(2,8)
      rdat%fq2(3,5)= rdat%fq2(3,5)*rdat%rab +rdat%fq2(3,8)
         rdat%r02(1,9)= rdat%r02(1,9)+(rdat%fq2(1,5)*rdat%aqx2-fqd11)*xmd3
         rdat%r02(2,9)= rdat%r02(2,9)+(rdat%fq2(1,5)*rdat%acy2-fqd11)*xmd3
         rdat%r02(3,9)= rdat%r02(3,9)+(rdat%fq2(3,5)     -fqd11)*xmd3
         rdat%r02(4,9)= rdat%r02(4,9)+ rdat%fq2(1,5)            *xmd3xy
         rdat%r02(5,9)= rdat%r02(5,9)+ rdat%fq2(2,5)            *xmd3x
         rdat%r02(6,9)= rdat%r02(6,9)+ rdat%fq2(2,5)            *xmd3y
         fwk(1,1)= rdat%fq2(1,9)*rdat%aqx2-rdat%fq1(1,9)
         fwk(2,1)= rdat%fq2(1,9)*rdat%acy2-rdat%fq1(1,9)
         fwk(3,1)= rdat%fq2(3,9)     -rdat%fq1(1,9)
         fwk(4,1)= rdat%fq2(1,9)*rdat%aqxy
         fwk(5,1)=-rdat%fq2(2,9)*rdat%aqx
         fwk(6,1)=-rdat%fq2(2,9)*rdat%acy
      do i=10,11
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,1)*work(i- 9)
         end do
      end do

      rdat%fq3(1,1)= rdat%fq3(1,1)*rdat%rab +rdat%fq3(1,4)
      rdat%fq3(2,1)= rdat%fq3(2,1)*rdat%rab +rdat%fq3(2,4)
      rdat%fq3(3,1)= rdat%fq3(3,1)*rdat%rab +rdat%fq3(3,4)
      rdat%fq3(4,1)= rdat%fq3(4,1)*rdat%rab +rdat%fq3(4,4)
      do i= 1, 4
         rdat%r03( 1,i)= rdat%r03( 1,i)+(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i+4)* 3 )*xmdtx
         rdat%r03( 2,i)= rdat%r03( 2,i)+(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i+4)    )*xmdty
         rdat%r03( 3,i)= rdat%r03( 3,i)+(rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i+4)    )*xmdt
         rdat%r03( 4,i)= rdat%r03( 4,i)+(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i+4)    )*xmdtx
         rdat%r03( 5,i)= rdat%r03( 5,i)+ rdat%fq3(2,i)                      *xmdtxy
         rdat%r03( 6,i)= rdat%r03( 6,i)+(rdat%fq3(3,i)     -rdat%fq2(1,i+4)    )*xmdtx
         rdat%r03( 7,i)= rdat%r03( 7,i)+(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i+4)* 3 )*xmdty
         rdat%r03( 8,i)= rdat%r03( 8,i)+(rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i+4)    )*xmdt
         rdat%r03( 9,i)= rdat%r03( 9,i)+(rdat%fq3(3,i)     -rdat%fq2(1,i+4)    )*xmdty
         rdat%r03(10,i)= rdat%r03(10,i)+(rdat%fq3(4,i)     -rdat%fq2(2,i+4)* 3 )*xmdt
      enddo
         rdat%r03( 1,5)= rdat%r03( 1,5)+(rdat%fq3(1,5)*rdat%aqx2-rdat%fq2(1,9)* 3 )*xmd3x
         rdat%r03( 2,5)= rdat%r03( 2,5)+(rdat%fq3(1,5)*rdat%aqx2-rdat%fq2(1,9)    )*xmd3y
         rdat%r03( 3,5)= rdat%r03( 3,5)+(rdat%fq3(2,5)*rdat%aqx2-rdat%fq2(2,9)    )*xmd3
         rdat%r03( 4,5)= rdat%r03( 4,5)+(rdat%fq3(1,5)*rdat%acy2-rdat%fq2(1,9)    )*xmd3x
         rdat%r03( 5,5)= rdat%r03( 5,5)+ rdat%fq3(2,5)                    *xmd3xy
         rdat%r03( 6,5)= rdat%r03( 6,5)+(rdat%fq3(3,5)     -rdat%fq2(1,9)    )*xmd3x
         rdat%r03( 7,5)= rdat%r03( 7,5)+(rdat%fq3(1,5)*rdat%acy2-rdat%fq2(1,9)* 3 )*xmd3y
         rdat%r03( 8,5)= rdat%r03( 8,5)+(rdat%fq3(2,5)*rdat%acy2-rdat%fq2(2,9)    )*xmd3
         rdat%r03( 9,5)= rdat%r03( 9,5)+(rdat%fq3(3,5)     -rdat%fq2(1,9)    )*xmd3y
         rdat%r03(10,5)= rdat%r03(10,5)+(rdat%fq3(4,5)     -rdat%fq2(2,9)* 3 )*xmd3

      aqx4= rdat%aqx2*rdat%aqx2
      acy4= rdat%acy2*rdat%acy2
      x2y2= rdat%aqx2*rdat%acy2
      q2c2= rdat%aqx2+rdat%acy2
      rdat%r04( 1,1)= rdat%r04( 1,1)+(rdat%fq4(1,1)*aqx4-rdat%fq3(1,5)* 6 *rdat%aqx2                 &
     &                                    +rdat%fq2(1,9)* 3            )*xmdt
      rdat%r04( 2,1)= rdat%r04( 2,1)+(rdat%fq4(1,1)*rdat%aqx2-rdat%fq3(1,5)* 3            )*xmdtxy
      rdat%r04( 3,1)= rdat%r04( 3,1)+(rdat%fq4(2,1)*rdat%aqx2-rdat%fq3(2,5)* 3            )*xmdtx
      rdat%r04( 4,1)= rdat%r04( 4,1)+(rdat%fq4(1,1)*x2y2-rdat%fq3(1,5)*q2c2+rdat%fq2(1,9))*xmdt
      rdat%r04( 5,1)= rdat%r04( 5,1)+(rdat%fq4(2,1)*rdat%aqx2-rdat%fq3(2,5)               )*xmdty
      rdat%r04( 6,1)= rdat%r04( 6,1)+(rdat%fq4(3,1)*rdat%aqx2-rdat%fq3(1,5)*rdat%aqx2                     &
     &                                    -rdat%fq3(3,5)+rdat%fq2(1,9)     )*xmdt
      rdat%r04( 7,1)= rdat%r04( 7,1)+(rdat%fq4(1,1)*rdat%acy2-rdat%fq3(1,5)* 3            )*xmdtxy
      rdat%r04( 8,1)= rdat%r04( 8,1)+(rdat%fq4(2,1)*rdat%acy2-rdat%fq3(2,5)               )*xmdtx
      rdat%r04( 9,1)= rdat%r04( 9,1)+(rdat%fq4(3,1)     -rdat%fq3(1,5)               )*xmdtxy
      rdat%r04(10,1)= rdat%r04(10,1)+(rdat%fq4(4,1)     -rdat%fq3(2,5)* 3            )*xmdtx
      rdat%r04(11,1)= rdat%r04(11,1)+(rdat%fq4(1,1)*acy4-rdat%fq3(1,5)* 6 *rdat%acy2                 &
     &                                    +rdat%fq2(1,9)* 3            )*xmdt
      rdat%r04(12,1)= rdat%r04(12,1)+(rdat%fq4(2,1)*rdat%acy2-rdat%fq3(2,5)* 3            )*xmdty
      rdat%r04(13,1)= rdat%r04(13,1)+(rdat%fq4(3,1)*rdat%acy2-rdat%fq3(1,5)*rdat%acy2-rdat%fq3(3,5)           &
     &                                    +rdat%fq2(1,9)               )*xmdt
      rdat%r04(14,1)= rdat%r04(14,1)+(rdat%fq4(4,1)     -rdat%fq3(2,5)* 3            )*xmdty
      rdat%r04(15,1)= rdat%r04(15,1)+(rdat%fq4(5,1)  -rdat%fq3(3,5)* 6 +rdat%fq2(1,9)* 3 )*xmdt
      end subroutine intk_13

! >
! >    @brief   ddps case
! >
! >    @details integration of the ddps case
! >
      subroutine intk_14(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd2, xmd3, xmd4, xmd6
      real(kind=dp) :: xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: y33, y34, y334, y344, y44
      real(kind=dp) :: work, fwk
      integer :: i, j

! Generate jtype=14 integrals

      dimension  work(15),fwk(15,3)

      if(ikl == 0) then
         rdat%r00(:,1:2)= 0.0_dp
         rdat%r01(:,1:13)= 0.0_dp
         rdat%r02(:,1:12)= 0.0_dp
         rdat%r03(:,1:8)= 0.0_dp
         rdat%r04(:,1:4)= 0.0_dp
         rdat%r05(:,1)= 0.0_dp
         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd3= xmd2
      xmd4= xmd3*xmd2
      xmd6= xmd4*xmd2
      xmdt= xmd6*xmd2
      xmd2= xmd3

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      y33 = rdat%y03 *rdat%y03
      y34 =-rdat%y03 *rdat%y04
      y44 = rdat%y04 *rdat%y04
      y334= y33 *rdat%y04
      y344= y34 *rdat%y04
      work( 1)= xmd4
      work( 2)= xmd2*y33
      work( 3)= xmd2*y34
      work( 4)= xmd2*y44
      work( 5)= y33 *y44

      work( 6)=-xmd4*rdat%y03
      work( 7)= xmd4*rdat%y04
      work( 8)= xmd2*y334
      work( 9)= xmd2*y344

      work(10)= xmd6
      work(11)= xmd4*y33
      work(12)= xmd4*y34
      work(13)= xmd4*y44

      work(14)=-xmd6*rdat%y03
      work(15)= xmd6*rdat%y04

      do i= 1, 2
         do j= 1, 5
            rdat%r00(j,i)= rdat%r00(j,i)+rdat%fq0(i)*work(j)
         end do
      end do

      do i= 1, 3
         fwk(1,i)=-rdat%fq1(1,i)*rdat%aqx
         fwk(2,i)=-rdat%fq1(1,i)*rdat%acy
         fwk(3,i)= rdat%fq1(2,i)
      enddo
      do i= 1, 4
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,1)*work(i+ 5)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,1)*work(i+ 5)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,1)*work(i+ 5)
      enddo
      do i= 5, 8
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,2)*work(i+ 1)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,2)*work(i+ 1)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,2)*work(i+ 1)
      enddo
      do i= 9,13
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,3)*work(i- 8)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,3)*work(i- 8)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,3)*work(i- 8)
      enddo

      do i= 1, 3
         fwk(1,i)= rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i)
         fwk(2,i)= rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i)
         fwk(3,i)= rdat%fq2(3,i)     -rdat%fq1(1,i)
         fwk(4,i)= rdat%fq2(1,i)*rdat%aqxy
         fwk(5,i)=-rdat%fq2(2,i)*rdat%aqx
         fwk(6,i)=-rdat%fq2(2,i)*rdat%acy
      enddo
      do i= 1, 4
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,1)*work(i+ 9)
         enddo
      enddo
      do i= 5, 8
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,2)*work(i+ 5)
         enddo
      enddo
      do i= 9,12
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,3)*work(i- 3)
         enddo
      enddo

      do i= 1, 3
         fwk( 1,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)* 3 )*rdat%aqx
         fwk( 2,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)    )*rdat%acy
         fwk( 3,i)=  rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i)
         fwk( 4,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)    )*rdat%aqx
         fwk( 5,i)=  rdat%fq3(2,i)                    *rdat%aqxy
         fwk( 6,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%aqx
         fwk( 7,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)* 3 )*rdat%acy
         fwk( 8,i)=  rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i)
         fwk( 9,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%acy
         fwk(10,i)=  rdat%fq3(4,i)     -rdat%fq2(2,i)* 3
      enddo
      do i= 1, 2
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,1)*work(i+13)
         enddo
      enddo
      do i= 3, 4
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,2)*work(i+11)
         enddo
      enddo
      do i= 5, 8
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,3)*work(i+ 5)
         enddo
      enddo

      aqx4= rdat%aqx2*rdat%aqx2
      acy4= rdat%acy2*rdat%acy2
      x2y2= rdat%aqx2*rdat%acy2
      q2c2= rdat%aqx2+rdat%acy2
      do i= 1, 2
         rdat%r04( 1,i)= rdat%r04( 1,i)+(rdat%fq4(1,i)*aqx4-rdat%fq3(1,i)* 6 *rdat%aqx2        &
     &                        +rdat%fq2(1,i)* 3                  )*xmdt
         rdat%r04( 2,i)= rdat%r04( 2,i)+(rdat%fq4(1,i)*rdat%aqx2-rdat%fq3(1,i)* 3   )*xmdtxy
         rdat%r04( 3,i)= rdat%r04( 3,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)* 3   )*xmdtx
         rdat%r04( 4,i)= rdat%r04( 4,i)+(rdat%fq4(1,i)*x2y2-rdat%fq3(1,i)*q2c2            &
     &                        +rdat%fq2(1,i)                     )*xmdt
         rdat%r04( 5,i)= rdat%r04( 5,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)      )*xmdty
         rdat%r04( 6,i)= rdat%r04( 6,i)+(rdat%fq4(3,i)*rdat%aqx2-rdat%fq3(1,i)*rdat%aqx2            &
     &                        -rdat%fq3(3,i)+rdat%fq2(1,i)           )*xmdt
         rdat%r04( 7,i)= rdat%r04( 7,i)+(rdat%fq4(1,i)*rdat%acy2-rdat%fq3(1,i)* 3   )*xmdtxy
         rdat%r04( 8,i)= rdat%r04( 8,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)      )*xmdtx
         rdat%r04( 9,i)= rdat%r04( 9,i)+(rdat%fq4(3,i)     -rdat%fq3(1,i)      )*xmdtxy
         rdat%r04(10,i)= rdat%r04(10,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3   )*xmdtx
         rdat%r04(11,i)= rdat%r04(11,i)+(rdat%fq4(1,i)*acy4-rdat%fq3(1,i)* 6 *rdat%acy2        &
     &                        +rdat%fq2(1,i)* 3                  )*xmdt
         rdat%r04(12,i)= rdat%r04(12,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)* 3   )*xmdty
         rdat%r04(13,i)= rdat%r04(13,i)+(rdat%fq4(3,i)*rdat%acy2-rdat%fq3(1,i)*rdat%acy2            &
     &                        -rdat%fq3(3,i)+rdat%fq2(1,i)           )*xmdt
         rdat%r04(14,i)= rdat%r04(14,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3   )*xmdty
         rdat%r04(15,i)= rdat%r04(15,i)+(rdat%fq4(5,i)-rdat%fq3(3,i)* 6                   &
     &                        +rdat%fq2(1,i)* 3                  )*xmdt
      enddo
      fwk( 1,3)=  rdat%fq4(1,3)*aqx4-rdat%fq3(1,3)* 6 *rdat%aqx2+rdat%fq2(1,3)* 3
      fwk( 2,3)= (rdat%fq4(1,3)*rdat%aqx2-rdat%fq3(1,3)* 3                )*rdat%aqxy
      fwk( 3,3)=-(rdat%fq4(2,3)*rdat%aqx2-rdat%fq3(2,3)* 3                )*rdat%aqx
      fwk( 4,3)=  rdat%fq4(1,3)*x2y2-rdat%fq3(1,3)*q2c2+rdat%fq2(1,3)
      fwk( 5,3)=-(rdat%fq4(2,3)*rdat%aqx2-rdat%fq3(2,3)                   )*rdat%acy
      fwk( 6,3)=  rdat%fq4(3,3)*rdat%aqx2-rdat%fq3(1,3)*rdat%aqx2-rdat%fq3(3,3)+rdat%fq2(1,3)
      fwk( 7,3)= (rdat%fq4(1,3)*rdat%acy2-rdat%fq3(1,3)* 3                )*rdat%aqxy
      fwk( 8,3)=-(rdat%fq4(2,3)*rdat%acy2-rdat%fq3(2,3)                   )*rdat%aqx
      fwk( 9,3)= (rdat%fq4(3,3)     -rdat%fq3(1,3)                   )*rdat%aqxy
      fwk(10,3)=-(rdat%fq4(4,3)     -rdat%fq3(2,3)* 3                )*rdat%aqx
      fwk(11,3)=  rdat%fq4(1,3)*acy4-rdat%fq3(1,3)* 6 *rdat%acy2+rdat%fq2(1,3)* 3
      fwk(12,3)=-(rdat%fq4(2,3)*rdat%acy2-rdat%fq3(2,3)* 3                )*rdat%acy
      fwk(13,3)=  rdat%fq4(3,3)*rdat%acy2-rdat%fq3(1,3)*rdat%acy2-rdat%fq3(3,3)+rdat%fq2(1,3)
      fwk(14,3)=-(rdat%fq4(4,3)     -rdat%fq3(2,3)* 3                )*rdat%acy
      fwk(15,3)=  rdat%fq4(5,3)     -rdat%fq3(3,3)* 6 +rdat%fq2(1,3)* 3
      do i= 3, 4
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,3)*work(i+11)
         end do
      end do

      rdat%r05( 1,1)= rdat%r05( 1,1)+(rdat%fq5(1,1)*aqx4-rdat%fq4(1,3)*10*rdat%aqx2                 &
     &                 +rdat%fq3(1,3)*15                        )*xmdtx
      rdat%r05( 2,1)= rdat%r05( 2,1)+(rdat%fq5(1,1)*aqx4-rdat%fq4(1,3)* 6 *rdat%aqx2                 &
     &                 +rdat%fq3(1,3)* 3                         )*xmdty
      rdat%r05( 3,1)= rdat%r05( 3,1)+(rdat%fq5(2,1)*aqx4-rdat%fq4(2,3)* 6 *rdat%aqx2                 &
     &                 +rdat%fq3(2,3)* 3                         )*xmdt
      rdat%r05( 4,1)= rdat%r05( 4,1)+(rdat%fq5(1,1)*x2y2-rdat%fq4(1,3)*rdat%aqx2-rdat%fq4(1,3)* 3 *rdat%acy2  &
     &                 +rdat%fq3(1,3)* 3                         )*xmdtx
      rdat%r05( 5,1)= rdat%r05( 5,1)+(rdat%fq5(2,1)*rdat%aqx2-rdat%fq4(2,3)* 3            )*xmdtxy
      rdat%r05( 6,1)= rdat%r05( 6,1)+(rdat%fq5(3,1)*rdat%aqx2-rdat%fq4(1,3)*rdat%aqx2-rdat%fq4(3,3)* 3        &
     &                 +rdat%fq3(1,3)* 3                         )*xmdtx
      rdat%r05( 7,1)= rdat%r05( 7,1)+(rdat%fq5(1,1)*x2y2-rdat%fq4(1,3)* 3 *rdat%aqx2-rdat%fq4(1,3)*rdat%acy2  &
     &                 +rdat%fq3(1,3)* 3                         )*xmdty
      rdat%r05( 8,1)= rdat%r05( 8,1)+(rdat%fq5(2,1)*x2y2-rdat%fq4(2,3)*q2c2+rdat%fq3(2,3))*xmdt
      rdat%r05( 9,1)= rdat%r05( 9,1)+(rdat%fq5(3,1)*rdat%aqx2-rdat%fq4(1,3)*rdat%aqx2-rdat%fq4(3,3)           &
     &                 +rdat%fq3(1,3)                            )*xmdty
      rdat%r05(10,1)= rdat%r05(10,1)+(rdat%fq5(4,1)*rdat%aqx2-rdat%fq4(2,3)* 3 *rdat%aqx2-rdat%fq4(4,3)       &
     &                 +rdat%fq3(2,3)* 3                         )*xmdt
      rdat%r05(11,1)= rdat%r05(11,1)+(rdat%fq5(1,1)*acy4-rdat%fq4(1,3)* 6 *rdat%acy2                 &
     &                 +rdat%fq3(1,3)* 3                         )*xmdtx
      rdat%r05(12,1)= rdat%r05(12,1)+(rdat%fq5(2,1)*rdat%acy2-rdat%fq4(2,3)* 3            )*xmdtxy
      rdat%r05(13,1)= rdat%r05(13,1)+(rdat%fq5(3,1)*rdat%acy2-rdat%fq4(1,3)*rdat%acy2-rdat%fq4(3,3)           &
     &                 +rdat%fq3(1,3)                            )*xmdtx
      rdat%r05(14,1)= rdat%r05(14,1)+(rdat%fq5(4,1)-rdat%fq4(2,3)* 3                 )*xmdtxy
      rdat%r05(15,1)= rdat%r05(15,1)+(rdat%fq5(5,1)-rdat%fq4(3,3)* 6 +rdat%fq3(1,3)* 3   )*xmdtx
      rdat%r05(16,1)= rdat%r05(16,1)+(rdat%fq5(1,1)*acy4-rdat%fq4(1,3)*10*rdat%acy2                 &
     &                 +rdat%fq3(1,3)*15                        )*xmdty
      rdat%r05(17,1)= rdat%r05(17,1)+(rdat%fq5(2,1)*acy4-rdat%fq4(2,3)* 6 *rdat%acy2                 &
     &                 +rdat%fq3(2,3)* 3                         )*xmdt
      rdat%r05(18,1)= rdat%r05(18,1)+(rdat%fq5(3,1)*rdat%acy2-rdat%fq4(1,3)*rdat%acy2-rdat%fq4(3,3)* 3        &
     &                 +rdat%fq3(1,3)* 3                         )*xmdty
      rdat%r05(19,1)= rdat%r05(19,1)+(rdat%fq5(4,1)*rdat%acy2-rdat%fq4(2,3)* 3 *rdat%acy2-rdat%fq4(4,3)       &
     &                 +rdat%fq3(2,3)* 3                         )*xmdt
      rdat%r05(20,1)= rdat%r05(20,1)+(rdat%fq5(5,1)-rdat%fq4(3,3)* 6 +rdat%fq3(1,3)* 3   )*xmdty
      rdat%r05(21,1)= rdat%r05(21,1)+(rdat%fq5(6,1)-rdat%fq4(4,3)*10+rdat%fq3(2,3)*15  )*xmdt

      end subroutine intk_14

! >
! >    @brief   dpds case
! >
! >    @details integration of the dpds case
! >
      subroutine intk_15(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd1, xmd2, xmd3, xmd4, xmd5
      real(kind=dp) :: xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: y33, y34
      real(kind=dp) :: work, fwk
      integer :: i, j

! Generate jtype=15 integrals

      dimension  work(12),fwk(15,4)

      if(ikl == 0) then
         rdat%r00(:,1:2)= 0.0_dp
         rdat%r01(:,1:13)= 0.0_dp
         rdat%r02(:,1:15)= 0.0_dp
         rdat%r03(:,1:9)= 0.0_dp
         rdat%r04(:,1:4)= 0.0_dp
         rdat%r05(:,1)= 0.0_dp
         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd1= xmd2
      xmd3= xmd2
      xmd4= xmd2*xmd2
      xmd5= xmd4
      xmdt= xmd4*xmd3

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      y33 = rdat%y03 *rdat%y03
      y34 =-rdat%y03 *rdat%y04
      work( 1)= xmd1
      work( 2)= y33
      work( 3)=-xmd3*rdat%y03
      work( 4)= xmd3*rdat%y04
      work( 5)= y33 *rdat%y04

      work( 6)=-xmd1*rdat%y03
      work( 7)= xmd5
      work( 8)= xmd3*y33
      work( 9)= xmd3*y34

      work(10)= xmd4
      work(11)=-xmd5*rdat%y03
      work(12)= xmd5*rdat%y04

      do i= 1, 2
         do j= 1, 5
            rdat%r00(j,i)= rdat%r00(j,i)+rdat%fq0(i)*work(j)
         end do
      end do

      do i= 1, 3
         fwk(1,i)=-rdat%fq1(1,i)*rdat%aqx
         fwk(2,i)=-rdat%fq1(1,i)*rdat%acy
         fwk(3,i)= rdat%fq1(2,i)
      enddo
      do i= 1, 4
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,1)*work(i+ 5)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,1)*work(i+ 5)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,1)*work(i+ 5)
      enddo
      do i= 5, 8
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,2)*work(i+ 1)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,2)*work(i+ 1)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,2)*work(i+ 1)
      enddo
      do i= 9,13
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,3)*work(i- 8)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,3)*work(i- 8)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,3)*work(i- 8)
      enddo

      do i= 1, 4
         fwk(1,i)= rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i)
         fwk(2,i)= rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i)
         fwk(3,i)= rdat%fq2(3,i)     -rdat%fq1(1,i)
         fwk(4,i)= rdat%fq2(1,i)*rdat%aqxy
         fwk(5,i)=-rdat%fq2(2,i)*rdat%aqx
         fwk(6,i)=-rdat%fq2(2,i)*rdat%acy
      enddo
      do i= 1, 3
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,1)*work(i+ 9)
         enddo
      enddo
      do i= 4, 6
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,2)*work(i+ 6)
         enddo
      enddo
      do i= 7,10
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,3)*work(i- 1)
         enddo
      enddo
      do i=11,15
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,4)*work(i-10)
         enddo
      enddo

      do i= 1, 2
         rdat%r03( 1,i)= rdat%r03( 1,i)+(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)* 3 )*xmdtx
         rdat%r03( 2,i)= rdat%r03( 2,i)+(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)    )*xmdty
         rdat%r03( 3,i)= rdat%r03( 3,i)+(rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i)    )*xmdt
         rdat%r03( 4,i)= rdat%r03( 4,i)+(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)    )*xmdtx
         rdat%r03( 5,i)= rdat%r03( 5,i)+ rdat%fq3(2,i)                    *xmdtxy
         rdat%r03( 6,i)= rdat%r03( 6,i)+(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*xmdtx
         rdat%r03( 7,i)= rdat%r03( 7,i)+(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)* 3 )*xmdty
         rdat%r03( 8,i)= rdat%r03( 8,i)+(rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i)    )*xmdt
         rdat%r03( 9,i)= rdat%r03( 9,i)+(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*xmdty
         rdat%r03(10,i)= rdat%r03(10,i)+(rdat%fq3(4,i)     -rdat%fq2(2,i)* 3 )*xmdt
      enddo
      do i= 3, 4
         fwk( 1,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)* 3 )*rdat%aqx
         fwk( 2,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)    )*rdat%acy
         fwk( 3,i)=  rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i)
         fwk( 4,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)    )*rdat%aqx
         fwk( 5,i)=  rdat%fq3(2,i)                    *rdat%aqxy
         fwk( 6,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%aqx
         fwk( 7,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)* 3 )*rdat%acy
         fwk( 8,i)=  rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i)
         fwk( 9,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%acy
         fwk(10,i)=  rdat%fq3(4,i)     -rdat%fq2(2,i)* 3
      enddo
      do i= 3, 5
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,3)*work(i+ 7)
         enddo
      enddo
      do i= 6, 9
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,4)*work(i)
         enddo
      enddo

      aqx4= rdat%aqx2*rdat%aqx2
      acy4= rdat%acy2*rdat%acy2
      x2y2= rdat%aqx2*rdat%acy2
      q2c2= rdat%aqx2+rdat%acy2
      rdat%r04( 1,1)= rdat%r04( 1,1)+(rdat%fq4(1,1)*aqx4-rdat%fq3(1,3)* 6 *rdat%aqx2           &
     &                                    +rdat%fq2(1,3)* 3      )*xmdt
      rdat%r04( 2,1)= rdat%r04( 2,1)+(rdat%fq4(1,1)*rdat%aqx2-rdat%fq3(1,3)* 3      )*xmdtxy
      rdat%r04( 3,1)= rdat%r04( 3,1)+(rdat%fq4(2,1)*rdat%aqx2-rdat%fq3(2,3)* 3      )*xmdtx
      rdat%r04( 4,1)= rdat%r04( 4,1)+(rdat%fq4(1,1)*x2y2-rdat%fq3(1,3)*q2c2               &
     &                                    +rdat%fq2(1,3)         )*xmdt
      rdat%r04( 5,1)= rdat%r04( 5,1)+(rdat%fq4(2,1)*rdat%aqx2-rdat%fq3(2,3)         )*xmdty
      rdat%r04( 6,1)= rdat%r04( 6,1)+(rdat%fq4(3,1)*rdat%aqx2-rdat%fq3(1,3)*rdat%aqx2               &
     &                     -rdat%fq3(3,3)+rdat%fq2(1,3)              )*xmdt
      rdat%r04( 7,1)= rdat%r04( 7,1)+(rdat%fq4(1,1)*rdat%acy2-rdat%fq3(1,3)* 3      )*xmdtxy
      rdat%r04( 8,1)= rdat%r04( 8,1)+(rdat%fq4(2,1)*rdat%acy2-rdat%fq3(2,3)         )*xmdtx
      rdat%r04( 9,1)= rdat%r04( 9,1)+(rdat%fq4(3,1)     -rdat%fq3(1,3)         )*xmdtxy
      rdat%r04(10,1)= rdat%r04(10,1)+(rdat%fq4(4,1)     -rdat%fq3(2,3)* 3      )*xmdtx
      rdat%r04(11,1)= rdat%r04(11,1)+(rdat%fq4(1,1)*acy4-rdat%fq3(1,3)* 6 *rdat%acy2           &
     &                     +rdat%fq2(1,3)* 3                     )*xmdt
      rdat%r04(12,1)= rdat%r04(12,1)+(rdat%fq4(2,1)*rdat%acy2-rdat%fq3(2,3)* 3      )*xmdty
      rdat%r04(13,1)= rdat%r04(13,1)+(rdat%fq4(3,1)*rdat%acy2-rdat%fq3(1,3)*rdat%acy2               &
     &                     -rdat%fq3(3,3)+rdat%fq2(1,3)              )*xmdt
      rdat%r04(14,1)= rdat%r04(14,1)+(rdat%fq4(4,1)     -rdat%fq3(2,3)* 3      )*xmdty
      rdat%r04(15,1)= rdat%r04(15,1)+(rdat%fq4(5,1)-rdat%fq3(3,3)* 6                      &
     &                     +rdat%fq2(1,3)* 3                     )*xmdt
      fwk( 1,2)=  rdat%fq4(1,2)*aqx4-rdat%fq3(1,4)* 6 *rdat%aqx2+rdat%fq2(1,4)* 3
      fwk( 2,2)= (rdat%fq4(1,2)*rdat%aqx2-rdat%fq3(1,4)* 3                )*rdat%aqxy
      fwk( 3,2)=-(rdat%fq4(2,2)*rdat%aqx2-rdat%fq3(2,4)* 3                )*rdat%aqx
      fwk( 4,2)=  rdat%fq4(1,2)*x2y2-rdat%fq3(1,4)*q2c2    +rdat%fq2(1,4)
      fwk( 5,2)=-(rdat%fq4(2,2)*rdat%aqx2-rdat%fq3(2,4)                   )*rdat%acy
      fwk( 6,2)=  rdat%fq4(3,2)*rdat%aqx2-rdat%fq3(1,4)*rdat%aqx2-rdat%fq3(3,4)+rdat%fq2(1,4)
      fwk( 7,2)= (rdat%fq4(1,2)*rdat%acy2-rdat%fq3(1,4)* 3                )*rdat%aqxy
      fwk( 8,2)=-(rdat%fq4(2,2)*rdat%acy2-rdat%fq3(2,4)                   )*rdat%aqx
      fwk( 9,2)= (rdat%fq4(3,2)     -rdat%fq3(1,4)                   )*rdat%aqxy
      fwk(10,2)=-(rdat%fq4(4,2)     -rdat%fq3(2,4)* 3                )*rdat%aqx
      fwk(11,2)=  rdat%fq4(1,2)*acy4-rdat%fq3(1,4)* 6 *rdat%acy2+rdat%fq2(1,4)* 3
      fwk(12,2)=-(rdat%fq4(2,2)*rdat%acy2-rdat%fq3(2,4)* 3                )*rdat%acy
      fwk(13,2)=  rdat%fq4(3,2)*rdat%acy2-rdat%fq3(1,4)*rdat%acy2-rdat%fq3(3,4)+rdat%fq2(1,4)
      fwk(14,2)=-(rdat%fq4(4,2)     -rdat%fq3(2,4)* 3                )*rdat%acy
      fwk(15,2)=  rdat%fq4(5,2)     -rdat%fq3(3,4)* 6 +rdat%fq2(1,4)* 3
      do i= 2, 4
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,2)*work(i+ 8)
         end do
      end do

      rdat%r05( 1,1)= rdat%r05( 1,1)+(rdat%fq5(1,1)*aqx4-rdat%fq4(1,2)*10*rdat%aqx2                 &
     &                 +rdat%fq3(1,4)*15                        )*xmdtx
      rdat%r05( 2,1)= rdat%r05( 2,1)+(rdat%fq5(1,1)*aqx4-rdat%fq4(1,2)* 6 *rdat%aqx2                 &
     &                 +rdat%fq3(1,4)* 3                         )*xmdty
      rdat%r05( 3,1)= rdat%r05( 3,1)+(rdat%fq5(2,1)*aqx4-rdat%fq4(2,2)* 6 *rdat%aqx2                 &
     &                 +rdat%fq3(2,4)* 3                         )*xmdt
      rdat%r05( 4,1)= rdat%r05( 4,1)+(rdat%fq5(1,1)*x2y2-rdat%fq4(1,2)*rdat%aqx2-rdat%fq4(1,2)* 3 *rdat%acy2  &
     &                 +rdat%fq3(1,4)* 3                         )*xmdtx
      rdat%r05( 5,1)= rdat%r05( 5,1)+(rdat%fq5(2,1)*rdat%aqx2-rdat%fq4(2,2)* 3            )*xmdtxy
      rdat%r05( 6,1)= rdat%r05( 6,1)+(rdat%fq5(3,1)*rdat%aqx2-rdat%fq4(1,2)*rdat%aqx2-rdat%fq4(3,2)* 3        &
     &                 +rdat%fq3(1,4)* 3                         )*xmdtx
      rdat%r05( 7,1)= rdat%r05( 7,1)+(rdat%fq5(1,1)*x2y2-rdat%fq4(1,2)* 3 *rdat%aqx2-rdat%fq4(1,2)*rdat%acy2  &
     &                 +rdat%fq3(1,4)* 3                         )*xmdty
      rdat%r05( 8,1)= rdat%r05( 8,1)+(rdat%fq5(2,1)*x2y2-rdat%fq4(2,2)*q2c2+rdat%fq3(2,4))*xmdt
      rdat%r05( 9,1)= rdat%r05( 9,1)+(rdat%fq5(3,1)*rdat%aqx2-rdat%fq4(1,2)*rdat%aqx2-rdat%fq4(3,2)           &
     &                 +rdat%fq3(1,4)                            )*xmdty
      rdat%r05(10,1)= rdat%r05(10,1)+(rdat%fq5(4,1)*rdat%aqx2-rdat%fq4(2,2)* 3 *rdat%aqx2-rdat%fq4(4,2)       &
     &                 +rdat%fq3(2,4)* 3                         )*xmdt
      rdat%r05(11,1)= rdat%r05(11,1)+(rdat%fq5(1,1)*acy4-rdat%fq4(1,2)* 6 *rdat%acy2                 &
     &                 +rdat%fq3(1,4)* 3                         )*xmdtx
      rdat%r05(12,1)= rdat%r05(12,1)+(rdat%fq5(2,1)*rdat%acy2-rdat%fq4(2,2)* 3            )*xmdtxy
      rdat%r05(13,1)= rdat%r05(13,1)+(rdat%fq5(3,1)*rdat%acy2-rdat%fq4(1,2)*rdat%acy2-rdat%fq4(3,2)           &
     &                 +rdat%fq3(1,4)                            )*xmdtx
      rdat%r05(14,1)= rdat%r05(14,1)+(rdat%fq5(4,1)-rdat%fq4(2,2)* 3                 )*xmdtxy
      rdat%r05(15,1)= rdat%r05(15,1)+(rdat%fq5(5,1)-rdat%fq4(3,2)* 6 +rdat%fq3(1,4)* 3   )*xmdtx
      rdat%r05(16,1)= rdat%r05(16,1)+(rdat%fq5(1,1)*acy4-rdat%fq4(1,2)*10*rdat%acy2                 &
     &                 +rdat%fq3(1,4)*15                        )*xmdty
      rdat%r05(17,1)= rdat%r05(17,1)+(rdat%fq5(2,1)*acy4-rdat%fq4(2,2)* 6 *rdat%acy2                 &
     &                 +rdat%fq3(2,4)* 3                         )*xmdt
      rdat%r05(18,1)= rdat%r05(18,1)+(rdat%fq5(3,1)*rdat%acy2-rdat%fq4(1,2)*rdat%acy2-rdat%fq4(3,2)* 3        &
     &                 +rdat%fq3(1,4)* 3                         )*xmdty
      rdat%r05(19,1)= rdat%r05(19,1)+(rdat%fq5(4,1)*rdat%acy2-rdat%fq4(2,2)* 3 *rdat%acy2-rdat%fq4(4,2)       &
     &                 +rdat%fq3(2,4)* 3                         )*xmdt
      rdat%r05(20,1)= rdat%r05(20,1)+(rdat%fq5(5,1)-rdat%fq4(3,2)* 6 +rdat%fq3(1,4)* 3   )*xmdty
      rdat%r05(21,1)= rdat%r05(21,1)+(rdat%fq5(6,1)-rdat%fq4(4,2)*10+rdat%fq3(2,4)*15  )*xmdt

      end subroutine intk_15

! >
! >    @brief   dppp case
! >
! >    @details integration of the dppp case
! >
      subroutine intk_16(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd1, xmd2, xmd3, xmd4, xmd5
      real(kind=dp) :: xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: y33, y34
      real(kind=dp) :: fqd11, fqd12, fqd21, fqd22, fqd23
      real(kind=dp) :: work, fwk
      integer :: i, j

! Generate jtype=16 integrals

      dimension  work(12),fwk(15,10)

      if(ikl == 0) then
         rdat%r00= 0.0_dp
         rdat%r01= 0.0_dp
         rdat%r02(:,1:36)= 0.0_dp
         rdat%r03(:,1:21)= 0.0_dp
         rdat%r04(:,1:7)= 0.0_dp
         rdat%r05(:,1)= 0.0_dp
         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd1= xmd2
      xmd3= xmd2
      xmd4= xmd2*xmd2
      xmd5= xmd4
      xmdt= xmd4*xmd3

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      y33 = rdat%y03 *rdat%y03
      y34 =-rdat%y03 *rdat%y04
      work( 1)= xmd1
      work( 2)= y33
      work( 3)=-xmd3*rdat%y03
      work( 4)= xmd3*rdat%y04
      work( 5)= y33 *rdat%y04

      work( 6)=-xmd1*rdat%y03
      work( 7)= xmd5
      work( 8)= xmd3*y33
      work( 9)= xmd3*y34

      work(10)= xmd4
      work(11)=-xmd5*rdat%y03
      work(12)= xmd5*rdat%y04

      do i= 1, 5
         do j= 1, 5
            rdat%r00(j,i)= rdat%r00(j,i)+rdat%fq0(i)*work(j)
         end do
      end do

      do i= 1, 8
         fwk(1,i)=-rdat%fq1(1,i)*rdat%aqx
         fwk(2,i)=-rdat%fq1(1,i)*rdat%acy
         fwk(3,i)= rdat%fq1(2,i)
      enddo
      fqd11 = rdat%fq1(1,5)*rdat%rab +rdat%fq1(1,8)
      fqd12 = rdat%fq1(2,5)*rdat%rab +rdat%fq1(2,8)
         fwk(1,9)=-fqd11*rdat%aqx
         fwk(2,9)=-fqd11*rdat%acy
         fwk(3,9)= fqd12
      do i= 1, 4
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,1)*work(i+ 5)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,1)*work(i+ 5)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,1)*work(i+ 5)
      enddo
      do i= 5, 8
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,2)*work(i+ 1)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,2)*work(i+ 1)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,2)*work(i+ 1)
      enddo
      do i= 9,12
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,3)*work(i- 3)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,3)*work(i- 3)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,3)*work(i- 3)
      enddo
      do i=13,16
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,4)*work(i- 7)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,4)*work(i- 7)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,4)*work(i- 7)
      enddo
      do i=17,20
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,5)*work(i-11)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,5)*work(i-11)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,5)*work(i-11)
      enddo
      do i=21,25
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,6)*work(i-20)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,6)*work(i-20)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,6)*work(i-20)
      enddo
      do i=26,30
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,7)*work(i-25)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,7)*work(i-25)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,7)*work(i-25)
      enddo
      do i=31,35
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,8)*work(i-30)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,8)*work(i-30)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,8)*work(i-30)
      enddo
      do i=36,40
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,9)*work(i-35)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,9)*work(i-35)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,9)*work(i-35)
      enddo

      do i= 1, 9
         fwk(1,i)= rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i)
         fwk(2,i)= rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i)
         fwk(3,i)= rdat%fq2(3,i)     -rdat%fq1(1,i)
         fwk(4,i)= rdat%fq2(1,i)*rdat%aqxy
         fwk(5,i)=-rdat%fq2(2,i)*rdat%aqx
         fwk(6,i)=-rdat%fq2(2,i)*rdat%acy
      enddo
      fqd21 = rdat%fq2(1,5)*rdat%rab +rdat%fq2(1,8)
      fqd22 = rdat%fq2(2,5)*rdat%rab +rdat%fq2(2,8)
      fqd23 = rdat%fq2(3,5)*rdat%rab +rdat%fq2(3,8)
         fwk(1,10)= fqd21*rdat%aqx2-fqd11
         fwk(2,10)= fqd21*rdat%acy2-fqd11
         fwk(3,10)= fqd23     -fqd11
         fwk(4,10)= fqd21*rdat%aqxy
         fwk(5,10)=-fqd22*rdat%aqx
         fwk(6,10)=-fqd22*rdat%acy
      do i= 1, 3
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,1)*work(i+ 9)
         enddo
      enddo
      do i= 4, 6
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,2)*work(i+ 6)
         enddo
      enddo
      do i= 7, 9
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,3)*work(i+ 3)
         enddo
      enddo
      do i=10,12
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,4)*work(i)
         enddo
      enddo
      do i=13,15
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,5)*work(i- 3)
         enddo
      enddo
      do i=16,19
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,6)*work(i-10)
         enddo
      enddo
      do i=20,23
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,7)*work(i-14)
         enddo
      enddo
      do i=24,27
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,8)*work(i-18)
         enddo
      enddo
      do i=28,31
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,10)*work(i-22)
         enddo
      enddo
      do i=32,36
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,9)*work(i-31)
         enddo
      enddo

      do i= 1, 5
         rdat%r03( 1,i)= rdat%r03( 1,i)+(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)* 3 )*xmdtx
         rdat%r03( 2,i)= rdat%r03( 2,i)+(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)    )*xmdty
         rdat%r03( 3,i)= rdat%r03( 3,i)+(rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i)    )*xmdt
         rdat%r03( 4,i)= rdat%r03( 4,i)+(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)    )*xmdtx
         rdat%r03( 5,i)= rdat%r03( 5,i)+ rdat%fq3(2,i)                    *xmdtxy
         rdat%r03( 6,i)= rdat%r03( 6,i)+(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*xmdtx
         rdat%r03( 7,i)= rdat%r03( 7,i)+(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)* 3 )*xmdty
         rdat%r03( 8,i)= rdat%r03( 8,i)+(rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i)    )*xmdt
         rdat%r03( 9,i)= rdat%r03( 9,i)+(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*xmdty
         rdat%r03(10,i)= rdat%r03(10,i)+(rdat%fq3(4,i)     -rdat%fq2(2,i)* 3 )*xmdt
      enddo
      rdat%fq2(1,5)= fqd21
      rdat%fq2(2,5)= fqd22
      rdat%fq3(1,5)= rdat%fq3(1,5)*rdat%rab+rdat%fq3(1,8)
      rdat%fq3(2,5)= rdat%fq3(2,5)*rdat%rab+rdat%fq3(2,8)
      rdat%fq3(3,5)= rdat%fq3(3,5)*rdat%rab+rdat%fq3(3,8)
      rdat%fq3(4,5)= rdat%fq3(4,5)*rdat%rab+rdat%fq3(4,8)
      do i= 5, 9
         fwk( 1,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)* 3 )*rdat%aqx
         fwk( 2,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)    )*rdat%acy
         fwk( 3,i)=  rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i)
         fwk( 4,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)    )*rdat%aqx
         fwk( 5,i)=  rdat%fq3(2,i)                    *rdat%aqxy
         fwk( 6,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%aqx
         fwk( 7,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)* 3 )*rdat%acy
         fwk( 8,i)=  rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i)
         fwk( 9,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%acy
         fwk(10,i)=  rdat%fq3(4,i)     -rdat%fq2(2,i)* 3
      enddo
      do i= 6, 8
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,5)*work(i+ 4)
         enddo
      enddo
      do i= 9,11
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,6)*work(i+ 1)
         enddo
      enddo
      do i=12,14
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,7)*work(i- 2)
         enddo
      enddo
      do i=15,17
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,8)*work(i- 5)
         enddo
      enddo
      do i=18,21
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,9)*work(i-12)
         enddo
      enddo

      do j= 1, 5
         rdat%fq4(j,1)= rdat%fq4(j,1)*rdat%rab +rdat%fq4(j,4)
      enddo
      aqx4= rdat%aqx2*rdat%aqx2
      acy4= rdat%acy2*rdat%acy2
      x2y2= rdat%aqx2*rdat%acy2
      q2c2= rdat%aqx2+rdat%acy2
      do i= 1, 4
         rdat%r04( 1,i)= rdat%r04( 1,i)+(rdat%fq4(1,i)*aqx4-rdat%fq3(1,i+4)* 6 *rdat%aqx2      &
     &                        +rdat%fq2(1,i+4)* 3                 )*xmdt
         rdat%r04( 2,i)= rdat%r04( 2,i)+(rdat%fq4(1,i)*rdat%aqx2-rdat%fq3(1,i+4)* 3  )*xmdtxy
         rdat%r04( 3,i)= rdat%r04( 3,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i+4)* 3  )*xmdtx
         rdat%r04( 4,i)= rdat%r04( 4,i)+(rdat%fq4(1,i)*x2y2-rdat%fq3(1,i+4)*q2c2          &
     &                                       +rdat%fq2(1,i+4)     )*xmdt
         rdat%r04( 5,i)= rdat%r04( 5,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i+4)     )*xmdty
         rdat%r04( 6,i)= rdat%r04( 6,i)+(rdat%fq4(3,i)*rdat%aqx2-rdat%fq3(1,i+4)*rdat%aqx2          &
     &                        -rdat%fq3(3,i+4)+rdat%fq2(1,i+4)        )*xmdt
         rdat%r04( 7,i)= rdat%r04( 7,i)+(rdat%fq4(1,i)*rdat%acy2-rdat%fq3(1,i+4)* 3  )*xmdtxy
         rdat%r04( 8,i)= rdat%r04( 8,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i+4)     )*xmdtx
         rdat%r04( 9,i)= rdat%r04( 9,i)+(rdat%fq4(3,i)     -rdat%fq3(1,i+4)     )*xmdtxy
         rdat%r04(10,i)= rdat%r04(10,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i+4)* 3  )*xmdtx
         rdat%r04(11,i)= rdat%r04(11,i)+(rdat%fq4(1,i)*acy4-rdat%fq3(1,i+4)* 6 *rdat%acy2      &
     &                        +rdat%fq2(1,i+4)* 3                 )*xmdt
         rdat%r04(12,i)= rdat%r04(12,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i+4)* 3  )*xmdty
         rdat%r04(13,i)= rdat%r04(13,i)+(rdat%fq4(3,i)*rdat%acy2-rdat%fq3(1,i+4)*rdat%acy2          &
     &                        -rdat%fq3(3,i+4)+rdat%fq2(1,i+4)        )*xmdt
         rdat%r04(14,i)= rdat%r04(14,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i+4)* 3  )*xmdty
         rdat%r04(15,i)= rdat%r04(15,i)+(rdat%fq4(5,i)-rdat%fq3(3,i+4)* 6                 &
     &                        +rdat%fq2(1,i+4)* 3                 )*xmdt
      enddo
      fwk( 1,5)=  rdat%fq4(1,5)*aqx4-rdat%fq3(1,9)* 6 *rdat%aqx2+rdat%fq2(1,9)* 3
      fwk( 2,5)= (rdat%fq4(1,5)*rdat%aqx2-rdat%fq3(1,9)* 3                 )*rdat%aqxy
      fwk( 3,5)=-(rdat%fq4(2,5)*rdat%aqx2-rdat%fq3(2,9)* 3                 )*rdat%aqx
      fwk( 4,5)=  rdat%fq4(1,5)*x2y2-rdat%fq3(1,9)*q2c2    +rdat%fq2(1,9)
      fwk( 5,5)=-(rdat%fq4(2,5)*rdat%aqx2-rdat%fq3(2,9)                    )*rdat%acy
      fwk( 6,5)=  rdat%fq4(3,5)*rdat%aqx2-rdat%fq3(1,9)*rdat%aqx2-rdat%fq3(3,9)+rdat%fq2(1,9)
      fwk( 7,5)= (rdat%fq4(1,5)*rdat%acy2-rdat%fq3(1,9)* 3                 )*rdat%aqxy
      fwk( 8,5)=-(rdat%fq4(2,5)*rdat%acy2-rdat%fq3(2,9)                    )*rdat%aqx
      fwk( 9,5)= (rdat%fq4(3,5)     -rdat%fq3(1,9)                    )*rdat%aqxy
      fwk(10,5)=-(rdat%fq4(4,5)     -rdat%fq3(2,9)* 3                 )*rdat%aqx
      fwk(11,5)=  rdat%fq4(1,5)*acy4-rdat%fq3(1,9)* 6 *rdat%acy2+rdat%fq2(1,9)* 3
      fwk(12,5)=-(rdat%fq4(2,5)*rdat%acy2-rdat%fq3(2,9)* 3                 )*rdat%acy
      fwk(13,5)=  rdat%fq4(3,5)*rdat%acy2-rdat%fq3(1,9)*rdat%acy2-rdat%fq3(3,9)+rdat%fq2(1,9)
      fwk(14,5)=-(rdat%fq4(4,5)     -rdat%fq3(2,9)* 3                 )*rdat%acy
      fwk(15,5)=  rdat%fq4(5,5)     -rdat%fq3(3,9)* 6 +rdat%fq2(1,9)* 3
      do i= 5, 7
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,5)*work(i+ 5)
         end do
      end do

      rdat%r05( 1,1)= rdat%r05( 1,1)+(rdat%fq5(1,1)*aqx4-rdat%fq4(1,5)*10*rdat%aqx2                 &
     &                 +rdat%fq3(1,9)*15                         )*xmdtx
      rdat%r05( 2,1)= rdat%r05( 2,1)+(rdat%fq5(1,1)*aqx4-rdat%fq4(1,5)* 6 *rdat%aqx2                 &
     &                 +rdat%fq3(1,9)* 3                          )*xmdty
      rdat%r05( 3,1)= rdat%r05( 3,1)+(rdat%fq5(2,1)*aqx4-rdat%fq4(2,5)* 6 *rdat%aqx2                 &
     &                 +rdat%fq3(2,9)* 3                          )*xmdt
      rdat%r05( 4,1)= rdat%r05( 4,1)+(rdat%fq5(1,1)*x2y2-rdat%fq4(1,5)*rdat%aqx2-rdat%fq4(1,5)* 3 *rdat%acy2  &
     &                 +rdat%fq3(1,9)* 3                          )*xmdtx
      rdat%r05( 5,1)= rdat%r05( 5,1)+(rdat%fq5(2,1)*rdat%aqx2-rdat%fq4(2,5)* 3             )*xmdtxy
      rdat%r05( 6,1)= rdat%r05( 6,1)+(rdat%fq5(3,1)*rdat%aqx2-rdat%fq4(1,5)*rdat%aqx2-rdat%fq4(3,5)* 3        &
     &                 +rdat%fq3(1,9)* 3                          )*xmdtx
      rdat%r05( 7,1)= rdat%r05( 7,1)+(rdat%fq5(1,1)*x2y2-rdat%fq4(1,5)* 3 *rdat%aqx2-rdat%fq4(1,5)*rdat%acy2  &
     &                 +rdat%fq3(1,9)* 3                          )*xmdty
      rdat%r05( 8,1)= rdat%r05( 8,1)+(rdat%fq5(2,1)*x2y2-rdat%fq4(2,5)*q2c2+rdat%fq3(2,9) )*xmdt
      rdat%r05( 9,1)= rdat%r05( 9,1)+(rdat%fq5(3,1)*rdat%aqx2-rdat%fq4(1,5)*rdat%aqx2-rdat%fq4(3,5)           &
     &                 +rdat%fq3(1,9)                             )*xmdty
      rdat%r05(10,1)= rdat%r05(10,1)+(rdat%fq5(4,1)*rdat%aqx2-rdat%fq4(2,5)* 3 *rdat%aqx2-rdat%fq4(4,5)       &
     &                 +rdat%fq3(2,9)* 3                          )*xmdt
      rdat%r05(11,1)= rdat%r05(11,1)+(rdat%fq5(1,1)*acy4-rdat%fq4(1,5)* 6 *rdat%acy2                 &
     &                 +rdat%fq3(1,9)* 3                          )*xmdtx
      rdat%r05(12,1)= rdat%r05(12,1)+(rdat%fq5(2,1)*rdat%acy2-rdat%fq4(2,5)* 3             )*xmdtxy
      rdat%r05(13,1)= rdat%r05(13,1)+(rdat%fq5(3,1)*rdat%acy2-rdat%fq4(1,5)*rdat%acy2-rdat%fq4(3,5)           &
     &                 +rdat%fq3(1,9)                             )*xmdtx
      rdat%r05(14,1)= rdat%r05(14,1)+(rdat%fq5(4,1)-rdat%fq4(2,5)* 3                  )*xmdtxy
      rdat%r05(15,1)= rdat%r05(15,1)+(rdat%fq5(5,1)-rdat%fq4(3,5)* 6 +rdat%fq3(1,9)* 3    )*xmdtx
      rdat%r05(16,1)= rdat%r05(16,1)+(rdat%fq5(1,1)*acy4-rdat%fq4(1,5)*10*rdat%acy2                 &
     &                 +rdat%fq3(1,9)*15                         )*xmdty
      rdat%r05(17,1)= rdat%r05(17,1)+(rdat%fq5(2,1)*acy4-rdat%fq4(2,5)* 6 *rdat%acy2                 &
     &                 +rdat%fq3(2,9)* 3                          )*xmdt
      rdat%r05(18,1)= rdat%r05(18,1)+(rdat%fq5(3,1)*rdat%acy2-rdat%fq4(1,5)*rdat%acy2-rdat%fq4(3,5)* 3        &
     &                 +rdat%fq3(1,9)* 3                          )*xmdty
      rdat%r05(19,1)= rdat%r05(19,1)+(rdat%fq5(4,1)*rdat%acy2-rdat%fq4(2,5)* 3 *rdat%acy2-rdat%fq4(4,5)       &
     &                 +rdat%fq3(2,9)* 3                          )*xmdt
      rdat%r05(20,1)= rdat%r05(20,1)+(rdat%fq5(5,1)-rdat%fq4(3,5)* 6 +rdat%fq3(1,9)* 3    )*xmdty
      rdat%r05(21,1)= rdat%r05(21,1)+(rdat%fq5(6,1)-rdat%fq4(4,5)*10+rdat%fq3(2,9)*15   )*xmdt

      end subroutine intk_16

! >
! >    @brief   ddds case
! >
! >    @details integration of the ddds case
! >
      subroutine intk_17(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd2, xmd3, xmd4, xmd6
      real(kind=dp) :: xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: y33, y34, y44
      real(kind=dp) :: work, fwk, fw6, fcu, fcc
      integer :: i, j

! Generate jtype=17 integrals

      dimension  work(15),fwk(21,4),fw6(28)
      dimension  fcu(45,8),fcc(45,8)

      if(ikl == 0) then
         rdat%r00(:,1:2)= 0.0_dp
         rdat%r01(:,1:13)= 0.0_dp
         rdat%r02(:,1:17)= 0.0_dp
         rdat%r03(:,1:12)= 0.0_dp
         rdat%r04(:,1:8)= 0.0_dp
         rdat%r05(:,1:3)= 0.0_dp
         rdat%r06(1:28,1)= 0.0_dp

         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd3= xmd2
      xmd4= xmd3*xmd2
      xmd6= xmd4*xmd2
      xmdt= xmd6*xmd2
      xmd2= xmd3

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      y33 = rdat%y03 *rdat%y03
      y34 =-rdat%y03 *rdat%y04
      y44 = rdat%y04 *rdat%y04
      work( 1)= xmd4
      work( 2)= xmd2*y33
      work( 3)= xmd2*y34
      work( 4)= xmd2*y44
      work( 5)= y33 *y44

      work( 6)=-xmd4*rdat%y03
      work( 7)= xmd4*rdat%y04
      work( 8)= xmd2*y33*rdat%y04
      work( 9)= xmd2*y34*rdat%y04

      work(10)= xmd6
      work(11)= xmd4*y33
      work(12)= xmd4*y34
      work(13)= xmd4*y44

      work(14)=-xmd6*rdat%y03
      work(15)= xmd6*rdat%y04

      call fcufcc(rdat,6,xmdt,fcu,fcc)

      do i= 1, 2
         do j= 1, 5
            rdat%r00(j,i)= rdat%r00(j,i)+rdat%fq0(i)*work(j)
         end do
      end do

      do i= 1, 3
         fwk(1,i)=-rdat%fq1(1,i)*rdat%aqx
         fwk(2,i)=-rdat%fq1(1,i)*rdat%acy
         fwk(3,i)= rdat%fq1(2,i)
      enddo
      do i= 1, 4
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,1)*work(i+ 5)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,1)*work(i+ 5)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,1)*work(i+ 5)
      enddo
      do i= 5, 8
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,2)*work(i+ 1)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,2)*work(i+ 1)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,2)*work(i+ 1)
      enddo
      do i= 9,13
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,3)*work(i- 8)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,3)*work(i- 8)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,3)*work(i- 8)
      enddo

      do i= 1, 4
         fwk(1,i)= rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i)
         fwk(2,i)= rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i)
         fwk(3,i)= rdat%fq2(3,i)     -rdat%fq1(1,i)
         fwk(4,i)= rdat%fq2(1,i)*rdat%aqxy
         fwk(5,i)=-rdat%fq2(2,i)*rdat%aqx
         fwk(6,i)=-rdat%fq2(2,i)*rdat%acy
      enddo
      do i= 1, 4
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,1)*work(i+ 9)
         enddo
      enddo
      do i= 5, 8
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,2)*work(i+ 5)
         enddo
      enddo
      do i= 9,12
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,3)*work(i- 3)
         enddo
      enddo
      do i=13,17
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,4)*work(i-12)
         enddo
      enddo

      do i= 1, 4
         fwk( 1,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)* 3 )*rdat%aqx
         fwk( 2,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)    )*rdat%acy
         fwk( 3,i)=  rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i)
         fwk( 4,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)    )*rdat%aqx
         fwk( 5,i)=  rdat%fq3(2,i)                    *rdat%aqxy
         fwk( 6,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%aqx
         fwk( 7,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)* 3 )*rdat%acy
         fwk( 8,i)=  rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i)
         fwk( 9,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%acy
         fwk(10,i)=  rdat%fq3(4,i)     -rdat%fq2(2,i)* 3
      enddo
      do i= 1, 2
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,1)*work(i+13)
         enddo
      enddo
      do i= 3, 4
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,2)*work(i+11)
         enddo
      enddo
      do i= 5, 8
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,3)*work(i+ 5)
         enddo
      enddo
      do i= 9,12
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,4)*work(i- 3)
         enddo
      enddo

      aqx4= rdat%aqx2*rdat%aqx2
      acy4= rdat%acy2*rdat%acy2
      x2y2= rdat%aqx2*rdat%acy2
      q2c2= rdat%aqx2+rdat%acy2
      do i= 1, 2
         rdat%r04( 1,i)= rdat%r04( 1,i)+(rdat%fq4(1,i)*aqx4-rdat%fq3(1,i)* 6 *rdat%aqx2        &
     &                        +rdat%fq2(1,i)* 3                 )*xmdt
         rdat%r04( 2,i)= rdat%r04( 2,i)+(rdat%fq4(1,i)*rdat%aqx2-rdat%fq3(1,i)* 3  )*xmdtxy
         rdat%r04( 3,i)= rdat%r04( 3,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)* 3  )*xmdtx
         rdat%r04( 4,i)= rdat%r04( 4,i)+(rdat%fq4(1,i)*x2y2-rdat%fq3(1,i)*q2c2            &
     &                                       +rdat%fq2(1,i)     )*xmdt
         rdat%r04( 5,i)= rdat%r04( 5,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)     )*xmdty
         rdat%r04( 6,i)= rdat%r04( 6,i)+(rdat%fq4(3,i)*rdat%aqx2-rdat%fq3(1,i)*rdat%aqx2            &
     &                        -rdat%fq3(3,i)+rdat%fq2(1,i)          )*xmdt
         rdat%r04( 7,i)= rdat%r04( 7,i)+(rdat%fq4(1,i)*rdat%acy2-rdat%fq3(1,i)* 3  )*xmdtxy
         rdat%r04( 8,i)= rdat%r04( 8,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)     )*xmdtx
         rdat%r04( 9,i)= rdat%r04( 9,i)+(rdat%fq4(3,i)     -rdat%fq3(1,i)     )*xmdtxy
         rdat%r04(10,i)= rdat%r04(10,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3  )*xmdtx
         rdat%r04(11,i)= rdat%r04(11,i)+(rdat%fq4(1,i)*acy4-rdat%fq3(1,i)* 6 *rdat%acy2        &
     &                        +rdat%fq2(1,i)* 3                 )*xmdt
         rdat%r04(12,i)= rdat%r04(12,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)* 3  )*xmdty
         rdat%r04(13,i)= rdat%r04(13,i)+(rdat%fq4(3,i)*rdat%acy2-rdat%fq3(1,i)*rdat%acy2            &
     &                        -rdat%fq3(3,i)+rdat%fq2(1,i)          )*xmdt
         rdat%r04(14,i)= rdat%r04(14,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3  )*xmdty
         rdat%r04(15,i)= rdat%r04(15,i)+(rdat%fq4(5,i)-rdat%fq3(3,i)* 6                   &
     &                        +rdat%fq2(1,i)* 3                 )*xmdt
      enddo
      do i= 3, 4
         fwk( 1,i)=  rdat%fq4(1,i)*aqx4-rdat%fq3(1,i)* 6 *rdat%aqx2+rdat%fq2(1,i)* 3
         fwk( 2,i)= (rdat%fq4(1,i)*rdat%aqx2-rdat%fq3(1,i)* 3            )*rdat%aqxy
         fwk( 3,i)=-(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)* 3            )*rdat%aqx
         fwk( 4,i)=  rdat%fq4(1,i)*x2y2-rdat%fq3(1,i)*q2c2    +rdat%fq2(1,i)
         fwk( 5,i)=-(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)               )*rdat%acy
         fwk( 6,i)=  rdat%fq4(3,i)*rdat%aqx2-rdat%fq3(1,i)*rdat%aqx2-rdat%fq3(3,i)+rdat%fq2(1,i)
         fwk( 7,i)= (rdat%fq4(1,i)*rdat%acy2-rdat%fq3(1,i)* 3            )*rdat%aqxy
         fwk( 8,i)=-(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)               )*rdat%aqx
         fwk( 9,i)= (rdat%fq4(3,i)     -rdat%fq3(1,i)               )*rdat%aqxy
         fwk(10,i)=-(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3            )*rdat%aqx
         fwk(11,i)=  rdat%fq4(1,i)*acy4-rdat%fq3(1,i)* 6 *rdat%acy2+rdat%fq2(1,i)* 3
         fwk(12,i)=-(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)* 3            )*rdat%acy
         fwk(13,i)=  rdat%fq4(3,i)*rdat%acy2-rdat%fq3(1,i)*rdat%acy2-rdat%fq3(3,i)+rdat%fq2(1,i)
         fwk(14,i)=-(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3            )*rdat%acy
         fwk(15,i)=  rdat%fq4(5,i)     -rdat%fq3(3,i)* 6 +rdat%fq2(1,i)* 3
      enddo
      do i= 3, 4
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,3)*work(i+11)
         enddo
      enddo
      do i= 5, 8
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,4)*work(i+ 5)
         enddo
      enddo

      rdat%r05( 1,1)= rdat%r05( 1,1)+(rdat%fq5(1,1)*aqx4-rdat%fq4(1,3)*10*rdat%aqx2           &
     &                     +rdat%fq3(1,3)*15                    )*xmdtx
      rdat%r05( 2,1)= rdat%r05( 2,1)+(rdat%fq5(1,1)*aqx4-rdat%fq4(1,3)* 6 *rdat%aqx2           &
     &                     +rdat%fq3(1,3)* 3                     )*xmdty
      rdat%r05( 3,1)= rdat%r05( 3,1)+(rdat%fq5(2,1)*aqx4-rdat%fq4(2,3)* 6 *rdat%aqx2           &
     &                     +rdat%fq3(2,3)* 3                     )*xmdt
      rdat%r05( 4,1)= rdat%r05( 4,1)+(rdat%fq5(1,1)*x2y2-rdat%fq4(1,3)*rdat%aqx2               &
     &                     -rdat%fq4(1,3)* 3 *rdat%acy2+rdat%fq3(1,3)* 3  )*xmdtx
      rdat%r05( 5,1)= rdat%r05( 5,1)+(rdat%fq5(2,1)*rdat%aqx2-rdat%fq4(2,3)* 3      )*xmdtxy
      rdat%r05( 6,1)= rdat%r05( 6,1)+(rdat%fq5(3,1)*rdat%aqx2-rdat%fq4(1,3)*rdat%aqx2               &
     &                     -rdat%fq4(3,3)* 3 +rdat%fq3(1,3)* 3       )*xmdtx
      rdat%r05( 7,1)= rdat%r05( 7,1)+(rdat%fq5(1,1)*x2y2-rdat%fq4(1,3)* 3 *rdat%aqx2           &
     &                     -rdat%fq4(1,3)*rdat%acy2+rdat%fq3(1,3)* 3      )*xmdty
      rdat%r05( 8,1)= rdat%r05( 8,1)+(rdat%fq5(2,1)*x2y2-rdat%fq4(2,3)*q2c2               &
     &                     +rdat%fq3(2,3)                        )*xmdt
      rdat%r05( 9,1)= rdat%r05( 9,1)+(rdat%fq5(3,1)*rdat%aqx2-rdat%fq4(1,3)*rdat%aqx2-rdat%fq4(3,3)     &
     &                     +rdat%fq3(1,3)                        )*xmdty
      rdat%r05(10,1)= rdat%r05(10,1)+(rdat%fq5(4,1)*rdat%aqx2-rdat%fq4(2,3)* 3 *rdat%aqx2           &
     &                     -rdat%fq4(4,3)+rdat%fq3(2,3)* 3           )*xmdt
      rdat%r05(11,1)= rdat%r05(11,1)+(rdat%fq5(1,1)*acy4-rdat%fq4(1,3)* 6 *rdat%acy2           &
     &                     +rdat%fq3(1,3)* 3                     )*xmdtx
      rdat%r05(12,1)= rdat%r05(12,1)+(rdat%fq5(2,1)*rdat%acy2-rdat%fq4(2,3)* 3      )*xmdtxy
      rdat%r05(13,1)= rdat%r05(13,1)+(rdat%fq5(3,1)*rdat%acy2-rdat%fq4(1,3)*rdat%acy2-rdat%fq4(3,3)     &
     &                     +rdat%fq3(1,3)                        )*xmdtx
      rdat%r05(14,1)= rdat%r05(14,1)+(rdat%fq5(4,1)-rdat%fq4(2,3)* 3           )*xmdtxy
      rdat%r05(15,1)= rdat%r05(15,1)+(rdat%fq5(5,1)-rdat%fq4(3,3)* 6                      &
     &                     +rdat%fq3(1,3)* 3                     )*xmdtx
      rdat%r05(16,1)= rdat%r05(16,1)+(rdat%fq5(1,1)*acy4-rdat%fq4(1,3)*10*rdat%acy2           &
     &                     +rdat%fq3(1,3)*15                    )*xmdty
      rdat%r05(17,1)= rdat%r05(17,1)+(rdat%fq5(2,1)*acy4-rdat%fq4(2,3)* 6 *rdat%acy2           &
     &                     +rdat%fq3(2,3)* 3                     )*xmdt
      rdat%r05(18,1)= rdat%r05(18,1)+(rdat%fq5(3,1)*rdat%acy2-rdat%fq4(1,3)*rdat%acy2               &
     &                     -rdat%fq4(3,3)* 3 +rdat%fq3(1,3)* 3       )*xmdty
      rdat%r05(19,1)= rdat%r05(19,1)+(rdat%fq5(4,1)*rdat%acy2-rdat%fq4(2,3)* 3 *rdat%acy2           &
     &                     -rdat%fq4(4,3)+rdat%fq3(2,3)* 3           )*xmdt
      rdat%r05(20,1)= rdat%r05(20,1)+(rdat%fq5(5,1)-rdat%fq4(3,3)* 6                      &
     &                     +rdat%fq3(1,3)* 3                     )*xmdty
      rdat%r05(21,1)= rdat%r05(21,1)+(rdat%fq5(6,1)-rdat%fq4(4,3)*10                     &
     &                     +rdat%fq3(2,3)*15                    )*xmdt

      fwk( 1,2)=-(rdat%fq5(1,2)*aqx4-rdat%fq4(1,4)*10*rdat%aqx2+rdat%fq3(1,4)*15)*rdat%aqx
      fwk( 2,2)=-(rdat%fq5(1,2)*aqx4-rdat%fq4(1,4)* 6 *rdat%aqx2+rdat%fq3(1,4)* 3 )*rdat%acy
      fwk( 3,2)=  rdat%fq5(2,2)*aqx4-rdat%fq4(2,4)* 6 *rdat%aqx2+rdat%fq3(2,4)* 3
      fwk( 4,2)=-(rdat%fq5(1,2)*x2y2-rdat%fq4(1,4)*rdat%aqx2-rdat%fq4(1,4)* 3 *rdat%acy2      &
     &                          +rdat%fq3(1,4)* 3                    )*rdat%aqx
      fwk( 5,2)= (rdat%fq5(2,2)*rdat%aqx2-rdat%fq4(2,4)* 3                   )*rdat%aqxy
      fwk( 6,2)=-(rdat%fq5(3,2)*rdat%aqx2-rdat%fq4(1,4)*rdat%aqx2-rdat%fq4(3,4)* 3            &
     &                          +rdat%fq3(1,4)* 3                    )*rdat%aqx
      fwk( 7,2)=-(rdat%fq5(1,2)*x2y2-rdat%fq4(1,4)* 3 *rdat%aqx2-rdat%fq4(1,4)*rdat%acy2      &
     &                          +rdat%fq3(1,4)* 3                    )*rdat%acy
      fwk( 8,2)=  rdat%fq5(2,2)*x2y2-rdat%fq4(2,4)*q2c2+rdat%fq3(2,4)
      fwk( 9,2)=-(rdat%fq5(3,2)*rdat%aqx2-rdat%fq4(1,4)*rdat%aqx2-rdat%fq4(3,4)               &
     &                          +rdat%fq3(1,4)                       )*rdat%acy
      fwk(10,2)=  rdat%fq5(4,2)*rdat%aqx2-rdat%fq4(2,4)* 3 *rdat%aqx2-rdat%fq4(4,4)           &
     &                          +rdat%fq3(2,4)* 3
      fwk(11,2)=-(rdat%fq5(1,2)*acy4-rdat%fq4(1,4)* 6 *rdat%acy2                     &
     &                          +rdat%fq3(1,4)* 3                    )*rdat%aqx
      fwk(12,2)= (rdat%fq5(2,2)*rdat%acy2-rdat%fq4(2,4)* 3                   )*rdat%aqxy
      fwk(13,2)=-(rdat%fq5(3,2)*rdat%acy2-rdat%fq4(1,4)*rdat%acy2-rdat%fq4(3,4)               &
     &                          +rdat%fq3(1,4)                       )*rdat%aqx
      fwk(14,2)= (rdat%fq5(4,2)     -rdat%fq4(2,4)* 3                   )*rdat%aqxy
      fwk(15,2)=-(rdat%fq5(5,2)     -rdat%fq4(3,4)* 6 +rdat%fq3(1,4)* 3      )*rdat%aqx
      fwk(16,2)=-(rdat%fq5(1,2)*acy4-rdat%fq4(1,4)*10*rdat%acy2+rdat%fq3(1,4)*15)*rdat%acy
      fwk(17,2)=  rdat%fq5(2,2)*acy4-rdat%fq4(2,4)* 6 *rdat%acy2+rdat%fq3(2,4)* 3
      fwk(18,2)=-(rdat%fq5(3,2)*rdat%acy2-rdat%fq4(1,4)*rdat%acy2-rdat%fq4(3,4)* 3            &
     &                          +rdat%fq3(1,4)* 3                    )*rdat%acy
      fwk(19,2)=  rdat%fq5(4,2)*rdat%acy2-rdat%fq4(2,4)* 3 *rdat%acy2-rdat%fq4(4,4)           &
     &                          +rdat%fq3(2,4)* 3
      fwk(20,2)=-(rdat%fq5(5,2)     -rdat%fq4(3,4)* 6 +rdat%fq3(1,4)* 3      )*rdat%acy
      fwk(21,2)=  rdat%fq5(6,2)     -rdat%fq4(4,4)*10+rdat%fq3(2,4)*15
      do i= 2, 3
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,2)*work(i+12)
         end do
      end do

      call frikr6(rdat, 1, 1,fw6,rdat%fq6, 1,rdat%fq5, 3,rdat%fq4, 3,rdat%fq3)

      do j= 1,28
         rdat%r06(j,1)= rdat%r06(j,1)+fw6(j)*fcc(j,6)
      enddo

      end subroutine intk_17

! >
! >    @brief   ddpp case
! >
! >    @details integration of the ddpp case
! >
      subroutine intk_18(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd2, xmd3, xmd4, xmd6
      real(kind=dp) :: xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: fqd11, fqd12, fqd21, fqd22, fqd23, &
                       fqd31, fqd32, fqd33, fqd34
      real(kind=dp) :: y33, y34, y44
      real(kind=dp) :: work, fwk, fw6, fcu, fcc
      integer :: i, j

! Generate jtype=18 integrals

      dimension  work(15),fwk(21,10),fw6(28)
      dimension  fcu(45,8),fcc(45,8)

      if(ikl == 0) then
         rdat%r00= 0.0_dp
         rdat%r01= 0.0_dp
         rdat%r02(:,1:41)= 0.0_dp
         rdat%r03(:,1:30)= 0.0_dp
         rdat%r04(:,1:17)= 0.0_dp
         rdat%r05(:,1:6)= 0.0_dp
         rdat%r06(1:28,1)= 0.0_dp
         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd3= xmd2
      xmd4= xmd3*xmd2
      xmd6= xmd4*xmd2
      xmdt= xmd6*xmd2
      xmd2= xmd3

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      y33 = rdat%y03 *rdat%y03
      y34 =-rdat%y03 *rdat%y04
      y44 = rdat%y04 *rdat%y04
      work( 1)= xmd4
      work( 2)= xmd2*y33
      work( 3)= xmd2*y34
      work( 4)= xmd2*y44
      work( 5)= y33 *y44

      work( 6)=-xmd4*rdat%y03
      work( 7)= xmd4*rdat%y04
      work( 8)= xmd2*y33*rdat%y04
      work( 9)= xmd2*y34*rdat%y04

      work(10)= xmd6
      work(11)= xmd4*y33
      work(12)= xmd4*y34
      work(13)= xmd4*y44

      work(14)=-xmd6*rdat%y03
      work(15)= xmd6*rdat%y04

      call fcufcc(rdat,6,xmdt,fcu,fcc)

      do i= 1, 5
         do j= 1, 5
            rdat%r00(j,i)= rdat%r00(j,i)+rdat%fq0(i)*work(j)
         end do
      end do

      do i= 1, 8
         fwk(1,i)=-rdat%fq1(1,i)*rdat%aqx
         fwk(2,i)=-rdat%fq1(1,i)*rdat%acy
         fwk(3,i)= rdat%fq1(2,i)
      enddo
      fqd11 = rdat%fq1(1,5)*rdat%rab +rdat%fq1(1,8)
      fqd12 = rdat%fq1(2,5)*rdat%rab +rdat%fq1(2,8)
      fwk(1,9)=-fqd11*rdat%aqx
      fwk(2,9)=-fqd11*rdat%acy
      fwk(3,9)= fqd12
      do i= 1, 4
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,1)*work(i+ 5)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,1)*work(i+ 5)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,1)*work(i+ 5)
      enddo
      do i= 5, 8
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,2)*work(i+ 1)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,2)*work(i+ 1)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,2)*work(i+ 1)
      enddo
      do i= 9,12
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,3)*work(i- 3)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,3)*work(i- 3)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,3)*work(i- 3)
      enddo
      do i=13,16
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,4)*work(i- 7)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,4)*work(i- 7)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,4)*work(i- 7)
      enddo
      do i=17,20
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,5)*work(i-11)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,5)*work(i-11)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,5)*work(i-11)
      enddo
      do i=21,25
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,6)*work(i-20)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,6)*work(i-20)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,6)*work(i-20)
      enddo
      do i=26,30
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,7)*work(i-25)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,7)*work(i-25)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,7)*work(i-25)
      enddo
      do i=31,35
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,8)*work(i-30)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,8)*work(i-30)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,8)*work(i-30)
      enddo
      do i=36,40
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,9)*work(i-35)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,9)*work(i-35)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,9)*work(i-35)
      enddo

      do i= 1, 9
         fwk(1,i)= rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i)
         fwk(2,i)= rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i)
         fwk(3,i)= rdat%fq2(3,i)     -rdat%fq1(1,i)
         fwk(4,i)= rdat%fq2(1,i)*rdat%aqxy
         fwk(5,i)=-rdat%fq2(2,i)*rdat%aqx
         fwk(6,i)=-rdat%fq2(2,i)*rdat%acy
      enddo
      fqd21 = rdat%fq2(1,5)*rdat%rab +rdat%fq2(1,8)
      fqd22 = rdat%fq2(2,5)*rdat%rab +rdat%fq2(2,8)
      fqd23 = rdat%fq2(3,5)*rdat%rab +rdat%fq2(3,8)
         fwk(1,10)= fqd21*rdat%aqx2-fqd11
         fwk(2,10)= fqd21*rdat%acy2-fqd11
         fwk(3,10)= fqd23     -fqd11
         fwk(4,10)= fqd21*rdat%aqxy
         fwk(5,10)=-fqd22*rdat%aqx
         fwk(6,10)=-fqd22*rdat%acy
      do i= 1, 4
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,1)*work(i+ 9)
         enddo
      enddo
      do i= 5, 8
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,2)*work(i+ 5)
         enddo
      enddo
      do i= 9,12
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,3)*work(i+ 1)
         enddo
      enddo
      do i=13,16
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,4)*work(i- 3)
         enddo
      enddo
      do i=17,20
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,5)*work(i- 7)
         enddo
      enddo
      do i=21,24
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,6)*work(i-15)
         enddo
      enddo
      do i=25,28
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,7)*work(i-19)
         enddo
      enddo
      do i=29,32
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,8)*work(i-23)
         enddo
      enddo
      do i=33,36
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,10)*work(i-27)
         enddo
      enddo
      do i=37,41
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,9)*work(i-36)
         enddo
      enddo

      do i= 1, 9
         fwk( 1,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)* 3 )*rdat%aqx
         fwk( 2,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)    )*rdat%acy
         fwk( 3,i)=  rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i)
         fwk( 4,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)    )*rdat%aqx
         fwk( 5,i)=  rdat%fq3(2,i)                    *rdat%aqxy
         fwk( 6,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%aqx
         fwk( 7,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)* 3 )*rdat%acy
         fwk( 8,i)=  rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i)
         fwk( 9,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%acy
         fwk(10,i)=  rdat%fq3(4,i)     -rdat%fq2(2,i)* 3
      enddo
      fqd31 = rdat%fq3(1,5)*rdat%rab +rdat%fq3(1,8)
      fqd32 = rdat%fq3(2,5)*rdat%rab +rdat%fq3(2,8)
      fqd33 = rdat%fq3(3,5)*rdat%rab +rdat%fq3(3,8)
      fqd34 = rdat%fq3(4,5)*rdat%rab +rdat%fq3(4,8)
      fwk( 1,10)=-(fqd31*rdat%aqx2-fqd21* 3 )*rdat%aqx
      fwk( 2,10)=-(fqd31*rdat%aqx2-fqd21    )*rdat%acy
      fwk( 3,10)=  fqd32*rdat%aqx2-fqd22
      fwk( 4,10)=-(fqd31*rdat%acy2-fqd21    )*rdat%aqx
      fwk( 5,10)=  fqd32                *rdat%aqxy
      fwk( 6,10)=-(fqd33     -fqd21    )*rdat%aqx
      fwk( 7,10)=-(fqd31*rdat%acy2-fqd21* 3 )*rdat%acy
      fwk( 8,10)=  fqd32*rdat%acy2-fqd22
      fwk( 9,10)=-(fqd33     -fqd21    )*rdat%acy
      fwk(10,10)=  fqd34     -fqd22* 3
      do i= 1, 2
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,1)*work(i+13)
         enddo
      enddo
      do i= 3, 4
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,2)*work(i+11)
         enddo
      enddo
      do i= 5, 6
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,3)*work(i+ 9)
         enddo
      enddo
      do i= 7, 8
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,4)*work(i+ 7)
         enddo
      enddo
      do i= 9,10
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,5)*work(i+ 5)
         enddo
      enddo
      do i=11,14
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,6)*work(i- 1)
         enddo
      enddo
      do i=15,18
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,7)*work(i- 5)
         enddo
      enddo
      do i=19,22
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,8)*work(i- 9)
         enddo
      enddo
      do i=23,26
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,10)*work(i-13)
         enddo
      enddo
      do i=27,30
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,9)*work(i-21)
         enddo
      enddo

      aqx4= rdat%aqx2*rdat%aqx2
      acy4= rdat%acy2*rdat%acy2
      x2y2= rdat%aqx2*rdat%acy2
      q2c2= rdat%aqx2+rdat%acy2
      do i= 1, 5
         rdat%r04( 1,i)= rdat%r04( 1,i)+(rdat%fq4(1,i)*aqx4-rdat%fq3(1,i)* 6 *rdat%aqx2        &
     &                                       +rdat%fq2(1,i)* 3   )*xmdt
         rdat%r04( 2,i)= rdat%r04( 2,i)+(rdat%fq4(1,i)*rdat%aqx2-rdat%fq3(1,i)* 3   )*xmdtxy
         rdat%r04( 3,i)= rdat%r04( 3,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)* 3   )*xmdtx
         rdat%r04( 4,i)= rdat%r04( 4,i)+(rdat%fq4(1,i)*x2y2-rdat%fq3(1,i)*q2c2            &
     &                                       +rdat%fq2(1,i)      )*xmdt
         rdat%r04( 5,i)= rdat%r04( 5,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)      )*xmdty
         rdat%r04( 6,i)= rdat%r04( 6,i)+(rdat%fq4(3,i)*rdat%aqx2-rdat%fq3(1,i)*rdat%aqx2            &
     &                             -rdat%fq3(3,i)+rdat%fq2(1,i)      )*xmdt
         rdat%r04( 7,i)= rdat%r04( 7,i)+(rdat%fq4(1,i)*rdat%acy2-rdat%fq3(1,i)* 3   )*xmdtxy
         rdat%r04( 8,i)= rdat%r04( 8,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)      )*xmdtx
         rdat%r04( 9,i)= rdat%r04( 9,i)+(rdat%fq4(3,i)     -rdat%fq3(1,i)      )*xmdtxy
         rdat%r04(10,i)= rdat%r04(10,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3   )*xmdtx
         rdat%r04(11,i)= rdat%r04(11,i)+(rdat%fq4(1,i)*acy4-rdat%fq3(1,i)* 6 *rdat%acy2        &
     &                                       +rdat%fq2(1,i)* 3   )*xmdt
         rdat%r04(12,i)= rdat%r04(12,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)* 3   )*xmdty
         rdat%r04(13,i)= rdat%r04(13,i)+(rdat%fq4(3,i)*rdat%acy2-rdat%fq3(1,i)*rdat%acy2            &
     &                             -rdat%fq3(3,i)+rdat%fq2(1,i)      )*xmdt
         rdat%r04(14,i)= rdat%r04(14,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3   )*xmdty
         rdat%r04(15,i)= rdat%r04(15,i)+(rdat%fq4(5,i)-rdat%fq3(3,i)* 6                   &
     &                                       +rdat%fq2(1,i)* 3   )*xmdt
      enddo
      rdat%fq2(1,5)= fqd21
      rdat%fq3(1,5)= fqd31
      rdat%fq3(2,5)= fqd32
      rdat%fq3(3,5)= fqd33
      do j= 1, 5
         rdat%fq4(j,5)= rdat%fq4(j,5)*rdat%rab +rdat%fq4(j,8)
      enddo

      do i= 5, 9
         fwk( 1,i)=  rdat%fq4(1,i)*aqx4-rdat%fq3(1,i)* 6 *rdat%aqx2+rdat%fq2(1,i)* 3
         fwk( 2,i)= (rdat%fq4(1,i)*rdat%aqx2-rdat%fq3(1,i)* 3                )*rdat%aqxy
         fwk( 3,i)=-(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)* 3                )*rdat%aqx
         fwk( 4,i)=  rdat%fq4(1,i)*x2y2-rdat%fq3(1,i)*q2c2    +rdat%fq2(1,i)
         fwk( 5,i)=-(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)                   )*rdat%acy
         fwk( 6,i)=  rdat%fq4(3,i)*rdat%aqx2-rdat%fq3(1,i)*rdat%aqx2-rdat%fq3(3,i)+rdat%fq2(1,i)
         fwk( 7,i)= (rdat%fq4(1,i)*rdat%acy2-rdat%fq3(1,i)* 3                )*rdat%aqxy
         fwk( 8,i)=-(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)                   )*rdat%aqx
         fwk( 9,i)= (rdat%fq4(3,i)     -rdat%fq3(1,i)                   )*rdat%aqxy
         fwk(10,i)=-(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3                )*rdat%aqx
         fwk(11,i)=  rdat%fq4(1,i)*acy4-rdat%fq3(1,i)* 6 *rdat%acy2+rdat%fq2(1,i)* 3
         fwk(12,i)=-(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)* 3                )*rdat%acy
         fwk(13,i)=  rdat%fq4(3,i)*rdat%acy2-rdat%fq3(1,i)*rdat%acy2-rdat%fq3(3,i)+rdat%fq2(1,i)
         fwk(14,i)=-(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3                )*rdat%acy
         fwk(15,i)=  rdat%fq4(5,i)     -rdat%fq3(3,i)* 6 +rdat%fq2(1,i)* 3
      enddo
      do i= 6, 7
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,5)*work(i+ 8)
         enddo
      enddo
      do i= 8, 9
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,6)*work(i+ 6)
         enddo
      enddo
      do i=10,11
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,7)*work(i+ 4)
         enddo
      enddo
      do i=12,13
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,8)*work(i+ 2)
         enddo
      enddo
      do i=14,17
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,9)*work(i- 4)
         enddo
      enddo

      do j= 1, 6
         rdat%fq5(j,1)= rdat%fq5(j,1)*rdat%rab +rdat%fq5(j,4)
      enddo
      do i= 1, 4
        rdat%r05( 1,i)= rdat%r05( 1,i)+(rdat%fq5(1,i)*aqx4-rdat%fq4(1,i+4)*10*rdat%aqx2       &
     &                       +rdat%fq3(1,i+4)*15                  )*xmdtx
        rdat%r05( 2,i)= rdat%r05( 2,i)+(rdat%fq5(1,i)*aqx4-rdat%fq4(1,i+4)* 6 *rdat%aqx2       &
     &                       +rdat%fq3(1,i+4)* 3                   )*xmdty
        rdat%r05( 3,i)= rdat%r05( 3,i)+(rdat%fq5(2,i)*aqx4-rdat%fq4(2,i+4)* 6 *rdat%aqx2       &
     &                       +rdat%fq3(2,i+4)* 3                   )*xmdt
        rdat%r05( 4,i)= rdat%r05( 4,i)+(rdat%fq5(1,i)*x2y2-rdat%fq4(1,i+4)* 3 *rdat%acy2       &
     &                       -rdat%fq4(1,i+4)*rdat%aqx2+rdat%fq3(1,i+4)* 3  )*xmdtx
        rdat%r05( 5,i)= rdat%r05( 5,i)+(rdat%fq5(2,i)*rdat%aqx2-rdat%fq4(2,i+4)* 3    )*xmdtxy
        rdat%r05( 6,i)= rdat%r05( 6,i)+(rdat%fq5(3,i)*rdat%aqx2-rdat%fq4(1,i+4)*rdat%aqx2           &
     &                       -rdat%fq4(3,i+4)* 3 +rdat%fq3(1,i+4)* 3   )*xmdtx
        rdat%r05( 7,i)= rdat%r05( 7,i)+(rdat%fq5(1,i)*x2y2-rdat%fq4(1,i+4)* 3 *rdat%aqx2       &
     &                       -rdat%fq4(1,i+4)*rdat%acy2+rdat%fq3(1,i+4)* 3  )*xmdty
        rdat%r05( 8,i)= rdat%r05( 8,i)+(rdat%fq5(2,i)*x2y2-rdat%fq4(2,i+4)*q2c2           &
     &                       +rdat%fq3(2,i+4)                      )*xmdt
        rdat%r05( 9,i)= rdat%r05( 9,i)+(rdat%fq5(3,i)*rdat%aqx2-rdat%fq4(1,i+4)*rdat%aqx2           &
     &                       -rdat%fq4(3,i+4)+rdat%fq3(1,i+4)          )*xmdty
        rdat%r05(10,i)= rdat%r05(10,i)+(rdat%fq5(4,i)*rdat%aqx2-rdat%fq4(2,i+4)* 3 *rdat%aqx2       &
     &                       -rdat%fq4(4,i+4)+rdat%fq3(2,i+4)* 3       )*xmdt
        rdat%r05(11,i)= rdat%r05(11,i)+(rdat%fq5(1,i)*acy4-rdat%fq4(1,i+4)* 6 *rdat%acy2       &
     &                       +rdat%fq3(1,i+4)* 3                   )*xmdtx
        rdat%r05(12,i)= rdat%r05(12,i)+(rdat%fq5(2,i)*rdat%acy2-rdat%fq4(2,i+4)* 3    )*xmdtxy
        rdat%r05(13,i)= rdat%r05(13,i)+(rdat%fq5(3,i)*rdat%acy2-rdat%fq4(1,i+4)*rdat%acy2           &
     &                       -rdat%fq4(3,i+4)+rdat%fq3(1,i+4)          )*xmdtx
        rdat%r05(14,i)= rdat%r05(14,i)+(rdat%fq5(4,i)-rdat%fq4(2,i+4)* 3         )*xmdtxy
        rdat%r05(15,i)= rdat%r05(15,i)+(rdat%fq5(5,i)-rdat%fq4(3,i+4)* 6                  &
     &                       +rdat%fq3(1,i+4)* 3                   )*xmdtx
        rdat%r05(16,i)= rdat%r05(16,i)+(rdat%fq5(1,i)*acy4-rdat%fq4(1,i+4)*10*rdat%acy2       &
     &                       +rdat%fq3(1,i+4)*15                  )*xmdty
        rdat%r05(17,i)= rdat%r05(17,i)+(rdat%fq5(2,i)*acy4-rdat%fq4(2,i+4)* 6 *rdat%acy2       &
     &                       +rdat%fq3(2,i+4)* 3                   )*xmdt
        rdat%r05(18,i)= rdat%r05(18,i)+(rdat%fq5(3,i)*rdat%acy2-rdat%fq4(1,i+4)*rdat%acy2           &
     &                       -rdat%fq4(3,i+4)* 3 +rdat%fq3(1,i+4)* 3   )*xmdty
        rdat%r05(19,i)= rdat%r05(19,i)+(rdat%fq5(4,i)*rdat%acy2-rdat%fq4(2,i+4)* 3 *rdat%acy2       &
     &                       -rdat%fq4(4,i+4)+rdat%fq3(2,i+4)* 3       )*xmdt
        rdat%r05(20,i)= rdat%r05(20,i)+(rdat%fq5(5,i)-rdat%fq4(3,i+4)* 6                  &
     &                       +rdat%fq3(1,i+4)* 3                   )*xmdty
        rdat%r05(21,i)= rdat%r05(21,i)+(rdat%fq5(6,i)-rdat%fq4(4,i+4)*10                 &
     &                       +rdat%fq3(2,i+4)*15                  )*xmdt
      enddo
      fwk( 1,5)=-(rdat%fq5(1,5)*aqx4-rdat%fq4(1,9)*10*rdat%aqx2+rdat%fq3(1,9)*15)*rdat%aqx
      fwk( 2,5)=-(rdat%fq5(1,5)*aqx4-rdat%fq4(1,9)* 6 *rdat%aqx2+rdat%fq3(1,9)* 3 )*rdat%acy
      fwk( 3,5)=  rdat%fq5(2,5)*aqx4-rdat%fq4(2,9)* 6 *rdat%aqx2+rdat%fq3(2,9)* 3
      fwk( 4,5)=-(rdat%fq5(1,5)*x2y2-rdat%fq4(1,9)*rdat%aqx2-rdat%fq4(1,9)* 3 *rdat%acy2      &
     &                          +rdat%fq3(1,9)* 3                    )*rdat%aqx
      fwk( 5,5)= (rdat%fq5(2,5)*rdat%aqx2-rdat%fq4(2,9)* 3                    )*rdat%aqxy
      fwk( 6,5)=-(rdat%fq5(3,5)*rdat%aqx2-rdat%fq4(1,9)*rdat%aqx2-rdat%fq4(3,9)* 3            &
     &                          +rdat%fq3(1,9)* 3                    )*rdat%aqx
      fwk( 7,5)=-(rdat%fq5(1,5)*x2y2-rdat%fq4(1,9)* 3 *rdat%aqx2-rdat%fq4(1,9)*rdat%acy2      &
     &                          +rdat%fq3(1,9)* 3                    )*rdat%acy
      fwk( 8,5)=  rdat%fq5(2,5)*x2y2-rdat%fq4(2,9)*q2c2+rdat%fq3(2,9)
      fwk( 9,5)=-(rdat%fq5(3,5)*rdat%aqx2-rdat%fq4(1,9)*rdat%aqx2-rdat%fq4(3,9)               &
     &                          +rdat%fq3(1,9)                       )*rdat%acy
      fwk(10,5)=  rdat%fq5(4,5)*rdat%aqx2-rdat%fq4(2,9)* 3 *rdat%aqx2-rdat%fq4(4,9)           &
     &                          +rdat%fq3(2,9)* 3
      fwk(11,5)=-(rdat%fq5(1,5)*acy4-rdat%fq4(1,9)* 6 *rdat%acy2+rdat%fq3(1,9)* 3 )*rdat%aqx
      fwk(12,5)= (rdat%fq5(2,5)*rdat%acy2-rdat%fq4(2,9)* 3                    )*rdat%aqxy
      fwk(13,5)=-(rdat%fq5(3,5)*rdat%acy2-rdat%fq4(1,9)*rdat%acy2-rdat%fq4(3,9)               &
     &                          +rdat%fq3(1,9)                       )*rdat%aqx
      fwk(14,5)= (rdat%fq5(4,5)     -rdat%fq4(2,9)* 3                    )*rdat%aqxy
      fwk(15,5)=-(rdat%fq5(5,5)     -rdat%fq4(3,9)* 6 +rdat%fq3(1,9)* 3      )*rdat%aqx
      fwk(16,5)=-(rdat%fq5(1,5)*acy4-rdat%fq4(1,9)*10*rdat%acy2+rdat%fq3(1,9)*15)*rdat%acy
      fwk(17,5)=  rdat%fq5(2,5)*acy4-rdat%fq4(2,9)* 6 *rdat%acy2+rdat%fq3(2,9)* 3
      fwk(18,5)=-(rdat%fq5(3,5)*rdat%acy2-rdat%fq4(1,9)*rdat%acy2-rdat%fq4(3,9)* 3            &
     &                          +rdat%fq3(1,9)* 3                    )*rdat%acy
      fwk(19,5)=  rdat%fq5(4,5)*rdat%acy2-rdat%fq4(2,9)* 3 *rdat%acy2-rdat%fq4(4,9)           &
     &                          +rdat%fq3(2,9)* 3
      fwk(20,5)=-(rdat%fq5(5,5)     -rdat%fq4(3,9)* 6 +rdat%fq3(1,9)* 3      )*rdat%acy
      fwk(21,5)=  rdat%fq5(6,5)     -rdat%fq4(4,9)*10+rdat%fq3(2,9)*15
      do i= 5, 6
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,5)*work(i+ 9)
         end do
      end do

      call frikr6(rdat, 1, 1,fw6,rdat%fq6, 4,rdat%fq5, 8,rdat%fq4, 8,rdat%fq3)

      do j= 1,28
         rdat%r06(j,1)= rdat%r06(j,1)+fw6(j)*fcc(j,6)
      enddo

      end subroutine intk_18

! >
! >    @brief   dpdp case
! >
! >    @details integration of the dpdp case
! >
      subroutine intk_19(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd1, xmd2, xmd3, xmd4, xmd5
      real(kind=dp) :: xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: fqd11, fqd12, fqd13, fqd21, &
                       fqd22, fqd23, fqd24, fqd25, fqd26, &
                       fqd31, fqd32, fqd33, fqd34
      real(kind=dp) :: y33, y34
      real(kind=dp) :: work, fwk, fw6, fcu, fcc
      integer :: i, j

! Generate jtype=19 integrals

      dimension  work(12),fwk(21,12),fw6(28)
      dimension  fcu(45,8),fcc(45,8)

      if(ikl == 0) then
         rdat%r00= 0.0_dp
         rdat%r01= 0.0_dp
         rdat%r02(:,1:46)= 0.0_dp
         rdat%r03(:,1:34)= 0.0_dp
         rdat%r04(:,1:17)= 0.0_dp
         rdat%r05(:,1:6)= 0.0_dp
         rdat%r06(:,1)= 0.0_dp
         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd1= xmd2
      xmd3= xmd2
      xmd4= xmd2*xmd2
      xmd5= xmd4
      xmdt= xmd4*xmd3

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      y33 = rdat%y03 *rdat%y03
      y34 =-rdat%y03 *rdat%y04
      work( 1)= xmd1
      work( 2)= y33
      work( 3)=-xmd3*rdat%y03
      work( 4)= xmd3*rdat%y04
      work( 5)= y33*rdat%y04

      work( 6)=-xmd1*rdat%y03
      work( 7)= xmd5
      work( 8)= xmd3*y33
      work( 9)= xmd3*y34

      work(10)= xmd4
      work(11)=-xmd5*rdat%y03
      work(12)= xmd5*rdat%y04

      call fcufcc(rdat,6,xmdt,fcu,fcc)

      do i= 1, 5
         do j= 1, 5
            rdat%r00(j,i)= rdat%r00(j,i)+rdat%fq0(i)*work(j)
         end do
      end do

      do i= 1, 8
         fwk(1,i)=-rdat%fq1(1,i)*rdat%aqx
         fwk(2,i)=-rdat%fq1(1,i)*rdat%acy
         fwk(3,i)= rdat%fq1(2,i)
      enddo
      fqd11 = rdat%fq1(1,5)*rdat%rab +rdat%fq1(1, 7)
      fqd12 = rdat%fq1(2,5)*rdat%rab +rdat%fq1(2, 7)
         fwk(1,9)=-fqd11*rdat%aqx
         fwk(2,9)=-fqd11*rdat%acy
         fwk(3,9)= fqd12
      do i= 1, 4
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,1)*work(i+ 5)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,1)*work(i+ 5)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,1)*work(i+ 5)
      enddo
      do i= 5, 8
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,2)*work(i+ 1)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,2)*work(i+ 1)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,2)*work(i+ 1)
      enddo
      do i= 9,12
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,3)*work(i- 3)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,3)*work(i- 3)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,3)*work(i- 3)
      enddo
      do i=13,16
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,4)*work(i- 7)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,4)*work(i- 7)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,4)*work(i- 7)
      enddo
      do i=17,20
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,5)*work(i-11)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,5)*work(i-11)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,5)*work(i-11)
      enddo
      do i=21,25
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,6)*work(i-20)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,6)*work(i-20)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,6)*work(i-20)
      enddo
      do i=26,30
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,7)*work(i-25)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,7)*work(i-25)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,7)*work(i-25)
      enddo
      do i=31,35
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,8)*work(i-30)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,8)*work(i-30)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,8)*work(i-30)
      enddo
      do i=36,40
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,9)*work(i-35)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,9)*work(i-35)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,9)*work(i-35)
      enddo

      do i= 1,10
         fwk(1,i)= rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i)
         fwk(2,i)= rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i)
         fwk(3,i)= rdat%fq2(3,i)     -rdat%fq1(1,i)
         fwk(4,i)= rdat%fq2(1,i)*rdat%aqxy
         fwk(5,i)=-rdat%fq2(2,i)*rdat%aqx
         fwk(6,i)=-rdat%fq2(2,i)*rdat%acy
      enddo
      fqd21 = rdat%fq2(1,5)*rdat%rab +rdat%fq2(1, 7)
      fqd22 = rdat%fq2(2,5)*rdat%rab +rdat%fq2(2, 7)
      fqd23 = rdat%fq2(3,5)*rdat%rab +rdat%fq2(3, 7)
         fwk(1,11)= fqd21*rdat%aqx2-fqd11
         fwk(2,11)= fqd21*rdat%acy2-fqd11
         fwk(3,11)= fqd23     -fqd11
         fwk(4,11)= fqd21*rdat%aqxy
         fwk(5,11)=-fqd22*rdat%aqx
         fwk(6,11)=-fqd22*rdat%acy
      fqd13 = rdat%fq1(1,8)*rdat%rab +rdat%fq1(1,10)
      fqd24 = rdat%fq2(1,8)*rdat%rab +rdat%fq2(1,10)
      fqd25 = rdat%fq2(2,8)*rdat%rab +rdat%fq2(2,10)
      fqd26 = rdat%fq2(3,8)*rdat%rab +rdat%fq2(3,10)
         fwk(1,12)= fqd24*rdat%aqx2-fqd13
         fwk(2,12)= fqd24*rdat%acy2-fqd13
         fwk(3,12)= fqd26     -fqd13
         fwk(4,12)= fqd24*rdat%aqxy
         fwk(5,12)=-fqd25*rdat%aqx
         fwk(6,12)=-fqd25*rdat%acy
      do i= 1, 3
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,1)*work(i+ 9)
         enddo
      enddo
      do i= 4, 6
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,2)*work(i+ 6)
         enddo
      enddo
      do i= 7, 9
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,3)*work(i+ 3)
         enddo
      enddo
      do i=10,12
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,4)*work(i)
         enddo
      enddo
      do i=13,15
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,5)*work(i- 3)
         enddo
      enddo
      do i=16,19
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,6)*work(i-10)
         enddo
      enddo
      do i=20,23
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,7)*work(i-14)
         enddo
      enddo
      do i=24,27
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,8)*work(i-18)
         enddo
      enddo
      do i=28,31
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,11)*work(i-22)
         enddo
      enddo
      do i=32,36
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,9)*work(i-31)
         enddo
      enddo
      do i=37,41
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,10)*work(i-36)
         enddo
      enddo
      do i=42,46
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,12)*work(i-41)
         enddo
      enddo

      do i= 1, 5
         rdat%r03( 1,i)= rdat%r03( 1,i)+(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)* 3 )*xmdtx
         rdat%r03( 2,i)= rdat%r03( 2,i)+(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)    )*xmdty
         rdat%r03( 3,i)= rdat%r03( 3,i)+(rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i)    )*xmdt
         rdat%r03( 4,i)= rdat%r03( 4,i)+(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)    )*xmdtx
         rdat%r03( 5,i)= rdat%r03( 5,i)+ rdat%fq3(2,i)                    *xmdtxy
         rdat%r03( 6,i)= rdat%r03( 6,i)+(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*xmdtx
         rdat%r03( 7,i)= rdat%r03( 7,i)+(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)* 3 )*xmdty
         rdat%r03( 8,i)= rdat%r03( 8,i)+(rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i)    )*xmdt
         rdat%r03( 9,i)= rdat%r03( 9,i)+(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*xmdty
         rdat%r03(10,i)= rdat%r03(10,i)+(rdat%fq3(4,i)     -rdat%fq2(2,i)* 3 )*xmdt
      enddo
      rdat%fq2(1,5)= fqd21
      rdat%fq2(2,5)= fqd22
      rdat%fq3(1,5)= rdat%fq3(1,5)*rdat%rab +rdat%fq3(1, 7)
      rdat%fq3(2,5)= rdat%fq3(2,5)*rdat%rab +rdat%fq3(2, 7)
      rdat%fq3(3,5)= rdat%fq3(3,5)*rdat%rab +rdat%fq3(3, 7)
      rdat%fq3(4,5)= rdat%fq3(4,5)*rdat%rab +rdat%fq3(4, 7)
      do i= 5,11
         fwk( 1,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)* 3 )*rdat%aqx
         fwk( 2,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)    )*rdat%acy
         fwk( 3,i)=  rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i)
         fwk( 4,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)    )*rdat%aqx
         fwk( 5,i)=  rdat%fq3(2,i)                    *rdat%aqxy
         fwk( 6,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%aqx
         fwk( 7,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)* 3 )*rdat%acy
         fwk( 8,i)=  rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i)
         fwk( 9,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%acy
         fwk(10,i)=  rdat%fq3(4,i)     -rdat%fq2(2,i)* 3
      enddo
      fqd31 = rdat%fq3(1,8)*rdat%rab +rdat%fq3(1,10)
      fqd32 = rdat%fq3(2,8)*rdat%rab +rdat%fq3(2,10)
      fqd33 = rdat%fq3(3,8)*rdat%rab +rdat%fq3(3,10)
      fqd34 = rdat%fq3(4,8)*rdat%rab +rdat%fq3(4,10)
         fwk( 1,12)=-(fqd31*rdat%aqx2-fqd24* 3 )*rdat%aqx
         fwk( 2,12)=-(fqd31*rdat%aqx2-fqd24    )*rdat%acy
         fwk( 3,12)=  fqd32*rdat%aqx2-fqd25
         fwk( 4,12)=-(fqd31*rdat%acy2-fqd24    )*rdat%aqx
         fwk( 5,12)=  fqd32                *rdat%aqxy
         fwk( 6,12)=-(fqd33     -fqd24    )*rdat%aqx
         fwk( 7,12)=-(fqd31*rdat%acy2-fqd24* 3 )*rdat%acy
         fwk( 8,12)=  fqd32*rdat%acy2-fqd25
         fwk( 9,12)=-(fqd33     -fqd24    )*rdat%acy
         fwk(10,12)=  fqd34     -fqd25* 3
      do i= 6, 8
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,5)*work(i+ 4)
         enddo
      enddo
      do i= 9,11
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,6)*work(i+ 1)
         enddo
      enddo
      do i=12,14
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,7)*work(i- 2)
         enddo
      enddo
      do i=15,17
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,8)*work(i- 5)
         enddo
      enddo
      do i=18,21
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,9)*work(i-12)
         enddo
      enddo
      do i=22,25
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,10)*work(i-16)
         enddo
      enddo
      do i=26,29
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,12)*work(i-20)
         enddo
      enddo
      do i=30,34
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,11)*work(i-29)
         enddo
      enddo

      do j= 1, 5
         rdat%fq4(j,1)= rdat%fq4(j,1)*rdat%rab +rdat%fq4(j,3)
      enddo
      aqx4= rdat%aqx2*rdat%aqx2
      acy4= rdat%acy2*rdat%acy2
      x2y2= rdat%aqx2*rdat%acy2
      q2c2= rdat%aqx2+rdat%acy2
      do i= 1, 4
         rdat%r04( 1,i)= rdat%r04( 1,i)+(rdat%fq4(1,i)*aqx4-rdat%fq3(1,i+4)* 6 *rdat%aqx2      &
     &                                       +rdat%fq2(1,i+4)* 3   )*xmdt
         rdat%r04( 2,i)= rdat%r04( 2,i)+(rdat%fq4(1,i)*rdat%aqx2-rdat%fq3(1,i+4)* 3   )*xmdtxy
         rdat%r04( 3,i)= rdat%r04( 3,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i+4)* 3   )*xmdtx
         rdat%r04( 4,i)= rdat%r04( 4,i)+(rdat%fq4(1,i)*x2y2-rdat%fq3(1,i+4)*q2c2          &
     &                                       +rdat%fq2(1,i+4)      )*xmdt
         rdat%r04( 5,i)= rdat%r04( 5,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i+4)      )*xmdty
         rdat%r04( 6,i)= rdat%r04( 6,i)+(rdat%fq4(3,i)*rdat%aqx2-rdat%fq3(1,i+4)*rdat%aqx2          &
     &                           -rdat%fq3(3,i+4)+rdat%fq2(1,i+4)      )*xmdt
         rdat%r04( 7,i)= rdat%r04( 7,i)+(rdat%fq4(1,i)*rdat%acy2-rdat%fq3(1,i+4)* 3   )*xmdtxy
         rdat%r04( 8,i)= rdat%r04( 8,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i+4)      )*xmdtx
         rdat%r04( 9,i)= rdat%r04( 9,i)+(rdat%fq4(3,i)     -rdat%fq3(1,i+4)      )*xmdtxy
         rdat%r04(10,i)= rdat%r04(10,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i+4)* 3   )*xmdtx
         rdat%r04(11,i)= rdat%r04(11,i)+(rdat%fq4(1,i)*acy4-rdat%fq3(1,i+4)* 6 *rdat%acy2      &
     &                                       +rdat%fq2(1,i+4)* 3   )*xmdt
         rdat%r04(12,i)= rdat%r04(12,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i+4)* 3   )*xmdty
         rdat%r04(13,i)= rdat%r04(13,i)+(rdat%fq4(3,i)*rdat%acy2-rdat%fq3(1,i+4)*rdat%acy2          &
     &                           -rdat%fq3(3,i+4)+rdat%fq2(1,i+4)      )*xmdt
         rdat%r04(14,i)= rdat%r04(14,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i+4)* 3   )*xmdty
         rdat%r04(15,i)= rdat%r04(15,i)+(rdat%fq4(5,i)-rdat%fq3(3,i+4)* 6                 &
     &                                       +rdat%fq2(1,i+4)* 3   )*xmdt
      enddo
      rdat%fq2(1,8)= fqd24
      rdat%fq3(1,8)= fqd31
      rdat%fq3(2,8)= fqd32
      rdat%fq3(3,8)= fqd33
      do j= 1, 5
         rdat%fq4(j,4)= rdat%fq4(j,4)*rdat%rab +rdat%fq4(j,6)
      enddo
      do i= 4, 7
         fwk( 1,i)=  rdat%fq4(1,i)*aqx4-rdat%fq3(1,i+4)* 6 *rdat%aqx2                &
     &                                              +rdat%fq2(1,i+4)* 3
         fwk( 2,i)= (rdat%fq4(1,i)*rdat%aqx2-rdat%fq3(1,i+4)* 3              )*rdat%aqxy
         fwk( 3,i)=-(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i+4)* 3              )*rdat%aqx
         fwk( 4,i)=  rdat%fq4(1,i)*x2y2-rdat%fq3(1,i+4)*q2c2+rdat%fq2(1,i+4)
         fwk( 5,i)=-(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i+4)                 )*rdat%acy
         fwk( 6,i)=  rdat%fq4(3,i)*rdat%aqx2-rdat%fq3(1,i+4)*rdat%aqx2-rdat%fq3(3,i+4)        &
     &                                              +rdat%fq2(1,i+4)
         fwk( 7,i)= (rdat%fq4(1,i)*rdat%acy2-rdat%fq3(1,i+4)* 3              )*rdat%aqxy
         fwk( 8,i)=-(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i+4)                 )*rdat%aqx
         fwk( 9,i)= (rdat%fq4(3,i)     -rdat%fq3(1,i+4)                 )*rdat%aqxy
         fwk(10,i)=-(rdat%fq4(4,i)     -rdat%fq3(2,i+4)* 3              )*rdat%aqx
         fwk(11,i)=  rdat%fq4(1,i)*acy4-rdat%fq3(1,i+4)* 6 *rdat%acy2                &
     &                                              +rdat%fq2(1,i+4)* 3
         fwk(12,i)=-(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i+4)* 3              )*rdat%acy
         fwk(13,i)=  rdat%fq4(3,i)*rdat%acy2-rdat%fq3(1,i+4)*rdat%acy2-rdat%fq3(3,i+4)        &
     &                                              +rdat%fq2(1,i+4)
         fwk(14,i)=-(rdat%fq4(4,i)     -rdat%fq3(2,i+4)* 3              )*rdat%acy
         fwk(15,i)=  rdat%fq4(5,i)     -rdat%fq3(3,i+4)* 6  +rdat%fq2(1,i+4)* 3
      enddo
      do i= 5, 7
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,5)*work(i+ 5)
         enddo
      enddo
      do i= 8,10
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,6)*work(i+ 2)
         enddo
      enddo
      do i=11,13
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,4)*work(i- 1)
         enddo
      enddo
      do i=14,17
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,7)*work(i- 8)
         enddo
      enddo

      do j= 1, 6
         rdat%fq5(j,1)= rdat%fq5(j,1)*rdat%rab +rdat%fq5(j,3)
      enddo
      do i= 1, 3
         rdat%r05( 1,i)= rdat%r05( 1,i)+(rdat%fq5(1,i)*aqx4-rdat%fq4(1,i+3)*10*rdat%aqx2      &
     &                        +rdat%fq3(1,i+7)*15                  )*xmdtx
         rdat%r05( 2,i)= rdat%r05( 2,i)+(rdat%fq5(1,i)*aqx4-rdat%fq4(1,i+3)* 6 *rdat%aqx2      &
     &                        +rdat%fq3(1,i+7)* 3                   )*xmdty
         rdat%r05( 3,i)= rdat%r05( 3,i)+(rdat%fq5(2,i)*aqx4-rdat%fq4(2,i+3)* 6 *rdat%aqx2      &
     &                        +rdat%fq3(2,i+7)* 3                   )*xmdt
         rdat%r05( 4,i)= rdat%r05( 4,i)+(rdat%fq5(1,i)*x2y2-rdat%fq4(1,i+3)* 3 *rdat%acy2      &
     &                        -rdat%fq4(1,i+3)*rdat%aqx2+rdat%fq3(1,i+7)* 3  )*xmdtx
         rdat%r05( 5,i)= rdat%r05( 5,i)+(rdat%fq5(2,i)*rdat%aqx2-rdat%fq4(2,i+3)* 3    )*xmdtxy
         rdat%r05( 6,i)= rdat%r05( 6,i)+(rdat%fq5(3,i)*rdat%aqx2-rdat%fq4(1,i+3)*rdat%aqx2          &
     &                        -rdat%fq4(3,i+3)* 3 +rdat%fq3(1,i+7)* 3   )*xmdtx
         rdat%r05( 7,i)= rdat%r05( 7,i)+(rdat%fq5(1,i)*x2y2-rdat%fq4(1,i+3)* 3 *rdat%aqx2      &
     &                        -rdat%fq4(1,i+3)*rdat%acy2+rdat%fq3(1,i+7)* 3  )*xmdty
         rdat%r05( 8,i)= rdat%r05( 8,i)+(rdat%fq5(2,i)*x2y2-rdat%fq4(2,i+3)*q2c2          &
     &                        +rdat%fq3(2,i+7)                      )*xmdt
         rdat%r05( 9,i)= rdat%r05( 9,i)+(rdat%fq5(3,i)*rdat%aqx2-rdat%fq4(1,i+3)*rdat%aqx2          &
     &                        -rdat%fq4(3,i+3)+rdat%fq3(1,i+7)          )*xmdty
         rdat%r05(10,i)= rdat%r05(10,i)+(rdat%fq5(4,i)*rdat%aqx2-rdat%fq4(2,i+3)* 3 *rdat%aqx2      &
     &                        -rdat%fq4(4,i+3)+rdat%fq3(2,i+7)* 3       )*xmdt
         rdat%r05(11,i)= rdat%r05(11,i)+(rdat%fq5(1,i)*acy4-rdat%fq4(1,i+3)* 6 *rdat%acy2      &
     &                        +rdat%fq3(1,i+7)* 3                   )*xmdtx
         rdat%r05(12,i)= rdat%r05(12,i)+(rdat%fq5(2,i)*rdat%acy2-rdat%fq4(2,i+3)* 3    )*xmdtxy
         rdat%r05(13,i)= rdat%r05(13,i)+(rdat%fq5(3,i)*rdat%acy2-rdat%fq4(1,i+3)*rdat%acy2          &
     &                        -rdat%fq4(3,i+3)+rdat%fq3(1,i+7)          )*xmdtx
         rdat%r05(14,i)= rdat%r05(14,i)+(rdat%fq5(4,i)-rdat%fq4(2,i+3)* 3         )*xmdtxy
         rdat%r05(15,i)= rdat%r05(15,i)+(rdat%fq5(5,i)-rdat%fq4(3,i+3)* 6                 &
     &                        +rdat%fq3(1,i+7)* 3                   )*xmdtx
         rdat%r05(16,i)= rdat%r05(16,i)+(rdat%fq5(1,i)*acy4-rdat%fq4(1,i+3)*10*rdat%acy2      &
     &                        +rdat%fq3(1,i+7)*15                  )*xmdty
         rdat%r05(17,i)= rdat%r05(17,i)+(rdat%fq5(2,i)*acy4-rdat%fq4(2,i+3)* 6 *rdat%acy2      &
     &                        +rdat%fq3(2,i+7)* 3                   )*xmdt
         rdat%r05(18,i)= rdat%r05(18,i)+(rdat%fq5(3,i)*rdat%acy2-rdat%fq4(1,i+3)*rdat%acy2          &
     &                        -rdat%fq4(3,i+3)* 3 +rdat%fq3(1,i+7)* 3   )*xmdty
         rdat%r05(19,i)= rdat%r05(19,i)+(rdat%fq5(4,i)*rdat%acy2-rdat%fq4(2,i+3)* 3 *rdat%acy2      &
     &                        -rdat%fq4(4,i+3)+rdat%fq3(2,i+7)* 3       )*xmdt
         rdat%r05(20,i)= rdat%r05(20,i)+(rdat%fq5(5,i)-rdat%fq4(3,i+3)* 6                 &
     &                        +rdat%fq3(1,i+7)* 3                   )*xmdty
         rdat%r05(21,i)= rdat%r05(21,i)+(rdat%fq5(6,i)-rdat%fq4(4,i+3)*10                &
     &                        +rdat%fq3(2,i+7)*15                  )*xmdt
      enddo
      fwk( 1,4)=-(rdat%fq5(1,4)*aqx4-rdat%fq4(1,7)*10*rdat%aqx2+rdat%fq3(1,11)*15)*rdat%aqx
      fwk( 2,4)=-(rdat%fq5(1,4)*aqx4-rdat%fq4(1,7)* 6 *rdat%aqx2+rdat%fq3(1,11)* 3 )*rdat%acy
      fwk( 3,4)=  rdat%fq5(2,4)*aqx4-rdat%fq4(2,7)* 6 *rdat%aqx2+rdat%fq3(2,11)* 3
      fwk( 4,4)=-(rdat%fq5(1,4)*x2y2-rdat%fq4(1,7)*rdat%aqx2-rdat%fq4(1,7)* 3 *rdat%acy2      &
     &                                         +rdat%fq3(1,11)* 3     )*rdat%aqx
      fwk( 5,4)= (rdat%fq5(2,4)*rdat%aqx2-rdat%fq4(2,7)* 3                     )*rdat%aqxy
      fwk( 6,4)=-(rdat%fq5(3,4)*rdat%aqx2-rdat%fq4(1,7)*rdat%aqx2-rdat%fq4(3,7)* 3            &
     &            +rdat%fq3(1,11)* 3                                  )*rdat%aqx
      fwk( 7,4)=-(rdat%fq5(1,4)*x2y2-rdat%fq4(1,7)* 3 *rdat%aqx2-rdat%fq4(1,7)*rdat%acy2      &
     &            +rdat%fq3(1,11)* 3                                  )*rdat%acy
      fwk( 8,4)=  rdat%fq5(2,4)*x2y2-rdat%fq4(2,7)*q2c2+rdat%fq3(2,11)
      fwk( 9,4)=-(rdat%fq5(3,4)*rdat%aqx2-rdat%fq4(1,7)*rdat%aqx2-rdat%fq4(3,7)               &
     &                                             +rdat%fq3(1,11)    )*rdat%acy
      fwk(10,4)=  rdat%fq5(4,4)*rdat%aqx2-rdat%fq4(2,7)* 3 *rdat%aqx2-rdat%fq4(4,7)           &
     &                                             +rdat%fq3(2,11)* 3
      fwk(11,4)=-(rdat%fq5(1,4)*acy4-rdat%fq4(1,7)* 6 *rdat%acy2                     &
     &                                             +rdat%fq3(1,11)* 3 )*rdat%aqx
      fwk(12,4)= (rdat%fq5(2,4)*rdat%acy2-rdat%fq4(2,7)* 3                     )*rdat%aqxy
      fwk(13,4)=-(rdat%fq5(3,4)*rdat%acy2-rdat%fq4(1,7)*rdat%acy2-rdat%fq4(3,7)               &
     &                                             +rdat%fq3(1,11)    )*rdat%aqx
      fwk(14,4)= (rdat%fq5(4,4)     -rdat%fq4(2,7)* 3                     )*rdat%aqxy
      fwk(15,4)=-(rdat%fq5(5,4)     -rdat%fq4(3,7)* 6 +rdat%fq3(1,11)* 3      )*rdat%aqx
      fwk(16,4)=-(rdat%fq5(1,4)*acy4-rdat%fq4(1,7)*10*rdat%acy2                     &
     &                                             +rdat%fq3(1,11)*15)*rdat%acy
      fwk(17,4)=  rdat%fq5(2,4)*acy4-rdat%fq4(2,7)* 6 *rdat%acy2+rdat%fq3(2,11)* 3
      fwk(18,4)=-(rdat%fq5(3,4)*rdat%acy2-rdat%fq4(1,7)*rdat%acy2-rdat%fq4(3,7)* 3            &
     &                                             +rdat%fq3(1,11)* 3 )*rdat%acy
      fwk(19,4)=  rdat%fq5(4,4)*rdat%acy2-rdat%fq4(2,7)* 3 *rdat%acy2-rdat%fq4(4,7)           &
     &                                             +rdat%fq3(2,11)* 3
      fwk(20,4)=-(rdat%fq5(5,4)     -rdat%fq4(3,7)* 6 +rdat%fq3(1,11)* 3      )*rdat%acy
      fwk(21,4)=  rdat%fq5(6,4)     -rdat%fq4(4,7)*10+rdat%fq3(2,11)*15
      do i= 4, 6
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,4)*work(i+ 6)
         end do
      end do

      call frikr6(rdat, 1, 1,fw6,rdat%fq6, 3,rdat%fq5, 6,rdat%fq4,10,rdat%fq3)

      do j= 1,28
         rdat%r06(j,1)= rdat%r06(j,1)+fw6(j)*fcc(j,6)
      enddo

      end subroutine intk_19

! >
! >    @brief   dddp case
! >
! >    @details integration of the dddp case
! >
      subroutine intk_20(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd2, xmd3, xmd4, xmd6
      real(kind=dp) :: xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: y33, y34, y44
      real(kind=dp) :: work, fwk, fw6, fw7, fcu, fcc
      integer :: i, j

! Generate jtype=20 integrals

      dimension  work(15),fwk(28,13),fw6(28, 4),fw7(36)
      dimension  fcu(45,8),fcc(45,8)

      if(ikl == 0) then
         rdat%r00= 0.0_dp
         rdat%r01= 0.0_dp
         rdat%r02(:,1:51)= 0.0_dp
         rdat%r03(:,1:43)= 0.0_dp
         rdat%r04(:,1:29)= 0.0_dp
         rdat%r05(:,1:14)= 0.0_dp
         rdat%r06(:,1:5)= 0.0_dp
         rdat%r07(:,1)= 0.0_dp
         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd3= xmd2
      xmd4= xmd3*xmd2
      xmd6= xmd4*xmd2
      xmdt= xmd6*xmd2
      xmd2= xmd3

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      y33 = rdat%y03 *rdat%y03
      y34 =-rdat%y03 *rdat%y04
      y44 = rdat%y04 *rdat%y04
      work( 1)= xmd4
      work( 2)= xmd2*y33
      work( 3)= xmd2*y34
      work( 4)= xmd2*y44
      work( 5)= y33 *y44

      work( 6)=-xmd4*rdat%y03
      work( 7)= xmd4*rdat%y04
      work( 8)= xmd2*y33*rdat%y04
      work( 9)= xmd2*y34*rdat%y04

      work(10)= xmd6
      work(11)= xmd4*y33
      work(12)= xmd4*y34
      work(13)= xmd4*y44

      work(14)=-xmd6*rdat%y03
      work(15)= xmd6*rdat%y04

      call fcufcc(rdat,7,xmdt,fcu,fcc)

      do i= 1, 5
         do j= 1, 5
            rdat%r00(j,i)= rdat%r00(j,i)+rdat%fq0(i)*work(j)
         end do
      end do

      do i= 1, 8
         fwk(1,i)=-rdat%fq1(1,i)*rdat%aqx
         fwk(2,i)=-rdat%fq1(1,i)*rdat%acy
         fwk(3,i)= rdat%fq1(2,i)
      enddo
      rdat%fq1(1,12)= rdat%fq1(1, 8)*rdat%rab +rdat%fq1(1,10)
      rdat%fq1(1,13)= rdat%fq1(1, 5)*rdat%rab +rdat%fq1(1, 7)
      rdat%fq1(2,13)= rdat%fq1(2, 5)*rdat%rab +rdat%fq1(2, 7)
         fwk(1,9)=-rdat%fq1(1,13)*rdat%aqx
         fwk(2,9)=-rdat%fq1(1,13)*rdat%acy
         fwk(3,9)= rdat%fq1(2,13)
      do i= 1, 4
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,1)*work(i+ 5)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,1)*work(i+ 5)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,1)*work(i+ 5)
      enddo
      do i= 5, 8
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,2)*work(i+ 1)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,2)*work(i+ 1)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,2)*work(i+ 1)
      enddo
      do i= 9,12
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,3)*work(i- 3)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,3)*work(i- 3)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,3)*work(i- 3)
      enddo
      do i=13,16
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,4)*work(i- 7)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,4)*work(i- 7)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,4)*work(i- 7)
      enddo
      do i=17,20
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,5)*work(i-11)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,5)*work(i-11)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,5)*work(i-11)
      enddo
      do i=21,25
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,6)*work(i-20)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,6)*work(i-20)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,6)*work(i-20)
      enddo
      do i=26,30
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,7)*work(i-25)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,7)*work(i-25)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,7)*work(i-25)
      enddo
      do i=31,35
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,8)*work(i-30)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,8)*work(i-30)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,8)*work(i-30)
      enddo
      do i=36,40
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,9)*work(i-35)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,9)*work(i-35)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,9)*work(i-35)
      enddo

      do i= 1,10
         fwk(1,i)= rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i)
         fwk(2,i)= rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i)
         fwk(3,i)= rdat%fq2(3,i)     -rdat%fq1(1,i)
         fwk(4,i)= rdat%fq2(1,i)*rdat%aqxy
         fwk(5,i)=-rdat%fq2(2,i)*rdat%aqx
         fwk(6,i)=-rdat%fq2(2,i)*rdat%acy
      enddo
      do i= 1, 3
         rdat%fq2(i,12)= rdat%fq2(i,8)*rdat%rab +rdat%fq2(i,10)
         rdat%fq2(i,13)= rdat%fq2(i,5)*rdat%rab +rdat%fq2(i, 7)
      enddo
      do i=12,13
         fwk(1,i)= rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i)
         fwk(2,i)= rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i)
         fwk(3,i)= rdat%fq2(3,i)     -rdat%fq1(1,i)
         fwk(4,i)= rdat%fq2(1,i)*rdat%aqxy
         fwk(5,i)=-rdat%fq2(2,i)*rdat%aqx
         fwk(6,i)=-rdat%fq2(2,i)*rdat%acy
      enddo
      do i= 1, 4
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,1)*work(i+ 9)
         enddo
      enddo
      do i= 5, 8
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,2)*work(i+ 5)
         enddo
      enddo
      do i= 9,12
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,3)*work(i+ 1)
         enddo
      enddo
      do i=13,16
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,4)*work(i- 3)
         enddo
      enddo
      do i=17,20
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,5)*work(i- 7)
         enddo
      enddo
      do i=21,24
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,6)*work(i-15)
         enddo
      enddo
      do i=25,28
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,7)*work(i-19)
         enddo
      enddo
      do i=29,32
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,8)*work(i-23)
         enddo
      enddo
      do i=33,36
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,13)*work(i-27)
         enddo
      enddo
      do i=37,41
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,9)*work(i-36)
         enddo
      enddo
      do i=42,46
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,10)*work(i-41)
         enddo
      enddo
      do i=47,51
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,12)*work(i-46)
         enddo
      enddo

      do j= 1, 4
         rdat%fq3(j,12)= rdat%fq3(j,8)*rdat%rab +rdat%fq3(j,10)
         rdat%fq3(j,13)= rdat%fq3(j,5)*rdat%rab +rdat%fq3(j, 7)
      enddo
      do i= 1,13
         fwk( 1,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)* 3 )*rdat%aqx
         fwk( 2,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)    )*rdat%acy
         fwk( 3,i)=  rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i)
         fwk( 4,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)    )*rdat%aqx
         fwk( 5,i)=  rdat%fq3(2,i)                    *rdat%aqxy
         fwk( 6,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%aqx
         fwk( 7,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)* 3 )*rdat%acy
         fwk( 8,i)=  rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i)
         fwk( 9,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%acy
         fwk(10,i)=  rdat%fq3(4,i)     -rdat%fq2(2,i)* 3
      enddo
      do i= 1, 2
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,1)*work(i+13)
         enddo
      enddo
      do i= 3, 4
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,2)*work(i+11)
         enddo
      enddo
      do i= 5, 6
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,3)*work(i+ 9)
         enddo
      enddo
      do i= 7, 8
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,4)*work(i+ 7)
         enddo
      enddo
      do i= 9,10
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,5)*work(i+ 5)
         enddo
      enddo
      do i=11,14
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,6)*work(i- 1)
         enddo
      enddo
      do i=15,18
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,7)*work(i- 5)
         enddo
      enddo
      do i=19,22
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,8)*work(i- 9)
         enddo
      enddo
      do i=23,26
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,13)*work(i-13)
         enddo
      enddo
      do i=27,30
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,9)*work(i-21)
         enddo
      enddo
      do i=31,34
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,10)*work(i-25)
         enddo
      enddo
      do i=35,38
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,12)*work(i-29)
         enddo
      enddo
      do i=39,43
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,11)*work(i-38)
         enddo
      enddo

      aqx4= rdat%aqx2*rdat%aqx2
      acy4= rdat%acy2*rdat%acy2
      x2y2= rdat%aqx2*rdat%acy2
      q2c2= rdat%aqx2+rdat%acy2
      do i= 1, 5
         rdat%r04( 1,i)= rdat%r04( 1,i)+(rdat%fq4(1,i)*aqx4-rdat%fq3(1,i)* 6 *rdat%aqx2        &
     &                                       +rdat%fq2(1,i)* 3    )*xmdt
         rdat%r04( 2,i)= rdat%r04( 2,i)+(rdat%fq4(1,i)*rdat%aqx2-rdat%fq3(1,i)* 3    )*xmdtxy
         rdat%r04( 3,i)= rdat%r04( 3,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)* 3    )*xmdtx
         rdat%r04( 4,i)= rdat%r04( 4,i)+(rdat%fq4(1,i)*x2y2-rdat%fq3(1,i)*q2c2            &
     &                                       +rdat%fq2(1,i)       )*xmdt
         rdat%r04( 5,i)= rdat%r04( 5,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)       )*xmdty
         rdat%r04( 6,i)= rdat%r04( 6,i)+(rdat%fq4(3,i)*rdat%aqx2-rdat%fq3(1,i)*rdat%aqx2            &
     &                                    -rdat%fq3(3,i)+rdat%fq2(1,i))*xmdt
         rdat%r04( 7,i)= rdat%r04( 7,i)+(rdat%fq4(1,i)*rdat%acy2-rdat%fq3(1,i)* 3    )*xmdtxy
         rdat%r04( 8,i)= rdat%r04( 8,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)       )*xmdtx
         rdat%r04( 9,i)= rdat%r04( 9,i)+(rdat%fq4(3,i)     -rdat%fq3(1,i)       )*xmdtxy
         rdat%r04(10,i)= rdat%r04(10,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3    )*xmdtx
         rdat%r04(11,i)= rdat%r04(11,i)+(rdat%fq4(1,i)*acy4-rdat%fq3(1,i)* 6 *rdat%acy2        &
     &                                       +rdat%fq2(1,i)* 3    )*xmdt
         rdat%r04(12,i)= rdat%r04(12,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)* 3    )*xmdty
         rdat%r04(13,i)= rdat%r04(13,i)+(rdat%fq4(3,i)*rdat%acy2-rdat%fq3(1,i)*rdat%acy2            &
     &                                    -rdat%fq3(3,i)+rdat%fq2(1,i))*xmdt
         rdat%r04(14,i)= rdat%r04(14,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3    )*xmdty
         rdat%r04(15,i)= rdat%r04(15,i)+(rdat%fq4(5,i)     -rdat%fq3(3,i)* 6              &
     &                                       +rdat%fq2(1,i)* 3    )*xmdt
      enddo
      rdat%fq2(1,5)= rdat%fq2(1,13)
      rdat%fq3(1,5)= rdat%fq3(1,13)
      rdat%fq3(2,5)= rdat%fq3(2,13)
      rdat%fq3(3,5)= rdat%fq3(3,13)
      do j= 1, 5
         rdat%fq4(j, 5)= rdat%fq4(j,5)*rdat%rab +rdat%fq4(j, 7)
         rdat%fq4(j,12)= rdat%fq4(j,8)*rdat%rab +rdat%fq4(j,10)
      enddo
      do i= 5,12
         fwk( 1,i)=  rdat%fq4(1,i)*aqx4-rdat%fq3(1,i)* 6 *rdat%aqx2+rdat%fq2(1,i)* 3
         fwk( 2,i)= (rdat%fq4(1,i)*rdat%aqx2-rdat%fq3(1,i)* 3                 )*rdat%aqxy
         fwk( 3,i)=-(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)* 3                 )*rdat%aqx
         fwk( 4,i)=  rdat%fq4(1,i)*x2y2-rdat%fq3(1,i)*q2c2+rdat%fq2(1,i)
         fwk( 5,i)=-(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)                    )*rdat%acy
         fwk( 6,i)=  rdat%fq4(3,i)*rdat%aqx2-rdat%fq3(1,i)*rdat%aqx2-rdat%fq3(3,i)+rdat%fq2(1,i)
         fwk( 7,i)= (rdat%fq4(1,i)*rdat%acy2-rdat%fq3(1,i)* 3                 )*rdat%aqxy
         fwk( 8,i)=-(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)                    )*rdat%aqx
         fwk( 9,i)= (rdat%fq4(3,i)     -rdat%fq3(1,i)                    )*rdat%aqxy
         fwk(10,i)=-(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3                 )*rdat%aqx
         fwk(11,i)=  rdat%fq4(1,i)*acy4-rdat%fq3(1,i)* 6 *rdat%acy2+rdat%fq2(1,i)* 3
         fwk(12,i)=-(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)* 3                 )*rdat%acy
         fwk(13,i)=  rdat%fq4(3,i)*rdat%acy2-rdat%fq3(1,i)*rdat%acy2-rdat%fq3(3,i)+rdat%fq2(1,i)
         fwk(14,i)=-(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3                 )*rdat%acy
         fwk(15,i)=  rdat%fq4(5,i)     -rdat%fq3(3,i)* 6  +rdat%fq2(1,i)* 3
      enddo
      do i= 6, 7
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,5)*work(i+ 8)
         enddo
      enddo
      do i= 8, 9
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,6)*work(i+ 6)
         enddo
      enddo
      do i=10,11
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,7)*work(i+ 4)
         enddo
      enddo
      do i=12,13
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,8)*work(i+ 2)
         enddo
      enddo
      do i=14,17
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,9)*work(i- 4)
         enddo
      enddo
      do i=18,21
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,10)*work(i- 8)
         enddo
      enddo
      do i=22,25
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,12)*work(i-12)
         enddo
      enddo
      do i=26,29
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,11)*work(i-20)
         enddo
      enddo

      do j= 1, 6
         rdat%fq5(j,1)= rdat%fq5(j,1)*rdat%rab +rdat%fq5(j,3)
      enddo
      do i= 1, 4
         rdat%r05( 1,i)= rdat%r05( 1,i)+(rdat%fq5(1,i)*aqx4-rdat%fq4(1,i+4)*10*rdat%aqx2      &
     &                                       +rdat%fq3(1,i+4)*15  )*xmdtx
         rdat%r05( 2,i)= rdat%r05( 2,i)+(rdat%fq5(1,i)*aqx4-rdat%fq4(1,i+4)* 6 *rdat%aqx2      &
     &                                       +rdat%fq3(1,i+4)* 3   )*xmdty
         rdat%r05( 3,i)= rdat%r05( 3,i)+(rdat%fq5(2,i)*aqx4-rdat%fq4(2,i+4)* 6 *rdat%aqx2      &
     &                                       +rdat%fq3(2,i+4)* 3   )*xmdt
         rdat%r05( 4,i)= rdat%r05( 4,i)+(rdat%fq5(1,i)*x2y2-rdat%fq4(1,i+4)* 3 *rdat%acy2      &
     &                        -rdat%fq4(1,i+4)*rdat%aqx2+rdat%fq3(1,i+4)* 3 )*xmdtx
         rdat%r05( 5,i)= rdat%r05( 5,i)+(rdat%fq5(2,i)*rdat%aqx2-rdat%fq4(2,i+4)* 3   )*xmdtxy
         rdat%r05( 6,i)= rdat%r05( 6,i)+(rdat%fq5(3,i)*rdat%aqx2-rdat%fq4(1,i+4)*rdat%aqx2          &
     &                        -rdat%fq4(3,i+4)* 3 +rdat%fq3(1,i+4)* 3  )*xmdtx
         rdat%r05( 7,i)= rdat%r05( 7,i)+(rdat%fq5(1,i)*x2y2-rdat%fq4(1,i+4)* 3 *rdat%aqx2      &
     &                        -rdat%fq4(1,i+4)*rdat%acy2+rdat%fq3(1,i+4)* 3 )*xmdty
         rdat%r05( 8,i)= rdat%r05( 8,i)+(rdat%fq5(2,i)*x2y2-rdat%fq4(2,i+4)*q2c2          &
     &                                       +rdat%fq3(2,i+4)      )*xmdt
         rdat%r05( 9,i)= rdat%r05( 9,i)+(rdat%fq5(3,i)*rdat%aqx2-rdat%fq4(1,i+4)*rdat%aqx2          &
     &                        -rdat%fq4(3,i+4)+rdat%fq3(1,i+4)         )*xmdty
         rdat%r05(10,i)= rdat%r05(10,i)+(rdat%fq5(4,i)*rdat%aqx2-rdat%fq4(2,i+4)* 3 *rdat%aqx2      &
     &                        -rdat%fq4(4,i+4)+rdat%fq3(2,i+4)* 3      )*xmdt
         rdat%r05(11,i)= rdat%r05(11,i)+(rdat%fq5(1,i)*acy4-rdat%fq4(1,i+4)* 6 *rdat%acy2      &
     &                                       +rdat%fq3(1,i+4)* 3   )*xmdtx
         rdat%r05(12,i)= rdat%r05(12,i)+(rdat%fq5(2,i)*rdat%acy2-rdat%fq4(2,i+4)* 3   )*xmdtxy
         rdat%r05(13,i)= rdat%r05(13,i)+(rdat%fq5(3,i)*rdat%acy2-rdat%fq4(1,i+4)*rdat%acy2          &
     &                        -rdat%fq4(3,i+4)+rdat%fq3(1,i+4)         )*xmdtx
         rdat%r05(14,i)= rdat%r05(14,i)+(rdat%fq5(4,i)-rdat%fq4(2,i+4)* 3        )*xmdtxy
         rdat%r05(15,i)= rdat%r05(15,i)+(rdat%fq5(5,i)-rdat%fq4(3,i+4)* 6                 &
     &                                       +rdat%fq3(1,i+4)* 3   )*xmdtx
         rdat%r05(16,i)= rdat%r05(16,i)+(rdat%fq5(1,i)*acy4-rdat%fq4(1,i+4)*10*rdat%acy2      &
     &                                       +rdat%fq3(1,i+4)*15  )*xmdty
         rdat%r05(17,i)= rdat%r05(17,i)+(rdat%fq5(2,i)*acy4-rdat%fq4(2,i+4)* 6 *rdat%acy2      &
     &                                       +rdat%fq3(2,i+4)* 3   )*xmdt
         rdat%r05(18,i)= rdat%r05(18,i)+(rdat%fq5(3,i)*rdat%acy2-rdat%fq4(1,i+4)*rdat%acy2          &
     &                        -rdat%fq4(3,i+4)* 3 +rdat%fq3(1,i+4)* 3  )*xmdty
         rdat%r05(19,i)= rdat%r05(19,i)+(rdat%fq5(4,i)*rdat%acy2-rdat%fq4(2,i+4)* 3 *rdat%acy2      &
     &                        -rdat%fq4(4,i+4)+rdat%fq3(2,i+4)* 3      )*xmdt
         rdat%r05(20,i)= rdat%r05(20,i)+(rdat%fq5(5,i)-rdat%fq4(3,i+4)* 6                 &
     &                                  +rdat%fq3(1,i+4)* 3        )*xmdty
         rdat%r05(21,i)= rdat%r05(21,i)+(rdat%fq5(6,i)-rdat%fq4(4,i+4)*10                &
     &                                  +rdat%fq3(2,i+4)*15       )*xmdt
      enddo
      rdat%fq3(1,8)= rdat%fq3(1,12)
      rdat%fq3(2,8)= rdat%fq3(2,12)
      rdat%fq4(1,8)= rdat%fq4(1,12)
      rdat%fq4(2,8)= rdat%fq4(2,12)
      rdat%fq4(3,8)= rdat%fq4(3,12)
      rdat%fq4(4,8)= rdat%fq4(4,12)
      do j= 1, 6
         rdat%fq5(j,4)= rdat%fq5(j,4)*rdat%rab +rdat%fq5(j,6)
      enddo
      do i= 4, 7
         fwk( 1,i)=-(rdat%fq5(1,i)*aqx4-rdat%fq4(1,i+4)*10*rdat%aqx2                &
     &                             +rdat%fq3(1,i+4)*15              )*rdat%aqx
         fwk( 2,i)=-(rdat%fq5(1,i)*aqx4-rdat%fq4(1,i+4)* 6 *rdat%aqx2                &
     &                             +rdat%fq3(1,i+4)* 3               )*rdat%acy
         fwk( 3,i)=  rdat%fq5(2,i)*aqx4-rdat%fq4(2,i+4)* 6 *rdat%aqx2                &
     &                             +rdat%fq3(2,i+4)* 3
         fwk( 4,i)=-(rdat%fq5(1,i)*x2y2-rdat%fq4(1,i+4)*rdat%aqx2                    &
     &              -rdat%fq4(1,i+4)* 3 *rdat%acy2+rdat%fq3(1,i+4)* 3         )*rdat%aqx
         fwk( 5,i)= (rdat%fq5(2,i)*rdat%aqx2-rdat%fq4(2,i+4)* 3               )*rdat%aqxy
         fwk( 6,i)=-(rdat%fq5(3,i)*rdat%aqx2-rdat%fq4(1,i+4)*rdat%aqx2-rdat%fq4(3,i+4)* 3     &
     &                             +rdat%fq3(1,i+4)* 3               )*rdat%aqx
         fwk( 7,i)=-(rdat%fq5(1,i)*x2y2-rdat%fq4(1,i+4)* 3 *rdat%aqx2                &
     &              -rdat%fq4(1,i+4)*rdat%acy2+rdat%fq3(1,i+4)* 3             )*rdat%acy
         fwk( 8,i)=  rdat%fq5(2,i)*x2y2-rdat%fq4(2,i+4)*q2c2+rdat%fq3(2,i+4)
         fwk( 9,i)=-(rdat%fq5(3,i)*rdat%aqx2-rdat%fq4(1,i+4)*rdat%aqx2-rdat%fq4(3,i+4)        &
     &                             +rdat%fq3(1,i+4)                  )*rdat%acy
         fwk(10,i)=  rdat%fq5(4,i)*rdat%aqx2-rdat%fq4(2,i+4)* 3 *rdat%aqx2-rdat%fq4(4,i+4)    &
     &                             +rdat%fq3(2,i+4)* 3
         fwk(11,i)=-(rdat%fq5(1,i)*acy4-rdat%fq4(1,i+4)* 6 *rdat%acy2                &
     &                             +rdat%fq3(1,i+4)* 3               )*rdat%aqx
         fwk(12,i)= (rdat%fq5(2,i)*rdat%acy2-rdat%fq4(2,i+4)* 3               )*rdat%aqxy
         fwk(13,i)=-(rdat%fq5(3,i)*rdat%acy2-rdat%fq4(1,i+4)*rdat%acy2-rdat%fq4(3,i+4)        &
     &                             +rdat%fq3(1,i+4)                  )*rdat%aqx
         fwk(14,i)= (rdat%fq5(4,i)     -rdat%fq4(2,i+4)* 3               )*rdat%aqxy
         fwk(15,i)=-(rdat%fq5(5,i)     -rdat%fq4(3,i+4)* 6                      &
     &                             +rdat%fq3(1,i+4)* 3               )*rdat%aqx
         fwk(16,i)=-(rdat%fq5(1,i)*acy4-rdat%fq4(1,i+4)*10*rdat%acy2                &
     &                             +rdat%fq3(1,i+4)*15              )*rdat%acy
         fwk(17,i)=  rdat%fq5(2,i)*acy4-rdat%fq4(2,i+4)* 6 *rdat%acy2                &
     &                             +rdat%fq3(2,i+4)* 3
         fwk(18,i)=-(rdat%fq5(3,i)*rdat%acy2-rdat%fq4(1,i+4)*rdat%acy2-rdat%fq4(3,i+4)* 3     &
     &                             +rdat%fq3(1,i+4)* 3               )*rdat%acy
         fwk(19,i)=  rdat%fq5(4,i)*rdat%acy2-rdat%fq4(2,i+4)* 3 *rdat%acy2-rdat%fq4(4,i+4)    &
     &                             +rdat%fq3(2,i+4)* 3
         fwk(20,i)=-(rdat%fq5(5,i)-rdat%fq4(3,i+4)* 6 +rdat%fq3(1,i+4)* 3    )*rdat%acy
         fwk(21,i)=  rdat%fq5(6,i)-rdat%fq4(4,i+4)*10+rdat%fq3(2,i+4)*15
      enddo
      do i= 5, 6
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,5)*work(i+ 9)
         enddo
      enddo
      do i= 7, 8
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,6)*work(i+ 7)
         enddo
      enddo
      do i= 9,10
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,4)*work(i+ 5)
         enddo
      enddo
      do i=11,14
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,7)*work(i- 1)
         enddo
      enddo

      do j= 1, 7
         rdat%fq6(j,1)=rdat%fq6(j,1)*rdat%rab +rdat%fq6(j,3)
      enddo

      call frikr6(rdat, 1, 4,fw6,rdat%fq6, 3,rdat%fq5, 7,rdat%fq4, 7,rdat%fq3)

      do i= 1, 3
         do j= 1,28
            rdat%r06(j,i)= rdat%r06(j,i)+fw6(j,i)*fcc(j,6)
         end do
      end do
      i= 4
! Fw6( 1,i)= fw6( 1,i)*fcu( 1,6)
         fw6( 2,i)= fw6( 2,i)*fcu( 2,6)
         fw6( 3,i)= fw6( 3,i)*fcu( 3,6)
! Fw6( 4,i)= fw6( 4,i)*fcu( 4,6)
         fw6( 5,i)= fw6( 5,i)*fcu( 5,6)
! Fw6( 6,i)= fw6( 6,i)*fcu( 6,6)
         fw6( 7,i)= fw6( 7,i)*fcu( 7,6)
         fw6( 8,i)= fw6( 8,i)*fcu( 8,6)
         fw6( 9,i)= fw6( 9,i)*fcu( 9,6)
         fw6(10,i)= fw6(10,i)*fcu(10,6)
! Fw6(11,i)= fw6(11,i)*fcu(11,6)
         fw6(12,i)= fw6(12,i)*fcu(12,6)
! Fw6(13,i)= fw6(13,i)*fcu(13,6)
         fw6(14,i)= fw6(14,i)*fcu(14,6)
! Fw6(15,i)= fw6(15,i)*fcu(15,6)
         fw6(16,i)= fw6(16,i)*fcu(16,6)
         fw6(17,i)= fw6(17,i)*fcu(17,6)
         fw6(18,i)= fw6(18,i)*fcu(18,6)
         fw6(19,i)= fw6(19,i)*fcu(19,6)
         fw6(20,i)= fw6(20,i)*fcu(20,6)
         fw6(21,i)= fw6(21,i)*fcu(21,6)
! Fw6(22,i)= fw6(22,i)*fcu(22,6)
         fw6(23,i)= fw6(23,i)*fcu(23,6)
! Fw6(24,i)= fw6(24,i)*fcu(24,6)
         fw6(25,i)= fw6(25,i)*fcu(25,6)
! Fw6(26,i)= fw6(26,i)*fcu(26,6)
         fw6(27,i)= fw6(27,i)*fcu(27,6)
! Fw6(28,i)= fw6(28,i)*fcu(28,6)
      do i= 4, 5
         do j= 1,28
            rdat%r06(j,i)= rdat%r06(j,i)+fw6(j,4)*work(i+10)
         end do
      end do

      call frikr7(rdat, 1, 1,fw7,rdat%fq7, 3,rdat%fq6, 6,rdat%fq5,10,rdat%fq4)

      rdat%r07(:,1)= rdat%r07(:,1)+fw7(:)*fcc(1:36,7)

      end subroutine intk_20

! >
! >    @brief   dddd case
! >
! >    @details integration of the dddd case
! >
      subroutine intk_21(rdat,ikl)

      implicit none
      type(rotaxis_data_t) :: rdat
      integer :: ikl
      real(kind=dp) :: xmd2, xmd3, xmd4, xmd6
      real(kind=dp) :: xmdt, xmdtx, xmdtxy, xmdty
      real(kind=dp) :: aqx4, acy4, q2c2, x2y2
      real(kind=dp) :: y33, y34, y44
      real(kind=dp) :: work, fwk, fw6, fw7, fw8, fcu, fcc
      integer :: i, j

! Generate jtype=21 integrals

      dimension  work(15),fwk(36,16),fw6(28, 7),fw7(36, 3),fw8(45)
      dimension  fcu(45,8),fcc(45,8)

      if(ikl == 0) then
         rdat%r00= 0.0_dp
         rdat%r01= 0.0_dp
         rdat%r02= 0.0_dp
         rdat%r03= 0.0_dp
         rdat%r04= 0.0_dp
         rdat%r05= 0.0_dp
         rdat%r06= 0.0_dp
         rdat%r07= 0.0_dp
         rdat%r08= 0.0_dp
         return
      endif

      xmd2= rdat%x43 *0.5d+00
      xmd3= xmd2
      xmd4= xmd3*xmd2
      xmd6= xmd4*xmd2
      xmdt= xmd6*xmd2
      xmd2= xmd3

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

      y33 = rdat%y03 *rdat%y03
      y34 =-rdat%y03 *rdat%y04
      y44 = rdat%y04 *rdat%y04
      work( 1)= xmd4
      work( 2)= xmd2*y33
      work( 3)= xmd2*y34
      work( 4)= xmd2*y44
      work( 5)= y33 *y44

      work( 6)=-xmd4*rdat%y03
      work( 7)= xmd4*rdat%y04
      work( 8)= xmd2*y33*rdat%y04
      work( 9)= xmd2*y34*rdat%y04

      work(10)= xmd6
      work(11)= xmd4*y33
      work(12)= xmd4*y34
      work(13)= xmd4*y44

      work(14)=-xmd6*rdat%y03
      work(15)= xmd6*rdat%y04

      call fcufcc(rdat,8,xmdt,fcu,fcc)

      do i= 1, 5
         do j= 1, 5
            rdat%r00(j,i)= rdat%r00(j,i)+rdat%fq0(i)*work(j)
         end do
      end do

      do i= 1, 9
         fwk(1,i)=-rdat%fq1(1,i)*rdat%aqx
         fwk(2,i)=-rdat%fq1(1,i)*rdat%acy
         fwk(3,i)= rdat%fq1(2,i)
      enddo
      do i= 1, 4
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,1)*work(i+ 5)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,1)*work(i+ 5)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,1)*work(i+ 5)
      enddo
      do i= 5, 8
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,2)*work(i+ 1)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,2)*work(i+ 1)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,2)*work(i+ 1)
      enddo
      do i= 9,12
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,3)*work(i- 3)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,3)*work(i- 3)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,3)*work(i- 3)
      enddo
      do i=13,16
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,4)*work(i- 7)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,4)*work(i- 7)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,4)*work(i- 7)
      enddo
      do i=17,20
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,5)*work(i-11)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,5)*work(i-11)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,5)*work(i-11)
      enddo
      do i=21,25
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,6)*work(i-20)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,6)*work(i-20)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,6)*work(i-20)
      enddo
      do i=26,30
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,7)*work(i-25)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,7)*work(i-25)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,7)*work(i-25)
      enddo
      do i=31,35
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,8)*work(i-30)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,8)*work(i-30)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,8)*work(i-30)
      enddo
      do i=36,40
         rdat%r01(1,i)= rdat%r01(1,i)+fwk(1,9)*work(i-35)
         rdat%r01(2,i)= rdat%r01(2,i)+fwk(2,9)*work(i-35)
         rdat%r01(3,i)= rdat%r01(3,i)+fwk(3,9)*work(i-35)
      enddo

      do i= 1,13
         fwk(1,i)= rdat%fq2(1,i)*rdat%aqx2-rdat%fq1(1,i)
         fwk(2,i)= rdat%fq2(1,i)*rdat%acy2-rdat%fq1(1,i)
         fwk(3,i)= rdat%fq2(3,i)     -rdat%fq1(1,i)
         fwk(4,i)= rdat%fq2(1,i)*rdat%aqxy
         fwk(5,i)=-rdat%fq2(2,i)*rdat%aqx
         fwk(6,i)=-rdat%fq2(2,i)*rdat%acy
      enddo
      do i= 1, 4
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,1)*work(i+ 9)
         enddo
      enddo
      do i= 5, 8
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,2)*work(i+ 5)
         enddo
      enddo
      do i= 9,12
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,3)*work(i+ 1)
         enddo
      enddo
      do i=13,16
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,4)*work(i- 3)
         enddo
      enddo
      do i=17,20
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,5)*work(i- 7)
         enddo
      enddo
      do i=21,24
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,6)*work(i-15)
         enddo
      enddo
      do i=25,28
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,7)*work(i-19)
         enddo
      enddo
      do i=29,32
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,8)*work(i-23)
         enddo
      enddo
      do i=33,36
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,9)*work(i-27)
         enddo
      enddo
      do i=37,41
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,10)*work(i-36)
         enddo
      enddo
      do i=42,46
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,11)*work(i-41)
         enddo
      enddo
      do i=47,51
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,12)*work(i-46)
         enddo
      enddo
      do i=52,56
         do j= 1, 6
            rdat%r02(j,i)= rdat%r02(j,i)+fwk(j,13)*work(i-51)
         enddo
      enddo

      do i= 1,15
         fwk( 1,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)* 3 )*rdat%aqx
         fwk( 2,i)=-(rdat%fq3(1,i)*rdat%aqx2-rdat%fq2(1,i)    )*rdat%acy
         fwk( 3,i)=  rdat%fq3(2,i)*rdat%aqx2-rdat%fq2(2,i)
         fwk( 4,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)    )*rdat%aqx
         fwk( 5,i)=  rdat%fq3(2,i)                    *rdat%aqxy
         fwk( 6,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%aqx
         fwk( 7,i)=-(rdat%fq3(1,i)*rdat%acy2-rdat%fq2(1,i)* 3 )*rdat%acy
         fwk( 8,i)=  rdat%fq3(2,i)*rdat%acy2-rdat%fq2(2,i)
         fwk( 9,i)=-(rdat%fq3(3,i)     -rdat%fq2(1,i)    )*rdat%acy
         fwk(10,i)=  rdat%fq3(4,i)     -rdat%fq2(2,i)* 3
      enddo
      do i= 1, 2
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,1)*work(i+13)
         enddo
      enddo
      do i= 3, 4
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,2)*work(i+11)
         enddo
      enddo
      do i= 5, 6
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,3)*work(i+ 9)
         enddo
      enddo
      do i= 7, 8
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,4)*work(i+ 7)
         enddo
      enddo
      do i= 9,10
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,5)*work(i+ 5)
         enddo
      enddo
      do i=11,14
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,6)*work(i- 1)
         enddo
      enddo
      do i=15,18
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,7)*work(i- 5)
         enddo
      enddo
      do i=19,22
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,8)*work(i- 9)
         enddo
      enddo
      do i=23,26
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,9)*work(i-13)
         enddo
      enddo
      do i=27,30
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,10)*work(i-21)
         enddo
      enddo
      do i=31,34
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,11)*work(i-25)
         enddo
      enddo
      do i=35,38
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,12)*work(i-29)
         enddo
      enddo
      do i=39,42
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,13)*work(i-33)
         enddo
      enddo
      do i=43,47
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,14)*work(i-42)
         enddo
      enddo
      do i=48,52
         do j= 1,10
            rdat%r03(j,i)= rdat%r03(j,i)+fwk(j,15)*work(i-47)
         enddo
      enddo

      aqx4= rdat%aqx2*rdat%aqx2
      acy4= rdat%acy2*rdat%acy2
      x2y2= rdat%aqx2*rdat%acy2
      q2c2= rdat%aqx2+rdat%acy2
      do i= 1, 5
         rdat%r04( 1,i)= rdat%r04( 1,i)+(rdat%fq4(1,i)*aqx4-rdat%fq3(1,i)* 6 *rdat%aqx2        &
     &                                       +rdat%fq2(1,i)* 3    )*xmdt
         rdat%r04( 2,i)= rdat%r04( 2,i)+(rdat%fq4(1,i)*rdat%aqx2-rdat%fq3(1,i)* 3    )*xmdtxy
         rdat%r04( 3,i)= rdat%r04( 3,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)* 3    )*xmdtx
         rdat%r04( 4,i)= rdat%r04( 4,i)+(rdat%fq4(1,i)*x2y2-rdat%fq3(1,i)*q2c2            &
     &                                       +rdat%fq2(1,i)       )*xmdt
         rdat%r04( 5,i)= rdat%r04( 5,i)+(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)       )*xmdty
         rdat%r04( 6,i)= rdat%r04( 6,i)+(rdat%fq4(3,i)*rdat%aqx2-rdat%fq3(1,i)*rdat%aqx2            &
     &                             -rdat%fq3(3,i)+rdat%fq2(1,i)       )*xmdt
         rdat%r04( 7,i)= rdat%r04( 7,i)+(rdat%fq4(1,i)*rdat%acy2-rdat%fq3(1,i)* 3    )*xmdtxy
         rdat%r04( 8,i)= rdat%r04( 8,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)       )*xmdtx
         rdat%r04( 9,i)= rdat%r04( 9,i)+(rdat%fq4(3,i)     -rdat%fq3(1,i)       )*xmdtxy
         rdat%r04(10,i)= rdat%r04(10,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3    )*xmdtx
         rdat%r04(11,i)= rdat%r04(11,i)+(rdat%fq4(1,i)*acy4-rdat%fq3(1,i)* 6 *rdat%acy2        &
     &                                       +rdat%fq2(1,i)* 3    )*xmdt
         rdat%r04(12,i)= rdat%r04(12,i)+(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)* 3    )*xmdty
         rdat%r04(13,i)= rdat%r04(13,i)+(rdat%fq4(3,i)*rdat%acy2-rdat%fq3(1,i)*rdat%acy2            &
     &                             -rdat%fq3(3,i)+rdat%fq2(1,i)       )*xmdt
         rdat%r04(14,i)= rdat%r04(14,i)+(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3    )*xmdty
         rdat%r04(15,i)= rdat%r04(15,i)+(rdat%fq4(5,i)     -rdat%fq3(3,i)* 6              &
     &                                       +rdat%fq2(1,i)* 3    )*xmdt
      enddo
      do i= 6,16
         fwk( 1,i)=  rdat%fq4(1,i)*aqx4-rdat%fq3(1,i)* 6 *rdat%aqx2+rdat%fq2(1,i)* 3
         fwk( 2,i)= (rdat%fq4(1,i)*rdat%aqx2-rdat%fq3(1,i)* 3                )*rdat%aqxy
         fwk( 3,i)=-(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)* 3                )*rdat%aqx
         fwk( 4,i)=  rdat%fq4(1,i)*x2y2-rdat%fq3(1,i)*q2c2+rdat%fq2(1,i)
         fwk( 5,i)=-(rdat%fq4(2,i)*rdat%aqx2-rdat%fq3(2,i)                   )*rdat%acy
         fwk( 6,i)=  rdat%fq4(3,i)*rdat%aqx2-rdat%fq3(1,i)*rdat%aqx2-rdat%fq3(3,i)+rdat%fq2(1,i)
         fwk( 7,i)= (rdat%fq4(1,i)*rdat%acy2-rdat%fq3(1,i)* 3                )*rdat%aqxy
         fwk( 8,i)=-(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)                   )*rdat%aqx
         fwk( 9,i)= (rdat%fq4(3,i)     -rdat%fq3(1,i)                   )*rdat%aqxy
         fwk(10,i)=-(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3                )*rdat%aqx
         fwk(11,i)=  rdat%fq4(1,i)*acy4-rdat%fq3(1,i)* 6 *rdat%acy2+rdat%fq2(1,i)* 3
         fwk(12,i)=-(rdat%fq4(2,i)*rdat%acy2-rdat%fq3(2,i)* 3                )*rdat%acy
         fwk(13,i)=  rdat%fq4(3,i)*rdat%acy2-rdat%fq3(1,i)*rdat%acy2-rdat%fq3(3,i)+rdat%fq2(1,i)
         fwk(14,i)=-(rdat%fq4(4,i)     -rdat%fq3(2,i)* 3                )*rdat%acy
         fwk(15,i)=  rdat%fq4(5,i)     -rdat%fq3(3,i)* 6 +rdat%fq2(1,i)* 3
      enddo
      do i= 6, 7
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,6)*work(i+ 8)
         enddo
      enddo
      do i= 8, 9
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,7)*work(i+ 6)
         enddo
      enddo
      do i=10,11
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,8)*work(i+ 4)
         enddo
      enddo
      do i=12,13
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,9)*work(i+ 2)
         enddo
      enddo
      do i=14,17
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,10)*work(i- 4)
         enddo
      enddo
      do i=18,21
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,11)*work(i- 8)
         enddo
      enddo
      do i=22,25
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,12)*work(i-12)
         enddo
      enddo
      do i=26,29
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,13)*work(i-16)
         enddo
      enddo
      do i=30,33
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,14)*work(i-24)
         enddo
      enddo
      do i=34,37
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,15)*work(i-28)
         enddo
      enddo
      do i=38,42
         do j= 1,15
            rdat%r04(j,i)= rdat%r04(j,i)+fwk(j,16)*work(i-37)
         enddo
      enddo

      do i= 1, 4
         rdat%r05( 1,i)= rdat%r05( 1,i)+(rdat%fq5(1,i)*aqx4-rdat%fq4(1,i+5)*10*rdat%aqx2      &
     &                                       +rdat%fq3(1,i+5)*15  )*xmdtx
         rdat%r05( 2,i)= rdat%r05( 2,i)+(rdat%fq5(1,i)*aqx4-rdat%fq4(1,i+5)* 6 *rdat%aqx2      &
     &                                       +rdat%fq3(1,i+5)* 3   )*xmdty
         rdat%r05( 3,i)= rdat%r05( 3,i)+(rdat%fq5(2,i)*aqx4-rdat%fq4(2,i+5)* 6 *rdat%aqx2      &
     &                                       +rdat%fq3(2,i+5)* 3   )*xmdt
         rdat%r05( 4,i)= rdat%r05( 4,i)+(rdat%fq5(1,i)*x2y2-rdat%fq4(1,i+5)* 3 *rdat%acy2      &
     &                        -rdat%fq4(1,i+5)*rdat%aqx2+rdat%fq3(1,i+5)* 3 )*xmdtx
         rdat%r05( 5,i)= rdat%r05( 5,i)+(rdat%fq5(2,i)*rdat%aqx2-rdat%fq4(2,i+5)* 3   )*xmdtxy
         rdat%r05( 6,i)= rdat%r05( 6,i)+(rdat%fq5(3,i)*rdat%aqx2-rdat%fq4(1,i+5)*rdat%aqx2          &
     &                        -rdat%fq4(3,i+5)* 3 +rdat%fq3(1,i+5)* 3  )*xmdtx
         rdat%r05( 7,i)= rdat%r05( 7,i)+(rdat%fq5(1,i)*x2y2-rdat%fq4(1,i+5)* 3 *rdat%aqx2      &
     &                        -rdat%fq4(1,i+5)*rdat%acy2+rdat%fq3(1,i+5)* 3 )*xmdty
         rdat%r05( 8,i)= rdat%r05( 8,i)+(rdat%fq5(2,i)*x2y2-rdat%fq4(2,i+5)*q2c2          &
     &                                       +rdat%fq3(2,i+5)      )*xmdt
         rdat%r05( 9,i)= rdat%r05( 9,i)+(rdat%fq5(3,i)*rdat%aqx2-rdat%fq4(1,i+5)*rdat%aqx2          &
     &                           -rdat%fq4(3,i+5)+rdat%fq3(1,i+5)      )*xmdty
         rdat%r05(10,i)= rdat%r05(10,i)+(rdat%fq5(4,i)*rdat%aqx2-rdat%fq4(2,i+5)* 3 *rdat%aqx2      &
     &                           -rdat%fq4(4,i+5)+rdat%fq3(2,i+5)* 3   )*xmdt
         rdat%r05(11,i)= rdat%r05(11,i)+(rdat%fq5(1,i)*acy4-rdat%fq4(1,i+5)* 6 *rdat%acy2      &
     &                                  +rdat%fq3(1,i+5)* 3        )*xmdtx
         rdat%r05(12,i)= rdat%r05(12,i)+(rdat%fq5(2,i)*rdat%acy2-rdat%fq4(2,i+5)* 3   )*xmdtxy
         rdat%r05(13,i)= rdat%r05(13,i)+(rdat%fq5(3,i)*rdat%acy2-rdat%fq4(1,i+5)*rdat%acy2          &
     &                           -rdat%fq4(3,i+5)+rdat%fq3(1,i+5)      )*xmdtx
         rdat%r05(14,i)= rdat%r05(14,i)+(rdat%fq5(4,i)-rdat%fq4(2,i+5)* 3        )*xmdtxy
         rdat%r05(15,i)= rdat%r05(15,i)+(rdat%fq5(5,i)-rdat%fq4(3,i+5)* 6                 &
     &                                  +rdat%fq3(1,i+5)* 3        )*xmdtx
         rdat%r05(16,i)= rdat%r05(16,i)+(rdat%fq5(1,i)*acy4-rdat%fq4(1,i+5)*10*rdat%acy2      &
     &                                  +rdat%fq3(1,i+5)*15       )*xmdty
         rdat%r05(17,i)= rdat%r05(17,i)+(rdat%fq5(2,i)*acy4-rdat%fq4(2,i+5)* 6 *rdat%acy2      &
     &                                  +rdat%fq3(2,i+5)* 3        )*xmdt
         rdat%r05(18,i)= rdat%r05(18,i)+(rdat%fq5(3,i)*rdat%acy2-rdat%fq4(1,i+5)*rdat%acy2          &
     &                        -rdat%fq4(3,i+5)* 3 +rdat%fq3(1,i+5)* 3  )*xmdty
         rdat%r05(19,i)= rdat%r05(19,i)+(rdat%fq5(4,i)*rdat%acy2-rdat%fq4(2,i+5)* 3 *rdat%acy2      &
     &                           -rdat%fq4(4,i+5)+rdat%fq3(2,i+5)* 3   )*xmdt
         rdat%r05(20,i)= rdat%r05(20,i)+(rdat%fq5(5,i)-rdat%fq4(3,i+5)* 6                 &
     &                                  +rdat%fq3(1,i+5)* 3        )*xmdty
         rdat%r05(21,i)= rdat%r05(21,i)+(rdat%fq5(6,i)-rdat%fq4(4,i+5)*10                &
     &                                  +rdat%fq3(2,i+5)*15       )*xmdt
      enddo
      do i= 5,11
         fwk( 1,i)=-(rdat%fq5(1,i)*aqx4-rdat%fq4(1,i+5)*10*rdat%aqx2                &
     &                             +rdat%fq3(1,i+5)*15              )*rdat%aqx
         fwk( 2,i)=-(rdat%fq5(1,i)*aqx4-rdat%fq4(1,i+5)* 6 *rdat%aqx2                &
     &                             +rdat%fq3(1,i+5)* 3               )*rdat%acy
         fwk( 3,i)=  rdat%fq5(2,i)*aqx4-rdat%fq4(2,i+5)* 6 *rdat%aqx2                &
     &                             +rdat%fq3(2,i+5)* 3
         fwk( 4,i)=-(rdat%fq5(1,i)*x2y2-rdat%fq4(1,i+5)*rdat%aqx2                    &
     &                      -rdat%fq4(1,i+5)* 3 *rdat%acy2+rdat%fq3(1,i+5)* 3 )*rdat%aqx
         fwk( 5,i)= (rdat%fq5(2,i)*rdat%aqx2-rdat%fq4(2,i+5)* 3               )*rdat%aqxy
         fwk( 6,i)=-(rdat%fq5(3,i)*rdat%aqx2-rdat%fq4(1,i+5)*rdat%aqx2-rdat%fq4(3,i+5)* 3     &
     &                             +rdat%fq3(1,i+5)* 3               )*rdat%aqx
         fwk( 7,i)=-(rdat%fq5(1,i)*x2y2-rdat%fq4(1,i+5)* 3 *rdat%aqx2                &
     &                    -rdat%fq4(1,i+5)*rdat%acy2+rdat%fq3(1,i+5)* 3       )*rdat%acy
         fwk( 8,i)=  rdat%fq5(2,i)*x2y2-rdat%fq4(2,i+5)*q2c2+rdat%fq3(2,i+5)
         fwk( 9,i)=-(rdat%fq5(3,i)*rdat%aqx2-rdat%fq4(1,i+5)*rdat%aqx2-rdat%fq4(3,i+5)        &
     &                             +rdat%fq3(1,i+5)                  )*rdat%acy
         fwk(10,i)=  rdat%fq5(4,i)*rdat%aqx2-rdat%fq4(2,i+5)* 3 *rdat%aqx2-rdat%fq4(4,i+5)    &
     &                             +rdat%fq3(2,i+5)* 3
         fwk(11,i)=-(rdat%fq5(1,i)*acy4-rdat%fq4(1,i+5)* 6 *rdat%acy2                &
     &                             +rdat%fq3(1,i+5)* 3               )*rdat%aqx
         fwk(12,i)= (rdat%fq5(2,i)*rdat%acy2-rdat%fq4(2,i+5)* 3               )*rdat%aqxy
         fwk(13,i)=-(rdat%fq5(3,i)*rdat%acy2-rdat%fq4(1,i+5)*rdat%acy2-rdat%fq4(3,i+5)        &
     &                             +rdat%fq3(1,i+5)                  )*rdat%aqx
         fwk(14,i)= (rdat%fq5(4,i)     -rdat%fq4(2,i+5)* 3               )*rdat%aqxy
         fwk(15,i)=-(rdat%fq5(5,i)     -rdat%fq4(3,i+5)* 6                      &
     &                             +rdat%fq3(1,i+5)* 3               )*rdat%aqx
         fwk(16,i)=-(rdat%fq5(1,i)*acy4-rdat%fq4(1,i+5)*10*rdat%acy2                &
     &                             +rdat%fq3(1,i+5)*15              )*rdat%acy
         fwk(17,i)=  rdat%fq5(2,i)*acy4-rdat%fq4(2,i+5)* 6 *rdat%acy2                &
     &                             +rdat%fq3(2,i+5)* 3
         fwk(18,i)=-(rdat%fq5(3,i)*rdat%acy2-rdat%fq4(1,i+5)*rdat%acy2-rdat%fq4(3,i+5)* 3     &
     &                             +rdat%fq3(1,i+5)* 3               )*rdat%acy
         fwk(19,i)=  rdat%fq5(4,i)*rdat%acy2-rdat%fq4(2,i+5)* 3 *rdat%acy2-rdat%fq4(4,i+5)    &
     &                             +rdat%fq3(2,i+5)* 3
         fwk(20,i)=-(rdat%fq5(5,i)-rdat%fq4(3,i+5)* 6 +rdat%fq3(1,i+5)* 3    )*rdat%acy
         fwk(21,i)=  rdat%fq5(6,i)-rdat%fq4(4,i+5)*10+rdat%fq3(2,i+5)*15
      enddo
      do i= 5, 6
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,5)*work(i+ 9)
         enddo
      enddo
      do i= 7, 8
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,6)*work(i+ 7)
         enddo
      enddo
      do i= 9,10
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,7)*work(i+ 5)
         enddo
      enddo
      do i=11,12
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,8)*work(i+ 3)
         enddo
      enddo
      do i=13,16
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,9)*work(i- 3)
         enddo
      enddo
      do i=17,20
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,10)*work(i- 7)
         enddo
      enddo
      do i=21,24
         do j= 1,21
            rdat%r05(j,i)= rdat%r05(j,i)+fwk(j,11)*work(i-15)
         enddo
      enddo

      call frikr6(rdat, 1, 7,fw6,rdat%fq6, 4,rdat%fq5, 9,rdat%fq4, 9,rdat%fq3)

      do i= 1, 4
         do j= 1,28
            rdat%r06(j,i)= rdat%r06(j,i)+fw6(j,i)*fcc(j,6)
         end do
      end do
      do i= 5, 7
! Fw6( 1,i)= fw6( 1,i)*fcu( 1,6)
         fw6( 2,i)= fw6( 2,i)*fcu( 2,6)
         fw6( 3,i)= fw6( 3,i)*fcu( 3,6)
! Fw6( 4,i)= fw6( 4,i)*fcu( 4,6)
         fw6( 5,i)= fw6( 5,i)*fcu( 5,6)
! Fw6( 6,i)= fw6( 6,i)*fcu( 6,6)
         fw6( 7,i)= fw6( 7,i)*fcu( 7,6)
         fw6( 8,i)= fw6( 8,i)*fcu( 8,6)
         fw6( 9,i)= fw6( 9,i)*fcu( 9,6)
         fw6(10,i)= fw6(10,i)*fcu(10,6)
! Fw6(11,i)= fw6(11,i)*fcu(11,6)
         fw6(12,i)= fw6(12,i)*fcu(12,6)
! Fw6(13,i)= fw6(13,i)*fcu(13,6)
         fw6(14,i)= fw6(14,i)*fcu(14,6)
! Fw6(15,i)= fw6(15,i)*fcu(15,6)
         fw6(16,i)= fw6(16,i)*fcu(16,6)
         fw6(17,i)= fw6(17,i)*fcu(17,6)
         fw6(18,i)= fw6(18,i)*fcu(18,6)
         fw6(19,i)= fw6(19,i)*fcu(19,6)
         fw6(20,i)= fw6(20,i)*fcu(20,6)
         fw6(21,i)= fw6(21,i)*fcu(21,6)
! Fw6(22,i)= fw6(22,i)*fcu(22,6)
         fw6(23,i)= fw6(23,i)*fcu(23,6)
! Fw6(24,i)= fw6(24,i)*fcu(24,6)
         fw6(25,i)= fw6(25,i)*fcu(25,6)
! Fw6(26,i)= fw6(26,i)*fcu(26,6)
         fw6(27,i)= fw6(27,i)*fcu(27,6)
! Fw6(28,i)= fw6(28,i)*fcu(28,6)
      enddo
      do i= 5, 6
         do j= 1,28
            rdat%r06(j,i)= rdat%r06(j,i)+fw6(j,5)*work(i+ 9)
         enddo
      enddo
      do i= 7, 8
         do j= 1,28
            rdat%r06(j,i)= rdat%r06(j,i)+fw6(j,6)*work(i+ 7)
         enddo
      enddo
      do i= 9,12
         do j= 1,28
            rdat%r06(j,i)= rdat%r06(j,i)+fw6(j,7)*work(i+ 1)
         enddo
      enddo

      call frikr7(rdat, 1, 3,fw7,rdat%fq7, 4,rdat%fq6, 8,rdat%fq5,13,rdat%fq4)

      do i= 1, 2
         do j= 1,36
            rdat%r07(j,i)= rdat%r07(j,i)+fw7(j,i)*fcc(j,7)
         end do
      end do
      i= 3
         fw7( 1,i)= fw7( 1,i)*fcu( 1,7)
         fw7( 2,i)= fw7( 2,i)*fcu( 2,7)
! Fw7( 3,i)= fw7( 3,i)*fcu( 3,7)
         fw7( 4,i)= fw7( 4,i)*fcu( 4,7)
         fw7( 5,i)= fw7( 5,i)*fcu( 5,7)
         fw7( 6,i)= fw7( 6,i)*fcu( 6,7)
         fw7( 7,i)= fw7( 7,i)*fcu( 7,7)
! Fw7( 8,i)= fw7( 8,i)*fcu( 8,7)
         fw7( 9,i)= fw7( 9,i)*fcu( 9,7)
! Fw7(10,i)= fw7(10,i)*fcu(10,7)
         fw7(11,i)= fw7(11,i)*fcu(11,7)
         fw7(12,i)= fw7(12,i)*fcu(12,7)
         fw7(13,i)= fw7(13,i)*fcu(13,7)
         fw7(14,i)= fw7(14,i)*fcu(14,7)
         fw7(15,i)= fw7(15,i)*fcu(15,7)
         fw7(16,i)= fw7(16,i)*fcu(16,7)
! Fw7(17,i)= fw7(17,i)*fcu(17,7)
         fw7(18,i)= fw7(18,i)*fcu(18,7)
! Fw7(19,i)= fw7(19,i)*fcu(19,7)
         fw7(20,i)= fw7(20,i)*fcu(20,7)
! Fw7(21,i)= fw7(21,i)*fcu(21,7)
         fw7(22,i)= fw7(22,i)*fcu(22,7)
         fw7(23,i)= fw7(23,i)*fcu(23,7)
         fw7(24,i)= fw7(24,i)*fcu(24,7)
         fw7(25,i)= fw7(25,i)*fcu(25,7)
         fw7(26,i)= fw7(26,i)*fcu(26,7)
         fw7(27,i)= fw7(27,i)*fcu(27,7)
         fw7(28,i)= fw7(28,i)*fcu(28,7)
         fw7(29,i)= fw7(29,i)*fcu(29,7)
! Fw7(30,i)= fw7(30,i)*fcu(30,7)
         fw7(31,i)= fw7(31,i)*fcu(31,7)
! Fw7(32,i)= fw7(32,i)*fcu(32,7)
         fw7(33,i)= fw7(33,i)*fcu(33,7)
! Fw7(34,i)= fw7(34,i)*fcu(34,7)
         fw7(35,i)= fw7(35,i)*fcu(35,7)
! Fw7(36,i)= fw7(36,i)*fcu(36,7)
      do i= 3, 4
         do j= 1,36
            rdat%r07(j,i)= rdat%r07(j,i)+fw7(j,3)*work(i+11)
         end do
      end do

      call frikr8(rdat, 1, 1,fw8,rdat%fq8, 2,rdat%fq7, 6,rdat%fq6,10,rdat%fq5,15,rdat%fq4)

      do j= 1,45
         rdat%r08(j)= rdat%r08(j)+fw8(j)*fcc(j,8)
      enddo

      end subroutine intk_21

! >
! >    @brief   psss case
! >
! >    @details integration of a psss case
! >
      subroutine mcdv_02 (f, rdat, qx, qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(3,1,1,*)
      real(kind=dp) :: qx, qz

      f (1, 1, 1, 1) = + rdat%r01 (1, 1) + rdat%r00 (2, 1) * qx
      f (2, 1, 1, 1) = + rdat%r01 (2, 1)
      f (3, 1, 1, 1) = + rdat%r01 (3, 1) + rdat%r00 (2, 1) * qz

      end subroutine mcdv_02

! >
! >    @brief   ppss case
! >
! >    @details integration of a ppss case
! >
      subroutine mcdv_03 (f, rdat, qx, qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(3,3,1,*)
      real(kind=dp) :: qx, qz

      f (1,1,1,1) = rdat%r02(1,1) + rdat%r00(4,1) + (rdat%r01(1,3) + rdat%r01(1,4) + rdat%r00(5,1) * qx) * qx
      f (2,1,1,1) = rdat%r02(4,1) + rdat%r01(2,4) * qx
      f (3,1,1,1) = rdat%r02(5,1) + rdat%r01(3,4) * qx + (rdat%r01(1,3) + rdat%r00 (5,1) * qx) * qz

      f (1,2,1,1) = rdat%r02(4,1) + rdat%r01(2,3) * qx
      f (2,2,1,1) = rdat%r02(2,1) + rdat%r00(4,1)
      f (3,2,1,1) = rdat%r02(6,1) + rdat%r01(2,3) * qz

      f (1,3,1,1) = rdat%r02(5,1) + rdat%r01(1,4) * qz + (rdat%r01(3,3) + rdat%r00(5,1) * qz) * qx
      f (2,3,1,1) = rdat%r02(6,1) + rdat%r01(2,4) * qz
      f (3,3,1,1) = rdat%r02(3,1) + rdat%r00(4,1) + (rdat%r01(3,3) + rdat%r01(3,4) + rdat%r00(5,1) * qz) * qz

      end subroutine mcdv_03

! >
! >    @brief   psps case
! >
! >    @details integration of a psps case
! >
      subroutine mcdv_04 (f, rdat, qx, qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(3,1,3,*)
      real(kind=dp) :: qx, qz

      f (1,1,1,1) = -rdat%r02(1,1) - rdat%r01(1,4) * qx
      f (2,1,1,1) = -rdat%r02(4,1)
      f (3,1,1,1) = -rdat%r02(5,1) - rdat%r01(1,4) * qz

      f (1,1,2,1) = -rdat%r02(4,1) - rdat%r01(2,4) * qx
      f (2,1,2,1) = -rdat%r02(2,1)
      f (3,1,2,1) = -rdat%r02(6,1) - rdat%r01(2,4) * qz

      f (1,1,3,1) = -rdat%r02(5,1) - rdat%r01(3,4) * qx + rdat%r01(1,3) + rdat%r00(4,1) * qx
      f (2,1,3,1) = -rdat%r02(6,1) + rdat%r01(2,3)
      f (3,1,3,1) = -rdat%r02(3,1) - rdat%r01(3,4) * qz + rdat%r01(3,3) + rdat%r00(4,1) * qz

      end subroutine mcdv_04

! >
! >    @brief   ppps case
! >
! >    @details integration of a ppps case
! >
      subroutine mcdv_05 (f, rdat, qx, qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(3,3,3,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 4, lx = 1
      real(kind=dp) :: b1 (5, kx, lx), b2 (4, 3, kx, lx), b3 (6, kx, lx)

      integer, parameter :: ind (6, kx, lx) = reshape(&
              [ 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 2, 4, 5, 7, 8, 9, 3, 5, 6, 8, 9, 10 ], &
              shape(ind))
      integer :: i, j, k, l, m

      do j = 1, 6
      b3(j,2,1) = -rdat%r03(j,1)
      enddo

      do j = 1, 3
      m = in6(j)
      b2(1,j,2,1) = -rdat%r02(m,3)
      b2(2,j,2,1) = -rdat%r02(m,4)
      b2(3,j,2,1) = -rdat%r02(m,5)
      b2(4,j,2,1) = -rdat%r02(m,6)
      enddo

      j = 1
      do i = 1, 5
      b1(i,2,1) = -rdat%r01(j,i+8)
      enddo

      do j = 1, 6
      k = ind(j,3,1)
      b3(j,3,1) = -rdat%r03(k,1)
      enddo

      do j = 1, 3
      k = ind(j,3,1)
      m = in6(k)
      b2(1,j,3,1) = -rdat%r02(m,3)
      b2(2,j,3,1) = -rdat%r02(m,4)
      b2(3,j,3,1) = -rdat%r02(m,5)
      b2(4,j,3,1) = -rdat%r02(m,6)
      enddo

      j = 1
      k = ind(j,3,1)
      do i = 1, 5
      b1(i,3,1) = -rdat%r01(k,i+8)
      enddo

      do j = 1, 6
      k = ind(j,4,1)
      m = in6(j)
      b3(j,4,1) = -rdat%r03(k,1) + rdat%r02(m,2)
      enddo

      do j = 1, 3
      k = ind(j,4,1)
      m = in6(k)
      b2(1,j,4,1) = -rdat%r02(m,3) + rdat%r01(j,5)
      b2(2,j,4,1) = -rdat%r02(m,4) + rdat%r01(j,6)
      b2(3,j,4,1) = -rdat%r02(m,5) + rdat%r01(j,7)
      b2(4,j,4,1) = -rdat%r02(m,6) + rdat%r01(j,8)
      enddo

      j = 1
      k = ind(j,4,1)
      do i = 1, 5
      b1(i,4,1) = -rdat%r01(k,i+8) + rdat%r00(i,2)
      enddo

      do l = 1, lx
      do k = 2, kx

      f (1,1,k-1,l) = b3(1,k,l) + b1 (5,k,l) + (b2(3,1,k,l) + b2(4,1,k,l) + b1(4,k,l)*qx) * qx
      f (2,1,k-1,l) = b3(2,k,l) + b2 (4,2,k,l) * qx
      f (3,1,k-1,l) = b3(3,k,l) + b2 (4,3,k,l) * qx + (b2(3,1,k,l) + b1(4,k,l)*qx) * qz

      f (1,2,k-1,l) = b3(2,k,l) + b2 (3,2,k,l) * qx
      f (2,2,k-1,l) = b3(4,k,l) + b1 (5,k,l)
      f (3,2,k-1,l) = b3(5,k,l) + b2 (3,2,k,l) * qz

      f (1,3,k-1,l) = b3(3,k,l) + b2 (4,1,k,l) * qz + (b2(3,3,k,l) + b1(4,k,l)*qz) * qx
      f (2,3,k-1,l) = b3(5,k,l) + b2 (4,2,k,l) * qz
      f (3,3,k-1,l) = b3(6,k,l) + b1 (5,k,l) + (b2(3,3,k,l) + b2(4,3,k,l) + b1(4,k,l)*qz) * qz

      enddo
      enddo

      end subroutine mcdv_05

! >
! >    @brief   pppp case
! >
! >    @details integration of a pppp case
! >
      subroutine mcdv_06 (f, rdat, qx, qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(3,3,3,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 4, lx = 4
      real(kind=dp) :: b1 (5, kx, lx), b2 (4, 3, kx, lx), b3 (6, kx, lx)

      integer, parameter :: ind (6, kx, lx) = reshape(&
              [ 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 2, 4, 5, 7, 8, 9, 3, 5, 6, &
                8, 9, 10, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 2, 4, 5, 7, 8, 9, &
                3, 5, 6, 8, 9, 10, 2, 4, 5, 7, 8, 9, 2, 4, 5, 7, 8, 9, 4, 7, 8, &
                11, 12, 13, 5, 8, 9, 12, 13, 14, 3, 5, 6, 8, 9, 10, 3, 5, 6, 8, &
                9, 10, 5, 8, 9, 12, 13, 14, 6, 9, 10, 13, 14, 15], &
              shape(ind))
      integer :: i, j, k, l, m

      do j = 1, 6
      m = in6 (j)
      b3 (j, 1, 1) = + rdat%r02 (m, 1)
      enddo

      do j = 1, 3
      b2 (1, j, 1, 1) = + rdat%r01( j,1)
      b2 (2, j, 1, 1) = + rdat%r01( j,2)
      b2 (3, j, 1, 1) = + rdat%r01( j,3)
      b2 (4, j, 1, 1) = + rdat%r01( j,4)
      enddo

      do i = 1, 5
      b1 (i, 1, 1) = + rdat%r00 (i,1)
      enddo

      do j = 1, 6
      b3 (j, 2, 1) = - rdat%r03 (j, 2)
      enddo

      do j = 1, 3
      m = in6 (j)
      b2 (1, j, 2, 1) = - rdat%r02( m,6)
      b2 (2, j, 2, 1) = - rdat%r02( m,7)
      b2 (3, j, 2, 1) = - rdat%r02( m,8)
      b2 (4, j, 2, 1) = - rdat%r02( m,9)
      enddo

      j = 1
      do i = 1, 5
      b1 (i, 2, 1) = - rdat%r01( j,i + 20)
      enddo

      do j = 1, 6
      k = ind (j, 3, 1)
      b3 (j, 3, 1) = - rdat%r03 (k, 2)
      enddo

      do j = 1, 3
      k = ind (j, 3, 1)
      m = in6 (k)
      b2 (1, j, 3, 1) = - rdat%r02( m,6)
      b2 (2, j, 3, 1) = - rdat%r02( m,7)
      b2 (3, j, 3, 1) = - rdat%r02( m,8)
      b2 (4, j, 3, 1) = - rdat%r02( m,9)
      enddo

      j = 1
      k = ind (j, 3, 1)
      do i = 1, 5
      b1 (i, 3, 1) = - rdat%r01( k,i + 20)
      enddo

      do j = 1, 6
      k = ind (j, 4, 1)
      m = in6 (j)
      b3 (j, 4, 1) = - rdat%r03 (k, 2) + rdat%r02 (m, 2)
      enddo

      do j = 1, 3
      k = ind (j, 4, 1)
      m = in6 (k)
      b2 (1, j, 4, 1) = - rdat%r02( m,6) + rdat%r01( j,5)
      b2 (2, j, 4, 1) = - rdat%r02( m,7) + rdat%r01( j,6)
      b2 (3, j, 4, 1) = - rdat%r02( m,8) + rdat%r01( j,7)
      b2 (4, j, 4, 1) = - rdat%r02( m,9) + rdat%r01( j,8)
      enddo

      j = 1
      k = ind (j, 4, 1)
      do i = 1, 5
      b1 (i, 4, 1) = - rdat%r01( k,i + 20) + rdat%r00 (i,2)
      enddo

      do j = 1, 6
      b3 (j, 1, 2) = - rdat%r03 (j, 3)
      enddo

      do j = 1, 3
      m = in6 (j)
      b2 (1, j, 1, 2) = - rdat%r02( m,10)
      b2 (2, j, 1, 2) = - rdat%r02( m,11)
      b2 (3, j, 1, 2) = - rdat%r02( m,12)
      b2 (4, j, 1, 2) = - rdat%r02( m,13)
      enddo

      j = 1
      do i = 1, 5
      b1 (i, 1, 2) = - rdat%r01( j,i + 25)
      enddo

      do j = 1, 6
      m = in6 (j)
      b3 (j, 2, 2) = + rdat%r04 (j,1) + rdat%r02 (m, 5)
      enddo

      do j = 1, 3
      b2 (1, j, 2, 2) = + rdat%r03( j,5) + rdat%r01( j,17)
      b2 (2, j, 2, 2) = + rdat%r03( j,6) + rdat%r01( j,18)
      b2 (3, j, 2, 2) = + rdat%r03( j,7) + rdat%r01( j,19)
      b2 (4, j, 2, 2) = + rdat%r03( j,8) + rdat%r01( j,20)
      enddo

      j = 1
      do i = 1, 5
      b1 (i, 2, 2) = + rdat%r02( j,i + 21) + rdat%r00 (i,5)
      enddo

      do j = 1, 6
      k = ind (j, 3, 2)
      b3 (j, 3, 2) = + rdat%r04 (k,1)
      enddo

      do j = 1, 3
      k = ind (j, 3, 2)
      b2 (1, j, 3, 2) = + rdat%r03( k,5)
      b2 (2, j, 3, 2) = + rdat%r03( k,6)
      b2 (3, j, 3, 2) = + rdat%r03( k,7)
      b2 (4, j, 3, 2) = + rdat%r03( k,8)
      enddo

      j = 1
      k = ind (j, 3, 2)
      m = in6 (k)
      do i = 1, 5
      b1 (i, 3, 2) = + rdat%r02( m,i + 21)
      enddo

      do j = 1, 6
      k = ind (j, 4, 2)
      b3 (j, 4, 2) = + rdat%r04 (k,1) - rdat%r03 (j, 4)
      enddo

      do j = 1, 3
      k = ind (j, 4, 2)
      m = in6 (j)
      b2 (1, j, 4, 2) = + rdat%r03( k,5) - rdat%r02( m,14)
      b2 (2, j, 4, 2) = + rdat%r03( k,6) - rdat%r02( m,15)
      b2 (3, j, 4, 2) = + rdat%r03( k,7) - rdat%r02( m,16)
      b2 (4, j, 4, 2) = + rdat%r03( k,8) - rdat%r02( m,17)
      enddo

      j = 1
      k = ind (j, 4, 2)
      m = in6 (k)
      do i = 1, 5
      b1 (i, 4, 2) = + rdat%r02( m,i + 21) - rdat%r01( j,i + 30)
      enddo

      do j = 1, 6
      k = ind (j, 1, 3)
      b3 (j, 1, 3) = - rdat%r03 (k, 3)
      enddo

      do j = 1, 3
      k = ind (j, 1, 3)
      m = in6 (k)
      b2 (1, j, 1, 3) = - rdat%r02( m,10)
      b2 (2, j, 1, 3) = - rdat%r02( m,11)
      b2 (3, j, 1, 3) = - rdat%r02( m,12)
      b2 (4, j, 1, 3) = - rdat%r02( m,13)
      enddo

      j = 1
      k = ind (j, 1, 3)
      do i = 1, 5
      b1 (i, 1, 3) = - rdat%r01( k,i + 25)
      enddo

      do j = 1, 6
      k = ind (j, 3, 3)
      m = in6 (j)
      b3 (j, 3, 3) = + rdat%r04 (k,1) + rdat%r02 (m, 5)
      enddo

      do j = 1, 3
      k = ind (j, 3, 3)
      b2 (1, j, 3, 3) = + rdat%r03( k,5) + rdat%r01( j,17)
      b2 (2, j, 3, 3) = + rdat%r03( k,6) + rdat%r01( j,18)
      b2 (3, j, 3, 3) = + rdat%r03( k,7) + rdat%r01( j,19)
      b2 (4, j, 3, 3) = + rdat%r03( k,8) + rdat%r01( j,20)
      enddo

      j = 1
      k = ind (j, 3, 3)
      m = in6 (k)
      do i = 1, 5
      b1 (i, 3, 3) = + rdat%r02( m,i + 21) + rdat%r00 (i,5)
      enddo

      b3 (1, 4, 3) = + rdat%r04 (5,1) - rdat%r03 (2, 4)

      b3 (2, 4, 3) = + rdat%r04 (8,1) - rdat%r03 (4, 4)
      b3 (3, 4, 3) = + rdat%r04 (9,1) - rdat%r03 (5, 4)

      b3 (4, 4, 3) = + rdat%r04 (12,1) - rdat%r03 (7, 4)
      b3 (5, 4, 3) = + rdat%r04 (13,1) - rdat%r03 (8, 4)
      b3 (6, 4, 3) = + rdat%r04 (14,1) - rdat%r03 (9, 4)

      b2 (1, 1, 4, 3) = + rdat%r03( 5,5) - rdat%r02( 4,14)
      b2 (2, 1, 4, 3) = + rdat%r03( 5,6) - rdat%r02( 4,15)
      b2 (3, 1, 4, 3) = + rdat%r03( 5,7) - rdat%r02( 4,16)
      b2 (4, 1, 4, 3) = + rdat%r03( 5,8) - rdat%r02( 4,17)

      b2 (1, 2, 4, 3) = + rdat%r03( 8,5) - rdat%r02( 2,14)
      b2 (2, 2, 4, 3) = + rdat%r03( 8,6) - rdat%r02( 2,15)
      b2 (3, 2, 4, 3) = + rdat%r03( 8,7) - rdat%r02( 2,16)
      b2 (4, 2, 4, 3) = + rdat%r03( 8,8) - rdat%r02( 2,17)
      b2 (1, 3, 4, 3) = + rdat%r03( 9,5) - rdat%r02( 6,14)
      b2 (2, 3, 4, 3) = + rdat%r03( 9,6) - rdat%r02( 6,15)
      b2 (3, 3, 4, 3) = + rdat%r03( 9,7) - rdat%r02( 6,16)
      b2 (4, 3, 4, 3) = + rdat%r03( 9,8) - rdat%r02( 6,17)

      do i = 1, 5
      b1 (i, 4, 3) = + rdat%r02( 6,i + 21) - rdat%r01( 2,i + 30)
      enddo

      do j = 1, 6
      k = ind (j, 1, 4)
      m = in6 (j)
      b3 (j, 1, 4) = - rdat%r03 (k, 3) + rdat%r02 (m, 3)
      enddo

      do j = 1, 3
      k = ind (j, 1, 4)
      m = in6 (k)
      b2 (1, j, 1, 4) = - rdat%r02( m,10) + rdat%r01( j,9)
      b2 (2, j, 1, 4) = - rdat%r02( m,11) + rdat%r01( j,10)
      b2 (3, j, 1, 4) = - rdat%r02( m,12) + rdat%r01( j,11)
      b2 (4, j, 1, 4) = - rdat%r02( m,13) + rdat%r01( j,12)
      enddo

      j = 1
      k = ind (j, 1, 4)
      do i = 1, 5
      b1 (i, 1, 4) = - rdat%r01( k,i + 25) + rdat%r00 (i,3)
      enddo

      do j = 1, 6
      k = ind (j, 2, 4)
      b3 (j, 2, 4) = + rdat%r04 (k,1) - rdat%r03 (j, 1)
      enddo

      do j = 1, 3
      k = ind (j, 2, 4)
      m = in6 (j)
      b2 (1, j, 2, 4) = + rdat%r03( k,5) - rdat%r02( m,18)
      b2 (2, j, 2, 4) = + rdat%r03( k,6) - rdat%r02( m,19)
      b2 (3, j, 2, 4) = + rdat%r03( k,7) - rdat%r02( m,20)
      b2 (4, j, 2, 4) = + rdat%r03( k,8) - rdat%r02( m,21)
      enddo

      j = 1
      k = ind (j, 2, 4)
      m = in6 (k)
      do i = 1, 5
      b1 (i, 2, 4) = + rdat%r02( m,i + 21) - rdat%r01( j,i + 35)
      enddo

      b3 (1, 3, 4) = + rdat%r04 (5,1) - rdat%r03 (2, 1)

      b3 (2, 3, 4) = + rdat%r04 (8,1) - rdat%r03 (4, 1)
      b3 (3, 3, 4) = + rdat%r04 (9,1) - rdat%r03 (5, 1)

      b3 (4, 3, 4) = + rdat%r04 (12,1) - rdat%r03 (7, 1)
      b3 (5, 3, 4) = + rdat%r04 (13,1) - rdat%r03 (8, 1)
      b3 (6, 3, 4) = + rdat%r04 (14,1) - rdat%r03 (9, 1)

      b2 (1, 1, 3, 4) = + rdat%r03( 5,5) - rdat%r02( 4,18)
      b2 (2, 1, 3, 4) = + rdat%r03( 5,6) - rdat%r02( 4,19)
      b2 (3, 1, 3, 4) = + rdat%r03( 5,7) - rdat%r02( 4,20)
      b2 (4, 1, 3, 4) = + rdat%r03( 5,8) - rdat%r02( 4,21)

      b2 (1, 2, 3, 4) = + rdat%r03( 8,5) - rdat%r02( 2,18)
      b2 (2, 2, 3, 4) = + rdat%r03( 8,6) - rdat%r02( 2,19)
      b2 (3, 2, 3, 4) = + rdat%r03( 8,7) - rdat%r02( 2,20)
      b2 (4, 2, 3, 4) = + rdat%r03( 8,8) - rdat%r02( 2,21)
      b2 (1, 3, 3, 4) = + rdat%r03( 9,5) - rdat%r02( 6,18)
      b2 (2, 3, 3, 4) = + rdat%r03( 9,6) - rdat%r02( 6,19)
      b2 (3, 3, 3, 4) = + rdat%r03( 9,7) - rdat%r02( 6,20)
      b2 (4, 3, 3, 4) = + rdat%r03( 9,8) - rdat%r02( 6,21)

      do i = 1, 5
      b1 (i, 3, 4) = + rdat%r02( 6,i + 21) - rdat%r01( 2,i + 35)
      enddo

      b3 (1, 4, 4) = rdat%r04 (6,1) - rdat%r03( 3,1) - rdat%r03( 3,4) + rdat%r02( 1,4)     &
      + rdat%r02( 1,5)

      b3 (2, 4, 4) = rdat%r04 (9,1) - rdat%r03( 5,1) - rdat%r03( 5,4) + rdat%r02( 4,4)     &
      + rdat%r02( 4,5)
      b3 (3, 4, 4) = rdat%r04 (10,1) - rdat%r03( 6,1) - rdat%r03( 6,4) + rdat%r02( 5,4)    &
      + rdat%r02( 5,5)

      b3 (4, 4, 4) = rdat%r04 (13,1) - rdat%r03( 8,1) - rdat%r03( 8,4) + rdat%r02( 2,4)    &
      + rdat%r02( 2,5)
      b3 (5, 4, 4) = rdat%r04 (14,1) - rdat%r03( 9,1) - rdat%r03( 9,4) + rdat%r02( 6,4)    &
      + rdat%r02( 6,5)
      b3 (6, 4, 4) = rdat%r04 (15,1) - rdat%r03( 10,1) - rdat%r03( 10,4) + rdat%r02( 3,4)  &
      + rdat%r02( 3,5)

      b2 (1, 1, 4, 4) = rdat%r03( 6,5) - rdat%r02( 5,14) - rdat%r02( 5,18) + &
      rdat%r01( 1,  13) + rdat%r01( 1,17)
      b2 (2, 1, 4, 4) = rdat%r03( 6,6) - rdat%r02( 5,15) - rdat%r02( 5,19) + &
      rdat%r01( 1,  14) + rdat%r01( 1,18)
      b2 (3, 1, 4, 4) = rdat%r03( 6,7) - rdat%r02( 5,16) - rdat%r02( 5,20) + &
      rdat%r01( 1,  15) + rdat%r01( 1,19)
      b2 (4, 1, 4, 4) = rdat%r03( 6,8) - rdat%r02( 5,17) - rdat%r02( 5,21) + &
      rdat%r01( 1,  16) + rdat%r01( 1,20)

      b2 (1, 2, 4, 4) = rdat%r03( 9,5) - rdat%r02( 6,14) - rdat%r02( 6,18) + &
      rdat%r01( 2,  13) + rdat%r01( 2,17)
      b2 (2, 2, 4, 4) = rdat%r03( 9,6) - rdat%r02( 6,15) - rdat%r02( 6,19) + &
      rdat%r01( 2,  14) + rdat%r01( 2,18)
      b2 (3, 2, 4, 4) = rdat%r03( 9,7) - rdat%r02( 6,16) - rdat%r02( 6,20) + &
      rdat%r01( 2,  15) + rdat%r01( 2,19)
      b2 (4, 2, 4, 4) = rdat%r03( 9,8) - rdat%r02( 6,17) - rdat%r02( 6,21) + &
      rdat%r01( 2,  16) + rdat%r01( 2,20)
      b2 (1, 3, 4, 4) = rdat%r03( 10,5) - rdat%r02( 3,14) - rdat%r02( 3,18) + &
      rdat%r01( 3, 13) + rdat%r01( 3,17)
      b2 (2, 3, 4, 4) = rdat%r03( 10,6) - rdat%r02( 3,15) - rdat%r02( 3,19) + &
      rdat%r01( 3, 14) + rdat%r01( 3,18)
      b2 (3, 3, 4, 4) = rdat%r03( 10,7) - rdat%r02( 3,16) - rdat%r02( 3,20) + &
      rdat%r01( 3, 15) + rdat%r01( 3,19)
      b2 (4, 3, 4, 4) = rdat%r03( 10,8) - rdat%r02( 3,17) - rdat%r02( 3,21) + &
      rdat%r01( 3, 16) + rdat%r01( 3,20)

      do i = 1, 5
      b1 (i, 4, 4) = + rdat%r02( 3,i + 21) - rdat%r01( 3,i + 30) - &
      rdat%r01( 3,i + 35) + rdat%r00 (i,4) + rdat%r00 (i,5)
      enddo

      do l = 2, lx
      do k = 2, kx

      if (k == 2.and.l == 3) cycle

      f(1,1,k-1,l-1) = b3(1,k,l) + b1(5,k,l) + (b2(3,1,k,l) + b2(4,1,k,l) + b1(4,k,l)*qx) * qx
      f(2,1,k-1,l-1) = b3(2,k,l) + b2(4,2,k,l)*qx
      f(3,1,k-1,l-1) = b3(3,k,l) + b2(4,3,k,l) * qx + (b2(3,1,k,l) + b1(4,k,l)*qx) * qz

      f(1,2,k-1,l-1) = b3(2,k,l) + b2(3,2,k,l) * qx
      f(2,2,k-1,l-1) = b3(4,k,l) + b1(5,k,l)
      f(3,2,k-1,l-1) = b3(5,k,l) + b2(3,2,k,l) * qz

      f(1,3,k-1,l-1) = b3(3,k,l) + b2(4,1,k,l) * qz + (b2(3,3,k,l) + b1(4,k,l)*qz) * qx
      f(2,3,k-1,l-1) = b3(5,k,l) + b2(4,2,k,l) * qz
      f(3,3,k-1,l-1) = b3(6,k,l) + b1(5,k,l) + (b2(3,3,k,l) + b2(4,3,k,l) + b1(4,k,l)*qz) * qz

      enddo
      enddo

      f (:,:,1,2) = f (:,:,2,1)

      end subroutine mcdv_06

! >
! >    @brief   dsss case
! >
! >    @details integration of a dsss case
! >             simplified calculation of f(i,j,k,l) for cases where
! >                 i = 1..6,  j = 1..1,  k = 1..kx,  and  l = 1..lx
! >             using auxiliary arrays c1, c2 and c3.
! >
      subroutine mcdv_07(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,1,1,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 1, lx = 1
      real(kind=dp) :: c1(  2,kx,lx),c2(  3,kx,lx),c3(  6,kx,lx)
      integer :: j, k, l, m

! Auxiliary arrays to simplify the formulation of f(i,j,k,l)
! Where:  i = 1..6,  j = 1..1,  k = 1  and  l = 1

      do j=1,6
         m= in6(j)
         c3(j,1,1)=+rdat%r02(m,1)
      enddo

      c2(1,1,1)=+rdat%r01(1,1)

      c2(2,1,1)=+rdat%r01(2,1)
      c2(3,1,1)=+rdat%r01(3,1)

      c1(1,1,1)=+rdat%r00(1,1)
      c1(2,1,1)=+rdat%r00(2,1)

! Do l=1,lx
! Do k=1,kx

         l=1
            k=1

            f(1,1,k,l)=+c3(  1,k,l)+c1(  1,k,l)                         &
     &               +(+c2(  1,k,l)+c2(  1,k,l)+c1(  2,k,l)*qx)*qx
            f(2,1,k,l)=+c3(  4,k,l)+c1(  1,k,l)
            f(3,1,k,l)=+c3(  6,k,l)+c1(  1,k,l)                         &
     &               +(+c2(  3,k,l)+c2(  3,k,l)+c1(  2,k,l)*qz)*qz
            f(4,1,k,l)=+c3(  2,k,l)+c2(  2,k,l)*qx
            f(5,1,k,l)=+c3(  3,k,l)+c2(  3,k,l)*qx                      &
     &               +(+c2(  1,k,l)+c1(  2,k,l)*qx)*qz
            f(6,1,k,l)=+c3(  5,k,l)+c2(  2,k,l)*qz

! Enddo
! Enddo

      end subroutine mcdv_07

! >
! >    @brief   dpss case
! >
! >    @details integration of a dpss case
! >
      subroutine mcdv_08(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,3,1,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 1, lx = 1
      real(kind=dp) :: d1(  5,kx,lx),d2(4,3,kx,lx),d3(3,6,kx,lx)
      real(kind=dp) :: d4( 10,kx,lx)

      real(kind=dp) :: xx, zz, xz, xxx, xxz, xzz, zzz
      real(kind=dp) :: qxd, qzd, xzd
      integer :: i, j, k, l, m

      do j=1,10
         d4(  j,1,1)=+rdat%r03(j,1)
      enddo

      do j=1,6
         m= in6(j)
         d3(1,j,1,1)=+rdat%r02(m,1)
         d3(2,j,1,1)=+rdat%r02(m,2)
         d3(3,j,1,1)=+rdat%r02(m,3)
      enddo

      do j=1,3
         d2(1,j,1,1)=+rdat%r01(j,1)
         d2(2,j,1,1)=+rdat%r01(j,2)
         d2(3,j,1,1)=+rdat%r01(j,3)
         d2(4,j,1,1)=+rdat%r01(j,4)
      enddo

      do i=1,5
         d1(i  ,1,1)=+rdat%r00(i,1)
      enddo

      xx= qx*qx
      zz= qz*qz
      xz= qx*qz
      xxx= xx*qx
      xxz= xx*qz
      xzz= zz*qx
      zzz= zz*qz

      qxd= qx+qx
      qzd= qz+qz

      xzd= xz+xz

      l=1
      k=1

      f(1,1,k,l) = d4(  1,k,l)+d2(2,1,k,l)* 3 +(+d3(2,1,k,l)* 2 +d3(3,1,k,l) +d1(  3,k,l)* 2 +d1(  4,k,l))*qx +(+d2(3,1,k,l)+d2(4,1,k,l)* 2 )*xx +d1(  5,k,l)*xxx
      f(2,1,k,l) = d4(  4,k,l)+d2(2,1,k,l) +(+d3(3,4,k,l)+d1(  4,k,l))*qx
      f(3,1,k,l) = d4(  6,k,l)+d2(2,1,k,l) +(+d3(3,6,k,l)+d1(  4,k,l))*qx +d3(2,3,k,l)*qzd +d2(4,3,k,l)*xzd +d2(3,1,k,l)*zz +d1(  5,k,l)*xzz
      f(4,1,k,l) = d4(  2,k,l)+d2(2,2,k,l) +(+d3(2,2,k,l)+d3(3,2,k,l))*qx +d2(4,2,k,l)*xx
      f(5,1,k,l) = d4(  3,k,l)+d2(2,3,k,l) +(+d3(2,3,k,l)+d3(3,3,k,l))*qx +(+d3(2,1,k,l)+d1(  3,k,l))*qz +d2(4,3,k,l)*xx +(+d2(3,1,k,l)+d2(4,1,k,l))*xz +d1(  5,k,l)*xxz
      f(6,1,k,l) = d4(  5,k,l) +d3(3,5,k,l)*qx +d3(2,2,k,l)*qz +d2(4,2,k,l)*xz

      f(1,2,k,l) = d4(  2,k,l)+d2(2,2,k,l) +d3(2,2,k,l)*qxd +d2(3,2,k,l)*xx
      f(2,2,k,l) = d4(  7,k,l)+d2(2,2,k,l)* 3
      f(3,2,k,l) = d4(  9,k,l)+d2(2,2,k,l) +d3(2,5,k,l)*qzd +d2(3,2,k,l)*zz
      f(4,2,k,l) = d4(  4,k,l)+d2(2,1,k,l) +(+d3(2,4,k,l)+d1(  3,k,l))*qx
      f(5,2,k,l) = d4(  5,k,l) +d3(2,5,k,l)*qx +d3(2,2,k,l)*qz +d2(3,2,k,l)*xz
      f(6,2,k,l) = d4(  8,k,l)+d2(2,3,k,l) +(+d3(2,4,k,l)+d1(  3,k,l))*qz

      f(1,3,k,l) = d4(  3,k,l)+d2(2,3,k,l) +d3(2,3,k,l)*qxd +(+d3(3,1,k,l)+d1(  4,k,l))*qz +d2(3,3,k,l)*xx +d2(4,1,k,l)*xzd +d1(  5,k,l)*xxz
      f(2,3,k,l) = d4(  8,k,l)+d2(2,3,k,l) +(+d3(3,4,k,l)+d1(  4,k,l))*qz
      f(3,3,k,l) = d4( 10,k,l)+d2(2,3,k,l)* 3 +(+d3(2,6,k,l)* 2 +d3(3,6,k,l) +d1(  3,k,l)* 2 +d1(  4,k,l))*qz +(+d2(3,3,k,l)+d2(4,3,k,l)* 2 )*zz +d1(  5,k,l)*zzz
      f(4,3,k,l) = d4(  5,k,l) +d3(2,5,k,l)*qx +d3(3,2,k,l)*qz +d2(4,2,k,l)*xz
      f(5,3,k,l) = d4(  6,k,l)+d2(2,1,k,l) +(+d3(2,6,k,l)+d1(  3,k,l))*qx +(+d3(2,3,k,l)+d3(3,3,k,l))*qz +(+d2(3,3,k,l)+d2(4,3,k,l))*xz +d2(4,1,k,l)*zz +d1(  5,k,l)*xzz
      f(6,3,k,l) = d4(  9,k,l)+d2(2,2,k,l) +(+d3(2,5,k,l)+d3(3,5,k,l))*qz +d2(4,2,k,l)*zz

      end subroutine mcdv_08

! >
! >    @brief   dsps case
! >
! >    @details integration of a dsps case
! >
      subroutine mcdv_09(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,1,3,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 4, lx = 1
      real(kind=dp) :: c1(  2,kx,lx),c2(  3,kx,lx),c3(  6,kx,lx)

      integer, parameter :: ind(6, kx, lx) = reshape(&
              [ 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 2, 4, 5, 7, 8, 9, 3, 5, 6, 8, 9, 10] &
              , shape(ind))
      integer :: j, k, l, m

      do j=1,6
         m= in6(j)
         c3(j,1,1)=+rdat%r02(m, 1)
      enddo

      c2(  1,1,1)=+rdat%r01( 1, 1)

      c2(  2,1,1)=+rdat%r01( 2, 1)
      c2(  3,1,1)=+rdat%r01( 3, 1)

      c1(  1,1,1)=+rdat%r00(1,1)
      c1(  2,1,1)=+rdat%r00(2,1)

      do j=1,6
         c3(j,2,1)=-rdat%r03(j, 1)
      enddo

      c2(  1,2,1)=-rdat%r02( 1, 3)

      c2(  2,2,1)=-rdat%r02( 4, 3)
      c2(  3,2,1)=-rdat%r02( 5, 3)

      c1(  1,2,1)=-rdat%r01( 1, 3)
      c1(  2,2,1)=-rdat%r01( 1, 4)

      do j=1,6
         k= ind(j,3,1)
         c3(j,3,1)=-rdat%r03(k, 1)
      enddo

      c2(  1,3,1)=-rdat%r02( 4, 3)

      c2(  2,3,1)=-rdat%r02( 2, 3)
      c2(  3,3,1)=-rdat%r02( 6, 3)

      c1(  1,3,1)=-rdat%r01( 2, 3)
      c1(  2,3,1)=-rdat%r01( 2, 4)

      do j=1,6
         k= ind(j,4,1)
         m= in6(j)
         c3(j,4,1)=-rdat%r03(k, 1)+rdat%r02( m, 2)
      enddo

      c2(1,4,1)=-rdat%r02(5,3)+rdat%r01(1,2)

      c2(2,4,1)=-rdat%r02(6,3)+rdat%r01(2,2)
      c2(3,4,1)=-rdat%r02(3,3)+rdat%r01(3,2)

      c1(1,4,1)=-rdat%r01(3,3)+rdat%r00(3,1)
      c1(2,4,1)=-rdat%r01(3,4)+rdat%r00(4,1)

      l=1
      do k=2,kx
        f(1,1,k-1,l) = c3(1,k,l)+c1(1,k,l) + (c2(1,k,l)+c2(1,k,l)+c1(2,k,l)*qx)*qx
        f(2,1,k-1,l) = c3(4,k,l)+c1(1,k,l)
        f(3,1,k-1,l) = c3(6,k,l)+c1(1,k,l) + (c2(3,k,l)+c2(3,k,l)+c1(2,k,l)*qz)*qz
        f(4,1,k-1,l) = c3(2,k,l)+c2(2,k,l)*qx
        f(5,1,k-1,l) = c3(3,k,l)+c2(3,k,l)*qx +(c2(1,k,l)+c1(2,k,l)*qx)*qz
        f(6,1,k-1,l) = c3(5,k,l)+c2(2,k,l)*qz
      enddo

      end subroutine mcdv_09

! >
! >    @brief   ddss case
! >
! >    @details integration of a ddss case
! >
      subroutine mcdv_10(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,6,1,*)
      real(kind=dp) :: qx, qz

      !integer, parameter :: kx = 4, lx = 1
      integer, parameter :: kx = 1, lx = 1
      real(kind=dp) :: e1(  5,kx,lx),e2(4,3,kx,lx),e3(4,6,kx,lx), &
                       e4(2,10,kx,lx),e5( 15,kx,lx)

      integer :: i, j, k, l, m
      real(kind=dp) :: xx, zz, xz, xxx, xxz, xzz, zzz
      real(kind=dp) :: xxxx, xxzz, xxxz, xzzz, zzzz
      real(kind=dp) :: xxxd, xxzd, xzzd, zzzd
      real(kind=dp) :: xzq
      real(kind=dp) :: qxd, qzd, xzd

      do j=1,15
         e5(  j,1,1)=+rdat%r04(j,1)
      enddo

      do j=1,10
         e4(1,j,1,1)=+rdat%r03(j,1)
         e4(2,j,1,1)=+rdat%r03(j,2)
      enddo

      do j=1,6
         m= in6(j)
         e3(1,j,1,1)=+rdat%r02(m,1)
         e3(2,j,1,1)=+rdat%r02(m,2)
         e3(3,j,1,1)=+rdat%r02(m,3)
         e3(4,j,1,1)=+rdat%r02(m,4)
      enddo

      do j=1,3
         e2(1,j,1,1)=+rdat%r01(j,1)
         e2(2,j,1,1)=+rdat%r01(j,2)
         e2(3,j,1,1)=+rdat%r01(j,3)
         e2(4,j,1,1)=+rdat%r01(j,4)
      enddo

      do i=1,5
         e1(i  ,1,1)=+rdat%r00(i,1)
      enddo

      xx= qx*qx
      zz= qz*qz
      xz= qx*qz
      xxx= xx*qx
      xxz= xx*qz
      xzz= zz*qx
      zzz= zz*qz
      xxxx= xx*xx
      xxxz= xx*xz
      xxzz= xx*zz
      xzzz= xz*zz
      zzzz= zz*zz

      qxd= qx+qx
      qzd= qz+qz
      xzd= xz+xz
      xxxd= xxx+xxx
      xxzd= xxz+xxz
      xzzd= xzz+xzz
      zzzd= zzz+zzz

      xzq= xz+xz+xz+xz

      f(1,1,1,1) = e5(  1,1,1)+e3(1,1,1,1)* 6 +e1(  1,1,1)* 3 +(+e4(1,1,1,1)+e4(2,1,1,1)+ (+e2(1,1,1,1)+e2(2,1,1,1))* 3 )*qxd +(+e3(2,1,1,1)+e3(3,1,1,1)* 4 +e3(4,1,1,1) +e1(  2,1,1)+e1(  3,1,1)* 4 +e1(  4,1,1))*xx +(+e2(3,1,1,1)+e2(4,1,1,1))*xxxd +e1(  5,1,1)*xxxx
      f(2,1,1,1) = e5(  4,1,1)+e3(1,4,1,1)+e3(1,1,1,1)+e1(  1,1,1) +(+e4(2,4,1,1)+e2(2,1,1,1))*qxd +(+e3(4,4,1,1)+e1(  4,1,1))*xx
      f(3,1,1,1) = e5(  6,1,1)+e3(1,6,1,1)+e3(1,1,1,1)+e1(  1,1,1) +(+e4(2,6,1,1)+e2(2,1,1,1))*qxd +(+e4(1,3,1,1)+e2(1,3,1,1))*qzd +(+e3(4,6,1,1)+e1(  4,1,1))*xx +e3(3,3,1,1)*xzq +(+e3(2,1,1,1)+e1(  2,1,1))*zz +e2(4,3,1,1)*xxzd +e2(3,1,1,1)*xzzd +e1(  5,1,1)*xxzz
      f(4,1,1,1) = e5(  2,1,1)+e3(1,2,1,1)* 3 +(+e4(1,2,1,1)+e4(2,2,1,1)* 2 +e2(1,2,1,1)+e2(2,2,1,1)* 2 )*qx +(+e3(3,2,1,1)* 2 +e3(4,2,1,1))*xx +e2(4,2,1,1)*xxx
      f(5,1,1,1) = e5(  3,1,1)+e3(1,3,1,1)* 3 +(+e4(1,3,1,1)+e4(2,3,1,1)* 2 +e2(1,3,1,1)+e2(2,3,1,1)* 2 )*qx +(+e4(1,1,1,1)+e2(1,1,1,1)* 3 )*qz +(+e3(3,3,1,1)* 2 +e3(4,3,1,1))*xx +(+e3(2,1,1,1)+e3(3,1,1,1)* 2 +e1(  2,1,1)+e1(  3,1,1)* 2 )*xz +e2(4,3,1,1)*xxx +(+e2(3,1,1,1)* 2 +e2(4,1,1,1))*xxz +e1(  5,1,1)*xxxz
      f(6,1,1,1) = e5(  5,1,1)+e3(1,5,1,1) +e4(2,5,1,1)*qxd +(+e4(1,2,1,1)+e2(1,2,1,1))*qz +e3(4,5,1,1)*xx +e3(3,2,1,1)*xzd +e2(4,2,1,1)*xxz

      f(1,2,1,1) = e5(  4,1,1)+e3(1,4,1,1)+e3(1,1,1,1)+e1(  1,1,1) +(+e4(1,4,1,1)+e2(1,1,1,1))*qxd +(+e3(2,4,1,1)+e1(  2,1,1))*xx
      f(2,2,1,1) = e5( 11,1,1)+e3(1,4,1,1)* 6 +e1(  1,1,1)* 3
      f(3,2,1,1) = e5( 13,1,1)+e3(1,6,1,1)+e3(1,4,1,1)+e1(  1,1,1) +(+e4(1,8,1,1)+e2(1,3,1,1))*qzd +(+e3(2,4,1,1)+e1(  2,1,1))*zz
      f(4,2,1,1) = e5(  7,1,1)+e3(1,2,1,1)* 3 +(+e4(1,7,1,1)+e2(1,2,1,1)* 3 )*qx
      f(5,2,1,1) = e5(  8,1,1)+e3(1,3,1,1) +(+e4(1,8,1,1)+e2(1,3,1,1))*qx +(+e4(1,4,1,1)+e2(1,1,1,1))*qz +(+e3(2,4,1,1)+e1(  2,1,1))*xz
      f(6,2,1,1) = e5( 12,1,1)+e3(1,5,1,1)* 3 +(+e4(1,7,1,1)+e2(1,2,1,1)* 3 )*qz

      f(1,3,1,1) = e5(  6,1,1)+e3(1,6,1,1)+e3(1,1,1,1)+e1(  1,1,1) +(+e4(1,6,1,1)+e2(1,1,1,1))*qxd +(+e4(2,3,1,1)+e2(2,3,1,1))*qzd +(+e3(2,6,1,1)+e1(  2,1,1))*xx +e3(3,3,1,1)*xzq +(+e3(4,1,1,1)+e1(  4,1,1))*zz +e2(3,3,1,1)*xxzd +e2(4,1,1,1)*xzzd +e1(  5,1,1)*xxzz
      f(2,3,1,1) = e5( 13,1,1)+e3(1,6,1,1)+e3(1,4,1,1)+e1(  1,1,1) +(+e4(2,8,1,1)+e2(2,3,1,1))*qzd +(+e3(4,4,1,1)+e1(  4,1,1))*zz
      f(3,3,1,1) = e5( 15,1,1)+e3(1,6,1,1)* 6 +e1(  1,1,1)* 3 +(+e4(1,10,1,1)+e4(2,10,1,1)+ (+e2(1,3,1,1)+e2(2,3,1,1))* 3 )*qzd +(+e3(2,6,1,1)+e3(3,6,1,1)* 4 +e3(4,6,1,1) +e1(  2,1,1)+e1(  3,1,1)* 4 +e1(  4,1,1))*zz +(+e2(3,3,1,1)+e2(4,3,1,1))*zzzd +e1(  5,1,1)*zzzz
      f(4,3,1,1) = e5(  9,1,1)+e3(1,2,1,1) +(+e4(1,9,1,1)+e2(1,2,1,1))*qx +e4(2,5,1,1)*qzd +e3(3,5,1,1)*xzd +e3(4,2,1,1)*zz +e2(4,2,1,1)*xzz
      f(5,3,1,1) = e5( 10,1,1)+e3(1,3,1,1)* 3 +(+e4(1,10,1,1)+e2(1,3,1,1)* 3 )*qx +(+e4(1,6,1,1)+e4(2,6,1,1)* 2 +e2(1,1,1,1)+e2(2,1,1,1)* 2 )*qz +(+e3(2,6,1,1)+e3(3,6,1,1)* 2 +e1(  2,1,1)+e1(  3,1,1)* 2 )*xz +(+e3(3,3,1,1)* 2 +e3(4,3,1,1))*zz +(+e2(3,3,1,1)* 2 +e2(4,3,1,1))*xzz +e2(4,1,1,1)*zzz +e1(  5,1,1)*xzzz
      f(6,3,1,1) = e5( 14,1,1)+e3(1,5,1,1)* 3 +(+e4(1,9,1,1)+e4(2,9,1,1)* 2 +e2(1,2,1,1)+e2(2,2,1,1)* 2 )*qz +(+e3(3,5,1,1)* 2 +e3(4,5,1,1))*zz +e2(4,2,1,1)*zzz

      f(1,4,1,1) = e5(  2,1,1)+e3(1,2,1,1)* 3 +(+e4(1,2,1,1)* 2 +e4(2,2,1,1) +e2(1,2,1,1)* 2 +e2(2,2,1,1))*qx +(+e3(2,2,1,1)+e3(3,2,1,1)* 2 )*xx +e2(3,2,1,1)*xxx
      f(2,4,1,1) = e5(  7,1,1)+e3(1,2,1,1)* 3 +(+e4(2,7,1,1)+e2(2,2,1,1)* 3 )*qx
      f(3,4,1,1) = e5(  9,1,1)+e3(1,2,1,1) +(+e4(2,9,1,1)+e2(2,2,1,1))*qx +e4(1,5,1,1)*qzd +e3(3,5,1,1)*xzd +e3(2,2,1,1)*zz +e2(3,2,1,1)*xzz
      f(4,4,1,1) = e5(  4,1,1)+e3(1,4,1,1)+e3(1,1,1,1)+e1(  1,1,1) +(+e4(1,4,1,1)+e4(2,4,1,1) +e2(1,1,1,1)+e2(2,1,1,1))*qx +(+e3(3,4,1,1)+e1(  3,1,1))*xx
      f(5,4,1,1) = e5(  5,1,1)+e3(1,5,1,1) +(+e4(1,5,1,1)+e4(2,5,1,1))*qx +(+e4(1,2,1,1)+e2(1,2,1,1))*qz +e3(3,5,1,1)*xx +(+e3(2,2,1,1)+e3(3,2,1,1))*xz +e2(3,2,1,1)*xxz
      f(6,4,1,1) = e5(  8,1,1)+e3(1,3,1,1) +(+e4(2,8,1,1)+e2(2,3,1,1))*qx +(+e4(1,4,1,1)+e2(1,1,1,1))*qz +(+e3(3,4,1,1)+e1(  3,1,1))*xz

      f(1,5,1,1) = e5(  3,1,1)+e3(1,3,1,1)* 3 +(+e4(1,3,1,1)* 2 +e4(2,3,1,1) +e2(1,3,1,1)* 2 +e2(2,3,1,1))*qx +(+e4(2,1,1,1)+e2(2,1,1,1)* 3 )*qz +(+e3(2,3,1,1)+e3(3,3,1,1)* 2 )*xx +(+e3(3,1,1,1)* 2 +e3(4,1,1,1) +e1(  3,1,1)* 2 +e1(  4,1,1))*xz +e2(3,3,1,1)*xxx +(+e2(3,1,1,1)+e2(4,1,1,1)* 2 )*xxz +e1(  5,1,1)*xxxz
      f(2,5,1,1) = e5(  8,1,1)+e3(1,3,1,1) +(+e4(2,8,1,1)+e2(2,3,1,1))*qx +(+e4(2,4,1,1)+e2(2,1,1,1))*qz +(+e3(4,4,1,1)+e1(  4,1,1))*xz
      f(3,5,1,1) = e5( 10,1,1)+e3(1,3,1,1)* 3 +(+e4(2,10,1,1)+e2(2,3,1,1)* 3 )*qx +(+e4(1,6,1,1)* 2 +e4(2,6,1,1) +e2(1,1,1,1)* 2 +e2(2,1,1,1))*qz +(+e3(3,6,1,1)* 2 +e3(4,6,1,1) +e1(  3,1,1)* 2 +e1(  4,1,1))*xz +(+e3(2,3,1,1)+e3(3,3,1,1)* 2 )*zz +(+e2(3,3,1,1)+e2(4,3,1,1)* 2 )*xzz +e2(3,1,1,1)*zzz +e1(  5,1,1)*xzzz
      f(4,5,1,1) = e5(  5,1,1)+e3(1,5,1,1) +(+e4(1,5,1,1)+e4(2,5,1,1))*qx +(+e4(2,2,1,1)+e2(2,2,1,1))*qz +e3(3,5,1,1)*xx +(+e3(3,2,1,1)+e3(4,2,1,1))*xz +e2(4,2,1,1)*xxz
      f(5,5,1,1) = e5(  6,1,1)+e3(1,6,1,1)+e3(1,1,1,1)+e1(  1,1,1) +(+e4(1,6,1,1)+e4(2,6,1,1) +e2(1,1,1,1)+e2(2,1,1,1))*qx +(+e4(1,3,1,1)+e4(2,3,1,1) +e2(1,3,1,1)+e2(2,3,1,1))*qz +(+e3(3,6,1,1)+e1(  3,1,1))*xx +(+e3(2,3,1,1)+e3(3,3,1,1)* 2 +e3(4,3,1,1))*xz +(+e3(3,1,1,1)+e1(  3,1,1))*zz +(+e2(3,3,1,1)+e2(4,3,1,1))*xxz +(+e2(3,1,1,1)+e2(4,1,1,1))*xzz +e1(  5,1,1)*xxzz
      f(6,5,1,1) = e5(  9,1,1)+e3(1,2,1,1) +(+e4(2,9,1,1)+e2(2,2,1,1))*qx +(+e4(1,5,1,1)+e4(2,5,1,1))*qz +(+e3(3,5,1,1)+e3(4,5,1,1))*xz +e3(3,2,1,1)*zz +e2(4,2,1,1)*xzz

      f(1,6,1,1) = e5(  5,1,1)+e3(1,5,1,1) +e4(1,5,1,1)*qxd +(+e4(2,2,1,1)+e2(2,2,1,1))*qz +e3(2,5,1,1)*xx +e3(3,2,1,1)*xzd +e2(3,2,1,1)*xxz
      f(2,6,1,1) = e5( 12,1,1)+e3(1,5,1,1)* 3 +(+e4(2,7,1,1)+e2(2,2,1,1)* 3 )*qz
      f(3,6,1,1) = e5( 14,1,1)+e3(1,5,1,1)* 3 +(+e4(1,9,1,1)* 2 +e4(2,9,1,1) +e2(1,2,1,1)* 2 +e2(2,2,1,1))*qz +(+e3(2,5,1,1)+e3(3,5,1,1)* 2 )*zz +e2(3,2,1,1)*zzz
      f(4,6,1,1) = e5(  8,1,1)+e3(1,3,1,1) +(+e4(1,8,1,1)+e2(1,3,1,1))*qx +(+e4(2,4,1,1)+e2(2,1,1,1))*qz +(+e3(3,4,1,1)+e1(  3,1,1))*xz
      f(5,6,1,1) = e5(  9,1,1)+e3(1,2,1,1) +(+e4(1,9,1,1)+e2(1,2,1,1))*qx +(+e4(1,5,1,1)+e4(2,5,1,1))*qz +(+e3(2,5,1,1)+e3(3,5,1,1))*xz +e3(3,2,1,1)*zz +e2(3,2,1,1)*xzz
      f(6,6,1,1) = e5( 13,1,1)+e3(1,6,1,1)+e3(1,4,1,1)+e1(  1,1,1) +(+e4(1,8,1,1)+e4(2,8,1,1) +e2(1,3,1,1)+e2(2,3,1,1))*qz +(+e3(3,4,1,1)+e1(  3,1,1))*zz

      end subroutine mcdv_10

! >
! >    @brief   dpps case
! >
! >    @details integration of a dpps case
! >
      subroutine mcdv_11(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,3,3,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 4, lx = 1
      real(kind=dp) :: d1(  5,kx,lx),d2(4,3,kx,lx),d3(3,6,kx,lx), &
                       d4( 10,kx,lx)

      integer, parameter :: ind (10, kx, lx) = reshape(&
              [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, &
                4, 5, 7, 8, 9, 11, 12, 13, 14, 3, 5, 6, 8, 9, 10, 12, 13, 14, 15] &
              , shape(ind))

      real(kind=dp) :: xx, zz, xz, xxx, xxz, xzz, zzz
      real(kind=dp) :: qxd, qzd, xzd
      integer :: i, j, k, l, m


      do j=1,10
         d4(  j,1,1)=+rdat%r03( j,1)
      enddo

      do j=1,6
         m= in6(j)
         d3(1,j,1,1)=+rdat%r02(m, 1)
         d3(2,j,1,1)=+rdat%r02(m, 2)
         d3(3,j,1,1)=+rdat%r02(m, 3)
      enddo

      do j=1,3
         d2(1,j,1,1)=+rdat%r01(j, 1)
         d2(2,j,1,1)=+rdat%r01(j, 2)
         d2(3,j,1,1)=+rdat%r01(j, 3)
         d2(4,j,1,1)=+rdat%r01(j, 4)
      enddo

      do i=1,5
         d1(i  ,1,1)=+rdat%r00(i,1)
      enddo

      do j=1,10
         d4(  j,2,1)=-rdat%r04( j,1)
      enddo

      do j=1,6
         d3(1,j,2,1)=-rdat%r03(j, 3)
         d3(2,j,2,1)=-rdat%r03(j, 4)
         d3(3,j,2,1)=-rdat%r03(j, 5)
      enddo

      do j=1,3
         m= in6(j)
         d2(1,j,2,1)=-rdat%r02(m, 7)
         d2(2,j,2,1)=-rdat%r02(m, 8)
         d2(3,j,2,1)=-rdat%r02(m, 9)
         d2(4,j,2,1)=-rdat%r02(m,10)
      enddo

         j=1
      do i=1,5
         d1(i  ,2,1)=-rdat%r01(j,i+ 8)
      enddo

      do j=1,10
         k= ind(j,3,1)
         d4(  j,3,1)=-rdat%r04( k,1)
      enddo

      do j=1,6
         k= ind(j,3,1)
         d3(1,j,3,1)=-rdat%r03(k, 3)
         d3(2,j,3,1)=-rdat%r03(k, 4)
         d3(3,j,3,1)=-rdat%r03(k, 5)
      enddo

      do j=1,3
         k= ind(j,3,1)
         m= in6(k)
         d2(1,j,3,1)=-rdat%r02(m, 7)
         d2(2,j,3,1)=-rdat%r02(m, 8)
         d2(3,j,3,1)=-rdat%r02(m, 9)
         d2(4,j,3,1)=-rdat%r02(m,10)
      enddo

         j=1
         k= ind(j,3,1)
      do i=1,5
         d1(i  ,3,1)=-rdat%r01(k,i+ 8)
      enddo

      do j=1,10
         k= ind(j,4,1)
         d4(  j,4,1)=-rdat%r04( k,1)+rdat%r03( j,2)
      enddo

      do j=1,6
         k= ind(j,4,1)
         m= in6(j)
         d3(1,j,4,1)=-rdat%r03(k, 3)+rdat%r02(m, 4)
         d3(2,j,4,1)=-rdat%r03(k, 4)+rdat%r02(m, 5)
         d3(3,j,4,1)=-rdat%r03(k, 5)+rdat%r02(m, 6)
      enddo

      do j=1,3
         k= ind(j,4,1)
         m= in6(k)
         d2(1,j,4,1)=-rdat%r02(m, 7)+rdat%r01(j, 5)
         d2(2,j,4,1)=-rdat%r02(m, 8)+rdat%r01(j, 6)
         d2(3,j,4,1)=-rdat%r02(m, 9)+rdat%r01(j, 7)
         d2(4,j,4,1)=-rdat%r02(m,10)+rdat%r01(j, 8)
      enddo

         j=1
         k= ind(j,4,1)
      do i=1,5
         d1(i  ,4,1)=-rdat%r01(k,i+ 8)+rdat%r00(i,2)
      enddo

      xx= qx*qx
      zz= qz*qz
      xz= qx*qz
      xxx= xx*qx
      xxz= xx*qz
      xzz= zz*qx
      zzz= zz*qz

      qxd= qx+qx
      qzd= qz+qz

      xzd= xz+xz

      l=1
      do k=2,kx

         f(1,1,k-1,l) = d4(  1,k,l)+d2(2,1,k,l)* 3 +(+d3(2,1,k,l)* 2 +d3(3,1,k,l) +d1(  3,k,l)* 2 +d1(  4,k,l))*qx +(+d2(3,1,k,l)+d2(4,1,k,l)* 2 )*xx +d1(  5,k,l)*xxx
         f(2,1,k-1,l) = d4(  4,k,l)+d2(2,1,k,l) +(+d3(3,4,k,l)+d1(  4,k,l))*qx
         f(3,1,k-1,l) = d4(  6,k,l)+d2(2,1,k,l) +(+d3(3,6,k,l)+d1(  4,k,l))*qx +d3(2,3,k,l)*qzd +d2(4,3,k,l)*xzd +d2(3,1,k,l)*zz +d1(  5,k,l)*xzz
         f(4,1,k-1,l) = d4(  2,k,l)+d2(2,2,k,l) +(+d3(2,2,k,l)+d3(3,2,k,l))*qx +d2(4,2,k,l)*xx
         f(5,1,k-1,l) = d4(  3,k,l)+d2(2,3,k,l) +(+d3(2,3,k,l)+d3(3,3,k,l))*qx +(+d3(2,1,k,l)+d1(  3,k,l))*qz +d2(4,3,k,l)*xx +(+d2(3,1,k,l)+d2(4,1,k,l))*xz +d1(  5,k,l)*xxz
         f(6,1,k-1,l) = d4(  5,k,l) +d3(3,5,k,l)*qx +d3(2,2,k,l)*qz +d2(4,2,k,l)*xz

         f(1,2,k-1,l) = d4(  2,k,l)+d2(2,2,k,l) +d3(2,2,k,l)*qxd +d2(3,2,k,l)*xx
         f(2,2,k-1,l) = d4(  7,k,l)+d2(2,2,k,l)* 3
         f(3,2,k-1,l) = d4(  9,k,l)+d2(2,2,k,l) +d3(2,5,k,l)*qzd +d2(3,2,k,l)*zz
         f(4,2,k-1,l) = d4(  4,k,l)+d2(2,1,k,l) +(+d3(2,4,k,l)+d1(  3,k,l))*qx
         f(5,2,k-1,l) = d4(  5,k,l) +d3(2,5,k,l)*qx +d3(2,2,k,l)*qz +d2(3,2,k,l)*xz
         f(6,2,k-1,l) = d4(  8,k,l)+d2(2,3,k,l) +(+d3(2,4,k,l)+d1(  3,k,l))*qz

         f(1,3,k-1,l) = d4(  3,k,l)+d2(2,3,k,l) +d3(2,3,k,l)*qxd +(+d3(3,1,k,l)+d1(  4,k,l))*qz +d2(3,3,k,l)*xx +d2(4,1,k,l)*xzd +d1(  5,k,l)*xxz
         f(2,3,k-1,l) = d4(  8,k,l)+d2(2,3,k,l) +(+d3(3,4,k,l)+d1(  4,k,l))*qz
         f(3,3,k-1,l) = d4( 10,k,l)+d2(2,3,k,l)* 3 +(+d3(2,6,k,l)* 2 +d3(3,6,k,l) +d1(  3,k,l)* 2 +d1(  4,k,l))*qz +(+d2(3,3,k,l)+d2(4,3,k,l)* 2 )*zz +d1(  5,k,l)*zzz
         f(4,3,k-1,l) = d4(  5,k,l) +d3(2,5,k,l)*qx +d3(3,2,k,l)*qz +d2(4,2,k,l)*xz
         f(5,3,k-1,l) = d4(  6,k,l)+d2(2,1,k,l) +(+d3(2,6,k,l)+d1(  3,k,l))*qx +(+d3(2,3,k,l)+d3(3,3,k,l))*qz +(+d2(3,3,k,l)+d2(4,3,k,l))*xz +d2(4,1,k,l)*zz +d1(  5,k,l)*xzz
         f(6,3,k-1,l) = d4(  9,k,l)+d2(2,2,k,l) +(+d3(2,5,k,l)+d3(3,5,k,l))*qz +d2(4,2,k,l)*zz

      enddo

      end subroutine mcdv_11

! >
! >    @brief   dsds case
! >
! >    @details integration of a dsds case
! >
      subroutine mcdv_12(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,1,6,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 6, lx = 1
      real(kind=dp) :: c1(  2,kx,lx),c2(  3,kx,lx),c3(  6,kx,lx)

      integer, parameter :: ind (6, kx, lx) = reshape(&
              [ 1, 2, 3, 4, 5, 6, 4, 7, 8, 11, 12, 13, 6, 9, 10, 13, 14, 15, 2, &
                4, 5, 7, 8, 9, 3, 5, 6, 8, 9, 10, 5, 8, 9, 12, 13, 14 ] &
              , shape(ind))

      integer :: j, k, l, m

      do j=1,6
         m= in6(j)
         c3(j,1,1)=+rdat%r04( j, 1)+rdat%r02( m, 1)
      enddo

      c2(  1,1,1)=+rdat%r03( 1, 2)+rdat%r01( 1, 1)

      c2(  2,1,1)=+rdat%r03( 2, 2)+rdat%r01( 2, 1)
      c2(  3,1,1)=+rdat%r03( 3, 2)+rdat%r01( 3, 1)

      c1(  1,1,1)=+rdat%r02( 1, 4)+rdat%r00( 1, 1)
      c1(  2,1,1)=+rdat%r02( 1, 5)+rdat%r00( 2, 1)
      do j=1,6
         k= ind(j,2,1)
         m= in6(j)
         c3(j,2,1)=+rdat%r04( k, 1)+rdat%r02( m, 1)
      enddo

      c2(  1,2,1)=+rdat%r03( 4, 2)+rdat%r01( 1, 1)

      c2(  2,2,1)=+rdat%r03( 7, 2)+rdat%r01( 2, 1)
      c2(  3,2,1)=+rdat%r03( 8, 2)+rdat%r01( 3, 1)

      c1(  1,2,1)=+rdat%r02( 2, 4)+rdat%r00( 1, 1)
      c1(  2,2,1)=+rdat%r02( 2, 5)+rdat%r00( 2, 1)
      c3(  1,3,1)=+rdat%r04( 6, 1)-rdat%r03( 3, 1)* 2 +rdat%r02( 1, 1)+rdat%r02( 1, 2)

      c3(  2,3,1)=+rdat%r04( 9, 1)-rdat%r03( 5, 1)* 2 +rdat%r02( 4, 1)+rdat%r02( 4, 2)
      c3(  3,3,1)=+rdat%r04(10, 1)-rdat%r03( 6, 1)* 2 +rdat%r02( 5, 1)+rdat%r02( 5, 2)

      c3(  4,3,1)=+rdat%r04(13, 1)-rdat%r03( 8, 1)* 2 +rdat%r02( 2, 1)+rdat%r02( 2, 2)
      c3(  5,3,1)=+rdat%r04(14, 1)-rdat%r03( 9, 1)* 2 +rdat%r02( 6, 1)+rdat%r02( 6, 2)
      c3(  6,3,1)=+rdat%r04(15, 1)-rdat%r03(10, 1)* 2 +rdat%r02( 3, 1)+rdat%r02( 3, 2)

      c2(  1,3,1)=+rdat%r03( 6, 2)-rdat%r02( 5, 3)* 2 +rdat%r01( 1, 1)+rdat%r01( 1, 2)

      c2(  2,3,1)=+rdat%r03( 9, 2)-rdat%r02( 6, 3)* 2 +rdat%r01( 2, 1)+rdat%r01( 2, 2)
      c2(  3,3,1)=+rdat%r03(10, 2)-rdat%r02( 3, 3)* 2 +rdat%r01( 3, 1)+rdat%r01( 3, 2)

      c1(  1,3,1)=+rdat%r02( 3, 4)-rdat%r01( 3, 3)* 2 +rdat%r00( 1, 1)+rdat%r00( 3, 1)
      c1(  2,3,1)=+rdat%r02( 3, 5)-rdat%r01( 3, 4)* 2 +rdat%r00( 2, 1)+rdat%r00( 4, 1)
      do j=1,6
         k= ind(j,4,1)
         c3(j,4,1)=+rdat%r04( k, 1)
      enddo

      c2(  1,4,1)=+rdat%r03( 2, 2)

      c2(  2,4,1)=+rdat%r03( 4, 2)
      c2(  3,4,1)=+rdat%r03( 5, 2)

      c1(  1,4,1)=+rdat%r02( 4, 4)
      c1(  2,4,1)=+rdat%r02( 4, 5)
      do j=1,6
         k= ind(j,5,1)
         c3(j,5,1)=+rdat%r04( k, 1)-rdat%r03( j, 1)
      enddo

      c2(  1,5,1)=+rdat%r03( 3, 2)-rdat%r02( 1, 3)

      c2(  2,5,1)=+rdat%r03( 5, 2)-rdat%r02( 4, 3)
      c2(  3,5,1)=+rdat%r03( 6, 2)-rdat%r02( 5, 3)

      c1(  1,5,1)=+rdat%r02( 5, 4)-rdat%r01( 1, 3)
      c1(  2,5,1)=+rdat%r02( 5, 5)-rdat%r01( 1, 4)
      c3(  1,6,1)=+rdat%r04( 5, 1)-rdat%r03( 2, 1)

      c3(  2,6,1)=+rdat%r04( 8, 1)-rdat%r03( 4, 1)
      c3(  3,6,1)=+rdat%r04( 9, 1)-rdat%r03( 5, 1)

      c3(  4,6,1)=+rdat%r04(12, 1)-rdat%r03( 7, 1)
      c3(  5,6,1)=+rdat%r04(13, 1)-rdat%r03( 8, 1)
      c3(  6,6,1)=+rdat%r04(14, 1)-rdat%r03( 9, 1)

      c2(  1,6,1)=+rdat%r03( 5, 2)-rdat%r02( 4, 3)

      c2(  2,6,1)=+rdat%r03( 8, 2)-rdat%r02( 2, 3)
      c2(  3,6,1)=+rdat%r03( 9, 2)-rdat%r02( 6, 3)

      c1(  1,6,1)=+rdat%r02( 6, 4)-rdat%r01( 2, 3)
      c1(  2,6,1)=+rdat%r02( 6, 5)-rdat%r01( 2, 4)

      do l=1,lx
         do k=1,kx

            f(1,1,k,l)=+c3(  1,k,l)+c1(  1,k,l)                         &
     &               +(+c2(  1,k,l)+c2(  1,k,l)+c1(  2,k,l)*qx)*qx
            f(2,1,k,l)=+c3(  4,k,l)+c1(  1,k,l)
            f(3,1,k,l)=+c3(  6,k,l)+c1(  1,k,l)                         &
     &               +(+c2(  3,k,l)+c2(  3,k,l)+c1(  2,k,l)*qz)*qz
            f(4,1,k,l)=+c3(  2,k,l)+c2(  2,k,l)*qx
            f(5,1,k,l)=+c3(  3,k,l)+c2(  3,k,l)*qx                      &
     &               +(+c2(  1,k,l)+c1(  2,k,l)*qx)*qz
            f(6,1,k,l)=+c3(  5,k,l)+c2(  2,k,l)*qz

         enddo
      enddo

      end subroutine mcdv_12

! >
! >    @brief   dspp case
! >
! >    @details integration of a dspp case
! >
      subroutine mcdv_13(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,1,3,*)
      real(kind=dp) :: qx, qz

      logical :: lsym13

      integer, parameter :: kx = 4, lx = 4
      real(kind=dp) :: c1(  2,kx,lx),c2(  3,kx,lx),c3(  6,kx,lx)

      integer, parameter :: ind (6, kx, lx) = reshape(&
        [ 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 2, 4, 5, 7, 8, 9, 3, 5, 6, 8, &
          9, 10, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 2, 4, 5, 7, 8, 9, 3, 5, &
          6, 8, 9, 10, 2, 4, 5, 7, 8, 9, 2, 4, 5, 7, 8, 9, 4, 7, 8, 11, 12, &
          13, 5, 8, 9, 12, 13, 14, 3, 5, 6, 8, 9, 10, 3, 5, 6, 8, 9, 10, 5, &
          8, 9, 12, 13, 14, 6, 9, 10, 13, 14, 15 ] &
        , shape(ind))
      integer :: i, j, k, l, m
      do j=1,6
         m= in6(j)
         c3(j,1,1)=+rdat%r02( m, 1)
      enddo

      c2(  1,1,1)=+rdat%r01( 1, 1)

      c2(  2,1,1)=+rdat%r01( 2, 1)
      c2(  3,1,1)=+rdat%r01( 3, 1)

      c1(  1,1,1)=+rdat%r00(1,1)
      c1(  2,1,1)=+rdat%r00(1,2)
      do j=1,6
         c3(j,2,1)=-rdat%r03( j, 2)
      enddo

      c2(  1,2,1)=-rdat%r02( 1, 6)

      c2(  2,2,1)=-rdat%r02( 4, 6)
      c2(  3,2,1)=-rdat%r02( 5, 6)

      c1(  1,2,1)=-rdat%r01( 1, 6)
      c1(  2,2,1)=-rdat%r01( 1, 7)
      do j=1,6
         k= ind(j,3,1)
         c3(j,3,1)=-rdat%r03( k, 2)
      enddo

      c2(  1,3,1)=-rdat%r02( 4, 6)

      c2(  2,3,1)=-rdat%r02( 2, 6)
      c2(  3,3,1)=-rdat%r02( 6, 6)

      c1(  1,3,1)=-rdat%r01( 2, 6)
      c1(  2,3,1)=-rdat%r01( 2, 7)
      do j=1,6
         k= ind(j,4,1)
         m= in6(j)
         c3(j,4,1)=-rdat%r03( k, 2)+rdat%r02( m, 2)
      enddo

      c2(  1,4,1)=-rdat%r02( 5, 6)+rdat%r01( 1, 2)

      c2(  2,4,1)=-rdat%r02( 6, 6)+rdat%r01( 2, 2)
      c2(  3,4,1)=-rdat%r02( 3, 6)+rdat%r01( 3, 2)

      c1(  1,4,1)=-rdat%r01( 3, 6)+rdat%r00( 2, 1)
      c1(  2,4,1)=-rdat%r01( 3, 7)+rdat%r00( 2, 2)
      do j=1,6
         c3(j,1,2)=-rdat%r03( j, 3)
      enddo

      c2(  1,1,2)=-rdat%r02( 1, 7)

      c2(  2,1,2)=-rdat%r02( 4, 7)
      c2(  3,1,2)=-rdat%r02( 5, 7)

      c1(  1,1,2)=-rdat%r01( 1, 8)
      c1(  2,1,2)=-rdat%r01( 1, 9)
      do j=1,6
         m= in6(j)
         c3(j,2,2)=+rdat%r04( j, 1)+rdat%r02( m, 5)
      enddo

      c2(  1,2,2)=+rdat%r03( 1, 5)+rdat%r01( 1, 5)

      c2(  2,2,2)=+rdat%r03( 2, 5)+rdat%r01( 2, 5)
      c2(  3,2,2)=+rdat%r03( 3, 5)+rdat%r01( 3, 5)

      c1(  1,2,2)=+rdat%r02( 1,10)+rdat%r00(5,1)
      c1(  2,2,2)=+rdat%r02( 1,11)+rdat%r00(5,2)
      do j=1,6
         k= ind(j,3,2)
         c3(j,3,2)=+rdat%r04( k, 1)
      enddo

      c2(  1,3,2)=+rdat%r03( 2, 5)

      c2(  2,3,2)=+rdat%r03( 4, 5)
      c2(  3,3,2)=+rdat%r03( 5, 5)

      c1(  1,3,2)=+rdat%r02( 4,10)
      c1(  2,3,2)=+rdat%r02( 4,11)
      do j=1,6
         k= ind(j,4,2)
         m= in6(j)
         c3(j,4,2)=+rdat%r04( k, 1)-rdat%r03( j, 4)
      enddo

      c2(  1,4,2)=+rdat%r03( 3, 5)-rdat%r02( 1, 8)

      c2(  2,4,2)=+rdat%r03( 5, 5)-rdat%r02( 4, 8)
      c2(  3,4,2)=+rdat%r03( 6, 5)-rdat%r02( 5, 8)

      c1(  1,4,2)=+rdat%r02( 5,10)-rdat%r01( 1,10)
      c1(  2,4,2)=+rdat%r02( 5,11)-rdat%r01( 1,11)
      do j=1,6
         k= ind(j,1,3)
         c3(j,1,3)=-rdat%r03( k, 3)
      enddo

      c2(  1,1,3)=-rdat%r02( 4, 7)

      c2(  2,1,3)=-rdat%r02( 2, 7)
      c2(  3,1,3)=-rdat%r02( 6, 7)

      c1(  1,1,3)=-rdat%r01( 2, 8)
      c1(  2,1,3)=-rdat%r01( 2, 9)
      do j=1,6
         k= ind(j,3,3)
         m= in6(j)
         c3(j,3,3)=+rdat%r04( k, 1)+rdat%r02( m, 5)
      enddo

      c2(  1,3,3)=+rdat%r03( 4, 5)+rdat%r01( 1, 5)

      c2(  2,3,3)=+rdat%r03( 7, 5)+rdat%r01( 2, 5)
      c2(  3,3,3)=+rdat%r03( 8, 5)+rdat%r01( 3, 5)

      c1(  1,3,3)=+rdat%r02( 2,10)+rdat%r00(5,1)
      c1(  2,3,3)=+rdat%r02( 2,11)+rdat%r00(5,2)
      c3(  1,4,3)=+rdat%r04( 5, 1)-rdat%r03( 2, 4)

      c3(  2,4,3)=+rdat%r04( 8, 1)-rdat%r03( 4, 4)
      c3(  3,4,3)=+rdat%r04( 9, 1)-rdat%r03( 5, 4)

      c3(  4,4,3)=+rdat%r04(12, 1)-rdat%r03( 7, 4)
      c3(  5,4,3)=+rdat%r04(13, 1)-rdat%r03( 8, 4)
      c3(  6,4,3)=+rdat%r04(14, 1)-rdat%r03( 9, 4)

      c2(  1,4,3)=+rdat%r03( 5, 5)-rdat%r02( 4, 8)

      c2(  2,4,3)=+rdat%r03( 8, 5)-rdat%r02( 2, 8)
      c2(  3,4,3)=+rdat%r03( 9, 5)-rdat%r02( 6, 8)

      c1(  1,4,3)=+rdat%r02( 6,10)-rdat%r01( 2,10)
      c1(  2,4,3)=+rdat%r02( 6,11)-rdat%r01( 2,11)
      do j=1,6
         k= ind(j,1,4)
         m= in6(j)
         c3(j,1,4)=-rdat%r03( k, 3)+rdat%r02( m, 3)
      enddo

      c2(  1,1,4)=-rdat%r02( 5, 7)+rdat%r01( 1, 3)

      c2(  2,1,4)=-rdat%r02( 6, 7)+rdat%r01( 2, 3)
      c2(  3,1,4)=-rdat%r02( 3, 7)+rdat%r01( 3, 3)

      c1(  1,1,4)=-rdat%r01( 3, 8)+rdat%r00(3,1)
      c1(  2,1,4)=-rdat%r01( 3, 9)+rdat%r00(3,2)
      do j=1,6
         k= ind(j,2,4)
         c3(j,2,4)=+rdat%r04( k, 1)-rdat%r03( j, 1)
      enddo

      c2(  1,2,4)=+rdat%r03( 3, 5)-rdat%r02( 1, 9)

      c2(  2,2,4)=+rdat%r03( 5, 5)-rdat%r02( 4, 9)
      c2(  3,2,4)=+rdat%r03( 6, 5)-rdat%r02( 5, 9)

      c1(  1,2,4)=+rdat%r02( 5,10)-rdat%r01( 1,12)
      c1(  2,2,4)=+rdat%r02( 5,11)-rdat%r01( 1,13)
      c3(  1,3,4)=+rdat%r04( 5, 1)-rdat%r03( 2, 1)

      c3(  2,3,4)=+rdat%r04( 8, 1)-rdat%r03( 4, 1)
      c3(  3,3,4)=+rdat%r04( 9, 1)-rdat%r03( 5, 1)

      c3(  4,3,4)=+rdat%r04(12, 1)-rdat%r03( 7, 1)
      c3(  5,3,4)=+rdat%r04(13, 1)-rdat%r03( 8, 1)
      c3(  6,3,4)=+rdat%r04(14, 1)-rdat%r03( 9, 1)

      c2(  1,3,4)=+rdat%r03( 5, 5)-rdat%r02( 4, 9)

      c2(  2,3,4)=+rdat%r03( 8, 5)-rdat%r02( 2, 9)
      c2(  3,3,4)=+rdat%r03( 9, 5)-rdat%r02( 6, 9)

      c1(  1,3,4)=+rdat%r02( 6,10)-rdat%r01( 2,12)
      c1(  2,3,4)=+rdat%r02( 6,11)-rdat%r01( 2,13)
      c3(  1,4,4)=+rdat%r04( 6, 1)-rdat%r03( 3, 1)-rdat%r03( 3, 4)+rdat%r02( 1, 4)+rdat%r02( 1, 5)

      c3(  2,4,4)=+rdat%r04( 9, 1)-rdat%r03( 5, 1)-rdat%r03( 5, 4)+rdat%r02( 4, 4)+rdat%r02( 4, 5)
      c3(  3,4,4)=+rdat%r04(10, 1)-rdat%r03( 6, 1)-rdat%r03( 6, 4)+rdat%r02( 5, 4)+rdat%r02( 5, 5)

      c3(  4,4,4)=+rdat%r04(13, 1)-rdat%r03( 8, 1)-rdat%r03( 8, 4)+rdat%r02( 2, 4)+rdat%r02( 2, 5)
      c3(  5,4,4)=+rdat%r04(14, 1)-rdat%r03( 9, 1)-rdat%r03( 9, 4)+rdat%r02( 6, 4)+rdat%r02( 6, 5)
      c3(  6,4,4)=+rdat%r04(15, 1)-rdat%r03(10, 1)-rdat%r03(10, 4)+rdat%r02( 3, 4)+rdat%r02( 3, 5)

      c2(  1,4,4)=+rdat%r03( 6, 5)-rdat%r02( 5, 8)-rdat%r02( 5, 9)+rdat%r01( 1, 4)+rdat%r01( 1, 5)

      c2(  2,4,4)=+rdat%r03( 9, 5)-rdat%r02( 6, 8)-rdat%r02( 6, 9)+rdat%r01( 2, 4)+rdat%r01( 2, 5)
      c2(  3,4,4)=+rdat%r03(10, 5)-rdat%r02( 3, 8)-rdat%r02( 3, 9)+rdat%r01( 3, 4)+rdat%r01( 3, 5)

      c1(  1,4,4)=+rdat%r02( 3,10)-rdat%r01( 3,10)-rdat%r01( 3,12)+rdat%r00(4,1)+rdat%r00(5,1)
      c1(  2,4,4)=+rdat%r02( 3,11)-rdat%r01( 3,11)-rdat%r01( 3,13)+rdat%r00(4,2)+rdat%r00(5,2)

      do l=2,lx
         do k=2,kx

           if(k == 2 .and. l == 3) cycle

            f(1,1,k-1,l-1)=+c3(  1,k,l)+c1(  1,k,l) +(+c2(  1,k,l)+c2(  1,k,l)+c1(  2,k,l)*qx)*qx
            f(2,1,k-1,l-1)=+c3(  4,k,l)+c1(  1,k,l)
            f(3,1,k-1,l-1)=+c3(  6,k,l)+c1(  1,k,l) +(+c2(  3,k,l)+c2(  3,k,l)+c1(  2,k,l)*qz)*qz
            f(4,1,k-1,l-1)=+c3(  2,k,l)+c2(  2,k,l)*qx
            f(5,1,k-1,l-1)=+c3(  3,k,l)+c2(  3,k,l)*qx +(+c2(  1,k,l)+c1(  2,k,l)*qx)*qz
            f(6,1,k-1,l-1)=+c3(  5,k,l)+c2(  2,k,l)*qz

         enddo
      enddo

     do i=1,6
        f(i,1,1,2)= f(i,1,2,1)
     enddo

      end subroutine mcdv_13

! >
! >    @brief   ddps case
! >
! >    @details integration of a ddps case
! >
      subroutine mcdv_14(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,6,3,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 4, lx = 1
      real(kind=dp) :: e1(  5,kx,lx),e2(4,3,kx,lx),e3(4,6,kx,lx), &
                       e4(2,10,kx,lx),e5( 15,kx,lx)

      integer, parameter :: ind (15, kx, lx) = reshape(&
        [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1, 2, 3, 4, 5, &
          6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 2, 4, 5, 7, 8, 9, 11, 12, 13, &
          14, 16, 17, 18, 19, 20, 3, 5, 6, 8, 9, 10, 12, 13, 14, 15, 17, 18, &
          19, 20, 21] &
        , shape(ind))
      integer :: i, j, k, l, m
      real(kind=dp) :: xx, zz, xz, xxx, xxz, xzz, zzz
      real(kind=dp) :: xxxx, xxzz, xxxz, xzzz, zzzz
      real(kind=dp) :: xxxd, xxzd, xzzd, zzzd
      real(kind=dp) :: xzq
      real(kind=dp) :: qxd, qzd, xzd

      do j=1,15
         e5(  j,1,1)=+rdat%r04(j,1)
      enddo

      do j=1,10
         e4(1,j,1,1)=+rdat%r03(j,1)
         e4(2,j,1,1)=+rdat%r03(j,2)
      enddo

      do j=1,6
         m= in6(j)
         e3(1,j,1,1)=+rdat%r02(m,1)
         e3(2,j,1,1)=+rdat%r02(m,2)
         e3(3,j,1,1)=+rdat%r02(m,3)
         e3(4,j,1,1)=+rdat%r02(m,4)
      enddo

      do j=1,3
         e2(1,j,1,1)=+rdat%r01(j,1)
         e2(2,j,1,1)=+rdat%r01(j,2)
         e2(3,j,1,1)=+rdat%r01(j,3)
         e2(4,j,1,1)=+rdat%r01(j,4)
      enddo

      do i=1,5
         e1(i  ,1,1)=+rdat%r00(i,1)
      enddo
      do j=1,15
         e5(  j,2,1)=-rdat%r05( j,1)
      enddo

      do j=1,10
         e4(1,j,2,1)=-rdat%r04(j, 3)
         e4(2,j,2,1)=-rdat%r04(j, 4)
      enddo

      do j=1,6
         e3(1,j,2,1)=-rdat%r03(j, 5)
         e3(2,j,2,1)=-rdat%r03(j, 6)
         e3(3,j,2,1)=-rdat%r03(j, 7)
         e3(4,j,2,1)=-rdat%r03(j, 8)
      enddo

      do j=1,3
         m= in6(j)
         e2(1,j,2,1)=-rdat%r02(m, 9)
         e2(2,j,2,1)=-rdat%r02(m,10)
         e2(3,j,2,1)=-rdat%r02(m,11)
         e2(4,j,2,1)=-rdat%r02(m,12)
      enddo

         j=1
      do i=1,5
         e1(i  ,2,1)=-rdat%r01(j,i+ 8)
      enddo
      do j=1,15
         k= ind(j,3,1)
         e5(  j,3,1)=-rdat%r05( k,1)
      enddo

      do j=1,10
         k= ind(j,3,1)
         e4(1,j,3,1)=-rdat%r04(k, 3)
         e4(2,j,3,1)=-rdat%r04(k, 4)
      enddo

      do j=1,6
         k= ind(j,3,1)
         e3(1,j,3,1)=-rdat%r03(k, 5)
         e3(2,j,3,1)=-rdat%r03(k, 6)
         e3(3,j,3,1)=-rdat%r03(k, 7)
         e3(4,j,3,1)=-rdat%r03(k, 8)
      enddo

      do j=1,3
         k= ind(j,3,1)
         m= in6(k)
         e2(1,j,3,1)=-rdat%r02(m, 9)
         e2(2,j,3,1)=-rdat%r02(m,10)
         e2(3,j,3,1)=-rdat%r02(m,11)
         e2(4,j,3,1)=-rdat%r02(m,12)
      enddo

         j=1
         k= ind(j,3,1)
      do i=1,5
         e1(i  ,3,1)=-rdat%r01(k,i+ 8)
      enddo
      do j=1,15
         k= ind(j,4,1)
         e5(  j,4,1)=-rdat%r05( k,1)+rdat%r04(j,2)
      enddo

      do j=1,10
         k= ind(j,4,1)
         e4(1,j,4,1)=-rdat%r04(k, 3)+rdat%r03(j,3)
         e4(2,j,4,1)=-rdat%r04(k, 4)+rdat%r03(j,4)
      enddo

      do j=1,6
         k= ind(j,4,1)
         m= in6(j)
         e3(1,j,4,1)=-rdat%r03(k, 5)+rdat%r02(m, 5)
         e3(2,j,4,1)=-rdat%r03(k, 6)+rdat%r02(m, 6)
         e3(3,j,4,1)=-rdat%r03(k, 7)+rdat%r02(m, 7)
         e3(4,j,4,1)=-rdat%r03(k, 8)+rdat%r02(m, 8)
      enddo

      do j=1,3
         k= ind(j,4,1)
         m= in6(k)
         e2(1,j,4,1)=-rdat%r02(m, 9)+rdat%r01(j, 5)
         e2(2,j,4,1)=-rdat%r02(m,10)+rdat%r01(j, 6)
         e2(3,j,4,1)=-rdat%r02(m,11)+rdat%r01(j, 7)
         e2(4,j,4,1)=-rdat%r02(m,12)+rdat%r01(j, 8)
      enddo

         j=1
         k= ind(j,4,1)
      do i=1,5
         e1(i  ,4,1)=-rdat%r01(k,i+ 8)+rdat%r00(i,2)
      enddo

      xx= qx*qx
      zz= qz*qz
      xz= qx*qz
      xxx= xx*qx
      xxz= xx*qz
      xzz= zz*qx
      zzz= zz*qz
      xxxx= xx*xx
      xxxz= xx*xz
      xxzz= xx*zz
      xzzz= xz*zz
      zzzz= zz*zz

      qxd= qx+qx
      qzd= qz+qz
      xzd= xz+xz
      xxxd= xxx+xxx
      xxzd= xxz+xxz
      xzzd= xzz+xzz
      zzzd= zzz+zzz

      xzq= xz+xz+xz+xz

      do k=2,kx

         f(1,1,k-1,1) = e5(  1,k,1)+e3(1,1,k,1)* 6 +e1(  1,k,1)* 3 +(+e4(1,1,k,1)+e4(2,1,k,1)+ (+e2(1,1,k,1)+e2(2,1,k,1))* 3 )*qxd +(+e3(2,1,k,1)+e3(3,1,k,1)* 4 +e3(4,1,k,1) +e1(  2,k,1)+e1(  3,k,1)* 4 +e1(  4,k,1))*xx +(+e2(3,1,k,1)+e2(4,1,k,1))*xxxd +e1(  5,k,1)*xxxx
         f(2,1,k-1,1) = e5(  4,k,1)+e3(1,4,k,1)+e3(1,1,k,1)+e1(  1,k,1) +(+e4(2,4,k,1)+e2(2,1,k,1))*qxd +(+e3(4,4,k,1)+e1(  4,k,1))*xx
         f(3,1,k-1,1) = e5(  6,k,1)+e3(1,6,k,1)+e3(1,1,k,1)+e1(  1,k,1) +(+e4(2,6,k,1)+e2(2,1,k,1))*qxd +(+e4(1,3,k,1)+e2(1,3,k,1))*qzd +(+e3(4,6,k,1)+e1(  4,k,1))*xx +e3(3,3,k,1)*xzq +(+e3(2,1,k,1)+e1(  2,k,1))*zz +e2(4,3,k,1)*xxzd +e2(3,1,k,1)*xzzd +e1(  5,k,1)*xxzz
         f(4,1,k-1,1) = e5(  2,k,1)+e3(1,2,k,1)* 3 +(+e4(1,2,k,1)+e4(2,2,k,1)* 2 +e2(1,2,k,1)+e2(2,2,k,1)* 2 )*qx +(+e3(3,2,k,1)* 2 +e3(4,2,k,1))*xx +e2(4,2,k,1)*xxx
         f(5,1,k-1,1) = e5(  3,k,1)+e3(1,3,k,1)* 3 +(+e4(1,3,k,1)+e4(2,3,k,1)* 2 +e2(1,3,k,1)+e2(2,3,k,1)* 2 )*qx +(+e4(1,1,k,1)+e2(1,1,k,1)* 3 )*qz +(+e3(3,3,k,1)* 2 +e3(4,3,k,1))*xx +(+e3(2,1,k,1)+e3(3,1,k,1)* 2 +e1(  2,k,1)+e1(  3,k,1)* 2 )*xz +e2(4,3,k,1)*xxx +(+e2(3,1,k,1)* 2 +e2(4,1,k,1))*xxz +e1(  5,k,1)*xxxz
         f(6,1,k-1,1) = e5(  5,k,1)+e3(1,5,k,1) +e4(2,5,k,1)*qxd +(+e4(1,2,k,1)+e2(1,2,k,1))*qz +e3(4,5,k,1)*xx +e3(3,2,k,1)*xzd +e2(4,2,k,1)*xxz

         f(1,2,k-1,1) = e5(  4,k,1)+e3(1,4,k,1)+e3(1,1,k,1)+e1(  1,k,1) +(+e4(1,4,k,1)+e2(1,1,k,1))*qxd +(+e3(2,4,k,1)+e1(  2,k,1))*xx
         f(2,2,k-1,1) = e5( 11,k,1)+e3(1,4,k,1)* 6 +e1(  1,k,1)* 3
         f(3,2,k-1,1) = e5( 13,k,1)+e3(1,6,k,1)+e3(1,4,k,1)+e1(  1,k,1) +(+e4(1,8,k,1)+e2(1,3,k,1))*qzd +(+e3(2,4,k,1)+e1(  2,k,1))*zz
         f(4,2,k-1,1) = e5(  7,k,1)+e3(1,2,k,1)* 3 +(+e4(1,7,k,1)+e2(1,2,k,1)* 3 )*qx
         f(5,2,k-1,1) = e5(  8,k,1)+e3(1,3,k,1) +(+e4(1,8,k,1)+e2(1,3,k,1))*qx +(+e4(1,4,k,1)+e2(1,1,k,1))*qz +(+e3(2,4,k,1)+e1(  2,k,1))*xz
         f(6,2,k-1,1) = e5( 12,k,1)+e3(1,5,k,1)* 3 +(+e4(1,7,k,1)+e2(1,2,k,1)* 3 )*qz

         f(1,3,k-1,1) = e5(  6,k,1)+e3(1,6,k,1)+e3(1,1,k,1)+e1(  1,k,1) +(+e4(1,6,k,1)+e2(1,1,k,1))*qxd +(+e4(2,3,k,1)+e2(2,3,k,1))*qzd +(+e3(2,6,k,1)+e1(  2,k,1))*xx +e3(3,3,k,1)*xzq +(+e3(4,1,k,1)+e1(  4,k,1))*zz +e2(3,3,k,1)*xxzd +e2(4,1,k,1)*xzzd +e1(  5,k,1)*xxzz
         f(2,3,k-1,1) = e5( 13,k,1)+e3(1,6,k,1)+e3(1,4,k,1)+e1(  1,k,1) +(+e4(2,8,k,1)+e2(2,3,k,1))*qzd +(+e3(4,4,k,1)+e1(  4,k,1))*zz
         f(3,3,k-1,1) = e5( 15,k,1)+e3(1,6,k,1)* 6 +e1(  1,k,1)* 3 +(+e4(1,10,k,1)+e4(2,10,k,1)+ (+e2(1,3,k,1)+e2(2,3,k,1))* 3 )*qzd +(+e3(2,6,k,1)+e3(3,6,k,1)* 4 +e3(4,6,k,1) +e1(  2,k,1)+e1(  3,k,1)* 4 +e1(  4,k,1))*zz +(+e2(3,3,k,1)+e2(4,3,k,1))*zzzd +e1(  5,k,1)*zzzz
         f(4,3,k-1,1) = e5(  9,k,1)+e3(1,2,k,1) +(+e4(1,9,k,1)+e2(1,2,k,1))*qx +e4(2,5,k,1)*qzd +e3(3,5,k,1)*xzd +e3(4,2,k,1)*zz +e2(4,2,k,1)*xzz
         f(5,3,k-1,1) = e5( 10,k,1)+e3(1,3,k,1)* 3 +(+e4(1,10,k,1)+e2(1,3,k,1)* 3 )*qx +(+e4(1,6,k,1)+e4(2,6,k,1)* 2 +e2(1,1,k,1)+e2(2,1,k,1)* 2 )*qz +(+e3(2,6,k,1)+e3(3,6,k,1)* 2 +e1(  2,k,1)+e1(  3,k,1)* 2 )*xz +(+e3(3,3,k,1)* 2 +e3(4,3,k,1))*zz +(+e2(3,3,k,1)* 2 +e2(4,3,k,1))*xzz +e2(4,1,k,1)*zzz +e1(  5,k,1)*xzzz
         f(6,3,k-1,1) = e5( 14,k,1)+e3(1,5,k,1)* 3 +(+e4(1,9,k,1)+e4(2,9,k,1)* 2 +e2(1,2,k,1)+e2(2,2,k,1)* 2 )*qz +(+e3(3,5,k,1)* 2 +e3(4,5,k,1))*zz +e2(4,2,k,1)*zzz

         f(1,4,k-1,1) = e5(  2,k,1)+e3(1,2,k,1)* 3 +(+e4(1,2,k,1)* 2 +e4(2,2,k,1) +e2(1,2,k,1)* 2 +e2(2,2,k,1))*qx +(+e3(2,2,k,1)+e3(3,2,k,1)* 2 )*xx +e2(3,2,k,1)*xxx
         f(2,4,k-1,1) = e5(  7,k,1)+e3(1,2,k,1)* 3 +(+e4(2,7,k,1)+e2(2,2,k,1)* 3 )*qx
         f(3,4,k-1,1) = e5(  9,k,1)+e3(1,2,k,1) +(+e4(2,9,k,1)+e2(2,2,k,1))*qx +e4(1,5,k,1)*qzd +e3(3,5,k,1)*xzd +e3(2,2,k,1)*zz +e2(3,2,k,1)*xzz
         f(4,4,k-1,1) = e5(  4,k,1)+e3(1,4,k,1)+e3(1,1,k,1)+e1(  1,k,1) +(+e4(1,4,k,1)+e4(2,4,k,1) +e2(1,1,k,1)+e2(2,1,k,1))*qx +(+e3(3,4,k,1)+e1(  3,k,1))*xx
         f(5,4,k-1,1) = e5(  5,k,1)+e3(1,5,k,1) +(+e4(1,5,k,1)+e4(2,5,k,1))*qx +(+e4(1,2,k,1)+e2(1,2,k,1))*qz +e3(3,5,k,1)*xx +(+e3(2,2,k,1)+e3(3,2,k,1))*xz +e2(3,2,k,1)*xxz
         f(6,4,k-1,1) = e5(  8,k,1)+e3(1,3,k,1) +(+e4(2,8,k,1)+e2(2,3,k,1))*qx +(+e4(1,4,k,1)+e2(1,1,k,1))*qz +(+e3(3,4,k,1)+e1(  3,k,1))*xz

         f(1,5,k-1,1) = e5(  3,k,1)+e3(1,3,k,1)* 3 +(+e4(1,3,k,1)* 2 +e4(2,3,k,1) +e2(1,3,k,1)* 2 +e2(2,3,k,1))*qx +(+e4(2,1,k,1)+e2(2,1,k,1)* 3 )*qz +(+e3(2,3,k,1)+e3(3,3,k,1)* 2 )*xx +(+e3(3,1,k,1)* 2 +e3(4,1,k,1) +e1(  3,k,1)* 2 +e1(  4,k,1))*xz +e2(3,3,k,1)*xxx +(+e2(3,1,k,1)+e2(4,1,k,1)* 2 )*xxz +e1(  5,k,1)*xxxz
         f(2,5,k-1,1) = e5(  8,k,1)+e3(1,3,k,1) +(+e4(2,8,k,1)+e2(2,3,k,1))*qx +(+e4(2,4,k,1)+e2(2,1,k,1))*qz +(+e3(4,4,k,1)+e1(  4,k,1))*xz
         f(3,5,k-1,1) = e5( 10,k,1)+e3(1,3,k,1)* 3 +(+e4(2,10,k,1)+e2(2,3,k,1)* 3 )*qx +(+e4(1,6,k,1)* 2 +e4(2,6,k,1) +e2(1,1,k,1)* 2 +e2(2,1,k,1))*qz +(+e3(3,6,k,1)* 2 +e3(4,6,k,1) +e1(  3,k,1)* 2 +e1(  4,k,1))*xz +(+e3(2,3,k,1)+e3(3,3,k,1)* 2 )*zz +(+e2(3,3,k,1)+e2(4,3,k,1)* 2 )*xzz +e2(3,1,k,1)*zzz +e1(  5,k,1)*xzzz
         f(4,5,k-1,1) = e5(  5,k,1)+e3(1,5,k,1) +(+e4(1,5,k,1)+e4(2,5,k,1))*qx +(+e4(2,2,k,1)+e2(2,2,k,1))*qz +e3(3,5,k,1)*xx +(+e3(3,2,k,1)+e3(4,2,k,1))*xz +e2(4,2,k,1)*xxz
         f(5,5,k-1,1) = e5(  6,k,1)+e3(1,6,k,1)+e3(1,1,k,1)+e1(  1,k,1) +(+e4(1,6,k,1)+e4(2,6,k,1) +e2(1,1,k,1)+e2(2,1,k,1))*qx +(+e4(1,3,k,1)+e4(2,3,k,1) +e2(1,3,k,1)+e2(2,3,k,1))*qz +(+e3(3,6,k,1)+e1(  3,k,1))*xx +(+e3(2,3,k,1)+e3(3,3,k,1)* 2 +e3(4,3,k,1))*xz +(+e3(3,1,k,1)+e1(  3,k,1))*zz +(+e2(3,3,k,1)+e2(4,3,k,1))*xxz +(+e2(3,1,k,1)+e2(4,1,k,1))*xzz +e1(  5,k,1)*xxzz
         f(6,5,k-1,1) = e5(  9,k,1)+e3(1,2,k,1) +(+e4(2,9,k,1)+e2(2,2,k,1))*qx +(+e4(1,5,k,1)+e4(2,5,k,1))*qz +(+e3(3,5,k,1)+e3(4,5,k,1))*xz +e3(3,2,k,1)*zz +e2(4,2,k,1)*xzz

         f(1,6,k-1,1) = e5(  5,k,1)+e3(1,5,k,1) +e4(1,5,k,1)*qxd +(+e4(2,2,k,1)+e2(2,2,k,1))*qz +e3(2,5,k,1)*xx +e3(3,2,k,1)*xzd +e2(3,2,k,1)*xxz
         f(2,6,k-1,1) = e5( 12,k,1)+e3(1,5,k,1)* 3 +(+e4(2,7,k,1)+e2(2,2,k,1)* 3 )*qz
         f(3,6,k-1,1) = e5( 14,k,1)+e3(1,5,k,1)* 3 +(+e4(1,9,k,1)* 2 +e4(2,9,k,1) +e2(1,2,k,1)* 2 +e2(2,2,k,1))*qz +(+e3(2,5,k,1)+e3(3,5,k,1)* 2 )*zz +e2(3,2,k,1)*zzz
         f(4,6,k-1,1) = e5(  8,k,1)+e3(1,3,k,1) +(+e4(1,8,k,1)+e2(1,3,k,1))*qx +(+e4(2,4,k,1)+e2(2,1,k,1))*qz +(+e3(3,4,k,1)+e1(  3,k,1))*xz
         f(5,6,k-1,1) = e5(  9,k,1)+e3(1,2,k,1) +(+e4(1,9,k,1)+e2(1,2,k,1))*qx +(+e4(1,5,k,1)+e4(2,5,k,1))*qz +(+e3(2,5,k,1)+e3(3,5,k,1))*xz +e3(3,2,k,1)*zz +e2(3,2,k,1)*xzz
         f(6,6,k-1,1) = e5( 13,k,1)+e3(1,6,k,1)+e3(1,4,k,1)+e1(  1,k,1) +(+e4(1,8,k,1)+e4(2,8,k,1) +e2(1,3,k,1)+e2(2,3,k,1))*qz +(+e3(3,4,k,1)+e1(  3,k,1))*zz

      enddo

      end subroutine mcdv_14

! >
! >    @brief   dpds case
! >
! >    @details integration of a dpds case
! >
      subroutine mcdv_15(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,3,6,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 6, lx = 1
      real(kind=dp) :: d1(  5,kx,lx),d2(4,3,kx,lx),d3(3,6,kx,lx), &
                       d4( 10,kx,lx)

      integer, parameter :: ind (10, kx, lx) = reshape(&
        [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 4, 7, 8, 11, 12, 13, 16, 17, 18, &
          19, 6, 9, 10, 13, 14, 15, 18, 19, 20, 21, 2, 4, 5, 7, 8, 9, 11, &
          12, 13, 14, 3, 5, 6, 8, 9, 10, 12, 13, 14, 15, 5, 8, 9, 12, 13, &
          14, 17, 18, 19, 20] &
        , shape(ind))
      integer :: i, j, k, l, m
      logical :: lsym16, lsym19

      real(kind=dp) :: xx, zz, xz, xxx, xxz, xzz, zzz
      real(kind=dp) :: qxd, qzd, xzd
      do j=1,10
         d4(  j,1,1)=+rdat%r05( j,1)+rdat%r03( j,1)
      enddo

      do j=1,6
         m= in6(j)
         d3(1,j,1,1)=+rdat%r04(j, 2)+rdat%r02(m, 1)
         d3(2,j,1,1)=+rdat%r04(j, 3)+rdat%r02(m, 2)
         d3(3,j,1,1)=+rdat%r04(j, 4)+rdat%r02(m, 3)
      enddo

      do j=1,3
         d2(1,j,1,1)=+rdat%r03(j, 6)+rdat%r01(j, 1)
         d2(2,j,1,1)=+rdat%r03(j, 7)+rdat%r01(j, 2)
         d2(3,j,1,1)=+rdat%r03(j, 8)+rdat%r01(j, 3)
         d2(4,j,1,1)=+rdat%r03(j, 9)+rdat%r01(j, 4)
      enddo

         j=1
      do i=1,5
         d1(i  ,1,1)=+rdat%r02(j,i+10)+rdat%r00(i,1)
      enddo
      do j=1,10
         k= ind(j,2,1)
         d4(  j,2,1)=+rdat%r05( k,1)+rdat%r03( j,1)
      enddo

      do j=1,6
         k= ind(j,2,1)
         m= in6(j)
         d3(1,j,2,1)=+rdat%r04(k, 2)+rdat%r02(m, 1)
         d3(2,j,2,1)=+rdat%r04(k, 3)+rdat%r02(m, 2)
         d3(3,j,2,1)=+rdat%r04(k, 4)+rdat%r02(m, 3)
      enddo

      do j=1,3
         k= ind(j,2,1)
         d2(1,j,2,1)=+rdat%r03(k, 6)+rdat%r01(j, 1)
         d2(2,j,2,1)=+rdat%r03(k, 7)+rdat%r01(j, 2)
         d2(3,j,2,1)=+rdat%r03(k, 8)+rdat%r01(j, 3)
         d2(4,j,2,1)=+rdat%r03(k, 9)+rdat%r01(j, 4)
      enddo

         j=1
         k= ind(j,2,1)
         m= in6(k)
      do i=1,5
         d1(i  ,2,1)=+rdat%r02(m,i+10)+rdat%r00(i,1)
      enddo
      d4(  1,3,1)=+rdat%r05( 6, 1)-rdat%r04( 3, 1)* 2 +rdat%r03( 1, 1)+rdat%r03( 1, 2)

      d4(  2,3,1)=+rdat%r05( 9, 1)-rdat%r04( 5, 1)* 2 +rdat%r03( 2, 1)+rdat%r03( 2, 2)
      d4(  3,3,1)=+rdat%r05(10, 1)-rdat%r04( 6, 1)* 2 +rdat%r03( 3, 1)+rdat%r03( 3, 2)

      d4(  4,3,1)=+rdat%r05(13, 1)-rdat%r04( 8, 1)* 2 +rdat%r03( 4, 1)+rdat%r03( 4, 2)
      d4(  5,3,1)=+rdat%r05(14, 1)-rdat%r04( 9, 1)* 2 +rdat%r03( 5, 1)+rdat%r03( 5, 2)
      d4(  6,3,1)=+rdat%r05(15, 1)-rdat%r04(10, 1)* 2 +rdat%r03( 6, 1)+rdat%r03( 6, 2)

      d4(  7,3,1)=+rdat%r05(18, 1)-rdat%r04(12, 1)* 2 +rdat%r03( 7, 1)+rdat%r03( 7, 2)
      d4(  8,3,1)=+rdat%r05(19, 1)-rdat%r04(13, 1)* 2 +rdat%r03( 8, 1)+rdat%r03( 8, 2)
      d4(  9,3,1)=+rdat%r05(20, 1)-rdat%r04(14, 1)* 2 +rdat%r03( 9, 1)+rdat%r03( 9, 2)
      d4( 10,3,1)=+rdat%r05(21, 1)-rdat%r04(15, 1)* 2 +rdat%r03(10, 1)+rdat%r03(10, 2)

      d3(1,1,3,1)=+rdat%r04( 6, 2)-rdat%r03( 3, 3)* 2 +rdat%r02( 1, 1)+rdat%r02( 1, 4)
      d3(2,1,3,1)=+rdat%r04( 6, 3)-rdat%r03( 3, 4)* 2 +rdat%r02( 1, 2)+rdat%r02( 1, 5)
      d3(3,1,3,1)=+rdat%r04( 6, 4)-rdat%r03( 3, 5)* 2 +rdat%r02( 1, 3)+rdat%r02( 1, 6)

      d3(1,2,3,1)=+rdat%r04( 9, 2)-rdat%r03( 5, 3)* 2 +rdat%r02( 4, 1)+rdat%r02( 4, 4)
      d3(2,2,3,1)=+rdat%r04( 9, 3)-rdat%r03( 5, 4)* 2 +rdat%r02( 4, 2)+rdat%r02( 4, 5)
      d3(3,2,3,1)=+rdat%r04( 9, 4)-rdat%r03( 5, 5)* 2 +rdat%r02( 4, 3)+rdat%r02( 4, 6)
      d3(1,3,3,1)=+rdat%r04(10, 2)-rdat%r03( 6, 3)* 2 +rdat%r02( 5, 1)+rdat%r02( 5, 4)
      d3(2,3,3,1)=+rdat%r04(10, 3)-rdat%r03( 6, 4)* 2 +rdat%r02( 5, 2)+rdat%r02( 5, 5)
      d3(3,3,3,1)=+rdat%r04(10, 4)-rdat%r03( 6, 5)* 2 +rdat%r02( 5, 3)+rdat%r02( 5, 6)

      d3(1,4,3,1)=+rdat%r04(13, 2)-rdat%r03( 8, 3)* 2 +rdat%r02( 2, 1)+rdat%r02( 2, 4)
      d3(2,4,3,1)=+rdat%r04(13, 3)-rdat%r03( 8, 4)* 2 +rdat%r02( 2, 2)+rdat%r02( 2, 5)
      d3(3,4,3,1)=+rdat%r04(13, 4)-rdat%r03( 8, 5)* 2 +rdat%r02( 2, 3)+rdat%r02( 2, 6)
      d3(1,5,3,1)=+rdat%r04(14, 2)-rdat%r03( 9, 3)* 2 +rdat%r02( 6, 1)+rdat%r02( 6, 4)
      d3(2,5,3,1)=+rdat%r04(14, 3)-rdat%r03( 9, 4)* 2 +rdat%r02( 6, 2)+rdat%r02( 6, 5)
      d3(3,5,3,1)=+rdat%r04(14, 4)-rdat%r03( 9, 5)* 2 +rdat%r02( 6, 3)+rdat%r02( 6, 6)
      d3(1,6,3,1)=+rdat%r04(15, 2)-rdat%r03(10, 3)* 2 +rdat%r02( 3, 1)+rdat%r02( 3, 4)
      d3(2,6,3,1)=+rdat%r04(15, 3)-rdat%r03(10, 4)* 2 +rdat%r02( 3, 2)+rdat%r02( 3, 5)
      d3(3,6,3,1)=+rdat%r04(15, 4)-rdat%r03(10, 5)* 2 +rdat%r02( 3, 3)+rdat%r02( 3, 6)

      d2(1,1,3,1)=+rdat%r03( 6, 6)-rdat%r02( 5, 7)* 2 +rdat%r01( 1, 1)+rdat%r01( 1, 5)
      d2(2,1,3,1)=+rdat%r03( 6, 7)-rdat%r02( 5, 8)* 2 +rdat%r01( 1, 2)+rdat%r01( 1, 6)
      d2(3,1,3,1)=+rdat%r03( 6, 8)-rdat%r02( 5, 9)* 2 +rdat%r01( 1, 3)+rdat%r01( 1, 7)
      d2(4,1,3,1)=+rdat%r03( 6, 9)-rdat%r02( 5,10)* 2 +rdat%r01( 1, 4)+rdat%r01( 1, 8)

      d2(1,2,3,1)=+rdat%r03( 9, 6)-rdat%r02( 6, 7)* 2 +rdat%r01( 2, 1)+rdat%r01( 2, 5)
      d2(2,2,3,1)=+rdat%r03( 9, 7)-rdat%r02( 6, 8)* 2 +rdat%r01( 2, 2)+rdat%r01( 2, 6)
      d2(3,2,3,1)=+rdat%r03( 9, 8)-rdat%r02( 6, 9)* 2 +rdat%r01( 2, 3)+rdat%r01( 2, 7)
      d2(4,2,3,1)=+rdat%r03( 9, 9)-rdat%r02( 6,10)* 2 +rdat%r01( 2, 4)+rdat%r01( 2, 8)
      d2(1,3,3,1)=+rdat%r03(10, 6)-rdat%r02( 3, 7)* 2 +rdat%r01( 3, 1)+rdat%r01( 3, 5)
      d2(2,3,3,1)=+rdat%r03(10, 7)-rdat%r02( 3, 8)* 2 +rdat%r01( 3, 2)+rdat%r01( 3, 6)
      d2(3,3,3,1)=+rdat%r03(10, 8)-rdat%r02( 3, 9)* 2 +rdat%r01( 3, 3)+rdat%r01( 3, 7)
      d2(4,3,3,1)=+rdat%r03(10, 9)-rdat%r02( 3,10)* 2 +rdat%r01( 3, 4)+rdat%r01( 3, 8)

      do i=1,5
         d1(i  ,3,1)=+rdat%r02(3,i+10)-rdat%r01(3,i+ 8)* 2 +rdat%r00(i,1)+rdat%r00(i,2)
      enddo
      do j=1,10
         k= ind(j,4,1)
         d4(  j,4,1)=+rdat%r05( k,1)
      enddo

      do j=1,6
         k= ind(j,4,1)
         d3(1,j,4,1)=+rdat%r04(k, 2)
         d3(2,j,4,1)=+rdat%r04(k, 3)
         d3(3,j,4,1)=+rdat%r04(k, 4)
      enddo

      do j=1,3
         k= ind(j,4,1)
         d2(1,j,4,1)=+rdat%r03(k, 6)
         d2(2,j,4,1)=+rdat%r03(k, 7)
         d2(3,j,4,1)=+rdat%r03(k, 8)
         d2(4,j,4,1)=+rdat%r03(k, 9)
      enddo

         j=1
         k= ind(j,4,1)
         m= in6(k)
      do i=1,5
         d1(i  ,4,1)=+rdat%r02(m,i+10)
      enddo
      do j=1,10
         k= ind(j,5,1)
         d4(  j,5,1)=+rdat%r05( k,1)-rdat%r04( j,1)
      enddo

      do j=1,6
         k= ind(j,5,1)
         d3(1,j,5,1)=+rdat%r04(k, 2)-rdat%r03(j, 3)
         d3(2,j,5,1)=+rdat%r04(k, 3)-rdat%r03(j, 4)
         d3(3,j,5,1)=+rdat%r04(k, 4)-rdat%r03(j, 5)
      enddo

      do j=1,3
         k= ind(j,5,1)
         m= in6(j)
         d2(1,j,5,1)=+rdat%r03(k, 6)-rdat%r02(m, 7)
         d2(2,j,5,1)=+rdat%r03(k, 7)-rdat%r02(m, 8)
         d2(3,j,5,1)=+rdat%r03(k, 8)-rdat%r02(m, 9)
         d2(4,j,5,1)=+rdat%r03(k, 9)-rdat%r02(m,10)
      enddo

         j=1
         k= ind(j,5,1)
         m= in6(k)
      do i=1,5
         d1(i  ,5,1)=+rdat%r02(m,i+10)-rdat%r01(j,i+ 8)
      enddo
      d4(  1,6,1)=+rdat%r05( 5, 1)-rdat%r04( 2, 1)

      d4(  2,6,1)=+rdat%r05( 8, 1)-rdat%r04( 4, 1)
      d4(  3,6,1)=+rdat%r05( 9, 1)-rdat%r04( 5, 1)

      d4(  4,6,1)=+rdat%r05(12, 1)-rdat%r04( 7, 1)
      d4(  5,6,1)=+rdat%r05(13, 1)-rdat%r04( 8, 1)
      d4(  6,6,1)=+rdat%r05(14, 1)-rdat%r04( 9, 1)

      d4(  7,6,1)=+rdat%r05(17, 1)-rdat%r04(11, 1)
      d4(  8,6,1)=+rdat%r05(18, 1)-rdat%r04(12, 1)
      d4(  9,6,1)=+rdat%r05(19, 1)-rdat%r04(13, 1)
      d4( 10,6,1)=+rdat%r05(20, 1)-rdat%r04(14, 1)

      d3(1,1,6,1)=+rdat%r04( 5, 2)-rdat%r03( 2, 3)
      d3(2,1,6,1)=+rdat%r04( 5, 3)-rdat%r03( 2, 4)
      d3(3,1,6,1)=+rdat%r04( 5, 4)-rdat%r03( 2, 5)

      d3(1,2,6,1)=+rdat%r04( 8, 2)-rdat%r03( 4, 3)
      d3(2,2,6,1)=+rdat%r04( 8, 3)-rdat%r03( 4, 4)
      d3(3,2,6,1)=+rdat%r04( 8, 4)-rdat%r03( 4, 5)
      d3(1,3,6,1)=+rdat%r04( 9, 2)-rdat%r03( 5, 3)
      d3(2,3,6,1)=+rdat%r04( 9, 3)-rdat%r03( 5, 4)
      d3(3,3,6,1)=+rdat%r04( 9, 4)-rdat%r03( 5, 5)

      d3(1,4,6,1)=+rdat%r04(12, 2)-rdat%r03( 7, 3)
      d3(2,4,6,1)=+rdat%r04(12, 3)-rdat%r03( 7, 4)
      d3(3,4,6,1)=+rdat%r04(12, 4)-rdat%r03( 7, 5)
      d3(1,5,6,1)=+rdat%r04(13, 2)-rdat%r03( 8, 3)
      d3(2,5,6,1)=+rdat%r04(13, 3)-rdat%r03( 8, 4)
      d3(3,5,6,1)=+rdat%r04(13, 4)-rdat%r03( 8, 5)
      d3(1,6,6,1)=+rdat%r04(14, 2)-rdat%r03( 9, 3)
      d3(2,6,6,1)=+rdat%r04(14, 3)-rdat%r03( 9, 4)
      d3(3,6,6,1)=+rdat%r04(14, 4)-rdat%r03( 9, 5)

      d2(1,1,6,1)=+rdat%r03( 5, 6)-rdat%r02( 4, 7)
      d2(2,1,6,1)=+rdat%r03( 5, 7)-rdat%r02( 4, 8)
      d2(3,1,6,1)=+rdat%r03( 5, 8)-rdat%r02( 4, 9)
      d2(4,1,6,1)=+rdat%r03( 5, 9)-rdat%r02( 4,10)

      d2(1,2,6,1)=+rdat%r03( 8, 6)-rdat%r02( 2, 7)
      d2(2,2,6,1)=+rdat%r03( 8, 7)-rdat%r02( 2, 8)
      d2(3,2,6,1)=+rdat%r03( 8, 8)-rdat%r02( 2, 9)
      d2(4,2,6,1)=+rdat%r03( 8, 9)-rdat%r02( 2,10)
      d2(1,3,6,1)=+rdat%r03( 9, 6)-rdat%r02( 6, 7)
      d2(2,3,6,1)=+rdat%r03( 9, 7)-rdat%r02( 6, 8)
      d2(3,3,6,1)=+rdat%r03( 9, 8)-rdat%r02( 6, 9)
      d2(4,3,6,1)=+rdat%r03( 9, 9)-rdat%r02( 6,10)

      do i=1,5
         d1(i  ,6,1)=+rdat%r02(6,i+10)-rdat%r01(2,i+ 8)
      enddo

      xx= qx*qx
      zz= qz*qz
      xz= qx*qz
      xxx= xx*qx
      xxz= xx*qz
      xzz= zz*qx
      zzz= zz*qz

      qxd= qx+qx
      qzd= qz+qz

      xzd= xz+xz

      l=1
      do k=1,kx

        f(1,1,k,l) = d4(  1,k,l)+d2(2,1,k,l)* 3 +(+d3(2,1,k,l)* 2 +d3(3,1,k,l) +d1(  3,k,l)* 2 +d1(  4,k,l))*qx +(+d2(3,1,k,l)+d2(4,1,k,l)* 2 )*xx +d1(  5,k,l)*xxx
        f(2,1,k,l) = d4(  4,k,l)+d2(2,1,k,l) +(+d3(3,4,k,l)+d1(  4,k,l))*qx
        f(3,1,k,l) = d4(  6,k,l)+d2(2,1,k,l) +(+d3(3,6,k,l)+d1(  4,k,l))*qx +d3(2,3,k,l)*qzd +d2(4,3,k,l)*xzd +d2(3,1,k,l)*zz +d1(  5,k,l)*xzz
        f(4,1,k,l) = d4(  2,k,l)+d2(2,2,k,l) +(+d3(2,2,k,l)+d3(3,2,k,l))*qx +d2(4,2,k,l)*xx
        f(5,1,k,l) = d4(  3,k,l)+d2(2,3,k,l) +(+d3(2,3,k,l)+d3(3,3,k,l))*qx +(+d3(2,1,k,l)+d1(  3,k,l))*qz +d2(4,3,k,l)*xx +(+d2(3,1,k,l)+d2(4,1,k,l))*xz +d1(  5,k,l)*xxz
        f(6,1,k,l) = d4(  5,k,l) +d3(3,5,k,l)*qx +d3(2,2,k,l)*qz +d2(4,2,k,l)*xz

        f(1,2,k,l) = d4(  2,k,l)+d2(2,2,k,l) +d3(2,2,k,l)*qxd +d2(3,2,k,l)*xx
        f(2,2,k,l) = d4(  7,k,l)+d2(2,2,k,l)* 3
        f(3,2,k,l) = d4(  9,k,l)+d2(2,2,k,l) +d3(2,5,k,l)*qzd +d2(3,2,k,l)*zz
        f(4,2,k,l) = d4(  4,k,l)+d2(2,1,k,l) +(+d3(2,4,k,l)+d1(  3,k,l))*qx
        f(5,2,k,l) = d4(  5,k,l) +d3(2,5,k,l)*qx +d3(2,2,k,l)*qz +d2(3,2,k,l)*xz
        f(6,2,k,l) = d4(  8,k,l)+d2(2,3,k,l) +(+d3(2,4,k,l)+d1(  3,k,l))*qz

        f(1,3,k,l) = d4(  3,k,l)+d2(2,3,k,l) +d3(2,3,k,l)*qxd +(+d3(3,1,k,l)+d1(  4,k,l))*qz +d2(3,3,k,l)*xx +d2(4,1,k,l)*xzd +d1(  5,k,l)*xxz
        f(2,3,k,l) = d4(  8,k,l)+d2(2,3,k,l) +(+d3(3,4,k,l)+d1(  4,k,l))*qz
        f(3,3,k,l) = d4( 10,k,l)+d2(2,3,k,l)* 3 +(+d3(2,6,k,l)* 2 +d3(3,6,k,l) +d1(  3,k,l)* 2 +d1(  4,k,l))*qz +(+d2(3,3,k,l)+d2(4,3,k,l)* 2 )*zz +d1(  5,k,l)*zzz
        f(4,3,k,l) = d4(  5,k,l) +d3(2,5,k,l)*qx +d3(3,2,k,l)*qz +d2(4,2,k,l)*xz
        f(5,3,k,l) = d4(  6,k,l)+d2(2,1,k,l) +(+d3(2,6,k,l)+d1(  3,k,l))*qx +(+d3(2,3,k,l)+d3(3,3,k,l))*qz +(+d2(3,3,k,l)+d2(4,3,k,l))*xz +d2(4,1,k,l)*zz +d1(  5,k,l)*xzz
        f(6,3,k,l) = d4(  9,k,l)+d2(2,2,k,l) +(+d3(2,5,k,l)+d3(3,5,k,l))*qz +d2(4,2,k,l)*zz

      enddo

      end subroutine mcdv_15

! >
! >    @brief   dppp case
! >
! >    @details integration of a dppp case
! >
      subroutine mcdv_16(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,3,3,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 4, lx = 4
      real(kind=dp) :: d1(  5,kx,lx),d2(4,3,kx,lx),d3(3,6,kx,lx), &
                       d4( 10,kx,lx)

      integer, parameter :: ind (10, kx, lx) = reshape(&
        [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, &
          4, 5, 7, 8, 9, 11, 12, 13, 14, 3, 5, 6, 8, 9, 10, 12, 13, 14, &
          15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, &
          10, 2, 4, 5, 7, 8, 9, 11, 12, 13, 14, 3, 5, 6, 8, 9, 10, 12, 13, &
          14, 15, 2, 4, 5, 7, 8, 9, 11, 12, 13, 14, 2, 4, 5, 7, 8, 9, 11, &
          12, 13, 14, 4, 7, 8, 11, 12, 13, 16, 17, 18, 19, 5, 8, 9, 12, &
          13, 14, 17, 18, 19, 20, 3, 5, 6, 8, 9, 10, 12, 13, 14, 15, 3, 5, &
          6, 8, 9, 10, 12, 13, 14, 15, 5, 8, 9, 12, 13, 14, 17, 18, 19, &
          20, 6, 9, 10, 13, 14, 15, 18, 19, 20, 21 ] &
        , shape(ind))
      integer :: i, j, k, l, m
      logical :: lsym16, lsym19

      real(kind=dp) :: xx, zz, xz, xxx, xxz, xzz, zzz
      real(kind=dp) :: qxd, qzd, xzd

      do j=1,10
         d4(  j,1,1)=+rdat%r03( j,1)
      enddo

      do j=1,6
         m= in6(j)
         d3(1,j,1,1)=+rdat%r02(m, 1)
         d3(2,j,1,1)=+rdat%r02(m, 2)
         d3(3,j,1,1)=+rdat%r02(m, 3)
      enddo

      do j=1,3
         d2(1,j,1,1)=+rdat%r01(j, 1)
         d2(2,j,1,1)=+rdat%r01(j, 2)
         d2(3,j,1,1)=+rdat%r01(j, 3)
         d2(4,j,1,1)=+rdat%r01(j, 4)
      enddo

      do i=1,5
         d1(i  ,1,1)=+rdat%r00(i,1)
      enddo
      do j=1,10
         d4(  j,2,1)=-rdat%r04( j,2)
      enddo

      do j=1,6
         d3(1,j,2,1)=-rdat%r03(j, 9)
         d3(2,j,2,1)=-rdat%r03(j,10)
         d3(3,j,2,1)=-rdat%r03(j,11)
      enddo

      do j=1,3
         m= in6(j)
         d2(1,j,2,1)=-rdat%r02(m,16)
         d2(2,j,2,1)=-rdat%r02(m,17)
         d2(3,j,2,1)=-rdat%r02(m,18)
         d2(4,j,2,1)=-rdat%r02(m,19)
      enddo

         j=1
      do i=1,5
         d1(i  ,2,1)=-rdat%r01(j,i+20)
      enddo
      do j=1,10
         k= ind(j,3,1)
         d4(  j,3,1)=-rdat%r04( k,2)
      enddo

      do j=1,6
         k= ind(j,3,1)
         d3(1,j,3,1)=-rdat%r03(k, 9)
         d3(2,j,3,1)=-rdat%r03(k,10)
         d3(3,j,3,1)=-rdat%r03(k,11)
      enddo

      do j=1,3
         k= ind(j,3,1)
         m= in6(k)
         d2(1,j,3,1)=-rdat%r02(m,16)
         d2(2,j,3,1)=-rdat%r02(m,17)
         d2(3,j,3,1)=-rdat%r02(m,18)
         d2(4,j,3,1)=-rdat%r02(m,19)
      enddo

         j=1
         k= ind(j,3,1)
      do i=1,5
         d1(i  ,3,1)=-rdat%r01(k,i+20)
      enddo
      do j=1,10
         k= ind(j,4,1)
         d4(  j,4,1)=-rdat%r04( k,2)+rdat%r03( j,2)
      enddo

      do j=1,6
         k= ind(j,4,1)
         m= in6(j)
         d3(1,j,4,1)=-rdat%r03(k, 9)+rdat%r02(m, 4)
         d3(2,j,4,1)=-rdat%r03(k,10)+rdat%r02(m, 5)
         d3(3,j,4,1)=-rdat%r03(k,11)+rdat%r02(m, 6)
      enddo

      do j=1,3
         k= ind(j,4,1)
         m= in6(k)
         d2(1,j,4,1)=-rdat%r02(m,16)+rdat%r01(j, 5)
         d2(2,j,4,1)=-rdat%r02(m,17)+rdat%r01(j, 6)
         d2(3,j,4,1)=-rdat%r02(m,18)+rdat%r01(j, 7)
         d2(4,j,4,1)=-rdat%r02(m,19)+rdat%r01(j, 8)
      enddo

         j=1
         k= ind(j,4,1)
      do i=1,5
         d1(i  ,4,1)=-rdat%r01(k,i+20)+rdat%r00(i,2)
      enddo
      do j=1,10
         d4(  j,1,2)=-rdat%r04( j,3)
      enddo

      do j=1,6
         d3(1,j,1,2)=-rdat%r03(j,12)
         d3(2,j,1,2)=-rdat%r03(j,13)
         d3(3,j,1,2)=-rdat%r03(j,14)
      enddo

      do j=1,3
         m= in6(j)
         d2(1,j,1,2)=-rdat%r02(m,20)
         d2(2,j,1,2)=-rdat%r02(m,21)
         d2(3,j,1,2)=-rdat%r02(m,22)
         d2(4,j,1,2)=-rdat%r02(m,23)
      enddo

         j=1
      do i=1,5
         d1(i  ,1,2)=-rdat%r01(j,i+25)
      enddo
      do j=1,10
         d4(  j,2,2)=+rdat%r05( j,1)+rdat%r03( j,5)
      enddo

      do j=1,6
         m= in6(j)
         d3(1,j,2,2)=+rdat%r04(j, 5)+rdat%r02(m,13)
         d3(2,j,2,2)=+rdat%r04(j, 6)+rdat%r02(m,14)
         d3(3,j,2,2)=+rdat%r04(j, 7)+rdat%r02(m,15)
      enddo

      do j=1,3
         d2(1,j,2,2)=+rdat%r03(j,18)+rdat%r01(j,17)
         d2(2,j,2,2)=+rdat%r03(j,19)+rdat%r01(j,18)
         d2(3,j,2,2)=+rdat%r03(j,20)+rdat%r01(j,19)
         d2(4,j,2,2)=+rdat%r03(j,21)+rdat%r01(j,20)
      enddo

         j=1
      do i=1,5
         d1(i  ,2,2)=+rdat%r02(j,i+31)+rdat%r00(i,5)
      enddo
      do j=1,10
         k= ind(j,3,2)
         d4(  j,3,2)=+rdat%r05( k,1)
      enddo

      do j=1,6
         k= ind(j,3,2)
         d3(1,j,3,2)=+rdat%r04(k, 5)
         d3(2,j,3,2)=+rdat%r04(k, 6)
         d3(3,j,3,2)=+rdat%r04(k, 7)
      enddo

      do j=1,3
         k= ind(j,3,2)
         d2(1,j,3,2)=+rdat%r03(k,18)
         d2(2,j,3,2)=+rdat%r03(k,19)
         d2(3,j,3,2)=+rdat%r03(k,20)
         d2(4,j,3,2)=+rdat%r03(k,21)
      enddo

         j=1
         k= ind(j,3,2)
         m= in6(k)
      do i=1,5
         d1(i  ,3,2)=+rdat%r02(m,i+31)
      enddo
      do j=1,10
         k= ind(j,4,2)
         d4(  j,4,2)=+rdat%r05( k,1)-rdat%r04( j,4)
      enddo

      do j=1,6
         k= ind(j,4,2)
         d3(1,j,4,2)=+rdat%r04(k, 5)-rdat%r03(j,15)
         d3(2,j,4,2)=+rdat%r04(k, 6)-rdat%r03(j,16)
         d3(3,j,4,2)=+rdat%r04(k, 7)-rdat%r03(j,17)
      enddo

      do j=1,3
         k= ind(j,4,2)
         m= in6(j)
         d2(1,j,4,2)=+rdat%r03(k,18)-rdat%r02(m,24)
         d2(2,j,4,2)=+rdat%r03(k,19)-rdat%r02(m,25)
         d2(3,j,4,2)=+rdat%r03(k,20)-rdat%r02(m,26)
         d2(4,j,4,2)=+rdat%r03(k,21)-rdat%r02(m,27)
      enddo

         j=1
         k= ind(j,4,2)
         m= in6(k)
      do i=1,5
         d1(i  ,4,2)=+rdat%r02(m,i+31)-rdat%r01(j,i+30)
      enddo
      do j=1,10
         k= ind(j,1,3)
         d4(  j,1,3)=-rdat%r04( k,3)
      enddo

      do j=1,6
         k= ind(j,1,3)
         d3(1,j,1,3)=-rdat%r03(k,12)
         d3(2,j,1,3)=-rdat%r03(k,13)
         d3(3,j,1,3)=-rdat%r03(k,14)
      enddo

      do j=1,3
         k= ind(j,1,3)
         m= in6(k)
         d2(1,j,1,3)=-rdat%r02(m,20)
         d2(2,j,1,3)=-rdat%r02(m,21)
         d2(3,j,1,3)=-rdat%r02(m,22)
         d2(4,j,1,3)=-rdat%r02(m,23)
      enddo

         j=1
         k= ind(j,1,3)
      do i=1,5
         d1(i  ,1,3)=-rdat%r01(k,i+25)
      enddo
      do j=1,10
         k= ind(j,3,3)
         d4(  j,3,3)=+rdat%r05( k,1)+rdat%r03( j,5)
      enddo

      do j=1,6
         k= ind(j,3,3)
         m= in6(j)
         d3(1,j,3,3)=+rdat%r04(k, 5)+rdat%r02(m,13)
         d3(2,j,3,3)=+rdat%r04(k, 6)+rdat%r02(m,14)
         d3(3,j,3,3)=+rdat%r04(k, 7)+rdat%r02(m,15)
      enddo

      do j=1,3
         k= ind(j,3,3)
         d2(1,j,3,3)=+rdat%r03(k,18)+rdat%r01(j,17)
         d2(2,j,3,3)=+rdat%r03(k,19)+rdat%r01(j,18)
         d2(3,j,3,3)=+rdat%r03(k,20)+rdat%r01(j,19)
         d2(4,j,3,3)=+rdat%r03(k,21)+rdat%r01(j,20)
      enddo

         j=1
         k= ind(j,3,3)
         m= in6(k)
      do i=1,5
         d1(i  ,3,3)=+rdat%r02(m,i+31)+rdat%r00(i,5)
      enddo
      d4(  1,4,3)=+rdat%r05( 5, 1)-rdat%r04( 2, 4)

      d4(  2,4,3)=+rdat%r05( 8, 1)-rdat%r04( 4, 4)
      d4(  3,4,3)=+rdat%r05( 9, 1)-rdat%r04( 5, 4)

      d4(  4,4,3)=+rdat%r05(12, 1)-rdat%r04( 7, 4)
      d4(  5,4,3)=+rdat%r05(13, 1)-rdat%r04( 8, 4)
      d4(  6,4,3)=+rdat%r05(14, 1)-rdat%r04( 9, 4)

      d4(  7,4,3)=+rdat%r05(17, 1)-rdat%r04(11, 4)
      d4(  8,4,3)=+rdat%r05(18, 1)-rdat%r04(12, 4)
      d4(  9,4,3)=+rdat%r05(19, 1)-rdat%r04(13, 4)
      d4( 10,4,3)=+rdat%r05(20, 1)-rdat%r04(14, 4)

      d3(1,1,4,3)=+rdat%r04( 5, 5)-rdat%r03( 2,15)
      d3(2,1,4,3)=+rdat%r04( 5, 6)-rdat%r03( 2,16)
      d3(3,1,4,3)=+rdat%r04( 5, 7)-rdat%r03( 2,17)

      d3(1,2,4,3)=+rdat%r04( 8, 5)-rdat%r03( 4,15)
      d3(2,2,4,3)=+rdat%r04( 8, 6)-rdat%r03( 4,16)
      d3(3,2,4,3)=+rdat%r04( 8, 7)-rdat%r03( 4,17)
      d3(1,3,4,3)=+rdat%r04( 9, 5)-rdat%r03( 5,15)
      d3(2,3,4,3)=+rdat%r04( 9, 6)-rdat%r03( 5,16)
      d3(3,3,4,3)=+rdat%r04( 9, 7)-rdat%r03( 5,17)

      d3(1,4,4,3)=+rdat%r04(12, 5)-rdat%r03( 7,15)
      d3(2,4,4,3)=+rdat%r04(12, 6)-rdat%r03( 7,16)
      d3(3,4,4,3)=+rdat%r04(12, 7)-rdat%r03( 7,17)
      d3(1,5,4,3)=+rdat%r04(13, 5)-rdat%r03( 8,15)
      d3(2,5,4,3)=+rdat%r04(13, 6)-rdat%r03( 8,16)
      d3(3,5,4,3)=+rdat%r04(13, 7)-rdat%r03( 8,17)
      d3(1,6,4,3)=+rdat%r04(14, 5)-rdat%r03( 9,15)
      d3(2,6,4,3)=+rdat%r04(14, 6)-rdat%r03( 9,16)
      d3(3,6,4,3)=+rdat%r04(14, 7)-rdat%r03( 9,17)

      d2(1,1,4,3)=+rdat%r03( 5,18)-rdat%r02( 4,24)
      d2(2,1,4,3)=+rdat%r03( 5,19)-rdat%r02( 4,25)
      d2(3,1,4,3)=+rdat%r03( 5,20)-rdat%r02( 4,26)
      d2(4,1,4,3)=+rdat%r03( 5,21)-rdat%r02( 4,27)

      d2(1,2,4,3)=+rdat%r03( 8,18)-rdat%r02( 2,24)
      d2(2,2,4,3)=+rdat%r03( 8,19)-rdat%r02( 2,25)
      d2(3,2,4,3)=+rdat%r03( 8,20)-rdat%r02( 2,26)
      d2(4,2,4,3)=+rdat%r03( 8,21)-rdat%r02( 2,27)
      d2(1,3,4,3)=+rdat%r03( 9,18)-rdat%r02( 6,24)
      d2(2,3,4,3)=+rdat%r03( 9,19)-rdat%r02( 6,25)
      d2(3,3,4,3)=+rdat%r03( 9,20)-rdat%r02( 6,26)
      d2(4,3,4,3)=+rdat%r03( 9,21)-rdat%r02( 6,27)

      do i=1,5
         d1(i  ,4,3)=+rdat%r02(6,i+31)-rdat%r01(2,i+30)
      enddo
      do j=1,10
         k= ind(j,1,4)
         d4(  j,1,4)=-rdat%r04( k,3)+rdat%r03( j,3)
      enddo

      do j=1,6
         k= ind(j,1,4)
         m= in6(j)
         d3(1,j,1,4)=-rdat%r03(k,12)+rdat%r02(m, 7)
         d3(2,j,1,4)=-rdat%r03(k,13)+rdat%r02(m, 8)
         d3(3,j,1,4)=-rdat%r03(k,14)+rdat%r02(m, 9)
      enddo

      do j=1,3
         k= ind(j,1,4)
         m= in6(k)
         d2(1,j,1,4)=-rdat%r02(m,20)+rdat%r01(j, 9)
         d2(2,j,1,4)=-rdat%r02(m,21)+rdat%r01(j,10)
         d2(3,j,1,4)=-rdat%r02(m,22)+rdat%r01(j,11)
         d2(4,j,1,4)=-rdat%r02(m,23)+rdat%r01(j,12)
      enddo

         j=1
         k= ind(j,1,4)
      do i=1,5
         d1(i  ,1,4)=-rdat%r01(k,i+25)+rdat%r00(i,3)
      enddo
      do j=1,10
         k= ind(j,2,4)
         d4(  j,2,4)=+rdat%r05( k,1)-rdat%r04( j,1)
      enddo

      do j=1,6
         k= ind(j,2,4)
         d3(1,j,2,4)=+rdat%r04(k, 5)-rdat%r03(j, 6)
         d3(2,j,2,4)=+rdat%r04(k, 6)-rdat%r03(j, 7)
         d3(3,j,2,4)=+rdat%r04(k, 7)-rdat%r03(j, 8)
      enddo

      do j=1,3
         k= ind(j,2,4)
         m= in6(j)
         d2(1,j,2,4)=+rdat%r03(k,18)-rdat%r02(m,28)
         d2(2,j,2,4)=+rdat%r03(k,19)-rdat%r02(m,29)
         d2(3,j,2,4)=+rdat%r03(k,20)-rdat%r02(m,30)
         d2(4,j,2,4)=+rdat%r03(k,21)-rdat%r02(m,31)
      enddo

         j=1
         k= ind(j,2,4)
         m= in6(k)
      do i=1,5
         d1(i  ,2,4)=+rdat%r02(m,i+31)-rdat%r01(j,i+35)
      enddo
      d4(  1,3,4)=+rdat%r05( 5, 1)-rdat%r04( 2, 1)

      d4(  2,3,4)=+rdat%r05( 8, 1)-rdat%r04( 4, 1)
      d4(  3,3,4)=+rdat%r05( 9, 1)-rdat%r04( 5, 1)

      d4(  4,3,4)=+rdat%r05(12, 1)-rdat%r04( 7, 1)
      d4(  5,3,4)=+rdat%r05(13, 1)-rdat%r04( 8, 1)
      d4(  6,3,4)=+rdat%r05(14, 1)-rdat%r04( 9, 1)

      d4(  7,3,4)=+rdat%r05(17, 1)-rdat%r04(11, 1)
      d4(  8,3,4)=+rdat%r05(18, 1)-rdat%r04(12, 1)
      d4(  9,3,4)=+rdat%r05(19, 1)-rdat%r04(13, 1)
      d4( 10,3,4)=+rdat%r05(20, 1)-rdat%r04(14, 1)

      d3(1,1,3,4)=+rdat%r04( 5, 5)-rdat%r03( 2, 6)
      d3(2,1,3,4)=+rdat%r04( 5, 6)-rdat%r03( 2, 7)
      d3(3,1,3,4)=+rdat%r04( 5, 7)-rdat%r03( 2, 8)

      d3(1,2,3,4)=+rdat%r04( 8, 5)-rdat%r03( 4, 6)
      d3(2,2,3,4)=+rdat%r04( 8, 6)-rdat%r03( 4, 7)
      d3(3,2,3,4)=+rdat%r04( 8, 7)-rdat%r03( 4, 8)
      d3(1,3,3,4)=+rdat%r04( 9, 5)-rdat%r03( 5, 6)
      d3(2,3,3,4)=+rdat%r04( 9, 6)-rdat%r03( 5, 7)
      d3(3,3,3,4)=+rdat%r04( 9, 7)-rdat%r03( 5, 8)

      d3(1,4,3,4)=+rdat%r04(12, 5)-rdat%r03( 7, 6)
      d3(2,4,3,4)=+rdat%r04(12, 6)-rdat%r03( 7, 7)
      d3(3,4,3,4)=+rdat%r04(12, 7)-rdat%r03( 7, 8)
      d3(1,5,3,4)=+rdat%r04(13, 5)-rdat%r03( 8, 6)
      d3(2,5,3,4)=+rdat%r04(13, 6)-rdat%r03( 8, 7)
      d3(3,5,3,4)=+rdat%r04(13, 7)-rdat%r03( 8, 8)
      d3(1,6,3,4)=+rdat%r04(14, 5)-rdat%r03( 9, 6)
      d3(2,6,3,4)=+rdat%r04(14, 6)-rdat%r03( 9, 7)
      d3(3,6,3,4)=+rdat%r04(14, 7)-rdat%r03( 9, 8)

      d2(1,1,3,4)=+rdat%r03( 5,18)-rdat%r02( 4,28)
      d2(2,1,3,4)=+rdat%r03( 5,19)-rdat%r02( 4,29)
      d2(3,1,3,4)=+rdat%r03( 5,20)-rdat%r02( 4,30)
      d2(4,1,3,4)=+rdat%r03( 5,21)-rdat%r02( 4,31)

      d2(1,2,3,4)=+rdat%r03( 8,18)-rdat%r02( 2,28)
      d2(2,2,3,4)=+rdat%r03( 8,19)-rdat%r02( 2,29)
      d2(3,2,3,4)=+rdat%r03( 8,20)-rdat%r02( 2,30)
      d2(4,2,3,4)=+rdat%r03( 8,21)-rdat%r02( 2,31)
      d2(1,3,3,4)=+rdat%r03( 9,18)-rdat%r02( 6,28)
      d2(2,3,3,4)=+rdat%r03( 9,19)-rdat%r02( 6,29)
      d2(3,3,3,4)=+rdat%r03( 9,20)-rdat%r02( 6,30)
      d2(4,3,3,4)=+rdat%r03( 9,21)-rdat%r02( 6,31)

      do i=1,5
         d1(i  ,3,4)=+rdat%r02(6,i+31)-rdat%r01(2,i+35)
      enddo
      d4(  1,4,4)=+rdat%r05( 6, 1)-rdat%r04( 3, 1)-rdat%r04( 3, 4)+rdat%r03( 1, 4)+rdat%r03( 1, 5)

      d4(  2,4,4)=+rdat%r05( 9, 1)-rdat%r04( 5, 1)-rdat%r04( 5, 4)+rdat%r03( 2, 4)+rdat%r03( 2, 5)
      d4(  3,4,4)=+rdat%r05(10, 1)-rdat%r04( 6, 1)-rdat%r04( 6, 4)+rdat%r03( 3, 4)+rdat%r03( 3, 5)

      d4(  4,4,4)=+rdat%r05(13, 1)-rdat%r04( 8, 1)-rdat%r04( 8, 4)+rdat%r03( 4, 4)+rdat%r03( 4, 5)
      d4(  5,4,4)=+rdat%r05(14, 1)-rdat%r04( 9, 1)-rdat%r04( 9, 4)+rdat%r03( 5, 4)+rdat%r03( 5, 5)
      d4(  6,4,4)=+rdat%r05(15, 1)-rdat%r04(10, 1)-rdat%r04(10, 4)+rdat%r03( 6, 4)+rdat%r03( 6, 5)

      d4(  7,4,4)=+rdat%r05(18, 1)-rdat%r04(12, 1)-rdat%r04(12, 4)+rdat%r03( 7, 4)+rdat%r03( 7, 5)
      d4(  8,4,4)=+rdat%r05(19, 1)-rdat%r04(13, 1)-rdat%r04(13, 4)+rdat%r03( 8, 4)+rdat%r03( 8, 5)
      d4(  9,4,4)=+rdat%r05(20, 1)-rdat%r04(14, 1)-rdat%r04(14, 4)+rdat%r03( 9, 4)+rdat%r03( 9, 5)
      d4( 10,4,4)=+rdat%r05(21, 1)-rdat%r04(15, 1)-rdat%r04(15, 4)+rdat%r03(10, 4)+rdat%r03(10, 5)

      d3(1,1,4,4)=+rdat%r04( 6, 5)-rdat%r03( 3, 6)-rdat%r03( 3,15)+rdat%r02( 1,10)+rdat%r02( 1,13)
      d3(2,1,4,4)=+rdat%r04( 6, 6)-rdat%r03( 3, 7)-rdat%r03( 3,16)+rdat%r02( 1,11)+rdat%r02( 1,14)
      d3(3,1,4,4)=+rdat%r04( 6, 7)-rdat%r03( 3, 8)-rdat%r03( 3,17)+rdat%r02( 1,12)+rdat%r02( 1,15)

      d3(1,2,4,4)=+rdat%r04( 9, 5)-rdat%r03( 5, 6)-rdat%r03( 5,15)+rdat%r02( 4,10)+rdat%r02( 4,13)
      d3(2,2,4,4)=+rdat%r04( 9, 6)-rdat%r03( 5, 7)-rdat%r03( 5,16)+rdat%r02( 4,11)+rdat%r02( 4,14)
      d3(3,2,4,4)=+rdat%r04( 9, 7)-rdat%r03( 5, 8)-rdat%r03( 5,17)+rdat%r02( 4,12)+rdat%r02( 4,15)
      d3(1,3,4,4)=+rdat%r04(10, 5)-rdat%r03( 6, 6)-rdat%r03( 6,15)+rdat%r02( 5,10)+rdat%r02( 5,13)
      d3(2,3,4,4)=+rdat%r04(10, 6)-rdat%r03( 6, 7)-rdat%r03( 6,16)+rdat%r02( 5,11)+rdat%r02( 5,14)
      d3(3,3,4,4)=+rdat%r04(10, 7)-rdat%r03( 6, 8)-rdat%r03( 6,17)+rdat%r02( 5,12)+rdat%r02( 5,15)

      d3(1,4,4,4)=+rdat%r04(13, 5)-rdat%r03( 8, 6)-rdat%r03( 8,15)+rdat%r02( 2,10)+rdat%r02( 2,13)
      d3(2,4,4,4)=+rdat%r04(13, 6)-rdat%r03( 8, 7)-rdat%r03( 8,16)+rdat%r02( 2,11)+rdat%r02( 2,14)
      d3(3,4,4,4)=+rdat%r04(13, 7)-rdat%r03( 8, 8)-rdat%r03( 8,17)+rdat%r02( 2,12)+rdat%r02( 2,15)
      d3(1,5,4,4)=+rdat%r04(14, 5)-rdat%r03( 9, 6)-rdat%r03( 9,15)+rdat%r02( 6,10)+rdat%r02( 6,13)
      d3(2,5,4,4)=+rdat%r04(14, 6)-rdat%r03( 9, 7)-rdat%r03( 9,16)+rdat%r02( 6,11)+rdat%r02( 6,14)
      d3(3,5,4,4)=+rdat%r04(14, 7)-rdat%r03( 9, 8)-rdat%r03( 9,17)+rdat%r02( 6,12)+rdat%r02( 6,15)
      d3(1,6,4,4)=+rdat%r04(15, 5)-rdat%r03(10, 6)-rdat%r03(10,15)+rdat%r02( 3,10)+rdat%r02( 3,13)
      d3(2,6,4,4)=+rdat%r04(15, 6)-rdat%r03(10, 7)-rdat%r03(10,16)+rdat%r02( 3,11)+rdat%r02( 3,14)
      d3(3,6,4,4)=+rdat%r04(15, 7)-rdat%r03(10, 8)-rdat%r03(10,17)+rdat%r02( 3,12)+rdat%r02( 3,15)

      d2(1,1,4,4)=+rdat%r03( 6,18)-rdat%r02( 5,24)-rdat%r02( 5,28)+rdat%r01( 1,13)+rdat%r01( 1,17)
      d2(2,1,4,4)=+rdat%r03( 6,19)-rdat%r02( 5,25)-rdat%r02( 5,29)+rdat%r01( 1,14)+rdat%r01( 1,18)
      d2(3,1,4,4)=+rdat%r03( 6,20)-rdat%r02( 5,26)-rdat%r02( 5,30)+rdat%r01( 1,15)+rdat%r01( 1,19)
      d2(4,1,4,4)=+rdat%r03( 6,21)-rdat%r02( 5,27)-rdat%r02( 5,31)+rdat%r01( 1,16)+rdat%r01( 1,20)

      d2(1,2,4,4)=+rdat%r03( 9,18)-rdat%r02( 6,24)-rdat%r02( 6,28)+rdat%r01( 2,13)+rdat%r01( 2,17)
      d2(2,2,4,4)=+rdat%r03( 9,19)-rdat%r02( 6,25)-rdat%r02( 6,29)+rdat%r01( 2,14)+rdat%r01( 2,18)
      d2(3,2,4,4)=+rdat%r03( 9,20)-rdat%r02( 6,26)-rdat%r02( 6,30)+rdat%r01( 2,15)+rdat%r01( 2,19)
      d2(4,2,4,4)=+rdat%r03( 9,21)-rdat%r02( 6,27)-rdat%r02( 6,31)+rdat%r01( 2,16)+rdat%r01( 2,20)
      d2(1,3,4,4)=+rdat%r03(10,18)-rdat%r02( 3,24)-rdat%r02( 3,28)+rdat%r01( 3,13)+rdat%r01( 3,17)
      d2(2,3,4,4)=+rdat%r03(10,19)-rdat%r02( 3,25)-rdat%r02( 3,29)+rdat%r01( 3,14)+rdat%r01( 3,18)
      d2(3,3,4,4)=+rdat%r03(10,20)-rdat%r02( 3,26)-rdat%r02( 3,30)+rdat%r01( 3,15)+rdat%r01( 3,19)
      d2(4,3,4,4)=+rdat%r03(10,21)-rdat%r02( 3,27)-rdat%r02( 3,31)+rdat%r01( 3,16)+rdat%r01( 3,20)

      do i=1,5
         d1(i  ,4,4)=+rdat%r02(3,i+31)-rdat%r01(3,i+30)-rdat%r01(3,i+35)+rdat%r00(i,4)+rdat%r00(i,5)
      enddo

      xx= qx*qx
      zz= qz*qz
      xz= qx*qz
      xxx= xx*qx
      xxz= xx*qz
      xzz= zz*qx
      zzz= zz*qz

      qxd= qx+qx
      qzd= qz+qz

      xzd= xz+xz

      do l=2,lx
         do k=2,kx

           if(k == 2 .and. l == 3) cycle

            f(1,1,k-1,l-1) = d4(  1,k,l)+d2(2,1,k,l)* 3 +(+d3(2,1,k,l)* 2 +d3(3,1,k,l) +d1(  3,k,l)* 2 +d1(  4,k,l))*qx +(+d2(3,1,k,l)+d2(4,1,k,l)* 2 )*xx +d1(  5,k,l)*xxx
            f(2,1,k-1,l-1) = d4(  4,k,l)+d2(2,1,k,l) +(+d3(3,4,k,l)+d1(  4,k,l))*qx
            f(3,1,k-1,l-1) = d4(  6,k,l)+d2(2,1,k,l) +(+d3(3,6,k,l)+d1(  4,k,l))*qx +d3(2,3,k,l)*qzd +d2(4,3,k,l)*xzd +d2(3,1,k,l)*zz +d1(  5,k,l)*xzz
            f(4,1,k-1,l-1) = d4(  2,k,l)+d2(2,2,k,l) +(+d3(2,2,k,l)+d3(3,2,k,l))*qx +d2(4,2,k,l)*xx
            f(5,1,k-1,l-1) = d4(  3,k,l)+d2(2,3,k,l) +(+d3(2,3,k,l)+d3(3,3,k,l))*qx +(+d3(2,1,k,l)+d1(  3,k,l))*qz +d2(4,3,k,l)*xx +(+d2(3,1,k,l)+d2(4,1,k,l))*xz +d1(  5,k,l)*xxz
            f(6,1,k-1,l-1) = d4(  5,k,l) +d3(3,5,k,l)*qx +d3(2,2,k,l)*qz +d2(4,2,k,l)*xz

            f(1,2,k-1,l-1) = d4(  2,k,l)+d2(2,2,k,l) +d3(2,2,k,l)*qxd +d2(3,2,k,l)*xx
            f(2,2,k-1,l-1) = d4(  7,k,l)+d2(2,2,k,l)* 3
            f(3,2,k-1,l-1) = d4(  9,k,l)+d2(2,2,k,l) +d3(2,5,k,l)*qzd +d2(3,2,k,l)*zz
            f(4,2,k-1,l-1) = d4(  4,k,l)+d2(2,1,k,l) +(+d3(2,4,k,l)+d1(  3,k,l))*qx
            f(5,2,k-1,l-1) = d4(  5,k,l) +d3(2,5,k,l)*qx +d3(2,2,k,l)*qz +d2(3,2,k,l)*xz
            f(6,2,k-1,l-1) = d4(  8,k,l)+d2(2,3,k,l) +(+d3(2,4,k,l)+d1(  3,k,l))*qz

            f(1,3,k-1,l-1) = d4(  3,k,l)+d2(2,3,k,l) +d3(2,3,k,l)*qxd +(+d3(3,1,k,l)+d1(  4,k,l))*qz +d2(3,3,k,l)*xx +d2(4,1,k,l)*xzd +d1(  5,k,l)*xxz
            f(2,3,k-1,l-1) = d4(  8,k,l)+d2(2,3,k,l) +(+d3(3,4,k,l)+d1(  4,k,l))*qz
            f(3,3,k-1,l-1) = d4( 10,k,l)+d2(2,3,k,l)* 3 +(+d3(2,6,k,l)* 2 +d3(3,6,k,l) +d1(  3,k,l)* 2 +d1(  4,k,l))*qz +(+d2(3,3,k,l)+d2(4,3,k,l)* 2 )*zz +d1(  5,k,l)*zzz
            f(4,3,k-1,l-1) = d4(  5,k,l) +d3(2,5,k,l)*qx +d3(3,2,k,l)*qz +d2(4,2,k,l)*xz
            f(5,3,k-1,l-1) = d4(  6,k,l)+d2(2,1,k,l) +(+d3(2,6,k,l)+d1(  3,k,l))*qx +(+d3(2,3,k,l)+d3(3,3,k,l))*qz +(+d2(3,3,k,l)+d2(4,3,k,l))*xz +d2(4,1,k,l)*zz +d1(  5,k,l)*xzz
            f(6,3,k-1,l-1) = d4(  9,k,l)+d2(2,2,k,l) +(+d3(2,5,k,l)+d3(3,5,k,l))*qz +d2(4,2,k,l)*zz

         enddo
      enddo

      f(:,:,1,2)= f(:,:,2,1)

      end subroutine mcdv_16

! >
! >    @brief   ddds case
! >
! >    @details integration of a ddds case
! >
      subroutine mcdv_17(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,6,6,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 6, lx = 1
      real(kind=dp) :: e1(  5,kx,lx),e2(4,3,kx,lx),e3(4,6,kx,lx), &
                       e4(2,10,kx,lx),e5( 15,kx,lx)

      integer, parameter :: ind (15, kx, lx) = reshape(&
        [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 4, 7, 8, 11, &
          12, 13, 16, 17, 18, 19, 22, 23, 24, 25, 26, 6, 9, 10, 13, 14, &
          15, 18, 19, 20, 21, 24, 25, 26, 27, 28, 2, 4, 5, 7, 8, 9, 11, &
          12, 13, 14, 16, 17, 18, 19, 20, 3, 5, 6, 8, 9, 10, 12, 13, 14, &
          15, 17, 18, 19, 20, 21, 5, 8, 9, 12, 13, 14, 17, 18, 19, 20, 23, &
          24, 25, 26, 27] &
        , shape(ind))
      integer :: i, j, k, l, m, ii, jj
      real(kind=dp) :: xx, zz, xz, xxx, xxz, xzz, zzz
      real(kind=dp) :: xxxx, xxzz, xxxz, xzzz, zzzz
      real(kind=dp) :: xxxd, xxzd, xzzd, zzzd
      real(kind=dp) :: xzq
      real(kind=dp) :: qxd, qzd, xzd

      do j=1,15
         e5(  j,1,1)=+rdat%r06( j,1)+rdat%r04( j,1)
      enddo

      do j=1,10
         e4(1,j,1,1)=+rdat%r05(j, 2)+rdat%r03(j, 1)
         e4(2,j,1,1)=+rdat%r05(j, 3)+rdat%r03(j, 2)
      enddo

      do j=1,6
         m= in6(j)
         e3(1,j,1,1)=+rdat%r04(j, 5)+rdat%r02(m, 1)
         e3(2,j,1,1)=+rdat%r04(j, 6)+rdat%r02(m, 2)
         e3(3,j,1,1)=+rdat%r04(j, 7)+rdat%r02(m, 3)
         e3(4,j,1,1)=+rdat%r04(j, 8)+rdat%r02(m, 4)
      enddo

      do j=1,3
         e2(1,j,1,1)=+rdat%r03(j, 9)+rdat%r01(j, 1)
         e2(2,j,1,1)=+rdat%r03(j,10)+rdat%r01(j, 2)
         e2(3,j,1,1)=+rdat%r03(j,11)+rdat%r01(j, 3)
         e2(4,j,1,1)=+rdat%r03(j,12)+rdat%r01(j, 4)
      enddo

         j=1
      do i=1,5
         e1(i  ,1,1)=+rdat%r02(j,i+12)+rdat%r00(i,1)
      enddo
      do j=1,15
         k= ind(j,2,1)
         e5(  j,2,1)=+rdat%r06( k,1)+rdat%r04( j,1)
      enddo

      do j=1,10
         k= ind(j,2,1)
         e4(1,j,2,1)=+rdat%r05(k, 2)+rdat%r03(j, 1)
         e4(2,j,2,1)=+rdat%r05(k, 3)+rdat%r03(j, 2)
      enddo

      do j=1,6
         k= ind(j,2,1)
         m= in6(j)
         e3(1,j,2,1)=+rdat%r04(k, 5)+rdat%r02(m, 1)
         e3(2,j,2,1)=+rdat%r04(k, 6)+rdat%r02(m, 2)
         e3(3,j,2,1)=+rdat%r04(k, 7)+rdat%r02(m, 3)
         e3(4,j,2,1)=+rdat%r04(k, 8)+rdat%r02(m, 4)
      enddo

      do j=1,3
         k= ind(j,2,1)
         e2(1,j,2,1)=+rdat%r03(k, 9)+rdat%r01(j, 1)
         e2(2,j,2,1)=+rdat%r03(k,10)+rdat%r01(j, 2)
         e2(3,j,2,1)=+rdat%r03(k,11)+rdat%r01(j, 3)
         e2(4,j,2,1)=+rdat%r03(k,12)+rdat%r01(j, 4)
      enddo

         j=1
         k= ind(j,2,1)
         m= in6(k)
      do i=1,5
         e1(i  ,2,1)=+rdat%r02(m,i+12)+rdat%r00(i,1)
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,3,1)
            l= k-2-jj

      e5(  j,3,1)=+rdat%r06(k,1)-rdat%r05(l,1)* 2 +rdat%r04(j,1)+rdat%r04(j,2)

            if(jj > 4) cycle
      e4(1,j,3,1)=+rdat%r05(k,2)-rdat%r04(l,3)* 2 +rdat%r03(j,1)+rdat%r03(j,3)
      e4(2,j,3,1)=+rdat%r05(k,3)-rdat%r04(l,4)* 2 +rdat%r03(j,2)+rdat%r03(j,4)

            if(jj > 3) cycle
            m= in6(j)
            do i=1,4
      e3(i,j,3,1)=+rdat%r04(k,i+ 4)-rdat%r03(l,i+ 4)* 2 +rdat%r02(m,i   )+rdat%r02(m,i+ 4)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      e2(i,j,3,1)=+rdat%r03(k,i+ 8)-rdat%r02(m,i+ 8)* 2 +rdat%r01(j,i   )+rdat%r01(j,i+ 4)
            enddo

            if(jj > 1) cycle
            m= in6(k)
            do i=1,5
      e1(i  ,3,1)=+rdat%r02(m,i+12)-rdat%r01(l,i+ 8)* 2 +rdat%r00(i,1)+rdat%r00(i,2)
            enddo
         enddo
      enddo
      do j=1,15
         k= ind(j,4,1)
         e5(  j,4,1)=+rdat%r06( k,1)
      enddo

      do j=1,10
         k= ind(j,4,1)
         e4(1,j,4,1)=+rdat%r05(k, 2)
         e4(2,j,4,1)=+rdat%r05(k, 3)
      enddo

      do j=1,6
         k= ind(j,4,1)
         e3(1,j,4,1)=+rdat%r04(k, 5)
         e3(2,j,4,1)=+rdat%r04(k, 6)
         e3(3,j,4,1)=+rdat%r04(k, 7)
         e3(4,j,4,1)=+rdat%r04(k, 8)
      enddo

      do j=1,3
         k= ind(j,4,1)
         e2(1,j,4,1)=+rdat%r03(k, 9)
         e2(2,j,4,1)=+rdat%r03(k,10)
         e2(3,j,4,1)=+rdat%r03(k,11)
         e2(4,j,4,1)=+rdat%r03(k,12)
      enddo

         j=1
         k= ind(j,4,1)
         m= in6(k)
      do i=1,5
         e1(i  ,4,1)=+rdat%r02(m,i+12)
      enddo
      do j=1,15
         k= ind(j,5,1)
         e5(  j,5,1)=+rdat%r06( k,1)-rdat%r05( j,1)
      enddo

      do j=1,10
         k= ind(j,5,1)
         e4(1,j,5,1)=+rdat%r05(k, 2)-rdat%r04(j, 3)
         e4(2,j,5,1)=+rdat%r05(k, 3)-rdat%r04(j, 4)
      enddo

      do j=1,6
         k= ind(j,5,1)
         e3(1,j,5,1)=+rdat%r04(k, 5)-rdat%r03(j, 5)
         e3(2,j,5,1)=+rdat%r04(k, 6)-rdat%r03(j, 6)
         e3(3,j,5,1)=+rdat%r04(k, 7)-rdat%r03(j, 7)
         e3(4,j,5,1)=+rdat%r04(k, 8)-rdat%r03(j, 8)
      enddo

      do j=1,3
         k= ind(j,5,1)
         m= in6(j)
         e2(1,j,5,1)=+rdat%r03(k, 9)-rdat%r02(m, 9)
         e2(2,j,5,1)=+rdat%r03(k,10)-rdat%r02(m,10)
         e2(3,j,5,1)=+rdat%r03(k,11)-rdat%r02(m,11)
         e2(4,j,5,1)=+rdat%r03(k,12)-rdat%r02(m,12)
      enddo

         j=1
         k= ind(j,5,1)
         m= in6(k)
      do i=1,5
         e1(i  ,5,1)=+rdat%r02(m,i+12)-rdat%r01(j,i+ 8)
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,6,1)
            l= k-2-jj

      e5(  j,6,1)=+rdat%r06(k,1)-rdat%r05(l,1)

            if(jj > 4) cycle
      e4(1,j,6,1)=+rdat%r05(k,2)-rdat%r04(l,3)
      e4(2,j,6,1)=+rdat%r05(k,3)-rdat%r04(l,4)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,6,1)=+rdat%r04(k,i+ 4)-rdat%r03(l,i+ 4)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      e2(i,j,6,1)=+rdat%r03(k,i+ 8)-rdat%r02(m,i+ 8)
            enddo

            if(jj > 1) cycle
            m= in6(k)
            do i=1,5
      e1(i  ,6,1)=+rdat%r02(m,i+12)-rdat%r01(l,i+ 8)
            enddo
         enddo
      enddo

      xx= qx*qx
      zz= qz*qz
      xz= qx*qz
      xxx= xx*qx
      xxz= xx*qz
      xzz= zz*qx
      zzz= zz*qz
      xxxx= xx*xx
      xxxz= xx*xz
      xxzz= xx*zz
      xzzz= xz*zz
      zzzz= zz*zz

      qxd= qx+qx
      qzd= qz+qz
      xzd= xz+xz
      xxxd= xxx+xxx
      xxzd= xxz+xxz
      xzzd= xzz+xzz
      zzzd= zzz+zzz

      xzq= xz+xz+xz+xz

      do k=1,kx
         f(1,1,k,1) = e5(  1,k,1)+e3(1,1,k,1)* 6 +e1(  1,k,1)* 3 +(+e4(1,1,k,1)+e4(2,1,k,1)+ (+e2(1,1,k,1)+e2(2,1,k,1))* 3 )*qxd +(+e3(2,1,k,1)+e3(3,1,k,1)* 4 +e3(4,1,k,1) +e1(  2,k,1)+e1(  3,k,1)* 4 +e1(  4,k,1))*xx +(+e2(3,1,k,1)+e2(4,1,k,1))*xxxd +e1(  5,k,1)*xxxx
         f(2,1,k,1) = e5(  4,k,1)+e3(1,4,k,1)+e3(1,1,k,1)+e1(  1,k,1) +(+e4(2,4,k,1)+e2(2,1,k,1))*qxd +(+e3(4,4,k,1)+e1(  4,k,1))*xx
         f(3,1,k,1) = e5(  6,k,1)+e3(1,6,k,1)+e3(1,1,k,1)+e1(  1,k,1) +(+e4(2,6,k,1)+e2(2,1,k,1))*qxd +(+e4(1,3,k,1)+e2(1,3,k,1))*qzd +(+e3(4,6,k,1)+e1(  4,k,1))*xx +e3(3,3,k,1)*xzq +(+e3(2,1,k,1)+e1(  2,k,1))*zz +e2(4,3,k,1)*xxzd +e2(3,1,k,1)*xzzd +e1(  5,k,1)*xxzz
         f(4,1,k,1) = e5(  2,k,1)+e3(1,2,k,1)* 3 +(+e4(1,2,k,1)+e4(2,2,k,1)* 2 +e2(1,2,k,1)+e2(2,2,k,1)* 2 )*qx +(+e3(3,2,k,1)* 2 +e3(4,2,k,1))*xx +e2(4,2,k,1)*xxx
         f(5,1,k,1) = e5(  3,k,1)+e3(1,3,k,1)* 3 +(+e4(1,3,k,1)+e4(2,3,k,1)* 2 +e2(1,3,k,1)+e2(2,3,k,1)* 2 )*qx +(+e4(1,1,k,1)+e2(1,1,k,1)* 3 )*qz +(+e3(3,3,k,1)* 2 +e3(4,3,k,1))*xx +(+e3(2,1,k,1)+e3(3,1,k,1)* 2 +e1(  2,k,1)+e1(  3,k,1)* 2 )*xz +e2(4,3,k,1)*xxx +(+e2(3,1,k,1)* 2 +e2(4,1,k,1))*xxz +e1(  5,k,1)*xxxz
         f(6,1,k,1) = e5(  5,k,1)+e3(1,5,k,1) +e4(2,5,k,1)*qxd +(+e4(1,2,k,1)+e2(1,2,k,1))*qz +e3(4,5,k,1)*xx +e3(3,2,k,1)*xzd +e2(4,2,k,1)*xxz

         f(1,2,k,1) = e5(  4,k,1)+e3(1,4,k,1)+e3(1,1,k,1)+e1(  1,k,1) +(+e4(1,4,k,1)+e2(1,1,k,1))*qxd +(+e3(2,4,k,1)+e1(  2,k,1))*xx
         f(2,2,k,1) = e5( 11,k,1)+e3(1,4,k,1)* 6 +e1(  1,k,1)* 3
         f(3,2,k,1) = e5( 13,k,1)+e3(1,6,k,1)+e3(1,4,k,1)+e1(  1,k,1) +(+e4(1,8,k,1)+e2(1,3,k,1))*qzd +(+e3(2,4,k,1)+e1(  2,k,1))*zz
         f(4,2,k,1) = e5(  7,k,1)+e3(1,2,k,1)* 3 +(+e4(1,7,k,1)+e2(1,2,k,1)* 3 )*qx
         f(5,2,k,1) = e5(  8,k,1)+e3(1,3,k,1) +(+e4(1,8,k,1)+e2(1,3,k,1))*qx +(+e4(1,4,k,1)+e2(1,1,k,1))*qz +(+e3(2,4,k,1)+e1(  2,k,1))*xz
         f(6,2,k,1) = e5( 12,k,1)+e3(1,5,k,1)* 3 +(+e4(1,7,k,1)+e2(1,2,k,1)* 3 )*qz

         f(1,3,k,1) = e5(  6,k,1)+e3(1,6,k,1)+e3(1,1,k,1)+e1(  1,k,1) +(+e4(1,6,k,1)+e2(1,1,k,1))*qxd +(+e4(2,3,k,1)+e2(2,3,k,1))*qzd +(+e3(2,6,k,1)+e1(  2,k,1))*xx +e3(3,3,k,1)*xzq +(+e3(4,1,k,1)+e1(  4,k,1))*zz +e2(3,3,k,1)*xxzd +e2(4,1,k,1)*xzzd +e1(  5,k,1)*xxzz
         f(2,3,k,1) = e5( 13,k,1)+e3(1,6,k,1)+e3(1,4,k,1)+e1(  1,k,1) +(+e4(2,8,k,1)+e2(2,3,k,1))*qzd +(+e3(4,4,k,1)+e1(  4,k,1))*zz
         f(3,3,k,1) = e5( 15,k,1)+e3(1,6,k,1)* 6 +e1(  1,k,1)* 3 +(+e4(1,10,k,1)+e4(2,10,k,1)+ (+e2(1,3,k,1)+e2(2,3,k,1))* 3 )*qzd +(+e3(2,6,k,1)+e3(3,6,k,1)* 4 +e3(4,6,k,1) +e1(  2,k,1)+e1(  3,k,1)* 4 +e1(  4,k,1))*zz +(+e2(3,3,k,1)+e2(4,3,k,1))*zzzd +e1(  5,k,1)*zzzz
         f(4,3,k,1) = e5(  9,k,1)+e3(1,2,k,1) +(+e4(1,9,k,1)+e2(1,2,k,1))*qx +e4(2,5,k,1)*qzd +e3(3,5,k,1)*xzd +e3(4,2,k,1)*zz +e2(4,2,k,1)*xzz
         f(5,3,k,1) = e5( 10,k,1)+e3(1,3,k,1)* 3 +(+e4(1,10,k,1)+e2(1,3,k,1)* 3 )*qx +(+e4(1,6,k,1)+e4(2,6,k,1)* 2 +e2(1,1,k,1)+e2(2,1,k,1)* 2 )*qz +(+e3(2,6,k,1)+e3(3,6,k,1)* 2 +e1(  2,k,1)+e1(  3,k,1)* 2 )*xz +(+e3(3,3,k,1)* 2 +e3(4,3,k,1))*zz +(+e2(3,3,k,1)* 2 +e2(4,3,k,1))*xzz +e2(4,1,k,1)*zzz +e1(  5,k,1)*xzzz
         f(6,3,k,1) = e5( 14,k,1)+e3(1,5,k,1)* 3 +(+e4(1,9,k,1)+e4(2,9,k,1)* 2 +e2(1,2,k,1)+e2(2,2,k,1)* 2 )*qz +(+e3(3,5,k,1)* 2 +e3(4,5,k,1))*zz +e2(4,2,k,1)*zzz

         f(1,4,k,1) = e5(  2,k,1)+e3(1,2,k,1)* 3 +(+e4(1,2,k,1)* 2 +e4(2,2,k,1) +e2(1,2,k,1)* 2 +e2(2,2,k,1))*qx +(+e3(2,2,k,1)+e3(3,2,k,1)* 2 )*xx +e2(3,2,k,1)*xxx
         f(2,4,k,1) = e5(  7,k,1)+e3(1,2,k,1)* 3 +(+e4(2,7,k,1)+e2(2,2,k,1)* 3 )*qx
         f(3,4,k,1) = e5(  9,k,1)+e3(1,2,k,1) +(+e4(2,9,k,1)+e2(2,2,k,1))*qx +e4(1,5,k,1)*qzd +e3(3,5,k,1)*xzd +e3(2,2,k,1)*zz +e2(3,2,k,1)*xzz
         f(4,4,k,1) = e5(  4,k,1)+e3(1,4,k,1)+e3(1,1,k,1)+e1(  1,k,1) +(+e4(1,4,k,1)+e4(2,4,k,1) +e2(1,1,k,1)+e2(2,1,k,1))*qx +(+e3(3,4,k,1)+e1(  3,k,1))*xx
         f(5,4,k,1) = e5(  5,k,1)+e3(1,5,k,1) +(+e4(1,5,k,1)+e4(2,5,k,1))*qx +(+e4(1,2,k,1)+e2(1,2,k,1))*qz +e3(3,5,k,1)*xx +(+e3(2,2,k,1)+e3(3,2,k,1))*xz +e2(3,2,k,1)*xxz
         f(6,4,k,1) = e5(  8,k,1)+e3(1,3,k,1) +(+e4(2,8,k,1)+e2(2,3,k,1))*qx +(+e4(1,4,k,1)+e2(1,1,k,1))*qz +(+e3(3,4,k,1)+e1(  3,k,1))*xz

         f(1,5,k,1) = e5(  3,k,1)+e3(1,3,k,1)* 3 +(+e4(1,3,k,1)* 2 +e4(2,3,k,1) +e2(1,3,k,1)* 2 +e2(2,3,k,1))*qx +(+e4(2,1,k,1)+e2(2,1,k,1)* 3 )*qz +(+e3(2,3,k,1)+e3(3,3,k,1)* 2 )*xx +(+e3(3,1,k,1)* 2 +e3(4,1,k,1) +e1(  3,k,1)* 2 +e1(  4,k,1))*xz +e2(3,3,k,1)*xxx +(+e2(3,1,k,1)+e2(4,1,k,1)* 2 )*xxz +e1(  5,k,1)*xxxz
         f(2,5,k,1) = e5(  8,k,1)+e3(1,3,k,1) +(+e4(2,8,k,1)+e2(2,3,k,1))*qx +(+e4(2,4,k,1)+e2(2,1,k,1))*qz +(+e3(4,4,k,1)+e1(  4,k,1))*xz
         f(3,5,k,1) = e5( 10,k,1)+e3(1,3,k,1)* 3 +(+e4(2,10,k,1)+e2(2,3,k,1)* 3 )*qx +(+e4(1,6,k,1)* 2 +e4(2,6,k,1) +e2(1,1,k,1)* 2 +e2(2,1,k,1))*qz +(+e3(3,6,k,1)* 2 +e3(4,6,k,1) +e1(  3,k,1)* 2 +e1(  4,k,1))*xz +(+e3(2,3,k,1)+e3(3,3,k,1)* 2 )*zz +(+e2(3,3,k,1)+e2(4,3,k,1)* 2 )*xzz +e2(3,1,k,1)*zzz +e1(  5,k,1)*xzzz
         f(4,5,k,1) = e5(  5,k,1)+e3(1,5,k,1) +(+e4(1,5,k,1)+e4(2,5,k,1))*qx +(+e4(2,2,k,1)+e2(2,2,k,1))*qz +e3(3,5,k,1)*xx +(+e3(3,2,k,1)+e3(4,2,k,1))*xz +e2(4,2,k,1)*xxz
         f(5,5,k,1) = e5(  6,k,1)+e3(1,6,k,1)+e3(1,1,k,1)+e1(  1,k,1) +(+e4(1,6,k,1)+e4(2,6,k,1) +e2(1,1,k,1)+e2(2,1,k,1))*qx +(+e4(1,3,k,1)+e4(2,3,k,1) +e2(1,3,k,1)+e2(2,3,k,1))*qz +(+e3(3,6,k,1)+e1(  3,k,1))*xx +(+e3(2,3,k,1)+e3(3,3,k,1)* 2 +e3(4,3,k,1))*xz +(+e3(3,1,k,1)+e1(  3,k,1))*zz +(+e2(3,3,k,1)+e2(4,3,k,1))*xxz +(+e2(3,1,k,1)+e2(4,1,k,1))*xzz +e1(  5,k,1)*xxzz
         f(6,5,k,1) = e5(  9,k,1)+e3(1,2,k,1) +(+e4(2,9,k,1)+e2(2,2,k,1))*qx +(+e4(1,5,k,1)+e4(2,5,k,1))*qz +(+e3(3,5,k,1)+e3(4,5,k,1))*xz +e3(3,2,k,1)*zz +e2(4,2,k,1)*xzz

         f(1,6,k,1) = e5(  5,k,1)+e3(1,5,k,1) +e4(1,5,k,1)*qxd +(+e4(2,2,k,1)+e2(2,2,k,1))*qz +e3(2,5,k,1)*xx +e3(3,2,k,1)*xzd +e2(3,2,k,1)*xxz
         f(2,6,k,1) = e5( 12,k,1)+e3(1,5,k,1)* 3 +(+e4(2,7,k,1)+e2(2,2,k,1)* 3 )*qz
         f(3,6,k,1) = e5( 14,k,1)+e3(1,5,k,1)* 3 +(+e4(1,9,k,1)* 2 +e4(2,9,k,1) +e2(1,2,k,1)* 2 +e2(2,2,k,1))*qz +(+e3(2,5,k,1)+e3(3,5,k,1)* 2 )*zz +e2(3,2,k,1)*zzz
         f(4,6,k,1) = e5(  8,k,1)+e3(1,3,k,1) +(+e4(1,8,k,1)+e2(1,3,k,1))*qx +(+e4(2,4,k,1)+e2(2,1,k,1))*qz +(+e3(3,4,k,1)+e1(  3,k,1))*xz
         f(5,6,k,1) = e5(  9,k,1)+e3(1,2,k,1) +(+e4(1,9,k,1)+e2(1,2,k,1))*qx +(+e4(1,5,k,1)+e4(2,5,k,1))*qz +(+e3(2,5,k,1)+e3(3,5,k,1))*xz +e3(3,2,k,1)*zz +e2(3,2,k,1)*xzz
         f(6,6,k,1) = e5( 13,k,1)+e3(1,6,k,1)+e3(1,4,k,1)+e1(  1,k,1) +(+e4(1,8,k,1)+e4(2,8,k,1) +e2(1,3,k,1)+e2(2,3,k,1))*qz +(+e3(3,4,k,1)+e1(  3,k,1))*zz

      enddo

      end subroutine mcdv_17

! >
! >    @brief   ddpp case
! >
! >    @details integration of a ddpp case
! >
      subroutine mcdv_18(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,6,3,3)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 4, lx = 4
      real(kind=dp) :: e1(  5,kx,lx),e2(4,3,kx,lx),e3(4,6,kx,lx), &
                       e4(2,10,kx,lx),e5( 15,kx,lx)

      integer, parameter :: ind (15, kx, lx) = reshape(&
        [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1, 2, 3, 4, &
          5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 2, 4, 5, 7, 8, 9, 11, 12, &
          13, 14, 16, 17, 18, 19, 20, 3, 5, 6, 8, 9, 10, 12, 13, 14, 15, &
          17, 18, 19, 20, 21, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, &
          14, 15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 2, 4, &
          5, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19, 20, 3, 5, 6, 8, 9, &
          10, 12, 13, 14, 15, 17, 18, 19, 20, 21, 2, 4, 5, 7, 8, 9, 11, &
          12, 13, 14, 16, 17, 18, 19, 20, 2, 4, 5, 7, 8, 9, 11, 12, 13, &
          14, 16, 17, 18, 19, 20, 4, 7, 8, 11, 12, 13, 16, 17, 18, 19, 22, &
          23, 24, 25, 26, 5, 8, 9, 12, 13, 14, 17, 18, 19, 20, 23, 24, 25, &
          26, 27, 3, 5, 6, 8, 9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21, &
          3, 5, 6, 8, 9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21, 5, 8, 9, &
          12, 13, 14, 17, 18, 19, 20, 23, 24, 25, 26, 27, 6, 9, 10, 13, &
          14, 15, 18, 19, 20, 21, 24, 25, 26, 27, 28] &
        , shape(ind))
      integer :: i, j, k, l, m, ii, jj
      real(kind=dp) :: xx, zz, xz, xxx, xxz, xzz, zzz
      real(kind=dp) :: xxxx, xxzz, xxxz, xzzz, zzzz
      real(kind=dp) :: xxxd, xxzd, xzzd, zzzd
      real(kind=dp) :: xzq
      real(kind=dp) :: qxd, qzd, xzd

      do j=1,15
         e5(  j,1,1)=+rdat%r04( j,1)
      enddo

      do j=1,10
         e4(1,j,1,1)=+rdat%r03(j, 1)
         e4(2,j,1,1)=+rdat%r03(j, 2)
      enddo

      do j=1,6
         m= in6(j)
         e3(1,j,1,1)=+rdat%r02(m, 1)
         e3(2,j,1,1)=+rdat%r02(m, 2)
         e3(3,j,1,1)=+rdat%r02(m, 3)
         e3(4,j,1,1)=+rdat%r02(m, 4)
      enddo

      do j=1,3
         e2(1,j,1,1)=+rdat%r01(j, 1)
         e2(2,j,1,1)=+rdat%r01(j, 2)
         e2(3,j,1,1)=+rdat%r01(j, 3)
         e2(4,j,1,1)=+rdat%r01(j, 4)
      enddo

      do i=1,5
         e1(i  ,1,1)=+rdat%r00(i,1)
      enddo
      do j=1,15
         e5(  j,2,1)=-rdat%r05( j,2)
      enddo

      do j=1,10
         e4(1,j,2,1)=-rdat%r04(j, 8)
         e4(2,j,2,1)=-rdat%r04(j, 9)
      enddo

      do j=1,6
         e3(1,j,2,1)=-rdat%r03(j,11)
         e3(2,j,2,1)=-rdat%r03(j,12)
         e3(3,j,2,1)=-rdat%r03(j,13)
         e3(4,j,2,1)=-rdat%r03(j,14)
      enddo

      do j=1,3
         m= in6(j)
         e2(1,j,2,1)=-rdat%r02(m,21)
         e2(2,j,2,1)=-rdat%r02(m,22)
         e2(3,j,2,1)=-rdat%r02(m,23)
         e2(4,j,2,1)=-rdat%r02(m,24)
      enddo

         j=1
      do i=1,5
         e1(i  ,2,1)=-rdat%r01(j,i+20)
      enddo
      do j=1,15
         k= ind(j,3,1)
         e5(  j,3,1)=-rdat%r05( k,2)
      enddo

      do j=1,10
         k= ind(j,3,1)
         e4(1,j,3,1)=-rdat%r04(k, 8)
         e4(2,j,3,1)=-rdat%r04(k, 9)
      enddo

      do j=1,6
         k= ind(j,3,1)
         e3(1,j,3,1)=-rdat%r03(k,11)
         e3(2,j,3,1)=-rdat%r03(k,12)
         e3(3,j,3,1)=-rdat%r03(k,13)
         e3(4,j,3,1)=-rdat%r03(k,14)
      enddo

      do j=1,3
         k= ind(j,3,1)
         m= in6(k)
         e2(1,j,3,1)=-rdat%r02(m,21)
         e2(2,j,3,1)=-rdat%r02(m,22)
         e2(3,j,3,1)=-rdat%r02(m,23)
         e2(4,j,3,1)=-rdat%r02(m,24)
      enddo

         j=1
         k= ind(j,3,1)
      do i=1,5
         e1(i  ,3,1)=-rdat%r01(k,i+20)
      enddo
      do j=1,15
         k= ind(j,4,1)
         e5(  j,4,1)=-rdat%r05( k,2)+rdat%r04( j,2)
      enddo

      do j=1,10
         k= ind(j,4,1)
         e4(1,j,4,1)=-rdat%r04(k, 8)+rdat%r03(j, 3)
         e4(2,j,4,1)=-rdat%r04(k, 9)+rdat%r03(j, 4)
      enddo

      do j=1,6
         k= ind(j,4,1)
         m= in6(j)
         e3(1,j,4,1)=-rdat%r03(k,11)+rdat%r02(m, 5)
         e3(2,j,4,1)=-rdat%r03(k,12)+rdat%r02(m, 6)
         e3(3,j,4,1)=-rdat%r03(k,13)+rdat%r02(m, 7)
         e3(4,j,4,1)=-rdat%r03(k,14)+rdat%r02(m, 8)
      enddo

      do j=1,3
         k= ind(j,4,1)
         m= in6(k)
         e2(1,j,4,1)=-rdat%r02(m,21)+rdat%r01(j, 5)
         e2(2,j,4,1)=-rdat%r02(m,22)+rdat%r01(j, 6)
         e2(3,j,4,1)=-rdat%r02(m,23)+rdat%r01(j, 7)
         e2(4,j,4,1)=-rdat%r02(m,24)+rdat%r01(j, 8)
      enddo

         j=1
         k= ind(j,4,1)
      do i=1,5
         e1(i  ,4,1)=-rdat%r01(k,i+20)+rdat%r00(i,2)
      enddo
      do j=1,15
         e5(  j,1,2)=-rdat%r05( j,3)
      enddo

      do j=1,10
         e4(1,j,1,2)=-rdat%r04(j,10)
         e4(2,j,1,2)=-rdat%r04(j,11)
      enddo

      do j=1,6
         e3(1,j,1,2)=-rdat%r03(j,15)
         e3(2,j,1,2)=-rdat%r03(j,16)
         e3(3,j,1,2)=-rdat%r03(j,17)
         e3(4,j,1,2)=-rdat%r03(j,18)
      enddo

      do j=1,3
         m= in6(j)
         e2(1,j,1,2)=-rdat%r02(m,25)
         e2(2,j,1,2)=-rdat%r02(m,26)
         e2(3,j,1,2)=-rdat%r02(m,27)
         e2(4,j,1,2)=-rdat%r02(m,28)
      enddo

         j=1
      do i=1,5
         e1(i  ,1,2)=-rdat%r01(j,i+25)
      enddo
      do j=1,15
         e5(  j,2,2)=+rdat%r06( j,1)+rdat%r04( j,5)
      enddo

      do j=1,10
         e4(1,j,2,2)=+rdat%r05(j, 5)+rdat%r03(j, 9)
         e4(2,j,2,2)=+rdat%r05(j, 6)+rdat%r03(j,10)
      enddo

      do j=1,6
         m= in6(j)
         e3(1,j,2,2)=+rdat%r04(j,14)+rdat%r02(m,17)
         e3(2,j,2,2)=+rdat%r04(j,15)+rdat%r02(m,18)
         e3(3,j,2,2)=+rdat%r04(j,16)+rdat%r02(m,19)
         e3(4,j,2,2)=+rdat%r04(j,17)+rdat%r02(m,20)
      enddo

      do j=1,3
         e2(1,j,2,2)=+rdat%r03(j,27)+rdat%r01(j,17)
         e2(2,j,2,2)=+rdat%r03(j,28)+rdat%r01(j,18)
         e2(3,j,2,2)=+rdat%r03(j,29)+rdat%r01(j,19)
         e2(4,j,2,2)=+rdat%r03(j,30)+rdat%r01(j,20)
      enddo

         j=1
      do i=1,5
         e1(i  ,2,2)=+rdat%r02(j,i+36)+rdat%r00(i,5)
      enddo
      do j=1,15
         k= ind(j,3,2)
         e5(  j,3,2)=+rdat%r06( k,1)
      enddo

      do j=1,10
         k= ind(j,3,2)
         e4(1,j,3,2)=+rdat%r05(k, 5)
         e4(2,j,3,2)=+rdat%r05(k, 6)
      enddo

      do j=1,6
         k= ind(j,3,2)
         e3(1,j,3,2)=+rdat%r04(k,14)
         e3(2,j,3,2)=+rdat%r04(k,15)
         e3(3,j,3,2)=+rdat%r04(k,16)
         e3(4,j,3,2)=+rdat%r04(k,17)
      enddo

      do j=1,3
         k= ind(j,3,2)
         e2(1,j,3,2)=+rdat%r03(k,27)
         e2(2,j,3,2)=+rdat%r03(k,28)
         e2(3,j,3,2)=+rdat%r03(k,29)
         e2(4,j,3,2)=+rdat%r03(k,30)
      enddo

         j=1
         k= ind(j,3,2)
         m= in6(k)
      do i=1,5
         e1(i  ,3,2)=+rdat%r02(m,i+36)
      enddo
      do j=1,15
         k= ind(j,4,2)
         e5(  j,4,2)=+rdat%r06( k,1)-rdat%r05( j,4)
      enddo

      do j=1,10
         k= ind(j,4,2)
         e4(1,j,4,2)=+rdat%r05(k, 5)-rdat%r04(j,12)
         e4(2,j,4,2)=+rdat%r05(k, 6)-rdat%r04(j,13)
      enddo

      do j=1,6
         k= ind(j,4,2)
         e3(1,j,4,2)=+rdat%r04(k,14)-rdat%r03(j,19)
         e3(2,j,4,2)=+rdat%r04(k,15)-rdat%r03(j,20)
         e3(3,j,4,2)=+rdat%r04(k,16)-rdat%r03(j,21)
         e3(4,j,4,2)=+rdat%r04(k,17)-rdat%r03(j,22)
      enddo

      do j=1,3
         k= ind(j,4,2)
         m= in6(j)
         e2(1,j,4,2)=+rdat%r03(k,27)-rdat%r02(m,29)
         e2(2,j,4,2)=+rdat%r03(k,28)-rdat%r02(m,30)
         e2(3,j,4,2)=+rdat%r03(k,29)-rdat%r02(m,31)
         e2(4,j,4,2)=+rdat%r03(k,30)-rdat%r02(m,32)
      enddo

         j=1
         k= ind(j,4,2)
         m= in6(k)
      do i=1,5
         e1(i  ,4,2)=+rdat%r02(m,i+36)-rdat%r01(j,i+30)
      enddo
      do j=1,15
         k= ind(j,1,3)
         e5(  j,1,3)=-rdat%r05( k,3)
      enddo

      do j=1,10
         k= ind(j,1,3)
         e4(1,j,1,3)=-rdat%r04(k,10)
         e4(2,j,1,3)=-rdat%r04(k,11)
      enddo

      do j=1,6
         k= ind(j,1,3)
         e3(1,j,1,3)=-rdat%r03(k,15)
         e3(2,j,1,3)=-rdat%r03(k,16)
         e3(3,j,1,3)=-rdat%r03(k,17)
         e3(4,j,1,3)=-rdat%r03(k,18)
      enddo

      do j=1,3
         k= ind(j,1,3)
         m= in6(k)
         e2(1,j,1,3)=-rdat%r02(m,25)
         e2(2,j,1,3)=-rdat%r02(m,26)
         e2(3,j,1,3)=-rdat%r02(m,27)
         e2(4,j,1,3)=-rdat%r02(m,28)
      enddo

         j=1
         k= ind(j,1,3)
      do i=1,5
         e1(i  ,1,3)=-rdat%r01(k,i+25)
      enddo
      do j=1,15
         k= ind(j,3,3)
         e5(  j,3,3)=+rdat%r06( k,1)+rdat%r04( j,5)
      enddo

      do j=1,10
         k= ind(j,3,3)
         e4(1,j,3,3)=+rdat%r05(k, 5)+rdat%r03(j, 9)
         e4(2,j,3,3)=+rdat%r05(k, 6)+rdat%r03(j,10)
      enddo

      do j=1,6
         k= ind(j,3,3)
         m= in6(j)
         e3(1,j,3,3)=+rdat%r04(k,14)+rdat%r02(m,17)
         e3(2,j,3,3)=+rdat%r04(k,15)+rdat%r02(m,18)
         e3(3,j,3,3)=+rdat%r04(k,16)+rdat%r02(m,19)
         e3(4,j,3,3)=+rdat%r04(k,17)+rdat%r02(m,20)
      enddo

      do j=1,3
         k= ind(j,3,3)
         e2(1,j,3,3)=+rdat%r03(k,27)+rdat%r01(j,17)
         e2(2,j,3,3)=+rdat%r03(k,28)+rdat%r01(j,18)
         e2(3,j,3,3)=+rdat%r03(k,29)+rdat%r01(j,19)
         e2(4,j,3,3)=+rdat%r03(k,30)+rdat%r01(j,20)
      enddo

         j=1
         k= ind(j,3,3)
         m= in6(k)
      do i=1,5
         e1(i  ,3,3)=+rdat%r02(m,i+36)+rdat%r00(i,5)
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,4,3)
            l= k-2-jj

      e5(  j,4,3)=+rdat%r06(k,1)-rdat%r05(l,4)

            if(jj > 4) cycle
      e4(1,j,4,3)=+rdat%r05(k, 5)-rdat%r04(l,12)
      e4(2,j,4,3)=+rdat%r05(k, 6)-rdat%r04(l,13)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,4,3)=+rdat%r04(k,i+13)-rdat%r03(l,i+18)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      e2(i,j,4,3)=+rdat%r03(k,i+26)-rdat%r02(m,i+28)
            enddo

            if(jj > 1) cycle
            m= in6(k)
            do i=1,5
      e1(i  ,4,3)=+rdat%r02(m,i+36)-rdat%r01(l,i+30)
            enddo
         enddo
      enddo
      do j=1,15
         k= ind(j,1,4)
         e5(  j,1,4)=-rdat%r05( k,3)+rdat%r04( j,3)
      enddo

      do j=1,10
         k= ind(j,1,4)
         e4(1,j,1,4)=-rdat%r04(k,10)+rdat%r03(j, 5)
         e4(2,j,1,4)=-rdat%r04(k,11)+rdat%r03(j, 6)
      enddo

      do j=1,6
         k= ind(j,1,4)
         m= in6(j)
         e3(1,j,1,4)=-rdat%r03(k,15)+rdat%r02(m, 9)
         e3(2,j,1,4)=-rdat%r03(k,16)+rdat%r02(m,10)
         e3(3,j,1,4)=-rdat%r03(k,17)+rdat%r02(m,11)
         e3(4,j,1,4)=-rdat%r03(k,18)+rdat%r02(m,12)
      enddo

      do j=1,3
         k= ind(j,1,4)
         m= in6(k)
         e2(1,j,1,4)=-rdat%r02(m,25)+rdat%r01(j, 9)
         e2(2,j,1,4)=-rdat%r02(m,26)+rdat%r01(j,10)
         e2(3,j,1,4)=-rdat%r02(m,27)+rdat%r01(j,11)
         e2(4,j,1,4)=-rdat%r02(m,28)+rdat%r01(j,12)
      enddo

         j=1
         k= ind(j,1,4)
      do i=1,5
         e1(i  ,1,4)=-rdat%r01(k,i+25)+rdat%r00(i,3)
      enddo
      do j=1,15
         k= ind(j,2,4)
         e5(  j,2,4)=+rdat%r06( k,1)-rdat%r05( j,1)
      enddo

      do j=1,10
         k= ind(j,2,4)
         e4(1,j,2,4)=+rdat%r05(k, 5)-rdat%r04(j, 6)
         e4(2,j,2,4)=+rdat%r05(k, 6)-rdat%r04(j, 7)
      enddo

      do j=1,6
         k= ind(j,2,4)
         e3(1,j,2,4)=+rdat%r04(k,14)-rdat%r03(j,23)
         e3(2,j,2,4)=+rdat%r04(k,15)-rdat%r03(j,24)
         e3(3,j,2,4)=+rdat%r04(k,16)-rdat%r03(j,25)
         e3(4,j,2,4)=+rdat%r04(k,17)-rdat%r03(j,26)
      enddo

      do j=1,3
         k= ind(j,2,4)
         m= in6(j)
         e2(1,j,2,4)=+rdat%r03(k,27)-rdat%r02(m,33)
         e2(2,j,2,4)=+rdat%r03(k,28)-rdat%r02(m,34)
         e2(3,j,2,4)=+rdat%r03(k,29)-rdat%r02(m,35)
         e2(4,j,2,4)=+rdat%r03(k,30)-rdat%r02(m,36)
      enddo

         j=1
         k= ind(j,2,4)
         m= in6(k)
      do i=1,5
         e1(i  ,2,4)=+rdat%r02(m,i+36)-rdat%r01(j,i+35)
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,3,4)
            l= k-2-jj

      e5(  j,3,4)=+rdat%r06(k,1)-rdat%r05(l,1)

            if(jj > 4) cycle
      e4(1,j,3,4)=+rdat%r05(k, 5)-rdat%r04(l, 6)
      e4(2,j,3,4)=+rdat%r05(k, 6)-rdat%r04(l, 7)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,3,4)=+rdat%r04(k,i+13)-rdat%r03(l,i+22)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      e2(i,j,3,4)=+rdat%r03(k,i+26)-rdat%r02(m,i+32)
            enddo

            if(jj > 1) cycle
            m= in6(k)
            do i=1,5
      e1(i  ,3,4)=+rdat%r02(m,i+36)-rdat%r01(l,i+35)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,4,4)
            l= k-2-jj

      e5(  j,4,4)=+rdat%r06(k,1)-rdat%r05(l,1)-rdat%r05(l,4)+rdat%r04(j,4)+rdat%r04(j,5)

            if(jj > 4) cycle
      e4(1,j,4,4)=+rdat%r05(k, 5)-rdat%r04(l, 6)-rdat%r04(l,12)+rdat%r03(j, 7)+rdat%r03(j, 9)
      e4(2,j,4,4)=+rdat%r05(k, 6)-rdat%r04(l, 7)-rdat%r04(l,13)+rdat%r03(j, 8)+rdat%r03(j,10)

            if(jj > 3) cycle
            m= in6(j)
            do i=1,4
      e3(i,j,4,4)=+rdat%r04(k,i+13)-rdat%r03(l,i+18)-rdat%r03(l,i+22)                     &
     &                       +rdat%r02(m,i+12)+rdat%r02(m,i+16)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      e2(i,j,4,4)=+rdat%r03(k,i+26)-rdat%r02(m,i+28)-rdat%r02(m,i+32)                     &
     &                       +rdat%r01(j,i+12)+rdat%r01(j,i+16)
            enddo

            if(jj > 1) cycle
            m= in6(k)
            do i=1,5
      e1(i  ,4,4)=+rdat%r02(m,i+36)-rdat%r01(l,i+30)-rdat%r01(l,i+35)                     &
     &                       +rdat%r00(i,4)+rdat%r00(i,5)
            enddo
         enddo
      enddo

      xx= qx*qx
      zz= qz*qz
      xz= qx*qz
      xxx= xx*qx
      xxz= xx*qz
      xzz= zz*qx
      zzz= zz*qz
      xxxx= xx*xx
      xxxz= xx*xz
      xxzz= xx*zz
      xzzz= xz*zz
      zzzz= zz*zz

      qxd= qx+qx
      qzd= qz+qz
      xzd= xz+xz
      xxxd= xxx+xxx
      xxzd= xxz+xxz
      xzzd= xzz+xzz
      zzzd= zzz+zzz

      xzq= xz+xz+xz+xz

      do l=2,lx
         do k=2,kx

           if(k == 2 .and. l == 3) cycle

            f(1,1,k-1,l-1) = e5(  1,k,l)+e3(1,1,k,l)* 6 +e1(  1,k,l)* 3 +(+e4(1,1,k,l)+e4(2,1,k,l)+ (+e2(1,1,k,l)+e2(2,1,k,l))* 3 )*qxd +(+e3(2,1,k,l)+e3(3,1,k,l)* 4 +e3(4,1,k,l) +e1(  2,k,l)+e1(  3,k,l)* 4 +e1(  4,k,l))*xx +(+e2(3,1,k,l)+e2(4,1,k,l))*xxxd +e1(  5,k,l)*xxxx
            f(2,1,k-1,l-1) = e5(  4,k,l)+e3(1,4,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(2,4,k,l)+e2(2,1,k,l))*qxd +(+e3(4,4,k,l)+e1(  4,k,l))*xx
            f(3,1,k-1,l-1) = e5(  6,k,l)+e3(1,6,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(2,6,k,l)+e2(2,1,k,l))*qxd +(+e4(1,3,k,l)+e2(1,3,k,l))*qzd +(+e3(4,6,k,l)+e1(  4,k,l))*xx +e3(3,3,k,l)*xzq +(+e3(2,1,k,l)+e1(  2,k,l))*zz +e2(4,3,k,l)*xxzd +e2(3,1,k,l)*xzzd +e1(  5,k,l)*xxzz
            f(4,1,k-1,l-1) = e5(  2,k,l)+e3(1,2,k,l)* 3 +(+e4(1,2,k,l)+e4(2,2,k,l)* 2 +e2(1,2,k,l)+e2(2,2,k,l)* 2 )*qx +(+e3(3,2,k,l)* 2 +e3(4,2,k,l))*xx +e2(4,2,k,l)*xxx
            f(5,1,k-1,l-1) = e5(  3,k,l)+e3(1,3,k,l)* 3 +(+e4(1,3,k,l)+e4(2,3,k,l)* 2 +e2(1,3,k,l)+e2(2,3,k,l)* 2 )*qx +(+e4(1,1,k,l)+e2(1,1,k,l)* 3 )*qz +(+e3(3,3,k,l)* 2 +e3(4,3,k,l))*xx +(+e3(2,1,k,l)+e3(3,1,k,l)* 2 +e1(  2,k,l)+e1(  3,k,l)* 2 )*xz +e2(4,3,k,l)*xxx +(+e2(3,1,k,l)* 2 +e2(4,1,k,l))*xxz +e1(  5,k,l)*xxxz
            f(6,1,k-1,l-1) = e5(  5,k,l)+e3(1,5,k,l) +e4(2,5,k,l)*qxd +(+e4(1,2,k,l)+e2(1,2,k,l))*qz +e3(4,5,k,l)*xx +e3(3,2,k,l)*xzd +e2(4,2,k,l)*xxz

            f(1,2,k-1,l-1) = e5(  4,k,l)+e3(1,4,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(1,4,k,l)+e2(1,1,k,l))*qxd +(+e3(2,4,k,l)+e1(  2,k,l))*xx
            f(2,2,k-1,l-1) = e5( 11,k,l)+e3(1,4,k,l)* 6 +e1(  1,k,l)* 3
            f(3,2,k-1,l-1) = e5( 13,k,l)+e3(1,6,k,l)+e3(1,4,k,l)+e1(  1,k,l) +(+e4(1,8,k,l)+e2(1,3,k,l))*qzd +(+e3(2,4,k,l)+e1(  2,k,l))*zz
            f(4,2,k-1,l-1) = e5(  7,k,l)+e3(1,2,k,l)* 3 +(+e4(1,7,k,l)+e2(1,2,k,l)* 3 )*qx
            f(5,2,k-1,l-1) = e5(  8,k,l)+e3(1,3,k,l) +(+e4(1,8,k,l)+e2(1,3,k,l))*qx +(+e4(1,4,k,l)+e2(1,1,k,l))*qz +(+e3(2,4,k,l)+e1(  2,k,l))*xz
            f(6,2,k-1,l-1) = e5( 12,k,l)+e3(1,5,k,l)* 3 +(+e4(1,7,k,l)+e2(1,2,k,l)* 3 )*qz

            f(1,3,k-1,l-1) = e5(  6,k,l)+e3(1,6,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(1,6,k,l)+e2(1,1,k,l))*qxd +(+e4(2,3,k,l)+e2(2,3,k,l))*qzd +(+e3(2,6,k,l)+e1(  2,k,l))*xx +e3(3,3,k,l)*xzq +(+e3(4,1,k,l)+e1(  4,k,l))*zz +e2(3,3,k,l)*xxzd +e2(4,1,k,l)*xzzd +e1(  5,k,l)*xxzz
            f(2,3,k-1,l-1) = e5( 13,k,l)+e3(1,6,k,l)+e3(1,4,k,l)+e1(  1,k,l) +(+e4(2,8,k,l)+e2(2,3,k,l))*qzd +(+e3(4,4,k,l)+e1(  4,k,l))*zz
            f(3,3,k-1,l-1) = e5( 15,k,l)+e3(1,6,k,l)* 6 +e1(  1,k,l)* 3 +(+e4(1,10,k,l)+e4(2,10,k,l)+ (+e2(1,3,k,l)+e2(2,3,k,l))* 3 )*qzd +(+e3(2,6,k,l)+e3(3,6,k,l)* 4 +e3(4,6,k,l) +e1(  2,k,l)+e1(  3,k,l)* 4 +e1(  4,k,l))*zz +(+e2(3,3,k,l)+e2(4,3,k,l))*zzzd +e1(  5,k,l)*zzzz
            f(4,3,k-1,l-1) = e5(  9,k,l)+e3(1,2,k,l) +(+e4(1,9,k,l)+e2(1,2,k,l))*qx +e4(2,5,k,l)*qzd +e3(3,5,k,l)*xzd +e3(4,2,k,l)*zz +e2(4,2,k,l)*xzz
            f(5,3,k-1,l-1) = e5( 10,k,l)+e3(1,3,k,l)* 3 +(+e4(1,10,k,l)+e2(1,3,k,l)* 3 )*qx +(+e4(1,6,k,l)+e4(2,6,k,l)* 2 +e2(1,1,k,l)+e2(2,1,k,l)* 2 )*qz +(+e3(2,6,k,l)+e3(3,6,k,l)* 2 +e1(  2,k,l)+e1(  3,k,l)* 2 )*xz +(+e3(3,3,k,l)* 2 +e3(4,3,k,l))*zz +(+e2(3,3,k,l)* 2 +e2(4,3,k,l))*xzz +e2(4,1,k,l)*zzz +e1(  5,k,l)*xzzz
            f(6,3,k-1,l-1) = e5( 14,k,l)+e3(1,5,k,l)* 3 +(+e4(1,9,k,l)+e4(2,9,k,l)* 2 +e2(1,2,k,l)+e2(2,2,k,l)* 2 )*qz +(+e3(3,5,k,l)* 2 +e3(4,5,k,l))*zz +e2(4,2,k,l)*zzz

            f(1,4,k-1,l-1) = e5(  2,k,l)+e3(1,2,k,l)* 3 +(+e4(1,2,k,l)* 2 +e4(2,2,k,l) +e2(1,2,k,l)* 2 +e2(2,2,k,l))*qx +(+e3(2,2,k,l)+e3(3,2,k,l)* 2 )*xx +e2(3,2,k,l)*xxx
            f(2,4,k-1,l-1) = e5(  7,k,l)+e3(1,2,k,l)* 3 +(+e4(2,7,k,l)+e2(2,2,k,l)* 3 )*qx
            f(3,4,k-1,l-1) = e5(  9,k,l)+e3(1,2,k,l) +(+e4(2,9,k,l)+e2(2,2,k,l))*qx +e4(1,5,k,l)*qzd +e3(3,5,k,l)*xzd +e3(2,2,k,l)*zz +e2(3,2,k,l)*xzz
            f(4,4,k-1,l-1) = e5(  4,k,l)+e3(1,4,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(1,4,k,l)+e4(2,4,k,l) +e2(1,1,k,l)+e2(2,1,k,l))*qx +(+e3(3,4,k,l)+e1(  3,k,l))*xx
            f(5,4,k-1,l-1) = e5(  5,k,l)+e3(1,5,k,l) +(+e4(1,5,k,l)+e4(2,5,k,l))*qx +(+e4(1,2,k,l)+e2(1,2,k,l))*qz +e3(3,5,k,l)*xx +(+e3(2,2,k,l)+e3(3,2,k,l))*xz +e2(3,2,k,l)*xxz
            f(6,4,k-1,l-1) = e5(  8,k,l)+e3(1,3,k,l) +(+e4(2,8,k,l)+e2(2,3,k,l))*qx +(+e4(1,4,k,l)+e2(1,1,k,l))*qz +(+e3(3,4,k,l)+e1(  3,k,l))*xz

            f(1,5,k-1,l-1) = e5(  3,k,l)+e3(1,3,k,l)* 3 +(+e4(1,3,k,l)* 2 +e4(2,3,k,l) +e2(1,3,k,l)* 2 +e2(2,3,k,l))*qx +(+e4(2,1,k,l)+e2(2,1,k,l)* 3 )*qz +(+e3(2,3,k,l)+e3(3,3,k,l)* 2 )*xx +(+e3(3,1,k,l)* 2 +e3(4,1,k,l) +e1(  3,k,l)* 2 +e1(  4,k,l))*xz +e2(3,3,k,l)*xxx +(+e2(3,1,k,l)+e2(4,1,k,l)* 2 )*xxz +e1(  5,k,l)*xxxz
            f(2,5,k-1,l-1) = e5(  8,k,l)+e3(1,3,k,l) +(+e4(2,8,k,l)+e2(2,3,k,l))*qx +(+e4(2,4,k,l)+e2(2,1,k,l))*qz +(+e3(4,4,k,l)+e1(  4,k,l))*xz
            f(3,5,k-1,l-1) = e5( 10,k,l)+e3(1,3,k,l)* 3 +(+e4(2,10,k,l)+e2(2,3,k,l)* 3 )*qx +(+e4(1,6,k,l)* 2 +e4(2,6,k,l) +e2(1,1,k,l)* 2 +e2(2,1,k,l))*qz +(+e3(3,6,k,l)* 2 +e3(4,6,k,l) +e1(  3,k,l)* 2 +e1(  4,k,l))*xz +(+e3(2,3,k,l)+e3(3,3,k,l)* 2 )*zz +(+e2(3,3,k,l)+e2(4,3,k,l)* 2 )*xzz +e2(3,1,k,l)*zzz +e1(  5,k,l)*xzzz
            f(4,5,k-1,l-1) = e5(  5,k,l)+e3(1,5,k,l) +(+e4(1,5,k,l)+e4(2,5,k,l))*qx +(+e4(2,2,k,l)+e2(2,2,k,l))*qz +e3(3,5,k,l)*xx +(+e3(3,2,k,l)+e3(4,2,k,l))*xz +e2(4,2,k,l)*xxz
            f(5,5,k-1,l-1) = e5(  6,k,l)+e3(1,6,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(1,6,k,l)+e4(2,6,k,l) +e2(1,1,k,l)+e2(2,1,k,l))*qx +(+e4(1,3,k,l)+e4(2,3,k,l) +e2(1,3,k,l)+e2(2,3,k,l))*qz +(+e3(3,6,k,l)+e1(  3,k,l))*xx +(+e3(2,3,k,l)+e3(3,3,k,l)* 2 +e3(4,3,k,l))*xz +(+e3(3,1,k,l)+e1(  3,k,l))*zz +(+e2(3,3,k,l)+e2(4,3,k,l))*xxz +(+e2(3,1,k,l)+e2(4,1,k,l))*xzz +e1(  5,k,l)*xxzz
            f(6,5,k-1,l-1) = e5(  9,k,l)+e3(1,2,k,l) +(+e4(2,9,k,l)+e2(2,2,k,l))*qx +(+e4(1,5,k,l)+e4(2,5,k,l))*qz +(+e3(3,5,k,l)+e3(4,5,k,l))*xz +e3(3,2,k,l)*zz +e2(4,2,k,l)*xzz

            f(1,6,k-1,l-1) = e5(  5,k,l)+e3(1,5,k,l) +e4(1,5,k,l)*qxd +(+e4(2,2,k,l)+e2(2,2,k,l))*qz +e3(2,5,k,l)*xx +e3(3,2,k,l)*xzd +e2(3,2,k,l)*xxz
            f(2,6,k-1,l-1) = e5( 12,k,l)+e3(1,5,k,l)* 3 +(+e4(2,7,k,l)+e2(2,2,k,l)* 3 )*qz
            f(3,6,k-1,l-1) = e5( 14,k,l)+e3(1,5,k,l)* 3 +(+e4(1,9,k,l)* 2 +e4(2,9,k,l) +e2(1,2,k,l)* 2 +e2(2,2,k,l))*qz +(+e3(2,5,k,l)+e3(3,5,k,l)* 2 )*zz +e2(3,2,k,l)*zzz
            f(4,6,k-1,l-1) = e5(  8,k,l)+e3(1,3,k,l) +(+e4(1,8,k,l)+e2(1,3,k,l))*qx +(+e4(2,4,k,l)+e2(2,1,k,l))*qz +(+e3(3,4,k,l)+e1(  3,k,l))*xz
            f(5,6,k-1,l-1) = e5(  9,k,l)+e3(1,2,k,l) +(+e4(1,9,k,l)+e2(1,2,k,l))*qx +(+e4(1,5,k,l)+e4(2,5,k,l))*qz +(+e3(2,5,k,l)+e3(3,5,k,l))*xz +e3(3,2,k,l)*zz +e2(3,2,k,l)*xzz
            f(6,6,k-1,l-1) = e5( 13,k,l)+e3(1,6,k,l)+e3(1,4,k,l)+e1(  1,k,l) +(+e4(1,8,k,l)+e4(2,8,k,l) +e2(1,3,k,l)+e2(2,3,k,l))*qz +(+e3(3,4,k,l)+e1(  3,k,l))*zz

         enddo
      enddo

      f(:,:,1,2)= f(:,:,2,1)

      end subroutine mcdv_18

! >
! >    @brief   dpdp case
! >
! >    @details integration of a dpdp case
! >
      subroutine mcdv_19(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,3,6,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 6, lx = 4
      real(kind=dp) :: d1(  5,kx,lx),d2(4,3,kx,lx),d3(3,6,kx,lx), &
                       d4( 10,kx,lx)

      integer, parameter :: ind (10, kx, lx) = reshape(&
        [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 4, 7, 8, 11, 12, 13, 16, 17, 18, &
          19, 6, 9, 10, 13, 14, 15, 18, 19, 20, 21, 2, 4, 5, 7, 8, 9, 11, &
          12, 13, 14, 3, 5, 6, 8, 9, 10, 12, 13, 14, 15, 5, 8, 9, 12, 13, &
          14, 17, 18, 19, 20, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 4, 7, 8, 11, &
          12, 13, 16, 17, 18, 19, 6, 9, 10, 13, 14, 15, 18, 19, 20, 21, 2, &
          4, 5, 7, 8, 9, 11, 12, 13, 14, 3, 5, 6, 8, 9, 10, 12, 13, 14, &
          15, 5, 8, 9, 12, 13, 14, 17, 18, 19, 20, 2, 4, 5, 7, 8, 9, 11, &
          12, 13, 14, 7, 11, 12, 16, 17, 18, 22, 23, 24, 25, 9, 13, 14, &
          18, 19, 20, 24, 25, 26, 27, 4, 7, 8, 11, 12, 13, 16, 17, 18, 19, &
          5, 8, 9, 12, 13, 14, 17, 18, 19, 20, 8, 12, 13, 17, 18, 19, 23, &
          24, 25, 26, 3, 5, 6, 8, 9, 10, 12, 13, 14, 15, 8, 12, 13, 17, &
          18, 19, 23, 24, 25, 26, 10, 14, 15, 19, 20, 21, 25, 26, 27, 28, &
          5, 8, 9, 12, 13, 14, 17, 18, 19, 20, 6, 9, 10, 13, 14, 15, 18, &
          19, 20, 21, 9, 13, 14, 18, 19, 20, 24, 25, 26, 27] &
        , shape(ind))
      integer :: i, j, k, l, m, n, ii, jj
      real(kind=dp) :: xx, zz, xz, xxx, xxz, xzz, zzz
      real(kind=dp) :: qxd, qzd, xzd


      do j=1,10
         d4(  j,1,1)=+rdat%r05( j,2)+rdat%r03( j,1)
      enddo

      do j=1,6
         m= in6(j)
         d3(1,j,1,1)=+rdat%r04(j, 5)+rdat%r02(m, 1)
         d3(2,j,1,1)=+rdat%r04(j, 6)+rdat%r02(m, 2)
         d3(3,j,1,1)=+rdat%r04(j, 7)+rdat%r02(m, 3)
      enddo

      do j=1,3
         d2(1,j,1,1)=+rdat%r03(j,18)+rdat%r01(j, 1)
         d2(2,j,1,1)=+rdat%r03(j,19)+rdat%r01(j, 2)
         d2(3,j,1,1)=+rdat%r03(j,20)+rdat%r01(j, 3)
         d2(4,j,1,1)=+rdat%r03(j,21)+rdat%r01(j, 4)
      enddo

         j=1
      do i=1,5
         d1(i  ,1,1)=+rdat%r02(j,i+31)+rdat%r00(i,1)
      enddo
      do j=1,10
         k= ind(j,2,1)
         d4(  j,2,1)=+rdat%r05( k,2)+rdat%r03( j,1)
      enddo

      do j=1,6
         k= ind(j,2,1)
         m= in6(j)
         d3(1,j,2,1)=+rdat%r04(k, 5)+rdat%r02(m, 1)
         d3(2,j,2,1)=+rdat%r04(k, 6)+rdat%r02(m, 2)
         d3(3,j,2,1)=+rdat%r04(k, 7)+rdat%r02(m, 3)
      enddo

      do j=1,3
         k= ind(j,2,1)
         d2(1,j,2,1)=+rdat%r03(k,18)+rdat%r01(j, 1)
         d2(2,j,2,1)=+rdat%r03(k,19)+rdat%r01(j, 2)
         d2(3,j,2,1)=+rdat%r03(k,20)+rdat%r01(j, 3)
         d2(4,j,2,1)=+rdat%r03(k,21)+rdat%r01(j, 4)
      enddo

         j=1
         k= ind(j,2,1)
         m= in6(k)
      do i=1,5
         d1(i  ,2,1)=+rdat%r02(m,i+31)+rdat%r00(i,1)
      enddo
      j= 0
      do jj=1,4
         do ii=1,jj
            j= j+1
            k= ind(j,3,1)
            l= k-2-jj

      d4(  j,3,1)=+rdat%r05(k,2)-rdat%r04(l,2)* 2 +rdat%r03(j,1)+rdat%r03(j,2)

            if(jj > 3) cycle
            m= in6(j)
            do i=1,3
      d3(i,j,3,1)=+rdat%r04(k,i+ 4)-rdat%r03(l,i+ 8)* 2 +rdat%r02(m,i   )+rdat%r02(m,i+ 3)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      d2(i,j,3,1)=+rdat%r03(k,i+17)-rdat%r02(m,i+15)* 2 +rdat%r01(j,i   )+rdat%r01(j,i+ 4)
            enddo

            if(jj > 1) cycle
            m= in6(k)
            do i=1,5
      d1(i  ,3,1)=+rdat%r02(m,i+31)-rdat%r01(l,i+20)* 2 +rdat%r00(i,1)+rdat%r00(i,2)
            enddo
         enddo
      enddo
      do j=1,10
         k= ind(j,4,1)
         d4(  j,4,1)=+rdat%r05( k,2)
      enddo

      do j=1,6
         k= ind(j,4,1)
         d3(1,j,4,1)=+rdat%r04(k, 5)
         d3(2,j,4,1)=+rdat%r04(k, 6)
         d3(3,j,4,1)=+rdat%r04(k, 7)
      enddo

      do j=1,3
         k= ind(j,4,1)
         d2(1,j,4,1)=+rdat%r03(k,18)
         d2(2,j,4,1)=+rdat%r03(k,19)
         d2(3,j,4,1)=+rdat%r03(k,20)
         d2(4,j,4,1)=+rdat%r03(k,21)
      enddo

         j=1
         k= ind(j,4,1)
         m= in6(k)
      do i=1,5
         d1(i  ,4,1)=+rdat%r02(m,i+31)
      enddo
      do j=1,10
         k= ind(j,5,1)
         d4(  j,5,1)=+rdat%r05( k,2)-rdat%r04( j,2)
      enddo

      do j=1,6
         k= ind(j,5,1)
         d3(1,j,5,1)=+rdat%r04(k, 5)-rdat%r03(j, 9)
         d3(2,j,5,1)=+rdat%r04(k, 6)-rdat%r03(j,10)
         d3(3,j,5,1)=+rdat%r04(k, 7)-rdat%r03(j,11)
      enddo

      do j=1,3
         k= ind(j,5,1)
         m= in6(j)
         d2(1,j,5,1)=+rdat%r03(k,18)-rdat%r02(m,16)
         d2(2,j,5,1)=+rdat%r03(k,19)-rdat%r02(m,17)
         d2(3,j,5,1)=+rdat%r03(k,20)-rdat%r02(m,18)
         d2(4,j,5,1)=+rdat%r03(k,21)-rdat%r02(m,19)
      enddo

         j=1
         k= ind(j,5,1)
         m= in6(k)
      do i=1,5
         d1(i  ,5,1)=+rdat%r02(m,i+31)-rdat%r01(j,i+20)
      enddo
      j= 0
      do jj=1,4
         do ii=1,jj
            j= j+1
            k= ind(j,6,1)
            l= k-2-jj

      d4(  j,6,1)=+rdat%r05(k,2)-rdat%r04(l,2)

            if(jj > 3) cycle
            do i=1,3
      d3(i,j,6,1)=+rdat%r04(k,i+ 4)-rdat%r03(l,i+ 8)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      d2(i,j,6,1)=+rdat%r03(k,i+17)-rdat%r02(m,i+15)
            enddo

            if(jj > 1) cycle
            m= in6(k)
            do i=1,5
      d1(i  ,6,1)=+rdat%r02(m,i+31)-rdat%r01(l,i+20)
            enddo
         enddo
      enddo
      do j=1,10
         d4(  j,1,2)=-rdat%r06( j,1)-rdat%r04( j,4)* 3
      enddo

      do j=1,6
         d3(1,j,1,2)=-rdat%r05(j, 4)-rdat%r03(j,15)* 3
         d3(2,j,1,2)=-rdat%r05(j, 5)-rdat%r03(j,16)* 3
         d3(3,j,1,2)=-rdat%r05(j, 6)-rdat%r03(j,17)* 3
      enddo

      do j=1,3
         m= in6(j)
         d2(1,j,1,2)=-rdat%r04(j,14)-rdat%r02(m,24)* 3
         d2(2,j,1,2)=-rdat%r04(j,15)-rdat%r02(m,25)* 3
         d2(3,j,1,2)=-rdat%r04(j,16)-rdat%r02(m,26)* 3
         d2(4,j,1,2)=-rdat%r04(j,17)-rdat%r02(m,27)* 3
      enddo

         j=1
      do i=1,5
         d1(i  ,1,2)=-rdat%r03(j,i+29)-rdat%r01(j,i+30)* 3
      enddo
      do j=1,10
         k= ind(j,2,2)
         d4(  j,2,2)=-rdat%r06( k,1)-rdat%r04( j,4)
      enddo

      do j=1,6
         k= ind(j,2,2)
         d3(1,j,2,2)=-rdat%r05(k, 4)-rdat%r03(j,15)
         d3(2,j,2,2)=-rdat%r05(k, 5)-rdat%r03(j,16)
         d3(3,j,2,2)=-rdat%r05(k, 6)-rdat%r03(j,17)
      enddo

      do j=1,3
         k= ind(j,2,2)
         m= in6(j)
         d2(1,j,2,2)=-rdat%r04(k,14)-rdat%r02(m,24)
         d2(2,j,2,2)=-rdat%r04(k,15)-rdat%r02(m,25)
         d2(3,j,2,2)=-rdat%r04(k,16)-rdat%r02(m,26)
         d2(4,j,2,2)=-rdat%r04(k,17)-rdat%r02(m,27)
      enddo

         j=1
         k= ind(j,2,2)
      do i=1,5
         d1(i  ,2,2)=-rdat%r03(k,i+29)-rdat%r01(j,i+30)
      enddo
      j= 0
      do jj=1,4
         do ii=1,jj
            j= j+1
            k= ind(j,3,2)
            l= k-2-jj

      d4(  j,3,2)=-rdat%r06(k,1)+rdat%r05(l,3)* 2 -rdat%r04(j,3)-rdat%r04(j,4)

            if(jj > 3) cycle
            do i=1,3
      d3(i,j,3,2)=-rdat%r05(k,i+ 3)+rdat%r04(l,i+ 7)* 2 -rdat%r03(j,i+11)-rdat%r03(j,i+14)
            enddo

            if(jj > 2) cycle
            m= in6(j)
            do i=1,4
      d2(i,j,3,2)=-rdat%r04(k,i+13)+rdat%r03(l,i+21)* 2 -rdat%r02(m,i+19)-rdat%r02(m,i+23)
            enddo

            if(jj > 1) cycle
            m= in6(l)
            do i=1,5
      d1(i  ,3,2)=-rdat%r03(k,i+29)+rdat%r02(m,i+36)* 2 -rdat%r01(j,i+25)-rdat%r01(j,i+30)
            enddo
         enddo
      enddo
      do j=1,10
         k= ind(j,4,2)
         d4(  j,4,2)=-rdat%r06( k,1)-rdat%r04( k,4)
      enddo

      do j=1,6
         k= ind(j,4,2)
         d3(1,j,4,2)=-rdat%r05(k, 4)-rdat%r03(k,15)
         d3(2,j,4,2)=-rdat%r05(k, 5)-rdat%r03(k,16)
         d3(3,j,4,2)=-rdat%r05(k, 6)-rdat%r03(k,17)
      enddo

      do j=1,3
         k= ind(j,4,2)
         m= in6(k)
         d2(1,j,4,2)=-rdat%r04(k,14)-rdat%r02(m,24)
         d2(2,j,4,2)=-rdat%r04(k,15)-rdat%r02(m,25)
         d2(3,j,4,2)=-rdat%r04(k,16)-rdat%r02(m,26)
         d2(4,j,4,2)=-rdat%r04(k,17)-rdat%r02(m,27)
      enddo

         j=1
         k= ind(j,4,2)
      do i=1,5
         d1(i  ,4,2)=-rdat%r03(k,i+29)-rdat%r01(k,i+30)
      enddo
      do j=1,10
         k= ind(j,5,2)
         d4(  j,5,2)=-rdat%r06( k,1)-rdat%r04( k,4)+rdat%r05( j,3)+rdat%r03( j,5)
      enddo

      do j=1,6
         k= ind(j,5,2)
         m= in6(j)
         d3(1,j,5,2)=-rdat%r05(k, 4)-rdat%r03(k,15)+rdat%r04(j, 8)+rdat%r02(m,13)
         d3(2,j,5,2)=-rdat%r05(k, 5)-rdat%r03(k,16)+rdat%r04(j, 9)+rdat%r02(m,14)
         d3(3,j,5,2)=-rdat%r05(k, 6)-rdat%r03(k,17)+rdat%r04(j,10)+rdat%r02(m,15)
      enddo

      do j=1,3
         k= ind(j,5,2)
         m= in6(k)
         d2(1,j,5,2)=-rdat%r04(k,14)-rdat%r02(m,24)+rdat%r03(j,22)+rdat%r01(j,17)
         d2(2,j,5,2)=-rdat%r04(k,15)-rdat%r02(m,25)+rdat%r03(j,23)+rdat%r01(j,18)
         d2(3,j,5,2)=-rdat%r04(k,16)-rdat%r02(m,26)+rdat%r03(j,24)+rdat%r01(j,19)
         d2(4,j,5,2)=-rdat%r04(k,17)-rdat%r02(m,27)+rdat%r03(j,25)+rdat%r01(j,20)
      enddo

         j=1
         k= ind(j,5,2)
      do i=1,5
         d1(i  ,5,2)=-rdat%r03(k,i+29)-rdat%r01(k,i+30)+rdat%r02(j,i+36)+rdat%r00(i,5)
      enddo
      j= 0
      do jj=1,4
         do ii=1,jj
            j= j+1
            k= ind(j,6,2)
            l= k-2-jj

      d4(  j,6,2)=-rdat%r06(k,1)+rdat%r05(l,3)

            if(jj > 3) cycle
            do i=1,3
      d3(i,j,6,2)=-rdat%r05(k,i+ 3)+rdat%r04(l,i+ 7)
            enddo

            if(jj > 2) cycle
            do i=1,4
      d2(i,j,6,2)=-rdat%r04(k,i+13)+rdat%r03(l,i+21)
            enddo

            if(jj > 1) cycle
            m= in6(l)
            do i=1,5
      d1(i  ,6,2)=-rdat%r03(k,i+29)+rdat%r02(m,i+36)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,4
         do ii=1,jj
            j= j+1
            k= ind(j,2,3)
            l= k-3-jj-jj

      d4(  j,2,3)=-rdat%r06(k,1)-rdat%r04(l,4)* 3

            if(jj > 3) cycle
            do i=1,3
      d3(i,j,2,3)=-rdat%r05(k,i+ 3)-rdat%r03(l,i+14)* 3
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      d2(i,j,2,3)=-rdat%r04(k,i+13)-rdat%r02(m,i+23)* 3
            enddo

            if(jj > 1) cycle
            do i=1,5
      d1(i  ,2,3)=-rdat%r03(k,i+29)-rdat%r01(l,i+30)* 3
            enddo
         enddo
      enddo
      j= 0
      do jj=1,4
         do ii=1,jj
            j= j+1
            k= ind(j,3,3)
            l= k-3-jj
            m= l-2-jj

      d4(  j,3,3)=-rdat%r06(k,1)+rdat%r05(l,3)* 2 -rdat%r04(m,3)-rdat%r04(m,4)

            if(jj > 3) cycle
            do i=1,3
      d3(i,j,3,3)=-rdat%r05(k,i+ 3)+rdat%r04(l,i+ 7)* 2 -rdat%r03(m,i+11)-rdat%r03(m,i+14)
            enddo

            if(jj > 2) cycle
            n= in6(m)
            do i=1,4
      d2(i,j,3,3)=-rdat%r04(k,i+13)+rdat%r03(l,i+21)* 2 -rdat%r02(n,i+19)-rdat%r02(n,i+23)
            enddo

            if(jj > 1) cycle
            n= in6(l)
            do i=1,5
      d1(i  ,3,3)=-rdat%r03(k,i+29)+rdat%r02(n,i+36)* 2 -rdat%r01(m,i+25)-rdat%r01(m,i+30)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,4
         do ii=1,jj
            j= j+1
            k= ind(j,6,3)
            l= k-3-jj
            m= l-jj

      d4(  j,6,3)=-rdat%r06(k,1)-rdat%r04(m,4)+rdat%r05(l,3)+rdat%r03(j,5)

            if(jj > 3) cycle
            n= in6(j)
            do i=1,3
      d3(i,j,6,3)=-rdat%r05(k,i+ 3)-rdat%r03(m,i+14)+rdat%r04(l,i+ 7)+rdat%r02(n,i+12)
            enddo

            if(jj > 2) cycle
            n= in6(m)
            do i=1,4
      d2(i,j,6,3)=-rdat%r04(k,i+13)-rdat%r02(n,i+23)+rdat%r03(l,i+21)+rdat%r01(j,i+16)
            enddo

            if(jj > 1) cycle
            n= in6(l)
            do i=1,5
      d1(i  ,6,3)=-rdat%r03(k,i+29)-rdat%r01(m,i+30)+rdat%r02(n,i+36)+rdat%r00(i,5)
            enddo
         enddo
      enddo
      do j=1,10
         k= ind(j,1,4)
         d4(  j,1,4)=-rdat%r06( k,1)-rdat%r04( k,4)+rdat%r05( j,1)+rdat%r03( j,4)
      enddo

      do j=1,6
         k= ind(j,1,4)
         m= in6(j)
         d3(1,j,1,4)=-rdat%r05(k, 4)-rdat%r03(k,15)+rdat%r04(j,11)+rdat%r02(m,10)
         d3(2,j,1,4)=-rdat%r05(k, 5)-rdat%r03(k,16)+rdat%r04(j,12)+rdat%r02(m,11)
         d3(3,j,1,4)=-rdat%r05(k, 6)-rdat%r03(k,17)+rdat%r04(j,13)+rdat%r02(m,12)
      enddo

      do j=1,3
         k= ind(j,1,4)
         m= in6(k)
         d2(1,j,1,4)=-rdat%r04(k,14)-rdat%r02(m,24)+rdat%r03(j,26)+rdat%r01(j,13)
         d2(2,j,1,4)=-rdat%r04(k,15)-rdat%r02(m,25)+rdat%r03(j,27)+rdat%r01(j,14)
         d2(3,j,1,4)=-rdat%r04(k,16)-rdat%r02(m,26)+rdat%r03(j,28)+rdat%r01(j,15)
         d2(4,j,1,4)=-rdat%r04(k,17)-rdat%r02(m,27)+rdat%r03(j,29)+rdat%r01(j,16)
      enddo

         j=1
         k= ind(j,1,4)
      do i=1,5
         d1(i  ,1,4)=-rdat%r03(k,i+29)-rdat%r01(k,i+30)+rdat%r02(j,i+41)+rdat%r00(i,4)
      enddo
      j= 0
      do jj=1,4
         do ii=1,jj
            j= j+1
            k= ind(j,2,4)
            l= k-3-jj
            m= l-jj

      d4(  j,2,4)=-rdat%r06(k,1)+rdat%r05(l,1)-rdat%r04(m,4)+rdat%r03(j,4)

            if(jj > 3) cycle
            n= in6(j)
            do i=1,3
      d3(i,j,2,4)=-rdat%r05(k,i+ 3)+rdat%r04(l,i+10)-rdat%r03(m,i+14)+rdat%r02(n,i+ 9)
            enddo

            if(jj > 2) cycle
            n= in6(m)
            do i=1,4
      d2(i,j,2,4)=-rdat%r04(k,i+13)+rdat%r03(l,i+25)-rdat%r02(n,i+23)+rdat%r01(j,i+12)
            enddo

            if(jj > 1) cycle
            n= in6(l)
            do i=1,5
      d1(i  ,2,4)=-rdat%r03(k,i+29)+rdat%r02(n,i+41)-rdat%r01(m,i+30)+rdat%r00(i,4)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,4
         do ii=1,jj
            j= j+1
            k= ind(j,3,4)
            l= k-3-jj
            m= l-2-jj

      d4(  j,3,4)=-rdat%r06(k,1)    +rdat%r05(l,1)+rdat%r05(l,3)* 2                       &
     &            -rdat%r04(m,1)* 2 -rdat%r04(m,3)-rdat%r04(m,4)* 3                       &
     &            +rdat%r03(j,3)    +rdat%r03(j,4)+rdat%r03(j,5)* 2

            if(jj > 3) cycle
            n= in6(j)
            do i=1,3
      d3(i,j,3,4)=-rdat%r05(k,i+ 3)    +rdat%r04(l,i+ 7)* 2 +rdat%r04(l,i+10)             &
     &            -rdat%r03(m,i+ 5)* 2 -rdat%r03(m,i+11)    -rdat%r03(m,i+14)* 3          &
     &            +rdat%r02(n,i+ 6)    +rdat%r02(n,i+ 9)    +rdat%r02(n,i+12)* 2
            enddo

            if(jj > 2) cycle
            n= in6(m)
            do i=1,4
      d2(i,j,3,4)=-rdat%r04(k,i+13)+rdat%r03(l,i+21)* 2 +rdat%r03(l,i+25)                 &
     &            -rdat%r02(n,i+19)-rdat%r02(n,i+23)* 3 -rdat%r02(n,i+27)* 2              &
     &            +rdat%r01(j,i+ 8)+rdat%r01(j,i+12)    +rdat%r01(j,i+16)* 2
            enddo

            if(jj > 1) cycle
            n= in6(l)
            do i=1,5
      d1(i  ,3,4)=-rdat%r03(k,i+29)+rdat%r02(n,i+36)* 2 +rdat%r02(n,i+41)                 &
     &            -rdat%r01(m,i+25)-rdat%r01(m,i+30)* 3 -rdat%r01(m,i+35)* 2              &
     &            +rdat%r00(i,3)+rdat%r00(i,4)    +rdat%r00(i,5)* 2
            enddo
         enddo
      enddo
      j= 0
      do jj=1,4
         do ii=1,jj
            j= j+1
            k= ind(j,4,4)
            l= k-2-jj

      d4(  j,4,4)=-rdat%r06(k,1)+rdat%r05(l,1)

            if(jj > 3) cycle
            do i=1,3
      d3(i,j,4,4)=-rdat%r05(k,i+ 3)+rdat%r04(l,i+10)
            enddo

            if(jj > 2) cycle
            do i=1,4
      d2(i,j,4,4)=-rdat%r04(k,i+13)+rdat%r03(l,i+25)
            enddo

            if(jj > 1) cycle
            m= in6(l)
            do i=1,5
      d1(i  ,4,4)=-rdat%r03(k,i+29)+rdat%r02(m,i+41)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,4
         do ii=1,jj
            j= j+1
            k= ind(j,5,4)
            l= k-2-jj

      d4(  j,5,4)=-rdat%r06(k,1)+rdat%r05(l,1)+rdat%r05(l,3)                              &
     &                    -rdat%r04(j,1)-rdat%r04(j,4)

            if(jj > 3) cycle
            do i=1,3
      d3(i,j,5,4)=-rdat%r05(k,i+ 3)+rdat%r04(l,i+ 7)+rdat%r04(l,i+10)                     &
     &                       -rdat%r03(j,i+ 5)-rdat%r03(j,i+14)
            enddo

            if(jj > 2) cycle
            m= in6(j)
            do i=1,4
      d2(i,j,5,4)=-rdat%r04(k,i+13)+rdat%r03(l,i+21)+rdat%r03(l,i+25)                     &
     &                       -rdat%r02(m,i+23)-rdat%r02(m,i+27)
            enddo

            if(jj > 1) cycle
            m= in6(l)
            do i=1,5
      d1(i  ,5,4)=-rdat%r03(k,i+29)+rdat%r02(m,i+36)+rdat%r02(m,i+41)                     &
     &                       -rdat%r01(j,i+30)-rdat%r01(j,i+35)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,4
         do ii=1,jj
            j= j+1
            k= ind(j,6,4)
            l= k-3-jj
            m= l-2-jj

      d4(  j,6,4)=-rdat%r06(k,1)+rdat%r05(l,1)+rdat%r05(l,3)                              &
     &                    -rdat%r04(m,1)-rdat%r04(m,4)

            if(jj > 3) cycle
            do i=1,3
      d3(i,j,6,4)=-rdat%r05(k,i+ 3)+rdat%r04(l,i+ 7)+rdat%r04(l,i+10)                     &
     &                       -rdat%r03(m,i+ 5)-rdat%r03(m,i+14)
            enddo

            if(jj > 2) cycle
            n= in6(m)
            do i=1,4
      d2(i,j,6,4)=-rdat%r04(k,i+13)+rdat%r03(l,i+21)+rdat%r03(l,i+25)                     &
     &                       -rdat%r02(n,i+23)-rdat%r02(n,i+27)
            enddo

            if(jj > 1) cycle
            n= in6(l)
            do i=1,5
      d1(i  ,6,4)=-rdat%r03(k,i+29)+rdat%r02(n,i+36)+rdat%r02(n,i+41)                     &
     &                       -rdat%r01(m,i+30)-rdat%r01(m,i+35)
            enddo
         enddo
      enddo

      xx= qx*qx
      zz= qz*qz
      xz= qx*qz
      xxx= xx*qx
      xxz= xx*qz
      xzz= zz*qx
      zzz= zz*qz

      qxd= qx+qx
      qzd= qz+qz

      xzd= xz+xz

      do l=2,lx
         do k=1,kx

           if(k == 1 .and. l == 3) cycle
           if(k == 4 .and. l == 3) cycle
           if(k == 5 .and. l == 3) cycle

            f(1,1,k,l-1) = d4(  1,k,l)+d2(2,1,k,l)* 3 +(+d3(2,1,k,l)* 2 +d3(3,1,k,l) +d1(  3,k,l)* 2 +d1(  4,k,l))*qx +(+d2(3,1,k,l)+d2(4,1,k,l)* 2 )*xx +d1(  5,k,l)*xxx
            f(2,1,k,l-1) = d4(  4,k,l)+d2(2,1,k,l) +(+d3(3,4,k,l)+d1(  4,k,l))*qx
            f(3,1,k,l-1) = d4(  6,k,l)+d2(2,1,k,l) +(+d3(3,6,k,l)+d1(  4,k,l))*qx +d3(2,3,k,l)*qzd +d2(4,3,k,l)*xzd +d2(3,1,k,l)*zz +d1(  5,k,l)*xzz
            f(4,1,k,l-1) = d4(  2,k,l)+d2(2,2,k,l) +(+d3(2,2,k,l)+d3(3,2,k,l))*qx +d2(4,2,k,l)*xx
            f(5,1,k,l-1) = d4(  3,k,l)+d2(2,3,k,l) +(+d3(2,3,k,l)+d3(3,3,k,l))*qx +(+d3(2,1,k,l)+d1(  3,k,l))*qz +d2(4,3,k,l)*xx +(+d2(3,1,k,l)+d2(4,1,k,l))*xz +d1(  5,k,l)*xxz
            f(6,1,k,l-1) = d4(  5,k,l) +d3(3,5,k,l)*qx +d3(2,2,k,l)*qz +d2(4,2,k,l)*xz

            f(1,2,k,l-1) = d4(  2,k,l)+d2(2,2,k,l) +d3(2,2,k,l)*qxd +d2(3,2,k,l)*xx
            f(2,2,k,l-1) = d4(  7,k,l)+d2(2,2,k,l)* 3
            f(3,2,k,l-1) = d4(  9,k,l)+d2(2,2,k,l) +d3(2,5,k,l)*qzd +d2(3,2,k,l)*zz
            f(4,2,k,l-1) = d4(  4,k,l)+d2(2,1,k,l) +(+d3(2,4,k,l)+d1(  3,k,l))*qx
            f(5,2,k,l-1) = d4(  5,k,l) +d3(2,5,k,l)*qx +d3(2,2,k,l)*qz +d2(3,2,k,l)*xz
            f(6,2,k,l-1) = d4(  8,k,l)+d2(2,3,k,l) +(+d3(2,4,k,l)+d1(  3,k,l))*qz

            f(1,3,k,l-1) = d4(  3,k,l)+d2(2,3,k,l) +d3(2,3,k,l)*qxd +(+d3(3,1,k,l)+d1(  4,k,l))*qz +d2(3,3,k,l)*xx +d2(4,1,k,l)*xzd +d1(  5,k,l)*xxz
            f(2,3,k,l-1) = d4(  8,k,l)+d2(2,3,k,l) +(+d3(3,4,k,l)+d1(  4,k,l))*qz
            f(3,3,k,l-1) = d4( 10,k,l)+d2(2,3,k,l)* 3 +(+d3(2,6,k,l)* 2 +d3(3,6,k,l) +d1(  3,k,l)* 2 +d1(  4,k,l))*qz +(+d2(3,3,k,l)+d2(4,3,k,l)* 2 )*zz +d1(  5,k,l)*zzz
            f(4,3,k,l-1) = d4(  5,k,l) +d3(2,5,k,l)*qx +d3(3,2,k,l)*qz +d2(4,2,k,l)*xz
            f(5,3,k,l-1) = d4(  6,k,l)+d2(2,1,k,l) +(+d3(2,6,k,l)+d1(  3,k,l))*qx +(+d3(2,3,k,l)+d3(3,3,k,l))*qz +(+d2(3,3,k,l)+d2(4,3,k,l))*xz +d2(4,1,k,l)*zz +d1(  5,k,l)*xzz
            f(6,3,k,l-1) = d4(  9,k,l)+d2(2,2,k,l) +(+d3(2,5,k,l)+d3(3,5,k,l))*qz +d2(4,2,k,l)*zz

         enddo
      enddo

      f(:,:,1,2)= f(:,:,4,1)
      f(:,:,4,2)= f(:,:,2,1)
      f(:,:,5,2)= f(:,:,6,1)

      end subroutine mcdv_19

! >
! >    @brief   dddp case
! >
! >    @details integration of a dddp case
! >
      subroutine mcdv_20(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,6,6,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 6, lx = 4
      real(kind=dp) :: e1(  5,kx,lx),e2(4,3,kx,lx),e3(4,6,kx,lx), &
                       e4(2,10,kx,lx),e5( 15,kx,lx)

      integer :: ind (15, kx, lx) = reshape(&
        [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 4, 7, 8, 11, &
          12, 13, 16, 17, 18, 19, 22, 23, 24, 25, 26, 6, 9, 10, 13, 14, 15, &
          18, 19, 20, 21, 24, 25, 26, 27, 28, 2, 4, 5, 7, 8, 9, 11, 12, 13, &
          14, 16, 17, 18, 19, 20, 3, 5, 6, 8, 9, 10, 12, 13, 14, 15, 17, 18, &
          19, 20, 21, 5, 8, 9, 12, 13, 14, 17, 18, 19, 20, 23, 24, 25, 26, &
          27, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 4, 7, 8, &
          11, 12, 13, 16, 17, 18, 19, 22, 23, 24, 25, 26, 6, 9, 10, 13, 14, &
          15, 18, 19, 20, 21, 24, 25, 26, 27, 28, 2, 4, 5, 7, 8, 9, 11, 12, &
          13, 14, 16, 17, 18, 19, 20, 3, 5, 6, 8, 9, 10, 12, 13, 14, 15, 17, &
          18, 19, 20, 21, 5, 8, 9, 12, 13, 14, 17, 18, 19, 20, 23, 24, 25, &
          26, 27, 2, 4, 5, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19, 20, 7, &
          11, 12, 16, 17, 18, 22, 23, 24, 25, 29, 30, 31, 32, 33, 9, 13, 14, &
          18, 19, 20, 24, 25, 26, 27, 31, 32, 33, 34, 35, 4, 7, 8, 11, 12, &
          13, 16, 17, 18, 19, 22, 23, 24, 25, 26, 5, 8, 9, 12, 13, 14, 17, &
          18, 19, 20, 23, 24, 25, 26, 27, 8, 12, 13, 17, 18, 19, 23, 24, 25, &
          26, 30, 31, 32, 33, 34, 3, 5, 6, 8, 9, 10, 12, 13, 14, 15, 17, 18, &
          19, 20, 21, 8, 12, 13, 17, 18, 19, 23, 24, 25, 26, 30, 31, 32, 33, &
          34, 10, 14, 15, 19, 20, 21, 25, 26, 27, 28, 32, 33, 34, 35, 36, 5, &
          8, 9, 12, 13, 14, 17, 18, 19, 20, 23, 24, 25, 26, 27, 6, 9, 10, &
          13, 14, 15, 18, 19, 20, 21, 24, 25, 26, 27, 28, 9, 13, 14, 18, 19, &
          20, 24, 25, 26, 27, 31, 32, 33, 34, 35] &
        , shape(ind))
      integer :: i, j, k, l, m, n, ii, jj
      real(kind=dp) :: xx, zz, xz, xxx, xxz, xzz, zzz
      real(kind=dp) :: xxxx, xxzz, xxxz, xzzz, zzzz
      real(kind=dp) :: xxxd, xxzd, xzzd, zzzd
      real(kind=dp) :: xzq
      real(kind=dp) :: qxd, qzd, xzd

      do j=1,15
         e5(  j,1,1)=+rdat%r06( j,2)+rdat%r04( j,1)
      enddo

      do j=1,10
         e4(1,j,1,1)=+rdat%r05(j, 5)+rdat%r03(j, 1)
         e4(2,j,1,1)=+rdat%r05(j, 6)+rdat%r03(j, 2)
      enddo

      do j=1,6
         m= in6(j)
         e3(1,j,1,1)=+rdat%r04(j,14)+rdat%r02(m, 1)
         e3(2,j,1,1)=+rdat%r04(j,15)+rdat%r02(m, 2)
         e3(3,j,1,1)=+rdat%r04(j,16)+rdat%r02(m, 3)
         e3(4,j,1,1)=+rdat%r04(j,17)+rdat%r02(m, 4)
      enddo

      do j=1,3
         e2(1,j,1,1)=+rdat%r03(j,27)+rdat%r01(j, 1)
         e2(2,j,1,1)=+rdat%r03(j,28)+rdat%r01(j, 2)
         e2(3,j,1,1)=+rdat%r03(j,29)+rdat%r01(j, 3)
         e2(4,j,1,1)=+rdat%r03(j,30)+rdat%r01(j, 4)
      enddo

         j=1
      do i=1,5
         e1(i  ,1,1)=+rdat%r02(j,i+36)+rdat%r00(i,1)
      enddo
      do j=1,15
         k= ind(j,2,1)
         e5(  j,2,1)=+rdat%r06( k,2)+rdat%r04( j,1)
      enddo

      do j=1,10
         k= ind(j,2,1)
         e4(1,j,2,1)=+rdat%r05(k, 5)+rdat%r03(j, 1)
         e4(2,j,2,1)=+rdat%r05(k, 6)+rdat%r03(j, 2)
      enddo

      do j=1,6
         k= ind(j,2,1)
         m= in6(j)
         e3(1,j,2,1)=+rdat%r04(k,14)+rdat%r02(m, 1)
         e3(2,j,2,1)=+rdat%r04(k,15)+rdat%r02(m, 2)
         e3(3,j,2,1)=+rdat%r04(k,16)+rdat%r02(m, 3)
         e3(4,j,2,1)=+rdat%r04(k,17)+rdat%r02(m, 4)
      enddo

      do j=1,3
         k= ind(j,2,1)
         e2(1,j,2,1)=+rdat%r03(k,27)+rdat%r01(j, 1)
         e2(2,j,2,1)=+rdat%r03(k,28)+rdat%r01(j, 2)
         e2(3,j,2,1)=+rdat%r03(k,29)+rdat%r01(j, 3)
         e2(4,j,2,1)=+rdat%r03(k,30)+rdat%r01(j, 4)
      enddo

         j=1
         k= ind(j,2,1)
         m= in6(k)
      do i=1,5
         e1(i  ,2,1)=+rdat%r02(m,i+36)+rdat%r00(i,1)
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,3,1)
            l= k-2-jj

      e5(  j,3,1)=+rdat%r06(k,2)-rdat%r05(l,2)* 2 +rdat%r04(j,1)+rdat%r04(j,2)

            if(jj > 4) cycle
      e4(1,j,3,1)=+rdat%r05(k,5)-rdat%r04(l,8)* 2 +rdat%r03(j,1)+rdat%r03(j,3)
      e4(2,j,3,1)=+rdat%r05(k,6)-rdat%r04(l,9)* 2 +rdat%r03(j,2)+rdat%r03(j,4)

            if(jj > 3) cycle
            m= in6(j)
            do i=1,4
      e3(i,j,3,1)=+rdat%r04(k,i+13)-rdat%r03(l,i+10)* 2 +rdat%r02(m,i   )+rdat%r02(m,i+ 4)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      e2(i,j,3,1)=+rdat%r03(k,i+26)-rdat%r02(m,i+20)* 2 +rdat%r01(j,i   )+rdat%r01(j,i+ 4)
            enddo

            if(jj > 1) cycle
            m= in6(k)
            do i=1,5
      e1(i  ,3,1)=+rdat%r02(m,i+36)-rdat%r01(l,i+20)* 2 +rdat%r00(i,1)+rdat%r00(i,2)
            enddo
         enddo
      enddo
      do j=1,15
         k= ind(j,4,1)
         e5(  j,4,1)=+rdat%r06( k,2)
      enddo

      do j=1,10
         k= ind(j,4,1)
         e4(1,j,4,1)=+rdat%r05(k, 5)
         e4(2,j,4,1)=+rdat%r05(k, 6)
      enddo

      do j=1,6
         k= ind(j,4,1)
         e3(1,j,4,1)=+rdat%r04(k,14)
         e3(2,j,4,1)=+rdat%r04(k,15)
         e3(3,j,4,1)=+rdat%r04(k,16)
         e3(4,j,4,1)=+rdat%r04(k,17)
      enddo

      do j=1,3
         k= ind(j,4,1)
         e2(1,j,4,1)=+rdat%r03(k,27)
         e2(2,j,4,1)=+rdat%r03(k,28)
         e2(3,j,4,1)=+rdat%r03(k,29)
         e2(4,j,4,1)=+rdat%r03(k,30)
      enddo

         j=1
         k= ind(j,4,1)
         m= in6(k)
      do i=1,5
         e1(i  ,4,1)=+rdat%r02(m,i+36)
      enddo
      do j=1,15
         k= ind(j,5,1)
         e5(  j,5,1)=+rdat%r06( k,2)-rdat%r05( j,2)
      enddo

      do j=1,10
         k= ind(j,5,1)
         e4(1,j,5,1)=+rdat%r05(k, 5)-rdat%r04(j, 8)
         e4(2,j,5,1)=+rdat%r05(k, 6)-rdat%r04(j, 9)
      enddo

      do j=1,6
         k= ind(j,5,1)
         e3(1,j,5,1)=+rdat%r04(k,14)-rdat%r03(j,11)
         e3(2,j,5,1)=+rdat%r04(k,15)-rdat%r03(j,12)
         e3(3,j,5,1)=+rdat%r04(k,16)-rdat%r03(j,13)
         e3(4,j,5,1)=+rdat%r04(k,17)-rdat%r03(j,14)
      enddo

      do j=1,3
         k= ind(j,5,1)
         m= in6(j)
         e2(1,j,5,1)=+rdat%r03(k,27)-rdat%r02(m,21)
         e2(2,j,5,1)=+rdat%r03(k,28)-rdat%r02(m,22)
         e2(3,j,5,1)=+rdat%r03(k,29)-rdat%r02(m,23)
         e2(4,j,5,1)=+rdat%r03(k,30)-rdat%r02(m,24)
      enddo

         j=1
         k= ind(j,5,1)
         m= in6(k)
      do i=1,5
         e1(i  ,5,1)=+rdat%r02(m,i+36)-rdat%r01(j,i+20)
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,6,1)
            l= k-2-jj

      e5(  j,6,1)=+rdat%r06(k,2)-rdat%r05(l,2)

            if(jj > 4) cycle
      e4(1,j,6,1)=+rdat%r05(k,5)-rdat%r04(l,8)
      e4(2,j,6,1)=+rdat%r05(k,6)-rdat%r04(l,9)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,6,1)=+rdat%r04(k,i+13)-rdat%r03(l,i+10)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      e2(i,j,6,1)=+rdat%r03(k,i+26)-rdat%r02(m,i+20)
            enddo

            if(jj > 1) cycle
            m= in6(k)
            do i=1,5
      e1(i  ,6,1)=+rdat%r02(m,i+36)-rdat%r01(l,i+20)
            enddo
         enddo
      enddo
      do j=1,15
         e5(  j,1,2)=-rdat%r07( j,1)-rdat%r05( j,4)* 3
      enddo

      do j=1,10
         e4(1,j,1,2)=-rdat%r06(j, 4)-rdat%r04(j,12)* 3
         e4(2,j,1,2)=-rdat%r06(j, 5)-rdat%r04(j,13)* 3
      enddo

      do j=1,6
         e3(1,j,1,2)=-rdat%r05(j,11)-rdat%r03(j,19)* 3
         e3(2,j,1,2)=-rdat%r05(j,12)-rdat%r03(j,20)* 3
         e3(3,j,1,2)=-rdat%r05(j,13)-rdat%r03(j,21)* 3
         e3(4,j,1,2)=-rdat%r05(j,14)-rdat%r03(j,22)* 3
      enddo

      do j=1,3
         m= in6(j)
         e2(1,j,1,2)=-rdat%r04(j,26)-rdat%r02(m,29)* 3
         e2(2,j,1,2)=-rdat%r04(j,27)-rdat%r02(m,30)* 3
         e2(3,j,1,2)=-rdat%r04(j,28)-rdat%r02(m,31)* 3
         e2(4,j,1,2)=-rdat%r04(j,29)-rdat%r02(m,32)* 3
      enddo

         j=1
      do i=1,5
         e1(i  ,1,2)=-rdat%r03(j,i+38)-rdat%r01(j,i+30)* 3
      enddo
      do j=1,15
         k= ind(j,2,2)
         e5(  j,2,2)=-rdat%r07( k,1)-rdat%r05( j,4)
      enddo

      do j=1,10
         k= ind(j,2,2)
         e4(1,j,2,2)=-rdat%r06(k, 4)-rdat%r04(j,12)
         e4(2,j,2,2)=-rdat%r06(k, 5)-rdat%r04(j,13)
      enddo

      do j=1,6
         k= ind(j,2,2)
         e3(1,j,2,2)=-rdat%r05(k,11)-rdat%r03(j,19)
         e3(2,j,2,2)=-rdat%r05(k,12)-rdat%r03(j,20)
         e3(3,j,2,2)=-rdat%r05(k,13)-rdat%r03(j,21)
         e3(4,j,2,2)=-rdat%r05(k,14)-rdat%r03(j,22)
      enddo

      do j=1,3
         k= ind(j,2,2)
         m= in6(j)
         e2(1,j,2,2)=-rdat%r04(k,26)-rdat%r02(m,29)
         e2(2,j,2,2)=-rdat%r04(k,27)-rdat%r02(m,30)
         e2(3,j,2,2)=-rdat%r04(k,28)-rdat%r02(m,31)
         e2(4,j,2,2)=-rdat%r04(k,29)-rdat%r02(m,32)
      enddo

         j=1
         k= ind(j,2,2)
      do i=1,5
         e1(i  ,2,2)=-rdat%r03(k,i+38)-rdat%r01(j,i+30)
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,3,2)
            l= k-2-jj

      e5(  j,3,2)=-rdat%r07(k,1)+rdat%r06(l,3)* 2 -rdat%r05(j,3)-rdat%r05(j,4)

            if(jj > 4) cycle
      e4(1,j,3,2)=-rdat%r06(k,4)+rdat%r05(l,7)* 2 -rdat%r04(j,10)-rdat%r04(j,12)
      e4(2,j,3,2)=-rdat%r06(k,5)+rdat%r05(l,8)* 2 -rdat%r04(j,11)-rdat%r04(j,13)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,3,2)=-rdat%r05(k,i+10)+rdat%r04(l,i+17)* 2 -rdat%r03(j,i+14)-rdat%r03(j,i+18)
            enddo

            if(jj > 2) cycle
            m= in6(j)
            do i=1,4
      e2(i,j,3,2)=-rdat%r04(k,i+25)+rdat%r03(l,i+30)* 2 -rdat%r02(m,i+24)-rdat%r02(m,i+28)
            enddo

            if(jj > 1) cycle
            m= in6(l)
            do i=1,5
      e1(i  ,3,2)=-rdat%r03(k,i+38)+rdat%r02(m,i+41)* 2 -rdat%r01(j,i+25)-rdat%r01(j,i+30)
            enddo
         enddo
      enddo
      do j=1,15
         k= ind(j,4,2)
         e5(  j,4,2)=-rdat%r07( k,1)-rdat%r05( k,4)
      enddo

      do j=1,10
         k= ind(j,4,2)
         e4(1,j,4,2)=-rdat%r06(k, 4)-rdat%r04(k,12)
         e4(2,j,4,2)=-rdat%r06(k, 5)-rdat%r04(k,13)
      enddo

      do j=1,6
         k= ind(j,4,2)
         e3(1,j,4,2)=-rdat%r05(k,11)-rdat%r03(k,19)
         e3(2,j,4,2)=-rdat%r05(k,12)-rdat%r03(k,20)
         e3(3,j,4,2)=-rdat%r05(k,13)-rdat%r03(k,21)
         e3(4,j,4,2)=-rdat%r05(k,14)-rdat%r03(k,22)
      enddo

      do j=1,3
         k= ind(j,4,2)
         m= in6(k)
         e2(1,j,4,2)=-rdat%r04(k,26)-rdat%r02(m,29)
         e2(2,j,4,2)=-rdat%r04(k,27)-rdat%r02(m,30)
         e2(3,j,4,2)=-rdat%r04(k,28)-rdat%r02(m,31)
         e2(4,j,4,2)=-rdat%r04(k,29)-rdat%r02(m,32)
      enddo

         j=1
         k= ind(j,4,2)
      do i=1,5
         e1(i  ,4,2)=-rdat%r03(k,i+38)-rdat%r01(k,i+30)
      enddo
      do j=1,15
         k= ind(j,5,2)
         e5(  j,5,2)=-rdat%r07( k,1)+rdat%r06( j,3)-rdat%r05( k,4)+rdat%r04( j,5)
      enddo

      do j=1,10
         k= ind(j,5,2)
         e4(1,j,5,2)=-rdat%r06(k, 4)+rdat%r05(j, 7)-rdat%r04(k,12)+rdat%r03(j, 9)
         e4(2,j,5,2)=-rdat%r06(k, 5)+rdat%r05(j, 8)-rdat%r04(k,13)+rdat%r03(j,10)
      enddo

      do j=1,6
         k= ind(j,5,2)
         m= in6(j)
         e3(1,j,5,2)=-rdat%r05(k,11)+rdat%r04(j,18)-rdat%r03(k,19)+rdat%r02(m,17)
         e3(2,j,5,2)=-rdat%r05(k,12)+rdat%r04(j,19)-rdat%r03(k,20)+rdat%r02(m,18)
         e3(3,j,5,2)=-rdat%r05(k,13)+rdat%r04(j,20)-rdat%r03(k,21)+rdat%r02(m,19)
         e3(4,j,5,2)=-rdat%r05(k,14)+rdat%r04(j,21)-rdat%r03(k,22)+rdat%r02(m,20)
      enddo

      do j=1,3
         k= ind(j,5,2)
         m= in6(k)
         e2(1,j,5,2)=-rdat%r04(k,26)+rdat%r03(j,31)-rdat%r02(m,29)+rdat%r01(j,17)
         e2(2,j,5,2)=-rdat%r04(k,27)+rdat%r03(j,32)-rdat%r02(m,30)+rdat%r01(j,18)
         e2(3,j,5,2)=-rdat%r04(k,28)+rdat%r03(j,33)-rdat%r02(m,31)+rdat%r01(j,19)
         e2(4,j,5,2)=-rdat%r04(k,29)+rdat%r03(j,34)-rdat%r02(m,32)+rdat%r01(j,20)
      enddo

         j=1
         k= ind(j,5,2)
      do i=1,5
         e1(i  ,5,2)=-rdat%r03(k,i+38)+rdat%r02(j,i+41)-rdat%r01(k,i+30)+rdat%r00(i,5)
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,6,2)
            l= k-2-jj

      e5(  j,6,2)=-rdat%r07(k,1)+rdat%r06(l,3)

            if(jj > 4) cycle
      e4(1,j,6,2)=-rdat%r06(k,4)+rdat%r05(l,7)
      e4(2,j,6,2)=-rdat%r06(k,5)+rdat%r05(l,8)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,6,2)=-rdat%r05(k,i+10)+rdat%r04(l,i+17)
            enddo

            if(jj > 2) cycle
            do i=1,4
      e2(i,j,6,2)=-rdat%r04(k,i+25)+rdat%r03(l,i+30)
            enddo

            if(jj > 1) cycle
            m= in6(l)
            do i=1,5
      e1(i  ,6,2)=-rdat%r03(k,i+38)+rdat%r02(m,i+41)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,2,3)
            l= k-3-jj-jj

      e5(  j,2,3)=-rdat%r07(k,1)-rdat%r05(l,4)* 3

            if(jj > 4) cycle
      e4(1,j,2,3)=-rdat%r06(k,4)-rdat%r04(l,12)* 3
      e4(2,j,2,3)=-rdat%r06(k,5)-rdat%r04(l,13)* 3

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,2,3)=-rdat%r05(k,i+10)-rdat%r03(l,i+18)* 3
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      e2(i,j,2,3)=-rdat%r04(k,i+25)-rdat%r02(m,i+28)* 3
            enddo

            if(jj > 1) cycle
            do i=1,5
      e1(i  ,2,3)=-rdat%r03(k,i+38)-rdat%r01(l,i+30)* 3
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,3,3)
            l= k-3-jj
            m= l-2-jj

      e5(  j,3,3)=-rdat%r07(k,1)+rdat%r06(l,3)* 2 -rdat%r05(m,3)-rdat%r05(m,4)

            if(jj > 4) cycle
      e4(1,j,3,3)=-rdat%r06(k,4)+rdat%r05(l,7)* 2 -rdat%r04(m,10)-rdat%r04(m,12)
      e4(2,j,3,3)=-rdat%r06(k,5)+rdat%r05(l,8)* 2 -rdat%r04(m,11)-rdat%r04(m,13)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,3,3)=-rdat%r05(k,i+10)+rdat%r04(l,i+17)* 2 -rdat%r03(m,i+14)-rdat%r03(m,i+18)
            enddo

            if(jj > 2) cycle
            n= in6(m)
            do i=1,4
      e2(i,j,3,3)=-rdat%r04(k,i+25)+rdat%r03(l,i+30)* 2 -rdat%r02(n,i+24)-rdat%r02(n,i+28)
            enddo

            if(jj > 1) cycle
            n= in6(l)
            do i=1,5
      e1(i  ,3,3)=-rdat%r03(k,i+38)+rdat%r02(n,i+41)* 2 -rdat%r01(m,i+25)-rdat%r01(m,i+30)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,6,3)
            l= k-3-jj
            m= l-jj

      e5(  j,6,3)=-rdat%r07(k,1)+rdat%r06(l,3)-rdat%r05(m,4)+rdat%r04(j,5)

            if(jj > 4) cycle
      e4(1,j,6,3)=-rdat%r06(k,4)+rdat%r05(l,7)-rdat%r04(m,12)+rdat%r03(j, 9)
      e4(2,j,6,3)=-rdat%r06(k,5)+rdat%r05(l,8)-rdat%r04(m,13)+rdat%r03(j,10)

            if(jj > 3) cycle
            n= in6(j)
            do i=1,4
      e3(i,j,6,3)=-rdat%r05(k,i+10)+rdat%r04(l,i+17)-rdat%r03(m,i+18)+rdat%r02(n,i+16)
            enddo

            if(jj > 2) cycle
            n= in6(m)
            do i=1,4
      e2(i,j,6,3)=-rdat%r04(k,i+25)+rdat%r03(l,i+30)-rdat%r02(n,i+28)+rdat%r01(j,i+16)
            enddo

            if(jj > 1) cycle
            n= in6(l)
            do i=1,5
      e1(i  ,6,3)=-rdat%r03(k,i+38)+rdat%r02(n,i+41)-rdat%r01(m,i+30)+rdat%r00(i,5)
            enddo
         enddo
      enddo
      do j=1,15
         k= ind(j,1,4)
         e5(  j,1,4)=-rdat%r07( k,1)+rdat%r06( j,1)-rdat%r05( k,4)+rdat%r04( j,4)
      enddo

      do j=1,10
         k= ind(j,1,4)
         e4(1,j,1,4)=-rdat%r06(k, 4)+rdat%r05(j, 9)-rdat%r04(k,12)+rdat%r03(j, 7)
         e4(2,j,1,4)=-rdat%r06(k, 5)+rdat%r05(j,10)-rdat%r04(k,13)+rdat%r03(j, 8)
      enddo

      do j=1,6
         k= ind(j,1,4)
         m= in6(j)
         e3(1,j,1,4)=-rdat%r05(k,11)+rdat%r04(j,22)-rdat%r03(k,19)+rdat%r02(m,13)
         e3(2,j,1,4)=-rdat%r05(k,12)+rdat%r04(j,23)-rdat%r03(k,20)+rdat%r02(m,14)
         e3(3,j,1,4)=-rdat%r05(k,13)+rdat%r04(j,24)-rdat%r03(k,21)+rdat%r02(m,15)
         e3(4,j,1,4)=-rdat%r05(k,14)+rdat%r04(j,25)-rdat%r03(k,22)+rdat%r02(m,16)
      enddo

      do j=1,3
         k= ind(j,1,4)
         m= in6(k)
         e2(1,j,1,4)=-rdat%r04(k,26)+rdat%r03(j,35)-rdat%r02(m,29)+rdat%r01(j,13)
         e2(2,j,1,4)=-rdat%r04(k,27)+rdat%r03(j,36)-rdat%r02(m,30)+rdat%r01(j,14)
         e2(3,j,1,4)=-rdat%r04(k,28)+rdat%r03(j,37)-rdat%r02(m,31)+rdat%r01(j,15)
         e2(4,j,1,4)=-rdat%r04(k,29)+rdat%r03(j,38)-rdat%r02(m,32)+rdat%r01(j,16)
      enddo

         j=1
         k= ind(j,1,4)
      do i=1,5
         e1(i  ,1,4)=-rdat%r03(k,i+38)+rdat%r02(j,i+46)-rdat%r01(k,i+30)+rdat%r00(i,4)
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,2,4)
            l= k-3-jj
            m= l-jj

      e5(  j,2,4)=-rdat%r07(k,1)+rdat%r06(l,1)-rdat%r05(m,4)+rdat%r04(j,4)

            if(jj > 4) cycle
      e4(1,j,2,4)=-rdat%r06(k,4)+rdat%r05(l, 9)-rdat%r04(m,12)+rdat%r03(j,7)
      e4(2,j,2,4)=-rdat%r06(k,5)+rdat%r05(l,10)-rdat%r04(m,13)+rdat%r03(j,8)

            if(jj > 3) cycle
            n= in6(j)
            do i=1,4
      e3(i,j,2,4)=-rdat%r05(k,i+10)+rdat%r04(l,i+21)-rdat%r03(m,i+18)+rdat%r02(n,i+12)
            enddo

            if(jj > 2) cycle
            n= in6(m)
            do i=1,4
      e2(i,j,2,4)=-rdat%r04(k,i+25)+rdat%r03(l,i+34)-rdat%r02(n,i+28)+rdat%r01(j,i+12)
            enddo

            if(jj > 1) cycle
            n= in6(l)
            do i=1,5
      e1(i  ,2,4)=-rdat%r03(k,i+38)+rdat%r02(n,i+46)-rdat%r01(m,i+30)+rdat%r00(i,4)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,3,4)
            l= k-3-jj
            m= l-2-jj

      e5(  j,3,4)=-rdat%r07(k,1)    +rdat%r06(l,1)+rdat%r06(l,3)* 2                       &
     &            -rdat%r05(m,1)* 2 -rdat%r05(m,3)-rdat%r05(m,4)* 3                       &
     &            +rdat%r04(j,3)    +rdat%r04(j,4)+rdat%r04(j,5)* 2

            if(jj > 4) cycle
            do i=1,2
      e4(i,j,3,4)=-rdat%r06(k,i+ 3)    +rdat%r05(l,i+ 6)* 2 +rdat%r05(l,i+ 8)             &
     &            -rdat%r04(m,i+ 5)* 2 -rdat%r04(m,i+ 9)    -rdat%r04(m,i+11)* 3          &
     &            +rdat%r03(j,i+ 4)    +rdat%r03(j,i+ 6)    +rdat%r03(j,i+ 8)* 2
            enddo

            if(jj > 3) cycle
            n= in6(j)
            do i=1,4
      e3(i,j,3,4)=-rdat%r05(k,i+10)+rdat%r04(l,i+17)* 2 +rdat%r04(l,i+21)                 &
     &            -rdat%r03(m,i+14)-rdat%r03(m,i+18)* 3 -rdat%r03(m,i+22)* 2              &
     &            +rdat%r02(n,i+ 8)+rdat%r02(n,i+12)    +rdat%r02(n,i+16)* 2
            enddo

            if(jj > 2) cycle
            n= in6(m)
            do i=1,4
      e2(i,j,3,4)=-rdat%r04(k,i+25)+rdat%r03(l,i+30)* 2 +rdat%r03(l,i+34)                 &
     &            -rdat%r02(n,i+24)-rdat%r02(n,i+28)* 3 -rdat%r02(n,i+32)* 2              &
     &            +rdat%r01(j,i+ 8)+rdat%r01(j,i+12)    +rdat%r01(j,i+16)* 2
            enddo

            if(jj > 1) cycle
            n= in6(l)
            do i=1,5
      e1(i  ,3,4)=-rdat%r03(k,i+38)+rdat%r02(n,i+41)* 2 +rdat%r02(n,i+46)                 &
     &            -rdat%r01(m,i+25)-rdat%r01(m,i+30)* 3 -rdat%r01(m,i+35)* 2              &
     &            +rdat%r00(i,3)+rdat%r00(i,4)    +rdat%r00(i,5)* 2
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,4,4)
            l= k-2-jj

      e5(  j,4,4)=-rdat%r07(k,1)+rdat%r06(l,1)

            if(jj > 4) cycle
      e4(1,j,4,4)=-rdat%r06(k,4)+rdat%r05(l, 9)
      e4(2,j,4,4)=-rdat%r06(k,5)+rdat%r05(l,10)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,4,4)=-rdat%r05(k,i+10)+rdat%r04(l,i+21)
            enddo

            if(jj > 2) cycle
            do i=1,4
      e2(i,j,4,4)=-rdat%r04(k,i+25)+rdat%r03(l,i+34)
            enddo

            if(jj > 1) cycle
            m= in6(l)
            do i=1,5
      e1(i  ,4,4)=-rdat%r03(k,i+38)+rdat%r02(m,i+46)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,5,4)
            l= k-2-jj

      e5(  j,5,4)=-rdat%r07(k,1)+rdat%r06(l,1)+rdat%r06(l,3)                              &
     &                    -rdat%r05(j,1)-rdat%r05(j,4)

            if(jj > 4) cycle
      e4(1,j,5,4)=-rdat%r06(k,4)+rdat%r05(l,7)+rdat%r05(l, 9)                             &
     &                    -rdat%r04(j,6)-rdat%r04(j,12)
      e4(2,j,5,4)=-rdat%r06(k,5)+rdat%r05(l,8)+rdat%r05(l,10)                             &
     &                    -rdat%r04(j,7)-rdat%r04(j,13)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,5,4)=-rdat%r05(k,i+10)+rdat%r04(l,i+17)+rdat%r04(l,i+21)                     &
     &                       -rdat%r03(j,i+18)-rdat%r03(j,i+22)
            enddo

            if(jj > 2) cycle
            m= in6(j)
            do i=1,4
      e2(i,j,5,4)=-rdat%r04(k,i+25)+rdat%r03(l,i+30)+rdat%r03(l,i+34)                     &
     &                       -rdat%r02(m,i+28)-rdat%r02(m,i+32)
            enddo

            if(jj > 1) cycle
            m= in6(l)
            do i=1,5
      e1(i  ,5,4)=-rdat%r03(k,i+38)+rdat%r02(m,i+41)+rdat%r02(m,i+46)                     &
     &                       -rdat%r01(j,i+30)-rdat%r01(j,i+35)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,6,4)
            l= k-3-jj
            m= l-2-jj

      e5(  j,6,4)=-rdat%r07(k,1)+rdat%r06(l,1)+rdat%r06(l,3)                              &
     &                    -rdat%r05(m,1)-rdat%r05(m,4)

            if(jj > 4) cycle
      e4(1,j,6,4)=-rdat%r06(k,4)+rdat%r05(l,7)+rdat%r05(l, 9)                             &
     &                    -rdat%r04(m,6)-rdat%r04(m,12)
      e4(2,j,6,4)=-rdat%r06(k,5)+rdat%r05(l,8)+rdat%r05(l,10)                             &
     &                    -rdat%r04(m,7)-rdat%r04(m,13)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,6,4)=-rdat%r05(k,i+10)+rdat%r04(l,i+17)+rdat%r04(l,i+21)                     &
     &                       -rdat%r03(m,i+18)-rdat%r03(m,i+22)
            enddo

            if(jj > 2) cycle
            n= in6(m)
            do i=1,4
      e2(i,j,6,4)=-rdat%r04(k,i+25)+rdat%r03(l,i+30)+rdat%r03(l,i+34)                     &
     &                       -rdat%r02(n,i+28)-rdat%r02(n,i+32)
            enddo

            if(jj > 1) cycle
            n= in6(l)
            do i=1,5
      e1(i  ,6,4)=-rdat%r03(k,i+38)+rdat%r02(n,i+41)+rdat%r02(n,i+46)                     &
     &                       -rdat%r01(m,i+30)-rdat%r01(m,i+35)
            enddo
         enddo
      enddo

      xx= qx*qx
      zz= qz*qz
      xz= qx*qz
      xxx= xx*qx
      xxz= xx*qz
      xzz= zz*qx
      zzz= zz*qz
      xxxx= xx*xx
      xxxz= xx*xz
      xxzz= xx*zz
      xzzz= xz*zz
      zzzz= zz*zz

      qxd= qx+qx
      qzd= qz+qz
      xzd= xz+xz
      xxxd= xxx+xxx
      xxzd= xxz+xxz
      xzzd= xzz+xzz
      zzzd= zzz+zzz

      xzq= xz+xz+xz+xz

      do l=2,lx
         do k=1,kx

            if(k == 1 .and. l == 3) cycle
            if(k == 4 .and. l == 3) cycle
            if(k == 5 .and. l == 3) cycle

            f(1,1,k,l-1) = e5(  1,k,l)+e3(1,1,k,l)* 6 +e1(  1,k,l)* 3 +(+e4(1,1,k,l)+e4(2,1,k,l)+ (+e2(1,1,k,l)+e2(2,1,k,l))* 3 )*qxd +(+e3(2,1,k,l)+e3(3,1,k,l)* 4 +e3(4,1,k,l) +e1(  2,k,l)+e1(  3,k,l)* 4 +e1(  4,k,l))*xx +(+e2(3,1,k,l)+e2(4,1,k,l))*xxxd +e1(  5,k,l)*xxxx
            f(2,1,k,l-1) = e5(  4,k,l)+e3(1,4,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(2,4,k,l)+e2(2,1,k,l))*qxd +(+e3(4,4,k,l)+e1(  4,k,l))*xx
            f(3,1,k,l-1) = e5(  6,k,l)+e3(1,6,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(2,6,k,l)+e2(2,1,k,l))*qxd +(+e4(1,3,k,l)+e2(1,3,k,l))*qzd +(+e3(4,6,k,l)+e1(  4,k,l))*xx +e3(3,3,k,l)*xzq +(+e3(2,1,k,l)+e1(  2,k,l))*zz +e2(4,3,k,l)*xxzd +e2(3,1,k,l)*xzzd +e1(  5,k,l)*xxzz
            f(4,1,k,l-1) = e5(  2,k,l)+e3(1,2,k,l)* 3 +(+e4(1,2,k,l)+e4(2,2,k,l)* 2 +e2(1,2,k,l)+e2(2,2,k,l)* 2 )*qx +(+e3(3,2,k,l)* 2 +e3(4,2,k,l))*xx +e2(4,2,k,l)*xxx
            f(5,1,k,l-1) = e5(  3,k,l)+e3(1,3,k,l)* 3 +(+e4(1,3,k,l)+e4(2,3,k,l)* 2 +e2(1,3,k,l)+e2(2,3,k,l)* 2 )*qx +(+e4(1,1,k,l)+e2(1,1,k,l)* 3 )*qz +(+e3(3,3,k,l)* 2 +e3(4,3,k,l))*xx +(+e3(2,1,k,l)+e3(3,1,k,l)* 2 +e1(  2,k,l)+e1(  3,k,l)* 2 )*xz +e2(4,3,k,l)*xxx +(+e2(3,1,k,l)* 2 +e2(4,1,k,l))*xxz +e1(  5,k,l)*xxxz
            f(6,1,k,l-1) = e5(  5,k,l)+e3(1,5,k,l) +e4(2,5,k,l)*qxd +(+e4(1,2,k,l)+e2(1,2,k,l))*qz +e3(4,5,k,l)*xx +e3(3,2,k,l)*xzd +e2(4,2,k,l)*xxz

            f(1,2,k,l-1) = e5(  4,k,l)+e3(1,4,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(1,4,k,l)+e2(1,1,k,l))*qxd +(+e3(2,4,k,l)+e1(  2,k,l))*xx
            f(2,2,k,l-1) = e5( 11,k,l)+e3(1,4,k,l)* 6 +e1(  1,k,l)* 3
            f(3,2,k,l-1) = e5( 13,k,l)+e3(1,6,k,l)+e3(1,4,k,l)+e1(  1,k,l) +(+e4(1,8,k,l)+e2(1,3,k,l))*qzd +(+e3(2,4,k,l)+e1(  2,k,l))*zz
            f(4,2,k,l-1) = e5(  7,k,l)+e3(1,2,k,l)* 3 +(+e4(1,7,k,l)+e2(1,2,k,l)* 3 )*qx
            f(5,2,k,l-1) = e5(  8,k,l)+e3(1,3,k,l) +(+e4(1,8,k,l)+e2(1,3,k,l))*qx +(+e4(1,4,k,l)+e2(1,1,k,l))*qz +(+e3(2,4,k,l)+e1(  2,k,l))*xz
            f(6,2,k,l-1) = e5( 12,k,l)+e3(1,5,k,l)* 3 +(+e4(1,7,k,l)+e2(1,2,k,l)* 3 )*qz

            f(1,3,k,l-1) = e5(  6,k,l)+e3(1,6,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(1,6,k,l)+e2(1,1,k,l))*qxd +(+e4(2,3,k,l)+e2(2,3,k,l))*qzd +(+e3(2,6,k,l)+e1(  2,k,l))*xx +e3(3,3,k,l)*xzq +(+e3(4,1,k,l)+e1(  4,k,l))*zz +e2(3,3,k,l)*xxzd +e2(4,1,k,l)*xzzd +e1(  5,k,l)*xxzz
            f(2,3,k,l-1) = e5( 13,k,l)+e3(1,6,k,l)+e3(1,4,k,l)+e1(  1,k,l) +(+e4(2,8,k,l)+e2(2,3,k,l))*qzd +(+e3(4,4,k,l)+e1(  4,k,l))*zz
            f(3,3,k,l-1) = e5( 15,k,l)+e3(1,6,k,l)* 6 +e1(  1,k,l)* 3 +(+e4(1,10,k,l)+e4(2,10,k,l)+ (+e2(1,3,k,l)+e2(2,3,k,l))* 3 )*qzd +(+e3(2,6,k,l)+e3(3,6,k,l)* 4 +e3(4,6,k,l) +e1(  2,k,l)+e1(  3,k,l)* 4 +e1(  4,k,l))*zz +(+e2(3,3,k,l)+e2(4,3,k,l))*zzzd +e1(  5,k,l)*zzzz
            f(4,3,k,l-1) = e5(  9,k,l)+e3(1,2,k,l) +(+e4(1,9,k,l)+e2(1,2,k,l))*qx +e4(2,5,k,l)*qzd +e3(3,5,k,l)*xzd +e3(4,2,k,l)*zz +e2(4,2,k,l)*xzz
            f(5,3,k,l-1) = e5( 10,k,l)+e3(1,3,k,l)* 3 +(+e4(1,10,k,l)+e2(1,3,k,l)* 3 )*qx +(+e4(1,6,k,l)+e4(2,6,k,l)* 2 +e2(1,1,k,l)+e2(2,1,k,l)* 2 )*qz +(+e3(2,6,k,l)+e3(3,6,k,l)* 2 +e1(  2,k,l)+e1(  3,k,l)* 2 )*xz +(+e3(3,3,k,l)* 2 +e3(4,3,k,l))*zz +(+e2(3,3,k,l)* 2 +e2(4,3,k,l))*xzz +e2(4,1,k,l)*zzz +e1(  5,k,l)*xzzz
            f(6,3,k,l-1) = e5( 14,k,l)+e3(1,5,k,l)* 3 +(+e4(1,9,k,l)+e4(2,9,k,l)* 2 +e2(1,2,k,l)+e2(2,2,k,l)* 2 )*qz +(+e3(3,5,k,l)* 2 +e3(4,5,k,l))*zz +e2(4,2,k,l)*zzz

            f(1,4,k,l-1) = e5(  2,k,l)+e3(1,2,k,l)* 3 +(+e4(1,2,k,l)* 2 +e4(2,2,k,l) +e2(1,2,k,l)* 2 +e2(2,2,k,l))*qx +(+e3(2,2,k,l)+e3(3,2,k,l)* 2 )*xx +e2(3,2,k,l)*xxx
            f(2,4,k,l-1) = e5(  7,k,l)+e3(1,2,k,l)* 3 +(+e4(2,7,k,l)+e2(2,2,k,l)* 3 )*qx
            f(3,4,k,l-1) = e5(  9,k,l)+e3(1,2,k,l) +(+e4(2,9,k,l)+e2(2,2,k,l))*qx +e4(1,5,k,l)*qzd +e3(3,5,k,l)*xzd +e3(2,2,k,l)*zz +e2(3,2,k,l)*xzz
            f(4,4,k,l-1) = e5(  4,k,l)+e3(1,4,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(1,4,k,l)+e4(2,4,k,l) +e2(1,1,k,l)+e2(2,1,k,l))*qx +(+e3(3,4,k,l)+e1(  3,k,l))*xx
            f(5,4,k,l-1) = e5(  5,k,l)+e3(1,5,k,l) +(+e4(1,5,k,l)+e4(2,5,k,l))*qx +(+e4(1,2,k,l)+e2(1,2,k,l))*qz +e3(3,5,k,l)*xx +(+e3(2,2,k,l)+e3(3,2,k,l))*xz +e2(3,2,k,l)*xxz
            f(6,4,k,l-1) = e5(  8,k,l)+e3(1,3,k,l) +(+e4(2,8,k,l)+e2(2,3,k,l))*qx +(+e4(1,4,k,l)+e2(1,1,k,l))*qz +(+e3(3,4,k,l)+e1(  3,k,l))*xz

            f(1,5,k,l-1) = e5(  3,k,l)+e3(1,3,k,l)* 3 +(+e4(1,3,k,l)* 2 +e4(2,3,k,l) +e2(1,3,k,l)* 2 +e2(2,3,k,l))*qx +(+e4(2,1,k,l)+e2(2,1,k,l)* 3 )*qz +(+e3(2,3,k,l)+e3(3,3,k,l)* 2 )*xx +(+e3(3,1,k,l)* 2 +e3(4,1,k,l) +e1(  3,k,l)* 2 +e1(  4,k,l))*xz +e2(3,3,k,l)*xxx +(+e2(3,1,k,l)+e2(4,1,k,l)* 2 )*xxz +e1(  5,k,l)*xxxz
            f(2,5,k,l-1) = e5(  8,k,l)+e3(1,3,k,l) +(+e4(2,8,k,l)+e2(2,3,k,l))*qx +(+e4(2,4,k,l)+e2(2,1,k,l))*qz +(+e3(4,4,k,l)+e1(  4,k,l))*xz
            f(3,5,k,l-1) = e5( 10,k,l)+e3(1,3,k,l)* 3 +(+e4(2,10,k,l)+e2(2,3,k,l)* 3 )*qx +(+e4(1,6,k,l)* 2 +e4(2,6,k,l) +e2(1,1,k,l)* 2 +e2(2,1,k,l))*qz +(+e3(3,6,k,l)* 2 +e3(4,6,k,l) +e1(  3,k,l)* 2 +e1(  4,k,l))*xz +(+e3(2,3,k,l)+e3(3,3,k,l)* 2 )*zz +(+e2(3,3,k,l)+e2(4,3,k,l)* 2 )*xzz +e2(3,1,k,l)*zzz +e1(  5,k,l)*xzzz
            f(4,5,k,l-1) = e5(  5,k,l)+e3(1,5,k,l) +(+e4(1,5,k,l)+e4(2,5,k,l))*qx +(+e4(2,2,k,l)+e2(2,2,k,l))*qz +e3(3,5,k,l)*xx +(+e3(3,2,k,l)+e3(4,2,k,l))*xz +e2(4,2,k,l)*xxz
            f(5,5,k,l-1) = e5(  6,k,l)+e3(1,6,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(1,6,k,l)+e4(2,6,k,l) +e2(1,1,k,l)+e2(2,1,k,l))*qx +(+e4(1,3,k,l)+e4(2,3,k,l) +e2(1,3,k,l)+e2(2,3,k,l))*qz +(+e3(3,6,k,l)+e1(  3,k,l))*xx +(+e3(2,3,k,l)+e3(3,3,k,l)* 2 +e3(4,3,k,l))*xz +(+e3(3,1,k,l)+e1(  3,k,l))*zz +(+e2(3,3,k,l)+e2(4,3,k,l))*xxz +(+e2(3,1,k,l)+e2(4,1,k,l))*xzz +e1(  5,k,l)*xxzz
            f(6,5,k,l-1) = e5(  9,k,l)+e3(1,2,k,l) +(+e4(2,9,k,l)+e2(2,2,k,l))*qx +(+e4(1,5,k,l)+e4(2,5,k,l))*qz +(+e3(3,5,k,l)+e3(4,5,k,l))*xz +e3(3,2,k,l)*zz +e2(4,2,k,l)*xzz

            f(1,6,k,l-1) = e5(  5,k,l)+e3(1,5,k,l) +e4(1,5,k,l)*qxd +(+e4(2,2,k,l)+e2(2,2,k,l))*qz +e3(2,5,k,l)*xx +e3(3,2,k,l)*xzd +e2(3,2,k,l)*xxz
            f(2,6,k,l-1) = e5( 12,k,l)+e3(1,5,k,l)* 3 +(+e4(2,7,k,l)+e2(2,2,k,l)* 3 )*qz
            f(3,6,k,l-1) = e5( 14,k,l)+e3(1,5,k,l)* 3 +(+e4(1,9,k,l)* 2 +e4(2,9,k,l) +e2(1,2,k,l)* 2 +e2(2,2,k,l))*qz +(+e3(2,5,k,l)+e3(3,5,k,l)* 2 )*zz +e2(3,2,k,l)*zzz
            f(4,6,k,l-1) = e5(  8,k,l)+e3(1,3,k,l) +(+e4(1,8,k,l)+e2(1,3,k,l))*qx +(+e4(2,4,k,l)+e2(2,1,k,l))*qz +(+e3(3,4,k,l)+e1(  3,k,l))*xz
            f(5,6,k,l-1) = e5(  9,k,l)+e3(1,2,k,l) +(+e4(1,9,k,l)+e2(1,2,k,l))*qx +(+e4(1,5,k,l)+e4(2,5,k,l))*qz +(+e3(2,5,k,l)+e3(3,5,k,l))*xz +e3(3,2,k,l)*zz +e2(3,2,k,l)*xzz
            f(6,6,k,l-1) = e5( 13,k,l)+e3(1,6,k,l)+e3(1,4,k,l)+e1(  1,k,l) +(+e4(1,8,k,l)+e4(2,8,k,l) +e2(1,3,k,l)+e2(2,3,k,l))*qz +(+e3(3,4,k,l)+e1(  3,k,l))*zz

         enddo
      enddo

      f(:,:,1,2) = f(:,:,4,1)
      f(:,:,4,2) = f(:,:,2,1)
      f(:,:,5,2) = f(:,:,6,1)

      end subroutine mcdv_20

! >
! >    @brief   dddd case
! >
! >    @details integration of a dddd case
! >
      subroutine mcdv_21(f,rdat, qx,qz)

      implicit none
      type(rotaxis_data_t) :: rdat

      real(kind=dp) :: f(6,6,6,*)
      real(kind=dp) :: qx, qz

      integer, parameter :: kx = 6, lx = 6
      real(kind=dp) :: e1(  5,kx,lx),e2(4,3,kx,lx),e3(4,6,kx,lx), &
                       e4(2,10,kx,lx),e5( 15,kx,lx)

      integer, parameter :: ind (15, kx, lx) = reshape(&
        [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 4, 7, 8, 11, &
          12, 13, 16, 17, 18, 19, 22, 23, 24, 25, 26, 6, 9, 10, 13, 14, 15, &
          18, 19, 20, 21, 24, 25, 26, 27, 28, 2, 4, 5, 7, 8, 9, 11, 12, 13, &
          14, 16, 17, 18, 19, 20, 3, 5, 6, 8, 9, 10, 12, 13, 14, 15, 17, 18, &
          19, 20, 21, 5, 8, 9, 12, 13, 14, 17, 18, 19, 20, 23, 24, 25, 26, 27, &
          4, 7, 8, 11, 12, 13, 16, 17, 18, 19, 22, 23, 24, 25, 26, 11, 16, 17, &
          22, 23, 24, 29, 30, 31, 32, 37, 38, 39, 40, 41, 13, 18, 19, 24, 25, &
          26, 31, 32, 33, 34, 39, 40, 41, 42, 43, 7, 11, 12, 16, 17, 18, 22, &
          23, 24, 25, 29, 30, 31, 32, 33, 8, 12, 13, 17, 18, 19, 23, 24, 25, &
          26, 30, 31, 32, 33, 34, 12, 17, 18, 23, 24, 25, 30, 31, 32, 33, 38, &
          39, 40, 41, 42, 6, 9, 10, 13, 14, 15, 18, 19, 20, 21, 24, 25, 26, &
          27, 28, 13, 18, 19, 24, 25, 26, 31, 32, 33, 34, 39, 40, 41, 42, 43, &
          15, 20, 21, 26, 27, 28, 33, 34, 35, 36, 41, 42, 43, 44, 45, 9, 13, &
          14, 18, 19, 20, 24, 25, 26, 27, 31, 32, 33, 34, 35, 10, 14, 15, 19, &
          20, 21, 25, 26, 27, 28, 32, 33, 34, 35, 36, 14, 19, 20, 25, 26, 27, &
          32, 33, 34, 35, 40, 41, 42, 43, 44, 2, 4, 5, 7, 8, 9, 11, 12, 13, &
          14, 16, 17, 18, 19, 20, 7, 11, 12, 16, 17, 18, 22, 23, 24, 25, 29, &
          30, 31, 32, 33, 9, 13, 14, 18, 19, 20, 24, 25, 26, 27, 31, 32, 33, &
          34, 35, 4, 7, 8, 11, 12, 13, 16, 17, 18, 19, 22, 23, 24, 25, 26, 5, &
          8, 9, 12, 13, 14, 17, 18, 19, 20, 23, 24, 25, 26, 27, 8, 12, 13, 17, &
          18, 19, 23, 24, 25, 26, 30, 31, 32, 33, 34, 3, 5, 6, 8, 9, 10, 12, &
          13, 14, 15, 17, 18, 19, 20, 21, 8, 12, 13, 17, 18, 19, 23, 24, 25, &
          26, 30, 31, 32, 33, 34, 10, 14, 15, 19, 20, 21, 25, 26, 27, 28, 32, &
          33, 34, 35, 36, 5, 8, 9, 12, 13, 14, 17, 18, 19, 20, 23, 24, 25, 26, &
          27, 6, 9, 10, 13, 14, 15, 18, 19, 20, 21, 24, 25, 26, 27, 28, 9, 13, &
          14, 18, 19, 20, 24, 25, 26, 27, 31, 32, 33, 34, 35, 5, 8, 9, 12, 13, &
          14, 17, 18, 19, 20, 23, 24, 25, 26, 27, 12, 17, 18, 23, 24, 25, 30, &
          31, 32, 33, 38, 39, 40, 41, 42, 14, 19, 20, 25, 26, 27, 32, 33, 34, &
          35, 40, 41, 42, 43, 44, 8, 12, 13, 17, 18, 19, 23, 24, 25, 26, 30, &
          31, 32, 33, 34, 9, 13, 14, 18, 19, 20, 24, 25, 26, 27, 31, 32, 33, &
          34, 35, 13, 18, 19, 24, 25, 26, 31, 32, 33, 34, 39, 40, 41, 42, 43] &
        , shape(ind))
      integer, parameter :: jnd1(15) = &
        [3, 5, 6, 8, 9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21]
      integer, parameter :: jnd2(15) = &
        [4, 7, 8, 11, 12, 13, 16, 17, 18, 19, 22, 23, 24, 25, 26]
      integer, parameter :: jnd3(15) = &
        [8, 12, 13, 17, 18, 19, 23, 24, 25, 26, 30, 31, 32, 33, 34]
      integer :: i, j, k, l, m, n, ii, jj
      integer :: i4, i5, j1, jn, m1
      real(kind=dp) :: xx, zz, xz, xxx, xxz, xzz, zzz
      real(kind=dp) :: xxxx, xxzz, xxxz, xzzz, zzzz
      real(kind=dp) :: xxxd, xxzd, xzzd, zzzd
      real(kind=dp) :: xzq
      real(kind=dp) :: qxd, qzd, xzd

      real(kind=dp) :: t531e1(  5),t531e2(4,3),t531e3(4,6),t531e4(2,10), &
        t531e5( 15), v531e1(  5),v531e2(4,3),v531e3(4,6),v531e4(2,10), &
        v531e5( 15), u5e1(  5,4),u5e2(4,3,4),u5e3(4,6,4),u5e4(2,10,4), &
        u5e5( 15,4), t632e1(  5),t632e2(4,3),t632e3(4,6),t632e4(2,10), &
        t632e5( 15), v632e1(  5),v632e2(4,3),v632e3(4,6),v632e4(2,10), &
        v632e5( 15), u6e1(  5,4),u6e2(4,3,4),u6e3(4,6,4),u6e4(2,10,4), &
        u6e5( 15,4)
      do j=1,15
         e5(  j,1,1)=+rdat%r08(   j)+rdat%r06( j,1)* 6 +rdat%r04( j,1)* 3
      enddo

      do j=1,10
         e4(1,j,1,1)=+rdat%r07(j, 3)+rdat%r05(j, 5)* 6 +rdat%r03(j, 1)* 3
         e4(2,j,1,1)=+rdat%r07(j, 4)+rdat%r05(j, 6)* 6 +rdat%r03(j, 2)* 3
      enddo

      do j=1,6
         m= in6(j)
         e3(1,j,1,1)=+rdat%r06(j, 9)+rdat%r04(j,14)* 6 +rdat%r02(m, 1)* 3
         e3(2,j,1,1)=+rdat%r06(j,10)+rdat%r04(j,15)* 6 +rdat%r02(m, 2)* 3
         e3(3,j,1,1)=+rdat%r06(j,11)+rdat%r04(j,16)* 6 +rdat%r02(m, 3)* 3
         e3(4,j,1,1)=+rdat%r06(j,12)+rdat%r04(j,17)* 6 +rdat%r02(m, 4)* 3
      enddo

      do j=1,3
         e2(1,j,1,1)=+rdat%r05(j,21)+rdat%r03(j,27)* 6 +rdat%r01(j, 1)* 3
         e2(2,j,1,1)=+rdat%r05(j,22)+rdat%r03(j,28)* 6 +rdat%r01(j, 2)* 3
         e2(3,j,1,1)=+rdat%r05(j,23)+rdat%r03(j,29)* 6 +rdat%r01(j, 3)* 3
         e2(4,j,1,1)=+rdat%r05(j,24)+rdat%r03(j,30)* 6 +rdat%r01(j, 4)* 3
      enddo

         j=1
      do i=1,5
         e1(i  ,1,1)=+rdat%r04(j,i+37)+rdat%r02(j,i+36)* 6 +rdat%r00(i,1)* 3
      enddo
      do j=1,15
         k= ind(j,2,1)
         e5(  j,2,1)=+rdat%r08(   k)+rdat%r06( k,1)+rdat%r06( j,1)+rdat%r04( j,1)
      enddo

      do j=1,10
         k= ind(j,2,1)
         e4(1,j,2,1)=+rdat%r07(k, 3)+rdat%r05(k, 5)+rdat%r05(j, 5)+rdat%r03(j, 1)
         e4(2,j,2,1)=+rdat%r07(k, 4)+rdat%r05(k, 6)+rdat%r05(j, 6)+rdat%r03(j, 2)
      enddo

      do j=1,6
         k= ind(j,2,1)
         m= in6(j)
         e3(1,j,2,1)=+rdat%r06(k, 9)+rdat%r04(k,14)+rdat%r04(j,14)+rdat%r02(m, 1)
         e3(2,j,2,1)=+rdat%r06(k,10)+rdat%r04(k,15)+rdat%r04(j,15)+rdat%r02(m, 2)
         e3(3,j,2,1)=+rdat%r06(k,11)+rdat%r04(k,16)+rdat%r04(j,16)+rdat%r02(m, 3)
         e3(4,j,2,1)=+rdat%r06(k,12)+rdat%r04(k,17)+rdat%r04(j,17)+rdat%r02(m, 4)
      enddo

      do j=1,3
         k= ind(j,2,1)
         e2(1,j,2,1)=+rdat%r05(k,21)+rdat%r03(k,27)+rdat%r03(j,27)+rdat%r01(j, 1)
         e2(2,j,2,1)=+rdat%r05(k,22)+rdat%r03(k,28)+rdat%r03(j,28)+rdat%r01(j, 2)
         e2(3,j,2,1)=+rdat%r05(k,23)+rdat%r03(k,29)+rdat%r03(j,29)+rdat%r01(j, 3)
         e2(4,j,2,1)=+rdat%r05(k,24)+rdat%r03(k,30)+rdat%r03(j,30)+rdat%r01(j, 4)
      enddo

         j=1
         k= ind(j,2,1)
         m= in6(k)
      do i=1,5
         e1(i  ,2,1)=+rdat%r04(k,i+37)+rdat%r02(m,i+36)+rdat%r02(j,i+36)+rdat%r00(i,1)
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,3,1)
            l= k-2-jj

      t531e5(  j)=+rdat%r08(k  )-rdat%r07(l,1)-rdat%r07(l,2)+rdat%r06(j,1)                      &
     &            +rdat%r06(k,1)-rdat%r05(l,1)-rdat%r05(l,2)+rdat%r04(j,1)

            if(jj > 4) cycle
      t531e4(1,j)=+rdat%r07(k,3)-rdat%r06(l,5)-rdat%r06(l,7)+rdat%r05(j,5)                      &
     &            +rdat%r05(k,5)-rdat%r04(l,6)-rdat%r04(l,8)+rdat%r03(j,1)
      t531e4(2,j)=+rdat%r07(k,4)-rdat%r06(l,6)-rdat%r06(l,8)+rdat%r05(j,6)                      &
     &            +rdat%r05(k,6)-rdat%r04(l,7)-rdat%r04(l,9)+rdat%r03(j,2)

            if(jj > 3) cycle
            m= in6(j)
            do i=1,4
      t531e3(i,j)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+12)-rdat%r05(l,i+16)+rdat%r04(j,i+13)          &
     &            +rdat%r04(k,i+13)-rdat%r03(l,i+10)-rdat%r03(l,i+14)+rdat%r02(m,i)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      t531e2(i,j)=+rdat%r05(k,i+20)-rdat%r04(l,i+29)-rdat%r04(l,i+33)+rdat%r03(j,i+26)          &
     &            +rdat%r03(k,i+26)-rdat%r02(m,i+20)-rdat%r02(m,i+24)+rdat%r01(j,i)
            enddo

            if(jj > 1) cycle
            do i=1,5
      t531e1(i  )=+rdat%r04(k,i+37)-rdat%r03(l,i+42)-rdat%r03(l,i+47)+rdat%r02(j,i+36)          &
     &            +rdat%r02(l,i+36)-rdat%r01(l,i+20)-rdat%r01(l,i+25)+rdat%r00(i,1)
            enddo
         enddo
      enddo

      do j1=2,4
         do j=1,15
               u5e5(  j,j1)= rdat%r06(  j,j1)+rdat%r04(  j,j1)
         enddo
         m1=j1+j1+ 2
         do j=1,10
               u5e4(1,j,j1)= rdat%r05(j,1+m1)+rdat%r03(j,1+m1- 4)
               u5e4(2,j,j1)= rdat%r05(j,2+m1)+rdat%r03(j,2+m1- 4)
         enddo
         m1=m1+m1+ 5
         do j=1,6
            m= in6(j)
            do i=1,4
               u5e3(i,j,j1)=+rdat%r04(j,i+m1)+rdat%r02(m,i+m1-13)
            enddo
         enddo
         m1=m1+13
         do j=1,3
            do i=1,4
               u5e2(i,j,j1)= rdat%r03(j,i+m1)+rdat%r01(j,i+m1-26)
            enddo
         enddo
         m1=m1+j1+ 9
            do i=1,5
               u5e1(i  ,j1)= rdat%r02(1,i+m1)+rdat%r00(i,j1)
            enddo
      enddo

      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= jnd1(j)

      v531e5(  j)=+rdat%r07(k,1)-rdat%r07(k,2)+rdat%r05(k,1)-rdat%r05(k,2)

            if(jj > 4) cycle
      v531e4(1,j)=+rdat%r06(k,5)-rdat%r06(k,7)+rdat%r04(k,6)-rdat%r04(k,8)
      v531e4(2,j)=+rdat%r06(k,6)-rdat%r06(k,8)+rdat%r04(k,7)-rdat%r04(k,9)

            if(jj > 3) cycle
            do i=1,4
      v531e3(i,j)=+rdat%r05(k,i+12)-rdat%r05(k,i+16)+rdat%r03(k,i+10)-rdat%r03(k,i+14)
            enddo

            if(jj > 2) cycle
            m= in6(k)
            do i=1,4
      v531e2(i,j)=+rdat%r04(k,i+29)-rdat%r04(k,i+33)+rdat%r02(m,i+20)-rdat%r02(m,i+24)
            enddo

            if(jj > 1) cycle
            do i=1,5
      v531e1(i  )=+rdat%r03(k,i+42)-rdat%r03(k,i+47)+rdat%r01(k,i+20)-rdat%r01(k,i+25)
            enddo
         enddo
      enddo

      j1=2
      do j=1,15
         e5(  j,3,1)= t531e5(  j)+u5e5(  j,j1)-v531e5(  j)
      enddo
      do j=1,10
         e4(1,j,3,1)= t531e4(1,j)+u5e4(1,j,j1)-v531e4(1,j)
         e4(2,j,3,1)= t531e4(2,j)+u5e4(2,j,j1)-v531e4(2,j)
      enddo
      do j=1,6
         e3(1,j,3,1)= t531e3(1,j)+u5e3(1,j,j1)-v531e3(1,j)
         e3(2,j,3,1)= t531e3(2,j)+u5e3(2,j,j1)-v531e3(2,j)
         e3(3,j,3,1)= t531e3(3,j)+u5e3(3,j,j1)-v531e3(3,j)
         e3(4,j,3,1)= t531e3(4,j)+u5e3(4,j,j1)-v531e3(4,j)
      enddo
      do j=1,3
         e2(1,j,3,1)= t531e2(1,j)+u5e2(1,j,j1)-v531e2(1,j)
         e2(2,j,3,1)= t531e2(2,j)+u5e2(2,j,j1)-v531e2(2,j)
         e2(3,j,3,1)= t531e2(3,j)+u5e2(3,j,j1)-v531e2(3,j)
         e2(4,j,3,1)= t531e2(4,j)+u5e2(4,j,j1)-v531e2(4,j)
      enddo
      do i=1,5
         e1(i  ,3,1)= t531e1(i  )+u5e1(i  ,j1)-v531e1(i  )
      enddo
      do j=1,15
         k= ind(j,4,1)
         e5(  j,4,1)=+rdat%r08(   k)+rdat%r06( k,1)* 3
      enddo

      do j=1,10
         k= ind(j,4,1)
         e4(1,j,4,1)=+rdat%r07(k, 3)+rdat%r05(k, 5)* 3
         e4(2,j,4,1)=+rdat%r07(k, 4)+rdat%r05(k, 6)* 3
      enddo

      do j=1,6
         k= ind(j,4,1)
         e3(1,j,4,1)=+rdat%r06(k, 9)+rdat%r04(k,14)* 3
         e3(2,j,4,1)=+rdat%r06(k,10)+rdat%r04(k,15)* 3
         e3(3,j,4,1)=+rdat%r06(k,11)+rdat%r04(k,16)* 3
         e3(4,j,4,1)=+rdat%r06(k,12)+rdat%r04(k,17)* 3
      enddo

      do j=1,3
         k= ind(j,4,1)
         e2(1,j,4,1)=+rdat%r05(k,21)+rdat%r03(k,27)* 3
         e2(2,j,4,1)=+rdat%r05(k,22)+rdat%r03(k,28)* 3
         e2(3,j,4,1)=+rdat%r05(k,23)+rdat%r03(k,29)* 3
         e2(4,j,4,1)=+rdat%r05(k,24)+rdat%r03(k,30)* 3
      enddo

         j=1
         k= ind(j,4,1)
         m= in6(k)
      do i=1,5
         e1(i  ,4,1)=+rdat%r04(k,i+37)+rdat%r02(m,i+36)* 3
      enddo
      do j=1,15
         k= ind(j,5,1)
         e5(  j,5,1)=+rdat%r08(   k)-rdat%r07( j,1)+(rdat%r06( k,1)-rdat%r05( j,1))* 3
      enddo

      do j=1,10
         k= ind(j,5,1)
         e4(1,j,5,1)=+rdat%r07(k, 3)-rdat%r06(j, 5)+(rdat%r05(k, 5)-rdat%r04(j, 6))* 3
         e4(2,j,5,1)=+rdat%r07(k, 4)-rdat%r06(j, 6)+(rdat%r05(k, 6)-rdat%r04(j, 7))* 3
      enddo

      do j=1,6
         k= ind(j,5,1)
         e3(1,j,5,1)=+rdat%r06(k, 9)-rdat%r05(j,13)+(rdat%r04(k,14)-rdat%r03(j,11))* 3
         e3(2,j,5,1)=+rdat%r06(k,10)-rdat%r05(j,14)+(rdat%r04(k,15)-rdat%r03(j,12))* 3
         e3(3,j,5,1)=+rdat%r06(k,11)-rdat%r05(j,15)+(rdat%r04(k,16)-rdat%r03(j,13))* 3
         e3(4,j,5,1)=+rdat%r06(k,12)-rdat%r05(j,16)+(rdat%r04(k,17)-rdat%r03(j,14))* 3
      enddo

      do j=1,3
         k= ind(j,5,1)
         m= in6(j)
         e2(1,j,5,1)=+rdat%r05(k,21)-rdat%r04(j,30)+(rdat%r03(k,27)-rdat%r02(m,21))* 3
         e2(2,j,5,1)=+rdat%r05(k,22)-rdat%r04(j,31)+(rdat%r03(k,28)-rdat%r02(m,22))* 3
         e2(3,j,5,1)=+rdat%r05(k,23)-rdat%r04(j,32)+(rdat%r03(k,29)-rdat%r02(m,23))* 3
         e2(4,j,5,1)=+rdat%r05(k,24)-rdat%r04(j,33)+(rdat%r03(k,30)-rdat%r02(m,24))* 3
      enddo

         j=1
         k= ind(j,5,1)
         m= in6(k)
      do i=1,5
         e1(i  ,5,1)=+rdat%r04(k,i+37)-rdat%r03(j,i+42)+(rdat%r02(m,i+36)-rdat%r01(j,i+20))* 3
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,6,1)
            l= k-2-jj

      e5(  j,6,1)=+rdat%r08(k  )-rdat%r07(l,1)+rdat%r06(k,1)-rdat%r05(l,1)

            if(jj > 4) cycle
      e4(1,j,6,1)=+rdat%r07(k,3)-rdat%r06(l,5)+rdat%r05(k,5)-rdat%r04(l,6)
      e4(2,j,6,1)=+rdat%r07(k,4)-rdat%r06(l,6)+rdat%r05(k,6)-rdat%r04(l,7)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,6,1)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+12)+rdat%r04(k,i+13)-rdat%r03(l,i+10)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      e2(i,j,6,1)=+rdat%r05(k,i+20)-rdat%r04(l,i+29)+rdat%r03(k,i+26)-rdat%r02(m,i+20)
            enddo

            if(jj > 1) cycle
            m= in6(k)
            do i=1,5
      e1(i  ,6,1)=+rdat%r04(k,i+37)-rdat%r03(l,i+42)+rdat%r02(m,i+36)-rdat%r01(l,i+20)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,2,2)
            l= k-5-jj-jj

      e5(  j,2,2)=+rdat%r08(k  )+rdat%r06(l,1)* 6 +rdat%r04(j,1)* 3

            if(jj > 4) cycle
      e4(1,j,2,2)=+rdat%r07(k,3)+rdat%r05(l,5)* 6 +rdat%r03(j,1)* 3
      e4(2,j,2,2)=+rdat%r07(k,4)+rdat%r05(l,6)* 6 +rdat%r03(j,2)* 3

            if(jj > 3) cycle
            m= in6(j)
            do i=1,4
      e3(i,j,2,2)=+rdat%r06(k,i+ 8)+rdat%r04(l,i+13)* 6 +rdat%r02(m,i)* 3
            enddo

            if(jj > 2) cycle
            do i=1,4
      e2(i,j,2,2)=+rdat%r05(k,i+20)+rdat%r03(l,i+26)* 6 +rdat%r01(j,i)* 3
            enddo

            if(jj > 1) cycle
            m= in6(l)
            do i=1,5
      e1(i  ,2,2)=+rdat%r04(k,i+37)+rdat%r02(m,i+36)* 6 +rdat%r00(i,1)* 3
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,3,2)
            l= k-4-jj
            m= l-3-jj
            n= m-jj

      t632e5(  j)=+rdat%r08(k  )-rdat%r07(l,1)-rdat%r07(l,2)+rdat%r06(m,1)                      &
     &          +rdat%r06(m+2,1)-rdat%r05(n,1)-rdat%r05(n,2)+rdat%r04(j,1)

            if(jj > 4) cycle
      t632e4(1,j)=+rdat%r07(k,3)-rdat%r06(l,5)-rdat%r06(l,7)+rdat%r05(m,5)                      &
     &          +rdat%r05(m+2,5)-rdat%r04(n,6)-rdat%r04(n,8)+rdat%r03(j,1)
      t632e4(2,j)=+rdat%r07(k,4)-rdat%r06(l,6)-rdat%r06(l,8)+rdat%r05(m,6)                      &
     &          +rdat%r05(m+2,6)-rdat%r04(n,7)-rdat%r04(n,9)+rdat%r03(j,2)

            if(jj > 3) cycle
            jn= in6(j)
            do i=1,4
      t632e3(i,j)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+12)-rdat%r05(l,i+16)+rdat%r04(m,i+13)          &
     &          +rdat%r04(m+2,i+13)-rdat%r03(n,i+10)-rdat%r03(n,i+14)+rdat%r02(jn,i)
            enddo

            if(jj > 2) cycle
            n= in6(n)
            do i=1,4
      t632e2(i,j)=+rdat%r05(k,i+20)-rdat%r04(l,i+29)-rdat%r04(l,i+33)+rdat%r03(m,i+26)          &
     &          +rdat%r03(m+2,i+26)-rdat%r02(n,i+20)-rdat%r02(n,i+24)+rdat%r01(j,i)
            enddo
            n= m-jj

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
      t632e1(i  )=+rdat%r04(k,i+37)-rdat%r03(l,i+42)-rdat%r03(l,i+47)+rdat%r02(m,i+36)          &
     &          +rdat%r02(m+1,i+36)-rdat%r01(n,i+20)-rdat%r01(n,i+25)+rdat%r00(i,1)
            enddo
         enddo
      enddo

      do j1=2,4
         do j=1,15
            k= jnd2(j)
               u6e5(  j,j1)=+rdat%r06(k  ,j1)+rdat%r04(  j,j1)
         enddo
         m1=j1+j1+ 2
         do j=1,10
            k= jnd2(j)
               u6e4(1,j,j1)=+rdat%r05(k,1+m1)+rdat%r03(j,1+m1- 4)
               u6e4(2,j,j1)=+rdat%r05(k,2+m1)+rdat%r03(j,2+m1- 4)
         enddo
         m1=m1+m1+ 5
         do j=1,6
            k= jnd2(j)
            m= in6(j)
            do i=1,4
               u6e3(i,j,j1)=+rdat%r04(k,i+m1)+rdat%r02(m,i+m1-13)
            enddo
         enddo
         m1=m1+13
         do j=1,3
            k= jnd2(j)
            do i=1,4
               u6e2(i,j,j1)=+rdat%r03(k,i+m1)+rdat%r01(j,i+m1-26)
            enddo
         enddo
         m1=m1+j1+ 9
            do i=1,5
               u6e1(i  ,j1)=+rdat%r02(2,i+m1)+rdat%r00(i,j1)
            enddo
      enddo

      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= jnd3(j)
            l= k-3-jj-jj

      v632e5(  j)=+rdat%r07(k,1)-rdat%r07(k,2)+rdat%r05(l,1)-rdat%r05(l,2)

            if(jj > 4) cycle
      v632e4(1,j)=+rdat%r06(k,5)-rdat%r06(k,7)+rdat%r04(l,6)-rdat%r04(l,8)
      v632e4(2,j)=+rdat%r06(k,6)-rdat%r06(k,8)+rdat%r04(l,7)-rdat%r04(l,9)

            if(jj > 3) cycle
            do i=1,4
      v632e3(i,j)=+rdat%r05(k,i+12)-rdat%r05(k,i+16)+rdat%r03(l,i+10)-rdat%r03(l,i+14)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      v632e2(i,j)=+rdat%r04(k,i+29)-rdat%r04(k,i+33)+rdat%r02(m,i+20)-rdat%r02(m,i+24)
            enddo

            if(jj > 1) cycle
            do i=1,5
      v632e1(i  )=+rdat%r03(k,i+42)-rdat%r03(k,i+47)+rdat%r01(l,i+20)-rdat%r01(l,i+25)
            enddo
         enddo
      enddo

      j1=2
      do j=1,15
         e5(  j,3,2)= t632e5(  j)+u6e5(  j,j1)-v632e5(  j)
      enddo
      do j=1,10
         e4(1,j,3,2)= t632e4(1,j)+u6e4(1,j,j1)-v632e4(1,j)
         e4(2,j,3,2)= t632e4(2,j)+u6e4(2,j,j1)-v632e4(2,j)
      enddo
      do j=1,6
         e3(1,j,3,2)= t632e3(1,j)+u6e3(1,j,j1)-v632e3(1,j)
         e3(2,j,3,2)= t632e3(2,j)+u6e3(2,j,j1)-v632e3(2,j)
         e3(3,j,3,2)= t632e3(3,j)+u6e3(3,j,j1)-v632e3(3,j)
         e3(4,j,3,2)= t632e3(4,j)+u6e3(4,j,j1)-v632e3(4,j)
      enddo
      do j=1,3
         e2(1,j,3,2)= t632e2(1,j)+u6e2(1,j,j1)-v632e2(1,j)
         e2(2,j,3,2)= t632e2(2,j)+u6e2(2,j,j1)-v632e2(2,j)
         e2(3,j,3,2)= t632e2(3,j)+u6e2(3,j,j1)-v632e2(3,j)
         e2(4,j,3,2)= t632e2(4,j)+u6e2(4,j,j1)-v632e2(4,j)
      enddo
      do i=1,5
         e1(i  ,3,2)= t632e1(i  )+u6e1(i  ,j1)-v632e1(i  )
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,4,2)
            l= k-3-jj-jj

      e5(  j,4,2)=+rdat%r08(k  )+rdat%r06(l,1)* 3

            if(jj > 4) cycle
      e4(1,j,4,2)=+rdat%r07(k,3)+rdat%r05(l,5)* 3
      e4(2,j,4,2)=+rdat%r07(k,4)+rdat%r05(l,6)* 3

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,4,2)=+rdat%r06(k,i+ 8)+rdat%r04(l,i+13)* 3
            enddo

            if(jj > 2) cycle
            do i=1,4
      e2(i,j,4,2)=+rdat%r05(k,i+20)+rdat%r03(l,i+26)* 3
            enddo

            if(jj > 1) cycle
            m= in6(l)
            do i=1,5
      e1(i  ,4,2)=+rdat%r04(k,i+37)+rdat%r02(m,i+36)* 3
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,5,2)
            l= k-3-jj
            m= l-jj

      e5(  j,5,2)=+rdat%r08(k  )-rdat%r07(l,1)+rdat%r06(m,1)-rdat%r05(j,1)

            if(jj > 4) cycle
      e4(1,j,5,2)=+rdat%r07(k,3)-rdat%r06(l,5)+rdat%r05(m,5)-rdat%r04(j,6)
      e4(2,j,5,2)=+rdat%r07(k,4)-rdat%r06(l,6)+rdat%r05(m,6)-rdat%r04(j,7)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,5,2)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+12)+rdat%r04(m,i+13)-rdat%r03(j,i+10)
            enddo

            if(jj > 2) cycle
            n= in6(j)
            do i=1,4
      e2(i,j,5,2)=+rdat%r05(k,i+20)-rdat%r04(l,i+29)+rdat%r03(m,i+26)-rdat%r02(n,i+20)
            enddo

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
      e1(i  ,5,2)=+rdat%r04(k,i+37)-rdat%r03(l,i+42)+rdat%r02(m,i+36)-rdat%r01(j,i+20)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,6,2)
            l= k-4-jj
            m= l-1-jj
            n= m-2-jj

      e5(  j,6,2)=+rdat%r08(k  )-rdat%r07(l,1)+rdat%r06(m,1)* 3 -rdat%r05(n,1)* 3

            if(jj > 4) cycle
      e4(1,j,6,2)=+rdat%r07(k,3)-rdat%r06(l,5)+rdat%r05(m,5)* 3 -rdat%r04(n,6)* 3
      e4(2,j,6,2)=+rdat%r07(k,4)-rdat%r06(l,6)+rdat%r05(m,6)* 3 -rdat%r04(n,7)* 3

            if(jj > 3) cycle
            do i=1,4
            i4= i+10
      e3(i,j,6,2)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+12)+rdat%r04(m,i+13)* 3 -rdat%r03(n,i4)* 3
            enddo

            if(jj > 2) cycle
            n= in6(n)
            do i=1,4
            i4= i+20
      e2(i,j,6,2)=+rdat%r05(k,i+20)-rdat%r04(l,i+29)+rdat%r03(m,i+26)* 3 -rdat%r02(n,i4)* 3
            enddo
            n= m-2-jj

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
            i4= i+20
      e1(i  ,6,2)=+rdat%r04(k,i+37)-rdat%r03(l,i+42)+rdat%r02(m,i+36)* 3 -rdat%r01(n,i4)* 3
            enddo
         enddo
      enddo
      j1=4
      do j=1,15
         e5(  j,1,3)= t531e5(  j)+u5e5(  j,j1)+v531e5(  j)
      enddo
      do j=1,10
         e4(1,j,1,3)= t531e4(1,j)+u5e4(1,j,j1)+v531e4(1,j)
         e4(2,j,1,3)= t531e4(2,j)+u5e4(2,j,j1)+v531e4(2,j)
      enddo
      do j=1,6
         e3(1,j,1,3)= t531e3(1,j)+u5e3(1,j,j1)+v531e3(1,j)
         e3(2,j,1,3)= t531e3(2,j)+u5e3(2,j,j1)+v531e3(2,j)
         e3(3,j,1,3)= t531e3(3,j)+u5e3(3,j,j1)+v531e3(3,j)
         e3(4,j,1,3)= t531e3(4,j)+u5e3(4,j,j1)+v531e3(4,j)
      enddo
      do j=1,3
         e2(1,j,1,3)= t531e2(1,j)+u5e2(1,j,j1)+v531e2(1,j)
         e2(2,j,1,3)= t531e2(2,j)+u5e2(2,j,j1)+v531e2(2,j)
         e2(3,j,1,3)= t531e2(3,j)+u5e2(3,j,j1)+v531e2(3,j)
         e2(4,j,1,3)= t531e2(4,j)+u5e2(4,j,j1)+v531e2(4,j)
      enddo
      do i=1,5
         e1(i  ,1,3)= t531e1(i  )+u5e1(i  ,j1)+v531e1(i  )
      enddo
      j1=4
      do j=1,15
         e5(  j,2,3)= t632e5(  j)+u6e5(  j,j1)+v632e5(  j)
      enddo
      do j=1,10
         e4(1,j,2,3)= t632e4(1,j)+u6e4(1,j,j1)+v632e4(1,j)
         e4(2,j,2,3)= t632e4(2,j)+u6e4(2,j,j1)+v632e4(2,j)
      enddo
      do j=1,6
         e3(1,j,2,3)= t632e3(1,j)+u6e3(1,j,j1)+v632e3(1,j)
         e3(2,j,2,3)= t632e3(2,j)+u6e3(2,j,j1)+v632e3(2,j)
         e3(3,j,2,3)= t632e3(3,j)+u6e3(3,j,j1)+v632e3(3,j)
         e3(4,j,2,3)= t632e3(4,j)+u6e3(4,j,j1)+v632e3(4,j)
      enddo
      do j=1,3
         e2(1,j,2,3)= t632e2(1,j)+u6e2(1,j,j1)+v632e2(1,j)
         e2(2,j,2,3)= t632e2(2,j)+u6e2(2,j,j1)+v632e2(2,j)
         e2(3,j,2,3)= t632e2(3,j)+u6e2(3,j,j1)+v632e2(3,j)
         e2(4,j,2,3)= t632e2(4,j)+u6e2(4,j,j1)+v632e2(4,j)
      enddo
      do i=1,5
         e1(i  ,2,3)= t632e1(i  )+u6e1(i  ,j1)+v632e1(i  )
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,3,3)
            l= k-4-jj
            m= l-3-jj
            n= m-2-jj

      e5(  j,3,3)=    +rdat%r08(k  )    -rdat%r07(l,1)* 2 -rdat%r07(l,2)* 2               &
     &    +rdat%r06(m,1)* 6 +rdat%r06(m,2)    +rdat%r06(m,3)* 4 +rdat%r06(m,4)+                 &
     &   (-rdat%r05(n,1)* 3 -rdat%r05(n,2)* 3 -rdat%r05(n,3)    -rdat%r05(n,4))* 2              &
     &    +rdat%r04(j,1)* 3 +rdat%r04(j,2)    +rdat%r04(j,3)* 4 +rdat%r04(j,4)+rdat%r04(j,5)

            if(jj > 4) cycle
            do i=1,2
               i4= i+ 4
               i5= i+ 5
      e4(i,j,3,3)=    +rdat%r07(k,i+ 2)    -rdat%r06(l,i+ 4)* 2 -rdat%r06(l,i+ 6)* 2      &
     &   +rdat%r05(m,i4)* 6 +rdat%r05(m,i4+2)    +rdat%r05(m,i4+4)* 4 +rdat%r05(m,i4+6)+        &
     &  (-rdat%r04(n,i5)* 3 -rdat%r04(n,i5+2)* 3 -rdat%r04(n,i5+4)    -rdat%r04(n,i5+6))* 2     &
     &   +rdat%r03(j,i   )* 3 +rdat%r03(j,i+ 2)                                     &
     &   +rdat%r03(j,i+ 4)* 4 +rdat%r03(j,i+ 6)+rdat%r03(j,i+ 8)
            enddo

            if(jj > 3) cycle
            jn= in6(j)
            do i=1,4
               i4= i+13
               i5= i+10
      e3(i,j,3,3)=    +rdat%r06(k,i+ 8)    -rdat%r05(l,i+12)* 2 -rdat%r05(l,i+16)* 2      &
     &   +rdat%r04(m,i4)* 6 +rdat%r04(m,i4+4)    +rdat%r04(m,i4+8)* 4 +rdat%r04(m,i4+12)+       &
     &  (-rdat%r03(n,i5)* 3 -rdat%r03(n,i5+4)* 3 -rdat%r03(n,i5+8)    -rdat%r03(n,i5+12))* 2    &
     &   +rdat%r02(jn,i   )* 3 +rdat%r02(jn,i+ 4)                                   &
     &   +rdat%r02(jn,i+ 8)* 4 +rdat%r02(jn,i+12)+rdat%r02(jn,i+16)
            enddo

            if(jj > 2) cycle
            n= in6(n)
            do i=1,4
               i4= i+26
               i5= i+20
      e2(i,j,3,3)=    +rdat%r05(k,i+20)    -rdat%r04(l,i+29)* 2 -rdat%r04(l,i+33)* 2      &
     &   +rdat%r03(m,i4)* 6 +rdat%r03(m,i4+4)    +rdat%r03(m,i4+8)* 4 +rdat%r03(m,i4+12)+       &
     &  (-rdat%r02(n,i5)* 3 -rdat%r02(n,i5+4)* 3 -rdat%r02(n,i5+8)    -rdat%r02(n,i5+12))* 2    &
     &   +rdat%r01(j,i   )* 3 +rdat%r01(j,i+ 4)                                     &
     &   +rdat%r01(j,i+ 8)* 4 +rdat%r01(j,i+12)+rdat%r01(j,i+16)
            enddo
            n= m-2-jj

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
               i4= i+36
               i5= i+20
      e1(i  ,3,3)=    +rdat%r04(k,i+37)    -rdat%r03(l,i+42)* 2 -rdat%r03(l,i+47)* 2      &
     &   +rdat%r02(m,i4)* 6 +rdat%r02(m,i4+5)    +rdat%r02(m,i4+10)* 4 +rdat%r02(m,i4+15)+      &
     &  (-rdat%r01(n,i5)* 3 -rdat%r01(n,i5+5)* 3 -rdat%r01(n,i5+10)    -rdat%r01(n,i5+15))* 2   &
     &   +rdat%r00(i,1)* 3 +rdat%r00(i,2)      +rdat%r00(i,3)* 4    +rdat%r00(i,4)+rdat%r00(i,5)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,4,3)
            l= k-3-jj
            m= l-2-jj

      e5(  j,4,3)=+rdat%r08(k  )-rdat%r07(l,2)* 2 +rdat%r06(m,1)+rdat%r06(m,4)

            if(jj > 4) cycle
      e4(1,j,4,3)=+rdat%r07(k,3)-rdat%r06(l,7)* 2 +rdat%r05(m,5)+rdat%r05(m,11)
      e4(2,j,4,3)=+rdat%r07(k,4)-rdat%r06(l,8)* 2 +rdat%r05(m,6)+rdat%r05(m,12)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,4,3)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+16)* 2 +rdat%r04(m,i+13)+rdat%r04(m,i+25)
            enddo

            if(jj > 2) cycle
            do i=1,4
      e2(i,j,4,3)=+rdat%r05(k,i+20)-rdat%r04(l,i+33)* 2 +rdat%r03(m,i+26)+rdat%r03(m,i+38)
            enddo

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
      e1(i  ,4,3)=+rdat%r04(k,i+37)-rdat%r03(l,i+47)* 2 +rdat%r02(m,i+36)+rdat%r02(m,i+51)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,5,3)
            l= k-3-jj
            m= l-2-jj

      e5(  j,5,3)=+rdat%r08(k  )-rdat%r07(l,1)    -rdat%r07(l,2)* 2                       &
     &                    +rdat%r06(m,1)* 3 +rdat%r06(m,3)* 2 +rdat%r06(m,4)              &
     &                    -rdat%r05(j,1)    -rdat%r05(j,2)* 2 -rdat%r05(j,4)

            if(jj > 4) cycle
            do i=1,2
      e4(i,j,5,3)=+rdat%r07(k,i+ 2)-rdat%r06(l,i+ 4)    -rdat%r06(l,i+ 6)* 2              &
     &                       +rdat%r05(m,i+ 4)* 3 +rdat%r05(m,i+ 8)* 2 +rdat%r05(m,i+10)  &
     &                       -rdat%r04(j,i+ 5)    -rdat%r04(j,i+ 7)* 2 -rdat%r04(j,i+11)
            enddo

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,5,3)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+12)    -rdat%r05(l,i+16)* 2              &
     &                       +rdat%r04(m,i+13)* 3 +rdat%r04(m,i+21)* 2 +rdat%r04(m,i+25)  &
     &                       -rdat%r03(j,i+10)    -rdat%r03(j,i+14)* 2 -rdat%r03(j,i+22)
            enddo

            if(jj > 2) cycle
            n= in6(j)
            do i=1,4
      e2(i,j,5,3)=+rdat%r05(k,i+20)-rdat%r04(l,i+29)    -rdat%r04(l,i+33)* 2              &
     &                       +rdat%r03(m,i+26)* 3 +rdat%r03(m,i+34)* 2 +rdat%r03(m,i+38)  &
     &                       -rdat%r02(n,i+20)    -rdat%r02(n,i+24)* 2 -rdat%r02(n,i+32)
            enddo

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
      e1(i  ,5,3)=+rdat%r04(k,i+37)-rdat%r03(l,i+42)    -rdat%r03(l,i+47)* 2              &
     &                       +rdat%r02(m,i+36)* 3 +rdat%r02(m,i+46)* 2 +rdat%r02(m,i+51)  &
     &                       -rdat%r01(j,i+20)    -rdat%r01(j,i+25)* 2 -rdat%r01(j,i+35)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,6,3)
            l= k-4-jj
            m= l-3-jj
            n= m-2-jj

            if(jj < 5) then
      e5(  j,6,3)= e5(j+jj,5,3)
            else
      e5(  j,6,3)=+rdat%r08(k  )-rdat%r07(l,1)    -rdat%r07(l,2)* 2                       &
     &                    +rdat%r06(m,1)* 3 +rdat%r06(m,3)* 2 +rdat%r06(m,4)              &
     &                    -rdat%r05(n,1)    -rdat%r05(n,2)* 2 -rdat%r05(n,4)
            endif

            if(jj > 4) cycle
            if(jj < 4) then
      e4(1,j,6,3)= e4(1,j+jj,5,3)
      e4(2,j,6,3)= e4(2,j+jj,5,3)
            else
               do i=1,2
      e4(i,j,6,3)=+rdat%r07(k,i+ 2)-rdat%r06(l,i+ 4)    -rdat%r06(l,i+ 6)* 2              &
     &                       +rdat%r05(m,i+ 4)* 3 +rdat%r05(m,i+ 8)* 2 +rdat%r05(m,i+10)  &
     &                       -rdat%r04(n,i+ 5)    -rdat%r04(n,i+ 7)* 2 -rdat%r04(n,i+11)
               enddo
            endif

            if(jj > 3) cycle
            if(jj < 3) then
               do i=1,4
      e3(i,j,6,3)= e3(i,j+jj,5,3)
               enddo
            else
               do i=1,4
      e3(i,j,6,3)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+12)    -rdat%r05(l,i+16)* 2              &
     &                       +rdat%r04(m,i+13)* 3 +rdat%r04(m,i+21)* 2 +rdat%r04(m,i+25)  &
     &                       -rdat%r03(n,i+10)    -rdat%r03(n,i+14)* 2 -rdat%r03(n,i+22)
               enddo
            endif

            if(jj > 2) cycle
            if(jj < 2) then
               do i=1,4
      e2(i,j,6,3)= e2(i,j+jj,5,3)
               enddo
            else
               n= in6(n)
               do i=1,4
      e2(i,j,6,3)=+rdat%r05(k,i+20)-rdat%r04(l,i+29)    -rdat%r04(l,i+33)* 2              &
     &                       +rdat%r03(m,i+26)* 3 +rdat%r03(m,i+34)* 2 +rdat%r03(m,i+38)  &
     &                       -rdat%r02(n,i+20)    -rdat%r02(n,i+24)* 2 -rdat%r02(n,i+32)
               enddo
               n= m-2-jj
            endif

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
      e1(i  ,6,3)=+rdat%r04(k,i+37)-rdat%r03(l,i+42)    -rdat%r03(l,i+47)* 2              &
     &                       +rdat%r02(m,i+36)* 3 +rdat%r02(m,i+46)* 2 +rdat%r02(m,i+51)  &
     &                       -rdat%r01(n,i+20)    -rdat%r01(n,i+25)* 2 -rdat%r01(n,i+35)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,3,4)
            l= k-3-jj
            m= l-2-jj

      e5(  j,3,4)=+rdat%r08(k  )-rdat%r07(l,1)* 2 +rdat%r06(m,1)+rdat%r06(m,2)

            if(jj > 4) cycle
      e4(1,j,3,4)=+rdat%r07(k,3)-rdat%r06(l,5)* 2 +rdat%r05(m,5)+rdat%r05(m,7)
      e4(2,j,3,4)=+rdat%r07(k,4)-rdat%r06(l,6)* 2 +rdat%r05(m,6)+rdat%r05(m,8)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,3,4)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+12)* 2 +rdat%r04(m,i+13)+rdat%r04(m,i+17)
            enddo

            if(jj > 2) cycle
            do i=1,4
      e2(i,j,3,4)=+rdat%r05(k,i+20)-rdat%r04(l,i+29)* 2 +rdat%r03(m,i+26)+rdat%r03(m,i+30)
            enddo

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
      e1(i  ,3,4)=+rdat%r04(k,i+37)-rdat%r03(l,i+42)* 2 +rdat%r02(m,i+36)+rdat%r02(m,i+41)
            enddo
         enddo
      enddo
      do j=1,15
         k= ind(j,1,5)
         e5(  j,1,5)=+rdat%r08(   k)-rdat%r07( j,2)+(rdat%r06( k,1)-rdat%r05( j,2))* 3
      enddo

      do j=1,10
         k= ind(j,1,5)
         e4(1,j,1,5)=+rdat%r07(k, 3)-rdat%r06(j, 7)+(rdat%r05(k, 5)-rdat%r04(j, 8))* 3
         e4(2,j,1,5)=+rdat%r07(k, 4)-rdat%r06(j, 8)+(rdat%r05(k, 6)-rdat%r04(j, 9))* 3
      enddo

      do j=1,6
         k= ind(j,1,5)
         e3(1,j,1,5)=+rdat%r06(k, 9)-rdat%r05(j,17)+(rdat%r04(k,14)-rdat%r03(j,15))* 3
         e3(2,j,1,5)=+rdat%r06(k,10)-rdat%r05(j,18)+(rdat%r04(k,15)-rdat%r03(j,16))* 3
         e3(3,j,1,5)=+rdat%r06(k,11)-rdat%r05(j,19)+(rdat%r04(k,16)-rdat%r03(j,17))* 3
         e3(4,j,1,5)=+rdat%r06(k,12)-rdat%r05(j,20)+(rdat%r04(k,17)-rdat%r03(j,18))* 3
      enddo

      do j=1,3
         k= ind(j,1,5)
         m= in6(j)
         e2(1,j,1,5)=+rdat%r05(k,21)-rdat%r04(j,34)+(rdat%r03(k,27)-rdat%r02(m,25))* 3
         e2(2,j,1,5)=+rdat%r05(k,22)-rdat%r04(j,35)+(rdat%r03(k,28)-rdat%r02(m,26))* 3
         e2(3,j,1,5)=+rdat%r05(k,23)-rdat%r04(j,36)+(rdat%r03(k,29)-rdat%r02(m,27))* 3
         e2(4,j,1,5)=+rdat%r05(k,24)-rdat%r04(j,37)+(rdat%r03(k,30)-rdat%r02(m,28))* 3
      enddo

         j=1
         k= ind(j,1,5)
         m= in6(k)
      do i=1,5
         e1(i  ,1,5)=+rdat%r04(k,i+37)-rdat%r03(j,i+47)+(rdat%r02(m,i+36)-rdat%r01(j,i+25))* 3
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,2,5)
            l= k-3-jj
            m= l-jj

      e5(  j,2,5)=+rdat%r08(k  )-rdat%r07(l,2)+rdat%r06(m,1)-rdat%r05(j,2)

            if(jj > 4) cycle
      e4(1,j,2,5)=+rdat%r07(k,3)-rdat%r06(l,7)+rdat%r05(m,5)-rdat%r04(j,8)
      e4(2,j,2,5)=+rdat%r07(k,4)-rdat%r06(l,8)+rdat%r05(m,6)-rdat%r04(j,9)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,2,5)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+16)+rdat%r04(m,i+13)-rdat%r03(j,i+14)
            enddo

            if(jj > 2) cycle
            n= in6(j)
            do i=1,4
      e2(i,j,2,5)=+rdat%r05(k,i+20)-rdat%r04(l,i+33)+rdat%r03(m,i+26)-rdat%r02(n,i+24)
            enddo

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
      e1(i  ,2,5)=+rdat%r04(k,i+37)-rdat%r03(l,i+47)+rdat%r02(m,i+36)-rdat%r01(j,i+25)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,3,5)
            l= k-3-jj
            m= l-2-jj

      e5(  j,3,5)=+rdat%r08(k  )-rdat%r07(l,1)* 2 -rdat%r07(l,2)                          &
     &                    +rdat%r06(m,1)* 3 +rdat%r06(m,2)+rdat%r06(m,3)* 2               &
     &                    -rdat%r05(j,1)* 2 -rdat%r05(j,2)-rdat%r05(j,3)

            if(jj > 4) cycle
            do i=1,2
      e4(i,j,3,5)=+rdat%r07(k,i+ 2)-rdat%r06(l,i+ 4)* 2 -rdat%r06(l,i+ 6)                 &
     &                       +rdat%r05(m,i+ 4)* 3 +rdat%r05(m,i+ 6)+rdat%r05(m,i+ 8)* 2   &
     &                       -rdat%r04(j,i+ 5)* 2 -rdat%r04(j,i+ 7)-rdat%r04(j,i+ 9)
            enddo

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,3,5)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+12)* 2 -rdat%r05(l,i+16)                 &
     &                       +rdat%r04(m,i+13)* 3 +rdat%r04(m,i+17)+rdat%r04(m,i+21)* 2   &
     &                       -rdat%r03(j,i+10)* 2 -rdat%r03(j,i+14)-rdat%r03(j,i+18)
            enddo

            if(jj > 2) cycle
            n= in6(j)
            do i=1,4
      e2(i,j,3,5)=+rdat%r05(k,i+20)-rdat%r04(l,i+29)* 2 -rdat%r04(l,i+33)                 &
     &                       +rdat%r03(m,i+26)* 3 +rdat%r03(m,i+30)+rdat%r03(m,i+34)* 2   &
     &                       -rdat%r02(n,i+20)* 2 -rdat%r02(n,i+24)-rdat%r02(n,i+28)
            enddo

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
      e1(i  ,3,5)=+rdat%r04(k,i+37)-rdat%r03(l,i+42)* 2 -rdat%r03(l,i+47)                 &
     &                       +rdat%r02(m,i+36)* 3 +rdat%r02(m,i+41)+rdat%r02(m,i+46)* 2   &
     &                       -rdat%r01(j,i+20)* 2 -rdat%r01(j,i+25)-rdat%r01(j,i+30)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,4,5)
            l= k-2-jj

      e5(  j,4,5)=+rdat%r08(k  )-rdat%r07(l,2)+rdat%r06(k,1)-rdat%r05(l,2)

            if(jj > 4) cycle
      e4(1,j,4,5)=+rdat%r07(k,3)-rdat%r06(l,7)+rdat%r05(k,5)-rdat%r04(l,8)
      e4(2,j,4,5)=+rdat%r07(k,4)-rdat%r06(l,8)+rdat%r05(k,6)-rdat%r04(l,9)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,4,5)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+16)+rdat%r04(k,i+13)-rdat%r03(l,i+14)
            enddo

            if(jj > 2) cycle
            m= in6(l)
            do i=1,4
      e2(i,j,4,5)=+rdat%r05(k,i+20)-rdat%r04(l,i+33)+rdat%r03(k,i+26)-rdat%r02(m,i+24)
            enddo

            if(jj > 1) cycle
            m= in6(k)
            do i=1,5
      e1(i  ,4,5)=+rdat%r04(k,i+37)-rdat%r03(l,i+47)+rdat%r02(m,i+36)-rdat%r01(l,i+25)
            enddo
         enddo
      enddo

      j1=3
      do j=1,15
         e5(  j,5,5)= t531e5(  j)+u5e5(  j,j1)
      enddo
      do j=1,10
         e4(1,j,5,5)= t531e4(1,j)+u5e4(1,j,j1)
         e4(2,j,5,5)= t531e4(2,j)+u5e4(2,j,j1)
      enddo
      do j=1,6
         e3(1,j,5,5)= t531e3(1,j)+u5e3(1,j,j1)
         e3(2,j,5,5)= t531e3(2,j)+u5e3(2,j,j1)
         e3(3,j,5,5)= t531e3(3,j)+u5e3(3,j,j1)
         e3(4,j,5,5)= t531e3(4,j)+u5e3(4,j,j1)
      enddo
      do j=1,3
         e2(1,j,5,5)= t531e2(1,j)+u5e2(1,j,j1)
         e2(2,j,5,5)= t531e2(2,j)+u5e2(2,j,j1)
         e2(3,j,5,5)= t531e2(3,j)+u5e2(3,j,j1)
         e2(4,j,5,5)= t531e2(4,j)+u5e2(4,j,j1)
      enddo
      do i=1,5
         e1(i  ,5,5)= t531e1(i  )+u5e1(i  ,j1)
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,6,5)
            l= k-3-jj
            m= l-2-jj

      e5(  j,6,5)=+rdat%r08(k  )-rdat%r07(l,1)-rdat%r07(l,2)+rdat%r06(m,1)+rdat%r06(m,3)

            if(jj > 4) cycle
      e4(1,j,6,5)=+rdat%r07(k,3)-rdat%r06(l,5)-rdat%r06(l,7)+rdat%r05(m,5)+rdat%r05(m, 9)
      e4(2,j,6,5)=+rdat%r07(k,4)-rdat%r06(l,6)-rdat%r06(l,8)+rdat%r05(m,6)+rdat%r05(m,10)

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,6,5)=+rdat%r06(k,i+ 8)-rdat%r05(l,i+12)-rdat%r05(l,i+16)                     &
     &                       +rdat%r04(m,i+13)+rdat%r04(m,i+21)
            enddo

            if(jj > 2) cycle
            do i=1,4
      e2(i,j,6,5)=+rdat%r05(k,i+20)-rdat%r04(l,i+29)-rdat%r04(l,i+33)                     &
     &                       +rdat%r03(m,i+26)+rdat%r03(m,i+34)
            enddo

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
      e1(i  ,6,5)=+rdat%r04(k,i+37)-rdat%r03(l,i+42)-rdat%r03(l,i+47)                     &
     &                       +rdat%r02(m,i+36)+rdat%r02(m,i+46)
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,2,6)
            l= k-4-jj
            m= l-1-jj
            n= m-2-jj

      e5(  j,2,6)=+rdat%r08(k  )-rdat%r07(l,2)+rdat%r06(m,1)* 3 -rdat%r05(n,2)* 3

            if(jj > 4) cycle
      e4(1,j,2,6)=+rdat%r07(k,3)-rdat%r06(l,7)+rdat%r05(m,5)* 3 -rdat%r04(n,8)* 3
      e4(2,j,2,6)=+rdat%r07(k,4)-rdat%r06(l,8)+rdat%r05(m,6)* 3 -rdat%r04(n,9)* 3

            if(jj > 3) cycle
            do i=1,4
      e3(i,j,2,6)=+rdat%r06(k,i+8)-rdat%r05(l,i+16)+rdat%r04(m,i+13)* 3 -rdat%r03(n,i+14)* 3
            enddo

            if(jj > 2) cycle
            n= in6(n)
            do i=1,4
            i4= i+20
      e2(i,j,2,6)=+rdat%r05(k,i4 )-rdat%r04(l,i+33)+rdat%r03(m,i+26)* 3 -rdat%r02(n,i+24)* 3
            enddo
            n= m-2-jj

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
            i4= i+37
      e1(i  ,2,6)=+rdat%r04(k,i4 )-rdat%r03(l,i+47)+rdat%r02(m,i+36)* 3 -rdat%r01(n,i+25)* 3
            enddo
         enddo
      enddo
      j= 0
      do jj=1,5
         do ii=1,jj
            j= j+1
            k= ind(j,3,6)
            l= k-4-jj
            m= l-3-jj
            n= m-2-jj

            if(jj < 5) then
      e5(  j,3,6)= e5(j+jj,3,5)
            else
      e5(  j,3,6)=+rdat%r08(k  )    -rdat%r07(l,1)* 2 -rdat%r07(l,2)                      &
     &            +rdat%r06(m,1)* 3 +rdat%r06(m,2)    +rdat%r06(m,3)* 2                   &
     &            -rdat%r05(n,1)* 2 -rdat%r05(n,2)    -rdat%r05(n,3)
            endif

            if(jj > 4) cycle
            if(jj < 4) then
      e4(1,j,3,6)= e4(1,j+jj,3,5)
      e4(2,j,3,6)= e4(2,j+jj,3,5)
            else
               do i=1,2
      e4(i,j,3,6)=+rdat%r07(k,i+ 2)    -rdat%r06(l,i+ 4)* 2 -rdat%r06(l,i+ 6)             &
     &            +rdat%r05(m,i+ 4)* 3 +rdat%r05(m,i+ 6)    +rdat%r05(m,i+ 8)* 2          &
     &            -rdat%r04(n,i+ 5)* 2 -rdat%r04(n,i+ 7)    -rdat%r04(n,i+ 9)
               enddo
            endif

            if(jj > 3) cycle
            if(jj < 3) then
               do i=1,4
      e3(i,j,3,6)= e3(i,j+jj,3,5)
               enddo
            else
               do i=1,4
               i4= i+10
      e3(i,j,3,6)=+rdat%r06(k,i+ 8)    -rdat%r05(l,i+12)* 2 -rdat%r05(l,i+16)             &
     &            +rdat%r04(m,i+13)* 3 +rdat%r04(m,i+17)    +rdat%r04(m,i+21)* 2          &
     &            -rdat%r03(n,i+10)* 2 -rdat%r03(n,i+14)    -rdat%r03(n,i+18)
               enddo
            endif

            if(jj > 2) cycle
            if(jj < 2) then
               do i=1,4
      e2(i,j,3,6)= e2(i,j+jj,3,5)
               enddo
            else
               n= in6(n)
               do i=1,4
      e2(i,j,3,6)=+rdat%r05(k,i+20)    -rdat%r04(l,i+29)* 2 -rdat%r04(l,i+33)             &
     &            +rdat%r03(m,i+26)* 3 +rdat%r03(m,i+30)    +rdat%r03(m,i+34)* 2          &
     &            -rdat%r02(n,i+20)* 2 -rdat%r02(n,i+24)    -rdat%r02(n,i+28)
               enddo
               n= m-2-jj
            endif

            if(jj > 1) cycle
            m= in6(m)
            do i=1,5
      e1(i  ,3,6)=+rdat%r04(k,i+37)    -rdat%r03(l,i+42)* 2 -rdat%r03(l,i+47)             &
     &            +rdat%r02(m,i+36)* 3 +rdat%r02(m,i+41)    +rdat%r02(m,i+46)* 2          &
     &            -rdat%r01(n,i+20)* 2 -rdat%r01(n,i+25)    -rdat%r01(n,i+30)
            enddo
         enddo
      enddo
      j1=3
      do j=1,15
         e5(  j,6,6)= t632e5(  j)+u6e5(  j,j1)
      enddo
      do j=1,10
         e4(1,j,6,6)= t632e4(1,j)+u6e4(1,j,j1)
         e4(2,j,6,6)= t632e4(2,j)+u6e4(2,j,j1)
      enddo
      do j=1,6
         e3(1,j,6,6)= t632e3(1,j)+u6e3(1,j,j1)
         e3(2,j,6,6)= t632e3(2,j)+u6e3(2,j,j1)
         e3(3,j,6,6)= t632e3(3,j)+u6e3(3,j,j1)
         e3(4,j,6,6)= t632e3(4,j)+u6e3(4,j,j1)
      enddo
      do j=1,3
         e2(1,j,6,6)= t632e2(1,j)+u6e2(1,j,j1)
         e2(2,j,6,6)= t632e2(2,j)+u6e2(2,j,j1)
         e2(3,j,6,6)= t632e2(3,j)+u6e2(3,j,j1)
         e2(4,j,6,6)= t632e2(4,j)+u6e2(4,j,j1)
      enddo
      do i=1,5
         e1(i  ,6,6)= t632e1(i  )+u6e1(i  ,j1)
      enddo

      xx= qx*qx
      zz= qz*qz
      xz= qx*qz
      xxx= xx*qx
      xxz= xx*qz
      xzz= zz*qx
      zzz= zz*qz
      xxxx= xx*xx
      xxxz= xx*xz
      xxzz= xx*zz
      xzzz= xz*zz
      zzzz= zz*zz

      qxd= qx+qx
      qzd= qz+qz
      xzd= xz+xz
      xxxd= xxx+xxx
      xxzd= xxz+xxz
      xzzd= xzz+xzz
      zzzd= zzz+zzz

      xzq= xz+xz+xz+xz

      do l=1,lx
         do k=1,kx

            if(k == 1 .and. l == 2) cycle

            if(k /= 3 .and. l == 4) cycle

            if(k == 1 .and. l == 6) cycle
            if(k == 4 .and. l == 6) cycle
            if(k == 5 .and. l == 6) cycle

            f(1,1,k,l) = e5(  1,k,l)+e3(1,1,k,l)* 6 +e1(  1,k,l)* 3 +(+e4(1,1,k,l)+e4(2,1,k,l)+ (+e2(1,1,k,l)+e2(2,1,k,l))* 3 )*qxd +(+e3(2,1,k,l)+e3(3,1,k,l)* 4 +e3(4,1,k,l) +e1(  2,k,l)+e1(  3,k,l)* 4 +e1(  4,k,l))*xx +(+e2(3,1,k,l)+e2(4,1,k,l))*xxxd +e1(  5,k,l)*xxxx
            f(2,1,k,l) = e5(  4,k,l)+e3(1,4,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(2,4,k,l)+e2(2,1,k,l))*qxd +(+e3(4,4,k,l)+e1(  4,k,l))*xx
            f(3,1,k,l) = e5(  6,k,l)+e3(1,6,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(2,6,k,l)+e2(2,1,k,l))*qxd +(+e4(1,3,k,l)+e2(1,3,k,l))*qzd +(+e3(4,6,k,l)+e1(  4,k,l))*xx +e3(3,3,k,l)*xzq +(+e3(2,1,k,l)+e1(  2,k,l))*zz +e2(4,3,k,l)*xxzd +e2(3,1,k,l)*xzzd +e1(  5,k,l)*xxzz
            f(4,1,k,l) = e5(  2,k,l)+e3(1,2,k,l)* 3 +(+e4(1,2,k,l)+e4(2,2,k,l)* 2 +e2(1,2,k,l)+e2(2,2,k,l)* 2 )*qx +(+e3(3,2,k,l)* 2 +e3(4,2,k,l))*xx +e2(4,2,k,l)*xxx
            f(5,1,k,l) = e5(  3,k,l)+e3(1,3,k,l)* 3 +(+e4(1,3,k,l)+e4(2,3,k,l)* 2 +e2(1,3,k,l)+e2(2,3,k,l)* 2 )*qx +(+e4(1,1,k,l)+e2(1,1,k,l)* 3 )*qz +(+e3(3,3,k,l)* 2 +e3(4,3,k,l))*xx +(+e3(2,1,k,l)+e3(3,1,k,l)* 2 +e1(  2,k,l)+e1(  3,k,l)* 2 )*xz +e2(4,3,k,l)*xxx +(+e2(3,1,k,l)* 2 +e2(4,1,k,l))*xxz +e1(  5,k,l)*xxxz
            f(6,1,k,l) = e5(  5,k,l)+e3(1,5,k,l) +e4(2,5,k,l)*qxd +(+e4(1,2,k,l)+e2(1,2,k,l))*qz +e3(4,5,k,l)*xx +e3(3,2,k,l)*xzd +e2(4,2,k,l)*xxz

            f(1,2,k,l) = e5(  4,k,l)+e3(1,4,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(1,4,k,l)+e2(1,1,k,l))*qxd +(+e3(2,4,k,l)+e1(  2,k,l))*xx
            f(2,2,k,l) = e5( 11,k,l)+e3(1,4,k,l)* 6 +e1(  1,k,l)* 3
            f(3,2,k,l) = e5( 13,k,l)+e3(1,6,k,l)+e3(1,4,k,l)+e1(  1,k,l) +(+e4(1,8,k,l)+e2(1,3,k,l))*qzd +(+e3(2,4,k,l)+e1(  2,k,l))*zz
            f(4,2,k,l) = e5(  7,k,l)+e3(1,2,k,l)* 3 +(+e4(1,7,k,l)+e2(1,2,k,l)* 3 )*qx
            f(5,2,k,l) = e5(  8,k,l)+e3(1,3,k,l) +(+e4(1,8,k,l)+e2(1,3,k,l))*qx +(+e4(1,4,k,l)+e2(1,1,k,l))*qz +(+e3(2,4,k,l)+e1(  2,k,l))*xz
            f(6,2,k,l) = e5( 12,k,l)+e3(1,5,k,l)* 3 +(+e4(1,7,k,l)+e2(1,2,k,l)* 3 )*qz

            f(1,3,k,l) = e5(  6,k,l)+e3(1,6,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(1,6,k,l)+e2(1,1,k,l))*qxd +(+e4(2,3,k,l)+e2(2,3,k,l))*qzd +(+e3(2,6,k,l)+e1(  2,k,l))*xx +e3(3,3,k,l)*xzq +(+e3(4,1,k,l)+e1(  4,k,l))*zz +e2(3,3,k,l)*xxzd +e2(4,1,k,l)*xzzd +e1(  5,k,l)*xxzz
            f(2,3,k,l) = e5( 13,k,l)+e3(1,6,k,l)+e3(1,4,k,l)+e1(  1,k,l) +(+e4(2,8,k,l)+e2(2,3,k,l))*qzd +(+e3(4,4,k,l)+e1(  4,k,l))*zz
            f(3,3,k,l) = e5( 15,k,l)+e3(1,6,k,l)* 6 +e1(  1,k,l)* 3 +(+e4(1,10,k,l)+e4(2,10,k,l)+ (+e2(1,3,k,l)+e2(2,3,k,l))* 3 )*qzd +(+e3(2,6,k,l)+e3(3,6,k,l)* 4 +e3(4,6,k,l) +e1(  2,k,l)+e1(  3,k,l)* 4 +e1(  4,k,l))*zz +(+e2(3,3,k,l)+e2(4,3,k,l))*zzzd +e1(  5,k,l)*zzzz
            f(4,3,k,l) = e5(  9,k,l)+e3(1,2,k,l) +(+e4(1,9,k,l)+e2(1,2,k,l))*qx +e4(2,5,k,l)*qzd +e3(3,5,k,l)*xzd +e3(4,2,k,l)*zz +e2(4,2,k,l)*xzz
            f(5,3,k,l) = e5( 10,k,l)+e3(1,3,k,l)* 3 +(+e4(1,10,k,l)+e2(1,3,k,l)* 3 )*qx +(+e4(1,6,k,l)+e4(2,6,k,l)* 2 +e2(1,1,k,l)+e2(2,1,k,l)* 2 )*qz +(+e3(2,6,k,l)+e3(3,6,k,l)* 2 +e1(  2,k,l)+e1(  3,k,l)* 2 )*xz +(+e3(3,3,k,l)* 2 +e3(4,3,k,l))*zz +(+e2(3,3,k,l)* 2 +e2(4,3,k,l))*xzz +e2(4,1,k,l)*zzz +e1(  5,k,l)*xzzz
            f(6,3,k,l) = e5( 14,k,l)+e3(1,5,k,l)* 3 +(+e4(1,9,k,l)+e4(2,9,k,l)* 2 +e2(1,2,k,l)+e2(2,2,k,l)* 2 )*qz +(+e3(3,5,k,l)* 2 +e3(4,5,k,l))*zz +e2(4,2,k,l)*zzz

            f(1,4,k,l) = e5(  2,k,l)+e3(1,2,k,l)* 3 +(+e4(1,2,k,l)* 2 +e4(2,2,k,l) +e2(1,2,k,l)* 2 +e2(2,2,k,l))*qx +(+e3(2,2,k,l)+e3(3,2,k,l)* 2 )*xx +e2(3,2,k,l)*xxx
            f(2,4,k,l) = e5(  7,k,l)+e3(1,2,k,l)* 3 +(+e4(2,7,k,l)+e2(2,2,k,l)* 3 )*qx
            f(3,4,k,l) = e5(  9,k,l)+e3(1,2,k,l) +(+e4(2,9,k,l)+e2(2,2,k,l))*qx +e4(1,5,k,l)*qzd +e3(3,5,k,l)*xzd +e3(2,2,k,l)*zz +e2(3,2,k,l)*xzz
            f(4,4,k,l) = e5(  4,k,l)+e3(1,4,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(1,4,k,l)+e4(2,4,k,l) +e2(1,1,k,l)+e2(2,1,k,l))*qx +(+e3(3,4,k,l)+e1(  3,k,l))*xx
            f(5,4,k,l) = e5(  5,k,l)+e3(1,5,k,l) +(+e4(1,5,k,l)+e4(2,5,k,l))*qx +(+e4(1,2,k,l)+e2(1,2,k,l))*qz +e3(3,5,k,l)*xx +(+e3(2,2,k,l)+e3(3,2,k,l))*xz +e2(3,2,k,l)*xxz
            f(6,4,k,l) = e5(  8,k,l)+e3(1,3,k,l) +(+e4(2,8,k,l)+e2(2,3,k,l))*qx +(+e4(1,4,k,l)+e2(1,1,k,l))*qz +(+e3(3,4,k,l)+e1(  3,k,l))*xz

            f(1,5,k,l) = e5(  3,k,l)+e3(1,3,k,l)* 3 +(+e4(1,3,k,l)* 2 +e4(2,3,k,l) +e2(1,3,k,l)* 2 +e2(2,3,k,l))*qx +(+e4(2,1,k,l)+e2(2,1,k,l)* 3 )*qz +(+e3(2,3,k,l)+e3(3,3,k,l)* 2 )*xx +(+e3(3,1,k,l)* 2 +e3(4,1,k,l) +e1(  3,k,l)* 2 +e1(  4,k,l))*xz +e2(3,3,k,l)*xxx +(+e2(3,1,k,l)+e2(4,1,k,l)* 2 )*xxz +e1(  5,k,l)*xxxz
            f(2,5,k,l) = e5(  8,k,l)+e3(1,3,k,l) +(+e4(2,8,k,l)+e2(2,3,k,l))*qx +(+e4(2,4,k,l)+e2(2,1,k,l))*qz +(+e3(4,4,k,l)+e1(  4,k,l))*xz
            f(3,5,k,l) = e5( 10,k,l)+e3(1,3,k,l)* 3 +(+e4(2,10,k,l)+e2(2,3,k,l)* 3 )*qx +(+e4(1,6,k,l)* 2 +e4(2,6,k,l) +e2(1,1,k,l)* 2 +e2(2,1,k,l))*qz +(+e3(3,6,k,l)* 2 +e3(4,6,k,l) +e1(  3,k,l)* 2 +e1(  4,k,l))*xz +(+e3(2,3,k,l)+e3(3,3,k,l)* 2 )*zz +(+e2(3,3,k,l)+e2(4,3,k,l)* 2 )*xzz +e2(3,1,k,l)*zzz +e1(  5,k,l)*xzzz
            f(4,5,k,l) = e5(  5,k,l)+e3(1,5,k,l) +(+e4(1,5,k,l)+e4(2,5,k,l))*qx +(+e4(2,2,k,l)+e2(2,2,k,l))*qz +e3(3,5,k,l)*xx +(+e3(3,2,k,l)+e3(4,2,k,l))*xz +e2(4,2,k,l)*xxz
            f(5,5,k,l) = e5(  6,k,l)+e3(1,6,k,l)+e3(1,1,k,l)+e1(  1,k,l) +(+e4(1,6,k,l)+e4(2,6,k,l) +e2(1,1,k,l)+e2(2,1,k,l))*qx +(+e4(1,3,k,l)+e4(2,3,k,l) +e2(1,3,k,l)+e2(2,3,k,l))*qz +(+e3(3,6,k,l)+e1(  3,k,l))*xx +(+e3(2,3,k,l)+e3(3,3,k,l)* 2 +e3(4,3,k,l))*xz +(+e3(3,1,k,l)+e1(  3,k,l))*zz +(+e2(3,3,k,l)+e2(4,3,k,l))*xxz +(+e2(3,1,k,l)+e2(4,1,k,l))*xzz +e1(  5,k,l)*xxzz
            f(6,5,k,l) = e5(  9,k,l)+e3(1,2,k,l) +(+e4(2,9,k,l)+e2(2,2,k,l))*qx +(+e4(1,5,k,l)+e4(2,5,k,l))*qz +(+e3(3,5,k,l)+e3(4,5,k,l))*xz +e3(3,2,k,l)*zz +e2(4,2,k,l)*xzz

            f(1,6,k,l) = e5(  5,k,l)+e3(1,5,k,l) +e4(1,5,k,l)*qxd +(+e4(2,2,k,l)+e2(2,2,k,l))*qz +e3(2,5,k,l)*xx +e3(3,2,k,l)*xzd +e2(3,2,k,l)*xxz
            f(2,6,k,l) = e5( 12,k,l)+e3(1,5,k,l)* 3 +(+e4(2,7,k,l)+e2(2,2,k,l)* 3 )*qz
            f(3,6,k,l) = e5( 14,k,l)+e3(1,5,k,l)* 3 +(+e4(1,9,k,l)* 2 +e4(2,9,k,l) +e2(1,2,k,l)* 2 +e2(2,2,k,l))*qz +(+e3(2,5,k,l)+e3(3,5,k,l)* 2 )*zz +e2(3,2,k,l)*zzz
            f(4,6,k,l) = e5(  8,k,l)+e3(1,3,k,l) +(+e4(1,8,k,l)+e2(1,3,k,l))*qx +(+e4(2,4,k,l)+e2(2,1,k,l))*qz +(+e3(3,4,k,l)+e1(  3,k,l))*xz
            f(5,6,k,l) = e5(  9,k,l)+e3(1,2,k,l) +(+e4(1,9,k,l)+e2(1,2,k,l))*qx +(+e4(1,5,k,l)+e4(2,5,k,l))*qz +(+e3(2,5,k,l)+e3(3,5,k,l))*xz +e3(3,2,k,l)*zz +e2(3,2,k,l)*xzz
            f(6,6,k,l) = e5( 13,k,l)+e3(1,6,k,l)+e3(1,4,k,l)+e1(  1,k,l) +(+e4(1,8,k,l)+e4(2,8,k,l) +e2(1,3,k,l)+e2(2,3,k,l))*qz +(+e3(3,4,k,l)+e1(  3,k,l))*zz

         enddo
      enddo

      f(:,:,1,2)= f(:,:,2,1)

      f(:,:,1,4)= f(:,:,4,1)
      f(:,:,2,4)= f(:,:,4,2)
      f(:,:,4,4)= f(:,:,2,1)
      f(:,:,5,4)= f(:,:,6,1)
      f(:,:,6,4)= f(:,:,5,2)

      f(:,:,1,6)= f(:,:,4,5)
      f(:,:,4,6)= f(:,:,2,5)
      f(:,:,5,6)= f(:,:,6,5)

      end subroutine mcdv_21

! >
! >    @brief   auxiliary routine internal to rot.axis integrations
! >
! >    @details auxiliary routine internal to rot.axis integrations
! >
      subroutine fcufcc(rdat,n,xmdt,fcu,fcc)

      implicit none
      type(rotaxis_data_t) :: rdat

      integer :: n
      real(kind=dp) :: xmdt
      real(kind=dp) :: fcu(45,8),fcc(45,8)

      integer :: i, j, k
      real(kind=dp) :: xmdtx, xmdty, xmdtxy

      k= 1
      do i= 1,n
         k= k+i+1
         do j= 1, k
            fcu(j,i)= 1.0_dp
            fcc(j,i)= xmdt
         enddo
      enddo

      xmdty=-xmdt*rdat%acy
      xmdtx=-xmdt*rdat%aqx
      xmdtxy=xmdt*rdat%aqxy

         fcu( 1,1)=-rdat%aqx
         fcu( 2,1)=-rdat%acy

         fcc( 1,1)= xmdtx
         fcc( 2,1)= xmdty
! If(n<=1) return
         fcu( 4,2)= rdat%aqxy
         fcu( 5,2)=-rdat%aqx
         fcu( 6,2)=-rdat%acy

         fcc( 4,2)= xmdtxy
         fcc( 5,2)= xmdtx
         fcc( 6,2)= xmdty
! If(n<=2) return
         fcu( 1,3)=-rdat%aqx
         fcu( 2,3)=-rdat%acy
         fcu( 4,3)=-rdat%aqx
         fcu( 5,3)= rdat%aqxy
         fcu( 6,3)=-rdat%aqx
         fcu( 7,3)=-rdat%acy
         fcu( 9,3)=-rdat%acy

         fcc( 1,3)= xmdtx
         fcc( 2,3)= xmdty
         fcc( 4,3)= xmdtx
         fcc( 5,3)= xmdtxy
         fcc( 6,3)= xmdtx
         fcc( 7,3)= xmdty
         fcc( 9,3)= xmdty
! If(n<=3) return
         fcu( 2,4)= rdat%aqxy
         fcu( 3,4)=-rdat%aqx
         fcu( 5,4)=-rdat%acy
         fcu( 7,4)= rdat%aqxy
         fcu( 8,4)=-rdat%aqx
         fcu( 9,4)= rdat%aqxy
         fcu(10,4)=-rdat%aqx
         fcu(12,4)=-rdat%acy
         fcu(14,4)=-rdat%acy

         fcc( 2,4)= xmdtxy
         fcc( 3,4)= xmdtx
         fcc( 5,4)= xmdty
         fcc( 7,4)= xmdtxy
         fcc( 8,4)= xmdtx
         fcc( 9,4)= xmdtxy
         fcc(10,4)= xmdtx
         fcc(12,4)= xmdty
         fcc(14,4)= xmdty
! If(n<=4) return
         do j= 1,10
            fcu(j,5)= fcu(j,3)
            fcc(j,5)= fcc(j,3)
         enddo
         fcu(11,5)=-rdat%aqx
         fcu(12,5)= rdat%aqxy
         fcu(13,5)=-rdat%aqx
         fcu(14,5)= rdat%aqxy
         fcu(15,5)=-rdat%aqx
         fcu(16,5)=-rdat%acy
         fcu(18,5)=-rdat%acy
         fcu(20,5)=-rdat%acy

         fcc(11,5)= xmdtx
         fcc(12,5)= xmdtxy
         fcc(13,5)= xmdtx
         fcc(14,5)= xmdtxy
         fcc(15,5)= xmdtx
         fcc(16,5)= xmdty
         fcc(18,5)= xmdty
         fcc(20,5)= xmdty
! If(n<=5) return
         do j= 1,15
            fcu(j,6)= fcu(j,4)
            fcc(j,6)= fcc(j,4)
         enddo
         fcu(16,6)= rdat%aqxy
         fcu(17,6)=-rdat%aqx
         fcu(18,6)= rdat%aqxy
         fcu(19,6)=-rdat%aqx
         fcu(20,6)= rdat%aqxy
         fcu(21,6)=-rdat%aqx
         fcu(23,6)=-rdat%acy
         fcu(25,6)=-rdat%acy
         fcu(27,6)=-rdat%acy

         fcc(16,6)= xmdtxy
         fcc(17,6)= xmdtx
         fcc(18,6)= xmdtxy
         fcc(19,6)= xmdtx
         fcc(20,6)= xmdtxy
         fcc(21,6)= xmdtx
         fcc(23,6)= xmdty
         fcc(25,6)= xmdty
         fcc(27,6)= xmdty
      if(n<=6) return
         do j= 1,21
            fcu(j,7)= fcu(j,5)
            fcc(j,7)= fcc(j,5)
         enddo
         fcu(22,7)=-rdat%aqx
         fcu(23,7)= rdat%aqxy
         fcu(24,7)=-rdat%aqx
         fcu(25,7)= rdat%aqxy
         fcu(26,7)=-rdat%aqx
         fcu(27,7)= rdat%aqxy
         fcu(28,7)=-rdat%aqx
         fcu(29,7)=-rdat%acy
         fcu(31,7)=-rdat%acy
         fcu(33,7)=-rdat%acy
         fcu(35,7)=-rdat%acy

         fcc(22,7)= xmdtx
         fcc(23,7)= xmdtxy
         fcc(24,7)= xmdtx
         fcc(25,7)= xmdtxy
         fcc(26,7)= xmdtx
         fcc(27,7)= xmdtxy
         fcc(28,7)= xmdtx
         fcc(29,7)= xmdty
         fcc(31,7)= xmdty
         fcc(33,7)= xmdty
         fcc(35,7)= xmdty
      if(n<=7) return
         do j= 1,28
            fcu(j,8)= fcu(j,6)
            fcc(j,8)= fcc(j,6)
         enddo
         fcu(29,8)= rdat%aqxy
         fcu(30,8)=-rdat%aqx
         fcu(31,8)= rdat%aqxy
         fcu(32,8)=-rdat%aqx
         fcu(33,8)= rdat%aqxy
         fcu(34,8)=-rdat%aqx
         fcu(35,8)= rdat%aqxy
         fcu(36,8)=-rdat%aqx
         fcu(38,8)=-rdat%acy
         fcu(40,8)=-rdat%acy
         fcu(42,8)=-rdat%acy
         fcu(44,8)=-rdat%acy

         fcc(29,8)= xmdtxy
         fcc(30,8)= xmdtx
         fcc(31,8)= xmdtxy
         fcc(32,8)= xmdtx
         fcc(33,8)= xmdtxy
         fcc(34,8)= xmdtx
         fcc(35,8)= xmdtxy
         fcc(36,8)= xmdtx
         fcc(38,8)= xmdty
         fcc(40,8)= xmdty
         fcc(42,8)= xmdty
         fcc(44,8)= xmdty

      end subroutine fcufcc

! >
! >    @brief   auxiliary routine of order 6 for rot.axis integrations
! >
! >    @details auxiliary routine of order 6 for rot.axis integrations
! >
      subroutine frikr6(rdat,i1,i2,wrk,qd6,j0,qd5,k0,qd4,l0,qd3)

      implicit none
      type(rotaxis_data_t) :: rdat

      integer :: i1, i2, j0, k0, l0
      real(kind=dp) :: wrk(28,*), &
              qd6( 7,*), qd5( 6,*), qd4( 5,*), qd3( 4,*)

      integer :: i, j, k, l
      real(kind=dp) :: a11, b11, b13, b22, b23, b33, d11, f31l03, f31l15, &
            f41k03, f41k06, f41k15, f41k45, f42k03, f42k15, f43k03, &
            f43k06, f43k45, f51j03, f51j06, f51j10, f51j15, &
            f52j03, f52j06, f52j10, f53j03, f53j06, &
            f54j03, f54j10, f55j15, r11, &
            s13, s22, s23, s33, u11

      do i=i1,i2
         j=i+j0
         k=i+k0
         l=i+l0
         f31l03= qd3(1,l)* 3
         f31l15= qd3(1,l)*15

         f41k03= qd4(1,k)* 3
         f41k06= f41k03+f41k03
         f41k15= qd4(1,k)*15
         f41k45= f41k15  * 3
         f42k03= qd4(2,k)* 3
         f42k15= qd4(2,k)*15
         f43k03= qd4(3,k)* 3
         f43k06= f43k03+f43k03
         f43k45= f43k03  *15

         f51j03= qd5(1,j)* 3
         f51j06= f51j03+f51j03
         f51j10= qd5(1,j)*10
         f51j15= qd5(1,j)*15
         f52j03= qd5(2,j)* 3
         f52j06= f52j03+f52j03
         f52j10= qd5(2,j)*10
         f53j03= qd5(3,j)* 3
         f53j06= f53j03+f53j03
         f54j03= qd5(4,j)* 3
         f54j10= qd5(4,j)*10
         f55j15= qd5(5,j)*15

         a11   =     qd4(1,k)*rdat%aqx2-qd3(1,l)
         r11   =     f41k03  *rdat%acy2-f31l03

         b11   =     qd5(1,j)*rdat%aqx2-qd4(1,k)
         b13   =     qd5(1,j)*rdat%aqx2-f41k03
         b22   =     f52j03  *rdat%aqx2-f42k03
         b23   =     qd5(2,j)*rdat%aqx2-f42k03
         b33   =     qd5(3,j)*rdat%aqx2-qd4(3,k)
         s13   =     qd5(1,j)*rdat%acy2-f41k03
         s22   =     f52j03  *rdat%acy2-f42k03
         s23   =     f52j03  *rdat%acy2-f42k03  * 3
         s33   =     f53j06  *rdat%acy2-f43k06

         d11   =    (qd5(1,j)*rdat%aqx2-f41k06  )*rdat%aqx2+f31l03
         u11   =    (qd5(1,j)*rdat%acy2-f41k06  )*rdat%acy2+f31l03

         wrk( 1,i)=((qd6(1,i)*rdat%aqx2-f51j15  )*rdat%aqx2+f41k45)*rdat%aqx2-f31l15
         wrk( 2,i)= (qd6(1,i)*rdat%aqx2-f51j10  )*rdat%aqx2+f41k15
         wrk( 3,i)= (qd6(2,i)*rdat%aqx2-f52j10  )*rdat%aqx2+f42k15
         wrk( 4,i)=((qd6(1,i)*rdat%aqx2-f51j06  )*rdat%aqx2+f41k03)*rdat%acy2-d11
         wrk( 5,i)= (qd6(2,i)*rdat%aqx2-f52j06  )*rdat%aqx2+f42k03
         wrk( 6,i)= (qd6(3,i)*rdat%aqx2-f53j06  )*rdat%aqx2+f43k03      -d11
         wrk( 7,i)= (qd6(1,i)*rdat%aqx2-f51j03  )*rdat%acy2-b13* 3
         wrk( 8,i)= (qd6(2,i)*rdat%aqx2-f52j03  )*rdat%acy2-b23
         wrk( 9,i)=  qd6(3,i)*rdat%aqx2-f53j03        -b13
         wrk(10,i)=  qd6(4,i)*rdat%aqx2-f54j03        -b23* 3
         wrk(11,i)=((qd6(1,i)*rdat%acy2-f51j06  )*rdat%acy2+f41k03)*rdat%aqx2-u11
         wrk(12,i)= (qd6(2,i)*rdat%aqx2-qd5(2,j))*rdat%acy2-b22
         wrk(13,i)= (qd6(3,i)*rdat%aqx2-qd5(3,j)      -b11   )*rdat%acy2-b33+a11
         wrk(14,i)=  qd6(4,i)*rdat%aqx2-qd5(4,j)      -b22
         wrk(15,i)=  qd6(5,i)*rdat%aqx2-qd5(5,j)      -b33* 6  +a11* 3
         wrk(16,i)= (qd6(1,i)*rdat%acy2-f51j10  )*rdat%acy2+f41k15
         wrk(17,i)= (qd6(2,i)*rdat%acy2-f52j06  )*rdat%acy2+f42k03
         wrk(18,i)=  qd6(3,i)*rdat%acy2-f53j03        -s13
         wrk(19,i)=  qd6(4,i)*rdat%acy2-qd5(4,j)      -s22
         wrk(20,i)=  qd6(5,i)     -f53j06        +f41k03
         wrk(21,i)=  qd6(6,i)     -f54j10        +f42k15
         wrk(22,i)=((qd6(1,i)*rdat%acy2-f51j15  )*rdat%acy2+f41k45)*rdat%acy2-f31l15
         wrk(23,i)= (qd6(2,i)*rdat%acy2-f52j10  )*rdat%acy2+f42k15
         wrk(24,i)= (qd6(3,i)*rdat%acy2-f53j06  )*rdat%acy2+f43k03      -u11
         wrk(25,i)=  qd6(4,i)*rdat%acy2-f54j03        -s23
         wrk(26,i)=  qd6(5,i)*rdat%acy2-qd5(5,j)      -s33         +r11
         wrk(27,i)=  qd6(6,i)     -f54j10        +f42k15
         wrk(28,i)=  qd6(7,i)     -f55j15        +f43k45      -f31l15
      enddo

      end subroutine frikr6

! >
! >    @brief   auxiliary routine of order 7 for rot.axis integrations
! >
! >    @details auxiliary routine of order 7 for rot.axis integrations
! >
      subroutine frikr7(rdat,i1,i2,wrk,qd7,j0,qd6,k0,qd5,l0,qd4)

      implicit none
      type(rotaxis_data_t) :: rdat

      integer :: i1, i2, j0, k0, l0
      real(kind=dp) :: wrk(36,*), qd7(8,*),qd6(7,*),qd5(6,*),qd4(5,*)

      integer :: i, j, k, l
      real(kind=dp) :: a11, a13, a22, b23, b33, b3t, b44, d1s, d1t, d2s, &
        f41l03, f41l15, f41l1h, f42l03, f42l15, f42l1h, &
        f51k03, f51k06, f51k10, f51k15, f51k1h, f51k45, &
        f52k03, f52k06, f52k09, f52k15, f52k45, &
        f53k03, f53k15, f53k45, f54k03, f54k1h, &
        f61j06, f61j10, f61j15, f61j21, f62j03, f62j06, f62j10, f62j15, &
        f63j03, f63j06, f63j10, f64j03, f64j06, f64j10, &
        f65j03, f65j15, f66j21, &
        r11, r13, r22, s11, s13, s22, s23, s33, s3t, s44, u1s, u1t, u2s


      do i=i1,i2
         j=i+j0
         k=i+k0
         l=i+l0
         f41l03= qd4(1,l)* 3
         f41l15= qd4(1,l)*15
         f41l1h= f41l15  *7
         f42l03= qd4(2,l)* 3
         f42l15= qd4(2,l)*15
         f42l1h= f42l15  *7

         f51k03= qd5(1,k)* 3
         f51k06= f51k03+f51k03
         f51k10= qd5(1,k)*10
         f51k15= qd5(1,k)*15
         f51k45= f51k15  * 3
         f51k1h= f51k15  *7
         f52k03= qd5(2,k)* 3
         f52k06= f52k03+f52k03
         f52k09= f52k03+f52k06
         f52k15= qd5(2,k)*15
         f52k45= f52k15  * 3
         f53k03= qd5(3,k)* 3
         f53k15= qd5(3,k)*15
         f53k45= f53k15  * 3
         f54k03= qd5(4,k)* 3
         f54k1h= f54k03  * 5 *7

         f61j06= qd6(1,j)* 6
         f61j10= qd6(1,j)*10
         f61j15= qd6(1,j)*15
         f61j21= qd6(1,j)*21
         f62j03= qd6(2,j)* 3
         f62j06= f62j03+f62j03
         f62j10= qd6(2,j)*10
         f62j15= qd6(2,j)*15
         f63j03= qd6(3,j)* 3
         f63j06= f63j03+f63j03
         f63j10= qd6(3,j)*10
         f64j03= qd6(4,j)* 3
         f64j06= f64j03+f64j03
         f64j10= qd6(4,j)*10
         f65j03= qd6(5,j)* 3
         f65j15= qd6(5,j)*15
         f66j21= qd6(6,j)*21

         a11   =     qd5(1,k)*rdat%aqx2-qd4(1,l)
         a13   =     qd5(1,k)*rdat%aqx2-f41l03
         a22   =     qd5(2,k)*rdat%aqx2-qd4(2,l)
         r11   =     f51k03  *rdat%acy2-f41l03
         r13   =     qd5(1,k)*rdat%acy2-f41l03
         r22   =     f52k03  *rdat%acy2-f42l03

         b23   =     f62j03  *rdat%aqx2-f52k09
         b33   =     qd6(3,j)*rdat%aqx2-qd5(3,k)
         b3t   =     qd6(3,j)*rdat%aqx2-f53k03
         b44   =     qd6(4,j)*rdat%aqx2-qd5(4,k)
         s11   =     qd6(1,j)*rdat%acy2-qd5(1,k)
         s13   =     qd6(1,j)*rdat%acy2-f51k03
         s22   =     f62j03  *rdat%acy2-f52k03
         s23   =     f62j03  *rdat%acy2-f52k09
         s33   =     f63j03  *rdat%acy2-f53k03
         s3t   =     qd6(3,j)*rdat%acy2-f53k03
         s44   =     qd6(4,j)*rdat%acy2-qd5(4,k)

         d1s   =    (qd6(1,j)*rdat%aqx2-f51k06  )*rdat%aqx2+f41l03
         d1t   =    (qd6(1,j)*rdat%aqx2-f51k10  )*rdat%aqx2+f41l15
         d2s   =    (qd6(2,j)*rdat%aqx2-f52k06  )*rdat%aqx2+f42l03
         u1s   =    (qd6(1,j)*rdat%acy2-f51k06  )*rdat%acy2+f41l03
         u1t   =    (qd6(1,j)*rdat%acy2-f51k10  )*rdat%acy2+f41l15
         u2s   =    (qd6(2,j)*rdat%acy2-f52k06  )*rdat%acy2+f42l03

         wrk( 1,i)=((qd7(1,i)*rdat%aqx2-f61j21  )*rdat%aqx2+f51k1h)*rdat%aqx2-f41l1h
         wrk( 2,i)=((qd7(1,i)*rdat%aqx2-f61j15  )*rdat%aqx2+f51k45)*rdat%aqx2-f41l15
         wrk( 3,i)=((qd7(2,i)*rdat%aqx2-f62j15  )*rdat%aqx2+f52k45)*rdat%aqx2-f42l15
         wrk( 4,i)=((qd7(1,i)*rdat%aqx2-f61j10  )*rdat%aqx2+f51k15)*rdat%acy2-d1t
         wrk( 5,i)= (qd7(2,i)*rdat%aqx2-f62j10  )*rdat%aqx2+f52k15
         wrk( 6,i)= (qd7(3,i)*rdat%aqx2-f63j10  )*rdat%aqx2+f53k15      -d1t
         wrk( 7,i)=((qd7(1,i)*rdat%aqx2-f61j06  )*rdat%aqx2+f51k03)*rdat%acy2-d1s* 3
         wrk( 8,i)=((qd7(2,i)*rdat%aqx2-f62j06  )*rdat%aqx2+f52k03)*rdat%acy2-d2s
         wrk( 9,i)= (qd7(3,i)*rdat%aqx2-f63j06  )*rdat%aqx2+f53k03      -d1s
         wrk(10,i)= (qd7(4,i)*rdat%aqx2-f64j06  )*rdat%aqx2+f54k03      -d2s* 3
         wrk(11,i)=((qd7(1,i)*rdat%acy2-f61j06  )*rdat%acy2+f51k03)*rdat%aqx2-u1s* 3
         wrk(12,i)= (qd7(2,i)*rdat%acy2-f62j03  )*rdat%aqx2-s23
         wrk(13,i)= (qd7(3,i)*rdat%acy2-qd6(3,j)      -s11   )*rdat%aqx2-s33+r11
         wrk(14,i)=  qd7(4,i)*rdat%aqx2-f64j03        -b23
         wrk(15,i)=  qd7(5,i)*rdat%aqx2-f65j03        -b3t* 6      +a13* 3
         wrk(16,i)=((qd7(1,i)*rdat%acy2-f61j10  )*rdat%acy2+f51k15)*rdat%aqx2-u1t
         wrk(17,i)=((qd7(2,i)*rdat%acy2-f62j06  )*rdat%acy2+f52k03)*rdat%aqx2-u2s
         wrk(18,i)= (qd7(3,i)*rdat%acy2-f63j03        -s13   )*rdat%aqx2-s3t+r13
         wrk(19,i)= (qd7(4,i)*rdat%acy2-qd6(4,j)      -s22   )*rdat%aqx2-s44+r22
         wrk(20,i)=  qd7(5,i)*rdat%aqx2-qd6(5,j)      -b33* 6      +a11* 3
         wrk(21,i)=  qd7(6,i)*rdat%aqx2-qd6(6,j)      -b44*10     +a22*15
         wrk(22,i)=((qd7(1,i)*rdat%acy2-f61j15  )*rdat%acy2+f51k45)*rdat%acy2-f41l15
         wrk(23,i)= (qd7(2,i)*rdat%acy2-f62j10  )*rdat%acy2+f52k15
         wrk(24,i)= (qd7(3,i)*rdat%acy2-f63j06  )*rdat%acy2+f53k03      -u1s
         wrk(25,i)=  qd7(4,i)*rdat%acy2-f64j03        -s23
         wrk(26,i)=  qd7(5,i)*rdat%acy2-qd6(5,j)      -s33-s33     +r11
         wrk(27,i)=  qd7(6,i)     -f64j10        +f52k15
         wrk(28,i)=  qd7(7,i)     -f65j15        +f53k45      -f41l15
         wrk(29,i)=((qd7(1,i)*rdat%acy2-f61j21  )*rdat%acy2+f51k1h)*rdat%acy2-f41l1h
         wrk(30,i)=((qd7(2,i)*rdat%acy2-f62j15  )*rdat%acy2+f52k45)*rdat%acy2-f42l15
         wrk(31,i)= (qd7(3,i)*rdat%acy2-f63j10  )*rdat%acy2+f53k15      -u1t
         wrk(32,i)= (qd7(4,i)*rdat%acy2-f64j06  )*rdat%acy2+f54k03      -u2s* 3
         wrk(33,i)=  qd7(5,i)*rdat%acy2-f65j03        -s3t* 6      +r13* 3
         wrk(34,i)=  qd7(6,i)*rdat%acy2-qd6(6,j)      -s44*10     +r22* 5
         wrk(35,i)=  qd7(7,i)     -f65j15        +f53k45      -f41l15
         wrk(36,i)=  qd7(8,i)     -f66j21        +f54k1h      -f42l1h
      enddo

      end subroutine frikr7

! >
! >    @brief   auxiliary routine of order 8 for rot.axis integrations
! >
! >    @details auxiliary routine of order 8 for rot.axis integrations
! >
      subroutine frikr8(rdat,i1,i2,wrk,qd8,j0,qd7,k0,qd6,l0,qd5,m0,qd4)

      implicit none
      type(rotaxis_data_t) :: rdat

      integer :: i1, i2, j0, k0, l0, m0
      real(kind=dp) :: wrk(45,*), qd8(9,*),qd7(8,*),qd6(7,*),qd5(6,*),qd4(5,*)

      integer :: i, j, k, l, m
      real(kind=dp) :: cy, qx
      real(kind=dp) :: a11, b13, b22, b23, b33, c33, c44, c4t, c55, &
          d11, e16, e1t, e26, e2t, e36, f41m03, f41m15, f41m1h, &
          f51l03, f51l06, f51l09, f51l15, f51l1h, f51l45, f51l4h, &
          f52l03, f52l15, f52l1h, f53l03, f53l15, f53l4h, f61k03, &
          f61k06, f61k10, f61k15, f61k1h, f61k2h, f61k45, f62k03, &
          f62k06, f62k09, f62k10, f62k15, f62k1h, f62k45, f63k03, &
          f63k06, f63k15, f63k45, f64k03, f64k15, f64k1h, f65k03, &
          f65k2h, f71j06, f71j10, f71j15, f71j21, f71j28, f72j03, &
          f72j06, f72j10, f72j15, f72j21, f73j03, f73j06, f73j10, &
          f73j15, f74j03, f74j06, f74j10, f75j03, f75j06, f75j15, &
          f76j03, f76j21, f77j28, g11, r11, s11, s13, s22, s23, &
          s33, t13, t22, t23, t33, t3t, t43, t44, t4d, t55, u11, &
          v16, v1t, v26, v2t, v36, w11

      cy= rdat%acy2
      qx= rdat%aqx2
      do i=i1,i2
         j=i+j0
         k=i+k0
         l=i+l0
         m=i+m0
         f41m03= qd4(1,m)* 3
         f41m15= qd4(1,m)*15
         f41m1h= f41m15  *7

         f51l03= qd5(1,l)* 3
         f51l06= f51l03+f51l03
         f51l09= f51l03+f51l06
         f51l15= qd5(1,l)*15
         f51l45= f51l15  * 3
         f51l1h= f51l15  *7
         f51l4h= qd5(1,l)*420
         f52l03= qd5(2,l)* 3
         f52l15= qd5(2,l)*15
         f52l1h= f52l15  *7
         f53l03= qd5(3,l)* 3
         f53l15= qd5(3,l)*15
         f53l4h= qd5(3,l)*420

         f61k03= qd6(1,k)* 3
         f61k06= f61k03+f61k03
         f61k10= qd6(1,k)*10
         f61k15= qd6(1,k)*15
         f61k45= f61k15  * 3
         f61k1h= f61k15  *7
         f61k2h= qd6(1,k)*210
         f62k03= qd6(2,k)* 3
         f62k06= f62k03+f62k03
         f62k09= f62k03+f62k06
         f62k10= qd6(2,k)*10
         f62k15= qd6(2,k)*15
         f62k45= f62k15  * 3
         f62k1h= f62k15  *7
         f63k03= qd6(3,k)* 3
         f63k06= f63k03+f63k03
         f63k15= qd6(3,k)*15
         f63k45= f63k15  * 3
         f64k03= qd6(4,k)* 3
         f64k15= qd6(4,k)*15
         f64k1h= f64k15  *7
         f65k03= qd6(5,k)* 3
         f65k2h= qd6(5,k)*210

         f71j06= qd7(1,j)* 6
         f71j10= qd7(1,j)*10
         f71j15= qd7(1,j)*15
         f71j21= qd7(1,j)*21
         f71j28= qd7(1,j)*28
         f72j03= qd7(2,j)* 3
         f72j06= f72j03+f72j03
         f72j10= qd7(2,j)*10
         f72j15= qd7(2,j)*15
         f72j21= qd7(2,j)*21
         f73j03= qd7(3,j)* 3
         f73j06= f73j03+f73j03
         f73j10= qd7(3,j)*10
         f73j15= qd7(3,j)*15
         f74j03= qd7(4,j)* 3
         f74j06= f74j03+f74j03
         f74j10= qd7(4,j)*10
         f75j03= qd7(5,j)* 3
         f75j06= f75j03+f75j03
         f75j15= qd7(5,j)*15
         f76j03= qd7(6,j)* 3
         f76j21= qd7(6,j)*21
         f77j28= qd7(7,j)*28

         a11   =     qd5(1,l)*qx-qd4(1,m)
         r11   =     qd5(1,l)*cy-qd4(1,m)

         b13   =     f61k03  *qx-f51l09
         b22   =     f62k03  *qx-f52l03
         b23   =     f62k15  *qx-f52l15  * 3
         b33   =     f63k03  *qx-f53l03
         s11   =     f61k03  *cy-f51l03
         s13   =     f61k03  *cy-f51l09
         s22   =     f62k03  *cy-f52l03
         s23   =     f62k03  *cy-f52l03  * 3
         s33   =     f63k03  *cy-f53l03

         c33   =     f73j06  *qx-f63k03  * 6
         c44   =     qd7(4,j)*qx-qd6(4,k)
         c4t   =     f74j10  *qx-f64k03  *10
         c55   =     qd7(5,j)*qx-qd6(5,k)
         t13   =     qd7(1,j)*cy-f61k03
         t22   =     f72j03  *cy-f62k03
         t23   =     f72j03  *cy-f62k09
         t33   =     qd7(3,j)*cy-qd6(3,k)
         t3t   =     qd7(3,j)*cy-f63k03
         t44   =     qd7(4,j)*cy-qd6(4,k)
         t43   =     qd7(4,j)*cy-f64k03
         t4d   =     f74j10  *cy-f64k03  *10
         t55   =     qd7(5,j)*cy-qd6(5,k)

         d11   =    (qd6(1,k)*qx-f51l06  )*qx+f41m03
         u11   =    (qd6(1,k)*cy-f51l06  )*cy+f41m03

         e16   =    (qd7(1,j)*qx-f61k06  )*qx+f51l03
         e1t   =    (qd7(1,j)*qx-f61k10  )*qx+f51l15
         e26   =   ((qd7(2,j)*qx-f62k06  )*qx+f52l03)* 3
         e2t   =    (qd7(2,j)*qx-f62k10  )*qx+f52l15
         e36   =    (qd7(3,j)*qx-f63k06  )*qx+f53l03
         v16   =    (qd7(1,j)*cy-f61k06  )*cy+f51l03
         v1t   =    (qd7(1,j)*cy-f61k10  )*cy+f51l15
         v26   =   ((qd7(2,j)*cy-f62k06  )*cy+f52l03)* 3
         v2t   =    (qd7(2,j)*cy-f62k10  )*cy+f52l15
         v36   =    (qd7(3,j)*cy-f63k06  )*cy+f53l03

         g11   =   ((qd7(1,j)*qx-f61k15  )*qx+f51l45)*qx-f41m15
         w11   =   ((qd7(1,j)*cy-f61k15  )*cy+f51l45)*cy-f41m15

         wrk( 1,i)=((qd8(1,i)*qx-f71j28  )*qx+f61k2h)*qx-f51l4h
         wrk( 1,i)=              wrk( 1,i)*qx+f41m1h
         wrk( 2,i)=((qd8(1,i)*qx-f71j21  )*qx+f61k1h)*qx-f51l1h
         wrk( 3,i)=((qd8(2,i)*qx-f72j21  )*qx+f62k1h)*qx-f52l1h
         wrk( 4,i)=((qd8(1,i)*qx-f71j15  )*qx+f61k45)*qx-f51l15
         wrk( 4,i)=              wrk( 4,i)*cy-g11
         wrk( 5,i)=((qd8(2,i)*qx-f72j15  )*qx+f62k45)*qx-f52l15
         wrk( 6,i)=((qd8(3,i)*qx-f73j15  )*qx+f63k45)*qx-f53l15 -g11
         wrk( 7,i)=((qd8(1,i)*qx-f71j10  )*qx+f61k15)*cy -e1t* 3
         wrk( 8,i)=((qd8(2,i)*qx-f72j10  )*qx+f62k15)*cy -e2t
         wrk( 9,i)= (qd8(3,i)*qx-f73j10  )*qx+f63k15     -e1t
         wrk(10,i)= (qd8(4,i)*qx-f74j10  )*qx+f64k15     -e2t* 3
         wrk(11,i)= (qd8(1,i)*cy-f71j06  )*cy+f61k03
         wrk(11,i)=             (wrk(11,i)*qx-v16* 6 )*qx+u11* 3
         wrk(12,i)=((qd8(2,i)*qx-f72j06  )*qx+f62k03)*cy -e26
         wrk(13,i)=((qd8(3,i)*qx-f73j06  )*qx+f63k03-e16)*cy-e36+d11
         wrk(14,i)= (qd8(4,i)*qx-f74j06  )*qx+f64k03-e26
         wrk(15,i)= (qd8(5,i)*qx-f75j06  )*qx+f65k03-e36* 6 +d11* 3
         wrk(16,i)=((qd8(1,i)*cy-f71j10  )*cy+f61k15)*qx-v1t* 3
         wrk(17,i)=((qd8(2,i)*cy-f72j06  )*cy+f62k03)*qx-v26
         wrk(18,i)= (qd8(3,i)*cy-f73j03  -t13)*qx-t3t* 3 +s13
         wrk(19,i)= (qd8(4,i)*cy-qd7(4,j)-t22)*qx-t44* 3 +s22* 3
         wrk(20,i)=  qd8(5,i)*qx-f75j03          -c33+b13
         wrk(21,i)=  qd8(6,i)*qx-f76j03          -c4t+b23
         wrk(22,i)=((qd8(1,i)*cy-f71j15  )*cy+f61k45)*cy-f51l15
         wrk(22,i)=              wrk(22,i)*qx-w11
         wrk(23,i)=((qd8(2,i)*cy-f72j10  )*cy+f62k15)*qx -v2t
         wrk(24,i)=((qd8(3,i)*cy-f73j06  )*cy+f63k03-v16)*qx-v36+u11
         wrk(25,i)= (qd8(4,i)*cy-f74j03  -t23)*qx-t43+s23
         wrk(26,i)=  qd8(5,i)*cy-qd7(5,j)-t33* 6 +s11
         wrk(26,i)=              wrk(26,i)*qx-(t55-s33-s33+r11* 3 )
         wrk(27,i)=  qd8(6,i)*qx-qd7(6,j)    -(c44*10-b22* 5 )
         wrk(28,i)=  qd8(7,i)*qx-qd7(7,j)    -(c55-b33+a11)*15
         wrk(29,i)=((qd8(1,i)*cy-f71j21  )*cy+f61k1h)*cy-f51l1h
         wrk(30,i)=((qd8(2,i)*cy-f72j15  )*cy+f62k45)*cy-f52l15
         wrk(31,i)= (qd8(3,i)*cy-f73j10  )*cy+f63k15-v1t
         wrk(32,i)= (qd8(4,i)*cy-f74j06  )*cy+f64k03-v26
         wrk(33,i)=  qd8(5,i)*cy-f75j03      -t3t* 6 +s13
         wrk(34,i)=  qd8(6,i)*cy-qd7(6,j)    -t44*10+s22* 5
         wrk(35,i)=  qd8(7,i)   -f75j15  +f63k45-f51l15
         wrk(36,i)=  qd8(8,i)   -f76j21  +f64k1h-f52l1h
         wrk(37,i)=((qd8(1,i)*cy-f71j28  )*cy+f61k2h)*cy-f51l4h
         wrk(37,i)=              wrk(37,i)*cy+f41m1h
         wrk(38,i)=((qd8(2,i)*cy-f72j21  )*cy+f62k1h)*cy-f52l1h
         wrk(39,i)=((qd8(3,i)*cy-f73j15  )*cy+f63k45)*cy-f53l15-w11
         wrk(40,i)= (qd8(4,i)*cy-f74j10  )*cy+f64k15   -v2t* 3
         wrk(41,i)= (qd8(5,i)*cy-f75j06  )*cy+f65k03   -v36* 6 +u11* 3
         wrk(42,i)=  qd8(6,i)*cy-f76j03                -t4d+s23* 5
         wrk(43,i)=  qd8(7,i)*cy-qd7(7,j)         -(t55-s33+r11)*15
         wrk(44,i)=  qd8(8,i)   -f76j21      +f64k1h   -f52l1h
         wrk(45,i)=  qd8(9,i)   -f77j28      +f65k2h   -f53l4h+f41m1h
      enddo

      end subroutine frikr8

      subroutine r30s1d_02(f,p)

      implicit none

      real(kind=dp) :: f(3,1,1,*), p(3,3)
      real(kind=dp) :: t(6)

      t (1:3) = f (1:3,1,1,1)
      f (1,1,1,1) = t (1)*p (1,1) + t (2)*p (2,1) + t (3)*p (3,1)
      f (2,1,1,1) = t (1)*p (1,2) + t (2)*p (2,2) + t (3)*p (3,2)
      f (3,1,1,1) = t (1)*p (1,3) + t (2)*p (2,3) + t (3)*p (3,3)

      end subroutine r30s1d_02

      subroutine r30s1d_03(f,p)

      implicit none

      real(kind=dp) :: f(3,3,1,*), p(3,3)
      real(kind=dp) :: t(6)
      integer :: i, j, k, l

     do i = 1, 3
        t (1) = f (i,1,1,1)
        t (2) = f (i,2,1,1)
        t (3) = f (i,3,1,1)
        f (i,1,1,1) = t(1) * p(1,1) + t(2) * p(2,1) + t(3) * p(3,1)
        f (i,2,1,1) = t(1) * p(1,2) + t(2) * p(2,2) + t(3) * p(3,2)
        f (i,3,1,1) = t(1) * p(1,3) + t(2) * p(2,3) + t(3) * p(3,3)
     end do
     do j = 1, 3
        t(1) = f(1,j,1,1)
        t(2) = f(2,j,1,1)
        t(3) = f(3,j,1,1)
        f(1,j,1,1) = t(1) * p(1,1) + t(2) * p(2,1) + t(3) * p(3,1)
        f(2,j,1,1) = t(1) * p(1,2) + t(2) * p(2,2) + t(3) * p(3,2)
        f(3,j,1,1) = t(1) * p(1,3) + t(2) * p(2,3) + t(3) * p(3,3)
     end do

      end subroutine r30s1d_03

      subroutine r30s1d_04(f,p)

      implicit none

      real(kind=dp) :: f(3,1,3,*), p(3,3)
      real(kind=dp) :: t(6)
      integer :: i, j, k, l

         do i = 1, 3
            t(1) = f(i,1,1,1)
            t(2) = f(i,1,2,1)
            t(3) = f(i,1,3,1)
            f(i,1,1,1) = t(1) * p(1,1) + t(2) * p(2,1) + t(3) * p(3,1)
            f(i,1,2,1) = t(1) * p(1,2) + t(2) * p(2,2) + t(3) * p(3,2)
            f(i,1,3,1) = t(1) * p(1,3) + t(2) * p(2,3) + t(3) * p(3,3)
         end do
         do k = 1, 3
            t(1) = f(1,1,k,1)
            t(2) = f(2,1,k,1)
            t(3) = f(3,1,k,1)
            f(1,1,k,1) = t(1) * p(1,1) + t(2) * p(2,1) + t(3) * p(3,1)
            f(2,1,k,1) = t(1) * p(1,2) + t(2) * p(2,2) + t(3) * p(3,2)
            f(3,1,k,1) = t(1) * p(1,3) + t(2) * p(2,3) + t(3) * p(3,3)
         end do

      end subroutine r30s1d_04

      subroutine r30s1d_05(f,p)

      implicit none

      real(kind=dp) :: f(3,3,3,*), p(3,3)
      real(kind=dp) :: t(6)
      integer :: i, j, k, l

         do j = 1, 3
            do i = 1, 3
               t(1) = f(i,j,1,1)
               t(2) = f(i,j,2,1)
               t(3) = f(i,j,3,1)
               f(i,j,1,1) = t(1) * p(1,1) + t(2) * p(2,1) + t(3) * p(3,1)
               f(i,j,2,1) = t(1) * p(1,2) + t(2) * p(2,2) + t(3) * p(3,2)
               f(i,j,3,1) = t(1) * p(1,3) + t(2) * p(2,3) + t(3) * p(3,3)
            end do
         end do
         do k = 1, 3
            do i = 1, 3
               t(1) = f(i,1,k,1)
               t(2) = f(i,2,k,1)
               t(3) = f(i,3,k,1)
               f(i,1,k,1) = t(1) * p(1,1) + t(2) * p(2,1) + t(3) * p(3,1)
               f(i,2,k,1) = t(1) * p(1,2) + t(2) * p(2,2) + t(3) * p(3,2)
               f(i,3,k,1) = t(1) * p(1,3) + t(2) * p(2,3) + t(3) * p(3,3)
            end do
         end do
         do k = 1, 3
            do j = 1, 3
               t(1) = f(1,j,k,1)
               t(2) = f(2,j,k,1)
               t(3) = f(3,j,k,1)
               f(1,j,k,1) = t(1) * p(1,1) + t(2) * p(2,1) + t(3) * p(3,1)
               f(2,j,k,1) = t(1) * p(1,2) + t(2) * p(2,2) + t(3) * p(3,2)
               f(3,j,k,1) = t(1) * p(1,3) + t(2) * p(2,3) + t(3) * p(3,3)
            end do
         end do

      end subroutine r30s1d_05

      subroutine r30s1d_06(f,p)

      implicit none

      real(kind=dp) :: f(3,3,3,*), p(3,3)
      real(kind=dp) :: t(6)
      integer :: i, j, k, l

         do k = 1, 3
            do j = 1, 3
               do i = 1, 3
                  t(1) = f(i,j,k,1)
                  t(2) = f(i,j,k,2)
                  t(3) = f(i,j,k,3)
                  f(i,j,k,1) = t(1) * p(1,1) + t(2) * p(2,1) + t(3) * p(3,1)
                  f(i,j,k,2) = t(1) * p(1,2) + t(2) * p(2,2) + t(3) * p(3,2)
                  f(i,j,k,3) = t(1) * p(1,3) + t(2) * p(2,3) + t(3) * p(3,3)
               end do
            end do
         end do
         do l = 1, 3
            do j = 1, 3
               do i = 1, 3
                  t(1) = f(i,j,1,l)
                  t(2) = f(i,j,2,l)
                  t(3) = f(i,j,3,l)
                  f(i,j,1,l) = t(1) * p(1,1) + t(2) * p(2,1) + t(3) * p(3,1)
                  f(i,j,2,l) = t(1) * p(1,2) + t(2) * p(2,2) + t(3) * p(3,2)
                  f(i,j,3,l) = t(1) * p(1,3) + t(2) * p(2,3) + t(3) * p(3,3)
               end do
            end do
         end do
         do l = 1, 3
            do k = 1, 3
               do i = 1, 3
                  t(1) = f(i,1,k,l)
                  t(2) = f(i,2,k,l)
                  t(3) = f(i,3,k,l)
                  f(i,1,k,l) = t(1) * p(1,1) + t(2) * p(2,1) + t(3) * p(3,1)
                  f(i,2,k,l) = t(1) * p(1,2) + t(2) * p(2,2) + t(3) * p(3,2)
                  f(i,3,k,l) = t(1) * p(1,3) + t(2) * p(2,3) + t(3) * p(3,3)
               end do
            end do
         end do
         do l = 1, 3
            do k = 1, 3
               do j = 1, 3
                  t(1) = f(1,j,k,l)
                  t(2) = f(2,j,k,l)
                  t(3) = f(3,j,k,l)
                  f(1,j,k,l) = t(1) * p(1,1) + t(2) * p(2,1) + t(3) * p(3,1)
                  f(2,j,k,l) = t(1) * p(1,2) + t(2) * p(2,2) + t(3) * p(3,2)
                  f(3,j,k,l) = t(1) * p(1,3) + t(2) * p(2,3) + t(3) * p(3,3)
               end do
            end do
         end do

      end subroutine r30s1d_06


      subroutine r30s1d_07(f,p)

      implicit none

      real(kind=dp) :: f(6,1,1,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      t(1:6)= f(1:6,1,1,1)
      do i=1,6
         f(i,1,1,1) = t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)             &
                     +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
      end do

      end subroutine r30s1d_07

      subroutine r30s1d_08(f,p)

      implicit none

      real(kind=dp) :: f(6,3,1,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do i=1,6
         t(1)= f(i,1,1,1)
         t(2)= f(i,2,1,1)
         t(3)= f(i,3,1,1)
         f(i,1,1,1)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
         f(i,2,1,1)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
         f(i,3,1,1)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
      end do
      do j=1,3
         t(1:6)= f(1:6,j,1,1)
         do i=1,6
            f(i,j,1,1)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)          &
                       +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
         end do
      end do

      end subroutine r30s1d_08

      subroutine r30s1d_09(f,p)

      implicit none

      real(kind=dp) :: f(6,1,3,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do i=1,6
         t(1)= f(i,1,1,1)
         t(2)= f(i,1,2,1)
         t(3)= f(i,1,3,1)
         f(i,1,1,1)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
         f(i,1,2,1)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
         f(i,1,3,1)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
      end do
      do j=1,3
        t(1:6)= f(1:6,1,j,1)
         do i=1,6
            f(i,1,j,1)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)          &
                       +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
         end do
      end do

      end subroutine r30s1d_09

      subroutine r30s1d_10(f,p)

      implicit none

      real(kind=dp) :: f(6,6,1,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do i=1,6
         t(1:6)= f(i,1:6,1,1)
         do j=1,6
            f(i,j,1,1)= t(1)*q(1,j)+t(2)*q(2,j)+t(3)*q(3,j)          &
                       +t(4)*q(4,j)+t(5)*q(5,j)+t(6)*q(6,j)
         end do
      end do
      do j=1,6
         t(1:6)= f(1:6,j,1,1)
         do i=1,6
            f(i,j,1,1)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)          &
                       +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
         end do
      end do

      end subroutine r30s1d_10

      subroutine r30s1d_11(f,p)

      implicit none

      real(kind=dp) :: f(6,3,3,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do j=1,3
         do i=1,6
            t(1)= f(i,j,1,1)
            t(2)= f(i,j,2,1)
            t(3)= f(i,j,3,1)
            f(i,j,1,1)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
            f(i,j,2,1)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
            f(i,j,3,1)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
         end do
      end do
      do k=1,3
         do i=1,6
            t(1)= f(i,1,k,1)
            t(2)= f(i,2,k,1)
            t(3)= f(i,3,k,1)
            f(i,1,k,1)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
            f(i,2,k,1)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
            f(i,3,k,1)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
         end do
      end do
      do k=1,3
         do j=1,3
            do i=1,6
               t(i)= f(i,j,k,1)
            end do
            do i=1,6
               f(i,j,k,1)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)       &
                          +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
            end do
         end do
      end do

      end subroutine r30s1d_11

      subroutine r30s1d_12(f,p)

      implicit none

      real(kind=dp) :: f(6,1,6,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do i=1,6
         do k=1,6
            t(k)= f(i,1,k,1)
         end do
         do k=1,6
            f(i,1,k,1)= t(1)*q(1,k)+t(2)*q(2,k)+t(3)*q(3,k)          &
                       +t(4)*q(4,k)+t(5)*q(5,k)+t(6)*q(6,k)
         end do
      end do
      do k=1,6
         do i=1,6
            t(i)= f(i,1,k,1)
         end do
         do i=1,6
            f(i,1,k,1)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)          &
                       +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
         end do
      end do

      end subroutine r30s1d_12

      subroutine r30s1d_13(f,p)

      implicit none

      real(kind=dp) :: f(6,1,3,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do k=1,3
         do i=1,6
            t(1)= f(i,1,k,1)
            t(2)= f(i,1,k,2)
            t(3)= f(i,1,k,3)
            f(i,1,k,1)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
            f(i,1,k,2)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
            f(i,1,k,3)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
         end do
      end do
      do l=1,3
         do i=1,6
            t(1)= f(i,1,1,l)
            t(2)= f(i,1,2,l)
            t(3)= f(i,1,3,l)
            f(i,1,1,l)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
            f(i,1,2,l)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
            f(i,1,3,l)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
         end do
      end do
      do l=1,3
         do k=1,3
            do i=1,6
               t(i)= f(i,1,k,l)
            end do
            do i=1,6
               f(i,1,k,l)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)       &
                          +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
            end do
         end do
      end do
      end subroutine r30s1d_13

      subroutine r30s1d_14(f,p)

      implicit none

      real(kind=dp) :: f(6,6,3,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do j=1,6
         do i=1,6
            t(1)= f(i,j,1,1)
            t(2)= f(i,j,2,1)
            t(3)= f(i,j,3,1)
            f(i,j,1,1)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
            f(i,j,2,1)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
            f(i,j,3,1)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
         end do
      end do
      do k=1,3
         do i=1,6
            do j=1,6
               t(j)= f(i,j,k,1)
            end do
            do j=1,6
               f(i,j,k,1)= t(1)*q(1,j)+t(2)*q(2,j)+t(3)*q(3,j)       &
                          +t(4)*q(4,j)+t(5)*q(5,j)+t(6)*q(6,j)
            end do
         end do
      end do
      do k=1,3
         do j=1,6
            do i=1,6
               t(i)= f(i,j,k,1)
            end do
            do i=1,6
               f(i,j,k,1)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)       &
                          +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
            end do
         end do
      end do
      end subroutine r30s1d_14

      subroutine r30s1d_15(f,p)

      implicit none

      real(kind=dp) :: f(6,3,6,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do j=1,3
         do i=1,6
            do k=1,6
               t(k)= f(i,j,k,1)
            end do
            do k=1,6
               f(i,j,k,1)= t(1)*q(1,k)+t(2)*q(2,k)+t(3)*q(3,k)       &
                          +t(4)*q(4,k)+t(5)*q(5,k)+t(6)*q(6,k)
            end do
         end do
      end do
      do k=1,6
         do i=1,6
            t(1)= f(i,1,k,1)
            t(2)= f(i,2,k,1)
            t(3)= f(i,3,k,1)
            f(i,1,k,1)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
            f(i,2,k,1)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
            f(i,3,k,1)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
         end do
      end do
      do k=1,6
         do j=1,3
            do i=1,6
               t(i)= f(i,j,k,1)
            end do
            do i=1,6
               f(i,j,k,1)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)       &
                          +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
            end do
         end do
      end do
      end subroutine r30s1d_15

      subroutine r30s1d_16(f,p)

      implicit none

      real(kind=dp) :: f(6,3,3,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do k=1,3
         do j=1,3
            do i=1,6
               t(1)= f(i,j,k,1)
               t(2)= f(i,j,k,2)
               t(3)= f(i,j,k,3)
               f(i,j,k,1)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
               f(i,j,k,2)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
               f(i,j,k,3)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
            end do
         end do
      end do
      do l=1,3
         do j=1,3
            do i=1,6
               t(1)= f(i,j,1,l)
               t(2)= f(i,j,2,l)
               t(3)= f(i,j,3,l)
               f(i,j,1,l)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
               f(i,j,2,l)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
               f(i,j,3,l)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
            end do
         end do
      end do
      do l=1,3
         do k=1,3
            do i=1,6
               t(1)= f(i,1,k,l)
               t(2)= f(i,2,k,l)
               t(3)= f(i,3,k,l)
               f(i,1,k,l)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
               f(i,2,k,l)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
               f(i,3,k,l)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
            end do
         end do
      end do
      do l=1,3
         do k=1,3
            do j=1,3
               do i=1,6
                  t(i)= f(i,j,k,l)
               end do
               do i=1,6
                  f(i,j,k,l)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)    &
                             +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
               end do
            end do
         end do
      end do
      end subroutine r30s1d_16

      subroutine r30s1d_17(f,p)

      implicit none

      real(kind=dp) :: f(6,6,6,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do j=1,6
         do i=1,6
            do k=1,6
               t(k)= f(i,j,k,1)
            end do
            do k=1,6
               f(i,j,k,1)= t(1)*q(1,k)+t(2)*q(2,k)+t(3)*q(3,k)       &
                          +t(4)*q(4,k)+t(5)*q(5,k)+t(6)*q(6,k)
            end do
         end do
      end do
      do k=1,6
         do i=1,6
            do j=1,6
               t(j)= f(i,j,k,1)
            end do
            do j=1,6
               f(i,j,k,1)= t(1)*q(1,j)+t(2)*q(2,j)+t(3)*q(3,j)       &
                          +t(4)*q(4,j)+t(5)*q(5,j)+t(6)*q(6,j)
            end do
         end do
      end do
      do k=1,6
         do j=1,6
            do i=1,6
               t(i)= f(i,j,k,1)
            end do
            do i=1,6
               f(i,j,k,1)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)       &
                          +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
            end do
         end do
      end do
      end subroutine r30s1d_17

      subroutine r30s1d_18(f,p)

      implicit none

      real(kind=dp) :: f(6,6,3,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do k=1,3
         do j=1,6
            do i=1,6
               t(1)= f(i,j,k,1)
               t(2)= f(i,j,k,2)
               t(3)= f(i,j,k,3)
               f(i,j,k,1)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
               f(i,j,k,2)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
               f(i,j,k,3)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
            end do
         end do
      end do
      do l=1,3
         do j=1,6
            do i=1,6
               t(1)= f(i,j,1,l)
               t(2)= f(i,j,2,l)
               t(3)= f(i,j,3,l)
               f(i,j,1,l)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
               f(i,j,2,l)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
               f(i,j,3,l)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
            end do
         end do
      end do
      do l=1,3
         do k=1,3
            do i=1,6
               do j=1,6
                  t(j)= f(i,j,k,l)
               end do
               do j=1,6
                  f(i,j,k,l)= t(1)*q(1,j)+t(2)*q(2,j)+t(3)*q(3,j)    &
                             +t(4)*q(4,j)+t(5)*q(5,j)+t(6)*q(6,j)
               end do
            end do
         end do
      end do
      do l=1,3
         do k=1,3
            do j=1,6
               do i=1,6
                  t(i)= f(i,j,k,l)
               end do
               do i=1,6
                  f(i,j,k,l)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)    &
                             +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
               end do
            end do
         end do
      end do
      end subroutine r30s1d_18

      subroutine r30s1d_19(f,p)

      implicit none

      real(kind=dp) :: f(6,3,6,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do k=1,6
         do j=1,3
            do i=1,6
               t(1)= f(i,j,k,1)
               t(2)= f(i,j,k,2)
               t(3)= f(i,j,k,3)
               f(i,j,k,1)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
               f(i,j,k,2)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
               f(i,j,k,3)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
            end do
         end do
      end do
      do l=1,3
         do j=1,3
            do i=1,6
               do k=1,6
                  t(k)= f(i,j,k,l)
               end do
               do k=1,6
                  f(i,j,k,l)= t(1)*q(1,k)+t(2)*q(2,k)+t(3)*q(3,k)    &
                             +t(4)*q(4,k)+t(5)*q(5,k)+t(6)*q(6,k)
               end do
            end do
         end do
      end do
      do l=1,3
         do k=1,6
            do i=1,6
               t(1)= f(i,1,k,l)
               t(2)= f(i,2,k,l)
               t(3)= f(i,3,k,l)
               f(i,1,k,l)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
               f(i,2,k,l)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
               f(i,3,k,l)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
            end do
         end do
      end do
      do l=1,3
         do k=1,6
            do j=1,3
               do i=1,6
                  t(i)= f(i,j,k,l)
               end do
               do i=1,6
                  f(i,j,k,l)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)    &
                             +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
               end do
            end do
         end do
      end do
      end subroutine r30s1d_19

      subroutine r30s1d_20(f,p)

      implicit none

      real(kind=dp) :: f(6,6,6,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do k=1,6
         do j=1,6
            do i=1,6
               t(1)= f(i,j,k,1)
               t(2)= f(i,j,k,2)
               t(3)= f(i,j,k,3)
               f(i,j,k,1)= t(1)*p(1,1)+t(2)*p(2,1)+t(3)*p(3,1)
               f(i,j,k,2)= t(1)*p(1,2)+t(2)*p(2,2)+t(3)*p(3,2)
               f(i,j,k,3)= t(1)*p(1,3)+t(2)*p(2,3)+t(3)*p(3,3)
            end do
         end do
      end do
      do l=1,3
         do j=1,6
            do i=1,6
               do k=1,6
                  t(k)= f(i,j,k,l)
               end do
               do k=1,6
                  f(i,j,k,l)= t(1)*q(1,k)+t(2)*q(2,k)+t(3)*q(3,k)    &
                             +t(4)*q(4,k)+t(5)*q(5,k)+t(6)*q(6,k)
               end do
            end do
         end do
      end do
      do l=1,3
         do k=1,6
            do i=1,6
               do j=1,6
                  t(j)= f(i,j,k,l)
               end do
               do j=1,6
                  f(i,j,k,l)= t(1)*q(1,j)+t(2)*q(2,j)+t(3)*q(3,j)    &
                             +t(4)*q(4,j)+t(5)*q(5,j)+t(6)*q(6,j)
               end do
            end do
         end do
      end do
      do l=1,3
         do k=1,6
            do j=1,6
               do i=1,6
                  t(i)= f(i,j,k,l)
               end do
               do i=1,6
                  f(i,j,k,l)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)    &
                             +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
               end do
            end do
         end do
      end do
      end subroutine r30s1d_20

      subroutine r30s1d_21(f,p)

      implicit none

      real(kind=dp) :: f(6,6,6,*), p(3,3)
      real(kind=dp) :: t(6), q(6,6)
      integer :: i, j, k, l

      q(1,1:3) = p(1,1:3)*p(1,1:3)
      q(2,1:3) = p(2,1:3)*p(2,1:3)
      q(3,1:3) = p(3,1:3)*p(3,1:3)
      q(4,1:3) = p(1,1:3)*p(2,1:3)*2.0_dp
      q(5,1:3) = p(1,1:3)*p(3,1:3)*2.0_dp
      q(6,1:3) = p(2,1:3)*p(3,1:3)*2.0_dp

      q(1,4:5) = sqrt3 * (p(1,1)*p(1,2:3) )
      q(2,4:5) = sqrt3 * (p(2,1)*p(2,2:3) )
      q(3,4:5) = sqrt3 * (p(3,1)*p(3,2:3) )
      q(4,4:5) = sqrt3 * (p(1,1)*p(2,2:3)+p(2,1)*p(1,2:3) )
      q(5,4:5) = sqrt3 * (p(1,1)*p(3,2:3)+p(3,1)*p(1,2:3) )
      q(6,4:5) = sqrt3 * (p(2,1)*p(3,2:3)+p(3,1)*p(2,2:3) )

      q(1,6) = sqrt3 * ( p(1,2)*p(1,3))
      q(2,6) = sqrt3 * ( p(2,2)*p(2,3))
      q(3,6) = sqrt3 * ( p(3,2)*p(3,3))
      q(4,6) = sqrt3 * ( p(1,2)*p(2,3)+p(2,2)*p(1,3) )
      q(5,6) = sqrt3 * ( p(1,2)*p(3,3)+p(3,2)*p(1,3) )
      q(6,6) = sqrt3 * ( p(2,2)*p(3,3)+p(3,2)*p(2,3) )

      do k=1,6
         do j=1,6
            do i=1,6
               do l=1,6
                  t(l)= f(i,j,k,l)
               end do
               do l=1,6
                  f(i,j,k,l)= t(1)*q(1,l)+t(2)*q(2,l)+t(3)*q(3,l)    &
                             +t(4)*q(4,l)+t(5)*q(5,l)+t(6)*q(6,l)
               end do
            end do
         end do
      end do
      do l=1,6
         do j=1,6
            do i=1,6
               do k=1,6
                  t(k)= f(i,j,k,l)
               end do
               do k=1,6
                  f(i,j,k,l)= t(1)*q(1,k)+t(2)*q(2,k)+t(3)*q(3,k)    &
                             +t(4)*q(4,k)+t(5)*q(5,k)+t(6)*q(6,k)
               end do
            end do
         end do
      end do
      do l=1,6
         do k=1,6
            do i=1,6
               do j=1,6
                  t(j)= f(i,j,k,l)
               end do
               do j=1,6
                  f(i,j,k,l)= t(1)*q(1,j)+t(2)*q(2,j)+t(3)*q(3,j)    &
                             +t(4)*q(4,j)+t(5)*q(5,j)+t(6)*q(6,j)
               end do
            end do
         end do
      end do
      do l=1,6
         do k=1,6
            do j=1,6
               do i=1,6
                  t(i)= f(i,j,k,l)
               end do
               do i=1,6
                  f(i,j,k,l)= t(1)*q(1,i)+t(2)*q(2,i)+t(3)*q(3,i)    &
                             +t(4)*q(4,i)+t(5)*q(5,i)+t(6)*q(6,i)
               end do
            end do
         end do
      end do
      end subroutine r30s1d_21

end module
