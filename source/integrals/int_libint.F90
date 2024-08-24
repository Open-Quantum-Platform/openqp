module int2e_libint
  use precision, only: dp
  use int2_pairs, only: int2_pair_storage, int2_cutoffs_t
  use iso_c_binding, only: c_double, c_ptr, c_null_ptr, c_int, &
                           c_funptr, c_f_pointer, c_f_procpointer, c_size_t
  use libint_f
  use constants, only: pi

  implicit none

#ifdef OQP_LIBINT_ENABLE
#include "libint2/config.h"
#include "libint2/util/generated/libint2_params.h"
#include "fortran_incldefs.h"
#endif

  private
  public libint_compute_eri
  public libint_print_eri
  public libint_to_ghondo
  public libint_static_init
  public libint_static_cleanup

  public libint_t
  public libint2_active
  public libint2_init_eri
  public libint2_cleanup_eri
  public libint2_build

  real(dp), parameter :: halfsqrtpi = 0.5d0*sqrt(pi)

contains

  subroutine libint_static_init
    if (libint2_active) call libint2_static_init
  end subroutine

  subroutine libint_static_cleanup
    if (libint2_active) call libint2_static_cleanup
  end subroutine



  subroutine libint_compute_eri(basis, ppairs, cutoffs, &
                           shell_ids, deriv_order, &
                           erieval, flips, zero_shq)
      use boys_lut, only: tmax
      use basis_tools, only: basis_set
      implicit none
      type(basis_set), intent(in) :: basis
      type(int2_cutoffs_t), intent(in) :: cutoffs
      type(int2_pair_storage), intent(in) :: ppairs
      integer, intent(in) :: shell_ids(4)
      integer, intent(out) :: flips(4)
      integer, intent(in) :: deriv_order
      logical, intent(out) :: zero_shq
      real(kind=dp), dimension(3) :: a, b, c, d
      real(kind=dp), dimension(21) :: f
      type(libint_t), intent(out) :: erieval(*)
      real(kind=dp) :: pq2, rhoq, gammapq, pfac
      real(kind=dp), dimension(3) :: pq
      integer :: am(4), am_tot, p1234
      procedure(libint2_build), pointer :: build_eri
      integer :: shl_new(4)
      real(kind=dp) :: gpqinv
      real(kind=dp) :: tt, t2m1
      integer :: m

      integer :: p12, p34
      integer :: ppid_p, ppid_q, npp_p, npp_q
      integer :: id1, id2

      zero_shq = .true.

#ifdef OQP_LIBINT_ENABLE
      am = basis%ktype(shell_ids) - 1
      am_tot = sum(am) + deriv_order

      flips = [1,2,3,4]
      if (am(1)<am(2)) then
          flips(1:2) = [2,1]
          am(1:2) = am([2,1])
      end if
      if (am(3)<am(4)) then
          flips(3:4) = [4,3]
          am(3:4) = am([4,3])
      end if
      if (am(1)+am(2)>am(3)+am(4)) then
          flips = flips([3,4,1,2])
          am = am([3,4,1,2])
      end if
      shl_new = shell_ids(flips)

      id1 = maxval(shl_new([1,2]))
      id2 = minval(shl_new([1,2]))
      npp_p = ppairs%ppid(1,id1*(id1-1)/2+id2)
      ppid_p = ppairs%ppid(2,id1*(id1-1)/2+id2)

      id1 = maxval(shl_new([3,4]))
      id2 = minval(shl_new([3,4]))
      npp_q = ppairs%ppid(1,id1*(id1-1)/2+id2)
      ppid_q = ppairs%ppid(2,id1*(id1-1)/2+id2)

      if (npp_p == 0 .or. npp_q == 0) return

      a = basis%shell_centers(shl_new(1),1:3)
      b = basis%shell_centers(shl_new(2),1:3)
      c = basis%shell_centers(shl_new(3),1:3)
      d = basis%shell_centers(shl_new(4),1:3)

      p1234 = 0

      associate( &
          alpha1 => ppairs%alpha_a(ppid_p:), &
          alpha2 => ppairs%alpha_b(ppid_p:), &
          gammap => ppairs%g(ppid_p:),       &
          gpinv  => ppairs%ginv(ppid_p:),    &
          p      => ppairs%p(:,ppid_p:),     &
          k1     => ppairs%k(ppid_p:),       &

          alpha3 => ppairs%alpha_a(ppid_q:), &
          alpha4 => ppairs%alpha_b(ppid_q:), &
          gammaq => ppairs%g(ppid_q:),       &
          gqinv  => ppairs%ginv(ppid_q:),    &
          q      => ppairs%p(:,ppid_q:),     &
          k2     => ppairs%k(ppid_q:) )

      DO p12 = 1, npp_p
          DO p34 = 1, npp_q

            if ((k1(p12)*k2(p34))**2 < &
                    (gammap(p12)*gammaq(p34))**2 * (gammap(p12)+gammaq(p34)) * &
                    cutoffs%pair_cutoff_squared) cycle

            gpqinv = 1/(gammap(p12)+gammaq(p34))
            pfac = k1(p12)*k2(p34)*gpinv(p12)*gqinv(p34)*sqrt(gpqinv)
!            if (abs(pfac) < cutoffs%quartet_cutoff) cycle

            pq = p(:,p12) - q(:,p34)
            pq2 = dot_product(pq,pq)
            gammapq = gammap(p12)*gammaq(p34)*gpqinv

!           Boys function calculation
            tt = PQ2*gammapq

            if (tt > tmax) then
!             Case of a large argument value (asymptotic expansion + forward recursion)
              t2m1 = 0.5d0/tt
              F(1) = halfsqrtpi/sqrt(tt)
              do m = 1, am_tot
                F(m+1) = (2*m-1)*F(m)*t2m1
              end do

            else
!             Case of a small argument (interpolation + backward recursion)
              call boysf_nonasym(am_tot, tt, F)
            end if

            p1234 = p1234 + 1

#if LIBINT2_DEFINED_PA_x
            erieval(p1234)%PA_x(1) = P(1,p12) - A(1)
#endif
#if LIBINT2_DEFINED_PA_y
            erieval(p1234)%PA_y(1) = P(2,p12) - A(2)
#endif
#if LIBINT2_DEFINED_PA_z
            erieval(p1234)%PA_z(1) = P(3,p12) - A(3)
#endif
#if LIBINT2_DEFINED_AB_x
            erieval(p1234)%AB_x(1) = A(1) - B(1)
#endif
#if LIBINT2_DEFINED_AB_y
            erieval(p1234)%AB_y(1) = A(2) - B(2)
#endif
#if LIBINT2_DEFINED_AB_z
            erieval(p1234)%AB_z(1) = A(3) - B(3)
#endif
#if LIBINT2_DEFINED_oo2z
            erieval(p1234)%oo2z(1) = 0.5_dp*gpinv(p12)
#endif
#if LIBINT2_DEFINED_QC_x
            erieval(p1234)%QC_x(1) = Q(1,p34) - C(1)
#endif
#if LIBINT2_DEFINED_QC_y
            erieval(p1234)%QC_y(1) = Q(2,p34) - C(2)
#endif
#if LIBINT2_DEFINED_QC_z
            erieval(p1234)%QC_z(1) = Q(3,p34) - C(3)
#endif
#if LIBINT2_DEFINED_CD_x
            erieval(p1234)%CD_x(1) = C(1) - D(1)
#endif
#if LIBINT2_DEFINED_CD_y
            erieval(p1234)%CD_y(1) = C(2) - D(2)
#endif
#if LIBINT2_DEFINED_CD_z
            erieval(p1234)%CD_z(1) = C(3) - D(3)
#endif
#if LIBINT2_DEFINED_oo2e
            erieval(p1234)%oo2e(1) = 0.5_dp*gqinv(p34)
#endif
#if LIBINT2_DEFINED_WP_x
            erieval(p1234)%WP_x(1) = -gammaq(p34) * gpqinv * pq(1)
#endif
#if LIBINT2_DEFINED_WP_y
            erieval(p1234)%WP_y(1) = -gammaq(p34) * gpqinv * pq(2)
#endif
#if LIBINT2_DEFINED_WP_z
            erieval(p1234)%WP_z(1) = -gammaq(p34) * gpqinv * pq(3)
#endif
#if LIBINT2_DEFINED_WQ_x
            erieval(p1234)%WQ_x(1) = gammap(p12) * gpqinv * pq(1)
#endif
#if LIBINT2_DEFINED_WQ_y
            erieval(p1234)%WQ_y(1) = gammap(p12) * gpqinv * pq(2)
#endif
#if LIBINT2_DEFINED_WQ_z
            erieval(p1234)%WQ_z(1) = gammap(p12) * gpqinv * pq(3)
#endif
#if LIBINT2_DEFINED_oo2ze
            erieval(p1234)%oo2ze(1) = 0.5_dp*gpqinv
#endif

#if LIBINT2_DEFINED_roz
            erieval(p1234)%roz(1) = gammapq*gpinv(p12)
#endif
#if LIBINT2_DEFINED_roe
            erieval(p1234)%roe(1) = gammapq*gqinv(p34)
#endif
            IF (deriv_order > 0) THEN
#if LIBINT2_DEFINED_alpha1rho_over_zeta2
              erieval(p1234)%alpha1rho_over_zeta2(1) = alpha1(p12)*gammapq*gpinv(p12)*gpinv(p12)
#endif
#if LIBINT2_DEFINED_alpha2rho_over_zeta2
              erieval(p1234)%alpha2rho_over_zeta2(1) = alpha2(p12)*gammapq*gpinv(p12)*gpinv(p12)
#endif
#if LIBINT2_DEFINED_alpha3rho_over_eta2
              erieval(p1234)%alpha3rho_over_eta2(1) = alpha3(p34)*gammapq*gqinv(p34)*gqinv(p34)
#endif
#if LIBINT2_DEFINED_alpha4rho_over_eta2
              erieval(p1234)%alpha4rho_over_eta2(1) = alpha4(p34)*gammapq*gqinv(p34)*gqinv(p34)
#endif
#if LIBINT2_DEFINED_alpha1over_zetapluseta
              erieval(p1234)%alpha1over_zetapluseta(1) = alpha1(p12)*gpqinv
#endif
#if LIBINT2_DEFINED_alpha2over_zetapluseta
              erieval(p1234)%alpha2over_zetapluseta(1) = alpha2(p12)*gpqinv
#endif
#if LIBINT2_DEFINED_alpha3over_zetapluseta
              erieval(p1234)%alpha3over_zetapluseta(1) = alpha3(p34)*gpqinv
#endif
#if LIBINT2_DEFINED_alpha4over_zetapluseta
              erieval(p1234)%alpha4over_zetapluseta(1) = alpha4(p34)*gpqinv
#endif
#if LIBINT2_DEFINED_rho12_over_alpha1
              erieval(p1234)%rho12_over_alpha1(1) = alpha2(p12)*gpinv(p12)
#endif
#if LIBINT2_DEFINED_rho12_over_alpha2
              erieval(p1234)%rho12_over_alpha2(1) = alpha1(p12)*gpinv(p12)
#endif
              rhoq = alpha3(p34)*alpha4(p34)*gqinv(p34)
#if LIBINT2_DEFINED_rho34_over_alpha3
              erieval(p1234)%rho34_over_alpha3(1) = alpha4(p34)*gqinv(p34)
#endif
#if LIBINT2_DEFINED_rho34_over_alpha4
              erieval(p1234)%rho34_over_alpha4(1) = alpha3(p34)*gqinv(p34)
#endif
#if LIBINT2_DEFINED_two_alpha0_bra
              erieval(p1234)%two_alpha0_bra(1) = 2.0_dp*alpha1(p12)
#endif
#if LIBINT2_DEFINED_two_alpha0_ket
              erieval(p1234)%two_alpha0_ket(1) = 2.0_dp*alpha2(p12)
#endif
#if LIBINT2_DEFINED_two_alpha1ket
              erieval(p1234)%two_alpha1ket(1) = 2.0_dp*alpha4(p34)
#endif
#if LIBINT2_DEFINED_two_alpha1bra
              erieval(p1234)%two_alpha1bra(1) = 2.0_dp*alpha3(p34)
#endif
            END IF
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0
            IF (am_tot >= 0) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0(1) = &
              pfac*F(1)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1
            IF (am_tot >= 1) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1(1) = &
              pfac*F(2)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2
            IF (am_tot >= 2) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2(1) = &
              pfac*F(3)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3
            IF (am_tot >= 3) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3(1) = &
              pfac*F(4)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4
            IF (am_tot >= 4) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4(1) = &
              pfac*F(5)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5
            IF (am_tot >= 5) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5(1) = &
              pfac*F(6)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6
            IF (am_tot >= 6) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6(1) = &
              pfac*F(7)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7
            IF (am_tot >= 7) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7(1) = &
              pfac*F(8)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8
            IF (am_tot >= 8) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8(1) = &
              pfac*F(9)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9
            IF (am_tot >= 9) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9(1) = &
              pfac*F(10)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10
            IF (am_tot >= 10) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10(1) = &
              pfac*F(11)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11
            IF (am_tot >= 11) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11(1) = &
              pfac*F(12)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12
            IF (am_tot >= 12) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12(1) = &
              pfac*F(13)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13
            IF (am_tot >= 13) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13(1) = &
              pfac*F(14)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14
            IF (am_tot >= 14) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14(1) = &
              pfac*F(15)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15
            IF (am_tot >= 15) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15(1) = &
              pfac*F(16)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16
            IF (am_tot >= 16) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16(1) = &
              pfac*F(17)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17
            IF (am_tot >= 17) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17(1) = &
              pfac*F(18)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18
            IF (am_tot >= 18) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18(1) = &
              pfac*F(19)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19
            IF (am_tot >= 19) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19(1) = &
              pfac*F(20)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20
            IF (am_tot >= 20) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20(1) = &
              pfac*F(21)
#endif
         END DO
      END DO
      end associate

#if LIBINT2_CONTRACTED_INTS
      erieval(1)%contrdepth = int(p1234, C_INT)
#endif

      select case (deriv_order)
      case (0)
          CALL C_F_PROCPOINTER(libint2_build_eri (am(4), am(3), am(2), am(1)), build_eri)
#if INCLUDE_ERI >= 1
      case (1)
          CALL C_F_PROCPOINTER(libint2_build_eri1(am(4), am(3), am(2), am(1)), build_eri)
#endif
#if INCLUDE_ERI >= 2
      case (2)
          CALL C_F_PROCPOINTER(libint2_build_eri2(am(4), am(3), am(2), am(1)), build_eri)
#endif
      case default
          error stop "Unsupported `deriv_order` ! Please reconfigure and rebuild libint"
      end select

      if (p1234 == 0) return
      CALL build_eri(erieval)
      zero_shq = .false.
#endif
   end subroutine

   subroutine libint_to_ghondo(basis, shell_ids, deriv_order, erieval, flips, ghondo)
     use constants, only: shells_pnrm2
     use basis_tools, only: basis_set
     type(basis_set) :: basis
     integer, intent(in) :: shell_ids(4), deriv_order, flips(4)
     type(libint_t), intent(in) :: erieval(*)
     real(kind=dp), contiguous, intent(inout) :: ghondo(:,:,:,:)
     real(kind=dp), dimension(:,:,:,:), pointer :: eri_shell_set

     integer :: i
     integer :: n(4), n0(4), am0(4), ids(4), i_target, na, nb, nc, nd
     real(kind=dp), dimension(28) :: pnorma, pnormb, pnormc, pnormd
     real(kind=dp) :: scaleab, scalecd

     integer, parameter, dimension(3) :: n_targets = [1, 12, 78]
     integer, parameter :: angs(0:7) = [ ( (i+1)*(i+2)/2, i = 0, 7)]

     am0 = basis%ktype(shell_ids)-1
     n0 = angs(am0)
     n = n0(flips)
     pnorma = shells_pnrm2(:n0(1),am0(1))
     pnormb = shells_pnrm2(:n0(2),am0(2))
     pnormc = shells_pnrm2(:n0(3),am0(3))
     pnormd = shells_pnrm2(:n0(4),am0(4))

     do i_target = 1, n_targets(deriv_order + 1)
       call c_f_pointer(erieval(1)%targets(i_target), eri_shell_set, &
           shape=n([4,3,2,1]))

       do na = 1, n0(1)
         do nb = 1, n0(2)
           scaleab =  pnorma(na) * pnormb(nb)
           do nc = 1, n0(3)
             do nd = 1, n0(4)
               ids = [na, nb, nc, nd]
               ids = ids(flips)
               scalecd =  pnormc(nc) * pnormd(nd)
               ghondo(nd,nc,nb,na) = &
                 eri_shell_set(ids(4),ids(3),ids(2),ids(1)) * scaleab * scalecd
             end do

           end do

         end do
       end do

     end do

   end subroutine

   subroutine libint_print_eri(basis, shell_ids, deriv_order, erieval, flips)
     use constants, only: shells_pnrm2
     use basis_tools, only: basis_set
     type(basis_set) :: basis
     integer, intent(in) :: shell_ids(4), deriv_order, flips(4)
     type(libint_t), dimension(*), intent(in) :: erieval
     real(kind=dp), dimension(:,:,:,:), pointer :: eri_shell_set
     integer :: n(4), n0(4), i_target, na, nb, nc, nd, ishell
     integer :: am(4), am0(4)
     integer, parameter, dimension(3) :: n_targets = [1, 12, 78]
     integer :: mins(4), maxs(4), ids(4), shls(4)
     real(kind=dp) :: pnorms(36,4)

     shls = shell_ids(flips)

     am = basis%ktype(shls)-1
     n = (am + 1)*(am + 2)/2

     mins = shmax(am-1)+1
     maxs = shmax(am)

     pnorms(:,1) = shells_pnrm2(:n(1),am(1))
     pnorms(:,2) = shells_pnrm2(:n(2),am(2))
     pnorms(:,3) = shells_pnrm2(:n(3),am(3))
     pnorms(:,4) = shells_pnrm2(:n(4),am(4))

     am0 = basis%ktype(shell_ids)-1
     n0 = (am0 + 1)*(am0 + 2)/2

     do i_target = 1, n_targets(deriv_order + 1)
       call c_f_pointer(erieval(1)%targets(i_target), eri_shell_set, &
           shape=n([4,3,2,1]))

       ishell = 0
       do na = 1, n0(1)
       do nb = 1, n0(2)
       do nc = 1, n0(3)
       do nd = 1, n0(4)
         ids = [na, nb, nc, nd]
         ids = ids(flips)
         write (*, "(a6, 2i3, a, 2i3, a, es30.15)") &
           "elem (", na, nb, " |", nc, nd,  ") = ", &
           eri_shell_set(ids(4),ids(3),ids(2),ids(1)) * &
           pnorms(ids(1),1)*pnorms(ids(2),2)*&
           pnorms(ids(3),3)*pnorms(ids(4),4)
       end do
       end do
       end do
       end do
     end do

   contains

      integer elemental function shmax(i)
        integer, intent(in) :: i
        shmax = (i+1)*(i+2)*(i+3)/6
      end function

   end subroutine

!> @brief Non-asymptotic Boys function case
  pure subroutine boysf_nonasym(n, tt, ft)

    use boys_lut, only: rxinc, rfinc, nord, igrid, irgrd, fgrid, rmr
    implicit none

    real(kind=8), intent(in) :: tt
    integer, intent(in) :: n

    real(kind=8), intent(out) :: ft(0:*)

    real(kind=8) :: tv, tx, fx, et, t2
    integer :: m, ip, ifxgrd, iftgrd

    ifxgrd = igrid(n)
    iftgrd = irgrd(n)

    tv = tt*rfinc(ifxgrd)
    tx = tt*rxinc
    ip = nint(tv)

    fx = 0
    et = 0
    et = exp(-tt)
    do m = nord, 0, -1
      fx = (fx*tv+fgrid(m, ip, ifxgrd))
    end do

    ft(iftgrd) = fx

    t2 = 2*tt
    do m = iftgrd, 1, -1
      ft(m-1) = (t2*ft(m)+et)*rmr(m)
    end do

  end subroutine boysf_nonasym

end module
