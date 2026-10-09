!> @file mp2_gradient.F90
!>
!> @brief Analytic nuclear gradient for canonical closed-shell MP2 and its
!>        spin-component-scaled variants.
!>
!> The implementation is a stationary MP2 Lagrangian derivative.  It builds
!> the unrelaxed MP2 one-particle density and the integral (two-particle)
!> density from the doubles amplitudes, solves the RHF orbital-response
!> equation with the native CPHF solver, and contracts the resulting relaxed
!> one- and two-particle densities with the existing grd1/grd2 derivative
!> integral engines.  No finite differences enter this routine.
!>
!> The first implementation deliberately supports RHF only.  UHF needs the
!> coupled alpha/beta density blocks and ROHF additionally needs the derivative
!> of the nonzero second-order singles term used by mp2_lib.  Refusing those
!> references is safer than returning a restricted formula for an open-shell
!> energy.
!>
!> Memory note: the original code allocated a full n^4 MO 2-RDM (g2mo) and
!> then transformed it to the AO basis via mo2ao_4index which required 3
!> additional n^4 temporary arrays.  The current implementation skips the MO
!> 2-RDM entirely: the orbital-Lagrangian intermediate I_mat is formed
!> directly from the l2 amplitudes (saving one n^4), and the AO 2-RDM is
!> built from l2 with O(n^3) scratch arrays (saving 3 more n^4 temps).
!> Peak workspace is ~3 n^4 (eri_mo + erip + g2c) instead of ~6 n^4.
!> DEFAULT_MAX_NBF raised from 60 to 100.
module mp2_gradient_mod

  use precision, only: dp
  use types, only: information
  use basis_tools, only: basis_set, bas_norm_matrix, build_cart_density
  use constants, only: HARMONIC_ACTIVE, NUM_CART_BF
  use grd2, only: grd2_driver, grd2_compute_data_t

  implicit none
  private

  character(len=*), parameter :: module_name = 'mp2_gradient_mod'
  integer, parameter :: DEFAULT_MAX_NBF = 100

  public :: mp2_gradient

  ! The derivative-ERI consumer.  d0 is the converged RHF spin-summed density,
  ! dc is the relaxed MP2 correction to the spin-summed density, and g2c is the
  ! MP2 transition two-particle density, symmetrized in the AO basis.  grd2's
  ! unique-quartet convention consumes four times the physical symmetrized
  ! 2-RDM, hence the factor used in get_density below.
  type, extends(grd2_compute_data_t) :: grd2_mp2_compute_data_t
    real(dp), allocatable :: d0(:,:), dc(:,:), g2c(:,:,:,:)
    real(dp), allocatable :: d0_cart(:,:), dc_cart(:,:), g2c_cart(:,:,:,:)
    integer, allocatable :: cart_off(:)
    integer :: nbf = 0
    integer :: nbf_cart = 0
  contains
    procedure :: init => grd2_mp2_init
    procedure :: clean => grd2_mp2_clean
    procedure :: get_density => grd2_mp2_get_density
    procedure :: build_cart => grd2_mp2_build_cart
  end type grd2_mp2_compute_data_t

contains

!###############################################################################

  subroutine mp2_gradient_C(c_handle) bind(C, name='mp2_gradient')
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call mp2_gradient(inf)
  end subroutine mp2_gradient_C

!###############################################################################

  subroutine mp2_gradient(infos)
    use io_constants, only: iw
    use messages, only: show_message, with_abort
    use printing, only: print_module_info
    use grd1, only: eijden, print_gradient, grad_nn, grad_ee_overlap, &
                    grad_ee_kinetic, grad_en_hellman_feynman, &
                    grad_en_pulay, grad_1e_ecp
    use constants, only: tol_int
    use mathlib, only: unpack_matrix, pack_matrix
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_VEC_MO_A, &
                                   OQP_E_MO_A, OQP_AO_ERI
    use fci_integrals_mod, only: fci_ao_integrals
    use mo_transform_mod, only: mo_transform_eri
    use cphf_mod, only: cphf_solve
    use iso_c_binding, only: c_int32_t, c_int64_t

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis
    type(grd2_mp2_compute_data_t) :: gcomp

    real(dp), contiguous, pointer :: mo(:,:), eps(:), dm0p(:), erip(:)
    real(dp), allocatable :: mo_c(:,:), eri_mo(:), t2(:,:,:,:), l2(:,:,:,:)
    real(dp), allocatable :: doo(:,:), dvv(:,:), dm1mo(:,:), imat(:,:)
    real(dp), allocatable :: dm0(:,:), dc(:,:), vmat(:,:), xvo(:,:), rhs(:,:), u(:,:)
    real(dp), allocatable :: zeta(:,:), imao(:,:), wao(:,:), pocc(:,:), tmp(:,:)
    real(dp), allocatable :: dtotal(:,:), wpack(:), dpack(:), de2(:,:)
    integer :: n, no, nv, i, j, a, b, p, q, ia, ierr, max_nbf, ln
    integer(c_int64_t) :: trc
    real(dp) :: den, os_scale, ss_scale, tol, e_dense, e_expected, energy_tol
    character(len=32) :: env
    logical :: cphf_converged

    if (infos%control%scftype /= 1) then
      call show_message('MP2 analytic gradient currently supports RHF references only', with_abort)
    end if
    if (infos%control%hamilton >= 20) then
      call show_message('MP2 analytic gradient requires a Hartree-Fock reference', with_abort)
    end if

    basis => infos%basis
    basis%atoms => infos%atoms
    n = basis%nbf
    no = infos%mol_prop%nocc
    nv = n - no
    if (no <= 0 .or. nv <= 0) &
      call show_message('MP2 analytic gradient requires occupied and virtual orbitals', with_abort)

    max_nbf = DEFAULT_MAX_NBF
    call get_environment_variable('OQP_MP2_GRAD_MAX_NBF', env, ln)
    if (ln > 0) read(env, *, iostat=ierr) max_nbf
    if (n > max_nbf) then
      call show_message('MP2 analytic gradient dense 2-RDM guard exceeded; set OQP_MP2_GRAD_MAX_NBF explicitly', &
                        with_abort)
    end if

    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, eps)
    call tagarray_get_data(infos%dat, OQP_DM_A, dm0p)
    call fci_ao_integrals(infos)
    call tagarray_get_data(infos%dat, OQP_AO_ERI, erip)

    allocate(mo_c(0:n-1,0:n-1), eri_mo(n**4), t2(no,no,nv,nv), l2(no,no,nv,nv), &
             doo(no,no), dvv(nv,nv), dm1mo(n,n), imat(n,n), &
             dm0(n,n), dc(n,n), vmat(n,n), xvo(nv,no), rhs(no*nv,1), &
             u(no*nv,1), zeta(n,n), imao(n,n), wao(n,n), pocc(n,n), &
             tmp(n,n), dtotal(n,n), wpack(n*(n+1)/2), dpack(n*(n+1)/2), &
             source=0.0_dp, stat=ierr)
    if (ierr /= 0) call show_message('MP2 analytic gradient allocation failed', with_abort)

    mo_c = transpose(mo)
    trc = mo_transform_eri(int(n, c_int32_t), erip, mo_c, eri_mo)
    if (trc /= 0_c_int64_t) &
      call show_message('MP2 analytic gradient AO-to-MO ERI transformation failed', with_abort)

    os_scale = infos%dft%MP2OS_Scale
    ss_scale = infos%dft%MP2SS_Scale

    ! Canonical RHF doubles and the derivative (lambda) amplitude appropriate
    ! to E = cOS E_OS + cSS E_SS.  For ordinary MP2 l2 = 2*t2-t2(ab<->ba).
    do i = 1, no
      do j = 1, no
        do a = 1, nv
          do b = 1, nv
            den = eps(i) + eps(j) - eps(no+a) - eps(no+b)
            if (abs(den) < 1.0e-10_dp) cycle
            t2(i,j,a,b) = eri4(eri_mo, n, i, no+a, j, no+b) / den
          end do
        end do
      end do
    end do
    e_dense = 0.0_dp
    do i = 1, no
      do j = 1, no
        do a = 1, nv
          do b = 1, nv
            l2(i,j,a,b) = (os_scale + ss_scale)*t2(i,j,a,b) &
                         - ss_scale*t2(i,j,b,a)
            e_dense = e_dense + l2(i,j,a,b) &
                                * eri4(eri_mo, n, i, no+a, j, no+b)
          end do
        end do
      end do
    end do
    e_expected = infos%mol_energy%energy - (infos%mol_energy%ehf1 &
                                            + infos%mol_energy%vee &
                                            + infos%mol_energy%nenergy)
    energy_tol = max(1.0e-8_dp, 1.0e-7_dp*abs(e_expected))
    if (abs(e_dense - e_expected) > energy_tol) then
      call show_message('(A,ES24.15)', 'MP2 energy from dense transformed integrals = ', e_dense)
      call show_message('(A,ES24.15)', 'MP2 energy from the production energy kernel = ', e_expected)
      call show_message('MP2 analytic gradient integral transformation failed its energy check', with_abort)
    end if

    ! Unrelaxed MP2 1-RDM blocks.
    do p = 1, no
      do q = 1, no
        do i = 1, no
          do a = 1, nv
            do b = 1, nv
              doo(p,q) = doo(p,q) - l2(i,p,a,b)*t2(i,q,a,b)
            end do
          end do
        end do
      end do
    end do
    do p = 1, nv
      do q = 1, nv
        do i = 1, no
          do j = 1, no
            do a = 1, nv
              dvv(p,q) = dvv(p,q) + l2(i,j,a,p)*t2(i,j,a,q)
            end do
          end do
        end do
      end do
    end do
    dm1mo(1:no,1:no) = doo + transpose(doo)
    dm1mo(no+1:n,no+1:n) = dvv + transpose(dvv)

    ! The amplitude-linear 2-RDM has only ovov and vovo blocks.  Form the
    ! orbital Lagrangian intermediate I_mat directly from the l2 amplitudes,
    ! avoiding a full n^4 MO 2-RDM intermediate.
    do p = 1, n
      do i = 1, no
        do j = 1, no
          do a = 1, nv
            do b = 1, nv
              imat(p,no+a) = imat(p,no+a) &
                - 2.0_dp * eri4(eri_mo, n, i, p, j, no+b) * l2(i,j,a,b)
              imat(p,i) = imat(p,i) &
                - 2.0_dp * eri4(eri_mo, n, no+a, p, no+b, j) * l2(i,j,a,b)
            end do
          end do
        end do
      end do
    end do

    ! AO correction density and its closed-shell response Fock, V=2J-K.
    call mo2ao_matrix(n, mo, dm1mo, dc)
    call dense_rhf_response(n, erip, dc, vmat)

    do a = 1, nv
      do i = 1, no
        xvo(a,i) = dot_product(mo(:,no+a), matmul(vmat,mo(:,i))) &
                   + imat(i,no+a) - imat(no+a,i)
        ia = (a-1)*no + i
        rhs(ia,1) = -xvo(a,i)
      end do
    end do
    open(unit=iw, file=infos%log_filename, position='append')
    call cphf_solve(infos, 1, rhs, u, tol=1.0e-10_dp, &
                    maxit=max(100,no*nv+10), converged=cphf_converged)
    if (.not. cphf_converged) then
      call show_message('MP2 analytic gradient CPHF response did not converge', with_abort)
    end if
    do a = 1, nv
      do i = 1, no
        ia = (a-1)*no + i
        dm1mo(no+a,i) = dm1mo(no+a,i) + u(ia,1)
        dm1mo(i,no+a) = dm1mo(i,no+a) + u(ia,1)
      end do
    end do
    call mo2ao_matrix(n, mo, dm1mo, dc)

    ! Energy-weighted density and the remaining overlap terms.
    do p = 1, n
      do q = 1, n
        zeta(p,q) = 0.5_dp*(eps(p)+eps(q))*dm1mo(p,q)
      end do
    end do
    do a = 1, nv
      do i = 1, no
        zeta(no+a,i) = eps(i)*dm1mo(no+a,i)
        zeta(i,no+a) = eps(i)*dm1mo(i,no+a)
        imat(no+a,i) = imat(i,no+a)
      end do
    end do
    call mo2ao_matrix(n, mo, zeta, wao)
    call mo2ao_matrix(n, mo, imat, imao)
    call dense_rhf_response(n, erip, dc, vmat)
    pocc = matmul(mo(:,1:no), transpose(mo(:,1:no)))
    tmp = matmul(pocc, vmat)
    vmat = matmul(tmp, pocc)
    wao = 0.5_dp*(imao + transpose(imao)) - wao - vmat

    call unpack_matrix(dm0p, dm0, 'U')
    dtotal = dm0 + dc
    call pack_matrix(dtotal, dpack, 'U')
    call eijden(wpack, n, infos)
    call add_packed_symmetric(n, wao, wpack)

    ! Build the AO 2-RDM directly from the compact l2 amplitudes, avoiding
    ! the full n^4 MO intermediate and 3 n^4 transformation temps.
    allocate(gcomp%g2c(n,n,n,n), gcomp%d0(n,n), gcomp%dc(n,n), stat=ierr)
    if (ierr /= 0) call show_message('MP2 gradient 2-RDM allocation failed', with_abort)
    call l2_to_ao_2rdm(n, no, mo, l2, gcomp%g2c)
    gcomp%d0 = dm0
    gcomp%dc = dc
    gcomp%nbf = n

    call print_module_info('MP2_Gradient', 'Analytic RHF-MP2 Nuclear Gradient')
    write(iw,'(3X,A,F10.6)') 'same-spin scale     = ', ss_scale
    write(iw,'(3X,A,F10.6)') 'opposite-spin scale = ', os_scale
    write(iw,'(3X,A,I0)') 'dense AO 2-RDM dimension = ', n

    tol = tol_int*log(10.0_dp)
    associate(grad => infos%atoms%grad, xyz => infos%atoms%xyz, &
              zn => infos%atoms%zn - infos%basis%ecp_zn_num)
      grad = 0.0_dp
      call grad_nn(infos%atoms, infos%basis%ecp_zn_num)
      call grad_ee_overlap(basis, wpack, grad, logtol=tol)
      call grad_ee_kinetic(basis, dpack, grad, logtol=tol)
      call grad_en_hellman_feynman(basis, xyz, zn, dpack, grad, logtol=tol)
      call grad_en_pulay(basis, xyz, zn, dpack, grad, logtol=tol)
      call grad_1e_ecp(infos, basis, xyz, dpack, grad, logtol=tol)
    end associate

    allocate(de2(3,ubound(infos%atoms%zn,1)), source=0.0_dp)
    gcomp%coulscale = 1.0_dp
    gcomp%hfscale = 1.0_dp
    call gcomp%init()
    call gcomp%build_cart(basis)
    call grd2_driver(infos, basis, de2, gcomp, petite=.true.)
    infos%atoms%grad = infos%atoms%grad + de2
    call print_gradient(infos)
    close(iw)

    call gcomp%clean()
  end subroutine mp2_gradient

!###############################################################################

  pure real(dp) function eri4(g, n, i, j, k, l) result(v)
    real(dp), intent(in) :: g(:)
    integer, intent(in) :: n, i, j, k, l
    v = g(i + (j-1)*n + (k-1)*n*n + (l-1)*n*n*n)
  end function eri4

  subroutine mo2ao_matrix(n, c, mmo, mao)
    integer, intent(in) :: n
    real(dp), intent(in) :: c(n,n), mmo(n,n)
    real(dp), intent(out) :: mao(n,n)
    mao = matmul(c, matmul(mmo, transpose(c)))
  end subroutine mo2ao_matrix

  subroutine dense_rhf_response(n, eri, d, v)
    integer, intent(in) :: n
    real(dp), intent(in) :: eri(:), d(n,n)
    real(dp), intent(out) :: v(n,n)
    integer :: p, q, r, s
    v = 0.0_dp
    do p = 1, n
      do q = 1, n
        do r = 1, n
          do s = 1, n
            v(p,q) = v(p,q) + d(r,s)*(2.0_dp*eri4(eri,n,p,q,r,s) &
                                             - eri4(eri,n,p,r,q,s))
          end do
        end do
      end do
    end do
  end subroutine dense_rhf_response

  subroutine add_packed_symmetric(n, a, ap)
    use mathlib, only: pack_matrix
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n,n)
    real(dp), intent(inout) :: ap(:)
    real(dp), allocatable :: p(:), s(:,:)
    allocate(p(n*(n+1)/2), s(n,n))
    s = 0.5_dp*(a + transpose(a))
    call pack_matrix(s, p, 'U')
    ap = ap + p
  end subroutine add_packed_symmetric

  !> Build the AO 2-RDM directly from the compact l2 tensor, avoiding
  !> the full n^4 MO intermediate and the 3 n^4 scratch arrays that
  !> mo2ao_4index required.
  !>
  !> The nonzero MO 2-RDM blocks are:
  !>   ovov: g2mo(i,no+a,j,no+b) = 2*l2(i,j,a,b)
  !>   vovo: g2mo(no+a,i,no+b,j) = 2*l2(i,j,a,b)
  !>
  !> The four-index MO->AO transform gao(mu,nu,la,si) =
  !>   sum_{pqrs} C(mu,p) C(nu,q) C(la,r) C(si,s) * g2mo(p,q,r,s)
  !> is computed by batching over one occupied index at a time with
  !> O(n^3) working storage.  After filling gao the eight-fold
  !> permutational symmetry is restored in-place.
  subroutine l2_to_ao_2rdm(n, no, c, l2, gao)
    integer, intent(in) :: n, no
    real(dp), intent(in) :: c(n,n), l2(no,no,n-no,n-no)
    real(dp), intent(inout) :: gao(n,n,n,n)
    integer :: nv, j, a, b, mu, nu, la, si, i
    real(dp), allocatable :: w1(:,:,:), w2(:,:,:), w3(:,:,:)

    nv = n - no
    allocate(w1(n,nv,nv), w2(n,n,nv), w3(n,n,n))
    gao = 0.0_dp

    ! === ovov contribution ===
    ! gao(mu,nu,la,si) += sum_{i,j,a,b} 2*l2(i,j,a,b)
    !                      * C(mu,i) * C(nu,no+a) * C(la,j) * C(si,no+b)
    do j = 1, no
      ! w1(mu,a,b) = sum_i C(mu,i) * l2(i,j,a,b)
      do b = 1, nv
        do a = 1, nv
          w1(:,a,b) = 0.0_dp
          do i = 1, no
            do mu = 1, n
              w1(mu,a,b) = w1(mu,a,b) + c(mu,i) * l2(i,j,a,b)
            end do
          end do
        end do
      end do

      ! w2(mu,nu,b) = sum_a w1(mu,a,b) * C(nu,no+a)
      w2 = 0.0_dp
      do b = 1, nv
        do a = 1, nv
          do nu = 1, n
            do mu = 1, n
              w2(mu,nu,b) = w2(mu,nu,b) + w1(mu,a,b) * c(nu,no+a)
            end do
          end do
        end do
      end do

      ! w3(mu,nu,si) = sum_b w2(mu,nu,b) * C(si,no+b)
      w3 = 0.0_dp
      do b = 1, nv
        do si = 1, n
          do nu = 1, n
            do mu = 1, n
              w3(mu,nu,si) = w3(mu,nu,si) + w2(mu,nu,b) * c(si,no+b)
            end do
          end do
        end do
      end do

      ! gao(mu,nu,la,si) += 2*C(la,j)*w3(mu,nu,si)
      do si = 1, n
        do la = 1, n
          do nu = 1, n
            do mu = 1, n
              gao(mu,nu,la,si) = gao(mu,nu,la,si) &
                + 2.0_dp * c(la,j) * w3(mu,nu,si)
            end do
          end do
        end do
      end do
    end do

    ! === vovo contribution ===
    ! gao(mu,nu,la,si) += sum_{i,j,a,b} 2*l2(i,j,a,b)
    !                      * C(mu,no+a) * C(nu,i) * C(la,no+b) * C(si,j)
    do j = 1, no
      ! w1(nu,a,b) = sum_i C(nu,i) * l2(i,j,a,b)
      do b = 1, nv
        do a = 1, nv
          w1(:,a,b) = 0.0_dp
          do i = 1, no
            do nu = 1, n
              w1(nu,a,b) = w1(nu,a,b) + c(nu,i) * l2(i,j,a,b)
            end do
          end do
        end do
      end do

      ! w2(mu,nu,b) = sum_a C(mu,no+a) * w1(nu,a,b)
      w2 = 0.0_dp
      do b = 1, nv
        do a = 1, nv
          do nu = 1, n
            do mu = 1, n
              w2(mu,nu,b) = w2(mu,nu,b) + c(mu,no+a) * w1(nu,a,b)
            end do
          end do
        end do
      end do

      ! w3(mu,nu,la) = sum_b C(la,no+b) * w2(mu,nu,b)
      w3 = 0.0_dp
      do b = 1, nv
        do la = 1, n
          do nu = 1, n
            do mu = 1, n
              w3(mu,nu,la) = w3(mu,nu,la) + w2(mu,nu,b) * c(la,no+b)
            end do
          end do
        end do
      end do

      ! gao(mu,nu,la,si) += 2*C(si,j)*w3(mu,nu,la)
      do si = 1, n
        do la = 1, n
          do nu = 1, n
            do mu = 1, n
              gao(mu,nu,la,si) = gao(mu,nu,la,si) &
                + 2.0_dp * c(si,j) * w3(mu,nu,la)
            end do
          end do
        end do
      end do
    end do

    ! Restore the eight-fold permutational symmetry of the AO 2-RDM.
    call symmetrize_eri_density(n, gao)
    deallocate(w1, w2, w3)
  end subroutine l2_to_ao_2rdm

  subroutine symmetrize_eri_density(n, g)
    integer, intent(in) :: n
    real(dp), intent(inout) :: g(n,n,n,n)
    real(dp), allocatable :: h(:,:,:,:)
    integer :: i,j,k,l
    allocate(h(n,n,n,n))
    h = g
    do i=1,n; do j=1,n; do k=1,n; do l=1,n
      g(i,j,k,l) = 0.125_dp*(h(i,j,k,l)+h(j,i,k,l)+h(i,j,l,k)+h(j,i,l,k) &
                              +h(k,l,i,j)+h(l,k,i,j)+h(k,l,j,i)+h(l,k,j,i))
    end do; end do; end do; end do
    deallocate(h)
  end subroutine symmetrize_eri_density

!###############################################################################

  subroutine grd2_mp2_init(this)
    class(grd2_mp2_compute_data_t), target, intent(inout) :: this
    this%nbf_cart = 0
  end subroutine grd2_mp2_init

  subroutine grd2_mp2_clean(this)
    class(grd2_mp2_compute_data_t), target, intent(inout) :: this
    if (allocated(this%d0)) deallocate(this%d0)
    if (allocated(this%dc)) deallocate(this%dc)
    if (allocated(this%g2c)) deallocate(this%g2c)
    if (allocated(this%d0_cart)) deallocate(this%d0_cart)
    if (allocated(this%dc_cart)) deallocate(this%dc_cart)
    if (allocated(this%g2c_cart)) deallocate(this%g2c_cart)
    if (allocated(this%cart_off)) deallocate(this%cart_off)
  end subroutine grd2_mp2_clean

  subroutine grd2_mp2_build_cart(this, basis)
    class(grd2_mp2_compute_data_t), intent(inout) :: this
    type(basis_set), intent(in) :: basis
    real(dp), allocatable :: tmp(:,:), expanded(:,:), pair12(:,:,:,:)
    integer, allocatable :: off(:)
    integer :: nc, p, q, mu, nu
    if (.not. HARMONIC_ACTIVE) return
    tmp = this%d0
    call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
    call build_cart_density(basis, tmp, this%d0_cart, this%cart_off, this%nbf_cart)
    tmp = this%dc
    call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
    call build_cart_density(basis, tmp, this%dc_cart, off, nc)
    allocate(pair12(this%nbf_cart,this%nbf_cart,this%nbf,this%nbf), &
             this%g2c_cart(this%nbf_cart,this%nbf_cart,this%nbf_cart,this%nbf_cart), &
             source=0.0_dp)
    do p = 1, this%nbf
      do q = 1, this%nbf
        tmp = this%g2c(:,:,p,q)
        call bas_norm_matrix(tmp, basis%bfnrm, basis%nbf)
        tmp = tmp*basis%bfnrm(p)*basis%bfnrm(q)
        call build_cart_density(basis, tmp, expanded, off, nc)
        pair12(:,:,p,q) = expanded
        deallocate(expanded)
      end do
    end do
    do mu = 1, this%nbf_cart
      do nu = 1, this%nbf_cart
        tmp = pair12(mu,nu,:,:)
        call build_cart_density(basis, tmp, expanded, off, nc)
        this%g2c_cart(mu,nu,:,:) = expanded
        deallocate(expanded)
      end do
    end do
    deallocate(pair12)
  end subroutine grd2_mp2_build_cart

  subroutine grd2_mp2_get_density(this, basis, id, dab, dabmax)
    use messages, only: show_message, with_abort
    class(grd2_mp2_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(dp), target, intent(out) :: dab(*)
    real(dp), intent(out) :: dabmax
    real(dp), pointer :: ab(:,:,:,:), d0(:,:), dc(:,:)
    real(dp), pointer :: g2(:,:,:,:)
    integer :: nb(4), loc(4), i,j,k,l,ii,jj,kk,ll
    real(dp) :: base, cross, corr, bfn

    if (HARMONIC_ACTIVE) then
      d0 => this%d0_cart; dc => this%dc_cart; g2 => this%g2c_cart
      loc = this%cart_off(id)-1
      nb = NUM_CART_BF(basis%am(id))
    else
      d0 => this%d0; dc => this%dc; g2 => this%g2c
      loc = basis%ao_offset(id)-1
      nb = basis%naos(id)
    end if
    ab(1:nb(4),1:nb(3),1:nb(2),1:nb(1)) => dab(1:product(nb))
    dabmax = 0.0_dp
    do i=1,nb(1); ii=loc(1)+i
      do j=1,nb(2); jj=loc(2)+j
        do k=1,nb(3); kk=loc(3)+k
          do l=1,nb(4); ll=loc(4)+l
            base = 4.0_dp*d0(ii,jj)*d0(kk,ll) &
                 - d0(ii,kk)*d0(jj,ll) - d0(ii,ll)*d0(jj,kk)
            cross = 4.0_dp*(dc(jj,ii)*d0(ll,kk) + dc(ll,kk)*d0(jj,ii)) &
                  - (dc(ii,kk)*d0(ll,jj) + dc(jj,kk)*d0(ll,ii) &
                    +dc(ll,ii)*d0(kk,jj) + dc(ll,jj)*d0(kk,ii))
            corr = 4.0_dp*g2(ii,jj,kk,ll)
            bfn = 1.0_dp
            if (.not. HARMONIC_ACTIVE) bfn = product(basis%bfnrm([ii,jj,kk,ll]))
            ab(l,k,j,i) = (base + cross + corr)*bfn
            dabmax = max(dabmax, abs(base+cross+corr))
          end do
        end do
      end do
    end do
  end subroutine grd2_mp2_get_density

end module mp2_gradient_mod
