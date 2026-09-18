!> @brief  Spin-restricted ensemble-referenced Kohn-Sham method, REKS(2,2),
!>         with state averaging (SA-REKS) and state interaction (SSR).
!>
!> @details
!>   Port of the GAMESS-US REKS(2,2) implementation by M. Filatov and S. Lee
!>   (source/reks.src, October 2021) to OpenQP.  Theory:
!>     M. Filatov, WIREs Comput. Mol. Sci. 5, 146 (2015);
!>     M. Filatov, F. Liu, T. J. Martinez, J. Chem. Phys. 147, 034113 (2017);
!>     M. Filatov, S. Lee, C. H. Choi, J. Chem. Theory Comput. 17, 5123 (2021).
!>
!>   Active space: two orbitals r (MO nocc) and s (MO nocc+1) sharing two
!>   electrons with fractional occupations n_r = 2 x, n_s = 2 (1-x).
!>   Four microstates (spatial orbitals shared, integer occupations):
!>     L=1 (core)^2 r^2, L=2 (core)^2 s^2, L=3 r(a) s(b), L=4 r(a) s(a).
!>   Microstates 3 and 4 stand for their spin-flipped partners as well and
!>   therefore carry a factor of 2 in every ensemble sum.
!>
!>   E_PPS = x E1 + (1-x) E2 - f(x) (E3 - E4),
!>   f(x)  = Y^(1 - (Y+delta)/(2(1+delta))),  Y = 4x(1-x), delta = 0.4
!>   E_OSS = 2 E3 - E4,   E_SA = w_PPS E_PPS + w_OSS E_OSS.
!>
!>   The orbitals are obtained from the coupling-operator (effective Fock)
!>   matrix assembled in the MO basis (GAMESS routine REXFM2FE).
!>
!> @author  Claude (port), from M. Filatov and S. Lee (GAMESS reks.src)
!> @date    September 2026
module reks

  use precision, only: dp
  use int2_compute, only: int2_fock_data_t, int2_storage_t, int2_compute_t
  use basis_tools, only: basis_set

  implicit none

  private
  public :: reks_driver
  public :: int2_reks_data_t
  public :: reks_cft, reks_fon_solve

  !> Two-electron consumer producing separate Coulomb and exchange matrices
  !> for every packed density supplied in `d(:,1:nd)`:
  !>   f(:,n)    = J[d_n]           (n = 1..nd)
  !>   f(:,nd+n) = c_x K[d_n]       (n = 1..nd), c_x = scale_exchange
  !> The Coulomb scale is applied to J.  Both are for per-spin densities.
  type, extends(int2_fock_data_t) :: int2_reks_data_t
  contains
    procedure :: parallel_start => int2_reks_data_t_parallel_start
    procedure :: update => int2_reks_data_t_update
  end type int2_reks_data_t

  integer, parameter :: nmic = 4
  !> microstate -> density block index (1=core, 2=core+r, 3=core+s, 4=core+r+s)
  integer, parameter :: mic_a(nmic) = [2, 3, 2, 4]
  integer, parameter :: mic_b(nmic) = [2, 3, 3, 1]
  !> occupation (0/1) of r and s in the alpha and beta parts of each microstate
  integer, parameter :: nra(nmic) = [1, 0, 1, 1]
  integer, parameter :: nsa(nmic) = [0, 1, 0, 1]
  integer, parameter :: nrb(nmic) = [1, 0, 0, 0]
  integer, parameter :: nsb(nmic) = [0, 1, 1, 0]
  !> spin-partner multiplicity of each microstate
  real(kind=dp), parameter :: fact(nmic) = [1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp]

contains

!###############################################################################
!  Two-electron consumer
!###############################################################################

  subroutine int2_reks_data_t_parallel_start(this, basis, nthreads)
    implicit none
    class(int2_reks_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nthreads

    this%fockdim = basis%nbf*(basis%nbf+1) / 2
    this%nfocks = 2*ubound(this%d, 2)

    call this%int2_fock_data_t_parallel_start(basis, nthreads)

  end subroutine int2_reks_data_t_parallel_start

!-------------------------------------------------------------------------------

  subroutine int2_reks_data_t_update(this, buf)
    implicit none
    class(int2_reks_data_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    integer :: ii, jj, kk, ll, ij, ik, il, jk, jl, kl, n, ii2, jj2, kk2
    integer :: nd, id, ix
    real(kind=dp) :: xval2, xval4, val, val1, val4
    real(kind=dp) :: aij, akl, aik, ajl, ail, ajk
    integer :: mythread

    nd = ubound(this%d, 2)
    xval2 = 2 * this%scale_exchange
    xval4 = 4 * this%scale_coulomb

    mythread = buf%thread_id
    if (this%atomic_fock) mythread = 1

    do n = 1, buf%ncur
      ii = buf%ids(1,n)
      jj = buf%ids(2,n)
      kk = buf%ids(3,n)
      ll = buf%ids(4,n)
      val = buf%ints(n)

      ii2 = ii*(ii-1)/2
      jj2 = jj*(jj-1)/2
      kk2 = kk*(kk-1)/2

      ij = ii2+jj
      ik = ii2+kk
      il = ii2+ll
      jk = jj2+kk
      jl = jj2+ll
      kl = kk2+ll
      if (jj<kk) jk = kk2 + jj
      if (jj<ll) jl = ll*(ll-1)/2 + jj

      val1 = val*xval2
      val4 = val*xval4

      do id = 1, nd
        ix = nd + id
        if (this%atomic_fock) then
          aij = val4*this%d(kl,id); akl = val4*this%d(ij,id)
          aik = val1*this%d(jl,id); ajl = val1*this%d(ik,id)
          ail = val1*this%d(jk,id); ajk = val1*this%d(il,id)
          !$omp atomic update
          this%f(ij,id,1) = this%f(ij,id,1) + aij
          !$omp atomic update
          this%f(kl,id,1) = this%f(kl,id,1) + akl
          !$omp atomic update
          this%f(ik,ix,1) = this%f(ik,ix,1) + aik
          !$omp atomic update
          this%f(jl,ix,1) = this%f(jl,ix,1) + ajl
          !$omp atomic update
          this%f(il,ix,1) = this%f(il,ix,1) + ail
          !$omp atomic update
          this%f(jk,ix,1) = this%f(jk,ix,1) + ajk
        else
          this%f(ij,id,mythread) = this%f(ij,id,mythread) + val4*this%d(kl,id)
          this%f(kl,id,mythread) = this%f(kl,id,mythread) + val4*this%d(ij,id)
          this%f(ik,ix,mythread) = this%f(ik,ix,mythread) + val1*this%d(jl,id)
          this%f(jl,ix,mythread) = this%f(jl,ix,mythread) + val1*this%d(ik,id)
          this%f(il,ix,mythread) = this%f(il,ix,mythread) + val1*this%d(jk,id)
          this%f(jk,ix,mythread) = this%f(jk,ix,mythread) + val1*this%d(il,id)
        end if
      end do
    end do

    buf%ncur = 0

  end subroutine int2_reks_data_t_update

!###############################################################################
!  REKS interpolating function and its derivatives with respect to x
!###############################################################################

  !> f(x) = Y^p(Y), Y = 4x(1-x), p = 1 - (Y+delta)/(2(1+delta))
  !> Returns f, df/dx, d2f/dx2.
  subroutine reks_cft(x, delta, f, f1, f2)
    implicit none
    real(kind=dp), intent(in) :: x, delta
    real(kind=dp), intent(out) :: f, f1, f2
    real(kind=dp) :: y, p, pp, lny, g, gy, gyy, yx, yxx
    ! GAMESS caps z = x(1-x) at 1e-10, i.e. Y at 4e-10
    real(kind=dp), parameter :: ymin = 4.0e-10_dp

    y = 4.0_dp*x*(1.0_dp-x)
    y = max(y, ymin)
    yx = 4.0_dp*(1.0_dp-2.0_dp*x)
    yxx = -8.0_dp

    p = 1.0_dp - (y+delta)/(2.0_dp*(1.0_dp+delta))
    pp = -1.0_dp/(2.0_dp*(1.0_dp+delta))
    lny = log(y)
    f = exp(p*lny)
    ! g = d ln f / dY
    g = pp*lny + p/y
    gy = f*g                      ! df/dY
    gyy = f*(g*g + 2.0_dp*pp/y - p/(y*y))   ! d2f/dY2
    f1 = gy*yx
    f2 = gyy*yx*yx + gy*yxx

  end subroutine reks_cft

!-------------------------------------------------------------------------------

  !> Minimize E(x) = x E1 + (1-x) E2 - f(x) (E3 - E4) over x in [0,1].
  !> GAMESS REXSOLVER/RexSlv4x4 logic: CI-type initial guess, Newton-Raphson
  !> with backtracking line search, then a safeguarded golden-section polish
  !> if Newton fails.
  subroutine reks_fon_solve(em, delta, x, ereks, iout, dbg)
    use io_constants, only: iw
    implicit none
    real(kind=dp), intent(in) :: em(nmic), delta
    real(kind=dp), intent(inout) :: x
    real(kind=dp), intent(out) :: ereks
    integer, intent(in) :: iout
    logical, intent(in) :: dbg

    real(kind=dp) :: e1, e2, dbc, x0, e0, e_new, grad, hess, dx, factor
    real(kind=dp) :: c1, d, disc, xn
    real(kind=dp), parameter :: cnvlim = 1.0e-12_dp, golden = 1.6180339887_dp
    real(kind=dp), parameter :: xlo = 0.0_dp, xhi = 1.0_dp
    integer :: iter, miter
    integer, parameter :: maxit = 20

    e1 = em(1); e2 = em(2); dbc = em(3) - em(4)

    if (abs(dbc) <= 1.0e-8_dp) then
      ! no coupling: pure single-configuration limit
      if (e2 < e1) then
        x = 0.0_dp
      else
        x = 1.0_dp
      end if
      ereks = min(e1, e2)
      return
    end if

    ! Initial guess from the 2x2 CI [[E1, Dbc],[Dbc, E2]] (GAMESS 4x4 with zeros)
    d = 0.5_dp*(e1 - e2)
    disc = sqrt(d*d + dbc*dbc)
    ! lowest eigenvector component on configuration 1
    if (abs(dbc) > 0.0_dp) then
      c1 = -dbc/sqrt(dbc*dbc + (d + disc)**2)
      x = c1*c1
    else
      x = 1.0_dp
    end if
    if (abs(x - 1.0_dp) <= 0.1_dp) x = 0.9_dp
    if (abs(x) <= 0.1_dp) x = 0.1_dp
    x = min(max(x, xlo), xhi)

    e0 = reks_efon(x)
    if (dbg) write(iout,'(2x,a,f12.8,a,f16.10)') 'FON solver: x0 =', x, ', E =', e0

    do iter = 1, maxit
      x0 = x
      call reks_grad_hess(x, grad, hess)
      if (abs(hess) < cnvlim) exit
      dx = -grad/hess
      if (hess < 0.0_dp) dx = -dx   ! move downhill if not convex
      factor = 1.0_dp
      miter = 0
      do
        xn = min(max(x0 + factor*dx, xlo), xhi)
        e_new = reks_efon(xn)
        if (e_new <= e0 .or. abs(factor*dx) < cnvlim) exit
        factor = factor/(golden*golden)
        miter = miter + 1
        if (miter > 40) exit
      end do
      x = xn
      e0 = reks_efon(x)
      if (dbg) write(iout,'(2x,i4,4x,f12.8,4x,es14.6,4x,f18.12)') iter, x, abs(x-x0), e0
      if (abs(x - x0) <= cnvlim) exit
    end do

    ! Safeguard: golden-section polish on [xlo,xhi] if the gradient is not small
    call reks_grad_hess(x, grad, hess)
    if (abs(grad) > 1.0e-6_dp) call reks_golden(x, e0)

    ereks = e0

  contains

    function reks_efon(xx) result(e)
      real(kind=dp), intent(in) :: xx
      real(kind=dp) :: e, f, f1, f2
      call reks_cft(xx, delta, f, f1, f2)
      e = xx*e1 + (1.0_dp-xx)*e2 - f*dbc
    end function reks_efon

    subroutine reks_grad_hess(xx, g, h)
      real(kind=dp), intent(in) :: xx
      real(kind=dp), intent(out) :: g, h
      real(kind=dp) :: f, f1, f2
      call reks_cft(xx, delta, f, f1, f2)
      g = e1 - e2 - f1*dbc
      h = -f2*dbc
    end subroutine reks_grad_hess

    subroutine reks_golden(xx, ee)
      real(kind=dp), intent(inout) :: xx, ee
      real(kind=dp) :: a, b, c, dd, fc, fd
      real(kind=dp), parameter :: gr = 0.6180339887_dp
      integer :: k
      a = xlo; b = xhi
      c = b - gr*(b-a); dd = a + gr*(b-a)
      fc = reks_efon(c); fd = reks_efon(dd)
      do k = 1, 200
        if (fc < fd) then
          b = dd; dd = c; fd = fc
          c = b - gr*(b-a); fc = reks_efon(c)
        else
          a = c; c = dd; fc = fd
          dd = a + gr*(b-a); fd = reks_efon(dd)
        end if
        if (abs(b-a) < 1.0e-13_dp) exit
      end do
      if (reks_efon(0.5_dp*(a+b)) < ee) then
        xx = 0.5_dp*(a+b)
        ee = reks_efon(xx)
      end if
    end subroutine reks_golden

  end subroutine reks_fon_solve

!###############################################################################
!  Ensemble coefficients (GAMESS REXCM)
!###############################################################################

  subroutine reks_coefficients(dnr, dns, delta, wpps, woss, cm)
    implicit none
    real(kind=dp), intent(in) :: dnr, dns, delta, wpps, woss
    real(kind=dp), intent(out) :: cm(nmic)
    real(kind=dp) :: f, f1, f2, ff
    logical :: statavg

    statavg = abs(wpps - 1.0_dp) > 1.0e-8_dp
    call reks_cft(dnr, delta, f, f1, f2)
    ff = wpps*f
    if (statavg) then
      cm(1) = wpps*dnr
      cm(2) = wpps*dns
      cm(3) = woss - 0.5_dp*ff
      cm(4) = 0.5_dp*ff - 0.5_dp*woss
    else
      cm(1) = dnr
      cm(2) = dns
      cm(3) = -0.5_dp*f
      cm(4) = -cm(3)
    end if
  end subroutine reks_coefficients

!###############################################################################
!  Main driver
!###############################################################################

  subroutine reks_driver(basis, infos, molGrid)
    use io_constants, only: iw
    use types, only: information
    use mod_dft_molgrid, only: dft_grid_t
    use oqp_tagarray_driver
    use messages, only: show_message, with_abort
    use mathlib, only: orb_to_dens, unpack_matrix, pack_matrix, &
                       orthogonal_transform_sym, traceprod_sym_packed
    use eigen, only: diag_symm_full
    use util, only: measure_time, e_charge_repulsion
    use printing, only: print_mo_range
    use mod_dft_gridint_energy, only: dmatd_density_blk
    use oqp_linalg

    implicit none

    character(len=*), parameter :: module_name = "reks"
    character(len=*), parameter :: subroutine_name = "reks_driver"

    type(basis_set), intent(in) :: basis
    type(information), target, intent(inout) :: infos
    type(dft_grid_t), intent(in) :: molGrid

    ! tag arrays
    real(kind=dp), contiguous, pointer :: hcore(:), smat(:), tmat(:)
    real(kind=dp), contiguous, pointer :: fock_a(:), fock_b(:), dmat_a(:), dmat_b(:)
    real(kind=dp), contiguous, pointer :: mo_a(:,:), mo_b(:,:)
    real(kind=dp), contiguous, pointer :: mo_energy_a(:), mo_energy_b(:)

    character(len=*), parameter :: tags_general(3) = &
      (/ character(len=80) :: OQP_Hcore, OQP_SM, OQP_TM /)
    character(len=*), parameter :: tags_alpha(4) = &
      (/ character(len=80) :: OQP_FOCK_A, OQP_DM_A, OQP_E_MO_A, OQP_VEC_MO_A /)
    character(len=*), parameter :: tags_beta(4) = &
      (/ character(len=80) :: OQP_FOCK_B, OQP_DM_B, OQP_E_MO_B, OQP_VEC_MO_B /)

    ! integrals
    type(int2_compute_t) :: int2_driver
    type(int2_reks_data_t) :: int2_data

    ! dimensions / control
    integer :: nbf, nbf_tri, nelec, na, ncore, ir, is, nvir0, maxit, iter, l
    integer :: reks_type, reks_target, nang, i, j, info, maxdiis
    logical :: is_dft, do_diis, converged, cvging_prev, cvging, dbg, cam
    real(kind=dp) :: hfscale, wpps, woss, delta, shift, conv
    real(kind=dp) :: dnr, dns, fr, fs, etot, eold, dele, offfock, diff, erdiis
    real(kind=dp) :: enuc, ereks, wrs, epps, eoss, edes, etrp, etarget
    real(kind=dp) :: exc(nmic), em(nmic), cm(nmic), cmpps(nmic)
    real(kind=dp) :: totele, totkin, ssr2(2,2), ssr3(3,3), w2(2), w3(3)
    real(kind=dp) :: sq2, ea, eb, s2, s3

    ! work arrays
    real(kind=dp), allocatable, target :: dblk(:,:)
    real(kind=dp), allocatable :: dblk_new(:,:), jmat(:,:), kmat(:,:)
    real(kind=dp), allocatable :: fa(:,:), fb(:,:), ga(:,:), gb(:,:)
    real(kind=dp), allocatable :: vxa(:,:), vxb(:,:), edfa(:), edfb(:)
    real(kind=dp), allocatable :: da(:,:), db(:,:), occ(:)
    real(kind=dp), allocatable :: famo(:,:,:), fbmo(:,:,:), ftmp(:)
    real(kind=dp), allocatable :: frex(:,:), wrex(:,:), umat(:,:), eig(:)
    real(kind=dp), allocatable :: smat_full(:,:), sc(:,:), work(:,:), errmo(:,:)
    real(kind=dp), allocatable :: diis_f(:,:), diis_e(:,:), fao(:,:), eao(:,:)
    integer :: ndiis
    real(kind=dp) :: wc, wra, wrb, wsa, wsb, wrca, wrcb, wsca, wscb, wsra, wsrb
    real(kind=dp) :: signrs, cl, ta, tb
    real(kind=dp), parameter :: lowlim = 1.0e-8_dp

    !--------------------------------------------------------------------------
    ! Parameters
    !--------------------------------------------------------------------------
    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2
    nelec = infos%mol_prop%nelec
    na = infos%mol_prop%nelec_a
    if (infos%mol_prop%nelec_b /= na) then
      call show_message('REKS(2,2) requires a singlet (nelec_a == nelec_b) reference', with_abort)
    end if
    if (na < 1 .or. na+1 > nbf) then
      call show_message('REKS(2,2): active orbitals r=nocc, s=nocc+1 are out of range', with_abort)
    end if
    ncore = na - 1
    ir = na
    is = na + 1
    nvir0 = na + 2
    maxit = infos%control%maxit
    conv = infos%control%conv
    maxdiis = max(2, int(infos%control%maxdiis))
    is_dft = infos%control%hamilton >= 20
    hfscale = 1.0_dp
    if (is_dft) hfscale = infos%dft%HFscale
    cam = is_dft .and. infos%dft%cam_flag

    reks_type = int(infos%control%reks_type)
    reks_target = int(infos%control%reks_target)
    wpps = infos%control%reks_wpps
    shift = infos%control%reks_shift
    delta = infos%control%reks_delta
    do_diis = infos%control%reks_diis /= 0
    dbg = infos%control%verbose >= 3
    if (abs(wpps - 1.0_dp) <= 1.0e-6_dp) then
      reks_type = 0
      reks_target = 0
    end if
    woss = 1.0_dp - wpps
    nang = maxval(basis%am) + 1 + 1
    sq2 = sqrt(2.0_dp)

    !--------------------------------------------------------------------------
    ! Tag arrays
    !--------------------------------------------------------------------------
    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_Hcore, hcore)
    call tagarray_get_data(infos%dat, OQP_SM, smat)
    call tagarray_get_data(infos%dat, OQP_TM, tmat)
    call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)

    !--------------------------------------------------------------------------
    ! Work arrays
    !--------------------------------------------------------------------------
    allocate(dblk(nbf_tri,4), dblk_new(nbf_tri,4), jmat(nbf_tri,4), kmat(nbf_tri,4), &
             fa(nbf_tri,nmic), fb(nbf_tri,nmic), ga(nbf_tri,nmic), gb(nbf_tri,nmic), &
             vxa(nbf_tri,nmic), vxb(nbf_tri,nmic), edfa(nmic), edfb(nmic), &
             da(nbf,nbf), db(nbf,nbf), occ(nbf), ftmp(nbf_tri), &
             famo(nbf,nbf,nmic), fbmo(nbf,nbf,nmic), frex(nbf,nbf), wrex(nbf,nbf), &
             umat(nbf,nbf), eig(nbf), smat_full(nbf,nbf), sc(nbf,nbf), work(nbf,nbf), &
             errmo(nbf,nbf), fao(nbf,nbf), eao(nbf,nbf), &
             diis_f(nbf*nbf, maxdiis), diis_e(nbf*nbf, maxdiis), source=0.0_dp)
    ndiis = 0

    call unpack_matrix(smat, smat_full, nbf, 'U')

    !--------------------------------------------------------------------------
    ! Banner
    !--------------------------------------------------------------------------
    write(iw,'(/3x,64("="))')
    select case (reks_type)
    case (0)
      if (abs(wpps-1.0_dp) <= 1.0e-6_dp) then
        write(iw,'(10x,a)') 'REKS(2,2) calculation'
      else
        write(iw,'(10x,a)') 'SA-REKS(2,2) calculation'
      end if
    case (1)
      write(iw,'(10x,a)') '2SI-2SA-REKS(2,2) (SSR(2,2)) calculation'
    case (2)
      write(iw,'(10x,a)') '3SI-2SA-REKS(2,2) (SSR(3,2)) calculation'
    end select
    write(iw,'(10x,a,i5,a,i5)')  'Active orbitals: r =', ir, ',  s =', is
    write(iw,'(10x,a,f10.6,a,f10.6)') 'Weights: W_PPS =', wpps, ',  W_OSS =', woss
    write(iw,'(10x,a,f8.4,a,f8.4)') 'Level shift =', shift, ',  delta =', delta
    if (do_diis) then
      write(iw,'(10x,a,i3)') 'DIIS on, max vectors =', maxdiis
    else
      write(iw,'(10x,a)') 'DIIS off'
    end if
    write(iw,'(10x,a,i5,a,es10.2)') 'maxit =', maxit, ',  conv =', conv
    write(iw,'(3x,64("="))')

    call measure_time(print_total=1, log_unit=iw)

    enuc = e_charge_repulsion(infos%atoms%xyz, infos%atoms%zn - infos%basis%ecp_zn_num)

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    call flush(iw)

    int2_data%d => dblk
    int2_data%scale_exchange = hfscale
    int2_data%scale_coulomb = 1.0_dp

    ! initial FONs: closed shell
    dnr = 1.0_dp
    dns = 0.0_dp
    call reks_coefficients(dnr, dns, delta, wpps, woss, cm)
    fr = dnr*wpps + 0.5_dp*woss
    fs = dns*wpps + 0.5_dp*woss

    call build_blocks(mo_a, dblk)

    eold = 0.0_dp
    converged = .false.
    cvging_prev = .false.
    erdiis = 0.0_dp
    wrs = 0.0_dp
    etot = 0.0_dp

    write(iw,'(/3x,a)') ' ITER   NR/2         TOTAL ENERGY     OFFDIAG FOCK      E CHANGE   DENSITY CHANGE    DIIS ERROR'
    write(iw,'(3x,110("-"))')

    !==========================================================================
    ! SCF iterations
    !==========================================================================
    do iter = 1, maxit

      !--- two-electron J/K for the four density blocks
      call int2_driver%run(int2_data, cam=cam, alpha=infos%dft%cam_alpha, &
                           beta=infos%dft%cam_beta, mu=infos%dft%cam_mu)
      jmat = 0.5_dp*int2_data%f(:,1:4,1)
      kmat = 0.5_dp*int2_data%f(:,5:8,1)
      do i = 1, 4
        j = 0
        do l = 1, nbf
          j = j + l
          jmat(j,i) = 2.0_dp*jmat(j,i)
          kmat(j,i) = 2.0_dp*kmat(j,i)
        end do
      end do

      !--- microstate two-electron parts and XC
      exc = 0.0_dp
      vxa = 0.0_dp
      vxb = 0.0_dp
      do l = 1, nmic
        ga(:,l) = jmat(:,mic_a(l)) + jmat(:,mic_b(l)) - kmat(:,mic_a(l))
        gb(:,l) = jmat(:,mic_a(l)) + jmat(:,mic_b(l)) - kmat(:,mic_b(l))
        if (is_dft) then
          call unpack_matrix(dblk(:,mic_a(l)), da, nbf, 'U')
          call unpack_matrix(dblk(:,mic_b(l)), db, nbf, 'U')
          call dmatd_density_blk(basis, molGrid, da, db, vxa(:,l), vxb(:,l), &
                                 exc(l), totele, totkin, nang, nbf, &
                                 infos%dft%grid_density_cutoff, .true., infos)
          if (ncore == 0 .and. l == 4) vxb(:,l) = 0.0_dp
        end if
        fa(:,l) = hcore + ga(:,l) + vxa(:,l)
        fb(:,l) = hcore + gb(:,l) + vxb(:,l)
      end do

      !--- microstate energies
      do l = 1, nmic
        ea = traceprod_sym_packed(dblk(:,mic_a(l)), hcore, nbf) &
           + 0.5_dp*traceprod_sym_packed(dblk(:,mic_a(l)), ga(:,l), nbf)
        eb = traceprod_sym_packed(dblk(:,mic_b(l)), hcore, nbf) &
           + 0.5_dp*traceprod_sym_packed(dblk(:,mic_b(l)), gb(:,l), nbf)
        em(l) = ea + eb + exc(l) + enuc
      end do

      !--- fractional occupation numbers
      call reks_fon_solve(em, delta, dnr, ereks, iw, dbg)
      if (dnr < 0.5_dp .and. .not. do_diis) then
        ! swap the active orbitals so that r is the more occupied one
        dnr = 1.0_dp - dnr
        occ(1:nbf) = mo_a(:,ir); mo_a(:,ir) = mo_a(:,is); mo_a(:,is) = occ(1:nbf)
        ea = em(1); em(1) = em(2); em(2) = ea
        ftmp = fa(:,1); fa(:,1) = fa(:,2); fa(:,2) = ftmp
        ftmp = fb(:,1); fb(:,1) = fb(:,2); fb(:,2) = ftmp
        ftmp = fa(:,3); fa(:,3) = fb(:,3); fb(:,3) = ftmp
        ftmp = dblk(:,2); dblk(:,2) = dblk(:,3); dblk(:,3) = ftmp
      end if
      dns = 1.0_dp - dnr
      call reks_coefficients(dnr, dns, delta, wpps, woss, cm)
      fr = dnr*wpps + 0.5_dp*woss
      fs = dns*wpps + 0.5_dp*woss

      !--- microstate Fock matrices to the MO basis
      do l = 1, nmic
        call orthogonal_transform_sym(nbf, nbf, fa(:,l), mo_a, nbf, ftmp)
        call unpack_matrix(ftmp, famo(:,:,l), nbf, 'U')
        call orthogonal_transform_sym(nbf, nbf, fb(:,l), mo_a, nbf, ftmp)
        call unpack_matrix(ftmp, fbmo(:,:,l), nbf, 'U')
      end do

      !--- effective (coupling-operator) Fock matrix, GAMESS REXFM2FE
      frex = 0.0_dp
      wrex = 0.0_dp
      wrs = 0.0_dp
      signrs = 1.0_dp
      if (fr - fs < 0.0_dp) signrs = -1.0_dp
      do l = 1, nmic
        cl = cm(l)*fact(l)
        wc = 0.5_dp*cl
        wra = wc*nra(l); if (fr > lowlim) wra = wra/fr
        wrb = wc*nrb(l); if (fr > lowlim) wrb = wrb/fr
        wsa = wc*nsa(l); if (fs > lowlim) wsa = wsa/fs
        wsb = wc*nsb(l); if (fs > lowlim) wsb = wsb/fs
        wrca = wc*(1-nra(l)); if (1.0_dp-fr > lowlim) wrca = wrca/(1.0_dp-fr)
        wrcb = wc*(1-nrb(l)); if (1.0_dp-fr > lowlim) wrcb = wrcb/(1.0_dp-fr)
        wsca = wc*(1-nsa(l)); if (1.0_dp-fs > lowlim) wsca = wsca/(1.0_dp-fs)
        wscb = wc*(1-nsb(l)); if (1.0_dp-fs > lowlim) wscb = wscb/(1.0_dp-fs)
        wsra = cl*(nra(l)-nsa(l))*signrs
        wsrb = cl*(nrb(l)-nsb(l))*signrs

        ! core-core
        do i = 1, ncore
          do j = 1, i
            frex(i,j) = frex(i,j) + wc*(famo(i,j,l) + fbmo(i,j,l))
            wrex(i,j) = wrex(i,j) + cl*(famo(i,j,l) + fbmo(i,j,l))
          end do
        end do
        ! virt-virt
        do i = nvir0, nbf
          do j = nvir0, i
            frex(i,j) = frex(i,j) + wc*(famo(i,j,l) + fbmo(i,j,l))
          end do
        end do
        ! virt-core
        do i = nvir0, nbf
          do j = 1, ncore
            frex(i,j) = frex(i,j) + wc*(famo(i,j,l) + fbmo(i,j,l))
          end do
        end do
        ! r-r
        frex(ir,ir) = frex(ir,ir) + wra*famo(ir,ir,l) + wrb*fbmo(ir,ir,l)
        wrex(ir,ir) = wrex(ir,ir) + cl*(nra(l)*famo(ir,ir,l) + nrb(l)*fbmo(ir,ir,l))
        ! s-s
        frex(is,is) = frex(is,is) + wsa*famo(is,is,l) + wsb*fbmo(is,is,l)
        wrex(is,is) = wrex(is,is) + cl*(nsa(l)*famo(is,is,l) + nsb(l)*fbmo(is,is,l))
        ! r-core, s-core
        do j = 1, ncore
          frex(ir,j) = frex(ir,j) + wrca*famo(ir,j,l) + wrcb*fbmo(ir,j,l)
          wrex(ir,j) = wrex(ir,j) + cl*(nra(l)*famo(ir,j,l) + nrb(l)*fbmo(ir,j,l))
          frex(is,j) = frex(is,j) + wsca*famo(is,j,l) + wscb*fbmo(is,j,l)
          wrex(is,j) = wrex(is,j) + cl*(nsa(l)*famo(is,j,l) + nsb(l)*fbmo(is,j,l))
        end do
        ! virt-r, virt-s
        do i = nvir0, nbf
          frex(i,ir) = frex(i,ir) + wra*famo(i,ir,l) + wrb*fbmo(i,ir,l)
          frex(i,is) = frex(i,is) + wsa*famo(i,is,l) + wsb*fbmo(i,is,l)
        end do
        ! s-r
        frex(is,ir) = frex(is,ir) + wsra*famo(is,ir,l) + wsrb*fbmo(is,ir,l)
        wrs = wrs + cl*(nra(l)*famo(is,ir,l) + nrb(l)*fbmo(is,ir,l))
      end do
      wrex(is,ir) = wrs
      ! symmetrize (lower triangle was filled)
      do i = 1, nbf
        do j = 1, i-1
          frex(j,i) = frex(i,j)
          wrex(j,i) = wrex(i,j)
        end do
      end do

      !--- largest off-diagonal element of the effective Fock (GAMESS REXOFFFOC)
      offfock = 0.0_dp
      do j = 1, ncore
        offfock = max(offfock, abs(frex(ir,j)), abs(frex(is,j)))
      end do
      offfock = max(offfock, abs(frex(is,ir)))
      do i = nvir0, nbf
        do j = 1, is
          offfock = max(offfock, abs(frex(i,j)))
        end do
      end do

      !--- total (SA) energy
      etot = 0.0_dp
      do l = 1, nmic
        etot = etot + fact(l)*cm(l)*em(l)
      end do
      dele = etot - eold
      eold = etot

      !--- DIIS on the AO representation of the effective Fock
      erdiis = 0.0_dp
      if (do_diis) then
        errmo = 0.0_dp
        do j = 1, ncore
          errmo(ir,j) = frex(ir,j)*(1.0_dp-fr); errmo(j,ir) = errmo(ir,j)
          errmo(is,j) = frex(is,j)*(1.0_dp-fs); errmo(j,is) = errmo(is,j)
        end do
        do i = nvir0, nbf
          errmo(i,ir) = frex(i,ir)*fr; errmo(ir,i) = errmo(i,ir)
          errmo(i,is) = frex(i,is)*fs; errmo(is,i) = errmo(i,is)
        end do
        errmo(is,ir) = frex(is,ir)*(fr-fs); errmo(ir,is) = errmo(is,ir)
        erdiis = maxval(abs(errmo))
        ! S C
        call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, smat_full, nbf, mo_a, nbf, 0.0_dp, sc, nbf)
        call mo_to_ao_full(errmo, sc, eao, work, nbf)
        call mo_to_ao_full(frex, sc, fao, work, nbf)
        call diis_push(fao, eao)
        call diis_extrapolate(fao)
        ! back to the MO basis: C^T F_AO C
        call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mo_a, nbf, fao, nbf, 0.0_dp, work, nbf)
        call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, work, nbf, mo_a, nbf, 0.0_dp, frex, nbf)
      end if

      !--- level shift: r + shift, s + 2 shift, virtuals + 3 shift
      frex(ir,ir) = frex(ir,ir) + shift
      frex(is,is) = frex(is,is) + 2.0_dp*shift
      do i = nvir0, nbf
        frex(i,i) = frex(i,i) + 3.0_dp*shift
      end do

      !--- diagonalize in the MO basis and rotate the orbitals
      umat = frex
      call diag_symm_full(1, nbf, umat, nbf, eig, info)
      if (info /= 0) call show_message('REKS: diagonalization of the effective Fock failed', with_abort)
      do i = 1, nbf
        if (umat(i,i) < 0.0_dp) umat(:,i) = -umat(:,i)
      end do
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, mo_a, nbf, umat, nbf, 0.0_dp, work, nbf)
      mo_a = work
      eig(ir) = eig(ir) - shift
      eig(is) = eig(is) - 2.0_dp*shift
      do i = nvir0, nbf
        eig(i) = eig(i) - 3.0_dp*shift
      end do
      mo_energy_a = eig

      !--- new density blocks and convergence
      call build_blocks(mo_a, dblk_new)
      diff = maxval(abs(dblk_new - dblk))
      dblk = dblk_new

      write(iw,'(3x,i4,2x,f11.8,2x,f20.10,2x,f14.9,2x,es14.6,2x,es14.6,2x,es14.6)') &
        iter, dnr, etot, offfock, dele, diff, erdiis
      call flush(iw)

      cvging = (diff < conv) .or. (offfock < conv)
      if (do_diis) cvging = cvging .or. (erdiis < conv .and. diff < 2.0_dp*conv)
      if (cvging .and. cvging_prev) then
        converged = .true.
        exit
      end if
      cvging_prev = cvging
    end do

    call int2_driver%clean()

    !==========================================================================
    ! Final report
    !==========================================================================
    write(iw,'(3x,110("-"))')
    if (converged) then
      write(iw,"(10x,'REKS SCF convergence achieved ....')")
      infos%mol_energy%SCF_converged = .true.
    else
      write(iw,"(10x,'REKS SCF did not converge.')")
      infos%mol_energy%SCF_converged = .false.
      iter = maxit
    end if
    write(iw,"(/' Final SA-REKS(2,2) energy is',F20.10,' after',I4,' iterations')") etot, iter
    write(iw,"(' SA-REKS(2,2): FON(',I4,') =',F10.6,'  FON(',I4,') =',F10.6/)") ir, 2.0_dp*dnr, is, 2.0_dp*dns

    ! individual states with the pure-PPS coefficients
    call reks_coefficients(dnr, dns, delta, 1.0_dp, 0.0_dp, cmpps)
    etrp = em(4)
    epps = cmpps(1)*em(1) + cmpps(2)*em(2) + 2.0_dp*cmpps(3)*em(3) + 2.0_dp*cmpps(4)*em(4)
    edes = cmpps(2)*em(1) + cmpps(1)*em(2) - 2.0_dp*cmpps(3)*em(3) - 2.0_dp*cmpps(4)*em(4)
    eoss = 2.0_dp*em(3) - em(4)

    write(iw,'(3x,a)') 'REKS(2,2) individual state energies (SA orbitals and FONs):'
    write(iw,'(3x,a,f22.12)') 'Triplet:                        ', etrp
    write(iw,'(3x,a,f22.12)') 'Doubly excited singlet (DES):   ', edes
    write(iw,'(3x,a,f22.12)') 'Open-shell singlet (OSS):       ', eoss
    write(iw,'(3x,a,f22.12)') 'Perfectly paired singlet (PPS): ', epps
    write(iw,'(3x,a,f22.12/)') 'Lagrangian W_rs:                ', wrs

    etarget = etot
    select case (reks_type)
    case (0)
      if (reks_target == 1) etarget = epps
      if (reks_target == 2) etarget = eoss
    case (1)
      ssr2 = 0.0_dp
      ssr2(1,1) = epps
      ssr2(2,2) = eoss
      ssr2(1,2) = wrs*(sqrt(dnr) - sqrt(dns))*sq2
      ssr2(2,1) = ssr2(1,2)
      call diag_symm_full(1, 2, ssr2, 2, w2, info)
      write(iw,'(3x,a)') '2SI-2SA-REKS(2,2) (SSR) states:'
      write(iw,'(3x,a)') '                   E_k              C_PPS        C_OSS'
      do i = 1, 2
        write(iw,'(3x,a,i2,f22.12,2f13.8)') 'SSR state', i-1, w2(i), ssr2(1,i), ssr2(2,i)
      end do
      write(iw,*)
      if (reks_target == 1) etarget = w2(1)
      if (reks_target == 2) etarget = w2(2)
    case (2)
      ssr3 = 0.0_dp
      ssr3(1,1) = epps
      ssr3(2,2) = eoss
      ssr3(3,3) = edes
      ssr3(1,2) = wrs*(sqrt(dnr) - sqrt(dns))*sq2
      ssr3(2,1) = ssr3(1,2)
      ssr3(2,3) = wrs*(sqrt(dnr) + sqrt(dns))*sq2
      ssr3(3,2) = ssr3(2,3)
      call diag_symm_full(1, 3, ssr3, 3, w3, info)
      write(iw,'(3x,a)') '3SI-2SA-REKS(2,2) (SSR(3,2)) states:'
      write(iw,'(3x,a)') '                   E_k              C_PPS        C_OSS        C_DES'
      do i = 1, 3
        write(iw,'(3x,a,i2,f22.12,3f13.8)') 'SSR state', i-1, w3(i), ssr3(1,i), ssr3(2,i), ssr3(3,i)
      end do
      write(iw,*)
      if (reks_target == 1) etarget = w3(1)
      if (reks_target == 2) etarget = w3(2)
    end select

    select case (reks_target)
    case (0)
      write(iw,'(3x,a,f22.12/)') 'REKS: reporting the SA-averaged energy:', etarget
    case (1)
      write(iw,'(3x,a,f22.12/)') 'REKS: reporting the S0 (PPS) state energy:', etarget
    case (2)
      write(iw,'(3x,a,f22.12/)') 'REKS: reporting the S1 (OSS) state energy:', etarget
    end select

    !--------------------------------------------------------------------------
    ! Store results: SA-averaged spin densities and Fock matrices, orbitals
    !--------------------------------------------------------------------------
    dmat_a = 0.0_dp
    dmat_b = 0.0_dp
    fock_a = 0.0_dp
    fock_b = 0.0_dp
    do l = 1, nmic
      dmat_a = dmat_a + fact(l)*cm(l)*dblk(:,mic_a(l))
      dmat_b = dmat_b + fact(l)*cm(l)*dblk(:,mic_b(l))
      fock_a = fock_a + fact(l)*cm(l)*fa(:,l)
      fock_b = fock_b + fact(l)*cm(l)*fb(:,l)
    end do
    mo_b = mo_a
    mo_energy_b = mo_energy_a

    call print_mo_range(basis, infos, mostart=1, moend=nbf)

    infos%mol_energy%energy = etarget
    infos%mol_energy%etot = etarget
    infos%mol_energy%nenergy = enuc
    infos%mol_energy%enuc = enuc
    infos%mol_energy%vnn = enuc
    infos%mol_energy%psinrm = (traceprod_sym_packed(dmat_a, smat, nbf) &
                             + traceprod_sym_packed(dmat_b, smat, nbf))/nelec
    infos%mol_energy%tkin = traceprod_sym_packed(dmat_a, tmat, nbf) &
                          + traceprod_sym_packed(dmat_b, tmat, nbf)
    infos%mol_energy%ehf1 = traceprod_sym_packed(dmat_a, hcore, nbf) &
                          + traceprod_sym_packed(dmat_b, hcore, nbf)
    infos%mol_energy%vne = infos%mol_energy%ehf1 - infos%mol_energy%tkin
    infos%mol_energy%vee = etot - enuc - infos%mol_energy%ehf1
    infos%mol_energy%vtot = infos%mol_energy%vne + infos%mol_energy%vee + enuc
    if (infos%mol_energy%tkin /= 0.0_dp) &
      infos%mol_energy%virial = -infos%mol_energy%vtot/infos%mol_energy%tkin

    write(iw,'(/3x,a,f20.10)') 'ONE ELECTRON ENERGY     =', infos%mol_energy%ehf1
    write(iw,'(3x,a,f20.10)')  'TWO ELECTRON ENERGY     =', infos%mol_energy%vee
    write(iw,'(3x,a,f20.10)')  'NUCLEAR REPULSION ENERGY=', enuc
    write(iw,'(3x,a,f20.10/)') 'TOTAL ENERGY            =', etarget

    call measure_time(print_total=1, log_unit=iw)

  contains

    !> Density blocks: 1=core, 2=core+r, 3=core+s, 4=core+r+s (per spin, occ 1)
    subroutine build_blocks(c, d)
      real(kind=dp), intent(in) :: c(:,:)
      real(kind=dp), intent(out) :: d(:,:)
      occ = 0.0_dp
      occ(1:ncore) = 1.0_dp
      if (ncore > 0) then
        call orb_to_dens(d(:,1), c, occ, ncore, nbf, nbf)
      else
        d(:,1) = 0.0_dp
      end if
      occ(ir) = 1.0_dp; occ(is) = 0.0_dp
      call orb_to_dens(d(:,2), c, occ, is, nbf, nbf)
      occ(ir) = 0.0_dp; occ(is) = 1.0_dp
      call orb_to_dens(d(:,3), c, occ, is, nbf, nbf)
      occ(ir) = 1.0_dp; occ(is) = 1.0_dp
      call orb_to_dens(d(:,4), c, occ, is, nbf, nbf)
    end subroutine build_blocks

    !> A_AO = (S C) A_MO (S C)^T
    subroutine mo_to_ao_full(amo, scmat, aao, wrk, n)
      integer, intent(in) :: n
      real(kind=dp), intent(in) :: amo(n,n), scmat(n,n)
      real(kind=dp), intent(out) :: aao(n,n), wrk(n,n)
      call dgemm('n','n', n, n, n, 1.0_dp, scmat, n, amo, n, 0.0_dp, wrk, n)
      call dgemm('n','t', n, n, n, 1.0_dp, wrk, n, scmat, n, 0.0_dp, aao, n)
    end subroutine mo_to_ao_full

    subroutine diis_push(f, e)
      real(kind=dp), intent(in) :: f(:,:), e(:,:)
      integer :: k
      if (ndiis < maxdiis) then
        ndiis = ndiis + 1
      else
        do k = 1, maxdiis-1
          diis_f(:,k) = diis_f(:,k+1)
          diis_e(:,k) = diis_e(:,k+1)
        end do
      end if
      diis_f(:,ndiis) = reshape(f, [nbf*nbf])
      diis_e(:,ndiis) = reshape(e, [nbf*nbf])
    end subroutine diis_push

    subroutine diis_extrapolate(f)
      real(kind=dp), intent(out) :: f(:,:)
      real(kind=dp), allocatable :: bmat(:,:), rhs(:)
      integer, allocatable :: ipiv(:)
      integer :: k, m, n1, ierr
      if (ndiis < 2) then
        f = reshape(diis_f(:,ndiis), [nbf,nbf])
        return
      end if
      n1 = ndiis + 1
      allocate(bmat(n1,n1), rhs(n1), ipiv(n1))
      bmat = 0.0_dp
      do k = 1, ndiis
        do m = 1, k
          bmat(k,m) = dot_product(diis_e(:,k), diis_e(:,m))
          bmat(m,k) = bmat(k,m)
        end do
        bmat(k,n1) = -1.0_dp
        bmat(n1,k) = -1.0_dp
      end do
      rhs = 0.0_dp
      rhs(n1) = -1.0_dp
      call dgesv(n1, 1, bmat, n1, ipiv, rhs, n1, ierr)
      if (ierr /= 0) then
        ! singular system: drop the oldest vector and use the latest Fock
        f = reshape(diis_f(:,ndiis), [nbf,nbf])
        ndiis = 0
        return
      end if
      f = 0.0_dp
      do k = 1, ndiis
        f = f + rhs(k)*reshape(diis_f(:,k), [nbf,nbf])
      end do
    end subroutine diis_extrapolate

  end subroutine reks_driver

end module reks
