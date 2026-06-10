!> @brief Real-solid-harmonic multipole construction in the EXACT ddX
!>        convention, for building a full-density ddPCM adjoint source (Psi).
!>
!> @details The ddX library represents a (classical) solute as a set of
!> atom-centred real-spherical-harmonic multipoles M_l^m and maps them to the
!> adjoint source via
!>
!>     psi(lm,isph) = 4*pi/((2l+1) * rsph(isph)^l) * M_lm(isph)            (ddX)
!>
!> (see ddx_multipolar_solutes.f90::multipole_psi). For a QM solute the
!> physically correct full-density source therefore requires the SAME multipole
!> convention ddX uses internally; otherwise the per-l normalisation/ordering
!> silently disagrees and Psi is wrong.
!>
!> To guarantee an exact match this module replicates, verbatim (only constants
!> renamed to OpenQP's precision module), the ddX routines that build multipoles
!> from point charges:
!>   * ylmscale     -> ddx_harmonics.f90::ylmscale  (scaling factors vscales)
!>   * polleg_work  -> ddx_harmonics.f90::polleg_work
!>   * trgev        -> ddx_harmonics.f90::trgev
!>   * ylmbas       -> ddx_harmonics.f90::ylmbas
!>   * accumulate_point_multipole -> ddx_harmonics.f90::fmm_p2m_work with the
!>     output radius dst_r = 1 and beta = 1, i.e. it accumulates the BARE moment
!>         M_lm += q * ||c||^l * Y_lm(c/||c||)
!>     which is exactly the multipole array consumed by ddx_multipole_psi /
!>     ddx_multipole_electrostatics (the C adapter sets the l=0 monopole to
!>     charge/sqrt(4*pi), i.e. q*Y_00, which this routine reproduces).
!>
!> Because the harmonic evaluation is ddX's own code, the convention (ordering
!> i = l*l+l+1+m, m = -l..l; sign of the sine/cosine harmonics; orthonormal real
!> Y_lm) is identical by construction for all l, not just the hand-checked
!> l<=2 table in solvent_pcm.F90::pack_ddx_l2_multipoles.
module solvent_pcm_harmonics

  use precision, only: dp

  implicit none
  private

  public :: pcm_ylmscale
  public :: accumulate_point_multipole
  public :: accumulate_point_multipole_leak

  real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp
  real(dp), parameter :: three = 3.0_dp, pt5 = 0.5_dp
  real(dp), parameter :: sqrt2 = sqrt(two)
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: fourpi = 4.0_dp * pi
  real(dp), parameter :: sqrt4pi = sqrt(fourpi)

contains

  !> @brief Scaling factors of real normalised spherical harmonics.
  !> Verbatim replica of ddx_harmonics.f90::ylmscale (vscales output only).
  subroutine pcm_ylmscale(p, vscales)
    integer,  intent(in)  :: p
    real(dp), intent(out) :: vscales((p+1)**2)

    real(dp) :: tmp, twolp1
    integer  :: l, ind, m

    twolp1 = one
    do l = 0, p
      ind = l*l + l + 1
      tmp = fourpi / twolp1
      tmp = sqrt(tmp)
      vscales(ind) = one / tmp
      twolp1 = twolp1 + two
      tmp = vscales(ind) * sqrt2
      do m = 1, l
        tmp = -tmp / sqrt(dble((l-m+1)*(l+m)))
        vscales(ind+m) = tmp
        vscales(ind-m) = tmp
      end do
    end do
  end subroutine pcm_ylmscale

  !> @brief Accumulate the bare ddX-convention multipole of a point charge.
  !>
  !> dst_m(lm) += src_q * ||c||^l * Y_lm(c/||c||)   for l = 0..p.
  !>
  !> Equivalent to ddx_harmonics.f90::fmm_p2m_work with dst_r = 1, beta = 1.
  !> @param[in]    c        radius vector from the charge to the harmonic centre
  !> @param[in]    src_q    charge of the source point
  !> @param[in]    p        maximal degree
  !> @param[in]    vscales  scaling factors from pcm_ylmscale, dim (p+1)**2
  !> @param[inout] dst_m    multipole coefficients, dim (p+1)**2 (accumulated)
  subroutine accumulate_point_multipole(c, src_q, p, vscales, dst_m)
    integer,  intent(in)    :: p
    real(dp), intent(in)    :: c(3), src_q, vscales((p+1)**2)
    real(dp), intent(inout) :: dst_m((p+1)**2)

    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), vcos(p+1), vsin(p+1)
    real(dp) :: rho, ctheta, stheta, cphi, sphi, t
    integer  :: n, ind

    if (src_q == zero) return

    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
                vcos, vsin)

    if (rho /= zero) then
      t = src_q
      do n = 0, p
        ind = n*n + n + 1
        dst_m(ind-n:ind+n) = dst_m(ind-n:ind+n) + t*vylm(ind-n:ind+n)
        t = t * rho
      end do
    else
      dst_m(1) = dst_m(1) + src_q / sqrt4pi
    end if
  end subroutine accumulate_point_multipole

  !> @brief Accumulate the ddX-convention multipole of a point charge WITH the
  !>        ddCOSMO/ddPCM outside-sphere "leak" correction (PySCF make_psi_vmat /
  !>        cache_fake_multipoles equivalent).
  !>
  !> ddX's multipole interface (multipole_psi: psi = 4*pi/((2l+1) rsph^l) M_lm)
  !> assumes every source charge lies INSIDE its sphere, so it consumes the bare
  !> interior moment M_lm += q*r^l*Y_lm. For a QM density the tail extends beyond
  !> the (small) vdW sphere; for those points the bare interior moment q*r^l (with
  !> r>rsph) is the WRONG continuation and blows up with l. PySCF handles this by
  !> scaling the fake-multipole factor by (r_vdw/r)^(2l+1) for points outside the
  !> sphere, i.e. it switches to the EXTERIOR multipole. Expressed in the bare ddX
  !> moment that multipole_psi consumes, the same correction is:
  !>
  !>   r <= rsph :  M_lm += q * r^l * Y_lm                         (interior)
  !>   r >  rsph :  M_lm += q * Y_lm * rsph^(2l+1) / r^(l+1)       (exterior leak)
  !>                       = (interior term) * (rsph/r)^(2l+1)
  !>
  !> One verifies psi_lm = 4*pi/((2l+1) rsph^l) * M_lm then reproduces PySCF's
  !> leak-corrected psi up to the fixed per-sphere rsph factor that distinguishes
  !> ddX's psi normalisation (rsph^l) from PySCF's (rsph^(l+1)); that factor is a
  !> convention difference between the two libraries, not a physical one.
  !>
  !> @param[in]    c        radius vector from the charge to the sphere centre
  !> @param[in]    src_q    charge of the source point
  !> @param[in]    p        maximal degree
  !> @param[in]    vscales  scaling factors from pcm_ylmscale, dim (p+1)**2
  !> @param[in]    rsph     radius of the sphere the moment is accumulated on
  !> @param[inout] dst_m    multipole coefficients, dim (p+1)**2 (accumulated)
  subroutine accumulate_point_multipole_leak(c, src_q, p, vscales, rsph, dst_m)
    integer,  intent(in)    :: p
    real(dp), intent(in)    :: c(3), src_q, vscales((p+1)**2), rsph
    real(dp), intent(inout) :: dst_m((p+1)**2)

    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), vcos(p+1), vsin(p+1)
    real(dp) :: rho, ctheta, stheta, cphi, sphi, t, ratio2
    integer  :: n, ind

    if (src_q == zero) return

    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
                vcos, vsin)

    if (rho == zero) then
      ! Exactly at the sphere centre (e.g. a nuclear charge): always interior,
      ! only the l=0 harmonic survives.
      dst_m(1) = dst_m(1) + src_q / sqrt4pi
      return
    end if

    if (rho <= rsph .or. rsph <= zero) then
      ! Interior point: standard bare moment q*r^l*Y_lm.
      t = src_q
      do n = 0, p
        ind = n*n + n + 1
        dst_m(ind-n:ind+n) = dst_m(ind-n:ind+n) + t*vylm(ind-n:ind+n)
        t = t * rho
      end do
    else
      ! Exterior (leak) point: q * rsph^(2l+1)/rho^(l+1) * Y_lm
      ! (= the interior term q*rho^l scaled by (rsph/rho)^(2l+1)).
      !   t_0 = src_q * rsph/rho ;  t_l = t_{l-1} * rsph^2/rho
      ratio2 = rsph * (rsph/rho)
      t = src_q * (rsph/rho)
      do n = 0, p
        ind = n*n + n + 1
        dst_m(ind-n:ind+n) = dst_m(ind-n:ind+n) + t*vylm(ind-n:ind+n)
        t = t * ratio2
      end do
    end if
  end subroutine accumulate_point_multipole_leak

  !> @brief Real spherical harmonics at a point. Verbatim replica of
  !>        ddx_harmonics.f90::ylmbas.
  subroutine ylmbas(x, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, &
                    vplm, vcos, vsin)
    integer,  intent(in)  :: p
    real(dp), intent(in)  :: x(3)
    real(dp), intent(in)  :: vscales((p+1)**2)
    real(dp), intent(out) :: rho, ctheta, stheta, cphi, sphi
    real(dp), intent(out) :: vylm((p+1)**2), vplm((p+1)**2)
    real(dp), intent(out) :: vcos(p+1), vsin(p+1)

    integer  :: l, m, ind
    real(dp) :: max12, ssq12, tmp

    if (x(1) == zero) then
      max12 = abs(x(2))
      ssq12 = one
    else if (abs(x(2)) > abs(x(1))) then
      max12 = abs(x(2))
      ssq12 = one + (x(1)/x(2))**2
    else
      max12 = abs(x(1))
      ssq12 = one + (x(2)/x(1))**2
    end if
    if (x(3) == zero) then
      rho = max12 * sqrt(ssq12)
    else if (abs(x(3)) > max12) then
      rho = one + ssq12 * (max12/x(3))**2
      rho = abs(x(3)) * sqrt(rho)
    else
      rho = ssq12 + (x(3)/max12)**2
      rho = max12 * sqrt(rho)
    end if
    if (rho == zero) return

    stheta = max12 * sqrt(ssq12)
    if (stheta /= zero) then
      cphi = x(1) / stheta
      sphi = x(2) / stheta
      ctheta = x(3) / rho
      stheta = stheta / rho
      select case (p)
      case (0)
        vylm(1) = vscales(1)
        return
      case (1)
        vylm(1) = vscales(1)
        vylm(2) = -vscales(4) * stheta * sphi
        vylm(3) = vscales(3) * ctheta
        vylm(4) = -vscales(4) * stheta * cphi
        return
      end select
      call polleg_work(ctheta, stheta, p, vplm, vcos)
      call trgev(cphi, sphi, p, vcos, vsin)
      vylm(1) = vscales(1)
      vylm(2) = -vscales(4) * stheta * sphi
      vylm(3) = vscales(3) * ctheta
      vylm(4) = -vscales(4) * stheta * cphi
      ind = 3
      do l = 2, p
        ind = ind + 2*l
        vylm(ind) = vscales(ind) * vplm(ind)
        do m = 1, l
          tmp = vplm(ind+m) * vscales(ind+m)
          vylm(ind+m) = tmp * vcos(m+1)
          vylm(ind-m) = tmp * vsin(m+1)
        end do
      end do
    else
      cphi = one
      sphi = zero
      ctheta = sign(one, x(3))
      stheta = zero
      vcos = one
      vsin = zero
      vylm = zero
      vplm = zero
      do l = 0, p, 2
        ind = l**2 + l + 1
        vplm(ind) = one
        vylm(ind) = vscales(ind)
      end do
      do l = 1, p, 2
        ind = l**2 + l + 1
        vplm(ind) = ctheta
        vylm(ind) = ctheta * vscales(ind)
      end do
    end if
  end subroutine ylmbas

  !> @brief cos(m*phi)/sin(m*phi) arrays. Verbatim replica of
  !>        ddx_harmonics.f90::trgev.
  subroutine trgev(cphi, sphi, p, vcos, vsin)
    integer,  intent(in)  :: p
    real(dp), intent(in)  :: cphi, sphi
    real(dp), intent(out) :: vcos(p+1), vsin(p+1)

    integer  :: m
    real(dp) :: c4phi, s4phi

    select case (p)
    case (:-1)
      return
    case (0)
      vcos(1) = one
      vsin(1) = zero
      return
    case (1)
      vcos(1) = one
      vsin(1) = zero
      vcos(2) = cphi
      vsin(2) = sphi
      return
    case (2)
      vcos(1) = one
      vsin(1) = zero
      vcos(2) = cphi
      vsin(2) = sphi
      vcos(3) = cphi**2 - sphi**2
      vsin(3) = 2 * cphi * sphi
      return
    case (3)
      vcos(1) = one
      vsin(1) = zero
      vcos(2) = cphi
      vsin(2) = sphi
      vcos(3) = cphi**2 - sphi**2
      vsin(3) = 2 * cphi * sphi
      vcos(4) = vcos(3)*cphi - vsin(3)*sphi
      vsin(4) = vcos(3)*sphi + vsin(3)*cphi
      return
    case default
      vcos(1) = one
      vsin(1) = zero
      vcos(2) = cphi
      vsin(2) = sphi
      vcos(3) = cphi**2 - sphi**2
      vsin(3) = 2 * cphi * sphi
      vcos(4) = vcos(3)*cphi - vsin(3)*sphi
      vsin(4) = vcos(3)*sphi + vsin(3)*cphi
      vcos(5) = vcos(3)**2 - vsin(3)**2
      vsin(5) = 2 * vcos(3) * vsin(3)
      c4phi = vcos(5)
      s4phi = vsin(5)
    end select
    do m = 6, p-2, 4
      vcos(m:m+3) = vcos(m-4:m-1)*c4phi - vsin(m-4:m-1)*s4phi
      vsin(m:m+3) = vcos(m-4:m-1)*s4phi + vsin(m-4:m-1)*c4phi
    end do
    select case (m-p)
    case (-1)
      vcos(p-1) = vcos(p-2)*cphi - vsin(p-2)*sphi
      vsin(p-1) = vcos(p-2)*sphi + vsin(p-2)*cphi
      vcos(p) = vcos(p-2)*vcos(3) - vsin(p-2)*vsin(3)
      vsin(p) = vcos(p-2)*vsin(3) + vsin(p-2)*vcos(3)
      vcos(p+1) = vcos(p-2)*vcos(4) - vsin(p-2)*vsin(4)
      vsin(p+1) = vcos(p-2)*vsin(4) + vsin(p-2)*vcos(4)
    case (0)
      vcos(p) = vcos(p-1)*cphi - vsin(p-1)*sphi
      vsin(p) = vcos(p-1)*sphi + vsin(p-1)*cphi
      vcos(p+1) = vcos(p-1)*vcos(3) - vsin(p-1)*vsin(3)
      vsin(p+1) = vcos(p-1)*vsin(3) + vsin(p-1)*vcos(3)
    case (1)
      vcos(p+1) = vcos(p)*cphi - vsin(p)*sphi
      vsin(p+1) = vcos(p)*sphi + vsin(p)*cphi
    end select
  end subroutine trgev

  !> @brief Associated Legendre polynomials (non-negative m). Verbatim replica
  !>        of ddx_harmonics.f90::polleg_work.
  subroutine polleg_work(ctheta, stheta, p, vplm, work)
    integer,  intent(in)  :: p
    real(dp), intent(in)  :: ctheta, stheta
    real(dp), intent(out) :: vplm((p+1)**2)
    real(dp), intent(out) :: work(p+1)

    integer  :: m, l
    real(dp) :: tmp1, tmp2, tmp3, tmp4, pl2m, pl1m, plm, pmm

    select case (p)
    case (0)
      vplm(1) = one
      return
    case (1)
      vplm(1) = one
      vplm(3) = ctheta
      vplm(4) = -stheta
      return
    case (2)
      vplm(1) = one
      vplm(3) = ctheta
      vplm(4) = -stheta
      tmp1 = three * stheta
      tmp2 = ctheta * ctheta
      vplm(7) = 1.5d0*tmp2 - pt5
      vplm(8) = -tmp1 * ctheta
      vplm(9) = tmp1 * stheta
      return
    case (3)
      vplm(1) = one
      vplm(3) = ctheta
      vplm(4) = -stheta
      tmp1 = three * stheta
      tmp2 = ctheta * ctheta
      vplm(7) = 1.5d0*tmp2 - pt5
      vplm(8) = -tmp1 * ctheta
      vplm(9) = tmp1 * stheta
      tmp3 = 2.5d0 * tmp2 - pt5
      vplm(13) = tmp3*ctheta - ctheta
      vplm(14) = -tmp1 * tmp3
      tmp4 = -5d0 * stheta
      vplm(15) = tmp4 * vplm(8)
      vplm(16) = tmp4 * vplm(9)
      return
    end select
    do l = 1, p
      work(l) = dble(2*l-1) * ctheta
    end do
    pmm = one
    do m = 0, p-1
      pl2m = pmm
      vplm((m+1)**2) = pmm
      tmp1 = dble(2*m+1) * pmm
      pmm = -stheta * tmp1
      pl1m = ctheta * tmp1
      vplm((m+1)*(m+3)) = pl1m
      do l = m+2, p
        plm = work(l)*pl1m - dble(l+m-1)*pl2m
        plm = plm / dble(l-m)
        vplm(l*l+l+1+m) = plm
        pl2m = pl1m
        pl1m = plm
      end do
    end do
    vplm((p+1)**2) = pmm
  end subroutine polleg_work

end module solvent_pcm_harmonics
