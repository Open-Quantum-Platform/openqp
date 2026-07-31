!> @file cc_uhf_lib.F90
!>
!> @brief Spin-orbital foundation for open-shell coupled cluster.
!>
!> Assembles antisymmetrised spin-orbital integrals <pq||rs> and the
!> spin-orbital MP2 energy from spin-blocked spatial MO integrals.  This is the
!> layer the open-shell CCSD equations sit on; MP2 is here because it is the
!> cheapest end-to-end check that the assembly and the orbital ordering are
!> right, and it shares every input with the coupled-cluster code to come.
!>
!> Spin ordering is INTERLEAVED: spin orbital 2p-1 is the alpha component of
!> spatial orbital p and 2p its beta component.  Interleaving keeps the
!> occupied spin orbitals contiguous whenever the alpha and beta occupations
!> are equal, and puts them in ascending orbital-energy order for the closed
!> shell case, so a UHF reference on a closed-shell system reduces to the RHF
!> answer without any special casing.
!>
!> Memory: <pq||rs> is (2*nmo)^4, sixteen times the spatial tensor.  That is
!> acceptable for validation-sized systems and for small open-shell cases, but
!> it is NOT the production layout -- see docs/ccsd_t_open_shell_plan.md for
!> the spin-integrated blocking that the full implementation should use.
!> cc_uhf_spinorb_gb reports the cost so callers can refuse in advance.
module cc_uhf_lib

  use precision, only: dp
  implicit none

  private
  public :: cc_uhf_spinorb_build
  public :: cc_uhf_spinorb_gb
  public :: cc_uhf_mp2

contains

!###############################################################################

!> Gigabytes needed for the spin-orbital integral tensor of @p nmo spatial MOs.
pure function cc_uhf_spinorb_gb(nmo) result(gb)
  integer, intent(in) :: nmo
  real(dp) :: gb
  gb = (2.0_dp*real(nmo,dp))**4 * 8.0_dp / 1.073741824e9_dp
end function cc_uhf_spinorb_gb

!###############################################################################

!> Build antisymmetrised spin-orbital integrals in physicist notation,
!>   g(p,q,r,s) = <pq||rs> = (pr|qs) - (ps|qr),
!> from spatial MO integrals supplied in chemist notation for the three spin
!> cases.  @p eri_aa and @p eri_bb are the same-spin spatial tensors and
!> @p eri_ab(p,q,r,s) = (pq|rs) with (pq) alpha and (rs) beta.
!>
!> The Coulomb operator is spin free, so (pr|qs) survives only when p,r share a
!> spin and q,s share a spin; that is the whole content of the spin tests below.
subroutine cc_uhf_spinorb_build(nmo, eri_aa, eri_bb, eri_ab, g)

  integer, intent(in) :: nmo
  real(dp), intent(in) :: eri_aa(nmo,nmo,nmo,nmo)
  real(dp), intent(in) :: eri_bb(nmo,nmo,nmo,nmo)
  real(dp), intent(in) :: eri_ab(nmo,nmo,nmo,nmo)
  real(dp), intent(out) :: g(2*nmo,2*nmo,2*nmo,2*nmo)

  integer :: p, q, r, s, ps, qs, rs_, ss
  integer :: pp, qq, rr, ssp
  real(dp) :: coul, exch

  !$omp parallel do collapse(4) default(shared) &
  !$omp   private(p,q,r,s,ps,qs,rs_,ss,pp,qq,rr,ssp,coul,exch) schedule(static)
  do s = 1, 2*nmo
    do r = 1, 2*nmo
      do q = 1, 2*nmo
        do p = 1, 2*nmo
          ! spatial index and spin (1 = alpha, 2 = beta) of each spin orbital
          pp = (p+1)/2; ps  = 2 - mod(p,2)
          qq = (q+1)/2; qs  = 2 - mod(q,2)
          rr = (r+1)/2; rs_ = 2 - mod(r,2)
          ssp= (s+1)/2; ss  = 2 - mod(s,2)

          coul = 0.0_dp
          if (ps == rs_ .and. qs == ss) coul = chem(pp, rr, qq, ssp, ps, qs)

          exch = 0.0_dp
          if (ps == ss .and. qs == rs_) exch = chem(pp, ssp, qq, rr, ps, qs)

          g(p,q,r,s) = coul - exch
        end do
      end do
    end do
  end do
  !$omp end parallel do

contains

  !> Spatial chemist integral (ab|cd) where (a,b) carry spin @p s1 and (c,d)
  !> carry spin @p s2.  The mixed case is stored alpha-pair first, so an
  !> all-beta-first request is served by the (cd|ab) symmetry of the tensor.
  pure real(dp) function chem(a, b, c, d, s1, s2) result(v)
    integer, intent(in) :: a, b, c, d, s1, s2
    if (s1 == 1 .and. s2 == 1) then
      v = eri_aa(a,b,c,d)
    else if (s1 == 2 .and. s2 == 2) then
      v = eri_bb(a,b,c,d)
    else if (s1 == 1) then
      v = eri_ab(a,b,c,d)
    else
      v = eri_ab(c,d,a,b)
    end if
  end function chem

end subroutine cc_uhf_spinorb_build

!###############################################################################

!> Spin-orbital MP2 correlation energy,
!>   E = 1/4 sum_ijab |<ij||ab>|^2 / (e_i + e_j - e_a - e_b).
!>
!> @p eso holds the spin-orbital energies in the same interleaved order as
!> @p g, and @p nocc is the total number of occupied spin orbitals
!> (nelec_a + nelec_b of the correlated space).
!>
!> This exists to validate the assembly above: with a correct build it
!> reproduces UMP2 exactly, and on a closed-shell system with a UHF reference
!> it reproduces RMP2.  It is not a production MP2 -- source/mp2_lib.F90 is.
pure function cc_uhf_mp2(nso, nocc, eso, g) result(e_mp2)

  integer, intent(in) :: nso, nocc
  real(dp), intent(in) :: eso(nso)
  real(dp), intent(in) :: g(nso,nso,nso,nso)
  real(dp) :: e_mp2

  integer :: i, j, a, b
  real(dp) :: denom

  e_mp2 = 0.0_dp
  do i = 1, nocc
    do j = 1, nocc
      do a = nocc+1, nso
        do b = nocc+1, nso
          denom = eso(i) + eso(j) - eso(a) - eso(b)
          if (abs(denom) < 1.0e-12_dp) cycle
          e_mp2 = e_mp2 + 0.25_dp * g(i,j,a,b)**2 / denom
        end do
      end do
    end do
  end do

end function cc_uhf_mp2

end module cc_uhf_lib
