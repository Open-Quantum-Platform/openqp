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
  public :: cc_uhf_ccsd_t

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

!###############################################################################

!> @brief Spin-orbital UCCSD with the perturbative triples correction.
!>
!> Equations: Stanton, Gauss, Watts, Bartlett, J. Chem. Phys. 94, 4334 (1991),
!> in the spin-orbital form tabulated by Crawford and Schaefer,
!> Rev. Comput. Chem. 14, 33 (2000).  Written for a SEMICANONICAL reference,
!> where the Fock matrix is diagonal in both the occupied-occupied and
!> virtual-virtual blocks; every f_ae, f_mi and f_me term therefore drops out
!> and only the orbital energies survive in the denominators.
!>
!> Straightforward loop form.  This is the correctness reference for the
!> open-shell path, not the performance path -- it is O(n^6) in scalar loops
!> against the closed-shell code's DGEMMs, and the spin-orbital tensors are
!> sixteen times the spatial ones.  See docs/ccsd_t_open_shell_plan.md.
subroutine cc_uhf_ccsd_t(nso, nocc, eso, g, maxit, conv, do_triples, &
                         e_ccsd, e_t, converged, niter)

  integer,  intent(in)  :: nso            !< total spin orbitals
  integer,  intent(in)  :: nocc           !< occupied spin orbitals
  real(dp), intent(in)  :: eso(nso)       !< spin-orbital energies, occ first
  real(dp), intent(in)  :: g(nso,nso,nso,nso)   !< <pq||rs>
  integer,  intent(in)  :: maxit
  real(dp), intent(in)  :: conv
  logical,  intent(in)  :: do_triples
  real(dp), intent(out) :: e_ccsd, e_t
  logical,  intent(out) :: converged
  integer,  intent(out) :: niter

  real(dp), allocatable :: t1(:,:), t2(:,:,:,:), t1n(:,:), t2n(:,:,:,:)
  real(dp), allocatable :: tau(:,:,:,:), taut(:,:,:,:)
  real(dp), allocatable :: Fae(:,:), Fmi(:,:), Fme(:,:)
  real(dp), allocatable :: Wmnij(:,:,:,:), Wabef(:,:,:,:), Wmbej(:,:,:,:)
  real(dp), allocatable :: gmnef(:,:,:,:)
  integer :: no2, nv2
  integer  :: no, nv, i, j, k, a, b, c, m, n, e, f_, it
  real(dp) :: s_, d, eold, rms, dia, dijab

  no = nocc
  nv = nso - nocc
  no2 = no*no
  nv2 = nv*nv
  e_ccsd = 0.0_dp; e_t = 0.0_dp; converged = .false.; niter = 0
  if (no <= 0 .or. nv <= 0) return

  allocate(t1(no,nv), t2(no,no,nv,nv), t1n(no,nv), t2n(no,no,nv,nv))
  allocate(tau(no,no,nv,nv), taut(no,no,nv,nv))
  allocate(Fae(nv,nv), Fmi(no,no), Fme(no,nv))
  allocate(Wmnij(no,no,no,no), Wabef(nv,nv,nv,nv), Wmbej(no,nv,nv,no))

  ! <mn||ef> over the correlated blocks.  It is constant, it is what every
  ! O(n^6) contraction below multiplies, and slicing it out of g once turns
  ! those contractions into single DGEMMs on contiguous memory.
  allocate(gmnef(no,no,nv,nv))
  do b = 1, nv; do a = 1, nv; do j = 1, no; do i = 1, no
    gmnef(i,j,a,b) = g(i,j,no+a,no+b)
  end do; end do; end do; end do

  ! MP1 amplitudes: t1 vanishes for a semicanonical reference.
  t1 = 0.0_dp
  do b = 1, nv; do a = 1, nv; do j = 1, no; do i = 1, no
    d = eso(i) + eso(j) - eso(no+a) - eso(no+b)
    t2(i,j,a,b) = g(i,j,no+a,no+b) / d
  end do; end do; end do; end do

  eold = 0.0_dp

  do it = 1, maxit
    niter = it

    ! --- tau and tau-tilde ---------------------------------------------------
    do b = 1, nv; do a = 1, nv; do j = 1, no; do i = 1, no
      d = t1(i,a)*t1(j,b) - t1(i,b)*t1(j,a)
      tau (i,j,a,b) = t2(i,j,a,b) + d
      taut(i,j,a,b) = t2(i,j,a,b) + 0.5_dp*d
    end do; end do; end do; end do

    ! --- one-particle intermediates (f is diagonal: no f_ae/f_mi/f_me) -------
    do e = 1, nv; do a = 1, nv
      s_ = 0.0_dp
      do f_ = 1, nv; do m = 1, no
        s_ = s_ + t1(m,f_)*g(m,no+a,no+f_,no+e)
      end do; end do
      do f_ = 1, nv; do n = 1, no; do m = 1, no
        s_ = s_ - 0.5_dp*taut(m,n,a,f_)*g(m,n,no+e,no+f_)
      end do; end do; end do
      Fae(a,e) = s_
    end do; end do

    do i = 1, no; do m = 1, no
      s_ = 0.0_dp
      do e = 1, nv; do n = 1, no
        s_ = s_ + t1(n,e)*g(m,n,i,no+e)
      end do; end do
      do f_ = 1, nv; do e = 1, nv; do n = 1, no
        s_ = s_ + 0.5_dp*taut(i,n,e,f_)*g(m,n,no+e,no+f_)
      end do; end do; end do
      Fmi(m,i) = s_
    end do; end do

    do e = 1, nv; do m = 1, no
      s_ = 0.0_dp
      do f_ = 1, nv; do n = 1, no
        s_ = s_ + t1(n,f_)*g(m,n,no+e,no+f_)
      end do; end do
      Fme(m,e) = s_
    end do; end do

    ! --- two-particle intermediates -----------------------------------------
    do j = 1, no; do i = 1, no; do n = 1, no; do m = 1, no
      s_ = g(m,n,i,j)
      do e = 1, nv
        s_ = s_ + t1(j,e)*g(m,n,i,no+e) - t1(i,e)*g(m,n,j,no+e)
      end do
      Wmnij(m,n,i,j) = s_
    end do; end do; end do; end do
    ! Wmnij(mn,ij) += 1/4 <mn||ef> tau(ij,ef)
    call dgemm('n','t', no2, no2, nv2, 0.25_dp, gmnef, no2, tau, no2, &
               1.0_dp, Wmnij, no2)

    do f_ = 1, nv; do e = 1, nv; do b = 1, nv; do a = 1, nv
      s_ = g(no+a,no+b,no+e,no+f_)
      do m = 1, no
        s_ = s_ - t1(m,b)*g(no+a,m,no+e,no+f_) + t1(m,a)*g(no+b,m,no+e,no+f_)
      end do
      Wabef(a,b,e,f_) = s_
    end do; end do; end do; end do
    ! Wabef(ab,ef) += 1/4 tau(mn,ab)^T <mn||ef>
    call dgemm('t','n', nv2, nv2, no2, 0.25_dp, tau, no2, gmnef, no2, &
               1.0_dp, Wabef, nv2)

    do j = 1, no; do e = 1, nv; do b = 1, nv; do m = 1, no
      s_ = g(m,no+b,no+e,j)
      do f_ = 1, nv
        s_ = s_ + t1(j,f_)*g(m,no+b,no+e,no+f_)
      end do
      do n = 1, no
        s_ = s_ - t1(n,b)*g(m,n,no+e,j)
      end do
      do f_ = 1, nv; do n = 1, no
        s_ = s_ - (0.5_dp*t2(j,n,f_,b) + t1(j,f_)*t1(n,b))*g(m,n,no+e,no+f_)
      end do; end do
      Wmbej(m,b,e,j) = s_
    end do; end do; end do; end do

    ! --- T1 ------------------------------------------------------------------
    do a = 1, nv; do i = 1, no
      s_ = 0.0_dp
      do e = 1, nv
        s_ = s_ + t1(i,e)*Fae(a,e)
      end do
      do m = 1, no
        s_ = s_ - t1(m,a)*Fmi(m,i)
      end do
      do e = 1, nv; do m = 1, no
        s_ = s_ + t2(i,m,a,e)*Fme(m,e)
      end do; end do
      do f_ = 1, nv; do n = 1, no
        s_ = s_ - t1(n,f_)*g(n,no+a,i,no+f_)
      end do; end do
      do f_ = 1, nv; do e = 1, nv; do m = 1, no
        s_ = s_ - 0.5_dp*t2(i,m,e,f_)*g(m,no+a,no+e,no+f_)
      end do; end do; end do
      do e = 1, nv; do n = 1, no; do m = 1, no
        s_ = s_ - 0.5_dp*t2(m,n,a,e)*g(n,m,no+e,i)
      end do; end do; end do
      dia = eso(i) - eso(no+a)
      t1n(i,a) = s_ / dia
    end do; end do

    ! --- T2 ------------------------------------------------------------------
    do b = 1, nv; do a = 1, nv; do j = 1, no; do i = 1, no
      s_ = g(i,j,no+a,no+b)

      ! P(ab) sum_e t2_ijae (Fbe - 1/2 sum_m t1_mb Fme)
      do e = 1, nv
        d = 0.0_dp
        do m = 1, no
          d = d + t1(m,b)*Fme(m,e)
        end do
        s_ = s_ + t2(i,j,a,e)*(Fae(b,e) - 0.5_dp*d)
        d = 0.0_dp
        do m = 1, no
          d = d + t1(m,a)*Fme(m,e)
        end do
        s_ = s_ - t2(i,j,b,e)*(Fae(a,e) - 0.5_dp*d)
      end do

      ! -P(ij) sum_m t2_imab (Fmj + 1/2 sum_e t1_je Fme)
      do m = 1, no
        d = 0.0_dp
        do e = 1, nv
          d = d + t1(j,e)*Fme(m,e)
        end do
        s_ = s_ - t2(i,m,a,b)*(Fmi(m,j) + 0.5_dp*d)
        d = 0.0_dp
        do e = 1, nv
          d = d + t1(i,e)*Fme(m,e)
        end do
        s_ = s_ + t2(j,m,a,b)*(Fmi(m,i) + 0.5_dp*d)
      end do

      do n = 1, no; do m = 1, no
        s_ = s_ + 0.5_dp*tau(m,n,a,b)*Wmnij(m,n,i,j)
      end do; end do

      ! P(ij)P(ab) [ t2_imae Wmbej - t1_ie t1_ma <mb||ej> ]
      do e = 1, nv; do m = 1, no
        s_ = s_ + t2(i,m,a,e)*Wmbej(m,b,e,j) - t1(i,e)*t1(m,a)*g(m,no+b,no+e,j)
        s_ = s_ - t2(j,m,a,e)*Wmbej(m,b,e,i) + t1(j,e)*t1(m,a)*g(m,no+b,no+e,i)
        s_ = s_ - t2(i,m,b,e)*Wmbej(m,a,e,j) + t1(i,e)*t1(m,b)*g(m,no+a,no+e,j)
        s_ = s_ + t2(j,m,b,e)*Wmbej(m,a,e,i) - t1(j,e)*t1(m,b)*g(m,no+a,no+e,i)
      end do; end do

      do e = 1, nv
        s_ = s_ + t1(i,e)*g(no+a,no+b,no+e,j) - t1(j,e)*g(no+a,no+b,no+e,i)
      end do
      do m = 1, no
        s_ = s_ - t1(m,a)*g(m,no+b,i,j) + t1(m,b)*g(m,no+a,i,j)
      end do

      t2n(i,j,a,b) = s_
    end do; end do; end do; end do

    ! Ladder: t2n(ij,ab) += 1/2 tau(ij,ef) Wabef(ab,ef).  This is the
    ! O(o^2 v^4) bottleneck and the one term worth a DGEMM above all others.
    call dgemm('n','t', no2, nv2, nv2, 0.5_dp, tau, no2, Wabef, nv2, &
               1.0_dp, t2n, no2)

    do b = 1, nv; do a = 1, nv; do j = 1, no; do i = 1, no
      dijab = eso(i) + eso(j) - eso(no+a) - eso(no+b)
      t2n(i,j,a,b) = t2n(i,j,a,b) / dijab
    end do; end do; end do; end do

    rms = sqrt(sum((t1n-t1)**2) + sum((t2n-t2)**2))
    t1 = t1n; t2 = t2n

    e_ccsd = 0.0_dp
    do b = 1, nv; do a = 1, nv; do j = 1, no; do i = 1, no
      e_ccsd = e_ccsd + 0.25_dp*g(i,j,no+a,no+b)*t2(i,j,a,b) &
                      + 0.5_dp *g(i,j,no+a,no+b)*t1(i,a)*t1(j,b)
    end do; end do; end do; end do

    if (rms < conv .and. abs(e_ccsd-eold) < conv) then
      converged = .true.
      exit
    end if
    eold = e_ccsd
  end do

  ! --- (T) -----------------------------------------------------------------
  if (do_triples) call triples(e_t)

  deallocate(t1, t2, t1n, t2n, tau, taut, Fae, Fmi, Fme, Wmnij, Wabef, Wmbej, gmnef)

contains

  !> Perturbative triples in the spin-orbital form:
  !>   E = 1/36 sum t3c * D * (t3c + t3d)
  !> with t3c the connected and t3d the disconnected triple, both already
  !> divided by their denominator.
  subroutine triples(et)
    real(dp), intent(out) :: et
    integer  :: i, j, k, a, b, c
    real(dp) :: dd, tc, td

    et = 0.0_dp
    !$omp parallel do collapse(3) default(shared) &
    !$omp   private(i,j,k,a,b,c,dd,tc,td) reduction(+:et) schedule(dynamic)
    do i = 1, no
      do j = 1, no
        do k = 1, no
          do a = 1, nv
            do b = 1, nv
              do c = 1, nv
                dd = eso(i)+eso(j)+eso(k)-eso(no+a)-eso(no+b)-eso(no+c)
                if (abs(dd) < 1.0e-12_dp) cycle
                tc = conn(i,j,k,a,b,c) / dd
                td = disc(i,j,k,a,b,c) / dd
                et = et + tc*dd*(tc+td) / 36.0_dp
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine triples

  !> Disconnected triple numerator: P(i/jk) P(a/bc) t1_ia <jk||bc>.
  pure real(dp) function disc(i,j,k,a,b,c) result(v)
    integer, intent(in) :: i,j,k,a,b,c
    v =   p1(i,j,k,a,b,c) - p1(j,i,k,a,b,c) - p1(k,j,i,a,b,c)
  end function disc

  !> P(a/bc) t1_ia <jk||bc> for a fixed occupied ordering.
  pure real(dp) function p1(i,j,k,a,b,c) result(v)
    integer, intent(in) :: i,j,k,a,b,c
    v =   t1s(i,a)*gs(j,k,b,c) - t1s(i,b)*gs(j,k,a,c) - t1s(i,c)*gs(j,k,b,a)
  end function p1

  !> Connected triple numerator:
  !>   P(i/jk) P(a/bc) [ sum_e t2_jkae <ei||bc> - sum_m t2_imbc <ma||jk> ]
  pure real(dp) function conn(i,j,k,a,b,c) result(v)
    integer, intent(in) :: i,j,k,a,b,c
    v =   p2(i,j,k,a,b,c) - p2(j,i,k,a,b,c) - p2(k,j,i,a,b,c)
  end function conn

  pure real(dp) function p2(i,j,k,a,b,c) result(v)
    integer, intent(in) :: i,j,k,a,b,c
    v =   q2(i,j,k,a,b,c) - q2(i,j,k,b,a,c) - q2(i,j,k,c,b,a)
  end function p2

  pure real(dp) function q2(i,j,k,a,b,c) result(v)
    integer, intent(in) :: i,j,k,a,b,c
    integer :: e, m
    v = 0.0_dp
    do e = 1, nv
      v = v + t2s(j,k,a,e)*g(no+e,i,no+b,no+c)
    end do
    do m = 1, no
      v = v - t2s(i,m,b,c)*g(m,no+a,j,k)
    end do
  end function q2

  pure real(dp) function t1s(i,a) result(v)
    integer, intent(in) :: i,a
    v = t1(i,a)
  end function t1s

  pure real(dp) function t2s(i,j,a,b) result(v)
    integer, intent(in) :: i,j,a,b
    v = t2(i,j,a,b)
  end function t2s

  pure real(dp) function gs(j,k,b,c) result(v)
    integer, intent(in) :: j,k,b,c
    v = g(j,k,no+b,no+c)
  end function gs

end subroutine cc_uhf_ccsd_t

end module cc_uhf_lib
