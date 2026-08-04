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
!> it is NOT the production layout: a full implementation should keep the
!> spin-integrated blocks (aa, ab, bb) separately and block the ladder over
!> the virtual index, the way cc_lib does for the closed shell.
!> cc_uhf_spinorb_gb reports the cost so callers can refuse in advance.
module cc_uhf_lib

  use precision, only: dp
  use blas_thread, only: blas_thread_count, blas_thread_set
  use, intrinsic :: iso_c_binding, only: c_int64_t
  implicit none

  !> Collapsed trip count below which a parallel region costs more than the
  !> loop it wraps.  Measured on the open-shell solver: without this guard the
  !> small cases run several times slower on four threads than on one.
  integer, parameter :: PAR_MIN = 4096

  !> Same idea for the regions whose cost is a product of four or six extents
  !> rather than a plain trip count: below this many elementary operations the
  !> region is not worth opening.
  real(dp), parameter :: PAR_WORK = 1.0e5_dp

  private
  public :: cc_uhf_spinorb_build
  public :: cc_uhf_spinorb_gb
  public :: cc_uhf_peak_gb
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

!> Peak gigabytes for the whole open-shell path, not just its largest array.
!>
!> Two stages compete for the peak and the larger one wins:
!>
!>   * assembly -- the three spatial spin-block tensors (3 nmo^4) are still
!>     alive while the spin-orbital tensor (nso^4) is being filled;
!>   * solution -- the spin-orbital tensor stays, and the solver adds the
!>     ladder intermediate Wabef (nv^4), Wmbej (no^2 nv^2), Wmnij (no^4),
!>     gmnef and the six amplitude arrays (~7 no^2 nv^2), the ten (ov)^2
!>     panels the ring and Wmbej contractions are grouped into, and a DIIS
!>     history of 2*ndiis vectors of no*nv + no^2 nv^2 each;
!>   * triples -- the tensor plus the four reordered panels build_q reads,
!>     the largest being nv^3 no, PLUS one Q(nv,nv,nv) panel per OpenMP
!>     thread.  That last term is the only one that grows with the machine
!>     rather than with the problem, and on a virtual-heavy case run wide it
!>     is not a rounding error: at nv = 200 and 64 threads it is 3.8 GB.
!>     Leaving it out let a job pass the guard and still be OOM-killed once
!>     the triples opened their parallel region.
!>
!> Reporting only nso^4 understated the real requirement by roughly a factor
!> of three, which is worse than useless in a guard: it waves through jobs
!> that then die in the allocator.
pure function cc_uhf_peak_gb(nmo, nocc, ndiis, nthreads) result(gb)
  integer, intent(in) :: nmo    !< correlated spatial MOs
  integer, intent(in) :: nocc   !< occupied spin orbitals
  integer, intent(in) :: ndiis  !< DIIS subspace size (0 = no DIIS)
  !> OpenMP threads the (T) loop will run with; each takes its own nv^3 panel.
  integer, intent(in), optional :: nthreads
  real(dp) :: gb
  real(dp) :: nso, no, nv, assembly, solve, trip, amp
  integer :: nthr

  nso = 2.0_dp*real(nmo,dp)
  no  = real(nocc,dp)
  nv  = nso - no

  assembly = 3.0_dp*real(nmo,dp)**4 + nso**4

  amp   = no*nv + no**2*nv**2
  solve = nso**4 + nv**4 + no**4 + 2.0_dp*no**2*nv**2 &
        + 7.0_dp*no**2*nv**2 + 10.0_dp*no**2*nv**2 &
        + 2.0_dp*real(max(ndiis,0),dp)*amp

  nthr = 1
  if (present(nthreads)) nthr = max(1, nthreads)
  trip = nso**4 + nv**3*no + no**2*nv**2 + no**2*nv**2 + no**3*nv &
       + real(nthr,dp)*nv**3

  gb = max(assembly, max(solve, trip)) * 8.0_dp / 1.073741824e9_dp
end function cc_uhf_peak_gb

!###############################################################################

!> Build antisymmetrised spin-orbital integrals in physicist notation,
!>   g(p,q,r,s) = <pq||rs> = (pr|qs) - (ps|qr),
!> from spatial MO integrals supplied in chemist notation for the three spin
!> cases.  @p eri_aa and @p eri_bb are the same-spin spatial tensors and
!> @p eri_ab(p,q,r,s) = (pq|rs) with (pq) alpha and (rs) beta.
!>
!> The Coulomb operator is spin free, so (pr|qs) survives only when p,r share a
!> spin and q,s share a spin; that is the whole content of the spin tests below.
!>
!> @param[in] ord  optional permutation, ord(i) = the interleaved spin orbital
!>                 that lands at position i.  The tensor is written directly in
!>                 that order.  Building it permuted costs nothing here, while
!>                 permuting afterwards would mean `g = g(ord,ord,ord,ord)` --
!>                 a whole second copy of the largest array in the run, which
!>                 doubles the peak of the open-shell path.
subroutine cc_uhf_spinorb_build(nmo, eri_aa, eri_bb, eri_ab, g, ord)

  integer, intent(in) :: nmo
  real(dp), intent(in) :: eri_aa(nmo,nmo,nmo,nmo)
  real(dp), intent(in) :: eri_bb(nmo,nmo,nmo,nmo)
  real(dp), intent(in) :: eri_ab(nmo,nmo,nmo,nmo)
  real(dp), intent(out) :: g(2*nmo,2*nmo,2*nmo,2*nmo)
  integer, intent(in), optional :: ord(2*nmo)

  integer :: p, q, r, s, ps, qs, rs_, ss
  integer :: pp, qq, rr, ssp
  integer :: omap(2*nmo)
  real(dp) :: coul, exch

  ! Resolve the permutation once so the hot loop indexes an array rather than
  ! branching on present() every element.
  if (present(ord)) then
    omap = ord
  else
    do p = 1, 2*nmo
      omap(p) = p
    end do
  end if

  !$omp parallel do collapse(4) default(shared) &
  !$omp   private(p,q,r,s,ps,qs,rs_,ss,pp,qq,rr,ssp,coul,exch) schedule(static) &
  !$omp   if(real(2*nmo,dp)**4 > PAR_WORK)
  do s = 1, 2*nmo
    do r = 1, 2*nmo
      do q = 1, 2*nmo
        do p = 1, 2*nmo
          ! spatial index and spin (1 = alpha, 2 = beta) of each spin orbital,
          ! taken from the orbital that belongs at this position
          pp = (omap(p)+1)/2; ps  = 2 - mod(omap(p),2)
          qq = (omap(q)+1)/2; qs  = 2 - mod(omap(q),2)
          rr = (omap(r)+1)/2; rs_ = 2 - mod(omap(r),2)
          ssp= (omap(s)+1)/2; ss  = 2 - mod(omap(s),2)

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
!> sixteen times the spatial ones.  Use it for small open-shell systems.
subroutine cc_uhf_ccsd_t(nso, nocc, eso, g, maxit, conv, do_triples, &
                         e_ccsd, e_t, converged, niter, fov, &
                         time_ccsd, time_triples, ndiis_in)

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
  !> Occupied-virtual Fock block in the same spin-orbital basis.  Zero for a
  !> canonical UHF reference; for ROHF it survives semicanonicalisation and
  !> dropping it costs ~5e-3 Hartree.  Absent means zero.
  real(dp), optional, intent(in) :: fov(nocc, nso-nocc)
  !> Wall-clock seconds in the CCSD iterations and in the (T) correction,
  !> reported separately -- (T) is the more expensive half on anything but the
  !> smallest cases, and folding it into the iteration time hides that.
  real(dp), optional, intent(out) :: time_ccsd, time_triples
  !> DIIS subspace size; 0 disables DIIS.  Absent keeps the default of 8.
  integer, optional, intent(in) :: ndiis_in

  real(dp), allocatable :: t1(:,:), t2(:,:,:,:), t1n(:,:), t2n(:,:,:,:)
  real(dp), allocatable :: tau(:,:,:,:), taut(:,:,:,:)
  real(dp), allocatable :: Fae(:,:), Fmi(:,:), Fme(:,:)
  ! Fock intermediates with the t1-dependent inner sums folded in (see T2).
  real(dp), allocatable :: FaeT(:,:), FmiT(:,:)
  real(dp), allocatable :: Wmnij(:,:,:,:), Wabef(:,:,:,:), Wmbej(:,:,:,:)
  real(dp), allocatable :: gmnef(:,:,:,:), f_ov(:,:)
  ! Reordered panels for the (T) correction, packed once by triples() and read
  ! by build_q().  Declared here so both see them by host association:
  !   gvv(e,bc,i) = <ei||bc>       t2o(m,bc,i) = t2(i,m,b,c)
  !   t2v(a,e,jk) = t2(j,k,a,e)    gov(m,a,jk) = <ma||jk>
  real(dp), allocatable :: gvv(:,:,:), t2o(:,:,:), t2v(:,:,:), gov(:,:,:)
  real(dp), allocatable :: gvvp(:,:,:)
  ! Ring-term panels, rebuilt each iteration (see the T2 section).
  real(dp), allocatable :: t2ring(:,:), wring(:,:), zring(:,:)
  real(dp), allocatable :: wjbnf(:,:), gnfme(:,:), wjbme(:,:)
  real(dp), allocatable :: gring(:,:,:,:), hring(:,:,:,:), hringp(:,:,:,:)
  real(dp), allocatable :: yring(:,:,:,:)
  integer(c_int64_t) :: nblas_save
  ! DIIS history over the concatenated (t1,t2) vector and its residual.
  real(dp), allocatable :: dv(:,:), de(:,:), bmat(:,:), rhs(:)
  integer, allocatable :: ipiv(:)
  integer :: no2, nv2, namp, ndiis, ndim, nvec, pos, ii, jj, info
  integer  :: no, nv, i, j, k, a, b, c, m, n, e, f_, it
  real(dp) :: s_, d, eold, rms, dia, dijab
  real(dp) :: tw0, tw1

  no = nocc
  nv = nso - nocc
  no2 = no*no
  nv2 = nv*nv

  ! This routine drives its own OpenMP regions and calls BLAS from between
  ! them.  A threaded BLAS underneath opens a second pool that spin-waits
  ! through the other's regions, and its per-call overhead swamps the small
  ! GEMMs the open-shell equations make.  Measured on four threads: with
  ! threaded BLAS the small cases run SLOWER than serial (O2 triplet 0.17 ->
  ! 0.28 s) and the large one gains only 2.7x; pinned, O2 speeds up and
  ! BH/cc-pVDZ reaches 3.6x.  Restored on exit.
  nblas_save = blas_thread_count()
  if (nblas_save > 0) call blas_thread_set(1_c_int64_t)
  e_ccsd = 0.0_dp; e_t = 0.0_dp; converged = .false.; niter = 0
  if (present(time_ccsd)) time_ccsd = 0.0_dp
  if (present(time_triples)) time_triples = 0.0_dp
  if (no <= 0 .or. nv <= 0) return

  call cc_uhf_wall_time(tw0)

  allocate(t1(no,nv), t2(no,no,nv,nv), t1n(no,nv), t2n(no,no,nv,nv))
  allocate(tau(no,no,nv,nv), taut(no,no,nv,nv))
  allocate(Fae(nv,nv), Fmi(no,no), Fme(no,nv), FaeT(nv,nv), FmiT(no,no))
  allocate(Wmnij(no,no,no,no), Wabef(nv,nv,nv,nv), Wmbej(no,nv,nv,no))
  allocate(t2ring(no*nv, no*nv), wring(no*nv, no*nv), zring(no*nv, no*nv))
  allocate(wjbnf(no*nv, no*nv), gnfme(no*nv, no*nv), wjbme(no*nv, no*nv))
  allocate(gring(nv,no,nv,no), hring(no,no,nv,no), hringp(no,no,nv,no), &
           yring(nv,no,nv,no))
  allocate(f_ov(no,nv))
  f_ov = 0.0_dp
  if (present(fov)) f_ov = fov

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

  ! DIIS.  Without it these equations take 40-50 iterations; the closed-shell
  ! solver reaches the same tolerance in 15-20 with it.  Eight vectors is the
  ! same subspace size cc_lib defaults to, and [cc] ndiis selects it here too
  ! -- the keyword used to apply to the closed-shell path only.
  namp = no*nv + no*no*nv*nv
  ndiis = 8
  if (present(ndiis_in)) ndiis = max(0, ndiis_in)
  if (ndiis > 0) allocate(dv(namp,ndiis), de(namp,ndiis))
  nvec = 0

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
    !$omp parallel do collapse(2) default(shared) private(a,e,m,n,f_,s_) schedule(static) if(nv2 > PAR_MIN)
    do e = 1, nv; do a = 1, nv
      s_ = 0.0_dp
      do f_ = 1, nv; do m = 1, no
        s_ = s_ + t1(m,f_)*g(m,no+a,no+f_,no+e)
      end do; end do
      do f_ = 1, nv; do n = 1, no; do m = 1, no
        s_ = s_ - 0.5_dp*taut(m,n,a,f_)*g(m,n,no+e,no+f_)
      end do; end do; end do
      do m = 1, no
        s_ = s_ - 0.5_dp*f_ov(m,e)*t1(m,a)
      end do
      Fae(a,e) = s_
    end do; end do
    !$omp end parallel do

    !$omp parallel do collapse(2) default(shared) private(i,m,e,n,f_,s_) schedule(static) if(no2 > PAR_MIN)
    do i = 1, no; do m = 1, no
      s_ = 0.0_dp
      do e = 1, nv; do n = 1, no
        s_ = s_ + t1(n,e)*g(m,n,i,no+e)
      end do; end do
      do f_ = 1, nv; do e = 1, nv; do n = 1, no
        s_ = s_ + 0.5_dp*taut(i,n,e,f_)*g(m,n,no+e,no+f_)
      end do; end do; end do
      do e = 1, nv
        s_ = s_ + 0.5_dp*t1(i,e)*f_ov(m,e)
      end do
      Fmi(m,i) = s_
    end do; end do
    !$omp end parallel do

    !$omp parallel do collapse(2) default(shared) private(m,e,n,f_,s_) schedule(static) if(no*nv > PAR_MIN)
    do e = 1, nv; do m = 1, no
      s_ = f_ov(m,e)
      do f_ = 1, nv; do n = 1, no
        s_ = s_ + t1(n,f_)*g(m,n,no+e,no+f_)
      end do; end do
      Fme(m,e) = s_
    end do; end do
    !$omp end parallel do

    ! --- two-particle intermediates -----------------------------------------
    !$omp parallel do collapse(4) default(shared) private(i,j,m,n,e,s_) schedule(static) if(no2*no2 > PAR_MIN)
    do j = 1, no; do i = 1, no; do n = 1, no; do m = 1, no
      s_ = g(m,n,i,j)
      do e = 1, nv
        s_ = s_ + t1(j,e)*g(m,n,i,no+e) - t1(i,e)*g(m,n,j,no+e)
      end do
      Wmnij(m,n,i,j) = s_
    end do; end do; end do; end do
    !$omp end parallel do
    ! Wmnij(mn,ij) += 1/4 <mn||ef> tau(ij,ef)
    call dgemm('n','t', no2, no2, nv2, 0.25_dp, gmnef, no2, tau, no2, &
               1.0_dp, Wmnij, no2)

    !$omp parallel do collapse(4) default(shared) private(a,b,e,f_,m,n,s_) schedule(static) if(nv2*nv2 > PAR_MIN)
    do f_ = 1, nv; do e = 1, nv; do b = 1, nv; do a = 1, nv
      s_ = g(no+a,no+b,no+e,no+f_)
      do m = 1, no
        s_ = s_ - t1(m,b)*g(no+a,m,no+e,no+f_) + t1(m,a)*g(no+b,m,no+e,no+f_)
      end do
      Wabef(a,b,e,f_) = s_
    end do; end do; end do; end do
    !$omp end parallel do
    ! Wabef(ab,ef) += 1/4 tau(mn,ab)^T <mn||ef>
    call dgemm('t','n', nv2, nv2, no2, 0.25_dp, tau, no2, gmnef, no2, &
               1.0_dp, Wabef, nv2)

    !$omp parallel do collapse(4) default(shared) private(m,b,e,j,f_,n,s_) schedule(static) if(no2*nv2 > PAR_MIN)
    do j = 1, no; do e = 1, nv; do b = 1, nv; do m = 1, no
      s_ = g(m,no+b,no+e,j)
      do f_ = 1, nv
        s_ = s_ + t1(j,f_)*g(m,no+b,no+e,no+f_)
      end do
      do n = 1, no
        s_ = s_ - t1(n,b)*g(m,n,no+e,j)
      end do
      Wmbej(m,b,e,j) = s_
    end do; end do; end do; end do
    !$omp end parallel do

    ! The remaining Wmbej term, -sum_nf (1/2 t2_jnfb + t1_jf t1_nb) <mn||ef>,
    ! is an (o v)^3 contraction and was the second-largest scalar cost of the
    ! iteration.  Grouped as (jb,nf) x (nf,me) it is a single DGEMM.
    !$omp parallel do collapse(2) default(shared) private(j,b,n,f_,m,e) schedule(static)
    do f_ = 1, nv
      do n = 1, no
        do b = 1, nv
          do j = 1, no
            wjbnf((b-1)*no+j, (f_-1)*no+n) = 0.5_dp*t2(j,n,f_,b) + t1(j,f_)*t1(n,b)
          end do
        end do
        do e = 1, nv
          do m = 1, no
            gnfme((f_-1)*no+n, (e-1)*no+m) = gmnef(m,n,e,f_)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call dgemm('n','n', no*nv, no*nv, no*nv, -1.0_dp, wjbnf, no*nv, &
               gnfme, no*nv, 0.0_dp, wjbme, no*nv)

    !$omp parallel do collapse(4) default(shared) private(j,e,b,m) schedule(static)
    do j = 1, no; do e = 1, nv; do b = 1, nv; do m = 1, no
      Wmbej(m,b,e,j) = Wmbej(m,b,e,j) + wjbme((b-1)*no+j, (e-1)*no+m)
    end do; end do; end do; end do
    !$omp end parallel do

    ! --- T1 ------------------------------------------------------------------
    !$omp parallel do collapse(2) default(shared) private(i,a,e,m,n,f_,s_,dia) schedule(static) if(no*nv > PAR_MIN)
    do a = 1, nv; do i = 1, no
      s_ = f_ov(i,a)
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
    !$omp end parallel do

    ! --- T2 ------------------------------------------------------------------
    ! The two Fock-dressed terms carry an inner sum that does not depend on the
    ! amplitude indices being looped over -- sum_m t1(m,b) Fme(m,e) is a
    ! function of (b,e) alone, and sum_e t1(j,e) Fme(m,e) of (j,m) alone.
    ! Evaluated in place they were recomputed no^2 nv^2 times each; dressing
    ! Fae and Fmi once per iteration costs a single small DGEMM.
    FaeT = Fae
    call dgemm('t','n', nv, nv, no, -0.5_dp, t1, no, Fme, no, 1.0_dp, FaeT, nv)
    FmiT = Fmi
    call dgemm('n','t', no, no, nv, 0.5_dp, Fme, no, t1, no, 1.0_dp, FmiT, no)

    ! The ring term, P(ij)P(ab)[ t2_imae Wmbej - t1_ie t1_ma <mb||ej> ], is the
    ! dominant cost of the iteration: four permutations of an (o v)^3
    ! contraction.  Evaluated element by element inside the loop below it is
    ! ~80% of the CCSD time, so both halves are precomputed here instead.
    !
    !   Zring(ia,jb) = sum_me t2(i,m,a,e) Wmbej(m,b,e,j)   -- one DGEMM
    !   Yring(a,i,b,j) = sum_m t1(m,a) sum_e t1(i,e) <mb||ej>
    !
    ! Yring factorises into two thin DGEMMs rather than one (ov)^3, because the
    ! t1 t1 product separates: contract e first, then m.
    !$omp parallel do collapse(2) default(shared) private(i,j,a,b,e,m) schedule(static)
    do e = 1, nv
      do m = 1, no
        do a = 1, nv
          do i = 1, no
            t2ring((a-1)*no+i, (e-1)*no+m) = t2(i,m,a,e)
          end do
        end do
        do b = 1, nv
          do j = 1, no
            wring((e-1)*no+m, (b-1)*no+j) = Wmbej(m,b,e,j)
            gring(e, m, b, j) = g(m, no+b, no+e, j)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call dgemm('n','n', no*nv, no*nv, no*nv, 1.0_dp, t2ring, no*nv, &
               wring, no*nv, 0.0_dp, zring, no*nv)

    ! hring(i, m,b,j) = sum_e t1(i,e) <mb||ej>
    call dgemm('n','n', no, no*nv*no, nv, 1.0_dp, t1, no, &
               gring, nv, 0.0_dp, hring, no)
    !$omp parallel do collapse(3) default(shared) private(i,j,b,m) schedule(static)
    do j = 1, no
      do b = 1, nv
        do m = 1, no
          do i = 1, no
            hringp(m, i, b, j) = hring(i, m, b, j)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    ! yring(a, i,b,j) = sum_m t1(m,a) hringp(m, i,b,j)
    call dgemm('t','n', nv, no*nv*no, no, 1.0_dp, t1, no, &
               hringp, no, 0.0_dp, yring, nv)

    !$omp parallel do collapse(4) default(shared) private(i,j,a,b,e,m,n,f_,s_,d) schedule(static) if(no2*nv2 > PAR_MIN)
    do b = 1, nv; do a = 1, nv; do j = 1, no; do i = 1, no
      s_ = g(i,j,no+a,no+b)

      ! P(ab) sum_e t2_ijae FaeT(b,e),  FaeT = Fae - 1/2 t1^T Fme
      do e = 1, nv
        s_ = s_ + t2(i,j,a,e)*FaeT(b,e) - t2(i,j,b,e)*FaeT(a,e)
      end do

      ! -P(ij) sum_m t2_imab FmiT(m,j),  FmiT = Fmi + 1/2 Fme t1^T
      do m = 1, no
        s_ = s_ - t2(i,m,a,b)*FmiT(m,j) + t2(j,m,a,b)*FmiT(m,i)
      end do

      do n = 1, no; do m = 1, no
        s_ = s_ + 0.5_dp*tau(m,n,a,b)*Wmnij(m,n,i,j)
      end do; end do

      ! P(ij)P(ab) [ t2_imae Wmbej - t1_ie t1_ma <mb||ej> ], from the panels
      ! built above.  Zring is indexed (ia,jb), Yring as (a,i,b,j).
      s_ = s_ + zring((a-1)*no+i, (b-1)*no+j) - zring((a-1)*no+j, (b-1)*no+i) &
              - zring((b-1)*no+i, (a-1)*no+j) + zring((b-1)*no+j, (a-1)*no+i)
      s_ = s_ - yring(a,i,b,j) + yring(a,j,b,i) &
              + yring(b,i,a,j) - yring(b,j,a,i)

      do e = 1, nv
        s_ = s_ + t1(i,e)*g(no+a,no+b,no+e,j) - t1(j,e)*g(no+a,no+b,no+e,i)
      end do
      do m = 1, no
        s_ = s_ - t1(m,a)*g(m,no+b,i,j) + t1(m,b)*g(m,no+a,i,j)
      end do

      t2n(i,j,a,b) = s_
    end do; end do; end do; end do
    !$omp end parallel do

    ! Ladder: t2n(ij,ab) += 1/2 tau(ij,ef) Wabef(ab,ef).  This is the
    ! O(o^2 v^4) bottleneck and the one term worth a DGEMM above all others.
    call dgemm('n','t', no2, nv2, nv2, 0.5_dp, tau, no2, Wabef, nv2, &
               1.0_dp, t2n, no2)

    do b = 1, nv; do a = 1, nv; do j = 1, no; do i = 1, no
      dijab = eso(i) + eso(j) - eso(no+a) - eso(no+b)
      t2n(i,j,a,b) = t2n(i,j,a,b) / dijab
    end do; end do; end do; end do

    ! RMS, not the Euclidean norm: cc_lib divides by the amplitude count and
    ! [cc] conv is documented against that meaning.  Without the division the
    ! effective per-amplitude tolerance tightens as sqrt(namp) with system
    ! size, so a larger open-shell job could exhaust maxit on a threshold the
    ! closed-shell path would have called converged.
    rms = sqrt((sum((t1n-t1)**2) + sum((t2n-t2)**2)) / real(namp, dp))

    ! Push (amplitude, residual) onto the DIIS subspace, oldest evicted.
    if (ndiis > 0) then
      if (nvec < ndiis) then
        nvec = nvec + 1
        pos = nvec
      else
        dv(:,1:ndiis-1) = dv(:,2:ndiis)
        de(:,1:ndiis-1) = de(:,2:ndiis)
        pos = ndiis
      end if
      dv(1:no*nv, pos) = reshape(t1n, [no*nv])
      dv(no*nv+1:namp, pos) = reshape(t2n, [no*no*nv*nv])
      de(1:no*nv, pos) = reshape(t1n-t1, [no*nv])
      de(no*nv+1:namp, pos) = reshape(t2n-t2, [no*no*nv*nv])

      if (nvec > 1) then
        ndim = nvec + 1
        allocate(bmat(ndim,ndim), rhs(ndim), ipiv(ndim))
        bmat = 0.0_dp
        do ii = 1, nvec
          do jj = 1, nvec
            bmat(ii,jj) = dot_product(de(:,ii), de(:,jj))
          end do
        end do
        bmat(1:nvec, ndim) = -1.0_dp
        bmat(ndim, 1:nvec) = -1.0_dp
        rhs = 0.0_dp
        rhs(ndim) = -1.0_dp
        call dgesv(ndim, 1, bmat, ndim, ipiv, rhs, ndim, info)
        if (info == 0) then
          t1n = 0.0_dp; t2n = 0.0_dp
          do ii = 1, nvec
            t1n = t1n + rhs(ii)*reshape(dv(1:no*nv,ii), [no,nv])
            t2n = t2n + rhs(ii)*reshape(dv(no*nv+1:namp,ii), [no,no,nv,nv])
          end do
        end if
        deallocate(bmat, rhs, ipiv)
      end if
    end if

    t1 = t1n; t2 = t2n

    ! E = sum_ia f_ia t1_ia + 1/4 <ij||ab> t2 + 1/2 <ij||ab> t1 t1.
    ! The singles term is zero for a canonical reference -- Brillouin -- but
    ! not for ROHF, where it is the same size as the terms f_ov adds to the
    ! amplitude equations.
    e_ccsd = 0.0_dp
    do a = 1, nv; do i = 1, no
      e_ccsd = e_ccsd + f_ov(i,a)*t1(i,a)
    end do; end do
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

  call cc_uhf_wall_time(tw1)
  if (present(time_ccsd)) time_ccsd = tw1 - tw0

  ! --- (T) -----------------------------------------------------------------
  ! Same reasoning as cc_lib: a non-converged CCSD is going to abort, so do
  ! not spend the dominant cost of the job on triples built from it.
  if (do_triples .and. converged) then
    call cc_uhf_wall_time(tw0)
    call triples(e_t)
    call cc_uhf_wall_time(tw1)
    if (present(time_triples)) time_triples = tw1 - tw0
  end if

  deallocate(t1, t2, t1n, t2n, tau, taut, Fae, Fmi, Fme, Wmnij, Wabef, Wmbej, gmnef, f_ov)
  deallocate(t2ring, wring, zring, gring, hring, hringp, yring, FaeT, FmiT)
  deallocate(wjbnf, gnfme, wjbme)
  if (allocated(dv)) deallocate(dv, de)
  if (nblas_save > 0) call blas_thread_set(nblas_save)

contains

  !> Perturbative triples in the spin-orbital form:
  !>   E = sum_{i<j<k} sum_{a<b<c} t3c * D * (t3c + t3d)
  !> with t3c the connected and t3d the disconnected triple, both already
  !> divided by their denominator.
  !>
  !> Equivalently E = 1/36 sum over all i,j,k and all a,b,c: both triples are
  !> fully antisymmetric in each index set, so the 36 orderings are equal and
  !> repeated indices vanish.  Restricting the loops is the same sum evaluated
  !> once instead of 36 times.
  subroutine triples(et)
    real(dp), intent(out) :: et
    integer  :: i, j, k, a, b, c, e, m, jk, ik, ji
    integer  :: itr, ntrip
    integer, allocatable :: ti(:), tj(:), tk(:)
    real(dp) :: dd, tc, td, num, eijk, eab
    real(dp), allocatable :: Q(:,:,:,:)

    et = 0.0_dp

    allocate(gvv(nv, nv*nv, no), t2o(no, nv*nv, no), &
             t2v(nv, nv, no*no), gov(no, nv, no*no), gvvp(nv, nv, no*no))

    !$omp parallel do collapse(2) default(shared) private(i,b,c,e,m) schedule(static)
    do i = 1, no
      do c = 1, nv
        do b = 1, nv
          do e = 1, nv
            gvv(e, (c-1)*nv+b, i) = g(no+e, i, no+b, no+c)
          end do
          do m = 1, no
            t2o(m, (c-1)*nv+b, i) = t2(i, m, b, c)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do collapse(2) default(shared) private(j,k,a,e,m,b,c) schedule(static)
    do k = 1, no
      do j = 1, no
        do e = 1, nv
          do a = 1, nv
            t2v(a, e, (k-1)*no+j) = t2(j, k, a, e)
          end do
        end do
        do a = 1, nv
          do m = 1, no
            gov(m, a, (k-1)*no+j) = g(m, no+a, j, k)
          end do
        end do
        ! <jk||bc> as a contiguous (b,c) plane per occupied pair.  The
        ! disconnected numerator reads this six times per element; straight
        ! out of g those are strides of nso^2 and nso^3, which is what made
        ! the energy loop -- not the contractions -- the cost of (T) once the
        ! contractions became DGEMMs.
        do c = 1, nv
          do b = 1, nv
            gvvp(b, c, (k-1)*no+j) = g(j, k, no+b, no+c)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! Both t3c and t3d are fully antisymmetric in (i,j,k) and in (a,b,c), so of
    ! the 36 orderings of a given index set all carry the same contribution --
    ! which is exactly what the 1/36 below used to divide back out -- and any
    ! ordering with a repeated index vanishes.  Summing i<j<k and a<b<c once,
    ! with no 1/36, is therefore the same number for 1/36 of the arithmetic.
    ! It is an identity, not an approximation: the energies are unchanged.
    !
    ! The occupied triples are flattened into a list because an OpenMP
    ! collapse() needs rectangular bounds, and a triangular loop nest collapsed
    ! by hand is the same thing with worse arithmetic.
    ntrip = no*(no-1)*(no-2)/6
    allocate(ti(max(1,ntrip)), tj(max(1,ntrip)), tk(max(1,ntrip)))
    ntrip = 0
    do i = 1, no
      do j = i+1, no
        do k = j+1, no
          ntrip = ntrip + 1
          ti(ntrip) = i; tj(ntrip) = j; tk(ntrip) = k
        end do
      end do
    end do

    !$omp parallel default(shared) &
    !$omp   private(itr,i,j,k,a,b,c,e,m,dd,tc,td,num,Q,jk,ik,ji,eijk,eab) &
    !$omp   reduction(+:et)
    allocate(Q(nv,nv,nv,1))
    !$omp do schedule(dynamic)
    do itr = 1, ntrip
      i = ti(itr); j = tj(itr); k = tk(itr)
      ! The three occupied orderings enter the numerator with fixed signs,
      !   num = P(a,b,c) - P(b,a,c) - P(c,b,a),   P = Q(ijk) - Q(jik) - Q(kji),
      ! so accumulate them into ONE panel through the DGEMM scale factors
      ! rather than keeping three and combining per element.  Identical
      ! algebra, a third of the panel memory, and -- what matters at scale
      ! -- a third of the traffic: the panels are far past cache (27 MB per
      ! triple at cc-pVTZ on O2) so the energy loop streams them from DRAM,
      ! and it walks each one in three different index orders.
      call accum_q(i,j,k,  1.0_dp, 0.0_dp, Q(:,:,:,1))
      call accum_q(j,i,k, -1.0_dp, 1.0_dp, Q(:,:,:,1))
      call accum_q(k,j,i, -1.0_dp, 1.0_dp, Q(:,:,:,1))

      ! The three occupied pairs the disconnected term needs, matching
      ! disc = p1(i,j,k) - p1(j,i,k) - p1(k,j,i): each p1 reads the pair
      ! formed by its last two occupied arguments.
      jk = (k-1)*no + j
      ik = (k-1)*no + i
      ji = (i-1)*no + j
      eijk = eso(i) + eso(j) + eso(k)

      do a = 1, nv
        do b = a+1, nv
          eab = eijk - eso(no+a) - eso(no+b)
          do c = b+1, nv
            dd = eab - eso(no+c)
            if (abs(dd) < 1.0e-12_dp) cycle
            num = Q(a,b,c,1) - Q(b,a,c,1) - Q(c,b,a,1)
            tc = num / dd
            ! disc(i,j,k,a,b,c), inlined against the packed (b,c) planes.
            td = ( t1(i,a)*gvvp(b,c,jk) - t1(i,b)*gvvp(a,c,jk)      &
                 - t1(i,c)*gvvp(b,a,jk)                             &
                 + f_ov(i,a)*t2v(b,c,jk) - f_ov(i,b)*t2v(a,c,jk)    &
                 - f_ov(i,c)*t2v(b,a,jk) )                          &
               - ( t1(j,a)*gvvp(b,c,ik) - t1(j,b)*gvvp(a,c,ik)      &
                 - t1(j,c)*gvvp(b,a,ik)                             &
                 + f_ov(j,a)*t2v(b,c,ik) - f_ov(j,b)*t2v(a,c,ik)    &
                 - f_ov(j,c)*t2v(b,a,ik) )                          &
               - ( t1(k,a)*gvvp(b,c,ji) - t1(k,b)*gvvp(a,c,ji)      &
                 - t1(k,c)*gvvp(b,a,ji)                             &
                 + f_ov(k,a)*t2v(b,c,ji) - f_ov(k,b)*t2v(a,c,ji)    &
                 - f_ov(k,c)*t2v(b,a,ji) )
            td = td / dd
            et = et + tc*dd*(tc+td)
          end do
        end do
      end do
    end do
    !$omp end do
    deallocate(Q)
    !$omp end parallel
    deallocate(ti, tj, tk)

    deallocate(gvv, t2o, t2v, gov, gvvp)

  end subroutine triples

  !> Accumulate @p sgn times q2 for one occupied ordering into @p Q:
  !>   q2(a,b,c) = sum_e t2(j,k,a,e) <ei||bc> - sum_m t2(i,m,b,c) <ma||jk>
  !>
  !> Two DGEMMs on the panels triples() packed, in place of the loop nest this
  !> used to be.  Same flops; the strided reads of g and t2 now happen once per
  !> panel instead of once per element.  @p beta is the usual DGEMM meaning --
  !> 0 to start a fresh panel, 1 to add to the one already there -- which is
  !> how the three orderings fold into a single panel for free.
  subroutine accum_q(i, j, k, sgn, beta, Q)
    integer, intent(in) :: i, j, k
    real(dp), intent(in) :: sgn, beta
    real(dp), intent(inout) :: Q(nv,nv,nv)
    integer :: jk
    jk = (k-1)*no + j
    call dgemm('n','n', nv, nv*nv, nv, sgn, t2v(1,1,jk), nv, &
               gvv(1,1,i), nv, beta, Q, nv)
    call dgemm('t','n', nv, nv*nv, no, -sgn, gov(1,1,jk), no, &
               t2o(1,1,i), no, 1.0_dp, Q, nv)
  end subroutine accum_q


end subroutine cc_uhf_ccsd_t

!###############################################################################

!> Wall-clock seconds from the system clock.
subroutine cc_uhf_wall_time(t)
  real(dp), intent(out) :: t
  integer(8) :: c, r
  call system_clock(c, r)
  t = real(c, kind=dp)/real(r, kind=dp)
end subroutine cc_uhf_wall_time

end module cc_uhf_lib
