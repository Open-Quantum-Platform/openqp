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
  !$omp   private(p,q,r,s,ps,qs,rs_,ss,pp,qq,rr,ssp,coul,exch) schedule(static) &
  !$omp   if(real(2*nmo,dp)**4 > PAR_WORK)
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
!> sixteen times the spatial ones.  Use it for small open-shell systems.
subroutine cc_uhf_ccsd_t(nso, nocc, eso, g, maxit, conv, do_triples, &
                         e_ccsd, e_t, converged, niter, fov)

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

  real(dp), allocatable :: t1(:,:), t2(:,:,:,:), t1n(:,:), t2n(:,:,:,:)
  real(dp), allocatable :: tau(:,:,:,:), taut(:,:,:,:)
  real(dp), allocatable :: Fae(:,:), Fmi(:,:), Fme(:,:)
  real(dp), allocatable :: Wmnij(:,:,:,:), Wabef(:,:,:,:), Wmbej(:,:,:,:)
  real(dp), allocatable :: gmnef(:,:,:,:), f_ov(:,:)
  integer(c_int64_t) :: nblas_save
  ! DIIS history over the concatenated (t1,t2) vector and its residual.
  real(dp), allocatable :: dv(:,:), de(:,:), bmat(:,:), rhs(:)
  integer, allocatable :: ipiv(:)
  integer :: no2, nv2, namp, ndiis, ndim, nvec, pos, ii, jj, info
  integer  :: no, nv, i, j, k, a, b, c, m, n, e, f_, it
  real(dp) :: s_, d, eold, rms, dia, dijab

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
  if (no <= 0 .or. nv <= 0) return

  allocate(t1(no,nv), t2(no,no,nv,nv), t1n(no,nv), t2n(no,no,nv,nv))
  allocate(tau(no,no,nv,nv), taut(no,no,nv,nv))
  allocate(Fae(nv,nv), Fmi(no,no), Fme(no,nv))
  allocate(Wmnij(no,no,no,no), Wabef(nv,nv,nv,nv), Wmbej(no,nv,nv,no))
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
  ! same subspace size cc_lib uses.
  namp = no*nv + no*no*nv*nv
  ndiis = 8
  allocate(dv(namp,ndiis), de(namp,ndiis))
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
      do f_ = 1, nv; do n = 1, no
        s_ = s_ - (0.5_dp*t2(j,n,f_,b) + t1(j,f_)*t1(n,b))*g(m,n,no+e,no+f_)
      end do; end do
      Wmbej(m,b,e,j) = s_
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
    !$omp parallel do collapse(4) default(shared) private(i,j,a,b,e,m,n,f_,s_,d) schedule(static) if(no2*nv2 > PAR_MIN)
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
    !$omp end parallel do

    ! Ladder: t2n(ij,ab) += 1/2 tau(ij,ef) Wabef(ab,ef).  This is the
    ! O(o^2 v^4) bottleneck and the one term worth a DGEMM above all others.
    call dgemm('n','t', no2, nv2, nv2, 0.5_dp, tau, no2, Wabef, nv2, &
               1.0_dp, t2n, no2)

    do b = 1, nv; do a = 1, nv; do j = 1, no; do i = 1, no
      dijab = eso(i) + eso(j) - eso(no+a) - eso(no+b)
      t2n(i,j,a,b) = t2n(i,j,a,b) / dijab
    end do; end do; end do; end do

    rms = sqrt(sum((t1n-t1)**2) + sum((t2n-t2)**2))

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

  ! --- (T) -----------------------------------------------------------------
  if (do_triples) call triples(e_t)

  deallocate(t1, t2, t1n, t2n, tau, taut, Fae, Fmi, Fme, Wmnij, Wabef, Wmbej, gmnef, f_ov)
  if (allocated(dv)) deallocate(dv, de)
  if (nblas_save > 0) call blas_thread_set(nblas_save)

contains

  !> Perturbative triples in the spin-orbital form:
  !>   E = 1/36 sum t3c * D * (t3c + t3d)
  !> with t3c the connected and t3d the disconnected triple, both already
  !> divided by their denominator.
  subroutine triples(et)
    real(dp), intent(out) :: et
    integer  :: i, j, k, a, b, c, e, m
    real(dp) :: dd, tc, td, num
    real(dp), allocatable :: Q(:,:,:,:)

    et = 0.0_dp
    !$omp parallel default(shared) private(i,j,k,a,b,c,e,m,dd,tc,td,num,Q) &
    !$omp   reduction(+:et)
    allocate(Q(nv,nv,nv,3))
    !$omp do collapse(3) schedule(dynamic)
    do i = 1, no
      do j = 1, no
        do k = 1, no
          ! q2 depends on the occupied triple only through three orderings,
          ! and the connected numerator asks for each of them three times over
          ! permuted virtuals.  Building the three nv^3 panels once per (i,j,k)
          ! replaces nine recomputations per element -- each of which had its
          ! own loop over e and m -- with nine array reads.
          call build_q(i,j,k, Q(:,:,:,1))
          call build_q(j,i,k, Q(:,:,:,2))
          call build_q(k,j,i, Q(:,:,:,3))

          do a = 1, nv
            do b = 1, nv
              do c = 1, nv
                dd = eso(i)+eso(j)+eso(k)-eso(no+a)-eso(no+b)-eso(no+c)
                if (abs(dd) < 1.0e-12_dp) cycle
                num =  (Q(a,b,c,1) - Q(b,a,c,1) - Q(c,b,a,1)) &
                     - (Q(a,b,c,2) - Q(b,a,c,2) - Q(c,b,a,2)) &
                     - (Q(a,b,c,3) - Q(b,a,c,3) - Q(c,b,a,3))
                tc = num / dd
                td = disc(i,j,k,a,b,c) / dd
                et = et + tc*dd*(tc+td) / 36.0_dp
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end do
    deallocate(Q)
    !$omp end parallel
  end subroutine triples

  !> q2 over the whole virtual cube for one occupied ordering:
  !>   Q(a,b,c) = sum_e t2(j,k,a,e) <ei||bc> - sum_m t2(i,m,b,c) <ma||jk>
  subroutine build_q(i, j, k, Q)
    integer, intent(in) :: i, j, k
    real(dp), intent(out) :: Q(nv,nv,nv)
    integer :: a, b, c, e, m
    real(dp) :: s1, s2
    do c = 1, nv
      do b = 1, nv
        do a = 1, nv
          s1 = 0.0_dp
          do e = 1, nv
            s1 = s1 + t2(j,k,a,e)*g(no+e,i,no+b,no+c)
          end do
          s2 = 0.0_dp
          do m = 1, no
            s2 = s2 + t2(i,m,b,c)*g(m,no+a,j,k)
          end do
          Q(a,b,c) = s1 - s2
        end do
      end do
    end do
  end subroutine build_q

  !> Disconnected triple numerator: P(i/jk) P(a/bc) [ t1_ia <jk||bc>
  !> + f_ia t2_jkbc ], the second term non-zero only for a non-HF reference.
  pure real(dp) function disc(i,j,k,a,b,c) result(v)
    integer, intent(in) :: i,j,k,a,b,c
    v =   p1(i,j,k,a,b,c) - p1(j,i,k,a,b,c) - p1(k,j,i,a,b,c)
  end function disc

  pure real(dp) function p1(i,j,k,a,b,c) result(v)
    integer, intent(in) :: i,j,k,a,b,c
    v =   t1(i,a)*g(j,k,no+b,no+c) - t1(i,b)*g(j,k,no+a,no+c) &
        - t1(i,c)*g(j,k,no+b,no+a) &
        + f_ov(i,a)*t2(j,k,b,c) - f_ov(i,b)*t2(j,k,a,c) - f_ov(i,c)*t2(j,k,b,a)
  end function p1

end subroutine cc_uhf_ccsd_t

end module cc_uhf_lib
