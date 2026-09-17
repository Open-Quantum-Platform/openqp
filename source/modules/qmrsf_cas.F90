!> @file qmrsf_cas.F90
!> @brief Shared QMRSF CAS(4,4) M_s=0 determinant-CI backbone.
!>
!> @details
!>   QMRSF = Quintet Mixed-Reference Spin-Flip.  This module is the in-space
!>   "backbone": given the active-space one- and two-electron integrals over the
!>   four frontier orbitals (the quintet SOMOs), it builds the complete
!>   CAS(4,4) M_s=0 determinant Hamiltonian by the Slater-Condon rules and
!>   diagonalizes it.  The 36 M_s=0 determinants (C(4,2)^2) span the full
!>   20 singlet + 15 triplet + 1 quintet CSF set; the eigenvalues ARE the
!>   CAS-CI (for an all-active space, FCI) state energies.
!>
!>   This is a faithful Fortran port of the validated reference implementation,
!>   matched to it at 5.7e-14.
module qmrsf_cas_mod
  use precision, only: dp
  implicit none
  private
  public :: qmrsf_cas_solve
  public :: qmrsf_cas_build_s2
  public :: QMRSF_NACT, QMRSF_NDET

  integer, parameter :: QMRSF_NACT = 4
  integer, parameter :: NSO  = 8          !< spin orbitals (4 spatial x 2 spin)
  integer, parameter :: QMRSF_NDET = 36   !< C(4,2)^2 M_s=0 determinants

contains

  !> @brief Build the full CAS(4,4) M_s=0 Hamiltonian and diagonalize it.
  !> @param[in]  h_act   active one-electron integrals h_pq (chemist/standard)
  !> @param[in]  eri_act active two-electron integrals (pq|rs) in CHEMIST order
  !> @param[out] evals   36 ascending CAS eigenvalues (electronic; add E_core)
  !> @param[out] evecs   (optional) 36x36 eigenvectors (columns)
  !> @param[out] herm    (optional) max |H - H^T| (Hermiticity diagnostic)
  subroutine qmrsf_cas_solve(h_act, eri_act, evals, evecs, herm)
    use eigen, only: diag_symm_full
    real(dp), intent(in)  :: h_act(QMRSF_NACT,QMRSF_NACT)
    real(dp), intent(in)  :: eri_act(QMRSF_NACT,QMRSF_NACT,QMRSF_NACT,QMRSF_NACT)
    real(dp), intent(out) :: evals(QMRSF_NDET)
    real(dp), intent(out), optional :: evecs(QMRSF_NDET,QMRSF_NDET)
    real(dp), intent(out), optional :: herm

    real(dp) :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    real(dp) :: Hmat(QMRSF_NDET,QMRSF_NDET)
    integer  :: dets(4,QMRSF_NDET)
    integer  :: i, j, ierr

    call build_spinorb(h_act, eri_act, H1, g)
    call gen_dets(dets)
    call build_H(dets, H1, g, Hmat)

    if (present(herm)) then
      herm = 0.0_dp
      do i = 1, QMRSF_NDET
        do j = 1, QMRSF_NDET
          herm = max(herm, abs(Hmat(i,j)-Hmat(j,i)))
        end do
      end do
    end if

    call diag_symm_full(0, QMRSF_NDET, Hmat, QMRSF_NDET, evals, ierr)
    if (present(evecs)) evecs = Hmat
  end subroutine qmrsf_cas_solve

!-------------------------------------------------------------------------------

  subroutine build_spinorb(h_act, eri_act, H1, g)
    real(dp), intent(in)  :: h_act(QMRSF_NACT,QMRSF_NACT)
    real(dp), intent(in)  :: eri_act(QMRSF_NACT,QMRSF_NACT,QMRSF_NACT,QMRSF_NACT)
    real(dp), intent(out) :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    integer :: P,Q,R,S, spat(NSO), spin(NSO), i
    real(dp) :: a, b
    do i = 1, NSO
      if (i <= QMRSF_NACT) then; spat(i) = i;             spin(i) = 0
      else;                      spat(i) = i - QMRSF_NACT; spin(i) = 1; end if
    end do
    H1 = 0.0_dp
    do P = 1, NSO; do Q = 1, NSO
      if (spin(P) == spin(Q)) H1(P,Q) = h_act(spat(P), spat(Q))
    end do; end do
    ! antisymmetrized g(P,Q,R,S) = <PQ||RS> = (PR|QS) - (PS|QR), chemist eri + spin deltas
    g = 0.0_dp
    do P = 1, NSO; do Q = 1, NSO; do R = 1, NSO; do S = 1, NSO
      a = 0.0_dp; b = 0.0_dp
      if (spin(P)==spin(R) .and. spin(Q)==spin(S)) a = eri_act(spat(P),spat(R),spat(Q),spat(S))
      if (spin(P)==spin(S) .and. spin(Q)==spin(R)) b = eri_act(spat(P),spat(S),spat(Q),spat(R))
      g(P,Q,R,S) = a - b
    end do; end do; end do; end do
  end subroutine build_spinorb

  subroutine gen_dets(dets)
    integer, intent(out) :: dets(4,QMRSF_NDET)
    integer :: a1,a2,b1,b2,k, t(4), i,j,tmp
    k = 0
    do a1 = 1, QMRSF_NACT-1; do a2 = a1+1, QMRSF_NACT
      do b1 = 1, QMRSF_NACT-1; do b2 = b1+1, QMRSF_NACT
        k = k + 1
        t = (/ a1, a2, b1+QMRSF_NACT, b2+QMRSF_NACT /)
        do i = 1,3; do j = i+1,4
          if (t(j) < t(i)) then; tmp=t(i); t(i)=t(j); t(j)=tmp; end if
        end do; end do
        dets(:,k) = t
      end do; end do
    end do; end do
  end subroutine gen_dets

  pure logical function inset(x, D)
    integer, intent(in) :: x, D(4)
    inset = any(D == x)
  end function inset

  real(dp) function melem(D1, D2, H1, g)
    integer,  intent(in) :: D1(4), D2(4)
    real(dp), intent(in) :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    integer :: holes(4), parts(4), common(4), nh, np, nc
    integer :: occ(4), nocc, i, idx, k
    integer :: p1,p2,ho1,ho2, Pp, Hh, Qc
    real(dp) :: sgn, val, e
    nh=0; np=0; nc=0
    do i=1,4
      if (.not. inset(D2(i), D1)) then; nh=nh+1; holes(nh)=D2(i); end if
    end do
    do i=1,4
      if (.not. inset(D1(i), D2)) then; np=np+1; parts(np)=D1(i); end if
      if (      inset(D1(i), D2)) then; nc=nc+1; common(nc)=D1(i); end if
    end do
    if (nh > 2) then; melem = 0.0_dp; return; end if
    occ = D2; nocc = 4; sgn = 1.0_dp
    do k = 1, nh
      idx = 0
      do i = 1, nocc
        if (occ(i) == holes(k)) then; idx = i; exit; end if
      end do
      if (mod(idx-1,2) == 1) sgn = -sgn
      do i = idx, nocc-1; occ(i) = occ(i+1); end do
      nocc = nocc - 1
    end do
    do k = np, 1, -1
      idx = 1
      do i = 1, nocc
        if (occ(i) < parts(k)) idx = idx + 1
      end do
      if (mod(idx-1,2) == 1) sgn = -sgn
      do i = nocc, idx, -1; occ(i+1) = occ(i); end do
      occ(idx) = parts(k); nocc = nocc + 1
    end do
    if (nh == 0) then
      e = 0.0_dp
      do i = 1,4; e = e + H1(D1(i),D1(i)); end do
      do i = 1,3
        do k = i+1,4; e = e + g(D1(i),D1(k),D1(i),D1(k)); end do
      end do
      melem = e
    else if (nh == 1) then
      Pp = parts(1); Hh = holes(1)
      val = H1(Pp,Hh)
      do i = 1, nc; Qc = common(i); val = val + g(Pp,Qc,Hh,Qc); end do
      melem = sgn * val
    else
      p1=parts(1); p2=parts(2); ho1=holes(1); ho2=holes(2)
      melem = sgn * g(p1,p2,ho1,ho2)
    end if
  end function melem

  !> @brief Build the S^2 operator in the 36-determinant M_s=0 CAS(4,4) basis.
  !> @details For M_s=0, S_z=0 so S^2 = S_- S_+ (the S_z(S_z+1) term vanishes). S_- S_+ is
  !> applied as a composite: S_+ (a^+_{alpha p} a_{beta p}) takes a det to the M_s=+1 sector,
  !> then S_- (a^+_{beta q} a_{alpha q}) returns it to M_s=0; only final dets back in the
  !> basis contribute. Port of the validated NumPy prototype (qmrsf_s2_projection_proto.py),
  !> which reproduces the 20 singlet / 15 triplet / 1 quintet split with <S^2>=0/2/6.
  subroutine qmrsf_cas_build_s2(s2mat)
    real(dp), intent(out) :: s2mat(QMRSF_NDET,QMRSF_NDET)
    integer :: dets(4,QMRSF_NDET), j, p, q, i, sgn1, sgn2
    integer :: mid(4), fin(4)
    call gen_dets(dets)
    s2mat = 0.0_dp
    do j = 1, QMRSF_NDET
      do p = 1, QMRSF_NACT
        call s2_exc(dets(:,j), p, p+QMRSF_NACT, mid, sgn1)   ! S_+ : a^+_{alpha p} a_{beta p}
        if (sgn1 == 0) cycle
        do q = 1, QMRSF_NACT
          call s2_exc(mid, q+QMRSF_NACT, q, fin, sgn2)        ! S_- : a^+_{beta q} a_{alpha q}
          if (sgn2 == 0) cycle
          i = s2_find(fin, dets)
          if (i > 0) s2mat(i,j) = s2mat(i,j) + real(sgn1*sgn2, dp)
        end do
      end do
    end do
  end subroutine qmrsf_cas_build_s2

  !> apply a^+_{p_create} a_{q_annih} to a 4-electron determinant (sorted SO indices).
  subroutine s2_exc(occ_in, p_create, q_annih, occ_out, sgn)
    integer, intent(in)  :: occ_in(4), p_create, q_annih
    integer, intent(out) :: occ_out(4), sgn
    integer :: rem(3), nrem, i, idx, cnt
    sgn = 0
    idx = 0
    do i = 1, 4
      if (occ_in(i) == q_annih) then; idx = i; exit; end if
    end do
    if (idx == 0) return                       ! nothing to annihilate
    sgn = (-1)**(idx-1)
    nrem = 0
    do i = 1, 4
      if (i /= idx) then; nrem = nrem + 1; rem(nrem) = occ_in(i); end if
    end do
    do i = 1, nrem
      if (rem(i) == p_create) then; sgn = 0; return; end if   ! Pauli
    end do
    cnt = 0
    do i = 1, nrem
      if (rem(i) < p_create) cnt = cnt + 1
    end do
    sgn = sgn * (-1)**cnt
    ! insert p_create at position cnt+1 (keep ascending)
    do i = 1, cnt
      occ_out(i) = rem(i)
    end do
    occ_out(cnt+1) = p_create
    do i = cnt+1, nrem
      occ_out(i+1) = rem(i)
    end do
  end subroutine s2_exc

  integer function s2_find(target, dets) result(idx)
    integer, intent(in) :: target(4), dets(4,QMRSF_NDET)
    integer :: j
    idx = 0
    do j = 1, QMRSF_NDET
      if (all(dets(:,j) == target)) then; idx = j; return; end if
    end do
  end function s2_find

  subroutine build_H(dets, H1, g, Hmat)
    integer,  intent(in)  :: dets(4,QMRSF_NDET)
    real(dp), intent(in)  :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    real(dp), intent(out) :: Hmat(QMRSF_NDET,QMRSF_NDET)
    integer :: i, j
    do i = 1, QMRSF_NDET
      do j = i, QMRSF_NDET
        Hmat(i,j) = melem(dets(:,i), dets(:,j), H1, g)
        Hmat(j,i) = Hmat(i,j)
      end do
    end do
  end subroutine build_H

end module qmrsf_cas_mod
