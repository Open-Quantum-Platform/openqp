!> @file qmrsf_ao2mo.F90
!> @brief QMRSF active-space AO->MO integral transformation (reuses int2).
!>
!> @details
!>   OpenQP is integral-direct: there is NO general ao2mo / four-index MO-integral
!>   transform anywhere in the code (only ao2mo_soc, a 1-index SOC transform).  The
!>   QMRSF backbone needs the one- and two-electron integrals over the four frontier
!>   active orbitals (the quintet SOMOs).  This module produces them WITHOUT ever
!>   forming or storing the AO 2e tensor, by REUSING the validated `int2` Coulomb
!>   digestor (the same engine that builds the SCF/MRSF Fock):
!>
!>     1e:  h_act(p,q) = sum_{mu,nu} C(mu,p) Hcore(mu,nu) C(nu,q).
!>     2e:  for each active pair (r,s) build the AO density  D^{rs} = C(:,r) C(:,s)^T,
!>          run it through int2 (pure Coulomb) to get  G^{rs}_{mu,nu} = (mu nu|rs),
!>          then transform the remaining pair:  (pq|rs) = C^T G^{rs} C.
!>          All nact*(nact+1)/2 pair densities are batched into ONE int2 sweep
!>          (the RHF consumer's `nfocks` channel), giving O(nact^2) screened builds
!>          and NO AO tensor in memory.
!>
!>   Frozen-core (inactive doubly-occupied orbitals, ncore = nelec_B for the quintet):
!>     the active one-electron integral is dressed by the core mean field
!>     v_core = 2 J[D_core] - K[D_core], and a constant E_core is returned.
!>     (For an all-active reference, e.g. H4/STO-3G, ncore=0 and these vanish, so
!>      E_core = E_nuc and h_eff = h_act -- the FCI gate case.)
!>
!>   CONVENTION NOTE (validated, not assumed): at the RAW int2_rhf consumer level
!>   (before fock_jk's 0.5x / diagonal-doubling post-scaling) each Fock channel gets
!>       f_raw(ij) += 4*scale_coulomb*val*d(kl)   -  scale_exchange*val*d(jl) ...
!>   uniformly for diagonal and off-diagonal ij.  Hence a pure Coulomb build with
!>   scale_coulomb=0.25, scale_exchange=0 and an off-diagonal-DOUBLED packed density
!>   returns J[D] as a clean (un-scaled) packed matrix.  The full eri_act is gated
!>   element-by-element against an independent closed-form oracle (route_a_oracle.py).
module qmrsf_ao2mo_mod
  use precision, only: dp
  implicit none
  private
  public :: qmrsf_active_integrals

contains

  !> @brief Build active-space h_eff, (pq|rs) and E_core from the converged
  !>        (quintet ROHF) reference, reusing int2.
  !> @param[inout] infos    OpenQP run container (provides basis, MOs, Hcore, enuc)
  !> @param[in]    nact     number of active orbitals (=4 for QMRSF)
  !> @param[in]    act      active MO indices (1-based, into the alpha MO set)
  !> @param[in]    ncore    number of inactive doubly-occupied (frozen-core) MOs
  !> @param[out]   h_act    nact x nact effective one-electron integrals (incl. core MF)
  !> @param[out]   eri_act  nact^4 two-electron integrals (pq|rs), CHEMIST order
  !> @param[out]   ecore    E_nuc + frozen-core electronic constant
  !> @param[out]   kcore_act optional active projection of K[D_core]
  subroutine qmrsf_active_integrals(infos, nact, act, ncore, h_act, eri_act, ecore, kcore_act)
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A, OQP_Hcore
    use int2_compute, only: int2_compute_t, int2_rhf_data_t
    use mathlib, only: unpack_matrix

    type(information), target, intent(inout) :: infos
    integer, intent(in)  :: nact, ncore
    integer, intent(in)  :: act(nact)
    real(dp), intent(out) :: h_act(nact,nact)
    real(dp), intent(out) :: eri_act(nact,nact,nact,nact)
    real(dp), intent(out) :: ecore
    real(dp), intent(out), optional :: kcore_act(nact,nact)

    type(basis_set), pointer :: basis
    real(dp), contiguous, pointer :: mo_a(:,:), hcore_p(:)
    real(dp), allocatable :: hcore_sq(:,:), cact(:,:), ccore(:,:)
    real(dp), allocatable :: tmp(:,:), gsq(:,:), gmo(:,:)
    real(dp), allocatable, target :: dprobe(:,:)
    real(dp), allocatable, target :: dcore(:,:)
    real(dp), allocatable :: vcore_sq(:,:), kcore_sq(:,:)
    integer, allocatable  :: pidx(:), qidx(:)
    type(int2_compute_t)  :: int2_driver
    type(int2_rhf_data_t) :: int2_data

    integer :: nbf, ntri, npair, i, p, q, r, s, pr, ii, jj, kl
    real(dp) :: psq

    basis => infos%basis
    if (present(kcore_act)) kcore_act = 0.0_dp
    nbf  = basis%nbf
    ntri = nbf*(nbf+1)/2

    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_Hcore,    hcore_p)

    allocate(hcore_sq(nbf,nbf))
    call unpack_matrix(hcore_p, hcore_sq, nbf, 'U')

    allocate(cact(nbf,nact))
    do i = 1, nact
      cact(:,i) = mo_a(:, act(i))
    end do

    ! ----------------------------- 1-electron -------------------------------
    allocate(tmp(nbf,nact))
    tmp   = matmul(hcore_sq, cact)
    h_act = matmul(transpose(cact), tmp)
    deallocate(tmp)

    ! ----------------------- 2-electron (probe Coulomb) ---------------------
    npair = nact*(nact+1)/2
    allocate(dprobe(ntri, npair), pidx(npair), qidx(npair))
    dprobe = 0.0_dp
    pr = 0
    do p = 1, nact
      do q = 1, p
        pr = pr + 1
        pidx(pr) = p; qidx(pr) = q
        do ii = 1, nbf
          do jj = 1, ii
            kl  = ii*(ii-1)/2 + jj          ! int2 packing: lower-tri row-major (=upper col-major)
            ! symmetric orbital-pair density P = 1/2 (C_p C_q^T + C_q C_p^T); NO off-diagonal
            ! doubling (matches OpenQP's pack_matrix/dtrttp convention).
            dprobe(kl, pr) = 0.5_dp*( cact(ii,p)*cact(jj,q) + cact(ii,q)*cact(jj,p) )
          end do
        end do
      end do
    end do

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    int2_data = int2_rhf_data_t(nfocks=npair, d=dprobe, &
                                scale_exchange=0.0_dp, scale_coulomb=1.0_dp)
    call int2_driver%run(int2_data, cam=.false.)

    allocate(gsq(nbf,nbf), gmo(nact,nact), tmp(nbf,nact))
    do pr = 1, npair
      r = pidx(pr); s = qidx(pr)
      ! raw consumer -> true J[P] via fock_jk's post-scaling (halve all, double diagonal)
      call fock_postscale(int2_data%f(:,pr,1), nbf, ntri, gsq)
      tmp = matmul(gsq, cact)
      gmo = matmul(transpose(cact), tmp)
      eri_act(:,:,r,s) = gmo
      eri_act(:,:,s,r) = gmo
    end do
    call int2_driver%clean()
    deallocate(gsq, gmo, tmp, dprobe, pidx, qidx)

    ! --------------------------- frozen core --------------------------------
    ecore = infos%mol_energy%nenergy   ! nuclear repulsion (SCF stores it here, not %enuc)
    if (ncore > 0) then
      allocate(ccore(nbf,ncore))
      do i = 1, ncore
        ccore(:,i) = mo_a(:, i)            ! inactive doubly-occupied MOs are 1..ncore
      end do
      ! D_core = sum_i C_i C_i^T (one-particle), packed, off-diagonal doubled
      allocate(dcore(ntri,1)); dcore = 0.0_dp
      do ii = 1, nbf
        do jj = 1, ii
          kl  = ii*(ii-1)/2 + jj          ! D_core = sum_i C_i C_i^T (undoubled)
          psq = 0.0_dp
          do i = 1, ncore
            psq = psq + ccore(ii,i)*ccore(jj,i)
          end do
          dcore(kl,1) = psq
        end do
      end do
      ! v_core = 2 J[D_core] - K[D_core]. After fock_postscale the raw int2_rhf consumer
      ! yields  scale_c*J - 0.5*scale_e*K  (calibrated vs the active-pair J build and the
      ! closed-shell fock_jk identity J-0.5K) -> set scale_c=2, scale_e=2.
      call int2_driver%init(basis, infos)
      call int2_driver%set_screening()
      int2_data = int2_rhf_data_t(nfocks=1, d=dcore, &
                                  scale_exchange=2.0_dp, scale_coulomb=2.0_dp)
      call int2_driver%run(int2_data, cam=.false.)
      allocate(vcore_sq(nbf,nbf))
      call fock_postscale(int2_data%f(:,1,1), nbf, ntri, vcore_sq)
      call int2_driver%clean()

      ! dress active one-electron integrals: h_eff += C_act^T v_core C_act
      allocate(tmp(nbf,nact))
      tmp   = matmul(vcore_sq, cact)
      h_act = h_act + matmul(transpose(cact), tmp)
      deallocate(tmp)

      ! E_core = E_nuc + 2 Tr[Dc.Hcore] + Tr[Dc.v_core] , with Dc = sum_i C_i C_i^T (square)
      block
        real(dp) :: dc_sq(nbf,nbf), e2h, evc
        integer  :: a, b
        dc_sq = matmul(ccore, transpose(ccore))
        e2h = 0.0_dp; evc = 0.0_dp
        do a = 1, nbf
          do b = 1, nbf
            e2h = e2h + dc_sq(a,b)*hcore_sq(a,b)
            evc = evc + dc_sq(a,b)*vcore_sq(a,b)
          end do
        end do
        ecore = ecore + 2.0_dp*e2h + evc
      end block

      ! The F^(0) paper builder needs the reference-functional core exchange
      ! separately: h_eff = h_act + (1-c_ref) K_core.  At the post-scaled
      ! int2_rhf consumer level, scale_exchange=-2 and scale_coulomb=0 gives
      ! +K[D_core].  Compute it only for callers that request the optional block.
      if (present(kcore_act)) then
        call int2_driver%init(basis, infos)
        call int2_driver%set_screening()
        int2_data = int2_rhf_data_t(nfocks=1, d=dcore, &
                                    scale_exchange=-2.0_dp, scale_coulomb=0.0_dp)
        call int2_driver%run(int2_data, cam=.false.)
        allocate(kcore_sq(nbf,nbf),tmp(nbf,nact))
        call fock_postscale(int2_data%f(:,1,1), nbf, ntri, kcore_sq)
        call int2_driver%clean()
        tmp = matmul(kcore_sq,cact)
        kcore_act = matmul(transpose(cact),tmp)
        deallocate(kcore_sq,tmp)
      end if
      deallocate(ccore, dcore, vcore_sq)
    end if

    deallocate(hcore_sq, cact)
  end subroutine qmrsf_active_integrals

  !> @brief Convert a RAW int2_rhf consumer Fock channel (packed) into the true
  !>        matrix, then unpack to a full square.  Mirrors fock_jk's post-scaling:
  !>        halve every element, then double the diagonal.  (At the raw consumer
  !>        level off-diagonal contributions are accumulated at 2x; the diagonal
  !>        at 1x after the global 0.5 -- hence the diagonal is restored by x2.)
  subroutine fock_postscale(fraw, nbf, ntri, fsq)
    use mathlib, only: unpack_matrix
    integer,  intent(in)  :: nbf, ntri
    real(dp), intent(in)  :: fraw(ntri)
    real(dp), intent(out) :: fsq(nbf,nbf)
    real(dp) :: fpk(ntri)
    integer  :: id
    fpk = 0.5_dp * fraw
    do id = 1, nbf
      fpk(id*(id+1)/2) = 2.0_dp * fpk(id*(id+1)/2)   ! diagonal packed index = id*(id+1)/2
    end do
    call unpack_matrix(fpk, fsq, nbf, 'U')
  end subroutine fock_postscale

end module qmrsf_ao2mo_mod
