module nmr_shielding_mod

  implicit none

  character(len=*), parameter :: module_name = "nmr_shielding_mod"

  private
  public nmr_shielding

contains

  subroutine nmr_shielding_C(c_handle) bind(C, name="nmr_shielding")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call nmr_shielding(inf)
  end subroutine nmr_shielding_C

!> @brief NMR nuclear magnetic shielding tensors (common gauge origin).
!> @details v1 scope: RHF and closed-shell DFT, common gauge origin (CGO). Both
!>  the uncoupled and the coupled (CPHF/CPKS) paramagnetic responses are reported.
!>  The coupled response includes the exact-exchange response of the imaginary
!>  antisymmetric first-order magnetic density (Phase 0); the Coulomb and
!>  semi-local XC-kernel responses vanish by symmetry, so the coupling is exact
!>  exchange scaled by the exchange fraction c_x (zero for pure functionals, where
!>  coupled == uncoupled). The isotropic shielding is sigma = sigma_dia +
!>  sigma_para per nucleus.
!>
!> Phase-0 validation (H2O/STO-3G, CGO at COM; PySCF common-gauge oracle in
!> tests/fixtures/nmr/pyscf_cgo_reference.json):
!>   - HF coupled para matches the oracle exactly (O -230.63, H 3.506 ppm).
!>   - Coulomb response of P^B ~0; exact-exchange response nonzero (gates 1-2).
!>   - Pure PBE: coupled == uncoupled (gate 3). HF/hybrid coupled != uncoupled,
!>     with the coupling scaling with c_x (gates 4, 6).
!>
!> Validation (H2O/STO-3G, RHF, CGO at the center of mass; PySCF reference):
!>   - Diamagnetic term matches PySCF common-gauge `dia()` to ~6 significant
!>     figures (O 411.418, H 28.062 ppm).
!>   - Paramagnetic term (MO transform of the orbital-Zeeman and PSO operators,
!>     occupied-virtual sum-over-states, 2*alpha^2 prefactor) matches the PySCF
!>     uncoupled reference for BOTH atoms (O para -113.63, H para 1.785 ppm;
!>     totals 297.79 / 29.85 ppm).
!>   - The PSO operator is anti-Hermitian; `pso_integrals` returns it exactly
!>     antisymmetric (max|diag| and max|A+A^T| are reported below as a check).
  subroutine nmr_shielding(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use int1, only: angular_momentum_integrals, nmr_dia_shielding, pso_integrals

    implicit none

    character(len=*), parameter :: subroutine_name = "nmr_shielding"
    ! CODATA fine-structure constant and derived prefactors
    real(kind=8), parameter :: ALPHA = 1.0d0/137.035999084d0
    real(kind=8), parameter :: halfa2 = 0.5d0*ALPHA*ALPHA
    real(kind=8), parameter :: twoa2 = 2.0d0*ALPHA*ALPHA
    real(kind=8), parameter :: PPM = 1.0d6

    type(information), target, intent(inout) :: infos

    integer :: nbf, nbf2, ok
    logical :: urohf
    type(basis_set), pointer :: basis
    real(kind=8), allocatable :: amom(:,:)      ! packed Lx,Ly,Lz (lower triangle)
    real(kind=8), allocatable :: lfull(:,:,:)   ! full antisymmetric (nbf,nbf,3)
    real(kind=8), allocatable :: gdia(:,:,:)    ! diamagnetic contracted integrals (3,3,nat)
    real(kind=8), allocatable :: sig_dia(:,:,:) ! diamagnetic shielding tensor (3,3,nat)
    real(kind=8), allocatable :: coords(:,:)    ! nuclear coordinates (3,nat)
    real(kind=8), allocatable :: siso_dia(:)    ! isotropic diamagnetic shielding (ppm)
    real(kind=8), allocatable :: pso_full(:,:,:)! full antisymmetric PSO (nbf,nbf,3)
    real(kind=8), allocatable :: moL(:,:,:)     ! orbital Zeeman in MO basis (nmo,nmo,3)
    real(kind=8), allocatable :: moP(:,:,:)     ! PSO in MO basis (nmo,nmo,3)
    real(kind=8), allocatable :: sig_para(:,:,:)! paramagnetic shielding tensor (3,3,nat)
    real(kind=8), allocatable :: siso_para(:), siso_tot(:)
    ! Phase 0: coupled (CPHF/CPKS) magnetic response
    real(kind=8), allocatable :: rcoup(:,:)     ! coupled response vectors (lvir,3)
    real(kind=8), allocatable :: sig_para_c(:,:,:), siso_para_c(:), siso_tot_c(:)
    real(kind=8) :: scale_exch, pb_asym, jnorm, knorm
    integer :: nvir, lvir
    logical :: is_dft
    real(kind=8) :: o(3), com(3), trg, pso_diag_max, pso_asym_max
    integer :: nat, i, c, t, s, nocc, nmo

    real(kind=8), contiguous, pointer :: dmat_a(:)
    real(kind=8), contiguous, pointer :: mo_a(:,:)
    real(kind=8), contiguous, pointer :: mo_e_a(:)
    integer(4) :: status

    urohf = infos%control%scftype == 2 .or. infos%control%scftype == 3

    open (unit=IW, file=infos%log_filename, position="append")

    basis => infos%basis
    basis%atoms => infos%atoms

    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2

    write(iw,'(2/)')
    write(iw,'(4x,a)') '======================================'
    write(iw,'(4x,a)') 'NMR nuclear magnetic shielding (CGO)'
    write(iw,'(4x,a)') '======================================'
    call flush(iw)

    if (urohf) then
      call show_message('NMR shielding currently supports closed-shell &
        &(RHF / pure-DFT) references only', with_abort)
    end if

    ! Confirm the SCF density is present (used by later stages)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a, status)
    call check_status(status, module_name, subroutine_name, OQP_DM_A)

    ! Gauge origin: center of mass (default)
    nat = ubound(basis%atoms%zn,1)
    com = 0
    do i = 1, nat
      com = com + basis%atoms%xyz(:,i)*basis%atoms%mass(i)
    end do
    com = com / sum(basis%atoms%mass)
    o = com

    write(iw,'(/4x,a)') 'Gauge origin (Bohr):'
    write(iw,'(4x,a,3f15.8)') 'O = ', o

    ! Angular momentum integrals about the gauge origin (packed, antisymmetric)
    allocate(amom(nbf2,3), source=0.0d0, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    call angular_momentum_integrals(basis, amom, o)

    ! Expand each component to a full antisymmetric nbf x nbf matrix (the
    ! orbital-Zeeman / magnetic-field perturbation used by the paramagnetic term).
    allocate(lfull(nbf,nbf,3), source=0.0d0, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    do c = 1, 3
      call expand_antisym(amom(:,c), lfull(:,:,c), nbf)
    end do

    !-------------------------------------------------------------------------
    ! Diamagnetic term
    !   sigma^dia_{ts}(N) = (alpha^2/2) [ delta_ts*Tr(g^N) - g^N_{s,t} ]
    ! where g^N_{ab} = sum_{mu,nu} P_{mu,nu} <mu|(r-O)_a (r-R_N)_b/|r-R_N|^3|nu>
    !-------------------------------------------------------------------------
    allocate(gdia(3,3,nat), coords(3,nat), source=0.0d0, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)
    do i = 1, nat
      coords(:,i) = basis%atoms%xyz(:,i)
    end do

    call nmr_dia_shielding(basis, dmat_a, o, coords, nat, gdia)

    allocate(sig_dia(3,3,nat), siso_dia(nat), source=0.0d0)
    do i = 1, nat
      trg = gdia(1,1,i) + gdia(2,2,i) + gdia(3,3,i)
      do t = 1, 3
        do s = 1, 3
          sig_dia(t,s,i) = halfa2 * (merge(trg, 0.0d0, t==s) - gdia(s,t,i))
        end do
      end do
      siso_dia(i) = (sig_dia(1,1,i)+sig_dia(2,2,i)+sig_dia(3,3,i))/3.0d0 * PPM
    end do

    !-------------------------------------------------------------------------
    ! Paramagnetic term (uncoupled; exact CPKS for pure functionals/HF-uncoupled)
    !   sigma^para_{xy}(N) = 2 alpha^2 sum_{i occ, a vir}
    !                          MO_x(a,i) * PSO^N_y(i,a) / (eps_a - eps_i)
    ! MO_x = C^T A_O[x] C  (orbital Zeeman, angular momentum about O)
    ! PSO^N_y = C^T A_PSO^N[y] C
    !-------------------------------------------------------------------------
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a, status)
    call check_status(status, module_name, subroutine_name, OQP_VEC_MO_A)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_e_a, status)
    call check_status(status, module_name, subroutine_name, OQP_E_MO_A)

    nocc = int(infos%mol_prop%nocc)
    nmo  = size(mo_e_a)

    ! Orbital-Zeeman operator in MO basis (3 components)
    allocate(moL(nmo,nmo,3), source=0.0d0)
    do c = 1, 3
      call ao_to_mo(lfull(:,:,c), mo_a, moL(:,:,c), nbf, nmo)
    end do

    ! ----------------------------------------------------------------------
    ! Phase 0: coupled magnetic response (CPHF/CPKS).
    ! The first-order magnetic density is imaginary/antisymmetric, so the
    ! Coulomb and (semi-local) XC-kernel responses vanish; only the exact
    ! exchange response survives, scaled by the exchange fraction c_x. For
    ! pure functionals (c_x = 0) the coupled response equals the uncoupled one.
    ! ----------------------------------------------------------------------
    nvir = nmo - nocc
    lvir = nocc*nvir
    is_dft = infos%control%hamilton == 20
    scale_exch = 1.0d0
    if (is_dft) scale_exch = infos%dft%HFscale
    allocate(rcoup(lvir,3), source=0.0d0)
    call compute_coupled_para(infos, basis, mo_a, mo_e_a, nocc, nbf, nvir, &
                              moL, scale_exch, rcoup, pb_asym, jnorm, knorm)

    allocate(pso_full(nbf,nbf,3), moP(nmo,nmo,3), source=0.0d0)
    allocate(sig_para(3,3,nat), siso_para(nat), siso_tot(nat), source=0.0d0)
    allocate(sig_para_c(3,3,nat), siso_para_c(nat), siso_tot_c(nat), source=0.0d0)

    pso_diag_max = 0.0d0
    pso_asym_max = 0.0d0
    do i = 1, nat
      ! Full antisymmetric PSO matrices A_a = [(r-R_N) x grad]_a/|r-R_N|^3
      call pso_integrals(basis, coords(:,i), pso_full)
      do c = 1, 3
        ! Diagnostics: the PSO operator is anti-Hermitian, so the diagonal and
        ! the symmetric part must vanish.
        do t = 1, nbf
          pso_diag_max = max(pso_diag_max, abs(pso_full(t,t,c)))
          do s = 1, nbf
            pso_asym_max = max(pso_asym_max, abs(pso_full(t,s,c)+pso_full(s,t,c)))
          end do
        end do
        call ao_to_mo(pso_full(:,:,c), mo_a, moP(:,:,c), nbf, nmo)
      end do
      do t = 1, 3
        do s = 1, 3
          sig_para(t,s,i)   = twoa2 * sum_ov(moL(:,:,t), moP(:,:,s), mo_e_a, nocc, nmo)
          sig_para_c(t,s,i) = -twoa2 * sum_resp(rcoup(:,t), moP(:,:,s), nocc, nvir)
        end do
      end do
      siso_para(i)   = (sig_para(1,1,i)+sig_para(2,2,i)+sig_para(3,3,i))/3.0d0 * PPM
      siso_tot(i)    = siso_dia(i) + siso_para(i)
      siso_para_c(i) = (sig_para_c(1,1,i)+sig_para_c(2,2,i)+sig_para_c(3,3,i))/3.0d0 * PPM
      siso_tot_c(i)  = siso_dia(i) + siso_para_c(i)
    end do

    write(iw,'(/4x,a,f8.4,a)') 'Isotropic shielding (CGO, ppm)   [exact-exchange c_x =', &
           scale_exch, ']'
    write(iw,'(4x,a)')  '   Atom    Z   sigma_dia   para_uncoupled   para_coupled'// &
           '   total_uncoupled   total_coupled'
    do i = 1, nat
      write(iw,'(4x,i6,f6.1,5f16.6)') i, basis%atoms%zn(i), &
             siso_dia(i), siso_para(i), siso_para_c(i), siso_tot(i), siso_tot_c(i)
    end do

    ! ---- Phase-0 validation gates (reported as diagnostics) ----
    write(iw,'(/4x,a)') 'Phase-0 magnetic-response gates:'
    write(iw,'(4x,a,es12.3)') '  gate0  max|P^B + (P^B)^T|        = ', pb_asym
    write(iw,'(4x,a,es12.3)') '  gate1  ||J(P^B)|| (Coulomb)      = ', jnorm
    write(iw,'(4x,a,es12.3)') '  gate2  ||K(P^B)|| (exact exch.)  = ', knorm
    write(iw,'(4x,a,2es12.3)') '  PSO    max|diag|, max|A+A^T|    = ', &
           pso_diag_max, pso_asym_max
    call flush(iw)

    deallocate(amom, lfull, gdia, coords, sig_dia, siso_dia)
    deallocate(pso_full, moL, moP, sig_para, siso_para, siso_tot)
    deallocate(rcoup, sig_para_c, siso_para_c, siso_tot_c)
    close(iw)

  end subroutine nmr_shielding

!> @brief Expand a packed lower-triangular antisymmetric matrix to full form.
!> @details Packed storage holds the bra>=ket elements A(p,q) (p>=q) at index
!>  q + p*(p-1)/2. The full matrix satisfies A(q,p) = -A(p,q), zero diagonal.
  subroutine expand_antisym(packed, full, n)
    real(kind=8), intent(in)  :: packed(:)
    real(kind=8), intent(out) :: full(:,:)
    integer, intent(in) :: n
    integer :: p, q, idx
    full = 0.0d0
    do p = 1, n
      do q = 1, p
        idx = q + p*(p-1)/2
        full(p,q) =  packed(idx)
        full(q,p) = -packed(idx)
      end do
    end do
  end subroutine expand_antisym

!> @brief Transform a full AO matrix to the MO basis: M = C^T A C.
  subroutine ao_to_mo(a_ao, c, m_mo, nbf, nmo)
    real(kind=8), intent(in)  :: a_ao(:,:)   ! (nbf,nbf)
    real(kind=8), intent(in)  :: c(:,:)      ! (nbf,nmo)
    real(kind=8), intent(out) :: m_mo(:,:)   ! (nmo,nmo)
    integer, intent(in) :: nbf, nmo
    real(kind=8), allocatable :: tmp(:,:)
    allocate(tmp(nbf,nmo))
    tmp = matmul(a_ao, c(:,1:nmo))
    m_mo = matmul(transpose(c(:,1:nmo)), tmp)
    deallocate(tmp)
  end subroutine ao_to_mo

!> @brief Occupied-virtual sum-over-states contraction
!>   sum_{i occ, a vir} L(a,i) * P(i,a) / (eps_a - eps_i)
  function sum_ov(lmo, pmo, e, nocc, nmo) result(val)
    real(kind=8), intent(in) :: lmo(:,:), pmo(:,:), e(:)
    integer, intent(in) :: nocc, nmo
    real(kind=8) :: val
    integer :: i, a
    val = 0.0d0
    do i = 1, nocc
      do a = nocc+1, nmo
        val = val + lmo(a,i)*pmo(i,a)/(e(a)-e(i))
      end do
    end do
  end function sum_ov

!> @brief Contract a coupled response vector (occ-vir, length nocc*nvir, packed
!>  as k=(a-1)*nocc+i to match iatogen) with the PSO occ-vir block.
  function sum_resp(rt, pmo, nocc, nvir) result(val)
    real(kind=8), intent(in) :: rt(:)        ! (nocc*nvir)
    real(kind=8), intent(in) :: pmo(:,:)     ! PSO[s] in MO basis (nmo,nmo)
    integer, intent(in) :: nocc, nvir
    real(kind=8) :: val
    integer :: i, a, k
    val = 0.0d0
    k = 0
    do a = 1, nvir
      do i = 1, nocc
        k = k + 1
        val = val + rt(k)*pmo(i, nocc+a)
      end do
    end do
  end function sum_resp

!> @brief Build the antisymmetric AO first-order magnetic density from an
!>  occ-vir response vector: pa = C * (av - av^T) * C^T, av(occ,vir) = rin.
  subroutine magnetic_pb_density(mo, pa, av, nbf, nocc, rin)
    use tdhf_lib, only: iatogen
    use mathlib, only: orthogonal_transform
    real(kind=8), intent(in) :: mo(:,:)
    real(kind=8), intent(inout), target :: pa(:,:,:)
    real(kind=8), intent(inout) :: av(:,:)
    integer, intent(in) :: nbf, nocc
    real(kind=8), intent(in) :: rin(:)
    call iatogen(rin, av, nocc, nocc)
    av = av - transpose(av)
    call orthogonal_transform('t', nbf, mo, av, pa(:,:,1))
  end subroutine magnetic_pb_density

!> @brief Solve the coupled (CPHF/CPKS) magnetic response for the three field
!>  components and report the Phase-0 gate diagnostics.
!> @details Fixed-point solve of  (eps_a-eps_i) R + c_x*K[P^B(R)] = b,
!>  with b = orbital-Zeeman occ-vir block and K the exact-exchange image of the
!>  antisymmetric first-order density. Coulomb (J) and the semi-local XC kernel
!>  do not contribute (imaginary antisymmetric density). For c_x = 0 the loop is
!>  skipped and R = b/(eps_a-eps_i) (uncoupled).
  subroutine compute_coupled_para(infos, basis, mo, mo_e, nocc, nbf, nvir, &
                                  moL, scale_exch, rcoup, pb_asym, jnorm, knorm)
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_td_data_t, mntoia
    use types, only: information
    use basis_tools, only: basis_set
    implicit none
    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=8), intent(in) :: mo(:,:), mo_e(:), moL(:,:,:)
    integer, intent(in) :: nocc, nbf, nvir
    real(kind=8), intent(in) :: scale_exch
    real(kind=8), intent(out) :: rcoup(:,:)         ! (lvir,3)
    real(kind=8), intent(out) :: pb_asym, jnorm, knorm

    integer, parameter :: maxit = 100
    real(kind=8), parameter :: tol = 1.0d-9, half = 0.5d0
    real(kind=8), parameter :: kappa = 1.0d0        ! coupling sign (validated vs PySCF)
    type(int2_compute_t) :: int2_driver
    type(int2_td_data_t), target :: kdat, kdat1, jdat
    real(kind=8), allocatable, target :: pa(:,:,:)
    real(kind=8), allocatable :: av(:,:), bb(:,:), dd(:), gx(:), rprev(:), gao(:,:)
    integer :: lvir, t, i, a, k, it

    lvir = nocc*nvir
    allocate(pa(nbf,nbf,1), av(nbf,nbf), gao(nbf,nbf), &
             bb(lvir,3), dd(lvir), gx(lvir), rprev(lvir), source=0.0d0)

    ! RHS (orbital-Zeeman occ-vir block) and orbital-energy denominators
    do t = 1, 3
      k = 0
      do a = 1, nvir
        do i = 1, nocc
          k = k + 1
          bb(k,t) = moL(i, nocc+a, t)
        end do
      end do
    end do
    k = 0
    do a = 1, nvir
      do i = 1, nocc
        k = k + 1
        dd(k) = mo_e(nocc+a) - mo_e(i)
      end do
    end do

    do t = 1, 3
      rcoup(:,t) = bb(:,t)/dd                        ! uncoupled start
    end do

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    ! Coupled iteration (skipped for pure functionals, c_x = 0)
    if (abs(scale_exch) > 1.0d-12) then
      kdat = int2_td_data_t(d2=pa, int_apb=.false., int_amb=.true., &
                            tamm_dancoff=.false., scale_exchange=scale_exch)
      do t = 1, 3
        do it = 1, maxit
          rprev = rcoup(:,t)
          call magnetic_pb_density(mo, pa, av, nbf, nocc, rcoup(:,t))
          call int2_driver%run(kdat)
          gao = half*kdat%amb(:,:,1,1)
          call mntoia(gao, gx, mo, mo, nocc, nocc)
          rcoup(:,t) = (bb(:,t) - kappa*gx)/dd
          if (maxval(abs(rcoup(:,t)-rprev)) < tol) exit
        end do
      end do
    end if

    ! ---- Gate diagnostics on the converged z-component first-order density ----
    call magnetic_pb_density(mo, pa, av, nbf, nocc, rcoup(:,3))
    pb_asym = maxval(abs(pa(:,:,1) + transpose(pa(:,:,1))))      ! gate 0

    jdat = int2_td_data_t(d2=pa, int_apb=.true., int_amb=.false., &
                          tamm_dancoff=.false., scale_exchange=0.0d0)
    call int2_driver%run(jdat)
    jnorm = sqrt(sum((half*jdat%apb(:,:,1,1))**2))               ! gate 1 (Coulomb)

    kdat1 = int2_td_data_t(d2=pa, int_apb=.false., int_amb=.true., &
                           tamm_dancoff=.false., scale_exchange=1.0d0)
    call int2_driver%run(kdat1)
    knorm = sqrt(sum((half*kdat1%amb(:,:,1,1))**2))              ! gate 2 (exchange)

    deallocate(pa, av, gao, bb, dd, gx, rprev)
  end subroutine compute_coupled_para

end module nmr_shielding_mod
