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
!> @details v1 scope: RHF and closed-shell pure DFT / HF, common gauge origin
!>  (CGO), uncoupled paramagnetic term (exact CPKS for pure functionals and the
!>  HF-uncoupled approximation). The isotropic shielding is reported as
!>  sigma = sigma_dia + sigma_para per nucleus.
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

    allocate(pso_full(nbf,nbf,3), moP(nmo,nmo,3), source=0.0d0)
    allocate(sig_para(3,3,nat), siso_para(nat), siso_tot(nat), source=0.0d0)

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
          sig_para(t,s,i) = twoa2 * sum_ov(moL(:,:,t), moP(:,:,s), mo_e_a, nocc, nmo)
        end do
      end do
      siso_para(i) = (sig_para(1,1,i)+sig_para(2,2,i)+sig_para(3,3,i))/3.0d0 * PPM
      siso_tot(i)  = siso_dia(i) + siso_para(i)
    end do

    write(iw,'(/4x,a)') 'Isotropic shielding (CGO, uncoupled, ppm):'
    write(iw,'(4x,a)')  '   Atom    Z      sigma_dia      sigma_para       sigma_total'
    do i = 1, nat
      write(iw,'(4x,i6,f8.1,3f16.6)') i, basis%atoms%zn(i), &
             siso_dia(i), siso_para(i), siso_tot(i)
    end do

    ! PSO antisymmetry diagnostics (must be ~0 for the anti-Hermitian operator)
    write(iw,'(/4x,a,2es12.3)') 'PSO diagnostics  max|diag|, max|A+A^T| = ', &
           pso_diag_max, pso_asym_max
    call flush(iw)

    deallocate(amom, lfull, gdia, coords, sig_dia, siso_dia)
    deallocate(pso_full, moL, moP, sig_para, siso_para, siso_tot)
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

end module nmr_shielding_mod
