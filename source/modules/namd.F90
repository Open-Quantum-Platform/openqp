!> @brief  Nonadiabatic molecular dynamics (NAMD) — Tully fewest-switches
!>         surface hopping (FSSH) core kernels for MRSF-TDDFT.
!>
!> @details
!>   Faithful port of the surface-hopping numerics from the GAMESS `namd.src`
!>   module (S. Lee), restructured into clean, argument-based modern Fortran so
!>   the kernels are unit-testable and free of COMMON-block / dynamic-memory
!>   coupling.  The physics mirrors the original exactly:
!>
!>     - time-derivative couplings (TDC) from wavefunction overlaps
!>       between consecutive nuclear steps        [GAMESS NACVFD]
!>     - RK4 propagation of the electronic amplitudes
!>       i*hbar*\dot{c} = (E - i*sigma) c          [GAMESS PPTDECOE/NDDTCR/NDDTCC]
!>     - cumulative Tully hopping probabilities    [GAMESS FSSHPRST/FSSHPR]
!>     - fewest-switches hop decision + isotropic
!>       velocity rescaling (energy conservation)  [GAMESS FSSH/FSSHT/RESCALV]
!>     - kinetic energy                            [GAMESS MDQKIN]
!>
!>   Internal-conversion accuracy upgrades, added per a verified literature
!>   survey (see session RESEARCH_ic_isc_methods.md) and absent from the GAMESS
!>   reference:
!>     - energy-based decoherence correction (EDC) — Granucci & Persico,
!>       J. Chem. Phys. 126, 134114 (2007); the SHARC default decoherence scheme
!>     - trivial / unavoided-crossing detection with diabatic state following,
!>       in the spirit of SC-FSSH — Wang & Prezhdo, JPCL 5, 713 (2014)
!>
!>   Planned next (documented, not yet implemented here):
!>     - norm-preserving interpolation (NPI) time-derivative couplings
!>       (Meek & Levine, JPCL 5, 2351 (2014)): rigorous multistate form is the
!>       real antisymmetric matrix logarithm of the Loewdin-orthonormalised
!>       step overlap, T = logm(orth(S))/dt, which reduces to the exact 2-state
!>       identity T*dt = arcsin(S_10).  Will replace namd_state_tdc when wired.
!>     - intersystem crossing (ISC) via the SHARC spin-adiabatic representation:
!>       diagonalise H = H_MCH + H_SOC, hop on the diagonal states, propagate
!>       c_diag = U' . P_MCH . U . c_diag.  Requires MRSF Breit-Pauli SOC
!>       matrix elements as input.
!>
!>   Deliberate, documented deviations from the original ("the GAMESS code may
!>   not be perfect"):
!>     * Everything is in consistent atomic units (energies in Hartree,
!>       velocities in bohr/atomic-time, masses in electron masses).  The
!>       original mixed Hartree (FSSH) and kcal/mol (FSSHT, QM/MM) paths; here a
!>       single code path is used and the caller converts units once.
!>     * The O(nstate^2 * nsub) per-substep probability buffer of FSSHPR is
!>       dropped: probabilities are accumulated on the fly, then clamped and
!>       row-normalised once — numerically identical to the original sum.
!>
!> @author  Port: OpenQP NAMD; original algorithm: Seunghoon Lee (GAMESS)
!> @date    2026-06
module namd_mod

  use precision, only: dp

  implicit none

  private

  character(len=*), parameter :: module_name = "namd_mod"

  public :: namd_state_tdc
  public :: namd_coeff_deriv
  public :: namd_propagate_coeff
  public :: namd_accumulate_hop_prob
  public :: namd_finalize_hop_prob
  public :: namd_kinetic_energy
  public :: namd_rescale_velocities
  public :: namd_fssh_decision
  public :: namd_decoherence_edc
  public :: namd_trivial_crossing

  !> Default empirical decoherence constant C in the energy-based correction
  !> (Granucci & Persico, J. Chem. Phys. 126, 134114 (2007)), in Hartree.
  real(kind=dp), parameter, public :: NAMD_EDC_C_DEFAULT = 0.1_dp

contains

!> @brief Time-derivative (nonadiabatic) coupling from state overlaps.
!>        sigma(i,j) = ( S(i,j) - S(j,i) ) / (2 dt)          [GAMESS NACVFD]
!>
!> @param[in]  stas   nstate x nstate overlap <Phi_i(t-dt)|Phi_j(t)> between
!>                     the previous and current nuclear geometries
!> @param[in]  dt     nuclear time step (atomic time units)
!> @param[out] tdc    nstate x nstate antisymmetric time-derivative coupling
  subroutine namd_state_tdc(stas, dt, tdc)
    real(kind=dp), intent(in)  :: stas(:,:)
    real(kind=dp), intent(in)  :: dt
    real(kind=dp), intent(out) :: tdc(:,:)
    integer :: i, j, n
    n = size(stas, 1)
    do j = 1, n
      do i = 1, n
        tdc(i,j) = (stas(i,j) - stas(j,i)) / (2.0_dp*dt)
      end do
    end do
  end subroutine namd_state_tdc

!> @brief Right-hand side of the electronic equation of motion in the adiabatic
!>        basis (amplitudes c = cr + i*ci):
!>           \dot{cr}_k = - sum_i sigma(k,i) cr_i + E_k ci_k
!>           \dot{ci}_k = - sum_i sigma(k,i) ci_i - E_k cr_k
!>        i.e. \dot{c} = -(i E + sigma) c.        [GAMESS NDDTCR/NDDTCC]
!>
!>   The returned increments are pre-multiplied by the integration step `h`
!>   (matching the original convention where k1..k4 are h*f).
  subroutine namd_coeff_deriv(cr, ci, tdc, eig, h, dcr, dci)
    real(kind=dp), intent(in)  :: cr(:), ci(:)
    real(kind=dp), intent(in)  :: tdc(:,:)
    real(kind=dp), intent(in)  :: eig(:)
    real(kind=dp), intent(in)  :: h
    real(kind=dp), intent(out) :: dcr(:), dci(:)
    integer :: k, i, n
    real(kind=dp) :: sr, si
    n = size(cr)
    do k = 1, n
      sr = 0.0_dp
      si = 0.0_dp
      do i = 1, n
        sr = sr - tdc(k,i)*cr(i)
        si = si - tdc(k,i)*ci(i)
      end do
      dcr(k) = (sr + eig(k)*ci(k)) * h
      dci(k) = (si - eig(k)*cr(k)) * h
    end do
  end subroutine namd_coeff_deriv

!> @brief One RK4 sub-step of the electronic amplitudes, followed by
!>        renormalisation.                         [GAMESS PPTDECOE]
!>
!> @param[in,out] cr,ci  real/imaginary amplitudes (nstate)
!> @param[in]     tdc    time-derivative coupling (nstate x nstate, constant
!>                       over the nuclear step)
!> @param[in]     eig    absolute adiabatic state energies (Hartree)
!> @param[in]     h      electronic sub-step length (atomic time units)
  subroutine namd_propagate_coeff(cr, ci, tdc, eig, h)
    real(kind=dp), intent(inout) :: cr(:), ci(:)
    real(kind=dp), intent(in)    :: tdc(:,:)
    real(kind=dp), intent(in)    :: eig(:)
    real(kind=dp), intent(in)    :: h
    integer :: n, k
    real(kind=dp), allocatable :: k1r(:), k1i(:), k2r(:), k2i(:)
    real(kind=dp), allocatable :: k3r(:), k3i(:), k4r(:), k4i(:)
    real(kind=dp), allocatable :: tr(:), ti(:)
    real(kind=dp) :: dnorm

    n = size(cr)
    allocate(k1r(n), k1i(n), k2r(n), k2i(n), k3r(n), k3i(n), k4r(n), k4i(n), &
             tr(n), ti(n))

    call namd_coeff_deriv(cr, ci, tdc, eig, h, k1r, k1i)
    tr = cr + 0.5_dp*k1r;  ti = ci + 0.5_dp*k1i
    call namd_coeff_deriv(tr, ti, tdc, eig, h, k2r, k2i)
    tr = cr + 0.5_dp*k2r;  ti = ci + 0.5_dp*k2i
    call namd_coeff_deriv(tr, ti, tdc, eig, h, k3r, k3i)
    tr = cr + k3r;         ti = ci + k3i
    call namd_coeff_deriv(tr, ti, tdc, eig, h, k4r, k4i)

    do k = 1, n
      cr(k) = cr(k) + (k1r(k) + 2.0_dp*k2r(k) + 2.0_dp*k3r(k) + k4r(k))/6.0_dp
      ci(k) = ci(k) + (k1i(k) + 2.0_dp*k2i(k) + 2.0_dp*k3i(k) + k4i(k))/6.0_dp
    end do

    dnorm = sqrt(sum(cr*cr) + sum(ci*ci))
    if (dnorm > 0.0_dp) then
      cr = cr/dnorm
      ci = ci/dnorm
    end if

    deallocate(k1r, k1i, k2r, k2i, k3r, k3i, k4r, k4i, tr, ti)
  end subroutine namd_propagate_coeff

!> @brief Accumulate the Tully transition probability over one electronic
!>        sub-step into the running cumulative matrix.   [GAMESS FSSHPRST]
!>
!>        g(i,j) += 2 sigma(i,j) Re(c_i^* c_j) h / |c_i|^2
!>
!>   Call once per sub-step (after propagating the amplitudes), then finalise
!>   with namd_finalize_hop_prob.
  subroutine namd_accumulate_hop_prob(cmhp, cr, ci, tdc, h)
    real(kind=dp), intent(inout) :: cmhp(:,:)
    real(kind=dp), intent(in)    :: cr(:), ci(:)
    real(kind=dp), intent(in)    :: tdc(:,:)
    real(kind=dp), intent(in)    :: h
    integer :: i, j, n
    real(kind=dp) :: pii
    n = size(cr)
    do i = 1, n
      pii = cr(i)*cr(i) + ci(i)*ci(i)
      if (pii <= 0.0_dp) cycle
      do j = 1, n
        cmhp(i,j) = cmhp(i,j) &
          + 2.0_dp*tdc(i,j)*(cr(i)*cr(j) + ci(i)*ci(j))*h/pii
      end do
    end do
  end subroutine namd_accumulate_hop_prob

!> @brief Finalise cumulative hopping probabilities: clamp negatives to zero
!>        and renormalise any row whose total exceeds one.   [GAMESS FSSHPR]
  subroutine namd_finalize_hop_prob(cmhp)
    real(kind=dp), intent(inout) :: cmhp(:,:)
    integer :: i, j, n
    real(kind=dp) :: rowsum
    n = size(cmhp, 1)
    do i = 1, n
      do j = 1, n
        if (cmhp(i,j) < 0.0_dp) cmhp(i,j) = 0.0_dp
      end do
      rowsum = sum(cmhp(i,:))
      if (rowsum > 1.0_dp) cmhp(i,:) = cmhp(i,:)/rowsum
    end do
  end subroutine namd_finalize_hop_prob

!> @brief Classical kinetic energy  KE = 1/2 sum_a m_a |v_a|^2  (atomic units).
!>                                                          [GAMESS MDQKIN]
!> @param[in] vel   3 x natom velocities (bohr / atomic-time)
!> @param[in] mass  natom atomic masses (electron masses)
  pure function namd_kinetic_energy(vel, mass) result(ke)
    real(kind=dp), intent(in) :: vel(:,:)
    real(kind=dp), intent(in) :: mass(:)
    real(kind=dp) :: ke
    integer :: a, nat
    nat = size(mass)
    ke = 0.0_dp
    do a = 1, nat
      ke = ke + mass(a)*(vel(1,a)**2 + vel(2,a)**2 + vel(3,a)**2)
    end do
    ke = 0.5_dp*ke
  end function namd_kinetic_energy

!> @brief Isotropic velocity rescaling after a hop to conserve total energy.
!>        v <- v * sqrt(1 + dE/KE),  dE = E_old - E_new.    [GAMESS RESCALV]
!>
!>   Caller must already have verified the hop is energetically allowed
!>   (KE >= |dE| when dE < 0); otherwise the argument of sqrt is negative.
  subroutine namd_rescale_velocities(vel, ke, de)
    real(kind=dp), intent(inout) :: vel(:,:)
    real(kind=dp), intent(in)    :: ke   !< kinetic energy on the old surface
    real(kind=dp), intent(in)    :: de   !< E_old - E_new (Hartree)
    real(kind=dp) :: scale
    if (ke <= 0.0_dp) return
    scale = sqrt(max(0.0_dp, 1.0_dp + de/ke))
    vel = scale*vel
  end subroutine namd_rescale_velocities

!> @brief Fewest-switches hop decision and (on accept) isotropic velocity
!>        rescaling.                                  [GAMESS FSSH/FSSHT]
!>
!>   All energies in Hartree, velocities/masses in atomic units.
!>
!> @param[in]     cmhp     finalised cumulative hop probabilities (nstate^2);
!>                         row `active` is used
!> @param[in]     eabs     absolute adiabatic state energies (Hartree)
!> @param[in]     rand     random number in [0,1)
!> @param[in]     thrshe   energy-gap gate: hops with |dE| > thrshe are blocked
!>                         (Hartree). Use a large value (e.g. huge) to disable.
!> @param[in]     mass     atomic masses (natom)
!> @param[in,out] vel      3 x natom velocities; rescaled in place on a hop
!> @param[in,out] active   active state index (1..nstate); updated on a hop
!> @param[out]    hopped   .true. if a hop occurred
!> @param[out]    target   state hopped to (= active on no hop)
!> @param[out]    blocked  .true. if a candidate hop was rejected (frustrated
!>                         or gated)
  subroutine namd_fssh_decision(cmhp, eabs, rand, thrshe, mass, vel, &
                                active, hopped, target, blocked)
    real(kind=dp), intent(in)    :: cmhp(:,:)
    real(kind=dp), intent(in)    :: eabs(:)
    real(kind=dp), intent(in)    :: rand
    real(kind=dp), intent(in)    :: thrshe
    real(kind=dp), intent(in)    :: mass(:)
    real(kind=dp), intent(inout) :: vel(:,:)
    integer,       intent(inout) :: active
    logical,       intent(out)   :: hopped
    integer,       intent(out)   :: target
    logical,       intent(out)   :: blocked

    integer :: i, ncrst, n
    real(kind=dp) :: lower, upper, de, ke

    n = size(eabs)
    ncrst = active
    hopped = .false.
    blocked = .false.
    target = active

    ! Walk the cumulative probability ladder over candidate target states.
    ! Self-transition probability is identically zero (sigma(i,i)=0), so the
    ! cumulative sum can include the diagonal without effect.
    lower = 0.0_dp
    do i = 1, n
      if (i == ncrst) then
        lower = lower + cmhp(ncrst, i)   ! adds 0; keeps ladder aligned
        cycle
      end if
      upper = lower + cmhp(ncrst, i)
      if (rand > lower .and. rand < upper) then
        de = eabs(ncrst) - eabs(i)       ! E_old - E_new
        ke = namd_kinetic_energy(vel, mass)
        ! Frustrated hop: not enough kinetic energy to climb uphill.
        if (de < 0.0_dp .and. ke < abs(de)) then
          blocked = .true.
          lower = upper
          cycle
        end if
        ! Energy-gap gate.
        if (abs(de) > thrshe) then
          blocked = .true.
          lower = upper
          cycle
        end if
        ! Accept the hop.
        active = i
        target = i
        hopped = .true.
        call namd_rescale_velocities(vel, ke, de)
        return
      end if
      lower = upper
    end do
  end subroutine namd_fssh_decision

!> @brief Energy-based decoherence correction (EDC).
!>        Granucci & Persico, J. Chem. Phys. 126, 134114 (2007); the pragmatic
!>        default in SHARC. Damps the non-active amplitudes toward zero on the
!>        decoherence time scale and restores the total norm via the active
!>        state:
!>           tau_k = (1/|E_k - E_a|) (1 + C/E_kin)     (atomic units, hbar=1)
!>           c_k  <- c_k exp(-dt/tau_k)        for k /= a
!>           c_a  <- c_a sqrt( (1 - sum_{k/=a}|c_k|^2) / |c_a|^2 )
!>
!>   Apply once per nuclear step, after the electronic propagation.
!>
!> @param[in,out] cr,ci  amplitudes (nstate)
!> @param[in]     eabs   absolute adiabatic state energies (Hartree)
!> @param[in]     active active state index
!> @param[in]     ekin   nuclear kinetic energy (Hartree)
!> @param[in]     dt     nuclear time step (atomic time units)
!> @param[in]     cval   empirical constant C (Hartree); see NAMD_EDC_C_DEFAULT
  subroutine namd_decoherence_edc(cr, ci, eabs, active, ekin, dt, cval)
    real(kind=dp), intent(inout) :: cr(:), ci(:)
    real(kind=dp), intent(in)    :: eabs(:)
    integer,       intent(in)    :: active
    real(kind=dp), intent(in)    :: ekin
    real(kind=dp), intent(in)    :: dt
    real(kind=dp), intent(in)    :: cval
    integer :: k, n
    real(kind=dp) :: gap, tau, decay, pa, sum_others, scale
    real(kind=dp), parameter :: tiny = 1.0e-12_dp

    n = size(cr)
    if (ekin <= 0.0_dp) return        ! no kinetic energy -> no decoherence
    sum_others = 0.0_dp
    do k = 1, n
      if (k == active) cycle
      gap = abs(eabs(k) - eabs(active))
      if (gap < tiny) then            ! (near-)degenerate: skip damping
        sum_others = sum_others + cr(k)*cr(k) + ci(k)*ci(k)
        cycle
      end if
      tau   = (1.0_dp/gap)*(1.0_dp + cval/ekin)
      decay = exp(-dt/tau)
      cr(k) = cr(k)*decay
      ci(k) = ci(k)*decay
      sum_others = sum_others + cr(k)*cr(k) + ci(k)*ci(k)
    end do
    pa = cr(active)*cr(active) + ci(active)*ci(active)
    if (pa > tiny) then
      scale = sqrt(max(0.0_dp, 1.0_dp - sum_others)/pa)
      cr(active) = cr(active)*scale
      ci(active) = ci(active)*scale
    end if
  end subroutine namd_decoherence_edc

!> @brief Trivial / unavoided-crossing detection and diabatic state following.
!>        Practical local-diabatization fix in the spirit of SC-FSSH
!>        (Wang & Prezhdo, J. Phys. Chem. Lett. 5, 713 (2014)).
!>
!>   At a trivial (non-interacting) crossing two adiabatic labels swap between
!>   consecutive steps: the active state's self-overlap |S(a,a)| collapses while
!>   |S(a,j)| ~ 1 for the partner j.  Following the diabatic character (relabel
!>   active -> j) prevents the spurious "hop far from the crossing" that plain
!>   FSSH suffers on dense PES.  No velocity rescaling is applied: at a trivial
!>   crossing the energy is continuous along the diabatic state.
!>
!> @param[in]     stas      nstate x nstate state overlap S(i,j)=<i(t-dt)|j(t)>
!> @param[in]     thresh    self-overlap threshold below which a crossing is
!>                          flagged (e.g. 0.5)
!> @param[in,out] active    active state index; relabelled on a trivial crossing
!> @param[out]    swapped   .true. if a relabel occurred
  subroutine namd_trivial_crossing(stas, thresh, active, swapped)
    real(kind=dp), intent(in)    :: stas(:,:)
    real(kind=dp), intent(in)    :: thresh
    integer,       intent(inout) :: active
    logical,       intent(out)   :: swapped
    integer :: j, n, jmax
    real(kind=dp) :: amax

    n = size(stas, 1)
    swapped = .false.
    if (abs(stas(active, active)) >= thresh) return   ! no trivial crossing

    ! Partner = state with the largest |overlap| to the (old) active state.
    jmax = active
    amax = abs(stas(active, active))
    do j = 1, n
      if (j == active) cycle
      if (abs(stas(active, j)) > amax) then
        amax = abs(stas(active, j))
        jmax = j
      end if
    end do
    if (jmax /= active .and. amax >= thresh) then
      active = jmax
      swapped = .true.
    end if
  end subroutine namd_trivial_crossing

!> @brief C-interoperable entry: one FSSH surface-hopping step for MRSF-TDDFT.
!>        Driven from the Python NAMD trajectory loop after the per-step
!>        electronic structure (energies, response vectors, phase-corrected
!>        state overlap) has been computed.
  subroutine namd_hop_C(c_handle) bind(C, name="mrsf_namd_hop")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call namd_hop(inf)
  end subroutine namd_hop_C

!> @brief One Tully FSSH step: TDC from the state overlap, RK4 amplitude
!>        propagation over sub-steps, optional EDC decoherence, trivial-crossing
!>        following, hop decision and isotropic velocity rescaling.
!>
!>   Exchanges all NAMD state with the Python driver via flat tagarray records
!>   (1-D, layout-unambiguous):
!>     in : OQP_td_states_overlap (n x n), OQP_td_energies (n),
!>          OQP_namd_coef (2n: re1,im1,re2,im2,...), OQP_namd_velocity (3*nat),
!>          OQP_namd_params (>=12 packed scalars)
!>     out: OQP_namd_coef, OQP_namd_velocity (rescaled), OQP_namd_params(active,
!>          hopped, target), OQP_namd_results (n*n cumulative probs + flags)
  subroutine namd_hop(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use messages, only: show_message, with_abort

    implicit none

    type(information), target, intent(inout) :: infos

    integer :: n, nat, i, a, isub, nsub, active, target, decoherence, trivial_en
    real(kind=dp) :: dt_fs, dt_au, hsub, thrshe, rand, edc_c, triv_thr, ekin
    logical :: hopped, blocked, swapped

    real(kind=dp), allocatable :: tdc(:,:), cmhp(:,:), cr(:), ci(:), eabs(:), vel(:,:)
    real(kind=dp), allocatable :: mass_au(:)

    ! tagarray records
    real(kind=dp), contiguous, pointer :: stas(:,:), omega(:), coef(:), velf(:), &
                                          params(:), results(:)
    real(kind=dp), contiguous, pointer :: mass(:)

    ! 1 atomic mass unit (Dalton) in electron masses
    real(kind=dp), parameter :: AMU_TO_AU = 1822.888486209_dp

    character(len=*), parameter :: subroutine_name = "namd_hop"
    character(len=*), parameter :: tags_req(*) = (/ character(len=80) :: &
        OQP_td_states_overlap, OQP_td_energies, OQP_namd_coef, &
        OQP_namd_velocity, OQP_namd_params /)
    character(len=*), parameter :: tags_out(*) = (/ character(len=80) :: &
        OQP_namd_results /)

    real(kind=dp), parameter :: FS_TO_AU = 41.341374575751_dp

    open(unit=iw, file=infos%log_filename, position="append")

    n   = int(infos%tddft%nstate)
    mass => infos%atoms%mass
    nat = size(mass)

    call data_has_tags(infos%dat, tags_req, module_name, subroutine_name, with_abort)
    call tagarray_get_data(infos%dat, OQP_td_states_overlap, stas)
    call tagarray_get_data(infos%dat, OQP_td_energies, omega)
    call tagarray_get_data(infos%dat, OQP_namd_coef, coef)
    call tagarray_get_data(infos%dat, OQP_namd_velocity, velf)
    call tagarray_get_data(infos%dat, OQP_namd_params, params)

    ! (re)allocate the results record
    call infos%dat%remove_records(tags_out)
    call infos%dat%reserve_data(OQP_namd_results, ta_type_real64, &
         n*n + 8, (/ n*n + 8 /), comment=OQP_namd_results_comment)
    call tagarray_get_data(infos%dat, OQP_namd_results, results)

    ! unpack parameters
    dt_fs       = params(1)
    nsub        = max(1, nint(params(2)))
    thrshe      = params(3)
    rand        = params(4)
    active      = nint(params(5))
    decoherence = nint(params(6))
    edc_c       = params(7)
    ! params(8) = tdc scheme (0 finite-diff, 1 NPI) -- NPI pending
    trivial_en  = nint(params(9))
    triv_thr    = params(10)
    dt_au       = dt_fs*FS_TO_AU
    hsub        = dt_au/real(nsub, dp)

    allocate(tdc(n,n), cmhp(n,n), cr(n), ci(n), eabs(n), vel(3,nat), mass_au(nat))
    mass_au = mass*AMU_TO_AU        ! infos%atoms%mass is in amu; integrate in a.u.
    do i = 1, n
      cr(i) = coef(2*i-1)
      ci(i) = coef(2*i)
    end do
    eabs = omega(1:n)                      ! excitation energies; ref cancels
    do a = 1, nat
      vel(1,a) = velf(3*a-2)
      vel(2,a) = velf(3*a-1)
      vel(3,a) = velf(3*a)
    end do
    cmhp = 0.0_dp

    ! 1) follow diabatic character across trivial/unavoided crossings
    swapped = .false.
    if (trivial_en == 1) call namd_trivial_crossing(stas, triv_thr, active, swapped)

    ! 2) time-derivative couplings from the phase-corrected state overlap
    call namd_state_tdc(stas, dt_au, tdc)

    ! 3) propagate amplitudes over electronic sub-steps; accumulate hop flux
    do isub = 1, nsub
      call namd_propagate_coeff(cr, ci, tdc, eabs, hsub)
      call namd_accumulate_hop_prob(cmhp, cr, ci, tdc, hsub)
    end do
    call namd_finalize_hop_prob(cmhp)

    ! 4) decoherence (energy-based correction)
    ekin = namd_kinetic_energy(vel, mass_au)
    if (decoherence == 1) &
      call namd_decoherence_edc(cr, ci, eabs, active, ekin, dt_au, edc_c)

    ! 5) fewest-switches hop + isotropic velocity rescaling
    call namd_fssh_decision(cmhp, eabs, rand, thrshe, mass_au, vel, &
                            active, hopped, target, blocked)

    ! pack results back
    do i = 1, n
      coef(2*i-1) = cr(i)
      coef(2*i)   = ci(i)
    end do
    do a = 1, nat
      velf(3*a-2) = vel(1,a)
      velf(3*a-1) = vel(2,a)
      velf(3*a)   = vel(3,a)
    end do
    params(5)  = real(active, dp)
    params(11) = merge(1.0_dp, 0.0_dp, hopped)
    params(12) = real(target, dp)

    results = 0.0_dp
    results(1:n*n) = reshape(cmhp, (/ n*n /))
    results(n*n+1) = merge(1.0_dp, 0.0_dp, hopped)
    results(n*n+2) = real(target, dp)
    results(n*n+3) = merge(1.0_dp, 0.0_dp, blocked)
    results(n*n+4) = ekin
    results(n*n+5) = merge(1.0_dp, 0.0_dp, swapped)

    deallocate(tdc, cmhp, cr, ci, eabs, vel, mass_au)
    close(iw)
  end subroutine namd_hop

end module namd_mod
