!> @brief  Nonadiabatic molecular dynamics (NAMD) — Tully fewest-switches
!>         surface hopping (FSSH) core kernels for MRSF-TDDFT.
!>
!> @details
!>   Native OpenQP implementation of surface-hopping numerics, structured as
!>   clean, argument-based modern Fortran so the kernels are unit-testable and
!>   free of COMMON-block / dynamic-memory coupling.  The implementation covers:
!>
!>     - time-derivative couplings (TDC) supplied either from consecutive-step
!>       wavefunction overlaps or from the resident analytic derivative coupling
!>     - RK4 propagation of the electronic amplitudes
!>       i*hbar*\dot{c} = (E - i*sigma) c
!>     - cumulative Tully hopping probabilities
!>     - fewest-switches hop decision with selectable isotropic or analytic-NAC
!>       directional velocity rescaling (energy conservation)
!>     - kinetic energy
!>
!>   Internal-conversion accuracy upgrades, added per a verified literature
!>   survey (see session RESEARCH_ic_isc_methods.md):
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
!>   Deliberate implementation choices:
!>     * Everything is in consistent atomic units (energies in Hartree,
!>       velocities in bohr/atomic-time, masses in electron masses).  The
!>       FSSH and QM/MM paths use a single code path and the caller converts
!>       units once.
!>     * The O(nstate^2 * nsub) per-substep probability buffer of FSSHPR is
!>       dropped: probabilities are accumulated on the fly, then clamped and
!>       row-normalised once.
!>
!> @author  OpenQP development team
!> @date    2026-06
module namd_mod

  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_int64_t
  use precision, only: dp
  use state_tracking_mod, only: maximum_overlap_assignment

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
  public :: namd_rescale_velocities_directional
  public :: namd_fssh_decision
  public :: namd_decoherence_edc
  public :: namd_trivial_crossing
  public :: namd_counter_random
  public :: namd_counter_normal_fill
  public :: namd_baeck_an_tdc
  public :: namd_nacme_gate
  public :: namd_droplet_boundary
  public :: namd_com_restraint
  public :: namd_langevin_thermostat

  !> Default empirical decoherence constant C in the energy-based correction
  !> (Granucci & Persico, J. Chem. Phys. 126, 134114 (2007)), in Hartree.
  real(kind=dp), parameter, public :: NAMD_EDC_C_DEFAULT = 0.1_dp

contains

!> @brief Add two 64-bit bit patterns modulo 2^64 without signed overflow.
  pure function namd_add64_mod(a, b) result(value)
    integer(c_int64_t), intent(in) :: a, b
    integer(c_int64_t) :: value
    integer(c_int64_t) :: total, carry, limb
    integer(c_int64_t), parameter :: mask16 = int(z'FFFF', c_int64_t)
    integer(c_int64_t), parameter :: base16 = 65536_c_int64_t
    integer :: k

    value = 0_c_int64_t
    carry = 0_c_int64_t
    do k = 0, 3
      total = iand(shiftr(a, 16*k), mask16) + &
              iand(shiftr(b, 16*k), mask16) + carry
      limb = modulo(total, base16)
      carry = total/base16
      value = ior(value, shiftl(limb, 16*k))
    end do
  end function namd_add64_mod

!> @brief Multiply two 64-bit bit patterns modulo 2^64 without signed overflow.
  pure function namd_mul64_mod(a, b) result(value)
    integer(c_int64_t), intent(in) :: a, b
    integer(c_int64_t) :: value
    integer(c_int64_t) :: aa(0:3), bb(0:3)
    integer(c_int64_t) :: total, carry, limb
    integer(c_int64_t), parameter :: mask16 = int(z'FFFF', c_int64_t)
    integer(c_int64_t), parameter :: base16 = 65536_c_int64_t
    integer :: j, k

    do k = 0, 3
      aa(k) = iand(shiftr(a, 16*k), mask16)
      bb(k) = iand(shiftr(b, 16*k), mask16)
    end do
    value = 0_c_int64_t
    carry = 0_c_int64_t
    do k = 0, 3
      total = carry
      do j = 0, k
        total = total + aa(j)*bb(k - j)
      end do
      limb = modulo(total, base16)
      carry = total/base16
      value = ior(value, shiftl(limb, 16*k))
    end do
  end function namd_mul64_mod

!> @brief Stateless, counter-based uniform random number for NAMD.
!>
!> The returned value depends only on (seed, stream, step), not on call order,
!> process scheduling, or restart history.  The SplitMix64 finalizer is used as
!> a high-quality integer mixer; its upper 53 bits map exactly to an IEEE-754
!> double in [0,1).  The helper operations implement modulo-2^64 arithmetic
!> with 16-bit limbs, avoiding non-standard signed overflow under optimization.
!>
!> @param[in] seed    campaign/trajectory seed
!> @param[in] stream  independent trajectory stream identifier
!> @param[in] step    physical nuclear-step index
  pure function namd_counter_random(seed, stream, step) result(rand)
    integer(c_int64_t), intent(in) :: seed, stream, step
    real(kind=dp) :: rand
    integer(c_int64_t) :: z, mantissa
    integer(c_int64_t), parameter :: gamma = int(z'9E3779B97F4A7C15', c_int64_t)
    integer(c_int64_t), parameter :: stream_mix = int(z'D2B74407B1CE6E93', c_int64_t)
    integer(c_int64_t), parameter :: mix1 = int(z'BF58476D1CE4E5B9', c_int64_t)
    integer(c_int64_t), parameter :: mix2 = int(z'94D049BB133111EB', c_int64_t)
    real(kind=dp), parameter :: two_to_minus_53 = 1.0_dp/9007199254740992.0_dp

    z = namd_add64_mod(seed, namd_mul64_mod(gamma, &
        namd_add64_mod(step, 1_c_int64_t)))
    z = ieor(z, namd_mul64_mod(stream_mix, &
        namd_add64_mod(stream, 1_c_int64_t)))
    z = namd_mul64_mod(ieor(z, shiftr(z, 30)), mix1)
    z = namd_mul64_mod(ieor(z, shiftr(z, 27)), mix2)
    z = ieor(z, shiftr(z, 31))
    mantissa = shiftr(z, 11)
    rand = real(mantissa, dp)*two_to_minus_53
  end function namd_counter_random

!> @brief C ABI for the stateless NAMD counter RNG.
  function namd_counter_random_C(seed, stream, step) result(rand) &
      bind(C, name="oqp_namd_counter_random")
    integer(c_int64_t), value, intent(in) :: seed, stream, step
    real(c_double) :: rand
    rand = real(namd_counter_random(seed, stream, step), c_double)
  end function namd_counter_random_C

!> @brief Fill an array with deterministic standard-normal deviates.
!>
!> Negative counter values form a domain separate from the positive physical
!> MD steps used for hopping.  Box-Muller pairs are generated wholly in the
!> resident Fortran layer so rng_stream also separates Maxwell initial
!> velocities between trajectories.
  pure subroutine namd_counter_normal_fill(seed, stream, count, values) &
      bind(C, name="oqp_namd_counter_normal_fill")
    integer(c_int64_t), value, intent(in) :: seed, stream, count
    real(c_double), intent(out) :: values(*)
    integer :: i, nvalue
    real(kind=dp) :: u1, u2, radius, angle
    real(kind=dp), parameter :: two_pi = 6.28318530717958647692528676655900577_dp

    nvalue = int(count)
    do i = 1, nvalue, 2
      u1 = max(namd_counter_random(seed, stream, -int(i, c_int64_t)), tiny(1.0_dp))
      u2 = namd_counter_random(seed, stream, -int(i + 1, c_int64_t))
      radius = sqrt(-2.0_dp*log(u1))
      angle = two_pi*u2
      values(i) = radius*cos(angle)
      if (i + 1 <= nvalue) values(i + 1) = radius*sin(angle)
    end do
  end subroutine namd_counter_normal_fill

!> @brief Native finite spherical-droplet boundary energy and force.
!>
!> Atom membership is supplied once by the driver as group_index(a).  A zero
!> index excludes an atom; positive indices define restrained sites.  A
!> one-atom group gives an atom/oxygen wall, while a multi-atom group gives a
!> molecular-COM wall whose force is distributed by mass.  All arguments use
!> atomic units.  The zero-force region is r <= radius.  With x=r-radius and
!> t=x/buffer, the positive radial derivative is
!>
!>   dU/dr = k*x*(3*t**2 - 2*t**3),  0 < x < buffer
!>         = k*x,                     x >= buffer.
!>
!> Integrating that expression gives a C2 onset at the boundary and a
!> continuous match to a half-harmonic wall outside the buffer.  The returned
!> Cartesian array is force (-grad U), in atom-major xyz order.
  function namd_droplet_boundary(natom, ngroup, coordinates, masses, &
      group_index, center, radius, buffer, force_constant, &
      max_penetration_limit, energy, forces, max_penetration, active_count) &
      result(info) bind(C, name="oqp_namd_droplet_boundary")
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer(c_int64_t), value, intent(in) :: natom, ngroup
    real(c_double), intent(in) :: coordinates(*), masses(*), center(3)
    integer(c_int64_t), intent(in) :: group_index(*)
    real(c_double), value, intent(in) :: radius, buffer, force_constant
    real(c_double), value, intent(in) :: max_penetration_limit
    real(c_double), intent(out) :: energy, forces(*)
    real(c_double), intent(out) :: max_penetration
    integer(c_int64_t), intent(out) :: active_count
    integer(c_int) :: info
    real(kind=dp), allocatable :: group_mass(:), group_com(:,:), group_force(:,:)
    real(kind=dp) :: displacement(3), distance, penetration, t, switch
    real(kind=dp) :: radial_derivative, atom_weight
    integer :: a, g, k, na, ng

    info = -1_c_int
    energy = 0.0_dp
    max_penetration = 0.0_dp
    active_count = 0_c_int64_t
    if (natom <= 0_c_int64_t .or. ngroup <= 0_c_int64_t) return
    if (.not. ieee_is_finite(radius) .or. radius <= 0.0_dp .or. &
        .not. ieee_is_finite(buffer) .or. buffer < 0.0_dp .or. &
        .not. ieee_is_finite(force_constant) .or. force_constant <= 0.0_dp .or. &
        .not. ieee_is_finite(max_penetration_limit) .or. &
        max_penetration_limit < 0.0_dp .or. &
        any(.not. ieee_is_finite(center))) return

    na = int(natom)
    ng = int(ngroup)
    allocate(group_mass(ng), group_com(3, ng), group_force(3, ng))
    group_mass = 0.0_dp
    group_com = 0.0_dp
    group_force = 0.0_dp
    do a = 1, na
      do k = 1, 3
        forces(3*(a - 1) + k) = 0.0_dp
        if (.not. ieee_is_finite(coordinates(3*(a - 1) + k))) then
          info = -2_c_int
          return
        end if
      end do
      if (.not. ieee_is_finite(masses(a)) .or. masses(a) <= 0.0_dp) then
        info = -2_c_int
        return
      end if
      g = int(group_index(a))
      if (g < 0 .or. g > ng) then
        info = -3_c_int
        return
      end if
      if (g == 0) cycle
      group_mass(g) = group_mass(g) + masses(a)
      do k = 1, 3
        group_com(k, g) = group_com(k, g) + &
                          masses(a)*coordinates(3*(a - 1) + k)
      end do
    end do
    do g = 1, ng
      if (group_mass(g) <= 0.0_dp) then
        info = -3_c_int
        return
      end if
      group_com(:, g) = group_com(:, g)/group_mass(g)
      displacement = group_com(:, g) - center
      distance = norm2(displacement)
      if (.not. ieee_is_finite(distance)) then
        info = -2_c_int
        return
      end if
      penetration = max(0.0_dp, distance - radius)
      max_penetration = max(max_penetration, penetration)
      if (penetration > 0.0_dp) active_count = active_count + 1_c_int64_t
    end do
    ! Stop before evaluating an unphysically large polynomial.  The caller can
    ! checkpoint/report this configuration without propagating another step.
    if (max_penetration_limit > 0.0_dp .and. &
        max_penetration > max_penetration_limit) then
      info = 1_c_int
      return
    end if

    do g = 1, ng
      displacement = group_com(:, g) - center
      distance = norm2(displacement)
      penetration = distance - radius
      if (penetration <= 0.0_dp) cycle
      if (buffer > 0.0_dp .and. penetration < buffer) then
        t = penetration/buffer
        switch = 3.0_dp*t*t - 2.0_dp*t*t*t
        radial_derivative = force_constant*penetration*switch
        energy = energy + force_constant*buffer*buffer* &
                 (0.75_dp*t**4 - 0.4_dp*t**5)
      else
        radial_derivative = force_constant*penetration
        energy = energy + 0.5_dp*force_constant*penetration*penetration
        if (buffer > 0.0_dp) &
          energy = energy - 0.15_dp*force_constant*buffer*buffer
      end if
      if (.not. ieee_is_finite(energy) .or. &
          .not. ieee_is_finite(radial_derivative)) then
        info = -4_c_int
        return
      end if
      group_force(:, g) = -radial_derivative*displacement/distance
    end do
    do a = 1, na
      g = int(group_index(a))
      if (g == 0) cycle
      atom_weight = masses(a)/group_mass(g)
      do k = 1, 3
        forces(3*(a - 1) + k) = atom_weight*group_force(k, g)
      end do
    end do
    info = 0_c_int
  end function namd_droplet_boundary

!> @brief Harmonic restraint of a selected solute centre of mass.
!>
!> selected(a)=0 excludes atom a.  Force is distributed by mass so the
!> restrained coordinate and its analytic Cartesian derivative are identical.
  function namd_com_restraint(natom, coordinates, masses, selected, center, &
      force_constant, energy, forces, displacement_norm) result(info) &
      bind(C, name="oqp_namd_com_restraint")
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer(c_int64_t), value, intent(in) :: natom
    real(c_double), intent(in) :: coordinates(*), masses(*), center(3)
    integer(c_int64_t), intent(in) :: selected(*)
    real(c_double), value, intent(in) :: force_constant
    real(c_double), intent(out) :: energy, forces(*), displacement_norm
    integer(c_int) :: info
    real(kind=dp) :: total_mass, com(3), displacement(3), weight
    integer :: a, k, na

    info = -1_c_int
    energy = 0.0_dp
    displacement_norm = 0.0_dp
    if (natom <= 0_c_int64_t .or. .not. ieee_is_finite(force_constant) .or. &
        force_constant <= 0.0_dp .or. any(.not. ieee_is_finite(center))) return
    na = int(natom)
    total_mass = 0.0_dp
    com = 0.0_dp
    do a = 1, na
      do k = 1, 3
        forces(3*(a - 1) + k) = 0.0_dp
        if (.not. ieee_is_finite(coordinates(3*(a - 1) + k))) then
          info = -2_c_int
          return
        end if
      end do
      if (.not. ieee_is_finite(masses(a)) .or. masses(a) <= 0.0_dp) then
        info = -2_c_int
        return
      end if
      if (selected(a) == 0_c_int64_t) cycle
      total_mass = total_mass + masses(a)
      do k = 1, 3
        com(k) = com(k) + masses(a)*coordinates(3*(a - 1) + k)
      end do
    end do
    if (total_mass <= 0.0_dp) then
      info = -3_c_int
      return
    end if
    com = com/total_mass
    displacement = com - center
    displacement_norm = norm2(displacement)
    energy = 0.5_dp*force_constant*displacement_norm*displacement_norm
    if (.not. ieee_is_finite(energy)) then
      info = -4_c_int
      return
    end if
    do a = 1, na
      if (selected(a) == 0_c_int64_t) cycle
      weight = masses(a)/total_mass
      do k = 1, 3
        forces(3*(a - 1) + k) = -weight*force_constant*displacement(k)
      end do
    end do
    info = 0_c_int
  end function namd_com_restraint

!> @brief Apply one exact Langevin Ornstein-Uhlenbeck velocity step.
!>
!> The random vector is a pure function of seed, stream, physical MD step and
!> Cartesian component, so restart and process scheduling do not affect it.
!> heat is K_after-K_before: positive values are energy transferred from the
!> thermostat into the physical system.  The caller may subsequently project
!> constraints and should then replace heat by the measured constrained dK.
  function namd_langevin_thermostat(natom, dt, temperature, friction, seed, &
      stream, step, masses, velocities, heat) result(info) &
      bind(C, name="oqp_namd_langevin_thermostat")
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer(c_int64_t), value, intent(in) :: natom, seed, stream, step
    real(c_double), value, intent(in) :: dt, temperature, friction
    real(c_double), intent(in) :: masses(*)
    real(c_double), intent(inout) :: velocities(*)
    real(c_double), intent(out) :: heat
    integer(c_int) :: info
    integer :: a, i, ncomponent, na
    real(kind=dp) :: damping, thermal_scale, u1, u2, radius, angle
    real(kind=dp) :: kinetic_before, kinetic_after, normal1, normal2
    real(kind=dp), parameter :: kb_hartree = 3.166811563e-6_dp
    real(kind=dp), parameter :: two_pi = 6.28318530717958647692528676655900577_dp
    integer(c_int64_t) :: component

    info = -1_c_int
    heat = 0.0_dp
    if (natom <= 0_c_int64_t .or. .not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. ieee_is_finite(temperature) .or. temperature < 0.0_dp .or. &
        .not. ieee_is_finite(friction) .or. friction < 0.0_dp) return
    na = int(natom)
    ncomponent = 3*na
    kinetic_before = 0.0_dp
    do a = 1, na
      if (.not. ieee_is_finite(masses(a)) .or. masses(a) <= 0.0_dp) then
        info = -2_c_int
        return
      end if
      do i = 1, 3
        if (.not. ieee_is_finite(velocities(3*(a - 1) + i))) then
          info = -2_c_int
          return
        end if
        kinetic_before = kinetic_before + &
          0.5_dp*masses(a)*velocities(3*(a - 1) + i)**2
      end do
    end do
    damping = exp(-friction*dt)
    do i = 1, ncomponent, 2
      component = int(i, c_int64_t)
      u1 = max(namd_counter_component_random(seed, stream, step, component), &
               tiny(1.0_dp))
      u2 = namd_counter_component_random(seed, stream, step, component + 1_c_int64_t)
      radius = sqrt(-2.0_dp*log(u1))
      angle = two_pi*u2
      normal1 = radius*cos(angle)
      normal2 = radius*sin(angle)
      a = (i + 2)/3
      thermal_scale = sqrt(max(0.0_dp, (1.0_dp - damping*damping)* &
                           kb_hartree*temperature/masses(a)))
      velocities(i) = damping*velocities(i) + thermal_scale*normal1
      if (i + 1 <= ncomponent) then
        a = (i + 3)/3
        thermal_scale = sqrt(max(0.0_dp, (1.0_dp - damping*damping)* &
                             kb_hartree*temperature/masses(a)))
        velocities(i + 1) = damping*velocities(i + 1) + thermal_scale*normal2
      end if
    end do
    kinetic_after = 0.0_dp
    do a = 1, na
      do i = 1, 3
        kinetic_after = kinetic_after + &
          0.5_dp*masses(a)*velocities(3*(a - 1) + i)**2
      end do
    end do
    heat = kinetic_after - kinetic_before
    if (.not. ieee_is_finite(heat)) then
      info = -3_c_int
      return
    end if
    info = 0_c_int
  end function namd_langevin_thermostat

!> @brief Counter random variate keyed by a fourth component index.
  pure function namd_counter_component_random(seed, stream, step, component) &
      result(rand)
    integer(c_int64_t), intent(in) :: seed, stream, step, component
    real(kind=dp) :: rand
    integer(c_int64_t) :: z, mantissa
    integer(c_int64_t), parameter :: gamma = int(z'9E3779B97F4A7C15', c_int64_t)
    integer(c_int64_t), parameter :: stream_mix = int(z'D2B74407B1CE6E93', c_int64_t)
    integer(c_int64_t), parameter :: component_mix = int(z'CA5A826395121157', c_int64_t)
    integer(c_int64_t), parameter :: mix1 = int(z'BF58476D1CE4E5B9', c_int64_t)
    integer(c_int64_t), parameter :: mix2 = int(z'94D049BB133111EB', c_int64_t)
    real(kind=dp), parameter :: two_to_minus_53 = 1.0_dp/9007199254740992.0_dp

    z = namd_add64_mod(seed, namd_mul64_mod(gamma, &
        namd_add64_mod(step, 1_c_int64_t)))
    z = ieor(z, namd_mul64_mod(stream_mix, &
        namd_add64_mod(stream, 1_c_int64_t)))
    z = ieor(z, namd_mul64_mod(component_mix, &
        namd_add64_mod(component, 1_c_int64_t)))
    z = namd_mul64_mod(ieor(z, shiftr(z, 30)), mix1)
    z = namd_mul64_mod(ieor(z, shiftr(z, 27)), mix2)
    z = ieor(z, shiftr(z, 31))
    mantissa = shiftr(z, 11)
    rand = real(mantissa, dp)*two_to_minus_53
  end function namd_counter_component_random

!> @brief Time-dependent Baeck-An TDC from three consecutive energy points.
!>
!> The result is centred on energies_center.  dt_left spans old -> center and
!> dt_right spans center -> current, so the nonuniform three-point curvature is
!>
!>   f''(center) = 2 [dt_left*f(current)
!>                     -(dt_left+dt_right)*f(center)
!>                     +dt_right*f(old)]
!>                   / [dt_left*dt_right*(dt_left+dt_right)].
!>
!> For each pair, sigma_ij = sign(DeltaE_ij)/2 * sqrt(DeltaE''_ij/DeltaE_ij)
!> when the radicand is positive and the gap is inside gap_max.  The method is
!> a magnitude-only diagnostic: its sign convention cannot validate the
!> wavefunction gauge of an overlap-derived NACME.
  function namd_baeck_an_tdc(nstate, dt_left, dt_right, gap_max, &
      energies_old, energies_center, energies_current, tdc_row_major) &
      result(info) bind(C, name="oqp_namd_baeck_an_tdc")
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer(c_int64_t), value, intent(in) :: nstate
    real(c_double), value, intent(in) :: dt_left, dt_right, gap_max
    real(c_double), intent(in) :: energies_old(*), energies_center(*), &
                                  energies_current(*)
    real(c_double), intent(out) :: tdc_row_major(*)
    integer(c_int) :: info
    integer :: i, j, n
    real(kind=dp) :: denominator, gap_old, gap_center, gap_current
    real(kind=dp) :: curvature, radicand, sigma

    info = -1_c_int
    if (nstate <= 0_c_int64_t .or. dt_left <= 0.0_dp .or. &
        dt_right <= 0.0_dp .or. .not. ieee_is_finite(gap_max)) return
    n = int(nstate)
    denominator = dt_left*dt_right*(dt_left + dt_right)
    if (.not. ieee_is_finite(denominator) .or. denominator <= 0.0_dp) return

    do i = 1, n*n
      tdc_row_major(i) = 0.0_dp
    end do
    do i = 1, n
      if (.not. ieee_is_finite(energies_old(i)) .or. &
          .not. ieee_is_finite(energies_center(i)) .or. &
          .not. ieee_is_finite(energies_current(i))) then
        info = -2_c_int
        return
      end if
      do j = i + 1, n
        gap_old = energies_old(i) - energies_old(j)
        gap_center = energies_center(i) - energies_center(j)
        gap_current = energies_current(i) - energies_current(j)
        if (abs(gap_center) <= tiny(1.0_dp)) cycle
        if (gap_max > 0.0_dp .and. abs(gap_center) > gap_max) cycle
        curvature = 2.0_dp*(dt_left*gap_current &
                    -(dt_left + dt_right)*gap_center &
                    +dt_right*gap_old)/denominator
        radicand = curvature/gap_center
        if (.not. ieee_is_finite(radicand) .or. radicand <= 0.0_dp) cycle
        sigma = sign(0.5_dp, gap_center)*sqrt(radicand)
        tdc_row_major((i - 1)*n + j) = sigma
        tdc_row_major((j - 1)*n + i) = -sigma
      end do
    end do
    info = 0_c_int
  end function namd_baeck_an_tdc

!> @brief Provider-neutral validation gate for an MD time-derivative coupling.
!>
!> Exact matrix invariants (finite values, zero diagonal, antisymmetry) are
!> checked independently of the optional reference.  Masked state pairs are
!> then compared either by magnitude (compare_mode=0, appropriate for TD-BA)
!> or with their signed gauge (compare_mode=1, intended for phase-aligned
!> analytic d_ij dot velocity values).
!>
!> metrics = [candidate diagonal max, candidate antisymmetry max,
!>            reference diagonal max, reference antisymmetry max,
!>            pair RMS error, pair max error, max tolerance ratio]
!> counts  = [compared pairs, invariant failures, reference failures]
  function namd_nacme_gate(nstate, candidate, reference, reference_mask, &
      compare_mode, invariant_tol, abs_tol, rel_tol, metrics, counts) &
      result(info) bind(C, name="oqp_namd_nacme_gate")
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer(c_int64_t), value, intent(in) :: nstate
    real(c_double), intent(in) :: candidate(*), reference(*)
    integer(c_int), intent(in) :: reference_mask(*)
    integer(c_int), value, intent(in) :: compare_mode
    real(c_double), value, intent(in) :: invariant_tol, abs_tol, rel_tol
    real(c_double), intent(out) :: metrics(*)
    integer(c_int64_t), intent(out) :: counts(*)
    integer(c_int) :: info
    integer :: i, j, n, ij, ji
    real(kind=dp) :: cand_i, cand_j, ref_i, ref_j
    real(kind=dp) :: error, scale, sum_error2
    logical :: pair_valid

    info = -1_c_int
    if (nstate <= 0_c_int64_t .or. &
        .not. ieee_is_finite(invariant_tol) .or. &
        .not. ieee_is_finite(abs_tol) .or. &
        .not. ieee_is_finite(rel_tol) .or. &
        invariant_tol < 0.0_dp .or. abs_tol < 0.0_dp .or. &
        rel_tol < 0.0_dp) return
    if (compare_mode /= 0_c_int .and. compare_mode /= 1_c_int) return
    n = int(nstate)
    metrics(1:7) = 0.0_dp
    counts(1:3) = 0_c_int64_t
    sum_error2 = 0.0_dp

    do i = 1, n
      ij = (i - 1)*n + i
      if (.not. ieee_is_finite(candidate(ij))) then
        info = -2_c_int
        return
      end if
      metrics(1) = max(metrics(1), abs(candidate(ij)))
      if (abs(candidate(ij)) > invariant_tol) &
        counts(2) = counts(2) + 1_c_int64_t
      if (reference_mask(ij) /= 0_c_int) then
        if (.not. ieee_is_finite(reference(ij))) then
          info = -3_c_int
          return
        end if
        metrics(3) = max(metrics(3), abs(reference(ij)))
        if (abs(reference(ij)) > invariant_tol) &
          counts(2) = counts(2) + 1_c_int64_t
      end if
      do j = i + 1, n
        ij = (i - 1)*n + j
        ji = (j - 1)*n + i
        cand_i = candidate(ij)
        cand_j = candidate(ji)
        if (.not. ieee_is_finite(cand_i) .or. &
            .not. ieee_is_finite(cand_j)) then
          info = -2_c_int
          return
        end if
        metrics(2) = max(metrics(2), abs(cand_i + cand_j))
        if (abs(cand_i + cand_j) > invariant_tol) &
          counts(2) = counts(2) + 1_c_int64_t

        if ((reference_mask(ij) /= 0_c_int) .neqv. &
            (reference_mask(ji) /= 0_c_int)) then
          info = -4_c_int
          return
        end if
        pair_valid = reference_mask(ij) /= 0_c_int
        if (.not. pair_valid) cycle
        ref_i = reference(ij)
        ref_j = reference(ji)
        if (.not. ieee_is_finite(ref_i) .or. &
            .not. ieee_is_finite(ref_j)) then
          info = -3_c_int
          return
        end if
        metrics(4) = max(metrics(4), abs(ref_i + ref_j))
        if (abs(ref_i + ref_j) > invariant_tol) &
          counts(2) = counts(2) + 1_c_int64_t

        if (compare_mode == 0_c_int) then
          error = abs(abs(cand_i) - abs(ref_i))
        else
          error = abs(cand_i - ref_i)
        end if
        scale = abs_tol + rel_tol*abs(ref_i)
        counts(1) = counts(1) + 1_c_int64_t
        sum_error2 = sum_error2 + error*error
        metrics(6) = max(metrics(6), error)
        if (scale > 0.0_dp) then
          metrics(7) = max(metrics(7), error/scale)
        else if (error > 0.0_dp) then
          metrics(7) = huge(1.0_dp)
        end if
        if (error > scale) counts(3) = counts(3) + 1_c_int64_t
      end do
    end do
    if (counts(1) > 0_c_int64_t) &
      metrics(5) = sqrt(sum_error2/real(counts(1), dp))
    info = 0_c_int
  end function namd_nacme_gate

!> @brief Time-derivative (nonadiabatic) coupling from state overlaps.
!>        sigma(i,j) = ( S(i,j) - S(j,i) ) / (2 dt)
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
!>        i.e. \dot{c} = -(i E + sigma) c.
!>
!>   The returned increments are pre-multiplied by the integration step `h`
!>   so that k1..k4 are h*f.
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
!>        renormalisation.
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
!>        sub-step into the running cumulative matrix.
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
!>        and renormalise any row whose total exceeds one.
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
!>        v <- v * sqrt(1 + dE/KE),  dE = E_old - E_new.
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

!> @brief Conserve energy by changing velocity only along a coupling vector.
!>
!> For DeltaE = E_new - E_old, use
!>   v'_A = v_A + gamma*d_A/M_A,
!> where 0.5*A*gamma^2 + B*gamma + DeltaE = 0,
!> A = sum_A |d_A|^2/M_A and B = sum_A v_A.d_A.  The real root with
!> the smallest absolute gamma is selected.  A state-gauge sign or any nonzero
!> scalar multiplication of d therefore leaves the final velocity unchanged.
  subroutine namd_rescale_velocities_directional(vel, mass, direction, delta_e, &
      accepted, gamma, discriminant)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(kind=dp), intent(inout) :: vel(:,:)
    real(kind=dp), intent(in) :: mass(:), direction(:,:), delta_e
    logical, intent(out) :: accepted
    real(kind=dp), intent(out) :: gamma, discriminant
    integer :: a, nat
    real(kind=dp) :: avec, bvec, sqrt_disc, q, gamma1, gamma2, scale

    accepted = .false.
    gamma = 0.0_dp
    discriminant = -huge(1.0_dp)
    nat = size(mass)
    if (size(vel, 1) /= 3 .or. size(vel, 2) /= nat .or. &
        size(direction, 1) /= 3 .or. size(direction, 2) /= nat) return
    if (.not. ieee_is_finite(delta_e)) return

    avec = 0.0_dp
    bvec = 0.0_dp
    do a = 1, nat
      if (.not. ieee_is_finite(mass(a)) .or. mass(a) <= 0.0_dp) return
      if (.not. all(ieee_is_finite(vel(:,a))) .or. &
          .not. all(ieee_is_finite(direction(:,a)))) return
      avec = avec + sum(direction(:,a)**2)/mass(a)
      bvec = bvec + dot_product(vel(:,a), direction(:,a))
    end do
    if (.not. ieee_is_finite(avec) .or. avec <= tiny(1.0_dp)) return

    discriminant = bvec*bvec - 2.0_dp*avec*delta_e
    scale = bvec*bvec + abs(2.0_dp*avec*delta_e)
    if (discriminant < 0.0_dp .and. &
        abs(discriminant) <= 64.0_dp*epsilon(1.0_dp)*max(scale, tiny(1.0_dp))) &
      discriminant = 0.0_dp
    if (.not. ieee_is_finite(discriminant) .or. discriminant < 0.0_dp) return

    sqrt_disc = sqrt(discriminant)
    if (bvec >= 0.0_dp) then
      q = -bvec - sqrt_disc
    else
      q = -bvec + sqrt_disc
    end if
    if (abs(q) <= tiny(1.0_dp)) then
      if (abs(delta_e) > tiny(1.0_dp)) return
      gamma = 0.0_dp
    else
      gamma1 = q/avec
      gamma2 = 2.0_dp*delta_e/q
      if (abs(gamma1) <= abs(gamma2)) then
        gamma = gamma1
      else
        gamma = gamma2
      end if
    end if
    if (.not. ieee_is_finite(gamma)) return
    do a = 1, nat
      vel(:,a) = vel(:,a) + gamma*direction(:,a)/mass(a)
    end do
    accepted = .true.
  end subroutine namd_rescale_velocities_directional

!> @brief C ABI for the directional velocity-rescaling kernel.
  function namd_rescale_directional_C(natom, velocity, mass, direction, &
      delta_e, gamma, discriminant) result(info) &
      bind(C, name="oqp_namd_rescale_directional")
    integer(c_int64_t), value, intent(in) :: natom
    real(c_double), intent(inout) :: velocity(*)
    real(c_double), intent(in) :: mass(*), direction(*)
    real(c_double), value, intent(in) :: delta_e
    real(c_double), intent(out) :: gamma, discriminant
    integer(c_int) :: info
    integer :: a, nat
    real(kind=dp), allocatable :: vel2(:,:), dir2(:,:), mass2(:)
    logical :: accepted

    info = -1_c_int
    if (natom <= 0_c_int64_t) return
    nat = int(natom)
    allocate(vel2(3,nat), dir2(3,nat), mass2(nat))
    do a = 1, nat
      vel2(:,a) = velocity(3*a-2:3*a)
      dir2(:,a) = direction(3*a-2:3*a)
      mass2(a) = mass(a)
    end do
    call namd_rescale_velocities_directional(&
        vel2, mass2, dir2, real(delta_e, dp), accepted, gamma, discriminant)
    if (.not. accepted) then
      info = 1_c_int
      return
    end if
    do a = 1, nat
      velocity(3*a-2:3*a) = vel2(:,a)
    end do
    info = 0_c_int
  end function namd_rescale_directional_C

!> @brief Fewest-switches hop decision and (on accept) selectable isotropic or
!>        derivative-coupling directional velocity rescaling.
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
                                active, hopped, target, blocked, &
                                direction_vectors, rescale_mode, &
                                rescale_gamma, rescale_discriminant)
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
    real(kind=dp), intent(in), optional :: direction_vectors(:,:,:,:)
    integer, intent(in), optional :: rescale_mode
    real(kind=dp), intent(out), optional :: rescale_gamma, rescale_discriminant

    integer :: i, ncrst, n, mode
    real(kind=dp) :: lower, upper, de, ke, gamma, discriminant
    logical :: directional_ok

    n = size(eabs)
    ncrst = active
    hopped = .false.
    blocked = .false.
    target = active
    if (present(rescale_gamma)) rescale_gamma = 0.0_dp
    if (present(rescale_discriminant)) rescale_discriminant = -huge(1.0_dp)
    mode = 0
    if (present(rescale_mode)) mode = rescale_mode

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
        ! Energy-gap gate.
        if (abs(de) > thrshe) then
          blocked = .true.
          lower = upper
          cycle
        end if
        if (mode == 1) then
          if (.not. present(direction_vectors)) then
            blocked = .true.
            lower = upper
            cycle
          end if
          call namd_rescale_velocities_directional(&
              vel, mass, direction_vectors(:,:,i,ncrst), -de, &
              directional_ok, gamma, discriminant)
          if (present(rescale_gamma)) rescale_gamma = gamma
          if (present(rescale_discriminant)) &
            rescale_discriminant = discriminant
          if (.not. directional_ok) then
            blocked = .true.
            lower = upper
            cycle
          end if
        else
          ! Frustrated hop under isotropic rescaling: insufficient total KE.
          if (de < 0.0_dp .and. ke < abs(de)) then
            blocked = .true.
            lower = upper
            cycle
          end if
          call namd_rescale_velocities(vel, ke, de)
        end if
        ! Accept the hop after a successful energy-conserving velocity update.
        active = i
        target = i
        hopped = .true.
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
    integer, allocatable :: assignment(:)
    real(kind=dp), allocatable :: signs(:), matched(:), margins(:)
    integer :: info, n

    n = size(stas, 1)
    swapped = .false.
    if (abs(stas(active, active)) >= thresh) return   ! no trivial crossing

    ! Use the same global one-to-one state map as MO/X tracking.  Independent
    ! row maxima can assign several old states to one new state and lose the
    ! physical lineage at simultaneous or near-degenerate exchanges.
    allocate(assignment(n), signs(n), matched(n), margins(n))
    call maximum_overlap_assignment(stas, assignment, signs, matched, margins, info)
    if (info == 0 .and. assignment(active) /= active .and. &
        matched(active) >= thresh) then
      active = assignment(active)
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

!> @brief One Tully FSSH step: supplied TDC, RK4 amplitude propagation over
!>        sub-steps, optional EDC decoherence, trivial-crossing following, hop
!>        decision, and selectable isotropic or directional velocity rescaling.
!>
!>   Exchanges all NAMD state with the Python driver via flat tagarray records
!>   (1-D, layout-unambiguous):
!>     in : OQP_td_states_overlap (n x n), OQP_td_energies (n),
!>          OQP_namd_coef (2n: re1,im1,re2,im2,...), OQP_namd_velocity (3*nat),
!>          OQP_namd_params (>=15 packed scalars; slot 14 < 0 suppresses
!>          state changes while retaining coefficient propagation; slot 15
!>          selects isotropic or analytic-NAC velocity rescaling)
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
    integer :: tdc_scheme, rescale_mode
    real(kind=dp) :: dt_fs, dt_au, hsub, thrshe, rand, edc_c, triv_thr, ekin
    real(kind=dp) :: rescale_gamma, rescale_discriminant
    logical :: hopped, blocked, swapped, allow_hop

    real(kind=dp), allocatable :: tdc(:,:), cmhp(:,:), cr(:), ci(:), eabs(:), vel(:,:)
    real(kind=dp), allocatable :: mass_au(:), dcv(:,:,:,:)

    ! tagarray records
    real(kind=dp), contiguous, pointer :: stas_in(:), eabs_in(:), coef(:), velf(:), &
                                          params(:), results(:), tdc_in(:), dcv_in(:)
    real(kind=dp), contiguous, pointer :: mass(:)
    real(kind=dp), allocatable :: stas2(:,:)

    ! 1 atomic mass unit (Dalton) in electron masses
    real(kind=dp), parameter :: AMU_TO_AU = 1822.888486209_dp

    character(len=*), parameter :: subroutine_name = "namd_hop"
    ! NAMD state (energies/overlap/couplings) is supplied entirely through the
    ! namd_* tags, so the same kernel serves same-spin MRSF (n = tddft%nstate)
    ! and spin-adiabatic SOC NAMD (n = ns + 3*nt).
    character(len=*), parameter :: tags_req(*) = (/ character(len=80) :: &
        OQP_namd_coef, OQP_namd_velocity, OQP_namd_params, OQP_namd_tdc, &
        OQP_namd_eabs, OQP_namd_stas, OQP_namd_dcv /)
    character(len=*), parameter :: tags_out(*) = (/ character(len=80) :: &
        OQP_namd_results /)

    real(kind=dp), parameter :: FS_TO_AU = 41.341374575751_dp

    open(unit=iw, file=infos%log_filename, position="append")

    mass => infos%atoms%mass
    nat = size(mass)

    call data_has_tags(infos%dat, tags_req, module_name, subroutine_name, with_abort)
    call tagarray_get_data(infos%dat, OQP_namd_coef, coef)
    call tagarray_get_data(infos%dat, OQP_namd_velocity, velf)
    call tagarray_get_data(infos%dat, OQP_namd_params, params)
    call tagarray_get_data(infos%dat, OQP_namd_tdc, tdc_in)
    call tagarray_get_data(infos%dat, OQP_namd_eabs, eabs_in)
    call tagarray_get_data(infos%dat, OQP_namd_stas, stas_in)
    call tagarray_get_data(infos%dat, OQP_namd_dcv, dcv_in)

    ! number of states: from params (slot 13); fall back to tddft%nstate
    n = nint(params(13))
    if (n <= 0) n = int(infos%tddft%nstate)

    ! (re)allocate the results record. erase + alloc_or_die replaces the removed
    ! remove_records/reserve_data API (main's tagarray container refactor); erase
    ! drops any stale record so alloc_or_die always binds a fresh n*n+8 buffer.
    call infos%dat%erase(tags_out)
    call infos%dat%alloc_or_die(OQP_namd_results, (/ n*n + 8 /), results, &
         description=OQP_namd_results_comment)

    ! unpack parameters
    dt_fs       = params(1)
    nsub        = max(1, nint(params(2)))
    thrshe      = params(3)
    rand        = params(4)
    active      = nint(params(5))
    decoherence = nint(params(6))
    edc_c       = params(7)
    ! params(8) = TDC scheme (0 finite-diff, 1 NPI, 2 analytic v.d)
    tdc_scheme = nint(params(8))
    trivial_en  = nint(params(9))
    triv_thr    = params(10)
    allow_hop   = .true.
    if (size(params) >= 14) allow_hop = params(14) >= 0.0_dp
    ! params(15): 0 isotropic rescaling; 1 analytic derivative-coupling direction
    rescale_mode = 0
    if (size(params) >= 15) rescale_mode = nint(params(15))
    dt_au       = dt_fs*FS_TO_AU
    hsub        = dt_au/real(nsub, dp)

    allocate(tdc(n,n), cmhp(n,n), cr(n), ci(n), eabs(n), vel(3,nat), mass_au(nat), &
             stas2(n,n), dcv(3,nat,n,n))
    mass_au = mass*AMU_TO_AU        ! infos%atoms%mass is in amu; integrate in a.u.
    do i = 1, n
      cr(i) = coef(2*i-1)
      ci(i) = coef(2*i)
    end do
    eabs = eabs_in(1:n)                     ! absolute state energies (Hartree)
    do a = 1, nat
      vel(1,a) = velf(3*a-2)
      vel(2,a) = velf(3*a-1)
      vel(3,a) = velf(3*a)
    end do
    do i = 1, n
      do a = 1, n
        stas2(i,a) = stas_in((i-1)*n + a)   ! state overlap, flat row-major
      end do
    end do
    cmhp = 0.0_dp
    dcv = reshape(dcv_in(1:3*nat*n*n), shape(dcv))

    ! 1) follow diabatic character across trivial/unavoided crossings
    swapped = .false.
    if (allow_hop .and. trivial_en == 1) then
      call namd_trivial_crossing(stas2, triv_thr, active, swapped)
    end if

    ! 2) time-derivative couplings: supplied by the Python driver as a flat
    !    row-major (n x n) matrix (finite difference, norm-preserving
    !    interpolation, or analytic v.d). tdc(i,j) = tdc_in((i-1)*n + j).
    !    Preserve an exactly zero analytic/NPI matrix: zero can be the physical
    !    answer. The legacy in-Fortran FD fallback applies only to FD mode.
    do i = 1, n
      do a = 1, n
        tdc(i,a) = tdc_in((i-1)*n + a)
      end do
    end do
    if (tdc_scheme == 0 .and. all(abs(tdc) < 1.0e-30_dp)) &
      call namd_state_tdc(stas2, dt_au, tdc)

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

    ! 5) fewest-switches hop + selected energy-conserving velocity rescaling
    hopped = .false.
    blocked = .false.
    target = active
    rescale_gamma = 0.0_dp
    rescale_discriminant = -huge(1.0_dp)
    if (allow_hop) then
      call namd_fssh_decision(cmhp, eabs, rand, thrshe, mass_au, vel, &
                              active, hopped, target, blocked, dcv, rescale_mode, &
                              rescale_gamma, rescale_discriminant)
    end if

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
    results(n*n+6) = real(rescale_mode, dp)
    results(n*n+7) = rescale_gamma
    results(n*n+8) = rescale_discriminant

    deallocate(tdc, cmhp, cr, ci, eabs, vel, mass_au, stas2, dcv)
    close(iw)
  end subroutine namd_hop

end module namd_mod
