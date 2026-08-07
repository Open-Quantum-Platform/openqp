"""OQP frequency analysis """

import numpy as np
from oqp.utils.constants import SPEED_OF_LIGHT, ATMOS, BOHR
from oqp.utils.constants import FREQ_TO_INV_CM, AMU_to_KG, J_TO_AU
from oqp.utils.constants import GAS_CONSTANT, PLANCK_CONSTANT, BOLTZMANN_CONSTANT, AVOGADRO_CONSTANT

# Common atomic-unit IR conversion used by quantum-chemistry frequency analyses:
# (d dipole / d normal coordinate)^2 -> km/mol.  Dipole derivatives are in
# e*Bohr and normal coordinates use amu-normalized Cartesian modes.
IR_INTENSITY_CONVERSION_KM_MOL = 42.255


def infrared_intensities(dipole_derivatives, modes):
    """Project Cartesian dipole derivatives onto normal modes.

    Parameters
    ----------
    dipole_derivatives
        Array with shape ``(3, 3N)`` containing d(mu_x,mu_y,mu_z)/dR in
        atomic units.
    modes
        Normal modes with shape ``(nmode, 3N)`` as returned by ``normal_mode``.

    Returns
    -------
    intensities, mode_dipoles
        IR intensities in km/mol and the mode-projected dipole derivatives.
    """

    dipole_derivatives = np.asarray(dipole_derivatives, dtype=float)
    modes = np.asarray(modes, dtype=float)
    if dipole_derivatives.ndim != 2 or dipole_derivatives.shape[0] != 3:
        raise ValueError("dipole_derivatives must have shape (3, 3N)")
    if modes.ndim != 2 or modes.shape[1] != dipole_derivatives.shape[1]:
        raise ValueError(
            f"modes must have shape (nmode, {dipole_derivatives.shape[1]}), got {modes.shape}"
        )

    mode_dipoles = modes @ dipole_derivatives.T
    intensities = IR_INTENSITY_CONVERSION_KM_MOL * np.einsum('ij,ij->i', mode_dipoles, mode_dipoles)
    return intensities, mode_dipoles


def raman_activities(polarizability_derivatives, modes):
    """Project Cartesian polarizability derivatives and compute Raman activities.

    ``polarizability_derivatives`` has shape ``(3, 3, 3N)`` in atomic units.
    The returned activities follow the standard non-resonant expression
    ``45 * alpha_bar'^2 + 7 * gamma'^2`` for each normal mode.
    """

    polarizability_derivatives = np.asarray(polarizability_derivatives, dtype=float)
    modes = np.asarray(modes, dtype=float)
    if polarizability_derivatives.shape[:2] != (3, 3) or polarizability_derivatives.ndim != 3:
        raise ValueError("polarizability_derivatives must have shape (3, 3, 3N)")
    if modes.ndim != 2 or modes.shape[1] != polarizability_derivatives.shape[2]:
        raise ValueError(
            f"modes must have shape (nmode, {polarizability_derivatives.shape[2]}), got {modes.shape}"
        )

    mode_polarizabilities = np.tensordot(modes, polarizability_derivatives, axes=(1, 2))
    mode_polarizabilities = 0.5 * (mode_polarizabilities + np.swapaxes(mode_polarizabilities, 1, 2))
    axx = mode_polarizabilities[:, 0, 0]
    ayy = mode_polarizabilities[:, 1, 1]
    azz = mode_polarizabilities[:, 2, 2]
    axy = mode_polarizabilities[:, 0, 1]
    ayz = mode_polarizabilities[:, 1, 2]
    azx = mode_polarizabilities[:, 2, 0]
    alpha_bar = (axx + ayy + azz) / 3.0
    gamma2 = 0.5 * ((axx - ayy) ** 2 + (ayy - azz) ** 2 + (azz - axx) ** 2) + 3.0 * (
        axy ** 2 + ayz ** 2 + azx ** 2
    )
    activities = 45.0 * alpha_bar ** 2 + 7.0 * gamma2
    return activities, mode_polarizabilities


def normal_mode(coord, mass, hessian):
    # compute normal mode
    # coord in Bohr, mass in g/mol, hessian in Hartree/Bohr**2
    natom = len(mass)
    ncoord = len(coord)
    xyz = coord.reshape((natom, 3))

    # compute mass-weighted hessian
    mr = np.repeat(mass, 3)
    mw_hess = hessian/np.outer(mr, mr) ** 0.5

    # compute center of mass
    com = np.sum(xyz * mass.reshape((natom, 1)) / np.sum(mass), axis=0)

    # compute inertial
    it = np.zeros((3, 3))

    ## compute momentum of inertia
    for n, i in enumerate(xyz - com):
        it += mass[n] * (np.sum(i ** 2) * np.diag(np.ones(3)) - np.outer(i, i))

    # compute principal moments and axis
    p_mom, p_axis = np.linalg.eigh(it)

    # rotate coordinate to principal coordinate system
    p_xyz = np.dot((xyz - com), p_axis)

    # compute modes for transitions and rotations
    tr = np.zeros((ncoord, 6))
    tr[0::3, 0] = mass ** 0.5
    tr[1::3, 1] = mass ** 0.5
    tr[2::3, 2] = mass ** 0.5

    for i, mi in enumerate(mass):
        mij = mi ** 0.5
        for j in range(3):
            tr[3*i+j, 3] = + mij * (p_xyz[i, 1] * p_axis[j, 2] - p_xyz[i, 2] * p_axis[j, 1])
            tr[3*i+j, 4] = - mij * (p_xyz[i, 0] * p_axis[j, 2] + p_xyz[i, 2] * p_axis[j, 0])
            tr[3*i+j, 5] = + mij * (p_xyz[i, 0] * p_axis[j, 1] - p_xyz[i, 1] * p_axis[j, 0])

    u, s, v = np.linalg.svd(tr, full_matrices=True)
    n_external = 3
    if natom > 1:
        linear = np.count_nonzero(p_mom > 1.0e-8) < 3
        n_external = 5 if linear else 6
    b = u[:, n_external:]

    h, u3 = np.linalg.eigh(np.dot(b.T, np.dot(mw_hess, b)))
    q = np.dot(b, u3)

    # compute frequencies
    freqs = np.zeros_like(h)
    freqs[h >= 0] = h[h >= 0] ** 0.5
    freqs[h < 0] = -(-h[h < 0]) ** 0.5
    freqs *= FREQ_TO_INV_CM

    # compute normal modes
    modes = q / np.outer(mr ** 0.5, np.ones((q.shape[1],)))

    return freqs, modes.T, p_mom

def thermal_analysis(
        energy, atoms, mass, freqs, inertia,
        temperature=298.15,
        linear=False,
        sigma=1,
        mult=0,
        freq_scale_factor=1,
        freq_cutoff=100,
):
    # -------- Remove imaginary freqs ---------
    freqs = freqs[freqs > 0]

    # ---------------- Damping ----------------
    alpha = 4
    damp = 1 / (1 + (freq_cutoff / freqs) ** alpha)

    # ---------------- ZPE ----------------
    # zero-pont energy, R * Sum(0.5 * hv/k) in Hartree
    hv_kt = PLANCK_CONSTANT * freqs * SPEED_OF_LIGHT * freq_scale_factor / (BOLTZMANN_CONSTANT * temperature)
    zpe = 0.5 * np.sum(hv_kt) * GAS_CONSTANT * temperature * J_TO_AU

    # ---------------- Enthalpy ----------------
    # electronic energy in Hartree
    el = energy

    # translational energy, 3/2 RT, in Hartree
    u_trans = 1.5 * GAS_CONSTANT * temperature * J_TO_AU

    # rotational energy, 0 (atomic) ; RT (linear); 3/2 RT (non-linear) in Hartree
    if len(atoms) == 1:
        u_rot = 0.0
    else:
        if linear:
            u_rot = GAS_CONSTANT * temperature * J_TO_AU
        else:
            u_rot = 1.5 * GAS_CONSTANT * temperature * J_TO_AU

    # vibrational energy, R * Sum((hv/k)/(e^(hv/KT)-1)), in Hartree
    hv_kt_ehv_kt = hv_kt / (np.exp(hv_kt) - 1.0)
    u_vib = np.sum(hv_kt_ehv_kt) * GAS_CONSTANT * temperature * J_TO_AU

    # Quasi-rigid rotor harmonic oscillator energy
    # 1/2(Nhv) + RT(hv/kT)e^(-hv/kT)/(1-e^(-hv/kT))
    # u_rrho_vib = GAS_CONSTANT * temperature * hv_kt * np.exp(-hv_kt) / (1 - np.exp(-hv_kt))
    # u_vib = np.sum(damp * u_rrho_vib + (1 - damp) * 0.5 * GAS_CONSTANT * temperature)

    # pV = RT in Hartree
    pv = GAS_CONSTANT * temperature * J_TO_AU

    # ---------------- Entropy ----------------
    # electronic entropy * T in Hartree
    st_el = GAS_CONSTANT * np.log(mult) * temperature * J_TO_AU

    # translational entropy, R(Ln(2πmkT/h^2)^3/2(1/C)) + 1 + 3/2) * T in Hartree
    r = (2.0 * np.pi * np.sum(mass) * AMU_to_KG * BOLTZMANN_CONSTANT * temperature) ** 0.5 / PLANCK_CONSTANT
    free_space = 1000  # get_free_space(solv)
    conc = ATMOS / (GAS_CONSTANT * temperature)  # g/L
    den = conc * 1000 * AVOGADRO_CONSTANT / (free_space / 1000.0)
    st_trans = GAS_CONSTANT * (2.5 + np.log(r ** 3 / den)) * temperature * J_TO_AU

    # rotational entropy, 0 (atomic) ; R(Ln(q)+1) * T (linear); R(Ln(q)+3/2) * T (non-linear) in Hartree
    # A linear rotor has a vanishing principal moment, so this legitimately
    # produces a meaningless rotational constant; both consumers below select
    # against it explicitly. Silence the divide warning rather than let it
    # reach the user.
    with np.errstate(divide='ignore', invalid='ignore'):
        rc = PLANCK_CONSTANT / (8 * np.pi ** 2 * inertia * AMU_to_KG * BOHR ** 2)
    rt = rc * PLANCK_CONSTANT / BOLTZMANN_CONSTANT

    # Select the physically meaningful moments on the INERTIA, at the same
    # threshold normal_mode uses to decide `linear` -- not on rc/rt.
    #
    # eigh returns the vanishing moment as exactly 0.0 only when the molecular
    # axis happens to lie along a coordinate axis. Tilt it and you get a tiny
    # value of either sign, and testing rc/rt instead lets it through:
    # CO2 along (1,1,1) gives -5.4e-31, whose huge negative rotational constant
    # makes S_vib NaN; CO2 tilted in the xy plane gives +5.3e-15, which passes
    # an isfinite-and-positive test and drives S_rot to -62.26 cal/mol/K
    # instead of +13.07.
    significant = np.asarray(inertia, dtype=float) > 1.0e-8

    # sigma is the rotational symmetry number: the rotational partition function
    # counts each indistinguishable orientation once, so q_rot carries a 1/sigma
    # that was missing here. Omitting it inflates S_rot by R*ln(sigma) -- 1.38
    # cal/mol/K for water, 4.94 for benzene -- and every S and G derived from it.
    sigma = max(1, int(sigma))

    if len(atoms) == 1:
        st_rot = 0.0
    elif linear:
        # A linear rotor has one vanishing principal moment, so rt holds one
        # infinity; the old code indexed rt[0], which numpy sorts ascending and
        # is therefore exactly that infinity -- every diatomic printed -inf.
        # Use the doubly degenerate finite rotational temperature instead.
        rot_temps = rt[significant]
        theta = float(rot_temps[0]) if rot_temps.size else np.nan
        qr = temperature / (sigma * theta)
        st_rot = GAS_CONSTANT * (np.log(qr) + 1) * temperature * J_TO_AU
    else:
        qr = (np.pi * temperature ** 3 / (rt[0] * rt[1] * rt[2])) ** 0.5 / sigma
        st_rot = GAS_CONSTANT * (np.log(qr) + 1.5) * temperature * J_TO_AU

    # vibrational entropy,
    # quasi-harmonic model, Grimme, S. (2012), Chem. Eur. J., 18: 9955-9964. DOI.org/10.1002/chem.201200497
    # rigid-rotor Sv = RSum(hv/kT/(e^(hv/kT)-1) - ln(1-e^(-hv/kT))) * T in Hartree
    s_rrho_vib = hv_kt / (np.exp(hv_kt) - 1) - np.log(1 - np.exp(-hv_kt))
    st_rrho_vib = s_rrho_vib * GAS_CONSTANT * temperature * J_TO_AU

    # free-rotor Sv = R(1/2 + 1/2ln((8π^3u'kT/h^2)) * T in Hartree
    # Grimme's free-rotor term needs an average rotational constant. For a
    # linear rotor the vanishing principal moment contributes a meaningless
    # one -- infinite when the moment is exactly zero, enormous and possibly
    # negative when it is not -- and averaging it in gives bav ~ 0, hence
    # log(0) and st_vib = -inf or NaN. That is why linear molecules printed
    # -inf for the TOTAL entropy even once st_rot was fixed.
    rot_consts = np.asarray(rc, dtype=float)[significant]
    av_rc = float(np.sum(rot_consts) / rot_consts.size) if rot_consts.size else np.inf  # s-1
    bav = PLANCK_CONSTANT / av_rc  # kg m^2
    mu = PLANCK_CONSTANT / (8 * np.pi ** 2 * freqs * SPEED_OF_LIGHT * freq_scale_factor)
    mu_p = bav / (mu + bav)
    log_mu = 8 * np.pi ** 3 * mu_p * BOLTZMANN_CONSTANT * temperature / PLANCK_CONSTANT ** 2
    st_frho_vib = (0.5 + np.log(log_mu ** 0.5)) * GAS_CONSTANT * temperature * J_TO_AU

    # combine vibrational entropy
    st_vib = np.sum(damp * st_rrho_vib + (1 - damp) * st_frho_vib)

    thermo_data = {
        'temp': temperature,
        'mass': float(np.sum(mass)),
        'rc': rc/SPEED_OF_LIGHT,
        'rt': rt,
        'el': el,
        'zpe': zpe,
        'u_trans': u_trans,
        'u_rot': u_rot,
        'u_vib': u_vib,
        'pv': pv,
        'st_el': st_el,
        'st_trans': st_trans,
        'st_rot': st_rot,
        'st_vib': st_vib,
        'sigma': int(sigma),
        'linear': bool(linear),
        # The printer must mark the SAME moments the entropy code selected
        # against, or the two disagree about which rotational constant is
        # meaningful. Testing isfinite() there is not equivalent: a tilted
        # linear rotor's vanishing moment comes back as a tiny finite number,
        # which is exactly the case this PR exists to handle.
        'rot_significant': [bool(x) for x in significant],
        'monatomic': len(atoms) == 1,
    }

    return thermo_data
