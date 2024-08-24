"""OQP frequency analysis """

import numpy as np
from oqp.utils.constants import SPEED_OF_LIGHT, ATMOS, BOHR
from oqp.utils.constants import FREQ_TO_INV_CM, AMU_to_KG, J_TO_AU
from oqp.utils.constants import GAS_CONSTANT, PLANCK_CONSTANT, BOLTZMANN_CONSTANT, AVOGADRO_CONSTANT

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
    b = u[:, 6:]

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
    rc = PLANCK_CONSTANT / (8 * np.pi ** 2 * inertia * AMU_to_KG * BOHR ** 2)
    rt = rc * PLANCK_CONSTANT / BOLTZMANN_CONSTANT

    if len(atoms) == 1:
        st_rot = 0.0
    else:
        if len(atoms) == 2:
            qr = temperature / rt[0]
        else:
            qr = (np.pi * temperature ** 3 / (rt[0] * rt[1] * rt[2])) ** 0.5

        if linear:
            st_rot = GAS_CONSTANT * (np.log(qr) + 1) * temperature * J_TO_AU
        else:
            st_rot = GAS_CONSTANT * (np.log(qr) + 1.5) * temperature * J_TO_AU

    # vibrational entropy,
    # quasi-harmonic model, Grimme, S. (2012), Chem. Eur. J., 18: 9955-9964. DOI.org/10.1002/chem.201200497
    # rigid-rotor Sv = RSum(hv/kT/(e^(hv/kT)-1) - ln(1-e^(-hv/kT))) * T in Hartree
    s_rrho_vib = hv_kt / (np.exp(hv_kt) - 1) - np.log(1 - np.exp(-hv_kt))
    st_rrho_vib = s_rrho_vib * GAS_CONSTANT * temperature * J_TO_AU

    # free-rotor Sv = R(1/2 + 1/2ln((8π^3u'kT/h^2)) * T in Hartree
    av_rc = sum(rc) / len(rc)  # s-1
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
    }

    return thermo_data
