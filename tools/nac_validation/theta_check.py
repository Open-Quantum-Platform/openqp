"""Semi-numerical theta^x decomposition of the MRSF NAC orbital term.

For each Cartesian coordinate x, compute the inter-geometry MO overlap
M_pq(+/-) = <phi_p(R0)|phi_q(R0 +/- delta e_x)> (SCF at the displaced
geometry, phases aligned to R0), giving

    theta^x_pq = <phi_p|d_x phi_q> = (M+ - M-)/(2 delta)
    Sket^x_pq  = <phi_p| (d_x chi)-part of phi_q>  (from the cross AO overlap)
    U^x        = theta^x - Sket^x

The exact orbital coupling is d^orb_x = sum_pq gamma^IJ_pq theta^x_pq with NO
code conventions involved.  Compare:
  (1) d_CI + d^orb(semi-num)  vs  d_num        -> is the decomposition complete?
  (2) frozen+skeleton terms (Fortran mrsf_nac_overlap)  vs  their semi-num
      counterparts                              -> are the coded terms right?
  (3) U-contraction (semi-num) vs the code's z-vector D -> the exact scale.
"""
import sys
import copy
import oqp                       # import BEFORE numpy (ILP64 LAPACK)
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import NAC, SinglePoint
import numpy as np

inp = sys.argv[1] if len(sys.argv) > 1 else '/tmp/nactest/H2O_energy.inp'
DELTA = 1.0e-3
PAIR = (2, 3)

r = Runner(input_file=inp, log=inp.replace('.inp', '_theta.log'))
r.run()
mol = r.mol
nac = NAC(mol)
nstate, natom = nac.nstate, nac.natom
nbf = None

# ---------- reference (R0) data ----------
xyz0 = np.array(mol.get_system(), copy=True)            # flat 3N, bohr
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)  # C-wrap of F(nbf,nbf)
nbf = C0raw.shape[0]
# Fortran mo_a(mu, p): python wrap is its transpose -> C0py[p, mu]
C0py = C0raw.copy()
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
C0raw_b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0_b = np.array(mol.data['OQP::E_MO_B'], copy=True)

# numerical NAC + analytical pieces at R0 (subprocesses; parent state intact)
nacv_n, dcv_n, fn = nac.numerical_nac()
oqp.mrsf_nac_overlap(mol)
ovraw = np.array(mol.data['OQP::nac_overlap'], copy=True)
ovF = ovraw.reshape(-1).reshape((nstate, nstate, 3 * natom))
ov = np.transpose(ovF, (1, 0, 2)).reshape((nstate, nstate, natom, 3))
# MO transition densities, Fortran (nbf*nbf, nstate, nstate)
tdraw = np.array(mol.data['OQP::nac_trden_mo'], copy=True)
tdF = tdraw.reshape(-1).reshape((nstate, nstate, nbf * nbf))


def gamma_mo(i, j):
    # flat index k = (p-1) + (q-1)*nbf -> C reshape gives [q,p] -> transpose
    return tdF[j - 1, i - 1].reshape((nbf, nbf)).T.copy()


nacv_a, dcv_a, fa = nac.analytical_nac()

i, j = PAIR
dn = dcv_n[i - 1, j - 1].reshape(-1)
da = dcv_a[i - 1, j - 1].reshape(-1)
dov = ov[i - 1, j - 1].reshape(-1)
dci = da - dov
gam = gamma_mo(i, j)                       # gamma_f[p, q]


def grad_now():
    oqp.tdhf_mrsf_gradient(mol)
    return mol.get_grad().reshape((natom, 3)).copy()


def cphf_term(ii, jj, block=0):
    mol.data.set_tdhf_target(ii)
    oqp.set_mrsf_nac_cphf(mol, ii, jj)
    oqp.set_mrsf_nac_cphf_block(mol, block)
    oqp.tdhf_mrsf_z_vector(mol)
    if not mol.mol_energy.Z_Vector_converged:
        oqp.set_mrsf_nac_cphf(mol, 0, 0)
        oqp.set_mrsf_nac_cphf_block(mol, 0)
        return None
    gZ = grad_now()
    mol.data['OQP::td_p'] = np.zeros_like(np.array(mol.data['OQP::td_p'], copy=True))
    mol.data['OQP::WAO'] = np.zeros_like(np.array(mol.data['OQP::WAO'], copy=True))
    gS = grad_now()
    oqp.set_mrsf_nac_cphf(mol, 0, 0)
    oqp.set_mrsf_nac_cphf_block(mol, 0)
    return (gZ - gS).reshape(-1)


D_code = cphf_term(i, j)

# occupation-space boundaries
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
print(f'noca={noca} nocb={nocb} nbf={nbf}')


def spaces(p):                              # 0-based orbital index
    if p < nocb:
        return 0                            # doc
    if p < noca:
        return 1                            # socc
    return 2                                # virt


SP = np.array([spaces(p) for p in range(nbf)])

# ---------- displaced SCF + cross overlaps ----------
cfg = mol.config
json0 = mol.log.replace('.log', '.json')
mol.save_data()
cfg['guess']['type'] = 'json'
cfg['guess']['file'] = json0
cfg['guess']['continue_geom'] = False


def cross_overlaps(coord):
    """SCF+TD at geometry `coord`; return (M_f, Sao_f, X_f) vs R0:
    M_f[p,q] = <phi_p(R0)|phi_q(R)>, Sao_f[mu,nu] = <chi_mu(R0)|chi_nu(R)>,
    X_f[ij, st] = response amplitudes at R (MO-phase corrected)."""
    mol.update_system(coord)
    oqp.library.ints_1e(mol)
    oqp.library.guess(mol)
    SinglePoint(mol).energy()
    # "old" = R0 reference
    mol.data['OQP::xyz_old'] = xyz0.reshape((3, natom))
    mol.data['OQP::VEC_MO_A_old'] = C0raw
    mol.data['OQP::E_MO_A_old'] = e0
    mol.data['OQP::VEC_MO_B_old'] = C0raw_b
    mol.data['OQP::E_MO_B_old'] = e0_b
    oqp.get_structures_ao_overlap(mol)
    Mraw = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
    Araw = np.array(mol.data['OQP::overlap_ao_non_orthogonal'], copy=True)
    # python wrap of Fortran (nbf,nbf) = transpose; overlap_mo_f[p_old, q_new]
    M_f = Mraw.T.copy()
    Sao_f = Araw.T.copy()
    # phase-fix new MOs to R0: flip columns with negative diagonal overlap
    sgn = np.sign(np.diag(M_f))
    sgn[sgn == 0] = 1.0
    M_f = M_f * sgn[None, :]
    if np.min(np.abs(np.diag(M_f))) < 0.9:
        print(f'  WARNING: weak diagonal MO overlap '
              f'{np.min(np.abs(np.diag(M_f))):.3f} (orbital swap?)')
    # response amplitudes at the displaced geometry, in the SAME folded
    # layout as R0: bvF[ij, st].  Correct each amplitude element for the MO
    # phase flips (element ij pairs occ i with virt a: ij=(a-1)*noca+i).
    bvraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
    nij = noca * (nbf - nocb)
    nst = bvraw.size // nij
    X_f = bvraw.reshape(-1).reshape((nst, nij)).T.copy()
    for ij in range(nij):
        a = ij // noca + nocb              # 0-based virt orbital
        ii = ij % noca                     # 0-based occ orbital
        X_f[ij, :] *= sgn[a] * sgn[ii]
    return M_f, Sao_f, X_f


# R0 amplitudes in folded layout bv[ij, st]
bvraw0 = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
nij = noca * (nbf - nocb)
X0 = bvraw0.reshape(-1).reshape((nstate, nij)).T.copy()

ncoord = 3 * natom
all_theta = np.zeros((ncoord, nbf, nbf))
all_sket = np.zeros((ncoord, nbf, nbf))
d_amp_true = np.zeros(ncoord)
d_orb_true = np.zeros(ncoord)
d_frozen_true = np.zeros(ncoord)
d_skel_true = np.zeros(ncoord)     # -1/2 same-space - cross(lo,hi), gamma.S^[x]
d_U_true = np.zeros(ncoord)
d_Ucphf_true = np.zeros(ncoord)    # [gamma(hi,lo)-gamma(lo,hi)] U(hi,lo)
maxcon = 0.0

for k in range(ncoord):
    Mp, Sp_, Xp = cross_overlaps(xyz0 + DELTA * np.eye(ncoord)[k])
    Mm, Sm_, Xm = cross_overlaps(xyz0 - DELTA * np.eye(ncoord)[k])
    # align displaced state phases (and check assignment) against R0
    for X in (Xp, Xm):
        for st in range(nstate):
            o = np.dot(X0[:, st], X[:, st])
            if abs(o) < 0.9:
                print(f'  WARNING coord {k}: state {st+1} overlap {o:.3f}')
            if o < 0:
                X[:, st] *= -1.0
    dXj = (Xp[:, j - 1] - Xm[:, j - 1]) / (2 * DELTA)
    d_amp_true[k] = np.dot(X0[:, i - 1], dXj)
    theta = (Mp - Mm) / (2 * DELTA)
    dSao = (Sp_ - Sm_) / (2 * DELTA)
    Sket = C0py @ dSao @ C0py.T            # C0_f^T dS C0_f
    U = theta - Sket
    Sx = Sket + Sket.T                     # full S^[x] in MO
    con = np.max(np.abs(U + U.T + Sx))
    maxcon = max(maxcon, con)
    all_theta[k] = theta
    all_sket[k] = Sket

    d_orb_true[k] = np.sum(gam * theta)
    d_frozen_true[k] = np.sum(gam * Sket)
    # coded skeleton terms: same-space 1/2, cross lower->higher 1, x (-S^[x])
    w = np.zeros((nbf, nbf))
    same = SP[:, None] == SP[None, :]
    lohi = SP[:, None] < SP[None, :]
    w[same] = 0.5
    w[lohi] = 1.0
    d_skel_true[k] = -np.sum(gam * w * Sx)
    d_U_true[k] = np.sum(gam * U)
    # the part the z-vector is supposed to deliver:
    hilo = SP[:, None] > SP[None, :]
    gA = np.where(hilo, gam - gam.T, 0.0)  # [g(hi,lo)-g(lo,hi)] at (hi,lo)
    d_Ucphf_true[k] = np.sum(gA * U)

print('\nmax orthonormality-constraint residual |U+U^T+S^[x]| =', maxcon)

sgn = np.sign(np.dot(dn, da)) or 1.0
dn_al = sgn * dn


def rep(name, v, ref):
    c = np.dot(ref, v) / (np.linalg.norm(ref) * np.linalg.norm(v) + 1e-30)
    print(f'{name}: cos={c:+.6f}  |v|={np.linalg.norm(v):.5f}  '
          f'|ref|={np.linalg.norm(ref):.5f}  resid={np.linalg.norm(v - ref):.5f}')


print('\n===== completeness: d_CI + d_orb(semi-num) vs d_num =====')
rep('d_CI + d_orb_true vs d_num', dci + d_orb_true, dn_al)
print('d_orb_true   :', np.round(d_orb_true, 5))
print('d_num - d_CI :', np.round(dn_al - dci, 5))

print('\n===== TRUE amplitude term X_I . dX_J vs polarization d_CI =====')
print('d_amp_true:', np.round(d_amp_true, 5), ' |.|=%.5f' % np.linalg.norm(d_amp_true))
print('d_CI(pol) :', np.round(dci, 5), ' |.|=%.5f' % np.linalg.norm(dci))
rep('d_amp_true + d_orb_true vs d_num', d_amp_true + d_orb_true, dn_al)

print('\n===== term-by-term =====')
rep('frozen+skeleton coded vs semi-num', dov, d_frozen_true + d_skel_true)
print('  frozen semi-num  :', np.round(d_frozen_true, 5))
print('  skeleton semi-num:', np.round(d_skel_true, 5))
print('  coded d_ov       :', np.round(dov, 5))

print('\n===== the U-term and the code D =====')
print('d_U_true (all blocks)  :', np.round(d_U_true, 5))
print('d_Ucphf_true (antisym) :', np.round(d_Ucphf_true, 5))
# exact identity: d_U = d_Ucphf + d_skel_cross&same + same-space antisym resp.
samesp_resp = d_U_true - d_Ucphf_true - d_skel_true
print('same-space antisym U response (MISSING from coded assembly if big):')
print('  d_U - d_Ucphf - d_skel:', np.round(samesp_resp, 5),
      ' |.|=%.5f' % np.linalg.norm(samesp_resp))
if D_code is not None:
    rep('D_code vs d_Ucphf_true', D_code, d_Ucphf_true)
    nz = np.abs(d_Ucphf_true) > 1e-4
    ratios = D_code[nz] / d_Ucphf_true[nz]
    print('per-component D_code/d_Ucphf_true:', np.round(ratios, 5))
    s = np.dot(d_Ucphf_true, D_code) / np.dot(D_code, D_code)
    print(f'lstsq scale s (d_Ucphf ~ s*D_code): {s:+.7f}   1/s = {1/s:+.7f}')

tag = 'hf' if 'HF' in inp else 'bhh'
np.savez(f'/tmp/nactest/theta_data_{tag}.npz',
         theta=all_theta, sket=all_sket, tdF=tdF,
         dn=dn, dci=dci, dov=dov, D=(D_code if D_code is not None else np.zeros(ncoord)),
         bvec=bvraw0, damp=d_amp_true,
         energies=np.array(mol.energies), noca=noca, nocb=nocb, nbf=nbf,
         dn_al=dn_al)
print(f'\nsaved /tmp/nactest/theta_data_{tag}.npz')

print('\n===== reassembled total =====')
if D_code is not None:
    s = np.dot(d_Ucphf_true, D_code) / np.dot(D_code, D_code)
    total = dci + d_frozen_true + d_skel_true + s * D_code + samesp_resp
    rep('dci+frozen+skel+s*D+samespace vs num', total, dn_al)
