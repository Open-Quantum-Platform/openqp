"""Full analytical-NAC benchmark against the production numerical reference.

Assembly (Phase 8 recipe):
    d_pred = -[ 2*d_amp + d_orb ]
    d_orb  = (frozen+skeleton)(gamma^TLF)  -  D_cphf(gamma^TLF)
with gamma^TLF built in closed form (validated to 1e-12 against the
overlap-code oracle), pushed to the verified Fortran machinery via the
OQP::nac_gamma_tlf tag.  The amplitude term d_amp = X_I . dX_J is measured
semi-numerically IN THIS PROCESS (phase-coherent with everything else);
its analytic replacement (bilinear interstate z-vector) is the one
remaining implementation item.
"""
import sys
import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import NAC, SinglePoint
import numpy as np

inp = sys.argv[1] if len(sys.argv) > 1 else '/tmp/nactest/H2O_tight_dx0.001.inp'
DELTA = 1.0e-3
RS = 1.0 / np.sqrt(2.0)
SGS = (1.0, -1.0, -1.0)          # (s_ij, s_ab, s_ia) from the closed-form fit
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=inp, log=inp.replace('.inp', '_bench.log'))
r.run()
mol = r.mol
nac = NAC(mol)
nstate, natom = nac.nstate, nac.natom
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = C0raw.shape[0]
nvirb = nbf - nocb
nij = noca * nvirb
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
C0raw_b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0_b = np.array(mol.data['OQP::E_MO_B'], copy=True)
xyz0 = np.array(mol.get_system(), copy=True)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
X0 = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()


# ===================== gamma^TLF closed form =====================
def unfold(bv, st):
    ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
    ijlr2 = (noca - nocb - 1) * noca + noca
    x = np.zeros((noca, nvirb))
    for i in range(1, noca + 1):
        for jj in range(nocb + 1, nbf + 1):
            ij = (jj - nocb - 1) * noca + i
            if ij == ijlr1:
                x[i - 1, jj - nocb - 1] = bv[ijlr1 - 1, st - 1] * RS
            elif ij == ijlr2:
                x[i - 1, jj - nocb - 1] = -bv[ijlr1 - 1, st - 1] * RS
            else:
                x[i - 1, jj - nocb - 1] = bv[ij - 1, st - 1]
    return x


iocc = np.arange(noca)
avir = np.arange(nvirb)
OO = (iocc[:, None] >= nocb) & (avir[None, :] < 2)
GEN = ~OO
S_IA0 = np.zeros((noca, nvirb))
for a in range(2):
    S_IA0[nocb + a, a] = 1.0


def kernel_pair(cI, cJ):
    sg_ij, sg_ab, sg_ia = SGS
    G = np.zeros((nbf, nbf))
    co = cI * GEN
    cn = cJ * GEN
    pq = slice(nocb, noca)
    qv = slice(0, 2)
    # term A
    G[nocb:, nocb:] += sg_ab * (co.T @ cn)
    G[:noca, :noca] += sg_ij * (cn @ co.T)
    # term B
    g0 = co @ S_IA0.T
    d0 = (cn @ S_IA0.T).T
    G[:noca, nocb:] += sg_ia * (co.T @ d0).T
    G[:noca, nocb:] += sg_ia * (g0 @ cn)
    # term C (OO x OO)
    cIo = cI[pq, qv]
    cJo = cJ[pq, qv]
    G[nocb:noca, nocb:noca] += sg_ij * (cJo @ cIo.T)
    G[nocb:noca, nocb:noca] += sg_ab * (cIo.T @ cJo)
    # terms D/E (OO x generic, 1/sqrt2)
    cIo_full = np.zeros_like(cI); cIo_full[pq, qv] = cI[pq, qv]
    cJo_full = np.zeros_like(cJ); cJo_full[pq, qv] = cJ[pq, qv]
    cJg = cJ * GEN
    cIg = cI * GEN
    G[:noca, nocb:noca] += RS * sg_ij * (cJg[:, 0:2] @ cIo_full[pq, qv].T)
    G[nocb:nocb + 2, nocb:] += RS * sg_ab * (cIo_full[pq, qv].T @ cJg[pq, :])
    G[nocb:noca, nocb:] += RS * sg_ia * (cIo_full[pq, qv] @ cJg[pq, :])
    G[:noca, nocb:nocb + 2] += RS * sg_ia * (cJg[:, 0:2] @ cIo_full[pq, qv])
    G[nocb:noca, :noca] += RS * sg_ij * (cJo_full[pq, qv] @ cIg[:, 0:2].T)
    G[nocb:, nocb:nocb + 2] += RS * sg_ab * (cIg[pq, :].T @ cJo_full[pq, qv])
    G[:noca, nocb:nocb + 2] += RS * sg_ia * (cIg[:, 0:2] @ cJo_full[pq, qv])
    G[nocb:noca, nocb:] += RS * sg_ia * (cJo_full[pq, qv] @ cIg[pq, :])
    return G


C = [unfold(X0, s + 1) for s in range(nstate)]
sp = np.zeros(nbf, dtype=int); sp[nocb:noca] = 1; sp[noca:] = 2
cross = sp[:, None] != sp[None, :]

gam_tlf = np.zeros((nbf, nbf, nstate, nstate))
for I in range(1, nstate + 1):
    for J in range(1, nstate + 1):
        if I == J:
            continue
        G = kernel_pair(C[I - 1], C[J - 1]) - kernel_pair(C[J - 1], C[I - 1])
        gam_tlf[:, :, I - 1, J - 1] = np.where(cross, G, 0.0)  # same-space = 0
                                                               # (pure gauge)

# push to Fortran: F-flat[k + nbf^2*(ist) + nbf^2*nst*(jst)], k = p + q*nbf
flatF = np.zeros(nbf * nbf * nstate * nstate)
for I in range(nstate):
    for J in range(nstate):
        flatF[(I + J * nstate) * nbf * nbf:(I + J * nstate + 1) * nbf * nbf] = \
            gam_tlf[:, :, I, J].T.reshape(-1)   # C-flat of transpose = F-flat
mol.data['OQP::nac_gamma_tlf'] = flatF.reshape((nbf * nbf, nstate, nstate))

# ===================== production numerical reference =====================
nacv_n, dcv_n, fn = nac.numerical_nac()
dn_all = dcv_n.reshape((nstate, nstate, -1))

# ===================== analytic orbital pieces =====================
oqp.mrsf_nac_overlap(mol)                  # frozen+skeleton with gamma^TLF
ovraw = np.array(mol.data['OQP::nac_overlap'], copy=True)
ovF = ovraw.reshape(-1).reshape((nstate, nstate, 3 * natom))
dov = np.transpose(ovF, (1, 0, 2))         # dov[I-1, J-1, :]


def grad_now():
    oqp.tdhf_mrsf_gradient(mol)
    return mol.get_grad().reshape((natom, 3)).copy()


def cphf_term(i, j):
    mol.data.set_tdhf_target(i)
    oqp.set_mrsf_nac_cphf(mol, i, j)
    oqp.tdhf_mrsf_z_vector(mol)
    if not mol.mol_energy.Z_Vector_converged:
        oqp.set_mrsf_nac_cphf(mol, 0, 0)
        return None
    gZ = grad_now()
    mol.data['OQP::td_p'] = np.zeros_like(np.array(mol.data['OQP::td_p'], copy=True))
    mol.data['OQP::WAO'] = np.zeros_like(np.array(mol.data['OQP::WAO'], copy=True))
    gS = grad_now()
    oqp.set_mrsf_nac_cphf(mol, 0, 0)
    return (gZ - gS).reshape(-1)


D = {}
for (i, j) in ((1, 2), (1, 3), (2, 3)):
    D[(i, j)] = cphf_term(i, j)

# ===================== amplitude term (semi-numerical, in-process) ========
mol.save_data()
cfg = mol.config
json0 = mol.log.replace('.log', '.json')
cfg['guess']['type'] = 'json'
cfg['guess']['file'] = json0
cfg['guess']['continue_geom'] = False


def displaced_amps(coord):
    mol.update_system(coord)
    oqp.library.ints_1e(mol)
    oqp.library.guess(mol)
    SinglePoint(mol).energy()
    mol.data['OQP::xyz_old'] = xyz0.reshape((3, natom))
    mol.data['OQP::VEC_MO_A_old'] = C0raw
    mol.data['OQP::E_MO_A_old'] = e0
    mol.data['OQP::VEC_MO_B_old'] = C0raw_b
    mol.data['OQP::E_MO_B_old'] = e0_b
    oqp.get_structures_ao_overlap(mol)
    M_f = np.array(mol.data[OVTAG], copy=True).reshape(-1).reshape((nbf, nbf)).T
    # Procrustes parallel-transport gauge (consistent with the cross-only
    # analytic orbital term). The rotation must act on the UNFOLDED
    # determinant amplitudes and be refolded: the spin-adapted ijlr slots do
    # NOT transform like grid entries (rotating the folded grid was the
    # source of the former '(1,2) 25% gap'). The singlet sector is closed
    # under socc rotations (the triplet ijlr2 combo is trace-like and stays
    # empty), so the refold is exact.
    Q = np.zeros((nbf, nbf))
    for lo, hi in ((0, nocb), (nocb, noca), (noca, nbf)):
        W, s, Vt = np.linalg.svd(M_f[lo:hi, lo:hi])
        Q[lo:hi, lo:hi] = Vt.T @ W.T
    Xd_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
    Xd_mat = Xd_raw.reshape(-1).reshape((nstate, nij))
    Xd_al = np.zeros_like(Xd_mat)
    for st in range(nstate):
        c = unfold(Xd_mat.T, st + 1)               # determinant amplitudes
        cp = Q[:noca, :noca].T @ c @ Q[nocb:, nocb:]
        leak = abs(cp[nocb, 0] + cp[nocb + 1, 1])  # triplet (trace) leakage
        if leak > 1e-8:
            print(f'  WARNING: ijlr2 leakage {leak:.2e} state {st+1}')
        f = cp.copy()
        f[nocb, 0] = np.sqrt(2.0) * cp[nocb, 0]    # refold singlet slot
        f[nocb + 1, 1] = 0.0                       # ijlr2 slot stays empty
        Xd_al[st] = f.T.reshape(-1)
    return Xd_al.T


ncoord = 3 * natom
d_amp = np.zeros((ncoord, nstate, nstate))
for k in range(ncoord):
    Xp = displaced_amps(xyz0 + DELTA * np.eye(ncoord)[k])
    Xm = displaced_amps(xyz0 - DELTA * np.eye(ncoord)[k])
    for X in (Xp, Xm):
        for st in range(nstate):
            if np.dot(X0[:, st], X[:, st]) < 0:
                X[:, st] *= -1.0
    d_amp[k] = X0.T @ ((Xp - Xm) / (2 * DELTA))

# ===================== assembly and report =====================
print('\n================== FULL ANALYTICAL NAC BENCHMARK ==================')
print(f'input: {inp}')
res = {}
for (i, j) in ((1, 2), (1, 3), (2, 3)):
    dn = dn_all[i - 1, j - 1]
    dov_ij = dov[i - 1, j - 1]
    Dij = D[(i, j)]
    damp = d_amp[:, i - 1, j - 1]
    # sum gamma^TLF theta (cross blocks) = D - dov in the Fortran-contraction
    # sign conventions (one global sign vs the kernel extraction, verified
    # against the exact requirement -dn - 2*damp at cos +1.0000)
    d_orb = Dij - dov_ij
    # factor-2 fix (2026-06-18): numerical_nac() now uses the standard
    # (S - S^T)/(2*dt) HST convention, so dn is the physical <I|dJ> (half the
    # former value). The master decomposition dn_old = -(2*damp + d_orb) was an
    # identity against the 2x-too-large numerical, so the physical prediction is
    # d_phys = -(2*damp + d_orb)/2 = -(damp + d_orb/2). Halve to match.
    pred = -(2 * damp + d_orb) * 0.5
    s = np.sign(np.dot(dn, pred)) or 1.0
    cos = np.dot(dn, s * pred) / (np.linalg.norm(dn) * np.linalg.norm(pred) + 1e-30)
    print(f'\npair ({i},{j}):  gap = {mol.energies[j] - mol.energies[i]:.6f}')
    print(f'  |d_num| = {np.linalg.norm(dn):.5f}')
    print(f'  components: |2*amp| = {np.linalg.norm(2*damp):.5f}  '
          f'|orb(frozen+skel)| = {np.linalg.norm(dov_ij):.5f}  '
          f'|orb(CPHF)| = {np.linalg.norm(Dij):.5f}')
    print(f'  PRED vs NUM: cos = {cos:+.6f}  |pred| = {np.linalg.norm(pred):.5f}  '
          f'resid = {np.linalg.norm(s*pred - dn):.5f}  '
          f'({100*np.linalg.norm(s*pred - dn)/np.linalg.norm(dn):.2f}%)')
    print('  pred:', np.round(s * pred, 5))
    print('  num :', np.round(dn, 5))
    res[(i, j)] = dict(dn=dn, pred=s * pred, cos=cos, damp=damp,
                       dov=dov_ij, D=Dij)

tag = 'hf' if 'HF' in inp else 'bhh'
np.savez(f'/tmp/nactest/benchmark_{tag}.npz',
         **{f'{k[0]}{k[1]}_{n}': v for k, d_ in res.items() for n, v in d_.items()})
print(f'\nsaved /tmp/nactest/benchmark_{tag}.npz')
