"""Pin the TLF state-overlap derivative decomposition.

Per displaced geometry (+/- delta per coordinate), call the production
state-overlap code (oqp.get_states_overlap) twice:
  (a) production-like: old = X(R0), new = X(R_disp), raw overlap_mo
      -> rebuild the numerical coupling, sanity-check vs numerical_nac;
  (b) FROZEN amplitudes: old = new = X(R0), overlap_mo phase-fixed to R0
      -> isolates the TLF orbital term with the code's own gamma weights.
Then test  d_num =?= d_amp (measured) + d_orb(TLF frozen).
"""
import sys
import oqp                       # import BEFORE numpy (ILP64 LAPACK)
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import NAC, SinglePoint
import numpy as np

inp = sys.argv[1] if len(sys.argv) > 1 else '/tmp/nactest/H2O_tight_dx0.001.inp'
DELTA = 1.0e-3
PAIR = (2, 3)
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=inp, log=inp.replace('.inp', '_tlf.log'))
r.run()
mol = r.mol
nac = NAC(mol)
nstate, natom = nac.nstate, nac.natom
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2

xyz0 = np.array(mol.get_system(), copy=True)
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = C0raw.shape[0]
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
C0raw_b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0_b = np.array(mol.data['OQP::E_MO_B'], copy=True)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
nij = noca * (nbf - nocb)
X0 = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()

# production numerical NAC (tight), the ground truth
nacv_n, dcv_n, fn = nac.numerical_nac()
i, j = PAIR
dn = dcv_n[i - 1, j - 1].reshape(-1)
dn_all = dcv_n.reshape((nstate, nstate, -1))

mol.save_data()
cfg = mol.config
json0 = mol.log.replace('.log', '.json')
cfg['guess']['type'] = 'json'
cfg['guess']['file'] = json0
cfg['guess']['continue_geom'] = False


def read_fmat(tag, n):
    raw = np.array(mol.data[tag], copy=True)
    return raw.reshape(-1).reshape((n, n)).T.copy()   # Fortran [row, col]


def write_fmat(tag, M_f):
    mol.data[tag] = np.ascontiguousarray(M_f.T)


def states_overlap():
    oqp.get_states_overlap(mol)
    # td_states_overlap Fortran (nstates, nstates): s_st(old_state, new_state)
    return read_fmat('OQP::td_states_overlap', nstate)


def displaced(coord):
    """SCF+TD at coord; returns (S_real, S_frozen, dXcol) vs R0."""
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
    M_f = read_fmat(OVTAG, nbf)            # raw <phi_p(R0)|phi_q(Rd)>
    Xd_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
    # Gauge-align the displaced MOs to R0 per orbital subspace (doc, socc,
    # virt) via Procrustes/polar decomposition of the diagonal blocks of M.
    # This removes BOTH sign flips and continuous rotations (the displaced
    # ROHF SCF can rotate the near-degenerate socc pair nondeterministically;
    # the production overlap is gauge-covariant, but the frozen-amplitude
    # decomposition is not unless the gauge is pinned).
    import os
    if os.environ.get('NAC_NOROT'):
        # production-equivalent gauge: per-orbital sign fixes only
        sgn0 = np.sign(np.diag(M_f)); sgn0[sgn0 == 0] = 1.0
        Q = np.diag(sgn0)
    else:
        Q = np.zeros((nbf, nbf))
        for lo, hi in ((0, nocb), (nocb, noca), (noca, nbf)):
            B = M_f[lo:hi, lo:hi]
            W, s, Vt = np.linalg.svd(B)
            Q[lo:hi, lo:hi] = Vt.T @ W.T   # Q_sp = V W^T -> B Q_sp = SPD
    M_al = M_f @ Q                         # aligned: diag(M_al) ~ +1
    if np.min(np.diag(M_al)) < 0.98:
        print('  WARNING: weak aligned diag', np.min(np.diag(M_al)))
    write_fmat(OVTAG, M_al)
    # rotate the displaced amplitudes into the aligned gauge:
    # X~ = Q_occ^T X Q_vb  (hole index in occ spaces, particle in beta-virt);
    # improper rotations (det -1) only flip the global state phase, which the
    # state alignment step absorbs.
    Xd_mat = Xd_raw.reshape(-1).reshape((nstate, nij))  # [st, ij]
    Xd_al = np.zeros_like(Xd_mat)
    nvirb = nbf - nocb
    for st in range(nstate):
        Xm = Xd_mat[st].reshape((nvirb, noca)).T        # X(i, a) i=occ, a=virtb
        Xt = Q[:noca, :noca].T @ Xm @ Q[nocb:, nocb:]
        Xd_al[st] = Xt.T.reshape(-1)
    Xd = Xd_al.T.copy()                                  # [ij, st]
    # state-sign-align the displaced amplitudes to R0 so the rebuilt state
    # overlap has +1 diagonal (clean central differencing)
    for st in range(nstate):
        if np.dot(X0[:, st], Xd[:, st]) < 0:
            Xd[:, st] *= -1.0
    # (a) production-like: old=X0, new=Xd (gauge- and sign-aligned), M aligned
    mol.data['OQP::td_bvec_mo_old'] = X0_raw
    mol.data['OQP::td_bvec_mo'] = np.ascontiguousarray(Xd.T).reshape(Xd_raw.shape)
    S_real = states_overlap()
    # (b) frozen amplitudes: old = new = X(R0), same aligned M
    mol.data['OQP::td_bvec_mo'] = X0_raw
    S_frozen = states_overlap()
    return S_real, S_frozen, Xd, M_al


ncoord = 3 * natom
A_real = np.zeros((ncoord, nstate, nstate))
A_frozen = np.zeros((ncoord, nstate, nstate))
d_amp = np.zeros((ncoord, nstate, nstate))
theta_al = np.zeros((ncoord, nbf, nbf))    # gauge-consistent with A_frozen

for k in range(ncoord):
    Sp_r, Sp_f, Xp, Mp = displaced(xyz0 + DELTA * np.eye(ncoord)[k])
    Sm_r, Sm_f, Xm, Mm = displaced(xyz0 - DELTA * np.eye(ncoord)[k])
    theta_al[k] = (Mp - Mm) / (2 * DELTA)
    print(f'k={k}: |S+f - I|={np.abs(Sp_f - np.eye(nstate)).max():.2e} '
          f'|S-f - I|={np.abs(Sm_f - np.eye(nstate)).max():.2e} '
          f'|S+r - I|={np.abs(Sp_r - np.eye(nstate)).max():.2e}', flush=True)
    for X in (Xp, Xm):                      # state phase alignment to R0
        for st in range(nstate):
            if np.dot(X0[:, st], X[:, st]) < 0:
                X[:, st] *= -1.0
    # production arithmetic: dcme = (S - S^T)/dt per displacement,
    # then (forward - backward)/2
    dc_p = (Sp_r - Sp_r.T) / DELTA
    dc_m = (Sm_r - Sm_r.T) / DELTA
    A_real[k] = (dc_p - dc_m) / 2.0
    dc_pf = (Sp_f - Sp_f.T) / DELTA
    dc_mf = (Sm_f - Sm_f.T) / DELTA
    A_frozen[k] = (dc_pf - dc_mf) / 2.0
    dX = (Xp - Xm) / (2 * DELTA)
    d_amp[k] = X0.T @ dX                    # d_amp[k][I,J] = X_I . dX_J


def vec(A, I, J):
    return A[:, I - 1, J - 1]


def rep(name, v, ref):
    c = np.dot(ref, v) / (np.linalg.norm(ref) * np.linalg.norm(v) + 1e-30)
    print(f'{name}: cos={c:+.6f}  |v|={np.linalg.norm(v):.5f}  '
          f'|ref|={np.linalg.norm(ref):.5f}  resid={np.linalg.norm(v - ref):.5f}')


print(f'\n===== pair {PAIR} =====')
print('d_num (production)    :', np.round(dn, 5), ' |.|=%.5f' % np.linalg.norm(dn))
dr = vec(A_real, i, j)
sgn_glob = np.sign(np.dot(dn, dr)) or 1.0
dn_al = sgn_glob * dn
rep('rebuilt (a) vs d_num  ', dr, dn_al)
df = vec(A_frozen, i, j)
da = vec(d_amp, i, j)
print('frozen-X orbital (b)  :', np.round(df, 5), ' |.|=%.5f' % np.linalg.norm(df))
print('amplitude term        :', np.round(da, 5), ' |.|=%.5f' % np.linalg.norm(da))
rep('amp + frozen vs d_num ', da + df, dn_al)
rep('amp + frozen vs rebuilt', da + df, dr)

tag = 'hf' if 'HF' in inp.upper() and 'BHH' not in inp.upper() else 'bhh'
if 'hf' in inp.lower() and 'tight' in inp.lower():
    tag = 'hf'
np.savez(f'/tmp/nactest/tlf_deriv_{tag}.npz', A_real=A_real, A_frozen=A_frozen,
         d_amp=d_amp, dn=dn, dn_all=dn_all, X0=X0, noca=noca, nocb=nocb, nbf=nbf,
         theta=theta_al)
print(f'\nsaved /tmp/nactest/tlf_deriv_{tag}.npz')

print('\n===== the exact relation, all pairs: dn =?= -(2*damp + A_frozen) =====')
for (I, J) in ((1, 2), (1, 3), (2, 3)):
    t = dn_all[I - 1, J - 1]
    v = -(2 * d_amp[:, I - 1, J - 1] + A_frozen[:, I - 1, J - 1])
    s = np.sign(np.dot(t, v)) or 1.0
    c = np.dot(t, s * v) / (np.linalg.norm(t) * np.linalg.norm(v) + 1e-30)
    print(f'({I},{J}): cos={c:+.6f}  |pred|={np.linalg.norm(v):.5f}  '
          f'|dn|={np.linalg.norm(t):.5f}  resid={np.linalg.norm(s*v - t):.5f}')

print('\n===== in-process rebuild vs production, all pairs =====')
for (I, J) in ((1, 2), (1, 3), (2, 3)):
    t = dn_all[I - 1, J - 1]
    v = A_real[:, I - 1, J - 1]
    s = np.sign(np.dot(t, v)) or 1.0
    c = np.dot(t, s * v) / (np.linalg.norm(t) * np.linalg.norm(v) + 1e-30)
    print(f'({I},{J}): cos={c:+.6f}  |rebuild|={np.linalg.norm(v):.5f}  '
          f'|dn|={np.linalg.norm(t):.5f}  resid={np.linalg.norm(s*v - t):.5f}')
    # and the split consistency: rebuild =?= -(...)? sign conventions equal:
    w = 2 * d_amp[:, I - 1, J - 1] + A_frozen[:, I - 1, J - 1]
    sw = np.sign(np.dot(v, w)) or 1.0
    print(f'        split check |rebuild -+ (2damp+Afrozen)| = '
          f'{min(np.linalg.norm(v - w), np.linalg.norm(v + w)):.5f}')
