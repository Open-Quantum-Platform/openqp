"""Extract gamma^TLF exactly by probing the production overlap code.

The TLF state overlap is bilinear in (X_old, X_new) and, to first order,
linear in (M - I).  With frozen amplitudes X_old = X_new = X(R0) and
synthetic MO overlaps M = I +/- eps*E_pq, the central difference

    gamma~_pq[I,J] = d(S_IJ - S_JI)/dM_pq

is exactly the TLF orbital-coupling kernel: A_frozen_IJ = sum_pq
gamma~_pq theta^x_pq.  361 oracle calls, no SCF, no derivation guesswork.

Validates against A_frozen / theta saved by tlf_deriv_check.py /
theta_check.py, and compares gamma~ to the analytic-side gamma_code
(get_mrsf_transition_density) block by block.
"""
import sys
import oqp
import oqp.library
from oqp.pyoqp import Runner
import numpy as np

inp = sys.argv[1] if len(sys.argv) > 1 else '/tmp/nactest/H2O_tight_dx0.001.inp'
EPS = 1.0e-6
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=inp, log=inp.replace('.inp', '_oracle.log'))
r.run()
mol = r.mol
nstate = mol.config['tdhf']['nstate']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = C0raw.shape[0]
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)

# the overlap tag must exist before we can write it: create it via one
# real call of get_structures_ao_overlap against R0 itself
xyz0 = np.array(mol.get_system(), copy=True)
mol.data['OQP::xyz_old'] = xyz0.reshape((3, len(xyz0) // 3))
mol.data['OQP::VEC_MO_A_old'] = C0raw
mol.data['OQP::E_MO_A_old'] = np.array(mol.data['OQP::E_MO_A'], copy=True)
mol.data['OQP::VEC_MO_B_old'] = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
mol.data['OQP::E_MO_B_old'] = np.array(mol.data['OQP::E_MO_B'], copy=True)
oqp.get_structures_ao_overlap(mol)
# use the real R0-R0 MO overlap (identity up to ~1e-12 fuzz) as the probe
# base: an EXACT identity makes comp_det divide by exactly-zero pivots on
# exactly-singular minors -> NaN.
Mbase = np.array(mol.data[OVTAG], copy=True).reshape(-1).reshape((nbf, nbf)).T

mol.data['OQP::td_bvec_mo_old'] = X0_raw
mol.data['OQP::td_bvec_mo'] = X0_raw


def write_fmat(tag, M_f):
    mol.data[tag] = np.ascontiguousarray(M_f.T)


def s_overlap(M_f):
    write_fmat(OVTAG, M_f)
    oqp.get_states_overlap(mol)
    raw = np.array(mol.data['OQP::td_states_overlap'], copy=True)
    return raw.reshape(-1).reshape((nstate, nstate)).T.copy()  # [old, new]


# sanity: base M -> S = I
S0 = s_overlap(Mbase)
print('sanity |S(Mbase) - I|max =', np.abs(S0 - np.eye(nstate)).max())

gam = np.zeros((nbf, nbf, nstate, nstate))
for p in range(nbf):
    for q in range(nbf):
        E = np.zeros((nbf, nbf))
        E[p, q] = EPS
        Sp = s_overlap(Mbase + E)
        Sm = s_overlap(Mbase - E)
        D = (Sp - Sp.T - Sm + Sm.T) / (2 * EPS)
        gam[p, q] = D
print('oracle probes done')

# export gamma_code from THIS process (phase-consistent with gam)
oqp.mrsf_nac_overlap(mol)
tdraw0 = np.array(mol.data['OQP::nac_trden_mo'], copy=True)
tdF0 = tdraw0.reshape(-1).reshape((nstate, nstate, nbf * nbf))
np.savez('/tmp/nactest/tlf_gamma.npz', gam=gam, tdF=tdF0, X0=X0_raw,
         noca=noca, nocb=nocb, nbf=nbf)

# ===== validation: sum gamma~ theta vs A_frozen =====
tag = 'hf' if 'HF' in inp else 'bhh'
try:
    dd = np.load(f'/tmp/nactest/tlf_deriv_{tag}.npz')
    th = np.load(f'/tmp/nactest/theta_data_{tag}.npz')
    theta = th['theta']                      # (9, nbf, nbf)
    for (I, J) in ((1, 2), (1, 3), (2, 3)):
        pred = np.einsum('pqk->k', np.zeros((1, 1, 1)))  # placeholder
        pred = np.array([np.sum(gam[:, :, I - 1, J - 1] * theta[k])
                         for k in range(theta.shape[0])])
        ref = dd['A_frozen'][:, I - 1, J - 1]
        c = np.dot(pred, ref) / (np.linalg.norm(pred) * np.linalg.norm(ref) + 1e-30)
        print(f'({I},{J}): cos(pred,A_frozen)={c:+.6f}  |pred|={np.linalg.norm(pred):.5f}'
              f'  |A_frozen|={np.linalg.norm(ref):.5f}  '
              f'resid={np.linalg.norm(pred - ref):.5f}')
except FileNotFoundError as e:
    print('validation data missing:', e)

# ===== structure: gamma~ vs gamma_code =====
# gamma_code from the analytic side (exported by mrsf_nac_overlap)
oqp.mrsf_nac_overlap(mol)
tdraw = np.array(mol.data['OQP::nac_trden_mo'], copy=True)
tdF = tdraw.reshape(-1).reshape((nstate, nstate, nbf * nbf))


def gamma_code(I, J):
    return tdF[J - 1, I - 1].reshape((nbf, nbf)).T.copy()


sp = np.zeros(nbf, dtype=int)
sp[nocb:noca] = 1
sp[noca:] = 2
names = {0: 'd', 1: 's', 2: 'v'}
print('\nblock-by-block gamma~(2,3) vs gamma_code(2,3):')
g_t = gam[:, :, 1, 2]
g_c = gamma_code(2, 3)
for bi in range(3):
    for bj in range(3):
        m = (sp[:, None] == bi) & (sp[None, :] == bj)
        nt, nc = np.linalg.norm(g_t[m]), np.linalg.norm(g_c[m])
        if nt < 1e-10 and nc < 1e-10:
            continue
        ratio = nt / nc if nc > 1e-10 else np.inf
        cc = (np.dot(g_t[m], g_c[m]) / (nt * nc + 1e-30)) if nc > 1e-10 else 0
        print(f'  {names[bi]}{names[bj]}: |g~|={nt:.5f} |g_code|={nc:.5f} '
              f'ratio={ratio:.5f} cos={cc:+.4f}')
print('\nsaved /tmp/nactest/tlf_gamma.npz')
