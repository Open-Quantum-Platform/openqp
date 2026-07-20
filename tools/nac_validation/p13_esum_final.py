"""
Phase 13 step 3: the CORRECT closed form for R_esum, using the fact that at the
reference canonical orbitals fa = diag(eps).

E_esum = sum_rs Gam_rs fa_rs   (Gam = Gam_A on occ-occ + Gam_B on virt-virt).
With fa = C^T F_AO C rebuilt from the rotated C, the frozen-Fock derivative is
    d fa_rs / dU_pq = delta_sq fa_rp + delta_rq fa_ps.
So
    R_pq = dE_esum/dU_pq = sum_rs Gam_rs (delta_sq fa_rp + delta_rq fa_ps)
         = sum_r Gam_rq fa_rp + sum_s Gam_qs fa_ps
         = (Gam^T fa)_qp + (Gam fa)_qp            [matrix form, fa the FULL fa]
    R    = fa Gam + fa Gam^T   (as an (nbf x nbf) matrix, R_pq)   ... (form 1)
or equivalently R_pq = [(Gam+Gam^T) fa]_qp = [(fa)(Gam+Gam^T)]... check both.

The earlier Gate-2 failure used the DIAGONAL fa (eps) in the closed form while the
MEASUREMENT saw the FULL dense d fa/dU.  Here fa0 IS diagonal at the reference, so
that cannot be the difference -- meaning the real bug was the matrix ORDER /
transpose in assembling R_pq from the two delta contractions.  This script tests
all the index orderings against the direct Fock-scaling measurement so the right
one is settled numerically, once, at the reference geometry.
"""
import math
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
ETA = 1.0e-4
EPS = 1.0e-4
r = Runner(input_file=INP,
           log='/bighome/alireza/openqp-nac/tools/nac_validation/p13_final.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate = mol.config['tdhf']['nstate']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
C0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
nbf = C0.shape[0]
nvirb = nbf - nocb
nij = noca * nvirb
F0a = np.array(mol.data['OQP::FOCK_A'], copy=True)
F0b = np.array(mol.data['OQP::FOCK_B'], copy=True)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
SQ = 1.0 / math.sqrt(2.0)
ijlr1 = (noca - 1 - nocb - 1) * noca + (noca - 1) - 1
ijlr2 = (noca - nocb - 1) * noca + (noca) - 1
diag = [(I, I) for I in range(nstate)]


def set_bvec(col):
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)


def put_C(Ca, Cb):
    mol.data['OQP::VEC_MO_A'] = Ca
    mol.data['OQP::VEC_MO_B'] = Cb


def put_F(s):
    mol.data['OQP::FOCK_A'] = F0a * s
    mol.data['OQP::FOCK_B'] = F0b * s


def E_and_f(I, J):
    set_bvec(X0[J])
    oqp.mrsf_matvec_apply(mol)
    ax = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
    fa = np.array(mol.data['OQP::nac_fa'], copy=True).reshape(-1).reshape((nbf, nbf)).T
    fb = np.array(mol.data['OQP::nac_fb'], copy=True).reshape(-1).reshape((nbf, nbf)).T
    return float(X0[I] @ ax), fa, fb


def E_esum(I, J):
    put_F(1.0)
    e0, fa, fb = E_and_f(I, J)
    put_F(1.0 + ETA)
    e1, _, _ = E_and_f(I, J)
    put_F(1.0)
    return (e1 - e0) / ETA, fa, fb


def unfold(col):
    x = np.zeros((noca, nvirb))
    for i in range(noca):
        for a in range(nvirb):
            ij = a * noca + i
            if ij == ijlr1:
                x[i, a] = col[ijlr1] * SQ
            elif ij == ijlr2:
                x[i, a] = -col[ijlr1] * SQ
            else:
                x[i, a] = col[ij]
    return x


Xt = [unfold(X0[k]) for k in range(nstate)]
put_F(1.0)
_, fa0, fb0 = E_and_f(0, 0)


def GamAB(I, J):
    GA = np.zeros((nbf, nbf))
    GB = np.zeros((nbf, nbf))
    GA[0:noca, 0:noca] = -0.5 * (Xt[I] @ Xt[J].T + Xt[J] @ Xt[I].T)
    GB[nocb:nbf, nocb:nbf] = 0.5 * (Xt[I].T @ Xt[J] + Xt[J].T @ Xt[I])
    return GA, GB


# candidate closed forms for R_pq (Gam symmetric so several coincide; list anyway)
def candidates(GA, GB):
    def one(G, f):
        return {
            'f@G + f@G^T': f @ G + f @ G.T,
            'G@f + G^T@f': G @ f + G.T @ f,
            '(f@G + f@G^T)^T': (f @ G + f @ G.T).T,
            '(G@f)+(f@G)': G @ f + f @ G,
        }
    ca = one(GA, fa0)
    cb = one(GB, fb0)
    return {k: ca[k] + cb[k] for k in ca}


print('# CHECK A: d fa/dU_pq (measured) vs delta_sq fa_rp + delta_rq fa_ps, FULL fa')
for (p, q) in [(2, 8), (7, 3), (5, 5)]:
    Ca = C0.copy(); Ca[:, q] += EPS * C0[:, p]; put_C(Ca, C0b)
    put_F(1.0); _, fp, _ = E_and_f(0, 0)
    Ca = C0.copy(); Ca[:, q] -= EPS * C0[:, p]; put_C(Ca, C0b)
    put_F(1.0); _, fm, _ = E_and_f(0, 0)
    put_C(C0, C0b); put_F(1.0)
    dfa = (fp - fm) / (2 * EPS)
    pred = np.zeros((nbf, nbf))
    pred[:, q] += fa0[:, p]
    pred[q, :] += fa0[p, :]
    c = np.sum(dfa * pred) / (np.linalg.norm(dfa) * np.linalg.norm(pred) + 1e-300)
    print(f'  (p,q)=({p},{q}): cos={c:+.6f} ratio={np.linalg.norm(pred)/np.linalg.norm(dfa):.4f}')

print('\n# CHECK C: full R_esum (measured by Fock scaling) vs each candidate form')
names = None
for (I, J) in diag:
    Rm = np.zeros((nbf, nbf))
    for p in range(nbf):
        for q in range(nbf):
            for sgn in (+1, -1):
                Ca, Cb = C0.copy(), C0b.copy()
                Ca[:, q] += sgn * EPS * C0[:, p]
                Cb[:, q] += sgn * EPS * C0b[:, p]
                put_C(Ca, Cb)
                e, _, _ = E_esum(I, J)
                Rm[p, q] += sgn * e / (2.0 * EPS)
    put_C(C0, C0b); put_F(1.0)
    GA, GB = GamAB(I, J)
    cands = candidates(GA, GB)
    if names is None:
        names = list(cands)
        print('  state | ' + ' | '.join(f'{n:>18}' for n in names))
    row = []
    for n in names:
        c = np.sum(Rm * cands[n]) / (np.linalg.norm(Rm) * np.linalg.norm(cands[n]) + 1e-300)
        rr = np.linalg.norm(cands[n]) / (np.linalg.norm(Rm) + 1e-300)
        row.append(f'{c:+.4f}/{rr:.3f}')
    print(f'   {I+1:>4} | ' + ' | '.join(f'{x:>18}' for x in row))
