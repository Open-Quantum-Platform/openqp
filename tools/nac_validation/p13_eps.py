"""fa0 looks DIAGONAL = orbital energies. If so, the esum orbital gradient is
trivial and the whole R_esum derivation simplifies.

Check 1: is fa0 == diag(E_MO_A)?  (up to the ixcore core shift on MO 0)
Check 2: THE decisive test.  If fa = C^T F_AO C is REBUILT from the rotated MOs
         in the matvec, then rotating C must change fa off-diagonally.  But if the
         matvec instead just re-reads a stored MO-basis diagonal (orbital
         energies) and does NOT rebuild it from an AO Fock, then rotating C leaves
         fa UNCHANGED.  That determines the entire form of dfa/dU.
         Measure d fa/dU_pq directly and see whether it is (a) ~0 (fa frozen),
         (b) delta_sq eps_p + delta_rq eps_p type (diagonal transported), or
         (c) the full delta_sq fa_rp + delta_rq fa_ps.
"""
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
EPS = 1.0e-4
r = Runner(input_file=INP,
           log='/bighome/alireza/openqp-nac/tools/nac_validation/p13_eps.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
C0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = C0.shape[0]
nij = noca * (nbf - nocb)
eMO = np.array(mol.data['OQP::E_MO_A'], copy=True).ravel()
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((-1, nij))


def fa_now():
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = X0[0]
    mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)
    oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_fa'], copy=True).reshape(-1).reshape((nbf, nbf)).T


fa0 = fa_now()
print(f'# E_MO_A            = {np.array2string(eMO, precision=3, max_line_width=200)}')
print(f'# diag(fa0)         = {np.array2string(np.diag(fa0), precision=3, max_line_width=200)}')
print(f'# |fa0 - diag(fa0)| off-diagonal max = {np.abs(fa0 - np.diag(np.diag(fa0))).max():.3e}')
d = np.diag(fa0).copy()
# MO 0 gets the -1e6 ixcore shift; compare the rest to eMO
print(f'# |diag(fa0)[1:] - E_MO_A[1:]| max = {np.abs(d[1:] - eMO[1:]).max():.3e}')
print(f'# fa0[0,0] = {fa0[0,0]:.4f}  (E_MO_A[0] = {eMO[0]:.4f}; ixcore shift expected)')

# ---- Check 2: does rotating C change fa? ----
print('\n# d fa/dU_pq structure (rotate C[:,q]+=eps*C[:,p], Fock frozen):')
for (p, q) in [(2, 8), (7, 3), (5, 5), (1, 12)]:
    Ca = C0.copy()
    Ca[:, q] += EPS * C0[:, p]
    mol.data['OQP::VEC_MO_A'] = Ca
    fp = fa_now()
    Ca = C0.copy()
    Ca[:, q] -= EPS * C0[:, p]
    mol.data['OQP::VEC_MO_A'] = Ca
    fm = fa_now()
    mol.data['OQP::VEC_MO_A'] = C0
    dfa = (fp - fm) / (2 * EPS)
    # prediction if fa is truly rebuilt as C^T F_AO C with fa=diag(eps):
    #   d fa_rs/dU_pq = delta_sq fa_rp + delta_rq fa_ps
    #   with fa diagonal, fa_rp = eps_p delta_rp, so this = delta_sq delta_rp eps_p
    #                                                      + delta_rq delta_ps eps_p
    pred = np.zeros((nbf, nbf))
    pred[p, q] += d[p]         # r=p,s=q term
    pred[q, p] += d[p]         # r=q,s=p term
    c = np.sum(dfa * pred) / (np.linalg.norm(dfa) * np.linalg.norm(pred) + 1e-300)
    nz = np.argwhere(np.abs(dfa) > 1e-3)
    print(f'  (p,q)=({p},{q}): |dfa|={np.linalg.norm(dfa):.4f}, '
          f'nonzero entries={len(nz)}, cos(dfa,pred_diag)={c:+.5f}')
    for (i, j) in nz[:6]:
        print(f'      dfa[{i},{j}]={dfa[i,j]:10.4f}  '
              f'(eps_{i}={d[i]:.3f} eps_{j}={d[j]:.3f})')
