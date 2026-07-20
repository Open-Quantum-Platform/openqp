"""
Nail down d fa / dU with the cleanest possible probe, no amplitudes, no Gam.

fa = mo_a^T F_AO mo_a  (confirmed: orthogonal_transform_sym, no renormalization).
Apply an EXACT unitary rotation in the (p,q) plane to mo_a:
    C' = C R,  R = Givens(theta) in columns p,q
Then fa' = R^T fa R exactly, and
    d fa/dtheta = [K^T fa + fa K],  K = dR/dtheta = E_pq - E_qp (antisymmetric).
So d fa/dtheta should equal  fa @ (E_pq - E_qp) + (E_pq - E_qp)^T @ fa
                           = fa @ K - K @ fa   (K antisym, K^T = -K).
Test:
  (1) does the exported nac_fa transform as R^T fa R under a finite rotation?
      -> confirms the tag write reaches the matvec AND the fa read layout is right.
  (2) does d fa/dtheta match fa@K - K@fa ?
      -> confirms the analytic derivative.
If (1) holds but (2) fails, the read layout is transposed; if (1) fails, the
VEC_MO_A tag write is not consumed (need a different injection).
"""
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
TH = 1.0e-4
r = Runner(input_file=INP,
           log='/bighome/alireza/openqp-nac/tools/nac_validation/p13_fadir.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
C0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = C0.shape[0]
nij = noca * (nbf - nocb)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((-1, nij))
rr = X0_raw.copy().reshape(-1)
rr[0:nij] = X0[0]
mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)


def fa_at(C, layout='T'):
    mol.data['OQP::VEC_MO_A'] = C
    oqp.mrsf_matvec_apply(mol)
    raw = np.array(mol.data['OQP::nac_fa'], copy=True).reshape(-1).reshape((nbf, nbf))
    return raw.T if layout == 'T' else raw


fa0 = fa_at(C0)


def givens(p, q, th):
    R = np.eye(nbf)
    c, s = np.cos(th), np.sin(th)
    R[p, p] = c; R[q, q] = c; R[p, q] = -s; R[q, p] = s
    return R


for (p, q) in [(2, 8), (7, 3), (1, 12)]:
    R = givens(p, q, TH)
    faR = fa_at(C0 @ R)
    # (1) finite-rotation transform check, both read layouts
    predT = R.T @ fa0 @ R
    faR_N = fa_at(C0 @ R, layout='N')
    e_T = np.abs(faR - predT).max()
    e_N = np.abs(faR_N - predT).max()
    # (2) analytic derivative check
    Rm = givens(p, q, -TH)
    faM = fa_at(C0 @ Rm)
    dfa = (faR - faM) / (2 * TH)
    K = np.zeros((nbf, nbf)); K[p, q] = 1.0; K[q, p] = -1.0
    pred = fa0 @ K - K @ fa0
    c = np.sum(dfa * pred) / (np.linalg.norm(dfa) * np.linalg.norm(pred) + 1e-300)
    print(f'(p,q)=({p},{q}): transform |faR-R^T fa0 R| ={e_T:.2e} (layout T) '
          f'{e_N:.2e} (layout N) | d fa: cos(meas,fa@K-K@fa)={c:+.6f} '
          f'ratio={np.linalg.norm(pred)/np.linalg.norm(dfa):.4f}')
mol.data['OQP::VEC_MO_A'] = C0
