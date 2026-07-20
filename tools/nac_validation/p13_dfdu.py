"""
Phase 13 step 2: check the Fock-derivative formula BEFORE contracting with Gam.

The matvec builds the MO Fock as  fa = C^T F_AO C  and exports it as OQP::nac_fa.
Perturb the MOs  C[:,q] += eps*C[:,p]  (Fock frozen) and the claim is

    d fa_rs / dU_pq  =  delta_sq * fa_rp  +  delta_rq * fa_ps

i.e. adding column p into column q inserts (fa row/col p) into (fa row/col q).
This is the sole analytic ingredient behind R_esum; Gate 2 failed, so test THIS
in isolation, with no Gam anywhere, so a mismatch localises the error to the
Fock derivative rather than to Gam or the trace.

Two checks:
  A. exact-index spot check: does d fa/dU_pq match delta_sq fa_rp + delta_rq fa_ps?
  B. the contraction actually used: dE_esum/dU_pq = sum_rs Gam_sr d fa_rs/dU_pq.
     Under the formula this is (Gam^T fa)_qp + (Gam fa)_qp = [(Gam+Gam^T) fa]_qp
     ONLY if Gam is symmetric.  Report Gam's asymmetry too -- if Gam is NOT
     symmetric the correct contraction is (Gam^T fa + fa Gam^T? ) and the closed
     form must keep Gam and Gam^T distinct.
"""
import math
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
EPS = 1.0e-4

r = Runner(input_file=INP,
           log='/bighome/alireza/openqp-nac/tools/nac_validation/p13_dfdu.log')
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
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
SQ = 1.0 / math.sqrt(2.0)
ijlr1 = (noca - 1 - nocb - 1) * noca + (noca - 1) - 1
ijlr2 = (noca - nocb - 1) * noca + (noca) - 1


def set_bvec(col):
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)


def get_fa_fb():
    """Read the MO Fock fa,fb the matvec actually built at the current MOs."""
    set_bvec(X0[0])
    oqp.mrsf_matvec_apply(mol)
    fa = np.array(mol.data['OQP::nac_fa'], copy=True).reshape(-1).reshape((nbf, nbf)).T
    fb = np.array(mol.data['OQP::nac_fb'], copy=True).reshape(-1).reshape((nbf, nbf)).T
    return fa, fb


def put_C(Ca, Cb):
    mol.data['OQP::VEC_MO_A'] = Ca
    mol.data['OQP::VEC_MO_B'] = Cb


fa0, fb0 = get_fa_fb()

# ------------------------------------------------------- reference: fa = C^T F_AO C ?
# OQP::FOCK_A is packed lower-triangular (nbf*(nbf+1)/2); unpack to full.
faopack = np.array(mol.data['OQP::FOCK_A'], copy=True).ravel()
FAO = np.zeros((nbf, nbf))
idx = 0
for i in range(nbf):
    for j in range(i + 1):
        FAO[i, j] = FAO[j, i] = faopack[idx]
        idx += 1
fa_from_C = C0.T @ FAO @ C0
print(f'# |fa0 - C^T FAO C| = {np.abs(fa0 - fa_from_C).max():.3e}  '
      f'(if ~0, fa0 IS C^T F_AO C and the derivative formula applies)')
print(f'# fa asymmetry |fa0 - fa0^T| = {np.abs(fa0 - fa0.T).max():.3e}')

# ------------------------------------------------------- CHECK A: exact index
print('\n# CHECK A: d fa/dU_pq  vs  delta_sq fa_rp + delta_rq fa_ps')
print(f'# {"(p,q)":>8} {"|dfa_meas|":>12} {"|dfa_pred|":>12} {"cos":>11} {"ratio":>9}')
worst = 0.0
for (p, q) in [(2, 8), (0, 5), (7, 3), (5, 5), (1, 12)]:
    Ca = C0.copy()
    Ca[:, q] += EPS * C0[:, p]
    put_C(Ca, C0b)
    fp, _ = get_fa_fb()
    Ca = C0.copy()
    Ca[:, q] -= EPS * C0[:, p]
    put_C(Ca, C0b)
    fm, _ = get_fa_fb()
    put_C(C0, C0b)
    dfa = (fp - fm) / (2 * EPS)
    pred = np.zeros((nbf, nbf))
    pred[:, q] += fa0[:, p]       # delta_sq fa_rp
    pred[q, :] += fa0[p, :]       # delta_rq fa_ps
    c = np.sum(dfa * pred) / (np.linalg.norm(dfa) * np.linalg.norm(pred) + 1e-300)
    worst = max(worst, 1 - c)
    print(f'  {str((p, q)):>8} {np.linalg.norm(dfa):12.6f} {np.linalg.norm(pred):12.6f} '
          f'{c:+11.7f} {np.linalg.norm(pred)/(np.linalg.norm(dfa)+1e-300):9.5f}')
print(f'# worst (1-cos) over the spot checks = {worst:.3e}')

# ------------------------------------------------------- CHECK B: Gam symmetry
Xt = []
for k in range(nstate):
    x = np.zeros((noca, nvirb))
    for i in range(noca):
        for a in range(nvirb):
            ij = a * noca + i
            if ij == ijlr1:
                x[i, a] = X0[k][ijlr1] * SQ
            elif ij == ijlr2:
                x[i, a] = -X0[k][ijlr1] * SQ
            else:
                x[i, a] = X0[k][ij]
    Xt.append(x)

print('\n# CHECK B: is Gam symmetric? (decides (Gam+Gam^T) vs keeping Gam,Gam^T distinct)')
for I in range(nstate):
    GA = np.zeros((nbf, nbf))
    GA[0:noca, 0:noca] = -0.5 * (Xt[I] @ Xt[I].T + Xt[I] @ Xt[I].T)
    print(f'  state {I+1}: Gam_A occ-occ asym |G-G^T| = {np.abs(GA - GA.T).max():.3e}')
