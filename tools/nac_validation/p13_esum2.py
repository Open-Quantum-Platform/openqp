"""
Phase 13 step 1: derive and VALIDATE the esum half of the interstate orbital
gradient in closed form.

STRUCTURE.  The MRSF matvec is  A = A_2e + A_esum, and A_esum is LINEAR in the
frozen AO Fock, which mrsf_matvec_apply reads from the OQP::FOCK_A/B tags.  So
scaling the Fock isolates it EXACTLY:

    E_esum(I,J) = [ E_IJ(F*(1+eta)) - E_IJ(F) ] / eta          (exact, linear)
    E_2e  (I,J) =   E_IJ(F) - E_esum(I,J)

That is a clean handle on half of A that involves NO displaced SCF, so none of
the orbital-phase contamination that invalidated Phase 12 can reach it.

CLOSED FORM.  E_esum = Tr(Gam_A . fa) + Tr(Gam_B . fb) with fa = C^T F_AO C and
Gam the interstate density, which is independent of both R and C.  Perturbing
C[:,q] += eps*C[:,p]:

    d fa_rs / dU_pq = delta_sq fa_rp + delta_rq fa_ps
 => dE_esum/dU_pq   = (Gam^T fa)_qp + (Gam fa)_qp = [ (Gam+Gam^T) fa ]_qp

GATES (both sharp, both cheap, all at the reference geometry):
  G1  Tr(Gam fa)+Tr(Gam fb) must equal the Fock-scaling measurement of E_esum.
      This also PINS the Gam convention (the sqrt2 SOMO fold and the +-1/2
      factors), which is otherwise guesswork.
  G2  [(Gam+Gam^T) fa]_qp must equal the Fock-scaling measurement of R_esum,
      i.e. of dE_esum/dU_pq.
If G1 and G2 pass, the esum half of the orbital-gradient RHS is DONE in closed
form and only the 2e channel half remains.
"""
import math
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
OUT = '/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p13_esum2.npz'
ETA = 1.0e-4     # Fock scaling (exact in principle; small only to stay in float range)
EPS = 1.0e-4     # MO mixing

r = Runner(input_file=INP,
           log='/bighome/alireza/openqp-nac/tools/nac_validation/p13_esum2.log')
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
pairs = [(I, J) for I in range(nstate) for J in range(nstate) if I < J]
diag  = [(I, I) for I in range(nstate)]


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


def E_and_fock(I, J):
    """(E_IJ, fa, fb) at the CURRENT C and the CURRENT scaled Fock."""
    set_bvec(X0[J])
    oqp.mrsf_matvec_apply(mol)
    ax = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
    fa = np.array(mol.data['OQP::nac_fa'], copy=True).reshape(-1).reshape((nbf, nbf)).T
    fb = np.array(mol.data['OQP::nac_fb'], copy=True).reshape(-1).reshape((nbf, nbf)).T
    return float(X0[I] @ ax), fa, fb


def E_esum_measured(I, J):
    """Isolate the esum part by Fock scaling (A_esum is linear in F_AO)."""
    put_F(1.0)
    e0, _, _ = E_and_fock(I, J)
    put_F(1.0 + ETA)
    e1, _, _ = E_and_fock(I, J)
    put_F(1.0)
    return (e1 - e0) / ETA, e0


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


# ---------------------------------------------------------------- GATE 1
put_F(1.0)
_, fa0, fb0 = E_and_fock(0, 1)
Xt = [unfold(X0[k]) for k in range(nstate)]

print('# GATE 1: pin the Gam convention against the Fock-scaling measurement')
print(f'# {"pair":>7} {"E_esum(meas)":>14} {"best cand":>12} {"Tr(Gam f)":>14} {"rel.err":>11}')
# candidate conventions for (alpha occ-occ, beta virt-virt) prefactors
CANDS = {'-1/2,+1/2': (-0.5, 0.5), '-1,+1': (-1.0, 1.0),
         '-1/4,+1/4': (-0.25, 0.25), '+1/2,-1/2': (0.5, -0.5)}
best = {}
for (I, J) in diag:
    meas, etot = E_esum_measured(I, J)
    res = {}
    for nm, (ca, cb) in CANDS.items():
        GA = np.zeros((nbf, nbf))
        GB = np.zeros((nbf, nbf))
        mo = 0.5 * (Xt[I] @ Xt[J].T + Xt[J] @ Xt[I].T)      # (noca x noca)
        mv = 0.5 * (Xt[I].T @ Xt[J] + Xt[J].T @ Xt[I])      # (nvirb x nvirb)
        GA[0:noca, 0:noca] = ca * mo
        GB[nocb:nbf, nocb:nbf] = cb * mv
        res[nm] = np.trace(GA @ fa0) + np.trace(GB @ fb0)
    nm = min(res, key=lambda k: abs(res[k] - meas))
    best[(I, J)] = (nm, res[nm], meas)
    print(f'  {str((I+1,J+1)):>7} {meas:14.8f} {nm:>12} {res[nm]:14.8f} '
          f'{abs(res[nm]-meas)/(abs(meas)+1e-30):11.2e}')

conv = max(set(b[0] for b in best.values()),
           key=lambda c: sum(1 for b in best.values() if b[0] == c))
ca, cb = CANDS[conv]
print(f'\n# chosen convention: {conv}')


def Gam(I, J):
    GA = np.zeros((nbf, nbf))
    GB = np.zeros((nbf, nbf))
    GA[0:noca, 0:noca] = ca * 0.5 * (Xt[I] @ Xt[J].T + Xt[J] @ Xt[I].T)
    GB[nocb:nbf, nocb:nbf] = cb * 0.5 * (Xt[I].T @ Xt[J] + Xt[J].T @ Xt[I])
    return GA, GB


# ---------------------------------------------------------------- GATE 2
print('\n# GATE 2 (on the DIAGONAL, where E_esum is nonzero): '
      'closed-form R_esum = [(Gam+Gam^T) f]_qp vs measurement')
print(f'# {"pair":>7} {"|R_meas|":>12} {"|R_closed|":>12} {"cos":>11} {"ratio":>9}')
out = {}
for (I, J) in diag:
    Rm = np.zeros((nbf, nbf))
    for p in range(nbf):
        for q in range(nbf):
            for sgn in (+1, -1):
                Ca, Cb = C0.copy(), C0b.copy()
                Ca[:, q] += sgn * EPS * C0[:, p]
                Cb[:, q] += sgn * EPS * C0b[:, p]
                put_C(Ca, Cb)
                e, _ = E_esum_measured(I, J)
                Rm[p, q] += sgn * e / (2.0 * EPS)
    put_C(C0, C0b)
    put_F(1.0)
    GA, GB = Gam(I, J)
    Rc = ((GA + GA.T) @ fa0).T + ((GB + GB.T) @ fb0).T
    c = np.sum(Rm * Rc) / (np.linalg.norm(Rm) * np.linalg.norm(Rc) + 1e-300)
    print(f'  {str((I+1,J+1)):>7} {np.linalg.norm(Rm):12.6f} {np.linalg.norm(Rc):12.6f} '
          f'{c:+11.7f} {np.linalg.norm(Rc)/(np.linalg.norm(Rm)+1e-300):9.5f}')
    out[f'Rmeas_{I+1}{J+1}'] = Rm
    out[f'Rclosed_{I+1}{J+1}'] = Rc
np.savez(OUT, **out)
print(f'\nsaved -> {OUT}')
