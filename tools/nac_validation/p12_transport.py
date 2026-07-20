"""
Phase 12 step 3b: verify the TRANSPORT relation, independently of any orbital
gradient.

The production oracle differentiates the TRANSPORTED operator A_ref = T^T A T,
where T is the amplitude-space transport (per-block Loewdin of the ref x displaced
MO overlap + SOMO det-grid fold).  Using A X_J = Om_J X_J, A^T = A and the
antisymmetry of tau = dT/dR:

    X_I^T (dA_ref/dR) X_J  =  X_I^T (dA/dR) X_J  +  (Om_I - Om_J) * X_I.tau.X_J   (*)

If (*) holds numerically, then the moving-frame derivative and the oracle differ by
EXACTLY one closed-form-able term, and

    d_amp = X_I^T(dA/dR)X_J/(Om_J-Om_I)  -  X_I.tau.X_J

That second piece is analytically tractable: dQ = antisym(dM) for a symmetric
(Loewdin) orthogonalisation, and X_I.tau.X_J contracts the interstate densities
(exactly the tij/tab that mrsf_interstate_tden already builds) with the
antisymmetric MO-overlap derivative in the doc/socc/virt blocks.

This isolates the transport term ALONE - no R, no U^x, no density channel - so it
is cheap (no nbf^2 perturbation loop) and it either confirms or kills the
structural claim before any Fortran is written.
"""
import math
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint, NAC

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
OUT = '/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p12_transport.npz'
DX = 1.0e-3
SQ = 1.0 / math.sqrt(2.0)
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=INP,
           log='/bighome/alireza/openqp-nac/tools/nac_validation/p12_transport.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20

nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
C0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
nbf = C0.shape[0]
nvirb = nbf - nocb
nij = noca * nvirb
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
xyz0 = np.array(mol.get_system(), copy=True)
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]
ncoord = 3 * natom
ijlr1 = (noca - 1 - nocb - 1) * noca + (noca - 1) - 1
ijlr2 = (noca - nocb - 1) * noca + (noca) - 1
pairs = [(I, J) for I in range(nstate) for J in range(nstate) if I < J]


def set_bvec(col):
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)


def Acol(c):
    set_bvec(c)
    oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()


def unfold_det(col):
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


def refold_det(cp):
    g = cp.copy()
    i1, a1 = ijlr1 % noca, ijlr1 // noca
    g[i1, a1] = math.sqrt(2.0) * cp[i1, a1]
    i2, a2 = ijlr2 % noca, ijlr2 // noca
    g[i2, a2] = 0.0
    return g.T.reshape(-1)


def transport_T(Q):
    Qo = Q[0:noca, 0:noca]
    Qv = Q[nocb:nbf, nocb:nbf]
    T = np.zeros((nij, nij))
    e = np.zeros(nij)
    for j in range(nij):
        e[:] = 0.0
        e[j] = 1.0
        T[:, j] = refold_det(Qo.T @ unfold_det(e) @ Qv)
    return T


cfg = mol.config
json0 = mol.log.replace('.log', '.json')
mol.save_data()
guess_bak = dict(cfg['guess'])
cfg['guess']['type'] = 'json'
cfg['guess']['file'] = json0
cfg['guess']['continue_geom'] = False


def at_geom(coord):
    """Re-SCF at `coord`; return (A in the displaced/moving frame, transport T)."""
    mol.update_system(coord)
    oqp.library.ints_1e(mol)
    oqp.library.guess(mol)
    SinglePoint(mol).energy()
    mol.data['OQP::xyz_old'] = xyz0.reshape((3, natom))
    mol.data['OQP::VEC_MO_A_old'] = C0
    mol.data['OQP::E_MO_A_old'] = e0
    mol.data['OQP::VEC_MO_B_old'] = C0b
    mol.data['OQP::E_MO_B_old'] = e0b
    oqp.get_structures_ao_overlap(mol)
    M = np.array(mol.data[OVTAG], copy=True).reshape(-1).reshape((nbf, nbf)).T
    Q = np.zeros((nbf, nbf))
    for lo, hi in ((0, nocb), (nocb, noca), (noca, nbf)):
        sub = M[lo:hi, lo:hi]
        w, U = np.linalg.eigh(sub.T @ sub)
        Q[lo:hi, lo:hi] = sub @ (U @ np.diag(1.0 / np.sqrt(w)) @ U.T)
    A = np.zeros((nij, nij))
    e = np.zeros(nij)
    for j in range(nij):
        e[:] = 0.0
        e[j] = 1.0
        A[:, j] = Acol(e)
    return A, transport_T(Q)


dA_mov = np.zeros((ncoord, nij, nij))     # moving frame
dA_ref = np.zeros((ncoord, nij, nij))     # transported (what the oracle uses)
tau = np.zeros((ncoord, nij, nij))
for k in range(ncoord):
    d = DX * np.eye(ncoord)[k]
    Ap, Tp = at_geom(xyz0 + d)
    Am, Tm = at_geom(xyz0 - d)
    dA_mov[k] = (Ap - Am) / (2 * DX)
    dA_ref[k] = (Tp.T @ Ap @ Tp - Tm.T @ Am @ Tm) / (2 * DX)
    tau[k] = (Tp - Tm) / (2 * DX)

cfg['guess'].update(guess_bak)
mol.update_system(xyz0)
oqp.library.ints_1e(mol)
oqp.library.guess(mol)
SinglePoint(mol).energy()
mol.data['OQP::td_bvec_mo'] = X0_raw

print('# TRANSPORT RELATION:  X_I^T dA_ref X_J  ==  X_I^T dA_mov X_J '
      '+ (Om_I-Om_J) X_I.tau.X_J')
print(f'# {"pair":>7} {"|lhs|":>11} {"|mov|":>11} {"|gauge|":>11} '
      f'{"|resid|":>11} {"resid/|lhs|":>12} {"cos":>11}')
out = {}
for (I, J) in pairs:
    lhs = np.array([X0[I] @ dA_ref[k] @ X0[J] for k in range(ncoord)])
    mov = np.array([X0[I] @ dA_mov[k] @ X0[J] for k in range(ncoord)])
    gge = (Om[I] - Om[J]) * np.array([X0[I] @ tau[k] @ X0[J] for k in range(ncoord)])
    res = lhs - (mov + gge)
    c = lhs @ (mov + gge) / (np.linalg.norm(lhs) * np.linalg.norm(mov + gge) + 1e-300)
    print(f'  {str((I+1, J+1)):>7} {np.linalg.norm(lhs):11.6f} {np.linalg.norm(mov):11.6f} '
          f'{np.linalg.norm(gge):11.6f} {np.linalg.norm(res):11.6f} '
          f'{np.linalg.norm(res)/(np.linalg.norm(lhs)+1e-300):12.5f} {c:+11.6f}')
    out[f'lhs_{I+1}{J+1}'] = lhs
    out[f'mov_{I+1}{J+1}'] = mov
    out[f'gauge_{I+1}{J+1}'] = gge
np.savez(OUT, **out)
print(f'\nsaved -> {OUT}')
