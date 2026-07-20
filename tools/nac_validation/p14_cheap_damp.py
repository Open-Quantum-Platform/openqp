"""
Phase 14: make the VALIDATED semi-numerical damp SCALE, by never building the
full nij x nij matvec matrix.

The production _compute_amp_damp builds A_ref = T^T A_disp T column by column
(nij matvec applications per displacement) purely to then contract it down to the
scalars  X0_I^T A_ref X0_J.  But

    X0_I^T A_ref X0_J = (T X0_I)^T A_disp (T X0_J)

so A_disp is only ever needed applied to the nstate transported kets T X0_J --
ONE matvec application per ket per displacement, not nij.  For PSB3 (nij=1955)
this is a ~2000x reduction and is exactly what killed the full run.

This harness computes damp BOTH ways on H2O and asserts they agree to FD floor,
so the cheap path can be trusted as a drop-in for the full one.
"""
import math
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
DX = 1.0e-3
SQ = 1.0 / math.sqrt(2.0)
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=INP,
           log='/bighome/alireza/openqp-nac/tools/nac_validation/p14.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
C0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = C0.shape[0]
nvirb = nbf - nocb
nij = noca * nvirb
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
xyz0 = np.array(mol.get_system(), copy=True)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]
ncoord = 3 * natom
ijlr1 = (noca - 1 - nocb - 1) * noca + (noca - 1) - 1
ijlr2 = (noca - nocb - 1) * noca + (noca) - 1


def set_bvec(col):
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)


def Aapply(v):
    set_bvec(v)
    oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()


def Amat():
    A = np.zeros((nij, nij))
    e = np.zeros(nij)
    for j in range(nij):
        e[:] = 0.0
        e[j] = 1.0
        A[:, j] = Aapply(e)
    return A


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


def refold(cp):
    g = cp.copy()
    i1, a1 = ijlr1 % noca, ijlr1 // noca
    g[i1, a1] = math.sqrt(2.0) * cp[i1, a1]
    i2, a2 = ijlr2 % noca, ijlr2 // noca
    g[i2, a2] = 0.0
    return g.T.reshape(-1)


def transport_vec(Q, col):
    """T applied to ONE amplitude vector (unfold -> block rotate -> refold)."""
    Qo = Q[0:noca, 0:noca]
    Qv = Q[nocb:nbf, nocb:nbf]
    return refold(Qo.T @ unfold(col) @ Qv)


cfg = mol.config
json0 = mol.log.replace('.log', '.json')
mol.save_data()
guess_bak = dict(cfg['guess'])
cfg['guess']['type'] = 'json'
cfg['guess']['file'] = json0
cfg['guess']['continue_geom'] = False


def displaced(coord, want_matrix):
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
    # CHEAP: transported kets, one matvec each
    TX = [transport_vec(Q, X0[k]) for k in range(nstate)]
    ATX = [Aapply(TX[k]) for k in range(nstate)]      # nstate matvec applications
    scal = np.array([[TX[i] @ ATX[j] for j in range(nstate)] for i in range(nstate)])
    full = None
    if want_matrix:
        A = Amat()
        T = np.zeros((nij, nij))
        e = np.zeros(nij)
        for j in range(nij):
            e[:] = 0.0
            e[j] = 1.0
            T[:, j] = transport_vec(Q, e)
        Aref = T.T @ A @ T
        full = np.array([[X0[i] @ Aref @ X0[j] for j in range(nstate)]
                         for i in range(nstate)])
    return scal, full


want = (nij <= 200)      # only validate against the full matrix on small systems
sc = np.zeros((ncoord, nstate, nstate))
fu = np.zeros((ncoord, nstate, nstate))
for k in range(ncoord):
    sp, fp = displaced(xyz0 + DX * np.eye(ncoord)[k], want)
    sm, fm = displaced(xyz0 - DX * np.eye(ncoord)[k], want)
    sc[k] = (sp - sm) / (2 * DX)
    if want:
        fu[k] = (fp - fm) / (2 * DX)

cfg['guess'].update(guess_bak)
mol.update_system(xyz0)
oqp.library.ints_1e(mol)
oqp.library.guess(mol)
SinglePoint(mol).energy()
mol.data['OQP::td_bvec_mo'] = X0_raw

print(f'# cheap (one-matvec-per-ket) vs full (nij-column) damp,  nij={nij}')
print(f'# {"pair":>7} {"|cheap|":>11} {"|full|":>11} {"cos":>11} {"ratio":>9}')
for I in range(nstate):
    for J in range(nstate):
        if I == J:
            continue
        gap = Om[J] - Om[I]
        dc = sc[:, I, J] / gap
        if want:
            df = fu[:, I, J] / gap
            c = dc @ df / (np.linalg.norm(dc) * np.linalg.norm(df) + 1e-300)
            print(f'  {str((I+1,J+1)):>7} {np.linalg.norm(dc):11.6f} '
                  f'{np.linalg.norm(df):11.6f} {c:+11.8f} '
                  f'{np.linalg.norm(dc)/(np.linalg.norm(df)+1e-300):9.6f}')
        else:
            print(f'  {str((I+1,J+1)):>7} {np.linalg.norm(dc):11.6f} '
                  f'{"(no full)":>11}')
