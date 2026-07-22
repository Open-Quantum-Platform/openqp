"""
FULL closed-form model with the derived TRANSPORT term, tested against the oracle.

Derivation (on paper, 2026-07-21):
  oracle*gap = X_I^T A' X_J  +  Om_J (T'X_I)^T X_J + Om_I X_I^T (T'X_J)
             = [ana2e + esum + L:U]  +  transport
  Orthonormality (C^T S C = 1, dM = S^x_MO + U)  =>  U = -dM^T
  => L:U = -tr(L @ dM)
  (T'X_I)^T X_J = tr(dM_oo @ Xi Xj^T) + tr(dM_vv @ Xi^T Xj)   [Xi = SOMO-unfolded]

Tested as an ASSEMBLED TOTAL vs the oracle (not term-by-term vs a subtracted
target), so a missing piece shows up as a residual instead of being absorbed.
MRSF = ROHF + functional -> DFT inputs only.
"""
import sys
import numpy as np
import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint, NAC

INP = sys.argv[1]
TAG = INP.split('/')[-1].replace('.inp', '')
EPS = 1.0e-4
DX = 1.0e-3
OVTAG = 'OQP::overlap_mo_non_orthogonal'
SQ = 1.0 / np.sqrt(2.0)

r = Runner(input_file=INP, log=f'/bighome/alireza/openqp-nac/tools/nac_validation/p36_{TAG}.log')
r.run()
mol = r.mol
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
Craw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = Craw.shape[0]
mo0 = Craw.T.copy()
Crb = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
mo0b = Crb.T.copy()
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
nvirb = nbf - nocb
nij = noca * nvirb
bkey = 'OQP::td_bvec_mo'
X0_raw = np.array(mol.data[bkey], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
xyz0 = np.array(mol.get_system(), copy=True)
E = list(mol.energies)
Om = np.array([E[k + 1] - E[0] for k in range(nstate)])
nc = 3 * natom
pairs = [(i, j) for i in range(1, nstate + 1) for j in range(1, nstate + 1) if i < j]
print(f'# {TAG}: natom={natom} nbf={nbf} noca={noca} nocb={nocb} nstate={nstate}')


def unfold(xv):
    """Python replica of mrsfxvec (mult=1): SOMO fold of the amplitude vector."""
    t = xv.copy()
    ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1     # 1-based
    ijlr2 = (noca - nocb - 1) * noca + noca             # 1-based
    t[ijlr1 - 1] = xv[ijlr1 - 1] * SQ
    t[ijlr2 - 1] = -xv[ijlr1 - 1] * SQ
    return t


def as_mat(xv):
    """(noca, nvirb) with the Fortran ij = (j-nocb-1)*noca + i layout."""
    return xv.reshape((noca, nvirb), order='F')


def set_mo(moa, mob):
    mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(moa.T)
    mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(mob.T)


def set_bvec(col):
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = col
    mol.data[bkey] = rr.reshape(Xshape)


def Eij(I, J):
    set_bvec(X0[J])
    oqp.mrsf_matvec_apply(mol)
    return float(X0[I] @ np.array(mol.data['OQP::nac_mvax'], copy=True).ravel())


def cos(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-300))


# ---- (1) clean-orbital terms ----
mol.data[bkey] = X0_raw
oqp.mrsf_nac_amp(mol)
ana2e = np.array(mol.data['OQP::nac_amp'], copy=True).reshape(-1).reshape(
    (nstate, nstate, natom, 3))
esum = {}
for (i, j) in pairs:
    oqp.mrsf_nac_esum(mol, i, j)
    esum[(i, j)] = np.array(mol.data['OQP::nac_esum'], copy=True).reshape(-1)

# ---- (2) L (orbital-rotation FD of the matvec coupling) ----
L = {p: np.zeros((nbf, nbf)) for p in pairs}
for p in range(nbf):
    for q in range(nbf):
        if p == q:
            continue
        for sgn in (+1, -1):
            moa, mob = mo0.copy(), mo0b.copy()
            moa[:, q] += sgn * EPS * mo0[:, p]
            mob[:, q] += sgn * EPS * mo0b[:, p]
            set_mo(moa, mob)
            for (i, j) in pairs:
                L[(i, j)][p, q] += sgn * Eij(i - 1, j - 1) / (2.0 * EPS)
set_mo(mo0, mo0b)
mol.data[bkey] = X0_raw
print('# L done')

# ---- (3) dM = MO-overlap derivative (re-SCF) ----
cfg = mol.config
json0 = mol.log.replace('.log', '.json')
mol.save_data()
gb = dict(cfg['guess'])
cfg['guess']['type'] = 'json'
cfg['guess']['file'] = json0
cfg['guess']['continue_geom'] = False


def mo_overlap(coord):
    mol.update_system(coord)
    oqp.library.ints_1e(mol)
    oqp.library.guess(mol)
    SinglePoint(mol).energy()
    mol.data['OQP::xyz_old'] = xyz0.reshape((3, natom))
    mol.data['OQP::VEC_MO_A_old'] = Craw
    mol.data['OQP::E_MO_A_old'] = e0
    mol.data['OQP::VEC_MO_B_old'] = Crb
    mol.data['OQP::E_MO_B_old'] = e0b
    oqp.get_structures_ao_overlap(mol)
    return np.array(mol.data[OVTAG], copy=True).reshape(-1).reshape((nbf, nbf)).T


dM = np.zeros((nc, nbf, nbf))
for k in range(nc):
    dM[k] = (mo_overlap(xyz0 + DX * np.eye(nc)[k]) -
             mo_overlap(xyz0 - DX * np.eye(nc)[k])) / (2.0 * DX)
cfg['guess'].update(gb)
mol.update_system(xyz0)
oqp.library.ints_1e(mol)
oqp.library.guess(mol)
SinglePoint(mol).energy()
mol.data[bkey] = X0_raw
print('# dM done')

# ---- (4) oracle (last) ----
oracle = NAC(mol)._compute_amp_damp(dx=1.0e-3)
mol.data[bkey] = X0_raw

# ---- (5) assemble and compare ----
print()
SAVE = {}
print('# ASSEMBLED MODEL vs ORACLE   (all terms, no subtracted target)')
print(f'{"pair":>6} {"|orc*g|":>9} {"|ana2e|":>9} {"|esum|":>9} {"|L:U|":>9} '
      f'{"|transp|":>9} {"|model|":>9} {"cos":>8} {"ratio":>7} {"|resid|":>9}')
for (i, j) in pairs:
    gap = Om[j - 1] - Om[i - 1]
    orc = oracle[(i, j)].reshape(-1) * gap
    a2 = ana2e[i - 1, j - 1].reshape(-1)
    es = esum[(i, j)]
    Xi = as_mat(unfold(X0[i - 1]))
    Xj = as_mat(unfold(X0[j - 1]))
    Aoo = Xi @ Xj.T                     # (noca, noca)
    Avv = Xi.T @ Xj                     # (nvirb, nvirb)
    lu = np.zeros(nc)
    tr = np.zeros(nc)
    for k in range(nc):
        lu[k] = -np.trace(L[(i, j)] @ dM[k])
        doo = dM[k][0:noca, 0:noca]
        dvv = dM[k][nocb:nbf, nocb:nbf]
        gIJ = np.trace(doo @ Aoo) + np.trace(dvv @ Avv)
        gJI = np.trace(doo @ Aoo.T) + np.trace(dvv @ Avv.T)
        tr[k] = Om[j - 1] * gIJ + Om[i - 1] * gJI
    model = a2 + es + lu + tr
    resid = orc - model
    SAVE[f'orc_{i}{j}'] = orc; SAVE[f'ana2e_{i}{j}'] = a2
    SAVE[f'esum_{i}{j}'] = es; SAVE[f'lu_{i}{j}'] = lu; SAVE[f'tr_{i}{j}'] = tr
    print(f'{str((i,j)):>6} {np.linalg.norm(orc):9.5f} {np.linalg.norm(a2):9.5f} '
          f'{np.linalg.norm(es):9.5f} {np.linalg.norm(lu):9.5f} '
          f'{np.linalg.norm(tr):9.5f} {np.linalg.norm(model):9.5f} '
          f'{cos(model, orc):+8.4f} {np.linalg.norm(model)/(np.linalg.norm(orc)+1e-30):7.3f} '
          f'{np.linalg.norm(resid):9.5f}')
print('# cos~+1, ratio~1, |resid|<<|orc*g| => closed form COMPLETE.')
np.savez(f'/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p36_{TAG}.npz',
         Om=Om, natom=natom, nbf=nbf, noca=noca, nocb=nocb, **SAVE)
print(f'saved -> data_snapshots/p36_{TAG}.npz')
