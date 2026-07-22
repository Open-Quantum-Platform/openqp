"""
CORRECTED transport term, from reading _compute_amp_damp line by line.

Three errors fixed vs p36:
 (1) ★ The oracle's transport Q is a PER-BLOCK LOEWDIN orthogonalization of the MO
     overlap M (single_point.py:1961-1966), on blocks [0,nocb) [nocb,noca) [noca,nbf):
         Q_blk = sub (sub^T sub)^(-1/2)
     d/dR of R = A(A^T A)^(-1/2) at A=I is 1/2(dA - dA^T), so
         dQ = BLOCK-DIAGONAL ANTISYMMETRIC part of dM
     (the symmetric S^x_MO half is projected OUT, off-block elements discarded).
     p36 used the raw full dM -> wrong object.
 (2) transport_vec is refold( Qo^T @ unfold(x) @ Qv ): Qo TRANSPOSED, Qv not.
 (3) refold is NOT unfold^T. <Refold(Y),X_J> = <Y, Refold^T(X_J)>, so the KET side
     needs Refold^T (entries as-is, ijlr1 * sqrt2, ijlr2 zeroed), not unfold.

Only dM is recomputed here (36 re-SCF). L is NOT recomputed: L:U = -tr(L@dM) was
already correct in form, so the cached `lu` from p36_eth_ana.npz is reused (ethylene
is deterministic -- p36/p37 gave identical numbers).
"""
import sys
import numpy as np
import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint

INP = sys.argv[1]
CACHE = '/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p36_eth_ana.npz'
DX = 1.0e-3
OVTAG = 'OQP::overlap_mo_non_orthogonal'
SQ = 1.0 / np.sqrt(2.0)

r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p38.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
Craw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = Craw.shape[0]
Crb = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
nvirb = nbf - nocb
nij = noca * nvirb
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
xyz0 = np.array(mol.get_system(), copy=True)
E = list(mol.energies)
Om = np.array([E[k + 1] - E[0] for k in range(nstate)])
nc = 3 * natom
pairs = [(i, j) for i in range(1, nstate + 1) for j in range(1, nstate + 1) if i < j]
# 0-based, exactly as in _compute_amp_damp lines 1871-1872
ijlr1 = (noca - 1 - nocb - 1) * noca + (noca - 1) - 1
ijlr2 = (noca - nocb - 1) * noca + (noca) - 1

cache = np.load(CACHE)
print(f'# {INP.split("/")[-1]}  nbf={nbf} nstate={nstate}')
print(f'# determinism check: max|Om_now - Om_cached| = '
      f'{np.max(np.abs(Om - cache["Om"])):.3e}')


def unfold(col):
    """Exact replica of _compute_amp_damp.unfold_det."""
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


def refoldT(col):
    """Refold^T: entries as-is, ijlr1 scaled by sqrt2, ijlr2 zeroed."""
    m = col.reshape((noca, nvirb), order='F').copy()
    i1, a1 = ijlr1 % noca, ijlr1 // noca
    m[i1, a1] = np.sqrt(2.0) * col[ijlr1]
    i2, a2 = ijlr2 % noca, ijlr2 // noca
    m[i2, a2] = 0.0
    return m


def cos(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-300))


# ---- dM (MO-overlap derivative) ----
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
print('# dM done')

# ---- dQ = block-diagonal antisymmetric part of dM (Loewdin derivative) ----
BLOCKS = ((0, nocb), (nocb, noca), (noca, nbf))
dQ = np.zeros_like(dM)
for k in range(nc):
    for lo, hi in BLOCKS:
        sub = dM[k][lo:hi, lo:hi]
        dQ[k][lo:hi, lo:hi] = 0.5 * (sub - sub.T)

# ---- corrected transport + assembly ----
print()
print('# CORRECTED transport (Loewdin dQ, right index convention, Refold^T ket)')
print(f'{"pair":>6} {"|orc|":>8} {"|tr_old|":>9} {"|tr_new|":>9} '
      f'{"cos(base)":>10} {"r(base)":>8} {"cos(new)":>9} {"r(new)":>8} {"cos(old)":>9}')
OUT = {}
for (i, j) in pairs:
    orc = cache[f'orc_{i}{j}']
    a2 = cache[f'ana2e_{i}{j}']
    es = cache[f'esum_{i}{j}']
    lu = cache[f'lu_{i}{j}']
    tr_old = cache[f'tr_{i}{j}']
    Ui, Uj = unfold(X0[i - 1]), unfold(X0[j - 1])
    Zi, Zj = refoldT(X0[i - 1]), refoldT(X0[j - 1])
    tr_new = np.zeros(nc)
    for k in range(nc):
        dQo = dQ[k][0:noca, 0:noca]
        dQv = dQ[k][nocb:nbf, nocb:nbf]
        gIJ = np.sum((dQo.T @ Ui) * Zj) + np.sum((Ui @ dQv) * Zj)
        gJI = np.sum((dQo.T @ Uj) * Zi) + np.sum((Uj @ dQv) * Zi)
        tr_new[k] = Om[j - 1] * gIJ + Om[i - 1] * gJI
    base = a2 + es + lu                      # no transport at all
    new = base + tr_new
    old = base + tr_old
    OUT[f'trnew_{i}{j}'] = tr_new
    print(f'{str((i,j)):>6} {np.linalg.norm(orc):8.5f} {np.linalg.norm(tr_old):9.5f} '
          f'{np.linalg.norm(tr_new):9.5f} {cos(base,orc):+10.4f} '
          f'{np.linalg.norm(base)/np.linalg.norm(orc):8.3f} {cos(new,orc):+9.4f} '
          f'{np.linalg.norm(new)/np.linalg.norm(orc):8.3f} {cos(old,orc):+9.4f}')
print('# base = ana2e+esum+L:U (no transport); new = base+corrected; old = base+p36 transport')
np.savez('/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p38_eth.npz',
         dM=dM, dQ=dQ, X0=X0, Om=Om, noca=noca, nocb=nocb, nbf=nbf, natom=natom, **OUT)
print('saved -> data_snapshots/p38_eth.npz  (dM/dQ/X0 cached for offline work)')
