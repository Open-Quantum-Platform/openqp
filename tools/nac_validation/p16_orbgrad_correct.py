"""
With the CORRECT layout (mo_a = VEC_MO_A.T, fa read as raw.T), re-measure the
interstate orbital gradient of the FULL matvec and test it against the U^x oracle.

Layout established (p15_combos.py): fa_export.T = mo_a^T FAO mo_a EXACTLY
(1.9e-15), mo_a = C.T. So the matvec does NOT re-canonicalize; all Phase 12/13
failures were this transpose. Rotations now transform exactly.

E_IJ = X_I^T A X_J (A = full MRSF matvec). R_pq = dE_IJ/dU_pq measured by rotating
the MO index of mo_a (= columns of mo_a = rows of C). Then contract with the
oracle orbital response U^x and compare to the saved target
`missing = oracle - ana2e/gap`  (data_snapshots/p12_ux_target.npz).

This is the same test as p12_orbgrad, but with the layout fixed -- so if the only
bug was the transpose, this should now CLOSE.
"""
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
TGT = '/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p12_ux_target.npz'
EPS = 1.0e-4
DX = 1.0e-3
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p16.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
Craw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)     # numpy layout
nbf = Craw.shape[0]
mo0 = Craw.T.copy()                                        # mo_a[:,k] = MO k
Crb = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
mo0b = Crb.T.copy()
nvirb = nbf - nocb
nij = noca * nvirb
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
xyz0 = np.array(mol.get_system(), copy=True)
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]
ncoord = 3 * natom
pairs = [(I, J) for I in range(nstate) for J in range(nstate) if I != J]


def set_mo(mo_a, mo_b):
    mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(mo_a.T)
    mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(mo_b.T)


def set_bvec(col):
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)


def E_all():
    out, cache = {}, {}
    for (I, J) in pairs:
        if J not in cache:
            set_bvec(X0[J])
            oqp.mrsf_matvec_apply(mol)
            cache[J] = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
        out[(I, J)] = float(X0[I] @ cache[J])
    return out


# ---- R_pq by rotating the MO index of mo_a (columns of mo_a) ----
R = {p: np.zeros((nbf, nbf)) for p in pairs}
for p in range(nbf):
    for q in range(nbf):
        if p == q:
            continue
        for sgn in (+1, -1):
            moa = mo0.copy()
            mob = mo0b.copy()
            moa[:, q] += sgn * EPS * mo0[:, p]     # rotate MO q by MO p
            mob[:, q] += sgn * EPS * mo0b[:, p]
            set_mo(moa, mob)
            e = E_all()
            for pr in pairs:
                R[pr][p, q] += sgn * e[pr] / (2.0 * EPS)
set_mo(mo0, mo0b)

# ---- U^x from the ref-vs-displaced MO overlap (central diff) ----
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


Ux = np.zeros((ncoord, nbf, nbf))
for k in range(ncoord):
    Mp = mo_overlap(xyz0 + DX * np.eye(ncoord)[k])
    Mm = mo_overlap(xyz0 - DX * np.eye(ncoord)[k])
    Ux[k] = (Mp - Mm) / (2.0 * DX)

cfg['guess'].update(gb)
mol.update_system(xyz0)
oqp.library.ints_1e(mol)
oqp.library.guess(mol)
SinglePoint(mol).energy()
mol.data['OQP::td_bvec_mo'] = X0_raw

tgt = np.load(TGT)
print(f'# CORRECT-LAYOUT orbital gradient: sum R.U^x  vs  missing = oracle - ana2e/gap')
print(f'# {"pair":>7} {"|R.Ux|":>11} {"|missing|":>11} {"cos":>11} {"ratio":>9}')
for (I, J) in pairs:
    key = f'missing_{I+1}{J+1}'
    if key not in tgt:
        continue
    gap = Om[J] - Om[I]
    contrib = np.array([np.sum(R[(I, J)] * Ux[k]) for k in range(ncoord)]) / gap
    miss = tgt[key]
    c = contrib @ miss / (np.linalg.norm(contrib) * np.linalg.norm(miss) + 1e-300)
    print(f'  {str((I+1,J+1)):>7} {np.linalg.norm(contrib):11.6f} {np.linalg.norm(miss):11.6f} '
          f'{c:+11.6f} {np.linalg.norm(contrib)/(np.linalg.norm(miss)+1e-300):9.4f}')
