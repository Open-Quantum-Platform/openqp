"""
The decisive test, CORRECT LAYOUT: the moving-frame chain rule for E_IJ = X_I^T A X_J

    dE/dR = dE/dR|_C(fixed orbitals) + sum_pq R_pq (dU_pq/dR)

with everything in one consistent frame:
  * total dE/dR : re-SCF at R+-dk, realign the displaced orbitals to the reference
      gauge continuously via mo_disp = mo0 @ M (M = C_ref^T S C_disp, correct
      layout), compute E_IJ, central-difference.
  * explicit dE/dR|_C : move nuclei + rebuild integrals but KEEP mo = mo0.
  * orbital : sum_pq R_pq U^x_pq, R_pq = dE/dU_pq (mo0 -> mo0(I+U)),
      U^x = dM/dR.

Layout (p15_combos): mo_a = VEC_MO_A.T ; the matvec computes mo_a^T FAO mo_a
EXACTLY (no re-canonicalization). If this chain rule CLOSES to FD floor, the
decomposition is validated and U^x can be replaced by the analytic z-vector.
"""
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
EPS = 1.0e-4
DX = 1.0e-3
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p17.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
Craw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = Craw.shape[0]
mo0 = Craw.T.copy()
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
pairs = [(I, J) for I in range(nstate) for J in range(nstate) if I < J]


def set_mo(moa, mob):
    mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(moa.T)
    mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(mob.T)


def set_bvec(col):
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)


def E_all():
    out, cache = {}, {}
    for (I, J) in [(a, b) for a in range(nstate) for b in range(nstate)]:
        if J not in cache:
            set_bvec(X0[J])
            oqp.mrsf_matvec_apply(mol)
            cache[J] = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
        out[(I, J)] = float(X0[I] @ cache[J])
    return out


# ---------- R_pq (orbital gradient, fixed geometry) ----------
R = {p: np.zeros((nbf, nbf)) for p in pairs}
for p in range(nbf):
    for q in range(nbf):
        if p == q:
            continue
        for sgn in (+1, -1):
            moa, mob = mo0.copy(), mo0b.copy()
            moa[:, q] += sgn * EPS * mo0[:, p]
            mob[:, q] += sgn * EPS * mo0b[:, p]
            set_mo(moa, mob)
            e = E_all()
            for (I, J) in pairs:
                R[(I, J)][p, q] += sgn * e[(I, J)] / (2.0 * EPS)
set_mo(mo0, mo0b)

# ---------- displaced quantities ----------
cfg = mol.config
json0 = mol.log.replace('.log', '.json')
mol.save_data()
gb = dict(cfg['guess'])
cfg['guess']['type'] = 'json'
cfg['guess']['file'] = json0
cfg['guess']['continue_geom'] = False


def at_geom(coord, mode):
    mol.update_system(coord)
    oqp.library.ints_1e(mol)
    if mode == 'explicit':
        set_mo(mo0, mo0b)                 # keep reference orbitals
        return E_all(), None
    oqp.library.guess(mol)
    SinglePoint(mol).energy()
    Cd = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    Cdb = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    mol.data['OQP::xyz_old'] = xyz0.reshape((3, natom))
    mol.data['OQP::VEC_MO_A_old'] = Craw
    mol.data['OQP::E_MO_A_old'] = e0
    mol.data['OQP::VEC_MO_B_old'] = Crb
    mol.data['OQP::E_MO_B_old'] = e0b
    oqp.get_structures_ao_overlap(mol)
    M = np.array(mol.data[OVTAG], copy=True).reshape(-1).reshape((nbf, nbf)).T
    # beta MO overlap: rebuild from Cd_b? get_structures uses alpha; approximate
    # beta by its own overlap via the same tag is not available -> use M for both
    # (ROHF: mo_a == mo_b for the spatial part, so M applies to both).
    moa = mo0 @ M
    mob = mo0b @ M
    set_mo(moa, mob)
    return E_all(), M


tot = {p: np.zeros(ncoord) for p in pairs}
exp = {p: np.zeros(ncoord) for p in pairs}
Ux = np.zeros((ncoord, nbf, nbf))
for k in range(ncoord):
    d = DX * np.eye(ncoord)[k]
    ep, Mp = at_geom(xyz0 + d, 'total')
    em, Mm = at_geom(xyz0 - d, 'total')
    xp, _ = at_geom(xyz0 + d, 'explicit')
    xm, _ = at_geom(xyz0 - d, 'explicit')
    Ux[k] = (Mp - Mm) / (2.0 * DX)
    for (I, J) in pairs:
        tot[(I, J)][k] = (ep[(I, J)] - em[(I, J)]) / (2.0 * DX)
        exp[(I, J)][k] = (xp[(I, J)] - xm[(I, J)]) / (2.0 * DX)

cfg['guess'].update(gb)
mol.update_system(xyz0)
oqp.library.ints_1e(mol)
oqp.library.guess(mol)
SinglePoint(mol).energy()
mol.data['OQP::td_bvec_mo'] = X0_raw

print('# MOVING-FRAME CHAIN RULE (correct layout): total == explicit + sum R.U^x')
print(f'# {"pair":>7} {"|total|":>11} {"|explicit|":>11} {"|R.Ux|":>11} '
      f'{"|resid|":>11} {"resid/tot":>10} {"cos":>10}')
out = {}
for (I, J) in pairs:
    rux = np.array([np.sum(R[(I, J)] * Ux[k]) for k in range(ncoord)])
    rhs = exp[(I, J)] + rux
    res = tot[(I, J)] - rhs
    c = tot[(I, J)] @ rhs / (np.linalg.norm(tot[(I, J)]) * np.linalg.norm(rhs) + 1e-300)
    print(f'  {str((I+1,J+1)):>7} {np.linalg.norm(tot[(I,J)]):11.6f} '
          f'{np.linalg.norm(exp[(I,J)]):11.6f} {np.linalg.norm(rux):11.6f} '
          f'{np.linalg.norm(res):11.6f} '
          f'{np.linalg.norm(res)/(np.linalg.norm(tot[(I,J)])+1e-300):10.4f} {c:+10.6f}')
    out[f'tot_{I+1}{J+1}'] = tot[(I, J)]
    out[f'exp_{I+1}{J+1}'] = exp[(I, J)]
    out[f'rux_{I+1}{J+1}'] = rux
np.savez('/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p17_chainrule.npz', **out)
