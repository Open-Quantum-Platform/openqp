"""
Step 3 derisk: pin down the z-vector RHS for the orbital-response part of d_amp,
entirely in Python, correct layout (mo_a = VEC_MO_A.T).

Theory (standard Lagrangian/z-vector):
  d_amp(I,J) = X_I^T dX_J
             = [explicit]/gap  +  sum_pq L_pq (dC/dR)_pq /gap
  L_pq = dE_IJ/dU_pq   is the BARE orbital gradient (fixed Fock) -- the z-vector
         RHS.  The 2e/density Fock response is in the z-vector LHS (orbital
         Hessian), NOT in L.
  (dC/dR) is the orbital response.  M = C_ref^T S C_disp, U^x = dM/dR, but
  U^x = (MO-basis orbital response) + (C^T dS/dR C) -- the second (symmetric,
  = -S^x_MO/... ) piece is the overlap/Pulay term that belongs to d_ov, NOT here.
  So the orbital-response part of d_amp uses only the OCC-VIRT rotation blocks of
  L contracted with the response.

This test:
  1. explicit = mrsf_nac_amp (frozen-C explicit 2e).
  2. L = bare orbital gradient (correct layout).
  3. U^x = dM/dR ; split into antisymmetric (response) and symmetric (overlap).
  4. oracle = _compute_amp_damp (the validated production damp).
  5. Check which contraction of L with U^x, over which blocks, reproduces
     (oracle - explicit/gap).  Report block-resolved cos/ratio so the RHS is
     unambiguous before any Fortran.

Blocks (MRSF): doc=0..nocb, socc=nocb..noca, virt=noca..nbf.
Rotation (occ-virt-like) blocks: doc-socc, doc-virt, socc-virt.
"""
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint, NAC

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
EPS = 1.0e-4
DX = 1.0e-3
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p18.log')
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
pairs = [(I, J) for I in range(nstate) for J in range(nstate) if I != J]


def set_mo(moa, mob):
    mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(moa.T)
    mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(mob.T)


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


# ---- oracle (production damp) ----
oracle = NAC(mol)._compute_amp_damp(dx=1.0e-3)
set_mo(mo0, mo0b)

# ---- explicit 2e (mrsf_nac_amp) ----
oqp.mrsf_nac_amp(mol)
ana2e = np.array(mol.data['OQP::nac_amp'], copy=True).reshape(-1).reshape((nstate, nstate, natom, 3))

# ---- bare orbital gradient L (correct layout) ----
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
            e = E_all()
            for pr in pairs:
                L[pr][p, q] += sgn * e[pr] / (2.0 * EPS)
set_mo(mo0, mo0b)

# ---- U^x = dM/dR ----
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
    Ux[k] = (mo_overlap(xyz0 + DX * np.eye(ncoord)[k]) -
             mo_overlap(xyz0 - DX * np.eye(ncoord)[k])) / (2.0 * DX)
cfg['guess'].update(gb)
mol.update_system(xyz0)
oqp.library.ints_1e(mol)
oqp.library.guess(mol)
SinglePoint(mol).energy()
mol.data['OQP::td_bvec_mo'] = X0_raw

# ---- block masks ----
blocks = {'doc-socc': (slice(0, nocb), slice(nocb, noca)),
          'doc-virt': (slice(0, nocb), slice(noca, nbf)),
          'socc-virt': (slice(nocb, noca), slice(noca, nbf))}


def masked(M, sel, antisym):
    """Contract only the given rotation blocks (both orientations)."""
    out = np.zeros_like(M[0]) if M.ndim == 3 else np.zeros_like(M)
    return out


print('# oracle vs explicit2e/gap + orbital-response (block-resolved). '
      'Correct layout.')
for (I, J) in pairs:
    if I > J:
        continue
    gap = Om[J] - Om[I]
    orc = oracle[(I + 1, J + 1)].reshape(-1)
    exp2e = (ana2e[I, J] / gap).reshape(-1)
    # antisymmetric (response) part of U^x
    Ua = 0.5 * (Ux - np.transpose(Ux, (0, 2, 1)))
    # orbital response over occ-virt rotation blocks only, using antisym U^x
    resp = np.zeros(ncoord)
    for name, (so, sv) in blocks.items():
        for k in range(ncoord):
            resp[k] += (np.sum(L[(I, J)][so, sv] * Ua[k][so, sv]) +
                        np.sum(L[(I, J)][sv, so] * Ua[k][sv, so]))
    resp /= gap
    model = exp2e + resp
    c = model @ orc / (np.linalg.norm(model) * np.linalg.norm(orc) + 1e-300)
    print(f'  pair {(I+1, J+1)}: |oracle|={np.linalg.norm(orc):.5f} '
          f'|exp2e|={np.linalg.norm(exp2e):.5f} |resp|={np.linalg.norm(resp):.5f} '
          f'|model|={np.linalg.norm(model):.5f}  cos={c:+.5f} '
          f'ratio={np.linalg.norm(model)/(np.linalg.norm(orc)+1e-300):.4f}')

# save for offline block analysis
np.savez('/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p18_zvec.npz',
         Ux=Ux, **{f'L_{I+1}{J+1}': L[(I, J)] for (I, J) in pairs if I < J},
         **{f'oracle_{I+1}{J+1}': oracle[(I + 1, J + 1)] for (I, J) in pairs if I < J},
         **{f'ana2e_{I+1}{J+1}': ana2e[I, J] for (I, J) in pairs if I < J},
         Om=np.array(Om), nocb=nocb, noca=noca, nbf=nbf, natom=natom)
print('\nsaved -> data_snapshots/p18_zvec.npz (for offline block analysis)')
