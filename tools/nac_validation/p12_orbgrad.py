"""
Phase 12 step 1-2: measure the interstate ORBITAL GRADIENT and test whether it,
contracted with the numerical orbital response, reproduces the U^x target.

    E_IJ(C, R) = X_I^T A(C) X_J
    dE_IJ/dR   = [dE_IJ/dR]_C  +  sum_pq (dE_IJ/dU_pq) U^x_pq

R_pq = dE_IJ/dU_pq is the INTERSTATE ORBITAL GRADIENT -- the RHS that the new
Z-vector solve needs, and the one thing still missing from the closed-form
amplitude term.  It is measured here by perturbing the MO coefficients
    C[:,q] -> C[:,q] +/- eps * C[:,p]
and finite-differencing E_IJ.  Crucially E_IJ costs ONE matvec application
(A X_J then dot with X_I), not a full 90x90 rebuild, so the whole gradient is
~2*nbf^2 matvec calls.

U^x_pq is measured from the MO overlap between the reference and displaced
geometries (tag OQP::overlap_mo_non_orthogonal), central-differenced.

PASS = sum_pq R_pq U^x_pq reproduces  missing = oracle - ana2e/gap
       (data_snapshots/p12_ux_target.npz) in direction and magnitude.

A PASS means the RHS is identified and the remaining work is mechanical:
write R analytically in Fortran and route it through the existing z-vector
solver + nuclear contraction instead of computing U^x per coordinate.
A FAIL means the operator is still wrong and no Fortran should be written.
"""
import os
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
TGT = '/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p12_ux_target.npz'
OUT = '/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p12_orbgrad.npz'
EPS = 1.0e-4          # MO-mixing step
DX = 1.0e-3           # nuclear step for U^x
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p12_orbgrad.log')
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


def set_bvec(col):
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)


def E_IJ(I, J):
    """X_I^T A X_J with the CURRENT MO coefficients and the frozen AO Fock."""
    set_bvec(X0[J])
    oqp.mrsf_matvec_apply(mol)
    ax = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
    return float(X0[I] @ ax)


def put_C(Ca, Cb):
    mol.data['OQP::VEC_MO_A'] = Ca
    mol.data['OQP::VEC_MO_B'] = Cb


# ------------------------------------------------------------------ R_pq
pairs = [(I, J) for I in range(nstate) for J in range(nstate) if I != J]
R = {p: np.zeros((nbf, nbf)) for p in pairs}
for p in range(nbf):
    for q in range(nbf):
        if p == q:
            continue
        for sgn in (+1, -1):
            Ca = C0.copy()
            Cb = C0b.copy()
            Ca[:, q] += sgn * EPS * C0[:, p]
            Cb[:, q] += sgn * EPS * C0b[:, p]
            put_C(Ca, Cb)
            for (I, J) in pairs:
                R[(I, J)][p, q] += sgn * E_IJ(I, J) / (2.0 * EPS)
put_C(C0, C0b)

# ------------------------------------------------------------------ U^x_pq
cfg = mol.config
json0 = mol.log.replace('.log', '.json')
mol.save_data()
guess_bak = dict(cfg['guess'])
cfg['guess']['type'] = 'json'
cfg['guess']['file'] = json0
cfg['guess']['continue_geom'] = False


def mo_overlap(coord):
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
    return np.array(mol.data[OVTAG], copy=True).reshape(-1).reshape((nbf, nbf)).T


Ux = np.zeros((ncoord, nbf, nbf))
for k in range(ncoord):
    Mp = mo_overlap(xyz0 + DX * np.eye(ncoord)[k])
    Mm = mo_overlap(xyz0 - DX * np.eye(ncoord)[k])
    Ux[k] = (Mp - Mm) / (2.0 * DX)

cfg['guess'].update(guess_bak)
mol.update_system(xyz0)
oqp.library.ints_1e(mol)
oqp.library.guess(mol)
SinglePoint(mol).energy()
mol.data['OQP::td_bvec_mo'] = X0_raw

# ------------------------------------------------------------------ test
tgt = np.load(TGT)
print(f'# {"pair":>7} {"|R.Ux|":>11} {"|missing|":>11} {"cos":>11} {"ratio":>9}')
out = {}
for (I, J) in pairs:
    key = f'missing_{I+1}{J+1}'
    if key not in tgt:
        continue
    gap = Om[J] - Om[I]
    contrib = np.array([np.sum(R[(I, J)] * Ux[k]) for k in range(ncoord)]) / gap
    miss = tgt[key]
    c = contrib @ miss / (np.linalg.norm(contrib) * np.linalg.norm(miss) + 1e-300)
    print(f'  {str((I+1,J+1)):>7} {np.linalg.norm(contrib):11.6f} '
          f'{np.linalg.norm(miss):11.6f} {c:+11.6f} '
          f'{np.linalg.norm(contrib)/(np.linalg.norm(miss)+1e-300):9.4f}')
    out[f'RUx_{I+1}{J+1}'] = contrib
    out[f'R_{I+1}{J+1}'] = R[(I, J)]
out['Ux'] = Ux
np.savez(OUT, **out)
print(f'\nsaved -> {OUT}')
