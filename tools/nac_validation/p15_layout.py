"""
Settle the VEC_MO_A layout, which decides whether the 'matvec re-canonicalizes'
diagnosis was real or just a transpose error in my probe.

fa0 exported by the matvec is diag(eps) (canonical). fa0 = mo_a^T F_AO mo_a where
mo_a is the Fortran mo coefficient matrix. If C = np.array(VEC_MO_A) is stored
TRANSPOSED vs mo_a (Fortran column-major, likely), then:
    C @ FAO @ C.T  == diag(eps)   (correct contraction)   [not C.T @ FAO @ C]
and an orbital rotation must be applied as  C -> R @ C  (rotate ROWS of C = the
MO index), NOT C @ R.

Tests:
  1. which of C.T@F@C , C@F@C.T equals diag(eps)?  -> the layout.
  2. with the correct layout + correct rotation side, does the exported fa satisfy
     fa(rot) == (rotation)^T fa0 (rotation) to MACHINE precision? If YES, there is
     NO re-canonicalization and the orbital gradient IS measurable.
"""
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
r = Runner(input_file=INP,
           log='/bighome/alireza/openqp-nac/tools/nac_validation/p15.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
Craw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = Craw.shape[0]
eMO = np.array(mol.data['OQP::E_MO_A'], copy=True).ravel()
nij = noca * (nbf - nocb)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((-1, nij))

# unpack FOCK_A (lower-tri packed) two ways
pk = np.array(mol.data['OQP::FOCK_A'], copy=True).ravel()
FAO = np.zeros((nbf, nbf))
idx = 0
for i in range(nbf):
    for j in range(i + 1):
        FAO[i, j] = FAO[j, i] = pk[idx]
        idx += 1


def fa_export(C):
    mol.data['OQP::VEC_MO_A'] = C
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = X0[0]
    mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)
    oqp.mrsf_matvec_apply(mol)
    raw = np.array(mol.data['OQP::nac_fa'], copy=True).reshape(-1).reshape((nbf, nbf))
    return raw     # decide layout below


D = np.diag(eMO)
opts = {
    'C.T @ F @ C': Craw.T @ FAO @ Craw,
    'C @ F @ C.T': Craw @ FAO @ Craw.T,
}
print('# TEST 1: which contraction == diag(eps)?')
for name, M in opts.items():
    off = np.abs(M - np.diag(np.diag(M))).max()
    de = np.abs(np.diag(M)[1:] - eMO[1:]).max()
    print(f'  {name:>14}:  max off-diag = {off:.3e},  |diag-eps|[1:] = {de:.3e}')

# The layout that diagonalizes tells us mo_a. If C.T@F@C is diagonal, mo_a = C
# (columns are MOs). If C@F@C.T is diagonal, mo_a = C.T (rows are MOs).
col_is_mo = np.abs(opts['C.T @ F @ C'] - np.diag(np.diag(opts['C.T @ F @ C']))).max() < 1e-6
print(f'\n# => columns of VEC_MO_A are MOs: {col_is_mo}')

# nac_fa export layout: does raw or raw.T equal diag(eps)?
fa_raw = fa_export(Craw)
lay = 'raw' if np.abs(fa_raw - np.diag(np.diag(fa_raw))).max() < \
    np.abs(fa_raw.T - np.diag(np.diag(fa_raw.T))).max() else 'raw.T'
fa0 = fa_raw if lay == 'raw' else fa_raw.T
print(f'# nac_fa export layout that is diagonal: {lay}')

# TEST 2: exact rotation transform with the correct side
th = 1.0e-4


def givens(p, q, t):
    R = np.eye(nbf)
    c, s = np.cos(t), np.sin(t)
    R[p, p] = R[q, q] = c
    R[p, q] = -s
    R[q, p] = s
    return R


print('\n# TEST 2: does fa transform as an EXACT rotation? (machine-0 => no re-canon)')
for (p, q) in [(2, 8), (7, 3), (1, 12)]:
    R = givens(p, q, th)
    # rotate MOs: if columns are MOs, C' = C @ R; if rows are MOs, C' = R @ C
    Cp = Craw @ R if col_is_mo else R @ Craw
    fr = fa_export(Cp)
    fr = fr if lay == 'raw' else fr.T
    pred = R.T @ fa0 @ R
    print(f'  (p,q)=({p},{q}): |fa(rot) - R^T fa0 R| = {np.abs(fr - pred).max():.3e}')
mol.data['OQP::VEC_MO_A'] = Craw
