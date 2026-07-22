"""
Exhaustive layout test: does the exported nac_fa equal mo_a^T FOCK_A mo_a for SOME
consistent (mo_a-layout, fa-layout) pair, to MACHINE precision, and does it
transform as an EXACT rotation? If yes -> no re-canonicalization; the earlier
diagnosis was a transpose error and the orbital gradient IS measurable.

FOCK_A is the ROHF alpha Fock (NOT diagonal in the MO basis), so I compare the
export to mo_a^T FOCK_A mo_a directly, NOT to diag(eps).
"""
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p15c.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
Craw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = Craw.shape[0]
nij = noca * (nbf - nocb)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((-1, nij))
pk = np.array(mol.data['OQP::FOCK_A'], copy=True).ravel()
FAO = np.zeros((nbf, nbf))
idx = 0
for i in range(nbf):
    for j in range(i + 1):
        FAO[i, j] = FAO[j, i] = pk[idx]
        idx += 1


def export_raw(C):
    mol.data['OQP::VEC_MO_A'] = C
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = X0[0]
    mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)
    oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_fa'], copy=True).reshape(-1).reshape((nbf, nbf))


fa_raw0 = export_raw(Craw)

# mo_a candidates and fa-read candidates
mo_opts = {'C': Craw, 'C.T': Craw.T}
fa_opts = {'raw': fa_raw0, 'raw.T': fa_raw0.T}
print('# IDENTITY: |fa_export - mo^T FAO mo| for each (mo layout, fa read):')
best = None
for mn, mo in mo_opts.items():
    ref = mo.T @ FAO @ mo
    for fn, fa in fa_opts.items():
        err = np.abs(fa - ref).max()
        tag = '  <== MACHINE ZERO' if err < 1e-9 else ''
        print(f'    mo={mn:>4}, fa={fn:>5}: {err:.3e}{tag}')
        if err < 1e-9:
            best = (mn, fn, mo)

if best is None:
    print('\n# no exact identity -> the export is NOT mo^T FAO mo '
          '(genuine internal transform). Diagnosis stands.')
    sys.exit()

mn, fn, mo = best
print(f'\n# EXACT identity at mo={mn}, fa={fn}. '
      f'Now test rotation (rotate the MO index of mo by R):')
th = 1.0e-4


def givens(p, q, t):
    R = np.eye(nbf)
    c, s = np.cos(t), np.sin(t)
    R[p, p] = R[q, q] = c
    R[p, q] = -s
    R[q, p] = s
    return R


def fa_of(C):
    raw = export_raw(C)
    return raw if fn == 'raw' else raw.T


fa0 = fa_of(Craw)
for (p, q) in [(2, 8), (7, 3), (1, 12)]:
    R = givens(p, q, th)
    # mo' = mo @ R rotates columns (MO index) of mo.  Map back to the C we set.
    moR = mo @ R
    C_set = moR if mn == 'C' else moR.T
    frot = fa_of(C_set)
    pred = R.T @ fa0 @ R           # mo^T FAO mo -> (moR)^T FAO (moR) = R^T fa0 R
    print(f'  (p,q)=({p},{q}): |fa(rot) - R^T fa0 R| = {np.abs(frot - pred).max():.3e}')
mol.data['OQP::VEC_MO_A'] = Craw
