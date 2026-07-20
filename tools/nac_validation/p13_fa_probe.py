"""Locate the 21.9 discrepancy between the exported fa and C^T F_AO C.

Hypotheses:
  (i)  the ixcore level-shift sets fa[core,core] = -1e6 on some diagonal(s);
       21.9 is way too small for that, so it is NOT the whole story, but check.
  (ii) FOCK_A packing order (row-major lower vs upper) is transposed.
  (iii) orthogonal_transform_sym uses a different convention than C^T F C.
"""
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
r = Runner(input_file=INP,
           log='/bighome/alireza/openqp-nac/tools/nac_validation/p13_faprobe.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
C0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = C0.shape[0]
nij = noca * (nbf - nocb)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((-1, nij))
rr = X0_raw.copy().reshape(-1)
rr[0:nij] = X0[0]
mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)
oqp.mrsf_matvec_apply(mol)
fa0 = np.array(mol.data['OQP::nac_fa'], copy=True).reshape(-1).reshape((nbf, nbf)).T

pk = np.array(mol.data['OQP::FOCK_A'], copy=True).ravel()
FAO = np.zeros((nbf, nbf))
idx = 0
for i in range(nbf):
    for j in range(i + 1):
        FAO[i, j] = FAO[j, i] = pk[idx]
        idx += 1
fcc = C0.T @ FAO @ C0

D = fa0 - fcc
print(f'nbf={nbf} nocb(core+doc)={nocb} noca={noca}')
print(f'|fa0 - C^T FAO C| max = {np.abs(D).max():.4e}, frob = {np.linalg.norm(D):.4e}')
print(f'diag(fa0)  = {np.array2string(np.diag(fa0), precision=3, max_line_width=200)}')
print(f'diag(fcc)  = {np.array2string(np.diag(fcc), precision=3, max_line_width=200)}')
# where are the big differences?
big = np.argwhere(np.abs(D) > 1e-3)
print(f'\n# entries with |diff|>1e-3: {len(big)}')
for (i, j) in big[:20]:
    print(f'  ({i:2d},{j:2d})  fa0={fa0[i,j]:12.5f}  fcc={fcc[i,j]:12.5f}  diff={D[i,j]:12.5f}')
# excluding the core row/col (iter=1 -> index 0)?
mask = np.ones((nbf, nbf), bool)
mask[0, :] = False
mask[:, 0] = False
print(f'\n# |fa0 - C^T FAO C| max EXCLUDING MO 0 (core) = {np.abs(D[mask]).max():.4e}')
mask2 = np.ones((nbf, nbf), bool)
mask2[:nocb, :] = False
mask2[:, :nocb] = False
print(f'# |fa0 - C^T FAO C| max in the VIRTUAL-VIRTUAL block only = '
      f'{np.abs(D[noca:, noca:]).max():.4e}')
print(f'# |fa0 - C^T FAO C| max in the OCC-OCC block (>=nocb) = '
      f'{np.abs(D[nocb:noca, nocb:noca]).max():.4e}')
