"""
Decisive degree-2 test of the polarization premise.
mrsf_nac_polarize under env NAC_HOMOG dumps G(sX) for s=0,1,2,3 (istate dir)
into OQP::nac_homog (3*natom, 4). If G(X) = c + L(X) + X^T(dA)X (degree<=2),
then with D_s = G(sX)-G(0):
    D1 = L + Q,  D2 = 2L + 4Q,  D3 = 3L + 9Q
    (D2 - 2 D1) = 2Q ,  (D3 - 3 D1) = 6Q  ->  ratio == 3.0 componentwise.
A componentwise ratio != 3 proves a degree>2 / non-polynomial term -> the
ground-config channel breaks polarization exactly. Runs for each state.
"""
import os
import sys
import numpy as np
import oqp
import oqp.library
from oqp.pyoqp import Runner

INP = sys.argv[1] if len(sys.argv) > 1 else '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
os.environ['NAC_HOMOG'] = '1'
r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p25.log')
r.run()
mol = r.mol
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']

print('# Homogeneity (degree-2) test of the MRSF gradient assembly G(sX)')
print('# state   |D1|        |D2|        |D3|     |(D3-3D1)|  |(D2-2D1)|   ratio(med)  ratio(max-comp)')
for st in range(1, nstate + 1):
    oqp.mrsf_nac_polarize(mol, st, st)   # istate==jstate: dir = state st
    raw = np.array(mol.data['OQP::nac_homog'], copy=True)      # (3*natom,4) F-order buffer
    G = raw.reshape(-1).reshape((4, 3 * natom)).T              # columns s=0,1,2,3
    G0, G1, G2, G3 = G[:, 0], G[:, 1], G[:, 2], G[:, 3]
    D1, D2, D3 = G1 - G0, G2 - G0, G3 - G0
    num = D3 - 3.0 * D1      # = 6Q if degree<=2
    den = D2 - 2.0 * D1      # = 2Q if degree<=2
    m = np.abs(den) > 1e-9 * (np.abs(D1).max() + 1e-30)
    comp = num[m] / den[m]   # each == 3.0 if degree-2
    ratio_med = np.median(comp) if comp.size else float('nan')
    ratio_max = comp[np.argmax(np.abs(comp - 3.0))] if comp.size else float('nan')
    print(f'  {st:3d}   {np.linalg.norm(D1):10.6f}  {np.linalg.norm(D2):10.6f}  '
          f'{np.linalg.norm(D3):10.6f}  {np.linalg.norm(num):9.5f}  '
          f'{np.linalg.norm(den):9.5f}   {ratio_med:9.5f}   {ratio_max:9.5f}')
print('# EXPECT ratio == 3.00000 for every state if G is a clean degree-2 form.')
print('# A state whose ratio deviates from 3 is the one whose polarization is deficient.')
