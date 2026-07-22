"""
Closed-form NAC via the Fortran polarization driver mrsf_nac_polarize.
Computes X_i^T(dA/dR)X_j = 1/2[G(X_i+X_j)-G(X_i)-G(X_j)+G(0)] entirely in Fortran
(amplitude injected in the target column, no Python reload), then d_amp = that/DOm.
Compare to the trusted oracle (_compute_amp_damp) and print cos/ratio.
"""
import sys
import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import NAC
import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p24.log')
r.run()
mol = r.mol
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]
pairs = [(i, j) for i in range(1, nstate + 1) for j in range(1, nstate + 1) if i < j]

# oracle first (it perturbs geometry; do before the polarize calls)
oracle = NAC(mol)._compute_amp_damp(dx=1.0e-3)

print(f'# Fortran-polarization closed-form d_amp vs oracle')
print(f'# {"pair":>7} {"|d_polar|":>11} {"|oracle|":>11} {"cos":>11} {"ratio":>9}')
for (i, j) in pairs:
    oqp.mrsf_nac_polarize(mol, i, j)
    raw = np.array(mol.data['OQP::nac_amp_polar'], copy=True)   # (3,natom) F-order
    hij = raw.reshape(-1).reshape((natom, 3))                   # X_i^T dA X_j
    dpol = (hij / (Om[j - 1] - Om[i - 1])).reshape(-1)
    orc = oracle[(i, j)].reshape(-1)
    c = dpol @ orc / (np.linalg.norm(dpol) * np.linalg.norm(orc) + 1e-300)
    print(f'  {str((i,j)):>7} {np.linalg.norm(dpol):11.6f} {np.linalg.norm(orc):11.6f} '
          f'{c:+11.7f} {np.linalg.norm(dpol)/(np.linalg.norm(orc)+1e-300):9.5f}')
