"""Parallelogram test: is G(X) (z-vector+gradient) purely quadratic in X?"""
import math
import oqp
from oqp.pyoqp import Runner
import numpy as np

inp = '/tmp/nactest/H2O_tight_dx0.001.inp'
r = Runner(input_file=inp, log='/tmp/nactest/quad.log')
r.run()
mol = r.mol
nstate = 3
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
nij = X0_raw.size // nstate
X = X0_raw.reshape(-1).reshape((nstate, nij))
s2 = math.sqrt(2.0)

def G(col, target=1):
    raw = X0_raw.copy().reshape(-1)
    raw[(target-1)*nij:target*nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(X0_raw.shape)
    mol.data.set_tdhf_target(target)
    oqp.tdhf_mrsf_z_vector(mol)
    ok = mol.mol_energy.Z_Vector_converged
    oqp.tdhf_mrsf_gradient(mol)
    g = mol.get_grad().reshape(-1).copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    return g if ok else None

for (a, b, name) in ((0, 2, 'X1,X3 (ground/generic)'),
                     (1, 2, 'X2,X3 (ijlr/generic)'),
                     (0, 1, 'X1,X2 (ground/ijlr)')):
    Xa, Xb = X[a], X[b]
    lhs = G(Xa + Xb) + G(Xa - Xb)
    rhs = G(s2 * Xa) + G(s2 * Xb)
    diff = np.linalg.norm(lhs - rhs)
    scale = np.linalg.norm(lhs - 2*G(Xa) - 2*G(Xb) + 2*G(np.zeros_like(Xa)) if False else lhs)
    print(f'{name}: |parallelogram violation| = {diff:.6e}   (|LHS|={np.linalg.norm(lhs):.4f})')
