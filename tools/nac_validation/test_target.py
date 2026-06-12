"""Does B_code depend on WHICH state column hosts X+-? (eigenvector-identity probe)"""
import math
import oqp
from oqp.pyoqp import Runner
import numpy as np

inp = '/tmp/nactest/H2O_tight_dx0.001.inp'
r = Runner(input_file=inp, log='/tmp/nactest/tgt.log')
r.run()
mol = r.mol
nstate = 3
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
nij = X0_raw.size // nstate
X = X0_raw.reshape(-1).reshape((nstate, nij))

def G(col, target):
    raw = X0_raw.copy().reshape(-1)
    raw[(target-1)*nij:target*nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(X0_raw.shape)
    mol.data.set_tdhf_target(target)
    oqp.tdhf_mrsf_z_vector(mol)
    oqp.tdhf_mrsf_gradient(mol)
    g = mol.get_grad().reshape(-1).copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    return g

def B(a, b, target):
    u, v = X[a-1], X[b-1]
    return 0.25*(G(u+v, target) - G(u-v, target))

q = np.load('/tmp/nactest/qz_ball7.npz')
for (a, b, z) in ((1,3,'Z13'), (1,2,'Z12'), (2,3,'Z23')):
    Bt = q[f'{z}_numeric'] - 0.5*q[f'X{a}_numeric'] - 0.5*q[f'X{b}_numeric']
    for tgt in (1, 2, 3):
        Bc = B(a, b, tgt)
        s = 1.0 if np.dot(Bc, Bt) >= 0 else -1.0
        c = np.dot(Bc, Bt)/(np.linalg.norm(Bc)*np.linalg.norm(Bt)+1e-30)
        sc = np.dot(Bc, Bt)/np.dot(Bt, Bt)
        print(f'B({a},{b}) target={tgt}: |B|={np.linalg.norm(Bc):.5f} '
              f'cos(B,Btrue)={c:+.6f} scalar={sc:+.4f}')
    print()
