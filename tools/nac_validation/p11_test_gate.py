"""Recompute coded_cross with NAC_OUTPUT_FOLD gate ON/OFF and compare to numeric_cross
(=gap*oracle, the truth, from qz_full.npz *_numeric).  Decides empirically whether the
output-fold gate closes the coded_cross (1,3) sign flip / (1,2) doc-socc ratio.

coded(col) = Gfull(col) = oqp.tdhf_mrsf_gradient on amplitude col minus SCF gradient
(verbatim test_QZ.py Gfull).  coded_cross(I,J)=coded(Z_IJ)-1/2[coded(X_I)+coded(X_J)].
"""
import os
import math
import oqp
from oqp.pyoqp import Runner
import numpy as np

GATE = os.environ.get('NAC_OUTPUT_FOLD', '')
inp = '/tmp/nactest/H2O_tight_dx0.001.inp'
r = Runner(input_file=inp, log='/tmp/nactest/testgate.log'); r.run(); mol = r.mol
nstate, natom = 3, mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca-2
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True); nbf = C0raw.shape[0]
nvirb = nbf-nocb; nij = noca*nvirb
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
X = X0_raw.reshape(-1).reshape((nstate, nij))


def Gfull(col):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(X0_raw.shape)
    mol.data.set_tdhf_target(1)
    oqp.tdhf_mrsf_z_vector(mol)
    oqp.tdhf_mrsf_gradient(mol)
    gZ = mol.get_grad().reshape(-1).copy()
    for tag in ('OQP::td_p', 'OQP::WAO', 'OQP::td_abxc', 'OQP::td_mrsf_density'):
        mol.data[tag] = np.zeros_like(np.array(mol.data[tag], copy=True))
    oqp.tdhf_mrsf_gradient(mol)
    gS = mol.get_grad().reshape(-1).copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    return gZ - gS


Zs = {'Z12': (X[0]+X[1])/math.sqrt(2), 'Z13': (X[0]+X[2])/math.sqrt(2),
      'Z23': (X[1]+X[2])/math.sqrt(2),
      'X1': X[0].copy(), 'X2': X[1].copy(), 'X3': X[2].copy()}
GZ = {n: Gfull(z) for n, z in Zs.items()}


def cross(d, p, a, b):
    return d[f'Z{p}'] - 0.5*(d[f'X{a}'] + d[f'X{b}'])


qz = np.load('/tmp/nactest/qz_full.npz')
num = {p: cross({k: qz[f'{k}_numeric'] for k in Zs}, p, a, b)
       for p, a, b in [('12', '1', '2'), ('13', '1', '3'), ('23', '2', '3')]}
cod = {p: cross(GZ, p, a, b)
       for p, a, b in [('12', '1', '2'), ('13', '1', '3'), ('23', '2', '3')]}

print(f"\n===== NAC_OUTPUT_FOLD = '{GATE}'  (coded_cross vs numeric_cross) =====")
allpass = True
for p in ('12', '13', '23'):
    a, b = cod[p], num[p]
    c = np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30)
    rr = np.linalg.norm(a)/(np.linalg.norm(b)+1e-30)
    ok = abs(c) > 0.99 and 0.9 < rr < 1.1
    allpass = allpass and ok
    print(f"  ({p}): cos={c:+.6f} ratio={rr:.4f}  {'PASS' if ok else 'FAIL'}")
print(f"VERDICT: {'ALL PASS' if allpass else 'NOT all pass'}")
np.savez(f"/tmp/nactest/qz_gate_{GATE or 'off'}.npz", **{f'{n}_coded': GZ[n] for n in Zs})
