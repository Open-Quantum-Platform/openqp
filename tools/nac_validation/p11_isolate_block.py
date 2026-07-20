"""Localise the coded_cross operator defect: A_coded is symmetric & matches A_true on
the diagonal but differs off-diagonal (block-structured: (1,2) doc-socc ~9.66x, (1,3)
virt ~-0.5x).  Decompose Gfull into its gradient-density contributions and find which
one carries the off-diagonal error by zeroing each between z-vector and gradient:
  td_p          : relaxed CPHF z-vector density (orbital response)
  td_abxc       : unrelaxed difference density (sfdmat, quadratic in folded amp)
  td_mrsf_density: spin-coupled (spc) density
  WAO           : energy-weighted density (Wao)
For each, coded_cross is recomputed and compared to numeric_cross (truth).  The
component whose removal/zeroing aligns coded_cross with numeric_cross IS the defect.
"""
import math
import oqp
from oqp.pyoqp import Runner
import numpy as np

inp = '/tmp/nactest/H2O_tight_dx0.001.inp'
r = Runner(input_file=inp, log='/tmp/nactest/isolate.log'); r.run(); mol = r.mol
nstate = 3
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca-2
nbf = np.array(mol.data['OQP::VEC_MO_A'], copy=True).shape[0]
nvirb = nbf-nocb; nij = noca*nvirb
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
X = X0_raw.reshape(-1).reshape((nstate, nij))
ALLTAGS = ('OQP::td_p', 'OQP::WAO', 'OQP::td_abxc', 'OQP::td_mrsf_density')


def Gfull(col, zero_after_zvec=()):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(X0_raw.shape)
    mol.data.set_tdhf_target(1)
    oqp.tdhf_mrsf_z_vector(mol)
    for tag in zero_after_zvec:
        mol.data[tag] = np.zeros_like(np.array(mol.data[tag], copy=True))
    oqp.tdhf_mrsf_gradient(mol)
    gZ = mol.get_grad().reshape(-1).copy()
    for tag in ALLTAGS:
        mol.data[tag] = np.zeros_like(np.array(mol.data[tag], copy=True))
    oqp.tdhf_mrsf_gradient(mol)
    gS = mol.get_grad().reshape(-1).copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    return gZ - gS


def cross_set(zero_after):
    Zs = {'Z12': (X[0]+X[1])/math.sqrt(2), 'Z13': (X[0]+X[2])/math.sqrt(2),
          'Z23': (X[1]+X[2])/math.sqrt(2), 'X1': X[0], 'X2': X[1], 'X3': X[2]}
    GZ = {n: Gfull(z, zero_after) for n, z in Zs.items()}
    out = {}
    for p, a, b in [('12', '1', '2'), ('13', '1', '3'), ('23', '2', '3')]:
        out[p] = GZ[f'Z{p}'] - 0.5*(GZ[f'X{a}'] + GZ[f'X{b}'])
    return out


qz = np.load('/tmp/nactest/qz_full.npz')
num = {}
for p, a, b in [('12', '1', '2'), ('13', '1', '3'), ('23', '2', '3')]:
    num[p] = qz[f'Z{p}_numeric'] - 0.5*(qz[f'X{a}_numeric'] + qz[f'X{b}_numeric'])

configs = [('FULL', ()), ('zero td_p', ('OQP::td_p',)),
           ('zero td_abxc', ('OQP::td_abxc',)),
           ('zero td_mrsf_density', ('OQP::td_mrsf_density',)),
           ('zero WAO', ('OQP::WAO',))]
for name, za in configs:
    cc = cross_set(za)
    line = f"{name:22}"
    for p in ('12', '13', '23'):
        a, b = cc[p], num[p]
        c = np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30)
        rr = np.linalg.norm(a)/(np.linalg.norm(b)+1e-30)
        line += f"  ({p}) cos={c:+.3f} r={rr:.2f}"
    print(line)
