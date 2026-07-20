"""HYPOTHESIS (i): the off-diagonal coded_cross defect is the AB1/B-matrix term.

coded(col) = Gfull(col) = oqp.tdhf_mrsf_gradient(mol) on amplitude col, minus the
SCF gradient. coded_cross(I,J) = coded(Z_IJ) - 1/2[coded(X_I)+coded(X_J)].

We recompute coded WITH the existing runtime gate NAC_ZERO_AB1=1 (zeros ab1_mo_a/b
in the z-vector RHS, tdhf_mrsf_z_vector.F90:1752-1755) -- no rebuild. Then build
coded_cross_noAB1 and compare cos/ratio to numeric_cross (from qz_full.npz) and to
the oracle (gap*oracle == numeric_cross) for all 3 pairs.

DECISIVE: if NAC_ZERO_AB1 brings coded_cross to cos>0.99 ratio~1 on all pairs
(esp. fixes the (1,3) sign flip), AB1 IS the cause.
"""
import os
os.environ["NAC_ZERO_AB1"] = "1"          # gate BEFORE any oqp gradient call
import math
import oqp, oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint
import numpy as np

inp = '/tmp/nactest/H2O_tight_dx0.001.inp'
r = Runner(input_file=inp, log='/tmp/nactest/ab1iso.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate, natom = 3, mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True); nbf = C0raw.shape[0]
nvirb = nbf - nocb; nij = noca*nvirb
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
X = X0_raw.reshape(-1).reshape((nstate, nij))


def Gfull(col):
    """coded dOmega(Z): full production gradient chain on amplitude col, minus SCF.
    With NAC_ZERO_AB1=1 in env, the ab1_mo term is zeroed inside the z-vector RHS."""
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


Zs = {'Z13': (X[0]+X[2])/math.sqrt(2), 'Z23': (X[1]+X[2])/math.sqrt(2),
      'Z12': (X[0]+X[1])/math.sqrt(2),
      'X1': X[0].copy(), 'X2': X[1].copy(), 'X3': X[2].copy()}
GZ = {n: Gfull(z) for n, z in Zs.items()}
np.savez('/tmp/nactest/qz_noab1.npz', **{f'{n}_coded': GZ[n] for n in Zs})

# ---- production S2 gradient anchor self-check (col = X2) ----
# O dE/dZ self anchor: coded(X2) should equal the production excited gradient minus
# SCF. The anchor value -0.182993575 is the production S2 gradient component; here
# we just report coded(X2) norm and the dominant comp so the parent can sanity-check
# AB1=0 hasn't broken the chain catastrophically.
print(f"coded(X2) sample [first 3] = {GZ['X2'][:3]}")

# ---- baseline coded (AB1 ON) from qz_full.npz, numeric_cross, oracle ----
base = np.load('/tmp/nactest/qz_full.npz')
oracle = np.load('/tmp/nactest/p11_damp_oracle.npz')


def cross(coded_dict, i, j):
    return coded_dict[f'Z{i}{j}'] - 0.5*(coded_dict[f'X{i}'] + coded_dict[f'X{j}'])


coded_noab1 = {'Z12': GZ['Z12'], 'Z13': GZ['Z13'], 'Z23': GZ['Z23'],
               'X1': GZ['X1'], 'X2': GZ['X2'], 'X3': GZ['X3']}
coded_base = {k: base[f'{k}_coded'] for k in coded_noab1}
numeric = {k: base[f'{k}_numeric'] for k in coded_noab1}


def stats(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    cos = np.dot(a, b)/(na*nb+1e-30)
    return cos, na/(nb+1e-30)


print("\n=== coded_cross with NAC_ZERO_AB1=1 vs numeric_cross / oracle ===")
PAIRS = [(1, 2), (1, 3), (2, 3)]
for (i, j) in PAIRS:
    cc_noab1 = cross(coded_noab1, i, j)
    cc_base = cross(coded_base, i, j)
    nc = cross(numeric, i, j)              # numeric_cross == truth
    o = oracle[f'{i}{j}_damp']             # gap*oracle == numeric_cross
    ab1_only = cc_base - cc_noab1          # the AB1 contribution to the cross
    c_nb, r_nb = stats(cc_noab1, nc)
    c_b, r_b = stats(cc_base, nc)
    c_o, r_o = stats(cc_noab1, o)
    ok = abs(c_nb) > 0.99 and 0.9 < r_nb < 1.1
    print(f"\n  PAIR ({i},{j}):")
    print(f"    BASELINE  coded(AB1 on)  vs numeric: cos={c_b:+.5f} ratio={r_b:.4f}")
    print(f"    NAC_ZERO_AB1 coded       vs numeric: cos={c_nb:+.5f} ratio={r_nb:.4f}  {'PASS' if ok else 'FAIL'}")
    print(f"    NAC_ZERO_AB1 coded       vs oracle : cos={c_o:+.5f} ratio={r_o:.4f}")
    print(f"    |cc_noAB1|={np.linalg.norm(cc_noab1):.5f} |cc_base|={np.linalg.norm(cc_base):.5f} |numeric|={np.linalg.norm(nc):.5f}")
    print(f"    AB1-only (base-noAB1) |.|={np.linalg.norm(ab1_only):.5f}  cos(AB1only, numeric)={stats(ab1_only,nc)[0]:+.5f}")
    print(f"    coded_noAB1 = {np.array2string(cc_noab1, precision=5, suppress_small=True)}")
    print(f"    numeric     = {np.array2string(nc, precision=5, suppress_small=True)}")
    print(f"    AB1-only    = {np.array2string(ab1_only, precision=5, suppress_small=True)}")
print("\nsaved qz_noab1.npz")
