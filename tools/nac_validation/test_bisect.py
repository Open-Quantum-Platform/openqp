"""Component bisection: which RHS/gradient piece carries the G-deficiency?"""
import os, math
import oqp
from oqp.pyoqp import Runner
import numpy as np

inp = '/tmp/nactest/H2O_tight_dx0.001.inp'
r = Runner(input_file=inp, log='/tmp/nactest/bis.log')
r.run()
mol = r.mol
nstate = 3
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
nij = X0_raw.size // nstate
X = X0_raw.reshape(-1).reshape((nstate, nij))

def G(col, pygate=None):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(X0_raw.shape)
    mol.data.set_tdhf_target(1)
    oqp.tdhf_mrsf_z_vector(mol)
    if pygate:
        for tag in pygate:
            mol.data[tag] = np.zeros_like(np.array(mol.data[tag], copy=True))
    oqp.tdhf_mrsf_gradient(mol)
    g = mol.get_grad().reshape(-1).copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    return g

def B(a, b, pygate=None):
    u, v = X[a-1], X[b-1]
    return 0.25*(G(u+v, pygate) - G(u-v, pygate))

# known deficiency for (1,3): from qz audit npz
q = np.load('/tmp/nactest/qz_ball7.npz')
Bc13 = q['Z13_coded'] - 0.5*q['X1_coded'] - 0.5*q['X3_coded']
Bt13 = q['Z13_numeric'] - 0.5*q['X1_numeric'] - 0.5*q['X3_numeric']
DEF = Bc13 - Bt13                # cross-process! own-process B_all below fixes phase

gates_f = {'HX': 'NAC_ZERO_HX', 'SP': 'NAC_ZERO_SP', 'T': 'NAC_ZERO_T',
           'AB1': 'NAC_ZERO_AB1'}
B_all = B(1, 3)
# phase: align own-process B_all with cross-process Bc13
s = np.sign(np.dot(B_all, Bc13)) or 1.0
DEF = s * DEF
print(f'|B_all|={np.linalg.norm(B_all):.5f}  own-vs-npz cos={np.dot(B_all, s*Bc13)/(np.linalg.norm(B_all)*np.linalg.norm(Bc13)):.6f}')
print(f'deficiency |DEF|={np.linalg.norm(DEF):.5f}')
print()
for name, ev in gates_f.items():
    os.environ[ev] = '1'
    Bw = B(1, 3)
    del os.environ[ev]
    dBc = B_all - Bw             # component-c bilinear contribution
    if np.linalg.norm(dBc) < 1e-8:
        print(f'{name}: |dB|=0'); continue
    c = np.dot(dBc, DEF)/(np.linalg.norm(dBc)*np.linalg.norm(DEF)+1e-30)
    sc = np.dot(DEF, dBc)/np.dot(dBc, dBc)
    print(f'{name}: |dB_c|={np.linalg.norm(dBc):.5f}  cos(dB_c, DEF)={c:+.6f}  '
          f'DEF~s*dB_c: s={sc:+.4f} resid={np.linalg.norm(sc*dBc-DEF):.5f}')
# python-side gates (direct 2e densities)
for name, tags in (('V(td_abxc)', ['OQP::td_abxc']),
                   ('SPC(td_mrsf_den)', ['OQP::td_mrsf_density']),
                   ('WAO', ['OQP::WAO'])):
    Bw = B(1, 3, pygate=tags)
    dBc = B_all - Bw
    if np.linalg.norm(dBc) < 1e-8:
        print(f'{name}: |dB|=0'); continue
    c = np.dot(dBc, DEF)/(np.linalg.norm(dBc)*np.linalg.norm(DEF)+1e-30)
    sc = np.dot(DEF, dBc)/np.dot(dBc, dBc)
    print(f'{name}: |dB_c|={np.linalg.norm(dBc):.5f}  cos(dB_c, DEF)={c:+.6f}  '
          f'DEF~s*dB_c: s={sc:+.4f} resid={np.linalg.norm(sc*dBc-DEF):.5f}')
