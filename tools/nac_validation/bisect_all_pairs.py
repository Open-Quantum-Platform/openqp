"""DECISIVE: which gradient-chain component carries the off-diagonal NAC deficiency?
Adversarial check (workflow 2026-06-16): the deficiency D=cross_coded-cross_numeric
is the error of the FULL chain incl. AB1 (the (A+B)/B-matrix term, ABSENT from the TDA
matvec, ~100-200% of the RHS magnitude). It has NEVER been gated out, so attributing D
to the 2e back-transform is unproven. This gates each component and correlates its
contribution dB_c = B_all - B_gated with D for ALL pairs ((2,3)~0 = control).

If cos(AB1, D)~+/-1 & resid~0 for (1,3) -> AB1 is the culprit (design's 2e theory wrong).
If cos(SP/SPC, D)~+/-1 & cos(AB1,D)~0 -> 2e/spin-pair confirmed.
"""
import os
import oqp
from oqp.pyoqp import Runner
import numpy as np

r = Runner(input_file='/tmp/nactest/H2O_tight_dx0.001.inp', log='/tmp/nactest/bisall.log')
r.run(); mol = r.mol
nstate = 3
X0 = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
nij = X0.size // nstate
X = X0.reshape(-1).reshape((nstate, nij))
qz = np.load('/tmp/nactest/qz_nogate.npz')


def G(col, fgate=None, pygate=None):
    raw = X0.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(X0.shape)
    mol.data.set_tdhf_target(1)
    for k in ('NAC_ZERO_AB1', 'NAC_ZERO_HX', 'NAC_ZERO_SP', 'NAC_ZERO_T'):
        os.environ.pop(k, None)
    if fgate:
        os.environ[fgate] = '1'
    oqp.tdhf_mrsf_z_vector(mol)
    if fgate:
        os.environ.pop(fgate, None)
    if pygate:
        for tag in pygate:
            mol.data[tag] = np.zeros_like(np.array(mol.data[tag], copy=True))
    oqp.tdhf_mrsf_gradient(mol)
    g = mol.get_grad().reshape(-1).copy()
    mol.data['OQP::td_bvec_mo'] = X0
    return g


def Bcross(a, b, fgate=None, pygate=None):
    u, v = X[a-1], X[b-1]
    return 0.25*(G(u+v, fgate, pygate) - G(u-v, fgate, pygate))


def qzc(I, J, kind):
    return qz[f'Z{I}{J}_{kind}'] - 0.5*(qz[f'X{I}_{kind}'] + qz[f'X{J}_{kind}'])


fgates = [('AB1', 'NAC_ZERO_AB1'), ('HX', 'NAC_ZERO_HX'),
          ('SP', 'NAC_ZERO_SP'), ('T', 'NAC_ZERO_T')]
pygates = [('V_td_abxc', ['OQP::td_abxc']), ('SPC_mrsf_den', ['OQP::td_mrsf_density']),
           ('WAO', ['OQP::WAO'])]

for (I, J) in [(1, 2), (1, 3), (2, 3)]:
    B_all = Bcross(I, J)
    cod, num = qzc(I, J, 'coded'), qzc(I, J, 'numeric')
    s = np.sign(np.dot(B_all, cod)) or 1.0          # align this-process phase to qz
    DEF = s*(cod - num)
    nDEF = np.linalg.norm(DEF)
    print(f"\n===== pair ({I},{J}):  |B_all|={np.linalg.norm(B_all):.5f}  "
          f"own-vs-qz cos={np.dot(B_all, s*cod)/(np.linalg.norm(B_all)*np.linalg.norm(cod)+1e-30):+.4f}  "
          f"|DEF|={nDEF:.5f} =====")
    rows = []
    for nm, ev in fgates:
        dB = B_all - Bcross(I, J, fgate=ev)
        rows.append((nm, dB))
    for nm, tags in pygates:
        dB = B_all - Bcross(I, J, pygate=tags)
        rows.append((nm, dB))
    for nm, dB in rows:
        nd = np.linalg.norm(dB)
        if nd < 1e-8:
            print(f"  {nm:14} |dB|=0"); continue
        cos = np.dot(dB, DEF)/(nd*nDEF + 1e-30)
        sc = np.dot(DEF, dB)/np.dot(dB, dB)
        resid = np.linalg.norm(sc*dB - DEF)/(nDEF + 1e-30)
        print(f"  {nm:14} |dB|={nd:.5f}  cos(dB,DEF)={cos:+.4f}  scale={sc:+.3f}  "
              f"resid/|DEF|={resid:.3f}")
print("\nDECISIVE: a component with cos~+/-1 AND resid/|DEF|~0 on (1,2) AND (1,3) IS the")
print("deficiency. (2,3) |DEF| should be ~0 (control). AB1 high => design 2e theory wrong.")
