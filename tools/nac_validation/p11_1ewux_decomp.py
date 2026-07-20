"""BP#1: characterize the MISSING (1e/W + Ux) damp term = oracle - analytic_2e.
oracle (= full transported-matvec FD damp, p11_poc_gaugefree closes it to FD floor)
is loaded from p11_damp_oracle.npz; analytic 2e from Fortran mrsf_nac_amp (NOP2).
Reports per pair: |miss|, cos(miss, -ana2e) [does it cancel the 2e?],
cos(miss, oracle), and component patterns. No SCF FD -> no crash.
"""
import os
import oqp
import oqp.library
from oqp.pyoqp import Runner
import numpy as np

INP = os.environ.get('AMP_INP', '/tmp/nactest/H2O_tight_dx0.001.inp')
ORACLE = os.environ.get('AMP_ORACLE', '/tmp/nactest/p11_damp_oracle.npz')

r = Runner(input_file=INP, log='/tmp/nactest/p11_1ewux.log'); r.run(); mol = r.mol
nstate = 3
natom = int(mol.data['natom']); ncoord = 3*natom
E = list(mol.energies); Om = [E[k+1]-E[0] for k in range(nstate)]

os.environ['OQP_NAC_AMP_NOP2'] = '1'
oqp.mrsf_nac_amp(mol)
raw = np.array(mol.data['OQP::nac_amp'], copy=True)
ampF = raw.reshape(-1).reshape((nstate, nstate, ncoord))
Gana = np.transpose(ampF, (1, 0, 2))   # [I,J,k] = G_IJ

oracle = np.load(ORACLE)
def cs(a, b): return float(np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30))
print("\n=== MISSING term = oracle - analytic_2e (per pair, 9-comp) ===")
for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    gap = Om[J]-Om[I]
    ana2e = Gana[I, J]/gap
    o = oracle[f'{I+1}{J+1}_damp']
    miss = o - ana2e
    print(f"\n({I+1},{J+1}): gap={gap:+.5f}")
    print(f"  |oracle|={np.linalg.norm(o):.5f}  |ana2e|={np.linalg.norm(ana2e):.5f}  cos(ana2e,oracle)={cs(ana2e,o):+.5f}  ratio={np.linalg.norm(ana2e)/(np.linalg.norm(o)+1e-30):.3f}")
    print(f"  |miss|={np.linalg.norm(miss):.5f}  cos(miss,oracle)={cs(miss,o):+.5f}  cos(miss,-ana2e)={cs(miss,-ana2e):+.5f}  cos(miss,+ana2e)={cs(miss,ana2e):+.5f}")
    print(f"  oracle={np.array2string(o, precision=5, suppress_small=True)}")
    print(f"  ana2e ={np.array2string(ana2e, precision=5, suppress_small=True)}")
    print(f"  miss  ={np.array2string(miss, precision=5, suppress_small=True)}")
