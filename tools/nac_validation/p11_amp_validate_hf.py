"""Same explicit-ERI damp validation for the HF case (clean ROHF half-factor
magnitude check on HF (1,3) |damp|=0.036)."""
import os
import oqp
import oqp.library
from oqp.pyoqp import Runner
import numpy as np

INP = '/tmp/nactest/H2O_tightHF_dx0.001.inp'
ORACLE = '/tmp/nactest/p11_damp_oracle_hf.npz'
r = Runner(input_file=INP, log='/tmp/nactest/p11_amp_hf.log')
r.run()
mol = r.mol
nstate = 3
natom = int(mol.data['natom'])
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]
os.environ['OQP_NAC_AMP_NOP2'] = '1'
oqp.mrsf_nac_amp(mol)
raw = np.array(mol.data['OQP::nac_amp'], copy=True)
G = np.transpose(raw.reshape(-1).reshape((nstate, nstate, 3 * natom)), (1, 0, 2))
oracle = np.load(ORACLE)
print("=== HF explicit-ERI 2e damp vs oracle ===")
for (I, J) in [(1, 2), (1, 3), (2, 3)]:
    gap = Om[J - 1] - Om[I - 1]
    damp = G[I - 1, J - 1] / gap
    o = oracle[f'{I}{J}_damp']
    c = np.dot(damp, o) / (np.linalg.norm(damp) * np.linalg.norm(o) + 1e-30)
    rr = np.linalg.norm(damp) / (np.linalg.norm(o) + 1e-30)
    print(f"({I},{J}): |2e|={np.linalg.norm(damp):.5f} |oracle|={np.linalg.norm(o):.5f} "
          f"cos={c:+.5f} ratio={rr:.4f}")
