"""Validate the explicit-ERI mrsf_nac_amp Fortran driver against the gauge-free
damp oracle. damp(I,J) = G_IJ / (Om_J - Om_I), G_IJ = OQP::nac_amp(:,I,J).

Run env per the build recipe; OQP_NAC_AMP_NOP2 controls FIX-2 (relaxed p2):
  unset/empty -> use OQP::td_p as the relaxed p2 (whatever the seam left)
  set         -> p2=0 (transition^2-only measurement, the clean first cut)
"""
import os
import math
import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint
import numpy as np

INP = os.environ.get('AMP_INP', '/tmp/nactest/H2O_tight_dx0.001.inp')
ORACLE = os.environ.get('AMP_ORACLE', '/tmp/nactest/p11_damp_oracle.npz')
LOG = '/tmp/nactest/p11_amp_validate.log'

r = Runner(input_file=INP, log=LOG)
r.run()
mol = r.mol

nstate = 3
natom = int(mol.data['natom'])
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]      # excitation energies
print("energies:", np.round(E[:nstate + 1], 6))
print("Om      :", np.round(Om, 6))

# zero td_p if requested so the engine sees p2=0 (transition^2-only)
nop2 = os.environ.get('OQP_NAC_AMP_NOP2', '')
print("OQP_NAC_AMP_NOP2 =", repr(nop2))

oqp.mrsf_nac_amp(mol)

# OQP::nac_amp is the Fortran (3*natom, nstate, nstate) column-major buffer
raw = np.array(mol.data['OQP::nac_amp'], copy=True)
ampF = raw.reshape(-1).reshape((nstate, nstate, 3 * natom))   # [jst, ist, k]
# ampF[jst-1, ist-1, k] = nac_amp(k, ist, jst); swap so G[I-1,J-1] = G_IJ
G = np.transpose(ampF, (1, 0, 2))                              # [I, J, k]

oracle = np.load(ORACLE)
PAIRS = [(1, 2), (1, 3), (2, 3)]
print("\n=== explicit-ERI damp vs oracle (damp = G_IJ / (Om_J - Om_I)) ===")
for (I, J) in PAIRS:
    gij = G[I - 1, J - 1]
    gap = Om[J - 1] - Om[I - 1]
    damp = gij / gap
    o = oracle[f'{I}{J}_damp']
    c = np.dot(damp, o) / (np.linalg.norm(damp) * np.linalg.norm(o) + 1e-30)
    rr = np.linalg.norm(damp) / (np.linalg.norm(o) + 1e-30)
    print(f"  ({I},{J}): gap={gap:+.5f}  |G|={np.linalg.norm(gij):.5f}  "
          f"|damp|={np.linalg.norm(damp):.5f}  |oracle|={np.linalg.norm(o):.5f}")
    print(f"           cos={c:+.6f}  ratio(|damp|/|oracle|)={rr:.4f}")
    s = 1.0 if c >= 0 else -1.0
    print(f"           damp  ={np.array2string(s*damp, precision=5, suppress_small=True)}")
    print(f"           oracle={np.array2string(o, precision=5, suppress_small=True)}")
