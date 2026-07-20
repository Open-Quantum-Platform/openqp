"""Decompose the full damp oracle into (a) the explicit-ERI 2e part captured by
the Fortran mrsf_nac_amp engine, and (b) the residual that the 1e/W +
orbital-response (Ux) terms must supply.

   oracle damp = X_I^T dA X_J / gap            (full matvec, the PoC/benchmark)
   G_fortran   = X_I^T dA_2e^explicit X_J      (my engine, ERI-deriv, frozen C)
   residual    = oracle - G_fortran/gap

Report per-pair cos/ratio of G_fortran vs oracle, and the residual magnitude
(= how much a separate 1e/W + Ux term must contribute).
"""
import os
import oqp
import oqp.library
from oqp.pyoqp import Runner
import numpy as np

INP = '/tmp/nactest/H2O_tight_dx0.001.inp'
ORACLE = '/tmp/nactest/p11_damp_oracle.npz'

r = Runner(input_file=INP, log='/tmp/nactest/p11_amp_decomp.log')
r.run()
mol = r.mol
nstate = 3
natom = int(mol.data['natom'])
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]

os.environ['OQP_NAC_AMP_NOP2'] = '1'
oqp.mrsf_nac_amp(mol)
raw = np.array(mol.data['OQP::nac_amp'], copy=True)
ampF = raw.reshape(-1).reshape((nstate, nstate, 3 * natom))
G = np.transpose(ampF, (1, 0, 2))

oracle = np.load(ORACLE)
print("=== explicit-ERI 2e damp vs full oracle; residual = 1e/W + Ux ===")
for (I, J) in [(1, 2), (1, 3), (2, 3)]:
    gap = Om[J - 1] - Om[I - 1]
    damp2e = G[I - 1, J - 1] / gap
    o = oracle[f'{I}{J}_damp']
    c = np.dot(damp2e, o) / (np.linalg.norm(damp2e) * np.linalg.norm(o) + 1e-30)
    rr = np.linalg.norm(damp2e) / (np.linalg.norm(o) + 1e-30)
    resid = o - damp2e
    print(f"({I},{J}): |2e|={np.linalg.norm(damp2e):.5f} |oracle|={np.linalg.norm(o):.5f} "
          f"cos={c:+.5f} ratio={rr:.4f} |residual|={np.linalg.norm(resid):.5f} "
          f"(={100*np.linalg.norm(resid)/np.linalg.norm(o):.0f}% of oracle)")
