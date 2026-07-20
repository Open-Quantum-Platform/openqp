"""Internal consistency: the bilinear engine G_IJ must satisfy the polarization
identity  G_IJ = 1/2 [ G(X_I+X_J) - G(X_I) - G(X_J) ]  where G(X)=G_KK is the
engine's own diagonal (validated bit-for-bit vs production). This confirms the
off-diagonal bilinear 1/2[f_I g_J + f_J g_I] is the correct symmetric extension
of the explicit-ERI 2e gradient X^T dA_2e X (independent of any oracle).

We exercise it by overwriting td_bvec_mo with pseudo-state combinations and
re-calling mrsf_nac_amp at I=J (diagonal), comparing to the true off-diagonal.
"""
import os
import oqp
import oqp.library
from oqp.pyoqp import Runner
import numpy as np

INP = '/tmp/nactest/H2O_tight_dx0.001.inp'
r = Runner(input_file=INP, log='/tmp/nactest/p11_amp_pol.log')
r.run()
mol = r.mol
nstate = 3
natom = int(mol.data['natom'])
os.environ['OQP_NAC_AMP_NOP2'] = '1'

bkey = 'OQP::td_bvec_mo'
raw0 = np.array(mol.data[bkey], copy=True)
shape = raw0.shape
# state axis
if shape[0] == nstate:
    sa = 0
elif shape[1] == nstate:
    sa = 1
else:
    raise SystemExit('cannot find state axis')


def get_state(k):
    return (raw0[k - 1] if sa == 0 else raw0[:, k - 1]).copy()


def run_amp():
    oqp.mrsf_nac_amp(mol)
    raw = np.array(mol.data['OQP::nac_amp'], copy=True)
    ampF = raw.reshape(-1).reshape((nstate, nstate, 3 * natom))
    return np.transpose(ampF, (1, 0, 2))


def set_diag(vec):
    """put vec in ALL state columns so G_IJ for any I!=J is the diagonal of vec."""
    raw = raw0.copy()
    for k in range(nstate):
        if sa == 0:
            raw[k] = vec
        else:
            raw[:, k] = vec
    mol.data[bkey] = raw


# true off-diagonal G_12 from the real amplitudes
mol.data[bkey] = raw0
Gtrue = run_amp()

x1 = get_state(1)
x2 = get_state(2)
sq = np.sqrt(2.0)

# diagonal engine evals via pseudo-states (the engine builds spcI=spcJ at I=J
# only if both columns equal; here all columns equal => any (I,J) is diagonal)
set_diag((x1 + x2) / sq); Gp = run_amp()[0, 1]
set_diag((x1 - x2) / sq); Gm = run_amp()[0, 1]
set_diag(x1); G11 = run_amp()[0, 1]
set_diag(x2); G22 = run_amp()[0, 1]

# polarization extraction of the bilinear from diagonals
G12_pol = 0.5 * (Gp - Gm)                      # = X1^T dA2e X2 (cross only)
G12_engine = Gtrue[0, 1]

c = np.dot(G12_pol, G12_engine) / (np.linalg.norm(G12_pol) * np.linalg.norm(G12_engine) + 1e-30)
print("=== bilinear engine self-consistency (explicit-ERI, no oracle) ===")
print(f"  |G12 engine|={np.linalg.norm(G12_engine):.6f}  |G12 polarization|={np.linalg.norm(G12_pol):.6f}")
print(f"  cos(engine, polarization)={c:+.8f}  max|diff|={np.max(np.abs(G12_engine - G12_pol)):.2e}")
mol.data[bkey] = raw0
