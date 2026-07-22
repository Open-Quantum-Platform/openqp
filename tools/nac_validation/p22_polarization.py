"""
Closed-form test (from Lee2019 + Zhang-Herbert): the amplitude term is the
POLARIZATION of the analytic MRSF gradient.

The excited-state gradient assembly (Lee2019 Eq 3.21 / Zhang-Herbert Eq A20)
contains NO state eigenvalue Omega -- it is a pure quadratic form G(X)=X^T(dA)X in
the amplitude. Hence
    X_I^T(dA)X_J = 1/2 [ G(X_I+X_J) - G(X_I) - G(X_J) ]
    d_amp(I,J)  = X_I^T(dA)X_J / (Om_J - Om_I).
G(X) = the existing, validated analytic MRSF gradient run with amplitude X.

Sanity 1 (homogeneity): G(c X) should be c^2 G(X) if the assembly is a clean
degree-2 quadratic (the whole premise). Report G(2 X)/G(X).
Sanity 2: G(X_k) for an eigenvector = the state-k gradient dOm_k (finite; nonzero).
Main: polarization d_amp vs the oracle (_compute_amp_damp) and vs the numerical.
"""
import math
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import NAC

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p22.log')
r.run()
mol = r.mol
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
bkey = 'OQP::td_bvec_mo'
raw0 = np.array(mol.data[bkey], copy=True)
state_axis = 0 if raw0.shape[0] == nstate else 1
nij = raw0.shape[1] if state_axis == 0 else raw0.shape[0]
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]
pairs = [(i, j) for i in range(1, nstate + 1) for j in range(1, nstate + 1) if i < j]


def get_state(k):
    return (raw0[k - 1, :] if state_axis == 0 else raw0[:, k - 1]).copy()


def grad_for(target, amp):
    raw = raw0.copy()
    if state_axis == 0:
        raw[target - 1, :] = amp
    else:
        raw[:, target - 1] = amp
    mol.data[bkey] = raw
    mol.data.set_tdhf_target(target)
    oqp.tdhf_mrsf_z_vector(mol)
    if not mol.mol_energy.Z_Vector_converged:
        return None
    oqp.tdhf_mrsf_gradient(mol)
    return mol.get_grad().reshape((natom, 3)).copy()


# ---- oracle (trusted semi-numerical damp) ----
oracle = NAC(mol)._compute_amp_damp(dx=1.0e-3)
mol.data[bkey] = raw0

# ---- Sanity 1: homogeneity degree 2 ----
X1 = get_state(1)
g1 = grad_for(1, X1)
g2 = grad_for(1, 2.0 * X1)
if g1 is not None and g2 is not None:
    ratio = np.linalg.norm(g2) / (np.linalg.norm(g1) + 1e-300)
    print(f'# Sanity homogeneity: |G(2X)|/|G(X)| = {ratio:.4f}  (expect 4.0 if '
          f'G is a clean degree-2 quadratic)')

# ---- ground-state gradient G(0): the constant part E_0^xi ----
g0 = grad_for(1, np.zeros_like(X1))
print(f'# |G(0)| (ground-state gradient, the constant part) = '
      f'{np.linalg.norm(g0) if g0 is not None else float("nan"):.6f}')

# ---- Main: CORRECTED polarization  X_I^T dA X_J = 1/2[G(X+)-G(X_I)-G(X_J)+G(0)] ----
print(f'# {"pair":>7} {"|d_pol|":>11} {"|oracle|":>11} {"cos":>11} {"ratio":>9}')
for (I, J) in pairs:
    gI = grad_for(I, get_state(I))
    gJ = grad_for(J, get_state(J))
    gS = grad_for(I, get_state(I) + get_state(J))   # target I, amp X_I+X_J
    mol.data[bkey] = raw0
    if gI is None or gJ is None or gS is None or g0 is None:
        print(f'  {str((I,J)):>7}  z-vector not converged')
        continue
    inter = 0.5 * (gS - gI - gJ + g0)               # +G(0) removes the constant
    dpol = (inter / (Om[J - 1] - Om[I - 1])).reshape(-1)
    orc = oracle[(I, J)].reshape(-1)
    c = dpol @ orc / (np.linalg.norm(dpol) * np.linalg.norm(orc) + 1e-300)
    print(f'  {str((I,J)):>7} {np.linalg.norm(dpol):11.6f} {np.linalg.norm(orc):11.6f} '
          f'{c:+11.6f} {np.linalg.norm(dpol)/(np.linalg.norm(orc)+1e-300):9.4f}')
