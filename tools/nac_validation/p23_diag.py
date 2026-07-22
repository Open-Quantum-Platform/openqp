"""
Decisive diagnostic: how does the analytic gradient use the INJECTED amplitude?
Determines whether polarization can work at all.

  gA = grad_for(1, X_1)   gB = grad_for(2, X_2)   (the two state gradients)
  gC = grad_for(1, X_2)   (inject state-2 amplitude into target slot 1)
  -> if gC == gB : the gradient is a pure function of the injected amplitude
     (target index irrelevant) => Omega^xi(X) is well-defined for any X.
  -> if gC == gA : the code uses STORED state-1 data, ignores the injection
     => polarization impossible, need the direct interstate builder.

Also: quadratic check on Omega^xi(X) = G(X) - G(0):
  |Omega^xi(2 X_1)| / |Omega^xi(X_1)| should be 4 if Omega^xi is degree-2.
"""
import sys
import oqp
import oqp.library
from oqp.pyoqp import Runner
import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p23.log')
r.run()
mol = r.mol
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
bkey = 'OQP::td_bvec_mo'
raw0 = np.array(mol.data[bkey], copy=True)
state_axis = 0 if raw0.shape[0] == nstate else 1


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


X1, X2 = get_state(1), get_state(2)
gA = grad_for(1, X1)
gB = grad_for(2, X2)
gC = grad_for(1, X2)     # inject X2 into slot 1
g0 = grad_for(1, np.zeros_like(X1))
gD = grad_for(1, 2.0 * X1)
mol.data[bkey] = raw0


def nrm(g):
    return np.linalg.norm(g) if g is not None else float('nan')


print(f'# gA=grad(1,X1) |{nrm(gA):.5f}|   gB=grad(2,X2) |{nrm(gB):.5f}|   '
      f'gC=grad(1,X2) |{nrm(gC):.5f}|   g0=grad(1,0) |{nrm(g0):.5f}|')
if gC is not None and gB is not None and gA is not None:
    print(f'# |gC - gB| = {np.linalg.norm(gC-gB):.3e}   (0 => injection PURE, '
          f'target index irrelevant)')
    print(f'# |gC - gA| = {np.linalg.norm(gC-gA):.3e}   (0 => code uses STORED '
          f'state-1 data, injection IGNORED)')
if gD is not None and gA is not None and g0 is not None:
    o1 = gA - g0
    o2 = gD - g0
    print(f'# Omega^xi quadratic check: |O(2X1)|/|O(X1)| = '
          f'{np.linalg.norm(o2)/(np.linalg.norm(o1)+1e-300):.4f}  (expect 4.0)')
    print(f'# cos(O(2X1), O(X1)) = '
          f'{o1.ravel()@o2.ravel()/(np.linalg.norm(o1)*np.linalg.norm(o2)+1e-300):+.4f}')
