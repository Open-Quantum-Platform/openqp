"""Sanity check on the Fock-scaling handle before trusting anything built on it.

E_II = X_I^T A X_I must equal Omega_I.  Scaling the frozen AO Fock must therefore
move E_II by exactly the esum content of Omega_I.  If E_II does not move at all,
the OQP::FOCK_A/B tag write is not reaching the matvec and the whole handle is
void; if it moves, the handle is real and the off-diagonal E_esum(I!=J)=0 is a
genuine property, not an artefact.
"""
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
r = Runner(input_file=INP,
           log='/bighome/alireza/openqp-nac/tools/nac_validation/p13_check.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20

nstate = mol.config['tdhf']['nstate']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
nbf = np.array(mol.data['OQP::VEC_MO_A']).shape[0]
nij = noca * (nbf - nocb)
F0a = np.array(mol.data['OQP::FOCK_A'], copy=True)
F0b = np.array(mol.data['OQP::FOCK_B'], copy=True)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]


def Eij(I, J, s):
    mol.data['OQP::FOCK_A'] = F0a * s
    mol.data['OQP::FOCK_B'] = F0b * s
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = X0[J]
    mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)
    oqp.mrsf_matvec_apply(mol)
    ax = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
    return float(X0[I] @ ax)


print(f'# {"quantity":>12} {"s=1.0":>16} {"s=1.001":>16} {"delta":>14}')
for I in range(nstate):
    a, b = Eij(I, I, 1.0), Eij(I, I, 1.001)
    print(f'  E_{I+1}{I+1}{"":>8} {a:16.9f} {b:16.9f} {b-a:14.3e}   '
          f'(Omega_{I+1} = {Om[I]:.9f})')
for I in range(nstate):
    for J in range(nstate):
        if I >= J:
            continue
        a, b = Eij(I, J, 1.0), Eij(I, J, 1.001)
        print(f'  E_{I+1}{J+1}{"":>8} {a:16.9f} {b:16.9f} {b-a:14.3e}')
mol.data['OQP::FOCK_A'] = F0a
mol.data['OQP::FOCK_B'] = F0b
