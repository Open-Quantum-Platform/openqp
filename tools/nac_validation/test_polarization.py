import math
import oqp
from oqp.pyoqp import Runner
import numpy as np

inp = '/tmp/nactest/H2O_tight_dx0.001.inp'
r = Runner(input_file=inp, log='/tmp/nactest/pol.log')
r.run()
mol = r.mol
nstate = 3
natom = mol.data['natom']
bkey = 'OQP::td_bvec_mo'
raw0 = np.array(mol.data[bkey], copy=True)
nij = raw0.size // nstate
sqrt2 = math.sqrt(2.0)


def grad_for(target, amp_col):
    raw = raw0.copy()
    flat = raw.reshape(-1)
    flat[(target - 1) * nij:target * nij] = amp_col   # state-contiguous layout
    mol.data[bkey] = flat.reshape(raw0.shape)
    mol.data.set_tdhf_target(target)
    oqp.tdhf_mrsf_z_vector(mol)
    ok = mol.mol_energy.Z_Vector_converged
    oqp.tdhf_mrsf_gradient(mol)
    g = mol.get_grad().reshape(-1).copy()
    mol.data[bkey] = raw0                              # restore
    return g if ok else None


X = raw0.reshape(-1).reshape((nstate, nij))
bench = np.load('/tmp/nactest/benchmark_bhh.npz')
print('\n===== polarization identity vs measured d_amp =====')
for (i, j) in ((1, 2), (1, 3), (2, 3)):
    Xp = (X[i - 1] + X[j - 1]) / sqrt2
    Xm = (X[i - 1] - X[j - 1]) / sqrt2
    gp = grad_for(i, Xp)
    gm = grad_for(i, Xm)
    if gp is None or gm is None:
        print(f'({i},{j}): z-vector failed'); continue
    h = 0.5 * (gp - gm)                                # X_I^T (dA) X_J
    gap = mol.energies[j] - mol.energies[i]
    dci = h / gap
    damp = bench[f'{i}{j}_damp']
    c = np.dot(dci, damp) / (np.linalg.norm(dci) * np.linalg.norm(damp) + 1e-30)
    print(f'({i},{j}): |dci|={np.linalg.norm(dci):.5f}  |damp|={np.linalg.norm(damp):.5f}  '
          f'ratio={np.linalg.norm(dci)/np.linalg.norm(damp):.4f}  |cos|={abs(c):.6f}')
