"""Phase 11 SYNTHESIS: push the antisymmetric interstate sp cross through the
mrsf_nac_cphf nuclear seam and compare its DIRECTION to the (1,3) nuclear
deficiency D = coded - numeric (qz_nogate). Confounded in magnitude (seam is
orbital-response only) but DIRECTION is the discriminator: if cos(seam(anti), D)
is high for (1,3) and the (1,3) deficiency is the antisymmetric H-H pattern, the
antisym sp piece is the missing physics.
"""
import os
import oqp
import oqp.library
from oqp.pyoqp import Runner
import numpy as np

INP = '/tmp/nactest/H2O_en.inp'
SQ = 1.0 / np.sqrt(2.0)
r = Runner(input_file=INP, log='/tmp/nactest/p11_an.log'); r.run(); mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
moa = np.array(mol.data['OQP::VEC_MO_A'], copy=True); nbf = moa.shape[0]
nvirb = nbf - nocb; nij = noca * nvirb; nstate = 3; natom = mol.data['natom']
Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); Xshape = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))
prs = ([(i, j) for i in range(nocb, noca) for j in range(0, nocb)]
       + [(k, j) for k in range(noca, nbf) for j in range(0, nocb)]
       + [(k, i) for k in range(noca, nbf) for i in range(nocb, noca)])
nconf = len(prs)

anti = np.load('/tmp/nactest/p11_sp_antisym.npz')
mol.save_data()

def inject_wrk1(R):
    W = np.zeros((nbf, nbf))
    for n, (hi, lo) in enumerate(prs):
        W[hi, lo] = R[n]
    return W

def push_gamma(W):
    flat = np.zeros(nbf * nbf * nstate * nstate)
    for I in range(nstate):
        for J in range(nstate):
            flat[(I + J * nstate) * nbf * nbf:(I + J * nstate + 1) * nbf * nbf] = W.T.reshape(-1)
    mol.data['OQP::nac_gamma_tlf'] = flat.reshape((nbf * nbf, nstate, nstate))

def grad_now():
    oqp.tdhf_mrsf_gradient(mol); return mol.get_grad().reshape(-1).copy()

def nuclear_from_rhs(i, j, R):
    push_gamma(inject_wrk1(R)); mol.data.set_tdhf_target(i)
    oqp.set_mrsf_nac_cphf(mol, i, j); oqp.tdhf_mrsf_z_vector(mol)
    if not mol.mol_energy.Z_Vector_converged:
        oqp.set_mrsf_nac_cphf(mol, 0, 0); return None
    gZ = grad_now()
    mol.data['OQP::td_p'] = np.zeros_like(np.array(mol.data['OQP::td_p'], copy=True))
    mol.data['OQP::WAO'] = np.zeros_like(np.array(mol.data['OQP::WAO'], copy=True))
    gS = grad_now(); oqp.set_mrsf_nac_cphf(mol, 0, 0)
    return (gZ - gS).reshape(-1)

qz = np.load('/tmp/nactest/qz_full.npz')   # nuclear oracle (9-comp)
qzn = np.load('/tmp/nactest/qz_nogate.npz')
def Dnuc(I, J):
    cc = qzn[f'Z{I}{J}_coded'] - 0.5*(qzn[f'X{I}_coded']+qzn[f'X{J}_coded'])
    cn = qzn[f'Z{I}{J}_numeric'] - 0.5*(qzn[f'X{I}_numeric']+qzn[f'X{J}_numeric'])
    return cc - cn
def cos9(a, b):
    return np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30), np.linalg.norm(a)/(np.linalg.norm(b)+1e-30)

print("=== seam(antisym sp cross) direction vs nuclear deficiency D=coded-numeric ===")
for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    i, j = I+1, J+1
    a = anti[f'anti_{i}{j}']
    g = nuclear_from_rhs(i, j, a)
    D = Dnuc(i, j)
    if g is None:
        print(f"  ({i},{j}): Z not converged"); continue
    c, rr = cos9(g, D)
    print(f"  ({i},{j}): cos(seam_anti, D)={c:+.5f} ratio={rr:.4f} |seam_anti|={np.linalg.norm(g):.5e} |D|={np.linalg.norm(D):.5e}")
    print(f"           seam_anti = {np.array2string(g, precision=4, suppress_small=True)}")
    print(f"           D         = {np.array2string(D, precision=4, suppress_small=True)}")
