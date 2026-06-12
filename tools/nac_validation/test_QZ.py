"""Is the coded quadratic form G(Z)-G_SCF equal to d(Z^T A Z)/dR for mixed Z?"""
import math
import oqp, oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint
import numpy as np

inp = '/tmp/nactest/H2O_tight_dx0.001.inp'
DELTA = 1e-3
r = Runner(input_file=inp, log='/tmp/nactest/qz.log')
r.run()
mol = r.mol
nstate, natom = 3, mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True); nbf = C0raw.shape[0]
nvirb = nbf - nocb; nij = noca*nvirb
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
xyz0 = np.array(mol.get_system(), copy=True)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
X = X0_raw.reshape(-1).reshape((nstate, nij))
Om0 = [mol.energies[k+1] - mol.energies[0] + (mol.energies[1]-mol.energies[1]) for k in range(nstate)]
# excitation energies: energies[1..3] are total state energies; Omega_K relative
Etot0 = list(mol.energies)
RS = 1/math.sqrt(2.0)

def unfold(v):
    ijlr1 = (noca-1-nocb-1)*noca + noca-1; ijlr2 = (noca-nocb-1)*noca + noca
    c = np.zeros((noca, nvirb))
    for i in range(1, noca+1):
        for jj in range(nocb+1, nbf+1):
            ij = (jj-nocb-1)*noca + i
            if ij == ijlr1: c[i-1, jj-nocb-1] = v[ijlr1-1]*RS
            elif ij == ijlr2: c[i-1, jj-nocb-1] = -v[ijlr1-1]*RS
            else: c[i-1, jj-nocb-1] = v[ij-1]
    return c

def refold(c):
    out = c.copy(); out[nocb,0] = math.sqrt(2.0)*c[nocb,0]; out[nocb+1,1] = 0.0
    return out.T.reshape(-1)

def Gfull(col):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(X0_raw.shape)
    mol.data.set_tdhf_target(1)
    oqp.tdhf_mrsf_z_vector(mol)
    oqp.tdhf_mrsf_gradient(mol)
    gZ = mol.get_grad().reshape(-1).copy()
    for tag in ('OQP::td_p', 'OQP::WAO', 'OQP::td_abxc', 'OQP::td_mrsf_density'):
        mol.data[tag] = np.zeros_like(np.array(mol.data[tag], copy=True))
    oqp.tdhf_mrsf_gradient(mol)
    gS = mol.get_grad().reshape(-1).copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    return gZ - gS                       # the coded dOmega(Z)

Zs = {'Z13': (X[0]+X[2])/math.sqrt(2), 'Z23': (X[1]+X[2])/math.sqrt(2),
      'Z12': (X[0]+X[1])/math.sqrt(2),
      'X1': X[0].copy(), 'X2': X[1].copy(), 'X3': X[2].copy()}
GZ = {n: Gfull(z) for n, z in Zs.items()}

mol.save_data()
cfg = mol.config; cfg['guess']['type']='json'
cfg['guess']['file']='/tmp/nactest/qz.json'; cfg['guess']['continue_geom']=False
OV='OQP::overlap_mo_non_orthogonal'

def disp(coord):
    mol.update_system(coord); oqp.library.ints_1e(mol); oqp.library.guess(mol)
    SinglePoint(mol).energy()
    mol.data['OQP::xyz_old'] = xyz0.reshape((3, natom))
    mol.data['OQP::VEC_MO_A_old'] = C0raw; mol.data['OQP::E_MO_A_old'] = e0
    mol.data['OQP::VEC_MO_B_old'] = C0b;  mol.data['OQP::E_MO_B_old'] = e0b
    oqp.get_structures_ao_overlap(mol)
    M = np.array(mol.data[OV], copy=True).reshape(-1).reshape((nbf,nbf)).T
    Q = np.zeros((nbf,nbf))
    for lo,hi in ((0,nocb),(nocb,noca),(noca,nbf)):
        W,s,Vt = np.linalg.svd(M[lo:hi,lo:hi]); Q[lo:hi,lo:hi] = Vt.T@W.T
    Xd = np.array(mol.data['OQP::td_bvec_mo'], copy=True).reshape(-1).reshape((nstate,nij))
    Xal = np.zeros((nstate, nij))
    for st in range(nstate):
        c = unfold(Xd[st]); cp = Q[:noca,:noca].T @ c @ Q[nocb:,nocb:]
        Xal[st] = refold(cp)
    Om = [mol.energies[k+1] - mol.energies[1] + (Etot0[1]-Etot0[1]) for k in range(nstate)]
    # use total state energies minus displaced SCF reference? simplest: total energies
    Etot = list(mol.energies)
    return Xal, Etot

results = {n: {'coded': GZ[n], 'numeric': np.zeros(9)} for n in Zs}
print('\n===== coded Q(Z) vs exact d(Z^T A Z) per coordinate =====')
for k in range(9):
    Xp, Ep = disp(xyz0 + DELTA*np.eye(9)[k])
    Xm, Em = disp(xyz0 - DELTA*np.eye(9)[k])
    for n, z in Zs.items():
        # Omega(Z; d) = sum_K Omega_K(d) <z|X_K(d)>^2 with Omega_K = E_K - E_SCFref...
        # use total TD state energies relative to nothing: Z^T A Z + const*|Z|^2?
        # A's eigenvalues are the EXCITATION energies Omega_K = energies[K] - energies[ref?]
        # OpenQP mol.energies: [E_SCF? E_1 E_2 E_3]? energies[j] used as total before.
        def om(E, Xal):
            tot = 0.0
            for K in range(nstate):
                ov = np.dot(z, Xal[K])
                tot += (E[K+1] - E[0]) * ov*ov   # excitation energies
            return tot
        num = (om(Ep, Xp) - om(Em, Xm)) / (2*DELTA)
        results[n]['numeric'][k] = num
np.savez('/tmp/nactest/qz_full.npz',
         **{f'{n}_{w}': results[n][w] for n in results for w in ('coded','numeric')})
print('saved qz_full.npz')
