"""Definitive in-process dataset: class-resolved bilinears + damp_PT."""
import math
import oqp, oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint
import numpy as np

inp = '/tmp/nactest/H2O_tight_dx0.001.inp'
DELTA = 1e-3
r = Runner(input_file=inp, log='/tmp/nactest/cf.log')
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
X0 = X.T.copy()
RS = 1/math.sqrt(2.0)

ijlr1 = (noca-1-nocb-1)*noca + noca - 1
ijg = (noca-1-nocb-1)*noca + noca
ijd = (noca-nocb-1)*noca + noca - 1
Pg = np.zeros(nij); Pg[ijg-1] = 1; Pg[ijd-1] = 1
Pl = np.zeros(nij); Pl[ijlr1-1] = 1
Px = 1.0 - Pg - Pl

def Gfull(col):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(X0_raw.shape)
    mol.data.set_tdhf_target(1)
    oqp.tdhf_mrsf_z_vector(mol)
    oqp.tdhf_mrsf_gradient(mol)
    g = mol.get_grad().reshape(-1).copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    return g

def B(u, v):
    return 0.25*(Gfull(u+v) - Gfull(u-v))

classes = {'g': Pg, 'l': Pl, 'x': Px}
Bsub = {}
for (i, j) in ((1,2),(1,3),(2,3)):
    XI, XJ = X[i-1], X[j-1]
    for a in 'glx':
        for b in 'glx':
            key = (i, j, a, b)
            u, v = XI*classes[a], XJ*classes[b]
            if np.linalg.norm(u) < 1e-10 or np.linalg.norm(v) < 1e-10:
                Bsub[key] = np.zeros(3*natom)
            else:
                Bsub[key] = B(u, v)

# damp_PT loop
mol.save_data()
cfg = mol.config; cfg['guess']['type']='json'
cfg['guess']['file']='/tmp/nactest/cf.json'; cfg['guess']['continue_geom']=False
OV='OQP::overlap_mo_non_orthogonal'

def unfold(v):
    c = np.zeros((noca, nvirb))
    for i in range(1, noca+1):
        for jj in range(nocb+1, nbf+1):
            ij = (jj-nocb-1)*noca + i
            if ij == ijlr1: c[i-1, jj-nocb-1] = v[ijlr1-1]*RS
            elif ij == ijlr1 + noca + 1: c[i-1, jj-nocb-1] = -v[ijlr1-1]*RS
            else: c[i-1, jj-nocb-1] = v[ij-1]
    return c

def refold(c):
    out = c.copy(); out[nocb,0] = math.sqrt(2.0)*c[nocb,0]; out[nocb+1,1] = 0.0
    return out.T.reshape(-1)

QDOTS = []

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
    QDOTS.append(Q.copy())
    return Xal.T

damp_PT = np.zeros((9, nstate, nstate))
for k in range(9):
    Xp = disp(xyz0 + DELTA*np.eye(9)[k])
    Xm = disp(xyz0 - DELTA*np.eye(9)[k])
    for Xx in (Xp, Xm):
        for st in range(nstate):
            if np.dot(X0[:,st], Xx[:,st]) < 0: Xx[:,st] *= -1
    damp_PT[k] = X0.T @ ((Xp-Xm)/(2*DELTA))

qdots = np.array([(QDOTS[2*k] - QDOTS[2*k+1])/(2*DELTA) for k in range(9)])
np.savez('/tmp/nactest/classfit.npz', damp_PT=damp_PT, qdots=qdots,
         X0=X0, noca=noca, nocb=nocb, nbf=nbf,
         energies=np.array(mol.energies),
         **{f'B_{i}{j}_{a}{b}': Bsub[(i,j,a,b)]
            for (i,j,a,b) in Bsub})
print('saved classfit.npz')
