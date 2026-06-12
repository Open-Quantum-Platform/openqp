"""Verify eq:gaugetransport: damp_PT =?= X(dA)X/dOmega + transport."""
import math
import oqp, oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint
import numpy as np

inp = '/tmp/nactest/H2O_tight_dx0.001.inp'
DELTA = 1e-3
r = Runner(input_file=inp, log='/tmp/nactest/gt.log')
r.run()
mol = r.mol
nstate, natom = 3, mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True); nbf = C0raw.shape[0]
nvirb = nbf - nocb; nij = noca * nvirb
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
xyz0 = np.array(mol.get_system(), copy=True)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
X0 = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
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

# ---- dci via polarization (canonical-gauge amplitude term) ----
def grad_for(t, col):
    raw = X0_raw.copy().reshape(-1)
    raw[(t-1)*nij:t*nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(X0_raw.shape)
    mol.data.set_tdhf_target(t)
    oqp.tdhf_mrsf_z_vector(mol)
    oqp.tdhf_mrsf_gradient(mol)
    g = mol.get_grad().reshape(-1).copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    return g

dci = {}
for (i, j) in ((1,2),(1,3),(2,3)):
    gp = grad_for(i, (X0[:,i-1]+X0[:,j-1])/math.sqrt(2))
    gm = grad_for(i, (X0[:,i-1]-X0[:,j-1])/math.sqrt(2))
    dci[(i,j)] = 0.5*(gp-gm)/(mol.energies[j]-mol.energies[i])

# ---- displaced loop: damp_PT + rotation rates ----
mol.save_data()
cfg = mol.config; cfg['guess']['type']='json'
cfg['guess']['file']='/tmp/nactest/gt.json'; cfg['guess']['continue_geom']=False
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
    Xal = np.zeros_like(Xd)
    for st in range(nstate):
        c = unfold(Xd[st]); cp = Q[:noca,:noca].T @ c @ Q[nocb:,nocb:]
        f = cp.copy(); f[nocb,0] = math.sqrt(2)*cp[nocb,0]; f[nocb+1,1] = 0.0
        Xal[st] = f.T.reshape(-1)
    return Xal.T, Q

damp = np.zeros((9, nstate, nstate)); trans = np.zeros((9, nstate, nstate))
cU = [unfold(X0[:, s]) for s in range(nstate)]
for k in range(9):
    Xp, Qp = disp(xyz0 + DELTA*np.eye(9)[k])
    Xm, Qm = disp(xyz0 - DELTA*np.eye(9)[k])
    for Xx in (Xp, Xm):
        for st in range(nstate):
            if np.dot(X0[:,st], Xx[:,st]) < 0: Xx[:,st] *= -1
    damp[k] = X0.T @ ((Xp-Xm)/(2*DELTA))
    # rotation rate (antisym generator; Q ~ exp(q) so q' ~ (Qp-Qm)/2d)
    qdot = (Qp - Qm)/(2*DELTA)
    qo = qdot[:noca,:noca]; qv = qdot[nocb:,nocb:]
    # transport: d/d(gauge) of c_J = qo^T c + c qv ; term = sum c_I * that
    for I in range(nstate):
        for J in range(nstate):
            dcJ = qo.T @ cU[J] + cU[J] @ qv
            trans[k, I, J] = np.sum(cU[I]*dcJ)

print('\n===== gauge-transport verification =====')
for (i, j) in ((1,2),(1,3),(2,3)):
    da = damp[:, i-1, j-1]; tr = trans[:, i-1, j-1]; dc = dci[(i,j)]
    pred = dc + tr
    c = np.dot(pred, da)/(np.linalg.norm(pred)*np.linalg.norm(da)+1e-30)
    print(f'({i},{j}): |damp|={np.linalg.norm(da):.5f} |dci|={np.linalg.norm(dc):.5f} '
          f'|transport|={np.linalg.norm(tr):.5f}')
    print(f'        dci+transport vs damp: cos={c:+.6f} '
          f'ratio={np.linalg.norm(pred)/np.linalg.norm(da):.4f} '
          f'resid={np.linalg.norm(pred-da):.5f}')
