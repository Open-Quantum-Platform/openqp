"""Phase 10: exact generator-convention scan.
transport == damp_PT - damp_raw is a pure data identity; scan L-candidates
against it. Separately test damp_raw vs dci (polarization, raw gauge)."""
import math
import oqp, oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint
import numpy as np

inp = '/tmp/nactest/H2O_tight_dx0.001.inp'
DELTA = 1e-3
r = Runner(input_file=inp, log='/tmp/nactest/gt2.log')
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

def refold(c):
    f = c.copy(); v = np.zeros(nij)
    f2 = f.T.reshape(-1)  # placeholder
    out = c.copy()
    out[nocb, 0] = math.sqrt(2.0)*c[nocb, 0]
    out[nocb+1, 1] = 0.0
    return out.T.reshape(-1)

# dci via polarization
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

mol.save_data()
cfg = mol.config; cfg['guess']['type']='json'
cfg['guess']['file']='/tmp/nactest/gt2.json'; cfg['guess']['continue_geom']=False
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
    sgn = np.sign(np.diag(M)); sgn[sgn==0]=1
    Xd = np.array(mol.data['OQP::td_bvec_mo'], copy=True).reshape(-1).reshape((nstate,nij))
    XPT = np.zeros_like(Xd); XRW = np.zeros_like(Xd)
    for st in range(nstate):
        c = unfold(Xd[st])
        cp = Q[:noca,:noca].T @ c @ Q[nocb:,nocb:]
        XPT[st] = refold(cp)
        cr = np.diag(sgn[:noca]) @ c @ np.diag(sgn[nocb:])   # sign-only gauge
        XRW[st] = refold(cr)
    return XPT.T, XRW.T, Q

damp_PT = np.zeros((9, nstate, nstate)); damp_RW = np.zeros_like(damp_PT)
qdots = np.zeros((9, nbf, nbf))
for k in range(9):
    Pp, Rp, Qp = disp(xyz0 + DELTA*np.eye(9)[k])
    Pm, Rm, Qm = disp(xyz0 - DELTA*np.eye(9)[k])
    for Xx in (Pp, Pm, Rp, Rm):
        for st in range(nstate):
            if np.dot(X0[:,st], Xx[:,st]) < 0: Xx[:,st] *= -1
    damp_PT[k] = X0.T @ ((Pp-Pm)/(2*DELTA))
    damp_RW[k] = X0.T @ ((Rp-Rm)/(2*DELTA))
    qdots[k] = (Qp - Qm)/(2*DELTA)

np.savez('/tmp/nactest/gt2.npz', damp_PT=damp_PT, damp_RW=damp_RW,
         qdots=qdots, X0=X0,
         dci12=dci[(1,2)], dci13=dci[(1,3)], dci23=dci[(2,3)],
         noca=noca, nocb=nocb, nbf=nbf)
print('saved gt2.npz')

cU = [unfold(X0[:, s]) for s in range(nstate)]
print('\n===== A: damp_raw vs dci (polarization validity, raw gauge) =====')
for (i, j) in ((1,2),(1,3),(2,3)):
    da = damp_RW[:, i-1, j-1]; dc = dci[(i,j)]
    c = np.dot(dc, da)/(np.linalg.norm(dc)*np.linalg.norm(da)+1e-30)
    print(f'({i},{j}): |damp_raw|={np.linalg.norm(da):.5f} |dci|={np.linalg.norm(dc):.5f} '
          f'cos={c:+.6f} ratio={np.linalg.norm(dc)/np.linalg.norm(da):.4f}')

print('\n===== B: generator scan vs exact target (damp_PT - damp_raw) =====')
def transport(qsign, transpose_occ):
    tr = np.zeros((9, nstate, nstate))
    for k in range(9):
        q = qsign*qdots[k]
        qo = q[:noca,:noca]; qv = q[nocb:,nocb:]
        if transpose_occ: qo = qo.T; qv = qv.T
        for I in range(nstate):
            for J in range(nstate):
                dcJ = qo.T @ cU[J] + cU[J] @ qv
                tr[k, I, J] = np.sum(cU[I]*dcJ)
    return tr

target = damp_PT - damp_RW
for qsign in (1, -1):
    for tro in (False, True):
        tr = transport(qsign, tro)
        errs = []
        for (i, j) in ((1,2),(1,3),(2,3)):
            t = target[:, i-1, j-1]; v = tr[:, i-1, j-1]
            errs.append(np.linalg.norm(v-t))
        print(f'qsign={qsign:+d} transposed={tro}: resid per pair = '
              + '  '.join(f'{e:.5f}' for e in errs)
              + f'   |target(1,3)|={np.linalg.norm(target[:,0,2]):.5f}')
