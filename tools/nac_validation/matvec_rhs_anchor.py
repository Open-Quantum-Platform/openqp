"""Phase 11: with the matvec made LINEAR by tight screening, finite-difference its
orbital-rotation gradient L^matvec = d(X^T A X)/dtheta and identify what it equals on
the DIAGONAL. The matvec is TDA (no B-matrix), so the expectation is
  L^matvec(I,I)  ==  -(sfrorhs without AB1)
i.e. the 2e back-transform + 2FT relaxation, with the (A+B)/ab1 term ABSENT.
This anchors L^matvec before the interstate/qdots projection test.
"""
import os
import oqp
from oqp.pyoqp import Runner
import numpy as np

r = Runner(input_file='/tmp/nactest/H2O_en.inp', log='/tmp/nactest/anchor.log')
r.run(); mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20            # LINEAR matvec
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca-2
C0a = np.array(mol.data['OQP::VEC_MO_A'], copy=True); nbf = C0a.shape[0]
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
nvirb = nbf-nocb; nij = noca*nvirb; nstate = 3
Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); shp = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))

prs = ([(i, j, 'b') for i in range(nocb, noca) for j in range(0, nocb)]
       + [(k, j, 'ab') for k in range(noca, nbf) for j in range(0, nocb)]
       + [(k, i, 'a') for k in range(noca, nbf) for i in range(nocb, noca)])
nconf = len(prs)
nds, ndv = nocb*(noca-nocb), (nbf-noca)*nocb
blocks = [('doc-socc', 0, nds), ('doc-virt', nds, nds+ndv), ('soc-virt', nds+ndv, nconf)]
TH = 1e-4


def setb(c):
    raw = Xraw.copy().reshape(-1); raw[0:nij] = c
    mol.data['OQP::td_bvec_mo'] = raw.reshape(shp)


def Ax(c):
    setb(c); oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()


def givens(C, a, b, th):
    Cr = C.copy(); cc, ss = np.cos(th), np.sin(th)
    Cr[:, a] = cc*C[:, a] + ss*C[:, b]; Cr[:, b] = -ss*C[:, a] + cc*C[:, b]
    return Cr


def Lmatvec(col):
    """frozen-Fock orbital-rotation gradient d(col^T A col)/dtheta over prs."""
    L = np.zeros(nconf)
    for n, (a, b, sp) in enumerate(prs):
        def f(th):
            Ca = givens(C0a, a, b, th) if sp in ('a', 'ab') else C0a
            Cb = givens(C0b, a, b, th) if sp in ('b', 'ab') else C0b
            mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(Ca)
            mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(Cb)
            return float(np.dot(col, Ax(col)))
        L[n] = (f(TH) - f(-TH))/(2*TH)
    mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(C0a)
    mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(C0b)
    return L


def dump(target, gates):
    mol.data.set_tdhf_target(target)
    allg = ('NAC_ZERO_2FT', 'NAC_ZERO_AB1', 'NAC_ZERO_HX', 'NAC_ZERO_T', 'NAC_ZERO_SP')
    for k in allg:
        os.environ.pop(k, None)
    for k in gates:
        os.environ[k] = '1'
    os.environ['NAC_DUMP_RHS'] = '1'
    oqp.tdhf_mrsf_z_vector(mol)
    os.environ.pop('NAC_DUMP_RHS')
    for k in allg:
        os.environ.pop(k, None)
    return np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).ravel()[:nconf]


def cmp(a, b):
    return (np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30),
            np.linalg.norm(a)/(np.linalg.norm(b)+1e-30))


print("=== L^matvec(I,I) vs -sfrorhs variants (tight screening) ===")
variants = {'full': [], 'noAB1': ['NAC_ZERO_AB1'],
            'no2FT': ['NAC_ZERO_2FT'], 'A-only(2e)': ['NAC_ZERO_2FT', 'NAC_ZERO_AB1']}
for st in range(1, nstate+1):
    L = Lmatvec(X[st-1])
    print(f"\nstate {st}: |L^matvec|={np.linalg.norm(L):.4e}")
    for name, g in variants.items():
        R = dump(st, g)
        c, rr = cmp(L, -R)
        print(f"   vs -sfrorhs[{name:10}]: cos={c:+.6f} ratio={rr:.4f}  (|R|={np.linalg.norm(R):.4e})")
