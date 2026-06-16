"""Phase 11 Task 4 (Step-3 core): interstate bilinear 2e RHS from G_MO.

The amplitude term's z-vector RHS is the orbital-rotation gradient of X_I^T A X_J.
Its 2-electron back-transform part is the symmetric bilinear generalization of the
validated diagonal hxa=2 G_MO[X] X_f^T:

    G_MO[.] is LINEAR in the amplitude (agdlr = int2(mrsfcbc(x)) is linear in x),
    so by the polarization identity the bilinear 2e kernel is
        hxa_IJ = G_I wrk3_J^T + G_J wrk3_I^T
        hxb_IJ = G_I^T wrk3_J + G_J^T wrk3_I      (rows 0:noca of G, wrk3)
    with G_K = G_MO[X_K] (OQP::nac_gmo) and wrk3_K = iatogen(mrsfxvec(X_K)).
Then the sfrorhs block projection (2e part only; 2FT/AB1 relaxation deferred).

Checks: (1) G_MO linearity G[aX_I+bX_J]=aG_I+bG_J; (2) the bilinear reduces to the
validated diagonal RHS at I=J to ~1e-12. This is the genuinely-different
ground-config-channel object the matvec route provides (vs the deficient gradient
chain); the cross RHS for (1,2)/(1,3)/(2,3) is reported for the next step.
"""
import oqp
from oqp.pyoqp import Runner
import numpy as np

INP = '/tmp/nactest/H2O_en.inp'
SQ = 1.0 / np.sqrt(2.0)

r = Runner(input_file=INP, log='/tmp/nactest/gmo_inter.log')
r.run(); mol = r.mol
print("RUN OK", flush=True)

noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
nbf = np.array(mol.data['OQP::VEC_MO_A'], copy=True).shape[0]
nvirb = nbf - nocb; nij = noca * nvirb; nstate = 3
Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); Xshape = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))
Om = np.array(mol.data['OQP::td_energies'], copy=True).ravel()

prs = ([(i, j) for i in range(nocb, noca) for j in range(0, nocb)]
       + [(k, j) for k in range(noca, nbf) for j in range(0, nocb)]
       + [(k, i) for k in range(noca, nbf) for i in range(nocb, noca)])
nconf = len(prs)
nds, ndv = nocb*(noca-nocb), (nbf-noca)*nocb
blocks = [('doc-socc', 0, nds), ('doc-virt', nds, nds+ndv), ('soc-virt', nds+ndv, nconf)]
ijlr1 = (noca-1-nocb-1)*noca + (noca-1) - 1
ijlr2 = (noca  -nocb-1)*noca + (noca  ) - 1


def set_bvec(col):
    raw = Xraw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def get_gmo(col):
    set_bvec(col); oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_gmo'], copy=True).ravel().reshape((nbf, nbf), order='F')


def mrsfxvec(col):
    f = col.copy(); f[ijlr1] = col[ijlr1]*SQ; f[ijlr2] = -col[ijlr1]*SQ
    return f


def iatogen(colf):
    av = np.zeros((nbf, nbf)); av[0:noca, nocb:nbf] = colf.reshape((nvirb, noca)).T
    return av


def project(hxa, hxb):
    rhs = np.zeros(nconf); n = 0
    for i in range(nocb, noca):
        for j in range(0, nocb):
            rhs[n] = hxa[i, j] - hxa[j, i] - hxb[j, i]; n += 1
    for k in range(noca, nbf):
        for j in range(0, nocb):
            rhs[n] = hxa[k, j] - hxb[j, k]; n += 1
    for k in range(noca, nbf):
        for i in range(nocb, noca):
            rhs[n] = hxa[k, i] + hxb[k, i] - hxb[i, k]; n += 1
    return -rhs


def rhs_2e_diag(col):
    """validated diagonal 2e RHS (== sfrorhs SP+2FT+AB1 off)."""
    G = get_gmo(col); w = iatogen(mrsfxvec(col))
    hxa = 2.0 * G @ w[0:noca, :].T
    hxb = 2.0 * G[0:noca, :].T @ w[0:noca, :]
    return project(hxa, hxb)


def rhs_2e_inter(GI, GJ, wI, wJ):
    """symmetric bilinear 2e RHS for pair (I,J)."""
    hxa = GI @ wJ[0:noca, :].T + GJ @ wI[0:noca, :].T
    hxb = GI[0:noca, :].T @ wJ[0:noca, :] + GJ[0:noca, :].T @ wI[0:noca, :]
    return project(hxa, hxb)


# ---- (1) G_MO linearity ----------------------------------------------------
print("\n=== (1) G_MO linearity: G[aX1+bX3] vs aG1+bG3 ===")
G1, G3 = get_gmo(X[0]), get_gmo(X[2])
a, b = 0.7, -1.3
Gmix = get_gmo(a*X[0] + b*X[2])
print(f"  max|G[aX1+bX3] - (aG1+bG3)| = {np.max(np.abs(Gmix - (a*G1+b*G3))):.2e}  (linear => ~1e-12)")

# ---- (2) I=J reduction -----------------------------------------------------
print("\n=== (2) interstate bilinear reduces to validated diagonal at I=J ===")
for st in range(1, nstate+1):
    G = get_gmo(X[st-1]); w = iatogen(mrsfxvec(X[st-1]))
    Rii = rhs_2e_inter(G, G, w, w)
    Rdi = rhs_2e_diag(X[st-1])
    print(f"  state {st}: max|R_inter(I=I) - R_diag| = {np.max(np.abs(Rii - Rdi)):.2e}")

# ---- (3) cross RHS for the 3 pairs -----------------------------------------
print("\n=== (3) interstate 2e cross RHS  (the bilinear amplitude-term kernel) ===")
gmos = [get_gmo(X[k]) for k in range(nstate)]
wf = [iatogen(mrsfxvec(X[k])) for k in range(nstate)]
for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    R = rhs_2e_inter(gmos[I], gmos[J], wf[I], wf[J])
    line = f"  pair ({I+1},{J+1}): |R_2e|={np.linalg.norm(R):.4e}"
    for name, lo, hi in blocks:
        line += f"  {name}={np.linalg.norm(R[lo:hi]):.3e}"
    print(line)
print("\nNOTE: this is the 2e back-transform part only (common to matvec & gradient")
print("chain). The relaxation term (matvec mrsfesum-F.X vs gradient 2FT) + AB1 are")
print("the remaining Step-3 pieces; the oracle is test_QZ numeric_cross.")
