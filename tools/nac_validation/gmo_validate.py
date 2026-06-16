"""Phase 11 Task 3: validate the exported full-MO 2e kernel G_MO.

We reproduce the production z-vector RHS hxa/hxb (the 2-electron part, A-only)
ANALYTICALLY from the exported OQP::nac_gmo, with the sqrt2 SOMO fold applied on
the INPUT (mrsfxvec) -- NOT via finite difference, NOT via mrsfesum's output fold.

Production construction (tdhf_mrsf_z_vector.F90:1707-1727):
    G_MO = mo_a^T agdlr mo_a                       (== OQP::nac_gmo, exported)
    X_folded = iatogen(mrsfxvec(X))                (input sqrt2 SOMO fold)
    hxa(p,i) = 2 sum_k G_MO(p,k) X_folded(i,k)     i in occ_alpha
    hxb(p,q) = 2 sum_{k in occ_a} G_MO(k,p) X_folded(k,q)
then the sfrorhs block projection (tdhf_sf_lib.F90:405-432), 2FT+AB1 OFF.

Target = the dumped production RHS with NAC_ZERO_SP/2FT/AB1 all set (pure G_MO
2e projection). A match to ~1e-10 element-wise validates the G_MO export and the
input-fold reproduction.
"""
import oqp                       # MUST precede numpy (ILP64 LAPACK clash)
from oqp.pyoqp import Runner
import numpy as np
import os

INP = '/tmp/nactest/H2O_en.inp'
SQ = 1.0 / np.sqrt(2.0)

r = Runner(input_file=INP, log='/tmp/nactest/gmo_val.log')
r.run(); mol = r.mol
print("RUN OK", flush=True)

noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
C0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = C0.shape[0]
nvirb = nbf - nocb
nij = noca * nvirb
nstate = 3

Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); Xshape = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))
Om = np.array(mol.data['OQP::td_energies'], copy=True).ravel()
print(f"nbf={nbf} noca={noca} nocb={nocb} nvirb={nvirb} nij={nij}")

# rotation-space block layout (matches sfrogen / sfrorhs / harness prs order)
prs = ([(i, j) for i in range(nocb, noca) for j in range(0, nocb)]        # doc-socc
       + [(k, j) for k in range(noca, nbf) for j in range(0, nocb)]       # doc-virt
       + [(k, i) for k in range(noca, nbf) for i in range(nocb, noca)])   # soc-virt
nconf = len(prs)
nds, ndv = nocb*(noca-nocb), (nbf-noca)*nocb
blocks = [('doc-socc', 0, nds), ('doc-virt', nds, nds+ndv), ('soc-virt', nds+ndv, nconf)]

# slot indices (1-based ijlr1=5, ijlr2=12 for noca6/nocb4) -> 0-based
ijlr1 = (noca-1-nocb-1)*noca + (noca-1) - 1
ijlr2 = (noca  -nocb-1)*noca + (noca  ) - 1
print(f"ijlr1(0based)={ijlr1} ijlr2(0based)={ijlr2}")


def set_bvec(col):
    raw = Xraw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def get_gmo(col):
    """Run the matvec (current orbitals) and return G_MO as a Fortran-order matrix."""
    set_bvec(col)
    oqp.mrsf_matvec_apply(mol)
    g = np.array(mol.data['OQP::nac_gmo'], copy=True).ravel()
    return g.reshape((nbf, nbf), order='F')        # G[p,q] = gmo(p+1,q+1)


def mrsfxvec(col):
    """Input-side sqrt2 SOMO fold (mult=1)."""
    f = col.copy()
    f[ijlr1] = col[ijlr1] * SQ
    f[ijlr2] = -col[ijlr1] * SQ
    return f


def iatogen(colf):
    """Compressed (nij) -> nbf x nbf amplitude matrix, (occ_a, virt_b) block,
    column-major with i (occ_alpha) fastest."""
    av = np.zeros((nbf, nbf))
    av[0:noca, nocb:nbf] = colf.reshape((nvirb, noca)).T
    return av


def analytic_rhs(col):
    """Reproduce sfrorhs A-only (2FT+AB1+SP off) from the exported G_MO."""
    G = get_gmo(col)
    wrk3 = iatogen(mrsfxvec(col))
    hxa = 2.0 * G @ wrk3[0:noca, :].T              # (nbf, noca)
    hxb = 2.0 * G[0:noca, :].T @ wrk3[0:noca, :]   # (nbf, nbf)
    rhs = np.zeros(nconf)
    n = 0
    for i in range(nocb, noca):          # doc-socc: socc i, doc j
        for j in range(0, nocb):
            rhs[n] = hxa[i, j] - hxa[j, i] - hxb[j, i]; n += 1
    for k in range(noca, nbf):           # doc-virt: virt k, doc j
        for j in range(0, nocb):
            rhs[n] = hxa[k, j] - hxb[j, k]; n += 1
    for k in range(noca, nbf):           # soc-virt: virt k, socc i
        for i in range(nocb, noca):
            rhs[n] = hxa[k, i] + hxb[k, i] - hxb[i, k]; n += 1
    return -rhs                          # sfrorhs final *-1


def dump_sfrorhs(target, gates):
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
            np.linalg.norm(a)/(np.linalg.norm(b)+1e-30),
            np.max(np.abs(a-b)))


print("\n=== DIAGONAL G_MO validation: analytic hxa/hxb  vs  sfrorhs (SP+2FT+AB1 off) ===")
GATES = ['NAC_ZERO_SP', 'NAC_ZERO_2FT', 'NAC_ZERO_AB1']
for st in range(1, nstate+1):
    Ran = analytic_rhs(X[st-1])
    Rsf = dump_sfrorhs(st, GATES)
    c, rr, mx = cmp(Ran, Rsf)
    print(f"\nstate {st}: |Ran|={np.linalg.norm(Ran):.4e} |Rsf|={np.linalg.norm(Rsf):.4e}")
    print(f"  OVERALL cos={c:+.8f}  ratio={rr:.6f}  max|diff|={mx:.2e}")
    for name, lo, hi in blocks:
        bc, br, bmx = cmp(Ran[lo:hi], Rsf[lo:hi])
        print(f"    {name:9} cos={bc:+.8f} ratio={br:.6f} max|diff|={bmx:.2e}")
