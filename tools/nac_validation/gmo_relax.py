"""Phase 11: validate the 2FT (gradient-chain) interstate relaxation in Python,
using the matvec-exported MO Fock fa,fb (OQP::nac_fa/fb) -- the bit-identical Fock
mrsfesum/sfrorhs use. This is the DEFICIENT baseline the matvec esum-relaxation
must differ from off-diagonal. Prep harness: once the esum-relaxation formula is
in hand, drop it in `relax_esum()` and compare L_esum - L_2ft per channel.

2FT (sfrorhs, tdhf_sf_lib.F90:382-403): xhxa += 2 Fa Tij (Tij occ-occ, noca),
xhxb += 2 Fb Tab (Tab placed in the vir-vir [nocb:,nocb:] block), then the same
block projection and *(-1). Tij/Tab = mrsf_interstate_tden (input-folded amps):
Tij=-1/2(Xi Xj^T+Xj Xi^T), Tab=+1/2(Xi^T Xj+Xj^T Xi).

CHECK: L_2ft(I=I) must equal the Fortran 2FT contribution
  R_2ft_fortran = dump[SP,AB1 off, 2FT ON] - dump[SP,2FT,AB1 off]   (per state).
"""
import oqp
from oqp.pyoqp import Runner
import numpy as np
import os

INP = '/tmp/nactest/H2O_en.inp'
SQ = 1.0/np.sqrt(2.0)

r = Runner(input_file=INP, log='/tmp/nactest/gmo_relax.log')
r.run(); mol = r.mol
print("RUN OK", flush=True)

noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca-2
nbf = np.array(mol.data['OQP::VEC_MO_A'], copy=True).shape[0]
nvirb = nbf-nocb; nij = noca*nvirb; nstate = 3
Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); Xshape = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))

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


def get_fock(col):
    set_bvec(col); oqp.mrsf_matvec_apply(mol)
    fa = np.array(mol.data['OQP::nac_fa'], copy=True).ravel().reshape((nbf, nbf), order='F')
    fb = np.array(mol.data['OQP::nac_fb'], copy=True).ravel().reshape((nbf, nbf), order='F')
    return fa, fb


def mrsfxvec_mat(col):
    """input fold -> (noca x nvirb) amplitude matrix."""
    f = col.copy(); f[ijlr1] = col[ijlr1]*SQ; f[ijlr2] = -col[ijlr1]*SQ
    return f.reshape((nvirb, noca)).T          # (noca, nvirb)


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


def relax_2ft(I, J):
    """gradient-chain 2FT interstate relaxation, projected (no 2e, no ab1)."""
    fa, fb = get_fock(X[I])                      # fa,fb are state-independent (frozen ROHF Fock)
    xi, xj = mrsfxvec_mat(X[I]), mrsfxvec_mat(X[J])
    Tij = -0.5*(xi @ xj.T + xj @ xi.T)           # (noca, noca)
    Tab = 0.5*(xi.T @ xj + xj.T @ xi)            # (nvirb, nvirb)
    xhxa = 2.0 * fa[:, 0:noca] @ Tij             # (nbf, noca)
    wrk = np.zeros((nbf, nbf)); wrk[nocb:, nocb:] = Tab
    xhxb = 2.0 * fb @ wrk                         # (nbf, nbf)
    return project(xhxa, xhxb)


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


# fa/fb sanity
fa, fb = get_fock(X[0])
print(f"\nfa symmetric? max|fa-fa^T|={np.max(np.abs(fa-fa.T)):.2e}  "
      f"fb symmetric? {np.max(np.abs(fb-fb.T)):.2e}")
eA = np.array(mol.data['OQP::E_MO_A'], copy=True).ravel()
print(f"max|diag(fa) - E_MO_A| = {np.max(np.abs(np.diag(fa)-eA)):.2e}  "
      f"(ROHF: occ/virt canonical may differ in SOMO block)")

print("\n=== V1a: Python 2FT(I=I) vs Fortran 2FT contribution (per state) ===")
for st in range(1, nstate+1):
    L = relax_2ft(st-1, st-1)
    R_2ft = dump_sfrorhs(st, ['NAC_ZERO_SP', 'NAC_ZERO_AB1'])         # 2e + 2FT
    R_a = dump_sfrorhs(st, ['NAC_ZERO_SP', 'NAC_ZERO_2FT', 'NAC_ZERO_AB1'])  # 2e only
    Rf = R_2ft - R_a                                                  # pure 2FT
    c = np.dot(L, Rf)/(np.linalg.norm(L)*np.linalg.norm(Rf)+1e-30)
    print(f"  state {st}: |L_2ft|={np.linalg.norm(L):.4e} |Rf_2ft|={np.linalg.norm(Rf):.4e} "
          f"cos={c:+.8f} ratio={np.linalg.norm(L)/(np.linalg.norm(Rf)+1e-30):.6f} "
          f"max|diff|={np.max(np.abs(L-Rf)):.2e}")

print("\n=== interstate 2FT cross RHS (the deficient baseline) ===")
for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    L = relax_2ft(I, J)
    line = f"  pair ({I+1},{J+1}): |L_2ft|={np.linalg.norm(L):.4e}"
    for name, lo, hi in blocks:
        line += f"  {name}={np.linalg.norm(L[lo:hi]):.3e}"
    print(line)
