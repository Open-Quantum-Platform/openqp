"""Phase 11: frozen-Fock matvec evaluator + orbital-rotation FD of X_I^T(d A)X_J.

Step A (correctness): the standalone matvec mrsf_matvec_apply must reproduce
A.X_K so that X_K . (A X_K) = Omega_K (the excitation energy) for every state.

Step B (frozen-Fock orbital gradient): rotate the ROHF orbitals in a single
occ/virt plane, hold the AO Fock fixed (mrsf_matvec_apply rebuilds C^T F_AO C),
and finite-difference z^T A z to get R^matvec_pq = d(z^T A z)/d theta_pq for
z = X1 (diagonal, must agree with the production sfrorhs RHS up to sign) and a
ground mixture.
"""
import oqp                       # MUST precede numpy: numpy's ILP64 LAPACK
from oqp.pyoqp import Runner      # otherwise intercepts oqp's DSYEVD (Huckel fails)
import numpy as np

INP = '/tmp/nactest/H2O_en.inp'
TH = 1e-3

r = Runner(input_file=INP, log='/tmp/nactest/mvrhs.log')
r.run()
mol = r.mol
print("RUN OK", flush=True)

nstate = 3
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2                      # MRSF: 2 SOMOs
C0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = C0.shape[0]
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
nvirb = nbf - nocb
nij = noca * nvirb

Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))
Om = np.array(mol.data['OQP::td_energies'], copy=True).ravel()   # excitation energies

print(f"nbf={nbf} noca={noca} nocb={nocb} nvirb={nvirb} nij={nij}")
print("Omega (Hartree):", np.array2string(Om, precision=6))


def set_bvec(col):
    raw = Xraw.copy().reshape(-1)
    raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def Ax(col):
    """A . col  via the frozen-Fock standalone matvec (current orbitals)."""
    set_bvec(col)
    oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()


# ---- Step A: eigenvalue reproduction (validates the standalone matvec) ----
print("\n=== Step A: X_K . (A X_K) vs Omega_K  (matvec correctness) ===")
okA = True
for k in range(nstate):
    xk = X[k]
    axk = Ax(xk)
    nn = np.dot(xk, xk)
    expk = np.dot(xk, axk) / nn
    resid = np.linalg.norm(axk - Om[k]*xk) / np.sqrt(nn)   # eigenvector residual
    err = abs(expk - Om[k])
    print(f"  state {k+1}: <X|AX>/<X|X> = {expk:.6f}  Omega = {Om[k]:.6f}  "
          f"|diff| = {err:.2e}   ||AX-OmX||/||X|| = {resid:.2e}")
    if resid > 1e-3:
        okA = False
print("Step A:", "PASS -- standalone matvec reproduces the spectrum" if okA
      else "FAIL -- matvec mismatch, debug before FD")

# ---- Step B: frozen-Fock orbital-rotation FD -> R^matvec, vs sfrorhs ----
# ROHF rotation space, matching sfrorhs order: doc-socc, doc-virt, soc-virt.
pairs = ([(i, j) for i in range(nocb, noca) for j in range(0, nocb)]      # doc-socc
         + [(k, j) for k in range(noca, nbf) for j in range(0, nocb)]     # doc-virt
         + [(k, i) for k in range(noca, nbf) for i in range(nocb, noca)]) # soc-virt
nconf = len(pairs)
print(f"\n=== Step B: frozen-Fock FD R^matvec  (nconf={nconf}) ===")


def zAz(col, Crot_a, Crot_b):
    mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(Crot_a)
    mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(Crot_b)
    return float(np.dot(col, Ax(col)))


def Rmatvec(col, axis):
    """axis='row' rotates Python rows (=Fortran MO columns under transpose)."""
    R = np.zeros(nconf)
    for n, (a, b) in enumerate(pairs):
        def rot(C):
            Cr = C.copy()
            if axis == 'row':
                Cr[a, :] = np.cos(TH)*C[a, :] + np.sin(TH)*C[b, :]
                Cr[b, :] = -np.sin(TH)*C[a, :] + np.cos(TH)*C[b, :]
            else:
                Cr[:, a] = np.cos(TH)*C[:, a] + np.sin(TH)*C[:, b]
                Cr[:, b] = -np.sin(TH)*C[:, a] + np.cos(TH)*C[:, b]
            return Cr
        fp = zAz(col, rot(C0), rot(C0b))
        TH2 = -TH
        def rotm(C):
            Cr = C.copy()
            if axis == 'row':
                Cr[a, :] = np.cos(TH2)*C[a, :] + np.sin(TH2)*C[b, :]
                Cr[b, :] = -np.sin(TH2)*C[a, :] + np.cos(TH2)*C[b, :]
            else:
                Cr[:, a] = np.cos(TH2)*C[:, a] + np.sin(TH2)*C[:, b]
                Cr[:, b] = -np.sin(TH2)*C[:, a] + np.cos(TH2)*C[:, b]
            return Cr
        fm = zAz(col, rotm(C0), rotm(C0b))
        R[n] = (fp - fm) / (2*TH)
    mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(C0)   # restore
    mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(C0b)
    return R


import os
TARGET = 2
nds, ndv = nocb*(noca-nocb), (nbf-noca)*nocb     # block sizes: doc-socc, doc-virt
blocks = [('doc-socc', 0, nds), ('doc-virt', nds, nds+ndv), ('soc-virt', nds+ndv, nconf)]

# FD R^matvec FIRST (both axes), clean post-energy state
Rmv = {ax: Rmatvec(X[TARGET-1], ax) for ax in ('row', 'col')}

# dump the production sfrorhs RHS for state 2
mol.data.set_tdhf_target(TARGET)
os.environ['NAC_DUMP_RHS'] = '1'
oqp.tdhf_mrsf_z_vector(mol)
os.environ.pop('NAC_DUMP_RHS')
Rsfr = np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).ravel()[:nconf]
print(f"\nR^sfrorhs(state {TARGET})  |R|={np.linalg.norm(Rsfr):.4e}", flush=True)


def cosr(a, b):
    return (np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30),
            np.linalg.norm(a)/(np.linalg.norm(b)+1e-30))


for ax in ('row', 'col'):
    c, rr = cosr(Rmv[ax], Rsfr)
    print(f"\naxis={ax}: overall cos={c:+.4f} ratio={rr:.3f}")
    for name, lo, hi in blocks:
        bc, br = cosr(Rmv[ax][lo:hi], Rsfr[lo:hi])
        print(f"    {name:9} cos={bc:+.4f} ratio={br:.3f}  "
              f"|mv|={np.linalg.norm(Rmv[ax][lo:hi]):.3e} |sfr|={np.linalg.norm(Rsfr[lo:hi]):.3e}")

