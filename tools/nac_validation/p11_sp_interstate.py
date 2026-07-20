"""Phase 11 SYNTHESIS PoC: build the interstate mrsfsp (channels 1-6) bilinear in
Python from the exported nac_gchan + mo_a, validate it reduces to the production
diagonal mrsfsp at I=J, then add it to the channel-7 interstate bilinear and test
whether the FULL interstate 2e cross moves the rotation-space deficiency toward
the oracle (especially the (1,3) sign flip).

mrsfsp structure (tdhf_mrsf_lib.F90:1775-2002), ca=cb=mo_a:
  for each channel ich, scr2 = (+/-) ca^T @ ch^AO @ cb = (+/-) gchan[ich]
  scr = 2 * scr2 @ cb = 2 * gchan[ich] @ mo_a           (MO x AO-index j)
  then accumulate into xhxa/xhxb with explicit xv(.,.) contractions.
So we need only gchan (exported) and mo_a (exported). No AO channel needed.

INTERSTATE generalization (bilinear in the two amplitude slots: fmrst2 carries
amplitude-1, xv carries amplitude-2):
  hx_sp^{IJ} = 1/2 [ mrsfsp(F=gchan_I, xv=wrk3_J) + mrsfsp(F=gchan_J, xv=wrk3_I) ]
At I=J this is exactly the production diagonal mrsfsp(F=gchan_I, xv=wrk3_I).
"""
import os
import oqp
from oqp.pyoqp import Runner
import numpy as np

INP = '/tmp/nactest/H2O_en.inp'
SQ = 1.0 / np.sqrt(2.0)

r = Runner(input_file=INP, log='/tmp/nactest/p11_sp.log')
r.run(); mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
print("RUN OK (tight int2e_cutoff=1e-20)", flush=True)

noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
moa = np.array(mol.data['OQP::VEC_MO_A'], copy=True)   # (nbf,nbf) col-major: columns=MOs
nbf = moa.shape[0]
nvirb = nbf - nocb; nij = noca * nvirb; nstate = 3
Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); Xshape = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))

# 0-based SOMO/ground slots
ijlr1 = (noca - 1 - nocb - 1) * noca + (noca - 1) - 1
ijlr2 = (noca - nocb - 1) * noca + (noca) - 1
lr1 = nocb        # 0-based SOMO1 row index (Fortran nocb+1)
lr2 = noca - 1    # 0-based SOMO2 row index (Fortran noca)
print(f"nbf={nbf} noca={noca} nocb={nocb} nvirb={nvirb} lr1={lr1} lr2={lr2}")

prs = ([(i, j) for i in range(nocb, noca) for j in range(0, nocb)]
       + [(k, j) for k in range(noca, nbf) for j in range(0, nocb)]
       + [(k, i) for k in range(noca, nbf) for i in range(nocb, noca)])
nconf = len(prs)
nds, ndv = nocb * (noca - nocb), (nbf - noca) * nocb
blocks = [('doc-socc', 0, nds), ('doc-virt', nds, nds + ndv), ('soc-virt', nds + ndv, nconf)]


def set_bvec(col):
    raw = Xraw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def matvec(col):
    """run matvec, return (gmo (nbf,nbf), gchan (nbf,nbf,6))."""
    set_bvec(col); oqp.mrsf_matvec_apply(mol)
    g = np.array(mol.data['OQP::nac_gmo'], copy=True).ravel().reshape((nbf, nbf), order='F')
    gc = np.array(mol.data['OQP::nac_gchan'], copy=True).ravel().reshape((nbf, nbf, 6), order='F')
    return g, gc


def mrsfxvec(col):
    f = col.copy(); f[ijlr1] = col[ijlr1] * SQ; f[ijlr2] = -col[ijlr1] * SQ
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


def mrsfsp_py(gchan, xv):
    """Replicate mrsfsp(xhxa,xhxb, mo_a, mo_a, xv, fmrst2) using exported gchan.
    gchan[:,:,ich] = mo_a^T @ ch^AO @ mo_a   (ich 0..5 = ado2v,ado1v,adco1,adco2,ao21v,aco12).
    Fortran: scr2 = ca^T @ ch (half-transform, MO x AO); scr = 2*scr2 @ cb.
    With full cb -> scr = 2 * gchan[ich] (pure MO).
    With cb(:,lr) -> scr(:,1) = 2 * gchan[ich][:,lr]."""
    ado2v, ado1v, adco1, adco2, ao21v, aco12 = (gchan[:, :, k] for k in range(6))
    xhxa = np.zeros((nbf, nbf)); xhxb = np.zeros((nbf, nbf))

    def scr_of(ch, sign):
        return 2.0 * sign * ch   # 2*scr2@cb = 2*sign*gchan  (full cb)

    # ---- xhxa ----
    # o1v
    scr = scr_of(ao21v, +1.0)
    for j in range(noca, nbf):
        xhxa[:, lr2] += scr[:, j] * xv[lr1, j]
        xhxa[:, lr1] -= scr[:, j] * xv[lr2, j]
    # co1
    scr = scr_of(aco12, -1.0)
    for i in range(0, nocb):
        xhxa[:, i] += scr[:, lr2] * xv[i, lr1]
        xhxa[:, i] -= scr[:, lr1] * xv[i, lr2]
    # adco2
    scr = scr_of(adco2, +1.0)
    for j in range(noca, nbf):
        xhxa[:, lr1] += scr[:, j] * xv[lr1, j]
    # adco1
    scr = scr_of(adco1, +1.0)
    for j in range(noca, nbf):
        xhxa[:, lr2] += scr[:, j] * xv[lr2, j]
    # ado2v with cb(:,lr1): scr(:,1) = 2*gchan[ado2v][:,lr1]
    scrv = 2.0 * ado2v[:, lr1]
    for i in range(0, nocb):
        xhxa[:, i] += scrv * xv[i, lr1]
    # ado1v with cb(:,lr2)
    scrv = 2.0 * ado1v[:, lr2]
    for i in range(0, nocb):
        xhxa[:, i] += scrv * xv[i, lr2]

    # ---- xhxb ----
    scr = scr_of(ao21v, +1.0)
    for j in range(noca, nbf):
        xhxb[:, j] += scr[lr2, :] * xv[lr1, j]
    scr = scr_of(ao21v, -1.0)
    for j in range(noca, nbf):
        xhxb[:, j] += scr[lr1, :] * xv[lr2, j]
    scr = scr_of(aco12, -1.0)
    for i in range(0, nocb):
        xhxb[:, lr2] += scr[i, :] * xv[i, lr1]
    scr = scr_of(aco12, +1.0)
    for i in range(0, nocb):
        xhxb[:, lr1] += scr[i, :] * xv[i, lr2]
    scr = scr_of(ado2v, +1.0)
    for i in range(0, nocb):
        xhxb[:, lr1] += scr[i, :] * xv[i, lr1]
    scr = scr_of(ado1v, +1.0)
    for i in range(0, nocb):
        xhxb[:, lr2] += scr[i, :] * xv[i, lr2]
    # O1V: scr2 = ca(:,lr1)^T @ adco2 (1 x AO); scr = 2*scr2@cb (1 x MO) = 2*gchan[adco2][lr1,:]
    scrr = 2.0 * adco2[lr1, :]
    for j in range(noca, nbf):
        xhxb[:, j] += scrr * xv[lr1, j]
    scrr = 2.0 * adco1[lr2, :]
    for j in range(noca, nbf):
        xhxb[:, j] += scrr * xv[lr2, j]
    return xhxa, xhxb


def rhs_full_diag(col):
    """ch7 + mrsfsp, diagonal -- must == sfrorhs (2FT,AB1 off, SP ON)."""
    g, gc = matvec(col); w = iatogen(mrsfxvec(col))
    hxa = 2.0 * g @ w[0:noca, :].T
    hxb = 2.0 * g[0:noca, :].T @ w[0:noca, :]
    spa, spb = mrsfsp_py(gc, w)
    hxa[:, :noca] += spa[:, :noca]; hxb += spb
    return project(hxa, hxb)


def dump_sfrorhs(target, gates):
    mol.data.set_tdhf_target(target)
    allg = ('NAC_ZERO_2FT', 'NAC_ZERO_AB1', 'NAC_ZERO_HX', 'NAC_ZERO_T', 'NAC_ZERO_SP')
    for k in allg:
        os.environ.pop(k, None)
    for k in gates:
        os.environ[k] = '1'
    os.environ['NAC_DUMP_RHS'] = '1'
    oqp.tdhf_mrsf_z_vector(mol)
    os.environ.pop('NAC_DUMP_RHS', None)
    for k in allg:
        os.environ.pop(k, None)
    return np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).ravel()[:nconf]


def cmp(a, b):
    return (np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-30),
            np.linalg.norm(a) / (np.linalg.norm(b) + 1e-30),
            np.max(np.abs(a - b)))


# ============================================================================
print("\n=== VALIDATION 1: Python mrsfsp_py reproduces production diagonal (SP on, 2FT/AB1 off) ===")
for st in range(1, nstate + 1):
    Rpy = rhs_full_diag(X[st - 1])
    Rsf = dump_sfrorhs(st, ['NAC_ZERO_2FT', 'NAC_ZERO_AB1'])  # SP kept
    c, rr, mx = cmp(Rpy, Rsf)
    print(f"  state {st}: cos={c:+.8f} ratio={rr:.6f} max|diff|={mx:.2e}  |Rpy|={np.linalg.norm(Rpy):.4e}")
    for nm, lo, hi in blocks:
        bc, br, bmx = cmp(Rpy[lo:hi], Rsf[lo:hi])
        print(f"     {nm:9} cos={bc:+.8f} ratio={br:.6f} max|diff|={bmx:.2e}")


# ============================================================================
# Build interstate 2e cross (ch7 + sp), test I=J reduction + off-diagonal asymmetry
# ============================================================================
def hx_inter(gI, gcI, wI, gJ, gcJ, wJ):
    """full interstate 2e hxa/hxb = ch7 bilinear + 1/2[sp(F_I,w_J)+sp(F_J,w_I)]."""
    hxa = gI @ wJ[0:noca, :].T + gJ @ wI[0:noca, :].T
    hxb = gI[0:noca, :].T @ wJ[0:noca, :] + gJ[0:noca, :].T @ wI[0:noca, :]
    spa1, spb1 = mrsfsp_py(gcI, wJ)
    spa2, spb2 = mrsfsp_py(gcJ, wI)
    hxa[:, :noca] += 0.5 * (spa1[:, :noca] + spa2[:, :noca])
    hxb += 0.5 * (spb1 + spb2)
    return hxa, hxb


def chain_cross_2e(I, J):
    """gradient-chain 2e cross via polarization: A2e(X_I,X_J)=RZ-1/2(RI+RJ),
    where R_K = the production diagonal 2e RHS (ch7+sp) for amplitude K at target I.
    This is what the deficient gradient chain produces off-diagonal (symmetric in I,J
    only because each R uses mrsfsp(F(amp),wrk3(amp)) -- DIAGONAL amplitude)."""
    Z = (X[I] + X[J]) / np.sqrt(2.0)
    return rhs_full_diag(Z) - 0.5 * (rhs_full_diag(X[I]) + rhs_full_diag(X[J]))


print("\n=== VALIDATION 2: interstate (ch7+sp) reduces to diagonal at I=J ===")
g_, gc_, w_ = [], [], []
for st in range(nstate):
    g, gc = matvec(X[st]); g_.append(g); gc_.append(gc); w_.append(iatogen(mrsfxvec(X[st])))
for st in range(nstate):
    hxa, hxb = hx_inter(g_[st], gc_[st], w_[st], g_[st], gc_[st], w_[st])
    Rii = project(hxa, hxb)
    Rdi = rhs_full_diag(X[st])
    print(f"  state {st+1}: max|R_inter(I=I) - R_diag| = {np.max(np.abs(Rii-Rdi)):.2e}")

print("\n=== DISCRIMINATOR: interstate (ch7+sp) cross vs gradient-chain 2e cross ===")
print("    (the deficiency must show as cos<1 / sign flip; (1,3) is the key)")
for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    hxa, hxb = hx_inter(g_[I], gc_[I], w_[I], g_[J], gc_[J], w_[J])
    Rinter = project(hxa, hxb)
    Rchain = chain_cross_2e(I, J)
    c = np.dot(Rinter, Rchain) / (np.linalg.norm(Rinter) * np.linalg.norm(Rchain) + 1e-30)
    d = Rinter - Rchain
    print(f"  pair ({I+1},{J+1}): |inter|={np.linalg.norm(Rinter):.4e} |chain|={np.linalg.norm(Rchain):.4e} "
          f"cos={c:+.5f} |inter-chain|={np.linalg.norm(d):.4e}")
    for nm, lo, hi in blocks:
        cb = np.dot(Rinter[lo:hi], Rchain[lo:hi]) / (np.linalg.norm(Rinter[lo:hi]) * np.linalg.norm(Rchain[lo:hi]) + 1e-30)
        print(f"     {nm:9} |inter|={np.linalg.norm(Rinter[lo:hi]):.3e} |chain|={np.linalg.norm(Rchain[lo:hi]):.3e} cos={cb:+.5f}")

# isolate the sp contribution to the cross (ch7 is already symmetric->same as chain ch7)
print("\n=== sp-only interstate cross asymmetry: sp(F_I,w_J) vs sp(F_J,w_I) ===")
for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    spa1, spb1 = mrsfsp_py(gc_[I], w_[J])
    spa2, spb2 = mrsfsp_py(gc_[J], w_[I])
    R1 = project(np.pad(spa1[:, :noca], ((0,0),(0,nbf-noca))), spb1)
    R2 = project(np.pad(spa2[:, :noca], ((0,0),(0,nbf-noca))), spb2)
    c = np.dot(R1, R2) / (np.linalg.norm(R1) * np.linalg.norm(R2) + 1e-30)
    print(f"  pair ({I+1},{J+1}): |sp(F_I,w_J)|={np.linalg.norm(R1):.4e} |sp(F_J,w_I)|={np.linalg.norm(R2):.4e} "
          f"cos={c:+.5f}  asymmetry|R1-R2|={np.linalg.norm(R1-R2):.4e}")


# ============================================================================
# CRUX TEST: the SYMMETRIC interstate 2e cross == gradient chain (proven above).
# The only genuinely-new object is the ANTISYMMETRIC sp difference, which the
# polarization of the diagonal DISCARDS. Does it align with the deficiency?
# We compare in ROTATION space against the FULL chain 2e cross direction, and
# report the antisymmetric sp vector for the nuclear-projection seam.
# ============================================================================
print("\n=== ANTISYMMETRIC interstate sp (the object the gradient chain CANNOT produce) ===")
def sp_cross_pieces(I, J):
    spa1, spb1 = mrsfsp_py(gc_[I], w_[J])   # F_I, w_J
    spa2, spb2 = mrsfsp_py(gc_[J], w_[I])   # F_J, w_I
    R1 = project(np.pad(spa1[:, :noca], ((0,0),(0,nbf-noca))), spb1)
    R2 = project(np.pad(spa2[:, :noca], ((0,0),(0,nbf-noca))), spb2)
    return 0.5*(R1+R2), 0.5*(R1-R2)   # sym, antisym

for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    sym, anti = sp_cross_pieces(I, J)
    print(f"  pair ({I+1},{J+1}): |sp_sym|={np.linalg.norm(sym):.4e}  |sp_antisym|={np.linalg.norm(anti):.4e}")
    for nm, lo, hi in blocks:
        print(f"     {nm:9} |sym|={np.linalg.norm(sym[lo:hi]):.3e} |anti|={np.linalg.norm(anti[lo:hi]):.3e}")

# Also: does the FULL 2e bilinear (ch7+sp) symmetric reproduce qz coded 2e?
# Save the antisymmetric sp vectors for the nuclear seam test.
np.savez('/tmp/nactest/p11_sp_antisym.npz',
         **{f'sym_{I+1}{J+1}': sp_cross_pieces(I, J)[0] for (I, J) in [(0,1),(0,2),(1,2)]},
         **{f'anti_{I+1}{J+1}': sp_cross_pieces(I, J)[1] for (I, J) in [(0,1),(0,2),(1,2)]})
print("\nsaved /tmp/nactest/p11_sp_antisym.npz (sym + antisym rotation-space sp cross)")
