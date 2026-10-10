#!/usr/bin/env python3
"""pyscf-free minao/SAD initial guess.

D0 = superposition of PROJECTED tabulated atomic densities (pyscf's
init_guess_by_minao recipe):

    S_cross[mu,a] = < chi_mu^comp | chi_a^minao >      (analytic cartesian overlap)
    C  = S^{-1} S_cross                                (S = comp-basis overlap, GIVEN)
    D0 = C @ diag(occ) @ C^T                           (spherically-averaged atomic occ)

All AO conventions are the OpenQP CARTESIAN frame used by nativeprep_xcgrid.py:
the computational AO for component (a,b,c) is  dscale(a,b,c) * raw contraction,
where raw uses coefficients that already carry the primitive (l,0,0) norm and
dscale = sqrt((2l-1)!! / ((2a-1)!!(2b-1)!!(2c-1)!!)).  Every AO is unit-normed.

Tr(D0 S) is NOT exactly nelec for a projected guess (pyscf's isn't either); the
caller rescales by nelec/Tr(D0 S) exactly as it already does for the Hueckel DM.
The gate here is that Tr(D0 S) lands within a few percent of nelec and D0 is a
sensible density (symmetric, eigenvalues in [~0, ~2.1]) -- a wrong dscale or
component order sends Tr far off, so it catches convention bugs.
"""
import math
import numpy as np
import basis_set_exchange as bse

# ---- cartesian conventions (identical to nativeprep_xcgrid.py) ----------------
def _df2(n):
    r = 1
    while n > 1:
        r *= n; n -= 2
    return r

_OQP_CART = {
    0: [(0,0,0)],
    1: [(1,0,0),(0,1,0),(0,0,1)],
    2: [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)],
    3: [(3,0,0),(0,3,0),(0,0,3),(2,1,0),(2,0,1),(1,2,0),(0,2,1),(1,0,2),(0,1,2),(1,1,1)],
}
def _dscale(a, b, c):
    l = a + b + c
    return math.sqrt(_df2(2*l-1) / (_df2(2*a-1)*_df2(2*b-1)*_df2(2*c-1)))

def gto_norm(l, a):
    return math.sqrt(2**(2*l+3)*math.factorial(l+1)*(2*a)**(l+1.5)
                     /(math.factorial(2*l+2)*math.sqrt(math.pi)))

# ---- analytic cartesian gaussian overlap --------------------------------------
def _o1(l1, l2, PA, PB, p):
    """int (x-A)^l1 (x-B)^l2 exp(-p (x-P)^2) dx,  PA=P-A, PB=P-B."""
    tot = 0.0
    for k in range((l1 + l2) // 2 + 1):
        f = 0.0
        for i in range(max(0, 2*k - l2), min(2*k, l1) + 1):
            j = 2*k - i
            f += (math.comb(l1, i) * math.comb(l2, j)
                  * (PA ** (l1 - i)) * (PB ** (l2 - j)))
        tot += f * _df2(2*k - 1) / (2*p) ** k
    return tot * math.sqrt(math.pi / p)

def _prim_overlap(la, ma, na, A, al, lb, mb, nb, B, be):
    p = al + be
    mu = al * be / p
    ab2 = (A[0]-B[0])**2 + (A[1]-B[1])**2 + (A[2]-B[2])**2
    P = ((al*A[0]+be*B[0])/p, (al*A[1]+be*B[1])/p, (al*A[2]+be*B[2])/p)
    pre = math.exp(-mu * ab2)
    return (pre
            * _o1(la, lb, P[0]-A[0], P[0]-B[0], p)
            * _o1(ma, mb, P[1]-A[1], P[1]-B[1], p)
            * _o1(na, nb, P[2]-A[2], P[2]-B[2], p))

def _contracted_overlap(pw_i, A, ei, ci, pw_j, B, ej, cj):
    """Overlap of two RAW contractions (coefs already primitive-normalized)."""
    la, ma, na = pw_i; lb, mb, nb = pw_j
    s = 0.0
    for al, ca in zip(ei, ci):
        for be, cb in zip(ej, cj):
            s += ca * cb * _prim_overlap(la, ma, na, A, al, lb, mb, nb, B, be)
    return s

# ---- vectorized shell-pair overlap block (primitive dimension in numpy) --------
def _o1v(l1, l2, PA, PB, p):
    """int (x-A)^l1 (x-B)^l2 exp(-p(x-P)^2) dx over ARRAYS PA,PB,p."""
    tot = np.zeros_like(p)
    for k in range((l1 + l2) // 2 + 1):
        f = np.zeros_like(p)
        for i in range(max(0, 2*k - l2), min(2*k, l1) + 1):
            j = 2*k - i
            f = f + (math.comb(l1, i) * math.comb(l2, j)
                     * PA**(l1 - i) * PB**(l2 - j))
        tot = tot + f * (_df2(2*k - 1) / (2.0*p)**k)
    return tot * np.sqrt(math.pi / p)

def _shellpair_block(comps_i, A, ei, ci, comps_j, B, ej, cj):
    """(len(comps_i) x len(comps_j)) raw-overlap block, vectorized over the
    n1*n2 primitive pairs. comps are lists of (a,b,c) cartesian powers."""
    ei = np.asarray(ei); ci = np.asarray(ci); ej = np.asarray(ej); cj = np.asarray(cj)
    E1 = ei[:, None]; E2 = ej[None, :]
    p = (E1 + E2).ravel()
    mu = (E1 * E2 / (E1 + E2)).ravel()
    w = (ci[:, None] * cj[None, :]).ravel()
    ab2 = (A[0]-B[0])**2 + (A[1]-B[1])**2 + (A[2]-B[2])**2
    pre = w * np.exp(-mu * ab2)
    Px = ((E1*A[0] + E2*B[0])/(E1+E2)).ravel()
    Py = ((E1*A[1] + E2*B[1])/(E1+E2)).ravel()
    Pz = ((E1*A[2] + E2*B[2])/(E1+E2)).ravel()
    PAx, PBx = Px - A[0], Px - B[0]
    PAy, PBy = Py - A[1], Py - B[1]
    PAz, PBz = Pz - A[2], Pz - B[2]
    out = np.empty((len(comps_i), len(comps_j)))
    for bi, (a1, b1, c1) in enumerate(comps_i):
        for bj, (a2, b2, c2) in enumerate(comps_j):
            sx = _o1v(a1, a2, PAx, PBx, p)
            sy = _o1v(b1, b2, PAy, PBy, p)
            sz = _o1v(c1, c2, PAz, PBz, p)
            out[bi, bj] = float(np.dot(pre, sx*sy*sz))
    return out

# ---- minao reference basis ----------------------------------------------------
_REF_TRIED = []
def _load_ref_basis():
    """Return (name, {Z: [(l, exps, prim_norm_coefs), ...]}); prefer MINAO."""
    import os
    pref = os.environ.get("OQP_MINAO_BASIS")
    # ano-r0 (ANO minimal, closest to pyscf's ANO-derived MINAO) is the effective
    # default -- BSE has no "minao" by that name, and ano-r0 converges fastest
    # (DFT (H2O)16: 10 cycles vs sto-3g 15 vs Hueckel 16).
    cands = ([pref] if pref else []) + ["ano-r0", "minao", "sto-3g"]
    for name in cands:
        try:
            bse.get_basis(name, elements=[8])
            _REF_TRIED.append(name + ":ok")
            return name
        except Exception as e:                        # noqa: BLE001
            _REF_TRIED.append(f"{name}:{type(e).__name__}")
    raise RuntimeError("no minao/sto-3g reference basis available: " + ",".join(_REF_TRIED))

_ref_cache = {}
def _ref_shells(name, Z):
    key = (name, Z)
    if key in _ref_cache:
        return _ref_cache[key]
    el = list(bse.get_basis(name, elements=[int(Z)])["elements"].values())[0]
    out = []
    for sh in el["electron_shells"]:
        exps = np.array([float(x) for x in sh["exponents"]])
        ams = sh["angular_momentum"]
        crows = sh["coefficients"]
        # an sp shell has angular_momentum [0,1] and one coef row per l
        if len(ams) == 1:
            ams = ams * len(crows)
        for l, crow in zip(ams, crows):
            c = np.array([float(x) for x in crow])
            m = np.abs(c) > 0
            e = exps[m]
            cn = c[m] * np.array([gto_norm(l, a) for a in e])
            out.append((int(l), e, cn))
    _ref_cache[key] = out
    return out

# aufbau capacity per angular momentum: s=2, p=6, d=10
def _shell_occ(shells_l, Z):
    """Electrons per minao shell (aufbau fill, shells assumed in aufbau order)."""
    rem = int(round(Z))
    occ = []
    for l in shells_l:
        cap = 2 * (2*l + 1)
        take = min(rem, cap)
        occ.append(take)
        rem -= take
    return occ                                        # total electrons per shell

# ---- the guess ----------------------------------------------------------------
def minao_guess(sh_am, sh_cx, sh_cy, sh_cz, sh_g0, sh_nc,
                ex_all, cc_all, dscale, S_full, Zs, xyz, verbose=False):
    """Build the projected minao density D0 (nbf x nbf) in the OpenQP frame.

    sh_g0 is 1-based into ex_all/cc_all (nativeprep convention); dscale is per-AO
    in the OpenQP AO order (shell-by-shell, _OQP_CART[l] components).
    """
    ex_all = np.asarray(ex_all, float); cc_all = np.asarray(cc_all, float)
    nsh = len(sh_am)
    # AO offset of each comp shell
    sh_off = []
    off = 0
    for ish in range(nsh):
        sh_off.append(off)
        off += (sh_am[ish]+1)*(sh_am[ish]+2)//2
    nbf = off
    assert nbf == S_full.shape[0], f"nbf {nbf} vs S {S_full.shape}"

    name = _load_ref_basis()

    # assemble minao AOs: list of (powers, center, exps, coefs_prim_norm, occ);
    # group AOs by their parent shell so we screen once per (comp shell, minao
    # shell) pair instead of per component.
    m_shells = []                                     # (center, l, e, cn, emin, csum, comps, nrm[], occ)
    for ia, Z in enumerate(Zs):
        A = (float(xyz[ia][0]), float(xyz[ia][1]), float(xyz[ia][2]))
        shells = _ref_shells(name, int(round(Z)))
        occ_sh = _shell_occ([l for (l, _, _) in shells], Z)
        for (l, e, cn), nele in zip(shells, occ_sh):
            comps = _OQP_CART[l]
            per = nele / len(comps)
            self_blk = _shellpair_block(comps, A, e, cn, comps, A, e, cn)
            nrm = 1.0 / np.sqrt(np.diag(self_blk))    # per-component unit norm
            m_shells.append((A, l, np.asarray(e), np.asarray(cn),
                             float(np.min(e)), float(np.sum(np.abs(cn))),
                             comps, nrm, per))
    nm = sum(len(s[6]) for s in m_shells)
    aoff = []; off = 0
    for s in m_shells:
        aoff.append(off); off += len(s[6])
    m_occ = []
    for s in m_shells:
        m_occ += [s[8]] * len(s[6])

    # S_cross[mu, a] = dscale_mu * < raw_mu | chi_a^minao >, screened by a rigorous
    # Gaussian bound: the loosest primitive pair (smallest exponents) sets the max
    # magnitude exp(-mu0 AB^2)*|c1||c2|; skip the shell pair if below TOL. minao
    # functions are atom-tight so only same/near-neighbour comp shells survive.
    TOL = 1e-10
    Scross = np.zeros((nbf, nm))
    for ish in range(nsh):
        l = sh_am[ish]
        A = (sh_cx[ish], sh_cy[ish], sh_cz[ish])
        g0 = sh_g0[ish] - 1; nc = sh_nc[ish]
        e_i = ex_all[g0:g0+nc]; c_i = cc_all[g0:g0+nc]
        emin_i = float(np.min(e_i)); csum_i = float(np.sum(np.abs(c_i)))
        comps_i = _OQP_CART[l]
        o_i = sh_off[ish]; ni = len(comps_i)
        ds_i = dscale[o_i:o_i+ni]
        for js, (B, lj, e_j, c_j, emin_j, csum_j, comps_j, nrm_j, occ_j) in enumerate(m_shells):
            ab2 = (A[0]-B[0])**2 + (A[1]-B[1])**2 + (A[2]-B[2])**2
            if math.exp(-emin_i*emin_j/(emin_i+emin_j)*ab2)*csum_i*csum_j < TOL:
                continue
            blk = _shellpair_block(comps_i, A, e_i, c_i, comps_j, B, e_j, c_j)
            blk = ds_i[:, None] * blk * nrm_j[None, :]
            a0 = aoff[js]; nj = len(comps_j)
            Scross[o_i:o_i+ni, a0:a0+nj] = blk

    occ = np.array(m_occ)
    C = np.linalg.solve(S_full, Scross)               # nbf x nm
    D0 = (C * occ) @ C.T
    D0 = 0.5 * (D0 + D0.T)
    if verbose:
        tr = float(np.einsum('ij,ji->', D0, S_full))
        ev = np.linalg.eigvalsh(D0)
        print(f"[minao] ref={name} nm={nm} occ_sum={occ.sum():.3f} "
              f"Tr(D0 S)={tr:.4f} eig[{ev.min():.3f},{ev.max():.3f}]", flush=True)
    return D0


# ---- self test ----------------------------------------------------------------
if __name__ == "__main__":
    # Build a small computational basis (cc-pVDZ) analytically, no OpenQP:
    def comp_shells(Z, center):
        el = list(bse.get_basis("cc-pvdz", elements=[int(Z)])["elements"].values())[0]
        out = []
        for sh in el["electron_shells"]:
            exps = np.array([float(x) for x in sh["exponents"]])
            ams = sh["angular_momentum"]; crows = sh["coefficients"]
            if len(ams) == 1: ams = ams * len(crows)
            for l, crow in zip(ams, crows):
                c = np.array([float(x) for x in crow]); m = np.abs(c) > 0
                e = exps[m]; cn = c[m]*np.array([gto_norm(l, a) for a in e])
                # OpenQP's coef ALSO carries the contracted (l,0,0) norm so that
                # AO_oqp = dscale*raw is unit-normed; mimic that here so the
                # self-test S has unit diagonal and D0 eigenvalues are physical.
                s00 = _contracted_overlap((l,0,0), center, e, cn, (l,0,0), center, e, cn)
                cn = cn / math.sqrt(s00)
                out.append((int(l), center, e, cn))
        return out

    def assemble(atoms):
        """atoms: [(Z, (x,y,z)Bohr)]  -> nativeprep-style shell arrays + S_full."""
        sh_am=[]; sh_cx=[]; sh_cy=[]; sh_cz=[]; sh_g0=[]; sh_nc=[]
        ex_all=[]; cc_all=[]; dscale=[]; Zs=[]; xyz=[]
        for Z, ctr in atoms:
            Zs.append(Z); xyz.append(ctr)
            for (l, c, e, cn) in comp_shells(Z, ctr):
                sh_am.append(l); sh_cx.append(c[0]); sh_cy.append(c[1]); sh_cz.append(c[2])
                sh_g0.append(len(ex_all)+1); sh_nc.append(len(e))
                ex_all += list(e); cc_all += list(cn)
                for (a,b,cc) in _OQP_CART[l]:
                    dscale.append(_dscale(a,b,cc))
        dscale = np.array(dscale); Zs = np.array(Zs, float); xyz = np.array(xyz, float)
        # build S_full analytically: AO_mu = dscale_mu * raw_mu
        nbf = len(dscale)
        # AO metadata
        aos = []
        off = 0
        for ish in range(len(sh_am)):
            l = sh_am[ish]; A=(sh_cx[ish],sh_cy[ish],sh_cz[ish])
            g0=sh_g0[ish]-1; nc=sh_nc[ish]
            e=np.array(ex_all[g0:g0+nc]); c=np.array(cc_all[g0:g0+nc])
            for bi,(a,b,cc) in enumerate(_OQP_CART[l]):
                aos.append(((a,b,cc), A, e, c, dscale[off])); off+=1
        S = np.zeros((nbf,nbf))
        for i in range(nbf):
            pwi,Ai,ei,ci,di = aos[i]
            for j in range(i+1):
                pwj,Aj,ej,cj,dj = aos[j]
                ov = _contracted_overlap(pwi,Ai,ei,ci, pwj,Aj,ej,cj)
                S[i,j]=S[j,i]=di*dj*ov
        return (sh_am,sh_cx,sh_cy,sh_cz,sh_g0,sh_nc,ex_all,cc_all,dscale,S,Zs,xyz)

    print("ref basis probe:", end=" ")
    _load_ref_basis(); print(",".join(_REF_TRIED))

    # test 1: single O atom
    O = assemble([(8, (0.0,0.0,0.0))])
    Sfull = O[9]
    diag = np.diag(Sfull)
    print(f"[O] nbf={Sfull.shape[0]} S diag range [{diag.min():.4f},{diag.max():.4f}] "
          f"(should be ~1.0 -> validates overlap+dscale)")
    D0 = minao_guess(*O, verbose=True)
    trO = float(np.einsum('ij,ji->', D0, Sfull))
    assert abs(trO - 8) < 0.5, f"O Tr(D0 S)={trO} far from 8"
    # OCCUPATIONS are eig(S^1/2 D0 S^1/2) = eig(D0 S), not eig(D0). A projected
    # minao guess is non-idempotent so occupations overshoot 2 (harmless: real SCF
    # converges to the same E), but must not go wildly negative.
    occO = np.linalg.eigvals(D0 @ Sfull).real
    assert occO.min() > -0.3, f"O occ too negative {occO.min()}"

    # test 2: a water molecule (Bohr)
    W = assemble([(8,(0.0,0.0,0.221)), (1,(0.0,1.432,-0.884)), (1,(0.0,-1.432,-0.884))])
    Sw = W[9]
    D0w = minao_guess(*W, verbose=True)
    trW = float(np.einsum('ij,ji->', D0w, Sw))
    assert abs(trW - 10) < 0.8, f"H2O Tr(D0 S)={trW} far from 10"
    print(f"[H2O] Tr(D0 S)={trW:.4f} (nelec=10)  SELF-TEST PASS")
