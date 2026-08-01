"""nac_formula_kernel: the verified MRSF state-overlap formula kernel.

Library form of the campaign's machine-precision machinery
(MRSF_NAC_DERIVATION.md 7.9-7.21):
  build_context(mol) -> ctx with amplitudes, dims, and masks
  replica_S(ctx, M)  -> exact Python replica of compute_states_overlap
                        over exact tlf=0 minors (== Fortran to ~2e-16)
  fortran_S(mol, ctx, theta, K) -> Fortran S at a staged MO rotation
                        (transpose-corrected staging + read)
  gamma_closed(ctx)  -> gamma^formula[I,J,p,q] in ONE linear-algebra pass
                        (cofactor sensitivities; == generator sweep to 1e-13)

Conventions (proven; see the derivation log): numpy 2-D tagarrays are the
TRANSPOSE of the Fortran matrices; the MO rotation C -> C e^{tK} is staged
as W -> expm(-tK) @ W for antisymmetric K (expm(+tK) @ W for symmetric).
"""
import numpy as np
import numpy.linalg as la


def build_context(mol):
    import oqp  # noqa: F401
    nstate = mol.config['tdhf']['nstate']
    noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
    nocb = noca - 2
    W = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    nbf = W.shape[0]
    nvirb = nbf - nocb
    nij = noca * nvirb
    noc = noca - 1
    RS = 1.0 / np.sqrt(2.0)
    X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
    X0 = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()

    def unfold(st):
        ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
        ijlr2 = (noca - nocb - 1) * noca + noca
        x = np.zeros((noca, nvirb))
        for i in range(1, noca + 1):
            for jj in range(nocb + 1, nbf + 1):
                ij = (jj - nocb - 1) * noca + i
                if ij == ijlr1:
                    x[i - 1, jj - nocb - 1] = X0[ijlr1 - 1, st - 1] * RS
                elif ij == ijlr2:
                    x[i - 1, jj - nocb - 1] = -X0[ijlr1 - 1, st - 1] * RS
                else:
                    x[i - 1, jj - nocb - 1] = X0[ij - 1, st - 1]
        return x

    Xt = [unfold(s + 1) for s in range(nstate)]
    genmask = np.ones((noca, nvirb))
    genmask[nocb:noca, 0:2] = 0.0
    return dict(nstate=nstate, noca=noca, nocb=nocb, nbf=nbf, nvirb=nvirb,
                nij=nij, noc=noc, RS=RS, W=W, X0_raw=X0_raw, Xt=Xt,
                genmask=genmask,
                Wb=np.array(mol.data['OQP::VEC_MO_B'], copy=True))


# ---------------- exact minors (literal ov_exact replicas) -----------------
def s_ij_of(ctx, M):
    noca = ctx['noca']
    G = np.zeros((noca, noca))
    for i1 in range(1, noca + 1):
        for i2 in range(1, noca + 1):
            if i1 == i2:
                keep = [k for k in range(noca) if k != i1 - 1]
                G[i1 - 1, i2 - 1] = la.det(M[np.ix_(keep, keep)])
            else:
                imin, imax = min(i1, i2), max(i1, i2)
                rows = ([k for k in range(noca) if k + 1 not in (imin, imax)]
                        + [i2 - 1])
                cols = ([k for k in range(noca) if k + 1 not in (imin, imax)]
                        + [i1 - 1])
                G[i1 - 1, i2 - 1] = -la.det(M[np.ix_(rows, cols)])
    return G


def s_ab_of(ctx, M):
    nvirb, nocb = ctx['nvirb'], ctx['nocb']
    G = np.zeros((nvirb, nvirb))
    core = list(range(nocb))
    for j1 in range(nvirb):
        for j2 in range(nvirb):
            G[j1, j2] = la.det(M[np.ix_(core + [nocb + j1],
                                        core + [nocb + j2])])
    return G


def _ia_maps(ctx, i1, j1):
    """NET row/col index maps of the LITERAL ov_exact case(3) layout
    (blocks 3/4 always win ddet rows noc-1, noc; block2 shifts +1 past i1)."""
    noc, nocb = ctx['noc'], ctx['nocb']
    rows = [((r if r <= i1 - 1 else r + 1) - 1) for r in range(1, noc - 1)]
    rows += [noc - 1, noc]                # 0-based of s_mo rows noc, noc+1
    cols = [((c if c <= i1 - 1 else c + 1) - 1) for c in range(1, noc - 1)]
    cols += [i1 - 1, nocb + j1]
    return rows, cols


def s_ia_of(ctx, M):
    noca, nvirb = ctx['noca'], ctx['nvirb']
    G = np.zeros((noca, nvirb))
    for i1 in range(1, noca + 1):
        for j1 in range(nvirb):
            rows, cols = _ia_maps(ctx, i1, j1)
            G[i1 - 1, j1] = la.det(M[np.ix_(rows, cols)])
    return G


# ---------------- the contraction (exact replica) --------------------------
def contraction(ctx, s_ij, s_ab, s_ia, xo, xn):
    ns = len(xo)
    noca, nocb, nbf = ctx['noca'], ctx['nocb'], ctx['nbf']
    RS, genmask = ctx['RS'], ctx['genmask']
    S = np.zeros((ns, ns))
    for oi in range(ns):
        for ni in range(ns):
            co, cn = xo[oi], xn[ni]
            cog, cng = co * genmask, cn * genmask
            acc = float(np.sum((cog @ s_ab) * (s_ij @ cng)))
            acc += float(np.sum((cog @ s_ia.T) * (cng @ s_ia.T).T))
            for pi in range(nocb, noca):
                for qi in range(nocb, noca):
                    for ri in range(nocb, noca):
                        for si in range(nocb, noca):
                            acc += (co[pi, qi - nocb] * s_ij[pi, ri]
                                    * cn[ri, si - nocb]
                                    * s_ab[qi - nocb, si - nocb])
            for pi in range(nocb, noca):
                for qi in range(nocb, noca):
                    for ri in range(noca):
                        for si in range(nocb, nbf):
                            if (ri >= nocb) and (si < noca):
                                continue
                            acc += (co[pi, qi - nocb] * cn[ri, si - nocb]
                                    * (s_ij[pi, ri] * s_ab[qi - nocb, si - nocb]
                                       + s_ia[pi, si - nocb]
                                       * s_ia[ri, qi - nocb]) * RS)
            for pi in range(noca):
                for qi in range(nocb, nbf):
                    if (pi >= nocb) and (qi < noca):
                        continue
                    for ri in range(nocb, noca):
                        for si in range(nocb, noca):
                            acc += (co[pi, qi - nocb] * cn[ri, si - nocb]
                                    * (s_ij[pi, ri] * s_ab[qi - nocb, si - nocb]
                                       + s_ia[pi, si - nocb]
                                       * s_ia[ri, qi - nocb]) * RS)
            S[oi, ni] = acc
    for i in range(ns):
        S[:, i] /= np.linalg.norm(S[:, i])
    return S


def replica_S(ctx, M):
    Xt = ctx['Xt']
    return contraction(ctx, s_ij_of(ctx, M), s_ab_of(ctx, M),
                       s_ia_of(ctx, M), Xt, Xt)


def fortran_S(mol, ctx, theta, K):
    """Fortran state overlap at the MO rotation C -> C e^{theta K} (tlf=0)."""
    import oqp
    from scipy.linalg import expm
    W, Wb = ctx['W'], ctx['Wb']
    mol.data['OQP::VEC_MO_A_old'] = W
    mol.data['OQP::VEC_MO_B_old'] = Wb.copy()
    mol.data['OQP::E_MO_A_old'] = np.array(mol.data['OQP::E_MO_A'], copy=True)
    mol.data['OQP::E_MO_B_old'] = np.array(mol.data['OQP::E_MO_B'], copy=True)
    mol.data['OQP::td_bvec_mo_old'] = ctx['X0_raw'].copy()
    mol.data['OQP::xyz_old'] = np.array(mol.get_system(),
                                        copy=True).reshape((3, -1))
    mol.data.set_tdhf_tlf(0)
    mol.data['OQP::VEC_MO_A'] = expm(-theta * K) @ W
    mol.data['OQP::VEC_MO_B'] = expm(-theta * K) @ W
    mol.data['OQP::td_bvec_mo'] = ctx['X0_raw']
    oqp.get_structures_ao_overlap(mol)
    oqp.get_states_overlap(mol)
    S = np.array(mol.data['OQP::td_states_overlap'], copy=True)
    mol.data['OQP::VEC_MO_A'] = W
    mol.data['OQP::VEC_MO_B'] = Wb
    return S.T


# ---------------- closed-form gamma^formula --------------------------------
def _adjugate(A):
    n = A.shape[0]
    adj = np.zeros_like(A)
    for r in range(n):
        for c in range(n):
            Mm = np.delete(np.delete(A, r, axis=0), c, axis=1)
            adj[c, r] = ((-1.0) ** (r + c)) * la.det(Mm)
    return adj


def gamma_closed(ctx):
    """gamma[I,J,p,q] = dS_IJ/dtheta_pq of the exact formula, one pass.
    Convention: for the antisym generator K (K_pq=+1, K_qp=-1),
    dS/dtheta = 2 * gamma[..., p, q] (the derivative is split between the
    two antisymmetric slots)."""
    nstate, noca, nocb = ctx['nstate'], ctx['noca'], ctx['nocb']
    nvirb, nbf = ctx['nvirb'], ctx['nbf']
    Xt = ctx['Xt']
    I_mo = np.eye(nbf)

    minors_def = {}
    for i1 in range(1, noca + 1):
        for i2 in range(1, noca + 1):
            if i1 == i2:
                keep = [k for k in range(noca) if k != i1 - 1]
                minors_def[('ij', i1, i2)] = (keep, keep, 1.0)
            else:
                imin, imax = min(i1, i2), max(i1, i2)
                rows = ([k for k in range(noca) if k + 1 not in (imin, imax)]
                        + [i2 - 1])
                cols = ([k for k in range(noca) if k + 1 not in (imin, imax)]
                        + [i1 - 1])
                minors_def[('ij', i1, i2)] = (rows, cols, -1.0)
    core = list(range(nocb))
    for j1 in range(nvirb):
        for j2 in range(nvirb):
            minors_def[('ab', j1, j2)] = (core + [nocb + j1],
                                          core + [nocb + j2], 1.0)
    for i1 in range(1, noca + 1):
        for j1 in range(nvirb):
            rows, cols = _ia_maps(ctx, i1, j1)
            minors_def[('ia', i1 - 1, j1)] = (rows, cols, 1.0)

    W = {}
    for key, (rows, cols, sgn) in minors_def.items():
        A0 = I_mo[np.ix_(rows, cols)]
        adj0 = _adjugate(A0)
        Wm = {}
        n = len(rows)
        for a in range(n):
            for b in range(n):
                v = sgn * adj0[b, a]
                if v != 0.0:
                    Wm[(rows[a], cols[b])] = Wm.get((rows[a], cols[b]),
                                                    0.0) + v
        W[key] = Wm

    sij0 = s_ij_of(ctx, I_mo)
    sab0 = s_ab_of(ctx, I_mo)
    sia0 = s_ia_of(ctx, I_mo)
    eps = 1e-6
    gam = np.zeros((nstate, nstate, nbf, nbf))
    for key, Wm in W.items():
        kind = key[0]
        if kind == 'ij':
            i1, i2 = key[1] - 1, key[2] - 1
            s = sij0.copy(); s[i1, i2] += eps
            Sp = contraction(ctx, s, sab0, sia0, Xt, Xt)
            s[i1, i2] -= 2 * eps
            Sm = contraction(ctx, s, sab0, sia0, Xt, Xt)
        elif kind == 'ab':
            j1, j2 = key[1], key[2]
            s = sab0.copy(); s[j1, j2] += eps
            Sp = contraction(ctx, sij0, s, sia0, Xt, Xt)
            s[j1, j2] -= 2 * eps
            Sm = contraction(ctx, sij0, s, sia0, Xt, Xt)
        else:
            i1, j1 = key[1], key[2]
            s = sia0.copy(); s[i1, j1] += eps
            Sp = contraction(ctx, sij0, sab0, s, Xt, Xt)
            s[i1, j1] -= 2 * eps
            Sm = contraction(ctx, sij0, sab0, s, Xt, Xt)
        dSds = (Sp - Sm) / (2 * eps)
        for (p, q), w in Wm.items():
            gam[:, :, p, q] += 0.5 * dSds * w
            gam[:, :, q, p] -= 0.5 * dSds * w
    return gam


# ---------------- independent exact-overlap generator oracle --------------
def antisymmetric_generator(nbf, p, q):
    """Return K with K[p,q]=+1 and K[q,p]=-1 for one p>q coordinate."""
    if not (0 <= q < p < nbf):
        raise ValueError("an independent orbital generator requires 0 <= q < p < nbf")
    generator = np.zeros((nbf, nbf))
    generator[p, q] = 1.0
    generator[q, p] = -1.0
    return generator


def exact_generator_rotation(nbf, p, q, theta):
    """Return exp(theta*K_pq) as its exact two-orbital plane rotation."""
    antisymmetric_generator(nbf, p, q)  # shared bounds/order validation
    rotation = np.eye(nbf)
    cosine = np.cos(theta)
    sine = np.sin(theta)
    rotation[p, p] = cosine
    rotation[q, q] = cosine
    rotation[p, q] = sine
    rotation[q, p] = -sine
    return rotation


def exact_overlap_generator_derivative(ctx, p, q, step=1.0e-4):
    """Differentiate the exact state-overlap replica along one generator.

    A fourth-order centred Richardson stencil is applied to ``replica_S`` at
    ``exp(theta*K_pq)``.  The underlying overlaps are the literal exact-tlf
    determinant minors.  The returned array is dS[I,J]/dtheta for every
    ordered state pair; no state transpose or antisymmetry is imposed.
    """
    if step <= 0.0:
        raise ValueError("the generator derivative step must be positive")
    nbf = ctx['nbf']
    plus_1 = replica_S(ctx, exact_generator_rotation(nbf, p, q, step))
    minus_1 = replica_S(ctx, exact_generator_rotation(nbf, p, q, -step))
    plus_2 = replica_S(ctx, exact_generator_rotation(nbf, p, q, 2.0*step))
    minus_2 = replica_S(ctx, exact_generator_rotation(nbf, p, q, -2.0*step))
    return (8.0*(plus_1 - minus_1) - (plus_2 - minus_2))/(12.0*step)


def generator_derivative_sweep(ctx, step=1.0e-4, progress=None):
    """Evaluate dS/dtheta for every independent p>q orbital generator.

    The result has shape ``(nstate,nstate,nbf,nbf)``.  Only the independent
    lower-triangle slots ``[...,p,q]`` with p>q are populated.  In particular,
    this routine does not manufacture the upper triangle by orbital
    antisymmetrization and does not manufacture reverse state pairs.
    """
    nstate, nbf = ctx['nstate'], ctx['nbf']
    derivative = np.zeros((nstate, nstate, nbf, nbf))
    total = nbf*(nbf - 1)//2
    done = 0
    for q in range(nbf - 1):
        for p in range(q + 1, nbf):
            derivative[:, :, p, q] = exact_overlap_generator_derivative(
                ctx, p, q, step=step
            )
            done += 1
            if progress is not None:
                progress(done, total, p, q)
    return derivative
