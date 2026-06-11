"""Offline probe: which transition-density convention closes the NAC gap?

Uses theta_data_{hf,bhh}.npz from theta_check.py.  Reconstructs gamma^IJ from
the raw folded amplitudes with adjustable open-shell scalings and tests
    d = d_CI + sum_pq gamma_pq theta_pq   vs   d_num
for each convention, on both functionals.
"""
import numpy as np

RS = 1.0 / np.sqrt(2.0)


def load(tag):
    d = np.load(f'/tmp/nactest/theta_data_{tag}.npz')
    return d


def unfold(bv, st, noca, nocb, nbf):
    """xv12[i, a] (i=1..noca occ, a=1..nvirb virt-of-beta), folded -> unfolded.
    st is 1-based state index. Follows get_mrsf_transition_density (mult=1)."""
    nvirb = nbf - nocb
    ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
    ijlr2 = (noca - nocb - 1) * noca + noca
    x = np.zeros((noca, nvirb))
    for i in range(1, noca + 1):
        for jj in range(nocb + 1, nbf + 1):
            ij = (jj - nocb - 1) * noca + i
            if ij == ijlr1:
                x[i - 1, jj - nocb - 1] = bv[ijlr1 - 1, st - 1] * RS
            elif ij == ijlr2:
                x[i - 1, jj - nocb - 1] = -bv[ijlr1 - 1, st - 1] * RS
            else:
                x[i - 1, jj - nocb - 1] = bv[ij - 1, st - 1]
    return x


def trans_den(xi, xj, noca, nocb, nbf, scal_oo=np.sqrt(2.0), scal_vv=np.sqrt(2.0)):
    """gamma[p,q] following get_trans_den, with adjustable open-shell scalings.
    xi = amplitudes of state I (bra), xj = of state J (ket)."""
    nvirb = nbf - nocb
    g = np.zeros((nbf, nbf))
    # vv block: trden(i+nocb, j+nocb) += sum_k scal * xj[k,i] * xi[k,j]
    for i in range(nvirb):
        for j in range(nvirb):
            t = 0.0
            for k in range(noca):
                ioo = (i in (0, 1)) and (k in (noca - 2, noca - 1))
                joo = (j in (0, 1)) and (k in (noca - 2, noca - 1))
                s = scal_vv if (ioo != joo) else 1.0
                t += s * xj[k, i] * xi[k, j]
            g[i + nocb, j + nocb] += t
    # oo block: trden(i,j) -= sum_k scal * xj[j,k] * xi[i,k]
    for i in range(noca):
        for j in range(noca):
            t = 0.0
            for k in range(nvirb):
                ioo = (i in (noca - 2, noca - 1)) and (k in (0, 1))
                joo = (j in (noca - 2, noca - 1)) and (k in (0, 1))
                s = scal_oo if (ioo != joo) else 1.0
                t += s * xj[j, k] * xi[i, k]
            g[i, j] -= t
    return g


def assemble(d, gam):
    theta = d['theta']                  # (ncoord, nbf, nbf)
    dci = d['dci']
    return dci + np.einsum('pq,kpq->k', gam, theta)


def cmp(name, v, ref):
    c = np.dot(ref, v) / (np.linalg.norm(ref) * np.linalg.norm(v) + 1e-30)
    print(f'  {name}: cos={c:+.6f}  |v|={np.linalg.norm(v):.5f}  '
          f'|ref|={np.linalg.norm(ref):.5f}  ratio={np.linalg.norm(v)/np.linalg.norm(ref):.4f}')


for tag in ('hf', 'bhh'):
    d = load(tag)
    noca, nocb, nbf = int(d['noca']), int(d['nocb']), int(d['nbf'])
    bvraw = d['bvec']
    nij = noca * (nbf - nocb)
    bv = bvraw.reshape(-1)
    nst = bv.size // nij
    # figure out the layout: Fortran (nij, nstate) -> python flat in F order
    bvF = bv.reshape((nst, nij)).T      # bvF[ij-1, st-1]
    i, j = 2, 3
    xi = unfold(bvF, i, noca, nocb, nbf)
    xj = unfold(bvF, j, noca, nocb, nbf)
    dn = d['dn_al']

    # baseline: reproduce the Fortran gamma exactly?
    tdF = d['tdF']
    gam_code = tdF[j - 1, i - 1].reshape((nbf, nbf)).T
    gam_re = trans_den(xi, xj, noca, nocb, nbf)
    print(f'\n===== {tag}: pair (2,3) =====')
    print('  reconstruction check |gam_re - gam_code| =',
          np.abs(gam_re - gam_code).max())

    # antisymmetry diagnostic: d_orb(I,J) + d_orb(J,I) should be ~0
    gam_ji = trans_den(xj, xi, noca, nocb, nbf)
    aIJ = np.einsum('pq,kpq->k', gam_re, d['theta'])
    aJI = np.einsum('pq,kpq->k', gam_ji, d['theta'])
    print('  |d_orb(I,J)+d_orb(J,I)| =', np.linalg.norm(aIJ + aJI),
          '  vs |d_orb| =', np.linalg.norm(aIJ))

    cmp('code gamma        ', assemble(d, gam_re), dn)
    cmp('no sqrt2 (scal=1) ', assemble(d, trans_den(xi, xj, noca, nocb, nbf, 1.0, 1.0)), dn)
    cmp('scal=2            ', assemble(d, trans_den(xi, xj, noca, nocb, nbf, 2.0, 2.0)), dn)
    cmp('scal=1/sqrt2      ', assemble(d, trans_den(xi, xj, noca, nocb, nbf, RS, RS)), dn)
    cmp('antisym part only ', assemble(d, 0.5 * (gam_re - gam_re.T) * 2)[:], dn)
    # least-squares: best uniform factor f on gamma: dci + f*sum(gam theta) = dn
    base = np.einsum('pq,kpq->k', gam_re, d['theta'])
    f = np.dot(dn - d['dci'], base) / np.dot(base, base)
    print(f'  best uniform gamma factor f = {f:+.6f}')
    cmp(f'f*gamma           ', d['dci'] + f * base, dn)
