"""Step-2 gate: EXACT interstate one-particle transition density gamma^IJ
for MRSF SF determinants via Slater-Condon, with fermionic signs derived,
not scanned.

  gamma^{sigma}_pq = sum_mn Xt^I_m Xt^J_n <Phi_m| a+_{p sigma} a_{q sigma} |Phi_n>

Phi_(i,a) = ROHF triplet reference with one alpha electron removed from
spatial i (i in 1..noca) and one beta electron added to spatial a
(a in nocb+1..nbf). Xt = unfolded (determinant-grid) amplitudes.

Checks performed:
  A. exact gamma vs closed-form dgemm candidate (occ-occ alpha,
     virt-virt beta) -> localizes every socc/sign subtlety.
  B. antisymmetrized exact gamma_a vs the branch's sign-scanned TLF
     kernel (kernel_pair(I,J) - kernel_pair(J,I), cross-mask).
  C. (needs liboqp) rotation-FD of the TLF state overlap:
     d/dtheta <Psi_I(C) | Psi_J(C e^{theta K})> == sum_pq gamma_pq K_pq
     for random antisymmetric K per rotation block.

Run:  python gamma_gate.py H2O_energy.inp
"""
import sys
import numpy as np


# ---------------------------------------------------------------- Slater-Condon
def sc_tdm(dets_bra_coef, dets_ket_coef, dets, nbf):
    """Exact spin-resolved 1-TDM between two CI vectors over determinants.

    dets: list of (alpha_occ_tuple, beta_occ_tuple) with orbitals 0-based,
          tuples SORTED ascending (canonical order defines the sign).
    Returns (gamma_alpha, gamma_beta), gamma[p, q] = <I|a+_p a_q|J>.
    """
    ga = np.zeros((nbf, nbf))
    gb = np.zeros((nbf, nbf))

    def one_sector(occ_bra, occ_ket):
        """Return (kind, data) for a single spin sector.
        kind 'same' -> data None; kind 'single' -> data (p, q, sign);
        kind 'far' -> differ by 2+."""
        sb, sk = set(occ_bra), set(occ_ket)
        db = sorted(sb - sk)          # in bra only  (create p)
        dk = sorted(sk - sb)          # in ket only  (annihilate q)
        if not db and not dk:
            return 'same', None
        if len(db) == 1 and len(dk) == 1:
            p, q = db[0], dk[0]
            # sign: move q to the end of ket list, p to end of bra list
            iq = occ_ket.index(q)
            ip = occ_bra.index(p)
            sign = (-1) ** (len(occ_ket) - 1 - iq) * (-1) ** (len(occ_bra) - 1 - ip)
            return 'single', (p, q, sign)
        return 'far', None

    n = len(dets)
    for m in range(n):
        am, bm = dets[m]
        for k in range(n):
            an, bn = dets[k]
            w = dets_bra_coef[m] * dets_ket_coef[k]
            if w == 0.0:
                continue
            ka, da_ = one_sector(am, an)
            kb, db_ = one_sector(bm, bn)
            if ka == 'far' or kb == 'far':
                continue
            if ka == 'same' and kb == 'same':      # diagonal: occupation numbers
                for p in am:
                    ga[p, p] += w
                for p in bm:
                    gb[p, p] += w
            elif ka == 'single' and kb == 'same':
                p, q, s = da_
                ga[p, q] += s * w
            elif ka == 'same' and kb == 'single':
                p, q, s = db_
                gb[p, q] += s * w
            # single x single -> two-particle, not in 1-TDM


# ------------------------------------------------------------------- main gate
def main():
    import oqp                        # before numpy-heavy imports (ILP64)
    from oqp.pyoqp import Runner

    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_gg.log'))
    r.run()
    mol = r.mol

    nstate = mol.config['tdhf']['nstate']
    noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
    nocb = noca - 2
    C = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    nbf = C.shape[0]
    nvirb = nbf - nocb
    nij = noca * nvirb
    RS = 1.0 / np.sqrt(2.0)
    X0 = np.array(mol.data['OQP::td_bvec_mo'], copy=True
                  ).reshape(-1).reshape((nstate, nij)).T.copy()

    # ---- unfold (constant V; the mult=1 convention, validated in Phase 11)
    def unfold(bv, st):
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

    Xt = [unfold(X0, s + 1) for s in range(nstate)]

    # ---- determinant list: (alpha_occ, beta_occ), 0-based, sorted
    ref_a = list(range(noca))
    ref_b = list(range(nocb))
    dets, amp_index = [], {}
    for i in range(noca):
        for a in range(nvirb):
            aocc = tuple(sorted(set(ref_a) - {i}))
            bocc = tuple(sorted(set(ref_b) | {nocb + a}))
            amp_index[(i, a)] = len(dets)
            dets.append((aocc, bocc))

    def coefs(x):
        # Code amplitudes multiply |i->a> = a+_{a b} a_{i a} |ref>; converting
        # to the sorted-occupancy determinant kets used here costs the alpha
        # annihilation parity (-1)^(noca-1-i) (electrons after slot i in the
        # alpha string). The beta creation parity is det-independent (a lands
        # at the end of the sorted beta string) and cancels in all bilinears.
        c = np.zeros(len(dets))
        for (i, a), idx in amp_index.items():
            c[idx] = x[i, a] * ((-1.0) ** (noca - 1 - i))
        return c

    print(f'\nnbf={nbf} noca={noca} nocb={nocb} nvirb={nvirb} ndet={len(dets)}')

    # ---- exact TDMs for every pair
    def tdm_pair(I, J):
        ga = np.zeros((nbf, nbf))
        gb = np.zeros((nbf, nbf))
        cb, ck = coefs(Xt[I]), coefs(Xt[J])

        def one_sector(occ_bra, occ_ket):
            sb, sk = set(occ_bra), set(occ_ket)
            db = sorted(sb - sk)
            dk = sorted(sk - sb)
            if not db and not dk:
                return 'same', None
            if len(db) == 1 and len(dk) == 1:
                p, q = db[0], dk[0]
                iq = occ_ket.index(q)
                ip = occ_bra.index(p)
                sign = (-1) ** (len(occ_ket) - 1 - iq) * (-1) ** (len(occ_bra) - 1 - ip)
                return 'single', (p, q, sign)
            return 'far', None

        for m in range(len(dets)):
            if cb[m] == 0.0:
                continue
            am, bm = dets[m]
            for k in range(len(dets)):
                w = cb[m] * ck[k]
                if w == 0.0:
                    continue
                an, bn = dets[k]
                ka, da_ = one_sector(am, an)
                if ka == 'far':
                    continue
                kb, db_ = one_sector(bm, bn)
                if kb == 'far':
                    continue
                if ka == 'same' and kb == 'same':
                    for p in am:
                        ga[p, p] += w
                    for p in bm:
                        gb[p, p] += w
                elif ka == 'single' and kb == 'same':
                    p, q, s = da_
                    ga[p, q] += s * w
                elif ka == 'same' and kb == 'single':
                    p, q, s = db_
                    gb[p, q] += s * w
        return ga, gb

    # ---- closed-form candidate (to be CONFIRMED against exact, incl. signs)
    # alpha occ-occ: gamma_a[j,i] = -sum_a Xt^I[i,a] Xt^J[j,a] (hole density)
    # beta virt-virt: gamma_b[a,b] = +sum_i Xt^I[i,a] Xt^J[i,b] (particle)
    # diagonal reference terms included separately.
    def closed_form(I, J):
        ga = np.zeros((nbf, nbf))
        gb = np.zeros((nbf, nbf))
        xi, xj = Xt[I], Xt[J]
        ov = float(np.sum(xi * xj))
        # diagonal occupation part
        for p in range(noca):
            ga[p, p] += ov
        for p in range(nocb):
            gb[p, p] += ov
        ga[:noca, :noca] -= xj @ xi.T          # gamma_a[j,i] -= sum_a xj[j,a] xi[i,a]
        Gv = xi.T @ xj                          # [a,b] = sum_i xi[i,a] xj[i,b]
        gb[nocb:, nocb:] += Gv
        return ga, gb

    print('\n===== A. exact Slater-Condon vs closed-form candidate =====')
    for I in range(nstate):
        for J in range(nstate):
            if I >= J:
                continue
            ga_e, gb_e = tdm_pair(I, J)
            ga_c, gb_c = closed_form(I, J)
            print(f'pair ({I+1},{J+1}): '
                  f'max|ga_exact - ga_closed| = {np.abs(ga_e-ga_c).max():.3e}   '
                  f'max|gb_exact - gb_closed| = {np.abs(gb_e-gb_c).max():.3e}')

    print('\n===== B. antisymmetrized spatial gamma vs branch TLF kernel =====')
    # branch kernel (sign-scanned) rebuilt via the NAC class helper
    from oqp.library.single_point import NAC
    nac = NAC(mol)
    nac._build_nac_gamma_tlf()
    gt = np.array(mol.data['OQP::nac_gamma_tlf'], copy=True
                  ).reshape(-1)
    gt = gt.reshape((nstate * nstate, nbf * nbf))
    for I in range(nstate):
        for J in range(nstate):
            if I >= J:
                continue
            ga_e, gb_e = tdm_pair(I, J)
            g_sp = ga_e + gb_e                 # spatial total (contracts T_pq)
            g_a = 0.5 * (g_sp - g_sp.T)        # only antisym part survives
            # branch buffer: F-flat index (I + J*nstate), C-flat of G^T
            k = np.linalg.norm
            gtlf = gt[I + J * nstate].reshape((nbf, nbf)).T
            num = float(np.sum(g_a * gtlf))
            cos = num / (k(g_a) * k(gtlf) + 1e-300)
            print(f'pair ({I+1},{J+1}): |gamma_a|={k(g_a):.6f} |gamma_TLF|={k(gtlf):.6f} '
                  f'cos={cos:+.6f} ratio={k(gtlf)/(k(g_a)+1e-300):.6f}')
            # per-block comparison (doc/socc/virt)
            sl = {'d': slice(0, nocb), 's': slice(nocb, noca), 'v': slice(noca, nbf)}
            for bn1 in 'dsv':
                for bn2 in 'dsv':
                    b1, b2 = sl[bn1], sl[bn2]
                    na_, nt_ = k(g_a[b1, b2]), k(gtlf[b1, b2])
                    if na_ < 1e-12 and nt_ < 1e-12:
                        continue
                    cb_ = float(np.sum(g_a[b1, b2] * gtlf[b1, b2])) / (na_ * nt_ + 1e-300)
                    print(f'    block {bn1}{bn2}: |exact|={na_:.5f} |tlf|={nt_:.5f} cos={cb_:+.4f}')

    # ---- C runs twice: production tlf=2 and exact-minor tlf=0 ----------
    try:
        mol.data.set_tdhf_tlf(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
        print(f'\n(tlf set to {int(sys.argv[2]) if len(sys.argv) > 2 else 0} for gate C)')
    except Exception as e:
        print('(could not set tlf:', e, ')')
    print('\n===== C. rotation-FD of the TLF overlap vs gamma (per block) =====')
    OVTAG = 'OQP::overlap_mo_non_orthogonal'
    ao_S = None
    try:
        import oqp as _oqp
        mol.data['OQP::VEC_MO_A_old'] = C
        mol.data['OQP::VEC_MO_B_old'] = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
        mol.data['OQP::E_MO_A_old'] = np.array(mol.data['OQP::E_MO_A'], copy=True)
        mol.data['OQP::E_MO_B_old'] = np.array(mol.data['OQP::E_MO_B'], copy=True)
        mol.data['OQP::td_bvec_mo_old'] = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
        mol.data['OQP::xyz_old'] = np.array(mol.get_system(), copy=True).reshape((3, -1))
    except Exception as e:
        print('  (skipped: cannot stage old-geometry records:', e, ')')
        return

    rng = np.random.default_rng(7)
    th = 1e-5
    Cb = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    X_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)

    def overlap_matrix():
        import oqp as _o
        _o.get_structures_ao_overlap(mol)
        _o.get_states_overlap(mol)
        return np.array(mol.data['OQP::td_states_overlap'], copy=True)

    blocks = {'ds': (slice(0, nocb), slice(nocb, noca)),
              'dv': (slice(0, nocb), slice(noca, nbf)),
              'sv': (slice(nocb, noca), slice(noca, nbf)),
              'ss': (slice(nocb, noca), slice(nocb, noca))}
    gam_sp = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            ga_e, gb_e = tdm_pair(min(I, J), max(I, J))
            g = ga_e + gb_e
            if I > J:
                g = g.T
            gam_sp[(I, J)] = g

    for bname, (lo, hi) in blocks.items():
        K = np.zeros((nbf, nbf))
        blk = rng.standard_normal((hi.stop - hi.start, lo.stop - lo.start))
        K[hi, lo] = blk
        K[lo, hi] = -blk.T
        # rotate ONLY the current MOs; "old" stays at reference.
        # NOTE: numpy tagarray 2-D = TRANSPOSE of the Fortran matrix, so the
        # MO rotation C_f -> C_f e^{th K} is staged as W -> expm(-th K) @ W.
        from scipy.linalg import expm
        mol.data['OQP::VEC_MO_A'] = expm(-th * K) @ C
        mol.data['OQP::VEC_MO_B'] = expm(-th * K) @ C
        mol.data['OQP::td_bvec_mo'] = X_raw
        Sp = overlap_matrix()
        mol.data['OQP::VEC_MO_A'] = expm(th * K) @ C
        mol.data['OQP::VEC_MO_B'] = expm(th * K) @ C
        mol.data['OQP::td_bvec_mo'] = X_raw
        Sm = overlap_matrix()
        dS = (Sp - Sm) / (2 * th)
        print(f'  block {bname}:')
        for I in range(nstate):
            for J in range(nstate):
                if I >= J:
                    continue
                an = float(np.sum(gam_sp[(I, J)] * K))
                # branch TLF kernel (formula-representation candidate):
                # buffer F-layout [I + J*nstate], C-flat of G^T
                gtlf = gt[I + J * nstate].reshape((nbf, nbf)).T
                at = float(np.sum(gtlf * K))
                print(f'    ({I+1},{J+1}): FD={dS[I, J]:+.8f}  raw_gamma.K={an:+.8f} '
                      f'(diff {abs(dS[I, J]-an):.1e})  gammaTLF.K={at:+.8f} '
                      f'(diff {abs(dS[I, J]-at):.1e})')
    # restore
    mol.data['OQP::VEC_MO_A'] = C
    mol.data['OQP::VEC_MO_B'] = Cb
    mol.data['OQP::td_bvec_mo'] = X_raw

    print('\n===== D. EXACT biorthogonal overlap FD vs gamma (code-independent) =====')
    # <Psi_I(C) | Psi_J(C e^{th K})> = sum_mn c^I_m c^J_n det(M[am,an]) det(M[bm,bn]),
    # M = e^{th K} (MO cross-overlap at the same geometry). First order must
    # equal sum_pq gamma_pq K_pq EXACTLY -- this is the self-consistency gate
    # for the Slater-Condon machinery, independent of the Fortran TLF code.
    from scipy.linalg import expm as _expm
    cvecs = [coefs(Xt[s]) for s in range(nstate)]

    def exact_overlap(M):
        S = np.zeros((nstate, nstate))
        for m in range(len(dets)):
            am, bm = dets[m]
            wa = np.array([cv[m] for cv in cvecs])
            if np.all(wa == 0.0):
                continue
            for k in range(len(dets)):
                an, bn = dets[k]
                wb = np.array([cv[k] for cv in cvecs])
                if np.all(wb == 0.0):
                    continue
                ov = (np.linalg.det(M[np.ix_(am, an)])
                      * np.linalg.det(M[np.ix_(bm, bn)]))
                S += np.outer(wa, wb) * ov
        return S

    for bname, (lo, hi) in blocks.items():
        K = np.zeros((nbf, nbf))
        rngD = np.random.default_rng(11)
        blk = rngD.standard_normal((hi.stop - hi.start, lo.stop - lo.start))
        K[hi, lo] = blk
        K[lo, hi] = -blk.T
        Sp = exact_overlap(_expm(th * K))
        Sm = exact_overlap(_expm(-th * K))
        dS = (Sp - Sm) / (2 * th)
        print(f'  block {bname}:')
        for I in range(nstate):
            for J in range(nstate):
                if I >= J:
                    continue
                an_ = float(np.sum(gam_sp[(I, J)] * K))
                print(f'    ({I+1},{J+1}): exactFD={dS[I, J]:+.8f}  gamma.K={an_:+.8f}  '
                      f'diff={abs(dS[I, J]-an_):.2e}')


if __name__ == '__main__':
    main()
