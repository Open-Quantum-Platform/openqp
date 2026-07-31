"""ABSOLUTE orientation gate: does the Python NAC pipeline's dc[i,j] hold
d_ij = <I|dJ> or d_ji = <J|dI>?

Method: one real geometry displacement. Stage (R0 -> R0+h) exactly like a
nacme worker, then compute the state overlap TWO ways:
  (a) production Fortran get_states_overlap, read the tag the way
      NACME.nacme() does (NO transpose) -> dc_python
  (b) code-independent EXACT biorthogonal overlap over the 90 SF
      determinants with operator phases (validated machinery, gate D) using
      the same MO cross-overlap and the same amplitude pair -> F_exact with
      UNAMBIGUOUS orientation F_exact[o,n] = <O(R0)|N(R0+h)>.
Then d_true[i,j] = (F_exact[i,j] - F_exact[j,i])/(2h) ~= <I|dJ> and the
verdict is sign(dc_python[i,j]) vs sign(d_true[i,j]).

Run:  python orientation_gate.py H2O_energy.inp
"""
import sys
import numpy as np


def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint

    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_og.log'))
    r.run()
    mol = r.mol

    nstate = mol.config['tdhf']['nstate']
    noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
    nocb = noca - 2
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)     # numpy = C_f^T
    nbf = W0.shape[0]
    nvirb = nbf - nocb
    nij = noca * nvirb
    RS = 1.0 / np.sqrt(2.0)

    X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    E0 = list(mol.energies)

    # displace along H1-x (large d(2,3) component in the frozen reference)
    h = 1e-3
    coord = xyz0.copy()
    coord[3] += h

    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False

    mol.update_system(coord)
    oqp.library.ints_1e(mol)
    oqp.library.guess(mol)
    SinglePoint(mol).energy()

    # stage old = R0
    mol.data['OQP::xyz_old'] = xyz0.reshape((3, -1))
    mol.data['OQP::VEC_MO_A_old'] = W0
    mol.data['OQP::VEC_MO_B_old'] = Wb0
    mol.data['OQP::E_MO_A_old'] = e0a
    mol.data['OQP::E_MO_B_old'] = e0b
    mol.data['OQP::td_bvec_mo_old'] = X0_raw
    mol.data.set_tdhf_tlf(0)

    oqp.get_structures_ao_overlap(mol)
    oqp.get_states_overlap(mol)

    # (a) production-python read (exactly like NACME.nacme)
    S_py = np.array(mol.data['OQP::td_states_overlap'], copy=True)
    dc_py = (S_py - S_py.T) / (2 * h)

    # (b) exact biorthogonal overlap with unambiguous orientation
    M_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
    M_f = M_np.reshape(-1).reshape((nbf, nbf)).T   # Fortran matrix: rows=old MO

    Xd_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
    X0m = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    Xdm = Xd_raw.reshape(-1).reshape((nstate, nij)).T.copy()

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

    Xt0 = [unfold(X0m, s + 1) for s in range(nstate)]
    Xtd = [unfold(Xdm, s + 1) for s in range(nstate)]

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
        c = np.zeros(len(dets))
        for (i, a), idx in amp_index.items():
            c[idx] = x[i, a] * ((-1.0) ** (noca - 1 - i))
        return c

    cvecs_o = [coefs(Xt0[s]) for s in range(nstate)]
    cvecs_n = [coefs(Xtd[s]) for s in range(nstate)]

    F = np.zeros((nstate, nstate))
    for m, (am, bm) in enumerate(dets):
        wa = np.array([cv[m] for cv in cvecs_o])
        if np.all(wa == 0.0):
            continue
        for k, (an, bn) in enumerate(dets):
            wb = np.array([cv[k] for cv in cvecs_n])
            if np.all(wb == 0.0):
                continue
            ov = (np.linalg.det(M_f[np.ix_(am, an)])
                  * np.linalg.det(M_f[np.ix_(bm, bn)]))
            F += np.outer(wa, wb) * ov

    # phase-align the displaced states (columns) to R0: flip any column whose
    # diagonal overlap is negative (random per-run Davidson phases otherwise
    # corrupt the one-sided FD); apply the SAME alignment to both oracles.
    colsign = np.sign(np.diag(F))
    colsign[colsign == 0] = 1.0
    F = F * colsign[None, :]
    S_py_al = S_py * colsign[None, :] if np.all(np.diag(S_py * colsign[None, :])
                                                > 0) else S_py
    dc_py = (S_py_al - S_py_al.T) / (2 * h)
    d_true = (F - F.T) / (2 * h)

    np.set_printoptions(precision=6, suppress=True)
    print('S_python (production read, no transpose):')
    print(S_py)
    print('F_exact  (rows = R0 states, cols = displaced states):')
    print(F)
    print('\ndc_python = (S_py - S_py^T)/(2h):')
    print(dc_py)
    print('d_true[i,j] = <I|dJ> (exact, unambiguous orientation):')
    print(d_true)
    print('\nVERDICT per pair (sign comparison on the dominant component):')
    for i in range(nstate):
        for j in range(nstate):
            if i >= j:
                continue
            s_py = dc_py[i, j]
            s_tr = d_true[i, j]
            rel = s_py / s_tr if abs(s_tr) > 1e-10 else float('nan')
            tag = 'd_ij (canonical)' if rel > 0 else 'd_ji (TRANSPOSED!)'
            print(f'  ({i+1},{j+1}): dc_py={s_py:+.6f}  d_true={s_tr:+.6f}  '
                  f'ratio={rel:+.3f}  -> dc_python holds {tag}')


if __name__ == '__main__':
    main()
