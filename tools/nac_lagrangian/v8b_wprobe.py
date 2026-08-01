"""v8b: FIXED-REFEREE probe: push sign-fixed C(x')diag(sg) into VEC_MO
(instead of sg-transporting vectors) -- the fold-consistent gauge.
  lhs = [w_ref - w_skel](c0)          (displaced FD, in-process)
  rhs = staged response along Ux(c0)  (A(C(1+tU))X - A(C)X)/t
      + G-channel vector: the Fock part of A applied with dF = G[dD]
Compare per raw amplitude slot; unfold the worst slots to (i,a) labels.
Run: python v8_wprobe.py <inp> <coord-index>
"""
import os, sys
import numpy as np

def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK
    inp = sys.argv[1]
    c0 = int(sys.argv[2])
    os.environ['NAC_DUMP_DS'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v8b.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    ctx = FK.build_context(mol)
    nbf, nij = ctx['nbf'], ctx['nij']
    noca, nocb, nvirb, RS = ctx['noca'], ctx['nocb'], ctx['nvirb'], ctx['RS']
    X0_raw = ctx['X0_raw']
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    C = W0.T
    ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
    ijlr2 = (noca - nocb - 1) * noca + noca

    def unfold_vec(v):
        x = np.zeros((noca, nvirb))
        for i in range(1, noca + 1):
            for jj in range(nocb + 1, nbf + 1):
                ij = (jj - nocb - 1) * noca + i
                if ij == ijlr1:
                    x[i - 1, jj - nocb - 1] = v[ijlr1 - 1] * RS
                elif ij == ijlr2:
                    x[i - 1, jj - nocb - 1] = -v[ijlr1 - 1] * RS
                else:
                    x[i - 1, jj - nocb - 1] = v[ij - 1]
        return x

    def refold(x):
        v = np.zeros(nij)
        for i in range(1, noca + 1):
            for jj in range(nocb + 1, nbf + 1):
                ij = (jj - nocb - 1) * noca + i
                if ij == ijlr1:
                    v[ijlr1 - 1] = x[i - 1, jj - nocb - 1] / RS
                elif ij == ijlr2:
                    pass
                else:
                    v[ij - 1] = x[i - 1, jj - nocb - 1]
        return v

    def sg_apply(v, sg):
        x = unfold_vec(v)
        x = sg[:noca, None] * x * sg[None, nocb:]
        return refold(x)

    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    def matvec(v):
        rr = X0_raw.copy().reshape(-1)
        rr[0:nij] = v
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        oqp.mrsf_matvec_apply(mol)
        return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()[:nij].copy()

    J = 2                      # state 3 (the near-degenerate pair's ket)
    XJ = Xf[:, J]
    base = matvec(XJ)

    # displaced full and skel w at coordinate c0
    SAVE0 = {k: np.array(mol.data[k], copy=True) for k in
             ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
              'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A', 'OQP::E_MO_B']}
    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False
    HD = 1e-3

    def disp(coord, frozen):
        mol.update_system(coord)
        oqp.library.ints_1e(mol)
        oqp.library.guess(mol)
        SinglePoint(mol).energy()
        mol.data['OQP::xyz_old'] = xyz0.reshape((3, -1))
        mol.data['OQP::VEC_MO_A_old'] = W0
        mol.data['OQP::VEC_MO_B_old'] = Wb0
        mol.data['OQP::E_MO_A_old'] = e0a
        mol.data['OQP::E_MO_B_old'] = e0b
        oqp.get_structures_ao_overlap(mol)
        M_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        M_f = M_np.reshape(-1).reshape((nbf, nbf)).T
        sg = np.sign(np.diag(M_f))
        sg[sg == 0] = 1.0
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass
        if frozen:
            mol.data['OQP::VEC_MO_A'] = W0
            mol.data['OQP::VEC_MO_B'] = Wb0
            return matvec(XJ)
        else:
            # sign-fix the displaced MOs themselves: C' <- C' diag(sg)
            Wd = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
            Wd = sg[:, None] * Wd
            mol.data['OQP::VEC_MO_A'] = Wd
            mol.data['OQP::VEC_MO_B'] = Wd.copy()
            return matvec(XJ)

    cp = xyz0.copy(); cp[c0] += HD
    cm = xyz0.copy(); cm[c0] -= HD
    wr = (disp(cp, False) - disp(cm, False)) / (2 * HD)
    # restore between the two modes
    mol.update_system(xyz0); oqp.library.ints_1e(mol)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    ws = (disp(cp, True) - disp(cm, True)) / (2 * HD)
    mol.update_system(xyz0); oqp.library.ints_1e(mol)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw

    lhs = wr - ws

    # staged response along Ux(c0): need Ux(c0): from the same disp machinery
    def dispM(coord):
        mol.update_system(coord)
        oqp.library.ints_1e(mol)
        oqp.library.guess(mol)
        SinglePoint(mol).energy()
        mol.data['OQP::xyz_old'] = xyz0.reshape((3, -1))
        mol.data['OQP::VEC_MO_A_old'] = W0
        mol.data['OQP::VEC_MO_B_old'] = Wb0
        mol.data['OQP::E_MO_A_old'] = e0a
        mol.data['OQP::E_MO_B_old'] = e0b
        oqp.get_structures_ao_overlap(mol)
        M_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        M_f = M_np.reshape(-1).reshape((nbf, nbf)).T
        sg = np.sign(np.diag(M_f))
        sg[sg == 0] = 1.0
        M_f = M_f * sg[None, :]
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        oqp.get_structures_ao_overlap(mol)
        Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
        return M_f, Skf

    Mp, Skp = dispM(cp)
    Mm, Skm = dispM(cm)
    Ux = (Mp - Mm) / (2 * HD) - (Skp - Skm) / (2 * HD)
    mol.update_system(xyz0); oqp.library.ints_1e(mol)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    t = 1e-5
    Wp1 = (np.eye(nbf) + t * Ux.T) @ W0
    mol.data['OQP::VEC_MO_A'] = Wp1
    mol.data['OQP::VEC_MO_B'] = Wp1
    a_p = matvec(XJ)
    Wm1 = (np.eye(nbf) - t * Ux.T) @ W0
    mol.data['OQP::VEC_MO_A'] = Wm1
    mol.data['OQP::VEC_MO_B'] = Wm1
    a_m = matvec(XJ)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    staged = (a_p - a_m) / (2 * t)

    # G-channel vector: replace FOCK by FOCK + eps*G[dD(Ux)], frozen C
    occ_a = np.zeros(nbf); occ_a[:noca] = 1.0
    occ_b = np.zeros(nbf); occ_b[:nocb] = 1.0
    Ma = Ux * occ_a[None, :] + Ux.T * occ_a[:, None]
    Mb = Ux * occ_b[None, :] + Ux.T * occ_b[:, None]
    dDa = C @ Ma @ C.T
    dDb = C @ Mb @ C.T
    nbf_tri = nbf * (nbf + 1) // 2

    def unpack_sym(pk):
        M = np.zeros((nbf, nbf))
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                M[p, q] = pk[idx]; M[q, p] = pk[idx]; idx += 1
        return M

    def pack_sym(M):
        pk = np.zeros(nbf_tri)
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                pk[idx] = M[p, q]; idx += 1
        return pk

    D0a = unpack_sym(SAVE0['OQP::DM_A'].ravel())
    D0b = unpack_sym(SAVE0['OQP::DM_B'].ravel())
    F0a = unpack_sym(SAVE0['OQP::FOCK_A'].ravel())
    F0b = unpack_sym(SAVE0['OQP::FOCK_B'].ravel())
    for attr in ('maxit', 'scf_maxit'):
        try:
            setattr(mol.data._data.control, attr, 1)
            break
        except Exception:
            pass
    mol.config['scf']['maxit'] = 1
    EPSG = 1e-4
    mol.data['OQP::DM_A'] = pack_sym(D0a + EPSG * dDa).reshape(SAVE0['OQP::DM_A'].shape)
    mol.data['OQP::DM_B'] = pack_sym(D0b + EPSG * dDb).reshape(SAVE0['OQP::DM_B'].shape)
    oqp.hf_energy(mol)
    Fa = unpack_sym(np.array(mol.data['OQP::FOCK_A'], copy=True).ravel())
    Fb = unpack_sym(np.array(mol.data['OQP::FOCK_B'], copy=True).ravel())
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    Ga = (Fa - F0a) / EPSG
    Gb = (Fb - F0b) / EPSG
    # perturb FOCK records by eps*G and matvec (frozen C) -> G-channel vector
    epsF = 1e-4
    mol.data['OQP::FOCK_A'] = pack_sym(F0a + epsF * Ga).reshape(SAVE0['OQP::FOCK_A'].shape)
    mol.data['OQP::FOCK_B'] = pack_sym(F0b + epsF * Gb).reshape(SAVE0['OQP::FOCK_B'].shape)
    gf = matvec(XJ)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    gvec = (gf - base) / epsF

    rhs = staged + gvec
    diff = lhs - rhs
    print(f'c0={c0} J={J+1}: |lhs|={np.linalg.norm(lhs):.5f} |rhs|={np.linalg.norm(rhs):.5f} '
          f'|diff|={np.linalg.norm(diff):.5f} maxdiff={np.abs(diff).max():.4e}')
    idx = np.argsort(-np.abs(diff))[:12]
    print('worst raw slots (i,a 1-based; a offset by nocb):')
    for k in idx:
        i = k % noca + 1
        a = k // noca + 1
        tag = ' LR1' if k == ijlr1 - 1 else (' LR2' if k == ijlr2 - 1 else '')
        print(f'  slot {k}: (i={i}, a={a + nocb}) diff={diff[k]:+.5f} '
              f'lhs={lhs[k]:+.5f} rhs={rhs[k]:+.5f}{tag}')

if __name__ == '__main__':
    main()
