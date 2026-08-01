"""v19: THE F-CHANNEL MEASUREMENT. One coordinate c0.
  w_skel : separate-process 1-iter workers (skel_gate.py --worker)
  w_ref  : in-process displaced full-SCF matvec (sign-fixed C)
  stagedC: staged response along the FIXED-Sk Ux(c0)
  gvec   : measured G-vector for dD_model(Ux(c0))
  trueF = w_ref - w_skel - stagedC;  diff = trueF - gvec
Run: python v19_fchan.py <inp> <coord> <skeldir>
"""
import os, sys, subprocess
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
    skeldir = sys.argv[3]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v19.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    ncoord = 3 * natom
    ctx = FK.build_context(mol)
    nbf, nij = ctx['nbf'], ctx['nij']
    noca, nocb = ctx['noca'], ctx['nocb']
    X0_raw = ctx['X0_raw']
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    C0 = W0.T
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

    J = 2
    XJ = Xf[:, J]
    base = matvec(XJ)
    SAVE0 = {k: np.array(mol.data[k], copy=True) for k in
             ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
              'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A', 'OQP::E_MO_B']}

    # workers with THIS process's phases
    ref_npz = os.path.join(skeldir, 'ref.npz')
    np.savez(ref_npz, C_a=W0, C_b=Wb0, X0_raw=X0_raw,
             nstate=nstate, nij=nij)
    helper = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'skel_gate.py')
    env = dict(os.environ, OMP_NUM_THREADS='4')
    procs = []
    for idx in (c0, c0 + ncoord):
        pinp = os.path.join(skeldir, f'p{idx}.inp')
        out = os.path.join(skeldir, f'p{idx}_v19.npz')
        procs.append(subprocess.Popen(
            [sys.executable, helper, '--worker', pinp, ref_npz, out],
            env=env, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT))
    # while workers run: in-process displaced full sweep for c0
    mol.save_data()
    cfg = mol.config
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = mol.log.replace('.log', '.json')
    cfg['guess']['continue_geom'] = False
    HD = 1e-3

    def disp(coord):
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
        S_np = np.array(mol.data['OQP::overlap_ao_non_orthogonal'], copy=True)
        S_f = S_np.reshape(-1).reshape((nbf, nbf)).T
        sg = np.sign(np.diag(M_f))
        sg[sg == 0] = 1.0
        M_f = M_f * sg[None, :]
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass
        Wd = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
        Wd_sf = sg[:, None] * Wd
        mol.data['OQP::VEC_MO_A'] = Wd_sf
        mol.data['OQP::VEC_MO_B'] = Wd_sf.copy()
        Ax = matvec(XJ)
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        Skf = C0.T @ S_f @ C0            # 7.49 fixed Sk
        return M_f, Skf, Ax

    cp = xyz0.copy(); cp[c0] += HD
    Mp, Skp, Axp = disp(cp)
    cm = xyz0.copy(); cm[c0] -= HD
    Mm, Skm, Axm = disp(cm)
    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass
    w_ref = (Axp - Axm) / (2 * HD)
    Ux = (Mp - Mm) / (2 * HD) - (Skp - Skm) / (2 * HD)

    # stagedC along fixed Ux
    t = 1e-5
    acc = np.zeros(nij)
    for sgn in (1.0, -1.0):
        Wp = (np.eye(nbf) + sgn * t * Ux.T) @ W0
        mol.data['OQP::VEC_MO_A'] = Wp
        mol.data['OQP::VEC_MO_B'] = Wp
        acc += sgn * matvec(XJ) / (2 * t)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    stagedC = acc

    # gvec for dD_model(Ux)
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
            continue
    mol.config['scf']['maxit'] = 1
    occ_a = np.zeros(nbf); occ_a[:noca] = 1.0
    occ_b = np.zeros(nbf); occ_b[:nocb] = 1.0
    Ma = Ux * occ_a[None, :] + Ux.T * occ_a[:, None]
    Mb = Ux * occ_b[None, :] + Ux.T * occ_b[:, None]
    dDa = C0 @ Ma @ C0.T
    dDb = C0 @ Mb @ C0.T
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
    epsF = 1e-4
    mol.data['OQP::FOCK_A'] = pack_sym(F0a + epsF * Ga).reshape(SAVE0['OQP::FOCK_A'].shape)
    mol.data['OQP::FOCK_B'] = pack_sym(F0b + epsF * Gb).reshape(SAVE0['OQP::FOCK_B'].shape)
    gf = matvec(XJ)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    gvec = (gf - base) / epsF

    # collect workers
    for p in procs:
        p.wait()
    Axp_w = np.load(os.path.join(skeldir, f'p{c0}_v19.npz'))['Ax'][J]
    Axm_w = np.load(os.path.join(skeldir, f'p{c0 + ncoord}_v19.npz'))['Ax'][J]
    w_skel = (Axp_w - Axm_w) / (2 * HD)

    trueF = w_ref - w_skel - stagedC
    diff = trueF - gvec
    print('===== v19 F-channel verdict =====')
    print(f'|w_ref|={np.linalg.norm(w_ref):.5f} |w_skel|={np.linalg.norm(w_skel):.5f} '
          f'|stagedC|={np.linalg.norm(stagedC):.5f}')
    print(f'|trueF|={np.linalg.norm(trueF):.5f} |gvec|={np.linalg.norm(gvec):.5f} '
          f'|trueF - gvec|={np.linalg.norm(diff):.5f} maxdiff={np.abs(diff).max():.4e}')
    ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
    ijlr2 = (noca - nocb - 1) * noca + noca
    idx = np.argsort(-np.abs(diff))[:10]
    for k in idx:
        i = k % noca + 1
        a = k // noca + 1
        tag = ' LR1' if k == ijlr1 - 1 else (' LR2' if k == ijlr2 - 1 else '')
        print(f'  slot {k}: (i={i}, a={a + nocb}) diff={diff[k]:+.5f} '
              f'trueF={trueF[k]:+.5f} gvec={gvec[k]:+.5f}{tag}')

if __name__ == '__main__':
    main()
