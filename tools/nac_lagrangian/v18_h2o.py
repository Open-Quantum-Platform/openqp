"""v18: v7j WITH THE 7.49 ONE-LINE FIX (Sk from the raw product
C0^T S_cross C0). Expected: theory-level closure of J1/J3. Mt harvest, direct-injection
seams, S^x-contraction terms, displaced referees, and ALL judges in a
single process (no cross-frame pairing anywhere).
  exact_U = ytil.w_ref - T1 + gamma:Ux_FD
  prod_U  = -seam(Lt+gamma) + T3pp
plus the gamma-side split and per-piece dumps.
  seam(e_pq) = -U^x_FULL_pq  =>  T2 = -seam(Lt+gamma)
  elimination V-mask = ov-engine weights: 1/2 same-space, 1 on the
  (lo,hi) cross rows; T3 = staged(V) + Tr[dD(V) G[P~]] + gamma:V.

  d[I,J]^c = [amp2e+esum](ytil,X_J)^c            (T1, slot injection)
           + zB[Lt_a(ytil,X_J) + gamma_a]^c      (T2, seam)
           + T3^c                                 (FULL sym channel)
           + gamma:Sk_an^c
  T3^c = staged directional derivative of E(ytil,X_J;C) along
         V_c = -1/2 S^x_MO (frozen-Fock matvec, 2 calls/coord)
       + sum_sigma Tr[dD_sigma(V_c) . G_sigma[P~]]   (Fock density
         response; G built in-process: hf_energy at maxit=1 on D+eps*P~)

Known omission: same-space antisym canonical channel (T4, v7b).
Run: python v7g_h2o.py <input.inp> <dnum.npz>
"""
import os
import sys
import numpy as np

EPSA = 1e-5
TROT = 1e-5
EPSG = 1e-4


def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint
    from scipy.sparse.linalg import minres, LinearOperator
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp, dnum_npz = sys.argv[1], sys.argv[2]
    dn_f = np.load(dnum_npz)
    dcv_n = dn_f['dcv' if 'dcv' in dn_f.files else dn_f.files[0]]

    os.environ['NAC_DUMP_DS'] = '1'
    os.environ['NAC_DUMP_PIJ'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v18.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    ncoord = 3 * natom
    E0 = list(mol.energies)
    Om = [E0[k + 1] - E0[0] for k in range(nstate)]

    ctx = FK.build_context(mol)
    nbf, nij = ctx['nbf'], ctx['nij']
    noca, nocb, nvirb, RS = ctx['noca'], ctx['nocb'], ctx['nvirb'], ctx['RS']
    X0_raw = ctx['X0_raw']
    Xt0 = ctx['Xt']
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a_r = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b_r = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    C = W0.T
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
    ijlr2 = (noca - nocb - 1) * noca + noca
    nbf_tri = nbf * (nbf + 1) // 2

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

    def unpack_sym(pk):
        M = np.zeros((nbf, nbf))
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                M[p, q] = pk[idx]
                M[q, p] = pk[idx]
                idx += 1
        return M

    def pack_sym(M):
        pk = np.zeros(nbf_tri)
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                pk[idx] = M[p, q]
                idx += 1
        return pk

    # gamma + dS exports
    gam = FK.gamma_closed(ctx)
    flatF = np.zeros(nbf * nbf * nstate * nstate)
    for I in range(nstate):
        for J in range(nstate):
            flatF[(I + J * nstate) * nbf * nbf:(I + J * nstate + 1) * nbf * nbf] = \
                gam[I, J].T.reshape(-1)
    mol.data['OQP::nac_gamma_tlf'] = flatF.reshape((nbf * nbf, nstate, nstate))
    oqp.mrsf_nac_overlap(mol)
    nbf2 = nbf * nbf
    dsk_raw = np.array(mol.data['OQP::dbg_dsket'], copy=True).reshape(-1)
    dsf_raw = np.array(mol.data['OQP::dbg_dsfull'], copy=True).reshape(-1)
    Sk_an = np.zeros((ncoord, nbf, nbf))
    Sx_MO = np.zeros((ncoord, nbf, nbf))
    dsf_AO = np.zeros((ncoord, nbf, nbf))
    for c in range(ncoord):
        dsk = dsk_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T
        dsf_AO[c] = dsf_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T
        Sk_an[c] = C.T @ dsk @ C
        Sx_MO[c] = C.T @ dsf_AO[c] @ C

    # elimination direction V_c with the ov-engine weights:
    # V[r,s] = -S^x[r,s] * (1/2 if same-space else 1 if space(r)<space(s) else 0)
    def space_of(i):
        return 0 if i < nocb else (1 if i < noca else 2)
    Wmask = np.zeros((nbf, nbf))
    for rr_ in range(nbf):
        for ss_ in range(nbf):
            sr, sc = space_of(rr_), space_of(ss_)
            if sr == sc:
                Wmask[rr_, ss_] = 0.5
            elif sr < sc:
                Wmask[rr_, ss_] = 1.0
    Vmask = [ -Sx_MO[c] * Wmask for c in range(ncoord) ]

    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    def matvec(v):
        rr = X0_raw.copy().reshape(-1)
        rr[0:nij] = v
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        oqp.mrsf_matvec_apply(mol)
        return np.array(mol.data['OQP::nac_mvax'],
                        copy=True).ravel()[:nij].copy()

    sij0 = FK.s_ij_of(ctx, np.eye(nbf))
    sab0 = FK.s_ab_of(ctx, np.eye(nbf))
    sia0 = FK.s_ia_of(ctx, np.eye(nbf))

    def ampdir_unf(J, dXt):
        Xp = [x.copy() for x in Xt0]
        Xm = [x.copy() for x in Xt0]
        Xp[J] = Xt0[J] + EPSA * dXt
        Xm[J] = Xt0[J] - EPSA * dXt
        Sp = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xp)
        Sm = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xm)
        return (Sp[:, J] - Sm[:, J]) / (2 * EPSA)

    print('G_met sweeps...', flush=True)
    G_met = np.zeros((nstate, nstate, nij))
    for J in range(nstate):
        for k in range(nij):
            if k == ijlr2 - 1:
                continue
            e = np.zeros(nij)
            e[k] = 1.0
            g = ampdir_unf(J, unfold_vec(e))
            for I in range(nstate):
                if I != J:
                    G_met[I, J, k] = g[I]

    phys = [k for k in range(nij) if k != ijlr2 - 1]

    def ytil_minres(I, J):
        XJ = Xf[:, J]

        def op(v):
            vv = np.zeros(nij)
            vv[phys] = v
            vv -= XJ * float(np.dot(XJ, vv))
            av = Om[J] * vv - matvec(vv)
            av -= XJ * float(np.dot(XJ, av))
            return av[phys]

        rhs = G_met[I, J].copy()
        rhs -= XJ * float(np.dot(XJ, rhs))
        y, info = minres(LinearOperator((nij - 1, nij - 1), matvec=op),
                         rhs[phys], rtol=1e-9, maxiter=3000)
        out = np.zeros(nij)
        out[phys] = y
        out -= XJ * float(np.dot(XJ, out))
        return out

    print('MINRES ytil solves...', flush=True)
    ytil = {}
    for I in range(nstate):
        for J in range(nstate):
            if I != J:
                ytil[(I, J)] = ytil_minres(I, J)
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # ---- T1 engines + pij densities --------------------------------------
    def inject(I, vec):
        rr = X0_raw.copy().reshape(-1)
        rr[I * nij:(I + 1) * nij] = vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)

    T1 = {}
    AMP2E = {}
    ESUM = {}
    PIJ = {}
    requested_t1 = os.environ.get('NAC_T1_PAIR', '')
    requested_t1_pair = None
    if requested_t1:
        requested_t1_pair = tuple(
            int(value) - 1 for value in requested_t1.split(',')
        )
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            if requested_t1_pair is not None and (I, J) != requested_t1_pair:
                continue
            inject(I, ytil[(I, J)])
            oqp.mrsf_nac_amp(mol)
            a = np.array(mol.data['OQP::nac_amp'], copy=True
                         ).reshape((nstate, nstate, natom, 3))
            oqp.mrsf_nac_esum(mol, I + 1, J + 1)
            es = np.array(mol.data['OQP::nac_esum'], copy=True).reshape(-1)
            pa = np.array(mol.data['OQP::dbg_pij_a'], copy=True
                          ).reshape(-1).reshape(nbf, nbf).T
            pb = np.array(mol.data['OQP::dbg_pij_b'], copy=True
                          ).reshape(-1).reshape(nbf, nbf).T
            AMP2E[(I, J)] = a[J, I].reshape(-1).copy()
            ESUM[(I, J)] = es.copy()
            T1[(I, J)] = AMP2E[(I, J)] + ESUM[(I, J)]
            PIJ[(I, J)] = (pa.copy(), pb.copy())
            mol.data['OQP::td_bvec_mo'] = X0_raw
    print('T1 engines done.', flush=True)

    if os.environ.get('NAC_T1_ONLY'):
        snapshot = {'Xf': Xf.copy(), 'energies': np.asarray(E0)}
        for I, J in T1:
            snapshot[f'ytil_{I}{J}'] = ytil[(I, J)].copy()
            snapshot[f'amp2e_{I}{J}'] = AMP2E[(I, J)].copy()
            snapshot[f'esum_{I}{J}'] = ESUM[(I, J)].copy()
            snapshot[f'T1_{I}{J}'] = T1[(I, J)].copy()
        output = os.environ.get('NAC_T1_OUTPUT', 't1_snapshot.npz')
        np.savez(output, **snapshot)
        print(f'saved T1 snapshot: {output}', flush=True)
        return

    # Cheap closed-form wpair snapshot for cross-molecule regression.  This
    # exits before the displaced-geometry and O(nbf**2) reference sweeps below.
    # The saved ytil vectors make state/process phase alignment explicit.
    if os.environ.get('NAC_WPAIR_ONLY'):
        requested = os.environ.get('NAC_WPAIR_PAIR', '')
        if requested:
            requested_pair = tuple(
                int(value) - 1 for value in requested.split(',')
            )
            pairs = [requested_pair]
        else:
            pairs = [
                (I, J) for I in range(nstate) for J in range(nstate)
                if I != J
            ]
        snapshot = {'Xf': Xf.copy(), 'energies': np.asarray(E0)}
        for I, J in pairs:
            mol.data['OQP::nac_ytil'] = ytil[(I, J)].copy()
            mol.data['OQP::nac_xstate'] = Xf[:, J].copy()
            oqp.mrsf_nac_wpair(mol, I + 1, J + 1)
            snapshot[f'ytil_{I}{J}'] = ytil[(I, J)].copy()
            closed_mt = np.array(
                mol.data['OQP::nac_mt_frozen'], copy=True
            ).reshape(nbf, nbf).T
            snapshot[f'MT_{I}{J}'] = closed_mt
            if os.environ.get('NAC_WPAIR_DIAG_FD'):
                diagonal_fd = np.full(nbf, np.nan)
                for p in range(nocb, noca):
                    value = 0.0
                    for sign in (1.0, -1.0):
                        rotated = W0.copy()
                        rotated[p, :] += sign * TROT * W0[p, :]
                        mol.data['OQP::VEC_MO_A'] = rotated
                        mol.data['OQP::VEC_MO_B'] = rotated
                        value += sign * float(
                            np.dot(ytil[(I, J)], matvec(Xf[:, J]))
                        ) / (2.0 * TROT)
                    diagonal_fd[p] = value
                mol.data['OQP::VEC_MO_A'] = W0
                mol.data['OQP::VEC_MO_B'] = Wb0
                snapshot[f'MT_diag_fd_{I}{J}'] = diagonal_fd
                error = np.nanmax(np.abs(np.diag(closed_mt) - diagonal_fd))
                print(
                    f'wpair diagonal ({I + 1},{J + 1}) maxdiff={error:.8e}',
                    flush=True,
                )
        output = os.environ.get('NAC_WPAIR_OUTPUT', 'wpair_snapshot.npz')
        np.savez(output, **snapshot)
        print(f'saved closed-form wpair snapshot: {output}', flush=True)
        return

    # ---- displaced w_ref + Ux sweep (exact referee) -----------------------
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

    SAVE0 = {k: np.array(mol.data[k], copy=True) for k in
             ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
              'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A',
              'OQP::E_MO_B']}
    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False
    HD = float(os.environ.get('NAC_HD', '1.0e-3'))

    def displaced_all(coord):
        mol.update_system(coord)
        oqp.library.ints_1e(mol)
        oqp.library.guess(mol)
        SinglePoint(mol).energy()
        mol.data['OQP::xyz_old'] = xyz0.reshape((3, -1))
        mol.data['OQP::VEC_MO_A_old'] = W0
        mol.data['OQP::VEC_MO_B_old'] = Wb0
        mol.data['OQP::E_MO_A_old'] = e0a_r
        mol.data['OQP::E_MO_B_old'] = e0b_r
        oqp.get_structures_ao_overlap(mol)
        M_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        M_f = M_np.reshape(-1).reshape((nbf, nbf)).T
        sg = np.sign(np.diag(M_f))
        sg[sg == 0] = 1.0
        M_f = M_f * sg[None, :]
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass
        Ax_ref = np.zeros((nstate, nij))
        for s in range(nstate):
            Ax_ref[s] = sg_apply(matvec(sg_apply(Xf[:, s], sg)), sg)
        S_np = np.array(mol.data['OQP::overlap_ao_non_orthogonal'], copy=True)
        S_f = S_np.reshape(-1).reshape((nbf, nbf)).T
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        C0m = W0.T
        Skf = C0m.T @ S_f @ C0m          # 7.49 one-line fix
        return M_f, Skf, Ax_ref

    print('displaced w_ref/Ux sweep...', flush=True)
    Ux_FD = np.zeros((ncoord, nbf, nbf))
    w_ref = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        cp = xyz0.copy(); cp[c] += HD
        Mp, Skp, Axp = displaced_all(cp)
        cm = xyz0.copy(); cm[c] -= HD
        Mm, Skm, Axm = displaced_all(cm)
        Ux_FD[c] = (Mp - Mm) / (2 * HD) - (Skp - Skm) / (2 * HD)
        w_ref[c] = (Axp - Axm) / (2 * HD)
    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    # ---- FULL Mt_frozen sweep (single-element directions) -----------------
    print('FULL Mt sweep (single-element directions)...', flush=True)
    MT = {}
    if os.environ.get('NAC_CLOSED_WPAIR'):
        for I in range(nstate):
            for J in range(nstate):
                if I == J:
                    continue
                mol.data['OQP::nac_ytil'] = ytil[(I, J)].copy()
                mol.data['OQP::nac_xstate'] = Xf[:, J].copy()
                oqp.mrsf_nac_wpair(mol, I + 1, J + 1)
                MT[(I, J)] = np.array(
                    mol.data['OQP::nac_mt_frozen'], copy=True
                ).reshape(nbf, nbf).T
                print(
                    f'  closed Mt ({I+1},{J+1}) '
                    f'|Mt|max={np.abs(MT[(I,J)]).max():.3f}',
                    flush=True,
                )
    else:
        for I in range(nstate):
            for J in range(nstate):
                if I == J:
                    continue
                y = ytil[(I, J)]
                M = np.zeros((nbf, nbf))
                for p in range(nbf):
                    for q in range(nbf):
                        acc = 0.0
                        for sgn in (1.0, -1.0):
                            Wp = W0.copy()
                            Wp[q, :] += sgn * TROT * W0[p, :]
                            mol.data['OQP::VEC_MO_A'] = Wp
                            mol.data['OQP::VEC_MO_B'] = Wp
                            acc += sgn * float(
                                np.dot(y, matvec(Xf[:, J]))
                            ) / (2 * TROT)
                        M[p, q] = acc
                mol.data['OQP::VEC_MO_A'] = W0
                mol.data['OQP::VEC_MO_B'] = Wb0
                MT[(I, J)] = M
                print(
                    f'  Mt ({I+1},{J+1}) done '
                    f'|Mt|max={np.abs(M).max():.3f}',
                    flush=True,
                )
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # ---- T3 G-channel -----------------------------------------------------
    SAVE_KEYS = ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
                 'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A',
                 'OQP::E_MO_B']
    saved = {k: np.array(mol.data[k], copy=True) for k in SAVE_KEYS}
    F0a = unpack_sym(saved['OQP::FOCK_A'].ravel())
    F0b = unpack_sym(saved['OQP::FOCK_B'].ravel())
    D0a = unpack_sym(saved['OQP::DM_A'].ravel())
    D0b = unpack_sym(saved['OQP::DM_B'].ravel())
    ok_maxit = False
    for attr in ('maxit', 'scf_maxit', 'maxit_scf'):
        try:
            setattr(mol.data._data.control, attr, 1)
            ok_maxit = True
            print(f'set control.{attr} = 1')
            break
        except Exception:
            continue
    mol.config['scf']['maxit'] = 1
    print(f'maxit control set: {ok_maxit}', flush=True)

    def gbuild(Pa, Pb):
        if os.environ.get('NAC_FULL_RESPONSE'):
            mol.data['OQP::nac_dm1_a'] = pack_sym(Pa)
            mol.data['OQP::nac_dm1_b'] = pack_sym(Pb)
            oqp.mrsf_nac_response(mol)
            return (
                unpack_sym(np.array(
                    mol.data['OQP::nac_v1_a'], copy=True
                ).reshape(-1)),
                unpack_sym(np.array(
                    mol.data['OQP::nac_v1_b'], copy=True
                ).reshape(-1)),
            )
        mol.data['OQP::DM_A'] = pack_sym(D0a + EPSG * Pa).reshape(
            saved['OQP::DM_A'].shape)
        mol.data['OQP::DM_B'] = pack_sym(D0b + EPSG * Pb).reshape(
            saved['OQP::DM_B'].shape)
        oqp.hf_energy(mol)
        Fa = unpack_sym(np.array(mol.data['OQP::FOCK_A'],
                                 copy=True).ravel())
        Fb = unpack_sym(np.array(mol.data['OQP::FOCK_B'],
                                 copy=True).ravel())
        for k in SAVE_KEYS:
            mol.data[k] = saved[k].copy()
        return (Fa - F0a) / EPSG, (Fb - F0b) / EPSG

    print('G-channel matrix builds...', flush=True)
    GSAVE = {}
    MTG = {}
    occ_a = np.zeros(nbf)
    occ_a[:noca] = 1.0
    occ_b = np.zeros(nbf)
    occ_b[:nocb] = 1.0
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            Pa, Pb = PIJ[(I, J)]
            Ga, Gb = gbuild(Pa, Pb)
            GSAVE[(I, J)] = (Ga, Gb)
            GMa = C.T @ Ga @ C
            GMb = C.T @ Gb @ C
            MTG[(I, J)] = 2.0 * (GMa * occ_a[None, :] + GMb * occ_b[None, :])
            print(f'  G ({I+1},{J+1}) |MtG|max={np.abs(MTG[(I,J)]).max():.4f}',
                  flush=True)
    for attr in ('maxit', 'scf_maxit', 'maxit_scf'):
        try:
            setattr(mol.data._data.control, attr, 30)
            break
        except Exception:
            continue
    mol.config['scf']['maxit'] = 30

    # Same-process source/fold referee.  This bypasses both the nuclear CPHF
    # block and the adjoint interchange: ytil.w_ref is the direct derivative
    # of the displaced MRSF sigma action, while gamma:(Sk+U) supplies the
    # overlap/state metric.  It therefore decides whether a residual belongs
    # to the interstate MRSF source or to the ROHF/ROKS nuclear-response side.
    if os.environ.get('NAC_SOURCE_GATE_ONLY'):
        exact_dp = np.zeros((ncoord, nstate, nstate))
        closed_dp = np.zeros_like(exact_dp)
        source_delta = {}
        directional_frozen = {}
        directional_coord = int(os.environ.get(
            'NAC_SOURCE_DIRECTION_COORD', '0'
        ))
        directional_step = float(os.environ.get(
            'NAC_SOURCE_DIRECTION_STEP', '1.0e-6'
        ))
        directional_response = {}
        if os.environ.get('NAC_SOURCE_DIRECTIONAL_AUDIT'):
            direction = Ux_FD[directional_coord]
            dmo_a = direction * occ_a[None, :]
            dmo_a = dmo_a + dmo_a.T
            dmo_b = direction * occ_b[None, :]
            dmo_b = dmo_b + dmo_b.T
            response_a, response_b = gbuild(
                C @ dmo_a @ C.T,
                C @ dmo_b @ C.T,
            )
            for I in range(nstate):
                for J in range(nstate):
                    if I == J:
                        continue
                    value = 0.0
                    for sign in (1.0, -1.0):
                        # W0 stores MO coefficients by row on the Python side.
                        # U[p,q] means C_q <- C_q + theta*C_p, hence the
                        # simultaneous directional rotation is U.T @ W0.
                        rotated = W0 + (
                            sign * directional_step * direction.T @ W0
                        )
                        mol.data['OQP::VEC_MO_A'] = rotated
                        mol.data['OQP::VEC_MO_B'] = rotated
                        value += sign * float(np.dot(
                            ytil[(I, J)], matvec(Xf[:, J])
                        )) / (2.0 * directional_step)
                    directional_frozen[(I, J)] = value
                    predicted = float(np.sum(MT[(I, J)] * direction))
                    probe_a, probe_b = PIJ[(I, J)]
                    response_value = float(
                        np.sum(probe_a * response_a)
                        + np.sum(probe_b * response_b)
                    )
                    directional_response[(I, J)] = response_value
                    predicted_response = float(np.sum(
                        MTG[(I, J)] * direction
                    ))
                    print(
                        f'  frozen-direction ({I+1},{J+1}) coord '
                        f'{directional_coord}: direct={value:.12e} '
                        f'closed={predicted:.12e} '
                        f'diff={predicted-value:.8e}',
                        flush=True,
                    )
                    print(
                        f'  response-direction ({I+1},{J+1}) coord '
                        f'{directional_coord}: direct={response_value:.12e} '
                        f'closed={predicted_response:.12e} '
                        f'diff={predicted_response-response_value:.8e}',
                        flush=True,
                    )
            mol.data['OQP::VEC_MO_A'] = W0
            mol.data['OQP::VEC_MO_B'] = Wb0
            mol.data['OQP::td_bvec_mo'] = X0_raw
        for I in range(nstate):
            for J in range(nstate):
                if I == J:
                    continue
                frozen_direct = np.array([
                    float(np.dot(ytil[(I, J)], w_ref[c, J]))
                    for c in range(ncoord)
                ])
                gsk = np.array([
                    float(np.sum(gam[I, J] * Sk_an[c]))
                    for c in range(ncoord)
                ])
                gu = np.array([
                    float(np.sum(gam[I, J] * Ux_FD[c]))
                    for c in range(ncoord)
                ])
                mfull = MT[(I, J)] + MTG[(I, J)]
                mu = np.array([
                    float(np.sum(mfull * Ux_FD[c]))
                    for c in range(ncoord)
                ])
                exact_dp[:, I, J] = frozen_direct + gsk + gu
                closed_dp[:, I, J] = T1[(I, J)] + mu + gsk + gu
                source_delta[(I, J)] = mu - (frozen_direct - T1[(I, J)])
                print(
                    f'  source ({I+1},{J+1}) '
                    f'max|Mt.U-(y.w-T1)|='
                    f'{np.max(np.abs(source_delta[(I,J)])):.8e}',
                    flush=True,
                )

        def report(label, ordered):
            antisym = 0.5 * (ordered - ordered.transpose(0, 2, 1))
            print(f'===== {label} vs numerical NAC =====', flush=True)
            for I in range(nstate):
                for J in range(I + 1, nstate):
                    value = antisym[:, I, J]
                    reference = dcv_n[I, J].reshape(-1)
                    same = np.max(np.abs(value - reference))
                    flipped = np.max(np.abs(value + reference))
                    print(
                        f'  ({I+1},{J+1}) sign-res maxdiff='
                        f'{min(same, flipped):.8e}',
                        flush=True,
                    )
            return antisym

        exact_dcv = report('direct displaced MRSF source', exact_dp)
        closed_dcv = report('closed-form source', closed_dp)
        output = os.environ.get(
            'NAC_SOURCE_OUTPUT', inp.replace('.inp', '_source_gate.npz')
        )
        payload = {
            'Ux_FD': Ux_FD,
            'w_ref': w_ref,
            'exact_dcv': exact_dcv,
            'closed_dcv': closed_dcv,
            'dcv_reference': dcv_n,
            'energies': np.asarray(E0),
            'Xf': Xf,
        }
        for (I, J), delta in source_delta.items():
            payload[f'source_delta_{I}{J}'] = delta
            payload[f'MT_{I}{J}'] = MT[(I, J)]
            payload[f'MTG_{I}{J}'] = MTG[(I, J)]
            payload[f'T1_{I}{J}'] = T1[(I, J)]
            payload[f'ytil_{I}{J}'] = ytil[(I, J)]
            if (I, J) in directional_frozen:
                payload[f'frozen_directional_{I}{J}'] = np.asarray(
                    directional_frozen[(I, J)]
                )
                payload[f'response_directional_{I}{J}'] = np.asarray(
                    directional_response[(I, J)]
                )
        np.savez(output, **payload)
        print(f'saved source/fold gate: {output}', flush=True)
        return



    # ---- direct-injection seams (X = MT + MTG + gamma) --------------------
    def grad_now():
        oqp.tdhf_mrsf_gradient(mol)
        return np.array(mol.get_grad(), copy=True).reshape(-1)

    print('direct-injection seams...', flush=True)
    T2d = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            Xm = MT[(I, J)] + MTG[(I, J)] + gam[I, J]
            mol.data['OQP::nac_orbgrad_L'] = Xm.ravel(order='F').copy()
            mol.data['OQP::td_bvec_mo'] = X0_raw
            mol.data.set_tdhf_target(J + 1)
            oqp.set_mrsf_nac_cphf(mol, I + 1, J + 1)
            oqp.tdhf_mrsf_z_vector(mol)
            conv = bool(mol.mol_energy.Z_Vector_converged)
            gZ = grad_now()
            mol.data['OQP::td_p'] = np.zeros_like(
                np.array(mol.data['OQP::td_p'], copy=True))
            mol.data['OQP::WAO'] = np.zeros_like(
                np.array(mol.data['OQP::WAO'], copy=True))
            gS = grad_now()
            oqp.set_mrsf_nac_cphf(mol, 0, 0)
            T2d[(I, J)] = gZ - gS
            print(f'  ({I+1},{J+1}) conv={conv} |seam|={np.linalg.norm(gZ-gS):.5f}',
                  flush=True)

    # Diagnostic-only decomposition of the same-process interchange seam.
    # Never use X/U records from another process here: both the SCF orbital
    # gauge and Davidson state phases belong to this run.
    if os.environ.get('NAC_BLOCK_AUDIT'):
        requested = os.environ.get('NAC_BLOCK_PAIR', '3,2')
        pair = tuple(int(value) - 1 for value in requested.split(','))

        def orbital_space(iorb):
            if iorb < nocb:
                return 0
            if iorb < noca:
                return 1
            return 2

        def pack_block(matrix, block):
            result = np.zeros(ncoord)
            for coord, ucoord in enumerate(Ux_FD):
                for p in range(nbf):
                    for q in range(p):
                        sp, sq = orbital_space(p), orbital_space(q)
                        if sp == sq:
                            continue
                        hi, lo = (p, q) if sp > sq else (q, p)
                        this_block = {
                            (1, 0): 1,
                            (2, 0): 2,
                            (2, 1): 3,
                        }[(orbital_space(hi), orbital_space(lo))]
                        if this_block == block:
                            result[coord] += (
                                matrix[hi, lo] - matrix[lo, hi]
                            ) * ucoord[hi, lo]
            return result

        I, J = pair
        Xm = MT[(I, J)] + MTG[(I, J)] + gam[I, J]
        seam_blocks = []
        pack_blocks = []
        print(f'===== SAME-PROCESS SEAM BLOCKS ({I+1},{J+1}) =====',
              flush=True)
        for block, label in (
            (1, 'doc-socc'), (2, 'doc-virt'), (3, 'socc-virt')
        ):
            mol.data['OQP::nac_orbgrad_L'] = Xm.ravel(order='F').copy()
            mol.data['OQP::td_bvec_mo'] = X0_raw
            mol.data.set_tdhf_target(J + 1)
            oqp.set_mrsf_nac_cphf_block(mol, block)
            oqp.set_mrsf_nac_cphf(mol, I + 1, J + 1)
            oqp.tdhf_mrsf_z_vector(mol)
            gZ = grad_now()
            mol.data['OQP::td_p'] = np.zeros_like(
                np.array(mol.data['OQP::td_p'], copy=True))
            mol.data['OQP::WAO'] = np.zeros_like(
                np.array(mol.data['OQP::WAO'], copy=True))
            gS = grad_now()
            seam = gZ - gS
            packed = pack_block(Xm, block)
            seam_blocks.append(seam)
            pack_blocks.append(packed)
            delta = -seam - packed
            print(f'  {label:9s} |seam|={np.linalg.norm(seam):.8e} '
                  f'|pack|={np.linalg.norm(packed):.8e} '
                  f'maxdiff={np.max(np.abs(delta)):.8e} '
                  f'arg={np.argmax(np.abs(delta))}', flush=True)
        seam_sum = np.sum(seam_blocks, axis=0)
        pack_sum = np.sum(pack_blocks, axis=0)
        delta = -seam_sum - pack_sum
        print(f'  sum       |seam|={np.linalg.norm(seam_sum):.8e} '
              f'|pack|={np.linalg.norm(pack_sum):.8e} '
              f'maxdiff={np.max(np.abs(delta)):.8e} '
              f'arg={np.argmax(np.abs(delta))} '
              f'linear={np.max(np.abs(seam_sum-T2d[(I,J)])):.8e}',
              flush=True)
        oqp.set_mrsf_nac_cphf_block(mol, 0)
        oqp.set_mrsf_nac_cphf(mol, 0, 0)

    # First-principles ROHF Lagrangian gate.  This deliberately compares the
    # native response and its adjoint in one orbital/state frame; it never
    # feeds a native CPHF solution through the legacy sfropcal gradient metric.
    if os.environ.get('NAC_LAGRANGIAN_GATE'):
        requested = os.environ.get('NAC_BLOCK_PAIR', '3,2')
        I, J = (int(value) - 1 for value in requested.split(','))
        lmat = MT[(I, J)] + MTG[(I, J)] + gam[I, J]

        def pack_rohf(matrix):
            packed = []
            for ii in range(nocb, noca):
                for jj in range(nocb):
                    packed.append(matrix[ii, jj] - matrix[jj, ii])
            for kk in range(noca, nbf):
                for jj in range(nocb):
                    packed.append(matrix[kk, jj] - matrix[jj, kk])
            for kk in range(noca, nbf):
                for ii in range(nocb, noca):
                    packed.append(matrix[kk, ii] - matrix[ii, kk])
            return np.asarray(packed)

        def pack_u(matrix):
            packed = []
            for ii in range(nocb, noca):
                for jj in range(nocb):
                    packed.append(matrix[ii, jj])
            for kk in range(noca, nbf):
                for jj in range(nocb):
                    packed.append(matrix[kk, jj])
            for kk in range(noca, nbf):
                for ii in range(nocb, noca):
                    packed.append(matrix[kk, ii])
            return np.asarray(packed)

        os.environ['NAC_DUMP_ROHF_RESPONSE'] = '1'
        oqp.hf_hessian(mol)

        def rotation_by_cartesian(tag):
            raw = np.array(mol.data[tag], copy=True)
            if raw.ndim != 2:
                raise RuntimeError(f'{tag} is not a 2-D TagArray: {raw.shape}')
            if raw.size % ncoord:
                raise RuntimeError(f'{tag} size is not divisible by '
                                   f'{ncoord}: {raw.shape}')
            # OQPData exposes the Fortran buffer through a C-order reshape
            # using the dimensions recorded by Fortran.  For a non-square
            # Fortran array A(rotation,Cartesian), neither ``raw`` nor
            # ``raw.T`` is A: its flat storage must first be reshaped as the
            # reversed dimensions, then transposed.  (The usual bare ``.T``
            # happened to work only for square TagArrays.)
            nrotation = raw.size // ncoord
            return raw.reshape(ncoord, nrotation).T

        b_hf = rotation_by_cartesian('OQP::nac_rohf_bvec_hf_jk_pulay')
        b_full = rotation_by_cartesian('OQP::nac_rohf_bvec_full')
        u_native = rotation_by_cartesian('OQP::nac_rohf_uvec')
        u_fd = np.column_stack([pack_u(Ux_FD[c]) for c in range(ncoord)])

        nds = (noca - nocb) * nocb
        ndv = (nbf - noca) * nocb
        blocks = (
            ('doc-socc', slice(0, nds)),
            ('doc-virt', slice(nds, nds + ndv)),
            ('socc-virt', slice(nds + ndv, None)),
        )
        print(f'===== NATIVE ROHF LAGRANGIAN GATE ({I+1},{J+1}) =====',
              flush=True)
        for label, block in blocks:
            delta = u_native[block] - u_fd[block]
            print(f'  U {label:9s} maxdiff={np.max(np.abs(delta)):.8e} '
                  f'rms={np.sqrt(np.mean(delta**2)):.8e}', flush=True)

        lvec = pack_rohf(lmat)
        mol.data['OQP::nac_rohf_rhs'] = lvec.copy()
        oqp.mrsf_nac_rohf_zvector(mol)
        zvec = np.array(
            mol.data['OQP::nac_rohf_solution'], copy=True
        ).reshape(-1)
        direct = lvec @ u_native
        direct_fd = lvec @ u_fd
        adjoint = zvec @ b_full
        print(f'  native-U vs FD-U contraction maxdiff='
              f'{np.max(np.abs(direct-direct_fd)):.8e} '
              f'|FD|={np.linalg.norm(direct_fd):.8e}', flush=True)
        print(f'  adjoint maxdiff={np.max(np.abs(direct-adjoint)):.8e} '
              f'|direct|={np.linalg.norm(direct):.8e}', flush=True)

        # The old nuclear RHS obtains this term from dftexcor(R+/-h).  The new
        # path is the analytic moving-grid mixed derivative
        # -1/2 d_R Tr[Vxc delta P_z].
        mol.data['OQP::nac_rohf_z'] = zvec.copy()
        oqp.mrsf_nac_rohf_hf_adjoint(mol)
        hf_analytic = np.array(
            mol.data['OQP::nac_rohf_hf_adjoint'], copy=True
        ).reshape(-1)
        hf_forward = zvec @ b_hf
        print(f'  HF adjoint-vs-forward-RHS maxdiff='
              f'{np.max(np.abs(hf_analytic-hf_forward)):.8e} '
              f'|HF|={np.linalg.norm(hf_analytic):.8e}', flush=True)
        oqp.mrsf_nac_xc_adjoint(mol)
        xc_analytic = np.array(
            mol.data['OQP::nac_rohf_xc_adjoint'], copy=True
        ).reshape(-1)
        xc_fd_rhs = zvec @ (b_full - b_hf)
        print(f'  XC analytic-vs-FD-RHS maxdiff='
              f'{np.max(np.abs(xc_analytic-xc_fd_rhs)):.8e} '
              f'|XC|={np.linalg.norm(xc_analytic):.8e}', flush=True)
        print('  XC analytic = ' + np.array2string(
            xc_analytic, precision=9, suppress_small=False), flush=True)
        print('  XC FD RHS   = ' + np.array2string(
            xc_fd_rhs, precision=9, suppress_small=False), flush=True)
        full_analytic = hf_analytic + xc_analytic
        full_forward = zvec @ b_full
        print(f'  FULL Z-vector analytic-vs-forward-RHS maxdiff='
              f'{np.max(np.abs(full_analytic-full_forward)):.8e} '
              f'|FULL|={np.linalg.norm(full_analytic):.8e}', flush=True)
        del os.environ['NAC_DUMP_ROHF_RESPONSE']
        if os.environ.get('NAC_LAGRANGIAN_GATE_ONLY'):
            return
    try:
        del mol.data['OQP::nac_orbgrad_L']
    except Exception:
        pass
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # ---- displaced w_ref + Ux sweep (exact referee) -----------------------
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

    SAVE0 = {k: np.array(mol.data[k], copy=True) for k in
             ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
              'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A',
              'OQP::E_MO_B']}
    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False
    HD = float(os.environ.get('NAC_HD', '1.0e-3'))

    def displaced_all(coord):
        mol.update_system(coord)
        oqp.library.ints_1e(mol)
        oqp.library.guess(mol)
        SinglePoint(mol).energy()
        mol.data['OQP::xyz_old'] = xyz0.reshape((3, -1))
        mol.data['OQP::VEC_MO_A_old'] = W0
        mol.data['OQP::VEC_MO_B_old'] = Wb0
        mol.data['OQP::E_MO_A_old'] = e0a_r
        mol.data['OQP::E_MO_B_old'] = e0b_r
        oqp.get_structures_ao_overlap(mol)
        M_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        M_f = M_np.reshape(-1).reshape((nbf, nbf)).T
        sg = np.sign(np.diag(M_f))
        sg[sg == 0] = 1.0
        M_f = M_f * sg[None, :]
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass
        Ax_ref = np.zeros((nstate, nij))
        for s in range(nstate):
            Ax_ref[s] = sg_apply(matvec(sg_apply(Xf[:, s], sg)), sg)
        S_np = np.array(mol.data['OQP::overlap_ao_non_orthogonal'], copy=True)
        S_f = S_np.reshape(-1).reshape((nbf, nbf)).T
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        C0m = W0.T
        Skf = C0m.T @ S_f @ C0m          # 7.49 one-line fix
        return M_f, Skf, Ax_ref

    print('displaced w_ref/Ux sweep...', flush=True)
    Ux_FD = np.zeros((ncoord, nbf, nbf))
    w_ref = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        cp = xyz0.copy(); cp[c] += HD
        Mp, Skp, Axp = displaced_all(cp)
        cm = xyz0.copy(); cm[c] -= HD
        Mm, Skm, Axm = displaced_all(cm)
        Ux_FD[c] = (Mp - Mm) / (2 * HD) - (Skp - Skm) / (2 * HD)
        w_ref[c] = (Axp - Axm) / (2 * HD)
    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    # ---- FULL Mt_frozen sweep (single-element directions) -----------------
    print('FULL Mt sweep (single-element directions)...', flush=True)
    MT = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            y = ytil[(I, J)]
            M = np.zeros((nbf, nbf))
            for p in range(nbf):
                for q in range(nbf):
                    acc = 0.0
                    for sgn in (1.0, -1.0):
                        Wp = W0.copy()
                        Wp[q, :] += sgn * TROT * W0[p, :]
                        mol.data['OQP::VEC_MO_A'] = Wp
                        mol.data['OQP::VEC_MO_B'] = Wp
                        acc += sgn * float(np.dot(y, matvec(Xf[:, J]))) / (2 * TROT)
                    M[p, q] = acc
            mol.data['OQP::VEC_MO_A'] = W0
            mol.data['OQP::VEC_MO_B'] = Wb0
            MT[(I, J)] = M
            print(f'  Mt ({I+1},{J+1}) done |Mt|max={np.abs(M).max():.3f}',
                  flush=True)
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # ---- T3 G-channel -----------------------------------------------------
    SAVE_KEYS = ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
                 'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A',
                 'OQP::E_MO_B']
    saved = {k: np.array(mol.data[k], copy=True) for k in SAVE_KEYS}
    F0a = unpack_sym(saved['OQP::FOCK_A'].ravel())
    F0b = unpack_sym(saved['OQP::FOCK_B'].ravel())
    D0a = unpack_sym(saved['OQP::DM_A'].ravel())
    D0b = unpack_sym(saved['OQP::DM_B'].ravel())
    ok_maxit = False
    for attr in ('maxit', 'scf_maxit', 'maxit_scf'):
        try:
            setattr(mol.data._data.control, attr, 1)
            ok_maxit = True
            print(f'set control.{attr} = 1')
            break
        except Exception:
            continue
    mol.config['scf']['maxit'] = 1
    print(f'maxit control set: {ok_maxit}', flush=True)

    def gbuild(Pa, Pb):
        if os.environ.get('NAC_FULL_RESPONSE'):
            mol.data['OQP::nac_dm1_a'] = pack_sym(Pa)
            mol.data['OQP::nac_dm1_b'] = pack_sym(Pb)
            oqp.mrsf_nac_response(mol)
            return (
                unpack_sym(np.array(
                    mol.data['OQP::nac_v1_a'], copy=True
                ).reshape(-1)),
                unpack_sym(np.array(
                    mol.data['OQP::nac_v1_b'], copy=True
                ).reshape(-1)),
            )
        mol.data['OQP::DM_A'] = pack_sym(D0a + EPSG * Pa).reshape(
            saved['OQP::DM_A'].shape)
        mol.data['OQP::DM_B'] = pack_sym(D0b + EPSG * Pb).reshape(
            saved['OQP::DM_B'].shape)
        oqp.hf_energy(mol)
        Fa = unpack_sym(np.array(mol.data['OQP::FOCK_A'],
                                 copy=True).ravel())
        Fb = unpack_sym(np.array(mol.data['OQP::FOCK_B'],
                                 copy=True).ravel())
        for k in SAVE_KEYS:
            mol.data[k] = saved[k].copy()
        return (Fa - F0a) / EPSG, (Fb - F0b) / EPSG

    print('G-channel matrix builds...', flush=True)
    GSAVE = {}
    MTG = {}
    occ_a = np.zeros(nbf)
    occ_a[:noca] = 1.0
    occ_b = np.zeros(nbf)
    occ_b[:nocb] = 1.0
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            Pa, Pb = PIJ[(I, J)]
            Ga, Gb = gbuild(Pa, Pb)
            GSAVE[(I, J)] = (Ga, Gb)
            GMa = C.T @ Ga @ C
            GMb = C.T @ Gb @ C
            MTG[(I, J)] = 2.0 * (GMa * occ_a[None, :] + GMb * occ_b[None, :])
            print(f'  G ({I+1},{J+1}) |MtG|max={np.abs(MTG[(I,J)]).max():.4f}',
                  flush=True)
    for attr in ('maxit', 'scf_maxit', 'maxit_scf'):
        try:
            setattr(mol.data._data.control, attr, 30)
            break
        except Exception:
            continue
    mol.config['scf']['maxit'] = 30



    os.environ['NAC_DUMP_RHS'] = '1'

    def rhs_of(target, amp_vec):
        rr = X0_raw.copy().reshape(-1)
        rr[(target - 1) * nij:target * nij] = amp_vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        mol.data.set_tdhf_target(target)
        oqp.tdhf_mrsf_z_vector(mol)
        out = np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).reshape(-1)
        mol.data['OQP::td_bvec_mo'] = X0_raw
        return out

    def unpack_rot_to_mat(v):
        Lm = np.zeros((nbf, nbf))
        idx = 0
        for ii in range(nocb, noca):
            for jj in range(nocb):
                Lm[ii, jj] = v[idx]; idx += 1
        for kk in range(noca, nbf):
            for jj in range(nocb):
                Lm[kk, jj] = v[idx]; idx += 1
        for kk in range(noca, nbf):
            for ii in range(nocb, noca):
                Lm[kk, ii] = v[idx]; idx += 1
        return Lm

    def grad_now():
        oqp.tdhf_mrsf_gradient(mol)
        return np.array(mol.get_grad(), copy=True).reshape(-1)

    print('T2 seams...', flush=True)
    T2 = {}
    T2g_only = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            y = ytil[(I, J)]
            qy = rhs_of(J + 1, y)
            qx = rhs_of(J + 1, Xf[:, J])
            qm = rhs_of(J + 1, y + Xf[:, J])
            Lpol = 0.5 * (qm - qy - qx)
            Lmat = unpack_rot_to_mat(Lpol) + gam[I, J]
            mol.data['OQP::nac_orbgrad_L'] = Lmat.ravel(order='F').copy()
            mol.data['OQP::td_bvec_mo'] = X0_raw
            mol.data.set_tdhf_target(J + 1)
            oqp.set_mrsf_nac_cphf(mol, I + 1, J + 1)
            oqp.tdhf_mrsf_z_vector(mol)
            if not bool(mol.mol_energy.Z_Vector_converged):
                oqp.set_mrsf_nac_cphf(mol, 0, 0)
                T2[(I, J)] = None
                print(f'  ({I+1},{J+1}): z NOT converged', flush=True)
                continue
            gZ = grad_now()
            mol.data['OQP::td_p'] = np.zeros_like(
                np.array(mol.data['OQP::td_p'], copy=True))
            mol.data['OQP::WAO'] = np.zeros_like(
                np.array(mol.data['OQP::WAO'], copy=True))
            gS = grad_now()
            oqp.set_mrsf_nac_cphf(mol, 0, 0)
            T2[(I, J)] = gZ - gS
            # gamma-only seam
            mol.data['OQP::nac_orbgrad_L'] = gam[I, J].ravel(order='F').copy()
            mol.data['OQP::td_bvec_mo'] = X0_raw
            mol.data.set_tdhf_target(J + 1)
            oqp.set_mrsf_nac_cphf(mol, I + 1, J + 1)
            oqp.tdhf_mrsf_z_vector(mol)
            gZ2 = grad_now()
            mol.data['OQP::td_p'] = np.zeros_like(
                np.array(mol.data['OQP::td_p'], copy=True))
            mol.data['OQP::WAO'] = np.zeros_like(
                np.array(mol.data['OQP::WAO'], copy=True))
            gS2 = grad_now()
            oqp.set_mrsf_nac_cphf(mol, 0, 0)
            T2g_only[(I, J)] = gZ2 - gS2
    try:
        del mol.data['OQP::nac_orbgrad_L']
    except Exception:
        pass

    # ---- ALL JUDGES (one process, one frame) ------------------------------
    def space_of(i):
        return 0 if i < nocb else (1 if i < noca else 2)

    print('\n===== J1: Mt-completeness (in-process) =====')
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            Mfull = MT[(I, J)] + MTG[(I, J)]
            lhs = np.array([float(np.sum(Mfull * Ux_FD[c]))
                            for c in range(ncoord)])
            rhs = np.array([float(np.dot(ytil[(I, J)], w_ref[c, J]))
                            for c in range(ncoord)]) - T1[(I, J)]
            print(f'({I+1},{J+1}): |Mt:U|={np.linalg.norm(lhs):.5f} '
                  f'|exact|={np.linalg.norm(rhs):.5f} '
                  f'maxdiff={np.abs(lhs-rhs).max():.3e}')

    print('\n===== J4: seam identity -seam vs pack(X).U =====')
    tpack = {}
    telim = {}
    tsssym = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            X = MT[(I, J)] + MTG[(I, J)] + gam[I, J]
            tp = np.zeros(ncoord)
            te = np.zeros(ncoord)
            ts = np.zeros(ncoord)
            for c in range(ncoord):
                a = e = s1 = 0.0
                U = Ux_FD[c]
                for p in range(nbf):
                    for q in range(p):
                        sp, sq = space_of(p), space_of(q)
                        if sp != sq:
                            hi, lo = (p, q) if sp > sq else (q, p)
                            a += (X[hi, lo] - X[lo, hi]) * U[hi, lo]
                            e += -X[lo, hi] * Sx_MO[c][hi, lo]
                        else:
                            s1 += (X[p, q] + X[q, p]) * (-0.5) * Sx_MO[c][p, q]
                tp[c] = a; te[c] = e; ts[c] = s1
            tpack[(I, J)] = tp
            telim[(I, J)] = te
            tsssym[(I, J)] = ts
            print(f'({I+1},{J+1}): |seam|={np.linalg.norm(T2d[(I,J)]):.5f} '
                  f'|pack|={np.linalg.norm(tp):.5f} '
                  f'maxdiff={np.abs(-T2d[(I,J)]-tp).max():.3e}')

    dn_all = dcv_n
    print('\n===== J3: total T1 + X:Ux + gSk vs d_num =====')
    def judge(dp, tag):
        dpa = 0.5 * (dp - dp.transpose(0, 2, 1))
        for I in range(nstate):
            for J in range(I + 1, nstate):
                v = dpa[:, I, J]
                ref = dn_all[I, J].reshape(-1)
                cc = float(np.dot(ref, v)) / (np.linalg.norm(ref)
                                              * np.linalg.norm(v) + 1e-300)
                md = min(np.abs(v - ref).max(), np.abs(v + ref).max())
                print(f'[{tag}] ({I+1},{J+1}): |d_num|={np.linalg.norm(ref):.6f} '
                      f'|pred|={np.linalg.norm(v):.6f} cos={cc:+.8f} '
                      f'sign-res maxdiff={md:.3e}')

    dp3 = np.zeros((ncoord, nstate, nstate))
    dp2 = np.zeros((ncoord, nstate, nstate))
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            X = MT[(I, J)] + MTG[(I, J)] + gam[I, J]
            gsk = np.array([float(np.sum(gam[I, J] * Sk_an[c]))
                            for c in range(ncoord)])
            uch = np.array([float(np.sum(X * Ux_FD[c])) for c in range(ncoord)])
            dp3[:, I, J] = T1[(I, J)] + uch + gsk
            dp2[:, I, J] = (T1[(I, J)] - T2d[(I, J)] + telim[(I, J)]
                            + tsssym[(I, J)] + gsk)
    judge(dp3, 'FD-U total')
    print('\n===== J2: PRODUCTION FORM T1 - seam + elim + ss_sym + gSk =====')
    judge(dp2, 'production')

    out = dict(Ux_FD=Ux_FD, w_ref=w_ref, gam=gam, Sk_an=Sk_an, Sx_MO=Sx_MO,
               Om=np.array(Om), Xf=Xf,
               noca=np.array(noca), nocb=np.array(nocb))
    for (I, J) in MT:
        out[f'MT_{I}{J}'] = MT[(I, J)]
        out[f'MTG_{I}{J}'] = MTG[(I, J)]
        out[f'T1_{I}{J}'] = T1[(I, J)]
        out[f'T2d_{I}{J}'] = T2d[(I, J)]
        out[f'ytil_{I}{J}'] = ytil[(I, J)]
    np.savez(inp.replace('.inp', '_v18.npz'), **out)
    print('saved v7j npz.', flush=True)


if __name__ == '__main__':
    main()
