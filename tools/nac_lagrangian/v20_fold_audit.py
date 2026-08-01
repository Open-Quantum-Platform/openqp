"""v20: identify and gate the fold-sector remainder in the v19 F channel.

The MRSF sigma source has only one ground-state-Fock dependency:
``mrsfesum(fa, fb)``.  Therefore compare the v19 residual with the response
obtained by injecting the *actual* converged-minus-one-iteration Fock
derivative.  In parallel, compare that Fock derivative with G[dD] from both
the orbital model and the finite-difference SCF density.

Driver:
    python v20_fold_audit.py <inp> <coord> <skeldir>

The driver launches two private workers with ``--worker``.  Worker outputs use
the reference Davidson phases written by the driver; do not reuse them across
processes.
"""

import os
import subprocess
import sys

import numpy as np


HD = 1.0e-3
EPSG = 1.0e-4
EPSF = 1.0e-4


def _unpack_sym(pk, nbf):
    matrix = np.zeros((nbf, nbf))
    idx = 0
    for q in range(nbf):
        for p in range(q + 1):
            matrix[p, q] = pk[idx]
            matrix[q, p] = pk[idx]
            idx += 1
    return matrix


def _pack_sym(matrix):
    nbf = matrix.shape[0]
    packed = np.zeros(nbf * (nbf + 1) // 2)
    idx = 0
    for q in range(nbf):
        for p in range(q + 1):
            packed[idx] = matrix[p, q]
            idx += 1
    return packed


def worker(inp, ref_npz, out_npz):
    import oqp
    from oqp.pyoqp import Runner

    ref = np.load(ref_npz)
    runner = Runner(input_file=inp, log=inp.replace('.inp', '_v20.log'))
    runner.run()
    mol = runner.mol
    nstate = int(ref['nstate'])
    nij = int(ref['nij'])
    x0_raw = ref['X0_raw']
    xj = x0_raw.reshape(-1).reshape((nstate, nij)).T[:, int(ref['J'])]

    # For maxit=1 these Fock records were built from the JSON-guess density,
    # before the single orbital/density update at the end of the iteration.
    fock_a = np.array(mol.data['OQP::FOCK_A'], copy=True)
    fock_b = np.array(mol.data['OQP::FOCK_B'], copy=True)
    dm_a = np.array(mol.data['OQP::DM_A'], copy=True)
    dm_b = np.array(mol.data['OQP::DM_B'], copy=True)
    vec_a = np.array(mol.data['OQP::VEC_MO_A'], copy=True)

    mol.data['OQP::VEC_MO_A'] = ref['W0'].copy()
    mol.data['OQP::VEC_MO_B'] = ref['Wb0'].copy()
    rr = x0_raw.copy().reshape(-1)
    rr[:nij] = xj
    mol.data['OQP::td_bvec_mo'] = rr.reshape(x0_raw.shape)
    try:
        mol.data._data.control.int2e_cutoff = 1.0e-20
    except Exception:
        pass
    oqp.mrsf_matvec_apply(mol)
    ax = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()[:nij]
    np.savez(out_npz, Ax=ax, F_a=fock_a, F_b=fock_b,
             D_a=dm_a, D_b=dm_b, W_a=vec_a)


def main(inp, c0, skeldir):
    import oqp
    import oqp.library
    from oqp.library.single_point import SinglePoint
    from oqp.pyoqp import Runner

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    os.environ['NAC_DUMP_PIJ'] = '1'

    runner = Runner(input_file=inp, log=inp.replace('.inp', '_v20.log'))
    runner.run()
    mol = runner.mol
    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    ncoord = 3 * natom
    ctx = FK.build_context(mol)
    nbf, nij = ctx['nbf'], ctx['nij']
    noca, nocb = ctx['noca'], ctx['nocb']
    x0_raw = ctx['X0_raw']
    xt0 = ctx['Xt']
    rs = ctx['RS']
    xf = x0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    jstate = 2
    xj = xf[:, jstate]
    ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
    ijlr2 = (noca - nocb - 1) * noca + noca

    w0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    cmat0 = w0.T
    keys = ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
            'OQP::VEC_MO_A', 'OQP::VEC_MO_B',
            'OQP::E_MO_A', 'OQP::E_MO_B']
    save0 = {key: np.array(mol.data[key], copy=True) for key in keys}
    d0a = _unpack_sym(save0['OQP::DM_A'].ravel(), nbf)
    d0b = _unpack_sym(save0['OQP::DM_B'].ravel(), nbf)
    f0a = _unpack_sym(save0['OQP::FOCK_A'].ravel(), nbf)
    f0b = _unpack_sym(save0['OQP::FOCK_B'].ravel(), nbf)

    try:
        mol.data._data.control.int2e_cutoff = 1.0e-20
    except Exception:
        pass

    def restore():
        mol.update_system(xyz0)
        oqp.library.ints_1e(mol)
        for key, value in save0.items():
            mol.data[key] = value.copy()
        mol.data['OQP::td_bvec_mo'] = x0_raw
        try:
            mol.data._data.control.int2e_cutoff = 1.0e-20
        except Exception:
            pass

    def matvec(vector):
        rr = x0_raw.copy().reshape(-1)
        rr[:nij] = vector
        mol.data['OQP::td_bvec_mo'] = rr.reshape(x0_raw.shape)
        oqp.mrsf_matvec_apply(mol)
        return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()[:nij]

    base = matvec(xj)

    # Workers use this process's Davidson phases and the already prepared v19
    # one-iteration displaced inputs.
    ref_npz = os.path.join(skeldir, 'ref_v20.npz')
    np.savez(ref_npz, W0=w0, Wb0=wb0, X0_raw=x0_raw, nstate=nstate,
             nij=nij, J=jstate)
    env = dict(os.environ, OMP_NUM_THREADS='4')
    worker_outputs = {}
    processes = []
    for idx, label in ((c0, 'p'), (c0 + ncoord, 'm')):
        pinp = os.path.join(skeldir, f'p{idx}.inp')
        out = os.path.join(skeldir, f'p{idx}_v20.npz')
        worker_outputs[label] = out
        processes.append(subprocess.Popen(
            [sys.executable, os.path.abspath(__file__), '--worker',
             pinp, ref_npz, out], env=env,
            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT))

    mol.save_data()
    mol.config['guess']['type'] = 'json'
    mol.config['guess']['file'] = mol.log.replace('.log', '.json')
    mol.config['guess']['continue_geom'] = False

    def displaced(sign):
        coord = xyz0.copy()
        coord[c0] += sign * HD
        mol.update_system(coord)
        oqp.library.ints_1e(mol)
        oqp.library.guess(mol)
        SinglePoint(mol).energy()

        fock_a = np.array(mol.data['OQP::FOCK_A'], copy=True)
        fock_b = np.array(mol.data['OQP::FOCK_B'], copy=True)
        dm_a = np.array(mol.data['OQP::DM_A'], copy=True)
        dm_b = np.array(mol.data['OQP::DM_B'], copy=True)

        mol.data['OQP::xyz_old'] = xyz0.reshape((3, -1))
        mol.data['OQP::VEC_MO_A_old'] = w0
        mol.data['OQP::VEC_MO_B_old'] = wb0
        mol.data['OQP::E_MO_A_old'] = e0a
        mol.data['OQP::E_MO_B_old'] = e0b
        oqp.get_structures_ao_overlap(mol)
        mmat = np.array(
            mol.data['OQP::overlap_mo_non_orthogonal'], copy=True
        ).reshape((nbf, nbf)).T
        smat_cross = np.array(
            mol.data['OQP::overlap_ao_non_orthogonal'], copy=True
        ).reshape((nbf, nbf)).T
        gauge = np.sign(np.diag(mmat))
        gauge[gauge == 0.0] = 1.0
        mmat *= gauge[None, :]
        wd = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
        wd_sf = gauge[:, None] * wd
        mol.data['OQP::VEC_MO_A'] = wd_sf
        mol.data['OQP::VEC_MO_B'] = wd_sf.copy()
        ax = matvec(xj)
        skmat = cmat0.T @ smat_cross @ cmat0
        return dict(Ax=ax, M=mmat, Sk=skmat, F_a=fock_a, F_b=fock_b,
                    D_a=dm_a, D_b=dm_b, W_sf=wd_sf)

    plus = displaced(+1.0)
    minus = displaced(-1.0)
    restore()

    for process in processes:
        if process.wait() != 0:
            raise RuntimeError('v20 worker failed')
    worker_plus = np.load(worker_outputs['p'])
    worker_minus = np.load(worker_outputs['m'])

    w_ref = (plus['Ax'] - minus['Ax']) / (2.0 * HD)
    w_skel = (worker_plus['Ax'] - worker_minus['Ax']) / (2.0 * HD)
    ux = ((plus['M'] - minus['M']) -
          (plus['Sk'] - minus['Sk'])) / (2.0 * HD)

    # Frozen-integral orbital response at the reference geometry.
    step = 1.0e-5
    staged_c = np.zeros(nij)
    for sign in (+1.0, -1.0):
        wp = (np.eye(nbf) + sign * step * ux.T) @ w0
        mol.data['OQP::VEC_MO_A'] = wp
        mol.data['OQP::VEC_MO_B'] = wp
        staged_c += sign * matvec(xj) / (2.0 * step)
    restore()
    true_f = w_ref - w_skel - staged_c

    def inject_fock(dfa, dfb):
        mol.data['OQP::FOCK_A'] = _pack_sym(
            f0a + EPSF * dfa
        ).reshape(save0['OQP::FOCK_A'].shape)
        mol.data['OQP::FOCK_B'] = _pack_sym(
            f0b + EPSF * dfb
        ).reshape(save0['OQP::FOCK_B'].shape)
        response = (matvec(xj) - base) / EPSF
        restore()
        return response

    # This is the exact record-level F channel sampled by true_f.
    dfa_actual = (
        (_unpack_sym(plus['F_a'].ravel(), nbf) -
         _unpack_sym(worker_plus['F_a'].ravel(), nbf)) -
        (_unpack_sym(minus['F_a'].ravel(), nbf) -
         _unpack_sym(worker_minus['F_a'].ravel(), nbf))
    ) / (2.0 * HD)
    dfb_actual = (
        (_unpack_sym(plus['F_b'].ravel(), nbf) -
         _unpack_sym(worker_plus['F_b'].ravel(), nbf)) -
        (_unpack_sym(minus['F_b'].ravel(), nbf) -
         _unpack_sym(worker_minus['F_b'].ravel(), nbf))
    ) / (2.0 * HD)
    gvec_actual = inject_fock(dfa_actual, dfb_actual)

    # Compare the finite-difference converged density with the orbital formula.
    dda_fd = (_unpack_sym(plus['D_a'].ravel(), nbf) -
              _unpack_sym(minus['D_a'].ravel(), nbf)) / (2.0 * HD)
    ddb_fd = (_unpack_sym(plus['D_b'].ravel(), nbf) -
              _unpack_sym(minus['D_b'].ravel(), nbf)) / (2.0 * HD)
    occ_a = np.zeros(nbf)
    occ_b = np.zeros(nbf)
    occ_a[:noca] = 1.0
    occ_b[:nocb] = 1.0
    ma = ux * occ_a[None, :] + ux.T * occ_a[:, None]
    mb = ux * occ_b[None, :] + ux.T * occ_b[:, None]
    dda_model = cmat0 @ ma @ cmat0.T
    ddb_model = cmat0 @ mb @ cmat0.T

    def gbuild(dda, ddb):
        for attr in ('maxit', 'scf_maxit'):
            try:
                setattr(mol.data._data.control, attr, 1)
                break
            except Exception:
                continue
        mol.config['scf']['maxit'] = 1
        mol.data['OQP::DM_A'] = _pack_sym(
            d0a + EPSG * dda
        ).reshape(save0['OQP::DM_A'].shape)
        mol.data['OQP::DM_B'] = _pack_sym(
            d0b + EPSG * ddb
        ).reshape(save0['OQP::DM_B'].shape)
        oqp.hf_energy(mol)
        ga = (_unpack_sym(
            np.array(mol.data['OQP::FOCK_A'], copy=True).ravel(), nbf
        ) - f0a) / EPSG
        gb = (_unpack_sym(
            np.array(mol.data['OQP::FOCK_B'], copy=True).ravel(), nbf
        ) - f0b) / EPSG
        restore()
        return ga, gb

    def gbuild_consistent(dda, ddb):
        """Differentiate F[D(C), C] by moving DM and MO records together.

        ``calc_fock`` uses DM_A/B for J/K but ``calc_dft_xc`` passes VEC_MO
        to ``dftexcor``.  A DM-only perturbation therefore omits f_xc[dD].
        """
        responses = []
        for sign in (+1.0, -1.0):
            for attr in ('maxit', 'scf_maxit'):
                try:
                    setattr(mol.data._data.control, attr, 1)
                    break
                except Exception:
                    continue
            mol.config['scf']['maxit'] = 1
            mol.data['OQP::DM_A'] = _pack_sym(
                d0a + sign * EPSG * dda
            ).reshape(save0['OQP::DM_A'].shape)
            mol.data['OQP::DM_B'] = _pack_sym(
                d0b + sign * EPSG * ddb
            ).reshape(save0['OQP::DM_B'].shape)
            wp = (np.eye(nbf) + sign * EPSG * ux.T) @ w0
            mol.data['OQP::VEC_MO_A'] = wp
            mol.data['OQP::VEC_MO_B'] = wp.copy()
            oqp.hf_energy(mol)
            responses.append((
                _unpack_sym(np.array(
                    mol.data['OQP::FOCK_A'], copy=True
                ).ravel(), nbf),
                _unpack_sym(np.array(
                    mol.data['OQP::FOCK_B'], copy=True
                ).ravel(), nbf),
            ))
            restore()
        ga = (responses[0][0] - responses[1][0]) / (2.0 * EPSG)
        gb = (responses[0][1] - responses[1][1]) / (2.0 * EPSG)
        return ga, gb

    ga_model, gb_model = gbuild(dda_model, ddb_model)
    ga_fd, gb_fd = gbuild(dda_fd, ddb_fd)
    ga_full, gb_full = gbuild_consistent(dda_model, ddb_model)
    mol.data['OQP::nac_dm1_a'] = _pack_sym(dda_model)
    mol.data['OQP::nac_dm1_b'] = _pack_sym(ddb_model)
    oqp.mrsf_nac_response(mol)
    ga_analytic = _unpack_sym(
        np.array(mol.data['OQP::nac_v1_a'], copy=True).ravel(), nbf
    )
    gb_analytic = _unpack_sym(
        np.array(mol.data['OQP::nac_v1_b'], copy=True).ravel(), nbf
    )
    ga_xc_analytic = _unpack_sym(
        np.array(mol.data['OQP::nac_vxc_a'], copy=True).ravel(), nbf
    )
    gb_xc_analytic = _unpack_sym(
        np.array(mol.data['OQP::nac_vxc_b'], copy=True).ravel(), nbf
    )
    gvec_model = inject_fock(ga_model, gb_model)
    gvec_fd = inject_fock(ga_fd, gb_fd)
    gvec_full = inject_fock(ga_full, gb_full)
    gvec_analytic = inject_fock(ga_analytic, gb_analytic)
    gvec_xc = inject_fock(ga_full - ga_model, gb_full - gb_model)

    # Convert the vector-level XC closure into the actual pair/coordinate
    # correction ytil_IJ . Delta-sigma_xc for the two pairs ending at J=3.
    from scipy.sparse.linalg import LinearOperator, minres

    def unfold_vec(vector):
        matrix = np.zeros((noca, nbf - nocb))
        for iorb in range(1, noca + 1):
            for aorb in range(nocb + 1, nbf + 1):
                slot = (aorb - nocb - 1) * noca + iorb
                if slot == ijlr1:
                    matrix[iorb - 1, aorb - nocb - 1] = vector[ijlr1 - 1] * rs
                elif slot == ijlr2:
                    matrix[iorb - 1, aorb - nocb - 1] = -vector[ijlr1 - 1] * rs
                else:
                    matrix[iorb - 1, aorb - nocb - 1] = vector[slot - 1]
        return matrix

    sij0 = FK.s_ij_of(ctx, np.eye(nbf))
    sab0 = FK.s_ab_of(ctx, np.eye(nbf))
    sia0 = FK.s_ia_of(ctx, np.eye(nbf))
    om = [mol.energies[k + 1] - mol.energies[0] for k in range(nstate)]
    phys = [k for k in range(nij) if k != ijlr2 - 1]

    def ampdir(dxt):
        xp = [x.copy() for x in xt0]
        xm = [x.copy() for x in xt0]
        xp[jstate] = xt0[jstate] + 1.0e-6 * dxt
        xm[jstate] = xt0[jstate] - 1.0e-6 * dxt
        sp = FK.contraction(ctx, sij0, sab0, sia0, xt0, xp)
        sm = FK.contraction(ctx, sij0, sab0, sia0, xt0, xm)
        return (sp[:, jstate] - sm[:, jstate]) / (2.0e-6)

    gmet = np.zeros((nstate, nij))
    for slot in phys:
        direction = np.zeros(nij)
        direction[slot] = 1.0
        gmet[:, slot] = ampdir(unfold_vec(direction))

    xc_pair = {}
    xc_mtg = {}
    for istate in range(nstate):
        if istate == jstate:
            continue
        xstate = xf[:, jstate]

        def projected_op(vector):
            full = np.zeros(nij)
            full[phys] = vector
            full -= xstate * float(np.dot(xstate, full))
            result = om[jstate] * full - matvec(full)
            result -= xstate * float(np.dot(xstate, result))
            return result[phys]

        rhs = gmet[istate].copy()
        rhs -= xstate * float(np.dot(xstate, rhs))
        solution, _ = minres(
            LinearOperator((nij - 1, nij - 1), matvec=projected_op),
            rhs[phys], rtol=1.0e-9, maxiter=3000,
        )
        ytil = np.zeros(nij)
        ytil[phys] = solution
        ytil -= xstate * float(np.dot(xstate, ytil))
        xc_pair[istate] = float(np.dot(ytil, gvec_xc))

        rr = x0_raw.copy().reshape(-1)
        rr[istate * nij:(istate + 1) * nij] = ytil
        mol.data['OQP::td_bvec_mo'] = rr.reshape(x0_raw.shape)
        oqp.mrsf_nac_esum(mol, istate + 1, jstate + 1)
        pij_a = np.array(
            mol.data['OQP::dbg_pij_a'], copy=True
        ).reshape(nbf, nbf).T
        pij_b = np.array(
            mol.data['OQP::dbg_pij_b'], copy=True
        ).reshape(nbf, nbf).T
        mol.data['OQP::nac_dm1_a'] = _pack_sym(pij_a)
        mol.data['OQP::nac_dm1_b'] = _pack_sym(pij_b)
        oqp.mrsf_nac_response(mol)
        vxc_a = _unpack_sym(
            np.array(mol.data['OQP::nac_vxc_a'], copy=True).ravel(), nbf
        )
        vxc_b = _unpack_sym(
            np.array(mol.data['OQP::nac_vxc_b'], copy=True).ravel(), nbf
        )
        gma = cmat0.T @ vxc_a @ cmat0
        gmb = cmat0.T @ vxc_b @ cmat0
        mtg = 2.0 * (
            gma * occ_a[None, :] + gmb * occ_b[None, :]
        )
        xc_mtg[istate] = (
            float(np.sum(mtg * ux)),
            float(np.sum(mtg.T * ux)),
        )
    restore()

    print('===== v20 fold/F-channel audit =====')
    print(f'|w_ref|={np.linalg.norm(w_ref):.8f} '
          f'|w_skel|={np.linalg.norm(w_skel):.8f} '
          f'|stagedC|={np.linalg.norm(staged_c):.8f}')
    print(f'|trueF|={np.linalg.norm(true_f):.8f}')
    print(f'actual-F injection: |g|={np.linalg.norm(gvec_actual):.8f} '
          f'|trueF-g|={np.linalg.norm(true_f-gvec_actual):.8e} '
          f'max={np.max(np.abs(true_f-gvec_actual)):.8e}')
    print(f'dD orbital model: |dDa_fd-model|={np.linalg.norm(dda_fd-dda_model):.8e} '
          f'|dDb_fd-model|={np.linalg.norm(ddb_fd-ddb_model):.8e}')
    print(f'F model: |dFa_actual-G[dDmodel]|={np.linalg.norm(dfa_actual-ga_model):.8e} '
          f'|dFb_actual-G[dDmodel]|={np.linalg.norm(dfb_actual-gb_model):.8e}')
    print(f'F from FD density: |dFa_actual-G[dDfd]|={np.linalg.norm(dfa_actual-ga_fd):.8e} '
          f'|dFb_actual-G[dDfd]|={np.linalg.norm(dfb_actual-gb_fd):.8e}')
    print(f'output model: |trueF-g_model|={np.linalg.norm(true_f-gvec_model):.8e} '
          f'|trueF-g_dDfd|={np.linalg.norm(true_f-gvec_fd):.8e}')
    print(f'MO+DM response: |dFa_actual-Gfull|={np.linalg.norm(dfa_actual-ga_full):.8e} '
          f'|dFb_actual-Gfull|={np.linalg.norm(dfb_actual-gb_full):.8e}')
    print(f'XC correction: |g_xc|={np.linalg.norm(gvec_xc):.8e} '
          f'|v19_residual-g_xc|={np.linalg.norm((true_f-gvec_model)-gvec_xc):.8e} '
          f'|trueF-g_full|={np.linalg.norm(true_f-gvec_full):.8e}')
    print(f'analytic JK+XC: |Ga-Gfull|={np.linalg.norm(ga_analytic-ga_full):.8e} '
          f'|Gb-Gfull|={np.linalg.norm(gb_analytic-gb_full):.8e} '
          f'|trueF-g_analytic|={np.linalg.norm(true_f-gvec_analytic):.8e}')
    print(f'analytic XC-only: |dGa_xc-FD|='
          f'{np.linalg.norm(ga_xc_analytic-(ga_full-ga_model)):.8e} '
          f'|dGb_xc-FD|='
          f'{np.linalg.norm(gb_xc_analytic-(gb_full-gb_model)):.8e}')
    print('pair XC correction at this coordinate: ' + ' '.join(
        f'({istate + 1},{jstate + 1})=y.g {value:+.8e} '
        f'Mtg:U {xc_mtg[istate][0]:+.8e} '
        f'MtgT:U {xc_mtg[istate][1]:+.8e}'
        for istate, value in xc_pair.items()
    ))

    residual = true_f - gvec_model
    for slot in np.argsort(-np.abs(residual))[:12]:
        iorb = slot % noca + 1
        aorb = slot // noca + 1 + nocb
        tag = ' LR1' if slot == ijlr1 - 1 else (
            ' LR2' if slot == ijlr2 - 1 else '')
        print(f'  slot {slot}: (i={iorb},a={aorb}) '
              f'residual={residual[slot]:+.8f} '
              f'actual={true_f[slot]:+.8f} model={gvec_model[slot]:+.8f}{tag}')

    np.savez(inp.replace('.inp', f'_c{c0}_v20.npz'), ux=ux,
             w_ref=w_ref, w_skel=w_skel, staged_c=staged_c,
             true_f=true_f, dfa_actual=dfa_actual, dfb_actual=dfb_actual,
             dda_fd=dda_fd, ddb_fd=ddb_fd, dda_model=dda_model,
             ddb_model=ddb_model, ga_model=ga_model, gb_model=gb_model,
             ga_fd=ga_fd, gb_fd=gb_fd, gvec_actual=gvec_actual,
             ga_full=ga_full, gb_full=gb_full, gvec_model=gvec_model,
             ga_analytic=ga_analytic, gb_analytic=gb_analytic,
             ga_xc_analytic=ga_xc_analytic,
             gb_xc_analytic=gb_xc_analytic,
             gvec_fd=gvec_fd, gvec_full=gvec_full,
             gvec_analytic=gvec_analytic, gvec_xc=gvec_xc)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--worker':
        worker(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
