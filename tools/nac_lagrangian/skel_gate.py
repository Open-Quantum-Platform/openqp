"""FROZEN-MO SKELETON FD GATE.

Certifies the skeleton engines against their composite defining object:
with MOs, amplitudes AND density frozen at the reference, displace the
geometry and central-difference the bilinear E_IJ(x) = X_I^T A(x) X_J
(matvec with F(x)[D_ref] and frozen C).  The exact contract is

    dE_IJ/dx |skeleton  ==  amp2e + esum     (wsx/zL/zG/ov are response)

Any mismatch localizes a Gamma-convention bug INSIDE the skeleton pair;
a match proves the six-term scatter lives in the response bookkeeping.

Driver:  python skel_gate.py <input.inp>
Worker:  python skel_gate.py --worker <inp> <ref.npz> <out.npz>
"""
import os
import sys
import copy
import subprocess
import numpy as np

H = float(os.environ.get(
    'NAC_SKEL_H', '1.0e-3'
))                   # displacement in get_system units
NPAR = int(os.environ.get('NAC_SKEL_NPAR', '12'))
WOMP = os.environ.get('NAC_SKEL_OMP', '2')


def worker(inp, ref_npz, out_npz):
    import oqp
    from oqp.pyoqp import Runner
    ref = np.load(ref_npz)
    r = Runner(input_file=inp, log=inp.replace('.inp', '.log'))
    r.run()
    mol = r.mol
    nstate = int(ref['nstate'])
    nij = int(ref['nij'])
    X0_raw = ref['X0_raw']
    mol.data['OQP::VEC_MO_A'] = ref['C_a'].copy()
    mol.data['OQP::VEC_MO_B'] = ref['C_b'].copy()
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    E = np.zeros((nstate, nstate))
    Ax_all = np.zeros((nstate, nij))
    for s in range(nstate):
        rr = X0_raw.copy().reshape(-1)
        rr[0:nij] = Xf[:, s]
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        oqp.mrsf_matvec_apply(mol)
        Ax = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
        Ax_all[s] = Ax[:nij]
        for I in range(nstate):
            E[I, s] = float(np.dot(Xf[:, I], Ax))
    np.savez(out_npz, E=E, Ax=Ax_all)


def main():
    import oqp
    from oqp.pyoqp import Runner
    from oqp.utils.file_utils import write_config, write_xyz
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp = sys.argv[1]
    base = inp.replace('.inp', '')
    wdir = base + '_skel'
    os.makedirs(wdir, exist_ok=True)
    r = Runner(input_file=inp, log=base + '_sk.log')
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    E0 = list(mol.energies)
    Om = [E0[k + 1] - E0[0] for k in range(nstate)]

    ctx = FK.build_context(mol)
    nij = ctx['nij']
    X0_raw = ctx['X0_raw']
    C_a = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    C_b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)

    # skeleton engines at reference
    oqp.mrsf_nac_amp(mol)
    amp2e = np.array(mol.data['OQP::nac_amp'], copy=True
                     ).reshape((nstate, nstate, natom, 3))
    esum_v, wsx_v = {}, {}
    for i in range(1, nstate + 1):
        for j in range(1, nstate + 1):
            if i == j:
                continue
            oqp.mrsf_nac_esum(mol, i, j)
            esum_v[(i, j)] = np.array(mol.data['OQP::nac_esum'],
                                      copy=True).reshape(-1)
            wsx_v[(i, j)] = np.array(mol.data['OQP::nac_wsx'],
                                     copy=True).reshape(-1)

    # sanity: frozen-path bilinear at the reference must be diag(omega)
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    Eref = np.zeros((nstate, nstate))
    for s in range(nstate):
        rr = X0_raw.copy().reshape(-1)
        rr[0:nij] = Xf[:, s]
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        oqp.mrsf_matvec_apply(mol)
        Ax = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
        for I in range(nstate):
            Eref[I, s] = float(np.dot(Xf[:, I], Ax))
    mol.data['OQP::td_bvec_mo'] = X0_raw
    print('SANITY E_ref (should be diag(omega)):')
    print(np.array2string(Eref, precision=8, suppress_small=True))
    print('omega:', np.round(Om, 8), flush=True)

    ref_npz = os.path.join(wdir, 'ref.npz')
    np.savez(ref_npz, C_a=C_a, C_b=C_b, X0_raw=X0_raw,
             nstate=nstate, nij=nij)

    # displaced inputs: 1-iteration SCF from the reference-density json guess
    mol.save_data()
    guess_file = mol.log.replace('.log', '.json')
    origin = mol.get_system()
    ncoord = len(origin)
    atoms = mol.get_atoms()
    config = copy.deepcopy(mol.config)
    config['input']['runtype'] = 'energy'
    config['guess']['type'] = 'json'
    config['guess']['file'] = guess_file
    config['guess']['file2'] = guess_file
    config['guess']['continue_geom'] = 'false'
    config['scf']['maxit'] = 1
    config['scf']['conv'] = 0.5
    config['tests']['exception'] = 'false'

    jobs = []
    for idx in range(2 * ncoord):
        c = idx % ncoord
        sgn = 1.0 if idx < ncoord else -1.0
        coord = origin.copy()
        coord[c] += sgn * H
        xyz = os.path.join(wdir, f'p{idx}.xyz')
        pinp = os.path.join(wdir, f'p{idx}.inp')
        out = os.path.join(wdir, f'p{idx}.npz')
        cfg = copy.deepcopy(config)
        cfg['input']['system'] = xyz
        with open(xyz, 'w') as f:
            f.write(write_xyz(atoms, coord, [idx]))
        input_file, _ = write_config(cfg)
        with open(pinp, 'w') as f:
            f.write(input_file)
        jobs.append((pinp, out))

    print(f'launching {len(jobs)} frozen-MO displaced matvec workers...',
          flush=True)
    env = dict(os.environ, OMP_NUM_THREADS=WOMP)
    procs = []
    for pinp, out in jobs:
        while len([p for p in procs if p.poll() is None]) >= NPAR:
            import time
            time.sleep(0.5)
        procs.append(subprocess.Popen(
            [sys.executable, os.path.abspath(__file__), '--worker',
             pinp, ref_npz, out],
            env=env, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT))
    for p in procs:
        p.wait()

    FD = np.zeros((nstate, nstate, ncoord))
    bad = []
    for c in range(ncoord):
        fp = os.path.join(wdir, f'p{c}.npz')
        fm = os.path.join(wdir, f'p{c + ncoord}.npz')
        if not (os.path.exists(fp) and os.path.exists(fm)):
            bad.append(c)
            continue
        FD[:, :, c] = (np.load(fp)['E'] - np.load(fm)['E']) / (2 * H)
    if bad:
        print('MISSING worker outputs for coords:', bad, flush=True)

    print('\n===== FROZEN-MO SKELETON GATE =====')
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            fd = FD[I, J]
            sk = amp2e[J, I].reshape(-1) + esum_v[(I + 1, J + 1)]
            skw = sk + wsx_v[(I + 1, J + 1)]
            c1 = float(np.dot(fd, sk)) / (np.linalg.norm(fd)
                                          * np.linalg.norm(sk) + 1e-300)
            c2 = float(np.dot(fd, skw)) / (np.linalg.norm(fd)
                                           * np.linalg.norm(skw) + 1e-300)
            print(f'pair ({I+1},{J+1}): |FD|={np.linalg.norm(fd):.6f} '
                  f'|amp+esum|={np.linalg.norm(sk):.6f} cos={c1:+.6f} '
                  f'maxdiff={np.abs(fd-sk).max():.3e} | '
                  f'+wsx: cos={c2:+.6f} maxdiff={np.abs(fd-skw).max():.3e}')
    # FD symmetry sanity (A symmetric => E^x symmetric)
    asym = max(np.abs(FD[i, j] - FD[j, i]).max()
               for i in range(nstate) for j in range(nstate) if i != j)
    print(f'FD pair-symmetry max|E_IJ^x - E_JI^x| = {asym:.3e}')
    np.savez(os.path.join(wdir, 'skel_result.npz'), FD=FD,
             amp2e=amp2e, esum=np.array(
                 [[esum_v.get((i + 1, j + 1), np.zeros(ncoord))
                   for j in range(nstate)] for i in range(nstate)]),
             wsx=np.array(
                 [[wsx_v.get((i + 1, j + 1), np.zeros(ncoord))
                   for j in range(nstate)] for i in range(nstate)]))


if __name__ == '__main__':
    if sys.argv[1] == '--worker':
        worker(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        main()
