"""Same-process ROHF/ROKS forward-response gate.

This is diagnostic only.  It compares the resident Fortran nuclear CPHF
solution with a central derivative of independently re-converged orbitals in
the same process and orbital gauge.  Production NAC uses one adjoint Z-vector,
never this 3N forward solve.
"""
import os
import sys

import numpy as np


def main():
    import oqp
    import oqp.library
    from oqp.library.single_point import SinglePoint
    from oqp.pyoqp import Runner

    inp = sys.argv[1]
    step = float(sys.argv[2]) if len(sys.argv) > 2 else 2.5e-4
    output_npz = sys.argv[3] if len(sys.argv) > 3 else None
    os.environ['NAC_DUMP_ROHF_RESPONSE'] = '1'

    runner = Runner(input_file=inp, log=inp.replace('.inp', '_rohf_gate.log'))
    runner.run()
    mol = runner.mol
    natom = mol.data['natom']
    ncoord = 3 * natom
    nbf = np.array(mol.data['OQP::VEC_MO_A'], copy=True).shape[0]
    noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
    nocb = int(np.asarray(mol.data['nelec_B']).ravel()[0])

    w0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    saved = {
        key: np.array(mol.data[key], copy=True)
        for key in (
            'OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
            'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A',
            'OQP::E_MO_B',
        )
    }

    mol.save_data()
    mol.config['guess']['type'] = 'json'
    mol.config['guess']['file'] = mol.log.replace('.log', '.json')
    mol.config['guess']['continue_geom'] = False

    def displaced_orbitals(coords):
        mol.update_system(coords)
        oqp.library.ints_1e(mol)
        oqp.library.guess(mol)
        SinglePoint(mol).energy()
        mol.data['OQP::xyz_old'] = xyz0.reshape((3, -1))
        mol.data['OQP::VEC_MO_A_old'] = w0
        mol.data['OQP::VEC_MO_B_old'] = wb0
        mol.data['OQP::E_MO_A_old'] = e0a
        mol.data['OQP::E_MO_B_old'] = e0b
        oqp.get_structures_ao_overlap(mol)
        overlap_mo = np.array(
            mol.data['OQP::overlap_mo_non_orthogonal'], copy=True
        ).reshape(nbf, nbf).T
        phase = np.sign(np.diag(overlap_mo))
        phase[phase == 0] = 1.0
        overlap_mo *= phase[None, :]
        overlap_ao = np.array(
            mol.data['OQP::overlap_ao_non_orthogonal'], copy=True
        ).reshape(nbf, nbf).T
        overlap_fixed = w0 @ overlap_ao @ w0.T
        return overlap_mo, overlap_fixed

    u_fd = np.zeros((ncoord, nbf, nbf))
    for coord in range(ncoord):
        plus = xyz0.copy()
        minus = xyz0.copy()
        plus[coord] += step
        minus[coord] -= step
        m_plus, s_plus = displaced_orbitals(plus)
        m_minus, s_minus = displaced_orbitals(minus)
        u_fd[coord] = (
            (m_plus - m_minus) - (s_plus - s_minus)
        ) / (2.0 * step)

    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    for key, value in saved.items():
        mol.data[key] = value

    oqp.hf_hessian(mol)

    def rotation_by_cartesian(tag):
        raw = np.array(mol.data[tag], copy=True)
        nrotation = raw.size // ncoord
        return raw.reshape(ncoord, nrotation).T

    def pack_u(matrix):
        packed = []
        for i in range(nocb, noca):
            for j in range(nocb):
                packed.append(matrix[i, j])
        for a in range(noca, nbf):
            for j in range(nocb):
                packed.append(matrix[a, j])
        for a in range(noca, nbf):
            for i in range(nocb, noca):
                packed.append(matrix[a, i])
        return np.asarray(packed)

    u_native = rotation_by_cartesian('OQP::nac_rohf_uvec')
    b_hf = rotation_by_cartesian('OQP::nac_rohf_bvec_hf_jk_pulay')
    b_full = rotation_by_cartesian('OQP::nac_rohf_bvec_full')
    u_ref = np.column_stack([pack_u(u_fd[c]) for c in range(ncoord)])
    nds = (noca - nocb) * nocb
    ndv = (nbf - noca) * nocb
    blocks = (
        ('doc-socc', slice(0, nds)),
        ('doc-virt', slice(nds, nds + ndv)),
        ('socc-virt', slice(nds + ndv, None)),
    )
    print(f'ROHF response gate: step={step:.3e}')
    for label, block in blocks:
        delta = u_native[block] - u_ref[block]
        print(
            f'{label:9s} max={np.max(np.abs(delta)):.9e} '
            f'rms={np.sqrt(np.mean(delta**2)):.9e}'
        )
    print(f'all       max={np.max(np.abs(u_native-u_ref)):.9e}')
    if output_npz is not None:
        np.savez(
            output_npz,
            u_native=u_native,
            u_reference=u_ref,
            b_hf=b_hf,
            b_full=b_full,
            reference_mo=w0,
            noca=noca,
            nocb=nocb,
            step=step,
            scf_conv=float(mol.config['scf']['conv']),
        )
        print(f'saved {os.path.abspath(output_npz)}')


if __name__ == '__main__':
    main()
