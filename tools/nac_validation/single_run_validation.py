import sys
import oqp                       # import oqp (loads its ILP64 LAPACK) BEFORE numpy
from oqp.pyoqp import Runner
from oqp.library.single_point import NAC
import numpy as np

inp = sys.argv[1] if len(sys.argv) > 1 else '/tmp/nactest/H2O_energy.inp'

# Proper setup + R0 reference SCF/excitation via the normal run() path
r = Runner(input_file=inp, log=inp.replace('.inp', '_harness.log'))
r.run()                      # runtype=energy: SCF + excitation, saves R0 anchor json
mol = r.mol

nac = NAC(mol)
nstate = nac.nstate
natom = nac.natom

# numerical NAC for all configured pairs (each displacement aligns to R0)
nacv_n, dcv_n, fn = nac.numerical_nac()

# orbital-overlap (frozen + skeleton) term, same R0 amplitudes
oqp.mrsf_nac_overlap(mol)
ovraw = np.array(mol.data['OQP::nac_overlap'], copy=True)
# tagarray getter exposes the Fortran column-major (3natom,nstate,nstate)
# buffer with C strides; re-read in Fortran order: ovF[jst-1, ist-1, k-1]
ovF = ovraw.reshape(-1).reshape((nstate, nstate, 3 * natom))
ov = np.transpose(ovF, (1, 0, 2)).reshape((nstate, nstate, natom, 3))

# analytical NAC (CI + ov), same R0 amplitudes
nacv_a, dcv_a, fa = nac.analytical_nac()


def grad_now():
    oqp.tdhf_mrsf_gradient(mol)
    return mol.get_grad().reshape((natom, 3)).copy()


def cphf_term(i, j, block=0):
    mol.data.set_tdhf_target(i)
    oqp.set_mrsf_nac_cphf(mol, i, j)      # CPHF mode on
    oqp.set_mrsf_nac_cphf_block(mol, block)
    oqp.tdhf_mrsf_z_vector(mol)
    if not mol.mol_energy.Z_Vector_converged:
        oqp.set_mrsf_nac_cphf(mol, 0, 0)
        oqp.set_mrsf_nac_cphf_block(mol, 0)
        return None
    gZ = grad_now()                        # SCF + nuc + d^cphf
    # zero the relaxed/transition densities AND the response energy-weighted
    # density W(Z) (wao) -> pure SCF+nuc gradient; the difference then
    # includes BOTH P(Z).F^x and W(Z).S^x, i.e. the complete -Z.B^x.
    mol.data['OQP::td_p'] = np.zeros_like(np.array(mol.data['OQP::td_p'], copy=True))
    mol.data['OQP::WAO'] = np.zeros_like(np.array(mol.data['OQP::WAO'], copy=True))
    gS = grad_now()
    oqp.set_mrsf_nac_cphf(mol, 0, 0)      # CPHF mode off
    oqp.set_mrsf_nac_cphf_block(mol, 0)
    return (gZ - gS).reshape(-1)


def cos(u, v):
    return np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-30)


pairs = [tuple(p) for p in nac.nac_states]
print('\n================= per-pair analysis =================')
for (i, j) in pairs:
    dn = dcv_n[i - 1, j - 1].reshape(-1)
    da = dcv_a[i - 1, j - 1].reshape(-1)
    dov = ov[i - 1, j - 1].reshape(-1)
    dci = da - dov
    sgn = np.sign(np.dot(dn, da)) or 1.0
    dn_al = sgn * dn                       # align numerical to analytical sign
    t = dn_al - dci - dov                  # exact CPHF target
    D = cphf_term(i, j)
    print(f'\n--- pair ({i},{j})  gap={mol.energies[j]-mol.energies[i]:.6f} ---')
    print(f'|d_num|={np.linalg.norm(dn):.5f}  |d_CI|={np.linalg.norm(dci):.5f}  '
          f'|d_ov|={np.linalg.norm(dov):.5f}  |t_cphf|={np.linalg.norm(t):.5f}')
    if D is None:
        print('CPHF z-vector did NOT converge')
        continue
    scale = np.dot(t, D) / np.dot(D, D)
    total = dci + dov + scale * D
    print(f'CPHF: |D|={np.linalg.norm(D):.5f}  cos(t,D)={cos(t, D):+.6f}  '
          f'GLOBAL SCALE s={scale:+.6f}')
    print(f'TOTAL vs num: cos={cos(total, dn_al):+.6f}  '
          f'|tot|={np.linalg.norm(total):.5f}  |num|={np.linalg.norm(dn_al):.5f}  '
          f'resid={np.linalg.norm(total - dn_al):.5f}')
    # per-block split (doc-virt RHS is identically zero at TDA)
    D1 = cphf_term(i, j, block=1)
    D3 = cphf_term(i, j, block=3)
    if D1 is not None and D3 is not None:
        M = np.stack([D1, D3], axis=1)
        s2, *_ = np.linalg.lstsq(M, t, rcond=None)
        fit = M @ s2
        tot2 = dci + dov + fit
        print(f'2-scale: s1(doc-socc)={s2[0]:+.6f} (|D1|={np.linalg.norm(D1):.4f})  '
              f's3(socc-virt)={s2[1]:+.6f} (|D3|={np.linalg.norm(D3):.4f})')
        print(f'         fit resid={np.linalg.norm(fit - t):.5f}  '
              f'TOTAL2 cos={cos(tot2, dn_al):+.6f}')
