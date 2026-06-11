import oqp                       # import oqp (loads its ILP64 LAPACK) BEFORE numpy
from oqp.pyoqp import Runner
from oqp.library.single_point import NAC
import numpy as np

# Proper setup + R0 reference SCF/excitation via the normal run() path
r = Runner(input_file='/tmp/nactest/H2O_energy.inp', log='/tmp/nactest/harness.log')
r.run()                      # runtype=energy: SCF + excitation, saves R0 anchor json
mol = r.mol

nac = NAC(mol)

# numerical NAC (each displacement aligns to the saved R0 anchor)
nacv_n, dcv_n, fn = nac.numerical_nac()

# orbital-overlap term alone, from the SAME R0 amplitudes
oqp.mrsf_nac_overlap(mol)
ovraw = np.array(mol.data['OQP::nac_overlap'], copy=True)
nstate = nac.nstate; natom = nac.natom
ov = np.transpose(ovraw, (1, 2, 0)).reshape((nstate, nstate, natom, 3))

# analytical NAC (CI + ov), same R0 amplitudes
nacv_a, dcv_a, fa = nac.analytical_nac()

i, j = 2, 3
dn = dcv_n[i-1, j-1].reshape(-1)
da = dcv_a[i-1, j-1].reshape(-1)
dov = ov[i-1, j-1].reshape(-1)
dci = da - dov


def cmp(name, v):
    s = np.sign(np.dot(dn, v)) or 1.0
    cos = np.dot(dn, v) / (np.linalg.norm(dn) * np.linalg.norm(v) + 1e-30)
    print(f'{name}: |v|={np.linalg.norm(v):.5f}  cos(num,v)={cos:+.4f}  '
          f'best-sign resid={np.linalg.norm(s*v - dn):.5f}')


print('\n=== pair (2,3), shared R0 anchor ===')
print('numerical |d| =', np.linalg.norm(dn))
cmp('analytical TOTAL', da)
cmp('CI part only    ', dci)
cmp('overlap only    ', dov)
print('\nnumerical d:', np.round(dn, 5))
print('analyt  d:', np.round(da, 5))
print('CI part  :', np.round(dci, 5))
print('overlap  :', np.round(dov, 5))

# align numerical to the analytical sign convention (global per-pair sign)
s = np.sign(np.dot(dn, da)) or 1.0
dn_al = s * dn
# orbital part (must be supplied by frozen-overlap + CPHF) and the CPHF target
d_orb_target = dn_al - dci
d_cphf_target = dn_al - dci - dov
print('\n--- decomposition targets ---')
print('orbital target (dn-dci)      :', np.round(d_orb_target, 5),
      ' |.|=%.4f' % np.linalg.norm(d_orb_target))
print('frozen overlap (coded)       :', np.round(dov, 5),
      ' |.|=%.4f' % np.linalg.norm(dov))
print('CPHF target (dn-dci-dov)      :', np.round(d_cphf_target, 5),
      ' |.|=%.4f' % np.linalg.norm(d_cphf_target))
print('frozen/orbital ratio = %.3f' % (np.linalg.norm(dov)/np.linalg.norm(d_orb_target)))
print('CPHF/orbital   ratio = %.3f' % (np.linalg.norm(d_cphf_target)/np.linalg.norm(d_orb_target)))

# ===== analytical CPHF term via z-vector interchange =====
def grad_now():
    oqp.tdhf_mrsf_gradient(mol)
    return mol.get_grad().reshape((natom, 3)).copy()

def cphf_term(i, j):
    mol.data.set_tdhf_target(i)
    oqp.set_mrsf_nac_cphf(mol, i, j)      # CPHF mode on
    oqp.tdhf_mrsf_z_vector(mol)
    if not mol.mol_energy.Z_Vector_converged:
        oqp.set_mrsf_nac_cphf(mol, 0, 0)
        return None
    gZ = grad_now()                        # SCF + nuc + d^cphf
    # zero the relaxed/transition densities -> pure SCF+nuc gradient
    mol.data['OQP::td_p'] = np.zeros_like(np.array(mol.data['OQP::td_p'], copy=True))
    gS = grad_now()
    oqp.set_mrsf_nac_cphf(mol, 0, 0)      # CPHF mode off
    return (gZ - gS).reshape(-1)

d_cphf_ana = cphf_term(i, j)
print('\n--- analytical CPHF term vs target ---')
if d_cphf_ana is None:
    print('CPHF z-vector did NOT converge')
else:
    t = d_cphf_target
    s = np.sign(np.dot(t, d_cphf_ana)) or 1.0
    cos = np.dot(t, d_cphf_ana)/(np.linalg.norm(t)*np.linalg.norm(d_cphf_ana)+1e-30)
    print('CPHF target :', np.round(t, 5), ' |.|=%.4f' % np.linalg.norm(t))
    print('CPHF analyt :', np.round(d_cphf_ana, 5), ' |.|=%.4f' % np.linalg.norm(d_cphf_ana))
    print('cos(target, analyt) = %+.4f   best-sign resid = %.5f' %
          (cos, np.linalg.norm(s*d_cphf_ana - t)))

