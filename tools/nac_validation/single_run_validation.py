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

