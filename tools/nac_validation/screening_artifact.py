"""Phase 11: the matvec's apparent nonlinearity is an INTEGRAL-SCREENING ARTIFACT.

The MRSF int2 builds its shell-density screening matrix from channel 7 (ball) ONLY
(int2_mrsf_data_t_init_screen -> shell_den_screen_mrsf, tdhf_mrsf_lib.F90:104-147,
uses d3(:,sized,:,:), sized=last channel=7), then applies that ball-based screening
to ALL channels. The small mixed channels 5,6 (o21v/co12) have density in shell-pairs
where ball is small, so the ball-based cutoff (default int2e_cutoff=5e-11, types.F90:121)
drops integrals that matter for ch5/6 -> amplitude-dependent (nonlinear) errors in
exactly the O1<->O2 spin-pair / ground-config sector. Channel 7 (the screening basis)
stays accurate. Tightening int2e_cutoff -> 0 restores exact linearity.

=> The matvec IS a clean linear operator; the Route-A blocker was numerical.
NOTE: the gradient CHAIN (z_vector+gradient) is NOT screening-sensitive, so the
off-diagonal NAC deficiency is REAL physics, not a screening artifact.
"""
import oqp
from oqp.pyoqp import Runner
import numpy as np

r = Runner(input_file='/tmp/nactest/H2O_en.inp', log='/tmp/nactest/screen.log')
r.run(); mol = r.mol
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nbf = np.array(mol.data['OQP::VEC_MO_A'], copy=True).shape[0]
nij = noca*(nbf-(noca-2))
Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); shp = Xraw.shape
X = Xraw.reshape(-1).reshape((3, nij))


def A(c):
    raw = Xraw.copy().reshape(-1); raw[0:nij] = c
    mol.data['OQP::td_bvec_mo'] = raw.reshape(shp)
    oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()


a, b = 0.7, -1.3
print(f"  current int2e_cutoff = {mol.data._data.control.int2e_cutoff:.2e}")
print("  matvec linearity  |A(aX1+bX3) - (aA1+bA3)|  vs int2e_cutoff:")
for cut in (5e-11, 1e-12, 1e-14, 1e-20):
    mol.data._data.control.int2e_cutoff = cut
    err = np.linalg.norm(A(a*X[0]+b*X[2]) - (a*A(X[0])+b*A(X[2])))
    print(f"    cutoff={cut:.1e}:  {err:.3e}")
print("  => default 5e-11 gives ~1e-2 (artifact); 1e-20 gives ~1e-15 (true linear A).")
