"""Phase 11 decisive experiment: decompose the NUCLEAR cross deficiency
D(I,J) = coded_cross - numeric_cross  into its gradient-density constituents.

Gfull(Z) = grad with all td tags - grad with td tags zeroed (= the SCF ref part).
The "all td tags" gradient is built from:
  - OQP::td_p              : the RELAXED (z-vector) difference density  -> 2e + 1e gradient
  - OQP::WAO               : the energy-weighted Lagrangian             -> overlap (Pulay) gradient
  - OQP::td_mrsf_density   : the 7 MRSF spin-channel densities          -> 2e gradient
  - OQP::td_abxc           : the abxc piece

We re-run Gfull but zero ONE tag at a time (keep the rest) to see which density's
cross carries the (1,2) O-overshoot and which carries the (1,3) H-flip.

This works directly in the ORACLE nuclear space (no seam, no rotation projection).
"""
import math
import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint
import numpy as np

INP = '/tmp/nactest/H2O_tight_dx0.001.inp'
r = Runner(input_file=INP, log='/tmp/nactest/p11gd.log')
r.run(); mol = r.mol
nstate = 3
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True); nbf = C0raw.shape[0]
nvirb = nbf - nocb; nij = noca * nvirb
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
X = X0_raw.reshape(-1).reshape((nstate, nij))
RS = 1 / math.sqrt(2.0)

TD_TAGS = ['OQP::td_p', 'OQP::WAO', 'OQP::td_abxc', 'OQP::td_mrsf_density']


def Gfull(col, zero_tags_extra=()):
    """coded dOmega(Z) = grad(with td) - grad(td zeroed).
    zero_tags_extra: additionally zero these tags in the 'with td' branch so we
    see the gradient WITHOUT that density piece (isolating its contribution by
    differencing against the full)."""
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(X0_raw.shape)
    mol.data.set_tdhf_target(1)
    oqp.tdhf_mrsf_z_vector(mol)
    # snapshot the freshly built td tags
    snap = {t: np.array(mol.data[t], copy=True) for t in TD_TAGS}
    # zero the requested extra tags before gradient
    for t in zero_tags_extra:
        mol.data[t] = np.zeros_like(snap[t])
    oqp.tdhf_mrsf_gradient(mol)
    gZ = mol.get_grad().reshape(-1).copy()
    # restore all td tags then zero ALL -> SCF ref grad
    for t in TD_TAGS:
        mol.data[t] = np.zeros_like(snap[t])
    oqp.tdhf_mrsf_gradient(mol)
    gS = mol.get_grad().reshape(-1).copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    return gZ - gS


def cross(GZ, GI, GJ):
    return GZ - 0.5 * (GI + GJ)


PAIRS = [(0, 1), (0, 2), (1, 2)]
Zs = {(I, J): (X[I] + X[J]) / math.sqrt(2) for (I, J) in PAIRS}

# numeric oracle from qz_nogate
qz = np.load('/tmp/nactest/qz_nogate.npz')
def numcross(I, J):
    return (qz[f'Z{I+1}{J+1}_numeric']
            - 0.5 * (qz[f'X{I+1}_numeric'] + qz[f'X{J+1}_numeric']))


def show(tag, D):
    return np.array2string(D, precision=4, suppress_small=True)


# baseline full coded cross + each "drop one tag" coded cross
variants = [('FULL', ())] + [(f'no_{t.split("::")[1]}', (t,)) for t in TD_TAGS]
print("Decomposing coded cross by dropping each td density tag\n")
for (I, J) in PAIRS:
    GZ = {nm: Gfull(Zs[(I, J)], ex) for nm, ex in variants}
    GI = {nm: Gfull(X[I], ex) for nm, ex in variants}
    GJ = {nm: Gfull(X[J], ex) for nm, ex in variants}
    cfull = cross(GZ['FULL'], GI['FULL'], GJ['FULL'])
    cn = numcross(I, J)
    D = cfull - cn
    print(f"===== pair ({I+1},{J+1}) =====")
    print(f"  coded_FULL = {show('c', cfull)}")
    print(f"  numeric    = {show('n', cn)}")
    print(f"  D=full-num = {show('D', D)}  |D|={np.linalg.norm(D):.5f}")
    print("  contribution of each dropped piece (coded_FULL - coded_no_TAG = that tag's cross):")
    for nm, _ in variants[1:]:
        contrib = cfull - cross(GZ[nm], GI[nm], GJ[nm])
        # how much of D does removing this piece kill?
        dot = np.dot(contrib, D) / (np.linalg.norm(D)**2 + 1e-30)
        print(f"    {nm:18}: contrib={show(nm, contrib)}")
        print(f"    {'':18}  |contrib|={np.linalg.norm(contrib):.5f}  proj_onto_D={dot:+.4f}")
    print()
