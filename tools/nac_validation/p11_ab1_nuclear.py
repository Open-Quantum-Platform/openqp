"""Phase 11 decisive nuclear experiment: decompose the coded cross deficiency
D(I,J)=coded-numeric by GATING each z-vector RHS piece (AB1, T, HX) in the ACTUAL
gradient solve (the gates at tdhf_mrsf_z_vector.F90:1744-1755 modify the rhs that
feeds the solve + nuclear gradient, NOT just the dump).

For each pair we compute the coded cross with a gate ON vs OFF; the difference is
that piece's contribution to the nuclear cross. We test which piece's contribution
aligns with D for (1,2) [O-atom symmetric], (1,3) [H-H antisym sign flip], (2,3)[~0].

This works in the ORACLE nuclear space directly. AB1 was 'ruled out' on the diagonal
(scale -0.577) but never tested in the symmetric nuclear CROSS, which is the oracle.
"""
import os
import math
import oqp
import oqp.library
from oqp.pyoqp import Runner
import numpy as np

INP = '/tmp/nactest/H2O_tight_dx0.001.inp'
r = Runner(input_file=INP, log='/tmp/nactest/p11ab1.log'); r.run(); mol = r.mol
nstate = 3
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
nbf = np.array(mol.data['OQP::VEC_MO_A'], copy=True).shape[0]
nvirb = nbf - nocb; nij = noca * nvirb
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
X = X0_raw.reshape(-1).reshape((nstate, nij))
RS = 1 / math.sqrt(2.0)
TD_TAGS = ['OQP::td_p', 'OQP::WAO', 'OQP::td_abxc', 'OQP::td_mrsf_density']
ALLG = ('NAC_ZERO_2FT', 'NAC_ZERO_AB1', 'NAC_ZERO_HX', 'NAC_ZERO_T', 'NAC_ZERO_SP')


def Gfull(col, gates=()):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(X0_raw.shape)
    mol.data.set_tdhf_target(1)
    for k in ALLG:
        os.environ.pop(k, None)
    for k in gates:
        os.environ[k] = '1'
    oqp.tdhf_mrsf_z_vector(mol)
    for k in ALLG:
        os.environ.pop(k, None)
    snap = {t: np.array(mol.data[t], copy=True) for t in TD_TAGS}
    oqp.tdhf_mrsf_gradient(mol)
    gZ = mol.get_grad().reshape(-1).copy()
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
qz = np.load('/tmp/nactest/qz_nogate.npz')


def numcross(I, J):
    return (qz[f'Z{I+1}{J+1}_numeric']
            - 0.5 * (qz[f'X{I+1}_numeric'] + qz[f'X{J+1}_numeric']))


def sh(D):
    return np.array2string(D, precision=4, suppress_small=True)


def u(v):
    n = np.linalg.norm(v); return v/n if n > 1e-14 else v


# gate sets to probe (each modifies the actual solve)
GATESETS = [('FULL', ()), ('noAB1', ('NAC_ZERO_AB1',)),
            ('noHX', ('NAC_ZERO_HX',)), ('noT', ('NAC_ZERO_T',)),
            ('no2FT', ('NAC_ZERO_2FT',)), ('noSP', ('NAC_ZERO_SP',))]

print("Probing which z-vector RHS piece carries each pair's nuclear deficiency\n")
for (I, J) in PAIRS:
    G = {nm: {'Z': Gfull(Zs[(I, J)], g), 'I': Gfull(X[I], g), 'J': Gfull(X[J], g)}
         for nm, g in GATESETS}
    cfull = cross(G['FULL']['Z'], G['FULL']['I'], G['FULL']['J'])
    cn = numcross(I, J)
    D = cfull - cn
    print(f"===== pair ({I+1},{J+1}) =====")
    print(f"  coded_FULL = {sh(cfull)}")
    print(f"  numeric    = {sh(cn)}  |num|={np.linalg.norm(cn):.5f}")
    print(f"  D          = {sh(D)}   |D|={np.linalg.norm(D):.5f}")
    for nm, g in GATESETS[1:]:
        cg = cross(G[nm]['Z'], G[nm]['I'], G[nm]['J'])
        piece = cfull - cg   # contribution of the gated-out piece to the coded cross
        # would removing this piece from coded move it toward numeric?
        cfix = cfull - piece  # == cg
        cos_fix = np.dot(u(cfix), u(cn)); rat_fix = np.linalg.norm(cfix)/(np.linalg.norm(cn)+1e-30)
        cos_pD = np.dot(u(piece), u(D))
        print(f"    {nm:7}: |piece|={np.linalg.norm(piece):.5f} cos(piece,D)={cos_pD:+.4f}  "
              f"| if removed: cos(fix,num)={cos_fix:+.4f} ratio={rat_fix:.3f}")
    print()
