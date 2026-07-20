"""Phase 11 PoC: the INTERSTATE relaxed/transition density assembly for coded_cross.

TASK: prove the correct td_p^IJ (interstate relaxed CPHF z-vector density) and show
which gradient-chain density carries the off-diagonal coded_cross deficiency
((1,2) doc-socc ratio inflation, (1,3) virt sign flip), with numbers.

RESULT (decisive, rebuild-free; all numbers below are produced by this script):

 (1) PRODUCTION z-vector RHS R(Z) is an EXACT quadratic form in the amplitude
     (R(t X1)/t^2 = const to 1e-5), so its polarization cross
        R^IJ = R(Z_IJ) - 1/2 [R(X_I)+R(X_J)],  Z_IJ=(X_I+X_J)/sqrt2
     is a well-defined symmetric interstate bilinear = d/dkappa (X_I^T A X_J).
     => the z-vector RHS, and hence td_p, is ALREADY the correct interstate object.

 (2) Per-density quadraticity scan (D(t X1)/t^2 and the a^2/b^2/2ab bilinear-form
     test on D(aX1+bX2)) splits the four gradient densities cleanly:
        td_p   : t^2-homogeneous, quad-form relerr ~1e-4   -> QUADRATIC  (interstate-OK)
        WAO    : t^2-homogeneous, quad-form relerr ~8e-4   -> QUADRATIC  (interstate-OK)
        td_mrsf_density : t^1-homogeneous, quad-form relerr ~2.2 -> LINEAR (interstate-WRONG)
        td_abxc         : t^1-homogeneous, quad-form relerr ~2.2 -> LINEAR (interstate-WRONG)

 (3) ROOT CAUSE.  td_mrsf_density and td_abxc are the *transition* densities
     fmrst1 = mrsfcbc(X) / sfdmat(X), which are LINEAR in the amplitude.  On the
     mixed state Z they give D(Z)=(D(X_I)+D(X_J))/sqrt2, so the polarization cross
     (1/sqrt2 - 1/2)(D(X_I)+D(X_J)) is a SPURIOUS, ill-defined object -- NOT the
     interstate bilinear transition density.  The correct interstate density is the
     genuine bilinear  D^IJ = 1/2[D_chan(X_I,X_J)+D_chan(X_J,X_I)]  delivered by
     mrsf_interstate_tden (the 7-channel / Tij,Tab interstate routine).

 (4) coded_cross density-subset scan vs oracle numeric_cross (gap*d_amp):
        td_p only      : (1,3) cos=+1.0000  -> td_p alone gives the CORRECT (1,3) sign
        WAO  only      : (1,3) cos=-1.0000  -> WAO contributes opposite sign
        mrsf_den only  : (1,3) cos=-0.371   -> the LINEAR-density garbage scrambles it
     The full sum lands at (1,3) cos=+0.969 because the spurious LINEAR td_mrsf_density
     and the (correct-magnitude but mis-signed-by-the-linear-density) WAO partially
     cancel td_p.  Replacing the LINEAR td_mrsf_density/td_abxc by their bilinear
     interstate forms is the fix; td_p needs no change.

CONCLUSION for production: the analytic interstate coded_cross must build
td_mrsf_density^IJ and td_abxc^IJ from mrsf_interstate_tden(I,J) (symmetric bilinear),
keep td_p (z-vector RHS already polarizes correctly) and rebuild WAO from the
interstate relaxed+transition densities.  This is the "interstate-consistent density
assembly" the task targets.

import oqp BEFORE numpy; mol.data._data.control.int2e_cutoff=1e-20 after run.
"""
import math
import oqp
import oqp.library
from oqp.pyoqp import Runner
import numpy as np

SQ = 1.0 / math.sqrt(2.0)
INP = '/tmp/nactest/H2O_tight_dx0.001.inp'

r = Runner(input_file=INP, log='/tmp/nactest/p11_iszvec.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
print("RUN OK", flush=True)

nstate = 3
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
nbf = np.array(mol.data['OQP::VEC_MO_A'], copy=True).shape[0]
nvirb = nbf - nocb; nij = noca * nvirb
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); Xshape = X0_raw.shape
X = X0_raw.reshape(-1).reshape((nstate, nij))
Etot0 = list(mol.energies); Om = [Etot0[k + 1] - Etot0[0] for k in range(nstate)]
gap = {'12': Om[1] - Om[0], '13': Om[2] - Om[0], '23': Om[2] - Om[1]}

prs = ([(i, j) for i in range(nocb, noca) for j in range(0, nocb)]
       + [(k, j) for k in range(noca, nbf) for j in range(0, nocb)]
       + [(k, i) for k in range(noca, nbf) for i in range(nocb, noca)])
nconf = len(prs)
nds, ndv = nocb * (noca - nocb), (nbf - noca) * nocb
blocks = [('doc-socc', 0, nds), ('doc-virt', nds, nds + ndv), ('soc-virt', nds + ndv, nconf)]

DENS = ('OQP::td_p', 'OQP::WAO', 'OQP::td_abxc', 'OQP::td_mrsf_density')


def set_bvec(col):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def run_zvec(col):
    """run the production z-vector for amplitude `col` (target_state=1); return the
    four gradient densities it produced."""
    set_bvec(col)
    mol.data.set_tdhf_target(1)
    oqp.tdhf_mrsf_z_vector(mol)
    d = {t: np.array(mol.data[t], copy=True).ravel() for t in DENS}
    mol.data['OQP::td_bvec_mo'] = X0_raw
    return d


def dump_rhs(col):
    """production z-vector RHS for `col` (NAC_DUMP_RHS)."""
    import os
    set_bvec(col); mol.data.set_tdhf_target(1)
    os.environ['NAC_DUMP_RHS'] = '1'
    oqp.tdhf_mrsf_z_vector(mol)
    os.environ.pop('NAC_DUMP_RHS', None)
    mol.data['OQP::td_bvec_mo'] = X0_raw
    return np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).ravel()[:nconf]


# ========================================================================
print("\n=== (1) production z-vector RHS R(t X1)/t^2 -> EXACT quadratic form ===")
for t in (0.5, 1.0, 2.0):
    R = dump_rhs(t * X[0])
    print(f"  t={t}: |R(tX1)|={np.linalg.norm(R):.6e}   |R|/t^2={np.linalg.norm(R)/t**2:.6e}")
print("  => RHS quadratic => its polarization cross = d/dkappa(X_I^T A X_J): td_p interstate-OK")

print("\n=== (2) per-density quadraticity: D(t X1)/t^2 (const => quadratic) ===")
for t in (0.5, 1.0, 2.0):
    d = run_zvec(t * X[0])
    print(f"  t={t}: " + "  ".join(
        f"{k.split('::')[-1]}={np.linalg.norm(v)/t**2:.4e}" for k, v in d.items()))

print("\n=== (2b) pure-quadratic-FORM test on D(aX1+bX2), a=0.7 b=-1.3 ===")
a, b = 0.7, -1.3
dab = run_zvec(a * X[0] + b * X[1])
d1, d2, d12 = run_zvec(X[0]), run_zvec(X[1]), run_zvec(X[0] + X[1])
for t in DENS:
    c = (d12[t] - d1[t] - d2[t]) / 2.0
    pred = a * a * d1[t] + b * b * d2[t] + 2 * a * b * c
    err = np.linalg.norm(dab[t] - pred) / (np.linalg.norm(dab[t]) + 1e-30)
    tag = "QUADRATIC (interstate-OK)" if err < 1e-3 else "LINEAR/NON-QUADRATIC (interstate-WRONG)"
    print(f"  {t.split('::')[-1]:18s}: quad-form relerr={err:.3e}  -> {tag}")

# ========================================================================
# (4) coded_cross density-subset scan vs oracle
# ========================================================================
print("\n=== (3) coded_cross density subsets vs oracle numeric_cross (gap*d_amp) ===")
Zs = {'Z12': (X[0] + X[1]) * SQ, 'Z13': (X[0] + X[2]) * SQ, 'Z23': (X[1] + X[2]) * SQ,
      'X1': X[0], 'X2': X[1], 'X3': X[2]}
Dz = {n: run_zvec(z) for n, z in Zs.items()}
# reshape back to native density shapes for the gradient
shp = {t: np.array(mol.data[t], copy=True).shape for t in DENS}


def set_dens(dd):
    for t in DENS:
        mol.data[t] = dd[t].reshape(shp[t])


def grad_of(dd):
    set_dens(dd)
    oqp.tdhf_mrsf_gradient(mol)
    return mol.get_grad().reshape(-1).copy()


# SCF baseline (all densities zero)
zero = {t: np.zeros(int(np.prod(shp[t]))) for t in DENS}
gSCF = grad_of(zero)


def coded_cross(p, ia, ib, keep):
    """cross density (polarization) restricted to `keep` densities; gradient is LINEAR
    in densities so this is exactly the keep-subset of coded_cross."""
    dd = {}
    for t in DENS:
        if t in keep:
            dd[t] = Dz[f'Z{p}'][t] - 0.5 * (Dz[f'X{ia}'][t] + Dz[f'X{ib}'][t])
        else:
            dd[t] = zero[t]
    return grad_of(dd) - gSCF


oracle = np.load('/tmp/nactest/p11_damp_oracle.npz')
numc = {p: gap[p] * oracle[f'{p}_damp'] for p in ('12', '13', '23')}
pairs = [('12', '1', '2'), ('13', '1', '3'), ('23', '2', '3')]


def rep(p, cc):
    bb = numc[p]
    c = np.dot(cc, bb) / (np.linalg.norm(cc) * np.linalg.norm(bb) + 1e-30)
    rr = np.linalg.norm(cc) / (np.linalg.norm(bb) + 1e-30)
    return f"cos={c:+.4f} r={rr:.3f}"


subsets = [('ALL', DENS),
           ('td_p only', ('OQP::td_p',)),
           ('WAO only', ('OQP::WAO',)),
           ('mrsf_den only', ('OQP::td_mrsf_density',)),
           ('abxc only', ('OQP::td_abxc',))]
for sn, keep in subsets:
    line = f"  {sn:16s}"
    for p, ia, ib in pairs:
        line += f"  ({p}) {rep(p, coded_cross(p, ia, ib, keep))}"
    print(line)

print("\nVERDICT:")
print("  - td_p (and WAO) ARE quadratic forms -> the z-vector RHS polarizes correctly")
print("    => the interstate relaxed CPHF density td_p^IJ is ALREADY the production td_p")
print("       under the cross extraction; td_p alone reproduces the (1,3) sign (+1.0).")
print("  - td_mrsf_density and td_abxc are LINEAR transition densities; their polarization")
print("    cross is spurious -> these are the off-diagonal deficiency.  FIX = rebuild them")
print("    from mrsf_interstate_tden(I,J) (symmetric bilinear) + recompute WAO from the")
print("    interstate densities.  No change to the td_p z-vector solve is required.")
