#!/usr/bin/env python3
"""pyscf-FREE OpenQP-frame prep for the native DF builder.

Runs OpenQP once (its own basis, its own AO frame), dumps the CARTESIAN orbital
shells from mol.data.get_basis() and the auxiliary shells from basis_set_exchange
(def2-universal-jkfit), plus the frozen transform operands the builder loads.
The 2-center metric is NOT computed here -- build_df computes it natively
(ROUTEC_NATIVE_METRIC).  No pyscf anywhere.

Then build the tensor with, e.g.:
  ROUTEC_CART=1 ROUTEC_NATIVE_METRIC=1 ROUTEC_B_SIGMA=1 \
    openqp_gpu_build_df tables.bin shells_c_wN.txt aux_c_wN.txt d.bin 1e-13 <dir> wN 1 B.bin

Usage: native_prep.py NX NY NZ <workdir>
Writes into <workdir>: shells_c_w{N}.txt aux_c_w{N}.txt meta_w{N}.bin
  perm_w{N}.bin dscale_w{N}.bin aperm_w{N}.bin ascale_w{N}.bin cag_w{N}.bin
  linv_w{N}.bin iu_w{N}.bin ju_w{N}.bin   (aperm/ascale/cag/linv are dummies,
  overridden by ROUTEC_NATIVE_METRIC in the builder)."""
import sys, os, math, numpy as np
import basis_set_exchange as bse
from oqp.pyoqp import Runner

# ---- cartesian component convention: routec output -> OpenQP AO frame --------
# routec emits cartesian components in lexicographic order (x power descending,
# then y descending).  OpenQP's order is the GAMESS-style table in constants.F90
# (diagonals first, then off-diagonals).  The per-component scale is the standard
# cartesian GTO factor sqrt((2l-1)!! / (2a-1)!!(2b-1)!!(2c-1)!!) (OpenQP's
# shells_pnrm2).  This was verified for d by matching physical integrals against
# the reference tensor (all components cos=1.0); s/p/d/f/g all follow it.
def _df2(n):                       # double factorial with (-1)!! = (0)!! = 1
    r = 1
    while n > 1: r *= n; n -= 2
    return r
def _routec_cart(l):               # lexicographic: x desc, then y desc
    return [(ax, ay, l-ax-ay) for ax in range(l, -1, -1) for ay in range(l-ax, -1, -1)]
# OpenQP cartesian component order (exponent tuples), from source/constants.F90.
_OQP_CART = {
    0: [(0,0,0)],
    1: [(1,0,0),(0,1,0),(0,0,1)],
    2: [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)],
    3: [(3,0,0),(0,3,0),(0,0,3),(2,1,0),(2,0,1),(1,2,0),(0,2,1),(1,0,2),(0,1,2),(1,1,1)],
    4: [(4,0,0),(0,4,0),(0,0,4),(3,1,0),(3,0,1),(1,3,0),(0,3,1),(1,0,3),(0,1,3),
        (2,2,0),(2,0,2),(0,2,2),(2,1,1),(1,2,1),(1,1,2)],
}
def cart_transform(l):
    """Return (perm, dscale): OpenQP slot m reads routec slot perm[m], scaled by
    dscale[m], to land the AO in OpenQP's exact cartesian frame."""
    if l not in _OQP_CART:
        raise SystemExit(f"cartesian transform for l={l} not tabulated "
                         f"(add its constants.F90 order to _OQP_CART)")
    rc = _routec_cart(l); df2l = _df2(2*l-1)
    perm = [rc.index(t) for t in _OQP_CART[l]]
    dsc  = [math.sqrt(df2l / (_df2(2*a-1)*_df2(2*b-1)*_df2(2*c-1)))
            for (a,b,c) in _OQP_CART[l]]
    return perm, dsc

# ---- geometry (water grid, same as the builder/seam drivers) ----------------
NX, NY, NZ = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); wd = sys.argv[4]
N = NX*NY*NZ; os.makedirs(wd, exist_ok=True)
def water(ox, oy, oz):
    r = 0.9572; th = np.deg2rad(104.52)
    return [(8,(ox,oy,oz)), (1,(ox+r,oy,oz)), (1,(ox+r*np.cos(th),oy+r*np.sin(th),oz))]
atoms = []; aa = 2.8
for i in range(NX):
    for j in range(NY):
        for k in range(NZ): atoms += water(i*aa, j*aa, k*aa)
SYSTEM = "".join(f"\n   {z}   {x:.9f} {y:.9f} {zz:.9f}" for z,(x,y,zz) in atoms)
CFG = {'input':{'system':SYSTEM,'charge':'0','runtype':'energy','basis':'cc-pvdz',
                'method':'hf','functional':'bhhlyp','d4':'False'},
       'guess':{'type':'huckel'},
       'scf':{'type':'rohf','multiplicity':'3','maxit':'5','conv':'1.0e-4',
              'save_molden':'False','incremental':'False'}}
r = Runner(project='nprep', input_dict=CFG, log=f'{wd}/nprep.log', silent=1, usempi=False)
r.run(); mol = r.mol

# ---- OpenQP's own basis (pyscf-free) ----------------------------------------
bas = mol.data.get_basis()
xyz = np.asarray(mol.get_system(), dtype=float).reshape(-1, 3)      # Bohr
angs = np.asarray(bas['angs']); ncontr = np.asarray(bas['ncontr'])
centers = np.asarray(bas['centers']); alpha = np.asarray(bas['alpha']); coef = np.asarray(bas['coef'])

# orbital shells: get_basis coef is already primitive-normalized (matches the
# builder's expected convention); write it directly, cartesian component count.
lines = []; nsh = 0; ip = 0; nao = 0; orb_subs = []
for ish in range(len(angs)):
    l = int(angs[ish]); nc = int(ncontr[ish])
    e = alpha[ip:ip+nc]; c = coef[ip:ip+nc]; ip += nc
    x, y, z = xyz[int(centers[ish])]
    lines.append(f"{l} {x:.12f} {y:.12f} {z:.12f} {nc}")
    for a, cc in zip(e, c): lines.append(f"{a:.12g} {cc:.12g}")
    orb_subs.append((l, int(centers[ish]), list(e), list(c)))   # for the gradient deck
    nsh += 1; nao += (l+1)*(l+2)//2
open(f"{wd}/shells_c_w{N}.txt", "w").write(f"{nsh}\n" + "\n".join(lines) + "\n")

# ---- auxiliary shells (def2-universal-jkfit via bse; spherical) -------------
def gto_norm(l, a):
    return math.sqrt(2**(2*l+3)*math.factorial(l+1)*(2*a)**(l+1.5)
                     /(math.factorial(2*l+2)*math.sqrt(math.pi)))
auxcache = {}; alines = []; ansh = 0; naux = 0
def aux_shells_for(znum):
    if znum in auxcache: return auxcache[znum]
    el = list(bse.get_basis("def2-universal-jkfit", elements=[znum])["elements"].values())[0]
    out = []
    for sh in el["electron_shells"]:
        exps = np.array([float(x) for x in sh["exponents"]])
        for l, crow in zip(sh["angular_momentum"]*len(sh["coefficients"]), sh["coefficients"]):
            c = np.array([float(x) for x in crow]); m = np.abs(c) > 0
            e = exps[m]; cn = c[m]*np.array([gto_norm(l, a) for a in e])   # primitive-normalized
            out.append((l, e, cn))
    auxcache[znum] = out; return out
aux_subs = []
for ia, zn in enumerate([a[0] for a in atoms]):
    x, y, zz = xyz[ia]
    for l, e, cn in aux_shells_for(zn):
        alines.append(f"{l} {x:.12f} {y:.12f} {zz:.12f} {len(e)}")
        for a, cc in zip(e, cn): alines.append(f"{a:.12g} {cc:.12g}")
        aux_subs.append((l, ia, list(e), list(cn)))   # for the gradient deck
        ansh += 1; naux += 2*l+1
open(f"{wd}/aux_c_w{N}.txt", "w").write(f"{ansh}\n" + "\n".join(alines) + "\n")

# ---- transform operands (perm/dscale = routec cart -> OpenQP frame) ----------
npair = nao*(nao+1)//2
perm = np.arange(nao, dtype=np.int32); dscale = np.ones(nao); o = 0
for ish in range(len(angs)):
    l = int(angs[ish]); nc = (l+1)*(l+2)//2
    pm, ds = cart_transform(l)
    for b in range(nc): perm[o+b] = o + pm[b]; dscale[o+b] = ds[b]
    o += nc
assert o == nao
iu, ju = np.tril_indices(nao)

# aperm/ascale/cag/linv are placeholders; ROUTEC_NATIVE_METRIC overrides them.
np.array([nao, naux, npair], dtype=np.int32).tofile(f"{wd}/meta_w{N}.bin")
np.zeros((naux, naux)).tofile(f"{wd}/cag_w{N}.bin")
np.zeros((naux, naux)).tofile(f"{wd}/linv_w{N}.bin")
perm.tofile(f"{wd}/perm_w{N}.bin")
np.ascontiguousarray(dscale, dtype=np.float64).tofile(f"{wd}/dscale_w{N}.bin")
np.arange(naux, dtype=np.int32).tofile(f"{wd}/aperm_w{N}.bin")
np.ones(naux, dtype=np.float64).tofile(f"{wd}/ascale_w{N}.bin")
iu.astype(np.int32).tofile(f"{wd}/iu_w{N}.bin")
ju.astype(np.int32).tofile(f"{wd}/ju_w{N}.bin")

# ---- gradient deck (OQP_ROUTEC_GRAD_INP) so routec_grad2 has the basis --------
# The gradient's derivative integrals need the shells/aux + the routec->OpenQP
# AO map, which the seam does not pass (it hands over only density + xyz). Emit
# them here from the SAME data as the tensor: atom index per shell (coords come
# from the seam's xyz), cartesian AO count, primitive-normalized coefficients,
# and the map = (perm, dscale).
def _subs(tag, subs):
    out = [f"{tag} {len(subs)}"]
    for l, atom, e, c in subs:
        out.append(f"{l} {atom} {len(e)}")
        out += [f"{a:.12g} {cc:.12g}" for a, cc in zip(e, c)]
    return out
# the gradient counts BOTH orbital and aux in cartesian components (ncart),
# unlike the tensor whose aux index is spherical -- so the deck's naux is the
# cartesian aux count.
naux_cart = sum((l+1)*(l+2)//2 for (l, _a, _e, _c) in aux_subs)
deck = [f"routec_grad_inp v1", f"natm {len(atoms)}", f"nbf {nao}", f"naux {naux_cart}"]
deck += _subs("nsub", orb_subs)
deck += _subs("nauxsub", aux_subs)
deck.append("map")
deck += [f"{int(perm[i])} {dscale[i]:.12g}" for i in range(nao)]
open(f"{wd}/grad_w{N}.inp", "w").write("\n".join(deck) + "\n")

print(f"w{N}: nao(cart)={nao} naux={naux} npair={npair} nsh={nsh} ansh={ansh} -> {wd}", flush=True)
