#!/usr/bin/env python3
"""pyscf-FREE shells.txt / aux.txt generator for the DF-tensor builder.

Uses basis_set_exchange for the basis data and reimplements the standard
primitive+contraction normalization (verified to reproduce the June builder
inputs, which equal pyscf's _env values, to ~1e-9).  Format per shell:
  l  x y z  nprim
  exp coeff        (nprim lines; coeff = normalized-primitive contraction coeff)
Coordinates in Bohr.  General contractions are segmented.
Usage: make_shells.py <nwat> <shells_out> <aux_out> [--grid AxBxC]
"""
import sys, itertools, math
import numpy as np
import basis_set_exchange as bse

ANG2BOHR = 1.8897259886
Z = {"O": 8, "H": 1}

def water_cluster(n, grid=(4, 4, 4)):
    r=0.9572; th=np.deg2rad(104.52); a=2.8; at=[]; c=0
    for i,j,k in itertools.product(range(grid[0]),range(grid[1]),range(grid[2])):
        if c>=n: break
        o=(i*a,j*a,k*a)
        at+=[("O",o),("H",(o[0]+r,o[1],o[2])),("H",(o[0]+r*np.cos(th),o[1]+r*np.sin(th),o[2]))]; c+=1
    return at

def gto_norm(l, a):
    """Primitive normalization (same convention as pyscf CINTgto_norm)."""
    return math.sqrt(2**(2*l+3) * math.factorial(l+1) * (2*a)**(l+1.5)
                     / (math.factorial(2*l+2) * math.sqrt(math.pi)))

def segmented_shells(basis_name, elem):
    """[(l, exps, coeffs_normalized)] with general contractions segmented and
    each contracted function normalized (radial overlap = 1)."""
    data = bse.get_basis(basis_name, elements=[elem])
    eldata = list(data["elements"].values())[0]
    out = []
    for shell in eldata["electron_shells"]:
        exps = np.array([float(x) for x in shell["exponents"]])
        for l, coefrow in zip(shell["angular_momentum"] * len(shell["coefficients"]),
                              shell["coefficients"]):
            c = np.array([float(x) for x in coefrow])
            m = np.abs(c) > 0
            e, cr = exps[m], c[m]
            cn = cr * np.array([gto_norm(l, a) for a in e])   # fold primitive norm
            # normalize the contracted function: <f|f> = sum_ij cn_i cn_j S_ij,
            # S_ij for normalized primitives = (2 sqrt(ai aj)/(ai+aj))^(l+1.5)
            ee = np.add.outer(e, e); gg = 2*np.sqrt(np.multiply.outer(e, e))/ee
            s = float(cn @ (gg**(l+1.5)) @ cn)
            out.append((l, e, cn/np.sqrt(s)))
    return out

def write_shells(path, atoms, basis_name):
    cache = {}
    recs = []
    for sym, xyz in atoms:
        if sym not in cache: cache[sym] = segmented_shells(basis_name, Z[sym])
        for l, e, c in cache[sym]:
            recs.append((l, tuple(q*ANG2BOHR for q in xyz), e, c))
    with open(path, "w") as f:
        f.write(f"{len(recs)}\n")
        for l, (x, y, z), e, c in recs:
            f.write(f"{l} {x:.12f} {y:.12f} {z:.12f} {len(e)}\n")
            for a, cc in zip(e, c):
                f.write(f"{a:.10g} {cc:.12g}\n")
    return len(recs)

if __name__ == "__main__":
    n = int(sys.argv[1]); shells_out = sys.argv[2]; aux_out = sys.argv[3]
    grid = (4,4,4)
    if len(sys.argv) > 4 and sys.argv[4].startswith("--grid"):
        grid = tuple(int(t) for t in sys.argv[4].split("=")[1].split("x"))
    atoms = water_cluster(n, grid)
    ns = write_shells(shells_out, atoms, "cc-pvdz")
    na = write_shells(aux_out, atoms, "def2-universal-jkfit")
    print(f"[make_shells] nwat={n} grid={grid}: {ns} orbital shells -> {shells_out}, "
          f"{na} aux shells -> {aux_out} (pyscf-free)")
