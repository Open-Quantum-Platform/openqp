"""Self-contained molecular DFT grid generator (NO pyscf).

Standard Becke molecular quadrature:
  per atom:  radial Treutler-Ahlrichs (Gauss-Chebyshev-2 mapped) x Lebedev angular
  pruning:   NWChem-style (smaller angular orders in the core region)
  partition: original Becke fuzzy cells (k=3 smoothing), Bragg-radii size adjust
Angular nodes/weights come from lebedev_np.py (ported from OpenQP's
source/lebedev.F90 -- Lebedev-Laikov constants; no pyscf anywhere).

Output: grid.bin in the engine's format: int64 npts | coords (npts,3) f64 Bohr
        | weights (npts) f64.

Gate: swap this grid.bin under an otherwise identical prep and require the DFT
energy to match the reference-grid energy to ~1e-6 Ha.
"""
import sys, struct
import numpy as np
from lebedev_np import lebedev

ANG2BOHR = 1.0 / 0.52917721092

# Bragg-Slater radii (Angstrom), H..Ar (pyscf/most codes use the same table)
BRAGG_A = {1:0.35, 2:1.40, 3:1.45, 4:1.05, 5:0.85, 6:0.70, 7:0.65, 8:0.60,
           9:0.50, 10:1.50, 11:1.80, 12:1.50, 13:1.25, 14:1.10, 15:1.00,
           16:1.00, 17:1.00, 18:1.88}
# Treutler-Ahlrichs xi parameters (their Table 1), H..Ar
TA_XI = {1:0.8, 2:0.9, 3:1.8, 4:1.4, 5:1.3, 6:1.1, 7:0.9, 8:0.9, 9:0.9,
         10:0.9, 11:1.4, 12:1.3, 13:1.3, 14:1.2, 15:1.1, 16:1.0, 17:1.0, 18:1.0}
ZNUM = {"H":1, "He":2, "Li":3, "Be":4, "B":5, "C":6, "N":7, "O":8, "F":9,
        "Ne":10, "Na":11, "Mg":12, "Al":13, "Si":14, "P":15, "S":16,
        "Cl":17, "Ar":18}
# radial points per element row (H-He lighter, rest heavier) -- level-3-like
def n_rad(z):
    import os
    f = float(os.environ.get("GRID_FINE","1.0"))
    return int((50 if z <= 2 else 75)*f)


def radial_treutler(n, xi):
    """Treutler-Ahlrichs M4 radial grid on Gauss-Chebyshev-2 nodes.
    Returns r (Bohr) and radial weights INCLUDING the r^2 volume factor."""
    i = np.arange(1, n + 1)
    x = np.cos(i * np.pi / (n + 1))                    # Chebyshev-2 nodes (-1,1)
    wch = np.pi / (n + 1) * np.sin(i * np.pi / (n + 1))**2  # GC2 weights
    ln2 = np.log(2.0)
    a = 1.0
    r = xi / ln2 * (a + x)**0.6 * np.log(2.0 / (1.0 - x))
    # dr/dx of the M4 map
    drdx = xi / ln2 * (a + x)**0.6 * (0.6 * np.log(2.0 / (1.0 - x)) / (a + x)
                                      + 1.0 / (1.0 - x))
    # GC2 quadrature integrates f(x)*sqrt(1-x^2); undo the sqrt factor
    w = wch / np.sqrt(1.0 - x**2) * drdx * r**2
    return r[::-1], w[::-1]                             # ascending r


# NWChem-style pruning: angular order per radial shell by r / bragg_radius.
# Regions (fractions of Bragg radius) -> Lebedev sizes scaled off the max order.
def prune_order(nang_max, r, rb):
    if nang_max <= 50:          # tiny grids: no pruning
        return nang_max
    frac = r / rb
    if   frac < 0.25: return 26
    elif frac < 0.5:  return 50
    elif frac < 1.0:  return 110
    elif frac < 4.5:  return nang_max
    else:             return nang_max


def becke_weights(points, iat, coords, radii, k=3, rcut=15.0, chunk=8192):
    """Original Becke fuzzy-cell weights for `points` belonging to atom iat.
    Fully vectorized over (point, a, b) with a neighbor cutoff: atoms farther
    than rcut (Bohr) from the point batch have cell functions ~1 (toward iat)
    and populations ~0, so they drop out of both product and normalization."""
    ctr = points.mean(0)
    rad = float(np.max(np.linalg.norm(points - ctr, axis=1)))
    d_at = np.linalg.norm(coords - ctr, axis=1)
    loc = np.where(d_at <= rad + rcut)[0]
    if iat not in loc: loc = np.append(loc, iat)
    ii = int(np.where(loc == iat)[0][0])
    C = coords[loc]; R = radii[loc]; nat = len(loc)
    Rab = np.linalg.norm(C[:, None, :] - C[None, :, :], axis=2)
    np.fill_diagonal(Rab, 1.0)
    chi = R[:, None] / R[None, :]
    u = (chi - 1.0) / (chi + 1.0)
    aab = np.clip(u / (u * u - 1.0), -0.5, 0.5)
    try:
        import cupy as xp                      # GPU path when available (~100x)
    except ImportError:
        xp = np
    Cx = xp.asarray(C); Rabx = xp.asarray(Rab); aabx = xp.asarray(aab)
    out = np.empty(len(points))
    idx = xp.arange(nat)
    for p0 in range(0, len(points), chunk):
        pts = xp.asarray(points[p0:p0 + chunk])
        d = xp.linalg.norm(pts[:, None, :] - Cx[None, :, :], axis=2)     # (np,nat)
        mu = (d[:, :, None] - d[:, None, :]) / Rabx[None, :, :]          # (np,a,b)
        mu = mu + aabx[None, :, :] * (1.0 - mu * mu)
        f = mu
        for _ in range(k):
            f = 1.5 * f - 0.5 * f**3
        s = 0.5 * (1.0 - f)
        s[:, idx, idx] = 1.0
        P = s.prod(axis=2)                                               # (np,nat)
        r = P[:, ii] / P.sum(axis=1)
        out[p0:p0 + chunk] = r.get() if xp is not np else r
    return out


def build_grid(atoms_ang, nang_max=302):
    """atoms_ang: list of (symbol, (x,y,z) in ANGSTROM). Returns coords (Bohr), weights."""
    syms = [a[0] for a in atoms_ang]
    coords = np.array([a[1] for a in atoms_ang], float) * ANG2BOHR
    zs = np.array([ZNUM[s] for s in syms])
    radii = np.array([BRAGG_A[z] for z in zs]) * ANG2BOHR
    all_p, all_w = [], []
    leb_cache = {}
    for ia, z in enumerate(zs):
        r, wr = radial_treutler(n_rad(z), TA_XI[z])
        rb = radii[ia]
        # group radial shells by pruned angular order
        orders = np.array([prune_order(nang_max, ri, rb) for ri in r])
        for order in np.unique(orders):
            if order not in leb_cache:
                leb_cache[order] = lebedev(int(order))
            ang_p, ang_w = leb_cache[order]
            sel = orders == order
            rs, ws = r[sel], wr[sel]
            # points: (nrad_sel, nang, 3)
            pts = coords[ia] + rs[:, None, None] * ang_p[None, :, :]
            wts = 4.0 * np.pi * ws[:, None] * ang_w[None, :]
            pts = pts.reshape(-1, 3); wts = wts.reshape(-1)
            bw = becke_weights(pts, ia, coords, radii)
            all_p.append(pts); all_w.append(wts * bw)
    P = np.vstack(all_p); W = np.concatenate(all_w)
    keep = W > 1e-14                                   # drop numerically dead points
    return P[keep], W[keep]


def write_grid_bin(path, coords, weights):
    import os
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        f.write(struct.pack("<q", len(weights)))
        f.write(np.ascontiguousarray(coords, dtype="<f8").tobytes())
        f.write(np.ascontiguousarray(weights, dtype="<f8").tobytes())
    os.replace(tmp, path)                      # atomic: consumers poll for `path`


def water_cluster(n, a=2.8, r=0.9572, th=104.52):
    import itertools
    thr = np.deg2rad(th); at = []; c = 0
    for i, j, k in itertools.product(range(4), repeat=3):
        if c >= n: break
        o = np.array([i * a, j * a, k * a], float)
        at += [("O", tuple(o)), ("H", (o[0] + r, o[1], o[2])),
               ("H", (o[0] + r * np.cos(thr), o[1] + r * np.sin(thr), o[2]))]
        c += 1
    return at


def molecule(n):
    """The benchmark geometry: a custom molecule if OQP_SYSTEM_XYZ is set (a file
    of `SYM x y z` lines in Angstrom), else the (H2O)n cluster. Shared by the grid
    generator and nativeprep so both see the SAME atoms."""
    import os
    path = os.environ.get("OQP_SYSTEM_XYZ")
    if not path:
        return water_cluster(n)
    at = []
    for ln in open(path):
        p = ln.split()
        if len(p) >= 4:
            at.append((p[0], (float(p[1]), float(p[2]), float(p[3]))))
    return at


if __name__ == "__main__":
    # usage: grid_native.py NWATER OUT.bin   (water benchmark geometry)
    nw = int(sys.argv[1]); out = sys.argv[2]
    nang = int(sys.argv[3]) if len(sys.argv)>3 else 302
    atoms = molecule(nw)
    P, W = build_grid(atoms, nang_max=nang)
    write_grid_bin(out, P, W)
    print(f"[grid_native] w{nw}: npts={len(W)}  sumW={W.sum():.6f}  -> {out}", flush=True)
