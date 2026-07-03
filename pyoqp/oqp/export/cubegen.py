"""Gaussian .cube export for MRSF orbitals and excited-state densities.

Uses the self-contained GTO evaluator (``oqp.analysis.gto_grid``).  Density
cubes store rho(r) = sum_uv D_uv chi_u(r) chi_v(r); orbital cubes store the
amplitude phi(r).  Grid integrals are validated against analytic traces:
state density -> N, transition density -> 0, attachment/detachment ->
n_promoted, |orbital|^2 -> 1.
"""
import numpy as np
from oqp.analysis.gto_grid import AOBasis, make_box_grid

__all__ = ["CubeExporter"]

_CHUNK = 200000


def _density_values(ao, D_ao, points):
    """rho(r) = sum_uv D_uv chi_u chi_v on points (chunked)."""
    npts = points.shape[0]
    out = np.empty(npts)
    Dsym = 0.5 * (D_ao + D_ao.T)
    for i in range(0, npts, _CHUNK):
        A = ao.eval_ao(points[i:i + _CHUNK])
        out[i:i + _CHUNK] = np.einsum("pm,pm->p", A @ Dsym, A)
    return out


def _orbital_values(ao, c_ao, points):
    npts = points.shape[0]
    out = np.empty(npts)
    for i in range(0, npts, _CHUNK):
        out[i:i + _CHUNK] = ao.eval_ao(points[i:i + _CHUNK]) @ c_ao
    return out


def _write_cube(path, comment1, comment2, Z, coords, origin, n, dvec, values):
    """Write a Gaussian cube. ``values`` is flat in (ix slow, iy, iz fast) order."""
    nx, ny, nz = n
    dx, dy, dz = dvec
    natom = len(Z)
    with open(path, "w") as f:
        f.write(comment1.rstrip() + "\n")
        f.write(comment2.rstrip() + "\n")
        f.write(f"{natom:5d} {origin[0]:12.6f} {origin[1]:12.6f} {origin[2]:12.6f}\n")
        f.write(f"{nx:5d} {dx:12.6f} {0.0:12.6f} {0.0:12.6f}\n")
        f.write(f"{ny:5d} {0.0:12.6f} {dy:12.6f} {0.0:12.6f}\n")
        f.write(f"{nz:5d} {0.0:12.6f} {0.0:12.6f} {dz:12.6f}\n")
        for z, c in zip(Z, coords):
            f.write(f"{int(z):5d} {float(z):12.6f} {c[0]:12.6f} {c[1]:12.6f} {c[2]:12.6f}\n")
        v = values.reshape(nx, ny, nz)
        for ix in range(nx):
            for iy in range(ny):
                col = v[ix, iy, :]
                for k in range(0, nz, 6):
                    f.write("".join(f"{x:13.5e}" for x in col[k:k + 6]) + "\n")


class CubeExporter:
    def __init__(self, states, ao=None, padding=5.0, spacing=0.15):
        self.st = states
        self.ao = ao if ao is not None else AOBasis(states.mol)
        self.Z = np.asarray(states.mol.get_atoms(), dtype=int)
        self.coords = self.ao.coords
        self.origin, self.n, self.dvec, self.points = make_box_grid(
            self.coords, padding=padding, spacing=spacing)
        self.dV = self.dvec[0] * self.dvec[1] * self.dvec[2]

    # ---- integral (trace) checks; can use a finer dedicated grid ----
    def integrate_density(self, D_ao, spacing=None, padding=5.0):
        pts, dV = self.points, self.dV
        if spacing is not None:
            _, _, dvec, pts = make_box_grid(self.coords, padding=padding, spacing=spacing)
            dV = dvec[0] * dvec[1] * dvec[2]
        total = 0.0
        for i in range(0, pts.shape[0], _CHUNK):
            total += _density_values(self.ao, D_ao, pts[i:i + _CHUNK]).sum() * dV
        return total

    def integrate_orbital_sq(self, c_ao, spacing=None, padding=5.0):
        pts, dV = self.points, self.dV
        if spacing is not None:
            _, _, dvec, pts = make_box_grid(self.coords, padding=padding, spacing=spacing)
            dV = dvec[0] * dvec[1] * dvec[2]
        total = 0.0
        for i in range(0, pts.shape[0], _CHUNK):
            phi = _orbital_values(self.ao, c_ao, pts[i:i + _CHUNK])
            total += np.sum(phi * phi) * dV
        return total

    # ---- cube writers ----
    def _density_cube(self, path, D_ao, label):
        vals = _density_values(self.ao, D_ao, self.points)
        _write_cube(path, f"OQP MRSF density: {label}", "rho(r), a.u.",
                    self.Z, self.coords, self.origin, self.n, self.dvec, vals)
        return float(vals.sum() * self.dV)

    def _orbital_cube(self, path, c_ao, label):
        vals = _orbital_values(self.ao, c_ao, self.points)
        _write_cube(path, f"OQP MRSF orbital: {label}", "phi(r), a.u.",
                    self.Z, self.coords, self.origin, self.n, self.dvec, vals)
        return float(np.sum(vals * vals) * self.dV)

    def mo_cube(self, path, mo_index, spin="alpha"):
        c = self.st.C[:, mo_index]
        return self._orbital_cube(path, c, f"MO {mo_index} ({spin})")

    def state_density_cube(self, path, n):
        return self._density_cube(path, self.st.state_density_ao(n), f"state S{n} density")

    def transition_density_cube(self, path, i, j):
        return self._density_cube(path, self.st.tdm_ao(i, j), f"transition {i}->{j}")

    def attachment_detachment_cubes(self, path_attach, path_detach, ad):
        a = self._density_cube(path_attach, ad["A_ao"], f"attachment S{ad['state']}")
        d = self._density_cube(path_detach, ad["D_ao"], f"detachment S{ad['state']}")
        return a, d

    def nto_cube(self, path, c_ao, label):
        return self._orbital_cube(path, c_ao, label)
