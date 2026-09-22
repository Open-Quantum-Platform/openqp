"""Magnetically induced current density and ACID (AICD) maps.

Builds the first-order current-density tensor field from the GIAO CPHF/CPKS
density response published by ``nmr_giao_shielding`` as ``OQP::nmr_pdens``, and
writes it as Gaussian cubes: one ACID scalar plus the three components of the
current-density vector for a chosen field direction.

Conventions
-----------
Atomic units, with the magnetic perturbation written as ``H1 = (1/2) B . L_O``
so that the factor c^-1 is absorbed into the field; shieldings derived from this
field acquire the usual alpha^2 on conversion to ppm.  The gauge origin is the
coordinate origin, matching the GIAO shielding driver (``o0 = 0``).

The current density is a second-rank tensor field ``J_ij(r) = dJ_i/dB_j`` with
four contributions,

    J1  density response      - sum_uv P^j_uv (d_i chi_u) chi_v
    J2  diamagnetic           - (1/2) rho(r) eps_ijm r_m
    J3  London phase gradient   (1/4) eps_ijm sum_uv D_uv (R_u + R_v)_m chi_u chi_v
    J4  London phase difference (1/2) eps_ljq r_l sum_uv D_uv (R_u - R_v)_q
                                       (d_i chi_u) chi_v

J3 and J4 exist only because the London factor ``exp(-i/2 B.(R_u x r))`` makes
the basis field-dependent; they need the AO centres together with the
unperturbed density matrix, so they cannot be reconstructed from exported
orbitals alone.  Dropping them is not a small correction -- on a displaced
methane they shift the reconstructed shielding by 72 and 7 ppm respectively.

ACID(r) is the anisotropy invariant of J_ij(r) in the normalisation of Herges
and Geuenich (see ``acid_scalar``), so isosurface values are directly the ones
used in the literature -- 0.05 a.u. is the standard plotting value.  Because it
is a tensor invariant it does not depend on the direction of B, which is what
makes ACID isosurfaces comparable between molecules in a way NICS is not.

Note on gauge: ACID must be built on the GIAO path.  With a common gauge origin
the map is badly gauge-contaminated -- benzene NICS(0) comes out near -105 ppm
in 6-31G* and is still ~43 ppm off at 6-311++G**, against -10 ppm for GIAO.

Basis limits: the grid evaluator (``oqp.analysis.gto_grid.AOBasis``) walks
Cartesian components up to f, so a spherical-harmonic basis needs
``[input] ispher=false`` and a Cartesian basis containing g or higher shells is
not supported.  Both are rejected with an explicit message rather than a partial
map.  The native ``OQP::nmr_pdens`` response itself carries no such restriction.
"""
import numpy as np

from oqp.analysis.gto_grid import AOBasis, make_box_grid
from .cubegen import _write_cube

__all__ = ["AcidExporter"]

_CHUNK = 20000

_LEVI_CIVITA = np.zeros((3, 3, 3))
_LEVI_CIVITA[0, 1, 2] = _LEVI_CIVITA[1, 2, 0] = _LEVI_CIVITA[2, 0, 1] = 1.0
_LEVI_CIVITA[0, 2, 1] = _LEVI_CIVITA[2, 1, 0] = _LEVI_CIVITA[1, 0, 2] = -1.0


def _f_order(mol, tag, shape):
    """Reinterpret a tagarray buffer with the given Fortran shape."""
    raw = np.array(mol.data[tag], copy=True).ravel(order="C")
    return raw.reshape(shape, order="F")


def _unpack_lt(packed, n):
    m = np.zeros((n, n))
    rows, cols = np.tril_indices(n)
    m[rows, cols] = packed
    m[cols, rows] = packed
    return m


def acid_scalar(jten):
    """ACID on the standard scale, so 0.05 a.u. means what the literature means.

    Herges and Geuenich, J. Phys. Chem. A 105, 3214 (2001), in the form the
    GIMIC reference implementation uses (``src/libgimic/acid.f90``):

        ACID^2 = 1/3 [(t_xx-t_yy)^2 + (t_yy-t_zz)^2 + (t_zz-t_xx)^2]
               + 1/2 [(t_xy+t_yx)^2 + (t_xz+t_zx)^2 + (t_yz+t_zy)^2]

    Only the symmetric part of the tensor survives this contraction.  That is
    the point of the method rather than an approximation: the antisymmetric part
    carries the diamagnetic circulation, which is proportional to the electron
    density and so peaks at the nuclei, and it is exactly what would otherwise
    bury the delocalisation the map is for.
    """
    t = jten
    diag = ((t[..., 0, 0] - t[..., 1, 1]) ** 2
            + (t[..., 1, 1] - t[..., 2, 2]) ** 2
            + (t[..., 2, 2] - t[..., 0, 0]) ** 2)
    off = ((t[..., 0, 1] + t[..., 1, 0]) ** 2
           + (t[..., 0, 2] + t[..., 2, 0]) ** 2
           + (t[..., 1, 2] + t[..., 2, 1]) ** 2)
    return np.sqrt(diag / 3.0 + off / 2.0)


def current_vector(jten, bfield):
    """Current-density vector field for one magnetic-field direction."""
    b = np.asarray(bfield, dtype=float)
    norm = np.linalg.norm(b)
    if norm == 0.0:
        raise ValueError("magnetic field direction must be non-zero")
    return jten @ (b / norm)


class AcidExporter:
    """ACID / current-density cubes for a finished GIAO NMR calculation."""

    def __init__(self, mol, ao=None, padding=5.0, spacing=0.20):
        self.mol = mol
        self.ao = ao if ao is not None else AOBasis(mol)
        self.nbf = self.ao.nbf
        try:
            self.pmat = _f_order(mol, "OQP::nmr_pdens", (3, self.nbf, self.nbf))
        except AttributeError:
            # OQPData raises AttributeError for an absent tag; say what is
            # actually missing rather than surfacing the tagarray lookup.
            raise ValueError(
                "OQP::nmr_pdens is absent; ACID needs a completed GIAO NMR run "
                "(properties=nmr with nmr_gauge=giao).") from None
        self.dm = self._total_density()
        self.rcen = self.ao.coords[self.ao.ao_atom]
        self.Z = np.asarray(mol.get_atoms(), dtype=int)
        self.coords = self.ao.coords
        self.origin, self.n, self.dvec, self.points = make_box_grid(
            self.coords, padding=padding, spacing=spacing)

    def _total_density(self):
        """Total AO density, matching the GIAO shielding driver's convention.

        For RHF ``OQP::DM_A`` is already the closed-shell total; for UHF/ROHF
        the total is D_alpha + D_beta.  DM_B is allocated and non-zero even in
        an RHF run, so branch on the SCF type rather than on its contents.
        """
        nbf2 = self.nbf * (self.nbf + 1) // 2
        dm = _unpack_lt(np.asarray(self.mol.data["OQP::DM_A"]).ravel()[:nbf2],
                        self.nbf)
        if int(self.mol.data["scftype"]) != 1:
            dm = dm + _unpack_lt(
                np.asarray(self.mol.data["OQP::DM_B"]).ravel()[:nbf2], self.nbf)
        return 0.5 * (dm + dm.T)

    def current_density(self, points):
        """First-order current-density tensor J_ij(r) -> (npts, 3, 3)."""
        pts = np.asarray(points, dtype=float).reshape(-1, 3)
        out = np.empty((pts.shape[0], 3, 3))
        for lo in range(0, pts.shape[0], _CHUNK):
            block = pts[lo:lo + _CHUNK]
            ao = self.ao.eval_ao_deriv1(block)
            ao0, ao1 = ao[0], ao[1:]
            w = ao0 @ self.dm

            j = -np.einsum("ipm,jmn,pn->pij", ao1, self.pmat, ao0, optimize=True)

            rho = np.einsum("pm,pm->p", ao0, w, optimize=True)
            j -= 0.5 * rho[:, None, None] * np.einsum(
                "ijm,pm->pij", _LEVI_CIVITA, block, optimize=True)

            gvec = 2.0 * np.einsum("pm,mq->pq", ao0 * w, self.rcen,
                                   optimize=True)
            j += 0.25 * np.einsum("ijm,pm->pij", _LEVI_CIVITA, gvec,
                                  optimize=True)

            amat = np.einsum("ipm,pm,mq->piq", ao1, w, self.rcen, optimize=True)
            vstack = np.einsum("pm,mq,mn->qpn", ao0, self.rcen, self.dm,
                               optimize=True)
            bmat = np.einsum("ipm,qpm->piq", ao1, vstack, optimize=True)
            j += 0.5 * np.einsum("ljq,pl,piq->pij", _LEVI_CIVITA, block,
                                 amat - bmat, optimize=True)

            out[lo:lo + _CHUNK] = j
        return out

    def acid(self, points=None):
        pts = self.points if points is None else points
        return acid_scalar(self.current_density(pts))

    def shielding_from_current(self, atom_index, points, weights):
        """Biot-Savart reconstruction of a shielding tensor, in atomic units.

        sigma_ij = - int d3r [ (r - R_N) x J_.j(r) ]_i / |r - R_N|^3

        Multiply by alpha^2 * 1e6 for ppm.  Integrating the plotted field back
        into the shielding the SCF already knows is the check that validates the
        current density itself rather than only its ingredients.
        """
        pts = np.asarray(points, dtype=float).reshape(-1, 3)
        rvec = pts - self.coords[atom_index]
        dist = np.linalg.norm(rvec, axis=1)
        kernel = np.where(dist > 1e-10, 1.0 / np.maximum(dist, 1e-10) ** 3, 0.0)
        jten = self.current_density(pts)
        cross = np.einsum("ilm,pl,pmj->pij", _LEVI_CIVITA, rvec, jten,
                          optimize=True)
        return -np.einsum("p,pij->ij", np.asarray(weights) * kernel, cross,
                          optimize=True)

    def write_cubes(self, prefix, bfield=(0.0, 0.0, 1.0)):
        """Write ``<prefix>_acid.cube`` and ``<prefix>_j{x,y,z}.cube``."""
        jten = self.current_density(self.points)
        scalar = acid_scalar(jten)
        jvec = current_vector(jten, bfield)
        b = np.asarray(bfield, dtype=float)
        b = b / np.linalg.norm(b)
        head = "GIAO current density, atomic units"
        paths = [f"{prefix}_acid.cube"]
        _write_cube(paths[0], "OQP ACID: anisotropy of the induced current density",
                    head, self.Z, self.coords, self.origin, self.n, self.dvec,
                    scalar)
        for k, tag in enumerate("xyz"):
            path = f"{prefix}_j{tag}.cube"
            _write_cube(path, f"OQP induced current density J{tag}",
                        f"{head}; B = ({b[0]:.3f}, {b[1]:.3f}, {b[2]:.3f})",
                        self.Z, self.coords, self.origin, self.n, self.dvec,
                        jvec[:, k])
            paths.append(path)
        return paths
