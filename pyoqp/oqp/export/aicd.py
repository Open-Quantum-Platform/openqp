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

Provenance: ``OQP::nmr_pdens`` is a persistent record -- it stays on the
molecule and in the .oqp file until another GIAO run replaces it, and nothing
else invalidates it.  A same-size SCF, a moved geometry or a CGO NMR call
therefore leaves behind a response that still loads at the right shape.  The
driver stamps the basis size, atom count, geometry and density fingerprint it
built the response from into ``OQP::nmr_pdens_ref``; ``AcidExporter`` refuses
anything that does not match the molecule it is asked to plot, rather than
quietly mixing one geometry's response with another's coordinates.

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


def unit_field(bfield):
    """Validated unit vector for the magnetic-field direction.

    One place that turns whatever the caller passed into three finite numbers,
    so the field used to contract the tensor and the field written into the cube
    header cannot disagree about shape.
    """
    b = np.asarray(bfield, dtype=float).reshape(-1)
    if b.size != 3 or not np.all(np.isfinite(b)):
        raise ValueError(
            "magnetic field direction must be three finite components; got "
            f"{bfield!r}")
    norm = float(np.linalg.norm(b))
    # A non-finite component would give a non-finite norm, which is not zero,
    # so an equality test against zero alone lets NaN through into every cube.
    if not np.isfinite(norm) or norm == 0.0:
        raise ValueError("magnetic field direction must be non-zero")
    return b / norm


def current_vector(jten, bfield):
    """Current-density vector field for one magnetic-field direction."""
    return jten @ unit_field(bfield)


class AcidExporter:
    """ACID / current-density cubes for a finished GIAO NMR calculation."""

    def __init__(self, mol, ao=None, padding=5.0, spacing=0.20):
        # This is public API, so it validates its own arguments rather than
        # relying on the workflow that usually calls it: nan and inf survive
        # float() and slip past an ordering test, and make_box_grid would then
        # build a degenerate box of non-finite coordinates and write cubes
        # nothing can read.
        spacing = float(spacing)
        padding = float(padding)
        if not (np.isfinite(spacing) and np.isfinite(padding)):
            raise ValueError(
                f"grid spacing and padding must be finite; got spacing={spacing}, "
                f"padding={padding}.")
        if spacing <= 0.0 or padding < 0.0:
            raise ValueError(
                f"grid spacing must be positive and padding non-negative; got "
                f"spacing={spacing}, padding={padding}.")
        self.mol = mol
        self.ao = ao if ao is not None else AOBasis(mol)
        self.nbf = self.ao.nbf
        self.coords = self.ao.coords
        self.rcen = self.ao.coords[self.ao.ao_atom]
        self.Z = np.asarray(mol.get_atoms(), dtype=int)
        self.dm = self._total_density()
        try:
            mol.data["OQP::nmr_pdens"]
        except AttributeError:
            # OQPData raises AttributeError for an absent tag; say what is
            # actually missing rather than surfacing the tagarray lookup.
            raise ValueError(
                "OQP::nmr_pdens is absent; ACID needs a completed GIAO NMR run "
                "(properties=nmr with nmr_gauge=giao).") from None
        # Provenance before shape: a response from a different molecule can
        # still carry the right shape, and when it does not, "built for 7 basis
        # functions" is the message worth having rather than a reshape error.
        self._verify_provenance()
        self.pmat = _f_order(mol, "OQP::nmr_pdens", (3, self.nbf, self.nbf))
        self.origin, self.n, self.dvec, self.points = make_box_grid(
            self.coords, padding=padding, spacing=spacing)

    def _verify_provenance(self):
        """Refuse a response that does not belong to this molecule.

        ``OQP::nmr_pdens`` is persistent: it lives on in the molecule and in
        the .oqp file, and only another GIAO run replaces it.  A same-size SCF,
        a moved geometry or a CGO NMR call therefore leaves the old response
        in place, still the right shape and still loadable, and combining it
        with the current density and coordinates gives a map that is wrong
        without looking wrong.  The GIAO driver stamps what the response
        belongs to (``OQP::nmr_pdens_ref``); this is the other half of that
        contract.
        """
        try:
            ref = np.asarray(self.mol.data["OQP::nmr_pdens_ref"],
                             dtype=float).ravel()
        except AttributeError:
            raise ValueError(
                "OQP::nmr_pdens carries no provenance stamp, so it cannot be "
                "shown to belong to this molecule; it predates the ACID "
                "export or was written by another code path.  Re-run the GIAO "
                "NMR calculation.") from None
        nat = self.coords.shape[0]
        # Written invalid when the buffer is allocated and filled in only once
        # the response is complete, so a run that aborted mid-way reads as
        # stale here rather than as current.
        if ref.size < 4 or ref[0] < 0.0:
            raise ValueError(
                "the GIAO magnetic response is marked incomplete: the run that "
                "allocated it did not finish.  Re-run the GIAO NMR "
                "calculation.")
        if ref.size < 5 or int(round(ref[1])) != nat:
            raise ValueError(
                f"the GIAO magnetic response was built for a molecule with "
                f"{int(round(ref[1])) if ref.size >= 2 else '?'} atoms, but "
                f"this one has {nat}.  Re-run the GIAO NMR calculation.")
        if int(round(ref[0])) != self.nbf:
            raise ValueError(
                f"the GIAO magnetic response was built for "
                f"{int(round(ref[0]))} basis functions, but this molecule now "
                f"has {self.nbf}.  Re-run the GIAO NMR calculation.")
        nsh = int(round(ref[4]))
        nprim = int(round(ref[5]))
        head = 6 + 3 * nat
        if ref.size != head + 3 * nsh + 2 * nprim:
            raise ValueError(
                "the GIAO magnetic response carries a malformed provenance "
                f"stamp ({ref.size} entries for {nat} atoms, {nsh} shells and "
                f"{nprim} primitives).  Re-run the GIAO NMR calculation.")
        moved = float(np.abs(ref[6:head].reshape(nat, 3) - self.coords).max())
        if moved > 1.0e-8:
            raise ValueError(
                f"the geometry moved by {moved:.3e} bohr since the GIAO "
                f"magnetic response was computed, so the response no longer "
                f"matches these coordinates.  Re-run the GIAO NMR "
                f"calculation.")
        # The AO ordering itself: nbf and the density invariants below are
        # blind to a permutation of the basis (trace and sum of squares survive
        # D -> P D P^T), so a same-size basis that came back with its shells in
        # a different order would pass everything else while the response is
        # indexed in the old order.
        got_sh, got_prim = self._basis_signature()
        if got_sh.shape[0] != nsh or got_prim.shape[1] != nprim:
            raise ValueError(
                f"the GIAO magnetic response was built in a basis of {nsh} "
                f"shells and {nprim} primitives, but this molecule now has "
                f"{got_sh.shape[0]} and {got_prim.shape[1]}.  Re-run the GIAO "
                f"NMR calculation.")
        want_sh = ref[head:head + 3 * nsh].reshape(nsh, 3)
        bad = np.flatnonzero(np.any(got_sh != want_sh, axis=1))
        if bad.size:
            s = int(bad[0])
            raise ValueError(
                f"the basis changed since the GIAO magnetic response was "
                f"computed: shell {s} is now (centre, am, ncontr) = "
                f"({got_sh[s, 0]:.0f}, {got_sh[s, 1]:.0f}, {got_sh[s, 2]:.0f}) "
                f"against ({want_sh[s, 0]:.0f}, {want_sh[s, 1]:.0f}, "
                f"{want_sh[s, 2]:.0f}) when the response was built, so the "
                f"response is indexed in a different AO order.  Re-run the "
                f"GIAO NMR calculation.")
        # Primitives in full, not reduced: a per-shell summary is not
        # injective, so exponents [1, 2] and [1.5, 1.5] would share it.  These
        # are the same stored doubles on both sides, hence the tight tolerance
        # -- it absorbs serialization, not a real difference.
        want_prim = ref[head + 3 * nsh:].reshape(2, nprim)
        bad = np.flatnonzero(np.any(
            np.abs(got_prim - want_prim)
            > 1.0e-13 * np.maximum(np.abs(want_prim), 1.0), axis=0))
        if bad.size:
            k = int(bad[0])
            raise ValueError(
                f"the basis primitives changed since the GIAO magnetic "
                f"response was computed: primitive {k} is now (exponent, "
                f"coefficient) = ({got_prim[0, k]:.10g}, {got_prim[1, k]:.10g}) "
                f"against ({want_prim[0, k]:.10g}, {want_prim[1, k]:.10g}) when "
                f"the response was built, so the AO functions it is indexed in "
                f"are not these.  Re-run the GIAO NMR calculation.")
        # The density fingerprint is summed in a different order here than in
        # the driver, so it is compared at a tolerance -- it is a detector of
        # a changed wavefunction, not a checksum.  Any real change (basis,
        # charge, spin state, a re-converged SCF) moves it enormously.
        got = np.array([float(np.trace(self.dm)), float(np.sum(self.dm ** 2))])
        want = ref[2:4]
        scale = np.maximum(np.abs(want), 1.0)
        if np.any(np.abs(got - want) > 1.0e-8 * scale):
            raise ValueError(
                "the electronic density changed since the GIAO magnetic "
                "response was computed (density fingerprint "
                f"{got[0]:.10g}/{got[1]:.10g} against {want[0]:.10g}/"
                f"{want[1]:.10g}), so the response belongs to a different "
                "wavefunction.  Re-run the GIAO NMR calculation.")

    def _basis_signature(self):
        """The basis as the grid evaluator sees it: shells and primitives.

        Returns per-shell ``(centre, am, ncontr)`` in AO order and the full
        ``(exponents, coefficients)``, from the same ``oqp_get_basis`` arrays
        the evaluator builds its shells from.  This is the identity of the AO
        functions the response is indexed in, not a summary of them.
        """
        b = self.mol.data.get_basis()
        shells = np.column_stack((
            np.asarray(b["centers"], dtype=float),
            np.asarray(b["angs"], dtype=float),
            np.asarray(b["ncontr"], dtype=float)))
        prims = np.vstack((np.asarray(b["alpha"], dtype=float),
                           np.asarray(b["coef"], dtype=float)))
        return shells, prims

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
        # Normalise before anything is written: a field this cannot use must
        # fail with no files on disk, not after the scalar cube has landed.
        b = unit_field(bfield)
        jten = self.current_density(self.points)
        scalar = acid_scalar(jten)
        jvec = jten @ b
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
