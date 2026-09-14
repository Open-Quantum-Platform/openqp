import openmm.app as app
import openmm as mm
import openmm.unit as unit
import os

import numpy as np

from oqp.openqp import OPENQP
from oqp.library.single_point import (
    SinglePoint, Gradient, Hessian, LastStep,
    BasisOverlap, NACME, NAC
)
from oqp.utils.file_utils import dump_log, dump_data, write_config, write_xyz
from oqp.utils.tb_backends import is_tb_method
from oqp.utils.state_labels import is_mrsf, public_state_label
from oqp.library.qmmm_ewald import EwaldQMMM
from oqp.library.qmmm_connectivity import (
    detect_link_atoms, link_atom_position,
    redistribute_frontier_charges, assemble_embedding_sites,
    select_boundary_switching,
    SWSCALE_WHOLE_MOLECULE, SWSCALE_COVALENT_BOUNDARY,
)

import oqp
from oqp.library.ints_1e import ints_1e


def unpack_lower_tri_single(packed_atom, nbf):
    """
    Unpack a single lower-triangular packed array (length nbf*(nbf+1)/2)
    into a full symmetric (nbf x nbf) matrix.
    """
    packed_atom = np.asarray(packed_atom)
    nbf_tri = nbf * (nbf + 1) // 2
    if packed_atom.size != nbf_tri:
        raise ValueError(f"Size mismatch: got {packed_atom.size}, expected {nbf_tri}")
    full = np.zeros((nbf, nbf), dtype=packed_atom.dtype)
    idx = 0
    for i in range(nbf):
        for j in range(i + 1):
            val = packed_atom[idx]
            full[i, j] = val
            full[j, i] = val
            idx += 1
    return full

def unpack_lower_tri_multi(packed, nbf, natm):
    packed = np.asarray(packed)
    nbf_tri = nbf * (nbf + 1) // 2
    flat = packed.ravel(order="C")
    if flat.size != natm * nbf_tri:
        raise ValueError(f"Size mismatch: got {flat.size}, expected {natm * nbf_tri}")
    packed_by_atom = flat.reshape((natm, nbf_tri))
    full_all = np.zeros((natm, nbf, nbf), dtype=packed.dtype)
    for a in range(natm):
        full_all[a] = unpack_lower_tri_single(packed_by_atom[a], nbf)
    return full_all

def pack_lower_tri_single(full):
    full = np.asarray(full)
    if full.shape[0] != full.shape[1]:
        raise ValueError("Matrix must be square")
    nbf = full.shape[0]
    nbf_tri = nbf * (nbf + 1) // 2
    packed = np.zeros(nbf_tri, dtype=full.dtype)
    idx = 0
    for i in range(nbf):
        for j in range(i + 1):
            packed[idx] = full[i, j]
            idx += 1
    return packed


def read_xyz(filepath):
    """
    Read a standard XYZ file.

    Returns
    -------
    symbols : list of str
    coords : np.ndarray, shape (natom, 3)  –  Angstroms
    """
    symbols = []
    coords = []
    with open(filepath, "r") as f:
        natom = int(f.readline().strip())
        f.readline()  # comment line
        for _ in range(natom):
            parts = f.readline().split()
            symbols.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.array(coords)


#: OpenMM nonbonded methods that impose periodic boundary conditions.  The
#: QM/MM electrostatics go through the Ewald branch only for these; NoCutoff
#: and CutoffNonPeriodic are a finite cluster even when the topology carries
#: box vectors.  (LJPME is not accepted by the [qmmm] cutoff parser.)
PERIODIC_METHODS = (app.PME, app.Ewald, app.CutoffPeriodic)

def smeared_coulomb(r, mu):
    """Pair kernel of a QM point charge with a Gaussian-smeared MM charge of
    damping parameter ``mu`` (1/bohr): returns (phi_kernel, force_kernel) with
    phi_kernel = erf(mu r)/r and force_kernel = -d(phi_kernel)/dr / r, so that
    the force on the QM centre is q_A Q_M * force_kernel * d.  ``mu=None`` is
    the point-charge limit 1/r and 1/r^3.  Both kernels are finite at r = 0
    (phi -> 2 mu/sqrt(pi), force -> 0), which is the point of the smearing."""
    from scipy.special import erf
    r = np.asarray(r, dtype=float)
    if mu is None:
        return 1.0 / r, 1.0 / r ** 3
    tiny = r < 1e-8
    rs = np.where(tiny, 1.0, r)
    phi = np.where(tiny, 2.0 * mu / np.sqrt(np.pi), erf(mu * rs) / rs)
    force = np.where(tiny, 0.0, erf(mu * rs) / rs ** 3
                     - 2.0 * mu / np.sqrt(np.pi) * np.exp(-(mu * rs) ** 2) / rs ** 2)
    return phi, force


#: Accepted [qmmm] embedding spellings (ported from PR #274).
_VALID_EMBEDDINGS = {"mechanical", "electrostatic", "espf", "espf_full", "split"}


def _normalize_embedding(value):
    """Lower-cased, stripped embedding keyword; rejects unknown spellings so a
    typo cannot silently select a different scheme."""
    embedding = str(value).strip().lower()
    if embedding not in _VALID_EMBEDDINGS:
        choices = ", ".join(sorted(_VALID_EMBEDDINGS))
        raise ValueError(f"Unknown QM/MM embedding '{value}'. Choices: {choices}")
    return embedding


def anderson_step(q_hist, f_hist, m=3, beta=1.0, max_step=0.5):
    """Anderson-accelerated update for a charge fixed point q = g(q).

    ``q_hist`` holds the last iterates q_k and ``f_hist`` the residuals
    f_k = g(q_k) - q_k (same length, newest last).  Returns the next iterate
    q_k + beta f_k - (dQ + beta dF) gamma with gamma the least-squares
    combination of the last ``m`` residual differences; falls back to the
    plain damped step (beta/2) when the history is too short, the
    least-squares problem is degenerate, or the extrapolated move exceeds
    ``max_step`` electrons on any atom."""
    q_k, f_k = np.asarray(q_hist[-1], dtype=float), np.asarray(f_hist[-1], dtype=float)
    n = min(m, len(q_hist) - 1)
    if n >= 1:
        dQ = np.column_stack([np.asarray(q_hist[-1 - i]) - np.asarray(q_hist[-2 - i]) for i in range(n)])
        dF = np.column_stack([np.asarray(f_hist[-1 - i]) - np.asarray(f_hist[-2 - i]) for i in range(n)])
        gamma, *_ = np.linalg.lstsq(dF, f_k, rcond=1e-10)
        q_new = q_k + beta * f_k - (dQ + beta * dF) @ gamma
        if np.all(np.isfinite(q_new)) and float(np.abs(q_new - q_k).max()) <= max_step:
            return q_new
    return q_k + 0.5 * f_k


def is_periodic_method(cutoff):
    """True when ``cutoff`` (an OpenMM nonbonded-method constant) is periodic."""
    return any(cutoff is m for m in PERIODIC_METHODS)


def _periodic_nonbonded_cutoff(topology, cutoff_method):
    """Return a safe OpenMM nonbonded cutoff for the current periodic box."""
    if cutoff_method is app.NoCutoff:
        return 1.0 * unit.nanometer
    vectors = topology.getPeriodicBoxVectors()
    if vectors is None:
        return 1.0 * unit.nanometer
    lengths = []
    for vec in vectors:
        xyz = vec.value_in_unit(unit.nanometer)
        lengths.append(float(np.linalg.norm(xyz)))
    min_len = min(lengths)
    return min(1.0, 0.4 * min_len) * unit.nanometer


def fold_virtual_site_forces(system, forces, positions=None, sites=None):
    """Move the force on every virtual site of ``system`` onto the particles
    that place it and zero the site row: the chain rule of the site position,
    so the real particles carry the whole derivative.  ``forces`` is (N, 3) in
    any force unit; ``positions`` (nm) are needed only for out-of-plane sites.

    An average site r_v = sum_k w_k r_k gives F_k += w_k F_v.  An out-of-plane
    site r_v = r_1 + w12 r_12 + w13 r_13 + wc (r_12 x r_13) gives
    F_2 += w12 F_v + wc (r_13 x F_v), F_3 += w13 F_v - wc (r_12 x F_v) and
    F_1 += F_v - (the other two), so the total force is kept.  A
    local-coordinates site raises NotImplementedError."""
    q_unit = None
    if unit.is_quantity(forces):
        q_unit = forces.unit
        forces = forces.value_in_unit(q_unit)
    f = np.array(forces, dtype=float)
    if sites is None:
        sites = [i for i in range(system.getNumParticles()) if system.isVirtualSite(i)]
    X = None
    for i in sites:
        site = system.getVirtualSite(i)
        name = type(site).__name__
        fv = f[i].copy()
        if name in ("TwoParticleAverageSite", "ThreeParticleAverageSite"):
            for k in range(site.getNumParticles()):
                f[site.getParticle(k)] += site.getWeight(k) * fv
        elif name == "OutOfPlaneSite":
            if X is None:
                if positions is None:
                    raise ValueError("positions are required to fold forces on an out-of-plane site")
                X = np.asarray(positions.value_in_unit(unit.nanometer)
                               if unit.is_quantity(positions) else positions, dtype=float)
            p1, p2, p3 = (site.getParticle(k) for k in range(3))
            r12, r13 = X[p2] - X[p1], X[p3] - X[p1]
            wc = site.getWeightCross()
            f2 = site.getWeight12() * fv + wc * np.cross(r13, fv)
            f3 = site.getWeight13() * fv - wc * np.cross(r12, fv)
            f[p1] += fv - f2 - f3
            f[p2] += f2
            f[p3] += f3
        else:
            raise NotImplementedError(
                f"virtual site {i} is a {name}; QM/MM forces are folded for average "
                "and out-of-plane sites (e.g. TIP4P, TIP5P) only.")
        f[i] = 0.0
    return f * q_unit if q_unit is not None else f


class OpenQpQMMM:
    """
    Low-level QM/MM driver using OpenMM (MM) + OpenQP (QM).

    Two modes:
      1. **Config mode** – ``oqp_cfg`` (dict).
      2. **Mol mode**    – pre-built ``mol`` object.

    Exactly one of ``oqp_cfg`` / ``mol`` must be provided.
    """

    def __init__(
        self,
        positions,
        topology,
        forcefield,
        qm_atoms,
        oqp_cfg=None,
        mol=None,
        Cutoff=app.NoCutoff,
        Embedding='mechanical',
        frontier_scheme='none',
        ewald_tol=None,
        lj_switch=False,
        h_lj=False,
        mm_charge_width=None,
    ):
        # Gaussian-smeared MM charges for the QM-MM electrostatics (width in
        # Angstrom; None = point charges).  The QM-MM pair potential becomes
        # erf(mu r)/r with mu = 1/(sqrt(2) w), which removes the 1/r singularity
        # a bare MM charge presents to the QM density -- the reference Tinker
        # ESPF code applies the same erf damping (its ERFMU keyword).  Energy and
        # force are modified consistently (direct sum and Ewald real-space part);
        # MM-MM interactions and the QM-image term are untouched.
        if mm_charge_width in (None, 0, 0.0):
            self.mm_damp_mu = None
        else:
            w = float(mm_charge_width)
            if not np.isfinite(w) or w <= 0.0:
                raise ValueError(
                    f"[qmmm] mm_charge_width must be a finite positive width in "
                    f"Angstrom (or 0/unset for point charges); got {mm_charge_width!r}")
            self.mm_damp_mu = 1.0 / (np.sqrt(2.0) * w * 1.8897259886)
        # Give Lennard-Jones parameters to MM hydrogens that have none (TIP3P
        # water H: sigma 1 nm / epsilon 0 in the AMBER XML).  Without them
        # nothing keeps a water hydrogen from collapsing onto a QM oxygen (the
        # QM density has no Pauli wall against a bare point charge), which
        # produces unphysical 1.4-1.5 A contacts that destabilise the MRSF
        # response.  Uses the CHARMM TIP3P HT values (Rmin/2 = 0.2245 A,
        # eps = 0.046 kcal/mol).  Off by default; [qmmm] h_lj=true.
        self.h_lj = bool(h_lj)
        # Smooth (switched) Lennard-Jones truncation for the MM systems: a
        # plain cutoff makes the MM energy discontinuous when pairs cross it,
        # which shows up as a drift in NVE tests.  Off by default (OpenMM's
        # createSystem default); switched on by [qmmm] lj_switch=true.
        self.lj_switch = bool(lj_switch)
        # OpenMM PME/Ewald error tolerance for the MM-MM systems (None = OpenMM
        # default 5e-4).  Tighten (1e-6) for force/energy consistency checks and
        # NVE validation: the default's force error (~0.5 kJ/mol/nm) is the
        # floor of any finite-difference test on a periodic box.
        if ewald_tol is not None:
            ewald_tol = float(ewald_tol)
            if not np.isfinite(ewald_tol) or ewald_tol <= 0.0:
                raise ValueError(
                    f"[qmmm] ewald_tol must be a finite positive tolerance; got {ewald_tol!r}")
        self.ewald_tol = ewald_tol
        if oqp_cfg is None and mol is None:
            raise ValueError("Either 'oqp_cfg' or 'mol' must be provided.")
        if oqp_cfg is not None and mol is not None:
            raise ValueError(
                "'oqp_cfg' and 'mol' are mutually exclusive – provide only one."
            )

        self.positions = positions
        self.topology = topology
        self.forcefield = forcefield
        # Sort the QM selection into ascending (topology) order. The QM geometry
        # handed to the engine is built by iterating self.topology.atoms() filtered
        # by membership (i.e. topology order), so gqm / f_qm / pchg come back in
        # topology order. The force-scatter loops in _assemble_force / _assemble_
        # force_espf and the link-atom host_row projection index those arrays by
        # position in self.qm_atoms; if qm_atoms were given out of order (e.g.
        # "5,2,7") those positions would not match topology order and QM gradients,
        # coupling forces, and link-atom projections would land on the wrong atoms.
        # Sorting makes input order == topology order without changing the QM
        # calculation (the engine already sees the atoms in topology order).
        self.qm_atoms = np.array(sorted(int(i) for i in qm_atoms), dtype=int)
        self.Cutoff = Cutoff
        self.Embedding = _normalize_embedding(Embedding)

        self.use_mol = mol is not None

        if self.use_mol:
            self.mol = mol
            self.oqp_cfg_base = None
        else:
            self.oqp_cfg_base = oqp_cfg
            self.mol = None

        self.op = None

        # Route the *entire* QM-MM electrostatic coupling through ESPF (energy in
        # the embedded SCF + analytic gradient), with OpenMM reduced to pure
        # MM-MM. This gives a finite-difference-exact analytic gradient for both
        # whole-molecule and covalent-boundary QM regions, so it is the default
        # for electrostatic embedding. ("split" selects the legacy scheme that
        # routes QM charges through OpenMM point charges -- kept for reference.)
        self.espf_full = self.Embedding in (
            "espf", "espf_full", "electrostatic")
        if self.mm_damp_mu is not None and not self.espf_full:
            # Only _full_field_potmm / _coupling_forces apply the erf damping;
            # the split and mechanical schemes route the QM-MM electrostatics
            # through OpenMM point charges and would silently ignore it.
            raise ValueError(
                "[qmmm] mm_charge_width (Gaussian-smeared MM charges) is "
                "implemented for the full-ESPF electrostatic embedding only "
                f"(embedding=electrostatic); got embedding={self.Embedding!r}.")

        # QM/MM boundary connectivity: hydrogen link atoms capping any covalent
        # bond that the QM/MM partition cuts.  Empty when the QM region is a set
        # of whole molecules (e.g. a water box), in which case every code path
        # below is a no-op and behaviour is identical to the pre-link-atom
        # driver.
        self.link_atoms = self._detect_link_atoms()

        # Frontier (M1) charge treatment across a covalent QM/MM cut. Default
        # 'none' = full-field, the validated ESPF baseline: ESPF couples the MM
        # potential to QM atomic-charge operators (H += sum_A phi_A Q_A), which
        # already suppresses the spill-out that would justify redistribution in a
        # density-based embedding. 'rcd' (conserve deleted M1 charge + its dipole
        # via virtual midpoint charges), 'rc', and 'z1' are optional refinements.
        # A no-op when the QM region is whole molecules (no link atoms).
        self.frontier_scheme = str(frontier_scheme or 'none').lower()
        if self.link_atoms and self.espf_full:
            note = ("full-field ESPF baseline" if self.frontier_scheme == 'none'
                    else f"'{self.frontier_scheme}' redistribution (optional refinement)")
            print(
                f"[QM/MM] {len(self.link_atoms)} covalent boundary bond(s); "
                f"frontier-charge embedding = {note}."
            )
        self._select_boundary_switching()

        self.mm_systems = self.prepare_mm()

    # --- Internal helpers -------------------------------------------------

    SWSCALE_WHOLE_MOLECULE = SWSCALE_WHOLE_MOLECULE
    SWSCALE_COVALENT_BOUNDARY = SWSCALE_COVALENT_BOUNDARY

    def _select_boundary_switching(self):
        """Pick ESPF_SWSCALE from whether the QM/MM cut is covalent.

        The selection itself lives in :mod:`oqp.library.qmmm_connectivity`
        (OpenMM-free, so it stays unit-testable without the optional MM
        stack); this driver supplies the one fact only it knows -- whether
        link atoms were built -- and reports the outcome.
        """
        covalent = bool(self.link_atoms)
        self.espf_swscale = select_boundary_switching(covalent)
        if self.espf_swscale is not None and covalent:
            print(f"[QM/MM] covalent boundary detected; "
                  f"ESPF_SWSCALE={self.espf_swscale} "
                  f"(whole-molecule default is "
                  f"{self.SWSCALE_WHOLE_MOLECULE}). Set ESPF_SWSCALE to "
                  f"override.")
        return self.espf_swscale

    def _detect_link_atoms(self):
        """Find dangling QM–MM bonds in the topology and build link atoms."""
        # extra particles (virtual sites such as the TIP4P M site) have no
        # element and no bonds, so they can never end a cut bond
        z_by_index = {
            atom.index: (0 if atom.element is None else atom.element.atomic_number)
            for atom in self.topology.atoms()
        }
        bonds = [(b[0].index, b[1].index) for b in self.topology.bonds()]
        return detect_link_atoms(bonds, self.qm_atoms, lambda i: z_by_index[i])

    def _qm_bond_adjacency(self):
        """{qm_index: [bonded qm_index, ...]} from the topology (cached)."""
        adj = getattr(self, "_qm_adj", None)
        if adj is None:
            qm_set = set(int(i) for i in self.qm_atoms)
            adj = {i: [] for i in qm_set}
            for b in self.topology.bonds():
                i, j = int(b[0].index), int(b[1].index)
                if i in qm_set and j in qm_set:
                    adj[i].append(j); adj[j].append(i)
            self._qm_adj = adj
        return adj

    def unwrap_qm(self, get_xyz, box):
        """Positions of the QM atoms with every bonded QM fragment made whole:
        starting from the lowest-index atom of each connected fragment, each
        neighbour is placed at the minimum-image bond vector from the atom it
        was reached from.  ``get_xyz(i)`` returns the raw (possibly wrapped)
        coordinate of atom i and ``box`` the orthorhombic box in the same
        units (None: no imaging, raw coordinates are returned).  A periodic
        frame that stores bonded atoms on opposite sides of the cell would
        otherwise hand the QM code a bond stretched by a box length."""
        out = {}
        if box is None:
            return {int(i): np.asarray(get_xyz(int(i)), dtype=float) for i in self.qm_atoms}
        adj = self._qm_bond_adjacency()
        box = np.asarray(box, dtype=float)
        placed = []
        for root in sorted(int(i) for i in self.qm_atoms):
            if root in out:
                continue
            out[root] = np.asarray(get_xyz(root), dtype=float)
            members = [root]
            stack = [root]
            while stack:
                i = stack.pop()
                for j in adj[i]:
                    if j not in out:
                        out[j] = out[i] + self._min_image(np.asarray(get_xyz(j), dtype=float) - np.asarray(get_xyz(i), dtype=float), box)
                        stack.append(j)
                        members.append(j)
            # Disconnected QM fragments (several QM molecules): each one is
            # whole now, but its root kept the raw wrapped coordinate, so two
            # fragments that neighbour each other across a box face would be
            # handed to the QM code a box length apart.  Translate every
            # fragment after the first by the lattice vector that puts its
            # centroid at the minimum image of the NEAREST already placed
            # fragment (greedy spanning tree: imaging only against the first
            # fragment leaves two later fragments on opposite sides of it a
            # box length apart even when they are neighbours across a face).
            centroid = np.mean([out[k] for k in members], axis=0)
            if not placed:
                placed.append(centroid)
            else:
                best = None
                for c in placed:
                    d = centroid - c
                    dm = self._min_image(d, box)
                    r = float(np.linalg.norm(dm))
                    if best is None or r < best[0]:
                        best = (r, dm - d)
                shift = best[1]
                if np.any(shift != 0.0):
                    for k in members:
                        out[k] = out[k] + shift
                placed.append(centroid + shift)
        return out

    def _qm_xyz_angstrom(self, positions):
        """Unwrapped QM-atom positions (Angstrom) keyed by atom index."""
        box = self._box_lengths_bohr()
        box_ang = None if box is None else np.asarray(box) / self._ANG2BOHR
        return self.unwrap_qm(lambda i: np.asarray(positions[i].value_in_unit(unit.angstrom), dtype=float), box_ang)

    def _link_positions_angstrom(self, positions, qm_xyz=None):
        """Link-atom Cartesian positions (Angstrom) for the given frame.  In a
        periodic box the QM->MM bond vector is taken as the minimum image, so
        a frame whose bonded hosts are wrapped to opposite sides of the cell
        still places the link hydrogen on the short (bonded) image."""
        box = self._box_lengths_bohr()
        box_ang = None if box is None else np.asarray(box) / self._ANG2BOHR
        if qm_xyz is None:
            qm_xyz = self._qm_xyz_angstrom(positions)
        coords = []
        for link in self.link_atoms:
            qm_raw = np.asarray(positions[link.qm_index].value_in_unit(unit.angstrom), dtype=float)
            mm_p = np.asarray(positions[link.mm_index].value_in_unit(unit.angstrom), dtype=float)
            bond = self._min_image(mm_p - qm_raw, box_ang)
            qm_p = qm_xyz[int(link.qm_index)]           # unwrapped host
            coords.append(link_atom_position(qm_p, qm_p + bond, link.g))
        return coords

    def _build_xyz_string(self):
        xyz_atoms = []
        qm_xyz = self._qm_xyz_angstrom(self.positions)
        for atom in self.topology.atoms():
            at_index = atom.index
            if at_index in self.qm_atoms:
                sym = atom.element.symbol
                x, y, z = qm_xyz[at_index]
                xyz_atoms.append(f"{sym} {x:.12f} {y:.12f} {z:.12f}")
        # Cap severed QM–MM bonds with hydrogen link atoms (appended last so the
        # QM-atom ordering above is preserved).
        for pos in self._link_positions_angstrom(self.positions, qm_xyz):
            xyz_atoms.append(f"H {pos[0]:.12f} {pos[1]:.12f} {pos[2]:.12f}")
        return '; '.join(xyz_atoms)

    def _update_mol_positions(self):
        qm_xyz = self._qm_xyz_angstrom(self.positions)
        coords = [list(qm_xyz[atom.index]) for atom in self.topology.atoms()
                  if atom.index in self.qm_atoms]
        # Append hydrogen link atoms capping severed QM–MM bonds.
        for pos in self._link_positions_angstrom(self.positions, qm_xyz):
            coords.append([pos[0], pos[1], pos[2]])
        coords = np.array(coords)
        ang2bohr = 1.8897259886
        # Molecule has no set_atoms2 (this branch was unreachable until the
        # QM/MM optimiser used mol mode); update_system writes the bohr
        # coordinates straight into the Fortran xyz buffer.
        self.mol.update_system((coords * ang2bohr).ravel())

    def forces_qm_openqp(self, potmm=None, potqm=None):

        # PR #205 review (M1a): the separate QM-QM POTQM correction was folded into
        # the embedded energy (and the SCF Fock via OQP::POTQM) but had no matching
        # force term, an energy/force inconsistency. PME POTMM already captures the
        # periodic MM embedding with the QM self-image removed, so zero POTQM here to
        # match the resolution already applied in the NAMD driver.
        if potqm is not None:
            potqm = np.zeros_like(potqm)

        if self.use_mol:
            # ---- Mol mode ------------------------------------------------
            self._update_mol_positions()
            if is_tb_method(str(self.mol.config['input']['method'])):
                # Native AO-based methods need explicit zero POTMM/POTQM
                # records in mechanical QM/MM, but the tight-binding adapter
                # uses ``None`` as its gas-phase/mechanical contract.
                tb_potmm = None if self.Embedding == "mechanical" else potmm
                return self._forces_qm_dftb(self.mol, tb_potmm)
            sp = SinglePoint(self.mol)
            if getattr(self, "_image_warm", False) or getattr(self, "_reuse_orbitals", False):
                # image iteration > 1 (same geometry), or a caller that moves
                # the geometry in small steps (the QM/MM optimiser): keep the
                # converged orbitals as the guess and rebuild the bare
                # one-electron integrals.  At a new geometry the basis is
                # rebuilt first: ECP centres are copied into the native basis
                # only by set_basis, so an ECP would otherwise stay at the
                # previous geometry (an iodide step gave -3e15 Hartree).
                if getattr(self, "_reuse_orbitals", False) and not getattr(self, "_image_warm", False):
                    oqp.library.set_basis(self.mol)
                ints_1e(self.mol)
            else:
                sp._prep_guess()

            self.mol.data["OQP::POTMM"] = potmm
            self.mol.data["OQP::POTQM"] = potqm
            oqp.espf_op_corr(self.mol)
            espf_op_corr = self.mol.data["OQP::ESPF_CORR"]

            basis = self.mol.data.get_basis()
            nat = self.mol.data["natom"]
            nbf = basis["nbf"]
            espf_op_corr_f = unpack_lower_tri_multi(espf_op_corr, nbf, nat)

            if potmm is not None:
                hcore = self.mol.get_hcore()
                hcore_full = unpack_lower_tri_single(hcore, nbf)
                hcore_full += np.einsum("ijk,i->jk", espf_op_corr_f, potmm)
                self.mol.set_hcore(pack_lower_tri_single(hcore_full))

            self._embedded_scf(sp)
            self.eqm = self.mol.get_scf_energy()

            self._native_embedded_energy_gradient(self.mol, sp, potmm, potqm)
            self._sp = sp

        else:
            # ---- Config mode ---------------------------------------------
            if getattr(self, "_image_warm", False) and getattr(self, "op", None) is not None:
                # image iteration > 1: same geometry, keep the converged
                # orbitals and only rebuild the bare one-electron integrals
                ints_1e(self.op.mol)
            else:
                xyz_atoms = self._build_xyz_string()
                self.oqp_cfg_base["input.system"] = xyz_atoms
                # one log per run: the first geometry opens it, later ones append
                self.op = OPENQP(self.oqp_cfg_base, True,
                                 append_log=getattr(self, "_log_started", False))
                self._log_started = True
                if is_tb_method(str(self.op.mol.config['input']['method'])):
                    tb_potmm = None if self.Embedding == "mechanical" else potmm
                    return self._forces_qm_dftb(self.op.mol, tb_potmm)
                self.op.sp._prep_guess()

            self.op.mol.data["OQP::POTMM"] = potmm
            self.op.mol.data["OQP::POTQM"] = potqm
            oqp.espf_op_corr(self.op.mol)
            espf_op_corr = self.op.mol.data["OQP::ESPF_CORR"]

            basis = self.op.mol.data.get_basis()
            nat = self.op.mol.data["natom"]
            nbf = basis["nbf"]
            espf_op_corr_f = unpack_lower_tri_multi(espf_op_corr, nbf, nat)

            if potmm is not None:
                hcore = self.op.mol.get_hcore()
                hcore_full = unpack_lower_tri_single(hcore, nbf)
                hcore_full += np.einsum("ijk,i->jk", espf_op_corr_f, potmm)
                self.op.mol.set_hcore(pack_lower_tri_single(hcore_full))

            self._embedded_scf(self.op.sp)
            self.eqm = self.op.mol.get_scf_energy()

            self._native_embedded_energy_gradient(self.op.mol, self.op.sp, potmm, potqm)
            self.op.mol.save_data()

        return self.eqm, self.gqm, self.pchg_qm

    @staticmethod
    def _embedded_scf(sp):
        """Embedded SCF through the robustness ladder (primary converger, then
        SOSCF/TRAH escalation from the current orbitals); stop if it still
        does not converge rather than assemble forces on a partial SCF."""
        if not sp._run_scf():
            raise RuntimeError(
                "QM/MM: the embedded SCF did not converge (primary converger and "
                "the SOSCF/TRAH escalation).  Raise [scf] maxit, loosen [scf] "
                "conv, or check the QM/MM contacts.")

    def _native_embedded_energy_gradient(self, mol, sp, potmm, potqm):
        """Post-SCF part of the native (AO-based) embedded QM step, shared by
        the mol and config modes: ESP charges, QM-MM energy bookkeeping
        (full-ESPF: + sum_A Z_A phi_A; split: - (q - Z).phi) and the analytic
        gradient with the ESPF terms, in OpenMM units.  Sets self.pchg_qm,
        self.eqm and self.gqm."""
        oqp.form_esp_charges(mol)
        self.pchg_qm = mol.data["OQP::partial_charges"]

        # The embedded SCF contains only the electronic QM-MM coupling
        # (dEqm/dphi_A = -Q_A, verified by finite differences). In the
        # full-ESPF scheme OpenMM carries no QM charge, so add the
        # nuclear-MM interaction sum_A Z_A phi_A to complete the QM-MM
        # electrostatic energy; its field derivative Z_A dphi/dx together
        # with the electronic response gives the net-charge coupling force
        # already supplied by the analytic coupling term.
        if self.espf_full and potmm is not None:
            self.eqm += float(
                np.dot(mol.get_atoms2("charge"), potmm)
            )

        if potqm is not None and potmm is not None:
            potmm -= np.einsum(
                "ij,j->i", potqm,
                self.pchg_qm - mol.get_atoms2("charge")
            )
            mol.data["OQP::POTMM"] = potmm

        # In the full-ESPF scheme the QM-MM coupling lives entirely in the
        # embedded SCF energy (and OpenMM carries no QM charges), so there is
        # no double count to remove. The split scheme subtracts it here
        # because OpenMM re-adds the coupling via the QM point charges.
        if potmm is not None and not self.espf_full:
            self.eqm -= np.dot(
                self.pchg_qm - mol.get_atoms2("charge"), potmm
            )

        # --- Gradients: pure QM + ESPF contribution -----------------------
        gradient = Gradient(mol)
        if gradient.method == 'hf':
            # Use the common wrapper: an active petite-list build leaves a
            # skeleton in the native buffer, and Gradient.gradient()
            # reconstructs it in the correct frame before returning and
            # writing the projected result back. Reading get_grad()
            # directly here used to bypass both operations.
            gqm = np.asarray(gradient.gradient(), dtype=float).reshape(
                (1, mol.get_atoms2("natom"), 3))
            oqp.grad_esp_qmmm(mol)
            # OQP::ESPF_GRAD is declared Fortran (3, natom) but its flat
            # buffer is atom-major (a0x,a0y,a0z,a1x,...), matching the QM
            # gradient. Reshape the flat buffer to (natom, 3). Adding the
            # (3, natom) view directly only works when natom == 3 (a square
            # coincidence), which is why non-3-atom QM regions - e.g.
            # link-atom-capped fragments - previously broke.
            natom_qm = mol.get_atoms2("natom")
            esp_grad = np.asarray(
                mol.data["OQP::ESPF_GRAD"]
            ).reshape(natom_qm, 3)
            gqm += esp_grad
            # --- Unit conversion to OpenMM conventions ------------------------
            self.eqm *= 2625.499639 * unit.kilojoule_per_mole
            self.gqm = gqm[0]*49614.75  # Hartree/bohr -> kJ/mol/nm (PR #205 review M1b)
        if gradient.method == 'tdhf':
            energies = sp.excitation([self.eqm])
            grads = np.zeros(( gradient.nstate + 1,  gradient.natom, 3))
            for i in gradient.grads:
                target = (public_state_label(gradient.mol.config, i)
                          if is_mrsf(gradient.mol.config) else 'Root %s' % i)
                dump_log(gradient.mol, title='PyOQP: Gradient of %s' % target)
                gradient.mol.data.set_tdhf_target(i)
                gradient.zvec_func[gradient.td](gradient.mol)

                # check convergence
                z_flag = gradient.mol.mol_energy.Z_Vector_converged

                if not z_flag:
                    dump_log(gradient.mol, title='PyOQP: TD Z-vector is not converged', section='end')

                    if gradient.exception is True:
                        raise ZVnotConverged()
                    else:
                        exit()

                gradient.grad_func[gradient.td](gradient.mol)
                gqm = gradient.mol.get_grad().reshape((gradient.natom, 3))
                # This state-by-state QM/MM path cannot call the common
                # wrapper as a batch because ESPF_GRAD is state-specific.
                # Apply the same reconstruction here before the force is
                # assembled, and keep the public native buffer consistent.
                gqm = np.asarray(
                    gradient.mol.symmetrize_gradient(gqm), dtype=float
                ).reshape((gradient.natom, 3))
                gradient.mol.set_grad(gqm)
                oqp.grad_esp_qmmm_excited(mol)
                # ESPF_GRAD flat buffer is atom-major; reshape to (natom, 3).
                gqm += np.asarray(
                    mol.data["OQP::ESPF_GRAD"]
                ).reshape(gradient.natom, 3)
                grads[i] = gqm.copy()
                self.gqm = gqm*49614.75  # Hartree/bohr -> kJ/mol/nm (PR #205 review M1b)
                self.eqm = energies[i] * 2625.499639 * unit.kilojoule_per_mole
                # The Z-vector step has replaced OQP::partial_charges by the
                # RELAXED ESPF charges of this state: they drive the classical
                # coupling forces and, in a periodic box, the QM-image
                # self-consistency loop in compute_force (which then converges
                # the field to the propagated state, not to the reference
                # density).  Publish them explicitly rather than through the
                # live view taken above.
                self.pchg_qm = np.array(mol.data["OQP::partial_charges"], dtype=float)

    def _forces_qm_dftb(self, mol, potmm):
        """QM energy/gradient/charges for the TB backends (method=dftb/xtb).

        Electrostatic embedding contract (see oqp.library.openqp_dftb): setting
        ``mol.dftb_external_potential`` = POTMM (Hartree/e, one value per QM
        centre INCLUDING hydrogen link-atom caps, in QM-geometry order) makes
        the library fold the potential directly into the SCC Hamiltonian, so

          * the returned state energy is the COMPLETE embedded QM energy: it
            already contains E_ext = sum_A q_A phi_A with q_A the NET atomic
            charge (valence cores + electrons). Unlike the native ESPF path
            there is NO separate nuclear ``+ sum_A Z_A phi_A`` term to add
            (dE/dphi_A = +q_A, verified by finite differences);
          * the returned analytic gradient is d(E_embedded)/dR_QM at FIXED
            potential values -- the charge-response (Pulay-type) coupling term
            is already inside it, so the native ``oqp.grad_esp_qmmm(_excited)``
            call and the ``OQP::ESPF_GRAD`` addition are SKIPPED;
          * the classical dphi/dR forces (MM field gradients acting on the QM
            charges, and the reaction on the MM charges) remain the driver's
            job: the backend-agnostic ``_coupling_forces`` reads
            ``OQP::partial_charges``, which the adapter publishes (relaxed net
            atomic charges of the active state) after the gradient call -- the
            DFTB analog of the native ``form_esp_charges`` step.

        Only the full-ESPF electrostatic scheme (``espf_full``) and mechanical
        embedding (``potmm is None``) are supported: the legacy 'split' scheme
        routes QM point charges through OpenMM, which would double-count the
        coupling already inside the embedded DFTB energy.
        """
        if potmm is not None and not self.espf_full:
            raise NotImplementedError(
                "method=dftb/xtb QM/MM supports Embedding='electrostatic'/'espf' "
                "(full-ESPF scheme) or 'mechanical'; the legacy 'split' scheme "
                "would double-count the coupling embedded in the DFTB energy."
            )

        mol.dftb_external_potential = (
            None if potmm is None else np.asarray(potmm, dtype=float)
        )

        gradient = Gradient(mol)
        grads = gradient.gradient()   # dispatches to the openqp-dftb adapter
        # Active state: same [properties] grad selection the native TDHF branch
        # uses (ground state runs carry grad=[0]).
        active = max(gradient.grads)
        eqm_ha = float(mol.energies[active])
        self.gqm = np.asarray(grads[active], dtype=float) * 49614.75  # Ha/bohr -> kJ/mol/nm
        # Relaxed net atomic charges of the active state, published by the
        # adapter after the gradient call (authoritative for the classical
        # coupling forces).
        self.pchg_qm = np.array(mol.data["OQP::partial_charges"], dtype=float)
        self.eqm = eqm_ha * 2625.499639 * unit.kilojoule_per_mole
        if not self.use_mol:
            mol.save_data()
        return self.eqm, self.gqm, self.pchg_qm

    def _get_mol(self):
        if self.use_mol:
            return self.mol
        else:
            return self.op.mol

    def forces_mm(self, pchg_qm):
        system = self.mm_systems["sys0"]
        simulation = self.mm_systems["sim0"]

        forces = { force.__class__.__name__ : force for force in system.getForces() }
        nonbonded = forces['NonbondedForce']

        if self.Embedding == "mechanical":
            # Mechanical embedding: the QM atoms keep their FIXED force-field
            # charges in the MM electrostatics (the system is used as built by
            # prepare_mm, intra-QM pairs excluded).  Injecting the fitted ESPF
            # charges of the gas-phase QM density here would make E_MM depend
            # on the geometry through q(R) while OpenMM differentiates it at
            # fixed charges, so the force would not be the derivative of the
            # energy; ``pchg_qm`` is therefore ignored on this path.
            state = simulation.context.getState(getEnergy=True, getForces=True)
            return state.getPotentialEnergy(), state.getForces(asNumpy=True)

        if self.espf_full:
            # Pure MM-MM: QM atoms carry no charge; all QM-MM electrostatics are
            # handled analytically by ESPF + the coupling force. vdW and bonded
            # terms (incl. boundary-spanning) are retained.
            for iatom in self.qm_atoms:
                charge, sigma, epsilon = nonbonded.getParticleParameters(iatom)
                nonbonded.setParticleParameters(
                    iatom, 0.0 * unit.elementary_charge, sigma, epsilon)
            nonbonded.updateParametersInContext(simulation.context)
            state = simulation.context.getState(getEnergy=True, getForces=True)
            return state.getPotentialEnergy(), state.getForces(asNumpy=True)

        for k, iatom in enumerate(self.qm_atoms):
            charge, sigma, epsilon = nonbonded.getParticleParameters(iatom)
            charge = pchg_qm[k]*unit.elementary_charge
            nonbonded.setParticleParameters(iatom, charge, sigma, epsilon)

        for i in range(nonbonded.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = nonbonded.getExceptionParameters(i)
            # A QM-MM bonded exclusion (1-2/1-3, chargeProd==0) must stay a full
            # exclusion: turning it into a charged pair would change the set of
            # non-excluded exceptions (OpenMM forbids that in
            # updateParametersInContext) and would double-count QM->MM
            # electrostatics already handled by the ESPF embedding. Only genuine
            # scaled (1-4) exceptions carry a non-zero charge product.
            if chargeProd.value_in_unit(unit.elementary_charge ** 2) == 0.0:
               continue
            if p1 in self.qm_atoms and p2 not in self.qm_atoms:
               k1 = int(np.where(self.qm_atoms == p1)[0][0])
               charge1 = pchg_qm[k1]*unit.elementary_charge
               charge2, sigma2, epsilon2 = nonbonded.getParticleParameters(p2)
               chargeProd=charge1*charge2
               nonbonded.setExceptionParameters(i,p1,p2,chargeProd,sigma,epsilon)
            elif p2 in self.qm_atoms and p1 not in self.qm_atoms:
               k2 = int(np.where(self.qm_atoms == p2)[0][0])
               charge2 = pchg_qm[k2]*unit.elementary_charge
               charge1, sigma1, epsilon1 = nonbonded.getParticleParameters(p1)
               chargeProd=charge1*charge2
               nonbonded.setExceptionParameters(i,p1,p2,chargeProd,sigma,epsilon)
        nonbonded.updateParametersInContext(simulation.context)

        state = simulation.context.getState(getEnergy=True,getForces=True)
        return state.getPotentialEnergy(), state.getForces(asNumpy=True)

    def _virtual_site_rows(self):
        """Indices of the virtual sites of sys0 (cached per System); empty
        when no OpenMM system is attached."""
        systems = getattr(self, "mm_systems", None)
        sys0 = systems.get("sys0") if isinstance(systems, dict) else None
        if sys0 is None or not hasattr(sys0, "isVirtualSite"):
            return []
        cache = getattr(self, "_vsite_cache", None)
        if cache is None or cache[0] is not sys0:
            self._vsite_cache = (sys0, [i for i in range(sys0.getNumParticles())
                                        if sys0.isVirtualSite(i)])
        return self._vsite_cache[1]

    def _drop_openmm_site_rows(self, forces):
        """OpenMM's State forces already add each virtual site's force to the
        particles that place it and also keep a copy on the site row.  Zero
        that copy, so that ``_fold_virtual_site_forces`` moves only the forces
        this driver adds on a site (the ESPF coupling force on a charged TIP4P
        M site) and nothing is counted twice."""
        idx = self._virtual_site_rows()
        if idx:
            forces[idx] = 0.0
        return forces

    def _fold_virtual_site_forces(self, forces):
        """Fold every force left on a virtual site onto its parents (site rows
        end at zero).  The real particles then carry the whole derivative: the
        optimiser takes their rows as the gradient, and an integrator that
        distributes site forces itself (the QM/MM MD system) adds nothing
        twice.  No-op without virtual sites."""
        idx = self._virtual_site_rows()
        if not idx:
            return forces
        return fold_virtual_site_forces(self.mm_systems["sys0"], forces,
                                        positions=self.positions, sites=idx)

    def compute_force(self, positions, topology, mm_systems, qm_atoms):
        """QM/MM energy and the force on every particle.  Forces on virtual
        sites are folded onto the particles that place them (site rows zero)."""
        self.positions = positions
        self.topology = topology
        self.mm_systems = mm_systems
        # Re-normalize to ascending (topology) order on every call: callers
        # (e.g. QMMM_MD) pass their own config-order qm_atoms here, which would
        # otherwise overwrite the sorted copy set in __init__ and reintroduce the
        # force mis-scatter this fix prevents. See __init__ for the rationale.
        self.qm_atoms = np.array(sorted(int(i) for i in qm_atoms), dtype=int)

        potmm = potqm = None
        if self.Embedding in ("electrostatic", "split") or self.espf_full:
            potmm, potqm = self.electrostatic_potential()
        elif self.Embedding == "mechanical":
            # Mechanical embedding (PR #274): the QM subsystem sees no MM
            # field, so the SCF is gas-phase and the QM-MM electrostatics is
            # left to OpenMM with the fixed force-field charges of the QM
            # atoms (forces_mm ignores the ESP charges on this path, so the
            # MM energy is differentiated at the charges it was built with).
            # The embedding
            # arrays must still exist and be ZERO rather than absent:
            # scf.F90 calls add_potqm_contributions on every SCF iteration
            # whenever qmmm_flag is set and aborts on a missing OQP::POTQM
            # record ("Record `OQP::POTQM` not found!"); grad_esp_qmmm needs
            # OQP::POTMM the same way.  A zero field reproduces gas-phase QM
            # exactly (adds 0 to hcore, 0 to the ESPF gradient, 0 to eqm).
            potmm, potqm = self._zero_embedding()
        else:   # guarded in __init__; keep a local invariant
            raise ValueError(f"Unknown QM/MM embedding '{self.Embedding}'")

        if (self.espf_full and self._ewald() is not None
                and os.environ.get("OQP_EWALD_NO_IMAGE", "").strip() not in ("1", "on")):
            # Periodic full-ESPF: the QM charges also interact with their own
            # periodic images (paper eq 8).  Solve for the QM charges self-
            # consistently in the total field  phi_eff = Phi^MM + psi_img q :
            # with E = E_QM[phi_eff] + Z.phi_eff - 1/2 q psi_img q the charge
            # response terms cancel at self-consistency, so the force needs only
            # the explicit derivatives (embedded QM gradient at fixed phi_eff,
            # classical Ewald coupling force, image force at fixed q).
            psi_img, dpsi = self._ewald().qm_image_matrix(self._qm_center_positions_bohr())
            n = len(potmm)
            q_prev = (self._q_prev if getattr(self, "_q_prev", None) is not None
                      and len(self._q_prev) == n else np.zeros(n))
            converged = False
            delta, it = float("inf"), -1
            self._image_warm = False
            # Reference-density charges converge to IMAGE_TOL; the relaxed
            # charges of a TDHF/MRSF target state carry the Z-vector residual
            # and converge to IMAGE_TOL_ACTIVE (as in NAMD_QMMM).
            tol = (self.IMAGE_TOL_ACTIVE if self._image_uses_relaxed_charges()
                   else self.IMAGE_TOL)
            for it in range(int(self.IMAGE_MAXITER)):
                phi_eff = np.asarray(potmm, dtype=float) + psi_img @ q_prev
                try:
                    eqm, gqm, pchg_qm = self.forces_qm_openqp(potmm=phi_eff.copy(), potqm=potqm)
                finally:
                    self._image_warm = True      # later iterations reuse the orbitals
                q_new = np.array(pchg_qm, dtype=float)
                delta = float(np.abs(q_new - q_prev).max())
                if delta < tol:
                    converged = True
                    break
                q_prev = 0.5 * (q_new + q_prev) if it > 6 else q_new
            if not converged:
                raise RuntimeError(
                    f"Periodic ESPF QM/MM: the QM-image charge self-consistency "
                    f"did not converge in {it + 1} iterations "
                    f"(max |dq| = {delta:.2e} e > {tol:.0e}); the "
                    "energy/force would be inconsistent.  Tighten the SCF "
                    "convergence or check the QM/MM contacts.")
            self._image_warm = False
            self._q_prev = q_new.copy()
            self._image_iterations = it + 1
            e_img = 0.5 * float(q_new @ psi_img @ q_new)                    # Hartree
            f_img = -np.einsum("a,b,abc->ac", q_new, q_new, dpsi)           # Hartree/bohr
            return self._assemble_force_espf(
                eqm, q_new, f_qm_extra=f_img,
                e_extra=-e_img * 2625.499639 * unit.kilojoule_per_mole)

        eqm, gqm, pchg_qm = self.forces_qm_openqp(potmm=potmm, potqm=potqm)
        nqm = len(self.qm_atoms)

        if self.espf_full:
            return self._assemble_force_espf(eqm, pchg_qm)

        # Fold each link atom's ESP charge onto its QM host so the charge the QM
        # region presents to the MM electrostatics is conserved (link atoms are
        # not MM particles).  No-op when there are no link atoms.
        pchg_mm = np.array(pchg_qm, dtype=float).copy()
        for a, link in enumerate(self.link_atoms):
            pchg_mm[link.host_row] += pchg_mm[nqm + a]

        emm, gmm = self.forces_mm(pchg_mm)

        total_energy = eqm + emm
        total_forces = self._drop_openmm_site_rows(gmm.copy())

        # Redistribute link-atom gradients onto their real host atoms by the
        # chain rule of the scaled capping position R_L = R_QM + g(R_MM-R_QM).
        # The QM-host share is folded into the QM gradient row (distributed
        # below); the MM-host share is applied directly to the MM host atom.
        for a, link in enumerate(self.link_atoms):
            g_link = self.gqm[nqm + a]
            self.gqm[link.host_row] = self.gqm[link.host_row] + (1.0 - link.g) * g_link
            total_forces[link.mm_index] = total_forces[link.mm_index] - link.g * g_link

        for k, i in enumerate(self.qm_atoms):
            total_forces[i] = total_forces[i] - self.gqm[k]
        cmm=np.sum(total_forces,axis=0)
        for i in range(len(total_forces)):
            total_forces[i]-=cmm/float(len(total_forces))

        return total_energy, self._fold_virtual_site_forces(total_forces)

    def _assemble_force_espf(self, eqm, pchg_qm, f_qm_extra=None, e_extra=None):
        """Assemble the total force in the full-ESPF scheme:
          F = F_MM(pure)  -  grad_QM(HF+ESPF charge-fluctuation)  +  F_coupling
        where F_coupling is the analytic QM-charge <-> MM-charge Coulomb force
        (field-fluctuation term), applied to both QM and MM atoms.
        ``f_qm_extra`` (Hartree/bohr, per QM centre incl. link atoms) and
        ``e_extra`` (OpenMM energy) carry the periodic QM-image term.
        """
        FCONV = 49614.75  # Hartree/bohr -> kJ/mol/nm
        nqm = len(self.qm_atoms)

        emm, gmm = self.forces_mm(pchg_qm)   # pure MM-MM (QM charges zeroed)
        total_energy = eqm + emm
        if e_extra is not None:
            total_energy = total_energy + e_extra
        total_forces = self._drop_openmm_site_rows(gmm.copy())

        f_qm, f_mm, mm_idx = self._coupling_forces(np.asarray(pchg_qm, dtype=float))
        if f_qm_extra is not None:
            f_qm = f_qm + np.asarray(f_qm_extra, dtype=float)
        f_qm = f_qm * FCONV
        f_mm = f_mm * FCONV

        # Project link-atom contributions (QM gradient and coupling force) onto
        # the real host atoms.
        for a, link in enumerate(self.link_atoms):
            g_link = self.gqm[nqm + a]
            self.gqm[link.host_row] = self.gqm[link.host_row] + (1.0 - link.g) * g_link
            total_forces[link.mm_index] = total_forces[link.mm_index] - link.g * g_link
            fl = f_qm[nqm + a]
            f_qm[link.host_row] = f_qm[link.host_row] + (1.0 - link.g) * fl
            total_forces[link.mm_index] = total_forces[link.mm_index] + link.g * fl

        for k, i in enumerate(self.qm_atoms):
            total_forces[i] = total_forces[i] - self.gqm[k] + f_qm[k]
        for j, m in enumerate(mm_idx):
            total_forces[m] = total_forces[m] + f_mm[j]

        # the coupling force on a charged virtual site (TIP4P M) goes to its parents
        return total_energy, self._fold_virtual_site_forces(total_forces)


    def prepare_mm(self):
        positions=self.positions
        topology=self.topology
        forcefield=self.forcefield
        qm_atoms =self.qm_atoms
        qm_set = set(int(i) for i in qm_atoms)
        Cutoff=self.Cutoff
        nb_cutoff = _periodic_nonbonded_cutoff(topology, Cutoff)

        _ew = {} if (self.ewald_tol is None or not is_periodic_method(Cutoff)) else {
            "ewaldErrorTolerance": float(self.ewald_tol)}
        system=forcefield.createSystem(
            topology, nonbondedMethod=Cutoff, nonbondedCutoff=nb_cutoff,
            constraints=None, rigidWater=False, **_ew)
        nonbonded = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
        if self.lj_switch and Cutoff is not app.NoCutoff:
            nonbonded.setUseSwitchingFunction(True)
            nonbonded.setSwitchingDistance(0.85 * nb_cutoff)
        if self.h_lj:
            sig_h = 2.0 * 0.02245 / (2.0 ** (1.0 / 6.0)) * unit.nanometer      # CHARMM HT Rmin/2 = 0.2245 A
            eps_h = 0.046 * 4.184 * unit.kilojoule_per_mole
            n_set = 0
            for atom in topology.atoms():
                if atom.element is None or atom.element.atomic_number != 1:
                    continue
                if int(atom.index) in qm_set:
                    continue          # QM hydrogens keep their own LJ parameters
                q, sig, eps = nonbonded.getParticleParameters(atom.index)
                if eps.value_in_unit(unit.kilojoule_per_mole) == 0.0:
                    nonbonded.setParticleParameters(atom.index, q, sig_h, eps_h)
                    n_set += 1
            print(f"[QM/MM] h_lj: Lennard-Jones parameters assigned to {n_set} MM hydrogen(s) that had none")

        n_14 = 0
        for i in range(nonbonded.getNumExceptions()):
            p1, p2, chgProd, sigma, epsilon = nonbonded.getExceptionParameters(i)
            if (int(p1) in qm_set) or (int(p2) in qm_set):
               if self.espf_full and chgProd.value_in_unit(unit.elementary_charge ** 2) != 0.0:
                   # Full ESPF routes the ENTIRE QM-MM electrostatics through
                   # the embedded SCF + coupling force, and forces_mm zeroes
                   # the QM particle charges -- but OpenMM keeps the 1-4
                   # exception charge products independently of the particle
                   # charges, so a scaled QM-MM 1-4 Coulomb pair across a
                   # covalent boundary would stay in the "pure MM" energy and
                   # force on top of the ESPF term.  Drop it here; the LJ part
                   # of QM-involving exceptions is handled as before.
                   chgProd = 0.0 * unit.elementary_charge ** 2
                   n_14 += 1
               nonbonded.setExceptionParameters(i, p1, p2, chgProd, 0.0, 0.0)
        if n_14:
            print(f"[QM/MM] full ESPF: {n_14} QM-MM 1-4 exception charge product(s) removed "
                  "from the MM system (the QM-MM electrostatics is carried by ESPF)")

        for p1 in qm_atoms:
           for p2 in qm_atoms:
              if p1 != p2:
                 nonbonded.addException(p1,p2,0,0,0,replace=True)

        int0=mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
        simulation=app.Simulation(topology, system, int0)
        simulation.context.setPositions(positions)

        #Deactivate QM-QM interactions
        for f in system.getForces():
        #Deactivating non-bonded terms between QM atoms (but keep QM-MM)
           if isinstance(f, mm.NonbondedForce):
              pass
           elif isinstance(f, mm.CustomNonbondedForce):
              for p1 in qm_atoms:
                for p2 in qm_atoms:
                   f.addExclusion(p1,p2)
              f.updateParametersInContext(simulation.context)
        #Deactivating QM bonded terms
           elif isinstance(f, mm.HarmonicBondForce):
              for i in range(f.getNumBonds()):
                 p1, p2, length, k = f.getBondParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms)
                 if exclude: f.setBondParameters(i, p1, p2, length, 0)
              f.updateParametersInContext(simulation.context)
           elif isinstance(f, mm.CustomBondForce):
              for i in range(f.getNumBonds()):
                 p1, p2, parameters = f.getBondParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms)
                 if exclude: f.setBondParameters(i, p1, p2, (0,0))
        #Deactivating QM angle terms
           elif isinstance(f, mm.HarmonicAngleForce):
              for i in range(f.getNumAngles()):
                 p1, p2, p3, angle, k = f.getAngleParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms and p3 in qm_atoms)
                 if exclude: f.setAngleParameters(i, p1, p2, p3, angle, 0)
              f.updateParametersInContext(simulation.context)
           elif isinstance(f, mm.CustomAngleForce):
              for i in range(f.getNumAngles()):
                 p1, p2, p3, parameters = f.getAngleParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms and p3 in qm_atoms)
                 if exclude: f.setAngleParameters(i, p1, p2, p3, (0, 0))
              f.updateParametersInContext(simulation.context)
        #Deactivating QM torsion terms
           elif isinstance(f, mm.PeriodicTorsionForce):
              for i in range(f.getNumTorsions()):
                 p1, p2, p3, p4, periodicity, phase, k = f.getTorsionParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms and p3 in qm_atoms and p4 in qm_atoms)
                 if exclude: f.setTorsionParameters(i, p1, p2, p3, p4, periodicity, phase, 0)
              f.updateParametersInContext(simulation.context)
           elif isinstance(f, mm.CustomTorsionForce):
              for i in range(f.getNumTorsions()):
                 p1, p2, p3, p4, parameters = f.getTorsionParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms and p3 in qm_atoms and p4 in qm_atoms)
                 if exclude: f.setTorsionParameters(i, p1, p2, p3, p4, (0,0,0))
              f.updateParametersInContext(simulation.context)
           elif isinstance(f, (mm.CMAPTorsionForce)):
              for i in range(f.getNumTorsions()):
                 cmap, p1, p2, p3, p4, q1, q2, q3, q4 = f.getTorsionParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms and p3 in qm_atoms and p4 in qm_atoms)
                 exclude = exclude and (q1 in qm_atoms and q2 in qm_atoms and q3 in qm_atoms and q4 in qm_atoms)
                 if exclude: f.setMapParameters(i,cmap.size,0)
              if f.getNumTorsions() != 0: f.updateParametersInContext(simulation.context)
        #Exception, unless CMMotionRemover
           else:
              if not isinstance(f, mm.CMMotionRemover): exit(f"Force not found")

        if is_periodic_method(Cutoff):
           sysew=forcefield.createSystem(
               topology, nonbondedMethod=app.Ewald, nonbondedCutoff=nb_cutoff,
               constraints=None, rigidWater=False, **_ew)
           intew=mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
           simew=app.Simulation(topology, sysew, intew)
           simew.context.setPositions(positions)

           sysor=forcefield.createSystem(topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
           intor=mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
           simor=app.Simulation(topology, sysor, intor)
           simor.context.setPositions(positions)
        else:
           sysew = simew = sysor = simor = None
        return {
         "sys0": system,
         "sim0": simulation,
         "sysew": sysew,
         "simew": simew,
         "sysor": sysor,
         "simor": simor,
        }


    def _zero_embedding(self):
        """Zero MM potential over every QM centre (real QM atoms + link
        atoms), sized like the arrays the ESPF path builds (nqm + nlink)."""
        n = len(self.qm_atoms) + len(self.link_atoms)
        return np.zeros(n), np.zeros((n, n))

    def _pad_potential_for_link_atoms(self, potmm, potqm):
       """Extend the ESPF embedding arrays to cover hydrogen link atoms.

       Link atoms are additional QM centres appended after the real QM atoms.
       They are capping hydrogens rather than physical atoms, so they are not
       embedded in the MM electrostatic field (their MM potential is taken as
       zero) and carry no periodic QM–QM self-image term.  This keeps the
       array dimensions consistent with the QM geometry (natom = nqm + nlink)
       while leaving the whole-molecule (no-link) case untouched.
       """
       nlink = len(self.link_atoms)
       if nlink == 0:
           return potmm, potqm
       nqm = len(self.qm_atoms)
       potmm = np.concatenate([np.asarray(potmm, dtype=float), np.zeros(nlink)])
       padded = np.zeros((nqm + nlink, nqm + nlink), dtype=float)
       padded[:nqm, :nqm] = potqm
       return potmm, padded

    # ------------------------------------------------------------------
    #  Full-ESPF (non-split) QM-MM electrostatics
    # ------------------------------------------------------------------
    _ANG2BOHR = 1.8897259886

    def _qm_center_positions_bohr(self):
        """Cartesian positions (bohr) of every QM centre (real QM atoms in
        topology order, then hydrogen link atoms), matching the QM geometry
        order used to build the QM system and the POTMM array."""
        qm_xyz = self._qm_xyz_angstrom(self.positions)
        coords = [[c * self._ANG2BOHR for c in qm_xyz[atom.index]]
                  for atom in self.topology.atoms() if atom.index in self.qm_atoms]
        for pos in self._link_positions_angstrom(self.positions, qm_xyz):
            coords.append([c * self._ANG2BOHR for c in pos])
        return np.asarray(coords, dtype=float)

    def _frontier_hosts(self):
        """[(m1_abs_idx, [m2_abs_idx, ...]), ...] : each unique MM host atom
        (the MM end of a bond the QM/MM boundary cuts) with its MM neighbours
        M2, from the topology. Empty for a whole-molecule QM region."""
        if not self.link_atoms:
            return []
        qm_set = set(int(i) for i in self.qm_atoms)
        nbrs = {}
        for b in self.topology.bonds():
            i, j = int(b[0].index), int(b[1].index)
            nbrs.setdefault(i, set()).add(j)
            nbrs.setdefault(j, set()).add(i)
        hosts = {}
        for link in self.link_atoms:
            m1 = int(link.mm_index)
            if m1 in hosts:
                continue
            hosts[m1] = sorted(n for n in nbrs.get(m1, ()) if n not in qm_set)
        return [(m1, hosts[m1]) for m1 in sorted(hosts)]

    def _embedding_sites(self):
        """The MM embedding point-charge set the QM density sees, with the
        frontier-host (M1) charges redistributed per ``self.frontier_scheme``
        (RCD by default). Returns (charges (S,), positions (S,3) bohr, scatter),
        where scatter maps a per-site force back onto the real MM atoms (identity
        for a real atom, 0.5/0.5 for a virtual midpoint charge). With no link
        atoms this is the identity -> whole-molecule QM/MM is unchanged."""
        mmq, mm_xyz, mm_idx = self._mm_charges_positions_bohr()
        q_of = {int(mm_idx[k]): float(mmq[k]) for k in range(len(mm_idx))}
        deleted, delta_q, virtuals = redistribute_frontier_charges(
            self._frontier_hosts(), lambda a: q_of.get(a, 0.0),
            self.frontier_scheme,
        )
        box = self._box_lengths_bohr()
        return assemble_embedding_sites(
            mm_idx, mmq, mm_xyz, deleted, delta_q, virtuals,
            min_image=None if box is None else (lambda d: self._min_image(d, box)))

    def _image_uses_relaxed_charges(self):
        """True when the QM step publishes the RELAXED ESPF charges of a
        response (TDHF/MRSF) target state, which carry the Z-vector residual
        and converge the QM-image loop to IMAGE_TOL_ACTIVE instead of the
        reference-density IMAGE_TOL."""
        if self.use_mol:
            method = self.mol.config.get("input", {}).get("method", "hf")
        else:
            method = (self.oqp_cfg_base or {}).get("input.method", "hf")
        return str(method).strip().lower() == "tdhf"

    def _is_periodic(self):
        """True when the MM nonbonded method is a periodic one (PME, Ewald,
        CutoffPeriodic).  NoCutoff and CutoffNonPeriodic are treated as a
        finite cluster even if the topology carries box vectors."""
        return is_periodic_method(self.Cutoff)

    #: QM-image charge self-consistency loop (periodic full-ESPF): iteration
    #: cap and convergence threshold on the ESPF charges (e).
    IMAGE_MAXITER = 50
    IMAGE_TOL = 1e-7
    #: Active-state refinement in NAMD (each iteration costs a Z-vector
    #: gradient).  The relaxed charges carry the Z-vector residual (1e-5 to
    #: 1e-4 e at the default Z-vector convergence), so the field is taken as
    #: self-consistent at 1e-4 e; the state energy is then stable to ~1e-7 Ha
    #: and the loop typically needs 3 gradient evaluations.
    IMAGE_TOL_ACTIVE = 1e-4
    IMAGE_MAXITER_ACTIVE = 20
    IMAGE_ETOL_ACTIVE = 1e-7      # Hartree: energy-stagnation acceptance (with |dq| < 10 IMAGE_TOL_ACTIVE)

    def _box_lengths_bohr(self):
        """Orthorhombic periodic box lengths (bohr), or None when the QM/MM
        electrostatics are non-periodic (NoCutoff / CutoffNonPeriodic). Used
        for the Ewald / minimum-image QM-MM electrostatics under PBC."""
        if not self._is_periodic():
            return None
        vecs = self.topology.getPeriodicBoxVectors()
        if vecs is None:
            raise ValueError(
                f"[qmmm] cutoff={self.Cutoff} is periodic but the PDB topology "
                "carries no box vectors (CRYST1 record).")
        box = np.array([[c.value_in_unit(unit.angstrom) for c in v] for v in vecs])
        if np.abs(box - np.diag(np.diag(box))).max() > 1e-8:
            raise NotImplementedError(
                "Periodic ESPF QM/MM supports orthorhombic boxes only; the PDB "
                "box vectors are not diagonal.")
        return np.diag(box) * self._ANG2BOHR

    def _min_image(self, d, box):
        """Minimum-image displacement(s) for an (N,3) array under an
        orthorhombic box (bohr). No-op when box is None."""
        if box is None:
            return d
        return d - box * np.round(d / box)

    def _mm_charges_positions_bohr(self):
        """MM (non-QM) force-field charges (e), positions (bohr) and absolute
        atom indices. These carry the classical electrostatics the QM density
        is embedded in."""
        nb = next(f for f in self.mm_systems["sys0"].getForces()
                  if isinstance(f, mm.NonbondedForce))
        qm_set = set(int(i) for i in self.qm_atoms)
        q, xyz, idx = [], [], []
        for i in range(nb.getNumParticles()):
            if i in qm_set:
                continue
            charge, _, _ = nb.getParticleParameters(i)
            p = self.positions[i].value_in_unit(unit.angstrom)
            q.append(charge.value_in_unit(unit.elementary_charge))
            xyz.append([c * self._ANG2BOHR for c in p])
            idx.append(i)
        return np.asarray(q), np.asarray(xyz, dtype=float), np.asarray(idx, dtype=int)

    def _ewald(self):
        """Ewald summation object for the current orthorhombic box, or None
        when the QM/MM electrostatics are non-periodic (NoCutoff)."""
        box = self._box_lengths_bohr()
        if box is None:
            return None
        cached = getattr(self, "_ewald_obj", None)
        if cached is None or not np.allclose(cached.box, box):
            self._ewald_obj = EwaldQMMM(box)
            self._ewald_obj.check_damping(self.mm_damp_mu)
        return self._ewald_obj

    def _full_field_potmm(self):
        """MM electrostatic potential at every QM centre from the (frontier-
        redistributed) embedding charge set (Hartree/e): a direct Coulomb sum
        for a non-periodic system, the Ewald lattice sum (all images, tin-foil
        boundary) for a periodic one."""
        qm_xyz = self._qm_center_positions_bohr()
        q_s, xyz_s, _ = self._embedding_sites()
        ew = self._ewald()
        if ew is not None:
            phi, _ = ew.mm_potential(qm_xyz, xyz_s, q_s, mu=self.mm_damp_mu)
            return phi
        box = self._box_lengths_bohr()
        potmm = np.zeros(len(qm_xyz))
        from scipy.special import erf
        for a in range(len(qm_xyz)):
            d = self._min_image(qm_xyz[a] - xyz_s, box)
            r = np.linalg.norm(d, axis=1)
            phi_kernel, _ = smeared_coulomb(r, self.mm_damp_mu)
            potmm[a] = np.sum(q_s * phi_kernel)
        return potmm

    def _coupling_forces(self, pchg):
        """Classical Coulomb force (Hartree/bohr) between the QM ESP charges q_A
        and the (frontier-redistributed) MM embedding charges Q_s -- the field-
        fluctuation term Tr[q dphi/dx]. A force on a virtual midpoint charge is
        chain-ruled onto its real hosts (M1, M2) via the site scatter map.
        Returns (F on QM centres, F on real MM atoms, real MM idx)."""
        qm_xyz = self._qm_center_positions_bohr()
        _, _, mm_idx = self._mm_charges_positions_bohr()
        q_s, xyz_s, scatter = self._embedding_sites()
        box = self._box_lengths_bohr()
        f_qm = np.zeros_like(qm_xyz)
        f_site = np.zeros_like(xyz_s)
        ew = self._ewald()
        if ew is not None:
            # Ewald: F_A = -q_A dPhi^MM_A/dr_A ; F_s from the QM charges' images
            _, dphi = ew.mm_potential(qm_xyz, xyz_s, q_s, mu=self.mm_damp_mu)
            f_qm = -np.asarray(pchg, dtype=float)[:, None] * dphi
            f_site = ew.mm_forces(qm_xyz, pchg, xyz_s, q_s, mu=self.mm_damp_mu)
        else:
          from scipy.special import erf
          mu = self.mm_damp_mu
          for a in range(len(qm_xyz)):
            d = self._min_image(qm_xyz[a] - xyz_s, box)   # r_A - r_s
            r = np.linalg.norm(d, axis=1)
            _, force_kernel = smeared_coulomb(r, mu)      # 1/r^3, or the smeared analogue
            coeff = pchg[a] * q_s * force_kernel          # q_A Q_s * kernel
            f_qm[a] = np.sum(coeff[:, None] * d, axis=0)
            f_site -= coeff[:, None] * d                  # Newton's third law
        # Scatter each site force onto the real MM atoms it is built from
        # (identity for a real atom, chain rule 0.5/0.5 for a midpoint).
        row_of = {int(m): j for j, m in enumerate(mm_idx)}
        f_mm = np.zeros((len(mm_idx), 3))
        for s, contribs in enumerate(scatter):
            for atom, w in contribs:
                f_mm[row_of[atom]] += w * f_site[s]
        return f_qm, f_mm, mm_idx

    def electrostatic_potential(self):

       if self.espf_full:
           # Full-field embedding: phi from all MM charges (no exclusions),
           # QM-QM periodic correction unused (non-periodic).
           nqm_c = len(self.qm_atoms) + len(self.link_atoms)
           return self._full_field_potmm(), np.zeros((nqm_c, nqm_c))

       syspbc=self.mm_systems["sys0"]
       simpbc=self.mm_systems["sim0"]

    #Focus on non-bonded interactions
       forces = { force.__class__.__name__ : force for force in syspbc.getForces() }
       nonbonded = forces['NonbondedForce']
       if nonbonded is None: ValueError(f"Non-bonded interactions are not present, what shall I do?")

    #######################################################################
    #  1. Compute MM potential (needs QM-QM contributions to be removed)  #
    #######################################################################
       potmm=np.zeros((len(self.qm_atoms)))
       for i in range(len(self.qm_atoms)):
           charge, sigma, epsilon = nonbonded.getParticleParameters(self.qm_atoms[i])
           nonbonded.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, sigma, epsilon)
       nonbonded.updateParametersInContext(simpbc.context)
       state = simpbc.context.getState(getEnergy=True)
       e_pbc_no_qm_charge = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879

       for i in range(len(self.qm_atoms)):
           charge, sigma, epsilon = nonbonded.getParticleParameters(self.qm_atoms[i])
           nonbonded.setParticleParameters(self.qm_atoms[i], 1.0*unit.elementary_charge, sigma, epsilon)
           nonbonded.updateParametersInContext(simpbc.context)
           state = simpbc.context.getState(getEnergy=True)
           e_pbc_qm_charge=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879
           potmm[i]=e_pbc_qm_charge-e_pbc_no_qm_charge
           nonbonded.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, sigma, epsilon)

       # Non-periodic embedding has no Ewald QM-QM self-interaction, so the
       # QM-QM correction potential is identically zero. Return a zero matrix
       # (not None) so the Fortran add_potqm_contributions has a valid record.
       if not self._is_periodic():
           return self._pad_potential_for_link_atoms(
               potmm, np.zeros((len(self.qm_atoms), len(self.qm_atoms)))
           )

    #######################################################################
    #                 2. Compute QM pair potential                        #
    #######################################################################
       potqm = np.zeros((len(self.qm_atoms), len(self.qm_atoms)))

    # Create an Ewald system
       sysew=self.mm_systems["sysew"]
       simew=self.mm_systems["simew"]
       forces = { force.__class__.__name__ : force for force in sysew.getForces() }
       nonbondedew = forces['NonbondedForce']

    # Create a non-periodic system
       sysor=self.mm_systems["sysor"]
       simor=self.mm_systems["simor"]
       forcesor = { force.__class__.__name__ : force for force in sysor.getForces() }
       nonbondedor = forcesor['NonbondedForce']

    # Create a system with no charges, both QM and MM (Ewald & Original)
       for i in range(sysew.getNumParticles()):
           nonbondedew.setParticleParameters(i, 0.0*unit.elementary_charge, 0.0, 0.0)
           nonbondedor.setParticleParameters(i, 0.0*unit.elementary_charge, 0.0, 0.0)
       nonbondedew.updateParametersInContext(simew.context)
       nonbondedor.updateParametersInContext(simor.context)

       state = simew.context.getState(getEnergy=True)
       e_ew_no_qm_charge = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879

       state = simor.context.getState(getEnergy=True)
       e_or_no_qm_charge = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879

    #######################################################################
    # 2.1. Compute the QM-QM diagonal potential and remove QM-QM from MM  #
    #         Note 1: the 1/2 factor needs to be corrected later          #
    #         Note 2: here the PME/Ew method                          #
    #######################################################################
       for i in range(len(self.qm_atoms)):
           nonbondedew.setParticleParameters(self.qm_atoms[i], 1.0*unit.elementary_charge, 0.0, 0.0)
           nonbondedew.updateParametersInContext(simew.context)
           state = simew.context.getState(getEnergy=True)
           e_ew_qm_chargei = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879
           potqm[i,i] = e_ew_qm_chargei - e_ew_no_qm_charge
           potmm[i] -= potqm[i,i] #Remove QM-QM interactions from MM potential
           potqm[i,i] *= 2.0
           nonbondedew.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, 0.0, 0.0)

    #######################################################################
    #        2.2. Compute the QM-QM off-diagonal potential                #
    #######################################################################
       for i in range(len(self.qm_atoms)):

           nonbondedew.setParticleParameters(self.qm_atoms[i], 1.0*unit.elementary_charge, 0.0, 0.0)
           nonbondedor.setParticleParameters(self.qm_atoms[i], 1.0*unit.elementary_charge, 0.0, 0.0)

           for j in range(i+1,len(self.qm_atoms)):

               nonbondedew.setParticleParameters(self.qm_atoms[j], 1.0*unit.elementary_charge, 0.0, 0.0)
               nonbondedew.updateParametersInContext(simew.context)
               state = simew.context.getState(getEnergy=True)
               e_ew_qm_chargeij = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879

               nonbondedor.setParticleParameters(self.qm_atoms[j], 1.0*unit.elementary_charge, 0.0, 0.0)
               nonbondedor.updateParametersInContext(simor.context)
               state = simor.context.getState(getEnergy=True)
               e_or_qm_chargeij = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879
               ecorr = e_or_qm_chargeij - e_or_no_qm_charge

               potij = e_ew_qm_chargeij - e_ew_no_qm_charge - ecorr - 0.5*(potqm[j,j] + potqm[i,i])
               potqm[i,j] = potqm[j,i] = potij

               nonbondedew.setParticleParameters(self.qm_atoms[j], 0.0*unit.elementary_charge, 0.0, 0.0)
               nonbondedor.setParticleParameters(self.qm_atoms[j], 0.0*unit.elementary_charge, 0.0, 0.0)

           nonbondedew.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, 0.0, 0.0)
           nonbondedor.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, 0.0, 0.0)

       return self._pad_potential_for_link_atoms(potmm, potqm)
