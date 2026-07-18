"""Step-determination engine for the oqp OpenQP optimizer.

Backend-agnostic (NumPy only).  The :class:`OQPEngine` owns the optimization
loop and drives an energy/gradient callback supplied by the caller, so it can be
exercised on analytic test potentials without the compiled OQP core.

Algorithm (v1):

* working coordinates from :mod:`oqp.library.oqp_coords` (redundant internals
  with Cartesian fall-back);
* restricted-step Rational Function Optimization (RFO) for minima and
  partitioned-RFO (P-RFO, eigenvector following) for transition states;
* Schlegel-type diagonal model Hessian, or an optional supplied Cartesian
  Hessian transformed into working coordinates, updated with BFGS (minima) or
  the Bofill SR1/PSB mixture (TS);
* a predictive trust region in Cartesian RMSD, grown/shrunk on the ratio of
  actual to predicted energy change.

References: Banerjee, Adams, Simons, Shepard, J. Phys. Chem. 89, 52 (1985);
Bofill, J. Comput. Chem. 15, 1 (1994); Peng, Ayala, Schlegel, Frisch,
J. Comput. Chem. 17, 49 (1996).
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as sla
import re

from oqp.library.oqp_coords import build_coordinates, CartesianCoordinates

# All dense linear algebra in this module is routed through LAPACK
# (scipy.linalg / numpy.linalg) and BLAS-3 GEMM (the ``@`` operators on the
# B-matrix and Hessian), so it inherits the threaded, vendor-tuned backends
# OpenQP already links against.  The matrices here are O(3N), so the wall-clock
# is dominated by the electronic-structure energy/gradient, not these solves.


def _sym_eigh(a):
    """Symmetric eigensolve via LAPACK (divide-and-conquer 'evd' driver).

    Avoids scipy's default MRRR driver ('evr'): some OpenBLAS builds (e.g. the
    0.3.27 wheel shipped with scipy 1.13) raise ``LinAlgError("Internal Error.")``
    in the threaded dsyevr/MRRR path on large matrices with clustered eigenvalues
    (the diagonal model Hessian used here is exactly that case). The matrix is
    finite and well-conditioned -- numpy.linalg.eigh and the 'evd' driver solve
    it fine -- so route through 'evd', falling back to NumPy if it is unavailable.
    """
    a = 0.5 * (a + a.T)
    try:
        return sla.eigh(a, driver="evd")
    except (sla.LinAlgError, TypeError, ValueError):
        return np.linalg.eigh(a)


class ConvergenceSignal(Exception):
    """Raised by the convergence callback to stop the loop cleanly."""


def parse_frozen_distance_spec(value):
    """Return one-based atom pairs from ``distance(i,j);...`` syntax."""
    text = str(value or "").strip()
    if not text:
        return []
    pattern = re.compile(r"(?:distance|r)\(\s*(\d+)\s*,\s*(\d+)\s*\)", re.I)
    pairs = []
    for item in (part.strip() for part in text.split(";") if part.strip()):
        match = pattern.fullmatch(item)
        if not match:
            raise ValueError(
                "freeze must contain semicolon-separated distance(i,j) pairs"
            )
        pair = (int(match.group(1)), int(match.group(2)))
        if pair[0] == pair[1]:
            raise ValueError("frozen-distance atom indices must be distinct")
        pairs.append(pair)
    return pairs


class OQPEngine:
    def __init__(self, atoms, x0, mode="min", trust=0.2, trust_min=5.0e-3,
                 trust_max=0.5, maxiter=100, follow_mode=0, coordsys="auto",
                 logger=None, initial_hessian=None, initial_gradient=None,
                 masses=None, project_global_rigid_modes=False,
                 frozen_distances=None):
        """Create an optimization engine.

        Parameters
        ----------
        initial_hessian : array_like, optional
            Cartesian Hessian at ``x0`` with shape ``(3N, 3N)``.  When it is
            provided, it is transformed into the active working coordinates
            with the same generalized Wilson-B inverse and rank cutoff used by
            gradients and back-transformation.  Redundant-coordinate
            null directions retain their model-Hessian curvature so they do not
            introduce artificial zero modes.  If omitted, the existing model
            Hessian is used unchanged.
        initial_gradient : array_like, optional
            Cartesian gradient at ``x0``.  When supplied with a Cartesian
            Hessian, the internal-coordinate transformation includes the
            gradient/coordinate-curvature correction and is exact away from a
            stationary point (up to the finite-difference B-matrix derivative).
            Global translation and rotation directions are removed from a real
            Hessian and assigned positive model curvature so numerical zero-mode
            noise cannot become the P-RFO mode followed by a TS search.
        masses : array_like, optional
            Per-atom masses used to define the global translation/rotation
            subspace for a real Hessian. Runtime drivers pass the molecule's
            actual isotope/custom masses; unit masses are the standalone API
            default when no real Hessian projection is needed.
        project_global_rigid_modes : bool, optional
            Replace global translation/rotation curvature by positive model
            curvature for an isolated, invariant molecular Hamiltonian. Runtime
            OpenQP enables this for non-QM/MM jobs. The standalone default is
            false because an external field or fixed embedding can give genuine
            lab-frame rigid-motion curvature.
        frozen_distances : sequence of pair, optional
            One-based atom-index pairs whose initial distances are held fixed.
            The optimizer projects gradients into the linearized constraint
            tangent space and applies a SHAKE-like Cartesian correction after
            every trial step.
        """
        self.atoms = np.asarray(atoms, dtype=int).reshape(-1)
        self.x = np.asarray(x0, dtype=float).reshape(-1)
        if masses is None:
            self.masses = np.ones(self.atoms.size, dtype=float)
        else:
            self.masses = np.asarray(masses, dtype=float).reshape(-1)
        if (self.masses.shape != self.atoms.shape
                or not np.all(np.isfinite(self.masses))
                or np.any(self.masses <= 0.0)):
            raise ValueError("masses must be a finite positive value for every atom")
        self.sqrt_cartesian_mass = np.sqrt(np.repeat(self.masses, 3))
        self.project_global_rigid_modes = bool(project_global_rigid_modes)
        self.frozen_distances = self._validate_frozen_distances(
            frozen_distances or []
        )
        xyz0 = self.x.reshape(-1, 3)
        self._frozen_distance_targets = [
            float(np.linalg.norm(xyz0[first] - xyz0[second]))
            for first, second in self.frozen_distances
        ]
        self.mode = mode
        self.trust = float(trust)
        self.trust_min = float(trust_min)
        self.trust_max = float(trust_max)
        self.maxiter = int(maxiter)
        self.follow_mode = int(follow_mode)
        self.logger = logger

        self.coords = build_coordinates(self.atoms, self.x, coordsys=coordsys)
        # Report the requested coordinate system by name (TRIC and RIC are both
        # RedundantInternalCoordinates objects, so the class name is ambiguous),
        # and flag when it fell back to Cartesians (e.g. for a linear molecule).
        requested = (coordsys or "tric").lower()
        if requested in ("auto",):
            requested = "tric"
        label = {"cart": "CART", "cartesian": "CART"}.get(requested, requested.upper())
        if (type(self.coords).__name__ == "CartesianCoordinates"
                and requested not in ("cart", "cartesian")):
            label = f"{label}->CART(fallback)"
        elif (requested == "tric"
              and type(self.coords).__name__ == "RedundantInternalCoordinates"):
            _, coordinate_rank = self.coords._g_inverse(
                self.coords.b_matrix(self.x)
            )
            if coordinate_rank < self.x.size:
                label = "TRIC->RIC(fallback)"
        elif (requested == "dlc"
              and type(self.coords).__name__ == "RedundantInternalCoordinates"):
            label = "DLC->RIC(fallback)"
        self.coordsys = label
        model_hessian = self.coords.guess_hessian(self.x)
        if initial_hessian is None:
            self.H = model_hessian
        else:
            self.H = self._transform_initial_hessian(initial_hessian,
                                                     model_hessian,
                                                     initial_gradient)
        self._prev = None

    def _validate_frozen_distances(self, pairs):
        validated = []
        natom = self.atoms.size
        for pair in pairs:
            if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                raise ValueError("each frozen distance must contain two atom indices")
            try:
                # Public input uses the chemistry convention: atoms start at 1.
                first, second = int(pair[0]) - 1, int(pair[1]) - 1
            except (TypeError, ValueError) as exc:
                raise ValueError("frozen-distance atom indices must be integers") from exc
            if first == second or min(first, second) < 0 or max(first, second) >= natom:
                raise ValueError(
                    "frozen-distance atom indices must be distinct and between 1 and %d"
                    % natom
                )
            canonical = tuple(sorted((first, second)))
            if canonical not in validated:
                validated.append(canonical)
        return validated

    def _constraint_jacobian(self, x):
        """Cartesian Jacobian rows for the active frozen distances."""
        if not self.frozen_distances:
            return np.zeros((0, self.x.size), dtype=float)
        xyz = np.asarray(x, dtype=float).reshape(-1, 3)
        rows = []
        for first, second in self.frozen_distances:
            delta = xyz[first] - xyz[second]
            norm = float(np.linalg.norm(delta))
            if norm <= 1.0e-14:
                raise ValueError("cannot freeze a zero atom-pair distance")
            unit = delta / norm
            row = np.zeros(self.x.size, dtype=float)
            row[3 * first:3 * first + 3] = unit
            row[3 * second:3 * second + 3] = -unit
            rows.append(row)
        return np.asarray(rows)

    def _project_constraint_tangent(self, vector, x=None):
        """Remove Cartesian components normal to all distance constraints."""
        vector = np.asarray(vector, dtype=float).reshape(-1)
        jac = self._constraint_jacobian(self.x if x is None else x)
        if not jac.size:
            return vector
        gram = jac @ jac.T
        coeff = np.linalg.pinv(gram, rcond=1.0e-12) @ (jac @ vector)
        return vector - jac.T @ coeff

    def _enforce_frozen_distance_values(self, x):
        """Return coordinates corrected to the stored pair distances."""
        if not self.frozen_distances:
            return np.asarray(x, dtype=float).reshape(-1)
        xyz = np.asarray(x, dtype=float).reshape(-1, 3).copy()
        for _ in range(50):
            largest = 0.0
            for (first, second), target in zip(
                    self.frozen_distances, self._frozen_distance_targets):
                delta = xyz[first] - xyz[second]
                distance = float(np.linalg.norm(delta))
                if distance <= 1.0e-14:
                    raise ValueError("cannot enforce a zero atom-pair distance")
                error = distance - target
                largest = max(largest, abs(error))
                correction = 0.5 * error * (delta / distance)
                xyz[first] -= correction
                xyz[second] += correction
            if largest <= 1.0e-12:
                break
        return xyz.reshape(-1)

    def _transform_initial_hessian(self, initial_hessian, model_hessian,
                                   initial_gradient=None):
        """Validate and transform a Cartesian Hessian to working coordinates."""
        ncart = self.x.size
        try:
            h_cart = np.asarray(initial_hessian, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError("initial_hessian must be a numeric Cartesian matrix") from exc

        expected = (ncart, ncart)
        if h_cart.shape != expected:
            raise ValueError(
                f"initial_hessian must have shape {expected}, got {h_cart.shape}"
            )
        if not np.all(np.isfinite(h_cart)):
            raise ValueError("initial_hessian must contain only finite values")
        if not np.allclose(h_cart, h_cart.T, rtol=1.0e-8, atol=1.0e-10):
            raise ValueError("initial_hessian must be symmetric")

        # Remove harmless roundoff asymmetry before transforming.  For dq = B dx,
        # the minimum-norm Cartesian displacement is dx = B+ dq; therefore
        # Hq = (B+)T Hx B+.  In a redundant coordinate system B B+ is the
        # projector onto the physical internal-coordinate subspace.  Complete
        # the orthogonal null space with the established model Hessian rather
        # than leaving spurious zero-curvature modes in the RFO problem.
        h_cart = 0.5 * (h_cart + h_cart.T)
        bmat = np.asarray(self.coords.b_matrix(self.x), dtype=float)
        if isinstance(self.coords, CartesianCoordinates):
            b_pinv = bmat.T
            active_rank = ncart
        else:
            # Use exactly the optimizer's active internal-coordinate subspace.
            # np.linalg.pinv has a much looser default cutoff and can amplify a
            # near-null RIC/DLC direction that grad_to_q/back_transform discard,
            # which can reorder the P-RFO modes of a real TS Hessian.
            g_inverse, rank = self.coords._g_inverse(bmat)
            if rank == 0:
                raise ValueError(
                    "initial_hessian cannot be transformed through a rank-zero B matrix"
                )
            b_pinv = bmat.T @ g_inverse
            active_rank = rank

        # Hx = B.T Hq B + sum_i gq_i d2(q_i)/dx2.  The second term vanishes at
        # a stationary point, but a TS guess can retain a non-negligible
        # gradient.  Numerically differentiate the already available B matrix
        # so injected real Hessians preserve curvature/mode ordering there too.
        coordinate_curvature = np.zeros_like(h_cart)
        if initial_gradient is not None:
            try:
                g_cart = np.asarray(initial_gradient, dtype=float).reshape(-1)
            except (TypeError, ValueError) as exc:
                raise ValueError("initial_gradient must be a numeric Cartesian vector") from exc
            if g_cart.shape != (ncart,):
                raise ValueError(
                    "initial_gradient must have shape (%d,), got %s"
                    % (ncart, g_cart.shape)
                )
            if not np.all(np.isfinite(g_cart)):
                raise ValueError("initial_gradient must contain only finite values")
            if not isinstance(self.coords, CartesianCoordinates) \
                    and np.linalg.norm(g_cart) > 0.0:
                g_work = np.asarray(
                    self.coords.grad_to_q(self.x, g_cart), dtype=float
                ).reshape(-1)
                fd_step = 1.0e-4
                for column in range(ncart):
                    xp = self.x.copy()
                    xm = self.x.copy()
                    xp[column] += fd_step
                    xm[column] -= fd_step
                    dbdx = (
                        np.asarray(self.coords.b_matrix(xp), dtype=float)
                        - np.asarray(self.coords.b_matrix(xm), dtype=float)
                    ) / (2.0 * fd_step)
                    coordinate_curvature[:, column] = g_work @ dbdx
                coordinate_curvature = 0.5 * (
                    coordinate_curvature + coordinate_curvature.T
                )

        h_work = b_pinv.T @ (h_cart - coordinate_curvature) @ b_pinv

        projector = bmat @ b_pinv
        null_projector = np.eye(projector.shape[0]) - projector
        h_work += null_projector.T @ model_hessian @ null_projector

        # A molecular Cartesian Hessian has five or six exact global zero modes.
        # Numerical Hessians and finite integral thresholds commonly make those
        # modes slightly negative. TRIC/Cartesian working coordinates retain the
        # global motions, so those artifacts could sort below the physical TS
        # vibration and make follow=0 translate/rotate the molecule. Map the
        # physical rigid-motion displacement subspace through B, remove all
        # Hessian coupling to it, and restore its established positive model
        # curvature. Relative fragment motion remains untouched because only
        # whole-system translation/rotation vectors are included.
        # RIC/DLC already remove global motions (active_rank < 3N). Their B+
        # gauge can contain an arbitrary COM component, so applying a second
        # mass-metric projector would misclassify and erase physical internal
        # curvature. Only full-rank TRIC/Cartesian coordinates carry the global
        # modes that require regularization.
        external_work = np.empty((h_work.shape[0], 0))
        if self.project_global_rigid_modes and active_rank == ncart:
            external_mw = self._global_external_mass_weighted_basis()
            # b_pinv maps a working-coordinate displacement to Cartesian dx.
            # Left-multiplying by sqrt(M), then projecting onto the mass-weighted
            # external basis, identifies exactly those working directions whose
            # realized Cartesian step contains global rigid motion.
            mass_weighted_back_transform = self.sqrt_cartesian_mass[:, None] * b_pinv
            external_work = mass_weighted_back_transform.T @ external_mw
        if external_work.size:
            left, singular, _ = np.linalg.svd(external_work, full_matrices=False)
            scale = float(singular[0]) if singular.size else 0.0
            keep = singular > 1.0e-8 * max(scale, 1.0)
            if np.any(keep):
                external_vectors = left[:, keep]
                external_projector = external_vectors @ external_vectors.T
                internal_projector = np.eye(h_work.shape[0]) - external_projector
                h_work = (
                    internal_projector @ h_work @ internal_projector
                    + external_projector @ model_hessian @ external_projector
                )
        h_work = 0.5 * (h_work + h_work.T)
        if not np.all(np.isfinite(h_work)):
            raise ValueError("initial_hessian produced a non-finite working Hessian")
        return h_work

    def _global_external_mass_weighted_basis(self):
        """Mass-weighted basis for whole-system rigid motions."""
        xyz = self.x.reshape(-1, 3)
        natom = xyz.shape[0]
        external = np.zeros((3 * natom, 6), dtype=float)
        centered = xyz - np.average(xyz, axis=0, weights=self.masses)
        sqrt_mass = np.sqrt(self.masses)
        for atom, position in enumerate(centered):
            external[3 * atom:3 * atom + 3, :3] = np.eye(3) * sqrt_mass[atom]
            for axis in range(3):
                unit = np.zeros(3)
                unit[axis] = 1.0
                external[3 * atom:3 * atom + 3, 3 + axis] = np.cross(
                    unit, position
                ) * sqrt_mass[atom]
        left, singular, _ = np.linalg.svd(external, full_matrices=False)
        scale = float(singular[0]) if singular.size else 0.0
        rank = int(np.count_nonzero(singular > 1.0e-10 * max(scale, 1.0)))
        return left[:, :rank]

    # -- public driver ------------------------------------------------------ #
    def run(self, energy_gradient, on_converged=None):
        """Optimize until ``on_converged`` raises or ``maxiter`` is reached.

        ``energy_gradient(x_flat) -> (E, g_flat)`` returns energy (Hartree) and
        Cartesian gradient (Hartree/Bohr).  ``on_converged()`` is called after
        every evaluation and may raise to stop (mirrors OQP's ``StopIteration``
        convergence protocol).
        """
        try:
            for _ in range(self.maxiter):
                e, g_cart = energy_gradient(self.x)
                if on_converged is not None:
                    on_converged()
                self.x = self._take_step(e, np.asarray(g_cart, dtype=float).reshape(-1))
        except (StopIteration, ConvergenceSignal):
            pass
        return self.x

    def rebase_previous_objective(self, energy, gradient):
        """Re-evaluate the stored previous point for a changed objective.

        Adaptive objectives may change a scalar parameter between two geometry
        evaluations.  Their caller can recombine the previous geometry's
        parameter-independent energy and Cartesian gradient at the new
        parameter, then use this helper to keep the BFGS secant and trust-ratio
        comparison on one objective.  The previous geometry is retained with
        the accepted step so the rebased gradient is transformed in the same
        coordinate frame in which it was evaluated.

        Returns ``False`` before the first step, when no previous point exists.
        """

        if self._prev is None:
            return False

        energy = float(energy)
        gradient = np.asarray(gradient, dtype=float).reshape(-1)
        previous_x = np.asarray(self._prev["x"], dtype=float).reshape(-1)
        if not np.isfinite(energy):
            raise ValueError("rebased objective energy must be finite")
        if gradient.shape != previous_x.shape:
            raise ValueError(
                "rebased objective gradient must match the previous geometry"
            )
        if not np.all(np.isfinite(gradient)):
            raise ValueError("rebased objective gradient must be finite")

        gradient = self._project_constraint_tangent(gradient, previous_x)
        gradient_q = self.coords.grad_to_q(previous_x, gradient)
        displacement_q = self._prev["dq"]

        self._prev["e"] = energy
        self._prev["g_q"] = gradient_q
        self._prev["pred"] = float(
            gradient_q @ displacement_q
            + 0.5 * displacement_q @ self.H @ displacement_q
        )
        return True

    # -- one macro-iteration ------------------------------------------------ #
    def _bounded_cartesian_recovery(self, gradient):
        """Reject an unstable internal step and take one bounded Cartesian step.

        A failed DLC/RIC back-transformation must never contaminate the stored
        secant pair or Hessian.  The next macroiteration can build a normal
        internal step again from this finite geometry.
        """
        self.trust = max(self.trust * 0.5, self.trust_min)
        self.H = self.coords.guess_hessian(self.x)
        self._prev = None
        gradient = np.asarray(gradient, dtype=float).reshape(-1)
        if not np.all(np.isfinite(gradient)):
            return self.x.copy()
        # Normalize in two stages to avoid overflowing the norm of a very
        # large, yet finite, electronic gradient.
        scale = float(np.max(np.abs(gradient)))
        if scale <= 1.0e-300 or not np.isfinite(scale):
            return self.x.copy()
        direction = gradient / scale
        rms = float(np.sqrt(np.mean(direction * direction)))
        if rms <= 1.0e-14 or not np.isfinite(rms):
            return self.x.copy()
        return self._enforce_frozen_distance_values(
            self.x - direction * (self.trust / rms)
        )

    def _take_step(self, e, g_cart):
        b = self.coords.b_matrix(self.x)
        g_cart = self._project_constraint_tangent(g_cart, self.x)
        if (not np.isfinite(e) or not np.all(np.isfinite(b))
                or not np.all(np.isfinite(g_cart))
                or not np.all(np.isfinite(self.H))):
            return self._bounded_cartesian_recovery(g_cart)
        g_q = self.coords.grad_to_q(self.x, g_cart)
        if not np.all(np.isfinite(g_q)):
            return self._bounded_cartesian_recovery(g_cart)

        if self._prev is not None:
            s = self._prev["dq"]
            y = g_q - self._prev["g_q"]
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                self.H = self._update_hessian(self.H, s, y)
            if not np.all(np.isfinite(self.H)):
                return self._bounded_cartesian_recovery(g_cart)
            actual = e - self._prev["e"]
            pred = self._prev["pred"]
            self._update_trust(actual, pred, self._prev["cart_step"])

        dq = self._rfo_step(self.H, g_q)
        dq = self._restrict_to_trust(b, dq)
        if not np.all(np.isfinite(dq)):
            return self._bounded_cartesian_recovery(g_cart)

        x_new, ok = self.coords.back_transform(self.x, dq)
        if not ok:
            # Shrink and retry once; otherwise use a finite Cartesian recovery.
            self.trust = max(self.trust * 0.5, self.trust_min)
            dq = self._restrict_to_trust(b, dq)
            x_new, ok = self.coords.back_transform(self.x, dq)
            if not ok:
                return self._bounded_cartesian_recovery(g_cart)

        x_new = self._enforce_frozen_distance_values(x_new)
        if not np.all(np.isfinite(x_new)):
            return self._bounded_cartesian_recovery(g_cart)

        dq_taken = self.coords.q_displacement(self.coords.q(x_new),
                                              self.coords.q(self.x))
        cart_step = self.coords.cart_rmsd(self.x, x_new)
        if not np.all(np.isfinite(dq_taken)) or not np.isfinite(cart_step):
            return self._bounded_cartesian_recovery(g_cart)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            pred = float(g_q @ dq_taken + 0.5 * dq_taken @ self.H @ dq_taken)
        if not np.isfinite(pred):
            return self._bounded_cartesian_recovery(g_cart)
        self._prev = {"e": e, "g_q": g_q, "dq": dq_taken,
                      "pred": pred, "cart_step": cart_step,
                      "x": self.x.copy()}
        return x_new

    # -- RFO / P-RFO -------------------------------------------------------- #
    def _rfo_step(self, h, g):
        if self.mode == "ts":
            return self._prfo_step(h, g)
        return self._simple_rfo_step(h, g)

    @staticmethod
    def _simple_rfo_step(h, g):
        """Minimization RFO: lowest-eigenvalue augmented-Hessian step."""
        n = len(g)
        aug = np.zeros((n + 1, n + 1))
        aug[:n, :n] = h
        aug[:n, n] = g
        aug[n, :n] = g
        w, v = _sym_eigh(aug)
        vec = v[:, 0]
        if abs(vec[-1]) < 1.0e-10:
            # Degenerate scaling: fall back to a damped steepest-descent step.
            return -g / (np.linalg.norm(g) + 1.0)
        return vec[:n] / vec[n]

    def _prfo_step(self, h, g):
        """Partitioned RFO for TS: maximize along one mode, minimize the rest."""
        w, v = _sym_eigh(h)
        f = v.T @ g  # gradient in the Hessian eigenbasis
        t = self._select_follow_mode(w, v)

        step = np.zeros_like(g)
        # Maximized mode: take the upper root of the 2x2 augmented problem.
        lam = w[t]
        nu_max = 0.5 * (lam + np.sqrt(lam * lam + 4.0 * f[t] * f[t]))
        denom = lam - nu_max
        if abs(denom) > 1.0e-10:
            step += (-f[t] / denom) * v[:, t]

        # Minimized subspace: RFO over the remaining modes.
        idx = [i for i in range(len(w)) if i != t]
        if idx:
            lam_min = w[idx]
            f_min = f[idx]
            m = len(idx)
            aug = np.zeros((m + 1, m + 1))
            aug[:m, :m] = np.diag(lam_min)
            aug[:m, m] = f_min
            aug[m, :m] = f_min
            wa, va = _sym_eigh(aug)
            nu_min = wa[0]
            for jj, i in enumerate(idx):
                d = lam_min[jj] - nu_min
                if abs(d) > 1.0e-10:
                    step += (-f_min[jj] / d) * v[:, i]
        self._followed = v[:, t].copy()
        return step

    def _select_follow_mode(self, w, v):
        """Pick the mode to maximize: lowest curvature, with mode following."""
        prev = getattr(self, "_followed", None)
        if prev is not None:
            overlaps = np.abs(v.T @ prev)
            return int(np.argmax(overlaps))
        order = np.argsort(w)
        if self.follow_mode < 0 or self.follow_mode >= len(order):
            raise ValueError(
                "follow_mode=%d is out of range for %d working-coordinate modes"
                % (self.follow_mode, len(order))
            )
        return int(order[self.follow_mode])

    # -- trust region & step restriction ------------------------------------ #
    def _restrict_to_trust(self, b, dq):
        # Estimate the Cartesian motion of this internal step and scale so its
        # RMSD does not exceed the trust radius.
        if not np.all(np.isfinite(b)) or not np.all(np.isfinite(dq)):
            return np.full_like(dq, np.nan)
        ginv, _ = self.coords._g_inverse(b) if hasattr(self.coords, "_g_inverse") \
            else (np.eye(b.shape[0]), 0)
        if not np.all(np.isfinite(ginv)):
            return np.full_like(dq, np.nan)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            dx = b.T @ (ginv @ dq) if hasattr(self.coords, "_g_inverse") else dq
        if not np.all(np.isfinite(dx)):
            return np.full_like(dq, np.nan)
        rmsd = float(np.sqrt(np.mean(dx ** 2)))
        if rmsd > self.trust and rmsd > 1.0e-14:
            dq = dq * (self.trust / rmsd)
        return dq

    def _update_trust(self, actual, pred, cart_step):
        if abs(pred) < 1.0e-14:
            return
        ratio = actual / pred
        if ratio < 0.25:
            self.trust = max(self.trust * 0.25, self.trust_min)
        elif ratio > 0.75 and cart_step >= 0.8 * self.trust:
            self.trust = min(self.trust * 2.0, self.trust_max)

    # -- Hessian updates ---------------------------------------------------- #
    def _update_hessian(self, h, s, y):
        if self.mode == "ts":
            return self._bofill_update(h, s, y)
        return self._bfgs_update(h, s, y)

    @staticmethod
    def _bfgs_update(h, s, y):
        sy = float(s @ y)
        hs = h @ s
        shs = float(s @ hs)
        if shs <= 1.0e-12:
            return h
        # Powell damping keeps the minimum-search Hessian positive definite
        # when a noisy/flat electronic surface violates the raw BFGS curvature
        # condition.  Previously those updates were simply skipped, leaving a
        # stale model Hessian that could drive repeated oscillatory steps.
        if sy < 0.2 * shs:
            denominator = shs - sy
            if denominator <= 1.0e-14:
                return h
            theta = 0.8 * shs / denominator
            y = theta * y + (1.0 - theta) * hs
            sy = float(s @ y)
        if sy > 1.0e-12:
            h = h + np.outer(y, y) / sy - np.outer(hs, hs) / shs
        return h

    @staticmethod
    def _bofill_update(h, s, y):
        ss = float(s @ s)
        if ss < 1.0e-14:
            return h
        dg = y - h @ s            # gradient-difference residual
        dgs = float(dg @ s)
        dgdg = float(dg @ dg)
        if dgdg < 1.0e-14:
            return h
        phi = (dgs * dgs) / (ss * dgdg)
        # SR1 (Murtagh-Sargent) part.
        if abs(dgs) > 1.0e-12:
            h_sr1 = h + np.outer(dg, dg) / dgs
        else:
            h_sr1 = h
        # PSB part.
        h_psb = (h + (np.outer(dg, s) + np.outer(s, dg)) / ss
                 - dgs * np.outer(s, s) / (ss * ss))
        return phi * h_sr1 + (1.0 - phi) * h_psb
