"""Self-contained Cartesian-GTO evaluator for OQP basis sets.

Builds AO values on arbitrary points from the basis exposed by
``mol.data.get_basis()``, using OQP's normalization/ordering convention
(verified by reproducing OQP's overlap matrix ``OQP::SM`` to ~1e-6 with an
analytic McMurchie-Davidson overlap, and on a grid to <=1e-3).

Convention (from the OQP molden writer):
  * shells store OQP-internal contraction coefficients; the molden/standard
    normalized-primitive coefficient is ``d = coef * NORMS[L] * (2a)^-(L/2+3/4)``.
  * Cartesian component order: L=0 [s]; L=1 [x,y,z]; L=2 [xx,yy,zz,xy,xz,yz];
    L=3 [xxx,yyy,zzz,xyy,xxy,xxz,xzz,yzz,yyz,xyz].
"""
import numpy as np

# sqrt(pi^{3/2} * f[L]); f from the OQP molden writer (per-L primitive factor).
_NORMS = np.sqrt(np.pi * np.sqrt(np.pi) *
                 np.array([1.0, 0.5, 0.75, 1.875, 6.5625, 29.53125, 162.421875]))

# Cartesian components per L in OQP order (lx, ly, lz).
_CART = {
    0: [(0, 0, 0)],
    1: [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    2: [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)],
    3: [(3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0),
        (2, 0, 1), (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)],
}


def _dfac(n):
    """(2n-1)!! with the convention dfac(0)=1."""
    r = 1.0
    k = 2 * n - 1
    while k > 1:
        r *= k
        k -= 2
    return r


class AOBasis:
    def __init__(self, mol):
        basis = mol.data.get_basis()
        self.coords = np.asarray(mol.get_system(), dtype=float).reshape(-1, 3)  # Bohr
        centers = np.asarray(basis["centers"]).astype(int)
        angs = np.asarray(basis["angs"]).astype(int)
        ncontr = np.asarray(basis["ncontr"]).astype(int)
        alpha = np.asarray(basis["alpha"], dtype=float)
        coef = np.asarray(basis["coef"], dtype=float)
        self.nbf = int(basis["nbf"])

        self.shells = []          # (atom_xyz, L, exps, dcoef)
        self.ao_index = []        # (shell_id, (lx,ly,lz), scomp)
        self.ao_atom = []         # atom index per AO
        p0 = 0
        for sh in range(int(basis["nsh"])):
            L = int(angs[sh]); nc = int(ncontr[sh]); at = int(centers[sh])
            a = alpha[p0:p0 + nc]; c = coef[p0:p0 + nc]; p0 += nc
            d = c * _NORMS[L] * (2.0 * a) ** (-(L / 2.0 + 0.75))   # molden-norm coef
            sh_id = len(self.shells)
            self.shells.append((self.coords[at], L, a, d))
            dfL = _dfac(L)
            for (lx, ly, lz) in _CART[L]:
                scomp = np.sqrt(dfL / (_dfac(lx) * _dfac(ly) * _dfac(lz)))
                self.ao_index.append((sh_id, (lx, ly, lz), scomp))
                self.ao_atom.append(at)
        self.ao_atom = np.asarray(self.ao_atom, dtype=int)
        assert len(self.ao_index) == self.nbf, (len(self.ao_index), self.nbf)

    def eval_orbitals(self, mo_coeff_ao, points):
        """Evaluate orbitals (columns of ``mo_coeff_ao``, AO basis) on points."""
        return self.eval_ao(points) @ np.asarray(mo_coeff_ao)

    @staticmethod
    def _nrad(a, L):
        return (2.0 * a / np.pi) ** 0.75 * (4.0 * a) ** (L / 2.0) / np.sqrt(_dfac(L))

    def eval_ao(self, points):
        """Evaluate all AOs on ``points`` (npts, 3) in Bohr -> (npts, nbf)."""
        pts = np.asarray(points, dtype=float).reshape(-1, 3)
        npts = pts.shape[0]
        # radial part per shell (npts,)
        shell_rad = []
        shell_dx = []
        for (A, L, a, d) in self.shells:
            dx = pts - A[None, :]
            r2 = np.einsum("pi,pi->p", dx, dx)
            rad = np.zeros(npts)
            for ai, di in zip(a, d):
                rad += di * self._nrad(ai, L) * np.exp(-ai * r2)
            shell_rad.append(rad)
            shell_dx.append(dx)
        out = np.empty((npts, self.nbf))
        for mu, (sh_id, (lx, ly, lz), scomp) in enumerate(self.ao_index):
            dx = shell_dx[sh_id]
            ang = (dx[:, 0] ** lx) * (dx[:, 1] ** ly) * (dx[:, 2] ** lz)
            out[:, mu] = scomp * ang * shell_rad[sh_id]
        return out

    # ---- analytic overlap (McMurchie-Davidson), for evaluator validation ----
    @staticmethod
    def _E(i, j, Qx, a, b):
        """Hermite expansion coefficient E^{ij}_0 for 1D overlap."""
        p = a + b
        mu = a * b / p
        # E[i][j][t]
        E = np.zeros((i + 1, j + 1, i + j + 1))
        E[0, 0, 0] = np.exp(-mu * Qx * Qx)
        PA = -b * Qx / p
        PB = a * Qx / p
        for ii in range(i + 1):
            for jj in range(j + 1):
                if ii == 0 and jj == 0:
                    continue
                if ii > 0:
                    for t in range(ii + jj + 1):
                        val = 0.0
                        if t - 1 >= 0:
                            val += 1.0 / (2 * p) * E[ii - 1, jj, t - 1]
                        val += PA * E[ii - 1, jj, t]
                        if t + 1 <= (ii - 1) + jj + 1:
                            val += (t + 1) * E[ii - 1, jj, t + 1]
                        E[ii, jj, t] = val
                else:
                    for t in range(ii + jj + 1):
                        val = 0.0
                        if t - 1 >= 0:
                            val += 1.0 / (2 * p) * E[ii, jj - 1, t - 1]
                        val += PB * E[ii, jj - 1, t]
                        if t + 1 <= ii + (jj - 1) + 1:
                            val += (t + 1) * E[ii, jj - 1, t + 1]
                        E[ii, jj, t] = val
        return E[i, j, 0]

    def overlap_analytic(self):
        """Analytic AO overlap matrix in OQP ordering/normalization."""
        nbf = self.nbf
        S = np.zeros((nbf, nbf))
        idx = self.ao_index
        for mu in range(nbf):
            shA, (lax, lay, laz), scA = idx[mu]
            A, LA, aA, dA = self.shells[shA]
            for nu in range(mu, nbf):
                shB, (lbx, lby, lbz), scB = idx[nu]
                B, LB, aB, dB = self.shells[shB]
                Q = A - B
                s = 0.0
                for ai, di in zip(aA, dA):
                    nA = self._nrad(ai, LA)
                    for bi, ei in zip(aB, dB):
                        nB = self._nrad(bi, LB)
                        p = ai + bi
                        Ex = self._E(lax, lbx, Q[0], ai, bi)
                        Ey = self._E(lay, lby, Q[1], ai, bi)
                        Ez = self._E(laz, lbz, Q[2], ai, bi)
                        s += di * ei * nA * nB * (np.pi / p) ** 1.5 * Ex * Ey * Ez
                S[mu, nu] = S[nu, mu] = scA * scB * s
        return S


def make_box_grid(coords, padding=5.0, spacing=0.15):
    """Regular Cartesian grid around ``coords`` (Bohr). Returns
    (origin, (nx,ny,nz), (dx,dy,dz), points[N,3]) in Bohr."""
    coords = np.asarray(coords, dtype=float).reshape(-1, 3)
    lo = coords.min(axis=0) - padding
    hi = coords.max(axis=0) + padding
    n = np.maximum(np.ceil((hi - lo) / spacing).astype(int) + 1, 2)
    axes = [lo[d] + spacing * np.arange(n[d]) for d in range(3)]
    gx, gy, gz = np.meshgrid(axes[0], axes[1], axes[2], indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    return lo, tuple(int(x) for x in n), (spacing, spacing, spacing), pts
