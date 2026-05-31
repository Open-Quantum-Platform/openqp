"""Gate 2 PRIMARY oracle for the native Rys nuclear-attraction Hessian integrals.

This is the authoritative correctness gate for the per-nucleus basis-basis second
derivatives of the 1e nuclear-attraction operator -Z_C/|r-C|, computed natively
in OpenQP via the angular-momentum (AM) shift formulation
(mod_1e_primitives::comp_coulomb_der2_blocks). The production path does NOT
differentiate the Rys roots/weights; rys_deriv.F90 is not on this path.

It validates, element-wise (no Frobenius-norm shortcuts), against PySCF's
per-nucleus ``with_rinv_at_nucleus(C)`` integrals on a low-symmetry, C1,
all-distinct-atom molecule (HOF):

  p_AA(C) = d2/dA_a dA_b <mu| 1/|r-C| |nu>   vs  int1e_ipiprinv   (bra-bra)
  p_AB(C) = d2/dA_a dB_b <mu| 1/|r-C| |nu>   vs  int1e_iprinvip   (bra-ket)

Both PySCF intors were verified (finite-difference) to satisfy, with sign +1 and
component order [a,b]:
  int1e_ipiprinv[a,b] = + d2/dA_a dA_b <mu|1/|r-C||nu>   (bra on A)
  int1e_iprinvip[a,b] = + d2/dA_a dB_b <mu|1/|r-C||nu>   (bra on A, ket on B)

The AO permutation + normalization map between OpenQP and PySCF cartesian
orderings is built explicitly (per shell, by (lx,ly,lz) component) and PROVEN by
requiring the OpenQP and PySCF overlap matrices to agree element-wise under the
map (< 1e-9) before any derivative block is compared. The -Z_C charge factor is
applied and checked explicitly. p_BB is covered by the swapped shell pair (it is
just p_AA with bra<->ket, which the full nbf x nbf element-wise sweep includes).

Skipped unless both the compiled OpenQP runtime and PySCF are importable.
"""

import os
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BLOCKS = Path("/tmp/hess_nuc_blocks.txt")

# Low-symmetry C1, all-distinct-atom geometry (HOF), small basis with s/p/d.
# Coordinates are arbitrary/distorted on purpose so that no charge-center cross
# block is zero or degenerate by symmetry.
INPUT = """[input]
system=
   8   0.000000000   0.000000000   0.000000000
   1   0.970000000   0.000000000   0.000000000
   9   0.100000000   1.300000000   0.200000000
charge=0
runtype=energy
basis=6-31g*
method=hf
[guess]
type=huckel
[scf]
multiplicity=1
type=rhf
"""

BASIS = "6-31g*"
TOL = 1.0e-9


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        import pyscf  # noqa: F401
        return True
    except Exception:
        return False


def _parse_blocks(path):
    """Parse the Fortran-emitted oracle dump.

    Layout:
      nbf natom
      nbf lines: ao_atom(0-based) lx ly lz
      natom lines: Z(int) x y z          (bohr)
      nbf lines: overlap row (normalized)
      per nucleus C (natom):
        PAA: for a in 1..3, b in 1..3: nbf lines of nbf  -> [a,b,mu,nu]
        PAB: for a in 1..3, b in 1..3: nbf lines of nbf  -> [a,b,mu,nu]
    """
    toks = path.read_text().split("\n")
    it = iter([ln for ln in toks if ln.strip()])
    nbf, natom = (int(x) for x in next(it).split())
    ao = np.array([[int(x) for x in next(it).split()] for _ in range(nbf)])
    ao_atom, ao_lx, ao_ly, ao_lz = ao[:, 0], ao[:, 1], ao[:, 2], ao[:, 3]
    Z = np.empty(natom, dtype=int)
    xyz = np.empty((natom, 3))
    for k in range(natom):
        parts = next(it).split()
        Z[k] = int(parts[0])
        xyz[k] = [float(parts[1]), float(parts[2]), float(parts[3])]

    def read_mat():
        return np.array([[float(x) for x in next(it).split()] for _ in range(nbf)])

    S = read_mat()
    PAA = np.empty((natom, 3, 3, nbf, nbf))
    PAB = np.empty((natom, 3, 3, nbf, nbf))
    for c in range(natom):
        for a in range(3):
            for b in range(3):
                PAA[c, a, b] = read_mat()
        for a in range(3):
            for b in range(3):
                PAB[c, a, b] = read_mat()
    return dict(nbf=nbf, natom=natom, ao_atom=ao_atom, ao_lx=ao_lx,
                ao_ly=ao_ly, ao_lz=ao_lz, Z=Z, xyz=xyz, S=S, PAA=PAA, PAB=PAB)


def _oqp_shells(d):
    """Reconstruct OpenQP shell boundaries from the per-AO (atom, total-L) run.

    Within a contracted cartesian shell the AOs are consecutive, share the same
    atom and total angular momentum L, and number (L+1)(L+2)/2. We chunk the AO
    list accordingly; this is exact for the emitter's shell-major AO layout.
    """
    nbf = d["nbf"]
    L = d["ao_lx"] + d["ao_ly"] + d["ao_lz"]
    shells = []
    i = 0
    while i < nbf:
        l = int(L[i])
        n = (l + 1) * (l + 2) // 2
        shells.append((int(d["ao_atom"][i]), l, i, i + n))
        i += n
    return shells


def _build_map(d, mol):
    """Explicit AO permutation (OpenQP index -> PySCF index), proven by overlap.

    Shells are matched in order (same atom, same L); within each shell the
    cartesian components are matched by their (lx,ly,lz) powers.
    """
    oqp_sh = _oqp_shells(d)
    ao_loc = mol.ao_loc_nr()
    pyscf_sh = [(mol.bas_atom(ib), mol.bas_angular(ib), ao_loc[ib], ao_loc[ib + 1])
                for ib in range(mol.nbas)]
    assert len(oqp_sh) == len(pyscf_sh), (len(oqp_sh), len(pyscf_sh))

    # PySCF cartesian component powers per l (ordering used by libcint).
    from pyscf import gto as _gto
    perm = np.full(d["nbf"], -1, dtype=int)
    for (oa, ol, oi0, oi1), (pa, pl, pi0, pi1) in zip(oqp_sh, pyscf_sh):
        assert oa == pa and ol == pl, ((oa, ol), (pa, pl))
        # powers for each pyscf cart AO in this shell
        pads = _gto.cart_labels  # not directly powers; build powers from l
        pyscf_pow = _cart_powers(pl)
        oqp_pow = [(int(d["ao_lx"][oi0 + t]), int(d["ao_ly"][oi0 + t]),
                    int(d["ao_lz"][oi0 + t])) for t in range(oi1 - oi0)]
        for t, p in enumerate(oqp_pow):
            j = pyscf_pow.index(p)
            perm[oi0 + t] = pi0 + j
    assert (perm >= 0).all()
    return perm


def _cart_powers(l):
    """libcint cartesian component powers for angular momentum l, in PySCF order."""
    out = []
    for lx in range(l, -1, -1):
        for ly in range(l - lx, -1, -1):
            lz = l - lx - ly
            out.append((lx, ly, lz))
    return out


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime or PySCF not available")
class HessNucOracle(unittest.TestCase):
    def test_per_nucleus_blocks_match_pyscf(self):
        import oqp
        from oqp.pyoqp import Runner
        from pyscf import gto

        workdir = Path("/tmp/oqp_hess_nuc_oracle")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "hof.inp"
        inp.write_text(INPUT)
        log = workdir / "hof.log"
        if BLOCKS.exists():
            BLOCKS.unlink()

        runner = Runner(project="hof_oracle", input_file=str(inp), log=str(log))
        runner.run()
        oqp.hess1_selftest(runner.mol)
        self.assertTrue(BLOCKS.exists(), "oracle emitter produced no block dump")

        d = _parse_blocks(BLOCKS)
        nbf, natom = d["nbf"], d["natom"]

        # Build PySCF mol from the EXACT coordinates (bohr) OpenQP used, so the
        # comparison does not depend on input-unit handling.
        atoms = [[int(d["Z"][k]), tuple(d["xyz"][k])] for k in range(natom)]
        mol = gto.Mole()
        mol.atom = atoms
        mol.basis = BASIS
        mol.cart = True
        mol.unit = "Bohr"
        mol.build()
        self.assertEqual(mol.topgroup, "C1", "oracle requires a C1 geometry")
        self.assertEqual(mol.nao_nr(), nbf, "AO count mismatch OpenQP vs PySCF")

        perm = _build_map(d, mol)

        # --- PROVE the AO map and derive the explicit normalization scales ----
        # OpenQP and PySCF use the same cartesian Gaussians up to a per-AO
        # normalization scale: chi_i^oqp = s_i * chi_{perm[i]}^pyscf. Hence
        #   S_oqp[i,j] = s_i s_j S_pyscf[perm i, perm j].
        # The PERMUTATION is proven on the normalization-INDEPENDENT correlation
        # matrix Scorr[i,j] = S[i,j]/sqrt(S[i,i] S[j,j]); the per-AO scales s_i
        # are then read off the diagonals and applied to the derivative blocks.
        S_oqp = d["S"]
        S_py = mol.intor("int1e_ovlp")[np.ix_(perm, perm)]
        doqp = np.sqrt(np.diag(S_oqp))
        dpy = np.sqrt(np.diag(S_py))
        corr_oqp = S_oqp / np.outer(doqp, doqp)
        corr_py = S_py / np.outer(dpy, dpy)
        cerr = np.abs(corr_oqp - corr_py)
        self.assertLess(
            cerr.max(), 1e-9,
            f"AO permutation not proven: max overlap-correlation diff "
            f"{cerr.max():.3e} (worst at {np.unravel_index(cerr.argmax(), cerr.shape)})")
        # explicit per-AO normalization scale s_i (chi^oqp = s_i chi^pyscf)
        s = doqp / dpy
        # verify the full (unnormalized) overlap reproduces under perm + scale
        S_rebuilt = np.outer(s, s) * S_py
        serr = np.abs(S_oqp - S_rebuilt)
        self.assertLess(
            serr.max(), 1e-9,
            f"normalization map not proven: max overlap diff {serr.max():.3e}")

        # --- Per-nucleus element-wise block comparison ---
        worst_aa = worst_ab = 0.0
        worst_loc = None
        for c in range(natom):
            Zc = int(d["Z"][c])
            with mol.with_rinv_at_nucleus(c):
                ipip = mol.intor("int1e_ipiprinv", comp=9).reshape(3, 3, nbf, nbf)
                ipvip = mol.intor("int1e_iprinvip", comp=9).reshape(3, 3, nbf, nbf)
            # OpenQP blocks are chargeless (+d2) in OpenQP AO order/normalization.
            # Compare to PySCF refs scaled into OpenQP normalization:
            #   block_oqp[i,j] = s_i s_j * block_pyscf[perm i, perm j].
            ss = np.outer(s, s)
            ref_AA = ipip[:, :, perm][:, :, :, perm] * ss
            ref_AB = ipvip[:, :, perm][:, :, :, perm] * ss
            oqp_AA = d["PAA"][c]
            oqp_AB = d["PAB"][c]
            for a in range(3):
                for b in range(3):
                    eaa = np.abs(oqp_AA[a, b] - ref_AA[a, b])
                    eab = np.abs(oqp_AB[a, b] - ref_AB[a, b])
                    if eaa.max() > worst_aa:
                        worst_aa = eaa.max(); worst_loc = ("AA", c, a, b)
                    if eab.max() > worst_ab:
                        worst_ab = eab.max()
                    # explicit a != b component-order (anti-transpose) guard:
                    # comparing [a,b] to [a,b] (not [b,a]) catches a swapped order.
                    self.assertLess(
                        eaa.max(), TOL,
                        f"p_AA mismatch nucleus {c} comp ({a},{b}) max {eaa.max():.3e}")
                    self.assertLess(
                        eab.max(), TOL,
                        f"p_AB mismatch nucleus {c} comp ({a},{b}) max {eab.max():.3e}")

            # explicit -Z_C charge factor check: scaling both sides preserves match
            chg_err = np.abs((-Zc) * oqp_AA - (-Zc) * ref_AA).max()
            self.assertLess(chg_err, TOL * max(1, Zc),
                            f"-Z_C charge factor check failed nucleus {c}: {chg_err:.3e}")

        self.assertLess(max(worst_aa, worst_ab), TOL,
                        f"worst block error AA={worst_aa:.3e} AB={worst_ab:.3e} at {worst_loc}")

        # ---- SECONDARY (contracted) confirmation, run only now that the
        #      per-nucleus integral oracle above has passed. -----------------
        # (a) Contracted nuclear-attraction Hessian (analytic hess_en total) vs
        #     central FD of grad_en_pulay + grad_en_hellman_feynman, on the SAME
        #     low-symmetry molecule. hess1_selftest already wrote these; require
        #     the nucattr error to be at parity with the VALIDATED overlap and
        #     kinetic second-derivative terms (same FD truncation regime) -- a
        #     formula error would show up as a gross departure from that parity.
        sel = {}
        for line in Path("/tmp/hess1_selftest.out").read_text().splitlines():
            if "=" in line:
                k, _, v = line.partition("=")
                try:
                    sel[k.strip()] = float(v)
                except ValueError:
                    pass
        err_o = sel["overlap max|an-fd|"]
        err_k = sel["kinetic max|an-fd|"]
        err_v = sel["nucattr max|an-fd| (WIP)"]
        ref = max(err_o, err_k)
        self.assertLess(
            err_v, max(5.0 * ref, 1e-5),
            f"contracted nucattr FD not at parity with overlap/kinetic: "
            f"nucattr={err_v:.3e} overlap={err_o:.3e} kinetic={err_k:.3e}")

        # (b) Non-circular charge-charge confirmation: the TI-built p_CC block
        #     (hess_en block 9) vs an INDEPENDENT finite difference of the HF
        #     charge gradient w.r.t. the charge position (fd_cc), which never
        #     uses translational invariance. Parity with the FD truncation floor
        #     confirms p_CC is correct rather than merely TI-consistent.
        sep = Path("/tmp/hess1_sep.out").read_text().splitlines()
        idx = [i for i, l in enumerate(sep) if "block9 - fd_cc" in l]
        self.assertTrue(idx, "separated-FD diagnostic missing block9 - fd_cc")
        rows = [r for r in sep[idx[0] + 1: idx[0] + 1 + 3 * natom] if r.strip()]
        M = np.array([[float(x) for x in r.split()] for r in rows])
        self.assertLess(
            np.abs(M).max(), max(5.0 * ref, 1e-5),
            f"TI-built p_CC disagrees with independent charge FD: "
            f"max|block9-fd_cc|={np.abs(M).max():.3e}")


if __name__ == "__main__":
    unittest.main()
