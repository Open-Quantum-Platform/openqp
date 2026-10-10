#!/usr/bin/env python3
"""Prototype: consensus / ensemble SCF accelerator.

Premise: the Fock build dominates SCF cost and happens ONCE per iteration. From
that single Fock + history, several extrapolators each produce a cheap prediction
of the next step. Combining them (statistically) may converge in fewer Fock builds
(= fewer iterations) than any single accelerator.

We compare, on identical (mol, Fock engine, guess, conv) footing, iterations to
reach ||grad|| (orbital gradient) < tol:

  diis     : plain Pulay C-DIIS on the Fock              (baseline 1)
  newton   : damped diagonal-Hessian quasi-Newton in κ   (baseline 2, "SOSCF"-like)
  ensemble : per-iter convex blend of the DIIS step and the Newton step in
             ROTATION (κ) space (valid density preserved), with online,
             energy-decrease-weighted blend coefficients (the "statistical game")

Everything is built on PySCF integrals (get_hcore/get_ovlp/get_jk) so the Fock
engine is identical across methods; only the step prediction differs.

Run (zeus, pyscf in venv):
  python3 harness/ensemble_scf_proto.py --root .
"""
import argparse, os
import numpy as np
from pyscf import gto, scf

TOL = 1e-7          # convergence on orbital-gradient RMS
MAXIT = 100

_RHF = ["t0_h2o", "t0_ch4", "t0_nh3", "t0_benzene",
        "t1_glycine", "t1_c8", "t1_c16", "t1_caffeine", "t1_porphine"]
SUBSET = [(sid, 0, 0, xc) for sid in _RHF for xc in (None, "b3lyp5")]


def build(root, sid, charge, spin, xc):
    atoms = open(os.path.join(root, "geometries", f"{sid}.xyz")).read().splitlines()
    n = int(atoms[0].split()[0])
    mol = gto.M(atom="\n".join(atoms[2:2 + n]), charge=charge, spin=spin,
                basis="def2-svp", cart=True, verbose=0)
    mf = scf.RKS(mol) if xc else scf.RHF(mol)
    if xc:
        mf.xc = xc
        mf.grids.atom_grid = (96, 302)
        mf.grids.build()
    return mol, mf


def orb_grad(fock, dm, s, mo):
    """RHF orbital gradient g_ai = F_ai (MO basis), occ-vir block. RMS."""
    nocc = int(round(dm.trace().real @ np.eye(1)[0])) if False else None
    return None  # unused; gradient computed inline below


def occ_from(mo_energy, nocc):
    idx = np.argsort(mo_energy)
    occ = np.zeros_like(mo_energy)
    occ[idx[:nocc]] = 2.0
    return occ


def run(mf, mode, noise=0.0, seed=1):
    """One SCF to convergence using `mode`. Returns (iters, converged, E, gnorm)."""
    mol = mf.mol
    s = mf.get_ovlp()
    h = mf.get_hcore()
    nao = s.shape[0]
    nocc = mol.nelec[0]
    # symmetric orthogonalizer X = S^-1/2
    sval, svec = np.linalg.eigh(s)
    X = svec @ np.diag(sval ** -0.5) @ svec.T

    # initial guess density (core Hamiltonian)
    e, c = np.linalg.eigh(X.T @ h @ X); c = X @ c
    occ = occ_from(e, nocc); dm = (c * occ) @ c.T

    from scipy.linalg import expm
    diis = scf.diis.CDIIS()
    # orbitals that BUILT the current density (start: core guess)
    co = c[:, np.argsort(e)[:nocc]]
    cv = c[:, np.argsort(e)[nocc:]]
    # L-BFGS history for the SOSCF predictor (no response builds): (s=kappa step, y=grad change)
    Sh, Yh, MMEM = [], [], 8
    g_prev_flat = None; s_last = None
    # history-built approximate Fock response: free exact samples dF_k = G[dD_k]
    DH, FH, MEMR = [], [], 10
    D_prev = F_prev = None
    nao = s.shape[0]

    def rotate(co, cv, kap):
        C = np.hstack([co, cv]); no = co.shape[1]
        K = np.zeros((C.shape[1], C.shape[1]))
        K[no:, :no] = kap; K[:no, no:] = -kap.T
        Cn = C @ expm(K)
        return Cn[:, :no], Cn[:, no:]

    def lbfgs_step(gflat, eo, ev):
        """Two-loop L-BFGS inverse-Hessian step in kappa space (no response build).
        Initial Hessian = orbital-energy-difference diagonal (cheap, exact diagonal)."""
        q = gflat.copy(); alpha = []
        for s, y in zip(reversed(Sh), reversed(Yh)):
            rho = 1.0 / (y @ s + 1e-30); a = rho * (s @ q); alpha.append(a); q = q - a * y
        # H0: diagonal (ev - eo)
        d0 = (ev[:, None] - eo[None, :]).ravel()
        d0 = np.where(np.abs(d0) < 0.1, np.sign(d0) * 0.1 + (d0 == 0), d0)
        z = q / d0
        for (s, y), a in zip(zip(Sh, Yh), reversed(alpha)):
            rho = 1.0 / (y @ s + 1e-30); b = rho * (y @ z); z = z + (a - b) * s
        return -z

    for it in range(1, MAXIT + 1):
        fock = mf.get_fock(h1e=h, s1e=s, dm=dm)        # <-- the only expensive op
        E = mf.energy_tot(dm=dm, h1e=h, vhf=mf.get_veff(mol, dm))
        if noise > 0.0:                                 # simulate grid/screening noise floor
            rng = np.random.default_rng(seed + it)
            nz = rng.standard_normal(fock.shape); nz = noise * (nz + nz.T) / 2.0
            fock = fock + nz
        # harvest the FREE exact response sample dF = G[dD] (h is constant, G linear)
        if D_prev is not None:
            dDk = (dm - D_prev).ravel(); dFk = (fock - F_prev).ravel()
            if np.linalg.norm(dDk) > 1e-12:
                DH.append(dDk); FH.append(dFk)
                if len(DH) > MEMR:
                    DH.pop(0); FH.pop(0)
        D_prev = dm.copy(); F_prev = fock.copy()
        # orbital gradient measured in the basis that BUILT dm (nonzero until converged)
        g = cv.T @ fock @ co                            # (nvir, nocc)
        gnorm = np.sqrt(np.mean(g ** 2)) if g.size else 0.0
        if gnorm < TOL:
            return it, True, E, gnorm
        # orbital energies in the current basis (diagonal of F in co/cv)
        eo = np.einsum('pi,pq,qi->i', co, fock, co)
        ev = np.einsum('pa,pq,qa->a', cv, fock, cv)

        # ---- prediction 1: DIIS (Fock extrapolation) -> occupied orbitals ----
        f_diis = diis.update(s, dm, fock)
        ed, cd = np.linalg.eigh(X.T @ f_diis @ X); cd = X @ cd
        cod = cd[:, np.argsort(ed)[:nocc]]
        if mode == "diis":
            dm = 2.0 * cod @ cod.T
            co, cv = cod, cd[:, np.argsort(ed)[nocc:]]
            continue

        # ---- prediction 2: SOSCF = L-BFGS quasi-Newton in κ space (NO response build) ----
        g_flat = g.ravel()
        if g_prev_flat is not None and s_last is not None:
            y = g_flat - g_prev_flat
            if y @ s_last > 1e-12:                       # curvature condition -> store pair
                Sh.append(s_last.copy()); Yh.append(y.copy())
                if len(Sh) > MMEM:
                    Sh.pop(0); Yh.pop(0)
        g_prev_flat = g_flat.copy()
        kap_b = lbfgs_step(g_flat, eo, ev).reshape(g.shape)
        kn = np.linalg.norm(kap_b)
        if kn > 0.5:
            kap_b *= 0.5 / kn

        if mode == "newton":                            # pure SOSCF (L-BFGS), DIIS-free
            co, cv = rotate(co, cv, kap_b); dm = 2.0 * co @ co.T
            s_last = kap_b.ravel(); continue

        if mode == "histqn":                            # DIIS far out; history-response Newton in tail
            if gnorm > 1.0e-3:
                dm = 2.0 * cod @ cod.T
                co, cv = cod, cd[:, np.argsort(ed)[nocc:]]
            else:
                # approximate orbital-Hessian action using the reconstructed exact
                # Fock response G ~ sum_k dF_k <dD_k|.> (projection onto explored subspace)
                MD = np.array(DH).T if DH else None
                FD = np.array(FH).T if FH else None
                def Hvec(kv):
                    kp = kv.reshape(g.shape)
                    dD = 2.0 * (cv @ kp @ co.T + co @ kp.T @ cv.T)
                    if MD is not None:
                        cc = np.linalg.lstsq(MD, dD.ravel(), rcond=None)[0]
                        GdD = (FD @ cc).reshape(nao, nao)
                    else:
                        GdD = np.zeros((nao, nao))
                    return ((ev[:, None] - eo[None, :]) * kp + cv.T @ GdD @ co).ravel()
                # CG solve H kap = -g (diagonal-dominant near convergence)
                b = (-g).ravel(); x = np.zeros_like(b); r = b - Hvec(x); p = r.copy(); rs = r @ r
                for _ in range(20):
                    Hp = Hvec(p); a = rs / (p @ Hp + 1e-30); x = x + a * p; r = r - a * Hp
                    rs2 = r @ r
                    if rs2 < 1e-12 * len(b):
                        break
                    p = r + (rs2 / rs) * p; rs = rs2
                kap = x.reshape(g.shape)
                kn = np.linalg.norm(kap)
                if kn > 0.5:
                    kap *= 0.5 / kn
                co, cv = rotate(co, cv, kap); dm = 2.0 * co @ co.T
            continue

        if mode == "switch":                            # DIIS first, pure SOSCF later
            if gnorm > 1.0e-3:                          # DIIS handles the global part
                dm = 2.0 * cod @ cod.T
                co, cv = cod, cd[:, np.argsort(ed)[nocc:]]
                Sh.clear(); Yh.clear(); g_prev_flat = None; s_last = None
            else:                                       # SOSCF polishes the tail
                co, cv = rotate(co, cv, kap_b); dm = 2.0 * co @ co.T
                s_last = kap_b.ravel()
            continue

        # ---- ensemble: DIIS warmup far out, DIIS+SOSCF blend near convergence ----
        # Both predictors are 1 Fock/iter (no response builds). DIIS density used
        # directly while ||g|| large (robust, basin-finding); blend the DIIS and
        # SOSCF κ-steps once small, weight ramping toward SOSCF as ||g||->0.
        G_SWITCH = 1.0e-2
        if gnorm > G_SWITCH:
            dm = 2.0 * cod @ cod.T
            co, cv = cod, cd[:, np.argsort(ed)[nocc:]]
            Sh.clear(); Yh.clear(); g_prev_flat = None; s_last = None  # frame jumped -> reset BFGS
        else:
            kap_diis = cv.T @ s @ cod
            w = 1.0 / (1.0 + 10.0 ** (np.log10(max(gnorm, 1e-30)) + 4.0))  # ->1 (SOSCF) as g->0
            kap_mix = (1.0 - w) * kap_diis + w * kap_b
            co, cv = rotate(co, cv, kap_mix); dm = 2.0 * co @ co.T
            s_last = kap_mix.ravel()

    return MAXIT, False, E, gnorm


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", default="."); a = ap.parse_args()
    root = os.path.abspath(a.root)
    MODES = ("diis", "newton", "ensemble", "switch", "histqn")
    hdr = "%-26s" % "cell" + "".join("%12s" % m for m in MODES)
    print(hdr); print("-" * len(hdr))
    for sid, ch, sp, xc in SUBSET:
        cell = f"{sid}/{xc or 'hf'}"
        res = {}
        for mode in MODES:
            mol, mf = build(root, sid, ch, sp, xc)
            try:
                it, conv, E, g = run(mf, mode)
                res[mode] = f"{it}{'' if conv else '*'}"
            except Exception as ex:
                res[mode] = f"ERR:{type(ex).__name__}"
        print("%-26s" % cell + "".join("%12s" % res[m] for m in MODES), flush=True)
    print("\n(newton=pure SOSCF, switch=DIIS-then-SOSCF; * = not converged in %d iters; counts=Fock builds)" % MAXIT)


if __name__ == "__main__":
    main()
