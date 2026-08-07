"""Native CASPT2 / MS-CASPT2 / XMS-CASPT2 for OpenQP (determinant-based).

Clean rewrite of the CASPT2 path (the module filename is historical).  The
retired ``caspt2.py`` used an Epstein-Nesbet zeroth-order Hamiltonian which
over-correlates by 6-11 mEh versus OpenMolcas.  This module solves the first-order
amplitudes in the external determinant space Q (the complement of the CAS model
space) built from semicanonical orbitals:

    (H0_QQ - E0_I) C_I = -V_I,     V_I = (H |Psi_I^(0)>)_Q,     E_I^(2) = V_I . C_I.

The zeroth-order Hamiltonian H0 is selectable (``[pt2] h0``):

* ``fock`` (default) -- the one-electron generalized Fock (``diag(eps)`` in
  semicanonical orbitals, no two-electron part): the **CASPT2** H0.  With the
  default PT2 frozen core (``[pt2] frozen``, see below) it reproduces OpenMolcas
  CIonly CASPT2 (IPEA=0, no shift) to **<= 22 uEh** for both all-valence (H4, H6)
  and atomic-core (LiH, HF, H2O, BeH2) references.  ``[pt2] frozen=0`` instead
  correlates the deep 1s core (all-electron) and lies ~0.03-0.23 mEh below the
  frozen-core OpenMolcas value for atomic-core systems.
* ``dyall`` -- Dyall's H0 (core-dressed active one-electron + active two-electron),
  for which the reference is an exact H0 eigenstate: this is **NEVPT2** (Dyall H0),
  validatable against PySCF/ORCA (<= 5 uEh CAS(2,2) vs PySCF SC-NEVPT2).  With
  ``[pt2] contraction=strong`` this becomes RDM-based **SC-NEVPT2** (``nevpt2_sc``),
  reproducing PySCF SC-NEVPT2 to ~nEh.

Variants:
* **Single-state** (``method=caspt2``): E(CASPT2)_I = E(CASCI)_I + E_I^(2).
* **MS-CASPT2** (``method=ms-caspt2``): the multistate effective Hamiltonian
  ``H_eff[I,J] = <I|H|J> + 1/2(V_I.C_J + V_J.C_I)`` (symmetrized), diagonalized;
  the diagonal is the single-state CASPT2 energies.  For the default Fock H0 the
  **multi-set** construction (per-state orbitals = OpenMolcas state-specific Fock,
  matched to ~1.8 uEh) is used; the Dyall-H0 and XMS paths use the single-set
  (state-averaged-orbital) construction.
* **XMS-CASPT2** (``method=xms-caspt2``): rotate the references by diagonalizing
  the state-averaged Fock in the model space, then the same construction.

QDPT family (GAMESS-convention determinant MRPT2), same machinery re-labelled
and re-routed:
* **MRMP2** (``method=mrmp2``): single-state, diagonal-Fock (MP-like) H0 --
  identical to ``method=caspt2`` with the default ``h0=fock``.
* **MCQDPT2** (``method=mcqdpt2``): multistate QDPT with the SINGLE-SET
  state-averaged semicanonical orbitals (the GAMESS ``IFORB=1`` canonical Fock),
  one common ``diag(eps)`` H0 operator and STATE-SPECIFIC zeroth-order energies
  ``E0_I = <I|H0|I>`` (Nakano's ket-dependent denominators), symmetrized
  ``H_eff``.  This is the determinant analogue of GAMESS ``MRPT=DETMRPT`` /
  ``MRPT=MCQDPT`` (those two agree numerically for CAS references).
* **XMCQDPT2** (``method=xmcqdpt2``): Granovsky's extended QDPT (``XZERO=.TRUE.``):
  the references are first rotated to diagonalize the state-averaged Fock in the
  model space and the same common H0 is used -- the invariant variant.  This is
  the same rotation + single-set path as ``xms-caspt2`` here.
* ``[pt2] edshft`` (GAMESS ``EDSHFT``): intruder-state-avoidance denominator
  regularization ``d -> d + edshft/d`` (sign-preserving); mutually exclusive
  with the real/imaginary level shifts.

Real and imaginary level shifts regularize intruder states; the smallest
denominator is reported.  H and H0 are assembled in the native determinant basis
with :mod:`oqp.library.fci`; the reference comes from :mod:`oqp.library.casscf`
(or a CASCI in the supplied orbitals).

Scope: validation-grade.  Q is the full external determinant space, so this is
limited to small systems; ``[pt2] max_det`` guards the determinant count.  The
IPEA shift (``[pt2] ipea_shift``) and the PT2 frozen core (``[pt2] frozen``,
default = the standard deep cores) are applied.
"""
from __future__ import annotations

from dataclasses import dataclass
import time

import numpy as np

import oqp
from oqp.library.fci import (
    _as_f64c,
    _as_i64c,
    _build_dense_hamiltonian,
    _determinants,
    _spin_orbital_integrals,
    _transform_integrals,
    _unpack_lower_triangle,
    fci_spin_diagnostics,
    settings_from_casci_config,
    _symmetric_eigh,
)
from oqp.library.casscf import CASSCF, _solve_active
from oqp.utils.file_utils import print_module_banner

_BOOL_TRUE = {"true", "t", "yes", "y", "1", "on"}
_BOOL_FALSE = {"false", "f", "no", "n", "0", "off", ""}
_COUPLING_CUTOFF = 1.0e-12
_INTRUDER_DENOM = 0.1   # |denominator| below this flags a likely intruder state


def _intruder_note(min_denom, options):
    """Return a log annotation when a small denominator signals an intruder."""
    if min_denom >= _INTRUDER_DENOM:
        return ""
    if options.imaginary_shift or options.level_shift or options.edshft:
        return "   (small; a level shift is applied)"
    return "   <-- WARNING: likely intruder; set [pt2] imaginary_shift (e.g. 0.1)"


# --------------------------------------------------------------------------- options
@dataclass
class CASPT2Options:
    reference: str = "casscf"        # casscf | casci
    variant: str = "caspt2"          # caspt2 | ms | xms
    family: str = "caspt2"           # caspt2 | qdpt (MRMP2/MCQDPT2/XMCQDPT2 labels
                                     # + single-set routing + edshft ISA)
    h0: str = "fock"                 # fock (=CASPT2, default) | dyall (=NEVPT2)
    contraction: str = "none"        # none (uncontracted determinant space) |
                                     # strong (RDM-based SC-NEVPT2; dyall H0 only)
    root: int = 0                    # single-state reference root
    target_roots: tuple = ()         # multistate roots (defaults to range(nroot))
    nroot: int = 1                   # multistate root count when target_roots empty
    ipea_shift: float = 0.0          # 1e diagonal H0 IPEA shift (Ghigo-Roos-Malmqvist 2004)
    level_shift: float = 0.0         # real denominator shift
    imaginary_shift: float = 0.0     # imaginary (intruder) regularization
    edshft: float = 0.0              # GAMESS-style ISA: d -> d + edshft/d
    frozen: int = -1                 # PT2 frozen core: -1 = auto (standard deep cores,
                                     # matches OpenMolcas); 0 = correlate all; N = freeze N
    semi_canonical: bool = True
    max_det: int = 12000
    engine: str = "auto"             # auto (qdpt->direct) | direct | dense
    max_terms: int = 30_000_000      # direct-engine streamed-term guard
    nproc: int = 1                   # direct-engine process parallelism over
                                     # reference-determinant chunks


def _as_bool(value, label):
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in _BOOL_TRUE:
        return True
    if text in _BOOL_FALSE:
        return False
    raise ValueError(f"{label} must be a boolean")


def _as_int_tuple(value):
    if value is None:
        return ()
    if isinstance(value, (list, tuple)):
        items = value
    else:
        items = str(value).replace(",", " ").split()
    return tuple(int(v) for v in items if str(v).strip() != "")


def _caspt2_options(config: dict) -> CASPT2Options:
    raw = config.get("pt2", {}) or {}
    reference = str(raw.get("reference", "casscf")).strip().lower()
    if reference not in {"casci", "casscf"}:
        raise ValueError("pt2.reference must be 'casscf' or 'casci'")
    method = str(config.get("input", {}).get("method", "caspt2")).strip().lower()
    qdpt_methods = {"mrmp2": "caspt2", "mcqdpt2": "ms", "xmcqdpt2": "xms"}
    family = "qdpt" if method in qdpt_methods else "caspt2"
    variant_raw = str(raw.get("variant", "")).strip().lower()
    if variant_raw in {"", "auto"}:
        if method in qdpt_methods:
            variant = qdpt_methods[method]
        elif method in {"ms-caspt2", "mscaspt2"}:
            variant = "ms"
        elif method in {"xms-caspt2", "xmscaspt2"}:
            variant = "xms"
        else:
            variant = "caspt2"
    elif variant_raw in {"caspt2", "mrmp2"}:
        variant = "caspt2"
    elif variant_raw in {"ms", "ms-caspt2", "mscaspt2", "mcqdpt2"}:
        variant = "ms"
    elif variant_raw in {"xms", "xms-caspt2", "xmscaspt2", "xmcqdpt2"}:
        variant = "xms"
    else:
        raise ValueError("pt2.variant must be auto, caspt2, ms-caspt2, xms-caspt2, "
                         "mrmp2, mcqdpt2 or xmcqdpt2")
    h0_raw = str(raw.get("h0", "fock")).strip().lower()
    if h0_raw in {"fock", "caspt2", ""}:
        h0 = "fock"
    elif h0_raw in {"dyall", "nevpt2", "nevpt"}:
        h0 = "dyall"
    else:
        raise ValueError("pt2.h0 must be 'fock' (CASPT2) or 'dyall' (NEVPT2)")
    if family == "qdpt" and h0 != "fock":
        raise ValueError(
            "MRMP2/MCQDPT2/XMCQDPT2 are defined for the diagonal-Fock (MP-like) "
            "zeroth order only; drop [pt2] h0 or set h0=fock."
        )
    contraction_raw = str(raw.get("contraction", "none")).strip().lower()
    if contraction_raw in {"none", "uncontracted", "full", ""}:
        contraction = "none"
    elif contraction_raw in {"strong", "sc", "sc-nevpt2", "ic", "internally-contracted"}:
        contraction = "strong"
    else:
        raise ValueError("pt2.contraction must be 'none' (uncontracted) or 'strong' (SC-NEVPT2)")
    if contraction == "strong" and h0 != "dyall":
        raise ValueError(
            "pt2.contraction='strong' (SC-NEVPT2) requires h0='dyall'; "
            "the strong contraction is defined for Dyall's H0 only."
        )
    if family == "qdpt" and contraction != "none":
        raise ValueError("MRMP2/MCQDPT2/XMCQDPT2 are uncontracted; "
                         "drop [pt2] contraction.")
    edshft = float(raw.get("edshft", 0.0))
    if edshft < 0.0:
        raise ValueError("pt2.edshft (ISA denominator shift) must be >= 0")
    if edshft and (float(raw.get("level_shift", 0.0)) or float(raw.get("imaginary_shift", 0.0))):
        raise ValueError(
            "pt2.edshft (GAMESS-style ISA) is mutually exclusive with "
            "level_shift/imaginary_shift; choose one regularization."
        )
    engine = str(raw.get("engine", "auto")).strip().lower() or "auto"
    if engine == "python":
        engine = "direct"       # alias: NumPy streaming engine
    if engine not in {"auto", "direct", "fortran", "dense"}:
        raise ValueError("pt2.engine must be auto, fortran, direct (python), or dense")
    if engine in {"direct", "fortran"} and family != "qdpt":
        raise ValueError(
            "pt2.engine=direct/fortran requires the diagonal-Fock QDPT family "
            "(mrmp2/mcqdpt2/xmcqdpt2); the CASPT2/NEVPT2 H0 is not diagonal "
            "in the determinant basis."
        )
    return CASPT2Options(
        reference=reference,
        variant=variant,
        family=family,
        h0=h0,
        contraction=contraction,
        root=int(raw.get("root", raw.get("state", 0))),
        target_roots=_as_int_tuple(raw.get("target_roots", ())),
        # 0 is the documented default meaning "follow the CI"; taking it
        # literally made _reference_roots do max(2, 0) and silently solve two
        # states for any input with [ci] nroot > 2 -- including the compact
        # job.caspt2(variant="ms-caspt2", nroot=3) helpers, which write nroot
        # only to [ci].
        nroot=int(raw.get("nroot", 0) or 0) or int(
            (config.get("ci", {}) or {}).get("nroot", 0) or 0),
        ipea_shift=float(raw.get("ipea_shift", 0.0)),
        level_shift=float(raw.get("level_shift", 0.0)),
        imaginary_shift=float(raw.get("imaginary_shift", 0.0)),
        edshft=edshft,
        frozen=(-1 if str(raw.get("frozen", "auto")).strip().lower() in {"auto", "", "none", "-1"}
                else int(raw.get("frozen"))),
        semi_canonical=_as_bool(raw.get("semi_canonical", True), "pt2.semi_canonical"),
        # [pt2] max_det is NOT in OQP_CONFIG_SCHEMA, so an input naming it is
        # rejected by the parser -- yet _build_operators' own error message
        # tells the user to raise it.  Fall back to [cas] max_det, which is a
        # real key, so the guard can actually be lifted.
        max_det=int(raw.get("max_det", 0) or 0)
                or int((config.get("cas", {}) or {}).get("max_det", 0) or 0)
                or 12000,
        engine=engine,
        max_terms=int(raw.get("max_terms", 30_000_000)),
        nproc=max(1, int(raw.get("nproc", 1))),
    )


def _reference_roots(options) -> list:
    if options.variant == "caspt2":
        return [int(options.root)]
    if options.target_roots:
        return [int(r) for r in options.target_roots]
    nroot = max(2, int(options.nroot) or 0)
    return list(range(nroot))


# --------------------------------------------------------------------------- log
def _log(mol, text: str = "") -> None:
    with open(mol.log, "a") as handle:
        handle.write(text + "\n")


# --------------------------------------------------------------------------- native engine
def _pt2_lib():
    """(lib, ffi) of the loaded native core, or ``None`` when unavailable.

    Every caller below falls back to the NumPy/Python implementation that
    follows it when this returns ``None`` or the symbol is missing, so the
    module keeps working against a liboqp built before these kernels existed --
    and the Python stays the numerical pin.
    """
    try:
        from oqp.library.fci import _lib_backend
    except Exception:
        return None
    return _lib_backend()


def _dptr(ffi, arr):
    return ffi.cast("double *", arr.ctypes.data)


def _iptr(ffi, arr):
    return ffi.cast("int64_t *", arr.ctypes.data)


# --------------------------------------------------------------------------- Fock
def _effective_fock(h1e, eri, D):
    """Closed+active mean-field (generalized) Fock; its diagonal is the orbital
    energy used by Dyall's zeroth-order Hamiltonian.

    The liboqp engine (``pt2_effective_fock``) evaluates the same
    ``F = h + J - 1/2 K``; the NumPy einsums below remain as the fallback and
    the numerical pin."""
    backend = _pt2_lib()
    if backend is not None:
        lib, ffi = backend
        if hasattr(lib, "pt2_effective_fock"):
            nbf = int(h1e.shape[0])
            h = np.ascontiguousarray(h1e, dtype=np.float64)
            g = np.ascontiguousarray(eri, dtype=np.float64)
            d = np.ascontiguousarray(D, dtype=np.float64)
            fock = np.zeros((nbf, nbf), dtype=np.float64)
            lib.pt2_effective_fock(nbf, _dptr(ffi, h), _dptr(ffi, g),
                                   _dptr(ffi, d), _dptr(ffi, fock))
            return fock
    J = np.einsum("rs,pqrs->pq", D, eri, optimize=True)
    K = np.einsum("rs,prsq->pq", D, eri, optimize=True)
    return h1e + J - 0.5 * K


# --------------------------------------------------------------------------- semicanonical reference
def _semicanonicalize(hcore_ao, eri_ao, coeff, ncore, nact, active_nelec, enuc,
                      settings, roots, weights, max_iter=8):
    """Rotate the inactive/active/virtual blocks to diagonalize the (state-
    averaged) effective Fock, re-solving the CASCI each pass.  Returns the
    semicanonical orbitals, MO integrals, the converged CASCI (coeffs, energies,
    determinants), the orbital energies ``eps`` and the averaging density."""
    C = np.array(coeff, dtype=float, copy=True)
    nbf = C.shape[1]
    blocks = [range(ncore), range(ncore, ncore + nact), range(ncore + nact, nbf)]
    weights = np.asarray(weights, dtype=float)
    roots = [int(r) for r in roots]
    h1e = eri = energies = coeffs = dets = D = None
    for _it in range(max_iter):
        h1e, eri = _transform_integrals(hcore_ao, eri_ao, C)
        energies, coeffs, dets, D, _G = _solve_active(
            h1e, eri, ncore, nact, active_nelec, enuc, settings, weights, roots
        )
        F = _effective_fock(h1e, eri, D)
        U = np.eye(nbf)
        offdiag = 0.0
        for blk in blocks:
            idx = list(blk)
            if not idx:
                continue
            sub = F[np.ix_(idx, idx)]
            if len(idx) > 1:
                offdiag = max(offdiag, float(np.max(np.abs(sub - np.diag(np.diag(sub))))))
            _w, vec = _symmetric_eigh(0.5 * (sub + sub.T))
            U[np.ix_(idx, idx)] = vec
        C = C @ U
        if offdiag < 1.0e-10:
            break
    h1e, eri = _transform_integrals(hcore_ao, eri_ao, C)
    energies, coeffs, dets, D, _G = _solve_active(
        h1e, eri, ncore, nact, active_nelec, enuc, settings, weights, roots
    )
    eps = np.diag(_effective_fock(h1e, eri, D)).copy()
    return C, h1e, eri, coeffs, np.asarray(energies), dets, eps, D


# --------------------------------------------------------------------------- determinant operators
def _embed_reference(ci, act_dets, ncore, nact, norb, det_index):
    """Embed an active-space CI vector into the full (correlated) determinant
    space: core orbitals doubly occupied, active bits shifted into [ncore, ncore+nact)."""
    vec = np.zeros(len(det_index))
    core_mask = sum(1 << i for i in range(ncore))
    act_mask = (1 << nact) - 1
    for k, ad in enumerate(act_dets):
        a_full = core_mask | ((ad & act_mask) << ncore)
        b_full = core_mask | ((ad >> nact) << ncore)
        full = a_full | (b_full << norb)
        idx = det_index.get(full)
        if idx is None:
            raise ValueError("CASPT2 reference embedding failed (determinant missing)")
        vec[idx] = ci[k]
    return vec


def _external_indices(full_dets, ncore, nact, norb):
    """Indices of the external (first-order) determinants: every determinant that
    is not a pure CAS configuration (core doubly occupied, no virtual occupied).

    The liboqp engine (``pt2_external_indices``) walks the same test over the
    determinant keys; the Python loop below is the fallback and the pin."""
    backend = _pt2_lib() if 2 * int(norb) <= 62 else None
    if backend is not None:
        lib, ffi = backend
        if hasattr(lib, "pt2_external_indices"):
            dets = np.ascontiguousarray(np.asarray(full_dets, dtype=np.int64))
            ext = np.zeros(max(dets.size, 1), dtype=np.int64)
            n = int(lib.pt2_external_indices(
                int(norb), int(ncore), int(nact), int(dets.size),
                _iptr(ffi, dets), _iptr(ffi, ext)))
            return ext[:n].astype(int, copy=False)
    core_mask = sum(1 << i for i in range(ncore))
    virt_mask = sum(1 << i for i in range(ncore + nact, norb))
    full_mask = (1 << norb) - 1
    external = []
    for idx, det in enumerate(full_dets):
        a = det & full_mask
        b = det >> norb
        internal = (
            (a & core_mask) == core_mask and (b & core_mask) == core_mask
            and (a & virt_mask) == 0 and (b & virt_mask) == 0
        )
        if not internal:
            external.append(idx)
    return np.asarray(external, dtype=int)


def _h0_integrals(h1e, eri, eps, ncore, nact, h0, active_occ=None, ipea=0.0):
    """Zeroth-order one-/two-electron integrals for the chosen H0.

    ``fock`` -- the one-electron generalized Fock (diag(eps) in semicanonical
    orbitals, no two-electron part): this is the **CASPT2** zeroth-order
    Hamiltonian and matches OpenMolcas to <= 1e-4.

    ``dyall`` -- diag(eps) on core/virtual, the core-dressed active one-electron
    matrix ``h_tu + sum_i[2(tu|ii)-(ti|iu)]`` on the active block, plus the active
    two-electron operator: **Dyall's H0 (NEVPT2)**, for which the reference is an
    exact H0 eigenstate.
    """
    h0_1e = np.diag(eps).astype(float)
    g0_2e = np.zeros_like(eri)
    if h0 == "dyall":
        A = list(range(ncore, ncore + nact))
        h_act = None
        backend = _pt2_lib()
        if backend is not None:
            lib, ffi = backend
            if hasattr(lib, "pt2_h0_dyall_active"):
                nbf = int(h1e.shape[0])
                h = np.ascontiguousarray(h1e, dtype=np.float64)
                g = np.ascontiguousarray(eri, dtype=np.float64)
                h_act = np.zeros((nact, nact), dtype=np.float64)
                lib.pt2_h0_dyall_active(nbf, int(ncore), int(nact),
                                        _dptr(ffi, h), _dptr(ffi, g),
                                        _dptr(ffi, h_act))
        if h_act is None:               # fallback / numerical pin
            h_act = h1e[np.ix_(A, A)].copy()
            for p, P in enumerate(A):
                for q, Q in enumerate(A):
                    for i in range(ncore):
                        h_act[p, q] += 2.0 * eri[P, Q, i, i] - eri[P, i, i, Q]
        h0_1e[np.ix_(A, A)] = h_act
        g0_2e[np.ix_(A, A, A, A)] = eri[np.ix_(A, A, A, A)]
    # IPEA shift (Ghigo-Roos-Malmqvist 2004): bias each active diagonal of H0
    # toward ionization (n_P<1) or electron attachment (n_P>1) by 0.5*ipea per
    # half-occupation.  Early return is byte-identical when ipea == 0.
    if ipea:
        for p in range(nact):
            P = ncore + p
            n_P = float(active_occ[p])
            h0_1e[P, P] += 0.5 * ipea * (1.0 - n_P / 2.0)
    return h0_1e, g0_2e


def _core_orbitals_for_z(z):
    """Standard frozen-core orbital count for atomic number ``z`` (the deep, chemically
    inert shells OpenMolcas auto-freezes in CASPT2): 0 for H/He, 1 for Li-Ne (1s), ..."""
    z = int(round(float(z)))
    for limit, ncore in ((2, 0), (10, 1), (18, 5), (36, 9), (54, 18), (86, 27)):
        if z <= limit:
            return ncore
    return 43


def _standard_core_count(mol):
    """Number of deep atomic-core orbitals to freeze in the PT2 (sum over atoms).
    This matches OpenMolcas' default CASPT2 frozen core and removes the spurious deep
    1s core-hole modes that over-correlate atomic-core systems by 0.03-0.23 mEh."""
    try:
        atoms = np.asarray(mol.get_atoms(), dtype=float).ravel()
    except Exception:
        return 0
    return int(sum(_core_orbitals_for_z(z) for z in atoms))


def _freeze_core(h1e, eri, eps, D_sa, ncore, norb, enuc, nfrozen):
    """Fold the ``nfrozen`` deepest (semicanonical) core orbitals into a dressed closed
    shell -- exactly as ``casci`` dresses the inactive core -- and reduce the first-order
    problem to the remaining ``norb - nfrozen`` orbitals, so no external determinant can
    have a hole in a frozen orbital.  The reference energy is preserved (``enuc`` absorbs
    the frozen-core energy and ``h1e`` is dressed with the frozen-core mean field; the
    generalized-Fock diagonal ``eps`` is already correct for the remaining block); only
    the second-order space shrinks.  Returns reduced ``(h1e, eri, eps, D_sa, ncore, norb,
    enuc)``."""
    f = int(nfrozen)
    e_fc = 2.0 * float(np.einsum("ii->", h1e[:f, :f]))
    e_fc += float(np.einsum("iijj->", 2.0 * eri[:f, :f, :f, :f])
                  - np.einsum("ijji->", eri[:f, :f, :f, :f]))
    h1e_d = h1e[f:, f:].astype(float).copy()
    if f:
        h1e_d += np.einsum("pqii->pq", 2.0 * eri[f:, f:, :f, :f])
        h1e_d -= np.einsum("piiq->pq", eri[f:, :f, :f, f:])
    eri_d = eri[f:, f:, f:, f:].copy()
    eps_d = np.asarray(eps)[f:].copy()
    D_sa_d = None if D_sa is None else np.asarray(D_sa)[f:, f:].copy()
    return h1e_d, eri_d, eps_d, D_sa_d, ncore - f, norb - f, enuc + e_fc


def _diagonal_zeroth_order(full_dets, h0_1e, g0_2e, norb):
    """Determinant-basis diagonal of H0 when H0 is exactly diagonal, else ``None``.

    ``h0=fock`` (the default CASPT2/MRMP2/MCQDPT2/XMCQDPT2 zeroth order) makes H0
    a pure one-electron operator that is *already diagonal in the MO basis*:
    ``diag(eps)`` plus the IPEA bias, which only touches active diagonals, and no
    two-electron part.  Slater-Condon then gives ``<D|H0|D'> = 0`` for ``D != D'``
    -- a single excitation would need an off-diagonal ``h_pq`` and there is none --
    so the whole ``ndet x ndet`` matrix is carried by

        E0(D) = sum_p eps_p n_p(D).

    Measured on the shipped examples: ``max|offdiag(H0)| = 0.0`` exactly for
    ``h0=fock`` and ``2.2e-01`` for ``h0=dyall`` (which has an active two-electron
    part), so the test below is what selects the path, not the option string.
    """
    if 2 * norb > 62:                       # determinant encoding must fit int64
        return None
    if g0_2e is not None and np.any(g0_2e):
        return None
    eps0 = np.asarray(np.diag(h0_1e), dtype=float)
    if np.any(np.asarray(h0_1e) - np.diag(eps0)):
        return None
    dets = np.ascontiguousarray(np.asarray(full_dets, dtype=np.int64))
    # The two exactness guards above are what SELECT this path, and they stay in
    # Python; only the evaluation of the resulting diagonal moves down.  The
    # native sum runs over p in the same ascending order and adds eps_p once per
    # occupied spin, so it is bit-identical to the NumPy accumulation below
    # (eps+eps and 2*eps agree exactly in IEEE double).
    backend = _pt2_lib()
    if backend is not None:
        lib, ffi = backend
        if hasattr(lib, "pt2_diagonal_h0"):
            e0 = np.ascontiguousarray(eps0, dtype=np.float64)
            diag = np.zeros(max(dets.size, 1), dtype=np.float64)
            lib.pt2_diagonal_h0(int(norb), int(dets.size), _iptr(ffi, dets),
                                _dptr(ffi, e0), _dptr(ffi, diag))
            return diag[:dets.size]
    alpha = dets & np.int64((1 << norb) - 1)
    beta = dets >> np.int64(norb)
    diag = np.zeros(dets.size, dtype=float)
    for p in range(norb):
        bit = np.int64(1 << p)
        occ = ((alpha & bit) != 0).astype(np.float64)
        occ += ((beta & bit) != 0).astype(np.float64)
        diag += float(eps0[p]) * occ
    return diag


def _occupation_blocks(full_dets, external, ncore, nact, norb):
    """Partition the external determinants into H0-invariant occupation blocks.

    Dyall's H0 is ``sum_p eps_p n_p`` on the inactive/virtual orbitals plus an
    operator that lives entirely *inside* the active window, so it cannot move an
    electron in or out of a core or virtual orbital: it conserves the full
    core+virtual occupation pattern, per spin.  Determinants sharing that pattern
    therefore form a closed block and ``H0_ext`` is exactly block-diagonal over
    them -- measured on H4/6-31G NEVPT2: 120 blocks of at most 24, off-block
    ``||.||_F^2 = 0.000e+00`` exactly, and the denominator eigendecomposition
    drops from 4.2e8 to 1.8e5 flops.

    Returns a list of index arrays *into the external space*, or ``None`` when
    the determinant encoding does not fit int64.
    """
    if 2 * norb > 62:
        return None
    # liboqp engine; the NumPy partition below is the fallback and the pin.  The
    # native sort is a stable LSD radix sort on the same signature, so it
    # reproduces ``np.argsort(kind="stable")`` exactly -- the intra-block member
    # order, and hence each block's eigenproblem, is bit-for-bit identical.
    backend = _pt2_lib()
    if backend is not None:
        lib, ffi = backend
        if hasattr(lib, "pt2_occupation_blocks"):
            dets = np.ascontiguousarray(np.asarray(full_dets, dtype=np.int64))
            ext = np.ascontiguousarray(np.asarray(external, dtype=np.int64))
            nx = int(ext.size)
            order = np.zeros(max(nx, 1), dtype=np.int64)
            starts = np.zeros(nx + 1, dtype=np.int64)
            nblock = int(lib.pt2_occupation_blocks(
                int(norb), int(ncore), int(nact), nx, _iptr(ffi, dets),
                _iptr(ffi, ext), _iptr(ffi, order), _iptr(ffi, starts)))
            if nx == 0:
                return [order[:0]]
            return [order[starts[b]:starts[b + 1]] for b in range(nblock)]
    act_mask = ((1 << nact) - 1) << ncore
    frozen = ((1 << norb) - 1) & ~act_mask
    keep = np.int64(frozen | (frozen << norb))     # both spins, one mask
    sig = np.asarray(full_dets, dtype=np.int64)[external] & keep
    order = np.argsort(sig, kind="stable")
    edges = np.flatnonzero(np.diff(sig[order])) + 1
    return np.split(order, edges)


class _ZerothOrder:
    """The PT2 zeroth-order Hamiltonian, in the cheapest *exact* representation.

    Three facts turn what used to be one ``O(n_ext^3)`` LAPACK diagonalization
    per reference root into (almost) nothing:

    * A one-electron H0 that is diagonal in the MO basis is exactly diagonal in
      the determinant basis (see :func:`_diagonal_zeroth_order`).  Then the
      denominator operator needs no eigendecomposition at all and the dense
      ``ndet x ndet`` H0 never has to be built.
    * Dyall H0 carries an active two-electron operator and is not diagonal, but
      it *is* exactly block-diagonal over the core+virtual occupation patterns
      (see :func:`_occupation_blocks`), so the one big diagonalization becomes
      many tiny ones.  The partition is verified numerically before it is used.
    * Either way the denominator ``H0_ext - E0_I`` differs between roots only by
      a *scalar*, so it shares its eigenvectors with ``H0_ext`` and its
      eigenvalues are just shifted by ``-E0_I``.  One eigendecomposition
      therefore serves every root instead of one per root.

    The single full-block eigendecomposition survives as the fallback and the
    numerical pin: it is the original code path, and :func:`_first_order` runs
    the same arithmetic whichever representation is in force.
    """

    def __init__(self, external, diagonal=None, dense=None, blocks=None):
        self.external = external
        self.diagonal = diagonal
        self.dense = dense
        self.blocks = blocks
        self._eigenbasis = None

    def expectation(self, vec):
        """``<vec|H0|vec>`` for a normalized determinant-space vector."""
        if self.diagonal is not None:
            return float(vec @ (self.diagonal * vec))
        return float(vec @ self.dense @ vec)

    def external_eigenbasis(self):
        """``(w, transforms)`` for the *unshifted* external block.

        ``w`` holds the eigenvalues at their own external positions; each entry
        of ``transforms`` is ``(idx, U)`` mapping that block's eigenbasis back to
        those positions, with ``U is None`` meaning the identity.  All of it is
        root-independent, so it is built once and cached."""
        if self._eigenbasis is None:
            self._eigenbasis = self._build_eigenbasis()
        return self._eigenbasis

    def _build_eigenbasis(self):
        n = len(self.external)
        if self.diagonal is not None:                 # already diagonal
            return self.diagonal[self.external], ((slice(None), None),)

        blk = self.dense[np.ix_(self.external, self.external)]
        blk = 0.5 * (blk + blk.T)
        if self.blocks is not None and len(self.blocks) > 1:
            # The partition is only used once it is *shown* to close: the largest
            # element outside the claimed blocks must vanish.  Measured directly
            # (no cancelling sums, so this is exact) rather than assumed from the
            # h0 option string; anything else falls through to the full solve.
            scale = max(float(np.max(np.abs(blk))), 1.0) if blk.size else 1.0
            off = 0.0
            for g in self.blocks:
                rows = np.abs(blk[g, :])
                rows[:, g] = 0.0
                off = max(off, float(rows.max()) if rows.size else 0.0)
            if off <= 1.0e-10 * scale:
                w = np.empty(n, dtype=float)
                transforms = []
                for g in self.blocks:
                    wg, Ug = _symmetric_eigh(blk[np.ix_(g, g)])
                    w[g] = wg
                    transforms.append((g, Ug))
                return w, tuple(transforms)

        w, U = _symmetric_eigh(blk)
        return w, ((slice(None), U),)


def _build_operators(h1e, eri, eps, ncore, nact, active_nelec, norb, max_det, h0="fock",
                     active_occ=None, ipea=0.0):
    """Assemble the full electronic Hamiltonian, the H0 operator and the external space."""
    full_nelec = (ncore + active_nelec[0], ncore + active_nelec[1])
    full_dets = _determinants(norb, full_nelec)
    ndet = len(full_dets)
    if ndet > max_det:
        raise ValueError(
            f"CASPT2 uncontracted determinant space too large: ndet={ndet} > "
            f"[pt2] max_det={max_det}.  Reduce the basis/active space; an RDM-based "
            "contraction (larger systems) is a separate step."
        )
    det_index = {det: i for i, det in enumerate(full_dets)}

    hspin, gspin = _spin_orbital_integrals(h1e, eri)
    hfull = _build_dense_hamiltonian(full_dets, det_index, hspin, gspin, 2 * norb, 0.0)

    h0_1e, g0_2e = _h0_integrals(h1e, eri, eps, ncore, nact, h0,
                                 active_occ=active_occ, ipea=ipea)
    external = _external_indices(full_dets, ncore, nact, norb)

    # A diagonal H0 (the default Fock zeroth order) needs neither the dense
    # ndet x ndet build nor the O(n_ext^3) denominator diagonalization below.
    diag0 = _diagonal_zeroth_order(full_dets, h0_1e, g0_2e, norb)
    if diag0 is not None:
        h0op = _ZerothOrder(external, diagonal=diag0)
    else:
        hspin0, gspin0 = _spin_orbital_integrals(h0_1e, g0_2e)
        h0op = _ZerothOrder(
            external,
            dense=_build_dense_hamiltonian(full_dets, det_index, hspin0, gspin0,
                                           2 * norb, 0.0),
            blocks=_occupation_blocks(full_dets, external, ncore, nact, norb))

    return full_dets, det_index, hfull, h0op, external


def _first_order(hfull, h0op, vec, external, options):
    """First-order amplitudes and second-order energy for one reference vector.

    Works in the eigenbasis of the ``(H0 - E0)`` denominator operator so the real
    (Roos-Andersson) and imaginary (Forsberg-Malmqvist) level shifts and the
    intruder diagnostic share a single path.  ``h0op`` (a :class:`_ZerothOrder`)
    owns that eigenbasis: it is independent of the root, because ``H0_ext - E0_I``
    differs from ``H0_ext`` by a scalar, so it is built once and shared across
    states -- and it is block-diagonal (Dyall) or the identity (Fock) rather than
    a dense rotation.  Returns ``(V, C, e2, min_denom)`` where ``min_denom`` is
    the smallest absolute bare denominator (intruder probe)."""
    e0 = h0op.expectation(vec)
    V = (hfull @ vec)[external]
    w0, transforms = h0op.external_eigenbasis()
    w = w0 - e0                                             # bare CASPT2 denominators
    Vt = np.empty_like(V)
    for idx, U in transforms:
        Vt[idx] = V[idx] if U is None else U.T @ V[idx]
    d = w + options.level_shift                             # real level shift
    if options.edshft:                                     # GAMESS ISA: d -> d + edshft/d
        with np.errstate(divide="ignore", invalid="ignore"):
            d = d + options.edshft / d
    if options.imaginary_shift:                            # imaginary level shift d -> d/(d^2+s^2)
        factor = d / (d * d + options.imaginary_shift ** 2)
    else:
        factor = np.divide(1.0, d, out=np.zeros_like(d),
                           where=(np.abs(d) > 1.0e-300) & np.isfinite(d))
    Ct = -(Vt * factor)
    C = np.empty_like(Ct)
    for idx, U in transforms:
        C[idx] = Ct[idx] if U is None else U @ Ct[idx]
    min_denom = float(np.min(np.abs(w))) if w.size else float("inf")
    return V, C, float(V @ C), min_denom


# --------------------------------------------------------------------------- multi-set MS-CASPT2
def _occ_tuples(norb, nelec):
    """Occupied-orbital index tuples for each single-spin string, in the SAME
    order as :func:`oqp.library.fci._string_basis` (lexicographic ``combinations``)."""
    from itertools import combinations
    return [tuple(occ) for occ in combinations(range(norb), nelec)]


def _minor_transform(R, occ):
    """T[t, s] = det( R[ix_(t_occ, s_occ)] ): the single-spin orbital-rotation
    transform on the string space (Loewdin minors), ordered by ``occ``.

    The liboqp engine (``pt2_minor_transform``) evaluates the same minors with
    the determinant inlined instead of one ``np.linalg.det`` call per string
    PAIR, which is what made this the most expensive Python on the PT2 path;
    the loop below stays as the fallback and the numerical pin."""
    n = len(occ)
    backend = _pt2_lib()
    if backend is not None and n:
        lib, ffi = backend
        if hasattr(lib, "pt2_minor_transform"):
            nelec = len(occ[0])
            r = _as_f64c(R)
            table = _as_i64c(np.asarray(occ, dtype=np.int64).reshape(n, nelec))
            T = np.zeros((n, n), dtype=np.float64)
            lib.pt2_minor_transform(int(R.shape[0]), int(nelec), int(n),
                                    _dptr(ffi, r), _iptr(ffi, table),
                                    _dptr(ffi, T))
            return T
    T = np.zeros((n, n))
    for it, t in enumerate(occ):
        rows = list(t)
        for js, s in enumerate(occ):
            T[it, js] = np.linalg.det(R[np.ix_(rows, list(s))])
    return T


def _rotate_det_vector(vec, R, norb, nelec):
    """Transform a determinant-space vector built in MO set ``src`` into MO set
    ``dst`` where ``C_src = C_dst @ R`` (R unitary).  ``vec`` is indexed over
    :func:`_determinants` (alpha-major); the transform factorizes into the alpha
    and beta single-spin minor matrices of ``R``.  Determinant-only; reproduces
    pyscf.fci.addons.transform_ci_for_orbital_rotation to machine precision."""
    occ_a = _occ_tuples(norb, nelec[0])
    occ_b = _occ_tuples(norb, nelec[1])
    Ta = _minor_transform(R.T, occ_a)         # R.T: matches the C_src=C_dst@R convention
    Tb = _minor_transform(R.T, occ_b)
    na, nb = len(occ_a), len(occ_b)
    return (Ta @ vec.reshape(na, nb) @ Tb.T).reshape(-1)


def _h0_fock_matrix_dets(h1e, eri, D_state, full_dets, det_index, norb):
    """Full generalized-Fock MATRIX zeroth-order Hamiltonian for a per-state (multi-
    set) reference: F_pq = h_pq + sum_rs D_rs[(pq|rs)-1/2(pr|sq)] from the STATE'S
    OWN 1-RDM ``D_state``, assembled as a pure one-electron determinant operator.
    Unlike the diagonal Fock (``diag(eps)``), it keeps the off-diagonal inactive-
    active / active-virtual Fock couplings that are nonzero for a state-specific
    reference -- essential for strongly correlated roots (31 mEh error otherwise)."""
    F = _effective_fock(h1e, eri, D_state)
    hspin0, gspin0 = _spin_orbital_integrals(F, np.zeros_like(eri))
    return _build_dense_hamiltonian(full_dets, det_index, hspin0, gspin0, 2 * norb, 0.0)


def _active_ci_from_full(vec, det_index, act_dets, ncore, nact, norb):
    """Inverse of :func:`_embed_reference`: pull the active-space CI vector (flat,
    indexed by ``act_dets`` = ``_determinants(nact, active_nelec)``) out of a full
    determinant vector (core doubly occupied, no virtual occupied)."""
    core_mask = sum(1 << i for i in range(ncore))
    act_mask = (1 << nact) - 1
    out = np.zeros(len(act_dets))
    for k, ad in enumerate(act_dets):
        a_full = core_mask | ((ad & act_mask) << ncore)
        b_full = core_mask | ((ad >> nact) << ncore)
        out[k] = vec[det_index[a_full | (b_full << norb)]]
    return out


def _lift_to_full(vec_frozen, det_index_frozen, det_index_full, nfrozen,
                  norbf, norb, active_nelec, ncf, ncore):
    """Lift a determinant vector from the FROZEN orbital set (norbf, cores ncf) into
    the FULL determinant set (norb): prepend ``nfrozen`` doubly-occupied frozen
    orbitals (lowest indices), shift every orbital up by ``nfrozen``.  Frozen-core-
    hole determinants stay zero (frozen core not correlated)."""
    fcore = sum(1 << i for i in range(nfrozen))
    mask_f = (1 << norbf) - 1
    out = np.zeros(len(det_index_full))
    inv = {idx: det for det, idx in det_index_frozen.items()}
    for kf, det in inv.items():
        v = vec_frozen[kf]
        if v == 0.0:
            continue
        a = det & mask_f
        b = det >> norbf
        a_full = fcore | (a << nfrozen)
        b_full = fcore | (b << nfrozen)
        out[det_index_full[a_full | (b_full << norb)]] = v
    return out


def _multiset_spin_diagnostics(hcore_ao, eri_ao, coeff_ref, ncore, nact,
                               active_nelec, enuc, settings, roots):
    """<S^2>/multiplicity of the reference roots (CASCI in the reference orbitals)."""
    weights = np.full(len(roots), 1.0 / len(roots))
    h1e_r, eri_r = _transform_integrals(hcore_ao, eri_ao, coeff_ref)
    _e, coeffs, dets, _D, _G = _solve_active(
        h1e_r, eri_r, ncore, nact, active_nelec, enuc, settings, weights, roots)
    return fci_spin_diagnostics(coeffs, dets, nact, active_nelec)


def _multistate_multiset(mol, hcore_ao, eri_ao, coeff_ref, settings, ncore, nact,
                         active_nelec, nbf, enuc, roots, options, nfrozen):
    """Multi-set MS-CASPT2 (per-state orbitals = OpenMolcas state-specific Fock).

    For each root I: build the generalized Fock from I's OWN active 1-RDM,
    block-semicanonicalize -> per-state orbitals C_I / orbital energies; re-express
    I's CI vector in C_I; solve I's first-order CASPT2 in its own basis with the
    deep core frozen and the full generalized-Fock-MATRIX H0.  The effective
    Hamiltonian off-diagonal Heff[I,J] = 1/2(<Psi_I^(1)|H|Psi_J^(0)> +
    <Psi_J^(1)|H|Psi_I^(0)>) is evaluated in the unfrozen determinant space with
    J's reference rotated into basis I by R = C_I^T S_ao C_J (here U_I^T U_J, since
    both C_I = C_ref @ U_I share the AO metric)."""
    # Note: the inter-state orbital rotation R = C_I^T S_ao C_J reduces to U_I^T U_J
    # here, because C_I = C_ref @ U_I and C_ref is orthonormal in the AO metric S_ao
    # (so C_I^T S_ao C_J = U_I^T (C_ref^T S_ao C_ref) U_J = U_I^T U_J); no explicit
    # S_ao is needed.
    weights = np.full(len(roots), 1.0 / len(roots))

    # reference CI vectors in the supplied orbitals (CASCI in C_ref)
    h1e_r, eri_r = _transform_integrals(hcore_ao, eri_ao, coeff_ref)
    energies, coeffs, dets, _D, _G = _solve_active(
        h1e_r, eri_r, ncore, nact, active_nelec, enuc, settings, weights, roots)
    E_ref = np.array([float(energies[r]) for r in roots])

    from oqp.library.rdm import make_rdm1_spatial
    blocks = [range(ncore), range(ncore, ncore + nact), range(ncore + nact, nbf)]

    full_nelec_un = (ncore + active_nelec[0], ncore + active_nelec[1])
    dets_un = _determinants(nbf, full_nelec_un)
    didx_un = {d: i for i, d in enumerate(dets_un)}

    per = []
    for r in roots:
        gamma = make_rdm1_spatial(coeffs[:, r], dets, nact)
        D_state = np.zeros((nbf, nbf))
        for i in range(ncore):
            D_state[i, i] = 2.0
        D_state[np.ix_(range(ncore, ncore + nact), range(ncore, ncore + nact))] = gamma
        # block-semicanonicalize from this state's own density
        F = _effective_fock(h1e_r, eri_r, D_state)
        U = np.eye(nbf)
        for blk in blocks:
            idx = list(blk)
            if not idx:
                continue
            sub = F[np.ix_(idx, idx)]
            _w, vec = _symmetric_eigh(0.5 * (sub + sub.T))
            U[np.ix_(idx, idx)] = vec
        coeff_I = coeff_ref @ U
        h1e_I, eri_I = _transform_integrals(hcore_ao, eri_ao, coeff_I)
        Ua = U[np.ix_(range(ncore, ncore + nact), range(ncore, ncore + nact))]
        gamma_I = Ua.T @ gamma @ Ua
        # I's reference CI rotated into C_I (flat, indexed by the active dets ``dets``)
        ref_un = _embed_reference(coeffs[:, r], dets, ncore, nact, nbf, didx_un)
        ref_un_I = _rotate_det_vector(ref_un, U, nbf, full_nelec_un)
        ci_active_I = _active_ci_from_full(ref_un_I, didx_un, dets, ncore, nact, nbf)

        # --- diagonal: I's first-order CASPT2 in its own basis, frozen core, full-Fock H0
        h1f, erif, _eps, _Dsa_f, ncf, norbf, enucf = _freeze_core(
            h1e_I, eri_I, np.diag(_effective_fock(h1e_I, eri_I, D_state)),
            D_state, ncore, nbf, enuc, nfrozen)
        full_nelec = (ncf + active_nelec[0], ncf + active_nelec[1])
        dets_f = _determinants(norbf, full_nelec)
        didx_f = {d: i for i, d in enumerate(dets_f)}
        hspin, gspin = _spin_orbital_integrals(h1f, erif)
        Hf_f = _build_dense_hamiltonian(dets_f, didx_f, hspin, gspin, 2 * norbf, 0.0)
        Hf_f = Hf_f + enucf * np.eye(len(dets_f))   # constant: enuc + frozen-core energy
        D_active_f = np.zeros((norbf, norbf))
        for i in range(ncf):
            D_active_f[i, i] = 2.0
        D_active_f[np.ix_(range(ncf, ncf + nact), range(ncf, ncf + nact))] = gamma_I
        H0_f = _h0_fock_matrix_dets(h1f, erif, D_active_f, dets_f, didx_f, norbf)
        ext = _external_indices(dets_f, ncf, nact, norbf)
        vec_f = _embed_reference(ci_active_I, dets, ncf, nact, norbf, didx_f)
        e0 = float(vec_f @ H0_f @ vec_f)
        V = (Hf_f @ vec_f)[ext]
        denom = H0_f[np.ix_(ext, ext)] - e0 * np.eye(len(ext))
        w, Uq = _symmetric_eigh(0.5 * (denom + denom.T))
        Vt = Uq.T @ V
        d = w + options.level_shift
        if options.edshft:                                 # GAMESS ISA: d -> d + edshft/d
            with np.errstate(divide="ignore", invalid="ignore"):
                d = d + options.edshft / d
        if options.imaginary_shift:
            factor = d / (d * d + options.imaginary_shift ** 2)
        else:
            factor = np.divide(1.0, d, out=np.zeros_like(d),
                               where=(np.abs(d) > 1.0e-300) & np.isfinite(d))
        C1 = -(Uq @ (Vt * factor))
        e2 = float(V @ C1)
        min_denom = float(np.min(np.abs(w))) if w.size else float("inf")

        # --- full (unfrozen) objects for the cross-coupling (R is exactly unitary)
        hspin_u, gspin_u = _spin_orbital_integrals(h1e_I, eri_I)
        Hf_full = _build_dense_hamiltonian(dets_un, didx_un, hspin_u, gspin_u, 2 * nbf, 0.0)
        Hf_full = Hf_full + enuc * np.eye(len(dets_un))
        ref_full = _embed_reference(ci_active_I, dets, ncore, nact, nbf, didx_un)
        psi1_f = np.zeros(len(dets_f)); psi1_f[ext] = C1
        psi1_full = _lift_to_full(psi1_f, didx_f, didx_un, nfrozen, norbf, nbf,
                                  active_nelec, ncf, ncore)
        per.append(dict(r=r, U=U, e2=e2, min_denom=min_denom, ci_active=ci_active_I,
                        Hf_full=Hf_full, ref_full=ref_full, psi1_full=psi1_full))

    nstate = len(roots)
    heff = np.zeros((nstate, nstate))
    for i in range(nstate):
        heff[i, i] = E_ref[i] + per[i]["e2"]
    for i in range(nstate):
        for j in range(nstate):
            if i == j:
                continue
            Rij = per[i]["U"].T @ per[j]["U"]
            vecj_in_i = _rotate_det_vector(per[j]["ref_full"], Rij, nbf, full_nelec_un)
            cij = float(per[i]["psi1_full"] @ (per[i]["Hf_full"] @ vecj_in_i))
            Rji = per[j]["U"].T @ per[i]["U"]
            veci_in_j = _rotate_det_vector(per[i]["ref_full"], Rji, nbf, full_nelec_un)
            cji = float(per[j]["psi1_full"] @ (per[j]["Hf_full"] @ veci_in_j))
            heff[i, j] = 0.5 * (cij + cji)
    heff = 0.5 * (heff + heff.T)
    ms_energies, mixing = _symmetric_eigh(heff)
    ss_energies = np.array([heff[i, i] for i in range(nstate)])
    return {
        "ms_energies": ms_energies,
        "mixing": mixing,
        "heff": heff,
        "ss_energies": ss_energies,
        "ref_energies": E_ref,
        "e2": np.array([per[i]["e2"] for i in range(nstate)]),
        "min_denoms": np.array([per[i]["min_denom"] for i in range(nstate)]),
        "n_external": len(per[0]["psi1_full"]) if per else 0,
        "rotated": False,
        "ref_drift": 0.0,
        "multiset": True,
    }


def _xms_rotation(h1e, eri, D_sa, coeffs, dets, ncore, nact, norb, roots, det_index, hfull):
    """XMS pre-rotation: diagonalize the state-averaged Fock in the model space.
    Returns the rotation matrix R (columns = rotated states) over the roots."""
    fock = _effective_fock(h1e, eri, D_sa)
    fspin, gzero = _spin_orbital_integrals(fock, np.zeros_like(eri))
    full_dets = list(det_index.keys())
    fock_det = _build_dense_hamiltonian(full_dets, det_index, fspin, gzero, 2 * norb, 0.0)
    refs = [_embed_reference(coeffs[:, r], dets, ncore, nact, norb, det_index) for r in roots]
    nstate = len(roots)
    fmodel = np.zeros((nstate, nstate))
    for i in range(nstate):
        fi = fock_det @ refs[i]
        for j in range(nstate):
            fmodel[i, j] = refs[j] @ fi
    _w, R = np.linalg.eigh(0.5 * (fmodel + fmodel.T))
    return R


def _multistate(h1e, eri, coeffs, energies, dets, eps, D_sa, ncore, nact,
                active_nelec, norb, enuc, roots, options):
    """MS-/XMS-CASPT2 effective Hamiltonian over the reference roots."""
    active_occ = np.diag(D_sa)[ncore:ncore + nact]
    full_dets, det_index, hfull, h0op, external = _build_operators(
        h1e, eri, eps, ncore, nact, active_nelec, norb, options.max_det, options.h0,
        active_occ=active_occ, ipea=options.ipea_shift
    )

    refs = [_embed_reference(coeffs[:, r], dets, ncore, nact, norb, det_index) for r in roots]
    rotated = options.variant == "xms"
    if rotated:
        R = _xms_rotation(h1e, eri, D_sa, coeffs, dets, ncore, nact, norb, roots, det_index, hfull)
        refs = [sum(R[k, i] * refs[k] for k in range(len(roots))) for i in range(len(roots))]

    nstate = len(roots)
    # model-space Hamiltonian <Psi_I|H|Psi_J> (+ E_nuc on the diagonal)
    hmodel = np.zeros((nstate, nstate))
    for i in range(nstate):
        hi = hfull @ refs[i]
        for j in range(nstate):
            hmodel[i, j] = refs[j] @ hi
        hmodel[i, i] += enuc

    Vs, Cs, e2s, min_denoms = [], [], [], []
    ref_drift = 0.0
    for i in range(nstate):
        V, C, e2, min_denom = _first_order(hfull, h0op, refs[i], external, options)
        Vs.append(V)
        Cs.append(C)
        e2s.append(e2)
        min_denoms.append(min_denom)
        ref_drift = max(ref_drift, abs((refs[i] @ hfull @ refs[i]) + enuc - hmodel[i, i]))

    heff = hmodel.copy()
    for i in range(nstate):
        for j in range(nstate):
            heff[i, j] += 0.5 * (Vs[i] @ Cs[j] + Vs[j] @ Cs[i])
    heff = 0.5 * (heff + heff.T)
    ms_energies, mixing = _symmetric_eigh(heff)

    ss_energies = np.array([hmodel[i, i] + e2s[i] for i in range(nstate)])
    return {
        "ms_energies": ms_energies,
        "mixing": mixing,
        "heff": heff,
        "ss_energies": ss_energies,
        "ref_energies": np.array([hmodel[i, i] for i in range(nstate)]),
        "e2": np.array(e2s),
        "min_denoms": np.array(min_denoms),
        "n_external": len(external),
        "rotated": rotated,
        "ref_drift": ref_drift,
    }


# --------------------------------------------------------------------------- driver
def native_caspt2_energy(mol, ref_energy=None):
    """CASPT2 / MS-CASPT2 / XMS-CASPT2 energy (uncontracted Dyall H0)."""
    options = _caspt2_options(mol.config)
    settings = settings_from_casci_config(mol.config)

    if options.contraction == "strong" and options.variant != "caspt2":
        raise NotImplementedError(
            "[pt2] contraction='strong' (SC-NEVPT2) is single-state only; "
            "multistate (MS/XMS) strong contraction is future work."
        )
    # The IPEA shift is applied as a one-electron diagonal bias on the active H0
    # block (see _h0_integrals); the schema default is 0.0 so existing benchmarks
    # are unaffected unless the user explicitly requests a shift.
    if settings.integral_backend != "native":
        raise ValueError("caspt2 integral_backend must be native")
    if mol.config["input"]["functional"]:
        raise ValueError("CASPT2 requires HF integrals; unset input.functional")
    if mol.config["scf"]["type"] != "rhf":
        raise ValueError("CASPT2 currently supports a closed-shell RHF reference")
    if int(mol.data["nelec_A"]) != int(mol.data["nelec_B"]):
        raise ValueError("CASPT2 currently supports closed-shell singlets")

    nbf = int(mol.data.get_basis()["nbf"])
    ncore = int(settings.frozen_core)
    nact = int(settings.active_orbitals)
    nelec = (int(mol.data["nelec_A"]), int(mol.data["nelec_B"]))
    active_nelec = (nelec[0] - ncore, nelec[1] - ncore)
    if nact <= 0 or ncore + nact > nbf:
        raise ValueError("CASPT2 needs a valid [cas] active_orbitals / frozen_core")
    if min(active_nelec) < 0:
        raise ValueError("frozen_core exceeds the available electron count")

    roots = _reference_roots(options)
    weights = np.full(len(roots), 1.0 / len(roots))

    t0 = time.time()
    # Reference orbitals: optimize for a CASSCF reference (state-averaged when
    # several roots are requested), else use the supplied RHF orbitals.
    if options.reference == "casscf":
        _run_casscf_reference(mol, ref_energy, roots, weights)

    oqp.fci_ao_integrals(mol)
    hcore_ao = _unpack_lower_triangle(np.asarray(mol.data["OQP::Hcore"], dtype=float), nbf)
    # Honour [cas] orbital_source: the standalone CASCI driver routes through
    # load_cas_mo_coeff, and this path did not -- so a PT2 run asking for
    # orbital_source=json silently used the current RHF coefficients instead,
    # giving a different reference (and a different active-space energy) from
    # the one requested.
    from oqp.library.cas_orbitals import load_cas_mo_coeff
    _default = np.asarray(mol.data["OQP::VEC_MO_A"], dtype=float).reshape((nbf, nbf)).T
    coeff, _orb_source = load_cas_mo_coeff(mol.config, nbf, _default)
    coeff = np.asarray(coeff, dtype=float)
    eri_ao = np.asarray(mol.data["OQP::AO_ERI"], dtype=float).reshape(
        (nbf, nbf, nbf, nbf), order="F"
    )
    enuc = float(mol.mol_energy.nenergy)

    # MS-CASPT2 with the Fock H0 uses the MULTI-SET construction (per-state
    # orbitals = the OpenMolcas state-specific Fock): each root's first-order
    # CASPT2 is solved in semicanonical orbitals built from its OWN density, with
    # the full generalized-Fock-MATRIX H0.  This is what matches OpenMolcas MS-
    # CASPT2; the prior single-set (state-averaged orbitals) path over-correlated
    # strongly-mixed roots by up to ~12 mEh.  Multi-set does its OWN per-state
    # semicanonicalization, so it uses the reference orbitals (``coeff``) and the
    # full, pre-freeze integrals directly.  Dyall H0 and XMS keep the single-set
    # path below.  The QDPT family (MCQDPT2) is BY DEFINITION single-set: one
    # state-averaged canonical Fock (GAMESS IFORB=1), one common diag(eps) H0,
    # ket-state-specific E0 -- so it never takes the multi-set branch.
    if options.variant == "ms" and options.h0 == "fock" and options.family == "caspt2":
        nfrozen0 = _standard_core_count(mol) if options.frozen < 0 else int(options.frozen)
        nfrozen0 = max(0, min(nfrozen0, ncore))
        s2_ms, mult_ms = _multiset_spin_diagnostics(
            hcore_ao, eri_ao, coeff, ncore, nact, active_nelec, enuc, settings, roots)
        result = _multistate_multiset(mol, hcore_ao, eri_ao, coeff, settings, ncore,
                                      nact, active_nelec, nbf, enuc, roots, options, nfrozen0)
        _multistate_finish(mol, ref_energy, options, settings, ncore, nact, active_nelec,
                           roots, result, s2_ms, mult_ms, time.time() - t0)
        return mol.energies

    coeff_sc, h1e, eri, coeffs, energies, dets, eps, D_sa = _semicanonicalize(
        hcore_ao, eri_ao, coeff, ncore, nact, active_nelec, enuc, settings, roots, weights
    )
    s2, mult = fci_spin_diagnostics(coeffs, dets, nact, active_nelec)

    # commit semicanonical orbitals (the rotation leaves the reference energy invariant)
    target = np.asarray(mol.data["OQP::VEC_MO_A"], dtype=float)
    mol.data["OQP::VEC_MO_A"][...] = np.ascontiguousarray(coeff_sc.T.reshape(target.shape))

    # PT2 frozen core: fold the deepest atomic cores out of the first-order space.
    # Default (frozen<0) = the standard deep cores, matching OpenMolcas and removing the
    # spurious deep-core over-correlation; frozen=0 correlates all; frozen=N freezes N.
    nfrozen = _standard_core_count(mol) if options.frozen < 0 else int(options.frozen)
    nfrozen = max(0, min(nfrozen, ncore))
    options._pt2_nfrozen = nfrozen
    if nfrozen:
        h1e, eri, eps, D_sa, ncore, nbf, enuc = _freeze_core(
            h1e, eri, eps, D_sa, ncore, nbf, enuc, nfrozen)

    # QDPT family defaults to the matrix-free direct engine (diagonal H0):
    # no dense N x N Hamiltonian, no O(N^3) external diagonalization, only the
    # singles/doubles interacting space streamed from the reference support.
    # auto = the NumPy streaming path (measured fastest at scale); the liboqp
    # hash kernel is the explicit engine=fortran opt-in (see qdpt2_direct.py).
    use_direct = options.family == "qdpt" and options.engine in {"auto", "direct", "fortran"}

    if options.variant == "caspt2":
        if use_direct:
            from oqp.library.qdpt2_direct import direct_qdpt2
            res = direct_qdpt2(h1e, eri, coeffs, energies, dets, eps, D_sa, ncore,
                               nact, active_nelec, nbf, enuc, [options.root],
                               options, rotate=False)
            e_casci = float(res["ref_energies"][0])
            e2 = float(res["e2"][0])
            e_caspt2 = e_casci + e2
            mol.energies = [e_caspt2]
            mol.mol_energy.energy = e_caspt2
            mol.data["OQP::CASPT2_ENERGIES"] = np.ascontiguousarray([e_caspt2], dtype=np.float64)
            mol.data["OQP::CASPT2_REFERENCE_ENERGIES"] = np.ascontiguousarray([e_casci], dtype=np.float64)
            mol.data["OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS"] = np.ascontiguousarray([e2], dtype=np.float64)
            _write_log(mol, ref_energy, options, settings, ncore, nact, active_nelec,
                       e_casci, e2, e_caspt2, res["n_external"], e_casci,
                       float(res["min_denoms"][0]), float(s2[options.root]),
                       int(mult[options.root]), time.time() - t0)
        else:
            _single_state_finish(mol, ref_energy, options, settings, ncore, nact, active_nelec,
                                 nbf, enuc, h1e, eri, coeffs, energies, dets, eps, D_sa,
                                 s2, mult, time.time() - t0)
    else:
        if use_direct:
            from oqp.library.qdpt2_direct import direct_qdpt2
            result = direct_qdpt2(h1e, eri, coeffs, energies, dets, eps, D_sa, ncore,
                                  nact, active_nelec, nbf, enuc, roots, options,
                                  rotate=(options.variant == "xms"))
        else:
            result = _multistate(h1e, eri, coeffs, energies, dets, eps, D_sa, ncore, nact,
                                 active_nelec, nbf, enuc, roots, options)
        _multistate_finish(mol, ref_energy, options, settings, ncore, nact, active_nelec,
                           roots, result, s2, mult, time.time() - t0)
    return mol.energies


def _run_casscf_reference(mol, ref_energy, roots, weights):
    """Run the native (SA-)CASSCF to provide reference orbitals."""
    cfg = mol.config["input"]
    original_method = cfg.get("method")
    if len(roots) > 1:
        cfg["method"] = "sa-casscf"
        sa = mol.config.setdefault("state_average", {}) or {}
        mol.config["state_average"] = sa
        sa.setdefault("enabled", "true")
        sa["target_roots"] = ",".join(str(r) for r in roots)
        sa["nstate"] = str(len(roots))
        sa.setdefault("weights", ",".join(f"{w:.10f}" for w in weights))
    else:
        cfg["method"] = "casscf"
    try:
        CASSCF(mol).energy(ref_energy=ref_energy)
    finally:
        cfg["method"] = original_method


def _single_state_finish(mol, ref_energy, options, settings, ncore, nact, active_nelec,
                         norb, enuc, h1e, eri, coeffs, energies, dets, eps, D_sa, s2, mult, wall):
    root = options.root
    e_casci = float(energies[root])
    # active occupations for the IPEA shift = diagonal of the active 1-RDM block
    active_occ = np.diag(D_sa)[ncore:ncore + nact]

    if options.contraction == "strong":
        # RDM-based strongly contracted NEVPT2 (SC-NEVPT2); reproduces PySCF/ORCA
        # to <~1 uEh.  No external determinant space: the contracted perturbers
        # are formed directly from the active 1-/2-/3-/4-RDM and integral blocks.
        from oqp.library.nevpt2_sc import sc_nevpt2_energy
        # sc_nevpt2_energy takes no regularisation: the contracted denominators
        # are built inside it.  Accepting a shift and returning the unshifted
        # energy -- while the PT2 summary prints the shift back to the user --
        # is worse than refusing, so refuse.
        _unapplied = [name for name, value in (
            ("level_shift", options.level_shift),
            ("imaginary_shift", options.imaginary_shift),
            ("edshft", options.edshft)) if value]
        if _unapplied:
            raise ValueError(
                "[pt2] %s cannot be applied to strongly contracted NEVPT2 "
                "(h0=dyall, contraction=strong): the contracted denominators are "
                "formed internally and no shift reaches them. Remove the shift, "
                "or use contraction=none for a shifted NEVPT2."
                % ", ".join(_unapplied))
        e2, comp = sc_nevpt2_energy(h1e, eri, eps, ncore, nact, active_nelec,
                                    coeffs[:, root])
        e_caspt2 = e_casci + e2
        e_ref_check = e_casci
        min_denom = float("inf")
        n_external = sum(1 for _ in comp)
        _write_log(mol, ref_energy, options, settings, ncore, nact, active_nelec,
                   e_casci, e2, e_caspt2, n_external, e_ref_check, min_denom,
                   float(s2[root]), int(mult[root]), wall, sc_components=comp)
        mol.energies = [e_caspt2]
        mol.mol_energy.energy = e_caspt2
        mol.data["OQP::CASPT2_ENERGIES"] = np.ascontiguousarray([e_caspt2], dtype=np.float64)
        mol.data["OQP::CASPT2_REFERENCE_ENERGIES"] = np.ascontiguousarray([e_casci], dtype=np.float64)
        mol.data["OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS"] = np.ascontiguousarray([e2], dtype=np.float64)
        return

    full_dets, det_index, hfull, h0op, external = _build_operators(
        h1e, eri, eps, ncore, nact, active_nelec, norb, options.max_det, options.h0,
        active_occ=active_occ, ipea=options.ipea_shift
    )
    vec = _embed_reference(coeffs[:, root], dets, ncore, nact, norb, det_index)
    e_ref_check = float(vec @ hfull @ vec) + enuc
    _V, _C, e2, min_denom = _first_order(hfull, h0op, vec, external, options)
    e_caspt2 = e_casci + e2

    mol.energies = [e_caspt2]
    mol.mol_energy.energy = e_caspt2
    mol.data["OQP::CASPT2_ENERGIES"] = np.ascontiguousarray([e_caspt2], dtype=np.float64)
    mol.data["OQP::CASPT2_REFERENCE_ENERGIES"] = np.ascontiguousarray([e_casci], dtype=np.float64)
    mol.data["OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS"] = np.ascontiguousarray([e2], dtype=np.float64)

    _write_log(mol, ref_energy, options, settings, ncore, nact, active_nelec,
               e_casci, e2, e_caspt2, len(external), e_ref_check, min_denom,
               float(s2[root]), int(mult[root]), wall)


def _multistate_finish(mol, ref_energy, options, settings, ncore, nact, active_nelec,
                       roots, result, s2, mult, wall):
    ms = result["ms_energies"]
    mol.energies = [float(e) for e in ms]
    mol.mol_energy.energy = float(ms[0])
    mol.data["OQP::CASPT2_ENERGIES"] = np.ascontiguousarray(ms, dtype=np.float64)
    mol.data["OQP::CASPT2_REFERENCE_ENERGIES"] = np.ascontiguousarray(result["ref_energies"], dtype=np.float64)
    mol.data["OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS"] = np.ascontiguousarray(result["e2"], dtype=np.float64)
    mol.data["OQP::CASPT2_SS_ENERGIES"] = np.ascontiguousarray(result["ss_energies"], dtype=np.float64)
    mol.data["OQP::CASPT2_EFFECTIVE_HAMILTONIAN"] = np.ascontiguousarray(result["heff"], dtype=np.float64)
    mol.data["OQP::CASPT2_MIXING_VECTORS"] = np.ascontiguousarray(result["mixing"], dtype=np.float64)
    mol.data["OQP::CASPT2_TARGET_ROOTS"] = np.ascontiguousarray(roots, dtype=np.int64)

    _write_ms_log(mol, ref_energy, options, settings, ncore, nact, active_nelec,
                  roots, result, s2, mult, wall)


# --------------------------------------------------------------------------- log blocks
def _header(mol, options, ref_energy, settings, ncore, nact, active_nelec, ref_label, title):
    ne = active_nelec[0] + active_nelec[1]
    print_module_banner(mol, title.split(" (")[0], title)
    _log(mol)
    _log(mol, "   ==============================================")
    _log(mol, f"   PyOQP: {title}")
    _log(mol, "   ==============================================")
    _log(mol)
    _log(mol, f"   PyOQP method:                       {str(mol.config['input']['method']).lower()}")
    _log(mol, f"   PyOQP reference:                    {ref_label}")
    if ref_energy:
        _log(mol, f"   PyOQP RHF reference energy:         {ref_energy[0]:18.10f}")
    _log(mol, f"   PyOQP active space:                 CAS({ne},{nact})")
    _log(mol, f"   PyOQP inactive (core) orbitals:     {ncore}")
    _log(mol, f"   PyOQP PT2 frozen-core orbitals:     {int(getattr(options, '_pt2_nfrozen', 0))}")
    _log(mol, f"   PyOQP active orbitals:              {nact}")
    _log(mol, f"   PyOQP CI solver:                    {settings.solver}")
    if options.h0 == "dyall":
        h0_label = "Dyall (NEVPT2)"
    elif options.family == "qdpt":
        h0_label = "diagonal Fock (MP-like, QDPT / GAMESS convention)"
    else:
        h0_label = "Fock 1e (CASPT2)"
    _log(mol, f"   PyOQP zeroth-order H0:              {h0_label}")
    _log(mol, f"   PyOQP semicanonical orbitals:       {'yes' if options.semi_canonical else 'no'}")
    if options.contraction == "strong":
        _log(mol, "   PyOQP IPEA shift:                   not applicable (strong contraction)")
    else:
        ipea_note = "  (1e diagonal H0 shift; not bit-identical to Molcas' contracted IPEA)" if options.ipea_shift else ""
        _log(mol, f"   PyOQP IPEA shift:                   {options.ipea_shift:.4f}{ipea_note}")
    _log(mol, f"   PyOQP real level shift:             {options.level_shift:.4f}")
    _log(mol, f"   PyOQP imaginary shift:              {options.imaginary_shift:.4f}")
    _log(mol, f"   PyOQP ISA shift (edshft):           {options.edshft:.4f}")


def _write_log(mol, ref_energy, options, settings, ncore, nact, active_nelec,
               e_casci, e2, e_caspt2, n_external, e_ref_check, min_denom, s2, mult, wall,
               sc_components=None):
    ref_label = "CASSCF" if options.reference == "casscf" else "CASCI"
    if options.contraction == "strong":
        method_name, title = "SC-NEVPT2", "SC-NEVPT2 (strongly contracted, Dyall zeroth order)"
    elif options.h0 == "dyall":
        method_name, title = "NEVPT2", "NEVPT2 (uncontracted, Dyall zeroth order)"
    elif options.family == "qdpt":
        method_name, title = "MRMP2", "MRMP2 (determinant MRPT2, diagonal-Fock zeroth order)"
    else:
        method_name, title = "CASPT2", "CASPT2 (uncontracted, Fock zeroth order)"
    _header(mol, options, ref_energy, settings, ncore, nact, active_nelec, ref_label, title)
    _log(mol, f"   PyOQP reference root:               {options.root}")
    _log(mol)
    _log(mol, "   --- second order ---")
    if sc_components is not None:
        _log(mol, f"   PyOQP contraction:                  strong (RDM-based SC-NEVPT2)")
        _log(mol, "   PyOQP Dyall subspace energies:")
        for key in ("Sr", "Si", "Sijrs", "Sijr", "Srsi", "Srs", "Sij", "Sir"):
            _log(mol, f"      {key:7s} {sc_components[key]:18.10f}")
    else:
        _log(mol, f"   PyOQP external (first-order) space:  {n_external}")
        _log(mol, f"   PyOQP smallest denominator:         {min_denom:12.4e}{_intruder_note(min_denom, options)}")
    _log(mol, f"   PyOQP reference <S^2> / mult:        {s2:8.6f} / {mult}")
    ref_drift = abs(e_ref_check - e_casci)
    flag = "" if ref_drift < 1.0e-7 else "   <-- WARNING: reference self-check drift"
    _log(mol, f"   PyOQP reference self-check:         {e_ref_check:18.10f}{flag}")
    _log(mol)
    _log(mol, f"   PyOQP E(0) reference  ({ref_label:6s}):     {e_casci:18.10f}")
    _log(mol, f"   PyOQP E(2) correlation:             {e2:18.10f}")
    _log(mol, f"   PyOQP E({method_name:<8}) total:        {e_caspt2:18.10f}")
    _log(mol)
    _log(mol, f"   PyOQP timing:                       {wall:.2f} s")
    _log(mol)


def _write_ms_log(mol, ref_energy, options, settings, ncore, nact, active_nelec,
                  roots, result, s2, mult, wall):
    ref_label = "CASSCF" if options.reference == "casscf" else "CASCI"
    if options.family == "qdpt":
        label = "XMCQDPT2" if result["rotated"] else "MCQDPT2"
    else:
        label = "XMS-CASPT2" if result["rotated"] else "MS-CASPT2"
    h0_word = "Dyall" if options.h0 == "dyall" else "Fock"
    construction = "multi-set" if result.get("multiset") else "single-set"
    _header(mol, options, ref_energy, settings, ncore, nact, active_nelec, ref_label,
            f"{label} ({construction}, {h0_word} zeroth order)")
    _log(mol, f"   PyOQP reference roots:              {', '.join(str(r) for r in roots)}")
    _log(mol, f"   PyOQP XMS reference rotation:       {'yes' if result['rotated'] else 'no'}")
    _log(mol, f"   PyOQP external (first-order) space:  {result['n_external']}")
    flag = "" if result["ref_drift"] < 1.0e-7 else "   <-- WARNING: reference self-check drift"
    _log(mol, f"   PyOQP reference self-check drift:   {result['ref_drift']:.2e}{flag}")
    _log(mol)
    _log(mol, "   --- per-state second order (effective-Hamiltonian diagonal) ---")
    _log(mol, f"   {'root':>5} {'E(ref)':>18} {'E(2)':>16} {'E(SS-CASPT2)':>18} {'min|denom|':>12}")
    min_denoms = result.get("min_denoms", np.full(len(roots), float("inf")))
    for i, r in enumerate(roots):
        note = _intruder_note(float(min_denoms[i]), options)
        _log(mol, f"   {r:5d} {result['ref_energies'][i]:18.10f} {result['e2'][i]:16.10f} "
                  f"{result['ss_energies'][i]:18.10f} {float(min_denoms[i]):12.4e}{note}")
    _log(mol)
    _log(mol, f"   --- {label} (effective Hamiltonian diagonalized) ---")
    _log(mol, f"   {'state':>5} {'E('+label+')':>18}    {'shift vs SS-CASPT2':>20}")
    ss_sorted = np.sort(result["ss_energies"])
    for i, e in enumerate(result["ms_energies"]):
        shift = float(e - ss_sorted[i])
        _log(mol, f"   {i:5d} {float(e):18.10f}    {shift:20.2e}")
    _log(mol)
    _log(mol, f"   PyOQP timing:                       {wall:.2f} s")
    _log(mol)
