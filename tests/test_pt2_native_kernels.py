"""Equivalence pin for the Fortran PT2 engines (pt2_kernel.F90, nevpt2_koopmans.F90).

These kernels replace arithmetic that used to run in Python/NumPy on the PT2
execution path:

* ``pt2_kernel.F90``       -- the mean-field Fock build, the Dyall active
                              dressing, and the determinant-space bookkeeping
                              (external space, diagonal H0, occupation blocks)
                              of ``caspt2_dyall.py``;
* ``nevpt2_koopmans.F90``  -- the two eri-folded 4-pdm intermediates and the
                              ``a16`` Koopmans intermediate of ``nevpt2_sc.py``.

Each Python implementation stays in the tree precisely so it can serve as this
pin, and as the fallback when liboqp predates the symbol.  Every test below
evaluates the same function twice -- once through liboqp, once with the backend
probe forced to ``None`` -- so the two paths are compared directly rather than
against a hard-coded constant.

The integer-valued results (external indices, occupation blocks) must agree
EXACTLY: the block partition feeds a per-block eigendecomposition, so a
different intra-block member order would perturb the printed energies even
though it is mathematically equivalent.  The floating-point results are held to
round-off rather than to a loose tolerance that would hide a real discrepancy.
"""
import numpy as np
import pytest

import oqp.library.caspt2_dyall as cp
import oqp.library.nevpt2_sc as sc



# ---- Required Fortran symbol checks (Definition of Done Rule 1) ----
# Every new bind(C) procedure must be covered by a test that FAILS (not skips)
# when the symbol is absent.  These tests are the REQUIRED gate: they exist
# only to fail when a needed export is missing.
_REQUIRED_NEVPT2 = [
    'nevpt2_f3ca_f3ac', 'nevpt2_a16', 'nevpt2_a22',
    'nevpt2_hdm1', 'nevpt2_hdm2', 'nevpt2_hdm3',
    'nevpt2_a9', 'nevpt2_a7', 'nevpt2_a12', 'nevpt2_a13',
    'nevpt2_a17', 'nevpt2_a19', 'nevpt2_a23', 'nevpt2_a25',
    'nevpt2_a3', 'nevpt2_k27',
    'nevpt2_sijrs', 'nevpt2_srs', 'nevpt2_sij',
    'nevpt2_sijr', 'nevpt2_srsi', 'nevpt2_sir',
    'nevpt2_sr', 'nevpt2_si',
    'nevpt2_fold_h1e',
    'nevpt2_make_rdms', 'nevpt2_make_f3',
]
_REQUIRED_PT2 = [
    'pt2_effective_fock', 'pt2_h0_dyall_active',
    'pt2_external_indices', 'pt2_diagonal_h0',
    'pt2_occupation_blocks', 'pt2_minor_transform',
]


def _has(name, probe):
    backend = probe()
    return backend is not None and hasattr(backend[0], name)


def test_nevpt2_koopmans_symbols_required():
    """Fail when any NEVPT2 Koopmans Fortran symbol is absent (DoD Rule 1)."""
    nloaded = sum(1 for sym in _REQUIRED_NEVPT2 if _has(sym, sc._koopmans_lib))
    assert nloaded == len(_REQUIRED_NEVPT2), (
        f"Only {nloaded}/{len(_REQUIRED_NEVPT2)} NEVPT2 Koopmans symbols "
        f"loaded; need all for Fortran engine to work.  Rebuild liboqp."
    )


def test_pt2_kernel_symbols_required():
    """Fail when any PT2 kernel Fortran symbol is absent (DoD Rule 1)."""
    nloaded = 0
    for sym in _REQUIRED_PT2:
        if _has(sym, cp._pt2_lib):
            nloaded += 1
    assert nloaded == len(_REQUIRED_PT2), (
        f"Only {nloaded}/{len(_REQUIRED_PT2)} PT2 kernel symbols "
        f"loaded; need all for Fortran engine to work.  Rebuild liboqp."
    )


def _without_backend(mod, probe_name, fn):
    """Evaluate ``fn`` with the module's native-backend probe disabled."""
    saved = getattr(mod, probe_name)
    setattr(mod, probe_name, lambda: None)
    try:
        return fn()
    finally:
        setattr(mod, probe_name, saved)


def _both(mod, probe_name, fn):
    return fn(), _without_backend(mod, probe_name, fn)


def _symmetric_eri(nbf, rng):
    """Random ERIs carrying the full 8-fold permutational symmetry."""
    a = rng.standard_normal((nbf, nbf, nbf, nbf)) * 0.1
    a = a + a.transpose(1, 0, 2, 3)
    a = a + a.transpose(0, 1, 3, 2)
    return a + a.transpose(2, 3, 0, 1)


# --------------------------------------------------------------------- pt2_kernel
pt2_only = pytest.mark.skipif(
    not _has("pt2_effective_fock", cp._pt2_lib),
    reason="liboqp pt2_kernel engine unavailable",
)

_skip_if_no_koopmans = pytest.mark.skipif(
    not _has("nevpt2_fold_h1e", sc._koopmans_lib),
    reason="NEVPT2 Koopmans Fortran kernels not loaded",
)


@pytest.mark.parametrize("nbf,ncore,nact", [(8, 2, 4), (12, 3, 6)])
def test_effective_fock_matches_python(nbf, ncore, nact):
    rng = np.random.default_rng(11 + nbf)
    h1e = rng.standard_normal((nbf, nbf))
    h1e = h1e + h1e.T
    eri = _symmetric_eri(nbf, rng)
    D = rng.standard_normal((nbf, nbf))
    D = D + D.T

    native, python = _both(cp, "_pt2_lib", lambda: cp._effective_fock(h1e, eri, D))
    assert native.shape == python.shape == (nbf, nbf)
    assert np.max(np.abs(native - python)) < 1e-11


@pytest.mark.parametrize("nbf,ncore,nact", [(8, 2, 4), (12, 3, 6)])
def test_h0_integrals_dyall_matches_python(nbf, ncore, nact):
    rng = np.random.default_rng(23 + nbf)
    h1e = rng.standard_normal((nbf, nbf))
    h1e = h1e + h1e.T
    eri = _symmetric_eri(nbf, rng)
    eps = rng.standard_normal(nbf)

    native, python = _both(
        cp, "_pt2_lib",
        lambda: cp._h0_integrals(h1e, eri, eps, ncore, nact, "dyall")[0])
    assert np.max(np.abs(native - python)) < 1e-11
    # the inactive/virtual diagonal is untouched by the dressing
    assert np.allclose(np.diag(native)[ncore + nact:], eps[ncore + nact:])


@pytest.mark.parametrize("nbf,ncore,nact", [(8, 2, 4), (10, 2, 4)])
def test_external_indices_match_python_exactly(nbf, ncore, nact):
    from oqp.library.fci import _determinants

    dets = _determinants(nbf, (ncore + 2, ncore + 2))
    native, python = _both(
        cp, "_pt2_lib", lambda: cp._external_indices(dets, ncore, nact, nbf))
    assert np.array_equal(np.asarray(native), np.asarray(python))
    # the CAS configurations are exactly the complement
    assert 0 < len(native) < len(dets)


@pytest.mark.parametrize("nbf,ncore,nact", [(8, 2, 4), (10, 2, 4)])
def test_diagonal_zeroth_order_matches_python(nbf, ncore, nact):
    from oqp.library.fci import _determinants

    rng = np.random.default_rng(31 + nbf)
    dets = _determinants(nbf, (ncore + 2, ncore + 2))
    eps = rng.standard_normal(nbf)
    h0_1e = np.diag(eps)
    g0 = np.zeros((nbf,) * 4)

    native, python = _both(
        cp, "_pt2_lib",
        lambda: cp._diagonal_zeroth_order(dets, h0_1e, g0, nbf))
    assert native is not None and python is not None
    assert np.max(np.abs(native - python)) < 1e-11


# (12, 3, 4) is 627k determinants: it is the size at which the native sort was
# actually measured, and with 24 signature bits it exercises three radix passes
# (an odd count, so the ping-pong ends in the scratch buffers and the copy-back
# path runs) over an array where nearly every signature is a tie.
@pytest.mark.parametrize("nbf,ncore,nact", [(8, 2, 4), (10, 2, 4), (12, 3, 4)])
def test_occupation_blocks_match_python_exactly(nbf, ncore, nact):
    from oqp.library.fci import _determinants

    dets = _determinants(nbf, (ncore + 2, ncore + 2))
    ext = cp._external_indices(dets, ncore, nact, nbf)

    native, python = _both(
        cp, "_pt2_lib",
        lambda: cp._occupation_blocks(dets, ext, ncore, nact, nbf))
    assert len(native) == len(python)
    for a, b in zip(native, python):
        # exact, including the member order inside each block
        assert np.array_equal(np.asarray(a), np.asarray(b))
    # the blocks partition the external space
    assert sum(len(g) for g in native) == len(ext)
    assert np.array_equal(np.sort(np.concatenate(native)), np.arange(len(ext)))

    # ... and pin the stability property directly rather than only through the
    # fallback: reading the blocks out in order must reproduce numpy's stable
    # argsort of the signature, ties included.
    act_mask = ((1 << nact) - 1) << ncore
    keep = np.int64(((1 << nbf) - 1) & ~act_mask) | (
        np.int64(((1 << nbf) - 1) & ~act_mask) << nbf)
    sig = np.asarray(dets, dtype=np.int64)[np.asarray(ext)] & keep
    assert np.array_equal(np.concatenate(native),
                          np.argsort(sig, kind="stable"))
    assert len(native) > 1          # the partition is not trivial


@pytest.mark.parametrize("norb,nelec", [(4, 2), (6, 3), (8, 4)])
def test_minor_transform_matches_python(norb, nelec):
    rng = np.random.default_rng(101 + norb)
    # a genuine orbital rotation: the minors are then O(1) and the pivots never
    # reach dgetf2's safe-minimum branch, which is the regime that matters
    R = np.linalg.qr(rng.standard_normal((norb, norb)))[0]
    occ = cp._occ_tuples(norb, nelec)

    native, python = _both(cp, "_pt2_lib", lambda: cp._minor_transform(R, occ))
    assert native.shape == python.shape == (len(occ), len(occ))
    assert np.max(np.abs(native - python)) < 1e-11
    # the string transform of a unitary R is itself unitary; a wrong sign on
    # any single minor breaks this while staying inside a loose tolerance
    assert np.max(np.abs(native @ native.T - np.eye(len(occ)))) < 1e-10


def test_minor_transform_handles_a_singular_rotation():
    """An exactly zero column of R drives dgetf2's zero-pivot early exit, which
    a well-conditioned rotation never reaches.  Every minor selecting that
    column must come back as an EXACT zero, not merely a small number: the
    elimination leaves the pivot bit-zero, so an implementation that divided by
    it (or that accumulated the diagonal product past the exit) would show up
    here as a NaN or a denormal rather than as a tolerance failure."""
    norb, nelec = 5, 2
    R = np.linalg.qr(np.random.default_rng(7).standard_normal((norb, norb)))[0]
    R[:, 3] = 0.0
    occ = cp._occ_tuples(norb, nelec)
    native, python = _both(cp, "_pt2_lib", lambda: cp._minor_transform(R, occ))
    assert np.all(np.isfinite(native))
    assert np.max(np.abs(native - python)) < 1e-11
    # The singular minors are exact on both sides (the well-conditioned ones
    # only agree to an ulp, so this is asserted on the singular columns only).
    singular = [s for s, cols in enumerate(occ) if 3 in cols]
    assert singular
    assert np.all(native[:, singular] == 0.0)
    assert np.all(python[:, singular] == 0.0)


# ---------------------------------------------------------------- Koopmans
koopmans_only = pytest.mark.skipif(
    not _has("nevpt2_f3ca_f3ac", sc._koopmans_lib),
    reason="liboqp nevpt2_koopmans engine unavailable",
)


def _koopmans_operands(n, seed):
    rng = np.random.default_rng(seed)
    h1e = rng.standard_normal((n, n))
    h1e = 0.5 * (h1e + h1e.T)
    h2e = rng.standard_normal((n, n, n, n)) * 0.1
    h2e = h2e + h2e.transpose(3, 2, 1, 0)      # real physicist-ordered ERIs
    dm1 = rng.standard_normal((n, n)) * 0.1
    dm2 = rng.standard_normal((n,) * 4) * 0.1
    dm3 = rng.standard_normal((n,) * 6) * 0.1
    dm4 = rng.standard_normal((n,) * 8) * 0.1
    return h1e, h2e, dm1, dm2, dm3, dm4


@pytest.mark.parametrize("n", [3, 4, 5])
def test_f3ca_f3ac_match_python(n):
    _h1e, h2e, _d1, _d2, _dm3, dm4 = _koopmans_operands(n, 41 + n)
    native, python = _both(sc, "_koopmans_lib", lambda: sc._f3ca_f3ac(h2e, dm4))
    for a, b in zip(native, python):
        assert a.shape == (n,) * 6
        assert np.max(np.abs(a - b)) < 1e-11


@pytest.mark.parametrize("n", [3, 4, 5])
def test_a16_matches_python(n):
    h1e, h2e, _dm1, _dm2, dm3, dm4 = _koopmans_operands(n, 53 + n)
    f3ca, f3ac = sc._f3ca_f3ac(h2e, dm4)
    native, python = _both(
        sc, "_koopmans_lib", lambda: sc._a16(h1e, h2e, dm3, f3ca, f3ac))
    assert native.shape == (n,) * 6
    assert np.max(np.abs(native - python)) < 1e-11


@pytest.mark.parametrize("n", [3, 4, 5])
def test_a22_matches_python(n):
    h1e, h2e, _dm1, dm2, dm3, dm4 = _koopmans_operands(n, 67 + n)
    f3ca, f3ac = sc._f3ca_f3ac(h2e, dm4)
    native, python = _both(
        sc, "_koopmans_lib", lambda: sc._a22(h1e, h2e, dm2, dm3, f3ca, f3ac))
    assert native.shape == (n,) * 6
    assert np.max(np.abs(native - python)) < 1e-11


@pytest.mark.parametrize("n", [3, 4, 5])
def test_hdm3_matches_python(n):
    _h1e, _h2e, dm1, dm2, dm3, _dm4 = _koopmans_operands(n, 79 + n)
    hdm1 = sc._hdm1(dm1)
    hdm2 = sc._hdm2(dm1, dm2)
    native, python = _both(
        sc, "_koopmans_lib", lambda: sc._hdm3(dm1, dm2, dm3, hdm1, hdm2))
    assert native.shape == (n,) * 6
    assert np.max(np.abs(native - python)) < 1e-11


@pytest.mark.parametrize("n", [3, 4, 5])
def test_a9_matches_python(n):
    h1e, h2e, dm1, dm2, dm3, _dm4 = _koopmans_operands(n, 83 + n)
    hdm1 = sc._hdm1(dm1)
    hdm2 = sc._hdm2(dm1, dm2)
    hdm3 = sc._hdm3(dm1, dm2, dm3, hdm1, hdm2)
    native, python = _both(
        sc, "_koopmans_lib", lambda: sc._a9(h1e, h2e, hdm1, hdm2, hdm3))
    assert native.shape == (n,) * 4
    assert np.max(np.abs(native - python)) < 1e-11


@pytest.mark.parametrize("n", [3, 4, 5])
def test_a7_matches_python(n):
    h1e, h2e, dm1, dm2, dm3, _dm4 = _koopmans_operands(n, 89 + n)
    native, python = _both(
        sc, "_koopmans_lib", lambda: sc._a7(h1e, h2e, dm1, dm2, dm3))
    # _a7 returns (rm2, a7); both are part of the contract -- rm2 goes straight
    # into the Srs norm, so a drift there would move the energy without ever
    # touching a7.
    for a, b in zip(native, python):
        assert a.shape == (n,) * 4
        assert np.max(np.abs(a - b)) < 1e-11


@pytest.mark.parametrize("n", [3, 4, 5])
def test_a12_matches_python(n):
    h1e, h2e, dm1, dm2, dm3, _dm4 = _koopmans_operands(n, 97 + n)
    native, python = _both(
        sc, "_koopmans_lib", lambda: sc._a12(h1e, h2e, dm1, dm2, dm3))
    assert native.shape == (n,) * 4
    assert np.max(np.abs(native - python)) < 1e-11


@pytest.mark.parametrize("n", [3, 4, 5])
def test_a13_matches_python(n):
    h1e, h2e, dm1, dm2, dm3, _dm4 = _koopmans_operands(n, 103 + n)
    native, python = _both(
        sc, "_koopmans_lib", lambda: sc._a13(h1e, h2e, dm1, dm2, dm3))
    assert native.shape == (n,) * 4
    assert np.max(np.abs(native - python)) < 1e-11

# ---- Additional Koopmans intermediate kernel tests (post-PR470) ----

@_skip_if_no_koopmans
@pytest.mark.parametrize("n", [3, 4, 5])
def test_hdm1_matches_python(n):
    _h1e, _h2e, dm1, _dm2, _dm3, _dm4 = _koopmans_operands(n, 131 + n)
    native, python = _both(sc, "_koopmans_lib", lambda: sc._hdm1(dm1))
    assert native.shape == (n, n)
    assert np.max(np.abs(native - python)) < 1e-11


@_skip_if_no_koopmans
@pytest.mark.parametrize("n", [3, 4, 5])
def test_hdm2_matches_python(n):
    _h1e, _h2e, dm1, dm2, _dm3, _dm4 = _koopmans_operands(n, 139 + n)
    native, python = _both(sc, "_koopmans_lib", lambda: sc._hdm2(dm1, dm2))
    assert native.shape == (n,) * 4
    assert np.max(np.abs(native - python)) < 1e-11


@_skip_if_no_koopmans
@pytest.mark.parametrize("n", [3, 4, 5])
def test_a17_matches_python(n):
    h1e, h2e, _dm1, dm2, dm3, _dm4 = _koopmans_operands(n, 149 + n)
    native, python = _both(sc, "_koopmans_lib", lambda: sc._a17(h1e, h2e, dm2, dm3))
    assert native.shape == (n,) * 4
    assert np.max(np.abs(native - python)) < 1e-11


@_skip_if_no_koopmans
@pytest.mark.parametrize("n", [3, 4, 5])
def test_a19_matches_python(n):
    h1e, h2e, dm1, dm2, _dm3, _dm4 = _koopmans_operands(n, 157 + n)
    native, python = _both(sc, "_koopmans_lib", lambda: sc._a19(h1e, h2e, dm1, dm2))
    assert native.shape == (n, n)
    assert np.max(np.abs(native - python)) < 1e-11


@_skip_if_no_koopmans
@pytest.mark.parametrize("n", [3, 4, 5])
def test_a23_matches_python(n):
    h1e, h2e, dm1, dm2, dm3, _dm4 = _koopmans_operands(n, 163 + n)
    native, python = _both(sc, "_koopmans_lib", lambda: sc._a23(h1e, h2e, dm1, dm2, dm3))
    assert native.shape == (n,) * 4
    assert np.max(np.abs(native - python)) < 1e-11


@_skip_if_no_koopmans
@pytest.mark.parametrize("n", [3, 4, 5])
def test_a25_matches_python(n):
    h1e, h2e, dm1, dm2, _dm3, _dm4 = _koopmans_operands(n, 173 + n)
    native, python = _both(sc, "_koopmans_lib", lambda: sc._a25(h1e, h2e, dm1, dm2))
    assert native.shape == (n, n)
    assert np.max(np.abs(native - python)) < 1e-11


@_skip_if_no_koopmans
@pytest.mark.parametrize("n", [3, 4, 5])
def test_a3_matches_python(n):
    h1e, h2e, dm1, dm2, _dm3, _dm4 = _koopmans_operands(n, 179 + n)
    hdm1 = sc._hdm1(dm1)
    native, python = _both(sc, "_koopmans_lib", lambda: sc._a3(h1e, h2e, dm1, dm2, hdm1))
    assert native.shape == (n, n)
    assert np.max(np.abs(native - python)) < 1e-11


@_skip_if_no_koopmans
@pytest.mark.parametrize("n", [3, 4, 5])
def test_k27_matches_python(n):
    h1e, h2e, dm1, dm2, _dm3, _dm4 = _koopmans_operands(n, 191 + n)
    native, python = _both(sc, "_koopmans_lib", lambda: sc._k27(h1e, h2e, dm1, dm2))
    assert native.shape == (n, n)
    assert np.max(np.abs(native - python)) < 1e-11


@_skip_if_no_koopmans
@pytest.mark.parametrize("n", [3, 4, 5])
def test_fold_h1e_matches_python(n):
    rng = np.random.default_rng(197 + n)
    h1e = rng.standard_normal((n, n))
    h2e = rng.standard_normal((n, n, n, n)) * 0.1
    native, python = _both(sc, "_koopmans_lib", lambda: sc._fold_h1e(h1e, h2e))
    assert native.shape == (n, n)
    assert np.max(np.abs(native - python)) < 1e-11
