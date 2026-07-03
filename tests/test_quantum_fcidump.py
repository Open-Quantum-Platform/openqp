"""Tests for the oqp.quantum second-quantized Hamiltonian / FCIDUMP bridge.

These exercise the pure-Python integral transforms and FCIDUMP I/O, which do
not require a compiled ``liboqp``. When the compiled extension is unavailable
(e.g. a source checkout without a build), the modules are loaded directly from
their files so the math is still covered.
"""

import importlib.util
import os
import sys
import types

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_QDIR = os.path.join(_HERE, os.pardir, "pyoqp", "oqp", "quantum")


def _load(modname, filename):
    try:
        return __import__(f"oqp.quantum.{modname}", fromlist=[modname])
    except Exception:
        spec = importlib.util.spec_from_file_location(
            f"_oqp_quantum_{modname}", os.path.join(_QDIR, filename))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod


integrals = _load("integrals", "integrals.py")
fcidump = _load("fcidump", "fcidump.py")


# --------------------------------------------------------------------------
# Integral transforms
# --------------------------------------------------------------------------

def test_unpack_triangular_roundtrip():
    n = 4
    full = np.arange(n * n, dtype=float).reshape(n, n)
    sym = full + full.T
    packed = sym[np.tril_indices(n)]
    out = integrals.unpack_triangular(packed, n)
    assert np.allclose(out, sym)
    assert np.allclose(out, out.T)


def test_unpack_triangular_bad_size():
    with pytest.raises(ValueError):
        integrals.unpack_triangular(np.zeros(5), 4)


def test_ao_to_mo_1body_identity():
    h = np.array([[1.0, 0.2], [0.2, -0.5]])
    c = np.eye(2)
    assert np.allclose(integrals.ao_to_mo_1body(h, c), h)


def test_ao_to_mo_1body_matches_explicit():
    rng = np.random.default_rng(0)
    h = rng.standard_normal((3, 3))
    h = h + h.T
    c = rng.standard_normal((3, 3))
    ref = c.T @ h @ c
    assert np.allclose(integrals.ao_to_mo_1body(h, c), ref)


def test_ao_to_mo_2body_matches_full_einsum():
    rng = np.random.default_rng(1)
    n = 3
    eri = rng.standard_normal((n, n, n, n))
    # symmetrize to a physical chemist-notation tensor
    eri = eri + eri.transpose(1, 0, 2, 3)
    eri = eri + eri.transpose(0, 1, 3, 2)
    eri = eri + eri.transpose(2, 3, 0, 1)
    c = rng.standard_normal((n, n))  # (nao, nmo): AO rows, MO columns
    # Contract the AO (first) axis of c on every AO axis of the integral, the
    # same [ao, mo] convention used by ao_to_mo_1body.
    ref = np.einsum('ip,jq,kr,ls,ijkl->pqrs', c, c, c, c, eri, optimize=True)
    out = integrals.ao_to_mo_2body(eri, c)
    assert np.allclose(out, ref)


def test_ao_to_mo_1body_2body_share_coeff_convention():
    # Both transforms must consume the SAME (nao, nmo) coefficient matrix.
    # Build a one-electron-like AO tensor as an outer product so the two-body
    # transform of it equals the outer product of the one-body transforms.
    rng = np.random.default_rng(7)
    n = 4
    a = rng.standard_normal((n, n))
    s = a + a.T
    eri = np.einsum('pq,rs->pqrs', s, s)
    c = rng.standard_normal((n, n))
    s_mo = integrals.ao_to_mo_1body(s, c)
    g_mo = integrals.ao_to_mo_2body(eri, c)
    assert np.allclose(g_mo, np.einsum('pq,rs->pqrs', s_mo, s_mo))


def test_ao_to_mo_2body_preserves_symmetry():
    rng = np.random.default_rng(2)
    n = 4
    a = rng.standard_normal((n, n))
    s = a + a.T
    eri = np.einsum('pq,rs->pqrs', s, s)  # valid (pq|rs) symmetry pattern
    c = rng.standard_normal((n, n))
    mo = integrals.ao_to_mo_2body(eri, c)
    assert np.allclose(mo, mo.transpose(1, 0, 2, 3))
    assert np.allclose(mo, mo.transpose(0, 1, 3, 2))
    assert np.allclose(mo, mo.transpose(2, 3, 0, 1))


# --------------------------------------------------------------------------
# FCIDUMP round trip
# --------------------------------------------------------------------------

def _random_hamiltonian(n, seed):
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((n, n))
    h1 = a + a.T
    g = rng.standard_normal((n, n, n, n))
    # impose full 8-fold permutation symmetry of (pq|rs)
    g = g + g.transpose(1, 0, 2, 3)
    g = g + g.transpose(0, 1, 3, 2)
    g = g + g.transpose(2, 3, 0, 1)
    return h1, g


def test_fcidump_write_read_roundtrip(tmp_path):
    n = 4
    h1, h2 = _random_hamiltonian(n, 3)
    ecore = -1.2345
    path = str(tmp_path / "test.FCIDUMP")
    fcidump.write_fcidump(path, h1, h2, ecore, n_electrons=4, ms2=0)

    data = fcidump.read_fcidump(path)
    assert data["norb"] == n
    assert data["nelec"] == 4
    assert data["ms2"] == 0
    assert np.isclose(data["ecore"], ecore)
    assert np.allclose(data["h1"], h1)
    assert np.allclose(data["h2"], h2)


def test_fcidump_header_format(tmp_path):
    n = 2
    h1, h2 = _random_hamiltonian(n, 5)
    path = str(tmp_path / "h.FCIDUMP")
    fcidump.write_fcidump(path, h1, h2, 0.7, n_electrons=2, ms2=0,
                          orbsym=[1, 1])
    with open(path) as fh:
        head = fh.read(200)
    assert "&FCI" in head
    assert "NORB=2" in head
    assert "NELEC=2" in head
    assert "ORBSYM=" in head
    assert "&END" in head


def test_fcidump_two_electron_count(tmp_path):
    # For norb=2 with all integrals non-zero, the 8-fold-unique count of
    # (pq|rs) pairs is the number of unique (P>=Q) compound-index pairs:
    # npair = 3 -> npair*(npair+1)/2 = 6 unique two-electron entries.
    n = 2
    h1, h2 = _random_hamiltonian(n, 9)
    path = str(tmp_path / "c.FCIDUMP")
    fcidump.write_fcidump(path, h1, h2, 0.0, n_electrons=2, tol=-1.0)
    twoe = 0
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) == 5 and int(p[3]) != 0 and int(p[4]) != 0:
                twoe += 1
    assert twoe == 6


def test_fcidump_tolerance_drops_small(tmp_path):
    n = 2
    h1 = np.array([[1.0, 1e-15], [1e-15, 2.0]])
    h2 = np.zeros((n, n, n, n))
    path = str(tmp_path / "t.FCIDUMP")
    fcidump.write_fcidump(path, h1, h2, 0.0, n_electrons=2)
    data = fcidump.read_fcidump(path)
    # tiny off-diagonal dropped -> reads back as exact zero
    assert data["h1"][0, 1] == 0.0
    assert np.isclose(data["h1"][0, 0], 1.0)
    assert np.isclose(data["h1"][1, 1], 2.0)


# --------------------------------------------------------------------------
# OQP::ERI_AO storage convention (mirrors source/modules/int2e.F90)
# --------------------------------------------------------------------------

def _scatter_like_int2e(unique, nbf):
    """Reproduce the Fortran dump consumer: scatter the 8-fold-unique integral
    list into a dense tensor, written in Fortran column-major order, then
    return the flat buffer exactly as Python would read OQP::ERI_AO."""
    eri_f = np.zeros((nbf, nbf, nbf, nbf), order="F")
    for (i, j, k, l), v in unique.items():
        for a, b, c, d in (
            (i, j, k, l), (j, i, k, l), (i, j, l, k), (j, i, l, k),
            (k, l, i, j), (l, k, i, j), (k, l, j, i), (l, k, j, i),
        ):
            eri_f[a, b, c, d] = v
    # Fortran writes column-major; Python reads the raw bytes as a flat array.
    return eri_f.ravel(order="F")


def test_eri_ao_flat_reshape_is_chemist():
    # Build a reference chemist-notation tensor with full 8-fold symmetry.
    rng = np.random.default_rng(11)
    nbf = 3
    a = rng.standard_normal((nbf, nbf))
    s = a + a.T
    ref = np.einsum('pq,rs->pqrs', s, s)  # (pq|rs) with proper symmetry

    # Collect its 8-fold-unique entries the way int2_run emits them.
    unique = {}
    for i in range(nbf):
        for j in range(i + 1):
            for k in range(nbf):
                for l in range(k + 1):
                    if (i * (i + 1) // 2 + j) < (k * (k + 1) // 2 + l):
                        continue
                    unique[(i, j, k, l)] = ref[i, j, k, l]

    flat = _scatter_like_int2e(unique, nbf)
    # The Python side does flat.reshape(nbf,nbf,nbf,nbf) in C order.
    recovered = flat.reshape(nbf, nbf, nbf, nbf)
    # Thanks to 8-fold symmetry this equals chemist (pq|rs) directly.
    assert np.allclose(recovered, ref)


# --------------------------------------------------------------------------
# ints_2e cache behavior
# --------------------------------------------------------------------------

def test_eri_ao_recomputes_by_default(tmp_path, monkeypatch):
    calls = {"n": 0}

    fake_oqp = types.ModuleType("oqp")

    class FakeData(dict):
        def get_basis(self):
            return {"nbf": 2}

    class FakeMol:
        def __init__(self):
            self.data = FakeData()

    def int2e(mol):
        calls["n"] += 1
        mol.data["OQP::ERI_AO"] = np.full(16, calls["n"], dtype=float)

    setattr(fake_oqp, "int2e", int2e)
    monkeypatch.setitem(sys.modules, "oqp", fake_oqp)

    spec = importlib.util.spec_from_file_location(
        "_oqp_library_ints_2e_under_test",
        os.path.join(_HERE, os.pardir, "pyoqp", "oqp", "library", "ints_2e.py"))
    assert spec is not None
    assert spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    mol = FakeMol()
    mol.data["OQP::ERI_AO"] = np.zeros(16)

    first = mod.eri_ao(mol)
    second = mod.eri_ao(mol)

    assert calls["n"] == 2
    assert np.all(first == 1.0)
    assert np.all(second == 2.0)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
