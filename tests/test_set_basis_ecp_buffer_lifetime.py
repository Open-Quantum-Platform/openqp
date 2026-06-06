"""Regression tests for the dangling ecp_zn/ecp_am CFFI buffers (1BNA bug).

set_ecp_data() used to build the ecp_zn/ecp_am pointers from temporary numpy
arrays: ``ffi.cast("int*", ffi.from_buffer(np.array(...)))``.  ffi.from_buffer
does not keep the array alive, so the buffer was freed before
``oqp.append_ecp`` copied it in Fortran.  For large systems (758 atoms in the
1BNA benchmark) the freed block was reused, ecp_zn_num filled with garbage,
and the nuclear repulsion energy ``e_charge_repulsion(xyz, zn - ecp_zn_num)``
became a huge constant — SCF energies printed as ``*****`` from iteration 1
while DIIS still converged.  Introduced by commit d1398dfe (ECP, 2025-02-18).
"""

import gc
import importlib
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SET_BASIS = ROOT / "pyoqp" / "oqp" / "library" / "set_basis.py"


def import_oqp_or_skip():
    """Import oqp, skipping cleanly when native build artifacts are absent."""
    try:
        return importlib.import_module("oqp")
    except (ImportError, FileNotFoundError, OSError) as exc:
        pytest.skip(f"oqp native runtime is not available: {exc}")


def test_no_inline_temporary_buffers_in_set_basis():
    """ffi.from_buffer(np.array(...)) hands CFFI a temporary — forbidden."""
    source = SET_BASIS.read_text()
    assert "ffi.from_buffer(np.array(" not in source


def test_append_ecp_sees_live_no_ecp_zn_buffer(monkeypatch):
    """ecp_zn data must still be intact when append_ecp reads the no-ECP path."""
    oqp = import_oqp_or_skip()
    from oqp import ffi
    from oqp.library.set_basis import BasisData

    natom = 758
    captured = {}

    def fake_append_ecp(mol):
        # Emulate allocator reuse of any freed temporary before the Fortran
        # side would copy the buffer.
        gc.collect()
        churn = [np.full(natom, 7, dtype=np.int32) for _ in range(64)]
        buf = ffi.buffer(mol.data["ecp_zn"], natom * 4)
        captured["ecp_zn"] = np.frombuffer(bytes(buf), dtype=np.int32).copy()
        del churn

    monkeypatch.setattr(oqp, "append_ecp", fake_append_ecp)

    class FakeMol:
        data = {}

    basis_data = object.__new__(BasisData)
    basis_data.mol = FakeMol()
    basis_data.use_ecp = False
    basis_data.ecp = {
        "ang": [],
        "r_expo": [],
        "g_expo": [],
        "coef": [],
        "coord": [],
        "element_id": 0,
        "ecp_electron": [0] * natom,
        "num_expo": [],
    }

    basis_data.set_ecp_data()

    assert captured["ecp_zn"].shape == (natom,)
    # No ECP: every entry must be exactly zero, not recycled-heap garbage.
    assert np.all(captured["ecp_zn"] == 0)


def test_append_ecp_sees_live_ecp_am_buffer(monkeypatch):
    """ecp_am data must also survive until append_ecp reads active-ECP data."""
    oqp = import_oqp_or_skip()
    from oqp import ffi
    from oqp.library.set_basis import BasisData

    natom = 4
    expected_am = np.array([0.0, 1.0, 2.0], dtype=np.float64)
    captured = {}

    def fake_append_ecp(mol):
        gc.collect()
        churn_i = [np.full(natom, 7, dtype=np.int32) for _ in range(64)]
        churn_f = [np.full(expected_am.size, -9.0, dtype=np.float64) for _ in range(64)]
        ecp_zn_buf = ffi.buffer(mol.data["ecp_zn"], natom * 4)
        ecp_am_buf = ffi.buffer(mol.data["ecp_am"], expected_am.nbytes)
        captured["ecp_zn"] = np.frombuffer(bytes(ecp_zn_buf), dtype=np.int32).copy()
        captured["ecp_am"] = np.frombuffer(bytes(ecp_am_buf), dtype=np.float64).copy()
        del churn_i, churn_f

    monkeypatch.setattr(oqp, "append_ecp", fake_append_ecp)

    class FakeMol:
        data = {}

    basis_data = object.__new__(BasisData)
    basis_data.mol = FakeMol()
    # Keep this false so the test only exercises append_ecp buffer handoff;
    # the post-append electron-count update needs a full native molecule.
    basis_data.use_ecp = False
    basis_data.ecp = {
        "ang": expected_am.tolist(),
        "r_expo": [2.0, 2.0, 2.0],
        "g_expo": [1.1, 1.2, 1.3],
        "coef": [0.1, 0.2, 0.3],
        "coord": [0.0, 0.0, 0.0],
        "element_id": 1,
        "ecp_electron": [0, 2, 0, 0],
        "num_expo": [3],
    }

    basis_data.set_ecp_data()

    np.testing.assert_array_equal(captured["ecp_zn"], np.array([0, 2, 0, 0], dtype=np.int32))
    np.testing.assert_array_equal(captured["ecp_am"], expected_am)
