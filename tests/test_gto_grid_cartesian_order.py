"""The Cartesian component order of ``oqp.analysis.gto_grid`` must be the
engine's.  ``_CART[3]`` once held the Molden 10F label order
(xxx,yyy,zzz,xyy,xxy,...) instead of the engine order of CART_X/Y/Z in
source/constants.F90 (xxx,yyy,zzz,xxy,xxz,xyy,...); s/p/d were unaffected, but
every f AO was assigned the wrong monomial, so ``overlap_analytic()`` differed
from ``OQP::SM`` by up to 0.66 and every cube/AICD/descriptor evaluated on a
basis with f shells was wrong.  g and h use the same engine table and are
checked the same way.

The table and guard tests run without the compiled library; the overlap tests need it.
"""
import importlib.util
import os
import re

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def _load_gto_grid():
    spec = importlib.util.spec_from_file_location(
        "_cart_order_gto_grid", os.path.join(_ROOT, "pyoqp/oqp/analysis/gto_grid.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _engine_cart_table():
    """Parse CART_X/Y/Z from source/constants.F90 into {L: [(lx, ly, lz), ...]}."""
    with open(os.path.join(_ROOT, "source", "constants.F90")) as fh:
        text = fh.read()
    axes = []
    for name in ("CART_X", "CART_Y", "CART_Z"):
        block = re.search(name + r"\(BAS_MXCART,0:BAS_MXANG\) = reshape\(\[(.*?)\], shape\(" + name,
                          text, re.S).group(1)
        rows = []
        for line in block.splitlines():
            m = re.match(r"\s*\[([0-9,\s]+?)(?:,\s*\(0,|\])", line)
            if m:
                rows.append([int(v) for v in m.group(1).replace(" ", "").split(",") if v])
        axes.append(rows)
    return {L: list(zip(*(ax[L] for ax in axes))) for L in range(len(axes[0]))}


def test_cart_table_matches_engine_constants():
    gto_grid = _load_gto_grid()
    engine = _engine_cart_table()
    for L, comps in gto_grid._CART.items():
        assert comps == engine[L], f"L={L}: gto_grid {comps} != engine {engine[L]}"


def _unpack_lower(packed, n):
    S = np.zeros((n, n))
    S[np.tril_indices(n)] = packed
    return S + np.tril(S, -1).T


# (basis, highest L, number of Cartesian AOs with that L) for off-axis CO.
_OVERLAP_CASES = [("cc-pvtz", 3, 20), ("cc-pvqz", 4, 30), ("cc-pv5z", 5, 42)]


@pytest.mark.parametrize("basis,lmax,n_lmax", _OVERLAP_CASES,
                         ids=[c[0] for c in _OVERLAP_CASES])
def test_overlap_analytic_reproduces_oqp_sm(basis, lmax, n_lmax, tmp_path, monkeypatch):
    pytest.importorskip("oqp")
    from oqp.openqp import OPENQP
    from oqp.analysis.gto_grid import AOBasis
    import oqp.library

    monkeypatch.chdir(tmp_path)
    # Off-axis CO so no Cartesian component is spared by symmetry; each basis
    # puts an L=lmax shell on both atoms, so one-centre and two-centre blocks
    # of the highest shell are both tested.
    cfg = {"input.runtype": "energy", "input.method": "hf", "input.basis": basis,
           "input.ispher": "false", "input.charge": 0,
           "input.system": "C 0.10 -0.20 0.05; O 0.72 0.31 0.83",
           "guess.type": "huckel"}
    mol = OPENQP(cfg, True).mol
    oqp.library.set_basis(mol)
    oqp.library.ints_1e(mol)

    ao = AOBasis(mol)
    n = ao.nbf
    ls = [lx + ly + lz for _, (lx, ly, lz), _ in ao.ao_index]
    assert max(ls) == lmax and ls.count(lmax) == n_lmax
    S_ref = _unpack_lower(np.asarray(mol.data["OQP::SM"], dtype=float).ravel(), n)
    S = ao.overlap_analytic()
    err = np.abs(S - S_ref)
    assert err.max() < 1e-10, (
        f"{basis}: max |S_analytic - OQP::SM| = {err.max():.3e} "
        f"at {np.unravel_index(err.argmax(), err.shape)}")


class _FakeData:
    def __init__(self, basis):
        self._b = basis

    def get_basis(self):
        return self._b


class _FakeMol:
    def __init__(self, basis):
        self.data = _FakeData(basis)

    def get_system(self):
        return np.zeros(3)


def test_shell_beyond_table_is_rejected():
    gto_grid = _load_gto_grid()
    L = max(gto_grid._CART) + 1
    basis = {"nbf": (L + 1) * (L + 2) // 2, "nsh": 1, "centers": np.array([0]),
             "angs": np.array([L]), "ncontr": np.array([1]),
             "alpha": np.array([1.0]), "coef": np.array([1.0])}
    with pytest.raises(NotImplementedError, match=f"L={L} shell"):
        gto_grid.AOBasis(_FakeMol(basis))
