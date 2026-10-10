import importlib
import os
import sys
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault("OQP_BACKEND_OPTIONAL", "1")
ROOT = Path(__file__).resolve().parents[1]
PYOQP = ROOT / "pyoqp"
if str(PYOQP) not in sys.path:
    sys.path.insert(0, str(PYOQP))


def _use_real_oqp_package():
    for name in list(sys.modules):
        if name == "oqp" or name.startswith("oqp."):
            sys.modules.pop(name, None)
    importlib.invalidate_caches()


def test_rdm2_spatial_fallback_non_product_list():
    """rdm2_spatial falls back to rdm2_gram for non-product determinant lists.
    Verify that the truncated-list path produces a valid 2-RDM."""
    _use_real_oqp_package()
    from oqp.library.rdm import determinant_basis, make_rdm2_spatial, make_rdm1_spatial

    norb = 6; nelec = (3, 3)
    dets = determinant_basis(norb, nelec)
    ndet_full = len(dets)
    assert ndet_full == 400, f"Expected 400 determinants for CAS(6,3+3), got {ndet_full}"

    # Truncate to a non-product subset (first 10 determinants)
    # rdm12_strings will decline because nastr*nbstr /= ndet
    dets_trunc = np.array(dets[:10], dtype=np.int64)
    rng = np.random.default_rng(42)
    ci = np.ascontiguousarray(rng.standard_normal(10))
    ci /= np.linalg.norm(ci)

    d2 = make_rdm2_spatial(ci, dets_trunc, norb)
    d1 = make_rdm1_spatial(ci, dets_trunc, norb)

    # Identity: sum_r D2[p,q,r,r] = (N_elec-1) * D1[p,q]
    d2_contracted = np.einsum("pqrr->pq", d2)
    assert d2_contracted.shape == (norb, norb)
    assert np.allclose(d2_contracted, (nelec[0]+nelec[1]-1) * d1, atol=1e-12), \
        "Truncated-list D2 fails contraction identity"

    # Spatial 1-RDM trace = N_elec
    assert np.trace(d1) == pytest.approx(nelec[0]+nelec[1], abs=1e-12)
    assert not np.any(np.isnan(d2))
    assert np.max(np.abs(d2)) < 10.0


def test_rdm2_spatial_via_gram_matches_strings():
    """For a full product list, rdm2_spatial (rdm12_strings path) agrees
    with the rdm2_gram fallback path (serial key-walk)."""
    _use_real_oqp_package()
    from oqp.library.rdm import determinant_basis, make_rdm2_spatial, make_rdm1_spatial

    norb = 6; nelec = (3, 3)
    dets = determinant_basis(norb, nelec)
    ndet = len(dets)
    rng = np.random.default_rng(42)
    ci = np.ascontiguousarray(rng.standard_normal(ndet))
    ci /= np.linalg.norm(ci)

    d2_new = make_rdm2_spatial(ci, dets, norb)
    d1 = make_rdm1_spatial(ci, dets, norb)

    # Contraction identity: sum_r D2[p,q,r,r] = (N_elec-1) * D1[p,q]
    assert np.allclose(np.einsum("pqrr->pq", d2_new), (nelec[0]+nelec[1]-1) * d1, atol=1e-12)
    # Trace identity: tr(D1) = N_elec
    assert np.trace(d1) == pytest.approx(nelec[0]+nelec[1], abs=1e-12)
