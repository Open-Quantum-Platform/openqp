"""Algebra/layout gates for resident MRSF NAC pair contractions."""

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
GRADIENT = ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"
ENERGY = ROOT / "source" / "modules" / "tdhf_mrsf_energy.F90"
INTERCHANGE = ROOT / "source" / "modules" / "mrsf_nac_interchange.F90"
SCF_ADDONS = ROOT / "source" / "scf_addons.F90"
DRIVER = ROOT / "source" / "modules" / "mrsf_nac_driver.F90"


def _body(source, name):
    start = source.rindex(f"subroutine {name}(")
    end = source.index(f"end subroutine {name}", start)
    return source[start:end]


def test_wpair_restricts_products_to_the_nonzero_iatogen_blocks():
    body = _body(GRADIENT.read_text(), "mrsf_nac_wpair_impl")
    assert "fb(1,nocb+1)" in body
    assert "gamma_b(nocb+1,nocb+1)" in body
    assert "g(:,nocc_b+1:n)" in body
    assert "v(:,nocc_b+1:n)" in body
    assert "hb(:,nocc_b+1:n)" in body
    assert "hx_tmp(nbf,nbf), hx_g(nbf,nbf), hx_f7(nbf,nbf)" in body
    assert "allocate(tmp(n,n)" not in body

    rng = np.random.default_rng(12)
    nbf, noca, nocb = 11, 6, 4
    nvirb = nbf - nocb
    fa = rng.standard_normal((nbf, nbf))
    fb = rng.standard_normal((nbf, nbf))
    gamma_a = np.zeros((nbf, nbf))
    gamma_b = np.zeros((nbf, nbf))
    gamma_a[:noca, :noca] = rng.standard_normal((noca, noca))
    gamma_b[nocb:, nocb:] = rng.standard_normal((nvirb, nvirb))
    mt = rng.standard_normal((nbf, nbf))
    full = mt + 2.0 * (fa @ gamma_a + fb @ gamma_b)
    blocked = mt.copy()
    blocked[:, :noca] += 2.0 * fa[:, :noca] @ gamma_a[:noca, :noca]
    blocked[:, nocb:] += 2.0 * fb[:, nocb:] @ gamma_b[nocb:, nocb:]
    np.testing.assert_allclose(blocked, full, rtol=0.0, atol=3.0e-14)

    g = rng.standard_normal((nbf, nbf))
    v = np.zeros((nbf, nbf))
    v[:noca, nocb:] = rng.standard_normal((noca, nvirb))
    ha_full = 2.0 * g @ v.T[:, :noca]
    hb_full = 2.0 * g.T @ v
    ha_blocked = 2.0 * g[:, nocb:] @ v[:noca, nocb:].T
    hb_blocked = np.zeros((nbf, nbf))
    hb_blocked[:, nocb:] = 2.0 * g[:noca, :].T @ v[:noca, nocb:]
    np.testing.assert_allclose(ha_blocked, ha_full, rtol=0.0, atol=3.0e-14)
    np.testing.assert_allclose(hb_blocked, hb_full, rtol=0.0, atol=3.0e-14)


def test_response_exports_jk_from_its_existing_integral_pass():
    response = _body(ENERGY.read_text(), "mrsf_nac_response")
    addons = _body(SCF_ADDONS.read_text(), "get_response_packed")
    assert "call fock_jk" not in response
    assert "mo_b_work, vjk)" in response
    assert "v1_tri, mo_b, vjk_tri)" in addons
    assert addons.count("if (present(vjk_tri)) vjk_tri = v1_tri") == 2


def test_pair_overlap_reverse_transforms_weights_once_per_pair():
    body = _body(INTERCHANGE.read_text(), "mrsf_nac_rohf_pair_overlap")
    coordinate_loop = body.split("do atom = 1, natom", 1)[1]
    assert body.count("call dgemm") == 4
    assert "call dgemm" not in coordinate_loop
    assert "value = sum(overlap_weight_ao*dsfull" in coordinate_loop
    assert "gsk = sum(gamma_ao*dsket" in coordinate_loop

    rng = np.random.default_rng(31)
    nbf, nocb, noca = 9, 3, 5
    mo = np.linalg.qr(rng.standard_normal((nbf, nbf)))[0]
    xmat = rng.standard_normal((nbf, nbf))
    gamma = rng.standard_normal((nbf, nbf))
    norms = np.exp(rng.normal(scale=0.1, size=nbf))
    dsket = rng.standard_normal((nbf, nbf))
    dsfull = rng.standard_normal((nbf, nbf))
    normalized_ket = dsket * norms[:, None] * norms[None, :]
    normalized_full = dsfull * norms[:, None] * norms[None, :]
    skmo = mo.T @ normalized_ket @ mo
    sxmo = mo.T @ normalized_full @ mo

    def orbital_space(index):
        if index < nocb:
            return 1
        if index < noca:
            return 2
        return 3

    original_value = 0.0
    weight = np.zeros((nbf, nbf))
    for p in range(nbf):
        original_value -= 0.5 * xmat[p, p] * sxmo[p, p]
        weight[p, p] = -0.5 * xmat[p, p]
        for q in range(p):
            if orbital_space(p) == orbital_space(q):
                original_value -= 0.5 * (xmat[p, q] + xmat[q, p]) * sxmo[p, q]
                coefficient = -0.5 * (xmat[p, q] + xmat[q, p])
            else:
                original_value -= xmat[q, p] * sxmo[p, q]
                coefficient = -xmat[q, p]
            weight[p, q] = coefficient
    original_gsk = np.sum(gamma * skmo)

    weight_ao = mo @ weight @ mo.T * norms[:, None] * norms[None, :]
    gamma_ao = mo @ gamma @ mo.T * norms[:, None] * norms[None, :]
    np.testing.assert_allclose(
        np.sum(weight_ao * dsfull), original_value, rtol=0.0, atol=5.0e-14
    )
    np.testing.assert_allclose(
        np.sum(gamma_ao * dsket), original_gsk, rtol=0.0, atol=5.0e-14
    )


def test_metric_only_reverse_reproduces_ordered_pair_antisymmetry():
    overlap = _body(INTERCHANGE.read_text(), "mrsf_nac_rohf_pair_overlap")
    driver = _body(DRIVER.read_text(), "mrsf_nac_lagrangian")
    assert "if (.not. only_metric) then" in overlap
    assert "xmat = gamma" in overlap
    direct_record_reads = overlap.split("if (.not. only_metric) then", 1)[1]
    assert "tag_mt_frozen" in direct_record_reads
    assert "tag_mt_response" in direct_record_reads
    assert "metric_only=.true." in driver
    assert "gamma_pair = pair_sign*gamma_column(:,istate)" in driver

    # Let R and O stand for arbitrary linear RHS and coordinate projections.
    # The identity therefore proves the production folding independently of
    # the particular ROHF tangent and overlap contraction implementations.
    rng = np.random.default_rng(91)
    nmat, nrhs, ncoord = 13, 7, 9
    direct = rng.standard_normal(nmat)
    gamma_ij = rng.standard_normal(nmat)
    gamma_ji = rng.standard_normal(nmat)
    rhs_projection = rng.standard_normal((nrhs, nmat))
    overlap_projection = rng.standard_normal((ncoord, nmat))
    direct_coordinates = rng.standard_normal(ncoord)

    ordered_ij = direct + gamma_ij
    ordered_ji = -direct + gamma_ji
    old_rhs = 0.5 * rhs_projection @ (ordered_ij - ordered_ji)
    new_rhs = (
        rhs_projection @ (direct + 0.5 * gamma_ij)
        + rhs_projection @ (-0.5 * gamma_ji)
    )
    old_coordinates = 0.5 * (
        (direct_coordinates + overlap_projection @ ordered_ij)
        - (-direct_coordinates + overlap_projection @ ordered_ji)
    )
    new_coordinates = (
        direct_coordinates
        + overlap_projection @ (direct + 0.5 * gamma_ij)
        + overlap_projection @ (-0.5 * gamma_ji)
    )
    np.testing.assert_allclose(new_rhs, old_rhs, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(
        new_coordinates, old_coordinates, rtol=0.0, atol=2.0e-14
    )
