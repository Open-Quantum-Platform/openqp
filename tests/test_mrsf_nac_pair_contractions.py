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
    body = _body(GRADIENT.read_text(), "mrsf_nac_wpair_batch_impl")
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


def test_wpair_batches_three_pairs_in_one_bounded_eri_traversal():
    gradient = GRADIENT.read_text()
    batch = _body(gradient, "mrsf_nac_wpair_batch_impl")
    driver = _body(DRIVER.read_text(), "mrsf_nac_lagrangian")
    assert "integer, parameter :: max_batch_width = 3" in batch
    assert "mrsf_density(2*nrhs,7,nbf,nbf)" in batch
    assert "d3 = mrsf_density(:2*nrhs,:,:,:)" in batch
    assert batch.count("call int2_driver%run") == 1
    assert "source_x = 2*ipair - 1" in batch
    assert "source_y = source_x + 1" in batch
    assert "integer, parameter :: wpair_batch_width = 3" in driver
    assert "wpair_mt(nbf,nbf,wpair_batch_width)" in driver
    assert "wpair_last = min(npair, wpair_first + wpair_batch_width - 1)" in driver
    assert "call mrsf_nac_wpair_batch_impl(" in driver
    assert "mt_frozen_tag = reshape(wpair_mt(:,:,wpair_index)" in driver
    assert "call mrsf_nac_wpair_impl(infos, istate, jstate)" not in driver


def test_response_exports_jk_from_its_existing_integral_pass():
    response = _body(ENERGY.read_text(), "mrsf_nac_response")
    addons = _body(SCF_ADDONS.read_text(), "get_response_packed")
    assert "call fock_jk" not in response
    assert "mo_b_work, vjk)" in response
    assert "v1_tri, mo_b, vjk_tri)" in addons
    assert addons.count("if (present(vjk_tri)) vjk_tri = v1_tri") == 2


def test_amp_subtracts_the_reference_before_the_integral_sweep():
    source = GRADIENT.read_text()
    amp = _body(source, "mrsf_nac_amp")
    density = _body(source, "grd2_mrsf_nac_compute_data_t_get_density")
    assert "subtract_reference = .true." in amp
    assert amp.count("call grd2_driver") == 1
    assert "deSCF" not in amp
    assert "if (this%subtract_reference) then" in density

    rng = np.random.default_rng(87)
    d = rng.standard_normal((2, 8))
    p = rng.standard_normal((2, 8))
    # The first four slots represent the (ik,jl) and (il,jk) index pairs
    # used for the total/spin exchange channels.  Verify the scalar identity
    # implemented in the shell-quartet callback.
    full_c = (d[0, 0] + p[0, 0]) * d[0, 1] + d[0, 0] * p[0, 1]
    base_c = d[0, 0] * d[0, 1]
    direct_c = p[0, 0] * d[0, 1] + d[0, 0] * p[0, 1]
    np.testing.assert_allclose(full_c - base_c, direct_c, rtol=0.0, atol=1e-15)

    full_x = 0.0
    base_x = 0.0
    direct_x = 0.0
    for spin in range(2):
        for left, right in ((2, 3), (4, 5)):
            full_x += (d[spin, left] + p[spin, left]) * d[spin, right]
            full_x += d[spin, left] * p[spin, right]
            base_x += d[spin, left] * d[spin, right]
            direct_x += p[spin, left] * d[spin, right]
            direct_x += d[spin, left] * p[spin, right]
    np.testing.assert_allclose(full_x - base_x, direct_x, rtol=0.0, atol=3e-15)


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
