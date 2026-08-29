"""Guards against dangling TagArray views in the resident NAC traversal."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INTERCHANGE = ROOT / "source" / "modules" / "mrsf_nac_interchange.F90"
GRADIENT = ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"
DRIVER = ROOT / "source" / "modules" / "mrsf_nac_driver.F90"
ZVECTOR = ROOT / "source" / "modules" / "tdhf_mrsf_z_vector.F90"
FOCK_DERIV = ROOT / "source" / "modules" / "fock_deriv.F90"


def _body(source, name):
    return source.split(f"subroutine {name}", 1)[1].split(
        f"end subroutine {name}", 1
    )[0]


def test_pair_overlap_reacquires_mos_after_output_reservations():
    body = _body(
        INTERCHANGE.read_text(),
        "mrsf_nac_rohf_pair_overlap(infos, metric_only)",
    )
    last_reserve = body.index("tagarray_reserve_data(infos%dat, tag_gsk")
    reacquire = body.index("tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo)", last_reserve)
    first_transform = body.index("call der_overlap_matrix_ket")
    assert last_reserve < reacquire < first_transform


def test_pair_finalize_reacquires_inputs_after_output_reservations():
    body = _body(INTERCHANGE.read_text(), "mrsf_nac_pair_finalize(infos)")
    last_reserve = body.index("tagarray_reserve_data(infos%dat, tag_nacv")
    dp_reacquire = body.index("tagarray_get_data(infos%dat, tag_dp, dp_ordered)", last_reserve)
    energy_reacquire = body.index(
        "tagarray_get_data(infos%dat, OQP_td_energies, energies)", last_reserve
    )
    first_use = body.index("gap = energies(jstate) - energies(istate)")
    assert last_reserve < dp_reacquire < first_use
    assert last_reserve < energy_reacquire < first_use


def test_amp_and_esum_own_inputs_that_survive_record_mutation():
    source = GRADIENT.read_text()
    amp = _body(source, "mrsf_nac_amp(infos, only_istate, only_jstate)")
    esum = _body(source, "mrsf_nac_esum(infos, istate, jstate)")
    for body in (amp, esum):
        assert "bvec_mo_owned = bvec_mo" in body
        assert "mo_a_owned = mo_a" in body
        assert "bvec_mo => bvec_mo_owned" in body
        assert "mo_a => mo_a_owned" in body
    assert "dmat_a_owned = dmat_a" in amp
    assert "dmat_a => dmat_a_owned" in amp
    assert "mo_b_owned = mo_b" in esum
    assert "dmat_b_owned = dmat_b" in esum


def test_driver_owns_energy_record_and_guards_integer_products_before_cast():
    body = _body(DRIVER.read_text(), "mrsf_nac_lagrangian(infos)")
    assert "energies_saved = energies" in body
    assert "gap = energies_saved(pair_j(batch_pair))" in body
    assert "energies_saved(pair_i(batch_pair))" in body
    cast = body.index("nstate = int(nstate64)")
    for guard in (
        "default_int_limit64/state_pair_size64",
        "default_int_limit64/nij64",
        "default_int_limit64/nbfsq64",
    ):
        assert body.index(guard) < cast


def test_tagarray_contains_results_are_used_as_logicals():
    gradient = GRADIENT.read_text()
    zvector = ZVECTOR.read_text()
    assert "have_custom = infos%dat%contains(tags_gamma, gtag_id)" in gradient
    assert "gstat = infos%dat%contains" not in gradient
    assert "have_orbgrad = infos%dat%contains(tag_L, ltag_id)" in zvector
    assert "have_gamma = infos%dat%contains(tags_gamma, gtag_id)" in zvector
    assert "gstat = infos%dat%contains" not in zvector
    assert "lstat = infos%dat%contains" not in zvector
    assert "gstat = -999" not in zvector


def test_batched_open_shell_fock_probes_prepare_spherical_cartesian_views():
    source = FOCK_DERIV.read_text()
    body = source.split("subroutine fock_deriv_contract_os_batch(", 1)[1].split(
        "end subroutine fock_deriv_contract_os_batch", 1
    )[0]
    assert "if (HARMONIC_ACTIVE) then" in body
    for view in (
        "gcomps(ia)%pcoul_cart",
        "gcomps(ia)%pexch_cart",
        "gcomps(ia)%mmat_cart",
        "gcomps(ib)%pcoul_cart",
        "gcomps(ib)%pexch_cart",
        "gcomps(ib)%mmat_cart",
    ):
        assert view in body
    assert "gcomps(ia)%cart_off" in body
    assert "gcomps(ib)%cart_off" in body
