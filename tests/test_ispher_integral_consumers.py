from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(relpath):
    return (ROOT / relpath).read_text()


def test_ground_and_response_gradients_expand_spherical_density_for_cartesian_derivatives():
    basis = read("source/basis_tools.F90")
    gradient_modules = {
        "HF/DFT": read("source/modules/hf_gradient.F90"),
        "TDDFT": read("source/modules/tdhf_gradient.F90"),
        "SF-TDDFT": read("source/modules/tdhf_sf_gradient.F90"),
        "MRSF-TDDFT": read("source/modules/tdhf_mrsf_gradient.F90"),
    }

    assert "subroutine build_cart_density" in basis
    assert "call c2s_expand_block(" in basis
    assert "basis%am(ish), basis%harmonic(ish)" in basis
    assert "basis%am(jsh), basis%harmonic(jsh)" in basis

    for name, source in gradient_modules.items():
        assert "build_cart_density" in source, name
        assert "usecart = HARMONIC_ACTIVE" in source, name
        assert "nbf = NUM_CART_BF(basis%am(id))" in source, name
        assert "nbf = basis%naos(id)" in source, name


def test_dft_grid_ao_builders_reduce_values_gradients_and_hessians_to_spherical():
    basis = read("source/basis_tools.F90")
    dft_grid = read("source/dftlib/dft_gridint.F90")
    dft_grad = read("source/dftlib/dft_gridint_grad.F90")
    dft_tdxc_grad = read("source/dftlib/dft_gridint_tdxc_grad.F90")

    assert "call cart2sph_vec(cbuf, sbuf, am)" in basis
    assert "call cart2sph_vec(cv, sv, am); aov" in basis
    assert "call cart2sph_vec(cx, sx, am); aogx" in basis
    assert "call cart2sph_vec(cb(:,5),  sb, am); aog2xx" in basis
    assert "nao = basis%naos(ish)" in dft_grid
    assert "naos => basis%naos(j)" in dft_grad
    assert "naos => basis%naos(j)" in dft_tdxc_grad


def test_response_energy_nmr_pcm_and_ekt_use_active_spherical_ao_dimensions():
    response_energy = {
        "TDDFT": read("source/modules/tdhf_energy.F90"),
        "SF-TDDFT": read("source/modules/tdhf_sf_energy.F90"),
        "MRSF-TDDFT": read("source/modules/tdhf_mrsf_energy.F90"),
        "TDDFT Z-vector": read("source/modules/tdhf_z_vector.F90"),
        "SF-TDDFT Z-vector": read("source/modules/tdhf_sf_z_vector.F90"),
        "MRSF-TDDFT Z-vector": read("source/modules/tdhf_mrsf_z_vector.F90"),
        "EKT": read("source/modules/tdhf_mrsf_ekt.F90"),
    }
    for name, source in response_energy.items():
        assert "nbf = basis%nbf" in source, name

    pcm = read("source/solvent_pcm.F90")
    assert "allocate(ao_pop(basis%nbf)" in pcm
    assert "traceprod_sym_packed(dtot, vpcm, basis%nbf)" in pcm
    assert "i1 = basis%ao_offset(ish) + basis%naos(ish) - 1" in pcm

    nmr = read("source/modules/nmr_giao_shielding.F90")
    nmr_giao_twoe = read("source/modules/nmr_giao_debug.F90")
    int1 = read("source/integrals/int1.F90")
    assert "nbf = basis%nbf" in nmr
    assert "HARMONIC_ACTIVE .and. any(basis%harmonic == 1)" not in nmr
    assert "call build_cart_density(basis, dm_norm, dm_cart" in nmr_giao_twoe
    assert "call reduce_cart_giao_matrix(basis, cart_off" in nmr_giao_twoe
    assert "ni = NUM_CART_BF(basis%am(si))" in nmr_giao_twoe
    assert "mapr = raised_cart_index(ii, am0(1), axis)" in nmr_giao_twoe
    assert "CALL cart2sph_mat(blk(:,m), shj%ang, shj%harmonic, shi%ang, shi%harmonic" in int1
