"""Regression gates for the native spin-adapted MRSF sigma refactor.

The fast tests protect the production representation and the memory-bounded
dense-column driver.  The opt-in live test compares the RouteC-disabled native
calculation with the pre-refactor H2O/BHHLYP reference.  Run the latter against
the freshly built package with::

    OPENQP_RUN_MRSF_SIGMA_REGRESSION=1 OPENQP_ROOT=/path/to/installed/oqp \
      pytest -q tests/test_mrsf_sigma_refactor.py
"""

from __future__ import annotations

import json
import os
import platform
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SIGMA = ROOT / "source/modules/tdhf_mrsf_sigma.F90"
OPERATOR = ROOT / "source/modules/tdhf_mrsf_hessian_operator.F90"
ENERGY = ROOT / "source/modules/tdhf_mrsf_energy.F90"
DENSITY_CONTRACTION = ROOT / "source/modules/tdhf_mrsf_density_contraction.F90"
REFERENCE_INPUT = ROOT / "examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.inp"
REFERENCE_JSON = ROOT / "examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.json"
RUN_LIVE = os.environ.get("OPENQP_RUN_MRSF_SIGMA_REGRESSION") == "1"


def _compact(path: Path) -> str:
    return "".join(path.read_text().lower().split())


def _runtime_available() -> bool:
    if not RUN_LIVE:
        return False
    suffix = {
        "Darwin": "dylib",
        "Linux": "so",
        "Windows": "dll",
    }.get(platform.system(), "so")
    runtime = Path(os.environ.get("OPENQP_ROOT", str(ROOT)))
    return (runtime / "lib" / f"liboqp.{suffix}").exists()


def test_production_sigma_is_only_the_spin_adapted_seven_density_chain():
    sigma = _compact(SIGMA)
    assert "density(nvec,7,nbf,nbf)" in sigma
    for routine in ("mrsfcbc", "mrsfmntoia", "mrsfesum"):
        assert f"call{routine}" in sigma
    for forbidden in (
        "slater",
        "determinant",
        "fock_space",
        "fci_hamiltonian",
        "build_mrsf_determinant_transform",
    ):
        assert forbidden not in sigma


def test_scale_resolution_and_pure_semilocal_override_use_grouped_channels():
    sigma = _compact(SIGMA)
    energy = _compact(ENERGY)
    assert "callprepare_mrsf_response_scaling" in sigma
    assert "callprepare_mrsf_response_scaling" in energy
    assert "infos%tddft%hfscale=infos%tddft%cam_alpha" in sigma
    assert "infos%tddft%spc_coco=response_scale" in sigma
    assert "grouped_channels=abs(response_scale)<=epsilon(1.0_dp).and.nonzero_spc" in sigma
    assert "callcontract_mrsf_seven_density_batch" in sigma


def test_zero_alpha_cam_spc_override_is_rejected_separately():
    sigma = _compact(SIGMA)
    assert "if(grouped_channels.and.dft.and.infos%dft%cam_flag)then" in sigma
    assert "status=-3" in sigma
    assert "zero-hfscalespcfallbackdoesnotsupportarange-separatedcamresponse" in sigma


def test_grouped_contraction_releases_result_pointer_before_cleaning_target():
    contraction = _compact(DENSITY_CONTRACTION)
    result_association = contraction.index("result=>data%f3")
    nullify = contraction.index("nullify(result)", result_association)
    clean = contraction.index("calldata%clean()", result_association)
    assert result_association < nullify < clean
    assert ".not.allocated(data%f3)" in contraction
    assert "size(data%f3,2)/=7" in contraction


def test_nonzero_exchange_preserves_the_common_seven_density_screening_batch():
    contraction = _compact(DENSITY_CONTRACTION)
    coupled = contraction.index("if(abs(response_scale)>epsilon(1.0_dp))then")
    grouped = contraction.index("callrun_group(1,4,spc_coov)", coupled)
    assert contraction.index("callrun_coupled(response_scale)", coupled) < grouped
    for scaling in (
        "contracted(:,1:4,:,:)=contracted(:,1:4,:,:)*spc_coov/response_scale",
        "contracted(:,5:5,:,:)=contracted(:,5:5,:,:)*spc_ovov/response_scale",
        "contracted(:,6:6,:,:)=contracted(:,6:6,:,:)*spc_coco/response_scale",
    ):
        assert scaling in contraction


def test_semilocal_grouped_contractions_keep_the_full_screening_envelope():
    contraction = _compact(DENSITY_CONTRACTION)
    assert "envelope=maxval(abs(density),dim=2)" in contraction
    assert "group_density(:,7,:,:)=envelope" in contraction


def test_sigma_validates_dimensions_and_cleans_pointer_targets():
    sigma = _compact(SIGMA)
    for guard in (
        "infos%basis%nbf/=nbf",
        "noccb<0",
        "nocca<2",
        "nocca>nbf",
        "nocca-noccb/=2",
        ".not.associated(int2_driver%basis)",
        ".not.allocated(int2_driver%schwarz_ints_regular)",
    ):
        assert guard in sigma
    assert "stat=alloc_status" in sigma
    assert "nullify(contracted)" in sigma
    assert "callint2_data%clean()" in sigma
    assert sigma.index("callint2_data%clean()") < sigma.index("deallocate(density,work)", sigma.index("callint2_data%clean()"))


def test_energy_has_no_duplicate_base_mrsf_density_or_dead_legacy_sigma():
    energy = _compact(ENERGY)
    assert "int2_mrsf_data_t" not in energy
    assert "callmrsfcbc" not in energy
    assert "callmrsfmntoia" not in energy
    assert "mrsf_density(nvec,7,nbf,nbf)" not in energy
    assert "callapply_mrsf_tda_sigma" in energy
    assert "mrsf_density(nvec,11,nbf,nbf)" in energy


def test_dense_operator_batches_columns_with_a_bounded_default():
    operator = _compact(OPERATOR)
    assert "infos%basis%nbf/=nbf" in operator
    assert "batch=1" in operator
    assert "dofirst=1,physical,batch" in operator
    assert "packed_sigma(packed,batch)" in operator
    assert "packed_basis(:,first:last)" in operator
    assert "packed_sigma(:,1:ncol)" in operator
    assert "packed_sigma(packed,physical)" not in operator


@pytest.mark.skipif(
    not _runtime_available(),
    reason=(
        "set OPENQP_RUN_MRSF_SIGMA_REGRESSION=1 and OPENQP_ROOT to the "
        "freshly installed OpenQP package"
    ),
)
def test_live_routec_off_native_sigma_matches_pre_refactor_energy(tmp_path, monkeypatch):
    """Exercise the refactored native sigma against the committed old-path result."""
    # An absent variable selects the native path.  The RouteC loader otherwise
    # interprets the value as a library name, so the literal string ``off``
    # would test a failed dlopen before falling back rather than a clean native
    # execution.
    monkeypatch.delenv("OQP_ROUTEC_SIG", raising=False)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "1")
    monkeypatch.setenv("MKL_NUM_THREADS", "1")

    from oqp.pyoqp import Runner

    reference = np.asarray(json.loads(REFERENCE_JSON.read_text())["td_energies"])
    runner = Runner(
        project="mrsf_sigma_refactor",
        input_file=str(REFERENCE_INPUT),
        log=str(tmp_path / "mrsf_sigma_refactor.log"),
        usempi=False,
    )
    runner.run()
    actual = np.asarray(runner.mol.data["OQP::td_energies"], dtype=float)

    assert actual.shape == reference.shape
    # The committed deck uses SCF/TD convergence 1e-6; separate clean builds can
    # consequently move a converged root by O(1e-8) even with the same response
    # operator.  Keep this cross-build regression below 5% of that requested
    # threshold.  Algebraic equality is covered independently by the seven-
    # density oracle tests.
    np.testing.assert_allclose(actual, reference, rtol=0.0, atol=5.0e-8)


@pytest.mark.skipif(
    not _runtime_available(),
    reason=(
        "set OPENQP_RUN_MRSF_SIGMA_REGRESSION=1 and OPENQP_ROOT to the "
        "freshly installed OpenQP package"
    ),
)
def test_live_pure_semilocal_nonzero_spc_is_finite(tmp_path, monkeypatch):
    """Exercise the four independent channel groups without SPC/HFscale."""
    monkeypatch.delenv("OQP_ROUTEC_SIG", raising=False)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "1")
    monkeypatch.setenv("MKL_NUM_THREADS", "1")

    from oqp.pyoqp import Runner

    inp = tmp_path / "pure_semilocal_spc.inp"
    text = REFERENCE_INPUT.read_text().replace("functional=bhhlyp", "functional=pbe")
    text += "\nspc_coco=0.20\nspc_ovov=0.20\nspc_coov=0.20\n"
    inp.write_text(text)
    runner = Runner(
        project="pure_semilocal_spc",
        input_file=str(inp),
        log=str(tmp_path / "pure_semilocal_spc.log"),
        usempi=False,
    )
    runner.run()
    energies = np.asarray(runner.mol.data["OQP::td_energies"], dtype=float)
    assert energies.size == 3
    assert np.all(np.isfinite(energies))
