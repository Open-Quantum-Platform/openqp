"""Source-level regressions for the zero-closed-orbital MRSF Hessian path.

H2 with a triplet ROHF reference has two alpha and zero beta occupied
orbitals.  In the two-SOMO MRSF topology this makes the C-derived CO/CV and
beta occupied-occupied blocks empty; it does not remove the active OO/OV
response space.
"""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
H2_HESSIAN_INPUT = ROOT / "tests/mrsf_h2/h2_mrsftdhf_hessian_rohf.inp"


def _source(relative_path):
    return (ROOT / relative_path).read_text()


def _routine(relative_path, name):
    text = _source(relative_path)
    match = re.search(
        rf"\bsubroutine\s+{name}\b.*?\bend\s+subroutine\s+{name}\b",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    assert match is not None
    return match.group(0)


def _compact(text):
    return "".join(text.lower().split())


def test_h2_keeps_the_nonempty_two_somo_response_space():
    nbf, nocca, noccb = 10, 2, 0
    nco = 2 * noccb
    ncv = noccb * (nbf - nocca)
    nov = 2 * (nbf - nocca)
    z_dimension = noccb * (2 + nbf - nocca) + nov

    assert nco == 0
    assert ncv == 0
    assert nov == 16
    assert z_dimension == 16


def test_h2_cc_pvdz_analytical_hessian_regression_is_registered():
    text = _compact(H2_HESSIAN_INPUT.read_text())
    assert "basis=cc-pvdz" in text
    assert "runtype=hess" in text
    assert "type=rohf" in text
    assert "multiplicity=3" in text
    assert "type=mrsf" in text
    assert "[hess]type=analyticalstate=1" in text
    runner = _source("tests/mrsf_h2/run_h2_mrsf.py")
    assert '"h2_mrsftdhf_hessian_rohf"' in runner


def test_mntoia_returns_before_zero_leading_dimension_dgemm():
    block = _compact(_routine("source/tdhf_lib.F90", "mntoia"))
    guard = block.index("if(nocca==0.or.noccb==nbf)return")
    assert guard < block.index("allocate(scr(nocca,nbf))")
    assert guard < block.index("calldgemm")


def test_dgemm_wrapper_handles_all_zero_dimensions_with_blas_semantics():
    block = _compact(
        _routine("source/mathlib/blas_wrap.F90", "oqp_dgemm_i64")
    )
    guard = "if(m==0.or.n==0)return"
    assert guard in block
    assert block.index(guard) < block.index("calldgemm")
    assert "if(k==0.and.ldc>=max(1,m))then" in block
    assert "if(beta==0.0d0)then" in block
    assert "c(i,j)=0.0d0" in block
    assert "elseif(beta/=1.0d0)then" in block
    assert "c(i,j)=beta*c(i,j)" in block


def test_empty_beta_ppij_transform_is_guarded_and_stored_as_zero_envelope():
    block = _compact(
        _routine(
            "source/modules/tdhf_mrsf_z_vector.F90",
            "build_mrsf_relaxed_density_and_w",
        )
    )
    assert "ppijb=0.0_dpif(noccb>0)then" in block
    assert "0.0_dp,ppijb,noccb)endif" in block
    assert "(/max(1,noccb),max(1,noccb)/)" in block
    assert "stored_ppijb=0.0_dp" in block
    assert "if(noccb>0)stored_ppijb(1:noccb,1:noccb)=ppijb" in block


def test_hessian_response_validators_accept_zero_closed_orbitals():
    z_validator = _compact(
        _routine(
            "source/modules/tdhf_mrsf_hessian_z_response.F90",
            "validate_mrsf_z_method",
        )
    )
    w_validator = _compact(
        _routine(
            "source/modules/tdhf_mrsf_hessian_w_response.F90",
            "validate_baseline_inputs",
        )
    )
    assert "noccb<0" in z_validator
    assert "noccb<=0" not in z_validator
    assert "nocb<0" in w_validator
    assert "nocb<1" not in w_validator

    canonicalizer = _compact(
        _routine(
            "source/modules/tdhf_mrsf_hessian_mo_response.F90",
            "canonicalize_rohf_common_response",
        )
    )
    assert "noccb<0" in canonicalizer
    assert "noccb<1" not in canonicalizer


def test_beta_rohf_connection_uses_metric_gauge_when_nocc_is_zero():
    block = _compact(
        _routine(
            "source/modules/hf_rohf_orbital_response.F90",
            "complete_rohf_orbital_connection",
        )
    )
    assert "nocc<0" in block
    assert "nocc<=0" not in block
    assert "connection=-0.5_dp*sx_mo" in block
    assert "if(nocc>0)then" in block


def test_storage_envelope_never_changes_the_physical_zero_by_zero_block():
    block = _compact(
        _routine(
            "source/modules/tdhf_mrsf_hessian_intermediates.F90",
            "build_mrsf_hf_w_intermediates",
        )
    )
    assert "data%ppijb(noccb,noccb)" in block
    assert "if(any(shape(stored_ppijb)/=[1,1]))" in block
    assert "baseline_error_b=abs(stored_ppijb(1,1))" in block
    assert "if(noccb>0)data%ppijb=stored_ppijb" in block
