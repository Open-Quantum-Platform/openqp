"""Regression tests for the reusable full-space ROHF nuclear response."""

from __future__ import annotations

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESPONSE = ROOT / "source/modules/hf_rohf_orbital_response.F90"
HESSIAN = ROOT / "source/modules/hf_hessian.F90"


def _compact(path: Path) -> str:
    return "".join(path.read_text().lower().split())


def _complete_connection(
    mo: np.ndarray,
    sx_mo: np.ndarray,
    nocc: int,
    xvo: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    connection = -0.5 * sx_mo
    connection[nocc:, :nocc] = xvo
    connection[:nocc, nocc:] = -sx_mo[:nocc, nocc:] - xvo.T
    return mo @ connection, connection


def test_production_response_has_full_alpha_beta_orbital_connections():
    source = _compact(RESPONSE)
    assert "response%dmo_alpha(nbf,nbf,ncart)" in source
    assert "response%dmo_beta(nbf,nbf,ncart)" in source
    assert "response%ds_ao(nbf,nbf,ncart)" in source
    assert "response%dhcore_ao(nbf,nbf,ncart)" in source
    assert "connection(nocc+1:nbf,1:nocc)=xvo" in source
    assert (
        "connection(1:nocc,nocc+1:nbf)=&"
        "-sx_mo(1:nocc,nocc+1:nbf)-transpose(xvo)"
    ) in source
    assert "dmo=matmul(mo,connection)" in source


def test_two_somo_default_and_cam_are_fail_closed():
    source = _compact(RESPONSE)
    assert "two_somo_only=.true." in source
    assert "if(two_somo_only.and.offset/=2)thenstatus=-4" in source
    assert "if(dft.and.infos%dft%cam_flag)thenstatus=-3" in source
    assert "if(infos%control%hamilton>20)thenstatus=-7" in source


def test_hf_hessian_reuses_response_without_a_second_rohf_cphf_prepass():
    source = HESSIAN.read_text().lower()
    rohf = source.split("subroutine hf_hessian_rohf", 1)[1].split(
        "end subroutine hf_hessian_rohf", 1
    )[0]
    compact = "".join(rohf.split())
    assert "callbuild_rohf_nuclear_response" in compact
    assert "require_two_somo=.false." in compact
    assert "callcphf_solve_rohf" not in compact
    assert "callrohf_pack_trial" not in compact
    assert "callrohf_unpack_trial" not in compact
    assert "orbital_response%dmo_alpha(:,1:nocca,x)" in compact
    assert "orbital_response%dmo_beta(:,1:noccb,x)" in compact


def test_origin_and_scope_are_explicit_without_forbidden_representation():
    source = RESPONSE.read_text().lower()
    assert "hiroya nakata" in source
    assert "total\n    ! df/dr is deliberately not returned" in source
    for forbidden in ("slater", "determinant", "fock_space"):
        assert forbidden not in source


def test_full_connections_obey_the_moving_metric_for_both_spin_partitions():
    rng = np.random.default_rng(81273)
    nbf = 7
    nocca = 4
    noccb = 2

    raw = rng.normal(size=(nbf, nbf))
    overlap = raw.T @ raw + 2.0 * np.eye(nbf)
    eigval, eigvec = np.linalg.eigh(overlap)
    invsqrt = (eigvec / np.sqrt(eigval)) @ eigvec.T
    q, _ = np.linalg.qr(rng.normal(size=(nbf, nbf)))
    mo = invsqrt @ q
    np.testing.assert_allclose(mo.T @ overlap @ mo, np.eye(nbf), atol=2.0e-14)

    raw_ds = rng.normal(size=(nbf, nbf))
    ds_ao = 0.5 * (raw_ds + raw_ds.T)
    sx_mo = mo.T @ ds_ao @ mo
    xa = rng.normal(size=(nbf - nocca, nocca))
    xb = rng.normal(size=(nbf - noccb, noccb))

    dmo_a, ua = _complete_connection(mo, sx_mo, nocca, xa)
    dmo_b, ub = _complete_connection(mo, sx_mo, noccb, xb)
    for dmo, connection in ((dmo_a, ua), (dmo_b, ub)):
        metric_derivative = (
            dmo.T @ overlap @ mo + mo.T @ overlap @ dmo + sx_mo
        )
        np.testing.assert_allclose(metric_derivative, 0.0, atol=3.0e-14)
        np.testing.assert_allclose(mo.T @ overlap @ dmo, connection, atol=2.0e-14)

    np.testing.assert_allclose(ua[nocca:, :nocca], xa, atol=0.0)
    np.testing.assert_allclose(ub[noccb:, :noccb], xb, atol=0.0)
    assert not np.allclose(dmo_a[:, noccb:nocca], dmo_b[:, noccb:nocca])
