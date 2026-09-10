"""Regression tests for the reusable full-space ROHF nuclear response."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESPONSE = ROOT / "source/modules/hf_rohf_orbital_response.F90"
HESSIAN = ROOT / "source/modules/hf_hessian.F90"
CPHF = ROOT / "source/modules/cphf.F90"


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


def test_roks_rhs_uses_all_coordinate_analytic_xc_derivative():
    source = _compact(RESPONSE)
    assert "callmrsf_xc_fock_total_derivative" in source
    assert "d0a_all(nbf,nbf,ncart)" in source
    assert "d0b_all(nbf,nbf,ncart)" in source
    assert "add_semilocal_roks_rhs" not in source
    assert "basis%atoms%xyz(cc,kc)=" not in source
    assert "step=1.0d-3" not in source


def test_hf_hessian_reuses_response_without_a_second_rohf_cphf_prepass():
    source = HESSIAN.read_text().lower()
    rohf = _routine(source, "hf_hessian_rohf")
    compact = "".join(rohf.split())
    assert "callbuild_rohf_nuclear_response" in compact
    assert "require_two_somo=.false." in compact
    assert "callcphf_solve_rohf" not in compact
    assert "callrohf_pack_trial" not in compact
    assert "callrohf_unpack_trial" not in compact
    assert "orbital_response%dmo_alpha(:,1:nocca,x)" in compact
    assert "orbital_response%dmo_beta(:,1:noccb,x)" in compact


def _routine(source, name):
    """Body of `subroutine <name>(infos)` in lower-cased source, exact name only."""
    match = re.search(
        rf"subroutine {name}\(infos\)(.*?)end subroutine {name}[ \t]*\n",
        source, re.S)
    return match.group(1) if match else ""


def cam_fallback_problems(source):
    """Return how lower-cased hf_hessian.F90 breaks the CAM fallback contract.

    The analytic ROHF nuclear response declines CAM/range-separated CPKS with
    status -3.  origin/main runs those ROKS Hessians through its
    semi-numerical response, so hf_hessian_rohf must hand over to that routine
    instead of aborting -- and only for status -3: every other response
    status still aborts, the handover returns immediately, it happens in
    exactly one place, and the analytic routine never calls the legacy CPHF
    solver itself.  The fallback must pass status= to cphf_solve_rohf and
    abort when the solve fails: without status the solver only logs the
    failure and returns partial amplitudes.
    """
    rohf = "".join(_routine(source, "hf_hessian_rohf").split())
    legacy = "".join(_routine(source, "hf_hessian_rohf_semi_numerical").split())
    problems = []
    if not rohf:
        return ["hf_hessian_rohf not found"]
    if not legacy:
        problems.append("the semi-numerical fallback routine is missing")
    elif "callcphf_solve_rohf" not in legacy:
        problems.append("the fallback routine no longer solves the ROHF CPHF")
    elif "callcphf_solve_rohf(infos,ncart,bvec,uvec,status=cphf_status)" not in legacy:
        problems.append("the fallback ignores the CPHF solver status")
    else:
        # Bound the check to this one call: in whitespace-free text a regex
        # would run on to some later WITH_ABORT in the legacy routine.
        parts = legacy.split("if(cphf_status/=0)callshow_message(", 1)
        call = parts[1].split("endblock", 1)[0] if len(parts) == 2 else ""
        if not call.endswith("with_abort)"):
            problems.append("a failed fallback CPHF solve does not abort")
    guard = "if(response_status==-3)then"
    if guard not in rohf:
        problems.append("no fallback guarded on response_status==-3")
    else:
        block = rohf.split(guard, 1)[1].split("endif", 1)[0]
        if "callhf_hessian_rohf_semi_numerical(infos)" not in block:
            problems.append("the status -3 branch does not call the fallback")
        if "return" not in block:
            problems.append("the status -3 branch does not return after the fallback")
    if rohf.count("callhf_hessian_rohf_semi_numerical(") != 1:
        problems.append("the fallback is called outside the status -3 branch")
    if "if(response_status/=0)callshow_message(" not in rohf:
        problems.append("other response failures no longer abort")
    if "callcphf_solve_rohf" in rohf:
        problems.append("the analytic routine calls the legacy CPHF solver")
    return problems


def test_cam_rohf_hessian_falls_back_only_when_the_analytic_response_declines():
    assert cam_fallback_problems(HESSIAN.read_text().lower()) == []


def test_rohf_hessian_response_fails_closed_on_any_cphf_rhs_failure():
    response = _compact(RESPONSE)
    solver = _compact(CPHF)
    assert "tol=1.0d-13,status=status" in response
    assert "if(status/=0)thencallresponse%clean()returnendif" in response
    assert "integer,intent(out),optional::status" in solver
    assert "if(any(failed).or.any(active).or.&" in solver
    assert "any(.not.ieee_is_finite(uvec)))status=-1" in solver


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
