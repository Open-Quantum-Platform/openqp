import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_hessian_prepare.F90"
RESPONSE = ROOT / "source/modules/tdhf_hessian_response.F90"


def test_prepare_layer_connects_complete_first_response():
    text = SOURCE.read_text()
    for call in (
        "build_rohf_nuclear_response",
        "build_orbital_density_derivatives",
        "mrsf_xc_fock_total_derivative",
        "solve_mrsf_first_nuclear_response",
    ):
        assert call in text
    for result in ("dmo_a", "dmo_b", "dfock_a", "dfock_b", "dax", "dx"):
        assert result in text
    assert "require_two_somo=.true." in text


def test_prepare_layer_fails_closed_outside_initial_scope():
    compact = "".join(SOURCE.read_text().lower().split())
    for condition in (
        "infos%control%scftype/=3",
        "nocca-noccb/=2",
        "infos%tddft%umrsf",
        "infos%dft%cam_flag",
        "infos%functional%needtau",
        "infos%functional%needlapl",
    ):
        assert condition in compact
    assert "hiroya" in compact and "nakata" in compact


_REMOVED_FLOOR = re.compile(r"max\(1\.0e-8_dp,100\.0_dp\*(?:solve_tol|tolerance)\)")
_RAISED_FLOOR = re.compile(r"max\(1\.0e-6_dp,100\.0_dp\*(solve_tol|tolerance)\)")


def eigenpair_contract_problems(compact):
    """Return how a compacted tdhf_hessian_response.F90 breaks the contract.

    5f20bba4 (recovered chc4 production edits) removed the independent
    eigenpair re-adjudication gate -- the Davidson run already certified
    (omega, x0), and the extra absolute threshold rejected converged MRSF
    pairs whenever |omega| < 1 Eh -- and raised the response-certification
    floors from 1e-8 to 1e-6 relative.  What must still hold: the eigenpair
    residual is formed and reported on the single-vector and batch paths, an
    operator failure still aborts the single-vector response, no
    eigenpair_tol gate is left behind, no 1e-8 certification floor survives
    beside the raised 1e-6 ones on both the solve_tol and tolerance paths,
    and a stagnated Krylov solve is accepted only at the 1e-6 relative floor.
    """
    problems = []
    for needle, why in (
            ("residual_max=maxval(abs(ax-omega*x0))",
             "single-vector eigenpair residual is no longer formed"),
            ("residual_max=maxval(abs(applied(:,1)-omega*x0))",
             "batch eigenpair residual is no longer formed"),
            ("if(operator_status/=0)thenwrite(error_unit,'(a,i0)')"
             "'mrsfresponseoperatorfailedwithstatus',operator_status",
             "an operator failure no longer aborts the response"),
            ("stagnation_accept=1.0e-6_dp",
             "stagnation acceptance floor is not 1e-6"),
            ("if(stagnation_residual<=stagnation_accept)",
             "a stagnated solve is not decided on its residual")):
        if needle not in compact:
            problems.append(why)
    if "eigenpair_tol" in compact:
        problems.append("the eigenpair re-adjudication gate is back")
    if _REMOVED_FLOOR.search(compact):
        problems.append("a 1e-8 certification floor survived")
    if {m.group(1) for m in _RAISED_FLOOR.finditer(compact)} != {
            "solve_tol", "tolerance"}:
        problems.append("the raised 1e-6 floor is missing on a solver path")
    return problems


def test_eigenpair_criterion_respects_double_precision_and_solver_resolution():
    compact = "".join(RESPONSE.read_text().lower().split()).replace("&", "")
    assert eigenpair_contract_problems(compact) == []
