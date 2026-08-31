"""Contracts for the differentiated Gxc term in the TDDFT Z-vector RHS."""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
Z_RHS = (ROOT / "source" / "modules" / "tdhf_hessian_z_rhs.F90").read_text().lower()
DENSITY = (ROOT / "source" / "modules" / "tdhf_hessian_density.F90").read_text().lower()
HELPER = (ROOT / "source" / "modules" / "tdhf_hessian_gxc_derivative.F90").read_text().lower()


def test_gxc_derivative_has_one_reusable_owner():
    assert "module tdhf_hessian_gxc_derivative_mod" in HELPER
    assert "public :: build_gxc_and_derivative" in HELPER
    assert "subroutine build_gxc_and_derivative" not in DENSITY
    assert "use tdhf_hessian_gxc_derivative_mod" in DENSITY
    assert "use tdhf_hessian_gxc_derivative_mod" in Z_RHS


def test_z_rhs_builds_gxc_derivative_from_xpy_and_dxpy():
    assert re.search(
        r"call\s+build_gxc_and_derivative\s*\(\s*infos\s*,\s*mo\s*,\s*umat\s*,"
        r"\s*um\s*,\s*du\s*,\s*gxp\s*,\s*dgxp\s*\)",
        Z_RHS,
    )
    assert "if (infos%control%hamilton == 20)" in Z_RHS


def test_z_rhs_contains_exact_missing_minus_two_dgxc_ov_term():
    assert re.search(
        r"drhs\(:,k\)\s*=\s*drhs\(:,k\)\s*&\s*"
        r"-\s*2\.0_dp\s*\*\s*reshape\s*\(\s*dgxp\(1:nocc,nocc\+1:,k\)\s*,"
        r"\s*\[nexc\]\s*\)",
        Z_RHS,
        re.DOTALL,
    )


def test_helper_restores_geometry_after_each_central_difference():
    assert "basis%atoms%xyz(cc,aa)=basis%atoms%xyz(cc,aa)+step" in HELPER
    assert "basis%atoms%xyz(cc,aa)=basis%atoms%xyz(cc,aa)-2.0_dp*step" in HELPER
    assert HELPER.count("call basis%init_shell_centers()") >= 3


def test_every_gxc_accumulator_is_zeroed_immediately_before_use():
    calls = [
        "call tddft_gxc(basis,grid,.true.,c,gp,xao,1,0.0_dp,infos)",
        "call tddft_gxc(basis,grid,.true.,cp,gp,xao,1,0.0_dp,infos)",
        "call tddft_gxc(basis,grid,.true.,cm,gm,xao,1,0.0_dp,infos)",
    ]
    lines = [line.strip() for line in HELPER.splitlines()]
    for call in calls:
        index = lines.index(call)
        accumulator = "gm=0.0_dp" if ",cm,gm," in call else "gp=0.0_dp"
        assert lines[index-1] == accumulator
