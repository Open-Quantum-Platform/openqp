from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_z_vector.F90"


def test_mrsf_zvector_preserves_multiplier_for_analytic_hessian():
    source = SOURCE.read_text()
    assert "alloc_or_die(OQP_td_z, (/ lzdim /), td_z" in source
    assert "description=OQP_td_z_comment" in source
    assert "td_z=xk" in source
    assignment = source.index("td_z=xk")
    solve = source.index("select case (infos%tddft%z_solver)", source.index("Step 2:"))
    relaxed = source.index("call build_mrsf_relaxed_density_and_w()")
    assert solve < assignment < relaxed
