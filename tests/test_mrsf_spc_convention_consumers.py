from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONSUMERS = (
    ROOT / "source/modules/tdhf_mrsf_energy.F90",
    ROOT / "source/modules/tdhf_mrsf_gradient.F90",
    ROOT / "source/modules/tdhf_mrsf_z_vector.F90",
)


def test_all_mrsf_derivative_consumers_use_the_named_raw_spc_sign():
    for path in CONSUMERS:
        source = path.read_text().lower()
        assert "use tdhf_mrsf_conventions_mod" in source, path
        assert "mrsf_raw_spc_multiplier" in source, path


def test_raw_spc_sign_is_not_reintroduced_as_a_triplet_literal():
    for path in CONSUMERS:
        source = "".join(path.read_text().lower().split())
        assert "if(mrst==3)sgnk=-1" not in source, path
        assert "if(mrst==3)fmrst2(:,1:6,:,:)=-" not in source, path
