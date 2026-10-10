from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_native_trah_allocates_rhf_beta_shadow_before_copying_mos():
    source = (ROOT / "source" / "trah_converger.F90").read_text()

    guard = "if (.not. allocated(conv%mo_b)) allocate(conv%mo_b, source=conv%mo_a)"
    copy = "allocate(mo0_a, source=conv%mo_a); allocate(mo0_b, source=conv%mo_b)"

    assert guard in source
    assert source.index(guard) < source.index(copy)
