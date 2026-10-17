"""Native pin of UMRSF bra-image reuse, including empty orbital spaces."""
import numpy as np


def test_umrsf_bra_images_match_independent_contractions():
    import oqp

    error = oqp.ffi.new("double *")
    checks = oqp.ffi.new("int *")
    # Missing entry points must fail instead of skipping a stale build.
    oqp.lib.umrsf_bra_reuse_selftest(error, checks)
    assert checks[0] == 1056
    assert np.isfinite(error[0]) and error[0] < 2e-13
