"""Native dependence decisions must depend on direction, not residual size."""
import ctypes
from pathlib import Path
import re


def test_native_davidson_candidate_scaling_and_dependence():
    from oqp import runtime

    root, suffix = runtime.resolve_oqp_root()
    library = ctypes.CDLL(runtime.library_path(root, suffix))
    check = library.davidson_orthogonalization_selftest
    check.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int),
                      ctypes.POINTER(ctypes.c_int)]
    check.restype = None
    error, checks, failures = ctypes.c_double(), ctypes.c_int(), ctypes.c_int()
    check(ctypes.byref(error), ctypes.byref(checks), ctypes.byref(failures))
    print(f"DAVIDSON_ORTHOGONALIZATION cases={checks.value} "
          f"failures={failures.value} max_abs={error.value:.16g}")
    assert checks.value == 61
    assert failures.value == 0
    assert 0 <= error.value <= 1e-12


def test_davidson_existing_basis_and_relative_dependence_source_contract():
    source = (Path(__file__).resolve().parents[1] / "source/tdhf_lib.F90").read_text()
    routine = re.search(r"(?is)\n  subroutine rpanewb\(.*?\n  end subroutine rpanewb", source).group()
    code = "\n".join(line.split("!")[0] for line in routine.lower().splitlines())
    # intent(out) makes the old basis undefined on entry, although MGS reads it.
    assert re.search(r"intent\(inout\).*::\s*bvec", code)
    # Keep the overflow-safe scaling before the first projection; testing only
    # the post-projection cutoff would miss the original amplitude dependence.
    finite = code.index("ieee_is_finite(q(:,k))")
    scale = code.index("maxval(abs(q(:,k)))")
    unit = code.index("q(:,k) = q(:,k)/fnorm")
    project = code.index("dot_product(bvec(:,istat),q(:,k))")
    assert finite < scale < unit < project


def test_sf_mrsf_drivers_mark_converged_roots_optional():
    # The native self-test covers rpanewb's optional-candidate ordering; this
    # pins that both SF-type drivers actually pass their converged roots.
    root = Path(__file__).resolve().parents[1] / "source/modules"
    for name, nroot in (("tdhf_sf_energy.F90", "nstates"), ("tdhf_mrsf_energy.F90", "nsolve")):
        code = re.sub(r"!.*", "", (root / name).read_text().lower())
        code = re.sub(r"&\s*\n\s*", "", code)
        assert re.search(rf"call rpanewb\({nroot},.*optional_candidate=rnorm\(1:{nroot}\)<=cnvtol\)", code), name
