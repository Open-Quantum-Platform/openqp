from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_lp64_configure_path_is_forbidden():
    cmake = (ROOT / "CMakeLists.txt").read_text()

    assert "if(NOT LINALG_LIB_INT64)" in cmake
    assert "LP64 BLAS/LAPACK is no longer supported" in cmake
    assert "OpenQP requires ILP64 BLAS/LAPACK with OQP_BLAS_INT=8" in cmake


def test_blas_lapack_discovery_requests_ilp64():
    top = (ROOT / "CMakeLists.txt").read_text()
    functions = (ROOT / "cmake" / "oqp_functions.cmake").read_text()

    assert "set(BLA_SIZEOF_INTEGER 8)" in top
    assert "set(BLA_SIZEOF_INTEGER 8)" in functions
    assert "Selected BLAS reports" in functions
    assert "Selected LAPACK reports" in functions


def test_oqp_and_wrapper_sources_receive_ilp64_define():
    source_cmake = (ROOT / "source" / "CMakeLists.txt").read_text()
    mathlib_cmake = (ROOT / "source" / "mathlib" / "CMakeLists.txt").read_text()

    assert "target_compile_definitions(oqp PRIVATE OQP_BLAS_INT=8)" in source_cmake
    assert "SOURCES_mathlib_ilp64" in source_cmake
    assert "OQP_BLAS_INT=8" in source_cmake
    assert "blas_wrap.F90" in mathlib_cmake
    assert "lapack_wrap.F90" in mathlib_cmake


def test_removed_dlfind_files_and_dependencies_do_not_reappear():
    assert not (ROOT / "pyoqp" / "oqp" / "library" / "libdlfind.py").exists()
    assert not (ROOT / "pyoqp" / "patch.sh").exists()

    pyproject = (ROOT / "pyproject.toml").read_text().lower()
    setup_py = (ROOT / "pyoqp" / "setup.py").read_text().lower()
    requirements = (ROOT / "pyoqp" / "requirements.txt").read_text().lower()

    assert "libdlfind" not in pyproject
    assert "libdlfind" not in setup_py
    assert "libdlfind" not in requirements
