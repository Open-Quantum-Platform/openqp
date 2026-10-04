"""The libecpint source transformation in external/fix_libecpint_accuracy.py.

libecpint 1.0.7 truncates the generated type-2 angular coefficients to six
significant digits and stops its one-point Gauss-Chebyshev rule at the first
level that passes its convergence test.  Together they made the analytic ECP
gradient of SiO/LANL2DZ differ from finite differences by 1.2e-4 Ha/bohr
(examples/ECP/SiO_RHF-BHHLYP_ECP_GRADIENT.inp).  These tests run the patch on
verbatim excerpts of the upstream sources (trailing blanks included), so a
silent no-op or a partial application fails here rather than as a gradient
drift.
"""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "external" / "fix_libecpint_accuracy.py"

GENERATE_CPP = (
    "#include <algorithm>\n"
    "#include <set>\n"
    "\n"
    '#include "generate.hpp"\n'
    "\t// Create the code file\n"
    '\tstd::string ofname = "generated/Q" + std::to_string(LA) + '
    'std::to_string(LB) + std::to_string(lam) + ".cpp"; \n'
    "\tstd::ofstream outfile(ofname); \n"
    "\t\n"
    "\tif (!outfile.is_open())\n"
)

GAUSSQUAD_CPP = (
    "\t\t\tint p = (M+1) / 2; // M / 2^n \n"
    "\t\t\twhile (n < maxN && !converged) {\n"
    "\t\t\t\t// Compute T_{2n+1}\n"
    "\t\t\t\tT2n1 = Tn + sumTerms(f, params, n, start, end, p, 2);\n"
    "\t\t\t\n"
    "\t\t\t\t// Check convergence\n"
    "\t\t\t\tdT = T2n1 - 2.0*Tn;\n"
    "\t\t\t\tn = 2*n + 1;\n"
    "\t\t\t\tif (dT*dT <= fabs(T2n1 - Tn12)*tolerance) {\n"
    "\t\t\t\t\tconverged = true;  \n"
    "\t\t\t\t} else {\n"
    "\t\t\t\t\tTn12 = 4.0 * Tn; \n"
    "\t\t\t\t\tTn = T2n1;\n"
    "\t\t\t\t\tp /= 2; \n"
    "\t\t\t\t}\n"
    "\t\t\t}\n"
    "\t\t\t// Finalise the integral\n"
    "\t\t\tI = 16.0 * T2n1 / (3.0*(n + 1.0));\n"
)


def _load():
    spec = importlib.util.spec_from_file_location("fix_libecpint_accuracy", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _tree(tmp_path, generate=GENERATE_CPP, gaussquad=GAUSSQUAD_CPP):
    (tmp_path / "src" / "lib").mkdir(parents=True)
    (tmp_path / "src" / "generate.cpp").write_text(generate)
    (tmp_path / "src" / "lib" / "gaussquad.cpp").write_text(gaussquad)
    return tmp_path


def test_patch_applies_once_and_is_idempotent(tmp_path):
    fix = _load()
    root = _tree(tmp_path)
    fix.patch(root)
    gen = (root / "src" / "generate.cpp").read_text()
    quad = (root / "src" / "lib" / "gaussquad.cpp").read_text()

    assert "#include <iomanip>" in gen and "#include <limits>" in gen
    assert gen.count("std::setprecision(std::numeric_limits<double>::max_digits10)") == 1
    # the precision is set on the stream before anything is written to it
    assert gen.index("max_digits10") < gen.index("if (!outfile.is_open())")

    assert "!converged" not in quad
    assert "while (n < maxN) {" in quad
    assert "if (dT*dT <= fabs(T2n1 - Tn12)*tolerance) converged = true;" in quad
    # the refinement advances unconditionally, so the full grid is summed
    loop = quad[quad.index("while (n < maxN)"):quad.index("// Finalise")]
    assert "else" not in loop
    for update in ("Tn12 = 4.0 * Tn;", "Tn = T2n1;", "p /= 2;"):
        assert update in loop

    fix.patch(root)
    assert (root / "src" / "generate.cpp").read_text() == gen
    assert (root / "src" / "lib" / "gaussquad.cpp").read_text() == quad


def test_unrecognized_source_fails_loudly(tmp_path):
    fix = _load()
    root = _tree(tmp_path, gaussquad=GAUSSQUAD_CPP.replace("p /= 2; ", "p >>= 1;"))
    with pytest.raises(RuntimeError, match="Unrecognized libecpint source"):
        fix.patch(root)


def test_cmake_applies_patch_and_keys_the_libecpint_cache_on_it():
    external_cmake = (ROOT / "external" / "CMakeLists.txt").read_text()
    block = external_cmake[
        external_cmake.index("ExternalProject_Add(libecpint"):
        external_cmake.index("if(_LINALG_LIB_TYPE STREQUAL NetLib)")
    ]
    assert "fix_libecpint_accuracy.py ${LIBECPINT_SOURCE_DIR}" in block
    assert "PATCH_COMMAND" in block
    # a cached unpatched libecpint must never be reused
    assert 'file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/fix_libecpint_accuracy.py"' in external_cmake
    assert "oqp_set_external_paths(LIBECPINT libecpint-p${_OQP_LIBECPINT_PATCH_KEY}" in external_cmake
