import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class OpenTrustRegionLinalgConfigTests(unittest.TestCase):
    def test_linalg_none_is_rejected_before_external_build(self):
        top_cmake = (ROOT / "CMakeLists.txt").read_text()

        self.assertIn("LINALG_LIB=none is not supported", top_cmake)
        self.assertIn("OpenTrustRegion", top_cmake)
        self.assertIn("findLinearAlgebra()", top_cmake)
        self.assertLess(
            top_cmake.index("findLinearAlgebra()"),
            top_cmake.index("add_subdirectory(external)"),
        )

    def test_opentrustregion_receives_resolved_linalg_configuration(self):
        external_cmake = (ROOT / "external" / "CMakeLists.txt").read_text()

        self.assertIn("OTR_CMAKE_ARGS", external_cmake)
        self.assertIn("INTEGER_SIZE=${BLA_SIZEOF_INTEGER}", external_cmake)
        self.assertIn("OTR_DEFS", external_cmake)
        self.assertIn("otr_integer_flags", external_cmake)
        self.assertIn("fdefault-integer-8", external_cmake)
        self.assertIn("fallow-argument-mismatch", external_cmake)
        self.assertIn("CMAKE_Fortran_FLAGS=${otr_fortran_flags}", external_cmake)
        self.assertIn("OQP_EXTERNAL_LIST_SEPARATOR", external_cmake)
        self.assertIn("oqp_external_cmake_list_arg(_otr_blas_arg BLAS_LIBRARIES", external_cmake)
        self.assertIn("oqp_external_cmake_list_arg(_otr_lapack_arg LAPACK_LIBRARIES", external_cmake)
        self.assertIn("LIST_SEPARATOR ${OQP_EXTERNAL_LIST_SEPARATOR}", external_cmake)
        self.assertIn("add_dependencies(libopentrustregion LAPACK)", external_cmake)
        self.assertIn("CMAKE_POLICY_VERSION_MINIMUM=3.5", external_cmake)
        self.assertIn("PATCH_COMMAND /usr/bin/perl", external_cmake)
        self.assertIn("cmake_minimum_required(VERSION 3.5)", external_cmake)
        self.assertLess(
            external_cmake.index("ExternalProject_Add(LAPACK"),
            external_cmake.index("ExternalProject_Add(libopentrustregion"),
        )

    def test_find_linalg_can_run_before_oqp_target_exists(self):
        functions_cmake = (ROOT / "cmake" / "oqp_functions.cmake").read_text()

        self.assertIn("if(TARGET oqp)", functions_cmake)
        self.assertIn("add_dependencies(oqp LAPACK)", functions_cmake)
        self.assertIn("target_link_libraries(oqp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})", functions_cmake)

    def test_python_wheel_uses_native_trah_without_building_opentrah(self):
        pyproject = (ROOT / "pyproject.toml").read_text()

        self.assertIn('ENABLE_OPENTRAH = "OFF"', pyproject)


if __name__ == "__main__":
    unittest.main()
