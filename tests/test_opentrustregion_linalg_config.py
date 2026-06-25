import shutil
import subprocess
import tempfile
import textwrap
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

    def test_lp64_blas_is_macos_only(self):
        top_cmake = (ROOT / "CMakeLists.txt").read_text()

        self.assertIn("if(NOT LINALG_LIB_INT64 AND NOT APPLE)", top_cmake)
        self.assertIn("LP64 BLAS/LAPACK is supported only on macOS", top_cmake)

    def test_default_integer8_flags_are_gated_on_linalg_int64(self):
        top_cmake = (ROOT / "CMakeLists.txt").read_text()

        guard = top_cmake.index("if(LINALG_LIB_INT64)")
        before_guard = top_cmake[:guard]
        int64_guard_block = top_cmake[guard:top_cmake.index("if(ENABLE_OPENMP)")]

        self.assertNotIn("fdefault-integer-8", before_guard)
        self.assertNotIn(":-i8>", before_guard)
        self.assertIn("fdefault-integer-8", int64_guard_block)
        self.assertIn(":-i8>", int64_guard_block)
        self.assertIn("set(BLA_SIZEOF_INTEGER 4)", int64_guard_block)

    def test_oqp_blas_integer_width_follows_linalg_configuration(self):
        source_cmake = (ROOT / "source" / "CMakeLists.txt").read_text()
        functions_cmake = (ROOT / "cmake" / "oqp_functions.cmake").read_text()

        self.assertIn("OQP_BLAS_INT=${BLA_SIZEOF_INTEGER}", source_cmake)
        self.assertIn("OpenQP configured for ${BLA_SIZEOF_INTEGER}-byte BLAS integers", functions_cmake)

    def test_atomic_structure_init_uses_explicit_atom_count_kind(self):
        atomic_structure = (ROOT / "source" / "atomic_structure.F90").read_text()

        self.assertIn("c_int64_t", atomic_structure)
        self.assertIn("integer(c_int64_t), intent(in) :: natoms", atomic_structure)

    def test_libtagarray_accepts_default_integer_reserve_sizes(self):
        external_cmake = (ROOT / "external" / "CMakeLists.txt").read_text()
        patch = (ROOT / "cmake" / "patches" / "libtagarray-v0.0.6-default-integer-shapes.patch").read_text()

        self.assertIn("libtagarray-v0.0.6-default-integer-shapes.patch", external_cmake)
        self.assertIn("grep -q reserve_data_default", external_cmake)
        self.assertIn("generic, public    :: reserve_data => reserve_data_i64, reserve_data_default", patch)
        self.assertIn("integer,                              optional, intent(in) :: array_shape(:)", patch)
        self.assertIn("int(array_size, c_int64_t)", patch)

    def test_external_projects_receive_top_level_compilers(self):
        external_cmake = (ROOT / "external" / "CMakeLists.txt").read_text()

        libecpint_block = external_cmake[
            external_cmake.index("ExternalProject_Add(libecpint"):
            external_cmake.index("if(_LINALG_LIB_TYPE STREQUAL NetLib)")
        ]

        self.assertIn("CMAKE_C_COMPILER=${CMAKE_C_COMPILER}", libecpint_block)
        self.assertIn("CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}", libecpint_block)
        self.assertIn("CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}", libecpint_block)

    def test_external_projects_receive_make_program_when_set(self):
        external_cmake = (ROOT / "external" / "CMakeLists.txt").read_text()

        self.assertIn("OQP_EXTERNAL_GENERATOR_ARGS", external_cmake)
        self.assertIn("CMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}", external_cmake)
        self.assertIn('string(MD5 _oqp_external_make_program_hash "${CMAKE_MAKE_PROGRAM}")', external_cmake)
        self.assertIn('string(APPEND _oqp_external_generator "-make${_oqp_external_make_program_hash}")', external_cmake)
        self.assertIn("${OQP_EXTERNAL_GENERATOR_ARGS}", external_cmake)

    def test_reusable_external_cache_roots_are_created_before_downloads(self):
        external_cmake = (ROOT / "external" / "CMakeLists.txt").read_text()

        self.assertIn("file(MAKE_DIRECTORY", external_cmake)
        self.assertIn('"${OQP_EXTERNALS_SOURCE_ROOT}"', external_cmake)
        self.assertIn('"${OQP_EXTERNALS_BUILD_ROOT}"', external_cmake)
        self.assertIn('"${OQP_EXTERNALS_INSTALL_ROOT}"', external_cmake)

    def test_gradient_atom_count_uses_explicit_kind(self):
        grd1 = (ROOT / "source" / "integrals" / "grd1.F90").read_text()

        self.assertIn("use iso_c_binding, only: c_int64_t", grd1)
        self.assertIn("integer(c_int64_t), intent(in) :: n", grd1)

    def test_dft_grid_passes_default_integer_selectors_to_legacy_helpers(self):
        dft = (ROOT / "source" / "dftlib" / "dft.F90").read_text()

        self.assertIn("rad_grid_type = int(infos%dft%rad_grid_type)", dft)
        self.assertIn("dft_partfun = int(infos%dft%dft_partfun)", dft)
        self.assertIn("dft_bfc_algo = int(infos%dft%dft_bfc_algo)", dft)
        self.assertIn("nrad = int(infos%dft%grid_rad_size)", dft)
        self.assertIn("partFunType=dft_partfun", dft)

    def test_basis_overlap_print_helper_receives_default_integer_count(self):
        get_basis_overlap = (ROOT / "source" / "modules" / "get_basis_overlap.F90").read_text()

        self.assertIn("int(infos%mol_prop%nelec_a)", get_basis_overlap)

    def test_hf_hessian_passes_default_integer_scf_selector_to_dft(self):
        hf_hessian = (ROOT / "source" / "modules" / "hf_hessian.F90").read_text()

        self.assertIn("int(infos%control%scftype)", hf_hessian)

    def test_mrsf_gmres_controls_are_default_integer_solver_inputs(self):
        mrsf_z = (ROOT / "source" / "modules" / "tdhf_mrsf_z_vector.F90").read_text()

        self.assertIn("restart = min(int(infos%tddft%gmres_dim), lzdim)", mrsf_z)
        self.assertIn("max_iter = int(infos%control%maxit_zv)", mrsf_z)

    def test_scf_converger_controls_are_default_integer_inputs(self):
        scf = (ROOT / "source" / "scf.F90").read_text()

        self.assertIn("control_converger = int(infos%control%converger_type)", scf)
        self.assertIn("control_diis = int(infos%control%diis_type)", scf)
        self.assertIn("control_verbose = int(infos%control%verbose)", scf)
        self.assertIn("control_maxit = int(infos%control%maxit)", scf)
        self.assertIn("control_scftype = int(infos%control%scftype)", scf)

    def test_openblas_config_fallback_is_forwarded_to_subbuilds(self):
        functions_cmake = (ROOT / "cmake" / "oqp_functions.cmake").read_text()
        external_cmake = (ROOT / "external" / "CMakeLists.txt").read_text()

        self.assertIn("find_package(OpenBLAS CONFIG QUIET)", functions_cmake)
        self.assertIn("OPENBLAS_USE64BITINT", functions_cmake)
        self.assertIn("OpenBLAS_LIBRARY", functions_cmake)
        self.assertIn("IMPORTED_LOCATION", functions_cmake)
        self.assertIn("set(BLAS_LIBRARIES ${_openblas_libraries})", functions_cmake)
        self.assertIn("set(LAPACK_LIBRARIES ${_openblas_libraries})", functions_cmake)
        self.assertIn("oqp_external_cmake_list_arg(_otr_blas_arg BLAS_LIBRARIES ${BLAS_LIBRARIES})", external_cmake)
        self.assertIn("oqp_external_cmake_list_arg(_otr_lapack_arg LAPACK_LIBRARIES ${LAPACK_LIBRARIES})", external_cmake)
        self.assertIn("oqp_external_cmake_list_arg(_ddx_blas_arg BLAS_LIBRARIES ${BLAS_LIBRARIES})", external_cmake)
        self.assertIn("oqp_external_cmake_list_arg(_ddx_lapack_arg LAPACK_LIBRARIES ${LAPACK_LIBRARIES})", external_cmake)

    def test_openblas_config_fallback_uses_imported_library_before_transitive_deps(self):
        cmake = shutil.which("cmake")
        self.assertIsNotNone(cmake)

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            package_root = tmp_path / "openblas"
            include_dir = package_root / "include"
            lib_dir = package_root / "lib"
            config_dir = lib_dir / "cmake" / "OpenBLAS"
            source_dir = tmp_path / "source"
            build_dir = tmp_path / "build"
            include_dir.mkdir(parents=True)
            lib_dir.mkdir(parents=True)
            config_dir.mkdir(parents=True)
            source_dir.mkdir()

            openblas_library = lib_dir / "libopenblas.a"
            openblas_library.write_text("")
            (include_dir / "openblas_config.h").write_text("#define OPENBLAS_USE64BITINT 1\n")
            (config_dir / "OpenBLASConfig.cmake").write_text(textwrap.dedent(f"""
                set(OpenBLAS_FOUND TRUE)
                set(OpenBLAS_INCLUDE_DIRS [=[{include_dir.as_posix()}]=])
                set(OpenBLAS_LIBRARY [=[{openblas_library.as_posix()}]=])
                set(OpenBLAS_LIBRARIES pthread;m)
                add_library(OpenBLAS::OpenBLAS UNKNOWN IMPORTED)
                set_target_properties(OpenBLAS::OpenBLAS PROPERTIES
                  IMPORTED_LOCATION [=[{openblas_library.as_posix()}]=]
                  INTERFACE_LINK_LIBRARIES "${{OpenBLAS_LIBRARIES}}"
                  INTERFACE_INCLUDE_DIRECTORIES "${{OpenBLAS_INCLUDE_DIRS}}")
            """))
            (source_dir / "CMakeLists.txt").write_text(textwrap.dedent(f"""
                cmake_minimum_required(VERSION 3.16)
                project(openblas_config_probe NONE)
                set(CMAKE_PREFIX_PATH [=[{package_root.as_posix()}]=])
                set(BLA_SIZEOF_INTEGER 8)
                set(LINALG_LIB_INT64 ON)
                include([=[{(ROOT / "cmake" / "oqp_functions.cmake").as_posix()}]=])
                findOpenBLASConfig()
                if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
                  message(FATAL_ERROR "OpenBLAS config fallback did not mark BLAS/LAPACK found")
                endif()
                list(GET BLAS_LIBRARIES 0 _blas_first)
                if(NOT _blas_first STREQUAL [=[{openblas_library.as_posix()}]=])
                  message(FATAL_ERROR "BLAS_LIBRARIES did not start with the OpenBLAS library: ${{BLAS_LIBRARIES}}")
                endif()
                if(BLAS_LIBRARIES MATCHES "OpenBLAS::OpenBLAS")
                  message(FATAL_ERROR "BLAS_LIBRARIES leaked an imported target into subbuild args: ${{BLAS_LIBRARIES}}")
                endif()
                if(NOT BLAS_LIBRARIES MATCHES "pthread" OR NOT BLAS_LIBRARIES MATCHES "m")
                  message(FATAL_ERROR "BLAS_LIBRARIES lost transitive dependencies: ${{BLAS_LIBRARIES}}")
                endif()
                if(NOT BLAS_SIZEOF_INTEGER EQUAL 8 OR NOT LAPACK_SIZEOF_INTEGER EQUAL 8)
                  message(FATAL_ERROR "OpenBLAS config fallback did not preserve the requested integer width")
                endif()
            """))

            result = subprocess.run(
                [cmake, "-S", str(source_dir), "-B", str(build_dir)],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
