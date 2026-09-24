import configparser
import importlib.util
import re
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENERGY = ROOT / "source" / "modules" / "tdhf_mrsf_energy.F90"
UMRSF_GRAD = ROOT / "source" / "modules" / "tdhf_umrsf_gradient.F90"
UMRSF_ZVEC = ROOT / "source" / "modules" / "tdhf_umrsf_z_vector.F90"
LIB = ROOT / "source" / "tdhf_mrsf_lib.F90"
LIBXC = ROOT / "source" / "dftlib" / "libxc.F90"
TAGARRAY = ROOT / "source" / "tagarray_driver.F90"
SINGLE_POINT = ROOT / "pyoqp" / "oqp" / "library" / "single_point.py"
OQPDATA = ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"
INPUT_CHECKER = ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py"
MOLECULE = ROOT / "pyoqp" / "oqp" / "molecule" / "molecule.py"
UMRSF_BHHLYP_FD = ROOT / "tests" / "data" / "umrsf" / "H2CO_BHHLYP_UMRSF_GRADIENT.inp"
UMRSF_BLYP_FD = ROOT / "tests" / "data" / "umrsf" / "H2CO_BLYP_UMRSF_GRADIENT.inp"

# UMRSF supports energy, gradients, and gradient-driven searches. Other
# runtypes still need Hessians, NACs, or optimization-level reuse that has not
# been implemented for UMRSF yet.
UMRSF_BLOCKED_RUNTYPES = (
    "prop", "data", "hess", "nac", "nacme",
    "mep", "ts", "irc", "neb",
)


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text.lower())


def _fortran_subroutine(path, name):
    source = compact(path.read_text())
    start = source.index(f"subroutine{name}(")
    return source[start:source.index(f"endsubroutine{name}", start)]


def _libxc_aliases_by_feature():
    source = LIBXC.read_text().split("select case (funcname)", 1)[1].split(
        "end select", 1
    )[0]
    heads = list(re.finditer(
        r"(?mi)^[ \t]*case(?:[ \t]*\(.*\)|[ \t]+default).*$", source
    ))
    features = {name: set() for name in ("cam", "meta", "dh", "spc")}
    for index, head in enumerate(heads):
        aliases = {
            value.lower()
            for value in re.findall(
                r'"([^"]+)"', head.group(0).split("!", 1)[0]
            )
        }
        stop = heads[index + 1].start() if index + 1 < len(heads) else len(source)
        body = "\n".join(
            line.split("!", 1)[0] for line in source[head.end():stop].splitlines()
        )
        if re.search(r"dft_params%cam_flag\s*=\s*\.true\.", body, re.I):
            features["cam"] |= aliases
        if re.search(
                r"add_functional\s*\(\s*XC_(?:HYB_)?MGGA_", body, re.I):
            features["meta"] |= aliases
        if re.search(r"dft_params%dh_flag\s*=\s*\.true\.", body, re.I):
            features["dh"] |= aliases
        if re.search(
                r"tddft_params%spc_(?:coco|ovov|coov)\s*=", body, re.I):
            features["spc"] |= aliases
    return features


def _load_regression():
    """Load the regression registry without importing the ``oqp`` package.

    ``import oqp`` dlopens liboqp, which is not available in a source-only
    checkout; regression.py itself has no native dependency.
    """
    path = ROOT / "pyoqp" / "oqp" / "utils" / "regression.py"
    spec = importlib.util.spec_from_file_location("_oqp_regression", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_input_checker():
    """Import input_checker.py in isolation with a minimal MPI stub."""
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        size = 1
        use_mpi = False

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils

    spec = importlib.util.spec_from_file_location(
        "input_checker_umrsf_under_test", INPUT_CHECKER
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _umrsf_config(runtype):
    return {
        "input": {
            "runtype": runtype,
            "method": "tdhf",
            "functional": "bhhlyp",
            "basis": "6-31g",
            "system": "\nO 0.0 0.0 0.0\nH 0.0 0.0 0.95\nH 0.9 0.0 -0.3",
        },
        "guess": {},
        "scf": {"type": "uhf", "multiplicity": 3},
        "tdhf": {"type": "umrsf", "nstate": 2},
        "properties": {"grad": [1]},
        "optimize": {"lib": "geometric", "istate": 1, "jstate": 2},
        "nac": {"states": [[1, 2]]},
        "neb": {"product": "", "nimage": 3},
    }


def _umrsf_guard_errors(report):
    return [
        diag
        for diag in report.errors
        if diag.path == "tdhf.type"
        and "umrsf-tddft supports runtype" in diag.message.lower()
    ]


class UMRSFEnergyRegressionTests(unittest.TestCase):
    def test_umrsf_mixed_exchange_channels_use_gamess_compatible_permutation(self):
        source = compact(LIB.read_text())
        expected_updates = (
            "f3(:nf,9:10,i,l)=f3(:nf,9:10,i,l)-xval*d3(:nf,9:10,k,j)",
            "f3(:nf,9:10,l,i)=f3(:nf,9:10,l,i)-xval*d3(:nf,9:10,j,k)",
            "f3(:nf,9:10,k,j)=f3(:nf,9:10,k,j)-xval*d3(:nf,9:10,i,l)",
            "f3(:nf,9:10,j,k)=f3(:nf,9:10,j,k)-xval*d3(:nf,9:10,l,i)",
            "f3(:nf,9:10,i,k)=f3(:nf,9:10,i,k)-xval*d3(:nf,9:10,l,j)",
            "f3(:nf,9:10,k,i)=f3(:nf,9:10,k,i)-xval*d3(:nf,9:10,j,l)",
            "f3(:nf,9:10,l,j)=f3(:nf,9:10,l,j)-xval*d3(:nf,9:10,i,k)",
            "f3(:nf,9:10,j,l)=f3(:nf,9:10,j,l)-xval*d3(:nf,9:10,k,i)",
        )
        bad_head_updates = (
            "f3(:nf,9:10,i,k)=f3(:nf,9:10,i,k)-xval*d3(:nf,9:10,j,l)",
            "f3(:nf,9:10,k,i)=f3(:nf,9:10,k,i)-xval*d3(:nf,9:10,l,j)",
            "f3(:nf,9:10,i,l)=f3(:nf,9:10,i,l)-xval*d3(:nf,9:10,j,k)",
            "f3(:nf,9:10,l,i)=f3(:nf,9:10,l,i)-xval*d3(:nf,9:10,k,j)",
        )

        for update in expected_updates:
            self.assertIn(update, source)
        for update in bad_head_updates:
            self.assertNotIn(update, source)

    def test_umrsf_mixed_exchange_gradient_differentiates_the_energy_permutation(self):
        density = _fortran_subroutine(
            UMRSF_GRAD, "grd2_umrsf_resp_get_density"
        )
        branch = re.search(
            r"if\(this%transpose_exchange\(ch\)\)then(.*?)else(.*?)endif",
            density,
        )
        self.assertIsNotNone(branch)
        product_pattern = (
            r"(?<![%a-z_])bden\(ch,[^)]*\)\*dden\(ch,[^)]*\)"
        )
        expected_mixed_terms = (
            "bden(ch,i1,k1)*dden(ch,l1,j1)",
            "bden(ch,j1,k1)*dden(ch,l1,i1)",
            "bden(ch,i1,l1)*dden(ch,k1,j1)",
            "bden(ch,j1,l1)*dden(ch,k1,i1)",
            "bden(ch,k1,i1)*dden(ch,j1,l1)",
            "bden(ch,l1,i1)*dden(ch,j1,k1)",
            "bden(ch,k1,j1)*dden(ch,i1,l1)",
            "bden(ch,l1,j1)*dden(ch,i1,k1)",
        )
        expected_ordinary_terms = (
            "bden(ch,i1,k1)*dden(ch,j1,l1)",
            "bden(ch,j1,k1)*dden(ch,i1,l1)",
            "bden(ch,i1,l1)*dden(ch,j1,k1)",
            "bden(ch,j1,l1)*dden(ch,i1,k1)",
            "bden(ch,k1,i1)*dden(ch,l1,j1)",
            "bden(ch,l1,i1)*dden(ch,k1,j1)",
            "bden(ch,k1,j1)*dden(ch,l1,i1)",
            "bden(ch,l1,j1)*dden(ch,k1,i1)",
        )
        self.assertCountEqual(
            re.findall(product_pattern, branch.group(1)), expected_mixed_terms
        )
        self.assertCountEqual(
            re.findall(product_pattern, branch.group(2)), expected_ordinary_terms
        )
        self.assertIn("df1=df1-c*(", branch.group(1))
        self.assertIn("df1=df1-c*(", branch.group(2))

        fill = _fortran_subroutine(UMRSF_GRAD, "umrsf_resp_2pdm_fill")
        self.assertEqual(fill.count("gcomp%transpose_exchange=.false."), 1)
        self.assertEqual(fill.count("gcomp%transpose_exchange(9:10)=.true."), 1)
        split = _fortran_subroutine(
            UMRSF_GRAD, "umrsf_resp_2e_grad_split"
        )
        self.assertEqual(
            split.count(
                "one%transpose_exchange(1)=gcomp%transpose_exchange(ch)"
            ),
            1,
        )

    def test_umrsf_mixed_exchange_quartet_matches_k_of_transposed_density(self):
        # Use four distinct AO indices so the ERI's eight symmetry-related
        # permutations are all distinct.  The derivative of
        # -<B,K[D^T]> with respect to that unique integral is the coefficient
        # emitted by grd2_umrsf_resp_get_density for channels 9 and 10.
        i, j, k, l = range(4)
        permutations = (
            (i, j, k, l),
            (j, i, k, l),
            (i, j, l, k),
            (j, i, l, k),
            (k, l, i, j),
            (l, k, i, j),
            (k, l, j, i),
            (l, k, j, i),
        )
        bra = [[11.0 * (p + 1) + 3.0 * (q + 1) for q in range(4)]
               for p in range(4)]
        ket = [[7.0 * (p + 1) - 2.0 * (q + 1) for q in range(4)]
               for p in range(4)]

        mixed_energy_derivative = -sum(
            bra[p][r] * ket[s][q] for p, q, r, s in permutations
        )
        mixed_quartet_coefficient = -(
            bra[i][k] * ket[l][j]
            + bra[j][k] * ket[l][i]
            + bra[i][l] * ket[k][j]
            + bra[j][l] * ket[k][i]
            + bra[k][i] * ket[j][l]
            + bra[l][i] * ket[j][k]
            + bra[k][j] * ket[i][l]
            + bra[l][j] * ket[i][k]
        )
        ordinary_k_coefficient = -sum(
            bra[p][r] * ket[q][s] for p, q, r, s in permutations
        )

        self.assertEqual(mixed_quartet_coefficient, mixed_energy_derivative)
        self.assertNotEqual(mixed_quartet_coefficient, ordinary_k_coefficient)

    def test_umrsf_flag_is_scoped_to_umrsf_entry_point(self):
        source = compact(ENERGY.read_text())
        self.assertIn("subroutinetdhf_mrsf_energy_c", source)
        self.assertIn("inf%tddft%umrsf=.false.", source)
        self.assertIn("logical::previous_umrsf", source)
        self.assertIn("previous_umrsf=inf%tddft%umrsf", source)
        self.assertIn("inf%tddft%umrsf=previous_umrsf", source)

    def test_umrsf_jacobi_rotation_intent_and_diagonal_are_consistent(self):
        lib = LIB.read_text().lower()
        energy = compact(ENERGY.read_text())
        self.assertRegex(
            lib,
            r"real\(kind=dp\),\s*intent\(inout\),\s*dimension\(:,:\)\s*::\s*mo_a,\s*mo_b",
        )
        self.assertIn("mo_energy_work_a", energy)
        self.assertIn("mo_energy_work_a(i)=fa(i,i)", energy)
        self.assertIn("mo_energy_work_b(i)=fb(i,i)", energy)
        self.assertIn("callmrinivec(infos,mo_energy_work_a,mo_energy_work_b", energy)

    def test_umrsf_is_registered_as_a_tdhf_type(self):
        oqpdata = compact(OQPDATA.read_text())
        single = SINGLE_POINT.read_text().lower()

        self.assertIn("'umrsf'", oqpdata)
        self.assertIn("self._data.tddft.umrsf=td_type=='umrsf'", oqpdata)
        self.assertIn("'umrsf': oqp.tdhf_umrsf_z_vector", single)
        self.assertIn("'umrsf': oqp.tdhf_umrsf_gradient", single)

    def test_umrsf_energy_runtype_is_not_blocked(self):
        checker = _load_input_checker()
        report = checker.check_input_values(
            _umrsf_config("energy"), raise_error=False, emit=False
        )
        self.assertEqual(
            _umrsf_guard_errors(report),
            [],
            "UMRSF energy must not be rejected by the runtype guard:\n"
            + report.to_text(),
        )

    def test_umrsf_non_energy_gradient_runtypes_are_blocked_at_the_single_choke_point(self):
        checker = _load_input_checker()
        for runtype in UMRSF_BLOCKED_RUNTYPES:
            with self.subTest(runtype=runtype):
                report = checker.check_input_values(
                    _umrsf_config(runtype), raise_error=False, emit=False
                )
                guard_errors = _umrsf_guard_errors(report)
                self.assertEqual(
                    len(guard_errors),
                    1,
                    f"runtype={runtype} should raise exactly one UMRSF guard "
                    f"error, got {len(guard_errors)}:\n" + report.to_text(),
                )
                self.assertIn(runtype, guard_errors[0].value)

    def test_umrsf_grad_runtype_is_allowed(self):
        checker = _load_input_checker()
        report = checker.check_input_values(
            _umrsf_config("grad"), raise_error=False, emit=False
        )
        self.assertEqual(
            _umrsf_guard_errors(report),
            [],
            "UMRSF grad should not be rejected by the runtype guard:\n"
            + report.to_text(),
        )

    def test_umrsf_gradient_accepts_validated_global_hybrid_and_pure_gga(self):
        checker = _load_input_checker()
        for functional in ("bhhlyp", "blyp"):
            with self.subTest(functional=functional):
                config = _umrsf_config("grad")
                config["input"]["functional"] = functional
                report = checker.check_input_values(
                    config, raise_error=False, emit=False
                )
                self.assertFalse(
                    any("UMRSF analytic gradients" in error.message
                        for error in report.errors),
                    report.to_text(),
                )

    def test_umrsf_gradient_rejects_all_resolved_unsupported_functionals(self):
        checker = _load_input_checker()
        features = _libxc_aliases_by_feature()
        self.assertTrue({"camb3lyp", "lb07"} <= features["cam"])
        self.assertTrue({"tpss", "scan", "m062x"} <= features["meta"])
        self.assertTrue({"b2-plyp", "b2plyp"} <= features["dh"])
        self.assertIn("stg1x", features["spc"])

        gradient_runtypes = ("grad", "optimize", "meci", "mecp", "tci")
        for feature, aliases in features.items():
            for functional in aliases:
                for runtype in gradient_runtypes:
                    with self.subTest(
                            feature=feature, functional=functional,
                            runtype=runtype):
                        config = _umrsf_config(runtype)
                        config["input"]["functional"] = functional
                        report = checker.check_input_values(
                            config, raise_error=False, emit=False
                        )
                        self.assertTrue(
                            any("UMRSF analytic gradients" in error.message
                                for error in report.errors),
                            report.to_text(),
                        )

    def test_umrsf_gradient_rejects_explicit_cam_spc_and_ignored_solver(self):
        checker = _load_input_checker()
        gradient_runtypes = ("grad", "optimize", "meci", "mecp", "tci")
        for runtype in gradient_runtypes:
            config = _umrsf_config(runtype)
            config["dftgrid"] = {"cam_flag": True}
            report = checker.check_input_values(
                config, raise_error=False, emit=False
            )
            self.assertTrue(
                any("range-separated CAM/LRC" in error.message
                    for error in report.errors),
                report.to_text(),
            )

            for key in ("spc_coco", "spc_ovov", "spc_coov"):
                with self.subTest(runtype=runtype, key=key):
                    config = _umrsf_config(runtype)
                    config["tdhf"][key] = 0.35
                    report = checker.check_input_values(
                        config, raise_error=False, emit=False
                    )
                    self.assertTrue(
                        any("default spin-pair-coupling scales" in error.message
                            for error in report.errors),
                        report.to_text(),
                    )

            config = _umrsf_config(runtype)
            config["tdhf"]["z_solver"] = 3
            report = checker.check_input_values(
                config, raise_error=False, emit=False
            )
            self.assertTrue(
                any(error.path == "tdhf.z_solver" for error in report.errors),
                report.to_text(),
            )

    def test_umrsf_energy_still_allows_broader_parameterizations(self):
        checker = _load_input_checker()
        config = _umrsf_config("energy")
        config["input"]["functional"] = "cam-b3lyp"
        config["dftgrid"] = {"cam_flag": True}
        config["tdhf"].update({
            "spc_coco": 0.35,
            "spc_ovov": 0.35,
            "spc_coov": 0.35,
            "z_solver": 3,
        })
        report = checker.check_input_values(config, raise_error=False, emit=False)
        self.assertFalse(
            any("UMRSF analytic gradients" in error.message
                or error.path == "tdhf.z_solver"
                for error in report.errors),
            report.to_text(),
        )

    def test_umrsf_finite_difference_fixtures_are_root_one_only(self):
        for path, functional in (
            (UMRSF_BHHLYP_FD, "bhhlyp"),
            (UMRSF_BLYP_FD, "blyp"),
        ):
            with self.subTest(path=path.name):
                self.assertTrue(path.is_file())
                deck = configparser.ConfigParser(interpolation=None)
                deck.read(path)
                self.assertEqual(deck["input"]["runtype"], "grad")
                self.assertEqual(deck["input"]["method"], "tdhf")
                self.assertEqual(deck["input"]["functional"], functional)
                self.assertEqual(deck["input"]["basis"], "6-31g*")
                self.assertEqual(deck["input"]["omp_threads"], "1")
                self.assertEqual(deck["scf"]["type"], "uhf")
                self.assertEqual(deck["scf"]["multiplicity"], "3")
                self.assertEqual(deck["scf"]["converger_type"], "diis")
                self.assertEqual(deck["scf"]["conv"], "1e-10")
                self.assertEqual(deck["dftgrid"]["rad_npts"], "96")
                self.assertEqual(deck["dftgrid"]["ang_npts"], "302")
                self.assertEqual(deck["dftgrid"]["pruned"], "")
                self.assertEqual(deck["dftgrid"]["grid_ao_pruned"], "false")
                self.assertEqual(deck["tdhf"]["type"], "umrsf")
                self.assertEqual(deck["tdhf"]["nstate"], "3")
                self.assertEqual(deck["tdhf"]["target"], "1")
                self.assertEqual(deck["tdhf"]["multiplicity"], "1")
                self.assertEqual(deck["tdhf"]["maxit_zv"], "100")
                self.assertEqual(deck["tdhf"]["conv"], "1e-9")
                self.assertEqual(deck["tdhf"]["zvconv"], "1e-9")
                self.assertNotIn("z_solver", deck["tdhf"])
                self.assertEqual(deck["properties"]["grad"], "1")
                atoms = [
                    line.split()
                    for line in deck["input"]["system"].splitlines()
                    if line.strip()
                ]
                self.assertEqual([atom[0] for atom in atoms], ["H", "H", "C", "O"])
                self.assertTrue(
                    all(float(value) != 0.0 for atom in atoms for value in atom[1:])
                )

    def test_umrsf_zvector_entry_prepares_cached_response(self):
        zvec = compact(UMRSF_ZVEC.read_text())
        grad = compact(UMRSF_GRAD.read_text())
        tags = compact(TAGARRAY.read_text())

        self.assertNotIn("stub", zvec)
        self.assertIn("calltdhf_umrsf_build_response_gradient", zvec)
        self.assertIn("coupledalpha/betaresponse", zvec)
        self.assertIn("oqp_umrsf_response_gradient", zvec)
        self.assertIn("oqp_umrsf_response_gradient", grad)
        self.assertIn("alpha/betaz-vectorresponse", grad)
        self.assertIn("callhf_gradient(infos)", grad)
        self.assertIn("rtol=umrsf_z_requested_tolerance(infos)", grad)
        self.assertIn("mxit=min(ndofov,max(1,int(infos%control%maxit_zv)))", grad)
        self.assertIn("umrsfz-vectordidnotreachtherequestedrelativeresidual", grad)
        self.assertIn("callumrsf_zov_matvec(resov,rhsov,c_loc(ctx))", grad)
        self.assertIn("relres=errout/bnorm", grad)
        self.assertIn("remaining=max(0,mxit-pcg_iters)", grad)
        self.assertIn("z_full_rel=sqrt(z_rnorm2/z_bnorm2)", grad)
        self.assertIn("umrsfcoupledz-vectordidnotreachtherequestedrelativeresidual", grad)
        self.assertIn("if(infos%dft%cam_flag)then", grad)
        self.assertIn("if(infos%functional%needtau)then", grad)
        self.assertIn("if(infos%dft%dh_flag)then", grad)
        self.assertIn("abs(spc_coco-hfs)", grad)
        self.assertNotIn("callumrsf_grad_run_gates(inf,de2e_resp", grad)
        self.assertIn("oqp_umrsf_response_gradient", tags)

    def test_umrsf_energy_does_not_use_mrsf_transition_density_output_path(self):
        source = compact(ENERGY.read_text())
        self.assertIn("if(umrsf)then", source)
        self.assertIn("trden=0.0_dp", source)
        self.assertIn("else", source)
        self.assertIn("callget_mrsf_transition_density", source)

    def test_umrsf_energy_uses_spin_resolved_transition_dipoles(self):
        source = compact(ENERGY.read_text())
        lib = compact(LIB.read_text())
        self.assertIn("callget_umrsf_transition_dipole", source)
        self.assertIn("mo_a,mo_b,bvec_mo", source)
        self.assertIn("subroutineget_umrsf_transition_dipole", lib)
        # The substance of the routine is that the AO transition density is
        # built from BOTH orbital sets (C_alpha on the left, C_beta on the
        # right) rather than from mo_a alone. Assert that against the code --
        # not against the prose comment describing it, which would turn any
        # rewording of a comment into a CI failure.
        body = lib.split("subroutineget_umrsf_transition_dipole", 1)[1]
        body = body.split("endsubroutineget_umrsf_transition_dipole", 1)[0]
        self.assertIn("mo_a(:,nocb+1:)", body)
        self.assertIn("mo_b(:,nocb+1:)", body)
        self.assertIn("mo_b", body)

    def test_transition_dipoles_are_exposed_for_umrsf(self):
        """The UMRSF dipoles are real, so they must reach the tagarray.

        Only the state-interaction DENSITY is a UMRSF placeholder. Bundling the
        dipole export behind the same ``.not. umrsf`` guard hid the quantity
        from downstream analysis and from the regression references -- which is
        precisely why identically zero UMRSF dipoles went unnoticed.
        """
        source = compact(ENERGY.read_text())
        # trden stays MRSF-only ...
        density_guard = source.split("oqp_td_trans_density_mo", 1)[0]
        self.assertIn("if(.not.umrsf)then", density_guard)
        # ... but the dipole export must NOT sit inside a .not.umrsf block.
        between = source.split("oqp_td_trans_density_mo", 1)[1]
        between = between.split("oqp_td_trans_dipole", 1)[0]
        self.assertIn("endif", between,
                      "OQP_td_trans_dipole must be exported outside the "
                      ".not. umrsf guard that protects the trden placeholder")

    def test_transition_dipole_export_mirrors_reverse_pairs_for_umrsf(self):
        """Reverse transition pairs must not expose different magnitudes."""
        source = compact(ENERGY.read_text())
        export = source.split("oqp_td_trans_dipole", 1)[1]
        export = export.split("oqp_td_dip_ao", 1)[0]
        self.assertIn("dip_store(:,ist,jst)=dip_store(:,jst,ist)", export)
        self.assertNotIn("if(.not.umrsf)then", export)

    def test_transition_dipole_json_uses_fortran_order(self):
        """The tag is allocated as Fortran (3,nstates,nstates)."""
        source = compact(MOLECULE.read_text())
        block = source.split("oqp::td_trans_dipole", 1)[1]
        block = block.split("nmr_shielding", 1)[0]
        self.assertIn("ravel(order='c')", block)
        self.assertIn("reshape((3,nstates,nstates),order='f')", block)

    def test_transition_dipole_is_a_compared_regression_key(self):
        """Excitation energies alone did not catch all-zero UMRSF dipoles."""
        regression = _load_regression()

        entry = next((e for e in regression.REGISTRY
                      if e.key == "td_trans_dipole"), None)
        self.assertIsNotNone(
            entry, "td_trans_dipole must be in the regression registry so a "
                   "regression to zero transition dipoles fails the examples")
        self.assertTrue(entry.needs_excited)
        self.assertTrue(entry.phase_invariant)
        self.assertIn(entry, regression.keys_to_compare("energy", excited=True))

    def test_spin_pair_scaling_avoids_hfscale_division_by_zero(self):
        source = compact(ENERGY.read_text())
        self.assertIn("if(abs(infos%tddft%hfscale)>epsilon(1.0_dp))then", source)
        self.assertIn("spc_scale_coco", source)
        self.assertIn("spc_scale_ovov", source)
        self.assertIn("spc_scale_coov", source)


if __name__ == "__main__":
    unittest.main()
