import configparser
import importlib.util
import re
import sys
import tempfile
import types
import unittest
from functools import lru_cache
from pathlib import Path
from unittest import mock


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
GRAD_VALIDATOR = ROOT / "tools" / "validate_gradients.py"
UMRSF_BHHLYP_FD = ROOT / "tools" / "validation_inputs" / "H2CO_BHHLYP_UMRSF_GRADIENT.inp"
UMRSF_BLYP_FD = ROOT / "tools" / "validation_inputs" / "H2CO_BLYP_UMRSF_GRADIENT.inp"

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


@lru_cache(maxsize=1)
def _load_gradient_validator():
    oqp_module = types.ModuleType("oqp")
    oqp_module.__path__ = []
    oqp_module.oqp_banner = lambda mol: None
    pyoqp_module = types.ModuleType("oqp.pyoqp")
    pyoqp_module.Runner = object
    library_module = types.ModuleType("oqp.library")
    library_module.__path__ = []
    single_point_module = types.ModuleType("oqp.library.single_point")
    single_point_module.SinglePoint = object
    stubs = {
        "oqp": oqp_module,
        "oqp.pyoqp": pyoqp_module,
        "oqp.library": library_module,
        "oqp.library.single_point": single_point_module,
    }
    spec = importlib.util.spec_from_file_location(
        "umrsf_gradient_validator_under_test", GRAD_VALIDATOR
    )
    module = importlib.util.module_from_spec(spec)
    with mock.patch.dict(sys.modules, stubs):
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
            r"this%bden\(ch,[^)]*\)\*this%dden\(ch,[^)]*\)"
        )
        expected_mixed_terms = (
            "this%bden(ch,i1,k1)*this%dden(ch,l1,j1)",
            "this%bden(ch,j1,k1)*this%dden(ch,l1,i1)",
            "this%bden(ch,i1,l1)*this%dden(ch,k1,j1)",
            "this%bden(ch,j1,l1)*this%dden(ch,k1,i1)",
            "this%bden(ch,k1,i1)*this%dden(ch,j1,l1)",
            "this%bden(ch,l1,i1)*this%dden(ch,j1,k1)",
            "this%bden(ch,k1,j1)*this%dden(ch,i1,l1)",
            "this%bden(ch,l1,j1)*this%dden(ch,i1,k1)",
        )
        expected_ordinary_terms = (
            "this%bden(ch,i1,k1)*this%dden(ch,j1,l1)",
            "this%bden(ch,j1,k1)*this%dden(ch,i1,l1)",
            "this%bden(ch,i1,l1)*this%dden(ch,j1,k1)",
            "this%bden(ch,j1,l1)*this%dden(ch,i1,k1)",
            "this%bden(ch,k1,i1)*this%dden(ch,l1,j1)",
            "this%bden(ch,l1,i1)*this%dden(ch,k1,j1)",
            "this%bden(ch,k1,j1)*this%dden(ch,l1,i1)",
            "this%bden(ch,l1,j1)*this%dden(ch,k1,i1)",
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

    def test_umrsf_validator_overrides_only_exact_input_keys(self):
        validator = _load_gradient_validator()
        source_text = """[input]
runtype=energy
[scf]
conv=1e-6
converger_type=trah
conv_extra=keep-scf
[tdhf]
conv=2e-6
zvconv=3e-6
nstate=9
nstate_buffer=keep-td
[properties]
grad=2
gradient_mode=keep-properties
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "source.inp"
            output = Path(tmpdir) / "output.inp"
            source.write_text(source_text)
            validator._write_input(
                source, output, runtype="grad", grad_list=[1], nstate=3
            )
            parsed = configparser.ConfigParser(interpolation=None)
            parsed.read(output)

        self.assertEqual(parsed["input"]["runtype"], "grad")
        self.assertEqual(parsed["scf"]["conv"], validator.SCF_CONV)
        self.assertEqual(parsed["scf"]["converger_type"], "trah")
        self.assertEqual(parsed["scf"]["conv_extra"], "keep-scf")
        self.assertEqual(parsed["tdhf"]["conv"], validator.TD_CONV)
        self.assertEqual(parsed["tdhf"]["zvconv"], validator.ZV_CONV)
        self.assertEqual(parsed["tdhf"]["nstate"], "3")
        self.assertEqual(parsed["tdhf"]["nstate_buffer"], "keep-td")
        self.assertEqual(parsed["properties"]["grad"], "1")
        self.assertEqual(
            parsed["properties"]["gradient_mode"], "keep-properties"
        )

    def test_umrsf_validator_behavior_forces_separated_root_one(self):
        validator = _load_gradient_validator()
        self.assertEqual(validator.METHODS["umrsf-bhhlyp"][2], "td-root1")
        self.assertEqual(validator.METHODS["umrsf-blyp"][2], "td-root1")
        captured = {}

        class FakeMol:
            data = {"natom": 1}
            grads = validator.np.arange(12.0)

        class FakeRunner:
            mol = FakeMol()

            @staticmethod
            def run():
                return None

        def fake_write(src, dst, runtype, grad_list, nstate=None):
            captured.update(
                runtype=runtype, grad_list=list(grad_list), nstate=nstate
            )
            return dst

        with tempfile.TemporaryDirectory() as tmpdir, \
                mock.patch.object(validator, "_write_input", side_effect=fake_write), \
                mock.patch.object(validator, "_runner", return_value=FakeRunner()):
            gradients, natom, states = validator.analytical_gradients(
                "umrsf-bhhlyp", tmpdir, nstate=99, compute_nstate=3
            )

        self.assertEqual(captured, {
            "runtype": "grad", "grad_list": [1], "nstate": 3
        })
        self.assertEqual(natom, 1)
        self.assertEqual(states, [1])
        self.assertEqual(list(gradients), [1])
        self.assertTrue(
            validator.np.array_equal(
                gradients[1], validator.np.array([[3.0, 4.0, 5.0]])
            )
        )

        calls = {}

        def fake_analytical(key, workdir, nstate, compute_nstate):
            calls["analytical"] = (nstate, compute_nstate)
            return {1: validator.np.zeros((1, 3))}, 1, [1]

        def fake_numerical(key, workdir, states, dx, nstate, compute_nstate):
            calls["numerical"] = (list(states), nstate, compute_nstate)
            return {1: validator.np.zeros((1, 3))}

        with tempfile.TemporaryDirectory() as tmpdir, \
                mock.patch.object(
                    validator, "analytical_gradients", side_effect=fake_analytical
                ), mock.patch.object(
                    validator, "numerical_gradients", side_effect=fake_numerical
                ):
            result = validator.validate(
                "umrsf-bhhlyp", dx=1.0e-3, nstate=99, workdir=tmpdir
            )
        self.assertNotIn("error", result)
        self.assertEqual(calls["analytical"], (99, 3))
        self.assertEqual(calls["numerical"], ([1], 99, 3))

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

    def test_spin_pair_scaling_avoids_hfscale_division_by_zero(self):
        source = compact(ENERGY.read_text())
        self.assertIn("if(abs(infos%tddft%hfscale)>epsilon(1.0_dp))then", source)
        self.assertIn("spc_scale_coco", source)
        self.assertIn("spc_scale_ovov", source)
        self.assertIn("spc_scale_coov", source)


if __name__ == "__main__":
    unittest.main()
