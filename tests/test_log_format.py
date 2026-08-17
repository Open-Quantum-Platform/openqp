"""Source-only checks for the common OpenQP text-log grammar."""

import ast
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _load_log_format():
    path = ROOT / "pyoqp" / "oqp" / "utils" / "log_format.py"
    spec = importlib.util.spec_from_file_location("oqp_log_format_under_test", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


LOG = _load_log_format()


def _literal_set_assignment(path, name):
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == name
               for target in node.targets):
            return set(ast.literal_eval(node.value))
    raise AssertionError(f"assignment {name} not found in {path}")


def _runner_runtypes():
    path = ROOT / "pyoqp" / "oqp" / "pyoqp.py"
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign) or not isinstance(node.value, ast.Dict):
            continue
        if not any(isinstance(target, ast.Attribute) and target.attr == "run_func"
                   for target in node.targets):
            continue
        return {
            key.value for key in node.value.keys
            if isinstance(key, ast.Constant) and isinstance(key.value, str)
        }
    raise AssertionError("Runner.run_func mapping not found")


def _literal_dump_log_sections():
    sections = set()
    for path in (ROOT / "pyoqp" / "oqp").rglob("*.py"):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if not isinstance(node.func, ast.Name) or node.func.id != "dump_log":
                continue
            for keyword in node.keywords:
                if (keyword.arg == "section" and
                        isinstance(keyword.value, ast.Constant) and
                        isinstance(keyword.value.value, str)):
                    sections.add(keyword.value.value)
    return sections


def _implicit_progress_dump_log_paths():
    paths = set()
    for path in (ROOT / "pyoqp" / "oqp").rglob("*.py"):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if not isinstance(node.func, ast.Name) or node.func.id != "dump_log":
                continue
            section = next(
                (keyword.value for keyword in node.keywords
                 if keyword.arg == "section"),
                None,
            )
            if section is None or (
                isinstance(section, ast.Constant) and section.value in (None, "")
            ):
                paths.add(path.relative_to(ROOT / "pyoqp" / "oqp").as_posix())
    return paths


def test_section_heading_snapshot():
    assert LOG.format_log_section(
        "PyOQP: Calculation request", LOG.INPUT_REFERENCE
    ) == (
        "\n"
        "   ========================================================================\n"
        "   PyOQP LOG | INPUT AND REFERENCE\n"
        "   PyOQP: Calculation request\n"
        "   ========================================================================\n"
    )


def test_common_field_and_energy_format_snapshots():
    assert LOG.format_log_fields((
        ("Method", "MRSF-TDDFT"),
        ("Reference", "triplet ROHF (internal working reference)"),
    )) == (
        "   PyOQP Method:                      MRSF-TDDFT\n"
        "   PyOQP Reference:                   triplet ROHF (internal working reference)"
    )
    assert LOG.format_energy(-75.1234567890123) == "-75.1234567890  "
    assert LOG.format_unit("Gradient", "Hartree/Bohr") == (
        "   PyOQP Gradient unit:               Hartree/Bohr"
    )


def test_every_supported_and_dispatched_runtype_has_an_ordered_log_profile():
    supported = _literal_set_assignment(
        ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py",
        "SUPPORTED_RUNTYPES",
    )
    routes = supported | _runner_runtypes()
    missing = routes - set(LOG.RUNTYPE_LOG_CATEGORIES)
    assert not missing, f"run types without a log profile: {sorted(missing)}"
    for runtype in routes:
        categories = LOG.RUNTYPE_LOG_CATEGORIES[runtype]
        indices = [LOG.LOG_SECTION_ORDER.index(category) for category in categories]
        assert indices == sorted(indices), runtype
        assert categories[0] == LOG.INPUT_REFERENCE
        assert categories[-1] == LOG.TERMINATION
        assert LOG.PROGRESS in categories


def test_every_literal_dump_log_section_is_classified():
    sections = _literal_dump_log_sections() - {""}
    unclassified = {
        section for section in sections
        if LOG.section_category(section) == LOG.PROGRESS
    }
    assert not unclassified, sorted(unclassified)


def test_implicit_progress_sections_are_declared_across_calculation_flows():
    assert LOG.PROGRESS in LOG.LOG_SECTION_ORDER
    assert LOG.section_category(None) == LOG.PROGRESS
    assert LOG.section_category("") == LOG.PROGRESS

    paths = _implicit_progress_dump_log_paths()
    expected = {
        "library/single_point.py",  # gradients and derivative properties
        "library/libscipy.py",      # optimizations and reaction paths
        "library/liboqp.py",        # native optimizations, NEB, IRC, and MEP
        "library/namd.py",          # molecular dynamics
    }
    assert expected <= paths, sorted(expected - paths)


def test_stage_specific_sections_use_forward_log_categories():
    assert LOG.section_category("correlation") == LOG.CONVERGENCE
    assert LOG.section_category("QM/MM") == LOG.CONVERGENCE
    assert LOG.section_category("text") == LOG.CONVERGENCE
    assert LOG.section_category("dftb_state_summary") == LOG.ENERGIES
    assert LOG.section_category("dftd") == LOG.INPUT_REFERENCE
    assert LOG.section_category("fci") == LOG.ENERGIES


def test_native_banner_is_appended_after_log_initialization():
    source = (ROOT / "pyoqp" / "oqp" / "pyoqp.py").read_text()
    start = source.index("dump_log(self.mol, title='', section='start'")
    banner = source.index("oqp.oqp_banner(self.mol)", start)
    calculation = source.index("section='calculation'", banner)
    assert start < banner < calculation


def test_mpi_banner_collective_precedes_single_root_write():
    source = (ROOT / "source" / "modules" / "oqp_banner.F90").read_text()
    collective = source.index("call pe%get_hostnames(hostnames)")
    root_guard = source.index("if (pe%rank == 0) then", collective)
    opened = source.index("open (newunit=iw", root_guard)
    banner = source.index("OpenQP: Open Quantum Platform", opened)
    closed = source.index("close (iw)", banner)
    root_end = source.index("endif", closed)

    assert collective < root_guard < opened < banner < closed < root_end
    assert source.count("open (newunit=iw") == 1
    assert source.count("OpenQP: Open Quantum Platform") == 1
    assert source.count("close (iw)") == 1
    root_block = source[root_guard:root_end]
    assert root_block.count("write(iw") == source.count("write(iw")


def test_compatibility_document_names_stable_markers_and_units():
    text = (ROOT / "docs" / "logging.md").read_text()
    for category in LOG.LOG_SECTION_ORDER:
        assert f"`{category}`" in text
    for marker in (
        "PyOQP started at",
        "PyOQP method:",
        "PyOQP electronic energies",
        "PyOQP electronic gradients",
        "PyOQP terminated at",
    ):
        assert marker in text
    for unit in ("Hartree", "Hartree/Bohr", "cm^-1", "km/mol"):
        assert unit in text
