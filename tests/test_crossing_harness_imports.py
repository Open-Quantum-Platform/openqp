"""Keep the test_crossing_objectives stub harness in step with its targets.

test_crossing_objectives.py imports libscipy and input_checker against stub
``oqp`` packages.  When a target gains a module-level ``oqp.*`` import the
harness does not provide, that whole file fails at collection (#534 added
``from oqp.utils.oqp_input import ...`` to input_checker).  This file reads
the harness and the targets without importing them, so it still runs then
and names the missing module.  Imports inside try blocks (which carry their
own fallback) and inside functions are not needed at load time.
"""
import ast
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HARNESS = ROOT / "tests" / "test_crossing_objectives.py"


def _literal(tree, name):
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == name
                for target in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError(f"{HARNESS.name} defines no {name}")


def _harness_modules(tree):
    provided = {name for name, _ in _literal(tree, "SOURCE_MODULES")}
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id in ("_stub", "_load") and node.args
                and isinstance(node.args[0], ast.Constant)
                and isinstance(node.args[0].value, str)):
            provided.add(node.args[0].value)
    return provided


def _is_source_module(name):
    """True when ``name`` is a module file or package under pyoqp/."""
    base = ROOT / "pyoqp" / Path(*name.split("."))
    return base.with_suffix(".py").is_file() or (base / "__init__.py").is_file()


def _module_level_oqp_imports(path):
    for node in ast.parse(path.read_text()).body:
        if isinstance(node, ast.ImportFrom) and node.module:
            names = [node.module]
            # ``from oqp.utils import name`` imports the submodule
            # oqp.utils.name when one exists; that needs providing too
            names += [f"{node.module}.{alias.name}" for alias in node.names
                      if _is_source_module(f"{node.module}.{alias.name}")]
        elif isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        else:
            continue
        for name in names:
            if name == "oqp" or name.startswith("oqp."):
                yield name


class CrossingHarnessImportTest(unittest.TestCase):
    def test_harness_provides_every_module_level_oqp_import(self):
        tree = ast.parse(HARNESS.read_text())
        provided = _harness_modules(tree)
        targets = _literal(tree, "LOADED_FILES")
        self.assertTrue(targets)
        for relative in targets:
            for name in _module_level_oqp_imports(ROOT / relative):
                with self.subTest(file=relative, module=name):
                    self.assertIn(
                        name, provided,
                        f"{relative} imports {name} at module level; add it "
                        f"to SOURCE_MODULES (pure Python) or stub it in "
                        f"{HARNESS.name}")

    def test_from_package_import_of_a_submodule_is_seen(self):
        import tempfile

        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder) / "target.py"
            target.write_text(
                "from oqp.utils import oqp_input, not_a_module\n"
                "import numpy\n"
                "try:\n    from oqp.utils import spin_states\n"
                "except ModuleNotFoundError:\n    pass\n")
            found = list(_module_level_oqp_imports(target))
        # the package, its real submodule, nothing from the try block
        self.assertEqual(found, ["oqp.utils", "oqp.utils.oqp_input"])


if __name__ == "__main__":
    unittest.main()
