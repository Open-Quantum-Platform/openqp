"""Every bind(C) procedure that PyOQP can call must be exported by liboqp.

GCC 16.2.0 gives a bind(C) procedure in a module with a PRIVATE default
hidden visibility unless the procedure is listed PUBLIC (GCC PR126872).  The
library still links, but cffi then cannot find the entry point, so whole
methods disappear at import time.  The first test checks the source and
therefore also fails on compilers without the regression (the CI compilers);
the second checks the library that is actually installed.
"""

import ctypes
import importlib.util
import re
import unittest
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1] / "source"

MODULE = re.compile(r"^\s*module\s+(?!procedure\b)(\w+)\s*$", re.I)
END_MODULE = re.compile(r"^\s*end\s*module\b", re.I)
INTERFACE = re.compile(r"^\s*(abstract\s+)?interface\b", re.I)
END_INTERFACE = re.compile(r"^\s*end\s*interface\b", re.I)
DEFAULT_PRIVATE = re.compile(r"^\s*private\s*$", re.I)
ACCESS = re.compile(r"^\s*(public|private)\s*(?:::)?\s*(\w+(?:\s*,\s*\w+)*)\s*$", re.I)
PROCEDURE = re.compile(r"\b(?:subroutine|function)\s+(\w+)\s*\(", re.I)
BIND_C = re.compile(r"bind\s*\(\s*c\s*,\s*name\s*=\s*[\"'](\w+)[\"']", re.I)
OPENER = re.compile(r"^\s*(?:end\s+)", re.I)


def statements(text):
    """Free-form statements with continuations joined and comments removed."""
    buf = ""
    for raw in text.splitlines():
        line = raw.split("!")[0].strip()
        if line.startswith("&"):
            line = line[1:]
        if line.endswith("&"):
            buf += line[:-1] + " "
            continue
        yield buf + line
        buf = ""


def bind_c_definitions():
    """(file, module, procedure, label, is_public) for bind(C) procedures defined in modules."""
    found = []
    for path in sorted(SOURCE.rglob("*")):
        if path.suffix.lower() != ".f90" or not path.is_file():
            continue
        module = None
        for stmt in statements(path.read_text(errors="replace")):
            m = MODULE.match(stmt)
            if m and module is None:
                module, private_default, public, private, procs, depth = (
                    m.group(1), False, set(), set(), [], 0)
                continue
            if module is None:
                continue
            if END_MODULE.match(stmt):
                for name, label in procs:
                    key = name.lower()
                    is_public = key in public or (not private_default and key not in private)
                    found.append((path.relative_to(SOURCE), module, name, label, is_public))
                module = None
                continue
            if INTERFACE.match(stmt):
                depth += 1
                continue
            if END_INTERFACE.match(stmt):
                depth = max(0, depth - 1)
                continue
            if DEFAULT_PRIVATE.match(stmt):
                private_default = True
                continue
            a = ACCESS.match(stmt)
            if a:
                names = {n.strip().lower() for n in a.group(2).split(",")}
                (public if a.group(1).lower() == "public" else private).update(names)
                continue
            if depth or OPENER.match(stmt):
                continue
            p, b = PROCEDURE.search(stmt), BIND_C.search(stmt)
            if p and b:
                procs.append((p.group(1), b.group(1)))
    return found


class BindCExportTest(unittest.TestCase):

    def test_source_scan_finds_entry_points(self):
        labels = {d[3] for d in bind_c_definitions()}
        for expected in ("oqp_clean", "hf_energy", "mp2_energy", "soc_mrsf"):
            self.assertIn(expected, labels)

    def test_bind_c_procedures_are_public(self):
        hidden = [f"{f}: {m}::{n} (bind name {lbl})"
                  for f, m, n, lbl, is_public in bind_c_definitions() if not is_public]
        self.assertEqual(hidden, [], "bind(C) procedures must be PUBLIC, or GCC 16.2 "
                         "(PR126872) drops them from liboqp:\n" + "\n".join(hidden))

    @unittest.skipIf(importlib.util.find_spec("oqp") is None, "oqp is not installed")
    def test_installed_library_exports_entry_points(self):
        import oqp
        from oqp import runtime
        root, suffix = runtime.resolve_oqp_root()
        lib = ctypes.CDLL(str(runtime.library_path(root, suffix)))
        missing = {lbl for _, _, _, lbl, _ in bind_c_definitions() if not hasattr(lib, lbl)}
        # ENABLE_DFTD4=OFF leaves dftd4_interface.F90 out of liboqp on purpose
        # (source/CMakeLists.txt), and PyOQP tolerates exactly these entry
        # points being absent.  Only a backend that is missing as a whole is
        # excused; losing part of it is still a visibility defect.
        optional = set(oqp._OPTIONAL_ENTRY_POINTS)
        if optional and optional <= missing:
            missing -= optional
        missing = sorted(missing)
        self.assertEqual(missing, [], "liboqp does not export these bind(C) entry points: "
                         + ", ".join(missing))


if __name__ == "__main__":
    unittest.main()
