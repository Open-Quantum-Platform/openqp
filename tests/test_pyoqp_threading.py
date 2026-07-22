"""Pre-import thread-policy tests for native DFTB inputs."""

import ast
import os
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "pyoqp" / "oqp" / "pyoqp.py"


def _load_preimport_helpers():
    """Load only the dependency-free functions that execute before ``import oqp``."""
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"), filename=str(SOURCE))
    names = {
        "_strip_preimport_comments",
        "_apply_omp_threads_from_input",
        "_set_threading_defaults",
    }
    nodes = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    ]
    namespace = {
        "__name__": "_pyoqp_preimport_threading_test",
        "os": os,
        "sys": sys,
    }
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), "exec"), namespace)
    return namespace


def test_native_dftb_input_seeds_mkl_from_explicit_omp_request(monkeypatch):
    helpers = _load_preimport_helpers()
    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)
    monkeypatch.delenv("MKL_NUM_THREADS", raising=False)
    with tempfile.NamedTemporaryFile("w", suffix=".oqp", delete=False) as handle:
        handle.write("mrsf-tddftb(nstate=3)\nomp_threads=8\n")
        path = handle.name
    try:
        helpers["_apply_omp_threads_from_input"](["openqp", path])
    finally:
        os.unlink(path)
    assert os.environ["OMP_NUM_THREADS"] == "8"
    assert os.environ["MKL_NUM_THREADS"] == "8"


def test_non_dftb_input_and_explicit_mkl_policy_are_preserved(monkeypatch):
    helpers = _load_preimport_helpers()
    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)
    monkeypatch.setenv("MKL_NUM_THREADS", "3")
    with tempfile.NamedTemporaryFile("w", suffix=".oqp", delete=False) as handle:
        handle.write("dft/pbe0/def2-svp\nomp_threads=8\n")
        path = handle.name
    try:
        helpers["_apply_omp_threads_from_input"](["openqp", path])
    finally:
        os.unlink(path)
    assert os.environ["OMP_NUM_THREADS"] == "8"
    assert os.environ["MKL_NUM_THREADS"] == "3"
