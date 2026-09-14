"""Reject failed SCF before state overlaps, excited states, or nuclear forces."""
import ast
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock
import pytest

SOURCE = Path(__file__).resolve().parents[1] / 'pyoqp/oqp/library/namd.py'


def driver_class(namespace):
    # Execute the actual NAMD methods without importing a native runtime.
    tree = ast.parse(SOURCE.read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'NAMD')
    names = {'_electronic', '_require_converged_reference', '_active_gradient'}
    cls.body = [n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name in names]
    exec(compile(ast.Module(body=[cls], type_ignores=[]), str(SOURCE), 'exec'), namespace)
    return namespace['NAMD']


def make_driver(converged):
    sp = Mock()
    namespace = {'SinglePoint': Mock(return_value=sp), 'Gradient': Mock(),
                 'BasisOverlap': Mock(), 'LastStep': Mock()}
    cls = driver_class(namespace)
    d = cls()
    d.mol = SimpleNamespace(mol_energy=SimpleNamespace(SCF_converged=converged))
    d.mo_reuse = False
    d.ref_follow = 'off'
    d.scf_fail = 'escalate'
    return d, namespace, sp


def test_failed_reference_return_cannot_reach_excitation_or_overlap():
    d, ns, sp = make_driver(False)
    # Even an erroneous backend that returns an energy must not advance NAMD.
    sp.reference.return_value = -414.0
    with pytest.raises(RuntimeError, match='SCF did not converge'):
        d._electronic(with_overlap=True)
    ns['BasisOverlap'].assert_not_called()
    sp.excitation.assert_not_called()
    ns['LastStep'].assert_not_called()


def test_failed_reference_cannot_reach_gradient():
    d, ns, _ = make_driver(False)
    with pytest.raises(RuntimeError, match='no force'):
        d._active_gradient()
    ns['Gradient'].assert_not_called()


def test_successful_reference_still_reaches_excitation():
    d, ns, sp = make_driver(True)
    sp.reference.return_value = -414.0
    d._electronic(with_overlap=False)
    sp.excitation.assert_called_once_with(-414.0)
    ns['LastStep'].assert_called_once()
