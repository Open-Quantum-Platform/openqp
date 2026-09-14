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
    names = {'_electronic', '_require_converged_reference', '_active_gradient', '_reference_fresh_guess'}
    cls.body = [n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name in names]
    exec(compile(ast.Module(body=[cls], type_ignores=[]), str(SOURCE), 'exec'), namespace)
    return namespace['NAMD']


def make_driver(converged):
    sp = Mock()
    namespace = {'SinglePoint': Mock(return_value=sp), 'Gradient': Mock(),
                 'BasisOverlap': Mock(), 'LastStep': Mock(),
                 'SCFnotConverged': type('SCFnotConverged', (Exception,), {}),
                 'dump_log': Mock(), 'is_tb_method': lambda method: method in ('dftb', 'xtb')}
    cls = driver_class(namespace)
    d = cls()
    d.mol = SimpleNamespace(mol_energy=SimpleNamespace(SCF_converged=converged),
                            config={'input': {'method': 'tdhf'}})
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


def test_tight_binding_uses_its_own_scc_acceptance():
    d, _, _ = make_driver(False)
    d.mol.config['input']['method'] = 'dftb'
    d._require_converged_reference()


def recovery_driver():
    d, ns, sp = make_driver(True)
    d.mo_reuse = True
    d.scf_guess_retry = True
    d._scf_fallback_steps = 0
    d.prev_xyz, d.prev_data = [], {}
    d.mol.config.update(scf={'converger_type': 'soscf', 'escalation': 'soscf,trah'},
                        guess={'type': 'previous'})
    d.mol.data = SimpleNamespace(set_scf_converger_type=Mock())
    return d, ns, sp


def test_failed_scf_retries_fresh_guess_and_restores_settings():
    d, ns, sp = recovery_driver()
    guesses = []
    def reference():
        guesses.append(d.mol.config['guess']['type'])
        if len(guesses) == 1:
            raise ns['SCFnotConverged']('failed criterion')
        return [-414.0]
    sp.reference.side_effect = reference
    d._electronic(with_overlap=True)
    assert guesses == ['previous', 'huckel']
    assert d._scf_fallback_steps == 1
    assert d.mol.config['guess']['type'] == 'previous'
    assert d.mol.config['scf']['converger_type'] == 'soscf'
    sp.excitation.assert_called_once_with([-414.0])


def test_disabled_guess_retry_stops_after_failed_scf():
    d, ns, sp = recovery_driver()
    d.scf_guess_retry = False
    sp.reference.side_effect = ns['SCFnotConverged']('failed criterion')
    with pytest.raises(ns['SCFnotConverged']):
        d._electronic(with_overlap=True)
    assert sp.reference.call_count == 1
    sp.excitation.assert_not_called()


def test_unrelated_runtime_error_does_not_trigger_guess_retry():
    d, ns, sp = recovery_driver()
    sp.reference.side_effect = RuntimeError('invalid integral buffer')
    with pytest.raises(RuntimeError, match='invalid integral'):
        d._electronic(with_overlap=True)
    assert sp.reference.call_count == 1
    sp.excitation.assert_not_called()


def test_failed_fresh_guess_cannot_advance_and_restores_settings():
    d, ns, sp = recovery_driver()
    sp.reference.side_effect = ns['SCFnotConverged']('failed criterion')
    with pytest.raises(ns['SCFnotConverged']):
        d._electronic(with_overlap=True)
    assert sp.reference.call_count == 2
    assert d.mol.config['guess']['type'] == 'previous'
    assert d.mol.config['scf']['converger_type'] == 'soscf'
    sp.excitation.assert_not_called()
    ns['BasisOverlap'].assert_not_called()


def test_false_success_from_fresh_guess_cannot_advance():
    d, ns, sp = recovery_driver()
    d.mol.mol_energy.SCF_converged = False
    sp.reference.side_effect = [ns['SCFnotConverged']('failed criterion'), -414.0]
    with pytest.raises(RuntimeError, match='SCF did not converge'):
        d._electronic(with_overlap=True)
    sp.excitation.assert_not_called()
    ns['BasisOverlap'].assert_not_called()
