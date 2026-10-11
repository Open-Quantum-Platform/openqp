"""Reject missing, skipped or failing named NEVPT2 CI evidence."""
import importlib.util
import copy
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest


@pytest.fixture
def driver():
    path = Path(__file__).parent / 'native' / 'run_nevpt2_qualification.py'
    spec = importlib.util.spec_from_file_location('nevpt2_ci_receipts', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _report(tmp_path, names, outcome=None):
    suite = ET.Element('testsuite')
    for name in names:
        case = ET.SubElement(suite, 'testcase', classname='tests.test_nevpt2_subspaces', name=name)
        if outcome is not None:
            ET.SubElement(case, outcome)
    path = tmp_path / 'junit.xml'
    ET.ElementTree(suite).write(path)
    return path


def test_complete_named_passes_are_accepted(driver, tmp_path):
    # Parametrized node names must still identify the committed test function.
    # This independent manifest prevents a shortened driver list passing its own oracle.
    names = [
        'test_final_contraction_matches_numpy[Sr]', 'test_final_contraction_matches_numpy[Si]',
        'test_final_contraction_matches_numpy[Sijrs]', 'test_final_contraction_matches_numpy[Sijr]',
        'test_final_contraction_matches_numpy[Srsi]', 'test_final_contraction_matches_numpy[Srs]',
        'test_final_contraction_matches_numpy[Sij]', 'test_final_contraction_matches_numpy[Sir]',
        'test_retained_subspace_outputs_survive_later_calls[Sr]',
        'test_retained_subspace_outputs_survive_later_calls[Si]',
        'test_retained_subspace_outputs_survive_later_calls[Sijrs]',
        'test_retained_subspace_outputs_survive_later_calls[Sijr]',
        'test_retained_subspace_outputs_survive_later_calls[Srsi]',
        'test_retained_subspace_outputs_survive_later_calls[Srs]',
        'test_retained_subspace_outputs_survive_later_calls[Sij]',
        'test_retained_subspace_outputs_survive_later_calls[Sir]',
        'test_fused_rdms_and_f3_match_numpy[cas2]', 'test_fused_rdms_and_f3_match_numpy[cas3]',
        'test_fused_rdms_and_f3_match_numpy[cas4]', 'test_retained_fused_outputs_survive_later_calls',
        'test_actual_native_allocation_failures_and_recovery',
        'test_malloc_fault_shim_survives_late_tls_modules',
        'test_whole_energy_matches_backend_unavailable[core1-cas2-virt2]',
        'test_whole_energy_matches_backend_unavailable[core0-cas3-virt2]',
    ]
    assert len(names) == 24
    assert set(driver.TESTS) == set(names)
    assert len(driver.require_named_passes(_report(tmp_path, names))) == 24


@pytest.mark.parametrize('outcome', ['skipped', 'failure', 'error'])
def test_success_summary_cannot_hide_nonpassing_node(driver, tmp_path, outcome):
    with pytest.raises(RuntimeError, match='did not pass'):
        driver.require_named_passes(_report(tmp_path, driver.TESTS, outcome))


def test_success_summary_cannot_hide_missing_native_node(driver, tmp_path):
    with pytest.raises(RuntimeError, match='Missing or unexpected named'):
        driver.require_named_passes(_report(tmp_path, driver.TESTS[:-1]))


def test_missing_parameter_case_is_not_complete_execution(driver, tmp_path):
    names = list(driver.TESTS)
    with pytest.raises(RuntimeError, match='Missing or unexpected named'):
        driver.require_named_passes(_report(tmp_path, names[1:]))


def test_empty_report_is_not_native_execution(driver, tmp_path):
    with pytest.raises(RuntimeError, match='No named'):
        driver.require_named_passes(_report(tmp_path, []))



def test_duplicate_case_cannot_replace_another_parameter(driver, tmp_path):
    names = list(driver.TESTS)
    names[0] = names[1]
    with pytest.raises(RuntimeError, match='Duplicate named'):
        driver.require_named_passes(_report(tmp_path, names))



def test_unexpected_parameter_cannot_replace_required_case(driver, tmp_path):
    names = list(driver.TESTS)
    names[0] = 'test_final_contraction_matches_numpy[not_a_subspace]'
    with pytest.raises(RuntimeError, match='Missing or unexpected named'):
        driver.require_named_passes(_report(tmp_path, names))


@pytest.mark.parametrize('key,value', [('errors', '1'), ('failures', '1'),
                                      ('skipped', '1'), ('tests', '23'), ('tests', 'invalid')])
def test_suite_summary_cannot_contradict_complete_nodes(driver, tmp_path, key, value):
    path = _report(tmp_path, driver.TESTS)
    tree = ET.parse(path)
    tree.getroot().set(key, value)
    tree.write(path)
    with pytest.raises(RuntimeError, match='suite'):
        driver.require_named_passes(path)


def test_suite_error_cannot_hide_outside_testcases(driver, tmp_path):
    path = _report(tmp_path, driver.TESTS)
    tree = ET.parse(path)
    ET.SubElement(tree.getroot(), 'error', message='collection failed')
    tree.write(path)
    with pytest.raises(RuntimeError, match='did not pass'):
        driver.require_named_passes(path)


def test_nested_pytest_summary_is_accepted(driver, tmp_path):
    path = _report(tmp_path, driver.TESTS)
    suite = ET.parse(path).getroot()
    suite.attrib.update(tests='24', errors='0', failures='0', skipped='0')
    root = ET.Element('testsuites')
    root.append(suite)
    ET.ElementTree(root).write(path)
    assert len(driver.require_named_passes(path)) == 24


def test_malformed_xml_cannot_qualify(driver, tmp_path):
    path = tmp_path / 'junit.xml'
    path.write_text('<testsuites>')
    with pytest.raises(RuntimeError, match='Malformed'):
        driver.require_named_passes(path)


@pytest.fixture
def runtime_identity():
    return {'library': '/installed/lib/liboqp.so', 'library_sha256': 'a' * 64,
            'package': '/installed/oqp/__init__.py',
            'nevpt2_module': '/installed/oqp/library/nevpt2_sc.py',
            'omp_requested': 4, 'fci_thread_argument': 4, 'omp_threads_from_env': True,
            'thread_environment': {'OMP_NUM_THREADS': '4', 'OQP_FCI_THREADS': '4',
                'OMP_DYNAMIC': 'FALSE', 'OMP_MAX_ACTIVE_LEVELS': '1',
                'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'}}


@pytest.mark.parametrize('filename', ['pytest-runtime.json', 'fault-runtime.json'])
def test_matching_actual_runtime_is_accepted(driver, tmp_path, runtime_identity, filename):
    path = tmp_path / filename
    path.write_text(json.dumps(runtime_identity))
    assert driver.require_runtime_receipt(path, runtime_identity) == runtime_identity


@pytest.mark.parametrize('key', ['library', 'library_sha256', 'package', 'nevpt2_module',
                                'omp_requested', 'fci_thread_argument',
                                'omp_threads_from_env', 'thread_environment'])
@pytest.mark.parametrize('mutation', ['absent', 'different', 'wrong_type'])
def test_runtime_identity_mismatch_is_fatal(driver, tmp_path, runtime_identity, key, mutation):
    actual = copy.deepcopy(runtime_identity)
    if mutation == 'absent':
        del actual[key]
    elif mutation == 'wrong_type':
        actual[key] = None
    elif key == 'thread_environment':
        actual[key]['OMP_NUM_THREADS'] = '1'
    elif isinstance(actual[key], bool):
        actual[key] = False
    elif isinstance(actual[key], int):
        actual[key] = 1
    else:
        actual[key] = 'different'
    path = tmp_path / 'fault-runtime.json'
    path.write_text(json.dumps(actual))
    with pytest.raises(RuntimeError, match=key):
        driver.require_runtime_receipt(path, runtime_identity)


@pytest.mark.parametrize('key', ['OMP_NUM_THREADS', 'OQP_FCI_THREADS', 'OMP_DYNAMIC',
                                'OMP_MAX_ACTIVE_LEVELS', 'OPENBLAS_NUM_THREADS',
                                'MKL_NUM_THREADS'])
def test_child_thread_policy_is_consumed(driver, tmp_path, runtime_identity, key):
    actual = copy.deepcopy(runtime_identity)
    actual['thread_environment'][key] = 'unexpected'
    path = tmp_path / 'pytest-runtime.json'
    path.write_text(json.dumps(actual))
    with pytest.raises(RuntimeError, match='thread_environment'):
        driver.require_runtime_receipt(path, runtime_identity)


@pytest.mark.parametrize('payload', [None, '{', '[]'])
def test_missing_or_malformed_runtime_cannot_qualify(driver, tmp_path, runtime_identity, payload):
    path = tmp_path / 'fault-runtime.json'
    if payload is not None:
        path.write_text(payload)
    with pytest.raises(RuntimeError):
        driver.require_runtime_receipt(path, runtime_identity)


def test_boolean_is_not_a_native_thread_count(driver, tmp_path, runtime_identity):
    runtime_identity['omp_requested'] = 1
    actual = copy.deepcopy(runtime_identity)
    actual['omp_requested'] = True
    path = tmp_path / 'pytest-runtime.json'
    path.write_text(json.dumps(actual))
    with pytest.raises(RuntimeError, match='omp_requested'):
        driver.require_runtime_receipt(path, runtime_identity)
