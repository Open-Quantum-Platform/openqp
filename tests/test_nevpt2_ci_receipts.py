"""Reject missing, skipped or failing named NEVPT2 CI evidence."""
import importlib.util
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
    names = list(driver.TESTS)
    assert len(driver.require_named_passes(_report(tmp_path, names))) == len(names)


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
