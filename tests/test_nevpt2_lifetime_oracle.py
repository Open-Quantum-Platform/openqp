"""Synthetic controls of the retained-array oracle; no native qualification."""
import numpy as np
import pytest

from tests.native.nevpt2_retained_outputs import RetainedOutputs


def test_owned_results_survive_dimension_changes_and_errors():
    retained = RetainedOutputs()
    for size in (2, 7, 3, 2, 0):
        result = np.arange(size, dtype=float).reshape(1, size)
        retained.remember((result,))
        # Independent later-call/error scratch cannot alter earlier results.
        scratch = np.full((size + 1, size + 1), np.nan)
        del result, scratch
        retained.assert_unchanged()
    assert len(retained) == 5


@pytest.mark.parametrize('mutation', ['reuse', 'resize', 'error'])
def test_borrowed_results_detect_backing_store_mutation(mutation):
    backing = np.arange(12, dtype=float)
    original = backing[:4].reshape(2, 2)
    retained = RetainedOutputs()
    retained.remember((original,))
    if mutation == 'reuse':
        backing[:4] = np.arange(4) + 20
    elif mutation == 'resize':
        # Simulate a persistent buffer reused by a smaller current state.
        backing[:2] = -4
    else:
        # Simulate error cleanup incorrectly poisoning a borrowed output.
        backing[:4] = np.nan
    with pytest.raises(AssertionError):
        retained.assert_unchanged()


def test_changed_original_shape_cannot_pass_as_same_values():
    original = np.arange(4, dtype=float)
    retained = RetainedOutputs()
    retained.remember((original,))
    original.shape = (2, 2)
    with pytest.raises(AssertionError):
        retained.assert_unchanged()
