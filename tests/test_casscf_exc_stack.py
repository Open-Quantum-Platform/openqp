"""Equivalence pin for the Fortran excitation-matrix engine (casscf_exc_stack.F90).

``(E_tu)_{row,col} = <row| E_tu |col>`` over the determinant list.  The engine
replaces the per-element Python enumeration in
``casscf_hessian._excitation_matrices_python``, which stays as the pin and the
fallback.  Both produce integer-valued phases, so agreement here is exact, not
approximate -- these tests assert bit equality rather than a tolerance.
"""
import numpy as np
import pytest

from oqp.library.casscf_hessian import (
    _excitation_matrices,
    _excitation_matrices_python,
    _lib_excitation_stack,
)
from oqp.library.fci import _determinants

pytestmark = pytest.mark.skipif(
    _lib_excitation_stack(3, tuple(_determinants(3, (1, 1)))) is None,
    reason="liboqp casscf_excitation_stack unavailable",
)


@pytest.mark.parametrize(
    "nact,nalpha,nbeta",
    [
        (3, 1, 1),
        (4, 2, 2),
        (5, 3, 2),
        (6, 3, 3),
        (4, 2, 1),   # unequal alpha/beta
        (5, 1, 4),   # strongly unequal, beta-heavy
        (4, 0, 2),   # empty alpha string
        (4, 4, 4),   # closed shell, single determinant
    ],
)
def test_fortran_stack_is_bit_identical_to_the_python_enumeration(nact, nalpha, nbeta):
    dets = tuple(_determinants(nact, (nalpha, nbeta)))
    ref = _excitation_matrices_python(nact, nalpha, nbeta)[1]
    got = _lib_excitation_stack(nact, dets)

    assert got is not None
    assert got.shape == ref.shape
    assert np.array_equal(got, ref)


def test_cached_entry_point_uses_the_engine_and_stays_read_only():
    _excitation_matrices.cache_clear()
    dets, stack = _excitation_matrices(5, 3, 2)
    assert np.array_equal(stack, _excitation_matrices_python(5, 3, 2)[1])
    assert not stack.flags.writeable
    assert dets == tuple(_determinants(5, (3, 2)))


def test_excitation_matrices_reproduce_the_particle_number():
    """sum_t (E_tt) is the electron-number operator: diagonal, equal to nelec.

    An index or phase slip that still produced a plausible-looking stack would
    break this, since it fixes the diagonal absolutely rather than relatively.
    """
    nact, nalpha, nbeta = 5, 3, 2
    dets = tuple(_determinants(nact, (nalpha, nbeta)))
    stack = _lib_excitation_stack(nact, dets)
    number = sum(stack[t, t] for t in range(nact))
    expected = np.eye(len(dets)) * (nalpha + nbeta)
    assert np.array_equal(number, expected)


def test_excitation_matrices_are_mutually_transpose():
    """<row|E_tu|col> = <col|E_ut|row>: E_ut is the transpose of E_tu."""
    nact, nalpha, nbeta = 4, 2, 2
    dets = tuple(_determinants(nact, (nalpha, nbeta)))
    stack = _lib_excitation_stack(nact, dets)
    for t in range(nact):
        for u in range(nact):
            assert np.array_equal(stack[t, u], stack[u, t].T)
