"""A failed TD solve must not leave the previous solve's imaginary-root flags.

SinglePoint.excitation() sets Molecule.td_imaginary_states after a successful
RPA solve; get_results() masks the energies of those roots.  When the same
Molecule is reused and the next TD solve raises (or exhausts its retries),
the flags from the earlier run must already be gone.
"""
from types import SimpleNamespace

import pytest

from oqp.library.single_point import SinglePoint


def test_failed_td_solve_clears_previous_imaginary_roots():
    mol = SimpleNamespace(td_imaginary_states=(2,), symmetry_metadata=None)

    def failing_tddft():
        raise RuntimeError("TD solve failed")

    caller = SimpleNamespace(mol=mol, tddft=failing_tddft)
    with pytest.raises(RuntimeError, match="TD solve failed"):
        SinglePoint.excitation(caller, [0.0])
    assert mol.td_imaginary_states == ()
