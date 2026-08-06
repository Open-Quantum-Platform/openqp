"""The non-abelian 'full' integral tier must decline to run under DFT.

``use_integral_symmetry=full`` promotes the ERI reduction to the full point
group, while the XC grid deliberately stays on the abelian sign operations --
Lebedev angular grids are invariant under the axis-aligned octahedral
operations but not under C3/C6, so a full-group grid reduction would be
inexact. The two halves therefore reduce over different groups, which costs
nothing for HF (no grid) and is a measured error for DFT:
``benzene_full_dft`` is out by 3.14e-04 against the tier's own 5e-7 tolerance
(docs/planned/integral-symmetry.md, "Still open").

This test exists because the combination only became REACHABLE when per-shell
purity started being exported. Before that, staging bailed with
``skipped_basis_mismatch`` on every cc-pVXZ, def2 and 6-311G input -- i.e. on
essentially every basis DFT is run in -- so the defect was cornered by an
unrelated bug. Widening engagement to spherical bases without this guard would
have turned it into a live way to get a wrong energy.

The decision is a pure function of the configuration, so it is tested directly
rather than through a live SCF: a run-based test would need a non-abelian
molecule, a converged reference and the reorientation allow-list, and would
still only cover one tier/functional combination.
"""

import unittest

from test_symmetry_metadata import load_molecule_module


def _molecule(config):
    molecule_module = load_molecule_module()
    molecule = molecule_module.Molecule.__new__(molecule_module.Molecule)
    molecule.config = config
    molecule.symmetry_metadata = {}
    return molecule


class FullGroupTierDeclinesForDFTTests(unittest.TestCase):
    def test_hf_is_allowed_on_the_full_tier(self):
        """No functional, no grid, nothing to disagree about."""
        molecule = _molecule({'input': {'method': 'hf', 'functional': ''}})
        self.assertIsNone(molecule._full_group_decline_reason())

    def test_a_functional_declines_the_full_tier(self):
        molecule = _molecule({'input': {'method': 'hf', 'functional': 'bhhlyp'}})
        self.assertEqual(molecule._full_group_decline_reason(),
                         'dft_functional')

    def test_whitespace_only_functional_is_not_a_functional(self):
        """`functional=` with trailing spaces is how an unset key often reads."""
        molecule = _molecule({'input': {'method': 'hf', 'functional': '   '}})
        self.assertIsNone(molecule._full_group_decline_reason())

    def test_missing_input_section_does_not_raise(self):
        """The guard runs inside staging; it must never be the thing that
        aborts a run."""
        self.assertIsNone(_molecule({})._full_group_decline_reason())


if __name__ == '__main__':
    unittest.main()
