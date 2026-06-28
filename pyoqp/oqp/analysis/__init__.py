"""MRSF excited-state analysis (NTOs, attachment/detachment, descriptors).

Built on the MRSF 1-TDM/1-RDM exposed by the ``misc-excited-analysis`` patch
and validated at GATE 2 (transition-dipole reconstruction to ~1e-15)."""
from .transition_density import MRSFExcitedStates
from .nto import nto_excitation, nto_transition
from .density_diff import attachment_detachment
from .descriptors import participation_ratio, tozer_lambda, fragment_ct_matrix
from .gto_grid import AOBasis, make_box_grid

__all__ = [
    "MRSFExcitedStates",
    "nto_excitation", "nto_transition",
    "attachment_detachment",
    "participation_ratio", "tozer_lambda", "fragment_ct_matrix",
    "AOBasis", "make_box_grid",
]
