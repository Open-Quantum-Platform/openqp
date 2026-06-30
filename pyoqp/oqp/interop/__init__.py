"""Single public surface for MRSF excited-state analysis & interoperability.

One place to import from; the implementation lives in :mod:`oqp.analysis`,
:mod:`oqp.export`, and the in-tree :mod:`oqp.quantum` (FCIDUMP / second-quantized
Hamiltonian)::

    from oqp.interop import (
        MRSFExcitedStates, nto_excitation, attachment_detachment,   # analysis
        CubeExporter, to_qcschema,                                  # export
        dump_fcidump, from_openqp,                                  # FCIDUMP (oqp.quantum)
        parse_output, parse_pyscf_tddft, compare_results,           # external diff
    )
"""
# external-output parsing + tolerance diff (defined here)
from .parsers import parse_output, parse_pyscf_tddft, parse_oqp
from .compare import compare_results, format_table

# excited-state analysis
from oqp.analysis import (
    MRSFExcitedStates, nto_excitation, nto_transition, attachment_detachment,
    participation_ratio, tozer_lambda, fragment_ct_matrix, AOBasis, make_box_grid,
)
# export formats (cubes, QCSchema, FCIDUMP)
from oqp.export import (
    CubeExporter, to_qcschema, validate_qcschema, dump_fcidump, verify_fcidump_fci,
)
# second-quantized Hamiltonian / FCIDUMP implementation (reused, not duplicated)
from oqp.quantum import from_openqp, MolecularHamiltonian, write_fcidump, read_fcidump

__all__ = [
    # external diff
    "parse_output", "parse_pyscf_tddft", "parse_oqp", "compare_results", "format_table",
    # analysis
    "MRSFExcitedStates", "nto_excitation", "nto_transition", "attachment_detachment",
    "participation_ratio", "tozer_lambda", "fragment_ct_matrix", "AOBasis", "make_box_grid",
    # export
    "CubeExporter", "to_qcschema", "validate_qcschema", "dump_fcidump", "verify_fcidump_fci",
    # FCIDUMP / Hamiltonian
    "from_openqp", "MolecularHamiltonian", "write_fcidump", "read_fcidump",
]
