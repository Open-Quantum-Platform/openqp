"""Small solvent-coupling helpers used by staged PCM scaffolding."""

from __future__ import annotations


def _as_float_list(values, *, name: str):
    return [float(value) for value in values]


def reference_scf_total_density(density_blocks):
    """Return the packed total AO density for reference-SCF PCM prototypes.

    OpenQP's first PCM target uses the RHF/ROHF reference density to build the
    ddX cavity potential.  RHF-style inputs have one packed density block, while
    ROHF/UHF-style storage has separate alpha/beta blocks that must be summed
    before evaluating the scalar electrostatic solvent response.
    """
    blocks = [_as_float_list(block, name="density_block") for block in density_blocks]
    if not blocks:
        raise ValueError("reference_scf_total_density requires at least one density block")
    if len(blocks) > 2:
        raise ValueError("reference_scf_total_density accepts one or two density blocks")
    packed_length = len(blocks[0])
    if packed_length == 0:
        raise ValueError("packed density blocks must not be empty")
    if any(len(block) != packed_length for block in blocks):
        raise ValueError("all density blocks must have the same packed length")
    return [sum(block[index] for block in blocks) for index in range(packed_length)]


def provisional_ddx_external_charges(q_cav, *, allow_provisional: bool = False):
    """Return candidate OpenQP external-charge weights from ddX ``q_cav``.

    The current ddPCM finite-difference smoke suggests ``chg = -0.5*q_cav``
    for the future ``external_charge_potential`` AO-matrix seam.  The
    sign/scale is still provisional, so callers must opt in explicitly until
    it is cross-checked against PySCF/ddX/reference-package data.
    """
    if not allow_provisional:
        raise ValueError(
            "ddX q_cav sign/scale is provisional; pass allow_provisional=True "
            "only in guarded validation/prototype code"
        )
    return [-0.5 * value for value in _as_float_list(q_cav, name="q_cav")]


def provisional_ddx_reaction_field_inputs(q_cav, cavity_xyz, *, allow_provisional: bool = False):
    """Validate candidate ddX cavity data for the future AO reaction-field seam.

    ``cavity_xyz`` is the flat ``(3, ncav)`` coordinate buffer copied from ddX,
    and ``q_cav`` is the projected cavity quantity whose sign/scale remains
    provisional.  This helper prepares only guarded prototype inputs; it does
    not enable runtime PCM coupling.
    """
    q_values = _as_float_list(q_cav, name="q_cav")
    xyz_values = _as_float_list(cavity_xyz, name="cavity_xyz")
    expected_xyz = 3 * len(q_values)
    if not q_values:
        raise ValueError("q_cav must contain at least one cavity value")
    if len(xyz_values) != expected_xyz:
        raise ValueError(
            "cavity_xyz must contain 3 * len(q_cav) values "
            f"({expected_xyz} expected, got {len(xyz_values)})"
        )
    charges = provisional_ddx_external_charges(q_values, allow_provisional=allow_provisional)
    return {
        "ncav": len(q_values),
        "cavity_xyz": xyz_values,
        "charges": charges,
        "x": xyz_values[0::3],
        "y": xyz_values[1::3],
        "z": xyz_values[2::3],
        "chg": charges,
    }
