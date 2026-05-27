"""Small solvent-coupling helpers used by staged PCM scaffolding."""

from __future__ import annotations


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
    return [-0.5 * float(value) for value in q_cav]
