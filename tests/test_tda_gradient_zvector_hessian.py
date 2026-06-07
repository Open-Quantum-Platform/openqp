"""Regression guard for the TDA (Tamm-Dancoff) excited-state gradient.

Background
----------
The excited-state gradient solves a Z-vector (orbital-relaxation) equation
whose operator is the *ground-state* orbital Hessian (A+B).  That Hessian is
identical for full RPA and for the Tamm-Dancoff approximation - TDA only
changes the excitation vectors and the right-hand side, never the relaxation
operator.

A previous version built the Z-vector operator with ``tamm_dancoff=tda``.
For a TDA run that makes ``int2_td_data_t`` accumulate the A matrix into its
``%amb`` buffer, while ``compute_apbx`` only ever reads ``%apb``.  The
two-electron part of the operator was therefore silently dropped, the
Z-vector collapsed to its uncoupled (diagonal) solution, and the TDA analytic
gradient disagreed with finite differences by a large factor (~3x on the
dominant component) even though the TDA *excitation energy* stayed correct.

Validation at the time of the fix (H2O/6-31G*, central FD step 0.005 bohr):
analytic-vs-FD ``max|diff|`` dropped from ~1.2e-1 to <1e-5 for TDA, and the
fixed TDA gradient matched the GAMESS reference for B3LYP/root-3.  Full RPA
gradients were unaffected.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ZVEC_SRC = ROOT / "source" / "modules" / "tdhf_z_vector.F90"


class TdaZVectorHessianRegressionTests(unittest.TestCase):
    def _zvector_source(self):
        return ZVEC_SRC.read_text()

    def test_cg_operator_uses_ground_state_apb_hessian(self):
        """The Z-vector CG operator must build (A+B), not the TDA A matrix.

        ``compute_apbx`` reads ``int2_data%apb``; the operator must therefore
        be constructed with ``tamm_dancoff=.false.`` so the (A+B) contribution
        lands in ``%apb`` for both TDA and RPA targets.
        """
        src = self._zvector_source()

        # Locate the int2_td_data_t construction that feeds the CG operator
        # (the one with int_apb=.true., int_amb=.false.).
        m = re.search(
            r"int2_td_data_t\(d2\s*=\s*pa,.*?\)",
            src,
            re.S | re.I,
        )
        self.assertIsNotNone(
            m,
            "Could not find the Z-vector CG operator construction "
            "(int2_td_data_t(d2=pa, ...)) in tdhf_z_vector.F90.",
        )
        block = m.group(0)

        self.assertRegex(
            block,
            r"int_apb\s*=\s*\.true\.",
            "The Z-vector CG operator must request the (A+B) combination "
            "(int_apb=.true.).",
        )
        self.assertRegex(
            block,
            r"tamm_dancoff\s*=\s*\.false\.",
            "The Z-vector orbital Hessian is the ground-state (A+B) operator "
            "for both TDA and RPA. Building it with tamm_dancoff=tda routes "
            "the A matrix into %amb, which compute_apbx never reads, breaking "
            "the TDA analytic gradient.",
        )

    def test_cg_operator_reads_apb_buffer(self):
        """compute_apbx must read the %apb buffer the operator writes into."""
        src = self._zvector_source()
        # Within compute_apbx, the operator result is pulled from %apb.
        self.assertRegex(
            src,
            r"apb\s*=>\s*int2_data%apb",
            "compute_apbx is expected to read int2_data%apb; if this changes, "
            "the tamm_dancoff wiring of the Z-vector operator must be "
            "revisited (see test_cg_operator_uses_ground_state_apb_hessian).",
        )


if __name__ == "__main__":
    unittest.main()
