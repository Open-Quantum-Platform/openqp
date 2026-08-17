"""Regression guard for the MRSF Davidson trial-vector diagonal on the ROHF path.

Background
----------
``mrinivec`` builds ``xm``, which both orders the initial trial vectors and
preconditions the Davidson residuals.  What ``xm`` has to approximate is the
diagonal of the response operator, ``A(ij,ij)`` -- not the one-electron part of
it.  On the ROHF path ``tdhf_mrsf_energy`` passes the ROHF canonical
eigenvalues as BOTH ``ea`` and ``eb``.  Those eigenvalues are the
Guest-Saunders spin average, ``eps(p) = 0.5*(fa(p,p) + fb(p,p))``; fitting
``eps(j) - eps(i)`` to an additive ``a(i) + b(j)`` form over all 53 (H2O) and
134 (CH2O) non-open-open amplitudes reproduces that identity to 1e-10.

Issue #328 proposed replacing this with the alpha/beta Fock diagonals, i.e.
``fb(j,j) - fa(i,i)``, which is indeed the exact one-electron diagonal of
``mrsfesum``.  Measurement rejects the substitution:

* The full diagonal also carries the spin-flip exchange ``-c_H*(ij|ji)``
  (JCP 149, 104101, Eq. 2.25), which ``xm`` omits on every path, and which is
  O(1 Eh) here because hole and particle both sit on or next to the two SOMOs.
  The Guest-Saunders subtraction is a surrogate for it.
* Applying the production sigma to unit vectors gives the exact ``A(ij,ij)``.
  Over 45 of the 54 amplitudes of H2O/6-31G (triplet ROHF reference, triplet
  target) the shipped diagonal is closer to it on 45 amplitudes out of 45:
  MRSF-TDHF MAE 0.315 vs 0.542 Eh, MRSF/BHHLYP MAE 0.142 vs 0.272 Eh.
* At the folded open-open slot the shipped value is exact for pure HF.  The
  measured ``A(OO,OO)`` is 0.0000000000 (HF) and 0.0317175663 (BHHLYP) against
  ``xm`` = 0 and a one-electron value of 0.609 / 0.337.  The "identically zero"
  open-open entry is the right answer, not an artefact of ``ea == eb``.
* Substituting the one-electron diagonal loses the 2.039974 eV triplet of
  ``examples/MRSF-TDDFT/CH2O_MRSFTDDFT_SYMMETRY_BLOCK_COVERAGE`` (T1 is then
  reported as 4.782 eV, converged and wrong) and stops
  ``h2o_rohf_mrsf-t_6-31g_{bhhlyp,cam-b3lyp}`` converging -- they reach
  5.9e-08 / 8.5e-08 against a 1e-08 threshold and then exit on
  ``nvec = mxvec``, which the auto-restart cannot rescue because
  ``mxvec = xvec_dim - 3`` is capped by the space dimension itself.

Those behavioural consequences are already covered by the shipped examples.
This file guards the source shape, so the one-line substitution cannot be
reapplied silently -- and, if the call site is ever legitimately restructured,
forces whoever does it to read the measurement first.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENERGY_SRC = ROOT / "source" / "modules" / "tdhf_mrsf_energy.F90"
LIB_SRC = ROOT / "source" / "tdhf_mrsf_lib.F90"

CALL_RE = re.compile(
    r"call\s+mrinivec\s*\(\s*infos\s*,\s*(\w+)\s*,\s*(\w+)\s*,", re.IGNORECASE
)


def _mrinivec_calls():
    return CALL_RE.findall(ENERGY_SRC.read_text())


class MrsfRohfDavidsonDiagonalTests(unittest.TestCase):
    def test_two_mrinivec_call_sites(self):
        """One ROHF call and one UMRSF call, nothing else."""
        self.assertEqual(
            len(_mrinivec_calls()), 2,
            "expected exactly two mrinivec call sites (ROHF and UMRSF); the "
            "guard below assumes that layout",
        )

    def test_rohf_path_uses_the_spin_averaged_diagonal(self):
        """The non-UMRSF call must pass one array as both ea and eb.

        Passing fa/fb there is the change issue #328 asked for; it degrades the
        diagonal (45/45 amplitudes worse), loses a physical CH2O root and stops
        two shipped decks converging.  See this module's docstring.
        """
        calls = _mrinivec_calls()
        rohf = [(ea, eb) for ea, eb in calls if ea == eb]
        self.assertEqual(
            len(rohf), 1,
            "the ROHF MRSF path must pass the Guest-Saunders spin average as "
            "both ea and eb; found %r" % (calls,),
        )

    def test_umrsf_path_keeps_spin_resolved_diagonals(self):
        """UMRSF has a genuine UHF reference, so it keeps alpha/beta."""
        calls = _mrinivec_calls()
        split = [(ea, eb) for ea, eb in calls if ea != eb]
        self.assertEqual(
            len(split), 1,
            "the UMRSF path must keep the spin-resolved alpha/beta diagonals; "
            "found %r" % (calls,),
        )

    def test_umrsf_branch_still_fills_both_arrays_from_fa_fb(self):
        """UMRSF fills its work arrays from the alpha/beta Fock diagonals.

        This is what makes the two paths genuinely different, and it must not
        be extended to the ROHF path.
        """
        text = ENERGY_SRC.read_text()
        self.assertRegex(
            text,
            r"if\s*\(\s*umrsf\s*\)\s*then\s*\n\s*do\s+i\s*=\s*1\s*,\s*nbf\s*\n"
            r"\s*mo_energy_work_a\(i\)\s*=\s*fa\(i,i\)\s*\n"
            r"\s*mo_energy_work_b\(i\)\s*=\s*fb\(i,i\)",
            "the fa/fb fill must stay guarded by 'if (umrsf)' alone",
        )

    def test_measurement_is_recorded_next_to_the_code(self):
        """The reasoning must travel with the code, not only with the issue."""
        energy = ENERGY_SRC.read_text()
        lib = LIB_SRC.read_text()
        for needle in ("#328", "Guest-Saunders", "A(ij,ij)"):
            self.assertIn(needle, energy,
                          "call-site comment lost its reference to %r" % needle)
        self.assertIn("#328", lib,
                      "mrinivec lost its pointer to the measurement")


if __name__ == "__main__":
    unittest.main()
