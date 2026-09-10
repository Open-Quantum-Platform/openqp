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

How the guard identifies the two calls
--------------------------------------
The two ``mrinivec`` call sites cannot be told apart by their arguments alone.
A source file whose ``if (.not. umrsf)`` predicate has been inverted, or whose
two calls have been exchanged between the branches, still contains exactly one
call with equal arguments and one call with unequal arguments -- while ROHF
receives the spin-resolved diagonal and UMRSF receives the averaged one, which
is exactly the regression this file exists to prevent.

So the source is parsed instead: comments are removed, ``&`` continuations are
joined, the ``if``/``else if``/``else``/``end if`` tree is walked, and each call
is reported together with the value ``umrsf`` is *known* to have where that call
stands.  Anything the walk cannot pin down (a call outside any ``umrsf``
branch, a compound predicate, a contradictory nesting, an unbalanced
construct) is reported as a problem rather than silently accepted.
``test_inverting_the_branch_predicate_is_rejected`` and
``test_exchanging_the_two_calls_is_rejected`` apply those two mutations to the
shipped source and require the guard to reject each of them, and
``test_inverting_and_exchanging_together_is_accepted`` applies both -- which
restores the original meaning -- and requires the guard to accept it, so the
guard is demonstrably reading control flow rather than line order.
"""

import re
import unittest
from collections import namedtuple
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENERGY_SRC = ROOT / "source" / "modules" / "tdhf_mrsf_energy.F90"
LIB_SRC = ROOT / "source" / "tdhf_mrsf_lib.F90"

IF_THEN_RE = re.compile(r"^(?:\w+\s*:\s*)?if\s*\((?P<cond>.+)\)\s*then$", re.I)
ELSE_IF_RE = re.compile(r"^else\s*if\s*\((?P<cond>.+)\)\s*then(?:\s+\w+)?$", re.I)
ELSE_RE = re.compile(r"^else(?:\s+\w+)?$", re.I)
END_IF_RE = re.compile(r"^end\s*if(?:\s+\w+)?$", re.I)
MRINIVEC_RE = re.compile(r"^call\s+mrinivec\s*\((?P<args>.*)\)$", re.I)


class _Unknown:
    """``umrsf`` has a value here, but this parser refuses to guess which."""

    def __repr__(self):
        return "unknown"


UNKNOWN = _Unknown()

# One ``mrinivec`` call site: the physical lines it occupies, its actual
# arguments, the value ``umrsf`` has there, and the line of the ``if`` /
# ``else if`` / ``else`` that fixes that value.
Call = namedtuple("Call", "start end args umrsf branch_line")


def _strip_comment(line):
    """Drop a trailing Fortran ``!`` comment, respecting quoted text."""
    out = []
    quote = None
    for char in line:
        if quote is not None:
            out.append(char)
            if char == quote:
                quote = None
        elif char in "'\"":
            quote = char
            out.append(char)
        elif char == "!":
            break
        else:
            out.append(char)
    return "".join(out).strip()


def _logical_lines(text):
    """Yield ``(first_line, last_line, code)``, ``&`` continuations joined."""
    lines = []
    buf = ""
    start = None
    for number, raw in enumerate(text.splitlines(), 1):
        code = _strip_comment(raw)
        if not code:
            # A blank or comment-only line may sit inside a continuation.
            continue
        if buf:
            code = code.lstrip("&").strip()
        else:
            start = number
        if code.endswith("&"):
            buf += code[:-1].rstrip() + " "
            continue
        lines.append((start, number, (buf + code).strip()))
        buf = ""
        start = None
    if buf:
        raise ValueError("file ends inside a continued statement at line %d" % start)
    return lines


def _umrsf_says(condition, taken):
    """What ``condition`` -- or its negation, when ``taken`` is False -- fixes.

    Returns True/False when ``umrsf`` is pinned down, None when the condition
    has nothing to do with ``umrsf``, and UNKNOWN when it mentions ``umrsf`` in
    a form this parser will not reason about.
    """
    text = re.sub(r"\s+", "", condition).lower()
    if text == "umrsf":
        return taken
    if text == ".not.umrsf":
        return not taken
    if "umrsf" in text:
        return UNKNOWN
    return None


def _combine(values):
    """Fold the branch constraints into one value for ``umrsf``."""
    fixed = [value for value in values if value is not None]
    if not fixed:
        return None
    if any(value is UNKNOWN for value in fixed):
        return UNKNOWN
    if len(set(fixed)) != 1:
        return UNKNOWN  # contradictory nesting; refuse to pick a side
    return fixed[0]


class _IfConstruct:
    """One open ``if`` construct and what its current branch says about umrsf."""

    def __init__(self, line, condition):
        self.branch_line = line
        self.conditions = [condition]
        self.umrsf = _umrsf_says(condition, True)

    def else_if(self, line, condition):
        self.branch_line = line
        self.umrsf = _combine(
            [_umrsf_says(prior, False) for prior in self.conditions]
            + [_umrsf_says(condition, True)]
        )
        self.conditions.append(condition)

    def else_(self, line):
        self.branch_line = line
        self.umrsf = _combine(
            [_umrsf_says(prior, False) for prior in self.conditions]
        )


def _enclosing_umrsf(stack):
    """The innermost branch that fixes ``umrsf``, as ``(value, line)``."""
    for frame in reversed(stack):
        if frame.umrsf is not None:
            return frame.umrsf, frame.branch_line
    return None, None


def mrinivec_calls(text):
    """Every ``mrinivec`` call in ``text``, tagged with its enclosing branch."""
    calls = []
    stack = []
    for start, end, code in _logical_lines(text):
        match = ELSE_IF_RE.match(code)
        if match:
            if not stack:
                raise ValueError("'else if' outside an if construct at line %d" % start)
            stack[-1].else_if(start, match.group("cond"))
            continue
        match = IF_THEN_RE.match(code)
        if match:
            stack.append(_IfConstruct(start, match.group("cond")))
            continue
        if ELSE_RE.match(code):
            if not stack:
                raise ValueError("'else' outside an if construct at line %d" % start)
            stack[-1].else_(start)
            continue
        if END_IF_RE.match(code):
            if not stack:
                raise ValueError("'end if' without a matching 'if' at line %d" % start)
            stack.pop()
            continue
        match = MRINIVEC_RE.match(code)
        if match:
            umrsf, branch_line = _enclosing_umrsf(stack)
            calls.append(
                Call(start, end, _split_args(match.group("args")), umrsf, branch_line)
            )
    if stack:
        raise ValueError(
            "if construct opened at line %d is never closed" % stack[-1].branch_line
        )
    return calls


def _split_args(text):
    """Split an argument list on its top-level commas."""
    args = []
    current = ""
    depth = 0
    for char in text:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        if char == "," and depth == 0:
            args.append(current.strip())
            current = ""
        else:
            current += char
    args.append(current.strip())
    return args


def diagnose(text):
    """Problems with the ``mrinivec`` branch layout; empty means the shape is right."""
    try:
        calls = mrinivec_calls(text)
    except ValueError as exc:
        return ["could not parse the if-construct tree: %s" % exc]

    problems = []
    if len(calls) != 2:
        problems.append(
            "expected exactly two mrinivec call sites (ROHF and UMRSF), found %d"
            % len(calls)
        )
    for call in calls:
        if call.umrsf is None:
            problems.append(
                "the mrinivec call at line %d is not inside any branch that fixes "
                "umrsf, so its role cannot be established" % call.start
            )
        elif call.umrsf is UNKNOWN:
            problems.append(
                "the branch at line %d guarding the mrinivec call at line %d does "
                "not fix umrsf to a single value" % (call.branch_line, call.start)
            )
        elif len(call.args) < 3:
            problems.append(
                "the mrinivec call at line %d has too few arguments to carry ea and "
                "eb: %r" % (call.start, call.args)
            )

    usable = [call for call in calls if call.umrsf in (True, False) and len(call.args) >= 3]
    rohf = [call for call in usable if call.umrsf is False]
    umrsf = [call for call in usable if call.umrsf is True]

    if len(rohf) != 1:
        problems.append(
            "expected exactly one mrinivec call under '.not. umrsf' (the ROHF "
            "path), found %d" % len(rohf)
        )
    elif rohf[0].args[1] != rohf[0].args[2]:
        problems.append(
            "the ROHF call at line %d (branch at line %d) passes ea=%s, eb=%s; the "
            "Guest-Saunders spin average has to go in as BOTH -- see this module's "
            "docstring and issue #328"
            % (rohf[0].start, rohf[0].branch_line, rohf[0].args[1], rohf[0].args[2])
        )

    if len(umrsf) != 1:
        problems.append(
            "expected exactly one mrinivec call under 'umrsf' (the UMRSF path), "
            "found %d" % len(umrsf)
        )
    elif umrsf[0].args[1] == umrsf[0].args[2]:
        problems.append(
            "the UMRSF call at line %d (branch at line %d) passes %s as both ea and "
            "eb; a UHF reference has genuine alpha/beta diagonals and must keep them"
            % (umrsf[0].start, umrsf[0].branch_line, umrsf[0].args[1])
        )

    return problems


def _invert_the_branch_predicate(text):
    """Rewrite the ROHF call's guard from ``.not. umrsf`` to ``umrsf``."""
    rohf = [call for call in mrinivec_calls(text) if call.umrsf is False]
    if len(rohf) != 1:
        raise ValueError("cannot locate the single ROHF call to mutate")
    lines = text.splitlines(keepends=True)
    index = rohf[0].branch_line - 1
    mutated, count = re.subn(r"\.not\.\s*umrsf", "umrsf", lines[index], flags=re.I)
    if count != 1:
        raise ValueError("no '.not. umrsf' to invert on line %d" % rohf[0].branch_line)
    lines[index] = mutated
    return "".join(lines)


def _exchange_the_two_calls(text):
    """Swap the two ``mrinivec`` call statements between their branches."""
    calls = mrinivec_calls(text)
    if len(calls) != 2:
        raise ValueError("cannot exchange calls: found %d" % len(calls))
    first, second = calls
    lines = text.splitlines(keepends=True)
    return "".join(
        lines[: first.start - 1]
        + lines[second.start - 1 : second.end]
        + lines[first.end : second.start - 1]
        + lines[first.start - 1 : first.end]
        + lines[second.end :]
    )


class MrsfRohfDavidsonDiagonalTests(unittest.TestCase):
    def setUp(self):
        self.energy = ENERGY_SRC.read_text()

    def test_two_mrinivec_call_sites(self):
        """One ROHF call and one UMRSF call, nothing else."""
        self.assertEqual(
            len(mrinivec_calls(self.energy)), 2,
            "expected exactly two mrinivec call sites (ROHF and UMRSF); the "
            "guard below assumes that layout",
        )

    def test_rohf_branch_passes_the_spin_averaged_diagonal(self):
        """The call inside ``.not. umrsf`` passes one array as both ea and eb.

        Identified by its enclosing branch, not by its arguments: passing fa/fb
        there is the change issue #328 asked for; it degrades the diagonal
        (45/45 amplitudes worse), loses a physical CH2O root and stops two
        shipped decks converging.  See this module's docstring.
        """
        rohf = [c for c in mrinivec_calls(self.energy) if c.umrsf is False]
        self.assertEqual(
            len(rohf), 1,
            "expected exactly one mrinivec call under '.not. umrsf'; found %r" % (rohf,),
        )
        ea, eb = rohf[0].args[1], rohf[0].args[2]
        self.assertEqual(
            ea, eb,
            "the ROHF MRSF path (branch at line %d, call at line %d) must pass the "
            "Guest-Saunders spin average as both ea and eb; it passes %s and %s"
            % (rohf[0].branch_line, rohf[0].start, ea, eb),
        )

    def test_umrsf_branch_keeps_spin_resolved_diagonals(self):
        """UMRSF has a genuine UHF reference, so its branch keeps alpha/beta."""
        umrsf = [c for c in mrinivec_calls(self.energy) if c.umrsf is True]
        self.assertEqual(
            len(umrsf), 1,
            "expected exactly one mrinivec call under 'umrsf'; found %r" % (umrsf,),
        )
        ea, eb = umrsf[0].args[1], umrsf[0].args[2]
        self.assertNotEqual(
            ea, eb,
            "the UMRSF path (branch at line %d, call at line %d) must keep the "
            "spin-resolved alpha/beta diagonals; it passes %s twice"
            % (umrsf[0].branch_line, umrsf[0].start, ea),
        )

    def test_every_call_site_is_inside_a_resolved_umrsf_branch(self):
        """Neither call may drift out of -- or into a fuzzier -- umrsf branch."""
        for call in mrinivec_calls(self.energy):
            self.assertIn(
                call.umrsf, (True, False),
                "the mrinivec call at line %d is not guarded by a branch that fixes "
                "umrsf to a single value (got %r); the ROHF and UMRSF diagonals can "
                "no longer be told apart" % (call.start, call.umrsf),
            )

    def test_layouts_that_hide_the_branch_are_rejected(self):
        """What the walk cannot pin down must be reported, not accepted.

        These are the ways a future restructuring could stop the two calls from
        being distinguishable; each has to fail loudly rather than pass by
        default.
        """
        cases = {
            "call outside any umrsf branch": (
                "    call mrinivec(infos, ea, ea, bvec_mo, xm, nvec)\n"
                "    if (umrsf) then\n"
                "      call mrinivec(infos, ea, eb, bvec_mo, xm, nvec)\n"
                "    end if\n"
            ),
            "umrsf hidden in a compound predicate": (
                "    if ((mrst==1 .or. mrst==3) .and. .not. umrsf) then\n"
                "      call mrinivec(infos, ea, ea, bvec_mo, xm, nvec)\n"
                "    else\n"
                "      call mrinivec(infos, ea, eb, bvec_mo, xm, nvec)\n"
                "    end if\n"
            ),
            "unbalanced if construct": (
                "    if (.not. umrsf) then\n"
                "      call mrinivec(infos, ea, ea, bvec_mo, xm, nvec)\n"
            ),
        }
        for name, source in cases.items():
            with self.subTest(layout=name):
                self.assertTrue(
                    diagnose(source),
                    "%s was accepted; the guard must fail closed" % name,
                )

    def test_shipped_source_has_no_layout_problems(self):
        """The whole guard, run over the shipped file, must be silent."""
        self.assertEqual(diagnose(self.energy), [])

    def test_inverting_the_branch_predicate_is_rejected(self):
        """Flipping ``.not. umrsf`` hands ROHF the spin-resolved diagonal.

        The argument-only classification cannot see this: there is still one
        equal pair and one unequal pair.
        """
        mutated = _invert_the_branch_predicate(self.energy)
        self.assertNotEqual(mutated, self.energy, "the mutation changed nothing")
        problems = diagnose(mutated)
        self.assertTrue(
            problems,
            "inverting the branch predicate gives ROHF the alpha/beta diagonal and "
            "the guard accepted it",
        )
        self.assertTrue(
            any("ROHF call" in problem for problem in problems),
            "expected the ROHF call to be reported; got %r" % (problems,),
        )

    def test_exchanging_the_two_calls_is_rejected(self):
        """Swapping the call statements has the same effect and must also fail."""
        mutated = _exchange_the_two_calls(self.energy)
        self.assertNotEqual(mutated, self.energy, "the mutation changed nothing")
        problems = diagnose(mutated)
        self.assertTrue(
            problems,
            "exchanging the two calls gives ROHF the alpha/beta diagonal and the "
            "guard accepted it",
        )
        self.assertTrue(
            any("ROHF call" in problem for problem in problems),
            "expected the ROHF call to be reported; got %r" % (problems,),
        )

    def test_inverting_and_exchanging_together_is_accepted(self):
        """Both mutations at once restore the shipped meaning, so this must pass.

        Positive control: it shows the guard reads the control flow rather than
        the order of the two call statements.
        """
        mutated = _invert_the_branch_predicate(_exchange_the_two_calls(self.energy))
        self.assertNotEqual(mutated, self.energy, "the mutation changed nothing")
        self.assertEqual(
            diagnose(mutated), [],
            "inverting the predicate and exchanging the calls is the shipped "
            "assignment written the other way round; the guard must accept it",
        )

    def test_umrsf_branch_still_fills_both_arrays_from_fa_fb(self):
        """UMRSF fills its work arrays from the alpha/beta Fock diagonals.

        This is what makes the two paths genuinely different, and it must not
        be extended to the ROHF path.
        """
        self.assertRegex(
            self.energy,
            r"if\s*\(\s*umrsf\s*\)\s*then\s*\n\s*do\s+i\s*=\s*1\s*,\s*nbf\s*\n"
            r"\s*mo_energy_work_a\(i\)\s*=\s*fa\(i,i\)\s*\n"
            r"\s*mo_energy_work_b\(i\)\s*=\s*fb\(i,i\)",
            "the fa/fb fill must stay guarded by 'if (umrsf)' alone",
        )

    def test_measurement_is_recorded_next_to_the_code(self):
        """The reasoning must travel with the code, not only with the issue."""
        lib = LIB_SRC.read_text()
        for needle in ("#328", "Guest-Saunders", "A(ij,ij)"):
            self.assertIn(needle, self.energy,
                          "call-site comment lost its reference to %r" % needle)
        self.assertIn("#328", lib,
                      "mrinivec lost its pointer to the measurement")


if __name__ == "__main__":
    unittest.main()
