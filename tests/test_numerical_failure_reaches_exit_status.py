"""A finite-difference run whose displacements failed must not exit 0.

The failure was always DETECTED -- `flags` carries 'failed' and the code acts on
it, deliberately keeping the scratch directory instead of deleting it -- but it
reached only the log, never the exit status. A numerical Hessian in which every
one of the twelve displacements failed printed

    PyOQP step: 12  displacement: 11  failed
    PyOQP: numerical hessian calculations failed

and exited 0. Any CI job, queue script or geometry scan driving OpenQP by exit
code recorded that as a success.

Silent success is worse than a wrong number: a wrong number can be caught
downstream, a false success cannot.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SINGLE_POINT = ROOT / 'pyoqp/oqp/library/single_point.py'


def block_after(source, marker, span=900, strip_comments=True):
    """Source following `marker`, comments removed by default.

    Comment stripping is not cosmetic here: the first version of the assertion
    below matched the phrase 'return None' inside the explanatory COMMENT above
    the fix, so a test written to prove the fix reported the fix missing. The
    same shape -- a structural assertion satisfied or defeated by prose --
    already produced two false passes elsewhere in this branch family.
    """
    idx = source.index(marker)
    block = source[idx:idx + span]
    return re.sub(r'#.*', '', block) if strip_comments else block


class FailedDisplacementsRaise(unittest.TestCase):
    def test_numerical_hessian_raises_rather_than_returning_none(self):
        source = SINGLE_POINT.read_text()
        block = block_after(
            source, "title='PyOQP: numerical hessian calculations failed'")
        self.assertIn('raise RuntimeError(', block)
        self.assertNotIn('return None', block.split('raise RuntimeError(')[0])

    def test_numerical_nac_raises_too(self):
        """The NAC driver has the identical structure and the identical hole."""
        source = SINGLE_POINT.read_text()
        block = block_after(
            source, "title='PyOQP: numerical nac calculations failed'")
        self.assertIn('raise RuntimeError(', block)
        self.assertNotIn('return None', block.split('raise RuntimeError(')[0])

    def test_the_message_says_how_many_failed(self):
        """'some failed' is not actionable; the count and total are."""
        source = SINGLE_POINT.read_text()
        for marker in ('numerical Hessian: %d of %d', 'numerical NAC: %d of %d'):
            with self.subTest(driver=marker):
                self.assertIn(marker, source)


class ChildrenRunTheParentsBuild(unittest.TestCase):
    """Displacement children must not resolve through PATH.

    sys.executable is guaranteed consistent with the running parent; PATH is
    not. Invoking a venv's `openqp` by absolute path without putting that venv
    on PATH -- exactly what a side-by-side comparison of two builds does -- had
    the parent running one OpenQP while every child ran another. Reproduced by
    making the failure appear and disappear with nothing changed but PATH
    order.
    """

    def test_the_external_runner_uses_the_running_interpreter(self):
        source = SINGLE_POINT.read_text()
        body = block_after(source, 'def _run_oqp_external(', span=2000)
        self.assertIn("cmd = [sys.executable, '-m', 'oqp.pyoqp'", body)

    def test_it_does_not_resolve_openqp_through_path(self):
        source = SINGLE_POINT.read_text()
        body = block_after(source, 'def _run_oqp_external(', span=2000)
        code = re.sub(r'#.*', '', body)
        self.assertNotIn("shutil.which('openqp')", code,
                         'PATH is not guaranteed to agree with the parent')


if __name__ == '__main__':
    unittest.main()
