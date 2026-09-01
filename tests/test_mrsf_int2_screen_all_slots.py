"""Regression guard: the MRSF 2e Fock build must screen against ALL densities.

int2_mrsf_data_t digests seven densities per shell quartet (bo2v, bo1v, bco1,
bco2, o21v, co12, ball), but its Schwarz shell-density bound used to be built
from slot 7 (`ball`) alone.  `ball` is a sum whose internal cancellations can
be small exactly where the mixed spin-pair densities (co12/o21v) are large, so
the driver dropped real contributions to the aco12/ao21v Fock matrices.  The
energy stage was effectively unaffected (multi-trial-vector Davidson builds
give union-like bounds), but the single-density Z-vector RHS build was not:
its spin-pair Lagrangian (mrsfsp) went wrong, and MRSF analytic gradients
disagreed with finite differences by up to ~4e-3 Hartree/Bohr in a
state-dependent way (largest for states with closed->open amplitudes), while
all state energies stayed exact.  MRSF-EKT values drift through the same
relaxed-density machinery.

The screening bound must majorize every contracted density: init_screen now
takes the elementwise max of the shell-density bound over all slots.

Two guards:
  * a source-text check that init_screen builds the bound over all slots;
  * a live check on the original repro (C2v-twisted water, MRSF-TDHF/6-31G,
    S3 = the 4th singlet, the worst case at ~4.3e-3): the analytic gradient
    must match a central finite difference of the state energy along H1-x to
    the FD-truncation level.  With the ball-only bound this component is off
    by 4.3e-3; correct screening gives ~2e-7.

Skipped unless the compiled OpenQP runtime is importable.
"""
import os
import re
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Same molecule and settings as examples/other/h2o_rohf_mrsf-s_6-31g_* but with
# a pure-HF reference (no grid noise) and tight thresholds for a clean FD.
INPUT_TMPL = """[input]
system=
 8   0.000000000   0.000000000  -0.041061554
 1  {h1x:12.9f}   0.533194329  -0.614469223
 1   0.533194329  -0.533194329  -0.614469223
charge=0
method=tdhf
basis=6-31g
runtype={runtype}
functional=
d4=False

[guess]
type=huckel
save_mol=False

[scf]
type=rohf
maxit=50
multiplicity=3
conv=1.0e-11
save_molden=False
{properties}
[tdhf]
type=mrsf
maxit=100
multiplicity=1
conv=1.0e-11
nstate=6
zvconv=1.0e-11
"""

H1X0 = -0.533194329
H_ANG = 5.0e-4                       # FD step, Angstrom
ANGSTROM_PER_BOHR = 0.52917721090299996
TARGET = 3                           # S3 (input state 4): worst pre-fix error

STATE_RE = re.compile(r"^\s+S(\d+)\s+(-\d+\.\d{8,})\s")


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


class MrsfInt2ScreenSourceTests(unittest.TestCase):
    def test_init_screen_bounds_every_density_slot(self):
        """init_screen must fold the shell-density bound over ALL d3 slots."""
        src = (ROOT / "source" / "tdhf_mrsf_lib.F90").read_text()
        blk = re.search(
            r"subroutine int2_mrsf_data_t_init_screen(.*?)end subroutine",
            src, re.S)
        self.assertIsNotNone(blk, "init_screen not found")
        body = blk.group(1)
        self.assertRegex(
            body, r"do\s+c\s*=\s*1\s*,\s*sized",
            "init_screen must loop over every density slot; a bound built "
            "from one slot (e.g. ball) does not majorize co12/o21v and lets "
            "the driver drop real aco12/ao21v contributions.")
        self.assertRegex(
            body, r"max\s*\(\s*this%dsh\s*,",
            "init_screen must keep the elementwise max of the per-slot "
            "shell-density bounds.")


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class MrsfGradientFdScreeningRepro(unittest.TestCase):
    def _run(self, workdir, name, h1x, runtype, properties):
        from oqp.pyoqp import Runner

        inp = Path(workdir) / f"{name}.inp"
        log = Path(workdir) / f"{name}.log"
        inp.write_text(INPUT_TMPL.format(h1x=h1x, runtype=runtype,
                                         properties=properties))
        runner = Runner(project=name, input_file=str(inp), log=str(log))
        runner.run()
        return log

    @staticmethod
    def _state_energies(log):
        """Absolute state energies (Hartree) from the final state table."""
        out = {}
        for line in log.read_text().splitlines():
            m = STATE_RE.match(line)
            if m and float(m.group(2)) <= -70.0:
                out[int(m.group(1))] = float(m.group(2))
        return out

    @staticmethod
    def _final_gradient(log):
        """Last 'PyOQP S<k>' gradient block (atoms x 3, Hartree/Bohr)."""
        grads, cur = [], None
        for line in log.read_text().splitlines():
            if re.search(r"PyOQP S\d+\s*$", line.rstrip()):
                cur = []
                continue
            if cur is not None:
                f = line.split()
                if len(f) == 4 and f[0] in ("O", "H"):
                    cur.append([float(x) for x in f[1:]])
                    if len(cur) == 3:
                        grads.append(cur)
                        cur = None
        return grads[-1]

    def test_s3_analytic_gradient_matches_fd_along_h1x(self):
        with tempfile.TemporaryDirectory(prefix="oqp_mrsf_screen_") as td:
            glog = self._run(td, "grad", H1X0, "grad",
                             "\n[properties]\ngrad=4\n")
            grad = self._final_gradient(glog)
            ep = self._state_energies(
                self._run(td, "ep", H1X0 + H_ANG, "energy", ""))[TARGET]
            em = self._state_energies(
                self._run(td, "em", H1X0 - H_ANG, "energy", ""))[TARGET]

        fd = (ep - em) / (2.0 * H_ANG) * ANGSTROM_PER_BOHR
        diff = abs(fd - grad[1][0])
        self.assertLess(
            diff, 1.0e-4,
            f"MRSF S3 analytic dE/dx(H1) = {grad[1][0]:.8f} vs FD {fd:.8f} "
            f"(|diff| = {diff:.2e}).  A ~4e-3 disagreement here is the "
            "signature of the ball-only 2e screening bound dropping "
            "spin-pair-coupling (aco12) contributions in the Z-vector build.")


if __name__ == "__main__":
    unittest.main()
