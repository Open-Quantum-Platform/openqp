"""Regression guard: MRSF gradients must be FD-exact for range-separated functionals.

With a CAM-type functional grd2_driver runs the MRSF 2e-derivative digest twice --
pass 1 over the full ERI, pass 2 over the erf-attenuated short-range one -- and
re-points coulscale/hfscale/hfscale2/mu for each pass.  It does NOT touch
``spcscale``, which is the MRSF compute type's own field, so the six spin-pair
coupling contractions (co12, o21v, bco1*bo2v, bco2*bo1v) used to be added in BOTH
passes.  The energy adds them in neither: int2_mrsf_data_t_update's ``cur_pass==2``
branch digests slot 7 (ball) alone, so slots 1-6 carry a pass-1 contribution only.

The gradient therefore held a term with no counterpart in the energy, and MRSF
analytic gradients disagreed with central finite differences by up to ~4e-3
Hartree/Bohr for any range-separated functional -- state-dependently, in
proportion to the state's closed->open (spin-pair) character.  Global hybrids
never run pass 2 and were always exact, which is why this hid behind the
committed dtcam-* references.

The check uses stock ``cam-b3lyp`` rather than one of the dtcam-* functionals on
purpose: the defect is a property of the attenuated pass, not of the tuned
parameter sets.  S3 of the C2v-twisted water below moved by 2.7e-4 Hartree/Bohr.

Skipped unless the compiled OpenQP runtime is importable.
"""
import os
import re
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

INPUT_TMPL = """[input]
system=
 8   0.000000000   0.000000000  -0.041061554
 1  {h1x:12.9f}   0.533194329  -0.614469223
 1   0.533194329  -0.533194329  -0.614469223
charge=0
method=tdhf
basis=6-31g
runtype={runtype}
functional=cam-b3lyp
d4=False

[guess]
type=huckel
save_mol=False

[scf]
type=rohf
maxit=50
maxdiis=5
multiplicity=3
conv=1.0e-10
save_molden=False

[dftgrid]
rad_npts=96
ang_npts=302
pruned=
{properties}
[tdhf]
type=mrsf
maxit=50
multiplicity=1
conv=1.0e-10
nstate=6
zvconv=1.0e-10
"""

H1X0 = -0.533194329
H_ANG = 5.0e-4                       # FD step, Angstrom
ANGSTROM_PER_BOHR = 0.52917721090299996
TARGET = 3                           # S3 (input state 4)

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


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class MrsfCamGradientSpcPassTests(unittest.TestCase):
    def _run(self, workdir, name, h1x, runtype, properties):
        from oqp.pyoqp import Runner

        inp = Path(workdir) / f"{name}.inp"
        log = Path(workdir) / f"{name}.log"
        inp.write_text(INPUT_TMPL.format(h1x=h1x, runtype=runtype,
                                         properties=properties))
        Runner(project=name, input_file=str(inp), log=str(log)).run()
        return log

    @staticmethod
    def _state_energies(log):
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

    def test_s3_cam_gradient_matches_fd_along_h1x(self):
        with tempfile.TemporaryDirectory(prefix="oqp_mrsf_cam_") as td:
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
            diff, 1.0e-5,
            f"MRSF/cam-b3lyp S3 analytic dE/dx(H1) = {grad[1][0]:.8f} vs FD "
            f"{fd:.8f} (|diff| = {diff:.2e}).  A ~3e-4 disagreement here is the "
            "signature of the spin-pair-coupling terms being digested in the "
            "attenuated (short-range) CAM pass of the 2e-derivative gradient, "
            "where the response Fock build contributes none.")


if __name__ == "__main__":
    unittest.main()
