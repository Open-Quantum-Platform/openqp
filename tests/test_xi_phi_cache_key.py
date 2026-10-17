"""The SCF Phi cache must not replay D^alpha chi built for other xi settings.

The collocation cache persists across calculations in one process and is
keyed by the geometry, the grid and, for xi^alpha, the ingredient settings.
Those settings used to be XORed into the 64-bit key, so (alpha=0.5,
cutoff=1.0) and (alpha=1.0, cutoff=0.5) with the same auto p produced the same
key and the second calculation replayed the first one's fractional-derivative
vectors.  Here, in ONE process with xc_phi_cache=on, run A=(0.5, 1.0) and then
B=(1.0, 0.5); B must equal a fresh process running B without the cache.
Skipped without a compiled OpenQP runtime.
"""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

# C1 water (no petite-list reduction, which would disable the cache)
INPUT = """[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.553194329  -0.614469223
   1   0.513194329  -0.533194329  -0.624469223
charge=0
runtype=energy
basis=6-31g
method=hf
functional=m06-2x
d4=False
[guess]
type=huckel
[scf]
multiplicity=1
type=rhf
conv=1.0e-9
xc_phi_cache={cache}
[dftgrid]
rad_npts=48
ang_npts=110
pruned=
xi_mode=1
xi_alpha={alpha}
xi_cutoff={cutoff}
"""

SCRIPT = r"""
import json, sys
from pathlib import Path
from oqp.pyoqp import Runner
out = []
for tag, alpha, cutoff, cache in json.loads(sys.argv[1]):
    inp = Path(f"{tag}.inp")
    inp.write_text(sys.argv[2].format(alpha=alpha, cutoff=cutoff, cache=cache))
    r = Runner(project=tag, input_file=str(inp), log=f"{tag}.log")
    r.run()
    out.append(float(r.mol.get_scf_energy()))
print("ENERGIES " + json.dumps(out))
"""


def _runtime_available():
    try:
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


def _run(jobs):
    with tempfile.TemporaryDirectory(prefix="oqp_xi_phi_") as tmp:
        proc = subprocess.run(
            [sys.executable, "-c", SCRIPT, json.dumps(jobs), INPUT],
            cwd=tmp, capture_output=True, text=True, timeout=900,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
        )
    line = [l for l in proc.stdout.splitlines() if l.startswith("ENERGIES ")]
    if proc.returncode != 0 or not line:
        raise AssertionError(f"rc={proc.returncode}\n{proc.stdout[-3000:]}\n{proc.stderr[-3000:]}")
    return json.loads(line[-1][len("ENERGIES "):])


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class XiPhiCacheKeyTest(unittest.TestCase):
    def test_cached_run_does_not_replay_other_xi_settings(self):
        e_a, e_b_cached = _run([["a", 0.5, 1.0, "on"], ["b", 1.0, 0.5, "on"]])
        (e_b_fresh,) = _run([["b", 1.0, 0.5, "off"]])
        self.assertNotAlmostEqual(e_a, e_b_fresh, places=6)  # settings really differ
        self.assertAlmostEqual(e_b_cached, e_b_fresh, places=8)


if __name__ == "__main__":
    unittest.main()
