"""The ladder result must not depend on how it is blocked or how wide it runs.

The ladder blocks over a virtual index and hands blocks to threads. Both the
block size and the thread count are tuning knobs -- the thread count is also
reduced automatically when memory is short -- so neither may change the
answer. A wrong leading dimension in the block contraction is invisible until
the last block is short, which happens for some block sizes and not others,
so it shows up as an energy that depends on thread count.

Skipped unless a built OpenQP is importable.
"""

import os
import re
import subprocess
import sys

import pytest


def _have_native_oqp():
    try:
        import oqp
    except Exception:
        return False
    return hasattr(oqp, "ccsd_t_energy")


pytestmark = pytest.mark.skipif(
    not _have_native_oqp(), reason="native OpenQP library not importable"
)

INPUT = """[input]
system=
 8   0.000000000   0.000000000  -0.041061554
 1  -0.533194329   0.533194329  -0.614469223
 1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=cc-pvdz
method=ccsd(t)
functional=

[guess]
type=huckel

[scf]
type=rhf
multiplicity=1
conv=1.0e-10

[cc]
nfzc=1
conv=1.0e-9
"""


def _energy(tmp_path, name, env_extra):
    inp = tmp_path / (name + ".inp")
    inp.write_text(INPUT)
    log = inp.with_suffix(".log")
    if log.exists():
        log.unlink()
    env = dict(os.environ, **env_extra)
    proc = subprocess.run(
        [sys.executable, "-m", "oqp.pyoqp", str(inp)],
        capture_output=True, cwd=str(tmp_path), env=env, timeout=1800,
    )
    assert proc.returncode == 0, proc.stderr.decode(errors="replace")[-2000:]
    text = log.read_text(errors="replace")
    m = re.search(r"E\(CCSD\(T\), correlation\)\s*=\s*(-?[\d.]+)", text)
    assert m, text[-2000:]
    return float(m.group(1))


def test_energy_is_independent_of_the_ladder_block_size(tmp_path):
    """Block sizes chosen so the last block is short as well as exact."""
    ref = _energy(tmp_path, "b_ref", {"OQP_CC_LADDER_BLOCK": "1",
                                      "OMP_NUM_THREADS": "1"})
    for bs in ("2", "3", "5", "7", "19"):
        got = _energy(tmp_path, "b_%s" % bs,
                      {"OQP_CC_LADDER_BLOCK": bs, "OMP_NUM_THREADS": "1"})
        assert abs(got - ref) < 1.0e-9, (bs, got, ref)


def test_energy_is_independent_of_thread_count(tmp_path):
    ref = _energy(tmp_path, "t_ref", {"OMP_NUM_THREADS": "1"})
    for nt in ("2", "3", "8"):
        got = _energy(tmp_path, "t_%s" % nt, {"OMP_NUM_THREADS": nt})
        assert abs(got - ref) < 1.0e-9, (nt, got, ref)
