"""Exercise both MP2 density allocation failures in Linux subprocesses.

Inputs and the dense output are allocated and touched before limiting address
space. A successful test requires the allocation site's diagnostic, rather
than a generic out-of-memory exit. Run after installing the native package.
"""
import json
import os
from pathlib import Path
import subprocess
import sys


def child(site):
    import resource
    import numpy as np
    import oqp

    n, no = 40, 20
    c = np.asfortranarray(np.eye(n))
    l2 = np.ones((no, no, n-no, n-no), order="F")
    g = np.ones((n, n, n, n), order="F")
    # Resolve native support and prepare pointer arguments before the limit.
    native = oqp.mp2_gradient_density_test
    pointers = [oqp.ffi.cast("double *", a.ctypes.data) for a in (c, l2, g)]
    current = int(Path("/proc/self/statm").read_text().split()[0]) * os.sysconf("SC_PAGE_SIZE")
    # Scratch is 896,000 bytes; the later dense copy is 20,480,000 bytes.
    allowance = 128 * 1024 if site == "scratch" else 4 * 1024**2
    ceiling = current + allowance
    print(json.dumps({"site": site, "n": n, "no": no,
                      "current_vms": current, "address_space_limit": ceiling}), flush=True)
    resource.setrlimit(resource.RLIMIT_AS, (ceiling, ceiling))
    native(n, no, *pointers)
    print("UNEXPECTED_NATIVE_SUCCESS", flush=True)


def main():
    if sys.platform != "linux":
        raise RuntimeError("This address-space fault test requires Linux /proc and RLIMIT_AS")
    if len(sys.argv) == 3 and sys.argv[1] == "child":
        child(sys.argv[2])
        return
    for site in ("scratch", "symmetrization"):
        env = os.environ.copy()
        env.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
        run = subprocess.run([sys.executable, __file__, "child", site],
                             env=env, capture_output=True, text=True, timeout=120)
        output = run.stdout + run.stderr
        print(output, end="", flush=True)
        expected = f"MP2 gradient density {site} allocation failed"
        assert run.returncode != 0, (site, run.returncode, output)
        assert expected in output, (site, run.returncode, output)
        assert "UNEXPECTED_NATIVE_SUCCESS" not in output, (site, output)
        print(f"ALLOCATION_FAILURE_PASS {site} exit={run.returncode}", flush=True)


if __name__ == "__main__":
    main()
