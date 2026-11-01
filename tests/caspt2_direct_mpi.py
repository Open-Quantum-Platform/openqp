"""Two-rank CASPT2 integration, run by the MPI-enabled CI lane.

Usage: mpiexec -n 2 python tests/caspt2_direct_mpi.py OUTPUT_DIRECTORY
Each rank compares both direct engines with its serial dense calculation;
shared logs must contain exactly one engine and denominator summary.
"""
import json
import sys
from pathlib import Path

import numpy as np
from mpi4py import MPI

from test_caspt2_dyall import _run_caspt2


def values(runner):
    return np.concatenate([
        np.asarray(runner.mol.data[key], dtype=float).ravel()
        for key in ("OQP::CASPT2_ENERGIES", "OQP::CASPT2_REFERENCE_ENERGIES",
                    "OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS")
    ])


def main():
    comm = MPI.COMM_WORLD
    assert comm.size == 2, "this integration test requires two MPI ranks"
    out = Path(sys.argv[1]).resolve()
    own = out / f"rank-{comm.rank}"
    own.mkdir(parents=True, exist_ok=True)
    kwargs = dict(system="\nH 0 0 0\nH 0.03 0 0.74\nH 0 0.06 1.53\nH 0.04 0 2.31",
                  basis="sto-3g",
                  cas={"active_electrons": "2", "active_orbitals": "2",
                       "frozen_core": "1", "max_det": "5000", "max_memory": "512"})
    baseline = values(_run_caspt2(own, "dense_serial", **kwargs,
                                  pt2={"reference": "casci", "engine": "dense", "nproc": "1"}))
    records = {"dense": baseline.tolist()}
    comm.Barrier()
    for engine in ("direct", "fortran"):
        result = values(_run_caspt2(out, engine, **kwargs, usempi=True,
                                   pt2={"reference": "casci", "engine": engine, "nproc": "1"}))
        np.testing.assert_allclose(result, baseline, rtol=0, atol=1e-9)
        all_results = comm.allgather(result)
        for peer in all_results:
            np.testing.assert_allclose(peer, baseline, rtol=0, atol=1e-9)
        comm.Barrier()
        if comm.rank == 0:
            log = (out / f"{engine}.log").read_text()
            assert log.count("PyOQP PT2 engine:") == 1
            assert log.count("generated external singles/doubles") == 1
        records[engine] = result.tolist()
        comm.Barrier()
    (own / "results.json").write_text(json.dumps(records, indent=2))
    if comm.rank == 0:
        print("CASPT2_DIRECT_MPI_PASS: both ranks match serial dense; one log writer")


if __name__ == "__main__":
    main()
