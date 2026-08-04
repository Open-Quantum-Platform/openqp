#!/usr/bin/env python3
"""Build a one-dimensional PMF from packed OpenQP ODP window trajectories.

Run one short NAMD calculation per umbrella window, then pass every resulting
``*.namd.trj`` file together. The temperature is explicit because an NVE
trajectory's initial velocity temperature is not automatically an NVT sampling
temperature.

Example
-------
python odp_wham.py --temperature 300 --output odp-wham.npz window-*.namd.trj
"""

from __future__ import annotations

import argparse

import numpy as np

from oqp.library import odp_wham


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectories", nargs="+", help="packed ODP *.namd.trj files")
    parser.add_argument("--temperature", type=float, required=True, help="WHAM temperature in K")
    parser.add_argument("--bins", type=int, default=100, help="number of xi bins")
    parser.add_argument("--discard", type=int, default=0, help="records to discard per window")
    parser.add_argument("--stride", type=int, default=1, help="record stride per window")
    parser.add_argument("--output", default="odp-wham.npz", help="portable output archive")
    args = parser.parse_args(argv)

    result = odp_wham(
        args.trajectories,
        args.temperature,
        bins=args.bins,
        discard=args.discard,
        stride=args.stride,
        output=args.output,
    )
    occupied = np.isfinite(result["free_energy_hartree"])
    print(
        f"WHAM converged in {result['iterations']} iteration(s); "
        f"{int(np.count_nonzero(occupied))}/{occupied.size} bins occupied; "
        f"effective sample size={result['effective_sample_size']:.1f}; "
        f"wrote {args.output}"
    )
    return result


if __name__ == "__main__":
    main()
