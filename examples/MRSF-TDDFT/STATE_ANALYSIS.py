#!/usr/bin/env python3
"""Physical-root MRSF state analysis for a finished OpenQP calculation.

Each MRSF root pair ``S_ref -> S_n`` is reported as two independent things:

* **LE/CT** -- how much of the transition moves density between the fragments
  you supply, from the fragment-partitioned state-interaction 1-TDM.
* **orbital character** -- the ``n/pi/sigma -> pi*/sigma*`` fractions of the
  same root pair's NTOs, after a Loewdin population analysis.

The analysis needs the state-interaction densities, and those only live in the
tagarray of a *finished MRSF run in the current process* -- so the run and the
analysis have to happen in one script, as below.

Usage
-----
    python STATE_ANALYSIS.py H2O_BHHLYP-MRSFTDDFT_ENERGY.inp
    python STATE_ANALYSIS.py my.inp 0,1 2      # fragment "atoms 0,1" | "atom 2"

Atom and state indices are zero-based.
"""

import os
import sys

from oqp.interop import AOBasis, MRSFExcitedStates, analyze_mrsf_transition


def report_states(mol, fragments, plane_normal=None):
    """Print and return the analysis of every ``S0 -> Sn`` pair of ``mol``."""
    states = MRSFExcitedStates(mol)
    ao = AOBasis(mol)
    reports = []
    for final_state in range(1, states.nstates):
        report = analyze_mrsf_transition(
            states, ao, final_state, fragments, ref=0,
            plane_normal=plane_normal)
        character = report["orbital_character"]
        print(
            f"S0 -> S{final_state}: {character['label']}; "
            f"LE={report['le_fraction']:.3f}, CT={report['ct_fraction']:.3f}; "
            f"fractions={character['fractions']}")
        if not character["classified"]:
            print(f"    orbital character unavailable: {character['reason']}")
        reports.append(report)
    return reports


def main():
    if len(sys.argv) < 2:
        sys.exit(f"usage: {sys.argv[0]} <input.inp> [fragment ...]")
    input_file = sys.argv[1]

    from oqp.pyoqp import Runner

    project = os.path.splitext(input_file)[0]
    runner = Runner(project=project, input_file=input_file,
                    log=f"{project}.log", usempi=False)
    runner.run()

    # One fragment per command-line argument, e.g. "0,1" "2".  With no
    # arguments the whole molecule is a single fragment, so LE=1 by
    # construction -- pass donor/acceptor fragments to see any CT.
    if len(sys.argv) > 2:
        fragments = [[int(a) for a in group.split(",")] for group in sys.argv[2:]]
    else:
        fragments = [list(range(len(runner.mol.get_atoms())))]
    print(f"fragments = {fragments}")

    report_states(runner.mol, fragments)


if __name__ == "__main__":
    main()
