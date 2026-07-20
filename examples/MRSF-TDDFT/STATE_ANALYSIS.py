"""Physical-root MRSF state analysis for a finished OpenQP calculation.

Run the MRSF energy job in the current process, then pass its ``mol`` object to
``report_states``.  Atom and state indices below are zero-based.
"""
from oqp.interop import AOBasis, MRSFExcitedStates, analyze_mrsf_transition


def report_states(mol, fragments, plane_normal=None):
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
        reports.append(report)
    return reports


# Example after running formaldehyde: one chromophore fragment means LE=1 by
# construction.  Use chemically meaningful donor/acceptor fragments for CT.
# reports = report_states(runner.mol, fragments=[[0, 1, 2, 3]])
