"""Parsers that normalize external-code excited-state results for comparison.

* :func:`parse_output` wraps cclib to read Gaussian/ORCA/Q-Chem/... output files.
* :func:`parse_pyscf_tddft` reads a PySCF TDDFT object's native attributes.
* :func:`parse_oqp` reads an OQP run (in-process molecule or results dict).

All return a normalized dict::

    {'scf_energy_ha', 'excitation_energies_ev', 'oscillator_strengths',
     'mo_energies_ev'}
"""
import numpy as np

__all__ = ["parse_output", "parse_pyscf_tddft", "parse_oqp"]

_HARTREE2EV = 27.211386245988
_CM2EV = 1.0 / 8065.543937


_PROGRAMS = {
    "gaussian": "Gaussian", "orca": "ORCA", "qchem": "QChem",
    "psi4": "Psi4", "molpro": "Molpro", "nwchem": "NWChem", "gamess": "GAMESS",
}


def parse_output(path, program=None):
    """Parse a quantum-chemistry output file with cclib (Gaussian/ORCA/Q-Chem...).

    ``program`` (e.g. "gaussian") forces a specific cclib parser, bypassing
    cclib's file-type auto-detection (useful for partial/synthetic logs)."""
    import logging
    if program is not None:
        import cclib.parser as ccp
        cls = getattr(ccp, _PROGRAMS.get(program.lower(), program))
        data = cls(path, loglevel=logging.ERROR).parse()
    else:
        from cclib.io import ccread
        data = ccread(path)
    res = {"source": "cclib", "parser": type(data).__name__ if data else None}
    scf = getattr(data, "scfenergies", None)
    if scf is not None and len(scf):
        res["scf_energy_ha"] = float(scf[-1]) / _HARTREE2EV   # cclib stores eV
    ete = getattr(data, "etenergies", None)
    if ete is not None and len(ete):
        res["excitation_energies_ev"] = (np.asarray(ete, float) * _CM2EV).tolist()
    osc = getattr(data, "etoscs", None)
    if osc is not None and len(osc):
        res["oscillator_strengths"] = np.asarray(osc, float).tolist()
    mo = getattr(data, "moenergies", None)
    if mo is not None and len(mo):
        res["mo_energies_ev"] = np.asarray(mo[0], float).tolist()
    return res


def parse_pyscf_tddft(td, mf=None):
    """Read excitation energies / oscillator strengths from a PySCF TDDFT object."""
    res = {"source": "pyscf"}
    e = np.asarray(td.e, float)                       # excitation energies (Hartree)
    res["excitation_energies_ev"] = (e * _HARTREE2EV).tolist()
    try:
        res["oscillator_strengths"] = np.asarray(td.oscillator_strength(), float).tolist()
    except Exception:
        pass
    if mf is not None:
        res["scf_energy_ha"] = float(mf.e_tot)
        res["mo_energies_ev"] = (np.asarray(mf.mo_energy, float) * _HARTREE2EV).tolist()
    return res


def parse_oqp(mol, states=None):
    """Normalize an in-process OQP molecule (and optional MRSFExcitedStates)."""
    res = {"source": "oqp", "scf_energy_ha": float(mol.get_scf_energy())}
    if states is not None:
        de = (states.energies - states.energies[0]) * _HARTREE2EV
        res["excitation_energies_ev"] = de[1:].tolist()
        res["oscillator_strengths"] = [states.oscillator_strength(0, j)
                                       for j in range(1, states.nstates)]
    return res
