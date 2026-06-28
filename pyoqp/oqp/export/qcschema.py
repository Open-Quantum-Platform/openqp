"""QCSchema AtomicResult export for an OQP (MRSF) calculation.

Produces a payload that validates with ``qcelemental.models.AtomicResult`` and
round-trips OQP's own energies to machine precision.  Excited-state data
(excitation energies, oscillator strengths, transition dipoles) are carried in
``extras`` since they are not first-class AtomicResultProperties fields.
"""
import numpy as np
from oqp.periodic_table import ELEMENTS_NAME

__all__ = ["to_qcschema", "validate_qcschema"]

_HARTREE2EV = 27.211386245988


def _symbols(mol):
    Z = np.asarray(mol.get_atoms(), dtype=int)
    return [ELEMENTS_NAME[int(z)].strip() for z in Z]


def to_qcschema(mol, states=None, driver="energy"):
    """Build a QCSchema AtomicResult payload (dict) from a finished run.

    ``states`` is an optional :class:`MRSFExcitedStates` to attach excited-state
    properties."""
    symbols = _symbols(mol)
    geometry = np.asarray(mol.get_system(), dtype=float).reshape(-1, 3)  # Bohr
    natom = len(symbols)
    nbf = int(round(np.asarray(mol.data["OQP::VEC_MO_A"]).size ** 0.5))
    na = int(np.asarray(mol.data["nelec_A"]).ravel()[0])
    nb = int(np.asarray(mol.data["nelec_B"]).ravel()[0])

    scf_energy = float(mol.get_scf_energy())
    cfg = mol.config
    method = cfg.get("input", {}).get("method", "tdhf")
    functional = cfg.get("input", {}).get("functional", "")
    td_type = cfg.get("tdhf", {}).get("type", "")
    method_name = f"{td_type or method}".upper()
    if functional:
        method_name = f"{method_name}/{functional.upper()}"

    try:
        enuc = float(mol.get_scf_energy("vnn"))
    except Exception:
        enuc = None

    molecule = {
        "symbols": symbols,
        "geometry": geometry.ravel().tolist(),
        "molecular_charge": float(cfg.get("input", {}).get("charge", 0) or 0),
        "molecular_multiplicity": int(cfg.get("scf", {}).get("multiplicity", 1) or 1),
        "fix_com": True,
        "fix_orientation": True,
    }

    properties = {
        "calcinfo_natom": natom,
        "calcinfo_nbasis": nbf,
        "calcinfo_nmo": nbf,
        "calcinfo_nalpha": na,
        "calcinfo_nbeta": nb,
        "scf_total_energy": scf_energy,
        "return_energy": scf_energy,
    }
    if enuc is not None:
        properties["nuclear_repulsion_energy"] = enuc

    extras = {"oqp": {"scf_total_energy": scf_energy}}
    if states is not None:
        de_ev = ((states.energies - states.energies[0]) * _HARTREE2EV).tolist()
        osc = []
        tdip = []
        for j in range(1, states.nstates):
            osc.append(states.oscillator_strength(0, j))
            tdip.append(states.transition_dipole(0, j).tolist())
        extras["oqp"].update({
            "state_energies_hartree": states.energies.tolist(),
            "excitation_energies_ev": de_ev,
            "oscillator_strengths": osc,            # S0 -> Sn, n=1..
            "transition_dipoles_au": tdip,
            "n_states": int(states.nstates),
        })

    payload = {
        "schema_name": "qcschema_output",
        "schema_version": 1,
        "molecule": molecule,
        "driver": driver,
        "model": {"method": method_name, "basis": cfg.get("input", {}).get("basis", "")},
        "keywords": {
            "scf_type": cfg.get("scf", {}).get("type", ""),
            "tdhf_type": td_type,
            "nstate": cfg.get("tdhf", {}).get("nstate", 0),
        },
        "properties": properties,
        "return_result": scf_energy,
        "success": True,
        "provenance": {"creator": "OpenQP", "routine": "tdhf_mrsf_energy"},
        "extras": extras,
    }
    return payload


def validate_qcschema(payload):
    """Validate the payload with qcelemental; returns the AtomicResult model."""
    import qcelemental as qcel
    return qcel.models.AtomicResult(**payload)
