"""Small solvent-coupling helpers used by staged PCM scaffolding."""

from __future__ import annotations

from collections.abc import Mapping
from math import isfinite, isqrt


def _as_float_list(values, *, name: str):
    if any(isinstance(value, bool) for value in values):
        raise ValueError(f"{name} values must be numeric, not boolean")
    floats = [float(value) for value in values]
    if any(not isfinite(value) for value in floats):
        raise ValueError(f"{name} values must be finite")
    return floats


def _packed_nbf(packed_length: int) -> int:
    discriminant = 1 + 8 * packed_length
    root = isqrt(discriminant)
    if root * root != discriminant or (root - 1) % 2 != 0:
        raise ValueError("reaction-field data must use triangular packed AO length")
    return (root - 1) // 2


def reference_scf_total_density(density_blocks):
    """Return the packed total AO density for reference-SCF PCM prototypes.

    OpenQP's first PCM target uses the RHF/ROHF reference density to build the
    ddX cavity potential.  RHF-style inputs have one packed density block, while
    ROHF/UHF-style storage has separate alpha/beta blocks that must be summed
    before evaluating the scalar electrostatic solvent response.
    """
    blocks = [_as_float_list(block, name="density_block") for block in density_blocks]
    if not blocks:
        raise ValueError("reference_scf_total_density requires at least one density block")
    if len(blocks) > 2:
        raise ValueError("reference_scf_total_density accepts one or two density blocks")
    packed_length = len(blocks[0])
    if packed_length == 0:
        raise ValueError("packed density blocks must not be empty")
    if any(len(block) != packed_length for block in blocks):
        raise ValueError("all density blocks must have the same packed length")
    return [sum(block[index] for block in blocks) for index in range(packed_length)]


def reference_scf_reaction_field_contract(density_blocks, reaction_potential):
    """Validate the first-scope packed AO reaction-field matrix contract.

    This is a dependency-light guard for the future `calc_jk_xc` hook: the
    reaction-field AO matrix must already be a packed triangular matrix matching
    the RHF/ROHF reference density length.  The helper deliberately does not
    enable runtime PCM or define final polarization-energy bookkeeping.
    """
    total_density = reference_scf_total_density(density_blocks)
    potential = _as_float_list(reaction_potential, name="reaction_potential")
    if len(potential) != len(total_density):
        raise ValueError("reaction_potential must have the same packed length as the reference density")
    nbf = _packed_nbf(len(potential))
    return {
        "nbf": nbf,
        "nfocks": len(density_blocks),
        "total_density": total_density,
        "reaction_potential": potential,
    }


def reference_scf_phi_cav_inputs(density_blocks, cavity_xyz):
    """Validate inputs for reference-density electrostatic potential on a ddX cavity.

    This prepares the dependency-light contract for the future
    ``electrostatic_potential_unweighted`` call that supplies ddX ``phi_cav``
    from the RHF/ROHF reference density.  It returns the summed reference
    density as ``density_packed`` so future callers cannot accidentally pass raw
    spin-density blocks into the scalar MEP path.  It does not enable runtime
    PCM or define solvent energy bookkeeping.
    """
    total_density = reference_scf_total_density(density_blocks)
    nbf = _packed_nbf(len(total_density))
    xyz_values = _as_float_list(cavity_xyz, name="cavity_xyz")
    if not xyz_values:
        raise ValueError("cavity_xyz must contain at least one cavity point")
    if len(xyz_values) % 3 != 0:
        raise ValueError("cavity_xyz must contain 3 * ncav values")
    return {
        "nbf": nbf,
        "ncav": len(xyz_values) // 3,
        "total_density": total_density,
        "density_packed": total_density,
        "cavity_xyz": xyz_values,
        "x": xyz_values[0::3],
        "y": xyz_values[1::3],
        "z": xyz_values[2::3],
    }


def reference_scf_pcm_coupling_contract(density_blocks, cavity_xyz, reaction_potential):
    """Package the first-scope reference-SCF PCM coupling handoff.

    The future SCF hook needs two validated dependency-light inputs from the
    same RHF/ROHF reference density: ``phi_cav`` inputs for ddX and a packed AO
    reaction-field matrix for ``calc_jk_xc``.  This helper keeps those contracts
    synchronized without enabling runtime PCM, state-specific response solvent,
    gradients, or final polarization-energy bookkeeping.
    """
    phi_cav = reference_scf_phi_cav_inputs(density_blocks, cavity_xyz)
    reaction = reference_scf_reaction_field_contract(density_blocks, reaction_potential)
    return {
        "nbf": phi_cav["nbf"],
        "nfocks": reaction["nfocks"],
        "ncav": phi_cav["ncav"],
        "total_density": phi_cav["total_density"],
        "density_packed": phi_cav["density_packed"],
        "cavity_xyz": phi_cav["cavity_xyz"],
        "x": phi_cav["x"],
        "y": phi_cav["y"],
        "z": phi_cav["z"],
        "reaction_potential": reaction["reaction_potential"],
        "pcm_scope": "reference_scf_energy_only",
        "reference_target": "RHF/ROHF reference density",
        "response_solvent_coupling": "not enabled",
        "gradient_support": "not enabled",
        "runtime_pcm_enabled": False,
    }


def reference_scf_pcm_energy_terms(density_blocks, reaction_potential):
    """Return guarded reference-SCF PCM energy bookkeeping terms.

    This dependency-light helper keeps the candidate polarization-energy
    convention tied to the same packed RHF/ROHF reference density and reaction
    potential validated for the future SCF hook.  It records only the
    conventional ``0.5 * D dot V`` bookkeeping candidate; backend-specific
    energy conventions still need independent validation before runtime PCM is
    enabled.
    """
    reaction = reference_scf_reaction_field_contract(density_blocks, reaction_potential)
    density = reaction["total_density"]
    potential = reaction["reaction_potential"]
    density_reaction_dot = sum(d_value * v_value for d_value, v_value in zip(density, potential))
    return {
        "nbf": reaction["nbf"],
        "nfocks": reaction["nfocks"],
        "density_packed": density,
        "reaction_potential": potential,
        "density_reaction_dot": density_reaction_dot,
        "candidate_polarization_energy": 0.5 * density_reaction_dot,
        "energy_convention": "0.5 * dot(reference_density_packed, reaction_potential)",
        "pcm_scope": "reference_scf_energy_only",
        "reference_target": "RHF/ROHF reference density",
        "response_solvent_coupling": "not enabled",
        "gradient_support": "not enabled",
        "runtime_pcm_enabled": False,
        "backend_validation_status": "pending PySCF/ddX/reference cross-check",
    }


def reference_scf_reaction_fock_updates(density_blocks, reaction_potential):
    """Return guarded packed AO reaction-potential updates for reference Focks.

    The first runtime PCM seam will add a validated packed AO reaction-field
    matrix to each RHF/ROHF reference-SCF Fock block.  This helper records the
    intended block replication and bookkeeping without enabling production PCM,
    state-specific response solvent, gradients, or optimizer support.
    """
    terms = reference_scf_pcm_energy_terms(density_blocks, reaction_potential)
    potential = terms["reaction_potential"]
    return {
        "nbf": terms["nbf"],
        "nfocks": terms["nfocks"],
        "density_packed": terms["density_packed"],
        "reaction_potential": potential,
        "fock_updates": [list(potential) for _ in range(terms["nfocks"])],
        "density_reaction_dot": terms["density_reaction_dot"],
        "candidate_polarization_energy": terms["candidate_polarization_energy"],
        "application_scope": "add reaction_potential to each reference SCF Fock block",
        "pcm_scope": "reference_scf_energy_only",
        "reference_target": "RHF/ROHF reference density",
        "response_solvent_coupling": "not enabled",
        "gradient_support": "not enabled",
        "runtime_pcm_enabled": False,
        "backend_validation_status": "pending PySCF/ddX/reference cross-check",
    }


def reference_scf_pcm_energy_handoff(density_blocks, cavity_xyz, reaction_potential):
    """Package the current no-runtime reference-SCF PCM energy handoff.

    The staged PCM prototype needs the ddX ``phi_cav`` input, the packed AO
    reaction-potential Fock update, and the candidate reference-SCF energy
    bookkeeping to remain tied to the same RHF/ROHF density.  This helper is a
    dependency-light handoff summary only; it does not enable runtime PCM,
    response-solvent coupling, gradients, or state-specific corrections.
    """
    coupling = reference_scf_pcm_coupling_contract(density_blocks, cavity_xyz, reaction_potential)
    updates = reference_scf_reaction_fock_updates(density_blocks, reaction_potential)
    energy_terms = reference_scf_pcm_energy_terms(density_blocks, reaction_potential)
    return {
        "nbf": coupling["nbf"],
        "nfocks": coupling["nfocks"],
        "ncav": coupling["ncav"],
        "phi_cav_inputs": {
            "nbf": coupling["nbf"],
            "ncav": coupling["ncav"],
            "density_packed": coupling["density_packed"],
            "cavity_xyz": coupling["cavity_xyz"],
            "x": coupling["x"],
            "y": coupling["y"],
            "z": coupling["z"],
        },
        "reaction_potential": coupling["reaction_potential"],
        "fock_updates": updates["fock_updates"],
        "energy_terms": energy_terms,
        "handoff_stage": "reference_scf_pcm_energy_prototype",
        "pcm_scope": "reference_scf_energy_only",
        "reference_target": "RHF/ROHF reference density",
        "response_solvent_coupling": "not enabled",
        "gradient_support": "not enabled",
        "runtime_pcm_enabled": False,
        "backend_validation_status": "pending PySCF/ddX/reference cross-check",
    }


def reference_scf_pcm_runtime_payload(density_blocks, reaction_potential):
    """Return a guarded no-runtime payload for future PCM runtime storage.

    The future SCF caller needs a validated packed ``pcm_reaction_potential``
    and the matching candidate ``epcm`` bookkeeping value to travel together.
    This helper deliberately records only a compact OpenQP-tag-style payload for
    review/prototype storage; it does not enable runtime PCM or expose raw
    spin-density blocks/state densities.
    """
    terms = reference_scf_pcm_energy_terms(density_blocks, reaction_potential)
    return {
        "OQP::pcm_reaction_potential": terms["reaction_potential"],
        "OQP::pcm_epcm": terms["candidate_polarization_energy"],
        "nbf": terms["nbf"],
        "packed_ao_length": len(terms["reaction_potential"]),
        "expected_packed_ao_length": terms["nbf"] * (terms["nbf"] + 1) // 2,
        "packed_ao_shape_formula": "nbf * (nbf + 1) / 2",
        "pcm_runtime_payload_version": 1,
        "pcm_scope": "reference_scf_energy_only",
        "reference_target": "RHF/ROHF reference density",
        "response_solvent_coupling": "not enabled",
        "gradient_support": "not enabled",
        "runtime_pcm_enabled": False,
        "backend_validation_status": "pending PySCF/ddX/reference cross-check",
    }


def reference_scf_pcm_reaction_potential_from_payload(payload):
    """Extract a reviewed packed reference-PCM reaction potential for calc_fock.

    This is the dependency-light consumer-side gate for a future prototype call
    to ``calc_fock(..., pcm_reaction_potential_in=...)``.  It accepts only the
    reviewed disabled-runtime/reference-SCF payload emitted by
    ``reference_scf_pcm_runtime_payload()``; it does not enable ``[pcm]`` runtime
    coupling and it rejects state-specific or malformed payloads before a packed
    AO reaction potential can reach the SCF Fock builder.
    """
    if not isinstance(payload, Mapping):
        raise ValueError("PCM runtime payload must be a mapping")
    if payload.get("pcm_scope") != "reference_scf_energy_only":
        raise ValueError("PCM runtime payload must have pcm_scope=reference_scf_energy_only")
    if payload.get("reference_target") != "RHF/ROHF reference density":
        raise ValueError("PCM runtime payload reference_target must be RHF/ROHF reference density")
    if payload.get("runtime_pcm_enabled") is not False:
        raise ValueError("PCM runtime payload runtime_pcm_enabled must be False")
    if payload.get("response_solvent_coupling") != "not enabled":
        raise ValueError("PCM runtime payload response_solvent_coupling must be not enabled")
    if payload.get("gradient_support") != "not enabled":
        raise ValueError("PCM runtime payload gradient_support must be not enabled")
    if payload.get("pcm_runtime_payload_version") != 1:
        raise ValueError("PCM runtime payload version must be 1")
    if payload.get("backend_validation_status") != "pending PySCF/ddX/reference cross-check":
        raise ValueError(
            "PCM runtime payload backend_validation_status must be pending PySCF/ddX/reference cross-check"
        )
    if "density_blocks" in payload:
        raise ValueError("PCM runtime payload must not include raw density_blocks")
    if "state_density" in payload:
        raise ValueError("PCM runtime payload must not include state_density")
    if "OQP::pcm_reaction_potential" not in payload:
        raise ValueError("PCM runtime payload missing OQP::pcm_reaction_potential")
    if "OQP::pcm_epcm" not in payload:
        raise ValueError("PCM runtime payload missing OQP::pcm_epcm")

    potential = _as_float_list(payload["OQP::pcm_reaction_potential"], name="OQP::pcm_reaction_potential")
    if not potential:
        raise ValueError("OQP::pcm_reaction_potential must not be empty")
    nbf = _packed_nbf(len(potential))
    expected_packed_ao_length = nbf * (nbf + 1) // 2
    for required_shape_key in (
        "nbf",
        "packed_ao_length",
        "expected_packed_ao_length",
        "packed_ao_shape_formula",
    ):
        if required_shape_key not in payload:
            raise ValueError(f"PCM runtime payload missing {required_shape_key}")
    payload_nbf = payload["nbf"]
    if isinstance(payload_nbf, bool) or not isinstance(payload_nbf, int):
        raise ValueError("PCM runtime payload nbf must be an integer")
    if payload_nbf != nbf:
        raise ValueError("PCM runtime payload nbf must match OQP::pcm_reaction_potential packed length")
    if "packed_ao_length" in payload:
        payload_packed_ao_length = payload["packed_ao_length"]
        if isinstance(payload_packed_ao_length, bool) or not isinstance(payload_packed_ao_length, int):
            raise ValueError("PCM runtime payload packed_ao_length must be an integer")
        if payload_packed_ao_length != len(potential):
            raise ValueError("PCM runtime payload packed_ao_length must match OQP::pcm_reaction_potential length")
    if "expected_packed_ao_length" in payload:
        payload_expected_length = payload["expected_packed_ao_length"]
        if isinstance(payload_expected_length, bool) or not isinstance(payload_expected_length, int):
            raise ValueError("PCM runtime payload expected_packed_ao_length must be an integer")
        if payload_expected_length != expected_packed_ao_length:
            raise ValueError("PCM runtime payload expected_packed_ao_length must match nbf * (nbf + 1) / 2")
    if payload.get("packed_ao_shape_formula", "nbf * (nbf + 1) / 2") != "nbf * (nbf + 1) / 2":
        raise ValueError("PCM runtime payload packed_ao_shape_formula must be nbf * (nbf + 1) / 2")
    epcm_value = payload["OQP::pcm_epcm"]
    if isinstance(epcm_value, bool):
        raise ValueError("OQP::pcm_epcm must be numeric, not boolean")
    epcm = float(epcm_value)
    if not isfinite(epcm):
        raise ValueError("OQP::pcm_epcm must be finite")

    return {
        "reaction_potential": potential,
        "nbf": nbf,
        "packed_ao_length": len(potential),
        "expected_packed_ao_length": expected_packed_ao_length,
        "packed_ao_shape_formula": "nbf * (nbf + 1) / 2",
        "candidate_polarization_energy": epcm,
        "pcm_runtime_payload_version": 1,
        "pcm_scope": "reference_scf_energy_only",
        "reference_target": "RHF/ROHF reference density",
        "response_solvent_coupling": "not enabled",
        "gradient_support": "not enabled",
        "runtime_pcm_enabled": False,
        "backend_validation_status": payload.get(
            "backend_validation_status",
            "pending PySCF/ddX/reference cross-check",
        ),
        "handoff_target": "calc_fock pcm_reaction_potential_in",
    }


def reference_scf_pcm_calc_fock_handoff(payload):
    """Package a reviewed PCM payload as opt-in ``calc_fock`` kwargs.

    This helper is the final Python-side guard before a future prototype caller
    forwards a packed reference-SCF reaction potential into
    ``calc_fock(..., pcm_reaction_potential_in=...)``.  It intentionally exposes
    only the keyword argument and compact disabled-runtime metadata; raw density
    blocks and state/response densities must not pass this boundary.
    """
    reviewed = reference_scf_pcm_reaction_potential_from_payload(payload)
    return {
        "calc_fock_kwargs": {
            "pcm_reaction_potential_in": list(reviewed["reaction_potential"]),
        },
        "nbf": reviewed["nbf"],
        "packed_ao_length": reviewed["packed_ao_length"],
        "expected_packed_ao_length": reviewed["expected_packed_ao_length"],
        "packed_ao_shape_formula": reviewed["packed_ao_shape_formula"],
        "candidate_polarization_energy": reviewed["candidate_polarization_energy"],
        "pcm_runtime_payload_version": reviewed["pcm_runtime_payload_version"],
        "pcm_scope": reviewed["pcm_scope"],
        "reference_target": reviewed["reference_target"],
        "response_solvent_coupling": reviewed["response_solvent_coupling"],
        "gradient_support": reviewed["gradient_support"],
        "runtime_pcm_enabled": reviewed["runtime_pcm_enabled"],
        "backend_validation_status": reviewed["backend_validation_status"],
        "handoff_target": reviewed["handoff_target"],
    }


def reference_scf_pcm_calc_fock_handoff_from_molecule(mol):
    """Return opt-in ``calc_fock`` kwargs only when a reviewed payload exists.

    This is a no-runtime prototype call-site gate: a molecule with no restored
    PCM runtime payload produces no ``calc_fock`` keyword arguments, while a
    molecule carrying the reviewed reference-SCF payload must pass the same
    consumer validation as ``reference_scf_pcm_calc_fock_handoff()`` before the
    packed reaction potential can be forwarded.
    """
    payload = mol.get_pcm_runtime_payload()
    if not isinstance(payload, Mapping):
        raise ValueError("PCM runtime payload must be a mapping")
    if not payload:
        return {
            "calc_fock_kwargs": {},
            "payload_present": False,
            "pcm_scope": "reference_scf_energy_only",
            "reference_target": "RHF/ROHF reference density",
            "runtime_pcm_enabled": False,
            "response_solvent_coupling": "not enabled",
            "gradient_support": "not enabled",
            "backend_validation_status": "pending PySCF/ddX/reference cross-check",
            "handoff_target": "calc_fock pcm_reaction_potential_in",
        }
    handoff = reference_scf_pcm_calc_fock_handoff(payload)
    handoff["payload_present"] = True
    return handoff


def reference_scf_pcm_calc_fock_request(mol, *, incremental_fock: bool = False):
    """Return a guarded prototype request for ``calc_fock`` PCM handoff.

    The reviewed reference-SCF PCM payload may only be forwarded through the
    non-incremental Fock path.  If no payload is present this returns an explicit
    disabled/no-op request; if a payload is present and the caller is attempting
    the incremental-Fock shortcut, fail fast until PCM incremental-energy
    behavior is derived and validated.
    """
    if not isinstance(incremental_fock, bool):
        raise ValueError("incremental_fock must be boolean")
    handoff = reference_scf_pcm_calc_fock_handoff_from_molecule(mol)
    if not handoff["calc_fock_kwargs"]:
        handoff["call_mode"] = "disabled_no_payload"
        handoff["requires_non_incremental_fock"] = False
        handoff["incremental_fock_allowed"] = False
        return handoff
    if incremental_fock:
        raise ValueError("reference PCM incremental Fock is not validated")
    handoff["call_mode"] = "non_incremental_only"
    handoff["requires_non_incremental_fock"] = True
    handoff["incremental_fock_allowed"] = False
    return handoff


def reference_scf_pcm_incremental_fock_audit(*, dens_old=None, f_old=None):
    """Return explicit old-buffer metadata for the reference-PCM guard.

    OpenQP's incremental Fock shortcut can be triggered by either an old density
    or an old Fock buffer.  The first reference-SCF PCM path must preserve these
    as separate audit fields, not just a combined boolean, so blocked prototype
    call sites remain diagnosable.
    """
    dens_old_present = dens_old is not None
    f_old_present = f_old is not None
    triggers = []
    if dens_old_present:
        triggers.append("dens_old")
    if f_old_present:
        triggers.append("f_old")
    return {
        "scf_state_incremental_fock": bool(triggers),
        "scf_state_dens_old_present": dens_old_present,
        "scf_state_f_old_present": f_old_present,
        "incremental_trigger_fields": triggers,
        "scf_state_guard": "dens_old/f_old presence blocks reviewed reference PCM payloads",
    }


def reference_scf_pcm_calc_fock_request_from_scf_state(mol, *, dens_old=None, f_old=None):
    """Derive the PCM ``calc_fock`` request from current SCF buffer state.

    A future SCF caller can use this dependency-light shim to map OpenQP's
    old-density/old-Fock buffers onto the existing non-incremental-only request
    guard.  Presence of either buffer means the incremental-Fock shortcut is in
    play and reviewed reference-PCM payloads must fail fast with explicit old
    buffer provenance.
    """
    audit = reference_scf_pcm_incremental_fock_audit(dens_old=dens_old, f_old=f_old)
    try:
        request = reference_scf_pcm_calc_fock_request(mol, incremental_fock=audit["scf_state_incremental_fock"])
    except ValueError as exc:
        if audit["scf_state_incremental_fock"] and "reference PCM incremental Fock is not validated" in str(exc):
            raise ValueError(
                "reference PCM incremental Fock is not validated; "
                f"dens_old_present={str(audit['scf_state_dens_old_present']).lower()}; "
                f"f_old_present={str(audit['scf_state_f_old_present']).lower()}; "
                f"incremental_trigger_fields={','.join(audit['incremental_trigger_fields'])}"
            ) from exc
        raise
    request.update(audit)
    return request


def reference_scf_pcm_calc_fock_call_site_bridge(mol, *, dens_old=None, f_old=None):
    """Return the guarded SCF call-site payload for ``calc_fock``.

    This dependency-light bridge is the final staging helper before a future SCF
    caller forwards ``pcm_reaction_potential_in``.  It forwards a reviewed
    reference-PCM payload only when the request is non-incremental; a no-payload
    molecule remains an explicit disabled/no-op bridge even if old SCF buffers
    are present.
    """
    request = reference_scf_pcm_calc_fock_request_from_scf_state(mol, dens_old=dens_old, f_old=f_old)
    request["call_site_bridge"] = "reference_scf_calc_fock"
    request["forward_pcm_reaction_potential"] = bool(request["calc_fock_kwargs"])
    return request


def _require_provisional_opt_in(allow_provisional: bool):
    if allow_provisional is not True:
        raise ValueError(
            "allow_provisional must be the boolean True; ddX q_cav sign/scale is provisional "
            "and may only be used in guarded validation/prototype code"
        )


def provisional_ddx_external_charges(q_cav, *, allow_provisional: bool = False):
    """Return candidate OpenQP external-charge weights from ddX ``q_cav``.

    The current ddPCM finite-difference smoke suggests ``chg = -0.5*q_cav``
    for the future ``external_charge_potential`` AO-matrix seam.  The
    sign/scale is still provisional, so callers must opt in explicitly until
    it is cross-checked against PySCF/ddX/reference-package data.
    """
    _require_provisional_opt_in(allow_provisional)
    return [-0.5 * value for value in _as_float_list(q_cav, name="q_cav")]


def provisional_ddx_reaction_field_inputs(q_cav, cavity_xyz, *, allow_provisional: bool = False):
    """Validate candidate ddX cavity data for the future AO reaction-field seam.

    ``cavity_xyz`` is the flat ``(3, ncav)`` coordinate buffer copied from ddX,
    and ``q_cav`` is the projected cavity quantity whose sign/scale remains
    provisional.  This helper prepares only guarded prototype inputs; it does
    not enable runtime PCM coupling.
    """
    q_values = _as_float_list(q_cav, name="q_cav")
    xyz_values = _as_float_list(cavity_xyz, name="cavity_xyz")
    expected_xyz = 3 * len(q_values)
    if not q_values:
        raise ValueError("q_cav must contain at least one cavity value")
    if len(xyz_values) != expected_xyz:
        raise ValueError(
            "cavity_xyz must contain 3 * len(q_cav) values "
            f"({expected_xyz} expected, got {len(xyz_values)})"
        )
    charges = provisional_ddx_external_charges(q_values, allow_provisional=allow_provisional)
    return {
        "ncav": len(q_values),
        "cavity_xyz": xyz_values,
        "charges": charges,
        "x": xyz_values[0::3],
        "y": xyz_values[1::3],
        "z": xyz_values[2::3],
        "chg": charges,
        "provisional_sign_scale": True,
        "sign_scale_convention": "chg = -0.5 * q_cav",
        "validation_status": "requires PySCF/ddX/reference cross-check before runtime use",
        "pcm_scope": "reference_scf_energy_only",
        "reference_target": "RHF/ROHF reference density",
        "response_solvent_coupling": "not enabled",
        "gradient_support": "not enabled",
        "runtime_pcm_enabled": False,
        "backend_validation_status": "pending PySCF/ddX/reference cross-check",
    }
