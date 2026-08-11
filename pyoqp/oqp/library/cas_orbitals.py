"""Orbital-file and active-space helpers for CASCI/CASSCF drivers."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np


def _matrix_from_json_value(value, nbf: int, *, restart_orientation: bool = False) -> np.ndarray:
    """Return AO-rows x MO-columns coefficients.

    ``restart_orientation`` selects the convention of the OQP::VEC_MO_A /
    OQP::VEC_MO_B keys, which OpenQP writes MO-major; the ``mo_coeff*`` aliases
    are AO-rows x MO-columns as documented.

    The orientation used to be decided by ndim rather than by key, so only the
    FLAT form was transposed -- and ``Molecule.save_data`` writes OQP::VEC_MO_A
    as a 2-D matrix.  Feeding a saved restart file straight back through
    ``[cas] orbital_file`` therefore loaded the TRANSPOSE of the intended
    orbitals, silently: on H4/STO-3G it turned a -2.1147 singlet CASCI into a
    -3.7105 "triplet" below the RHF energy.  The convention now follows the key
    name, so both spellings of the same key round-trip.
    """
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 1:
        if arr.size != nbf * nbf:
            raise ValueError(f"MO coefficient vector has {arr.size} values, expected {nbf * nbf}")
        arr = arr.reshape((nbf, nbf))
        return arr.T if restart_orientation else arr
    if arr.ndim == 2:
        if arr.shape != (nbf, nbf):
            raise ValueError(f"MO coefficient matrix has shape {arr.shape}, expected {(nbf, nbf)}")
        return arr.T if restart_orientation else arr
    raise ValueError("MO coefficients must be a flat nbf*nbf array or an nbf x nbf matrix")


def load_mo_coeff_from_json(filename: str | Path, nbf: int, *, spin: str = "alpha") -> np.ndarray:
    """Load MO coefficients from an OpenQP-style JSON file.

    Accepted keys are ``OQP::VEC_MO_A``/``OQP::VEC_MO_B`` and convenience aliases
    ``mo_coeff``, ``mo_coeff_alpha``, and ``mo_coeff_beta``.  Matrices are
    expected to have AO rows and MO columns; flat OpenQP vectors use the existing
    OpenQP restart orientation.
    """

    with open(filename) as handle:
        data = json.load(handle)

    spin = spin.lower()
    keys = ["OQP::VEC_MO_A", "mo_coeff_alpha", "mo_coeff"] if spin.startswith("a") else ["OQP::VEC_MO_B", "mo_coeff_beta", "mo_coeff"]
    for key in keys:
        if key in data:
            return np.ascontiguousarray(
                _matrix_from_json_value(data[key], nbf,
                                        restart_orientation=key.startswith("OQP::")),
                dtype=np.float64)
    raise ValueError(f"No MO coefficient key found in {filename}; tried {', '.join(keys)}")


def check_mo_orthonormality(coeff: np.ndarray, overlap: np.ndarray, *,
                            label: str, tol: float = 1.0e-6) -> None:
    """Verify C^T S C = I for imported orbitals.

    Every consumer -- the integral transform, the determinant Hamiltonians, the
    CASSCF gradient -- assumes the coefficients are orthonormal in the CURRENT
    AO metric.  Shape agreement does not imply that: a restart from a different
    geometry or basis with the same nbf, a malformed matrix, or a transposed
    one all pass a shape check and then produce nonvariational nonsense with no
    diagnostic.  This is not hypothetical -- a transposed orbital file turned an
    H4/STO-3G CASCI into -3.7105, below that system's RHF energy, and nothing
    complained.
    """
    c = np.asarray(coeff, dtype=float)
    s_ao = np.asarray(overlap, dtype=float)
    if s_ao.shape != c.shape:
        return
    gram = c.T @ s_ao @ c
    dev = float(np.abs(gram - np.eye(gram.shape[0])).max())
    if not np.isfinite(dev) or dev > tol:
        raise ValueError(
            "orbitals loaded from %s are not orthonormal in the current AO "
            "basis: max |C^T S C - I| = %.3e (tolerance %.1e).  They are "
            "probably from a different geometry or basis, transposed, or "
            "otherwise malformed -- every CAS consumer assumes orthonormality, "
            "so continuing would give nonvariational energies with no other "
            "symptom." % (label, dev, tol))


def load_cas_mo_coeff(config: dict, nbf: int, default_coeff: np.ndarray,
                      overlap: np.ndarray | None = None,
                      input_dir: str | None = None) -> tuple[np.ndarray, str]:
    """Return the MO coefficient matrix requested by ``[cas] orbital_source``.

    ``overlap`` is the current AO overlap; when supplied, imported coefficients
    are checked against it before any consumer sees them.
    """

    cas = config.get("cas", {})
    source = str(cas.get("orbital_source", "rhf")).lower()
    orbital_file = str(cas.get("orbital_file", "") or "")
    if source in {"rhf", "canonical"}:
        return np.ascontiguousarray(default_coeff, dtype=np.float64), "rhf"
    if not orbital_file:
        raise ValueError(f"[cas] orbital_source={source} requires orbital_file")
    # A relative path in a legacy .inp means "beside the input", which is where
    # the user put it -- resolving against the process CWD instead made
    # `openqp examples/WF_methods/X.inp` from the repository root fail on a file
    # sitting next to X.inp.  The concise .oqp lowering already resolves this.
    if input_dir and not os.path.isabs(orbital_file):
        # Unconditionally, not "unless the CWD happens to have a file of the
        # same name": my first version made which file gets loaded depend on an
        # unrelated collision in the working directory.  A relative path in a
        # legacy deck means "beside the deck".
        orbital_file = os.path.join(input_dir, orbital_file)
    if source == "json":
        loaded = load_mo_coeff_from_json(orbital_file, nbf)
        if overlap is not None:
            check_mo_orthonormality(loaded, overlap, label=orbital_file)
        return loaded, f"json:{orbital_file}"
    raise ValueError(f"Unknown [cas] orbital_source={source}")


def _one_based_indices(raw, *, field: str) -> list[int]:
    indices = [int(x) for x in (raw or [])]
    if any(idx < 1 for idx in indices):
        raise ValueError(f"{field} uses 1-based orbital indices; all entries must be >= 1")
    if len(set(indices)) != len(indices):
        raise ValueError(f"{field} contains duplicate orbital indices")
    return indices


def active_core_lists(config: dict, norb: int) -> tuple[list[int], list[int]]:
    """Return zero-based explicit active/core lists from ``[cas]``.

    User input is 1-based because that matches printed
    orbital labels.  Empty lists mean sequential active-space inference.
    """

    cas = config.get("cas", {})
    active_1 = _one_based_indices(cas.get("active_orbital_indices", []), field="cas.active_orbital_indices")
    core_1 = _one_based_indices(cas.get("core_orbital_indices", []), field="cas.core_orbital_indices")
    for idx in active_1 + core_1:
        if idx > norb:
            raise ValueError(f"CAS orbital index {idx} exceeds available orbital count {norb}")
    if set(active_1) & set(core_1):
        raise ValueError("CAS active and core orbital index lists must not overlap")
    return [idx - 1 for idx in active_1], [idx - 1 for idx in core_1]
