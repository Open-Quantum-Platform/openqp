"""Post-run MRSF physical-root NTO and density analysis.

OpenQP writes state-interaction density matrices in its full JSON export.  The
matrices are analysed here after the engine process exits, using the matching
Molden file for AO ordering and cube-grid evaluation.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import molden

HARTREE_EV = 27.211386245988
KINDS = {
    "nto_hole", "nto_particle", "attachment", "detachment",
    "difference", "transition", "state_density",
}


class AnalysisUnavailable(ValueError):
    pass


def _json_payload(path: Path) -> dict | None:
    try:
        payload = json.loads(path.read_text())
    except (OSError, ValueError):
        return None
    if not isinstance(payload, dict):
        return None
    data = payload.get("data")
    if isinstance(data, dict) and "OQP::td_trans_density_mo" in data:
        return data
    return payload


def _pick_molden(paths: list[Path]) -> tuple[Path, molden.MoldenData]:
    candidates = [path for path in paths if path.suffix.lower() == ".molden"]
    candidates.sort(key=lambda path: (
        "dyson" in path.name.lower(),
        "freq" in path.name.lower(),
        "_scf_" not in path.name.lower(),
        path.name,
    ))
    for path in candidates:
        data = molden.parse_molden(path.read_text(errors="replace"))
        if data.supported and data.basis and data.orbitals:
            return path, data
    if candidates:
        raise AnalysisUnavailable("excited-state maps require a Cartesian-basis Molden file")
    raise AnalysisUnavailable("no SCF Molden file was found for the excited-state maps")


def _spin_orbitals(data: molden.MoldenData, nbf: int) -> tuple[np.ndarray, np.ndarray]:
    alpha = [orbital for orbital in data.orbitals if orbital.spin.lower().startswith("a")]
    beta = [orbital for orbital in data.orbitals if orbital.spin.lower().startswith("b")]
    if len(alpha) < nbf:
        alpha = data.orbitals[:nbf]
        beta = []
    if len(alpha) < nbf or any(orbital.coefficients.size != nbf for orbital in alpha[:nbf]):
        raise AnalysisUnavailable("Molden orbital dimensions do not match the saved MRSF density")
    coefficients = np.column_stack([orbital.coefficients for orbital in alpha[:nbf]])
    occupations = np.array([
        float(orbital.occupancy or 0.0) for orbital in alpha[:nbf]
    ])
    if len(beta) >= nbf:
        occupations += np.array([
            float(orbital.occupancy or 0.0) for orbital in beta[:nbf]
        ])
    return coefficients, occupations


@dataclass
class ExcitedStateData:
    energies: np.ndarray
    transition_densities: np.ndarray
    coefficients: np.ndarray
    reference_occupations: np.ndarray
    molden_data: molden.MoldenData
    json_name: str
    molden_name: str

    @property
    def nstates(self) -> int:
        return int(self.energies.size)

    @property
    def nbf(self) -> int:
        return int(self.coefficients.shape[0])

    @classmethod
    def load(cls, paths: list[Path]) -> ExcitedStateData:
        json_path = None
        payload = None
        for path in sorted(paths, key=lambda item: ("hess" in item.name.lower(), item.name)):
            if path.suffix.lower() != ".json":
                continue
            candidate = _json_payload(path)
            if candidate and "OQP::td_trans_density_mo" in candidate:
                json_path, payload = path, candidate
                break
        if payload is None or json_path is None:
            raise AnalysisUnavailable(
                "this calculation has no saved MRSF state-interaction density; "
                "run it with guess(save_mol=true) using an engine that exports that tag"
            )

        energies = np.asarray(
            payload.get("OQP::td_energies", payload.get("td_energies", [])), dtype=float
        ).ravel()
        if energies.size < 2:
            raise AnalysisUnavailable("the saved MRSF result has fewer than two physical states")
        vec = np.asarray(payload.get("OQP::VEC_MO_A", []), dtype=float)
        nbf = math.isqrt(vec.size)
        if nbf <= 0 or nbf * nbf != vec.size:
            raise AnalysisUnavailable("the saved alpha-MO coefficient matrix is missing or invalid")
        raw = np.asarray(payload["OQP::td_trans_density_mo"], dtype=float).ravel(order="C")
        expected = nbf * nbf * energies.size * energies.size
        if raw.size != expected:
            raise AnalysisUnavailable(
                f"saved MRSF density has {raw.size} values; expected {expected}"
            )
        densities = raw.reshape((nbf, nbf, energies.size, energies.size), order="F")
        molden_path, molden_data = _pick_molden(paths)
        if len(molden_data.basis) != nbf:
            raise AnalysisUnavailable(
                "Molden basis size does not match the saved MRSF density matrix"
            )
        coefficients, occupations = _spin_orbitals(molden_data, nbf)
        return cls(
            energies, densities, coefficients, occupations, molden_data,
            json_path.name, molden_path.name,
        )

    def _check_pair(self, ref: int, target: int) -> None:
        if not (0 <= ref < self.nstates and 0 <= target < self.nstates):
            raise ValueError("state index out of range")
        if ref == target:
            raise ValueError("source and target states must differ")

    def tdm_mo(self, ref: int, target: int) -> np.ndarray:
        self._check_pair(ref, target)
        if ref <= target:
            return self.transition_densities[:, :, ref, target].copy()
        return self.transition_densities[:, :, target, ref].T.copy()

    def difference_mo(self, ref: int, target: int) -> np.ndarray:
        self._check_pair(ref, target)
        delta = (self.transition_densities[:, :, target, target]
                 - self.transition_densities[:, :, ref, ref])
        return 0.5 * (delta + delta.T)

    def nto(self, ref: int, target: int) -> dict:
        particles, singular, holes_t = np.linalg.svd(self.tdm_mo(ref, target))
        weights = singular * singular
        total = float(weights.sum())
        fractions = weights / total if total > 0 else np.zeros_like(weights)
        return {
            "singular": singular,
            "weights": weights,
            "fractions": fractions,
            "holes_ao": self.coefficients @ holes_t.T,
            "particles_ao": self.coefficients @ particles,
            "sum_weights": total,
            "participation_ratio": (
                float(total * total / np.sum(weights * weights))
                if np.any(weights) else 0.0
            ),
        }

    def attachment_detachment(self, ref: int, target: int) -> dict:
        eigenvalues, vectors = np.linalg.eigh(self.difference_mo(ref, target))
        positive = eigenvalues > 0
        negative = eigenvalues < 0
        attach = ((vectors[:, positive] * eigenvalues[positive])
                  @ vectors[:, positive].T)
        detach = ((vectors[:, negative] * (-eigenvalues[negative]))
                  @ vectors[:, negative].T)
        n_attach = float(eigenvalues[positive].sum())
        n_detach = float(-eigenvalues[negative].sum())
        return {
            "attachment_mo": attach,
            "detachment_mo": detach,
            "n_attach": n_attach,
            "n_detach": n_detach,
            "n_promoted": 0.5 * (n_attach + n_detach),
        }

    def ao_density(self, matrix_mo: np.ndarray) -> np.ndarray:
        return self.coefficients @ matrix_mo @ self.coefficients.T

    def summary(self, ref: int = 0, target: int = 1) -> dict:
        self._check_pair(ref, target)
        nto = self.nto(ref, target)
        ad = self.attachment_detachment(ref, target)
        largest_weight = float(np.max(nto["weights"])) if nto["weights"].size else 0.0
        cutoff = max(largest_weight * 1.0e-6, 1.0e-14)
        pairs = [
            {
                "index": index,
                "singular_value": float(sigma),
                "weight": float(weight),
                "fraction": float(fraction),
            }
            for index, (sigma, weight, fraction) in enumerate(zip(
                nto["singular"], nto["weights"], nto["fractions"]
            ))
            if weight >= cutoff
        ]
        states = [
            {
                "index": index,
                "label": f"S{index}",
                "relative_ev": float((energy - self.energies[0]) * HARTREE_EV),
            }
            for index, energy in enumerate(self.energies)
        ]
        return {
            "available": True,
            "source_state": ref,
            "target_state": target,
            "transition_ev": float(abs(self.energies[target] - self.energies[ref]) * HARTREE_EV),
            "states": states,
            "nto_pairs": pairs,
            "nto_participation_ratio": nto["participation_ratio"],
            "transition_density_norm": nto["sum_weights"],
            "n_promoted": ad["n_promoted"],
            "n_attach": ad["n_attach"],
            "n_detach": ad["n_detach"],
            "json_file": self.json_name,
            "molden_file": self.molden_name,
        }

    def cube(self, kind: str, ref: int, target: int, rank: int = 0) -> str:
        if kind not in KINDS:
            raise ValueError(f"unknown excited-state map: {kind}")
        self._check_pair(ref, target)
        transition = f"S{ref} to S{target}"
        if kind.startswith("nto_"):
            nto = self.nto(ref, target)
            if not (0 <= rank < len(nto["weights"])):
                raise ValueError("NTO pair index out of range")
            key = "holes_ao" if kind == "nto_hole" else "particles_ao"
            role = "hole" if kind == "nto_hole" else "particle"
            return molden.orbital_coeff_cube(
                self.molden_data, nto[key][:, rank],
                f"OQP Studio {transition} {role} NTO {rank + 1}",
                f"weight={float(nto['fractions'][rank]):.8f}",
            )

        if kind in {"attachment", "detachment"}:
            ad = self.attachment_detachment(ref, target)
            matrix = ad[f"{kind}_mo"]
        elif kind == "difference":
            matrix = self.difference_mo(ref, target)
        elif kind == "transition":
            tdm = self.tdm_mo(ref, target)
            matrix = 0.5 * (tdm + tdm.T)
        else:
            matrix = (
                np.diag(self.reference_occupations)
                + self.transition_densities[:, :, target, target]
            )
        return molden.matrix_density_cube(
            self.molden_data, self.ao_density(matrix),
            f"OQP Studio {transition} {kind.replace('_', ' ')}",
            "MRSF physical-root density (e/bohr^3)",
        )


def summary(paths: list[Path], ref: int = 0, target: int = 1) -> dict:
    try:
        return ExcitedStateData.load(paths).summary(ref, target)
    except AnalysisUnavailable as exc:
        return {"available": False, "reason": str(exc)}
