"""Configuration and thin C-ABI bridge for native Fortran ODP umbrellas."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import warnings

import numpy as np

import oqp


_CV_TYPES = {
    "distance": (1, 2),
    "asymmetric_distance": (2, 4),
    "asymmetric-distance": (2, 4),
    "asym_distance": (2, 4),
    "angle": (3, 3),
}
_CV_CANONICAL_NAMES = {1: "distance", 2: "asymmetric_distance", 3: "angle"}
_CV_PATTERN = re.compile(r"([A-Za-z][A-Za-z_-]*)\s*\(([^()]*)\)")
KB_HARTREE_PER_KELVIN = 3.166811563e-6


def parse_odp_cv_specification(specification):
    """Parse ``distance(1,2); angle(3,4,5)`` using 1-based atom labels."""
    text = str(specification or "").strip()
    if not text:
        raise ValueError("[odp] cv must define at least one continuous CV")
    pieces = [piece.strip() for piece in text.split(";") if piece.strip()]
    types = []
    atoms = []
    labels = []
    for piece in pieces:
        match = _CV_PATTERN.fullmatch(piece)
        if match is None:
            raise ValueError(
                "[odp] cv entries must look like distance(1,2), "
                "asymmetric_distance(1,2,3,4), or angle(1,2,3)"
            )
        name = match.group(1).strip().lower()
        if name not in _CV_TYPES:
            raise ValueError(f"[odp] unsupported CV type {name!r}")
        cv_type, count = _CV_TYPES[name]
        try:
            indices = [
                int(token) for token in re.split(r"[\s,]+", match.group(2).strip())
                if token
            ]
        except ValueError as exc:
            raise ValueError(f"[odp] non-integer atom index in {piece!r}") from exc
        if len(indices) != count or min(indices, default=0) < 1:
            raise ValueError(
                f"[odp] {_CV_CANONICAL_NAMES[cv_type]} requires {count} "
                "positive, 1-based atom indices"
            )
        if cv_type == 1 and indices[0] == indices[1]:
            raise ValueError("[odp] a distance requires two distinct atoms")
        if cv_type == 2 and (indices[0] == indices[1] or indices[2] == indices[3]):
            raise ValueError("[odp] each asymmetric-distance pair must be distinct")
        if cv_type == 2 and frozenset(indices[:2]) == frozenset(indices[2:]):
            raise ValueError(
                "[odp] the two asymmetric-distance atom pairs must differ"
            )
        if cv_type == 3 and len(set(indices)) != 3:
            raise ValueError("[odp] an angle requires three distinct atoms")
        row = [-1, -1, -1, -1]
        row[:count] = [index - 1 for index in indices]
        types.append(cv_type)
        atoms.append(row)
        labels.append(f"{_CV_CANONICAL_NAMES[cv_type]}({','.join(map(str, indices))})")
    return (
        np.ascontiguousarray(types, dtype=np.int32),
        np.ascontiguousarray(atoms, dtype=np.int64),
        tuple(labels),
    )


class ODPUmbrella:
    """One configured ODP window evaluated by the resident Fortran kernel."""

    def __init__(self, section):
        self.enabled = self._as_bool(section.get("enabled", False))
        if not self.enabled:
            raise ValueError("ODPUmbrella requires [odp] enabled=true")
        self.cv_types, self.cv_atoms, self.cv_labels = parse_odp_cv_specification(
            section.get("cv", "")
        )
        self.ncv = len(self.cv_types)
        self.scales = self._vector(section.get("scale", ()), "scale")
        self.reference_r = self._vector(
            section.get("reference_r", ()), "reference_r"
        )
        self.reference_p = self._vector(
            section.get("reference_p", ()), "reference_p"
        )
        for name, value in (
            ("scale", self.scales),
            ("reference_r", self.reference_r),
            ("reference_p", self.reference_p),
        ):
            if value.size != self.ncv:
                raise ValueError(
                    f"[odp] {name} requires exactly {self.ncv} value(s), one per CV"
                )
        if np.any(self.scales <= 0.0):
            raise ValueError("[odp] every scale must be positive")
        self.center = self._finite_float(section.get("center", 0.0), "center")
        self.k_parallel = self._finite_float(
            section.get("k_parallel", 0.0), "k_parallel"
        )
        self.k_perpendicular = self._finite_float(
            section.get("k_perpendicular", 0.0), "k_perpendicular"
        )
        if self.k_parallel <= 0.0:
            raise ValueError("[odp] k_parallel must be positive when ODP is enabled")
        if self.k_perpendicular < 0.0:
            raise ValueError("[odp] k_perpendicular must be non-negative")
        self.window = int(section.get("window", 0))
        if self.window < 0:
            raise ValueError("[odp] window must be non-negative")
        direction = self.scales*(self.reference_p - self.reference_r)
        self.path_length = float(np.linalg.norm(direction))
        if not math.isfinite(self.path_length) or self.path_length <= 1.0e-12:
            raise ValueError("[odp] scaled reactant and product references must differ")
        self.last = None

    @staticmethod
    def _as_bool(value):
        return value is True or str(value).strip().lower() in {
            "true", "t", "1", "yes", "on", ".true.",
        }

    @staticmethod
    def _finite_float(value, name):
        result = float(value)
        if not math.isfinite(result):
            raise ValueError(f"[odp] {name} must be finite")
        return result

    def _vector(self, value, name):
        if isinstance(value, str):
            try:
                value = [item.strip() for item in value.split(",") if item.strip()]
            except AttributeError as exc:
                raise ValueError(f"[odp] invalid {name}") from exc
        try:
            array = np.ascontiguousarray(value, dtype=np.float64).reshape(-1)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"[odp] {name} must be a comma-separated numeric vector") from exc
        if not np.all(np.isfinite(array)):
            raise ValueError(f"[odp] {name} must contain only finite values")
        return array

    def validate_atom_count(self, natom):
        natom = int(natom)
        used = self.cv_atoms[self.cv_atoms >= 0]
        if natom <= 0 or used.size == 0 or int(used.max()) >= natom:
            largest = int(used.max()) + 1 if used.size else 0
            raise ValueError(
                f"[odp] CV atom {largest} is outside the {natom}-atom MD system"
            )

    def evaluate(self, coordinates):
        """Return native bias energy, force, CVs, and projection diagnostics."""
        xyz = np.ascontiguousarray(coordinates, dtype=np.float64)
        if xyz.ndim != 2 or xyz.shape[1] != 3:
            raise ValueError("ODP coordinates must have shape (natom, 3)")
        self.validate_atom_count(len(xyz))
        if not np.all(np.isfinite(xyz)):
            raise RuntimeError("ODP refused non-finite MD coordinates")
        energy = np.zeros(1, dtype=np.float64)
        force = np.zeros_like(xyz)
        xi = np.zeros(1, dtype=np.float64)
        raw = np.zeros(self.ncv, dtype=np.float64)
        scaled = np.zeros(self.ncv, dtype=np.float64)
        perpendicular = np.zeros(self.ncv, dtype=np.float64)
        perpendicular_norm = np.zeros(1, dtype=np.float64)
        energy_parallel = np.zeros(1, dtype=np.float64)
        energy_perpendicular = np.zeros(1, dtype=np.float64)
        status = oqp.oqp_odp_umbrella_evaluate(
            len(xyz), self.ncv,
            oqp.ffi.cast("int *", self.cv_types.ctypes.data),
            oqp.ffi.cast("int64_t *", self.cv_atoms.ctypes.data),
            oqp.ffi.cast("double *", self.scales.ctypes.data),
            oqp.ffi.cast("double *", self.reference_r.ctypes.data),
            oqp.ffi.cast("double *", self.reference_p.ctypes.data),
            self.center, self.k_parallel, self.k_perpendicular,
            oqp.ffi.cast("double *", xyz.ctypes.data),
            oqp.ffi.cast("double *", energy.ctypes.data),
            oqp.ffi.cast("double *", force.ctypes.data),
            oqp.ffi.cast("double *", xi.ctypes.data),
            oqp.ffi.cast("double *", raw.ctypes.data),
            oqp.ffi.cast("double *", scaled.ctypes.data),
            oqp.ffi.cast("double *", perpendicular.ctypes.data),
            oqp.ffi.cast("double *", perpendicular_norm.ctypes.data),
            oqp.ffi.cast("double *", energy_parallel.ctypes.data),
            oqp.ffi.cast("double *", energy_perpendicular.ctypes.data),
        )
        if status != 0:
            meanings = {
                1: "invalid dimensions or force constants",
                2: "invalid CV type or atom indices",
                3: "invalid scale or reference",
                4: "singular distance CV",
                5: "zero-length angle bond",
                6: "collinear angle CV",
                7: "degenerate scaled reactant/product path",
                8: "non-finite input or output",
            }
            raise RuntimeError(
                f"native ODP evaluation failed ({status}: "
                f"{meanings.get(status, 'unknown status')})"
            )
        result = {
            "energy": float(energy[0]),
            "force": force,
            "xi": float(xi[0]),
            "cv_raw": raw,
            "cv_scaled": scaled,
            "cv_perpendicular": perpendicular,
            "perpendicular_norm": float(perpendicular_norm[0]),
            "energy_parallel": float(energy_parallel[0]),
            "energy_perpendicular": float(energy_perpendicular[0]),
        }
        if not all(
            np.all(np.isfinite(value)) if isinstance(value, np.ndarray)
            else math.isfinite(value)
            for value in result.values()
        ):
            raise RuntimeError("native ODP evaluation returned non-finite data")
        self.last = result
        return result

    def provenance(self):
        """JSON-safe definition sufficient to reconstruct this umbrella bias."""
        return {
            "enabled": True,
            "window": self.window,
            "cv": list(self.cv_labels),
            "cv_atom_indexing": "1-based",
            "cv_native_units": [
                "radian" if cv_type == 3 else "bohr"
                for cv_type in self.cv_types
            ],
            "scale": self.scales.tolist(),
            "reference_r": self.reference_r.tolist(),
            "reference_p": self.reference_p.tolist(),
            "scaled_path_length": self.path_length,
            "center": self.center,
            "k_parallel_hartree": self.k_parallel,
            "k_perpendicular_hartree": self.k_perpendicular,
            "projection": "signed_scaled_dot_over_path_norm_squared",
            "perpendicular_restraint": self.k_perpendicular > 0.0,
        }

    def signature(self):
        return json.dumps(self.provenance(), sort_keys=True, separators=(",", ":"))


def odp_from_config(config):
    section = config.get("odp", {})
    if not ODPUmbrella._as_bool(section.get("enabled", False)):
        return None
    return ODPUmbrella(section)


def _logsumexp(values, axis=None):
    values = np.asarray(values, dtype=np.float64)
    maximum = np.max(values, axis=axis, keepdims=True)
    shifted = np.exp(values - maximum)
    result = maximum + np.log(np.sum(shifted, axis=axis, keepdims=True))
    if axis is None:
        return float(result.reshape(-1)[0])
    return np.squeeze(result, axis=axis)


def _file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _same_odp_coordinate(left, right):
    exact_keys = ("cv", "cv_atom_indexing", "cv_native_units", "projection")
    numeric_keys = (
        "scale", "reference_r", "reference_p", "scaled_path_length"
    )
    return all(left.get(key) == right.get(key) for key in exact_keys) and all(
        np.array_equal(
            np.asarray(left.get(key), dtype=np.float64),
            np.asarray(right.get(key), dtype=np.float64),
        )
        for key in numeric_keys
    )


def odp_wham(trajectory_paths, temperature_kelvin, bins=100, *, discard=0,
             stride=1, tolerance=1.0e-10, max_iterations=10000,
             output=None):
    """Recover a one-dimensional unbiased ODP PMF from packed window TRJs.

    The fixed-point equations use every saved sample rather than evaluating a
    window bias only at histogram-bin centres.  This preserves optional
    perpendicular-restraint reweighting.  ``temperature_kelvin`` is mandatory:
    an NVE trajectory's initial velocity temperature is never silently treated
    as a canonical thermostat temperature.
    """
    from oqp.library.namd import read_odp_wham_series

    paths = [os.fspath(path) for path in trajectory_paths]
    if not paths:
        raise ValueError("ODP WHAM requires at least one packed trajectory")
    temperature = float(temperature_kelvin)
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("ODP WHAM temperature_kelvin must be positive and finite")
    discard = int(discard)
    stride = int(stride)
    max_iterations = int(max_iterations)
    tolerance = float(tolerance)
    if discard < 0 or stride < 1:
        raise ValueError("ODP WHAM requires discard >= 0 and stride >= 1")
    if max_iterations < 1 or not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("ODP WHAM requires positive convergence controls")

    resolved_paths = {}
    trajectory_hashes = []
    hashed_paths = {}
    for path in paths:
        resolved = os.path.realpath(os.path.abspath(path))
        if resolved in resolved_paths:
            raise ValueError(
                "ODP WHAM duplicate trajectory input resolves to the same file: "
                f"{resolved_paths[resolved]!r} and {path!r}"
            )
        resolved_paths[resolved] = path
        digest = _file_sha256(path)
        if digest in hashed_paths:
            raise ValueError(
                "ODP WHAM duplicate trajectory content: "
                f"{hashed_paths[digest]!r} and {path!r}"
            )
        hashed_paths[digest] = path
        trajectory_hashes.append(digest)

    loaded = []
    coordinate_provenance = None
    system_identity = None
    ensembles = []
    window_definitions = {}
    for path in paths:
        series = read_odp_wham_series(path)
        provenance = series["provenance"]
        if coordinate_provenance is None:
            coordinate_provenance = provenance
        elif not _same_odp_coordinate(coordinate_provenance, provenance):
            raise ValueError(
                "ODP WHAM trajectories use incompatible CV ordering, metric, or R/P references"
            )
        if system_identity is None:
            system_identity = series["system_identity"]
        elif system_identity != series["system_identity"]:
            raise ValueError(
                "ODP WHAM trajectories come from different molecular systems "
                "or electronic Hamiltonians"
            )
        window = int(provenance["window"])
        definition = (
            float(provenance["center"]),
            float(provenance["k_parallel_hartree"]),
            float(provenance["k_perpendicular_hartree"]),
        )
        if window in window_definitions and window_definitions[window] != definition:
            raise ValueError(f"ODP WHAM window {window} has conflicting bias definitions")
        window_definitions[window] = definition
        selection = slice(discard, None, stride)
        xi = np.asarray(series["xi"][selection], dtype=np.float64)
        perpendicular_norm = np.asarray(
            series["perpendicular_norm"][selection], dtype=np.float64)
        recorded_bias = np.asarray(
            series["bias_hartree"][selection], dtype=np.float64)
        steps = np.asarray(series["step"][selection], dtype=np.int64)
        if xi.size == 0:
            raise ValueError(f"ODP WHAM trajectory {path!r} has no selected samples")
        if not (np.all(np.isfinite(xi))
                and np.all(np.isfinite(perpendicular_norm))
                and np.all(np.isfinite(recorded_bias))):
            raise ValueError(f"ODP WHAM trajectory {path!r} contains non-finite records")
        expected_bias = (
            0.5*definition[1]*(xi - definition[0])**2
            + 0.5*definition[2]*perpendicular_norm**2
        )
        if not np.allclose(recorded_bias, expected_bias, rtol=2.0e-12, atol=2.0e-12):
            raise ValueError(
                f"ODP WHAM trajectory {path!r} bias records disagree with provenance"
            )
        ensembles.append(series.get("ensemble"))
        loaded.append((path, window, xi, perpendicular_norm, steps))

    windows = np.asarray(sorted(window_definitions), dtype=np.int64)
    window_lookup = {window: index for index, window in enumerate(windows)}
    xi = np.concatenate([entry[2] for entry in loaded])
    perpendicular_norm = np.concatenate([entry[3] for entry in loaded])
    sample_steps = np.concatenate([entry[4] for entry in loaded])
    sample_window = np.concatenate([
        np.full(entry[2].size, entry[1], dtype=np.int64) for entry in loaded
    ])
    origins = np.asarray(
        [window_lookup[int(window)] for window in sample_window], dtype=np.int64)
    counts = np.bincount(origins, minlength=len(windows)).astype(np.float64)
    if np.any(counts <= 0.0):
        raise ValueError("ODP WHAM encountered an empty umbrella window")

    bias = np.empty((xi.size, len(windows)), dtype=np.float64)
    for column, window in enumerate(windows):
        center, k_parallel, k_perpendicular = window_definitions[int(window)]
        bias[:, column] = (
            0.5*k_parallel*(xi - center)**2
            + 0.5*k_perpendicular*perpendicular_norm**2
        )
    beta = 1.0/(KB_HARTREE_PER_KELVIN*temperature)
    log_counts = np.log(counts)
    reduced_free_energies = np.zeros(len(windows), dtype=np.float64)
    converged = False
    residual = math.inf
    for iteration in range(1, max_iterations + 1):
        log_denominator = _logsumexp(
            log_counts[None, :] + reduced_free_energies[None, :] - beta*bias,
            axis=1,
        )
        updated = np.empty_like(reduced_free_energies)
        for column in range(len(windows)):
            updated[column] = -_logsumexp(
                -beta*bias[:, column] - log_denominator)
        updated -= updated[0]
        residual = float(np.max(np.abs(updated - reduced_free_energies)))
        reduced_free_energies = updated
        if residual < tolerance:
            converged = True
            break
    if not converged:
        raise RuntimeError(
            f"ODP WHAM did not converge in {max_iterations} iterations "
            f"(residual={residual:.3e})"
        )

    log_denominator = _logsumexp(
        log_counts[None, :] + reduced_free_energies[None, :] - beta*bias,
        axis=1,
    )
    log_weights = -log_denominator
    log_weights -= _logsumexp(log_weights)
    weights = np.exp(log_weights)
    responsibilities = np.exp(
        log_counts[None, :] + reduced_free_energies[None, :]
        - beta*bias - log_denominator[:, None]
    )
    overlap = np.empty((len(windows), len(windows)), dtype=np.float64)
    for row in range(len(windows)):
        overlap[row] = responsibilities[origins == row].mean(axis=0)

    if np.isscalar(bins):
        bin_count = int(bins)
        if bin_count < 2:
            raise ValueError("ODP WHAM bins must be at least 2")
        lower = float(np.min(xi))
        upper = float(np.max(xi))
        if not upper > lower:
            raise ValueError("ODP WHAM requires a nonzero sampled xi range")
        edges = np.linspace(lower, upper, bin_count + 1)
    else:
        edges = np.asarray(bins, dtype=np.float64).reshape(-1)
        if (edges.size < 3 or not np.all(np.isfinite(edges))
                or np.any(np.diff(edges) <= 0.0)):
            raise ValueError("ODP WHAM bin edges must be finite and strictly increasing")
    probability, _ = np.histogram(xi, bins=edges, weights=weights)
    widths = np.diff(edges)
    density = probability/widths
    free_energy = np.full_like(density, np.inf)
    occupied = density > 0.0
    free_energy[occupied] = -np.log(density[occupied])/beta
    if np.any(occupied):
        free_energy[occupied] -= np.min(free_energy[occupied])

    ensemble_warning = None
    normalized_ensembles = sorted({str(value or "unknown") for value in ensembles})
    if normalized_ensembles != ["NVT"]:
        ensemble_warning = (
            "WHAM assumes canonical sampling at the supplied temperature; "
            f"trajectory ensemble metadata is {', '.join(normalized_ensembles)}"
        )
        warnings.warn(ensemble_warning, RuntimeWarning, stacklevel=2)
    result = {
        "temperature_kelvin": temperature,
        "beta_hartree_inverse": beta,
        "windows": windows,
        "window_centers": np.asarray([
            window_definitions[int(window)][0] for window in windows
        ], dtype=np.float64),
        "window_k_parallel_hartree": np.asarray([
            window_definitions[int(window)][1] for window in windows
        ], dtype=np.float64),
        "window_k_perpendicular_hartree": np.asarray([
            window_definitions[int(window)][2] for window in windows
        ], dtype=np.float64),
        "window_counts": counts.astype(np.int64),
        "window_free_energies_hartree": reduced_free_energies/beta,
        "window_overlap": overlap,
        "bin_edges": edges,
        "bin_centers": 0.5*(edges[:-1] + edges[1:]),
        "probability": probability,
        "probability_density": density,
        "free_energy_hartree": free_energy,
        "free_energy_kcal_mol": free_energy*627.5094740631,
        "sample_weights": weights,
        "sample_xi": xi,
        "sample_window": sample_window,
        "sample_step": sample_steps,
        "effective_sample_size": float(1.0/np.sum(weights*weights)),
        "iterations": iteration,
        "residual": residual,
        "converged": converged,
        "coordinate_provenance": {
            key: coordinate_provenance.get(key)
            for key in ("cv", "cv_atom_indexing", "cv_native_units", "scale",
                        "reference_r", "reference_p", "scaled_path_length",
                        "projection")
        },
        "system_identity": system_identity,
        "trajectory_paths": paths,
        "trajectory_sha256": trajectory_hashes,
        "discard": discard,
        "stride": stride,
        "tolerance": tolerance,
        "max_iterations": max_iterations,
        "ensemble_metadata": normalized_ensembles,
        "ensemble_warning": ensemble_warning,
    }
    if output is not None:
        write_odp_wham(output, result)
    return result


def write_odp_wham(path, result):
    """Write a portable NPZ containing WHAM arrays and JSON provenance."""
    array_keys = {
        key for key, value in result.items() if isinstance(value, np.ndarray)
    }
    metadata = {
        key: value for key, value in result.items() if key not in array_keys
    }
    payload = {key: result[key] for key in sorted(array_keys)}
    payload["metadata_json"] = np.array([
        json.dumps(metadata, sort_keys=True, separators=(",", ":"))
    ])
    np.savez_compressed(os.fspath(path), **payload)
