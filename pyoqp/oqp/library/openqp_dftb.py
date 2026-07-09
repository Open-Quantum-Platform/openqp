"""OpenQP runtime adapter for the optional OpenQP-DFTB backend."""

from __future__ import annotations

import ctypes
from dataclasses import dataclass
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import textwrap

import numpy as np

import oqp
from oqp.periodic_table import ELEMENTS_NAME
from oqp.utils.constants import ANGSTROM_TO_BOHR as BOHR_TO_ANGSTROM
from oqp.utils.file_utils import dump_log


_EXCITED_METHOD_BY_TDHF_TYPE = {
    "rpa": "tddftb",
    "tda": "tddftb",
    "sf": "sf",
    "mrsf": "mrsf",
}

_GROUND_TYPES = {"ground", "dftb", "dftb0", "ground_noscc", "noscc"}


@dataclass(frozen=True)
class _StateResult:
    state: int
    reference_energy: float
    state_energy: float
    gradient_bohr: np.ndarray
    excitation_energy: float | None = None
    spin_square: float | None = None
    stdout: str = ""


class OpenQPDFTBAdapter:
    """Make OpenQP-DFTB look like a normal OpenQP energy/gradient provider."""

    def __init__(self, mol):
        self.mol = mol
        self.config = mol.config
        self.dftb = self.config.get("dftb", {})
        self.natom = int(mol.data["natom"])
        self.nstate = self._effective_nstate()
        if not hasattr(mol, "_openqp_dftb_cache"):
            mol._openqp_dftb_cache = {}

    def energy(self):
        """Return an OpenQP-style state-energy array and update molecule data."""
        if self._resolved_method() in _GROUND_TYPES:
            states = [0]
        else:
            states = range(1, self.nstate + 1)
        energies, _ = self._evaluate(states, need_grad=False)
        self.mol.energies = energies
        self._store_td_energies(energies)
        return energies

    def reference(self):
        """Return the DFTB reference energy in the same shape as SCF reference()."""
        method = self._resolved_method()
        if method in _GROUND_TYPES:
            result = self._run_state("ground", 0, need_grad=False)
        else:
            result = self._run_state(method, 1, need_grad=False)
        energies = np.array([result.reference_energy], dtype=float)
        self.mol.energies = energies
        return energies

    def excitation(self, _ref_energy=None):
        """Return DFTB total state energies after a reference-style call."""
        return self.energy()

    def gradient(self, states):
        """Return an OpenQP-style gradient array indexed by state number."""
        states = sorted({int(state) for state in states})
        if not states:
            states = [0]
        energies, gradient_map = self._evaluate(states, need_grad=True)

        nrows = max(len(energies), max(states) + 1)
        grads = np.zeros((nrows, self.natom, 3))
        for state in states:
            grads[state] = gradient_map[state]

        self.mol.energies = energies
        self.mol.grads = grads
        self._store_td_energies(energies)
        return grads

    def _evaluate(self, states, *, need_grad):
        method = self._resolved_method()
        states = sorted({int(state) for state in states})

        if method in _GROUND_TYPES:
            requested = [0]
        else:
            requested = [state for state in states if state > 0]
            if not requested:
                requested = list(range(1, self.nstate + 1))
            requested = sorted(set(requested) | set(range(1, self.nstate + 1)))

        results = {}
        for state in requested:
            results[state] = self._run_state(method, state, need_grad=False)

        if method in _GROUND_TYPES or 0 in states:
            results[0] = self._run_state("ground", 0, need_grad=False)

        if method in _GROUND_TYPES:
            energies = np.array([results[0].state_energy], dtype=float)
        else:
            energies = np.zeros(self.nstate + 1, dtype=float)
            first = results[min(results)]
            energies[0] = first.reference_energy
            for state, result in results.items():
                if state > 0 and state <= self.nstate:
                    energies[state] = result.state_energy

        gradient_map = {}
        if need_grad:
            for state in states:
                run_method = "ground" if state == 0 else method
                gradient_result = self._run_state(run_method, state, need_grad=True)
                results[state] = gradient_result
                gradient_map[state] = gradient_result.gradient_bohr

        return energies, gradient_map

    def _run_state(self, method: str, state: int, *, need_grad: bool) -> _StateResult:
        key = self._cache_key(method, state, need_grad=need_grad)
        cache = self.mol._openqp_dftb_cache
        if key in cache:
            return cache[key]

        backend = str(self.dftb.get("backend", "native")).lower()
        if backend in {"native", "auto"}:
            result = self._run_native(method, state, need_grad=need_grad)
        elif backend == "probe":
            result = self._run_probe(method, state)
        else:
            raise ValueError(f"Unknown OpenQP-DFTB backend: {backend}")

        cache[key] = result
        return result

    def _run_native(self, method: str, state: int, *, need_grad: bool) -> _StateResult:
        """Call libopenqp_dftb_c (the standalone openqp-dftb shared C API) in-process.

        openqp-dftb stays a fully separate library: PyOQP hands over geometry and
        options as plain C scalars/arrays and receives energies and the analytic
        state gradient back. No OpenQP build coupling, no Fortran module seam.
        """
        lib = self._native_library()
        natom = self.natom
        atoms = np.ascontiguousarray(
            np.asarray(self.mol.get_atoms(), dtype=np.int64).reshape(-1)
        )
        coords_bohr = np.ascontiguousarray(
            np.asarray(self.mol.get_system(), dtype=np.float64).reshape(-1)
        )
        parameter = self._parameter_path().encode("utf-8")
        method_name = self._probe_method_name(method).encode("ascii")

        reference_energy = ctypes.c_double()
        state_energy = ctypes.c_double()
        excitation_energy = ctypes.c_double()
        spin_square = ctypes.c_double()
        gradient = (ctypes.c_double * (3 * natom))()
        status_message = ctypes.create_string_buffer(1024)
        status = ctypes.c_int64()

        lib.openqp_dftb_state_gradient(
            ctypes.c_int64(natom),
            atoms.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
            coords_bohr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            parameter,
            ctypes.c_int32(len(parameter)),
            method_name,
            ctypes.c_int32(len(method_name)),
            ctypes.c_int64(int(state)),
            ctypes.c_int64(int(max(1, self.nstate))),
            ctypes.c_int64(int(bool(need_grad))),
            ctypes.c_int64(int(self.config.get("input", {}).get("charge", 0))),
            ctypes.c_int64(self._scc_mixer_code()),
            ctypes.c_int64(int(self.dftb.get("scc_history", 12))),
            ctypes.c_int64(int(self.dftb.get("max_scc_iterations", 1200))),
            ctypes.c_int64(int(self.dftb.get("response_max_iterations", 50))),
            ctypes.c_int64(int(self.dftb.get("response_max_subspace", 100))),
            ctypes.c_int64(self._response_solver_code()),
            ctypes.c_int64(self._reference_multiplicity(method)),
            ctypes.c_int64(int(self.dftb.get("target_multiplicity", 1))),
            ctypes.c_int64(int(bool(self.dftb.get("spin_complete", True)))),
            ctypes.c_int64(int(bool(self.dftb.get("lc_ground_state", False)))),
            ctypes.c_int64(int(self._lc_gamma_is_erf())),
            ctypes.c_int64(int(bool(self.dftb.get("zvector", True)))),
            ctypes.c_double(float(self.dftb.get("scc_tolerance", 1.0e-8))),
            ctypes.c_double(float(self.dftb.get("scc_mixing", 0.35))),
            ctypes.c_double(float(self.dftb.get("scc_max_step", 0.5))),
            ctypes.c_double(float(self.dftb.get("response_tolerance", 1.0e-6))),
            ctypes.c_double(self._spc_channel("spc_coco")),
            ctypes.c_double(self._spc_channel("spc_ovov")),
            ctypes.c_double(self._spc_channel("spc_coov")),
            ctypes.c_double(float(self.dftb.get("omega", 0.3))),
            ctypes.c_double(float(self.dftb.get("cam_alpha", 0.0))),
            ctypes.c_double(float(self.dftb.get("cam_beta", 1.0))),
            ctypes.c_double(float(self.dftb.get("mrsf_shift_oo", 0.0))),
            ctypes.c_double(float(self.dftb.get("mrsf_shift_co", 0.0))),
            ctypes.c_double(float(self.dftb.get("mrsf_shift_ov", 0.0))),
            ctypes.c_double(float(self.dftb.get("mrsf_shift_cv", 0.0))),
            ctypes.byref(reference_energy),
            ctypes.byref(state_energy),
            ctypes.byref(excitation_energy),
            ctypes.byref(spin_square),
            gradient,
            status_message,
            ctypes.c_int32(1024),
            ctypes.byref(status),
        )
        if status.value != 0:
            detail = status_message.value.decode("utf-8", errors="replace").strip()
            suffix = f": {detail}" if detail else ""
            raise RuntimeError(
                f"openqp-dftb native library call failed for method={method}, "
                f"state={state}, status={status.value}{suffix}"
            )

        gradient_bohr = np.frombuffer(gradient, dtype=np.float64).reshape((natom, 3)).copy()
        return _StateResult(
            state=state,
            reference_energy=float(reference_energy.value),
            state_energy=float(state_energy.value),
            gradient_bohr=gradient_bohr,
            excitation_energy=float(excitation_energy.value) if state > 0 else None,
            spin_square=float(spin_square.value) if state > 0 else None,
        )

    def _native_library(self):
        cache = self.mol._openqp_dftb_cache
        lib = cache.get("__native_library__")
        if lib is not None:
            return lib
        path = self._native_library_path()
        try:
            lib = ctypes.CDLL(str(path))
        except OSError as exc:
            raise RuntimeError(f"Could not load the openqp-dftb library {path}: {exc}") from exc
        if not hasattr(lib, "openqp_dftb_state_gradient"):
            raise RuntimeError(
                f"{path} does not export openqp_dftb_state_gradient; "
                "rebuild openqp-dftb with OPENQP_DFTB_BUILD_SHARED=ON."
            )
        lib.openqp_dftb_state_gradient.restype = None
        cache["__native_library__"] = lib
        return lib

    def _native_library_path(self) -> Path:
        raw = self.dftb.get("library_path") or os.environ.get("OPENQP_DFTB_LIBRARY")
        if raw:
            path = self._resolve_user_path(raw)
            if Path(path).exists():
                return Path(path)
            raise FileNotFoundError(f"openqp-dftb library not found at {path}")
        # pip-installed openqp-dftb wheel bundles the library next to its
        # locator package: the most robust default when nothing is configured.
        try:
            import openqp_dftb  # noqa: PLC0415

            return Path(openqp_dftb.library_path())
        except (ImportError, FileNotFoundError):
            pass
        oqp_root = os.environ.get("OPENQP_ROOT", "")
        for name in ("libopenqp_dftb_c.dylib", "libopenqp_dftb_c.so"):
            if oqp_root:
                staged = Path(oqp_root) / "lib" / name
                if staged.exists():
                    return staged
            found = shutil.which(name)
            if found:
                return Path(found)
        message = (
            "openqp-dftb not found: could not locate libopenqp_dftb_c. "
            "Install it with `pip install openqp-dftb` (or `pip install "
            "git+https://github.com/Open-Quantum-Platform/openqp-dftb.git`), "
            "set [dftb] library_path / OPENQP_DFTB_LIBRARY, or build OpenQP "
            "with -DENABLE_OPENQP_DFTB=ON to stage it next to liboqp."
        )
        dump_log(
            self.mol,
            title="PyOQP: openqp-dftb NOT FOUND\n   "
            + "\n   ".join(textwrap.wrap(message, width=76)),
        )
        raise FileNotFoundError(message)

    def _run_probe(self, method: str, state: int) -> _StateResult:
        executable = self._probe_executable()
        parameter_path = self._parameter_path()
        with tempfile.TemporaryDirectory(prefix="openqp-dftb-") as tmpdir:
            xyz_path = Path(tmpdir) / "molecule.xyz"
            self._write_xyz(xyz_path)
            cmd = [
                executable,
                parameter_path,
                str(xyz_path),
                self._probe_method_name(method),
                str(max(1, int(state))),
                str(max(1, int(self.nstate))),
                self._fmt(self.dftb.get("scc_tolerance", 1.0e-8)),
                str(self.dftb.get("scc_mixer", "auto")).lower(),
                self._fmt(self.dftb.get("scc_mixing", 0.35)),
                str(int(self.dftb.get("scc_history", 12))),
                self._fmt(self.dftb.get("scc_max_step", 0.5)),
                str(int(self.dftb.get("max_scc_iterations", 1200))),
                self._fmt(self.dftb.get("spc", 0.5)),
                self._fmt(self.dftb.get("omega", 0.3)),
                self._fmt(self.dftb.get("cam_alpha", 0.0)),
                self._fmt(self.dftb.get("cam_beta", 1.0)),
                str(int(self.config.get("input", {}).get("charge", 0))),
                "erf" if self._lc_gamma_is_erf() else "yukawa",
            ]
            dump_log(
                self.mol,
                title="PyOQP: OpenQP-DFTB state calculation",
                section="dftb",
                info=[method, state, parameter_path],
            )
            completed = subprocess.run(
                cmd,
                check=False,
                capture_output=True,
                text=True,
                timeout=int(self.dftb.get("timeout", 300)),
            )

        if completed.returncode != 0:
            tail = "\n".join((completed.stdout + completed.stderr).splitlines()[-20:])
            raise RuntimeError(
                "OpenQP-DFTB probe failed for "
                f"method={method}, state={state}, rc={completed.returncode}\n{tail}"
            )
        return self._parse_probe_output(state, completed.stdout)

    def _parse_probe_output(self, state: int, stdout: str) -> _StateResult:
        reference_energy = None
        state_energy = None
        excitation_energy = None
        spin_square = None
        gradient = np.zeros((self.natom, 3), dtype=float)

        for raw_line in stdout.splitlines():
            fields = raw_line.split()
            if not fields:
                continue
            key = fields[0]
            if key == "reference_total_energy_hartree" and len(fields) >= 2:
                reference_energy = float(fields[1])
            elif key == "state_energy_hartree" and len(fields) >= 2:
                state_energy = float(fields[1])
            elif key == "state_excitation_energy_hartree" and len(fields) >= 2:
                excitation_energy = float(fields[1])
            elif key == "state_spin_square" and len(fields) >= 2:
                spin_square = float(fields[1])
            elif key == "gradient_hartree_per_angstrom" and len(fields) >= 5:
                atom_index = int(fields[1]) - 1
                gradient[atom_index, :] = [float(x) * BOHR_TO_ANGSTROM for x in fields[2:5]]

        if reference_energy is None:
            reference_energy = state_energy
        if state_energy is None:
            state_energy = reference_energy
        if state_energy is None:
            raise RuntimeError("OpenQP-DFTB probe output did not contain an energy.")

        return _StateResult(
            state=state,
            reference_energy=float(reference_energy),
            state_energy=float(state_energy),
            gradient_bohr=gradient,
            excitation_energy=excitation_energy,
            spin_square=spin_square,
            stdout=stdout,
        )

    def _resolved_method(self) -> str:
        explicit = str(self.dftb.get("type", "auto")).strip().lower()
        if explicit and explicit != "auto":
            return explicit

        runtype = str(self.config.get("input", {}).get("runtype", "energy")).lower()
        istate = int(self.config.get("optimize", {}).get("istate", 0))
        if runtype in {"optimize", "mep"} and istate == 0:
            return "ground"

        td_type = str(self.config.get("tdhf", {}).get("type", "tda")).lower()
        try:
            return _EXCITED_METHOD_BY_TDHF_TYPE[td_type]
        except KeyError as exc:
            raise ValueError(
                f"Cannot derive OpenQP-DFTB response type from tdhf.type={td_type!r}."
            ) from exc

    def _effective_nstate(self) -> int:
        config = self.config
        requested = [int(config.get("tdhf", {}).get("nstate", 1))]
        runtype = str(config.get("input", {}).get("runtype", "energy")).lower()
        if runtype in {"grad", "data"}:
            requested.extend(int(x) for x in config.get("properties", {}).get("grad", [0]))
        if runtype in {"optimize", "mep"}:
            requested.append(int(config.get("optimize", {}).get("istate", 0)))
        if runtype == "meci":
            requested.append(int(config.get("optimize", {}).get("istate", 1)))
            requested.append(int(config.get("optimize", {}).get("jstate", 2)))
        return max(0, max(requested))

    def _cache_key(self, method: str, state: int, *, need_grad: bool):
        coords = np.ascontiguousarray(self.mol.get_system(), dtype=np.float64)
        atoms = tuple(int(z) for z in np.asarray(self.mol.get_atoms()).reshape(-1))
        option_key = (
            self._parameter_path(),
            method,
            state,
            bool(need_grad),
            self.nstate,
            int(self.config.get("input", {}).get("charge", 0)),
            float(self.dftb.get("scc_tolerance", 1.0e-8)),
            str(self.dftb.get("scc_mixer", "auto")).lower(),
            float(self.dftb.get("scc_mixing", 0.35)),
            int(self.dftb.get("scc_history", 12)),
            float(self.dftb.get("scc_max_step", 0.5)),
            int(self.dftb.get("max_scc_iterations", 1200)),
            float(self.dftb.get("spc", 0.5)),
            self._spc_channel("spc_coco"),
            self._spc_channel("spc_ovov"),
            self._spc_channel("spc_coov"),
            float(self.dftb.get("omega", 0.3)),
            float(self.dftb.get("cam_alpha", 0.0)),
            float(self.dftb.get("cam_beta", 1.0)),
            float(self.dftb.get("mrsf_shift_oo", 0.0)),
            float(self.dftb.get("mrsf_shift_co", 0.0)),
            float(self.dftb.get("mrsf_shift_ov", 0.0)),
            float(self.dftb.get("mrsf_shift_cv", 0.0)),
            str(self.dftb.get("lc_gamma", "yukawa")).lower(),
            bool(self.dftb.get("lc_ground_state", False)),
            bool(self.dftb.get("zvector", True)),
            bool(self.dftb.get("spin_complete", True)),
            int(self.dftb.get("reference_multiplicity", 0)),
            int(self.dftb.get("target_multiplicity", 1)),
            str(self.dftb.get("response_solver", "auto")).lower(),
            int(self.dftb.get("response_max_subspace", 100)),
        )
        return atoms, coords.tobytes(), option_key

    def _parameter_path(self) -> str:
        raw = self.dftb.get("parameter_path") or os.environ.get("OPENQP_DFTB_PARAMETER_PATH")
        if not raw:
            raise ValueError("Set [dftb] parameter_path or OPENQP_DFTB_PARAMETER_PATH.")
        return str(self._resolve_user_path(raw))

    def _probe_executable(self) -> str:
        raw = self.dftb.get("executable") or os.environ.get("OPENQP_DFTB_STATE_GRADIENT_PROBE")
        candidates = []
        if raw:
            candidates.append(self._resolve_user_path(raw))
            found = shutil.which(str(raw))
            if found:
                candidates.append(Path(found))
        found = shutil.which("openqp_dftb_state_gradient_probe")
        if found:
            candidates.append(Path(found))

        for candidate in candidates:
            if candidate and Path(candidate).exists():
                return str(candidate)
        raise FileNotFoundError(
            "Could not find openqp_dftb_state_gradient_probe. Set [dftb] executable "
            "or OPENQP_DFTB_STATE_GRADIENT_PROBE for the probe fallback."
        )

    def _resolve_user_path(self, raw_path) -> Path:
        path = Path(str(raw_path)).expanduser()
        if path.is_absolute():
            return path
        input_file = getattr(self.mol, "input_file", "") or ""
        base = Path(input_file).resolve().parent if input_file else Path.cwd()
        resolved = base / path
        return resolved if resolved.exists() else Path.cwd() / path

    def _write_xyz(self, path: Path) -> None:
        atoms = np.asarray(self.mol.get_atoms(), dtype=int).reshape(-1)
        coords_angstrom = np.asarray(self.mol.get_system(), dtype=float).reshape((-1, 3)) * BOHR_TO_ANGSTROM
        lines = [str(len(atoms)), "OpenQP-DFTB geometry from OpenQP"]
        for atomic_number, xyz in zip(atoms, coords_angstrom):
            symbol = ELEMENTS_NAME[int(atomic_number)].strip()
            lines.append(f"{symbol:<2s} {xyz[0]: .12f} {xyz[1]: .12f} {xyz[2]: .12f}")
        path.write_text("\n".join(lines) + "\n", encoding="ascii")

    @staticmethod
    def _probe_method_name(method: str) -> str:
        method = method.lower()
        if method in {"dftb", "dftb0", "ground_noscc", "noscc"}:
            return method
        if method == "sf-tddftb":
            return "sf"
        if method == "mrsf-tddftb":
            return "mrsf"
        return method

    @staticmethod
    def _fmt(value) -> str:
        return f"{float(value):.16g}"

    def _response_solver_code(self) -> int:
        solver = str(self.dftb.get("response_solver", "auto")).lower()
        return {"dense": 0, "davidson": 1, "auto": 2}[solver]

    def _lc_gamma_is_erf(self) -> bool:
        kernel = str(self.dftb.get("lc_gamma", "yukawa")).lower()
        if kernel not in {"yukawa", "erf"}:
            raise ValueError(f"[dftb] lc_gamma must be yukawa or erf, got {kernel!r}")
        return kernel == "erf"

    def _reference_multiplicity(self, method: str) -> int:
        explicit = int(self.dftb.get("reference_multiplicity", 0))
        if explicit > 0:
            return explicit
        return 3 if method in {"sf", "sftddftb", "sf-tddftb", "mrsf", "mrsftddftb", "mrsf-tddftb"} else 1

    def _spc_channel(self, key: str) -> float:
        value = self.dftb.get(key, None)
        if value is None or str(value).strip() == "" or float(value) <= -900.0:
            return float(self.dftb.get("spc", 0.5))
        return float(value)

    def _scc_mixer_code(self) -> int:
        mixer = str(self.dftb.get("scc_mixer", "auto")).lower()
        return {
            "linear": 0,
            "anderson": 1,
            "pulay": 1,
            "broyden": 2,
            "auto": 3,
            "diis": 4,
            "trust": 5,
            "trah": 5,
        }[mixer]

    def _store_td_energies(self, energies) -> None:
        energies = np.asarray(energies, dtype=float)
        if energies.size > 1:
            self.mol.data["OQP::td_energies"] = energies[1:] - energies[0]
