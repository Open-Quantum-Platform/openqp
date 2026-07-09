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
    # C API v2 extras (native backend only).
    all_state_energies: np.ndarray | None = None
    all_spin_squares: np.ndarray | None = None
    relaxed_charges: np.ndarray | None = None
    nbf: int = 0
    noca: int = 0
    nocb: int = 0
    mo_energies: np.ndarray | None = None
    mo_coefficients: np.ndarray | None = None
    response_vectors: np.ndarray | None = None


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
        # SF/MRSF-TDDFTB uses a high-spin ROKS reference: mark scf.type
        # accordingly so generic bookkeeping (get_data tag selection for the
        # NAMD back_door carry, beta-set handling) treats both spin sets.
        if self._probe_method_name(self._resolved_method()) in {"sf", "mrsf"} and \
                str(self.config.get("scf", {}).get("type", "rhf")).lower() == "rhf":
            self.config.setdefault("scf", {})["type"] = "rohf"

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
            result = self._run_state(method, 0, need_grad=False)
        else:
            result = self._run_state(method, 1, need_grad=False)
            self._store_wavefunction_tags(result)
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
        active = max(states)
        active_result = getattr(self, "_last_results", {}).get(active)
        if active_result is not None:
            if active_result.relaxed_charges is not None:
                self.mol.data["OQP::partial_charges"] = np.ascontiguousarray(
                    active_result.relaxed_charges)
            self._store_wavefunction_tags(active_result)
        return grads

    def _evaluate(self, states, *, need_grad):
        method = self._resolved_method()
        # The ground reference for a state-0 call is the requested ground
        # variant (dftb0/noscc/ground_noscc run NoSCC), not always plain SCC.
        ground_method = method if method in _GROUND_TYPES else "ground"
        states = sorted({int(state) for state in states})

        if method in _GROUND_TYPES:
            requested = [0]
        else:
            requested = [state for state in states if state > 0]
            if not requested:
                requested = list(range(1, self.nstate + 1))
            requested = sorted(set(requested) | set(range(1, self.nstate + 1)))

        results = {}
        if method not in _GROUND_TYPES:
            # C API v2 returns every response root from ONE solve.
            base = self._run_state(method, min(requested), need_grad=False)
            for state in requested:
                results[state] = base

        if method in _GROUND_TYPES or 0 in states:
            results[0] = self._run_state(ground_method, 0, need_grad=False)

        if method in _GROUND_TYPES:
            energies = np.array([results[0].state_energy], dtype=float)
        else:
            energies = np.zeros(self.nstate + 1, dtype=float)
            energies[0] = base.reference_energy
            if base.all_state_energies is not None and len(base.all_state_energies) >= self.nstate:
                energies[1:] = base.all_state_energies[: self.nstate]
            else:
                # probe backend fallback: one call per state
                for state in requested:
                    result = self._run_state(method, state, need_grad=False)
                    results[state] = result
                    if 0 < state <= self.nstate:
                        energies[state] = result.state_energy
            self._store_wavefunction_tags(base)

        gradient_map = {}
        if need_grad:
            for state in states:
                run_method = ground_method if state == 0 else method
                gradient_result = self._run_state(run_method, state, need_grad=True)
                results[state] = gradient_result
                gradient_map[state] = gradient_result.gradient_bohr

        self._last_results = results
        return energies, gradient_map

    def _reject_probe_with_embedding(self, backend: str) -> None:
        if backend == "probe" and self._external_potential() is not None:
            raise ValueError(
                "OpenQP-DFTB QM/MM electrostatic embedding requires the native "
                "backend: the probe executable does not carry the per-atom "
                "external potential, so it would return unembedded energies and "
                "gradients while the driver assembles embedded coupling forces. "
                "Set [dftb] backend=native for QM/MM jobs."
            )

    def _run_state(self, method: str, state: int, *, need_grad: bool) -> _StateResult:
        key = self._cache_key(method, state, need_grad=need_grad)
        cache = self.mol._openqp_dftb_cache
        if key in cache:
            return cache[key]
        # One geometry (+ embedding potential) per cache generation: an MD
        # trajectory would otherwise grow the cache without bound.
        generation = (key[0], key[1], key[2][0])
        if cache.get("__generation__") != generation:
            library = cache.get("__native_library__")
            cache.clear()
            if library is not None:
                cache["__native_library__"] = library
            cache["__generation__"] = generation

        backend = str(self.dftb.get("backend", "native")).lower()
        self._reject_probe_with_embedding(backend)
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

        v_ext = self._external_potential()
        if v_ext is None:
            n_ext = 0
            ext_potential = (ctypes.c_double * 1)()
        else:
            n_ext = natom
            ext_potential = (ctypes.c_double * natom)(*np.asarray(v_ext, dtype=float))

        nstate = int(max(1, self.nstate))
        n_roots_out = ctypes.c_int64()
        all_state_energies = (ctypes.c_double * nstate)()
        all_spin_squares = (ctypes.c_double * nstate)()
        relaxed_charges = (ctypes.c_double * natom)()
        nbf_out = ctypes.c_int64()
        noca_out = ctypes.c_int64()
        nocb_out = ctypes.c_int64()
        vec_dim_out = ctypes.c_int64()
        # Minimal SK basis has at most 9 orbitals/atom (spd): over-allocate and
        # let the library fill what it needs (capacity-guarded on its side).
        nbf_max = 9 * natom
        mo_capacity = nbf_max * nbf_max
        mo_energies = (ctypes.c_double * nbf_max)()
        mo_coefficients = (ctypes.c_double * mo_capacity)()
        vec_capacity = mo_capacity * nstate
        response_vectors = (ctypes.c_double * vec_capacity)()

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
            ctypes.c_int64(n_ext),
            ext_potential,
            ctypes.byref(reference_energy),
            ctypes.byref(state_energy),
            ctypes.byref(excitation_energy),
            ctypes.byref(spin_square),
            gradient,
            ctypes.byref(n_roots_out),
            all_state_energies,
            all_spin_squares,
            relaxed_charges,
            ctypes.byref(nbf_out),
            ctypes.byref(noca_out),
            ctypes.byref(nocb_out),
            ctypes.byref(vec_dim_out),
            ctypes.c_int64(mo_capacity),
            mo_energies,
            mo_coefficients,
            ctypes.c_int64(vec_capacity),
            response_vectors,
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
        n_roots = int(n_roots_out.value)
        nbf = int(nbf_out.value)
        vec_dim = int(vec_dim_out.value)
        all_e = np.frombuffer(all_state_energies, dtype=np.float64)[:n_roots].copy() \
            if n_roots > 0 else None
        all_s2 = np.frombuffer(all_spin_squares, dtype=np.float64)[:n_roots].copy() \
            if n_roots > 0 else None
        mo_c = None
        mo_e = None
        if 0 < nbf <= nbf_max:
            mo_c = np.frombuffer(mo_coefficients, dtype=np.float64)[: nbf * nbf].reshape(
                (nbf, nbf), order="F").copy()
            mo_e = np.frombuffer(mo_energies, dtype=np.float64)[:nbf].copy()
        vecs = None
        if vec_dim > 0 and n_roots > 0:
            vecs = np.frombuffer(response_vectors, dtype=np.float64)[: vec_dim * n_roots].reshape(
                (vec_dim, n_roots), order="F").copy()
        return _StateResult(
            state=state,
            reference_energy=float(reference_energy.value),
            state_energy=float(state_energy.value),
            gradient_bohr=gradient_bohr,
            excitation_energy=float(excitation_energy.value) if state > 0 else None,
            spin_square=float(spin_square.value) if state > 0 else None,
            all_state_energies=all_e,
            all_spin_squares=all_s2,
            relaxed_charges=np.frombuffer(relaxed_charges, dtype=np.float64).copy(),
            nbf=nbf,
            noca=int(noca_out.value),
            nocb=int(nocb_out.value),
            mo_energies=mo_e,
            mo_coefficients=mo_c,
            response_vectors=vecs,
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
        # The -DENABLE_OPENQP_DFTB=ON hook stages libopenqp_dftb_c next to the
        # liboqp that the Python package already resolved. Self-locating installs
        # and source-tree runs do not set OPENQP_ROOT, so check the RESOLVED
        # runtime root's lib dir before the OPENQP_ROOT env fallback and PATH.
        lib_dirs = []
        resolved_root = getattr(oqp, "oqp_root", "")
        if resolved_root:
            lib_dirs.append(Path(resolved_root) / "lib")
        oqp_root = os.environ.get("OPENQP_ROOT", "")
        if oqp_root:
            lib_dirs.append(Path(oqp_root) / "lib")
        for name in ("libopenqp_dftb_c.dylib", "libopenqp_dftb_c.so"):
            for lib_dir in lib_dirs:
                staged = lib_dir / name
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

        # A plain energy/gradient run that targets no excited state and does not
        # explicitly request an open-shell (SF/MRSF) response is a ground-state
        # DFTB job. An explicit tdhf.type=sf/mrsf is honored even for
        # runtype=energy, so it is NOT collapsed to ground here.
        if runtype in {"energy", "grad"} and td_type in {"rpa", "tda"}:
            grad_states = [int(x) for x in self.config.get("properties", {}).get("grad", [])]
            if not any(s > 0 for s in grad_states):
                return "ground"

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

    def states_overlap(self, xyz_old, xyz_new, mo_old, mo_new, x_old, x_new, *,
                       noca, nocb, multiplicity, tlf_order=2):
        """Cross-geometry MRSF/SF state overlap via libopenqp_dftb_c.

        All matrix inputs are FORTRAN-FLAT 1-D arrays (i.e. tag.ravel() of the
        native-convention tags). Returns (overlap_mo_tag, states_overlap_tag)
        in the same native tag layout (transpose view of the Fortran matrices).
        """
        lib = self._native_library()
        natom = self.natom
        nbf = int(round(np.sqrt(mo_new.size)))
        nstate = x_new.size // (noca * (nbf - nocb))
        atoms = np.ascontiguousarray(
            np.asarray(self.mol.get_atoms(), dtype=np.int64).reshape(-1))
        parameter = self._parameter_path().encode("utf-8")
        overlap_mo = (ctypes.c_double * (nbf * nbf))()
        states_overlap = (ctypes.c_double * (nstate * nstate))()
        status_message = ctypes.create_string_buffer(1024)
        status = ctypes.c_int64()
        as_c = lambda a: np.ascontiguousarray(np.asarray(a, dtype=np.float64).reshape(-1))
        xyz_old = as_c(xyz_old); xyz_new = as_c(xyz_new)
        mo_old = as_c(mo_old); mo_new = as_c(mo_new)
        x_old = as_c(x_old); x_new = as_c(x_new)
        dptr = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        lib.openqp_dftb_states_overlap(
            ctypes.c_int64(natom),
            atoms.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
            dptr(xyz_old), dptr(xyz_new),
            parameter, ctypes.c_int32(len(parameter)),
            ctypes.c_int64(nbf), ctypes.c_int64(int(noca)), ctypes.c_int64(int(nocb)),
            ctypes.c_int64(int(nstate)), ctypes.c_int64(int(multiplicity)),
            ctypes.c_int64(int(tlf_order)),
            dptr(mo_old), dptr(mo_new), dptr(x_old), dptr(x_new),
            overlap_mo, states_overlap,
            status_message, ctypes.c_int32(1024), ctypes.byref(status),
        )
        if status.value != 0:
            detail = status_message.value.decode("utf-8", errors="replace").strip()
            raise RuntimeError(f"openqp-dftb states_overlap failed: {detail}")
        # C-order reshape of the Fortran-flat buffers = the native tag layout.
        s_mo = np.frombuffer(overlap_mo, dtype=np.float64).reshape((nbf, nbf)).copy()
        s_st = np.frombuffer(states_overlap, dtype=np.float64).reshape((nstate, nstate)).copy()
        return s_mo, s_st

    def soc_matrix(self, mo, x_singlet, x_triplet, *, noca, nocb):
        """Raw one-center MRSF SOC matrix (no alpha^2/2) via libopenqp_dftb_c.

        Inputs are Fortran-flat 1-D arrays; returns (hsoc_re, hsoc_im) as
        (dim, dim) numpy arrays in the GAMESS state ordering, where the
        row/column index layout matches the native soc_mrsf convention.
        """
        lib = self._native_library()
        natom = self.natom
        nbf = int(round(np.sqrt(mo.size)))
        n_dim = noca * (nbf - nocb)
        ns = x_singlet.size // n_dim
        nt = x_triplet.size // n_dim
        dim = ns + 3 * nt
        atoms = np.ascontiguousarray(
            np.asarray(self.mol.get_atoms(), dtype=np.int64).reshape(-1))
        coords_bohr = np.ascontiguousarray(
            np.asarray(self.mol.get_system(), dtype=np.float64).reshape(-1))
        parameter = self._parameter_path().encode("utf-8")
        hsoc_re = (ctypes.c_double * (dim * dim))()
        hsoc_im = (ctypes.c_double * (dim * dim))()
        status_message = ctypes.create_string_buffer(1024)
        status = ctypes.c_int64()
        as_c = lambda a: np.ascontiguousarray(np.asarray(a, dtype=np.float64).reshape(-1))
        mo = as_c(mo); x_singlet = as_c(x_singlet); x_triplet = as_c(x_triplet)
        dptr = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        lib.openqp_dftb_soc_matrix(
            ctypes.c_int64(natom),
            atoms.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
            dptr(coords_bohr),
            parameter, ctypes.c_int32(len(parameter)),
            ctypes.c_int64(nbf), ctypes.c_int64(int(noca)), ctypes.c_int64(int(nocb)),
            ctypes.c_int64(ns), ctypes.c_int64(nt),
            dptr(mo), dptr(x_singlet), dptr(x_triplet),
            hsoc_re, hsoc_im,
            status_message, ctypes.c_int32(1024), ctypes.byref(status),
        )
        if status.value != 0:
            detail = status_message.value.decode("utf-8", errors="replace").strip()
            raise RuntimeError(f"openqp-dftb soc_matrix failed: {detail}")
        # Hermitian: fortran-vs-C reshape only transposes; fix with explicit .T
        re = np.frombuffer(hsoc_re, dtype=np.float64).reshape((dim, dim)).T.copy()
        im = np.frombuffer(hsoc_im, dtype=np.float64).reshape((dim, dim)).T.copy()
        return re, im

    def _external_potential(self):
        """Per-atom external electrostatic potential (QM/MM embedding), or None.

        The QM/MM drivers set mol.dftb_external_potential = POTMM (Hartree/e,
        one value per QM atom incl. link-atom caps) before each QM call.
        """
        v_ext = getattr(self.mol, "dftb_external_potential", None)
        if v_ext is None:
            return None
        v_ext = np.asarray(v_ext, dtype=float).reshape(-1)
        if v_ext.size != self.natom:
            raise ValueError(
                f"dftb_external_potential must have one value per atom "
                f"({self.natom}), got {v_ext.size}"
            )
        return v_ext

    def _mo_tag_present(self, geom_key) -> bool:
        """True if MO coefficient tags are already published for this geometry.

        Used to avoid reverting a BasisOverlap phase/reorder alignment when
        excitation()/energy() re-stores tags at the same geometry.
        """
        if getattr(self.mol, "_dftb_mo_tag_geom", None) != geom_key:
            return False
        try:
            self.mol.data["OQP::VEC_MO_A"]
            return True
        except (KeyError, AttributeError):
            return False

    @staticmethod
    def _fortran_tag(array2d):
        """Store a true (rows, cols) matrix the way native OQP:: tags are laid out.

        Native tag arrays carry the FORTRAN column-major flat data under a
        numpy view with the Fortran shape, so python code that indexes them
        (align_mo, align_x, NACME) effectively sees the transpose. Emulate
        exactly that so all downstream consumers behave like the native path.
        """
        return np.ascontiguousarray(array2d.T).reshape(array2d.shape)

    def _store_wavefunction_tags(self, result: _StateResult) -> None:
        """Publish MO coefficients/energies and response vectors as OQP:: tags.

        Uses the SAME tag names and layout conventions as the native path
        (liboqp never sets them for method=dftb), so BasisOverlap, NACME, and
        the NAMD back_door carry work unchanged. Arrays carry DFTB dimensions.
        """
        if result.mo_coefficients is None:
            return
        data = self.mol.data
        # Preserve MO tags that were already published (and possibly
        # phase/reorder-aligned by BasisOverlap) for the CURRENT geometry: a
        # later excitation()/energy() call at the same geometry must not revert
        # the alignment the DFTB state-overlap path relies on. Response vectors
        # and dims are always refreshed (align_x re-aligns td_bvec_mo after).
        geom_key = np.ascontiguousarray(
            np.asarray(self.mol.get_system(), dtype=np.float64)).tobytes()
        if not self._mo_tag_present(geom_key):
            mo_tag = self._fortran_tag(result.mo_coefficients)
            data["OQP::VEC_MO_A"] = mo_tag
            data["OQP::VEC_MO_B"] = mo_tag.copy()
            data["OQP::E_MO_A"] = np.ascontiguousarray(result.mo_energies)
            data["OQP::E_MO_B"] = np.ascontiguousarray(result.mo_energies)
            self.mol._dftb_mo_tag_geom = geom_key
        if result.response_vectors is not None:
            data["OQP::td_bvec_mo"] = self._fortran_tag(result.response_vectors)
        data["OQP::dftb_wf_dims"] = np.array(
            [result.nbf, result.noca, result.nocb], dtype=float)

    def _cache_key(self, method: str, state: int, *, need_grad: bool):
        coords = np.ascontiguousarray(self.mol.get_system(), dtype=np.float64)
        atoms = tuple(int(z) for z in np.asarray(self.mol.get_atoms()).reshape(-1))
        v_ext = self._external_potential()
        v_ext_key = v_ext.tobytes() if v_ext is not None else b""
        option_key = (
            v_ext_key,
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
        # Normalize every accepted [dftb] type alias to one of the backend's
        # method names (ground variants, tddftb, sf, mrsf). The input checker
        # accepts spellings such as tda/td-dftb/sftddftb/mrsftddftb; the native
        # C API does not recognize all of them, so canonicalize here.
        method = method.lower()
        canon = {
            "dftb": "ground", "ground": "ground",
            "dftb0": "ground_noscc", "noscc": "ground_noscc",
            "ground_noscc": "ground_noscc",
            "tddftb": "tddftb", "tda": "tddftb", "td-dftb": "tddftb", "rpa": "tddftb",
            "sf": "sf", "sftddftb": "sf", "sf-tddftb": "sf",
            "mrsf": "mrsf", "mrsftddftb": "mrsf", "mrsf-tddftb": "mrsf",
        }
        return canon.get(method, method)

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


HA_TO_WAVENUMBER = 219474.6313708
FINE_STRUCTURE = 7.2973525693e-3


def dftb_soc(mol):
    """Standalone SOC driver for method=dftb (mirror of the native compute_soc).

    Runs the MRSF response for both target multiplicities from the same ROKS
    reference, builds the one-center SOC matrix, diagonalizes
    diag((E - e0) * ha2wn) + hsoc * (alpha^2/2 * ha2wn) in numpy, and fills the
    same OQP:: tags the native soc_mrsf produces (soc_eval in cm^-1 relative to
    the lowest MCH excitation; soc_hsoc_* raw, without the alpha^2/2 prefactor).
    """
    adapter = OpenQPDFTBAdapter(mol)
    data = mol.data
    dftb = mol.config['dftb']
    saved = dftb.get('target_multiplicity', 1)
    try:
        dftb['target_multiplicity'] = 1
        singlet_energies = np.array(OpenQPDFTBAdapter(mol).energy(), dtype=float)
        data['OQP::td_singlet_energies'] = np.asarray(data['OQP::td_energies']).copy()
        x_s = np.asarray(data['OQP::td_bvec_mo']).copy()
        data['OQP::td_bvec_mo_s'] = x_s.copy()
        mo_tag = np.asarray(data['OQP::VEC_MO_A']).copy()
        dims = np.asarray(data['OQP::dftb_wf_dims']).ravel()

        dftb['target_multiplicity'] = 3
        triplet_energies = np.array(OpenQPDFTBAdapter(mol).energy(), dtype=float)
        data['OQP::td_triplet_energies'] = np.asarray(data['OQP::td_energies']).copy()
        x_t = np.asarray(data['OQP::td_bvec_mo']).copy()
        data['OQP::td_bvec_mo_t'] = x_t.copy()
    finally:
        dftb['target_multiplicity'] = saved

    mol.singlet_energies = singlet_energies
    mol.triplet_energies = triplet_energies

    nbf, noca, nocb = (int(round(v)) for v in dims[:3])
    hsoc_re, hsoc_im = adapter.soc_matrix(
        mo_tag.ravel(), x_s.ravel(), x_t.ravel(), noca=noca, nocb=nocb)

    ns = len(singlet_energies) - 1
    nt = len(triplet_energies) - 1
    e_s = singlet_energies[1:]
    e_t = triplet_energies[1:]
    e0 = min(e_s[0], e_t[0])
    diag = np.concatenate([e_s - e0, np.repeat(e_t - e0, 3)]) * HA_TO_WAVENUMBER
    dfac = 0.5 * FINE_STRUCTURE ** 2 * HA_TO_WAVENUMBER
    h_total = np.diag(diag).astype(complex) + (hsoc_re + 1j * hsoc_im) * dfac
    eigenvalues, eigenvectors = np.linalg.eigh(h_total)

    fortran_tag = OpenQPDFTBAdapter._fortran_tag
    data['OQP::soc_eval'] = np.ascontiguousarray(eigenvalues.real)
    data['OQP::soc_evec_re'] = fortran_tag(np.ascontiguousarray(eigenvectors.real))
    data['OQP::soc_evec_im'] = fortran_tag(np.ascontiguousarray(eigenvectors.imag))
    data['OQP::soc_hsoc_re'] = fortran_tag(np.ascontiguousarray(hsoc_re))
    data['OQP::soc_hsoc_im'] = fortran_tag(np.ascontiguousarray(hsoc_im))
    mol.soc = eigenvalues.real
    return eigenvalues.real
