"""OpenQP runtime adapter for the optional OpenQP-DFTB backend."""

from __future__ import annotations

import ctypes
from dataclasses import dataclass
import itertools
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import textwrap
import threading
import time

import numpy as np

import oqp
from oqp.periodic_table import ELEMENTS_NAME
from oqp.utils.constants import ANGSTROM_TO_BOHR as BOHR_TO_ANGSTROM
from oqp.utils import dftb_trace
from oqp.utils.file_utils import dump_log
from oqp.utils.state_labels import (
    DFTB_CAP_STATE_SPECTRUM,
    DFTB_CAP_STRUCTURED_TRACE,
    canonical_dftb_type,
    dftb_method_name,
    public_state_label,
    resolved_dftb_type,
)


_GROUND_TYPES = {"ground", "dftb", "dftb0", "ground_noscc", "noscc"}
_NATIVE_DIAGNOSTIC_LOCK = threading.RLock()
_NATIVE_DIAGNOSTIC_ENV = (
    "OPENQP_DFTB_SCC_TRACE",
    "OPENQP_DFTB_SOLVER_TRACE",
    "OPENQP_ZVEC_TRACE",
    "OPENQP_DFTB_STATE_SPECTRUM",
)


# Warm-start seed paths whose stale files have already been cleared once in
# this process (see OpenQPDFTBAdapter._maybe_enable_optimization_solution_following).
_OPT_SOLUTION_FOLLOW_SEEDS = set()


def _call_with_native_diagnostics(
        callback, *, print_level, state_spectrum, structured_trace=True):
    """Run one native call while capturing flushed Fortran stdout.

    The C ABI predates structured diagnostics.  Its optional trace hooks write
    to the preconnected Fortran output unit, so capture file descriptor 1 for
    the duration of the call and return the text for the normal OpenQP log.
    Environment and descriptor changes are process-global; serialize them
    within this Python process and restore both even on failure.
    """
    requested_level = max(0, min(2, int(print_level)))
    # The per-iteration records the OpenQP-style tables are built from are
    # emitted at native trace level 2 only; request them whenever progress is
    # wanted at all and let the Python formatter decide how much to show.
    level = 2 if (structured_trace and requested_level > 0) else 0
    updates = {
        "OPENQP_DFTB_SCC_TRACE": str(level),
        "OPENQP_DFTB_SOLVER_TRACE": str(level),
        "OPENQP_ZVEC_TRACE": str(level),
        "OPENQP_DFTB_STATE_SPECTRUM": "1" if state_spectrum else "0",
    }

    with _NATIVE_DIAGNOSTIC_LOCK:
        previous = {name: os.environ.get(name) for name in _NATIVE_DIAGNOSTIC_ENV}
        for name, value in updates.items():
            os.environ[name] = value
        saved_stdout = None
        try:
            sys.stdout.flush()
            saved_stdout = os.dup(1)
            with tempfile.TemporaryFile(mode="w+b") as capture:
                os.dup2(capture.fileno(), 1)
                try:
                    callback()
                finally:
                    # Core trace and spectrum writers flush their Fortran
                    # output unit; restore stdout before decoding the capture.
                    os.dup2(saved_stdout, 1)
                capture.seek(0)
                return capture.read().decode("utf-8", errors="replace")
        finally:
            if saved_stdout is not None:
                os.close(saved_stdout)
            for name, value in previous.items():
                if value is None:
                    os.environ.pop(name, None)
                else:
                    os.environ[name] = value

# C ABI v1 includes embedding/exports and the first three DTCAM controls
# (c_mrsf, response_global_hybrid, onsite_exchange_scale). These later v2
# controls are the only current options its shorter argument list cannot carry.
_ABI1_UNSUPPORTED_OPTION_DEFAULTS = {
    "c_mrsf_oo": -1.0,
    "w_scale": 1.0,
    "response_w_scale": -1.0,
    "response_omega": -1.0,
    "response_cam_alpha": -1.0,
    "response_cam_beta": -1.0,
    "onsite_ss": 0.0,
    "onsite_sp": 0.0,
    "onsite_pp": 0.0,
}


def _state_gradient_argtypes(abi_version: int) -> list[object]:
    """Return the exact ctypes signature of ``openqp_dftb_state_gradient``.

    The C entry point is a fixed-arity Fortran ``bind(C)`` routine.  Leaving
    it untyped makes ctypes use its permissive foreign-call path; an argument
    layout drift can then corrupt native pointers instead of producing a
    Python exception.  The explicit signature also matters on arm64, where
    fixed and variadic foreign calls use different ABI classification rules.
    """
    i64 = ctypes.c_int64
    i32 = ctypes.c_int32
    f64 = ctypes.c_double
    p_i64 = ctypes.POINTER(i64)
    p_f64 = ctypes.POINTER(f64)

    # slots 1--37 in _state_gradient_common_arguments
    common = [
        i64, p_i64, p_f64, ctypes.c_char_p, i32, ctypes.c_char_p, i32,
        *([i64] * 16),
        *([f64] * 14),
    ]
    # c_mrsf, response_global_hybrid, onsite_exchange_scale
    controls = [f64, i64, f64]
    if abi_version != 1:
        # ABI v2/v3 extended DTCAM controls and the preset name.
        controls.extend([
            *([f64] * 9), ctypes.c_char_p, i32,
        ])
    if abi_version >= 4:
        # ABI v4 appended the reference_selector (0 default / 1 UKS / 2 ROKS).
        controls.append(i64)
    # n_ext_pot/ext_potential followed by the result/output tail.
    tail = [
        i64, p_f64,
        p_f64, p_f64, p_f64, p_f64, p_f64,
        p_i64, p_f64, p_f64, p_f64,
        p_i64, p_i64, p_i64, p_i64,
        i64, p_f64, p_f64, i64, p_f64,
        ctypes.c_char_p, i32, p_i64,
    ]
    return common + controls + tail


def _bundled_parameter_path() -> str | None:
    """Bundled OB2W0PT3 default of the installed openqp-dftb wheel, if any.

    The wheel ships the SKF set plus the official shell-resolved spinw.txt
    (required by every W kernel) next to its locator package; wheels that
    predate the bundled data simply lack the accessor.
    """
    try:
        import openqp_dftb  # noqa: PLC0415

        return openqp_dftb.default_parameter_path()
    except (ImportError, AttributeError, FileNotFoundError):
        return None



def _bundled_parameter_path() -> str | None:
    """Return the installed openqp-dftb default parameter set, if available."""

    try:
        import openqp_dftb  # noqa: PLC0415

        return openqp_dftb.default_parameter_path()
    except (ImportError, AttributeError, FileNotFoundError):
        return None


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
    """Make OpenQP-DFTB look like a normal OpenQP energy/gradient provider.

    The backend-identity knobs below are class attributes so the OpenQP-xTB
    adapter (oqp.library.openqp_xtb) can subclass this adapter and only swap
    the section name, C-symbol prefix, library basenames, and env vars. The
    tag names published by the store/response helpers and the
    mol.dftb_external_potential attribute are intentionally SHARED by both TB
    backends (see oqp.utils.tb_backends).
    """

    SECTION = "dftb"
    SYMBOL_PREFIX = "openqp_dftb"
    LIB_BASENAMES = ("libopenqp_dftb_c.dylib", "libopenqp_dftb_c.so")
    ENV_LIBRARY = "OPENQP_DFTB_LIBRARY"
    ENV_PARAMETER = "OPENQP_DFTB_PARAMETER_PATH"
    PIP_LOCATOR = "openqp_dftb"
    BACKEND_NAME = "openqp-dftb"      # pip package / diagnostics name
    DISPLAY_NAME = "OpenQP-DFTB"      # human-readable diagnostics name
    CACHE_ATTR = "_openqp_dftb_cache"
    # ABI negotiation: DFTB accepts C ABI v1/2/3 and treats an unversioned
    # library (no version symbol) as the stable v1 layout. OpenQP-xTB overrides
    # these (its first released layout carries the model block, so an
    # unversioned probe must resolve to the current layout, not legacy v1).
    _SUPPORTED_ABI = (1, 2, 3, 4)
    _UNVERSIONED_ABI_FALLBACK = 1

    def __init__(self, mol):
        self.mol = mol
        self.config = mol.config
        self.dftb = self.config.get(self.SECTION, {})
        # Open-shell ground-state reference selection ([dftb] reference=...).
        # Lets a plain `dftb`/`dftb0` ground SCC run an open-shell reference
        # (ROKS / CUKS / UKS) WITHOUT the MRSF response -- the reference energy
        # decomposition (H0/gamma/repulsive/spin) is what one compares to DFTB+.
        #   reference = rhf|rks   -> closed-shell (default)
        #               roks|rohf -> restricted-open (plain Guest-Saunders)
        #               cuks|cuhf -> constrained-UKS (default open-shell operator)
        #               uks|uhf   -> genuine unrestricted (DFTB+ udftb analogue)
        # `unpaired` (default 2) sets the number of unpaired electrons.
        # reference_selector crosses the C-API as an ABI-4 argument (0 default /
        # 1 UKS / 2 ROKS).  For an older (ABI < 4) native library the selector is
        # instead driven by the OPENQP_DFTB_UKS/CUHF env vars -- but those are set
        # *around the native call* (see _run_native) and restored afterwards, NOT
        # here: setting them in __init__ would leak this molecule's reference into
        # a later default-reference job sharing the process (adapters may be
        # constructed before any of them runs, and an ABI<4 default job never
        # re-enters this branch to clear a previous setting).
        self._reference = str(self.dftb.get("reference", "")).strip().lower()
        self._reference_selector = 0
        if self._reference:
            _open = self._reference in {
                "roks", "rohf", "cuks", "cuhf", "uks", "uhf"}
            if _open and int(self.dftb.get("reference_multiplicity", 0)) <= 1:
                self.dftb["reference_multiplicity"] = (
                    int(self.dftb.get("unpaired", 2)) + 1
                )
            if self._reference in {"uks", "uhf"}:
                self._reference_selector = 1
            elif self._reference in {"roks", "rohf"}:
                self._reference_selector = 2
        self.natom = int(mol.data["natom"])
        self.nstate = self._effective_nstate()
        # Named operator preset ([dftb] model=...). The published parameter
        # vector lives in openqp-dftb (single source of truth) and is applied
        # by the native library; the input checker forbids mixing a model
        # with individual operator keys.  'none' is the explicit opt-out of
        # the SF/MRSF dtcam default (materialized in Molecule.get_config).
        preset = str(self.dftb.get("model", "")).strip()
        if preset.lower() == "none":
            preset = ""
        self._preset_bytes = preset.encode("ascii")
        if not hasattr(mol, self.CACHE_ATTR):
            setattr(mol, self.CACHE_ATTR, {})
        # SF/MRSF-TDDFTB uses a high-spin ROKS reference: mark scf.type
        # accordingly so generic bookkeeping (get_data tag selection for the
        # NAMD back_door carry, beta-set handling) treats both spin sets.
        if self._probe_method_name(self._resolved_method()) in {"sf", "mrsf"} and \
                str(self.config.get("scf", {}).get("type", "rhf")).lower() == "rhf":
            self.config.setdefault("scf", {})["type"] = "rohf"
        self._maybe_enable_optimization_solution_following()

    def _maybe_enable_optimization_solution_following(self):
        """Steer the bistable LC-ROKS reference SCF onto a consistent branch
        during geometry optimizations and seam searches.

        The long-range-corrected restricted-open-shell reference SCF has
        multiple self-consistent solutions in near-degenerate regions (bond
        alternated interiors, near-CI geometries).  A geometry optimizer that
        solves every step from a cold guess lands on different branches at
        neighboring points, producing SCC non-convergence and spurious jumps in
        MECI / S1-min searches.  Because the optimizer runs every SCC in this
        Python process, a guess-only warm-start file -- dumped after each SCC
        and loaded before the next -- carries the converged reference orbitals
        from one geometry step to the next, following a single solution branch.
        No-op for single points; skipped if the user configured the warm-start
        env explicitly or set OPENQP_DFTB_OPT_SOLUTION_FOLLOW=0.
        """
        runtype = str(self.config.get("input", {}).get("runtype", "energy")).lower()
        if runtype not in {"optimize", "meci", "mecp", "mep", "ts"}:
            return
        if os.environ.get("OPENQP_DFTB_OPT_SOLUTION_FOLLOW", "1") == "0":
            return
        if os.environ.get("OPENQP_DFTB_WARMSTART_MO") or \
                os.environ.get("OPENQP_DFTB_DUMP_MO"):
            return  # respect an explicit user-supplied warm-start configuration
        seed = os.path.join(
            self.mol.log_path or ".",
            ".{0}_{1}_dftb_scc_seed.mo".format(self.mol.project_name, os.getpid()),
        )
        # Start each optimization from a cold first step: drop any stale seed
        # left by an earlier same-PID run the first time this path is used.
        if seed not in _OPT_SOLUTION_FOLLOW_SEEDS:
            _OPT_SOLUTION_FOLLOW_SEEDS.add(seed)
            try:
                os.remove(seed)
            except OSError:
                pass
        os.environ["OPENQP_DFTB_DUMP_MO"] = seed
        os.environ["OPENQP_DFTB_WARMSTART_MO"] = seed

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
            self._dump_excited_state_summary(base)

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
        if backend != "probe":
            return
        qmmm_flag = self.config.get("input", {}).get("qmmm_flag", False)
        is_qmmm = (qmmm_flag is True) or (
            str(qmmm_flag).strip().lower() in {"true", "1", "on", "yes"})
        if self._external_potential() is not None or is_qmmm:
            raise ValueError(
                "OpenQP-DFTB QM/MM requires the native backend: the probe "
                "executable neither carries the per-atom external potential "
                "(electrostatic embedding) nor publishes the relaxed atomic "
                "charges the QM/MM driver reads to assemble the MM coupling "
                "forces, so it would return unembedded energies/gradients or "
                "fail on missing OQP::partial_charges (mechanical embedding). "
                "Set [dftb] backend=native for QM/MM jobs."
            )

    def _resolve_model_default(self) -> None:
        """Materialize the production [dftb] model default, ABI permitting.

        The default is decided here -- after the native library is loaded --
        rather than at input parsing, because a C ABI v1 library cannot carry
        model presets at all: materializing the default unconditionally would
        abort previously working no-model ABI-v1 jobs.  On ABI v1 the default
        is simply skipped and the input keeps its explicit-keys meaning; an
        explicit user-provided model on ABI v1 still fails loudly in
        _validate_abi1_request.
        """
        if getattr(self.mol, "_dftb_model_default_resolved", False):
            return
        backend = str(self.dftb.get("backend", "native")).lower()
        if backend not in {"native", "auto"}:
            return
        lib = self._native_library()  # noqa: F841 -- caches the ABI version
        self.mol._dftb_model_default_resolved = True
        if int(getattr(self.mol, self.CACHE_ATTR)["__native_abi_version__"]) == 1:
            return
        from oqp.utils.input_checker import apply_dftb_model_default  # noqa: PLC0415
        model = apply_dftb_model_default(self.config)
        if model:
            self._preset_bytes = model.encode("ascii")

    def _run_state(self, method: str, state: int, *, need_grad: bool) -> _StateResult:
        self._resolve_model_default()
        key = self._cache_key(method, state, need_grad=need_grad)
        cache = getattr(self.mol, self.CACHE_ATTR)
        if key in cache:
            return cache[key]
        # One geometry (+ embedding potential) per cache generation: an MD
        # trajectory would otherwise grow the cache without bound.
        generation = (key[0], key[1], key[2][0])
        if cache.get("__generation__") != generation:
            library = cache.get("__native_library__")
            library_abi = cache.get("__native_abi_version__")
            library_capabilities = cache.get("__native_capabilities__")
            cache.clear()
            if library is not None:
                cache["__native_library__"] = library
            if library_abi is not None:
                cache["__native_abi_version__"] = library_abi
            if library_capabilities is not None:
                cache["__native_capabilities__"] = library_capabilities
            cache["__generation__"] = generation

        backend = str(self.dftb.get("backend", "native")).lower()
        self._reject_probe_with_embedding(backend)
        if backend in {"native", "auto"}:
            result = self._run_native_with_scc_recovery(
                method, state, need_grad=need_grad)
        elif backend == "probe":
            result = self._run_probe(method, state)
        else:
            raise ValueError(f"Unknown OpenQP-DFTB backend: {backend}")

        cache[key] = result
        return result

    @staticmethod
    def _is_scc_nonconvergence(error: RuntimeError) -> bool:
        """Return whether a native failure is specifically an SCC failure.

        Parameter, ABI, response, and gradient failures must not be hidden by
        changing an unrelated electronic optimizer.  The native ground-state
        drivers use the phrases below for both closed- and open-shell SCC
        exhaustion.
        """
        message = str(error).lower()
        return (
            "scc" in message
            and (
                "did not converge" in message
                or "not converged" in message
                or "nonconver" in message
            )
        )

    @staticmethod
    def _scc_recovery_ladder(primary: str) -> tuple[str, ...]:
        """Minimal per-geometry SCC ladder, ending in charge/orbital TRAH."""
        primary = str(primary).lower()
        if primary in {"trust", "trah"}:
            return (primary,)
        if primary in {"broyden", "diis"}:
            return (primary, "trust")
        if primary in {"anderson", "pulay"}:
            return (primary, "broyden", "trust")
        if primary == "linear":
            return (primary, "anderson", "broyden", "trust")
        # The native auto mixer starts from Anderson/Pulay.  Make the two
        # genuinely different fallbacks explicit after it is exhausted.
        return (primary, "broyden", "trust")

    def _run_native_with_scc_recovery(
            self, method: str, state: int, *, need_grad: bool) -> _StateResult:
        """Retry SCC failures at this geometry, then restore the primary mixer.

        Geometry optimizers construct a fresh adapter at each evaluation, but
        ``dftb`` is the molecule's shared configuration dictionary.  Restoring
        it in ``finally`` is therefore essential: a geometry that needed
        Broyden or TRAH must not force the next geometry to start there.
        """
        primary = str(self.dftb.get("scc_mixer", "auto")).lower()
        ladder = self._scc_recovery_ladder(primary)
        try:
            for attempt, mixer in enumerate(ladder):
                self.dftb["scc_mixer"] = mixer
                try:
                    return self._run_native(method, state, need_grad=need_grad)
                except RuntimeError as error:
                    if not self._is_scc_nonconvergence(error) or attempt + 1 == len(ladder):
                        raise
                    next_mixer = ladder[attempt + 1]
                    dump_log(
                        self.mol,
                        title=(
                            "PyOQP: DFTB SCC recovery at the current geometry\n"
                            f"   {mixer} did not converge; retrying with {next_mixer}.\n"
                            f"   The next geometry will restart with {primary}."
                        ),
                    )
        finally:
            self.dftb["scc_mixer"] = primary

        raise AssertionError("DFTB SCC recovery ladder exhausted unexpectedly")

    def _run_native(self, method: str, state: int, *, need_grad: bool) -> _StateResult:
        """Call libopenqp_dftb_c (the standalone openqp-dftb shared C API) in-process.

        openqp-dftb stays a fully separate library: PyOQP hands over geometry and
        options as plain C scalars/arrays and receives energies and the analytic
        state gradient back. No OpenQP build coupling, no Fortran module seam.
        """
        lib = self._native_library()
        abi_version = int(getattr(self.mol, self.CACHE_ATTR)["__native_abi_version__"])
        capabilities = int(
            getattr(self.mol, self.CACHE_ATTR).get("__native_capabilities__", 0)
        )
        parameter_path = self._parameter_path()
        self._log_settings_once(
            method,
            backend="native",
            parameter_path=parameter_path,
            library_path=str(getattr(lib, "_name", "")),
            abi_version=abi_version,
            capabilities=capabilities,
        )
        natom = self.natom
        atoms = np.ascontiguousarray(
            np.asarray(self.mol.get_atoms(), dtype=np.int64).reshape(-1)
        )
        coords_bohr = np.ascontiguousarray(
            np.asarray(self.mol.get_system(), dtype=np.float64).reshape(-1)
        )
        parameter = parameter_path.encode("utf-8")
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

        call_kwargs = {
            "natom": natom,
            "atoms": atoms,
            "coords_bohr": coords_bohr,
            "parameter": parameter,
            "method_name": method_name,
            "state": state,
            "need_grad": need_grad,
            "n_ext": n_ext,
            "ext_potential": ext_potential,
            "reference_energy": reference_energy,
            "state_energy": state_energy,
            "excitation_energy": excitation_energy,
            "spin_square": spin_square,
            "gradient": gradient,
            "n_roots_out": n_roots_out,
            "all_state_energies": all_state_energies,
            "all_spin_squares": all_spin_squares,
            "relaxed_charges": relaxed_charges,
            "nbf_out": nbf_out,
            "noca_out": noca_out,
            "nocb_out": nocb_out,
            "vec_dim_out": vec_dim_out,
            "mo_capacity": mo_capacity,
            "mo_energies": mo_energies,
            "mo_coefficients": mo_coefficients,
            "vec_capacity": vec_capacity,
            "response_vectors": response_vectors,
            "status_message": status_message,
            "status": status,
        }

        def native_call():
            fn = getattr(lib, f"{self.SYMBOL_PREFIX}_state_gradient")

            def _dispatch():
                if abi_version == 1:
                    self._call_state_gradient_abi1(fn, **call_kwargs)
                else:
                    self._call_state_gradient_current(fn, **call_kwargs)

            # ABI >= 4 carries reference_selector as an explicit C argument, so no
            # environment is involved.  Only an older library needs the env
            # fallback -- set it immediately around the call and restore the prior
            # values afterwards so a UKS/ROKS job cannot leak its reference into a
            # subsequent default-reference job that shares this process.
            if abi_version >= 4 or not self._reference:
                _dispatch()
                return
            with _NATIVE_DIAGNOSTIC_LOCK:
                saved = {
                    key: os.environ.get(key)
                    for key in ("OPENQP_DFTB_UKS", "OPENQP_DFTB_CUHF")
                }
                os.environ.pop("OPENQP_DFTB_UKS", None)
                os.environ.pop("OPENQP_DFTB_CUHF", None)
                if self._reference_selector == 1:
                    os.environ["OPENQP_DFTB_UKS"] = "1"
                elif self._reference_selector == 2:
                    os.environ["OPENQP_DFTB_CUHF"] = "0"
                try:
                    _dispatch()
                finally:
                    for key, value in saved.items():
                        if value is None:
                            os.environ.pop(key, None)
                        else:
                            os.environ[key] = value

        if abi_version == 1:
            self._validate_abi1_request()

        state_spectrum = (
            method not in _GROUND_TYPES
            and not need_grad
            and bool(self.dftb.get("state_to_state_spectrum", True))
            and bool(capabilities & DFTB_CAP_STATE_SPECTRUM)
        )
        step_wall_start = time.perf_counter()
        step_cpu_start = time.process_time()
        native_trace = _call_with_native_diagnostics(
            native_call,
            print_level=int(self.dftb.get("print_level", 1)),
            state_spectrum=state_spectrum,
            structured_trace=bool(capabilities & DFTB_CAP_STRUCTURED_TRACE),
        )
        step_times = (time.process_time() - step_cpu_start,
                      time.perf_counter() - step_wall_start)
        charge_info = None
        mo_info = None
        if status.value == 0:
            charge_info = self._stash_mulliken_charges(
                method, state, need_grad,
                np.frombuffer(relaxed_charges, dtype=np.float64)[:natom],
                atoms,
            )
            if not need_grad:
                mo_info = self._stash_orbital_info(
                    method, int(nbf_out.value), int(noca_out.value),
                    int(nocb_out.value), mo_energies, mo_coefficients)
        if native_trace.strip():
            self._log_native_progress(method, state, need_grad, native_trace,
                                      charge_info=charge_info,
                                      mo_info=mo_info,
                                      step_times=step_times)
        self._raise_native_status(
            method=method,
            state=state,
            status=status,
            status_message=status_message,
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
            stdout=native_trace,
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

    def _public_target_label(self, state):
        """Public label for a response root (S1, T0, ...; 'state N' fallback)."""
        label = public_state_label(self.config, state)
        if label.startswith("state "):
            label = "state %d" % int(state)
        return label

    def _stash_mulliken_charges(self, method, state, need_grad, charges, atoms):
        """Record net Mulliken charges + point-charge dipole on the molecule.

        On energy calls the native library returns the converged SCC reference
        charges; on excited-state gradient calls it returns the Z-vector
        relaxed charges of the target state.  Both are stored for the results
        JSON; the returned dict feeds the log block.
        """
        charges = [float(value) for value in charges]
        coords = np.asarray(self.mol.get_system(),
                            dtype=np.float64).reshape(-1, 3)
        dipole = dftb_trace.point_charge_dipole(charges, coords)
        symbols = [ELEMENTS_NAME[int(z)].strip() for z in atoms]
        if need_grad and int(state) > 0:
            title = ("Relaxed Mulliken charges (%s)"
                     % self._public_target_label(state))
            self.mol.dftb_relaxed_mulliken = {
                "state": int(state),
                "charges": charges,
                "dipole_au": dipole,
            }
        else:
            if canonical_dftb_type(method) in {"sf", "mrsf"}:
                reference = "high-spin ROKS reference"
            else:
                reference = "SCC reference"
            title = "Mulliken atomic charges (%s)" % reference
            self.mol.dftb_mulliken = {
                "reference": reference,
                "charges": charges,
                "dipole_au": dipole,
            }
        return {
            "title": title,
            "symbols": symbols,
            "charges": charges,
            "dipole": dipole,
        }

    def _timing_footer(self, step_times):
        """Native-OpenQP style step + cumulative CPU/wall timing lines."""
        if step_times is None:
            return ""
        step_cpu, step_wall = step_times
        start = getattr(self.mol, "start_time", None)
        total_wall = (time.time() - start) if start else step_wall
        return dftb_trace.format_step_timing(
            step_cpu, step_wall, time.process_time(), total_wall)

    def _dump_effective_options(self, parsed):
        """Report the preset-resolved native options once per distinct set."""
        options = parsed.get("resolved_options")
        if not options:
            return
        self.mol.dftb_resolved_options = options
        if getattr(self.mol, "_dftb_resolved_options_logged", None) == options:
            return
        self.mol._dftb_resolved_options_logged = options
        text = dftb_trace.format_resolved_options(options)
        if text:
            dump_log(self.mol, title="PyOQP: DFTB effective options",
                     section="text", info={"text": text})

    def _log_native_progress(self, method, state, need_grad, native_trace,
                             charge_info=None, mo_info=None, step_times=None):
        """Render the captured native trace as OpenQP-style iteration tables.

        openqp-dftb cannot append to the OpenQP log itself (the log unit
        belongs to liboqp), so its structured stdout records are parsed here
        and written in the same style as the native SCF/Davidson/Z-vector
        sections of an OpenQP log.  The state-to-state spectrum block is NOT
        printed here: it is stored and reported after the excited-state
        summary table.
        """
        print_level = max(0, int(self.dftb.get("print_level", 1)))
        if print_level == 0:
            return
        verbose = print_level >= 2
        parsed = dftb_trace.parse_native_trace(native_trace)
        if parsed.get("energy_components"):
            self.mol.dftb_energy_components = parsed["energy_components"]
        if not any(parsed[key] for key in (
                "scc_passes", "davidson", "zvector", "zvector_dense")):
            # Nothing structured (old library or probe-style text): keep the
            # raw capture so no diagnostic is silently dropped.
            dump_log(
                self.mol,
                title="PyOQP: OpenQP-DFTB native progress",
                section="dftb_runtime",
                info={
                    "method": canonical_dftb_type(method),
                    "state": int(state),
                    "gradient": bool(need_grad),
                    "text": native_trace,
                },
            )
            return

        extra = dftb_trace.format_other_lines(parsed) if verbose else ""
        charge_block = ""
        if charge_info is not None:
            charge_block = dftb_trace.format_charge_block(
                charge_info["symbols"], charge_info["charges"],
                charge_info["dipole"], charge_info["title"])
        if not need_grad:
            self._dump_effective_options(parsed)
            sections = []
            scc_blocks = []
            scc_text = dftb_trace.format_scc_block(parsed, verbose=verbose)
            if scc_text:
                scc_blocks.append(scc_text)
            components_text = dftb_trace.format_energy_components(
                parsed.get("energy_components"),
                reference_label="converged SCC reference")
            if components_text:
                scc_blocks.append(components_text)
            if charge_block:
                scc_blocks.append(charge_block)
            if scc_blocks:
                sections.append(("PyOQP: DFTB SCC steps", scc_blocks))
            if mo_info is not None:
                mo_text = dftb_trace.format_mo_table(mo_info)
                if mo_text:
                    sections.append(("PyOQP: DFTB molecular orbitals",
                                     [mo_text]))
            response_blocks = [
                dftb_trace.format_davidson_block(
                    solve, dftb_method_name(self.config))
                for solve in parsed["davidson"]
            ]
            if response_blocks:
                sections.append(("PyOQP: DFTB response steps",
                                 response_blocks))
            if extra:
                sections.append(("PyOQP: DFTB native diagnostics", [extra]))
            timing = self._timing_footer(step_times)
            if timing and sections:
                sections[-1][1].append(timing)
            for title, blocks in sections:
                dump_log(self.mol, title=title, section="text",
                         info={"text": "\n\n".join(blocks)})
            return

        # Gradient call: the reference/response are merely re-solved for the
        # requested state, so summarize them in one line each and give the
        # Z-vector solve the full iteration table.
        blocks = []
        if verbose:
            blocks.append("   Native call: method=%s  state=%s  gradient=yes"
                          % (canonical_dftb_type(method), int(state)))
            scc_text = dftb_trace.format_scc_block(parsed, verbose=True)
            if scc_text:
                blocks.append(scc_text)
            blocks.extend(
                dftb_trace.format_davidson_block(
                    solve, dftb_method_name(self.config))
                for solve in parsed["davidson"])
        else:
            summary_lines = []
            scc_lines = dftb_trace.format_scc_summary_lines(parsed)
            if scc_lines:
                summary_lines.extend(scc_lines.splitlines())
            summary_lines.extend(dftb_trace.collapse_repeats(
                [dftb_trace.format_davidson_summary_line(solve)
                 for solve in parsed["davidson"]]))
            if summary_lines:
                blocks.append("\n".join(summary_lines))
        blocks.extend(dftb_trace.format_zvector_block(solve)
                      for solve in parsed["zvector"])
        blocks.extend(dftb_trace.format_zvector_dense_block(solve)
                      for solve in parsed["zvector_dense"])
        if charge_block:
            blocks.append(charge_block)
        if extra:
            blocks.append(extra)
        timing = self._timing_footer(step_times)
        if timing:
            blocks.append(timing)
        if not blocks:
            return
        if state > 0:
            title = ("PyOQP: DFTB gradient steps (%s)"
                     % self._public_target_label(state))
        else:
            title = "PyOQP: DFTB gradient steps (ground state)"
        dump_log(self.mol, title=title, section="text",
                 info={"text": "\n\n".join(block for block in blocks if block)})

    # Configurations below this |coefficient| stay out of the per-state
    # excitation analysis (the native MRSF-TDDFT log uses the same cutoff).
    _CONFIGURATION_COEFF_CUTOFF = 0.05

    def _dftb_shells(self, nbf):
        """Reconstruct the minimal-SK shell list [(atom, l, pure)], or None.

        The C API does not export the per-AO map, but the minimal valence
        basis is deterministic per element (1 = s, 4 = s+p, 9 = s+p+d AOs)
        and shells are stored per atom in ascending l.  The assignment is
        accepted only when exactly ONE per-element size combination matches
        the exported nbf; d elements are declined because the SK d-component
        order has not been validated against the labeler's convention.
        """
        atoms = [int(z) for z in np.asarray(self.mol.get_atoms()).reshape(-1)]
        elements = sorted(set(atoms))
        counts = {z: atoms.count(z) for z in elements}
        matches = [
            dict(zip(elements, sizes))
            for sizes in itertools.product((1, 4, 9), repeat=len(elements))
            if sum(counts[z] * size for z, size in zip(elements, sizes)) == nbf
        ]
        if len(matches) != 1:
            return None
        assignment = matches[0]
        if any(size == 9 for size in assignment.values()):
            return None
        shells = []
        # The labeler indexes detection permutations with the shell's atom
        # index directly, so atoms are 0-based here.
        for atom_index, z in enumerate(atoms):
            for l in range({1: 1, 4: 2}[assignment[z]]):
                shells.append((atom_index, l, True))
        return shells

    def _mo_symmetry_labels(self, coefficients, nbf):
        """Abelian irrep label per MO when symmetry labeling is active.

        Reuses the geometry-based point-group detection and the metadata-only
        irrep assigner of the native path.  The DFTB s/p AO component order
        (y, z, x) matches the labeler's CCA m-ascending spherical convention,
        and the AO overlap follows from the complete orthonormal MO set as
        S = (C C^T)^-1.  Any failure returns None (labeling is optional).
        """
        meta = getattr(self.mol, "symmetry_metadata", None) or {}
        if not meta or meta.get("status", "disabled") == "disabled":
            return None
        if not meta.get("label_mo", True):
            return None
        detection = meta.get("detection")
        if not detection:
            return None
        shells = self._dftb_shells(nbf)
        if shells is None:
            return None
        try:
            from oqp.library.symmetry import assign_mo_irreps  # noqa: PLC0415

            overlap = np.linalg.inv(coefficients @ coefficients.T)
            tolerance = max(float(meta.get("tolerance", 1.0e-5)), 1.0e-4)
            result = assign_mo_irreps(
                coefficients, overlap, shells,
                detection["operations"], detection["character_table"],
                tolerance=tolerance,
                matrix_key="matrix_input_frame",
            )
            meta.setdefault("mo_labels", {})["dftb"] = result
            labels = list(result.get("labels", []))
            return labels or None
        except Exception as exc:
            # Metadata only: record why labeling was skipped, never fail.
            meta.setdefault("mo_labels", {})["dftb"] = {
                "status": "error", "error": str(exc),
            }
            return None

    def _stash_orbital_info(self, method, nbf, noca, nocb,
                            mo_energies, mo_coefficients):
        """Record MO energies/occupations (+irrep labels) for log and JSON."""
        if nbf <= 0:
            return None
        energies = np.frombuffer(mo_energies, dtype=np.float64)[:nbf].copy()
        coefficients = np.frombuffer(
            mo_coefficients, dtype=np.float64)[: nbf * nbf].reshape(
            (nbf, nbf), order="F").copy()
        occupations = [2 if index <= nocb else (1 if index <= noca else 0)
                       for index in range(1, nbf + 1)]
        labels = self._mo_symmetry_labels(coefficients, nbf)
        if canonical_dftb_type(method) in {"sf", "mrsf"}:
            reference = "high-spin ROKS reference"
        else:
            reference = "SCC reference"
        mo_info = {
            "reference": reference,
            "energies": [float(value) for value in energies],
            "occupations": occupations,
            "labels": labels,
            "noca": int(noca),
            "nocb": int(nocb),
        }
        self.mol.dftb_mo = mo_info
        return mo_info

    def _state_configurations(self, result, labels):
        """Dominant SF/MRSF amplitudes per state, keyed by public label."""
        vectors = result.response_vectors
        noca, nocb, nbf = int(result.noca), int(result.nocb), int(result.nbf)
        if vectors is None or noca <= 0 or nbf <= nocb:
            return {}
        if vectors.shape[0] != noca * (nbf - nocb):
            return {}
        configurations = {}
        for position, label in enumerate(labels):
            if position >= vectors.shape[1]:
                break
            column = vectors[:, position]
            entries = []
            for raw_index in np.nonzero(
                    np.abs(column) >= self._CONFIGURATION_COEFF_CUTOFF)[0]:
                index = int(raw_index) + 1
                occ, vir = dftb_trace.spin_flip_pair(index, noca, nocb)
                entries.append({
                    "index": index,
                    "coeff": float(column[raw_index]),
                    "occ": occ,
                    "vir": vir,
                })
            configurations[label] = entries
        return configurations

    def _state_irreps(self, configurations, mo_labels, noca, nocb):
        """State irrep from the dominant configuration (abelian groups).

        The SF/MRSF CSF is the reference determinant with one alpha-occupied
        electron moved to a beta-virtual orbital, so its irrep is the product
        of the reference irrep (the singly occupied orbitals) with the hole
        and particle irreps.  Metadata only; any failure returns {}.
        """
        meta = getattr(self.mol, "symmetry_metadata", None) or {}
        detection = meta.get("detection")
        if not detection:
            return {}
        try:
            from oqp.library.symmetry import product_irrep  # noqa: PLC0415

            table = detection["character_table"]
            somos = [mo_labels[index - 1]
                     for index in range(int(nocb) + 1, int(noca) + 1)]
            reference = product_irrep(somos, table) if somos else None
            irreps = {}
            for label, entries in configurations.items():
                if not entries:
                    continue
                dominant = max(entries, key=lambda entry: abs(entry["coeff"]))
                factors = [mo_labels[dominant["occ"] - 1],
                           mo_labels[dominant["vir"] - 1]]
                if reference:
                    factors.append(reference)
                irreps[label] = product_irrep(factors, table)
            return irreps
        except Exception:
            return {}

    def _excited_state_summary(self, result):
        """Build the excited-state summary payload from one response solve."""
        method = self._resolved_method()
        if method in _GROUND_TYPES or result.all_state_energies is None:
            return None
        canonical = canonical_dftb_type(method)
        energies = [float(value) for value in result.all_state_energies]
        n_roots = len(energies)
        if n_roots == 0:
            return None
        spins = ([float(value) for value in result.all_spin_squares]
                 if result.all_spin_squares is not None else [])
        parsed = dftb_trace.parse_native_trace(result.stdout or "")
        spectrum = parsed.get("spectrum")

        if canonical == "mrsf":
            labels = [public_state_label(self.config, root)
                      for root in range(1, n_roots + 1)]
            prefix = labels[0][0] if labels and labels[0][1:].isdigit() else "S"
            reference_label = "Reference: high-spin ROKS (internal)"
            reference_energy = result.reference_energy
        elif canonical == "sf":
            labels = ["root %d" % root for root in range(1, n_roots + 1)]
            prefix = None
            reference_label = "Reference: high-spin ROKS (internal)"
            reference_energy = result.reference_energy
        else:
            # Closed-shell TDDFTB: the reference IS S0; roots are S1..Sn.
            labels = ["S0"] + ["S%d" % root for root in range(1, n_roots + 1)]
            energies = [float(result.reference_energy)] + energies
            spins = [0.0] + spins
            prefix = "S"
            reference_label = None
            reference_energy = None

        oscillator = {}
        transitions = []
        approximation = ""
        if spectrum:
            approximation = spectrum.get("approximation", "")
            oscillator = {
                label: row for label, row in dftb_trace.spectrum_from_ground(
                    spectrum, prefix).items()
            }
            transitions = [{
                "initial": dftb_trace.map_spectrum_label(row["initial"], prefix),
                "final": dftb_trace.map_spectrum_label(row["final"], prefix),
                "de_ev": row["de_ev"],
                "dipole": row["dipole"],
                "dipole_xyz": row.get("dipole_xyz"),
                "oscillator": row["oscillator"],
            } for row in spectrum.get("rows", [])]

        configurations = {}
        if canonical in {"mrsf", "sf"}:
            configurations = self._state_configurations(result, labels)
        state_irreps = {}
        mo_labels = (getattr(self.mol, "dftb_mo", None) or {}).get("labels")
        if configurations and mo_labels:
            state_irreps = self._state_irreps(
                configurations, mo_labels, result.noca, result.nocb)

        return {
            "method": canonical,
            "state_labels": labels,
            "state_energies": energies,
            "spin_squares": spins,
            "oscillator": oscillator,
            "transitions": transitions,
            "configurations": configurations,
            "state_irreps": state_irreps,
            "approximation": approximation,
            "reference_label": reference_label,
            "reference_energy": reference_energy,
        }

    def _dump_excited_state_summary(self, base):
        """Report state energies + oscillator data once per response solve.

        This runs at the END of the energy stage, so the summary table (and
        the state-to-state spectrum after it) appear in the log BEFORE any
        'Entering Gradient Calculation' section, mirroring the native OpenQP
        TDDFT summary placement.
        """
        summary = self._excited_state_summary(base)
        if summary is None:
            return
        self.mol.dftb_excited_states = summary
        key = self._cache_key(self._resolved_method(), base.state,
                              need_grad=False)
        if getattr(self.mol, "_dftb_summary_logged_key", None) == key:
            return
        self.mol._dftb_summary_logged_key = key
        if max(0, int(self.dftb.get("print_level", 1))) == 0:
            return
        dump_log(
            self.mol,
            title="PyOQP: %s Excited States" % dftb_method_name(self.config),
            section="text",
            info={"text": dftb_trace.format_excited_state_summary(summary)},
        )

    def _state_gradient_common_arguments(self, values):
        """Arguments shared by every state-gradient ABI (slots 1--37)."""
        return [
            ctypes.c_int64(values["natom"]),
            values["atoms"].ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
            values["coords_bohr"].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            values["parameter"],
            ctypes.c_int32(len(values["parameter"])),
            values["method_name"],
            ctypes.c_int32(len(values["method_name"])),
            ctypes.c_int64(int(values["state"])),
            ctypes.c_int64(int(max(1, self.nstate))),
            ctypes.c_int64(int(bool(values["need_grad"]))),
            ctypes.c_int64(int(self.config.get("input", {}).get("charge", 0))),
            ctypes.c_int64(self._scc_mixer_code()),
            ctypes.c_int64(int(self.dftb.get("scc_history", 12))),
            ctypes.c_int64(int(self.dftb.get("max_scc_iterations", 1200))),
            ctypes.c_int64(int(self.dftb.get("response_max_iterations", 50))),
            ctypes.c_int64(int(self.dftb.get("response_max_subspace", 100))),
            ctypes.c_int64(self._response_solver_code()),
            ctypes.c_int64(
                self._reference_multiplicity(
                    values["method_name"].decode("ascii")
                )
            ),
            ctypes.c_int64(int(self.dftb.get("target_multiplicity", 1))),
            ctypes.c_int64(int(bool(self.dftb.get("spin_complete", True)))),
            ctypes.c_int64(int(bool(self.dftb.get("lc_ground_state", False)))),
            ctypes.c_int64(self._lc_gamma_code()),
            ctypes.c_int64(int(bool(self.dftb.get("zvector", True)))),
            # C ABI model block: empty for DFTB (no-op splice), the five
            # GFN1 model scalars for OpenQP-xTB. Sits between the zvector
            # flag and scc_tolerance, matching openqp_xtb_state_gradient.
            *self._model_args(),
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
        ]

    @staticmethod
    def _state_gradient_output_arguments(values):
        """Embedding and result tail shared by ABI v1 and current ABIs."""
        return [
            ctypes.c_int64(values["n_ext"]),
            values["ext_potential"],
            ctypes.byref(values["reference_energy"]),
            ctypes.byref(values["state_energy"]),
            ctypes.byref(values["excitation_energy"]),
            ctypes.byref(values["spin_square"]),
            values["gradient"],
            ctypes.byref(values["n_roots_out"]),
            values["all_state_energies"],
            values["all_spin_squares"],
            values["relaxed_charges"],
            ctypes.byref(values["nbf_out"]),
            ctypes.byref(values["noca_out"]),
            ctypes.byref(values["nocb_out"]),
            ctypes.byref(values["vec_dim_out"]),
            ctypes.c_int64(values["mo_capacity"]),
            values["mo_energies"],
            values["mo_coefficients"],
            ctypes.c_int64(values["vec_capacity"]),
            values["response_vectors"],
            values["status_message"],
            ctypes.c_int32(1024),
            ctypes.byref(values["status"]),
        ]

    def _call_state_gradient_abi1(self, function, **values) -> None:
        """Call the stable 63-argument openqp-dftb C ABI v1.

        ABI v1 includes the first DTCAM controls, QM/MM embedding, all-root
        results, relaxed charges, and wavefunction exports. ABI v2 inserted
        later DTCAM/preset arguments immediately before the shared result tail,
        so dispatching only on an accepted version number would shift pointers.
        """
        function(
            *self._state_gradient_common_arguments(values),
            ctypes.c_double(float(self.dftb.get("c_mrsf", -1.0))),
            ctypes.c_int64(
                int(bool(self.dftb.get("response_global_hybrid", False)))
            ),
            ctypes.c_double(
                float(self.dftb.get("onsite_exchange_scale", 0.0))
            ),
            *self._state_gradient_output_arguments(values),
        )

    def _call_state_gradient_current(self, function, **values) -> None:
        """Call the ABI-v2/v3 state-gradient layout."""
        function(
            *self._state_gradient_common_arguments(values),
            ctypes.c_double(float(self.dftb.get("c_mrsf", -1.0))),
            ctypes.c_int64(
                int(bool(self.dftb.get("response_global_hybrid", False)))
            ),
            ctypes.c_double(
                float(self.dftb.get("onsite_exchange_scale", 0.0))
            ),
            ctypes.c_double(float(self.dftb.get("w_scale", 1.0))),
            ctypes.c_double(float(self.dftb.get("response_w_scale", -1.0))),
            ctypes.c_double(float(self.dftb.get("response_omega", -1.0))),
            ctypes.c_double(float(self.dftb.get("response_cam_alpha", -1.0))),
            ctypes.c_double(float(self.dftb.get("response_cam_beta", -1.0))),
            ctypes.c_double(float(self.dftb.get("c_mrsf_oo", -1.0))),
            ctypes.c_double(float(self.dftb.get("onsite_ss", 0.0))),
            ctypes.c_double(float(self.dftb.get("onsite_sp", 0.0))),
            ctypes.c_double(float(self.dftb.get("onsite_pp", 0.0))),
            self._preset_bytes,
            ctypes.c_int32(len(self._preset_bytes)),
            *([ctypes.c_int64(getattr(self, "_reference_selector", 0))]
              if int(getattr(self.mol, self.CACHE_ATTR)[
                  "__native_abi_version__"]) >= 4 else []),
            *self._state_gradient_output_arguments(values),
        )

    def _validate_abi1_request(self) -> None:
        """Reject features that cannot be represented by the ABI-v1 call."""
        limitations = []

        model = str(self.dftb.get("model", "")).strip()
        if model:
            limitations.append(f"model={model!r}")

        for key, default in _ABI1_UNSUPPORTED_OPTION_DEFAULTS.items():
            value = self.dftb.get(key, default)
            if isinstance(default, bool):
                differs = bool(value) is not default
            else:
                differs = float(value) != default
            if differs:
                limitations.append(f"{key}={value!r}")

        if limitations:
            detail = ", ".join(limitations)
            raise RuntimeError(
                "The loaded openqp-dftb library provides C ABI v1, which cannot "
                f"represent: {detail}. Install openqp-dftb >= 0.2.0 (C ABI v3) "
                "or reset the unsupported option(s) to their ABI-v1 defaults."
            )

    @staticmethod
    def _raise_native_status(*, method, state, status, status_message) -> None:
        if status.value == 0:
            return
        detail = status_message.value.decode("utf-8", errors="replace").strip()
        suffix = f": {detail}" if detail else ""
        raise RuntimeError(
            f"openqp-dftb native library call failed for method={method}, "
            f"state={state}, status={status.value}{suffix}"
        )

    def _native_library(self):
        cache = getattr(self.mol, self.CACHE_ATTR)
        lib = cache.get("__native_library__")
        if lib is not None:
            if "__native_capabilities__" not in cache:
                cache["__native_capabilities__"] = \
                    self._native_capabilities(lib)
            return lib
        path = self._native_library_path()
        try:
            lib = ctypes.CDLL(str(path))
        except OSError as exc:
            raise RuntimeError(f"Could not load the openqp-dftb library {path}: {exc}") from exc
        if not hasattr(lib, f"{self.SYMBOL_PREFIX}_state_gradient"):
            raise RuntimeError(
                f"{path} does not export {self.SYMBOL_PREFIX}_state_gradient; "
                "rebuild openqp-dftb with OPENQP_DFTB_BUILD_SHARED=ON."
            )
        # Detect the ABI before the first call and install its exact fixed
        # ctypes signature. The released v1 predates the version symbol.
        abi_probe = getattr(lib, self._abi_version_symbol(), None)
        if abi_probe is not None:
            abi_probe.restype = ctypes.c_int64
            abi_version = int(abi_probe())
        else:
            # Stable ABI v1 shipped alongside these two C exports. Earlier
            # unversioned development snapshots used other argument layouts
            # and cannot be distinguished by inspecting state_gradient itself.
            abi1_markers = (
                f"{self.SYMBOL_PREFIX}_states_overlap",
                f"{self.SYMBOL_PREFIX}_soc_matrix",
            )
            missing_markers = [
                symbol for symbol in abi1_markers if not hasattr(lib, symbol)
            ]
            if missing_markers:
                missing = ", ".join(missing_markers)
                raise RuntimeError(
                    f"{path} has no openqp-dftb C ABI version symbol and lacks "
                    f"the stable ABI-v1 marker export(s): {missing}. Its "
                    "state-gradient argument layout cannot be identified "
                    "safely; install openqp-dftb >= 0.2.0."
                )
            abi_version = self._UNVERSIONED_ABI_FALLBACK
        if abi_version not in self._SUPPORTED_ABI:
            raise RuntimeError(
                f"{path} exports {self.BACKEND_NAME} C ABI version {abi_version}, "
                f"but this OpenQP adapter supports only versions {self._SUPPORTED_ABI}."
            )
        self._install_state_gradient_argtypes(lib, abi_version)
        cache["__native_library__"] = lib
        cache["__native_abi_version__"] = abi_version
        cache["__native_capabilities__"] = self._native_capabilities(lib)
        return lib

    def _install_state_gradient_argtypes(self, lib, abi_version) -> None:
        """Pin the exact ctypes signature of the state-gradient entry point.

        DFTB uses the fixed per-ABI argtypes table. OpenQP-xTB overrides this:
        its layout inserts the model block, so it leaves the symbol untyped
        (the adapter always passes fully-typed ctypes objects) rather than
        install the DFTB signature.
        """
        fn = getattr(lib, f"{self.SYMBOL_PREFIX}_state_gradient")
        fn.argtypes = _state_gradient_argtypes(abi_version)
        fn.restype = None

    def _native_capabilities(self, lib) -> int:
        """Return additive native capabilities; a missing symbol means none."""
        probe = getattr(lib, f"{self.SYMBOL_PREFIX}_capi_capabilities", None)
        if probe is None:
            return 0
        probe.restype = ctypes.c_int64
        capabilities = int(probe())
        if capabilities < 0:
            raise RuntimeError(
                "openqp-dftb returned a negative C capability mask: "
                f"{capabilities}"
            )
        return capabilities

    def _native_library_path(self) -> Path:
        raw = self.dftb.get("library_path") or os.environ.get(self.ENV_LIBRARY)
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
        for name in self.LIB_BASENAMES:
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
        self._log_settings_once(
            method,
            backend="probe",
            parameter_path=parameter_path,
            executable=executable,
        )
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

    def _log_settings_once(
        self,
        method: str,
        *,
        backend: str,
        parameter_path: str,
        library_path: str = "",
        executable: str = "",
        abi_version: int | None = None,
        capabilities: int | None = None,
    ) -> None:
        """Write one settings block per distinct DFTB request, not per state."""
        signature = (
            canonical_dftb_type(method),
            backend,
            parameter_path,
            library_path,
            executable,
            abi_version,
            capabilities,
            tuple(sorted((str(key), repr(value)) for key, value in self.dftb.items())),
            tuple(sorted(
                (str(key), repr(value))
                for key, value in self.config.get("tdhf", {}).items()
            )),
        )
        logged = getattr(self.mol, "_openqp_dftb_logged_settings", None)
        if logged is None:
            logged = set()
            self.mol._openqp_dftb_logged_settings = logged
        if signature in logged:
            return
        logged.add(signature)
        dump_log(
            self.mol,
            title="PyOQP: OpenQP-DFTB settings",
            section="dftb",
            info={
                "backend": backend,
                "parameter_path": parameter_path,
                "library_path": library_path,
                "executable": executable,
                "abi_version": abi_version,
                "capabilities": capabilities,
            },
        )

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
        method = resolved_dftb_type(self.config)
        if method not in {"ground", "ground_noscc", "tddftb", "sf", "mrsf"}:
            raise ValueError(f"Unknown OpenQP-DFTB calculation type: {method!r}.")
        return method

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
        getattr(lib, f"{self.SYMBOL_PREFIX}_states_overlap")(
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
        (dim, dim) numpy arrays in the OpenQP state ordering, where the
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
        getattr(lib, f"{self.SYMBOL_PREFIX}_soc_matrix")(
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
            # Only the CURRENT geometry (mol.idx == 1) owns the alignment
            # marker. BasisOverlap.load_previous_data() evaluates the previous
            # geometry with SinglePoint.energy() at mol.idx == 2: that call
            # stores its own MO tags but must NOT advance the marker, or the
            # later current-geometry excitation() would see a stale marker,
            # treat the phase/reorder-aligned tags written by dftb_overlap() as
            # absent, and overwrite them -- corrupting standalone system2
            # NACME/TDC state overlaps. (put_data() restores the current data
            # tags but not this plain attribute.)
            if int(getattr(self.mol, "idx", 1)) == 1:
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
            str(self.dftb.get("backend", "native")).lower(),
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
            # Open-shell ground-state reference: UKS/ROKS/CUKS can share the same
            # derived multiplicity, so the selector must key the cache too, or a
            # ROKS result would be reused for a UKS request at the same geometry.
            getattr(self, "_reference", ""),
            int(getattr(self, "_reference_selector", 0)),
            str(self.dftb.get("response_solver", "auto")).lower(),
            int(self.dftb.get("response_max_subspace", 100)),
            int(self.dftb.get("response_max_iterations", 50)),
            float(self.dftb.get("response_tolerance", 1.0e-6)),
            int(self.dftb.get("print_level", 1)),
            bool(self.dftb.get("state_to_state_spectrum", True)),
            # DTCAM operator surface + preset: a Python workflow may retune
            # these on the same molecule, so cached results must not outlive
            # them.
            str(self.dftb.get("model", "")).strip().lower(),
            float(self.dftb.get("c_mrsf", -1.0)),
            float(self.dftb.get("c_mrsf_oo", -1.0)),
            bool(self.dftb.get("response_global_hybrid", False)),
            float(self.dftb.get("onsite_exchange_scale", 0.0)),
            float(self.dftb.get("w_scale", 1.0)),
            float(self.dftb.get("response_w_scale", -1.0)),
            float(self.dftb.get("response_omega", -1.0)),
            float(self.dftb.get("response_cam_alpha", -1.0)),
            float(self.dftb.get("response_cam_beta", -1.0)),
            float(self.dftb.get("onsite_ss", 0.0)),
            float(self.dftb.get("onsite_sp", 0.0)),
            float(self.dftb.get("onsite_pp", 0.0)),
        )
        return atoms, coords.tobytes(), option_key

    def _parameter_path(self) -> str:
        raw = self.dftb.get("parameter_path") or os.environ.get("OPENQP_DFTB_PARAMETER_PATH")
        if raw:
            return str(self._resolve_user_path(raw))
        bundled = _bundled_parameter_path()
        if bundled:
            return bundled
        raise ValueError(
            "Set [dftb] parameter_path or OPENQP_DFTB_PARAMETER_PATH "
            "(this openqp-dftb installation ships no bundled parameter set)."
        )

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
        return canonical_dftb_type(method)

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

    # --- TB-backend hooks (overridden by the OpenQP-xTB adapter) ----------
    def _abi_version_symbol(self) -> str:
        """C symbol that reports the native ABI version for this backend."""
        return f"{self.SYMBOL_PREFIX}_capi_abi_version"

    def _lc_gamma_code(self) -> int:
        """Long-range gamma kernel as an integer for the C ABI.

        DFTB exposes only yukawa/erf, so this is exactly the historical
        ``lc_gamma_erf`` 0/1 flag (yukawa=0, erf=1); the OpenQP-xTB adapter
        overrides it to add the Ohno-Klopman ``ok`` kernel (=2).
        """
        return 1 if self._lc_gamma_is_erf() else 0

    def _model_args(self) -> list:
        """C ABI model-block scalars spliced in right after the zvector flag.

        Empty for DFTB (the DFTB C ABI has no model block); the OpenQP-xTB
        adapter overrides it to insert model/dispersion/halogen_bond/
        third_order/spin_scale.
        """
        return []

    def _reference_multiplicity(self, method: str) -> int:
        explicit = int(self.dftb.get("reference_multiplicity", 0))
        if canonical_dftb_type(method) in {"sf", "mrsf"}:
            return explicit if explicit > 1 else 3
        return explicit if explicit > 0 else 1

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


def _tb_soc(mol, adapter_cls):
    """Shared SOC driver for the tight-binding backends (dftb/xtb).

    Runs the MRSF response for both target multiplicities from the same ROKS
    reference, builds the one-center SOC matrix, diagonalizes
    diag((E - e0) * ha2wn) + hsoc * (alpha^2/2 * ha2wn) in numpy, and fills the
    same OQP:: tags the native soc_mrsf produces (soc_eval in cm^-1 relative to
    the lowest MCH excitation; soc_hsoc_* raw, without the alpha^2/2 prefactor).
    The backend (its section name and native library) is selected by
    ``adapter_cls`` -- OpenQPDFTBAdapter for dftb, OpenQPXTBAdapter for xtb.
    """
    adapter = adapter_cls(mol)
    data = mol.data
    dftb = mol.config[adapter_cls.SECTION]
    tdhf = mol.config['tdhf']
    saved = dftb.get('target_multiplicity', 1)
    saved_mult = int(tdhf.get('multiplicity', 1))
    saved_nstate = int(tdhf.get('nstate', 1))
    ns = int(tdhf.get('nstate_s', 0)) or saved_nstate
    nt = int(tdhf.get('nstate_t', 0)) or saved_nstate

    def select_manifold(mult, nstate):
        dftb['target_multiplicity'] = mult
        tdhf['multiplicity'] = mult
        tdhf['nstate'] = nstate
        data.set_tdhf_multiplicity(mult)
        data.set_tdhf_nstate(nstate)

    try:
        select_manifold(1, ns)
        singlet_energies = np.array(adapter_cls(mol).energy(), dtype=float)
        data['OQP::td_singlet_energies'] = np.asarray(data['OQP::td_energies']).copy()
        x_s = np.asarray(data['OQP::td_bvec_mo']).copy()
        data['OQP::td_bvec_mo_s'] = x_s.copy()
        mo_tag = np.asarray(data['OQP::VEC_MO_A']).copy()
        dims = np.asarray(data['OQP::dftb_wf_dims']).ravel()

        select_manifold(3, nt)
        triplet_energies = np.array(adapter_cls(mol).energy(), dtype=float)
        data['OQP::td_triplet_energies'] = np.asarray(data['OQP::td_energies']).copy()
        x_t = np.asarray(data['OQP::td_bvec_mo']).copy()
        data['OQP::td_bvec_mo_t'] = x_t.copy()
    finally:
        dftb['target_multiplicity'] = saved
        tdhf['multiplicity'] = saved_mult
        tdhf['nstate'] = saved_nstate
        data.set_tdhf_multiplicity(saved_mult)
        data.set_tdhf_nstate(saved_nstate)

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

    fortran_tag = adapter_cls._fortran_tag
    data['OQP::soc_eval'] = np.ascontiguousarray(eigenvalues.real)
    data['OQP::soc_evec_re'] = fortran_tag(np.ascontiguousarray(eigenvectors.real))
    data['OQP::soc_evec_im'] = fortran_tag(np.ascontiguousarray(eigenvectors.imag))
    data['OQP::soc_hsoc_re'] = fortran_tag(np.ascontiguousarray(hsoc_re))
    data['OQP::soc_hsoc_im'] = fortran_tag(np.ascontiguousarray(hsoc_im))
    mol.soc = eigenvalues.real
    return eigenvalues.real


def dftb_soc(mol):
    """Standalone SOC driver for method=dftb (see _tb_soc)."""
    return _tb_soc(mol, OpenQPDFTBAdapter)
