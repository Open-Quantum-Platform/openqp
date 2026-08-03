"""Production guards for SOC-NAMD-QMMM option selection and hop bookkeeping.

These tests avoid importing the compiled OpenQP runtime. They pin the Python
dispatch and source-level invariants that make ``soc_basis=mch`` the preferred
SOC-QMMM production path.
"""

from __future__ import annotations

import importlib.util
import re
import sys
import tempfile
import types
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RUNFUNC = ROOT / "pyoqp" / "oqp" / "library" / "runfunc.py"
NAMD = ROOT / "pyoqp" / "oqp" / "library" / "namd.py"
OQPDATA = ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"
PYOQP = ROOT / "pyoqp" / "oqp" / "pyoqp.py"
QMMM_DRIVER = ROOT / "pyoqp" / "oqp" / "library" / "qmmm_driver.py"
MRSF_Z_VECTOR = ROOT / "source" / "modules" / "tdhf_mrsf_z_vector.F90"
SOC_MRSF = ROOT / "source" / "modules" / "soc_mrsf.F90"


def load_runfunc_with_namd_stubs():
    """Import runfunc.py with fake NAMD classes and no compiled oqp module."""
    calls = []

    class _Runner:
        def __init__(self, mol):
            self.mol = mol

        def run(self):
            calls.append(type(self).__name__)

    namd = types.ModuleType("oqp.library.namd")
    for name in (
        "NAMD",
        "NAMD_QMMM",
        "NAMD_SOC",
        "NAMD_SOC_QMMM",
        "NAMD_SOC_MCH",
        "NAMD_SOC_MCH_QMMM",
    ):
        setattr(namd, name, type(name, (_Runner,), {}))

    oqp = types.ModuleType("oqp")
    oqp.__path__ = []
    for name in ("electric_moments", "mulliken", "lowdin", "resp_charges"):
        setattr(oqp, name, lambda *_args, **_kwargs: None)

    library = types.ModuleType("oqp.library")
    library.__path__ = []
    utils = types.ModuleType("oqp.utils")
    utils.__path__ = []
    tb_backends = types.ModuleType("oqp.utils.tb_backends")
    setattr(tb_backends, "is_tb_method", lambda *_args, **_kwargs: False)
    setattr(tb_backends, "tb_section_name", lambda *_args, **_kwargs: "dftb")

    single_point = types.ModuleType("oqp.library.single_point")
    noop = type("_Noop", (), {
        "__init__": lambda self, *_args, **_kwargs: None,
        "energy": lambda self: None,
        "gradient": lambda self: None,
        "hessian": lambda self: None,
        "compute": lambda self, *_args, **_kwargs: None,
        "overlap": lambda self: None,
        "nacme": lambda self: None,
        "nac": lambda self: None,
    })
    for name in ("SinglePoint", "Gradient", "Hessian", "LastStep", "BasisOverlap", "NACME", "NAC"):
        setattr(single_point, name, noop)

    modules = {
        "oqp": oqp,
        "oqp.library": library,
        "oqp.library.namd": namd,
        "oqp.library.single_point": single_point,
        "oqp.library.libscipy": types.ModuleType("oqp.library.libscipy"),
        "oqp.library.libgeometric": types.ModuleType("oqp.library.libgeometric"),
        "oqp.library.liboqp": types.ModuleType("oqp.library.liboqp"),
        "oqp.utils": utils,
        "oqp.utils.tb_backends": tb_backends,
    }
    for module_name, names in {
        "oqp.library.libscipy": ("StateSpecificOpt", "MECIOpt", "MECPOpt", "MEP", "QMMMOpt"),
        "oqp.library.libgeometric": (
            "GeometricIRCOpt",
            "GeometricMECIOpt",
            "GeometricMECPOpt",
            "GeometricNEBOpt",
            "GeometricOpt",
            "GeometricTSOpt",
            "GeometricQMMMOpt",
        ),
        "oqp.library.liboqp": (
            "OQPOpt",
            "OQPTSOpt",
            "OQPMECIOpt",
            "OQPMECPOpt",
            "OQPTCIOpt",
            "OQPNEBOpt",
            "OQPIRCOpt",
            "OQPMEPOpt",
        ),
    }.items():
        for name in names:
            setattr(modules[module_name], name, noop)

    saved = {name: sys.modules.get(name) for name in modules}
    sys.modules.update(modules)

    def cleanup():
        for name, value in saved.items():
            if value is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = value

    spec = importlib.util.spec_from_file_location("runfunc_soc_namd_under_test", RUNFUNC)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, calls, cleanup


def load_namd_with_stubs():
    """Import namd.py without the compiled extension for state-map tests."""
    logs = []
    oqp = types.ModuleType("oqp")
    oqp.__path__ = []
    setattr(
        oqp,
        "oqp_namd_counter_random",
        lambda seed, stream, step: ((int(seed) + 17 * int(stream) + 31 * int(step)) % 997) / 997.0,
    )
    setattr(oqp, "oqp_namd_counter_normal_fill", lambda *_args: None)
    setattr(oqp, "ffi", types.SimpleNamespace(cast=lambda _ctype, ptr: ptr))
    library = types.ModuleType("oqp.library")
    library.__path__ = []
    utils = types.ModuleType("oqp.utils")
    utils.__path__ = []
    single_point = types.ModuleType("oqp.library.single_point")
    noop = type("_Noop", (), {})
    for name in ("SinglePoint", "Gradient", "LastStep", "BasisOverlap", "NACME"):
        setattr(single_point, name, noop)
    file_utils = types.ModuleType("oqp.utils.file_utils")
    setattr(file_utils, "dump_log", lambda *_args, **kwargs: logs.append(kwargs.get("title")))
    nac_utils = types.ModuleType("oqp.library.nac_utils")
    setattr(
        nac_utils,
        "canonical_state_overlap",
        lambda overlap: np.array(np.asarray(overlap).T, copy=True),
    )
    tb_backends = types.ModuleType("oqp.utils.tb_backends")
    setattr(tb_backends, "is_tb_method", lambda *_args, **_kwargs: False)
    setattr(tb_backends, "make_tb_adapter", lambda *_args, **_kwargs: None)
    setattr(tb_backends, "tb_section_name", lambda *_args, **_kwargs: "dftb")
    odp = types.ModuleType("oqp.library.odp")
    setattr(odp, "odp_from_config", lambda *_args, **_kwargs: None)

    modules = {
        "oqp": oqp,
        "oqp.library": library,
        "oqp.library.nac_utils": nac_utils,
        "oqp.library.odp": odp,
        "oqp.library.single_point": single_point,
        "oqp.utils": utils,
        "oqp.utils.file_utils": file_utils,
        "oqp.utils.tb_backends": tb_backends,
    }
    saved = {name: sys.modules.get(name) for name in modules}
    sys.modules.update(modules)

    def cleanup():
        for name, value in saved.items():
            if value is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = value

    spec = importlib.util.spec_from_file_location("namd_state_map_under_test", NAMD)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, logs, cleanup


def load_pyoqp_with_stubs():
    """Import pyoqp.py without native/OpenMM dependencies for Runner tests."""
    oqp = types.ModuleType("oqp")
    oqp.__path__ = []
    library = types.ModuleType("oqp.library")
    library.__path__ = []
    utils = types.ModuleType("oqp.utils")
    utils.__path__ = []

    file_utils = types.ModuleType("oqp.utils.file_utils")
    setattr(file_utils, "dump_log", lambda *_args, **_kwargs: None)
    input_checker = types.ModuleType("oqp.utils.input_checker")
    setattr(input_checker, "check_input_values", lambda *_args, **_kwargs: None)
    molecule = types.ModuleType("oqp.molecule")
    setattr(molecule, "Molecule", type("Molecule", (), {}))
    runfunc = types.ModuleType("oqp.library.runfunc")
    for name in (
        "compute_energy", "compute_grad", "compute_nac", "compute_soc",
        "compute_geom", "compute_md", "compute_nacme", "compute_properties",
        "compute_data", "compute_hess", "compute_thermo", "compute_namd",
    ):
        setattr(runfunc, name, lambda *_args, **_kwargs: None)
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")
    setattr(mpi_utils, "MPIManager", type("MPIManager", (), {"use_mpi": 0, "rank": 0}))

    modules = {
        "oqp": oqp,
        "oqp.library": library,
        "oqp.library.runfunc": runfunc,
        "oqp.utils": utils,
        "oqp.utils.file_utils": file_utils,
        "oqp.utils.input_checker": input_checker,
        "oqp.utils.mpi_utils": mpi_utils,
        "oqp.molecule": molecule,
    }
    saved = {name: sys.modules.get(name) for name in modules}
    sys.modules.update(modules)

    def cleanup():
        for name, value in saved.items():
            if value is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = value

    spec = importlib.util.spec_from_file_location("pyoqp_runner_under_test", PYOQP)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, cleanup


class DummyMol:
    def __init__(self, *, qmmm, soc, soc_basis="adiabatic"):
        self.config = {
            "input": {"qmmm_flag": qmmm},
            "md": {"soc": soc, "soc_basis": soc_basis},
        }


class SOCNAMDQMMMProductionTests(unittest.TestCase):
    def test_rank_zero_io_failure_is_broadcast_to_every_mpi_rank(self):
        namd, _logs, cleanup = load_namd_with_stubs()
        try:
            root = namd.NAMD.__new__(namd.NAMD)

            class RootManager:
                rank = 0
                use_mpi = 1

                @staticmethod
                def bcast(value, root=0, barrier=True):
                    self.assertEqual(root, 0)
                    self.assertFalse(barrier)
                    return value

            root.mol = types.SimpleNamespace(mpi_manager=RootManager())
            with self.assertRaisesRegex(
                RuntimeError, "checkpoint write failed on MPI rank 0: OSError: disk full"
            ):
                root._run_io_rank(
                    lambda: (_ for _ in ()).throw(OSError("disk full")),
                    "checkpoint write",
                )

            follower = namd.NAMD.__new__(namd.NAMD)

            class FollowerManager:
                rank = 1
                use_mpi = 1

                @staticmethod
                def bcast(value, root=0, barrier=True):
                    self.assertIsNone(value)
                    self.assertEqual(root, 0)
                    self.assertFalse(barrier)
                    return {"ok": False, "error": "OSError: disk full"}

            follower.mol = types.SimpleNamespace(mpi_manager=FollowerManager())
            action_called = False

            def unexpected_action():
                nonlocal action_called
                action_called = True

            with self.assertRaisesRegex(
                RuntimeError, "checkpoint write failed on MPI rank 0: OSError: disk full"
            ):
                follower._run_io_rank(unexpected_action, "checkpoint write")
            self.assertFalse(action_called)
        finally:
            cleanup()

    def test_soc_basis_mch_dispatches_to_mch_qmmm_class(self):
        runfunc, calls, cleanup = load_runfunc_with_namd_stubs()
        try:
            runfunc.compute_namd(DummyMol(qmmm=True, soc=True, soc_basis="mch"))
        finally:
            cleanup()

        self.assertEqual(calls, ["NAMD_SOC_MCH_QMMM"])

    def test_adiabatic_soc_qmmm_dispatch_remains_available(self):
        runfunc, calls, cleanup = load_runfunc_with_namd_stubs()
        try:
            runfunc.compute_namd(DummyMol(qmmm=True, soc=True, soc_basis="adiabatic"))
        finally:
            cleanup()

        self.assertEqual(calls, ["NAMD_SOC_QMMM"])

    def test_soc_basis_schema_documents_mch_and_correction_switches(self):
        src = OQPDATA.read_text()

        self.assertIn("'soc_basis'", src)
        self.assertIn("'adiabatic' (SHARC) | 'mch' (spin-pure exact-gradient)", src)
        self.assertIn("'soc_du_dt_corr'", src)
        self.assertIn("'soc_tdc_grad_corr'", src)

    def test_mrsf_triplet_labels_are_zero_based_in_soc_and_namd(self):
        namd = NAMD.read_text()
        soc = SOC_MRSF.read_text()

        self.assertIn("base = self.ns + n * 3", namd)
        self.assertIn("return f'S{target - 1}' if mult == 1 else f'T{target - 1}'", namd)
        self.assertIn("active = self.ns + target * 3 + 1", namd)
        self.assertIn("'S', ist-1, 'T', jst-1", soc)
        self.assertIn("'T', ist-1, trim(trip(ims_i)), 'T', jst-1", soc)
        self.assertIn("'T', (i-ns-1)/3, '( 0)'", soc)

    def test_soc_namd_preserves_legacy_t1_and_public_t0_state_maps(self):
        namd, logs, cleanup = load_namd_with_stubs()
        try:
            self.assertEqual(namd._parse_soc_init_state("T1", 1, 1), (3, 0, "T0"))
            self.assertEqual(
                namd._parse_soc_init_state("T0", 1, 1, public_labels=True),
                (3, 0, "T0"),
            )
            self.assertEqual(
                namd._parse_soc_init_state("T1", 2, 2, public_labels=True),
                (3, 1, "T1"),
            )
            with self.assertRaisesRegex(ValueError, "outside the SOC MCH basis"):
                namd._parse_soc_init_state("T1", 1, 1, public_labels=True)

            # Legacy T1 must select the first triplet block in both the
            # spin-adiabatic character resolver and the explicit MCH path.
            adiabatic = namd.NAMD_SOC.__new__(namd.NAMD_SOC)
            adiabatic.init_state = "T1"
            adiabatic.ns = adiabatic.nt = 1
            adiabatic.nstate_soc = 4
            adiabatic.mol = types.SimpleNamespace(oqp_public_state_labels=False)
            adiabatic.coef = np.zeros(4, dtype=complex)
            adiabatic._resolve_initial_active(np.eye(4, dtype=complex))
            self.assertEqual(adiabatic.active, 2)
            self.assertEqual(np.flatnonzero(adiabatic.coef).tolist(), [1])

            mch = namd.NAMD_SOC_MCH.__new__(namd.NAMD_SOC_MCH)
            mch.init_state = "T1"
            mch.ns = mch.nt = 1
            mch.nstate_soc = 4
            mch.mol = types.SimpleNamespace(oqp_public_state_labels=False)
            mch.coef = np.zeros(4, dtype=complex)
            mch._resolve_initial_mch_active()
            self.assertEqual(mch.active, 2)
            self.assertEqual(np.flatnonzero(mch.coef).tolist(), [1])
            self.assertTrue(any("T0 (legacy T1)" in message for message in logs))
        finally:
            cleanup()

    def test_programmatic_runner_dispatches_qmmm_md_to_openmm_driver(self):
        pyoqp, cleanup = load_pyoqp_with_stubs()
        calls = []

        class QMMM_MD:
            def __init__(self, *, mol):
                calls.append(("init", mol))

            def run(self):
                calls.append(("run", None))

        qmmm_md = types.ModuleType("oqp.library.qmmm_md")
        setattr(qmmm_md, "QMMM_MD", QMMM_MD)
        previous = sys.modules.get("oqp.library.qmmm_md")
        sys.modules["oqp.library.qmmm_md"] = qmmm_md
        try:
            mol = types.SimpleNamespace(
                config={
                    "input": {"runtype": "md", "qmmm_flag": True},
                    "tests": {"exception": False},
                },
                usempi=False,
            )
            runner = pyoqp.Runner.__new__(pyoqp.Runner)
            runner.mol = mol
            runner.mpi_manager = types.SimpleNamespace(use_mpi=0)
            runner.run_func = {"md": lambda _mol: self.fail("compute_md was called")}

            runner.run()

            self.assertEqual(calls, [("init", mol), ("run", None)])
            self.assertIsInstance(runner.qmmm_md, QMMM_MD)
        finally:
            if previous is None:
                sys.modules.pop("oqp.library.qmmm_md", None)
            else:
                sys.modules["oqp.library.qmmm_md"] = previous
            cleanup()

    def test_soc_qmmm_hops_rescale_for_espf_energy_change(self):
        src = NAMD.read_text()

        self.assertGreaterEqual(src.count("de_espf = ((epot_old - epot) +"), 2)
        self.assertIn("eval_ha[self.active - 1] - eval_ha[active_old - 1]", src)
        self.assertIn("e_mch[self.active - 1] - e_mch[active_old - 1]", src)
        self.assertIn("self.v_all *= np.sqrt(max(0.0, 1.0 + de_espf / ekin_all))", src)

    def test_mch_qmmm_recomputes_exact_gradient_after_accepted_hop(self):
        src = NAMD.read_text()
        block = re.search(
            r"class NAMD_SOC_MCH_QMMM\(.*?dump_log\(mol, title='PyOQP: SOC-MCH-QMMM-NAMD trajectory complete'\)",
            src,
            re.S,
        )
        self.assertIsNotNone(block)
        body = block.group(0)

        self.assertIn("hopped = self._mch_propagate_and_hop(h_mch, e_mch)", body)
        self.assertIn("if hopped:", body)
        self.assertIn("g_qm, e_pure, mult, state, pchg = self._mch_exact_gradient_qmmm(self.active)", body)
        self.assertIn("f_all, epot = self._total_force_soc(potmm, g_qm, e_pure, pchg)", body)

    def test_overlap_tracking_and_rng_are_reproducible_contracts(self):
        src = NAMD.read_text()

        self.assertIn("oqp.oqp_namd_counter_random(", src)
        self.assertIn("oqp.oqp_namd_counter_normal_fill(", src)
        self.assertIn("self._rng_step >= self.first_hop_step", src)
        self.assertNotIn("np.random.default_rng", src)
        self.assertIn("cfg['properties']['back_door'] = True", src)
        self.assertIn("NAMD_SOC._store_prev(self, self.r_all[self.qm_atoms].reshape((self.natom, 3)), u, eval_ha)", src)
        self.assertIn("BasisOverlap(mol).overlap()", src)
        self.assertIn("s_mch = NAMD_SOC._mch_overlap(self)", src)

    def test_hop_rng_is_step_indexed_and_supports_exact_replay(self):
        namd, _logs, cleanup = load_namd_with_stubs()
        try:
            driver = namd.NAMD.__new__(namd.NAMD)
            driver.seed = 11
            driver.rng_stream = 7
            driver.first_hop_step = 2
            driver._rng_step = 0
            driver._last_hop_random = np.nan

            self.assertFalse(driver._prepare_hop_step(1))
            self.assertTrue(np.isnan(driver._last_hop_random))
            self.assertTrue(driver._prepare_hop_step(2))
            expected = ((11 + 17 * 7 + 31 * 2) % 997) / 997.0
            self.assertEqual(driver._hop_random(), expected)
            self.assertEqual(driver._hop_random(), expected)

            class Replay:
                calls = 0

                def random(self):
                    self.calls += 1
                    return 0.0605

            replay = Replay()
            driver._hop_random_override = replay
            driver._prepare_hop_step(2)
            self.assertEqual(driver._hop_random(), 0.0605)
            self.assertEqual(driver._hop_random(), 0.0605)
            self.assertEqual(replay.calls, 1)
            self.assertEqual(
                driver._hop_rng_log(),
                "rng_step=2 random=0.060499999999999998",
            )
        finally:
            cleanup()

    def test_hop_rng_configuration_rejects_invalid_domains(self):
        src = NAMD.read_text()

        self.assertIn("seed must fit in a signed 64-bit integer", src)
        self.assertIn("rng_stream must be a non-negative signed 64-bit integer", src)
        self.assertIn("first_hop_step must be at least 1", src)

    def test_gate_tolerances_reject_nan_inf_and_negative_values(self):
        namd, _logs, cleanup = load_namd_with_stubs()
        try:
            namd._validate_gate_tolerances("NVE", (0.0, 1.0e-6, 2.0))
            for invalid in (np.nan, np.inf, -np.inf, -1.0e-6):
                with self.assertRaisesRegex(ValueError, "finite and non-negative"):
                    namd._validate_gate_tolerances("NVE", (0.0, invalid, 1.0))
        finally:
            cleanup()

    def test_restart_skips_fresh_velocity_source_initialization(self):
        namd, _logs, cleanup = load_namd_with_stubs()
        try:
            driver = namd.NAMD.__new__(namd.NAMD)
            driver.restart_requested = True
            driver.natom = 2
            driver.velocity_source = "/missing/or/moved/velocity.dat"
            np.testing.assert_array_equal(
                driver._init_velocities(), np.zeros((2, 3)))

            driver.restart_requested = False
            with self.assertRaisesRegex(ValueError, "not zero/maxwell"):
                driver._init_velocities()
        finally:
            cleanup()

    def test_restart_manifests_are_job_specific(self):
        namd, _logs, cleanup = load_namd_with_stubs()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                first = namd._restart_manifest_path(
                    str(Path(tmpdir) / "window-01.log"))
                second = namd._restart_manifest_path(
                    str(Path(tmpdir) / "window-02.log"))
                self.assertEqual(first, str(Path(tmpdir) / "window-01.restart.oqp"))
                self.assertEqual(second, str(Path(tmpdir) / "window-02.restart.oqp"))
                self.assertNotEqual(first, second)
        finally:
            cleanup()

    def test_periodic_odp_is_rejected_until_minimum_images_exist(self):
        namd, _logs, cleanup = load_namd_with_stubs()
        try:
            self.assertIsNone(
                namd._validate_odp_boundary_conditions(object(), False))
            with self.assertRaisesRegex(NotImplementedError, "minimum-image"):
                namd._validate_odp_boundary_conditions(object(), True)
            src = NAMD.read_text()
            self.assertIn(
                "_validate_odp_boundary_conditions(self.odp, self.periodic)",
                src,
            )
        finally:
            cleanup()

    def test_soc_nve_gate_is_explicitly_rejected(self):
        src = NAMD.read_text()

        self.assertIn("if soc_requested and self.nve_gate != 'off':", src)
        self.assertIn("nve_gate currently supports same-spin NAMD only", src)

    def test_soc_dense_trajectory_uses_full_electronic_state_space(self):
        namd, _logs, cleanup = load_namd_with_stubs()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                mol = types.SimpleNamespace(
                    config={
                        "input": {
                            "method": "tdhf", "functional": "bhhlyp",
                            "basis": "sto-3g", "charge": 0,
                        },
                        "scf": {"type": "rohf", "multiplicity": 3},
                        "tdhf": {
                            "type": "mrsf", "multiplicity": 3,
                            "nstate": 1, "tlf": 2,
                        },
                        "md": {"soc": True, "soc_basis": "adiabatic"},
                    },
                    data={"OQP::td_energies": np.array([-1.0])},
                    get_state_tracking=lambda: {
                        "order": [0], "phase_step": [1.0],
                        "matched_overlap": [1.0], "margin": [1.0],
                    },
                )
                driver = namd.NAMD_SOC.__new__(namd.NAMD_SOC)
                driver.mol = mol
                driver.nstate = 1
                driver.coef = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
                driver._trajectory_state_energies = np.array(
                    [-1.0, -0.9, -0.9, -0.9])
                driver._trajectory_representation = "soc_adiabatic"
                driver.econs = True
                driver._restart_system_identity = {
                    "kind": "test", "sha256": "soc-system",
                }
                driver.trajectory_file = str(Path(tmpdir) / "soc.namd.trj")
                driver.trajectory_interval = 1
                driver.dt_fs = 0.5
                driver.dt_adaptive = False
                driver._t_fs = 0.0
                driver.seed = 1
                driver.rng_stream = 0
                driver.init_temp = 300.0
                driver.active = 1
                driver.vel = np.zeros((2, 3))
                driver.odp = None
                driver._last_hop_random = np.nan
                driver._unbiased_potential_energy = -1.0
                driver._last_state_overlap = None
                driver._last_overlap_tdc = None
                driver._nacme_reference_tdc = None
                driver._nacme_reference_mask = None
                driver._nacme_reference_source = 0
                driver._nacme_gate_last = None
                driver._nve_gate_last = None
                driver._pending_nve_gate_error = None
                driver._electronic_config_identity = (
                    namd._electronic_config_identity(mol.config))
                frozen_signature = driver._restart_signature()
                mol.config["tdhf"]["multiplicity"] = 1
                self.assertEqual(
                    driver._restart_signature(), frozen_signature)

                driver._write_md_trajectory(
                    0, np.zeros((2, 3)), -1.0, 0.1, False)
                header, records = namd.read_namd_trajectory(
                    driver.trajectory_file)

                self.assertEqual(header["nstate"], 4)
                self.assertEqual(
                    header["electronic_representation"], "soc_adiabatic")
                self.assertEqual(
                    header["ensemble"],
                    "ENERGY_CONSTRAINED_VELOCITY_RESCALING")
                self.assertTrue(
                    header["ensemble_provenance"][
                        "per_step_velocity_rescaling"])
                self.assertEqual(
                    header["ensemble_provenance"]["velocity_rescaling_mode"],
                    "restore_initial_total_energy")
                self.assertEqual(header["signature"], frozen_signature)
                np.testing.assert_allclose(
                    records["state_energies"][0], driver._trajectory_state_energies)
                np.testing.assert_allclose(
                    records["populations"][0], [1.0, 0.0, 0.0, 0.0])
                self.assertEqual(int(records["tracking_valid"][0]), 0)
        finally:
            cleanup()

    def test_all_soc_run_paths_prepare_and_write_dense_trajectories(self):
        src = NAMD.read_text()

        self.assertEqual(src.count("self._prepare_md_outputs()"), 6)
        for name in ("_log_soc", "_log_mch", "_log_soc_qmmm", "_log_mch_qmmm"):
            block = re.search(
                rf"    def {name}\(.*?(?=\n    def |\n\nclass |\n\ndef )",
                src, re.S)
            self.assertIsNotNone(block, name)
            self.assertIn("self._write_md_trajectory(", block.group(0), name)

    def test_qmmm_single_point_inputs_stay_on_runner_path(self):
        src = PYOQP.read_text()

        self.assertIn("if qmmm_flag and runtype_l == 'md':", src)
        self.assertNotIn("if qmmm_flag and runtype_l != 'namd':", src)

    def test_mm_mm_nonbonded_exceptions_are_preserved(self):
        src = QMMM_DRIVER.read_text()

        self.assertIn("qm_set = set(int(i) for i in qm_atoms)", src)
        self.assertIn("if (int(p1) in qm_set) or (int(p2) in qm_set):", src)
        self.assertNotIn("if (p1 not in qm_atoms) or (p2 not in qm_atoms):", src)

    def test_trivial_crossing_active_state_updates_are_applied(self):
        src = NAMD.read_text()

        self.assertGreaterEqual(src.count("active_changed = new_active != active_old"), 2)
        self.assertIn("kernel can update ACTIVE without marking HOPPED", src)

    def test_espf_charge_fit_is_qmmm_only_in_mrsf_z_vector(self):
        src = MRSF_Z_VECTOR.read_text()

        block = re.search(
            r"if \(infos%control%qmmm_flag\) then\s+call mulliken_excited\(infos\)\s+call form_esp_charges_excited\(infos\)\s+end if",
            src,
            re.S,
        )
        self.assertIsNotNone(block)


if __name__ == "__main__":
    unittest.main()
