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
                    self.assertTrue(barrier)
                    return value

            root.mol = types.SimpleNamespace(mpi_manager=RootManager())
            with self.assertRaisesRegex(OSError, "disk full"):
                root._run_io_collective(
                    lambda: (_ for _ in ()).throw(OSError("disk full")))

            follower = namd.NAMD.__new__(namd.NAMD)

            class FollowerManager:
                rank = 1
                use_mpi = 1

                @staticmethod
                def bcast(value, root=0, barrier=True):
                    self.assertIsNone(value)
                    self.assertEqual(root, 0)
                    self.assertTrue(barrier)
                    return (False, "OSError", "disk full")

            follower.mol = types.SimpleNamespace(mpi_manager=FollowerManager())
            action_called = False

            def unexpected_action():
                nonlocal action_called
                action_called = True

            with self.assertRaisesRegex(OSError, "disk full"):
                follower._run_io_collective(unexpected_action)
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

        self.assertIn("hopped = self._mch_propagate_and_hop(", body)
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

    def test_soc_restart_round_trip_restores_complex_gauge_and_response_history(self):
        namd, _logs, cleanup = load_namd_with_stubs()
        try:
            driver = namd.NAMD_SOC.__new__(namd.NAMD_SOC)
            driver.nstate_soc = driver.nstate = 4
            driver.prev_u = np.diag([1.0, 1j, -1.0, -1j]).astype(complex)
            driver.prev_eval = np.array([0.0, 0.01, 0.02, 0.03])
            driver.prev_sbvec = np.arange(6.0).reshape(2, 3)
            driver.prev_tbvec = np.arange(6.0, 12.0).reshape(2, 3)
            driver.sbvec = driver.prev_sbvec.copy()
            driver.tbvec = driver.prev_tbvec.copy()
            prev_data = {
                'OQP::td_bvec_mo_s': driver.prev_sbvec.copy(),
                'OQP::td_bvec_mo_t': driver.prev_tbvec.copy(),
            }

            payload = driver._restart_extra_payload()
            restored = driver._load_restart_extra(payload, prev_data)
            target = namd.NAMD_SOC.__new__(namd.NAMD_SOC)
            target._restore_restart_extra(restored)

            np.testing.assert_allclose(target.prev_u, driver.prev_u)
            np.testing.assert_allclose(target.prev_eval, driver.prev_eval)
            np.testing.assert_allclose(target.prev_sbvec, driver.prev_sbvec)
            np.testing.assert_allclose(target.prev_tbvec, driver.prev_tbvec)

            broken = dict(payload)
            broken['soc_prev_u_real'] = np.ones((4, 4))
            with self.assertRaisesRegex(RuntimeError, "not unitary"):
                driver._load_restart_extra(broken, prev_data)

            wrong_shape = dict(payload)
            wrong_shape['soc_prev_sbvec'] = np.ones((3, 2))
            with self.assertRaisesRegex(RuntimeError, "singlet.*shape or dtype"):
                driver._load_restart_extra(wrong_shape, prev_data)

            wrong_dtype = dict(payload)
            wrong_dtype['soc_prev_tbvec'] = np.asarray(
                payload['soc_prev_tbvec'], dtype=np.float32)
            with self.assertRaisesRegex(RuntimeError, "triplet.*shape or dtype"):
                driver._load_restart_extra(wrong_dtype, prev_data)
        finally:
            cleanup()

    def test_soc_overlap_generator_is_complex_antihermitian(self):
        from scipy.linalg import expm

        namd, _logs, cleanup = load_namd_with_stubs()
        try:
            driver = namd.NAMD_SOC.__new__(namd.NAMD_SOC)
            driver.dt = 2.5
            generator = np.array(
                [[0.2j, 0.03 + 0.04j], [-0.03 + 0.04j, -0.1j]],
                dtype=complex,
            )
            overlap = expm(generator)
            unitary, recovered = driver._soc_unitary_overlap(overlap)

            np.testing.assert_allclose(
                unitary.conj().T @ unitary, np.eye(2), atol=1.0e-12)
            np.testing.assert_allclose(
                recovered + recovered.conj().T, 0.0, atol=1.0e-12)
            np.testing.assert_allclose(driver._last_state_overlap, overlap)
            np.testing.assert_allclose(driver._last_overlap_tdc,
                                       recovered / driver.dt)
            self.assertGreater(
                np.max(np.abs(driver._last_overlap_tdc.imag)), 0.0)
        finally:
            cleanup()

    def test_soc_packed_trajectory_keeps_complex_overlap_components(self):
        namd, _logs, cleanup = load_namd_with_stubs()
        try:
            dtype = namd._namd_trajectory_dtype(4, 3)
            self.assertIn('state_overlap_imag', dtype.names)
            self.assertIn('overlap_tdc_imag_au', dtype.names)
            self.assertEqual(namd.NAMD_TRAJECTORY_SCHEMA_VERSION, 7)
            self.assertEqual(namd.NAMD_RESTART_SCHEMA_VERSION, 8)
        finally:
            cleanup()

    def test_all_soc_drivers_wire_output_restart_and_nve_gate(self):
        src = NAMD.read_text()
        for class_name, end_name in (
            ('NAMD_SOC', 'NAMD_SOC_MCH'),
            ('NAMD_SOC_MCH', 'NAMD_SOC_QMMM'),
            ('NAMD_SOC_QMMM', 'NAMD_SOC_MCH_QMMM'),
        ):
            block = src.split(f'class {class_name}', 1)[1].split(
                f'class {end_name}', 1)[0]
            self.assertIn('self._prepare_md_outputs()', block)
            self.assertIn('restart = self._load_restart()', block)
            self.assertIn('self._save_restart(', block)
            self.assertIn('self._update_nve_gate(', block)
            self.assertIn('self._write_md_trajectory(', block)
        final_block = src.split('class NAMD_SOC_MCH_QMMM', 1)[1]
        self.assertIn('self._prepare_md_outputs()', final_block)
        self.assertIn('restart = self._load_restart()', final_block)
        self.assertIn('self._save_restart(', final_block)
        self.assertIn('self._update_nve_gate(', final_block)
        self.assertIn('self._write_md_trajectory(', final_block)

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

    def test_delayed_hop_still_propagates_all_electronic_coefficients(self):
        namd, _logs, cleanup = load_namd_with_stubs()
        try:
            # Same-spin path: the native call must still occur with a negative
            # allow-hop flag, while no stochastic random number is consumed.
            driver = namd.NAMD.__new__(namd.NAMD)
            driver.nstate = 2
            driver.natom = 1
            driver.dt_fs = 0.5
            driver.dt = 0.5 * namd.FS_TO_AU
            driver.substep = 4
            driver.thrshe = 0.1
            driver.active = 1
            driver.decoherence = 0
            driver.edc_c = 0.1
            driver.tdc_scheme = 0
            driver.trivial = 0
            driver.trivial_thresh = 0.5
            driver.coef = np.array([1.0 + 0.0j, 0.0 + 0.0j])
            driver.vel = np.zeros((1, 3))
            driver._hop_random_override = lambda: self.fail(
                "delayed hop consumed a random number")
            driver._last_hop_random = np.nan
            driver.mol = types.SimpleNamespace(data={
                "OQP::td_states_overlap": np.eye(2),
                "OQP::td_energies": np.array([0.0, 0.1]),
            })
            driver._compute_tdc = lambda _overlap: np.zeros((2, 2))

            native_calls = []

            def propagate_only(mol):
                native_calls.append(np.array(
                    mol.data["OQP::namd_params"], copy=True))
                mol.data["OQP::namd_coef"] = np.array(
                    [np.sqrt(0.75), 0.0, 0.0, 0.5])

            namd.oqp.mrsf_namd_hop = propagate_only
            new_active, hopped = driver._hop(allow_hop=False)
            self.assertEqual(len(native_calls), 1)
            self.assertLess(native_calls[0][namd._P_ALLOW_HOP], 0.0)
            np.testing.assert_allclose(
                driver.coef, [np.sqrt(0.75), 0.5j])
            self.assertEqual(new_active, 1)
            self.assertFalse(hopped)
            self.assertTrue(np.isnan(driver._last_hop_random))

            # SOC adiabatic path: coefficient propagation likewise precedes
            # the allow-hop return and must not consume the RNG.
            soc = namd.NAMD_SOC.__new__(namd.NAMD_SOC)
            soc.nstate_soc = 2
            soc.active = 1
            soc.dt = 0.2
            soc.substep = 1
            soc.coef = np.array([1.0 + 0.0j, 0.0 + 0.0j])
            soc.decoherence = 0
            soc.mass = np.ones(1)
            soc.vel = np.zeros((1, 3))
            soc.thrshe = 0.1
            soc._hop_random_override = lambda: self.fail(
                "delayed SOC hop consumed a random number")
            soc._last_hop_random = np.nan
            theta = 0.2
            overlap = np.array([
                [np.cos(theta), -np.sin(theta)],
                [np.sin(theta), np.cos(theta)],
            ], dtype=complex)
            hopped = soc._propagate_and_hop(
                np.zeros(2), np.zeros(2), overlap, allow_hop=False)
            self.assertFalse(hopped)
            self.assertEqual(soc.active, 1)
            self.assertGreater(abs(soc.coef[1]), 0.1)
            self.assertTrue(np.isnan(soc._last_hop_random))

            # MCH-basis path follows the same contract.
            mch = namd.NAMD_SOC_MCH.__new__(namd.NAMD_SOC_MCH)
            mch.nstate_soc = 2
            mch.active = 1
            mch.dt = 0.2
            mch.coef = np.array([1.0 + 0.0j, 0.0 + 0.0j])
            mch.decoherence = 0
            mch.mass = np.ones(1)
            mch.vel = np.zeros((1, 3))
            mch.thrshe = 0.1
            mch._hop_random_override = lambda: self.fail(
                "delayed MCH hop consumed a random number")
            mch._last_hop_random = np.nan
            hopped = mch._mch_propagate_and_hop(
                np.array([[0.0, 0.2], [0.2, 0.0]]),
                np.zeros(2), allow_hop=False)
            self.assertFalse(hopped)
            self.assertEqual(mch.active, 1)
            self.assertGreater(abs(mch.coef[1]), 0.01)
            self.assertTrue(np.isnan(mch._last_hop_random))
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
                self.assertEqual(
                    first, str(Path(tmpdir) / "window-01.namd.restart.oqp"))
                self.assertEqual(
                    second, str(Path(tmpdir) / "window-02.namd.restart.oqp"))
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

    def test_soc_nve_gate_and_restart_are_enabled(self):
        src = NAMD.read_text()

        self.assertNotIn("nve_gate currently supports same-spin NAMD only", src)
        self.assertNotIn("restart currently supports same-spin NAMD only", src)
        self.assertIn("'dt_adaptive', 'dt_min', 'dx_max', 'econs'", src)

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
                        "md": {
                            "soc": True, "soc_basis": "adiabatic",
                            "econs": False,
                        },
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
                driver._restart_molecular_identity = {
                    "kind": "test", "sha256": "soc-molecule",
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
                mol.config["md"]["econs"] = True
                self.assertNotEqual(
                    driver._restart_signature(), frozen_signature)
                mol.config["md"]["econs"] = False
                mol.config["tdhf"]["multiplicity"] = 1
                self.assertNotEqual(
                    driver._restart_signature(), frozen_signature)
                mol.config["tdhf"]["multiplicity"] = 3
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

    def test_native_qmmm_reconstructs_petite_gradients_before_force_assembly(self):
        src = QMMM_DRIVER.read_text()

        self.assertIn("gqm = np.asarray(gradient.gradient(), dtype=float)", src)
        self.assertIn("gradient.mol.symmetrize_gradient(gqm)", src)
        self.assertIn("gradient.mol.set_grad(gqm)", src)

        native_start = src.index("# --- Gradients: pure QM + ESPF contribution")
        native_end = src.index("self.op.mol.save_data()", native_start)
        native_path = src[native_start:native_end]
        self.assertNotIn("oqp.hf_gradient(self.op.mol)", native_path)

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
