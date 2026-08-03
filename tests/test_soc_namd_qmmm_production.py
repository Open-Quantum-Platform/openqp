"""Production guards for SOC-NAMD-QMMM option selection and hop bookkeeping.

These tests avoid importing the compiled OpenQP runtime. They pin the Python
dispatch and source-level invariants that make ``soc_basis=mch`` the preferred
SOC-QMMM production path.
"""

from __future__ import annotations

import importlib.util
import re
import sys
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

    modules = {
        "oqp": oqp,
        "oqp.library": library,
        "oqp.library.nac_utils": nac_utils,
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

    def test_soc_restart_round_trip_restores_complex_gauge_and_response_history(self):
        namd, _logs, cleanup = load_namd_with_stubs()
        try:
            driver = namd.NAMD_SOC.__new__(namd.NAMD_SOC)
            driver.nstate_soc = driver.nstate = 4
            driver.prev_u = np.diag([1.0, 1j, -1.0, -1j]).astype(complex)
            driver.prev_eval = np.array([0.0, 0.01, 0.02, 0.03])
            driver.prev_sbvec = np.arange(6.0).reshape(2, 3)
            driver.prev_tbvec = np.arange(6.0, 12.0).reshape(2, 3)

            payload = driver._restart_extra_payload()
            restored = driver._load_restart_extra(payload)
            target = namd.NAMD_SOC.__new__(namd.NAMD_SOC)
            target._restore_restart_extra(restored)

            np.testing.assert_allclose(target.prev_u, driver.prev_u)
            np.testing.assert_allclose(target.prev_eval, driver.prev_eval)
            np.testing.assert_allclose(target.prev_sbvec, driver.prev_sbvec)
            np.testing.assert_allclose(target.prev_tbvec, driver.prev_tbvec)

            broken = dict(payload)
            broken['soc_prev_u_real'] = np.ones((4, 4))
            with self.assertRaisesRegex(RuntimeError, "not unitary"):
                driver._load_restart_extra(broken)
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
            self.assertEqual(namd.NAMD_TRAJECTORY_SCHEMA_VERSION, 6)
            self.assertEqual(namd.NAMD_RESTART_SCHEMA_VERSION, 7)
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

    def test_hop_rng_configuration_rejects_invalid_domains(self):
        src = NAMD.read_text()

        self.assertIn("seed must fit in a signed 64-bit integer", src)
        self.assertIn("rng_stream must be a non-negative signed 64-bit integer", src)
        self.assertIn("first_hop_step must be at least 1", src)

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
