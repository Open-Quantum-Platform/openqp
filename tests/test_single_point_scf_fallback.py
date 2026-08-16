import importlib.util
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def install_single_point_stubs():
    oqp = types.ModuleType("oqp")
    oqp.__path__ = []
    setattr(oqp, "lib", types.SimpleNamespace())
    setattr(oqp, "hf_energy", lambda mol: None)
    setattr(oqp, "tdhf_energy", lambda mol: None)
    setattr(oqp, "tdhf_sf_energy", lambda mol: None)
    setattr(oqp, "tdhf_mrsf_energy", lambda mol: None)
    setattr(oqp, "tdhf_umrsf_energy", lambda mol: None)
    setattr(oqp, "library", types.SimpleNamespace(
        set_basis=lambda mol: None,
        ints_1e=lambda mol: None,
        guess=lambda mol: None,
        project_basis=lambda mol: None,
    ))
    sys.modules["oqp"] = oqp

    molecule = types.ModuleType("oqp.molecule")
    setattr(molecule, "Molecule", type("Molecule", (), {}))
    sys.modules["oqp.molecule"] = molecule

    utils = types.ModuleType("oqp.utils")
    utils.__path__ = []
    sys.modules["oqp.utils"] = utils

    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")
    setattr(mpi_utils, "MPIManager", type("MPIManager", (), {}))
    setattr(mpi_utils, "MPIPool", type("MPIPool", (), {}))
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils

    matrix = types.ModuleType("oqp.utils.matrix")
    setattr(matrix, "DampingParam", type("DampingParam", (), {}))
    setattr(matrix, "DispersionModel", type("DispersionModel", (), {}))
    sys.modules["oqp.utils.matrix"] = matrix

    file_utils = types.ModuleType("oqp.utils.file_utils")
    setattr(file_utils, "dump_log", lambda *args, **kwargs: None)
    setattr(file_utils, "dump_data", lambda *args, **kwargs: None)
    setattr(file_utils, "write_config", lambda *args, **kwargs: None)
    setattr(file_utils, "write_xyz", lambda *args, **kwargs: None)
    sys.modules["oqp.utils.file_utils"] = file_utils

    state_labels = types.ModuleType("oqp.utils.state_labels")
    setattr(state_labels, "is_mrsf", lambda config: False)
    setattr(
        state_labels,
        "public_state_label",
        lambda config, state, **kwargs: f"state {int(state)}",
    )
    sys.modules["oqp.utils.state_labels"] = state_labels

    qmmm = types.ModuleType("oqp.utils.qmmm")
    sys.modules["oqp.utils.qmmm"] = qmmm

    # Stubbed like the rest so this file runs on its own.  Without it the
    # import only resolved when some earlier test in the same session had
    # already cached the real package, which made a single-file run fail.
    tb_backends = types.ModuleType("oqp.utils.tb_backends")
    setattr(tb_backends, "is_tb_method", lambda method: False)
    setattr(tb_backends, "make_tb_adapter", lambda mol: None)
    setattr(tb_backends, "tb_config", lambda *args, **kwargs: {})
    sys.modules["oqp.utils.tb_backends"] = tb_backends

    library = types.ModuleType("oqp.library")
    library.__path__ = []
    sys.modules["oqp.library"] = library

    state_tracking = types.ModuleType("oqp.library.state_tracking")
    setattr(state_tracking, "diagonal_phase_tracking", lambda *args, **kwargs: None)
    setattr(state_tracking, "maximum_overlap_assignment", lambda *args, **kwargs: None)
    sys.modules["oqp.library.state_tracking"] = state_tracking

    nac_utils = types.ModuleType("oqp.library.nac_utils")
    for name in (
        "canonical_state_overlap",
        "hst_derivative_coupling",
        "interstate_coupling",
        "load_numerical_nac_cache",
        "write_numerical_nac_cache_marker",
    ):
        setattr(nac_utils, name, lambda *args, **kwargs: None)
    sys.modules["oqp.library.nac_utils"] = nac_utils

    tb_backends = types.ModuleType("oqp.utils.tb_backends")
    setattr(tb_backends, "is_tb_method", lambda method: False)
    setattr(tb_backends, "make_tb_adapter", lambda mol: None)
    setattr(tb_backends, "tb_config", lambda config: {})
    sys.modules["oqp.utils.tb_backends"] = tb_backends

    frequency = types.ModuleType("oqp.library.frequency")
    setattr(frequency, "normal_mode", lambda *args, **kwargs: None)
    setattr(frequency, "thermal_analysis", lambda *args, **kwargs: None)
    sys.modules["oqp.library.frequency"] = frequency

    openqp_dftb = types.ModuleType("oqp.library.openqp_dftb")
    setattr(openqp_dftb, "OpenQPDFTBAdapter", type("OpenQPDFTBAdapter", (), {}))
    sys.modules["oqp.library.openqp_dftb"] = openqp_dftb


class FakeData(dict):
    def __init__(self):
        super().__init__()
        self.convergers = []
        self.diis_types = []
        self.sd_scf_flags = []

    def set_scf_converger_type(self, converger):
        self.convergers.append(converger)

    def set_scf_diis_type(self, diis_type):
        self.diis_types.append(diis_type)

    def set_sd_scf(self, flag):
        self.sd_scf_flags.append(flag)

    def set_trah_stability(self, flag):
        self.trah_stability = flag


class FakeMol:
    def __init__(self):
        self.config = {
            "scf": {"converger_type": "diis", "trh_stab": False},
        }
        self.data = FakeData()
        self.mol_energy = types.SimpleNamespace(energy=-1.0, SCF_converged=False)
        self.energies = []

    def write_molden(self, filename):
        raise AssertionError("save_molden is disabled in this test")


class SlottedMolEnergy:
    __slots__ = ("energy", "SCF_converged")

    def __init__(self, energy=-1.0, converged=False):
        self.energy = energy
        self.SCF_converged = converged


class TestSinglePointScfFallback(unittest.TestCase):
    def setUp(self):
        install_single_point_stubs()
        self.single_point = load_module(
            "single_point_under_test",
            "pyoqp/oqp/library/single_point.py",
        )

    def make_calculator(self):
        calc = self.single_point.SinglePoint.__new__(self.single_point.SinglePoint)
        calc.mol = FakeMol()
        calc.method = "hf"
        calc.td = "rpa"
        calc.init_scf = "no"
        calc.forced_attempt = 2
        calc.converger_type = "diis"
        calc.alternative_scf = "trah"
        calc.stability = False
        calc.save_molden = False
        calc.exception = False
        calc._prep_guess = lambda: None
        calc.swapmo = lambda: None
        calc.ixcore_shift = lambda: None
        calc.scf_calls = 0

        def scf():
            calc.scf_calls += 1
            calc.mol.mol_energy.energy = -1.0 - calc.scf_calls
            calc.mol.mol_energy.SCF_converged = calc.scf_calls == 2

        calc.scf = scf
        return calc

    def test_petite_state_requires_the_live_native_enable_flag(self):
        calc = self.make_calculator()
        calc.mol.symmetry_metadata = {
            'integral_symmetry': {'status': 'active'},
        }
        calc.mol.data['OQP::sym_petite_enable'] = [0]
        self.assertFalse(calc._petite_is_staged())

        calc.mol.data['OQP::sym_petite_enable'] = [1]
        self.assertTrue(calc._petite_is_staged())

    def test_energy_dispatches_native_caspt2_after_reference(self):
        caspt2 = types.ModuleType("oqp.library.caspt2_dyall")
        calls = []

        def fake_native_caspt2_energy(mol, ref_energy):
            calls.append((mol, ref_energy))
            return ["caspt2"]

        setattr(caspt2, "native_caspt2_energy", fake_native_caspt2_energy)
        sys.modules["oqp.library.caspt2_dyall"] = caspt2

        calc = self.make_calculator()
        calc.method = "caspt2"
        calc.reference = lambda *args, **kwargs: [-2.0]

        self.assertEqual(calc.energy(), ["caspt2"])
        self.assertEqual(calls, [(calc.mol, [-2.0])])

    def test_mp2_ignores_tdhf_ixcore_before_correlation(self):
        calc = self.make_calculator()
        calc.method = "mp2"
        calc.mol.config["tdhf"] = {"ixcore": "1"}
        calc.reference = lambda *args, **kwargs: [-2.0]
        calc.correlation = lambda ref_energy: [ref_energy[0] - 0.1]
        calc.ixcore_shift = lambda: self.fail("MP2 must not apply the TD ixcore shift")
        logs = []
        self.single_point.dump_log = lambda *args, **kwargs: logs.append(kwargs)

        self.assertEqual(calc.energy(), [-2.1])
        self.assertTrue(any("ignoring [tdhf] ixcore for mp2" in item["title"] for item in logs))

    def test_energy_dispatches_native_ms_caspt2_after_reference(self):
        caspt2 = types.ModuleType("oqp.library.caspt2_dyall")
        calls = []

        def fake_native_caspt2_energy(mol, ref_energy):
            calls.append((mol, ref_energy))
            return ["ms-caspt2"]

        setattr(caspt2, "native_caspt2_energy", fake_native_caspt2_energy)
        sys.modules["oqp.library.caspt2_dyall"] = caspt2

        calc = self.make_calculator()
        calc.method = "mscaspt2"
        calc.reference = lambda *args, **kwargs: [-2.0]

        self.assertEqual(calc.energy(), ["ms-caspt2"])
        self.assertEqual(calls, [(calc.mol, [-2.0])])

    def test_energy_dispatches_native_xms_caspt2_after_reference(self):
        caspt2 = types.ModuleType("oqp.library.caspt2_dyall")
        calls = []

        def fake_native_caspt2_energy(mol, ref_energy):
            calls.append((mol, ref_energy))
            return ["xms-caspt2"]

        setattr(caspt2, "native_caspt2_energy", fake_native_caspt2_energy)
        sys.modules["oqp.library.caspt2_dyall"] = caspt2

        calc = self.make_calculator()
        calc.method = "xmscaspt2"
        calc.reference = lambda *args, **kwargs: [-2.0]

        self.assertEqual(calc.energy(), ["xms-caspt2"])
        self.assertEqual(calls, [(calc.mol, [-2.0])])

    def test_energy_dispatches_sa_casscf_after_reference(self):
        casscf = types.ModuleType("oqp.library.casscf")
        calls = []

        class FakeCASSCF:
            def __init__(self, mol):
                calls.append(("init", mol))

            def energy(self, ref_energy):
                calls.append(("energy", ref_energy))
                return ["sa-casscf"]

        setattr(casscf, "CASSCF", FakeCASSCF)
        sys.modules["oqp.library.casscf"] = casscf

        calc = self.make_calculator()
        calc.method = "sacasscf"
        calc.reference = lambda *args, **kwargs: [-3.0]

        self.assertEqual(calc.energy(), ["sa-casscf"])
        self.assertEqual(calls, [("init", calc.mol), ("energy", [-3.0])])

    def test_energy_still_rejects_unknown_methods(self):
        calc = self.make_calculator()
        calc.method = "definitely-not-a-method"

        with self.assertRaisesRegex(ValueError, "Unknown method type"):
            calc.energy()

    def test_energy_restores_fallback_converger_after_successful_recovery(self):
        calc = self.make_calculator()

        energy = calc.energy(do_init_scf=False)

        self.assertEqual(energy, [-3.0])
        # Default escalation ladder: primary DIIS -> SOSCF (recovers here) -> ...,
        # then the primary converger is restored (in _run_scf and again in energy()).
        self.assertEqual(calc.mol.data.convergers, ["diis", "soscf", "diis", "diis"])
        self.assertEqual(calc.mol.data.sd_scf_flags, [False])

    def test_energy_restores_primary_converger_after_failed_recovery(self):
        calc = self.make_calculator()

        def scf_never_converges():
            calc.scf_calls += 1
            calc.mol.mol_energy.energy = -1.0 - calc.scf_calls
            calc.mol.mol_energy.SCF_converged = False

        calc.scf = scf_never_converges

        with self.assertRaises(RuntimeError):
            calc.energy(do_init_scf=False)

        # Full default ladder is walked (DIIS -> SOSCF -> TRAH) before giving up,
        # then the primary converger is restored.
        self.assertEqual(
            calc.mol.data.convergers, ["diis", "soscf", "trah", "diis", "diis"])

    def test_reference_restores_primary_converger_for_matching_gradient(self):
        calc = self.make_calculator()

        energy = calc.reference(do_init_scf=False)

        self.assertEqual(energy, [-3.0])
        self.assertEqual(calc.mol.data.convergers, ["diis", "soscf", "diis"])

    def test_next_geometry_restarts_from_diis_after_trah_recovery(self):
        calc = self.make_calculator()

        def scf_needs_trah():
            calc.scf_calls += 1
            # Each geometry requires DIIS, SOSCF, then TRAH.
            calc.mol.mol_energy.energy = -1.0 - calc.scf_calls
            calc.mol.mol_energy.SCF_converged = calc.scf_calls % 3 == 0

        calc.scf = scf_needs_trah

        calc.reference(do_init_scf=False)
        calc.reference(do_init_scf=False)

        self.assertEqual(
            calc.mol.data.convergers,
            [
                "diis", "soscf", "trah", "diis",
                "diis", "soscf", "trah", "diis",
            ],
        )

    def test_escalation_override_replaces_default_ladder(self):
        # scf.escalation overrides the default DIIS->SOSCF->TRAH chain with an
        # explicit comma-separated list (here: straight to TRAH, the old behavior).
        calc = self.make_calculator()
        calc.mol.config["scf"]["escalation"] = "trah"

        energy = calc.reference(do_init_scf=False)

        self.assertEqual(energy, [-3.0])
        self.assertEqual(calc.mol.data.convergers, ["diis", "trah", "diis"])

    def test_escalation_override_multi_stage_chain(self):
        # An explicit multi-stage chain is honored in order; the primary (diis)
        # is dropped from it. SOSCF recovers on the second SCF call.
        calc = self.make_calculator()
        calc.mol.config["scf"]["escalation"] = "diis,soscf,trah"

        energy = calc.reference(do_init_scf=False)

        self.assertEqual(energy, [-3.0])
        self.assertEqual(calc.mol.data.convergers, ["diis", "soscf", "diis"])

    def test_rstctmo_never_escalates_to_trah(self):
        calc = self.make_calculator()
        calc.mol.config["scf"]["rstctmo"] = True

        def scf_never_converges():
            calc.scf_calls += 1
            calc.mol.mol_energy.SCF_converged = False

        calc.scf = scf_never_converges

        self.assertFalse(calc._run_scf())
        self.assertEqual(calc.mol.data.convergers, ["diis", "soscf", "diis"])

    def test_rstctmo_disables_requested_trah_stability_pass(self):
        calc = self.make_calculator()
        calc.mol.config["scf"]["rstctmo"] = True
        calc.stability = True

        self.assertTrue(calc._run_scf())
        self.assertEqual(calc.mol.data.convergers, ["diis", "soscf", "diis"])

    def test_ml_selector_maps_diis_subtype_to_native_controls(self):
        # The distilled model labels DIIS variants, while native control uses a
        # top-level DIIS converger plus a DIIS subtype.
        model = types.ModuleType("oqp.library.scf_selector_model")
        setattr(model, "predict", lambda features: "ediis")
        sys.modules["oqp.library.scf_selector_model"] = model

        calc = self.make_calculator()
        calc.converger_type = "ml"
        calc.mol.config["scf"]["converger_type"] = "ml"

        converged = calc._run_scf()

        self.assertTrue(converged)
        self.assertEqual(calc.mol.data.diis_types, ["ediis"])
        self.assertEqual(calc.mol.data.convergers, ["diis", "soscf", "diis"])

    def test_stability_noop_restores_pre_trah_energy_metadata(self):
        calc = self.make_calculator()
        calc.stability = True

        def scf_stable_after_trah_check():
            calc.scf_calls += 1
            if calc.scf_calls == 1:
                calc.mol.mol_energy.energy = -2.0
                calc.mol.mol_energy.SCF_converged = True
            else:
                calc.mol.mol_energy.energy = -1.99999999
                calc.mol.mol_energy.SCF_converged = True

        calc.scf = scf_stable_after_trah_check

        energy = calc.reference(do_init_scf=False)

        self.assertEqual(energy, [-2.0])
        self.assertEqual(calc.mol.mol_energy.energy, -2.0)

    def test_stability_noop_restores_pre_trah_energy_metadata_without_dict(self):
        calc = self.make_calculator()
        calc.mol.mol_energy = SlottedMolEnergy()
        calc.stability = True

        def scf_stable_after_trah_check():
            calc.scf_calls += 1
            if calc.scf_calls == 1:
                calc.mol.mol_energy.energy = -2.0
                calc.mol.mol_energy.SCF_converged = True
            else:
                calc.mol.mol_energy.energy = -1.99999999
                calc.mol.mol_energy.SCF_converged = True

        calc.scf = scf_stable_after_trah_check

        energy = calc.reference(do_init_scf=False)

        self.assertEqual(energy, [-2.0])
        self.assertEqual(calc.mol.mol_energy.energy, -2.0)
        self.assertTrue(calc.mol.mol_energy.SCF_converged)

    def test_default_scf_stability_is_opt_in(self):
        # The post-SCF TRAH stability safeguard is expensive and can rotate
        # orbitals; it must not run unless the user writes [scf]
        # stability=true.
        oqpdata = (ROOT / "pyoqp/oqp/molecule/oqpdata.py").read_text()
        self.assertIn("'stability': {'type': bool, 'default': 'False'}", oqpdata)

    def test_stability_safeguard_relaxes_unstable_tdhf_reference(self):
        # The SCF stability safeguard can also run for the excited-state
        # reference SCF (method='tdhf', e.g. MRSF) when the user explicitly opts
        # in with [scf] stability=true.  A DIIS-converged but unstable
        # open-shell reference can otherwise make MRSF disagree with the
        # standalone SCF along a PES.  Before the fix this path was gated to
        # method=='hf', so no TRAH stability pass ran here and the unstable -2.0
        # solution would have been returned.
        calc = self.make_calculator()
        calc.method = "tdhf"
        calc.td = "mrsf"
        calc.stability = True

        def scf_unstable_then_relaxes():
            calc.scf_calls += 1
            if calc.scf_calls == 1:
                # primary DIIS converges to an unstable, higher solution
                calc.mol.mol_energy.energy = -2.0
                calc.mol.mol_energy.SCF_converged = True
            else:
                # the TRAH stability pass escapes to a genuinely lower solution
                calc.mol.mol_energy.energy = -2.5
                calc.mol.mol_energy.SCF_converged = True

        calc.scf = scf_unstable_then_relaxes

        converged = calc._run_scf()

        self.assertTrue(converged)
        # The safeguard ran for tdhf (a TRAH stability pass was invoked) ...
        self.assertIn("trah", calc.mol.data.convergers)
        # ... and the lower (relaxed) solution was kept.
        self.assertEqual(calc.mol.mol_energy.energy, -2.5)

    def test_stability_safeguard_skips_ordinary_tdhf_reference(self):
        # Ordinary closed-shell TDHF/TDA/RPA references should not run the extra
        # TRAH stability pass even when the user explicitly opts in for SCF
        # stability; the safeguard is only meaningful here for spin-flip/open-
        # shell references.
        calc = self.make_calculator()
        calc.method = "tdhf"
        calc.td = "rpa"
        calc.stability = True

        def scf_converges_once():
            calc.scf_calls += 1
            calc.mol.mol_energy.energy = -2.0
            calc.mol.mol_energy.SCF_converged = True

        calc.scf = scf_converges_once

        converged = calc._run_scf()

        self.assertTrue(converged)
        self.assertEqual(calc.mol.data.convergers, ["diis", "diis"])
        self.assertEqual(calc.mol.mol_energy.energy, -2.0)

    def test_stability_safeguard_relaxes_unstable_coupled_cluster_reference(self):
        # CCSD and CCSD(T) correlate the converged ground-state determinant, so
        # an unstable UHF/ROHF reference is as wrong for them as it is for the
        # HF energy itself.  The gate used to read method=='hf', which silently
        # dropped the safeguard the user had explicitly opted into and fed the
        # unstable -2.0 solution straight into the correlation step.
        for method in ("ccsd", "ccsd(t)"):
            with self.subTest(method=method):
                calc = self.make_calculator()
                calc.method = method
                calc.stability = True

                def scf_unstable_then_relaxes(calc=calc):
                    calc.scf_calls += 1
                    if calc.scf_calls == 1:
                        calc.mol.mol_energy.energy = -2.0
                        calc.mol.mol_energy.SCF_converged = True
                    else:
                        calc.mol.mol_energy.energy = -2.5
                        calc.mol.mol_energy.SCF_converged = True

                calc.scf = scf_unstable_then_relaxes

                converged = calc._run_scf()

                self.assertTrue(converged)
                self.assertIn("trah", calc.mol.data.convergers)
                self.assertEqual(calc.mol.mol_energy.energy, -2.5)


if __name__ == "__main__":
    unittest.main()
