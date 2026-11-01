import importlib.util
import sys
import types
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_file_utils():
    stub_names = (
        "oqp",
        "oqp.molden",
        "oqp.molden.moldenwriter",
        "oqp.periodic_table",
        "oqp.runtime",
        "oqp.utils",
        "oqp.utils.constants",
        "oqp.utils.dftb_trace",
        "oqp.utils.log_format",
        "oqp.utils.mpi_utils",
        "oqp.utils.qmmm",
        "oqp.utils.state_labels",
    )
    saved_modules = {name: sys.modules.get(name) for name in stub_names}

    try:
        oqp = types.ModuleType("oqp")
        oqp.__path__ = []
        sys.modules["oqp"] = oqp

        molden = types.ModuleType("oqp.molden")
        molden.__path__ = []
        sys.modules["oqp.molden"] = molden

        moldenwriter = types.ModuleType("oqp.molden.moldenwriter")
        setattr(moldenwriter, "write_frequency", lambda *args, **kwargs: "")
        sys.modules["oqp.molden.moldenwriter"] = moldenwriter

        periodic_table = types.ModuleType("oqp.periodic_table")
        setattr(periodic_table, "SYMBOL_MAP", {1: 1, 8: 8, "1": 1, "8": 8})
        elements = [""] * 9
        elements[1] = "H"
        elements[8] = "O"
        setattr(periodic_table, "ELEMENTS_NAME", elements)
        sys.modules["oqp.periodic_table"] = periodic_table

        utils = types.ModuleType("oqp.utils")
        utils.__path__ = []
        sys.modules["oqp.utils"] = utils

        constants = types.ModuleType("oqp.utils.constants")
        setattr(constants, "ANGSTROM_TO_BOHR", 0.529177210903)
        sys.modules["oqp.utils.constants"] = constants

        runtime = types.ModuleType("oqp.runtime")
        setattr(runtime, "basis_search_paths", lambda: [])
        sys.modules["oqp.runtime"] = runtime

        dftb_trace = types.ModuleType("oqp.utils.dftb_trace")
        setattr(dftb_trace, "final_energy_annotation", lambda *args, **kwargs: "")
        setattr(dftb_trace, "final_energy_header", lambda *args, **kwargs: "")
        setattr(dftb_trace, "final_energy_label", lambda _summary, _root, label: label)
        sys.modules["oqp.utils.dftb_trace"] = dftb_trace

        log_format_spec = importlib.util.spec_from_file_location(
            "oqp.utils.log_format",
            ROOT / "pyoqp/oqp/utils/log_format.py",
        )
        assert log_format_spec is not None and log_format_spec.loader is not None
        log_format = importlib.util.module_from_spec(log_format_spec)
        sys.modules["oqp.utils.log_format"] = log_format
        log_format_spec.loader.exec_module(log_format)

        mpi_utils = types.ModuleType("oqp.utils.mpi_utils")
        setattr(mpi_utils, "mpi_dump", lambda func: func)
        sys.modules["oqp.utils.mpi_utils"] = mpi_utils

        qmmm = types.ModuleType("oqp.utils.qmmm")
        setattr(qmmm, "gradient_qmmm", lambda *args, **kwargs: None)
        sys.modules["oqp.utils.qmmm"] = qmmm

        state_labels = types.ModuleType("oqp.utils.state_labels")
        setattr(state_labels, "format_calculation_request", lambda *args, **kwargs: "")
        setattr(state_labels, "format_dftb_settings", lambda *args, **kwargs: "")
        setattr(state_labels, "is_mrsf", lambda *args, **kwargs: False)
        setattr(state_labels, "public_method_name", lambda *args, **kwargs: "HF")
        setattr(state_labels, "public_state_label", lambda _cfg, root, **kwargs: f"state {root}")
        setattr(state_labels, "spin_name", lambda mult: str(mult))
        sys.modules["oqp.utils.state_labels"] = state_labels

        spec = importlib.util.spec_from_file_location(
            "file_utils_under_test",
            ROOT / "pyoqp/oqp/utils/file_utils.py",
        )
        assert spec is not None
        assert spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        sys.modules["file_utils_under_test"] = module
        spec.loader.exec_module(module)
        return module
    finally:
        for name, module in saved_modules.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


def minimal_log_mol(log):
    native = types.SimpleNamespace(
        mol_prop=types.SimpleNamespace(natom=2, charge=0, mult=1),
        control=types.SimpleNamespace(scftype=1, maxit=50),
        tddft=types.SimpleNamespace(mult=1),
    )
    return types.SimpleNamespace(
        log=str(log),
        silent=True,
        config={
            "input": {
                "method": "hf", "basis": "sto-3g", "functional": "",
                "qmmm_flag": False, "runtype": "grad",
            },
            "scf": {
                "forced_attempt": 0, "conv": 1.0e-8,
                "incremental": False, "diis_type": "cdiis",
                "cdiis_switch": 0.0, "vdiis_vshift_switch": 0.0,
                "vshift": 0.0,
            },
            "tdhf": {
                "type": "tda", "maxit": 50, "maxit_zv": 50,
                "conv": 1.0e-6, "nstate": 1,
                "zvconv": 1.0e-8, "nvdav": 20,
            },
        },
        data=types.SimpleNamespace(_data=native),
        energies=[-1.0],
        grads=np.zeros((1, 2, 3)),
        get_atoms=lambda: np.array([1, 1]),
        get_system=lambda: np.zeros(6),
    )


class TestWriteXYZ(unittest.TestCase):
    def test_accepts_numpy_2d_atomic_number_arrays(self):
        file_utils = load_file_utils()
        atoms = np.array([[8], [1], [1]], dtype=np.int64)
        coord = np.array(
            [
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.8897261254576558,
                0.0,
                1.8897261254576558,
                0.0,
            ],
            dtype=float,
        )

        xyz = file_utils.write_xyz(atoms, coord, (3, -76.0))

        lines = xyz.splitlines()
        self.assertEqual(lines[0], "3")
        self.assertEqual(lines[1], "Geom 3 -76.0")
        self.assertTrue(lines[2].startswith("O"), lines[2])
        self.assertTrue(lines[3].startswith("H"), lines[3])
        self.assertTrue(lines[4].startswith("H"), lines[4])

    def test_representative_run_uses_ordered_sections_units_and_stable_markers(self):
        file_utils = load_file_utils()
        with self.subTest("render"):
            from tempfile import TemporaryDirectory
            import time

            with TemporaryDirectory() as folder:
                log = Path(folder) / "representative.log"
                mol = minimal_log_mol(log)
                mol.start_time = time.time() - 1.0

                file_utils.dump_log(mol, section="start", info={"build": "git HEAD abc123"})
                file_utils.print_module_banner(mol, "HF", "Hartree-Fock reference")
                file_utils.dump_log(mol, title="PyOQP: Calculation request",
                                    section="calculation")
                file_utils.dump_log(mol, title="PyOQP: Entering Gradient Calculation")
                file_utils.dump_log(mol, title="PyOQP: Normal SCF steps", section="scf")
                file_utils.dump_log(
                    mol, title="PyOQP: Final energy", section="energy",
                    info={"el": [-1.0], "d4": 0.0},
                )
                file_utils.dump_log(
                    mol, title="PyOQP: Final gradient", section="grad",
                    info={"el": np.zeros((1, 2, 3)), "d4": np.zeros((2, 3)),
                          "grad_list": [0]},
                )
                file_utils.dump_log(mol, section="end")

                text = log.read_text()

        headings = [
            "PyOQP LOG | RUN",
            "PyOQP LOG | INPUT AND REFERENCE",
            "PyOQP LOG | CALCULATION PROGRESS",
            "PyOQP LOG | CONVERGENCE AND ITERATIONS",
            "PyOQP LOG | ENERGIES AND STATES",
            "PyOQP LOG | GRADIENTS AND PROPERTIES",
            "PyOQP LOG | TIMING AND TERMINATION",
        ]
        positions = [text.index(heading) for heading in headings]
        self.assertEqual(positions, sorted(positions))
        self.assertIn("MODULE: HF", text)
        self.assertIn("PyOQP method:", text)
        self.assertIn("PyOQP Energy unit:                 Hartree", text)
        self.assertIn("PyOQP state 0", text)
        self.assertIn("-1.0000000000", text)
        self.assertIn("PyOQP Gradient unit:               Hartree/Bohr", text)
        self.assertIn("PyOQP electronic gradients", text)
        self.assertIn("PyOQP terminated at", text)

    def test_post_scf_correlation_never_returns_to_input_category(self):
        file_utils = load_file_utils()
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as folder:
            log = Path(folder) / "post-scf.log"
            mol = minimal_log_mol(log)
            file_utils.dump_log(
                mol, title="PyOQP: Calculation request", section="calculation")
            file_utils.dump_log(
                mol, title="PyOQP: Normal SCF steps", section="scf")
            for method in ("MP2", "CCSD", "CCSD(T)"):
                file_utils.dump_log(
                    mol, title=f"PyOQP: {method} correlation steps",
                    section="correlation",
                )
            text = log.read_text()

        input_heading = "PyOQP LOG | INPUT AND REFERENCE"
        convergence_heading = "PyOQP LOG | CONVERGENCE AND ITERATIONS"
        self.assertEqual(text.count(input_heading), 1)
        self.assertEqual(text.count(convergence_heading), 4)
        self.assertLess(text.index(input_heading), text.index(convergence_heading))
        self.assertNotIn(input_heading, text[text.index(convergence_heading):])
        for method in ("MP2", "CCSD", "CCSD(T)"):
            self.assertIn(f"PyOQP: {method} correlation steps", text)

    def test_dftd_optimization_settings_precede_method_stages(self):
        file_utils = load_file_utils()
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as folder:
            log = Path(folder) / "dftd-optimization.log"
            mol = minimal_log_mol(log)
            mol.config["input"]["functional"] = "pbe"
            file_utils.dump_log(
                mol, title="PyOQP: Calculation request", section="calculation")
            file_utils.dump_log(
                mol,
                title="PyOQP: Dispersion Correction",
                section="dftd",
                info={"type": "dftd4", "d4": True},
            )
            file_utils.dump_log(
                mol, title="PyOQP: Entering Geometry Optimization")
            file_utils.dump_log(
                mol, title="PyOQP: Normal SCF steps", section="scf")
            file_utils.dump_log(
                mol,
                title="PyOQP: Final Energy",
                section="energy",
                info={"el": [-1.0], "d4": 0.0},
            )
            text = log.read_text()

        categories = [
            line.split("| ", 1)[1]
            for line in text.splitlines()
            if "PyOQP LOG | " in line
        ]
        self.assertEqual(
            categories,
            [
                "INPUT AND REFERENCE",
                "INPUT AND REFERENCE",
                "CALCULATION PROGRESS",
                "CONVERGENCE AND ITERATIONS",
                "ENERGIES AND STATES",
            ],
        )
        self.assertIn("PyOQP dftd correction:", text)
        self.assertIn("PyOQP dftd method:", text)
        self.assertIn("PyOQP dftd functional:", text)

    def test_fci_results_use_energy_state_category_and_keep_markers(self):
        file_utils = load_file_utils()
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as folder:
            log = Path(folder) / "fci-results.log"
            mol = minimal_log_mol(log)
            file_utils.dump_log(
                mol,
                title="PyOQP: Complete Active Space Configuration Interaction",
                section="fci",
                info={
                    "method": "casci",
                    "ci_label": "CASCI",
                    "active_electrons": 2,
                    "active_orbitals": 2,
                    "frozen_core": 0,
                    "determinants": 4,
                    "orbital_source": "rhf",
                    "orbital_selection": "sequential",
                    "active_orbital_indices": "1,2",
                    "core_orbital_indices": "",
                    "energies": [-1.0, -0.8],
                    "s2": [0.0, 2.0],
                    "multiplicity": [1, 3],
                    "state_average": {
                        "energy": -0.9,
                        "weights": [0.5, 0.5],
                        "roots": [0, 1],
                        "root_indices": [0, 1],
                    },
                    "ci_vector_log": {
                        "threshold": 0.05,
                        "root_indices": (0,),
                        "entries": [
                            {
                                "index": 0,
                                "alpha": "1",
                                "beta": "1",
                                "coefficients": [0.95],
                            },
                        ],
                    },
                },
            )
            text = log.read_text()

        self.assertIn("PyOQP LOG | ENERGIES AND STATES", text)
        for marker in (
            "PyOQP CASCI energies",
            "<S^2>",
            "multiplicity",
            "PyOQP state-average roots:",
            "PyOQP state-average weights:",
            "PyOQP state-average energy:",
            "PyOQP CASCI CI vectors",
        ):
            self.assertIn(marker, text)

    def test_numerical_sections_preserve_category_heading_and_payload(self):
        file_utils = load_file_utils()
        from tempfile import TemporaryDirectory

        cases = (
            ("num_hess", [0, 6, 0.005, False, 6, 2, 4],
             "PyOQP hessian type"),
            ("hess_worker", [1, 2, "completed", (0.0, 2.0, 0, 4, "node0")],
             "displacement: 2"),
            ("num_nacv", [6, 0.005, False, 6, 2, 4],
             "PyOQP nac type"),
            ("nacv_worker", [1, 3, "completed", (0.0, 3.0, 0, 4, "node0")],
             "displacement: 3"),
            ("read_hess", None, "PyOQP read hessian file"),
        )
        with TemporaryDirectory() as folder:
            for section, info, payload in cases:
                with self.subTest(section=section):
                    log = Path(folder) / f"{section}.log"
                    mol = minimal_log_mol(log)
                    file_utils.dump_log(
                        mol, title=f"PyOQP: {section}", section=section, info=info)
                    text = log.read_text()
                    self.assertIn(
                        "PyOQP LOG | " + (
                            "GRADIENTS AND PROPERTIES"
                            if section == "read_hess"
                            else "CONVERGENCE AND ITERATIONS"
                        ),
                        text,
                    )
                    self.assertIn(payload, text)

    def test_qmmm_optimization_metrics_use_convergence_category(self):
        file_utils = load_file_utils()
        from tempfile import TemporaryDirectory

        metrics = {
            "istate": 0,
            "de": 0.001,
            "energy_shift": 0.0001,
            "rmsd_step": 0.01,
            "target_rmsd_step": 0.001,
            "max_step": 0.02,
            "target_max_step": 0.002,
            "rmsd_grad": 0.003,
            "target_rmsd_grad": 0.0003,
            "max_grad": 0.004,
            "target_max_grad": 0.0004,
        }
        with TemporaryDirectory() as folder:
            log = Path(folder) / "qmmm-opt.log"
            mol = minimal_log_mol(log)
            file_utils.dump_log(
                mol,
                title="Geometry Optimization Convergence 1",
                section="QM/MM",
                info=metrics,
            )
            text = log.read_text()

        self.assertIn("PyOQP LOG | CONVERGENCE AND ITERATIONS", text)
        self.assertIn("PyOQP energy shift:", text)
        self.assertIn("PyOQP rmsd grad:", text)


class TestMECPStateLabels(unittest.TestCase):
    """MECP joins two multiplicity solves into one energy array.

    The final summary must label each block in its own multiplicity.  Before
    the fix the triplet roots of ``mecp(S0,T0)`` with nstate=5 printed as
    S5..S9 and the T0 gradient as S5.
    """

    def _render(self, info_energy, info_grad):
        from tempfile import TemporaryDirectory

        file_utils = load_file_utils()
        # use the real labelling rules rather than the stub; state_labels
        # imports only the pure log_format sibling, provided here so the
        # compiled package is never imported
        names = ("oqp", "oqp.utils", "oqp.utils.log_format")
        saved = {name: sys.modules.get(name) for name in names}
        try:
            for name in names[:2]:
                package = types.ModuleType(name)
                package.__path__ = []
                sys.modules[name] = package
            format_spec = importlib.util.spec_from_file_location(
                "oqp.utils.log_format", ROOT / "pyoqp/oqp/utils/log_format.py")
            log_format = importlib.util.module_from_spec(format_spec)
            sys.modules["oqp.utils.log_format"] = log_format
            format_spec.loader.exec_module(log_format)
            spec = importlib.util.spec_from_file_location(
                "state_labels_for_mecp_test",
                ROOT / "pyoqp/oqp/utils/state_labels.py",
            )
            labels = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(labels)
        finally:
            for name, module in saved.items():
                if module is None:
                    sys.modules.pop(name, None)
                else:
                    sys.modules[name] = module
        file_utils.is_mrsf = labels.is_mrsf
        file_utils.public_state_label = labels.public_state_label

        with TemporaryDirectory() as folder:
            log = Path(folder) / "mecp.log"
            mol = minimal_log_mol(log)
            mol.config["input"].update(method="tdhf", runtype="mecp")
            mol.config["scf"].update(type="rohf", multiplicity=3)
            mol.config["tdhf"].update(type="mrsf", multiplicity=1, nstate=5)
            energies = info_energy["el"]
            mol.energies = list(energies)
            mol.grads = info_grad["el"]
            file_utils.dump_log(mol, title="PyOQP: Final Energy",
                                section="energy", info=info_energy)
            file_utils.dump_log(mol, title="PyOQP: Final Gradient",
                                section="grad", info=info_grad)
            return log.read_text()

    def test_triplet_block_is_labelled_t0_to_t4(self):
        singlets = [-78.4066, -78.3101, -78.3081, -78.1360, -78.1345]
        triplets = [-78.4066, -78.1515, -78.1472, -78.1470, -78.1417]
        energies = [-78.4318] + singlets + triplets
        blocks = ((1, 5), (3, 5))
        grads = np.arange(11 * 2 * 3, dtype=float).reshape(11, 2, 3)
        text = self._render(
            {"el": energies, "d4": 0.0, "state_blocks": blocks},
            {"el": grads, "d4": np.zeros((2, 3)), "grad_list": [1, 6],
             "state_blocks": blocks},
        )
        energy_block = text.split("ENERGIES AND STATES", 1)[1].split(
            "GRADIENTS AND PROPERTIES", 1)[0]
        labels = [line.split()[1] for line in energy_block.splitlines()
                  if line.strip().startswith("PyOQP ")
                  and line.split()[1][:1] in "ST"
                  and line.split()[1][1:].isdigit()]
        self.assertEqual(
            labels, ["S0", "S1", "S2", "S3", "S4", "T0", "T1", "T2", "T3", "T4"])
        self.assertIn("triplet ROHF reference (internal)", energy_block)
        gradient_block = text.split("GRADIENTS AND PROPERTIES", 1)[1]
        self.assertIn("PyOQP T0", gradient_block)
        self.assertNotIn("PyOQP S5", text)

    def test_clipped_first_block_moves_the_triplet_start(self):
        # nstate=5 but only 3 singlet roots were solved
        energies = [-1.0, -0.9, -0.8, -0.7, -0.95, -0.85]
        blocks = ((1, 3), (3, 2))
        grads = np.zeros((6, 2, 3))
        text = self._render(
            {"el": energies, "d4": 0.0, "state_blocks": blocks},
            {"el": grads, "d4": np.zeros((2, 3)), "grad_list": [1, 4],
             "state_blocks": blocks},
        )
        self.assertIn("PyOQP S2 ", text)
        self.assertIn("PyOQP T0 ", text)
        self.assertIn("PyOQP T1 ", text)
        self.assertNotIn("PyOQP S3", text)

    def test_index_outside_the_blocks_is_an_error(self):
        file_utils = load_file_utils()
        with self.assertRaisesRegex(IndexError, "outside the 4 roots"):
            file_utils.stacked_state_label({}, 5, ((1, 2), (3, 2)))


class TestVerbosityInPythonLog(unittest.TestCase):
    def _energy_log(self, d4_flag, verbose=None):
        file_utils = load_file_utils()
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as folder:
            log = Path(folder) / "verbosity.log"
            mol = minimal_log_mol(log)
            mol.config["input"]["d4"] = d4_flag
            if verbose is not None:
                mol.config["input"]["verbose"] = verbose
            file_utils.dump_log(
                mol,
                title="PyOQP: Dispersion Correction",
                section="dftd",
                info={"type": "dftd4", "d4": d4_flag},
            )
            file_utils.dump_log(
                mol,
                title="PyOQP: Final Energy",
                section="energy",
                info={"el": [-1.0], "d4": 0.0},
            )
            return log.read_text()

    def test_dispersion_records_only_when_dispersion_is_requested(self):
        off = self._energy_log(False)
        self.assertNotIn("PyOQP: Dispersion Correction", off)
        self.assertNotIn("dispersion corrected energies", off)
        self.assertIn("PyOQP electronic energies", off)
        on = self._energy_log(True)
        self.assertEqual(on.split("PyOQP dftd correction:", 1)[1].split()[0], "yes")
        self.assertIn("PyOQP dispersion corrected energies", on)

    def test_detailed_log_keeps_the_disabled_dispersion_block(self):
        text = self._energy_log(False, verbose=2)
        self.assertIn("PyOQP: Dispersion Correction", text)
        self.assertEqual(text.split("PyOQP dftd correction:", 1)[1].split()[0], "no")

    def test_python_log_lines_carry_no_trailing_whitespace(self):
        text = self._energy_log(True)
        self.assertEqual([line for line in text.splitlines() if line != line.rstrip()], [])


class TestMECPStatusTable(unittest.TestCase):
    """opt_status.txt keeps the mean energy and the search objective apart."""

    def test_energy_column_is_the_mean_and_objective_has_its_own(self):
        from tempfile import TemporaryDirectory

        file_utils = load_file_utils()
        atoms = np.array([[8], [1]])
        coordinates = np.zeros(6)
        # step 4 of the twisted-ethylene SQP MECP: mean -78.40055722,
        # Lagrangian -78.52136313 at multiplier -53
        rows = [
            (1, atoms, coordinates, -78.39592156, -78.39592156, -78.39592156,
             -0.00038109, 0.0, 0.0, 0.034120, 0.102155, 0.0, 0.0),
            (4, atoms, coordinates, -78.40055722, -78.52136313, -0.12303000,
             0.00228043, 0.038897, 0.080952, 0.355466, 0.734927,
             0.353892, 0.730535),
        ]
        with TemporaryDirectory() as folder:
            for row in rows:
                file_utils.dump_data(types.SimpleNamespace(), row,
                                     title="MECP", fpath=folder)
            lines = (Path(folder) / "opt_status.txt").read_text().splitlines()
            comment = (Path(folder) / "opt.xyz").read_text().splitlines()[1]

        header = lines[0].split()
        self.assertEqual(header[:5],
                         ["Step", "Energy", "Objective", "Obj", "Shift"])
        step4 = lines[-1].split()
        self.assertEqual(step4[0], "4")
        self.assertAlmostEqual(float(step4[1]), -78.40055722, places=8)
        self.assertAlmostEqual(float(step4[2]), -78.52136313, places=8)
        self.assertAlmostEqual(float(step4[3]), -0.12303000, places=8)
        self.assertAlmostEqual(float(step4[4]), 0.00228043, places=8)
        self.assertEqual(len(step4), 11)
        # the geometry comment carries the energy, not the objective
        self.assertIn("-78.40055722", comment)

    def test_convergence_block_names_the_objective(self):
        from tempfile import TemporaryDirectory

        file_utils = load_file_utils()
        metrics = {
            "istate": 1, "jstate": 1, "nstate": 5, "mecp_search": "sqp",
            "energy": -78.40055722, "objective": -78.52136313,
            "de": -0.12303, "energy_shift": 1.0e-6,
            "gap": 0.00228043, "energy_gap": 1.0e-5,
            "rmsd_step": 0.038897, "target_rmsd_step": 1.0e-3,
            "max_step": 0.080952, "target_max_step": 2.0e-3,
            "rmsd_grad": 0.355466, "target_rmsd_grad": 1.0e-4,
            "max_grad": 0.734927, "target_max_grad": 3.0e-4,
        }
        with TemporaryDirectory() as folder:
            log = Path(folder) / "mecp-convergence.log"
            mol = minimal_log_mol(log)
            mol.config["optimize"] = {"imult": 1, "jmult": 3}
            file_utils.dump_log(mol, title="Geometry Optimization Convergence 4",
                                section="mecp", info=metrics)
            text = log.read_text()

        self.assertIn("PyOQP mean energy (Hartree):          -78.40055722", text)
        self.assertIn("PyOQP search objective (Hartree):     -78.52136313", text)
        self.assertIn("PyOQP objective shift:", text)
        self.assertNotIn("PyOQP energy shift:", text)

    def test_convergence_block_without_the_energies_is_an_error(self):
        # record() sets both before check_convergence logs; a caller that
        # skips it must fail rather than print nan
        from tempfile import TemporaryDirectory

        file_utils = load_file_utils()
        metrics = {
            "istate": 1, "jstate": 1, "nstate": 5, "mecp_search": "sqp",
            "objective": -78.5, "de": -0.1, "energy_shift": 1.0e-6,
            "gap": 0.002, "energy_gap": 1.0e-5,
            "rmsd_step": 0.03, "target_rmsd_step": 1.0e-3,
            "max_step": 0.08, "target_max_step": 2.0e-3,
            "rmsd_grad": 0.3, "target_rmsd_grad": 1.0e-4,
            "max_grad": 0.7, "target_max_grad": 3.0e-4,
        }
        with TemporaryDirectory() as folder:
            mol = minimal_log_mol(Path(folder) / "mecp.log")
            mol.config["optimize"] = {"imult": 1, "jmult": 3}
            with self.assertRaisesRegex(KeyError, "energy"):
                file_utils.dump_log(mol, title="Geometry Optimization Convergence 1",
                                    section="mecp", info=metrics)


if __name__ == "__main__":
    unittest.main()
