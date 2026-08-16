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


if __name__ == "__main__":
    unittest.main()
