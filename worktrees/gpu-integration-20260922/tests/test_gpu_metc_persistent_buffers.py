import importlib.util
import sys
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


class GpuMetcPersistentBufferTests(unittest.TestCase):
    def test_buffer_plan_accounts_for_reusable_metc_device_arrays(self):
        buffers = load_module(
            "gpu_metc_buffers_under_test", "python/openqp_gpu/gpu_metc_buffers.py"
        )

        plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=7,
            nf=5,
            nmatrix=3,
            max_integrals=11,
            dtype_bytes=8,
        )

        self.assertEqual(plan.reuse_key, (7, 5, 3, 11, 8))
        self.assertEqual(plan.bytes_for("ids"), 11 * 4 * 4)
        self.assertEqual(plan.bytes_for("integrals"), 11 * 8)
        self.assertEqual(plan.bytes_for("density"), 3 * 5 * 7 * 7 * 8)
        self.assertEqual(plan.bytes_for("fock"), 3 * 5 * 7 * 7 * 8)
        self.assertEqual(
            plan.total_bytes,
            plan.bytes_for("ids")
            + plan.bytes_for("integrals")
            + plan.bytes_for("density")
            + plan.bytes_for("fock"),
        )

    def test_allocation_manifest_is_stable_for_fortran_cuda_wiring(self):
        buffers = load_module(
            "gpu_metc_buffers_manifest_under_test", "python/openqp_gpu/gpu_metc_buffers.py"
        )

        plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=4,
            nf=3,
            nmatrix=2,
            max_integrals=5,
        )

        self.assertEqual(
            plan.allocation_manifest(),
            (
                {"name": "ids", "bytes": 5 * 4 * 4, "role": "eri_index"},
                {"name": "integrals", "bytes": 5 * 8, "role": "eri_value"},
                {"name": "density", "bytes": 2 * 3 * 4 * 4 * 8, "role": "input_matrix"},
                {"name": "fock", "bytes": 2 * 3 * 4 * 4 * 8, "role": "output_matrix"},
            ),
        )

    def test_fortran_allocation_table_uses_stable_one_based_slots(self):
        buffers = load_module(
            "gpu_metc_buffers_fortran_table_under_test",
            "python/openqp_gpu/gpu_metc_buffers.py",
        )

        plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=4,
            nf=3,
            nmatrix=2,
            max_integrals=5,
        )

        self.assertEqual(
            plan.fortran_allocation_table(),
            (
                (1, "ids", 5 * 4 * 4, "eri_index"),
                (2, "integrals", 5 * 8, "eri_value"),
                (3, "density", 2 * 3 * 4 * 4 * 8, "input_matrix"),
                (4, "fock", 2 * 3 * 4 * 4 * 8, "output_matrix"),
            ),
        )

    def test_fortran_allocation_table_validator_catches_abi_drift(self):
        buffers = load_module(
            "gpu_metc_buffers_abi_validation_under_test",
            "python/openqp_gpu/gpu_metc_buffers.py",
        )

        plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=4,
            nf=3,
            nmatrix=2,
            max_integrals=5,
        )
        table = plan.fortran_allocation_table()

        self.assertEqual(plan.validate_fortran_allocation_table(table), table)
        with self.assertRaisesRegex(ValueError, "slot 3"):
            plan.validate_fortran_allocation_table(
                table[:2] + ((3, "scratch", 1, "temporary"),) + table[3:]
            )

    def test_runtime_stub_validates_fortran_table_before_registering_plan(self):
        buffers = load_module(
            "gpu_metc_buffers_runtime_plan_under_test",
            "python/openqp_gpu/gpu_metc_buffers.py",
        )
        runtime = load_module(
            "gpu_metc_persistent_runtime_under_test",
            "python/openqp_gpu/gpu_metc_persistent_runtime.py",
        )
        plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=6,
            nf=2,
            nmatrix=4,
            max_integrals=17,
        )
        stub = runtime.PersistentMetcAllocationRegistry()

        record = stub.register_plan(plan, plan.fortran_allocation_table())

        self.assertEqual(record.reuse_key, plan.reuse_key)
        self.assertEqual(record.total_bytes, plan.total_bytes)
        self.assertEqual(record.table, plan.fortran_allocation_table())
        self.assertEqual(stub.lookup(plan.reuse_key), record)
        with self.assertRaisesRegex(ValueError, "slot 2"):
            stub.register_plan(
                plan,
                ((1, "ids", 17 * 4 * 4, "eri_index"), (2, "bad", 1, "temporary")),
            )

    def test_runtime_stub_releases_records_by_reuse_key(self):
        buffers = load_module(
            "gpu_metc_buffers_runtime_release_plan_under_test",
            "python/openqp_gpu/gpu_metc_buffers.py",
        )
        runtime = load_module(
            "gpu_metc_persistent_runtime_release_under_test",
            "python/openqp_gpu/gpu_metc_persistent_runtime.py",
        )
        plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=5,
            nf=2,
            nmatrix=3,
            max_integrals=13,
        )
        registry = runtime.PersistentMetcAllocationRegistry()
        record = registry.register_plan(plan, plan.fortran_allocation_table())

        self.assertEqual(registry.release(plan.reuse_key), record)
        self.assertIsNone(registry.lookup(plan.reuse_key))
        self.assertIsNone(registry.release(plan.reuse_key))

    def test_buffer_plan_rejects_nonpositive_dimensions(self):
        buffers = load_module(
            "gpu_metc_buffers_validation_under_test", "python/openqp_gpu/gpu_metc_buffers.py"
        )

        with self.assertRaisesRegex(ValueError, "nbf"):
            buffers.PersistentMetcBufferPlan.from_problem(
                nbf=0,
                nf=5,
                nmatrix=3,
                max_integrals=11,
            )

    def test_buffer_plan_rejects_boolean_dimensions(self):
        buffers = load_module(
            "gpu_metc_buffers_bool_validation_under_test",
            "python/openqp_gpu/gpu_metc_buffers.py",
        )

        for field in ("nbf", "nf", "nmatrix", "max_integrals", "dtype_bytes"):
            kwargs = {
                "nbf": 4,
                "nf": 3,
                "nmatrix": 2,
                "max_integrals": 5,
                "dtype_bytes": 8,
            }
            kwargs[field] = True
            with self.subTest(field=field):
                with self.assertRaisesRegex(ValueError, field):
                    buffers.PersistentMetcBufferPlan.from_problem(**kwargs)

    def test_runtime_registry_rejects_boolean_reuse_keys(self):
        runtime = load_module(
            "gpu_metc_persistent_runtime_bool_key_under_test",
            "python/openqp_gpu/gpu_metc_persistent_runtime.py",
        )
        registry = runtime.PersistentMetcAllocationRegistry()

        with self.assertRaisesRegex(ValueError, "reuse_key"):
            registry.lookup((True, 3, 2, 5, 8))
        with self.assertRaisesRegex(ValueError, "reuse_key"):
            registry.release((4, False, 2, 5, 8))

    def test_runtime_registry_rejects_plan_with_boolean_reuse_key(self):
        runtime = load_module(
            "gpu_metc_persistent_runtime_bool_plan_key_under_test",
            "python/openqp_gpu/gpu_metc_persistent_runtime.py",
        )

        class BadPlan:
            reuse_key = (4, True, 2, 5, 8)
            total_bytes = 1

            def validate_fortran_allocation_table(self, table):
                return table

        registry = runtime.PersistentMetcAllocationRegistry()
        with self.assertRaisesRegex(ValueError, "reuse_key"):
            registry.register_plan(BadPlan(), ((1, "ids", 1, "eri_index"),))


if __name__ == "__main__":
    unittest.main()
