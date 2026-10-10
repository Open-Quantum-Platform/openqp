"""Unit + parity tests for the unified GpuWorkspaceManager (Stage 1).

These tests are dependency-light (stdlib unittest, modules loaded by path) so they run
without CUDA or a built OpenQP shared library.
"""

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def ws():
    return load_module("gpu_workspace_under_test", "python/openqp_gpu/workspace_legacy.py")


def metc_buffers():
    return load_module("gpu_metc_buffers_for_ws", "python/openqp_gpu/gpu_metc_buffers.py")


def xc_cache():
    return load_module("tdhf_xc_response_cache_for_ws", "python/openqp_gpu/tdhf_xc_response_cache.py")


class ResidencyAndAccumulatorTests(unittest.TestCase):
    def test_all_four_residency_classes_exist(self):
        m = ws()
        self.assertEqual(
            {c.value for c in m.ResidencyClass},
            {"host_only", "device_resident", "host_device_mirrored", "borrowed_external"},
        )

    def test_per_thread_accumulator_must_be_device_resident(self):
        m = ws()
        with self.assertRaisesRegex(ValueError, "device-resident"):
            m.BufferSpec(
                "fock", 8, m.ResidencyClass.MIRRORED, "out",
                accumulator=m.AccumulatorMode.PER_THREAD, nthreads=4,
            )

    def test_per_thread_total_scales_by_nthreads_single_atomic_does_not(self):
        m = ws()
        per_thread = m.BufferSpec(
            "fock", 100, m.ResidencyClass.DEVICE_RESIDENT, "out",
            accumulator=m.AccumulatorMode.PER_THREAD, nthreads=8,
        )
        single = m.BufferSpec(
            "fock", 100, m.ResidencyClass.DEVICE_RESIDENT, "out",
            accumulator=m.AccumulatorMode.SINGLE_ATOMIC,
        )
        self.assertEqual(per_thread.total_bytes, 800)
        self.assertEqual(single.total_bytes, 100)

    def test_host_only_buffer_is_supported(self):
        m = ws()
        spec = m.BufferSpec("ids_repack", 64, m.ResidencyClass.HOST_ONLY, "staging")
        self.assertEqual(spec.residency, m.ResidencyClass.HOST_ONLY)
        self.assertEqual(spec.total_bytes, 64)


class ManagerCoreTests(unittest.TestCase):
    def make_schema(self, target="metc", shape=(1,), fock_bytes=100, nthreads=4):
        m = ws()
        buffers = (
            m.BufferSpec("ids", 16, m.ResidencyClass.MIRRORED, "eri_index"),
            m.BufferSpec(
                "fock", fock_bytes, m.ResidencyClass.DEVICE_RESIDENT, "out",
                accumulator=m.AccumulatorMode.PER_THREAD, nthreads=nthreads,
            ),
        )
        return m.WorkspaceSchema(target=target, shape_key=shape, buffers=buffers)

    def test_reuse_key_is_namespaced_by_target(self):
        schema = self.make_schema(target="metc", shape=(7, 5, 3, 11, 8))
        self.assertEqual(schema.reuse_key, ("metc", 7, 5, 3, 11, 8))

    def test_allocate_then_borrow_returns_pointer_equivalent_handle(self):
        m = ws()
        mgr = m.GpuWorkspaceManager()
        schema = self.make_schema(shape=(2,), fock_bytes=100, nthreads=4)
        record = mgr.allocate(schema)

        ids = mgr.borrow(schema.reuse_key, "ids")
        fock = mgr.borrow(schema.reuse_key, "fock")
        self.assertEqual((ids.offset, ids.total_bytes), (0, 16))
        # fock follows ids; per-thread expands to 4 * 100
        self.assertEqual((fock.offset, fock.total_bytes), (16, 400))
        self.assertEqual(fock.accumulator, m.AccumulatorMode.PER_THREAD)
        self.assertEqual(record.total_bytes, 16 + 400)

    def test_allocate_is_idempotent_for_identical_schema(self):
        m = ws()
        mgr = m.GpuWorkspaceManager()
        a = mgr.allocate(self.make_schema(shape=(3,)))
        b = mgr.allocate(self.make_schema(shape=(3,)))
        self.assertIs(a, b)

    def test_reuse_key_collision_with_different_layout_is_rejected(self):
        m = ws()
        mgr = m.GpuWorkspaceManager()
        mgr.allocate(self.make_schema(shape=(4,), fock_bytes=100))
        with self.assertRaisesRegex(ValueError, "different layout"):
            mgr.allocate(self.make_schema(shape=(4,), fock_bytes=999))

    def test_borrow_before_allocate_raises(self):
        m = ws()
        mgr = m.GpuWorkspaceManager()
        with self.assertRaisesRegex(KeyError, "not allocated"):
            mgr.borrow(("metc", 4), "ids")

    def test_borrow_unknown_buffer_raises(self):
        m = ws()
        mgr = m.GpuWorkspaceManager()
        schema = self.make_schema(shape=(5,))
        mgr.allocate(schema)
        with self.assertRaisesRegex(KeyError, "nope"):
            mgr.borrow(schema.reuse_key, "nope")

    def test_release_forgets_record(self):
        m = ws()
        mgr = m.GpuWorkspaceManager()
        schema = self.make_schema(shape=(6,))
        mgr.allocate(schema)
        self.assertIsNotNone(mgr.release(schema.reuse_key))
        self.assertIsNone(mgr.lookup(schema.reuse_key))
        self.assertIsNone(mgr.release(schema.reuse_key))

    def test_allocation_table_validator_catches_abi_drift(self):
        m = ws()
        schema = self.make_schema(shape=(7,))
        table = schema.allocation_table()
        self.assertEqual(schema.validate_table(table), table)
        bad = table[:1] + ((2, ("metc", "fock"), 12345, "out", "device_resident", "per_thread"),)
        with self.assertRaisesRegex(ValueError, "slot 2"):
            schema.validate_table(bad)

    def test_density_name_does_not_collide_across_targets(self):
        m = ws()
        metc = m.WorkspaceSchema(
            target="metc", shape_key=(1,),
            buffers=(m.BufferSpec("density", 8, m.ResidencyClass.MIRRORED, "input_matrix"),),
        )
        xc = m.WorkspaceSchema(
            target="xc_response", shape_key=(1,),
            buffers=(m.BufferSpec("density", 8, m.ResidencyClass.MIRRORED, "grid_density"),),
        )
        self.assertNotEqual(metc.namespaced_name("density"), xc.namespaced_name("density"))
        self.assertNotEqual(metc.reuse_key, xc.reuse_key)
        # both coexist in one manager
        mgr = m.GpuWorkspaceManager()
        mgr.allocate(metc)
        mgr.allocate(xc)
        metc_density = mgr.borrow(metc.reuse_key, "density")
        xc_density = mgr.borrow(xc.reuse_key, "density")
        self.assertEqual(metc_density.target, "metc")
        self.assertEqual(xc_density.target, "xc_response")
        self.assertNotEqual(metc_density, xc_density)

    def test_namespace_for_config_reads_target(self):
        m = ws()
        class FakeConfig:
            target = "XC_Response"
        self.assertEqual(m.namespace_for_config(FakeConfig()), "xc_response")


class MetcParityTests(unittest.TestCase):
    """Prove the legacy METC reuse key + byte sizes survive into the unified schema."""

    def make_plan(self):
        return metc_buffers().PersistentMetcBufferPlan.from_problem(
            nbf=7, nf=5, nmatrix=3, max_integrals=11, dtype_bytes=8,
        )

    def test_reuse_key_is_legacy_key_namespaced_by_metc(self):
        m = ws()
        plan = self.make_plan()
        schema = m.metc_schema_from_plan(plan, fock_accumulator=m.AccumulatorMode.SINGLE_ATOMIC)
        self.assertEqual(schema.shape_key, plan.reuse_key)
        self.assertEqual(schema.reuse_key, ("metc",) + plan.reuse_key)

    def test_base_bytes_match_legacy_bytes_for(self):
        m = ws()
        plan = self.make_plan()
        schema = m.metc_schema_from_plan(plan, fock_accumulator=m.AccumulatorMode.SINGLE_ATOMIC)
        for name in ("ids", "integrals", "density", "fock"):
            self.assertEqual(schema.buffer(name).base_bytes, plan.bytes_for(name), name)

    def test_single_atomic_total_equals_legacy_plan_total(self):
        m = ws()
        plan = self.make_plan()
        schema = m.metc_schema_from_plan(plan, fock_accumulator=m.AccumulatorMode.SINGLE_ATOMIC)
        self.assertEqual(schema.total_bytes, plan.total_bytes)

    def test_per_thread_fock_scales_total_by_nthreads(self):
        m = ws()
        plan = self.make_plan()
        schema = m.metc_schema_from_plan(
            plan, fock_accumulator=m.AccumulatorMode.PER_THREAD, nthreads=4,
        )
        # only fock replicates; the other three are single replicas
        expected = (
            plan.bytes_for("ids")
            + plan.bytes_for("integrals")
            + plan.bytes_for("density")
            + plan.bytes_for("fock") * 4
        )
        self.assertEqual(schema.total_bytes, expected)
        self.assertEqual(schema.buffer("fock").total_bytes, plan.bytes_for("fock") * 4)

    def test_fock_accumulator_policy_must_be_explicit(self):
        m = ws()
        plan = self.make_plan()
        # required keyword -> cannot silently assume a global accumulator
        with self.assertRaises(TypeError):
            m.metc_schema_from_plan(plan)
        with self.assertRaisesRegex(ValueError, "PER_THREAD or SINGLE_ATOMIC"):
            m.metc_schema_from_plan(plan, fock_accumulator=m.AccumulatorMode.NONE)

    def test_fock_is_device_resident_density_is_mirrored(self):
        m = ws()
        plan = self.make_plan()
        schema = m.metc_schema_from_plan(plan, fock_accumulator=m.AccumulatorMode.SINGLE_ATOMIC)
        self.assertEqual(schema.buffer("fock").residency, m.ResidencyClass.DEVICE_RESIDENT)
        self.assertEqual(schema.buffer("density").residency, m.ResidencyClass.MIRRORED)


class XcParityTests(unittest.TestCase):
    """Prove the legacy XC reuse key, byte sizes, and offsets survive unification."""

    def make_plan(self):
        return xc_cache().XcResponseCachePlan(
            nbf=10, ngrid=50, functional="b3lypv5", basis="3-21g",
            scf_type="rohf", response_type="rpa", spin_channels=2,
        )

    def test_reuse_key_is_legacy_key_namespaced_by_xc_response(self):
        m = ws()
        plan = self.make_plan()
        schema = m.xc_response_schema_from_plan(plan)
        self.assertEqual(schema.shape_key, plan.reuse_key())
        self.assertEqual(schema.reuse_key, ("xc_response",) + plan.reuse_key())

    def test_base_bytes_match_legacy_scalar_counts_times_dtype(self):
        m = ws()
        plan = self.make_plan()
        schema = m.xc_response_schema_from_plan(plan, dtype_bytes=8)
        self.assertEqual(schema.buffer("density").base_bytes, plan.density_values * 8)
        self.assertEqual(schema.buffer("potential").base_bytes, plan.potential_values * 8)
        self.assertEqual(schema.buffer("weights").base_bytes, plan.weight_values * 8)
        self.assertEqual(schema.buffer("ao_grid").base_bytes, plan.ao_grid_values * 8)
        self.assertEqual(schema.total_bytes, plan.total_workspace_bytes(8))

    def test_arena_offsets_match_legacy_layout_scaled_by_dtype(self):
        m = ws()
        plan = self.make_plan()
        schema = m.xc_response_schema_from_plan(plan, dtype_bytes=8)
        legacy = plan.workspace_layout()  # (name, scalar_offset, scalar_len)
        unified = schema.arena_layout()    # (name, byte_offset, byte_len)
        for (lname, loff, llen), (uname, uoff, ulen) in zip(legacy, unified):
            self.assertEqual(uname, lname)
            self.assertEqual(uoff, loff * 8)
            self.assertEqual(ulen, llen * 8)

    def test_borrowed_and_device_resident_classes_are_used(self):
        m = ws()
        plan = self.make_plan()
        schema = m.xc_response_schema_from_plan(plan)
        self.assertEqual(schema.buffer("ao_grid").residency, m.ResidencyClass.BORROWED)
        self.assertEqual(schema.buffer("weights").residency, m.ResidencyClass.BORROWED)
        self.assertEqual(schema.buffer("potential").residency, m.ResidencyClass.DEVICE_RESIDENT)
        self.assertEqual(schema.buffer("density").residency, m.ResidencyClass.MIRRORED)


if __name__ == "__main__":
    unittest.main()
