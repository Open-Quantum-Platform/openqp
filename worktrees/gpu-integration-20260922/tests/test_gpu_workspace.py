"""Tests for the Stage-1 unified GPU workspace manager.

These cover the manager's own contract (residency classes, the four C-ABI
verbs, target namespacing, the explicit f3 accumulator policy) plus parity
tests proving the legacy METC and XC-response schemas map cleanly into the
unified manager without colliding.

As with the legacy planning tests, modules are loaded by path so the suite runs
without importing the full OpenQP package or any CUDA/Fortran runtime.
"""

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


def load_workspace():
    return load_module("gpu_workspace_under_test", "python/openqp_gpu/gpu_workspace.py")


def load_metc_buffers():
    return load_module(
        "gpu_metc_buffers_ws_under_test", "python/openqp_gpu/gpu_metc_buffers.py"
    )


def load_xc_cache():
    return load_module(
        "tdhf_xc_response_cache_ws_under_test",
        "python/openqp_gpu/tdhf_xc_response_cache.py",
    )


class ResidencyAndVerbTests(unittest.TestCase):
    def test_residency_classes_present(self):
        ws = load_workspace()
        names = {member.name for member in ws.Residency}
        self.assertEqual(
            names, {"HOST_ONLY", "DEVICE_RESIDENT", "MIRRORED", "BORROWED"}
        )

    def test_valid_targets_are_metc_and_xc_response(self):
        ws = load_workspace()
        self.assertEqual(set(ws.VALID_TARGETS), {"metc", "xc_response"})

    def test_validate_table_normalizes_tuple_rows(self):
        ws = load_workspace()
        manager = ws.GpuWorkspaceManager()
        table = manager.validate_table(
            "metc",
            [
                ("ids", 176, "eri_index", ws.Residency.MIRRORED),
                ("fock", 4704, "output_matrix", ws.Residency.DEVICE_RESIDENT),
            ],
        )
        self.assertEqual(table[0].name, "ids")
        self.assertEqual(table[1].residency, ws.Residency.DEVICE_RESIDENT)

    def test_validate_table_rejects_unknown_target(self):
        ws = load_workspace()
        manager = ws.GpuWorkspaceManager()
        with self.assertRaisesRegex(ValueError, "unknown workspace target"):
            manager.validate_table("eri", [("x", 1, "r", ws.Residency.MIRRORED)])

    def test_validate_table_rejects_boolean_bytes(self):
        ws = load_workspace()
        manager = ws.GpuWorkspaceManager()
        with self.assertRaisesRegex(ValueError, "bytes must be a positive integer"):
            manager.validate_table(
                "metc", [("ids", True, "eri_index", ws.Residency.MIRRORED)]
            )

    def test_validate_table_rejects_non_residency(self):
        ws = load_workspace()
        manager = ws.GpuWorkspaceManager()
        with self.assertRaisesRegex(ValueError, "residency must be a Residency"):
            manager.validate_table("metc", [("ids", 16, "eri_index", "mirrored")])

    def test_validate_table_rejects_duplicate_names(self):
        ws = load_workspace()
        manager = ws.GpuWorkspaceManager()
        with self.assertRaisesRegex(ValueError, "duplicate workspace buffer name"):
            manager.validate_table(
                "metc",
                [
                    ("ids", 16, "eri_index", ws.Residency.MIRRORED),
                    ("ids", 16, "eri_index", ws.Residency.MIRRORED),
                ],
            )

    def test_allocate_lookup_and_total_bytes(self):
        ws = load_workspace()
        manager = ws.GpuWorkspaceManager()
        record = manager.allocate(
            "metc",
            (7, 5, 3, 11, 8),
            [
                ("ids", 176, "eri_index", ws.Residency.MIRRORED),
                ("fock", 4704, "output_matrix", ws.Residency.MIRRORED),
            ],
        )
        self.assertEqual(record.total_bytes, 176 + 4704)
        self.assertEqual(record.key.target, "metc")
        self.assertEqual(record.key.reuse_key, (7, 5, 3, 11, 8))
        self.assertEqual(manager.lookup("metc", (7, 5, 3, 11, 8)), record)

    def test_allocate_rejects_conflicting_reuse(self):
        ws = load_workspace()
        manager = ws.GpuWorkspaceManager()
        manager.allocate(
            "metc", (1, 1, 1, 1, 8), [("ids", 16, "eri_index", ws.Residency.MIRRORED)]
        )
        with self.assertRaisesRegex(ValueError, "reuse conflict"):
            manager.allocate(
                "metc",
                (1, 1, 1, 1, 8),
                [("ids", 32, "eri_index", ws.Residency.MIRRORED)],
            )

    def test_borrow_marks_buffers_borrowed_without_freeing_owner(self):
        ws = load_workspace()
        manager = ws.GpuWorkspaceManager()
        manager.allocate(
            "xc_response",
            ("bhhlyp", "6-31g*", "rhf", "rpa", 24, 1024, 1),
            [("density", 8192, "density", ws.Residency.DEVICE_RESIDENT)],
        )
        alias = manager.borrow("xc_response", ("bhhlyp", "6-31g*", "rhf", "rpa", 24, 1024, 1))
        self.assertTrue(alias.borrowed)
        self.assertTrue(all(b.residency is ws.Residency.BORROWED for b in alias.buffers))
        # Owner is still registered after a borrow.
        self.assertIsNotNone(
            manager.lookup("xc_response", ("bhhlyp", "6-31g*", "rhf", "rpa", 24, 1024, 1))
        )

    def test_borrow_unallocated_raises(self):
        ws = load_workspace()
        manager = ws.GpuWorkspaceManager()
        with self.assertRaises(KeyError):
            manager.borrow("metc", (9, 9, 9, 9, 8))

    def test_release_balances_borrow_then_frees_owner(self):
        ws = load_workspace()
        manager = ws.GpuWorkspaceManager()
        key = (5, 2, 3, 13, 8)
        manager.allocate(
            "metc", key, [("ids", 16, "eri_index", ws.Residency.MIRRORED)]
        )
        manager.borrow("metc", key)
        # First release returns the outstanding borrow; owner survives.
        released_borrow = manager.release("metc", key)
        self.assertTrue(released_borrow.borrowed)
        self.assertIsNotNone(manager.lookup("metc", key))
        # Second release (no borrows left) frees the owner.
        owner = manager.release("metc", key)
        self.assertFalse(owner.borrowed)
        self.assertIsNone(manager.lookup("metc", key))
        # Releasing again returns None.
        self.assertIsNone(manager.release("metc", key))


class F3AccumulatorPolicyTests(unittest.TestCase):
    def test_per_thread_is_the_default_policy(self):
        ws = load_workspace()
        buffers = load_metc_buffers()
        manager = ws.GpuWorkspaceManager()
        plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=4, nf=3, nmatrix=2, max_integrals=5
        )
        record = manager.from_metc_plan(plan, nthreads=8)
        self.assertEqual(record.f3_policy, ws.F3AccumulatorPolicy.PER_THREAD)

    def test_per_thread_replicates_the_fock_accumulator(self):
        ws = load_workspace()
        buffers = load_metc_buffers()
        manager = ws.GpuWorkspaceManager()
        plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=4, nf=3, nmatrix=2, max_integrals=5
        )
        single_fock = plan.bytes_for("fock")

        record = manager.from_metc_plan(
            plan, f3_policy=ws.F3AccumulatorPolicy.PER_THREAD, nthreads=8
        )
        fock = next(b for b in record.buffers if b.name == "fock")
        self.assertEqual(fock.bytes, single_fock * 8)

    def test_single_atomic_keeps_one_accumulator(self):
        ws = load_workspace()
        buffers = load_metc_buffers()
        manager = ws.GpuWorkspaceManager()
        plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=4, nf=3, nmatrix=2, max_integrals=5
        )
        record = manager.from_metc_plan(
            plan, f3_policy=ws.F3AccumulatorPolicy.SINGLE_ATOMIC, nthreads=8
        )
        fock = next(b for b in record.buffers if b.name == "fock")
        self.assertEqual(fock.bytes, plan.bytes_for("fock"))
        self.assertEqual(record.f3_policy, ws.F3AccumulatorPolicy.SINGLE_ATOMIC)


class MetcParityTests(unittest.TestCase):
    def test_metc_reuse_key_maps_into_namespaced_unified_key(self):
        ws = load_workspace()
        buffers = load_metc_buffers()
        manager = ws.GpuWorkspaceManager()
        plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=7, nf=5, nmatrix=3, max_integrals=11
        )
        record = manager.from_metc_plan(plan, nthreads=1)
        # Namespace prefix + legacy 5-int reuse key.
        self.assertEqual(record.key.target, "metc")
        self.assertEqual(record.key.reuse_key, plan.reuse_key)

    def test_metc_manifest_bytes_and_roles_survive_mapping(self):
        ws = load_workspace()
        buffers = load_metc_buffers()
        manager = ws.GpuWorkspaceManager()
        plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=7, nf=5, nmatrix=3, max_integrals=11
        )
        record = manager.from_metc_plan(plan, nthreads=1)
        mapped = {(b.name, b.bytes, b.role) for b in record.buffers}
        expected = {
            (item["name"], item["bytes"], item["role"])
            for item in plan.allocation_manifest()
        }
        # With nthreads=1 PER_THREAD the bytes are identical to the legacy manifest.
        self.assertEqual(mapped, expected)


class XcResponseParityTests(unittest.TestCase):
    def test_xc_reuse_key_maps_into_namespaced_unified_key(self):
        ws = load_workspace()
        cache = load_xc_cache()
        manager = ws.GpuWorkspaceManager()
        plan = cache.XcResponseCachePlan(
            nbf=24, ngrid=1024, functional="bhhlyp", basis="6-31g*",
            scf_type="rhf", response_type="rpa", spin_channels=1,
        )
        record = manager.from_xc_response_plan(plan)
        self.assertEqual(record.key.target, "xc_response")
        self.assertEqual(record.key.reuse_key, plan.reuse_key())

    def test_xc_scalar_offsets_survive_mapping(self):
        ws = load_workspace()
        cache = load_xc_cache()
        manager = ws.GpuWorkspaceManager()
        plan = cache.XcResponseCachePlan(
            nbf=10, ngrid=50, functional="b3lypv5", basis="3-21g",
            scf_type="rohf", response_type="rpa", spin_channels=2,
        )
        record = manager.from_xc_response_plan(plan)
        # The legacy scalar (name, offset, length) layout is preserved verbatim.
        self.assertEqual(record.scalar_layout, plan.workspace_layout())
        # And byte sizing applies the dtype width on top of the scalar lengths.
        density = next(b for b in record.buffers if b.name == "density")
        self.assertEqual(density.bytes, plan.density_values * 8)


class TargetNamespaceParityTests(unittest.TestCase):
    def test_metc_and_xc_density_do_not_collide(self):
        ws = load_workspace()
        buffers = load_metc_buffers()
        cache = load_xc_cache()
        manager = ws.GpuWorkspaceManager()

        metc_plan = buffers.PersistentMetcBufferPlan.from_problem(
            nbf=10, nf=2, nmatrix=1, max_integrals=4
        )
        xc_plan = cache.XcResponseCachePlan(
            nbf=10, ngrid=50, functional="b3lypv5", basis="3-21g",
            scf_type="rohf", response_type="rpa", spin_channels=2,
        )
        metc_record = manager.from_metc_plan(metc_plan, nthreads=1)
        xc_record = manager.from_xc_response_plan(xc_plan)

        # Both subsystems contribute a buffer named "density" ...
        metc_density = next(b for b in metc_record.buffers if b.name == "density")
        xc_density = next(b for b in xc_record.buffers if b.name == "density")
        self.assertNotEqual(metc_density.bytes, xc_density.bytes)

        # ... but the namespaced keys keep both records independently addressable.
        self.assertIsNotNone(manager.lookup("metc", metc_plan.reuse_key))
        self.assertIsNotNone(manager.lookup("xc_response", xc_plan.reuse_key()))
        self.assertNotEqual(metc_record.key, xc_record.key)
        self.assertEqual(len(manager._records), 2)


if __name__ == "__main__":
    unittest.main()
