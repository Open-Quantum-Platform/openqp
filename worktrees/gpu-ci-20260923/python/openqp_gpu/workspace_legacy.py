"""Unified pure-Python GPU workspace manager for OpenQP experimental backends.

Stage 1 of the GPU HF/DFT roadmap (see ``docs/gpu_workspace_design.md``): collapse
the two independently-grown planning registries

* ``gpu_metc_buffers.PersistentMetcBufferPlan`` / ``PersistentMetcAllocationRegistry``
  (MRSF/UMRSF METC contraction buffers), and
* ``tdhf_xc_response_cache.XcResponseCachePlan``
  (TDHF/TDDFT XC-response cache workspace)

into one :class:`GpuWorkspaceManager`.

The module is deliberately pure Python and free of CUDA / built-OpenQP imports so the
ownership, namespacing, residency, and reuse-key contracts can be tested at source level
before any C-ABI or kernel wiring exists.  The verbs (:meth:`allocate`,
:meth:`validate_table`, :meth:`borrow`, :meth:`release`) are the exact surface a later
``oqp_gpu_ws_*`` C ABI is expected to export.

NOTE: no CUDA buffers are touched here.  A "device pointer" is modelled as a
:class:`BufferHandle` carrying a byte offset into a single notional workspace arena --
exactly what ``base_ptr + offset`` would resolve to once real allocation is wired.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


# --------------------------------------------------------------------------------------
# Namespace authority (mirrors ``GpuConfig.target`` from oqp.utils.gpu on the
# gpu-xc-response branch).  Buffers are namespaced by target because METC and XC both use
# names like ``density`` with different meanings; the namespace keeps them distinct in one
# manager.  Kept as a duck-typed string so the core does not hard-depend on the xc-response
# scaffold file; ``namespace_for_config`` extracts it from any GpuConfig-like object.
# --------------------------------------------------------------------------------------

class GpuTarget(str, Enum):
    """Known GPU workspace namespaces (one per ``GpuConfig.target``)."""

    METC = "metc"
    XC_RESPONSE = "xc_response"


def namespace_for_config(gpu_config: Any) -> str:
    """Return the workspace namespace for a ``GpuConfig``-like object.

    Duck-typed on a ``.target`` attribute so source-level tests (and the eventual
    ``oqp.utils.gpu.GpuConfig``) can drive namespacing without importing that module.
    """

    return str(gpu_config.target).lower()


# --------------------------------------------------------------------------------------
# Residency + accumulator taxonomy
# --------------------------------------------------------------------------------------

class ResidencyClass(str, Enum):
    """Where a workspace buffer physically lives and who owns its memory."""

    HOST_ONLY = "host_only"            # host RAM only; never uploaded (e.g. index repack staging)
    DEVICE_RESIDENT = "device_resident"  # device only; host copy not maintained during compute
    MIRRORED = "host_device_mirrored"  # exists on both host and device, kept in sync
    BORROWED = "borrowed_external"     # reference to memory the manager neither allocates nor frees


class AccumulatorMode(str, Enum):
    """How a device-resident accumulator handles OpenMP thread concurrency.

    METC's Fortran ``this%f3`` accumulator is *per OpenMP thread*
    (``tdhf_mrsf_lib.F90``).  A device-resident accumulator must therefore declare its
    policy explicitly -- the manager never silently assumes a single global accumulator.
    """

    NONE = "none"                  # not an accumulator
    PER_THREAD = "per_thread"      # nthreads independent replicas; total = nthreads * base
    SINGLE_ATOMIC = "single_atomic"  # one device buffer, atomic cross-thread accumulation


# --------------------------------------------------------------------------------------
# Buffer + schema
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class BufferSpec:
    """One named device/host buffer within a workspace schema.

    ``base_bytes`` is the size of a single replica.  Per-thread accumulators expand to
    ``base_bytes * nthreads`` (see :meth:`total_bytes`); every other buffer is a single
    replica.
    """

    name: str
    base_bytes: int
    residency: ResidencyClass
    role: str
    accumulator: AccumulatorMode = AccumulatorMode.NONE
    nthreads: int = 1

    def __post_init__(self) -> None:
        if self.base_bytes <= 0:
            raise ValueError(f"{self.name}: base_bytes must be positive, got {self.base_bytes}")
        if self.nthreads <= 0:
            raise ValueError(f"{self.name}: nthreads must be positive, got {self.nthreads}")
        if self.accumulator is AccumulatorMode.PER_THREAD and self.residency is not ResidencyClass.DEVICE_RESIDENT:
            raise ValueError(f"{self.name}: per-thread accumulators must be device-resident")
        if self.accumulator is AccumulatorMode.NONE and self.nthreads != 1:
            raise ValueError(f"{self.name}: nthreads is only meaningful for accumulators")

    @property
    def total_bytes(self) -> int:
        """Total device bytes including per-thread replication."""

        if self.accumulator is AccumulatorMode.PER_THREAD:
            return self.base_bytes * self.nthreads
        return self.base_bytes


# C-ABI export row: (slot, namespaced_name, total_bytes, role, residency, accumulator)
AllocationRow = tuple
AllocationTable = tuple


@dataclass(frozen=True)
class WorkspaceSchema:
    """A namespaced, ordered set of buffers identified by a shape reuse key.

    ``shape_key`` is the legacy planner's ``reuse_key`` (the problem-shape identity);
    :attr:`reuse_key` prepends ``target`` so two targets that share an inner key (and a
    buffer name like ``density``) never collide in one manager.
    """

    target: str
    shape_key: tuple
    buffers: tuple[BufferSpec, ...]

    def __post_init__(self) -> None:
        if not self.buffers:
            raise ValueError("schema must declare at least one buffer")
        names = [b.name for b in self.buffers]
        if len(names) != len(set(names)):
            raise ValueError(f"duplicate buffer names within target {self.target!r}: {names}")

    @property
    def reuse_key(self) -> tuple:
        """Namespaced reuse key: ``(target, *shape_key)``."""

        return (self.target, *self.shape_key)

    def namespaced_name(self, name: str) -> tuple[str, str]:
        """Return ``(target, name)`` -- the collision-proof identity of a buffer."""

        return (self.target, name)

    def buffer(self, name: str) -> BufferSpec:
        for b in self.buffers:
            if b.name == name:
                return b
        raise KeyError(f"no buffer {name!r} in target {self.target!r}")

    @property
    def total_bytes(self) -> int:
        return sum(b.total_bytes for b in self.buffers)

    def allocation_table(self) -> AllocationTable:
        """Stable one-based C-ABI export table.

        Each row is ``(slot, namespaced_name, total_bytes, role, residency, accumulator)``.
        This is the artifact a future Fortran/C ABI validates against, analogous to the
        legacy ``fortran_allocation_table``.
        """

        return tuple(
            (
                slot,
                self.namespaced_name(b.name),
                b.total_bytes,
                b.role,
                b.residency.value,
                b.accumulator.value,
            )
            for slot, b in enumerate(self.buffers, start=1)
        )

    def validate_table(self, table: AllocationTable) -> AllocationTable:
        """Validate an ABI table against this schema, pinpointing the first drift."""

        expected = self.allocation_table()
        if tuple(table) != expected:
            for exp_row, act_row in zip(expected, tuple(table)):
                if act_row != exp_row:
                    raise ValueError(
                        f"allocation table mismatch at slot {exp_row[0]}: "
                        f"expected {exp_row!r}, got {act_row!r}"
                    )
            raise ValueError(
                f"allocation table length mismatch: expected {len(expected)}, got {len(tuple(table))}"
            )
        return expected

    def arena_layout(self) -> tuple[tuple[str, int, int], ...]:
        """Contiguous byte layout as ``(name, offset, total_bytes)`` in declaration order.

        Offsets are the future ``base_ptr + offset`` for each buffer.
        """

        layout = []
        offset = 0
        for b in self.buffers:
            layout.append((b.name, offset, b.total_bytes))
            offset += b.total_bytes
        return tuple(layout)


# --------------------------------------------------------------------------------------
# Allocation record + borrowed handle
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class BufferHandle:
    """A borrowed reference to one buffer -- the stand-in for a device pointer.

    Carries everything an ``oqp_gpu_ws_ptr`` call would need to resolve a real pointer:
    namespace, name, ABI slot, byte offset into the arena, size, and residency/accumulator.
    """

    target: str
    name: str
    slot: int
    offset: int
    total_bytes: int
    residency: ResidencyClass
    accumulator: AccumulatorMode


@dataclass(frozen=True)
class WorkspaceRecord:
    """Validated metadata for one allocated workspace shape."""

    reuse_key: tuple
    target: str
    total_bytes: int
    table: AllocationTable
    _handles: dict[str, BufferHandle] = field(default_factory=dict, repr=False)

    def handle(self, name: str) -> BufferHandle:
        try:
            return self._handles[name]
        except KeyError:
            raise KeyError(f"no buffer {name!r} in target {self.target!r}") from None


# --------------------------------------------------------------------------------------
# The unified manager
# --------------------------------------------------------------------------------------

class GpuWorkspaceManager:
    """Single owner of all (notional) GPU workspace allocations.

    Verbs map one-to-one onto the planned C ABI:

    * :meth:`allocate` -> ``oqp_gpu_ws_acquire`` (alloc-or-reuse by reuse key, validates table)
    * :meth:`validate_table` -> ``oqp_gpu_ws_validate`` (ABI drift guard for export)
    * :meth:`borrow` -> ``oqp_gpu_ws_ptr`` (borrow a device pointer; never transfers ownership)
    * :meth:`release` -> ``oqp_gpu_ws_release``

    No CUDA memory is allocated; records are pure metadata.
    """

    def __init__(self) -> None:
        self._records: dict[tuple, WorkspaceRecord] = {}

    def allocate(self, schema: WorkspaceSchema) -> WorkspaceRecord:
        """Allocate-or-reuse the workspace for ``schema`` by its namespaced reuse key.

        Idempotent: re-allocating an identical schema returns the existing record.  A
        reuse-key collision with a *different* table is rejected (shape drift).
        """

        key = schema.reuse_key
        existing = self._records.get(key)
        table = schema.validate_table(schema.allocation_table())
        if existing is not None:
            if existing.table != table:
                raise ValueError(
                    f"reuse key {key!r} already allocated with a different layout"
                )
            return existing

        handles: dict[str, BufferHandle] = {}
        for (slot, (target, name), total_bytes, role, residency, accumulator), (
            _, offset, _
        ) in zip(table, schema.arena_layout()):
            handles[name] = BufferHandle(
                target=target,
                name=name,
                slot=slot,
                offset=offset,
                total_bytes=total_bytes,
                residency=ResidencyClass(residency),
                accumulator=AccumulatorMode(accumulator),
            )
        record = WorkspaceRecord(
            reuse_key=key,
            target=schema.target,
            total_bytes=schema.total_bytes,
            table=table,
            _handles=handles,
        )
        self._records[key] = record
        return record

    def validate_table(self, schema: WorkspaceSchema, table: AllocationTable) -> AllocationTable:
        """Validate an externally-supplied ABI table against ``schema`` (C-ABI export guard)."""

        return schema.validate_table(table)

    def lookup(self, reuse_key: tuple) -> WorkspaceRecord | None:
        return self._records.get(reuse_key)

    def borrow(self, reuse_key: tuple, name: str) -> BufferHandle:
        """Borrow a pointer-equivalent handle to one buffer of an allocated workspace."""

        record = self._records.get(reuse_key)
        if record is None:
            raise KeyError(f"workspace {reuse_key!r} is not allocated")
        return record.handle(name)

    def release(self, reuse_key: tuple) -> WorkspaceRecord | None:
        """Forget and return the record for ``reuse_key`` if one is allocated."""

        return self._records.pop(reuse_key, None)


# --------------------------------------------------------------------------------------
# Adapters: legacy plans -> unified schema (these are the parity bridge)
# --------------------------------------------------------------------------------------

def metc_schema_from_plan(
    plan: Any,
    *,
    fock_accumulator: AccumulatorMode,
    nthreads: int = 1,
    target: str = GpuTarget.METC.value,
) -> WorkspaceSchema:
    """Map a ``PersistentMetcBufferPlan`` onto the unified schema.

    ``fock_accumulator`` is **required** (no default): the per-thread ``this%f3`` policy
    must be chosen explicitly -- pass :attr:`AccumulatorMode.PER_THREAD` (with
    ``nthreads``) or :attr:`AccumulatorMode.SINGLE_ATOMIC`.  Per-buffer ``base_bytes`` are
    taken verbatim from the legacy ``bytes_for`` so byte sizes stay in parity.
    """

    if fock_accumulator not in (AccumulatorMode.PER_THREAD, AccumulatorMode.SINGLE_ATOMIC):
        raise ValueError("fock_accumulator must be PER_THREAD or SINGLE_ATOMIC")
    fock_threads = nthreads if fock_accumulator is AccumulatorMode.PER_THREAD else 1

    buffers = (
        BufferSpec("ids", plan.bytes_for("ids"), ResidencyClass.MIRRORED, "eri_index"),
        BufferSpec("integrals", plan.bytes_for("integrals"), ResidencyClass.MIRRORED, "eri_value"),
        BufferSpec("density", plan.bytes_for("density"), ResidencyClass.MIRRORED, "input_matrix"),
        BufferSpec(
            "fock",
            plan.bytes_for("fock"),
            ResidencyClass.DEVICE_RESIDENT,
            "output_matrix",
            accumulator=fock_accumulator,
            nthreads=fock_threads,
        ),
    )
    return WorkspaceSchema(target=target, shape_key=tuple(plan.reuse_key), buffers=buffers)


def xc_response_schema_from_plan(
    plan: Any,
    *,
    dtype_bytes: int = 8,
    target: str = GpuTarget.XC_RESPONSE.value,
) -> WorkspaceSchema:
    """Map an ``XcResponseCachePlan`` onto the unified schema.

    Buffer order and scalar counts match the legacy ``workspace_layout``; ``base_bytes``
    are ``scalar_values * dtype_bytes`` so the byte arena is in parity with the legacy
    ``total_workspace_bytes``.  Residency reflects ownership: AO-on-grid and quadrature
    weights are produced and owned by the grid machinery (borrowed); grid density is
    mirrored; the XC potential is pure device scratch.
    """

    if dtype_bytes <= 0:
        raise ValueError(f"dtype_bytes must be positive, got {dtype_bytes}")

    buffers = (
        BufferSpec("density", plan.density_values * dtype_bytes, ResidencyClass.MIRRORED, "grid_density"),
        BufferSpec("potential", plan.potential_values * dtype_bytes, ResidencyClass.DEVICE_RESIDENT, "xc_potential"),
        BufferSpec("weights", plan.weight_values * dtype_bytes, ResidencyClass.BORROWED, "quadrature_weight"),
        BufferSpec("ao_grid", plan.ao_grid_values * dtype_bytes, ResidencyClass.BORROWED, "ao_on_grid"),
    )
    return WorkspaceSchema(target=target, shape_key=tuple(plan.reuse_key()), buffers=buffers)
