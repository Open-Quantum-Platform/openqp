"""Stage-1 unified GPU workspace manager (residency + reuse cache).

This module is the single source of truth for *where* a planned GPU buffer
lives (its residency class) and *when* a previously planned allocation can be
reused.  It deliberately subsumes the two legacy planning experiments:

- ``gpu_metc_buffers`` / ``gpu_metc_persistent_runtime`` (persistent METC
  device-array planning + reuse registry), and
- ``tdhf_xc_response_cache`` (TDHF/TDDFT XC-response cache layout planning).

Stage-1 scope is pure planning and bookkeeping.  As with the legacy modules,
nothing here allocates CUDA memory or imports the built OpenQP runtime; it is a
dependency-light, source-testable ABI contract that future Fortran/CUDA wiring
can validate against.

Two design rules keep the unified schema honest:

1. **Buffers are namespaced by target.**  Every record is keyed by
   ``(target, reuse_key)`` so a METC ``density`` buffer and an XC-response
   ``density`` buffer never collide even though they share a role name.

2. **The f3 accumulator policy is explicit.**  The METC Fock/``f3`` accumulator
   defaults to :attr:`F3AccumulatorPolicy.PER_THREAD` (mirroring the
   per-OpenMP-thread ``this%f3`` host model).  ``SINGLE_ATOMIC`` is supported
   but must be requested deliberately -- the manager never silently collapses
   the accumulator to one global buffer.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from typing import Any, Iterable, Optional


# Concrete target namespaces understood by the Stage-1 manager.  The workspace
# manager namespaces buffers per concrete target string (handoff decision 4);
# the GpuConfig may carry a multi-valued ``target`` but residency/cache records
# are always stored under one concrete target at a time.
VALID_TARGETS = ("metc", "xc_response")


class Residency(Enum):
    """Where a planned workspace buffer physically lives.

    ``HOST_ONLY``       -- staged on the host; never uploaded.
    ``DEVICE_RESIDENT`` -- lives on the device for the lifetime of the plan.
    ``MIRRORED``        -- exists on both host and device (upload/download).
    ``BORROWED``        -- aliases memory owned by another allocation; the
                           borrower must not free it.
    """

    HOST_ONLY = "host_only"
    DEVICE_RESIDENT = "device_resident"
    MIRRORED = "mirrored"
    BORROWED = "borrowed"


class F3AccumulatorPolicy(Enum):
    """How the METC Fock/``f3`` accumulator is laid out across threads.

    ``PER_THREAD``    -- one accumulator per OpenMP thread (default; mirrors the
                         host ``this%f3`` residency model).
    ``SINGLE_ATOMIC`` -- one shared accumulator updated with atomics.  Supported
                         but never assumed implicitly.
    """

    PER_THREAD = "per_thread"
    SINGLE_ATOMIC = "single_atomic"


# --- ABI codes shared with the source-side workspace manifest bridge ---------
#
# These integer codes are the contract between this Python control plane and
# ``source/gpu_workspace_bridge.F90``.  The bridge's parity test asserts that
# the Fortran ``GPU_WS_*`` parameters carry exactly these values, so the two
# sides never silently disagree on residency, target namespace, or f3 policy.
WORKSPACE_MANIFEST_SCHEMA = 1
WORKSPACE_ELEMENT_BYTES = 8  # float64-only ABI

TARGET_ABI_CODES = {
    "metc": 0,
    "xc_response": 1,
}

RESIDENCY_ABI_CODES = {
    Residency.HOST_ONLY: 0,
    Residency.DEVICE_RESIDENT: 1,
    Residency.MIRRORED: 2,
    Residency.BORROWED: 3,
}

F3_POLICY_ABI_CODES = {
    F3AccumulatorPolicy.PER_THREAD: 0,
    F3AccumulatorPolicy.SINGLE_ATOMIC: 1,
}


@dataclass(frozen=True)
class WorkspaceBuffer:
    """One planned buffer in a target's workspace table.

    ``bytes`` is the device byte count for the buffer.  ``role`` carries the
    legacy ABI role string (e.g. ``"eri_index"``, ``"input_matrix"``) so the
    unified table stays comparable with the legacy allocation manifests.
    """

    name: str
    bytes: int
    role: str
    residency: Residency


@dataclass(frozen=True)
class WorkspaceKey:
    """Namespaced cache key: a concrete target plus the subsystem reuse key."""

    target: str
    reuse_key: tuple


@dataclass(frozen=True)
class WorkspaceAllocation:
    """A validated, reusable workspace record for one ``(target, reuse_key)``.

    ``scalar_layout`` preserves the legacy XC-response ``(name, offset, length)``
    scalar slices so callers can apply their own dtype width while keeping the
    cache regions non-overlapping.  It is empty for byte-oriented plans (METC).
    """

    key: WorkspaceKey
    buffers: tuple[WorkspaceBuffer, ...]
    total_bytes: int
    f3_policy: Optional[F3AccumulatorPolicy] = None
    nthreads: int = 1
    scalar_layout: tuple = ()
    borrowed: bool = False


def _reject_bool_int(name: str, value: Any) -> int:
    """Return ``int(value)`` rejecting bool aliasing and non-positive values."""

    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive integer, not bool")
    if not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


class GpuWorkspaceManager:
    """Unified residency + reuse registry, namespaced per concrete target.

    The four C-ABI-ready verbs are :meth:`allocate`, :meth:`validate_table`,
    :meth:`borrow`, and :meth:`release`.  Plans are duck-typed so source-level
    tests can drive the manager without importing the full OpenQP package.
    """

    def __init__(self) -> None:
        self._records: dict[WorkspaceKey, WorkspaceAllocation] = {}
        self._borrows: dict[WorkspaceKey, int] = {}

    # -- helpers ---------------------------------------------------------

    @staticmethod
    def _check_target(target: str) -> str:
        if target not in VALID_TARGETS:
            raise ValueError(
                f"unknown workspace target {target!r}; expected one of {VALID_TARGETS}"
            )
        return target

    def _make_key(self, target: str, reuse_key: Any) -> WorkspaceKey:
        self._check_target(target)
        if not isinstance(reuse_key, tuple) or not reuse_key:
            raise ValueError("reuse_key must be a non-empty tuple")
        return WorkspaceKey(target=target, reuse_key=tuple(reuse_key))

    # -- verb: validate_table -------------------------------------------

    def validate_table(
        self, target: str, table: Iterable[Any]
    ) -> tuple[WorkspaceBuffer, ...]:
        """Normalize and validate a workspace table for ``target``.

        Accepts rows that are either :class:`WorkspaceBuffer` instances or
        ``(name, bytes, role, residency)`` tuples.  Returns the validated,
        ordered tuple of buffers, rejecting drift in types, byte counts, or
        residency classes before any allocation path exists.
        """

        self._check_target(target)
        validated: list[WorkspaceBuffer] = []
        seen: set[str] = set()
        for index, row in enumerate(table):
            if isinstance(row, WorkspaceBuffer):
                name, nbytes, role, residency = (
                    row.name,
                    row.bytes,
                    row.role,
                    row.residency,
                )
            else:
                try:
                    name, nbytes, role, residency = row
                except (TypeError, ValueError):
                    raise ValueError(
                        f"workspace row {index} must be (name, bytes, role, residency)"
                    )
            if not isinstance(name, str) or not name:
                raise ValueError(f"workspace row {index} name must be a non-empty str")
            if name in seen:
                raise ValueError(f"duplicate workspace buffer name {name!r}")
            seen.add(name)
            nbytes = _reject_bool_int(f"{name}.bytes", nbytes)
            if not isinstance(role, str) or not role:
                raise ValueError(f"{name}.role must be a non-empty str")
            if not isinstance(residency, Residency):
                raise ValueError(
                    f"{name}.residency must be a Residency, got {residency!r}"
                )
            validated.append(
                WorkspaceBuffer(name=name, bytes=nbytes, role=role, residency=residency)
            )
        if not validated:
            raise ValueError("workspace table must contain at least one buffer")
        return tuple(validated)

    # -- verb: allocate --------------------------------------------------

    def allocate(
        self,
        target: str,
        reuse_key: Any,
        table: Iterable[Any],
        *,
        f3_policy: Optional[F3AccumulatorPolicy] = None,
        nthreads: int = 1,
        scalar_layout: tuple = (),
    ) -> WorkspaceAllocation:
        """Validate ``table`` and record a reusable allocation under ``target``.

        Re-allocating the same ``(target, reuse_key)`` returns the cached record
        if the validated table matches, mirroring the legacy reuse behavior.
        """

        key = self._make_key(target, reuse_key)
        buffers = self.validate_table(target, table)
        nthreads = _reject_bool_int("nthreads", nthreads)
        if f3_policy is not None and not isinstance(f3_policy, F3AccumulatorPolicy):
            raise ValueError("f3_policy must be an F3AccumulatorPolicy or None")

        record = WorkspaceAllocation(
            key=key,
            buffers=buffers,
            total_bytes=sum(buf.bytes for buf in buffers),
            f3_policy=f3_policy,
            nthreads=nthreads,
            scalar_layout=tuple(scalar_layout),
        )

        existing = self._records.get(key)
        if existing is not None and existing.buffers != record.buffers:
            raise ValueError(
                f"workspace reuse conflict for {key.target}:{key.reuse_key}: "
                "table differs from the cached allocation"
            )
        self._records[key] = record
        return record

    # -- verb: borrow ----------------------------------------------------

    def borrow(self, target: str, reuse_key: Any) -> WorkspaceAllocation:
        """Return a BORROWED alias of an existing owned allocation.

        The underlying record keeps ownership of the memory; the returned alias
        marks every buffer (and the record) as :attr:`Residency.BORROWED` so a
        consumer cannot accidentally free buffers it does not own.  Each borrow
        must be matched by a :meth:`release`.
        """

        key = self._make_key(target, reuse_key)
        owner = self._records.get(key)
        if owner is None:
            raise KeyError(
                f"cannot borrow {key.target}:{key.reuse_key}: not allocated"
            )
        self._borrows[key] = self._borrows.get(key, 0) + 1
        borrowed_buffers = tuple(
            replace(buf, residency=Residency.BORROWED) for buf in owner.buffers
        )
        return replace(owner, buffers=borrowed_buffers, borrowed=True)

    # -- verb: release ---------------------------------------------------

    def release(self, target: str, reuse_key: Any) -> Optional[WorkspaceAllocation]:
        """Release one outstanding borrow, or free the owning allocation.

        If borrows are outstanding for ``(target, reuse_key)``, one borrow is
        returned and the count decremented.  Once no borrows remain, releasing
        forgets and returns the owning record (or ``None`` if absent).
        """

        key = self._make_key(target, reuse_key)
        outstanding = self._borrows.get(key, 0)
        if outstanding > 0:
            self._borrows[key] = outstanding - 1
            if self._borrows[key] == 0:
                del self._borrows[key]
            owner = self._records.get(key)
            if owner is None:
                return None
            borrowed_buffers = tuple(
                replace(buf, residency=Residency.BORROWED) for buf in owner.buffers
            )
            return replace(owner, buffers=borrowed_buffers, borrowed=True)
        return self._records.pop(key, None)

    # -- lookup ----------------------------------------------------------

    def lookup(self, target: str, reuse_key: Any) -> Optional[WorkspaceAllocation]:
        """Return the owned record for ``(target, reuse_key)`` if registered."""

        return self._records.get(self._make_key(target, reuse_key))

    # -- parity adapters -------------------------------------------------

    def from_metc_plan(
        self,
        plan: Any,
        *,
        f3_policy: F3AccumulatorPolicy = F3AccumulatorPolicy.PER_THREAD,
        nthreads: int = 1,
    ) -> WorkspaceAllocation:
        """Map a legacy ``PersistentMetcBufferPlan`` into the unified schema.

        The legacy allocation manifest (``ids``, ``integrals``, ``density``,
        ``fock``) becomes a namespaced ``metc`` workspace table.  The ``fock``
        buffer is the ``f3`` accumulator: under ``PER_THREAD`` it is replicated
        across ``nthreads`` (mirroring per-OpenMP-thread ``this%f3``); under
        ``SINGLE_ATOMIC`` a single shared accumulator is kept.
        """

        if not isinstance(f3_policy, F3AccumulatorPolicy):
            raise ValueError("f3_policy must be an F3AccumulatorPolicy")
        nthreads = _reject_bool_int("nthreads", nthreads)

        residency_by_role = {
            "eri_index": Residency.MIRRORED,
            "eri_value": Residency.MIRRORED,
            "input_matrix": Residency.MIRRORED,
            "output_matrix": Residency.MIRRORED,
        }
        table: list[WorkspaceBuffer] = []
        for item in plan.allocation_manifest():
            name = str(item["name"])
            nbytes = int(item["bytes"])
            role = str(item["role"])
            # The Fock/f3 accumulator is replicated per thread under PER_THREAD.
            if role == "output_matrix" and f3_policy is F3AccumulatorPolicy.PER_THREAD:
                nbytes *= nthreads
            table.append(
                WorkspaceBuffer(
                    name=name,
                    bytes=nbytes,
                    role=role,
                    residency=residency_by_role.get(role, Residency.MIRRORED),
                )
            )
        return self.allocate(
            "metc",
            plan.reuse_key,
            table,
            f3_policy=f3_policy,
            nthreads=nthreads,
        )

    def from_metc_session(
        self,
        plan: Any,
        *,
        nthreads: int,
        max_ncur: int,
        f3_policy: F3AccumulatorPolicy = F3AccumulatorPolicy.PER_THREAD,
    ) -> WorkspaceAllocation:
        """Build the METC *residency-session* workspace (METC-C1b).

        Unlike :meth:`from_metc_plan` (the pre-refit MIRRORED planning view),
        this models the resident arena the source-side session actually owns:

            [ density(d3) | fock(f3) | ids_scratch | ints_scratch ]

        - ``density`` (d3): shared, read-only input tensor; DEVICE_RESIDENT.
        - ``fock`` (f3): PER_THREAD accumulator (``nthreads x base``);
          DEVICE_RESIDENT.  ``SINGLE_ATOMIC`` keeps a single base accumulator
          (schema-only; not the default refit path).
        - ``ids_scratch`` / ``ints_scratch``: PER_THREAD scratch sized for
          ``max_ncur`` integrals, so concurrent OpenMP ``update`` calls never
          share-write a resident buffer.

        Buffer order and byte sizes mirror ``oqp_gpu_metc_layout`` in
        ``source/gpu_workspace_runtime.c`` so the Python manifest and the
        source-side arena agree on offsets.
        """

        if not isinstance(f3_policy, F3AccumulatorPolicy):
            raise ValueError("f3_policy must be an F3AccumulatorPolicy")
        nthreads = _reject_bool_int("nthreads", nthreads)
        max_ncur = _reject_bool_int("max_ncur", max_ncur)

        bytes_tensor = int(plan.bytes_for("density"))  # nf*nmatrix*nbf*nbf*8
        ids_cap = 4 * max_ncur * 4   # four int32 indices per integral
        ints_cap = max_ncur * 8      # one float64 per integral
        fock_factor = nthreads if f3_policy is F3AccumulatorPolicy.PER_THREAD else 1

        table = [
            WorkspaceBuffer("density", bytes_tensor, "input_matrix",
                            Residency.DEVICE_RESIDENT),
            WorkspaceBuffer("fock", bytes_tensor * fock_factor, "output_matrix",
                            Residency.DEVICE_RESIDENT),
            WorkspaceBuffer("ids_scratch", ids_cap * nthreads, "index_scratch",
                            Residency.DEVICE_RESIDENT),
            WorkspaceBuffer("ints_scratch", ints_cap * nthreads, "integral_scratch",
                            Residency.DEVICE_RESIDENT),
        ]
        return self.allocate(
            "metc",
            tuple(plan.reuse_key) + (nthreads, max_ncur),
            table,
            f3_policy=f3_policy,
            nthreads=nthreads,
        )

    def from_xc_response_plan(
        self, plan: Any, *, dtype_bytes: int = 8
    ) -> WorkspaceAllocation:
        """Map a legacy ``XcResponseCachePlan`` into the unified schema.

        The legacy scalar workspace layout (``density``, ``potential``,
        ``weights``, ``ao_grid``) is preserved verbatim in ``scalar_layout`` so
        the cache offsets survive the mapping, while each slice also becomes a
        byte-sized :class:`WorkspaceBuffer` namespaced under ``xc_response``.
        """

        dtype_bytes = _reject_bool_int("dtype_bytes", dtype_bytes)
        layout = plan.workspace_layout()
        table = [
            WorkspaceBuffer(
                name=name,
                bytes=length * dtype_bytes,
                role=name,
                residency=Residency.MIRRORED,
            )
            for name, _offset, length in layout
        ]
        return self.allocate(
            "xc_response",
            plan.reuse_key(),
            table,
            scalar_layout=tuple(layout),
        )


# --- ABI manifest export (source-side bridge contract) -----------------------


def workspace_manifest_rows(allocation: WorkspaceAllocation) -> tuple[dict, ...]:
    """Return ABI manifest rows with contiguous byte offsets.

    Each row mirrors one ``gpu_ws_buffer_t`` on the Fortran/CUDA side: slot,
    target (name + code), logical name, role, residency (name + code), byte
    size, contiguous byte offset, and element size.  Offsets are assigned in
    buffer order so the source-side bridge can validate the identical layout.
    """

    rows: list[dict] = []
    offset = 0
    for slot, buf in enumerate(allocation.buffers, start=1):
        rows.append(
            {
                "slot": slot,
                "target": allocation.key.target,
                "target_code": TARGET_ABI_CODES[allocation.key.target],
                "name": buf.name,
                "role": buf.role,
                "residency": buf.residency.name,
                "residency_code": RESIDENCY_ABI_CODES[buf.residency],
                "bytes": buf.bytes,
                "offset": offset,
                "element_bytes": WORKSPACE_ELEMENT_BYTES,
            }
        )
        offset += buf.bytes
    return tuple(rows)


def workspace_manifest(allocation: WorkspaceAllocation) -> dict:
    """Return the full ABI manifest for a workspace allocation.

    The manifest is the source-side bridge contract: it carries the schema
    version, the namespaced target, the reuse key, the (optional) f3
    accumulator policy with its ABI code, the per-thread fan-out, the ordered
    buffer rows, and the total byte size.  ``f3_policy`` is ``None`` for targets
    that do not carry an accumulator policy (e.g. ``xc_response``); it is never
    silently defaulted to a single global accumulator.
    """

    f3_policy = allocation.f3_policy
    return {
        "schema": WORKSPACE_MANIFEST_SCHEMA,
        "target": allocation.key.target,
        "target_code": TARGET_ABI_CODES[allocation.key.target],
        "reuse_key": allocation.key.reuse_key,
        "f3_policy": f3_policy.name if f3_policy is not None else None,
        "f3_policy_code": (
            F3_POLICY_ABI_CODES[f3_policy] if f3_policy is not None else None
        ),
        "nthreads": allocation.nthreads,
        "element_bytes": WORKSPACE_ELEMENT_BYTES,
        "rows": workspace_manifest_rows(allocation),
        "total_bytes": allocation.total_bytes,
    }
