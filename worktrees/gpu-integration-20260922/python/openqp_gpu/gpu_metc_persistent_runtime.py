"""Pure-Python persistent METC allocation runtime scaffolding.

This module deliberately does not allocate CUDA memory.  It is a source-level
contract for future Fortran/CUDA wiring: a caller must provide the allocation
table generated from :mod:`gpu_metc_buffers`, and the registry validates that
ABI table before remembering the reusable plan metadata.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


AllocationTable = tuple[tuple[int, str, int, str], ...]
ReuseKey = tuple[int, int, int, int, int]


def _validate_reuse_key(reuse_key: Any) -> ReuseKey:
    """Return a strict METC reuse key, rejecting bool-as-int aliasing."""

    if not isinstance(reuse_key, tuple) or len(reuse_key) != 5:
        raise ValueError("reuse_key must be a five-integer METC allocation key")
    if any(not isinstance(value, int) or isinstance(value, bool) for value in reuse_key):
        raise ValueError("reuse_key must contain non-boolean integers")
    return reuse_key


@dataclass(frozen=True)
class PersistentMetcAllocationRecord:
    """Validated metadata for one reusable METC allocation shape."""

    reuse_key: ReuseKey
    total_bytes: int
    table: AllocationTable


class PersistentMetcAllocationRegistry:
    """Record validated persistent-buffer plans without touching CUDA."""

    def __init__(self) -> None:
        self._records: dict[ReuseKey, PersistentMetcAllocationRecord] = {}

    def register_plan(self, plan: Any, table: AllocationTable) -> PersistentMetcAllocationRecord:
        """Validate and remember a persistent METC allocation plan.

        ``plan`` is intentionally duck-typed so source-level tests can load this
        module without importing the full OpenQP package.  The real runtime path
        is expected to pass a ``PersistentMetcBufferPlan`` instance.
        """

        validated_table = plan.validate_fortran_allocation_table(table)
        reuse_key = _validate_reuse_key(plan.reuse_key)
        record = PersistentMetcAllocationRecord(
            reuse_key=reuse_key,
            total_bytes=plan.total_bytes,
            table=validated_table,
        )
        self._records[record.reuse_key] = record
        return record

    def lookup(self, reuse_key: ReuseKey) -> PersistentMetcAllocationRecord | None:
        """Return the validated record for ``reuse_key`` if one is registered."""

        return self._records.get(_validate_reuse_key(reuse_key))

    def release(self, reuse_key: ReuseKey) -> PersistentMetcAllocationRecord | None:
        """Forget and return the record for ``reuse_key`` if one is registered."""

        return self._records.pop(_validate_reuse_key(reuse_key), None)
