"""Persistent GPU buffer planning helpers for METC experiments.

The persistent-buffer branch is about removing repeated device allocation and
host/device setup overhead from GPU METC contractions.  Keep this helper pure
Python so sizing and cache-key behavior can be tested without CUDA or a built
OpenQP shared library.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PersistentMetcBufferPlan:
    """Byte-size plan for reusable METC device arrays."""

    nbf: int
    nf: int
    nmatrix: int
    max_integrals: int
    dtype_bytes: int = 8

    @classmethod
    def from_problem(
        cls,
        *,
        nbf: int,
        nf: int,
        nmatrix: int,
        max_integrals: int,
        dtype_bytes: int = 8,
    ) -> "PersistentMetcBufferPlan":
        """Build a validated persistent-buffer plan for a METC problem shape."""

        values = {
            "nbf": nbf,
            "nf": nf,
            "nmatrix": nmatrix,
            "max_integrals": max_integrals,
            "dtype_bytes": dtype_bytes,
        }
        for name, value in values.items():
            if isinstance(value, bool):
                raise ValueError(f"{name} must be a positive integer, not bool")
            if int(value) <= 0:
                raise ValueError(f"{name} must be positive")
        return cls(
            nbf=int(nbf),
            nf=int(nf),
            nmatrix=int(nmatrix),
            max_integrals=int(max_integrals),
            dtype_bytes=int(dtype_bytes),
        )

    @property
    def reuse_key(self) -> tuple[int, int, int, int, int]:
        """Stable cache key for deciding when device buffers can be reused."""

        return (self.nbf, self.nf, self.nmatrix, self.max_integrals, self.dtype_bytes)

    def bytes_for(self, name: str) -> int:
        """Return the allocation size for one named persistent device array."""

        if name == "ids":
            # Four shell/basis indices per integral, stored as 32-bit integers.
            return self.max_integrals * 4 * 4
        if name == "integrals":
            return self.max_integrals * self.dtype_bytes
        if name in {"density", "fock"}:
            return self.nmatrix * self.nf * self.nbf * self.nbf * self.dtype_bytes
        raise KeyError(name)

    def allocation_manifest(self) -> tuple[dict[str, int | str], ...]:
        """Ordered ABI-facing allocation manifest for planned device buffers.

        The manifest is intentionally simple (name, byte count, semantic role) so
        future Fortran/CUDA wiring can validate allocation order and reuse without
        importing a built OpenQP runtime in source-level tests.
        """

        roles = {
            "ids": "eri_index",
            "integrals": "eri_value",
            "density": "input_matrix",
            "fock": "output_matrix",
        }
        return tuple(
            {"name": name, "bytes": self.bytes_for(name), "role": roles[name]}
            for name in ("ids", "integrals", "density", "fock")
        )

    def fortran_allocation_table(self) -> tuple[tuple[int, str, int, str], ...]:
        """Return a stable one-based allocation table for Fortran/C ABI stubs."""

        return tuple(
            (slot, str(item["name"]), int(item["bytes"]), str(item["role"]))
            for slot, item in enumerate(self.allocation_manifest(), start=1)
        )

    def validate_fortran_allocation_table(
        self, table: tuple[tuple[int, str, int, str], ...]
    ) -> tuple[tuple[int, str, int, str], ...]:
        """Validate an ABI allocation table against this plan.

        Future Fortran/CUDA stubs can use this pure-Python contract to catch
        drift in slot order, names, byte counts, or semantic roles before any
        device allocation path exists.
        """

        expected = self.fortran_allocation_table()
        if tuple(table) != expected:
            for expected_row, actual_row in zip(expected, tuple(table)):
                if actual_row != expected_row:
                    raise ValueError(
                        f"allocation table mismatch at slot {expected_row[0]}: "
                        f"expected {expected_row!r}, got {actual_row!r}"
                    )
            raise ValueError(
                f"allocation table length mismatch: expected {len(expected)}, got {len(tuple(table))}"
            )
        return expected

    @property
    def total_bytes(self) -> int:
        """Total bytes needed by all currently planned persistent buffers."""

        return sum(
            self.bytes_for(name) for name in ("ids", "integrals", "density", "fock")
        )
