"""Planning helpers for TDHF/TDDFT XC-response cache workspaces.

This module is intentionally dependency-light so source-level tests can define
cache contracts before the Fortran/CUDA response kernels are wired in.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class XcResponseCachePlan:
    """Describe reusable XC-response intermediates for one TDHF job shape.

    The plan records dimensions and identity fields that must remain stable for
    cached quadrature intermediates to be reused safely.  It does not allocate
    arrays or imply that production response code is already cached.
    """

    nbf: int
    ngrid: int
    functional: str
    basis: str
    scf_type: str
    response_type: str
    spin_channels: int = 1

    def __post_init__(self):
        for name in ("nbf", "ngrid", "spin_channels"):
            value = getattr(self, name)
            if value <= 0:
                raise ValueError(f"{name} must be positive, got {value}")

    def reuse_key(self):
        """Return the stable identity key for a compatible cache entry."""

        return (
            self.functional.lower(),
            self.basis.lower(),
            self.scf_type.lower(),
            self.response_type.lower(),
            self.nbf,
            self.ngrid,
            self.spin_channels,
        )

    @property
    def density_values(self) -> int:
        """Scalar density slots over the integration grid."""

        return self.ngrid * self.spin_channels

    @property
    def potential_values(self) -> int:
        """Scalar XC-potential slots over the integration grid."""

        return self.ngrid * self.spin_channels

    @property
    def weight_values(self) -> int:
        """Quadrature weights needed by the cached response contraction."""

        return self.ngrid

    @property
    def ao_grid_values(self) -> int:
        """AO-on-grid slots reused by the TDHF/TDDFT XC response."""

        return self.nbf * self.ngrid

    def total_scalar_values(self) -> int:
        """Total scalar workspace values covered by this cache plan."""

        return (
            self.density_values
            + self.potential_values
            + self.weight_values
            + self.ao_grid_values
        )

    def total_workspace_bytes(self, dtype_bytes: int = 8) -> int:
        """Total workspace size in bytes for one scalar dtype width."""

        if dtype_bytes <= 0:
            raise ValueError(f"dtype_bytes must be positive, got {dtype_bytes}")
        return self.total_scalar_values() * dtype_bytes

    def invalidation_reasons(self, other):
        """Return cache-identity fields that differ from another plan.

        String fields are compared case-insensitively to match ``reuse_key``.
        The returned field names are ordered for stable diagnostics.
        """

        fields = (
            "functional",
            "basis",
            "scf_type",
            "response_type",
            "nbf",
            "ngrid",
            "spin_channels",
        )
        changed = []
        for field in fields:
            lhs = getattr(self, field)
            rhs = getattr(other, field)
            if isinstance(lhs, str) and isinstance(rhs, str):
                lhs = lhs.lower()
                rhs = rhs.lower()
            if lhs != rhs:
                changed.append(field)
        return tuple(changed)

    def workspace_layout(self):
        """Return ordered scalar-workspace slices as ``(name, offset, length)``.

        The layout is a stable call-site contract for future Fortran/CUDA wiring:
        offsets are scalar indices, not byte offsets, so callers can apply their
        own dtype width while preserving non-overlapping cache regions.
        """

        layout = []
        offset = 0
        for name, length in (
            ("density", self.density_values),
            ("potential", self.potential_values),
            ("weights", self.weight_values),
            ("ao_grid", self.ao_grid_values),
        ):
            layout.append((name, offset, length))
            offset += length
        return tuple(layout)
