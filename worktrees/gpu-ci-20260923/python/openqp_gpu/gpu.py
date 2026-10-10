"""Runtime helpers for experimental OpenQP GPU backends.

The GPU branches are intentionally split by target.  This helper stays pure
Python so input parsing and branch-level unit tests do not require CUDA
libraries or a built OpenQP shared library.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any


XC_RESPONSE_STATUS_OK = 0
XC_RESPONSE_STATUS_DISABLED = 1
XC_RESPONSE_STATUS_INVALID_INPUT = 2
XC_RESPONSE_STATUS_OVERFLOW = 3
GPU_FALLBACK_POLICIES = {"cpu", "error"}

# Concrete GPU acceleration targets understood by the unified config surface.
GPU_TARGETS = ("metc", "xc_response")

# Environment variables that gate the experimental GPU METC residency path.
# These names are shared verbatim with the Fortran backend
# (source/gpu_backend.F90) and the C runtime (source/gpu_workspace_runtime.c)
# so the Python tooling, the runtime, and the validation/benchmark harness all
# agree without a dedicated binding.
OQP_GPU_METC_ENV = "OQP_GPU_METC"
OQP_GPU_METC_STRICT_ENV = "OQP_GPU_METC_STRICT"
OQP_GPU_METC_TIMING_ENV = "OQP_GPU_METC_TIMING"
OQP_GPU_STRICT_ENV = "OQP_GPU_STRICT"
OQP_GPU_F3_CHECK_ENV = "OQP_GPU_F3_CHECK"
OQP_GPU_PROFILE_ENV = "OQP_GPU_PROFILE"
OQP_GPU_METC_PROFILE_ENV = "OQP_GPU_METC_PROFILE"

_GPU_METC_TRUTHY = {"1", "true", "on", "yes"}


def _gpu_metc_env_flag(name: str, environ: Any = None) -> bool:
    """Return True when ``name`` is set to a recognized truthy token."""

    source = os.environ if environ is None else environ
    value = source.get(name)
    if value is None:
        return False
    return str(value).strip().lower() in _GPU_METC_TRUTHY


def gpu_metc_enabled_from_env(environ: Any = None) -> bool:
    """Return True when the GPU METC path is requested via the environment."""

    return _gpu_metc_env_flag(OQP_GPU_METC_ENV, environ)


def gpu_metc_strict(environ: Any = None) -> bool:
    """Return True when strict (no-CPU-fallback) METC mode is requested.

    Strict mode makes a GPU METC *begin-time* failure fatal instead of allowing
    a CPU fallback; it mirrors the Fortran ``OQP_GPU_METC_STRICT`` handling.
    (A failure once the resident session is ACTIVE is always fatal, regardless
    of this flag, because partial device results cannot be safely combined with
    a host fallback.)
    """

    return _gpu_metc_env_flag(OQP_GPU_METC_STRICT_ENV, environ) or _gpu_metc_env_flag(
        OQP_GPU_STRICT_ENV, environ
    )


def gpu_metc_timing_enabled(environ: Any = None) -> bool:
    """Return True when the runtime should emit separated-region METC timings."""

    return _gpu_metc_env_flag(OQP_GPU_METC_TIMING_ENV, environ)


def gpu_metc_f3_check(environ: Any = None) -> bool:
    """Return True when validation-only CPU mirror f3 checking is requested."""

    return _gpu_metc_env_flag(OQP_GPU_F3_CHECK_ENV, environ) or _gpu_metc_env_flag(
        "OQP_GPU_METC_F3_CHECK", environ
    )


def gpu_metc_profile_enabled(environ: Any = None) -> bool:
    """Return True when detailed GPU METC profile counters should be emitted."""

    return _gpu_metc_env_flag(OQP_GPU_PROFILE_ENV, environ) or _gpu_metc_env_flag(
        OQP_GPU_METC_PROFILE_ENV, environ
    )


def gpu_metc_runtime_config(environ: Any = None) -> dict:
    """Return a snapshot of the env-driven GPU METC runtime configuration."""

    return {
        "enabled": gpu_metc_enabled_from_env(environ),
        "strict": gpu_metc_strict(environ),
        "timing": gpu_metc_timing_enabled(environ),
        "f3_check": gpu_metc_f3_check(environ),
        "profile": gpu_metc_profile_enabled(environ),
    }


def normalize_gpu_targets(value: Any) -> tuple[str, ...]:
    """Normalize a GPU target specification into an immutable, ordered tuple.

    Implements the locked multi-valued target policy with backward-compatible
    singleton coercion::

        "metc"                     -> ("metc",)
        "xc_response"              -> ("xc_response",)
        ("metc", "xc_response")    -> ("metc", "xc_response")
        ["metc", "xc_response"]    -> ("metc", "xc_response")

    First-occurrence order is preserved and duplicates are collapsed so the
    result is deterministic.  Booleans, empty specifications, non-string
    entries, and unknown targets are rejected.
    """

    if isinstance(value, bool):
        raise ValueError("gpu.target must not be a boolean")
    if isinstance(value, str):
        items = [value]
    elif isinstance(value, (list, tuple)):
        items = list(value)
    else:
        raise ValueError(
            f"gpu.target must be a string or sequence of strings, got {type(value).__name__}"
        )
    if not items:
        raise ValueError("gpu.target must not be empty")
    normalized: list[str] = []
    for item in items:
        if isinstance(item, bool):
            raise ValueError("gpu.target entries must not be booleans")
        if not isinstance(item, str):
            raise ValueError(
                f"gpu.target entries must be strings, got {type(item).__name__}"
            )
        name = item.strip().lower()
        if name not in GPU_TARGETS:
            raise ValueError(
                f"unknown gpu.target {item!r}; expected a subset of {GPU_TARGETS}"
            )
        if name not in normalized:
            normalized.append(name)
    return tuple(normalized)


def _normalize_xc_response_status(status: int) -> int:
    """Return an integer ABI status while rejecting bool-as-int ambiguity."""

    if isinstance(status, bool):
        raise ValueError("XC-response runtime status must be a non-boolean integer")
    return int(status)


def xc_response_status_message(status: int) -> str:
    """Return Python-side text matching the Fortran XC-response ABI statuses."""

    status = _normalize_xc_response_status(status)
    messages = {
        XC_RESPONSE_STATUS_OK: "CUDA XC-response contract completed",
        XC_RESPONSE_STATUS_DISABLED: "CUDA XC-response backend disabled; using CPU fallback",
        XC_RESPONSE_STATUS_INVALID_INPUT: "CUDA XC-response invalid input dimensions or null buffer",
        XC_RESPONSE_STATUS_OVERFLOW: "CUDA XC-response dimension overflow before launch",
    }
    return messages.get(status, f"CUDA XC-response runtime status {status}")


def xc_response_status_report(status: int) -> dict[str, Any]:
    """Classify an XC-response ABI status for future Python/Fortran handoffs."""

    status = _normalize_xc_response_status(status)
    known_statuses = {
        XC_RESPONSE_STATUS_OK,
        XC_RESPONSE_STATUS_DISABLED,
        XC_RESPONSE_STATUS_INVALID_INPUT,
        XC_RESPONSE_STATUS_OVERFLOW,
    }
    return {
        "status": status,
        "message": xc_response_status_message(status),
        "ok": status == XC_RESPONSE_STATUS_OK,
        "fallback": status == XC_RESPONSE_STATUS_DISABLED,
        "blocking": status not in {XC_RESPONSE_STATUS_OK, XC_RESPONSE_STATUS_DISABLED},
        "known": status in known_statuses,
    }


def xc_response_handoff_validation_error(handoff: Any) -> str | None:
    """Return a diagnostic when a CUDA XC-response handoff payload is malformed."""

    if not isinstance(handoff, dict):
        return "Malformed XC-response CUDA handoff: payload must be a dictionary"
    required = {"nbasis", "nstate", "slots", "total_elements", "total_bytes", "launch", "fortran_table", "layout"}
    missing = sorted(required.difference(handoff))
    if missing:
        return "Malformed XC-response CUDA handoff: missing keys " + ", ".join(missing)
    launch = handoff.get("launch")
    if not isinstance(launch, dict):
        return "Malformed XC-response CUDA handoff: launch must be a dictionary"
    launch_missing = sorted({"blocks", "threads_per_block", "elements"}.difference(launch))
    if launch_missing:
        return "Malformed XC-response CUDA handoff: launch missing keys " + ", ".join(launch_missing)
    launch_values = (launch.get("blocks"), launch.get("threads_per_block"), launch.get("elements"))
    if not all(isinstance(value, int) and not isinstance(value, bool) and value > 0 for value in launch_values):
        return (
            "Malformed XC-response CUDA handoff: "
            "launch blocks, threads_per_block, and elements must be non-boolean positive integers"
        )
    slots = handoff.get("slots")
    if launch.get("elements") != slots:
        return "Malformed XC-response CUDA handoff: launch elements must match slots"
    nbasis = handoff.get("nbasis")
    nstate = handoff.get("nstate")
    fortran_table = handoff.get("fortran_table")
    if (
        not isinstance(nbasis, int)
        or not isinstance(nstate, int)
        or isinstance(nbasis, bool)
        or isinstance(nstate, bool)
    ):
        return "Malformed XC-response CUDA handoff: nbasis and nstate must be non-boolean integers"
    if not isinstance(fortran_table, (list, tuple)):
        return "Malformed XC-response CUDA handoff: fortran_table must be a sequence"
    try:
        plan = XcResponseGpuPlan(nbasis=nbasis, nstate=nstate)
        plan.validate_fortran_buffer_table(fortran_table)
    except (TypeError, ValueError) as exc:
        return "Malformed XC-response CUDA handoff: " + str(exc)
    if handoff.get("total_elements") != plan.total_elements:
        return "Malformed XC-response CUDA handoff: total_elements must match validated XC-response plan"
    if handoff.get("total_bytes") != plan.total_bytes:
        return "Malformed XC-response CUDA handoff: total_bytes must match validated XC-response plan"
    layout = handoff.get("layout")
    if not isinstance(layout, dict):
        return "Malformed XC-response CUDA handoff: layout must be a dictionary"
    layout_missing = sorted({"density", "kernel", "response"}.difference(layout))
    if layout_missing:
        return "Malformed XC-response CUDA handoff: layout missing buffers " + ", ".join(layout_missing)
    expected_offset = 0
    for name in ("density", "kernel", "response"):
        entry = layout.get(name)
        if not isinstance(entry, dict):
            return f"Malformed XC-response CUDA handoff: layout entry {name} must be a dictionary"
        elements = entry.get("elements")
        offset = entry.get("offset")
        stop = entry.get("stop")
        if not all(isinstance(value, int) and not isinstance(value, bool) for value in (elements, offset, stop)):
            return f"Malformed XC-response CUDA handoff: layout entry {name} must use non-boolean integer elements/offset/stop"
        if elements != slots:
            return f"Malformed XC-response CUDA handoff: layout elements for {name} must match slots"
        if offset != expected_offset:
            return f"Malformed XC-response CUDA handoff: layout offset for {name} must be contiguous"
        if stop != offset + elements:
            return f"Malformed XC-response CUDA handoff: layout stop for {name} must equal offset plus elements"
        expected_offset = stop
    if expected_offset != handoff.get("total_elements"):
        return "Malformed XC-response CUDA handoff: layout stop must match total_elements"
    return None


def xc_response_call_site_action(decision: dict[str, Any]) -> str:
    """Return the safe action for a future XC-response Python call site."""

    report = decision.get("status_report", {}) or {}
    fallback_policy = str(decision.get("fallback_policy", "cpu")).lower()
    if fallback_policy not in GPU_FALLBACK_POLICIES:
        return "blocked"
    if report.get("ok") and xc_response_handoff_validation_error(decision.get("handoff")) is None:
        return "cuda"
    if report.get("fallback") and fallback_policy != "error":
        return "cpu_fallback"
    return "blocked"


def xc_response_dispatch_plan(decision: dict[str, Any]) -> dict[str, Any]:
    """Translate an XC-response runtime decision into guarded call-site flags."""

    action = xc_response_call_site_action(decision)
    status = _normalize_xc_response_status(decision.get("status", XC_RESPONSE_STATUS_INVALID_INPUT))
    report = decision.get("status_report")
    if not isinstance(report, dict) or _normalize_xc_response_status(report.get("status", status)) != status:
        report = xc_response_status_report(status)
    if action == "cuda":
        reason = str(decision.get("reason", ""))
        handoff = decision.get("handoff")
    elif action == "cpu_fallback":
        reason = str(decision.get("status_message", decision.get("reason", "")))
        handoff = None
    else:
        fallback_policy = str(decision.get("fallback_policy", "cpu")).lower()
        handoff_error = xc_response_handoff_validation_error(decision.get("handoff")) if report.get("ok") else None
        if fallback_policy not in GPU_FALLBACK_POLICIES:
            reason = f"Unsupported GPU fallback policy: {fallback_policy}"
        elif handoff_error is not None:
            reason = handoff_error
        elif report.get("fallback"):
            reason = str(decision.get("status_message", decision.get("reason", "")))
        else:
            reason = str(decision.get("reason", decision.get("status_message", "")))
        handoff = None
    return {
        "action": action,
        "call_cuda": action == "cuda",
        "use_cpu_fallback": action == "cpu_fallback",
        "blocked": action == "blocked",
        "enabled": bool(decision.get("enabled", False)),
        "supported": bool(decision.get("supported", False)),
        "status": status,
        "status_report": report,
        "fallback": bool(report.get("fallback", False)),
        "blocking": bool(report.get("blocking", False) or action == "blocked"),
        "fallback_policy": str(decision.get("fallback_policy", "cpu")),
        "reason": reason,
        "handoff": handoff,
    }


def xc_response_dispatch_audit_record(plan: dict[str, Any]) -> dict[str, Any]:
    """Return a compact JSON-friendly audit record for an XC-response call site."""

    handoff = plan.get("handoff")
    if isinstance(handoff, dict) and xc_response_handoff_validation_error(handoff) is None:
        handoff_valid = True
        launch = handoff.get("launch", {})
        launch_elements = launch.get("elements") if isinstance(launch, dict) else None
        total_bytes = handoff.get("total_bytes")
    else:
        handoff_valid = False
        launch_elements = None
        total_bytes = None
    status_report = plan.get("status_report", {}) or {}
    action = str(plan.get("action", "blocked"))
    status = _normalize_xc_response_status(plan.get("status", XC_RESPONSE_STATUS_INVALID_INPUT))
    diagnostic_code = {
        "cuda": "cuda_ready",
        "cpu_fallback": "cpu_fallback",
    }.get(action, "blocked")
    return {
        "action": action,
        "diagnostic_code": diagnostic_code,
        "status": status,
        "status_known": bool(status_report.get("known", False)),
        "status_blocking": bool(status_report.get("blocking", False)),
        "enabled": bool(plan.get("enabled", False)),
        "supported": bool(plan.get("supported", False)),
        "fallback_policy": str(plan.get("fallback_policy", "cpu")),
        "handoff_valid": handoff_valid,
        "launch_elements": launch_elements,
        "total_bytes": total_bytes,
        "reason": str(plan.get("reason", "")),
    }


def _xc_response_audit_log_value(value: Any) -> str:
    """Return a stable single-token value for XC-response audit log lines."""

    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "none"
    return "_".join(str(value).strip().split())


def xc_response_dispatch_audit_log_line(audit: dict[str, Any]) -> str:
    """Format a compact schema-versioned XC-response audit line.

    The line is intentionally limited to the public audit-record fields and does
    not expose raw CUDA handoff payloads, layout maps, or Fortran buffer tables.
    """

    ordered_fields = [
        "action",
        "diagnostic_code",
        "status",
        "status_known",
        "status_blocking",
        "enabled",
        "supported",
        "fallback_policy",
        "handoff_valid",
        "launch_elements",
        "total_bytes",
        "reason",
    ]
    parts = ["OQP_GPU_XC_RESPONSE_AUDIT", "schema=1"]
    for field in ordered_fields:
        parts.append(f"{field}={_xc_response_audit_log_value(audit.get(field))}")
    return " ".join(parts)


def _parse_xc_response_audit_int_field(fields: dict[str, Any], key: str, *, allow_none: bool = False) -> None:
    """Parse an integer audit field with a field-specific diagnostic."""

    value = fields[key]
    if allow_none and value == "none":
        fields[key] = None
        return
    try:
        fields[key] = int(value)
    except (TypeError, ValueError) as exc:
        suffix = " or none" if allow_none else ""
        raise ValueError(f"XC-response audit field {key} must be an integer{suffix}") from exc


def parse_xc_response_dispatch_audit_log_line(line: str) -> dict[str, Any]:
    """Parse an ``OQP_GPU_XC_RESPONSE_AUDIT`` line for log consumers."""

    prefix = "OQP_GPU_XC_RESPONSE_AUDIT"
    tokens = str(line).strip().split()
    if not tokens or tokens[0] != prefix:
        raise ValueError("line is not an OQP_GPU_XC_RESPONSE_AUDIT record")
    fields: dict[str, Any] = {}
    for token in tokens[1:]:
        if "=" not in token:
            raise ValueError(f"malformed XC-response audit token: {token}")
        key, value = token.split("=", 1)
        if key in fields:
            raise ValueError(f"duplicate XC-response audit field: {key}")
        fields[key] = value
    required_fields = {
        "schema",
        "action",
        "diagnostic_code",
        "status",
        "status_known",
        "status_blocking",
        "enabled",
        "supported",
        "fallback_policy",
        "handoff_valid",
        "launch_elements",
        "total_bytes",
        "reason",
    }
    missing = sorted(required_fields.difference(fields))
    if missing:
        raise ValueError("missing required audit fields: " + ", ".join(missing))
    unknown = sorted(set(fields).difference(required_fields))
    if unknown:
        raise ValueError("unknown XC-response audit fields: " + ", ".join(unknown))
    for key in ("schema", "status"):
        _parse_xc_response_audit_int_field(fields, key)
    for key in ("status_known", "status_blocking", "enabled", "supported", "handoff_valid"):
        if key in fields:
            if fields[key] not in {"true", "false"}:
                raise ValueError(f"XC-response audit field {key} must be true or false")
            fields[key] = fields[key] == "true"
    if fields["schema"] != 1:
        raise ValueError(f"unsupported XC-response audit schema: {fields['schema']}")
    expected_diagnostic_code = {
        "cuda": "cuda_ready",
        "cpu_fallback": "cpu_fallback",
        "blocked": "blocked",
    }.get(str(fields["action"]))
    if expected_diagnostic_code is None:
        raise ValueError(f"unknown XC-response audit action: {fields['action']}")
    if fields["diagnostic_code"] != expected_diagnostic_code:
        raise ValueError(
            f"diagnostic_code {fields['diagnostic_code']} does not match action {fields['action']}"
        )
    for key in ("launch_elements", "total_bytes"):
        _parse_xc_response_audit_int_field(fields, key, allow_none=True)
    return fields


def xc_response_audit_call_site_guard(parsed: dict[str, Any]) -> dict[str, Any]:
    """Return execution flags for a parsed XC-response audit record.

    This guard is for downstream log/audit consumers that no longer have access
    to the raw CUDA handoff payload.  It therefore enforces semantic consistency
    between the parsed action, status, and public handoff summary before a
    consumer treats an audit line as evidence that the CUDA call site was ready.
    """

    if not isinstance(parsed, dict):
        raise ValueError("parsed XC-response audit record must be a dictionary")
    action = str(parsed.get("action", ""))
    status = _normalize_xc_response_status(parsed.get("status", XC_RESPONSE_STATUS_INVALID_INPUT))
    diagnostic_code = str(parsed.get("diagnostic_code", ""))
    expected_diagnostic_code = {
        "cuda": "cuda_ready",
        "cpu_fallback": "cpu_fallback",
        "blocked": "blocked",
    }.get(action)
    if expected_diagnostic_code is None:
        raise ValueError(f"unknown XC-response audit action: {action}")
    if diagnostic_code != expected_diagnostic_code:
        raise ValueError(f"diagnostic_code {diagnostic_code} does not match action {action}")
    if action == "cuda":
        if status != XC_RESPONSE_STATUS_OK or parsed.get("status_blocking"):
            raise ValueError("cuda audit record requires OK non-blocking status")
        if not parsed.get("handoff_valid") or parsed.get("launch_elements") is None or parsed.get("total_bytes") is None:
            raise ValueError("cuda audit record requires a valid handoff summary")
        launch_elements = parsed.get("launch_elements")
        total_bytes = parsed.get("total_bytes")
        if (
            not isinstance(launch_elements, int)
            or isinstance(launch_elements, bool)
            or not isinstance(total_bytes, int)
            or isinstance(total_bytes, bool)
        ):
            raise ValueError("cuda audit record launch_elements and total_bytes must be non-boolean integers")
        if launch_elements <= 0 or total_bytes <= 0:
            raise ValueError("cuda audit record requires positive launch_elements and total_bytes")
        dispatch = "cuda"
        execute_cuda = True
        execute_cpu_fallback = False
        raise_error = False
    elif action == "cpu_fallback":
        if status != XC_RESPONSE_STATUS_DISABLED or parsed.get("status_blocking"):
            raise ValueError("cpu_fallback audit record requires disabled non-blocking status")
        if parsed.get("handoff_valid") or parsed.get("launch_elements") is not None or parsed.get("total_bytes") is not None:
            raise ValueError("cpu_fallback audit record must not advertise CUDA handoff data")
        dispatch = "cpu_fallback"
        execute_cuda = False
        execute_cpu_fallback = True
        raise_error = False
    elif action == "blocked":
        if not parsed.get("status_blocking") and status in {XC_RESPONSE_STATUS_OK, XC_RESPONSE_STATUS_DISABLED}:
            raise ValueError("blocked audit record requires a blocking status or diagnostic")
        if parsed.get("handoff_valid") or parsed.get("launch_elements") is not None or parsed.get("total_bytes") is not None:
            raise ValueError("blocked audit record must not advertise CUDA handoff data")
        dispatch = "blocked"
        execute_cuda = False
        execute_cpu_fallback = False
        raise_error = True
    else:
        raise ValueError(f"unknown XC-response audit action: {action}")
    return {
        "dispatch": dispatch,
        "execute_cuda": execute_cuda,
        "execute_cpu_fallback": execute_cpu_fallback,
        "raise_error": raise_error,
        "diagnostic_code": diagnostic_code,
        "status": status,
        "reason": str(parsed.get("reason", "")),
    }


def xc_response_call_site_guard(plan: dict[str, Any]) -> dict[str, Any]:
    """Return safe execution flags plus audit metadata for an XC-response call site.

    The guard deliberately does not return the raw CUDA handoff payload.  A
    future production call site should use this for logging/control flow, while
    the lower-level dispatch plan remains the only object carrying validated
    launch metadata for a CUDA invocation.
    """

    audit = xc_response_dispatch_audit_record(plan)
    expected_flags = {
        "cuda": (True, False),
        "cpu_fallback": (False, True),
        "blocked": (False, False),
    }.get(audit["action"])
    actual_flags = (bool(plan.get("call_cuda", False)), bool(plan.get("use_cpu_fallback", False)))
    if expected_flags is not None and actual_flags != expected_flags:
        raise ValueError(
            "XC-response dispatch flags drift from audited action: "
            f"expected call_cuda/use_cpu_fallback={expected_flags}, got {actual_flags}"
        )
    execute_cuda = bool(plan.get("call_cuda", False) and audit["handoff_valid"] and audit["status"] == XC_RESPONSE_STATUS_OK)
    execute_cpu_fallback = bool(plan.get("use_cpu_fallback", False) and not execute_cuda)
    raise_error = not execute_cuda and not execute_cpu_fallback
    if execute_cuda:
        dispatch = "cuda"
    elif execute_cpu_fallback:
        dispatch = "cpu_fallback"
    else:
        dispatch = "blocked"
    return {
        "dispatch": dispatch,
        "execute_cuda": execute_cuda,
        "execute_cpu_fallback": execute_cpu_fallback,
        "raise_error": raise_error,
        "audit": audit,
        "audit_log_line": xc_response_dispatch_audit_log_line(audit),
    }


@dataclass(frozen=True)
class XcResponseGpuPlan:
    """Shape contract for the first CUDA XC-response scaffold.

    The current CUDA ABI smoke path treats density, kernel, and response as
    equally sized ``nbasis * nstate`` slots.  Keeping that contract in a pure
    Python helper gives the Fortran/CUDA wiring a stable manifest to validate
    before production quadrature/cache integration is added.
    """

    nbasis: int
    nstate: int

    def __post_init__(self) -> None:
        """Reject nonsensical buffer dimensions before CUDA allocation planning."""

        if not isinstance(self.nbasis, int) or isinstance(self.nbasis, bool):
            raise ValueError(f"nbasis must be a non-boolean integer for XC-response GPU buffers, got {self.nbasis}")
        if not isinstance(self.nstate, int) or isinstance(self.nstate, bool):
            raise ValueError(f"nstate must be a non-boolean integer for XC-response GPU buffers, got {self.nstate}")
        if self.nbasis <= 0:
            raise ValueError(f"nbasis must be positive for XC-response GPU buffers, got {self.nbasis}")
        if self.nstate <= 0:
            raise ValueError(f"nstate must be positive for XC-response GPU buffers, got {self.nstate}")
        max_fortran_c_int = 2_147_483_647
        if self.nbasis > max_fortran_c_int // self.nstate:
            raise ValueError(
                "XC-response GPU buffer dimension overflow: "
                f"nbasis*nstate={self.nbasis * self.nstate} exceeds Fortran C int limit"
            )

    @property
    def slots(self) -> int:
        """Return the number of scalar response slots per buffer."""

        return self.nbasis * self.nstate

    @property
    def total_elements(self) -> int:
        """Return total scalar elements required by the ABI smoke buffers."""

        return 3 * self.slots

    @property
    def bytes_per_element(self) -> int:
        """Return scalar storage width for the current float64-only ABI."""

        return 8

    @property
    def total_bytes(self) -> int:
        """Return total float64 workspace bytes required by the ABI smoke buffers."""

        return self.total_elements * self.bytes_per_element

    def buffer_manifest(self) -> list[dict[str, Any]]:
        """Return a stable ordered buffer contract for source-level tests."""

        return [
            {"name": "density", "elements": self.slots, "role": "input transition-density slots"},
            {"name": "kernel", "elements": self.slots, "role": "input XC kernel slots"},
            {"name": "response", "elements": self.slots, "role": "output contracted XC-response slots"},
        ]

    def buffer_byte_manifest(self) -> list[dict[str, Any]]:
        """Return stable per-buffer float64 byte counts for allocation planning."""

        return [
            {"name": item["name"], "bytes": item["elements"] * self.bytes_per_element, "role": item["role"]}
            for item in self.buffer_manifest()
        ]

    def fortran_buffer_table(self) -> list[tuple[int, str, int, str]]:
        """Return one-based buffer rows for future Fortran/C ABI checks."""

        return [
            (slot, item["name"], item["elements"], item["role"])
            for slot, item in enumerate(self.buffer_manifest(), start=1)
        ]

    def validate_fortran_buffer_table(
        self, table: list[tuple[int, str, int, str]] | tuple[tuple[int, str, int, str], ...]
    ) -> list[tuple[int, str, int, str]]:
        """Validate a Fortran/C ABI buffer table against this XC plan.

        This pure-Python guard catches slot/name/element/role drift before the
        experimental XC-response branch wires the Fortran table to real CUDA
        runtime allocation or contraction code.  Rows may arrive as lists after
        JSON decoding, so normalize row containers while preserving strict row
        arity and value checks.
        """

        expected = self.fortran_buffer_table()
        actual = []
        for idx, row in enumerate(table, start=1):
            if not isinstance(row, (list, tuple)) or len(row) != 4:
                raise ValueError(
                    f"XC-response buffer table row {idx} must contain slot, name, elements, and role"
                )
            slot, name, elements, role = row
            actual.append((slot, name, elements, role))
        if actual != expected:
            for expected_row, actual_row in zip(expected, actual):
                if actual_row != expected_row:
                    raise ValueError(
                        f"XC-response buffer table mismatch at slot {expected_row[0]}: "
                        f"expected {expected_row!r}, got {actual_row!r}"
                    )
            raise ValueError(
                f"XC-response buffer table length mismatch: expected {len(expected)}, got {len(actual)}"
            )
        return expected

    def consumer_layout_from_fortran_table(
        self, table: list[tuple[int, str, int, str]] | tuple[tuple[int, str, int, str], ...]
    ) -> dict[str, dict[str, int | str]]:
        """Return contiguous buffer offsets after validating the ABI table."""

        layout: dict[str, dict[str, int | str]] = {}
        offset = 0
        for slot, name, elements, role in self.validate_fortran_buffer_table(table):
            stop = offset + elements
            layout[name] = {
                "slot": slot,
                "elements": elements,
                "offset": offset,
                "stop": stop,
                "role": role,
            }
            offset = stop
        return layout

    def cuda_launch_shape(self, *, block_size: int = 128) -> dict[str, int]:
        """Return a deterministic launch-shape contract for the smoke kernel."""

        if not isinstance(block_size, int) or isinstance(block_size, bool) or block_size <= 0:
            raise ValueError(
                f"block_size must be a non-boolean positive integer for XC-response CUDA launch planning, got {block_size}"
            )
        return {
            "blocks": (self.slots + block_size - 1) // block_size,
            "threads_per_block": block_size,
            "elements": self.slots,
        }

    def runtime_handoff_payload(self, *, block_size: int = 128) -> dict[str, Any]:
        """Return a typed metadata payload for the future Python/Fortran call site."""

        table = self.fortran_buffer_table()
        return {
            "nbasis": self.nbasis,
            "nstate": self.nstate,
            "slots": self.slots,
            "total_elements": self.total_elements,
            "total_bytes": self.total_bytes,
            "launch": self.cuda_launch_shape(block_size=block_size),
            "fortran_table": table,
            "layout": self.consumer_layout_from_fortran_table(table),
        }


@dataclass(frozen=True)
class GpuConfig:
    """Normalized GPU runtime configuration."""

    enabled: bool = False
    backend: str = "cuda"
    target: Any = "metc"
    device: int = 0
    precision: str = "float64"
    fallback: str = "cpu"

    def __post_init__(self) -> None:
        """Normalize ``target`` into the locked multi-valued schema.

        A single requested target is preserved as a plain string for backward
        compatibility; multiple targets are stored as an immutable tuple.  Use
        :attr:`targets` for the always-normalized tuple and :meth:`has_target`
        for membership checks.
        """

        normalized = normalize_gpu_targets(self.target)
        coerced: Any = normalized[0] if len(normalized) == 1 else normalized
        object.__setattr__(self, "target", coerced)

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> "GpuConfig":
        """Create a normalized GPU config from an OpenQP config dictionary."""

        gpu = config.get("gpu", {}) or {}
        return cls(
            enabled=bool(gpu.get("enabled", False)),
            backend=str(gpu.get("backend", "cuda")).lower(),
            target=gpu.get("target", "metc"),
            device=int(gpu.get("device", 0)),
            precision=str(gpu.get("precision", "float64")).lower(),
            fallback=str(gpu.get("fallback", "cpu")).lower(),
        )

    @property
    def targets(self) -> tuple[str, ...]:
        """Return the normalized, immutable tuple of requested GPU targets."""

        if isinstance(self.target, str):
            return (self.target,)
        return tuple(self.target)

    def has_target(self, name: str) -> bool:
        """Return True when ``name`` is among the requested GPU targets."""

        return name in self.targets

    @property
    def targets_metc(self) -> bool:
        """Return True when METC contractions are among the requested targets."""

        return self.has_target("metc")

    @property
    def targets_xc_response(self) -> bool:
        """Return True when TDHF/TDDFT XC response is among the requested targets."""

        return self.has_target("xc_response")

    def supports_metc(self, config: dict[str, Any]) -> bool:
        """Return True when this config can use the scoped GPU METC path.

        METC is scoped to MRSF/UMRSF TDHF Davidson contractions.  Membership is
        used so a multi-target config that includes ``"metc"`` still qualifies.
        """

        input_section = config.get("input", {}) or {}
        tdhf_section = config.get("tdhf", {}) or {}
        method = str(input_section.get("method", "hf")).lower()
        td_type = str(tdhf_section.get("type", "rpa")).lower()
        return (
            self.enabled
            and self.backend == "cuda"
            and self.has_target("metc")
            and self.precision == "float64"
            and method == "tdhf"
            and td_type in {"mrsf", "umrsf"}
        )

    def xc_response_support_report(self, config: dict[str, Any]) -> dict[str, bool | str]:
        """Explain whether this config can use the scoped XC-response CUDA path."""

        input_section = config.get("input", {}) or {}
        tdhf_section = config.get("tdhf", {}) or {}
        method = str(input_section.get("method", "hf")).lower()
        functional = str(input_section.get("functional", "")).strip()
        td_type = str(tdhf_section.get("type", "rpa")).lower()

        if not self.enabled:
            return {"supported": False, "reason": "gpu.enabled is false"}
        if self.backend != "cuda":
            return {"supported": False, "reason": "gpu.backend must be cuda for XC-response GPU support"}
        if not self.has_target("xc_response"):
            return {"supported": False, "reason": "gpu.target must include xc_response"}
        if self.precision != "float64":
            return {"supported": False, "reason": "gpu.precision must be float64 for validated XC-response comparisons"}
        if method != "tdhf":
            return {"supported": False, "reason": "input.method must be tdhf for XC-response GPU support"}
        if td_type not in {"rpa", "tda"}:
            return {"supported": False, "reason": "tdhf.type must be rpa or tda for XC-response GPU support"}
        if not functional:
            return {"supported": False, "reason": "input.functional is required for XC-response GPU support"}
        return {"supported": True, "reason": "TDHF/TDDFT RPA/TDA XC-response CUDA path is enabled"}

    def xc_response_runtime_decision(
        self,
        config: dict[str, Any],
        *,
        nbasis: int,
        nstate: int,
        block_size: int = 128,
        runtime_status: int = XC_RESPONSE_STATUS_DISABLED,
    ) -> dict[str, Any]:
        """Package the safe XC-response runtime handoff decision.

        The branch is still scaffolding-only by default, so supported inputs
        produce a validated plan plus the current disabled/fallback ABI status
        instead of implying production CUDA execution.
        """

        support = self.xc_response_support_report(config)
        if not support["supported"]:
            status = xc_response_status_report(XC_RESPONSE_STATUS_INVALID_INPUT)
            return {
                "enabled": self.enabled,
                "supported": False,
                "reason": support["reason"],
                "status": status["status"],
                "status_message": status["message"],
                "fallback": status["fallback"],
                "blocking": status["blocking"],
                "fallback_policy": self.fallback,
                "status_report": status,
                "action": "blocked",
                "handoff": None,
                "plan": None,
            }

        try:
            plan = XcResponseGpuPlan(nbasis=nbasis, nstate=nstate)
        except ValueError as exc:
            message = str(exc)
            status_code = XC_RESPONSE_STATUS_OVERFLOW if "overflow" in message.lower() else XC_RESPONSE_STATUS_INVALID_INPUT
            status = xc_response_status_report(status_code)
            return {
                "enabled": self.enabled,
                "supported": True,
                "reason": message,
                "status": status["status"],
                "status_message": status["message"],
                "fallback": status["fallback"],
                "blocking": status["blocking"],
                "fallback_policy": self.fallback,
                "status_report": status,
                "action": "blocked",
                "handoff": None,
                "plan": None,
            }

        try:
            payload = plan.runtime_handoff_payload(block_size=block_size)
        except ValueError as exc:
            status = xc_response_status_report(XC_RESPONSE_STATUS_INVALID_INPUT)
            return {
                "enabled": self.enabled,
                "supported": True,
                "reason": str(exc),
                "status": status["status"],
                "status_message": status["message"],
                "fallback": status["fallback"],
                "blocking": status["blocking"],
                "fallback_policy": self.fallback,
                "status_report": status,
                "action": "blocked",
                "handoff": None,
                "plan": None,
            }

        status = xc_response_status_report(runtime_status)
        action = "cuda" if status["ok"] else "cpu_fallback" if status["fallback"] else "blocked"
        return {
            "enabled": self.enabled,
            "supported": True,
            "reason": support["reason"],
            "status": status["status"],
            "status_message": status["message"],
            "fallback": status["fallback"],
            "blocking": status["blocking"],
            "fallback_policy": self.fallback,
            "status_report": status,
            "action": action,
            "handoff": payload,
            "plan": {
                "elements": payload["launch"]["elements"],
                "blocks": payload["launch"]["blocks"],
                "threads_per_block": payload["launch"]["threads_per_block"],
                "total_bytes": payload["total_bytes"],
                "buffers": plan.buffer_byte_manifest(),
                "layout": payload["layout"],
            },
        }

    def xc_response_runtime_dispatch_plan(
        self,
        config: dict[str, Any],
        *,
        nbasis: int,
        nstate: int,
        block_size: int = 128,
        runtime_status: int = XC_RESPONSE_STATUS_DISABLED,
    ) -> dict[str, Any]:
        """Return the guarded XC-response runtime dispatch plan for call sites."""

        decision = self.xc_response_runtime_decision(
            config,
            nbasis=nbasis,
            nstate=nstate,
            block_size=block_size,
            runtime_status=runtime_status,
        )
        return xc_response_dispatch_plan(decision)

    def xc_response_runtime_dispatch_audit_record(
        self,
        config: dict[str, Any],
        *,
        nbasis: int,
        nstate: int,
        block_size: int = 128,
        runtime_status: int = XC_RESPONSE_STATUS_DISABLED,
    ) -> dict[str, Any]:
        """Return a compact call-site audit record without exposing handoff payloads."""

        plan = self.xc_response_runtime_dispatch_plan(
            config,
            nbasis=nbasis,
            nstate=nstate,
            block_size=block_size,
            runtime_status=runtime_status,
        )
        return xc_response_dispatch_audit_record(plan)

    def xc_response_runtime_call_site_guard(
        self,
        config: dict[str, Any],
        *,
        nbasis: int,
        nstate: int,
        block_size: int = 128,
        runtime_status: int = XC_RESPONSE_STATUS_DISABLED,
    ) -> dict[str, Any]:
        """Return final XC-response call-site control flags without raw CUDA payloads."""

        plan = self.xc_response_runtime_dispatch_plan(
            config,
            nbasis=nbasis,
            nstate=nstate,
            block_size=block_size,
            runtime_status=runtime_status,
        )
        return xc_response_call_site_guard(plan)

    def supports_xc_response(self, config: dict[str, Any]) -> bool:
        """Return True for the first scoped XC-response GPU validation path."""

        return bool(self.xc_response_support_report(config)["supported"])
