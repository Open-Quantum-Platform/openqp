"""Tolerance harness: diff an OQP result against a parsed external result."""
import numpy as np

__all__ = ["compare_results", "format_table"]


def _delta(a, b):
    """Return ``(delta, note)``.

    ``note`` is non-empty when the two quantities cannot be compared
    element-for-element (one missing, or arrays of unequal length); such cases
    are treated as failures by :func:`compare_results` rather than silently
    comparing only the shared prefix."""
    if a is None or b is None:
        return None, "missing"
    if np.isscalar(a) and np.isscalar(b):
        return abs(float(a) - float(b)), ""
    a = np.atleast_1d(np.asarray(a, dtype=float))
    b = np.atleast_1d(np.asarray(b, dtype=float))
    note = "" if a.size == b.size else f"length {a.size} vs {b.size}"
    n = min(a.size, b.size)
    if n == 0:
        return None, (note or "empty")
    return float(np.max(np.abs(a[:n] - b[:n]))), note


def compare_results(ref, other, tolerances, ref_label="OQP", other_label="external"):
    """Compare two normalized result dicts under per-key tolerances.

    Returns ``(rows, all_pass)`` where each row is
    ``(quantity, delta, tolerance, status, note)``. A length mismatch (e.g. the
    external code reported fewer states) is a FAIL, not a prefix-only PASS."""
    rows = []
    all_pass = True
    for key, tol in tolerances.items():
        d, note = _delta(ref.get(key), other.get(key))
        if d is None:
            rows.append((key, None, tol, "N/A", note))
            continue
        if note:                       # length mismatch -> fail outright
            all_pass = False
            rows.append((key, d, tol, "FAIL", note))
            continue
        ok = d <= tol
        all_pass &= ok
        rows.append((key, d, tol, "PASS" if ok else "FAIL", ""))
    return rows, all_pass


def format_table(rows, ref_label="OQP", other_label="external", title=None):
    lines = []
    if title:
        lines.append(title)
    header = f"{'quantity':30s} {'|Δ| ('+ref_label+' vs '+other_label+')':26s} {'tol':>10s}  status  note"
    lines.append(header)
    lines.append("-" * len(header))
    for row in rows:
        key, d, tol, status = row[0], row[1], row[2], row[3]
        note = row[4] if len(row) > 4 else ""
        ds = "   N/A   " if d is None else f"{d:.3e}"
        lines.append(f"{key:30s} {ds:26s} {tol:10.1e}  {status:6s}  {note}")
    return "\n".join(lines)
