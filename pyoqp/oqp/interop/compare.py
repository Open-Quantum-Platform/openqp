"""Tolerance harness: diff an OQP result against a parsed external result."""
import numpy as np

__all__ = ["compare_results", "format_table"]


def _delta(a, b):
    """Max abs difference over the common length (arrays) or scalar diff."""
    if a is None or b is None:
        return None
    if np.isscalar(a) and np.isscalar(b):
        return abs(float(a) - float(b))
    a = np.atleast_1d(np.asarray(a, float))
    b = np.atleast_1d(np.asarray(b, float))
    n = min(a.size, b.size)
    if n == 0:
        return None
    return float(np.max(np.abs(a[:n] - b[:n])))


def compare_results(ref, other, tolerances, ref_label="OQP", other_label="external"):
    """Compare two normalized result dicts under per-key tolerances.

    Returns (rows, all_pass) where each row is
    (quantity, delta, tolerance, status)."""
    rows = []
    all_pass = True
    for key, tol in tolerances.items():
        d = _delta(ref.get(key), other.get(key))
        if d is None:
            rows.append((key, None, tol, "N/A"))
            continue
        ok = d <= tol
        all_pass &= ok
        rows.append((key, d, tol, "PASS" if ok else "FAIL"))
    return rows, all_pass


def format_table(rows, ref_label="OQP", other_label="external", title=None):
    lines = []
    if title:
        lines.append(title)
    header = f"{'quantity':32s} {'|Δ| ('+ref_label+' vs '+other_label+')':28s} {'tol':>10s}  status"
    lines.append(header)
    lines.append("-" * len(header))
    for key, d, tol, status in rows:
        ds = "   N/A   " if d is None else f"{d:.3e}"
        lines.append(f"{key:32s} {ds:28s} {tol:10.1e}  {status}")
    return "\n".join(lines)
