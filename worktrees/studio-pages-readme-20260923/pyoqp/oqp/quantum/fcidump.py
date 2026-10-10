"""Read and write the FCIDUMP interchange format.

FCIDUMP is the de-facto standard exchange format for the second-quantized
molecular electronic Hamiltonian

    H = E_core
        + sum_{pq,s}        h_{pq}      a^+_{p,s} a_{q,s}
        + 1/2 sum_{pqrs,st} (pq|rs)     a^+_{p,s} a^+_{r,t} a_{s,t} a_{q,s}

expressed over a set of (real, orthonormal) molecular spatial orbitals. It is
consumed directly by Qiskit Nature, OpenFermion (``MolecularData`` via
``from_pyscf``/FCIDUMP loaders), PySCF (``pyscf.tools.fcidump``), Block2/DMRG,
and most active-space / quantum-computing front ends, which makes it the
shortest path from an OpenQP mean-field calculation to a quantum-computing
workflow.

This module has no dependency on the compiled ``oqp`` extension.
"""

import numpy as np

# Default threshold below which integrals are treated as numerical zero and
# omitted from the file (matches common quantum-chemistry practice).
DEFAULT_TOLERANCE = 1.0e-12


def _iter_unique_2e(h2, tol):
    """Yield ``(value, p, q, r, s)`` for the 8-fold-unique two-electron block.

    Indices are 0-based here; the writer converts to 1-based on output.
    ``h2`` is in chemist notation ``(pq|rs)``.
    """
    norb = h2.shape[0]
    for p in range(norb):
        for q in range(p + 1):
            for r in range(norb):
                for s in range(r + 1):
                    # Restrict to the canonical pair ordering pq >= rs so each
                    # symmetry-equivalent integral is emitted exactly once.
                    if (p * (p + 1) // 2 + q) < (r * (r + 1) // 2 + s):
                        continue
                    val = h2[p, q, r, s]
                    if abs(val) > tol:
                        yield val, p, q, r, s


def write_fcidump(filename, h1, h2, ecore, n_electrons, ms2=0,
                  orbsym=None, isym=1, tol=DEFAULT_TOLERANCE):
    """Write a FCIDUMP file.

    Parameters
    ----------
    filename : str
        Output path.
    h1 : array_like, shape (norb, norb)
        One-electron MO integrals ``h_pq``.
    h2 : array_like, shape (norb, norb, norb, norb)
        Two-electron MO integrals ``(pq|rs)`` in chemist notation.
    ecore : float
        Scalar (core / nuclear-repulsion + frozen-core) energy.
    n_electrons : int
        Number of correlated electrons.
    ms2 : int
        ``2 * S_z`` (number of alpha minus beta electrons). Default 0.
    orbsym : sequence of int, optional
        Per-orbital symmetry labels. Defaults to all-totally-symmetric (1).
    isym : int
        Spatial symmetry of the target state.
    tol : float
        Integrals with magnitude at or below this are omitted.
    """
    h1 = np.asarray(h1, dtype=float)
    h2 = np.asarray(h2, dtype=float)
    norb = h1.shape[0]
    if h1.shape != (norb, norb):
        raise ValueError("h1 must be square (norb, norb)")
    if h2.shape != (norb, norb, norb, norb):
        raise ValueError("h2 must have shape (norb, norb, norb, norb)")
    if orbsym is None:
        orbsym = [1] * norb
    if len(orbsym) != norb:
        raise ValueError("orbsym must have length norb")

    with open(filename, "w") as fh:
        fh.write(f" &FCI NORB={norb},NELEC={int(n_electrons)},MS2={int(ms2)},\n")
        fh.write("  ORBSYM=" + ",".join(str(int(s)) for s in orbsym) + ",\n")
        fh.write(f"  ISYM={int(isym)},\n")
        fh.write(" &END\n")

        fmt = "{:28.20E}{:4d}{:4d}{:4d}{:4d}\n"
        # Two-electron integrals (8-fold unique), 1-based indices.
        for val, p, q, r, s in _iter_unique_2e(h2, tol):
            fh.write(fmt.format(val, p + 1, q + 1, r + 1, s + 1))
        # One-electron integrals (lower triangle), k = l = 0.
        for p in range(norb):
            for q in range(p + 1):
                val = h1[p, q]
                if abs(val) > tol:
                    fh.write(fmt.format(val, p + 1, q + 1, 0, 0))
        # Core energy, all indices zero.
        fh.write(fmt.format(float(ecore), 0, 0, 0, 0))


def read_fcidump(filename):
    """Read a FCIDUMP file produced by :func:`write_fcidump` (or PySCF).

    Returns a dict with keys ``norb``, ``nelec``, ``ms2``, ``isym``,
    ``orbsym``, ``h1`` (norb, norb), ``h2`` (norb^4, chemist notation, fully
    symmetrized) and ``ecore``.
    """
    header = {}
    with open(filename) as fh:
        lines = fh.readlines()

    # Parse the namelist header up to &END (or /).
    body_start = 0
    header_text = ""
    for i, line in enumerate(lines):
        stripped = line.strip()
        header_text += " " + stripped
        if stripped.upper().endswith("&END") or stripped == "/" \
                or stripped.upper() == "&END":
            body_start = i + 1
            break

    # Extract key=value tokens from the namelist.
    cleaned = header_text.replace("&FCI", " ").replace("&END", " ")
    cleaned = cleaned.replace("/", " ")
    tokens = [t for t in cleaned.replace("\n", " ").split() if t]
    joined = " ".join(tokens)
    import re
    for key in ("NORB", "NELEC", "MS2", "ISYM"):
        m = re.search(rf"{key}\s*=\s*(-?\d+)", joined, re.IGNORECASE)
        if m:
            header[key.lower()] = int(m.group(1))
    m = re.search(r"ORBSYM\s*=\s*([0-9,\s]+)", joined, re.IGNORECASE)
    if m:
        header["orbsym"] = [int(x) for x in m.group(1).split(",") if x.strip()]

    norb = header["norb"]
    h1 = np.zeros((norb, norb))
    h2 = np.zeros((norb, norb, norb, norb))
    ecore = 0.0

    for line in lines[body_start:]:
        parts = line.split()
        if len(parts) != 5:
            continue
        val = float(parts[0])
        p, q, r, s = (int(x) for x in parts[1:])
        if p == 0 and q == 0 and r == 0 and s == 0:
            ecore = val
        elif r == 0 and s == 0:
            h1[p - 1, q - 1] = val
            h1[q - 1, p - 1] = val
        else:
            i, j, k, l = p - 1, q - 1, r - 1, s - 1
            # Restore all 8 permutation symmetries of (ij|kl).
            for a, b, c, d in (
                (i, j, k, l), (j, i, k, l), (i, j, l, k), (j, i, l, k),
                (k, l, i, j), (l, k, i, j), (k, l, j, i), (l, k, j, i),
            ):
                h2[a, b, c, d] = val

    return {
        "norb": norb,
        "nelec": header.get("nelec"),
        "ms2": header.get("ms2", 0),
        "isym": header.get("isym", 1),
        "orbsym": header.get("orbsym", [1] * norb),
        "h1": h1,
        "h2": h2,
        "ecore": ecore,
    }
