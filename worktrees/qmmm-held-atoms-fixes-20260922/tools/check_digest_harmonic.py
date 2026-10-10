#!/usr/bin/env python3
"""Rule 6: every grd2 two-particle-density digest must handle a spherical basis.

grd2 drives `get_density` with CARTESIAN shell extents. A digest that indexes
its densities with `basis%ao_offset` / `basis%naos`, which count the ACTUAL
AOs, therefore reads the wrong elements whenever the basis is spherical -- and
only then. In a Cartesian basis the two index spaces coincide, so the mistake
is completely invisible: it passes every example, every regression test and
every CI job that uses a Pople basis.

That is not hypothetical. The analytic MRSF NAC digest shipped without the
branch and was wrong for every spherical basis: with d functions (5 vs 6) it
read neighbouring AO elements, broke molecular symmetry, and amplified the
1e-13 run-to-run SCF reduction-order noise into an O(1) change in the answer;
with f functions (7 vs 10) it wrote past the end of the array and aborted. It
was added to a file whose OTHER digest handled the spherical case correctly,
so a check at file granularity would have missed it. This one is per type.

The gate is deliberately dumb: it asks only that the `get_density` bound to
each concrete `grd2_compute_data_t` extension mentions HARMONIC_ACTIVE. It
cannot tell a correct branch from an incorrect one. What it does guarantee is
that nobody adds a digest having never considered the question.
"""

import argparse
import pathlib
import re
import sys

BASE = "grd2_compute_data_t"
# Attributes may appear on either side of extends(...), so take the whole
# attribute list and pull the parent out of it.
TYPE_RE = re.compile(
    r"^[ \t]*type[ \t]*,(?P<attrs>[^:\n]*)::[ \t]*(?P<name>\w+)",
    re.IGNORECASE | re.MULTILINE)
EXTENDS_RE = re.compile(r"extends\s*\(\s*(?P<parent>\w+)\s*\)", re.IGNORECASE)
BINDING_RE = re.compile(
    r"^\s*procedure\s*(?:\([^)]*\))?\s*(?:,[^:]*)?::\s*get_density\s*=>\s*(?P<impl>\w+)",
    re.IGNORECASE | re.MULTILINE)


def _type_block(text, start):
    """Return the source of one derived-type definition, from `start`."""
    end = re.compile(r"^\s*end\s+type\b", re.IGNORECASE | re.MULTILINE)
    match = end.search(text, start)
    return text[start:match.end()] if match else text[start:]


def _procedure_body(text, name):
    """Return the body of `subroutine name`, or None if it is not here."""
    opener = re.compile(
        rf"^\s*(?:pure\s+|elemental\s+|recursive\s+)*subroutine\s+{re.escape(name)}\s*\(",
        re.IGNORECASE | re.MULTILINE)
    start = opener.search(text)
    if start is None:
        return None
    closer = re.compile(
        rf"^\s*end\s+subroutine\s+{re.escape(name)}\b", re.IGNORECASE | re.MULTILINE)
    stop = closer.search(text, start.end())
    return text[start.start():stop.end() if stop else len(text)]


def _derived_types(root):
    """Every declaration of a type transitively deriving from BASE.

    A digest that extends an intermediate subtype rather than BASE itself is
    still a digest. hf_gradient.F90 has exactly that shape -- an abstract
    grd2_hf_compute_data_t with grd2_rhf/uhf_compute_data_t under it, each
    binding its own get_density -- and matching only the literal base name left
    both of them unexamined.
    """
    # A name may legitimately appear more than once across the tree, so keep
    # every declaration rather than letting the last one seen win. A silently
    # shadowed record would make this gate examine the wrong procedure body,
    # which is the exact class of failure it exists to prevent.
    found = {}
    for path in sorted(pathlib.Path(root).rglob("*.F90")):
        text = path.read_text(errors="replace")
        for match in TYPE_RE.finditer(text):
            attrs = match.group("attrs")
            parent = EXTENDS_RE.search(attrs)
            if parent is None:
                continue
            found.setdefault(match.group("name").lower(), []).append({
                "name": match.group("name"),
                "parent": parent.group("parent").lower(),
                "abstract": "abstract" in attrs.lower(),
                "path": path,
                "text": text,
                "start": match.start(),
            })

    derived, changed = {BASE.lower()}, True
    while changed:
        changed = False
        for name, records in found.items():
            if name in derived:
                continue
            if any(r["parent"] in derived for r in records):
                derived.add(name)
                changed = True
    return [r for name, records in found.items() if name in derived
            for r in records]


def check(root):
    root = pathlib.Path(root)
    failures = []
    checked = 0
    records = sorted(_derived_types(root),
                     key=lambda r: (str(r["path"]), r["start"]))
    for record in records:
        if record["abstract"]:
            continue
        block = _type_block(record["text"], record["start"])
        binding = BINDING_RE.search(block)
        if binding is None:
            # Inherits its parent's digest; nothing of its own to check.
            continue
        impl = binding.group("impl")
        body = _procedure_body(record["text"], impl)
        rel = record["path"].relative_to(root)
        if body is None:
            failures.append(
                f"{rel}: type {record['name']} binds get_density to "
                f"{impl}, which is not defined in this file")
            continue
        checked += 1
        if "HARMONIC_ACTIVE" not in body:
            line = record["text"][:record["start"]].count("\n") + 1
            failures.append(
                f"{rel}:{line}: type {record['name']} digests a "
                f"two-particle density in {impl} without ever consulting "
                f"HARMONIC_ACTIVE.\n"
                f"    grd2 drives get_density with CARTESIAN shell "
                f"extents. Index the densities with basis%ao_offset / "
                f"basis%naos and the digest is silently wrong for every "
                f"spherical basis.\n"
                f"    See grd2_mrsf_compute_data_t_get_density and "
                f"grd2_mrsf_build_cart in source/modules/"
                f"tdhf_mrsf_gradient.F90 for the pattern.")
    return checked, failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", nargs="?", default="source",
                        help="directory to scan (default: source)")
    args = parser.parse_args()
    root = pathlib.Path(args.root).resolve()
    if not root.is_dir():
        print(f"grd2 spherical-basis gate (rule 6): no such directory {root}")
        return 2
    checked, failures = check(root)
    if failures:
        print("grd2 spherical-basis gate (rule 6): FAIL")
        for failure in failures:
            print(f"  {failure}")
        return 1
    print(f"grd2 spherical-basis gate (rule 6): PASS, {checked} digest(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
