"""Borrowed CASSCF descriptors must belong to each procedure invocation.

Declaration initialization gives a Fortran local implicit SAVE storage. This
source guard catches that cause independently of compiler and scheduling; it
does not establish reentrancy of the complete native driver.
"""
from pathlib import Path
import re

import pytest


SOURCE = Path(__file__).resolve().parents[1] / "source/modules/casscf_driver.F90"


def _check_borrowed_descriptors(source, routine, pointers):
    match = re.search(
        rf"^\s*(?:function|subroutine) {routine}\(.*?"
        rf"^\s*end (?:function|subroutine) {routine}\b",
        source, re.M | re.S | re.I,
    )
    assert match, f"missing {routine}"
    # Join continued statements and ignore comments before checking attributes.
    text = re.sub(r"!.*", "", match.group()).lower()
    text = re.sub(r"&\s*\n\s*&?", " ", text)
    statements = [line.strip() for line in text.splitlines() if line.strip()]
    assert not any(re.match(r"save\b", line) for line in statements), (
        f"{routine}: borrowed locals must not have explicit SAVE"
    )
    initialization = set()
    for line in statements:
        nullify = re.fullmatch(r"nullify\s*\(([^)]*)\)", line)
        if nullify:
            initialization.update(name.strip() for name in nullify[1].split(","))
        if re.match(r"status\s*=", line):
            break
    for pointer in pointers:
        declarations = [line for line in statements if "::" in line
                        and re.search(rf"\b{pointer}\s*\(", line)]
        assert len(declarations) == 1, f"{routine}: missing {pointer} declaration"
        declaration = declarations[0]
        attributes, _ = declaration.split("::", 1)
        assert re.search(r"\bpointer\b", attributes)
        assert not re.search(r"\bsave\b", attributes)
        assert "=>" not in declaration, (
            f"{routine}::{pointer}: declaration initialization implies SAVE"
        )
        assert pointer in initialization, (
            f"{routine}::{pointer}: initialize before status/validation work"
        )


@pytest.mark.parametrize("routine,pointers", [
    ("casscf_energy", ("mo_p",)),
    ("cas_context_init", ("hcore_p", "eri_p", "mo_p")),
])
def test_borrowed_pointer_descriptors_are_per_call(routine, pointers):
    _check_borrowed_descriptors(SOURCE.read_text(), routine, pointers)
