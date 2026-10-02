"""One definition of the QM region for a QM/MM job.

Two parts of the code need the QM atoms and count them differently:
``[qmmm] qm_atoms`` is a 0-based selection in the OpenMM topology, and the
1-based indices after the PDB name in ``[input] system`` build the QM
molecule.  Asking the user for both invites a silent mismatch, and leaving
out the second used to stop the run with "Atom list not defined!".

:func:`reconcile_qm_selection` makes either one sufficient: the missing list
is derived from the given one, and when both are given they must name the
same atoms.  Pure Python; no OpenMM, no OpenQP runtime.
"""
import re


def parse_index_list(value):
    """Indices from ``"0-3, 8 9"``, ``"9 10 17-19"`` or a list of ints."""
    if value is None:
        return []
    if not isinstance(value, (str, bytes)):
        # a list, a tuple, a NumPy array (a selection built in a script) or a
        # single integer -- anything but text is taken element by element
        try:
            return [int(v) for v in value]
        except TypeError:
            return [int(value)]
    out = []
    for item in re.split(r"[,\s]+", str(value).strip()):
        if not item:
            continue
        if "-" in item[1:]:
            first, last = item.split("-", 1)
            out.extend(range(int(first), int(last) + 1))
        else:
            out.append(int(item))
    return out


def split_pdb_reference(system):
    """``"ala.pdb 9 10"`` -> ``("ala.pdb", [9, 10])``; None when the geometry
    is not a PDB reference (an xyz file or an inline atom table)."""
    text = "" if system is None else str(system)
    if "\n" in text.strip():
        return None
    reference = text.strip()
    end = reference.lower().find(".pdb")
    if end < 0:
        return None
    return reference[:end + 4].strip(), parse_index_list(reference[end + 4:])


def reconcile_qm_selection(system, qm_atoms, pdb_file=None):
    """Return ``(system, qm_atoms)`` with both QM-atom lists filled in.

    * only ``qm_atoms`` (0-based)        -> the 1-based list is appended to the
      PDB reference; with no ``system`` at all the PDB is ``pdb_file``;
    * only the indices in ``system``     -> ``qm_atoms`` is derived from them;
    * both                               -> they must be the same atoms;
    * neither, or a non-PDB geometry     -> returned unchanged.
    """
    zero_based = parse_index_list(qm_atoms)
    reference = split_pdb_reference(system)
    if reference is None:
        blank = system is None or not str(system).strip()
        if blank and zero_based and pdb_file and str(pdb_file).strip():
            reference = (str(pdb_file).strip(), [])
        else:
            return system, qm_atoms
    path, one_based = reference

    if zero_based and not one_based:
        derived = " ".join(str(i + 1) for i in sorted(set(zero_based)))
        return f"{path} {derived}", qm_atoms
    if one_based and not zero_based:
        if min(one_based) < 1:
            raise ValueError(
                "QM atoms after the PDB name are 1-based; index 0 is not an atom")
        return system, " ".join(str(i - 1) for i in sorted(set(one_based)))
    if one_based and zero_based:
        if sorted(set(i - 1 for i in one_based)) != sorted(set(zero_based)):
            raise ValueError(
                "The QM region is defined twice and the two lists differ: "
                f"[qmmm] qm_atoms = {sorted(set(zero_based))} (0-based) but the "
                f"indices after {path} select "
                f"{sorted(set(i - 1 for i in one_based))} (given 1-based). "
                "Give qm_atoms alone; the QM molecule is built from it.")
    return system, qm_atoms
