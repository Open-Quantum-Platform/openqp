"""FCIDUMP export for OQP reference orbitals.

This delegates to the in-tree :mod:`oqp.quantum` package, which builds the
second-quantized Hamiltonian from OpenQP's *own* AO integrals (``OQP::Hcore``
and the native two-electron integrals ``OQP::ERI_AO``) in OpenQP's MO basis --
so the AO ordering/normalization is guaranteed consistent and the file is the
genuine OpenQP Hamiltonian (no foreign integral engine, nothing silently
mixed). It does not reimplement the writer.

:func:`verify_fcidump_fci` is an optional cross-check: read the FCIDUMP, run a
PySCF FCI, and compare to PySCF's native FCI for the same system (FCI is
invariant to the orbital choice, so this is a clean ground truth). It honours
the spin sector (``MS2``) recorded in the file, so open-shell dumps verify
correctly.
"""
import numpy as np
from oqp.periodic_table import ELEMENTS_NAME

__all__ = ["dump_fcidump", "verify_fcidump_fci", "build_pyscf_mol"]


def dump_fcidump(path, mol, tol=1.0e-12):
    """Write a FCIDUMP for the OQP reference MOs via :mod:`oqp.quantum`.

    Returns a small metadata dict (norb / nelec / ms2 / core energy)."""
    from oqp.quantum import from_openqp
    ham = from_openqp(mol)                 # OQP native Hcore + ERI in OQP MO basis
    ham.to_fcidump(path, tol=tol)
    return {
        "path": path,
        "norb": int(ham.n_orbitals),
        "nelec": int(ham.n_electrons),
        "ms2": int(ham.ms2),
        "core_energy": float(ham.core_energy),
        "engine": "oqp.quantum (native OQP integrals)",
    }


def build_pyscf_mol(mol):
    """Construct a PySCF Mole matching the OQP molecule/basis (for the FCI
    cross-check only)."""
    from pyscf import gto
    Z = np.asarray(mol.get_atoms(), dtype=int)
    xyz = np.asarray(mol.get_system(), dtype=float).reshape(-1, 3)   # Bohr
    atom = [[ELEMENTS_NAME[int(z)].strip(), tuple(xyz[i])] for i, z in enumerate(Z)]
    pm = gto.Mole()
    pm.atom = atom
    pm.unit = "Bohr"
    pm.basis = mol.config.get("input", {}).get("basis", "sto-3g")
    pm.charge = int(mol.config.get("input", {}).get("charge", 0) or 0)
    pm.spin = int(mol.config.get("scf", {}).get("multiplicity", 1) or 1) - 1
    pm.cart = True
    pm.build()
    return pm


def verify_fcidump_fci(path, mol):
    """Read FCIDUMP -> PySCF FCI, honouring MS2; compare to native PySCF FCI."""
    from pyscf import scf, fci, ao2mo
    from pyscf.tools import fcidump as fcidump_tools

    d = fcidump_tools.read(path)
    norb = d["NORB"]
    nelec_tot = d["NELEC"]
    ms2 = d.get("MS2", 0)
    # split the total electron count into the (nalpha, nbeta) sector recorded
    # in the file so open-shell dumps are verified in the right spin state.
    na = (nelec_tot + ms2) // 2
    nb = (nelec_tot - ms2) // 2
    h1 = d["H1"]
    h2 = ao2mo.restore(1, d["H2"], norb)
    ecore = d["ECORE"]
    e_dump, _ = fci.direct_spin1.kernel(h1, h2, norb, (na, nb), ecore=ecore)

    pm = build_pyscf_mol(mol)
    mf = (scf.RHF(pm) if pm.spin == 0 else scf.ROHF(pm)).run(verbose=0)
    e_ref = fci.FCI(mf).kernel()[0]
    return {"e_fcidump": float(e_dump), "e_pyscf_fci": float(e_ref),
            "diff": float(abs(e_dump - e_ref)), "nelec_ab": (na, nb)}
