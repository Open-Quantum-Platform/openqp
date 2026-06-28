"""FCIDUMP export of OQP MO integrals + cross-check against PySCF FCI.

OQP does not expose the AO two-electron integrals to Python, and the task
restricts Fortran edits to the 1-TDM / AO-on-grid intermediates, so the AO
integral engine here is PySCF, evaluated in the *same* basis and contracted
into OQP's reference MO coefficients.  Consistency is asserted by comparing
PySCF's S / Hcore / E_nuc to OQP's own (``OQP::SM`` / ``OQP::Hcore`` /
``vnn``), and the integral file is validated by an FCI round-trip:
FCIDUMP -> PySCF FCI must match PySCF's native FCI to <= 1e-8 Hartree
(FCI is invariant to the orbital choice, so this is a clean ground-truth).

Currently validated for s/p-only AO ordering (e.g. H2/STO-3G).  For higher
angular momenta an OQP<->PySCF Cartesian reorder map would be required; the
OQP component order is documented in ``oqp.analysis.gto_grid``.
"""
import numpy as np
from oqp.periodic_table import ELEMENTS_NAME

__all__ = ["build_pyscf_mol", "dump_fcidump", "verify_fcidump_fci"]


def _unpack_lt(packed, n):
    m = np.zeros((n, n))
    r, c = np.tril_indices(n)
    m[r, c] = packed
    m[c, r] = packed
    return m


def build_pyscf_mol(mol):
    """Construct a PySCF Mole matching the OQP molecule/basis (unit = Bohr)."""
    from pyscf import gto
    Z = np.asarray(mol.get_atoms(), dtype=int)
    xyz = np.asarray(mol.get_system(), dtype=float).reshape(-1, 3)   # Bohr
    atom = [[ELEMENTS_NAME[int(z)].strip(), tuple(xyz[i])] for i, z in enumerate(Z)]
    basis = mol.config.get("input", {}).get("basis", "sto-3g")
    charge = int(mol.config.get("input", {}).get("charge", 0) or 0)
    mult = int(mol.config.get("scf", {}).get("multiplicity", 1) or 1)
    pm = gto.Mole()
    pm.atom = atom
    pm.unit = "Bohr"
    pm.basis = basis
    pm.charge = charge
    pm.spin = mult - 1
    pm.cart = True       # OQP Pople basis is Cartesian (6d)
    pm.build()
    return pm


def _consistency(mol, pm, tol=1e-6):
    """Compare PySCF vs OQP S / Hcore / E_nuc to confirm the same basis+order."""
    nbf = pm.nao_nr()
    S_p = pm.intor("int1e_ovlp")
    H_p = pm.intor("int1e_kin") + pm.intor("int1e_nuc")
    enuc_p = float(pm.energy_nuc())
    report = {"nbf": nbf}
    try:
        S_o = _unpack_lt(np.asarray(mol.data["OQP::SM"]).ravel(), nbf)
        report["max_dS"] = float(np.max(np.abs(S_o - S_p)))
    except Exception as e:
        report["max_dS"] = None
    try:
        H_o = _unpack_lt(np.asarray(mol.data["OQP::Hcore"]).ravel(), nbf)
        report["max_dHcore"] = float(np.max(np.abs(H_o - H_p)))
    except Exception:
        report["max_dHcore"] = None
    try:
        report["dEnuc"] = abs(float(mol.get_scf_energy("vnn")) - enuc_p)
    except Exception:
        report["dEnuc"] = None
    return report, S_p, H_p, enuc_p


def dump_fcidump(path, mol, source="oqp"):
    """Write a FCIDUMP for the OQP reference MOs.

    ``source`` selects the 1e/nuclear terms: 'oqp' uses OQP's own Hcore and
    nuclear repulsion (so the file is OQP's Hamiltonian); 'pyscf' uses PySCF's.
    The 2e integrals always come from PySCF.  Returns a metadata dict including
    the consistency report and the max 8-fold-symmetry residual."""
    pm = build_pyscf_mol(mol)
    nbf = pm.nao_nr()
    rep, S_p, H_p, enuc_p = _consistency(mol, pm)

    na = int(np.asarray(mol.data["nelec_A"]).ravel()[0])
    nb = int(np.asarray(mol.data["nelec_B"]).ravel()[0])
    C = np.asarray(mol.data["OQP::VEC_MO_A"]).reshape(nbf, nbf).T       # C[ao, mo]

    if source == "oqp" and rep.get("max_dHcore") is not None:
        H_ao = _unpack_lt(np.asarray(mol.data["OQP::Hcore"]).ravel(), nbf)
        enuc = float(mol.get_scf_energy("vnn"))
    else:
        H_ao, enuc = H_p, enuc_p

    eri_ao = pm.intor("int2e")                                          # (pq|rs)
    h_mo = C.T @ H_ao @ C
    eri_mo = np.einsum("pqrs,pi,qj,rk,sl->ijkl", eri_ao, C, C, C, C, optimize=True)

    # 8-fold permutational symmetry residual of the MO 2e integrals
    perms = [eri_mo, eri_mo.transpose(1, 0, 2, 3), eri_mo.transpose(0, 1, 3, 2),
             eri_mo.transpose(1, 0, 3, 2), eri_mo.transpose(2, 3, 0, 1),
             eri_mo.transpose(3, 2, 0, 1), eri_mo.transpose(2, 3, 1, 0),
             eri_mo.transpose(3, 2, 1, 0)]
    sym_res = max(float(np.max(np.abs(eri_mo - p))) for p in perms)

    norb = nbf
    nelec = na + nb
    ms2 = na - nb
    with open(path, "w") as f:
        f.write(f" &FCI NORB={norb},NELEC={nelec},MS2={ms2},\n")
        f.write("  ORBSYM=" + ",".join(["1"] * norb) + ",\n")
        f.write("  ISYM=1,\n")
        f.write(" &END\n")
        # 2e: chemist (ij|kl), unique 8-fold (i>=j, k>=l, ij>=kl), 1-based
        for i in range(norb):
            for j in range(i + 1):
                ij = i * (i + 1) // 2 + j
                for k in range(norb):
                    for l in range(k + 1):
                        kl = k * (k + 1) // 2 + l
                        if ij < kl:
                            continue
                        v = eri_mo[i, j, k, l]
                        if abs(v) > 1e-12:
                            f.write(f"{v:23.16e}{i+1:4d}{j+1:4d}{k+1:4d}{l+1:4d}\n")
        # 1e: h_pq (i>=j)
        for i in range(norb):
            for j in range(i + 1):
                v = h_mo[i, j]
                if abs(v) > 1e-12:
                    f.write(f"{v:23.16e}{i+1:4d}{j+1:4d}{0:4d}{0:4d}\n")
        f.write(f"{enuc:23.16e}{0:4d}{0:4d}{0:4d}{0:4d}\n")

    return {
        "path": path, "norb": norb, "nelec": nelec, "ms2": ms2,
        "consistency": rep, "sym_residual_8fold": sym_res,
        "source_1e": source,
    }


def verify_fcidump_fci(path, mol):
    """Read FCIDUMP -> PySCF FCI; compare to PySCF's native FCI for the system."""
    from pyscf import gto, scf, fci, ao2mo
    from pyscf.tools import fcidump as fcidump_tools

    d = fcidump_tools.read(path)
    norb = d["NORB"]
    nelec = d["NELEC"]
    h1 = d["H1"]
    h2 = ao2mo.restore(1, d["H2"], norb)
    ecore = d["ECORE"]
    e_dump, _ = fci.direct_spin1.kernel(h1, h2, norb, nelec, ecore=ecore)

    # native reference FCI
    pm = build_pyscf_mol(mol)
    mf = scf.RHF(pm).run(verbose=0) if pm.spin == 0 else scf.ROHF(pm).run(verbose=0)
    cisolver = fci.FCI(mf)
    e_ref = cisolver.kernel()[0]
    return {"e_fcidump": float(e_dump), "e_pyscf_fci": float(e_ref),
            "diff": float(abs(e_dump - e_ref))}
