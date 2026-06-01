"""OQP external quantum chemical program"""
import os
import sys
import json
import warnings
from pathlib import Path

import numpy as np
import oqp
from oqp.utils.constants import ANGSTROM_TO_BOHR
from oqp.periodic_table import SYMBOL_MAP, ELEMENTS_NAME

try:
    from pyscf import gto, scf, dft, lib

except ModuleNotFoundError:
    gto = None
    scf = None
    dft = None
    lib = None

pyscf_functional = {
    "pbe0": "pbe0",
    "b3lyp": "b3lyp",
    "bhhlyp": "bhandhlyp",
    "cam-b3lyp": "camb3lyp",
    "camb3lyp": "camb3lyp",
}


def _json_array(value):
    """Convert NumPy/scalar values into JSON-serializable Python values."""

    return np.asarray(value).tolist()


def _spin_pair(value):
    """Return alpha/beta arrays from PySCF RHF/ROHF/UHF-style data."""

    if isinstance(value, (tuple, list)) and len(value) == 2:
        return np.asarray(value[0]), np.asarray(value[1])
    array = np.asarray(value)
    if array.ndim >= 3 and array.shape[0] == 2:
        return array[0], array[1]
    return array, array


def _pack_lower_triangle(matrix):
    """Pack a square matrix in OpenQP's lower-triangular restart layout."""

    mat = np.asarray(matrix)
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise ValueError(f"Expected a square density matrix, got shape={mat.shape}")
    idx = np.tril_indices(mat.shape[0])
    return mat[idx].tolist()


def _density_pair_from_mf(mf):
    """Return alpha/beta AO density matrices from a PySCF mean-field object."""

    mo_coeff = mf.mo_coeff
    mo_occ = getattr(mf, 'mo_occ', None)
    if isinstance(mo_coeff, (tuple, list)) and len(mo_coeff) == 2:
        return _spin_pair(mf.make_rdm1())
    if isinstance(mo_occ, (tuple, list)) and len(mo_occ) == 2:
        return _spin_pair(mf.make_rdm1())

    coeff = np.asarray(mo_coeff)
    occ = np.asarray(mo_occ) if mo_occ is not None else None
    if coeff.ndim == 2 and occ is not None and occ.ndim == 1:
        occ_a = np.where(occ > 0.0, 1.0, 0.0)
        occ_b = np.where(occ > 1.0, 1.0, 0.0)
        dm_a = np.dot(coeff * occ_a, coeff.T)
        dm_b = np.dot(coeff * occ_b, coeff.T)
        return dm_a, dm_b

    return _spin_pair(mf.make_rdm1())


def _pyscf_total_energy(mf):
    """Read the total energy from a PySCF mean-field object if available."""

    try:
        return float(mf.energy_tot())
    except Exception:
        pass
    try:
        return float(mf.e_tot)
    except Exception:
        pass
    try:
        return float(mf.mol.energy_tot())
    except Exception:
        return 0.0


def export_pyscf_guess_to_openqp_json(mf, mol, filename='pyscf_wfn.json'):
    """Export a PySCF mean-field object as an OpenQP JSON restart guess.

    This native exporter replaces the previous MOKIT ``py2openqp`` bridge for
    PySCF-backed guess modes. It writes the OpenQP fields consumed by
    ``Molecule.load_data``: MO coefficients, orbital energies, packed density
    matrices, atoms, coordinates, and basis/library metadata.
    """

    filename = Path(filename)
    coeff_a, coeff_b = _spin_pair(mf.mo_coeff)
    energy_a, energy_b = _spin_pair(mf.mo_energy)
    dm_a, dm_b = _density_pair_from_mf(mf)

    data = {
        'OQP::DM_A': _pack_lower_triangle(dm_a),
        'OQP::DM_B': _pack_lower_triangle(dm_b),
        'OQP::E_MO_A': _json_array(energy_a),
        'OQP::E_MO_B': _json_array(energy_b),
        'OQP::VEC_MO_A': _json_array(coeff_a),
        'OQP::VEC_MO_B': _json_array(coeff_b),
        'atoms': _json_array(mol.get_atoms()),
        'coord': _json_array(mol.get_system()),
        'energy': _pyscf_total_energy(mf),
        'td_energies': [0.0],
        'grad': [],
        'nac': [],
        'soc': [],
        'hess': [],
        'json': {
            'scf_type': mol.config['scf']['type'],
            'basis': mol.config['input']['basis'],
            'library': mol.config['input']['library'],
        },
    }

    with open(filename, 'w') as outdata:
        json.dump(data, outdata, indent=2)
    return filename


def load_pyscf_guess(mf, mol, filename='pyscf_wfn.json'):
    """Export a PySCF guess, load it into OpenQP, and project the basis."""

    target_basis = mol.config['input']['basis']
    target_library = mol.config['input']['library']
    guess_file = export_pyscf_guess_to_openqp_json(mf, mol, filename)
    mol.config['guess']['file'] = str(guess_file)
    mol.load_data()

    mol.data.set_scf_active_basis(1)
    mol.config['input']['basis'] = mol.config['json']['basis']
    mol.config['input']['library'] = mol.config['json']['library']
    oqp.library.set_basis(mol)

    mol.data.set_scf_active_basis(0)
    mol.config['input']['basis'] = target_basis
    mol.config['input']['library'] = target_library
    oqp.library.project_basis(mol)


def guess_from_pyscf(mol):
    """Set up pySCF calculation"""

    if gto is None:
        print(f'\n   PyOQP cannot find PySCF\n')
        ## pyscf run mpi with a different lib, which might cause conflict with mpi4py, disabling mpi now for future dev
        sys.exit(1)

    if mol.usempi:
        print(f'\n   PyOQP cannot run pySCF with MPI\n')
        ## pyscf run mpi with a different lib, which might cause conflict with mpi4py, disabling mpi now for future dev
        sys.exit(1)

    atoms = []
    coord = mol.get_system().reshape((-1, 3)) * ANGSTROM_TO_BOHR
    for n, at in enumerate(mol.get_atoms()):
        atoms.append([ELEMENTS_NAME[SYMBOL_MAP[int(at)]], coord[n]])

    mole = gto.Mole()
    mole.atom = atoms
    mole.unit = 'Bohr'
    mole.basis = mol.config["input"]["basis"]
    mole.charge = mol.config["input"]["charge"]
    mole.spin = mol.config["scf"]["multiplicity"] - 1
    mole.output = '%s.pyscf' % mol.project_name
    mole.verbose = 4
    mole.build(cart=True)
    functional = mol.config["input"]["functional"]

    if mol.config["scf"]["type"] == 'rohf':
        mf = dft.ROKS(mole)
    elif mol.config["scf"]["type"] == 'uhf':
        mf = dft.UKS(mole)
    elif mol.config["scf"]["type"] == 'rhf':
        mf = dft.RKS(mole)
    else:
        print(f'\n   PyOQP only support pySCF with rhf, uhf, and rohf\n')
        sys.exit(1)

    try:
        xcfun = pyscf_functional[functional]
    except KeyError:
        xcfun = 'pbe0'

    mf.xc = xcfun
    mf.max_cycle = mol.config["scf"]["maxit"]
    mf.level_shift = mol.config["scf"]["vshift"]

    try:
        os.environ["PYSCF_MAX_MEMORY"]
    except KeyError:
        mf.max_memory = 1000  # MB
    try:
        os.environ["OMP_NUM_THREADS"]
    except KeyError:
        lib.num_threads(1)

    mf.kernel()

    if mf.converged is False:
        mf = mf.newton()
        mf.kernel()

    load_pyscf_guess(mf, mol)

    return

def guess_from_pyscf_initial_density(mol, guess_type):
    """Build a PySCF atomic-density/potential initial guess for OpenQP.

    The resulting PySCF orbitals are exported through OpenQP's native JSON
    restart writer. ``sad`` uses PySCF's atomic-density guess. ``sap`` requests
    PySCF's superposition-of-atomic-potentials guess when available.
    """

    if gto is None:
        print(f'\n   PyOQP cannot find PySCF\n')
        sys.exit(1)

    if mol.usempi:
        print(f'\n   PyOQP cannot run pySCF with MPI\n')
        sys.exit(1)

    guess_type = guess_type.lower()
    if guess_type not in {'sad', 'sap'}:
        raise ValueError(f'Unsupported PySCF initial-density guess={guess_type}')

    atoms = []
    coord = mol.get_system().reshape((-1, 3)) * ANGSTROM_TO_BOHR
    for n, at in enumerate(mol.get_atoms()):
        atoms.append([ELEMENTS_NAME[SYMBOL_MAP[int(at)]], coord[n]])

    mole = gto.Mole()
    mole.atom = atoms
    mole.unit = 'Bohr'
    mole.basis = mol.config["input"]["basis"]
    mole.charge = mol.config["input"]["charge"]
    mole.spin = mol.config["scf"]["multiplicity"] - 1
    mole.output = '%s.%s.pyscf' % (mol.project_name, guess_type)
    mole.verbose = 4
    mole.build(cart=True)

    if mol.config["scf"]["type"] == 'rohf':
        mf = dft.ROKS(mole)
    elif mol.config["scf"]["type"] == 'uhf':
        mf = dft.UKS(mole)
    elif mol.config["scf"]["type"] == 'rhf':
        mf = dft.RKS(mole)
    else:
        print(f'\n   PyOQP only support pySCF with rhf, uhf, and rohf\n')
        sys.exit(1)

    functional = mol.config["input"]["functional"]
    try:
        xcfun = pyscf_functional[functional]
    except KeyError:
        xcfun = 'pbe0'

    mf.xc = xcfun
    mf.max_cycle = 0
    mf.level_shift = mol.config["scf"]["vshift"]

    try:
        os.environ["PYSCF_MAX_MEMORY"]
    except KeyError:
        mf.max_memory = 1000  # MB
    try:
        os.environ["OMP_NUM_THREADS"]
    except KeyError:
        lib.num_threads(1)

    pyscf_guess_key = 'atom' if guess_type == 'sad' else 'sap'
    try:
        dm0 = mf.get_init_guess(mole, key=pyscf_guess_key)
    except Exception as exc:
        if guess_type == 'sap':
            print(f'\n   PyOQP/PySCF SAP guess failed ({exc}); falling back to SAD\n')
            dm0 = mf.get_init_guess(mole, key='atom')
        else:
            raise

    # max_cycle=0 performs the initial Fock build/diagonalization from dm0,
    # giving the native exporter a consistent MO coefficient set to pass to OpenQP.
    mf.kernel(dm0=dm0)

    load_pyscf_guess(mf, mol)

    return


def _compact_hessian_summary(flat_hessian, max_asymmetry):
    """Return compact Hessian metadata without embedding the matrix payload."""

    return {
        "schema_version": "analytic_hessian_validation.v1",
        "report_type": "runtime_hessian_bridge",
        "matrix_payload": "omitted",
        "shape": list(flat_hessian.shape),
        "max_asymmetry_before_symmetrization": float(max_asymmetry),
    }


def _flatten_pyscf_hessian(raw_hessian, return_metadata=False):
    """Convert PySCF's (nat, nat, 3, 3) Hessian to OpenQP (3N, 3N)."""

    hess4 = np.asarray(raw_hessian, dtype=float)
    if hess4.ndim != 4 or hess4.shape[2:] != (3, 3) or hess4.shape[0] != hess4.shape[1]:
        raise ValueError(f"Expected PySCF Hessian shape (nat, nat, 3, 3), got {hess4.shape}")
    flat = hess4.transpose(0, 2, 1, 3).reshape(hess4.shape[0] * 3, hess4.shape[1] * 3)
    max_asymmetry = float(np.max(np.abs(flat - flat.T))) if flat.size else 0.0
    sym_flat = 0.5 * (flat + flat.T)
    if return_metadata:
        return sym_flat, _compact_hessian_summary(sym_flat, max_asymmetry)
    return sym_flat


def _build_pyscf_hessian_mf(mol, coord_bohr=None):
    """Build a guarded PySCF mean-field object for HF/RHF analytic Hessians."""

    if gto is None or scf is None:
        raise NotImplementedError(
            "PySCF is required for the external HF/RHF analytic Hessian bridge; "
            "no numerical fallback will be used."
        )
    # Only disable under *real* multi-rank MPI (mol.usempi is hard-coded True even
    # in serial runs; the actual indicator is the MPI manager's use_mpi flag).
    mpi_mgr = getattr(mol, "mpi_manager", None)
    if mpi_mgr is not None and getattr(mpi_mgr, "use_mpi", 0):
        raise NotImplementedError(
            "PySCF external analytic Hessian bridge is disabled under MPI; no numerical fallback will be used."
        )

    scf_type = mol.config.get("scf", {}).get("type", "rhf").lower()
    if scf_type not in ("rhf", "uhf"):
        raise NotImplementedError(
            "PySCF external analytic Hessian bridge supports RHF and UHF, got "
            f"scf.type={scf_type} (ROHF analytic Hessian is not available in PySCF; "
            "no numerical fallback will be used)."
        )

    atoms = []
    # get_system() already returns coordinates in Bohr; PySCF is built with
    # unit='Bohr' below, so no conversion is applied (note ANGSTROM_TO_BOHR is
    # actually the Bohr->Angstrom factor, 0.529177).
    if coord_bohr is None:
        coord = np.asarray(mol.get_system(), dtype=float).reshape((-1, 3))
    else:
        coord = np.asarray(coord_bohr, dtype=float).reshape((-1, 3))
    for n, at in enumerate(mol.get_atoms()):
        atoms.append([ELEMENTS_NAME[SYMBOL_MAP[int(at)]], coord[n]])

    mole = gto.Mole()
    mole.atom = atoms
    mole.unit = 'Bohr'
    mole.basis = mol.config["input"]["basis"]
    mole.charge = mol.config["input"].get("charge", 0)
    mole.spin = mol.config["scf"].get("multiplicity", 1) - 1
    mole.output = '%s.analytic_hess.pyscf' % mol.project_name
    mole.verbose = 4
    mole.build(cart=True)

    functional = mol.config.get("input", {}).get("functional", "hf").lower()
    is_hf = functional in {"", "hf"}
    if not is_hf and dft is None:
        raise NotImplementedError(
            "PySCF DFT support is required for external DFT analytic Hessians; no numerical fallback will be used."
        )
    if scf_type == "rhf":
        mf = scf.RHF(mole) if is_hf else dft.RKS(mole)
    else:  # uhf
        mf = scf.UHF(mole) if is_hf else dft.UKS(mole)
    if not is_hf:
        mf.xc = pyscf_functional.get(functional, functional)

    try:
        os.environ["PYSCF_MAX_MEMORY"]
    except KeyError:
        mf.max_memory = 1000
    try:
        os.environ["OMP_NUM_THREADS"]
    except KeyError:
        lib.num_threads(1)
    return mf


def _polarizability_driver(mf):
    """Return the PySCF polarizability driver module for an SCF object."""

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Module .* is under testing")
        if dft is not None and isinstance(mf, dft.rks.RKS):
            from pyscf.prop.polarizability import rks as polar_rks
            return polar_rks
        if dft is not None and isinstance(mf, dft.uks.UKS):
            from pyscf.prop.polarizability import uks as polar_uks
            return polar_uks
        if isinstance(mf, scf.uhf.UHF):
            from pyscf.prop.polarizability import uhf as polar_uhf
            return polar_uhf
        from pyscf.prop.polarizability import rhf as polar_rhf
        return polar_rhf


def _pyscf_polarizability(mf):
    driver = _polarizability_driver(mf)
    if driver is None:
        raise NotImplementedError("PySCF polarizability module is unavailable")
    return np.asarray(driver.Polarizability(mf).polarizability(), dtype=float)


def _pyscf_dipole(mf):
    return np.asarray(mf.dip_moment(unit='AU', verbose=0), dtype=float)


def vibrational_intensities_from_pyscf(mol, modes, displacement=1.0e-3):
    """Finite-difference IR intensities and Raman activities for HF/DFT modes.

    Dipoles and static polarizabilities are evaluated at +/- Cartesian nuclear
    displacements with the same PySCF SCF settings used by the Hessian bridge.
    """

    modes = np.asarray(modes, dtype=float)
    coord0 = np.asarray(mol.get_system(), dtype=float).reshape((-1, 3))
    ncoord = coord0.size
    if modes.ndim != 2 or modes.shape[1] != ncoord:
        raise ValueError(f"Expected modes with shape (nmode, {ncoord}), got {modes.shape}")

    from oqp.library.frequency import infrared_intensities, raman_activities

    dipole_derivs = np.zeros((3, ncoord), dtype=float)
    polar_derivs = np.zeros((3, 3, ncoord), dtype=float)
    for idx in range(ncoord):
        disp = np.zeros(ncoord, dtype=float)
        disp[idx] = displacement
        plus = (coord0.reshape(-1) + disp).reshape(coord0.shape)
        minus = (coord0.reshape(-1) - disp).reshape(coord0.shape)

        mf_plus = _build_pyscf_hessian_mf(mol, coord_bohr=plus)
        mf_plus.kernel()
        dip_plus = _pyscf_dipole(mf_plus)
        pol_plus = _pyscf_polarizability(mf_plus)

        mf_minus = _build_pyscf_hessian_mf(mol, coord_bohr=minus)
        mf_minus.kernel()
        dip_minus = _pyscf_dipole(mf_minus)
        pol_minus = _pyscf_polarizability(mf_minus)

        dipole_derivs[:, idx] = (dip_plus - dip_minus) / (2.0 * displacement)
        polar_derivs[:, :, idx] = (pol_plus - pol_minus) / (2.0 * displacement)

    ir, mode_dipoles = infrared_intensities(dipole_derivs, modes)
    raman, mode_polarizabilities = raman_activities(polar_derivs, modes)
    return {
        "infrared_intensities": ir,
        "infrared_mode_dipole_derivatives": mode_dipoles,
        "raman_activities": raman,
        "raman_mode_polarizability_derivatives": mode_polarizabilities,
        "dipole_derivatives": dipole_derivs,
        "polarizability_derivatives": polar_derivs,
        "metadata": {
            "backend": "external_pyscf_finite_difference",
            "displacement_bohr": float(displacement),
            "ir_units": "km/mol",
            "raman_units": "a.u.",
        },
    }


def analytic_hessian_from_pyscf(mol, mf_factory=None):
    """Return a conservative PySCF-backed analytic Hessian for HF/DFT only.

    This is a dispatch scaffold, not TDDFT/SF/MRSF support. Unsupported
    response-theory methods raise NotImplementedError explicitly so callers do
    not silently fall back to numerical Hessians.
    """

    method = mol.config.get("input", {}).get("method", "hf").lower()
    td_type = mol.config.get("tdhf", {}).get("type", "rpa").lower()
    if method == "tdhf":
        if td_type == "mrsf":
            raise NotImplementedError("MRSF-TDDFT analytic Hessian is not implemented")
        if td_type in {"sf", "umrsf"}:
            raise NotImplementedError(f"{td_type.upper()} analytic Hessian is not implemented")
        raise NotImplementedError(f"TDHF/TDDFT analytic Hessian is not implemented for tdhf.type={td_type}")
    if method != "hf":
        raise NotImplementedError(f"Analytic Hessian is not implemented for method={method}")

    if mf_factory is None:
        mf_factory = _build_pyscf_hessian_mf
    mf = mf_factory(mol)
    mf.kernel()
    hessian, compact_summary = _flatten_pyscf_hessian(mf.Hessian().kernel(), return_metadata=True)
    metadata = {
        "backend": "external_pyscf",
        "native_openqp_kernel": False,
        "no_numerical_fallback": True,
        "shape": list(hessian.shape),
        "max_asymmetry_before_symmetrization": compact_summary["max_asymmetry_before_symmetrization"],
        "compact_validation_summary": compact_summary,
    }
    setattr(mol, "hessian_metadata", metadata)
    return hessian, ["computed", "external_pyscf"]

