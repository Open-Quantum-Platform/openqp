import importlib
import os
import sys
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault("OQP_BACKEND_OPTIONAL", "1")
ROOT = Path(__file__).resolve().parents[1]
PYOQP = ROOT / "pyoqp"
if str(PYOQP) not in sys.path:
    sys.path.insert(0, str(PYOQP))


def _use_real_oqp_package():
    for name in list(sys.modules):
        if name == "oqp" or name.startswith("oqp."):
            sys.modules.pop(name, None)
    importlib.invalidate_caches()


def _model_integrals():
    h1e = np.array([[-1.0, -0.12], [-0.12, -0.45]], dtype=float)
    eri = np.zeros((2, 2, 2, 2), dtype=float)

    def put(p, q, r, s, value):
        for idx in {
            (p, q, r, s),
            (q, p, r, s),
            (p, q, s, r),
            (q, p, s, r),
            (r, s, p, q),
            (s, r, p, q),
            (r, s, q, p),
            (s, r, q, p),
        }:
            eri[idx] = value

    put(0, 0, 0, 0, 0.70)
    put(1, 1, 1, 1, 0.60)
    put(0, 0, 1, 1, 0.25)
    put(0, 1, 0, 1, 0.10)
    return h1e, eri


def _scalar(value):
    return float(np.real_if_close(value))


def test_spinorb_and_spatial_rdms_for_single_determinant_are_exact():
    _use_real_oqp_package()
    from oqp.library.rdm import (
        make_rdm1_spatial,
        make_rdm1_spinorb,
        make_rdm2_spatial,
        make_rdm2_spinorb,
    )

    determinants = [(1 << 0) | (1 << 2)]
    ci = np.array([1.0])

    rdm1 = make_rdm1_spinorb(ci, determinants, 4)
    rdm2 = make_rdm2_spinorb(ci, determinants, 4)

    assert np.diag(rdm1).tolist() == pytest.approx([1.0, 0.0, 1.0, 0.0])
    assert rdm2[0, 2, 0, 2] == pytest.approx(1.0)
    assert rdm2[2, 0, 0, 2] == pytest.approx(-1.0)
    assert rdm2[0, 2, 2, 0] == pytest.approx(-1.0)
    assert rdm2[2, 0, 2, 0] == pytest.approx(1.0)
    assert np.allclose(np.einsum("pqrq->pr", rdm2), rdm1)

    spatial1 = make_rdm1_spatial(ci, determinants, 2)
    spatial2 = make_rdm2_spatial(ci, determinants, 2)

    assert spatial1 == pytest.approx(np.diag([2.0, 0.0]))
    assert spatial2[0, 0, 0, 0] == pytest.approx(2.0)
    assert np.allclose(np.einsum("pqrr->pq", spatial2), spatial1)


def test_rdms_reconstruct_native_fci_energy_and_identities():
    _use_real_oqp_package()
    from oqp.library.fci import _spin_orbital_integrals, solve_fci
    from oqp.library.rdm import (
        determinant_basis,
        make_rdm1_spatial,
        make_rdm1_spinorb,
        make_rdm2_spatial,
        make_rdm2_spinorb,
    )

    h1e, eri = _model_integrals()
    nelec = (1, 1)
    energies, coeffs = solve_fci(h1e, eri, nelec, nroot=1, max_det=1000, solver="dense")
    ci = coeffs[:, 0]
    determinants = determinant_basis(h1e.shape[0], nelec)

    spin1 = make_rdm1_spinorb(ci, determinants, 2 * h1e.shape[0])
    spin2 = make_rdm2_spinorb(ci, determinants, 2 * h1e.shape[0])
    spatial1 = make_rdm1_spatial(ci, determinants, h1e.shape[0])
    spatial2 = make_rdm2_spatial(ci, determinants, h1e.shape[0])

    assert np.trace(spin1) == pytest.approx(sum(nelec), abs=1.0e-12)
    assert np.allclose(spin1, spin1.T.conj(), atol=1.0e-12)
    assert np.allclose(np.einsum("pqrq->pr", spin2), (sum(nelec) - 1) * spin1, atol=1.0e-12)
    assert np.trace(spatial1) == pytest.approx(sum(nelec), abs=1.0e-12)
    assert np.allclose(spatial1, spatial1.T.conj(), atol=1.0e-12)
    assert np.allclose(np.einsum("pqrr->pq", spatial2), (sum(nelec) - 1) * spatial1, atol=1.0e-12)

    hspin, gspin = _spin_orbital_integrals(h1e, eri)
    spin_energy = np.einsum("pq,pq", hspin, spin1) + 0.5 * np.einsum("pqrs,pqrs", gspin, spin2)
    spatial_energy = np.einsum("pq,pq", h1e, spatial1) + 0.5 * np.einsum("pqrs,pqrs", eri, spatial2)

    assert _scalar(spin_energy) == pytest.approx(energies[0], abs=1.0e-11)
    assert _scalar(spatial_energy) == pytest.approx(energies[0], abs=1.0e-11)


def test_natural_orbital_occupations_are_sorted_eigenvalues():
    _use_real_oqp_package()
    from oqp.library.rdm import natural_orbital_occupations

    np.testing.assert_allclose(
        natural_orbital_occupations(np.diag([1.75, 0.25])),
        [1.75, 0.25],
    )
    np.testing.assert_allclose(
        natural_orbital_occupations(np.asarray([[1.0, 0.2], [0.2, 1.0]], dtype=float)),
        [1.2, 0.8],
        atol=1.0e-12,
    )


@pytest.mark.parametrize(
    "rdm1",
    [
        np.asarray([[True, False], [False, True]], dtype=bool),
        [["1.0", "0.0"], ["0.0", "1.0"]],
        np.ones((2, 3), dtype=float),
        np.asarray([[1.0, 0.1], [0.2, 1.0]], dtype=float),
        np.asarray([[1.0, np.nan], [np.nan, 1.0]], dtype=float),
    ],
)
def test_natural_orbital_occupations_reject_malformed_rdm1(rdm1):
    _use_real_oqp_package()
    from oqp.library.rdm import natural_orbital_occupations

    with pytest.raises(ValueError, match="spatial 1-RDM"):
        natural_orbital_occupations(rdm1)


@pytest.mark.parametrize(
    ("atom", "active_electrons", "active_orbitals"),
    [
        ("H 0 0 0; H 0 0 0.740", 2, 2),
        (
            "O 0.000000 0.000000 0.000000;"
            " H 0.000000 -0.757000 0.587000;"
            " H 0.000000 0.757000 0.587000",
            2,
            2,
        ),
    ],
)
def test_spatial_rdms_match_pyscf_casci_reference(atom, active_electrons, active_orbitals):
    pytest.importorskip("pyscf")
    _use_real_oqp_package()
    from pyscf import gto, mcscf, scf

    from oqp.library.rdm import determinant_basis, make_rdm1_spatial, make_rdm2_spatial

    mol = gto.Mole()
    mol.atom = atom
    mol.unit = "Angstrom"
    mol.basis = "sto-3g"
    mol.verbose = 0
    mol.build(cart=True)

    mf = scf.RHF(mol)
    mf.conv_tol = 1.0e-12
    mf.kernel()
    assert mf.converged

    mc = mcscf.CASCI(mf, active_orbitals, active_electrons)
    mc.fcisolver.conv_tol = 1.0e-12
    mc.kernel()

    ci = np.asarray(mc.ci).reshape(-1)
    determinants = determinant_basis(active_orbitals, mc.nelecas)
    rdm1 = make_rdm1_spatial(ci, determinants, active_orbitals)
    rdm2 = make_rdm2_spatial(ci, determinants, active_orbitals)
    rdm1_ref, rdm2_ref = mc.fcisolver.make_rdm12(mc.ci, active_orbitals, mc.nelecas)

    assert np.allclose(rdm1, rdm1_ref, atol=1.0e-10)
    assert np.allclose(rdm2, rdm2_ref, atol=1.0e-10)
