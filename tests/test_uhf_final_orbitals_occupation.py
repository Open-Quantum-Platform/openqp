"""Published SCF orbitals must span the converged MOM/restricted-MO density.

Two converged exits re-diagonalise the Fock with an aufbau fill: the UHF
branch after the SCF loop, and the SOSCF/TRAH handler for RHF and UHF.
A non-aufbau solution selected by restricted-MO reordering (rstctmo, the
delta-CHP/XAS core-hole route) or by MOM was therefore published as the
aufbau state while DM_A/DM_B kept the converged density.

H2O/6-31G UHF, swapmo=1,6 (alpha O 1s -> LUMO), rstctmo, before the fix:
the stored occupied alpha orbitals held 4.000 of the 5 alpha electrons of
DM_A, max |DM_A - C_occ C_occ^T| = 1.45, and the occupied particle orbital
(eps = -0.091 Eh) was printed and stored as a virtual.  A JSON restart from
that run then started from the ground-state occupation and converged to the
ground state (-75.985 Eh instead of -56.162 Eh).  A SOSCF restart from
correct orbitals kept the state but published the aufbau orbitals again.
"""

import json
import os
import subprocess
import sys

import pytest

np = pytest.importorskip("numpy")

# Distorted water, no symmetry; 13 basis functions in 6-31G.
H2O = """\
[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.633194329  -0.483194329  -0.554469223
charge=0
runtype=energy
basis=6-31g
method=hf

[guess]
{guess}
save_mol=True

[scf]
type={scf}
multiplicity={mult}
maxit=100
conv=1.0e-8
{extra}
"""

# scf type -> (multiplicity, swapmo, n_alpha, n_beta)
STATES = {
    "uhf": (1, "1,6", 5, 5),    # alpha O 1s -> LUMO core excitation
    "rhf": (1, "5,6", 5, 5),    # HOMO^2 -> LUMO^2 double excitation
    "rohf": (3, "1,5", 6, 4),   # O 1s singly occupied core-hole triplet
}


def _have_native_oqp():
    try:
        import oqp  # noqa: F401
    except Exception:
        return False
    return True


pytestmark = pytest.mark.skipif(
    not _have_native_oqp(),
    reason="native OpenQP library not importable",
)


def _run(tmp_path, name, deck, threads=1):
    inp = tmp_path / f"{name}.inp"
    inp.write_text(deck)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    proc = subprocess.run(
        [sys.executable, "-m", "oqp.pyoqp", str(inp)],
        capture_output=True, text=True, cwd=str(tmp_path), env=env, timeout=900,
    )
    log_path = inp.with_suffix(".log")
    log = log_path.read_text() if log_path.exists() else ""
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    assert "SCF convergence achieved" in log, log[-3000:]
    data = json.loads((tmp_path / f"{name}.json").read_text())
    return data, tmp_path / f"{name}.json"


def _unpack(packed, n):
    m = np.zeros((n, n))
    m[np.tril_indices(n)] = packed
    return m + np.tril(m, -1).T


def _check_spin(data, spin, nocc, factor=1.0):
    """Stored occupied orbitals reproduce the stored density of one spin."""
    vec = np.asarray(data[f"OQP::VEC_MO_{spin}"], dtype=float)   # rows are MOs
    n = vec.shape[1]
    dens = _unpack(np.asarray(data[f"OQP::DM_{spin}"], dtype=float).ravel(), n) / factor
    smat = _unpack(np.asarray(data["OQP::SM"], dtype=float).ravel(), n)
    occ = vec[:nocc]
    # Electrons of the density captured by the stored occupied orbitals:
    # nocc when they span the converged occupied space, nocc-1 per lost orbital.
    captured = np.trace(occ @ smat @ dens @ smat @ occ.T)
    assert captured == pytest.approx(nocc, abs=1e-6), (spin, captured)
    assert np.abs(dens - occ.T @ occ).max() < 1e-5, spin


def _check_state(data, scf):
    _, _, na, nb = STATES[scf]
    if scf == "rhf":
        _check_spin(data, "A", na, factor=2.0)
    else:
        _check_spin(data, "A", na)
        _check_spin(data, "B", nb)


@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("mult,swap,nalpha,nbeta", [(1, "1,6", 5, 5), (3, "1,7", 6, 4)])
def test_uhf_rstctmo_core_hole_orbitals_span_density(tmp_path, threads, mult, swap, nalpha, nbeta):
    data, _ = _run(tmp_path, f"ch_m{mult}_t{threads}",
                   H2O.format(guess=f"type=huckel\nswapmo={swap}", scf="uhf",
                              mult=mult, extra="rstctmo=True"), threads)
    _check_spin(data, "A", nalpha)
    _check_spin(data, "B", nbeta)
    # The alpha 1s stays a hole: the lowest alpha orbital energy is unoccupied.
    e_a = np.asarray(data["OQP::E_MO_A"], dtype=float).ravel()
    assert e_a.min() < -15.0
    assert e_a.min() not in e_a[:nalpha]


# TRAH is left out: the input checker refuses it together with rstctmo/MOM.
@pytest.mark.parametrize("scf,converger", [
    ("uhf", "diis"), ("uhf", "soscf"),
    ("rhf", "diis"), ("rhf", "soscf"),
    ("rohf", "diis"), ("rohf", "soscf"),
])
def test_json_restart_keeps_restricted_mo_state(tmp_path, scf, converger):
    mult, swap, _, _ = STATES[scf]
    first, js = _run(tmp_path, "first",
                     H2O.format(guess=f"type=huckel\nswapmo={swap}", scf=scf,
                                mult=mult, extra="rstctmo=True"))
    _check_state(first, scf)
    again, _ = _run(tmp_path, "again",
                    H2O.format(guess=f"type=json\nfile={js}", scf=scf, mult=mult,
                               extra=f"rstctmo=True\nconverger_type={converger}"))
    _check_state(again, scf)
    assert again["energy"] == pytest.approx(first["energy"], abs=1e-6)


def test_uhf_mom_soscf_restart_keeps_core_hole(tmp_path):
    first, js = _run(tmp_path, "first",
                     H2O.format(guess="type=huckel\nswapmo=1,6", scf="uhf", mult=1,
                                extra="rstctmo=True"))
    again, _ = _run(tmp_path, "again",
                    H2O.format(guess=f"type=json\nfile={js}", scf="uhf", mult=1,
                               extra="mom=True\nmom_switch=1.0\nconverger_type=soscf"))
    _check_state(again, "uhf")
    assert again["energy"] == pytest.approx(first["energy"], abs=1e-6)


# A JSON guess whose alpha O 1s and LUMO rows are exchanged in the file.  For
# UHF the JSON guess rebuilds the densities from the stored orbitals, and
# swapmo now rebuilds a supplied density as well, so orbitals and
# density of a supplied guess always agree: both name the core hole.  With and
# without MOM the SCF must keep that state and publish orbitals that build the
# stored density -- the same state as the swapmo route.
@pytest.mark.parametrize("mom", [False, True])
def test_uhf_reordered_json_orbitals_keep_their_state(tmp_path, mom):
    ground, js = _run(tmp_path, "ground",
                      H2O.format(guess="type=huckel", scf="uhf", mult=1, extra=""))
    # Reference: the same swap made with swapmo, which rebuilds the density.
    hole, _ = _run(tmp_path, "hole",
                   H2O.format(guess=f"type=json\nfile={js}\nswapmo=1,6", scf="uhf", mult=1,
                              extra="mom=True\nmom_switch=1.0\nconverger_type=soscf"))
    data = json.loads(js.read_text())
    vec = np.array(data["OQP::VEC_MO_A"], dtype=float)          # rows are MOs
    vec[[0, 5]] = vec[[5, 0]]
    data["OQP::VEC_MO_A"] = vec.tolist()
    e_a = np.array(data["OQP::E_MO_A"], dtype=float)
    flat = e_a.reshape(-1)
    flat[[0, 5]] = flat[[5, 0]]
    data["OQP::E_MO_A"] = flat.reshape(e_a.shape).tolist()
    mixed = tmp_path / "mixed.json"
    mixed.write_text(json.dumps(data))
    extra = "mom=True\nmom_switch=1.0\n" if mom else ""
    again, _ = _run(tmp_path, "again",
                    H2O.format(guess=f"type=json\nfile={mixed}", scf="uhf", mult=1,
                               extra=extra + "converger_type=soscf"))
    _check_state(again, "uhf")
    assert again["energy"] == pytest.approx(hole["energy"], abs=1e-6)
    assert again["energy"] > ground["energy"] + 10.0
