"""Tests for oqp.interop.ccpy (external CCpy coupled-cluster library).

The unit tests need neither a compiled OpenQP nor CCpy. The end-to-end cases
run an OpenQP SCF, build the CCpy Hamiltonian from OpenQP's own integrals and
compare with reference energies from the independent PySCF -> CCpy route
(tests/data/ccpy_openqp_reference.json). They run wherever a built OpenQP is
importable; CCpy is then downloaded on demand by ensure_ccpy(), and a failed
download or build fails the test.
"""

import errno
import importlib.util
import json
import os
import subprocess
import sys
import types
from pathlib import Path

import numpy as np
import pytest

_HERE = Path(__file__).parent
DATA = _HERE / "data" / "ccpy_openqp_reference.json"
_MODULE = _HERE.parent / "pyoqp" / "oqp" / "interop" / "ccpy.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("_oqp_interop_ccpy", _MODULE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ccpy_interop = _load_module()


# --------------------------------------------------------------------------
# On-demand installation
# --------------------------------------------------------------------------

def _hide_ccpy(monkeypatch, tmp_path, failures=2):
    """Make the first ``failures`` imports of ccpy fail, later ones succeed."""
    real = ccpy_interop.importlib.import_module
    calls = {"n": 0}
    monkeypatch.setattr(ccpy_interop.sys, "prefix", str(tmp_path))
    # CCpy's floor; the version gate has its own test below
    monkeypatch.setattr(ccpy_interop.sys, "version_info", (3, 11, 0, "final", 0))

    def fake(name, *args, **kwargs):
        if name == "ccpy":
            calls["n"] += 1
            if calls["n"] <= failures:
                raise ImportError("No module named 'ccpy'")
            return types.SimpleNamespace(__name__="ccpy")
        return real(name, *args, **kwargs)

    monkeypatch.setattr(ccpy_interop.importlib, "import_module", fake)
    return calls


def test_ensure_ccpy_disabled_raises_with_pinned_command(monkeypatch, tmp_path):
    _hide_ccpy(monkeypatch, tmp_path)
    monkeypatch.setenv("OPENQP_CCPY_AUTO_INSTALL", "0")
    with pytest.raises(ImportError) as err:
        ccpy_interop.ensure_ccpy()
    assert ccpy_interop.CCPY_COMMIT in str(err.value)


def test_ensure_ccpy_installs_pinned_commit(monkeypatch, tmp_path):
    calls = _hide_ccpy(monkeypatch, tmp_path)
    monkeypatch.delenv("OPENQP_CCPY_AUTO_INSTALL", raising=False)
    seen = {}

    def fake_run(cmd, **kwargs):
        seen["cmd"] = cmd
        seen["constraints"] = open(cmd[cmd.index("--constraint") + 1]).read().split()
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(ccpy_interop.subprocess, "run", fake_run)
    mod = ccpy_interop.ensure_ccpy()
    assert mod.__name__ == "ccpy"
    assert calls["n"] == 3          # before the lock, under the lock, after pip
    assert seen["cmd"][:4] == [sys.executable, "-m", "pip", "install"]
    assert seen["cmd"][-1] == (f"coupled-cluster-py @ git+{ccpy_interop.CCPY_REPOSITORY}"
                               f"@{ccpy_interop.CCPY_COMMIT}")
    # the loaded numpy (and scipy, when installed) must stay as installed
    from importlib import metadata
    assert f"numpy=={metadata.version('numpy')}" in seen["constraints"]
    try:
        assert f"scipy=={metadata.version('scipy')}" in seen["constraints"]
    except metadata.PackageNotFoundError:
        pass


def test_build_env_adds_homebrew_openblas_last(monkeypatch, tmp_path):
    keg = tmp_path / "opt" / "openblas"
    (keg / "lib" / "pkgconfig").mkdir(parents=True)
    (keg / "lib" / "pkgconfig" / "openblas.pc").write_text("")
    monkeypatch.setattr(ccpy_interop, "_HOMEBREW_OPENBLAS", (str(tmp_path / "none"), str(keg)))
    monkeypatch.setattr(ccpy_interop.sys, "platform", "darwin")
    monkeypatch.setenv("PKG_CONFIG_PATH", "/user/first")
    paths = ccpy_interop._build_env()["PKG_CONFIG_PATH"].split(os.pathsep)
    assert paths == ["/user/first", str(keg / "lib" / "pkgconfig")]

    monkeypatch.setattr(ccpy_interop.sys, "platform", "linux")
    assert ccpy_interop._build_env()["PKG_CONFIG_PATH"] == "/user/first"


def test_ensure_ccpy_installed_by_another_process_while_waiting(monkeypatch, tmp_path):
    calls = _hide_ccpy(monkeypatch, tmp_path, failures=1)
    monkeypatch.setattr(ccpy_interop.subprocess, "run",
                        lambda *a, **k: pytest.fail("pip must not run twice"))
    assert ccpy_interop.ensure_ccpy(install=True).__name__ == "ccpy"
    assert calls["n"] == 2


def test_ensure_ccpy_pip_failure_raises(monkeypatch, tmp_path):
    _hide_ccpy(monkeypatch, tmp_path)
    monkeypatch.setattr(
        ccpy_interop.subprocess, "run",
        lambda cmd, **kw: subprocess.CompletedProcess(cmd, 1, "", "meson: no fortran"))
    with pytest.raises(RuntimeError, match="no fortran"):
        ccpy_interop.ensure_ccpy(install=True)


def test_ensure_ccpy_rejects_numpy_1_before_pip(monkeypatch, tmp_path):
    _hide_ccpy(monkeypatch, tmp_path)
    monkeypatch.setattr(ccpy_interop.np, "__version__", "1.26.4")
    monkeypatch.setattr(ccpy_interop.subprocess, "run",
                        lambda *a, **k: pytest.fail("pip must not replace a loaded numpy"))
    with pytest.raises(ImportError, match="numpy >= 2.0"):
        ccpy_interop.ensure_ccpy(install=True)


def test_ensure_ccpy_rejects_old_python_before_pip(monkeypatch, tmp_path):
    _hide_ccpy(monkeypatch, tmp_path)
    monkeypatch.setattr(ccpy_interop.sys, "version_info", (3, 9, 18, "final", 0))
    monkeypatch.setattr(ccpy_interop.subprocess, "run",
                        lambda *a, **k: pytest.fail("pip must not run on Python 3.9"))
    with pytest.raises(ImportError, match="Python >= 3.10"):
        ccpy_interop.ensure_ccpy(install=True)


def test_install_lock_falls_back_when_prefix_is_not_writable(monkeypatch, tmp_path):
    monkeypatch.setattr(ccpy_interop.sys, "prefix", str(tmp_path / "missing" / "env"))
    import tempfile
    monkeypatch.setattr(tempfile, "gettempdir", lambda: str(tmp_path))
    with ccpy_interop._install_lock():
        locks = list(tmp_path.glob("openqp-ccpy-install-*.lock"))
    assert len(locks) == 1


@pytest.mark.skipif(os.name == "nt", reason="exercises the fcntl branch")
def test_install_lock_closes_file_when_locking_fails(monkeypatch, tmp_path):
    import fcntl
    monkeypatch.setattr(ccpy_interop.sys, "prefix", str(tmp_path))
    opened = []
    real_open = open

    def tracking_open(*a, **k):
        fh = real_open(*a, **k)
        opened.append(fh)
        return fh

    monkeypatch.setattr(ccpy_interop, "open", tracking_open, raising=False)
    monkeypatch.setattr(fcntl, "flock", lambda *a: (_ for _ in ()).throw(OSError("busy")))
    with pytest.raises(OSError, match="busy"):
        with ccpy_interop._install_lock():
            pass
    assert opened and all(fh.closed for fh in opened)


@pytest.mark.parametrize("config, kwargs", [
    ({"scf": {"type": "rhf"}, "input": {"functional": "bhhlyp"}}, {}),
    ({"scf": {"type": "rks"}, "input": {}}, {}),
    ({"scf": {"type": "rhf"}, "input": {}}, {"nfrozen": -1}),
])
def test_driver_rejects_invalid_requests_before_installing(monkeypatch, config, kwargs):
    monkeypatch.setattr(ccpy_interop, "ensure_ccpy",
                        lambda *a, **k: pytest.fail("installed before validation"))
    with pytest.raises(ValueError):
        ccpy_interop.driver_from_openqp(types.SimpleNamespace(config=config), **kwargs)


def _fake_msvcrt(monkeypatch, tmp_path, failures):
    """msvcrt whose LK_LOCK raises the given OSErrors, then succeeds."""
    attempts = []
    fake = types.ModuleType("msvcrt")
    fake.LK_LOCK, fake.LK_UNLCK = 1, 0
    pending = list(failures)

    def locking(fd, mode, nbytes):
        attempts.append(mode)
        if mode == fake.LK_LOCK and pending:
            raise pending.pop(0)

    fake.locking = locking
    monkeypatch.setitem(sys.modules, "msvcrt", fake)
    monkeypatch.setattr(ccpy_interop, "_WINDOWS", True)
    monkeypatch.setattr(ccpy_interop.sys, "prefix", str(tmp_path))
    return fake, attempts


def test_windows_install_lock_waits_past_lk_lock_timeouts(monkeypatch, tmp_path):
    """msvcrt.LK_LOCK raises EDEADLOCK after ten one-second attempts (EACCES
    for a locking violation); the lock must keep waiting, because a CCpy build
    takes much longer."""
    edeadlock = getattr(errno, "EDEADLOCK", errno.EDEADLK)
    fake, attempts = _fake_msvcrt(monkeypatch, tmp_path, [
        OSError(edeadlock, "Resource deadlock avoided"),
        OSError(errno.EACCES, "Permission denied"),
        OSError(edeadlock, "Resource deadlock avoided")])
    with ccpy_interop._install_lock():
        pass
    assert attempts == [fake.LK_LOCK] * 4 + [fake.LK_UNLCK]


@pytest.mark.parametrize("err", [errno.EINVAL, errno.EBADF, None])
def test_windows_install_lock_raises_permanent_errors(monkeypatch, tmp_path, err):
    """A permanent locking error (e.g. no byte-range locks on the filesystem)
    propagates after one attempt instead of retrying forever."""
    exc = OSError(err, "locking not supported") if err else OSError("no errno")
    fake, attempts = _fake_msvcrt(monkeypatch, tmp_path, [exc] * 1000)
    with pytest.raises(OSError, match="locking not supported|no errno"):
        with ccpy_interop._install_lock():
            pytest.fail("lock acquired despite a permanent error")
    assert attempts == [fake.LK_LOCK]


class _RecordingData:
    """mol.data stand-in that records which native records were deleted."""

    def __init__(self, present=True):
        self.deleted, self.present = [], present

    def __delitem__(self, key):
        if not self.present:
            raise KeyError(key)
        self.deleted.append(key)


@pytest.mark.parametrize("fail", [False, True])
def test_eri_record_is_released_on_success_and_failure(monkeypatch, fail):
    seen = {}
    _fake_ccpy_models(monkeypatch, seen)
    if fail:
        def reject(*a):
            raise RuntimeError("reference energy mismatch")
        monkeypatch.setattr(ccpy_interop, "_check_reference", reject)
    data = _RecordingData()
    mol = types.SimpleNamespace(config={"scf": {"type": "rohf"}, "input": {}}, data=data)
    if fail:
        with pytest.raises(RuntimeError, match="mismatch"):
            ccpy_interop.hamiltonian_from_openqp(mol)
    else:
        ccpy_interop.hamiltonian_from_openqp(mol)
    assert data.deleted == ["OQP::ERI_AO"]


def test_missing_eri_record_is_not_an_error(monkeypatch):
    _fake_ccpy_models(monkeypatch, {})
    mol = types.SimpleNamespace(config={"scf": {"type": "uhf"}, "input": {}},
                                data=_RecordingData(present=False))
    ccpy_interop.hamiltonian_from_openqp(mol)


class _CountData:
    """mol.data stand-in with only the counts: no integrals, no orbitals."""

    def __init__(self, nbf, na, nb):
        self._basis, self._items = {"nbf": nbf}, {"nelec_A": [na], "nelec_B": [nb]}

    def get_basis(self):
        return self._basis

    def __getitem__(self, key):
        return self._items[key]


@pytest.mark.parametrize("kwargs, match", [({"nfrozen": 2}, "nfrozen=2 exceeds"),
                                           ({"ndelete": 2}, "ndelete=2 exceeds")])
def test_orbital_count_limits_are_checked_before_installing(monkeypatch, kwargs, match):
    monkeypatch.setattr(ccpy_interop, "ensure_ccpy",
                        lambda *a, **k: pytest.fail("installed before validation"))
    monkeypatch.setattr(ccpy_interop, "_openqp_ao_data",
                        lambda mol: pytest.fail("built integrals before validation"))
    mol = types.SimpleNamespace(config={"scf": {"type": "rohf"}, "input": {}},
                                data=_CountData(nbf=3, na=2, nb=1))
    with pytest.raises(ValueError, match=match):
        ccpy_interop.hamiltonian_from_openqp(mol, **kwargs)


def test_negative_frozen_or_deleted_orbitals_are_rejected(monkeypatch):
    monkeypatch.setattr(ccpy_interop, "ensure_ccpy",
                        lambda *a, **k: pytest.fail("checked after ensure_ccpy"))
    mol = types.SimpleNamespace(config={"scf": {"type": "rhf"}, "input": {}})
    with pytest.raises(ValueError, match="must be >= 0"):
        ccpy_interop.hamiltonian_from_openqp(mol, nfrozen=-1)


def _fake_ccpy_models(monkeypatch, seen):
    """Minimal stand-ins for the CCpy modules hamiltonian_from_openqp imports."""
    def module(name, **attrs):
        mod = types.ModuleType(name)
        mod.__dict__.update(attrs)
        monkeypatch.setitem(sys.modules, name, mod)

    def system(kind):
        def make(*args, **kwargs):
            seen[kind] = kwargs
            return types.SimpleNamespace()
        return make

    for name in ("ccpy", "ccpy.energy", "ccpy.interfaces", "ccpy.models"):
        module(name)
    module("ccpy.energy.hf_energy", calc_hf_energy=lambda *a: 0.0,
           calc_hf_frozen_core_energy=lambda *a: 0.0)
    module("ccpy.interfaces.fcidump_tools", build_rohf_fock=lambda *a: (None, None))
    module("ccpy.models.integrals", getHamiltonian=lambda *a: "H",
           getUHFHamiltonian=lambda *a: "H")
    module("ccpy.models.system", System=system("System"), UHFSystem=system("UHFSystem"))
    monkeypatch.setattr(ccpy_interop, "ensure_ccpy", lambda *a, **k: None)
    monkeypatch.setattr(ccpy_interop, "_check_reference", lambda *a: None)
    monkeypatch.setattr(ccpy_interop, "_scf_energy", lambda mol: 0.0)
    monkeypatch.setattr(ccpy_interop, "_counts", lambda mol: (3, 1, 0))
    n = 3
    monkeypatch.setattr(ccpy_interop, "_openqp_ao_data", lambda mol: {
        "nbf": n, "hcore": np.eye(n), "v_ao": np.zeros((n, n, n, n)),
        "ca": np.eye(n), "cb": np.eye(n), "ea": np.zeros(n), "eb": np.zeros(n),
        "na": 1, "nb": 0, "vnn": 0.0})


@pytest.mark.parametrize("scf_type, kind", [("rohf", "System"), ("uhf", "UHFSystem")])
def test_charge_reaches_the_ccpy_system(monkeypatch, scf_type, kind):
    seen = {}
    _fake_ccpy_models(monkeypatch, seen)
    mol = types.SimpleNamespace(config={"scf": {"type": scf_type},
                                        "input": {"charge": 1}},
                                data=_RecordingData())
    ccpy_interop.hamiltonian_from_openqp(mol)
    assert seen[kind]["charge"] == 1


class _FakeDriver:
    """CCpy's rule: RHF_symmetry = (noccupied_alpha == noccupied_beta)."""

    def __init__(self, system, hamiltonian):
        self.system = system
        self.options = {"RHF_symmetry": system.noccupied_alpha == system.noccupied_beta}


def _fake_driver_module(monkeypatch):
    for name in ("ccpy", "ccpy.drivers"):
        monkeypatch.setitem(sys.modules, name, types.ModuleType(name))
    mod = types.ModuleType("ccpy.drivers.driver")
    mod.Driver = _FakeDriver
    monkeypatch.setitem(sys.modules, "ccpy.drivers.driver", mod)
    monkeypatch.setattr(ccpy_interop, "ensure_ccpy", lambda *a, **k: None)
    system = types.SimpleNamespace(noccupied_alpha=4, noccupied_beta=4)
    monkeypatch.setattr(ccpy_interop, "hamiltonian_from_openqp",
                        lambda mol, **k: (system, object()))


@pytest.mark.parametrize("scf_type, expected", [("uhf", False), ("rhf", True)])
def test_uhf_singlet_driver_is_unrestricted(monkeypatch, scf_type, expected):
    _fake_driver_module(monkeypatch)
    mol = types.SimpleNamespace(config={"scf": {"type": scf_type}, "input": {}})
    assert ccpy_interop.driver_from_openqp(mol).options["RHF_symmetry"] is expected


# --------------------------------------------------------------------------
# UHF reference energy from AO integrals
# --------------------------------------------------------------------------

def _random_eri_chemist(n, rng):
    """(pq|rs) with the 8-fold permutational symmetry of real orbitals."""
    x = rng.standard_normal((n, n, n, n))
    x = x + x.transpose(1, 0, 2, 3)
    x = x + x.transpose(0, 1, 3, 2)
    return x + x.transpose(2, 3, 0, 1)


def test_uhf_energy_matches_mo_basis_formula():
    rng = np.random.default_rng(7)
    n, na, nb = 6, 3, 2
    h = rng.standard_normal((n, n))
    h = h + h.T
    g = _random_eri_chemist(n, rng)
    ca = rng.standard_normal((n, n))          # alpha and beta differ, not orthogonal:
    cb = rng.standard_normal((n, n))          # a swapped index or spin changes E
    v_ao = np.transpose(g, (0, 2, 1, 3))      # physicist, as the adapter passes it

    got = ccpy_interop._uhf_energy(h, v_ao, ca, cb, na, nb)

    def mo(c1, c2):
        return np.einsum("pqrs,pi,qj,rk,sl->ijkl", g, c1, c1, c2, c2)
    haa, hbb = ca.T @ h @ ca, cb.T @ h @ cb
    gaa, gbb, gab = mo(ca, ca), mo(cb, cb), mo(ca, cb)
    oa, ob = range(na), range(nb)
    ref = sum(haa[i, i] for i in oa) + sum(hbb[i, i] for i in ob)
    ref += 0.5 * sum(gaa[i, i, j, j] - gaa[i, j, j, i] for i in oa for j in oa)
    ref += 0.5 * sum(gbb[i, i, j, j] - gbb[i, j, j, i] for i in ob for j in ob)
    ref += sum(gab[i, i, j, j] for i in oa for j in ob)
    assert got == pytest.approx(ref, abs=1e-10, rel=1e-12)


# --------------------------------------------------------------------------
# OpenQP -> CCpy end to end
# --------------------------------------------------------------------------

def _have_native_oqp():
    try:
        import oqp
    except Exception:
        return False
    return hasattr(oqp, "int2e")


native = pytest.mark.skipif(not _have_native_oqp(),
                            reason="native OpenQP library with int2e not importable")

_RUNNER = r"""
import json, sys
from oqp.pyoqp import Runner
from oqp.interop.ccpy import driver_from_openqp

def run(case, k):
    atoms = "\n".join("   %d %.10f %.10f %.10f" % tuple(a) for a in case["atoms"])
    ispher = "ispher=true\n" if case.get("ispher") else ""
    inp = "case%d.inp" % k
    open(inp, "w").write(
        "[input]\nsystem=\n%s\ncharge=0\nruntype=energy\nbasis=%s\nmethod=hf\n%s\n"
        "[guess]\ntype=huckel\n\n[scf]\ntype=%s\nmultiplicity=%d\nconv=1.0e-11\nmaxit=200\n"
        % (atoms, case["basis"], ispher, case["scf_type"], case["multiplicity"]))
    runner = Runner(project="case%d" % k, input_file=inp, log="case%d.log" % k)
    runner.run()
    drv = driver_from_openqp(runner.mol, nfrozen=case["nfrozen"])
    drv.options["energy_convergence"] = 1.0e-10
    drv.options["amp_convergence"] = 1.0e-8
    drv.options["maximum_iterations"] = 200
    return case_tasks(case, drv)

def case_tasks(case, drv):
    out = {}
    drv.run_cc(method="ccsd")
    out["e_ref"] = drv.system.reference_energy
    out["e_ccsd"] = drv.system.reference_energy + drv.correlation_energy
    if "crcc23" in case["tasks"]:
        drv.run_hbar(method="ccsd")
        drv.run_leftcc(method="left_ccsd")
        drv.run_ccp3(method="crcc23")
        out["e_crcc23_d"] = (drv.system.reference_energy + drv.correlation_energy
                             + drv.deltap3[0]["D"])
    if "sfeomccsd" in case["tasks"]:
        drv.run_hbar(method="ccsd")
        n = case["nroot"]
        drv.run_guess(method="sfcis", multiplicity=1, roots_per_irrep={"A": n},
                      nact_occupied=0, nact_unoccupied=0, use_symmetry=False)
        drv.run_sfeomcc(method="sfeomccsd", state_index=list(range(n)))
        out["omega_singlets"] = [float(drv.vertical_excitation_energy[i]) for i in range(n)]
    return out

results = [run(case, k) for k, case in enumerate(json.loads(sys.argv[1]))]
print("CCPY_RESULT " + json.dumps(results))
"""


def _cases():
    return json.loads(DATA.read_text())["cases"]


def _run_cases(cases, tmp_path, threads):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    proc = subprocess.run([sys.executable, "-c", _RUNNER, json.dumps(cases)],
                          cwd=tmp_path, capture_output=True, text=True,
                          timeout=3600, env=env)
    lines = [l for l in proc.stdout.splitlines() if l.startswith("CCPY_RESULT ")]
    assert proc.returncode == 0 and lines, proc.stdout[-4000:] + proc.stderr[-4000:]
    return json.loads(lines[-1].split(" ", 1)[1])


def _check(got, case):
    tol = json.loads(DATA.read_text())["tolerance"]
    ref = case["reference"]
    assert got["e_ref"] == pytest.approx(ref["e_ref"], abs=tol)
    assert got["e_ccsd"] == pytest.approx(ref["e_ccsd"], abs=tol)
    if "e_crcc23_d" in ref:
        assert got["e_crcc23_d"] == pytest.approx(ref["e_crcc23_d"], abs=tol)
    if "omega_singlets" in ref:
        assert got["omega_singlets"] == pytest.approx(ref["omega_singlets"], abs=tol)


@native
@pytest.mark.parametrize("case", _cases(), ids=lambda c: c["name"])
def test_openqp_orbitals_reproduce_pyscf_ccpy(case, tmp_path):
    [got] = _run_cases([case], tmp_path, case.get("threads", 2))
    _check(got, case)


@native
def test_repeated_calls_with_changing_basis_size(tmp_path):
    """Several molecules in one process, nbf up, down and up again."""
    by_name = {c["name"]: c for c in _cases()}
    seq = [by_name[n] for n in ("h2o_rhf_ccpvdz_spherical", "ch2_uhf_631gs",
                                "h2o_rhf_631gs_fc1", "h2o_rhf_ccpvdz_spherical")]
    for got, case in zip(_run_cases(seq, tmp_path, 2), seq):
        _check(got, case)


def test_dft_orbitals_are_rejected_before_ccpy_is_needed(monkeypatch):
    monkeypatch.setattr(ccpy_interop, "ensure_ccpy",
                        lambda *a, **k: pytest.fail("checked after ensure_ccpy"))
    mol = types.SimpleNamespace(config={"scf": {"type": "rhf"},
                                        "input": {"functional": "bhhlyp"}})
    with pytest.raises(ValueError, match="Hartree-Fock"):
        ccpy_interop.hamiltonian_from_openqp(mol)


def test_correlated_run_is_rejected_before_ccpy_is_needed(monkeypatch):
    monkeypatch.setattr(ccpy_interop, "ensure_ccpy",
                        lambda *a, **k: pytest.fail("checked after ensure_ccpy"))
    mol = types.SimpleNamespace(config={"scf": {"type": "rhf"},
                                        "input": {"method": "ccsd"}})
    with pytest.raises(ValueError, match="method=hf"):
        ccpy_interop.hamiltonian_from_openqp(mol)


def test_unsupported_reference_is_rejected(monkeypatch):
    monkeypatch.setattr(ccpy_interop, "ensure_ccpy",
                        lambda *a, **k: pytest.fail("checked after ensure_ccpy"))
    mol = types.SimpleNamespace(config={"scf": {"type": "rks"}, "input": {}})
    with pytest.raises(ValueError, match="RHF, ROHF or UHF"):
        ccpy_interop.hamiltonian_from_openqp(mol)
