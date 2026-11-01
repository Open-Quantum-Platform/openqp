"""Coupled-cluster calculations with the external CCpy library on OpenQP orbitals.

CCpy (https://github.com/piecuch-group/ccpy, GPL-3.0) provides CCSD, CR-CC(2,3),
CC(P;Q), EOM-CC, IP/EA-EOM-CC and spin-flip EOM-CC. It is an external library:
OpenQP does not vendor or redistribute it. :func:`ensure_ccpy` downloads and
builds the pinned commit into the running Python environment the first time it
is needed, unless ``OPENQP_CCPY_AUTO_INSTALL=0``.

The Hamiltonian is built directly from OpenQP's own integrals, so the AO
ordering and normalization are those of the SCF that produced the orbitals:

* ``OQP::Hcore``                    packed lower triangle -> AO core Hamiltonian
* :func:`oqp.library.ints_2e.eri_ao`  AO ERIs ``(mu nu|la si)``, chemist order
* ``OQP::VEC_MO_A`` / ``OQP::VEC_MO_B``  MO coefficients
* ``OQP::E_MO_A`` / ``OQP::E_MO_B``      orbital energies
* ``nelec_A``, ``nelec_B``, ``vnn``       electron counts, nuclear repulsion

CCpy works in physicist order ``<pq|rs> = (pr|qs)``, hence
``transpose(0, 2, 1, 3)``. Point-group symmetry is not passed on (C1). The
two-electron integrals are held in core, as in ``oqp.quantum``. Building the
Hamiltonian peaks at about 6.6 times ``8 * nbf**4`` bytes for RHF and UHF
alike: 0.7 GiB at nbf = 58, 3.9 GiB at nbf = 92 and 9.3 GiB at nbf = 115
(H2O, spherical cc-pVTZ / aug-cc-pVTZ / cc-pVQZ, peak RSS on Apple M3 Ultra,
2026-10-19), before CCpy allocates its amplitudes. The build itself took
0.3 / 1.5 / 3.7 s on one thread there. These come from a single-host timing
record, not from the OpenQP verification gate.

Usage::

    from oqp.interop.ccpy import driver_from_openqp
    driver = driver_from_openqp(mol, nfrozen=1)    # mol: converged OpenQP SCF
    driver.run_cc(method="ccsd")
    driver.run_hbar(method="ccsd")
    driver.run_leftcc(method="left_ccsd")
    driver.run_ccp3(method="crcc23")
"""
import errno
import importlib
import os
import subprocess
import sys
import traceback

import numpy as np

__all__ = [
    "CCPY_REPOSITORY", "CCPY_COMMIT", "ensure_ccpy",
    "hamiltonian_from_openqp", "driver_from_openqp",
]

CCPY_REPOSITORY = "https://github.com/piecuch-group/ccpy"
#: Commit validated against OpenQP RHF/UHF/ROHF references (2026-10-17).
CCPY_COMMIT = "cf72da4865e054f7a935e2df27e920cf5e335435"

#: The pinned CCpy commit declares requires-python >= 3.10.
CCPY_MIN_PYTHON = (3, 10)

#: ... and numpy >= 2.0. OpenQP itself accepts numpy >= 1.20.
CCPY_MIN_NUMPY = (2, 0)

#: Runtime dependencies this process has already imported (numpy) or that
#: OpenQP needs at its own version (scipy). pip must not replace them under
#: the running interpreter, so the install pins them to what is installed.
_PINNED_DISTRIBUTIONS = ("numpy", "scipy")

#: Same agreement threshold CCpy's own PySCF interface applies.
REFERENCE_ENERGY_TOLERANCE = 1.0e-6


def _auto_install_enabled():
    value = os.environ.get("OPENQP_CCPY_AUTO_INSTALL", "1").strip().lower()
    return value not in ("0", "false", "no", "off")


def ensure_ccpy(install=None):
    """Import CCpy, downloading and building the pinned commit when it is absent.

    ``install=None`` follows ``OPENQP_CCPY_AUTO_INSTALL`` (default on). The
    build needs git, a Fortran compiler and OpenBLAS (CCpy's meson defaults).
    Returns the imported ``ccpy`` module; raises ``ImportError`` when CCpy is
    absent and installation is disabled, ``RuntimeError`` when pip fails.
    """
    try:
        return importlib.import_module("ccpy")
    except ImportError as exc:
        missing = exc
    if install is None:
        install = _auto_install_enabled()
    spec = f"coupled-cluster-py @ git+{CCPY_REPOSITORY}@{CCPY_COMMIT}"
    if tuple(sys.version_info[:2]) < CCPY_MIN_PYTHON:
        raise ImportError(
            "CCpy needs Python >= %d.%d; this OpenQP runs on Python %d.%d."
            % (CCPY_MIN_PYTHON + tuple(sys.version_info[:2]))) from missing
    numpy_version = tuple(int(x) for x in np.__version__.split(".")[:2])
    if numpy_version < CCPY_MIN_NUMPY:
        raise ImportError(
            f"CCpy needs numpy >= {CCPY_MIN_NUMPY[0]}.{CCPY_MIN_NUMPY[1]}; this "
            f"environment has numpy {np.__version__}, already loaded by OpenQP. "
            "Upgrade numpy (and scipy) in a separate step, restart, and retry.") from missing
    if not install:
        raise ImportError(
            "CCpy is not installed and OPENQP_CCPY_AUTO_INSTALL=0. Install it with: "
            f"{sys.executable} -m pip install '{spec}'") from missing

    # Concurrent calculations in one environment must not run two pip builds
    # into the same site-packages; the second waits and then finds CCpy.
    with _install_lock():
        importlib.invalidate_caches()
        try:
            return importlib.import_module("ccpy")
        except ImportError:
            pass
        print(f"   CCpy not found; downloading and building {CCPY_REPOSITORY}"
              f"@{CCPY_COMMIT} (GPL-3.0, external library) into {sys.prefix}. "
              "Set OPENQP_CCPY_AUTO_INSTALL=0 to disable.", flush=True)
        with _pin_constraints() as constraints:
            proc = subprocess.run(
                [sys.executable, "-m", "pip", "install", "--constraint", constraints, spec],
                capture_output=True, text=True, env=_build_env())
        if proc.returncode != 0:
            tail = "\n".join((proc.stdout + proc.stderr).splitlines()[-40:])
            raise RuntimeError(
                f"pip could not install CCpy (exit {proc.returncode}). The build "
                "needs git, a Fortran compiler and OpenBLAS, and must not change "
                f"the installed {', '.join(_PINNED_DISTRIBUTIONS)}.\n{tail}")
        importlib.invalidate_caches()
        return importlib.import_module("ccpy")


class _pin_constraints:
    """Temporary pip constraints file pinning the installed runtime packages."""

    def __enter__(self):
        import tempfile
        from importlib import metadata
        lines = []
        for name in _PINNED_DISTRIBUTIONS:
            try:
                lines.append(f"{name}=={metadata.version(name)}")
            except metadata.PackageNotFoundError:
                continue
        fd, self.path = tempfile.mkstemp(prefix="openqp-ccpy-constraints-", suffix=".txt")
        with os.fdopen(fd, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        return self.path

    def __exit__(self, *exc):
        os.unlink(self.path)
        return False


#: Homebrew installs OpenBLAS keg-only, so its openblas.pc is not on the
#: default pkg-config path that CCpy's meson build searches.
_HOMEBREW_OPENBLAS = ("/opt/homebrew/opt/openblas", "/usr/local/opt/openblas")


def _build_env():
    """Environment for the CCpy pip build.

    On macOS, append Homebrew's OpenBLAS pkg-config directory (lowest
    priority, so an OpenBLAS the user already exposes still wins).
    """
    env = dict(os.environ)
    if sys.platform != "darwin":
        return env
    for prefix in _HOMEBREW_OPENBLAS:
        pcdir = os.path.join(prefix, "lib", "pkgconfig")
        if os.path.isfile(os.path.join(pcdir, "openblas.pc")):
            paths = [p for p in env.get("PKG_CONFIG_PATH", "").split(os.pathsep) if p]
            if pcdir not in paths:
                env["PKG_CONFIG_PATH"] = os.pathsep.join(paths + [pcdir])
            break
    return env


def _lock_paths():
    """Lock file in the environment, else in the temp directory.

    The temp-directory name is derived from sys.prefix, so processes that
    install into the same environment still share one lock.
    """
    import hashlib
    import tempfile
    tag = hashlib.sha256(os.path.realpath(sys.prefix).encode()).hexdigest()[:16]
    return (os.path.join(sys.prefix, ".openqp-ccpy-install.lock"),
            os.path.join(tempfile.gettempdir(), f"openqp-ccpy-install-{tag}.lock"))


_WINDOWS = os.name == "nt"

#: errno values msvcrt.locking reports when another process holds the range:
#: EDEADLOCK after LK_LOCK's ten one-second attempts, EACCES for a locking
#: violation (Microsoft _locking documentation; CPython raises OSError from errno).
_LOCK_CONTENTION_ERRNOS = frozenset({getattr(errno, "EDEADLOCK", errno.EDEADLK), errno.EACCES})


class _install_lock:
    """Exclusive lock serializing the CCpy install for one environment."""

    def __enter__(self):
        errors = []
        for path in _lock_paths():
            try:
                self._fh = open(path, "a+")
                break
            except OSError as exc:
                errors.append(f"{path}: {exc}")
        else:
            raise OSError("cannot create a CCpy install lock: " + "; ".join(errors))
        try:
            if _WINDOWS:
                import msvcrt
                self._fh.seek(0)
                self._unlock = lambda: msvcrt.locking(self._fh.fileno(), msvcrt.LK_UNLCK, 1)
                # LK_LOCK gives up after ten one-second attempts, far shorter
                # than a CCpy build; keep waiting while another process holds
                # the lock. Any other error (EBADF, EINVAL, a filesystem
                # without byte-range locks) is permanent and propagates.
                while True:
                    try:
                        msvcrt.locking(self._fh.fileno(), msvcrt.LK_LOCK, 1)
                        break
                    except OSError as exc:
                        if exc.errno not in _LOCK_CONTENTION_ERRNOS:
                            raise
            else:
                import fcntl
                self._unlock = lambda: fcntl.flock(self._fh, fcntl.LOCK_UN)
                fcntl.flock(self._fh, fcntl.LOCK_EX)
        except BaseException:
            self._fh.close()
            raise
        return self

    def __exit__(self, *exc):
        try:
            self._unlock()
        finally:
            self._fh.close()
        return False


def _scalar(mol, key):
    return float(np.asarray(mol.data[key]).ravel()[0])


def _count(mol, key):
    return int(np.asarray(mol.data[key]).ravel()[0])


def _counts(mol):
    """(nbf, n_alpha, n_beta) of a finished OpenQP run, without any integrals."""
    basis = mol.data.get_basis()
    if not basis:
        raise RuntimeError("No basis available; run the OpenQP SCF first.")
    return int(basis["nbf"]), _count(mol, "nelec_A"), _count(mol, "nelec_B")


def _openqp_ao_data(mol):
    from oqp.library.ints_2e import eri_ao
    from oqp.quantum.integrals import unpack_triangular

    nbf, _, _ = _counts(mol)
    return {
        "nbf": nbf,
        "hcore": unpack_triangular(mol.data["OQP::Hcore"], nbf),
        "v_ao": np.transpose(eri_ao(mol), (0, 2, 1, 3)),     # chemist -> physicist
        "ca": mol._mo_coefficients("OQP::VEC_MO_A", nbf),
        "cb": mol._mo_coefficients("OQP::VEC_MO_B", nbf),
        "ea": np.asarray(mol.data["OQP::E_MO_A"], dtype=float).ravel()[:nbf],
        "eb": np.asarray(mol.data["OQP::E_MO_B"], dtype=float).ravel()[:nbf],
        "na": _count(mol, "nelec_A"),
        "nb": _count(mol, "nelec_B"),
        "vnn": _scalar(mol, "vnn"),
    }


def _uhf_energy(hcore, v_ao, ca, cb, na, nb):
    """Electronic UHF energy of the determinant (ca[:, :na], cb[:, :nb])."""
    da = ca[:, :na] @ ca[:, :na].T
    db = cb[:, :nb] @ cb[:, :nb].T
    dt = da + db
    # v_ao is physicist <mu la|nu si> = (mu nu|la si)
    j = np.einsum("mlns,ls->mn", v_ao, dt, optimize=True)
    ka = np.einsum("mlsn,ls->mn", v_ao, da, optimize=True)
    kb = np.einsum("mlsn,ls->mn", v_ao, db, optimize=True)
    return (np.einsum("mn,mn->", dt, hcore)
            + 0.5 * np.einsum("mn,mn->", dt, j)
            - 0.5 * np.einsum("mn,mn->", da, ka)
            - 0.5 * np.einsum("mn,mn->", db, kb))


def _check_reference(name, e_integrals, e_scf):
    if abs(e_integrals - e_scf) > REFERENCE_ENERGY_TOLERANCE:
        raise RuntimeError(
            f"{name} energy rebuilt from OpenQP integrals ({e_integrals:.10f}) "
            f"differs from the stored OpenQP energy ({e_scf:.10f}). Use the "
            "molecule of a converged method=hf run: a correlated method "
            "overwrites the stored energy, and other orbitals do not match it.")


def _scf_energy(mol):
    if not bool(mol.data._data.mol_energy.SCF_converged):
        raise RuntimeError("The OpenQP SCF is not converged; CCpy needs a "
                           "converged RHF, ROHF or UHF reference.")
    return float(mol.get_scf_energy())


def hamiltonian_from_openqp(mol, nfrozen=0, ndelete=0, reference=None,
                            normal_ordered=True, sorted=True,
                            rohf_canonicalization="Guest-Saunders"):
    """Return ``(system, hamiltonian)`` for ``ccpy.drivers.driver.Driver``.

    ``reference`` defaults to the OpenQP ``[scf] type`` (``rhf``, ``rohf`` or
    ``uhf``). The reference energy is rebuilt from the integrals and must agree
    with the OpenQP SCF energy to ``REFERENCE_ENERGY_TOLERANCE``.
    """
    ref = (reference or mol.config["scf"]["type"]).lower()
    if ref not in ("rhf", "rohf", "uhf"):
        raise ValueError(
            f"CCpy needs an RHF, ROHF or UHF reference; got [scf] type={ref!r}.")
    method = str(mol.config["input"].get("method", "hf")).lower()
    if method != "hf":
        raise ValueError(
            f"CCpy needs the molecule of a method=hf run; got method={method!r}, "
            "whose stored energy is not the SCF energy.")
    if mol.config["input"].get("functional"):
        raise ValueError(
            "CCpy needs Hartree-Fock orbitals; the OpenQP SCF used functional="
            f"{mol.config['input']['functional']!r}.")
    if nfrozen < 0 or ndelete < 0:
        raise ValueError(f"nfrozen={nfrozen} and ndelete={ndelete} must be >= 0.")
    # Everything that can reject the request is checked here, from cheap
    # counts, before CCpy is downloaded or the nbf**4 integrals are built.
    norb, na, nb = _counts(mol)
    if nfrozen > nb:
        raise ValueError(f"nfrozen={nfrozen} exceeds the {nb} doubly occupied "
                         "(beta) orbitals.")
    if ndelete > norb - na:
        raise ValueError(f"ndelete={ndelete} exceeds the {norb - na} virtual "
                         "(alpha) orbitals.")
    e_scf = _scf_energy(mol)
    ensure_ccpy()
    from ccpy.energy.hf_energy import calc_hf_energy, calc_hf_frozen_core_energy
    from ccpy.interfaces.fcidump_tools import build_rohf_fock
    from ccpy.models.integrals import getHamiltonian, getUHFHamiltonian
    from ccpy.models.system import System, UHFSystem

    # eri_ao() stores the nbf**4 AO tensor in the native OQP::ERI_AO record
    # and returns views of it; release the record once CCpy has its own MO
    # integrals, on success and on failure, after every view is gone.
    d = None
    try:
        d = _openqp_ao_data(mol)
        return _build_hamiltonian(
            d, ref, na, nb, norb, nfrozen, ndelete,
            int(mol.config["input"].get("charge", 0) or 0), e_scf,
            normal_ordered, sorted, rohf_canonicalization,
            (calc_hf_energy, calc_hf_frozen_core_energy, build_rohf_fock,
             getHamiltonian, getUHFHamiltonian, System, UHFSystem))
    except BaseException as exc:
        # finished frames in the traceback still hold the borrowed views
        traceback.clear_frames(exc.__traceback__)
        raise
    finally:
        if d is not None:
            d.clear()
        _release_eri_ao(mol)


def _release_eri_ao(mol):
    try:
        del mol.data["OQP::ERI_AO"]
    except KeyError:
        pass                                    # not allocated (failed before eri_ao)


def _build_hamiltonian(d, ref, na, nb, norb, nfrozen, ndelete, charge, e_scf,
                       normal_ordered, sorted, rohf_canonicalization, ccpy_api):
    (calc_hf_energy, calc_hf_frozen_core_energy, build_rohf_fock,
     getHamiltonian, getUHFHamiltonian, System, UHFSystem) = ccpy_api
    nelec, mult = na + nb, na - nb + 1

    if ref == "uhf":
        occ = {"a": np.array([1.0] * na + [0.0] * (norb - na)),
               "b": np.array([1.0] * nb + [0.0] * (norb - nb))}
        system = UHFSystem(
            nelec, norb, mult, nfrozen, ndelete=ndelete, point_group="C1",
            orbital_symmetries={"a": ["A"] * norb, "b": ["A"] * norb},
            charge=charge, nuclear_repulsion=d["vnn"],
            mo_energies={"a": d["ea"], "b": d["eb"]}, mo_occupation=occ)
        e_ref = _uhf_energy(d["hcore"], d["v_ao"], d["ca"], d["cb"], na, nb) + d["vnn"]
        _check_reference("UHF", e_ref, e_scf)
        system.reference_energy = e_ref
        hamiltonian = getUHFHamiltonian(d["ca"], d["cb"], d["hcore"], d["v_ao"],
                                        system, normal_ordered, sorted)
        return system, hamiltonian

    c = d["ca"]                       # RHF and ROHF share one set of spatial MOs
    e1int = np.asfortranarray(c.T @ d["hcore"] @ c)
    e2int = np.asfortranarray(np.einsum("ijkl,ip,jq,kr,ls->pqrs", d["v_ao"], c, c, c, c,
                                        optimize=True))
    system = System(nelec, norb, mult, nfrozen, ndelete=ndelete, point_group="C1",
                    orbital_symmetries=["A"] * norb, charge=charge,
                    nuclear_repulsion=d["vnn"],
                    mo_energies=d["ea"])
    if ref == "rohf":
        system.mo_energies, _ = build_rohf_fock(e1int, e2int, system,
                                                rohf_canonicalization)
    e_ref = calc_hf_energy(e1int, e2int, system) + d["vnn"]
    _check_reference(ref.upper(), e_ref, e_scf)
    system.reference_energy = e_ref
    system.frozen_energy = calc_hf_frozen_core_energy(e1int, e2int, system)
    return system, getHamiltonian(e1int, e2int, system, normal_ordered, sorted)


def driver_from_openqp(mol, nfrozen=0, ndelete=0, reference=None, **kwargs):
    """Return a ``ccpy.drivers.driver.Driver`` for a converged OpenQP SCF.

    For a UHF reference the driver runs with ``RHF_symmetry=False``. Build the
    ``Driver`` through this function rather than from
    :func:`hamiltonian_from_openqp` directly, or set that option yourself.
    """
    ref = (reference or mol.config["scf"]["type"]).lower()
    # hamiltonian_from_openqp validates the request before it installs CCpy,
    # so an invalid request never downloads anything.
    system, hamiltonian = hamiltonian_from_openqp(
        mol, nfrozen=nfrozen, ndelete=ndelete, reference=reference, **kwargs)
    from ccpy.drivers.driver import Driver
    driver = Driver(system, hamiltonian)
    if ref == "uhf":
        # CCpy turns RHF_symmetry on whenever the alpha and beta occupied
        # counts are equal and then copies alpha amplitudes into beta ones.
        # A UHF singlet (broken symmetry) has different alpha and beta
        # orbitals, so that shortcut would be wrong.
        driver.options["RHF_symmetry"] = False
    return driver
