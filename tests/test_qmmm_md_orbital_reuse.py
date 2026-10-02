"""QM/MM MD must apply the public orbital-reuse control after success."""

import pytest

try:
    from oqp.library.qmmm_md import QMMM_MD
    from oqp.library.qmmm_driver import OpenQpQMMM
    import oqp.library.qmmm_driver as qmmm_driver_module
except (ImportError, RuntimeError) as error:
    pytest.skip(f"OpenMM/OpenQP runtime required: {error}", allow_module_level=True)


class _ForceDriver:
    def __init__(self, *, fail=False):
        self._reuse_orbitals = False
        self.fail = fail
        self.reuse_at_call = []

    def compute_force(self, positions, topology, mm_systems, qm_atoms):
        self.reuse_at_call.append(self._reuse_orbitals)
        if self.fail:
            raise RuntimeError("SCF failure")
        return "energy", "force"


def _bare_qmmm_md(driver, *, mo_reuse):
    dynamics = object.__new__(QMMM_MD)
    dynamics.oqp_driver = driver
    dynamics.mo_reuse = mo_reuse
    dynamics.pdb = type("PDB", (), {"topology": "topology"})()
    dynamics.mm_systems = "systems"
    dynamics.qm_atoms = "qm-atoms"
    return dynamics


@pytest.mark.parametrize("mo_reuse", [True, False])
def test_orbital_reuse_starts_only_after_a_successful_force(mo_reuse):
    force_driver = _ForceDriver()
    dynamics = _bare_qmmm_md(force_driver, mo_reuse=mo_reuse)

    assert dynamics._compute_qmmm_force("positions") == ("energy", "force")
    assert force_driver.reuse_at_call == [False]
    assert force_driver._reuse_orbitals is mo_reuse

    dynamics._compute_qmmm_force("next-positions")
    assert force_driver.reuse_at_call == [False, mo_reuse]


def test_failed_force_disables_orbital_reuse_before_resume():
    force_driver = _ForceDriver(fail=True)
    force_driver._reuse_orbitals = True
    dynamics = _bare_qmmm_md(force_driver, mo_reuse=True)

    with pytest.raises(RuntimeError, match="SCF failure"):
        dynamics._compute_qmmm_force("positions")
    assert force_driver.reuse_at_call == [True]
    assert force_driver._reuse_orbitals is False


def test_config_mode_moves_retained_molecule_and_rebuilds_integrals(monkeypatch):
    calls = []
    mol = type("Mol", (), {"config": {"input": {"method": "hf"}}})()
    driver = object.__new__(OpenQpQMMM)
    driver.op = type("Op", (), {"mol": mol})()
    driver._image_warm = False
    driver._reuse_orbitals = True
    driver._update_mol_positions = lambda target: calls.append(("move", target))
    monkeypatch.setattr(
        qmmm_driver_module.oqp.library,
        "set_basis",
        lambda target: calls.append(("basis", target)),
    )
    monkeypatch.setattr(
        qmmm_driver_module,
        "ints_1e",
        lambda target: calls.append(("ints", target)),
    )

    assert driver._prepare_config_continuation() is True
    assert calls == [("move", mol), ("basis", mol), ("ints", mol)]


def test_config_mode_fresh_guess_path_does_not_claim_continuation():
    driver = object.__new__(OpenQpQMMM)
    driver.op = type(
        "Op", (),
        {"mol": type("Mol", (), {"config": {"input": {"method": "hf"}}})()},
    )()
    driver._image_warm = False
    driver._reuse_orbitals = False
    assert driver._prepare_config_continuation() is False
