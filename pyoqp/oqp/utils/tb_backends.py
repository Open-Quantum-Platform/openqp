"""Shared dispatch helpers for the tight-binding (TB) backends.

OpenQP currently ships two decoupled TB backends behind the same PyOQP
plumbing:

* ``method=dftb`` -> the [dftb] section and the openqp-dftb library
  (``oqp.library.openqp_dftb.OpenQPDFTBAdapter``);
* ``method=xtb``  -> the [xtb] section and the openqp-xtb library
  (``oqp.library.openqp_xtb.OpenQPXTBAdapter``).

Every ``method == 'dftb'`` dispatch site in single_point.py / namd.py /
qmmm_driver.py / runfunc.py goes through the helpers below so both backends
share one code path.

INTENTIONALLY SHARED between the two backends (do not fork per backend):

* the wavefunction/response tag names the adapters publish and the drivers
  consume -- ``OQP::dftb_wf_dims``, ``OQP::VEC_MO_A``, ``OQP::VEC_MO_B``,
  ``OQP::E_MO_A``, ``OQP::E_MO_B``, ``OQP::td_bvec_mo``, ``OQP::td_energies``,
  ``OQP::partial_charges``, the ``OQP::soc_*`` family, ...;
* the ``mol.dftb_external_potential`` attribute that carries the per-atom
  QM/MM embedding potential into the SCC Hamiltonian.

Both libraries implement the same C API family, so the drivers stay agnostic:
only the section name, the shared-library/symbol prefix, and a few model
options differ (see the adapter classes).
"""

from __future__ import annotations


TB_METHODS = {"dftb", "xtb"}


def is_tb_method(method_or_mol) -> bool:
    """True when the input method is one of the tight-binding backends.

    Accepts either the method name itself (any case) or an object carrying a
    ``config`` dict (a Molecule), so both ``is_tb_method(self.method)`` and
    ``is_tb_method(mol)`` read naturally at the dispatch sites.
    """
    if isinstance(method_or_mol, str):
        return method_or_mol.strip().lower() in TB_METHODS
    config = getattr(method_or_mol, "config", None)
    if config is None:
        return False
    return str(config.get("input", {}).get("method", "")).strip().lower() in TB_METHODS


def tb_section_name(config) -> str:
    """Name of the TB options section for this config: 'dftb' or 'xtb'.

    The section name equals the input method by construction ([dftb] for
    method=dftb, [xtb] for method=xtb).
    """
    method = str(config.get("input", {}).get("method", "")).strip().lower()
    if method not in TB_METHODS:
        raise ValueError(
            f"tb_section_name called for non-TB method {method!r}; "
            f"expected one of {sorted(TB_METHODS)}"
        )
    return method


def tb_config(config) -> dict:
    """The TB options section ([dftb] or [xtb]) of this config, {} by default."""
    return config.get(tb_section_name(config), {})


def make_tb_adapter(mol):
    """Instantiate the adapter matching mol's input method (dftb or xtb)."""
    section = tb_section_name(mol.config)
    if section == "xtb":
        from oqp.library.openqp_xtb import OpenQPXTBAdapter  # noqa: PLC0415
        return OpenQPXTBAdapter(mol)
    from oqp.library.openqp_dftb import OpenQPDFTBAdapter  # noqa: PLC0415
    return OpenQPDFTBAdapter(mol)
