"""OpenQP runtime adapter for the optional OpenQP-xTB (LC-GFN1-xTB) backend.

Thin subclass of the OpenQP-DFTB adapter: the two standalone libraries share
the same C API family, tag names (OQP::dftb_wf_dims, OQP::VEC_MO_A,
OQP::td_bvec_mo, ...), and the mol.dftb_external_potential QM/MM embedding
attribute (see oqp.utils.tb_backends). Differences, per the pinned C ABI v3 of
docs/design/gfn1-implementation-spec.md (sections 4, 8a, 8b of the openqp-xtb
design spec):

* [xtb] input section, openqp_xtb_* symbols in libopenqp_xtb_c.{dylib,so};
* lc_gamma gains the 'ok' (Ohno-Klopman) kind: {yukawa:0, erf:1, ok:2},
  default 'ok' (response exchange defaults: omega=0.3, cam_alpha=0,
  cam_beta=1);
* five model scalars inserted by value immediately AFTER the zvector flag and
  BEFORE scc_tolerance: int64 model (0=dftb2, 1=gfn1), int64 dispersion,
  int64 halogen_bond, int64 third_order, double spin_scale;
* openqp_xtb_states_overlap / openqp_xtb_soc_matrix keep the v2 signatures.
"""

from __future__ import annotations

import ctypes

from oqp.library.openqp_dftb import OpenQPDFTBAdapter, _tb_soc


class OpenQPXTBAdapter(OpenQPDFTBAdapter):
    """Make OpenQP-xTB look like a normal OpenQP energy/gradient provider."""

    SECTION = "xtb"
    SYMBOL_PREFIX = "openqp_xtb"
    LIB_BASENAMES = ("libopenqp_xtb_c.dylib", "libopenqp_xtb_c.so")
    ENV_LIBRARY = "OPENQP_XTB_LIBRARY"
    ENV_PARAMETER = "OPENQP_XTB_PARAMETER_PATH"
    PIP_LOCATOR = "openqp_xtb"
    BACKEND_NAME = "openqp-xtb"
    DISPLAY_NAME = "OpenQP-xTB"
    CACHE_ATTR = "_openqp_xtb_cache"

    # C ABI model codes ([xtb] model key; only gfn1 is exposed for now).
    _MODEL_CODES = {"gfn1": 1}
    _LC_GAMMA_CODES = {"yukawa": 0, "erf": 1, "ok": 2}

    # openqp-xtb C ABI. Its first released layout (openqp_xtb_abi_version == 4)
    # already carries the model block, so an unversioned probe resolves to the
    # current layout rather than the legacy DFTB v1. The state-gradient
    # argument order is: the shared common block (with the five model scalars
    # spliced after the zvector flag via _model_args), then the DTCAM-TB
    # response knobs c_mrsf/response_global_hybrid/onsite_exchange_scale, then
    # the shared output tail -- i.e. the same trailing block as DFTB ABI v1.
    _SUPPORTED_ABI = (1, 2, 3, 4)
    _UNVERSIONED_ABI_FALLBACK = 4

    def _abi_version_symbol(self) -> str:
        # openqp-xtb exports openqp_xtb_abi_version (no "capi" infix).
        return "openqp_xtb_abi_version"

    def _call_state_gradient_current(self, function, **values) -> None:
        # openqp_xtb_state_gradient == common(+model block) + the three
        # DTCAM-TB response knobs + output tail, which is exactly the DFTB
        # ABI-v1 argument shape. Reuse it rather than the DFTB "current" layout
        # (which appends preset/onsite/w_scale scalars openqp-xtb does not take).
        self._call_state_gradient_abi1(function, **values)

    def _install_state_gradient_argtypes(self, lib, abi_version) -> None:
        # The xTB layout inserts the model block, so the DFTB argtypes table
        # does not apply. Leave the symbol untyped -- the adapter always passes
        # fully-typed ctypes objects -- and only pin the void return type.
        getattr(lib, f"{self.SYMBOL_PREFIX}_state_gradient").restype = None

    def _lc_gamma_code(self) -> int:
        """lc_gamma_kind for the C ABI: 0=yukawa, 1=erf, 2=ok (default)."""
        kernel = str(self.dftb.get("lc_gamma", "ok")).lower()
        if kernel not in self._LC_GAMMA_CODES:
            raise ValueError(
                f"[xtb] lc_gamma must be yukawa, erf, or ok, got {kernel!r}")
        return self._LC_GAMMA_CODES[kernel]

    def _model_code(self) -> int:
        model = str(self.dftb.get("model", "gfn1")).strip().lower()
        if model not in self._MODEL_CODES:
            raise ValueError(
                f"[xtb] model must be gfn1 for now, got {model!r}")
        return self._MODEL_CODES[model]

    def _model_args(self) -> list:
        """The five C ABI v3 scalars, spliced right BEFORE scc_tolerance.

        Order pinned by the spec: int64 model (0=dftb2, 1=gfn1), int64
        dispersion, int64 halogen_bond, int64 third_order, double spin_scale.
        """
        spin_scale = float(self.dftb.get("spin_scale", 1.0))
        return [
            ctypes.c_int64(self._model_code()),
            ctypes.c_int64(int(bool(self.dftb.get("dispersion", True)))),
            ctypes.c_int64(int(bool(self.dftb.get("halogen_bond", True)))),
            ctypes.c_int64(int(bool(self.dftb.get("third_order", True)))),
            ctypes.c_double(spin_scale),
        ]

    def _run_probe(self, method: str, state: int):
        raise ValueError(
            "[xtb] backend=probe is not supported: the probe executable is a "
            "DFTB-only developer fallback. Use [xtb] backend=native (loads the "
            "standalone libopenqp_xtb_c shared library)."
        )

    def _cache_key(self, method: str, state: int, *, need_grad: bool):
        atoms, coords_bytes, option_key = super()._cache_key(
            method, state, need_grad=need_grad)
        option_key = option_key + (
            str(self.dftb.get("model", "gfn1")).strip().lower(),
            bool(self.dftb.get("dispersion", True)),
            bool(self.dftb.get("halogen_bond", True)),
            bool(self.dftb.get("third_order", True)),
            float(self.dftb.get("spin_scale", 1.0)),
            # base key records lc_gamma with the dftb default; re-record with
            # the xtb default 'ok' so unset and explicit 'ok' hash identically.
            str(self.dftb.get("lc_gamma", "ok")).lower(),
        )
        return atoms, coords_bytes, option_key


def xtb_soc(mol):
    """Standalone SOC driver for method=xtb (mirror of dftb_soc).

    Runs the MRSF response for both target multiplicities from the same ROKS
    reference, builds the one-center SOC matrix (xi_l sourced from the `soc`
    records of the .opxtb parameter file), diagonalizes in numpy, and fills
    the same OQP::soc_* tags the native soc_mrsf produces.
    """
    return _tb_soc(mol, OpenQPXTBAdapter)
