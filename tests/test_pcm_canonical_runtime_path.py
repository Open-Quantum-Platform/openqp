"""Source-level invariants for the single canonical PCM (ddX) runtime path.

These checks were previously recorded as prose in docs/solvent_pcm_*.md and
docs/solvent_ddx_scf_integration_seam.md (a "validation matrix" / "integration
seam" describing the intended design). The documents have been removed; the
enforceable claims they carried are kept here as tests that assert against the
*actual source code* instead of against markdown — so the invariants are
executable rather than narrative.

Pure text/AST-free assertions over the committed sources: no build, no ddX, no
Python OpenQP import required, so this runs in every environment.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def _read(*parts):
    return (ROOT.joinpath(*parts)).read_text(encoding="utf-8")


class PcmCanonicalRuntimePathTests(unittest.TestCase):
    """One gate, one Fock hook, one energy field, one backend."""

    def test_single_canonical_runtime_path(self):
        # Gate: a single control flag toggles the whole feature.
        types = _read("source", "types.F90")
        self.assertIn("pcm_enabled", types)

        # Hook + energy field: the public reaction-field entry point sets e_pcm.
        pcm = _read("source", "solvent_pcm.F90")
        self.assertIn("public :: add_pcm_reaction_field", pcm)
        self.assertIn("subroutine add_pcm_reaction_field", pcm)
        self.assertIn("e_pcm", pcm)

        # SCF wiring: gated by pcm_enabled, applied via add_pcm_reaction_field,
        # reported in E%e_pcm.
        scf = _read("source", "scf_addons.F90")
        self.assertIn("infos%control%pcm_enabled", scf)
        self.assertIn("add_pcm_reaction_field", scf)
        self.assertIn("E%e_pcm", scf)

    def test_no_second_runtime_reaction_field_hook(self):
        # The retired "reference-supplied-potential" prototype was removed: the
        # Fock reaction field is applied at exactly one site in the SCF builder.
        scf = _read("source", "scf_addons.F90")
        self.assertEqual(
            scf.count("call add_pcm_reaction_field"),
            1,
            "expected exactly one PCM reaction-field call site (single canonical path)",
        )

    def test_pcm_energy_convention(self):
        # e_pcm = -1/2 <phi_cav, q_cav>: the surface charges contracted with the
        # exact total solute potential, with the canonical -0.5 polarization factor.
        pcm = _read("source", "solvent_pcm.F90")
        self.assertIn("PCM_QCAV_TO_FOCK_SCALE = -0.5_dp", pcm)
        self.assertIn(
            "e_pcm = PCM_QCAV_TO_FOCK_SCALE * dot_product(phi_cav, q_cav)", pcm
        )

    def test_provisional_conventions_documented_in_source(self):
        # The sign/scale conventions that are not yet validated against a trusted
        # PCM reference are recorded next to the code that relies on them.
        pcm = _read("source", "solvent_pcm.F90")
        self.assertIn("PROVISIONAL CONVENTIONS", pcm)
        self.assertIn("phi_cav sign", pcm)
        self.assertIn("ddx_get_xi", pcm)         # q_cav from ddX adjoint charge
        self.assertIn("q_cav sign/scale", pcm)

    def test_ddx_model_discretization_defaults(self):
        # ddPCM (model=2), lmax=8, 302-point Lebedev grid, eta=0.1, water eps.
        adapter = _read("source", "solvent_ddx_adapter.c")
        self.assertIn("const int model_pcm = 2;", adapter)   # ddPCM
        self.assertIn("const int lmax = 8;", adapter)
        self.assertIn("const int n_lebedev = 302;", adapter)
        self.assertIn("const double eta = 0.1;", adapter)
        self.assertIn("const double epsilon = 78.3553;", adapter)

    def test_energy_only_rhf_rohf_scope_recorded(self):
        # First-scope: RHF/ROHF reference single-point energy; gradients and
        # state-specific MRSF PCM are out of scope. Enforced by the input checker.
        checker = _read("pyoqp", "oqp", "utils", "input_checker.py")
        self.assertIn("energy-only solvent backend", checker)
        self.assertIn("RHF/ROHF", checker)
        self.assertIn("PCM first scope supports RHF/ROHF", checker)


if __name__ == "__main__":
    unittest.main()
