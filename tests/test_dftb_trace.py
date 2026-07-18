"""Unit tests for the openqp-dftb structured-trace parser and formatters."""

import importlib.util
from pathlib import Path

MODULE = Path(__file__).parents[1] / "pyoqp" / "oqp" / "utils" / "dftb_trace.py"
SPEC = importlib.util.spec_from_file_location("oqp_dftb_trace", MODULE)
dftb_trace = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dftb_trace)


LEVEL2_TRACE = """
openqp_dftb_scc_pass rhf_seed pass 1 stage seed mixer broyden elapsed_seconds  0.000000E+000
openqp_dftb_scc_iter rhf_seed 1 residual  2.500000E-002 tol  1.000000E-005 pass 1 elapsed_seconds  1.000000E-004
openqp_dftb_scc_iter rhf_seed 2 residual  3.100000E-006 tol  1.000000E-005 pass 1 elapsed_seconds  2.000000E-004
openqp_dftb_scc_pass_end rhf_seed pass 1 outcome converged iterations 2 elapsed_seconds  2.000000E-003
openqp_dftb_scc_pass rohf pass 1 stage initial mixer broyden elapsed_seconds  0.000000E+000
openqp_dftb_scc_iter rohf 1 residual  1.200000E-002 charge  1.200000E-002 spin  3.000000E-003 tol  1.000000E-008
openqp_dftb_scc_iter rohf_trah 2 residual  4.000000E-006 charge  4.000000E-006 spin  1.000000E-007 tol  1.000000E-008
openqp_dftb_scc_iter rohf 3 residual  5.000000E-009 charge  5.000000E-009 spin  2.000000E-011 tol  1.000000E-008
openqp_dftb_scc_pass_end rohf pass 1 outcome converged iterations 3 elapsed_seconds  9.000000E-003
openqp_dftb_scc_post_update rohf residual  6.434884E-010 charge  6.434884E-010 spin  4.727729E-011 density  3.197395E-010 tol  1.000000E-008 converged T elapsed_seconds  1.000000E-002
openqp_dftb_scc_recovery rohf kind level_shift attempted F succeeded F elapsed_seconds  1.000000E-002
openqp_dftb_scc_recovery rohf kind charge_spin_trah attempted F succeeded F elapsed_seconds  1.000000E-002
openqp_dftb_davidson_start dimension 143 roots 3 max_subspace 100 max_iterations 50 tolerance  1.000000E-006 elapsed_seconds  0.000000E+000
openqp_dftb_davidson_iter iteration 1 basis 3 converged_roots 0 residual  1.500000E-003 matvecs 3 restarts 0 elapsed_seconds  1.000000E-004
openqp_dftb_davidson_iter iteration 2 basis 6 converged_roots 2 residual  8.000000E-006 matvecs 6 restarts 0 elapsed_seconds  3.000000E-004
openqp_dftb_davidson_end outcome converged iterations 7 residual  6.480064E-007 matvecs 20 restarts 0 elapsed_seconds  1.000000E-003

DFTB state-to-state oscillator strengths
  Approximation: unrelaxed TDA/state interaction, length gauge, DFTB Mulliken-monopole transition dipoles
  Initial  Final      Delta E (eV)      |mu| (a.u.)                f
  root 1   root 2         3.984014         0.716700     9.990000E-002
  root 1   root 3         7.942290         0.124000     3.300000E-003
  root 2   root 3         3.958276         0.500000     2.500000E-002
"""

GRADIENT_TRACE = """
openqp_dftb_scc_pass rohf pass 1 stage initial mixer broyden elapsed_seconds  0.000000E+000
openqp_dftb_scc_iter rohf 1 residual  1.200000E-002 charge  1.200000E-002 spin  3.000000E-003 tol  1.000000E-008
openqp_dftb_scc_pass_end rohf pass 1 outcome converged iterations 12 elapsed_seconds  9.000000E-003
openqp_dftb_zvector_start dimension 286 max_iterations 50 tolerance  1.000000E-010 elapsed_seconds  0.000000E+000
[zvec-pair] using pair-space single solve, gmres_iters 11 fixed_point_max_abs  3.200000E-012 preconditioned_l2  1.100000E-013
openqp_dftb_zvector_iter iteration 1 residual  1.000000E-004 inner_solves 2 inner_matvecs 14 elapsed_seconds  1.000000E-002
openqp_dftb_zvector_iter iteration 2 residual  9.000000E-011 inner_solves 2 inner_matvecs 12 elapsed_seconds  2.000000E-002
openqp_dftb_zvector_end outcome pair_space_converged pair_iterations 11 anderson_iterations 2 fallback_gmres_iterations 0 residual_kind fixed_point_max_abs residual  3.200000E-012 elapsed_seconds  8.000000E-001
"""

PAIR_SOLVE_TRACE = """
openqp_dftb_zvector_start dimension 968 max_iterations 2000 tolerance  1.000000E-009 elapsed_seconds  0.000000E+000
openqp_dftb_gmres_start dimension 240 max_iterations 2000 max_subspace 100 tolerance  1.000000E-009 elapsed_seconds  0.000000E+000
openqp_dftb_gmres_iter iteration 0 cycle 1 residual  4.100000E-002 matvecs 1 elapsed_seconds  1.000000E-004
openqp_dftb_gmres_iter iteration 1 cycle 1 residual  6.300000E-004 matvecs 2 elapsed_seconds  2.000000E-004
openqp_dftb_gmres_iter iteration 2 cycle 1 residual  8.800000E-007 matvecs 3 elapsed_seconds  3.000000E-004
openqp_dftb_gmres_iter iteration 3 cycle 1 residual  2.100000E-010 matvecs 4 elapsed_seconds  4.000000E-004
openqp_dftb_gmres_end outcome converged iterations 3 residual  2.100000E-010 matvecs 5 elapsed_seconds  5.000000E-004
openqp_dftb_zvector_end outcome pair_space_converged pair_iterations 3 anderson_iterations 0 fallback_gmres_iterations 0 residual_kind fixed_point_max_abs residual  2.650000E-013 elapsed_seconds  1.000000E-003
"""

COMPONENT_SPECTRUM_TRACE = """
openqp_dftb_energy_components total -1.04652311259000E+001 electronic -1.07000000000000E+001 h0 -1.05000000000000E+001 scc -2.00000000000000E-001 spin  0.00000000000000E+000 exchange  0.00000000000000E+000 repulsive  2.34768874100000E-001 band -5.60000000000000E+000

DFTB state-to-state oscillator strengths
  Approximation: unrelaxed TDA/state interaction, length gauge, DFTB Mulliken-monopole transition dipoles
  Initial  Final      Delta E (eV)      |mu| (a.u.)                f    mu_x (a.u.)    mu_y (a.u.)    mu_z (a.u.)
  root 1   root 2         6.951845         2.504000     1.067900E+000       0.000000       2.504000       0.000000
  root 1   root 3         7.376022         0.000000     0.000000E+000       0.000000       0.000000       0.000000
"""


def test_parse_level2_trace_collects_iterations_and_spectrum():
    trace = dftb_trace.parse_native_trace(LEVEL2_TRACE)

    assert [entry["label"] for entry in trace["scc_passes"]] == \
        ["rhf_seed", "rohf"]
    seed, rohf = trace["scc_passes"]
    assert len(seed["iterations"]) == 2
    assert seed["outcome"] == "converged"
    assert len(rohf["iterations"]) == 3
    assert rohf["iterations"][1]["mode"] == "TRAH"
    assert rohf["iterations"][2]["charge"] == 5.0e-9

    assert trace["post_updates"][0]["converged"] is True
    assert all(not entry["attempted"] for entry in trace["recoveries"])

    davidson = trace["davidson"][0]
    assert davidson["dimension"] == 143
    assert len(davidson["iterations"]) == 2
    assert davidson["outcome"] == "converged"
    assert davidson["iteration_count"] == 7

    spectrum = trace["spectrum"]
    assert spectrum["rows"][0]["initial"] == "root 1"
    assert spectrum["rows"][0]["oscillator"] == 9.99e-2
    assert len(spectrum["rows"]) == 3


def test_scc_table_shows_iterations_and_convergence():
    trace = dftb_trace.parse_native_trace(LEVEL2_TRACE)
    text = dftb_trace.format_scc_block(trace)

    assert "SCC iterations: restricted charge seed" in text
    assert "SCC iterations: open-shell (ROKS) reference" in text
    assert "Iter" in text and "Residual" in text
    assert "TRAH" in text
    assert "SCC converged in 3 iterations" in text
    assert "Final self-consistency" in text
    # Recoveries that were never attempted stay out of the default log.
    assert "recovery" not in text


def test_davidson_table_shows_iterations():
    trace = dftb_trace.parse_native_trace(LEVEL2_TRACE)
    text = dftb_trace.format_davidson_block(
        trace["davidson"][0], "MRSF-TDDFTB")

    assert "Davidson algorithm for MRSF-TDDFTB response states" in text
    assert "dimension = 143" in text
    assert "Converged roots" in text
    assert "Davidson converged in 7 iterations" in text
    assert "20 matvecs" in text


def test_zvector_table_shows_iterations_and_outcome():
    trace = dftb_trace.parse_native_trace(GRADIENT_TRACE)
    text = dftb_trace.format_zvector_block(trace["zvector"][0])

    assert "Z-vector relaxed-density solve" in text
    assert "dimension = 286" in text
    assert "Inner solves" in text
    assert "pair-space single solve" in text
    assert "residual 3.20E-12" in text
    # The bracketed debug line is not a structured record.
    assert "[zvec-pair]" in "\n".join(trace["other"])


def test_gradient_summary_lines_are_compact():
    trace = dftb_trace.parse_native_trace(GRADIENT_TRACE)
    text = dftb_trace.format_scc_summary_lines(trace)

    assert "SCC open-shell (ROKS) reference: converged in 12 iterations" in text


def test_spectrum_root_labels_map_to_public_states():
    assert dftb_trace.map_spectrum_label("root 1", "S") == "S0"
    assert dftb_trace.map_spectrum_label("root 3", "T") == "T2"
    assert dftb_trace.map_spectrum_label("S2", "S") == "S2"
    assert dftb_trace.map_spectrum_label("root 2", None) == "root 2"


def test_excited_state_summary_table_and_transitions():
    trace = dftb_trace.parse_native_trace(LEVEL2_TRACE)
    spectrum = trace["spectrum"]
    oscillator = dftb_trace.spectrum_from_ground(spectrum, "S")
    summary = {
        "method": "mrsf",
        "state_labels": ["S0", "S1", "S2"],
        "state_energies": [-114.43, -114.28, -114.13],
        "spin_squares": [0.0, 0.0, 0.0],
        "oscillator": oscillator,
        "transitions": [{
            "initial": dftb_trace.map_spectrum_label(row["initial"], "S"),
            "final": dftb_trace.map_spectrum_label(row["final"], "S"),
            "de_ev": row["de_ev"],
            "dipole": row["dipole"],
            "oscillator": row["oscillator"],
        } for row in spectrum["rows"]],
        "approximation": spectrum["approximation"],
        "reference_label": "Reference: high-spin ROKS (internal)",
        "reference_energy": -114.32,
    }
    text = dftb_trace.format_excited_state_summary(summary)

    assert "Summary table" in text
    assert "Excitation" in text and "Oscillator" in text
    lines = [line for line in text.splitlines() if line.strip().startswith("S1")]
    assert lines and "0.7167" in lines[0] and "0.0999" in lines[0]
    assert "Reference: high-spin ROKS (internal)" in text
    assert "State-to-state transitions" in text
    assert text.index("Summary table") < text.index("State-to-state transitions")
    transition_lines = [line for line in text.splitlines()
                        if line.strip().startswith("S1 ")
                        or "S1        S2" in line.replace("   ", " ")]
    assert "S2" in text


def test_final_energy_annotation_alignment():
    summary = {
        "method": "mrsf",
        "state_labels": ["S0", "S1"],
        "state_energies": [-114.43, -114.28],
        "oscillator": {
            "S1": {"dipole": 0.7167, "oscillator": 0.0999},
        },
    }
    # Row 0 is the internal reference: no annotation.
    assert dftb_trace.final_energy_annotation(summary, 0) == ""
    # Column names live in the header, not in every row.
    header = dftb_trace.final_energy_header(summary)
    assert "VEE (eV)" in header and "|mu| (a.u.)" in header
    s0 = dftb_trace.final_energy_annotation(summary, 1)
    assert "VEE" not in s0 and s0.strip() == "0.000000"
    s1 = dftb_trace.final_energy_annotation(summary, 2)
    assert "VEE" not in s1
    assert "0.7167" in s1 and "0.0999" in s1
    assert dftb_trace.final_energy_annotation(summary, 3) == ""
    assert dftb_trace.final_energy_annotation(None, 1) == ""
    assert dftb_trace.final_energy_header(None) == ""


def test_final_energy_labels_for_closed_shell_tddftb():
    summary = {"method": "tddftb", "state_labels": ["S0", "S1", "S2"]}
    assert dftb_trace.final_energy_label(summary, 0, "state 0") == "S0"
    assert dftb_trace.final_energy_label(summary, 1, "state 1") == "S1"
    assert dftb_trace.final_energy_label(summary, 9, "state 9") == "state 9"
    mrsf = {"method": "mrsf", "state_labels": ["S0"]}
    # MRSF rows keep public_state_label (row 0 is the internal reference).
    assert dftb_trace.final_energy_label(mrsf, 0, "reference") == "reference"


def test_zvector_pair_solve_iterations_are_rendered():
    trace = dftb_trace.parse_native_trace(PAIR_SOLVE_TRACE)
    solve = trace["zvector"][0]

    assert not solve["iterations"]
    assert len(solve["krylov_solves"]) == 1
    assert len(solve["krylov_solves"][0]["iterations"]) == 4

    text = dftb_trace.format_zvector_block(solve)
    assert "Pair-space GMRES iterations" in text
    assert "2.100000E-10" in text
    assert "pair-space single solve" in text


def test_gmres_outside_zvector_stays_generic():
    trace = dftb_trace.parse_native_trace(
        "openqp_dftb_gmres_start dimension 10 max_iterations 5 "
        "max_subspace 5 tolerance  1.000000E-008 elapsed_seconds  0.0\n"
        "openqp_dftb_gmres_end outcome converged iterations 2 residual "
        " 1.000000E-009 matvecs 3 elapsed_seconds  1.000000E-003\n")
    assert len(trace["linear"]) == 1
    assert not trace["zvector"]


def test_energy_components_block_and_spectrum_components():
    trace = dftb_trace.parse_native_trace(COMPONENT_SPECTRUM_TRACE)

    components = trace["energy_components"]
    assert components["total"] == -10.4652311259
    text = dftb_trace.format_energy_components(components)
    assert "Core Hamiltonian (H0) energy" in text
    assert "TOTAL energy" in text
    # Zero-valued optional terms stay out of the block.
    assert "Spin polarization" not in text
    assert "Long-range exchange" not in text

    rows = trace["spectrum"]["rows"]
    assert rows[0]["dipole_xyz"] == [0.0, 2.504, 0.0]
    assert rows[0]["dipole"] == 2.504
    assert rows[0]["oscillator"] == 1.0679


RESOLVED_OPTIONS_TRACE = (
    "openqp_dftb_resolved_options preset dtcam-tb lc_ground_state T "
    "lc_gamma_erf T omega  3.00000000E-001 cam_alpha  0.00000000E+000 "
    "cam_beta  4.00000000E-002 response_omega  2.62500000E-001 "
    "response_cam_alpha  0.00000000E+000 response_cam_beta  1.01250000E+000 "
    "w_scale  1.00000000E+000 response_w_scale  6.37500000E-001 "
    "spc_coco  1.02500000E+000 spc_ovov  2.50000000E-001 "
    "spc_coov  2.62500000E-001 onsite_ss  0.00000000E+000 "
    "onsite_sp  0.00000000E+000 onsite_pp -1.25000000E-002 "
    "onsite_exchange_scale  0.00000000E+000 c_mrsf -1.00000000E+000 "
    "c_mrsf_oo -1.00000000E+000 response_global_hybrid F "
    "global_hybrid_exchange F builtin_spin_constants F scc_mixer 2 "
    "scc_history 12 scc_mixing  3.50000000E-001 scc_max_step "
    " 1.00000000E+000 scc_tolerance  1.00000000E-008 max_scc_iterations "
    "4000 response_solver 1 response_tolerance  1.00000000E-006 "
    "response_max_iterations 50 response_max_subspace 100 zvector T\n"
)


def test_resolved_options_record_expands_preset_values():
    trace = dftb_trace.parse_native_trace(RESOLVED_OPTIONS_TRACE)
    options = trace["resolved_options"]

    assert options["preset"] == "dtcam-tb"
    assert options["lc_ground_state"] is True
    assert options["spc_coco"] == 1.025
    assert options["scc_mixer"] == 2
    assert options["max_scc_iterations"] == 4000

    text = dftb_trace.format_resolved_options(options)
    assert "preset: dtcam-tb" in text
    assert "erf kernel" in text
    assert "0.3000 / 0.0000 / 0.0400" in text
    assert "0.2625 / 0.0000 / 1.0125" in text
    assert "1.0250 / 0.2500 / 0.2625" in text
    assert "broyden" in text and "davidson" in text
    assert "Analytic Z-vector gradient:     yes" in text


def test_timing_footer_format():
    text = dftb_trace.format_step_timing(0.012, 0.013, 0.5, 0.7)
    assert "Step  CPU time (seconds):" in text
    assert "Total CPU time (seconds):" in text
    assert "0.013" in text and "0.700" in text


def test_spin_flip_pair_mapping_is_column_major():
    # amplitudes(noca, n_virt_beta) flattened column-major:
    # index = (vir - nocb - 1) * noca + occ.
    noca, nocb = 12, 10
    assert dftb_trace.spin_flip_pair(1, noca, nocb) == (1, 11)
    assert dftb_trace.spin_flip_pair(12, noca, nocb) == (12, 11)
    assert dftb_trace.spin_flip_pair(13, noca, nocb) == (1, 12)
    assert dftb_trace.spin_flip_pair(26, noca, nocb) == (2, 13)


def test_state_configurations_render_before_summary_table():
    summary = {
        "method": "mrsf",
        "state_labels": ["S0", "S1"],
        "state_energies": [-10.4652311259, -10.2097555403],
        "spin_squares": [0.009, 0.0],
        "oscillator": {},
        "transitions": [],
        "configurations": {
            "S0": [{"index": 131, "coeff": 0.981013, "occ": 11, "vir": 20}],
            "S1": [
                {"index": 11, "coeff": -0.171425, "occ": 11, "vir": 11},
                {"index": 132, "coeff": 0.953585, "occ": 12, "vir": 20},
            ],
        },
        "reference_energy": -10.0480223817,
        "reference_label": "Reference: high-spin ROKS (internal)",
    }
    text = dftb_trace.format_excited_state_summary(summary)

    assert "Spin-adapted spin-flip excitations" in text
    assert text.index("Spin-adapted spin-flip excitations") \
        < text.index("Summary table")
    assert "0.981013" in text and "11  ->    20" in text
    assert "-0.171425" in text
    assert "<S^2> =" in text and "VEE =" in text
    # States with no dominant configurations are simply omitted.
    summary["configurations"] = {"S0": [], "S1": []}
    assert "Spin-adapted" not in dftb_trace.format_excited_state_summary(summary)


def test_charge_block_lists_atoms_and_dipole():
    text = dftb_trace.format_charge_block(
        ["C", "H"], [-0.1, 0.1], [0.0, 0.2, 0.0],
        "Mulliken atomic charges (SCC reference)")
    assert "Mulliken atomic charges (SCC reference)" in text
    assert "C" in text and "H" in text
    assert "Sum" in text
    assert "Debye" in text


def test_unstructured_text_yields_no_tables():
    trace = dftb_trace.parse_native_trace("hello world\nsome text\n")
    assert not trace["scc_passes"]
    assert not trace["davidson"]
    assert trace["other"] == ["hello world", "some text"]
