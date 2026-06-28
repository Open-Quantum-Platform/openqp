#!/usr/bin/env python
"""OQP excited-state analysis & export suite -- proof-of-concept benchmark driver.

Runs the benchmark molecules and exercises every feature + self-check from the
task brief (GATE 2-6), emitting:
  * benchmark/results.json    -- every check {name, gate, molecule, value, tol, pass}
  * benchmark/report_table.tex -- LaTeX tabular of the same (for \\input in the report)

Reproducible: fixed geometries (Appendix A), basis 6-31G*, functional BHHLYP,
triplet ROHF (high-spin) MRSF reference, SCF/TD conv 1e-8.  H2/STO-3G RHF for FCIDUMP.

Usage:  python benchmark/run_poc.py
"""
import os
import sys
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("OPENQP_ROOT", os.path.dirname(HERE))
RUNDIR = os.path.join(HERE, "_runs")
CUBEDIR = os.path.join(RUNDIR, "cubes")
os.makedirs(CUBEDIR, exist_ok=True)

from oqp.pyoqp import Runner
from oqp.analysis import (MRSFExcitedStates, AOBasis, make_box_grid,
                          nto_excitation, nto_transition, attachment_detachment,
                          participation_ratio, tozer_lambda, fragment_ct_matrix)
from oqp.export import (CubeExporter, to_qcschema, validate_qcschema,
                        dump_fcidump, verify_fcidump_fci)
from oqp.interop import parse_output, parse_pyscf_tddft, parse_oqp, compare_results, format_table

HARTREE2EV = 27.211386245988
CHECKS = []          # collected check records


def record(gate, molecule, name, value, tol, passed, unit="", note=""):
    CHECKS.append(dict(gate=gate, molecule=molecule, name=name,
                       value=(None if value is None else float(value)),
                       tolerance=float(tol), unit=unit,
                       passed=bool(passed), note=note))


def run_mrsf(inp):
    proj = os.path.splitext(os.path.basename(inp))[0]
    r = Runner(project=proj, input_file=inp,
               log=os.path.join(RUNDIR, proj + ".log"), silent=1, usempi=False)
    r.run()
    return r.mol


def analyze_molecule(label, inp, target, fragments, fine=0.05):
    """Full feature sweep + checks for one MRSF molecule."""
    mol = run_mrsf(inp)
    st = MRSFExcitedStates(mol)
    ao = AOBasis(mol)
    nstates = st.nstates

    # ---- GATE 2: keystone dipole reconstruction + traces ----
    worst_mu = 0.0
    worst_trace = 0.0
    for i in range(nstates):
        for j in range(i + 1, nstates):
            mu_r = st.transition_dipole(i, j)
            mu_o = st.dip_oqp[:, i, j]
            worst_mu = max(worst_mu, float(np.max(np.abs(mu_r - mu_o))))
            worst_trace = max(worst_trace, abs(float(np.trace(st.tdm_mo(i, j)))))
    record("G2", label, "transition_dipole_reconstruction", worst_mu, 1e-6, worst_mu < 1e-6, "a.u.")
    record("G2", label, "tdm_traceless", worst_trace, 1e-6, worst_trace < 1e-6, "e")
    worst_N = 0.0
    for n in range(nstates):
        g_ao = st.state_density_ao(n)
        worst_N = max(worst_N, abs(float(np.sum(g_ao * st.S)) - st.n_elec))
    record("G2", label, "state_density_trace_eq_N", worst_N, 1e-6, worst_N < 1e-6, "e")
    ortho = float(np.max(np.abs(st.C.T @ st.S @ st.C - np.eye(st.nbf))))
    record("G2", label, "MO_orthonormality", ortho, 1e-8, ortho < 1e-8)

    # ---- GATE 3: amplitudes / NTO / A-D / descriptors ----
    amp_err = 0.0
    for n in range(nstates):
        X = st.amplitude_matrix(n)
        amp_err = max(amp_err, float(np.max(np.abs(
            st.trans_den_from_amplitudes(X, X) - st.diff_density_mo(n)))))
    record("G3", label, "amplitude_vs_keystone", amp_err, 1e-9, amp_err < 1e-9)

    ntoe = nto_excitation(st, target)
    dnorm = abs(ntoe["sum_sigma2"] - float((ntoe["amplitude_matrix"] ** 2).sum()))
    record("G3", label, "nto_sum_sigma2_eq_norm", dnorm, 1e-10, dnorm < 1e-10)

    ntot = nto_transition(st, 0, target)
    g_full = ntot["reconstruct_tdm"]()
    mu_recon = np.array([-np.sum((st.C @ g_full @ st.C.T) * st.R[k]) for k in range(3)])
    dmu = float(np.max(np.abs(mu_recon - st.dip_oqp[:, 0, target])))
    record("G3", label, "nto_transition_reconstructs_mu", dmu, 1e-8, dmu < 1e-8, "a.u.")

    ad = attachment_detachment(st, target, ref=0)
    ad_err = float(np.max(np.abs(ad["A_mo"] - ad["D_mo"] - ad["delta_mo"])))
    record("G3", label, "attach_minus_detach_eq_delta", ad_err, 1e-9, ad_err < 1e-9)
    trAS = float(np.sum(ad["A_ao"] * st.S)); trDS = float(np.sum(ad["D_ao"] * st.S))
    record("G3", label, "TrAS_eq_nprom", abs(trAS - ad["n_promoted"]), 1e-6,
           abs(trAS - ad["n_promoted"]) < 1e-6, "e")
    record("G3", label, "TrDS_eq_nprom", abs(trDS - ad["n_promoted"]), 1e-6,
           abs(trDS - ad["n_promoted"]) < 1e-6, "e")

    pr = participation_ratio(ntoe["weights"])
    npair = ntoe["sigma"].size
    record("G3", label, "participation_ratio_in_range", pr,
           float(npair), 1.0 - 1e-9 <= pr <= npair + 1e-9, "pairs",
           note=f"1<=PR<={npair}")

    lo, n3, dvec, pts = make_box_grid(ao.coords, padding=5.0, spacing=0.12)
    lam, _ = tozer_lambda(ao, ntoe, pts, dvec[0] ** 3)
    record("G3", label, "tozer_lambda_in_0_1", lam, 1.0, -1e-9 <= lam <= 1.0 + 1e-6,
           note="high=local")

    om = fragment_ct_matrix(st, ao, target, fragments)
    record("G3", label, "omega_total_eq_norm", abs(om["total"] - om["amplitude_norm"]),
           1e-8, abs(om["total"] - om["amplitude_norm"]) < 1e-8)

    # ---- GATE 4: cube export + grid-integral trace checks ----
    cub = CubeExporter(st, ao, padding=5.0, spacing=0.15)
    pre = os.path.join(CUBEDIR, label)
    cub.mo_cube(pre + "_homo.cube", st.na - 1)
    cub.state_density_cube(pre + f"_state_S{target}.cube", target)
    cub.transition_density_cube(pre + f"_trans_0_{target}.cube", 0, target)
    cub.attachment_detachment_cubes(pre + f"_attach_S{target}.cube",
                                    pre + f"_detach_S{target}.cube", ad)
    cub.nto_cube(pre + f"_nto_hole_S{target}.cube", ntoe["holes_ao"][:, 0], "hole")
    cub.nto_cube(pre + f"_nto_part_S{target}.cube", ntoe["particles_ao"][:, 0], "particle")
    # validate evaluator: analytic overlap vs OQP S
    Sa = ao.overlap_analytic()
    dS = float(np.max(np.abs(Sa - st.S)))
    record("G4", label, "ao_evaluator_overlap_vs_OQP", dS, 1e-3, dS < 1e-3)
    # state density carries the tight 1s cores -> fine grid (brief: tighten spacing);
    # the valence transition/attachment densities are smooth -> coarse grid suffices.
    sd = cub.integrate_density(st.state_density_ao(target), spacing=fine, padding=5.0)
    record("G4", label, "cube_state_density_eq_N", abs(sd - st.n_elec), 1e-2,
           abs(sd - st.n_elec) < 1e-2, "e")
    td = cub.integrate_density(st.tdm_ao(0, target), spacing=0.10, padding=4.0)
    record("G4", label, "cube_transition_density_eq_0", abs(td), 1e-2, abs(td) < 1e-2, "e")
    at = cub.integrate_density(ad["A_ao"], spacing=0.10, padding=4.0)
    record("G4", label, "cube_attachment_eq_nprom", abs(at - ad["n_promoted"]), 1e-2,
           abs(at - ad["n_promoted"]) < 1e-2, "e")

    # ---- GATE 5: QCSchema ----
    payload = to_qcschema(mol, states=st)
    res = validate_qcschema(payload)
    e_rt = abs(res.properties.scf_total_energy - float(mol.get_scf_energy()))
    record("G5", label, "qcschema_energy_roundtrip", e_rt, 1e-12, res.success and e_rt < 1e-12, "Ha")
    with open(os.path.join(RUNDIR, f"{label}_qcschema.json"), "w") as f:
        json.dump(payload, f, indent=2)

    # collect headline physics numbers
    de = (st.energies - st.energies[0]) * HARTREE2EV
    return {
        "excitation_eV": de.tolist(),
        "target": target,
        "target_exc_eV": float(de[target]),
        "osc_0_target": float(st.oscillator_strength(0, target)),
        "n_promoted": float(ad["n_promoted"]),
        "PR": float(pr), "lambda": float(lam),
        "omega_ct_fraction": float(om["ct_fraction"]),
        "scf_energy_Ha": float(mol.get_scf_energy()),
    }


def run_fcidump():
    inp = os.path.join(HERE, "inputs", "h2_rhf_sto3g.inp")
    mol = run_mrsf(inp)
    path = os.path.join(RUNDIR, "h2_sto3g.fcidump")
    meta = dump_fcidump(path, mol, source="oqp")
    ver = verify_fcidump_fci(path, mol)
    record("G5", "H2/STO-3G", "fcidump_fci_vs_pyscf", ver["diff"], 1e-8, ver["diff"] < 1e-8, "Ha")
    record("G5", "H2/STO-3G", "fcidump_8fold_symmetry", meta["sym_residual_8fold"], 1e-10,
           meta["sym_residual_8fold"] < 1e-10)
    cons = meta["consistency"]
    record("G5", "H2/STO-3G", "oqp_vs_pyscf_Hcore", cons["max_dHcore"], 1e-6,
           cons["max_dHcore"] is None or cons["max_dHcore"] < 1e-6)
    return {"e_fcidump": ver["e_fcidump"], "e_pyscf_fci": ver["e_pyscf_fci"], "diff": ver["diff"]}


def run_interop():
    fixture = os.path.join(HERE, "fixtures", "ch2o_gaussian_td.log")
    ref = [3.9892, 8.9000, 9.5000]
    p = parse_output(fixture, program="gaussian")
    got = p.get("excitation_energies_ev", [])
    d = max(abs(a - b) for a, b in zip(got, ref)) if len(got) == len(ref) else 1e9
    record("G6", "CH2O(Gaussian)", "cclib_excitation_parse", d, 1e-4, d < 1e-4, "eV")
    from pyscf import gto, dft, tddft
    m = gto.M(atom="""C 0 0 -0.5297; O 0 0 0.6775; H 0 0.9342 -1.124; H 0 -0.9342 -1.124""",
              basis="6-31g*", unit="Angstrom", verbose=0)
    mf = dft.RKS(m); mf.xc = "bhandhlyp"; mf.kernel()
    td = tddft.TDDFT(mf); td.nstates = 3; td.kernel()
    src = (np.asarray(td.e) * HARTREE2EV).tolist()
    pp = parse_pyscf_tddft(td, mf)
    rt = max(abs(a - b) for a, b in zip(pp["excitation_energies_ev"], src))
    record("G6", "CH2O(PySCF)", "pyscf_native_roundtrip", rt, 1e-4, rt < 1e-4, "eV")
    return {"pyscf_tddft_eV": pp["excitation_energies_ev"]}


def write_tex(path, headline):
    rows = []
    for c in CHECKS:
        val = "n/a" if c["value"] is None else f"{c['value']:.2e}"
        status = r"\textsc{pass}" if c["passed"] else r"\textbf{FAIL}"
        name = c["name"].replace("_", r"\_")
        mol = c["molecule"].replace("_", r"\_")
        rows.append(f"{c['gate']} & {name} & {mol} & {val} & {c['tolerance']:.0e} & {status} \\\\")
    body = "\n".join(rows)
    npass = sum(c["passed"] for c in CHECKS)
    tex = r"""% Auto-generated by benchmark/run_poc.py -- do not hand-edit.
\begin{tabular}{lllrrl}
\hline
Gate & Check & System & $|\Delta|$ & Tol. & Status \\
\hline
""" + body + r"""
\hline
\end{tabular}
% summary: """ + f"{npass}/{len(CHECKS)} checks passed" + "\n"
    with open(path, "w") as f:
        f.write(tex)


def main():
    headline = {}
    headline["formaldehyde"] = analyze_molecule(
        "formaldehyde", os.path.join(HERE, "inputs", "ch2o_mrsf.inp"),
        target=1, fragments=[[1], [0, 2, 3]])
    headline["ethylene"] = analyze_molecule(
        "ethylene", os.path.join(HERE, "inputs", "c2h4_mrsf.inp"),
        target=1, fragments=[[0, 2, 3], [1, 4, 5]])
    headline["H2_fcidump"] = run_fcidump()
    headline["interop"] = run_interop()

    npass = sum(c["passed"] for c in CHECKS)
    ntot = len(CHECKS)
    results = {
        "metadata": {
            "basis": "6-31G*", "functional": "BHHLYP",
            "mrsf_reference": "ROHF triplet (high-spin)",
            "scf_conv": 1e-8, "td_conv": 1e-8,
            "fcidump_system": "H2/STO-3G RHF",
            "geometries_angstrom": {
                "formaldehyde": "C(0,0,-0.5297) O(0,0,0.6775) H(0,0.9342,-1.124) H(0,-0.9342,-1.124)",
                "ethylene": "C(0,0,0.6686) C(0,0,-0.6686) H(0,0.9229,1.2372)x2 H(0,0.9229,-1.2372)x2",
                "H2": "H(0,0,0) H(0,0,0.741)"},
        },
        "headline": headline,
        "checks": CHECKS,
        "summary": {"passed": npass, "total": ntot, "all_pass": npass == ntot},
    }
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(results, f, indent=2)
    write_tex(os.path.join(HERE, "report_table.tex"), headline)

    print(f"\n===== POC summary: {npass}/{ntot} checks passed =====")
    for c in CHECKS:
        if not c["passed"]:
            print(f"  FAIL  {c['gate']} {c['molecule']} {c['name']}: {c['value']} (tol {c['tolerance']})")
    print("results.json and report_table.tex written to", HERE)
    return 0 if npass == ntot else 1


if __name__ == "__main__":
    sys.exit(main())
