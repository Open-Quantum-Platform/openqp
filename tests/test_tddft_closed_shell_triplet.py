"""Closed-shell triplet RPA/TDA (tdhf_energy, [tdhf] multiplicity=3).

For a closed-shell reference the triplet response has no Coulomb term and the
f_aa - f_ab XC kernel.  Checks:

* H2O/6-31G*, HF: the triplet TDA roots equal the alpha->beta spin-flip TDA
  roots of the same closed-shell reference (identical matrix
  eps_a - eps_i - (ij|ab), independent code), and the printed triplet
  oscillator strengths vanish.
* H2O/6-31G*, B3LYP: the singlet routed through the unrestricted XC contraction
  (OQP_TDDFT_XC_UNRESTRICTED=1, dxa = dxb = D) equals the restricted singlet,
  which fixes the normalisation that the triplet (dxa = D, dxb = -D) relies on;
  every triplet lies below the singlet of the same order.

Before the triplet branch existed, [tdhf] multiplicity=3 was rejected by the
input checker, so every case here fails without it.
"""
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

EV = 27.211386245988


def _runtime_available():
    try:
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


H2O = ("   O 0.000000 0.000000 0.117300\n   H 0.000000 0.757200 -0.469200\n"
       "   H 0.000000 -0.757200 -0.469200\n")
# C1 (no symmetry element): unequal OH bonds, out-of-plane displacement, so a
# transposed index or a sign error in the response cannot cancel by symmetry.
H2O_C1 = ("   O 0.010000 -0.020000 0.117300\n   H 0.080000 0.770200 -0.449200\n"
          "   H -0.050000 -0.731200 -0.509200\n")


def deck(geom, basis, functional, td_type, mult, nstate):
    s = f"[input]\nsystem=\n{geom}charge=0\nruntype=energy\nbasis={basis}\nmethod=tdhf\n"
    if functional:
        s += f"functional={functional}\n"
    s += ("\n[guess]\ntype=huckel\n\n[scf]\ntype=rhf\nmultiplicity=1\nmaxit=100\nconv=1e-10\n\n"
          f"[tdhf]\ntype={td_type}\nnstate={nstate}\nmaxit=100\nconv=1e-10\nmultiplicity={mult}\n")
    if functional:
        s += "\n[dftgrid]\nrad_npts=96\nang_npts=302\n"
    return s


def run(text, env_extra=None):
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "job.inp").write_text(text)
        env = dict(os.environ, OMP_NUM_THREADS="2")
        env.pop("OQP_TDDFT_XC_UNRESTRICTED", None)
        env.update(env_extra or {})
        proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", "job.inp"], cwd=tmp, env=env,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600)
        log = (Path(tmp) / "job.log").read_text(errors="ignore") if (Path(tmp) / "job.log").exists() else ""
        summary = log.split("Summary of the TD-DFT calculation")[-1]
        ev = [float(x) for x in re.findall(r"^\s*0 ->\s*\d+\s+(-?[\d.]+)", summary, re.M)]
        return proc.returncode, ev, proc.stdout.decode(errors="ignore")[-1500:] + log[-3000:]


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestClosedShellTripletHF(unittest.TestCase):
    """HF limit: the closed-shell triplet TDA matrix, eps_a - eps_i - (ij|ab),
    is the alpha->beta spin-flip TDA matrix of the same closed-shell reference,
    computed by the independent SF-TDA code (tdhf_sf_energy)."""

    def test_triplet_tda_equals_closed_shell_spin_flip_tda(self):
        for geom in (H2O, H2O_C1):
            with self.subTest(c1=geom is H2O_C1):
                self._compare_with_spin_flip(geom)

    def _compare_with_spin_flip(self, geom):
        rct, et, dt = run(deck(geom, "6-31g*", None, "tda", 3, 4))
        sf = deck(geom, "6-31g*", None, "sf", 3, 4).replace("type=rhf", "type=rohf")
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "job.inp").write_text(sf)
            env = dict(os.environ, OMP_NUM_THREADS="2")
            proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", "job.inp"], cwd=tmp, env=env,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600)
            log = (Path(tmp) / "job.log").read_text(errors="ignore")
        esf = [float(x) for x in re.findall(r"State #\s+\d+\s+Energy =\s+(-?[\d.]+) eV", log)]
        self.assertEqual(rct, 0, dt)
        self.assertEqual(proc.returncode, 0, log[-3000:])
        self.assertEqual(len(et), 4, dt)
        self.assertEqual(len(esf), 4, log[-3000:])
        for a, b in zip(et, esf):
            self.assertAlmostEqual(a, b, delta=2e-6)

    def test_triplet_oscillator_strengths_vanish(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "job.inp").write_text(deck(H2O, "6-31g*", None, "tda", 3, 4))
            env = dict(os.environ, OMP_NUM_THREADS="2")
            subprocess.run([sys.executable, "-m", "oqp.pyoqp", "job.inp"], cwd=tmp, env=env,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600)
            log = (Path(tmp) / "job.log").read_text(errors="ignore")
        summary = log.split("Summary of the TD-DFT calculation")[-1]
        osc = [float(x) for x in re.findall(r"^\s*0 ->\s*\d+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)", summary, re.M)]
        self.assertEqual(len(osc), 4, summary[:1500])
        self.assertTrue(all(abs(f) < 1e-12 for f in osc), osc)


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestClosedShellTripletDFT(unittest.TestCase):
    def test_unrestricted_xc_path_reproduces_the_restricted_singlet(self):
        for td in ("tda", "rpa"):
            with self.subTest(td=td):
                rc0, e0, d0 = run(deck(H2O_C1, "6-31g*", "b3lypv5", td, 1, 4))
                rc1, e1, d1 = run(deck(H2O_C1, "6-31g*", "b3lypv5", td, 1, 4), {"OQP_TDDFT_XC_UNRESTRICTED": "1"})
                self.assertEqual(rc0, 0, d0)
                self.assertEqual(rc1, 0, d1)
                self.assertEqual(len(e0), len(e1))
                for a, b in zip(e0, e1):
                    self.assertAlmostEqual(a, b, delta=1e-6)

    def test_triplets_lie_below_singlets(self):
        for td in ("tda", "rpa"):
            with self.subTest(td=td):
                rcs, es, ds = run(deck(H2O, "6-31g*", "b3lypv5", td, 1, 4))
                rct, et, dt = run(deck(H2O, "6-31g*", "b3lypv5", td, 3, 4))
                self.assertEqual(rcs, 0, ds)
                self.assertEqual(rct, 0, dt)
                self.assertLess(et[0], es[0])
                self.assertNotEqual([round(x, 4) for x in et], [round(x, 4) for x in es], dt)


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestClosedShellTripletAntisymmetricPath(unittest.TestCase):
    """The triplet skips the beta half of the unrestricted XC contraction
    (dxb = -dxa makes it the negated alpha half).  It must equal the full
    two-spin contraction (OQP_TDDFT_TRIPLET_FULL_UXC=1), for an LDA (no
    gradient contractions), a GGA hybrid and a meta-GGA, whose tau response
    takes a separate branch."""

    def test_fast_triplet_xc_equals_full_two_spin_contraction(self):
        for func in ("svwn", "b3lypv5", "m06-2x"):
            for td in ("tda", "rpa"):
                with self.subTest(functional=func, td=td):
                    r0, e0, d0 = run(deck(H2O_C1, "6-31g*", func, td, 3, 4))
                    r1, e1, d1 = run(deck(H2O_C1, "6-31g*", func, td, 3, 4),
                                     {"OQP_TDDFT_TRIPLET_FULL_UXC": "1"})
                    self.assertEqual(r0, 0, d0)
                    self.assertEqual(r1, 0, d1)
                    self.assertEqual(len(e0), 4, d0)
                    self.assertEqual(len(e1), 4, d1)
                    for a, b in zip(e0, e1):
                        self.assertAlmostEqual(a, b, delta=1e-6)


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestClosedShellTripletCAM(unittest.TestCase):
    """Range-separated (CAM) exchange runs two integral passes, and
    int2_run_cam sets each pass's Coulomb scale from its arguments.  With
    [dftgrid] cam_alpha=0.2, cam_beta=0 the CAM route reproduces B3LYP's 20 %
    full-range exchange, so its roots must equal the plain B3LYP roots.  A
    triplet that leaks the Coulomb term through the CAM passes differs by
    ~0.8 eV (RPA, C1 H2O)."""

    CAM_AS_B3LYP = "cam_flag=true\ncam_alpha=0.2\ncam_beta=0.0\ncam_mu=0.33\n"

    def test_cam_route_with_full_range_exchange_equals_plain_hybrid(self):
        for mult in (1, 3):
            for td in ("tda", "rpa"):
                with self.subTest(mult=mult, td=td):
                    plain = deck(H2O_C1, "6-31g*", "b3lypv5", td, mult, 4)
                    r0, e0, d0 = run(plain)
                    r1, e1, d1 = run(plain + self.CAM_AS_B3LYP)
                    self.assertEqual(r0, 0, d0)
                    self.assertEqual(r1, 0, d1)
                    self.assertEqual(len(e0), 4, d0)
                    self.assertEqual(len(e1), 4, d1)
                    for a, b in zip(e0, e1):
                        self.assertAlmostEqual(a, b, delta=1e-6)

    def test_cam_b3lyp_triplets_lie_below_singlets(self):
        for td in ("tda", "rpa"):
            with self.subTest(td=td):
                rcs, es, ds = run(deck(H2O_C1, "6-31g*", "cam-b3lyp", td, 1, 4))
                rct, et, dt = run(deck(H2O_C1, "6-31g*", "cam-b3lyp", td, 3, 4))
                self.assertEqual(rcs, 0, ds)
                self.assertEqual(rct, 0, dt)
                self.assertEqual(len(et), 4, dt)
                self.assertLess(et[0], es[0] - 0.1)


C2H4 = ("   C 0.000000 0.000000 0.669500\n   C 0.000000 0.000000 -0.669500\n"
        "   H 0.000000 0.928900 1.232100\n   H 0.000000 -0.928900 1.232100\n"
        "   H 0.000000 0.928900 -1.232100\n   H 0.000000 -0.928900 -1.232100\n")


def _full_log(text):
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "job.inp").write_text(text)
        env = dict(os.environ, OMP_NUM_THREADS="2")
        proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", "job.inp"], cwd=tmp, env=env,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600)
        log = (Path(tmp) / "job.log").read_text(errors="ignore") if (Path(tmp) / "job.log").exists() else ""
        return proc.returncode, log


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestImaginaryRpaRoots(unittest.TestCase):
    """RHF ethylene is triplet-unstable: the lowest HF RPA triplet has
    omega^2 < 0.  rpaechk keeps |omega|, so the root must be reported as
    imaginary, not as a 0.58 eV excitation (NWChem 7.0.2 aborts on it,
    ORCA 5.0.3 prints a different positive number).  TDA (real eigenvalue
    3.55 eV) and a stable reference must stay unmarked."""

    def _deck(self, geom, td_type):
        return deck(geom, "cc-pvdz", None, td_type, 3, 3).replace(
            "method=tdhf\n", "method=tdhf\nispher=true\n")

    def test_unstable_rpa_triplet_root_is_marked_imaginary(self):
        rc, log = _full_log(self._deck(C2H4, "rpa"))
        self.assertEqual(rc, 0, log[-3000:])
        summary = log.split("Summary of the TD-DFT calculation")[-1]
        rows = re.findall(r"^\s*0 ->\s*(\d+)\s+(\S+)", summary, re.M)
        self.assertEqual(len(rows), 3, summary[:1500])
        # the imaginary root is not converged, so no magnitude is reported
        self.assertEqual(rows[0][1], "imaginary", summary[:1500])
        self.assertTrue(all(float(r[1]) > 8.0 for r in rows[1:]), summary[:1500])
        self.assertIn("root(s) imaginary (reference unstable)", summary)
        state1 = re.search(r"PyOQP state 1\s.*", log).group(0)
        state2 = re.search(r"PyOQP state 2\s.*", log).group(0)
        self.assertIn("imaginary", state1)
        self.assertNotRegex(state1, r"-\d+\.\d+")
        self.assertNotIn("imaginary", state2)

    def test_all_requested_roots_imaginary_is_a_breakdown(self):
        # nstate=1: the only requested root has omega^2 < 0, so nothing can
        # be converged; the run must fail instead of reporting it.
        rc, log = _full_log(self._deck(C2H4, "rpa").replace("nstate=3", "nstate=1"))
        self.assertNotEqual(rc, 0, log[-3000:])
        self.assertIn("all 1 requested RPA roots have omega^2 < 0", log)
        self.assertNotIn("TD-DFT energies converged", log)

    def test_tda_and_stable_reference_are_not_marked(self):
        for geom, td_type in ((C2H4, "tda"), (H2O, "rpa")):
            with self.subTest(td=td_type, c2h4=geom is C2H4):
                rc, log = _full_log(self._deck(geom, td_type))
                self.assertEqual(rc, 0, log[-3000:])
                self.assertNotIn("imaginary", log)


H2O_STRETCHED = ("   O 0.000000 0.000000 0.230000\n   H 0.000000 1.550000 -0.920000\n"
                 "   H 0.000000 -1.550000 -0.920000\n")


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestNegativeTdaTripletRoot(unittest.TestCase):
    """H2O with both O-H bonds stretched to ~1.8 A: the RHF reference is
    triplet-unstable and the lowest HF triplet TDA roots are negative.  A
    negative TDA eigenvalue is a real
    root: it must be converged and agree with the independent closed-shell
    SF-TDA code (same matrix eps_a - eps_i - (ij|ab)), also when it is the
    only requested root.  They used to be excluded from the convergence
    test, so nstate=1 reported the first Ritz value (-1.3636 eV instead of
    -1.4793 eV) as converged."""

    def test_negative_triplet_tda_root_equals_spin_flip_tda(self):
        for nstate in (1, 3):
            with self.subTest(nstate=nstate):
                rct, et, dt = run(deck(H2O_STRETCHED, "cc-pvdz", None, "tda", 3, nstate))
                sf = deck(H2O_STRETCHED, "cc-pvdz", None, "sf", 3, nstate).replace("type=rhf", "type=rohf")
                with tempfile.TemporaryDirectory() as tmp:
                    (Path(tmp) / "job.inp").write_text(sf)
                    env = dict(os.environ, OMP_NUM_THREADS="2")
                    proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", "job.inp"], cwd=tmp, env=env,
                                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600)
                    log = (Path(tmp) / "job.log").read_text(errors="ignore")
                esf = [float(x) for x in re.findall(r"State #\s+\d+\s+Energy =\s+(-?[\d.]+) eV", log)]
                self.assertEqual(rct, 0, dt)
                self.assertEqual(proc.returncode, 0, log[-3000:])
                self.assertEqual(len(et), nstate, dt)
                self.assertLess(et[0], 0.0, dt)
                for a, b in zip(et, esf[:nstate]):
                    self.assertAlmostEqual(a, b, delta=2e-6)


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestClosedShellTripletThreads(unittest.TestCase):
    def test_triplet_energies_independent_of_thread_count(self):
        for td in ("tda", "rpa"):
            with self.subTest(td=td):
                r1, e1, d1 = run(deck(H2O_C1, "6-31g*", "b3lypv5", td, 3, 4), {"OMP_NUM_THREADS": "1"})
                r4, e4, d4 = run(deck(H2O_C1, "6-31g*", "b3lypv5", td, 3, 4), {"OMP_NUM_THREADS": "4"})
                self.assertEqual(r1, 0, d1)
                self.assertEqual(r4, 0, d4)
                self.assertEqual(len(e1), 4, d1)
                for a, b in zip(e1, e4):
                    self.assertAlmostEqual(a, b, delta=1e-6)


if __name__ == "__main__":
    unittest.main()
