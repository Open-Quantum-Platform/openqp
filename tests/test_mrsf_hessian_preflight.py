"""Input preflight for analytic MRSF Hessians.

The preflight must agree with the native MRSF Hessian gates: accept the LibXC
spellings OpenQP maps to a verified functional, and reject range separation
before any SCF or response work is done.  For runtype=hess it must also
reject property requests that cannot run at all, and only warn about ones
that complete partially.
"""

import unittest

from oqp.utils.input_checker import analytic_hessian_capability, check_input_values


def _mrsf_request(functional, **dftgrid):
    config = {
        "input": {"method": "tdhf", "functional": functional},
        "scf": {"type": "rohf", "multiplicity": 3},
        "tdhf": {"type": "mrsf", "multiplicity": 1, "nstate": 6},
        "hess": {"state": 3},
    }
    if dftgrid:
        config["dftgrid"] = dict(dftgrid)
    return config



def _checked_request(nstate, **hess):
    # A complete runtype=hess request that passes the whole input checker.
    return {
        "input": {"method": "tdhf", "runtype": "hess", "basis": "sto-3g",
                  "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3"},
        "scf": {"type": "rohf", "multiplicity": 3},
        "tdhf": {"type": "mrsf", "nstate": nstate, "multiplicity": 3},
        "hess": dict({"type": "analytical", "state": 1, "nproc": 1,
                      "temperature": [298.15]}, **hess),
    }

class MrsfHessianPreflight(unittest.TestCase):
    def test_libxc_aliases_of_verified_functionals_are_accepted(self):
        for alias, canonical in (("b3lypv5", "b3lyp5"), ("pbepbe", "pbe")):
            with self.subTest(alias=alias):
                self.assertEqual(
                    analytic_hessian_capability(_mrsf_request(canonical))[0], "supported")
                self.assertEqual(
                    analytic_hessian_capability(_mrsf_request(alias))[0], "supported")

    def test_unverified_functionals_stay_rejected(self):
        status, _ = analytic_hessian_capability(_mrsf_request("tpss"))
        self.assertEqual(status, "unsupported_feature")

    def test_range_separation_is_rejected_before_the_calculation(self):
        status, reason = analytic_hessian_capability(_mrsf_request("pbe", cam_flag="true"))
        self.assertEqual(status, "unsupported_feature")
        self.assertIn("cam_flag", reason)
        self.assertEqual(
            analytic_hessian_capability(_mrsf_request("pbe", cam_flag="false"))[0],
            "supported")

    def test_intensities_need_a_solved_root_above_the_target(self):
        request = _mrsf_request("bhhlyp")
        request["input"]["runtype"] = "hess"
        request["hess"]["state"] = request["tdhf"]["nstate"]
        status, reason = analytic_hessian_capability(request)
        self.assertEqual(status, "unsupported_feature")
        self.assertIn("tdhf.nstate", reason)
        request["hess"]["vibrational_intensities"] = False
        self.assertEqual(analytic_hessian_capability(request)[0], "supported")

    def test_short_sos_tail_only_warns_because_ir_still_runs(self):
        # truncated_sos_polarizability's ValueError is caught per displacement:
        # IR intensities are published and Raman is marked unavailable.
        request = _checked_request(nstate=4)
        self.assertEqual(analytic_hessian_capability(request)[0], "supported")
        report = check_input_values(request, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())
        self.assertIn("raman_sos_tail_states", report.to_text())

    def test_finite_field_raman_is_rejected_because_it_cannot_run(self):
        request = _mrsf_request("bhhlyp")
        request["input"]["runtype"] = "hess"
        request["hess"]["raman_backend"] = "finite_field"
        status, reason = analytic_hessian_capability(request)
        self.assertEqual(status, "unsupported_feature")
        self.assertIn("finite_field", reason)
        request["hess"]["vibrational_intensities"] = False
        self.assertEqual(analytic_hessian_capability(request)[0], "supported")

    def test_cached_hessian_reads_skip_the_property_requirements(self):
        request = _mrsf_request("bhhlyp")
        request["input"]["runtype"] = "hess"
        request["hess"]["state"] = request["tdhf"]["nstate"]
        request["hess"]["read"] = True
        self.assertEqual(analytic_hessian_capability(request)[0], "supported")

    def test_numerical_mrsf_hessians_get_the_same_root_check(self):
        request = _checked_request(nstate=1, type="numerical")
        report = check_input_values(request, raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertIn("tdhf.nstate", report.to_text())
        request["tdhf"]["nstate"] = 6
        report = check_input_values(request, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())

    def test_bare_b3lyp_is_rejected_as_ambiguous(self):
        # source/dftlib/libxc.F90 aborts on bare B3LYP and asks for a variant.
        status, reason = analytic_hessian_capability(_mrsf_request("b3lyp"))
        self.assertEqual(status, "unsupported_feature")
        self.assertIn("ambiguous", reason)
        self.assertEqual(analytic_hessian_capability(_mrsf_request("b3lyp5"))[0], "supported")

    def test_numerical_root_check_does_not_depend_on_intensities(self):
        # numerical_hess tracks the target root at every displacement anyway.
        request = _checked_request(nstate=1, type="numerical", vibrational_intensities=False)
        report = check_input_values(request, raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertIn("tdhf.nstate", report.to_text())

    def test_invalid_mrsf_property_options_fail_preflight(self):
        # MRSFPropertyFDRequest.create rejects these, but only after the Hessian.
        for options in ({"raman_backend": "truncated-SOS"}, {"raman_backend": "TRUNCATED_SOS"},
                        {"raman_sos_tail_states": 0}, {"property_dx": -1.0e-3},
                        {"property_min_overlap": 0.5}):
            with self.subTest(**options):
                request = _mrsf_request("bhhlyp")
                request["input"]["runtype"] = "hess"
                request["hess"].update(options)
                status, reason = analytic_hessian_capability(request)
                self.assertEqual(status, "unsupported_feature")
                self.assertIn("property options are invalid", reason)
                request["hess"]["vibrational_intensities"] = False
                self.assertEqual(analytic_hessian_capability(request)[0], "supported")

    def test_rejected_requests_get_no_raman_tail_warning(self):
        report = check_input_values(_checked_request(nstate=1), raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertNotIn("raman_sos_tail_states", report.to_text())

    def test_shipped_analytic_mrsf_decks_still_pass(self):
        # examples/HESS/*MRSF_ANALYTIC_HESSIAN*.inp: runtype=hess, state=3, nstate=6.
        for functional in ("bhhlyp", ""):
            request = _mrsf_request(functional)
            request["input"]["runtype"] = "hess"
            self.assertEqual(analytic_hessian_capability(request)[0], "supported")

    def test_ts_and_irc_initial_hessians_skip_the_property_requirements(self):
        request = _mrsf_request("bhhlyp")
        request["input"]["runtype"] = "ts"
        request["tdhf"]["nstate"] = 3
        request["hess"]["state"] = 1
        self.assertEqual(analytic_hessian_capability(request)[0], "supported")


if __name__ == "__main__":
    unittest.main()
