"""Input preflight for analytic MRSF Hessians.

The preflight must agree with the native MRSF Hessian gates: accept the LibXC
spellings OpenQP maps to a verified functional, and reject range separation
before any SCF or response work is done.
"""

import unittest

from oqp.utils.input_checker import analytic_hessian_capability


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

    def test_truncated_sos_raman_needs_its_tail_roots(self):
        request = _mrsf_request("bhhlyp")
        request["input"]["runtype"] = "hess"
        request["tdhf"]["nstate"] = 4
        request["hess"]["state"] = 1
        status, reason = analytic_hessian_capability(request)
        self.assertEqual(status, "unsupported_feature")
        self.assertIn("raman_sos_tail_states", reason)
        request["hess"]["raman_backend"] = "finite_field"
        self.assertEqual(analytic_hessian_capability(request)[0], "supported")

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
