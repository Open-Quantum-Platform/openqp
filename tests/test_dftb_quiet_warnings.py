"""A quiet DFTB log still carries native warnings and keeps the energy data.

``[dftb] print_level`` resolves to 0 at ``[input] verbose = 0``.  The trace
used to be discarded before parsing on that path, which dropped the state
spectrum warnings and the energy components the final report reads.
"""
import types
import unittest
from unittest import mock

try:
    from oqp.library import openqp_dftb
    _HAVE = True
except Exception:  # pragma: no cover - runtime unavailable
    _HAVE = False

TRACE = "\n".join([
    "openqp_dftb_energy_components total -1.25 repulsive 0.5",
    "openqp_dftb_state_spectrum_warning pair=1-2 reason=near_degenerate",
    "note: warning from the native backend about the reference",
])


@unittest.skipUnless(_HAVE, "OpenQP-DFTB adapter unavailable")
class QuietDftbWarnings(unittest.TestCase):
    def _adapter(self, verbose):
        adapter = openqp_dftb.OpenQPDFTBAdapter.__new__(openqp_dftb.OpenQPDFTBAdapter)
        adapter.config = {"input": {"verbose": verbose}, "dftb": {"print_level": 1}}
        adapter.mol = types.SimpleNamespace()
        return adapter

    def test_quiet_level_keeps_warnings_and_energy_components(self):
        adapter = self._adapter(0)
        with mock.patch.object(openqp_dftb, "dump_log") as dump:
            adapter._log_native_progress("dftb", 0, False, TRACE)
        self.assertEqual(adapter.mol.dftb_energy_components.get("total"), -1.25)
        self.assertEqual(dump.call_count, 1)
        text = dump.call_args.kwargs["info"]["text"]
        self.assertIn("openqp_dftb_state_spectrum_warning", text)
        self.assertIn("warning from the native backend", text)


if __name__ == "__main__":
    unittest.main()
