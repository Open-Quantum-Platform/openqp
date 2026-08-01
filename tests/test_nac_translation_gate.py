"""Regressions for the raw analytic-NAC translation structural gate."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from tools.nac_lagrangian.translation_gate import (
    TranslationGateError,
    analyze_file,
    analyze_payload,
)


def synthetic_debug_payload() -> dict[str, np.ndarray]:
    nstate = 3
    natom = 2
    payload: dict[str, np.ndarray] = {
        "Xf": np.zeros((5, nstate)),
        "energies": np.array([-10.0, -9.8, -9.6, -9.5]),
    }
    components = ("t1", "z_hf", "z_xc", "vmask")
    for istate in range(nstate):
        for jstate in range(istate + 1, nstate):
            scale = float(1 + istate + 2 * jstate)
            nonoverlap: dict[str, np.ndarray] = {}
            for index, component in enumerate(components, start=1):
                atom = scale * index * np.array([0.01, -0.02, 0.03])
                values = np.vstack((atom, -atom))
                nonoverlap[component] = values
                payload[f"{component}_{istate}{jstate}"] = values
                payload[f"{component}_{jstate}{istate}"] = -values

            # Deliberately retain a large nonzero raw translation component.
            # Passing this payload proves the gate does not impose an ETF/zero
            # sum on a raw electronic derivative coupling.
            gamma = np.array(
                [[0.20 * scale, -0.10 * scale, 0.05 * scale],
                 [0.03 * scale, 0.04 * scale, -0.02 * scale]]
            )
            total = gamma + sum(nonoverlap.values(), np.zeros((natom, 3)))
            payload[f"gsk_{istate}{jstate}"] = gamma
            payload[f"gsk_{jstate}{istate}"] = -gamma
            payload[f"dp_{istate}{jstate}"] = total
            payload[f"dp_{jstate}{istate}"] = -total
    return payload


class NACTranslationGateTests(unittest.TestCase):
    def test_nonzero_raw_atom_sum_passes_structural_gate(self):
        result = analyze_payload(synthetic_debug_payload(), term_atol=1.0e-14)
        self.assertGreater(result.max_raw_atom_sum, 0.1)
        self.assertLessEqual(result.max_nonoverlap_atom_sum, 1.0e-14)
        self.assertLessEqual(result.max_identity_error, 1.0e-14)

    def test_nonoverlap_translation_leak_is_rejected(self):
        payload = synthetic_debug_payload()
        payload["t1_01"] = np.array(payload["t1_01"], copy=True)
        payload["t1_01"][0, 0] += 1.0e-4
        with self.assertRaisesRegex(TranslationGateError, "t1 atom sum"):
            analyze_payload(payload, term_atol=1.0e-8)

    def test_dp_gamma_sk_identity_failure_is_rejected(self):
        payload = synthetic_debug_payload()
        payload["dp_12"] = np.array(payload["dp_12"], copy=True)
        payload["dp_12"][1, 2] += 2.0e-4
        with self.assertRaisesRegex(TranslationGateError, "differs from gamma:Sk"):
            analyze_payload(payload, identity_atol=1.0e-8)

    def test_npz_loader_uses_pickle_free_payload(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "analytic_debug.npz"
            np.savez(path, **synthetic_debug_payload())
            result = analyze_file(path, term_atol=1.0e-14)
        self.assertEqual(len(result.pair_metrics), 3)


if __name__ == "__main__":
    unittest.main()
