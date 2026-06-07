#!/usr/bin/env python3
"""Probe the pyddx API with the upstream point-charge example.

This is intentionally standalone and disposable. It answers whether a local
Python environment can import pyddx and run the basic host-code workflow:
model -> solute electrostatics/psi -> state solve -> energy/forces.
"""

from __future__ import annotations

import importlib
import sys


def main() -> int:
    try:
        import numpy as np
        pyddx = importlib.import_module("pyddx")
    except Exception as exc:  # pragma: no cover - diagnostic script
        print("IMPORT_FAILED")
        print(f"python={sys.executable}")
        print(f"error={type(exc).__name__}: {exc}")
        return 2

    tobohr = 1 / 0.52917721092
    charges = np.array([
        -0.04192, -0.04192, -0.04198, -0.04192, -0.04192, -0.04198,
        0.04193, 0.04193, 0.04197, 0.04193, 0.04193, 0.04197,
    ])
    radii = tobohr * np.array([
        4.00253, 4.00253, 4.00253, 4.00253, 4.00253, 4.00253,
        2.99956, 2.99956, 2.99956, 2.99956, 2.99956, 2.99956,
    ])
    centres = tobohr * np.array([
        [0.00000, 2.29035, 1.32281],
        [0.00000, 2.29035, -1.32281],
        [0.00000, 0.00000, -2.64562],
        [0.00000, -2.29035, -1.32281],
        [0.00000, -2.29035, 1.32281],
        [0.00000, 0.00000, 2.64562],
        [0.00103, 4.05914, 2.34326],
        [0.00103, 4.05914, -2.34326],
        [0.00000, 0.00000, -4.68652],
        [-0.00103, -4.05914, -2.34326],
        [-0.00103, -4.05914, 2.34326],
        [0.00000, 0.00000, 4.68652],
    ]).T

    model = pyddx.Model("pcm", centres, radii, solvent_epsilon=78.3553)
    solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
    solute_field = model.multipole_electrostatics(solute_multipoles)
    solute_psi = model.multipole_psi(solute_multipoles)

    state = pyddx.State(model, solute_psi, solute_field["phi"])
    state.fill_guess()
    state.solve()
    state.fill_guess_adjoint()
    state.solve_adjoint()

    energy = 0.5 * np.sum(state.x * solute_psi)
    force = state.solvation_force_terms(solute_field)
    force += state.multipole_force_terms(solute_multipoles)

    print("OK")
    print(f"energy={energy:.16f}")
    print(f"force_shape={force.shape}")
    print(f"force_max_abs={np.max(np.abs(force)):.8e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
