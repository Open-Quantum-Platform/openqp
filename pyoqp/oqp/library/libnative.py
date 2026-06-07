"""Native geometry optimizer backend for OpenQP.

A NumPy/SciPy-only alternative to the external ``geomeTRIC`` driver.  The
electronic-structure work, convergence test and logging are reused verbatim from
:class:`oqp.library.libscipy.StateSpecificOpt`; this module only supplies the
step-determination loop through :class:`oqp.library.native_engine.NativeEngine`
(redundant internal coordinates + restricted-step RFO / P-RFO with model-Hessian
BFGS/Bofill updates).

Selected with ``[optimize] lib=native``.  Supported runtypes:

* ``optimize`` -> :class:`NativeOpt`   (state-specific minimum)
* ``ts``       -> :class:`NativeTSOpt` (transition state, eigenvector following)

Options are read from the ``[native]`` input section (see
``oqp.molecule.oqpdata``).
"""

from __future__ import annotations

import numpy as np

from oqp.library.libscipy import StateSpecificOpt, MECIOpt, MECPOpt
from oqp.library.native_engine import NativeEngine
from oqp.utils.file_utils import dump_log


class _NativeRunner:
    """Mixin that drives a PyOQP ``one_step`` objective with NativeEngine."""

    mode = "min"

    def _native_config(self):
        cfg = self.mol.config.get("native", {})
        return {
            "coordsys": cfg.get("coordsys", "auto"),
            "trust": float(cfg.get("trust", 0.2)),
            "trust_max": float(cfg.get("trust_max", 0.5)),
            "follow_mode": int(cfg.get("follow", 0)),
        }

    def optimize(self):
        opts = self._native_config()
        atoms = np.asarray(self.mol.get_atoms(), dtype=int).reshape(-1)
        x0 = np.asarray(self.pre_coord, dtype=float).reshape(-1)

        engine = NativeEngine(
            atoms, x0,
            mode=self.mode,
            trust=opts["trust"],
            trust_max=opts["trust_max"],
            follow_mode=opts["follow_mode"],
            coordsys=opts["coordsys"],
            # NativeEngine's own loop is a backstop; the OQP convergence test
            # (which raises StopIteration at maxit) governs termination.
            maxiter=self.maxit,
        )
        dump_log(
            self.mol,
            title="PyOQP: Native optimizer [mode=%s, coordsys=%s, trust=%.3f]"
            % (self.mode, engine.coordsys, opts["trust"]),
        )

        def energy_gradient(coords):
            # one_step updates metrics/pre_coord/pre_energy and returns
            # (energy [Hartree], gradient [Hartree/Bohr]) for the active state.
            return self.one_step(np.asarray(coords, dtype=float).reshape(-1))

        try:
            engine.run(energy_gradient, on_converged=self.check_convergence)
        except StopIteration:
            pass


class NativeOpt(_NativeRunner, StateSpecificOpt):
    """State-specific minimum via redundant internals + restricted-step RFO."""

    mode = "min"

    def __init__(self, mol):
        super().__init__(mol)


class NativeTSOpt(_NativeRunner, StateSpecificOpt):
    """Transition-state search via partitioned RFO (eigenvector following)."""

    mode = "ts"

    def __init__(self, mol):
        super().__init__(mol)


class NativeMECIOpt(_NativeRunner, MECIOpt):
    """MECI search: minimize the penalty/UBP objective with the native engine.

    Reuses ``MECIOpt.one_step`` (which returns the penalty objective and its
    gradient) and ``MECIOpt.check_convergence`` (which adds the energy-gap
    criterion) verbatim -- only the step determination is native.
    """

    mode = "min"

    def __init__(self, mol):
        MECIOpt.__init__(self, mol)


class NativeMECPOpt(_NativeRunner, MECPOpt):
    """MECP search: minimize the gap-penalty objective with the native engine."""

    mode = "min"

    def __init__(self, mol):
        MECPOpt.__init__(self, mol)
