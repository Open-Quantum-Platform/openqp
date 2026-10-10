"""Compatibility imports for the former PT2-only numerical-gradient driver."""

from oqp.library.wf_numgrad import (  # noqa: F401
    DEFAULT_GAP_WARN,
    DEFAULT_GRAD_GUESS,
    DEFAULT_GRAD_STEP_BOHR,
    PT2_NUMGRAD_METHODS,
    wavefunction_numerical_gradient,
)


pt2_numerical_gradient = wavefunction_numerical_gradient
