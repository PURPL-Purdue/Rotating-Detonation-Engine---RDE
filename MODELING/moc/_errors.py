"""Exceptions raised by the MoC solver (MATLAB ``error(...)`` equivalents)."""


class MocConvergenceError(RuntimeError):
    """A predictor-corrector or point-location iteration did not converge."""
