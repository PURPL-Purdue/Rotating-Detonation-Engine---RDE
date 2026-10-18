"""Detonation triple point: oblique shock in the bounding gas vs expansion of the products.

Python translation of ``RDE_Toolbox_2/src/MoC_Resolve_Triple_Point.m``
(Fievisohn eqn. 8; Sousa eqns. 1-4).  MATLAB ``fsolve`` (trust-region-dogleg,
FunctionTolerance = StepTolerance = 1e-8) -> ``scipy.optimize.fsolve``
(MINPACK hybrd, Powell dogleg) with ``xtol = 1e-8``.
"""

import math
import warnings

import numpy as np
from scipy.optimize import fsolve


def MoC_Resolve_Triple_Point(iTripleParam):
    """Return ``(beta_shock, delta_slip, P_matched, M3, M1_prime)``; angles in degrees.

    Unknowns: M3 (post-detonation Mach after the Prandtl-Meyer expansion) and
    the shock angle theta such that pressure and flow direction match across
    the slip line.
    """
    P1 = iTripleParam.P1
    P2 = iTripleParam.P2
    T2 = iTripleParam.T2
    Vcj = iTripleParam.Vcj
    gamma2 = iTripleParam.gamma2
    R1 = iTripleParam.R1

    # 1. Bounding gas: detonation products isentropically expanded to the reactant pressure
    T1_prime = T2 * (P1 / P2) ** ((gamma2 - 1) / gamma2)

    # Acoustic speed and Mach number of the bounding gas
    a1_prime = math.sqrt(gamma2 * R1 * T1_prime)
    M1_prime = Vcj / a1_prime

    # 2. Nonlinear system for the triple point
    M3_guess = 1.5
    theta_shock_guess = 45 * (math.pi / 180)  # radians
    x0 = np.array([M3_guess, theta_shock_guess])

    sol, info, ier, msg = fsolve(
        lambda x: triple_point_residuals(x, P1, P2, gamma2, gamma2, M1_prime),
        x0, xtol=1e-8, full_output=True,
    )
    if ier != 1:
        warnings.warn("Triple point solver did not converge.", RuntimeWarning, stacklevel=2)

    # 3. Outputs
    M3 = float(sol[0])
    beta_shock = float(sol[1])  # radians

    # Final matched slip-line angle (delta_3) (Sousa eqn. 3)
    term1 = math.sqrt((gamma2 + 1) / (gamma2 - 1))
    term2 = math.sqrt((gamma2 - 1) / (gamma2 + 1) * (M3**2 - 1))
    delta_slip = term1 * math.atan(term2) - math.atan(math.sqrt(M3**2 - 1))

    # Final matched static pressure (P3) (Sousa eqn. 4)
    P_matched = P2 * ((1 + (gamma2 - 1) / 2) / (1 + (gamma2 - 1) / 2 * M3**2)) ** (gamma2 / (gamma2 - 1))

    # Angles in degrees
    beta_shock = beta_shock * (180 / math.pi)
    delta_slip = delta_slip * (180 / math.pi)
    return beta_shock, delta_slip, P_matched, M3, M1_prime


def triple_point_residuals(x, P1, P2, gamma1, gamma2, M1_prime):
    """Residuals [P3 - P2', delta_3 - delta_2'] of the triple-point matching."""
    M3 = x[0]
    theta = x[1]

    # Prevent unphysical values during iteration
    if M3 <= 1 or theta <= 0 or theta >= math.pi / 2:
        return np.array([1e6, 1e6])

    # Right side: Prandtl-Meyer expansion (assuming Mw2 = 1) (Sousa eqn. 3)
    term1_exp = math.sqrt((gamma2 + 1) / (gamma2 - 1))
    term2_exp = math.sqrt((gamma2 - 1) / (gamma2 + 1) * (M3**2 - 1))
    delta_3 = term1_exp * math.atan(term2_exp) - math.atan(math.sqrt(M3**2 - 1))

    # (Sousa eqn. 4)
    P3 = P2 * ((1 + (gamma2 - 1) / 2) / (1 + (gamma2 - 1) / 2 * M3**2)) ** (gamma2 / (gamma2 - 1))

    # Left side: oblique shock (Sousa eqn. 2), pressure ratio across the shock
    P2_prime = P1 * (1 + (2 * gamma1 / (gamma1 + 1)) * (M1_prime**2 * math.sin(theta) ** 2 - 1))

    # Flow deflection across the oblique shock (Sousa eqn. 1)
    num = M1_prime**2 * math.sin(theta) ** 2 - 1
    den = M1_prime**2 * (gamma1 + math.cos(2 * theta)) + 2
    delta_2_prime = math.atan(2 * (1.0 / math.tan(theta)) * (num / den))

    # Residuals to drive to zero
    return np.array([P3 - P2_prime, delta_3 - delta_2_prime])
