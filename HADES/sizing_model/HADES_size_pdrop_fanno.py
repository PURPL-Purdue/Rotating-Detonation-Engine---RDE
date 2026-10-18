"""Port of ``HADES_size_pdrop_fanno.m`` (choked Fanno flow, Tkachenko friction factor)."""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from rdelib.gasdynamics import fanno_fLD, tkachenko_friction_factor

__all__ = ["HADES_size_pdrop_fanno"]


def HADES_size_pdrop_fanno(P0, D, L, epsilon, m_dot, gamma, mu):
    """Inlet/exit pressures of a choked Fanno tube.

    Parameters
    ----------
    P0 : stagnation pressure [Pa]
    D : diameter [m]
    L : length [m]
    epsilon : surface roughness [m]
    m_dot : mass flow rate [kg/s]
    gamma : ratio of specific heats [-]
    mu : viscosity [Pa s]

    Returns
    -------
    (P_inlet [Pa], P_exit [Pa], Ma_inlet [-], f [-, Darcy])
    """
    Re = (4 * m_dot) / (np.pi * D * mu)
    rel_rough = epsilon / D
    f = float(tkachenko_friction_factor(Re, rel_rough))

    fanno_param = (f * L) / D
    f_to_solve = lambda Ma: float(fanno_fLD(Ma, gamma)) - fanno_param  # noqa: E731
    # fzero(f, [0.001, 0.999]) -- MATLAB errors if the bracket has no sign change
    Ma_inlet = brentq(f_to_solve, 0.001, 0.999, xtol=1e-15, rtol=4 * np.finfo(float).eps)

    P_inlet = P0 * (1 + (gamma - 1) / 2 * Ma_inlet**2) ** (-gamma / (gamma - 1))
    P_ratio_star = (1 / Ma_inlet) * np.sqrt((gamma + 1) / (2 + (gamma - 1) * Ma_inlet**2))
    P_exit = P_inlet / P_ratio_star
    return P_inlet, P_exit, Ma_inlet, f
