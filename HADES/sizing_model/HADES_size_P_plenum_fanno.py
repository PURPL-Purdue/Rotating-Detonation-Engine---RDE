"""Port of ``HADES_size_P_plenum_fanno.m`` (H2 plenum pressure via Fanno flow)."""

from __future__ import annotations

import math
import warnings

import numpy as np
from scipy.optimize import brentq, fsolve

from rdelib.gasdynamics import colebrook

__all__ = ["HADES_size_P_plenum_fanno"]


def HADES_size_P_plenum_fanno(Pc, m_dot, A_inj, L, D, T, Cd, gamma, R, mu_ref, T_ref, S_ref,
                              epsilon, M_G):
    """Required plenum stagnation pressure P01 [Pa] for an H2 injector orifice.

    Works backwards: chamber pressure -> Fanno flow in the orifice tube -> Cd
    contraction -> plenum pressure.

    Parameters
    ----------
    Pc : chamber pressure [Pa]
    m_dot : H2 mass flow rate per orifice [kg/s]
    A_inj : orifice area [m^2]
    L : orifice tube length [mm]
    D : orifice tube diameter [mm]
    T : static inlet temperature [K]
    Cd : discharge coefficient [-]
    gamma : specific heat ratio [-]
    R : gas constant [J/(kg K)]
    mu_ref, T_ref, S_ref : Sutherland reference viscosity [Pa s], temperature [K], constant [K]
    epsilon : surface roughness [mm]
    M_G : initial (subsonic) Mach number guess
    """
    tol = 1e-6
    max_iter = 100
    err = 1.0
    it = 0

    P2 = Pc * 1.5  # initial downstream pressure guess [Pa]
    T_stag = T  # initial stagnation temperature guess [K]
    f = M_inlet = math.nan

    while err > tol and it < max_iter:
        it += 1
        # Sutherland viscosity (as written in the MATLAB)
        mu = mu_ref * (T_ref + S_ref) / (T_stag + S_ref) * (T / T_ref) ** (3 / 2)
        # Reynolds number
        Re = (P2 * M_G * math.sqrt(gamma * R * T) * (D / 1000)) / (R * T * mu)
        # Colebrook (Darcy) friction factor, fsolve from 0.02
        f = colebrook(Re, epsilon / D, 0.02)

        # Fanno: solve fL/D(M) = f L / D for the inlet Mach number (fsolve from M_G)
        fLD = f * L / D

        def fanno_fun(M):
            M = M[0] if np.ndim(M) else M
            return (1 - M**2) / (gamma * M**2) + ((gamma + 1) / (2 * gamma)) * math.log(
                ((gamma + 1) * M**2) / (2 + (gamma - 1) * M**2)
            ) - fLD

        M_inlet = float(fsolve(fanno_fun, M_G)[0])
        M_G = M_inlet
        if M_inlet <= 0:
            raise RuntimeError("HADES: Fanno solver returned nonphysical Mach number.")

        T_stag_new = T * (1 + ((gamma - 1) / 2) * M_inlet**2)
        P2_new = Pc * (1 / M_inlet) * (1 / math.sqrt((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M_inlet**2)))

        err = abs((P2_new - P2) / P2)
        P2 = P2_new
        T_stag = T_stag_new

    if it == max_iter:
        warnings.warn("HADES: Fanno iteration did not fully converge.")

    # Cd contraction -- subcritical compressible orifice equation solved for P01
    def P01_fun(P0):
        return Cd * A_inj * P0 * math.sqrt(
            (2 * gamma / (R * T * (gamma - 1))) * (P2 / P0) ** (2 / gamma) * (1 - (P2 / P0) ** ((gamma - 1) / gamma))
        ) - m_dot

    a = P2 * 1.01
    b = a + 1e8  # 1000 bar search window
    if P01_fun(a) * P01_fun(b) < 0:
        P01 = brentq(P01_fun, a, b, xtol=1e-12, rtol=4 * np.finfo(float).eps)
    else:
        raise RuntimeError("HADES: P01 solver could not find a root in search window. Check inputs.")

    print("")
    print("HADES P_PLENUM FANNO - H2 RESULTS")
    print("Darcy friction factor: %.4f" % f)
    print("Orifice Inlet Mach Number : %.4f" % M_inlet)
    print("Orifice Inlet Static Pressure : %.4f bar" % (P2 / 1e5))
    print("Required Plenum Pressure (P01) : %.4f bar" % (P01 / 1e5))
    print("Required Plenum Pressure (P01) : %.4f psi" % (P01 * 14.504 / 1e5))
    print("Pressure Drop (P01 - Pc) : %.4f bar" % ((P01 - Pc) / 1e5))
    return P01
