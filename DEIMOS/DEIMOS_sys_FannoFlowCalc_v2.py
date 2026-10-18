"""The big FFC (Fanno Flow Calculator).

Canonical Python port of ``DEIMOS_sys_FannoFlowCalc_v2.m``; supersedes
``DEIMOS_sys_FannoFlowCalc_v1.m`` (see ``PORTING_NOTES.md``).

NOTES (from the MATLAB header)
    - Oxygen is first column, methane second
    - Don't guess 1 for M_G, calc needs either subsonic or supersonic guess
    - All other guesses used only for Darcy friction calculations,
      iterative convergence with Fanno fixes output values
    - P01 Solver Diagnostic is toggled with ``main(diagnostic=...)``
      (the MATLAB toggled it with a ``%{`` block comment; it was ON)
    - The MATLAB needed the Optimization Toolbox (fsolve). Here the scalar roots
      are found with bracketed ``scipy.optimize.brentq`` (same roots; the
      Colebrook root is unique and the Fanno branch is picked by ``M_G``).

Model, per propellant (oxygen, methane), iterated on the tube-inlet static
pressure ``P2`` until ``max|dP2/P2| <= tol``:

1. Sutherland viscosity at the current stagnation temperature ``T_stag``.
2. Reynolds number ``P2 * M_G * sqrt(gamma R T) * D / (R T mu)`` (uses the
   Mach *guess* ``M_G``, exactly as the MATLAB does).
3. Darcy friction factor from Colebrook-White.
4. Tube inlet Mach ``M_TubeIN`` from the Fanno relation f L*/D = f L/D
   (tube choked at its exit into the chamber).
5. Upstream (manifold) stagnation pressure ``P01`` from the subcritical
   isentropic orifice equation with downstream pressure ``P2``.
6. ``T_stag = T (1 + (gamma-1)/2 M^2)`` and ``P2 = Pc * (P/P*)_Fanno(M)``.

Units: inj_OD, D, L, epsilon, Rmid in mm; pressures in Pa; temperatures in K.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.optimize import brentq

from rdelib import gasdynamics as _gd
from rdelib.gasdynamics import (
    colebrook_darcy,
    critical_pressure_ratio,
    fanno_fLD_max,
    subcritical_orifice_mass_flow,
)
from rdelib.gasdynamics import fanno_T_over_Tstar as fanno_T_Tstar  # noqa: F401 (public DEIMOS name)
from rdelib.units import PA_PER_BAR

_EPS = np.finfo(float).eps
_BRENT_XTOL = 1e-12
_BRENT_RTOL = 4 * _EPS  # brentq minimum; comparable to MATLAB fzero's TolX = eps
PSI_PER_BAR = 14.504  # conversion used in the MATLAB printout (rdelib.units.PSI_PER_BAR is 14.5038)
N_PROPELLANTS = 2  # [oxygen, methane]; the MATLAB loops `for i = 1:2`

# ---------------------------------------------------------------------------
# Flow-relation helpers
#
# Shared with rdelib.gasdynamics (same operation order, so bit-identical):
# sutherland_viscosity (a wrapper keeping the DEIMOS keyword S_ref),
# colebrook_darcy (bracketed Colebrook), fanno_T_Tstar
# (= fanno_T_over_Tstar), fanno_fLD_max, critical_pressure_ratio, and
# orifice_mdot (a wrapper around subcritical_orifice_mass_flow).
# fanno_fLD, fanno_P_Pstar and fanno_mach_from_fLD stay local: rdelib's
# versions use a different operation order (np.log, another P/P* form) or a
# different bracket, so they agree only to round-off.
# ---------------------------------------------------------------------------


def sutherland_viscosity(T, mu_ref, T_ref, S_ref):
    """Sutherland's law: mu = mu_ref (T_ref + S)/(T + S) (T/T_ref)^(3/2).

    DEIMOS signature (``S_ref``); the formula is
    :func:`rdelib.gasdynamics.sutherland_viscosity`.
    """
    return _gd.sutherland_viscosity(T, mu_ref, T_ref, S_ref)


def colebrook_residual(ff: float, Re: float, epsilon: float, D: float) -> float:
    """MATLAB ``colebrook`` handle: 1/sqrt(f) + 2 log10(eps/(3.7 D) + 2.51/(Re sqrt f))."""
    return 1 / math.sqrt(ff) + 2 * math.log10(epsilon / (3.7 * D) + 2.51 / (Re * math.sqrt(ff)))


def fanno_fLD(M: float, gamma: float) -> float:
    """Darcy Fanno parameter f L*/D (= 4 f_Fanning L*/D) at Mach M."""
    return (1 - M**2) / (gamma * M**2) + ((gamma + 1) / (2 * gamma)) * math.log(
        ((gamma + 1) * M**2) / (2 + (gamma - 1) * M**2)
    )


def fanno_P_Pstar(M, gamma):
    """Fanno static-pressure ratio P/P* = (1/M) / sqrt(2/(gamma+1) (1 + (gamma-1)/2 M^2))."""
    return (1 / M) * (1 / np.sqrt((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M**2)))


def fanno_mach_from_fLD(fLD: float, gamma: float, M_guess: float) -> float:
    """Mach number with Darcy Fanno parameter f L*/D = ``fLD``.

    The branch is chosen by ``M_guess`` (subsonic if < 1, supersonic if > 1),
    mirroring the role of the fsolve starting guess in the MATLAB.
    """

    def fun(M):
        return fanno_fLD(M, gamma) - fLD

    if M_guess == 1:
        raise ValueError("Don't guess 1 for M_G: give either a subsonic or a supersonic guess.")
    if M_guess < 1:
        lo, hi = 1e-6, 1.0
    else:
        if fLD >= fanno_fLD_max(gamma):
            raise ValueError(
                f"No supersonic Fanno solution: f*L/D = {fLD:.6g} >= limit {fanno_fLD_max(gamma):.6g} "
                f"for gamma = {gamma}."
            )
        lo, hi = 1.0, 10.0
        while fun(hi) < 0:
            hi *= 10.0
            if hi > 1e12:
                raise ValueError("Supersonic Fanno root could not be bracketed.")
    return brentq(fun, lo, hi, xtol=1e-15, rtol=_BRENT_RTOL, maxiter=500)


def orifice_mdot(P0, P, Cd, A, gamma, R, T0):
    """Subcritical isentropic orifice mass flow (kg/s).

    mdot = Cd A P0 sqrt( 2 gamma / (R T0 (gamma-1)) (P/P0)^(2/gamma) (1 - (P/P0)^((gamma-1)/gamma)) )
    Valid for P/P0 >= critical ratio (2/(gamma+1))^(gamma/(gamma-1)).
    Argument order of the MATLAB handle; the formula is
    :func:`rdelib.gasdynamics.subcritical_orifice_mass_flow`.
    """
    return subcritical_orifice_mass_flow(Cd, A, P0, P, T0, gamma, R)


def solve_P01(m_dot: float, P2: float, Cd: float, A: float, gamma: float, R: float, T: float) -> float:
    """Upstream stagnation pressure that drives ``m_dot`` through the orifice.

    Same solution window as the MATLAB: [1.01 P2, 1.01 P2 + 1e8 Pa] (1000 bar).
    """

    def P01_fun(P0):
        return orifice_mdot(P0, P2, Cd, A, gamma, R, T) - m_dot

    a = P2 * 1.01
    b = a + 1e8  # 1000 bar (1e8 Pa) Solution Window
    if P01_fun(a) * P01_fun(b) < 0:
        return brentq(P01_fun, a, b, xtol=_BRENT_XTOL, rtol=_BRENT_RTOL, maxiter=500)
    raise ValueError("ERROR: P01 Solver Solution Window Does Not Include A Valid Root")


# ---------------------------------------------------------------------------
# Main calculation
# ---------------------------------------------------------------------------


@dataclass
class FannoFlowResult:
    """Outputs of :func:`DEIMOS_sys_FannoFlowCalc_v2` (names follow the MATLAB variables).

    Two-element arrays are ordered [oxygen, methane].
    """

    OF_ratio: float  # oxygen-to-fuel mass ratio
    Y: np.ndarray  # propellant mass fractions
    m_dot: np.ndarray  # propellant mass flow rates (kg/s)
    m_dotT: float  # total mass flow rate (kg/s)
    A_inj: np.ndarray  # total injector orifice area per propellant (m^2)
    A_t: float  # throat area (mm^2)
    A_gap: float  # annulus gap (mm)
    Rinner: float  # inner throat radius (mm)
    Router: float  # outer throat radius (mm)
    iter: int  # outer iterations performed
    err: float  # final max relative change in P2
    mu: np.ndarray  # Sutherland viscosity of the last iteration (Pa-s)
    Re: np.ndarray  # Reynolds number of the last iteration
    f: np.ndarray  # Darcy friction factor of the last iteration
    M_TubeIN: np.ndarray  # tube inlet Mach number
    P01: np.ndarray  # upstream stagnation pressure (Pa) (solved with the P2 of the last iteration's start)
    P2: np.ndarray  # tube-inlet ("downstream") static pressure after the last update (Pa)
    T_stag: np.ndarray  # stagnation temperature after the last update (K)
    P01_residual: np.ndarray  # orifice MFR(P01; final P2) - m_dot (kg/s); MATLAB printed only the last entry
    r_crit: np.ndarray  # critical orifice pressure ratio per propellant
    orifice_choked: np.ndarray  # True where P2/P01 < r_crit (subcritical orifice equation out of range)
    # inputs needed by the P01 diagnostic
    Cd: float
    gamma: np.ndarray
    R: np.ndarray
    T: float


def DEIMOS_sys_FannoFlowCalc_v2(
    Thrust: float = 1780.0,  # (N)
    Isp: float = 233.37,  # Specific impulse, from CEA (s)
    Pc: float = 1e6,  # Chamber pressure (Pa)
    cstar: float = 1806.2,  # Characteristic velocity, from CEA (m/s)
    phi: float = 1.3,  # Equivalence ratio
    Cd: float = 0.6,  # Discharge coefficient
    elements: float = 36,  # Number of injector elements
    Imp_angle: float = 45.0,  # Injector impingement angle, off axial (deg) -- unused, as in MATLAB
    ope: Sequence[float] = (1, 2),  # Orifices per injector element
    inj_OD: Sequence[float] = (2.5, 1.0),  # Injector diameter (mm)
    T: float = 283.15,  # Static inlet temperature (K)
    Rmid: float = 22.3524953,  # Mean annulus radius (mm) (MATLAB comment says "diameter"; used as a radius)
    M_G: Sequence[float] = (0.5, 0.5),  # Speed regime guess for FFC, used in Re calc as well
    P_G: Sequence[float] = (3e6, 3e6),  # Pressure guess for Re (Pa)
    T_stagG: Sequence[float] = (300.0, 300.0),  # Temperature guess for Re (K)
    f_G: Sequence[float] = (0.02, 0.02),  # Darcy friction factor guess (MATLAB fsolve start; unneeded by brentq)
    L: Sequence[float] = (9.48944, 9.56818),  # Tube length (mm, same units as D)
    D: Sequence[float] | None = None,  # Tube diameter (mm); defaults to inj_OD
    g: float = 9.81,  # (m/s^2)
    gamma: Sequence[float] = (1.4, 1.31),  # specific heat ratio
    R: Sequence[float] = (259.8, 518.3),  # Gas constants (J/kg-K)
    mu_ref: Sequence[float] = (2.07e-5, 1.1e-5),  # Reference viscosity @ reference temp (Pa-s)
    T_ref: Sequence[float] = (291.0, 291.0),  # Reference temperature (K)
    S_ref: Sequence[float] = (127.0, 148.0),  # Reference Sutherland Constant (K)
    epsilon: float = 0.015,  # Stainless Steel surface roughness (mm)
    tol: float = 1e-6,
    max_iter: int = 100,
) -> FannoFlowResult:
    """Injector Fanno-flow / feed-pressure calculation (defaults = MATLAB v2 inputs).

    Per-propellant inputs are [oxygen, methane]. ``ope``, ``T_stagG``, ``mu_ref``,
    ``T_ref`` and ``S_ref`` may also be scalars, which are applied to both
    propellants, as MATLAB broadcasting does. The other per-propellant inputs are
    indexed per propellant in the MATLAB and must have exactly two entries.
    """
    # Imp_angle was declared but never used in the MATLAB. f_G was the MATLAB
    # Colebrook fsolve start point; the bracketed solver here does not need it.
    del Imp_angle

    # Iteration controls: the MATLAB loop starts with err = 1 and runs while
    # err > tol && iter < max_iter, so tol >= 1 or max_iter <= 0 runs zero
    # iterations and the MATLAB then fails on undefined P01. A fractional
    # max_iter > 0 is valid (0.5 -> one iteration), as in MATLAB.
    if not tol < 1:
        raise ValueError(f"tol must be < 1 (the iteration starts with err = 1); got {tol!r}")
    if not max_iter > 0:
        raise ValueError(f"max_iter must be > 0; got {max_iter!r}")

    def per_propellant(name, v, scalar_ok=False):
        a = np.array(v, dtype=float).ravel()
        if scalar_ok and a.size == 1:  # MATLAB scalar broadcasting (checked against Octave)
            return np.full(N_PROPELLANTS, a[0])
        if a.size != N_PROPELLANTS:
            raise ValueError(
                f"{name} must have exactly {N_PROPELLANTS} entries [oxygen, methane]"
                f"{' or be a scalar' if scalar_ok else ''}; got {a.size}"
            )
        return a

    # Scalars allowed: used only in element-wise expressions in the MATLAB.
    ope = per_propellant("ope", ope, scalar_ok=True)
    T_stagG = per_propellant("T_stagG", T_stagG, scalar_ok=True)
    mu_ref = per_propellant("mu_ref", mu_ref, scalar_ok=True)
    T_ref = per_propellant("T_ref", T_ref, scalar_ok=True)
    S_ref = per_propellant("S_ref", S_ref, scalar_ok=True)
    # Exactly two entries: indexed X(i) in the MATLAB (a scalar fails with "X(2): out of bound 1").
    inj_OD = per_propellant("inj_OD", inj_OD)
    M_G = per_propellant("M_G", M_G)
    P_G = per_propellant("P_G", P_G)
    per_propellant("f_G", f_G)  # validated for interface parity; not needed by the solver
    L = per_propellant("L", L)
    D = inj_OD.copy() if D is None else per_propellant("D", D)
    gamma = per_propellant("gamma", gamma)
    R = per_propellant("R", R)
    n = N_PROPELLANTS

    # ---------------- calc ----------------
    OF_ratio = 4.0 / phi  # Oxygen to fuel ratio
    Y = np.array([OF_ratio / (1 + OF_ratio), 1 / (1 + OF_ratio)])  # Propellant mass fraction
    m_dot = (Y * Thrust) / (Isp * g)  # Propellant Mass flow rate (kg/s)
    m_dotT = Thrust / (Isp * g)  # Total mass flow rate (kg/s)
    A_inj = (inj_OD / 2000) ** 2 * np.pi * elements * ope  # Injector area (m^2)
    A_t = ((m_dotT * cstar) / Pc) * 1e6  # Throat area (mm^2)
    A_gap = A_t / (Rmid * 2 * np.pi)  # Annulus gap (mm)
    Rinner = Rmid - A_gap / 2  # Inner throat radius (mm)
    Router = Rmid + A_gap / 2  # Outer throat radius (mm)

    # ---------------- Iteration Setup ----------------
    err = 1.0
    it = 0

    P2 = P_G.copy()  # initial guess (Pa)
    T_stag = T_stagG.copy()  # initial guess (K)
    f = np.zeros(n)
    M_TubeIN = np.zeros(n)
    P01 = np.zeros(n)
    mu = Re = None

    while err > tol and it < max_iter:
        it += 1

        # ---- Update viscosity (Sutherland) ----
        mu = sutherland_viscosity(T_stag, mu_ref, T_ref, S_ref)

        # ---- Reynolds number ----
        Re = (P2 * M_G * np.sqrt(gamma * R * T) * (D / 1000)) / (R * T * mu)

        # ---- Loop over propellants ----
        for i in range(n):
            # ---- Friction factor via Colebrook ----
            f[i] = colebrook_darcy(Re[i], epsilon, inj_OD[i])

            # ---- Solve Fanno for inlet Mach ----
            M_TubeIN[i] = fanno_mach_from_fLD(f[i] * L[i] / D[i], gamma[i], M_G[i])

            if M_TubeIN[i] <= 0:
                raise ValueError("ERROR: Fanno Solver returned nonphysical Mach number.")

            # ---- Solve for upstream pressure P01 ----
            P01[i] = solve_P01(m_dot[i], P2[i], Cd, A_inj[i], gamma[i], R[i], T)

        # ---- Update stagnation properties ----
        T_stag_new = T * (1 + ((gamma - 1) / 2) * M_TubeIN**2)

        P_new = Pc * fanno_P_Pstar(M_TubeIN, gamma)

        # ---- Convergence check (pressure only) ----
        err = float(np.max(np.abs((P_new - P2) / P2)))

        # Update guesses
        P2 = P_new
        T_stag = T_stag_new

    if it == 0:  # unreachable after the checks above; never return unsolved placeholders
        raise RuntimeError("No Fanno iteration executed: P01 / M_TubeIN were never solved.")

    if it == max_iter:
        warnings.warn("Iteration did not fully converge", stacklevel=2)

    # P01 solver residual, evaluated as in the MATLAB diagnostic (final P2, solved P01)
    P01_residual = orifice_mdot(P01, P2, Cd, A_inj, gamma, R, T) - m_dot

    # Orifice validity check (added in the port; the numbers are unchanged, see PORTING_NOTES.md)
    r_crit = critical_pressure_ratio(gamma)
    orifice_choked = (P2 / P01) < r_crit
    for i in np.flatnonzero(orifice_choked):
        warnings.warn(
            f"Propellant {i + 1}: orifice pressure ratio P2/P01 = {P2[i] / P01[i]:.4f} is below the "
            f"critical ratio {r_crit[i]:.4f}; the subcritical orifice equation used for P01 is "
            "outside its range (the orifice would be choked).",
            stacklevel=2,
        )

    return FannoFlowResult(
        OF_ratio=OF_ratio,
        Y=Y,
        m_dot=m_dot,
        m_dotT=m_dotT,
        A_inj=A_inj,
        A_t=A_t,
        A_gap=A_gap,
        Rinner=Rinner,
        Router=Router,
        iter=it,
        err=err,
        mu=mu,
        Re=Re,
        f=f,
        M_TubeIN=M_TubeIN,
        P01=P01,
        P2=P2,
        T_stag=T_stag,
        P01_residual=P01_residual,
        r_crit=r_crit,
        orifice_choked=orifice_choked,
        Cd=Cd,
        gamma=gamma,
        R=R,
        T=T,
    )


# ---------------------------------------------------------------------------
# P01 solver diagnostic (plot + residual) and console output
# ---------------------------------------------------------------------------


def P01_diagnostic(res: FannoFlowResult, n_points: int = 3000):
    """Curves of the MATLAB diagnostic: x (Pa), MFR(x) and required MFR, one row per propellant."""
    x = res.P2[:, None] * (1 + np.logspace(-6, 0.25, n_points))[None, :]
    ym = np.repeat(res.m_dot[:, None], n_points, axis=1)
    y = orifice_mdot(x, res.P2[:, None], res.Cd, res.A_inj[:, None], res.gamma[:, None], res.R[:, None], res.T)
    return x, y, ym


def plot_P01_diagnostic(res: FannoFlowResult):
    """Injector mass flow rate vs upstream stagnation pressure (MATLAB figure). Returns the figure."""
    import matplotlib.pyplot as plt

    x, y, ym = P01_diagnostic(res)
    fig, ax = plt.subplots()
    ax.plot(x[0] / 1e5, y[0], "b", linewidth=2, label="Oxygen MFR($P_0$)")
    ax.plot(x[0] / 1e5, ym[0], "b--", linewidth=1.8, label="Oxygen Required MFR")
    ax.plot(x[1] / 1e5, y[1], "r", linewidth=2, label="Methane MFR($P_0$)")
    ax.plot(x[1] / 1e5, ym[1], "r--", linewidth=1.8, label="Methane Required MFR")
    ax.plot(res.P01[0] / 1e5, res.m_dot[0], "bo", markerfacecolor="b", markersize=7, label="Solved P01 (Oxygen)")
    ax.plot(res.P01[1] / 1e5, res.m_dot[1], "ro", markerfacecolor="r", markersize=7, label="Solved P01 (Methane)")
    ax.grid(True)
    ax.set_xlabel("Upstream Stagnation Pressure, $P_0$ (bar)", fontsize=12)
    ax.set_ylabel("Mass Flow Rate (MFR) (kg/s)", fontsize=12)
    ax.set_title("Injector Mass Flow Rate vs Upstream Stagnation Pressure", fontsize=14)
    ax.legend(loc="upper right")
    return fig


_RULE = "=" * 52


def format_P01_residual(res: FannoFlowResult) -> str:
    """MATLAB 'P01 Solver Residual' block (prints the last propellant's residual, as the MATLAB did)."""
    return (
        f"\n{_RULE}\n                P01 Solver Residual                \n{_RULE}\n\n"
        f"Residual                        : {res.P01_residual[-1]:18.15f} \n\n"
    )


def format_summary(res: FannoFlowResult) -> str:
    """MATLAB 'FLOW CALCULATION SUMMARY' block, character for character."""

    def bar_psi(p):
        return f"{p / PA_PER_BAR:.2f} bar | {(p * PSI_PER_BAR) / PA_PER_BAR:.2f} psi"

    return (
        f"\n{_RULE}\n              FLOW CALCULATION SUMMARY              \n{_RULE}\n\n"
        "--- Mass Flow Rates ---\n"
        f"Total Mass Flow Rate            : {res.m_dotT:.3f} kg/s\n"
        f"Oxygen Mass Flow Rate           : {res.m_dot[0]:.3f} kg/s\n"
        f"Methane Mass Flow Rate          : {res.m_dot[1]:.3f} kg/s\n\n"
        "--- Annular Throat Geometry ---\n"
        f"Outer Radius                    : {res.Router:.2f} mm\n"
        f"Inner Radius                    : {res.Rinner:.2f} mm\n"
        f"Annulus Gap                     : {res.A_gap:.2f} mm\n\n"
        "--- Fanno Output & Pressure Correlations ---\n"
        f"Oxygen Inlet Mach               : {res.M_TubeIN[0]:.2f} \n"
        f"Methane Inlet Mach              : {res.M_TubeIN[1]:.2f} \n"
        f"P2 Oxygen Downstream Pressure   : {bar_psi(res.P2[0])}\n"
        f"P2 Methane Downstream Pressure  : {bar_psi(res.P2[1])}\n"
        f"P01 Oxygen Upstream Pressure    : {bar_psi(res.P01[0])}\n"
        f"P01 Methane Upstream Pressure   : {bar_psi(res.P01[1])}\n\n"
        f"{_RULE}\n\n"
    )


def main(diagnostic: bool = True, show: bool = True) -> FannoFlowResult:
    """Run with the MATLAB inputs, print the MATLAB output and (optionally) show the diagnostic plot."""
    res = DEIMOS_sys_FannoFlowCalc_v2()
    if diagnostic:
        plot_P01_diagnostic(res)
        print(format_P01_residual(res), end="")
    print(format_summary(res), end="")
    if diagnostic and show:
        import matplotlib.pyplot as plt

        plt.show()
    return res


if __name__ == "__main__":
    main()
