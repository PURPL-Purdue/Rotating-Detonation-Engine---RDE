"""Compressible Subcritical Orifice Flow Calculator.

Python port of ``DEIMOS_sys_CdCalc_v1.m``.

WARNING - This is an empirical formulation, as D1 --> D2 error explodes.

Method (unchanged from the MATLAB):

1. Check choking with the critical pressure ratio ((gamma+1)/2)^(gamma/(gamma-1)).
2. Upstream stagnation pressure P0 from isentropic continuity between the pipe
   (D1, P1) and the orifice bore (D2, P2).
3. Fixed-point iteration of the empirical corner-tap discharge coefficient
   correlation (pipe diameter in inches inside the correlation), with the
   incompressible mass flow feeding the pipe Reynolds number.
4. Loss coefficient f = 1/Cd - 1/(2 Cd^2), isentropic flow function Kn(P2/P0)
   and the compressible discharge coefficient Cd_comp.
5. Compressible mass flow  mdot_comp = Kn * Cd_comp * A2 * sqrt(P0 * rho0).

Units: SI (m, Pa, K, kg/s, Pa-s).

Run ``uv run python -m DEIMOS.DEIMOS_sys_CdCalc_v1`` to print the MATLAB results.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

IN_TO_M = 0.0254  # (m/in)


@dataclass
class CdCalcResult:
    """Outputs of :func:`DEIMOS_sys_CdCalc_v1` (names follow the MATLAB variables)."""

    Cd_inc: float  # converged incompressible Cd (corner tap)
    Cd_comp: float  # compressible Cd
    mdot_comp: float  # compressible mass flow rate (kg/s)
    Re_D1: float  # upstream-pipe Reynolds number from the last Cd iteration
    A1: float  # upstream pipe area (m^2)
    A2: float  # orifice bore area (m^2)
    beta: float  # D2/D1
    D1_in: float  # upstream pipe diameter (in)
    PR_critical: float  # critical pressure ratio P1/P2
    PR_actual: float  # actual pressure ratio P1/P2
    P0: float  # upstream stagnation pressure (Pa)
    rho0: float  # upstream stagnation density (kg/m^3)
    mdot: float  # incompressible mass flow of the last Cd iteration (kg/s)
    f: float  # loss coefficient
    r: float  # P2/P0
    Kn: float  # isentropic flow function
    term: float  # discriminant of the Cd_comp quadratic
    iter: int  # number of Cd iterations performed (1-based, as in MATLAB)


def corner_tap_Cd(beta: float, D1_in: float, Re_D1: float) -> float:
    """Empirical corner-tap orifice discharge coefficient (D1 in inches).

    Same expression as the MATLAB ``Cv / sqrt(1 - beta^4)``.
    """
    Cv = (
        0.5991 + 0.0044 / D1_in + (0.3155 + 0.0175 / D1_in) * (beta**4 + 2 * beta**16)
    ) * math.sqrt(1 - beta**4) + (
        0.52 / D1_in - 0.192 + (16.48 - 1.16 / D1_in) * (beta**4 + 4 * beta**16)
    ) * math.sqrt((1 - beta**4) / Re_D1)
    return Cv / math.sqrt(1 - beta**4)


def DEIMOS_sys_CdCalc_v1(
    D1: float = 0.01125,  # Upstream pipe diameter (m)
    D2: float = 0.00225,  # Orifice bore diameter (m)
    P1: float = 19e5,  # Upstream static pressure (Pa)
    P2: float = 13e5,  # Downstream static pressure (Pa)
    T1: float = 283.15,  # Upstream temperature (K)
    gamma: float = 1.4,  # Ratio of specific heats
    Rg: float = 259.8,  # Gas constant (J/kg-K)
    mu: float = 2.35e-5,  # Dynamic viscosity (Pa-s)
    Cd_guess: float = 0.61,  # Initial Cd guess (look at slide 17)
    tol: float = 1e-6,
    maxiter: float = 100,
) -> CdCalcResult:
    """Compressible subcritical orifice Cd and mass flow (defaults = MATLAB inputs)."""
    if not maxiter >= 1:
        raise ValueError("maxiter must be >= 1")
    # MATLAB's `for iter = 1:maxiter` runs floor(maxiter) times, so accept non-integer limits.
    n_iter = math.floor(maxiter)

    # ---------------- Geometry ----------------
    A2 = math.pi / 4 * D2**2
    A1 = math.pi / 4 * D1**2
    beta = D2 / D1

    D1_in = D1 / IN_TO_M

    # ---------------- Check choking ----------------
    PR_critical = ((gamma + 1) / 2) ** (gamma / (gamma - 1))
    PR_actual = P1 / P2

    if PR_actual >= PR_critical:
        warnings.warn("Flow is CHOKED - this model is invalid.", stacklevel=2)

    # ---------------- Stagnation pressure ----------------
    P0 = (
        ((D1 / D2) ** 4 * P1 ** ((gamma + 1) / gamma) - P2 ** ((gamma + 1) / gamma))
        / ((D1 / D2) ** 4 * P1 ** (2 / gamma) - P2 ** (2 / gamma))
    ) ** (gamma / (gamma - 1))

    rho0 = P0 / (Rg * T1)

    # ---------------- Iterate Cd (corner tap) ----------------
    Cd_inc = Cd_guess
    for it in range(1, n_iter + 1):
        mdot = Cd_inc * A2 * math.sqrt(2 * (P1 / (Rg * T1)) * (P1 - P2))

        Re_D1 = 4 * mdot / (math.pi * D1 * mu)

        Cd_new = corner_tap_Cd(beta, D1_in, Re_D1)

        if abs(Cd_new - Cd_inc) / Cd_inc < tol:
            Cd_inc = Cd_new
            break

        Cd_inc = Cd_new

    # ---------------- Loss coefficient ----------------
    f = 1 / Cd_inc - 1 / (2 * Cd_inc**2)

    # ---------------- Flow coefficient ----------------
    r = P2 / P0

    Kn = math.sqrt((2 * gamma / (gamma - 1)) * r ** (2 / gamma) * (1 - r ** ((gamma - 1) / gamma)))

    # ---------------- Compressible Cd ----------------
    term = 1 - f * (2 * r ** (1 / gamma)) ** 2 * (1 - r) / Kn**2

    if term < 0:
        raise ValueError("Negative discriminant - check inputs.")

    Cd_comp = (1 - math.sqrt(term)) / (2 * f * r ** (1 / gamma))

    # ---------------- Mass flow ----------------
    mdot_comp = Kn * Cd_comp * A2 * math.sqrt(P0 * rho0)

    return CdCalcResult(
        Cd_inc=Cd_inc,
        Cd_comp=Cd_comp,
        mdot_comp=mdot_comp,
        Re_D1=Re_D1,
        A1=A1,
        A2=A2,
        beta=beta,
        D1_in=D1_in,
        PR_critical=PR_critical,
        PR_actual=PR_actual,
        P0=P0,
        rho0=rho0,
        mdot=mdot,
        f=f,
        r=r,
        Kn=Kn,
        term=term,
        iter=it,
    )


def format_results(res: CdCalcResult) -> str:
    """Console output identical to the MATLAB ``fprintf`` block."""
    return (
        "\n--- Results ---\n"
        f"Cd (incompressible) = {res.Cd_inc:.6f}\n"
        f"Cd (compressible)   = {res.Cd_comp:.6f}\n"
        f"Mass flow rate      = {res.mdot_comp:.6f} kg/s\n"
        f"Re_D1               = {res.Re_D1:.0f}\n"
    )


def main() -> CdCalcResult:
    res = DEIMOS_sys_CdCalc_v1()
    print(format_results(res), end="")
    return res


if __name__ == "__main__":
    main()
