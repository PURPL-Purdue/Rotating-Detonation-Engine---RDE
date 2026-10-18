"""Port of ``HADES_size_thermal_v1.m`` (chamber hot-fire wall-temperature limit).

Author of the MATLAB: Ansh Patel. 1-D transient radial conduction in the
outer wall with gas-side convection (Ruan et al. recovery temperature and a
Bartz-type Nusselt correlation) and ambient convection outside. MATLAB's
``pdepe`` is replaced by :func:`rdelib.ode.pdepe` (same Skeel-Berzins
spatial discretisation, stiff BDF time integration).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from rdelib.ode import pdepe

__all__ = ["HADES_size_thermal_v1", "ThermalResult", "heatPDE", "icFun", "bcFun", "main",
           "PHI_VEC", "CJ_T_VEC", "RHO_G_VEC", "CP_G_KJ_VEC", "GAMMA_G_VEC", "PR_T_G_VEC",
           "VISC_MP_VEC", "POINT8_CJ"]

# CEA data: 4 phi values, CJ state for each (columns of the MATLAB tables)
PHI_VEC = [0.8, 0.9, 1.0, 1.1]
CJ_T_VEC = [2790.67, 2914.30, 2990.66, 3019.33]  # K
RHO_G_VEC = [3.3327, 3.1250, 2.9747, 2.8724]  # kg/m^3
CP_G_KJ_VEC = [2.3336, 2.7540, 3.1354, 3.2067]  # kJ/(kg K), equilibrium
GAMMA_G_VEC = [1.1924, 1.1780, 1.1693, 1.1700]
PR_T_G_VEC = [0.5242, 0.4850, 0.4988, 0.5123]  # equilibrium Prandtl
VISC_MP_VEC = [0.89074, 0.91898, 0.93556, 0.94054]  # millipoise
POINT8_CJ = [0.6, 0.57, 0.57, 0.60]  # gas temperature as a fraction of CJ T


@dataclass
class ThermalResult:
    phi: float
    T_gasK: float
    T_filmK: float
    T_rec_g: float
    Re_g: float
    Nu_g: float
    h_inner: float
    h_outer: float
    r: np.ndarray  # radial grid [m]
    t: np.ndarray  # time grid [s]
    T: np.ndarray  # T(t_index, r_index) [C]
    T_inner: np.ndarray  # inner-wall temperature vs time [C]
    T_limits_C: np.ndarray
    time_to_limit: np.ndarray  # [s] (NaN if never reached / already above)
    params: dict = field(default_factory=dict)


def heatPDE(r, t, u, dudr, params):
    """rho cp dT/dt = (1/r) d/dr (k r dT/dr): c, f, s."""
    return params["rho"] * params["cp"], params["k"] * dudr, 0.0


def icFun(r, params):
    """Uniform initial temperature [C]."""
    return params["T_init"]


def bcFun(rl, ul, rr, ur, t, params):
    """Inner: -k dT/dr = h_in (T_gas - T); outer: -k dT/dr = h_out (T - T_amb)."""
    pl = params["h_inner"] * (params["T_gas"] - ul)
    ql = 1.0
    pr = params["h_outer"] * (ur - params["T_amb"])
    qr = 1.0
    return pl, ql, pr, qr


def HADES_size_thermal_v1(phiIndex: int = 2, *, plot: bool = True, rtol: float = 1e-6,
                          atol: float = 1e-8) -> ThermalResult:
    """Run the hot-fire wall-temperature script.

    ``phiIndex`` is 1-based like the MATLAB (1 -> phi 0.8, 2 -> 0.9 [default],
    3 -> 1.0, 4 -> 1.1). ``rtol``/``atol`` are the time-integration
    tolerances (see :func:`rdelib.ode.pdepe` for why they are tighter than
    MATLAB's 1e-3/1e-6 defaults).
    """
    # Flow annulus geometry (hot gas region)
    ri_flow = 41.275e-3  # m (1.625")
    ro_flow = 50.800e-3  # m (2.000")

    # Solid wall geometry
    wall_extra = 0.5 * 0.0254  # m (0.5")
    ri = ro_flow
    ro = ro_flow + wall_extra

    T_init = 25  # C
    T_ambC = 25  # C

    T_limits_C = np.array([700, 926.67, 1188.12])
    nLimits = T_limits_C.size
    time_to_limit = np.full(nLimits, np.nan)

    T_wallC_design = 1188.12
    T_wallK_design = T_wallC_design + 273.15

    rho_wall = 7930
    cp_wall = 500
    k_wall = 15

    if isinstance(phiIndex, bool) or not isinstance(phiIndex, (int, np.integer)) or not 1 <= phiIndex <= len(PHI_VEC):
        raise ValueError(f"phiIndex must be an integer in 1..{len(PHI_VEC)} (MATLAB 1-based), got {phiIndex!r}")
    i = int(phiIndex) - 1
    phi = PHI_VEC[i]
    T_CJ = CJ_T_VEC[i]
    rho_g = RHO_G_VEC[i]
    cp_g = 1000 * CP_G_KJ_VEC[i]  # J/(kg K)
    gamma_g = GAMMA_G_VEC[i]
    Pr_t_g = PR_T_G_VEC[i]
    mu_g = VISC_MP_VEC[i] * 1e-4  # Pa s

    T_gasK = POINT8_CJ[i] * T_CJ

    print("\n=== Gas state (phi = %.2f) using its CJ T ===" % phi)
    print("CJ T(phi) = %.2f K" % T_CJ)
    print("T_gasK (0.6CJ) = %.2f K (%.2f °C)" % (T_gasK, T_gasK - 273.15))
    print("rho_g  = %.4f kg/m^3" % rho_g)
    print("cp_g   = %.2f J/(kg·K)" % cp_g)
    print("gamma  = %.4f" % gamma_g)
    print("Pr_t,g = %.4f" % Pr_t_g)
    print("mu_g   = %.4e Pa·s" % mu_g)

    # Film temperature and gas-side h (Ruan / Bartz)
    T_filmK = 0.5 * (T_gasK + T_wallK_design)
    rho_am_g = rho_g * (T_gasK / T_filmK)
    mu_am_g = mu_g

    U_g = 1952.6 * 0.6  # m/s (placeholder in the MATLAB)
    M_g = 0.3  # (placeholder in the MATLAB)

    theta_rec = Pr_t_g ** (1 / 3)
    T_rec_g = T_gasK * (1 + theta_rec * (gamma_g - 1) / 2 * M_g**2)

    k_gas = mu_g * cp_g / Pr_t_g
    D_h = 2 * (ro_flow - ri_flow)
    Re_g = rho_g * U_g * D_h / mu_g
    Nu_g = 0.026 * Re_g**0.8 * Pr_t_g**0.4 * (rho_am_g / rho_g) ** 0.8 * (mu_am_g / mu_g) ** 0.2
    h_inner = Nu_g * k_gas / D_h
    h_outer = 15

    print("T_wall_design = %.1f K (%.1f °C)" % (T_wallK_design, T_wallC_design))
    print("T_filmK = %.1f K" % T_filmK)
    print("T_rec_g  = %.1f K (%.1f °C)" % (T_rec_g, T_rec_g - 273.15))
    print("Re_g = %.3e, Nu_g = %.1f, h_inner = %.1f W/m^2-K" % (Re_g, Nu_g, h_inner))

    params = dict(rho=rho_wall, cp=cp_wall, k=k_wall, h_inner=h_inner, h_outer=h_outer,
                  T_gas=T_rec_g - 273.15, T_amb=T_ambC, T_init=T_init)

    Nr = 120
    r = np.linspace(ri, ro, Nr)
    tEnd = 2
    Nt = 400
    t = np.linspace(0, tEnd, Nt)
    m = 1

    T = pdepe(m, lambda rr, tt, u, dudr: heatPDE(rr, tt, u, dudr, params),
              lambda rr: icFun(rr, params),
              lambda rl, ul, rr_, ur, tt: bcFun(rl, ul, rr_, ur, tt, params),
              r, t, rtol=rtol, atol=atol)

    T_inner = T[:, 0]
    T_inner_F = T_inner * 9 / 5 + 32
    T_limits_F = T_limits_C * 9 / 5 + 32

    for j in range(nLimits):
        Tcrit = T_limits_C[j]
        if (T_inner[0] >= Tcrit) or np.all(T_inner < Tcrit):
            time_to_limit[j] = np.nan
        else:
            # interp1(T_inner, t, Tcrit, 'linear'); T_inner rises monotonically
            time_to_limit[j] = np.interp(Tcrit, T_inner, t)

    print("\n=== Max runtime based on INNER WALL temperature (phi = %.2f) ===" % phi)
    for j in range(nLimits):
        if np.isnan(time_to_limit[j]):
            print("Inner wall never reaches %.2f °C within %.2f s simulated." % (T_limits_C[j], tEnd))
        else:
            print("Time for inner wall to reach %.2f °C (%.2f °F): %.5f s"
                  % (T_limits_C[j], T_limits_F[j], time_to_limit[j]))

    if plot:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(t, T_inner_F, color="#9100FF", linewidth=5)
        for j in range(nLimits):
            plt.axhline(T_limits_F[j], linestyle="--", color="k")
            plt.text(t[-1], T_limits_F[j], "%.2f °F" % T_limits_F[j], ha="right", va="bottom")
        plt.xlabel("Time [s]")
        plt.ylabel("Inner Wall Temperature [°F]")
        plt.title("Inner Wall Temperature vs Time ($\\phi$ = %.2f)" % phi)
        plt.grid(True)

        plt.figure()
        for tp in [0.01, 0.05, 0.1, 0.2, 0.5]:
            idx = int(np.argmin(np.abs(t - tp)))
            plt.plot(r * 1000, T[idx, :] * 9 / 5 + 32, linewidth=1.5, label="t = %.3f s" % t[idx])
        plt.xlabel("Radius [mm]")
        plt.ylabel("Temperature [°F]")
        plt.title("Radial temperature profiles in wall ($\\phi$ = %.2f)" % phi)
        plt.legend(loc="best")
        plt.grid(True)

    return ThermalResult(phi=phi, T_gasK=T_gasK, T_filmK=T_filmK, T_rec_g=T_rec_g, Re_g=Re_g,
                         Nu_g=Nu_g, h_inner=h_inner, h_outer=h_outer, r=r, t=t, T=T,
                         T_inner=T_inner, T_limits_C=T_limits_C, time_to_limit=time_to_limit,
                         params=params)


def main():
    import matplotlib.pyplot as plt

    HADES_size_thermal_v1()
    plt.show()


if __name__ == "__main__":
    main()
