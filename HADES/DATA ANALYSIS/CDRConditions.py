import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from scipy.optimize import brentq

# === Fixed Inputs ===
m_dot_air  = 0.907                        # kg/s
phi        = np.array([0.8, 0.9, 1.0, 1.1, 1.2])
OF_stoich  = 34.3
m_dot_h2   = m_dot_air / OF_stoich * phi  # kg/s

# === Geometry ===
D_outer    = 4 * 0.0254                   # m (4 inch outer diameter)
w_choked   = 1e-3                         # m (1mm choked slot width)
w_annulus  = 9.525e-3                     # m (9.525mm annulus width)
n_h2       = 32                           # number of H2 orifices
D_h2_orif  = 1e-3                         # m (1mm orifice diameter)

R_outer    = D_outer / 2
r_inner_i  = R_outer - w_choked
r_inner_w  = R_outer - w_annulus

A_air      = np.pi * (R_outer**2 - r_inner_i**2)
A_w        = np.pi * (R_outer**2 - r_inner_w**2)
A_h2_orif  = n_h2 * np.pi * (D_h2_orif / 2)**2

Ai_Aw_air  = A_air / A_w
Ai_Aw_h2   = A_h2_orif / A_w

A_h2       = 25.13e-6                     # m^2
Cd_air     = 0.62
Cd_h2      = 0.71
T_0_C      = 20.0
T_0_K      = T_0_C + 273.15
P_REF      = 101325.0

# === Fanno Inputs ===
L          = 0.5 * 0.0254                 # m (0.5 inches)
Dh_air     = 2e-3                         # m
Dh_h2      = 1e-3                         # m
epsilon    = 1.5e-6                       # m
M_0        = 0.65

# === RDE Inputs ===
cj_vel     = 1920.0                       # m/s
wave_modes = np.arange(1, 6)             # 1 to 5 waves
p_p1       = 16.0                         # pressure ratio from NASA CEA
D_mid      = D_outer - w_annulus
pathlength = np.pi * D_mid               # m

# === Fluid Properties ===
def get_props(fluid, T_K, P):
    Cp  = CP.PropsSI('Cpmass',   'T', T_K, 'P', P, fluid)
    Cv  = CP.PropsSI('Cvmass',   'T', T_K, 'P', P, fluid)
    mu  = CP.PropsSI('viscosity','T', T_K, 'P', P, fluid)
    return Cp / Cv, Cp - Cv, mu

g_air, R_air, mu_air = get_props('Air',      T_0_K, P_REF)
g_h2,  R_h2,  mu_h2  = get_props('Hydrogen', T_0_K, P_REF)

# === Choked Flow ===
def MFP(gamma):
    return np.sqrt(gamma) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))

def P0_plenum(m_dot, R, T_K, Cd, A, gamma):
    return m_dot * np.sqrt(R * T_K) / (Cd * A * MFP(gamma)) / 1e5  # bar

# === Fanno Flow ===
def fanno_pdrop(P0_bar, D, L, epsilon, m_dot, gamma, mu):
    P0        = P0_bar * 1e5
    Re        = (4 * m_dot) / (np.pi * D * mu)
    rel_rough = epsilon / D
    A0        = -0.79638 * np.log((rel_rough / 8.208) + (7.3357 / Re))
    A1        = Re * rel_rough + 9.3120665 * A0
    num       = 8.128943 + A1
    den       = 8.128943 * A0 - 0.86859209 * A1 * np.log(A1 / (3.7099535 * Re))
    f         = (num / den) ** 2
    fanno_param  = (f * L) / D
    func         = lambda Ma: ((1 - Ma**2) / (gamma * Ma**2) + (gamma + 1) / (2 * gamma)
                               * np.log((gamma + 1) * Ma**2 / (2 + (gamma - 1) * Ma**2))
                               - fanno_param)
    Ma_inlet     = brentq(func, 1e-9, 1.0 - 1e-9)
    P_inlet      = P0 * (1 + (gamma - 1) / 2 * Ma_inlet**2) ** (-gamma / (gamma - 1))
    P_ratio_star = (1 / Ma_inlet) * np.sqrt((gamma + 1) / (2 + (gamma - 1) * Ma_inlet**2))
    P_exit       = P_inlet / P_ratio_star
    return P_inlet / 1e5, P_exit / 1e5, Ma_inlet, f

# === Kessler Chamber Pressure ===
def kessler_P0(P_p_bar, Ai_Aw, M_0, gamma):
    P_c_k    = (2 / (gamma + 1)) ** (gamma / (gamma - 1)) * P_p_bar
    P_exit_k = Ai_Aw * (1 / M_0) * np.sqrt((gamma + 1) / (2 + (gamma - 1) * M_0**2)) * P_c_k
    return P_c_k, P_exit_k

# === RDE Average Chamber Pressure ===
def rde_chamber_pressure(p_p1, P_i_arr, phi_arr, pathlength, cj_vel, wave_modes):
    """
    P_c(t) = P_CJ * (P_i / P_CJ)^(t / t_cycle)
    Rows = phi, Cols = wave modes
    """
    fig, axes = plt.subplots(len(phi_arr), len(wave_modes),
                             figsize=(18, 14), constrained_layout=True)
    fig.suptitle('RDE Chamber Pressure vs Time', fontsize=14, fontweight='medium')

    avg_pressures = np.zeros((len(phi_arr), len(wave_modes)))

    for i, (ph, P_i) in enumerate(zip(phi_arr, P_i_arr)):
        P_CJ = p_p1 * P_i
        for j, wm in enumerate(wave_modes):
            ax         = axes[i, j]
            cycle_time = (pathlength / cj_vel) * 1e3 / wm            # ms
            time       = np.linspace(0, cycle_time, 1000)            # ms
            P_c        = P_CJ * (P_i / P_CJ) ** (time / cycle_time)  # bar
            avg_p      = np.mean(P_c)
            avg_pressures[i, j] = avg_p

            ax.plot(time, P_c, linewidth=1.5)
            ax.axhline(avg_p, color='red', linestyle='--', linewidth=1.2,
                       label=f'Avg = {avg_p:.2f} bar')
            ax.axhline(P_i, color='green', linestyle=':', linewidth=1.0,
                       label=f'P_i = {P_i:.2f} bar')
            ax.set_title(f'$\\phi$ = {ph:.1f}  |  {wm} Wave{"s" if wm > 1 else ""}', fontsize=9)
            ax.set_xlabel('Time (ms)', fontsize=8)
            ax.set_ylabel('$P_c$ (bar)', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)

    plt.show()
    return avg_pressures

# === CDR Conditions ===
P0_air = P0_plenum(m_dot_air, R_air, T_0_K, Cd_air, A_air, g_air)
P0_h2  = P0_plenum(m_dot_h2,  R_h2,  T_0_K, Cd_h2,  A_h2,  g_h2)

P_in_air, P_ex_air, Ma_air, f_air = fanno_pdrop(P0_air, Dh_air, L, epsilon, m_dot_air, g_air, mu_air)

fanno_h2 = np.array([fanno_pdrop(p, Dh_h2, L, epsilon, m, g_h2, mu_h2)
                      for p, m in zip(P0_h2, m_dot_h2)])
P_in_h2  = fanno_h2[:, 0]
P_ex_h2  = fanno_h2[:, 1]
Ma_h2    = fanno_h2[:, 2]
f_h2     = fanno_h2[:, 3]

# Kessler uses plenum stagnation pressure (P0) as input
P_c_k_air, P_exit_k_air = kessler_P0(P0_air, Ai_Aw_air, M_0, g_air)
P_c_k_h2,  P_exit_k_h2  = np.array([kessler_P0(p, Ai_Aw_h2, M_0, g_h2) for p in P0_h2]).T

# P_i = average of air and H2 kessler exit pressures for each phi
P_i_arr = (m_dot_air * P_exit_k_air + m_dot_h2 * P_exit_k_h2) / (m_dot_air + m_dot_h2)

avg_chamber_pressures = rde_chamber_pressure(p_p1, P_i_arr, phi, pathlength, cj_vel, wave_modes)

# === Print Results ===
print(f"=== CDR Design Point Conditions (T0 = {T_0_C} °C) ===\n")
print(f"  Geometry:  L = {L*1000:.1f} mm  |  Dh_air = {Dh_air*1000:.1f} mm  |  Dh_H2 = {Dh_h2*1000:.1f} mm  |  epsilon = {epsilon*1e6:.1f} um")
print(f"             A_air = {A_air*1e6:.2f} mm^2  |  A_w = {A_w*1e6:.2f} mm^2  |  A_H2 = {A_h2_orif*1e6:.2f} mm^2")
print(f"             Ai_Aw_air = {Ai_Aw_air:.4f}  |  Ai_Aw_H2 = {Ai_Aw_h2:.4f}")
print(f"             Pathlength = {pathlength*1000:.1f} mm  |  CJ vel = {cj_vel} m/s  |  p_p1 = {p_p1}\n")
print(f"  Air:  m_dot = {m_dot_air:.4f} kg/s  |  Re = {4*m_dot_air/(np.pi*Dh_air*mu_air):.0f}  |  f = {f_air:.4f}")
print(f"        P0 = {P0_air:.2f} bar  |  P_inlet = {P_in_air:.2f} bar  |  P_exit_c = {P_ex_air:.2f} bar  |  Ma_inlet = {Ma_air:.4f}")
print(f"        P_c_k_air = {P_c_k_air:.2f} bar  |  P_exit_k_air = {P_exit_k_air:.2f} bar\n")
print(f"  {'phi':>5}  {'m_dot [kg/s]':>13}  {'P0 [bar]':>10}  {'P_in [bar]':>11}  {'P_exit_c [bar]':>15}  {'Ma_in':>8}  {'f':>8}  {'P_c_k [bar]':>12}  {'P_exit_k_h2 [bar]':>18}  {'P_i [bar]':>10}  {'P_CJ [bar]':>11}")
print(f"  {'-'*130}")
for p, m, P0h, Pi, Pe, Ma, fh, Pck, Pek, P_i, P_CJ in zip(
        phi, m_dot_h2, P0_h2, P_in_h2, P_ex_h2, Ma_h2, f_h2,
        P_c_k_h2, P_exit_k_h2, P_i_arr, p_p1 * P_i_arr):
    print(f"  {p:>5.1f}  {m:>13.5f}  {P0h:>10.2f}  {Pi:>11.2f}  {Pe:>15.2f}  {Ma:>8.4f}  {fh:>8.4f}  {Pck:>12.2f}  {Pek:>18.2f}  {P_i:>10.2f}  {P_CJ:>11.2f}")
print(f"\n  RDE Average Chamber Pressures (bar):")
print(f"  {'phi':<8}", end="")
for wm in wave_modes:
    print(f"  {wm} wave{'s' if wm > 1 else ''}  ", end="")
print()
print(f"  {'-'*60}")
for i, p in enumerate(phi):
    print(f"  {p:<8.1f}", end="")
    for j in range(len(wave_modes)):
        print(f"  {avg_chamber_pressures[i,j]:>8.2f}  ", end="")
    print()