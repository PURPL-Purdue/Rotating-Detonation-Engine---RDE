import numpy as np
import math
import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq

# =============================================================================
# CONVERSIONS
# =============================================================================
psi_to_Pa = 6894.76         # psi -> Pa
in_to_m   = 0.0254          # in -> m
m_to_ft   = 3.28084         # m/s -> ft/s
R_univ    = 8.3144          # J/mol-K; 
P_atm     = 101325.0        # Pa

# =============================================================================
# INPUTS  <- edit operating case values in this block
# =============================================================================

# --- Gas properties ---
MW_O2 = 31.9988 / 1000      # kg/mol
MW_H2 = 2.01588 / 1000      # kg/mol
o_f_target = 8              # target O/F mass ratio
Y_O2 = o_f_target / (o_f_target + 1)   # reference mixture mass fraction
Y_H2 = 1 / (o_f_target + 1)

# --- Operating conditions ---
T     = 10 + 273.15         # K - feed temperature
p1_O2 = 400 * psi_to_Pa     # Pa - illustrative GOx operating case
Cd    = 0.61                # discharge coefficient for both injectors
# GH2 pressure is calculated below to meet the verified O/F = 8 condition.

# --- Injector orifice diameters ---
d_inj_O2 = 0.01 * in_to_m   # m
d_inj_H2 = 0.01 * in_to_m   # m

# --- Supply line geometry ---
d_line = 0.250 * in_to_m    # m - line OD
t_line = 0.049 * in_to_m    # m - line wall thickness

# --- Pre-detonator tube geometry ---
do1     = 0.250 * in_to_m  # m - 1/4 in tube OD (tube 2 in FINAL)
do2     = 0.375 * in_to_m  # m - 3/8 in tube OD (tube 1 in FINAL)
wall_t1 = 0.049            # in - 1/4 in tube wall thickness
wall_t2 = 0.065            # in - 3/8 in tube wall thickness
l1      = 1.5 * in_to_m    # m - 1/4 in tube length
l2      = 2.5 * in_to_m    # m - 3/8 in tube length
Vt      = 0.1407 * in_to_m**3  # m^3 - original assumed CAD volume; NOT in FINAL

# --- Shchelkin spiral ---
OD_spiral = 0.30            # in
ID_spiral = 0.21            # in

# --- Hoop stress ---
P_i           = 3004 / 2   # psia - half of CJ detonation pressure (FINAL)
P_o           = 14.7       # psia - external pressure
allow_1_4     = 7500       # psi - 1/4 in reference stress limit
allow_3_8     = 6500       # psi - 3/8 in reference stress limit

# --- Blowdown simulation ---
dt      = 1e-6             # s - time step
t_final = 0.1              # s - max simulation time

# =============================================================================
# DERIVED GEOMETRY 
# =============================================================================
di1    = do1 - 2 * wall_t1 * in_to_m
di2    = do2 - 2 * wall_t2 * in_to_m
a1     = math.pi * (di1 / 2) ** 2       # m^2 - 1/4 in tube flow area
a2     = math.pi * (di2 / 2) ** 2       # m^2 - 3/8 in tube flow area
A_line = math.pi * ((d_line - 2 * t_line) / 2) ** 2
A_inj_O2 = math.pi * (d_inj_O2 / 2) ** 2
A_inj_H2 = math.pi * (d_inj_H2 / 2) ** 2

# Hoop stress radii in inches
r_o  = 0.250 / 2
r_i  = r_o - wall_t1
r2_o = 0.375 / 2
r2_i = r2_o - wall_t2

# =============================================================================
# SUPPLY PRESSURES 
# =============================================================================
R_O2 = R_univ / MW_O2
R_H2 = R_univ / MW_H2

def choked_mdot(Cd, A, P, g, R, Tt):
    """Isentropic choked mass flow rate, kg/s."""
    return (Cd * A * P / math.sqrt(Tt)) * math.sqrt(g / R) * \
           ((g + 1) / 2) ** (-(g + 1) / (2 * (g - 1)))

def mdot_at_pressure(P, fluid, R, A):
    """Choked-flow relationship; P in Pa."""
    cp_gas = PropsSI('CPMASS', 'T', T, 'P', P, fluid)
    cv_gas = PropsSI('CVMASS', 'T', T, 'P', P, fluid)
    g_gas = cp_gas / cv_gas
    return choked_mdot(Cd, A, P, g_gas, R, T)

def required_GH2_pressure(p_GOx_psia):
    """Same O/F=8 pressure-solving method. Returns psia."""
    md_GOx = mdot_at_pressure(
        p_GOx_psia * psi_to_Pa, 'Oxygen', R_O2, A_inj_O2
    )
    md_GH2_target = md_GOx / o_f_target

    def residual(p_GH2_psia):
        md_GH2 = mdot_at_pressure(
            p_GH2_psia * psi_to_Pa, 'Hydrogen', R_H2, A_inj_H2
        )
        return md_GH2 - md_GH2_target

    return brentq(residual, 1, 10000)

p1_H2 = required_GH2_pressure(p1_O2 / psi_to_Pa) * psi_to_Pa

# =============================================================================
# MIXTURE COMPOSITION  (computed once, shared everywhere)
# =============================================================================
nO2    = Y_O2 / MW_O2
nH2    = Y_H2 / MW_H2
n_sum  = nO2 + nH2
xO2    = nO2 / n_sum
xH2    = nH2 / n_sum
MW_mix = xO2 * MW_O2 + xH2 * MW_H2
R_mix  = R_univ / MW_mix    # J/kg-K - mixture specific gas constant

# =============================================================================
# THERMODYNAMIC PROPERTIES  (computed at feed conditions)
# =============================================================================
cp_O2 = PropsSI('CPMASS', 'T', T, 'P', p1_O2, 'Oxygen')
cv_O2 = PropsSI('CVMASS', 'T', T, 'P', p1_O2, 'Oxygen')
cp_H2 = PropsSI('CPMASS', 'T', T, 'P', p1_H2, 'Hydrogen')
cv_H2 = PropsSI('CVMASS', 'T', T, 'P', p1_H2, 'Hydrogen')
g_O2  = cp_O2 / cv_O2
g_H2  = cp_H2 / cv_H2

# Critical pressure ratios (per species)
cf_O2 = (2 / (g_O2 + 1)) ** (g_O2 / (g_O2 - 1))
cf_H2 = (2 / (g_H2 + 1)) ** (g_H2 / (g_H2 - 1))

# Sonic-throat temperatures (not CoolProp TCRIT, which means critical point)
Tc_O2 = T * 2 / (g_O2 + 1)
Tc_H2 = T * 2 / (g_H2 + 1)

# Feed densities (for line velocity)
rho_O2 = PropsSI('DMASS', 'T', T, 'P', p1_O2, 'Oxygen')
rho_H2 = PropsSI('DMASS', 'T', T, 'P', p1_H2, 'Hydrogen')

# Throat conditions
p2_O2 = p1_O2 * cf_O2
p2_H2 = p1_H2 * cf_H2
P0_O2 = p2_O2
P0_H2 = p2_H2
# Original weighted-throat pressure approximation; not independently verified.
P0    = xO2 * P0_O2 + xH2 * P0_H2

rho2_O2 = PropsSI('DMASS', 'T', Tc_O2, 'P', p2_O2, 'Oxygen')
rho2_H2 = PropsSI('DMASS', 'T', Tc_H2, 'P', p2_H2, 'Hydrogen')

# Mixture Cp/Cv/gamma at initial blowdown pressure
Cp_O2_b      = PropsSI('CPMASS', 'T', T, 'P', P0_O2, 'Oxygen')
Cp_H2_b      = PropsSI('CPMASS', 'T', T, 'P', P0_H2, 'Hydrogen')
Cp_mix_molar = xO2 * Cp_O2_b * MW_O2 + xH2 * Cp_H2_b * MW_H2
Cp_mix       = Cp_mix_molar / MW_mix
Cv_mix       = Cp_mix - R_mix
gamma        = Cp_mix / Cv_mix

# Mixture viscosity via Wilke mixing rule
mu_O2 = PropsSI('VISCOSITY', 'T', T, 'P', P0_O2, 'Oxygen')
mu_H2 = PropsSI('VISCOSITY', 'T', T, 'P', P0_H2, 'Hydrogen')

def phi_wilke(mu_i, mu_j, MW_i, MW_j):
    return (1 + np.sqrt(mu_i / mu_j) * (MW_j / MW_i) ** 0.25) ** 2 / \
           np.sqrt(8 * (1 + MW_i / MW_j))

mu_mix = (xO2 * mu_O2 / (xO2 + xH2 * phi_wilke(mu_O2, mu_H2, MW_O2, MW_H2)) +
          xH2 * mu_H2 / (xH2 + xO2 * phi_wilke(mu_H2, mu_O2, MW_H2, MW_O2)))

# =============================================================================
# CALCULATIONS
# =============================================================================

# --- Supply pressures ---
md_O2      = choked_mdot(Cd, A_inj_O2, p1_O2, g_O2, R_O2, T)
md_H2      = choked_mdot(Cd, A_inj_H2, p1_H2, g_H2, R_H2, T)
v_choke_O2 = math.sqrt(g_O2 * R_O2 * Tc_O2)
v_choke_H2 = math.sqrt(g_H2 * R_H2 * Tc_H2)
v_line_O2  = md_O2 / (rho_O2 * A_line)
v_line_H2  = md_H2 / (rho_H2 * A_line)
q_O2       = md_O2 / rho2_O2
q_H2       = md_H2 / rho2_H2
q_tot      = q_O2 + q_H2
t_fill     = Vt / q_tot                 # original fill-time approximation
q_O2_feed  = md_O2 / rho_O2
q_H2_feed  = md_H2 / rho_H2
q_tot_feed = q_O2_feed + q_H2_feed

# --- Structural ---
BR = (OD_spiral**2 - ID_spiral**2) / OD_spiral**2

def lame_hoop(Pi, Po, ri, ro):
    return ((Pi * ri**2 - Po * ro**2) / (ro**2 - ri**2) +
            ri**2 * ro**2 * (Pi - Po) / (ri**2 * (ro**2 - ri**2)))

hoop1   = lame_hoop(P_i, P_o, r_i, r_o)
hoop2   = lame_hoop(P_i, P_o, r2_i, r2_o)
safety1 = allow_1_4 / hoop1
safety2 = allow_3_8 / hoop2

# --- Blowdown (original simplified isothermal model) ---
steps          = int(t_final / dt)
time           = np.zeros(steps)
P_tube         = np.zeros(steps)
P              = P0
m              = P * Vt / (R_mix * T)
lastIndex      = steps
unChokeTime    = np.nan
critical_ratio = (2 / (gamma + 1)) ** (gamma / (gamma - 1))

for i in range(steps):
    if P_atm / P <= critical_ratio:
        m_dot = (Cd * a1 * P * np.sqrt(gamma / (R_mix * T)) *
                 (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1))))
    else:
        if np.isnan(unChokeTime):
            unChokeTime = i * dt
        pr    = P_atm / P
        m_dot = (Cd * a1 * P * np.sqrt(
            (2 * gamma / (R_mix * T * (gamma - 1))) *
            (pr ** (2 / gamma) - pr ** ((gamma + 1) / gamma))
        ))
    m -= m_dot * dt
    m  = max(m, 0)
    P  = m * R_mix * T / Vt
    time[i]   = (i + 1) * dt
    P_tube[i] = P
    if P <= P_atm:
        lastIndex = i + 1
        break

time   = time[:lastIndex]
P_tube = P_tube[:lastIndex]

# =============================================================================
# OUTPUTS
# =============================================================================
print()
print("=" * 60)
print("  PRE-DETONATOR ANALYSIS — RESULTS")
print("=" * 60)

print("\n--- SUPPLY PRESSURES ---")
print(f"  GOx Feed Pressure    : {p1_O2 / psi_to_Pa:.3f} psia")
print(f"  GH2 Feed Pressure    : {p1_H2 / psi_to_Pa:.3f} psia")
print(f"  O2 Mass Flow Rate    : {md_O2 * 1000:.3f} g/s")
print(f"  H2 Mass Flow Rate    : {md_H2 * 1000:.3f} g/s")
print(f"  Mass Ratio (O2/H2)   : {md_O2 / md_H2:.4f}")
print(f"  O2 Choked Velocity   : {v_choke_O2 * m_to_ft:.1f} ft/s")
print(f"  O2 Line Velocity     : {v_line_O2 * m_to_ft:.1f} ft/s")
print(f"  H2 Choked Velocity   : {v_choke_H2 * m_to_ft:.1f} ft/s")
print(f"  H2 Line Velocity     : {v_line_H2 * m_to_ft:.1f} ft/s")
print(f"  Fill Time (model)    : {t_fill * 1000:.3f} ms")
print(f"  Chamber Pressure (*) : {P0 / psi_to_Pa:.3f} psia")

print("\n--- TESTING CONDITIONS ---")
print(f"  Total Volume         : {Vt * 1e6:.4f} cm^3")
print(f"  Blockage Ratio       : {BR:.4f}")
print(f"  Hoop Stress 1/4 in   : {hoop1:.2f} psi  (SF = {safety1:.2f})")
print(f"  Hoop Stress 3/8 in   : {hoop2:.2f} psi  (SF = {safety2:.2f})")

print("\n--- BLOWDOWN ---")
print(f"  Initial Pressure (*) : {P0 / psi_to_Pa:.2f} psia")
print(f"  Blowdown Duration    : {time[-1] * 1e3:.3f} ms")
if not np.isnan(unChokeTime):
    print(f"  Flow Unchokes At     : {unChokeTime * 1e3:.3f} ms")
print("\n(*) Weighted sonic-throat pressure approximation, not a validated chamber pressure.")

# =============================================================================
# PLOT
# =============================================================================
plt.figure(figsize=(8, 5))
plt.plot(time * 1e3, P_tube / psi_to_Pa, linewidth=2, label="Tube Pressure")
if not np.isnan(unChokeTime):
    plt.axvline(unChokeTime * 1e3, linestyle="--", linewidth=1.5,
                color="orange", label="Unchoked")
plt.axhline(14.7, linestyle="--", linewidth=1, color="gray", label="Atmospheric")
plt.xlabel("Time (ms)")
plt.ylabel("Pressure (psia)")
plt.title("Pre-Detonator Pressure Decay")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
