import math
import CoolProp.CoolProp as cp

#--------------------------------- FUNCTIONS ---------------------------------#
def chokedFlow(C_D, orificeSize, gas, pressure, temp):

    area = math.pi * (orificeSize/2 * in_to_m) ** 2

    gamma = cp.PropsSI("CPMOLAR", "T", temp, "P", pressure, gas) / cp.PropsSI("CVMOLAR", "T", temp, "P", pressure, gas)
    rho = cp.PropsSI("D", "T", temp, "P", pressure, gas)

    m_dot = C_D * area * math.sqrt(gamma * rho * pressure * (2/(gamma + 1))**((gamma+1)/(gamma-1)))
    Q_gas = m_dot / rho

    M = cp.PropsSI("M", gas)  # molar mass [kg/mol]
    R = 8.314 / M # J/kg-K

    return m_dot, Q_gas, R

#--------------------------------- CONVERSIONS -------------------------------#
in_to_m = 0.0254
psi_to_Pa = 6894.76
kg_to_lbm = 2.20462
choke_factor = 0.528

#--------------------------------- MAIN CODE ---------------------------------#

C_D = 0.6

# H2 Conditions
# Orifice Link: https://www.mcmaster.com/2275N39/
h2_orifice = 0.018 # in
h2_temp = 10 + 273.15 # K
h2_press = 55.125 * psi_to_Pa
h2_dpress = h2_press * choke_factor / psi_to_Pa # psia
[h2_m_dot, h2_Q, h2_R] = chokedFlow(C_D, h2_orifice, "H2", h2_press, h2_temp) # kg/s, m^3/s

# O2 Conditions
# Orifice Link: https://www.mcmaster.com/2275N39/
o2_orifice = 0.035 # in
o2_temp = 10 + 273.15 # K
o2_press = 29.4 * psi_to_Pa
o2_dpress = o2_press * choke_factor / psi_to_Pa # psia
[o2_m_dot, o2_Q, o2_R] = chokedFlow(C_D, o2_orifice, "O2", o2_press, o2_temp) # kg/s, m^3/s

# Pre-Det Dimensions
quarter_length = (4.3955 + 0.74) * in_to_m # m
threeeight_length = 1.5 * in_to_m # m
wall_thickness = 0.035 # in
spring_volume = 0.0216 * in_to_m**3 # m^3
manifold_volume = math.pi * ((0.3873*in_to_m)/2)**2 * (1.19 * in_to_m) # m^3
tube_volume = math.pi * (((0.25-wall_thickness*2)*in_to_m)/2)**2 * quarter_length + math.pi * ((0.305 * in_to_m)/2)**2 * threeeight_length - spring_volume # m^3
total_volume = manifold_volume + tube_volume # m^3

# Testing Values
total_Q = h2_Q + o2_Q
total_mdot = h2_m_dot + o2_m_dot
R_mix = 8/9 * o2_R + 1/9 * h2_R # J/kg-K
fill_time = total_volume / total_Q # s
# print("Total Vol. Flow Rate:", str(total_Q), "m^3/s")
print("H2 Mass Flow Rate:", str(h2_m_dot*1000), "g/s")
print("H2 Feed Pressure: ", str(h2_dpress), "psia")
print("O2 Mass Flow Rate:", str(o2_m_dot*1000), "g/s")
print("H2 Feed Pressure: ", str(o2_dpress), "psia")
print("O/F Ratio:", str(o2_m_dot/h2_m_dot))
# print("Total Volume:", str(total_volume), "m^3")
print("Total Fill Time:", str(fill_time*1000), "ms")
print("Pre-Det Chamber Pressure:", str((total_mdot * R_mix * h2_temp) / total_Q / psi_to_Pa), "psia")

# Shchelkin Spiral Calculations
OD = 0.3 # inches
ID = 0.21 # inches
BR = (OD**2 - ID**2) / (OD**2)
print("Blockage Ratio:", str(BR))

# Hoop Stress Calculations
P_i = 80 * 14.7 # psia
P_o = 14.7 # psia
r_o = 0.25/2 # inches
r2_o = 0.375/2 # inches
r2_i = r2_o - wall_thickness # inches
r_i = r_o - wall_thickness # inches
hoop_stress = ((P_i * r_i**2 - P_o * r_o**2)/(r_o**2 - r_i**2)) + (r_i**2 * r_o**2 * (P_i - P_o) / (r_i**2 * (r_o**2 - r_i**2))) # psia
hoop_stress2 = ((P_i * r2_i**2 - P_o * r2_o**2)/(r2_o**2 - r2_i**2)) + (r2_i**2 * r2_o**2 * (P_i - P_o) / (r2_i**2 * (r2_o**2 - r2_i**2))) # psia
print("Hoop Stress 1/4 in:", str(hoop_stress), "psia")
print("Hoop Stress 3/8 in:", str(hoop_stress2), "psia")
hoop_stress_safety = 5100 # psia
hoop_stress_safety2 = 3300 # psia
print("Pressure Safety Factor 1/4 in:", str(hoop_stress_safety/hoop_stress))
print("Pressure Safety Factor 3/8 in:", str(hoop_stress_safety2/hoop_stress2))