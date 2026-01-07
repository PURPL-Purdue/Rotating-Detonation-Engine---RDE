import math
import CoolProp.CoolProp as cp

#--------------------------------- FUNCTIONS ---------------------------------#
def chokedFlow(C_D, orificeSize, gas, pressure, temp):

    area = math.pi * (orificeSize/2 * in_to_m) ** 2

    gamma = cp.PropsSI("CPMOLAR", "T", temp, "P", pressure, gas) / cp.PropsSI("CVMOLAR", "T", temp, "P", pressure, gas)
    rho = cp.PropsSI("D", "T", temp, "P", pressure, gas)

    m_dot = C_D * area * math.sqrt(gamma * rho * pressure * (2/(gamma + 1))**((gamma+1)/(gamma-1)))
    Q_gas = m_dot / rho

    return m_dot, Q_gas

#--------------------------------- CONVERSIONS -------------------------------#
in_to_m = 0.0254
psi_to_Pa = 6894.76
kg_to_lbm = 2.20462

#--------------------------------- MAIN CODE ---------------------------------#

C_D = 0.6

# H2 Conditions
# Orifice Link: https://www.mcmaster.com/2275N39/
h2_orifice = 0.018 # in
h2_temp = 25 + 273.15 # K
h2_press = 200 * psi_to_Pa
[h2_m_dot, h2_Q] = chokedFlow(C_D, h2_orifice, "H2", h2_press, h2_temp) # kg/s, m^3/s

# O2 Conditions
# Orifice Link: https://www.mcmaster.com/2275N39/
o2_orifice = 0.035 # in
o2_temp = 25 + 273.15 # K
o2_press = 110 * psi_to_Pa
[o2_m_dot, o2_Q] = chokedFlow(C_D, o2_orifice, "O2", o2_press, o2_temp) # kg/s, m^3/s

# Pre-Det Dimensions
length = 6 * in_to_m # m
wall_thickness = 0.035 # in
manifold_volume = math.pi * ((0.25*in_to_m)/2)**2 * (2 * in_to_m) # m^3
PT_port_volume = math.pi * (((0.25-wall_thickness*2)*in_to_m)/2)**2 * (4 * in_to_m) # m^3
tube_volume = math.pi * (((0.25-wall_thickness*2)*in_to_m)/2)**2 * length # m^3
total_volume = manifold_volume + tube_volume + PT_port_volume # m^3

# Testing Values
total_Q = h2_Q + o2_Q
fill_time = total_volume / total_Q # s
# print("Total Vol. Flow Rate:", str(total_Q), "m^3/s")
print("H2 Mass Flow Rate:", str(h2_m_dot*1000), "g/s")
print("O2 Mass Flow Rate:", str(o2_m_dot*1000), "g/s")
print("O/F Ratio:", str(o2_m_dot/h2_m_dot))
# print("Total Volume:", str(total_volume), "m^3")
print("Total Fill Time:", str(fill_time*1000), "ms")
print("Pre-Det Chamber Pressure:", str(min([o2_press, h2_press])/psi_to_Pa/14.7), "bar")

# Shchelkin Spiral Calculations
OD = 0.3 # inches
ID = 0.21 # inches
BR = (OD**2 - ID**2) / (OD**2)
print("Blockage Ratio:", str(BR))

# Hoop Stress Calculations
P_i = 2531.25 # psia
P_o = 14.7 # psia
r_o = 0.25/2 # inches
r_i = r_o - wall_thickness # inches
hoop_stress = ((P_i * r_i**2 - P_o * r_o**2)/(r_o**2 - r_i**2)) + (r_i**2 * r_o**2 * (P_i - P_o) / (r_i**2 * (r_o**2 - r_i**2))) # psia
print("Hoop Stress:", str(hoop_stress), "psia")
P_i_safety = 4910.1 # psia
hoop_stress_safety = ((P_i_safety * r_i**2 - P_o * r_o**2)/(r_o**2 - r_i**2)) + (r_i**2 * r_o**2 * (P_i_safety - P_o) / (r_i**2 * (r_o**2 - r_i**2))) # psia
print("Pressure Safety Factor:", str(hoop_stress_safety/hoop_stress))

P_i = 2276.22 # psia
P_o = 14.7 # psia
r_i = 2 # in
r_o = 2.5 # in
chamber_stress = ((P_i * r_i**2 - P_o * r_o**2)/(r_o**2 - r_i**2)) + (r_i**2 * r_o**2 * (P_i - P_o) / (r_i**2 * (r_o**2 - r_i**2))) # psia
print(chamber_stress, "psia")