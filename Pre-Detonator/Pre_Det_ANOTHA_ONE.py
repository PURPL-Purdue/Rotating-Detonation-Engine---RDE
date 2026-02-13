print("##### PURPL RDE Pre-Det Propellant Orifice Sizing #####\n")
# Import stuff
import math
import CoolProp.CoolProp as CoolProp
from CoolProp.CoolProp import PropsSI

## Define Parameters
m_to_in = 39.3701; # in/m
m2_to_in2 = 1550; # in^2/m^2
kgs_to_lbm = 2.20462; # lbm/kg
psia_to_Pa = 6894.76; # Pa/psia
in_to_ft = 1/12; # ft/in
g = 9.81; #m/s^2
lbs_to_slug = 0.031081; #slug*ft/s^2/lb
bar_to_Pa = 100000; #Pa/bar
estimated_discharge_coefficient_thick_plate = 0.86; #unitless; assumed based on thick plate, square edge
#actual discharge coefficient will be determined experimentally - c_d_actual = mdot_measured/mdot_calculated

### PRE-DETONATOR PARAMETERS FROM DEEPESH
#min_chamber_pressure = 14.7*psia_to_Pa; #14.7psia to Pa
#print(f"Minimum Manifold Chamber Pressure: {round(min_chamber_pressure*(1/psia_to_Pa), 5)} psia")
#max_chamber_pressure = 40*psia_to_Pa; #40psia to Pa
#print(f"Maximum Manifold Chamber Pressure: {round(max_chamber_pressure*(1/psia_to_Pa), 5)} psia")
pre_det_manifold_pressure = 14.30*psia_to_Pa; #14.30 psia to Pa (pressure at Purdue elevation of 750ft above sea level)

#GOx Parameters
print("\n----GOx Calculations----")
num_orifices_GOx = 1; #unitless
print(f"Number of GOx Orifices: {num_orifices_GOx}")
m_dot_GOx = .274/1000; #.274 g/s to kg/s (MAX MASS FLOW OF GOx from CEA)
print(f"Total GOx Mass Flowrate: {round(m_dot_GOx*kgs_to_lbm, 5)} lbm/s")
inlet_temp_GOx_cold = 20 + 273.15; #K (10 C -- 50 degF)
print(f"GOx Inlet Temperature: {inlet_temp_GOx_cold} K")
min_GOx_feed_pressure = 13.05*1.893*(pre_det_manifold_pressure); #Pa
print(f"Minimum GOx Feed Pressure: {round(min_GOx_feed_pressure*(1/psia_to_Pa), 5)} psia")
density_GOx = PropsSI('D', 'T', inlet_temp_GOx_cold, 'P', min_GOx_feed_pressure, 'Oxygen') #kg/m^3
print(f"GOx Density at {inlet_temp_GOx_cold} K and {round(min_GOx_feed_pressure*(1/psia_to_Pa), 5)} psia: {round(density_GOx, 2)} kg/m^3")
dynamic_viscosity_GOx = PropsSI('V', 'T', inlet_temp_GOx_cold, 'P', min_GOx_feed_pressure, 'Oxygen') #Pa*s
print(f"GOx Dynamic Viscosity at {inlet_temp_GOx_cold} K and {round(min_GOx_feed_pressure*(1/psia_to_Pa), 5)} psia: {round(dynamic_viscosity_GOx, 5)} Pa*s")
cp_GOx = PropsSI('CPMASS', 'T', inlet_temp_GOx_cold, 'P', min_GOx_feed_pressure, 'Oxygen') #J/kg-K
cv_GOx = PropsSI('CVMASS', 'T', inlet_temp_GOx_cold, 'P', min_GOx_feed_pressure, 'Oxygen') #J/kg-K
gamma_GOx = cp_GOx/cv_GOx; #unitless
print(f"GOx Specific Heat Ratio at {inlet_temp_GOx_cold} K and {round(min_GOx_feed_pressure*(1/psia_to_Pa), 5)} psia: {round(gamma_GOx, 5)}")

#GH2 Parameters
print("\n----GH2 Calculations----")
num_orifices_GH2 = 1; #unitless
print(f"Number of GH2 Orifices: {num_orifices_GH2}")
m_dot_GH2 = .034/1000; #.034g/s to kg/s (MAX MASS FLOW OF GH2 from CEA)
print(f"Total GH2 Mass Flowrate: {round(m_dot_GH2*kgs_to_lbm, 5)} lbm/s")
inlet_temp_GH2_cold = 20 + 273.15; #K (10 C -- 50 degF)
print(f"GH2 Inlet Temperature: {inlet_temp_GH2_cold} K")
min_GH2_feed_pressure = 6.58*1.899*(pre_det_manifold_pressure); #Pa
print(f"Minimum GH2 Feed Pressure: {round(min_GH2_feed_pressure*(1/psia_to_Pa), 5)} psia")
density_GH2 = PropsSI('D', 'T', inlet_temp_GH2_cold, 'P', min_GH2_feed_pressure, 'Hydrogen') #kg/m^3
print(f"GH2 Density at {inlet_temp_GH2_cold} K and {round(min_GH2_feed_pressure*(1/psia_to_Pa), 5)} psia: {round(density_GH2, 5)} kg/m^3")
dynamic_viscosity_GH2 = PropsSI('V', 'T', inlet_temp_GH2_cold, 'P', min_GH2_feed_pressure, 'Hydrogen') #Pa*s
print(f"GH2 Dynamic Viscosity at {inlet_temp_GH2_cold} K and {round(min_GH2_feed_pressure*(1/psia_to_Pa), 5)} psia: {round(dynamic_viscosity_GH2, 7)} Pa*s")
cp_GH2 = PropsSI('CPMASS', 'T', inlet_temp_GH2_cold, 'P', min_GH2_feed_pressure, 'Hydrogen') #J/kg-K
cv_GH2 = PropsSI('CVMASS', 'T', inlet_temp_GH2_cold, 'P', min_GH2_feed_pressure, 'Hydrogen') #J/kg-K
gamma_GH2 = cp_GH2/cv_GH2; #unitless
print(f"GH2 Specific Heat Ratio at {inlet_temp_GH2_cold} K and {round(min_GH2_feed_pressure*(1/psia_to_Pa), 5)} psia: {round(gamma_GH2, 5)}")


### C_d*A Method Calculations with Square Shoulder Orifice
print("\n----C_d*A Method Calculations----")
print(f"Assuming C = Cd & Estimated Discharge Coefficient of {estimated_discharge_coefficient_thick_plate} for a Thick Plate Orifice")
# GOx Orifice Diameter for Choked Flow Equation
GOx_CdA_individual_orifice_C_equals_Cd = (m_dot_GOx/num_orifices_GOx)*math.sqrt(1/(gamma_GOx*density_GOx*min_GOx_feed_pressure*((2/(gamma_GOx+1))**((gamma_GOx+1)/(gamma_GOx-1))))); #m^2
GOx_individual_orifice_area_C_equals_Cd = GOx_CdA_individual_orifice_C_equals_Cd/estimated_discharge_coefficient_thick_plate; #m^2
print(f"GOx Individual Orifice Area (C = Cd): {round(GOx_individual_orifice_area_C_equals_Cd*m2_to_in2, 5)} in^2")
GOx_individual_orifice_dia_C_equals_Cd = math.sqrt((4*GOx_individual_orifice_area_C_equals_Cd)/math.pi); #m
print(f"GOx Individual Orifice Diameter (C = Cd): {round(GOx_individual_orifice_dia_C_equals_Cd*m_to_in, 5)} in")

# GH2 Orifice Diameter for Choked Flow Equation
GH2_CdA_individual_orifice_C_equals_Cd = (m_dot_GH2/num_orifices_GH2)*math.sqrt(1/(gamma_GH2*density_GH2*min_GH2_feed_pressure*((2/(gamma_GH2+1))**((gamma_GH2+1)/(gamma_GH2-1))))); #m^2
GH2_individual_orifice_area_C_equals_Cd = GH2_CdA_individual_orifice_C_equals_Cd/estimated_discharge_coefficient_thick_plate; #m^2
print(f"GH2 Individual Orifice Area (C = Cd): {round(GH2_individual_orifice_area_C_equals_Cd*m2_to_in2, 5)} in^2")
GH2_individual_orifice_dia_C_equals_Cd = math.sqrt((4*GH2_individual_orifice_area_C_equals_Cd)/math.pi); #m
print(f"GH2 Individual Orifice Diameter (C = Cd): {round(GH2_individual_orifice_dia_C_equals_Cd*m_to_in, 5)} in")


# GOx Line Velocity Check
GOx_orifice_velocity = m_dot_GOx/(density_GOx*GOx_individual_orifice_area_C_equals_Cd); #m/s
print(f"\nGOx Line Velocity through Orifice: {round(GOx_orifice_velocity*m_to_in**in_to_ft, 2)} ft/s (Note: must be below 80-100 ft/s)")