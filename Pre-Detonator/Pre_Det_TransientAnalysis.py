import numpy as np
from math import pi
import matplotlib.pyplot as plt
import sys

#DEFINING A BUNCH OF CONSTANTS

# conversions
in_m = 0.0254 # in -> m
psi_pa = 6894.76 # psi -> Pa
mpa_to_psi = 145.038 # MPa -> psi conversion factor for outputs

# conditions

P_inital = 110 #[psig]
P_standard = 14.7 # 1 atm
T_inital_gas = 293.15 # 68 F
T_ambient_wall = 298.0 #25 C 

# detparameters (NASA CEA)

det_width = 0.0000010 # source
CJ = 2952.3 #[m/s]
flame_vel = 0.8 * CJ #[m/s]
T_ratio = 14.010 #[N/A]
P_ratio = 20.537 #[N/A]
T_gas = 4106.92 #[K]

# geometry
tube_length = 6 * in_m #[m]
tube_thickness = 0.035 * in_m #[m]
D2 = 0.25 * in_m #[m]
D1 = D2 - 2 * tube_thickness #[m]
mean_diameter = (D1 + D2) / 2

# ss316 properties: https://www.azom.com/properties.aspx?ArticleID=863 using lower bounds
k_wall = 13.0 #thermal conductivity [W / (m*K)]
rho_wall = 7870 # density [Kg/m^3]
c_wall = 490 # specific heat [J / (Kg*K)]

# gas properties (FROM NASA CEA --> converted into units for equations)
k_gas = 2.814 # thermal conductivity (w/ (m * K))
mu_gas = 0.00012366 # dynamic viscosity (kg / (m * s))
cp_gas = 11656.9 # specific heat (J / (kg * K))
rho_gas = 7.7717 # density (kg/m^3)


h_ambient = 10 # convection outside the tube [W / (m^2 * K)]


#yeild strength model

temp_points_K = [
    300.15,  # 27 C
    422.15,  # 149 C
    533.15,  # 260 C
    700.15,  # 371 C
    755.15,  # 482 C
    866.15,  # 593 C
    977.15,  # 704 C
    1088.15, # 816 C
    1200.15, # 927 C
    1311.15, # 1038 C
    1366.15,  # 1093 C
    1648.00  # melting 
]

yeild_points_MPa = [
    290,
    201,
    172,
    159,
    148,
    140,
    131,
    110,
    80,  # From 1700F, MPa=80
    39,  # From 1900F, MPa=39
    28,   # From 2000F, MPa=28
    1.0   # Melting point strength
]

def get_yield_strength(T_kelvin):
    strength = np.interp(T_kelvin, temp_points_K, yeild_points_MPa)
    return strength

# simulation model

def run_simulation(t_total, num_layers=1000, capture_history=False):
    # Set up Mesh 
    radii = np.linspace(D1 / 2, D2 / 2, num_layers + 1)
    dr = radii[1] - radii[0]
    
    # Automatic Stability Check (CFL Condition)
    alpha = k_wall / (rho_wall * c_wall) 
    dt_critical = 0.5 * (dr**2) / alpha  
    
    required_steps = int(t_total / dt_critical) + 100 
    num_steps = max(required_steps, 5000) 
    dt = t_total / num_steps

    #calculate h (convective coeeficient)
    # Dittus-Boelter
    Re = (rho_gas * flame_vel * D1) / mu_gas 
    Pr = (cp_gas * mu_gas) / k_gas              
    Nu = 0.023 * (Re ** 0.8) * (Pr ** 0.4)      
    h_inner = (Nu * k_gas) / D1

    A_interfaces = 2 * pi * radii * tube_length
    mass_layers = pi * (radii[1:]**2 - radii[:-1]**2) * tube_length * rho_wall

    #LOOPING
    T_profile = np.full(num_layers, T_ambient_wall)
    T_new = np.copy(T_profile)

    #history for single pulse
    single_pulse_temp = []
    single_pulse_time = []

    for t_step in range(num_steps):
        #Inner Wall (Layer 0)
        #convection from gas: q = h * Area * (t_gas - t_wall)
        q_conv = h_inner * (pi * D1 * tube_length) * (T_gas - T_profile[0])
        #conduction to next layer: fourier's law: q_out = k * Area * (dT / dr)
        q_cond_out = k_wall * A_interfaces[1] * (T_profile[0] - T_profile[1]) / dr
        
        # temp left over stays in new temperature array
        T_new[0] += (q_conv - q_cond_out) * dt / (mass_layers[0] * c_wall)

        # Middle Layers (Vectorized for 50x Speedup)
        # We calculate flux for layers 1 to N-2
        area_in = A_interfaces[1 : num_layers-1]
        area_out = A_interfaces[2 : num_layers]
        
        #conduction in
        q_in = k_wall * area_in * (T_profile[:-2] - T_profile[1:-1]) / dr
        #conduction out
        q_out = k_wall * area_out * (T_profile[1:-1] - T_profile[2:]) / dr
        
        # update new temp array
        T_new[1:-1] += (q_in - q_out) * dt / (mass_layers[1:-1] * c_wall)

        #Outer Wall (Layer num_layers-1)
        # conduction in
        q_cond_in = k_wall * A_interfaces[-2] * (T_profile[-2] - T_profile[-1]) / dr
        # convection out to ambient air
        q_conv_out = h_ambient * (pi * D2 * tube_length) * (T_profile[-1] - T_ambient_wall)
        
        # update new temp array
        T_new[-1] += (q_cond_in - q_conv_out) * dt / (mass_layers[-1] * c_wall)

        T_profile[:] = T_new[:]

        if capture_history and t_step % 100 == 0:
            current_time = t_step * dt
            single_pulse_time.append(current_time * 1000)
            single_pulse_temp.append(T_profile[0])
    
    # Return the final temperature profile
    if capture_history:
        return T_profile, h_inner, single_pulse_time, single_pulse_temp
    else:
        return T_profile, h_inner, [], []


# 4. MAIN ITERATION LOOP

print("--- Starting Failure Analysis ---")
print("Note: Calculation performed in Metric, Output converted to Imperial.")

# Set up failure parameters
safety_factor = 1.0
pulse_coefficient = 1 # pulse number
max_pulses = 100 # how many pulses to check

# Hoop Stress
# constant so calculate it once

P_initial_abs_psi = P_inital + P_standard
P_peak_Pa = (P_initial_abs_psi * psi_pa) * P_ratio # P_internal (Pi)
P_atm_Pa = P_standard * psi_pa                     # P_external (Po)

# Using the thick walled hoop stress calculation

r_i = D1 / 2  # Inner Radius
r_o = D2 / 2  # Outer Radius
r = r_i
# Term 1: (Pi*ri^2 - Po*ro^2) / (ro^2 - ri^2)

term1 = (P_peak_Pa * r_i**2 - P_atm_Pa * r_o**2) / (r_o**2 - r_i**2)

# Term 2: (ri^2 * ro^2 * (Pi - Po)) / (r^2 * (ro^2 - ri^2))
term2 = (r_i**2 * r_o**2 * (P_peak_Pa - P_atm_Pa)) / (r**2 * (r_o**2 - r_i**2))

# Total Hoop Stress
hoop_stress_Pa = term1 + term2
hoop_stress_MPa = hoop_stress_Pa / 1e6
hoop_stress_psi = hoop_stress_MPa * mpa_to_psi

print(f"Constant Hoop Stress (Load): {hoop_stress_psi:.0f} psi ({hoop_stress_MPa:.2f} MPa)")


# THE loop

# lists for plotting
list_pulse_count = []
list_peak_temp_F = []
list_yield_psi = []
list_hoop_psi = []

# get data for one pulse
t_one_pulse = 1 * (det_width / flame_vel)
_, _, hist_time, hist_temp_K = run_simulation(t_one_pulse, num_layers=1000, capture_history=True)

# Convert Single Pulse History to F
hist_temp_F = [(t - 273.15) * 9/5 + 32 for t in hist_temp_K]

pulse_coefficient = 1
while pulse_coefficient <= max_pulses:
    t_tot = pulse_coefficient * (det_width / flame_vel)
    
    # Run Sim
    T_final, _, _, _ = run_simulation(t_tot, num_layers=1000, capture_history=False)
    
    # Check Failure
    peak_T_K = T_final[0]
    strength_MPa = get_yield_strength(peak_T_K)
    limit_MPa = strength_MPa / safety_factor
    
    # Convert for Output
    peak_T_F = (peak_T_K - 273.15) * 9/5 + 32
    strength_psi = strength_MPa * mpa_to_psi
    limit_psi = limit_MPa * mpa_to_psi
    
    # Store Data
    list_pulse_count.append(pulse_coefficient)
    list_peak_temp_F.append(peak_T_F)
    list_yield_psi.append(limit_psi)
    list_hoop_psi.append(hoop_stress_psi)
    
    print(f"Pulse {pulse_coefficient}: T={peak_T_F:.0f} F | Stress={hoop_stress_psi:.0f} psi vs Limit={limit_psi:.0f} psi")
    
    if hoop_stress_MPa > limit_MPa:
        print(f"\n!!! FAILURE at Pulse {pulse_coefficient} !!!")
        break
    pulse_coefficient += 1

if pulse_coefficient > max_pulses:
    print(f"\n--- Analysis Stopped ---")
    print(f"Reached max iterations ({max_pulses}) without failure.")


# Create a figure with two distinct graphs (subplots)
fig = plt.figure(figsize=(12, 10))
gs = fig.add_gridspec(2, 2) # Create a 2x2 grid

# Graph 1: Single Pulse Transient (Top Left) 
ax3 = fig.add_subplot(gs[0, 0])
ax3.plot(hist_time, hist_temp_F, 'k-o', linewidth=1.5, markersize=4, markevery=100)
ax3.set_title('Single Pulse Transient Rise')
ax3.set_xlabel('Time (ms)')
ax3.set_ylabel('Inner Wall Temp (F)')
ax3.grid(True)

# Graph 2: Accumulation vs Pulse Number (Top Right) 
ax1 = fig.add_subplot(gs[0, 1])
ax1.plot(list_pulse_count, list_peak_temp_F, 'b-o')
ax1.set_title('Temp Accumulation vs. Detonation Number')
ax1.set_xlabel('Detonation Number')
ax1.set_ylabel('Peak Wall Temp (F)')
ax1.grid(True)

# Graph 3: Failure Analysis (Bottom - Full Width)
ax2 = fig.add_subplot(gs[1, :])
ax2.plot(list_pulse_count, list_yield_psi, 'g-o', label='Material Strength Limit (psi)')
ax2.plot(list_pulse_count, list_hoop_psi, 'r-', linewidth=2, label='Hoop Stress Load (psi)')
ax2.set_title('Failure Analysis: Load vs. Strength')
ax2.set_xlabel('Detonation Number')
ax2.set_ylabel('Stress (psi)')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()

# --- Graph 4: Temperature Gradient (The "Heat Soak" View) ---
# 1. Re-calculate radii (because the variable was local to the function)
num_layers = 1000
radii = np.linspace(D1 / 2, D2 / 2, num_layers + 1)

# 2. Convert to mm for readability
radii_in = radii / 0.0254 # Meters -> Inches

# 3. Plot
plt.figure(figsize=(8, 6))

# We use T_final from the last run of the loop (Pulse 11 or 12)
layer_centers_in = (radii_in[:-1] + radii_in[1:]) / 2 

# Convert final profile to F
T_final_F = [(t - 273.15) * 9/5 + 32 for t in T_final]

plt.plot(layer_centers_in, T_final_F, 'r-', linewidth=2)

plt.title(f'Temperature Gradient Through Wall at Failure (Pulse {len(list_pulse_count)})')
plt.xlabel('Distance from Center (in)')
plt.ylabel('Temperature (F)')
plt.grid(True)

# Add dashed lines for the wall boundaries
plt.axvline(x=radii_in[0], color='k', linestyle='--', label='Inner Surface')
plt.axvline(x=radii_in[-1], color='k', linestyle='--', label='Outer Surface')
plt.legend()

plt.show()

print("\n--- Analysis Complete ---")