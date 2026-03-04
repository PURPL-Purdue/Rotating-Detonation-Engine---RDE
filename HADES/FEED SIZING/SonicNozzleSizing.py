import math

# Conversions
in_to_m = 0.0254
lb_to_kg = 0.453592

# Feed requirements
m_dot = 3 * lb_to_kg # kg/s
P_m = 2 * 10**6 # Pa
R = 287 # J/kg-K
gamma = 1.4
T = 283 # K
mass_flow_margin = 0.85
crit_ratio = (2 / (gamma+1))**(gamma/(gamma-1))
P_t_u_req = P_m / crit_ratio

# Solve for the feed area
wall_thickness = 0.095
feed_diameter = (16/16 - wall_thickness*2) * in_to_m
A_2 = math.pi * feed_diameter**2 / 4

# Solve for the downstream Mach number required
M_2 = (m_dot * mass_flow_margin) / (P_m * A_2) * math.sqrt(R * T / gamma)

# Solve for throat area
A_factor = 2/(gamma+1) * (1 + (gamma-1)/2 * M_2**2) ** (0.5 * (gamma+1)/(gamma-1))
A_t = A_2 / A_factor * M_2
D_t = 2 * math.sqrt(A_t/math.pi)

# Solve for total upstream pressure
T_t = T * (1 + (gamma-1)/2 * M_2**2)
A_1 = A_2
P_t_u = m_dot/A_t * math.sqrt(T_t * R / gamma) * ((gamma+1)/2)**((gamma+1)/(2*(gamma-1)))

# Output statements
print()
print("M_2 =", str(M_2))
print("A_2 =", str(A_2), "m^2")
print("A_t =", str(A_t), "m^2")
print("Throat Diameter:", str(D_t * 1000), "mm")
print("Upstream Pressure:", str(P_t_u * 10**(-5)), "bar")
print("Required Upstream Pressure:", str(P_t_u_req * 10**(-5)), "bar")

# Print if it chokes
if P_t_u > P_t_u_req:
    print("Flow is choked!")
else:
    print("Flow does not choke :(")

print()