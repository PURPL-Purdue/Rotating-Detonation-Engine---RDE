# Global
max_phi = 1.2
stoic_OF = 34.3
T_amb = 293 # K
T_amb_F = 67.73 # degF

# Air
m_dot_air = 0.907185 # kg/s
p_air = 5.171e+6 # Pa
R_air = 287 # J/kg-K
rho_air = p_air / (R_air * T_amb) # kg/m^3
Q_air = m_dot_air / rho_air # m^3/s
Q_air_acfm = Q_air * 2118.8799727597 # acfm
Q_air_scfm = Q_air_acfm * ((p_air * 0.000145038)/14.7) * (519/(460 + (T_amb_F)))

# H2
m_dot_h2 = m_dot_air / stoic_OF * max_phi # kg/s
p_h2 = 5.171e+6 # Pa
R_h2 = 4126 # J/kg-K
rho_h2 = p_h2 / (R_h2 * T_amb) # kg/m^3
Q_h2 = m_dot_h2 / rho_h2 # m^3/s
Q_h2_acfm = Q_h2 * 2118.8799727597 # acfm
Q_h2_scfm = Q_h2_acfm * ((p_h2 * 0.000145038)/14.7) * (519/(460 + (T_amb_F)))

# Output
print("Air:", str(Q_air_scfm), "scfm")
print("H2:", str(Q_h2_scfm), "scfm")