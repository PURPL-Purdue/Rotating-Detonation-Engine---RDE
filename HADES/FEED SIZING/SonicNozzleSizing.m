% Conversions
in_to_m = 0.0254;
lb_to_kg = 0.453592;

% Feed requirements
m_dot = 3 * lb_to_kg; % kg/s
P_m = 2 * 10^6; % Pa
R = 287; % J/kg-K
gamma = 1.4;
T = 283; % K
mass_flow_margin = 0.85;
crit_ratio = (2 / (gamma+1))^(gamma/(gamma-1));
P_t_u_req = P_m / crit_ratio;

% Solve for the feed area
wall_thickness = 0.095;
feed_diameter = (16/16 - wall_thickness*2) * in_to_m;
A_2 = pi * feed_diameter^2 / 4;

% Solve for the downstream Mach number required
M_2 = (m_dot * mass_flow_margin) / (P_m * A_2) * sqrt(R * T / gamma);

% Solve for throat area
A_factor = 2/(gamma+1) * (1 + (gamma-1)/2 * M_2^2) ^ (0.5 * (gamma+1)/(gamma-1));
A_t = A_2 / A_factor * M_2;
D_t = 2 * sqrt(A_t/pi);

% Solve for total upstream pressure
T_t = T * (1 + (gamma-1)/2 * M_2^2);
A_1 = A_2;
P_t_u = m_dot/A_t * sqrt(T_t * R / gamma) * ((gamma+1)/2)^((gamma+1)/(2*(gamma-1)));