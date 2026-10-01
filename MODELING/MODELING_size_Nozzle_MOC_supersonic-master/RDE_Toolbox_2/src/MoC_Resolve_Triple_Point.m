function [beta, delta, P_match, M3, M1_p] = MoC_Resolve_Triple_Point(iTripleParam)
% Inputs:
% P1     - Injection static pressure (assumed equal to bounding gas initial pressure) [Pa]
% P2     - Post-detonation static pressure [Pa]
% T2     - Post-detonation static temperature [K]
% Vcj    - Chapman-Jouguet velocity [m/s]
% gamma2 - Specific heat ratio of post-detonation gas
% R1     - Specific gas constant of bounding gas [J/kg-K]

% Read the inputs passed from main into this function's local workspace.
P1 = iTripleParam.P1;
P2 = iTripleParam.P2;
T2 = iTripleParam.T2;
Vcj = iTripleParam.Vcj;
gamma2 = iTripleParam.gamma2;
R1 = iTripleParam.R1;

% 1. Calculate bounding gas properties 
% The bounding gas consists of detonation products isentropically 
% expanded to the initial reactant pressure. (Fievisohn eqn. 8)
T1_prime = T2 * (P1 / P2)^((gamma2 - 1) / gamma2);

% Calculate the acoustic speed and Mach number of the bounding gas
a1_prime = sqrt(gamma2 * R1 * T1_prime);
M1_prime = Vcj / a1_prime;

% 2. Setup the nonlinear system for the triple point
% We must find M3 and theta_shock such that P3 = P2_prime and delta_3 = delta_2_prime.
% Flow direction and static pressures must be equal on both sides of the slip line

% Initial guesses
M3_guess = 1.5; 
theta_shock_guess = 45 * (pi/180); % Radians
x0 = [M3_guess, theta_shock_guess];

% fsolve options
options = optimoptions('fsolve', 'Display', 'none', 'FunctionTolerance', 1e-8, 'StepTolerance', 1e-8);

% Solve the system
[sol, fval, exitflag] = fsolve(@(x) triple_point_residuals(x, P1, P2, gamma2, M1_prime), x0, options);

if exitflag <= 0
    warning('Triple point solver did not converge.');
end

% 3. Extract outputs
M3 = sol(1);
beta_shock = sol(2); % In radians

% Calculate final matched slip-line angle (delta_3) (Sousa eqn. 3)
term1 = sqrt((gamma2 + 1) / (gamma2 - 1));
term2 = sqrt((gamma2 - 1) / (gamma2 + 1) * (M3^2 - 1));
delta_slip = term1 * atan(term2) - atan(sqrt(M3^2 - 1));

% Calculate final matched static pressure (P3) (Sousa eqn. 4)
P_matched = P2 * ( (1 + (gamma2 - 1)/2) / (1 + (gamma2 - 1)/2 * M3^2) )^(gamma2 / (gamma2 - 1));

% Convert angles to degrees for user readability if desired
beta_shock = beta_shock * (180/pi);
delta_slip = delta_slip * (180/pi);

% Return the names requested by main (angles in degrees).
beta = beta_shock;
delta = delta_slip;
P_match = P_matched;
M1_p = M1_prime;
end

function res = triple_point_residuals(x, P1, P2, gamma2, M1_prime)
M3 = x(1);
theta = x(2);

% Prevent unphysical values during iteration
if M3 <= 1 || theta <= 0 || theta >= pi/2
    res = [1e6, 1e6];
    return;
end

% Right Side: Prandtl-Meyer Expansion (assuming Mw2 = 1) (Sousa eqn. 3)
term1_exp = sqrt((gamma2 + 1) / (gamma2 - 1));
term2_exp = sqrt((gamma2 - 1) / (gamma2 + 1) * (M3^2 - 1));
delta_3 = term1_exp * atan(term2_exp) - atan(sqrt(M3^2 - 1));

% (Sousa eqn. 4)
P3 = P2 * ( (1 + (gamma2 - 1)/2) / (1 + (gamma2 - 1)/2 * M3^2) )^(gamma2 / (gamma2 - 1));

% Left Side: Oblique Shock (Sousa eqn. 2)
% Pressure ratio across the oblique shock
P2_prime = P1 * (1 + (2 * gamma2 / (gamma2 + 1)) * (M1_prime^2 * sin(theta)^2 - 1));

% Flow deflection angle across the oblique shock (Sousa eqn. 1)
num = M1_prime^2 * sin(theta)^2 - 1;
den = M1_prime^2 * (gamma2 + cos(2 * theta)) + 2;
delta_2_prime = atan(2 * cot(theta) * (num / den));

% Residuals to drive to zero[cite: 2]
res(1) = P3 - P2_prime;
res(2) = delta_3 - delta_2_prime;
end
