% =========================================================
%  Compressible Subcritical Orifice Flow Calculator
% =========================================================
%% WARNING - This is an empirical formulation, as D1 --> D2 error explodes. 

clear; clc;

% ---------------------------------------------------------
%  USER INPUTS
% ---------------------------------------------------------
D1      = 0.01125;   % Upstream pipe diameter (m)
D2      = 0.00225;       % Orifice bore diameter (m)
P1      = 19e5;       % Upstream static pressure (Pa)
P2      = 13e5;       % Downstream static pressure (Pa)
T1      = 283.15;           % Upstream temperature (K)
gamma   = 1.4;          % Ratio of specific heats
Rg      = 259.8;          % Gas constant (J/kg-K)
mu      = 2.35e-5;     % Dynamic viscosity (Pa-s)

Cd_guess = 0.61;        % Initial Cd guess (look at slide 17)

tol     = 1e-6;
maxiter = 100;

% ---------------------------------------------------------
% Geometry
% ---------------------------------------------------------
A2   = pi/4 * D2^2;
A1   = pi/4 * D1^2;
beta = D2 / D1;

D1_in = D1 / 0.0254;

% ---------------------------------------------------------
% Check choking
% ---------------------------------------------------------
PR_critical = ((gamma + 1) / 2)^(gamma / (gamma - 1));
PR_actual   = P1 / P2;

if PR_actual >= PR_critical
    warning('Flow is CHOKED - this model is invalid.');
end

% ---------------------------------------------------------
% Stagnation pressure
% ---------------------------------------------------------
P0 = (( (D1/D2)^4 * P1^((gamma+1)/gamma) - P2^((gamma+1)/gamma) ) / ...
      ( (D1/D2)^4 * P1^(2/gamma) - P2^(2/gamma) ))^(gamma/(gamma-1));


rho0 = P0 / (Rg * T1);

% ---------------------------------------------------------
% Iterate Cd (corner tap)
% ---------------------------------------------------------
Cd_inc = Cd_guess;

for iter = 1:maxiter

    mdot = Cd_inc * A2 * sqrt(2 * (P1/(Rg*T1)) * (P1 - P2));

    Re_D1 = 4 * mdot / (pi * D1 * mu);

    Cv = ( 0.5991 + 0.0044/D1_in ...
         + (0.3155 + 0.0175/D1_in)*(beta^4 + 2*beta^16) ) ...
         * sqrt(1 - beta^4) + ( 0.52/D1_in - 0.192 ...
         + (16.48 - 1.16/D1_in)*(beta^4 + 4*beta^16) ) ...
         * sqrt((1 - beta^4)/Re_D1);

    Cd_new = Cv / sqrt(1 - beta^4);

    if abs(Cd_new - Cd_inc)/Cd_inc < tol
        Cd_inc = Cd_new;
        break
    end

    Cd_inc = Cd_new;
end

% ---------------------------------------------------------
% Loss coefficient
% ---------------------------------------------------------
f = 1/Cd_inc - 1/(2*Cd_inc^2);

% ---------------------------------------------------------
% Flow coefficient
% ---------------------------------------------------------
r  = P2 / P0;

Kn = sqrt((2*gamma/(gamma-1)) * r^(2/gamma) * (1 - r^((gamma-1)/gamma)));

% ---------------------------------------------------------
% Compressible Cd
% ---------------------------------------------------------
term = 1 - f * (2*r^(1/gamma))^2 * (1 - r) / Kn^2;

if term < 0
    error('Negative discriminant - check inputs.');
end

Cd_comp = (1 - sqrt(term)) / (2 * f * r^(1/gamma));

% ---------------------------------------------------------
% Mass flow
% ---------------------------------------------------------
mdot_comp = Kn * Cd_comp * A2 * sqrt(P0 * rho0);

% ---------------------------------------------------------
%  OUTPUT
% ---------------------------------------------------------
fprintf('\n--- Results ---\n');
fprintf('Cd (incompressible) = %.6f\n', Cd_inc);
fprintf('Cd (compressible)   = %.6f\n', Cd_comp);
fprintf('Mass flow rate      = %.6f kg/s\n', mdot_comp);
fprintf('Re_D1               = %.0f\n', Re_D1);
