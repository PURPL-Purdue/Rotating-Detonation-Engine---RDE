%% Choked Flow Nozzle Analysis — Multiple Equivalence Ratios
% Formula: mdot = Cd * A * sqrt(gamma * rho0 * P0 * (2/(gamma+1))^((gamma+1)/(gamma-1)))

clear; clc; close all;

%% Given Constants
mdot_air  = 2 * 0.453592;       % Convert lbm/s to kg/s
T0        = 283;                 % Stagnation temperature [K]

% Stoichiometric fuel-to-air ratio for H2
FAR_stoich = 0.0292;             % kg H2 / kg Air (stoichiometric)

% Equivalence ratios to plot
phi_values = 0.8:0.1:1.2;        % [0.8, 0.9, 1.0, 1.1, 1.2]

% H2 mass flow rates for each equivalence ratio
% mdot_h2 = phi * FAR_stoich * mdot_air
mdot_h2_vec = phi_values * FAR_stoich * mdot_air;  % [kg/s]

% Throat areas: convert in^2 to m^2 (1 in^2 = 6.4516e-4 m^2)
A_air = 0.472   * 6.4516e-4;    % Air throat area [m^2]
A_h2  = 0.03875 * 6.4516e-4;    % H2 throat area  [m^2]

%% Gas Properties
% Air
gamma_air = 1.4;
R_air     = 287;                 % J/(kg·K)

% Hydrogen (H2)
gamma_h2  = 1.405;
R_h2      = 4157;                % J/(kg·K)

%% Cd Range
Cd = linspace(0.5, 0.9, 200);

%% Compute mass flow factors
K_air   = (2 / (gamma_air + 1))^((gamma_air + 1) / (gamma_air - 1));
MFF_air = sqrt(gamma_air * K_air / (R_air * T0));

K_h2    = (2 / (gamma_h2  + 1))^((gamma_h2  + 1) / (gamma_h2  - 1));
MFF_h2  = sqrt(gamma_h2  * K_h2  / (R_h2  * T0));

%% Solve for P0 — Air (single curve, mdot_air is fixed)
P0_air     = mdot_air ./ (Cd * A_air * MFF_air);
P0_air_bar = P0_air * 1e-5;

%% Color map for equivalence ratio curves
colors = lines(length(phi_values));

%% Plotting
figure('Position', [100, 100, 1200, 500]);

% --- Air Plot ---
subplot(1, 2, 1);
plot(Cd, P0_air_bar, 'b-', 'LineWidth', 2.5);
xlabel('Discharge Coefficient C_D',       'FontSize', 13);
ylabel('Stagnation Pressure P_0 [bar]',   'FontSize', 13);
title('Air: C_D vs. P_0',                 'FontSize', 14, 'FontWeight', 'bold');
subtitle(sprintf('\\dot{m} = %.3f kg/s,  T_0 = %d K,  A = %.6f m^2', mdot_air, T0, A_air));
grid on; grid minor;
set(gca, 'FontSize', 11);

% --- H2 Plot ---
subplot(1, 2, 2);
hold on;
for i = 1:length(phi_values)
    P0_h2     = mdot_h2_vec(i) ./ (Cd * A_h2 * MFF_h2);
    P0_h2_bar = P0_h2 * 1e-5;
    plot(Cd, P0_h2_bar, 'Color', colors(i,:), 'LineWidth', 2.5, ...
         'DisplayName', sprintf('\\phi = %.1f  (\\dot{m} = %.4f kg/s)', ...
         phi_values(i), mdot_h2_vec(i)));
end
hold off;
xlabel('Discharge Coefficient C_D',       'FontSize', 13);
ylabel('Stagnation Pressure P_0 [bar]',   'FontSize', 13);
title('H_2: C_D vs. P_0',                 'FontSize', 14, 'FontWeight', 'bold');
subtitle(sprintf('T_0 = %d K,  A = %.8f m^2,  FAR_{stoich} = %.4f', T0, A_h2, FAR_stoich));
legend('Location', 'best', 'FontSize', 10);
grid on; grid minor;
set(gca, 'FontSize', 11);

sgtitle('Choked Nozzle Flow: C_D vs. Stagnation Pressure', 'FontSize', 15, 'FontWeight', 'bold');