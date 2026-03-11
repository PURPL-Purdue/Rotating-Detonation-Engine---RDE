%% PURPL RDE Chamber Cooldown Time v3
% Jon DeSimone

%% Solid wall geometry
r_outer = 0.066675; %outer wall radius [m], 2.625" 
r_inner = 0.0508; %inner wall radius [m], 2" 
t_wall = r_outer - r_inner; % wall thickness [m]
length = 0.093345; % axial length [m], 3.675"
rho_wall = 7930; % SS316 density [kg/m^3]
cp_wall = 500; % SS316 spec heat cap [J/kg-K]
k_wall = 15; % SS316 thermal cond [W/m-K]

%% Temperatures
t_wall_i_c = 1188.12; %Initial wall temp after firing [C]
t_surr_c = 10; %Ambient temp [C]
t_wall_i = t_wall_i_c + 273.15; %Initial wall temp [K]
t_surr = t_surr_c + 273.15; %Ambient temp [K]

%% Air properties

h_air_outer_base = 10; %Baseline external convective coeff of air [W/m^2-K]
h_air_inner_max = 25; %Maximum internal convective coeff at high temp [W/m^2-K]
h_air_inner_min = 10; %Minimum internal convective coeff at low temp [W/m^2-K]
rho_air = 0.3; %Density of air at ~1200C [kg/m^3]
cp_air = 1150; %Specific heat of air at high temp [J/kg-K]

%% Enhanced cooling scenarios
h_air_outer_rough = h_air_outer_base * 1.8; %Rough surface (increased turbulence)
h_air_n2 = 60; %N2 forced convection cooling [W/m^2-K]
t_n2_start = 0.5*3600; %Start N2 cooling after 30 minutes [s]

%% Radiation properties
epsilon = 0.6; %Emissivity of oxidized SS316 (typical range: 0.5-0.85)
sigma = 5.67e-8; %Stefan-Boltzmann constant [W/m^2-K^4]

%% Wall calculations
v_wall = pi*length*(r_outer^2 - r_inner^2); %Wall volume [m^3]
m_wall = rho_wall*v_wall; %Wall mass [kg]
as_wall_outer = 2*pi*r_outer*length; %Outer wall surface area [m^2]
as_wall_inner = 2*pi*r_inner*length; %Inner wall surface area [m^2]
lc = v_wall/as_wall_outer; %Wall characteristic length [m]

%% Air calculations
v_air = pi*r_inner^2*length; %Air volume inside chamber [m^3]
m_air = rho_air*v_air; %Air mass [kg]

%% Biot number validation
bi = h_air_outer_base*lc/k_wall;
fprintf('Bi = %.4f (lumped capacitance %s)\n', bi, ...
    ternary(bi < 0.1, 'valid', 'questionable'));

%% Solve four scenarios
tspan = [0 3*3600]; % [s] 0 to 3 hours
Y0 = [t_wall_i; t_wall_i]; % Initial conditions

fprintf('\n=== RUNNING SIMULATIONS ===\n\n');

% Scenario 1: Baseline
fprintf('Scenario 1: Baseline\n');
cooldownODE_base = @(t, Y) twoNodeModel_tempDependent(t, Y, t_surr, t_wall_i, ...
    m_air, cp_air, m_wall, cp_wall, h_air_inner_min, h_air_inner_max, ...
    h_air_outer_base, as_wall_inner, as_wall_outer, epsilon, sigma);
[t_base, Y_base] = ode45(cooldownODE_base, tspan, Y0);

% Scenario 2: With N2 cooling starting at t_n2_start
fprintf('Scenario 2: N2 cooling after %.1f minutes\n', t_n2_start/60);
cooldownODE_N2 = @(t, Y) twoNodeModel_N2(t, Y, t_surr, t_wall_i, ...
    m_air, cp_air, m_wall, cp_wall, h_air_inner_min, h_air_inner_max, ...
    h_air_outer_base, h_air_n2, t_n2_start, as_wall_inner, as_wall_outer, epsilon, sigma);
[t_N2, Y_N2] = ode45(cooldownODE_N2, tspan, Y0);

% Scenario 3: Rough surface (higher h_outer)
fprintf('Scenario 3: Rough surface (h_outer = %.1f W/m²-K)\n', h_air_outer_rough);
cooldownODE_rough = @(t, Y) twoNodeModel_tempDependent(t, Y, t_surr, t_wall_i, ...
    m_air, cp_air, m_wall, cp_wall, h_air_inner_min, h_air_inner_max, ...
    h_air_outer_rough, as_wall_inner, as_wall_outer, epsilon, sigma);
[t_rough, Y_rough] = ode45(cooldownODE_rough, tspan, Y0);

% Scenario 4: Rough surface + N2 cooling
fprintf('Scenario 4: Rough surface + N2 cooling\n');
cooldownODE_combo = @(t, Y) twoNodeModel_N2(t, Y, t_surr, t_wall_i, ...
    m_air, cp_air, m_wall, cp_wall, h_air_inner_min, h_air_inner_max, ...
    h_air_outer_rough, h_air_n2, t_n2_start, as_wall_inner, as_wall_outer, epsilon, sigma);
[t_combo, Y_combo] = ode45(cooldownODE_combo, tspan, Y0);

%% Extract temperatures for all scenarios
t_hours_base = t_base/3600;
T_wall_celsius_base = Y_base(:,2) - 273.15;

t_hours_N2 = t_N2/3600;
T_wall_celsius_N2 = Y_N2(:,2) - 273.15;

t_hours_rough = t_rough/3600;
T_wall_celsius_rough = Y_rough(:,2) - 273.15;

t_hours_combo = t_combo/3600;
T_wall_celsius_combo = Y_combo(:,2) - 273.15;

%% Find cooling times for 95% and 90% targets
t_target_95 = t_surr + 0.05*(t_wall_i - t_surr);
t_target_95_c = t_target_95 - 273.15;

t_target_90 = t_surr + 0.10*(t_wall_i - t_surr);
t_target_90_c = t_target_90 - 273.15;

% Baseline scenario
[idx_95_base, t_95_base, T_95_base] = findCoolingTime(t_base, Y_base(:,2), t_target_95);
[idx_90_base, t_90_base, T_90_base] = findCoolingTime(t_base, Y_base(:,2), t_target_90);

% N2 cooling scenario
[idx_95_N2, t_95_N2, T_95_N2] = findCoolingTime(t_N2, Y_N2(:,2), t_target_95);
[idx_90_N2, t_90_N2, T_90_N2] = findCoolingTime(t_N2, Y_N2(:,2), t_target_90);

% Rough surface scenario
[idx_95_rough, t_95_rough, T_95_rough] = findCoolingTime(t_rough, Y_rough(:,2), t_target_95);
[idx_90_rough, t_90_rough, T_90_rough] = findCoolingTime(t_rough, Y_rough(:,2), t_target_90);

% Combo scenario (rough + N2)
[idx_95_combo, t_95_combo, T_95_combo] = findCoolingTime(t_combo, Y_combo(:,2), t_target_95);
[idx_90_combo, t_90_combo, T_90_combo] = findCoolingTime(t_combo, Y_combo(:,2), t_target_90);

%% Display results
fprintf('\n=== COOLING TIME COMPARISON ===\n\n');

fprintf('To 95%% Cooling (%.1f°C):\n', t_target_95_c);
fprintf('  Baseline:      %.2f hrs (%.1f min)\n', t_95_base/3600, t_95_base/60);
fprintf('  N2 Cooling:    %.2f hrs (%.1f min)\n', t_95_N2/3600, t_95_N2/60);
fprintf('  Rough Surface: %.2f hrs (%.1f min)\n', t_95_rough/3600, t_95_rough/60);
fprintf('  Rough + N2:    %.2f hrs (%.1f min)\n\n', t_95_combo/3600, t_95_combo/60);

fprintf('To 90%% Cooling (%.1f°C):\n', t_target_90_c);
fprintf('  Baseline:      %.2f hrs (%.1f min)\n', t_90_base/3600, t_90_base/60);
fprintf('  N2 Cooling:    %.2f hrs (%.1f min)\n', t_90_N2/3600, t_90_N2/60);
fprintf('  Rough Surface: %.2f hrs (%.1f min)\n', t_90_rough/3600, t_90_rough/60);
fprintf('  Rough + N2:    %.2f hrs (%.1f min)\n\n', t_90_combo/3600, t_90_combo/60);

fprintf('Time Savings vs Baseline (to 90%%):\n');
fprintf('  N2 Cooling:    %.1f min (%.1f%% faster)\n', ...
    (t_90_base - t_90_N2)/60, 100*(t_90_base - t_90_N2)/t_90_base);
fprintf('  Rough Surface: %.1f min (%.1f%% faster)\n', ...
    (t_90_base - t_90_rough)/60, 100*(t_90_base - t_90_rough)/t_90_base);
fprintf('  Rough + N2:    %.1f min (%.1f%% faster)\n\n', ...
    (t_90_base - t_90_combo)/60, 100*(t_90_base - t_90_combo)/t_90_base);

%% Plotting
figure('Position', [100, 100, 1200, 700]);

% Plot all four scenarios with distinct colors
plot(t_hours_base, T_wall_celsius_base, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Baseline');
hold on;
plot(t_hours_N2, T_wall_celsius_N2, 'r-', 'LineWidth', 2.5, 'DisplayName', 'N₂ Cooling');
plot(t_hours_rough, T_wall_celsius_rough, 'Color', [0 0.6 0], 'LineWidth', 2.5, 'DisplayName', 'Rough Surface');
plot(t_hours_combo, T_wall_celsius_combo, 'm-', 'LineWidth', 2.5, 'DisplayName', 'Rough + N₂');

% Plot target lines (thinner and less prominent)
plot([0 3], [t_target_95_c t_target_95_c], 'k--', 'LineWidth', 1, 'DisplayName', '95% Target', 'Color', [0.5 0.5 0.5]);
plot([0 3], [t_target_90_c t_target_90_c], 'k:', 'LineWidth', 1.5, 'DisplayName', '90% Target', 'Color', [0.3 0.3 0.3]);

% Mark 95% cooling points (circles) - smaller and less prominent
% Use different vertical offsets to prevent overlap
offset_95_base = 100;
offset_95_N2 = 80;
offset_95_rough = 60;
offset_95_combo = 40;

if ~isnan(idx_95_base)
    plot(t_95_base/3600, T_95_base - 273.15, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
    text(t_95_base/3600, T_95_base - 273.15 + offset_95_base, ...
        sprintf('%.0f min', t_95_base/60), 'FontSize', 8, 'Color', 'b', 'FontWeight', 'normal', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'b', 'Margin', 1);
end
if ~isnan(idx_95_N2)
    plot(t_95_N2/3600, T_95_N2 - 273.15, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    text(t_95_N2/3600, T_95_N2 - 273.15 + offset_95_N2, ...
        sprintf('%.0f min', t_95_N2/60), 'FontSize', 8, 'Color', 'r', 'FontWeight', 'normal', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'r', 'Margin', 1);
end
if ~isnan(idx_95_rough)
    plot(t_95_rough/3600, T_95_rough - 273.15, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.6 0], 'MarkerEdgeColor', [0 0.6 0], 'HandleVisibility', 'off');
    text(t_95_rough/3600, T_95_rough - 273.15 + offset_95_rough, ...
        sprintf('%.0f min', t_95_rough/60), 'FontSize', 8, 'Color', [0 0.6 0], 'FontWeight', 'normal', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', [0 0.6 0], 'Margin', 1);
end
if ~isnan(idx_95_combo)
    plot(t_95_combo/3600, T_95_combo - 273.15, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'HandleVisibility', 'off');
    text(t_95_combo/3600, T_95_combo - 273.15 + offset_95_combo, ...
        sprintf('%.0f min', t_95_combo/60), 'FontSize', 8, 'Color', 'm', 'FontWeight', 'normal', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'm', 'Margin', 1);
end

% Mark 90% cooling points (squares) with time labels
% Use different vertical offsets to prevent overlap
offset_base = 80;
offset_N2 = 60;
offset_rough = 40;
offset_combo = 20;

if ~isnan(idx_90_base)
    plot(t_90_base/3600, T_90_base - 273.15, 'bs', 'MarkerSize', 11, 'MarkerFaceColor', 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    text(t_90_base/3600, T_90_base - 273.15 + offset_base, ...
        sprintf('%.0f min', t_90_base/60), 'FontSize', 9, 'Color', 'b', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'b', 'Margin', 2);
end
if ~isnan(idx_90_N2)
    plot(t_90_N2/3600, T_90_N2 - 273.15, 'rs', 'MarkerSize', 11, 'MarkerFaceColor', 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    text(t_90_N2/3600, T_90_N2 - 273.15 + offset_N2, ...
        sprintf('%.0f min', t_90_N2/60), 'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'r', 'Margin', 2);
end
if ~isnan(idx_90_rough)
    plot(t_90_rough/3600, T_90_rough - 273.15, 's', 'MarkerSize', 11, 'MarkerFaceColor', [0 0.6 0], 'MarkerEdgeColor', [0 0.6 0], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    text(t_90_rough/3600, T_90_rough - 273.15 + offset_rough, ...
        sprintf('%.0f min', t_90_rough/60), 'FontSize', 9, 'Color', [0 0.6 0], 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', [0 0.6 0], 'Margin', 2);
end
if ~isnan(idx_90_combo)
    plot(t_90_combo/3600, T_90_combo - 273.15, 'ms', 'MarkerSize', 11, 'MarkerFaceColor', 'm', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    text(t_90_combo/3600, T_90_combo - 273.15 + offset_combo, ...
        sprintf('%.0f min', t_90_combo/60), 'FontSize', 9, 'Color', 'm', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'm', 'Margin', 2);
end

% Mark N2 start time (subtle vertical line)
xline(t_n2_start/3600, '--', 'LineWidth', 1, 'HandleVisibility', 'off', 'Alpha', 0.4, 'Color', [0.7 0.3 0.3]);
text(t_n2_start/3600, 1150, 'N₂ Start', 'FontSize', 9, 'Color', [0.7 0.3 0.3], ...
    'HorizontalAlignment', 'center', 'FontWeight', 'normal', 'Rotation', 90);

% Initial temperature marker (subtle)
plot(0, t_wall_i_c, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.3 0.3], 'HandleVisibility', 'off');

grid on;
xlabel('Time (hours)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Wall Temperature (°C)', 'FontSize', 13, 'FontWeight', 'bold');
title('RDE Chamber Cooldown Strategies', 'FontSize', 15, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11);
ylim([0 1250]);
xlim([0 3]);


%% Helper function for finding cooling time
function [idx, t_cool, T_cool] = findCoolingTime(t, T, T_target)
    idx = find(T <= T_target, 1, 'first');
    if ~isempty(idx)
        t_cool = t(idx);
        T_cool = T(idx);
    else
        idx = NaN;
        t_cool = NaN;
        T_cool = NaN;
    end
end

%% Helper function for ternary operator
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

%% ODE Function: Temperature-dependent h_inner (Baseline and Rough Surface)
function dYdt = twoNodeModel_tempDependent(~, Y, t_surr, t_wall_i, m_air, cp_air, ...
    m_wall, cp_wall, h_air_inner_min, h_air_inner_max, h_air_outer, as_wall_inner, ...
    as_wall_outer, epsilon, sigma)
    
    % Extract state variables
    t_air = Y(1);   % Air temperature [K]
    t_wall = Y(2);  % Wall temperature [K]
    
    % Temperature-dependent internal convection coefficient
    % Higher convection at high temperatures (stronger buoyancy)
    temp_ratio = max(0, (t_air - t_surr) / (t_wall_i - t_surr));
    h_air_inner = h_air_inner_min + (h_air_inner_max - h_air_inner_min) * temp_ratio^0.5;
    
    % Heat transfer terms
    q_air_to_wall = h_air_inner * as_wall_inner * (t_air - t_wall);
    q_wall_conv = h_air_outer * as_wall_outer * (t_wall - t_surr);
    q_wall_rad = epsilon * sigma * as_wall_outer * (t_wall^4 - t_surr^4);
    
    % Rate equations
    dt_air_dt = -q_air_to_wall / (m_air * cp_air);
    dt_wall_dt = (q_air_to_wall - q_wall_conv - q_wall_rad) / (m_wall * cp_wall);
    
    dYdt = [dt_air_dt; dt_wall_dt];
end

%% ODE Function: With N2 cooling (switches h_outer at t_n2_start)
function dYdt = twoNodeModel_N2(t, Y, t_surr, t_wall_i, m_air, cp_air, ...
    m_wall, cp_wall, h_air_inner_min, h_air_inner_max, h_air_outer_base, h_air_n2, ...
    t_n2_start, as_wall_inner, as_wall_outer, epsilon, sigma)
    
    % Extract state variables
    t_air = Y(1);   % Air temperature [K]
    t_wall = Y(2);  % Wall temperature [K]
    
    % Temperature-dependent internal convection coefficient
    temp_ratio = max(0, (t_air - t_surr) / (t_wall_i - t_surr));
    h_air_inner = h_air_inner_min + (h_air_inner_max - h_air_inner_min) * temp_ratio^0.5;
    
    % Switch to N2 cooling after t_n2_start
    if t >= t_n2_start
        h_air_outer = h_air_n2;
    else
        h_air_outer = h_air_outer_base;
    end
    
    % Heat transfer terms
    q_air_to_wall = h_air_inner * as_wall_inner * (t_air - t_wall);
    q_wall_conv = h_air_outer * as_wall_outer * (t_wall - t_surr);
    q_wall_rad = epsilon * sigma * as_wall_outer * (t_wall^4 - t_surr^4);
    
    % Rate equations
    dt_air_dt = -q_air_to_wall / (m_air * cp_air);
    dt_wall_dt = (q_air_to_wall - q_wall_conv - q_wall_rad) / (m_wall * cp_wall);
    
    dYdt = [dt_air_dt; dt_wall_dt];
end