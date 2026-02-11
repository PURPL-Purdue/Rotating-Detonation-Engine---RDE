%% PURPL RDE Chamber Cooldown Time with Two-Node Model
% Jon DeSimone

%% Solid wall geometry
ro = 0.066675; %outer wall radius [m], 2.625" 
ri = 0.0508; %inner wall radius [m], 2" 
tWall = ro - ri; % wall thickness [m]
L  = 0.093345; % axial length [m], 3.675"
rhoWall = 7930; % SS316 density [kg/m^3]
cpWall = 500; % SS316 spec heat cap [J/kg-K]
kWall = 15; % SS316 thermal cond [W/m-K]

%% Temperatures (convert to Kelvin for radiation calculations)
TWalli_C = 1188.12; %Initial wall temp after firing [C]
TSurr_C = 10; %Ambient temp [C]
TWalli = TWalli_C + 273.15; %Initial wall temp [K]
TSurr = TSurr_C + 273.15; %Ambient temp [K]

%% Air properties
% Note: Using properties at elevated temperature for better accuracy
hAir_outer = 10; %External convective coeff of air [W/m^2-K]
hAir_inner = 10; %Internal convective coeff of air [W/m^2-K] (natural convection)
rhoAir = 0.3; %Density of air at ~1200C [kg/m^3] (much lower than sea level due to temp)
cpAir = 1150; %Specific heat of air at high temp [J/kg-K]

%% Radiation properties
epsilon = 0.6; %Emissivity of oxidized SS316 (typical range: 0.5-0.85)
sigma = 5.67e-8; %Stefan-Boltzmann constant [W/m^2-K^4]

%% Wall calculations
VWall = pi*L*(ro^2 - ri^2); %Wall volume [m^3]
mWall = rhoWall*VWall; %Wall mass [kg]
AsWall_outer = 2*pi*ro*L; %Outer wall surface area [m^2]
AsWall_inner = 2*pi*ri*L; %Inner wall surface area [m^2]
Lc = VWall/AsWall_outer; %Wall characteristic length [m]

%% Air calculations
VAir = pi*ri^2*L; %Air volume inside chamber [m^3]
mAir = rhoAir*VAir; %Air mass [kg]

%% Thermal mass comparison
C_wall = mWall * cpWall; %Wall thermal capacitance [J/K]
C_air = mAir * cpAir; %Air thermal capacitance [J/K]
fprintf('\n=== THERMAL MASS ANALYSIS ===\n');
fprintf('Wall thermal mass: %.2f J/K\n', C_wall);
fprintf('Air thermal mass: %.2f J/K\n', C_air);
fprintf('Ratio (C_wall/C_air): %.1f\n', C_wall/C_air);
fprintf('Wall will cool %.1fx slower than air\n\n', C_wall/C_air);

%% Biot number validation
Bi = hAir_outer*Lc/kWall; %Biot number at given conditions
if Bi < 0.1
    fprintf('Bi = %.4f < 0.1, lumped capacitance valid for wall\n', Bi);
else
    fprintf('Bi = %.4f > 0.1, lumped capacitance may be invalid\n', Bi);
end

% Check air equilibration time
tau_air = (mAir * cpAir) / (hAir_inner * AsWall_inner);
fprintf('Air equilibration time constant: %.2f seconds\n\n', tau_air);

%% Solve the coupled ODEs for two-node model
% Node 1: Air inside chamber (T_air)
% Node 2: Wall (T_wall)
%
% Energy balance for air:
%   m_air * cp_air * dT_air/dt = -h_inner * A_inner * (T_air - T_wall)
%
% Energy balance for wall:
%   m_wall * cp_wall * dT_wall/dt = h_inner * A_inner * (T_air - T_wall)
%                                    - h_outer * A_outer * (T_wall - T_surr)
%                                    - epsilon * sigma * A_outer * (T_wall^4 - T_surr^4)

% Define the coupled ODE system
% State vector: Y = [T_air; T_wall]
cooldownODE = @(t, Y) twoNodeModel(t, Y, TSurr, mAir, cpAir, mWall, cpWall, ...
                                   hAir_inner, hAir_outer, AsWall_inner, ...
                                   AsWall_outer, epsilon, sigma);

% Time span (simulate for 3 hours)
tspan = [0 3*3600]; % [s] 0 to 3 hours

% Initial conditions: both air and wall start at same temperature
Y0 = [TWalli; TWalli];

% Solve ODE system
[t, Y] = ode45(cooldownODE, tspan, Y0);

% Extract temperatures
T_air = Y(:,1); % Air temperature [K]
T_wall = Y(:,2); % Wall temperature [K]

% Convert time to hours and temperatures to Celsius
t_hours = t/3600;
T_air_celsius = T_air - 273.15;
T_wall_celsius = T_wall - 273.15;

%% Find time to reach 95% of original temperature (based on wall temp)
% Target temperature: TSurr + 0.05*(TWalli - TSurr)
TTarget = TSurr + 0.05*(TWalli - TSurr);
TTarget_C = TTarget - 273.15;

% Find when wall temperature drops to target
idx = find(T_wall <= TTarget, 1, 'first');
if ~isempty(idx)
    tCool95 = t(idx);
    tCool95_hours = tCool95/3600;
    tCool95_minutes = tCool95/60;
    T95_wall_celsius = T_wall_celsius(idx);
    T95_air_celsius = T_air_celsius(idx);
    fprintf('\n=== COOLING TIME TO 95%% OF ORIGINAL TEMPERATURE ===\n');
    fprintf('Time to cool wall to %.2f°C: %.2f hours (%.2f minutes)\n', ...
            TTarget_C, tCool95_hours, tCool95_minutes);
    fprintf('Wall temperature at this time: %.2f°C\n', T95_wall_celsius);
    fprintf('Air temperature at this time: %.2f°C\n\n', T95_air_celsius);
else
    fprintf('Wall did not reach target temperature within 3 hours.\n');
    fprintf('Final wall temperature after 3 hours: %.2f°C\n', T_wall_celsius(end));
    fprintf('Final air temperature after 3 hours: %.2f°C\n\n', T_air_celsius(end));
    tCool95_hours = NaN;
end

%% Calculate heat transfer contributions at initial conditions
q_conv_outer_initial = hAir_outer*AsWall_outer*(TWalli - TSurr);
q_rad_initial = epsilon*sigma*AsWall_outer*(TWalli^4 - TSurr^4);
q_conv_inner_initial = hAir_inner*AsWall_inner*(TWalli - TWalli); % Zero initially
q_total_initial = q_conv_outer_initial + q_rad_initial;

fprintf('=== HEAT TRANSFER ANALYSIS AT INITIAL TEMPERATURE ===\n');
fprintf('Initial wall temperature: %.2f°C\n', TWalli_C);
fprintf('Initial air temperature: %.2f°C\n', TWalli_C);
fprintf('Ambient temperature: %.2f°C\n', TSurr_C);
fprintf('External convective heat loss: %.2f W (%.1f%%)\n', q_conv_outer_initial, ...
        100*q_conv_outer_initial/q_total_initial);
fprintf('Radiative heat loss: %.2f W (%.1f%%)\n', q_rad_initial, ...
        100*q_rad_initial/q_total_initial);
fprintf('Total initial heat loss: %.2f W\n\n', q_total_initial);

%% Plotting
figure('Position', [100, 100, 1000, 700]);

% Main temperature plot - wall temperature only
plot(t_hours, T_wall_celsius, 'b-', 'LineWidth', 2.5);
hold on;
plot([0 max(t_hours)], [TTarget_C TTarget_C], 'r--', 'LineWidth', 1.5);

% Mark the 95% cooldown point
if ~isnan(tCool95_hours)
    plot(tCool95_hours, T95_wall_celsius, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    text(tCool95_hours*1.05, T95_wall_celsius, ...
         sprintf('  %.2f hrs (%.1f min)\n  %.1f°C', tCool95_hours, tCool95_minutes, T95_wall_celsius), ...
         'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');
end

% Plot initial temperature point
plot(0, TWalli_C, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
text(0.05, TWalli_C, sprintf('  Initial: %.1f°C', TWalli_C), ...
     'FontSize', 11, 'Color', 'g', 'FontWeight', 'bold');

% Plot final temperature point
plot(max(t_hours), T_wall_celsius(end), 'mo', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
text(max(t_hours)*0.95, T_wall_celsius(end), sprintf('%.1f°C  ', T_wall_celsius(end)), ...
     'FontSize', 10, 'Color', 'm', 'HorizontalAlignment', 'right', 'FontWeight', 'bold');

grid on;
xlabel('Time (hours)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Wall Temperature (°C)', 'FontSize', 13, 'FontWeight', 'bold');
title('RDE Chamber Cooldown: Two-Node Model', 'FontSize', 15, 'FontWeight', 'bold');
legend('Wall Temperature', '95% Cooldown Target', 'Target Reached', 'Initial Temp', 'Final Temp', ...
       'Location', 'northeast', 'FontSize', 11);
set(gca, 'FontSize', 11);
ylim([0 max(TWalli_C*1.05, 1300)]);
xlim([0 3]);

fprintf('=== SIMULATION COMPLETE ===\n');
fprintf('Final wall temperature after 3 hours: %.2f°C\n', T_wall_celsius(end));
fprintf('Final air temperature after 3 hours: %.2f°C\n', T_air_celsius(end));
fprintf('Wall temperature drop: %.2f°C\n', TWalli_C - T_wall_celsius(end));
fprintf('Air temperature drop: %.2f°C\n', TWalli_C - T_air_celsius(end));

%% ODE Function Definition
function dYdt = twoNodeModel(t, Y, TSurr, mAir, cpAir, mWall, cpWall, ...
                             hAir_inner, hAir_outer, AsWall_inner, ...
                             AsWall_outer, epsilon, sigma)
    % Extract state variables
    T_air = Y(1);   % Air temperature [K]
    T_wall = Y(2);  % Wall temperature [K]
    
    % Heat transfer terms
    % Air to wall (internal convection)
    q_air_to_wall = hAir_inner * AsWall_inner * (T_air - T_wall);
    
    % Wall to surroundings (external convection)
    q_wall_conv = hAir_outer * AsWall_outer * (T_wall - TSurr);
    
    % Wall to surroundings (radiation)
    q_wall_rad = epsilon * sigma * AsWall_outer * (T_wall^4 - TSurr^4);
    
    % Rate equations
    % Air energy balance
    dT_air_dt = -q_air_to_wall / (mAir * cpAir);
    
    % Wall energy balance
    dT_wall_dt = (q_air_to_wall - q_wall_conv - q_wall_rad) / (mWall * cpWall);
    
    % Return derivatives
    dYdt = [dT_air_dt; dT_wall_dt];
end