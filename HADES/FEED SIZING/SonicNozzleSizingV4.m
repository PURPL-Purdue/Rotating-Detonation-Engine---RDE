%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Purpose:
%
%   Calculate Required RDE Manifold Pressure from
%   Injector Area and Desired Mass Flow Rates in a CSV
%
%   Assumes each injector array is choked.
%
%   Input CSV columns:
%       H2, Air, CH4, GOx
%
%   Each column contains total required mass flow [kg/s].
%
%   Output:
%       manifold_pressure_results.csv
%
%   Programmers: Deepesh Balwani, Noah Ha
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;

%% =====================================================
%  USER INPUTS
% ======================================================

input_file = 'flow_requirements_manifold_input.csv';

% Mean chamber pressure [bar absolute]
%
% This is used only to CHECK whether the injector is
% actually choked at the calculated manifold pressure.

P_chamber_bar = 10.0;

%% =====================================================
%  TOTAL INJECTOR AREAS
% ======================================================

% Total physical injector flow area [mm^2]

A_H2_mm2  = 25.13;
A_Air_mm2 = 304.77;
A_CH4_mm2 = 56.55;
A_GOx_mm2 = 86.59;

% Convert to m^2

A_H2  = A_H2_mm2  * 1e-6;
A_Air = A_Air_mm2 * 1e-6;
A_CH4 = A_CH4_mm2 * 1e-6;
A_GOx = A_GOx_mm2 * 1e-6;

%% =====================================================
%  EFFECTIVE INJECTOR DISCHARGE COEFFICIENTS
% ======================================================

% Start with Cd = 1.0 for an ideal choked-flow calculation.
%
% Replace these values later with experimentally
% determined effective injector discharge coefficients.

Cd_H2  = 0.71;
Cd_Air = 0.62;
Cd_CH4 = 0.6;
Cd_GOx = 0.6;

%% =====================================================
%  GAS PROPERTIES
% ======================================================

% Manifold stagnation temperature [K]

T0_H2  = 283;
T0_Air = 283;
T0_CH4 = 283;
T0_GOx = 283;

% Specific gas constants [J/(kg*K)]

R_H2  = 4124.0;
R_Air = 287.05;
R_CH4 = 518.3;
R_GOx = 259.84;

% Specific heat ratios [-]

gamma_H2  = 1.405;
gamma_Air = 1.400;
gamma_CH4 = 1.310;
gamma_GOx = 1.400;

%% =====================================================
%  READ CSV
% ======================================================

data = readtable( ...
    input_file, ...
    'VariableNamingRule', 'preserve');

required_columns = {'H2', 'Air', 'CH4', 'GOx'};

for i = 1:length(required_columns)

    if ~any(strcmpi( ...
            data.Properties.VariableNames, ...
            required_columns{i}))

        error( ...
            'CSV must contain a column named "%s".', ...
            required_columns{i});

    end

end

column_names = data.Properties.VariableNames;

idx_H2  = find(strcmpi(column_names, 'H2'), 1);
idx_Air = find(strcmpi(column_names, 'Air'), 1);
idx_CH4 = find(strcmpi(column_names, 'CH4'), 1);
idx_GOx = find(strcmpi(column_names, 'GOx'), 1);

m_dot_H2  = data{:, idx_H2};
m_dot_Air = data{:, idx_Air};
m_dot_CH4 = data{:, idx_CH4};
m_dot_GOx = data{:, idx_GOx};

%% =====================================================
%  CHOKED MASS-FLOW PARAMETERS
% ======================================================

% For choked flow:
%
% mdot =
%
% Cd * A * P0 / sqrt(T0)
%
% *
%
% sqrt(gamma/R)
%
% *
%
% (2/(gamma+1))^((gamma+1)/(2*(gamma-1)))
%
%
% Define:
%
% K =
%
% sqrt(gamma/R)
%
% *
%
% (2/(gamma+1))^((gamma+1)/(2*(gamma-1)))

K_H2 = ...
    sqrt(gamma_H2 / R_H2) * ...
    (2 / (gamma_H2 + 1)) ^ ...
    ((gamma_H2 + 1) / ...
    (2 * (gamma_H2 - 1)));

K_Air = ...
    sqrt(gamma_Air / R_Air) * ...
    (2 / (gamma_Air + 1)) ^ ...
    ((gamma_Air + 1) / ...
    (2 * (gamma_Air - 1)));

K_CH4 = ...
    sqrt(gamma_CH4 / R_CH4) * ...
    (2 / (gamma_CH4 + 1)) ^ ...
    ((gamma_CH4 + 1) / ...
    (2 * (gamma_CH4 - 1)));

K_GOx = ...
    sqrt(gamma_GOx / R_GOx) * ...
    (2 / (gamma_GOx + 1)) ^ ...
    ((gamma_GOx + 1) / ...
    (2 * (gamma_GOx - 1)));

%% =====================================================
%  REQUIRED MANIFOLD PRESSURES
% ======================================================

% Rearranged choked-flow equation:
%
% P0_manifold =
%
% mdot * sqrt(T0)
% --------------------------------
% Cd * A * K

Pman_H2 = ...
    (m_dot_H2 .* sqrt(T0_H2)) ./ ...
    (Cd_H2 * A_H2 * K_H2);

Pman_Air = ...
    (m_dot_Air .* sqrt(T0_Air)) ./ ...
    (Cd_Air * A_Air * K_Air);

Pman_CH4 = ...
    (m_dot_CH4 .* sqrt(T0_CH4)) ./ ...
    (Cd_CH4 * A_CH4 * K_CH4);

Pman_GOx = ...
    (m_dot_GOx .* sqrt(T0_GOx)) ./ ...
    (Cd_GOx * A_GOx * K_GOx);

% Convert Pa -> bar absolute

Pman_H2_bar  = Pman_H2  / 1e5;
Pman_Air_bar = Pman_Air / 1e5;
Pman_CH4_bar = Pman_CH4 / 1e5;
Pman_GOx_bar = Pman_GOx / 1e5;

%% =====================================================
%  CRITICAL PRESSURE RATIOS
% ======================================================

% Choked if:
%
% P_chamber / P0_manifold <= critical pressure ratio

crit_H2 = ...
    (2 / (gamma_H2 + 1)) ^ ...
    (gamma_H2 / (gamma_H2 - 1));

crit_Air = ...
    (2 / (gamma_Air + 1)) ^ ...
    (gamma_Air / (gamma_Air - 1));

crit_CH4 = ...
    (2 / (gamma_CH4 + 1)) ^ ...
    (gamma_CH4 / (gamma_CH4 - 1));

crit_GOx = ...
    (2 / (gamma_GOx + 1)) ^ ...
    (gamma_GOx / (gamma_GOx - 1));

%% =====================================================
%  MINIMUM MANIFOLD PRESSURE REQUIRED FOR CHOKING
% ======================================================

Pman_min_choke_H2_bar = ...
    P_chamber_bar / crit_H2;

Pman_min_choke_Air_bar = ...
    P_chamber_bar / crit_Air;

Pman_min_choke_CH4_bar = ...
    P_chamber_bar / crit_CH4;

Pman_min_choke_GOx_bar = ...
    P_chamber_bar / crit_GOx;

%% =====================================================
%  CHECK CHOKING
% ======================================================

H2_Choked = ...
    (P_chamber_bar ./ Pman_H2_bar) <= crit_H2;

Air_Choked = ...
    (P_chamber_bar ./ Pman_Air_bar) <= crit_Air;

CH4_Choked = ...
    (P_chamber_bar ./ Pman_CH4_bar) <= crit_CH4;

GOx_Choked = ...
    (P_chamber_bar ./ Pman_GOx_bar) <= crit_GOx;

% Preserve blank / NaN rows

H2_Choked(isnan(m_dot_H2)) = false;
Air_Choked(isnan(m_dot_Air)) = false;
CH4_Choked(isnan(m_dot_CH4)) = false;
GOx_Choked(isnan(m_dot_GOx)) = false;

%% =====================================================
%  CREATE OUTPUT TABLE
% ======================================================

output = data;

output.H2_ManifoldPressure_bar = ...
    Pman_H2_bar;

output.H2_Choked = ...
    H2_Choked;

output.Air_ManifoldPressure_bar = ...
    Pman_Air_bar;

output.Air_Choked = ...
    Air_Choked;

output.CH4_ManifoldPressure_bar = ...
    Pman_CH4_bar;

output.CH4_Choked = ...
    CH4_Choked;

output.GOx_ManifoldPressure_bar = ...
    Pman_GOx_bar;

output.GOx_Choked = ...
    GOx_Choked;

%% =====================================================
%  WRITE OUTPUT CSV
% ======================================================

output_file = ...
    'manifold_pressure_results.csv';

writetable( ...
    output, ...
    output_file);

%% =====================================================
%  COMMAND WINDOW SUMMARY
% ======================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('              RDE MANIFOLD PRESSURE MODEL\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf('Chamber Pressure: %.2f bar absolute\n', ...
    P_chamber_bar);

fprintf('\n');
fprintf('Minimum manifold pressures for ideal choking:\n');

fprintf('H2  : %.3f bar\n', ...
    Pman_min_choke_H2_bar);

fprintf('Air : %.3f bar\n', ...
    Pman_min_choke_Air_bar);

fprintf('CH4 : %.3f bar\n', ...
    Pman_min_choke_CH4_bar);

fprintf('GOx : %.3f bar\n', ...
    Pman_min_choke_GOx_bar);

fprintf('\n');
fprintf('Injector areas:\n');

fprintf('H2  : %.2f mm^2\n', A_H2_mm2);
fprintf('Air : %.2f mm^2\n', A_Air_mm2);
fprintf('CH4 : %.2f mm^2\n', A_CH4_mm2);
fprintf('GOx : %.2f mm^2\n', A_GOx_mm2);

fprintf('\n');
fprintf('Effective injector Cd values:\n');

fprintf('H2  : %.3f\n', Cd_H2);
fprintf('Air : %.3f\n', Cd_Air);
fprintf('CH4 : %.3f\n', Cd_CH4);
fprintf('GOx : %.3f\n', Cd_GOx);

fprintf('\n');
fprintf('Results written to:\n%s\n', output_file);

fprintf('============================================================\n');