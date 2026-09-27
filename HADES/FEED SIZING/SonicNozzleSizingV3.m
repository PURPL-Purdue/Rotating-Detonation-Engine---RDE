%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Purpose:
%
%   Two-Stage Choked Flow Sizing from CSV
%
%   Stage 1:
%       Size fixed injector effective area using:
%
%       Chamber Pressure -> Injector -> Manifold Pressure
%
%   Stage 2:
%       Size fixed upstream sonic orifice using:
%
%       Manifold Pressure -> Sonic Orifice -> Supply P0
%
%   Input CSV columns:
%       H2, Air, CH4, GOx
%
%   Each column contains required mass flow rates [kg/s].
%
%   The lowest positive mass flow for each gas is used
%   to size:
%
%       1. One fixed injector effective flow area
%       2. One fixed upstream sonic orifice
%
%   The code then calculates:
%
%       - Required manifold pressure for every flow point
%       - Required upstream P0 for every flow point
%       - Injector choking status
%       - Sonic orifice choking status
%
%   Programmers: Deepesh Balwani, Noah Ha
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;

%% =====================================================
%  USER INPUTS
% ======================================================

% Input CSV filename

input_file = 'flow_requirements.csv';

%% CHAMBER PRESSURE

% Both engines assumed to operate at the same
% chamber pressure.

P_chamber_bar = 10;         % Chamber pressure [bar absolute]

%% INJECTOR DISCHARGE COEFFICIENTS

% Effective discharge coefficients for the engine
% injectors.

Cd_inj_H2  = 0.85;
Cd_inj_Air = 0.85;
Cd_inj_CH4 = 0.85;
Cd_inj_GOx = 0.85;

%% SONIC ORIFICE DISCHARGE COEFFICIENTS

% Discharge coefficients for the upstream metering
% orifices.

Cd_orifice_H2  = 0.85;
Cd_orifice_Air = 0.85;
Cd_orifice_CH4 = 0.85;
Cd_orifice_GOx = 0.85;

%% PRESSURE MARGINS

% Injector pressure margin:
%
% At the minimum-flow operating point, the manifold
% pressure is this factor above the theoretical minimum
% pressure required to choke the injector into the
% chamber.
%
% Example:
% 1.15 = 15% above theoretical choking pressure.

injector_pressure_margin = 1.15;

% Orifice pressure margin:
%
% At the minimum-flow operating point, the upstream
% pressure is this factor above the theoretical minimum
% pressure required to choke the metering orifice into
% the manifold.

orifice_pressure_margin = 1.15;

%% =====================================================
%  CONVERSIONS
% ======================================================

in_to_m   = 0.0254;
psi_to_Pa = 6894.76;

%% =====================================================
%  CHAMBER PRESSURE
% ======================================================

P_chamber = P_chamber_bar * 1e5;     % [Pa]

%% =====================================================
%  GAS PROPERTIES
% ======================================================

% Approximate gas properties near room temperature.
%
% Update if more accurate values at actual operating
% conditions are available.

% ---------------- H2 ----------------

R_H2     = 4124;        % [J/(kg*K)]
gamma_H2 = 1.405;       % [-]
T0_H2    = 283;         % [K]

% ---------------- AIR ----------------

R_Air     = 287.05;     % [J/(kg*K)]
gamma_Air = 1.400;      % [-]
T0_Air    = 283;        % [K]

% ---------------- CH4 ----------------

R_CH4     = 518.3;      % [J/(kg*K)]
gamma_CH4 = 1.310;      % [-]
T0_CH4    = 283;        % [K]

% ---------------- GOx ----------------

R_GOx     = 259.84;     % [J/(kg*K)]
gamma_GOx = 1.400;      % [-]
T0_GOx    = 283;        % [K]

%% =====================================================
%  READ CSV
% ======================================================

data = readtable( ...
    input_file, ...
    'VariableNamingRule', 'preserve');

%% =====================================================
%  CHECK REQUIRED COLUMNS
% ======================================================

required_columns = { ...
    'H2', ...
    'Air', ...
    'CH4', ...
    'GOx'};

for i = 1:length(required_columns)

    if ~any(strcmpi( ...
            data.Properties.VariableNames, ...
            required_columns{i}))

        error( ...
            'CSV must contain a column named "%s".', ...
            required_columns{i});

    end

end

%% =====================================================
%  FIND COLUMN INDICES
% ======================================================

column_names = data.Properties.VariableNames;

idx_H2 = ...
    find(strcmpi(column_names, 'H2'), 1);

idx_Air = ...
    find(strcmpi(column_names, 'Air'), 1);

idx_CH4 = ...
    find(strcmpi(column_names, 'CH4'), 1);

idx_GOx = ...
    find(strcmpi(column_names, 'GOx'), 1);

%% =====================================================
%  EXTRACT MASS FLOW REQUIREMENTS
% ======================================================

m_dot_H2 = ...
    data{:, idx_H2};

m_dot_Air = ...
    data{:, idx_Air};

m_dot_CH4 = ...
    data{:, idx_CH4};

m_dot_GOx = ...
    data{:, idx_GOx};

%% =====================================================
%  STORE GAS INFORMATION
% ======================================================

gas_names = { ...
    'H2', ...
    'Air', ...
    'CH4', ...
    'GOx'};

m_dot_all = { ...
    m_dot_H2, ...
    m_dot_Air, ...
    m_dot_CH4, ...
    m_dot_GOx};

R_all = [ ...
    R_H2, ...
    R_Air, ...
    R_CH4, ...
    R_GOx];

gamma_all = [ ...
    gamma_H2, ...
    gamma_Air, ...
    gamma_CH4, ...
    gamma_GOx];

T0_all = [ ...
    T0_H2, ...
    T0_Air, ...
    T0_CH4, ...
    T0_GOx];

Cd_inj_all = [ ...
    Cd_inj_H2, ...
    Cd_inj_Air, ...
    Cd_inj_CH4, ...
    Cd_inj_GOx];

Cd_orifice_all = [ ...
    Cd_orifice_H2, ...
    Cd_orifice_Air, ...
    Cd_orifice_CH4, ...
    Cd_orifice_GOx];

%% =====================================================
%  PREALLOCATE RESULTS
% ======================================================

num_rows = height(data);

% Required manifold pressures

P_manifold_bar_all = ...
    nan(num_rows, 4);

P_manifold_psi_all = ...
    nan(num_rows, 4);

% Required supply/orifice upstream pressures

P0_bar_all = ...
    nan(num_rows, 4);

P0_psi_all = ...
    nan(num_rows, 4);

% Choking checks

injector_choked_all = ...
    false(num_rows, 4);

orifice_choked_all = ...
    false(num_rows, 4);

% Injector geometry

injector_area = ...
    nan(1,4);

injector_diameter_m = ...
    nan(1,4);

% Sonic orifice geometry

orifice_area = ...
    nan(1,4);

orifice_diameter_m = ...
    nan(1,4);

% Critical pressures

Pman_critical_all = ...
    nan(1,4);

Pman_design_min_all = ...
    nan(1,4);

P0_orifice_critical_all = ...
    nan(1,4);

P0_orifice_design_min_all = ...
    nan(1,4);

crit_ratio_all = ...
    nan(1,4);

%% =====================================================
%  CALCULATE EACH GAS
% ======================================================

for g = 1:4

    %% -------------------------------------------------
    %  GAS PROPERTIES
    % --------------------------------------------------

    m_dot = ...
        m_dot_all{g};

    R = ...
        R_all(g);

    gamma = ...
        gamma_all(g);

    T0 = ...
        T0_all(g);

    Cd_inj = ...
        Cd_inj_all(g);

    Cd_orifice = ...
        Cd_orifice_all(g);

    %% -------------------------------------------------
    %  VALID FLOW REQUIREMENTS
    % --------------------------------------------------

    % Only positive mass flow values are considered.

    valid = ...
        (m_dot > 0) & ...
        ~isnan(m_dot);

    if ~any(valid)

        warning( ...
            'No positive flow rates found for %s.', ...
            gas_names{g});

        continue;

    end

    %% -------------------------------------------------
    %  MINIMUM MASS FLOW
    % --------------------------------------------------

    % Lowest test-condition mass flow is used to size
    % both the fixed injector area and fixed metering
    % orifice.

    m_dot_min = ...
        min(m_dot(valid));

    %% =================================================
    %  COMMON CHOKED FLOW PARAMETERS
    % ==================================================

    % Critical pressure ratio:
    %
    % P* / P0 =
    %
    % (2/(gamma+1))^(gamma/(gamma-1))

    crit_ratio = ...
        (2 / (gamma + 1)) ^ ...
        (gamma / (gamma - 1));

    crit_ratio_all(g) = ...
        crit_ratio;

    % Choked mass-flow parameter:
    %
    % sqrt(gamma/R)
    %
    % *
    %
    % (2/(gamma+1))^
    % ((gamma+1)/(2*(gamma-1)))

    choked_param = ...
        sqrt(gamma / R) * ...
        (2 / (gamma + 1)) ^ ...
        ((gamma + 1) / ...
        (2 * (gamma - 1)));

    %% =================================================
    %  STAGE 1:
    %
    %  CHAMBER -> INJECTOR -> MANIFOLD
    % ==================================================

    %% MINIMUM MANIFOLD PRESSURE FOR INJECTOR CHOKING

    % Choking condition:
    %
    % P_chamber / P_manifold <= crit_ratio
    %
    % Therefore:
    %
    % P_manifold_critical =
    %
    % P_chamber / crit_ratio

    Pman_critical = ...
        P_chamber / crit_ratio;

    Pman_critical_all(g) = ...
        Pman_critical;

    %% DESIGN MINIMUM MANIFOLD PRESSURE

    % Apply margin above theoretical injector
    % choking boundary.

    Pman_design_min = ...
        injector_pressure_margin * ...
        Pman_critical;

    Pman_design_min_all(g) = ...
        Pman_design_min;

    %% SIZE FIXED INJECTOR EFFECTIVE AREA

    % Choked injector equation:
    %
    % mdot =
    %
    % Cd_inj * A_inj * Pman / sqrt(T0)
    %
    % *
    %
    % choked_param
    %
    % Therefore:
    %
    % A_inj =
    %
    % mdot_min * sqrt(T0)
    % ---------------------------------------
    % Cd_inj * Pman_design_min * choked_param

    A_inj = ...
        (m_dot_min * sqrt(T0)) / ...
        (Cd_inj * ...
        Pman_design_min * ...
        choked_param);

    % Equivalent single circular diameter.
    %
    % NOTE:
    % This is the total effective flow area expressed
    % as one equivalent circular diameter.
    %
    % If multiple injector holes are used, divide the
    % total area among the desired number of holes.

    D_inj = ...
        sqrt(4 * A_inj / pi);

    injector_area(g) = ...
        A_inj;

    injector_diameter_m(g) = ...
        D_inj;

    %% REQUIRED MANIFOLD PRESSURE AT EVERY FLOW POINT

    % Rearranging choked injector equation:
    %
    % Pman =
    %
    % mdot * sqrt(T0)
    % ---------------------------------
    % Cd_inj * A_inj * choked_param

    P_manifold_required = ...
        nan(size(m_dot));

    P_manifold_required(valid) = ...
        (m_dot(valid) .* sqrt(T0)) ./ ...
        (Cd_inj * ...
        A_inj * ...
        choked_param);

    %% CHECK INJECTOR CHOKING

    injector_pressure_ratio = ...
        nan(size(m_dot));

    injector_pressure_ratio(valid) = ...
        P_chamber ./ ...
        P_manifold_required(valid);

    injector_choked = ...
        false(size(m_dot));

    injector_choked(valid) = ...
        injector_pressure_ratio(valid) <= ...
        crit_ratio;

    %% =================================================
    %  STAGE 2:
    %
    %  MANIFOLD -> SONIC ORIFICE -> SUPPLY
    % ==================================================

    % At the minimum-flow condition, the manifold
    % pressure calculated above becomes the downstream
    % pressure of the upstream sonic orifice.

    P_down_orifice_min = ...
        Pman_design_min;

    %% MINIMUM UPSTREAM PRESSURE FOR ORIFICE CHOKING

    % Choking condition:
    %
    % P_manifold / P0 <= crit_ratio
    %
    % Therefore:
    %
    % P0_orifice_critical =
    %
    % P_manifold / crit_ratio

    P0_orifice_critical = ...
        P_down_orifice_min / ...
        crit_ratio;

    P0_orifice_critical_all(g) = ...
        P0_orifice_critical;

    %% DESIGN MINIMUM ORIFICE UPSTREAM PRESSURE

    P0_orifice_design_min = ...
        orifice_pressure_margin * ...
        P0_orifice_critical;

    P0_orifice_design_min_all(g) = ...
        P0_orifice_design_min;

    %% SIZE FIXED SONIC ORIFICE

    % A_orifice =
    %
    % mdot_min * sqrt(T0)
    % ---------------------------------------------
    % Cd_orifice * P0_design_min * choked_param

    A_orifice = ...
        (m_dot_min * sqrt(T0)) / ...
        (Cd_orifice * ...
        P0_orifice_design_min * ...
        choked_param);

    D_orifice = ...
        sqrt(4 * A_orifice / pi);

    orifice_area(g) = ...
        A_orifice;

    orifice_diameter_m(g) = ...
        D_orifice;

    %% REQUIRED UPSTREAM P0 FOR EVERY FLOW POINT

    % Fixed orifice:
    %
    % P0 =
    %
    % mdot * sqrt(T0)
    % --------------------------------------
    % Cd_orifice * A_orifice * choked_param

    P0_required = ...
        nan(size(m_dot));

    P0_required(valid) = ...
        (m_dot(valid) .* sqrt(T0)) ./ ...
        (Cd_orifice * ...
        A_orifice * ...
        choked_param);

    %% CHECK ORIFICE CHOKING

    % IMPORTANT:
    %
    % Downstream pressure of the orifice changes with
    % mass flow because manifold pressure changes with
    % mass flow.
    %
    % Therefore choking is checked row-by-row using
    % each operating point's manifold pressure.

    orifice_pressure_ratio = ...
        nan(size(m_dot));

    orifice_pressure_ratio(valid) = ...
        P_manifold_required(valid) ./ ...
        P0_required(valid);

    orifice_choked = ...
        false(size(m_dot));

    orifice_choked(valid) = ...
        orifice_pressure_ratio(valid) <= ...
        crit_ratio;

    %% -------------------------------------------------
    %  SAVE RESULTS
    % --------------------------------------------------

    P_manifold_bar_all(:,g) = ...
        P_manifold_required / 1e5;

    P_manifold_psi_all(:,g) = ...
        P_manifold_required / psi_to_Pa;

    P0_bar_all(:,g) = ...
        P0_required / 1e5;

    P0_psi_all(:,g) = ...
        P0_required / psi_to_Pa;

    injector_choked_all(:,g) = ...
        injector_choked;

    orifice_choked_all(:,g) = ...
        orifice_choked;

end

%% =====================================================
%  ADD MANIFOLD PRESSURES TO OUTPUT TABLE
% ======================================================

% ---------------- H2 ----------------

data.H2_Manifold_bar = ...
    P_manifold_bar_all(:,1);

data.H2_Manifold_psi = ...
    P_manifold_psi_all(:,1);

data.H2_Injector_Choked = ...
    injector_choked_all(:,1);

% ---------------- AIR ----------------

data.Air_Manifold_bar = ...
    P_manifold_bar_all(:,2);

data.Air_Manifold_psi = ...
    P_manifold_psi_all(:,2);

data.Air_Injector_Choked = ...
    injector_choked_all(:,2);

% ---------------- CH4 ----------------

data.CH4_Manifold_bar = ...
    P_manifold_bar_all(:,3);

data.CH4_Manifold_psi = ...
    P_manifold_psi_all(:,3);

data.CH4_Injector_Choked = ...
    injector_choked_all(:,3);

% ---------------- GOx ----------------

data.GOx_Manifold_bar = ...
    P_manifold_bar_all(:,4);

data.GOx_Manifold_psi = ...
    P_manifold_psi_all(:,4);

data.GOx_Injector_Choked = ...
    injector_choked_all(:,4);

%% =====================================================
%  ADD SUPPLY PRESSURES TO OUTPUT TABLE
% ======================================================

% ---------------- H2 ----------------

data.H2_P0_bar = ...
    P0_bar_all(:,1);

data.H2_P0_psi = ...
    P0_psi_all(:,1);

data.H2_Orifice_Choked = ...
    orifice_choked_all(:,1);

% ---------------- AIR ----------------

data.Air_P0_bar = ...
    P0_bar_all(:,2);

data.Air_P0_psi = ...
    P0_psi_all(:,2);

data.Air_Orifice_Choked = ...
    orifice_choked_all(:,2);

% ---------------- CH4 ----------------

data.CH4_P0_bar = ...
    P0_bar_all(:,3);

data.CH4_P0_psi = ...
    P0_psi_all(:,3);

data.CH4_Orifice_Choked = ...
    orifice_choked_all(:,3);

% ---------------- GOx ----------------

data.GOx_P0_bar = ...
    P0_bar_all(:,4);

data.GOx_P0_psi = ...
    P0_psi_all(:,4);

data.GOx_Orifice_Choked = ...
    orifice_choked_all(:,4);

%% =====================================================
%  ADD INJECTOR EFFECTIVE SIZES TO OUTPUT TABLE
% ======================================================

data.H2_InjectorArea_mm2 = ...
    repmat( ...
    injector_area(1) * 1e6, ...
    num_rows, 1);

data.Air_InjectorArea_mm2 = ...
    repmat( ...
    injector_area(2) * 1e6, ...
    num_rows, 1);

data.CH4_InjectorArea_mm2 = ...
    repmat( ...
    injector_area(3) * 1e6, ...
    num_rows, 1);

data.GOx_InjectorArea_mm2 = ...
    repmat( ...
    injector_area(4) * 1e6, ...
    num_rows, 1);

data.H2_InjectorEqDiameter_mm = ...
    repmat( ...
    injector_diameter_m(1) * 1000, ...
    num_rows, 1);

data.Air_InjectorEqDiameter_mm = ...
    repmat( ...
    injector_diameter_m(2) * 1000, ...
    num_rows, 1);

data.CH4_InjectorEqDiameter_mm = ...
    repmat( ...
    injector_diameter_m(3) * 1000, ...
    num_rows, 1);

data.GOx_InjectorEqDiameter_mm = ...
    repmat( ...
    injector_diameter_m(4) * 1000, ...
    num_rows, 1);

%% =====================================================
%  ADD SONIC ORIFICE SIZES TO OUTPUT TABLE
% ======================================================

data.H2_Orifice_mm = ...
    repmat( ...
    orifice_diameter_m(1) * 1000, ...
    num_rows, 1);

data.Air_Orifice_mm = ...
    repmat( ...
    orifice_diameter_m(2) * 1000, ...
    num_rows, 1);

data.CH4_Orifice_mm = ...
    repmat( ...
    orifice_diameter_m(3) * 1000, ...
    num_rows, 1);

data.GOx_Orifice_mm = ...
    repmat( ...
    orifice_diameter_m(4) * 1000, ...
    num_rows, 1);

data.H2_Orifice_in = ...
    repmat( ...
    orifice_diameter_m(1) / in_to_m, ...
    num_rows, 1);

data.Air_Orifice_in = ...
    repmat( ...
    orifice_diameter_m(2) / in_to_m, ...
    num_rows, 1);

data.CH4_Orifice_in = ...
    repmat( ...
    orifice_diameter_m(3) / in_to_m, ...
    num_rows, 1);

data.GOx_Orifice_in = ...
    repmat( ...
    orifice_diameter_m(4) / in_to_m, ...
    num_rows, 1);

%% =====================================================
%  OUTPUT CSV FILENAME
% ======================================================

[input_path, ...
 input_name, ...
 ~] = ...
    fileparts(input_file);

if isempty(input_path)

    input_path = pwd;

end

output_file = ...
    fullfile( ...
    input_path, ...
    [input_name ...
    '_two_stage_choked_results.csv']);

%% =====================================================
%  WRITE OUTPUT CSV
% ======================================================

writetable( ...
    data, ...
    output_file);

%% =====================================================
%  COMMAND WINDOW SUMMARY
% ======================================================

fprintf('\n');

fprintf( ...
    '==============================================================================================\n');

fprintf( ...
    '                        TWO-STAGE CHOKED FLOW SIZING SUMMARY\n');

fprintf( ...
    '==============================================================================================\n');

fprintf('\n');

fprintf( ...
    'Chamber Pressure              : %.2f bar absolute\n', ...
    P_chamber / 1e5);

fprintf( ...
    'Injector Pressure Margin      : %.3f\n', ...
    injector_pressure_margin);

fprintf( ...
    'Orifice Pressure Margin       : %.3f\n', ...
    orifice_pressure_margin);

fprintf('\n');

fprintf( ...
    ['Gas      Pman,min    Inj Area     Inj Eq D      ', ...
     'Orifice D     P0,min\n']);

fprintf( ...
    ['         [bar]       [mm^2]       [mm]          ', ...
     '[mm]          [bar]\n']);

fprintf( ...
    '----------------------------------------------------------------------------------------------\n');

for g = 1:4

    fprintf( ...
        '%-5s   %9.2f   %10.3f   %10.3f   %11.3f   %10.2f\n', ...
        gas_names{g}, ...
        Pman_design_min_all(g) / 1e5, ...
        injector_area(g) * 1e6, ...
        injector_diameter_m(g) * 1000, ...
        orifice_diameter_m(g) * 1000, ...
        P0_orifice_design_min_all(g) / 1e5);

end

fprintf( ...
    '==============================================================================================\n');

fprintf('\n');

fprintf( ...
    'Results written to:\n%s\n\n', ...
    output_file);