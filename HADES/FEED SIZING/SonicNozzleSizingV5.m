%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Purpose:
%
%   Combined RDE Injector + Fixed Sonic Orifice Model
%
%   1. Read desired mass flow rates from CSV
%
%   2. Calculate required manifold pressure from:
%          - Injector area
%          - Injector Cd
%          - Desired mass flow
%
%   3. Check whether injector is choked
%
%   4. Size ONE fixed sonic orifice for each gas
%
%   5. Calculate required upstream stagnation pressure
%      for that fixed sonic orifice at every test point
%
%   Input CSV columns:
%
%       H2, Air, CH4, GOx
%
%   Mass flow units:
%
%       kg/s
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

P_chamber_bar = 10.0;

% Atmospheric pressure used only for converting
% absolute supply pressure -> gauge pressure

P_atm = 101325;       % [Pa]

%% =====================================================
%  PRESSURE MARGIN FOR SONIC ORIFICE
% ======================================================

% Sonic orifice will be sized at the minimum requested
% mass flow such that its upstream pressure is this
% factor above the theoretical choking boundary.
%
% Example:
%
% 1.15 = 15% pressure margin

pressure_margin = 1.15;

%% =====================================================
%  TOTAL INJECTOR AREAS
% ======================================================

% Total physical injector flow areas [mm^2]

A_H2_mm2  = 25.13;
A_Air_mm2 = 304.77;
A_CH4_mm2 = 56.55;
A_GOx_mm2 = 86.59;

% Convert mm^2 -> m^2

A_H2 = ...
    A_H2_mm2 * 1e-6;

A_Air = ...
    A_Air_mm2 * 1e-6;

A_CH4 = ...
    A_CH4_mm2 * 1e-6;

A_GOx = ...
    A_GOx_mm2 * 1e-6;

%% =====================================================
%  EFFECTIVE INJECTOR DISCHARGE COEFFICIENTS
% ======================================================

Cd_inj_H2  = 0.71;
Cd_inj_Air = 0.62;
Cd_inj_CH4 = 0.60;
Cd_inj_GOx = 0.60;

%% =====================================================
%  SONIC ORIFICE DISCHARGE COEFFICIENTS
% ======================================================

% These are separate from the injector Cd values.
%
% Replace with experimentally determined sonic-orifice
% discharge coefficients when available.

Cd_so_H2  = 0.85;
Cd_so_Air = 0.85;
Cd_so_CH4 = 0.85;
Cd_so_GOx = 0.85;

%% =====================================================
%  GAS PROPERTIES
% ======================================================

% Stagnation temperature [K]

T0_H2  = 283;
T0_Air = 283;
T0_CH4 = 283;
T0_GOx = 283;

% Specific gas constant [J/(kg*K)]

R_H2  = 4124.0;
R_Air = 287.05;
R_CH4 = 518.3;
R_GOx = 259.84;

% Specific heat ratio [-]

gamma_H2  = 1.405;
gamma_Air = 1.400;
gamma_CH4 = 1.310;
gamma_GOx = 1.400;

%% =====================================================
%  CONVERSIONS
% ======================================================

in_to_m = ...
    0.0254;

psi_to_Pa = ...
    6894.757293;

%% =====================================================
%  READ INPUT CSV
% ======================================================

data = readtable( ...
    input_file, ...
    'VariableNamingRule', ...
    'preserve');

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

column_names = ...
    data.Properties.VariableNames;

idx_H2 = ...
    find(strcmpi( ...
    column_names, ...
    'H2'), ...
    1);

idx_Air = ...
    find(strcmpi( ...
    column_names, ...
    'Air'), ...
    1);

idx_CH4 = ...
    find(strcmpi( ...
    column_names, ...
    'CH4'), ...
    1);

idx_GOx = ...
    find(strcmpi( ...
    column_names, ...
    'GOx'), ...
    1);

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

A_inj_all = [ ...
    A_H2, ...
    A_Air, ...
    A_CH4, ...
    A_GOx];

Cd_inj_all = [ ...
    Cd_inj_H2, ...
    Cd_inj_Air, ...
    Cd_inj_CH4, ...
    Cd_inj_GOx];

Cd_so_all = [ ...
    Cd_so_H2, ...
    Cd_so_Air, ...
    Cd_so_CH4, ...
    Cd_so_GOx];

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

%% =====================================================
%  PREALLOCATE RESULTS
% ======================================================

num_rows = ...
    height(data);

Pman_bar_all = ...
    nan(num_rows, 4);

injector_choked_all = ...
    false(num_rows, 4);

P0_supply_bar_all = ...
    nan(num_rows, 4);

P0_supply_psia_all = ...
    nan(num_rows, 4);

P0_supply_psig_all = ...
    nan(num_rows, 4);

sonic_choked_all = ...
    false(num_rows, 4);

model_valid_all = ...
    false(num_rows, 4);

sonic_area_all = ...
    nan(1, 4);

sonic_diameter_all = ...
    nan(1, 4);

critical_ratio_all = ...
    nan(1, 4);

min_flow_all = ...
    nan(1, 4);

Pman_minflow_all = ...
    nan(1, 4);

P0_critical_minflow_all = ...
    nan(1, 4);

P0_design_minflow_all = ...
    nan(1, 4);

%% =====================================================
%  LOOP THROUGH EACH GAS
% ======================================================

for g = 1:4

    %% -------------------------------------------------
    %  GAS DATA
    % --------------------------------------------------

    m_dot = ...
        m_dot_all{g};

    A_inj = ...
        A_inj_all(g);

    Cd_inj = ...
        Cd_inj_all(g);

    Cd_so = ...
        Cd_so_all(g);

    R = ...
        R_all(g);

    gamma = ...
        gamma_all(g);

    T0 = ...
        T0_all(g);

    %% -------------------------------------------------
    %  VALID MASS FLOW POINTS
    % --------------------------------------------------

    valid = ...
        (m_dot > 0) & ...
        ~isnan(m_dot);

    if ~any(valid)

        warning( ...
            'No valid positive mass flow rates for %s.', ...
            gas_names{g});

        continue;

    end

    %% -------------------------------------------------
    %  CHOKED FLOW PARAMETER
    % --------------------------------------------------

    % K =
    %
    % sqrt(gamma/R)
    %
    % *
    %
    % (2/(gamma+1))^
    %
    % ((gamma+1)/(2*(gamma-1)))

    K = ...
        sqrt(gamma / R) * ...
        (2 / (gamma + 1)) ^ ...
        ((gamma + 1) / ...
        (2 * (gamma - 1)));

    %% -------------------------------------------------
    %  CRITICAL PRESSURE RATIO
    % --------------------------------------------------

    % P* / P0

    critical_ratio = ...
        (2 / (gamma + 1)) ^ ...
        (gamma / ...
        (gamma - 1));

    critical_ratio_all(g) = ...
        critical_ratio;

    %% =================================================
    %  STEP 1:
    %  CALCULATE REQUIRED MANIFOLD PRESSURE
    % ==================================================

    % Choked injector equation:
    %
    % mdot =
    %
    % Cd_inj * A_inj * Pman / sqrt(T0) * K
    %
    % Therefore:
    %
    % Pman =
    %
    % mdot * sqrt(T0)
    % -------------------------
    % Cd_inj * A_inj * K

    Pman = ...
        nan(size(m_dot));

    Pman(valid) = ...
        (m_dot(valid) .* sqrt(T0)) ./ ...
        (Cd_inj * ...
        A_inj * ...
        K);

    Pman_bar = ...
        Pman / 1e5;

    %% -------------------------------------------------
    %  CHECK INJECTOR CHOKING
    % --------------------------------------------------

    % Injector choking condition:
    %
    % P_chamber / Pman <= critical_ratio

    injector_pressure_ratio = ...
        nan(size(m_dot));

    injector_pressure_ratio(valid) = ...
        P_chamber_bar ./ ...
        Pman_bar(valid);

    injector_choked = ...
        false(size(m_dot));

    injector_choked(valid) = ...
        injector_pressure_ratio(valid) <= ...
        critical_ratio;

    %% -------------------------------------------------
    %  WARN IF INJECTOR IS NOT CHOKED
    % --------------------------------------------------

    bad_injector = ...
        valid & ...
        ~injector_choked;

    if any(bad_injector)

        warning( ...
            ['%s injector is NOT choked at %d test point(s). ' ...
            'The choked-injector manifold pressure equation ' ...
            'is not valid at those points.'], ...
            gas_names{g}, ...
            sum(bad_injector));

    end

    %% =================================================
    %  STEP 2:
    %  FIND MINIMUM MASS FLOW POINT
    % ==================================================

    valid_indices = ...
        find(valid);

    [m_dot_min, ...
     local_min_index] = ...
        min(m_dot(valid));

    min_index = ...
        valid_indices(local_min_index);

    Pman_min = ...
        Pman(min_index);

    min_flow_all(g) = ...
        m_dot_min;

    Pman_minflow_all(g) = ...
        Pman_min;

    %% =================================================
    %  STEP 3:
    %  SONIC ORIFICE MINIMUM CHOKING PRESSURE
    % ==================================================

    % Sonic orifice downstream pressure is the manifold
    % pressure.
    %
    % Choking condition:
    %
    % Pman / P0 <= critical_ratio
    %
    % Therefore:
    %
    % P0_critical =
    %
    % Pman / critical_ratio

    P0_critical_min = ...
        Pman_min / ...
        critical_ratio;

    P0_critical_minflow_all(g) = ...
        P0_critical_min;

    %% =================================================
    %  STEP 4:
    %  APPLY PRESSURE MARGIN
    % ==================================================

    P0_design_min = ...
        pressure_margin * ...
        P0_critical_min;

    P0_design_minflow_all(g) = ...
        P0_design_min;

    %% =================================================
    %  STEP 5:
    %  SIZE ONE FIXED SONIC ORIFICE
    % ==================================================

    % Choked sonic-orifice equation:
    %
    % mdot =
    %
    % Cd_so * A_so * P0 / sqrt(T0) * K
    %
    %
    % Therefore:
    %
    % A_so =
    %
    % mdot * sqrt(T0)
    % -------------------------
    % Cd_so * P0 * K

    A_so = ...
        (m_dot_min * sqrt(T0)) / ...
        (Cd_so * ...
        P0_design_min * ...
        K);

    D_so = ...
        sqrt(4 * A_so / pi);

    sonic_area_all(g) = ...
        A_so;

    sonic_diameter_all(g) = ...
        D_so;

    %% =================================================
    %  STEP 6:
    %  REQUIRED SONIC-ORIFICE UPSTREAM PRESSURE
    % ==================================================

    % Fixed A_so for every test point.
    %
    % Rearranging:
    %
    % P0 =
    %
    % mdot * sqrt(T0)
    % -------------------------
    % Cd_so * A_so * K

    P0_supply = ...
        nan(size(m_dot));

    P0_supply(valid) = ...
        (m_dot(valid) .* sqrt(T0)) ./ ...
        (Cd_so * ...
        A_so * ...
        K);

    %% -------------------------------------------------
    %  CONVERT SUPPLY PRESSURES
    % --------------------------------------------------

    P0_supply_bar = ...
        P0_supply / 1e5;

    P0_supply_psia = ...
        P0_supply / ...
        psi_to_Pa;

    P0_supply_psig = ...
        (P0_supply - P_atm) / ...
        psi_to_Pa;

    %% =================================================
    %  STEP 7:
    %  CHECK SONIC ORIFICE CHOKING
    % ==================================================

    % Sonic-orifice choking condition:
    %
    % Pman / P0_supply <= critical_ratio

    sonic_pressure_ratio = ...
        nan(size(m_dot));

    sonic_pressure_ratio(valid) = ...
        Pman(valid) ./ ...
        P0_supply(valid);

    sonic_choked = ...
        false(size(m_dot));

    sonic_choked(valid) = ...
        sonic_pressure_ratio(valid) <= ...
        critical_ratio;

    %% =================================================
    %  STEP 8:
    %  OVERALL MODEL VALIDITY
    % ==================================================

    % This simplified model is valid only when BOTH:
    %
    % 1. Injector is choked
    % 2. Sonic orifice is choked

    model_valid = ...
        injector_choked & ...
        sonic_choked;

    %% -------------------------------------------------
    %  SAVE RESULTS
    % --------------------------------------------------

    Pman_bar_all(:, g) = ...
        Pman_bar;

    injector_choked_all(:, g) = ...
        injector_choked;

    P0_supply_bar_all(:, g) = ...
        P0_supply_bar;

    P0_supply_psia_all(:, g) = ...
        P0_supply_psia;

    P0_supply_psig_all(:, g) = ...
        P0_supply_psig;

    sonic_choked_all(:, g) = ...
        sonic_choked;

    model_valid_all(:, g) = ...
        model_valid;

end

%% =====================================================
%  ADD H2 RESULTS TO OUTPUT TABLE
% ======================================================

data.H2_ManifoldPressure_bar = ...
    Pman_bar_all(:,1);

data.H2_InjectorChoked = ...
    injector_choked_all(:,1);

data.H2_SupplyP0_bar_abs = ...
    P0_supply_bar_all(:,1);

data.H2_SupplyP0_psia = ...
    P0_supply_psia_all(:,1);

data.H2_SupplyP0_psig = ...
    P0_supply_psig_all(:,1);

data.H2_SonicOrificeChoked = ...
    sonic_choked_all(:,1);

data.H2_ModelValid = ...
    model_valid_all(:,1);

data.H2_SonicOrificeArea_mm2 = ...
    repmat( ...
    sonic_area_all(1) * 1e6, ...
    num_rows, ...
    1);

data.H2_SonicOrificeDiameter_mm = ...
    repmat( ...
    sonic_diameter_all(1) * 1000, ...
    num_rows, ...
    1);

data.H2_SonicOrificeDiameter_in = ...
    repmat( ...
    sonic_diameter_all(1) / in_to_m, ...
    num_rows, ...
    1);

%% =====================================================
%  ADD AIR RESULTS TO OUTPUT TABLE
% ======================================================

data.Air_ManifoldPressure_bar = ...
    Pman_bar_all(:,2);

data.Air_InjectorChoked = ...
    injector_choked_all(:,2);

data.Air_SupplyP0_bar_abs = ...
    P0_supply_bar_all(:,2);

data.Air_SupplyP0_psia = ...
    P0_supply_psia_all(:,2);

data.Air_SupplyP0_psig = ...
    P0_supply_psig_all(:,2);

data.Air_SonicOrificeChoked = ...
    sonic_choked_all(:,2);

data.Air_ModelValid = ...
    model_valid_all(:,2);

data.Air_SonicOrificeArea_mm2 = ...
    repmat( ...
    sonic_area_all(2) * 1e6, ...
    num_rows, ...
    1);

data.Air_SonicOrificeDiameter_mm = ...
    repmat( ...
    sonic_diameter_all(2) * 1000, ...
    num_rows, ...
    1);

data.Air_SonicOrificeDiameter_in = ...
    repmat( ...
    sonic_diameter_all(2) / in_to_m, ...
    num_rows, ...
    1);

%% =====================================================
%  ADD CH4 RESULTS TO OUTPUT TABLE
% ======================================================

data.CH4_ManifoldPressure_bar = ...
    Pman_bar_all(:,3);

data.CH4_InjectorChoked = ...
    injector_choked_all(:,3);

data.CH4_SupplyP0_bar_abs = ...
    P0_supply_bar_all(:,3);

data.CH4_SupplyP0_psia = ...
    P0_supply_psia_all(:,3);

data.CH4_SupplyP0_psig = ...
    P0_supply_psig_all(:,3);

data.CH4_SonicOrificeChoked = ...
    sonic_choked_all(:,3);

data.CH4_ModelValid = ...
    model_valid_all(:,3);

data.CH4_SonicOrificeArea_mm2 = ...
    repmat( ...
    sonic_area_all(3) * 1e6, ...
    num_rows, ...
    1);

data.CH4_SonicOrificeDiameter_mm = ...
    repmat( ...
    sonic_diameter_all(3) * 1000, ...
    num_rows, ...
    1);

data.CH4_SonicOrificeDiameter_in = ...
    repmat( ...
    sonic_diameter_all(3) / in_to_m, ...
    num_rows, ...
    1);

%% =====================================================
%  ADD GOx RESULTS TO OUTPUT TABLE
% ======================================================

data.GOx_ManifoldPressure_bar = ...
    Pman_bar_all(:,4);

data.GOx_InjectorChoked = ...
    injector_choked_all(:,4);

data.GOx_SupplyP0_bar_abs = ...
    P0_supply_bar_all(:,4);

data.GOx_SupplyP0_psia = ...
    P0_supply_psia_all(:,4);

data.GOx_SupplyP0_psig = ...
    P0_supply_psig_all(:,4);

data.GOx_SonicOrificeChoked = ...
    sonic_choked_all(:,4);

data.GOx_ModelValid = ...
    model_valid_all(:,4);

data.GOx_SonicOrificeArea_mm2 = ...
    repmat( ...
    sonic_area_all(4) * 1e6, ...
    num_rows, ...
    1);

data.GOx_SonicOrificeDiameter_mm = ...
    repmat( ...
    sonic_diameter_all(4) * 1000, ...
    num_rows, ...
    1);

data.GOx_SonicOrificeDiameter_in = ...
    repmat( ...
    sonic_diameter_all(4) / in_to_m, ...
    num_rows, ...
    1);

%% =====================================================
%  OUTPUT CSV FILENAME
% ======================================================

[input_path, ...
 input_name, ...
 ~] = ...
    fileparts(input_file);

if isempty(input_path)

    input_path = ...
        pwd;

end

output_file = ...
    fullfile( ...
    input_path, ...
    [input_name ...
    '_combined_results.csv']);

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
    '========================================================================================\n');

fprintf( ...
    '             RDE INJECTOR + FIXED SONIC ORIFICE SIZING SUMMARY\n');

fprintf( ...
    '========================================================================================\n');

fprintf('\n');

fprintf( ...
    'Chamber Pressure             : %.3f bar absolute\n', ...
    P_chamber_bar);

fprintf( ...
    'Sonic Orifice Pressure Margin: %.3f\n', ...
    pressure_margin);

fprintf('\n');

fprintf( ...
    ['Gas      Cd_inj   Cd_SO    Min Flow     Pman@Min    ', ...
    'SO D [mm]   SO D [in]   P0@Min [bar]\n']);

fprintf( ...
    ['                           [kg/s]        [bar abs]    ', ...
    '                        [bar abs]\n']);

fprintf( ...
    '----------------------------------------------------------------------------------------\n');

for g = 1:4

    fprintf( ...
        '%-5s    %6.3f   %6.3f    %8.4f     %8.3f     %8.3f    %8.4f     %8.3f\n', ...
        gas_names{g}, ...
        Cd_inj_all(g), ...
        Cd_so_all(g), ...
        min_flow_all(g), ...
        Pman_minflow_all(g) / 1e5, ...
        sonic_diameter_all(g) * 1000, ...
        sonic_diameter_all(g) / in_to_m, ...
        P0_design_minflow_all(g) / 1e5);

end

fprintf( ...
    '========================================================================================\n');

fprintf('\n');

%% =====================================================
%  CHOKING SUMMARY
% ======================================================

for g = 1:4

    valid = ...
        (m_dot_all{g} > 0) & ...
        ~isnan(m_dot_all{g});

    fprintf( ...
        '%s:\n', ...
        gas_names{g});

    fprintf( ...
        '  Injector choked at %d of %d valid points\n', ...
        sum(injector_choked_all(valid,g)), ...
        sum(valid));

    fprintf( ...
        '  Sonic orifice choked at %d of %d valid points\n', ...
        sum(sonic_choked_all(valid,g)), ...
        sum(valid));

    fprintf( ...
        '  Full simplified model valid at %d of %d valid points\n', ...
        sum(model_valid_all(valid,g)), ...
        sum(valid));

    fprintf('\n');

end

fprintf( ...
    'Results written to:\n%s\n\n', ...
    output_file);