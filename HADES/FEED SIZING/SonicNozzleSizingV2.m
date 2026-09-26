%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Purpose: Fixed Sonic Orifice Sizing from CSV
%
%   Input CSV columns:
%       H2, Air, CH4, GOx
%
%   Each column contains required mass flow rates [kg/s].
%
%   The code sizes ONE fixed sonic orifice for each gas
%   using the minimum positive mass flow requirement.
%
%   Each gas can have its own target downstream pressure.
%
%   It then calculates the upstream stagnation pressure
%   P0 required to achieve every requested mass flow
%   through that fixed orifice.
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

%% DISCHARGE COEFFICIENTS

Cd_H2  = 0.85;
Cd_Air = 0.85;
Cd_CH4 = 0.85;
Cd_GOx = 0.85;

%% PRESSURE MARGIN ABOVE CHOKING BOUNDARY

% Example:
%
% pressure_margin = 1.15 means the minimum-flow
% operating point uses an upstream stagnation pressure
% 15% above the theoretical minimum pressure required
% for choking.

pressure_margin = 1.15;

%% =====================================================
%  CONVERSIONS
% ======================================================

in_to_m   = 0.0254;
psi_to_Pa = 6894.76;

%% =====================================================
%  TARGET DOWNSTREAM PRESSURES
% ======================================================

% Target downstream pressures [bar absolute]
%
% These may be independently adjusted for each gas.

P_down_H2_bar  = 20;
P_down_Air_bar = 20;
P_down_CH4_bar = 24;
P_down_GOx_bar = 36;

%% ADDITIONAL DOWNSTREAM PRESSURE MARGINS

% Additional pressure margin added to each target
% downstream pressure [psi].
%
% Set to zero if no additional margin is desired.

P_margin_H2_psi  = 20;
P_margin_Air_psi = 20;
P_margin_CH4_psi = 20;
P_margin_GOx_psi = 20;

%% CONVERT DOWNSTREAM PRESSURES TO Pa

P_down_H2 = ...
    P_down_H2_bar * 1e5 + ...
    P_margin_H2_psi * psi_to_Pa;

P_down_Air = ...
    P_down_Air_bar * 1e5 + ...
    P_margin_Air_psi * psi_to_Pa;

P_down_CH4 = ...
    P_down_CH4_bar * 1e5 + ...
    P_margin_CH4_psi * psi_to_Pa;

P_down_GOx = ...
    P_down_GOx_bar * 1e5 + ...
    P_margin_GOx_psi * psi_to_Pa;

%% =====================================================
%  GAS PROPERTIES
% ======================================================

% Approximate gas properties near room temperature.
%
% Update these if more accurate values at the actual
% operating conditions are available.

% ---------------- H2 ----------------

R_H2     = 4124;        % Specific gas constant [J/(kg*K)]
gamma_H2 = 1.405;       % Specific heat ratio [-]
T0_H2    = 283;         % Upstream stagnation temperature [K]

% ---------------- AIR ----------------

R_Air     = 287.05;     % Specific gas constant [J/(kg*K)]
gamma_Air = 1.400;      % Specific heat ratio [-]
T0_Air    = 283;        % Upstream stagnation temperature [K]

% ---------------- CH4 ----------------

R_CH4     = 518.3;      % Specific gas constant [J/(kg*K)]
gamma_CH4 = 1.310;      % Specific heat ratio [-]
T0_CH4    = 283;        % Upstream stagnation temperature [K]

% ---------------- GOx ----------------

R_GOx     = 259.84;     % Specific gas constant [J/(kg*K)]
gamma_GOx = 1.400;      % Specific heat ratio [-]
T0_GOx    = 283;        % Upstream stagnation temperature [K]

%% =====================================================
%  READ CSV
% ======================================================

data = readtable(input_file, ...
    'VariableNamingRule', 'preserve');

%% =====================================================
%  CHECK REQUIRED COLUMNS
% ======================================================

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

Cd_all = [ ...
    Cd_H2, ...
    Cd_Air, ...
    Cd_CH4, ...
    Cd_GOx];

P_down_all = [ ...
    P_down_H2, ...
    P_down_Air, ...
    P_down_CH4, ...
    P_down_GOx];

%% =====================================================
%  PREALLOCATE RESULTS
% ======================================================

num_rows = height(data);

P0_bar_all = ...
    nan(num_rows, 4);

P0_psi_all = ...
    nan(num_rows, 4);

choked_all = ...
    false(num_rows, 4);

orifice_area = ...
    nan(1, 4);

orifice_diameter_m = ...
    nan(1, 4);

P0_critical_all = ...
    nan(1, 4);

P0_design_min_all = ...
    nan(1, 4);

crit_ratio_all = ...
    nan(1, 4);

%% =====================================================
%  SIZE EACH GAS ORIFICE
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

    Cd = ...
        Cd_all(g);

    P_down = ...
        P_down_all(g);

    %% -------------------------------------------------
    %  VALID FLOW REQUIREMENTS
    % --------------------------------------------------

    % Only positive flow rates are used.
    %
    % Zero and NaN entries are ignored.

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

    % The minimum positive mass flow is used to size
    % the fixed orifice.

    m_dot_min = ...
        min(m_dot(valid));

    %% -------------------------------------------------
    %  CRITICAL PRESSURE RATIO
    % --------------------------------------------------

    % At M = 1:
    %
    % P* / P0 =
    %
    % (2/(gamma+1))^(gamma/(gamma-1))

    crit_ratio = ...
        (2 / (gamma + 1)) ^ ...
        (gamma / (gamma - 1));

    crit_ratio_all(g) = ...
        crit_ratio;

    %% -------------------------------------------------
    %  MINIMUM PRESSURE FOR CHOKING
    % --------------------------------------------------

    % Choking requirement:
    %
    % P_down / P0 <= critical pressure ratio
    %
    % Therefore:
    %
    % P0_critical =
    %
    % P_down / critical pressure ratio

    P0_critical = ...
        P_down / crit_ratio;

    P0_critical_all(g) = ...
        P0_critical;

    %% -------------------------------------------------
    %  DESIGN MINIMUM UPSTREAM PRESSURE
    % --------------------------------------------------

    % Apply pressure margin above theoretical choking
    % boundary.

    P0_design_min = ...
        pressure_margin * ...
        P0_critical;

    P0_design_min_all(g) = ...
        P0_design_min;

    %% -------------------------------------------------
    %  CHOKED FLOW PARAMETER
    % --------------------------------------------------

    % Choked flow equation:
    %
    % mdot =
    %
    % Cd * At * P0 / sqrt(T0)
    %
    % *
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

    %% -------------------------------------------------
    %  SIZE ONE FIXED ORIFICE
    % --------------------------------------------------

    % Rearranging:
    %
    % At =
    %
    % mdot * sqrt(T0)
    % ------------------------------
    % Cd * P0 * choked_param

    A_t = ...
        (m_dot_min * sqrt(T0)) / ...
        (Cd * ...
        P0_design_min * ...
        choked_param);

    %% ORIFICE DIAMETER

    D_t = ...
        sqrt(4 * A_t / pi);

    %% SAVE FIXED ORIFICE RESULTS

    orifice_area(g) = ...
        A_t;

    orifice_diameter_m(g) = ...
        D_t;

    %% -------------------------------------------------
    %  REQUIRED P0 FOR EVERY FLOW POINT
    % --------------------------------------------------

    % Rearranging the choked-flow equation:
    %
    % P0 =
    %
    % mdot * sqrt(T0)
    % -----------------------------
    % Cd * At * choked_param

    P0_required = ...
        nan(size(m_dot));

    P0_required(valid) = ...
        (m_dot(valid) .* sqrt(T0)) ./ ...
        (Cd * ...
        A_t * ...
        choked_param);

    %% -------------------------------------------------
    %  CHECK CHOKING
    % --------------------------------------------------

    % Choking occurs if:
    %
    % P_down / P0 <= critical pressure ratio

    pressure_ratio = ...
        nan(size(m_dot));

    pressure_ratio(valid) = ...
        P_down ./ ...
        P0_required(valid);

    is_choked = ...
        false(size(m_dot));

    is_choked(valid) = ...
        pressure_ratio(valid) <= ...
        crit_ratio;

    %% -------------------------------------------------
    %  SAVE RESULTS
    % --------------------------------------------------

    P0_bar_all(:, g) = ...
        P0_required / 1e5;

    P0_psi_all(:, g) = ...
        P0_required / psi_to_Pa;

    choked_all(:, g) = ...
        is_choked;

end

%% =====================================================
%  ADD RESULTS TO OUTPUT TABLE
% ======================================================

% ---------------- H2 ----------------

data.H2_P0_bar = ...
    P0_bar_all(:,1);

data.H2_P0_psi = ...
    P0_psi_all(:,1);

data.H2_Choked = ...
    choked_all(:,1);

% ---------------- AIR ----------------

data.Air_P0_bar = ...
    P0_bar_all(:,2);

data.Air_P0_psi = ...
    P0_psi_all(:,2);

data.Air_Choked = ...
    choked_all(:,2);

% ---------------- CH4 ----------------

data.CH4_P0_bar = ...
    P0_bar_all(:,3);

data.CH4_P0_psi = ...
    P0_psi_all(:,3);

data.CH4_Choked = ...
    choked_all(:,3);

% ---------------- GOx ----------------

data.GOx_P0_bar = ...
    P0_bar_all(:,4);

data.GOx_P0_psi = ...
    P0_psi_all(:,4);

data.GOx_Choked = ...
    choked_all(:,4);

%% =====================================================
%  ADD DOWNSTREAM PRESSURES TO OUTPUT TABLE
% ======================================================

data.H2_Pdown_bar = ...
    repmat( ...
    P_down_H2 / 1e5, ...
    num_rows, 1);

data.Air_Pdown_bar = ...
    repmat( ...
    P_down_Air / 1e5, ...
    num_rows, 1);

data.CH4_Pdown_bar = ...
    repmat( ...
    P_down_CH4 / 1e5, ...
    num_rows, 1);

data.GOx_Pdown_bar = ...
    repmat( ...
    P_down_GOx / 1e5, ...
    num_rows, 1);

%% =====================================================
%  ADD FIXED ORIFICE SIZES TO OUTPUT TABLE
% ======================================================

% These values are constant for every row because
% each gas uses ONE fixed orifice.

% ---------------- H2 ----------------

data.H2_Orifice_mm = ...
    repmat( ...
    orifice_diameter_m(1) * 1000, ...
    num_rows, 1);

data.H2_Orifice_in = ...
    repmat( ...
    orifice_diameter_m(1) / in_to_m, ...
    num_rows, 1);

% ---------------- AIR ----------------

data.Air_Orifice_mm = ...
    repmat( ...
    orifice_diameter_m(2) * 1000, ...
    num_rows, 1);

data.Air_Orifice_in = ...
    repmat( ...
    orifice_diameter_m(2) / in_to_m, ...
    num_rows, 1);

% ---------------- CH4 ----------------

data.CH4_Orifice_mm = ...
    repmat( ...
    orifice_diameter_m(3) * 1000, ...
    num_rows, 1);

data.CH4_Orifice_in = ...
    repmat( ...
    orifice_diameter_m(3) / in_to_m, ...
    num_rows, 1);

% ---------------- GOx ----------------

data.GOx_Orifice_mm = ...
    repmat( ...
    orifice_diameter_m(4) * 1000, ...
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

    input_path = ...
        pwd;

end

output_file = ...
    fullfile( ...
    input_path, ...
    [input_name ...
    '_sonic_orifice_results.csv']);

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
    '===============================================================================\n');

fprintf( ...
    '                   FIXED SONIC ORIFICE SIZING SUMMARY\n');

fprintf( ...
    '===============================================================================\n');

fprintf('\n');

fprintf( ...
    'Pressure Margin Factor : %.3f\n', ...
    pressure_margin);

fprintf('\n');

fprintf( ...
    ['Gas       Cd      Pdown [bar]     D [mm]       D [in]      ', ...
    'P0crit [bar]   P0,min [bar]\n']);

fprintf( ...
    '-------------------------------------------------------------------------------\n');

for g = 1:4

    fprintf( ...
        '%-5s   %6.3f     %9.2f     %9.3f     %9.4f     %11.2f     %11.2f\n', ...
        gas_names{g}, ...
        Cd_all(g), ...
        P_down_all(g) / 1e5, ...
        orifice_diameter_m(g) * 1000, ...
        orifice_diameter_m(g) / in_to_m, ...
        P0_critical_all(g) / 1e5, ...
        P0_design_min_all(g) / 1e5);

end

fprintf( ...
    '===============================================================================\n');

fprintf('\n');

fprintf( ...
    'Results written to:\n%s\n\n', ...
    output_file);