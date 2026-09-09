%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Purpose: Sonic Nozzle Sizing from Downstream
%            Pressure and Mass Flow Requirements
%
%   Programmers: Deepesh Balwani, Noah Ha
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;

%% CONVERSIONS

in_to_m   = 0.0254;
lb_to_kg  = 0.453592;
psi_to_Pa = 6894.76;

%% DESIGN REQUIREMENTS

% Mass flow requirement
m_dot_req   = 2 * lb_to_kg;      % Required mass flow [kg/s]
flow_margin = 1.15;              % 15% design margin
m_dot       = m_dot_req * flow_margin;

% Downstream pressure requirement
P_down = 20e5;                   % Downstream/back pressure [Pa] = 20 bar

% Air properties
R       = 287;                   % Specific gas constant [J/(kg*K)]
gamma   = 1.4;                   % Specific heat ratio [-]
T_tu    = 283;                   % Upstream total temperature [K]

%% CHOKING CONDITION

% Critical pressure ratio for choking:
%
% P_t / P_tu = (2/(gamma+1))^(gamma/(gamma-1))
%
% For air, this is approximately 0.5283.

crit_ratio = (2/(gamma + 1))^(gamma/(gamma - 1));

% Minimum upstream total pressure required for choking.
%
% At the limiting choking condition:
%
% P_down = P_t = crit_ratio * P_tu
%
% Therefore:
%
% P_tu_min = P_down / crit_ratio

P_tu_min = P_down / crit_ratio;

%% CHOKED FLOW FUNCTION

% Choked mass-flow parameter:
%
% mdot = (P_tu * A_t / sqrt(T_tu)) *
%        sqrt(gamma/R) *
%        (2/(gamma+1))^((gamma+1)/(2*(gamma-1)))

choked_param = sqrt(gamma / R) * ...
               (2 / (gamma + 1))^((gamma + 1) / ...
               (2 * (gamma - 1)));

%% THROAT AREA

% Rearranging the choked-flow equation:
%
% A_t = mdot * sqrt(T_tu) / (P_tu * choked_param)

A_t = (m_dot * sqrt(T_tu)) / ...
      (P_tu_min * choked_param);

%% THROAT DIAMETER

D_t_m = 2 * sqrt(A_t / pi);

%% THROAT CONDITIONS

% At M = 1:
%
% T_t / T_tu = 2/(gamma+1)
% P_t / P_tu = crit_ratio
%
% Calculate throat static conditions.

T_throat = T_tu * (2 / (gamma + 1));
P_throat = P_tu_min * crit_ratio;

%% FEED LINE GEOMETRY

wall_thickness = 0.049 * in_to_m;

feed_OD = (16 / 16) * in_to_m;       % 1.000 in OD
feed_ID = feed_OD - 2 * wall_thickness;

%% OUTPUTS

fprintf('\n');
fprintf('=============================================\n');
fprintf('       SONIC NOZZLE SIZING RESULTS\n');
fprintf('=============================================\n');

fprintf('\n--- DESIGN REQUIREMENTS ---\n');

fprintf('Required Mass Flow          : %.3f lb/s\n', ...
        m_dot_req / lb_to_kg);

fprintf('Design Mass Flow            : %.3f lb/s\n', ...
        m_dot / lb_to_kg);

fprintf('Design Mass Flow            : %.3f kg/s\n', ...
        m_dot);

fprintf('Downstream Pressure         : %.2f psi\n', ...
        P_down / psi_to_Pa);

fprintf('Downstream Pressure         : %.2f bar\n', ...
        P_down / 1e5);

fprintf('Upstream Total Temperature  : %.2f K\n', ...
        T_tu);

fprintf('\n--- CHOKING CONDITION ---\n');

fprintf('Critical Pressure Ratio     : %.4f\n', ...
        crit_ratio);

fprintf('Minimum Upstream Total P    : %.2f psi\n', ...
        P_tu_min / psi_to_Pa);

fprintf('Minimum Upstream Total P    : %.2f bar\n', ...
        P_tu_min / 1e5);

fprintf('\n--- THROAT ---\n');

fprintf('Throat Area                 : %.6f m^2\n', ...
        A_t);

fprintf('Throat Area                 : %.4f in^2\n', ...
        A_t / in_to_m^2);

fprintf('Throat Diameter             : %.2f thou\n', ...
        (D_t_m / in_to_m) * 1000);

fprintf('Throat Diameter             : %.4f in\n', ...
        D_t_m / in_to_m);

fprintf('Throat Diameter             : %.3f mm\n', ...
        D_t_m * 1000);

fprintf('\n--- THROAT STATIC CONDITIONS ---\n');

fprintf('Throat Pressure (M=1)       : %.2f psi\n', ...
        P_throat / psi_to_Pa);

fprintf('Throat Pressure (M=1)       : %.2f bar\n', ...
        P_throat / 1e5);

fprintf('Throat Temperature (M=1)    : %.2f K\n', ...
        T_throat);

fprintf('\n--- FEED LINE ---\n');

fprintf('Feed OD                     : %.3f in\n', ...
        feed_OD / in_to_m);

fprintf('Feed ID                     : %.3f in\n', ...
        feed_ID / in_to_m);

fprintf('=============================================\n');