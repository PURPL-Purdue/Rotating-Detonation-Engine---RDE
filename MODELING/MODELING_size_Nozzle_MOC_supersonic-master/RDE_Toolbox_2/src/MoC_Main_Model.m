%% MoC_Main_Model - starting framework for the hand-written RDE model
% This is a script: inputs and ceaOut remain visible in the MATLAB workspace.
% CEA results are passed to the triple-point resolver below.

clear
clc

%% Paths - relative to this file, not MATLAB's current folder
modelDir = fileparts(mfilename('fullpath'));
toolboxDir = fileparts(modelDir);
resultsDir = fullfile(toolboxDir, 'results');
ceaResultsDir = fullfile(resultsDir, 'CEA');

addpath(modelDir);
% Plotting placeholders now live in src, so no post_processing path is needed.

%% CEA Chemistry inputs
ox_type = 'O2';
fuel_type = 'CH4';

% Choose one: 'of', 'pctFuel', 'phi', or 'r'.
mixture_spec = 'phi';
mixture_spec_value = 1.3;

% These are the initial UNBURNED reactant conditions supplied to CEA.
P0Units = 'bar';
P0 = 2.5;
T0Units = 'K';
T0 = 283;

%% Run CEA 
fprintf('\n--- RUNNING CEA ---\n');
ceaPath = getCEAPath();

% CEA inputs:
ceaOut = HADES_size_ceaDet( ...
    'ox', ox_type, ...
    'fuel', fuel_type, ...
    mixture_spec, mixture_spec_value, ...
    'P0', P0, ...
    'P0Units', P0Units, ...
    'T0', T0, ...
    'T0Units', T0Units, ...
    'ceaExe', ceaPath, ...
    'outputDir', ceaResultsDir);

% Prints CEA output into console. 
disp(ceaOut);

fprintf('CEA input and output saved in:\n%s\n', ceaResultsDir);

%% Chamber dimensions and injector port inputs
ChamberDim.outer_chamber_diameter = 52.243;  % [mm]
ChamberDim.inner_chamber_diameter = 37.167;  % [mm]
ChamberDim.chamber_length = 29.21;           % [mm], injector face to exit

InjV.injector_area = 3.7138;  % [mm^2] summed area of one injector port set
InjV.N_holes = 36;           % number of port sets; one for a single opening

%% Injection inputs - original values, not yet used by a model
InjV.P3 = 1e6;        % [Pa] original pre-injection pressure input
InjV.T3 = 398;        % [K] original pre-injection temperature input
InjV.mdot = 0.59;     % [kg/s]
InjV.gamma1 = 1.3569;

InjV.Pa = 1e6;       % [Pa] chamber-side pressure input
InjV.R = 315;        % [J/(kg K)] gas constant
InjV.Cp = 315;       % [J/(kg K)] ORIGINAL placeholder: verify before use
% Cp = R is not a consistent ideal-gas property pair. Kept for review, unused.
% Decide whether P3/T3 are static or stagnation values before deriving injection.

%% Triple-point inputs - defined after CEA has returned its results
% P1 and R1 are chosen bounding-gas inputs; the other inputs come from CEA.
iTripleParam.P1 = P0 * 1e5;  % [Pa] original injection/bounding pressure
iTripleParam.P2 = ceaOut.P_burned_bar * 1e5; % [Pa] CEA pressure: bar -> Pa
iTripleParam.T2 = ceaOut.T_cj;    % [K] post-detonation temperature
iTripleParam.Vcj = ceaOut.cjVel;   % [m/s] detonation propagation velocity
iTripleParam.gamma2 = ceaOut.gamma_burned;
iTripleParam.R1 = ceaOut.R_specific;  % [J/(kg K)] original bounding-gas constant

%% 1. Cell sizing and unwrapped chamber dimensions
% Write and verify the geometry here or in MoC_Calculate_Domain_Size.m.
[xMax, yMax] = MoC_Calculate_Domain_Size(ChamberDim);
MoC_Plot_Field([], [], xMax, yMax); %  no field is solved yet so we're just plotting the chamber domain

%% 2. Injection velocity and the first velocity triangle
% Write and verify the injection model before connecting it to CEA.
% injection = MoC_Injection_Velocity(InjV, ChamberDim);
% TODO: define the laboratory and wave reference frames explicitly.
% TODO: calculate the wave angle and fill height from the chosen triangle.

%% 3. Detonation line and triple point
% TODO: place the angled detonation line and its triple point.
% TODO: define all angle conventions before using beta and delta.
[beta, delta, P_match, M3, M1_p] = MoC_Resolve_Triple_Point(iTripleParam);
fprintf('Beta Shock: %.2f° | Slip Angle: %.2f° | P Matched: %.2f Pa | M3: %.2f | M1'': %.2f\n\n', beta, delta, P_match, M3, M1_p);

%% 4. Initial-value line (IVLine)
% TODO: interpolate seed points along the detonation line.
% TODO: assign and verify the state at each seed point.
% No startup Mach offset or expansion-fan assumption is supplied here.

%% 5. Product-region C+ and C- characteristics
% TODO: implement and verify one internal point before marching a mesh.
% mesh = MoC_Solve_Field(seedPoints, geometry);

%% 6. Slip line, oblique shock, and post-shock characteristics
% TODO: calculate the chosen straight boundary geometry.
% TODO: seed and solve the post-shock region separately.

%% 7. Refill state and region-by-region interpolation
% TODO: assign primitive variables in refill using the verified injection model.
% TODO: document which frame each Mach number uses.
% field = MoC_Interpolate_Field(mesh, geometry);

%% 8. Plot the unwrapped field and characteristic net
% paths = MoC_Trace_Characteristics(field, geometry);
% [fieldFigure, netFigure] = MoC_Plot_Field(field, geometry);

fprintf('\nSolution complete.\n');

