%% MoC_Main_Model - starting framework for the hand-written RDE model
% This is a script: inputs and ceaOut remain visible in the MATLAB workspace.
% CEA results are passed to the triple-point resolver below.


%% TODO: FIGURE OUT WHAT IS GOING ON WITH THE SLIP
%% SHOCK REGION HAVING CHAR LIES BE SO STRAIGHT AND REFLECTING AT 90 DEGREES.
% POSSIBLE IMMEDIATE IDEAS FOR CAUSES
% 1. Coordinate system adjustment making it look like the char lines
%  reflecting at 90 degrees. 
% 2. PM angle calc could be wrong. 
% 3. The lines coming from the points themselves are very straight as well
% so it's prob some bs to do with the existing calcs for angles.  
% 90 degree angle lines means that region is likely slow. 
% WRite the oblique shock jump conditions part of the paper. 
% Work on delta_h calc. 
clear
clc
close all

tic;

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

%% Chamber dimensions and injector port inputs
ChamberDim.outer_chamber_diameter = 52.243;  % [mm]
ChamberDim.inner_chamber_diameter = 37.167;  % [mm]
ChamberDim.chamber_length = 29.21;           % [mm], injector face to exit

InjV.injector_area = 3.7138;  % [mm^2] summed area of one injector port set
InjV.N_holes = 36;           % number of port sets; one for a single opening

%% Plot Selection
% Select which output graphs to display by providing a vector of numbers:
%   1 = Method of Characteristics Net
%   2 = Static Temperature Field [K]
%   3 = Static Pressure Field [bar]
%   4 = Mach Number Field
%   5 = Flow Angle (Theta) Field [deg]
%   6 = Density Field [kg/m^3]
%   7 = Velocity Magnitude Field [m/s]
%   8 = Characteristic Lines (products: C+ red, C- blue; slip-shock: green)
% Examples:
%   plots_requested = [1, 2];       % MoC Net and Temperature
%   plots_requested = [3];          % Only Pressure
%   plots_requested = [2, 3, 5];    % Temperature, Pressure, and Flow Angle
%   plots_requested = [8];          % Characteristic lines instead of points

plots_requested = [2, 4, 5, 3, 8];

%% Plot Options
% Select which plot options to use by providing a vector of numbers:
%   1 = Show Slip Line 
%   2 = Show Shock Line
%   3 = Show Infill Region
%   4 = Show Detonation Triangle

plot_options = [1, 2];

%% Injection inputs 

InjV.Pa = 1e6;       % [Pa] chamber-side pressure input

%% Solver inputs

% Set the max and min number of points along the various lines, if you want to use the non-optimized version, toggle is below.
num_det_points = 250; % Max number of points along the detonation line (from triple point to floor)
min_det_points = 50;  % Minimum number of points along the detonation line (from triple point to floor)
num_fan_points = 300; % Max number of points along the Prandtl-Meyer expansion fan (from triple point to floor)
min_fan_points = 50;  % Minimum number of points along the Prandtl-Meyer expansion fan (from triple point to floor)
slipShockOptions.maxSteps = 100000;

% delta_h optimization toggle: If off, the solver will use the max number of points for every iteration.
delta_h_opt_toggle = 1; % 1 = on, 0 = off

if delta_h_opt_toggle == 0
    min_det_points = num_det_points;
    min_fan_points = num_fan_points;
    fprintf('Delta_h optimization is OFF. Using max points for every iteration.\n');

elseif delta_h_opt_toggle == 1
    fprintf('Delta_h optimization is ON. Using min points for first iteration, then scaling to max points as convergence is approached.\n');
end 


tol = 1e-4; % Tolerance for the y_inj iteration convergence
max_iter = 50; % Maximum number of iterations for the y_inj convergence loop
y_inj_update_method = 2; % Method for updating y_inj: 1 = simple update, 2 = weighted average



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
%disp(ceaOut);

%fprintf('CEA input and output saved in:\n%s\n', ceaResultsDir);

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

%% 2. Detonation line and triple point
% Calculate the triple point properties. This depends only on CEA/gas parameters.
[beta, delta, P_match, M3, M1_p] = MoC_Resolve_Triple_Point(iTripleParam);
fprintf('Beta Shock: %.2f degrees | Slip Angle: %.2f degrees | P Matched: %.2f Pa | M3: %.2f | M1'': %.2f\n\n', beta, delta, P_match, M3, M1_p);

%% 3. Target Pressure for Infill Triangle
% Find the target Prandtl-Meyer angle (nu_target) where the chamber pressure equals injection pressure (InjV.Pa)
% Because we are assuming isentropic expansion, we calculate Nu on every line that hits the floor, and when it equals
% the target Nu, we can determine where to start the infill triangle.

gamma = iTripleParam.gamma2;
% gp1 and gm1 are just used to make the code cleaner
gp1 = gamma + 1; 
gm1 = gamma - 1;

% Stagnation pressure of the characteristic field (constant for this MOC)
P0_mesh = P_match * (1 + (gm1 / 2) * M3.^2)^(gamma / gm1);
P_target = InjV.Pa;

% Target Mach number
M_target = sqrt((2 / gm1) * ((P0_mesh / P_target)^(gm1 / gamma) - 1));

% Target Nu
nu_target = (sqrt(gp1/gm1)*atan(sqrt((gm1/gp1)*(M_target^2 - 1))) - atan(sqrt(M_target^2 - 1))) * 180/pi;

% Initial Nu based on post-detonation Mach M3
M2_m1 = (1.001)^2 - 1; % Used to clean up the code a bit, just means M^2 - 1
if M2_m1 < 0 % Checks for mach numbers below 0 and sets them to 0 to avoid complex numbers in the atan function
              % This could happen due to rounding errors in the MOC solver, but it should not happen in a real flow.
    M2_m1 = 0; 
end

% Calculates the nu value for the post-detonation Mach number M3
nu_initial = (sqrt(gp1/gm1)*atan(sqrt((gm1/gp1)*M2_m1)) - atan(sqrt(M2_m1))) * 180/pi; 

%% 4. Iterative Loop to Find y_inj (Convergence Condition)
v_inj = 150; % TODO: placeholder injection velocity [m/s]
y_inj = .3 * ChamberDim.chamber_length; % Initial guess for injection triangle height [mm] (Converged value is independent of this)

delta_h = Inf; % Initialize delta_h for mesh scaling on iteration 1

fprintf('--- Starting y_inj iteration ---\n');
for iter = 1:max_iter
    % Determine mesh resolution as a function of delta_h
    if iter == 1 || isinf(delta_h)
        cur_det_points = min_det_points;
        cur_fan_points = min_fan_points;
    else
        % Scale points based on |delta_h|: as error drops toward tol, points increase to max
        delta_h_scale = 0.05; % [mm] threshold below which mesh progressively refines
        refine_factor = min(1.0, max(0.0, 1.0 - abs(delta_h) / delta_h_scale));
        
        cur_det_points = round(min_det_points + refine_factor * (num_det_points - min_det_points));
        cur_fan_points = round(min_fan_points + refine_factor * (num_fan_points - min_fan_points));
    end

    % Solves for y_inj assuming that infill does not exist, then it will 
    % be added after, based on the converged infill height.

    [phi, ~, x_inj] = MoC_Infill_Triangle(v_inj, iTripleParam.Vcj, y_inj);
    
    % 1. Seed points along the actual angled Detonation Line
    x_det_floor = y_inj * tand(phi); % Calculate the x-coordinate where the detonation line hits the floor
    X_det = linspace(x_det_floor, 0, cur_det_points); % Calculate the x-coordinates along the detonation line
    Y_det = linspace(0, y_inj - 1e-4, cur_det_points); % Calculate the y-coordinates along the detonation line
    
    Theta_det = zeros(1, cur_det_points); % Pre-allocate vector for detonation line angles
    Nu_det = repmat(nu_initial, 1, cur_det_points); % Pre-allocate vector for detonation line Prandtl-Meyer angles
    
    % 2. Prandtl-Meyer Expansion Fan at the Triple Point (0, y_inj)
    Theta_fan = linspace(0, delta + phi, cur_fan_points); % Pre-allocate vector for fan angles
    Theta_fan = Theta_fan(2:end); % Remove the first point to avoid duplication at the triple point
    Nu_fan = nu_initial + Theta_fan; % Calculate the Prandtl-Meyer angles for the fan points
    
    X_fan = zeros(1, length(Theta_fan)); % x-coordinates for the fan points (all at x=0)
    Y_fan = linspace(y_inj - 0.8e-4, y_inj, length(Theta_fan)); % y-coordinates for the fan points
    
    seedPoints.X = [X_det, X_fan]; % Combine detonation line and fan x-coordinates
    seedPoints.Y = [Y_det, Y_fan]; % Combine detonation line and fan y-coordinates
    seedPoints.M = zeros(1, cur_det_points + length(Theta_fan)); % Pre-allocate Mach numbers
    seedPoints.M(1:cur_det_points) = 1.001; % Set Mach numbers for detonation line points
    seedPoints.Theta = [Theta_det, Theta_fan]; % Combine detonation line and fan angle flow angle (theta) vectors
    seedPoints.Nu = [Nu_det, Nu_fan]; % Combine detonation line and fan Prandtl-Meyer angle (Nu) vectors
    
    for i = 1:length(Theta_fan)
        seedPoints.M(cur_det_points + i) = MoC_PM_Inv(Nu_fan(i), gamma); % Calculate Mach number for the fan point based on its PM angle
    end 
    
    % 3. Set up geometry struct to test expansion along a solid floor
    geometry.xMax = xMax;
    geometry.y_engine_top = 160; % Calculate higher for a fuller mesh (prevents mesh from being too sparse near the top and cutting off)
    geometry.x_infill_start = xMax * 10; % Push infill way out so floor is entirely solid 
                                         % for the first run to get a good initial mesh without considering infill
    geometry.infill_nu = nu_target; 
    geometry.infill_phi = phi; 
    geometry.slip_angle = delta + phi;
    geometry.shock_angle = beta + phi; 
    geometry.x_slip0 = 0;
    geometry.y_slip0 = y_inj;
    geometry.y_inj = y_inj;
    geometry.yMax = ChamberDim.chamber_length; % Set the maximum y-coordinate for the geometry to the chamber length
    geometry.beta = beta; % Set the shock angle for the geometry to the calculated beta
    geometry.delta = delta; % Set the slip angle for the geometry to the calculated delta
    
    gas.gamma = gamma;
    
    % Run solver to find floor pressure crossing
    field = MoC_Solve_Field(seedPoints, geometry, gas);
    
    % Find X on the floor where Nu crosses nu_target
    idx_floor = find(field.Type == 1); % Filter for all floor points (Type 1)
    x_floor = field.X(idx_floor); % Extract all x values of the floor points
    nu_floor = field.Nu(idx_floor); % Extract all Nu values of the floor points
    
    [x_floor, sort_idx] = sort(x_floor); % Sort the x values of the floor points in ascending order
    nu_floor = nu_floor(sort_idx); % Sort the Nu values of the floor points based on the sorted x values
    
    if max(nu_floor) < nu_target % Check if the maximum Nu on the floor is less than the target Nu
        warning('Pressure never dropped to Inj.Pa. Infill cannot start.');
        x_infill = xMax;
    else
        [unique_nu, unique_idx] = unique(nu_floor); % Find unique Nu values and their corresponding indices
        x_infill = interp1(unique_nu, x_floor(unique_idx), nu_target); % Interpolate to find the x-coordinate where Nu equals nu_target
    end
    
    y_inj_new = (xMax - x_infill) * tand(phi); % Calculate the new y_inj based on the x-coordinate of the infill start and the angle phi
    delta_h = y_inj_new - y_inj; % Calculate the difference between the new and old y_inj values
    
    fprintf('Iter %d [det: %d, fan: %d]: y_inj = %.3f, x_infill = %.2f, y_inj_new = %.5f, delta_h = %.3f\n', ...
        iter, cur_det_points, cur_fan_points, y_inj, x_infill, y_inj_new, delta_h);
    
    if abs(y_inj_new - y_inj) < tol % Check for convergence based on tolerance
        y_inj = y_inj_new;
        x_inj = xMax - x_infill;
        fprintf('Converged!\n\n');
        break;
    end
    
    % Choose whether to use a simple update or a weighted average for y_inj.
    if y_inj_update_method == 1
        y_inj = y_inj_new; % Simple update
    elseif y_inj_update_method == 2
        y_inj = 0.5 * y_inj + 0.5 * y_inj_new; % Weighted average update
    end
   
end

%% 5. Final Geometry and Solved Field
% Regenerate seed points at final converged y_inj using MAX number of points
x_det_floor = y_inj * tand(phi);
X_det = linspace(x_det_floor, 0, num_det_points);
Y_det = linspace(0, y_inj - 1e-4, num_det_points);
Theta_det = zeros(1, num_det_points);
Nu_det = repmat(nu_initial, 1, num_det_points);

Theta_fan = linspace(0, delta + phi, num_fan_points);
Theta_fan = Theta_fan(2:end);
Nu_fan = nu_initial + Theta_fan;
X_fan = zeros(1, length(Theta_fan));
Y_fan = linspace(y_inj - 0.8e-4, y_inj, length(Theta_fan));

seedPoints.X = [X_det, X_fan];
seedPoints.Y = [Y_det, Y_fan];
seedPoints.M = zeros(1, num_det_points + length(Theta_fan));
seedPoints.M(1:num_det_points) = 1.001;
seedPoints.Theta = [Theta_det, Theta_fan];
seedPoints.Nu = [Nu_det, Nu_fan];

for i = 1:length(Theta_fan)
    seedPoints.M(num_det_points + i) = MoC_PM_Inv(Nu_fan(i), gamma);
end

% Re-solve the final characteristic field with the infill boundary at max points
geometry.x_infill_start = xMax - x_inj;
geometry.y_inj = y_inj;
geometry.shock_angle = beta + phi;
fprintf('Solving final characteristic field with infill boundary at max resolution (%d det, %d fan points)...\n', num_det_points, num_fan_points);
[field, shockSeeds] = MoC_Solve_Field(seedPoints, geometry, gas);
fprintf('Solved %d points.\n\n', length(field.X));

%% 5.1 Compute Thermodynamic and Kinematic Field Properties
T0_mesh = ceaOut.T_cj * (1 + 0.5 * (gamma - 1)); % Stagnation Temperature [K]
R_gas   = ceaOut.R_specific;                     % Specific gas constant [J/(kg*K)]

% 1. Temperature [K]
field.T = T0_mesh ./ (1 + 0.5 * (gamma - 1) .* field.M.^2);

% 2. Static Pressure [Pa]
field.P = P0_mesh ./ ((1 + 0.5 * (gamma - 1) .* field.M.^2).^(gamma / (gamma - 1)));

% 3. Density [kg/m^3]
field.rho = field.P ./ (R_gas .* field.T);

% 4. Velocity Components [m/s]
field.a   = sqrt(gamma * R_gas .* field.T);     % Local speed of sound [m/s]
field.V   = field.M .* field.a;                 % Total velocity magnitude [m/s]
field.v_x = field.V .* cosd(field.Theta);       % Horizontal / circumferential velocity [m/s]
field.v_y = field.V .* sind(field.Theta);       % Vertical / axial velocity [m/s]

%% 5.2 Irrotational MoC in Shock Slip Region
% Shock seeds already contain the jump in Mach, angle, and total-pressure ratio.
% Solve Mach/angle with MoC; map P0/T0 along straight lines parallel to the slip.
gas.R = R_gas;
gas.P0 = P0_mesh;
gas.T0 = T0_mesh;
shockSlipField = MoC_Solve_Slip_Shock(shockSeeds, geometry, field, gas, slipShockOptions);
fprintf('Shock-slip region: %s (%d points).\n', ...
    shockSlipField.Diagnostics.status, numel(shockSlipField.X));

% Append solved region points to the same field used by the existing plot.
% Shock seeds already exist in field: update their thermodynamics, do not
% duplicate them. Convert the new solver's local parent indices to field IDs.
old_count = numel(field.X);
new_points = find(shockSlipField.Type ~= 5);
seed_points = find(shockSlipField.Type == 5);
existing_shock = find(field.Type == 5);
point_map = nan(numel(shockSlipField.X),1);
point_map(new_points) = old_count + (1:numel(new_points));
point_map(seed_points) = existing_shock(shockSlipField.ShockSeed(seed_points));

field.P0 = repmat(P0_mesh,old_count,1);
field.T0 = repmat(T0_mesh,old_count,1);
field.ShockSeed = nan(old_count,1);
point_names = {'X','Y','M','Theta','Nu','P','T','rho','P0','T0', ...
    'a','V','v_x','v_y','Type','Region','ShockSeed'};
for k = 1:numel(point_names)
    name = point_names{k};
    field.(name)(point_map(seed_points)) = shockSlipField.(name)(seed_points);
    field.(name) = [field.(name); shockSlipField.(name)(new_points)];
end
for name = {'Parent_C_Minus','Parent_C_Plus'}
    parents = shockSlipField.(name{1})(new_points);
    valid_parent = ~isnan(parents);
    parents(valid_parent) = point_map(parents(valid_parent));
    field.(name{1}) = [field.(name{1}); parents];
end

field.SlipShockDiagnostics = shockSlipField.Diagnostics;

%% 5.1.1 Infill Region Point Creation

num_infill_points = round(0.11 * (length(field.X))); % Number of points to create inside the infill region
point_spacing = sqrt((xMax - geometry.x_infill_start) * y_inj / (2 * num_infill_points)); % Spacing between points in the infill region [mm]

% Define horizontal velocity component v_cj 
v_cj = iTripleParam.Vcj;

% Infill static pressure in Pascals [Pa] 
P_infill = P0 * 1e5;

% Generate evenly spaced 2D grid of dots inside the infill triangle
dx = point_spacing;
dy = point_spacing;
x_infill_vals = (geometry.x_infill_start + dx/2) : dx : xMax;
y_infill_vals = (dy/2) : dy : y_inj;
[X_grid, Y_grid] = meshgrid(x_infill_vals, y_infill_vals);

% Hypotenuse boundary from (geometry.x_infill_start, 0) to (xMax, y_inj)
slope_infill = y_inj / (xMax - geometry.x_infill_start);
inside_infill = (Y_grid <= slope_infill * (X_grid - geometry.x_infill_start)) & ...
                (X_grid <= xMax) & (Y_grid >= 0);

infill_X = X_grid(inside_infill);
infill_Y = Y_grid(inside_infill);
N_infill = numel(infill_X);

% Append infill points to field struct with specified properties
field.X      = [field.X; infill_X];
field.Y      = [field.Y; infill_Y];
field.Region = [field.Region; ones(N_infill, 1)];
field.T      = [field.T; repmat(T0, N_infill, 1)];
field.P      = [field.P; repmat(P_infill, N_infill, 1)];
field.Theta  = [field.Theta; repmat(phi, N_infill, 1)];
field.v_x    = [field.v_x; repmat(v_cj, N_infill, 1)];
field.v_y    = [field.v_y; repmat(v_inj, N_infill, 1)];

% Mark any unspecified point-data properties as NaN
unspecified_fields = {'M', 'Nu', 'rho', 'a', 'V', 'Type', 'Parent_C_Minus', 'Parent_C_Plus'};
for k_fn = 1:numel(unspecified_fields)
    fn_name = unspecified_fields{k_fn};
    if isfield(field, fn_name)
        field.(fn_name) = [field.(fn_name); nan(N_infill, 1)];
    end
end

%% 6. Plotting - Method of Characteristics Net and RDE Flow Field
geometry.y_inj = y_inj;
geometry.x_inj = x_inj;
geometry.phi = phi;
geometry.beta = beta;
geometry.delta = delta;
geometry.shock_angle = beta + phi;
geometry.slip_angle = delta + phi;
geometry.xMax = xMax;
geometry.yMax = ChamberDim.chamber_length;
geometry.chamber_length = ChamberDim.chamber_length;
geometry.x_infill_start = xMax - x_inj;

% Moves the x-coordinates of the field points to be within the unwrapped chamber circumference
% Retain marching coordinates so characteristic segments wrap at the seam.
field.X_unwrapped = field.X;

for i = 1:length(field.X)
    if field.X(i) >= geometry.xMax 
        field.X(i) = field.X(i) - geometry.xMax;
    end

    if field.Y(i) > geometry.yMax
        field.Region(i) = 4; % Mark points above the chamber length as Region 4 (outside the chamber)
    end
end

field.plot_logic = (field.Region ~= 4) & ~isnan(field.Region); % Logical index for plotting only points within the chamber (Region 1, 2, and 3)
field.datapoint_indices = (1:length(field.X))'; % Index numbers for all points in the field

MoC_Plot_Field(field, geometry, plots_requested, plot_options);

toc;
