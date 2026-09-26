clear
clc
close all

%% Description:
% In addition to running the entire program this main file initializes by
% computing the 2D injector velocity CJ speed. It constructs the first
% velocity triangle to determine the angle and height of the detonation
% wave. Before interpolating starting points along the detonation wave and
% drawing characteristic lines from them.

%% Definitions:
% IVLine = Initial value line = detonation wave front

% ChamberDim -> MoC_Calculate_Domain.m
% iTripleParam -> MoC_Resolve_Triple_Point.m
% InjV -> velocity_wave_calculation.m


%%% INPUTS %%%

%% Cell sizer from CEA should go here to render the unwrapping hand inputs obsolete.

%% Annulus mean diameter and getting circumference for X and Y max length. (unwrapping RDE)

ChamberDim.outer_chamber_diameter = 50; % [mm]

ChamberDim.inner_chamber_diameter = 35; % [mm]

ChamberDim.chamber_length         = 100; % [mm] Axial length of the chamber measured from the injector face to the exit. 

InjV.injector_area                = 40; % [mm^2] Area of a single set of injector orifices. Used to calculate injector velocity. 
% So total area of a doublet injector would be the sum of the two orifice
% areas. 


%% Triple Point Inputs

iTripleParam.P1     = 200000;  % Injection static pressure (assumed equal to bounding gas initial pressure) [Pa]
iTripleParam.P2     = 150000;  % Post-detonation static pressure [Pa]
iTripleParam.T2     = 800;  % Post-detonation static temperature [K]
iTripleParam.Vcj    = 3486;  % Chapman-Jouguet *detonation wave velocity* [m/s]
iTripleParam.gamma2 = 1.11;  % Specific heat ratio of post-detonation gas
iTripleParam.R1     = 105;  % Specific gas constant of bounding gas [J/kg-K]



%% Drawing the first velocity triangle and getting phi.

%% Delta and beta angles are relative to phi angle. delta + beta = phi angle. 
% calculate retrieve injector velocity and determine phi angle. 

InjV.P3 = 3e6; % [pa] % Initial manifold pressure before entering orifice
InjV.T3 = 398; % [K] % Initial gas temp
InjV.mdot = 0.59; % [kg/s]

% Post injection inputs
InjV.Pa = 101325 ; % Ambient air pressure
InjV.R = 287; % Gas Constant
InjV.Cp = 1005; % Idk what this is

% plug in phi into triple point equation to put all angles into phi


%% Draw Detonation wave line and put triple point at the end of it.

%% Parameters for detonation wave (IVLine) + Fill Height
% The IVLine is derived from a 1D model to calculate injector velocity.
% Use CEA "tp" mode and assign a temperature and pressure to determine
% chamber pressure?


fprintf('\n--- RUNNING WAVE CALCULATIONS ---\n');
out = velocity_wave_calculations(InjV, ChamberDim);

fprintf('\n--- RUNNING TRIPLE POINT RESOLUTION ---\n');
[beta, delta, P_match, M3, M1_p] = MoC_Resolve_Triple_Point(iTripleParam);

% Single line formatted printout for Function 2 outputs
fprintf('Beta Shock: %.2f° | Slip Angle: %.2f° | P Matched: %.2f Pa | M3: %.2f | M1'': %.2f\n\n', beta, delta, P_match, M3, M1_p);


addpath('./RDE_Toolbox_2/src/');


% Logic process -> Draw first velocity triangle -> Draw main lines,
% detonation line, slip line, injector fill line, shock line -> Define
% conditions in each zone -> Draw characteristic lines -> Draw more
% characteristic lines in the slip line and oblique shock. 