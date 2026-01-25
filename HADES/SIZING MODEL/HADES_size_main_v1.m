clear; clc; close all;

% put vals in this 
phi = 0.9;
outer_radius = 2; %in
wall_thickness = 0.625; % in
wave_modes = 2; %our fav assumption

%% CEA inputs (wrapper)
initial_temp = 283.15; %K (10C)
temp_units = 'K';
initial_pressure = 26; %psia
pressure_units = 'psia';
fuel_type = 'H2';
ox_type = 'Air';

%% print givens
fprintf("Given: \n");
fprintf("Equivalence Ratio = %.1f\n", phi);
fprintf("Initial Pressure: %.3f %s \n\n", initial_pressure, pressure_units);


%% calcs
fprintf("Outputs: \n");
[cellWidth, annulus_gap, fillHeight] = HADES_size_annulusgap_fillheight(phi);
[det_wave_path_length] = HADES_size_geometry(outer_radius, wall_thickness, annulus_gap);
ceaDet_results = HADES_size_ceaDet('ox',ox_type,'fuel',fuel_type,'phi', phi,'P0', initial_pressure,'P0Units',pressure_units,'T0', initial_temp,'T0Units', temp_units);
[avg_chamber_p] = HADES_size_chamberPressure(wave_modes, det_wave_path_length, ceaDet_results.cjVel, ceaDet_results.P_ratio, initial_pressure);


