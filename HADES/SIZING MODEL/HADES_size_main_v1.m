clear; clc; close all;

% put vals in this 
phi = 0.9;
outer_radius = 2; %in
wall_thickness = 0.625; % in
wave_modes = 1; %our fav assumption

ambient_pressure_Mpa = 0.101325; %ambient pressure
ox_mdot = 2; %lbm/s
fuel_mdot = ox_mdot * phi / 34.3; %also lbm/s
total_mdot = ox_mdot + fuel_mdot;


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
ceaRock_results = HADES_size_ceaRocket('ox', ox_type,'fuel',fuel_type,'phi', phi,'Pc',avg_chamber_p,'PcUnits', pressure_units);
[thrust] = HADES_size_thrust(ceaRock_results.isp, total_mdot);
Failure_temps = HADES_size_HoopStressTemps(outer_radius - wall_thickness, outer_radius, ambient_pressure_Mpa, HADES_size_convertPressure(initial_pressure * ceaDet_results.P_ratio, pressure_units, 'bar'));
[P_0] = HADES_size_P0_v1();

[h2_area] = HADES_size_h2Inj(fuel_mdot, initial_temp, ambient_pressure_Mpa);
