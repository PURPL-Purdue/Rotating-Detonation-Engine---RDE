clear; clc; close all;

% put vals in this 
phi = 0.9;
outer_radius = 2; %in
wall_thickness = 0.625; % in
wave_modes = 1; %our fav assumption

azi_annulus_length = 0.2893; %Azimuthal annulus length [m]
det_wave_num = 1; %Detonation wave number [-]
CJ_det_speed = 1931.5; %Chapman-Jouguet detonation speed [m/s]
burned_gas_mach = 1075.6; %Sound speed in burned gas [m/s]
burned_gas_p = 32.33; %Pressure of burned gas region [Pa or bar, consistent with Pc]
inj_crit_p = 2; %Injector critical pressure [same units as P2]
unburned_axial_vel = 391; %Axial velocity of unburned propellant [m/s]
surface_roughness_304 = 3.2e-6; %Stainless steel 304 surface roughness [m]
injector_dia =1e-3; % [m]
injector_length = 0.25; %[in]



ambient_pressure_Mpa = 0.101325; %ambient pressure
ox_mdot = 2; %lbm/s
fuel_mdot = ox_mdot * phi / 34.3; %lbm/s
total_mdot = ox_mdot + fuel_mdot;




%% CEA inputs (wrapper)
initial_temp = 283.15; %K (10C)
temp_units = 'K';
initial_pressure = 26; %psia
pressure_units = 'psia';
fuel_type = 'H2';
ox_type = 'Air';

%% print givens
fprintf("Inputs: \n");
fprintf("Equivalence Ratio = %.1f\n", phi);
fprintf("Initial Pressure: %.3f %s \n", initial_pressure, pressure_units);
fprintf("Initial Temp: %.2f %s \n\n",initial_temp, temp_units);

%% calcs
fprintf("Outputs: \n");
[cellWidth, annulus_gap, fillHeight] = HADES_size_annulusgap_fillheight(phi, azi_annulus_length, det_wave_num, CJ_det_speed, burned_gas_mach, burned_gas_p, inj_crit_p, unburned_axial_vel);

[det_wave_path_length] = HADES_size_geometry(outer_radius, wall_thickness, annulus_gap);

ceaDet_results = HADES_size_ceaDet('ox',ox_type,'fuel',fuel_type,'phi', phi,'P0', initial_pressure,'P0Units',pressure_units,'T0', initial_temp,'T0Units', temp_units);

[P_inlet, P_exit, Ma_inlet, f] = HADES_size_pdrop_fanno(20e-5, injector_dia, injector_length, surface_roughness_304, HADES_size_convert(fuel_mdot,'lbms', 'kgs') , 1.41, ceaDet_results.Mu * 1e-4);

[avg_chamber_p] = HADES_size_chamberPressure(wave_modes, det_wave_path_length, ceaDet_results.cjVel, ceaDet_results.P_ratio, HADES_size_convert(P_exit * 1e5, 'bar', 'psia'));

ceaRock_results = HADES_size_ceaRocket('ox', ox_type,'fuel',fuel_type,'phi', phi,'Pc',avg_chamber_p,'PcUnits', pressure_units);

[thrust] = HADES_size_thrust(ceaRock_results.isp, total_mdot);
Failure_temps = HADES_size_HoopStressTemps(outer_radius - wall_thickness, outer_radius, ambient_pressure_Mpa, HADES_size_convert(initial_pressure * ceaDet_results.P_ratio, pressure_units, 'bar'));

[P_0] = HADES_size_P0_v1();

[h2_area] = HADES_size_h2Inj(fuel_mdot, initial_temp, ambient_pressure_Mpa);

[air_area] = HADES_size_airInj(h2_area, annulus_gap, outer_radius, initial_temp, initial_pressure);
