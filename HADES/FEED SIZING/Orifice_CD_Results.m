%% WASTE PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%



clear;
clc;
close all;

%% CONVERSION FACTORS

psi_to_Pa = 6894.75729;

%% INPUTS

gas = ""

%% READ CSV FILE

filename = 'orifice_blowdown_input.csv';

data = readtable(filename, ...
    'VariableNamingRule', 'preserve');


%% EXTRACTING DATA

time_ms = data.("time_ms");

supply_pressure = data.("supply_pressure_psi");
upstream_pressure = data.("upstream_pressure_psi");
downstream_pressure = data.("downstream_pressure_psi");

supply_temperature = data.("supply_temperature_K");
upstream_temperature = data.("upstream_temperature_K");
downstream_temperature = data.("downstream_temperature_K");

%% CALCULATIONS

% Starting Supply Conditions
start_pressure = supply_pressure(1) .* psi_to_Pa;
start_temp = supply_temperature(1);

end_pressure = supply_pressure(end);
