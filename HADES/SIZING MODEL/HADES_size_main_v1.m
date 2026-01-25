clear; clc; close;

% put vals in this 
phi = 0.9;

% CEA inputs (wrapper)

initial_temp = 283.15; %K (10C)
temp_units = "K";
initial_pressure = 26; %psia
pressure_units = "psia";

% print everythign
fprintf("Given: \n");
fprintf("Equivalence Ratio = %.1f\n", phi);
fprintf("Initial Pressure: %.3f %s \n\n", initial_pressure, pressure_units);


%calcs
fprintf("Outputs: \n");
[cellWidth, annulus_gap, fillHeight] = HADES_size_annulusgap_fillheight(phi);
