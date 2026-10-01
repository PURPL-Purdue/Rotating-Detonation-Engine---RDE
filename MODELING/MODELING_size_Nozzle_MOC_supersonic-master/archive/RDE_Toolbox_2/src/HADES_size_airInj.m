%%air sizing

%%% needed function inputs
% hydrogen area estimation

function [air_area, Area_Annulus] = HADES_size_airInj(ox_mdot, an_gap, an_rad_outter, temp)

% givens
R = 8.314; % R value for Ideal Gas Law
Cd = 0.6; %Coefficient of Discharge
Delta_P = 1000000; % Pa
ambient_pressure_air = 1000000; %Pa
molarmass_air = 28.97; %g/mol


% Total Area Air Calculations
Air_mass_FR = ox_mdot * 0.453592; % convert lbm/s to kg/s
density = (ambient_pressure_air  * molarmass_air)/ (R * temp * 1000); %kg/m^3
air_area = (Air_mass_FR / (Cd * sqrt(2 * density * Delta_P))) * 1e6; %calculate air area in mm^2
fprintf("\nTotal Area Air: %.5f mm^2\n", air_area);

% Annulus Area Calculations
an_rad_outter_mm = an_rad_outter * 25.4; %convert to mm
Area_Annulus = pi * (an_rad_outter_mm^2) - pi * (an_rad_outter_mm - an_gap)^2; % mm^2
fprintf("\nArea Annulus: %.5f mm^2\n", Area_Annulus);

% Air injector width
air_inj_width = an_rad_outter_mm - sqrt(an_rad_outter_mm^2 - (air_area / pi)); % mm
fprintf("\nAir Injector Width: %.5f mm\n", air_inj_width);



