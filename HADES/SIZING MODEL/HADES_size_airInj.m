%%air sizing

%%% needed function inputs
% hydrogen area estimation

function [air_area] = HADES_size_airInj(h2_area, an_length, an_rad_1, temp, psia)

% givens
pressure = psia * 6894.76; % converts psia to Pa
R_m = 8.314 / 28.97; %% ****TBH IDK WHERE THIS CAME FROM, I JUST PUT IT IN HERE
Cd = 0.6; %% ****TBH IDK WHERE THIS CAME FROM, I JUST PUT IT IN HERE
Delta_P = 1000000; % Pa

% given in sheets but not used in function
% mach = .8; % assumed
% air_temp = temp;
% gamma = 1.4;
% R_air = .287; % kj/kg*K

% Totoal Area Air Calculations
Air_mass_FR = h2_area * 0.453592; % kg/s
density = pressure / (R_m * temp) / 1000; 
air_area = Air_mass_FR / (Cd * sqrt(2 * density * Delta_P)) * (10 ^ 6);
fprintf("\nTotal Area Air: %.5f mm^2\n", air_area);

% Annulus Area Calculations
an_rad_2 = an_rad_1 * 25.4; % mm
Area_Annulus = pi * (an_rad_2^2) - pi * (an_rad_2 - an_length)^2; % mm^2
fprintf("\nArea Annulus: %.5f mm^2\n", Area_Annulus);

% Air injector width
air_inj_width = an_rad_2 - sqrt(an_rad_2^2 - (air_area / pi)); % mm
fprintf("\nAir Injector Width: %.5f mm\n", air_inj_width);



