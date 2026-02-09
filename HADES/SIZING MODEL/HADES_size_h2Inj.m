%%H2 inj sizing

% make function define inputs ONLY FROM MAIN
function [h2_area] = HADES_size_h2Inj(fuel_mdot, initial_temp, h2_inj_density)

% to check inputs
fprintf("\n-------------\nh2 mdot:%.4f\ninitial temp: %.3f\nambient pressure: %.3f\n", fuel_mdot,initial_temp, h2_inj_density)

%defining constants
%gamma = 1.4; % specific heat ratio for hydrogen
%R = 4124; % specific gas constant for hydrogen in J/(kg*K)

% make print statemnts
fprintf("Given: \n")
fprintf("gamma (specific heat ratio) = 1.4 \n")
fprintf("R (specific gas constant for hydrogen in J/(kg*K) = 4124\n")

%% this was another test output that works
mani_p = 1.9995 * 10^6; % 290 converted to Pa
cham_p = 999739.807; % 145 converted to Pa
C_d = 0.71;
h2_area = (fuel_mdot / 2.205) / (C_d * sqrt(2 * h2_inj_density *(mani_p - cham_p)));
h2_area = h2_area * 1e6; % Converting m^2 to mm^2
fprintf("\nh2 area: %.4f mm^2 \n", h2_area);

