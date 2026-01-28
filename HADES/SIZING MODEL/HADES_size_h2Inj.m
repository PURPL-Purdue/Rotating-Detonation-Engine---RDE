%%H2 inj sizing

% make function define inputs ONLY FROM MAIN
function [h2_area] = HADES_size_h2Inj(fuel_mdot, initial_temp, ...
ambient_pressure)

%defining constants
gamma = 1.4; % specific heat ratio for hydrogen
R = 4124; % specific gas constant for hydrogen in J/(kg*K)

% make print statemnts
fprintf("Given: \n")
fprintf("gamma (specific heat ratio) = 1.4 \n")
fprintf("R (specific gas constant for hydrogen in J/(kg*K) = 4124")

%choked flow equation
h2_area = (fuel_mdot * sqrt(initial_temp) * sqrt(R/gamma))/...
    (ambient_pressure * ((gamma + 1)/2) .^ (-1 / 2 * (gamma + 1)/...
    (gamma - 1)));
