%%H2 inj sizing

% make function define inputs ONLY FROM MAIN
%<<<<<<< Updated upstream
%<<<<<<< Updated upstream
%function [h2_area] = HADES_size_h2Inj(fuel_mdot, initial_temp, ambient_pressure)
%=======
%function [h2_area] = HADES_size_h2Inj(fuel_mdot, chamber_pressure)
%>>>>>>> Stashed changes
%=======
function [h2_area] = HADES_size_h2Inj(fuel_mdot, chamber_pressure)
%>>>>>>> Stashed changes

%defining constants
stiffness = 100;
fuel_mdot = fuel_mdot * .453592; %%fuel mass flow rate converted to kg/s
dischargeCoef = .7; %Discharge coefficient
density = 1.012; %%UNSURE OF WHERE THIS NUM IS FROM
manifold_pressure = chamber_pressure * (1 + (stiffness / 100)); %%manifold pressure (psia)
pressureDrop = (manifold_pressure - chamber_pressure) * 6895; %%pressure drop converted to PA


% make print statemnts
fprintf("Given: \n")
fprintf("gamma (specific heat ratio) = 1.4 \n")
fprintf("R (specific gas constant for hydrogen in J/(kg*K) = 4124\n")



h2_area = (fuel_mdot / dischargeCoef) * sqrt(1 / (2 * density * pressureDrop));
h2_area = h2_area * 1000000;

fprintf("\nHydrogen Area = " + h2_area + " mm^2 \n");

orificeSize = 1; %%diameter mm
orificeArea = pi * (orificeSize / 2) ^ 2;
orificeCount = h2_area / orificeArea;

fprintf("Orifice count: " + orificeCount);
end