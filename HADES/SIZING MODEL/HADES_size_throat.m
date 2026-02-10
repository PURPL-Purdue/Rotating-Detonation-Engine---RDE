function [At] = HADES_size_throat(mdot, P0, T0, gamma, R)
% Calculates the minimum throat area for choked RDE nozzle.
%
% Inputs:
% mdot - Mass flow rate (kg/s)
% P0 - Stagnation pressure (Pa)
% T0 - Stagnation temperature (K)
% gamma - Ratio of specific heats
% R - Specific gas constant (J/kg-K)
%
% Output:
% At - Minimum throat area (m^2)

    term1 = gamma / R;
    term2 = (2 / (gamma + 1))^((gamma + 1) / (gamma - 1));
    
    flow_factor = sqrt(term1 * term2);

    % Calculate Throat Area
    At = (mdot * sqrt(T0)) / (P0 * flow_factor);
end