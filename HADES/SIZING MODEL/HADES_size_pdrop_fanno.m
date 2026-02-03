function [P_inlet, P_exit, Ma_inlet, f] = HADES_size_pdrop_fanno(P0, D, L, epsilon, m_dot, gamma, mu)
% HADES_size_pdrop_fanno: Choked Fanno flow with Tkachenko friction factor
% 
% Inputs:
% P0 - Stagnation Pressure (Pa)
% D - Diameter (m)
% L - Length (m)
% epsilon - Surface roughness (m)
% m_dot - Mass flow rate (kg/s)
% gamma - Ratio of specific heats
% mu - viscosity (Pa-s)

    %% Reynolds Number and Friction Factor
    Re = (4 * m_dot) / (pi * D * mu);
    rel_rough = epsilon / D;
    
    A0 = -0.79638 * log( (rel_rough / 8.208) + (7.3357 / Re) );
    A1 = Re * rel_rough + 9.3120665 * A0;
    
    num = 8.128943 + A1;
    den = 8.128943 * A0 - 0.86859209 * A1 * log( A1 / (3.7099535 * Re) );
    f = (num / den)^2;

    %% Fanno Flow Parameters
    fanno_param = (4 * f * L) / D; 
    
    options = optimset('Display','off');
    f_to_solve = @(Ma) ((1 - Ma.^2)./(gamma * Ma.^2)) + ((gamma + 1)/(2 * gamma)) * log(((gamma + 1) * Ma.^2) ./ (2 + (gamma - 1) * Ma.^2)) - fanno_param;

    Ma_inlet = fzero(f_to_solve, [0.001, 0.999], options);
    
    %% Pressure Calculations
    P_inlet = P0 * (1 + (gamma-1)/2 * Ma_inlet^2)^(-gamma/(gamma-1));
    P_ratio_star = (1/Ma_inlet) * sqrt((gamma+1) / (2 + (gamma-1)*Ma_inlet^2));
    P_exit = P_inlet / P_ratio_star;
end