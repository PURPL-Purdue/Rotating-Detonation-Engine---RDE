function [M2, theta2, nu2, shockProps] = MoC_Oblique_Shock_Jump(M1, theta1, beta_shock, gamma)
% MOC_OBLIQUE_SHOCK_JUMP Calculates post-shock flow properties across an oblique shock.
%
% Inputs:
%   M1         - Upstream Mach number
%   theta1     - Upstream flow angle [degrees] (relative to horizontal)
%   beta_shock - Oblique shock wave angle [degrees] (relative to horizontal)
%   gamma      - Ratio of specific heats
%
% Outputs:
%   M2         - Downstream Mach number
%   theta2     - Downstream flow angle [degrees] (relative to horizontal)
%   nu2        - Downstream Prandtl-Meyer angle [degrees]
%   shockProps - Struct containing:
%                .Mn1       : Upstream normal Mach number
%                .Mn2       : Downstream normal Mach number
%                .theta_def : Deflection angle across shock [degrees]
%                .P02_P01   : Stagnation pressure ratio P02 / P01
%                .P2_P1     : Static pressure ratio P2 / P1
%                .T2_T1     : Static temperature ratio T2 / T1

    % Shock angle relative to oncoming flow vector
    beta_rel = beta_shock - theta1;
    
    % Upstream normal Mach number
    Mn1 = M1 * sind(beta_rel);
    
    % Check if shock condition is met (Mn1 > 1)
    if Mn1 <= 1.0001
        % No shock or sound wave limit
        M2 = M1;
        theta2 = theta1;
        if M2 >= 1
            pm = MoC_PM(M2, gamma);
            nu2 = pm.nu;
        else
            nu2 = 0;
        end
        shockProps.Mn1 = Mn1;
        shockProps.Mn2 = Mn1;
        shockProps.theta_def = 0;
        shockProps.P02_P01 = 1.0;
        shockProps.P2_P1 = 1.0;
        shockProps.T2_T1 = 1.0;
        return;
    end
    
    % 1. Deflection Angle theta_def (Theta-Beta-M relation)
    num = 2 * cotd(beta_rel) * (Mn1^2 - 1);
    den = M1^2 * (gamma + cosd(2 * beta_rel)) + 2;
    theta_def = atand(num / den);
    
    % Downstream flow angle (deflected upwards toward slip line)
    theta2 = theta1 + theta_def;
    
    % 2. Downstream Normal Mach Number Mn2
    gm1 = gamma - 1;
    gp1 = gamma + 1;
    Mn2_sq = (2 + gm1 * Mn1^2) / (2 * gamma * Mn1^2 - gm1);
    if Mn2_sq < 0
        Mn2_sq = 0.01;
    end
    Mn2 = sqrt(Mn2_sq);
    
    % 3. Downstream Total Mach Number M2
    sin_term = sind(beta_rel - theta_def);
    if sin_term <= 1e-4
        sin_term = 1e-4;
    end
    M2 = Mn2 / sin_term;
    
    % 4. Downstream Prandtl-Meyer Angle nu2
    if M2 >= 1.0
        pm = MoC_PM(M2, gamma);
        nu2 = pm.nu;
    else
        nu2 = 0.0;
    end
    
    % 5. Stagnation Pressure Ratio P02 / P01 (Rayleigh pitot / shock total pressure loss)
    term1 = ((gp1 / 2) * Mn1^2) / (1 + (gm1 / 2) * Mn1^2);
    term2 = (2 * gamma * Mn1^2 - gm1) / gp1;
    P02_P01 = (term1^(gamma / gm1)) * (term2^(-1 / gm1));
    
    % 6. Static Ratios
    P2_P1 = 1 + (2 * gamma / gp1) * (Mn1^2 - 1);
    T2_T1 = P2_P1 * ((2 + gm1 * Mn1^2) / (gp1 * Mn1^2));
    
    % Pack properties struct
    shockProps.Mn1 = Mn1;
    shockProps.Mn2 = Mn2;
    shockProps.theta_def = theta_def;
    shockProps.P02_P01 = P02_P01;
    shockProps.P2_P1 = P2_P1;
    shockProps.T2_T1 = T2_T1;
end

