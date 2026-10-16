function M = MoC_PM_Inv(nu_deg, gamma)
% MoC_PM_Inv Returns the Mach number for a given Prandtl-Meyer angle (in degrees).
% Uses an efficient Newton-Raphson solver for high performance.

    nu = nu_deg * pi / 180; % Convert to radians
    
    % Initialize Mach number arrays
    M = zeros(size(nu));
    
    % Handle Mach 1 exactly
    M(nu <= 0) = 1.0;
    
    % For nu > 0, solve using Newton-Raphson
    idx = nu > 0;
    if any(idx)
        nu_target = nu(idx);
        
        % Initial guess: a linear approximation M = 1 + c * nu
        M_guess = 1.0 + 1.2 * nu_target;
        
        tol = 1e-6;
        max_iter = 20;
        
        gp1 = gamma + 1;
        gm1 = gamma - 1;
        coeff1 = sqrt(gp1 / gm1);
        coeff2 = gm1 / gp1;
        
        for i = 1:max_iter
            M2 = M_guess.^2;
            M2_minus_1 = M2 - 1;
            
            % Check for invalid Mach numbers
            M2_minus_1(M2_minus_1 < 0) = 0.0001; 
            
            sq_M2_m1 = sqrt(M2_minus_1);
            
            % Prandtl-Meyer function
            nu_calc = coeff1 * atan(sqrt(coeff2 * M2_minus_1)) - atan(sq_M2_m1);
            
            % Derivative d(nu)/dM
            % Formula: d(nu)/dM = sqrt(M^2 - 1) / (M * (1 + (gamma-1)/2 * M^2))
            deriv = sq_M2_m1 ./ (M_guess .* (1 + (gm1 / 2) * M2));
            
            % Avoid divide by zero
            deriv(deriv == 0) = 1e-6;
            
            % Newton step
            delta_M = (nu_calc - nu_target) ./ deriv;
            M_guess = M_guess - delta_M;
            
            if max(abs(delta_M)) < tol
                break;
            end
        end
        
        M(idx) = M_guess;
    end
end

