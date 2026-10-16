function output = MoC_PM(Mach, gamma)
% Evaluate the Prandtl-Meyer function if your formulation requires it.
% mu is the Mach angle in degrees
% nu is the Prandtl-Meyer function in degrees

% Returns the Prandtl-Meyer function nu and the Mach angle mu
  % Assert ensures Mach number is greater than 1 (supersonic), and throws an error if not
  assert(Mach>=1, 'Error: Mach is under 1 (subsonic)'); 
  mu = asind(1./Mach);
  
  sqroot = sqrt((gamma+1)/(gamma-1));
  nu = sqroot * atand(sqrt(Mach.^2-1)/sqroot) - atand(sqrt(Mach.^2-1));
  
  output.nu = nu;
  output.mu = mu;
end