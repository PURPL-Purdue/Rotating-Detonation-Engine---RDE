function [xMax, yMax] = MoC_Calculate_Domain_Size(ChamberDim)

% Calculates 2D unwrapped domain of engine from given chamber size values
% and outputs the domain of the graph that will be used to display the sim.  


% Convert all variables to SI units (meters)

OC_dia = ChamberDim.outer_chamber_diameter * 1e-3;
IC_dia = ChamberDim.inner_chamber_diameter * 1e-3;
Chamb_L = ChamberDim.chamber_length * 1e-3;

% Calculate annulus circumference and define the unwrapped domain.

annulusMeanDiameter = (OC_dia + IC_dia) / 2;
annulusCircumference = pi * annulusMeanDiameter;
xMax = annulusCircumference;
yMax = Chamb_L;

end