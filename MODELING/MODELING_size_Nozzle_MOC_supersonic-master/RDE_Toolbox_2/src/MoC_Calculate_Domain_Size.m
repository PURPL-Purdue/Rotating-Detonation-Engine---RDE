function [xMax, yMax] = MoC_Calculate_Domain_Size(ChamberDim)
    
% Convert the chamber dimensions in mm into an unwrapped domain.

OC_dia = ChamberDim.outer_chamber_diameter * 1e-3;
IC_dia = ChamberDim.inner_chamber_diameter * 1e-3;
Chamb_L = ChamberDim.chamber_length * 1e-3;

% Calculate annulus circumference and define the unwrapped domain.

annulusMeanDiameter = (OC_dia + IC_dia) / 2;
annulusCircumference = pi * annulusMeanDiameter;
xMax = annulusCircumference * 1000;
yMax = ChamberDim.chamber_length;


end
