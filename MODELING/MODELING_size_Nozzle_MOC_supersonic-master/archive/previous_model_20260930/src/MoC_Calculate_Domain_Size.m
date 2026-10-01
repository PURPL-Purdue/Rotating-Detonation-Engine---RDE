function [circumference,height] = MoC_Calculate_Domain_Size(c)
% Input dimensions mm; output meters.
validateattributes(c.outer_chamber_diameter,{'numeric'},{'scalar','finite','positive'});
validateattributes(c.inner_chamber_diameter,{'numeric'},{'scalar','finite','positive'});
validateattributes(c.chamber_length,{'numeric'},{'scalar','finite','positive'});
if c.outer_chamber_diameter<=c.inner_chamber_diameter
    error('MoC:Geometry','Outer diameter must exceed inner diameter.');
end
circumference = pi*(c.outer_chamber_diameter+c.inner_chamber_diameter)*0.5e-3;
height = c.chamber_length*1e-3;
end
