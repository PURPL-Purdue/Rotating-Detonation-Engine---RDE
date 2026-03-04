
function [max_flux] = HADES_size_max_flux_range(total_mdot, Area_Annulus)

max_flux = total_mdot / Area_Annulus;

fprintf('Max Flux Range: %.5f kg/(s*mm^2)\n', max_flux)