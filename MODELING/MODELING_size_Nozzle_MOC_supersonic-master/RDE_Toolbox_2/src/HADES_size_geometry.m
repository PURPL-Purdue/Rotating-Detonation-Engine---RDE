function [det_wave_path_length] = HADES_size_geometry(outer_r, wall_thickness, annulus_gap)

%calcs
annulus_gap_in = annulus_gap / 25.4; %convert to in (is inputted as mm)
inner_r = outer_r - annulus_gap_in;
annulus_r = (outer_r + inner_r) / 2;
det_wave_path_length = 2 * pi() * annulus_r;

% Print Values
fprintf("\nGeometry:\n");
fprintf("Outer Radius: %.3f in.\n", outer_r);
fprintf("Wall Thickness: %.3f in.\n", wall_thickness);
fprintf("Annulus Gap: %.3f mm (%.3f in). \n", annulus_gap, annulus_gap_in);
fprintf("Inner Radius: %.3f in.\n", inner_r);
fprintf("Annulus Radius: %.3f in.\n", annulus_r);
fprintf("Detonation Wave Path Length: %.3f in.\n", det_wave_path_length);