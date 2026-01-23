clear; clc; close;
% Author: anjali
% vals from caltech database
% (https://shepherd.caltech.edu/detn_db/html/H2-Air2.html)

%inputs
data = readmatrix("HADES_size_cellsizedata.xlsx");
chosen_phi = 0.9;

phis = data(:, 1);
cell_width = data(:, 2);

log_phis = log10(phis);
log_cell_width = log10(cell_width);

coeffs = polyfit(log_phis, log_cell_width, 3);
log_phis_fit = linspace(min(log_phis), max(log_phis), 100);
log_cw_fit = polyval(coeffs, log_phis_fit);

figure();
scatter(log_phis, log_cell_width);
hold on;
plot(log_phis_fit, log_cw_fit);
xlabel("log(equivalence ratio)");
ylabel("log(detonation cell width)");
title("detonation cell width vs equivalence ratio");
grid on;


log_chosen_phi = log10(chosen_phi);
log_chosen_cellwidth = polyval(coeffs, log_chosen_phi);
chosen_cellwidth = 10 ^ log_chosen_cellwidth;

% annulus gap correlation from karasu continuous spin detonation
annulus_gap = 1.4 * chosen_cellwidth; 

fprintf("Chosen equivalence ratio: %.1f. \n", chosen_phi);
fprintf("Calculated Cell Width: %.3f mm \n", chosen_cellwidth);
fprintf("Annulus Gap: %.3f mm \n", annulus_gap);

chosen_annulus_gap = 9.575; %mm
%for Hades chose annulus gap of 3/8 in or 9.575 mm 


%% fill height

