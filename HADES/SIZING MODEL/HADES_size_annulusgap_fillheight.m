function [chosen_cellwidth, chosen_annulus_gap, chosen_fillheight] = HADES_size_annulusgap_fillheight(chosen_phi, L_theta, N_det, D_CJ, a2, P2, Pc, u0)
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

chosen_annulus_gap = 9.525; %mm
%for Hades chose annulus gap of 3/8 in or 9.525 mm

fprintf("Calculated Cell Width: %.3f mm \n", chosen_cellwidth);
fprintf("Calculated Annulus Gap: %.3f mm \n", annulus_gap);



%% fill height

%fill_HEIGHT Compute detonation wave height h_det from Kawashima Eq. (15).
%
%   fill_h = detonation_height(L_theta, N_det, D_CJ, a2, P2, Pc, u0)
%
%   Inputs:
%     L_theta : azimuthal annulus length [m]
%     N_det   : detonation wave number [-]
%     D_CJ    : Chapman-Jouguet detonation speed [m/s]
%     a2      : sound speed in burned gas [m/s]
%     P2      : pressure of burned gas region [Pa or bar, consistent with Pc]
%     Pc      : injector critical pressure [same units as P2]
%     u0      : axial velocity of unburned propellant [m/s]
%
%   Output:
%     h_det   : detonation (fill) height [m]

%Getting value for function (eq 9 in paper)
 x = P2 ./ Pc;

 f_val = (231/2^10)*x.^(1/15) + ...
            ( 63/ 2^9)*x.^(3/15) + ...
            (105/2^10)*x.^(5/15) + ...
            ( 25/ 2^8)*x.^(7/15) + ...
            (105/2^10)*x.^(9/15) + ...
            ( 63/ 2^9)*x.^(11/15) + ...
            (231/2^10)*x.^(13/15);

% Calculations for eq15

% Calculations for the bracket
bracket = (1./a2) .* (f_val) + 1./u0;

chosen_fillheight = (L_theta ./ (N_det .* D_CJ)) .* (1 ./ bracket);

% Display the calculated fill height
fprintf("Fill Height: %.3f m (%.3f mm). \n", chosen_fillheight, chosen_fillheight * 1000);



