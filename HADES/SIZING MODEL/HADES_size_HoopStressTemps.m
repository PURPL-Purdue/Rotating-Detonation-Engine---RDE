function [Fail_temp_C, Fail_temp_F] = HADES_size_HoopStressTemps (r_i_in, r_o_in, p_o_MPa, P_cj_bar, fos)
% Find failure temperature based on cj Pressure and hoop stress
% calculations
% 
% Example Usage
% HADES_size_HoopStressTemps( 2, 2.5, 0.101325, 30.7);
%
% Inputs
% r_i_in - inner radius in inches
% r_o_in - outer radius in inches
% p_o_Mpa - outside pressure in Mpa
% P_cj_bar - CJ pressure in bar
% fos - factor of safety
%
% Outputs
% Fail_temp_C - failure temperature in celcius
% Fail_temp_F

%% Constants
% Temperature data in degrees Fahrenheit
temp_F = [1500, 1440, 1330, 1240, 1180, 1140, 1100, 1020, ...
          980, 922, 835, 719, 661, 589, 502, 459, 401, ...
          357, 285, 242, 213, 184, 155, 126, 68];

yeild_points_ksi = [25.8, 29.5, 37.4, 42.7, 45.9, 48.1, 49.9, 52.6, ...
              53.7, 54.8, 55.7, 56.3, 56.3, 56.6, 57, 57.4, 58.3, ...
              59.2, 62.1, 64, 65.5, 67.5, 69.8, 72.5, 79] ./  fos;

temp_points_K = ((temp_F - 32) .* 5/9) + 273.15;

yeild_points_MPa = yeild_points_ksi .* 6.89476;


temp_points_C = temp_points_K - 273.15;


%% Calculations

%convert to mm for calculations
r_i_mm = r_i_in * 25.4;
r_o_mm = r_o_in * 25.4;
r_mm = r_i_mm;

u_mm = (r_o_mm ^ 2) - (r_i_mm ^ 2);
a_mm = (r_i_mm ^ 2 / u_mm) + (r_i_mm ^ 2 * r_o_mm ^ 2 / (r_mm ^ 2 * u_mm));
b_mm = (r_o_mm ^ 2 / u_mm) + (r_i_mm ^ 2 * r_o_mm ^ 2 / (r_mm ^ 2 * u_mm));

p_i_MPa = ((p_o_MPa * b_mm) + yeild_points_MPa) ./ a_mm;
p_i_psi = p_i_MPa * 145.03773773;
P_cj_psi = P_cj_bar * 14.5038;

%Fitting function
p = polyfit(temp_points_C, p_i_psi, 3);

temp_points_ideal_C = linspace(min(temp_points_C), max(temp_points_C), 100);
p_i_ideal_psi = polyval(p, temp_points_C);
p_i_curve_psi = polyval(p, temp_points_ideal_C);
p_i_psi_bar = mean(p_i_psi);

SSE = sum((p_i_psi - p_i_ideal_psi) .^ 2);
SST = sum((p_i_psi - p_i_psi_bar) .^ 2);
R_squared = 1 - (SSE/SST);

fprintf("The R squared value is %.3f\n", R_squared);

%find intersection point using fzero method
intersection_func = @(T) polyval(p, T) - P_cj_psi;

T_iniial_guess_C = mean(temp_points_C);
T_intersection_C = fzero(intersection_func, T_iniial_guess_C);

fprintf("Failure Temperature occurs at %.2f C\n", T_intersection_C);
fprintf("Failure Temperature occurs at %.2f F\n", (T_intersection_C * 1.8) + 32);

Fail_temp_C = T_intersection_C;
Fail_temp_F = (T_intersection_C * 1.8) + 32;


%% Plots

figure;
hold on;
plot((temp_points_C * 1.8) + 32, p_i_psi, "r*");
plot((temp_points_ideal_C * 1.8) + 32, p_i_curve_psi, LineWidth=1.6,Color = [0.5 0 0.8]);
yline(P_cj_psi,"-", "CJ Presure");
title("Maximum Pressure of SS316 vs Temperature")
xlabel("Temperature [F]");
ylabel("Allowed Chamber Pressure [psi]")
grid on;

