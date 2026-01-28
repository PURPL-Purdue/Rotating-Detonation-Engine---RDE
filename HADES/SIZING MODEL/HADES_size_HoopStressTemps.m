function [Fail_temp_C, Fail_temp_F] = HADES_size_HoopStressTemps (r_i_in, r_o_in, p_o_MPa, P_cj_bar)


%% Constants
temp_points_K = [
    300.15, 
    422.15,  
    533.15,  
    700.15,  
    755.15, 
    866.15, 
    977.15,  
    1088.15, 
    1200.15, 
    1311.15, 
    1366.15,  
    1648.00
];

yeild_points_MPa = [
    290,
    201,
    172,
    159,
    148,
    140,
    131,
    110,
    80, 
    39, 
    28,  
    1.0
];

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
p = polyfit(temp_points_C, p_i_psi, 2);

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

