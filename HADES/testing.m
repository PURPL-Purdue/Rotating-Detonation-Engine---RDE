clc; clear;

d = ceaDet('ox','O2','fuel','H2','phi',1,'P0',1000,'P0Units','psia','T0',300,'T0Units','K' );

fprintf('\nCEA DETONATION RESULTS\n');

fprintf('Detonation Mach Number = %.4f\n', d.detMach);
fprintf('Detonation Velocity = %.2f m/s\n', d.detVel);
fprintf('P/P1 = %.4f\n', d.P_ratio);
fprintf('T/T1 = %.4f\n', d.T_ratio);
fprintf('Burned Pressure = %.2f bar\n', d.P_burned_bar);
fprintf('Burned Temp = %.2f K\n', d.T_burned_K);
