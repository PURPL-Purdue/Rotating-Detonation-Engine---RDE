
function [fake_plenum_pressure] = PURPL_Fake_Plenum_pressure()

base_plenum_pressure = 20; %20 bar expected static pressure
base_plenum_pa = base_plenum_pressure * 100000; %base static pressure
f_samp = 1000; %sampling frequency
dt  = 1/f_samp; %time step
time = (0: f_samp-1) * dt; %time vector containing data points for one second

%Below some random frequencies are defined that will cause minor
%oscillations in the plenum pressures
f1 = 300 + (600 - 300) * rand;
f2 = 600 + (1200 - 600) * rand;
f3 = 1000 + (1300 - 1000) * rand;

%sets the baseline pressure for the hydrogen and air plenums
pressure_h2 = ones(1, length(time)) * base_plenum_pa;
pressure_air = ones(1, length(time)) * base_plenum_pa;

%adds oscillations that are scaled for each frequency
plenum_pressure_h2 = pressure_h2 + abs((0.003 * pressure_h2 .* sin(time .* f1)) + (0.0025 * pressure_h2 .* sin(time .* f2)) + (0.0015 * pressure_h2 .* sin(time .* f2)));
plenum_pressure_air = pressure_air + abs((0.002 * pressure_air .* sin(time .* f1 * 2)) + (0.00095 * pressure_air .* sin(time .* f2 * 3)) + (0.00045 * pressure_air .* sin(time .* f2 * 1.5)));


%combines outputs into one matrix
fake_plenum_pressure = [time; plenum_pressure_h2; plenum_pressure_air];

%plots the outputs
figure(1)
plot(time, plenum_pressure_h2/100000)
hold on
plot(time, plenum_pressure_air/100000)
ylabel('Pressure (bar)')
xlabel('Time (s)')
title('Static Plenum Pressure Over Time')
legend('Hydrogen Plenum', 'Air Plenum')
hold off
grid on

