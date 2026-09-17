function [data] = PURPL_ITP_data_fft_new()

clc;
close all;

% Parameters
fs = 100000; % Sampling frequency (100 kHz to satisfy Nyquist for >10 kHz signals)
duration = 1; % Total duration in seconds
n_points = fs * duration; % 100,000 total points
t = (0:n_points-1)' / fs; % Time vector

% 1. Pick 3 random dominant frequencies between 10,000 Hz and 25,000 Hz
f_min = 10000;
f_max = 25000;
freqs = f_min + (f_max - f_min) * rand(1, 3);

% Random phases (0 to 2*pi)
phases1 = 2*pi * rand(1, 3);
phases2 = 2*pi * rand(1, 3);

% Random relative amplitudes for the 3 tones
amp1 = 0.2 + 0.8 * rand(1, 3);
amp2 = 0.2 + 0.8 * rand(1, 3);

% 2. Combine high-frequency tones with noise
raw_signal1 = amp1(1)*sin(2*pi*freqs(1)*t + phases1(1)) + ...
              amp1(2)*sin(2*pi*freqs(2)*t + phases1(2)) + ...
              amp1(3)*sin(2*pi*freqs(3)*t + phases1(3)) + ...
              0.1 * randn(n_points, 1);

raw_signal2 = amp2(1)*sin(2*pi*freqs(1)*t + phases2(1)) + ...
              amp2(2)*sin(2*pi*freqs(2)*t + phases2(2)) + ...
              amp2(3)*sin(2*pi*freqs(3)*t + phases2(3)) + ...
              0.1 * randn(n_points, 1);

% 3. Scale signals into the [1, 5] bar range
scale_to_range = @(x, target_min, target_max) ...
    target_min + (x - min(x)) * (target_max - target_min) / (max(x) - min(x));

pressure1 = scale_to_range(raw_signal1, 1.0, 5.0);
pressure2 = scale_to_range(raw_signal2, 1.0, 5.0);

% 4. Save to CSV
data = [t, pressure1, pressure2];
filename = 'synthetic_pressure_data.csv';
writematrix(data, filename);


% Plot Frequencies over time (confirm randomness)
figure(1);
plot(data(:,1), data(:,2), 'DisplayName', 'Sensor 1');
hold on;
plot(data(:,1), data(:,3), 'DisplayName', 'Sensor 2');
grid on;
xlabel('Time (s)');
ylabel('Pressure (bar)');
title('Fake ITP Pressure Data');
legend('Location', 'best');

end