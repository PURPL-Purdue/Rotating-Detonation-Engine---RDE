function [fake_pressure_data] = PURPL_Fake_ITP_Pressure_Data()
%% High-Magnitude Equalized Correlated Pressure Data with Time Delay
n_points = 100000;
fs = 100000;           % 100 kHz (10 microseconds per sample)
dt = 1/fs;
t = (0:n_points-1)' * dt;

%% 1. Define Equal-Magnitude Dominant Frequencies
f_dom = 1000 + (2000 - 1000) * rand;
f2 = 2000 + (4000 - 2000) * rand;
f3 = 4000 + (8000 - 4000) * rand;

% The "Source" signal
shared_signal = 300 * (sin(2*pi*f_dom*t) + sin(2*pi*f2*t) + sin(2*pi*f3*t));

%% 2. Apply Time Delay to the second dataset
% Let's delay it by 3 samples (30 microseconds)
delay_samples = 3; 

% Sensor 1 uses the original signal
phys_signal1 = shared_signal;

% Sensor 2 uses a shifted version of the signal
% We pad the beginning with zeros and truncate the end to keep the length same
phys_signal2 = [zeros(delay_samples, 1); shared_signal(1:end-delay_samples)];

%% 3. Normalize and Scale to PSI (30-70 range)
% We normalize based on Sensor 1 to keep them consistent
s_min = min(phys_signal1);
s_max = max(phys_signal1);

norm1 = (phys_signal1 - s_min) / (s_max - s_min);
norm2 = (phys_signal2 - s_min) / (s_max - s_min);

pressure1 = norm1 * 170 + 30 + 0.5*randn(n_points,1);
pressure2 = norm2 * 188 + 30 + 0.5*randn(n_points,1); % Slightly different scaling

%% 4. Add Synchronized Spikes (also delayed)
num_events = 150;
for i = 1:num_events
    idx = randi([1 n_points-10]);
    width = randi([1 3]);
    mag = 1000;
    pulse = mag * exp(-3*(0:width-1)');
    L = length(pulse);
    
    % Sensor 1 gets it at idx
    pressure1(idx:idx+L-1) = pressure1(idx:idx+L-1) + pulse;
    
    % Sensor 2 gets it at idx + delay
    idx2 = idx + delay_samples;
    if idx2 + L - 1 <= n_points
        pressure2(idx2:idx2+L-1) = pressure2(idx2:idx2+L-1) + 0.95*pulse;
    end
end

fake_pressure_data = [t pressure1 pressure2];

figure(2)
plot(t,pressure1)
hold on
plot(t,pressure2)
xlabel('Time (t)')
ylabel('Pressure (kPa)')
legend('ITP 1', 'ITP 2', 'Location', 'best')
title('Fake ITP Data')

fprintf("Dominant Frequency: %.5f (kHz)\n", f_dom/1000)
fprintf("Prominent Frequency 2: %.5f (kHz)\n", f2/1000)
fprintf("Prominent Frequency 3: %.5f (kHz)\n", f3/1000)

end