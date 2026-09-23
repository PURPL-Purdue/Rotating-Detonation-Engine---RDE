function [] = PURPL_test_FFT_v2()

clc;
close all;

[fake_data] = PURPL_test_fake_ITP_data_v2(); %Calls fake ITP Data

%Assigns vector names for each column in the data matrix
time = fake_data(:,1);
pressure_1  = fake_data(:,2);
pressure_2 = fake_data(:,3);

time_step = time(2)-time(1); %finds the time step between samples (1/100000 s)
num_samples = length(time); %Counts the number of samples
sampling_freq = 1/time_step; %Finds the sampling frequency (Should be 100k hz)

f_1 = fft(pressure_1); %Runs FFT on first pressure data set
f_2 = fft(pressure_2); %Runs FFT on second pressure data set

p1_two_sided = abs(f_1/num_samples); %creates the two sided fft on the first data set
p2_two_sided = abs(f_2/num_samples); %creates one sided fft on the second data set

half_N = floor(num_samples / 2); %finds half the number of samples
p1_single_sided = p1_two_sided(1 : half_N + 1); %takes only the first half of the two sided fft
p2_single_sided = p2_two_sided(1 : half_N + 1); %takes only the first half of the two sided fft


p1_single_sided(2:end-1) = 2 * p1_single_sided(2:end-1); %doubles the strength of the one sided fft to account for lost power from only taking positive frequencies
p2_single_sided(2:end-1) = 2 * p2_single_sided(2:end-1); %doubles the strength of the one sided fft to account for lost power from only taking positive frequencies

%averages the fft and deletes the first point (always super high for some
%reason at 0 frequency)
average_fft = (p1_single_sided + p2_single_sided) ./ 2;
average_fft(1) = 0;

all_f = (0:half_N)' * (sampling_freq / num_samples); %gets all the frequencies found in a vector


%plots the fft
figure(2)
plot(all_f , average_fft, 'LineWidth', 1.25)
grid on
xlabel('Frequencies (Hz)')
ylabel('Relative Strength')
title('FFT Plot')
