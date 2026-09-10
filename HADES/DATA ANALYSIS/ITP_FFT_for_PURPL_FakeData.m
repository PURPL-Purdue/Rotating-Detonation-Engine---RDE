function [] = ITP_FFT_for_PURPL_FakeData()

%imports fake pressure data
fake_pressure_data = PURPL_Fake_ITP_Pressure_Data();

%converts the data matrix into individual columns for easier use later
time = fake_pressure_data(:,1); %time data
itp_pt_1 = fake_pressure_data(:,2); %first ITP data
itp_pt_2 = fake_pressure_data(:,3); %second ITP data

dt = time(2) - time(1);  %Time step
fs = 1/dt;               %sampling frequency
f = linspace(1,fs,fs+1) /1000; %frequencies vector
n = length(time); %length of time term or number of samples

% Compute the magnitude of each FFT
fft1 = abs(fft(itp_pt_1));
fft2 = abs(fft(itp_pt_2));

% Average them
combined_fft = (fft1 + fft2) / 2;


fft_final = zeros(1, length(combined_fft) - 1);

%shifts the FFT back one entry (trust me its needed)
for index = 1: (length(combined_fft) - 1)
    fft_final(index) = combined_fft(index + 1);
end


% Plots the single-sided spectrum
figure(1)
plot(f(1:n/2), fft_final(1:n/2))
title('Fast Fourier Transform of Data set')
ylabel('|a_{k}|')
xlabel('Frequency (kHz)')










