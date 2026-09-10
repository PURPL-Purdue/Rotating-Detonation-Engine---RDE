%% CONSTANTS
x1 = 0.25;         % sensor 1 depth
x2 = 0.575;        % sensor 2 depth

k = 390;           % thermal conductivity (W/m-K)
alpha = 1.11e-4;   % thermal diffusivity (m^2/s)

sgolay_order = 2;
sgolay_window = 21;

%% DATA READING
data = readmatrix('/Users/Abby/Desktop/RDE/sample_thermocouple_data.csv');

t  = data(:,1);
T1 = data(:,2);
T2 = data(:,3);

dt = mean(diff(t));

%% SAVISTSKY-GOLAY FILTER
T1_s = sgolayfilt(T1, sgolay_order, sgolay_window);
T2_s = sgolayfilt(T2, sgolay_order, sgolay_window);

%% TEMPERATURE DERIVATIVES
T1_dot = T1_s / dt;
T2_dot = T2_s / dt;

%% SURFACE TEMPERATURE: EQ. 3
T0 = ( (x2/(x2-x1)) * T1_s ...
     - (x1/(x2-x1)) * T2_s ...
     + alpha * ((x2 * (x2-2*x1))/(6*(x2-x1))) * T1_dot ...
     + alpha * ((x1 * (2*x2-x1))/(6*(x2-x1))) * T2_dot );

%% HEAT FLUX: EQ. 4
q0 = -k * ( ...
      (T2_s - T1_s)/(x2-x1) ...
    + alpha * ((x2^2 + x1 * x2 - 2 * x1^2) / (6 * (x2-x1))) * T1_dot ...
    + alpha * ((2 * x2^2 - x1 * x2 - x1^2) / (6 * (x2-x1))) * T2_dot );

%% --- PLOTS ---
figure
plot(t,T1,'b',t,T2,'r')
xlabel('Time (s)')
ylabel('Temperature (K)')
legend('T1','T2')
title('Measured Thermocouple Temperatures')

figure
plot(t,T0,'c')
xlabel('Time (s)')
ylabel('Surface Temperature (K)')
title('Calculated Surface Temperature')

figure
plot(t,q0)
xlabel('Time (s)')
ylabel('Heat Flux (W/m^2)')
title('Calculated Surface Heat Flux')