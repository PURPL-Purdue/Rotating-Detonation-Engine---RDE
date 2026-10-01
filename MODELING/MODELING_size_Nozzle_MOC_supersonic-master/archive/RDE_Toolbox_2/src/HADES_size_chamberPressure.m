function [avg_chamber_p] = HADES_size_chamberPressure(wave_mode, detwave_pathlength, cj_vel_ms, p_p1, p_0)

cj_vel = HADES_size_convert(cj_vel_ms,'m/s','ft/s');

cycle_time = (detwave_pathlength /(cj_vel*12)) * 10 ^ 3 / wave_mode;
exp_decay_const = log(p_p1)/(cycle_time/1000);

time = 0:0.001:cycle_time;

chamber_pressure = p_p1 * p_0 * exp (-1 * exp_decay_const * time ./ 1000); %psia
avg_chamber_p = mean(chamber_pressure); %psia

% print sol and plot
fprintf('\nAverage Chamber Pressure: %.2f psia (%.3f bar). \n', avg_chamber_p, avg_chamber_p / 14.504);

figure();
plot(time, chamber_pressure);
title("Chamber Pressure vs Time");
xlabel("Time (ms)");
ylabel("Chmaber Pressure (psia)");
grid on;
