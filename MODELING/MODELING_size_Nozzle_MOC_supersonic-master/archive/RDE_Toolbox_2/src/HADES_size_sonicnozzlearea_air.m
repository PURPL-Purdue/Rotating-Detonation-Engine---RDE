clc;

%%NOT DONE NEED TO INTEGRATE TO MAIN
air_mdot_gay = 2; %(lb/s)
air_mdot = air_mdot_gay / 2.205; % kg/s

% Guesses
p_feed_gay = 20:1:60; % feed pressue (bar)
p_feed = p_feed_gay * 100000; % Pa
density_air = 1.225; %m^3/kg

c_d = 0.9;
c_f = 0.6989; % critical flow
R = 8.134; % J/mol * K
M = 0.02897; %kg/mol
T_0 = 283.15; %K

d_feed_in = 0.61;%in
d_feed = 0.61 * 25.4; %mm
A_feed_gay = pi * (d_feed / 2) ^2; % (mm^2) -16 : assuming wall thickness of 0.095"
A_feed = A_feed_gay / 10^6; % m^2

A_choke_better = air_mdot .* sqrt(R ./ M .* T_0) ./ (c_d .* c_f .* p_feed);
dia_mm = 2 * sqrt(A_choke_better * 10^6 / pi());

%fprintf("The area of choke is %.3f m^2 (%.3f mm^2)\n", A_choke_better, A_choke_better * 10^6);
%fprintf("Area Ratio: %.3f\n", A_feed./A_choke_better);
%fprintf("Diameter = %.3f mm ", dia_mm);

figure();
plot(p_feed_gay, dia_mm);
xlabel('Feed Pressure (bar)');
ylabel('Choke Diameter (mm)');
yline(11,'r-', "Minimum Diameter");
title('Choke Diameter vs Feed Pressure');
grid on;

