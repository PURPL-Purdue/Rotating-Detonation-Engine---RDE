


C_d_h2  = 0.71; %Cd for hydrogen injector
C_d_air = 0.62; %Cd for air injector

temp = 293.15; %room temperature (293.15K)

R_air = 287.05; %(J/kgK)
R_h2 = 4124.2; %(J/kgK)

A_h2 = pi * (0.001/2) ^2 * 32; %area of the hydrogen injetcor (m^2)
A_air = 304.77 / 1000 /1000; %area of the air injector (m^2)

gamma_h2 = 1.405; %specific heat capacity for hydrogen
gamma_air = 1.4; %specific heat capacity for air

stag_density_h2 = 0.08988; %stagnation density of hydrogen at 293K(kg/m^3)
stag_density_air  = 1.225; %density of air (kg/m^3)


fake_plenum_pressure = PURPL_Fake_Plenum_pressure(); %imports plenum pressure data 

time = fake_plenum_pressure(1,:); %renames and isolates time vector
air_plenum_pressure = fake_plenum_pressure(2,:); %renames and isolates air pressure vector
h2_plenum_pressure = fake_plenum_pressure(3,:); %renames and isolates hydrogen pressure vector

density_air = air_plenum_pressure ./ (R_air * temp); %calculates the denisty of air using ideal gas law
density_h2  = h2_plenum_pressure ./ (R_h2 * temp); %calculates the denisty of hydrogen using ideal gas law


%uses choked mass flow equation to calculate mass flow rate from pressure 
mass_flow_air = C_d_air .* A_air .* sqrt(gamma_air .* density_air .* air_plenum_pressure .* ((2 / (gamma_air + 1)) ^ ((gamma_air + 1) / (gamma_air - 1))));
mass_flow_h2 = C_d_h2 .* A_h2 .* sqrt(gamma_h2 .* density_h2 .* h2_plenum_pressure .* ((2 / (gamma_h2 + 1)) ^ ((gamma_h2 + 1) / (gamma_h2 - 1))));


%plots mass flow rate in lbs/2
figure(2)
plot(time, mass_flow_h2 .* 2.20462)
hold on
plot(time, mass_flow_air .* 2.20462)
xlabel('Time (s)')
ylabel('Mass Flow Rate (lbs/s)')
title('Static Plenum Pressure Over Time')
legend('Hydrogen Plenum', 'Air Plenum')
hold off
grid on