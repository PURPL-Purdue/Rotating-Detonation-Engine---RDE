function [thrust] = HADES_size_thrust(isp, total_mdot_lb)

total_mdot = total_mdot_lb * 0.453592; % Convert total mass flow rate from lb/s to kg/s
thrust = isp * total_mdot;

fprintf("Thrust: %.3f N (%.3f lbf) \n", thrust, thrust / 0.453592);