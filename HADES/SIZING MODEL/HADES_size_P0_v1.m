function [P_0] = HADES_size_P0_v1()

Ai_Aw = .1135; % unitless
M_0 = 0.65; % unitless
P_p = 20; % bar
gamma_unburned = 1.4; %unitless
P_c = (2/(gamma_unburned+1))^(gamma_unburned/(gamma_unburned-1)) * P_p; %bar
P_0 = Ai_Aw * (1/M_0) * ((gamma_unburned + 1)/(2+(gamma_unburned-1) * M_0^2))^(1/2) * P_c; %bar


