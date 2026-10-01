%testing
% Output structure:


d = HADES_size_ceaDet('ox','Air','fuel','H2','phi',0.9,'P0',26,'P0Units','psia','T0',283,'T0Units','K');

r_i_in = 2; %internal radius
r_o_in = 2.25; %outer radius 
p_o_MPa = 0.101325; %ambient pressure
P_cj_bar = 30.7; %CJ pressure

res = HADES_size_HoopStressTemps(r_i_in, r_o_in, p_o_MPa, P_cj_bar);
