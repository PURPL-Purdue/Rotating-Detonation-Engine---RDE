function out = MoC_Velocity_Wave_Calculations(InjV, ChamberDim)

% Use SI units; zetaDeg is in degrees.
% Uses Blake's "Derivation of Momentum Flux-Area Change Relations"

P3 = InjV.P3;
T3 = InjV.T3;
mdot = InjV.mdot;
Pa = InjV.Pa;
R = InjV.R;
Cp = InjV.Cp;
N_holes = InjV.N_holes;

% Converting to SI and calling ChamberDim variables to calculate A4
A3 = InjV.injector_area * N_holes * 1e-6; % Injector orifice area in from mm^2 to m^2
InjV.A4 = (pi / 4) * (ChamberDim.outer_chamber_diameter^2 - ChamberDim.inner_chamber_diameter^2) * 1e-6;
A4 = InjV.A4;

V3 = mdot*R*T3/(P3*A3);

T03   = T3 + V3^2/(2*Cp);
beta  = mdot/A4;
alpha = (P3 - Pa)*A3/A4 + Pa + beta*V3;

a = beta*(R/(2*Cp) - 1);
b = alpha;
c = -beta*R*T03;

D = b^2 - 4*a*c;
if D < 0
    error('No real velocity solutions for these inputs.');
end

out.V4 = [(-b + sqrt(D))/(2*a);
          (-b - sqrt(D))/(2*a)];

out.P4 = alpha - beta*out.V4;
out.T4 = T03 - out.V4.^2/(2*Cp);
out.Vwave = out.V4;

% NaN marks a root with nonpositive velocity, pressure, or temperature.
invalid = out.V4 <= 0 | out.P4 <= 0 | out.T4 <= 0;
out.Vwave(invalid) = NaN;

disp('Root     V4 (m/s)       P4 (Pa)        T4 (K)      Vwave (m/s)')
for k = 1:2
    fprintf('%d     %11.3f   %12.3f   %11.3f   %13.3f\n', ...
        k, out.V4(k), out.P4(k), out.T4(k), out.Vwave(k));
end
end