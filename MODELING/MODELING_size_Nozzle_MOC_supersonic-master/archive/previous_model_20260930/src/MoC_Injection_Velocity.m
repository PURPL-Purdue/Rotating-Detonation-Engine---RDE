function out = MoC_Injection_Velocity(inj,chamber)
% Adiabatic momentum/area-change model; physical subsonic root only.
fields = {'P3','T3','mdot','Pa','R','gamma1','N_holes','injector_area'};
for k=1:numel(fields)
    validateattributes(inj.(fields{k}),{'numeric'},{'scalar','real','finite','positive'});
end
if inj.gamma1<=1, error('MoC:Injection','gamma must exceed one.'); end
validateattributes(inj.N_holes,{'numeric'},{'integer'});
if isfield(inj,'port_diameter') && ~isempty(inj.port_diameter)
    validateattributes(inj.port_diameter,{'numeric'},{'scalar','positive','finite'});
    inj.injector_area = pi*inj.port_diameter^2/4;
end
MoC_Calculate_Domain_Size(chamber);
A3 = inj.injector_area*inj.N_holes*1e-6;
A4 = pi/4*(chamber.outer_chamber_diameter^2-chamber.inner_chamber_diameter^2)*1e-6;
if A3>=A4, error('MoC:Injection','Total port area must be smaller than annulus area.'); end
cp = inj.gamma1*inj.R/(inj.gamma1-1);
V3 = inj.mdot*inj.R*inj.T3/(inj.P3*A3);
T03 = inj.T3+V3^2/(2*cp);
flux = inj.mdot/A4;
alpha = (inj.P3-inj.Pa)*A3/A4+inj.Pa+flux*V3;
a = flux*(inj.R/(2*cp)-1); b=alpha; c=-flux*inj.R*T03;
D = b^2-4*a*c;
if D<0, error('MoC:Injection','No real injector momentum/energy solution.'); end
q = -0.5*(b+sqrt(D)); % stable quadratic
V = [q/a;c/q]; P=alpha-flux*V; T=T03-V.^2/(2*cp);
valid = V>0 & P>0 & T>0;
M = nan(size(V)); M(valid)=V(valid)./sqrt(inj.gamma1*inj.R*T(valid));
index = find(valid & M<1);
if numel(index)~=1
    error('MoC:Injection','Expected one physical subsonic injection root; found %d.',numel(index));
end
k=index(1);
out = struct('V',V(k),'P',P(k),'T',T(k),'M',M(k), ...
    'rho',P(k)/(inj.R*T(k)),'gamma',inj.gamma1,'R',inj.R,'cp',cp, ...
    'Aports',A3,'Aannulus',A4,'Vport',V3,'roots',[V P T M]);
out.massResidual = out.rho*out.V*A4-inj.mdot;
out.energyResidual = cp*out.T+out.V^2/2-cp*T03;
end
