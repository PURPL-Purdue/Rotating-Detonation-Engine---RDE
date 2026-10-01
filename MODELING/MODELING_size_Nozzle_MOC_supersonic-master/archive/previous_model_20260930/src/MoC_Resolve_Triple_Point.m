function t = MoC_Resolve_Triple_Point(p,inj,geom,seedMach)
% Scalar weak-shock match, Grunenwald Eqs. 2.22--2.28.
% Compatibility invariants are used for matching; no centered fan is seeded.
% Closure: uniform bounding products at injection pressure with the same
% stagnation enthalpy/entropy as the regularized post-detonation seed.
g=p.gamma; f=1+(g-1)*seedMach^2/2;
t.gas=struct('gamma',g,'R',p.R,'T0',p.T*f,'P0',p.P*f^(g/(g-1)));
if inj.P>=p.P, error('MoC:Triple','Product pressure must exceed refill pressure.'); end
M5=sqrt(2/(g-1)*((t.gas.P0/inj.P)^((g-1)/g)-1));
mu=asin(1/M5);
beta=linspace(mu+1e-7,pi/2-1e-7,500);
delta=atan(2*cot(beta).*(M5^2*sin(beta).^2-1)./(M5^2*(g+cos(2*beta))+2));
[~,peak]=max(delta); beta=beta(1:peak);
res=arrayfun(@residual,beta);
k=find(res(1:end-1).*res(2:end)<=0,1);
if isempty(k), error('MoC:Triple','No attached weak shock pressure match for these inputs.'); end
b=fzero(@residual,beta(k:k+1));
[d,P4,M4,ratioRho]=shock(b);
if M4<=1, error('MoC:Triple','Post-shock flow is subsonic; supersonic MOC is inapplicable.'); end
M3=MoC_PM(MoC_PM(seedMach,g)+d,g,true);
t.M3=M3; t.M4=M4; t.M5=M5; t.beta=b; t.delta=d;
t.slipAngle=geom.zeta+d; t.shockAngle=geom.zeta+b;
t.P=P4; t.T5=t.gas.T0/(1+(g-1)*M5^2/2);
t.T4=t.T5*(P4/inj.P)/ratioRho;
t.shockGas=struct('gamma',g,'R',p.R,'T0',t.gas.T0, ...
    'P0',P4*(1+(g-1)*M4^2/2)^(g/(g-1)));
t.pressureResidual=abs(residual(b))/P4;
t.seedMach=seedMach;
t.boundingPressure=inj.P;
t.betaMax=b;
for j=1:numel(beta)
    [~,~,candidateMach]=shock(beta(j));
    if candidateMach>1.001, t.betaMax=beta(j); end
end
    function r=residual(b)
        [d,ps]=shock(b);
        m=MoC_PM(MoC_PM(seedMach,g)+d,g,true);
        r=t.gas.P0/(1+(g-1)*m^2/2)^(g/(g-1))-ps;
    end
    function [d,ps,m,rho]=shock(b)
        n=M5^2*sin(b)^2;
        d=atan(2*cot(b)*(n-1)/(M5^2*(g+cos(2*b))+2));
        ps=inj.P*(1+2*g/(g+1)*(n-1));
        rho=(g+1)*n/((g-1)*n+2);
        m=sqrt((1+(g-1)*n/2)/(g*n-(g-1)/2))/sin(b-d);
    end
end
